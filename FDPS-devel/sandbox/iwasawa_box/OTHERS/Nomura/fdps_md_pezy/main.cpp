#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<iostream>
#include<fstream>
#include<sys/time.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>
#include<unistd.h>
#include "main.hpp"

#include "timer.h"
Profile timer;

#include "kernel.hpp"


#define MULTI_WALK

const int unit_num = 40;

template <class Tpsys>
void integrate_vel(Tpsys &system, const PS::F64 dt){
  const int np = system.getNumberOfParticleLocal();
  for(int i=0;i<np;i++) system[i].vel += system[i].acc * dt;
}

template <class Tpsys>
void integrate_pos(Tpsys &system, const PS::F64 dt){
  const int np = system.getNumberOfParticleLocal();
  for(int i=0;i<np;i++) system[i].pos += system[i].vel * dt;
}

// This is the program of MD, LJ particles, with using FDPS ///
int main(int argc, char *argv[]){
  PS::Initialize(argc,argv);
  int processes = PS::Comm::getNumberOfProc();
  int threads = PS::Comm::getNumberOfThread();
  int myrank = PS::Comm::getRank();
  if(myrank == 0) {
    fprintf(stdout, "Number of processes: %d\n", processes);
    fprintf(stdout, "Number of threads per process: %d\n", threads);
  }
  double density_coexist = 0.8; //Argon

  PS::S32 np, n_tot = 4 * unit_num * unit_num * unit_num;

  double cutoff = 4.5;
  //double cutoff = 9.0;
  //double cutoff = 13.5;

  double density = density_coexist;
  PS::F32vec side,sideh;
  /* Cell sizes */
  side.x = pow(n_tot/density, 1./3.);
  side.y = side.x;
  side.z = side.x;
  sideh  = side*.5;
  double unit_length = side.x/unit_num;

  PS::F32 cut_off2 = cutoff*cutoff;
  PS::F64 theta = 0.5;
  const PS::S32 n_leaf_limit = 8;
  //PS::S32 n_group_limit = 64;
  //PS::S32 n_group_limit = 128;
  //PS::S32 n_group_limit = 256;
  PS::S32 n_group_limit = 512;
  char dir_name[1024];

  PS::ParticleSystem<FPLJ> system_lj;
  system_lj.initialize();

  /* Initial Pos */
  /*for read at rank 0 only */
  if( myrank == 0 ){
    system_lj.setNumberOfParticleLocal(n_tot);
    set_initial_pos(system_lj,unit_num,unit_length,sideh);
    set_initial_vel(system_lj);
    /*
    for(int i=0;i<n_tot;i++) printf("%d %f %f %f %f %f %f\n",i
				    ,system_lj[i].pos.x,system_lj[i].pos.y,system_lj[i].pos.z
				    ,system_lj[i].vel.x,system_lj[i].vel.y,system_lj[i].vel.z);
    //*/
  }else{
    system_lj.setNumberOfParticleLocal(0);
  }
  np = system_lj.getNumberOfParticleLocal();

  /* FDPS domein and tree force */
  const PS::F64 coef_ema = 0.3; //defolt 1.0
  PS::DomainInfo dinfo;
  dinfo.initialize(coef_ema);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F32vec(-sideh.x,-sideh.y,-sideh.z),
  			 PS::F32vec( sideh.x, sideh.y, sideh.z));
  dinfo.collectSampleParticle(system_lj);
  dinfo.decomposeDomain();
  system_lj.exchangeParticle(dinfo);
  np = system_lj.getNumberOfParticleLocal();
  for(int i = 0; i<np; i++){
    system_lj[i].search_radius = cutoff;
    system_lj[i].acc = 0.0d;
  }
  ForceLJ *pezy = new ForceLJ[np];
  ForceLJ *cpu  = new ForceLJ[np];

  /* pezy test part */
#ifdef ENABLE_PEZY
  InitializeDEVICE();
#endif
  PS::TreeForForceShort<ForceLJ, EPILJ, EPJLJ>::Scatter tree_pezy;
  PS::TreeForForceShort<ForceLJ, EPILJ, EPJLJ>::Scatter tree_cpu;
  tree_pezy.initialize(np, theta, n_leaf_limit, n_group_limit);
  tree_cpu.initialize(np, theta, n_leaf_limit, n_group_limit);
  //n_group_limit

  const PS::S32 n_walk_limit = N_WALK_LIMIT;
  const PS::S32 tag_max = 1;
  std::cerr<<"n_walk_limit="<<n_walk_limit<<std::endl;

  if(myrank == 0){
#ifdef ENABLE_PEZY
    printf(" PEZY is ON\n");
#else
    printf(" PEZY is OFF\n");
#endif
  }

  dinfo.collectSampleParticle(system_lj);
  dinfo.decomposeDomainAll(system_lj);
  system_lj.exchangeParticle(dinfo);
  tree_cpu.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSPcpu,
					      RetrieveKernelcpu,
					      tag_max,
					      system_lj,
					      dinfo,
					      n_walk_limit,
					      true);
  PS::F64 pot,kin,tot0;
  calc_energy(system_lj,pot,kin);
  tot0 = pot + kin;

  const PS::F64 dt =0.0005;
  const int nstep = 5;

  timer.flush();
  //for(int s=-2000;s<nstep;s++){
  tree_pezy.clearTimeProfile();
  tree_cpu.clearTimeProfile();
  for(int s=0;s<nstep;s++){
    if(s<0) velocity_scaling(system_lj,1.0);
    PS::Comm::barrier();
    timer.beg(Profile::INTEG);
    integrate_vel(system_lj,0.5*dt);
    integrate_pos(system_lj,dt);
    PBC(system_lj,sideh,side);
    timer.end(Profile::INTEG);

    /////* Calc Force *//////
    PS::Comm::barrier();
    timer.beg(Profile::DECOMP);
    dinfo.collectSampleParticle(system_lj);
    dinfo.decomposeDomainAll(system_lj);
    timer.end(Profile::DECOMP);

    PS::Comm::barrier();
    timer.beg(Profile::EXCHANGE);
    system_lj.exchangeParticle(dinfo);
    timer.end(Profile::EXCHANGE);

    np = system_lj.getNumberOfParticleLocal();
    // calc force by pezy-sc
    PS::Comm::barrier();
    timer.beg(Profile::PEZY);
#ifdef MULTI_WALK
  tree_pezy.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
						RetrieveKernel,
						tag_max,
						system_lj,
						dinfo,
						n_walk_limit,
						true);
#else
    tree_peay.calcForceAllAndWriteBack(CalcForceEpEpPEZY(cut_off2),
				       system_lj, dinfo);
#endif
    timer.end(Profile::PEZY);
    for(int i=0;i<np;i++) pezy[i].acc = system_lj[i].acc;
    //calc force by host cpu
    PS::Comm::barrier();
    timer.beg(Profile::CPU);
#ifdef MULTI_WALK
    tree_cpu.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSPcpu,
						RetrieveKernelcpu,
						tag_max,
						system_lj,
						dinfo,
						n_walk_limit,
						true);
#else
    tree_cpu.calcForceAllAndWriteBack(CalcForceEpEp(cut_off2),
				       system_lj, dinfo);
#endif
    timer.end(Profile::CPU);
    for(int i=0;i<np;i++) cpu[i].acc = system_lj[i].acc;

    PS::Comm::barrier();
    timer.beg(Profile::INTEG);
    integrate_vel(system_lj,0.5*dt);
    timer.end(Profile::INTEG);

    PS::Comm::barrier();
    timer.beg(Profile::OTHER);
    calc_energy(system_lj,pot,kin);
    if(s==0) tot0 = pot + kin;
    if(myrank == 0){
      PS::F64vec max_error = 0.0;
      PS::F64vec min_error = 100.0;
      int err_count = 0;
      for(int i=0;i<np;i++){
	PS::F64vec error;
	error.x = fabs((cpu[i].acc.x - pezy[i].acc.x) / cpu[i].acc.x);
	error.y = fabs((cpu[i].acc.y - pezy[i].acc.y) / cpu[i].acc.y);
	error.z = fabs((cpu[i].acc.z - pezy[i].acc.z) / cpu[i].acc.z);
	if(max_error.x < error.x) max_error.x = error.x;
	if(max_error.y < error.y) max_error.y = error.y;
	if(max_error.z < error.z) max_error.z = error.z;
	if(min_error.x > error.x) min_error.x = error.x;
	if(min_error.y > error.y) min_error.y = error.y;
	if(min_error.z > error.z) min_error.z = error.z;
	if(error.x > 1e-6 || error.y > 1e-6 || error.z > 1e-6){
	  if(err_count<32)
	    printf("%d: %e %e %e, %e %e %e\n",i,
		   cpu[i].acc.x,cpu[i].acc.y,cpu[i].acc.z,
		   pezy[i].acc.x,pezy[i].acc.y,pezy[i].acc.z);
	  err_count++;
	}
      }
      printf("%d %e %e %e %e %e %e %e\n",
	     s,max_error.x,max_error.y,max_error.z,pot,kin,pot+kin,(pot+kin-tot0)/tot0);
#if 0
      if(max_error.x > 1e-5 || max_error.y > 1e-5 || max_error.z > 1e-5){
	for(int i=0;i<np;i++)
	  printf("force %d cpu, pezy: %e %e %e, %e %e %e\n",
		 i,cpu[i].acc.x,cpu[i].acc.y,cpu[i].acc.z
		 ,pezy[i].acc.x,pezy[i].acc.y,pezy[i].acc.z);
      }
#endif
    }
    timer.end(Profile::OTHER);
  }
#ifdef TUNING
  PS::TimeProfile time_profile = tree_pezy.getTimeProfile();
  printf("Total time of TreeForForceShort for pezy is %e\n",time_profile.getTotalTime());
  printf("calc_force:\t\t\t\t\t%e \n",time_profile.calc_force);
  printf("calc_force__core:\t\t\t\t%e \n",time_profile.calc_force__core);
  printf("calc_force__core__walk_tree:\t\t\t%e \n",time_profile.calc_force__core__walk_tree);
  printf("time_profile_.calc_force__copy_original_order:\t%e \n",time_profile.calc_force__copy_original_order);
  
  double time = timer.time[Profile::KERNEL];
  printf("%lu flop  / %lf s = %lf Gflops\n",flop,time, (double)flop / time / 1000000000.0);
#endif
  timer.show();
  PS::Finalize();
  return 0;
}
