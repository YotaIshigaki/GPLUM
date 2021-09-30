#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>
#include <float.h>
#include <cstdio>
#include <cstdlib>
#include <particle_simulator.hpp>
#include <particle_mesh.hpp>

#include "treepm.hpp"
#include "run_param.hpp"
#include "cosmology.hpp"
#include "mpi.h"

#include "prototype.h"

// unit_m = 1.212e48[kg/h] = 6.097e17[Msun/h]
// unit_L = 64[Mpc/h]
// unit_T = 1/(100h)[sec*Mpc/km]
// unit_v = 6400[km/sec]

template <class Tp>
void setup_random_particle(Tp &ptcl, 
                           run_param &this_run, 
                           const PS::S32 npart_total)
{
  PS::S32 rank = PS::Comm::getRank();
  PS::S64 npart_local = (rank == 0) ? npart_total : 0;

  ptcl.setNumberOfParticleLocal(npart_local);

  this_run.npart_total = npart_total;
  this_run.mpi_nproc = PS::Comm::getNumberOfProc();
  this_run.mpi_rank = rank;

  PS::MT::init_genrand(0);

  for(PS::S32 i=0;i<npart_local;i++) {
    ptcl[i].mass = 3.0*this_run.cosm.omegam/(8.0*AC::PI*(PS::F32)npart_total);

    ptcl[i].pos[0] = PS::MT::genrand_res53();
    ptcl[i].pos[1] = PS::MT::genrand_res53();
    ptcl[i].pos[2] = PS::MT::genrand_res53();
    ptcl[i].vel[0] = 0.0;
    ptcl[i].vel[1] = 0.0;
    ptcl[i].vel[2] = 0.0;
    ptcl[i].eps = 0.1/pow(npart_total,1.0/3.0);
  }

}

int main(int argc, char **argv)
{
  std::cout << std::setprecision(15);
  std::cerr << std::setprecision(15);

  PS::Initialize(argc, argv);

  PS::PM::ParticleMesh pm;
  PS::ParticleSystem<FPtreepm> ptcl;
  PS::DomainInfo domain_info;

  run_param this_run;

  //this_run.npart_total = 128*128*128;
  //std::cerr<<"this_run.npart_total="<<this_run.npart_total<<std::endl;
  //this_run.npart_total = 2097152;
  //  this_run.npart_total = 262144;
  //  this_run.npart_total = 32768;
  this_run.cosm.omegam = 1.0;
  this_run.cosm.omegav = 0.0;
  this_run.cosm.omegab = 0.0;
  this_run.cosm.omeganu = 0.0;

  this_run.step = 0;
  this_run.znow = 63.0;
  this_run.anow = 1.0 / (1.0 + this_run.znow);
  this_run.tnow = this_run.cosm.atotime(this_run.anow);
  std::cerr<<"this_run.znow="<<this_run.znow
	   <<" this_run.anow="<<this_run.anow
      	   <<" this_run.tnow="<<this_run.tnow<<std::endl;
  //this_run.tnow = 0.01;
  this_run.update_expansion(this_run.tnow);
  std::cerr<<"this_run.znow="<<this_run.znow
           <<" this_run.anow="<<this_run.anow
      	   <<" this_run.tnow="<<this_run.tnow<<std::endl;

#if 1
  ptcl.initialize();
  domain_info.initialize();

  //setup_random_particle(ptcl, this_run, this_run.npart_total);
  this_run.mpi_nproc = PS::Comm::getNumberOfProc();
  this_run.mpi_rank = PS::Comm::getRank();
  std::cerr<<"read particle"<<std::endl;
  //ptcl.readParticleBinary("./particles_ic_sb128");
  ptcl.readParticleBinary("./particles_ic_sb256");

  this_run.npart_total = ptcl.getNumberOfParticleGlobal();
  this_run.npart_local = ptcl.getNumberOfParticleLocal();
  ptcl[0].eps = 0.1 / pow(this_run.npart_total, 1.0/3.0);



  if(PS::Comm::getRank()==0){
      std::cerr<<"ptcl.getNumberOfParticleLocal0()="<<ptcl.getNumberOfParticleLocal()<<std::endl;
      for(PS::S32 i=0; i<10; i++){
	  std::cerr<<"ptcl0[i].mass="<<ptcl[i].mass<<std::endl;
	  std::cerr<<"ptcl0[i].pos="<<ptcl[i].pos<<std::endl;
	  std::cerr<<"ptcl0[i].vel="<<ptcl[i].vel<<std::endl;
      std::cerr<<"ptcl0[i].eps="<<ptcl[i].eps<<std::endl;
	  std::cerr<<"ptcl0[i].id="<<ptcl[i].id<<std::endl;
      }
  }
  if(PS::Comm::getRank() == PS::Comm::getNumberOfProc()-1 ){
      std::cerr<<"ptcl.getNumberOfParticleLocal1()="<<ptcl.getNumberOfParticleLocal()<<std::endl;
      for(PS::S32 i=ptcl.getNumberOfParticleLocal()-1; i>=ptcl.getNumberOfParticleLocal()-10; i--){
          std::cerr<<"ptcl1[i].mass="<<ptcl[i].mass<<std::endl;
          std::cerr<<"ptcl1[i].pos="<<ptcl[i].pos<<std::endl;
          std::cerr<<"ptcl1[i].vel="<<ptcl[i].vel<<std::endl;
          std::cerr<<"ptcl1[i].eps="<<ptcl[i].eps<<std::endl;
          std::cerr<<"ptcl1[i].id="<<ptcl[i].id<<std::endl;
      }
  }

  this_run.input_params(argv[1]);

  domain_info.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  domain_info.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0),
                               PS::F64vec(1.0, 1.0, 1.0));

  domain_info.decomposeDomainAll(ptcl);
  ptcl.adjustPositionIntoRootDomain(domain_info);
  ptcl.exchangeParticle(domain_info);
  this_run.npart_local = ptcl.getNumberOfParticleLocal();

  PS::TreeForForceLong<Result_treepm, EPItreepm, EPJtreepm>::MonopoleWithCutoff treepm_tree;

  treepm_tree.initialize(3*ptcl.getNumberOfParticleGlobal());

  treepm_tree.calcForceAllAndWriteBack
    (calc_pp_force<EPJtreepm>(),
     calc_pp_force<PS::SPJMonopoleCutoff>(),
     ptcl,
     domain_info);

  pm.calcForceAllAndWriteBack(ptcl, domain_info);
  PS::F64 dtime = calc_dtime(ptcl, this_run);
  PS::F64 dtime_prev, dtime_mid;

  std::cerr<<"dtime="<<dtime<<std::endl;

  this_run.output_diag(dtime);

  drift_ptcl(ptcl, domain_info, 0.5*dtime);



  domain_info.decomposeDomainAll(ptcl);
  ptcl.exchangeParticle(domain_info);
  this_run.npart_local = ptcl.getNumberOfParticleLocal();
  PS::S64 n_loop = 0;
  //while(this_run.tnow < 0.5) {
  while(this_run.tnow < 0.663) {
      std::cerr<<"n_loop="<<n_loop<<std::endl;
      treepm_tree.calcForceAllAndWriteBack
          (calc_pp_force<EPJtreepm>(),
           calc_pp_force<PS::SPJMonopoleCutoff>(),
           ptcl,
           domain_info);
    
      std::cerr<<"this_run.tnow="<<this_run.tnow<<std::endl;

      pm.calcForceAllAndWriteBack(ptcl, domain_info);

      std::cerr<<"this_run.znow="<<this_run.znow<<std::endl;

      this_run.tnow += 0.5*dtime;

      std::cerr<<"this_run.anow="<<this_run.anow<<std::endl;

      this_run.update_expansion(this_run.tnow);

      std::cerr<<"this_run.hnow="<<this_run.anow<<std::endl;

      kick_ptcl(ptcl, dtime, this_run);

      this_run.tnow += 0.5*dtime;

      this_run.update_expansion(this_run.tnow);

      dtime_prev = dtime;
      dtime = calc_dtime(ptcl, this_run);

      std::cerr<<"dtime="<<dtime<<std::endl;

      dtime_mid = 0.5*(dtime_prev + dtime);

      drift_ptcl(ptcl, domain_info, dtime_mid);

      domain_info.decomposeDomainAll(ptcl);
      ptcl.exchangeParticle(domain_info);
      this_run.npart_local = ptcl.getNumberOfParticleLocal();

      output_data_in_run(ptcl, this_run);

      this_run.step++;
      this_run.output_diag(dtime);

      n_loop++;

  }
#endif
  
  PS::Finalize();

  return(0);
}
