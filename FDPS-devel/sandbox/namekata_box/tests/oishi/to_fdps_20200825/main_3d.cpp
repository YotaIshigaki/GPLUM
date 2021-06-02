/*****************************************************************
* main.cpp
* using FDPS, OpenMP
* compiler: g++
* Date:
* Author: saki nonaka.
******************************************************************/

#include <particle_simulator.hpp>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <iomanip>
#include <time.h>
#include <cstdlib>// ->exit(0)

#include "operator.h"
#include "constant.h"
#include "fp_ep.h"
#include "force.h" // -> kernelfunc
#include "func.h"
#include "preprocess.h"





/******************************************************************************************/
/* MAIN                                                                                   */
/******************************************************************************************/

int main(int argc,char* argv[]){
  clock_t ex;
  ex = clock();

  //Initialize FDPS
  PS::Initialize(argc,argv);

  //make  directory
  char dir_name[1024];
  sprintf(dir_name, "./bin");
  makeOutputDirectory(dir_name);
  sprintf(dir_name, "./result");
  makeOutputDirectory(dir_name);

  //display # of MPI process and threads
  PS::S32 nprocs = PS::Comm::getNumberOfProc();
  PS::S32 nthrds = PS::Comm::getNumberOfThread();
  std::cout <<"process: "<<nprocs <<"   thread: "<<nthrds <<std::endl;
  //Make an instance of ParticleSystem and initialize it
  PS::ParticleSystem<FP> sph_system;
  sph_system.initialize();

  /***** make an initial condition *****/
  generate_ptcls(sph_system);

  //Make an instance of DomainInfo and initialize it
  PS::DomainInfo dinfo;
  dinfo.initialize();
  //set the boundary consition
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
  //dinfo.setPosRootDomain(PS::F64vec( -0.5,-0.5, ),PS::F64vec( 0.6,0.32));
  //perform domain decomposition
  dinfo.decomposeDomainAll(sph_system);
  //exchange the SPH particles between the (MPI) process
  sph_system.exchangeParticle(dinfo);
  //Make 2 tree structures(density,force calculation)
	// PS::TreeForForceShort<Dens,EP,EP>::Gather dens_tree;
	// dens_tree.initialize(3*sph_system.getNumberOfParticleGlobal());
	PS::TreeForForceShort<Hydro, EP, EP>::Gather hydr_tree;
  PS::F64 theta = 0.0;
  PS::U32 n_leaf_limit = 1;
  PS::U32 n_group_limit = 4;
	hydr_tree.initialize(3.0 * sph_system.getNumberOfParticleGlobal(), theta, n_leaf_limit, n_group_limit);
	// hydr_tree.initialize(3.0 * sph_system.getNumberOfParticleGlobal());
  //compute density, pressure, eleration at initial condition
  /*dens_tree.calcForceAllAndWriteBack(CalcDensity(),sph_system,dinfo);//density
  hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(),sph_system,dinfo);//force*/

  PS::S32 step = 0;
  ene << "time,potential,kinetic, elastic,plastic,total,H+K\n";

  clock_t pre;
  pre = clock() - ex;
  double pre_t;
  pre_t = (double) pre / CLOCKS_PER_SEC;

    FileHeader header0;
    header0.time  = 0.0;
    header0.Nbody = sph_system.getNumberOfParticleGlobal();
    char ascii0_name[256];
    sprintf( ascii0_name, "result/initial.csv");
    sph_system.writeParticleAscii( ascii0_name, header0);


  PS::F64 wgh = 1.0;


/********************************************************************/
/* MAIN LOOP                                                        */
/********************************************************************/

  for(PS::F64 time = 0 ; time < end_time ; time += dt, ++step){
    double time0, time1, time_fin;
    PS::Comm::barrier();
    time0 = MPI_Wtime();

    InitialKick(sph_system,dt);
    time_fin = MPI_Wtime();
    if(PS::Comm::getRank() == 0) std::cout<< "step="<<step<<" :initinal kick =  "<<time_fin-time0<<std::endl; 

    time1 = MPI_Wtime();
    sph_system.adjustPositionIntoRootDomain(dinfo);
    dinfo.decomposeDomainAll(sph_system);
    dinfo.decomposeDomainAll(sph_system, wgh);
    time_fin = MPI_Wtime();
    if(PS::Comm::getRank() == 0) std::cout<< "step="<<step<<" :domain =  "<<time_fin-time1<<std::endl; 

    time1 = MPI_Wtime();
    sph_system.exchangeParticle(dinfo);
    time_fin = MPI_Wtime();
    if(PS::Comm::getRank() == 0) std::cout<< "step="<<step<<" :exchange =  "<<time_fin-time1<<std::endl; 
    
    time1 = MPI_Wtime();
    EigenSolver(sph_system);
    time_fin = MPI_Wtime();
    if(PS::Comm::getRank() == 0) std::cout<< "step="<<step<<" :eigen solver =  "<<time_fin-time1<<std::endl;
    time1 = MPI_Wtime();
    hydr_tree.calcForceAllAndWriteBack(CalcHydroForce_3d(),sph_system,dinfo);
    time_fin = MPI_Wtime();
    if(PS::Comm::getRank() == 0) std::cout<< "step="<<step<<" :force =  "<<time_fin-time1<<std::endl; 

    time1 = MPI_Wtime();
    FinalKick(sph_system,dt);
    time_fin = MPI_Wtime();
    if(PS::Comm::getRank() == 0) std::cout<< "step="<<step<<" :final kick= "<<time_fin-time1<<std::endl; 
    if(PS::Comm::getRank() == 0) std::cout<< "step="<<step<<" :total= "<<time_fin-time0<<std::endl; 

    PS::TimeProfile tp = hydr_tree.getTimeProfile();
    wgh = tp.calc_force;
    hydr_tree.clearCounterAll();

    if(PS::Comm::getRank() == 0) std::cout<<"step = "<<step<<std::endl;

    //output result files
    if(step % OUTPUT_INTERVAL == 0){
      GetPRank(sph_system);
      FileHeader header;
      header.time  = time;
      header.Nbody = sph_system.getNumberOfParticleGlobal();
      char ascii_name[256], bin_name[256];
      //sprintf( ascii_name, "result/%04d.csv", step );
      sprintf( bin_name, "bin/%05d.bin", step );
      //sph_system.writeParticleAscii( ascii_name, header);
      sph_system.writeParticleBinary( bin_name, header);

      if(PS::Comm::getRank() == 0){
        //std::cout<<"==============="<<std::endl;
        std::cout<<"output  " <<bin_name <<"." <<std::endl;
        //std::cout<<"==============="<<std::endl;
      }
      CalcSum(sph_system , time);
    }
    //std::cout<<"time = "<<time<<std::endl;
  }

  ene.close();

  PS::Finalize();
  ex = clock() - ex;
  double cal_t;
  cal_t = (double)ex/CLOCKS_PER_SEC;
  std::cout<<"program excution time is (s)...."<<cal_t<<std::endl;
  return 0;
}
