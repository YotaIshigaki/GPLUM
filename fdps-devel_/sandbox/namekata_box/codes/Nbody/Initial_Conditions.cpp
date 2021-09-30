/* Standard headers */
#include <cstdio>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <sys/time.h>
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#endif
/* Numerical libraries */
#include <gsl/gsl_rng.h>
/* FDPS headers */
#include <particle_simulator.hpp>
#include <particle_mesh.hpp>
/* User-defined headers */
#include "preprocess_keywords.h"
#include "Particle_Class.h"
#include "Physical_Constants.h"
#include "Simulation_Units.h"
#include "Calculation_Conditions.h"
#include "IO.h"
#include "Initial_Conditions.h"

/*-------------------------------------------------------------------*/
////////////////////      S U B R O U T I N E      ////////////////////
//////////////////// < M A K E _ G L A S S _ I C > ////////////////////
/*-------------------------------------------------------------------*/
void PM_test_IC(PS::ParticleSystem<Nbody_FP>& system,
                PS::DomainInfo& dinfo) {
   using namespace physical_constants;
   using namespace simulation_units;
   using namespace calculation_conditions;
   //* Local variables
   PS::S32 nprocs,myrank;
   //-(lattice)
   PS::F64 Lx,Ly,Lz;
   PS::F64 xmin,xmax,ymin,ymax,zmin,zmax;
   PS::F64 x,y,z,eps;
   PS::F64 totalMass;
   PS::S32 numPtcl,numPtclLocal,numPtclRem;
   PS::S32 i_start,i_end;
   //-(GSL)
   gsl_rng *r;
   //-(IO)
   std::string filename;
   std::ofstream output_file;
   //-(Time measurement)
   struct timeval start_time,end_time;
   double e_time;

   //* Get # of MPI processes and Rank number and
   //  initialize gsl_rng
   nprocs = PS::Comm::getNumberOfProc();
   myrank = PS::Comm::getRank();
   r = gsl_rng_alloc(gsl_rng_ranlux389);

   //* Define the parameters
   //numPtcl = std::pow(static_cast<double>(2),18);
   //numPtcl = std::pow(static_cast<double>(2),24);
   //numPtcl = std::pow(static_cast<double>(2),27);
   //numPtcl = std::pow(static_cast<double>(2),30);
   numPtcl = std::pow(static_cast<double>(2),33);
   Lx = Ly = Lz = 1.0e0 * (pc / unit_leng);
   totalMass = 1.0e0 * Msolar / unit_mass;
   eps = 1.0/256.0; // gravitatioal softening
   IO_ctrl.set_config(0.0, 0.0); // tentatively
   //** error check for PM calculation
   if ((Lx > 1.0) || (Ly > 1.0) || (Lz > 1.0)) {
      if (PS::Comm::getRank() == 0)
         std::cout << "Invalid values for Lx,Ly,Lz !" << std::endl;
      PS::Finalize();
      std::exit(0);
   }

   //* Compute numPtclList
   std::vector<PS::S32> numPtclList(nprocs);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
   //* MPI-parallel
   numPtclLocal = numPtcl/nprocs;
   if ((numPtclRem = numPtcl % nprocs) != 0) {
      if ((myrank+1) <= numPtclRem) {
         numPtclLocal++;
      }
   }
   std::vector<PS::S32> sendbuf(nprocs);
   for (PS::S32 i=0; i<nprocs; i++) sendbuf[i]=numPtclLocal;
   MPI::COMM_WORLD.Alltoall(&sendbuf[0],1,MPI::INT, 
                            &numPtclList[0],1,MPI::INT);
   i_start = 0;
   if (myrank > 0) {
      for (PS::S32 irank=0; irank<myrank; irank++)
         i_start += numPtclList[irank];
   }
   i_end = i_start + numPtclLocal - 1;
#else
   numPtclList.push_back(numPtcl);
   numPtclLocal = numPtcl;
   i_start = 0;
   i_end   = i_start + numPtclLocal - 1;
#endif
   system.setNumberOfParticleLocal(numPtclLocal);

   //* Make a particle distribution
   xmin = 0.1e0; xmax = 0.9e0*Lz;
   ymin = 0.1e0; ymax = 0.9e0*Ly;
   zmin = 0.1e0; zmax = 0.9e0*Lz;
   for (PS::S32 i=0; i<numPtcl; i++) {
      if (myrank == 0) {
         if (i % 10000000 == 0) 
            std::cout << "i = " << i << std::endl;
      }
      //* Create a position
      x = xmin + (xmax-xmin)*gsl_rng_uniform(r);
      y = ymin + (ymax-ymin)*gsl_rng_uniform(r);
      z = zmin + (zmax-zmin)*gsl_rng_uniform(r);
      //* Make paticles be within the box
      if (x <  xmin) x += Lx;
      if (x >= xmax) x -= Lx;
      if (y <  ymin) y += Ly;
      if (y >= ymax) y -= Ly;
      if (z <  zmin) z += Lz;
      if (z >= zmax) z -= Lz;
      //* Double check
      if ((x < xmin) || (xmax <= x) ||
          (y < ymin) || (ymax <= y) ||
          (z < zmin) || (zmax <= z)) {
         if (myrank == 0) {
            std::cout << "A particle is out of the box !" 
                      << "(particle num. = " << i << ")" 
                      << std::endl;
         }
         std::exit(0);
      } 
      //* Substitute
      if ((i_start <= i) && (i <= i_end)) {
         PS::S32 ii = i - i_start;
         system[ii].id     = i;
         system[ii].m      = totalMass/numPtcl;
         system[ii].eps    = eps;
         system[ii].x.x    = x;
         system[ii].x.y    = y;
         system[ii].x.z    = z;
         system[ii].v      = 0.0;
      }
   }

   //* Initialize domain info and apply periodic boundary conditions
   if (myrank == 0) std::cout << "OK1" << std::endl;
   dinfo.initialize();
   dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
   dinfo.setPosRootDomain(PS::F64vec(0.0,0.0,0.0),
                          PS::F64vec(Lx, Ly, Lz));
   if (myrank == 0) std::cout << "OK2" << std::endl;
   dinfo.decomposeDomainAll(system);
   if (myrank == 0) std::cout << "OK3" << std::endl;
   //* Exchange particle
   system.exchangeParticle(dinfo);
   if (myrank == 0) std::cout << "OK4" << std::endl;
   //* Check load balance
   numPtclLocal = system.getNumberOfParticleLocal();
   {
      PS::S32 Min = PS::Comm::getMinValue(numPtclLocal);
      PS::S32 Max = PS::Comm::getMaxValue(numPtclLocal);
      PS::S32 Sum = PS::Comm::getSum(numPtclLocal);
      if (myrank == 0) {
         std::cout << "Max. " << Max << std::endl;
         std::cout << "Min. " << Min << std::endl;
         std::cout << "Sum. " << Sum << std::endl;
      }
   }
  
   //* Check exchanged particle distribution
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
   std::ostringstream ss;
   ss << "IC" << std::setfill('0') << std::setw(5) << myrank << ".txt";
   filename = ss.str();
#else
   filename = "IC.txt";
#endif
   //output_file.open(filename.c_str(),std::ios::trunc);
   //output_file.setf(std::ios_base::scientific,
   //                 std::ios_base::floatfield);
   //for (PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++) {
   //   output_file << std::setprecision(15)
   //               << std::showpos
   //               << system[i].x.x * unit_leng << " "
   //               << system[i].x.y * unit_leng << " "
   //               << system[i].x.z * unit_leng << " "
   //               << std::endl;
   //}
   //output_file.close();
   if (myrank == 0) 
      std::cout << "complete to make I.C.!" << std::endl;

   //===================================================================
   //* Test for Particle Mesh method
   //===================================================================
   //* Notification
   if (myrank == 0) 
      std::cout << "start to PM test !" << std::endl;
   //* Compute force
   numPtclLocal = system.getNumberOfParticleLocal();
   PS::PM::ParticleMesh pm;
   pm.setDomainInfoParticleMesh(dinfo);
   pm.setParticleParticleMesh(system);
   gettimeofday(&start_time, NULL);
   pm.calcMeshForceOnly();
   gettimeofday(&end_time, NULL);
   e_time = (end_time.tv_sec  - start_time.tv_sec)
          + (end_time.tv_usec - start_time.tv_usec)*1.0e-6;
   e_time = PS::Comm::getSum(e_time)/nprocs;
   if (myrank == 0)
      std::cout << "[etime] calcMeshForceOnly = " << e_time << std::endl;

   gettimeofday(&start_time, NULL);
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      PS::F32vec x32 = system[i].x;
      system[i].agrv = pm.getForce(x32);
   }
   gettimeofday(&end_time, NULL);
   e_time = (end_time.tv_sec  - start_time.tv_sec)
          + (end_time.tv_usec - start_time.tv_usec)*1.0e-6;
   e_time = PS::Comm::getSum(e_time)/nprocs;
   if (myrank == 0)
      std::cout << "[etime] getForce = " << e_time << std::endl;
   //* Compute potential
   std::vector<PS::F64vec> aFD;
   PS::F64 ds = (1.0/32.0);
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      PS::F32vec xm,xp; 
      PS::F32 phim,phip;
      PS::F64vec a;
      //** X
      xm = xp = system[i].x;
      xm.x -= 0.5 * ds;
      xp.x += 0.5 * ds;
      phim = pm.getPotential(xm);
      phip = pm.getPotential(xp);
      a.x = - (phip - phim)/(xp.x - xm.x);
      //** Y
      xm = xp = system[i].x;
      xm.y -= 0.5 * ds;
      xp.y += 0.5 * ds;
      phim = pm.getPotential(xm);
      phip = pm.getPotential(xp);
      a.y = - (phip - phim)/(xp.y - xm.y);
      //** Z
      xm = xp = system[i].x;
      xm.z -= 0.5 * ds;
      xp.z += 0.5 * ds;
      phim = pm.getPotential(xm);
      phip = pm.getPotential(xp);
      a.z = - (phip - phim)/(xp.z - xm.z);
      //** push back
      aFD.push_back(a);
   }
   //** Compute RMS force error
   PS::F64 RMS_error = 0.0;
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      PS::F64 scale = system[i].agrv * system[i].agrv;
      PS::F64vec da = system[i].agrv - aFD[i];
      RMS_error += (da * da)/scale;
   }
   std::cout << std::sqrt(RMS_error)/numPtclLocal 
             << " (RMS, " << myrank << ")" 
             << std::endl;
   RMS_error = std::sqrt(PS::Comm::getSum(RMS_error)/numPtcl);
   if (myrank == 0) 
      std::cout << "RMS error = " << RMS_error << std::endl;
   //** Measure elapsed time for getPotential
   gettimeofday(&start_time, NULL);
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      PS::F32vec x32 = system[i].x;
      system[i].pot = pm.getPotential(x32);
   }
   gettimeofday(&end_time, NULL);
   e_time = (end_time.tv_sec  - start_time.tv_sec)
          + (end_time.tv_usec - start_time.tv_usec)*1.0e-6;
   e_time = PS::Comm::getSum(e_time)/nprocs;
   if (myrank == 0)
      std::cout << "[etime] getPotential = " << e_time << std::endl;
   //===================================================================

   //* Free objects
   gsl_rng_free(r);

   //* Terminate the program
   PS::Finalize();
   std::exit(0);
}

