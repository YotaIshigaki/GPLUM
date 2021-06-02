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
#include <param.h>
#include <param_fdps.h>
/* User-defined headers */
#include "preprocess_keywords.h"
#include "Particle_Class.h"
#include "Misc_Class.h"
#include "Physical_Constants.h"
#include "Simulation_Units.h"
#include "Calculation_Conditions.h"
#include "Initial_Conditions.h"

/*-------------------------------------------------------------------*/
///////////////////////// S U B R O U T I N E /////////////////////////
/////////////////////////  < N A C L _ I C >  /////////////////////////
/*-------------------------------------------------------------------*/
void NaCl_IC(PS::ParticleSystem<Nbody_FP>& system,
             PS::DomainInfo& dinfo,
             Crystal_Parameters& params) {
   using namespace physical_constants;
   using namespace simulation_units;
   using namespace calculation_conditions;
   //* Local variables
   PS::S32 nprocs,myrank;
   //-(lattice)
   PS::F64 Lx,Ly,Lz;
   PS::F64 xmin,xmax,ymin,ymax,zmin,zmax;
   PS::F64 m,x,y,z,eps;
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

   //* Define the parameters
   numPtcl = params.numPtcl_per_side
           * params.numPtcl_per_side
           * params.numPtcl_per_side;
   if (params.numPtcl_per_side % 2 != 0) {
      std::cout << "numPtcl_per_side is an invalid value: "
                << params.numPtcl_per_side 
                << std::endl;
      PS::Finalize();
      std::exit(1);
   }
   Lx = Ly = Lz = 1.0e0; 
   eps = 1.0/256.0; // gravitatioal softening
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
   PS::S32 id {0};
   for(int i=0; i<params.numPtcl_per_side; i++) {
      for(int j=0; j<params.numPtcl_per_side; j++) {
         for(int k=0; k<params.numPtcl_per_side; k++){
            //* Substitute
            if ((i_start <= id) && (id <= i_end)) {
               PS::S32 ii = id - i_start;
               system[ii].id  = id;
               system[ii].m   = ((i+j+k)%2 ? 1.0 : -1.0)
                              / (params.numPtcl_per_side
                                *params.numPtcl_per_side);
               system[ii].eps = eps;
               system[ii].rc  = PS::ParticleMesh::CUTOFF_RADIUS / SIZE_OF_MESH; // see param.h & param_fdps.h
               system[ii].x.x = (params.pos_vertex.x + i) * (1.0/params.numPtcl_per_side);
               system[ii].x.y = (params.pos_vertex.y + j) * (1.0/params.numPtcl_per_side);
               system[ii].x.z = (params.pos_vertex.z + k) * (1.0/params.numPtcl_per_side);
               system[ii].v   = 0.0;
            } // END of i_start & i_end

            //* Update id
            id++;
         }
      }
   }

   //* Make the system charge-neutral
#if (Charge_Neutrality == 1)
   PS::F64 msum {0.0};
   for (PS::S32 i=0; i<numPtclLocal; i++) msum+=system[i].m;
   msum = PS::Comm::getSum(msum);
   for (PS::S32 i=0; i<numPtclLocal; i++) system[i].m -= msum/numPtcl;
   //* re-check
   msum = 0.0;
   for (PS::S32 i=0; i<numPtclLocal; i++) msum+=system[i].m;
   msum = PS::Comm::getSum(msum);
   if (myrank == 0) std::cout << "[check] msum = " << msum << std::endl;
#endif

   //* Initialize domain info and apply periodic boundary conditions
   dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
   dinfo.setPosRootDomain(PS::F64vec(0.0,0.0,0.0),
                          PS::F64vec(Lx, Ly, Lz));
   dinfo.decomposeDomainAll(system);
   //* Exchange particle
   system.exchangeParticle(dinfo);
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
//#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
//   std::ostringstream ss;
//   ss << "IC" << std::setfill('0') << std::setw(5) << myrank << ".txt";
//   filename = ss.str();
//#else
//   filename = "IC.txt";
//#endif
//   output_file.open(filename.c_str(),std::ios::trunc);
//   output_file.setf(std::ios_base::scientific,
//                    std::ios_base::floatfield);
//   for (PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++) {
//      output_file << std::setprecision(15)
//                  << std::showpos
//                  << system[i].x.x << " "
//                  << system[i].x.y << " "
//                  << system[i].x.z << " "
//                  << std::endl;
//   }
//   output_file.close();

   if (myrank == 0) 
      std::cout << "complete to make I.C.!" << std::endl;

}
/*-------------------------------------------------------------------*/
////////////////////      S U B R O U T I N E      ////////////////////
//////////////////// < M A K E _ G L A S S _ I C > ////////////////////
/*-------------------------------------------------------------------*/
void PM_test_IC(PS::ParticleSystem<Nbody_FP>& system,
                PS::DomainInfo& dinfo,
                Crystal_Parameters& params) {
   using namespace physical_constants;
   using namespace simulation_units;
   using namespace calculation_conditions;
   //* Local variables
   PS::S32 nprocs,myrank;
   //-(lattice)
   PS::F64 Lx,Ly,Lz;
   PS::F64 xmin,xmax,ymin,ymax,zmin,zmax;
   PS::F64 m,x,y,z,eps;
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
   numPtcl = params.numPtcl_per_side
           * params.numPtcl_per_side
           * params.numPtcl_per_side;
   Lx = Ly = Lz = 1.0e0; 
   eps = 1.0/256.0; // gravitatioal softening
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
   xmin = 0.0; xmax = Lx;
   ymin = 0.0; ymax = Ly;
   zmin = 0.0; zmax = Lz;
   for (PS::S32 i=0; i<numPtcl; i++) {
      //* Create a charge and a position
      // (1) Debug-purpose
      //if (i == 0) {
      //   //* Charged central particle
      //   m = 1.0;
      //} else {
      //   //* charge-less particle
      //   m = 0.0;
      //}
      //PS::F64vec end_point{1.0,1.0,1.0};
      //x = (static_cast<PS::F64>(i)/numPtcl) * end_point.x;
      //y = (static_cast<PS::F64>(i)/numPtcl) * end_point.y;
      //z = (static_cast<PS::F64>(i)/numPtcl) * end_point.z;
      // (2) Random charge + 1D distribution
      //PS::F64vec end_point{1.0,1.0,1.0};
      //m = gsl_rng_uniform(r);
      //x = (static_cast<PS::F64>(i)/numPtcl) * end_point.x;
      //y = (static_cast<PS::F64>(i)/numPtcl) * end_point.y;
      //z = (static_cast<PS::F64>(i)/numPtcl) * end_point.z;
      // (3) Random distribution
      m = gsl_rng_uniform(r);
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
            std::cout << "A particle is out of the box !" << std::endl; 
            std::cout << "  - ID = " << i << std::endl;
            std::cout << "  - x  = " << x << std::endl;
            std::cout << "  - y  = " << y << std::endl;
            std::cout << "  - z  = " << z << std::endl;
         }
         std::exit(0);
      } 
      //* Substitute
      if ((i_start <= i) && (i <= i_end)) {
         PS::S32 ii = i - i_start;
         system[ii].id     = i;
         system[ii].m      = m;
         system[ii].eps    = eps;
         system[ii].rc     = PS::ParticleMesh::CUTOFF_RADIUS / SIZE_OF_MESH; // see param.h & param_fdps.h
         system[ii].x.x    = x;
         system[ii].x.y    = y;
         system[ii].x.z    = z;
         system[ii].v      = 0.0;
      } // END of i_start & i_end
   }

   //* Make the system charge-neutral
#if (Charge_Neutrality == 1)
   PS::F64 msum {0.0};
   for (PS::S32 i=0; i<numPtclLocal; i++) msum+=system[i].m;
   msum = PS::Comm::getSum(msum);
   for (PS::S32 i=0; i<numPtclLocal; i++) system[i].m -= msum/numPtcl;
   //* re-check
   msum = 0.0;
   for (PS::S32 i=0; i<numPtclLocal; i++) msum+=system[i].m;
   msum = PS::Comm::getSum(msum);
   if (myrank == 0) std::cout << "[check] msum = " << msum << std::endl;
#endif

   //* Initialize domain info and apply periodic boundary conditions
   dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
   dinfo.setPosRootDomain(PS::F64vec(0.0,0.0,0.0),
                          PS::F64vec(Lx, Ly, Lz));
   dinfo.decomposeDomainAll(system);
   //* Exchange particle
   system.exchangeParticle(dinfo);
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
   output_file.open(filename.c_str(),std::ios::trunc);
   output_file.setf(std::ios_base::scientific,
                    std::ios_base::floatfield);
   for (PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++) {
      output_file << std::setprecision(15)
                  << std::showpos
                  << system[i].x.x << " "
                  << system[i].x.y << " "
                  << system[i].x.z << " "
                  << std::endl;
   }
   output_file.close();
   if (myrank == 0) 
      std::cout << "complete to make I.C.!" << std::endl;

   //* Free objects
   gsl_rng_free(r);

}

