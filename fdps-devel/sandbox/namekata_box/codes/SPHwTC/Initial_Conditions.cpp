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
#include "Misc_Class.h"
#include "Physical_Constants.h"
#include "Simulation_Units.h"
#include "Calculation_Conditions.h"
#include "IO.h"
#include "Initial_Conditions.h"

/*-------------------------------------------------------------------*/
////////////////////      S U B R O U T I N E      ////////////////////
//////////////////// < M A K E _ G L A S S _ I C > ////////////////////
/*-------------------------------------------------------------------*/
void PM_test_IC(PS::ParticleSystem<SPH_FP>& system,
                PS::DomainInfo& dinfo) {
   using namespace physical_constants;
   using namespace simulation_units;
   using namespace calculation_conditions;
   //* Local variables
   PS::S32 nprocs,myrank;
   //-(lattice)
   PS::F64 Lx,Ly,Lz;
   PS::F64 xmin,xmax,ymin,ymax,zmin,zmax;
   PS::F64 mu,gameff,rho,Tgas,pres,cs,tSC;
   PS::F64 totalMass;
   PS::S32 numPtcl,numPtclLocal,numPtclRem;
   std::vector<SPH_FP> ptcl;
   SPH_FP tmp;
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
   numPtcl = std::pow(static_cast<double>(2),18);
   //numPtcl = std::pow(static_cast<double>(2),24);
   Lx = Ly = Lz = 1.0e0 * (pc / unit_leng);
   mu     = 1.0e0;
   gameff = 5.0e0/3.0e0;
   rho    = 1.0e4 * M_H / unit_dens;
   Tgas   = 1.0e3 / unit_temp;
   pres   = rho * Tgas / mu ;
   cs     = std::sqrt(gameff * Tgas / mu);
   //** output time, etc.
   tSC = std::sqrt(Lx*Lx + Ly*Ly + Lz*Lz) / cs;
   IO_ctrl.set_config(10.0e0*tSC, tSC);
   if (PS::Comm::getRank() == 0) {
      std::cout << "tSC = " << tSC << std::endl;
   }
   //** error check for PM calculation
   if ((Lx > 1.0) || (Ly > 1.0) || (Lz > 1.0)) {
      if (PS::Comm::getRank() == 0)
         std::cout << "Invalid values for Lx,Ly,Lz !" << std::endl;
      PS::Finalize();
      std::exit(0);
   }

   //* Make a particle distribution
   xmin = 0.1e0; xmax = 0.9e0*Lz;
   ymin = 0.1e0; ymax = 0.9e0*Ly;
   zmin = 0.1e0; zmax = 0.9e0*Lz;
   for (PS::S32 i=0; i<numPtcl; i++) {
      tmp.x.x = xmin + (xmax-xmin)*gsl_rng_uniform(r);
      tmp.x.y = ymin + (ymax-ymin)*gsl_rng_uniform(r);
      tmp.x.z = zmin + (zmax-zmin)*gsl_rng_uniform(r);
      //* Make paticles be within the box
      if (tmp.x.x <  xmin) tmp.x.x += Lx;
      if (tmp.x.x >= xmax) tmp.x.x -= Lx;
      if (tmp.x.y <  ymin) tmp.x.y += Ly;
      if (tmp.x.y >= ymax) tmp.x.y -= Ly;
      if (tmp.x.z <  zmin) tmp.x.z += Lz;
      if (tmp.x.z >= zmax) tmp.x.z -= Lz;
      //* Double check
      if ((tmp.x.x < xmin) || (xmax <= tmp.x.x) ||
          (tmp.x.y < ymin) || (ymax <= tmp.x.y) ||
          (tmp.x.z < zmin) || (zmax <= tmp.x.z)) {
         if (myrank == 0) {
            std::cout << "A particle is out of the box !" 
                      << "(particle num. = " << i << ")" 
                      << std::endl;
         }
         std::exit(0);
      } 
      tmp.rho    = rho;
      tmp.pres   = pres;
      tmp.Tgas   = Tgas;
      tmp.mu     = mu;
      tmp.gameff = gameff;
      tmp.u      = Tgas/(mu*(gameff-1.0e0));
      tmp.id     = i;
      //tmp.h      = 3.0e0 * (3.0e0/(4.0e0*M_PI))
      //           * std::cbrt(numPtclNeighb*Lx*Ly*Lz/numPtcl);
      tmp.h      = 3.0e0 * (3.0e0/(4.0e0*M_PI))
                 * std::pow(numPtclNeighb*Lx*Ly*Lz/numPtcl,1.0e0/3.0e0);
      tmp.alpha  = 1.0e0;
      tmp.fAV    = 1.0e0;
      ptcl.push_back(tmp);
   }
   totalMass = rho*Lx*Ly*Lz;
   for (PS::S32 i=0; i<ptcl.size(); i++) 
      ptcl[i].m = totalMass/(PS::F64)(ptcl.size());
   if (myrank == 0) 
      std::cout << "# of ptcls = " << ptcl.size() << std::endl;

   //* [3] Scatter particles
   numPtclLocal = ptcl.size()/nprocs;
   if ((numPtclRem = ptcl.size() % nprocs) != 0) {
      if ((myrank+1) <= numPtclRem) {
         numPtclLocal++;
      }
   }
   std::cout << "# of local ptcls = " << numPtclLocal
             << "(" << myrank << ")" << std::endl;
   system.setNumberOfParticleLocal(numPtclLocal);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
   //* MPI-parallel
   std::vector<PS::S32> sendbuf(nprocs), ptcl_nums(nprocs);
   for (PS::S32 i=0; i<nprocs; i++) sendbuf[i]=numPtclLocal;
   MPI::COMM_WORLD.Alltoall(&sendbuf[0],1,MPI::INT, 
                            &ptcl_nums[0],1,MPI::INT);
   PS::S32 offset = 0;
   if (myrank > 0) {
      for (PS::S32 irank=0; irank<myrank; irank++)
         offset += ptcl_nums[irank];
   }
   for (PS::S32 i=offset; i<offset+numPtclLocal; i++) {
      PS::S32 ii = i - offset;
      system[ii] = ptcl[i];
   }
#else
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      system[i] = ptcl[i];
   }
#endif

   //* Initialize domain info and apply periodic boundary conditions
   dinfo.initialize();
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
   if (myrank == 0)
      std::cout << "[etime] PM = " << e_time << std::endl;
   gettimeofday(&start_time, NULL);
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      PS::F32vec x32 = system[i].x;
      PS::F32vec a32 = pm.getForce(x32);
      system[i].atot = a32;
   }
   gettimeofday(&end_time, NULL);
   e_time = (end_time.tv_sec  - start_time.tv_sec)
          + (end_time.tv_usec - start_time.tv_usec)*1.0e-6;
   if (myrank == 0)
      std::cout << "[etime] getForce = " << e_time << std::endl;
   //* Compute potential
   std::vector<PS::F64vec> aFD;
   PS::F64 ds = (1.0/32.0);
   gettimeofday(&start_time, NULL);
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
   gettimeofday(&end_time, NULL);
   e_time = (end_time.tv_sec  - start_time.tv_sec)
          + (end_time.tv_usec - start_time.tv_usec)*1.0e-6;
   if (myrank == 0)
      std::cout << "[etime] getPotential = " << e_time/6.0 << std::endl;
   //** Compute RMS force error
   PS::F64 RMS_error = 0.0;
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      PS::F64 scale = system[i].atot * system[i].atot;
      PS::F64vec da = system[i].atot - aFD[i];
      RMS_error += (da * da)/scale;
   }
   std::cout << std::sqrt(RMS_error)/numPtclLocal 
             << " (RMS, " << myrank << ")" 
             << std::endl;
   RMS_error = std::sqrt(PS::Comm::getSum(RMS_error)/numPtcl);
   if (myrank == 0) 
      std::cout << "RMS error = " << RMS_error << std::endl;
   //===================================================================

   //* Free objects
   gsl_rng_free(r);

   //* Terminate the program
   PS::Finalize();
   std::exit(0);
}

/*-------------------------------------------------------------------*/
////////////////////      S U B R O U T I N E      ////////////////////
//////////////////// < M A K E _ G L A S S _ I C > ////////////////////
/*-------------------------------------------------------------------*/
void make_glass_IC(PS::ParticleSystem<SPH_FP>& system,
                   PS::DomainInfo& dinfo,
                   Fluctuation_Monitor& f_monitor) {
   using namespace physical_constants;
   using namespace simulation_units;
   using namespace calculation_conditions;
   //* Local variables
   PS::S32 nprocs,myrank;
   //-(lattice)
   PS::F64 Lx,Ly,Lz;
   PS::F64 xmin,xmax,ymin,ymax,zmin,zmax;
   PS::F64 mu,gameff,rho,Tgas,pres,cs,tSC;
   PS::F64 totalMass;
   PS::S32 numPtcl,numPtclLocal,numPtclRem;
   std::vector<SPH_FP> ptcl;
   SPH_FP tmp;
   //-(GSL)
   gsl_rng *r;
   //-(IO)
   std::string filename;
   std::ofstream output_file;

   //* Get # of MPI processes and Rank number and
   //  initialize gsl_rng
   nprocs = PS::Comm::getNumberOfProc();
   myrank = PS::Comm::getRank();
   r = gsl_rng_alloc(gsl_rng_ranlux389);

   //* Define the parameters
   numPtcl = std::pow(static_cast<double>(2),18);
   Lx = Ly = Lz = 1.0e0 * pc / unit_leng;
   mu     = 1.0e0;
   gameff = 5.0e0/3.0e0;
   rho    = 1.0e4 * M_H / unit_dens;
   Tgas   = 1.0e3 / unit_temp;
   pres   = rho * Tgas / mu ;
   cs     = std::sqrt(gameff * Tgas / mu);
   //** output time, etc.
   tSC = std::sqrt(Lx*Lx + Ly*Ly + Lz*Lz) / cs;
   IO_ctrl.set_config(10.0e0*tSC, tSC);
   f_monitor.set_config(0.02e0*tSC, 0.02e0*tSC);
   if (PS::Comm::getRank() == 0) {
      std::cout << "tSC = " << tSC << std::endl;
   }

   //* Make a particle distribution
   xmin = -0.5e0*Lx; xmax = 0.5e0*Lz;
   ymin = -0.5e0*Ly; ymax = 0.5e0*Ly;
   zmin = -0.5e0*Lz; zmax = 0.5e0*Lz;
   for (PS::S32 i=0; i<numPtcl; i++) {
      tmp.x.x = xmin + (xmax-xmin)*gsl_rng_uniform(r);
      tmp.x.y = ymin + (ymax-ymin)*gsl_rng_uniform(r);
      tmp.x.z = zmin + (zmax-zmin)*gsl_rng_uniform(r);
      //* Make paticles be within the box
      if (tmp.x.x <  xmin) tmp.x.x += Lx;
      if (tmp.x.x >= xmax) tmp.x.x -= Lx;
      if (tmp.x.y <  ymin) tmp.x.y += Ly;
      if (tmp.x.y >= ymax) tmp.x.y -= Ly;
      if (tmp.x.z <  zmin) tmp.x.z += Lz;
      if (tmp.x.z >= zmax) tmp.x.z -= Lz;
      //* Double check
      if ((tmp.x.x < xmin) || (xmax <= tmp.x.x) ||
          (tmp.x.y < ymin) || (ymax <= tmp.x.y) ||
          (tmp.x.z < zmin) || (zmax <= tmp.x.z)) {
         if (myrank == 0) {
            std::cout << "A particle is out of the box !" 
                      << "(particle num. = " << i << ")" 
                      << std::endl;
         }
         std::exit(0);
      } 
      tmp.rho    = rho;
      tmp.pres   = pres;
      tmp.Tgas   = Tgas;
      tmp.mu     = mu;
      tmp.gameff = gameff;
      tmp.u      = Tgas/(mu*(gameff-1.0e0));
      tmp.id     = i;
      //tmp.h      = 3.0e0 * (3.0e0/(4.0e0*M_PI))
      //           * std::cbrt(numPtclNeighb*Lx*Ly*Lz/numPtcl);
      tmp.h      = 3.0e0 * (3.0e0/(4.0e0*M_PI))
                 * std::pow(numPtclNeighb*Lx*Ly*Lz/numPtcl,1.0e0/3.0e0);
      tmp.alpha  = 1.0e0;
      tmp.fAV    = 1.0e0;
      ptcl.push_back(tmp);
   }
   totalMass = rho*Lx*Ly*Lz;
   for (PS::S32 i=0; i<ptcl.size(); i++) 
      ptcl[i].m = totalMass/(PS::F64)(ptcl.size());
   if (myrank == 0) 
      std::cout << "# of ptcls = " << ptcl.size() << std::endl;

   //* [3] Scatter particles
   numPtclLocal = ptcl.size()/nprocs;
   if ((numPtclRem = ptcl.size() % nprocs) != 0) {
      if ((myrank+1) <= numPtclRem) {
         numPtclLocal++;
      }
   }
   std::cout << "# of local ptcls = " << numPtclLocal
             << "(" << myrank << ")" << std::endl;
   system.setNumberOfParticleLocal(numPtclLocal);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
   //* MPI-parallel
   std::vector<PS::S32> sendbuf(nprocs), ptcl_nums(nprocs);
   for (PS::S32 i=0; i<nprocs; i++) sendbuf[i]=numPtclLocal;
   MPI::COMM_WORLD.Alltoall(&sendbuf[0],1,MPI::INT, 
                            &ptcl_nums[0],1,MPI::INT);
   PS::S32 offset = 0;
   if (myrank > 0) {
      for (PS::S32 irank=0; irank<myrank; irank++)
         offset += ptcl_nums[irank];
   }
   for (PS::S32 i=offset; i<offset+numPtclLocal; i++) {
      PS::S32 ii = i - offset;
      system[ii] = ptcl[i];
   }
#else
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      system[i] = ptcl[i];
   }
#endif

   //* Initialize domain info and apply periodic boundary conditions
   dinfo.initialize();
   dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
   dinfo.setPosRootDomain(PS::F64vec(-0.5e0*Lx,-0.5e0*Ly,-0.5e0*Lz),
                          PS::F64vec(+0.5e0*Lx,+0.5e0*Ly,+0.5e0*Lz));
   //dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
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

   //* Free objects
   gsl_rng_free(r);
}

/*-------------------------------------------------------------------*/
////////////////////      S U B R O U T I N E      ////////////////////
//////////////////// < S H O C K _ T U B E _ I C > ////////////////////
/*-------------------------------------------------------------------*/
void shock_tube_IC(PS::ParticleSystem<SPH_FP>& system,
                   PS::DomainInfo& dinfo) {
   using namespace physical_constants;
   using namespace simulation_units;
   //* Local parameters
   const PS::F64 mu_fixed=28.8e0;
   const PS::F64 gam_fixed=1.4e0;
   //* Local variables
   PS::S32 nprocs,myrank;
   PS::F64 rho_L,pres_L,Tgas_L;
   PS::F64 rho_R,pres_R,Tgas_R;
   PS::F64 totalMass;
   //-(lattice)
   PS::F64 Lx,Ly,Lz; // boxsize
   PS::S32 nx_L,ny_L,nz_L,numPtcl_L;
   PS::F64 dx_L,dy_L,dz_L;
   PS::S32 nx_R,ny_R,nz_R,numPtcl_R;
   PS::F64 dx_R,dy_R,dz_R;
   PS::F64 xmin,xmax,ymin,ymax,zmin,zmax;
   PS::S32 numPtcl,numPtclLocal,numPtclRem;
   std::vector<SPH_FP> ptcl;
   SPH_FP tmp;
   //-(GSL)
   gsl_rng *r;
   //-(IO)
   std::string filename;
   std::ofstream output_file;

   //* Get # of MPI processes and Rank number and
   //  initialize gsl_rng
   nprocs = PS::Comm::getNumberOfProc();
   myrank = PS::Comm::getRank();
   r = gsl_rng_alloc(gsl_rng_ranlux389);

   //* Define the left- and right-states
   //* [A] Sod problem (Test1 in Toro's book; tstop=0.25)
   rho_L   = 1.0e0;
   pres_L  = 1.0e0;
   rho_R   = 0.125e0;
   pres_R  = 0.1e0;
   IO_ctrl.set_config(0.251e0/unit_time,0.025e0/unit_time);
   //* [B] Shock tube problem in Saitoh & Makino (2012) (tstop=0.1)
   //rho_L   = 1.0e0;
   //pres_L  = 1.0e0;
   //rho_R   = 0.25e0;
   //pres_R  = 0.1795e0;
   //IO_ctrl.set_config(0.1e0/unit_time,0.01e0/unit_time);
   //* [C] Strong shock tube problem in Saitoh & Makino (2012) (tstop=0.012)
   //rho_L   = 1.0e0;
   //pres_L  = 1000.0e0;
   //rho_R   = 1.0e0;
   //pres_R  = 0.01e0;
   //IO_ctrl.set_config(0.012e0/unit_time,0.0012e0/unit_time);
   //* [Common] Common settings
   Tgas_L  = (mu_fixed*M_H*pres_L)/(rho_L*kBoltz);
   Tgas_R  = (mu_fixed*M_H*pres_R)/(rho_R*kBoltz);

   //* Make a lattice of particles
   //* Define boxsize
   Lx = 2.0e0/unit_leng;
   Ly = Lz = Lx/16.0e0;
   //* Define # of particles in the left ([-0.5,0.0])
   nx_L = 128;
   dx_L = (0.5e0*Lx)/(PS::F64)nx_L;
   ny_L = nz_L = (PS::S32)(Ly/dx_L);
   dy_L = dz_L = Ly/(PS::F64)ny_L;
   numPtcl_L = nx_L * ny_L * nz_L;
   if (myrank  == 0) {
      std::cout << "*** The left region ***"      << std::endl;
      std::cout << "   numPtcl_L = " << numPtcl_L << std::endl;
      std::cout << "   nx_L      = " << nx_L      << std::endl;
      std::cout << "   ny_L      = " << ny_L      << std::endl;
      std::cout << "   nz_L      = " << nz_L      << std::endl;
      std::cout << "   dx_L      = " << dx_L      << std::endl;
      std::cout << "   dy_L      = " << dy_L      << std::endl;
      std::cout << "   dz_L      = " << dz_L      << std::endl;
   }
   //* Determine # of particles in the right ([0.0,0.5])
   numPtcl_R = (PS::S32)((rho_R/rho_L)*(PS::F64)numPtcl_L);
   nx_R = (PS::S32) std::pow((0.5e0*Lx/Ly)*(0.5e0*Lx/Lz)*numPtcl_R,1.0e0/3.0e0);
   dx_R = (0.5e0*Lx)/(PS::F64)nx_R;
   ny_R = nz_R = (PS::S32)(Ly/dx_R);
   dy_R = dz_R = Ly/(PS::F64)ny_R;
   numPtcl_R = nx_R * ny_R * nz_R;
   if (myrank  == 0) {
      std::cout << "*** The right region ***"     << std::endl;
      std::cout << "   numPtcl_R = " << numPtcl_R << std::endl;
      std::cout << "   nx_R      = " << nx_R      << std::endl;
      std::cout << "   ny_R      = " << ny_R      << std::endl;
      std::cout << "   nz_R      = " << nz_R      << std::endl;
      std::cout << "   dx_R      = " << dx_R      << std::endl;
      std::cout << "   dy_R      = " << dy_R      << std::endl;
      std::cout << "   dz_R      = " << dz_R      << std::endl;
   }
   //* Create a lattice
   xmin = -0.5e0*Lx; xmax = 0.0e0;
   ymin = -0.5e0*Ly; ymax = 0.5e0*Ly;
   zmin = -0.5e0*Lz; zmax = 0.5e0*Lz;
   numPtcl=0;
   for(PS::S32 i=0; i<nx_L; i++){
      for(PS::S32 j=0; j<ny_L; j++){
         for(PS::S32 k=0; k<nz_L; k++){
            numPtcl++;
            tmp.x.x    = (xmin + dx_L*((PS::F64)i + 0.5e0))
                       + (xmax - dx_L*((PS::F64)(nx_L-i) - 0.5e0));
            tmp.x.x    *= 0.5e0;
            tmp.x.y    = (ymin + dy_L*((PS::F64)j + 0.5e0))
                       + (ymax - dy_L*((PS::F64)(ny_L-j) - 0.5e0));
            tmp.x.y    *= 0.5e0;
            tmp.x.z    = (zmin + dz_L*((PS::F64)k + 0.5e0))
                       + (zmax - dz_L*((PS::F64)(nz_L-k) - 0.5e0));
            tmp.x.z    *= 0.5e0;
            tmp.rho    = rho_L/unit_dens;
            tmp.pres   = pres_L/unit_pres;
            tmp.Tgas   = Tgas_L/unit_temp;
            tmp.mu     = mu_fixed;
            tmp.gameff = gam_fixed;
            tmp.u      = tmp.Tgas/(tmp.mu*(tmp.gameff-1.0e0));
            tmp.id     = numPtcl;
            tmp.h      = 2.5*(dx_L+dx_R)*0.5e0;
            tmp.alpha  = 1.0e0;
            tmp.fAV    = 1.0e0;
            ptcl.push_back(tmp);
         }
      }
   }
   xmin = 0.0e0;     xmax = 0.5e0*Lx;
   ymin = -0.5e0*Ly; ymax = 0.5e0*Ly;
   zmin = -0.5e0*Lz; zmax = 0.5e0*Lz;
   for(PS::S32 i=0; i<nx_R; i++){
      for(PS::S32 j=0; j<ny_R; j++){
         for(PS::S32 k=0; k<nz_R; k++){
            numPtcl++;
            tmp.x.x    = (xmin + dx_R*((PS::F64)i + 0.5e0))
                       + (xmax - dx_R*((PS::F64)(nx_R-i) - 0.5e0));
            tmp.x.x    *= 0.5e0;
            tmp.x.y    = (ymin + dy_R*((PS::F64)j + 0.5e0))
                       + (ymax - dy_R*((PS::F64)(ny_R-j) - 0.5e0));
            tmp.x.y    *= 0.5e0;
            tmp.x.z    = (zmin + dz_R*((PS::F64)k + 0.5e0))
                       + (zmax - dz_R*((PS::F64)(nz_R-k) - 0.5e0));
            tmp.x.z    *= 0.5e0;
            tmp.rho    = rho_R/unit_dens;
            tmp.pres   = pres_R/unit_pres;
            tmp.Tgas   = Tgas_R/unit_temp;
            tmp.mu     = mu_fixed;
            tmp.gameff = gam_fixed;
            tmp.u      = tmp.Tgas/(tmp.mu*(tmp.gameff-1.0e0));
            tmp.id     = numPtcl;
            tmp.h      = 2.5*(dx_L+dx_R)*0.5e0;
            tmp.alpha  = 1.0e0;
            tmp.fAV    = 1.0e0;
            ptcl.push_back(tmp);
         }
      }
   }
   totalMass = (rho_L+rho_R)*0.5e0*Lx*Ly*Lz
             * (unit_leng*unit_leng*unit_leng)/unit_mass;
   for (PS::S32 i=0; i<ptcl.size(); i++) {
      ptcl[i].m = totalMass/(PS::F64)(ptcl.size());
      // Add position perturbation
      ptcl[i].x.x += 1.0e-8*Lx*(2.0e0*gsl_rng_uniform(r)-1.0e0);
      ptcl[i].x.y += 1.0e-8*Ly*(2.0e0*gsl_rng_uniform(r)-1.0e0);
      ptcl[i].x.z += 1.0e-8*Lz*(2.0e0*gsl_rng_uniform(r)-1.0e0);
   }
   if (myrank == 0) 
      std::cout << "# of ptcls = " << ptcl.size() << std::endl;
   //* Check particle distribution
   filename = "initial_distribution1.txt";
   output_file.open(filename.c_str(),std::ios::trunc);
   output_file.setf(std::ios_base::scientific,
                    std::ios_base::floatfield);
      for (PS::S32 i=0; i<ptcl.size(); i++) {
         output_file << std::setprecision(15)
                     << std::showpos
                     << ptcl[i].x.x * unit_leng << " "
                     << ptcl[i].x.y * unit_leng << " "
                     << ptcl[i].x.z * unit_leng << " "
                     << std::endl;
      }
   output_file.close();

   //* [3] Scatter particles
   numPtclLocal = ptcl.size()/nprocs;
   if ((numPtclRem = ptcl.size() % nprocs) != 0) {
      if ((myrank+1) <= numPtclRem) {
         numPtclLocal++;
      }
   }
   std::cout << "# of local ptcls = " << numPtclLocal
             << "(" << myrank << ")" << std::endl;
   system.setNumberOfParticleLocal(numPtclLocal);
   // [!!!Caution!!!] The following is valid only for serial processing.
   PS::S32 i_head = 0;
   PS::S32 i_tail = numPtclLocal;
   for (PS::S32 i=i_head; i<i_tail; i++) {
      PS::S32 ii = i - i_head;
      system[ii] = ptcl[i];
   }

   //* Initialize domain info and apply periodic boundary conditions
   dinfo.initialize();
   dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
   dinfo.setPosRootDomain(PS::F64vec(-0.5e0*Lx,-0.5e0*Ly,-0.5e0*Lz),
                          PS::F64vec(+0.5e0*Lx,+0.5e0*Ly,+0.5e0*Lz));
   dinfo.decomposeDomainAll(system);
   //* Exchange particle
   system.exchangeParticle(dinfo);
  
   //* Check exchanged particle distribution
   filename = "initial_distribution2.txt";
   output_file.open(filename.c_str(),std::ios::trunc);
   output_file.setf(std::ios_base::scientific,
                    std::ios_base::floatfield);
      for (PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++) {
         output_file << std::setprecision(15)
                     << std::showpos
                     << system[i].x.x * unit_leng << " "
                     << system[i].x.y * unit_leng << " "
                     << system[i].x.z * unit_leng << " "
                     << std::endl;
      }
   output_file.close();
   if (myrank == 0) {
      std::cout << "complete to make I.C.!" << std::endl;
   }

   //* Free objects
   gsl_rng_free(r);
}

