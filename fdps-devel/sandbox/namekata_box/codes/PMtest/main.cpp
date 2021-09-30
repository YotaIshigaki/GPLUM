/* Standard headers */
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
/* Numerical libraries */
#include <gsl/gsl_rng.h>
/* FDPS headers */
#include <particle_simulator.hpp>
#include <param.h>
#include <param_fdps.h>
/* User-defined headers */
#include "preprocess_keywords.h"
#include "Simulation_Units.h"
#include "Calculation_Conditions.h"
#include "Nbody_Objects_Class.h"
#include "Misc_Class.h"
#include "Initialize.h"
#include "Error_Analysis.h"

/*-------------------------------------------------------------------*/
////////////////////// M A I N   F U N C T I O N //////////////////////
/*-------------------------------------------------------------------*/
int main(int argc, char* argv[]) {
   using namespace simulation_units;
   using namespace calculation_conditions;

   //* Local parameters
   const PS::S32 numSample=512;
   const bool flag_test1=false;
   const bool flag_test2=false;
   const bool flag_test3=true;
   //* Local variables
   Nbody_Objects Nbody_objs;
   Crystal_Parameters NaCl_params;
   //-(GSL)
   gsl_rng *r;
   r = gsl_rng_alloc(gsl_rng_ranlux389);

   //* Initialize 
   PS::Initialize(argc,argv);
   setup_sim_units();

   //==================================================================
   //* [Test1] Compute relative energy errors of the Madelung energy
   //          due to the P^{3}M method for different # of particles
   //          and for different configurations.
   //==================================================================
   if (flag_test1 == true) {
   for (PS::S32 numPtcl1D=4; numPtcl1D <= 32; numPtcl1D += 2) {
      NaCl_params.numPtcl_per_side = numPtcl1D;
      PS::F64 relerr {0.0};
      for (PS::S32 nstep=1; nstep <= numSample; nstep++) {
         if (PS::Comm::getRank() == 0)
            std::cout << "nstep = " << nstep << std::endl;

         //* [1] Randomly choose a configuration of the grid
         NaCl_params.pos_vertex.x = gsl_rng_uniform(r);
         NaCl_params.pos_vertex.y = gsl_rng_uniform(r);
         NaCl_params.pos_vertex.z = gsl_rng_uniform(r);

         //* [2] Compute force and potential with P^{3}M method
         Initialize(Nbody_objs,NaCl_params,nstep);

         //* [3] Compare with the result of the Ewald summation
#if (NaCl_Crystal_Mode == 1)
         relerr += calc_NaCl_error(Nbody_objs.system,nstep);
#else
         //error_analysis(Nbody_objs.system);
#endif
      }
      relerr /= numSample;

      //* Output relative error
      if (PS::Comm::getRank() == 0) {
         //* Output to file
         std::string filename {"EnergyError_test1.dat"};
         std::ofstream output_file;
         output_file.open(filename.c_str(),std::ios::out | std::ios::app);
         output_file << std::setprecision(15)
                     << numPtcl1D << " "
                     << relerr << std::endl;
         output_file.close();
         //* Output to STDOUT
         std::cout << "********** Result of this experiment **********" << std::endl;
         std::cout << "   numSample        = " << numSample << std::endl;
         std::cout << "   numPtcl_per_side = " << numPtcl1D << std::endl;
         std::cout << "   Relative Error   = " << relerr    << std::endl;
         std::cout << "***********************************************" << std::endl;
      }

   } // END of numPtcl1D
   } // END of if statement

   //==================================================================
   //* [Test2] Compute relative energy errors of the Madelung energy
   //          due to the P^{3}M method for different SIZE_OF_MESH
   //          and for different CUTOFF_RADIUS.
   //==================================================================
   std::vector<PS::S32> numPtclList;
   numPtclList.push_back(4);
   numPtclList.push_back(10);
   if (flag_test2 == true) {
   for (std::vector<PS::S32>::iterator iter=numPtclList.begin(); iter != numPtclList.end(); ++iter) {
      PS::S32 numPtcl1D = *iter;
      NaCl_params.numPtcl_per_side = numPtcl1D;
      PS::F64 relerr {0.0};
      for (PS::S32 nstep=1; nstep <= numSample; nstep++) {
         if (PS::Comm::getRank() == 0)
            std::cout << "nstep = " << nstep << std::endl;

         //* [1] Randomly choose a configuration of the grid
         NaCl_params.pos_vertex.x = gsl_rng_uniform(r);
         NaCl_params.pos_vertex.y = gsl_rng_uniform(r);
         NaCl_params.pos_vertex.z = gsl_rng_uniform(r);

         //* [2] Compute force and potential with P^{3}M method
         Initialize(Nbody_objs,NaCl_params,nstep);

         //* [3] Compare with the result of the Ewald summation
#if (NaCl_Crystal_Mode == 1)
         relerr += calc_NaCl_error(Nbody_objs.system,nstep);
#else
         //error_analysis(Nbody_objs.system);
#endif
      }
      relerr /= numSample;

      //* Output relative error
      if (PS::Comm::getRank() == 0) {
         //* Output to file
         std::string filename;
         std::ostringstream ss;
         ss << "EnergyError_np" 
            << std::setfill('0') << std::setw(2) << numPtcl1D
            << "_test2.dat";
         filename = ss.str();
         std::ofstream output_file;
         output_file.open(filename.c_str(),std::ios::out | std::ios::app);
         output_file << std::setprecision(15)
                     << SIZE_OF_MESH << " "
                     << ParticleSimulator::ParticleMesh::CUTOFF_RADIUS << " "
                     << relerr << std::endl;
         output_file.close();
      }

   } // END of numPtcl1D
   } // END of if statement

   //==================================================================
   //* [Test3] Compute relative energy errors of the Madelung energy
   //          due to the P^{3}M method in the case of
   //          CUTOFF_RADIUS = const. x numPtcl1D x SIZE_OF_MESH.
   //==================================================================
   if (flag_test3 == true) {
      PS::S32 numPtcl1D = 128;
      NaCl_params.numPtcl_per_side = numPtcl1D;
      PS::F64 relerr {0.0};
      for (PS::S32 nstep=1; nstep <= numSample; nstep++) {
         if (PS::Comm::getRank() == 0)
            std::cout << "nstep = " << nstep << std::endl;

         //* [1] Randomly choose a configuration of the grid
         NaCl_params.pos_vertex.x = gsl_rng_uniform(r);
         NaCl_params.pos_vertex.y = gsl_rng_uniform(r);
         NaCl_params.pos_vertex.z = gsl_rng_uniform(r);

         //* [2] Compute force and potential with P^{3}M method
         Initialize(Nbody_objs,NaCl_params,nstep);

         //* [3] Compare with the result of the Ewald summation
#if (NaCl_Crystal_Mode == 1)
         relerr += calc_NaCl_error(Nbody_objs.system,nstep);
#else
         //error_analysis(Nbody_objs.system);
#endif
      }
      relerr /= numSample;

      //* Output relative error
      if (PS::Comm::getRank() == 0) {
         //* Compute rcut_coef
         PS::F64 rcut_coef = (PS::ParticleMesh::CUTOFF_RADIUS * numPtcl1D)/(SIZE_OF_MESH);
         //* Output to file
         std::string filename;
         std::ostringstream ss;
         ss << "EnergyError_SoM" 
            << std::setfill('0') << std::setw(2) << SIZE_OF_MESH
            << "_test3.dat";
         filename = ss.str();
         std::ofstream output_file;
         output_file.open(filename.c_str(),std::ios::out | std::ios::app);
         output_file << std::setprecision(15)
                     << numPtcl1D << " "
                     << rcut_coef << " "
                     << relerr << std::endl;
         output_file.close();
      }

   } // END of if statement

   //* Finalize main
   gsl_rng_free(r);
   Finalize();
   return 0;
}
