/* Standard headers */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <limits>
#include <vector>
/* FDPS headers */
#include <particle_simulator.hpp>
#include <particle_mesh.hpp>
/* User-defined headers */
#include "preprocess_keywords.h"
#include "Particle_Class.h"
#include "Nbody_Objects_Class.h"
#include "Calculation_Conditions.h"

//================================
//* Class Defs.: Nbody_Objects
//================================
void Nbody_Objects::init_tree() {
   PS::S32 numPtclLocal = system.getNumberOfParticleLocal();
   PS::U64 n_glb_tot = 3 * numPtclLocal;
   pp_tree.initialize(n_glb_tot,0.0);
}

void Nbody_Objects::calc_gravity() {
   //* Local variables
   PS::S32 numPtclLocal = system.getNumberOfParticleLocal();

   //* Reset potential and accelerations
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      system[i].pot  = 0.0;
      system[i].agrv = 0.0;
   }

   //=================
   //* [1] PM part 
   //=================
   pm.setDomainInfoParticleMesh(dinfo);
   pm.setParticleParticleMesh(system);
   pm.calcMeshForceOnly();
   for (PS::S32 i=0; i<numPtclLocal; i++) { 
      PS::F32vec x32 = system[i].x; 
      system[i].pot  += pm.getPotential(x32);
      system[i].agrv += pm.getForce(x32);
   }

   //=================
   //* [1] PP part 
   //=================
   pp_tree.calcForceAll(Calc_gravity(), system, dinfo);
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      Nbody_PP_Results result = pp_tree.getForce(i);
      system[i].pot  += result.pot;
      system[i].agrv += result.agrv;
   }

   //* Check
   //std::ostringstream string_stream;
   //string_stream << "check"
   //              << std::setfill('0')
   //              << std::setw(5)
   //              << PS::Comm::getRank()
   //              << ".dat";
   //std::string filename = string_stream.str();
   //std::ofstream output_file;
   //output_file.open(filename.c_str(),std::ios::trunc);
   //for (PS::S32 i=0; i<numPtclLocal; i++) {
   //   output_file << system[i].x.x << " "
   //               << system[i].x.y << " "
   //               << system[i].x.z << " "
   //               << system[i].pot << " "
   //               << std::endl;
   //}
   //output_file.close();
   //PS::Finalize();
   //std::exit(0);

}

