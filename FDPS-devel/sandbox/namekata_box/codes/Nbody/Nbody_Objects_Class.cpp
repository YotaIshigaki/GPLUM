/* Standard headers */
#include <cmath>
#include <limits>
#include <vector>
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "preprocess_keywords.h"
#include "Particle_Class.h"
#include "Nbody_Objects_Class.h"
#include "Calculation_Conditions.h"

//================================
//* Class Defs.: Nbody_Objects
//================================
void Nbody_Objects::init_tree() {
   //PS::S32 numPtclLocal = system.getNumberOfParticleLocal();
   //PS::U64 n_glb_tot = 3 * numPtclLocal;
   //grv_tree.initialize(n_glb_tot);
}
void Nbody_Objects::calc_gravity() { }
PS::F64 Nbody_Objects::calc_timestep() { } 
