/* Standard headers */
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cfloat>
/* Libraries */
//#ifndef __GNUC__
//#  define  __attribute__(x)  /*NOTHING*/
//#endif
#include <Eigen/Dense>
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#include <omp.h>
#endif
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "preprocess_keywords.h"
#include "Math_Functions.hpp"
#include "Kernel_Functions.h"
#include "Particle_Class.h"
#include "Calculation_Conditions.h"

//* Constants
#if (Simulation_Dimension == 1)
static const PS::S32 dimSpace=1;
#elif (Simulation_Dimension == 2)
static const PS::S32 dimSpace=2;
#elif (Simulation_Dimension == 3)
static const PS::S32 dimSpace=3;
#endif

//======================================
//* Class Defs.: Nbody_FP 
//======================================
// Constructor & destructor
Nbody_FP::Nbody_FP() {};
Nbody_FP::~Nbody_FP() {};
// Member functions required by FDPS
PS::F64 Nbody_FP::getCharge() const{
   return m;
}

PS::F64 Nbody_FP::getChargeParticleMesh() const{
   return m;
}

PS::F64vec Nbody_FP::getPos() const{
   return x;
}

void Nbody_FP::setPos(const PS::F64vec& x){
   this->x = x;
}

void Nbody_FP::copyFromForceParticleMesh(const PS::F64 apm){}
// The other member functions 

//=========================================
//* Class Defs.: Nbody_EP 
//=========================================
// Constructor & destructor
Nbody_EP::Nbody_EP() {};
Nbody_EP::~Nbody_EP() {};
// Member functions required by FDPS
PS::F64vec Nbody_EP::getPos() const{
   return x;
}
void Nbody_EP::setPos(const PS::F64vec& x){
   this->x = x;
}
void Nbody_EP::copyFromFP(const Nbody_FP& FP){
   id  = FP.id;
   m   = FP.m;
   eps = FP.eps;
   x   = FP.x;
}

//=========================================
//* Class Defs.: Nbody_IO
//=========================================
void Nbody_IO::copyFromFP(const Nbody_FP& FP) {
   id     = FP.id;
   m      = FP.m;
   eps    = FP.eps;
   x      = FP.x;
   v      = FP.v;
   agrv   = FP.agrv;
   pot    = FP.pot;
}

void Nbody_IO::copyToFP(Nbody_FP& FP) const {
   FP.id     = id;
   FP.m      = m;
   FP.eps    = eps;
   FP.x      = x;
   FP.v      = v;
   FP.agrv   = agrv;
   FP.pot    = pot;
}


//==================================
//* Class Defs.: Calc_gravity
//==================================
//void Calc_gravity::operator ()(const Nbody_EP* const ep_i,
//                               const PS::S32 Nip,
//                               const Nbody_EP* const ep_j,
//                               const PS::S32 Njp,
//                               Nbody_Symmetry_Results* const result){
//}

