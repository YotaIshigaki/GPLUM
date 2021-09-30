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
#include "Particle_Class.h"
#include "Calculation_Conditions.h"

PS::F64 S2_pcut(const PS::F64 xi);
PS::F64 S2_fcut(const PS::F64 xi);

//======================================
//* Class Defs.: Nbody_PP_Results
//======================================
void Nbody_PP_Results::clear() {
   pot  = 0.0;
   agrv = 0.0;
}

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

PS::F64 Nbody_FP::getRSearch() const {
   return rc;
}

void Nbody_FP::setPos(const PS::F64vec& x){
   this->x = x;
}

void Nbody_FP::copyFromForce(const Nbody_PP_Results& result){}

void Nbody_FP::copyFromForceParticleMesh(const PS::F64 apm){}
// The other member functions 

//=========================================
//* Class Defs.: Nbody_EP 
//=========================================
// Constructor & destructor
Nbody_EP::Nbody_EP() {};
Nbody_EP::~Nbody_EP() {};
// Member functions required by FDPS
PS::F64 Nbody_EP::getCharge() const {
   return m;
}

PS::F64vec Nbody_EP::getPos() const{
   return x;
}

PS::F64 Nbody_EP::getRSearch() const{
   return rc;
}

void Nbody_EP::setPos(const PS::F64vec& x){
   this->x = x;
}

void Nbody_EP::copyFromFP(const Nbody_FP& FP){
   id  = FP.id;
   m   = FP.m;
   eps = FP.eps;
   rc  = FP.rc;
   x   = FP.x;
}

//==================================
//* Class Defs.: Calc_gravity
//==================================
void Calc_gravity::operator ()(const Nbody_EP* const ep_i,
                               const PS::S32 Nip,
                               const Nbody_EP* const ep_j,
                               const PS::S32 Njp,
                               Nbody_PP_Results* const result){

   for (PS::S32 i=0; i<Nip; i++) {
      for (PS::S32 j=0; j<Njp; j++) {
         PS::F64vec dx = ep_i[i].x - ep_j[j].x;
         PS::F64 rij = std::sqrt(dx * dx);
         if ((ep_i[i].id == ep_j[j].id) && (rij == 0.0)) continue;
         PS::F64 rinv = 1.0/rij;
         PS::F64 rinv3 = rinv*rinv*rinv;
         PS::F64 xi = 2.0*rij/ep_i[i].rc;
         result[i].pot  -= ep_j[j].m * S2_pcut(xi) * rinv;
         result[i].agrv -= ep_j[j].m * S2_fcut(xi) * rinv3 * dx;
      }
      //* Self-interaction term
#if (Self_Interaction_Term == 1)
      result[i].pot += ep_i[i].m * (208.0/(70.0*ep_i[i].rc));
      // Note that += instead of -= (remind 1-S2_pcut(xi) for xi=0)
#endif
   }

}

/*--------------------------------------------------------------------*/
/////////////////////////  F U N C T I O N  ////////////////////////////
///////////////////////// < S 2 _ P C U T > ////////////////////////////
/*--------------------------------------------------------------------*/
PS::F64 S2_pcut(const PS::F64 xi) {
   // This is the potential cutoff function where we used Eq.(8.75)
   // in Hockney & Eastwood (1987).

   if (xi <= 1.0) {
      return 1.0 - xi*(208.0 
                      +(xi*xi)*(-112.0 
                               +(xi*xi)*(56.0 
                                        +xi*(-14.0 
                                            +xi*(-8.0
                                                +3.0*xi)))))/140.0;
   } else if ((1.0 < xi) && (xi < 2.0)) {
      return 1.0 - (12.0 
                   +xi*(128.0
                       +xi*(224.0 
                           +xi*(-448.0 
                               +xi*(280.0 
                                   +xi*(-56.0 
                                       +xi*(-14.0 
                                           +xi*(8.0
                                               -xi))))))))/140.0;
   } else {
      return 0.0;
   }
}

/*--------------------------------------------------------------------*/
/////////////////////////  F U N C T I O N  ////////////////////////////
///////////////////////// < S 2 _ F C U T > ////////////////////////////
/*--------------------------------------------------------------------*/
PS::F64 S2_fcut(const PS::F64 xi) {
   // This function returns 1 - R(\xi), where \xi is r/(a/2), a is the
   // scale length of the cutoff function, and R(\xi) is almost the same
   // as the function defined as Eq.(8-72) in Hockney & Eastwood (1987).
   // The only difference is that [1/(r/(a/2))]^2 is factored out 
   // in this function from Eq.(8-72).

   if (xi <= 1.0) {
      return 1.0 - (xi*xi*xi)*(224.0 
                              +(xi*xi)*(-224.0
                                       +xi*(70.0
                                           +xi*(48.0-21.0*xi))))/140.0;
   } else if ((1.0 < xi) && (xi < 2.0)) {
      return 1.0 - (12.0 
                   +(xi*xi)*(-224.0 
                            +xi*(896.0 
                                +xi*(-840.0 
                                    +xi*(224.0 
                                        +xi*(70.0 
                                            +xi*(-48.0+7.0*xi)))))))/140.0;
   } else {
      return 0.0;
   }
}

