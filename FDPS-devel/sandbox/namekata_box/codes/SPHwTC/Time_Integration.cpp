/* Standard headers */
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "preprocess_keywords.h"
#include "Time_Integration.h"

namespace time_integrator {

/*-------------------------------------------------------------------*/
///////////////          S U B R O U T I N E            ///////////////
/////////////// < L E A P F R O G _ K I C K D R I F T > ///////////////
/*-------------------------------------------------------------------*/
void Leapfrog_KickDrift(PS::ParticleSystem<SPH_FP>& sph_system,
                        PS::DomainInfo& dinfo, const PS::F64 dt) {
   //* Local variables
   PS::F64 dt_half = 0.5e0*dt;
   PS::F64vec atot;
   PS::F64 udot;

   //* [1] Kick-Drift-Prediction
   for (PS::S32 i=0; i<sph_system.getNumberOfParticleLocal(); i++) {
      // Total acceleration
      atot = sph_system[i].ap
           + sph_system[i].avis; 
      udot = sph_system[i].gamad
           + sph_system[i].gamvis;
      // Error check
      if (isnan(atot.x) || isnan(atot.y) || 
          isnan(atot.z) || isnan(udot)   ||
          isnan(sph_system[i].cs)) {
         std::cout << "[LF1] NaN detected (i= " << i << ")" << std::endl;
         std::cout << "atot.x = " << atot.x << std::endl;
         std::cout << "atot.y = " << atot.y << std::endl;
         std::cout << "atot.z = " << atot.z << std::endl;
         std::cout << "udot   = " << atot   << std::endl;
         std::cout << "cs     = " << sph_system[i].cs << std::endl;
         std::exit(1);
      }
      sph_system[i].v_half = sph_system[i].v + atot * dt_half;
      sph_system[i].u_half = sph_system[i].u + udot * dt_half;
      sph_system[i].v += atot * dt;
      sph_system[i].u += udot * dt;
      sph_system[i].x += sph_system[i].v_half * dt;
      // where v_half and u_half are the velocity and the specific internal
      // energy at t^{n+1/2}, respectively. v and u are the predicted velocity
      // the predicted internal energy at t^{n+1}.
   }

   //* [2] Update domain info
   sph_system.adjustPositionIntoRootDomain(dinfo);
   dinfo.decomposeDomain();
   sph_system.exchangeParticle(dinfo);
}


/*-------------------------------------------------------------------*/
///////////////          S U B R O U T I N E            ///////////////
/////////////// < L E A P F R O G _ F I N A L K I C K > ///////////////
/*-------------------------------------------------------------------*/
void Leapfrog_FinalKick(PS::ParticleSystem<SPH_FP>& sph_system,
                        const PS::F64 dt){
   //* Local variables
   PS::F64 dt_half = 0.5e0*dt;
   PS::F64vec atot;
   PS::F64 udot;

   for (PS::S32 i=0; i<sph_system.getNumberOfParticleLocal(); i++){
      // Total acceleration
      atot = sph_system[i].ap 
           + sph_system[i].avis;
      udot = sph_system[i].gamad 
           + sph_system[i].gamvis;
      // Error check
      if (isnan(atot.x) || isnan(atot.y) || 
          isnan(atot.z) || isnan(udot)   ||
          isnan(sph_system[i].cs)) {
         std::cout << "[LF2] NaN detected (i= " << i << ")" << std::endl;
         std::cout << "atot.x = " << atot.x << std::endl;
         std::cout << "atot.y = " << atot.y << std::endl;
         std::cout << "atot.z = " << atot.z << std::endl;
         std::cout << "udot   = " << atot   << std::endl;
         std::cout << "cs     = " << sph_system[i].cs << std::endl;
         std::exit(1);
      }
      sph_system[i].v = sph_system[i].v_half + atot * dt_half;
      sph_system[i].u = sph_system[i].u_half + udot * dt_half;
   }

}

} // END of time_integrator

