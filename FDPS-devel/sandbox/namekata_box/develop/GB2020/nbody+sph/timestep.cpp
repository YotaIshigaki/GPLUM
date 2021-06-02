/* C++ header */
#include "common.h"
/* FDPS header */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "debug_utilities.h"
#include "run_parameters.h"
#include "user_defined.h"

PS::F64 getTimeStep(const PS::ParticleSystem<FP_nbody>& psys_nbody,
                    const PS::ParticleSystem<FP_star>& psys_star,
                    const PS::ParticleSystem<FP_gas>& psys_gas) {
    PS::F64 dt = DBL_MAX; 
    if (run_param::sim::dt_max > 0.0) dt = run_param::sim::dt_max;

    // Timescale for N-body system
    for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
        const PS::F64 acc = std::sqrt(psys_nbody[i].acc * psys_nbody[i].acc);
        if (acc > 0.0)
            dt = std::min(dt, run_param::sim::CFL_dyn * std::sqrt(run_param::sim::eps_grav / acc));
    }

   // Timescale for stellar system
    for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
        const PS::F64 acc = std::sqrt(psys_star[i].acc * psys_star[i].acc);
        if (acc > 0.0)
            dt = std::min(dt, run_param::sim::CFL_dyn * std::sqrt(run_param::sim::eps_grav / acc));
    }

   // Timescale for SPH system
   for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
      const PS::F64 acc = std::sqrt((psys_gas[i].acc_grav + psys_gas[i].acc_hydro)
                                   *(psys_gas[i].acc_grav + psys_gas[i].acc_hydro));
      if (acc > 0.0)
          dt = std::min(dt, run_param::sim::CFL_dyn * std::sqrt(run_param::sim::eps_grav / acc));
      dt = std::min(dt, psys_gas[i].dt);
   }
   return PS::Comm::getMinValue(dt);
}

