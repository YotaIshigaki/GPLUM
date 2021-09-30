/* C++ header */
#include "common.h"
/* FDPS header */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "debug_utilities.h"
#include "run_parameters.h"
#include "user_defined.h"

PS::F64 getTimeStep(const PS::ParticleSystem<FP_dm>& psys_dm,
                    const PS::ParticleSystem<FP_star>& psys_star,
                    const PS::ParticleSystem<FP_gas>& psys_gas) {
    PS::F64 dt = DBL_MAX; 
    if (run_param::grav::soft::dt_fid > 0.0) dt = run_param::grav::soft::dt_fid;

    // Timescale for N-body system
    for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
        const PS::F64 acc = std::sqrt(psys_dm[i].acc * psys_dm[i].acc);
        if (acc > 0.0)
            dt = std::min(dt, run_param::grav::soft::CFL * std::sqrt(run_param::grav::soft::eps_dm / acc));
    }

    // Timescale for stellar system
    for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
        const PS::F64 acc = std::sqrt(psys_star[i].acc * psys_star[i].acc);
        if (acc > 0.0)
            dt = std::min(dt, run_param::grav::soft::CFL * std::sqrt(run_param::grav::soft::eps_star / acc));
    }

    // Timescale for SPH system
    for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
        const PS::F64 acc = std::sqrt((psys_gas[i].acc_grav + psys_gas[i].acc_hydro)
                                     *(psys_gas[i].acc_grav + psys_gas[i].acc_hydro));
        if (acc > 0.0)
            dt = std::min(dt, run_param::grav::soft::CFL * std::sqrt(run_param::grav::soft::eps_gas / acc));
        dt = std::min(dt, psys_gas[i].dt);
    }
    return PS::Comm::getMinValue(dt);
}

