#pragma once
/* C++ headers */
#include "common.h"
/* FDPS header */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "user_defined.h"

PS::F64 getTimeStep(const PS::ParticleSystem<FP_nbody>& psys_nbody,
                    const PS::ParticleSystem<FP_star>& psys_star,
                    const PS::ParticleSystem<FP_gas>& psys_gas);
