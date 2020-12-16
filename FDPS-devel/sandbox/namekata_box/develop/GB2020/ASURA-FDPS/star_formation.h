#pragma once
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "user_defined.h"

bool StarFormation(PS::ParticleSystem<FP_gas> & psys_gas,
                   PS::ParticleSystem<FP_star> & psys_star,
                   const PS::F64 t, const PS::F64 dt);
