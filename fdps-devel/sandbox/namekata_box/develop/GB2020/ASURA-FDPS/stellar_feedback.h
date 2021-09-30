#pragma once
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "user_defined.h"

bool StellarFeedback(PS::ParticleSystem<FP_star> & psys_star,
                     PS::ParticleSystem<FP_gas> & psys_gas,
                     PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather & tree,
                     const PS::F64 time, const PS::F64 dt);
