#pragma once
/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "user_defined.h"

void CoolingHeating(PS::ParticleSystem<FP_gas> & psys,
                    const PS::F64 dt_normalized);
void CalcEquilibriumState(const PS::F64 Zmetal,
                          const PS::F64 epsilon_FUV,
                          const PS::F64 G0_FUV);
