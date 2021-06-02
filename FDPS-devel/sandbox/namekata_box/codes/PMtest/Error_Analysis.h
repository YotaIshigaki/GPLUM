#pragma once
#include <particle_simulator.hpp>
#include "Particle_Class.h"

void error_analysis(PS::ParticleSystem<Nbody_FP>& system);

PS::F64 calc_NaCl_error(PS::ParticleSystem<Nbody_FP>& system,
                        PS::S32 nstep);
