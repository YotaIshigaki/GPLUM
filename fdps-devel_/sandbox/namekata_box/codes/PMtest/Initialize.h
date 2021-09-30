#pragma once
#include <particle_simulator.hpp>
#include "Particle_Class.h"
#include "Nbody_Objects_Class.h"
#include "Misc_Class.h"

extern void Initialize(Nbody_Objects& Nbody_objs,
                       Crystal_Parameters& NaCl_params,
                       PS::S32 nstep);
extern void Finalize();
