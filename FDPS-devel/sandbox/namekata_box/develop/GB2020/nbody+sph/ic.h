#pragma once
/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "user_defined.h"

void readTipsyFile(std::string& input_file_name,
                   PS::ParticleSystem<FP_nbody>& psys);

void GalaxyIC(PS::ParticleSystem<FP_nbody>& psys_nbody,
              PS::ParticleSystem<FP_star>& psys_star,
              PS::ParticleSystem<FP_gas>& psys_gas);
