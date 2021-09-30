#pragma once
#include <particle_simulator.hpp>
#include "Particle_Class.h"
#include "Misc_Class.h"

extern void NaCl_IC(PS::ParticleSystem<Nbody_FP>& nbody_system,
                    PS::DomainInfo& dinfo,
                    Crystal_Parameters& params);
extern void PM_test_IC(PS::ParticleSystem<Nbody_FP>& nbody_system,
                       PS::DomainInfo& dinfo,
                       Crystal_Parameters& params);
