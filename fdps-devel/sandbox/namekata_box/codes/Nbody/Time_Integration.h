#pragma once
#include <particle_simulator.hpp>
#include "Particle_Class.h"

namespace time_integrator {

void Leapfrog_KickDrift(PS::ParticleSystem<Nbody_FP>& nbody_system,
                        PS::DomainInfo& dinfo, const PS::F64 dt);
void Leapfrog_FinalKick(PS::ParticleSystem<Nbody_FP>& nbody_system,
                        const PS::F64 dt);

} // END of time_integrator
