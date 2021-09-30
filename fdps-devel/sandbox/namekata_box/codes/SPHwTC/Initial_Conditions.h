#pragma once
#include <particle_simulator.hpp>
#include "Particle_Class.h"
#include "Misc_Class.h"

extern void PM_test_IC(PS::ParticleSystem<SPH_FP>& sph_system,
                       PS::DomainInfo& dinfo);
extern void make_glass_IC(PS::ParticleSystem<SPH_FP>& sph_system,
                          PS::DomainInfo& dinfo,
                          Fluctuation_Monitor& f_monitor);
extern void shock_tube_IC(PS::ParticleSystem<SPH_FP>& sph_system,
                          PS::DomainInfo& dinfo);
