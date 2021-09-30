#pragma once
/* C++ headers */
#include "common.h"
/* FDPS header */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "user_defined.h"

void calcDensity(PS::ParticleSystem<FP_gas> & psys_gas,
                 PS::ParticleSystem<FP_star> & psys_star,
                 PS::DomainInfo & dinfo, 
                 PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather & tree,
                 const PS::F64 t, const PS::F64 dt);

void setEntropy(PS::ParticleSystem<FP_gas>& psys);

void setPressure(PS::ParticleSystem<FP_gas>& psys);
