#pragma once
/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "user_defined.h"

void InitialKick(PS::ParticleSystem<FP_dm> & psys, const PS::F64 dt);
void InitialKick(PS::ParticleSystem<FP_star> & psys, const PS::F64 dt);
void InitialKick(PS::ParticleSystem<FP_gas> & psys, const PS::F64 dt);

void FullDrift(PS::ParticleSystem<FP_dm>& psys, const PS::F64 dt);
void FullDrift(PS::ParticleSystem<FP_star>& psys, const PS::F64 dt);
void FullDrift(PS::ParticleSystem<FP_gas>& psys, const PS::F64 dt);

void Predict(PS::ParticleSystem<FP_gas>& psys, const PS::F64 dt);

void FinalKick(PS::ParticleSystem<FP_dm>& psys, const PS::F64 dt);
void FinalKick(PS::ParticleSystem<FP_star>& psys, const PS::F64 dt);
void FinalKick(PS::ParticleSystem<FP_gas>& psys, const PS::F64 dt);
