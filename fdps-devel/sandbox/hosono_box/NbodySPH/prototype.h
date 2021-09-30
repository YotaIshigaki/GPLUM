#pragma once

void DisplayInfo(void);
void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, system_t&, PS::DomainInfo& dinfo);

void InitialKick(PS::ParticleSystem<RealPtcl>& sph_system, const system_t sysinfo);
void FullDrift(PS::ParticleSystem<RealPtcl>& sph_system, const system_t sysinfo);
void Predict(PS::ParticleSystem<RealPtcl>& sph_system, const system_t sysinfo);
void FinalKick(PS::ParticleSystem<RealPtcl>& sph_system, const system_t sysinfo);

void Initialize(PS::ParticleSystem<RealPtcl>& sph_system);
void CalcPressure(PS::ParticleSystem<RealPtcl>& sph_system);

PS::F64 getTimeStepGlobal(const PS::ParticleSystem<RealPtcl>& sph_system);

void CheckConservativeVariables(const PS::ParticleSystem<RealPtcl>& sph_system);
void CalcExternalForce(PS::ParticleSystem<RealPtcl>& sph_system, const system_t sysinfo);
