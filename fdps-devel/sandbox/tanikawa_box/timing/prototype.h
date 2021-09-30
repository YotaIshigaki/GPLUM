#pragma once

void DisplayInfo(void);
void SetupIC(PS::ParticleSystem<RealPtcl>&, system_t&, PS::DomainInfo&);

void InitialKick(PS::ParticleSystem<RealPtcl>&, const system_t);
void FullDrift(PS::ParticleSystem<RealPtcl>&, const system_t);
void Predict(PS::ParticleSystem<RealPtcl>&, const system_t);
void FinalKick(PS::ParticleSystem<RealPtcl>&, const system_t);

void Initialize(PS::ParticleSystem<RealPtcl>&);
void CalcPressure(PS::ParticleSystem<RealPtcl>&);

void ShiftOrigin(PS::ParticleSystem<RealPtcl>&);

PS::F64 getTimeStepGlobal(const PS::ParticleSystem<RealPtcl>&);

void CheckConservativeVariables(const PS::ParticleSystem<RealPtcl>&);
void CalcExternalForce(PS::ParticleSystem<RealPtcl>&, const system_t);
