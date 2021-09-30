#pragma once

void DisplayInfo(void);
void SetupIC(PS::ParticleSystem<RealPtcl>&, system_t&);
void SetupICWithCMB(PS::ParticleSystem<RealPtcl>&, system_t&);

void InitialKick(PS::ParticleSystem<RealPtcl>&, const system_t);
void FullDrift(PS::ParticleSystem<RealPtcl>&, const system_t);
void Predict(PS::ParticleSystem<RealPtcl>&, const system_t);
void FinalKick(PS::ParticleSystem<RealPtcl>&, const system_t);

void Initialize(PS::ParticleSystem<RealPtcl>&);
void CalcPressure(PS::ParticleSystem<RealPtcl>&);
void CalcDensityAndSmoothing(PS::ParticleSystem<RealPtcl>&, PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather&, PS::DomainInfo&);

void ShiftOrigin(PS::ParticleSystem<RealPtcl>&);

PS::F64 getTimeStepGlobal(const PS::ParticleSystem<RealPtcl>&);

void CheckConservativeVariables(const PS::ParticleSystem<RealPtcl>&, system_t);
void CalcExternalForce(PS::ParticleSystem<RealPtcl>&, const system_t);

void OutputFile(PS::ParticleSystem<RealPtcl>&, const system_t&);
void OutputFileWithTimeInterval(PS::ParticleSystem<RealPtcl>&, const system_t&);
void InputFile (PS::ParticleSystem<RealPtcl>&, system_t&);
void InputFileWithTimeInterval(PS::ParticleSystem<RealPtcl>&, system_t&);
