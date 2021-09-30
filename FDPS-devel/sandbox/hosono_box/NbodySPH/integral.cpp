#include "header.h"

PS::F64 getTimeStepGlobal(const PS::ParticleSystem<RealPtcl>& sph_system){
	PS::F64 dt = 1.0e+30;//set VERY LARGE VALUE
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		dt = std::min(dt, sph_system[i].dt);
	}
	//return 0.01;
	return PS::Comm::getMinValue(dt);
}

void InitialKick(PS::ParticleSystem<RealPtcl>& sph_system, const system_t sysinfo){
	#pragma omp parallel for
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].vel_half = sph_system[i].vel + 0.5 * sysinfo.dt * (sph_system[i].acc + sysinfo.Grav * sph_system[i].grav);
		sph_system[i].eng_half = sph_system[i].eng + 0.5 * sysinfo.dt * sph_system[i].eng_dot;
	}
}

void FullDrift(PS::ParticleSystem<RealPtcl>& sph_system, const system_t sysinfo){
	//time becomes t + dt;
	#pragma omp parallel for
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].pos += sysinfo.dt * sph_system[i].vel_half;
	}
}

void Predict(PS::ParticleSystem<RealPtcl>& sph_system, const system_t sysinfo){
	#pragma omp parallel for
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].vel += sysinfo.dt * (sph_system[i].acc + sysinfo.Grav * sph_system[i].grav);
		sph_system[i].eng += sysinfo.dt * sph_system[i].eng_dot;
	}
}

void FinalKick(PS::ParticleSystem<RealPtcl>& sph_system, const system_t sysinfo){
	#pragma omp parallel for
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].vel = sph_system[i].vel_half + 0.5 * sysinfo.dt * (sph_system[i].acc + sysinfo.Grav * sph_system[i].grav);
		sph_system[i].eng = sph_system[i].eng_half + 0.5 * sysinfo.dt * sph_system[i].eng_dot;
	}
}
