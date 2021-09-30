#include "header.h"

void Initialize(PS::ParticleSystem<RealPtcl>& sph_system){
	#pragma omp parallel for
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].smth = PARAM::SMTH * pow(sph_system[i].mass / sph_system[i].dens, 1.0/(PS::F64)(PARAM::Dim));
	}
}
