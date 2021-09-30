#pragma once
template <class ThisPtcl> void CheckConservativeVariables(const PS::ParticleSystem<ThisPtcl>& sph_system){
	PS::F64vec Mom = 0;//total momentum
	PS::F64    Eng = 0;//total enegry
	PS::F64    Mass = 0;//total enegry
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		Mom += sph_system[i].vel * sph_system[i].mass;
		Eng += (sph_system[i].eng + 0.5 * sph_system[i].vel * sph_system[i].vel) * sph_system[i].mass;
		Mass += sph_system[i].mass;
	}
	Eng  = PS::Comm::getSum(Eng);
	Mom  = PS::Comm::getSum(Mom);
	Mass = PS::Comm::getSum(Mass);
	printf("%.16e\n", Mass);
	printf("%.16e\n", Eng);
	printf("%.16e\n", Mom.x);
	printf("%.16e\n", Mom.y);
	#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
	printf("%.16e\n", Mom.z);
	#endif
}

