#pragma once
template <class ThisPtcl> PS::F64 getTimeStepGlobal(const PS::ParticleSystem<ThisPtcl>& sph_system){
	PS::F64 dt = 1.0e+30;//set VERY LARGE VALUE
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		PS::F64 dt_tmp = 1.0e+30;
		#if 1
		const PS::F64 acc = sqrt((sph_system[i].grav + sph_system[i].acc) * (sph_system[i].grav + sph_system[i].acc));
		dt_tmp = std::min(dt_tmp, sqrt(sph_system[i].smth / (acc + 1.0e-16)));
		//dt_tmp = std::min(dt_tmp, std::abs(sph_system[i].eng) / std::abs(sph_system[i].eng_dot + 1.0e-16));
		//dt_tmp = std::min(dt_tmp, 1.0 / std::abs(sph_system[i].div_v + 1.0e-16));
		dt = std::min(dt, std::min(sph_system[i].dt, dt_tmp));
		#else
		dt_tmp = std::max(sph_system[i].dt, PARAM::C_CFL * sph_system[i].smth / (sph_system[i].snds + 1.2 * (sph_system[i].snds + 2.0 * sph_system[i].smth * std::min(sph_system[i].div_v, 0.0))));
		dt = std::min(dt, dt_tmp);
		#endif
	}
	return PS::Comm::getMinValue(dt);
}

