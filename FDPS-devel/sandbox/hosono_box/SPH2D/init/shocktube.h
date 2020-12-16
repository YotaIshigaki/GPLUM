#pragma once


#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
template <class PTCL> void SetupIC(PS::ParticleSystem<PTCL>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
	/////////
	//place ptcls
	/////////
	std::vector<PTCL> ptcl;
	const PS::F64 dx = 1.0 / 256;
	const PS::F64 box_x = 1.0;
	const PS::F64 box_y = box_x / 8.0;
	PS::S32 i = 0;
	for(PS::F64 x = 0 ; x < box_x * 0.5 ; x += dx){
		for(PS::F64 y = 0 ; y < box_y ; y += dx){
			PTCL ith;
			ith.pos.x = x;
			ith.pos.y = y;
			ith.dens = 1.0;
			ith.mass = 0.75;
			ith.eng  = 2.5;
			ith.id   = i++;
			ith.setPressure(&Diatomic);
			ptcl.push_back(ith);
		}
	}
	for(PS::F64 x = box_x * 0.5 ; x < box_x * 1.0 ; x += dx * 2.0){
		for(PS::F64 y = 0 ; y < box_y ; y += dx){
			PTCL ith;
			ith.pos.x = x;
			ith.pos.y = y;
			ith.dens = 0.5;
			ith.mass = 0.75;
			ith.eng  = 1.0;
			ith.id   = i++;
			ith.setPressure(&Diatomic);
			ptcl.push_back(ith);
		}
	}
	for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
		ptcl[i].mass = ptcl[i].mass * box_x * box_y / (PS::F64)(ptcl.size());
	}
	std::cout << "# of ptcls is... " << ptcl.size() << std::endl;
	//
	dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XY);
	dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0), PS::F64vec(box_x, box_y));
	if(PS::Comm::getRank() == 0){
		const PS::S32 numPtclLocal = ptcl.size();
		sph_system.setNumberOfParticleLocal(numPtclLocal);
		for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
			sph_system[i] = ptcl[i];
		}
	}else{
		sph_system.setNumberOfParticleLocal(0);
	}
	/////////
	sysinfo.end_time = 0.2;
	//Fin.
	std::cout << "setup..." << std::endl;
}

#else
void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
	/////////
	//place ptcls
	/////////
	std::vector<RealPtcl> ptcl;
	const PS::F64 dx = 1.0 / 256;
	const PS::F64 box_x = 1.0;
	const PS::F64 box_y = box_x / 8.0;
	const PS::F64 box_z = box_x / 8.0;
	PS::S32 i = 0;
	for(PS::F64 x = 0 ; x < box_x * 0.5 ; x += dx){
		for(PS::F64 y = 0 ; y < box_y ; y += dx){
			for(PS::F64 z = 0 ; z < box_z ; z += dx){
				RealPtcl ith;
				ith.pos.x = x;
				ith.pos.y = y;
				ith.pos.z = z;
				ith.dens = 1.0;
				ith.mass = 0.75;
				ith.eng  = 2.5;
				ith.id   = i++;
				ith.mat = EoS::Diatomic;
				ptcl.push_back(ith);
			}
		}
	}
	for(PS::F64 x = box_x * 0.5 ; x < box_x * 1.0 ; x += dx * 2.0){
		for(PS::F64 y = 0 ; y < box_y ; y += dx){
			for(PS::F64 z = 0 ; z < box_z ; z += dx){
				RealPtcl ith;
				ith.pos.x = x;
				ith.pos.y = y;
				ith.pos.z = z;
				ith.dens = 0.5;
				ith.mass = 0.75;
				ith.eng  = 1.0;
				ith.id   = i++;
				ith.mat = EoS::Diatomic;
				ptcl.push_back(ith);
			}
		}
	}
	for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
		ptcl[i].mass = ptcl[i].mass * box_x * box_y * box_z / (PS::F64)(ptcl.size());
	}
	std::cout << "# of ptcls is... " << ptcl.size() << std::endl;
	//
	dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
	dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), PS::F64vec(box_x, box_y, box_z));
	if(PS::Comm::getRank() == 0){
		const PS::S32 numPtclLocal = ptcl.size();
		sph_system.setNumberOfParticleLocal(numPtclLocal);
		for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
			sph_system[i] = ptcl[i];
		}
	}else{
		sph_system.setNumberOfParticleLocal(0);
	}
	/////////
	sysinfo.end_time = 0.2;
	//Fin.
	std::cout << "setup..." << std::endl;
}
#endif

template <class ThisPtcl> void SetEoS(PS::ParticleSystem<ThisPtcl>& sph_system){
	for(PS::U64 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].setPressure(&Diatomic);
	}
}

