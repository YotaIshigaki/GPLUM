#pragma once


#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
template <class Ptcl> class RTI: public Problem<Ptcl>{
	static double Density(const double y){
		return 1.5 - 1.0 * y;
	}
	public:
	static const double END_TIME = 5.0;
	static void setupIC(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
		/////////
		//place ptcls
		/////////
		std::vector<Ptcl> ptcl;
		const PS::F64 dx = 1.0 / 128;
		const PS::F64 box_x = 0.5;
		const PS::F64 box_y = 1.0;
		PS::S32 i = 0;
		for(PS::F64 x = 0 ; x < box_x ; x += dx){
			for(PS::F64 y = 0 ; y < box_y ; y += dx / Density(y)){
				std::cout << Density(y) << std::endl;
				Ptcl ith;
				ith.type = HYDRO;
				//if(y < 0.5){
					ith.pos.x = x;
					ith.pos.y = y;
					ith.dens = Density(y);
					ith.mass = 1.0;
					ith.eng  = 2.5 / (0.4);
					ith.id   = i++;
				//}else{
				//	ith.pos.x = x;
				//	ith.pos.y = y;
				//	ith.dens = 2.0;
				//	ith.mass = 1.0;
				//	ith.eng  = 2.5 / (0.4 * ith.dens);
				//	ith.id   = i++;
				//}
				ith.setPressure(&Diatomic);
				if(y > 0.9 || y < 0.1){
					ith.type = FREEZE;
				}
				ptcl.push_back(ith);
			}
		}
		for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
			ptcl[i].mass = 0.5 * (Density(0.0) + Density(1.0)) * box_x * box_y / (PS::F64)(ptcl.size());
		}
		std::cout << "# of ptcls is... " << ptcl.size() << std::endl;
		//
		dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_X);
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
		//Fin.
		std::cout << "setup..." << std::endl;
	}
	static void setEoS(PS::ParticleSystem<Ptcl>& sph_system){
		for(PS::U64 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].setPressure(&Diatomic);
		}
	}
	static void addExternalForce(PS::ParticleSystem<Ptcl>& sph_system, system_t& system){
		for(PS::U64 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].acc.y -= 1.0;
		}
	}
};
#else
#error NO
#endif

