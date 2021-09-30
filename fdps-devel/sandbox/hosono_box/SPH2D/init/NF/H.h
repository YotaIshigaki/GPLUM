#pragma once

#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
template <class Ptcl> class KHI: public Problem<Ptcl>{
	public:
	static const double END_TIME = 3.0;
	static void setupIC(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
		/////////
		//place ptcls
		/////////
		std::vector<Ptcl> ptcl;
		const PS::F64 dx = 1.0 / 512.0;
		const PS::F64 box_x = 1.5;
		const PS::F64 box_y = 1.0;
		PS::S32 cnt = 0;
		for(PS::F64 x = 0 ; x < box_x ; x += dx){
			for(PS::F64 y = 0 ; y < box_y ; y += dx){
				++cnt;
			}
		}
		const int NptclIn1Node = cnt / PS::Comm::getNumberOfProc();
		PS::S32 i = 0;
		for(PS::F64 x = 0 ; x < box_x ; x += dx){
			for(PS::F64 y = 0 ; y < box_y ; y += dx){
				Ptcl ith;
				ith.type = HYDRO;
				ith.pos.x = x;
				ith.pos.y = y;
				if(box_y / 3.0 > y){
					ith.vel.x = 0.5;
					ith.dens = 4.0;
				}else if(box_y * 2.0 / 3.0 > y){
					ith.vel.x = 0.0;
					ith.dens = 1.0;
				}else{
					ith.vel.x = -0.5;
					ith.dens = 2.0;
				}
				ith.mass = ith.dens;
				ith.eng  = 2.5 / (0.4 * ith.dens);
				ith.setPressure(&Diatomic);
				ith.vel.y = 0.1 * sin(4.0 * math::pi * x) * (exp(- powf((y/0.05), 2.0)) + exp(- powf(((y - box_y)/0.05), 2.0)) + exp(- powf((y - box_y / 3.0) / 0.05, 2.0)) + exp(- powf((y - box_y * 2.0 / 3.0) / 0.05, 2.0)));
				ith.id   = i++;
				if(ith.id / NptclIn1Node == PS::Comm::getRank()) ptcl.push_back(ith);
			}
		}
		for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
			ptcl[i].mass = ptcl[i].mass / (PS::F64)(cnt);
		}
		std::cout << "# of ptcls is... " << ptcl.size() << std::endl;
		//
		dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XY);
		dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0), PS::F64vec(box_x, box_y));
		const PS::S32 numPtclLocal = ptcl.size();
		sph_system.setNumberOfParticleLocal(numPtclLocal);
		for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
			sph_system[i] = ptcl[i];
		}
		/////////
		//Fin.
		std::cout << "setup..." << std::endl;
	}
	static void setEoS(PS::ParticleSystem<Ptcl>& sph_system){
		for(PS::U64 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].setPressure(&Diatomic);
		}
	}
};

#else
#error NO
#endif

