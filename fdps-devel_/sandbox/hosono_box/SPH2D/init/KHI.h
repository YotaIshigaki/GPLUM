#pragma once

#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
template <class Ptcl> class KHI: public Problem<Ptcl>{
	public:
	static const double END_TIME = 3.0;
	static const PS::F64 dx = 1.0 / 256.0;
	static const PS::F64 box_x = 1.0;
	static const PS::F64 box_y = 1.0;
	static void setupIC(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
		/////////
		//place ptcls
		/////////
		std::vector<Ptcl> ptcl;
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
				if(0.25 < y && y < 0.75){
					ith.vel.x = 0.5;
					ith.dens = 2.0;
					ith.mass = 2.0;
				}else{
					ith.vel.x = -0.5;
					ith.dens = 1.0;
					ith.mass = 1.0;
				}
				ith.vel.y = (fabs(ith.pos.y - 0.25) < 0.025 || fabs(ith.pos.y - 0.75) < 0.025) ? 0.025 * sin(2.0 * M_PI * (ith.pos.x + 0.5) * 6) : 0;
				ith.eng  = 2.5 / (ith.dens * (2./3.));
				ith.setPressure(&Monoatomic);
				//ith.vel.y = 0.1 * sin(4.0 * math::pi * x) * (exp(- (y - 0.25) * (y - 0.25) / (0.05 * 0.05)) + exp(- (y - 0.75) * (y - 0.75) / (0.05 * 0.05)));
				//ith.eng  = 2.5 / (0.4 * ith.dens);
				//ith.setPressure(&Diatomic);
				ith.id   = i++;
				if(ith.id / NptclIn1Node == PS::Comm::getRank()) ptcl.push_back(ith);
			}
		}
		for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
			ptcl[i].mass = ptcl[i].mass / (PS::F64)(cnt);
		}
		std::cout << "# of ptcls is... " << ptcl.size() << std::endl;
		//
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
			sph_system[i].setPressure(&Monoatomic);
		}
	}
	static void setDomain(PS::DomainInfo& dinfo){
		dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XY);
		dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0), PS::F64vec(box_x, box_y));
	}
};

#else
#error NO
#endif

