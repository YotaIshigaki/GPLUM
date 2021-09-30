#pragma once

#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
template <class Ptcl> class TEST: public Problem<Ptcl>{
	public:
	static const double END_TIME = 1.0;
	static const PS::F64 dx = 1.0 / 64.0;
	static const PS::F64 box_x = 1.0;
	static const PS::F64 box_y = 1.0;
	static const PS::F64 box_z = 1.0;
	static void setupIC(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
		/////////
		//place ptcls
		/////////
		std::vector<Ptcl> ptcl;
		PS::S32 cnt = 0;
		for(PS::F64 x = 0 ; x < box_x ; x += dx){
			for(PS::F64 y = 0 ; y < box_y ; y += dx){
				for(PS::F64 z = 0 ; z < box_y ; z += dx){
					++cnt;
				}
			}
		}
		const int NptclIn1Node = cnt / PS::Comm::getNumberOfProc();
		PS::S32 i = 0;
		for(PS::F64 x = 0 ; x < box_x ; x += dx){
			for(PS::F64 y = 0 ; y < box_y ; y += dx){
				for(PS::F64 z = 0 ; z < box_z ; z += dx){
					Ptcl ith;
					ith.type = HYDRO;
					ith.pos.x = x;
					ith.pos.y = y;
					ith.pos.z = z;
					ith.dens = 1.0;
					ith.mass = 1.0;
					ith.eng  = 2.5 / (0.4 * ith.dens);
					ith.setPressure(&Diatomic);
					ith.id   = i++;
					if(ith.id / NptclIn1Node == PS::Comm::getRank()) ptcl.push_back(ith);
				}
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
			sph_system[i].setPressure(&Diatomic);
		}
	}
	static void setDomain(PS::DomainInfo& dinfo){
		dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
		dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), PS::F64vec(box_x, box_y, box_z));
	}
};

#else
#error NO
#endif

