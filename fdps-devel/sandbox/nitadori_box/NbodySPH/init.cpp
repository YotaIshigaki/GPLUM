#include "header.h"

/*
void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
	/////////
	//place ptcls
	/////////
	std::vector<RealPtcl> ptcl;
	const PS::F64 dx = 1.0 / 128.0;
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
				ith.eng  = 2.5;
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
	sysinfo.end_time = 0.11;
	//Fin.
	std::cout << "setup..." << std::endl;
}
*/

#if 0
void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
	/////////
	//place ptcls
	/////////
	std::vector<RealPtcl> ptcl;
	const PS::S32 NumPtcl = 256;
	const PS::F64 box_x = 1.0;
	const PS::F64 box_y = 1.0;
	const PS::F64 box_z = 1.0;
	for(int i = 0 ; i < NumPtcl ; ++ i){
		const PS::F64 x = rand() / (double)(RAND_MAX);
		const PS::F64 y = rand() / (double)(RAND_MAX);
		const PS::F64 z = rand() / (double)(RAND_MAX);
		const PS::F64 u = rand() / (double)(RAND_MAX);
		const PS::F64 v = rand() / (double)(RAND_MAX);
		const PS::F64 w = rand() / (double)(RAND_MAX);
		RealPtcl ith;
		ith.pos.x = x;
		ith.pos.y = y;
		ith.pos.z = z;
		ith.vel.x = u;
		ith.vel.y = v;
		ith.vel.z = w;
		ith.dens = 1.0;
		ith.mass = 1.0;
		ith.eng  = 1.0;
		ith.id   = i;
		ith.mat = EoS::Diatomic;
		ptcl.push_back(ith);
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
	sysinfo.end_time = 0.0;
	//Fin.
	std::cout << "setup..." << std::endl;
}
#endif

#if 1
void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
	/////////
	//place ptcls
	/////////
	std::vector<RealPtcl> ptcl;
	/////////
	const PS::F64 Mass = 1.0;
	const PS::F64 Radi = 1.0;
	sysinfo.Grav = 1.0;
	sysinfo.end_time = 2.2;
	const PS::F64 dx = 1.0 / 16.0;
	PS::S32 id = 0;
	for(PS::F64 x = -Radi ; x <= Radi ; x += dx){
		for(PS::F64 y = -Radi ; y <= Radi ; y += dx){
			for(PS::F64 z = -Radi ; z <= Radi ; z += dx){
				const PS::F64 r = sqrt(x*x + y*y + z*z);
				if(r > Radi) continue;
				RealPtcl ith;
				ith.pos.x = x;
				ith.pos.y = y;
				ith.pos.z = z;
				ith.pos *= sqrt(r);

				ith.dens = (r <= 1.0e-10) ? 1.0 / dx : Mass / (2.0 * math::pi * Radi * Radi * r);
				ith.mass = Mass;
				ith.eng  = 0.05 * sysinfo.Grav * Mass / Radi;
				ith.id   = id++;
				ith.mat = EoS::Monoatomic;
				ptcl.push_back(ith);
			}
		}
	}
	for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
		ptcl[i].mass = ptcl[i].mass / (PS::F64)(ptcl.size());
	}

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
#endif

void Initialize(PS::ParticleSystem<RealPtcl>& sph_system){
	#pragma omp parallel for
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].smth = PARAM::SMTH * pow(sph_system[i].mass / sph_system[i].dens, 1.0/(PS::F64)(PARAM::Dim));
	}
}





