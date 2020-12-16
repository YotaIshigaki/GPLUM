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

/*
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
*/

#ifdef INITIAL_CONDITION
void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
	/////////
	//place ptcls
	/////////
	std::vector<RealPtcl> ptcl;
	std::vector<RealPtcl> tar;//Target
	std::vector<RealPtcl> imp;//Impactor
	/////////
	const PS::F64 EarthMass = 6.0e+24;
	const PS::F64 EarthRadi = 6400e+3;
	const PS::F64 tarMass = EarthMass;
	const PS::F64 tarRadi = EarthRadi;
	const PS::F64 impMass = 0.1 * EarthMass;
	const PS::F64 impRadi = 0.5 * EarthRadi;
	sysinfo.Grav = 6.67e-11;
	sysinfo.end_time = 1.0e+4;

	const PS::F64 dx = 1.0 / 128.0;
	PS::S32 id = 0;

	//Put Tar.
	for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
		for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
			for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
				const PS::F64 r = sqrt(x*x + y*y + z*z);
				if(r > 1.0) continue;
				RealPtcl ith;
				ith.pos.x = tarRadi * x;
				ith.pos.y = tarRadi * y;
				ith.pos.z = tarRadi * z;
				ith.dens = tarMass / (4.0 / 3.0 * math::pi * tarRadi * tarRadi * tarRadi);
				ith.mass = tarMass + impMass;
				ith.eng  = 0.1 * sysinfo.Grav * tarMass / tarRadi;
				ith.id   = id++;
				ith.mat = EoS::Granite;
				ith.tag = 0;
				tar.push_back(ith);
			}
		}
	}
	//Put Imp.
	for(PS::F64 x = -1.0 ; x <= 1.0 ; x += 2.0 * dx){
		for(PS::F64 y = -1.0 ; y <= 1.0 ; y += 2.0 * dx){
			for(PS::F64 z = -1.0 ; z <= 1.0 ; z += 2.0 * dx){
				const PS::F64 r = sqrt(x*x + y*y + z*z);
				if(r > 1.0) continue;
				RealPtcl ith;
				ith.pos.x = impRadi * x + 5.0 * tarRadi;
				ith.pos.y = impRadi * y;
				ith.pos.z = impRadi * z;
				ith.dens = impMass / (4.0 / 3.0 * math::pi * impRadi * impRadi * impRadi);
				ith.mass = tarMass + impMass;
				ith.eng  = 0.1 * sysinfo.Grav * impMass / impRadi;
				ith.id   = id++;
				ith.mat = EoS::Granite;
				ith.tag = 1;
				imp.push_back(ith);
			}
		}
	}
	for(PS::U32 i = 0 ; i < tar.size() ; ++ i){
		tar[i].mass /= (PS::F64)(tar.size() + imp.size());
	}
	for(PS::U32 i = 0 ; i < imp.size() ; ++ i){
		imp[i].mass /= (PS::F64)(tar.size() + imp.size());
	}

	for(PS::U32 i = 0 ; i < tar.size() ; ++ i){
		ptcl.push_back(tar[i]);
	}	
	for(PS::U32 i = 0 ; i < imp.size() ; ++ i){
		ptcl.push_back(imp[i]);
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
	std::cout << "# of ptcls = " << ptcl.size() << std::endl;
	std::cout << "setup..." << std::endl;
}
#else

void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
	FileHeader header;
	sph_system.readParticleAscii("init.dat", header);
	//Fin.
	sysinfo.Grav = 6.67e-11;
	sysinfo.end_time = 1.0e+5;

	//Shift positions
	std::vector<RealPtcl> tar, imp;
	PS::F64vec pos_tar = 0;
	PS::F64vec pos_imp = 0;
	PS::F64 mass_tar = 0;
	PS::F64 mass_imp = 0;
	PS::F64 radi_tar = 0;
	PS::F64 radi_imp = 0;
	for(PS::U32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		if(sph_system[i].tag == 0){
			//target
			pos_tar  += sph_system[i].mass * sph_system[i].pos;
			mass_tar += sph_system[i].mass;
		}else{
			//impactor
			pos_imp  += sph_system[i].mass * sph_system[i].pos;
			mass_imp += sph_system[i].mass;
		}
	}
	//accumurate
	pos_tar  = PS::Comm::getSum( pos_tar);
	pos_imp  = PS::Comm::getSum( pos_imp);
	mass_tar = PS::Comm::getSum(mass_tar);
	mass_imp = PS::Comm::getSum(mass_imp);
	pos_tar /= mass_tar;
	pos_imp /= mass_imp;
	for(PS::U32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		if(sph_system[i].tag == 0){
			//target
			radi_tar = std::max(radi_tar, sqrt((pos_tar - sph_system[i].pos) * (pos_tar - sph_system[i].pos)));
		}else{
			//impactor
			radi_imp = std::max(radi_imp, sqrt((pos_imp - sph_system[i].pos) * (pos_imp - sph_system[i].pos)));
		}
	}
	radi_tar = PS::Comm::getMaxValue(radi_tar);
	radi_imp = PS::Comm::getMaxValue(radi_imp);
	std::cout << radi_tar << std::endl;
	std::cout << radi_imp << std::endl;

	const double L_EM = 3.5e+34;
	const double v_esc = sqrt(2.0 * sysinfo.Grav * (mass_tar + mass_imp) / (radi_tar + radi_imp));
	const double v_inf = 0.0;
	const double L_init = L_EM * 1.2;
	const double x_init = 3.0 * radi_tar;

	double y_init = radi_tar;//Initial guess.
	double v_init;
	std::cout << "v_esc = " << v_esc << std::endl;
	for(int it = 0 ; it < 10 ; ++ it){
		v_init = sqrt(v_inf * v_inf + 2.0 * sysinfo.Grav * mass_tar / sqrt(x_init * x_init + y_init * y_init));
		y_init = L_init / (mass_imp * v_init);
	}
	const double v_imp = sqrt(v_inf * v_inf + 2.0 * sysinfo.Grav * (mass_tar + mass_imp) / (radi_tar + radi_imp));
	std::cout << "v_init = "  << v_init <<  std::endl;
	std::cout << "y_init / Rtar = "  << y_init / radi_tar <<  std::endl;
	std::cout << "v_imp  = "  << v_imp <<  std::endl;
	//shift'em

	for(PS::U32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		if(sph_system[i].tag == 0){
			//Do nothing
		}else{
			sph_system[i].pos = sph_system[i].pos - pos_imp;
			sph_system[i].pos.x += x_init;
			sph_system[i].pos.y += y_init;
			sph_system[i].vel.x -= v_init;
		}
	}
	std::cout << "setup..." << std::endl;
}
#endif

void Initialize(PS::ParticleSystem<RealPtcl>& sph_system){
	#pragma omp parallel for
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].smth = PARAM::SMTH * pow(sph_system[i].mass / sph_system[i].dens, 1.0/(PS::F64)(PARAM::Dim));
	}
}





