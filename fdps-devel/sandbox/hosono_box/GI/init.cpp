#include "header.h"

#ifdef INITIAL_CONDITION
#ifdef WITH_CORE
void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, system_t& sysinfo){
	/////////
	//place ptcls
	/////////
	std::vector<RealPtcl> ptcl;
	std::vector<RealPtcl> tar;//Target
	std::vector<RealPtcl> imp;//Impactor
	/////////
	const PS::F64 UnitMass = 6.0e+24; // Earth mass
	const PS::F64 UnitRadi = 6400e+3; // Earth radii
	//total-to-core frac.
	const PS::F64 coreFracRadi = 3500.0e+3 / 6400.0e+3;//Earth
	const PS::F64 coreFracMass = 0.3;//Earth
	/////////
	const PS::F64 tarMass = UnitMass;
	const PS::F64 tarRadi = UnitRadi;
	const PS::F64 tarCoreMass = tarMass * coreFracMass;
	const PS::F64 tarCoreRadi = tarRadi * coreFracRadi;
	const PS::F64 impMass = 0.1 * tarMass;
	const PS::F64 impFracRadi = cbrt(impMass / tarMass);
	const PS::F64 impRadi = impFracRadi * UnitRadi;
	const PS::F64 impCoreMass = impMass * coreFracMass;
	sysinfo.Grav = 6.67e-11;
	sysinfo.end_time = 1.0e+4;
	const double offset = 5.0 * UnitRadi;
	const PS::F64 dx = 1.0 / 512.0;
	///////////////////
	//Dummy put to determine # of ptcls
	///////////////////
	//target
	int tarNmntl = 0;
	for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
		for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
			for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
				const PS::F64 r = sqrt(x*x + y*y + z*z);
				if(r >= 1.0 || r <= coreFracRadi) continue;
				++ tarNmntl;
			}
		}
	}
	int tarNcore;
	double tarCoreShrinkFactor = 1.0;
	while(tarCoreShrinkFactor *= 0.99){
		tarNcore = 0;
		for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
			for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
				for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
					const PS::F64 r = tarCoreShrinkFactor * sqrt(x*x + y*y + z*z);
					if(r + dx >= coreFracRadi) continue;
					++ tarNcore;
				}
			}
		}
		if((double)(tarNcore) / (double)(tarNcore + tarNmntl) > coreFracMass) break;
	}
	//imp
	int impNmntl = 0;
	for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
		for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
			for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
				const PS::F64 r = sqrt(x*x + y*y + z*z);
				if(r >= impFracRadi || r <= impFracRadi * coreFracRadi) continue;
				++ impNmntl;
			}
		}
	}
	double impCoreShrinkFactor = 1.0;
	int impNcore;
	while(impCoreShrinkFactor *= 0.99){
		impNcore = 0;
		for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
			for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
				for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
					const PS::F64 r = impCoreShrinkFactor * sqrt(x*x + y*y + z*z);
					if(r + dx >= impFracRadi * coreFracRadi) continue;
					++ impNcore;
				}
			}
		}
		if((double)(impNcore) / (double)(impNcore + impNmntl) > coreFracMass) break;
	}
	///////////////////
	//Dummy end
	///////////////////
	const int tarNptcl = tarNcore + tarNmntl;
	const int impNptcl = impNcore + impNmntl;
	const int Nptcl    = tarNptcl + impNptcl;
	std::cout << "Target  :" << tarNptcl << std::endl;
	std::cout << "    total-to-core    : " << (double)(tarNcore) / (double)(tarNptcl) << std::endl;
	std::cout << "    # of core ptcls  : " << tarNcore << std::endl;
	std::cout << "    # of mantle ptcls: " << tarNmntl << std::endl;
	std::cout << "    total-to-core    : " << (double)(tarNcore) / (double)(tarNptcl) << std::endl;
	std::cout << "Impactor:" << impNptcl << std::endl;
	std::cout << "    total-to-core    : " << (double)(impNcore) / (double)(impNptcl) << std::endl;
	std::cout << "    # of core ptcls  : " << impNcore << std::endl;
	std::cout << "    # of mantle ptcls: " << impNmntl << std::endl;
	std::cout << "    total-to-core    : " << (double)(impNcore) / (double)(impNptcl) << std::endl;
	std::cout << "Total:" << Nptcl << std::endl;
	std::cout << "Tar-to-Imp mass ratio: " << (double)(impNmntl) / (double)(tarNmntl) << std::endl;
	//exit(0);
	const int NptclIn1Node = Nptcl / PS::Comm::getNumberOfProc();
	///////////////////
	//Real put
	///////////////////
	PS::S32 id = 0;
	//Put Tar.
	for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
		for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
			for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
				const PS::F64 r = sqrt(x*x + y*y + z*z);
				if(r >= 1.0 || r <= coreFracRadi) continue;
				RealPtcl ith;
				ith.pos.x = UnitRadi * x;
				ith.pos.y = UnitRadi * y;
				ith.pos.z = UnitRadi * z;
				ith.dens = tarMass / (4.0 / 3.0 * math::pi * tarRadi * tarRadi * tarRadi);
				ith.mass = tarMass + impMass;
				ith.eng  = 0.1 * sysinfo.Grav * tarMass / tarRadi;
				ith.id   = id++;
				ith.mat = EoS::Granite;
				ith.tag = 0;
				if(ith.id / NptclIn1Node == PS::Comm::getRank()) tar.push_back(ith);
			}
		}
	}
	for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
		for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
			for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
				const PS::F64 r = tarCoreShrinkFactor * sqrt(x*x + y*y + z*z);
				if(r + dx >= coreFracRadi) continue;
				RealPtcl ith;
				ith.pos.x = tarCoreShrinkFactor * UnitRadi * x;
				ith.pos.y = tarCoreShrinkFactor * UnitRadi * y;
				ith.pos.z = tarCoreShrinkFactor * UnitRadi * z;
				ith.dens = tarCoreMass / (4.0 / 3.0 * math::pi * tarCoreRadi * tarCoreRadi * tarCoreRadi);
				ith.mass = tarMass + impMass;
				ith.eng  = 0.1 * sysinfo.Grav * tarMass / tarRadi;
				ith.id   = id++;
				ith.mat = EoS::Iron;
				ith.tag = 0;
				if(ith.id / NptclIn1Node == PS::Comm::getRank()) tar.push_back(ith);
			}
		}
	}
	//imp
	for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
		for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
			for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
				const PS::F64 r = sqrt(x*x + y*y + z*z);
				if(r >= impFracRadi || r <= impFracRadi * coreFracRadi) continue;
				RealPtcl ith;
				ith.pos.x = UnitRadi * x + offset;
				ith.pos.y = UnitRadi * y;
				ith.pos.z = UnitRadi * z;
				ith.dens = impMass / (4.0 / 3.0 * math::pi * impRadi * impRadi * impRadi);
				ith.mass = tarMass + impMass;
				ith.eng  = 0.1 * sysinfo.Grav * tarMass / tarRadi;
				ith.id   = id++;
				ith.mat = EoS::Granite;
				ith.tag = 1;
				if(ith.id / NptclIn1Node == PS::Comm::getRank()) imp.push_back(ith);
			}
		}
	}
	for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
		for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
			for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
				const PS::F64 r = impCoreShrinkFactor * sqrt(x*x + y*y + z*z);
				if(r + dx >= impFracRadi * coreFracRadi) continue;
				RealPtcl ith;
				const double impCoreRadi = impCoreShrinkFactor * UnitRadi;
				ith.pos.x = impCoreShrinkFactor * UnitRadi * x + offset;
				ith.pos.y = impCoreShrinkFactor * UnitRadi * y;
				ith.pos.z = impCoreShrinkFactor * UnitRadi * z;
				ith.dens = impCoreMass / (4.0 / 3.0 * math::pi * impCoreRadi * impCoreRadi * impCoreRadi);
				ith.mass = tarMass + impMass;
				ith.eng  = 0.1 * sysinfo.Grav * tarMass / tarRadi;
				ith.id   = id++;
				ith.mat = EoS::Iron;
				ith.tag = 1;
				if(ith.id / NptclIn1Node == PS::Comm::getRank()){
					imp.push_back(ith);
				}
			}
		}
	}
	for(PS::U32 i = 0 ; i < tar.size() ; ++ i){
		tar[i].mass /= (PS::F64)(Nptcl);
	}
	for(PS::U32 i = 0 ; i < imp.size() ; ++ i){
		imp[i].mass /= (PS::F64)(Nptcl);
	}
	
	for(PS::U32 i = 0 ; i < tar.size() ; ++ i){
		ptcl.push_back(tar[i]);
	}
	for(PS::U32 i = 0 ; i < imp.size() ; ++ i){
		//ptcl.push_back(imp[i]);
	}

	const PS::S32 numPtclLocal = ptcl.size();
	sph_system.setNumberOfParticleLocal(numPtclLocal);
	for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
		sph_system[i] = ptcl[i];
	}
	//Fin.
	std::cout << "# of ptcls = " << ptcl.size() << std::endl;
	std::cout << "setup..." << std::endl;
}
#else
void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, system_t& sysinfo){
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
	const PS::F64 impMass = 0.125 * EarthMass;
	const PS::F64 impRadi = 0.5 * EarthRadi;
	sysinfo.Grav = 6.67e-11;
	sysinfo.end_time = 1.0e+4;

	const PS::F64 dx = 1.0 / 16.0;
	PS::S32 id = 0;
	//Dummy put to determine # of ptcls
	int Nptcl = 0;
	for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
		for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
			for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
				const PS::F64 r = sqrt(x*x + y*y + z*z);
				if(r >= 1.0) continue;
				++ Nptcl;
			}
		}
	}
	for(PS::F64 x = -1.0 ; x <= 1.0 ; x += 2.0 * dx){
		for(PS::F64 y = -1.0 ; y <= 1.0 ; y += 2.0 * dx){
			for(PS::F64 z = -1.0 ; z <= 1.0 ; z += 2.0 * dx){
				const PS::F64 r = sqrt(x*x + y*y + z*z);
				if(r >= 1.0) continue;
				++ Nptcl;
			}
		}
	}
	const int NptclIn1Node = Nptcl / PS::Comm::getNumberOfProc();
	//Dummy end

	//Put Tar.
	for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
		for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
			for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
				const PS::F64 r = sqrt(x*x + y*y + z*z);
				if(r >= 1.0) continue;
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
				if(ith.id / NptclIn1Node == PS::Comm::getRank()) tar.push_back(ith);
			}
		}
	}
	//Put Imp.
	for(PS::F64 x = -1.0 ; x <= 1.0 ; x += 2.0 * dx){
		for(PS::F64 y = -1.0 ; y <= 1.0 ; y += 2.0 * dx){
			for(PS::F64 z = -1.0 ; z <= 1.0 ; z += 2.0 * dx){
				const PS::F64 r = sqrt(x*x + y*y + z*z);
				if(r >= 1.0) continue;
				RealPtcl ith;
				ith.pos.x = impRadi * x + 5.0 * tarRadi;
				ith.pos.y = impRadi * y;
				ith.pos.z = impRadi * z;
				ith.dens = impMass / (4.0 / 3.0 * math::pi * impRadi * impRadi * impRadi);
				ith.mass = tarMass + impMass;
				//ith.eng  = 0.1 * sysinfo.Grav * impMass / impRadi;
				ith.eng  = 0.1 * sysinfo.Grav * tarMass / tarRadi;
				ith.id   = id++;
				ith.mat = EoS::Granite;
				ith.tag = 1;
				if(ith.id / NptclIn1Node == PS::Comm::getRank()) imp.push_back(ith);
			}
		}
	}
	for(PS::U32 i = 0 ; i < tar.size() ; ++ i){
		tar[i].mass /= (PS::F64)(Nptcl);
	}
	for(PS::U32 i = 0 ; i < imp.size() ; ++ i){
		imp[i].mass /= (PS::F64)(Nptcl);
	}

	for(PS::U32 i = 0 ; i < tar.size() ; ++ i){
		ptcl.push_back(tar[i]);
	}
	for(PS::U32 i = 0 ; i < imp.size() ; ++ i){
		ptcl.push_back(imp[i]);
	}
	const PS::S32 numPtclLocal = ptcl.size();
	sph_system.setNumberOfParticleLocal(numPtclLocal);
	for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
		sph_system[i] = ptcl[i];
	}
	//Fin.
	std::cout << "# of ptcls = " << ptcl.size() << std::endl;
	std::cout << "setup..." << std::endl;
}
#endif
#else

void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, system_t& sysinfo){
	FileHeader header;
	//sph_system.readParticleAscii("init.dat", header);
	char filename[256];
	sprintf(filename, "init");
	sph_system.readParticleAscii(filename, "%s_%05d_%05d.dat", header);
	std::cout << "TEST" << std::endl;
	std::cout << "TEST2" << std::endl;

	//Fin.
	sysinfo.Grav = 6.67e-11;
	sysinfo.end_time = 1.0e+5;

	//Shift positions
	std::vector<RealPtcl> tar, imp;
	PS::F64vec pos_tar = 0;
	PS::F64vec pos_imp = 0;
	PS::F64vec vel_imp = 0;
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
			vel_imp  += sph_system[i].mass * sph_system[i].vel;
			mass_imp += sph_system[i].mass;
		}
	}
	//accumurate
	pos_tar  = PS::Comm::getSum( pos_tar);
	pos_imp  = PS::Comm::getSum( pos_imp);
	vel_imp  = PS::Comm::getSum( vel_imp);
	mass_tar = PS::Comm::getSum(mass_tar);
	mass_imp = PS::Comm::getSum(mass_imp);
	pos_tar /= mass_tar;
	pos_imp /= mass_imp;
	vel_imp /= mass_imp;
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
	const double L_init = L_EM * 1.21;
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
			sph_system[i].vel -= vel_imp;
		}
	}
	std::cout << "setup..." << std::endl;
}
#endif

