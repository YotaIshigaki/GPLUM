#define SELF_GRAVITY
#define FLAG_GI
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
#error
#endif

template <class Ptcl> class GI_init : public Problem<Ptcl>{
	public:
	static const PS::F64 UnitMass = 6.0e+24 * 0.55; // Earth mass
	static const PS::F64 UnitRadi = 6400e+3 * 1.; // Earth radii
	static const PS::F64 Grav = 6.67e-11;
	
	static PS::F64 UnitTime(){
		return sqrt(UnitRadi * UnitRadi * UnitRadi / (Grav * UnitMass));
	}
	static PS::F64 DynamicalTime(){
		return math::pi * UnitTime();
	}
	static PS::F64 EscapeVelocity(){
		return sqrt(2.0 * Grav * UnitMass / UnitRadi);
	}
	static const double END_TIME = 1.0e+4;

	static void setupIC(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
		int putImpactor = 1;
		const PS::F64 dx = 1.0 / 50.0;
		//std::cin >> putImpactor;

		std::cout << "putImpactor = " << putImpactor << std::endl;
		/////////
		//place ptcls
		/////////
		std::vector<Ptcl> ptcl;
		std::vector<Ptcl> tar;//Target
		std::vector<Ptcl> imp;//Impactor
		/////////
		/////////
		const PS::F64 Expand = 1.0;//1.1 for EM
		const PS::F64 tarMass = UnitMass;
		const PS::F64 tarRadi = UnitRadi;
		const PS::F64 tarCoreMass = 0.3 * tarMass;
		const PS::F64 tarCoreRadi = cbrt(tarCoreMass / tarMass) * tarRadi;

		const PS::F64 impMass = tarMass;
		const PS::F64 impRadi = Expand * cbrt(impMass / tarMass) * UnitRadi;
		const PS::F64 impCoreMass = 0.3 * impMass;
		const PS::F64 impCoreRadi = cbrt(impCoreMass / impMass) * impRadi;

		const double offset = 0.0 * UnitRadi;
		///////////////////
		//Dummy put to determine # of ptcls
		///////////////////
		//target
		int tarNmntl = 0;
		for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
			for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
				for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
					const PS::F64 r = sqrt(x*x + y*y + z*z) * UnitRadi;
					if(r >= tarRadi) continue;
					++ tarNmntl;
				}
			}
		}
		//imp
		int impNmntl = 0;
		for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
			for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
				for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
					const PS::F64 r = Expand * sqrt(x*x + y*y + z*z) * UnitRadi;
					if(r >= impRadi) continue;
					++ impNmntl;
				}
			}
		}
		///////////////////
		//Dummy end
		///////////////////
		const int tarNptcl = tarNmntl;
		const int impNptcl = impNmntl;
		const int Nptcl    = tarNptcl + impNptcl;
		std::cout << "Target  :" << tarNptcl << std::endl;
		std::cout << "    radius           : " << tarRadi << std::endl;
		std::cout << "    # of mantle ptcls: " << tarNmntl << std::endl;
		std::cout << "    mantle density   : " << (tarMass) / (4.0 * math::pi / 3.0 * (tarRadi * tarRadi * tarRadi)) << std::endl;
		std::cout << "    core radius      : " << tarCoreRadi << std::endl;
		std::cout << "Impactor:" << impNptcl << std::endl;
		std::cout << "    radius           : " << impRadi << std::endl;
		std::cout << "    # of mantle ptcls: " << impNmntl << std::endl;
		std::cout << "    mantle density   : " << (impMass) / (4.0 * math::pi / 3.0 * (impRadi * impRadi * impRadi)) << std::endl;
		std::cout << "Total:" << Nptcl << std::endl;
		std::cout << "Tar-to-Imp mass ratio: " << (double)(impNmntl) / (double)(tarNmntl) << std::endl;
		const int NptclIn1Node = Nptcl / PS::Comm::getNumberOfProc();
		///////////////////
		//Real put
		///////////////////
		PS::S32 id = 0;
		//Put Tar.
		for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
			for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
				for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
					const PS::F64 r = sqrt(x*x + y*y + z*z) * UnitRadi;
					if(r >= tarRadi) continue;
					Ptcl ith;
					ith.pos.x = UnitRadi * x;
					ith.pos.y = UnitRadi * y;
					ith.pos.z = UnitRadi * z;
					/*
					ith.vel.x = EscapeVelocity() * 0.1 * (double)rand() / (double)RAND_MAX;
					ith.vel.y = EscapeVelocity() * 0.1 * (double)rand() / (double)RAND_MAX;
					ith.vel.z = EscapeVelocity() * 0.1 * (double)rand() / (double)RAND_MAX;
					*/
					ith.dens = (tarMass) / (4.0 / 3.0 * math::pi * (tarRadi * tarRadi * tarRadi));
					ith.mass = tarMass + impMass;
					ith.eng  = 0.1 * Grav * tarMass / tarRadi;
					ith.id   = id++;
					ith.setPressure(&Iron);
					ith.tag = (sqrt(ith.pos * ith.pos) < tarCoreRadi) ? 1 : 0;
					if(ith.id / NptclIn1Node == PS::Comm::getRank()) tar.push_back(ith);
				}
			}
		}
		//imp
		for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
			for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
				for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
					const PS::F64 r = Expand * sqrt(x*x + y*y + z*z) * UnitRadi;
					if(r >= impRadi) continue;
					Ptcl ith;
					ith.pos.x = Expand * UnitRadi * x;
					ith.pos.y = Expand * UnitRadi * y;
					ith.pos.z = Expand * UnitRadi * z;
					ith.dens = (impMass) / (4.0 / 3.0 * math::pi * (impRadi * impRadi * impRadi));
					ith.mass = tarMass + impMass;
					ith.eng  = 0.1 * Grav * impMass / impRadi;
					ith.id   = id++;
					ith.setPressure(&Granite);
					ith.tag = (sqrt(ith.pos * ith.pos) < impCoreRadi) ? 3 : 2;
					ith.pos.x += offset;
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

		if(putImpactor == 0){
			for(PS::U32 i = 0 ; i < tar.size() ; ++ i){
				ptcl.push_back(tar[i]);
			}
		}else{
			for(PS::U32 i = 0 ; i < imp.size() ; ++ i){
				ptcl.push_back(imp[i]);
			}
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
	static void setEoS(PS::ParticleSystem<Ptcl>& sph_system){
		#pragma omp parallel for
		for(PS::U64 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			if(sph_system[i].tag >= 2){
				sph_system[i].setPressure(&Granite);
			}else{
				sph_system[i].setPressure(&Iron);
			}
		}
	}
	static void addExternalForce(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo){
		if(sysinfo.time < 5000.0){
			std::cout << "Add Frictional Force." << std::endl;
			#pragma omp parallel for
			for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
				sph_system[i].acc += - sph_system[i].vel * 0.05 / sph_system[i].dt;
			}
		}
	}
};

template <class Ptcl> class GI_init2 : public Problem<Ptcl>{
	public:
	static const double END_TIME = 1.0e+4;
	static void setupIC(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
		FileHeader header;
		char filename[256];
		sprintf(filename, "init");
		//sph_system.readParticleAscii(filename, "init/%s_%05d_%05d.dat", header);
		sph_system.readParticleAscii(filename, "%s_%05d_%05d.dat", header);
	}
	static void setEoS(PS::ParticleSystem<Ptcl>& sph_system){
		#pragma omp parallel for
		for(PS::U64 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			if(sph_system[i].tag % 2 == 1){
				sph_system[i].setPressure(&Iron);
			}else{
				sph_system[i].setPressure(&Granite);
			}
		}
	}
	static void addExternalForce(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo){
		if(sysinfo.time < 5000.0){
			std::cout << "Add Frictional Force." << std::endl;
			#pragma omp parallel for
			for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
				sph_system[i].acc += - sph_system[i].vel * 0.05 / sph_system[i].dt;
			}
		}
	}
};

template <class Ptcl> class GI_init3 : public Problem<Ptcl>{
	public:
	static const double END_TIME = 4.0e+4;
	static const double RotationPeriod = 9.0;//hours
	static void setupIC(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
		FileHeader header;
		char filename[256];
		sprintf(filename, "init");
		sph_system.readParticleAscii(filename, "%s_%05d_%05d.dat", header);
		const PS::F64 ang_vel = 2.0 * math::pi / (RotationPeriod * 3600.0);
		PS::F64vec AM_loc = 0.0;
		for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			const PS::F64 theta = atan2(sph_system[i].pos.y, sph_system[i].pos.x);
			const PS::F64vec n = PS::F64vec(-sin(theta), cos(theta), 0.0);
			sph_system[i].vel += sqrt(sph_system[i].pos * sph_system[i].pos) * ang_vel * n;
			AM_loc += sph_system[i].mass * sph_system[i].pos ^ sph_system[i].vel;
		}
		PS::F64vec AM = PS::Comm::getSum(AM_loc);
		std::cout << "Spin AM: "<< AM.z / (3.5e+34) << std::endl;
	}
	static void setEoS(PS::ParticleSystem<Ptcl>& sph_system){
		#pragma omp parallel for
		for(PS::U64 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			if(sph_system[i].tag % 2 == 1){
				sph_system[i].setPressure(&Iron);
			}else{
				sph_system[i].setPressure(&Granite);
			}
		}
	}
	static void addExternalForce(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo){
		if(sysinfo.time < 0.5 * END_TIME){
			std::cout << "Add Frictional Force. to rad." << std::endl;
			#pragma omp parallel for
			for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
				sph_system[i].acc += - 0.1 * (sph_system[i].pos * sph_system[i].vel) / (sph_system[i].pos * sph_system[i].pos) / sysinfo.dt * sph_system[i].pos;
			}
		}
	}
};


