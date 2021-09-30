#define SELF_GRAVITY
#define FLAG_GI
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
#error
#endif
template <class Ptcl> class GI : public Problem<Ptcl>{
	public:
	static const double END_TIME = 1.0e+4;
	static void setupIC(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
		const double Corr = .98;//Correction Term
		/////////
		//place ptcls
		/////////
		std::vector<Ptcl> ptcl;
		std::vector<Ptcl> tar;//Target
		std::vector<Ptcl> imp;//Impactor
		/////////
		const PS::F64 UnitMass = 6.0e+24; // Earth mass
		const PS::F64 UnitRadi = 6400e+3; // Earth radii
		//total-to-core frac.
		const PS::F64 coreFracRadi = 3500.0e+3 / 6400.0e+3;//Earth
		const PS::F64 coreFracMass = 0.3;//Earth
		/////////
		const PS::F64 Expand = 1.1;
		const PS::F64 tarMass = UnitMass;
		const PS::F64 tarRadi = UnitRadi;
		const PS::F64 tarCoreMass = tarMass * coreFracMass;
		const PS::F64 tarCoreRadi = tarRadi * coreFracRadi;
		const PS::F64 impMass = 0.1 * tarMass;
		const PS::F64 impRadi = Expand * cbrt(impMass / tarMass) * UnitRadi;
		const PS::F64 impCoreMass = impMass * coreFracMass;
		const PS::F64 impCoreRadi = impRadi * coreFracRadi;

		const double offset = 5.0 * UnitRadi;
		const PS::F64 dx = 1.0 / 39;
		const PS::F64 Grav = 6.67e-11;
		std::cout << impRadi / tarRadi << std::endl;
		std::cout << impCoreRadi / impRadi << std::endl;
		///////////////////
		//Dummy put to determine # of ptcls
		///////////////////
		//target
		int tarNmntl = 0;
		for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
			for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
				for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
					const PS::F64 r = sqrt(x*x + y*y + z*z) * UnitRadi;
					if(r >= tarRadi || r <= tarCoreRadi) continue;
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
						const PS::F64 r = tarCoreShrinkFactor * sqrt(x*x + y*y + z*z) * UnitRadi;
						if(r >= Corr * tarCoreRadi) continue;
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
					const PS::F64 r = Expand * sqrt(x*x + y*y + z*z) * UnitRadi;
					if(r >= impRadi || r <= impCoreRadi) continue;
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
						const PS::F64 r = Expand * impCoreShrinkFactor * sqrt(x*x + y*y + z*z) * UnitRadi;
						if(r >= Corr * impCoreRadi) continue;
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
		std::cout << "    radius           : " << tarRadi << std::endl;
		std::cout << "    total-to-core    : " << (double)(tarNcore) / (double)(tarNptcl) << std::endl;
		std::cout << "    # of core ptcls  : " << tarNcore << std::endl;
		std::cout << "    # of mantle ptcls: " << tarNmntl << std::endl;
		std::cout << "    core density     : " << tarCoreMass / (4.0 * math::pi / 3.0 * tarCoreRadi * tarCoreRadi * tarCoreRadi * Corr * Corr * Corr) << std::endl;
		std::cout << "    mantle density   : " << (tarMass - tarCoreMass) / (4.0 * math::pi / 3.0 * (tarRadi * tarRadi * tarRadi - tarCoreRadi * tarCoreRadi * tarCoreRadi)) << std::endl;
		std::cout << "    mean density     : " << tarMass / (4.0 * math::pi / 3.0 * tarRadi * tarRadi * tarRadi) << std::endl;
		std::cout << "Impactor:" << impNptcl << std::endl;
		std::cout << "    radius           : " << impRadi << std::endl;
		std::cout << "    total-to-core    : " << (double)(impNcore) / (double)(impNptcl) << std::endl;
		std::cout << "    # of core ptcls  : " << impNcore << std::endl;
		std::cout << "    # of mantle ptcls: " << impNmntl << std::endl;
		std::cout << "    core density     : " << impCoreMass / (4.0 * math::pi / 3.0 * impCoreRadi * impCoreRadi * impCoreRadi * Corr * Corr * Corr) << std::endl;
		std::cout << "    mantle density   : " << (impMass - impCoreMass) / (4.0 * math::pi / 3.0 * (impRadi * impRadi * impRadi - impCoreRadi * impCoreRadi * impCoreRadi)) << std::endl;
		std::cout << "    mean density     : " << impMass / (4.0 * math::pi / 3.0 * impRadi * impRadi * impRadi) << std::endl;
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
					if(r >= tarRadi || r <= tarCoreRadi) continue;
					Ptcl ith;
					ith.pos.x = UnitRadi * x;
					ith.pos.y = UnitRadi * y;
					ith.pos.z = UnitRadi * z;
					ith.dens = (tarMass - tarCoreMass) / (4.0 / 3.0 * math::pi * (tarRadi * tarRadi * tarRadi - tarCoreRadi * tarCoreRadi * tarCoreRadi));
					ith.mass = tarMass + impMass;
					ith.eng  = 0.1 * Grav * tarMass / tarRadi;
					ith.id   = id++;
					ith.setPressure(&Granite);
					ith.tag = 0;
					if(ith.id / NptclIn1Node == PS::Comm::getRank()) tar.push_back(ith);
				}
			}
		}
		for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
			for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
				for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
					const PS::F64 r = tarCoreShrinkFactor * sqrt(x*x + y*y + z*z) * UnitRadi;
					if(r >= Corr * tarCoreRadi) continue;
					Ptcl ith;
					ith.pos.x = tarCoreShrinkFactor * UnitRadi * x;
					ith.pos.y = tarCoreShrinkFactor * UnitRadi * y;
					ith.pos.z = tarCoreShrinkFactor * UnitRadi * z;
					ith.dens = tarCoreMass / (4.0 / 3.0 * math::pi * tarCoreRadi * tarCoreRadi * tarCoreRadi * Corr * Corr * Corr);
					ith.mass = tarMass + impMass;
					ith.eng  = 0.1 * Grav * tarMass / tarRadi;
					ith.id   = id++;
					ith.setPressure(&Iron);
					ith.tag = 1;
					if(ith.id / NptclIn1Node == PS::Comm::getRank()) tar.push_back(ith);
				}
			}
		}
		//imp
		for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
			for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
				for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
					const PS::F64 r = Expand * sqrt(x*x + y*y + z*z) * UnitRadi;
					if(r >= impRadi || r <= impCoreRadi) continue;
					Ptcl ith;
					ith.pos.x = Expand * UnitRadi * x + offset;
					ith.pos.y = Expand * UnitRadi * y;
					ith.pos.z = Expand * UnitRadi * z;
					ith.dens = (impMass - impCoreMass) / (4.0 / 3.0 * math::pi * (impRadi * impRadi * impRadi - impCoreRadi * impCoreRadi * impCoreRadi));
					ith.mass = tarMass + impMass;
					ith.eng  = 0.1 * Grav * tarMass / tarRadi;
					ith.id   = id++;
					ith.setPressure(&Granite);
					ith.tag = 2;
					if(ith.id / NptclIn1Node == PS::Comm::getRank()) imp.push_back(ith);
				}
			}
		}
		for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
			for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
				for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
					const PS::F64 r = Expand * impCoreShrinkFactor * sqrt(x*x + y*y + z*z) * UnitRadi;
					if(r >= impCoreRadi) continue;
					Ptcl ith;
					ith.pos.x = Expand * impCoreShrinkFactor * UnitRadi * x + offset;
					ith.pos.y = Expand * impCoreShrinkFactor * UnitRadi * y;
					ith.pos.z = Expand * impCoreShrinkFactor * UnitRadi * z;
					ith.dens = impCoreMass / (4.0 / 3.0 * math::pi * impCoreRadi * impCoreRadi * impCoreRadi * Corr * Corr * Corr);
					ith.mass = tarMass + impMass;
					ith.eng  = 0.1 * Grav * tarMass / tarRadi;
					ith.id   = id++;
					ith.setPressure(&Iron);
					ith.tag = 3;
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

	static void setEoS(PS::ParticleSystem<Ptcl>& sph_system){
		for(PS::U64 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			if(sph_system[i].tag % 2 == 0){
				sph_system[i].setPressure(&Granite);
			}else{
				sph_system[i].setPressure(&Iron);
			}
		}
	}

	static void addExternalForce(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo){
		if(sysinfo.time >= 5000) return;
		std::cout << "Add Ext. Force!!!" << std::endl;
		#pragma omp parallel for
		for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].acc += - sph_system[i].vel * 0.05 / sph_system[i].dt;
		}
	}
};

