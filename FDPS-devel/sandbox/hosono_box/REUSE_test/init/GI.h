#define SELF_GRAVITY
#define FLAG_GI
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
#error
#endif
#define GRAVITY
template <class Ptcl> class GI_init : public Problem<Ptcl>{
	public:
	static const PS::F64 UnitMass = 6.0e+24; // Earth mass
	static const PS::F64 UnitRadi = 6400e+3; // Earth radii
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

	static void setupIC(PS::ParticleSystem<Ptcl>& sph_system, PS::DomainInfo& dinfo){
		int putImpactor = 0;
		const PS::F64 dx = 1.0 / 100.0;
		//const PS::F64 dx = 1.0 / 39.0;
		//std::cin >> putImpactor;

		/////////
		//place ptcls
		/////////
		std::vector<Ptcl> ptcl;
		std::vector<Ptcl> tar;//Target
		std::vector<Ptcl> imp;//Impactor
		/////////
		/////////
		const PS::F64 Expand = 1.1;
		const PS::F64 tarMass = UnitMass;
		const PS::F64 tarRadi = UnitRadi;
		const PS::F64 tarCoreMass = 0.3 * tarMass;
		const PS::F64 tarCoreRadi = cbrt(tarCoreMass / tarMass) * tarRadi;

		const PS::F64 impMass = 0.1 * tarMass;
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
					ith.mat = (sqrt(ith.pos * ith.pos) < tarCoreRadi) ? 1 : 0;
					if(ith.mat == 1){
						ith.EoS  = &Iron;
					}else{
						ith.EoS  = &Granite;
					}
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
					ith.mat = (sqrt(ith.pos * ith.pos) < impCoreRadi) ? 3 : 2;
					if(ith.mat == 3){
						ith.EoS  = &Iron;
					}else{
						ith.EoS  = &Granite;
					}
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
	}
};


