#define SELF_GRAVITY
#define FLAG_GI
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
#error
#endif
template <class Ptcl> class GI : public Problem<Ptcl>{
	public:
	static const double END_TIME = 1.0e+5;
	static const PS::F64 R = 6400.0e+3;
	static const double Grav = 6.67e-11;
	static const double L_EM = 3.5e+34;
	//static const PS::F64 AngMomEarthSpin = 7.1e+33;???
	static const PS::F64 InertiaEarth    = 8.1e+37;

	static void setupIC(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
		FileHeader header;
		char filename[256];
		sprintf(filename, "init");
		//std::cout << "TEST" << std::endl;
		sph_system.readParticleAscii(filename, "%s_%05d_%05d.dat", header);
		//Shift positions
		std::vector<Ptcl> tar, imp;
		PS::F64vec pos_tar = 0;
		PS::F64vec pos_imp = 0;
		PS::F64vec vel_imp = 0;
		PS::F64 mass_tar = 0;
		PS::F64 mass_imp = 0;
		PS::F64 radi_tar = 0;
		PS::F64 radi_imp = 0;
		for(PS::U32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			if(sph_system[i].tag <= 1){
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
			if(sph_system[i].tag <= 1){
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

		const double v_esc = sqrt(2.0 * Grav * (mass_tar + mass_imp) / (radi_tar + radi_imp));
		const double x_init = 3.0 * radi_tar;
		double input = 0;
		std::cout << "Input L_init" << std::endl;
		std::cin >> input;
		const double L_init = L_EM * input;
		std::cout << "Input v_imp" << std::endl;
		std::cin >> input;
		const double v_imp = v_esc * input;
		std::cout << "Input rotation angle on yz-plane (1 = prograde; -1 = retrograde; 0 = ...)" << std::endl;
		std::cin >> input;
		const double angle = acos(input);

		const double v_inf = sqrt(std::max(v_imp * v_imp - v_esc * v_esc, 0.0));
		double y_init = radi_tar;//Initial guess.
		double v_init;
		std::cout << "v_esc = " << v_esc << std::endl;
		for(int it = 0 ; it < 10 ; ++ it){
			v_init = sqrt(v_inf * v_inf + 2.0 * Grav * mass_tar / sqrt(x_init * x_init + y_init * y_init));
			y_init = L_init / (mass_imp * v_init);
		}

		std::cout << "v_init = "  << v_init <<  std::endl;
		std::cout << "y_init / Rtar = "  << y_init / radi_tar <<  std::endl;
		std::cout << "v_imp  = "  << v_imp <<  std::endl;
		std::cout << "m_imp  = " << mass_imp / 6.0e+24 << std::endl;
		//shift'em
		for(PS::U32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			if(sph_system[i].tag % 2 == 0){
				sph_system[i].setPressure(&Granite);
			}else{
				sph_system[i].setPressure(&Iron);
			}
			if(sph_system[i].tag <= 1){
			}else{
				sph_system[i].pos = sph_system[i].pos - pos_imp;
				sph_system[i].pos.x += x_init;
				sph_system[i].pos.y += y_init;
				sph_system[i].vel.x -= v_init;
				sph_system[i].vel -= vel_imp;

				const double y =   cos(angle) * sph_system[i].pos.y + sin(angle) * sph_system[i].pos.z;
				const double z = - sin(angle) * sph_system[i].pos.y + cos(angle) * sph_system[i].pos.z;
				sph_system[i].pos.y = y;
				sph_system[i].pos.z = z;
			}
		}
		std::cout << "setup..." << std::endl;
	}
	static void postTimestepProcess(PS::ParticleSystem<Ptcl>& sph_system, system_t& sys){
		//SANITY Check
		for(PS::U64 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].eng  = std::max(sph_system[i].eng , 0.0);
			sph_system[i].dens = std::max(sph_system[i].dens, 5.0);
		}
		//Shift Origin
		PS::F64vec com_loc = 0;//center of mass of target core
		PS::F64vec mom_loc = 0;//moment of target core
		PS::F64 mass_loc = 0;//mass of target core
		for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			if(sph_system[i].tag != 1) continue;
			com_loc += sph_system[i].pos * sph_system[i].mass;
			mom_loc += sph_system[i].vel * sph_system[i].mass;
			mass_loc += sph_system[i].mass;
		}
		PS::F64vec com = PS::Comm::getSum(com_loc);
		PS::F64vec mom = PS::Comm::getSum(mom_loc);
		PS::F64 mass = PS::Comm::getSum(mass_loc);
		com /= mass;
		mom /= mass;
		for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].pos -= com;
			sph_system[i].vel -= mom;
		}
		#if 0
		std::size_t Nptcl = sph_system.getNumberOfParticleLocal();
		for(PS::S32 i = 0 ; i < Nptcl ; ++ i){
			if(sqrt(sph_system[i].pos * sph_system[i].pos) / R > 40.0){
				//bounded particles should not be killed.
				if(0.5 * sph_system[i].vel * sph_system[i].vel + sph_system[i].pot < 0) continue;
				std::cout << "KILL" << std::endl;
				sph_system[i] = sph_system[-- Nptcl];
			}
		}
		sph_system.setNumberOfParticleLocal(Nptcl);
		#endif
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
};

