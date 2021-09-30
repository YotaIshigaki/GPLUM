#pragma once


#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
template <class Ptcl> class Imp: public Problem<Ptcl>{
	static double pres(const double z){
		return 1000.0 * 9.8 * (550000.0 * 0.05 - z);
	}
	public:
	static const double END_TIME = 2.0e-4;
	static void setupIC(PS::ParticleSystem<Ptcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
		/////////
		//place ptcls
		/////////
		std::vector<Ptcl> ptcl;
		const PS::F64 box_x = 0.05;
		const PS::F64 dx = box_x / 1024;
		const double dens_tar = 1000.0;
		const double radi = 1.0e-3;
		const double mass = box_x * box_x * dens_tar + math::pi * radi * radi * 2.0 * dens_tar;

		PS::S32 i = 0;
		for(PS::F64 x = 0 ; x < box_x - 0.5 * dx ; x += dx){
			for(PS::F64 y = dx ; y <= box_x ; y += dx){
				Ptcl ith;
				ith.type = HYDRO;
				ith.pos.x = x;
				ith.pos.y = y;
				ith.dens = dens_tar;
				ith.mass = mass;
				ith.setPressure(&Water);
				ith.tag = 0;
				double eng;
				const double pressure = pres(y);
				{
					double u_prev, u_curr, u_next;
					u_curr = 1.0e+6;
					u_prev = 1.0e+5;
					for(int i = 0 ; i < 100 ; ++ i){
						u_next = u_curr - (ith.EoS->Pressure(ith.dens, u_curr) - pressure) * (u_curr - u_prev) / (ith.EoS->Pressure(ith.dens, u_curr) - ith.EoS->Pressure(ith.dens, u_prev));
						if(std::abs(u_next / u_curr - 1.0) < 1.0e-5) break;
						u_prev = u_curr;
						u_curr = u_next;
					}
					eng = u_next;
				}

				ith.eng  = eng;
				ith.id   = i++;
				if(y <= 0.1 * box_x){
					ith.type = FREEZE;
				}
				ptcl.push_back(ith);
			}
		}

		for(PS::F64 x = 0 ; x < 2.0 * box_x ; x += dx){
			for(PS::F64 y = box_x + dx ; y < 2.0 * box_x ; y += dx){
				Ptcl ith;
				ith.type = HYDRO;
				PS::F64vec pos = PS::F64vec(x, y);
				pos -= PS::F64vec(0.5 * box_x, 1.5 * box_x - (0.5 - 4.0 * radi / box_x) * box_x);
				if(sqrt(pos * pos) > radi - 0.9 * dx){
					continue;
				}
				ith.pos.x = x;
				ith.pos.y = y + 0.0025;
				ith.dens = 2.0 * dens_tar;
				ith.mass = mass;
				ith.setPressure(&WetTuff);
				ith.tag = 1;
				ith.vel.y = -4.64e+3;
				double eng;
				const double pressure = pres(box_x) * 4.8;
				{
					double u_prev, u_curr, u_next;
					u_curr = 1.0e+6;
					u_prev = 1.0e+5;
					for(int i = 0 ; i < 100 ; ++ i){
						u_next = u_curr - (ith.EoS->Pressure(ith.dens, u_curr) - pressure) * (u_curr - u_prev) / (ith.EoS->Pressure(ith.dens, u_curr) - ith.EoS->Pressure(ith.dens, u_prev));
						if(std::abs(u_next / u_curr - 1.0) < 1.0e-5) break;
						u_prev = u_curr;
						u_curr = u_next;
					}
					eng = u_next;
				}

				ith.eng  = eng;
				ith.id   = i++;
				if(y <= 0.1 * box_x){
					ith.type = FREEZE;
				}
				ptcl.push_back(ith);

			}
		}

		for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
			ptcl[i].mass = mass / (PS::F64)(ptcl.size());
		}
		std::cout << "# of ptcls is... " << ptcl.size() << std::endl;
		//
		dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_X);
		dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0), PS::F64vec(box_x, box_x));
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
	static void setEoS(PS::ParticleSystem<Ptcl>& sph_system){
		for(PS::U64 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			if(sph_system[i].tag == 0){
				sph_system[i].setPressure(&Water);
			}else{
				sph_system[i].setPressure(&WetTuff);
			}
		}
	}
	static void addExternalForce(PS::ParticleSystem<Ptcl>& sph_system, system_t& system){
		for(PS::U64 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].acc.y -= 9.8;
		}
	}
};
#else
#error NO
#endif

