#pragma once

template <typename real> real pres(const real y){
	return 2700.0 * 9.8 * (65000.0 - y);
}

#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
template <class RealPtcl> void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
	const PS::F64 imp_radi = 10.0;
	const PS::F64 ref_dens = 2700.0;
	const PS::F64 dx = imp_radi / 64;
	const PS::F64 box_x = 10.0 * imp_radi;
	const PS::F64 box_y = box_x;
	/////////
	//place ptcls
	/////////
	int Nptcl = 0;
	for(PS::F64 x = 0 ; x < box_x ; x += dx){
		for(PS::F64 y = 0 ; y < box_y ; y += dx){
			++ Nptcl;
		}
	}

	for(PS::F64 x = - imp_radi ; x <= imp_radi ; x += dx){
		for(PS::F64 y = - imp_radi ; y <= imp_radi ; y += dx){
			++ Nptcl;
		}
	}
	const int NptclIn1Node = Nptcl / PS::Comm::getNumberOfProc();

	//Real put
	int id = 0;
	std::vector<RealPtcl> ptcl;
	for(PS::F64 x = 0 ; x < box_x ; x += dx){
		for(PS::F64 y = 0 ; y < box_y ; y += dx){
			RealPtcl ith;
			ith.type = HYDRO;
			ith.pos.x = x;
			ith.pos.y = y;
			ith.dens = ref_dens;
			ith.setPressure(&Granite);
			if(y < 0.1 * box_y){
				ith.type = FREEZE;
			}
			ith.tag = 0;
			ith.id = id++;
			if(ith.id / NptclIn1Node == PS::Comm::getRank()) ptcl.push_back(ith);
		}
	}

	for(PS::F64 x = - imp_radi ; x <= imp_radi ; x += dx){
		for(PS::F64 y = - imp_radi ; y <= imp_radi ; y += dx){
			RealPtcl ith;
			ith.type = HYDRO;
			const PS::F64 r = sqrt(x*x + y*y);
			if(r >= imp_radi) continue;
			ith.pos.x = x + 0.5 * box_x;
			ith.pos.y = y + box_y + imp_radi + 5.0 * dx;
			ith.vel.y = -10.0e+3;
			ith.dens = 2700.0;
			ith.setPressure(&Granite);
			ith.tag = 1;
			ith.id = id++;
			if(ith.id / NptclIn1Node == PS::Comm::getRank()) ptcl.push_back(ith);
		}
	}

	for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
		ptcl[i].mass = (ref_dens * box_x * box_y + 4.0 / 3.0 * math::pi * imp_radi * imp_radi * ref_dens) / (PS::F64)(Nptcl);
		ptcl[i].grav.y = -9.8;

		double eng;
		const double pressure = pres(ptcl[i].pos.y);
		const EoS::Tillotson<PS::F64> Granite(2680.0, 16.0e+6, 3.5e+6, 18.00e+6, 18.0e+9, 18.0e+9, 0.5, 1.3, 5.0, 5.0);
		{
			double u_prev, u_curr, u_next;
			u_curr = 1.0e+6;
			u_prev = 1.0e+5;
			for(int i = 0 ; i < 100 ; ++ i){
				u_next = u_curr - (Granite.Pressure(ptcl[i].dens, u_curr) - pressure) * (u_curr - u_prev) / (Granite.Pressure(ptcl[i].dens, u_curr) - Granite.Pressure(ptcl[i].dens, u_prev));
				if(std::abs(u_next / u_curr - 1.0) < 1.0e-5) break;
				u_prev = u_curr;
				u_curr = u_next;
			}
			eng = u_next;
		}
		ptcl[i].eng = eng;
	}
	std::cout << "# of ptcls is... " << ptcl.size() << std::endl;
	//
	dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_X);
	dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0), PS::F64vec(box_x, box_y));
	const PS::S32 numPtclLocal = ptcl.size();
	sph_system.setNumberOfParticleLocal(numPtclLocal);
	for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){
		sph_system[i] = ptcl[i];
	}
	/////////
	sysinfo.end_time = 1.0e-2;
	//Fin.
	std::cout << "setup..." << std::endl;
}

#else
#error NO
#endif

template <class ThisPtcl> void SetEoS(PS::ParticleSystem<ThisPtcl>& sph_system){
	for(PS::U64 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].setPressure(&Granite);
	}
}

