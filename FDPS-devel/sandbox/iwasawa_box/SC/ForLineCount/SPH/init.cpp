#include "header.h"
void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, system_t& sysinfo, PS::DomainInfo& dinfo){
    std::vector<RealPtcl> ptcl;
    std::vector<RealPtcl> tar;//Target
    std::vector<RealPtcl> imp;//Impactor
    const PS::F64 EarthMass = 6.0e+24;
    const PS::F64 EarthRadi = 6400e+3;
    const PS::F64 tarMass = EarthMass;
    const PS::F64 tarRadi = EarthRadi;
    const PS::F64 impMass = 0.1 * EarthMass;
    const PS::F64 impRadi = 0.5 * EarthRadi;
    sysinfo.Grav = 6.67e-11;
    sysinfo.end_time = 1.0e+4;
    const int Nx = 160;
    const PS::F64 dx = 1.0 / (double)(Nx);// for 64 nodes
    int Nptcl = 0;
    for(PS::F64 x = -1.0 ; x <= 1.0 ; x += dx){
	for(PS::F64 y = -1.0 ; y <= 1.0 ; y += dx){
	    for(PS::F64 z = -1.0 ; z <= 1.0 ; z += dx){
		const PS::F64 r = sqrt(x*x + y*y + z*z);
		if(r > 1.0) continue;
		++ Nptcl;
	    }
	}
    }
    for(PS::F64 x = -1.0 ; x <= 1.0 ; x += 2.0 * dx){
	for(PS::F64 y = -1.0 ; y <= 1.0 ; y += 2.0 * dx){
	    for(PS::F64 z = -1.0 ; z <= 1.0 ; z += 2.0 * dx){
		const PS::F64 r = sqrt(x*x + y*y + z*z);
		if(r > 1.0) continue;
		++ Nptcl;
	    }
	}
    }
    const int NptclIn1Node = Nptcl / PS::Comm::getNumberOfProc();
    PS::S32 id = 0;
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
		if(ith.id / NptclIn1Node == PS::Comm::getRank()) tar.push_back(ith);
	    }
	}
    }
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
		if(ith.id / NptclIn1Node == PS::Comm::getRank()) imp.push_back(ith);
	    }
	}
    }
    for(PS::U32 i = 0 ; i < tar.size() ; ++ i){	tar[i].mass /= (double)(Nptcl); }
    for(PS::U32 i = 0 ; i < imp.size() ; ++ i){ imp[i].mass /= (double)(Nptcl); }
    for(PS::U32 i = 0 ; i < tar.size() ; ++ i){	ptcl.push_back(tar[i]);}	
    for(PS::U32 i = 0 ; i < imp.size() ; ++ i){ptcl.push_back(imp[i]);}
    const PS::S32 numPtclLocal = ptcl.size();
    sph_system.setNumberOfParticleLocal(numPtclLocal);
    for(PS::U32 i = 0 ; i < ptcl.size() ; ++ i){sph_system[i] = ptcl[i];}
}
void Initialize(PS::ParticleSystem<RealPtcl>& sph_system){
#pragma omp parallel for
    for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
	sph_system[i].smth = PARAM::SMTH * pow(sph_system[i].mass / sph_system[i].dens, 1.0/(PS::F64)(PARAM::Dim));
    }
}
