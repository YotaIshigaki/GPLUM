#include "header.h"

namespace CalcDensity{
	int DispatchKernel(const PS::S32, const int, const EPI::Dens**, const int*, const EPJ::Dens**, const int*);
	int RetrieveKernel(const PS::S32, const PS::S32, const PS::S32*, RESULT::Dens**);
}

namespace CalcHydroForce{
	int DispatchKernel(const PS::S32, const int, const EPI::Hydro**, const int*, const EPJ::Hydro**, const int*);
	int RetrieveKernel(const PS::S32, const PS::S32, const PS::S32*, RESULT::Hydro**);
}

void OutputFileWithTimeInterval(PS::ParticleSystem<RealPtcl>&, const system_t&);
void InputFileWithTimeInterval(PS::ParticleSystem<RealPtcl>&, system_t&);

PS::F32 getTimeStepGlobal(const PS::ParticleSystem<RealPtcl>&);
void InitialKick(PS::ParticleSystem<RealPtcl>&, const system_t);
void FullDrift(PS::ParticleSystem<RealPtcl>&, const system_t);
void Predict(PS::ParticleSystem<RealPtcl>&, const system_t);
void FinalKick(PS::ParticleSystem<RealPtcl>&, const system_t);

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
	sysinfo.end_time = 0.2;
	//Fin.
	std::cout << "setup..." << std::endl;
}

void Initialize(PS::ParticleSystem<RealPtcl>& sph_system){
	#pragma omp parallel for
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].smth = PARAM::SMTH * pow(sph_system[i].mass / sph_system[i].dens, 1.0/(PS::F64)(PARAM::Dim));
		sph_system[i].grad_smth = 1.0;
	}
}

void SetPressure(PS::ParticleSystem<RealPtcl>& sph_system){
	const double hcr = 1.4;
	#pragma omp parallel for
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].pres = (hcr - 1.0) * sph_system[i].dens * sph_system[i].eng;
		sph_system[i].snds = sqrt(hcr * sph_system[i].pres / sph_system[i].dens);
		sph_system[i].grad_smth = 1.0;
	}
}


int main(int argc, char* argv[]){
	const PS::S64 n_walk_limit  = 64;

	PS::Initialize(argc, argv);

	PS::ParticleSystem<RealPtcl> sph_system;
	sph_system.initialize();
	PS::DomainInfo dinfo;
	dinfo.initialize();
	system_t sysinfo;

	//////////////////
	//Setup Initial
	//////////////////
	SetupIC(sph_system, sysinfo, dinfo);
	Initialize(sph_system);
	//Dom. info
	//dinfo.setDomain(PS::Comm::getNumberOfProc(), 1, 1);
	dinfo.decomposeDomain();
	sph_system.exchangeParticle(dinfo);
	//plant tree
	PS::TreeForForceShort<RESULT::Dens , EPI::Dens , EPJ::Dens >::Gather   dens_tree;
	//PS::TreeForForceShort<RESULT::Drvt , EPI::Drvt , EPJ::Drvt >::Gather   drvt_tree;
	PS::TreeForForceShort<RESULT::Hydro, EPI::Hydro, EPJ::Hydro>::Symmetry hydr_tree;

	dens_tree.initialize(sph_system.getNumberOfParticleGlobal());
	//drvt_tree.initialize(sph_system.getNumberOfParticleGlobal());
	hydr_tree.initialize(sph_system.getNumberOfParticleGlobal());

	for(int loop = 0 ; loop < 2 ; ++ loop){
		dens_tree.calcForceAllAndWriteBackMultiWalk(CalcDensity::DispatchKernel, CalcDensity::RetrieveKernel, 0, sph_system, dinfo, n_walk_limit, true);
	}
	SetPressure(sph_system);
	hydr_tree.calcForceAllAndWriteBackMultiWalk(CalcHydroForce::DispatchKernel, CalcHydroForce::RetrieveKernel, 0, sph_system, dinfo, n_walk_limit, true);
	sysinfo.dt = getTimeStepGlobal(sph_system);
	OutputFileWithTimeInterval(sph_system, sysinfo);
	if(PS::Comm::getRank() == 0){
		std::cout << "//================================" << std::endl;
		std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt << std::endl;
		std::cout << "step = " << sysinfo.step << std::endl;
		std::cout << "//================================" << std::endl;
	}

	while(sysinfo.time < sysinfo.end_time){
		InitialKick(sph_system, sysinfo);
		FullDrift(sph_system, sysinfo);
		sph_system.adjustPositionIntoRootDomain(dinfo);
		Predict(sph_system, sysinfo);
		dinfo.decomposeDomainAll(sph_system);
		sph_system.exchangeParticle(dinfo);
		OutputFileWithTimeInterval(sph_system, sysinfo);
		for(int loop = 0 ; loop < 2 ; ++ loop){
			dens_tree.calcForceAllAndWriteBackMultiWalk(CalcDensity::DispatchKernel, CalcDensity::RetrieveKernel, 0, sph_system, dinfo, n_walk_limit, true);
		}
		SetPressure(sph_system);
		hydr_tree.calcForceAllAndWriteBackMultiWalk(CalcHydroForce::DispatchKernel, CalcHydroForce::RetrieveKernel, 0, sph_system, dinfo, n_walk_limit, true);
		sysinfo.dt = getTimeStepGlobal(sph_system);
		FinalKick(sph_system, sysinfo);
		OutputFileWithTimeInterval(sph_system, sysinfo);
		sysinfo.time += sysinfo.dt;
		++ sysinfo.step;

		if(PS::Comm::getRank() == 0){
			std::cout << "//================================" << std::endl;
			std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt << std::endl;
			std::cout << "step = " << sysinfo.step << std::endl;
			std::cout << "//================================" << std::endl;
		}
		sysinfo.time += sysinfo.dt;
	}

	return 0;
}

