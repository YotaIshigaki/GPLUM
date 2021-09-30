//#define SANITY_CHECK_REALLOCATABLE_ARRAY

#include <particle_simulator.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "param.h"
#include "mathfunc.h"
#include "kernel.h"
#include "EoS.h"
#include "class.h"
//=============
//#include "init/shock_tube.h"
#include "init/KHI.h"
//#include "init/RTI.h"
//#include "init/conv.h"
//#include "init/GI_singlecomp.h"
//#include "init/GI_Maindl.h"
//#include "init/GI.h"
//#include "init/GI_imp.h"
//#include "init/GoW.h"
//=============
#include "force.h"
#include "test.h"
#include "io.h"
#include "integral.h"

int main(int argc, char* argv[]){
	namespace PTCL = STD;
	//typedef ShockTube<PTCL::RealPtcl> PROBLEM;
	typedef KHI<PTCL::RealPtcl> PROBLEM;
	//typedef RTI<PTCL::RealPtcl> PROBLEM;
	//typedef GI_init<PTCL::RealPtcl> PROBLEM;
	//typedef GI<PTCL::RealPtcl> PROBLEM;
	//typedef Imp<PTCL::RealPtcl> PROBLEM;
	//////////////////
	//Create vars.
	//////////////////
	PS::Initialize(argc, argv);
	PS::ParticleSystem<PTCL::RealPtcl> sph_system;
	sph_system.initialize();
	PS::DomainInfo dinfo;
	dinfo.initialize();
	system_t sysinfo;

	sph_system.setAverageTargetNumberOfSampleParticlePerProcess(200);
	//////////////////
	//Setup Initial
	//////////////////
	if(argc == 1){
		PROBLEM::setupIC(sph_system, sysinfo, dinfo);
		PROBLEM::setEoS(sph_system);
		PTCL::CalcPressure(sph_system);
	}else{
		sysinfo.step = atoi(argv[1]);
		InputFileWithTimeInterval<PTCL::RealPtcl>(sph_system, sysinfo);
		PROBLEM::setEoS(sph_system);
	}

	#pragma omp parallel for
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].initialize();
	}

	//Dom. info
	dinfo.decomposeDomainAll(sph_system);
	sph_system.exchangeParticle(dinfo);
	//plant tree
	PS::TreeForForceShort<PTCL::RESULT::Dens , PTCL::EPI::Dens , PTCL::EPJ::Dens >::Gather   dens_tree;
	PS::TreeForForceShort<PTCL::RESULT::Drvt , PTCL::EPI::Drvt , PTCL::EPJ::Drvt >::Gather   drvt_tree;
	PS::TreeForForceShort<PTCL::RESULT::Hydro, PTCL::EPI::Hydro, PTCL::EPJ::Hydro>::Symmetry hydr_tree;
	#ifdef SELF_GRAVITY
	PS::TreeForForceLong <PTCL::RESULT::Grav , PTCL::EPI::Grav , PTCL::EPJ::Grav >::Monopole grav_tree;
	#endif

	dens_tree.initialize(sph_system.getNumberOfParticleLocal(), 0.5, 1, 1);
	drvt_tree.initialize(sph_system.getNumberOfParticleLocal(), 0.5, 1, 1);
	hydr_tree.initialize(sph_system.getNumberOfParticleLocal(), 0.5, 1, 1);
	#ifdef SELF_GRAVITY
	grav_tree.initialize(sph_system.getNumberOfParticleLocal(), 0.5, 8, 256);
	#endif
	std::cout << "Calc Dens" << std::endl;
	for(short int loop = 0 ; loop <= PARAM::NUMBER_OF_DENSITY_SMOOTHING_LENGTH_LOOP ; ++ loop){
		dens_tree.calcForceAllAndWriteBack(PTCL::CalcDensity(), sph_system, dinfo, PS::MAKE_LIST_FOR_REUSE);
	}

	std::cout << "Calc Pres" << std::endl;
	PTCL::CalcPressure(sph_system);
	std::cout << "Calc Drvt" << std::endl;
	drvt_tree.calcForceAllAndWriteBack(PTCL::CalcDerivative(), sph_system, dinfo, PS::MAKE_LIST_FOR_REUSE);
	std::cout << "Calc Hydro" << std::endl;
	hydr_tree.calcForceAllAndWriteBack(PTCL::CalcHydroForce(), sph_system, dinfo, PS::MAKE_LIST_FOR_REUSE);
	#ifdef SELF_GRAVITY
	std::cout << "Calc Grav" << std::endl;
	grav_tree.calcForceAllAndWriteBack(PTCL::CalcGravityForce<PTCL::EPJ::Grav>(), PTCL::CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo, PS::MAKE_LIST_FOR_REUSE);
	#endif
	PROBLEM::addExternalForce(sph_system, sysinfo);
	OutputFileWithTimeInterval(sph_system, sysinfo, PROBLEM::END_TIME);

	std::cout << "get dt" << std::endl;
	sysinfo.dt = getTimeStepGlobal<PTCL::RealPtcl>(sph_system);

	/////////////
	if(PS::Comm::getRank() == 0){
		std::cout << "//================================" << std::endl;
		std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt << std::endl;
		std::cout << "step = " << sysinfo.step << std::endl;
		std::cout << "//================================" << std::endl;
	}
	/////////////

	while(sysinfo.time < PROBLEM::END_TIME){
		#pragma omp parallel for
		for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].initialKick(sysinfo.dt);
			sph_system[i].fullDrift(sysinfo.dt);
			sph_system[i].predict(sysinfo.dt);
		}
		sysinfo.time += sysinfo.dt;
		sph_system.adjustPositionIntoRootDomain(dinfo);
		#if 1
		if(sysinfo.step % 4 == 0){
			std::cout << "make lists" << std::endl;
			dinfo.decomposeDomainAll(sph_system);
			sph_system.exchangeParticle(dinfo);
			for(short int loop = 0 ; loop <= PARAM::NUMBER_OF_DENSITY_SMOOTHING_LENGTH_LOOP ; ++ loop){
				dens_tree.calcForceAllAndWriteBack(PTCL::CalcDensity(), sph_system, dinfo, PS::MAKE_LIST_FOR_REUSE);
			}
			PTCL::CalcPressure(sph_system);
			drvt_tree.calcForceAllAndWriteBack(PTCL::CalcDerivative(), sph_system, dinfo, PS::MAKE_LIST_FOR_REUSE);
			hydr_tree.calcForceAllAndWriteBack(PTCL::CalcHydroForce(), sph_system, dinfo, PS::MAKE_LIST_FOR_REUSE);
			#ifdef SELF_GRAVITY
			grav_tree.calcForceAllAndWriteBack(PTCL::CalcGravityForce<PTCL::EPJ::Grav>(), PTCL::CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo, PS::MAKE_LIST_FOR_REUSE);
			#endif
		}else{
			std::cout << "reuse lists" << std::endl;
			for(short int loop = 0 ; loop <= PARAM::NUMBER_OF_DENSITY_SMOOTHING_LENGTH_LOOP ; ++ loop){
				dens_tree.calcForceAllAndWriteBack(PTCL::CalcDensity(), sph_system, dinfo, PS::REUSE_LIST);
			}
			PTCL::CalcPressure(sph_system);
			drvt_tree.calcForceAllAndWriteBack(PTCL::CalcDerivative(), sph_system, dinfo, PS::REUSE_LIST);
			hydr_tree.calcForceAllAndWriteBack(PTCL::CalcHydroForce(), sph_system, dinfo, PS::REUSE_LIST);
			#ifdef SELF_GRAVITY
			grav_tree.calcForceAllAndWriteBack(PTCL::CalcGravityForce<PTCL::EPJ::Grav>(), PTCL::CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo, PS::REUSE_LIST);
			#endif
		}
		#else
		std::cout << "make lists" << std::endl;
		dinfo.decomposeDomainAll(sph_system);
		sph_system.exchangeParticle(dinfo);
		for(short int loop = 0 ; loop <= PARAM::NUMBER_OF_DENSITY_SMOOTHING_LENGTH_LOOP ; ++ loop){
			dens_tree.calcForceAllAndWriteBack(PTCL::CalcDensity(), sph_system, dinfo);
		}
		PTCL::CalcPressure(sph_system);
		drvt_tree.calcForceAllAndWriteBack(PTCL::CalcDerivative(), sph_system, dinfo);
		hydr_tree.calcForceAllAndWriteBack(PTCL::CalcHydroForce(), sph_system, dinfo);
		#ifdef SELF_GRAVITY
		grav_tree.calcForceAllAndWriteBack(PTCL::CalcGravityForce<PTCL::EPJ::Grav>(), PTCL::CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);
		#endif
		#endif
		PROBLEM::addExternalForce(sph_system, sysinfo);

		sysinfo.dt = getTimeStepGlobal<PTCL::RealPtcl>(sph_system);

		#pragma omp parallel for
		for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].finalKick(sysinfo.dt);
		}

		PROBLEM::postTimestepProcess(sph_system, sysinfo);
		OutputFileWithTimeInterval<PTCL::RealPtcl>(sph_system, sysinfo, PROBLEM::END_TIME);
		++ sysinfo.step;
		if(PS::Comm::getRank() == 0){
			std::cout << "//================================" << std::endl;
			std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt << std::endl;
			std::cout << "step = " << sysinfo.step << std::endl;
			std::cout << "//================================" << std::endl;
		}
		//CheckConservativeVariables(sph_system);
	}

	PS::Finalize();
	return 0;
}

