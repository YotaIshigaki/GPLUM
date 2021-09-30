//#define SANITY_CHECK_REALLOCATABLE_ARRAY
#include "header.h"
// #include "force.h"
// #include "force-k.h"
#include "force-k-sp.h"

int main(int argc, char* argv[]){
	//////////////////
	//Create vars.
	//////////////////
	PS::Initialize(argc, argv);
	PS::ParticleSystem<RealPtcl> sph_system;
	sph_system.initialize();
	PS::DomainInfo dinfo;
	dinfo.initialize();
	system_t sysinfo;

	//////////////////
	//Disp. Info
	//////////////////
	DisplayInfo();
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
	PS::TreeForForceShort<RESULT::Drvt , EPI::Drvt , EPJ::Drvt >::Gather   drvt_tree;
	PS::TreeForForceShort<RESULT::Hydro, EPI::Hydro, EPJ::Hydro>::Symmetry hydr_tree;
	PS::TreeForForceLong <RESULT::Grav , EPI::Grav , EPJ::Grav >::Monopole grav_tree;

	dens_tree.initialize(sph_system.getNumberOfParticleGlobal());
	drvt_tree.initialize(sph_system.getNumberOfParticleGlobal());
	hydr_tree.initialize(sph_system.getNumberOfParticleGlobal());
	grav_tree.initialize(sph_system.getNumberOfParticleGlobal());

	for(int loop = 0 ; loop <= 2 ; ++ loop){
		dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
	}

	CalcPressure(sph_system);
	drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo);
	hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
	if(sysinfo.Grav > 0.0){
		grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);
	}
	sysinfo.dt = getTimeStepGlobal(sph_system);
	CalcExternalForce(sph_system, sysinfo);

	std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt << std::endl;

	PS::S32 step = 0;
	if(step % PARAM::OUTPUT_INTERVAL == 0){
		FileHeader header;
		header.time = sysinfo.time;
		header.Nbody = sph_system.getNumberOfParticleGlobal();
		char filename[256];
		sprintf(filename, "result/%04d.dat", step / PARAM::OUTPUT_INTERVAL);
		sph_system.writeParticleAscii(filename, header);
		if(PS::Comm::getRank() == 0){
			std::cout << "//================================" << std::endl;
			std::cout << "output " << filename << "." << std::endl;
			std::cout << "//================================" << std::endl;
		}
	}

	for(sysinfo.time = 0 ; sysinfo.time < sysinfo.end_time ; sysinfo.time += sysinfo.dt, ++ step){
		InitialKick(sph_system, sysinfo);
		FullDrift(sph_system, sysinfo);
		sph_system.adjustPositionIntoRootDomain(dinfo);
		Predict(sph_system, sysinfo);
		dinfo.decomposeDomainAll(sph_system);
		sph_system.exchangeParticle(dinfo);
		for(int loop = 0 ; loop <= 2 ; ++ loop){
			dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
		}
		CalcPressure(sph_system);
		drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo);
		hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
		if(sysinfo.Grav > 0.0){
			grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);
		}
		sysinfo.dt = getTimeStepGlobal(sph_system);
		CalcExternalForce(sph_system, sysinfo);

		FinalKick(sph_system, sysinfo);
		if(step % PARAM::OUTPUT_INTERVAL == 0){
			FileHeader header;
			header.time = sysinfo.time;
			header.Nbody = sph_system.getNumberOfParticleLocal();
			char filename[256];
			sprintf(filename, "result/%04d.dat", step / PARAM::OUTPUT_INTERVAL);
			sph_system.writeParticleAscii(filename, header);
			if(PS::Comm::getRank() == 0){
				std::cout << "//================================" << std::endl;
				std::cout << "output " << filename << "." << std::endl;
				std::cout << "//================================" << std::endl;
			}
		}

		if(PS::Comm::getRank() == 0){
			std::cout << "//================================" << std::endl;
			std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt << std::endl;
			std::cout << "step = " << step << std::endl;
			std::cout << "//================================" << std::endl;
		}
		CheckConservativeVariables(sph_system);
	}

	PS::Finalize();
	return 0;
}

