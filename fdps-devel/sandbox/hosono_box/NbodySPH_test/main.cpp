//#define SANITY_CHECK_REALLOCATABLE_ARRAY
#include "header.h"

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
	//timer
	PS::Timer timer;
	timer.reset();
	timer.start();
	
	std::ofstream fout_tcal;
	char timerfile[256];
	sprintf(timerfile, "time_%04d.dat", PS::Comm::getRank());
	fout_tcal.open(timerfile);
	//////////////////
	//Disp. Info
	//////////////////
	DisplayInfo();
	//////////////////
	//Setup Initial
	//////////////////
	SetupIC(sph_system, sysinfo, dinfo);
	std::cout << "Init." << std::endl;
	Initialize(sph_system);
	//Dom. info
	std::cout << "decomp. All" << std::endl;
	dinfo.decomposeDomainAll(sph_system);
	std::cout << "ex ptcl" << std::endl;
	sph_system.exchangeParticle(dinfo);
	//plant tree
	std::cout << "dens tree" << std::endl;
	PS::TreeForForceShort<RESULT::Dens , EPI::Dens , EPJ::Dens >::Gather   dens_tree;
	std::cout << "drvt tree" << std::endl;
	PS::TreeForForceShort<RESULT::Drvt , EPI::Drvt , EPJ::Drvt >::Gather   drvt_tree;
	std::cout << "hydro tree" << std::endl;
	PS::TreeForForceShort<RESULT::Hydro, EPI::Hydro, EPJ::Hydro>::Symmetry hydr_tree;
	std::cout << "grav tree" << std::endl;
	PS::TreeForForceLong <RESULT::Grav , EPI::Grav , EPJ::Grav >::Monopole grav_tree;

	std::cout << "Init. dens" << std::endl;
	dens_tree.initialize(sph_system.getNumberOfParticleGlobal());
	std::cout << "Init. drvt" << std::endl;
	drvt_tree.initialize(sph_system.getNumberOfParticleGlobal());
	std::cout << "Init. hydro" << std::endl;
	hydr_tree.initialize(sph_system.getNumberOfParticleGlobal());
	std::cout << "Init. grav" << std::endl;
	grav_tree.initialize(sph_system.getNumberOfParticleGlobal());

	std::cout << "Calc Dens" << std::endl;
	for(int loop = 0 ; loop <= 2 ; ++ loop){
		dens_tree.calcForceAllAndWriteBackWithTimer(CalcDensity(), sph_system, dinfo, timer, true);
		dens_tree.dump_calc_cost(1.0, fout_tcal);
		std::cout << loop << std::endl;
	}
	timer.dump(fout_tcal);

	std::cout << "Calc Pres" << std::endl;
	CalcPressure(sph_system);
	std::cout << "Calc Drvt" << std::endl;
	drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo);
	std::cout << "Calc Hydro" << std::endl;
	hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
	std::cout << "Calc Grav" << std::endl;
	grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);
	
	std::cout << "get dt" << std::endl;
	sysinfo.dt = getTimeStepGlobal(sph_system);
	std::cout << "calc EXT" << std::endl;
	CalcExternalForce(sph_system, sysinfo);

	std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt << std::endl;

	PS::S32 step = 0;
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
		grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);

		sysinfo.dt = getTimeStepGlobal(sph_system);
		CalcExternalForce(sph_system, sysinfo);

		FinalKick(sph_system, sysinfo);
		if(step % PARAM::OUTPUT_INTERVAL == 0){
			FileHeader header;
			header.time = sysinfo.time;
			header.Nbody = sph_system.getNumberOfParticleLocal();
			char filename[256];
			sprintf(filename, "result/%05d", step / PARAM::OUTPUT_INTERVAL);
			sph_system.writeParticleAscii(filename, "%s_%05d_%05d.dat", header);
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
		ShiftOrigin(sph_system);
	}

	PS::Finalize();
	return 0;
}

