//#define SANITY_CHECK_REALLOCATABLE_ARRAY
//#define FAST
#define FLAG_TREE_COMB
#include "header.h"

#include "phantomquad.hpp"
#include "force-k.h"
//#include "force.h"
const PS::F64 EXPAND = 1.2;
//const PS::F64 EXPAND = 3.0;

int main(int argc, char* argv[]){
	//////////////////
	//Create vars.
	//////////////////
	PS::Initialize(argc, argv);
	PS::ParticleSystem<RealPtcl> sph_system;
	sph_system.initialize();
	sph_system.setAverageTargetNumberOfSampleParticlePerProcess(100);
	PS::DomainInfo dinfo;
	dinfo.initialize();
	system_t sysinfo;
	//timer
	PS::Timer density_timer, derivative_timer, hydro_timer, gravity_timer;
	density_timer.reset();
	density_timer.start();
	derivative_timer.reset();
	derivative_timer.start();
	hydro_timer.reset();
	hydro_timer.start();
	gravity_timer.reset();
	gravity_timer.start();

	std::ofstream density_fout_tcal;
	char density_timerfile[256];
	sprintf(density_timerfile, "./result/density_time_%04d.dat", PS::Comm::getRank());
	density_fout_tcal.open(density_timerfile);

	std::ofstream derivative_fout_tcal;
	char derivative_timerfile[256];
	sprintf(derivative_timerfile, "./result/derivative_time_%04d.dat", PS::Comm::getRank());
	derivative_fout_tcal.open(derivative_timerfile);

	std::ofstream hydro_fout_tcal;
	char hydro_timerfile[256];
	sprintf(hydro_timerfile, "./result/hydro_time_%04d.dat", PS::Comm::getRank());
	hydro_fout_tcal.open(hydro_timerfile);

	std::ofstream gravity_fout_tcal;
	char gravity_timerfile[256];
	sprintf(gravity_timerfile, "./result/gravity_time_%04d.dat", PS::Comm::getRank());
	gravity_fout_tcal.open(gravity_timerfile);

	//////////////////
	//Disp. Info
	//////////////////
	DisplayInfo();
	//Timer
	PS::F64 TbeginInit = PS::GetWtime();
	//////////////////
	//Setup Initial
	//////////////////
	SetupIC(sph_system, sysinfo, dinfo);
	std::cout << "Init." << TbeginInit - PS::GetWtime() << std::endl;
	Initialize(sph_system);
	//Dom. info
	std::cout << "decomp. All" << TbeginInit - PS::GetWtime() << std::endl;
	dinfo.decomposeDomainAll(sph_system);
	std::cout << "ex ptcl" <<  TbeginInit - PS::GetWtime() << std::endl;
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
	const int Ngrp = 128;
	std::cout << "Init. dens" << std::endl;
	dens_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, 4, Ngrp);
	std::cout << "Init. drvt" << std::endl;
	drvt_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, 4, Ngrp);
	std::cout << "Init. hydro" << std::endl;
	hydr_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, 4, Ngrp);
	std::cout << "Init. grav" << std::endl;
	grav_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, 8, 256);

	std::cout << "Calc Dens" << TbeginInit - PS::GetWtime() << std::endl;
	#ifdef FAST
	for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].Rsearch = kernel_t::supportRadius() * sph_system[i].smth * EXPAND;
	}
    int cnt = 0;
    for(bool repeat = true; repeat==true;){
        std::cout<<"cnt="<<cnt<<std::endl;
        bool repeat_loc = false;
        repeat = false;
        //dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
        dens_tree.calcForceAll(CalcDensity(), sph_system, dinfo);
        for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
            if( dens_tree.getForce(i).itr == true){
                repeat_loc = true;
                sph_system[i].Rsearch *= EXPAND;
            }
        }
        repeat = PS::Comm::synchronizeConditionalBranchOR(repeat_loc);
        if(!repeat){
            for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
                sph_system[i].copyFromForce(dens_tree.getForce(i));
            }
        }
        cnt++;
    }

    for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
        sph_system[i].Rsearch = sph_system[i].smth * kernel_t::supportRadius();
    }

	#else
	for(int loop = 0 ; loop <= 2 ; ++ loop){
		dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
		std::cout << loop << std::endl;
	}
	#endif




	std::cout << "Calc Pres" << std::endl;
	CalcPressure(sph_system);
	std::cout << "Calc Drvt" << TbeginInit - PS::GetWtime() << std::endl;
	#ifdef FLAG_TREE_COMB
	dens_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo, false);
	#else
	drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo);
	std::cout << "Calc Hydro" << TbeginInit - PS::GetWtime() << std::endl;
	hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
	std::cout << "Calc Grav" << TbeginInit - PS::GetWtime() << std::endl;
	grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);

	
	std::cout << "get dt" << TbeginInit - PS::GetWtime() << std::endl;
	sysinfo.dt = getTimeStepGlobal(sph_system);
	std::cout << "calc EXT" << TbeginInit - PS::GetWtime() << std::endl;
	CalcExternalForce(sph_system, sysinfo);

	std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt << std::endl;

	PS::S32 step = 0;
	for(sysinfo.time = 0 ; sysinfo.time < sysinfo.end_time ; sysinfo.time += sysinfo.dt, ++ step){
		PS::F64 Tbegin = PS::GetWtime();
		InitialKick(sph_system, sysinfo);
		FullDrift(sph_system, sysinfo);
		std::cout << "IK & FD ... " << PS::GetWtime() - Tbegin << std::endl;
		sph_system.adjustPositionIntoRootDomain(dinfo);
		std::cout << "AdjustPos... " << PS::GetWtime() - Tbegin << std::endl;
		Predict(sph_system, sysinfo);
		std::cout << "Predict ... " << PS::GetWtime() - Tbegin << std::endl;
		dinfo.collectSampleParticle(sph_system);
		dinfo.decomposeDomainMultiStep();
		std::cout << "decomose Dom. all ... " << PS::GetWtime() - Tbegin << std::endl;
		sph_system.exchangeParticle(dinfo);
		std::cout << "Ex.Ptcl ... " << PS::GetWtime() - Tbegin << std::endl;
		density_timer.reset();
		density_timer.start();
#ifdef FAST
        for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
            sph_system[i].Rsearch = kernel_t::supportRadius() * sph_system[i].smth * EXPAND;
        }
        int cnt = 0;
        for(bool repeat = true; repeat==true;){
            //std::cout<<"cnt="<<cnt<<std::endl;
            bool repeat_loc = false;
            repeat = false;
            //dens_tree.calcForceAll(CalcDensity(), sph_system, dinfo);
            dens_tree.calcForceAllWithTimer(CalcDensity(), sph_system, dinfo, density_timer, true);
            for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
                if( dens_tree.getForce(i).itr == true){
                    repeat_loc = true;
                    sph_system[i].Rsearch *= EXPAND;
                }
            }
            repeat = PS::Comm::synchronizeConditionalBranchOR(repeat_loc);
            if(!repeat){
                for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
                    sph_system[i].copyFromForce(dens_tree.getForce(i));
                }
            }
            cnt++;
        }
        for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
            sph_system[i].Rsearch = sph_system[i].smth * kernel_t::supportRadius();
        }
#else

		for(int loop = 0 ; loop <= 2 ; ++ loop){
			dens_tree.calcForceAllAndWriteBackWithTimer(CalcDensity(), sph_system, dinfo, density_timer, true);
		}
#endif
		std::cout << "calc Dens... ... " << PS::GetWtime() - Tbegin << std::endl;
		dens_tree.dump_calc_cost(1.0, density_fout_tcal);
		density_timer.dump(density_fout_tcal);




		CalcPressure(sph_system);
		std::cout << "set pres " << PS::GetWtime() - Tbegin << std::endl;
		derivative_timer.reset();
		derivative_timer.start();
		drvt_tree.calcForceAllAndWriteBackWithTimer(CalcDerivative(), sph_system, dinfo, derivative_timer, true);
		drvt_tree.dump_calc_cost(1.0, derivative_fout_tcal);
		derivative_timer.dump(derivative_fout_tcal);
		#endif
		std::cout << "calc Deriv ... " << PS::GetWtime() - Tbegin << std::endl;
		hydro_timer.reset();
		hydro_timer.start();
		hydr_tree.calcForceAllAndWriteBackWithTimer(CalcHydroForce(), sph_system, dinfo, hydro_timer, true);
		hydr_tree.dump_calc_cost(1.0, hydro_fout_tcal);
		hydro_timer.dump(hydro_fout_tcal);
		std::cout << "calc Hydro ... " << PS::GetWtime() - Tbegin << std::endl;
		gravity_timer.reset();
		gravity_timer.start();
		grav_tree.calcForceAllAndWriteBackWithTimer(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo, gravity_timer, true);
		grav_tree.dump_calc_cost(1.0, gravity_fout_tcal);
		gravity_timer.dump(gravity_fout_tcal);
		std::cout << "calc Grav ... " << PS::GetWtime() - Tbegin << std::endl;

		sysinfo.dt = getTimeStepGlobal(sph_system);
		std::cout << "get dt ... " << PS::GetWtime() - Tbegin << std::endl;
		CalcExternalForce(sph_system, sysinfo);

		FinalKick(sph_system, sysinfo);
		std::cout << "Ext & FK ... " << PS::GetWtime() - Tbegin << std::endl;
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
		std::cout << "Output ... " << PS::GetWtime() - Tbegin << std::endl;
		if(PS::Comm::getRank() == 0){
			std::cout << "//================================" << std::endl;
			std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt << std::endl;
			std::cout << "step = " << step << std::endl;
			std::cout << "//================================" << std::endl;
		}
		//CheckConservativeVariables(sph_system);
		ShiftOrigin(sph_system);
		std::cout << "Shift Origin ... " << PS::GetWtime() - Tbegin << std::endl;
		std::cout << "//================================" << std::endl;
	}

	PS::Finalize();
	return 0;
}

