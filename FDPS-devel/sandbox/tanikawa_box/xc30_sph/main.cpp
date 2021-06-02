//#define SANITY_CHECK_REALLOCATABLE_ARRAY
#define FAST
#include "header.h"
#include "force.h"
//#include "force-k.h"
//#include "force-k2.h"
const PS::F64 EXPAND = 1.3;
//const PS::F64 EXPAND = 1.5;

int main(int argc, char* argv[]){

    PS::F64 flops_per_proc = 1.28e11; // for K computer (8thread)
    PS::F64 n_op_dens = 42.0*3.0 + 74.0;
    PS::F64 n_op_hydr = 119.0;
    PS::F64 n_op_ep_ep_grav = 28.0;
    PS::F64 n_op_ep_sp_grav = 28.0;
    PS::F64 wtime_dens = 0.0;
    PS::F64 wtime_hydr = 0.0;
    PS::F64 wtime_grav = 0.0;

	//////////////////
	//Create vars.
	//////////////////
	PS::Initialize(argc, argv);
	PS::ParticleSystem<RealPtcl> sph_system;
	sph_system.initialize();
	//sph_system.setAverageTargetNumberOfSampleParticlePerProcess(100);
    sph_system.setAverageTargetNumberOfSampleParticlePerProcess(200);
	PS::DomainInfo dinfo;
    const PS::F32 coef_ema = 0.2;
	dinfo.initialize(coef_ema);
	system_t sysinfo;
	//timer
    PS::Timer timer;
	std::ofstream fout_tcal;
	char timerfile[256];
	sprintf(timerfile, "./result/time_%04d.dat", PS::Comm::getRank());
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
	std::cout << "hydro tree" << std::endl;
	PS::TreeForForceShort<RESULT::Hydro, EPI::Hydro, EPJ::Hydro>::Symmetry hydr_tree;
	std::cout << "grav tree" << std::endl;
	PS::TreeForForceLong <RESULT::Grav , EPI::Grav , EPJ::Grav >::Monopole grav_tree;
	const int Ngrp = 128;
	std::cout << "Init. dens" << std::endl;
	dens_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, 4, Ngrp);

    std::cout<<"dens_tree.getMemSize()="<<dens_tree.getMemSizeUsed()*1e-9<<std::endl;


	std::cout << "Init. hydro" << std::endl;
	hydr_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, 4, Ngrp);

    std::cout<<"hydr_tree.getMemSize()="<<hydr_tree.getMemSizeUsed()*1e-9<<std::endl;

	std::cout << "Init. grav" << std::endl;
	grav_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, 8, 256);

    std::cout<<"grav_tree.getMemSize()="<<grav_tree.getMemSizeUsed()*1e-9<<std::endl;

	std::cout << "Calc Dens" << std::endl;

	for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].Rsearch = kernel_t::supportRadius() * sph_system[i].smth * EXPAND;
	}
    int cnt = 0;
#if 0
    for(bool repeat = true; repeat==true;){
        std::cout<<"cnt="<<cnt<<std::endl;
        cnt++;
        bool repeat_loc = false;
        repeat = false;
        dens_tree.calcForceAllWithTimer(CalcDensity(), sph_system, dinfo, timer, true);
        for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
            if(dens_tree.getForce(i).itr == true){
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

    }
    CalcPressure(sph_system);
    for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
        //sph_system[i].copyFromForce(dens_tree.getForce(i));
        sph_system[i].setBalsalaSwitch();
        sph_system[i].Rsearch = sph_system[i].smth * kernel_t::supportRadius();
    }
#else
    for(bool repeat = true; repeat==true;){
        std::cout<<"cnt="<<cnt<<std::endl;
        cnt++;
        bool repeat_loc = false;
        repeat = false;
        dens_tree.calcForceAllWithTimer(CalcDensity(), sph_system, dinfo, timer, true);
        for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
            if(sph_system[i].Rsearch != 0.0){
                if(dens_tree.getForce(i).itr == true){
                    repeat_loc = true;
                    sph_system[i].Rsearch *= EXPAND;
                }
                else{
                    sph_system[i].Rsearch = 0.0;
                    sph_system[i].copyFromForce(dens_tree.getForce(i));
                }
            }
        }
        repeat = PS::Comm::synchronizeConditionalBranchOR(repeat_loc);
    }
    CalcPressure(sph_system);
    for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
        //sph_system[i].copyFromForce(dens_tree.getForce(i));
        sph_system[i].setBalsalaSwitch();
        sph_system[i].Rsearch = sph_system[i].smth * kernel_t::supportRadius();
    }
#endif

    std::cout<<"dens_tree.getMemSize()="<<dens_tree.getMemSizeUsed()*1e-9<<std::endl;

/*
    for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
        sph_system[i].Rsearch = sph_system[i].smth * kernel_t::supportRadius();
    }
*/

	//std::cout << "Calc Pres" << std::endl;
	//CalcPressure(sph_system);

	std::cout << "Calc Hydro" << std::endl;
	hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);

    std::cout<<"hydr_tree.getMemSize()="<<hydr_tree.getMemSizeUsed()*1e-9<<std::endl;

	std::cout << "Calc Grav" << std::endl;
	grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);

    std::cout<<"grav_tree.getMemSize()="<<grav_tree.getMemSizeUsed()*1e-9<<std::endl;
	
	std::cout << "get dt" << std::endl;
	sysinfo.dt = getTimeStepGlobal(sph_system);
	std::cout << "calc EXT" << std::endl;
	CalcExternalForce(sph_system, sysinfo);

	std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt << std::endl;

/*
    if(1){
        FileHeader header;
        header.time = sysinfo.time;
        header.Nbody = sph_system.getNumberOfParticleLocal();
        char filename[256];
        sprintf(filename, "result/%s_%05d.dat", "init", PS::Comm::getRank());
        std::ofstream fout;
        fout.open(filename);
        fout<<std::setprecision(15);
        for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
            sph_system[i].dump(fout);
        }
        fout.close();
    }

    sph_system[0].dens = 2998.02878243082;
    sph_system[0].eng = 6253125;
    sph_system[1].dens = 2998.02878243083;
    sph_system[1].eng = 6253125;
    sph_system[2].dens = 2998.02878243082;
    sph_system[2].eng = 6253125;
    sph_system[3].dens = 2998;
    sph_system[3].eng = 6253125;

    CalcPressure(sph_system);

    if(PS::Comm::getRank() == 0){
        std::cout<<"sph_system[0].snds="<<sph_system[0].snds<<std::endl;
        std::cout<<"sph_system[1].snds="<<sph_system[1].snds<<std::endl;
        std::cout<<"sph_system[2].snds="<<sph_system[2].snds<<std::endl;
        std::cout<<"sph_system[3].snds="<<sph_system[3].snds<<std::endl;
        std::cout<<"sph_system[0].pres="<<sph_system[0].pres<<std::endl;
        std::cout<<"sph_system[1].pres="<<sph_system[1].pres<<std::endl;
        std::cout<<"sph_system[2].pres="<<sph_system[2].pres<<std::endl;
        std::cout<<"sph_system[3].pres="<<sph_system[3].pres<<std::endl;
    }
*/

    if(1){
        FileHeader header;
        header.time = sysinfo.time;
        header.Nbody = sph_system.getNumberOfParticleLocal();
        char filename[256];
        sprintf(filename, "result/%s", "init");
        sph_system.writeParticleAscii(filename, "%s_%05d_%05d.dat", header);
        if(PS::Comm::getRank() == 0){
            std::cout << "//================================" << std::endl;
            std::cout << "output " << filename << "." << std::endl;
            std::cout << "//================================" << std::endl;
        }
    }

    //PS::Finalize();
    //return 0;



	PS::S64 step = 0;
    PS::F64 n_op_tot = 1.0;
	for(sysinfo.time = 0 ; sysinfo.time < sysinfo.end_time && step < 64; sysinfo.time += sysinfo.dt, ++ step){

        timer.reset();
        timer.start();

		//PS::F64 Tbegin = PS::GetWtime();
        PS::F64 wtime_tot = PS::GetWtime();
		InitialKick(sph_system, sysinfo);

        timer.restart("InitialKick");

		FullDrift(sph_system, sysinfo);

        timer.restart("FullDrift");

		Predict(sph_system, sysinfo);

        timer.restart("Predict");

        if(step < 2){
            dinfo.collectSampleParticle(sph_system);
            timer.restart("collectSampleParticle");
            dinfo.decomposeDomainMultiStep();
        }
        else if(step < 12 || step % 4 == 0){
            //dinfo.collectSampleParticle(sph_system, true, wtime_dens+wtime_hydr+wtime_grav);
            dinfo.collectSampleParticle(sph_system, true, n_op_tot);
            timer.restart("collectSampleParticle");
            dinfo.decomposeDomainMultiStep();
        }
        else{
            timer.restart("collectSampleParticle");

        }
        timer.restart("decomposeDomainMultiStep");

		sph_system.exchangeParticle(dinfo);

        timer.restart("exchangeParticle");

        for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
            sph_system[i].Rsearch = kernel_t::supportRadius() * sph_system[i].smth * EXPAND;
        }

        timer.restart("setRsearch");

/*
        int cnt = 0;
        for(bool repeat = true; repeat==true;){
            std::cout<<"cnt="<<cnt<<std::endl;
            cnt++;
            bool repeat_loc = false;
            repeat = false;
            dens_tree.calcForceAllWithTimer(CalcDensity(), sph_system, dinfo, timer, true);
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
        }
*/

        PS::S64 n_int_ep_ep = 0;
        wtime_dens = PS::GetWtime();
#if 0
        cnt = 0;
        for(bool repeat = true; repeat==true;){

            std::cout<<"cnt="<<cnt<<std::endl;
            cnt++;
            bool repeat_loc = false;
            repeat = false;
            dens_tree.calcForceAllWithTimer(CalcDensity(), sph_system, dinfo, timer, true);
            for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
                if(dens_tree.getForce(i).itr == true){
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
            n_int_ep_ep += dens_tree.getNInteractionEPEP();
        }
/*
        CalcPressure(sph_system);
        for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
            //sph_system[i].copyFromForce(dens_tree.getForce(i));
            sph_system[i].setBalsalaSwitch();
            sph_system[i].Rsearch = sph_system[i].smth * kernel_t::supportRadius();
        }
*/
#else
        cnt = 0;
        for(bool repeat = true; repeat==true;){
            //std::cout<<"cnt="<<cnt<<std::endl;

            cnt++;
            bool repeat_loc = false;
            repeat = false;
            dens_tree.calcForceAllWithTimer(CalcDensity(), sph_system, dinfo, timer, true);
            for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
                if(sph_system[i].Rsearch != 0.0){
                    if(dens_tree.getForce(i).itr == true){
                        repeat_loc = true;
                        sph_system[i].Rsearch *= EXPAND;
                    }
                    else{
                        sph_system[i].Rsearch = 0.0;
                        sph_system[i].copyFromForce(dens_tree.getForce(i));
                    }
                }
            }
            repeat = PS::Comm::synchronizeConditionalBranchOR(repeat_loc);
            n_int_ep_ep += dens_tree.getNInteractionEPEP();
        }
#endif
        dens_tree.setNInteractionEPEP(n_int_ep_ep);
        wtime_dens = PS::GetWtime() - wtime_dens;
        timer.restart("misc_and_copy");


		CalcPressure(sph_system);

        timer.restart("CalcPressure");

        for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
            //sph_system[i].copyFromForce(dens_tree.getForce(i));
            sph_system[i].setBalsalaSwitch();
            sph_system[i].Rsearch = sph_system[i].smth * kernel_t::supportRadius();
        }

        timer.restart("setRsearch");

        wtime_hydr = PS::GetWtime();
		hydr_tree.calcForceAllAndWriteBackWithTimer(CalcHydroForce(), sph_system, dinfo, timer, true);
        wtime_hydr = PS::GetWtime() - wtime_hydr;

        timer.restart("copyHydroForce");

        wtime_grav = PS::GetWtime();
		grav_tree.calcForceAllAndWriteBackWithTimer(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo, timer, true);
        wtime_grav = PS::GetWtime() - wtime_grav;

        timer.restart("copyGravityForce");

		sysinfo.dt = getTimeStepGlobal(sph_system);

        timer.restart("setDt");

		CalcExternalForce(sph_system, sysinfo);

        timer.restart("CalcExternalForce");

		FinalKick(sph_system, sysinfo);

        timer.restart("FinalKick");

        wtime_tot = PS::GetWtime() - wtime_tot;

		//dens_tree.dump_calc_cost(1.0, fout_tcal);
        PS::S64 n_int_dens_glb = PS::Comm::getSum(dens_tree.getNInteractionEPEP());
        PS::S64 n_int_hydr_glb = PS::Comm::getSum(hydr_tree.getNInteractionEPEP());
        PS::S64 n_int_grav_ep_ep_glb = PS::Comm::getSum(grav_tree.getNInteractionEPEP());
        PS::S64 n_int_grav_ep_sp_glb = PS::Comm::getSum(grav_tree.getNInteractionEPSP());
        //PS::F64 n_op_tot = (PS::F64)()*n_int_dens + (PS::F64)(hydr_tree.getNInteractionEPEP())*n_int_hydr + (PS::F64)(grav_tree.getNInteractionEPEP())*n_int_grav;
        n_op_tot = (PS::F64)(n_int_dens_glb)*n_op_dens + (PS::F64)(n_int_hydr_glb)*n_op_hydr 
            + (PS::F64)(n_int_grav_ep_ep_glb)*n_op_ep_ep_grav + (PS::F64)(n_int_grav_ep_sp_glb)*n_op_ep_sp_grav;
        PS::F64 total_speed = (PS::F64)(n_op_tot)/wtime_tot;
        PS::F64 total_speed_max = flops_per_proc*((PS::F64)(PS::Comm::getNumberOfProc()));
        fout_tcal<<"total_speed="<<total_speed*1e-12<<"[Tflops] efficiency="<<total_speed/total_speed_max<<" n_op_tot="<<n_op_tot<<" wtime_tot="<<wtime_tot<<std::endl;
/*
                 <<" (PS::F64)(dens_tree.getNInteractionEPEP())="<<(PS::F64)(dens_tree.getNInteractionEPEP())
                 <<" (PS::F64)(dens_tree.getNInteractionEPEP())*n_int_dens="<<(PS::F64)(dens_tree.getNInteractionEPEP())*n_int_dens<<std::endl;
*/
        dens_tree.dump_calc_cost(wtime_dens, fout_tcal, flops_per_proc, n_op_dens, 0.0);
		hydr_tree.dump_calc_cost(wtime_hydr, fout_tcal, flops_per_proc, n_op_hydr, 0.0);
		grav_tree.dump_calc_cost(wtime_grav, fout_tcal, flops_per_proc, n_op_ep_ep_grav, n_op_ep_sp_grav);
		timer.dump(fout_tcal);
        fout_tcal<<std::endl;

/*
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
*/
		//CheckConservativeVariables(sph_system);
		ShiftOrigin(sph_system);
		//std::cout << "Shift Origin ... " << PS::GetWtime() - Tbegin << std::endl;
		//std::cout << "//================================" << std::endl;
	}

	PS::Finalize();
	return 0;
}

