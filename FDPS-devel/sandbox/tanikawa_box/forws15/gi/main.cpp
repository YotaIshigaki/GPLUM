//#define SANITY_CHECK_REALLOCATABLE_ARRAY
#define FAST
#include "header.h"
#include "force.h"
//#include "force-k.h"
//#include "force-k2.h"
const PS::F64 EXPAND = 1.3;
//const PS::F64 EXPAND = 1.5;

void DumpTimeProfile(const PS::TimeProfile & tp, const PS::TimeProfile & tp_max, const PS::S32 rank_max[], std::ostream & fout){
    PS::S32 id = 0;
    fout<<"collect_sample_particle= "<<tp.collect_sample_particle<<", max= "<<tp_max.collect_sample_particle<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"decompose_domain= "<<tp.decompose_domain<<", max= "<<tp_max.decompose_domain<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"exchange_particle= "<<tp.exchange_particle<<", max= "<<tp_max.exchange_particle<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"set_particle_local_tree= "<<tp.set_particle_local_tree<<", max= "<<tp_max.set_particle_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"set_root_cell= "<<tp.set_root_cell<<", max= "<<tp_max.set_root_cell<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_local_tree= "<<tp.make_local_tree<<", max= "<<tp_max.make_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"calc_moment_local_tree= "<<tp.calc_moment_local_tree<<", max= "<<tp_max.calc_moment_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_LET_1st= "<<tp.make_LET_1st<<", max= "<<tp_max.make_LET_1st<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"exchange_LET_1st= "<<tp.exchange_LET_1st<<", max= "<<tp_max.exchange_LET_1st<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_LET_2nd= "<<tp.make_LET_2nd<<", max= "<<tp_max.make_LET_2nd<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"exchange_LET_2nd= "<<tp.exchange_LET_2nd<<", max= "<<tp_max.exchange_LET_2nd<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"set_particle_global_tree= "<<tp.set_particle_global_tree<<", max= "<<tp_max.set_particle_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_global_tree= "<<tp.make_global_tree<<", max= "<<tp_max.make_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"calc_moment_global_tree= "<<tp.calc_moment_global_tree<<", max= "<<tp_max.calc_moment_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"calc_force = "<<tp.calc_force<<", max= "<<tp_max.calc_force<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<std::endl;
} 

void DumpTimeProfile0(const PS::TimeProfile & tp, const PS::TimeProfile & tp_max, const PS::S32 rank_max[], std::ostream & fout){
    PS::S32 id = 0;
    fout<<"collect_sample_particle= "<<tp.collect_sample_particle<<", max= "<<tp_max.collect_sample_particle<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"decompose_domain= "<<tp.decompose_domain<<", max= "<<tp_max.decompose_domain<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<std::endl;
} 

void DumpTimeProfile1(const PS::TimeProfile & tp, const PS::TimeProfile & tp_max, const PS::S32 rank_max[], std::ostream & fout){
    PS::S32 id = 2;
    fout<<"exchange_particle= "<<tp.exchange_particle<<", max= "<<tp_max.exchange_particle<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<std::endl;
} 

void DumpTimeProfile2(const PS::TimeProfile & tp, const PS::TimeProfile & tp_max, const PS::S32 rank_max[], std::ostream & fout){
    PS::S32 id = 3;
    fout<<"set_particle_local_tree= "<<tp.set_particle_local_tree<<", max= "<<tp_max.set_particle_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"set_root_cell= "<<tp.set_root_cell<<", max= "<<tp_max.set_root_cell<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_local_tree= "<<tp.make_local_tree<<", max= "<<tp_max.make_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"calc_moment_local_tree= "<<tp.calc_moment_local_tree<<", max= "<<tp_max.calc_moment_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_LET_1st= "<<tp.make_LET_1st<<", max= "<<tp_max.make_LET_1st<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"exchange_LET_1st= "<<tp.exchange_LET_1st<<", max= "<<tp_max.exchange_LET_1st<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_LET_2nd= "<<tp.make_LET_2nd<<", max= "<<tp_max.make_LET_2nd<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"exchange_LET_2nd= "<<tp.exchange_LET_2nd<<", max= "<<tp_max.exchange_LET_2nd<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"set_particle_global_tree= "<<tp.set_particle_global_tree<<", max= "<<tp_max.set_particle_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_global_tree= "<<tp.make_global_tree<<", max= "<<tp_max.make_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"calc_moment_global_tree= "<<tp.calc_moment_global_tree<<", max= "<<tp_max.calc_moment_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"calc_force = "<<tp.calc_force<<", max= "<<tp_max.calc_force<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<std::endl;
} 


void getTimeProfileMax(const PS::TimeProfile & tp, const PS::S32 rank, PS::TimeProfile & tp_max, PS::S32 rank_max[]){
    PS::S32 id = 0;
    PS::Comm::getMaxValue(tp.collect_sample_particle, rank, tp_max.collect_sample_particle, rank_max[id++]);
    PS::Comm::getMaxValue(tp.decompose_domain, rank, tp_max.decompose_domain, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_particle, rank, tp_max.exchange_particle, rank_max[id++]);
    PS::Comm::getMaxValue(tp.set_particle_local_tree, rank, tp_max.set_particle_local_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.set_root_cell, rank, tp_max.set_root_cell, rank_max[id++]);
    PS::Comm::getMaxValue(tp.make_local_tree, rank, tp_max.make_local_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.calc_moment_local_tree, rank, tp_max.calc_moment_local_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.make_LET_1st, rank, tp_max.make_LET_1st, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_LET_1st, rank, tp_max.exchange_LET_1st, rank_max[id++]);
    PS::Comm::getMaxValue(tp.make_LET_2nd, rank, tp_max.make_LET_2nd, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_LET_2nd, rank, tp_max.exchange_LET_2nd, rank_max[id++]);
    PS::Comm::getMaxValue(tp.set_particle_global_tree, rank, tp_max.set_particle_global_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.make_global_tree, rank, tp_max.make_global_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.calc_moment_global_tree, rank, tp_max.calc_moment_global_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.calc_force, rank, tp_max.calc_force, rank_max[id++]);
    //if(PS::Comm::getRank() == 0){ std::cout<<"tp_max.calc_force="<<tp_max.calc_force<<std::endl; }
}



int main(int argc, char* argv[]){
    PS::F64 flops_per_proc = 1.28e11; // for K computer (8thread)
    PS::F64 n_op_dens = 42.0*3.0 + 74.0;
    PS::F64 n_op_hydr = 119.0;
    PS::F64 n_op_ep_ep_grav = 29.0;
    PS::F64 n_op_ep_sp_grav = 29.0;
    PS::F64 wtime_dens = 0.0;
    PS::F64 wtime_hydr = 0.0;
    PS::F64 wtime_grav = 0.0;

	//////////////////
	//Create vars.
	//////////////////
    std::cout<<"initialize"<<std::endl;
	PS::Initialize(argc, argv);
    std::cout<<"end of initialize"<<std::endl;
    PS::F64 wtime0 = PS::GetWtime();
    std::cout<<"start"<<std::endl;
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
    std::cout<<"befor SetupIC: GetWimte()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	SetupIC(sph_system, sysinfo, dinfo);
    std::cout<<"after SetupIC: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	std::cout << "Init: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	Initialize(sph_system);
	//Dom. info
	std::cout << "decomp. All: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	dinfo.decomposeDomainAll(sph_system);
	std::cout << "ex ptcl: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	sph_system.exchangeParticle(dinfo);
	//plant tree
	std::cout << "dens tree: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	PS::TreeForForceShort<RESULT::Dens , EPI::Dens , EPJ::Dens >::Gather   dens_tree;
	std::cout << "hydro tree: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	PS::TreeForForceShort<RESULT::Hydro, EPI::Hydro, EPJ::Hydro>::Symmetry hydr_tree;
	std::cout << "grav tree: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	PS::TreeForForceLong <RESULT::Grav , EPI::Grav , EPJ::Grav >::Monopole grav_tree;
	const int n_grp_dens = 128;
    //const int n_grp_dens = 64;
	//const int n_grp_dens = 32;
	const int n_leaf_dens = 1;
	std::cout << "Init. dens: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	dens_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, n_leaf_dens, n_grp_dens);
	std::cout<<"dens_tree.getMemSize()="<<dens_tree.getMemSizeUsed()*1e-9<<std::endl;

	std::cout << "Init. hydro: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	hydr_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, n_leaf_dens, n_grp_dens);
	std::cout<<"hydr_tree.getMemSize()="<<hydr_tree.getMemSizeUsed()*1e-9<<std::endl;

	std::cout << "Init. grav: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	grav_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, 8, 512);
	std::cout<<"grav_tree.getMemSize()="<<grav_tree.getMemSizeUsed()*1e-9<<std::endl;

	std::cout << "Calc Dens: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
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
        //dens_tree.calcForceAllWithTimer(CalcDensity(), sph_system, dinfo, timer, true);
        dens_tree.calcForceAll(CalcDensity(), sph_system, dinfo, true);
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
        //dens_tree.calcForceAllWithTimer(CalcDensity(), sph_system, dinfo, timer, true);
        dens_tree.calcForceAll(CalcDensity(), sph_system, dinfo, true);
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
        sph_system[i].setBalsalaSwitch();
        sph_system[i].Rsearch = sph_system[i].smth * kernel_t::supportRadius();
    }
#endif
    std::cout<<"sph_system.getMemSize()="<<sph_system.getMemSizeUsed()*1e-9<<std::endl;
    std::cout<<"dens_tree.getMemSize()="<<dens_tree.getMemSizeUsed()*1e-9<<std::endl;

	//std::cout << "Calc Pres" << std::endl;
	//CalcPressure(sph_system);

	std::cout << "Calc Hydro: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
    std::cout<<"hydr_tree.getMemSize()="<<hydr_tree.getMemSizeUsed()*1e-9<<std::endl;

	std::cout << "Calc Grav: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);
    std::cout<<"grav_tree.getMemSize()="<<grav_tree.getMemSizeUsed()*1e-9<<std::endl;
	
	std::cout << "get dt: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	sysinfo.dt = getTimeStepGlobal(sph_system);
	std::cout << "calc EXT: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	CalcExternalForce(sph_system, sysinfo);

	std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt << std::endl;

	PS::S64 step = 0;
	PS::F64 n_op_tot = 1.0;
	for(sysinfo.time = 0 ; sysinfo.time < sysinfo.end_time ; sysinfo.time += sysinfo.dt, ++ step){
        if(step >= 100) break;
	    std::cout<<"step="<<step<<std::endl;
	    PS::Comm::barrier();
	    dinfo.clearTimeProfile();
	    sph_system.clearTimeProfile();
	    grav_tree.clearTimeProfile();
	    dens_tree.clearTimeProfile();
	    hydr_tree.clearTimeProfile();
	    grav_tree.clearTimeProfile();
	    dens_tree.clearNumberOfInteraction();
	    hydr_tree.clearNumberOfInteraction();
	    grav_tree.clearNumberOfInteraction();


	    PS::F64 wtime_offset = PS::GetWtime();
	    timer.reset();
	    timer.start();

	    //std::cout<<"check 1"<<std::endl;

	    InitialKick(sph_system, sysinfo);
	    timer.restart("InitialKick");

	    //std::cout<<"check 2"<<std::endl;

	    FullDrift(sph_system, sysinfo);
	    timer.restart("FullDrift");

	    //std::cout<<"check 3"<<std::endl;

	    Predict(sph_system, sysinfo);
	    timer.restart("Predict");

	    //std::cout<<"check 4"<<std::endl;

	    if(step < 2){
	      dinfo.collectSampleParticle(sph_system);
	      //timer.restart("collectSampleParticle");
	      dinfo.decomposeDomainMultiStep();
	    }
	    else if(step < 12 || step % 4 == 0){
            //dinfo.collectSampleParticle(sph_system, true, n_op_tot);
            dinfo.collectSampleParticle(sph_system, true);
	      //timer.restart("collectSampleParticle");
	      dinfo.decomposeDomainMultiStep();
	    }
	    else{
            //timer.restart("collectSampleParticle");
	    }
	    timer.restart("collect_and_decompose_domain");

	    //std::cout<<"check 5"<<std::endl;

	    sph_system.exchangeParticle(dinfo);
	    timer.restart("exchangeParticle");

	    //std::cout<<"sph_system.getMemSize()="<<sph_system.getMemSizeUsed()*1e-9<<std::endl;

	    //std::cout<<"check 6"<<std::endl;

#pragma omp parallel for
	    for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].Rsearch = kernel_t::supportRadius() * sph_system[i].smth * EXPAND;
	    }

	    //std::cout<<"check 7"<<std::endl;

	    timer.restart("setRsearch");

	    PS::S64 n_int_ep_ep = 0;
	    wtime_dens = PS::GetWtime();
#if 0
	    cnt = 0;
	    for(bool repeat = true; repeat==true;){

		//std::cout<<"cnt="<<cnt<<std::endl;
		cnt++;
		bool repeat_loc = false;
		repeat = false;
		//dens_tree.calcForceAllWithTimer(CalcDensity(), sph_system, dinfo, timer, true);
		dens_tree.calcForceAll(CalcDensity(), sph_system, dinfo, true);
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
#else
        cnt = 0;
        for(bool repeat = true; repeat==true;){
	    cnt++;
            bool repeat_loc = false;
            repeat = false;
            //dens_tree.calcForceAllWithTimer(CalcDensity(), sph_system, dinfo, timer, true);
            dens_tree.calcForceAll(CalcDensity(), sph_system, dinfo, true);
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
        //std::cout<<"check 8"<<std::endl;
        dens_tree.setNInteractionEPEP(n_int_ep_ep);
        //std::cout<<"dens_tree.getMemSize()="<<dens_tree.getMemSizeUsed()*1e-9<<std::endl;
        wtime_dens = PS::GetWtime() - wtime_dens;
        timer.restart("calc_ens");
	

        //std::cout<<"check 9"<<std::endl;
        CalcPressure(sph_system);
        //std::cout<<"check 10"<<std::endl;
        timer.restart("CalcPressure");

        for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
            sph_system[i].setBalsalaSwitch();
            sph_system[i].Rsearch = sph_system[i].smth * kernel_t::supportRadius();
        }
        timer.restart("setRsearch");
        //std::cout<<"check 11"<<std::endl;
        //wtime_hydr = PS::GetWtime();
		//hydr_tree.calcForceAllAndWriteBackWithTimer(CalcHydroForce(), sph_system, dinfo, timer, true);
        hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo, true);
        //std::cout<<"hydr_tree.getMemSize()="<<hydr_tree.getMemSizeUsed()*1e-9<<std::endl;
        //wtime_hydr = PS::GetWtime() - wtime_hydr;
        timer.restart("calc_hydr");

        //wtime_grav = PS::GetWtime();
		//grav_tree.calcForceAllAndWriteBackWithTimer(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo, timer, true);
        grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo, true);
        //std::cout<<"grav_tree.getMemSize()="<<grav_tree.getMemSizeUsed()*1e-9<<std::endl;
        //wtime_grav = PS::GetWtime() - wtime_grav;
        timer.restart("calc_grav");

        //std::cout<<"check 12"<<std::endl;
		sysinfo.dt = getTimeStepGlobal(sph_system);
        timer.restart("set_dt");

        //std::cout<<"check 13"<<std::endl;
		CalcExternalForce(sph_system, sysinfo);
        timer.restart("CalcExternalForce");

        //std::cout<<"check 14"<<std::endl;
		FinalKick(sph_system, sysinfo);
        timer.restart("FinalKick");

        //std::cout<<"check 15"<<std::endl;
        PS::Comm::barrier();
        PS::F64 wtime_tot = PS::GetWtime() - wtime_offset;

        PS::TimeProfile tp_dinfo = dinfo.getTimeProfile();
        PS::TimeProfile tp_system = sph_system.getTimeProfile();
        PS::TimeProfile tp_dens = dens_tree.getTimeProfile();
        PS::TimeProfile tp_hydr = hydr_tree.getTimeProfile();
        PS::TimeProfile tp_grav = grav_tree.getTimeProfile();
        PS::TimeProfile tp_dinfo_max, tp_system_max, tp_dens_max, tp_hydr_max, tp_grav_max;
        PS::S32 rank_dinfo_max[15], rank_system_max[15], rank_dens_max[15], 
            rank_hydr_max[15], rank_grav_max[15];
        getTimeProfileMax(tp_dinfo, PS::Comm::getRank(), tp_dinfo_max, rank_dinfo_max);
        getTimeProfileMax(tp_system, PS::Comm::getRank(), tp_system_max, rank_system_max);
        getTimeProfileMax(tp_dens, PS::Comm::getRank(), tp_dens_max, rank_dens_max);
        getTimeProfileMax(tp_hydr, PS::Comm::getRank(), tp_hydr_max, rank_hydr_max);
        getTimeProfileMax(tp_grav, PS::Comm::getRank(), tp_grav_max, rank_grav_max);
        PS::CountT n_int_ep_ep_dens = dens_tree.getNumberOfInteractionEPEPGlobal();
        PS::CountT n_int_ep_ep_hydr = hydr_tree.getNumberOfInteractionEPEPGlobal();
        PS::CountT n_int_ep_ep_grav = grav_tree.getNumberOfInteractionEPEPGlobal();
        PS::CountT n_int_ep_sp_grav = grav_tree.getNumberOfInteractionEPSPGlobal();
        PS::CountT n_op_tot = n_int_ep_ep_dens*n_op_dens + n_int_ep_ep_hydr*n_op_hydr + n_int_ep_ep_grav * n_op_ep_ep_grav + n_int_ep_sp_grav * n_op_ep_sp_grav;
        fout_tcal<<"speed="<<(PS::F64)(n_op_tot)/(wtime_tot)*1e-12<<"[Tflops]"<<std::endl;
        fout_tcal<<"efficiency="<<(PS::F64)(n_op_tot)/(wtime_tot)/(flops_per_proc*PS::Comm::getNumberOfProc())<<std::endl;
        fout_tcal<<"wtime_tot="<<wtime_tot<<std::endl;
        fout_tcal<<"n_op_tot="<<n_op_tot<<std::endl;

        timer.dump(fout_tcal);

        DumpTimeProfile0(tp_dinfo, tp_dinfo_max, rank_dinfo_max, fout_tcal);
        DumpTimeProfile1(tp_system, tp_system_max, rank_system_max, fout_tcal);
        fout_tcal<<"n_int_ep_ep_dens="<<n_int_ep_ep_dens<<" cnt="<<cnt<<std::endl;
        fout_tcal<<"ni_ave="<<(PS::F64)(sph_system.getNumberOfParticleGlobal())/dens_tree.getNumberOfWalkGlobal()
                 <<" nj_ave="<<(PS::F64)(n_int_ep_ep_dens)/sph_system.getNumberOfParticleGlobal()<<std::endl;
        DumpTimeProfile2(tp_dens, tp_dens_max, rank_dens_max, fout_tcal);
        fout_tcal<<"n_int_ep_ep_hydr="<<n_int_ep_ep_hydr<<std::endl;
        fout_tcal<<"ni_ave="<<(PS::F64)(sph_system.getNumberOfParticleGlobal())/hydr_tree.getNumberOfWalkGlobal()
                 <<" nj_ave="<<(PS::F64)(n_int_ep_ep_hydr)/sph_system.getNumberOfParticleGlobal()<<std::endl;
        DumpTimeProfile2(tp_hydr, tp_hydr_max, rank_hydr_max, fout_tcal);
        fout_tcal<<"n_int_ep_ep_grav="<<n_int_ep_ep_grav<<", n_int_ep_sp_grav="<<n_int_ep_sp_grav<<std::endl;
        fout_tcal<<"ni_ave="<<(PS::F64)(sph_system.getNumberOfParticleGlobal())/grav_tree.getNumberOfWalkGlobal()
                 <<" nj_ave(EP-EP)="<<(PS::F64)(n_int_ep_ep_grav)/sph_system.getNumberOfParticleGlobal()
                 <<" nj_ave(EP-SP)="<<(PS::F64)(n_int_ep_sp_grav)/sph_system.getNumberOfParticleGlobal()<<std::endl;
        DumpTimeProfile2(tp_grav, tp_grav_max, rank_grav_max, fout_tcal);

        fout_tcal<<std::endl;
        fout_tcal<<std::endl;
        fout_tcal<<std::endl;

        //std::cout<<"check 16"<<std::endl;
	ShiftOrigin(sph_system);
        //std::cout<<"check 17"<<std::endl;
	}

	PS::Finalize();
	return 0;
}

