//#define SANITY_CHECK_REALLOCATABLE_ARRAY
#define FAST
#include "header.h"

#ifdef INTRINSIC_K
#include "force-k2.h"
#else 
#include "force.h"
#endif
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
    fout<<"calc_force= "<<tp.calc_force<<", max= "<<tp_max.calc_force<<", rank= "<<rank_max[id++]<<std::endl;
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
    fout<<"calc_force= "<<tp.calc_force<<", max= "<<tp_max.calc_force<<", rank= "<<rank_max[id++]<<std::endl;

    fout<<"make_local_tree_tot= "<<tp.make_local_tree_tot<<", max= "<<tp_max.make_local_tree_tot<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_global_tree_tot= "<<tp.make_global_tree_tot<<", max= "<<tp_max.make_global_tree_tot<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"exchange_LET_tot= "<<tp.exchange_LET_tot<<", max= "<<tp_max.exchange_LET_tot<<", rank= "<<rank_max[id++]<<std::endl;

    fout<<"morton_sort_local_tree= "<<tp.morton_sort_local_tree<<", max= "<<tp_max.morton_sort_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"link_cell_local_tree= "<<tp.link_cell_local_tree<<", max= "<<tp_max.link_cell_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"morton_sort_global_tree= "<<tp.morton_sort_global_tree<<", max= "<<tp_max.morton_sort_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"link_cell_global_tree= "<<tp.link_cell_global_tree<<", max= "<<tp_max.link_cell_global_tree<<", rank= "<<rank_max[id++]<<std::endl;

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

    PS::Comm::getMaxValue(tp.make_local_tree_tot, rank, tp_max.make_local_tree_tot, rank_max[id++]);
    PS::Comm::getMaxValue(tp.make_global_tree_tot, rank, tp_max.make_global_tree_tot, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_LET_tot, rank, tp_max.exchange_LET_tot, rank_max[id++]);

    PS::Comm::getMaxValue(tp.morton_sort_local_tree, rank, tp_max.morton_sort_local_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.link_cell_local_tree, rank, tp_max.link_cell_local_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.morton_sort_global_tree, rank, tp_max.morton_sort_global_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.link_cell_global_tree, rank, tp_max.link_cell_global_tree, rank_max[id++]);
    //if(PS::Comm::getRank() == 0){ std::cout<<"tp_max.calc_force="<<tp_max.calc_force<<std::endl; }
}



int main(int argc, char* argv[]){

    std::cout<<"initialize"<<std::endl;
    PS::Initialize(argc, argv);
    std::cout<<"end of initialize"<<std::endl;
    PS::F64 wtime0 = PS::GetWtime();

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
	PS::TreeForForceShort<RESULT::Dens , EPI::Dens , EPJ::Dens>::Gather   dens_tree;
    dens_tree.setPrefixOfProfile("dens");

	std::cout << "hydro tree: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	PS::TreeForForceShort<RESULT::Hydro, EPI::Hydro, EPJ::Hydro>::Symmetry hydr_tree;
    hydr_tree.setPrefixOfProfile("hydr");

	std::cout << "grav tree: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	PS::TreeForForceLong <RESULT::Grav, EPI::Grav , EPJ::Grav>::Monopole grav_tree;
    grav_tree.setPrefixOfProfile("grav");


    const int n_grp_dens = 128;
	const int n_leaf_dens = 8;
	const int n_grp_hydr = 128;
	const int n_leaf_hydr = 8;
	const int n_grp_grav = 512;
	const int n_leaf_grav = 8;

/*
	const int n_grp_dens = 32;
	const int n_leaf_dens = 8;
	const int n_grp_hydr = 32;
	const int n_leaf_hydr = 8;
	const int n_grp_grav = 32;
	const int n_leaf_grav = 8;
*/
	//const int n_leaf_dens = 1;
	std::cout << "Init. dens: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	dens_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, n_leaf_dens, n_grp_dens);
	std::cout<<"dens_tree.getMemSize()="<<dens_tree.getMemSizeUsed()*1e-9<<std::endl;

	std::cout << "Init. hydro: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	hydr_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, n_leaf_hydr, n_grp_hydr);
	std::cout<<"hydr_tree.getMemSize()="<<hydr_tree.getMemSizeUsed()*1e-9<<std::endl;

	std::cout << "Init. grav: PS::GetWtime()-wtime0<<std::endl="<<PS::GetWtime()-wtime0<<std::endl;
	grav_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, n_leaf_grav, n_grp_grav);
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
	PS::F64 weight = 0.0;
	for(sysinfo.time = 0 ; sysinfo.time < sysinfo.end_time ; sysinfo.time += sysinfo.dt, step++){
        //if(step >= 1) break;
	    PS::Comm::barrier();
	    dinfo.clearTimeProfile();
	    sph_system.clearTimeProfile();
/*
	    dens_tree.clearTimeProfile();
	    hydr_tree.clearTimeProfile();
	    grav_tree.clearTimeProfile();
	    dens_tree.clearNumberOfInteraction();
	    hydr_tree.clearNumberOfInteraction();
	    grav_tree.clearNumberOfInteraction();
*/
	    dens_tree.clearCounterAll();
	    hydr_tree.clearCounterAll();
	    grav_tree.clearCounterAll();

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

	    if(step < 4){
	      dinfo.collectSampleParticle(sph_system);
	      //timer.restart("collectSampleParticle");
	      dinfo.decomposeDomainMultiStep();
	    }
	    else if(step < 12 || step % 4 == 0){
            //dinfo.collectSampleParticle(sph_system, true, n_op_tot);
            //dinfo.collectSampleParticle(sph_system, true);
            dinfo.collectSampleParticle(sph_system, true, weight);
            //timer.restart("collectSampleParticle");
            dinfo.decomposeDomainMultiStep();
            weight = 0.0;
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
        PS::S32 n_remain[100]; 
        for(bool repeat = true; repeat==true;){
            bool repeat_loc = false;
            repeat = false;
            n_remain[cnt] = 0;
            //dens_tree.calcForceAllWithTimer(CalcDensity(), sph_system, dinfo, timer, true);
            dens_tree.calcForceAll(CalcDensity(), sph_system, dinfo, true);
            for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
                if(sph_system[i].Rsearch != 0.0){
                    if(dens_tree.getForce(i).itr == true){
                        repeat_loc = true;
                        sph_system[i].Rsearch *= EXPAND;
                        n_remain[cnt]++;
                    }
                    else{
                        sph_system[i].Rsearch = 0.0;
                        sph_system[i].copyFromForce(dens_tree.getForce(i));
                    }
                }
            }
            repeat = PS::Comm::synchronizeConditionalBranchOR(repeat_loc);
            n_int_ep_ep += dens_tree.getNInteractionEPEP();
            cnt++;
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
        hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo, true);
        //hydr_tree.calcForceAllWalkOnlyAndWriteBack(CalcHydroForce(), sph_system, dinfo, true);

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

        weight += dens_tree.getTimeProfile().calc_force + hydr_tree.getTimeProfile().calc_force;
        //std::cout<<"check 15"<<std::endl;
        PS::Comm::barrier();
        PS::F64 wtime_tot = PS::GetWtime() - wtime_offset;

        PS::TimeProfile tp_dinfo = dinfo.getTimeProfile();
        PS::TimeProfile tp_system = sph_system.getTimeProfile();
        PS::TimeProfile tp_dens = dens_tree.getTimeProfile();
        PS::TimeProfile tp_hydr = hydr_tree.getTimeProfile();
        PS::TimeProfile tp_grav = grav_tree.getTimeProfile();
        PS::TimeProfile tp_dinfo_max, tp_system_max, tp_dens_max, tp_hydr_max, tp_grav_max;
        const PS::S32 n_profile = 30;
        PS::S32 rank_dinfo_max[n_profile], rank_system_max[n_profile], rank_dens_max[n_profile], 
            rank_hydr_max[n_profile], rank_grav_max[n_profile];
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
        fout_tcal<<"sysinfo.time= "<<sysinfo.time<<" step= "<<step<<std::endl;
        fout_tcal<<"speed="<<(PS::F64)(n_op_tot)/(wtime_tot)*1e-12<<"[Tflops]"<<std::endl;
        fout_tcal<<"efficiency="<<(PS::F64)(n_op_tot)/(wtime_tot)/(flops_per_proc*PS::Comm::getNumberOfProc())<<std::endl;
        fout_tcal<<"dinfo.getPosDomain(PS::Comm::getRank()).getCenter()= "<<dinfo.getPosDomain(PS::Comm::getRank()).getCenter()
                 <<" dinfo.getPosDomain(PS::Comm::getRank()).getHalfLength()= "<<dinfo.getPosDomain(PS::Comm::getRank()).getHalfLength()<<std::endl;
        fout_tcal<<"dinfo.getPosDomain(PS::Comm::getRank()).low_= "<<dinfo.getPosDomain(PS::Comm::getRank()).low_
                 <<" dinfo.getPosDomain(PS::Comm::getRank()).high_= "<<dinfo.getPosDomain(PS::Comm::getRank()).high_<<std::endl;
        fout_tcal<<"wtime_tot= "<<wtime_tot<<std::endl;
        fout_tcal<<"weight= "<<weight<<std::endl;
        fout_tcal<<"sph_system.getNumberOfParticleLocal()= "<<sph_system.getNumberOfParticleLocal()
                 <<" sph_system.getNumberOfParticleGlobal()= "<<sph_system.getNumberOfParticleGlobal()<<std::endl;
        fout_tcal<<"n_op_tot= "<<n_op_tot<<std::endl;

        fout_tcal<<"n_ptcl_send_local= "<<sph_system.getNumberOfParticleSendLocal()
                 <<" n_ptcl_recv_local= "<<sph_system.getNumberOfParticleRecvLocal()
                 <<" n_ptcl_send_golbal= "<<sph_system.getNumberOfParticleSendGlobal()
                 <<" n_ptcl_recv_golbal= "<<sph_system.getNumberOfParticleRecvGlobal()<<std::endl;

        fout_tcal<<"n_let_ep_send_1st_local(grav)= "<<grav_tree.getNumberOfLETEPSend1stLocal()
                 <<" n_let_ep_recv_1st_local(grav)= "<<grav_tree.getNumberOfLETEPRecv1stLocal()
                 <<" n_let_sp_send_1st_local(grav)= "<<grav_tree.getNumberOfLETEPSend1stLocal()
                 <<" n_let_sp_recv_1st_local(grav)= "<<grav_tree.getNumberOfLETEPRecv1stLocal()<<std::endl;
        fout_tcal<<"n_let_ep_send_1st_global(grav)= "<<grav_tree.getNumberOfLETEPSend1stGlobal()
                 <<" n_let_ep_recv_1st_global(grav)= "<<grav_tree.getNumberOfLETEPRecv1stGlobal()
                 <<" n_let_sp_send_1st_global(grav)= "<<grav_tree.getNumberOfLETSPSend1stGlobal()
                 <<" n_let_sp_recv_1st_global(grav)= "<<grav_tree.getNumberOfLETSPRecv1stGlobal()<<std::endl;
        fout_tcal<<"n_cell_openl_local(grav)= "<<grav_tree.getNumberOfCellOpenLocal()
                 <<" n_cell_openl_global(grav)= "<<grav_tree.getNumberOfCellOpenGlobal()<<std::endl;

        fout_tcal<<"n_let_ep_send_1st_local(dens)= "<<dens_tree.getNumberOfLETEPSend1stLocal()
                 <<" n_let_ep_recv_1st_local(dens)= "<<dens_tree.getNumberOfLETEPRecv1stLocal()
                 <<" n_let_ep_send_2nd_local(dens)= "<<dens_tree.getNumberOfLETEPSend2ndLocal()
                 <<" n_let_ep_recv_2nd_local(dens)= "<<dens_tree.getNumberOfLETEPRecv2ndLocal()<<std::endl;
        fout_tcal<<"n_let_ep_send_1st_global(dens)= "<<dens_tree.getNumberOfLETEPSend1stGlobal()
                 <<" n_let_ep_recv_1st_global(dens)= "<<dens_tree.getNumberOfLETEPRecv1stGlobal()
                 <<" n_let_sp_send_2nd_global(dens)= "<<dens_tree.getNumberOfLETEPSend2ndGlobal()
                 <<" n_let_sp_recv_2nd_global(dens)= "<<dens_tree.getNumberOfLETEPRecv2ndGlobal()<<std::endl;
        fout_tcal<<"n_cell_open_local(dens)= "<<dens_tree.getNumberOfCellOpenLocal()
                 <<" n_cell_open_global(dens)= "<<dens_tree.getNumberOfCellOpenGlobal()<<std::endl;

        fout_tcal<<"n_let_ep_send_1st_local(hydr)= "<<hydr_tree.getNumberOfLETEPSend1stLocal()
                 <<" n_let_ep_recv_1st_local(hydr)= "<<hydr_tree.getNumberOfLETEPRecv1stLocal()
                 <<" n_let_ep_send_2nd_local(hydr)= "<<hydr_tree.getNumberOfLETEPSend2ndLocal()
                 <<" n_let_ep_recv_2nd_local(hydr)= "<<hydr_tree.getNumberOfLETEPRecv2ndLocal()<<std::endl;
        fout_tcal<<"n_let_ep_send_1st_global(hydr)= "<<hydr_tree.getNumberOfLETEPSend1stGlobal()
                 <<" n_let_ep_recv_1st_global(hydr)= "<<hydr_tree.getNumberOfLETEPRecv1stGlobal()
                 <<" n_let_sp_send_2nd_global(hydr)= "<<hydr_tree.getNumberOfLETEPSend2ndGlobal()
                 <<" n_let_sp_recv_2nd_global(hydr)= "<<hydr_tree.getNumberOfLETEPRecv2ndGlobal()<<std::endl;
        fout_tcal<<"n_cell_open_local(hydr)= "<<hydr_tree.getNumberOfCellOpenLocal()
                 <<" n_cell_open_global(hydr)= "<<hydr_tree.getNumberOfCellOpenGlobal()<<std::endl;


        timer.dump(fout_tcal);

        DumpTimeProfile0(tp_dinfo, tp_dinfo_max, rank_dinfo_max, fout_tcal);
        DumpTimeProfile1(tp_system, tp_system_max, rank_system_max, fout_tcal);
        fout_tcal<<"n_int_ep_ep_dens="<<n_int_ep_ep_dens<<" cnt="<<cnt<<std::endl;
        fout_tcal<<"n_remain[0]"<<n_remain[0];
        for(PS::S32 i=1; i<cnt; i++) fout_tcal<<" n_remain[i]="<<n_remain[i];
        fout_tcal<<std::endl;
        fout_tcal<<"ni_ave(local)= "<<(PS::F64)(sph_system.getNumberOfParticleLocal())/dens_tree.getNumberOfWalkLocal()
                 <<" nj_ave(local)= "<<(PS::F64)(dens_tree.getNumberOfInteractionEPEPLocal())/sph_system.getNumberOfParticleLocal()
                 <<" ni_ave(global)= "<<(PS::F64)(sph_system.getNumberOfParticleGlobal())/dens_tree.getNumberOfWalkGlobal()
                 <<" nj_ave(global)= "<<(PS::F64)(n_int_ep_ep_dens)/sph_system.getNumberOfParticleGlobal()
                 <<" dens_tree.getNumberOfWalkLocal()= "<<dens_tree.getNumberOfWalkLocal()
                 <<" dens_tree.getNumberOfWalkGlobal()= "<<dens_tree.getNumberOfWalkGlobal()<<std::endl;
        DumpTimeProfile2(tp_dens, tp_dens_max, rank_dens_max, fout_tcal);
        fout_tcal<<"n_int_ep_ep_hydr="<<n_int_ep_ep_hydr<<std::endl;
        fout_tcal<<"ni_ave(local)= "<<(PS::F64)(sph_system.getNumberOfParticleLocal())/hydr_tree.getNumberOfWalkLocal()
                 <<" nj_ave(local)= "<<(PS::F64)(hydr_tree.getNumberOfInteractionEPEPLocal())/sph_system.getNumberOfParticleLocal()
                 <<" ni_ave(Global)= "<<(PS::F64)(sph_system.getNumberOfParticleGlobal())/hydr_tree.getNumberOfWalkGlobal()
                 <<" nj_ave(Global)= "<<(PS::F64)(n_int_ep_ep_hydr)/sph_system.getNumberOfParticleGlobal()
                 <<" hydr_tree.getNumberOfWalkLocal()= "<<hydr_tree.getNumberOfWalkLocal()
                 <<" hydr_tree.getNumberOfWalkGlobal()= "<<hydr_tree.getNumberOfWalkGlobal()<<std::endl;
        DumpTimeProfile2(tp_hydr, tp_hydr_max, rank_hydr_max, fout_tcal);
        fout_tcal<<"n_int_ep_ep_grav="<<n_int_ep_ep_grav<<", n_int_ep_sp_grav="<<n_int_ep_sp_grav<<std::endl;
        fout_tcal<<"ni_ave(local)= "<<(PS::F64)(sph_system.getNumberOfParticleLocal())/grav_tree.getNumberOfWalkLocal()
                 <<" nj_ave(EP-EP)(local)= "<<(PS::F64)(grav_tree.getNumberOfInteractionEPEPLocal())/sph_system.getNumberOfParticleLocal()
                 <<" nj_ave(EP-SP)(local)= "<<(PS::F64)(grav_tree.getNumberOfInteractionEPSPLocal())/sph_system.getNumberOfParticleLocal()
                 <<" ni_ave(global)= "<<(PS::F64)(sph_system.getNumberOfParticleGlobal())/grav_tree.getNumberOfWalkGlobal()
                 <<" nj_ave(EP-EP)(global)= "<<(PS::F64)(n_int_ep_ep_grav)/sph_system.getNumberOfParticleGlobal()
                 <<" nj_ave(EP-SP)(global)= "<<(PS::F64)(n_int_ep_sp_grav)/sph_system.getNumberOfParticleGlobal()
                 <<" grav_tree.getNumberOfWalkLocal()= "<<grav_tree.getNumberOfWalkLocal()
                 <<" grav_tree.getNumberOfWalkGlobal()= "<<grav_tree.getNumberOfWalkGlobal()<<std::endl;
        DumpTimeProfile2(tp_grav, tp_grav_max, rank_grav_max, fout_tcal);

        fout_tcal<<std::endl;
        fout_tcal<<"dens_tree.getMemSize()= "<<dens_tree.getMemSizeUsed()*1e-9<<"[GB]"
                 <<" hydr_tree.getMemSize()= "<<hydr_tree.getMemSizeUsed()*1e-9<<"[GB]"
                 <<" grav_tree.getMemSize()= "<<grav_tree.getMemSizeUsed()*1e-9<<"[GB]"
                 <<" sph_system.getMemSize()= "<<sph_system.getMemSizeUsed()*1e-9<<"[GB]"<<std::endl;
        fout_tcal<<std::endl;
        fout_tcal<<std::endl;

        //std::cout<<"check 16"<<std::endl;
        ShiftOrigin(sph_system);
        //std::cout<<"check 17"<<std::endl;
	}

	PS::Finalize();
	return 0;
}

