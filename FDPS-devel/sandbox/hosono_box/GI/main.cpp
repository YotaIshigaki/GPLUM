//#define SANITY_CHECK_REALLOCATABLE_ARRAY
#include "header.h"

#ifdef WORKING_ON_K
	#warning Working on K.
	#ifdef FAST
		#warning FUSION
		#include "force-k2.h"
	#else
		#warning NOT FUSION
		#include "force-k.h"
	#endif
#else
	#warning NOT K
	#include "force.h"
#endif

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

void CalcDensityAndSmoothing(PS::ParticleSystem<RealPtcl>& sph_system, PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather& dens_tree, PS::DomainInfo& dinfo){
	#ifdef FAST
	PS::F64 EXPAND = 1.3;
	for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
	    sph_system[i].Rsearch = kernel_t::supportRadius() * sph_system[i].smth * EXPAND;
	}
	int cnt = 0;
	for(bool repeat = true; repeat==true;){
		//
		dens_tree.clearTimeProfile();
		dens_tree.clearNumberOfInteraction();
		int miss = 0;
		double st = MPI_Wtime();
		//
		std::cout << "cnt= " << cnt <<std::endl;
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
					++ miss;
				}else{
					sph_system[i].Rsearch = 0.0;
					sph_system[i].copyFromForce(dens_tree.getForce(i));
				}
			}
		}
		repeat = PS::Comm::synchronizeConditionalBranchOR(repeat_loc);
		//
		PS::TimeProfile tp_dens = dens_tree.getTimeProfile();
		PS::CountT n_int_ep_ep_dens = dens_tree.getNumberOfInteractionEPEPGlobal();
		PS::TimeProfile tp_dens_max;
		PS::S32 rank_dens_max[15];
		//getTimeProfileMax(tp_dens, PS::Comm::getRank(), tp_dens_max, rank_dens_max);
		//DumpTimeProfile2(tp_dens, tp_dens_max, rank_dens_max, std::cout);

		std::cout << "WTIME... " << MPI_Wtime() - st << std::endl;
		miss = PS::Comm::getSum(miss);
		const double missrate = static_cast<double>(miss) / static_cast<double>(sph_system.getNumberOfParticleGlobal());
		std::cout << "MISS ... " << missrate << "(" << miss << ")" << std::endl;
		if(missrate < 0.1e-2) EXPAND *= 2.0;
	}
	CalcPressure(sph_system);
	for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].setBalsalaSwitch();
		sph_system[i].Rsearch = sph_system[i].smth * kernel_t::supportRadius();
	}

	#else
	for(int loop = 0 ; loop <= 2 ; ++ loop){
		dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
	}
	#endif
}

void removeParticle(PS::ParticleSystem<RealPtcl>& sph_system, const std::size_t i){
	const std::size_t Nptcl = sph_system.getNumberOfParticleLocal();
	sph_system[i] = sph_system[Nptcl - 1];
	sph_system.setNumberOfParticleLocal(Nptcl - 1);
}


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
	//Utilize
	sph_system.setAverageTargetNumberOfSampleParticlePerProcess(200);
	const int Ngrp = 128;
	//////////////////
	//Setup Initial
	//////////////////
	if(argc == 1){
		SetupIC(sph_system, sysinfo);
	}else{
		sysinfo.step = atoi(argv[1]);
		//InputFile(sph_system, sysinfo);
		InputFileWithTimeInterval(sph_system, sysinfo);
	}
	std::cout << "Init." << std::endl;
	Initialize(sph_system);
	ShiftOrigin(sph_system);
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
	dens_tree.initialize(sph_system.getNumberOfParticleLocal(), 0.5, 4, Ngrp);
	std::cout << "Init. drvt" << std::endl;
	drvt_tree.initialize(sph_system.getNumberOfParticleLocal(), 0.5, 4, Ngrp);
	std::cout << "Init. hydro" << std::endl;
	hydr_tree.initialize(sph_system.getNumberOfParticleLocal(), 0.5, 4, Ngrp);
	std::cout << "Init. grav" << std::endl;
	grav_tree.initialize(sph_system.getNumberOfParticleLocal(), 0.5, 8, Ngrp * 2);

	std::cout << "Calc Dens" << std::endl;
	CalcDensityAndSmoothing(sph_system, dens_tree, dinfo);
	std::cout << "Calc Pres" << std::endl;
	CalcPressure(sph_system);
	#ifndef FAST
	std::cout << "Calc Drvt" << std::endl;
	drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo);
	#endif
	std::cout << "Calc Hydro" << std::endl;
	hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
	std::cout << "Calc Grav" << std::endl;
	grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);
	
	std::cout << "get dt" << std::endl;
	sysinfo.dt = getTimeStepGlobal(sph_system);
	std::cout << "calc EXT" << std::endl;
	CalcExternalForce(sph_system, sysinfo);

	std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt << std::endl;

	for(; sysinfo.time < sysinfo.end_time ; ++ sysinfo.step){
		//DEBUG
		/*
		if(sysinfo.time > 5500.0){
			std::cout << "DEBUG" << std::endl;
			break;
		}
		*/
		InitialKick(sph_system, sysinfo);
		FullDrift(sph_system, sysinfo);
		//sph_system.adjustPositionIntoRootDomain(dinfo);
		Predict(sph_system, sysinfo);
		dinfo.decomposeDomainAll(sph_system);
		sph_system.exchangeParticle(dinfo);
		CalcDensityAndSmoothing(sph_system, dens_tree, dinfo);
		CalcPressure(sph_system);
		#ifndef FAST
		drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo);
		#endif
		hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
		grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);
		CalcExternalForce(sph_system, sysinfo);
		FinalKick(sph_system, sysinfo);
		OutputFileWithTimeInterval(sph_system, sysinfo);

		sysinfo.time += sysinfo.dt;
		sysinfo.dt = getTimeStepGlobal(sph_system);
		
		if(PS::Comm::getRank() == 0){
			std::cout << "//================================" << std::endl;
			std::cout << std::scientific << std::setprecision(16) << "time = " << sysinfo.time << ", dt = " << sysinfo.dt << std::endl;
			std::cout << "step = " << sysinfo.step << std::endl;
			std::cout << "//================================" << std::endl;
		}
		CheckConservativeVariables(sph_system, sysinfo);
		ShiftOrigin(sph_system);
	}

	OutputFileWithTimeInterval(sph_system, sysinfo);
	PS::Finalize();
	return 0;
}

