//#define SANITY_CHECK_REALLOCATABLE_ARRAY
#include<sys/stat.h>
#include "header.h"

void makeOutputDirectory(const char * dir_name) {
    struct stat st;
    if (stat(dir_name, &st) != 0) {
        PS::S32 ret = -1;
        if (PS::Comm::getRank() == 0)
            ret = mkdir(dir_name, 0777);
        PS::Comm::broadcast(&ret, 1);
        if (ret == 0) {
            if(PS::Comm::getRank() == 0)
                fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
        } else {
            fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
            PS::Abort();
        }
    }
}

int main(int argc, char* argv[]){
	//////////////////
	//Create vars.
	//////////////////
	PS::Initialize(argc, argv);

#if defined(MULTI_WALK) || defined(MULTI_WALK_INDEX)
        const int tag_max = 1;
        const int n_walk_limit = 200;
#endif

	makeOutputDirectory("result");
	PS::ParticleSystem<RealPtcl> sph_system;
	sph_system.initialize();
	PS::DomainInfo dinfo;
	dinfo.initialize();

	PS::F64 dt, end_time;
	//////////////////
	//Disp. Info
	//////////////////
	DisplayInfo();
	//////////////////
	//Setup Initial
	//////////////////
	SetupIC(sph_system, &end_time, dinfo);
	Initialize(sph_system);
	//Dom. info
	dinfo.setDomain(PS::Comm::getNumberOfProc(), 1, 1);
	dinfo.decomposeDomainAll(sph_system);
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

	for(int loop = 0 ; loop <= 5 ; ++ loop){
		dens_tree.calcForceAllAndWriteBack(CalcDensity()   , sph_system, dinfo);
	}
	for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
		sph_system[i].setPressure();
	}
	drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo);
	hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
	grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);

	dt = getTimeStepGlobal(sph_system);

	PS::F64 time = 0.0;
	std::cout << std::scientific << std::setprecision(16) << "time = " << time << ", dt = " << dt << std::endl;

	PS::S32 step = 0;
        bool last_loop = false;
        PS::F64 wtime_loop = 0.0;
        PS::F64 wtime_loop_beg = 0.0;
	//for( ; time < end_time ; time += dt, ++ step){
        for( ; ; ++ step){
            if(time + dt >= end_time){
                last_loop = true;
                dt = end_time - time;
            }
            wtime_loop_beg = PS::GetWtime();
            dens_tree.clearTimeProfile();
            drvt_tree.clearTimeProfile();
            hydr_tree.clearTimeProfile();
            grav_tree.clearTimeProfile();

            time += dt;
            InitialKick(sph_system, dt);
            FullDrift(sph_system, dt);
            sph_system.adjustPositionIntoRootDomain(dinfo);
            Predict(sph_system, dt);
#ifdef REUSE_LIST_MODE
    #ifdef MULTI_WALK_INDEX
                if(step % 4 == 0){
                    //std::cerr<<"make step"<<std::endl;
                    dinfo.decomposeDomainAll(sph_system);
                    sph_system.exchangeParticle(dinfo);
                    for(int loop = 0 ; loop <= 2 ; ++ loop){
			//dens_tree.calcForceAllAndWriteBack(CalcDensity()   , sph_system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
                        dens_tree.calcForceAllAndWriteBackMultiWalkIndex
                            (CalcDensDispatch,
                             RetrieveKernel<RESULT::Dens>,
                             tag_max, sph_system, dinfo, n_walk_limit, true, 
                             PS::MAKE_LIST_FOR_REUSE);
                    }
                    for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].setPressure();
                    }
                    //drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
                    drvt_tree.calcForceAllAndWriteBackMultiWalkIndex
                        (CalcDrvtDispatch,
                         RetrieveKernel<RESULT::Drvt>,
                         tag_max, sph_system, dinfo, n_walk_limit, true, 
                         PS::MAKE_LIST_FOR_REUSE);

                    //hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
                    hydr_tree.calcForceAllAndWriteBackMultiWalkIndex
                        (CalcHydroDispatch,
                         RetrieveKernel<RESULT::Hydro>,
                         tag_max, sph_system, dinfo, n_walk_limit, true, 
                         PS::MAKE_LIST_FOR_REUSE);

                    //grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
                    grav_tree.calcForceAllAndWriteBackMultiWalkIndex
                        (CalcGravDispatch,
                         RetrieveKernel<RESULT::Grav>,
                         tag_max, sph_system, dinfo, n_walk_limit, true, 
                         PS::MAKE_LIST_FOR_REUSE);
                }
                else{
                    //std::cerr<<"reuse step"<<std::endl;
                    //dens_tree.calcForceAllAndWriteBack(CalcDensity()   , sph_system, dinfo, true, PS::REUSE_LIST);
                    dens_tree.calcForceAllAndWriteBackMultiWalkIndex
                        (CalcDensDispatch,
                         RetrieveKernel<RESULT::Dens>,
                         tag_max, sph_system, dinfo, n_walk_limit, true, 
                         PS::REUSE_LIST);
                    for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].setPressure();
                    }
                    //drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo, true, PS::REUSE_LIST);
                    drvt_tree.calcForceAllAndWriteBackMultiWalkIndex
                        (CalcDrvtDispatch,
                         RetrieveKernel<RESULT::Drvt>,
                         tag_max, sph_system, dinfo, n_walk_limit, true, 
                         PS::REUSE_LIST);
                    //hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo, true, PS::REUSE_LIST);
                    hydr_tree.calcForceAllAndWriteBackMultiWalkIndex
                        (CalcHydroDispatch,
                         RetrieveKernel<RESULT::Hydro>,
                         tag_max, sph_system, dinfo, n_walk_limit, true, 
                         PS::REUSE_LIST);


                    //grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo, true, PS::REUSE_LIST);
                    grav_tree.calcForceAllAndWriteBackMultiWalkIndex
                        (CalcGravDispatch,
                         RetrieveKernel<RESULT::Grav>,
                         tag_max, sph_system, dinfo, n_walk_limit, true, 
                         PS::REUSE_LIST);
                }
    #else //MULTI_WALK_INDEX
                if(step % 4 == 0){
                    dinfo.decomposeDomainAll(sph_system);
                    sph_system.exchangeParticle(dinfo);
                    for(int loop = 0 ; loop <= 2 ; ++ loop){
			dens_tree.calcForceAllAndWriteBack(CalcDensity()   , sph_system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
                    }
                    for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].setPressure();
                    }
                    drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
                    hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
                    grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
                }
                else{
                    dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo, true, PS::REUSE_LIST);
                    for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].setPressure();
                    }
                    drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo, true, PS::REUSE_LIST);
                    hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo, true, PS::REUSE_LIST);
                    grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo, true, PS::REUSE_LIST);
                }
    #endif //MULTI_WALK_INDEX
#else //REUSE_LIST_MODE
                // NO REUSE
    #ifdef MULTI_WALK_INDEX
                if(step % 4 == 0){
                    dinfo.decomposeDomainAll(sph_system);
                    sph_system.exchangeParticle(dinfo);
                }
                for(int loop = 0 ; loop <= 2 ; ++ loop){
                    dens_tree.calcForceAllAndWriteBackMultiWalkIndex
                        (CalcDensDispatch,
                         RetrieveKernel<RESULT::Dens>,
                         tag_max, sph_system, dinfo, n_walk_limit, true, 
                         PS::MAKE_LIST);
                }
                for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
                    sph_system[i].setPressure();
                }
                drvt_tree.calcForceAllAndWriteBackMultiWalkIndex
                    (CalcDrvtDispatch,
                     RetrieveKernel<RESULT::Drvt>,
                     tag_max, sph_system, dinfo, n_walk_limit, true, 
                     PS::MAKE_LIST);

                hydr_tree.calcForceAllAndWriteBackMultiWalkIndex
                    (CalcHydroDispatch,
                     RetrieveKernel<RESULT::Hydro>,
                     tag_max, sph_system, dinfo, n_walk_limit, true, 
                     PS::MAKE_LIST);

                grav_tree.calcForceAllAndWriteBackMultiWalkIndex
                    (CalcGravDispatch,
                     RetrieveKernel<RESULT::Grav>,
                     tag_max, sph_system, dinfo, n_walk_limit, true, 
                     PS::MAKE_LIST);
    #else //MULTI_WALK_INDEX
		dinfo.decomposeDomainAll(sph_system);
		sph_system.exchangeParticle(dinfo);
		for(int loop = 0 ; loop <= 2 ; ++ loop){
			dens_tree.calcForceAllAndWriteBack(CalcDensity()   , sph_system, dinfo);
		}
		for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
			sph_system[i].setPressure();
		}
		drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo);
		hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
		grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);
    #endif //MULTI_WALK_INDEX
#endif //REUSE_LIST_MODE
		dt = getTimeStepGlobal(sph_system);

		FinalKick(sph_system, dt);

                PS::TimeProfile dens_tp = dens_tree.getTimeProfile();
                PS::TimeProfile drvt_tp = drvt_tree.getTimeProfile();
                PS::TimeProfile hydr_tp = hydr_tree.getTimeProfile();
                PS::TimeProfile grav_tp = grav_tree.getTimeProfile();
                PS::TimeProfile tree_tp = dens_tp + drvt_tp + hydr_tp + grav_tp;
                wtime_loop = PS::GetWtime() - wtime_loop_beg;

		//if(step % PARAM::OUTPUT_INTERVAL == 0){
                if(step % PARAM::OUTPUT_INTERVAL == 0 || last_loop){
			FileHeader header;
			header.time = time;
			header.Nbody = sph_system.getNumberOfParticleGlobal();
			char filename[256];
			sprintf(filename, "result/%04d.dat", step);
			sph_system.writeParticleAscii(filename, header);
                        /*
			if(PS::Comm::getRank() == 0){
				std::cout << "//================================" << std::endl;
				std::cout << "output " << filename << "." << std::endl;
				std::cout << "//================================" << std::endl;
			}
                        */
		}
                /*
		if(PS::Comm::getRank() == 0){
			std::cout << "//================================" << std::endl;
			std::cout << std::scientific << std::setprecision(16) << "time = " << time << ", dt = " << dt << std::endl;
			std::cout << "step = " << step << std::endl;
			std::cout << "//================================" << std::endl;
		}
                */
		if(PS::Comm::getRank() == 0){
			std::cerr << "//================================" << std::endl;
			std::cerr << std::scientific << std::setprecision(16) << "time = " << time << ", dt = " << dt << std::endl;
			std::cerr << "step = " << step << std::endl;
			std::cerr << "//================================" << std::endl;
		}
                for(PS::S32 i=0; i<sph_system.getNumberOfParticleLocal(); i++){
                    EPI::Grav p_tmp;
                    PS::F64 inv_eps = 1.0 / sqrt(p_tmp.getEps2());
                    sph_system[i].pot += sph_system[i].mass * inv_eps;
                }
		if(PS::Comm::getRank() == 0){
                    std::cout<<"time= "<<time<<" wtime_loop= "<<wtime_loop
                             <<" tree_tp.calc_force= "<<tree_tp.calc_force;
                    std::cerr<<"time= "<<time<<" wtime_loop= "<<wtime_loop
                             <<" tree_tp.calc_force= "<<tree_tp.calc_force
                             <<" dens_tp.calc_force= "<<dens_tp.calc_force
                             <<" drvt_tp.calc_force= "<<drvt_tp.calc_force
                             <<" hydr_tp.calc_force= "<<hydr_tp.calc_force
                             <<" grav_tp.calc_force= "<<grav_tp.calc_force
                             <<std::endl;
                }
		CheckConservativeVariables(sph_system);
                if(last_loop) break;
	}

	PS::Finalize();
	return 0;
}

