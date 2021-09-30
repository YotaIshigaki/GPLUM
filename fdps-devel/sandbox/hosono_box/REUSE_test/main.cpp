#include <iostream>
#include <iomanip>
#include <vector>
#include <assert.h>
#include <cmath>
#include <sys/time.h>
#include <particle_simulator.hpp>

#include "EoS.h"
#include "class_platform.hpp"
#include "init/GI.h"

#ifdef ENABLE_PHANTOM_GRAPE_X86
#include <gp5util.h>
#endif

#ifdef ENABLE_PEZY
#include "pezycl.h"
#endif
#ifdef ENABLE_GPU
#include "kernel.h"
#endif
PS::S32 DensDispatchKernel(const PS::S32, const PS::S32, const EPI::Dens**, const PS::S32*, const EPJ::Dens**, const PS::S32*);
PS::S32 DrvtDispatchKernel(const PS::S32, const PS::S32, const EPI::Drvt**, const PS::S32*, const EPJ::Drvt**, const PS::S32*);
PS::S32 HydrDispatchKernel(const PS::S32, const PS::S32, const EPI::Hydr**, const PS::S32*, const EPJ::Hydr**, const PS::S32*);
PS::S32 GravDispatchKernel(const PS::S32, const PS::S32, const EPI::Grav**, const PS::S32*, const EPJ::Grav**, const PS::S32*, const PS::SPJMonopole**, const PS::S32*);

PS::S32 DensRetrieveKernel(const PS::S32, const PS::S32, const PS::S32*, RESULT::Dens**);
PS::S32 DrvtRetrieveKernel(const PS::S32, const PS::S32, const PS::S32*, RESULT::Drvt**);
PS::S32 HydrRetrieveKernel(const PS::S32, const PS::S32, const PS::S32*, RESULT::Hydr**);
PS::S32 GravRetrieveKernel(const PS::S32, const PS::S32, const PS::S32*, RESULT::Grav**);

double get_dtime(void){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return ((double)(tv.tv_sec) + (double)(tv.tv_usec) * 0.001 * 0.001);
}

int main(int argc, char *argv[]){
	std::cout << std::scientific << std::setprecision(16);
	std::cerr << std::scientific << std::setprecision(16);
	double time[4];
	PS::Initialize(argc, argv);
	#ifdef ENABLE_KNL
	omp_set_num_threads(64);
	#endif
	#ifdef ENABLE_PEZY
	omp_set_num_threads(4);
	extern PezyDevice device;
	device.initialize();
	#endif
	PS::DomainInfo dinfo;
	//
	PS::ParticleSystem<RealPtcl> ptcl;
	ptcl.initialize();
	//
	typedef GI_init<RealPtcl> PROBLEM;
	//
	PROBLEM::setupIC(ptcl, dinfo);
	std::cerr << ptcl.getNumberOfParticleGlobal() << std::endl;
	for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
		ptcl[i].setPressure();
		ptcl[i].initialize();
	}

	int Nlf, Ngrp;
	Nlf = 16;
	Ngrp = 1024;
	const int loop = 16;

	dinfo.initialize(0.3);
	dinfo.collectSampleParticle(ptcl);
	PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather dens_tree;
	dens_tree.initialize(ptcl.getNumberOfParticleLocal(), 0.5, Nlf, Ngrp);
	PS::TreeForForceShort<RESULT::Drvt, EPI::Drvt, EPJ::Drvt>::Gather drvt_tree;
	drvt_tree.initialize(ptcl.getNumberOfParticleLocal(), 0.5, Nlf, Ngrp);
	PS::TreeForForceShort<RESULT::Hydr, EPI::Hydr, EPJ::Hydr>::Symmetry hydr_tree;
	hydr_tree.initialize(ptcl.getNumberOfParticleLocal(), 0.5, Nlf, Ngrp);

	for(int l = 0 ; l < loop ; ++ l){
		if(l % 4 == 0){
			time[0] = get_dtime();
			dinfo.decomposeDomainAll(ptcl);
			ptcl.exchangeParticle(dinfo);
			time[0] = PS::Comm::getMaxValue((double)(get_dtime() - time[0]));
			//Density
			time[1] = get_dtime();
			#ifdef ENABLE_PEZY
				dens_tree.calcForceAllAndWriteBackMultiWalk(DensDispatchKernel, DensRetrieveKernel, 1, ptcl, dinfo, N_WALK_LIMIT, true);
			#elif ENABLE_GPU
				dens_tree.calcForceAllAndWriteBackMultiWalk(DensDispatchKernel, DensRetrieveKernel, 1, ptcl, dinfo, N_WALK_LIMIT, true);
			#else
				dens_tree.calcForceAllAndWriteBack(CalcDensity(), ptcl, dinfo, PS::MAKE_LIST_FOR_REUSE);
			#endif
			time[1] = PS::Comm::getMaxValue((double)(get_dtime() - time[1]));
			for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
				ptcl[i].setPressure();
				ptcl[i].initialize();
			}
			//Drvt
			time[2] = get_dtime();
			#ifdef ENABLE_PEZY
				drvt_tree.calcForceAllAndWriteBackMultiWalk(DrvtDispatchKernel, DrvtRetrieveKernel, 1, ptcl, dinfo, N_WALK_LIMIT, true);
			#elif ENABLE_GPU
				drvt_tree.calcForceAllAndWriteBackMultiWalk(DrvtDispatchKernel, DrvtRetrieveKernel, 1, ptcl, dinfo, N_WALK_LIMIT, true);
			#else
				drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), ptcl, dinfo, PS::MAKE_LIST_FOR_REUSE);
			#endif
			time[2] = PS::Comm::getMaxValue((double)(get_dtime() - time[2]));
			//Force
			time[3] = get_dtime();
			#ifdef ENABLE_PEZY
				hydr_tree.calcForceAllAndWriteBackMultiWalk(HydrDispatchKernel, HydrRetrieveKernel, 1, ptcl, dinfo, N_WALK_LIMIT, true);
			#elif ENABLE_GPU
				hydr_tree.calcForceAllAndWriteBackMultiWalk(HydrDispatchKernel, HydrRetrieveKernel, 1, ptcl, dinfo, N_WALK_LIMIT, true);
			#else
				hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), ptcl, dinfo, PS::MAKE_LIST_FOR_REUSE);
			#endif
			time[3] = PS::Comm::getMaxValue((double)(get_dtime() - time[3]));
			//time[0] = PS::Comm::getMaxValue((double)(get_dtime() - time[0]));
			if(PS::Comm::getRank() == 0){
				std::cout << PS::Comm::getNumberOfProc() << " " << l << " " << time[0] << " " << time[1] << " " << time[2] << " " << time[3] << std::endl;
			}
		}else{
			time[0] = get_dtime();
			time[0] = PS::Comm::getMaxValue((double)(get_dtime() - time[0]));
			//Density
			time[1] = get_dtime();
			#ifdef ENABLE_PEZY
				dens_tree.calcForceAllAndWriteBackMultiWalk(DensDispatchKernel, DensRetrieveKernel, 1, ptcl, dinfo, N_WALK_LIMIT, true);
			#elif ENABLE_GPU
				dens_tree.calcForceAllAndWriteBackMultiWalk(DensDispatchKernel, DensRetrieveKernel, 1, ptcl, dinfo, N_WALK_LIMIT, true);
			#else
				dens_tree.calcForceAllAndWriteBack(CalcDensity(), ptcl, dinfo, PS::REUSE_LIST);
			#endif
			time[1] = PS::Comm::getMaxValue((double)(get_dtime() - time[1]));
			for(int i = 0 ; i < ptcl.getNumberOfParticleLocal() ; ++ i){
				ptcl[i].setPressure();
				ptcl[i].initialize();
			}
			//Drvt
			time[2] = get_dtime();
			#ifdef ENABLE_PEZY
				drvt_tree.calcForceAllAndWriteBackMultiWalk(DrvtDispatchKernel, DrvtRetrieveKernel, 1, ptcl, dinfo, N_WALK_LIMIT, true);
			#elif ENABLE_GPU
				drvt_tree.calcForceAllAndWriteBackMultiWalk(DrvtDispatchKernel, DrvtRetrieveKernel, 1, ptcl, dinfo, N_WALK_LIMIT, true);
			#else
				drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), ptcl, dinfo, PS::REUSE_LIST);
			#endif
			time[2] = PS::Comm::getMaxValue((double)(get_dtime() - time[2]));
			//Force
			time[3] = get_dtime();
			#ifdef ENABLE_PEZY
				hydr_tree.calcForceAllAndWriteBackMultiWalk(HydrDispatchKernel, HydrRetrieveKernel, 1, ptcl, dinfo, N_WALK_LIMIT, true);
			#elif ENABLE_GPU
				hydr_tree.calcForceAllAndWriteBackMultiWalk(HydrDispatchKernel, HydrRetrieveKernel, 1, ptcl, dinfo, N_WALK_LIMIT, true);
			#else
				hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), ptcl, dinfo, PS::REUSE_LIST);
			#endif
			time[3] = PS::Comm::getMaxValue((double)(get_dtime() - time[3]));
			//time[0] = PS::Comm::getMaxValue((double)(get_dtime() - time[0]));
			if(PS::Comm::getRank() == 0){
				std::cout << PS::Comm::getNumberOfProc() << " " << l << " " << time[0] << " " << time[1] << " " << time[2] << " " << time[3] << std::endl;
			}
		}
	}



	PS::Finalize();
	return 0;
}

