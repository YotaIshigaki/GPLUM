#include "header.h"
#include "force-k2.h"
const PS::F64 EXPAND = 1.3;
int main(int argc, char* argv[]){
    PS::Initialize(argc, argv);
    PS::ParticleSystem<RealPtcl> sph_system;
    sph_system.initialize();
    sph_system.setAverageTargetNumberOfSampleParticlePerProcess(200);
    PS::DomainInfo dinfo;
    const PS::F32 coef_ema = 0.2;
    dinfo.initialize(coef_ema);
    system_t sysinfo;
    SetupIC(sph_system, sysinfo, dinfo);
    Initialize(sph_system);
    dinfo.decomposeDomainAll(sph_system);
    sph_system.exchangeParticle(dinfo);
    PS::TreeForForceShort<RESULT::Dens , EPI::Dens , EPJ::Dens >::Gather   dens_tree;
    PS::TreeForForceShort<RESULT::Hydro, EPI::Hydro, EPJ::Hydro>::Symmetry hydr_tree;
    PS::TreeForForceLong <RESULT::Grav , EPI::Grav , EPJ::Grav >::Monopole grav_tree;
    const int Ngrp = 128;
    dens_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, 4, Ngrp);
    hydr_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, 4, Ngrp);
    grav_tree.initialize(sph_system.getNumberOfParticleGlobal(), 0.5, 8, 256);
    for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
	sph_system[i].Rsearch = kernel_t::supportRadius() * sph_system[i].smth * EXPAND;
    }
    for(bool repeat = true; repeat==true;){
        bool repeat_loc = false;
        repeat = false;
        dens_tree.calcForceAll(CalcDensity(), sph_system, dinfo);
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
    hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
    grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);
    sysinfo.dt = getTimeStepGlobal(sph_system);
    CalcExternalForce(sph_system, sysinfo);
    PS::S64 step = 0;
    for(sysinfo.time = 0 ; sysinfo.time < sysinfo.end_time ; sysinfo.time += sysinfo.dt, ++ step){
	InitialKick(sph_system, sysinfo);
	FullDrift(sph_system, sysinfo);
	Predict(sph_system, sysinfo);
    if(step < 12 || step % 4 == 0){
        dinfo.collectSampleParticle(sph_system);
        dinfo.decomposeDomainMultiStep();
    }
	sph_system.exchangeParticle(dinfo);
    for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
        sph_system[i].Rsearch = kernel_t::supportRadius() * sph_system[i].smth * EXPAND;
    }
    PS::S64 n_int_ep_ep = 0;
    for(bool repeat = true; repeat==true;){
        bool repeat_loc = false;
            repeat = false;
            dens_tree.calcForceAll(CalcDensity(), sph_system, dinfo);
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
        dens_tree.setNInteractionEPEP(n_int_ep_ep);
	CalcPressure(sph_system);
        for(int i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
            sph_system[i].setBalsalaSwitch();
            sph_system[i].Rsearch = sph_system[i].smth * kernel_t::supportRadius();
        }
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
	}
    }
    PS::Finalize();
    return 0;
}
