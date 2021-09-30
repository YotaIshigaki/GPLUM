#pragma once

#include<particle_simulator.hpp>
#include "user-defined.hpp"

#ifdef SUNWAY
extern "C"{
   #include <athread.h>
   #include "cpe_func.h"
   extern void SLAVE_FUN(ForceKernelSunWay1st)(void*);
   extern void SLAVE_FUN(ForceKernelSunWay2nd)(void*);
}
Force_Kernel_Pars cpe_pars; // Force_Kernel_Pars is defined in cpe_func.h
#endif
#ifdef MPE_FORCE_KERNEL
#include "mpe_func.h"
MPE_pars mpe_pars;
#endif

PS::F64 wtime_copy_all_epj   = 0.0;
PS::F64 wtime_copy_all_spj   = 0.0;
PS::F64 wtime_calc_pj    = 0.0;
PS::F64 wtime_calc_epj   = 0.0;
PS::F64 wtime_calc_spj   = 0.0;
//PS::F64 wtime_dispatch   = 0.0;
//PS::F64 wtime_retrieve   = 0.0;
PS::F64 wtime_copy_ptcl_to_mm = 0.0;
PS::F64 wtime_copy_force_from_mm = 0.0;
PS::F64 wtime_kick = 0.0;
PS::F64 wtime_drift = 0.0;
PS::F64 wtime_calc_energy = 0.0;

PS::F64 ENERGY_TOT, ENERGY_KIN, ENERGY_POT;
bool DISPATCH_MODE;


extern PS::F64 DT_TREE;
extern Planet PLANET;
extern PS::S32 N_SATELLITE;
extern PS::ReallocatableArray<Satellite> SATELLITE;
extern PS::ReallocatableArray<SatelliteForce> SATELLITE_FORCE;
extern PS::ReallocatableArray<SatelliteForce> SATELLITE_FORCE_LOC;
extern std::ofstream fout_debug;

TreeType * TREE_POINTER;

#if !defined(SUNWAY) || !defined(SUNWAY_FORCE_KERNEL)
const PS::S64 N_EPJ_LIMIT = 20000000;
const PS::S64 N_SPJ_LIMIT = 10000000;

//Epi EPI_ARRAY[N_EPI_LIMIT];
Epj EPJ_ARRAY[N_EPJ_LIMIT];
PS::SPJMonopoleScatter SPJ_ARRAY[N_SPJ_LIMIT];
#endif
const PS::S64 N_WALK_LIMIT = 500;
const PS::S64 N_EPI_LIMIT = 2000000;
void * EPI_POINTER;
Force FORCE_ARRAY[N_EPI_LIMIT];

template<class Tptcl, class Tptclforce, class Tsat, class Tsatforce>
void CalcEnergy(const Tptcl ptcl[],
                const Tptclforce ptcl_force[],
                const PS::S32 n_ptcl,
                const Tsat satellite[],
                const Tsatforce satellite_force[],
                const PS::S32 n_sat,
                PS::F64 & etot,
                PS::F64 & ekin,
                PS::F64 & epot){
    etot = ekin = epot = 0.0;
    PS::F64 etot_loc = 0.0;
    PS::F64 ekin_loc = 0.0;
    PS::F64 epot_loc = 0.0;
    for(PS::S32 i=0; i<n_ptcl; i++){
        ekin_loc += ptcl[i].mass * ptcl[i].vel * ptcl[i].vel;
        epot_loc += ptcl[i].mass * ptcl_force[i].pot;
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    etot_loc  = ekin_loc + epot_loc;
    etot = PS::Comm::getSum(etot_loc);
    epot = PS::Comm::getSum(epot_loc);
    ekin = PS::Comm::getSum(ekin_loc);
    //std::cerr<<"ekin= "<<ekin<<" epot= "<<epot<<" etot= "<<etot<<std::endl;
    for(PS::S32 i=0; i<n_sat; i++){
        ekin += 0.5 * satellite[i].mass * satellite[i].vel * satellite[i].vel;
        epot += 0.5 * satellite[i].mass * satellite_force[i].pot;
    }
    etot = ekin + epot;
}


template<class Tpi, class Tpj, class Tforce>
void CalcGravPair(const Tpi & pi,
                  const Tpj & pj,
                  Tforce & force,
                  const PS::F64 eps_sq){
    const PS::F64vec rij = pi.pos - pj.pos;
    const PS::F64 mj = pj.mass;
    const PS::F64 r_sq = rij * rij + eps_sq;
    const PS::F64 r_inv = 1.0 / sqrt(r_sq);
    const PS::F64 pij = mj * r_inv;
    const double mri3 = r_inv*r_inv*pij;
    force.acc -= mri3 * rij;
    force.pot -= pij;
}

template<class Tpi, class Tpj, class Tforce>
void CalcGravEx(const Tpi & pi,
                const Tpj & pj,
                Tforce & force,
                const PS::F64 eps_sq){
    const PS::F64vec rij = pi.pos - pj.pos;
    const PS::F64 mj = pj.mass;
    const PS::F64 r_sq = rij * rij + eps_sq;
    const PS::F64 r_inv = 1.0 / sqrt(r_sq);
    const PS::F64 pij = mj * r_inv;
    const double mri3 = r_inv*r_inv*pij;
    force.acc -= mri3 * rij;
    force.pot -= 2.0*pij; // factor 2 is just a trick
}

template<class Tpi, class Tpj, class Tforce>
void CalcGravFormArray(const Tpi & pi,
                       const Tpj & pj,
                       Tforce & force,
                       const PS::F64 eps_sq){
    const PS::S32 nj = pj.size();
    for(PS::S32 j=0; j<nj; j++){
        const PS::F64 mj = pj[j].mass;
        const PS::F64vec rij = pi.pos - pj[j].pos;
        const PS::F64 r_sq = rij * rij + eps_sq;
        const PS::F64 r_inv = 1.0 / sqrt(r_sq);
        const PS::F64 pij = mj * r_inv;
        const double mri3 = r_inv*r_inv*pij;
        force.acc -= mri3 * rij;
        force.pot -= pij;
    }
}


template<class Tepi, class Tepj, class Tspj>
PS::S32 DispatchKernelWithSP(
                             const PS::S32   n_walk,
                             Tepi      epi[],
                             const PS::S32   n_epi[],
                             const PS::S32   n_disp_epi[],                             
                             const PS::S32   adr_epj[],
                             const PS::S32   n_epj[],
                             const PS::S32   n_disp_epj[],
                             const PS::S32   adr_spj[],
                             const PS::S32   n_spj[],
                             const PS::S32   n_disp_spj[],
                             const Tepj      epj[],
                             const Tspj      spj[]){

    //******** debug [start] *********
    /*
    for(PS::S32 ig=0; ig<n_walk; ig++){
        PS::F64 mass = 0.0;
        for(PS::S32 jp=0; jp<n_epj[ig]; jp++){
            const PS::S32 adr_des = n_disp_epj[ig]+jp;
            const PS::S32 adr_src = adr_epj[adr_des];
            mass += epj[adr_src].mass;
        }
        for(PS::S32 jp=0; jp<n_spj[ig]; jp++){
            const PS::S32 adr_des = n_disp_spj[ig]+jp;
            const PS::S32 adr_src = adr_spj[adr_des];            
            mass += spj[adr_src].mass;
        }
        if(PS::Comm::getRank()==0)std::cerr<<"mass= "<<mass<<std::endl;
    } 
    */    
#ifdef CHECK_NAN
    PS::S32 num_nan_loc,num_nan;
    PS::S32 n_epi_tot,n_epj_tot,n_spj_tot;
    static int num_called=0;
    num_called++;
    if (PS::Comm::getRank() == 0)
        std::cout << "num_called = " << num_called << std::endl;

    n_epi_tot = 0;
    n_epj_tot = 0;
    n_spj_tot = 0;
    for (PS::S32 ig=0; ig<n_walk; ig++) {
       n_epi_tot += n_epi[ig];
       n_epj_tot += n_epj[ig];
       n_spj_tot += n_spj[ig];
    }
    PS::S32 n_epi_tot_min,n_epi_tot_max;
    PS::S32 n_epj_tot_min,n_epj_tot_max;
    PS::S32 n_spj_tot_min,n_spj_tot_max;
    PS::S32 n_jp_tot_min,n_jp_tot_max;
    PS::S32 ntmp;
    MPI::COMM_WORLD.Allreduce(&n_epi_tot,&n_epi_tot_min,1,MPI::INT,MPI::MIN);
    MPI::COMM_WORLD.Allreduce(&n_epi_tot,&n_epi_tot_max,1,MPI::INT,MPI::MAX);
    MPI::COMM_WORLD.Allreduce(&n_epj_tot,&n_epj_tot_min,1,MPI::INT,MPI::MIN);
    MPI::COMM_WORLD.Allreduce(&n_epj_tot,&n_epj_tot_max,1,MPI::INT,MPI::MAX);
    MPI::COMM_WORLD.Allreduce(&n_spj_tot,&n_spj_tot_min,1,MPI::INT,MPI::MIN);
    MPI::COMM_WORLD.Allreduce(&n_spj_tot,&n_spj_tot_max,1,MPI::INT,MPI::MAX);
    ntmp = n_epj_tot + n_spj_tot;
    MPI::COMM_WORLD.Allreduce(&ntmp,&n_jp_tot_min,1,MPI::INT,MPI::MIN);
    MPI::COMM_WORLD.Allreduce(&ntmp,&n_jp_tot_max,1,MPI::INT,MPI::MAX);
    if (PS::Comm::getRank() == 0) {
        std::cout << "n_epi_tot_min = " << n_epi_tot_min << std::endl;
        std::cout << "n_epi_tot_max = " << n_epi_tot_max << std::endl;
        std::cout << "n_epj_tot_min = " << n_epj_tot_min << std::endl;
        std::cout << "n_epj_tot_max = " << n_epj_tot_max << std::endl;
        std::cout << "n_spj_tot_min = " << n_spj_tot_min << std::endl;
        std::cout << "n_spj_tot_max = " << n_spj_tot_max << std::endl;
        std::cout << "n_jp_tot_min  = " << n_jp_tot_min  << std::endl;
        std::cout << "n_jp_tot_max  = " << n_jp_tot_max  << std::endl;
    }

    num_nan_loc = 0;
    for (PS::S32 i=0; i<n_epi_tot; i++) {
        if ((std::isfinite(epi[i].pos.x) != true) ||
            (std::isfinite(epi[i].pos.y) != true) ||
            (std::isfinite(epi[i].pos.z) != true)) {
            num_nan_loc++;
        }
    } 
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "[nan-epi1] " << num_nan << std::endl;

    num_nan_loc = 0;
    for (PS::S32 ig=0; ig<n_walk; ig++) {
        for (PS::S32 jp=0; jp<n_epj[ig]; jp++) {
            const PS::S32 adr_des = n_disp_epj[ig]+jp;
            const PS::S32 adr_src = adr_epj[adr_des];
            if ((std::isfinite(epj[adr_src].pos.x) != true) ||
                (std::isfinite(epj[adr_src].pos.y) != true) ||
                (std::isfinite(epj[adr_src].pos.z) != true)) {
                num_nan_loc++;
            }
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "[nan-epj] " << num_nan << std::endl;
    
    num_nan_loc = 0;
    for (PS::S32 ig=0; ig<n_walk; ig++) {
        for (PS::S32 jp=0; jp<n_spj[ig]; jp++) {
            const PS::S32 adr_des = n_disp_spj[ig]+jp;
            const PS::S32 adr_src = adr_spj[adr_des];
            if ((std::isfinite(spj[adr_src].pos.x) != true) ||
                (std::isfinite(spj[adr_src].pos.y) != true) ||
                (std::isfinite(spj[adr_src].pos.z) != true)) {
                num_nan_loc++;
            }
        } 
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "[nan-spj] " << num_nan << std::endl;
    //if (num_nan_loc > 0) {
    //    std::cout << "myrank = " << PS::Comm::getRank() << " "
    //              << "ntot = " << n_epi_tot + ntmp << std::endl;
    //}
#endif
    //******** debug [end] *********

#if defined(SUNWAY) && defined(SUNWAY_FORCE_KERNEL)
    // [*** Note ****]
    //   Long-Wang's force kernel is implemented here.
    //   cpe_pars below is a global variable defined at nbody.cpp.
    //* Set force kernel parameters
    cpe_pars.epi        = (EpiMM*) epi;
    cpe_pars.epj        = (EpjMM*) epj;
    cpe_pars.spj        = (SpjMM*) spj;
    cpe_pars.sat        = (EpiMM*) SATELLITE.getPointer();
    cpe_pars.force_sat  = (ForceMM*) SATELLITE_FORCE_LOC.getPointer();
    cpe_pars.n_sat      = N_SATELLITE;
    cpe_pars.n_disp_epi = const_cast<int*>(n_disp_epi);
    cpe_pars.n_disp_epj = const_cast<int*>(n_disp_epj);
    cpe_pars.n_disp_spj = const_cast<int*>(n_disp_spj);
    cpe_pars.n_epi      = const_cast<int*>(n_epi);
    cpe_pars.n_epj      = const_cast<int*>(n_epj);
    cpe_pars.n_spj      = const_cast<int*>(n_spj);
    cpe_pars.adr_epj    = const_cast<int*>(adr_epj);
    cpe_pars.adr_spj    = const_cast<int*>(adr_spj);
    cpe_pars.n_walk     = n_walk;
#ifdef SUNWAY_FORCE_KERNEL_CHECK
    cpe_pars.forcei     = (ForceMM*)force_mm;
#endif
    //* Compute force
    PS::F64 starttime = MPI::Wtime();
    __real_athread_spawn((void*)slave_ForceKernelSunWay1st,(void*)&cpe_pars);
    athread_join();
    PS::wtime_calc_force += MPI::Wtime() - starttime;
    //* Reduce satellite force 
    MPI_Allreduce( (double*)(SATELLITE_FORCE_LOC.getPointer()),
                   (double*)(SATELLITE_FORCE.getPointer()),
                   4*SATELLITE.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //* Compute satellite-satellite force and perform kick-drift for satellite
    cpe_pars.sat        = (EpiMM*) SATELLITE.getPointer();
    cpe_pars.force_sat  = (ForceMM*) SATELLITE_FORCE.getPointer();
    cpe_pars.n_sat      = N_SATELLITE; 
    __real_athread_spawn((void*)slave_ForceKernelSunWay2nd,(void*)&cpe_pars);
    athread_join();

    //* Copy
    EPI_POINTER = epi;

    //* For test-run
    //fout_debug.close();
    //athread_halt();
    //PS::Finalize();
    //std::exit(0);
#elif defined(MPE_FORCE_KERNEL)
    EpiMM* sat_mpe = (EpiMM*)SATELLITE.getPointer();
    ForceMM* force_sat_mpe = (ForceMM*)SATELLITE_FORCE_LOC.getPointer();
    memset(force_sat_mpe,0,N_SATELLITE*sizeof(Force));
    std::cerr<<"Force calculation on MPE\n";
    //      wtime_copy0 = PS::GetWtime();
    for(int iw=0; iw<n_walk; iw++){
      Force* force_mpe = &FORCE_ARRAY[n_disp_epi[iw]];
      memset(force_mpe,0,n_epi[iw]*sizeof(Force));
      //      std::cerr<<"iw="<<iw<<" I-J Ni="<<n_epi[iw]<<" Nj="<<n_epj[iw]<<"\n";
      calcDirectGrav((EpiMM*)&(epi[n_disp_epi[iw]]), n_epi[iw], 
                     (EpjMM*)epj, n_epj[iw],(int*)&adr_epj[n_disp_epj[iw]],
                     (ForceMM*)force_mpe, mpe_pars.r_ep, mpe_pars.r_ep, mpe_pars.kappa, mpe_pars.eta);
      //      std::cerr<<"iw="<<iw<<" I-SP\n";
      calcSPGrav((EpiMM*)&(epi[n_disp_epi[iw]]),n_epi[iw],
                 (SpjMM*)spj, n_spj[iw],(int*)&adr_spj[n_disp_spj[iw]],
                 (ForceMM*)force_mpe);
      //      std::cerr<<"iw="<<iw<<" I-SAT\n";
      calcDirectGrav((EpiMM*)&(epi[n_disp_epi[iw]]), n_epi[iw],
                     (EpjMM*)sat_mpe, N_SATELLITE, NULL,
                     (ForceMM*)force_mpe, mpe_pars.r_ep, mpe_pars.r_sat, mpe_pars.kappa, mpe_pars.eta, false);
      //      std::cerr<<"iw="<<iw<<" SAT-I\n";
      calcDirectGrav((EpiMM*)sat_mpe, N_SATELLITE,
                     (EpjMM*)&(epi[n_disp_epi[iw]]), n_epi[iw], NULL,
                     (ForceMM*)force_sat_mpe, mpe_pars.r_ep, mpe_pars.r_sat, mpe_pars.kappa, mpe_pars.eta, false);
      //      std::cerr<<"iw="<<iw<<" I-P\n";
      calcForcePlanet((EpiMM*)&(epi[n_disp_epi[iw]]), n_epi[iw], PLANET.mass, (ForceMM*)force_mpe);
    }
    //    std::cerr<<"Sat MPI Allreduce\n";
    MPI_Allreduce( (double*)(SATELLITE_FORCE_LOC.getPointer()),
                   (double*)(SATELLITE_FORCE.getPointer()),
                   4*SATELLITE.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    //    std::cerr<<"Sat-Sat\n";
    calcDirectGrav((EpiMM*)sat_mpe, N_SATELLITE,
                   (EpjMM*)sat_mpe, N_SATELLITE, NULL,
                   (ForceMM*)force_sat_mpe, mpe_pars.r_ep, mpe_pars.r_sat, mpe_pars.kappa, mpe_pars.eta,false);
    //    std::cerr<<"Sat-Planet\n";
    calcForcePlanet(sat_mpe, N_SATELLITE, PLANET.mass, force_sat_mpe);

    //      wtime_copy0 = PS::GetWtime() - wtime_copy0;
    //      fprintf(stderr,"host force time per loop=%e \n",wtime_copy0);
    
#else
    PS::F64 wtime_offset_out = PS::GetWtime();
    const PS::F64 eps_sq = Tepi::eps * Tepi::eps;
    assert(N_EPI_LIMIT > n_disp_epi[n_walk]);
    assert(N_EPJ_LIMIT > n_disp_epj[n_walk]);
    assert(N_SPJ_LIMIT > n_disp_spj[n_walk]);
    const PS::S32 n_epi_tot = n_disp_epi[n_walk];
#ifndef CANCEL_FORCE    
    /////////////////////
    // force on satellite
    for(PS::S32 ip=0; ip<SATELLITE.size(); ip++){
        //SATELLITE_FORCE[ip].clear();
        SATELLITE_FORCE_LOC[ip].clear();
        for(PS::S32 jp=0; jp<n_epi_tot; jp++){
            CalcGravPair(SATELLITE[ip], epi[jp], SATELLITE_FORCE_LOC[ip], eps_sq);
        }
    }
    MPI_Allreduce( (double*)(SATELLITE_FORCE_LOC.getPointer()),
                   (double*)(SATELLITE_FORCE.getPointer()),
                   4*SATELLITE.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for(PS::S32 ip=0; ip<SATELLITE.size(); ip++){
        for(PS::S32 jp=0; jp<SATELLITE.size(); jp++){
            if(ip == jp) continue;
            CalcGravPair(SATELLITE[ip], SATELLITE[jp], SATELLITE_FORCE[ip], eps_sq);
        }
        CalcGravEx(SATELLITE[ip], PLANET, SATELLITE_FORCE[ip], 0.0);
    }
    // force on satellite
    /////////////////////
#endif
    

    /////////////////////
    // force on particle
    for(PS::S32 ig=0; ig<n_walk; ig++){
        PS::F64 mass = 0.0;
        PS::F64vec cm_pos = 0.0;
        //PS::F64 wtime_offset = PS::GetWtime();
        for(PS::S32 jp=0; jp<n_epj[ig]; jp++){
            const PS::S32 adr_des = n_disp_epj[ig]+jp;
            const PS::S32 adr_src = adr_epj[adr_des];
            EPJ_ARRAY[adr_des].id = epj[adr_src].id;
            EPJ_ARRAY[adr_des].mass = epj[adr_src].mass;
            EPJ_ARRAY[adr_des].pos = epj[adr_src].pos;
            mass += EPJ_ARRAY[adr_des].mass;
            cm_pos += EPJ_ARRAY[adr_des].mass*EPJ_ARRAY[adr_des].pos;
        }
        //wtime_copy_all_epj += PS::GetWtime() - wtime_offset;
        //wtime_offset = PS::GetWtime();
        for(PS::S32 jp=0; jp<n_spj[ig]; jp++){
            const PS::S32 adr_des = n_disp_spj[ig]+jp;
            const PS::S32 adr_src = adr_spj[adr_des];
            SPJ_ARRAY[adr_des].mass = spj[adr_src].mass;
            SPJ_ARRAY[adr_des].pos  = spj[adr_src].pos;
            mass += SPJ_ARRAY[adr_des].mass;
            cm_pos += SPJ_ARRAY[adr_des].mass*SPJ_ARRAY[adr_des].pos;
        }
        //wtime_copy_all_spj += PS::GetWtime() - wtime_offset;        
        //if(PS::Comm::getRank() == 0) std::cerr<<"mass= "<<mass<<std::endl;
        //wtime_offset = PS::GetWtime();
        for(PS::S32 ip=0; ip<n_epi[ig]; ip++){
            mass = 0.0;
            const PS::S32 adr_i = n_disp_epi[ig] + ip;
            const Tepi & iptcl = epi[adr_i];
            Force & iforce = FORCE_ARRAY[adr_i];
            iforce.acc = 0.0;
            iforce.pot = 0.0;
#ifndef CANCEL_FORCE            
            for(PS::S32 jp=0; jp<n_epj[ig]; jp++){
                const PS::S32 adr_j = n_disp_epj[ig] + jp;
                if(iptcl.id == EPJ_ARRAY[adr_j].id) continue;
                //mass += EPJ_ARRAY[adr_j].mass;
                CalcGravPair(iptcl, EPJ_ARRAY[adr_j], iforce, eps_sq);
            }
            for(PS::S32 jp=0; jp<n_spj[ig]; jp++){
                const PS::S32 adr_j = n_disp_spj[ig] + jp;
                //mass += SPJ_ARRAY[adr_j].mass;

                CalcGravPair(iptcl, SPJ_ARRAY[adr_j], iforce, eps_sq);

            }
            for(PS::S32 jp=0; jp<SATELLITE.size(); jp++){
                CalcGravPair(iptcl, SATELLITE[jp], iforce, eps_sq);
            }
#endif
            CalcGravEx(iptcl, PLANET, iforce, 0.0);
        }
        //wtime_calc_pj += PS::GetWtime() - wtime_offset;
    }

    //wtime_dispatch += PS::GetWtime() - wtime_offset_out;
#endif // End of SUNWAY_FORCE_KERNEL

    //****** debug[start] *******
#ifdef CHECK_NAN
    num_nan_loc = 0;
    for (PS::S32 i=0; i<n_epi_tot; i++) {
        if ((std::isfinite(epi[i].pos.x) != true) ||
            (std::isfinite(epi[i].pos.y) != true) ||
            (std::isfinite(epi[i].pos.z) != true)) {
            num_nan_loc++;
        }
    } 
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0) 
        std::cout << "[nan-epi2] " << num_nan << std::endl;
#endif
    //******* debug[end] *******


}

template<class Tforce>
PS::S32 RetrieveKernel(const PS::S32 ni_tot,
                       Tforce force[]){
    /*
    //PS::F64 wtime_offset = PS::GetWtime();
    for(PS::S32 i=0; i<ni_tot; i++){
        force[i].acc = FORCE_ARRAY[i].acc;
        force[i].pot = FORCE_ARRAY[i].pot;
    }
    //wtime_retrieve += PS::GetWtime() - wtime_offset;
    */
    return 0;
}



template<class Tepi, class Tepj, class Tspj>
PS::S32 DispatchKernelWithSPKickDrift(
                                      const PS::S32   n_walk,
                                      Tepi      epi[],
                                      const PS::S32   n_epi[],
                                      const PS::S32   n_disp_epi[],
                                      const PS::S32   adr_epj[],
                                      const PS::S32   n_epj[],
                                      const PS::S32   n_disp_epj[],
                                      const PS::S32   adr_spj[],
                                      const PS::S32   n_spj[],
                                      const PS::S32   n_disp_spj[],
                                      const Tepj      epj[],
                                      const Tspj      spj[]){
    
    
#if defined(SUNWAY) && defined(SUNWAY_FORCE_KERNEL)
    DispatchKernelWithSP(n_walk, epi, n_epi, n_disp_epi,
                         adr_epj, n_epj, n_disp_epj,
                         adr_spj, n_spj, n_disp_spj,
                         epj, spj);
#else
    DispatchKernelWithSP(n_walk, epi, n_epi, n_disp_epi,
                         adr_epj, n_epj, n_disp_epj,
                         adr_spj, n_spj, n_disp_spj,
                         epj, spj);
    const PS::S32 n_epi_tot = n_disp_epi[n_walk];
    const PS::F64 dt_kick = (DISPATCH_MODE == 0) ? DT_TREE : DT_TREE*0.5;
    PS::F64 wtime_offset = PS::GetWtime();
    for(PS::S32 ip=0; ip<n_epi_tot; ip++){
        epi[ip].vel += dt_kick * FORCE_ARRAY[ip].acc;
    }
    for(PS::S32 ip=0; ip<SATELLITE.size(); ip++){
        SATELLITE[ip].vel += dt_kick * SATELLITE_FORCE[ip].acc;
    }
    wtime_kick += PS::GetWtime() - wtime_offset;
    if(DISPATCH_MODE == 1){
        wtime_offset = PS::GetWtime();
        CalcEnergy(epi, FORCE_ARRAY, n_epi_tot,
                   SATELLITE.getPointer(), SATELLITE_FORCE.getPointer(), N_SATELLITE,
                   ENERGY_TOT, ENERGY_KIN, ENERGY_POT);
        wtime_calc_energy += PS::GetWtime() - wtime_offset;
        wtime_offset = PS::GetWtime();
        for(PS::S32 ip=0; ip<n_epi_tot; ip++){
            epi[ip].vel += dt_kick * FORCE_ARRAY[ip].acc;
        }
        for(PS::S32 ip=0; ip<SATELLITE.size(); ip++){
            SATELLITE[ip].vel += dt_kick * SATELLITE_FORCE[ip].acc;
        }
        wtime_kick += PS::GetWtime() - wtime_offset;
    }
    // TO DO; collision procedure
    wtime_offset = PS::GetWtime();    
    for(PS::S32 ip=0; ip<n_epi_tot; ip++){
        epi[ip].pos += DT_TREE * epi[ip].vel;
        //EPI_ARRAY[ip].pos = epi[ip].pos;
    }
    EPI_POINTER = epi;
    for(PS::S32 ip=0; ip<SATELLITE.size(); ip++){
        SATELLITE[ip].pos += DT_TREE * SATELLITE[ip].vel;
    }
    wtime_drift += PS::GetWtime() - wtime_offset;
#endif
}
