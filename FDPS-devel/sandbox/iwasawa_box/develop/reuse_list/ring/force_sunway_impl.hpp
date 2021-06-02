#pragma once

#include<particle_simulator.hpp>

extern "C"{
#include"kernel_sunway.h"
}

const PS::S64 N_EPJ_RESERVE_LIMIT = 10000000;
const PS::S64 N_SPJ_RESERVE_LIMIT = 10000000;

int         disp_mm[N_WALK_LIMIT+2][3]; // 0:epi, 1:epj, 2:spj
EpiMM       epi_mm[NI_LIMIT];
EpjMM       epj_mm[NJ_LIMIT];
SpjMM       spj_mm[NJ_LIMIT];
ForceMM     force_mm[NI_LIMIT];
double      r_coll_sq_mm;

void      * EPJ_POINTER;

template<class Tepj, class Tspj>
class ReserverJptcl{
public:
    PS::ReallocatableArray<Tepj> epj_reserver;
    PS::ReallocatableArray<Tspj> spj_reserver;
    void copy(const Tepj _epj[], const PS::S32 _n_epj,
              const Tspj _spj[], const PS::S32 _n_spj){
        epj_reserver.resizeNoInitialize(_n_epj);
        spj_reserver.resizeNoInitialize(_n_spj);
        for(PS::S32 i=0; i<_n_epj; i++) epj_reserver[i] = _epj[i];
        for(PS::S32 i=0; i<_n_spj; i++) spj_reserver[i] = _spj[i];
    }
};

template<class Tepi, class Tepj, class Tspj>
PS::S32 DispatchKernelWithSP(
                             const PS::S32    tag,
                             const PS::S32    n_walk,
                             const Tepi      *epi[],
                             const PS::S32    n_epi[],
                             const PS::S32   *id_epj[],
                             const PS::S32    n_epj[],
                             const PS::S32   *id_spj[],
                             const PS::S32    n_spj[],
                             const Tepj       epj[],
                             const PS::S32    n_send_epj,
                             const Tspj       spj[],
                             const PS::S32    n_send_spj,
                             const bool       flag_send_ptcl){
    static ReserverJptcl<Tepj, Tspj> reserver;
    if(flag_send_ptcl){
        reserver.copy(epj, n_send_epj, spj, n_send_spj);
        return 0;
    }
    EPJ_POINTER =  (Tepj*)reserver.epj_reserver.getPointer();
    
    const float eps2 = Tepi::eps * Tepi::eps;
    r_coll_sq_mm = 2.0*Tepj::r_coll*2.0*Tepj::r_coll;
    disp_mm[0][0] = disp_mm[0][1] = disp_mm[0][2] = 0;
    int n_epi_cum = 0;
    int n_epj_cum = 0;
    int n_spj_cum = 0;
    for(int iw=0; iw<n_walk; iw++){
        PS::F64 mass = 0.0;        
        disp_mm[iw+1][0] = disp_mm[iw][0] + n_epi[iw];
        disp_mm[iw+1][1] = disp_mm[iw][1] + n_epj[iw];
        disp_mm[iw+1][2] = disp_mm[iw][2] + n_spj[iw];
        for(int ip=0; ip<n_epi[iw]; ip++){
            epi_mm[n_epi_cum].id = epi[iw][ip].id;
            epi_mm[n_epi_cum].x  = epi[iw][ip].pos.x;
            epi_mm[n_epi_cum].y  = epi[iw][ip].pos.y;
            epi_mm[n_epi_cum].z  = epi[iw][ip].pos.z;
            n_epi_cum++;
        }
        for(int jp=0; jp<n_epj[iw]; jp++){
            const PS::S32 adr = id_epj[iw][jp];
            epj_mm[n_epj_cum].id   = reserver.epj_reserver[adr].id;
            epj_mm[n_epj_cum].mass = reserver.epj_reserver[adr].mass;
            epj_mm[n_epj_cum].x    = reserver.epj_reserver[adr].pos.x;
            epj_mm[n_epj_cum].y    = reserver.epj_reserver[adr].pos.y;
            epj_mm[n_epj_cum].z    = reserver.epj_reserver[adr].pos.z;
            mass += epj_mm[n_epj_cum].mass;
            n_epj_cum++;
        }
        for(int jp=0; jp<n_spj[iw]; jp++){
            const PS::S32 adr = id_spj[iw][jp];
            spj_mm[n_spj_cum].mass = reserver.spj_reserver[adr].mass;
            spj_mm[n_spj_cum].x    = reserver.spj_reserver[adr].pos.x;
            spj_mm[n_spj_cum].y    = reserver.spj_reserver[adr].pos.y;
            spj_mm[n_spj_cum].z    = reserver.spj_reserver[adr].pos.z;
            mass += spj_mm[n_spj_cum].mass;
            //std::cerr<<"spj_mm[n_spj_cum].mass= "<<spj_mm[n_spj_cum].mass<<std::endl;
            n_spj_cum++;
        }
        //std::cerr<<"mass= "<<mass<<std::endl; // for debuging
    }
    disp_mm[n_walk+1][0] = disp_mm[n_walk][0];
    disp_mm[n_walk+1][1] = disp_mm[n_walk][1];
    disp_mm[n_walk+1][2] = disp_mm[n_walk][2];
    DispatchKernelSunWay(n_walk, eps2); // after puting them to disp_lm, we can use MPE
    return 0;
}

template<class Tforce>
PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 ni[],
                       Tforce    *force[]){
    int ni_tot = 0;
    for(int k=0; k<n_walk; k++) ni_tot += ni[k];
    RetrieveKernelSunWay(ni_tot);
    int n_cnt = 0;
    for(int iw=0; iw<n_walk; iw++){
        for(int i=0; i<ni[iw]; i++){
            force[iw][i].acc.x = force_mm[n_cnt].ax;
            force[iw][i].acc.y = force_mm[n_cnt].ay;
            force[iw][i].acc.z = force_mm[n_cnt].az;
            force[iw][i].pot   = force_mm[n_cnt].pot;
            force[iw][i].n_coll   = force_mm[n_cnt].n_coll;
            force[iw][i].r_ngb_sq = force_mm[n_cnt].r_ngb_sq;
            const int adr_ngb = force_mm[n_cnt].adr_ngb;
            force[iw][i].adr_ngb = force_mm[n_cnt].adr_ngb;
            n_cnt++;
        }
    }
    return 0;
}
