#pragma once

#include "gravity_kernel_epep.hpp"
#include "gravity_kernel_epsp.hpp"

struct calcForceEPEPWithSearch{
    calcForceEPEPWithSearch(){}
    
    void operator()(const EPI_t* epi,
                    const int    ni,
                    const EPJ_t* epj,
                    const int    nj,
                    Force_t*     force){
#ifdef USE_INDIVIDUAL_CUTOFF
        CalcForceLongEPEP func(FP_t::eps2);
#else
        CalcForceLongEPEP func(FP_t::eps2, FP_t::r_out, FP_t::r_search);
#endif
        for (PS::S32 i=0; i<ni; i++) force[i].neighbor.setNumber();
        func(epi, ni, epj, nj, force);
     
        for (PS::S32 i=0; i<ni; i++) {
            PS::S32 n_ngb = force[i].neighbor.getNumber();
            if ( n_ngb == 0 ) continue; 
            
            PS::S32 k = 1;
            while ( force[i].neighbor.id_local[k] != -1 && k<8 ) k++;
            if (k<8) {
                if ( n_ngb == 1 ) {
                    force[i].neighbor.id_local[k] = force[i].neighbor.id_local[0];
                    force[i].neighbor.rank[k]     = force[i].neighbor.rank[0];
                    k ++;
                    
                } else if ( n_ngb > 1 ) {
                    PS::F32vec xi = epi[i].pos - epi[0].pos;
                    
                    for (PS::S32 j=0; j<nj; j++) {
                        if ( epi[i].id_local == epj[j].id_local
                             && epi[i].myrank == epj[j].myrank ) continue;
                        
                        PS::F32vec xj = epj[j].pos - epi[0].pos;
                        
#ifdef USE_INDIVIDUAL_CUTOFF
                        PS::F32    rsearchi = epi[i].r_search;
                        PS::F32    rsearchj = epj[j].r_search;
                        PS::F32    rsearch  = std::max(rsearchi, rsearchj);             
#else
                        PS::F32    rsearch  = FP_t::r_search;
#endif
                        PS::F32    rsearch2 = rsearch * rsearch * SAFTY_FACTOR2;
                        
                        PS::F32vec rij = xj - xi;
                        PS::F32    dr2 = rij * rij + (PS::F32)FP_t::eps2;
                        
                        if ( dr2 < rsearch2 ) {
                            if (k<8) {
                                force[i].neighbor.id_local[k] = epj[j].id_local;
                                force[i].neighbor.rank[k]     = epj[j].myrank;
                            }
                            k ++;
                        }
                    }
                }
                force[i].neighbor.number = k-1;
                //assert( force[i].neighbor.number == k-1 );
            }
        }
    }
};

struct calcForceEPSP{
    calcForceEPSP(){}
    void operator()(const EPI_t* epi,
                    const int    ni,
                    const SPJ_t* epj,
                    const int    nj,
                    Force_t*     force){
        
        CalcForceLongEPSP func(FP_t::eps2);
        func(epi, ni, epj, nj, force);
    }
};



