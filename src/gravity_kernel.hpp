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
        func(epi, ni, epj, nj, force);

        
        for (PS::S32 i=0; i<ni; i++) if ( force[i].neighbor.number > 2 ) {
                PS::S32 ii=0;
                for (PS::S32 j=0; j<nj; j++) {
                    if ( epi[i].id_local == epj[j].id_local
                         && epi[i].myrank == epj[j].myrank ) continue;

                    PS::F32vec xi       = epi[i].pos - epi[0].pos;
                    PS::F32    rsearchi = epi[i].r_search;
                    
                    PS::F32vec xj       = epj[j].pos - epi[0].pos;
                    PS::F32    mj       = epj[j].mass;
                    PS::F32    rsearchj = epj[j].r_search;

                    PS::F32    rsearch  = std::max(rsearchi, rsearchj);
                    PS::F32    rsearch2 = rsearch * rsearch * 1.1025f;
                    
                    PS::F32vec rij = xj - xi;
                    PS::F32    dr2 = rij * rij + (PS::F32)FP_t::eps2;
                    
                    if ( dr2 < rsearch2 ) {
                        if (ii<8) {
                            force[i].neighbor.id_local[ii] = epj[j].id_local;
                            force[i].neighbor.rank[ii]     = epj[j].myrank;
                        }
                        ii ++;
                    }
                }
                assert ( force[i].neighbor.number == ii + 1 );
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



