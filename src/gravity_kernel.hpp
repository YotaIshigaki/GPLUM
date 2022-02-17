#pragma once

#ifdef USE_PIKG
#include "gravity_kernel_epep.hpp"
#include "gravity_kernel_epsp.hpp"
#endif

struct calcForceEPEPWithSearch{
    calcForceEPEPWithSearch(){}

#ifdef USE_PIKG
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
    }
#else
    void operator()(const EPI_t* epi,
                    const int    ni,
                    const EPJ_t* epj,
                    const int    nj,
                    Force_t*     force){
        
        const PS::F32 eps2 = FP_t::eps2;
#ifndef USE_INDIVIDUAL_CUTOFF
        PS::F32 rout    = FP_t::rout;
        PS::F32 rsearch = FP_t::rsearch;
#endif
        PS::F32 xj[nj];
        PS::F32 yj[nj];
        PS::F32 zj[nj];
        PS::F32 mj[nj];
#ifdef USE_INDIVIDUAL_CUTOFF
        PS::F32 routj[nj];
        PS::F32 rsearchj[nj];
#endif
        PS::S32 idj[nj];
        PS::S32 rankj[nj];
        for(int j=0; j<nj; j++){
            xj[j]    = epj[j].pos[0] - epi[0].pos[0];
            yj[j]    = epj[j].pos[1] - epi[0].pos[1];
            zj[j]    = epj[j].pos[2] - epi[0].pos[2];
            mj[j]    = epj[j].getCharge();
#ifdef USE_INDIVIDUAL_CUTOFF
            routj[j]    = epj[j].r_out;
            rsearchj[j] = epj[j].r_search;
#endif
            idj[j]   = epj[j].id_local;
            rankj[j] = epj[j].myrank;
        }
        
        for(int i=0; i<ni; i++){
            //const PS::F64vec xi = epi[i].pos;
            PS::F32    afx     = 0.0;
            PS::F32    afy     = 0.0;
            PS::F32    afz     = 0.0;
            PS::F32    potf    = 0.0;
            PS::S32    ngbf    = 0;
            PS::S32    rankf   = 0;
            PS::S32    id_maxf = -1;
            PS::S32    id_minf = std::numeric_limits<PS::S32>::max();
            
            PS::F32    xix   = epi[i].pos[0] - epi[0].pos[0];
            PS::F32    xiy   = epi[i].pos[1] - epi[0].pos[1];
            PS::F32    xiz   = epi[i].pos[2] - epi[0].pos[2];
#ifdef USE_INDIVIDUAL_CUTOFF
            PS::F32    routi    = epi[i].r_out;
            PS::F32    rsearchi = epi[i].r_search;
#endif
            PS::S32    idi   = epi[i].id_local;
            PS::S32    ranki = epi[i].myrank;
            
#pragma omp simd reduction(+:afx,afy,afz,potf,ngbf,rankf), reduction(max:id_maxf), reduction(min:id_minf)
            for(int j=0; j<nj; j++){
                PS::F32 rijx    = xix - xj[j];
                PS::F32 rijy    = xiy - yj[j];
                PS::F32 rijz    = xiz - zj[j];
#ifdef USE_INDIVIDUAL_CUTOFF
                PS::F32 rout    = std::max(routi,    routj[j]);
                PS::F32 rsearch = std::max(rsearchi, rsearchj[j]);
#endif
                PS::F32 rout2    = rout * rout;
                PS::F32 rsearch2 = rsearch * rsearch * 1.0201f;
                
                //if(idi == idj[j]) continue;
                PS::F32 r2_real = rijx*rijx + rijy*rijy + rijz*rijz + eps2;
                PS::F32 r2      = std::max(r2_real, rout2);
                PS::F32 r_inv   = 1.0f/sqrt(r2);
                PS::F32 tmp     = 3.0f - r2*(r_inv*r_inv);
                r_inv   *= (tmp * 0.5f);
                PS::F32 r2_inv  = r_inv * r_inv;
                PS::F32 pot     = r_inv * mj[j];

                afx    +=  pot * r2_inv * rijx;
                afy    +=  pot * r2_inv * rijy;
                afz    +=  pot * r2_inv * rijz;
                potf   +=  pot;                
                if ( r2_real < rsearch2 ) {
                    if ( idi != idj[j] || ranki != rankj[j] ) {
                        ngbf    += 1;
                        rankf   += std::abs(ranki - rankj[j]);
                        id_maxf  = std::max(id_maxf, idj[j]);
                        id_minf  = std::min(id_minf, idj[j]);
                    }
                }
            }
            force[i].acc -= PS::F64vec(afx, afy, afz);
            force[i].phi -= potf;
            force[i].neighbor.number += ngbf;
            force[i].neighbor.rank   += rankf;
            force[i].neighbor.id_max = std::max(id_maxf, force[i].neighbor.id_max);
            force[i].neighbor.id_min = std::min(id_minf, force[i].neighbor.id_min);
        }
    }
#endif
};

struct calcForceEPSP{
    calcForceEPSP(){}
#ifdef USE_PIKG
    void operator()(const EPI_t* epi,
                    const int    ni,
                    const SPJ_t* epj,
                    const int    nj,
                    Force_t*     force){
        
        CalcForceLongEPSP func(FP_t::eps2);
        func(epi, ni, epj, nj, force);
    }
#else
    void operator()(const EPI_t* epi,
                    const int    ni,
                    const SPJ_t* epj,
                    const int    nj,
                    Force_t*     force){
        
        const PS::F32 eps2 = FP_t::eps2;
        PS::F32 xjx[nj];
        PS::F32 xjy[nj];
        PS::F32 xjz[nj];
        PS::F32 mj[nj];
#ifdef USE_QUAD
        PS::F32 qjxx[nj];
        PS::F32 qjyy[nj];
        PS::F32 qjzz[nj];
        PS::F32 qjxy[nj];
        PS::F32 qjyz[nj];
        PS::F32 qjzx[nj];
        PS::F32 tr[nj];
#endif
        for(int j=0; j<nj; j++){
#ifdef USE_POLAR_COORDINATE
            xjx[j]  = epj[j].pos_car[0] - epi[0].pos[0];
            xjy[j]  = epj[j].pos_car[1] - epi[0].pos[1];
            xjz[j]  = epj[j].pos_car[2] - epi[0].pos[2];
#else
            xjx[j]  = epj[j].pos[0] - epi[0].pos[0];
            xjy[j]  = epj[j].pos[1] - epi[0].pos[1];
            xjz[j]  = epj[j].pos[2] - epi[0].pos[2];
#endif
            mj[j]   = epj[j].getCharge();
#ifdef USE_QUAD
            qjxx[j] = epj[j].quad.xx;
            qjyy[j] = epj[j].quad.yy;
            qjzz[j] = epj[j].quad.zz;
            qjxy[j] = epj[j].quad.xy;
            qjyz[j] = epj[j].quad.yz;
            qjzx[j] = epj[j].quad.xz;
            //tr[j]   = epj[j].quad.getTrace();
            tr[j]   = epj[j].quad.xx + epj[j].quad.yy + epj[j].quad.xx;
#endif
        }
        for(int i=0; i<ni; i++){
            PS::F32 xix = epi[i].pos[0] - epi[0].pos[0];
            PS::F32 xiy = epi[i].pos[1] - epi[0].pos[1];
            PS::F32 xiz = epi[i].pos[2] - epi[0].pos[2];
            
            PS::F32 afx  = 0.0;
            PS::F32 afy  = 0.0;
            PS::F32 afz  = 0.0;
            PS::F32 potf = 0.0;
#pragma clan loop unroll_count(6)
#pragma omp simd reduction(+:afx,afy,afz,potf)	    
            for(int j=0; j<nj; j++){
                PS::F32 rijx = xjx[j] - xix;
                PS::F32 rijy = xjy[j] - xiy;
                PS::F32 rijz = xjz[j] - xiz;
                PS::F32 r2   = rijx*rijx + rijy*rijy + rijz*rijz + eps2;
                PS::F32 r_inv = 1.0f/sqrt(r2);
                PS::F32 tmp    = 3.0f - r2*(r_inv*r_inv);
                r_inv *= (tmp * 0.5f);
                PS::F32 r2_inv = r_inv * r_inv;
#ifdef USE_QUAD
                PS::F32 r3_inv = r2_inv * r_inv;
                PS::F32 r4_inv = r2_inv * r2_inv;
                PS::F32 r5_inv = r2_inv * r3_inv;
                
                PS::F32 qxx = 3.0f * qjxx[j] - tr[j];
                PS::F32 qyy = 3.0f * qjyy[j] - tr[j];
                PS::F32 qzz = 3.0f * qjzz[j] - tr[j];
                PS::F32 qxy = 3.0f * qjxy[j];
                PS::F32 qyz = 3.0f * qjyz[j];
                PS::F32 qzx = 3.0f * qjzx[j];
                PS::F32 mtr = -(eps2 * tr[j]);
                    
                PS::F32 qrx =qxx*rijx + qxy*rijy + qzx*rijz;
                PS::F32 qry =qyy*rijy + qyz*rijz + qxy*rijx;
                PS::F32 qrz =qzz*rijz + qzx*rijx + qyz*rijy;
                PS::F32 rqr = mtr + qrx*rijx + qry*rijy + qrz*rijz;
                PS::F32 rqr_r4_inv = rqr * r4_inv;

                PS::F32 meff  =  mj[j] + 0.5f * rqr_r4_inv;
                PS::F32 meff3 = (mj[j] + 2.5f * rqr_r4_inv) * r3_inv;

                afx  += r5_inv * qrx - meff3 * rijx;
                afy  += r5_inv * qry - meff3 * rijy;
                afz  += r5_inv * qrz - meff3 * rijz;
                potf += meff * r_inv;
#else
                PS::F32 pot     = r_inv * mj[j];
                
                afx    +=  -pot * r2_inv * rijx;
                afy    +=  -pot * r2_inv * rijy;
                afz    +=  -pot * r2_inv * rijz;
                potf   +=  pot;
#endif
            }
            force[i].acc -= PS::F64vec(afx, afy, afz);
            force[i].phi -= potf;
        }
    }
#endif
};



