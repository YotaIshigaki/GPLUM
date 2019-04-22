#include <particle_simulator.hpp>
#include "avx.h"
#include "particle.h"

#ifndef PARALLEL_I2J8
    
template <class TParticleI, class TParticleJ, class TForce>
void CalcForceLongEP(const TParticleI * ep_i,
                     const PS::S32 n_ip,
                     const TParticleJ * ep_j,
                     const PS::S32 n_jp,
                     TForce * force)
{
#ifdef __AVX512F__
#ifdef PARALLEL_I4J4
    const PS::S32 n_iparallel = 4;
    const PS::S32 n_jparallel = 4;
    //const PS::S32 N = 16;
    typedef v16sf VEC0;
    typedef v4sf  VEC1;
    typedef PS::F32  REAL;
#else
    const PS::S32 n_iparallel = 8;
    const PS::S32 n_jparallel = 2;
    //const PS::S32 N = 32;
    typedef v16sf VEC0;
    typedef v8sf  VEC1;
    typedef PS::F32  REAL;
#endif
#else //__AVX512F__
#ifndef CALC_EP_64bit
    const PS::S32 n_iparallel = 4;
    const PS::S32 n_jparallel = 2;
    //const PS::S32 N = 16;
    typedef v8sf VEC0;
    typedef v4sf VEC1;
    typedef PS::F32  REAL;
#else
    const PS::S32 n_iparallel = 2;
    const PS::S32 n_jparallel = 2;
    //const PS::S32 N = 16;
    typedef v4df VEC0;
    typedef v2df VEC1;
    typedef PS::F64  REAL;
#endif
#endif //__AVX512F__
    
    VEC0 (*rsqrt)(VEC0) = VEC0::rsqrt_1st;
    const PS::S32 nvector = VEC0::getVectorLength();
    assert ( n_iparallel * n_jparallel == nvector );

    const VEC0 eps2((REAL)ep_i[0].getEps2());
#ifndef CUTOFF_0
    const VEC0 g2((REAL)ep_i[0].getGamma2());
    const VEC0 g2_1_inv((REAL)ep_i[0].getGamma2_1_inv());

    const VEC0 v1(1.);
#endif
    const VEC1 v0(0.);

    //#pragma omp parallel for 
    for (PS::S32 i = 0; i < n_ip; i += n_iparallel) {
#if !defined(__AVX512F__) || defined(PARALLEL_I4J4)
        REAL buf_px[n_iparallel] __attribute__((aligned(16)));
        REAL buf_py[n_iparallel] __attribute__((aligned(16)));
        REAL buf_pz[n_iparallel] __attribute__((aligned(16)));
        REAL buf_rc[n_iparallel] __attribute__((aligned(16)));
#else
        REAL buf_px[n_iparallel] __attribute__((aligned(32)));
        REAL buf_py[n_iparallel] __attribute__((aligned(32)));
        REAL buf_pz[n_iparallel] __attribute__((aligned(32)));
        REAL buf_rc[n_iparallel] __attribute__((aligned(32)));
#endif
        const PS::S32 nii = std::min(n_ip - i, n_iparallel);
        for (PS::S32 ii = 0; ii < nii; ii++) {
            buf_px[ii] = (REAL)ep_i[i+ii].pos[0];
            buf_py[ii] = (REAL)ep_i[i+ii].pos[1];
            buf_pz[ii] = (REAL)ep_i[i+ii].pos[2];
#ifdef CUTOFF_0
            buf_rc[ii] = (REAL)(ep_i[i+ii].getROut());
#else
            buf_rc[ii] = (REAL)(ep_i[i+ii].getROut_inv());
#endif
        }
        VEC1 px_i_1;
        VEC1 py_i_1;
        VEC1 pz_i_1;
        VEC1 rc_i_1;
        px_i_1.load(buf_px);
        py_i_1.load(buf_py);
        pz_i_1.load(buf_pz);
        rc_i_1.load(buf_rc);
        VEC0 px_i(px_i_1);
        VEC0 py_i(py_i_1);
        VEC0 pz_i(pz_i_1);
        VEC0 rc_i(rc_i_1);
        VEC1 ax_i(0.0);
        VEC1 ay_i(0.0);
        VEC1 az_i(0.0);
        VEC1 pt_i(0.0);

        VEC1 px_j_1[n_jparallel];
        VEC1 py_j_1[n_jparallel];
        VEC1 pz_j_1[n_jparallel];
        VEC1 rc_j_1[n_jparallel];
        VEC1 ms_j_1[n_jparallel];
        for (PS::S32 j = 0; j < n_jp; j += n_jparallel) { 
            const PS::S32 njj = std::min(n_jp - j, n_jparallel);
            for (PS::S32 jj = 0; jj < njj; jj++) {
                px_j_1[jj] = VEC1((REAL)ep_j[j+jj].pos[0]);
                py_j_1[jj] = VEC1((REAL)ep_j[j+jj].pos[1]);
                pz_j_1[jj] = VEC1((REAL)ep_j[j+jj].pos[2]);
#ifdef CUTOFF_0
                rc_j_1[jj] = VEC1((REAL)ep_j[j+jj].getROut());
#else
                rc_j_1[jj] = VEC1((REAL)ep_j[j+jj].getROut_inv());
#endif
                ms_j_1[jj] = VEC1((REAL)ep_j[j+jj].getCharge());
            }
            for (PS::S32 jj = njj; jj < n_jparallel; jj++)
                ms_j_1[jj] = v0;
#ifdef PARALLEL_I4J4
            VEC0 dpx_ij = px_i - VEC0(px_j_1[0], px_j_1[1], px_j_1[2], px_j_1[3]);
            VEC0 dpy_ij = py_i - VEC0(py_j_1[0], py_j_1[1], py_j_1[2], py_j_1[3]);
            VEC0 dpz_ij = pz_i - VEC0(pz_j_1[0], pz_j_1[1], pz_j_1[2], pz_j_1[3]);
#ifdef CUTOFF_0
            VEC0 rc_ij  = VEC0::max(rc_i, VEC0(rc_j_1[0], rc_j_1[1], rc_j_1[2], rc_j_1[3]));
#else
            VEC0 rc_ij  = VEC0::min(rc_i, VEC0(rc_j_1[0], rc_j_1[1], rc_j_1[2], rc_j_1[3]));
#endif
            VEC0 ms_j   = VEC0(ms_j_1[0], ms_j_1[1], ms_j_1[2], ms_j_1[3]);
#else
            VEC0 dpx_ij = px_i - VEC0(px_j_1[0], px_j_1[1]);
            VEC0 dpy_ij = py_i - VEC0(py_j_1[0], py_j_1[1]);
            VEC0 dpz_ij = pz_i - VEC0(pz_j_1[0], pz_j_1[1]);
#ifdef CUTOFF_0
            VEC0 rc_ij  = VEC0::max(rc_i, VEC0(rc_j_1[0], rc_j_1[1]));
#else
            VEC0 rc_ij  = VEC0::min(rc_i, VEC0(rc_j_1[0], rc_j_1[1]));
#endif
            VEC0 ms_j   = VEC0(ms_j_1[0], ms_j_1[1]);
#endif       

            VEC0 dr2 = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij + eps2;
#ifdef CUTOFF_0
            VEC0 dr_inv  = rsqrt(VEC0::max(dr2, rc_ij*rc_ij));
            VEC0 dr3_inv = dr_inv * dr_inv * dr_inv;
#else
            VEC0 dr_rc_2 = VEC0::max( dr2 * rc_ij*rc_ij, g2 );
            VEC0 dr_inv0 = rsqrt(dr_rc_2) * rc_ij;
            VEC0 dr_inv  = dr_inv0 * VEC0::min( (g2 - dr_rc_2) * g2_1_inv, v1 );
            VEC0 dr3_inv = dr_inv * dr_inv0 * dr_inv0;
#endif
            VEC0 dg2_ij = ms_j * dr3_inv;

#ifdef PARALLEL_I4J4
            pt_i -= VEC0::reduce16to4(ms_j   * dr_inv);
            ax_i -= VEC0::reduce16to4(dpx_ij * dg2_ij);
            ay_i -= VEC0::reduce16to4(dpy_ij * dg2_ij);
            az_i -= VEC0::reduce16to4(dpz_ij * dg2_ij);
#else
            pt_i -= VEC0::reduce(ms_j   * dr_inv);
            ax_i -= VEC0::reduce(dpx_ij * dg2_ij);
            ay_i -= VEC0::reduce(dpy_ij * dg2_ij);
            az_i -= VEC0::reduce(dpz_ij * dg2_ij);
#endif
        }

#if !defined(__AVX512F__) || defined(PARALLEL_I4J4)
        REAL buf_ax[n_iparallel] __attribute__((aligned(16)));
        REAL buf_ay[n_iparallel] __attribute__((aligned(16)));
        REAL buf_az[n_iparallel] __attribute__((aligned(16)));
        REAL buf_pt[n_iparallel] __attribute__((aligned(16)));
#else
        REAL buf_ax[n_iparallel] __attribute__((aligned(32)));
        REAL buf_ay[n_iparallel] __attribute__((aligned(32)));
        REAL buf_az[n_iparallel] __attribute__((aligned(32)));
        REAL buf_pt[n_iparallel] __attribute__((aligned(32)));
#endif
   
        ax_i.store(buf_ax);
        ay_i.store(buf_ay);
        az_i.store(buf_az);
        pt_i.store(buf_pt);

        for(PS::S32 ii = 0; ii < nii; ii++) {
            force[i+ii].acc[0] += buf_ax[ii];
            force[i+ii].acc[1] += buf_ay[ii];
            force[i+ii].acc[2] += buf_az[ii];
            force[i+ii].phi    += buf_pt[ii];
        }
    }
}
 
template <class TParticleI, class TParticleJ, class TForce>
void CalcForceLongSP(const TParticleI * ep_i,
                     const PS::S32 n_ip,
                     const TParticleJ * ep_j,
                     const PS::S32 n_jp,
                     TForce * force)
{
#ifdef __AVX512F__
#ifdef PARALLEL_I4J4
    const PS::S32 n_iparallel = 4;
    const PS::S32 n_jparallel = 4;
    //const PS::S32 N = 16;
    typedef v16sf VEC0;
    typedef v4sf  VEC1;
#else
    const PS::S32 n_iparallel = 8;
    const PS::S32 n_jparallel = 2;
    //const PS::S32 N = 32;
    typedef v16sf VEC0;
    typedef v8sf  VEC1;
#endif
#else //__AVX512F__
    const PS::S32 n_iparallel = 4;
    const PS::S32 n_jparallel = 2;
    //const PS::S32 N = 16;
    typedef v8sf VEC0;
    typedef v4sf VEC1;
#endif //__AVX512F__

    VEC0 (*rsqrt)(VEC0) = VEC0::rsqrt_1st;
    const PS::S32 nvector = VEC0::getVectorLength();
    const VEC0 eps2((PS::F32)ep_i[0].getEps2());
    assert ( n_iparallel * n_jparallel == nvector );
    const VEC1 v0(0.);

#ifdef USE_QUAD
    VEC0 v5(5.0);
    VEC0 mv2(-2.0);
    VEC0 v0p5(0.5);
    VEC0 v1p5(1.5);
#endif
    
    //#pragma omp parallel for 
    for(PS::S32 i = 0; i < n_ip; i += n_iparallel) {
#if !defined(__AVX512F__) || defined(PARALLEL_I4J4)
        PS::F32 buf_px[n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_py[n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_pz[n_iparallel] __attribute__((aligned(16)));
        //PS::F32 buf_e2[n_iparallel] __attribute__((aligned(16)));
#else
        PS::F32 buf_px[n_iparallel] __attribute__((aligned(32)));
        PS::F32 buf_py[n_iparallel] __attribute__((aligned(32)));
        PS::F32 buf_pz[n_iparallel] __attribute__((aligned(32)));
        //PS::F32 buf_e2[n_iparallel] __attribute__((aligned(32)));
#endif
        const PS::S32 nii = std::min(n_ip - i, n_iparallel);
        for(PS::S32 ii = 0; ii < nii; ii++) {
            buf_px[ii] = ep_i[i+ii].pos[0];
            buf_py[ii] = ep_i[i+ii].pos[1];
            buf_pz[ii] = ep_i[i+ii].pos[2];
            //buf_e2[ii] = eps2;
        }

        VEC1 px_i_1;
        VEC1 py_i_1;
        VEC1 pz_i_1;
        //VEC1 e2_i_1;
        px_i_1.load(buf_px);
        py_i_1.load(buf_py);
        pz_i_1.load(buf_pz);
        //e2_i_1.load(buf_e2);
        VEC0 px_i(px_i_1);
        VEC0 py_i(py_i_1);
        VEC0 pz_i(pz_i_1);
        //VEC0 e2_i(e2_i_1);
        VEC1 ax_i(0.0);
        VEC1 ay_i(0.0);
        VEC1 az_i(0.0);
        VEC1 pt_i(0.0);

        VEC1 px_j_1[n_jparallel];
        VEC1 py_j_1[n_jparallel];
        VEC1 pz_j_1[n_jparallel];
        VEC1 ms_j_1[n_jparallel];
#ifdef USE_QUAD
        VEC1 qj_xx_1[n_jparallel];
        VEC1 qj_xy_1[n_jparallel];
        VEC1 qj_xz_1[n_jparallel];
        VEC1 qj_yy_1[n_jparallel];
        VEC1 qj_yz_1[n_jparallel];
        VEC1 qj_zz_1[n_jparallel];
#endif
        for (PS::S32 j = 0; j < n_jp; j += n_jparallel) {
            const PS::S32 njj = std::min(n_jp - j, n_jparallel);
            for (PS::S32 jj = 0; jj < njj; jj++) {
                px_j_1[jj] = VEC1((PS::F32)ep_j[j+jj].pos[0]);
                py_j_1[jj] = VEC1((PS::F32)ep_j[j+jj].pos[1]);
                pz_j_1[jj] = VEC1((PS::F32)ep_j[j+jj].pos[2]);
                ms_j_1[jj] = VEC1((PS::F32)ep_j[j+jj].mass);
#ifdef USE_QUAD
                qj_xx_1[jj] = VEC1((PS::F32)ep_j[j+jj].quad.xx);
                qj_xy_1[jj] = VEC1((PS::F32)ep_j[j+jj].quad.xy);
                qj_xz_1[jj] = VEC1((PS::F32)ep_j[j+jj].quad.xz);
                qj_yy_1[jj] = VEC1((PS::F32)ep_j[j+jj].quad.yy);
                qj_yz_1[jj] = VEC1((PS::F32)ep_j[j+jj].quad.yz);
                qj_zz_1[jj] = VEC1((PS::F32)ep_j[j+jj].quad.zz);
#endif
            }
            for (PS::S32 jj = njj; jj < n_jparallel; jj++) {
                px_j_1[jj] = VEC1((PS::F32)ep_j[0].pos[0]);
                py_j_1[jj] = VEC1((PS::F32)ep_j[0].pos[1]);
                pz_j_1[jj] = VEC1((PS::F32)ep_j[0].pos[2]);
                ms_j_1[jj] = v0;
#ifdef USE_QUAD
                qj_xx_1[jj] = qj_xy_1[jj] = qj_xz_1[jj]
                    = qj_yy_1[jj] = qj_yz_1[jj] = qj_zz_1[jj] = v0;
#endif
            }
            
#ifdef PARALLEL_I4J4
            VEC0 dpx_ij = px_i - VEC0(px_j_1[0], px_j_1[1], px_j_1[2], px_j_1[3]);
            VEC0 dpy_ij = py_i - VEC0(py_j_1[0], py_j_1[1], py_j_1[2], py_j_1[3]);
            VEC0 dpz_ij = pz_i - VEC0(pz_j_1[0], pz_j_1[1], pz_j_1[2], pz_j_1[3]);
            VEC0 ms_j   = VEC0(ms_j_1[0], ms_j_1[1], ms_j_1[2], ms_j_1[3]);
#ifdef USE_QUAD
            VEC0 qj_xx = VEC0(qj_xx_1[0], qj_xx_1[1], qj_xx_1[2], qj_xx_1[3]);
            VEC0 qj_xy = VEC0(qj_xy_1[0], qj_xy_1[1], qj_xy_1[2], qj_xy_1[3]);
            VEC0 qj_xz = VEC0(qj_xz_1[0], qj_xz_1[1], qj_xz_1[2], qj_xz_1[3]);
            VEC0 qj_yy = VEC0(qj_yy_1[0], qj_yy_1[1], qj_yy_1[2], qj_yy_1[3]);
            VEC0 qj_yz = VEC0(qj_yz_1[0], qj_yz_1[1], qj_yz_1[2], qj_yz_1[3]);
            VEC0 qj_zz = VEC0(qj_zz_1[0], qj_zz_1[1], qj_zz_1[2], qj_zz_1[3]);
            VEC0 qj_tr = qj_xx + qj_yy + qj_zz;
#endif
#else //PARALLEL_I4J4
            VEC0 dpx_ij = px_i - VEC0(px_j_1[0], px_j_1[1]);
            VEC0 dpy_ij = py_i - VEC0(py_j_1[0], py_j_1[1]);
            VEC0 dpz_ij = pz_i - VEC0(pz_j_1[0], pz_j_1[1]);
            VEC0 ms_j   = VEC0(ms_j_1[0], ms_j_1[1]);
#ifdef USE_QUAD
            VEC0 qj_xx = VEC0(qj_xx_1[0], qj_xx_1[1]);
            VEC0 qj_xy = VEC0(qj_xy_1[0], qj_xy_1[1]);
            VEC0 qj_xz = VEC0(qj_xz_1[0], qj_xz_1[1]);
            VEC0 qj_yy = VEC0(qj_yy_1[0], qj_yy_1[1]);
            VEC0 qj_yz = VEC0(qj_yz_1[0], qj_yz_1[1]);
            VEC0 qj_zz = VEC0(qj_zz_1[0], qj_zz_1[1]);
            VEC0 qj_tr = qj_xx + qj_yy + qj_zz;
#endif
#endif //PARALLEL_I4J4

            VEC0 dr2 = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij + eps2;

            VEC0 dr_inv = rsqrt(dr2);
            VEC0 dr2_inv = dr_inv * dr_inv;
            VEC0 dr3_inv = dr2_inv * dr_inv;

#ifdef USE_QUAD
            VEC0 qr_x = qj_xx*dpx_ij + qj_xy*dpy_ij + qj_xz*dpz_ij;
            VEC0 qr_y = qj_xy*dpx_ij + qj_yy*dpy_ij + qj_yz*dpz_ij;
            VEC0 qr_z = qj_xz*dpx_ij + qj_yz*dpy_ij + qj_zz*dpz_ij;

            VEC0 rqr = dpx_ij*qr_x + dpy_ij*qr_y + dpz_ij*qr_z;

            VEC0 dr5_inv = dr3_inv * dr2_inv * v1p5;
            VEC0 rqr_r5  = dr5_inv * rqr;
            VEC0 rqr_r7  = dr2_inv * rqr_r5;            
            VEC0 A = ms_j*dr3_inv - qj_tr*dr5_inv + v5*rqr_r7;
            VEC0 B = mv2*dr5_inv;

#ifdef PARALLEL_I4J4 
            pt_i -= VEC0::reduce16to4(ms_j*dr_inv - v0p5*qj_tr*dr3_inv + rqr_r5);
            ax_i -= VEC0::reduce16to4(A*dpx_ij + B*qr_x);
            ay_i -= VEC0::reduce16to4(A*dpy_ij + B*qr_y);
            az_i -= VEC0::reduce16to4(A*dpz_ij + B*qr_z);
#else
            pt_i -= VEC0::reduce(ms_j*dr_inv - v0p5*qj_tr*dr3_inv + rqr_r5);
            ax_i -= VEC0::reduce(A*dpx_ij + B*qr_x);
            ay_i -= VEC0::reduce(A*dpy_ij + B*qr_y);
            az_i -= VEC0::reduce(A*dpz_ij + B*qr_z);
#endif
#else // USE_QUAD
#ifdef PARALLEL_I4J4
            pt_i -= VEC0::reduce16to4(ms_j    * dr_inv);
            ax_i -= VEC0::reduce16to4(dpx_ij * dg2_ij);
            ay_i -= VEC0::reduce16to4(dpy_ij * dg2_ij);
            az_i -= VEC0::reduce16to4(dpz_ij * dg2_ij);
#else
            pt_i -= VEC0::reduce(ms_j    * dr_inv);
            ax_i -= VEC0::reduce(dpx_ij * dg2_ij);
            ay_i -= VEC0::reduce(dpy_ij * dg2_ij);
            az_i -= VEC0::reduce(dpz_ij * dg2_ij);
#endif
#endif // USE_QUAD
        }

#if !defined(__AVX512F__) || defined(PARALLEL_I4J4)
        PS::F32 buf_ax[n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_ay[n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_az[n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_pt[n_iparallel] __attribute__((aligned(16)));
#else
        PS::F32 buf_ax[n_iparallel] __attribute__((aligned(32)));
        PS::F32 buf_ay[n_iparallel] __attribute__((aligned(32)));
        PS::F32 buf_az[n_iparallel] __attribute__((aligned(32)));
        PS::F32 buf_pt[n_iparallel] __attribute__((aligned(32)));
#endif

        ax_i.store(buf_ax);
        ay_i.store(buf_ay);
        az_i.store(buf_az);
        pt_i.store(buf_pt);

        for(PS::S32 ii = 0; ii < nii; ii++) {
            force[i+ii].acc[0] += buf_ax[ii];
            force[i+ii].acc[1] += buf_ay[ii];
            force[i+ii].acc[2] += buf_az[ii];
            force[i+ii].phi    += buf_pt[ii];
        }
    }
}

#else //PARALLEL_I2J8

#ifdef __AVX512F__

template <class TParticleI, class TParticleJ, class TForce>
void CalcForceLongEP(const TParticleI * ep_i,
                     const PS::S32 n_ip,
                     const TParticleJ * ep_j,
                     const PS::S32 n_jp,
                     TForce * force)
{
    const PS::S32 n_iparallel = 2;
    const PS::S32 n_jparallel = 8;
    
    v16sf (*rsqrt)(v16sf) = v16sf::rsqrt_1st;
    const PS::S32 nvector = v16sf::getVectorLength();
    assert ( n_iparallel * n_jparallel == nvector );

    const v16sf eps2((PS::F32)ep_i[0].getEps2());
#ifndef CUTOFF_0
    const v16sf g2((PS::F32)ep_i[0].getGamma2());
    const v16sf g2_1_inv((PS::F32)ep_i[0].getGamma2_1_inv());

    const v16sf v1(1.);
#endif

    v8sf px_i_v8[n_iparallel];
    v8sf py_i_v8[n_iparallel];
    v8sf pz_i_v8[n_iparallel];
    v8sf rc_i_v8[n_iparallel];
#pragma omp parallel for 
    for (PS::S32 i = 0; i < n_ip; i += n_iparallel) {
        const PS::S32 nii = std::min(n_ip - i, n_iparallel);
        for (PS::S32 ii = 0; ii < nii; ii++) {
            px_i_v8[ii] = v8sf((PS::F32)ep_i[i+ii].pos[0]);
            py_i_v8[ii] = v8sf((PS::F32)ep_i[i+ii].pos[1]);
            pz_i_v8[ii] = v8sf((PS::F32)ep_i[i+ii].pos[2]);
#ifdef CUTOFF_0
            rc_i_v8[ii] = v8sf((PS::F32)ep_i[i+ii].getROut());
#else
            rc_i_v8[ii] = v8sf((PS::F32)ep_i[i+ii].getROut_inv());
#endif
        }
        v16sf px_i = v16sf(px_i_v8[0], px_i_v8[1]);
        v16sf py_i = v16sf(py_i_v8[0], py_i_v8[1]);
        v16sf pz_i = v16sf(pz_i_v8[0], pz_i_v8[1]);
        v16sf rc_i = v16sf(rc_i_v8[0], rc_i_v8[1]);
        v4sf ax_i(0.0);
        v4sf ay_i(0.0);
        v4sf az_i(0.0);
        v4sf pt_i(0.0);

        
        PS::F32 buf_px[n_jparallel] __attribute__((aligned(32)));
        PS::F32 buf_py[n_jparallel] __attribute__((aligned(32)));
        PS::F32 buf_pz[n_jparallel] __attribute__((aligned(32)));
        PS::F32 buf_rc[n_jparallel] __attribute__((aligned(32)));
        PS::F32 buf_ms[n_jparallel] __attribute__((aligned(32)));
        for (PS::S32 j = 0; j < n_jp; j += n_jparallel) {
            const PS::S32 njj = std::min(n_jp - j, n_jparallel);
            for (PS::S32 jj = 0; jj < njj; jj++) {
                buf_px[jj] = (PS::F32)ep_j[j+jj].pos[0];
                buf_py[jj] = (PS::F32)ep_j[j+jj].pos[1];
                buf_pz[jj] = (PS::F32)ep_j[j+jj].pos[2];
#ifdef CUTOFF_0
                buf_rc[jj] = (PS::F32)ep_j[j+jj].getROut();
#else
                buf_rc[jj] = (PS::F32)ep_j[j+jj].getROut_inv();
#endif
                buf_ms[jj] = (PS::F32)ep_j[j+jj].getCharge();
            }
            for (PS::S32 jj = njj; jj < n_jparallel; jj++)
                buf_ms[jj] = 0.;
            v8sf px_j_v8;
            v8sf py_j_v8;
            v8sf pz_j_v8;
            v8sf rc_j_v8;
            v8sf ms_j_v8;
            px_j_v8.load(buf_px);
            py_j_v8.load(buf_py);
            pz_j_v8.load(buf_pz);
            rc_j_v8.load(buf_rc);
            ms_j_v8.load(buf_ms);
            v16sf dpx_ij = px_i - v16sf(px_j_v8);
            v16sf dpy_ij = py_i - v16sf(py_j_v8);
            v16sf dpz_ij = pz_i - v16sf(pz_j_v8);
#ifdef CUTOFF_0
            v16sf rc_ij = v16sf::max(rc_i, v16sf(rc_j_v8));
#else
            v16sf rc_ij = v16sf::min(rc_i, v16sf(rc_j_v8));
#endif
            v16sf ms_j  = v16sf(ms_j_v8);
            
            v16sf dr2 = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij + eps2;
#ifdef CUTOFF_0
            v16sf dr_inv  = rsqrt(v16sf::max(dr2, rc_ij*rc_ij));
            v16sf dr3_inv = dr_inv * dr_inv * dr_inv;
#else
            v16sf dr_rc_2 = v16sf::max( dr2 * rc_ij*rc_ij, g2 );
            v16sf dr_inv0 = rsqrt(dr_rc_2) * rc_ij;
            v16sf dr_inv  = dr_inv0 * v16sf::min( (g2 - dr_rc_2) * g2_1_inv, v1 );
            v16sf dr3_inv = dr_inv * dr_inv0 * dr_inv0;
#endif
            v16sf dg2_ij = ms_j * dr3_inv;
            
            pt_i -= v16sf::reduce16to4(ms_j  * dr_inv);
            ax_i -= v16sf::reduce16to4(dpx_ij * dg2_ij);
            ay_i -= v16sf::reduce16to4(dpy_ij * dg2_ij);
            az_i -= v16sf::reduce16to4(dpz_ij * dg2_ij);
        }

        PS::F32 buf_ax[2*n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_ay[2*n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_az[2*n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_pt[2*n_iparallel] __attribute__((aligned(16)));
   
        ax_i.store(buf_ax);
        ay_i.store(buf_ay);
        az_i.store(buf_az);
        pt_i.store(buf_pt);

        for(PS::S32 ii = 0; ii < nii; ii++) {
            force[i+ii].acc[0] += buf_ax[2*ii] + buf_ax[2*ii+1];
            force[i+ii].acc[1] += buf_ay[2*ii] + buf_ay[2*ii+1];
            force[i+ii].acc[2] += buf_az[2*ii] + buf_az[2*ii+1];
            force[i+ii].phi    += buf_pt[2*ii] + buf_pt[2*ii+1];
        }
    }
}
 
template <class TParticleI, class TParticleJ, class TForce>
void CalcForceLongSP(const TParticleI * ep_i,
                     const PS::S32 n_ip,
                     const TParticleJ * ep_j,
                     const PS::S32 n_jp,
                     TForce * force)
{
    const PS::S32 n_iparallel = 2;
    const PS::S32 n_jparallel = 8;

    PS::F32 eps2 = ep_i[0].getEps2();
    v16sf (*rsqrt)(v16sf) = v16sf::rsqrt_1st;
    PS::S32 nvector = v16sf::getVectorLength();
    assert ( n_iparallel * n_jparallel == nvector );

#ifdef USE_QUAD
    v16sf v5(5.0);
    v16sf mv2(-2.0);
    v16sf v0p5(0.5);
    v16sf v1p5(1.5);
#endif

    v8sf px_i_v8[n_iparallel];
    v8sf py_i_v8[n_iparallel];
    v8sf pz_i_v8[n_iparallel];
    v8sf e2_i_v8[n_iparallel];
#pragma omp parallel for 
    for (PS::S32 i = 0; i < n_ip; i += n_iparallel) {
        const PS::S32 nii = std::min(n_ip - i, n_iparallel);
        for (PS::S32 ii = 0; ii < nii; ii++) {
            px_i_v8[ii] = v8sf((PS::F32)ep_i[i+ii].pos[0]);
            py_i_v8[ii] = v8sf((PS::F32)ep_i[i+ii].pos[1]);
            pz_i_v8[ii] = v8sf((PS::F32)ep_i[i+ii].pos[2]);
            e2_i_v8[ii] = v8sf(eps2);
        }
        v16sf px_i = v16sf(px_i_v8[0], px_i_v8[1]);
        v16sf py_i = v16sf(py_i_v8[0], py_i_v8[1]);
        v16sf pz_i = v16sf(pz_i_v8[0], pz_i_v8[1]);
        v16sf e2_i = v16sf(e2_i_v8[0], e2_i_v8[1]);
        v4sf ax_i(0.0);
        v4sf ay_i(0.0);
        v4sf az_i(0.0);
        v4sf pt_i(0.0);

        PS::F32 buf_px[n_jparallel] __attribute__((aligned(32)));
        PS::F32 buf_py[n_jparallel] __attribute__((aligned(32)));
        PS::F32 buf_pz[n_jparallel] __attribute__((aligned(32)));
        PS::F32 buf_ms[n_jparallel] __attribute__((aligned(32)));
#ifdef USE_QUAD
        PS::F32 buf_qj_xx[n_jparallel] __attribute__((aligned(32)));
        PS::F32 buf_qj_xy[n_jparallel] __attribute__((aligned(32)));
        PS::F32 buf_qj_xz[n_jparallel] __attribute__((aligned(32)));
        PS::F32 buf_qj_yy[n_jparallel] __attribute__((aligned(32)));
        PS::F32 buf_qj_yz[n_jparallel] __attribute__((aligned(32)));
        PS::F32 buf_qj_zz[n_jparallel] __attribute__((aligned(32)));
#endif
        for (PS::S32 j = 0; j < n_jp; j += n_jparallel) {
            const PS::S32 njj = std::min(n_jp - j, n_jparallel);
            for (PS::S32 jj = 0; jj < njj; jj++) {
                buf_px[jj] = (PS::F32)ep_j[j+jj].pos[0];
                buf_py[jj] = (PS::F32)ep_j[j+jj].pos[1];
                buf_pz[jj] = (PS::F32)ep_j[j+jj].pos[2];
                buf_ms[jj] = (PS::F32)ep_j[j+jj].mass;
#ifdef USE_QUAD
                buf_qj_xx[jj] = (PS::F32)ep_j[j+jj].quad.xx;
                buf_qj_xy[jj] = (PS::F32)ep_j[j+jj].quad.xy;
                buf_qj_xz[jj] = (PS::F32)ep_j[j+jj].quad.xz;
                buf_qj_yy[jj] = (PS::F32)ep_j[j+jj].quad.yy;
                buf_qj_yz[jj] = (PS::F32)ep_j[j+jj].quad.yz;
                buf_qj_zz[jj] = (PS::F32)ep_j[j+jj].quad.zz;
#endif
            }
            for (PS::S32 jj = njj; jj < n_jparallel; jj++){
                buf_ms[jj] = 0.;
#ifdef USE_QUAD
                buf_qj_xx[jj] = 0.;
                buf_qj_xy[jj] = 0.;
                buf_qj_xz[jj] = 0.;
                buf_qj_yy[jj] = 0.;
                buf_qj_yz[jj] = 0.;
                buf_qj_zz[jj] = 0.;
#endif
            }
            v8sf px_j_v8;
            v8sf py_j_v8;
            v8sf pz_j_v8;
            v8sf ms_j_v8;
            px_j_v8.load(buf_px);
            py_j_v8.load(buf_py);
            pz_j_v8.load(buf_pz);
            ms_j_v8.load(buf_ms);
#ifdef USE_QUAD
            v8sf qj_xx_v8;
            v8sf qj_xy_v8;
            v8sf qj_xz_v8;
            v8sf qj_yy_v8;
            v8sf qj_yz_v8;
            v8sf qj_zz_v8;
            qj_xx_v8.load(buf_qj_xx);
            qj_xy_v8.load(buf_qj_xy);
            qj_xz_v8.load(buf_qj_xz);
            qj_yy_v8.load(buf_qj_yy);
            qj_yz_v8.load(buf_qj_yz);
            qj_zz_v8.load(buf_qj_zz);
#endif
            v16sf dpx_ij = px_i - v16sf(px_j_v8);
            v16sf dpy_ij = py_i - v16sf(py_j_v8);
            v16sf dpz_ij = pz_i - v16sf(pz_j_v8);
            v16sf ms_j   = v16sf(ms_j_v8);
#ifdef USE_QUAD
            v16sf qj_xx  = v16sf(qj_xx_v8);
            v16sf qj_xy  = v16sf(qj_xy_v8);
            v16sf qj_xz  = v16sf(qj_xz_v8);
            v16sf qj_yy  = v16sf(qj_yy_v8);
            v16sf qj_yz  = v16sf(qj_yz_v8);
            v16sf qj_zz  = v16sf(qj_zz_v8);
            v16sf qj_tr  = qj_xx + qj_yy + qj_zz;
#endif

            v16sf dr2 = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij + e2_i;

            v16sf dr_inv = rsqrt(dr2);
            v16sf dr2_inv = dr_inv * dr_inv;
            v16sf dr3_inv = dr2_inv * dr_inv;

#ifdef USE_QUAD
            v16sf qr_x = qj_xx*dpx_ij + qj_xy*dpy_ij + qj_xz*dpz_ij;
            v16sf qr_y = qj_xy*dpx_ij + qj_yy*dpy_ij + qj_yz*dpz_ij;
            v16sf qr_z = qj_xz*dpx_ij + qj_yz*dpy_ij + qj_zz*dpz_ij;

            v16sf rqr = dpx_ij*qr_x + dpy_ij*qr_y + dpz_ij*qr_z;

            v16sf dr5_inv = dr3_inv * dr2_inv * v1p5;
            v16sf rqr_r5  = dr5_inv * rqr;
            v16sf rqr_r7  = dr2_inv * rqr_r5;            
            v16sf A = ms_j*dr3_inv - qj_tr*dr5_inv + v5*rqr_r7;
            v16sf B = mv2*dr5_inv;

            pt_i -= v16sf::reduce16to4(ms_j*dr_inv - v0p5*qj_tr*dr3_inv + rqr_r5);
            ax_i -= v16sf::reduce16to4(A*dpx_ij + B*qr_x);
            ay_i -= v16sf::reduce16to4(A*dpy_ij + B*qr_y);
            az_i -= v16sf::reduce16to4(A*dpz_ij + B*qr_z);
#else
            pt_i -= v16sf::reduce16to4(ms_j    * dr_inv);
            ax_i -= v16sf::reduce16to4(dpx_ij * dg2_ij);
            ay_i -= v16sf::reduce16to4(dpy_ij * dg2_ij);
            az_i -= v16sf::reduce16to4(dpz_ij * dg2_ij);
#endif
        }

        PS::F32 buf_ax[2*n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_ay[2*n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_az[2*n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_pt[2*n_iparallel] __attribute__((aligned(16)));

        ax_i.store(buf_ax);
        ay_i.store(buf_ay);
        az_i.store(buf_az);
        pt_i.store(buf_pt);

        for(PS::S32 ii = 0; ii < nii; ii++) {
            force[i+ii].acc[0] += buf_ax[2*ii] + buf_ax[2*ii+1];
            force[i+ii].acc[1] += buf_ay[2*ii] + buf_ay[2*ii+1];
            force[i+ii].acc[2] += buf_az[2*ii] + buf_az[2*ii+1];
            force[i+ii].phi    += buf_pt[2*ii] + buf_pt[2*ii+1];
        }
    }
}

#else //__AVX512F__

template <class TParticleI, class TParticleJ, class TForce>
void CalcForceLongEP(const TParticleI * ep_i,
                     const PS::S32 n_ip,
                     const TParticleJ * ep_j,
                     const PS::S32 n_jp,
                     TForce * force)
{
    const PS::S32 n_iparallel = 2;
    const PS::S32 n_jparallel = 4;
    
    v8sf (*rsqrt)(v8sf) = v8sf::rsqrt_1st;
    const PS::S32 nvector = v8sf::getVectorLength();
    assert ( n_iparallel * n_jparallel == nvector );

    const v8sf eps2((PS::F32)ep_i[0].getEps2());
#ifndef CUTOFF_0
    const v8sf g2((PS::F32)ep_i[0].getGamma2());
    const v8sf g2_1_inv((PS::F32)ep_i[0].getGamma2_1_inv());

    const v8sf v1(1.);
#endif

    v4sf px_i_v4[n_iparallel];
    v4sf py_i_v4[n_iparallel];
    v4sf pz_i_v4[n_iparallel];
    v4sf rc_i_v4[n_iparallel];
#pragma omp parallel for 
    for (PS::S32 i = 0; i < n_ip; i += n_iparallel) {
        const PS::S32 nii = std::min(n_ip - i, n_iparallel);
        for (PS::S32 ii = 0; ii < nii; ii++) {
            px_i_v4[ii] = v4sf((PS::F32)ep_i[i+ii].pos[0]);
            py_i_v4[ii] = v4sf((PS::F32)ep_i[i+ii].pos[1]);
            pz_i_v4[ii] = v4sf((PS::F32)ep_i[i+ii].pos[2]);
#ifdef CUTOFF_0
            rc_i_v4[ii] = v4sf((PS::F32)ep_i[i+ii].getROut());
#else
            rc_i_v4[ii] = v4sf((PS::F32)ep_i[i+ii].getROut_inv());
#endif
        }
        v8sf px_i = v8sf(px_i_v4[0], px_i_v4[1]);
        v8sf py_i = v8sf(py_i_v4[0], py_i_v4[1]);
        v8sf pz_i = v8sf(pz_i_v4[0], pz_i_v4[1]);
        v8sf rc_i = v8sf(rc_i_v4[0], rc_i_v4[1]);
        v4sf ax_i(0.0);
        v4sf ay_i(0.0);
        v4sf az_i(0.0);
        v4sf pt_i(0.0);

        
        PS::F32 buf_px[n_jparallel] __attribute__((aligned(16)));
        PS::F32 buf_py[n_jparallel] __attribute__((aligned(16)));
        PS::F32 buf_pz[n_jparallel] __attribute__((aligned(16)));
        PS::F32 buf_rc[n_jparallel] __attribute__((aligned(16)));
        PS::F32 buf_ms[n_jparallel] __attribute__((aligned(16)));
        for (PS::S32 j = 0; j < n_jp; j += n_jparallel) {
            const PS::S32 njj = std::min(n_jp - j, n_jparallel);
            for (PS::S32 jj = 0; jj < njj; jj++) {
                buf_px[jj] = (PS::F32)ep_j[j+jj].pos[0];
                buf_py[jj] = (PS::F32)ep_j[j+jj].pos[1];
                buf_pz[jj] = (PS::F32)ep_j[j+jj].pos[2];
#ifdef CUTOFF_0
                buf_rc[jj] = (PS::F32)ep_j[j+jj].getROut();
#else
                buf_rc[jj] = (PS::F32)ep_j[j+jj].getROut_inv();
#endif
                buf_ms[jj] = (PS::F32)ep_j[j+jj].getCharge();
            }
            for (PS::S32 jj = njj; jj < n_jparallel; jj++)
                buf_ms[jj] = 0.;
            v4sf px_j_v4;
            v4sf py_j_v4;
            v4sf pz_j_v4;
            v4sf rc_j_v4;
            v4sf ms_j_v4;
            px_j_v4.load(buf_px);
            py_j_v4.load(buf_py);
            pz_j_v4.load(buf_pz);
            rc_j_v4.load(buf_rc);
            ms_j_v4.load(buf_ms);
            v8sf dpx_ij = px_i - v8sf(px_j_v4);
            v8sf dpy_ij = py_i - v8sf(py_j_v4);
            v8sf dpz_ij = pz_i - v8sf(pz_j_v4);
#ifdef CUTOFF_0
            v8sf rc_ij = v8sf::max(rc_i, v8sf(rc_j_v4));
#else
            v8sf rc_ij = v8sf::min(rc_i, v8sf(rc_j_v4));
#endif
            v8sf ms_j  = v8sf(ms_j_v4);
            
            v8sf dr2 = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij + eps2;
#ifdef CUTOFF_0
            v8sf dr_inv  = rsqrt(v8sf::max(dr2, rc_ij*rc_ij));
            v8sf dr3_inv = dr_inv * dr_inv * dr_inv;
#else
            v8sf dr_rc_2 = v8sf::max( dr2 * rc_ij*rc_ij, g2 );
            v8sf dr_inv0 = rsqrt(dr_rc_2) * rc_ij;
            v8sf dr_inv  = dr_inv0 * v8sf::min( (g2 - dr_rc_2) * g2_1_inv, v1 );
            v8sf dr3_inv = dr_inv * dr_inv0 * dr_inv0;
#endif
            v8sf dg2_ij = ms_j * dr3_inv;
            
            pt_i -= v8sf::reduce(ms_j  * dr_inv);
            ax_i -= v8sf::reduce(dpx_ij * dg2_ij);
            ay_i -= v8sf::reduce(dpy_ij * dg2_ij);
            az_i -= v8sf::reduce(dpz_ij * dg2_ij);
        }

        PS::F32 buf_ax[2*n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_ay[2*n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_az[2*n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_pt[2*n_iparallel] __attribute__((aligned(16)));
   
        ax_i.store(buf_ax);
        ay_i.store(buf_ay);
        az_i.store(buf_az);
        pt_i.store(buf_pt);

        for(PS::S32 ii = 0; ii < nii; ii++) {
            force[i+ii].acc[0] += buf_ax[2*ii] + buf_ax[2*ii+1];
            force[i+ii].acc[1] += buf_ay[2*ii] + buf_ay[2*ii+1];
            force[i+ii].acc[2] += buf_az[2*ii] + buf_az[2*ii+1];
            force[i+ii].phi    += buf_pt[2*ii] + buf_pt[2*ii+1];
        }
    }
}
 
template <class TParticleI, class TParticleJ, class TForce>TForce
void CalcForceLongSP(const TParticleI * ep_i,
                     const PS::S32 n_ip,
                     const TParticleJ * ep_j,
                     const PS::S32 n_jp,
                     TForce * force)
{
    const PS::S32 n_iparallel = 2;
    const PS::S32 n_jparallel = 4;
    
    PS::F32 eps2 = ep_i[0].getEps2();
    v8sf (*rsqrt)(v8sf) = v8sf::rsqrt_1st;
    PS::S32 nvector = v8sf::getVectorLength();
    assert ( n_iparallel * n_jparallel == nvector );

#ifdef USE_QUAD
    v8sf v5(5.0);
    v8sf mv2(-2.0);
    v8sf v0p5(0.5);
    v8sf v1p5(1.5);
#endif

    v4sf px_i_v4[n_iparallel];
    v4sf py_i_v4[n_iparallel];
    v4sf pz_i_v4[n_iparallel];
    v4sf e2_i_v4[n_iparallel];
#pragma omp parallel for 
    for (PS::S32 i = 0; i < n_ip; i += n_iparallel) {
        const PS::S32 nii = std::min(n_ip - i, n_iparallel);
        for (PS::S32 ii = 0; ii < nii; ii++) {
            px_i_v4[ii] = v4sf((PS::F32)ep_i[i+ii].pos[0]);
            py_i_v4[ii] = v4sf((PS::F32)ep_i[i+ii].pos[1]);
            pz_i_v4[ii] = v4sf((PS::F32)ep_i[i+ii].pos[2]);
            e2_i_v4[ii] = v4sf(eps2);
        }
        v8sf px_i = v8sf(px_i_v4[0], px_i_v4[1]);
        v8sf py_i = v8sf(py_i_v4[0], py_i_v4[1]);
        v8sf pz_i = v8sf(pz_i_v4[0], pz_i_v4[1]);
        v8sf e2_i = v8sf(e2_i_v4[0], e2_i_v4[1]);
        v4sf ax_i(0.0);
        v4sf ay_i(0.0);
        v4sf az_i(0.0);
        v4sf pt_i(0.0);

        PS::F32 buf_px[n_jparallel] __attribute__((aligned(16)));
        PS::F32 buf_py[n_jparallel] __attribute__((aligned(16)));
        PS::F32 buf_pz[n_jparallel] __attribute__((aligned(16)));
        PS::F32 buf_ms[n_jparallel] __attribute__((aligned(16)));
#ifdef USE_QUAD
        PS::F32 buf_qj_xx[n_jparallel] __attribute__((aligned(16)));
        PS::F32 buf_qj_xy[n_jparallel] __attribute__((aligned(16)));
        PS::F32 buf_qj_xz[n_jparallel] __attribute__((aligned(16)));
        PS::F32 buf_qj_yy[n_jparallel] __attribute__((aligned(16)));
        PS::F32 buf_qj_yz[n_jparallel] __attribute__((aligned(16)));
        PS::F32 buf_qj_zz[n_jparallel] __attribute__((aligned(16)));
#endif
        for (PS::S32 j = 0; j < n_jp; j += n_jparallel) {
            const PS::S32 njj = std::min(n_jp - j, n_jparallel);
            for (PS::S32 jj = 0; jj < njj; jj++) {
                buf_px[jj] = (PS::F32)ep_j[j+jj].pos[0];
                buf_py[jj] = (PS::F32)ep_j[j+jj].pos[1];
                buf_pz[jj] = (PS::F32)ep_j[j+jj].pos[2];
                buf_ms[jj] = (PS::F32)ep_j[j+jj].mass;
#ifdef USE_QUAD
                buf_qj_xx[jj] = (PS::F32)ep_j[j+jj].quad.xx;
                buf_qj_xy[jj] = (PS::F32)ep_j[j+jj].quad.xy;
                buf_qj_xz[jj] = (PS::F32)ep_j[j+jj].quad.xz;
                buf_qj_yy[jj] = (PS::F32)ep_j[j+jj].quad.yy;
                buf_qj_yz[jj] = (PS::F32)ep_j[j+jj].quad.yz;
                buf_qj_zz[jj] = (PS::F32)ep_j[j+jj].quad.zz;
#endif
            }
            for (PS::S32 jj = njj; jj < n_jparallel; jj++){
                buf_ms[jj] = 0.;
#ifdef USE_QUAD
                buf_qj_xx[jj] = 0.;
                buf_qj_xy[jj] = 0.;
                buf_qj_xz[jj] = 0.;
                buf_qj_yy[jj] = 0.;
                buf_qj_yz[jj] = 0.;
                buf_qj_zz[jj] = 0.;
#endif
            }
            v4sf px_j_v4;
            v4sf py_j_v4;
            v4sf pz_j_v4;
            v4sf ms_j_v4;
            px_j_v4.load(buf_px);
            py_j_v4.load(buf_py);
            pz_j_v4.load(buf_pz);
            ms_j_v4.load(buf_ms);
#ifdef USE_QUAD
            v4sf qj_xx_v4;
            v4sf qj_xy_v4;
            v4sf qj_xz_v4;
            v4sf qj_yy_v4;
            v4sf qj_yz_v4;
            v4sf qj_zz_v4;
            qj_xx_v4.load(buf_qj_xx);
            qj_xy_v4.load(buf_qj_xy);
            qj_xz_v4.load(buf_qj_xz);
            qj_yy_v4.load(buf_qj_yy);
            qj_yz_v4.load(buf_qj_yz);
            qj_zz_v4.load(buf_qj_zz);
#endif
            v8sf dpx_ij = px_i - v8sf(px_j_v4);
            v8sf dpy_ij = py_i - v8sf(py_j_v4);
            v8sf dpz_ij = pz_i - v8sf(pz_j_v4);
            v8sf ms_j   = v8sf(ms_j_v4);
#ifdef USE_QUAD
            v8sf qj_xx  = v8sf(qj_xx_v4);
            v8sf qj_xy  = v8sf(qj_xy_v4);
            v8sf qj_xz  = v8sf(qj_xz_v4);
            v8sf qj_yy  = v8sf(qj_yy_v4);
            v8sf qj_yz  = v8sf(qj_yz_v4);
            v8sf qj_zz  = v8sf(qj_zz_v4);
            v8sf qj_tr  = qj_xx + qj_yy + qj_zz;
#endif

            v8sf dr2 = dpx_ij * dpx_ij + dpy_ij * dpy_ij + dpz_ij * dpz_ij + e2_i;

            v8sf dr_inv = rsqrt(dr2);
            v8sf dr2_inv = dr_inv * dr_inv;
            v8sf dr3_inv = dr2_inv * dr_inv;

#ifdef USE_QUAD
            v8sf qr_x = qj_xx*dpx_ij + qj_xy*dpy_ij + qj_xz*dpz_ij;
            v8sf qr_y = qj_xy*dpx_ij + qj_yy*dpy_ij + qj_yz*dpz_ij;
            v8sf qr_z = qj_xz*dpx_ij + qj_yz*dpy_ij + qj_zz*dpz_ij;

            v8sf rqr = dpx_ij*qr_x + dpy_ij*qr_y + dpz_ij*qr_z;

            v8sf dr5_inv = dr3_inv * dr2_inv * v1p5;
            v8sf rqr_r5  = dr5_inv * rqr;
            v8sf rqr_r7  = dr2_inv * rqr_r5;            
            v8sf A = ms_j*dr3_inv - qj_tr*dr5_inv + v5*rqr_r7;
            v8sf B = mv2*dr5_inv;

            pt_i -= v8sf::reduce(ms_j*dr_inv - v0p5*qj_tr*dr3_inv + rqr_r5);
            ax_i -= v8sf::reduce(A*dpx_ij + B*qr_x);
            ay_i -= v8sf::reduce(A*dpy_ij + B*qr_y);
            az_i -= v8sf::reduce(A*dpz_ij + B*qr_z);
#else
            pt_i -= v8sf::reduce(ms_j    * dr_inv);
            ax_i -= v8sf::reduce(dpx_ij * dg2_ij);
            ay_i -= v8sf::reduce(dpy_ij * dg2_ij);
            az_i -= v8sf::reduce(dpz_ij * dg2_ij);
#endif
        }

        PS::F32 buf_ax[2*n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_ay[2*n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_az[2*n_iparallel] __attribute__((aligned(16)));
        PS::F32 buf_pt[2*n_iparallel] __attribute__((aligned(16)));

        ax_i.store(buf_ax);
        ay_i.store(buf_ay);
        az_i.store(buf_az);
        pt_i.store(buf_pt);

        for(PS::S32 ii = 0; ii < nii; ii++) {
            force[i+ii].acc[0] += buf_ax[2*ii] + buf_ax[2*ii+1];
            force[i+ii].acc[1] += buf_ay[2*ii] + buf_ay[2*ii+1];
            force[i+ii].acc[2] += buf_az[2*ii] + buf_az[2*ii+1];
            force[i+ii].phi    += buf_pt[2*ii] + buf_pt[2*ii+1];
        }
    }
}

#endif //__AVX512F__

#endif //PARALLEL_I2J8



template 
void CalcForceLongEP<EPGrav, EPGrav, ForceGrav>(const EPGrav * ep_i,
                                                const PS::S32 n_ip,
                                                const EPGrav * ep_j,
                                                const PS::S32 n_jp,
                                                ForceGrav * force);

#ifdef USE_INDIVIDUAL_RADII
#ifdef USE_QUAD
template 
void CalcForceLongSP<EPGrav, PS::SPJQuadrupoleInAndOut, ForceGrav>(const EPGrav * ep_i,
                                                                   const PS::S32 n_ip,
                                                                   const PS::SPJQuadrupoleInAndOut * ep_j,
                                                                   const PS::S32 n_jp,
                                                                   ForceGrav * force);
#else
template 
void CalcForceLongSP<EPGrav, PS::SPJMonopoleInAndOut, ForceGrav>(const EPGrav * ep_i,
                                                                   const PS::S32 n_ip,
                                                                   const PS::SPJMonopoleInAndOut * ep_j,
                                                                   const PS::S32 n_jp,
                                                                   ForceGrav * force);
#endif //USE_QUAD
#else
#ifdef USE_QUAD
template 
void CalcForceLongSP<EPGrav, PS::SPJQuadrupoleScatter, ForceGrav>(const EPGrav * ep_i,
                                                                   const PS::S32 n_ip,
                                                                   const PS::SPJQuadrupoleScatter * ep_j,
                                                                   const PS::S32 n_jp,
                                                                   ForceGrav * force);
#else
template 
void CalcForceLongSP<EPGrav, PS::SPJMonopoleScatter, ForceGrav>(const EPGrav * ep_i,
                                                                   const PS::S32 n_ip,
                                                                   const PS::SPJMonopoleScatter * ep_j,
                                                                   const PS::S32 n_jp,
                                                                   ForceGrav * force);
#endif //USE_QUAD
#endif //USE_INDIVIDUAL_RADII
