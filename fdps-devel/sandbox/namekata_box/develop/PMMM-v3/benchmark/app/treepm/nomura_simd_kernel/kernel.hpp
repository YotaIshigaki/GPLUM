//
// kernel.cpp
//
// PPPM/TreePM EPEP/EPSP kernel using AVX-512
//
// Kentaro NOMURA 2019/10/05
//
// Known problem:
// nothing so far

#include<particle_simulator.hpp>
#include<particle_mesh.hpp>
#include<param.h>
#include<param_fdps.h>

#include <cstdio>

#include "user_defined_class.h"

#include <immintrin.h>

/* Original Cutoff functions not used in this code */
PS::F64 S2_pcut(const PS::F64 xi) {
  // This is the potential cutoff function where we used Eq.(8.75)
  // in Hockney & Eastwood (1987).

  if (xi <= 1.0) {
    return 1.0 - xi*(208.0 
                     +(xi*xi)*(-112.0 
                               +(xi*xi)*(56.0 
                                         +xi*(-14.0 
                                              +xi*(-8.0
                                                   +3.0*xi)))))/140.0;
  } else if ((1.0 < xi) && (xi < 2.0)) {
    return 1.0 - (12.0 
                  +xi*(128.0
                       +xi*(224.0 
                            +xi*(-448.0 
                                 +xi*(280.0 
                                      +xi*(-56.0 
                                           +xi*(-14.0 
                                                +xi*(8.0
                                                     -xi))))))))/140.0;
  } else {
    return 0.0;
  }
}

PS::F64 S2_fcut(const PS::F64 xi) {
  // This function returns 1 - R(\xi), where \xi is r/(a/2), a is the
  // scale length of the cutoff function, and R(\xi) is almost the same
  // as the function defined as Eq.(8-72) in Hockney & Eastwood (1987).
  // The only difference is that [1/(r/(a/2))]^2 is factored out 
  // in this function from Eq.(8-72).

  if (xi <= 1.0) {
    return 1.0 - (xi*xi*xi)*(224.0 
                             +(xi*xi)*(-224.0
                                       +xi*(70.0
                                            +xi*(48.0-21.0*xi))))/140.0;
  } else if ((1.0 < xi) && (xi < 2.0)) {
    return 1.0 - (12.0 
                  +(xi*xi)*(-224.0 
                            +xi*(896.0 
                                 +xi*(-840.0 
                                      +xi*(224.0 
                                           +xi*(70.0 
                                                +xi*(-48.0+7.0*xi)))))))/140.0;
  } else {
    return 0.0;
  }
}

#ifdef ENABLE_AVX512
__attribute__((always_inline))
static inline __m512 _mm512_rsqrt14to28_ps(const __m512 x, const __m512 y){
  const __m512 half = _mm512_set1_ps(0.5f);
  const __m512 a    = _mm512_mul_ps(y,y);
  const __m512 b    = _mm512_mul_ps(half,x);
  const __m512 hh   = _mm512_fnmadd_ps(a,b,half);
  return _mm512_fmadd_ps(y,hh,y);
}
#endif

namespace PP{
  const PS::S32 G0_ORDER = 1;
  const PS::S32 G1_ORDER = 4;
  const PS::S32 G2_ORDER = 3;
  const PS::S32 G3_ORDER = 2;

  enum{
       CalcForce,
       CalcPot,
       CalcForcePot
  };

  /* original functions to compute g0, g1, g2, g3 */
  static inline PS::F32 calc_g0(const PS::F32 r){
    const PS::F32 s = std::max(2.0-r,0.0);
    const PS::F32 s2 = s*s;
    const PS::F32 s3 = s2*s;
    const PS::F32 t = std::max(1.0-r,0.0);
    const PS::F32 t2 = t*t;
    const PS::F32 t4 = t2*t2;
    PS::F32 s_g0 = s3*s3*(0.1 + s*(-8./140. + s/140.));
    PS::F32 t_g0 = t4*t*(t - 4.0)/35.;
    return s_g0 - t2*t_g0;
  }
  static inline PS::F32 calc_g1(const PS::F32 r){
    const PS::F32 s = std::max(2.0-r,0.0);
    const PS::F32 s2 = s*s;
    const PS::F32 s3 = s2*s;
    const PS::F32 t = std::max(1.0-r,0.0);
    const PS::F32 t2 = t*t;
    const PS::F32 t4 = t2*t2;
    const PS::F32 s_g1 = s3*s2*(-168./140. + s*(182./140. + s*(-64./140. + s*(7./140.))));
    const PS::F32 t_g1 = t4*(28.0/35. + t*(-32./35.0 + t * (7./35.)));
    return s_g1 - t2*t_g1;
  }
  static inline PS::F32 calc_g2(const PS::F32 r){
    const PS::F32 s = std::max(2.0-r,0.0);
    const PS::F32 s2 = s*s;
    const PS::F32 s3 = s2*s;
    const PS::F32 t = std::max(1.0-r,0.0);
    const PS::F32 t2 = t*t;
    const PS::F32 t4 = t2*t2;
    const PS::F32 s_g2 = s2*s2*(1680./420. + s*(-2520./420. + s*(1442./420. + s*(-368./420. + s*(35./420.)))));
    const PS::F32 t_g2 = t*t2*(-168./105. + t*(308./105. + t*(-184./105. + t*(35./105.))));
    return 3.0*(s_g2 - t2*t_g2);
  }
  static inline PS::F32 calc_g3(const PS::F32 r){
    const PS::F32 s = std::max(2.0-r,0.0);
    const PS::F32 s2 = s*s;
    const PS::F32 s3 = s2*s;
    const PS::F32 t = std::max(1.0-r,0.0);
    const PS::F32 t2 = t*t;
    const PS::F32 t4 = t2*t2;
    const PS::F32 s_g3 = s3*(-4480./700. + s*(7840./700. + s*(-5768/700. + s*(2198./700. + s*(-432./700. + s*(35./700.))))));
    const PS::F32 t_g3 = t2*(280./175. + t*(-616./175. + t*(532./175. + t*(-216./175. + t*(35./175.)))));

    return 15.0*(s_g3 - t2*t_g3);
  }

  /* function to compute g0, g1, g2, g3 using piecewise polynomial approximation */
  static inline PS::F32 calc_g0_2nd(const PS::F32 r){
    const PS::F32 tab[16][3] = {
#include "./inc/g0_2nd.h"
    };
    const PS::S32 i = r * 8;
    const PS::F32 x = r - 0.125*i;
    return tab[i][0] + x * (tab[i][1] + x * (tab[i][2]));
  }
  static inline PS::F32 calc_g1_2nd(const PS::F32 r){
    const PS::F32 tab[16][3] = {
#include "./inc/g1_2nd.h"
    };
    const PS::S32 i = r * 8;
    const PS::F32 x = r - 0.125*i;
    return tab[i][0] + x * (tab[i][1] + x * (tab[i][2]));
  }
  static inline PS::F32 calc_g2_2nd(const PS::F32 r){
    const PS::F32 tab[16][3] = {
#include "./inc/g2_2nd.h"
    };
    const PS::S32 i = r * 8;
    const PS::F32 x = r - 0.125*i;
    return tab[i][0] + x * (tab[i][1] + x * (tab[i][2]));
  }
  static inline PS::F32 calc_g3_2nd(const PS::F32 r){
    const PS::F32 tab[16][3] = {
#include "./inc/g3_2nd.h"
    };
    const PS::S32 i = r * 8;
    const PS::F32 x = r - 0.125*i;
    return tab[i][0] + x * (tab[i][1] + x * (tab[i][2]));
  }

  static inline PS::F32 calc_g0_3rd(const PS::F32 r){
    const PS::F32 tab[16][4] = {
#include "./inc/g0_3rd.h"
    };
    const PS::S32 i = r * 8;
    const PS::F32 x = r - 0.125*i;
    return tab[i][0] + x * (tab[i][1] + x * (tab[i][2] + x*(tab[i][3])));
  }
  static inline PS::F32 calc_g1_3rd(const PS::F32 r){
    const PS::F32 tab[16][4] = {
#include "./inc/g1_3rd.h"
    };
    const PS::S32 i = r * 8;
    const PS::F32 x = r - 0.125*i;
    return tab[i][0] + x * (tab[i][1] + x * (tab[i][2] + x*(tab[i][3])));
  }
  static inline PS::F32 calc_g2_3rd(const PS::F32 r){
    const PS::F32 tab[16][4] = {
#include "./inc/g2_3rd.h"
    };
    const PS::S32 i = r * 8;
    const PS::F32 x = r - 0.125*i;
    return tab[i][0] + x * (tab[i][1] + x * (tab[i][2] + x*(tab[i][3])));
  }
  static inline PS::F32 calc_g3_3rd(const PS::F32 r){
    const PS::F32 tab[16][4] = {
#include "./inc/g3_3rd.h"
    };
    const PS::S32 i = r * 8;
    const PS::F32 x = r - 0.125*i;
    return tab[i][0] + x * (tab[i][1] + x * (tab[i][2] + x*(tab[i][3])));
  }

  static inline PS::F32 calc_g0_4th(const PS::F32 r){
    const PS::F32 tab[16][5] = {
#include "./inc/g0_4th.h"
    };
    const PS::S32 i = r * 8;
    const PS::F32 x = r - 0.125*i;
    return tab[i][0] + x * (tab[i][1] + x * (tab[i][2] + x*(tab[i][3] + x*(tab[i][4]))));
  }
  static inline PS::F32 calc_g1_4th(const PS::F32 r){
    const PS::F32 tab[16][5] = {
#include "./inc/g1_4th.h"
    };
    const PS::S32 i = r * 8;
    const PS::F32 x = r - 0.125*i;
    return tab[i][0] + x * (tab[i][1] + x * (tab[i][2] + x*(tab[i][3] + x*(tab[i][4]))));
  }
  static inline PS::F32 calc_g2_4th(const PS::F32 r){
    const PS::F32 tab[16][5] = {
#include "./inc/g2_4th.h"
    };
    const PS::S32 i = r * 8;
    const PS::F32 x = r - 0.125*i;
    return tab[i][0] + x * (tab[i][1] + x * (tab[i][2] + x*(tab[i][3] + x*(tab[i][4]))));
  }
  static inline PS::F32 calc_g3_4th(const PS::F32 r){
    const PS::F32 tab[16][5] = {
#include "./inc/g3_4th.h"
    };
    const PS::S32 i = r * 8;
    const PS::F32 x = r - 0.125*i;
    return tab[i][0] + x * (tab[i][1] + x * (tab[i][2] + x*(tab[i][3] + x*(tab[i][4]))));
  }

  PS::F32 calc_g0(const PS::F32 r,const int i){
    switch(i){
    case 1:
      return calc_g0(r);
    case 2:
      return calc_g0_2nd(r);
    case 3:
      return calc_g0_3rd(r);
    case 4:
      return calc_g0_4th(r);
    default:
      fprintf(stderr,"error: unsuported order of polynomial approximation for g0 calculation\n");
      PS::Abort();
      return 0;
    }
  }
  PS::F32 calc_g1(const PS::F32 r,const int i){
    switch(i){
    case 1:
      return calc_g1(r);
    case 2:
      return calc_g1_2nd(r);
    case 3:
      return calc_g1_3rd(r);
    case 4:
      return calc_g1_4th(r);
    default:
      fprintf(stderr,"error: unsuported order of polynomial approximation for g1 calculation\n");
      PS::Abort();
      return 0;
    }
  }
  PS::F32 calc_g2(const PS::F32 r,const int i){
    switch(i){
    case 1:
      return calc_g2(r);
    case 2:
      return calc_g2_2nd(r);
    case 3:
      return calc_g2_3rd(r);
    case 4:
      return calc_g2_4th(r);
    default:
      fprintf(stderr,"error: unsuported order of polynomial approximation for g2 calculation\n");
      PS::Abort();
      return 0;
    }
  }
  PS::F32 calc_g3(const PS::F32 r,const int i){
    switch(i){
    case 1:
      return calc_g3(r);
    case 2:
      return calc_g3_2nd(r);
    case 3:
      return calc_g3_3rd(r);
    case 4:
      return calc_g3_4th(r);
    default:
      fprintf(stderr,"error: unsuported order of polynomial approximation for g3 calculation\n");
      PS::Abort();
      return 0;
    }
  }

#ifdef ENABLE_AVX512
  // for debug
  int __m512_to_int16(__m512i a,int index){
    long long tmp = a[index/2];
    if(index%2==0) return (int)(tmp & 0xffffffff);
    else           return (int)((tmp & 0xffffffff00000000)>>32);
  }
#endif

  template <bool enableAVX512>
  struct CalcForceEpEp{
    const PS::F32 rcut,rcut2;
    const PS::F32 eps;
    CalcForceEpEp(const PS::F32 _rcut,const PS::F32 _eps = 0.001)
      :        rcut(_rcut), rcut2(_rcut*_rcut),eps(_eps){}

    void operator () (const EP*     ep_i,
                      const PS::S32 n_ip,
                      const EP*     ep_j,
                      const PS::S32 n_jp,
                      Force*  force){
      int i = 0;
      if(enableAVX512){
      //printf("--- AVX-512 kernel (ni = %d) ---\n",n_ip);
      PS::F32 g1_tab[16][5] = {
#include "./inc/g1_4th.h"
      };
      assert(G1_ORDER == 4);
      const __m512 vrc_sq = _mm512_set1_ps(rcut2);
      const __m512 vrc_inv = _mm512_set1_ps(1.f / rcut);
      const __m512 zero = _mm512_set1_ps(0.f);
      const __m512 one  = _mm512_set1_ps(1.f);
      const __m512 two  = _mm512_set1_ps(2.f);
      const __m512 veps  = _mm512_set1_ps(eps);

      int index_4th[16] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75};
      const __m512 g1_tab0 = _mm512_i32gather_ps(_mm512_load_epi32(index_4th),(float*)(g1_tab),4);
      for(int j=0;j<16;j++) index_4th[j]++;
      const __m512 g1_tab1 = _mm512_i32gather_ps(_mm512_load_epi32(index_4th),(float*)(g1_tab),4);
      for(int j=0;j<16;j++) index_4th[j]++;
      const __m512 g1_tab2 = _mm512_i32gather_ps(_mm512_load_epi32(index_4th),(float*)(g1_tab),4);
      for(int j=0;j<16;j++) index_4th[j]++;
      const __m512 g1_tab3 = _mm512_i32gather_ps(_mm512_load_epi32(index_4th),(float*)(g1_tab),4);
      for(int j=0;j<16;j++) index_4th[j]++;
      const __m512 g1_tab4 = _mm512_i32gather_ps(_mm512_load_epi32(index_4th),(float*)(g1_tab),4);
      const __m512 dh = _mm512_set1_ps(0.125f);

      const int index[16] = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60};
      const __m512i i0 = _mm512_load_epi32(index);
      const __m512i i1 = _mm512_add_epi32(i0,_mm512_set1_epi32(1));
      const __m512i i2 = _mm512_add_epi32(i0,_mm512_set1_epi32(2));

      const int zero2fifteen[16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
      const __m512i vzero2fifteen = _mm512_load_epi32(zero2fifteen);
      const int n_ip_avx = ((n_ip+15)/16)*16;
      for(;i<n_ip_avx;i+=16){
        __mmask16 mask_i = _mm512_cmp_epi32_mask(_mm512_add_epi32(_mm512_set1_epi32(i),vzero2fifteen),_mm512_set1_epi32(n_ip),_MM_CMPINT_LT);
        const __m512 rix  = _mm512_mask_i32gather_ps(zero,mask_i,i0,&ep_i[i],4);
        const __m512 riy  = _mm512_mask_i32gather_ps(zero,mask_i,i1,&ep_i[i],4);
        const __m512 riz  = _mm512_mask_i32gather_ps(zero,mask_i,i2,&ep_i[i],4);

        __m512d fix_lo = _mm512_castps_pd(zero);
        __m512d fiy_lo = _mm512_castps_pd(zero);
        __m512d fiz_lo = _mm512_castps_pd(zero);
        __m512d fix_hi = _mm512_castps_pd(zero);
        __m512d fiy_hi = _mm512_castps_pd(zero);
        __m512d fiz_hi = _mm512_castps_pd(zero);

        for(int j=0;j<n_jp;j++){
          const PS::F32vec rj = ep_j[j].pos;
          const __m512 rjx = _mm512_set1_ps(rj.x);
          const __m512 rjy = _mm512_set1_ps(rj.y);
          const __m512 rjz = _mm512_set1_ps(rj.z);
          const __m512 mj  = _mm512_set1_ps(ep_j[j].charge);

          const __m512 rijx = _mm512_sub_ps(rix,rjx);
          const __m512 rijy = _mm512_sub_ps(riy,rjy);
          const __m512 rijz = _mm512_sub_ps(riz,rjz); // 3 sub, 3mul
          __m512 r2 = _mm512_fmadd_ps(rijx,rijx,veps);
          r2 = _mm512_fmadd_ps(rijy,rijy,r2);
          r2 = _mm512_fmadd_ps(rijz,rijz,r2); // 1 mul, 2 fmadd
          __m512 ri = _mm512_rsqrt14_ps(r2);
          ri = _mm512_rsqrt14to28_ps(r2,ri);
          const __m512 r2i = _mm512_mul_ps(ri,ri); // 1 rsqrt, 3 mul, 2 fmadd

          const __m512 r = r2 * ri;

          const __m512 eta = _mm512_mul_ps(_mm512_mul_ps(r,two),vrc_inv);
          const __m512i k = _mm512_cvt_roundps_epi32(_mm512_mul_ps(eta,_mm512_set1_ps(8.0)),(_MM_FROUND_TO_NEG_INF |_MM_FROUND_NO_EXC));
          const __m512 x = _mm512_mul_ps(dh,_mm512_cvt_roundepi32_ps(k,(_MM_FROUND_TO_NEG_INF |_MM_FROUND_NO_EXC)));
          const __m512 d = _mm512_sub_ps(eta,x); //4 mul, 1 sub

          __m512 g1 = _mm512_permutexvar_ps(k,g1_tab4);
          g1 = _mm512_fmadd_ps(d,g1,_mm512_permutexvar_ps(k,g1_tab3));
          g1 = _mm512_fmadd_ps(d,g1,_mm512_permutexvar_ps(k,g1_tab2));
          g1 = _mm512_fmadd_ps(d,g1,_mm512_permutexvar_ps(k,g1_tab1));
          g1 = _mm512_fmadd_ps(d,g1,_mm512_permutexvar_ps(k,g1_tab0));// 4 fmadd

          __m512 f = _mm512_mul_ps(mj,g1);
          f = _mm512_mul_ps(ri,f); // 2 mul

          const __m512 fix = _mm512_mul_ps(rijx,f);
          const __m512d xlo = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(fix),0)));
          const __m512d xhi = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(fix),1)));
          fix_lo = _mm512_add_pd(fix_lo, xlo);
          fix_hi = _mm512_add_pd(fix_hi, xhi);
          const __m512 fiy = _mm512_mul_ps(rijy,f);
          const __m512d ylo = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(fiy),0)));
          const __m512d yhi = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(fiy),1)));
          fiy_lo = _mm512_add_pd(fiy_lo, ylo);
          fiy_hi = _mm512_add_pd(fiy_hi, yhi);
          const __m512 fiz = _mm512_mul_ps(rijz,f);
          const __m512d zlo = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(fiz),0)));
          const __m512d zhi = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(fiz),1)));
          fiz_lo = _mm512_add_pd(fiz_lo, zlo);
          fiz_hi = _mm512_add_pd(fiz_hi, zhi);
        }
        const long findex[8] = {0,1,2,3,4,5,6,7};
        const __m512i vfindex = _mm512_load_epi64(findex);
        const __m512i i0f = _mm512_sll_epi64(vfindex,_mm_set1_epi64x(2));
        const __m512i i1f = _mm512_add_epi64(i0f,_mm512_set1_epi64(1));
        const __m512i i2f = _mm512_add_epi64(i0f,_mm512_set1_epi64(2));
        const __m512i vi_lo = _mm512_add_epi64(_mm512_set1_epi64(i),vfindex);
        const __m512i vi_hi = _mm512_add_epi64(_mm512_set1_epi64(i+8),vfindex);
        const __mmask8 mask_lo = _mm512_cmp_epi64_mask(vi_lo,_mm512_set1_epi64(n_ip),_MM_CMPINT_LT);
        const __mmask8 mask_hi = _mm512_cmp_epi64_mask(vi_hi,_mm512_set1_epi64(n_ip),_MM_CMPINT_LT);
          
        _mm512_mask_i64scatter_pd(&force[i],  mask_lo,i0f,fix_lo,8);
        _mm512_mask_i64scatter_pd(&force[i],  mask_lo,i1f,fiy_lo,8);
        _mm512_mask_i64scatter_pd(&force[i],  mask_lo,i2f,fiz_lo,8);
        _mm512_mask_i64scatter_pd(&force[i+8],mask_hi,i0f,fix_hi,8);
        _mm512_mask_i64scatter_pd(&force[i+8],mask_hi,i1f,fiy_hi,8);
        _mm512_mask_i64scatter_pd(&force[i+8],mask_hi,i2f,fiz_hi,8);
      } // i loop
      }// enableAVX512
      for(; i<n_ip; i++){
        const PS::F32vec ri = ep_i[i].pos;
        PS::F32vec fi = 0.0;
        PS::F32  poti = 0.0;
        for(PS::S32 j=0; j<n_jp; j++){
          PS::F32vec rij = ri - ep_j[j].pos;
          const PS::F32 r2 = eps + rij*rij;
          if(r2 < rcut2){
            const PS::F32 r_inv = 1.f / sqrt(r2);
            const PS::F32 r = r2 * r_inv;
            const PS::F32 xi = 2.f * r / rcut;
            const PS::F32 g1 = calc_g1(xi,G1_ORDER);
            fi += (ep_j[j].charge * g1 * r_inv * r_inv * r_inv)*rij;
          }
        }
        force[i].acc += fi;
      }
    }
  };

  template <bool enableAVX512>
  struct CalcForceEpSp{
    const PS::F32 rcut,rcut2;
    const PS::F32 eps;
    CalcForceEpSp(const PS::F32 _rcut,const PS::F32 _eps = 0.001)
      : rcut(_rcut), rcut2(_rcut*_rcut), eps(_eps){}

    void operator () (const EP*     ep_i,
                      const PS::S32 n_ip,
                      const SPJQuadrupoleCutoff* sp_j,
                      const PS::S32 n_jp,
                      Force*  force){
      int i = 0;
      if(enableAVX512){
      PS::F32 g1_tab[16][5] = {
#include "./inc/g1_4th.h"
      };
      PS::F32 g2_tab[16][4] = {
#include "./inc/g2_3rd.h"
      };
      PS::F32 g3_tab[16][3] = {
#include "./inc/g3_2nd.h"
      };
      assert(G1_ORDER==4 && G2_ORDER == 3 && G3_ORDER == 2);
      const __m512 vrc_sq = _mm512_set1_ps(rcut2);
      const __m512 vrc_inv = _mm512_set1_ps(1.f / rcut);
      const __m512 zero = _mm512_set1_ps(0.0f);
      const __m512 half = _mm512_set1_ps(0.5f);
      const __m512 two  = _mm512_set1_ps(2.0f);
      const __m512 veps = _mm512_set1_ps(eps);

      int index_4th[16] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75};
      const __m512 g1_tab0 = _mm512_i32gather_ps(_mm512_load_epi32(index_4th),(float*)(g1_tab),4);
      for(int j=0;j<16;j++) index_4th[j]++;
      const __m512 g1_tab1 = _mm512_i32gather_ps(_mm512_load_epi32(index_4th),(float*)(g1_tab),4);
      for(int j=0;j<16;j++) index_4th[j]++;
      const __m512 g1_tab2 = _mm512_i32gather_ps(_mm512_load_epi32(index_4th),(float*)(g1_tab),4);
      for(int j=0;j<16;j++) index_4th[j]++;
      const __m512 g1_tab3 = _mm512_i32gather_ps(_mm512_load_epi32(index_4th),(float*)(g1_tab),4);
      for(int j=0;j<16;j++) index_4th[j]++;
      const __m512 g1_tab4 = _mm512_i32gather_ps(_mm512_load_epi32(index_4th),(float*)(g1_tab),4);

      int index_3rd[16] = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60};
      const __m512 g2_tab0 = _mm512_i32gather_ps(_mm512_load_epi32(index_3rd),(float*)(g2_tab),4);
      for(int j=0;j<16;j++) index_3rd[j]++;
      const __m512 g2_tab1 = _mm512_i32gather_ps(_mm512_load_epi32(index_3rd),(float*)(g2_tab),4);
      for(int j=0;j<16;j++) index_3rd[j]++;
      const __m512 g2_tab2 = _mm512_i32gather_ps(_mm512_load_epi32(index_3rd),(float*)(g2_tab),4);
      for(int j=0;j<16;j++) index_3rd[j]++;
      const __m512 g2_tab3 = _mm512_i32gather_ps(_mm512_load_epi32(index_3rd),(float*)(g2_tab),4);

      int index_2nd[16] = {0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45};
      const __m512 g3_tab0 = _mm512_i32gather_ps(_mm512_load_epi32(index_2nd),(float*)(g3_tab),4);
      for(int j=0;j<16;j++) index_2nd[j]++;
      const __m512 g3_tab1 = _mm512_i32gather_ps(_mm512_load_epi32(index_2nd),(float*)(g3_tab),4);
      for(int j=0;j<16;j++) index_2nd[j]++;
      const __m512 g3_tab2 = _mm512_i32gather_ps(_mm512_load_epi32(index_2nd),(float*)(g3_tab),4);

      const __m512 dh = _mm512_set1_ps(0.125f);
      const int index[16] = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60};
      const __m512i i0 = _mm512_load_epi32(index);
      const __m512i i1 = _mm512_add_epi32(i0,_mm512_set1_epi32(1));
      const __m512i i2 = _mm512_add_epi32(i0,_mm512_set1_epi32(2));
      const int zero2fifteen[16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
      const __m512i vzero2fifteen = _mm512_load_epi32(zero2fifteen);
      const int n_ip_avx = ((n_ip+15)/16)*16;
      for(;i<n_ip_avx;i+=16){
        __mmask16 mask_i = _mm512_cmp_epi32_mask(_mm512_add_epi32(_mm512_set1_epi32(i),vzero2fifteen),_mm512_set1_epi32(n_ip),_MM_CMPINT_LT);
        const __m512 rix  = _mm512_mask_i32gather_ps(zero,mask_i,i0,&ep_i[i],4);
        const __m512 riy  = _mm512_mask_i32gather_ps(zero,mask_i,i1,&ep_i[i],4);
        const __m512 riz  = _mm512_mask_i32gather_ps(zero,mask_i,i2,&ep_i[i],4);
        __m512d fix_lo = _mm512_castps_pd(zero);
        __m512d fiy_lo = _mm512_castps_pd(zero);
        __m512d fiz_lo = _mm512_castps_pd(zero);
        __m512d fix_hi = _mm512_castps_pd(zero);
        __m512d fiy_hi = _mm512_castps_pd(zero);
        __m512d fiz_hi = _mm512_castps_pd(zero);
        for(int j=0;j<n_jp;j++){
          const PS::F32vec rj = sp_j[j].pos;
          const __m512 rjx = _mm512_set1_ps(rj.x);
          const __m512 rjy = _mm512_set1_ps(rj.y);
          const __m512 rjz = _mm512_set1_ps(rj.z);
          const PS::F32vec dj = sp_j[j].dipole;
          const __m512 djx = _mm512_set1_ps(dj.x);
          const __m512 djy = _mm512_set1_ps(dj.y);
          const __m512 djz = _mm512_set1_ps(dj.z);
          const PS::F32mat qj = sp_j[j].quadrupole;
          const __m512 qjxx = _mm512_set1_ps(qj.xx);
          const __m512 qjyy = _mm512_set1_ps(qj.yy);
          const __m512 qjzz = _mm512_set1_ps(qj.zz);
          const __m512 qjxy = _mm512_set1_ps(qj.xy);
          const __m512 qjxz = _mm512_set1_ps(qj.xz);
          const __m512 qjyz = _mm512_set1_ps(qj.yz);
          const __m512 trqh = _mm512_mul_ps(_mm512_add_ps(qjxx,_mm512_add_ps(qjyy,qjzz)),half);
          // 2 add, 1 mul

          const __m512 rijx = _mm512_sub_ps(rix,rjx);
          const __m512 rijy = _mm512_sub_ps(riy,rjy);
          const __m512 rijz = _mm512_sub_ps(riz,rjz);
          // 3 sub, 3 mul

          const __m512 dr  = _mm512_fmadd_ps(rijz,djz, _mm512_fmadd_ps(rijy,djy, _mm512_mul_ps(rijx,djx)));
          const __m512 qrx = _mm512_fmadd_ps(rijz,qjxz,_mm512_fmadd_ps(rijy,qjxy,_mm512_mul_ps(rijx,qjxx)));
          const __m512 qry = _mm512_fmadd_ps(rijz,qjyz,_mm512_fmadd_ps(rijy,qjyy,_mm512_mul_ps(rijx,qjxy)));
          const __m512 qrz = _mm512_fmadd_ps(rijz,qjzz,_mm512_fmadd_ps(rijy,qjyz,_mm512_mul_ps(rijx,qjxz)));
          const __m512 rqr = _mm512_fmadd_ps(rijz,qrz, _mm512_fmadd_ps(rijy,qry, _mm512_mul_ps(rijx,qrx)));
          const __m512 rqrh= _mm512_mul_ps(rqr,half);
          // 6 mul, 10 fmadd

          __m512 r2 = _mm512_fmadd_ps(rijx,rijx,veps);
          r2 = _mm512_fmadd_ps(rijy,rijy,r2);
          r2 = _mm512_fmadd_ps(rijz,rijz,r2);
          __m512 ri = _mm512_rsqrt14_ps(r2);
          ri = _mm512_rsqrt14to28_ps(r2,ri);
          //ri = _mm512_rsqrt14to28_ps(r2,ri);
          const __m512 r2i = _mm512_mul_ps(ri,ri);
          const __m512 r3i = _mm512_mul_ps(ri,r2i);
          // 5 mul, 4 fmadd, 1 rsqrt

          // coulomb
          __mmask16 mask_cl = _mm512_cmp_ps_mask(r2,vrc_sq,_CMP_LT_OS);
          const __m512 m = _mm512_set1_ps(sp_j[j].getCharge());
          const __m512 r = r2 * ri;
          const __m512 eta = _mm512_mul_ps(_mm512_mul_ps(r,two),vrc_inv);
          const __m512i k = _mm512_cvt_roundps_epi32(_mm512_mul_ps(eta,_mm512_set1_ps(8.0)),(_MM_FROUND_TO_NEG_INF |_MM_FROUND_NO_EXC));
          const __m512 x = _mm512_mul_ps(dh,_mm512_cvt_roundepi32_ps(k,(_MM_FROUND_TO_NEG_INF |_MM_FROUND_NO_EXC)));
          const __m512 d = _mm512_sub_ps(eta,x);
          // 1 sub, 4 mul

          __m512 g1 = _mm512_permutexvar_ps(k,g1_tab4);
          g1 = _mm512_fmadd_ps(d,g1,_mm512_permutexvar_ps(k,g1_tab3));
          g1 = _mm512_fmadd_ps(d,g1,_mm512_permutexvar_ps(k,g1_tab2));
          g1 = _mm512_fmadd_ps(d,g1,_mm512_permutexvar_ps(k,g1_tab1));
          g1 = _mm512_fmadd_ps(d,g1,_mm512_permutexvar_ps(k,g1_tab0));
          __m512 g2 = _mm512_permutexvar_ps(k,g2_tab3);
          g2 = _mm512_fmadd_ps(d,g2,_mm512_permutexvar_ps(k,g2_tab2));
          g2 = _mm512_fmadd_ps(d,g2,_mm512_permutexvar_ps(k,g2_tab1));
          g2 = _mm512_fmadd_ps(d,g2,_mm512_permutexvar_ps(k,g2_tab0));
          __m512 g3 = _mm512_permutexvar_ps(k,g3_tab2);
          g3 = _mm512_fmadd_ps(d,g3,_mm512_permutexvar_ps(k,g3_tab1));
          g3 = _mm512_fmadd_ps(d,g3,_mm512_permutexvar_ps(k,g3_tab0));
          // 9 fmadd

          //if(PS::Comm::getRank()==0)for(int x=0;x<16;x++) printf("%e %d %e\n",eta[x],__m512_to_int16(k,x),g1[x]);
          const __m512 tmp0 = _mm512_mask_mul_ps(zero,mask_cl,g3,rqrh);
          const __m512 tmp1 = _mm512_mask_mul_ps(zero,mask_cl,g2,_mm512_sub_ps(dr,trqh));
          const __m512 tmp2 = _mm512_mask_mul_ps(zero,mask_cl,g1,m);
          const __m512 fr = _mm512_fmsub_ps(_mm512_fnmadd_ps(tmp0,r2i,tmp1),r2i,tmp2);
          const __m512 fq = _mm512_mask_mul_ps(zero,mask_cl,g2,r2i);
          // 1 sub, 4 mul, 2 fmadd
          __m512 fx = _mm512_mask_mul_ps(zero,mask_cl,g1,djx);
          fx = _mm512_fmadd_ps(fr,rijx,fx);
          fx = _mm512_fnmadd_ps(fq,qrx,fx);
          __m512 fy = _mm512_mask_mul_ps(zero,mask_cl,g1,djy);
          fy = _mm512_fmadd_ps(fr,rijy,fy);
          fy = _mm512_fnmadd_ps(fq,qry,fy);
          __m512 fz = _mm512_mask_mul_ps(zero,mask_cl,g1,djz);
          fz = _mm512_fmadd_ps(fr,rijz,fz);
          fz = _mm512_fnmadd_ps(fq,qrz,fz);
          // 3 mul, 6 fmadd

          __m512 fix = _mm512_mul_ps(fx,r3i);
          __m512 fiy = _mm512_mul_ps(fy,r3i);
          __m512 fiz = _mm512_mul_ps(fz,r3i);
          // 3 fmadd
          const __m512d xlo = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(fix),0)));
          const __m512d xhi = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(fix),1)));
          fix_lo = _mm512_add_pd(fix_lo, xlo);
          fix_hi = _mm512_add_pd(fix_hi, xhi);
          const __m512d ylo = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(fiy),0)));
          const __m512d yhi = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(fiy),1)));
          fiy_lo = _mm512_add_pd(fiy_lo, ylo);
          fiy_hi = _mm512_add_pd(fiy_hi, yhi);
          const __m512d zlo = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(fiz),0)));
          const __m512d zhi = _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(fiz),1)));
          fiz_lo = _mm512_add_pd(fiz_lo, zlo);
          fiz_hi = _mm512_add_pd(fiz_hi, zhi);
        }
        const long findex[8] = {0,1,2,3,4,5,6,7};
        const __m512i vfindex = _mm512_load_epi64(findex);
        const __m512i i0f = _mm512_sll_epi64(vfindex,_mm_set1_epi64x(2));
        const __m512i i1f = _mm512_add_epi64(i0f,_mm512_set1_epi64(1));
        const __m512i i2f = _mm512_add_epi64(i0f,_mm512_set1_epi64(2));
        const __m512i vi_lo = _mm512_add_epi64(_mm512_set1_epi64(i),vfindex);
        const __m512i vi_hi = _mm512_add_epi64(_mm512_set1_epi64(i+8),vfindex);
        const __mmask8 mask_lo = _mm512_cmp_epi64_mask(vi_lo,_mm512_set1_epi64(n_ip),_MM_CMPINT_LT);
        const __mmask8 mask_hi = _mm512_cmp_epi64_mask(vi_hi,_mm512_set1_epi64(n_ip),_MM_CMPINT_LT);
        __m512d fxtmp_lo  = _mm512_mask_i64gather_pd(_mm512_castps_pd(zero),mask_lo,i0f,&force[i],8);
        __m512d fytmp_lo  = _mm512_mask_i64gather_pd(_mm512_castps_pd(zero),mask_lo,i1f,&force[i],8);
        __m512d fztmp_lo  = _mm512_mask_i64gather_pd(_mm512_castps_pd(zero),mask_lo,i2f,&force[i],8);
        fxtmp_lo = _mm512_add_pd(fix_lo,fxtmp_lo);
        fytmp_lo = _mm512_add_pd(fiy_lo,fytmp_lo);
        fztmp_lo = _mm512_add_pd(fiz_lo,fztmp_lo);
        _mm512_mask_i64scatter_pd(&force[i],  mask_lo,i0f,fxtmp_lo,8);
        _mm512_mask_i64scatter_pd(&force[i],  mask_lo,i1f,fytmp_lo,8);
        _mm512_mask_i64scatter_pd(&force[i],  mask_lo,i2f,fztmp_lo,8);

        __m512d fxtmp_hi  = _mm512_mask_i64gather_pd(_mm512_castps_pd(zero),mask_hi,i0f,&force[i+8],8);
        __m512d fytmp_hi  = _mm512_mask_i64gather_pd(_mm512_castps_pd(zero),mask_hi,i1f,&force[i+8],8);
        __m512d fztmp_hi  = _mm512_mask_i64gather_pd(_mm512_castps_pd(zero),mask_hi,i2f,&force[i+8],8);
        fxtmp_hi = _mm512_add_pd(fix_hi,fxtmp_hi);
        fytmp_hi = _mm512_add_pd(fiy_hi,fytmp_hi);
        fztmp_hi = _mm512_add_pd(fiz_hi,fztmp_hi);
        _mm512_mask_i64scatter_pd(&force[i+8],mask_hi,i0f,fxtmp_hi,8);
        _mm512_mask_i64scatter_pd(&force[i+8],mask_hi,i1f,fytmp_hi,8);
        _mm512_mask_i64scatter_pd(&force[i+8],mask_hi,i2f,fztmp_hi,8);
      } // loop i
      } // enableAVX512
      for(; i<n_ip; i++){
        const PS::F32vec ri = ep_i[i].pos;
        PS::F64vec fi = 0.0;
        for(PS::S32 j=0; j<n_jp; j++){
          PS::F32vec rij = ri - sp_j[j].pos;
          const PS::F32 r2 = rij*rij;
          if(r2 < rcut2){
            const PS::F32 r_inv = 1.0 / sqrt(r2);
            const PS::F32 r2_inv = r_inv * r_inv;
            const PS::F32 r = r2 * r_inv;
            const PS::F32 xi = 2.0 * r / rcut;

            const PS::F32  mono = sp_j[j].charge;
            const PS::F32vec di   = sp_j[j].dipole;
            PS::F32mat quad = sp_j[j].quadrupole;

            const PS::F32 dr = di * rij;
            const PS::F32vec qr = PS::F32vec(quad.xx*rij.x + quad.xy*rij.y + quad.xz*rij.z,
                                             quad.xy*rij.x + quad.yy*rij.y + quad.yz*rij.z,
                                             quad.xz*rij.x + quad.yz*rij.y + quad.zz*rij.z); // F64mat is symmetric
            const PS::F32 trq = quad.getTrace();

            PS::F32 g1,g2,g3;
            g1 = calc_g1(xi,G1_ORDER);
            g2 = calc_g2(xi,G2_ORDER);
            g3 = calc_g3(xi,G3_ORDER);
            const PS::F32vec fm = - g1 * mono * rij;
            const PS::F32vec fd = (g2 * dr * r2_inv) * rij + g1 * di;
            const PS::F32vec fq = - 0.5*r2_inv * ((g3*(rij*qr)*r2_inv + g2*trq)*rij + 2.0*g2*qr);
            fi += (r2_inv*r_inv)*(fm + fd + fq);
          } // cutoff
        }
        force[i].acc += fi;
      }
    }
  };
}; // namespace PP

class CalcForceAll{
public:
  PS::TreeForForceLong<Force,EP,EP,MomentQuadrupoleCutoff,SPJQuadrupoleCutoff>::WithCutoff pp;
  PS::PM::ParticleMesh pm;
  PS::F64 box_size, box_size_inv;
  bool doTimeCount = true;
  void initialize(const PS::S32 n, const PS::F64 t,const PS::U32 l,PS::U32 g){
    pp.initialize(n,t,l,g);
  }
  void clearTimeProfile() {
    pp.clearTimeProfile();
  }
  PS::TimeProfile getTimeProfile() {
    return pp.getTimeProfile();
  }
  template <class Tfunc_ep_ep,class Tfunc_ep_sp,class Tpsys,class Tdinfo>
  void calcForceAllAndWriteBack(Tfunc_ep_ep epep,
                                Tfunc_ep_sp epsp,
                                Tpsys& psys,
                                Tdinfo& dinfo,
                                const bool clear = true){
    pp.calcForceAllAndWriteBack(epep,epsp,psys,dinfo,clear);
#if 0
    pm.setDomainInfoParticleMesh(dinfo);
    pm.setParticleParticleMesh(psys);
    pm.calcMeshForceOnly();
    for (PS::S32 i=0; i<psys.getNumberOfParticleLocal(); i++) { 
      PS::F64vec pos    = psys[i].pos;
      PS::F64    charge = psys[i].charge;
      psys[i].acc -= charge * pm.getForce(pos);
      psys[i].pot -= charge * pm.getPotential(pos);
    }
#endif
  }
};


