#pragma once

#include <immintrin.h>

struct v4sf;
struct v8sf;
struct v2df;
struct v4df;
#ifdef __AVX512F__
struct v16sf;
struct v8df;
#endif

struct v4sf {
    typedef float _v4sf __attribute__((vector_size(16))) __attribute__((aligned(16)));
    _v4sf val;
    
    static inline int getVectorLength() {
        return 4;
    }
    
    v4sf(const v4sf & rhs) : val(rhs.val) {}
    v4sf operator = (const v4sf rhs) {
        val = rhs.val;
        return (*this);
    }
    v4sf() : val(_mm_setzero_ps()) {}
    v4sf(const float x) : val(_mm_set1_ps(x)) {}
    v4sf(const float x, const float y, const float z, const float w)
        : val(_mm_set_ps(w, z, y, x)) {}
    v4sf(const _v4sf _val) : val(_val) {}
    operator _v4sf() {return val;}
    
    v4sf operator + (const v4sf rhs) const {
        return v4sf(val + rhs.val);
    }
    v4sf operator - (const v4sf rhs) const {
        return v4sf(val - rhs.val);
    }
    v4sf operator * (const v4sf rhs) const {
        return v4sf(val * rhs.val);
    }
    v4sf operator / (const v4sf rhs) const {
        return v4sf(_mm_div_ps(val, rhs.val));
    }
    
    v4sf operator += (const v4sf rhs) {
        val = val + rhs.val;
        return (*this);
    }
    v4sf operator -= (const v4sf rhs) {
        val = val - rhs.val;
        return (*this);
    }
    v4sf operator *= (const v4sf rhs) {
        val = val * rhs.val;
        return (*this);
    }
    v4sf operator /= (const v4sf rhs) {
        val = _mm_div_ps(val, rhs.val);
        return (*this);
    }

    v4df cvtps2pd();

    static inline v4sf rcp_0th(const v4sf rhs) {
        return _mm_rcp_ps(rhs.val);
    }
    static inline v4sf rcp_1st(const v4sf rhs) {
        v4sf x0 = _mm_rcp_ps(rhs.val);
        v4sf h  = v4sf(1.) - rhs * x0;
        return x0 + h * x0;
    }
    
    static inline v4sf sqrt(const v4sf rhs) {
        return v4sf(_mm_sqrt_ps(rhs.val));
    }
    static inline v4sf rsqrt_0th(const v4sf rhs) {
        return v4sf(_mm_rsqrt_ps(rhs.val));
    }
    static inline v4sf rsqrt_1st(const v4sf rhs) {
        v4sf x0 = v4sf(_mm_rsqrt_ps(rhs.val));
        v4sf h  = v4sf(1.) - rhs * x0 * x0;
        return x0 + v4sf(0.5) * h * x0;
    }
    static inline v4sf rsqrt_1st_phantom(const v4sf rhs) {
        v4sf x0 = v4sf(_mm_rsqrt_ps(rhs.val));
        return v4sf(x0 * (rhs * x0 * x0 - v4sf(3.)));
    }
    static inline v4sf hadd(v4sf x0, v4sf x1) {
        return _mm_hadd_ps(x0, x1);
    }
    
    inline void store(float *p) const {
        _mm_store_ps(p, val);
    }
    inline void load(float const *p) {
        val = _mm_load_ps(p);
    }
    
    void print(FILE * fp = stdout,
               const char * fmt = "%+e %+e %+e %+e\n") const {
        float a[4] __attribute__((aligned(16)));
        store(a);
        fprintf(fp, fmt, a[0], a[1], a[2], a[3]);
    }    
};

struct v8sf {
    typedef float _v8sf __attribute__((vector_size(32))) __attribute__((aligned(32)));
    _v8sf val;
    
    static inline int getVectorLength() {
        return 8;
    }
    
    v8sf(const v8sf & rhs) : val(rhs.val) {}
    v8sf operator = (const v8sf rhs) {
        val = rhs.val;
        return (*this);
    }
    v8sf() : val(_mm256_setzero_ps()) {}
    v8sf(const float x) : val(_mm256_set1_ps(x)) {}
    v8sf(const float x0, const float y0, const float z0, const float w0,
         const float x1, const float y1, const float z1, const float w1)
        : val(_mm256_set_ps(w1, z1, y1, x1, w0, z0, y0, x0)) {}
    v8sf(const _v8sf _val) : val(_val) {}
    operator _v8sf() {return val;}
    v8sf(const v4sf rhs) {
        val = _mm256_broadcast_ps((const __m128*)&rhs);
    }
    v8sf(v4sf rhs0, v4sf rhs1) {
        v8sf temp(rhs1);
        val = _mm256_insertf128_ps(temp, rhs0, 0);
    }
    
    v8sf operator + (const v8sf rhs) const {
        return val + rhs.val;
    }
    v8sf operator - (const v8sf rhs) const {
        return val - rhs.val;
    }
    v8sf operator * (const v8sf rhs) const {
        return val * rhs.val;
    }
    v8sf operator / (const v8sf rhs) const {
        return v8sf(_mm256_div_ps(val, rhs.val));
    }
    
    static inline v8sf max(const v8sf a, const v8sf b) {
        return v8sf(_mm256_max_ps(a.val, b.val));
    }
    static inline v8sf min(const v8sf a, const v8sf b) {
        return v8sf(_mm256_min_ps(a.val, b.val));
    }
    
    v8sf operator += (const v8sf rhs) {
        val = val + rhs.val;
        return (*this);
    }
    v8sf operator -= (const v8sf rhs) {
        val = val - rhs.val;
        return (*this);
    }
    v8sf operator *= (const v8sf rhs) {
        val = val * rhs.val;
        return (*this);
    }
    v8sf operator /= (const v8sf rhs) {
        val = _mm256_div_ps(val, rhs.val);
        return (*this);
    }
    
    v8sf operator & (const v8sf rhs) {
        return v8sf(_mm256_and_ps(val, rhs.val));
    }
    v8sf operator != (const v8sf rhs) {
        return v8sf(_mm256_cmp_ps(val, rhs.val, _CMP_NEQ_UQ));
    }
    v8sf operator < (const v8sf rhs) {
        return v8sf(_mm256_cmp_ps(val, rhs.val, _CMP_LT_OS));
    }
    v8sf operator <= (const v8sf rhs) {
        return v8sf(_mm256_cmp_ps(val, rhs.val, _CMP_LE_OS));
    }
    v8sf operator > (const v8sf rhs) {
        return v8sf(_mm256_cmp_ps(val, rhs.val, _CMP_GT_OS));
    }
    v8sf operator >= (const v8sf rhs) {
        return v8sf(_mm256_cmp_ps(val, rhs.val, _CMP_GE_OS));
    }
    
    static inline v8sf rcp_0th(const v8sf rhs) {
        return _mm256_rcp_ps(rhs.val);
    }
    static inline v8sf rcp_1st(const v8sf rhs) {
        v8sf x0 = _mm256_rcp_ps(rhs.val);
        v8sf h  = v8sf(1.) - rhs * x0;
        return x0 + h * x0;
    }
    
    static inline v8sf sqrt(const v8sf rhs) {
        return v8sf(_mm256_sqrt_ps(rhs.val));
    }
    static inline v8sf rsqrt_0th(const v8sf rhs) {
        return v8sf(_mm256_rsqrt_ps(rhs.val));
    }
    static inline v8sf rsqrt_1st(const v8sf rhs) {
        v8sf x0 = v8sf(_mm256_rsqrt_ps(rhs.val));
        v8sf h  = v8sf(1.) - rhs * x0 * x0;
        return x0 + v8sf(0.5) * h * x0;
    }
    static inline v8sf rsqrt_1st_phantom(const v8sf rhs) {
        v8sf x0 = v8sf(_mm256_rsqrt_ps(rhs.val));
        return v8sf(x0 * (rhs * x0 * x0 - v8sf(3.)));
    }
    static inline v8sf hadd(v8sf x0, v8sf x1) {
        return _mm256_hadd_ps(x0, x1);
    }
    
    inline void store(float *p) const {
        _mm256_store_ps(p, val);
    }
    inline void load(float const *p) {
        val = _mm256_load_ps(p);
    }

    void cvtpd2ps(v4df & x0, v4df & x1);
    
    inline void extractf128(v4sf & x0, v4sf & x1) {
        x0 = _mm256_extractf128_ps(val, 0);
        x1 = _mm256_extractf128_ps(val, 1);
    }
    
    void print(FILE * fp = stdout,
               const char * fmt = "%+e %+e %+e %+e\n%+e %+e %+e %+e\n") const {
        float a[8] __attribute__((aligned(32)));
        store(a);
        fprintf(fp, fmt, a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]);
    }
    
    static inline v8sf shuffle0(v8sf rhs) {
        return _mm256_permute_ps(rhs, 0x00);
    }
    
    static inline v8sf shuffle1(v8sf rhs) {
        return _mm256_permute_ps(rhs, 0x55);
    }
    
    static inline v8sf shuffle2(v8sf rhs) {
        return _mm256_permute_ps(rhs, 0xaa);
    }
    
    static inline v8sf shuffle3(v8sf rhs) {
        return _mm256_permute_ps(rhs, 0xff);
    }

#ifndef PARALLEL_I2J8
    static inline v4sf reduce(v8sf rhs) {
        v4sf x0 = _mm256_extractf128_ps(rhs, 0);
        v4sf x1 = _mm256_extractf128_ps(rhs, 1);
        return v4sf(x0 + x1);
    }
#else
    static inline v4sf reduce(v8sf rhs) {
        v4sf x0 = _mm256_extractf128_ps(rhs, 0);
        v4sf x1 = _mm256_extractf128_ps(rhs, 1);
        return v4sf::hadd(x0, x1);
    }
#endif
    
};
struct v2df {
    typedef double _v2df __attribute__((vector_size(16)))
        __attribute__((aligned(16)));
    _v2df val;
    
    static inline int getVectorLength() {
        return 2;
    }
    
    v2df(const v2df & rhs) : val(rhs.val) {}
    v2df operator = (const v2df rhs) {
        val = rhs.val;
        return (*this);
    }
    v2df() : val(_mm_setzero_pd()) {}
    v2df(const double x) : val(_mm_set1_pd(x)) {}
    v2df(const double x0, const double y0)
        : val(_mm_set_pd(y0, x0)) {}
    v2df(const _v2df _val) : val(_val) {}
    operator _v2df() {return val;}
    
    v2df operator + (const v2df rhs) const {
        return v2df(val + rhs.val);
    }
    v2df operator - (const v2df rhs) const {
        return v2df(val - rhs.val);
    }
    v2df operator * (const v2df rhs) const {
        return v2df(val * rhs.val);
    }
    v2df operator / (const v2df rhs) const {
        return v2df(_mm_div_pd(val, rhs.val));
    }

    static v2df max(const v2df a, const v2df b) {
        return v2df(_mm_max_pd(a.val, b.val));
    }
    static v2df min(const v2df a, const v2df b) {
        return v2df(_mm_min_pd(a.val, b.val));
    }

    v2df operator += (const v2df rhs) {
        //val = _mm_add_pd(val, rhs.val);
        val = val + rhs.val;
        return (*this);
    }
    v2df operator -= (const v2df rhs) {
        //val = _mm_sub_pd(val, rhs.val);
        val = val - rhs.val;
        return (*this);
    }
    v2df operator *= (const v2df rhs) {
        //val = _mm_mul_pd(val, rhs.val);
        val = val * rhs.val;
        return (*this);
    }
    v2df operator /= (const v2df rhs) {
        val = _mm_div_pd(val, rhs.val);
        return (*this);
    }
    
    v2df operator & (const v2df rhs) {
        return v2df(_mm_and_pd(val, rhs.val));
    }
    v2df operator != (const v2df rhs) {
        return v2df(_mm_cmp_pd(val, rhs.val, _CMP_NEQ_UQ));
    }
    v2df operator < (const v2df rhs) {
        return v2df(_mm_cmp_pd(val, rhs.val, _CMP_LT_OS));
    }

    static v2df sqrt(const v2df rhs) {
        return v2df(_mm_sqrt_pd(rhs.val));
    }
    static v2df hadd(v2df x0, v2df x1) {
        return _mm_hadd_pd(x0, x1);
    }
    
    inline void store(double *p) const {
        _mm_store_pd(p, val);
    }
    inline void load(double const *p) {
        val = _mm_load_pd(p);
    }
    v4sf cvtpd2ps() {
        return _mm_cvtpd_ps(val);
    }
    
    void print(FILE * fp = stdout,
               const char * fmt = "%+e %+e\n") const {
        double a[2] __attribute__((aligned(16)));
        store(a);
        fprintf(fp, fmt, a[0], a[1]);
    }     
};


struct v4df {
    typedef double _v4df __attribute__((vector_size(32)))
        __attribute__((aligned(32)));
    _v4df val;
    
    static inline int getVectorLength() {
        return 4;
    }
    
    v4df(const v4df & rhs) : val(rhs.val) {}
    v4df operator = (const v4df rhs) {
        val = rhs.val;
        return (*this);
    }
    v4df() : val(_mm256_setzero_pd()) {}
    v4df(const double x) : val(_mm256_set1_pd(x)) {}
    v4df(const double x0, const double y0, const double z0, const double w0)
        : val(_mm256_set_pd(w0, z0, y0, x0)) {}
    v4df(const _v4df _val) : val(_val) {}
    operator _v4df() {return val;}
    v4df(const v2df rhs) {
        val = _mm256_broadcast_pd((const __m128d*)&rhs);
    }
    v4df(v2df rhs0, v2df rhs1) {
        v4df temp(rhs1);
        val = _mm256_insertf128_pd(temp, rhs0, 0);
    }
    
    v4df operator + (const v4df rhs) const {
        return v4df(val + rhs.val);
    }
    v4df operator - (const v4df rhs) const {
        return v4df(val - rhs.val);
    }
    v4df operator * (const v4df rhs) const {
        return v4df(val * rhs.val);
    }
    v4df operator / (const v4df rhs) const {
        return v4df(_mm256_div_pd(val, rhs.val));
    }

    static v4df max(const v4df a, const v4df b) {
        return v4df(_mm256_max_pd(a.val, b.val));
    }
    static v4df min(const v4df a, const v4df b) {
        return v4df(_mm256_min_pd(a.val, b.val));
    }
    
    v4df operator += (const v4df rhs) {
        val = val + rhs.val;
        return (*this);
    }
    v4df operator -= (const v4df rhs) {
        val = val - rhs.val;
        return (*this);
    }
    v4df operator *= (const v4df rhs) {
        val = val * rhs.val;
        return (*this);
    }
    v4df operator /= (const v4df rhs) {
        val = v4df(_mm256_div_pd(val, rhs.val));
        return (*this);
    }

    v4df operator & (const v4df rhs) {
        return v4df(_mm256_and_pd(val, rhs.val));
    }
    v4df operator != (const v4df rhs) {
        return v4df(_mm256_cmp_pd(val, rhs.val, _CMP_NEQ_UQ));
    }
    v4df operator < (const v4df rhs) {
        return v4df(_mm256_cmp_pd(val, rhs.val, _CMP_LT_OS));
    }
    
    inline void store(double *p) const {
        _mm256_store_pd(p, val);
    }
    inline void load(double const *p) {
        val = _mm256_load_pd(p);
    }
    v4sf cvtpd2ps() {
        return _mm256_cvtpd_ps(val);
    }
    void extractf128(v2df & x0, v2df & x1) {
        x0 = _mm256_extractf128_pd(val, 0);
        x1 = _mm256_extractf128_pd(val, 1);
    }

    static v4df sqrt(const v4df rhs) {
        return v4df(_mm256_sqrt_pd(rhs.val));
    }
    static v4df rsqrt_0th(v4df rhs) {
        v4df x0 = (v4sf::rsqrt_0th(rhs.cvtpd2ps())).cvtps2pd();
        return x0;
    }
    static v4df rsqrt_1st(v4df rhs) {
        v4df x0 = (v4sf::rsqrt_0th(rhs.cvtpd2ps())).cvtps2pd();
        v4df h  = v4df(1.) - rhs * x0 * x0;
        return x0 + v4df(0.5) * h * x0;
    }
    static v4df rsqrt_2nd(v4df rhs) {
        v4df x0 = (v4sf::rsqrt_0th(rhs.cvtpd2ps())).cvtps2pd();
        v4df h  = v4df(1.) - rhs * x0 * x0;
        return x0 + (v4df(0.5) + v4df(0.375) * h) * h * x0;
    }
    static v4df rsqrt_3rd(v4df rhs) {
        v4df x0 = (v4sf::rsqrt_0th(rhs.cvtpd2ps())).cvtps2pd();
        v4df h  = v4df(1.) - rhs * x0 * x0;
        return x0 + (v4df(0.5) + (v4df(0.375) + v4df(0.3125) * h) * h) * h * x0;
    }
    static v4df rsqrt_4th(v4df rhs) {
        v4df x0 = (v4sf::rsqrt_0th(rhs.cvtpd2ps())).cvtps2pd();
        v4df h  = v4df(1.) - rhs * x0 * x0;
        return x0
            + (v4df(0.5) + (v4df(0.375) + (v4df(0.3125) + v4df(0.2734375) * h) * h) * h) * h * x0;
    }
    static v4df rsqrt_1st_phantom(v4df rhs) {
        v4sf x1 = v4sf::rsqrt_1st_phantom(rhs.cvtpd2ps());
        return x1.cvtps2pd();
    }

    static v4df fabs(v4df rhs) {
        //v4df signmask(-0.0d);
        v4df signmask(-0.0);
        return _mm256_andnot_pd(signmask, rhs);
    }
    
    static v4df hadd(v4df x0, v4df x1) {
        return _mm256_hadd_pd(x0, x1);
    }
    
    void print(FILE * fp = stdout,
               const char * fmt = "%+e %+e %+e %+e\n") const {
        double a[4] __attribute__((aligned(32)));
        store(a);
        fprintf(fp, fmt, a[0], a[1], a[2], a[3]);
    }

    static inline v2df reduce(v4df rhs) {
        v2df x0 = _mm256_extractf128_pd(rhs, 0);
        v2df x1 = _mm256_extractf128_pd(rhs, 1);
        return v2df(x0 + x1);
    }
};

void v8sf::cvtpd2ps(v4df & x0, v4df & x1) {
    x0 = _mm256_cvtps_pd(_mm256_extractf128_ps(val, 0));
    x1 = _mm256_cvtps_pd(_mm256_extractf128_ps(val, 1));
}

v4df v4sf::cvtps2pd() {
    return _mm256_cvtps_pd(val);
}

#ifdef __AVX512F__

struct v16sf {
    typedef float _v16sf __attribute__((vector_size(64)))
        __attribute__((aligned(64)));
    _v16sf val;
    
    static inline int getVectorLength() {
        return 16;
    }
    
    v16sf(const v16sf & rhs) : val(rhs.val) {}
    v16sf operator = (const v16sf rhs) {
        val = rhs.val;
        return (*this);
    }
    v16sf() : val(_mm512_setzero_ps()) {}
    v16sf(const float x) : val(_mm512_set1_ps(x)) {}
    v16sf(const float x0, const float y0, const float z0, const float w0,
          const float x1, const float y1, const float z1, const float w1,
          const float x2, const float y2, const float z2, const float w2,
          const float x3, const float y3, const float z3, const float w3)
        : val(_mm512_set_ps(w3, z3, y3, x3,
                            w2, z2, y2, x2,
                            w1, z1, y1, x1,
                            w0, z0, y0, x0)) {}
    v16sf(const _v16sf _val) : val(_val) {}
    operator _v16sf() {return val;}
    v16sf(v4sf rhs) {
        val = _mm512_broadcast_f32x4(rhs);
    }
#ifndef __AVX512DQ__
    // terrible implementation (A. Tanikawa);
    v16sf(v4sf rhs0, v4sf rhs1, v4sf rhs2, v4sf rhs3) {
        val = v16sf(rhs0);
        val = _mm512_insertf32x4(val, rhs1, 1);
        val = _mm512_insertf32x4(val, rhs2, 2);
        val = _mm512_insertf32x4(val, rhs3, 3);
    }
    v16sf(v8sf rhs0, v8sf rhs1) {
        v4sf x0 = _mm256_extractf128_ps(rhs0, 0);
        v4sf x1 = _mm256_extractf128_ps(rhs0, 1);
        v4sf x2 = _mm256_extractf128_ps(rhs1, 0);
        v4sf x3 = _mm256_extractf128_ps(rhs1, 1);
        val = v16sf(x0, x1, x2, x3);
    }
    v16sf(v8sf rhs) {
        v4sf x0 = _mm256_extractf128_ps(rhs, 0);
        v4sf x1 = _mm256_extractf128_ps(rhs, 1);
        v4sf x2 = x0;
        v4sf x3 = x1;
        val = v16sf(x0, x1, x2, x3);
    }
#else
    v16sf(v8sf rhs) {
        val = _mm512_broadcast_f32x8(rhs);
    }
    v16sf(v8sf rhs0, v8sf rhs1) {
        val = _mm512_insertf32x8(v16sf(rhs0), rhs1, 1);
    }
    v16sf(v4sf rhs0, v4sf rhs1, v4sf rhs2, v4sf rhs3) {
        val = v16sf(v8sf(rhs0, rhs1), v8sf(rhs2, rhs3));
    }
#endif

    static inline v16sf max(const v16sf a, const v16sf b) {
        return v16sf(_mm512_max_ps(a.val, b.val));
    }
    static inline v16sf min(const v16sf a, const v16sf b) {
        return v16sf(_mm512_min_ps(a.val, b.val));
    }
        
    v16sf operator + (const v16sf rhs) const {
        return val + rhs.val;
    }
    v16sf operator - (const v16sf rhs) const {
        return val - rhs.val;
    }
    v16sf operator * (const v16sf rhs) const {
        return val * rhs.val;
    }
    v16sf operator / (const v16sf rhs) const {
        return v16sf(_mm512_div_ps(val, rhs.val));
    }
    
    v16sf operator += (const v16sf rhs) {
        val = val + rhs.val;
        return (*this);
    }
    v16sf operator -= (const v16sf rhs) {
        val = val - rhs.val;
        return (*this);
    }
    v16sf operator *= (const v16sf rhs) {
        val = val * rhs.val;
        return (*this);
    }
    v16sf operator /= (const v16sf rhs) {
        val = _mm512_div_ps(val, rhs.val);
        return (*this);
    }
    
    static inline v16sf rcp_0th(const v16sf rhs) {
        return _mm512_rcp14_ps(rhs.val);
    }
    static inline v16sf rcp_1st(const v16sf rhs) {
        /*
          v16sf x0 = _mm512_rcp_ps(rhs.val);
          v16sf h  = v16sf(1.) - rhs * x0;
          return x0 + h * x0;
        */
        return _mm512_rcp28_ps(rhs.val);
    }
    
    static inline v16sf sqrt(const v16sf rhs) {
        return _mm512_sqrt_ps(rhs.val);
    }
    static inline v16sf rsqrt_0th(const v16sf rhs) {
        return _mm512_rsqrt14_ps(rhs.val);
    }
    static inline v16sf rsqrt_1st(const v16sf rhs) {
        v16sf x0 = _mm512_rsqrt14_ps(rhs.val);
        v16sf h  = v16sf(1.) - rhs * x0 * x0;
        return x0 + v16sf(0.5) * h * x0;
    }
    static inline v16sf rsqrt_1st_phantom(const v16sf rhs) {
        v16sf x0 = _mm512_rsqrt14_ps(rhs.val);
        return x0 * (rhs * x0 * x0 - v16sf(3.));
    }
    
    inline void store(float *p) const {
        _mm512_store_ps(p, val);
    }
    inline void load(float const *p) {
        val = _mm512_load_ps(p);
    }

    void print(FILE * fp = stdout,
               const char * fmt =
               "%+e %+e %+e %+e\n%+e %+e %+e %+e\n%+e %+e %+e %+e\n%+e %+e %+e %+e\n")
        const {
        float a[16] __attribute__((aligned(64)));
        store(a);
        fprintf(fp, fmt,
                a[0],  a[1],  a[2],  a[3],
                a[4],  a[5],  a[6],  a[7],
                a[8],  a[9],  a[10], a[11],
                a[12], a[13], a[14], a[15]);
    }
    
    static inline v16sf shuffle0(v16sf rhs) {
        return _mm512_permute_ps(rhs, 0x00);
    }
    
    static inline v16sf shuffle1(v16sf rhs) {
        return _mm512_permute_ps(rhs, 0x55);
    }
    
    static inline v16sf shuffle2(v16sf rhs) {
        return _mm512_permute_ps(rhs, 0xaa);
    }
    
    static inline v16sf shuffle3(v16sf rhs) {
        return _mm512_permute_ps(rhs, 0xff);
    }

#ifndef PARALLEL_I2J8
    static inline v4sf reduce16to4(v16sf rhs) {
        v4sf x0 = _mm512_extractf32x4_ps(rhs, 0);
        v4sf x1 = _mm512_extractf32x4_ps(rhs, 1);
        v4sf x2 = _mm512_extractf32x4_ps(rhs, 2);
        v4sf x3 = _mm512_extractf32x4_ps(rhs, 3);
        return (x0 + x1 + x2 + x3);
    }
#else
    static inline v4sf reduce16to4(v16sf rhs) {
        v4sf x0 = _mm512_extractf32x4_ps(rhs, 0);
        v4sf x1 = _mm512_extractf32x4_ps(rhs, 1);
        v4sf x2 = _mm512_extractf32x4_ps(rhs, 2);
        v4sf x3 = _mm512_extractf32x4_ps(rhs, 3);
        return v4sf::hadd((x0 + x1), (x2 + x3));
    }
#endif

    static inline v8sf reduce(v16sf rhs) {
        v4sf x0 = _mm512_extractf32x4_ps(rhs, 0);
        v4sf x1 = _mm512_extractf32x4_ps(rhs, 1);
        v4sf x2 = _mm512_extractf32x4_ps(rhs, 2);
        v4sf x3 = _mm512_extractf32x4_ps(rhs, 3);
        v4sf x02 = x0 + x2;
        v4sf x13 = x1 + x3;
        return (v8sf(x02, x13));
    }

    /*
    static inline v8sf reduce16to08(v16sf rhs) {
        v4sf x0 = _mm512_extractf32x4_ps(rhs, 0);
        v4sf x1 = _mm512_extractf32x4_ps(rhs, 1);
        v4sf x2 = _mm512_extractf32x4_ps(rhs, 2);
        v4sf x3 = _mm512_extractf32x4_ps(rhs, 3);
        v4sf x01 = x0 + x1;
        v4sf x23 = x2 + x3;
        return (v8sf(x01, x23));
    }
    */
};


struct v8df {
    typedef double _v8df __attribute__((vector_size(64)))
        __attribute__((aligned(64)));
    _v8df val;
    
    static inline int getVectorLength() {
        return 8;
    }
    
    v8df(const v8df & rhs) : val(rhs.val) {}
    v8df operator = (const v8df rhs) {
        val = rhs.val;
        return (*this);
    }
    v8df() : val(_mm512_setzero_pd()) {}
    v8df(const double x) : val(_mm512_set1_pd(x)) {}
    v8df(const double x1, const double y1, const double z1, const double w1,
         const double x0, const double y0, const double z0, const double w0)
        : val(_mm512_set_pd(w1, z1, y1, x1, w0, z0, y0, x0)) {}
    v8df(const _v8df _val) : val(_val) {}
    operator _v8df() {return val;}
    v8df(v4df rhs0) {
        val = _mm512_broadcast_f64x4(rhs0);
    }
    v8df(v4df rhs0, v4df rhs1) {
        v8df temp(rhs0);
        val = _mm512_insertf64x4(temp, rhs1, 1);
    }
    
    v8df operator + (const v8df rhs) const {
        return v8df(val + rhs.val);
    }
    v8df operator - (const v8df rhs) const {
        return v8df(val - rhs.val);
    }
    v8df operator * (const v8df rhs) const {
        return v8df(val * rhs.val);
    }
    v8df operator / (const v8df rhs) const {
        return v8df(_mm512_div_pd(val, rhs.val));
    }
    v8df operator += (const v8df rhs) {
        val = val + rhs.val;
        return (*this);
    }
    v8df operator -= (const v8df rhs) {
        val = val - rhs.val;
        return (*this);
    }
    v8df operator *= (const v8df rhs) {
        val = val * rhs.val;
        return (*this);
    }
    v8df operator /= (const v8df rhs) {
        val = v8df(_mm512_div_pd(val, rhs.val));
        return (*this);
    }
    
    inline void store(double *p) const {
        _mm512_store_pd(p, val);
    }
    inline void load(double const *p) {
        val = _mm512_load_pd(p);
    }

    // terrible implementation (A. Tanikawa)
    static inline v8df hadd_pd(v8df rhs0, v8df rhs1) {
        v4df x0 = _mm512_extractf64x4_pd(rhs0, 0);;
        v4df x1 = _mm512_extractf64x4_pd(rhs0, 1);;
        v4df x2 = _mm512_extractf64x4_pd(rhs1, 0);;
        v4df x3 = _mm512_extractf64x4_pd(rhs1, 1);;
        v4df x01 = _mm256_hadd_pd(x0, x1);
        v4df x23 = _mm256_hadd_pd(x2, x3);
        return v8df(x01, x23);
    }
    
    void print(FILE * fp = stdout,
               const char * fmt = "%+e %+e %+e %+e\n%+e %+e %+e %+e\n") const {
        double a[8] __attribute__((aligned(64)));
        store(a);
        fprintf(fp, fmt,
                a[0], a[1], a[2], a[3],
                a[4], a[5], a[6], a[7]);
    }    
};

/*
  #ifndef __AVX512DQ__
  inline v8sf _mm512_extractf32x8_ps(v16sf rhs, int imm8) {
  v4sf x0, x1;
  if(imm8 == 0) {
  x0 = _mm512_extractf32x4_ps(rhs, 0);
  x1 = _mm512_extractf32x4_ps(rhs, 1);
  } else {
  x0 = _mm512_extractf32x4_ps(rhs, 2);
  x1 = _mm512_extractf32x4_ps(rhs, 3);
  }
  return v8sf(x0, x1);
  }
  #endif
*/

#endif
