#pragma once

#include <immintrin.h>

struct v4sf;
struct v8sf;
struct v2df;
struct v4df;

struct v4sf {
    typedef float _v4sf __attribute__((vector_size(16))) __attribute__((aligned(16)));
    _v4sf val;

    static int getVectorLength() {
        return 4;
    }

    v4sf(const v4sf & rhs) : val(rhs.val) {}
    v4sf operator = (const v4sf rhs) {
        val = rhs.val;
        return (*this);
    }

    v4sf() : val(_mm_setzero_ps()) {}
    v4sf(const float x) : val(_mm_set_ps(x, x, x, x)) {}
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
#ifdef __AVX2__
    static v4sf madd(const v4sf c, const v4sf a, const v4sf b) {
        return v4sf(_mm_fmadd_ps(a.val, b.val, c.val));
    }
    static v4sf nmadd(const v4sf c, const v4sf a, const v4sf b) {
        return v4sf(_mm_fnmadd_ps(a.val, b.val, c.val));
    }
#endif

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

    static v4sf rcp_0th(const v4sf rhs) {
        return _mm_rcp_ps(rhs.val);
    }
    static v4sf rcp_1st(const v4sf rhs) {
        v4sf x0 = _mm_rcp_ps(rhs.val);
        v4sf h  = v4sf(1.) - rhs * x0;
        return x0 + h * x0;
    }

    static v4sf sqrt(const v4sf rhs) {
        return v4sf(_mm_sqrt_ps(rhs.val));
    }
    static v4sf rsqrt_0th(const v4sf rhs) {
        return v4sf(_mm_rsqrt_ps(rhs.val));
    }
    static v4sf rsqrt_1st(const v4sf rhs) {
        v4sf x0 = v4sf(_mm_rsqrt_ps(rhs.val));
        v4sf h  = v4sf(1.) - rhs * x0 * x0;
        return x0 + v4sf(0.5) * h * x0;
    }
    static v4sf rsqrt_1st_phantom(const v4sf rhs) {
        v4sf x0 = v4sf(_mm_rsqrt_ps(rhs.val));
        return v4sf(x0 * (rhs * x0 * x0 - v4sf(3.)));
    }

    void store(float *p) const {
        _mm_store_ps(p, val);
    }
    void load(float const *p) {
        val = _mm_load_ps(p);
    }

    void print(FILE * fp = stdout,
               const char * fmt = "%+e %+e %+e %+e\n") const {
        int vl = getVectorLength();
        float a[vl];
        store(a);
        fprintf(fp, fmt, a[0], a[1], a[2], a[3]);
    }    
};

struct v8sf {
    typedef float _v8sf __attribute__((vector_size(32))) __attribute__((aligned(32)));
    _v8sf val;

    static int getVectorLength() {
        return 8;
    }

    v8sf(const v8sf & rhs) : val(rhs.val) {}
    v8sf operator = (const v8sf rhs) {
        val = rhs.val;
        return (*this);
    }

    v8sf() : val(_mm256_setzero_ps()) {}
    v8sf(const float x) : val(_mm256_set_ps(x, x, x, x, x, x, x, x)) {}
    v8sf(const float x0, const float y0, const float z0, const float w0,
         const float x1, const float y1, const float z1, const float w1)
        : val(_mm256_set_ps(w1, z1, y1, x1, w0, z0, y0, x0)) {}

    v8sf(const _v8sf _val) : val(_val) {}
    operator _v8sf() {return val;}

    v8sf(const v4sf rhs) {
        float buf[8];
        rhs.store(&buf[0]);
        rhs.store(&buf[4]);
        load(buf);
    }

    v8sf(const v4sf x0, const v4sf x1) {
        float buf[8];
        x0.store(&buf[0]);
        x1.store(&buf[4]);
        load(buf);
    }

    v8sf operator + (const v8sf rhs) const {
        //return v8sf(_mm256_add_ps(val, rhs.val));
        return val + rhs.val;
    }
    v8sf operator - (const v8sf rhs) const {
        //return v8sf(_mm256_sub_ps(val, rhs.val));
        return val - rhs.val;
    }
    v8sf operator * (const v8sf rhs) const {
        //return v8sf(_mm256_mul_ps(val, rhs.val));
        return val * rhs.val;
    }
    v8sf operator / (const v8sf rhs) const {
        return v8sf(_mm256_div_ps(val, rhs.val));
    }
#ifdef __AVX2__
    static v8sf madd(const v8sf c, const v8sf a, const v8sf b) {
        return v8sf(_mm256_fmadd_ps(a.val, b.val, c.val));
    }
    static v8sf nmadd(const v8sf c, const v8sf a, const v8sf b) {
        return v8sf(_mm256_fnmadd_ps(a.val, b.val, c.val));
    }
#endif

    static v8sf max(const v8sf a, const v8sf b) {
        return v8sf(_mm256_max_ps(a.val, b.val));
    }
    static v8sf min(const v8sf a, const v8sf b) {
        return v8sf(_mm256_min_ps(a.val, b.val));
    }

    v8sf operator += (const v8sf rhs) {
        //val = _mm256_add_ps(val, rhs.val);
        val = val + rhs.val;
        return (*this);
    }
    v8sf operator -= (const v8sf rhs) {
        //val = _mm256_sub_ps(val, rhs.val);
        val = val - rhs.val;
        return (*this);
    }
    v8sf operator *= (const v8sf rhs) {
        //val = _mm256_mul_ps(val, rhs.val);
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

    static v8sf rcp_0th(const v8sf rhs) {
        return _mm256_rcp_ps(rhs.val);
    }
    static v8sf rcp_1st(const v8sf rhs) {
        v8sf x0 = _mm256_rcp_ps(rhs.val);
        v8sf h  = v8sf(1.) - rhs * x0;
        return x0 + h * x0;
    }

    static v8sf sqrt(const v8sf rhs) {
        return v8sf(_mm256_sqrt_ps(rhs.val));
    }
    static v8sf rsqrt_0th(const v8sf rhs) {
        return v8sf(_mm256_rsqrt_ps(rhs.val));
    }
    static v8sf rsqrt_1st(const v8sf rhs) {
        v8sf x0 = v8sf(_mm256_rsqrt_ps(rhs.val));
        v8sf h  = v8sf(1.) - rhs * x0 * x0;
        return x0 + v8sf(0.5) * h * x0;
    }
    static v8sf rsqrt_1st_phantom(const v8sf rhs) {
        v8sf x0 = v8sf(_mm256_rsqrt_ps(rhs.val));
        return v8sf(x0 * (rhs * x0 * x0 - v8sf(3.)));
    }
    static v8sf hadd(v8sf x0, v8sf x1) {
        return _mm256_hadd_ps(x0, x1);
    }

    void store(float *p) const {
        _mm256_store_ps(p, val);
    }
    void load(float const *p) {
        val = _mm256_load_ps(p);
    }
    void cvtpd2ps(v4df & x0, v4df & x1);
    void extractf128(v4sf & x0, v4sf & x1) {
        x0 = _mm256_extractf128_ps(val, 0);
        x1 = _mm256_extractf128_ps(val, 1);
    }

    void print(FILE * fp = stdout,
               const char * fmt = "%+e %+e %+e %+e\n%+e %+e %+e %+e\n\n") const {
        int vl = getVectorLength();
        float a[vl];
        store(a);
        fprintf(fp, fmt, a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]);
    }

    static v8sf shuffle0(v8sf rhs) {
        return _mm256_permute_ps(rhs, 0x00);
    }

    static v8sf shuffle1(v8sf rhs) {
        return _mm256_permute_ps(rhs, 0x55);
    }

    static v8sf shuffle2(v8sf rhs) {
        return _mm256_permute_ps(rhs, 0xaa);
    }

    static v8sf shuffle3(v8sf rhs) {
        return _mm256_permute_ps(rhs, 0xff);
    }

    static v4sf reduce(const v8sf rhs) {
        int nv = getVectorLength();
        float buf[nv];
        rhs.store(buf);
        v4sf x0(buf[0], buf[1], buf[2], buf[3]);
        v4sf x1(buf[4], buf[5], buf[6], buf[7]);
        return v4sf(x0 + x1);
    }

};

struct v2df {
    typedef double _v2df __attribute__((vector_size(16))) __attribute__((aligned(16)));
    _v2df val;

    static int getVectorLength() {
        return 2;
    }

    v2df(const v2df & rhs) : val(rhs.val) {}
    v2df operator = (const v2df rhs) {
        val = rhs.val;
        return (*this);
    }

    v2df() : val(_mm_setzero_pd()) {}
    v2df(const double x) : val(_mm_set_pd(x, x)) {}
    v2df(const double x, const double y)
        : val(_mm_set_pd(y, x)) {}

    v2df(const _v2df _val) : val(_val) {}
    operator _v2df() {return val;}

    v2df operator + (const v2df rhs) const {
        //return v2df(_mm_add_pd(val, rhs.val));
        return val + rhs.val;
    }
    v2df operator - (const v2df rhs) const {
        //return v2df(_mm_sub_pd(val, rhs.val));
        return val - rhs.val;
    }
    v2df operator * (const v2df rhs) const {
        //return v2df(_mm_mul_pd(val, rhs.val));
        return val * rhs.val;
    }
    v2df operator / (const v2df rhs) const {
        return v2df(_mm_div_pd(val, rhs.val));
    }
#ifdef __AVX2__
    static v2df madd(const v2df c, const v2df a, const v2df b) {
        return v2df(_mm_fmadd_pd(a.val, b.val, c.val));
    }
    static v2df nmadd(const v2df c, const v2df a, const v2df b) {
        return v2df(_mm_fnmadd_pd(a.val, b.val, c.val));
    }
#endif
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

    void store(double *p) const {
        _mm_store_pd(p, val);
    }
    void storel(double *p) const {
        _mm_storel_pd(p, val);
    }
    void storeh(double *p) const {
        _mm_storeh_pd(p, val);
    }
    void load(double const *p) {
        val = _mm_load_pd(p);
    }
    v4sf cvtpd2ps() {
        return _mm_cvtpd_ps(val);
    }

    void print(FILE * fp = stdout,
               const char * fmt = "%+e %+e\n") const {
        int vl = getVectorLength();
        double a[vl];
        _mm_store_pd(a, val);
        fprintf(fp, fmt, a[0], a[1]);
    }    
};

struct v4df {
    typedef double _v4df __attribute__((vector_size(32))) __attribute__((aligned(32)));
    _v4df val;

    static int getVectorLength() {
        return 4;
    }

    v4df(const v4df & rhs) : val(rhs.val) {}
    v4df operator = (const v4df rhs) {
        val = rhs.val;
        return (*this);
    }

    v4df() : val(_mm256_setzero_pd()) {}
    v4df(const double x) : val(_mm256_set_pd(x, x, x, x)) {}
    v4df(const double x, const double y, const double z, const double w)
        : val(_mm256_set_pd(w, z, y, x)) {}

    v4df(const _v4df _val) : val(_val) {}
    operator _v4df() {return val;}

    //~v4df() {}

    v4df operator + (const v4df rhs) const {
        //return v4df(_mm256_add_pd(val, rhs.val));
        return val + rhs.val;
    }
    v4df operator - (const v4df rhs) const {
        //return v4df(_mm256_sub_pd(val, rhs.val));
        return val - rhs.val;
    }
    v4df operator * (const v4df rhs) const {
        //return v4df(_mm256_mul_pd(val, rhs.val));
        return val * rhs.val;
    }
    v4df operator / (const v4df rhs) const {
        return v4df(_mm256_div_pd(val, rhs.val));
    }

#ifdef __AVX2__
    static v4df madd(const v4df c, const v4df a, const v4df b) {
        return v4df(_mm256_fmadd_pd(a.val, b.val, c.val));
    }
    static v4df nmadd(const v4df c, const v4df a, const v4df b) {
        return v4df(_mm256_fnmadd_pd(a.val, b.val, c.val));
    }
#endif
    static v4df max(const v4df a, const v4df b) {
        return v4df(_mm256_max_pd(a.val, b.val));
    }
    static v4df min(const v4df a, const v4df b) {
        return v4df(_mm256_min_pd(a.val, b.val));
    }

    v4df operator += (const v4df rhs) {
        //val = _mm256_add_pd(val, rhs.val);
        val = val + rhs.val;
        return (*this);
    }
    v4df operator -= (const v4df rhs) {
        //val = _mm256_sub_pd(val, rhs.val);
        val = val - rhs.val;
        return (*this);
    }
    v4df operator *= (const v4df rhs) {
        //val = _mm256_mul_pd(val, rhs.val);
        val = val * rhs.val;
        return (*this);
    }
    v4df operator /= (const v4df rhs) {
        val = _mm256_div_pd(val, rhs.val);
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

    void store(double *p) const {
        _mm256_store_pd(p, val);
    }
    void load(double const *p) {
        val = _mm256_load_pd(p);
    }
    v4sf cvtpd2ps() {
        return _mm256_cvtpd_ps(val);
    }
    void extractf128(v2df & x0, v2df & x1) {
        x0 = _mm256_extractf128_pd(val, 0);
        x1 = _mm256_extractf128_pd(val, 1);
    }

    static v4df rcp_0th(v4df rhs) {
        v4df x0 = (v4sf::rcp_0th(rhs.cvtpd2ps())).cvtps2pd();
        return x0;
    }
    static v4df rcp_1st(v4df rhs) {
        v4df x0 = (v4sf::rcp_0th(rhs.cvtpd2ps())).cvtps2pd();
        v4df h  = v4df(1.) - rhs * x0;
        return x0 + h * x0;
    }
    static v4df rcp_4th(v4df rhs) {
        v4df x0 = (v4sf::rcp_0th(rhs.cvtpd2ps())).cvtps2pd();
        v4df h  = v4df(1.) - rhs * x0;
        return (v4df(1.) + h) * (v4df(1.) + h * h) * x0;
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
//        v4df signmask(-0.0d);
        v4df signmask(-0.0);
        return _mm256_andnot_pd(signmask, rhs);
    }

    static v4df hadd(v4df x0, v4df x1) {
        return _mm256_hadd_pd(x0, x1);
    }

    void print(FILE * fp = stdout,
               const char * fmt = "%+e %+e %+e %+e\n") const {
        int vl = getVectorLength();
        double a[vl];
        _mm256_store_pd(a, val);
        fprintf(fp, fmt, a[0], a[1], a[2], a[3]);
    }    
};

void v8sf::cvtpd2ps(v4df & x0, v4df & x1) {
    x0 = _mm256_cvtps_pd(_mm256_extractf128_ps(val, 0));
    x1 = _mm256_cvtps_pd(_mm256_extractf128_ps(val, 1));
}

v4df v4sf::cvtps2pd() {
    return _mm256_cvtps_pd(val);
}

