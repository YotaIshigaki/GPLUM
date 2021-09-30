#include <micvec.h>
//#include <immintrin.h>

struct v16sf {    
    typedef float _v16sf __attribute__((vector_size(64))) __attribute__((aligned(64)));
    //__m512 val;
    _v16sf val;

    static int getVectorLength() {
        return 16;
    }
    
    v16sf(const v16sf & rhs) : val(rhs.val) {}
    v16sf operator = (const v16sf rhs) {
        val = rhs.val;
        return (*this);
    }
    v16sf() : val(_mm512_setzero_ps()) {}
    v16sf(const float x) :
        val(_mm512_set_ps(x, x, x, x,
                          x, x, x, x,
                          x, x, x, x,
                          x, x, x, x)) {}
    v16sf(const float x, const float y) :
        val(_mm512_set_ps(y, y, y, y,
                          y, y, y, y,
                          x, x, x, x,
                          x, x, x, x)) {}
    v16sf(const float x, const float y, const float z, const float w) :
        val(_mm512_set_ps(w, w, w, w,
                          z, z, z, z,
                          y, y, y, y,
                          x, x, x, x)) {}
    v16sf(const float x0, const float y0, const float z0, const float w0,
          const float x1, const float y1, const float z1, const float w1) :
        val(_mm512_set_ps(w1, w1, z1, z1,
                          y1, y1, x1, x1,
                          w0, w0, z0, z0,
                          y0, y0, x0, x0)) {}
    v16sf(const float x0, const float y0, const float z0, const float w0,
          const float x1, const float y1, const float z1, const float w1,
          const float x2, const float y2, const float z2, const float w2,
          const float x3, const float y3, const float z3, const float w3) :
        val(_mm512_set_ps(w3, z3, y3, x3,
                          w2, z2, y2, x2,
                          w1, z1, y1, x1,
                          w0, z0, y0, x0)) {}
    v16sf(const __m512 _val) : val(_val) {}
    operator __m512() {return val;}

    v16sf operator + (const v16sf rhs) const {
        return _mm512_add_ps(val, rhs.val);
    }
    v16sf operator - (const v16sf rhs) const {
        return _mm512_sub_ps(val, rhs.val);
    }
    v16sf operator * (const v16sf rhs) const {
        return _mm512_mul_ps(val, rhs.val);
    }
    static v16sf madd(const v16sf c, const v16sf a, const v16sf b) {
        return v16sf(_mm512_fmadd_ps(a.val, b.val, c.val));
    }
    static v16sf nmadd(const v16sf c, const v16sf a, const v16sf b) {
        return v16sf(_mm512_fnmadd_ps(a.val, b.val, c.val));
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

    static v16sf rsqrt_0th(const v16sf rhs) {
        return v16sf(_mm512_rsqrt23_ps(rhs.val));
    }

    static v16sf permute(const v16sf rhs, const _MM_PERM_ENUM imm8) {
        return v16sf(_mm512_permute4f128_ps(rhs.val, imm8));
    }
    static v16sf swizzle(const v16sf rhs, const _MM_SWIZZLE_ENUM imm8) {
        return v16sf(_mm512_swizzle_ps(rhs.val, imm8));
    }
    static v16sf broadcast(float mt) {
        return v16sf(_mm512_extload_ps(&mt, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, 0));
    }
    static float mask_reduce_add(const __mmask16 imm, v16sf rhs) {
        return _mm512_mask_reduce_add_ps(imm, rhs);
    }

    void store(float *p) const {
        _mm512_store_ps(p, val);
    }
    void load(float const *p) {
        val = _mm512_load_ps(p);
    }

    void print(FILE * fp = stdout,
               const char * fmt = "%+e %+e %+e %+e\n%+e %+e %+e %+e\n%+e %+e %+e %+e\n%+e %+e %+e %+e\n\n") const {
        int vl = getVectorLength();
        float a[vl];
        store(a);
        fprintf(fp, fmt,
                a[0],  a[1] , a[2],  a[3],
                a[4],  a[5],  a[6],  a[7],
                a[8],  a[9],  a[10], a[11],
                a[12], a[13], a[14], a[15]);
    }    

};

