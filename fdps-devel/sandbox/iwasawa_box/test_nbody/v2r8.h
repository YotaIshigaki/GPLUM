#pragma once
#include <cstdio>
struct v2r8_mask{
	typedef __builtin_v2r8 _v2r8;
	_v2r8 val;

	v2r8_mask(const v2r8_mask &rhs) : val(rhs.val) {}
	v2r8_mask operator=(const v2r8_mask rhs){
		val = rhs.val;
		return (*this);
	}
	v2r8_mask(const _v2r8 _val) : val(_val) {}
};

struct v2r8{
	typedef __builtin_v2r8 _v2r8;
	_v2r8 val;

	v2r8(const v2r8 &rhs) : val(rhs.val) {}
	v2r8 operator=(const v2r8 rhs){
		val = rhs.val;
		return (*this);
	}

	v2r8() : val(__builtin_fj_setzero_v2r8()) {}
	v2r8(const double x) : val(__builtin_fj_set_v2r8(x, x)) {}
	v2r8(const double lo, const double hi) : val(__builtin_fj_set_v2r8(lo, hi)) {}
	v2r8(const _v2r8 _val) : val(_val) {}
	operator _v2r8() { return val; }

	v2r8 operator+(const v2r8 rhs) const {
		return v2r8(__builtin_fj_add_v2r8(val, rhs.val));
	}
	v2r8 operator-(const v2r8 rhs) const {
		return v2r8(__builtin_fj_sub_v2r8(val, rhs.val));
	}
	v2r8 operator*(const v2r8 rhs) const {
		return v2r8(__builtin_fj_mul_v2r8(val, rhs.val));
	}
	static v2r8 madd(const v2r8 a, const v2r8 b, const v2r8 c){
		return v2r8(__builtin_fj_madd_v2r8(a.val, b.val, c.val));
	}
	static v2r8 nmsub(const v2r8 a, const v2r8 b, const v2r8 c){
		return v2r8(__builtin_fj_nmsub_v2r8(a.val, b.val, c.val));
	}


	v2r8 operator+=(const v2r8 rhs){
		val = __builtin_fj_add_v2r8(val, rhs.val);
		return (*this);
	}
	v2r8 operator-=(const v2r8 rhs){
		val = __builtin_fj_sub_v2r8(val, rhs.val);
		return (*this);
	}
	v2r8 operator*=(const v2r8 rhs){
		val = __builtin_fj_mul_v2r8(val, rhs.val);
		return (*this);
	}

	v2r8 rsqrta() const {
		return v2r8(__builtin_fj_rsqrta_v2r8(val));
	}
	__attribute__ ((always_inline))
	v2r8 rsqrta_x3() const{
		const v2r8 x(val);
		const v2r8 y = v2r8(__builtin_fj_rsqrta_v2r8(val));
		const v2r8 y2 = y*y;
		const v2r8 h = v2r8(1.0) - x*y2;
#if 0
		const v2r8 p = (v2r8(3./8.) * h) * 
			// (v2r8(7./3.) - x*y2);
			v2r8::nmsub(x, y2, v2r8(7./3.)); // c - a*b
		return y + p*y; // 3-mul, 3-fma
#else
		const v2r8 p = v2r8(1.0) + h*(v2r8(0.5) + h*v2r8(3./8.));
		return p*y; // 2-mul, 3-fma
#endif
	}
	__attribute__ ((always_inline))
	v2r8 rsqrta_x8() const {
		const v2r8 x(val);
		const v2r8 y0 = v2r8(__builtin_fj_rsqrta_v2r8(val));
		const v2r8 half(0.5);
		const v2r8 xhalf = x * half;
		const v2r8 h0half = half - y0*(y0*xhalf);
		const v2r8 y1 = y0 + y0*h0half;           // x2
		const v2r8 h1half = half - y1*(y1*xhalf);
		const v2r8 y2 = y1 + y1*h1half;           // x4
		const v2r8 h2half = half - y2*(y2*xhalf);
		const v2r8 y3 = y2 + y2*h2half;           // x8
		return y3; // 4-mul, 6-fma
	}
	v2r8 rsqrta_x9() const {
		const v2r8 x50 = v2r8(0.50) * v2r8(val);
		const v2r8 x75 = v2r8(0.75) * v2r8(val);
		v2r8 y = v2r8(__builtin_fj_rsqrta_v2r8(val));
		y += y * ( (v2r8(0.5) - (y*y)*x50) * (v2r8(1.75) - (y*y)*x75) ); // 2-mul, 3-fma
		y += y * ( (v2r8(0.5) - (y*y)*x50) * (v2r8(1.75) - (y*y)*x75) ); // 2-mul, 3-fma
		return y; // 6-mul, 6-fma
	}
	__attribute__ ((always_inline))
	v2r8 rsqrta_x7() const{
		const v2r8 x(val);
		const v2r8 y0 = v2r8(__builtin_fj_rsqrta_v2r8(x.val));
		const v2r8 h  = v2r8(1.0) - x*(y0*y0);
		const v2r8 h2 = h*h;
		const v2r8 a  = v2r8(1./2.)    + h * v2r8(3./8.);
		const v2r8 b  = v2r8(5./16.)   + h * v2r8(35./128.);
		const v2r8 c  = v2r8(63./256.) + h * v2r8(231./1024.);
		const v2r8 p  = h*(a + h2*(b + h2*c));
		return y0 + p*y0; // 3-mul, 7-fma
	}

	void storel(double *p) const {
		__builtin_fj_storel_v2r8(p, val);
	}
	void storeh(double *p) const {
		__builtin_fj_storeh_v2r8(p, val);
	}
	// explicit store to generate std,s
	void store(double *p) const {
		__builtin_fj_store_v2r8(p, val);
	}
	// explicit load
	static v2r8 load(const double *p){
		return v2r8(__builtin_fj_load_v2r8(p));
	}
	v2r8 hadd() const {
		const v2r8 one(1.0);
		return v2r8(__builtin_fj_madd_sr1_v2r8(val, one.val, val));
	}

	void print(
			FILE *fp = stdout,
			const char *fmt = "lo = %f, hi = %f\n") const
	{
		double lo, hi;
		storel(&lo);
		storeh(&hi);
		fprintf(fp, fmt, lo, hi);
	}

	static v2r8 unpckl(const v2r8 a, const v2r8 b){
		return v2r8(__builtin_fj_unpacklo_v2r8(a.val, b.val));
	}
	static v2r8 unpckh(const v2r8 a, const v2r8 b){
		return v2r8(__builtin_fj_unpackhi_v2r8(a.val, b.val));
	}

	v2r8_mask operator<(const v2r8 rhs) const {
		return v2r8_mask(__builtin_fj_cmplt_v2r8(val, rhs.val));
	}
	v2r8 operator&(const v2r8_mask rhs){
		return v2r8(__builtin_fj_and_v2r8(val, rhs.val));
	}
};

// operators for broadcasting the lower element
struct v2r8_bcl{
	typedef __builtin_v2r8 _v2r8;
	_v2r8 val;

	v2r8_bcl(const v2r8_bcl &rhs) : val(rhs.val) {}
	v2r8_bcl operator=(const v2r8_bcl rhs){
		val = rhs.val;
		return (*this);
	}

	v2r8_bcl(const v2r8 src) : val(src.val) {}
	v2r8_bcl(const _v2r8 _val) : val(_val) {}

	v2r8 operator+(const v2r8 rhs) const {
		const v2r8 one(1.0);
		return v2r8(__builtin_fj_madd_cp_v2r8(val, one.val, rhs.val));
	}
	v2r8 operator-(const v2r8 rhs) const {
		const v2r8 one(1.0);
		return v2r8(__builtin_fj_msub_cp_v2r8(val, one.val, rhs.val));
	}
	v2r8 operator*(const v2r8 rhs) const {
		const v2r8 zero(0.0);
		return v2r8(__builtin_fj_madd_cp_v2r8(val, rhs.val, zero.val));
	}
	v2r8 madd(const v2r8 b, const v2r8 c) const {
		return v2r8(__builtin_fj_madd_cp_v2r8(val, b.val, c.val));
	}

	friend v2r8 operator+(const v2r8 lhs, const v2r8_bcl rhs){
		const v2r8 one(1.0);
		return v2r8(__builtin_fj_madd_cp_v2r8(rhs.val, one.val, lhs.val));
	}
	friend v2r8 operator-(const v2r8 lhs, const v2r8_bcl rhs){
		const v2r8 one(1.0);
		return v2r8(__builtin_fj_nmsub_cp_v2r8(rhs.val, one.val, lhs.val));
	}
	friend v2r8 operator*(const v2r8 lhs, const v2r8_bcl rhs){
		const v2r8 zero(0.0);
		return v2r8(__builtin_fj_madd_cp_v2r8(rhs.val, lhs.val, zero.val));
	}
};

// operators for broadcasting the higher element
struct v2r8_bch{
	typedef __builtin_v2r8 _v2r8;
	_v2r8 val;

	v2r8_bch(const v2r8_bch &rhs) : val(rhs.val) {}
	v2r8_bch operator=(const v2r8_bch rhs){
		val = rhs.val;
		return (*this);
	}

	v2r8_bch(const v2r8 src) : val(src.val) {}
	v2r8_bch(const _v2r8 _val) : val(_val) {}

	v2r8 operator+(const v2r8 rhs) const {
		const v2r8 one(1.0);
		return v2r8(__builtin_fj_madd_cp_sr1_v2r8(val, one.val, rhs.val));
	}
	v2r8 operator-(const v2r8 rhs) const {
		const v2r8 one(1.0);
		return v2r8(__builtin_fj_msub_cp_sr1_v2r8(val, one.val, rhs.val));
	}
	v2r8 operator*(const v2r8 rhs) const {
		const v2r8 zero(0.0);
		return v2r8(__builtin_fj_madd_cp_sr1_v2r8(val, rhs.val, zero.val));
	}
	v2r8 madd(const v2r8 b, const v2r8 c) const {
		return v2r8(__builtin_fj_madd_cp_sr1_v2r8(val, b.val, c.val));
	}
	v2r8 nmsub(const v2r8 b, const v2r8 c) const {
		return v2r8(__builtin_fj_nmsub_cp_sr1_v2r8(val, b.val, c.val));
	}

	friend v2r8 operator+(const v2r8 lhs, const v2r8_bch rhs){
		const v2r8 one(1.0);
		return v2r8(__builtin_fj_madd_cp_sr1_v2r8(rhs.val, one.val, lhs.val));
	}
	friend v2r8 operator-(const v2r8 lhs, const v2r8_bch rhs){
		const v2r8 one(1.0);
		return v2r8(__builtin_fj_nmsub_cp_sr1_v2r8(rhs.val, one.val, lhs.val));
	}
	friend v2r8 operator*(const v2r8 lhs, const v2r8_bch rhs){
		const v2r8 zero(0.0);
		return v2r8(__builtin_fj_madd_cp_sr1_v2r8(rhs.val, lhs.val, zero.val));
	}
};
