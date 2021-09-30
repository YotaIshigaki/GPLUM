// gcc -O2 -march=core-avx2
#include <cassert>

#define RSQRTPS_NR_EPJ

class PhantomGrapeQuad{
public:
	enum{
	    //NIMAX = 4096,
        NIMAX = 8192,
	    //NJMAX = 32768,
        NJMAX = 131072,
	};
    
private:
#if 1
	float xibuf  [NIMAX/8]  [3][8];   // x, y, z
	float accpbuf[NIMAX/8]  [4][8];   // ax, ay, az, pot
	float epjbuf [NJMAX]    [4];      // x, y, z, m
	float spjbuf [NJMAX]    [3][4];   // x, y, z, m, | xx, yy, zz, pad, | xy, yz, zx, tr
#else
        float *** xibuf;
        float *** accpbuf;
        float ** epjbuf;
        float *** spjbuf;
#endif
	float eps2;
	static double get_a_NaN(){
		union{ long   l; double d; } m;
		m.l = -1;
		return m.d;
	}

public:
  PhantomGrapeQuad() : eps2(get_a_NaN()) {} // default NaN

  /*
        PhantomGrapeQuad() : eps2(get_a_NaN()) {
	  xibuf = new float**[NIMAX/8];
	  accpbuf = new float**[NIMAX/8];
	  for(int i=0; i<NIMAX/8; i++){
	    xibuf[i] = new float*[3];
	      for(int j=0; j<3; j++){
		xibuf[i][j] = new float[8];
	      }
	      accpbuf[i] = new float*[4];
	      for(int j=0; j<4; j++){
		accpbuf[i][j] = new float[8];
	      }
	  }
	  epjbuf = new float*[NJMAX];
	  spjbuf = new float**[NJMAX];
	  for(int i=0; i<NJMAX; i++){
	    epjbuf[i] = new float[4];
	    spjbuf[i] = new float*[3];
	    for(int j=0; j<3; j++){
	      spjbuf[i][j] = new float[4];
	    }
	  }

	} // default NaN

        ~PhantomGrapeQuad(){
	  for(int i=0; i<NIMAX/8; i++){
	    for(int j=0; j<3; j++){
	      delete [] xibuf[i][j];
	    }
	    for(int j=0; j<4; j++){
	      delete [] accpbuf[i][j];
	    }
	    delete [] xibuf[i];
	    delete [] accpbuf[i];
	  }
	  delete [] xibuf;
	  delete [] accpbuf;

	  for(int i=0; i<NJMAX; i++){
	    for(int j=0; j<3; j++){
	      delete [] spjbuf[i][j];
	    }
	    delete [] spjbuf[i];
	    delete [] epjbuf[i];
	  }
	  delete [] spjbuf;
	  delete [] epjbuf;
        }
  */

	void set_eps2(const double _eps2){
		this->eps2 = _eps2;
	}

	void set_epj_one(const int addr, const double x, const double y, const double z, const double m) {
		epjbuf[addr][0] = x;
		epjbuf[addr][1] = y;
		epjbuf[addr][2] = z;
		epjbuf[addr][3] = m;
	}
	// please specialize the following
	template <typename EPJ_t>
	void set_epj(const int nj, const EPJ_t epj[]);

	void set_spj_one(
			const int addr, 
			const double x,   const double y,   const double z,   const double m,
			const double qxx, const double qyy, const double qzz,
			const double qxy, const double qyz, const double qzx)
	{
		const double tr = qxx + qyy + qzz;
		spjbuf[addr][0][0] = x;
		spjbuf[addr][0][1] = y;
		spjbuf[addr][0][2] = z;
		spjbuf[addr][0][3] = m;

		spjbuf[addr][1][0] = 3.0 * qxx - tr;
		spjbuf[addr][1][1] = 3.0 * qyy - tr;
		spjbuf[addr][1][2] = 3.0 * qzz - tr;
		spjbuf[addr][1][3] = m;

		spjbuf[addr][2][0] = 3.0 * qxy;
		spjbuf[addr][2][1] = 3.0 * qyz;
		spjbuf[addr][2][2] = 3.0 * qzx;
		spjbuf[addr][2][3] = -(eps2 * tr);
	}
	// please specialize the following
	template <typename SPJ_t>
	void set_spj(const int nj, const SPJ_t spj[]);

	void set_xi_one(const int addr, const double x, const double y, const double z){
		const int ah = addr / 8;
		const int al = addr % 8;
		xibuf[ah][0][al] = x;
		xibuf[ah][1][al] = y;
		xibuf[ah][2][al] = z;
	}
	// please specialize the following
	template <typename EPI_t>
	void set_xi(const int ni, const EPI_t spj[]);

	template <typename real_t>
	void get_accp_one(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot){
		const int ah = addr / 8;
		const int al = addr % 8;
		ax  = accpbuf[ah][0][al];
		ay  = accpbuf[ah][1][al];
		az  = accpbuf[ah][2][al];
		pot = accpbuf[ah][3][al];
	}
	template <typename real_t>
	void accum_accp_one(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot){
		const int ah = addr / 8;
		const int al = addr % 8;
		ax  += accpbuf[ah][0][al];
		ay  += accpbuf[ah][1][al];
		az  += accpbuf[ah][2][al];
		pot += accpbuf[ah][3][al];
	}

	void run_epj(const int ni, const int nj){
		assert(ni <= NIMAX);
		assert(nj <= NJMAX);

		kernel_epj_nounroll(ni, nj);
		// kernel_epj_unroll2(ni, nj);
	}
	void run_spj(const int ni, const int nj){
		assert(ni <= NIMAX);
		assert(nj <= NJMAX);

		kernel_spj_nounroll(ni, nj);
		// kernel_spj_unroll2(ni, nj);
	}
private:
	typedef float v4sf __attribute__((vector_size(16)));
	typedef float v8sf __attribute__((vector_size(32)));
	
	__attribute__ ((noinline))
	void kernel_epj_nounroll(const int ni, const int nj){
		const v8sf veps2 = {eps2, eps2, eps2, eps2, eps2, eps2, eps2, eps2};
		for(int i=0; i<ni; i+=8){
			const v8sf xi = *(v8sf *)(xibuf[i/8][0]);
			const v8sf yi = *(v8sf *)(xibuf[i/8][1]);
			const v8sf zi = *(v8sf *)(xibuf[i/8][2]);

			v8sf ax, ay, az, pot;
			ax = ay = az = pot = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

			v8sf jbuf = __builtin_ia32_vbroadcastf128_ps256((v4sf *)epjbuf);
			v8sf xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
			v8sf yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
			v8sf zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
			v8sf mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);
			for(int j=0; j<nj; j++){
				jbuf = __builtin_ia32_vbroadcastf128_ps256((v4sf *)(epjbuf + j+1));

				v8sf dx = xj - xi;
				v8sf dy = yj - yi;
				v8sf dz = zj - zi;

				v8sf r2   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
				v8sf ri1  = __builtin_ia32_rsqrtps256(r2);
#ifdef RSQRTPS_NR_EPJ
				v8sf v3p0 = {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f}; 
				ri1 *= (v3p0 - r2*(ri1*ri1));
#endif
				v8sf mri1 = mj * ri1;
				v8sf ri2  = ri1 * ri1;
				v8sf mri3 = mri1 * ri2;

				xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
				yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
				zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
				mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);

				pot -= mri1;
				ax += mri3 * dx;
				ay += mri3 * dy;
				az += mri3 * dz;
			}
#ifdef RSQRTPS_NR_EPJ
			v8sf v0p5 = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f}; 
			v8sf v0p125 = {0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f}; 
			pot *= v0p5;
			ax  *= v0p125;
			ay  *= v0p125;
			az  *= v0p125;
#endif

			*(v8sf *)(accpbuf[i/8][0]) = ax;
			*(v8sf *)(accpbuf[i/8][1]) = ay;
			*(v8sf *)(accpbuf[i/8][2]) = az;
			*(v8sf *)(accpbuf[i/8][3]) = pot;
		}
	}
	__attribute__ ((noinline))
	void kernel_spj_nounroll(const int ni, const int nj){
		const v8sf veps2 = {eps2, eps2, eps2, eps2, eps2, eps2, eps2, eps2};
		for(int i=0; i<ni; i+=8){
			const v8sf xi = *(v8sf *)(xibuf[i/8][0]);
			const v8sf yi = *(v8sf *)(xibuf[i/8][1]);
			const v8sf zi = *(v8sf *)(xibuf[i/8][2]);

			v8sf ax, ay, az, pot;
			ax = ay = az = pot = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

#define PRELOAD_SPJ

#ifdef PRELOAD_SPJ
			v8sf jbuf0 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[0][0]);
			v8sf jbuf1 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[0][1]);
			v8sf jbuf2 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[0][2]);
#else
			v8sf jbuf0, jbuf1, jbuf2;
#endif
			for(int j=0; j<nj; j++){
#ifndef PRELOAD_SPJ
				jbuf0 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+0][0]);
#endif
				v8sf xj  = __builtin_ia32_shufps256(jbuf0, jbuf0, 0x00);
				v8sf yj  = __builtin_ia32_shufps256(jbuf0, jbuf0, 0x55);
				v8sf zj  = __builtin_ia32_shufps256(jbuf0, jbuf0, 0xaa);
#ifdef PRELOAD_SPJ
				jbuf0 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+1][0]);
#endif

#ifndef PRELOAD_SPJ
				jbuf1 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+0][1]);
#endif
				v8sf qxx = __builtin_ia32_shufps256(jbuf1, jbuf1, 0x00);
				v8sf qyy = __builtin_ia32_shufps256(jbuf1, jbuf1, 0x55);
				v8sf qzz = __builtin_ia32_shufps256(jbuf1, jbuf1, 0xaa);
				v8sf mj  = __builtin_ia32_shufps256(jbuf1, jbuf1, 0xff);
#ifdef PRELOAD_SPJ
				jbuf1 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+1][1]);
#endif

#ifndef PRELOAD_SPJ
				jbuf2 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+0][2]);
#endif
				v8sf qxy = __builtin_ia32_shufps256(jbuf2, jbuf2, 0x00);
				v8sf qyz = __builtin_ia32_shufps256(jbuf2, jbuf2, 0x55);
				v8sf qzx = __builtin_ia32_shufps256(jbuf2, jbuf2, 0xaa);
				v8sf mtr = __builtin_ia32_shufps256(jbuf2, jbuf2, 0xff);
#ifdef PRELOAD_SPJ
				jbuf2 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+1][2]);
#endif

				v8sf dx = xj - xi;
				v8sf dy = yj - yi;
				v8sf dz = zj - zi;

				v8sf r2  = ((veps2 + dx*dx) + dy*dy) + dz*dz;
				v8sf ri1 = __builtin_ia32_rsqrtps256(r2);
				v8sf ri2 = ri1 * ri1;
				v8sf ri3 = ri1 * ri2;
				v8sf ri4 = ri2 * ri2;
				v8sf ri5 = ri2 * ri3;

				v8sf qr_x = (qxx*dx + qxy*dy) + qzx*dz;
				v8sf qr_y = (qyy*dy + qxy*dx) + qyz*dz;
				v8sf qr_z = (qzz*dz + qzx*dx) + qyz*dy;

				v8sf rqr = ((mtr + qr_x*dx) + qr_y*dy) + qr_z*dz;
				v8sf rqr_ri4 = rqr * ri4;

				v8sf v0p5 = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};
				v8sf v2p5 = {2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f}; 

				// rqr_ri4 -= rqr_ri4; // debug

				v8sf meff  =  mj + v0p5 * rqr_ri4;
				v8sf meff3 = (mj + v2p5 * rqr_ri4) * ri3;

				pot -= meff * ri1;

				ax = (ax - ri5*qr_x) + meff3*dx;
				ay = (ay - ri5*qr_y) + meff3*dy;
				az = (az - ri5*qr_z) + meff3*dz;
			}

			*(v8sf *)(accpbuf[i/8][0]) = ax;
			*(v8sf *)(accpbuf[i/8][1]) = ay;
			*(v8sf *)(accpbuf[i/8][2]) = az;
			*(v8sf *)(accpbuf[i/8][3]) = pot;
		}
	}
} __attribute__ ((aligned(128)));
