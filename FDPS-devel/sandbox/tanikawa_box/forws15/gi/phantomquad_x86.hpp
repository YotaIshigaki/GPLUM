// gcc -O2 -march=core-avx2
#include <cassert>

#define RSQRTPS_NR_EPJ

class PhantomGrapeQuad{
public:
	enum{
	    NIMAX = 4096,
//	    NJMAX = 32768,
//	    NJMAX = 65536,
	    NJMAX = 262144,
	};
    
private:
	float xibuf  [NIMAX/8]  [3][8];   // x, y, z
	float accpbuf[NIMAX/8]  [4][8];   // ax, ay, az, pot
	float epjbuf [NJMAX]    [4];      // x, y, z, m
#ifdef USE_AT_PHANTOM
	float spjbuf [NJMAX/2]  [3][8];   // (x,y,z,pad)*2 | (xx,yy,zz,m)*2 | (xy,yz,zx,tr)*2
#else
	float spjbuf [NJMAX]    [3][4];   // x, y, z, pad | xx, yy, zz, m, | xy, yz, zx, tr
#endif
	float eps2;
	static double get_a_NaN(){
		union{ long   l; double d; } m;
		m.l = -1;
		return m.d;
	}
    
public:
    /*
    PhantomGrapeQuad() : eps2(get_a_NaN()) {} // default NaN
    */
    PhantomGrapeQuad() : eps2(get_a_NaN()) {
        assert(0 == size_t(xibuf)   % 32);
        assert(0 == size_t(accpbuf) % 32);
        assert(0 == size_t(epjbuf)  % 16);
        assert(0 == size_t(spjbuf)  % 32);
    } // default NaN
    
	void set_eps2(const double _eps2){
		this->eps2 = _eps2;
	}
    
	void set_epj_one(const int addr,
                     const double x,
                     const double y,
                     const double z,
                     const double m) {
        
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
#ifdef USE_AT_PHANTOM
            const int ah = addr / 2;
            const int al = 4 * (addr % 2);
            const double tr = qxx + qyy + qzz;
            spjbuf[ah][0][al+0] = x;
            spjbuf[ah][0][al+1] = y;
            spjbuf[ah][0][al+2] = z;
            
            spjbuf[ah][1][al+0] = 3.0 * qxx - tr;
            spjbuf[ah][1][al+1] = 3.0 * qyy - tr;
            spjbuf[ah][1][al+2] = 3.0 * qzz - tr;
            spjbuf[ah][1][al+3] = m;
            
            spjbuf[ah][2][al+0] = 3.0 * qxy;
            spjbuf[ah][2][al+1] = 3.0 * qyz;
            spjbuf[ah][2][al+2] = 3.0 * qzx;
            spjbuf[ah][2][al+3] = -(eps2 * tr);
#else
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
#endif
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

    void check_range(const int ni, const int nj) {
        assert(ni <= NIMAX);
        assert(nj <= NJMAX);
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

#ifdef USE_AT_PHANTOM
	__attribute__ ((noinline))
	void kernel_epj_nounroll(const int ni, const int nj){
		const v8sf veps2 = {eps2, eps2, eps2, eps2, eps2, eps2, eps2, eps2};
        if(nj % 2 == 1) {
            epjbuf[nj+1][0] = 0.0;
            epjbuf[nj+1][1] = 0.0;
            epjbuf[nj+1][2] = 0.0;
            epjbuf[nj+1][3] = 0.0;
        }

		for(int i=0; i<ni; i+=4){
            const int  ii = i % 8;
			const v8sf xi = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&xibuf[i/8][0][ii]);
			const v8sf yi = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&xibuf[i/8][1][ii]);
			const v8sf zi = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&xibuf[i/8][2][ii]);

			v8sf ax, ay, az, pot;
			ax = ay = az = pot = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

            v8sf jbuf = *(v8sf *)&epjbuf[0];
			v8sf xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
			v8sf yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
			v8sf zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
			v8sf mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);
			for(int j=0; j<nj; j+=2){
                v8sf jbuf = *(v8sf *)&epjbuf[j+2];
                
				v8sf dx = xj - xi;
				v8sf dy = yj - yi;
				v8sf dz = zj - zi;

				xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
				yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
				zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
				mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);
                
				v8sf r2   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
				v8sf ri1  = __builtin_ia32_rsqrtps256(r2);
#ifdef RSQRTPS_NR_EPJ
				v8sf v3p0 = {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f}; 
				ri1 *= (v3p0 - r2*(ri1*ri1));
#endif
				v8sf mri1 = mj * ri1;
				v8sf ri2  = ri1 * ri1;
				v8sf mri3 = mri1 * ri2;
                
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

            *(v4sf *)(&accpbuf[i/8][0][ii]) = *(v4sf *)&ax[0] + *(v4sf *)&ax[4];
            *(v4sf *)(&accpbuf[i/8][1][ii]) = *(v4sf *)&ay[0] + *(v4sf *)&ay[4];
            *(v4sf *)(&accpbuf[i/8][2][ii]) = *(v4sf *)&az[0] + *(v4sf *)&az[4];
            *(v4sf *)(&accpbuf[i/8][3][ii]) = *(v4sf *)&pot[0] + *(v4sf *)&pot[4];

		}
	}
	__attribute__ ((noinline))
	void kernel_spj_nounroll(const int ni, const int nj){
		const v8sf veps2 = {eps2, eps2, eps2, eps2, eps2, eps2, eps2, eps2};
        if(nj % 2 == 1) {
            const int ah = nj / 2;
            spjbuf[ah][0][4] = spjbuf[ah][0][5] = spjbuf[ah][0][6] = spjbuf[ah][0][7] = 0.0;
            spjbuf[ah][1][4] = spjbuf[ah][1][5] = spjbuf[ah][1][6] = spjbuf[ah][1][7] = 0.0;
            spjbuf[ah][2][4] = spjbuf[ah][2][5] = spjbuf[ah][2][6] = spjbuf[ah][2][7] = 0.0;
        }
		for(int i=0; i<ni; i+=4){
            const int  ii = i % 8;
			const v8sf xi = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&xibuf[i/8][0][ii]);
			const v8sf yi = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&xibuf[i/8][1][ii]);
			const v8sf zi = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&xibuf[i/8][2][ii]);
            
			v8sf ax, ay, az, pot;
			ax = ay = az = pot = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
            
			v8sf jbuf0 = *(v8sf *)&spjbuf[0][0];
			v8sf jbuf1 = *(v8sf *)&spjbuf[0][1];
			v8sf jbuf2 = *(v8sf *)&spjbuf[0][2];
			for(int j=0, jj=0; j<nj; j+=2, jj++){
				v8sf xj  = __builtin_ia32_shufps256(jbuf0, jbuf0, 0x00);
				v8sf yj  = __builtin_ia32_shufps256(jbuf0, jbuf0, 0x55);
				v8sf zj  = __builtin_ia32_shufps256(jbuf0, jbuf0, 0xaa);
				jbuf0 = *(v8sf *)&spjbuf[jj+1][0];
                
				v8sf qxx = __builtin_ia32_shufps256(jbuf1, jbuf1, 0x00);
				v8sf qyy = __builtin_ia32_shufps256(jbuf1, jbuf1, 0x55);
				v8sf qzz = __builtin_ia32_shufps256(jbuf1, jbuf1, 0xaa);
				v8sf mj  = __builtin_ia32_shufps256(jbuf1, jbuf1, 0xff);
				jbuf1 = *(v8sf *)&spjbuf[jj+1][1];
                
				v8sf qxy = __builtin_ia32_shufps256(jbuf2, jbuf2, 0x00);
				v8sf qyz = __builtin_ia32_shufps256(jbuf2, jbuf2, 0x55);
				v8sf qzx = __builtin_ia32_shufps256(jbuf2, jbuf2, 0xaa);
				v8sf mtr = __builtin_ia32_shufps256(jbuf2, jbuf2, 0xff);
				jbuf2 = *(v8sf *)&spjbuf[jj+1][2];
                
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
            
			*(v4sf *)(&accpbuf[i/8][0][ii]) = *(v4sf *)&ax[0] + *(v4sf *)&ax[4];
			*(v4sf *)(&accpbuf[i/8][1][ii]) = *(v4sf *)&ay[0] + *(v4sf *)&ay[4];
			*(v4sf *)(&accpbuf[i/8][2][ii]) = *(v4sf *)&az[0] + *(v4sf *)&az[4];
			*(v4sf *)(&accpbuf[i/8][3][ii]) = *(v4sf *)&pot[0] + *(v4sf *)&pot[4];
		}
	}
#else
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
#endif
} __attribute__ ((aligned(128)));
