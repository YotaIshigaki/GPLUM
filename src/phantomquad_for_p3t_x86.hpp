// gcc -O2 -march=core-avx2

#pragma once

#include <cassert>
#include <iostream>

#ifdef __AVX2__
#include <immintrin.h>
#elif defined(__AVX512DQ__)
#include <immintrin.h>
//#include <zmmintrin.h>
#endif


class PhantomGrapeQuad{
    
#if 1
    
public:
    enum{
        //NIMAX = 32768,
        //NJMAX = 131072,
        NIMAX = 512,
        NJMAX = 2048,
    };
    
private:
#ifdef USE_INDIVIDUAL_CUTOFF
    float  xibuf   [NIMAX/8]  [5][8];   // x, y, z, r_out, r_in;
    double xibufd  [NIMAX/8]  [5][8];   // x, y, z, r_out, r_in;
    float  epjbuf  [NJMAX]    [2][4];   // x, y, z, m, | r_out, r_in
    double epjbufd [NJMAX]    [2][4];   // x, y, z, m, | r_out, r_in
#else
    float  xibuf   [NIMAX/8]  [3][8];   // x, y, z
    double xibufd  [NIMAX/8]  [3][8];   // x, y, z
    float  epjbuf  [NJMAX]    [4];      // x, y, z, m
    double epjbufd [NJMAX]    [4];      // x, y, z, m
#endif
    float  spjbuf  [NJMAX]    [3][4];   // x, y, z, m, | xx, yy, zz, pad, | xy, yz, zx, tr
    double spjbufd [NJMAX]    [3][4];   // x, y, z, m, | xx, yy, zz, pad, | xy, yz, zx, tr
    float  accpbuf [NIMAX/8]  [4][8];   // ax, ay, az, pot
    double accpbufd[NIMAX/8]  [4][8];   // ax, ay, az, pot
    
#else //if 0

private:
#ifdef USE_INDIVIDUAL_CUTOFF
    float  (*xibuf)    [5][8];  // x, y, z, r_out, r_in;
    double (*xibufd)   [5][8];  // x, y, z, r_out, r_in;
    float  (*epjbuf)   [2][4];  // x, y, z, m, | r_out, r_in
    double (*epjbufd)  [2][4];  // x, y, z, m, | r_out, r_in
#else
    float  (*xibuf)    [3][8];  // x, y, z
    double (*xibufd)   [3][8];  // x, y, z
    float  (*epjbuf)   [4];     // x, y, z, m
    double (*epjbufd)  [4];     // x, y, z, m
#endif
    float  (*spjbuf)   [3][4];  // x, y, z, m, | xx, yy, zz, pad, | xy, yz, zx, tr
    double (*spjbufd)  [3][4];  // x, y, z, m, | xx, yy, zz, pad, | xy, yz, zx, tr
    float  (*accpbuf)  [4][8];  // ax, ay, az, pot
    double (*accpbufd) [4][8];  // ax, ay, az, pot

#endif

    double eps2;
    static double get_a_NaN(){
        union{ long   l; double d; } m;
        m.l = -1;
        return m.d;
    }
    double r_out;
    double r_in;
    double denominator; // for cut off
public:

    PhantomGrapeQuad() : eps2(get_a_NaN()) {} // default NaN

    void set_cutoff(const double _r_out, const double _r_in){
        r_out = _r_out;
        r_in = _r_in;
        denominator = 1.0 / (r_out - r_in);
    }

    void set_eps2(const double _eps2){
        this->eps2 = _eps2;
    }

#ifdef USE_INDIVIDUAL_CUTOFF
    void set_epj_one(const int addr, const double x, const double y, const double z, const double m,
                     const double r_out, const double r_in) {
        epjbuf[addr][0][0] = x;
        epjbuf[addr][0][1] = y;
        epjbuf[addr][0][2] = z;
        epjbuf[addr][0][3] = m;

        epjbuf[addr][1][0] = r_out;
        epjbuf[addr][1][1] = r_in;
        epjbuf[addr][1][2] = 0.;
        epjbuf[addr][1][3] = 0.;
    }
    void set_epj_one_d(const int addr, const double x, const double y, const double z, const double m,
                       const double r_out, const double r_in){
        epjbufd[addr][0][0] = x;
        epjbufd[addr][0][1] = y;
        epjbufd[addr][0][2] = z;
        epjbufd[addr][0][3] = m;

        epjbufd[addr][1][0] = r_out;
        epjbufd[addr][1][1] = r_in;
        epjbufd[addr][1][2] = 0.;
        epjbufd[addr][1][3] = 0.;
    }
#else //USE_INDIVIDUAL_CUTOFF
    void set_epj_one(const int addr, const double x, const double y, const double z, const double m) {
        epjbuf[addr][0] = x;
        epjbuf[addr][1] = y;
        epjbuf[addr][2] = z;
        epjbuf[addr][3] = m;
    }
    void set_epj_one_d(const int addr, const double x, const double y, const double z, const double m){
        epjbufd[addr][0] = x;
        epjbufd[addr][1] = y;
        epjbufd[addr][2] = z;
        epjbufd[addr][3] = m;
    }
#endif //USE_INDIVIDUAL_CUTOFF

    void set_spj_one(const int addr, 
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
    void set_spj_one_d(const int addr, 
                       const double x,   const double y,   const double z,   const double m,
                       const double qxx, const double qyy, const double qzz,
                       const double qxy, const double qyz, const double qzx)
    {
        const double tr = qxx + qyy + qzz;
        spjbufd[addr][0][0] = x;
        spjbufd[addr][0][1] = y;
        spjbufd[addr][0][2] = z;
        spjbufd[addr][0][3] = m;

        spjbufd[addr][1][0] = 3.0 * qxx - tr;
        spjbufd[addr][1][1] = 3.0 * qyy - tr;
        spjbufd[addr][1][2] = 3.0 * qzz - tr;
        spjbufd[addr][1][3] = m;

        spjbufd[addr][2][0] = 3.0 * qxy;
        spjbufd[addr][2][1] = 3.0 * qyz;
        spjbufd[addr][2][2] = 3.0 * qzx;
        spjbufd[addr][2][3] = -(eps2 * tr);
    }

#ifdef USE_INDIVIDUAL_CUTOFF
    void set_xi_one(const int addr, const double x, const double y, const double z,
                    const double r_out, const double r_in){
        const int ah = addr / 8;
        const int al = addr % 8;
        xibuf[ah][0][al] = x;
        xibuf[ah][1][al] = y;
        xibuf[ah][2][al] = z;
        xibuf[ah][3][al] = r_out;
        xibuf[ah][4][al] = r_in;
    }
    void set_xi_one_d(const int addr, const double x, const double y, const double z,
                      const double r_out, const double r_in){
        const int ah = addr / 8;
        const int al = addr % 8;
        xibufd[ah][0][al] = x;
        xibufd[ah][1][al] = y;
        xibufd[ah][2][al] = z;
        xibufd[ah][3][al] = r_out;
        xibufd[ah][4][al] = r_in;
    }
#else //USE_INDIVIDUAL_CUTOFF
    void set_xi_one(const int addr, const double x, const double y, const double z){
        const int ah = addr / 8;
        const int al = addr % 8;
        xibuf[ah][0][al] = x;
        xibuf[ah][1][al] = y;
        xibuf[ah][2][al] = z;
    }
    void set_xi_one_d(const int addr, const double x, const double y, const double z){
        const int ah = addr / 8;
        const int al = addr % 8;
        xibufd[ah][0][al] = x;
        xibufd[ah][1][al] = y;
        xibufd[ah][2][al] = z;
    }
#endif //USE_INDIVIDUAL_CUTOFF

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
    void get_accp_one_d(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot){
        const int ah = addr / 8;
        const int al = addr % 8;
        ax  = accpbufd[ah][0][al];
        ay  = accpbufd[ah][1][al];
        az  = accpbufd[ah][2][al];
        pot = accpbufd[ah][3][al];
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
    template <typename real_t>
    void accum_accp_one_d(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot){
        const int ah = addr / 8;
        const int al = addr % 8;
        ax  += accpbufd[ah][0][al];
        ay  += accpbufd[ah][1][al];
        az  += accpbufd[ah][2][al];
        pot += accpbufd[ah][3][al];
    }

    void run_epj(const int ni, const int nj){
        if(ni > NIMAX || nj > NJMAX){
            std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
            for(int i=0; i<(ni-1)/8+1; i++){
                for(int j=0; j<8; j++){
                    for(int k=0; k<3; k++){
                        std::cout<<"i,j,k="<<i<<" "<<j<<" "<<k<<std::endl;
                        std::cout<<"xibuf[i][k][j]="<<xibuf[i][k][j]<<std::endl;
                    }
                    std::cout<<std::endl;
                }
            }
        }
        assert(ni <= NIMAX);
        assert(nj <= NJMAX);

        kernel_epj_nounroll(ni, nj);
    }

    void run_epj_d(const int ni, const int nj){
        if(ni > NIMAX || nj > NJMAX){
            std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
            for(int i=0; i<(ni-1)/8+1; i++){
                for(int j=0; j<8; j++){
                    for(int k=0; k<3; k++){
                        std::cout<<"i,j,k="<<i<<" "<<j<<" "<<k<<std::endl;
                        std::cout<<"xibuf[i][k][j]="<<xibuf[i][k][j]<<std::endl;
                    }
                    std::cout<<std::endl;
                }
            }
        }
        assert(ni <= NIMAX);
        assert(nj <= NJMAX);
        kernel_epj_nounroll_64bit(ni, nj);
    }

    //////////
    // include linear cutoff
    void run_epj_for_p3t_with_linear_cutoff(const int ni, const int nj){
        if(ni > NIMAX || nj > NJMAX){
            std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
            for(int i=0; i<(ni-1)/8+1; i++){
                for(int j=0; j<8; j++){
                    for(int k=0; k<3; k++){
                        std::cout<<"xibufd[i][k][j]="<<xibufd[i][k][j]<<std::endl;
                    }
                    std::cout<<std::endl;
                }
            }
        }
        assert(ni <= NIMAX);
        assert(nj <= NJMAX);
        kernel_epj_nounroll_for_p3t_with_linear_cutoff(ni, nj);
    }

    void run_epj_for_p3t_with_linear_cutoff_d(const int ni, const int nj){
        if(ni > NIMAX || nj > NJMAX){
            std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
            for(int i=0; i<(ni-1)/8+1; i++){
                for(int j=0; j<8; j++){
                    for(int k=0; k<3; k++){
                        std::cout<<"xibufd[i][k][j]="<<xibufd[i][k][j]<<std::endl;
                    }
                    std::cout<<std::endl;
                }
            }
        }
        assert(ni <= NIMAX);
        assert(nj <= NJMAX);
        kernel_epj_64bit_nounroll_for_p3t_with_linear_cutoff(ni, nj);
    }
    // include linear cutoff
    //////////

    void run_spj(const int ni, const int nj){
        if(ni > NIMAX || nj > NJMAX){
            std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
        }
        assert(ni <= NIMAX);
        assert(nj <= NJMAX);

        kernel_spj_nounroll(ni, nj);
        // kernel_spj_unroll2(ni, nj);
    }

    void run_spj_d(const int ni, const int nj){
        if(ni > NIMAX || nj > NJMAX){
            std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
        }
        assert(ni <= NIMAX);
        assert(nj <= NJMAX);

        kernel_spj_64bit_nounroll(ni, nj);
        // kernel_spj_unroll2(ni, nj);
    }

private:
	typedef float v4sf __attribute__((vector_size(16)));
	typedef float v8sf __attribute__((vector_size(32)));
	typedef double v4df __attribute__((vector_size(32)));

#ifdef __AVX2__
    
	__attribute__ ((noinline))
	void kernel_epj_nounroll(const int ni, const int nj){

	    const v8sf veps2 = {(float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2};
        
		for(int i=0; i<ni; i+=8){
            const v8sf xi = *(v8sf *)(xibuf[i/8][0]);
			const v8sf yi = *(v8sf *)(xibuf[i/8][1]);
			const v8sf zi = *(v8sf *)(xibuf[i/8][2]);

			v8sf ax, ay, az, pot;
			ax = ay = az = pot = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
			//v8sf jbuf = __builtin_ia32_vbroadcastf128_ps256((v4sf *)epjbuf);
			//v8sf xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
			//v8sf yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
			//v8sf zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
			//v8sf mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);
#ifdef USE_INDIVIDUAL_CUTOFF
            v8sf jbuf = _mm256_broadcast_ps((v4sf *)&epjbuf[0][0]);
#else
            v8sf jbuf = _mm256_broadcast_ps((v4sf *)epjbuf);
#endif
            v8sf xj =  _mm256_shuffle_ps(jbuf, jbuf, 0x00);
            v8sf yj =  _mm256_shuffle_ps(jbuf, jbuf, 0x55);
            v8sf zj =  _mm256_shuffle_ps(jbuf, jbuf, 0xaa);
            v8sf mj =  _mm256_shuffle_ps(jbuf, jbuf, 0xff);
			for(int j=0; j<nj; j++){
				//jbuf = __builtin_ia32_vbroadcastf128_ps256((v4sf *)(epjbuf + j+1));
#ifdef USE_INDIVIDUAL_CUTOFF
                jbuf = _mm256_broadcast_ps((v4sf *)&epjbuf[j+1][0]);
#else
                jbuf = _mm256_broadcast_ps((v4sf*)(epjbuf+j+1));
#endif

				v8sf dx = xj - xi;
				v8sf dy = yj - yi;
				v8sf dz = zj - zi;

				v8sf r2   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                //v8sf mask = _mm256_cmp_ps(veps2, r2, 0x4);
				//v8sf ri1  = __builtin_ia32_rsqrtps256(r2);
                v8sf ri1  = _mm256_rsqrt_ps(r2);
                //v8sf ri1  = _mm256_and_pd( __mm256_rsqrt_ps( r2), mask );
#ifdef RSQRT_NR_EPJ_X2
				v8sf v3p0 = {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f}; 
				ri1 *= (v3p0 - r2*(ri1*ri1));
#endif
				v8sf mri1 = mj * ri1;
				v8sf ri2  = ri1 * ri1;
				v8sf mri3 = mri1 * ri2;

				//xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
				//yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
				//zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
				//mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);
                xj =  _mm256_shuffle_ps(jbuf, jbuf, 0x00);
                yj =  _mm256_shuffle_ps(jbuf, jbuf, 0x55);
                zj =  _mm256_shuffle_ps(jbuf, jbuf, 0xaa);
                mj =  _mm256_shuffle_ps(jbuf, jbuf, 0xff);

				pot -= mri1;
				ax += mri3 * dx;
				ay += mri3 * dy;
				az += mri3 * dz;
			}
#ifdef RSQRT_NR_EPJ_X2
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
    void kernel_epj_nounroll_64bit(const int ni, const int nj){
        //const v4df vzero = {0.0, 0.0, 0.0, 0.0};

        const v4df veps2 = {eps2, eps2, eps2, eps2};
        for(int i=0; i<ni; i+=4){
            const int il = i%8;
            const v4df xi = *(v4df *)(&xibufd[i/8][0][il]);
            const v4df yi = *(v4df *)(&xibufd[i/8][1][il]);
            const v4df zi = *(v4df *)(&xibufd[i/8][2][il]);

            v4df ax, ay, az, pot;
            ax = ay = az = pot = (v4df){0.0, 0.0, 0.0, 0.0};

#ifdef USE_INDIVIDUAL_CUTOFF
            v4df jbuf = *((v4df *)&epjbufd[0][0]);
#else
            v4df jbuf = *((v4df*)epjbufd);
#endif
            v4df xj =  _mm256_permute4x64_pd(jbuf, 0x00);
            v4df yj =  _mm256_permute4x64_pd(jbuf, 0x55);
            v4df zj =  _mm256_permute4x64_pd(jbuf, 0xaa);
            v4df mj =  _mm256_permute4x64_pd(jbuf, 0xff);

            for(int j=0; j<nj; j++){
#ifdef USE_INDIVIDUAL_CUTOFF
                jbuf = *((v4df *)&epjbufd[j+1][0]);
#else
                jbuf = *((v4df*)(epjbufd+j+1));
#endif
                v4df dx = xj - xi;
                v4df dy = yj - yi;
                v4df dz = zj - zi;
                v4df r2   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                //v4df mask = _mm256_cmp_pd(vrcrit2, r2, 0x01); // vrcrit2 < r2
                //v4df mask = _mm256_cmp_pd(veps2, r2, 0x4); // veps2 != r2
                //v4df ri1  = _mm256_and_pd( __builtin_ia32_cvtps2pd256( __builtin_ia32_rsqrtps( __builtin_ia32_cvtpd2ps256(r2))), mask );
                //v4df ri1  = _mm256_and_pd( _mm256_cvtps_pd( _mm_rsqrt_ps( _mm256_cvtpd_ps(r2))), mask );
                v4df ri1  = _mm256_cvtps_pd( _mm_rsqrt_ps( _mm256_cvtpd_ps(r2)));
#ifdef RSQRT_NR_EPJ_X2
                //x2
                v4df v0p5 = {0.5, 0.5, 0.5, 0.5};
                v4df v3p0 = {3.0, 3.0, 3.0, 3.0};
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
                // x4
                v4df vone = {1.0, 1.0, 1.0, 1.0};
                v4df v8p0 = {8.0, 8.0, 8.0, 8.0};
                v4df v6p0 = {6.0, 6.0, 6.0, 6.0};
                v4df v5p0 = {5.0, 5.0, 5.0, 5.0};
                v4df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
                v4df h = vone - r2*(ri1*ri1);
                ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
                v4df mri1 = mj * ri1;
                v4df ri2  = ri1 * ri1;
                v4df mri3 = mri1 * ri2;
		
                xj =  _mm256_permute4x64_pd(jbuf, 0x00);
                yj =  _mm256_permute4x64_pd(jbuf, 0x55);
                zj =  _mm256_permute4x64_pd(jbuf, 0xaa);
                mj =  _mm256_permute4x64_pd(jbuf, 0xff);

                pot -= mri1;
                ax += mri3 * dx;
                ay += mri3 * dy;
                az += mri3 * dz;
            }
            *(v4df *)(&accpbufd[i/8][0][il]) = ax;
            *(v4df *)(&accpbufd[i/8][1][il]) = ay;
            *(v4df *)(&accpbufd[i/8][2][il]) = az;
            *(v4df *)(&accpbufd[i/8][3][il]) = pot;
        }
    }

    ///////////////////////
    // with linear cutoff
    __attribute__ ((noinline))
    void kernel_epj_nounroll_for_p3t_with_linear_cutoff(const int ni, const int nj){
        
        const v8sf vone = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
        const v8sf veps2 = {(float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2};
#ifndef USE_INDIVIDUAL_CUTOFF
        const v8sf vr_out  = {(float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out};
        const v8sf vr_out2  = vr_out * vr_out;
#endif
        for(int i=0; i<ni; i+=8){
            const v8sf xi = *(v8sf *)(xibuf[i/8][0]);
            const v8sf yi = *(v8sf *)(xibuf[i/8][1]);
            const v8sf zi = *(v8sf *)(xibuf[i/8][2]);
#ifdef USE_INDIVIDUAL_CUTOFF
            const v8sf r_outi = *(v8sf *)(xibuf[i/8][3]);
#endif
            
            v8sf ax, ay, az, pot;
            ax = ay = az = pot = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
            //v8sf jbuf = __builtin_ia32_vbroadcastf128_ps256((v4sf *)epjbuf);
            //v8sf xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
            //v8sf yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
            //v8sf zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
            //v8sf mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);
#ifdef USE_INDIVIDUAL_CUTOFF
            v8sf jbuf0 = _mm256_broadcast_ps((v4sf *)&epjbuf[0][0]);
            v8sf jbuf1 = _mm256_broadcast_ps((v4sf *)&epjbuf[0][1]);
            v8sf xj     =  _mm256_shuffle_ps(jbuf0, jbuf0, 0x00);
            v8sf yj     =  _mm256_shuffle_ps(jbuf0, jbuf0, 0x55);
            v8sf zj     =  _mm256_shuffle_ps(jbuf0, jbuf0, 0xaa);
            v8sf mj     =  _mm256_shuffle_ps(jbuf0, jbuf0, 0xff);
            v8sf r_outj =  _mm256_shuffle_ps(jbuf1, jbuf1, 0x00);
#else
            v8sf jbuf = _mm256_broadcast_ps((v4sf *)epjbuf);
            v8sf xj =  _mm256_shuffle_ps(jbuf, jbuf, 0x00);
            v8sf yj =  _mm256_shuffle_ps(jbuf, jbuf, 0x55);
            v8sf zj =  _mm256_shuffle_ps(jbuf, jbuf, 0xaa);
            v8sf mj =  _mm256_shuffle_ps(jbuf, jbuf, 0xff);
#endif
            
            for(int j=0; j<nj; j++){
                //jbuf = __builtin_ia32_vbroadcastf128_ps256((v4sf *)(epjbuf + j+1));
#ifdef USE_INDIVIDUAL_CUTOFF
                jbuf0 = _mm256_broadcast_ps((v4sf *)&epjbuf[j+1][0]);
                jbuf1 = _mm256_broadcast_ps((v4sf *)&epjbuf[j+1][1]);
                v8sf vr_out  = _mm256_max_ps(r_outi, r_outj);
                v8sf vr_out2 = vr_out * vr_out;
#else
                jbuf = _mm256_broadcast_ps((v4sf*)(epjbuf+j+1));
#endif
                v8sf dx = xj - xi;
                v8sf dy = yj - yi;
                v8sf dz = zj - zi;
                v8sf r2_real   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v8sf r2 = _mm256_max_ps(r2_real, vr_out2);
                //v8sf ri1  = __builtin_ia32_rsqrtps256(r2);
                v8sf ri1  = _mm256_rsqrt_ps(r2);
#ifdef RSQRT_NR_EPJ_X2
                v8sf v3p0 = {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f}; 
                v8sf v0p5 = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
                // x4
                v8sf v8p0 = {8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f};
                v8sf v6p0 = {6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f};
                v8sf v5p0 = {5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f};
                v8sf v0p0625 = {(float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0};
                v8sf h = vone - r2*(ri1*ri1);
                ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
                v8sf ri2 = ri1*ri1;
                v8sf mri1 = mj*ri1;
                v8sf mri3 = mri1 * ri2;
                //xj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x00);
                //yj =  __builtin_ia32_shufps256(jbuf, jbuf, 0x55);
                //zj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xaa);
                //mj =  __builtin_ia32_shufps256(jbuf, jbuf, 0xff);
#ifdef USE_INDIVIDUAL_CUTOFF
                xj     =  _mm256_shuffle_ps(jbuf0, jbuf0, 0x00);
                yj     =  _mm256_shuffle_ps(jbuf0, jbuf0, 0x55);
                zj     =  _mm256_shuffle_ps(jbuf0, jbuf0, 0xaa);
                mj     =  _mm256_shuffle_ps(jbuf0, jbuf0, 0xff);
                r_outj =  _mm256_shuffle_ps(jbuf1, jbuf1, 0x00);
#else
                xj =  _mm256_shuffle_ps(jbuf, jbuf, 0x00);
                yj =  _mm256_shuffle_ps(jbuf, jbuf, 0x55);
                zj =  _mm256_shuffle_ps(jbuf, jbuf, 0xaa);
                mj =  _mm256_shuffle_ps(jbuf, jbuf, 0xff);
#endif
                pot -= mri1;
                ax += mri3 * dx;
                ay += mri3 * dy;
                az += mri3 * dz;
            }
            *(v8sf *)(accpbuf[i/8][0]) = ax;
            *(v8sf *)(accpbuf[i/8][1]) = ay;
            *(v8sf *)(accpbuf[i/8][2]) = az;
            *(v8sf *)(accpbuf[i/8][3]) = pot;
        }
    }

    __attribute__ ((noinline))
    void kernel_epj_64bit_nounroll_for_p3t_with_linear_cutoff(const int ni, const int nj){
        const v4df vone = {1.0, 1.0, 1.0, 1.0};
        const v4df veps2 = {eps2, eps2, eps2, eps2};
#ifndef USE_INDIVIDUAL_CUTOFF
        const v4df vr_out  = {r_out,  r_out,  r_out,  r_out};
        const v4df vr_out2  = vr_out * vr_out;
#endif
        for(int i=0; i<ni; i+=4){
            const int il = i%8;
            const v4df xi = *(v4df *)(&xibufd[i/8][0][il]);
            const v4df yi = *(v4df *)(&xibufd[i/8][1][il]);
            const v4df zi = *(v4df *)(&xibufd[i/8][2][il]);
#ifdef USE_INDIVIDUAL_CUTOFF
            const v4df r_outi = *(v4df *)(&xibufd[i/8][3][il]);
#endif
            
            v4df ax, ay, az, pot;
            ax = ay = az = pot = (v4df){0.0, 0.0, 0.0, 0.0};
#ifdef USE_INDIVIDUAL_CUTOFF
            v4df jbuf0 = *((v4df *)&epjbufd[0][0]);
            v4df jbuf1 = *((v4df *)&epjbufd[0][1]);
            v4df xj     =  _mm256_permute4x64_pd(jbuf0, 0x00);
            v4df yj     =  _mm256_permute4x64_pd(jbuf0, 0x55);
            v4df zj     =  _mm256_permute4x64_pd(jbuf0, 0xaa);
            v4df mj     =  _mm256_permute4x64_pd(jbuf0, 0xff);
            v4df r_outj =  _mm256_permute4x64_pd(jbuf1, 0x00);
#else
            v4df jbuf = *((v4df*)epjbufd);
            v4df xj =  _mm256_permute4x64_pd(jbuf, 0x00);
            v4df yj =  _mm256_permute4x64_pd(jbuf, 0x55);
            v4df zj =  _mm256_permute4x64_pd(jbuf, 0xaa);
            v4df mj =  _mm256_permute4x64_pd(jbuf, 0xff);
#endif
            
            for(int j=0; j<nj; j++){
#ifdef USE_INDIVIDUAL_CUTOFF
                jbuf0 = *((v4df *)&epjbufd[j+1][0]);
                jbuf1 = *((v4df *)&epjbufd[j+1][1]);
                v4df vr_out  = _mm256_max_pd(r_outi, r_outj);
                v4df vr_out2 = vr_out * vr_out;
#else
                jbuf = *((v4df*)(epjbufd+j+1));
#endif
                v4df dx = xj - xi;
                v4df dy = yj - yi;
                v4df dz = zj - zi;
                v4df r2_real   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v4df r2 = _mm256_max_pd( r2_real, vr_out2);
                //v4df ri1  = __builtin_ia32_cvtps2pd256(__builtin_ia32_rsqrtps( __builtin_ia32_cvtpd2ps256(r2)));
                v4df ri1  = _mm256_cvtps_pd( _mm_rsqrt_ps( _mm256_cvtpd_ps(r2)));
#ifdef RSQRT_NR_EPJ_X2
                //x2
                v4df v0p5 = {0.5, 0.5, 0.5, 0.5};
                v4df v3p0 = {3.0, 3.0, 3.0, 3.0};
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
                // x4
                v4df v8p0 = {8.0, 8.0, 8.0, 8.0};
                v4df v6p0 = {6.0, 6.0, 6.0, 6.0};
                v4df v5p0 = {5.0, 5.0, 5.0, 5.0};
                v4df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
                v4df h = vone - r2*(ri1*ri1);
                ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
                v4df ri2 = ri1*ri1;
                v4df mri1 = mj*ri1;
                v4df mri3 = mri1 * ri2;
#ifdef USE_INDIVIDUAL_CUTOFF
                xj     =  _mm256_permute4x64_pd(jbuf0, 0x00);
                yj     =  _mm256_permute4x64_pd(jbuf0, 0x55);
                zj     =  _mm256_permute4x64_pd(jbuf0, 0xaa);
                mj     =  _mm256_permute4x64_pd(jbuf0, 0xff);
                r_outj =  _mm256_permute4x64_pd(jbuf1, 0x00);
#else
                xj =  _mm256_permute4x64_pd(jbuf, 0x00);
                yj =  _mm256_permute4x64_pd(jbuf, 0x55);
                zj =  _mm256_permute4x64_pd(jbuf, 0xaa);
                mj =  _mm256_permute4x64_pd(jbuf, 0xff);
#endif
                pot -= mri1;
                ax += mri3 * dx;
                ay += mri3 * dy;
                az += mri3 * dz;
            }
            *(v4df *)(&accpbufd[i/8][0][il]) = ax;
            *(v4df *)(&accpbufd[i/8][1][il]) = ay;
            *(v4df *)(&accpbufd[i/8][2][il]) = az;
            *(v4df *)(&accpbufd[i/8][3][il]) = pot;
        }
    }
    // linear cutoff
    //////////////
    

	__attribute__ ((noinline))
	void kernel_spj_nounroll(const int ni, const int nj){
        
	    const v8sf veps2 = {(float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2};
		for(int i=0; i<ni; i+=8){
			const v8sf xi = *(v8sf *)(xibuf[i/8][0]);
			const v8sf yi = *(v8sf *)(xibuf[i/8][1]);
			const v8sf zi = *(v8sf *)(xibuf[i/8][2]);

			v8sf ax, ay, az, pot;
			ax = ay = az = pot = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

#define PRELOAD_SPJ

#ifdef PRELOAD_SPJ
			//v8sf jbuf0 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[0][0]);
			//v8sf jbuf1 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[0][1]);
			//v8sf jbuf2 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[0][2]);
            v8sf jbuf0 = _mm256_broadcast_ps((v4sf *)&spjbuf[0][0]);
            v8sf jbuf1 = _mm256_broadcast_ps((v4sf *)&spjbuf[0][1]);
            v8sf jbuf2 = _mm256_broadcast_ps((v4sf *)&spjbuf[0][2]);
#else
			v8sf jbuf0, jbuf1, jbuf2;
#endif
			for(int j=0; j<nj; j++){
#ifndef PRELOAD_SPJ
				//jbuf0 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+0][0]);
                jbuf0 = _mm256_broadcast_ps((v4sf *)&spjbuf[j+0][0]);
#endif
				//v8sf xj  = __builtin_ia32_shufps256(jbuf0, jbuf0, 0x00);
				//v8sf yj  = __builtin_ia32_shufps256(jbuf0, jbuf0, 0x55);
				//v8sf zj  = __builtin_ia32_shufps256(jbuf0, jbuf0, 0xaa);
                v8sf xj =  _mm256_shuffle_ps(jbuf0, jbuf0, 0x00);
                v8sf yj =  _mm256_shuffle_ps(jbuf0, jbuf0, 0x55);
                v8sf zj =  _mm256_shuffle_ps(jbuf0, jbuf0, 0xaa);
#ifdef PRELOAD_SPJ
				//jbuf0 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+1][0]);
                jbuf0 = _mm256_broadcast_ps((v4sf *)&spjbuf[j+1][0]);
#endif

#ifndef PRELOAD_SPJ
				//jbuf1 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+0][1]);
                jbuf1 = _mm256_broadcast_ps((v4sf *)&spjbuf[j+0][1]);
#endif
				//v8sf qxx = __builtin_ia32_shufps256(jbuf1, jbuf1, 0x00);
				//v8sf qyy = __builtin_ia32_shufps256(jbuf1, jbuf1, 0x55);
				//v8sf qzz = __builtin_ia32_shufps256(jbuf1, jbuf1, 0xaa);
				//v8sf mj  = __builtin_ia32_shufps256(jbuf1, jbuf1, 0xff);
                v8sf qxx =  _mm256_shuffle_ps(jbuf1, jbuf1, 0x00);
                v8sf qyy =  _mm256_shuffle_ps(jbuf1, jbuf1, 0x55);
                v8sf qzz =  _mm256_shuffle_ps(jbuf1, jbuf1, 0xaa);
                v8sf mj  =  _mm256_shuffle_ps(jbuf1, jbuf1, 0xff);
#ifdef PRELOAD_SPJ
				//jbuf1 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+1][1]);
                jbuf1 = _mm256_broadcast_ps((v4sf *)&spjbuf[j+1][1]);
#endif

#ifndef PRELOAD_SPJ
				//jbuf2 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+0][2]);
                jbuf2 = _mm256_broadcast_ps((v4sf *)&spjbuf[j+0][2]);
#endif
				//v8sf qxy = __builtin_ia32_shufps256(jbuf2, jbuf2, 0x00);
				//v8sf qyz = __builtin_ia32_shufps256(jbuf2, jbuf2, 0x55);
				//v8sf qzx = __builtin_ia32_shufps256(jbuf2, jbuf2, 0xaa);
				//v8sf mtr = __builtin_ia32_shufps256(jbuf2, jbuf2, 0xff);
                v8sf qxy =  _mm256_shuffle_ps(jbuf2, jbuf2, 0x00);
                v8sf qyz =  _mm256_shuffle_ps(jbuf2, jbuf2, 0x55);
                v8sf qzx =  _mm256_shuffle_ps(jbuf2, jbuf2, 0xaa);
                v8sf mtr =  _mm256_shuffle_ps(jbuf2, jbuf2, 0xff);
#ifdef PRELOAD_SPJ
				//jbuf2 = __builtin_ia32_vbroadcastf128_ps256((v4sf *)&spjbuf[j+1][2]);
                jbuf2 = _mm256_broadcast_ps((v4sf *)&spjbuf[j+1][2]);
#endif

				v8sf dx = xj - xi;
				v8sf dy = yj - yi;
				v8sf dz = zj - zi;

				v8sf r2  = ((veps2 + dx*dx) + dy*dy) + dz*dz;
				//v8sf ri1 = __builtin_ia32_rsqrtps256(r2);
                v8sf ri1  = _mm256_rsqrt_ps(r2);
				v8sf v0p5 = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
#ifdef RSQRT_NR_SPJ_X2
                v8sf v3p0 = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#endif
				v8sf ri2 = ri1 * ri1;
				v8sf ri3 = ri1 * ri2;
				v8sf ri4 = ri2 * ri2;
				v8sf ri5 = ri2 * ri3;

				v8sf qr_x = (qxx*dx + qxy*dy) + qzx*dz;
				v8sf qr_y = (qyy*dy + qxy*dx) + qyz*dz;
				v8sf qr_z = (qzz*dz + qzx*dx) + qyz*dy;

				v8sf rqr = ((mtr + qr_x*dx) + qr_y*dy) + qr_z*dz;
				v8sf rqr_ri4 = rqr * ri4;

				//v8sf v0p5 = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};
				v8sf v2p5 = {2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f}; 

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

    __attribute__ ((noinline))
    void kernel_spj_64bit_nounroll(const int ni, const int nj){
        const v4df veps2 = {eps2, eps2, eps2, eps2};
        for(int i=0; i<ni; i+=4){
            const int il = i%8;
            const v4df xi = *(v4df *)(&xibufd[i/8][0][il]);
            const v4df yi = *(v4df *)(&xibufd[i/8][1][il]);
            const v4df zi = *(v4df *)(&xibufd[i/8][2][il]);

            v4df ax, ay, az, pot;
            ax = ay = az = pot = (v4df){0.0, 0.0, 0.0, 0.0};

            v4df jbuf0 = *((v4df*)spjbufd[0][0]);
            v4df jbuf1 = *((v4df*)spjbufd[0][1]);
            v4df jbuf2 = *((v4df*)spjbufd[0][2]);

            for(int j=0; j<nj; j++){
                v4df xj  = _mm256_permute4x64_pd(jbuf0, 0x00);
                v4df yj  = _mm256_permute4x64_pd(jbuf0, 0x55);
                v4df zj  = _mm256_permute4x64_pd(jbuf0, 0xaa);
                jbuf0 = *((v4df*)spjbufd[j+1][0]);

                v4df qxx = _mm256_permute4x64_pd(jbuf1, 0x00);
                v4df qyy = _mm256_permute4x64_pd(jbuf1, 0x55);
                v4df qzz = _mm256_permute4x64_pd(jbuf1, 0xaa);
                v4df mj  = _mm256_permute4x64_pd(jbuf1, 0xff);
                jbuf1 = *((v4df*)spjbufd[j+1][1]);
		
                v4df qxy = _mm256_permute4x64_pd(jbuf2, 0x00);
                v4df qyz = _mm256_permute4x64_pd(jbuf2, 0x55);
                v4df qzx = _mm256_permute4x64_pd(jbuf2, 0xaa);
                v4df mtr = _mm256_permute4x64_pd(jbuf2, 0xff);
                jbuf2 = *((v4df*)spjbufd[j+1][2]);
		
                v4df dx = xj - xi;
                v4df dy = yj - yi;
                v4df dz = zj - zi;
		
                v4df r2  = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                //v4df ri1 = __builtin_ia32_rsqrtps256(r2);
                //v4df ri1 = __builtin_ia32_cvtps2pd256( __builtin_ia32_rsqrtps( __builtin_ia32_cvtpd2ps256(r2)));
                v4df ri1  = _mm256_cvtps_pd( _mm_rsqrt_ps( _mm256_cvtpd_ps(r2)));
		
#ifdef RSQRT_NR_SPJ_X2
                //x2
                v4df v3p0 = {3.0, 3.0, 3.0, 3.0};
                ri1 *= (v3p0 - r2*(ri1*ri1));
#elif defined(RSQRT_NR_SPJ_X4)
                // x4
                v4df v8p0 = {8.0, 8.0, 8.0, 8.0};
                v4df v6p0 = {6.0, 6.0, 6.0, 6.0};
                v4df v5p0 = {5.0, 5.0, 5.0, 5.0};
                v4df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
                v4df v1p0 = {1.0, 1.0, 1.0, 1.0};
                v4df h = v1p0 - r2*(ri1*ri1);
                ri1 *= v1p0 + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
                v4df ri2 = ri1 * ri1;
                v4df ri3 = ri1 * ri2;
                v4df ri4 = ri2 * ri2;
                v4df ri5 = ri2 * ri3;
		
                v4df qr_x = (qxx*dx + qxy*dy) + qzx*dz;
                v4df qr_y = (qyy*dy + qxy*dx) + qyz*dz;
                v4df qr_z = (qzz*dz + qzx*dx) + qyz*dy;
		
                v4df rqr = ((mtr + qr_x*dx) + qr_y*dy) + qr_z*dz;
                v4df rqr_ri4 = rqr * ri4;
		
                v4df v0p5 = {0.5, 0.5, 0.5, 0.5};
                v4df v2p5 = {2.5, 2.5, 2.5, 2.5};
		
                v4df meff  =  mj + v0p5 * rqr_ri4;
                v4df meff3 = (mj + v2p5 * rqr_ri4) * ri3;

                pot -= meff * ri1;

                ax = (ax - ri5*qr_x) + meff3*dx;
                ay = (ay - ri5*qr_y) + meff3*dy;
                az = (az - ri5*qr_z) + meff3*dz;
            }

#ifdef RSQRT_NR_SPJ_X2
            //x2
            v4df v0p5 = {0.5, 0.5, 0.5, 0.5};
            v4df v0p125 = {0.125, 0.125, 0.125, 0.125};
            pot *= v0p5;
            ax  *= v0p125;
            ay  *= v0p125;
            az  *= v0p125;
#endif
	    
            *(v4df *)(&accpbufd[i/8][0][il]) = ax;
            *(v4df *)(&accpbufd[i/8][1][il]) = ay;
            *(v4df *)(&accpbufd[i/8][2][il]) = az;
            *(v4df *)(&accpbufd[i/8][3][il]) = pot;
        }
    }
    

#elif defined(__AVX512DQ__)
    //////////////////////
    ///   For AVX512   ///
    //////////////////////

    typedef float v16sf __attribute__((vector_size(64)));
    typedef double v8df __attribute__((vector_size(64)));
    
    __attribute__ ((noinline))
	void kernel_epj_nounroll(const int ni, const int nj){
        
	    const v16sf veps2 = {(float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2,
                             (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2};
        const v8sf vzero = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        
		for(int i=0; i<ni; i+=8){
            const v16sf xi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][0]));
			const v16sf yi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][1]));
			const v16sf zi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][2]));
            
			v16sf ax, ay, az, pot;
			ax = ay = az = pot = (v16sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                                         0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
#ifdef USE_INDIVIDUAL_CUTOFF
            v16sf jbuf = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&epjbuf[1][0]),
                                            _mm256_broadcast_ps((v4sf *)&epjbuf[0][0]), 0x00);
#else
            v16sf jbuf = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)(epjbuf+1)),
                                            _mm256_broadcast_ps((v4sf *)(epjbuf+0)), 0x00);
#endif      
            v16sf xj =  _mm512_shuffle_ps(jbuf, jbuf, 0x00);
            v16sf yj =  _mm512_shuffle_ps(jbuf, jbuf, 0x55);
            v16sf zj =  _mm512_shuffle_ps(jbuf, jbuf, 0xaa);
            v16sf mj =  _mm512_shuffle_ps(jbuf, jbuf, 0xff);
            if ( nj<=1 ) mj = _mm512_insertf32x8(mj, vzero, 0x01);
			for(int j=0; j<nj; j+=2){
#ifdef USE_INDIVIDUAL_CUTOFF
                jbuf = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&epjbuf[j+3][0]),
                                          _mm256_broadcast_ps((v4sf *)&epjbuf[j+2][0]), 0x00);
#else
                jbuf = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)(epjbuf+3)),
                                          _mm256_broadcast_ps((v4sf *)(epjbuf+2)), 0x00);
#endif

				v16sf dx = xj - xi;
				v16sf dy = yj - yi;
				v16sf dz = zj - zi;

				v16sf r2   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v16sf ri1  = _mm512_rsqrt14_ps(r2);
#ifdef RSQRT_NR_EPJ_X2
				v16sf v3p0 = {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f,
                              3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f}; 
				ri1 *= (v3p0 - r2*(ri1*ri1));
#endif
				v16sf mri1 = mj * ri1;
				v16sf ri2  = ri1 * ri1;
				v16sf mri3 = mri1 * ri2;

                xj =  _mm512_shuffle_ps(jbuf, jbuf, 0x00);
                yj =  _mm512_shuffle_ps(jbuf, jbuf, 0x55);
                zj =  _mm512_shuffle_ps(jbuf, jbuf, 0xaa);
                mj =  _mm512_shuffle_ps(jbuf, jbuf, 0xff);
                if ( nj<=j+3 ) mj = _mm512_insertf32x8(mj, vzero, 0x01);

				pot -= mri1;
				ax += mri3 * dx;
				ay += mri3 * dy;
				az += mri3 * dz;
			}
#ifdef RSQRT_NR_EPJ_X2
			v16sf v0p5 = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
                          0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f}; 
			v16sf v0p125 = {0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f,
                            0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f}; 
			pot *= v0p5;
			ax  *= v0p125;
			ay  *= v0p125;
			az  *= v0p125;
#endif

			*(v8sf *)(accpbuf[i/8][0]) = _mm512_extractf32x8_ps(ax,  0x00) + _mm512_extractf32x8_ps(ax,  0x01);
			*(v8sf *)(accpbuf[i/8][1]) = _mm512_extractf32x8_ps(ay,  0x00) + _mm512_extractf32x8_ps(ay,  0x01);
			*(v8sf *)(accpbuf[i/8][2]) = _mm512_extractf32x8_ps(az,  0x00) + _mm512_extractf32x8_ps(az,  0x01);
			*(v8sf *)(accpbuf[i/8][3]) = _mm512_extractf32x8_ps(pot, 0x00) + _mm512_extractf32x8_ps(pot, 0x01);
		}
	}

    __attribute__ ((noinline))
    void kernel_epj_nounroll_64bit(const int ni, const int nj){

        const v8df veps2 = {eps2, eps2, eps2, eps2, eps2, eps2, eps2, eps2};
        const v4df vzero = (v4df){0.0, 0.0, 0.0, 0.0};
        
        for(int i=0; i<ni; i+=4){
            const int il = i%8;
            const v8df xi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][0][il]));
            const v8df yi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][1][il]));
            const v8df zi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][2][il]));

            v8df ax, ay, az, pot;
            ax = ay = az = pot = (v8df){0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

#ifdef USE_INDIVIDUAL_CUTOFF
            v8df jbuf = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&epjbufd[1][0])),
                                           *((v4df *)&epjbufd[0][0]), 0x00);
#else
            v8df jbuf = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df*)(epjbufd+1))),
                                           *(v4df*)(epjbufd+0), 0x00);
#endif
            v8df xj =  _mm512_permutex_pd(jbuf, 0x00);
            v8df yj =  _mm512_permutex_pd(jbuf, 0x55);
            v8df zj =  _mm512_permutex_pd(jbuf, 0xaa);
            v8df mj =  _mm512_permutex_pd(jbuf, 0xff);
            if ( nj<=1 ) mj = _mm512_insertf64x4(mj, vzero, 0x01);

            for(int j=0; j<nj; j+=2){
#ifdef USE_INDIVIDUAL_CUTOFF
                jbuf = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&epjbufd[j+3][0])),
                                          *((v4df *)&epjbufd[j+2][0]), 0x00);
#else
                jbuf = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df*)(epjbufd+3))),
                                          *(v4df*)(epjbufd+2), 0x00);
#endif
                v8df dx = xj - xi;
                v8df dy = yj - yi;
                v8df dz = zj - zi;
                v8df r2   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v8df ri1  = _mm512_cvtps_pd( _mm256_rsqrt_ps( _mm512_cvtpd_ps(r2)));
#ifdef RSQRT_NR_EPJ_X2
                //x2
                v8df v0p5 = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
                v8df v3p0 = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
                // x4
                v8df vone = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
                v8df v8p0 = {8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0};
                v8df v6p0 = {6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0};
                v8df v5p0 = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
                v8df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
                v8df h = vone - r2*(ri1*ri1);
                ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
                v8df mri1 = mj * ri1;
                v8df ri2  = ri1 * ri1;
                v8df mri3 = mri1 * ri2;
		
                xj = _mm512_permutex_pd(jbuf, 0x00);
                yj = _mm512_permutex_pd(jbuf, 0x55);
                zj = _mm512_permutex_pd(jbuf, 0xaa);
                mj = _mm512_permutex_pd(jbuf, 0xff);
                if ( nj<=j+3 ) mj = _mm512_insertf64x4(mj, vzero, 0x01);

                pot -= mri1;
                ax += mri3 * dx;
                ay += mri3 * dy;
                az += mri3 * dz;
            }
            *(v4df *)(&accpbufd[i/8][0][il]) = _mm512_extractf64x4_pd(ax,  0x00) + _mm512_extractf64x4_pd(ax,  0x01);
			*(v4df *)(&accpbufd[i/8][1][il]) = _mm512_extractf64x4_pd(ay,  0x00) + _mm512_extractf64x4_pd(ay,  0x01);
			*(v4df *)(&accpbufd[i/8][2][il]) = _mm512_extractf64x4_pd(az,  0x00) + _mm512_extractf64x4_pd(az,  0x01);
			*(v4df *)(&accpbufd[i/8][3][il]) = _mm512_extractf64x4_pd(pot, 0x00) + _mm512_extractf64x4_pd(pot, 0x01);
        }
    }

    ///////////////////////
    // with linear cutoff
    __attribute__ ((noinline))
    void kernel_epj_nounroll_for_p3t_with_linear_cutoff(const int ni, const int nj){
        
        const v16sf vone = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                            1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
        const v16sf veps2 = {(float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2,
                             (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2};
#ifndef USE_INDIVIDUAL_CUTOFF
        const v16sf vr_out  = {(float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,
                              (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out,  (float)r_out};
        const v16sf vr_out2  = vr_out * vr_out;
#endif
        const v8sf vzero = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        
        for(int i=0; i<ni; i+=8){
            const v16sf xi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][0]));
            const v16sf yi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][1]));
            const v16sf zi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][2]));
#ifdef USE_INDIVIDUAL_CUTOFF
            const v16sf r_outi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][3]));
#endif
            
            v16sf ax, ay, az, pot;
            ax = ay = az = pot = (v16sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                                         0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
#ifdef USE_INDIVIDUAL_CUTOFF
            v16sf jbuf0 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&epjbuf[1][0]),
                                             _mm256_broadcast_ps((v4sf *)&epjbuf[0][0]), 0x00);
            v16sf jbuf1 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&epjbuf[1][1]),
                                             _mm256_broadcast_ps((v4sf *)&epjbuf[0][1]), 0x00);
            v16sf xj     =  _mm512_shuffle_ps(jbuf0, jbuf0, 0x00);
            v16sf yj     =  _mm512_shuffle_ps(jbuf0, jbuf0, 0x55);
            v16sf zj     =  _mm512_shuffle_ps(jbuf0, jbuf0, 0xaa);
            v16sf mj     =  _mm512_shuffle_ps(jbuf0, jbuf0, 0xff);
            v16sf r_outj =  _mm512_shuffle_ps(jbuf1, jbuf1, 0x00);
#else
            v16sf jbuf = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)(epjbuf+1)),
                                            _mm256_broadcast_ps((v4sf *)(epjbuf+0)), 0x00);
            v16sf xj    =  _mm512_shuffle_ps(jbuf, jbuf, 0x00);
            v16sf yj    =  _mm512_shuffle_ps(jbuf, jbuf, 0x55);
            v16sf zj    =  _mm512_shuffle_ps(jbuf, jbuf, 0xaa);
            v16sf mj    =  _mm512_shuffle_ps(jbuf, jbuf, 0xff);
#endif
            if ( nj<=1 ) mj = _mm512_insertf32x8(mj, vzero, 0x01);
            
            for(int j=0; j<nj; j+=2){
#ifdef USE_INDIVIDUAL_CUTOFF
                jbuf0 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&epjbuf[j+3][0]),
                                           _mm256_broadcast_ps((v4sf *)&epjbuf[j+2][0]), 0x00);
                jbuf1 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&epjbuf[j+3][1]),
                                           _mm256_broadcast_ps((v4sf *)&epjbuf[j+2][1]), 0x00);
                v16sf vr_out  = _mm512_max_ps(r_outi, r_outj);
                v16sf vr_out2 = vr_out * vr_out;
#else
                jbuf = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)(epjbuf+j+3)),
                                          _mm256_broadcast_ps((v4sf *)(epjbuf+j+2)), 0x00);
#endif
                v16sf dx = xj - xi;
                v16sf dy = yj - yi;
                v16sf dz = zj - zi;
                v16sf r2_real   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v16sf r2 = _mm512_max_ps(r2_real, vr_out2);
                v16sf ri1  = _mm512_rsqrt14_ps(r2);
#ifdef RSQRT_NR_EPJ_X2
                v16sf v3p0 = {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f,
                              3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f}; 
                v16sf v0p5 = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
                              0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
                // x4
                v16sf v8p0 = {8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f,
                              8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f};
                v16sf v6p0 = {6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f,
                              6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f};
                v16sf v5p0 = {5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f,
                              5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f};
                v16sf v0p0625 = {(float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0,
                                 (float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0};
                v16sf h = vone - r2*(ri1*ri1);
                ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
                v16sf ri2 = ri1*ri1;
                v16sf mri1 = mj*ri1;
                v16sf mri3 = mri1 * ri2;
#ifdef USE_INDIVIDUAL_CUTOFF
                xj     = _mm512_shuffle_ps(jbuf0, jbuf0, 0x00);
                yj     = _mm512_shuffle_ps(jbuf0, jbuf0, 0x55);
                zj     = _mm512_shuffle_ps(jbuf0, jbuf0, 0xaa);
                mj     = _mm512_shuffle_ps(jbuf0, jbuf0, 0xff);
                r_outj = _mm512_shuffle_ps(jbuf1, jbuf1, 0x00);
#else
                xj = _mm512_shuffle_ps(jbuf, jbuf, 0x00);
                yj = _mm512_shuffle_ps(jbuf, jbuf, 0x55);
                zj = _mm512_shuffle_ps(jbuf, jbuf, 0xaa);
                mj = _mm512_shuffle_ps(jbuf, jbuf, 0xff);
#endif
                if ( nj<=j+3 ) mj = _mm512_insertf32x8(mj, vzero, 0x01);
                
                pot -= mri1;
                ax += mri3 * dx;
                ay += mri3 * dy;
                az += mri3 * dz;
            }
            *(v8sf *)(accpbuf[i/8][0]) = _mm512_extractf32x8_ps(ax,  0x00) + _mm512_extractf32x8_ps(ax,  0x01);
			*(v8sf *)(accpbuf[i/8][1]) = _mm512_extractf32x8_ps(ay,  0x00) + _mm512_extractf32x8_ps(ay,  0x01);
			*(v8sf *)(accpbuf[i/8][2]) = _mm512_extractf32x8_ps(az,  0x00) + _mm512_extractf32x8_ps(az,  0x01);
			*(v8sf *)(accpbuf[i/8][3]) = _mm512_extractf32x8_ps(pot, 0x00) + _mm512_extractf32x8_ps(pot, 0x01);
        }
    }

    __attribute__ ((noinline))
    void kernel_epj_64bit_nounroll_for_p3t_with_linear_cutoff(const int ni, const int nj){
        const v8df vone = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        const v8df veps2 = {eps2, eps2, eps2, eps2, eps2, eps2, eps2, eps2};
        const v4df vzero = (v4df){0.0, 0.0, 0.0, 0.0};
#ifndef USE_INDIVIDUAL_CUTOFF
        const v8df vr_out  = {r_out, r_out, r_out, r_out, r_out, r_out, r_out, r_out};
        const v8df vr_out2  = vr_out * vr_out;
#endif
        for(int i=0; i<ni; i+=4){
            const int il = i%8;
            const v8df xi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][0][il]));
            const v8df yi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][1][il]));
            const v8df zi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][2][il]));
#ifdef USE_INDIVIDUAL_CUTOFF
            const v8df r_outi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][3][il]));
#endif
            
            v8df ax, ay, az, pot;
            ax = ay = az = pot = (v8df){0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#ifdef USE_INDIVIDUAL_CUTOFF
            v8df jbuf0 = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&epjbufd[1][0])),
                                            *((v4df *)&epjbufd[0][0]), 0x00);
            v8df jbuf1 = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&epjbufd[1][1])),
                                            *((v4df *)&epjbufd[0][1]), 0x00);
            v8df xj     =  _mm512_permutex_pd(jbuf0, 0x00);
            v8df yj     =  _mm512_permutex_pd(jbuf0, 0x55);
            v8df zj     =  _mm512_permutex_pd(jbuf0, 0xaa);
            v8df mj     =  _mm512_permutex_pd(jbuf0, 0xff);
            v8df r_outj =  _mm512_permutex_pd(jbuf1, 0x00);
#else
            v8df jbuf = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df*)(epjbufd+1))),
                                           *(v4df*)(epjbufd+0), 0x00);
            v8df xj =  _mm512_permutex_pd(jbuf, 0x00);
            v8df yj =  _mm512_permutex_pd(jbuf, 0x55);
            v8df zj =  _mm512_permutex_pd(jbuf, 0xaa);
            v8df mj =  _mm512_permutex_pd(jbuf, 0xff);
#endif
            if ( nj<=1 ) mj = _mm512_insertf64x4(mj, vzero, 0x01);
            
            for(int j=0; j<nj; j+=2){
#ifdef USE_INDIVIDUAL_CUTOFF
                jbuf0 = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&epjbufd[j+3][0])),
                                           *((v4df *)&epjbufd[j+2][0]), 0x00);
                jbuf1 = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&epjbufd[j+3][1])),
                                           *((v4df *)&epjbufd[j+2][1]), 0x00);
                v8df vr_out  = _mm512_max_pd(r_outi, r_outj);
                v8df vr_out2 = vr_out * vr_out;
#else
                jbuf = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df*)(epjbufd+j+3))),
                                          *(v4df*)(epjbufd+j+2), 0x00);
#endif
                v8df dx = xj - xi;
                v8df dy = yj - yi;
                v8df dz = zj - zi;
                v8df r2_real   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v8df r2  = _mm512_max_pd( r2_real, vr_out2);
                v8df ri1 = _mm512_cvtps_pd( _mm256_rsqrt_ps( _mm512_cvtpd_ps(r2)));
#ifdef RSQRT_NR_EPJ_X2
                //x2
                v8df v0p5 = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
                v8df v3p0 = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
                // x4
                v8df v8p0 = {8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0};
                v8df v6p0 = {6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0};
                v8df v5p0 = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
                v8df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
                v8df h = vone - r2*(ri1*ri1);
                ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
                v8df ri2 = ri1*ri1;
                v8df mri1 = mj*ri1;
                v8df mri3 = mri1 * ri2;
#ifdef USE_INDIVIDUAL_CUTOFF
                xj     = _mm512_permutex_pd(jbuf0, 0x00);
                yj     = _mm512_permutex_pd(jbuf0, 0x55);
                zj     = _mm512_permutex_pd(jbuf0, 0xaa);
                mj     = _mm512_permutex_pd(jbuf0, 0xff);
                r_outj = _mm512_permutex_pd(jbuf1, 0x00);
#else
                xj = _mm512_permutex_pd(jbuf, 0x00);
                yj = _mm512_permutex_pd(jbuf, 0x55);
                zj = _mm512_permutex_pd(jbuf, 0xaa);
                mj = _mm512_permutex_pd(jbuf, 0xff);
#endif
                if ( nj<=j+3 ) mj = _mm512_insertf64x4(mj, vzero, 0x01);
                                
                pot -= mri1;
                ax += mri3 * dx;
                ay += mri3 * dy;
                az += mri3 * dz;
            }
            *(v4df *)(&accpbufd[i/8][0][il]) = _mm512_extractf64x4_pd(ax,  0x00) + _mm512_extractf64x4_pd(ax,  0x01);
			*(v4df *)(&accpbufd[i/8][1][il]) = _mm512_extractf64x4_pd(ay,  0x00) + _mm512_extractf64x4_pd(ay,  0x01);
			*(v4df *)(&accpbufd[i/8][2][il]) = _mm512_extractf64x4_pd(az,  0x00) + _mm512_extractf64x4_pd(az,  0x01);
			*(v4df *)(&accpbufd[i/8][3][il]) = _mm512_extractf64x4_pd(pot, 0x00) + _mm512_extractf64x4_pd(pot, 0x01);
        }
    }
    // linear cutoff
    //////////////

    __attribute__ ((noinline))
	void kernel_spj_nounroll(const int ni, const int nj){
        
	    const v16sf veps2 = {(float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2,
                             (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2};
        const v8sf vzero = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        
		for(int i=0; i<ni; i+=8){
            const v16sf xi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][0]));
			const v16sf yi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][1]));
			const v16sf zi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][2]));

			v16sf ax, ay, az, pot;
			ax = ay = az = pot = (v16sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                                         0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

#define PRELOAD_SPJ

#ifdef PRELOAD_SPJ
            v16sf jbuf0 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&spjbuf[1][0]),
                                             _mm256_broadcast_ps((v4sf *)&spjbuf[0][0]), 0x00);
            v16sf jbuf1 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&spjbuf[1][1]),
                                             _mm256_broadcast_ps((v4sf *)&spjbuf[0][1]), 0x00);
            v16sf jbuf2 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&spjbuf[1][2]),
                                             _mm256_broadcast_ps((v4sf *)&spjbuf[0][2]), 0x00);
#else
			v16sf jbuf0, jbuf1, jbuf2;
#endif
			for(int j=0; j<nj; j+=2){
#ifndef PRELOAD_SPJ
                jbuf0 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&spjbuf[j+1][0]),
                                           _mm256_broadcast_ps((v4sf *)&spjbuf[j+0][0]), 0x00);
#endif
                v16sf xj =  _mm512_shuffle_ps(jbuf0, jbuf0, 0x00);
                v16sf yj =  _mm512_shuffle_ps(jbuf0, jbuf0, 0x55);
                v16sf zj =  _mm512_shuffle_ps(jbuf0, jbuf0, 0xaa);
#ifdef PRELOAD_SPJ
                jbuf0 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&spjbuf[j+3][0]),
                                           _mm256_broadcast_ps((v4sf *)&spjbuf[j+2][0]), 0x00);
#endif

#ifndef PRELOAD_SPJ
                jbuf1 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&spjbuf[j+1][1]),
                                           _mm256_broadcast_ps((v4sf *)&spjbuf[j+0][1]), 0x00);
#endif
                v16sf qxx =  _mm512_shuffle_ps(jbuf1, jbuf1, 0x00);
                v16sf qyy =  _mm512_shuffle_ps(jbuf1, jbuf1, 0x55);
                v16sf qzz =  _mm512_shuffle_ps(jbuf1, jbuf1, 0xaa);
                v16sf mj  =  _mm512_shuffle_ps(jbuf1, jbuf1, 0xff);
#ifdef PRELOAD_SPJ
                jbuf1 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&spjbuf[j+3][1]),
                                           _mm256_broadcast_ps((v4sf *)&spjbuf[j+2][1]), 0x00);
#endif

#ifndef PRELOAD_SPJ
                jbuf2 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&spjbuf[j+1][2]),
                                           _mm256_broadcast_ps((v4sf *)&spjbuf[j+0][2]), 0x00);
#endif
                v16sf qxy =  _mm512_shuffle_ps(jbuf2, jbuf2, 0x00);
                v16sf qyz =  _mm512_shuffle_ps(jbuf2, jbuf2, 0x55);
                v16sf qzx =  _mm512_shuffle_ps(jbuf2, jbuf2, 0xaa);
                v16sf mtr =  _mm512_shuffle_ps(jbuf2, jbuf2, 0xff);
                if ( nj<=j+1 ) {
                    qxx = _mm512_insertf32x8(qxx, vzero, 0x01);
                    qyy = _mm512_insertf32x8(qyy, vzero, 0x01);
                    qzz = _mm512_insertf32x8(qzz, vzero, 0x01);
                    mj  = _mm512_insertf32x8(mj,  vzero, 0x01);
                    qxy = _mm512_insertf32x8(qxy, vzero, 0x01);
                    qyz = _mm512_insertf32x8(qyz, vzero, 0x01);
                    qzx = _mm512_insertf32x8(qzx, vzero, 0x01);
                    mtr = _mm512_insertf32x8(mtr, vzero, 0x01);
                }
#ifdef PRELOAD_SPJ
                jbuf2 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&spjbuf[j+3][2]),
                                           _mm256_broadcast_ps((v4sf *)&spjbuf[j+2][2]), 0x00);
#endif

				v16sf dx = xj - xi;
				v16sf dy = yj - yi;
				v16sf dz = zj - zi;

				v16sf r2  = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v16sf ri1  = _mm512_rsqrt14_ps(r2);
				v16sf v0p5 = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                              0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
#ifdef RSQRT_NR_SPJ_X2
                v16sf v3p0 = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
                              3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#endif
				v16sf ri2 = ri1 * ri1;
				v16sf ri3 = ri1 * ri2;
				v16sf ri4 = ri2 * ri2;
				v16sf ri5 = ri2 * ri3;

				v16sf qr_x = (qxx*dx + qxy*dy) + qzx*dz;
				v16sf qr_y = (qyy*dy + qxy*dx) + qyz*dz;
				v16sf qr_z = (qzz*dz + qzx*dx) + qyz*dy;

				v16sf rqr = ((mtr + qr_x*dx) + qr_y*dy) + qr_z*dz;
				v16sf rqr_ri4 = rqr * ri4;

				v16sf v2p5 = {2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f,
                              2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f}; 

				v16sf meff  =  mj + v0p5 * rqr_ri4;
				v16sf meff3 = (mj + v2p5 * rqr_ri4) * ri3;

				pot -= meff * ri1;

				ax = (ax - ri5*qr_x) + meff3*dx;
				ay = (ay - ri5*qr_y) + meff3*dy;
				az = (az - ri5*qr_z) + meff3*dz;
			}
            *(v8sf *)(accpbuf[i/8][0]) = _mm512_extractf32x8_ps(ax,  0x00) + _mm512_extractf32x8_ps(ax,  0x01);
			*(v8sf *)(accpbuf[i/8][1]) = _mm512_extractf32x8_ps(ay,  0x00) + _mm512_extractf32x8_ps(ay,  0x01);
			*(v8sf *)(accpbuf[i/8][2]) = _mm512_extractf32x8_ps(az,  0x00) + _mm512_extractf32x8_ps(az,  0x01);
			*(v8sf *)(accpbuf[i/8][3]) = _mm512_extractf32x8_ps(pot, 0x00) + _mm512_extractf32x8_ps(pot, 0x01);
		}
	}

    __attribute__ ((noinline))
    void kernel_spj_64bit_nounroll(const int ni, const int nj){
        const v8df veps2 = {eps2, eps2, eps2, eps2, eps2, eps2, eps2, eps2};
        const v4df vzero = (v4df){0.0, 0.0, 0.0, 0.0};
        
        for(int i=0; i<ni; i+=4){
            const int il = i%8;
            const v8df xi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][0][il]));
            const v8df yi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][1][il]));
            const v8df zi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][2][il]));

            v8df ax, ay, az, pot;
            ax = ay = az = pot = (v8df){0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

            v8df jbuf0 = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&spjbufd[1][0])),
                                            *((v4df *)&spjbufd[0][0]), 0x00);
            v8df jbuf1 = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&spjbufd[1][1])),
                                            *((v4df *)&spjbufd[0][1]), 0x00);
            v8df jbuf2 = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&spjbufd[1][2])),
                                            *((v4df *)&spjbufd[0][2]), 0x00);

            for(int j=0; j<nj; j+=2){
                v8df xj =  _mm512_permutex_pd(jbuf0, 0x00);
                v8df yj =  _mm512_permutex_pd(jbuf0, 0x55);
                v8df zj =  _mm512_permutex_pd(jbuf0, 0xaa);
                jbuf0 = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&spjbufd[j+3][0])),
                                           *((v4df *)&spjbufd[j+2][0]), 0x00);

                v8df qxx = _mm512_permutex_pd(jbuf1, 0x00);
                v8df qyy = _mm512_permutex_pd(jbuf1, 0x55);
                v8df qzz = _mm512_permutex_pd(jbuf1, 0xaa);
                v8df mj  = _mm512_permutex_pd(jbuf1, 0xff);
                jbuf1 = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&spjbufd[j+3][1])),
                                           *((v4df *)&spjbufd[j+2][1]), 0x00);
		
                v8df qxy = _mm512_permutex_pd(jbuf2, 0x00);
                v8df qyz = _mm512_permutex_pd(jbuf2, 0x55);
                v8df qzx = _mm512_permutex_pd(jbuf2, 0xaa);
                v8df mtr = _mm512_permutex_pd(jbuf2, 0xff);
                jbuf2 = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&spjbufd[j+3][2])),
                                           *((v4df *)&spjbufd[j+2][2]), 0x00);
                if ( nj<=j+1 ) {
                    qxx = _mm512_insertf64x4(qxx, vzero, 0x01);
                    qyy = _mm512_insertf64x4(qyy, vzero, 0x01);
                    qzz = _mm512_insertf64x4(qzz, vzero, 0x01);
                    mj  = _mm512_insertf64x4(mj,  vzero, 0x01);
                    qxy = _mm512_insertf64x4(qxy, vzero, 0x01);
                    qyz = _mm512_insertf64x4(qyz, vzero, 0x01);
                    qzx = _mm512_insertf64x4(qzx, vzero, 0x01);
                    mtr = _mm512_insertf64x4(mtr, vzero, 0x01);
                }
		
                v8df dx = xj - xi;
                v8df dy = yj - yi;
                v8df dz = zj - zi;
		
                v8df r2  = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v8df ri1  = _mm512_cvtps_pd( _mm256_rsqrt_ps( _mm512_cvtpd_ps(r2)));
		
#ifdef RSQRT_NR_SPJ_X2
                //x2
                v8df v3p0 = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
                ri1 *= (v3p0 - r2*(ri1*ri1));
#elif defined(RSQRT_NR_SPJ_X4)
                // x4
                v8df v8p0 = {8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0};
                v8df v6p0 = {6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0};
                v8df v5p0 = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
                v8df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
                v8df v1p0 = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
                v8df h = v1p0 - r2*(ri1*ri1);
                ri1 *= v1p0 + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
                v8df ri2 = ri1 * ri1;
                v8df ri3 = ri1 * ri2;
                v8df ri4 = ri2 * ri2;
                v8df ri5 = ri2 * ri3;
		
                v8df qr_x = (qxx*dx + qxy*dy) + qzx*dz;
                v8df qr_y = (qyy*dy + qxy*dx) + qyz*dz;
                v8df qr_z = (qzz*dz + qzx*dx) + qyz*dy;
		
                v8df rqr = ((mtr + qr_x*dx) + qr_y*dy) + qr_z*dz;
                v8df rqr_ri4 = rqr * ri4;
		
                v8df v0p5 = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
                v8df v2p5 = {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5};
		
                v8df meff  =  mj + v0p5 * rqr_ri4;
                v8df meff3 = (mj + v2p5 * rqr_ri4) * ri3;

                pot -= meff * ri1;

                ax = (ax - ri5*qr_x) + meff3*dx;
                ay = (ay - ri5*qr_y) + meff3*dy;
                az = (az - ri5*qr_z) + meff3*dz;
            }

#ifdef RSQRT_NR_SPJ_X2
            //x2
            v8df v0p5 = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
            v8df v0p125 = {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125};
            pot *= v0p5;
            ax  *= v0p125;
            ay  *= v0p125;
            az  *= v0p125;
#endif

            *(v4df *)(&accpbufd[i/8][0][il]) = _mm512_extractf64x4_pd(ax,  0x00) + _mm512_extractf64x4_pd(ax,  0x01);
			*(v4df *)(&accpbufd[i/8][1][il]) = _mm512_extractf64x4_pd(ay,  0x00) + _mm512_extractf64x4_pd(ay,  0x01);
			*(v4df *)(&accpbufd[i/8][2][il]) = _mm512_extractf64x4_pd(az,  0x00) + _mm512_extractf64x4_pd(az,  0x01);
			*(v4df *)(&accpbufd[i/8][3][il]) = _mm512_extractf64x4_pd(pot, 0x00) + _mm512_extractf64x4_pd(pot, 0x01);
        }
    }

#endif //__AVX512DQ__


} __attribute__ ((aligned(128)));







