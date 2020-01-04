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

#define SAFETY_FACTOR 1.01


class PhantomGrapeQuad{
    
public:
    enum{
        //NIMAX = 32768,
        //NJMAX = 131072,
        NIMAX = 512,
        NJMAX = 2048,
    };
    
private:
    float  xibuf   [NIMAX/8]  [5][8];   // x, y, z, r_out, r_search
    double xibufd  [NIMAX/8]  [5][8];   // x, y, z, r_out, r_search
    float  epjbuf  [NJMAX]    [2][4];   // x, y, z, m, | r_out, r_search, idx
    double epjbufd [NJMAX]    [2][4];   // x, y, z, m, | r_out, r_search, idx
    float  spjbuf  [NJMAX]    [3][4];   // x, y, z, m, | xx, yy, zz, pad, | xy, yz, zx, tr
    double spjbufd [NJMAX]    [3][4];   // x, y, z, m, | xx, yy, zz, pad, | xy, yz, zx, tr
    float  accpbuf [NIMAX/8]  [6][8];   // ax, ay, az, pot, nngb, idxngb
    double accpbufd[NIMAX/8]  [6][8];   // ax, ay, az, pot, nngb, idxngb


    double eps2;
    double gamma;
    static double get_a_NaN(){
        union{ long   l; double d; } m;
        m.l = -1;
        return m.d;
    }
    double r_out;
    //double r_in;
    double r_search;
    double denominator; // for cut off
public:

    PhantomGrapeQuad() : eps2(get_a_NaN()) {} // default NaN

    void set_cutoff(const double _r_out, const double _r_search){
        r_out = _r_out;
        //r_in = _r_in;
        r_search = SAFETY_FACTOR * _r_search;
        denominator = 1.0 / ((1.-gamma)*r_out);
    }

    void set_eps2(const double _eps2, const double _gamma){
        this->eps2 = _eps2;
        this->gamma = _gamma;
    }

    void set_epj_one(const int addr, const double x, const double y, const double z, const double m,
#ifdef USE_INDIVIDUAL_CUTOFF
                     const double _r_out, const double _r_search,
#endif
                     const double idx) {
        epjbuf[addr][0][0] = x;
        epjbuf[addr][0][1] = y;
        epjbuf[addr][0][2] = z;
        epjbuf[addr][0][3] = m;

#ifdef USE_INDIVIDUAL_CUTOFF
        epjbuf[addr][1][0] = _r_out;
        epjbuf[addr][1][1] = SAFETY_FACTOR * _r_search;
#else
        epjbuf[addr][1][0] = r_out;
        epjbuf[addr][1][1] = SAFETY_FACTOR * r_search;
#endif
        epjbuf[addr][1][2] = idx;
        epjbuf[addr][1][3] = 0.;
    }
    void set_epj_one_d(const int addr, const double x, const double y, const double z, const double m,
#ifdef USE_INDIVIDUAL_CUTOFF
                       const double _r_out, const double _r_search,
#endif
                       const double idx){
        epjbufd[addr][0][0] = x;
        epjbufd[addr][0][1] = y;
        epjbufd[addr][0][2] = z;
        epjbufd[addr][0][3] = m;

#ifdef USE_INDIVIDUAL_CUTOFF
        epjbufd[addr][1][0] = _r_out;
        epjbufd[addr][1][1] = SAFETY_FACTOR * _r_search;
#else
        epjbufd[addr][1][0] = r_out;
        epjbufd[addr][1][1] = SAFETY_FACTOR * r_search;
#endif
        epjbufd[addr][1][2] = idx;
        epjbufd[addr][1][3] = 0.;
    }

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
    
    void copyz_epj_one(const int addr, const int addr0)
    {
        epjbuf[addr][0][0] = epjbuf[addr0][0][0];
        epjbuf[addr][0][1] = epjbuf[addr0][0][1];
        epjbuf[addr][0][2] = epjbuf[addr0][0][2];
        epjbuf[addr][0][3] = 0.;
        
        epjbuf[addr][1][0] = 0.;
        epjbuf[addr][1][1] = 0.;
        
        epjbuf[addr][1][2] = -1;
        epjbuf[addr][1][3] = 0.;
    }
    void copyz_epj_one_d(const int addr, const int addr0)
    {
        epjbufd[addr][0][0] = epjbufd[addr0][0][0];
        epjbufd[addr][0][1] = epjbufd[addr0][0][1];
        epjbufd[addr][0][2] = epjbufd[addr0][0][2];
        epjbufd[addr][0][3] = 0.;
        
        epjbufd[addr][1][0] = 0.;
        epjbufd[addr][1][1] = 0.;
        
        epjbufd[addr][1][2] = -1;
        epjbufd[addr][1][3] = 0.;
    }
    void copyz_spj_one(const int addr, const int addr0)
    {
        spjbuf[addr][0][0] = spjbuf[addr0][0][0];
        spjbuf[addr][0][1] = spjbuf[addr0][0][1];
        spjbuf[addr][0][2] = spjbuf[addr0][0][2];
        spjbuf[addr][0][3] = 0.;

        spjbuf[addr][1][0] = 0.;
        spjbuf[addr][1][1] = 0.;
        spjbuf[addr][1][2] = 0.;
        spjbuf[addr][1][3] = 0.;
        
        spjbuf[addr][2][0] = 0.;
        spjbuf[addr][2][1] = 0.;
        spjbuf[addr][2][2] = 0.;
        spjbuf[addr][2][3] = 0.;
    }
    void copyz_spj_one_d(const int addr, const int addr0)
    {
        spjbufd[addr][0][0] = spjbufd[addr0][0][0];
        spjbufd[addr][0][1] = spjbufd[addr0][0][1];
        spjbufd[addr][0][2] = spjbufd[addr0][0][2];
        spjbufd[addr][0][3] = 0.;

        spjbufd[addr][1][0] = 0.;
        spjbufd[addr][1][1] = 0.;
        spjbufd[addr][1][2] = 0.;
        spjbufd[addr][1][3] = 0.;

        spjbufd[addr][2][0] = 0.;
        spjbufd[addr][2][1] = 0.;
        spjbufd[addr][2][2] = 0.;
        spjbufd[addr][2][3] = 0.;
    }

    void set_xi_one(const int addr, const double x, const double y, const double z
#ifdef USE_INDIVIDUAL_CUTOFF
                    , const double r_out, const double r_search
#endif
                    ){
        const int ah = addr / 8;
        const int al = addr % 8;
        xibuf[ah][0][al] = x;
        xibuf[ah][1][al] = y;
        xibuf[ah][2][al] = z;
#ifdef USE_INDIVIDUAL_CUTOFF
        xibuf[ah][3][al] = r_out;
        xibuf[ah][4][al] = r_search;
#endif
    }
    void set_xi_one_d(const int addr, const double x, const double y, const double z
#ifdef USE_INDIVIDUAL_CUTOFF
                      , const double r_out, const double r_search
#endif
                      ){
        const int ah = addr / 8;
        const int al = addr % 8;
        xibufd[ah][0][al] = x;
        xibufd[ah][1][al] = y;
        xibufd[ah][2][al] = z;
#ifdef USE_INDIVIDUAL_CUTOFF
        xibufd[ah][3][al] = r_out;
        xibufd[ah][4][al] = r_search;
#endif
    }

    template <typename real_t>
    void get_accp_one(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot, real_t &nngb, real_t &idxngb){
        const int ah = addr / 8;
        const int al = addr % 8;
        ax     = accpbuf[ah][0][al];
        ay     = accpbuf[ah][1][al];
        az     = accpbuf[ah][2][al];
        pot    = accpbuf[ah][3][al];
        nngb   = accpbuf[ah][4][al];
        idxngb = accpbuf[ah][5][al];
    }
    template <typename real_t>
    void get_accp_one_d(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot, real_t &nngb, real_t &idxngb){
        const int ah = addr / 8;
        const int al = addr % 8;
        ax     = accpbufd[ah][0][al];
        ay     = accpbufd[ah][1][al];
        az     = accpbufd[ah][2][al];
        pot    = accpbufd[ah][3][al];
        nngb   = accpbufd[ah][4][al];
        idxngb = accpbufd[ah][5][al];
    }

    template <typename real_t>
    void accum_accp_one(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot, real_t &nngb, real_t &idxngb){
        const int ah = addr / 8;
        const int al = addr % 8;
        ax    += accpbuf[ah][0][al];
        ay    += accpbuf[ah][1][al];
        az    += accpbuf[ah][2][al];
        pot   += accpbuf[ah][3][al];
        nngb  += accpbuf[ah][4][al];
        idxngb = accpbuf[ah][5][al];
    }
    template <typename real_t>
    void accum_accp_one_d(const int addr, real_t &ax, real_t &ay, real_t &az, real_t &pot, real_t &nngb, real_t &idxngb){
        const int ah = addr / 8;
        const int al = addr % 8;
        ax    += accpbufd[ah][0][al];
        ay    += accpbufd[ah][1][al];
        az    += accpbufd[ah][2][al];
        pot   += accpbufd[ah][3][al];
        nngb  += accpbufd[ah][4][al];
        idxngb = accpbufd[ah][5][al];
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

#ifdef __AVX512DQ__
        if ( nj%2==1 ) copyz_epj_one(nj, nj-1);
#endif

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

#ifdef __AVX512DQ__
        if ( nj%2==1 ) copyz_epj_one_d(nj, nj-1);
#endif
        
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

#ifdef __AVX512DQ__
        if ( nj%2==1 ) copyz_epj_one(nj, nj-1);
#endif
        
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

#ifdef __AVX512DQ__
        if ( nj%2==1 ) copyz_epj_one_d(nj, nj-1);
#endif
        
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

#ifdef __AVX512DQ__
        if ( nj%2==1 ) copyz_spj_one(nj, nj-1);
#endif

        kernel_spj_nounroll(ni, nj);
        // kernel_spj_unroll2(ni, nj);
    }

    void run_spj_d(const int ni, const int nj){
        if(ni > NIMAX || nj > NJMAX){
            std::cout<<"ni= "<<ni<<" NIMAX= "<<NIMAX<<" nj= "<<nj<<" NJMAX= "<<NJMAX<<std::endl;
        }
        assert(ni <= NIMAX);
        assert(nj <= NJMAX);

#ifdef __AVX512DQ__
        if ( nj%2==1 ) copyz_spj_one_d(nj, nj-1);
#endif

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
        const v8sf vzero = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        const v8sf vone  = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
#ifdef RSQRT_NR_EPJ_X2
        const v8sf v3p0   = {  3.0f,   3.0f,   3.0f,   3.0f,   3.0f,   3.0f,   3.0f,   3.0f};
        const v8sf v0p5   = {  0.5f,   0.5f,   0.5f,   0.5f,   0.5f,   0.5f,   0.5f,   0.5f}; 
        const v8sf v0p125 = {0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f}; 
#endif
        
		for(int i=0; i<ni; i+=8){
            const v8sf xi = *(v8sf *)(xibuf[i/8][0]);
			const v8sf yi = *(v8sf *)(xibuf[i/8][1]);
			const v8sf zi = *(v8sf *)(xibuf[i/8][2]);

			v8sf ax, ay, az, pot;
			ax = ay = az = pot = vzero;
            
            v8sf jbuf = _mm256_broadcast_ps((v4sf *)&epjbuf[0][0]);
            
            v8sf xj =  _mm256_shuffle_ps(jbuf, jbuf, 0x00);
            v8sf yj =  _mm256_shuffle_ps(jbuf, jbuf, 0x55);
            v8sf zj =  _mm256_shuffle_ps(jbuf, jbuf, 0xaa);
            v8sf mj =  _mm256_shuffle_ps(jbuf, jbuf, 0xff);
			for(int j=0; j<nj; j++){
                jbuf = _mm256_broadcast_ps((v4sf *)&epjbuf[j+1][0]);

				v8sf dx = xj - xi;
				v8sf dy = yj - yi;
				v8sf dz = zj - zi;

				v8sf r2   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v8sf ri1  = _mm256_rsqrt_ps(r2);
#ifdef RSQRT_NR_EPJ_X2
				ri1 *= (v3p0 - r2*(ri1*ri1));
#endif
				v8sf mri1 = mj * ri1;
				v8sf ri2  = ri1 * ri1;
				v8sf mri3 = mri1 * ri2;

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
			pot *= v0p5;
			ax  *= v0p125;
			ay  *= v0p125;
			az  *= v0p125;
#endif

			*(v8sf *)(accpbuf[i/8][0]) = ax;
			*(v8sf *)(accpbuf[i/8][1]) = ay;
			*(v8sf *)(accpbuf[i/8][2]) = az;
			*(v8sf *)(accpbuf[i/8][3]) = pot;
            *(v8sf *)(accpbuf[i/8][4]) = vzero;
            *(v8sf *)(accpbuf[i/8][5]) = -vone;
		}
	}

    __attribute__ ((noinline))
    void kernel_epj_nounroll_64bit(const int ni, const int nj){

        const v4df veps2 = {eps2, eps2, eps2, eps2};
        const v4df vzero = {0.0, 0.0, 0.0, 0.0};
        const v4df vone  = {1.0, 1.0, 1.0, 1.0};
#ifdef RSQRT_NR_EPJ_X2
        const v4df v0p5 = {0.5, 0.5, 0.5, 0.5};
        const v4df v3p0 = {3.0, 3.0, 3.0, 3.0};
#elif defined(RSQRT_NR_EPJ_X4)
        const v4df v8p0 = {8.0, 8.0, 8.0, 8.0};
        const v4df v6p0 = {6.0, 6.0, 6.0, 6.0};
        const v4df v5p0 = {5.0, 5.0, 5.0, 5.0};
        const v4df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
#endif
        
        for(int i=0; i<ni; i+=4){
            const int il = i%8;
            const v4df xi = *(v4df *)(&xibufd[i/8][0][il]);
            const v4df yi = *(v4df *)(&xibufd[i/8][1][il]);
            const v4df zi = *(v4df *)(&xibufd[i/8][2][il]);

            v4df ax, ay, az, pot;
            ax = ay = az = pot = vzero;

            v4df jbuf = *((v4df *)&epjbufd[0][0]);

            v4df xj =  _mm256_permute4x64_pd(jbuf, 0x00);
            v4df yj =  _mm256_permute4x64_pd(jbuf, 0x55);
            v4df zj =  _mm256_permute4x64_pd(jbuf, 0xaa);
            v4df mj =  _mm256_permute4x64_pd(jbuf, 0xff);

            for(int j=0; j<nj; j++){
                jbuf = *((v4df *)&epjbufd[j+1][0]);
                
                v4df dx = xj - xi;
                v4df dy = yj - yi;
                v4df dz = zj - zi;
                v4df r2   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v4df ri1  = _mm256_cvtps_pd( _mm_rsqrt_ps( _mm256_cvtpd_ps(r2)));
#ifdef RSQRT_NR_EPJ_X2
                //x2
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
                // x4
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
            *(v4df *)(&accpbufd[i/8][4][il]) = vzero;
            *(v4df *)(&accpbufd[i/8][5][il]) = -vone;
        }
    }

    ///////////////////////
    // with linear cutoff
    __attribute__ ((noinline))
    void kernel_epj_nounroll_for_p3t_with_linear_cutoff(const int ni, const int nj){
        
        const v8sf veps2 = {(float)eps2, (float)eps2, (float)eps2, (float)eps2,
                            (float)eps2, (float)eps2, (float)eps2, (float)eps2};
#ifndef USE_INDIVIDUAL_CUTOFF
        const v8sf vr_out     = {(float)r_out,    (float)r_out,    (float)r_out,    (float)r_out,
                                 (float)r_out,    (float)r_out,    (float)r_out,    (float)r_out};
        const v8sf vr_search  = {(float)r_search, (float)r_search, (float)r_search, (float)r_search,
                                 (float)r_search, (float)r_search, (float)r_search, (float)r_search};
        const v8sf vr_out2    = vr_out * vr_out;
        const v8sf vr_search2 = vr_search * vr_search;
#endif
        const v8sf vzero = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        const v8sf vone  = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
#ifdef RSQRT_NR_EPJ_X2
        const v8sf v3p0 = {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f}; 
        const v8sf v0p5 = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};
#elif defined(RSQRT_NR_EPJ_X4)
        const v8sf v8p0 = {8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f};
        const v8sf v6p0 = {6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f};
        const v8sf v5p0 = {5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f};
        const v8sf v0p0625 = {(float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0};
#endif
        //const v8sf allbits = _mm256_cmp_ps(vzero, vzero, _CMP_EQ_OS);
        
        for(int i=0; i<ni; i+=8){
            const v8sf xi = *(v8sf *)(xibuf[i/8][0]);
            const v8sf yi = *(v8sf *)(xibuf[i/8][1]);
            const v8sf zi = *(v8sf *)(xibuf[i/8][2]);
#ifdef USE_INDIVIDUAL_CUTOFF
            const v8sf r_outi    = *(v8sf *)(xibuf[i/8][3]);
            const v8sf r_searchi = *(v8sf *)(xibuf[i/8][4]);
#endif
            
            v8sf ax, ay, az, pot, nngb, idxngb;
            ax = ay = az = pot = nngb = vzero;
            idxngb = -vone;
            
            v8sf jbuf0 = _mm256_broadcast_ps((v4sf *)&epjbuf[0][0]);
            v8sf jbuf1 = _mm256_broadcast_ps((v4sf *)&epjbuf[0][1]);
            
            v8sf xj        =  _mm256_shuffle_ps(jbuf0, jbuf0, 0x00);
            v8sf yj        =  _mm256_shuffle_ps(jbuf0, jbuf0, 0x55);
            v8sf zj        =  _mm256_shuffle_ps(jbuf0, jbuf0, 0xaa);
            v8sf mj        =  _mm256_shuffle_ps(jbuf0, jbuf0, 0xff);
#ifdef USE_INDIVIDUAL_CUTOFF
            v8sf r_outj    =  _mm256_shuffle_ps(jbuf1, jbuf1, 0x00);
            v8sf r_searchj =  _mm256_shuffle_ps(jbuf1, jbuf1, 0x55);
#endif
            v8sf idxj      =  _mm256_shuffle_ps(jbuf1, jbuf1, 0xaa);
            
            for(int j=0; j<nj; j++){
                jbuf0 = _mm256_broadcast_ps((v4sf *)&epjbuf[j+1][0]);
                jbuf1 = _mm256_broadcast_ps((v4sf *)&epjbuf[j+1][1]);

#ifdef USE_INDIVIDUAL_CUTOFF
                v8sf vr_out     = _mm256_max_ps(r_outi,    r_outj);
                v8sf vr_search  = _mm256_max_ps(r_searchi, r_searchj);
                v8sf vr_out2    = vr_out    * vr_out;
                v8sf vr_search2 = vr_search * vr_search;
#endif

                v8sf dx = xj - xi;
                v8sf dy = yj - yi;
                v8sf dz = zj - zi;
                v8sf r2_real = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v8sf r2      = _mm256_max_ps(r2_real, vr_out2);
                v8sf ri1     = _mm256_rsqrt_ps(r2);
#ifdef RSQRT_NR_EPJ_X2
                // x2
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
                // x4
                v8sf h = vone - r2*(ri1*ri1);
                ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
                v8sf ri2  = ri1*ri1;
                v8sf mri1 = mj*ri1;
                v8sf mri3 = mri1 * ri2;

                v8sf mask0 = _mm256_cmp_ps(r2_real, vr_search2, _CMP_LT_OS); // for neighbour search
                v8sf mask1 = _mm256_and_ps(_mm256_cmp_ps(veps2, r2_real, _CMP_LT_OS), mask0);

                pot -= mri1;
                ax += mri3 * dx;
                ay += mri3 * dy;
                az += mri3 * dz;
                nngb  += _mm256_and_ps(mask0, vone);
                idxngb = _mm256_and_ps(mask1, idxj) + _mm256_andnot_ps(mask1, idxngb);

                xj        =  _mm256_shuffle_ps(jbuf0, jbuf0, 0x00);
                yj        =  _mm256_shuffle_ps(jbuf0, jbuf0, 0x55);
                zj        =  _mm256_shuffle_ps(jbuf0, jbuf0, 0xaa);
                mj        =  _mm256_shuffle_ps(jbuf0, jbuf0, 0xff);
#ifdef USE_INDIVIDUAL_CUTOFF
                r_outj    =  _mm256_shuffle_ps(jbuf1, jbuf1, 0x00);
                r_searchj =  _mm256_shuffle_ps(jbuf1, jbuf1, 0x55);
#endif
                idxj      =  _mm256_shuffle_ps(jbuf1, jbuf1, 0xaa);
            }
            
            *(v8sf *)(accpbuf[i/8][0]) = ax;
            *(v8sf *)(accpbuf[i/8][1]) = ay;
            *(v8sf *)(accpbuf[i/8][2]) = az;
            *(v8sf *)(accpbuf[i/8][3]) = pot;
            *(v8sf *)(accpbuf[i/8][4]) = nngb;
            *(v8sf *)(accpbuf[i/8][5]) = idxngb;
        }
    }

    __attribute__ ((noinline))
    void kernel_epj_64bit_nounroll_for_p3t_with_linear_cutoff(const int ni, const int nj){
        const v4df veps2 = {eps2, eps2, eps2, eps2};
#ifndef USE_INDIVIDUAL_CUTOFF
        const v4df vr_out    = {r_out,  r_out,  r_out,  r_out};
        const v4df vr_search = {r_search, r_search, r_search, r_search};
        const v4df vr_out2    = vr_out * vr_out;
        const v4df vr_search2 = vr_search * vr_search;
#endif
        const v4df vzero = {0.0, 0.0, 0.0, 0.0};
        const v4df vone  = {1.0, 1.0, 1.0, 1.0};
#ifdef RSQRT_NR_EPJ_X2
        const v4df v0p5 = {0.5, 0.5, 0.5, 0.5};
        const v4df v3p0 = {3.0, 3.0, 3.0, 3.0};
#elif defined(RSQRT_NR_EPJ_X4)
        const v4df v8p0 = {8.0, 8.0, 8.0, 8.0};
        const v4df v6p0 = {6.0, 6.0, 6.0, 6.0};
        const v4df v5p0 = {5.0, 5.0, 5.0, 5.0};
        const v4df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
#endif
        //const v8sf allbits = _mm256_cmp_pd(vzero, vzero, _CMP_EQ_OS);
        
        for(int i=0; i<ni; i+=4){
            const int il = i%8;
            const v4df xi = *(v4df *)(&xibufd[i/8][0][il]);
            const v4df yi = *(v4df *)(&xibufd[i/8][1][il]);
            const v4df zi = *(v4df *)(&xibufd[i/8][2][il]);
#ifdef USE_INDIVIDUAL_CUTOFF
            const v4df r_outi    = *(v4df *)(&xibufd[i/8][3][il]);
            const v4df r_searchi = *(v4df *)(&xibufd[i/8][4][il]);
#endif
            
            v4df ax, ay, az, pot, nngb, idxngb;
            ax = ay = az = pot = nngb = vzero;
            idxngb = -vone;
            
            v4df jbuf0 = *((v4df *)&epjbufd[0][0]);
            v4df jbuf1 = *((v4df *)&epjbufd[0][1]);
            v4df xj        =  _mm256_permute4x64_pd(jbuf0, 0x00);
            v4df yj        =  _mm256_permute4x64_pd(jbuf0, 0x55);
            v4df zj        =  _mm256_permute4x64_pd(jbuf0, 0xaa);
            v4df mj        =  _mm256_permute4x64_pd(jbuf0, 0xff);
#ifdef USE_INDIVIDUAL_CUTOFF
            v4df r_outj    =  _mm256_permute4x64_pd(jbuf1, 0x00);
            v4df r_searchj =  _mm256_permute4x64_pd(jbuf1, 0x55);
#endif
            v4df idxj      =  _mm256_permute4x64_pd(jbuf1, 0xaa);
            
            for(int j=0; j<nj; j++){
                jbuf0 = *((v4df *)&epjbufd[j+1][0]);
                jbuf1 = *((v4df *)&epjbufd[j+1][1]);

#ifdef USE_INDIVIDUAL_CUTOFF
                v4df vr_out     = _mm256_max_pd(r_outi,    r_outj);
                v4df vr_search  = _mm256_max_pd(r_searchi, r_searchj);
                v4df vr_out2    = vr_out    * vr_out;
                v4df vr_search2 = vr_search * vr_search;
#endif
                
                v4df dx = xj - xi;
                v4df dy = yj - yi;
                v4df dz = zj - zi;
                v4df r2_real = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v4df r2      = _mm256_max_pd(r2_real, vr_out2);
                v4df ri1     = _mm256_cvtps_pd( _mm_rsqrt_ps( _mm256_cvtpd_ps(r2)));
#ifdef RSQRT_NR_EPJ_X2
                //x2
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
                // x4
                v4df h = vone - r2*(ri1*ri1);
                ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
                
                v4df ri2 = ri1*ri1;
                v4df mri1 = mj*ri1;
                v4df mri3 = mri1 * ri2;

                v4df mask0 = _mm256_cmp_pd(r2_real, vr_search2, _CMP_LT_OS); // for neighbour search
                v4df mask1 = _mm256_and_pd(_mm256_cmp_pd(veps2, r2_real, _CMP_LT_OS), mask0);

                pot -= mri1;
                ax += mri3 * dx;
                ay += mri3 * dy;
                az += mri3 * dz;
                nngb  += _mm256_and_pd(mask0, vone);
                idxngb = _mm256_and_pd(mask1, idxj) + _mm256_andnot_pd(mask1, idxngb);

                xj        =  _mm256_permute4x64_pd(jbuf0, 0x00);
                yj        =  _mm256_permute4x64_pd(jbuf0, 0x55);
                zj        =  _mm256_permute4x64_pd(jbuf0, 0xaa);
                mj        =  _mm256_permute4x64_pd(jbuf0, 0xff);
#ifdef USE_INDIVIDUAL_CUTOFF
                r_outj    =  _mm256_permute4x64_pd(jbuf1, 0x00);
                r_searchj =  _mm256_permute4x64_pd(jbuf1, 0x55);
#endif
                idxj      =  _mm256_permute4x64_pd(jbuf1, 0xaa);
            }
            
            *(v4df *)(&accpbufd[i/8][0][il]) = ax;
            *(v4df *)(&accpbufd[i/8][1][il]) = ay;
            *(v4df *)(&accpbufd[i/8][2][il]) = az;
            *(v4df *)(&accpbufd[i/8][3][il]) = pot;
            *(v4df *)(&accpbufd[i/8][4][il]) = nngb;
            *(v4df *)(&accpbufd[i/8][5][il]) = idxngb;
        }
    }
    // linear cutoff
    //////////////
    

	__attribute__ ((noinline))
	void kernel_spj_nounroll(const int ni, const int nj){
        
	    const v8sf veps2 = {(float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2};
        const v8sf vzero = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        const v8sf vone  = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};  
        const v8sf v0p5  = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};
        const v8sf v2p5  = {2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f}; 
#ifdef RSQRT_NR_SPJ_X2
        const v8sf v3p0  = {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f};
#endif
        
		for(int i=0; i<ni; i+=8){
			const v8sf xi = *(v8sf *)(xibuf[i/8][0]);
			const v8sf yi = *(v8sf *)(xibuf[i/8][1]);
			const v8sf zi = *(v8sf *)(xibuf[i/8][2]);

			v8sf ax, ay, az, pot;
			ax = ay = az = pot = vzero;

#define PRELOAD_SPJ

#ifdef PRELOAD_SPJ
            v8sf jbuf0 = _mm256_broadcast_ps((v4sf *)&spjbuf[0][0]);
            v8sf jbuf1 = _mm256_broadcast_ps((v4sf *)&spjbuf[0][1]);
            v8sf jbuf2 = _mm256_broadcast_ps((v4sf *)&spjbuf[0][2]);
#else
			v8sf jbuf0, jbuf1, jbuf2;
#endif
			for(int j=0; j<nj; j++){
#ifndef PRELOAD_SPJ
                jbuf0 = _mm256_broadcast_ps((v4sf *)&spjbuf[j+0][0]);
#endif
                v8sf xj =  _mm256_shuffle_ps(jbuf0, jbuf0, 0x00);
                v8sf yj =  _mm256_shuffle_ps(jbuf0, jbuf0, 0x55);
                v8sf zj =  _mm256_shuffle_ps(jbuf0, jbuf0, 0xaa);
#ifdef PRELOAD_SPJ
                jbuf0 = _mm256_broadcast_ps((v4sf *)&spjbuf[j+1][0]);
#endif

#ifndef PRELOAD_SPJ
                jbuf1 = _mm256_broadcast_ps((v4sf *)&spjbuf[j+0][1]);
#endif
                v8sf qxx =  _mm256_shuffle_ps(jbuf1, jbuf1, 0x00);
                v8sf qyy =  _mm256_shuffle_ps(jbuf1, jbuf1, 0x55);
                v8sf qzz =  _mm256_shuffle_ps(jbuf1, jbuf1, 0xaa);
                v8sf mj  =  _mm256_shuffle_ps(jbuf1, jbuf1, 0xff);
#ifdef PRELOAD_SPJ
                jbuf1 = _mm256_broadcast_ps((v4sf *)&spjbuf[j+1][1]);
#endif

#ifndef PRELOAD_SPJ
                jbuf2 = _mm256_broadcast_ps((v4sf *)&spjbuf[j+0][2]);
#endif
                v8sf qxy =  _mm256_shuffle_ps(jbuf2, jbuf2, 0x00);
                v8sf qyz =  _mm256_shuffle_ps(jbuf2, jbuf2, 0x55);
                v8sf qzx =  _mm256_shuffle_ps(jbuf2, jbuf2, 0xaa);
                v8sf mtr =  _mm256_shuffle_ps(jbuf2, jbuf2, 0xff);
#ifdef PRELOAD_SPJ
                jbuf2 = _mm256_broadcast_ps((v4sf *)&spjbuf[j+1][2]);
#endif

				v8sf dx = xj - xi;
				v8sf dy = yj - yi;
				v8sf dz = zj - zi;

				v8sf r2  = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v8sf ri1  = _mm256_rsqrt_ps(r2);
#ifdef RSQRT_NR_SPJ_X2
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
            *(v8sf *)(accpbuf[i/8][4]) = vzero;
            *(v8sf *)(accpbuf[i/8][5]) = -vone;
		}
	}

    __attribute__ ((noinline))
    void kernel_spj_64bit_nounroll(const int ni, const int nj){
        
        const v4df veps2 = {eps2, eps2, eps2, eps2};
        const v4df vzero = {0.0, 0.0, 0.0, 0.0};
        const v4df vone  = {1.0, 1.0, 1.0, 1.0};
        const v4df v0p5  = {0.5, 0.5, 0.5, 0.5};
        const v4df v2p5  = {2.5, 2.5, 2.5, 2.5};
#ifdef RSQRT_NR_SPJ_X2
        const v4df v3p0  = {3.0, 3.0, 3.0, 3.0};
        //const v4df v0p125 = {0.125, 0.125, 0.125, 0.125};
#elif defined(RSQRT_NR_SPJ_X4)
        const v4df v8p0  = {8.0, 8.0, 8.0, 8.0};
        const v4df v6p0  = {6.0, 6.0, 6.0, 6.0};
        const v4df v5p0  = {5.0, 5.0, 5.0, 5.0};
        const v4df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
#endif
        
        for(int i=0; i<ni; i+=4){
            const int il = i%8;
            const v4df xi = *(v4df *)(&xibufd[i/8][0][il]);
            const v4df yi = *(v4df *)(&xibufd[i/8][1][il]);
            const v4df zi = *(v4df *)(&xibufd[i/8][2][il]);

            v4df ax, ay, az, pot;
            ax = ay = az = pot = vzero;

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
                v4df ri1  = _mm256_cvtps_pd( _mm_rsqrt_ps( _mm256_cvtpd_ps(r2)));
#ifdef RSQRT_NR_SPJ_X2
                //x2
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_SPJ_X4)
                // x4
                v4df h = vone - r2*(ri1*ri1);
                ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
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
		
                v4df meff  =  mj + v0p5 * rqr_ri4;
                v4df meff3 = (mj + v2p5 * rqr_ri4) * ri3;

                pot -= meff * ri1;

                ax = (ax - ri5*qr_x) + meff3*dx;
                ay = (ay - ri5*qr_y) + meff3*dy;
                az = (az - ri5*qr_z) + meff3*dz;
            }
	    
            *(v4df *)(&accpbufd[i/8][0][il]) = ax;
            *(v4df *)(&accpbufd[i/8][1][il]) = ay;
            *(v4df *)(&accpbufd[i/8][2][il]) = az;
            *(v4df *)(&accpbufd[i/8][3][il]) = pot;
            *(v4df *)(&accpbufd[i/8][4][il]) = vzero;
            *(v4df *)(&accpbufd[i/8][5][il]) = -vone;
        }
    }
    

#elif defined(__AVX512DQ__)
    //////////////////////
    ///   For AVX512   ///
    //////////////////////

    typedef float v16sf __attribute__((vector_size(64)));
    typedef double v8df __attribute__((vector_size(64)));

    typedef unsigned char  mask8;
    typedef unsigned short mask16;
    
    __attribute__ ((noinline))
	void kernel_epj_nounroll(const int ni, const int nj){
        
	    const v16sf veps2 = {(float)eps2, (float)eps2, (float)eps2, (float)eps2,
                             (float)eps2, (float)eps2, (float)eps2, (float)eps2,
                             (float)eps2, (float)eps2, (float)eps2, (float)eps2,
                             (float)eps2, (float)eps2, (float)eps2, (float)eps2};
        const v8sf vzero8  = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        const v8sf vone8   = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
        const v8sf vzero   = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                              0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
#ifdef RSQRT_NR_EPJ_X2
        const v16sf v3p0 = {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f,
                            3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f};
        const v16sf v0p5 = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
                            0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f}; 
        const v16sf v0p125 = {0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f,
                              0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f}; 
#endif
        
		for(int i=0; i<ni; i+=8){
            const v16sf xi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][0]));
			const v16sf yi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][1]));
			const v16sf zi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][2]));
            
			v16sf ax, ay, az, pot;
			ax = ay = az = pot = vzero;

            v16sf jbuf = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&epjbuf[1][0]),
                                            _mm256_broadcast_ps((v4sf *)&epjbuf[0][0]), 0x00);
            
            v16sf xj =  _mm512_shuffle_ps(jbuf, jbuf, 0x00);
            v16sf yj =  _mm512_shuffle_ps(jbuf, jbuf, 0x55);
            v16sf zj =  _mm512_shuffle_ps(jbuf, jbuf, 0xaa);
            v16sf mj =  _mm512_shuffle_ps(jbuf, jbuf, 0xff);
            if ( nj<=1 ) mj = _mm512_insertf32x8(mj, vzero8, 0xff);
			for(int j=0; j<nj; j+=2){
                jbuf = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&epjbuf[j+3][0]),
                                          _mm256_broadcast_ps((v4sf *)&epjbuf[j+2][0]), 0x00);
                
				v16sf dx = xj - xi;
				v16sf dy = yj - yi;
				v16sf dz = zj - zi;

				v16sf r2   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v16sf ri1  = _mm512_rsqrt14_ps(r2);
#ifdef RSQRT_NR_EPJ_X2
				ri1 *= (v3p0 - r2*(ri1*ri1));
#endif
				v16sf mri1 = mj * ri1;
				v16sf ri2  = ri1 * ri1;
				v16sf mri3 = mri1 * ri2;

                xj =  _mm512_shuffle_ps(jbuf, jbuf, 0x00);
                yj =  _mm512_shuffle_ps(jbuf, jbuf, 0x55);
                zj =  _mm512_shuffle_ps(jbuf, jbuf, 0xaa);
                mj =  _mm512_shuffle_ps(jbuf, jbuf, 0xff);
                if ( nj<=j+3 ) mj = _mm512_insertf32x8(mj, vzero8, 0xff);

				pot -= mri1;
				ax += mri3 * dx;
				ay += mri3 * dy;
				az += mri3 * dz;
			}
#ifdef RSQRT_NR_EPJ_X2
			pot *= v0p5;
			ax  *= v0p125;
			ay  *= v0p125;
			az  *= v0p125;
#endif

			*(v8sf *)(accpbuf[i/8][0]) = _mm512_extractf32x8_ps(ax,  0x00) + _mm512_extractf32x8_ps(ax,  0xff);
			*(v8sf *)(accpbuf[i/8][1]) = _mm512_extractf32x8_ps(ay,  0x00) + _mm512_extractf32x8_ps(ay,  0xff);
			*(v8sf *)(accpbuf[i/8][2]) = _mm512_extractf32x8_ps(az,  0x00) + _mm512_extractf32x8_ps(az,  0xff);
			*(v8sf *)(accpbuf[i/8][3]) = _mm512_extractf32x8_ps(pot, 0x00) + _mm512_extractf32x8_ps(pot, 0xff);
            *(v8sf *)(accpbuf[i/8][4]) = vzero8;
            *(v8sf *)(accpbuf[i/8][5]) = -vone8;
		}
	}

    __attribute__ ((noinline))
    void kernel_epj_nounroll_64bit(const int ni, const int nj){

        const v8df veps2 = {eps2, eps2, eps2, eps2, eps2, eps2, eps2, eps2};
        const v4df vzero4 = {0.0, 0.0, 0.0, 0.0};
        const v4df vone4  = {1.0, 1.0, 1.0, 1.0};
        const v8df vzero = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#ifdef RSQRT_NR_EPJ_X2
        const v8df v0p5 = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
        const v8df v3p0 = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
#elif defined(RSQRT_NR_EPJ_X4)
        const v8df vone = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        const v8df v8p0 = {8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0};
        const v8df v6p0 = {6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0};
        const v8df v5p0 = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
        const v8df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0,
                              1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
#endif
        
        for(int i=0; i<ni; i+=4){
            const int il = i%8;
            const v8df xi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][0][il]));
            const v8df yi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][1][il]));
            const v8df zi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][2][il]));
            
            v8df ax, ay, az, pot;
            ax = ay = az = pot = vzero;
            
            v8df jbuf = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&epjbufd[1][0])),
                                           *((v4df *)&epjbufd[0][0]), 0x00);

            v8df xj =  _mm512_permutex_pd(jbuf, 0x00);
            v8df yj =  _mm512_permutex_pd(jbuf, 0x55);
            v8df zj =  _mm512_permutex_pd(jbuf, 0xaa);
            v8df mj =  _mm512_permutex_pd(jbuf, 0xff);
            if ( nj<=1 ) mj = _mm512_insertf64x4(mj, vzero4, 0xff);

            for(int j=0; j<nj; j+=2){
                jbuf = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&epjbufd[j+3][0])),
                                          *((v4df *)&epjbufd[j+2][0]), 0x00);
                
                v8df dx = xj - xi;
                v8df dy = yj - yi;
                v8df dz = zj - zi;
                v8df r2   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v8df ri1  = _mm512_cvtps_pd( _mm256_rsqrt_ps( _mm512_cvtpd_ps(r2)));
#ifdef RSQRT_NR_EPJ_X2
                //x2
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
                // x4
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
                if ( nj<=j+3 ) mj = _mm512_insertf64x4(mj, vzero4, 0xff);

                pot -= mri1;
                ax += mri3 * dx;
                ay += mri3 * dy;
                az += mri3 * dz;
            }
            *(v4df *)(&accpbufd[i/8][0][il]) = _mm512_extractf64x4_pd(ax,  0x00) + _mm512_extractf64x4_pd(ax,  0xff);
			*(v4df *)(&accpbufd[i/8][1][il]) = _mm512_extractf64x4_pd(ay,  0x00) + _mm512_extractf64x4_pd(ay,  0xff);
			*(v4df *)(&accpbufd[i/8][2][il]) = _mm512_extractf64x4_pd(az,  0x00) + _mm512_extractf64x4_pd(az,  0xff);
			*(v4df *)(&accpbufd[i/8][3][il]) = _mm512_extractf64x4_pd(pot, 0x00) + _mm512_extractf64x4_pd(pot, 0xff);
            *(v4df *)(&accpbufd[i/8][4][il]) = vzero4;
            *(v4df *)(&accpbufd[i/8][5][il]) = -vone4;
        }
    }

    ///////////////////////
    // with linear cutoff
    __attribute__ ((noinline))
    void kernel_epj_nounroll_for_p3t_with_linear_cutoff(const int ni, const int nj){
        
        const v16sf veps2 = {(float)eps2, (float)eps2, (float)eps2, (float)eps2,
                             (float)eps2, (float)eps2, (float)eps2, (float)eps2,
                             (float)eps2, (float)eps2, (float)eps2, (float)eps2,
                             (float)eps2, (float)eps2, (float)eps2, (float)eps2};
        const v8sf  vzero8 = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        const v8sf  vone8  = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
        const v16sf vone   = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                              1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
        const v16sf vzero  = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                              0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
#ifdef RSQRT_NR_EPJ_X2
        const v16sf v3p0 = {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f,
                            3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f}; 
        const v16sf v0p5 = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
                            0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};
#elif defined(RSQRT_NR_EPJ_X4)
        const v16sf v8p0 = {8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f,
                            8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f};
        const v16sf v6p0 = {6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f,
                            6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f, 6.0f};
        const v16sf v5p0 = {5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f,
                            5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f};
        const v16sf v0p0625 = {(float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0,
                               (float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0, (float)1.0/16.0};
#endif
#ifndef USE_INDIVIDUAL_CUTOFF
        const v16sf vr_out    = {(float)r_out,    (float)r_out,    (float)r_out,    (float)r_out,
                                 (float)r_out,    (float)r_out,    (float)r_out,    (float)r_out,
                                 (float)r_out,    (float)r_out,    (float)r_out,    (float)r_out,
                                 (float)r_out,    (float)r_out,    (float)r_out,    (float)r_out};
        const v16sf vr_search = {(float)r_search, (float)r_search, (float)r_search, (float)r_search,
                                 (float)r_search, (float)r_search, (float)r_search, (float)r_search,
                                 (float)r_search, (float)r_search, (float)r_search, (float)r_search,
                                 (float)r_search, (float)r_search, (float)r_search, (float)r_search};
        const v16sf vr_out2    = vr_out    * vr_out;
        const v16sf vr_search2 = vr_search * vr_search;
#endif

        for(int i=0; i<ni; i+=8){
            const v16sf xi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][0]));
            const v16sf yi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][1]));
            const v16sf zi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][2]));
#ifdef USE_INDIVIDUAL_CUTOFF
            const v16sf r_outi    = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][3]));
            const v16sf r_searchi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][4]));
#endif
            
            v16sf ax, ay, az, pot, nngb, idxngb;
            ax = ay = az = pot = nngb = vzero;
            idxngb = -vone;
            
            v16sf jbuf0 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&epjbuf[1][0]),
                                             _mm256_broadcast_ps((v4sf *)&epjbuf[0][0]), 0x00);
            v16sf jbuf1 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&epjbuf[1][1]),
                                             _mm256_broadcast_ps((v4sf *)&epjbuf[0][1]), 0x00);
            v16sf xj        =  _mm512_shuffle_ps(jbuf0, jbuf0, 0x00);
            v16sf yj        =  _mm512_shuffle_ps(jbuf0, jbuf0, 0x55);
            v16sf zj        =  _mm512_shuffle_ps(jbuf0, jbuf0, 0xaa);
            v16sf mj        =  _mm512_shuffle_ps(jbuf0, jbuf0, 0xff);
#ifdef USE_INDIVIDUAL_CUTOFF
            v16sf r_outj    =  _mm512_shuffle_ps(jbuf1, jbuf1, 0x00);
            v16sf r_searchj =  _mm512_shuffle_ps(jbuf1, jbuf1, 0x55);
#endif
            v16sf idxj      =  _mm512_shuffle_ps(jbuf1, jbuf1, 0xaa);
            //if ( nj<=1 ) mj = _mm512_insertf32x8(mj, vzero8, 0xff);
            
            for(int j=0; j<nj; j+=2){
                jbuf0 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&epjbuf[j+3][0]),
                                           _mm256_broadcast_ps((v4sf *)&epjbuf[j+2][0]), 0x00);
                jbuf1 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&epjbuf[j+3][1]),
                                           _mm256_broadcast_ps((v4sf *)&epjbuf[j+2][1]), 0x00);
#ifdef USE_INDIVIDUAL_CUTOFF
                v16sf vr_out     = _mm512_max_ps(r_outi,    r_outj);
                v16sf vr_search  = _mm512_max_ps(r_searchi, r_searchj);
                v16sf vr_out2    = vr_out    * vr_out;
                v16sf vr_search2 = vr_search * vr_search;
#endif
                
                v16sf dx = xj - xi;
                v16sf dy = yj - yi;
                v16sf dz = zj - zi;
                v16sf r2_real = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v16sf r2      = _mm512_max_ps(r2_real, vr_out2);
                v16sf ri1     = _mm512_rsqrt14_ps(r2);
#ifdef RSQRT_NR_EPJ_X2
                // x2
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
                // x4
                v16sf h = vone - r2*(ri1*ri1);
                ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
                v16sf ri2 = ri1*ri1;
                v16sf mri1 = mj*ri1;
                v16sf mri3 = mri1 * ri2;

                mask16 mask0 = _mm512_mask_cmp_ps_mask(_mm512_cmp_ps_mask(idxj, vzero, _CMP_GE_OQ),
                                                       r2_real, vr_search2, _CMP_LT_OS); // for neighbour search
                mask16 mask1 = _mm512_mask_cmp_ps_mask(mask0, veps2, r2_real, _CMP_LT_OS);

                pot -= mri1;
                ax += mri3 * dx;
                ay += mri3 * dy;
                az += mri3 * dz;
                nngb  += _mm512_maskz_shuffle_ps(mask0, vone, vone, 0x00);
                idxngb = _mm512_mask_shuffle_ps(idxngb, mask1, idxj, idxj, 0x00);

                xj        = _mm512_shuffle_ps(jbuf0, jbuf0, 0x00);
                yj        = _mm512_shuffle_ps(jbuf0, jbuf0, 0x55);
                zj        = _mm512_shuffle_ps(jbuf0, jbuf0, 0xaa);
                mj        = _mm512_shuffle_ps(jbuf0, jbuf0, 0xff);
#ifdef USE_INDIVIDUAL_CUTOFF
                r_outj    = _mm512_shuffle_ps(jbuf1, jbuf1, 0x00);
                r_searchj = _mm512_shuffle_ps(jbuf1, jbuf1, 0x55);
#endif
                idxj      = _mm512_shuffle_ps(jbuf1, jbuf1, 0xaa);
                //if ( nj<=j+3 ) mj = _mm512_insertf32x8(mj, vzero8, 0xff);
            }
            *(v8sf *)(accpbuf[i/8][0]) = _mm512_extractf32x8_ps(ax,   0x00) + _mm512_extractf32x8_ps(ax,   0xff);
			*(v8sf *)(accpbuf[i/8][1]) = _mm512_extractf32x8_ps(ay,   0x00) + _mm512_extractf32x8_ps(ay,   0xff);
			*(v8sf *)(accpbuf[i/8][2]) = _mm512_extractf32x8_ps(az,   0x00) + _mm512_extractf32x8_ps(az,   0xff);
			*(v8sf *)(accpbuf[i/8][3]) = _mm512_extractf32x8_ps(pot,  0x00) + _mm512_extractf32x8_ps(pot,  0xff);
            *(v8sf *)(accpbuf[i/8][4]) = _mm512_extractf32x8_ps(nngb, 0x00) + _mm512_extractf32x8_ps(nngb, 0xff);
            *(v8sf *)(accpbuf[i/8][5]) = _mm256_max_ps(_mm512_extractf32x8_ps(idxngb,0x00), _mm512_extractf32x8_ps(idxngb, 0xff));
        }
    }

    __attribute__ ((noinline))
    void kernel_epj_64bit_nounroll_for_p3t_with_linear_cutoff(const int ni, const int nj){
        
        const v8df veps2 = {eps2, eps2, eps2, eps2, eps2, eps2, eps2, eps2};
        const v8df vzero  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        const v8df vone   = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        const v4df vzero4 = {0.0, 0.0, 0.0, 0.0};
        const v4df vone4  = {1.0, 1.0, 1.0, 1.0};
#ifdef RSQRT_NR_EPJ_X2
        const v8df v0p5 = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
        const v8df v3p0 = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
#elif defined(RSQRT_NR_EPJ_X4)
        const v8df v8p0 = {8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0};
        const v8df v6p0 = {6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0};
        const v8df v5p0 = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
        const v8df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
#endif
#ifndef USE_INDIVIDUAL_CUTOFF
        const v8df vr_out     = {r_out, r_out, r_out, r_out, r_out, r_out, r_out, r_out};
        const v8df vr_search  = {r_search, r_search, r_search, r_search, r_search, r_search, r_search, r_search};
        const v8df vr_out2    = vr_out    * vr_out;
        const v8df vr_search2 = vr_search * vr_search;
#endif
        for(int i=0; i<ni; i+=4){
            const int il = i%8;
            const v8df xi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][0][il]));
            const v8df yi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][1][il]));
            const v8df zi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][2][il]));
#ifdef USE_INDIVIDUAL_CUTOFF
            const v8df r_outi    = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][3][il]));
            const v8df r_searchi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][4][il]));
#endif
            
            v8df ax, ay, az, pot, nngb, idxngb;
            ax = ay = az = pot = nngb = vzero;
            idxngb = -vone;
            
            v8df jbuf0 = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&epjbufd[1][0])),
                                            *((v4df *)&epjbufd[0][0]), 0x00);
            v8df jbuf1 = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&epjbufd[1][1])),
                                            *((v4df *)&epjbufd[0][1]), 0x00);
            v8df xj        =  _mm512_permutex_pd(jbuf0, 0x00);
            v8df yj        =  _mm512_permutex_pd(jbuf0, 0x55);
            v8df zj        =  _mm512_permutex_pd(jbuf0, 0xaa);
            v8df mj        =  _mm512_permutex_pd(jbuf0, 0xff);
#ifdef USE_INDIVIDUAL_CUTOFF
            v8df r_outj    =  _mm512_permutex_pd(jbuf1, 0x00);
            v8df r_searchj =  _mm512_permutex_pd(jbuf1, 0x55);
#endif
            v8df idxj      =  _mm512_permutex_pd(jbuf1, 0xaa);
            if ( nj<=1 ) mj = _mm512_insertf64x4(mj, vzero4, 0xff);
            
            for(int j=0; j<nj; j+=2){
                jbuf0 = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&epjbufd[j+3][0])),
                                           *((v4df *)&epjbufd[j+2][0]), 0x00);
                jbuf1 = _mm512_insertf64x4(_mm512_broadcast_f64x4(*((v4df *)&epjbufd[j+3][1])),
                                           *((v4df *)&epjbufd[j+2][1]), 0x00);
#ifdef USE_INDIVIDUAL_CUTOFF
                v8df vr_out     = _mm512_max_pd(r_outi,    r_outj);
                v8df vr_search  = _mm512_max_pd(r_searchi, r_searchj);
                v8df vr_out2    = vr_out    * vr_out;
                v8df vr_search2 = vr_search * vr_search;
#endif
                v8df dx = xj - xi;
                v8df dy = yj - yi;
                v8df dz = zj - zi;
                v8df r2_real   = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v8df r2  = _mm512_max_pd( r2_real, vr_out2);
                v8df ri1 = _mm512_cvtps_pd( _mm256_rsqrt_ps( _mm512_cvtpd_ps(r2)));
#ifdef RSQRT_NR_EPJ_X2
                //x2
                ri1 *= (v3p0 - r2*(ri1*ri1))*v0p5;
#elif defined(RSQRT_NR_EPJ_X4)
                // x4
                v8df h = vone - r2*(ri1*ri1);
                ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
#endif
                v8df ri2 = ri1*ri1;
                v8df mri1 = mj*ri1;
                v8df mri3 = mri1 * ri2;

                mask8 mask0 = _mm512_mask_cmp_pd_mask(_mm512_cmp_pd_mask(idxj, vzero, _CMP_GE_OQ),
                                                      r2_real, vr_search2, _CMP_LT_OS); // for neighbour search
                mask8 mask1 = _mm512_mask_cmp_pd_mask(mask0, veps2, r2_real, _CMP_LT_OS);

                pot -= mri1;
                ax += mri3 * dx;
                ay += mri3 * dy;
                az += mri3 * dz;
                nngb  += _mm512_maskz_shuffle_pd(mask0, vone, vone, 0x00);
                idxngb = _mm512_mask_shuffle_pd(idxngb, mask1, idxj, idxj, 0x00);
                
                xj        = _mm512_permutex_pd(jbuf0, 0x00);
                yj        = _mm512_permutex_pd(jbuf0, 0x55);
                zj        = _mm512_permutex_pd(jbuf0, 0xaa);
                mj        = _mm512_permutex_pd(jbuf0, 0xff);
#ifdef USE_INDIVIDUAL_CUTOFF
                r_outj    = _mm512_permutex_pd(jbuf1, 0x00);
                r_searchj = _mm512_permutex_pd(jbuf1, 0x55);
#endif
                idxj      = _mm512_permutex_pd(jbuf1, 0xaa);
                //if ( nj<=j+3 ) mj = _mm512_insertf64x4(mj, vzero, 0xff);
                                
            }
            *(v4df *)(&accpbufd[i/8][0][il]) = _mm512_extractf64x4_pd(ax,   0x00) + _mm512_extractf64x4_pd(ax,   0xff);
			*(v4df *)(&accpbufd[i/8][1][il]) = _mm512_extractf64x4_pd(ay,   0x00) + _mm512_extractf64x4_pd(ay,   0xff);
			*(v4df *)(&accpbufd[i/8][2][il]) = _mm512_extractf64x4_pd(az,   0x00) + _mm512_extractf64x4_pd(az,   0xff);
			*(v4df *)(&accpbufd[i/8][3][il]) = _mm512_extractf64x4_pd(pot,  0x00) + _mm512_extractf64x4_pd(pot,  0xff);
            *(v4df *)(&accpbufd[i/8][4][il]) = _mm512_extractf64x4_pd(nngb, 0x00) + _mm512_extractf64x4_pd(nngb, 0xff);
            *(v4df *)(&accpbufd[i/8][5][il]) = _mm256_max_pd(_mm512_extractf64x4_pd(idxngb,0x00), _mm512_extractf64x4_pd(idxngb,0xff));
        }
    }
    // linear cutoff
    //////////////

    __attribute__ ((noinline))
	void kernel_spj_nounroll(const int ni, const int nj){
        
	    const v16sf veps2  = {(float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2,
                             (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2, (float)eps2};
        const v16sf vzero  = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                              0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        const v16sf vone   = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
                              1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};  
        const v8sf  vzero8 = (v8sf){0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        
        const v16sf v0p5   = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
                              0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};
        
        const v16sf v2p5   = {2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f,
                              2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f, 2.5f}; 
#ifdef RSQRT_NR_SPJ_X2
        const v16sf v3p0   = {3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f,
                              3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f};
#endif
        
		for(int i=0; i<ni; i+=8){
            const v16sf xi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][0]));
			const v16sf yi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][1]));
			const v16sf zi = _mm512_broadcast_f32x8(*(v8sf *)(xibuf[i/8][2]));

			v16sf ax, ay, az, pot;
			ax = ay = az = pot = vzero;

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
                //if ( nj<=j+1 ) {
                //    qxx = _mm512_insertf32x8(qxx, vzero8, 0xff);
                //    qyy = _mm512_insertf32x8(qyy, vzero8, 0xff);
                //    qzz = _mm512_insertf32x8(qzz, vzero8, 0xff);
                //    mj  = _mm512_insertf32x8(mj,  vzero8, 0xff);
                //    qxy = _mm512_insertf32x8(qxy, vzero8, 0xff);
                //    qyz = _mm512_insertf32x8(qyz, vzero8, 0xff);
                //    qzx = _mm512_insertf32x8(qzx, vzero8, 0xff);
                //    mtr = _mm512_insertf32x8(mtr, vzero8, 0xff);
                //}
#ifdef PRELOAD_SPJ
                jbuf2 = _mm512_insertf32x8(_mm512_broadcast_f32x4(*(v4sf *)&spjbuf[j+3][2]),
                                           _mm256_broadcast_ps((v4sf *)&spjbuf[j+2][2]), 0x00);
#endif

				v16sf dx = xj - xi;
				v16sf dy = yj - yi;
				v16sf dz = zj - zi;

				v16sf r2  = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v16sf ri1  = _mm512_rsqrt14_ps(r2);
#ifdef RSQRT_NR_SPJ_X2
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

				v16sf meff  =  mj + v0p5 * rqr_ri4;
				v16sf meff3 = (mj + v2p5 * rqr_ri4) * ri3;

				pot -= meff * ri1;

				ax = (ax - ri5*qr_x) + meff3*dx;
				ay = (ay - ri5*qr_y) + meff3*dy;
				az = (az - ri5*qr_z) + meff3*dz;
			}
            *(v8sf *)(accpbuf[i/8][0]) = _mm512_extractf32x8_ps(ax,  0x00) + _mm512_extractf32x8_ps(ax,  0xff);
			*(v8sf *)(accpbuf[i/8][1]) = _mm512_extractf32x8_ps(ay,  0x00) + _mm512_extractf32x8_ps(ay,  0xff);
			*(v8sf *)(accpbuf[i/8][2]) = _mm512_extractf32x8_ps(az,  0x00) + _mm512_extractf32x8_ps(az,  0xff);
			*(v8sf *)(accpbuf[i/8][3]) = _mm512_extractf32x8_ps(pot, 0x00) + _mm512_extractf32x8_ps(pot, 0xff);
            *(v8sf *)(accpbuf[i/8][4]) = vzero;
            *(v8sf *)(accpbuf[i/8][5]) = -vone;
		}
	}

    __attribute__ ((noinline))
    void kernel_spj_64bit_nounroll(const int ni, const int nj){
        const v8df veps2  = {eps2, eps2, eps2, eps2, eps2, eps2, eps2, eps2};
        const v8df vzero  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        const v8df vone   = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        const v8df v0p5   = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
        const v8df v2p5   = {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5};
        const v4df vzero4 = (v4df){0.0, 0.0, 0.0, 0.0};
#ifdef RSQRT_NR_SPJ_X2
        const v8df v3p0 = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
#elif defined(RSQRT_NR_SPJ_X4)
        const v8df v8p0 = {8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0};
        const v8df v6p0 = {6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0};
        const v8df v5p0 = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
        const v8df v0p0625 = {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0};
#endif
        
        for(int i=0; i<ni; i+=4){
            const int il = i%8;
            const v8df xi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][0][il]));
            const v8df yi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][1][il]));
            const v8df zi = _mm512_broadcast_f64x4(*(v4df *)(&xibufd[i/8][2][il]));

            v8df ax, ay, az, pot;
            ax = ay = az = pot = vzero;

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
                //if ( nj<=j+1 ) {
                //    qxx = _mm512_insertf64x4(qxx, vzero4, 0xff);
                //    qyy = _mm512_insertf64x4(qyy, vzero4, 0xff);
                //    qzz = _mm512_insertf64x4(qzz, vzero4, 0xff);
                //    mj  = _mm512_insertf64x4(mj,  vzero4, 0xff);
                //    qxy = _mm512_insertf64x4(qxy, vzero4, 0xff);
                //    qyz = _mm512_insertf64x4(qyz, vzero4, 0xff);
                //    qzx = _mm512_insertf64x4(qzx, vzero4, 0xff);
                //    mtr = _mm512_insertf64x4(mtr, vzero4, 0xff);
                //}
		
                v8df dx = xj - xi;
                v8df dy = yj - yi;
                v8df dz = zj - zi;
		
                v8df r2  = ((veps2 + dx*dx) + dy*dy) + dz*dz;
                v8df ri1  = _mm512_cvtps_pd( _mm256_rsqrt_ps( _mm512_cvtpd_ps(r2)));
		
#ifdef RSQRT_NR_SPJ_X2
                //x2
                ri1 *= (v3p0 - r2*(ri1*ri1)) * v0p5;
#elif defined(RSQRT_NR_SPJ_X4)
                // x4
                v8df h = v1p0 - r2*(ri1*ri1);
                ri1 *= vone + h*(v8p0+h*(v6p0+v5p0*h))*v0p0625;
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
		
                v8df meff  =  mj + v0p5 * rqr_ri4;
                v8df meff3 = (mj + v2p5 * rqr_ri4) * ri3;

                pot -= meff * ri1;

                ax = (ax - ri5*qr_x) + meff3*dx;
                ay = (ay - ri5*qr_y) + meff3*dy;
                az = (az - ri5*qr_z) + meff3*dz;
            }

            *(v4df *)(&accpbufd[i/8][0][il]) = _mm512_extractf64x4_pd(ax,  0x00) + _mm512_extractf64x4_pd(ax,  0xff);
			*(v4df *)(&accpbufd[i/8][1][il]) = _mm512_extractf64x4_pd(ay,  0x00) + _mm512_extractf64x4_pd(ay,  0xff);
			*(v4df *)(&accpbufd[i/8][2][il]) = _mm512_extractf64x4_pd(az,  0x00) + _mm512_extractf64x4_pd(az,  0xff);
			*(v4df *)(&accpbufd[i/8][3][il]) = _mm512_extractf64x4_pd(pot, 0x00) + _mm512_extractf64x4_pd(pot, 0xff);
            *(v4df *)(&accpbufd[i/8][4][il]) = vzero;
            *(v4df *)(&accpbufd[i/8][5][il]) = -vone;
        }
    }

#endif //__AVX512DQ__


} __attribute__ ((aligned(128)));







