#pragma once
#include "../ps_defs.hpp"
#include "multidimensional_array.hpp"

namespace ParticleSimulator {
    namespace ParticleMeshMultipole {

        template <typename real_t=double,
                  typename cplx_t=std::complex<real_t> >
        class GreenFunction {
        public:
            S32 p, bc, NX, NY, NZ;
            MultidimensionalArray<real_t, 4> gf_r;
            MultidimensionalArray<cplx_t, 4> gf_k;
      
        private:
            inline void initOpenImpl() {
                const S32 LEN2 = (2*p+1)*(2*p+1);
                S32 sizes_r[4] = {2*NZ, 2*NY, 2*NX, LEN2};
                S32 sizes_k[4] = {2*NZ, 2*NY, 1+NX, LEN2};
                gf_r.initialize(sizes_r);
                gf_k.initialize(sizes_k);
            }
    
            inline void initPeriodicXYZImpl() {
                const S32 LEN2 = (2*p+1)*(2*p+1);
                S32 sizes_r[4] = {NZ, NY, NX,     LEN2};
                S32 sizes_k[4] = {NZ, NY, 1+NX/2, LEN2};
                gf_r.initialize(sizes_r);
                gf_k.initialize(sizes_k);
            }
    
            void setOpenImpl(const S32 icut, const F64vec cell_length) {
                const S32 n_thread = Comm::getNumberOfThread();
                Slm<real_t> * slm = new Slm<real_t>[n_thread];
                for (S32 i=0; i<n_thread; i++) slm[i].alloc(2*p);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                for(S32 k=0; k<2*NZ; k++) for(S32 j=0; j<2*NY; j++) for(S32 i=0; i<2*NX; i++)
                {
                    const S32 ith = Comm::getThreadNum();
                    const S32 LEN2 = (2*p+1)*(2*p+1);
                    if (k==NZ || j==NY || i==NX){
                        for (S32 lm=0; lm<LEN2; lm++){
                            gf_r(k,j,i,lm) = 0.0;
                        }
                        continue;
                    }
                    const S32 kk = (k>NZ) ? k - 2*NZ : k;
                    const S32 jj = (j>NY) ? j - 2*NY : j;
                    const S32 ii = (i>NX) ? i - 2*NX : i;
                    if(abs(kk)<=icut && abs(jj)<=icut && abs(ii)<=icut){
                        for(S32 lm=0; lm<LEN2; lm++){
                            gf_r(k,j,i,lm) = 0.0;
                        }
                        continue;
                    }
                    const F64 dx = cell_length.x * F64(ii);
                    const F64 dy = cell_length.y * F64(jj);
                    const F64 dz = cell_length.z * F64(kk);
                    slm[ith].eval_opt(-dx, dy, dz); // eval S_l^{-m}
                    for(int lm=0; lm<LEN2; lm++){
                        gf_r(k,j,i,lm) = slm[ith].buf[lm];
                    }
                }
                // Free memory
                for (S32 i=0; i<n_thread; i++) slm[i].freeMem();
                delete [] slm;
            }
    
            void setPeriodicXYZImpl(const S32 icut,
                                    const F64vec cell_length, 
                                    const F64 alpha,
                                    const S32 NMAX,
                                    const S32 MMAX) {
                const S32 n_thread = Comm::getNumberOfThread();   
                CutFunc<real_t, cplx_t> *cf;
                cf = new CutFunc<real_t, cplx_t>[n_thread];
                Slm<real_t> *slm, *rsum, *ksum;
                slm  = new Slm<real_t>[n_thread];
                rsum = new Slm<real_t>[n_thread];
                ksum = new Slm<real_t>[n_thread];
                for (S32 i=0; i < n_thread; i++) {
                    cf[i].init(2*p);
                    slm[i].alloc(2*p);
                    rsum[i].alloc(2*p);
                    ksum[i].alloc(2*p);
                }
                {
                    Slm<real_t> tmp;
                    tmp.dry_run(2*p);
                }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                for(S32 k=0; k<NZ; k++) for(S32 j=0; j<NY; j++) for(S32 i=0; i<NX; i++)
                {
                    const S32 ith = Comm::getThreadNum();
                    const S32 LEN2 = (2*p+1)*(2*p+1);

                    // real-space sum
                    rsum[ith].clear();
                    const S32 kk = (k>=NZ/2) ? k - NZ : k;
                    const S32 jj = (j>=NY/2) ? j - NY : j;
                    const S32 ii = (i>=NX/2) ? i - NX : i;
                    S32 nx, ny, nz;
                    for(nz=-NMAX; nz<=NMAX; nz++) for(ny=-NMAX; ny<=NMAX; ny++) for(nx=-NMAX; nx<=NMAX; nx++)
                    {
                        const S32 kkk = kk + nz*NZ;
                        const S32 jjj = jj + ny*NY;
                        const S32 iii = ii + nx*NX;
                        if( 0 == (iii|jjj|kkk) ) continue;
                        const F64 dx = cell_length.x * F64(iii);
                        const F64 dy = cell_length.y * F64(jjj);
                        const F64 dz = cell_length.z * F64(kkk);
                        const F64 dr2 = dx*dx + dy*dy + dz*dz;
                        cf[ith].eval_rcut(dr2 * (alpha*alpha));
                        slm[ith].eval_opt(-dx, dy, dz); // eval S_l^{-m}
                        // near cell correction
                        const bool near = (abs(kkk)<=icut && abs(jjj)<=icut && abs(iii)<=icut);
                        if(near){
                            for(S32 l=0; l<=2*p; l++){
                                cf[ith].rcut[l] -= 1.0;
                            }
                        }
                        for(S32 l=0; l<=2*p; l++){
                            for(S32 m=0; m<=l; m++){
                                const cplx_t val = cf[ith].rcut[l] * slm[ith].val_at(l, m);
                                rsum[ith].accum_at(l, m, val);
                            }
                        }
                    }
        
                    // wave-space sum
                    ksum[ith].clear();
                    S32 mx, my, mz;
                    for(mz=-MMAX; mz<=MMAX; mz++) for(my=-MMAX; my<=MMAX; my++) for(mx=-MMAX; mx<=MMAX; mx++)
                    {
                        if(0 == (mx|my|mz)) continue;
#if 1
                        const F64 dx = cell_length.x * F64(ii);
                        const F64 dy = cell_length.y * F64(jj);
                        const F64 dz = cell_length.z * F64(kk);
#else
                        // This code is the same as that of the original PMMM code.
                        // But, the use of (i,j,k) is probably WRONG.
                        const F64 dx = cell_length.x * F64(i);
                        const F64 dy = cell_length.y * F64(j);
                        const F64 dz = cell_length.z * F64(k);
#endif
                        const F64 theta = (+8.0 * atan(1.0)) * (dx*mx + dy*my + dz*mz);
                        const cplx_t phase(cos(theta), sin(theta));
                        slm[ith].eval_opt(-F64(mx), F64(my), F64(mz));
                        const F64 m2 = mx*mx + my*my + mz*mz;
                        cf[ith].eval_kcut(m2, alpha);
                        for(S32 l=0; l<=2*p; l++){
                            for(S32 m=0; m<=l; m++){
                                const cplx_t val = (cf[ith].kcut[l] * phase) * slm[ith].val_at(l, m);
                                ksum[ith].accum_at(l, m, val);
                            }
                        }
                    }
                    // store sum
#ifdef PARTICLE_SIMULATOR_PMM_IGNORE_RSPACE
                    rsum[ith].clear();
#endif
#ifdef PARTICLE_SIMULATOR_PMM_IGNORE_KSPACE
                    ksum[ith].clear();
#endif
                    for(S32 lm=0; lm<LEN2; lm++){
                        gf_r(k,j,i,lm) = rsum[ith].buf[lm] + ksum[ith].buf[lm];
                    }
#ifdef PARTICLE_SIMULATOR_PMM_SHOW_RSUM_AND_KSUM
                    if(0 == (i|j|k)){
                        rsum.show();
                        ksum.show();
                    };
#endif
                } // for(k, j, i)
   
                // Free memory
                for (S32 i=0; i<n_thread; i++) {
                    cf[i].freeMem();
                    slm[i].freeMem();
                    rsum[i].freeMem();
                    ksum[i].freeMem();
                }
                delete [] cf;
                delete [] slm;
                delete [] rsum;
                delete [] ksum;
            }
    
            void doFFTOpenImpl() {
                const S32 LEN2 = (2*p+1)*(2*p+1);
                MultidimensionalArray<fftw_real_t, 3> rbuf;
                MultidimensionalArray<fftw_cplx_t, 3> kbuf;
                S32 sizes_r[3] = {2*NZ, 2*NY, 2*NX};
                S32 sizes_k[3] = {2*NZ, 2*NY, 1+NX};
                rbuf.initialize(sizes_r, 1);
                kbuf.initialize(sizes_k, 1);
                fftw_plan plan_fwd = 
                    fftw_plan_dft_r2c_3d(
                        2*NZ, 2*NY, 2*NX, 
                        (fftw_real_t *)(rbuf.getPointer()),
                        (fftw_cplx_t *)(kbuf.getPointer()),
                        FFTW_ESTIMATE);
                for(S32 lm=0; lm<LEN2; lm++){
                    S32 i, j, k;
                    for(k=0; k<2*NZ; k++) for(int j=0; j<2*NY; j++) for(i=0; i<2*NX; i++)
                    {
                        rbuf(k,j,i) = gf_r(k,j,i,lm);
                    }
                    // CALL FFTW
                    fftw_execute(plan_fwd);
        
                    for(k=0; k<2*NZ; k++) for(j=0; j<2*NY; j++) for(i=0; i<1+NX; i++)
                    {
                        gf_k(k,j,i,lm).real(kbuf(k,j,i)[0]);
                        gf_k(k,j,i,lm).imag(kbuf(k,j,i)[1]);
                        // fftw_complex = double[2]
                    }
                }
                fftw_destroy_plan(plan_fwd);
            }
    
            void doFFTPeriodicXYZImpl() {
                const S32 LEN2 = (2*p+1)*(2*p+1);
                MultidimensionalArray<fftw_real_t, 3> rbuf;
                MultidimensionalArray<fftw_cplx_t, 3> kbuf;
                S32 sizes_r[3] = {NZ, NY, NX};
                S32 sizes_k[3] = {NZ, NY, 1+NX/2};
                rbuf.initialize(sizes_r, 1);
                kbuf.initialize(sizes_k, 1);
                fftw_plan plan_fwd = 
                    fftw_plan_dft_r2c_3d(
                        NZ, NY, NX, 
                        (fftw_real_t *)(rbuf.getPointer()),
                        (fftw_cplx_t *)(kbuf.getPointer()),
                        FFTW_ESTIMATE);
                for(S32 lm=0; lm<LEN2; lm++){
                    S32 i, j, k;
                    for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
                    {
                        rbuf(k,j,i) = gf_r(k,j,i,lm);
                    }
                    // CALL FFTW
                    fftw_execute(plan_fwd);
        
                    for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<1+NX/2; i++)
                    {
                        gf_k(k,j,i,lm).real(kbuf(k,j,i)[0]);
                        gf_k(k,j,i,lm).imag(kbuf(k,j,i)[1]);
                    }
                }
                fftw_destroy_plan(plan_fwd);
            }
    
        public:
            void init(const S32 p,
                      const S32 bc,
                      const S32 NX,
                      const S32 NY,
                      const S32 NZ) {
                this->p  = p;
                this->bc = bc;
                this->NX = NX;
                this->NY = NY;
                this->NZ = NZ;
                if (bc == BOUNDARY_CONDITION_OPEN) {
                    initOpenImpl();
                } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                    initPeriodicXYZImpl();
                }
            }

            void free() {
                p = bc = NX = NY = NZ = -1;
                gf_r.freeMem();
                gf_k.freeMem();
            }
      
            void set(const S32 icut,
                     const F64vec cell_length,
                     const F64 alpha = 2.4,
                     const S32 NMAX = 3,
                     const S32 MMAX = 5) {
                if (bc == BOUNDARY_CONDITION_OPEN) {
                    setOpenImpl(icut, cell_length);
                } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                    setPeriodicXYZImpl(icut, cell_length, alpha, NMAX, MMAX);
                }
            }
        
            void doFFT(){
                if (bc == BOUNDARY_CONDITION_OPEN) {
                    doFFTOpenImpl();
                } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                    doFFTPeriodicXYZImpl();
                }
    
            }
        
            void transform(const MultidimensionalArray<cplx_t, 4> & mm_k,
                           MultidimensionalArray<cplx_t, 4> & le_k)
            {
                const S32 n_thread = Comm::getNumberOfThread();
                const S32 kmax = gf_k.getSize(0);
                const S32 jmax = gf_k.getSize(1);
                const S32 imax = gf_k.getSize(2);
                Slm<cplx_t> * slm_tmp;
                MultipoleMoment <cplx_t> * mm_tmp;
                LocalExpansion  <cplx_t> * le_tmp;
                slm_tmp = new Slm<cplx_t>[n_thread];
                mm_tmp  = new MultipoleMoment<cplx_t>[n_thread];
                le_tmp  = new LocalExpansion<cplx_t>[n_thread];
                for (S32 i=0; i<n_thread; i++) {
                    slm_tmp[i].alloc(2*p);
                    mm_tmp[i].alloc(p);
                    le_tmp[i].alloc(p);
                }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                for(S32 k=0; k<kmax; k++){
                    for(S32 j=0; j<jmax; j++){
                        for(S32 i=0; i<imax; i++){
                            const S32 ith = Comm::getThreadNum();
                            const S32 LEN = (p+1)*(p+1);
                            const S32 LEN2 = (2*p+1)*(2*p+1);
                            for (S32 lm=0; lm<LEN; lm++) {
                                mm_tmp[ith].buf[lm] = mm_k(k,j,i,lm);
                            }
                            for (S32 lm=0; lm<LEN2; lm++) {
                                slm_tmp[ith].buf[lm] = gf_k(k,j,i,lm);
                            }
                            slm_tmp[ith].transform_M2L(mm_tmp[ith], le_tmp[ith], false);
                            for (S32 lm=0; lm<LEN; lm++) {
                                le_k(k,j,i,lm) = le_tmp[ith].buf[lm];
                            }
                        }
                    }
                }
                // Free memory
                for (S32 i=0; i<n_thread; i++) {
                    slm_tmp[i].freeMem();
                    mm_tmp[i].freeMem();
                    le_tmp[i].freeMem();
                }
                delete [] slm_tmp;
                delete [] mm_tmp;
                delete [] le_tmp;
            }
        };
        
    } // END of namespace of ParticleSimulator
} // END of namespace of ParticleSimulator
