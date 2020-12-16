#pragma once
#include "multidimensional_array.hpp"
#include "green_function.hpp"

namespace ParticleSimulator {
    namespace ParticleMeshMultipole {

        template<typename real_t, typename cplx_t>
        static void M2LConvolutionOpenImpl(GreenFunction<real_t, cplx_t> & gf,
                                           MultidimensionalArray<Cell_FMM<real_t, cplx_t>, 3> & cell) {
            const S32 p = gf.p;
            const S32 LEN = (p+1)*(p+1);
            // Get # of grid points 
            const S32 NX = gf.NX;
            const S32 NY = gf.NY;
            const S32 NZ = gf.NZ;
            // Multipole Moments
            MultidimensionalArray<cplx_t, 4> mm_k;
            // Local Expansions
            MultidimensionalArray<cplx_t, 4> le_k;
            // FFT buffer
            MultidimensionalArray<fftw_real_t, 3> rbuf;
            MultidimensionalArray<fftw_cplx_t, 3> kbuf;
            // Memory allocation
            S32 sizes_mm_k[4] = {2*NZ, 2*NY, 1+NX, LEN};
            S32 sizes_le_k[4] = {2*NZ, 2*NY, 1+NX, LEN};
            S32 sizes_rbuf[3] = {2*NZ, 2*NY, 2*NX};
            S32 sizes_kbuf[3] = {2*NZ, 2*NY, 1+NX};
            mm_k.initialize(sizes_mm_k);
            le_k.initialize(sizes_le_k);
            rbuf.initialize(sizes_rbuf, 1);
            kbuf.initialize(sizes_kbuf, 1);
            // Create plans of FFTW
            fftw_plan plan_fwd = 
                fftw_plan_dft_r2c_3d(
                    2*NZ, 2*NY, 2*NX, 
                    (fftw_real_t *)(rbuf.getPointer()),
                    (fftw_cplx_t *)(kbuf.getPointer()),
                    FFTW_ESTIMATE);
            fftw_plan plan_bkw  =
                fftw_plan_dft_c2r_3d(
                    2*NZ, 2*NY, 2*NX, 
                    (fftw_cplx_t *)(kbuf.getPointer()),
                    (fftw_real_t *)(rbuf.getPointer()),
                    FFTW_ESTIMATE);
        
            S32 i, j, k;
            // clear rbuf
            for(k=0; k<2*NZ; k++) for(j=0; j<2*NY; j++) for(i=0; i<2*NX; i++)
            {
                rbuf(k,j,i) = 0.0;
            }
            // forward multipole
            for(S32 lm=0; lm<LEN; lm++){
                for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
                {
                    rbuf(k,j,i) = cell(k,j,i).mm.buf[lm];
                }
        
                fftw_execute(plan_fwd);
        
                for(k=0; k<2*NZ; k++) for(j=0; j<2*NY; j++) for(i=0; i<1+NX; i++)
                {
                    mm_k(k,j,i,lm).real(kbuf(k,j,i)[0]);
                    mm_k(k,j,i,lm).imag(kbuf(k,j,i)[1]);
                }
            }
            // M2L transformation
            gf.transform(mm_k, le_k);
            // backward local expansion
            for(S32 lm=0; lm<LEN; lm++){
                for(k=0; k<2*NZ; k++) for(j=0; j<2*NY; j++) for(i=0; i<1+NX; i++)
                {
                    kbuf(k,j,i)[0] = le_k(k,j,i,lm).real();
                    kbuf(k,j,i)[1] = le_k(k,j,i,lm).imag();
                }
    
                fftw_execute(plan_bkw);
    
                const F64 norm = 1.0 / (8*NX*NY*NZ);
                for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
                {
                    cell(k,j,i).le.buf[lm] = norm * rbuf(k,j,i);
                }
            }
            fftw_destroy_plan(plan_fwd);
            fftw_destroy_plan(plan_bkw);
        }
    
        template <typename real_t, typename cplx_t>
        static void M2LConvolutionPeriodicXYZImpl(GreenFunction<real_t, cplx_t> & gf,
                                                  MultidimensionalArray<Cell_FMM<real_t, cplx_t>, 3> & cell) {
            const S32 p = gf.p;
            const S32 LEN = (p+1)*(p+1);
            // Get # of grid points
            const S32 NX = gf.NX;
            const S32 NY = gf.NY;
            const S32 NZ = gf.NZ; 
            // Multipole Moments
            MultidimensionalArray<cplx_t, 4> mm_k;
            // Local Expansions
            MultidimensionalArray<cplx_t, 4> le_k;
            // FFT buffer
            MultidimensionalArray<fftw_real_t, 3> rbuf;
            MultidimensionalArray<fftw_cplx_t, 3> kbuf;
            // Memory allocation
            S32 sizes_mm_k[4] = {NZ, NY, 1+NX/2, LEN};
            S32 sizes_le_k[4] = {NZ, NY, 1+NX/2, LEN};
            S32 sizes_rbuf[3] = {NZ, NY, NX};
            S32 sizes_kbuf[3] = {NZ, NY, 1+NX/2};
            mm_k.initialize(sizes_mm_k);
            le_k.initialize(sizes_le_k);
            rbuf.initialize(sizes_rbuf, 1);
            kbuf.initialize(sizes_kbuf, 1);
            // Create plans of FFTW 
            fftw_plan plan_fwd = 
                fftw_plan_dft_r2c_3d(
                    NZ, NY, NX, 
                    (fftw_real_t *)(rbuf.getPointer()),
                    (fftw_cplx_t *)(kbuf.getPointer()),
                    FFTW_ESTIMATE);
            fftw_plan plan_bkw  =
                fftw_plan_dft_c2r_3d(
                    NZ, NY, NX, 
                    (fftw_cplx_t *)(kbuf.getPointer()),
                    (fftw_real_t *)(rbuf.getPointer()),
                    FFTW_ESTIMATE);
            S32 i, j, k;
            // forward multipole
            for(S32 lm=0; lm<LEN; lm++){
                for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
                {
                    rbuf(k,j,i) = cell(k,j,i).mm.buf[lm];
                }
        
                fftw_execute(plan_fwd);
        
                for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<1+NX/2; i++)
                {
                    mm_k(k,j,i,lm).real(kbuf(k,j,i)[0]);
                    mm_k(k,j,i,lm).imag(kbuf(k,j,i)[1]);
                }
            }
            // M2L transformation
            gf.transform(mm_k, le_k);
            // backward local expansion
            for(S32 lm=0; lm<LEN; lm++){
                for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<1+NX/2; i++)
                {
                    kbuf(k,j,i)[0] = le_k(k,j,i,lm).real();
                    kbuf(k,j,i)[1] = le_k(k,j,i,lm).imag();
                }
        
                fftw_execute(plan_bkw);
        
                const F64 norm = 1.0 / (NX*NY*NZ);
                for(k=0; k<NZ; k++) for(j=0; j<NY; j++) for(i=0; i<NX; i++)
                {
                    cell(k,j,i).le.buf[lm] = norm * rbuf(k,j,i);
                }
        
            }
            fftw_destroy_plan(plan_fwd);
            fftw_destroy_plan(plan_bkw);
        }
    
        template<typename real_t, typename cplx_t>
        void M2LConvolution(GreenFunction<real_t, cplx_t> & gf,
                            MultidimensionalArray<Cell_FMM<real_t, cplx_t>, 3> & cell) {
            if (gf.bc == BOUNDARY_CONDITION_OPEN) {
                M2LConvolutionOpenImpl<real_t, cplx_t>(gf, cell);
            } else if (gf.bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                M2LConvolutionPeriodicXYZImpl<real_t, cplx_t>(gf, cell);
            }
        }
    } // END of namespace of ParticleMeshMultipole
} // END of namespace of ParticleSimulator
