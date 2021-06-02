#pragma once
#include <fftw3.h>

namespace ParticleSimulator {
    namespace ParticleMeshMultipole {

#if defined(PARTICLE_SIMULATOR_PMM_FFTW_FLOAT)
        using fftw_real_t = float;
        using fftw_cplx_t = fftwf_complex;
#define fftw_plan            fftwf_plan
#define fftw_plan_dft_r2c_3d fftwf_plan_dft_r2c_3d
#define fftw_plan_dft_c2r_3d fftwf_plan_dft_c2r_3d 
#define fftw_execute         fftwf_execute
#define fftw_destroy_plan    fftwf_destroy_plan
#elif defined(PARTICLE_SIMULATOR_PMM_FFTW_LONG_DOUBLE)
        using fftw_real_t = long double;
        using fftw_cplx_t = fftwl_complex;
#define fftw_plan            fftwl_plan
#define fftw_plan_dft_r2c_3d fftwl_plan_dft_r2c_3d
#define fftw_plan_dft_c2r_3d fftwl_plan_dft_c2r_3d 
#define fftw_execute         fftwl_execute
#define fftw_destroy_plan    fftwl_destroy_plan
#else
        using fftw_real_t = double;
        using fftw_cplx_t = fftw_complex;
#endif


    } // END of namespace of ParticleMeshMultipole
} // END of namespace of ParticleSimulator
