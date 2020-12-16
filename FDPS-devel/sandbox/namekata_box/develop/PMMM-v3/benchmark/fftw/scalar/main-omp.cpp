// Include the standard C++ headers
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <stdlib.h>
#include <omp.h>
// Include header files of FFTW3
#include <fftw3.h>

double get_time_of_day(){
    struct timeval t;
    gettimeofday(&t,NULL);
    return (double)(t.tv_sec + t.tv_usec*1.0e-6);
}

void measure(const int NC) {

    // Allocate buffers
    double *rbuf = (double *) fftw_malloc((size_t)  sizeof(double) * NC * NC * NC);
    fftw_complex *kbuf = (fftw_complex * ) fftw_malloc((size_t) sizeof(fftw_complex) * NC * NC * (1+NC/2));

    // Create plans
    fftw_plan plan_fwd = fftw_plan_dft_r2c_3d(NC, NC, NC, rbuf, kbuf, FFTW_EXHAUSTIVE);
    fftw_plan plan_bkw = fftw_plan_dft_c2r_3d(NC, NC, NC, kbuf, rbuf, FFTW_EXHAUSTIVE);

    // Make input data
    const long seed=19810614;
    srand48(seed);
    for (int k=0; k<NC; k++) for (int j=0; j<NC; j++) for (int i=0; i<NC; i++) {
        const int idx = NC*(NC*k + j) + i;
        rbuf[idx] = 2.0 * drand48() - 1.0;
    }

    // Measure the performance
    const int n_trials = 32;
    double time_offset = get_time_of_day();
    for (int n=0; n<n_trials; n++) {
        fftw_execute(plan_fwd);
        fftw_execute(plan_bkw);
    }
    const double etime_ave = (get_time_of_day() - time_offset)/n_trials;

    // Output the result
    std::cout << "------------ Summary ------------" << std::endl;
    std::cout << " NC              = " << NC << std::endl;
    std::cout << " etime (fwd+bkw) = " << etime_ave << " [s]" << std::endl;

    // Destroy plans
    fftw_destroy_plan(plan_fwd);
    fftw_destroy_plan(plan_bkw);

    // Free 
    fftw_free(rbuf);
    fftw_free(kbuf);
}

void measure(const int NC, const int nthreads) {

    // Allocate buffers
    double **rbuf;
    fftw_complex **kbuf;
    rbuf = (double **) malloc((size_t) sizeof(double *) * nthreads);
    kbuf = (fftw_complex **) malloc((size_t) sizeof(fftw_complex *) * nthreads);
    for (int i=0; i<nthreads; i++) {
        rbuf[i] = (double *) fftw_malloc((size_t)  sizeof(double) * NC * NC * NC);
        kbuf[i] = (fftw_complex * ) fftw_malloc((size_t) sizeof(fftw_complex) * NC * NC * (1+NC/2));
    }

    // Create plans
    fftw_plan *plan_fwd, *plan_bkw;
    plan_fwd = (fftw_plan *) malloc((size_t) sizeof(fftw_plan) * nthreads); 
    plan_bkw = (fftw_plan *) malloc((size_t) sizeof(fftw_plan) * nthreads); 
    for (int i=0; i<nthreads; i++) {
        plan_fwd[i] = fftw_plan_dft_r2c_3d(NC, NC, NC, rbuf[i], kbuf[i], FFTW_ESTIMATE);
        plan_bkw[i] = fftw_plan_dft_c2r_3d(NC, NC, NC, kbuf[i], rbuf[i], FFTW_ESTIMATE);
    }

    // Make input data
    const long seed=19810614;
    srand48(seed);
    for (int n=0; n<nthreads; n++) {
        for (int k=0; k<NC; k++) for (int j=0; j<NC; j++) for (int i=0; i<NC; i++) {
            const int idx = NC*(NC*k + j) + i;
            rbuf[n][idx] = 2.0 * drand48() - 1.0;
        }
    }

    // Measure the performance
    omp_set_num_threads(nthreads); // set # of threads used in the parallel region
    const int nFFTs = 64;
    double time_offset = get_time_of_day();
#pragma omp parallel for
    for (int n=0; n<nFFTs; n++) {
        const int th = omp_get_thread_num();
        fftw_execute(plan_fwd[th]);
        fftw_execute(plan_bkw[th]);
    }
    const double etime = (get_time_of_day() - time_offset);
    omp_set_num_threads(omp_get_max_threads()); // reset

    // Output the result
    std::cout << "------------ Summary ------------" << std::endl;
    std::cout << " NC              = " << NC << std::endl;
    std::cout << " nFFTs           = " << nFFTs << std::endl;
    std::cout << " etime (fwd+bkw) = " << etime << " [s]" << std::endl;

    // Destroy plans
    for (int i=0; i<nthreads; i++) {
        fftw_destroy_plan(plan_fwd[i]);
        fftw_destroy_plan(plan_bkw[i]);
    }
    free(plan_fwd);
    free(plan_bkw);

    // Free 
    for (int i=0; i<nthreads; i++) {
        fftw_free(rbuf[i]);
        fftw_free(kbuf[i]);
    }
    free(rbuf);
    free(kbuf);
}


int main(int argc, char argv[]) {

    fftw_init_threads();
    const int n_thrd_max = omp_get_max_threads();

    std::cout << "================================" << std::endl;
    std::cout << " Test1: multithreaded FFT"        << std::endl;
    std::cout << "================================" << std::endl;
    for (int n_thrd=2; n_thrd <= n_thrd_max; n_thrd*=2) {
        std::cout << "n_thrd = " << n_thrd << std::endl;
        fftw_plan_with_nthreads(n_thrd);
        for (int n=8; n<=256; n *= 2) {
            measure(n);
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "================================" << std::endl;
    std::cout << " Test2: multipole FFT"            << std::endl;
    std::cout << "================================" << std::endl;
    fftw_plan_with_nthreads(1); // reset
    for (int n_thrd=1; n_thrd <= n_thrd_max; n_thrd*=2) {
        std::cout << "n_thrd = " << n_thrd << std::endl;
        for (int n=8; n<=256; n *= 2) {
            measure(n,n_thrd);
        }
        std::cout << std::endl;
    }

    fftw_cleanup_threads();

    return 0;
}
