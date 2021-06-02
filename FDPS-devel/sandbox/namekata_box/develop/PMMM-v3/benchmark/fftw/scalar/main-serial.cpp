// Include the standard C++ headers
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <stdlib.h>
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
    std::cout << " n_trials        = " << n_trials << std::endl;
    std::cout << " etime (fwd+bkw) = " << etime_ave << " [s]" << std::endl;

    // Free 
    fftw_free(rbuf);
    fftw_free(kbuf);
}


int main(int argc, char argv[]) {

    for (int n=8; n<=512; n *= 2) {
        measure(n);
    }

    return 0;
}
