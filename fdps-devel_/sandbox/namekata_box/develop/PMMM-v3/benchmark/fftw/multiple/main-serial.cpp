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

double get_time_of_FFT(const int NC, const int NFFT, const int mode) {

    if (mode == 0) {
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
            const int idx = (2*(NC/2+1))*(NC*k + j) + i;
            rbuf[idx] = 2.0 * drand48() - 1.0;
        }

        // Measure the performance
        double time_offset = get_time_of_day();
        for (int n=0; n<NFFT; n++) {
            fftw_execute(plan_fwd);
            fftw_execute(plan_bkw);
        }
        const double etime_ave = (get_time_of_day() - time_offset);

        // Free 
        fftw_free(rbuf);
        fftw_free(kbuf);

        return etime_ave;

    } else if (mode == 1) {
        // Allocate buffers
        double *rbuf = (double *) fftw_malloc((size_t)  sizeof(double) * NFFT * NC * NC * NC);
        fftw_complex *kbuf = (fftw_complex * ) fftw_malloc((size_t) sizeof(fftw_complex) * NFFT * NC * NC * (1+NC/2));

        // Create plans
        const int rank = 3;
        const int sizes[] = {NC, NC, NC};
        const int howmany = NFFT;
        const int istride = 1;
        const int ostride = 1;
        const int idist = NC*NC*NC;
        const int odist = NC*NC*(1+NC/2);
        fftw_plan plan_fwd = fftw_plan_many_dft_r2c(rank, sizes, howmany,
                                                    rbuf, NULL, istride, idist,
                                                    kbuf, NULL, ostride, odist,
                                                    FFTW_EXHAUSTIVE);
        fftw_plan plan_bkw = fftw_plan_many_dft_c2r(rank, sizes, howmany,
                                                    kbuf, NULL, ostride, odist,
                                                    rbuf, NULL, istride, idist,
                                                    FFTW_EXHAUSTIVE);

        // Make input data
        const long seed=19810614;
        srand48(seed);
        for (int n=0; n<NFFT; n++)
        for (int k=0; k<NC; k++)
        for (int j=0; j<NC; j++)
        for (int i=0; i<NC; i++) {
            const int idx = NC*(NC*(NC*n + k) + j) + i;
            rbuf[idx] = 2.0 * drand48() - 1.0;
        }

        // Measure the performance
        double time_offset = get_time_of_day();
        fftw_execute(plan_fwd);
        fftw_execute(plan_bkw);
        const double etime_ave = (get_time_of_day() - time_offset);

        // Free 
        fftw_free(rbuf);
        fftw_free(kbuf);

        return etime_ave;

    } else if (mode == 2) {
        // Allocate buffers
        double *rbuf = (double *) fftw_malloc((size_t)  sizeof(double) * NFFT * NC * NC * NC);
        fftw_complex *kbuf = (fftw_complex * ) fftw_malloc((size_t) sizeof(fftw_complex) * NFFT * NC * NC * (1+NC/2));

        // Create plans
        const int rank = 3;
        const int sizes[3] = {NC, NC, NC};
        const int howmany = NFFT;
        const int istride = NFFT;
        const int ostride = NFFT;
        const int idist = 1;
        const int odist = 1;
        fftw_plan plan_fwd = fftw_plan_many_dft_r2c(rank, sizes, howmany,
                                                    rbuf, NULL, istride, idist,
                                                    kbuf, NULL, ostride, odist,
                                                    FFTW_EXHAUSTIVE);
        fftw_plan plan_bkw = fftw_plan_many_dft_c2r(rank, sizes, howmany,
                                                    kbuf, NULL, istride, idist,
                                                    rbuf, NULL, ostride, odist,
                                                    FFTW_EXHAUSTIVE);

        // Make input data
        const long seed=19810614;
        srand48(seed);
        for (int k=0; k<NC; k++)
        for (int j=0; j<NC; j++)
        for (int i=0; i<NC; i++) 
        for (int n=0; n<NFFT; n++) {
            const int idx = NFFT*(NC*(NC*k + j) + i) + n;
            rbuf[idx] = 2.0 * drand48() - 1.0;
        }

        // Measure the performance
        double time_offset = get_time_of_day();
        fftw_execute(plan_fwd);
        fftw_execute(plan_bkw);
        const double etime_ave = (get_time_of_day() - time_offset);

        // Free 
        fftw_free(rbuf);
        fftw_free(kbuf);

        return etime_ave;

    } else {
        return -1.0; // representing an error.
    }
}


int main(int argc, char *argv[]) {

    const int NC = 64;
    const int NFFT = 32;
    std::cout << "etime0 = " << get_time_of_FFT(NC,NFFT,0) << std::endl;
    std::cout << "etime1 = " << get_time_of_FFT(NC,NFFT,1) << std::endl;
    std::cout << "etime2 = " << get_time_of_FFT(NC,NFFT,2) << std::endl;
    

    return 0;
}
