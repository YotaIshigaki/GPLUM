// Include the standard C++ headers
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <vector>
#include <unordered_map>
#include <stdlib.h>
// Include header files of OpenMP
#include <omp.h>
// Include header files of MPI
#include "mpi.h"
// Include header files of FFTW3
#include <fftw3-mpi.h>

#define DEBUG

double get_time_of_day(){
    struct timeval t;
    gettimeofday(&t,NULL);
    return (double)(t.tv_sec + t.tv_usec*1.0e-6);
}

double get_time_of_FFT(const int NC, MPI_Comm & comm) {

#if !defined(DEBUG)
    // Allocate buffers
    ptrdiff_t alloc_local, local_n0, local_0_start;
    ptrdiff_t local_n1, local_1_start;
    alloc_local = fftw_mpi_local_size_3d_transposed(NC, NC, NC/2+1, comm,
                                                    &local_n0, &local_0_start,
                                                    &local_n1, &local_1_start);
    double *rbuf = fftw_alloc_real(2 * alloc_local);
    fftw_complex *kbuf = fftw_alloc_complex(alloc_local);

    // Create plans
    fftw_plan plan_fwd = fftw_mpi_plan_dft_r2c_3d(NC, NC, NC, rbuf, kbuf,
                                                  comm,FFTW_EXHAUSTIVE || FFTW_MPI_TRANSPOSED_OUT);
    fftw_plan plan_bkw = fftw_mpi_plan_dft_c2r_3d(NC, NC, NC, kbuf, rbuf,
                                                  comm,FFTW_EXHAUSTIVE || FFTW_MPI_TRANSPOSED_IN);

    // Make input data
    const long seed=19810614;
    srand48(seed);
    for (int k=0; k<local_n0; k++) for (int j=0; j<NC; j++) for (int i=0; i<NC; i++) {
        const int idx = (2*(NC/2+1))*(NC*k + j) + i;
        rbuf[idx] = 2.0 * drand48() - 1.0;
    }
#else // DEBUG
    // Get rank number
    int my_rank;
    MPI_Comm_rank(comm,&my_rank);

    // Set the problem size
    const int NX = 64, NY = 81, NZ = 125;

    // Allocate buffers
    ptrdiff_t alloc_local, local_n0, local_0_start;
    ptrdiff_t local_n1, local_1_start;
    alloc_local = fftw_mpi_local_size_3d_transposed(NZ, NY, NX/2+1, comm,
                                                    &local_n0, &local_0_start,
                                                    &local_n1, &local_1_start);
    std::cout << "my_rank = " << my_rank 
              << " local_n0 = " << local_n0 
              << " local_0_start = " << local_0_start
              << " local_n1 = " << local_n1
              << " local_1_start = " << local_1_start
              << std::endl;
    double *rbuf = fftw_alloc_real(2 * alloc_local);
    double *rbuf_org = fftw_alloc_real(2 * alloc_local);
    fftw_complex *kbuf = fftw_alloc_complex(alloc_local);

    // Create plans
    fftw_plan plan_fwd = fftw_mpi_plan_dft_r2c_3d(NZ, NY, NX, rbuf, kbuf,
                                                  comm,FFTW_EXHAUSTIVE || FFTW_MPI_TRANSPOSED_OUT);
    fftw_plan plan_bkw = fftw_mpi_plan_dft_c2r_3d(NZ, NY, NX, kbuf, rbuf,
                                                  comm,FFTW_EXHAUSTIVE || FFTW_MPI_TRANSPOSED_IN);

    // Make input data
    const long seed=19810614;
    srand48(seed);
    for (int k=0; k<local_n0; k++) for (int j=0; j<NY; j++) for (int i=0; i<NX; i++) {
        const int idx = (2*(NX/2+1))*(NY*k + j) + i;
        rbuf[idx] = 2.0 * drand48() - 1.0;
        rbuf_org[idx] = rbuf[idx];
    }
#endif // DEBUG

    // Measure the performance
#ifdef DEBUG
    const int n_trials = 1;
#else
    const int n_trials = 32;
#endif
    const double time_offset = get_time_of_day();
    for (int n=0; n<n_trials; n++) {
        fftw_execute(plan_fwd);
        fftw_execute(plan_bkw);
    }
    const double etime_ave = (get_time_of_day() - time_offset)/n_trials;

#ifdef DEBUG
    std::cout << "etime_ave = " << etime_ave << std::endl;
    for (int k=0; k<local_n0; k++) for (int j=0; j<NY; j++) for (int i=0; i<NX; i++) {
        const int idx = (2*(NX/2+1))*(NY*k + j) + i;
        rbuf[idx] /= (NX*NY*NZ);
        const double diff = (rbuf[idx] - rbuf_org[idx])/rbuf_org[idx];
        if (std::fabs(diff) > 1.0e-10) {
            std::cout << "i = " << i
                      << " j = " << j
                      << " k = " << k
                      << " diff = " << diff
                      << std::endl;
        }
    }
    MPI_Finalize();
    std::exit(0);
#endif

    // Free 
    fftw_free(rbuf);
    fftw_free(kbuf);

    // Return value
    return etime_ave;
}

int main(int argc, char *argv[]) {
   
    int n_proc_max, my_rank, provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    const bool threads_ok = (provided >= MPI_THREAD_FUNNELED);
    MPI_Comm_size(MPI_COMM_WORLD,&n_proc_max);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    if (my_rank == 0) std::cout << "threads_ok = " << threads_ok << std::endl;
    if (threads_ok) {
        fftw_init_threads();
        fftw_plan_with_nthreads(omp_get_max_threads());
    }
    fftw_mpi_init();

    std::unordered_map<int, std::unordered_map<int, double> > data; // NC, n_proc, etime

    for (int n_proc=2; n_proc <= n_proc_max; n_proc *= 2) {
        // Output
        if (my_rank == 0) {
            std::cout << "Measuring the case of n_proc = " << n_proc << std::endl;
        }
        // Get the parent group
        MPI_Group parent_group;
        MPI_Comm_group(MPI_COMM_WORLD,&parent_group);
        // Create a group of n_proc processes
        std::vector<int> ranks;
        ranks.reserve(n_proc);
        for (int i=0; i<n_proc; i++) ranks[i]=i;
        MPI_Group group;
        MPI_Group_incl(parent_group, n_proc, &ranks[0], &group);
        // Create a communicator corresponding to this group
        MPI_Comm comm;
        MPI_Comm_create(MPI_COMM_WORLD, group, &comm);
        // Measure the speed of FFT
        if (comm != MPI_COMM_NULL) { 
            for (int NC=32; NC<=512; NC *= 2) {
                if (n_proc <= NC) {
                    const double etime = get_time_of_FFT(NC,comm);
                    if (my_rank == 0) data[NC][n_proc] = etime;
                }
            }
        }
        // Free the communicator
        if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
        // Free the group
        MPI_Group_free(&group);
    }

    // Output the measured data
    if (my_rank == 0) {
        for (auto itr_outer = data.begin(); itr_outer != data.end(); ++itr_outer) {
            const int NC = itr_outer->first;
            const long long int NC3 = NC*NC*NC;
            std::cout << "# NC = " << NC << "    NC3 = " << NC3 << std::endl; // NC, NC3
            std::cout << "# (n_proc, etime)" << std::endl;
            for (auto itr_inner = itr_outer->second.begin(); itr_inner != itr_outer->second.end(); ++itr_inner) {
                std::cout << itr_inner->first << "    " << itr_inner->second << std::endl; // (n_proc, etime) 
            }
            std::cout << std::endl << std::endl << std::endl;
        }
    }

    // Finalize
    fftw_mpi_cleanup();
    MPI_Finalize();

    return 0;
}
