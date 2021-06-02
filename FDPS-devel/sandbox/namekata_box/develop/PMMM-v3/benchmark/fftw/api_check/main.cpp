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

int main(int argc, char *argv[]) {
   
    int n_proc_max, my_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&n_proc_max);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    fftw_mpi_init();

    // Set problem size
    const int NX = 64, NY = 81, NZ = 125;
    const int tgt_rank = 2;

    // Allocate buffers
    ptrdiff_t alloc_local, local_n0, local_0_start;
    ptrdiff_t local_n1, local_1_start;
    ptrdiff_t nelem_rbuf, nelem_kbuf;
    // (1) 
    alloc_local = fftw_mpi_local_size_3d_transposed(NZ, NY, NX/2+1, MPI_COMM_WORLD,
                                                    &local_n0, &local_0_start,
                                                    &local_n1, &local_1_start);
    nelem_rbuf = 2 * alloc_local;
    nelem_kbuf = alloc_local;
    if (my_rank == tgt_rank) {
        std::cout << "[1] my_rank = " << my_rank
                  << " local_0_start = " << local_0_start
                  << " local_n0 = " << local_n0
                  << " local_1_start = " << local_1_start
                  << " local_n1 = " << local_n1
                  << " nelem_rbuf = " << nelem_rbuf
                  << " nelem_kbuf = " << nelem_kbuf
                  << std::endl;
                
    }
    // (2)
    alloc_local = fftw_mpi_local_size_3d_transposed(NY, NZ, NX/2+1, MPI_COMM_WORLD,
                                                    &local_n0, &local_0_start,
                                                    &local_n1, &local_1_start);
    nelem_rbuf = 2 * alloc_local;
    nelem_kbuf = alloc_local;
    if (my_rank == tgt_rank) {
        std::cout << "[2] my_rank = " << my_rank
                  << " local_0_start = " << local_0_start
                  << " local_n0 = " << local_n0
                  << " local_1_start = " << local_1_start
                  << " local_n1 = " << local_n1
                  << " nelem_rbuf = " << nelem_rbuf
                  << " nelem_kbuf = " << nelem_kbuf
                  << std::endl;
                
    }

    // Finalize
    fftw_mpi_cleanup();
    MPI_Finalize();

    return 0;
}
