// Include the standard C++ headers
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <new>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <vector>
#include <unordered_map>
#include <stdlib.h>
#include <algorithm>
// Include header files of MPI
#include "mpi.h"

double get_time_of_day(){
    struct timeval t;
    gettimeofday(&t,NULL);
    return (double)(t.tv_sec + t.tv_usec*1.0e-6);
}

int main(int argc, char *argv[]) {
   
    int n_proc, my_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

    // Data size
    std::vector<int> data_size; // in bytes
    data_size.push_back(8);
    data_size.push_back(32);
    data_size.push_back(128);
    data_size.push_back(1024);
    data_size.push_back(4*1024);
    data_size.push_back(16*1024);
    data_size.push_back(64*1024);
    data_size.push_back(256*1024);
    data_size.push_back(1024*1024);

    std::unordered_map<int, double> etime; 

    // Set the seed
    const long seed=19810614;
    srand48(seed);

    for (int i=0; i < data_size.size(); ++i) {
        // Output
        if (my_rank == 0) {
            std::cout << "Measuring the case of data size of " << data_size[i] << " [B]" << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // Measure the performance for different configurations
        const int n_trials_outer=64, n_trials_inner=4;
        for (int nn=0; nn<n_trials_outer; ++nn) {
            // Set recvcounts
            double *fractions = new(std::nothrow) double[n_proc];
            if (fractions == nullptr) {
                std::cerr << "new for fractions[] is failed!" << std::endl;
                MPI_Abort(MPI_COMM_WORLD,-1);
                std::exit(-1);
            }
            const int recvcount_tot = data_size[i] * n_proc;
            int *recvcounts = new(std::nothrow) int[n_proc];
            if (recvcounts == nullptr) {
                std::cerr << "new for recvcounts[] is failed!" << std::endl;
                MPI_Abort(MPI_COMM_WORLD,-1);
                std::exit(-1);
            }
            // Set fractions using a pseudo-random number generator
            double fsum = 0.0;
            for (int k=0; k<n_proc; k++) {
                double f = drand48();
                fractions[k] = f;
                fsum += f;
            }
            for (int k=0; k<n_proc; k++) fractions[k] /= fsum; // scaling for \sum(fractions)=1
            // Set recvcounts
            int csum {0};
            for (int k=0; k<n_proc; k++) {
                int cnt = recvcount_tot * fractions[k];
                recvcounts[k] = cnt;
                csum += cnt;
            }
            if (csum != recvcount_tot) {
                int diff = recvcount_tot - csum;
                assert(0 < diff && diff <= n_proc);
                for (int k=0; k<diff; k++) recvcounts[k] += 1;
            }
            // Check if recvcounts[] is consistent with recvcount_tot
            csum = 0;
            for (int k=0; k<n_proc; k++) csum += recvcounts[k];
            assert(csum == recvcount_tot);
            // Set recvdispls
            int *recvdispls = new(std::nothrow) int[n_proc];
            if (recvdispls == nullptr) {
                std::cerr << "new for recvdispls[] is failed!" << std::endl;
                MPI_Abort(MPI_COMM_WORLD,-1);
                std::exit(-1);
            }
            recvdispls[0] = 0;
            for (int k=1; k<n_proc; k++) 
                recvdispls[k] = recvdispls[k-1] + recvcounts[k-1];
            // Set recvbuf
            char *recvbuf = new(std::nothrow) char[recvcount_tot];
            if (recvbuf == nullptr) {
                std::cerr << "new for recvbuf[] is failed!" << std::endl;
                MPI_Abort(MPI_COMM_WORLD,-1);
                std::exit(-1);
            }
            // Set sendbuf
            const int sendcount = recvcounts[my_rank];
            char *sendbuf;
            if (sendcount > 0) {
                sendbuf = new(std::nothrow) char[sendcount];
                if (sendbuf == nullptr) {
                    std::cerr << "new for sendbuf[] is failed!" << std::endl;
                    MPI_Abort(MPI_COMM_WORLD,-1);
                    std::exit(-1);
                }
            } else sendbuf = nullptr;
            // Measure the performance
            MPI_Request req;
            MPI_Status stat;
            const double time_offset = get_time_of_day();
            for (int n=0; n<n_trials_inner; ++n) {
                MPI_Igatherv(sendbuf, sendcount, MPI_CHAR,
                             recvbuf, recvcounts, recvdispls, MPI_CHAR,
                             0,MPI_COMM_WORLD,&req);
                MPI_Wait(&req, &stat);
            }
            const double et = get_time_of_day() - time_offset;
            if (my_rank == 0) {
                if (nn == 0) etime[data_size[i]] = et;
                else etime[data_size[i]] += et;
            }
            // Free
            delete [] fractions;
            delete [] recvcounts;
            delete [] recvdispls;
            delete [] recvbuf;
            if (sendbuf != nullptr) delete [] sendbuf;
        }
        if (my_rank == 0) etime[data_size[i]] /= (n_trials_inner * n_trials_outer);
    }

    // Output
    if (my_rank == 0) {
        std::cout << "# data size [B],     etime [s]" << std::endl;
        for (auto itr = etime.begin(); itr != etime.end(); itr++) {
            std::cout << itr->first << "    " << itr->second << std::endl;
        }
    }

    // Finalize
    MPI_Finalize();

    return 0;
}
