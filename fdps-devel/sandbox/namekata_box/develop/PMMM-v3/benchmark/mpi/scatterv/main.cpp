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
            // Set sendcounts
            double *fractions = new(std::nothrow) double[n_proc];
            if (fractions == nullptr) {
                std::cerr << "new for fractions[] is failed!" << std::endl;
                MPI_Abort(MPI_COMM_WORLD,-1);
                std::exit(-1);
            }
            const int sendcount_tot = data_size[i] * n_proc;
            int *sendcounts = new(std::nothrow) int[n_proc];
            if (sendcounts == nullptr) {
                std::cerr << "new for sendcounts[] is failed!" << std::endl;
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
            // Set sendcounts
            int csum {0};
            for (int k=0; k<n_proc; k++) {
                int cnt = sendcount_tot * fractions[k];
                sendcounts[k] = cnt;
                csum += cnt;
            }
            if (csum != sendcount_tot) {
                int diff = sendcount_tot - csum;
                assert(0 < diff && diff <= n_proc);
                for (int k=0; k<diff; k++) sendcounts[k] += 1;
            }
            // Check if sendcounts[] is consistent with sendcount_tot
            csum = 0;
            for (int k=0; k<n_proc; k++) csum += sendcounts[k];
            assert(csum == sendcount_tot);
            // Set senddispls
            int *senddispls = new(std::nothrow) int[n_proc];
            if (senddispls == nullptr) {
                std::cerr << "new for senddispls[] is failed!" << std::endl;
                MPI_Abort(MPI_COMM_WORLD,-1);
                std::exit(-1);
            }
            senddispls[0] = 0;
            for (int k=1; k<n_proc; k++) 
                senddispls[k] = senddispls[k-1] + sendcounts[k-1];
            // Set sendbuf
            char *sendbuf = new(std::nothrow) char[sendcount_tot];
            if (sendbuf == nullptr) {
                std::cerr << "new for sendbuf[] is failed!" << std::endl;
                MPI_Abort(MPI_COMM_WORLD,-1);
                std::exit(-1);
            }
            // Set recvbuf
            const int recvcount = sendcounts[my_rank];
            char *recvbuf;
            if (recvcount > 0) {
                recvbuf = new(std::nothrow) char[recvcount];
                if (recvbuf == nullptr) {
                    std::cerr << "new for recvbuf[] is failed!" << std::endl;
                    MPI_Abort(MPI_COMM_WORLD,-1);
                    std::exit(-1);
                }
            } else recvbuf = nullptr;
            // Measure the performance
            const double time_offset = get_time_of_day();
            for (int n=0; n<n_trials_inner; ++n) {
                MPI_Scatterv(sendbuf, sendcounts, senddispls, MPI_CHAR,
                             recvbuf, recvcount, MPI_CHAR,
                             0,MPI_COMM_WORLD);
            }
            const double et = get_time_of_day() - time_offset;
            if (my_rank == 0) {
                if (nn == 0) etime[data_size[i]] = et;
                else etime[data_size[i]] += et;
            }
            // Free
            delete [] fractions;
            delete [] sendcounts;
            delete [] senddispls;
            delete [] sendbuf;
            if (recvbuf != nullptr) delete [] recvbuf;
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
