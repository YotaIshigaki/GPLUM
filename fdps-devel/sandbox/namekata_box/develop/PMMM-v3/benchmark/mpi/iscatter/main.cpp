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

    for (int i=0; i < data_size.size(); ++i) {
        // Output
        if (my_rank == 0) {
            std::cout << "Measuring the case of data size of " << data_size[i] << " [B]" << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // Set sendbuf
        const int sendcount = data_size[i];
        char *sendbuf = new char[sendcount * n_proc];
        // Set recvbuf
        const int recvcount = sendcount;
        char *recvbuf = new char[recvcount];
        // Measure the performance
        MPI_Request req;
        MPI_Status stat;
        const int n_trials = 32;
        const double time_offset = get_time_of_day();
        for (int n=0; n<n_trials; ++n) {
            MPI_Iscatter(sendbuf, sendcount, MPI_CHAR,
                         recvbuf, recvcount, MPI_CHAR,
                         0,MPI_COMM_WORLD,&req);
            MPI_Wait(&req, &stat);
        }
        if (my_rank == 0) etime[sendcount] = (get_time_of_day() - time_offset)/n_trials;
        // Free
        delete [] sendbuf;
        delete [] recvbuf;
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
