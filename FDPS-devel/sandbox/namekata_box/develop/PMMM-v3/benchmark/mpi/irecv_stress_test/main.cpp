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
    
    // Set sendbuf
    //constexpr int sendcount = 100; // 100[B]
    //constexpr int sendcount = 10000; // 10[kB]
    constexpr int sendcount = 100000; // 100[kB]
    //constexpr int sendcount = 1000000; // 1[MB]
    char *sendbuf = new char[sendcount];
    // Set recvbuf
    constexpr int recvcount = sendcount;
    char *recvbuf = new char[recvcount * n_proc];
    // Output information
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == 0) std::cout << "Stress test of MPI_Irecv started." << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    const double time_offset = get_time_of_day();
    // Perform MPI_Irecv
    MPI_Request *req = new MPI_Request[n_proc];
    for (int i=0; i<n_proc; i++) {
        constexpr int recvtag = 0;
        const int adr = recvcount * i;
        MPI_Irecv(&recvbuf[adr], recvcount, MPI_CHAR, i, recvtag, MPI_COMM_WORLD, &req[i]);
    }
    // Perform MPI_Send
    for (int i=0; i<n_proc; i++) {
        constexpr int sendtag = 0;
        MPI_Send(sendbuf, sendcount, MPI_CHAR, i, sendtag, MPI_COMM_WORLD);
    }
    // Perform MPI_Waitall
    MPI_Status *stat = new MPI_Status[n_proc];
    MPI_Waitall(n_proc, req, stat);
    // Output information
    MPI_Barrier(MPI_COMM_WORLD);
    const double etime = get_time_of_day() - time_offset;
    if (my_rank == 0) {
        std::cout << "Stress test of MPI_Irecv ended." << std::endl;
        std::cout << "time required was " << etime << " [s]." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Free
    delete [] sendbuf;
    delete [] recvbuf;
    delete [] req;
    delete [] stat;

    // Finalize
    MPI_Finalize();

    return 0;
}
