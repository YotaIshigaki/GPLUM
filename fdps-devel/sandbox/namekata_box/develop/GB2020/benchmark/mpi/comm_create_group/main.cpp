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

    int n_proc_in_group = n_proc/2; 
    // Get the parent group
    MPI_Group parent_group;
    MPI_Comm_group(MPI_COMM_WORLD,&parent_group);
    // Create a group of n_proc processes
    std::vector<int> ranks;
    ranks.reserve(n_proc_in_group);
    for (int i=0; i<n_proc_in_group; i++) ranks[i]=i;
    MPI_Group group;
    if (my_rank < n_proc_in_group) {
        MPI_Group_incl(parent_group, n_proc_in_group, &ranks[0], &group);
    } else {
        group = MPI_GROUP_NULL;
    }
    // Measure a time to create and destroy a new comm
    const int n_trials = 32;
    MPI_Barrier(MPI_COMM_WORLD);
    double time_offset = get_time_of_day();
    if (group != MPI_GROUP_NULL) {
        for (int i = 0; i<n_trials; i++) {
            // Create a communicator corresponding to this group
            MPI_Comm comm;
            int tag = 0;
            int ret = MPI_Comm_create_group(MPI_COMM_WORLD, group, tag, &comm);
            // Free the communicator
            MPI_Comm_free(&comm); 
        }
        // Free the group
        MPI_Group_free(&group);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double etime = (get_time_of_day() - time_offset)/n_trials;

    // Output the measured data
    if (my_rank == 0) std::cout << "etime = " << etime << std::endl;

    // Finalize
    MPI_Finalize();
    return 0;
}
