// Include the standard C++ headers
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <vector>
#include <unordered_map>
#include <stdlib.h>
// Include header files of MPI
#include "mpi.h"

// DO NOT CHANGE THE VALUES OF THE MACROS BELOW:
#define CREATE_AND_FREE_SINGLE_COMM     (1)
#define CREATE_AND_FREE_MULTIPLE_COMM   (2)
#define ALLREDUCE                       (3)

#define MEASURE_MODE  (3)


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

    std::unordered_map<int,double> etime;

#if MEASURE_MODE == CREATE_AND_FREE_SINGLE_COMM
    const int n_proc_start = 2;
    const int n_proc_end   = n_proc/2;
    for (int n_proc_in_group = n_proc_start; n_proc_in_group <= n_proc_end; n_proc_in_group *= 2) { 
        // Get the parent group
        MPI_Group parent_group;
        MPI_Comm_group(MPI_COMM_WORLD,&parent_group);
        // Create a group of n_proc_in_group processes
        std::vector<int> ranks;
        ranks.reserve(n_proc_in_group);
        for (int i=0; i<n_proc_in_group; i++) ranks[i]=i;
        MPI_Group group;
        MPI_Group_incl(parent_group, n_proc_in_group, &ranks[0], &group);
        // Measure a time to create and destroy a new comm
        const int n_trials = 32;
        MPI_Barrier(MPI_COMM_WORLD);
        double time_offset = get_time_of_day();
        for (int i = 0; i<n_trials; i++) {
            // Create a communicator corresponding to this group
            MPI_Comm comm;
            MPI_Comm_create(MPI_COMM_WORLD, group, &comm);
            // Free the communicator
            if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        etime[n_proc_in_group] = (get_time_of_day() - time_offset)/n_trials;
        // Free the group
        MPI_Group_free(&group);
    }
#elif MEASURE_MODE == CREATE_AND_FREE_MULTIPLE_COMM
    const int n_proc_start = 2;
    const int n_proc_end   = n_proc/2;
    for (int n_proc_in_group = n_proc_start; n_proc_in_group <= n_proc_end; n_proc_in_group *= 2) { 
        // Get the parent group
        MPI_Group parent_group;
        MPI_Comm_group(MPI_COMM_WORLD,&parent_group);
        // Examine my group id
        const int n_group = n_proc/n_proc_in_group;
        assert(n_proc % n_proc_in_group == 0);
        const int group_id = my_rank/n_proc_in_group;
        // Create groups of n_proc_in_group processes
        std::vector<int> ranks;
        ranks.reserve(n_proc_in_group);
        const int offset = n_proc_in_group * group_id;
        for (int i=0; i<n_proc_in_group; i++) ranks[i]=offset+i;
        MPI_Group group;
        MPI_Group_incl(parent_group, n_proc_in_group, &ranks[0], &group);
        // Measure a time to create and destroy a new comm
        const int n_trials = 32;
        MPI_Barrier(MPI_COMM_WORLD);
        double time_offset = get_time_of_day();
        for (int i = 0; i<n_trials; i++) {
            // Create a communicator corresponding to this group
            MPI_Comm comm;
            MPI_Comm_create(MPI_COMM_WORLD, group, &comm);
            // Free the communicator
            if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        etime[n_proc_in_group] = (get_time_of_day() - time_offset)/n_trials;
        // Free the group
        MPI_Group_free(&group);
    }
#elif MEASURE_MODE == ALLREDUCE
    const int n_proc_start = 2;
    const int n_proc_end   = n_proc/2;
    for (int n_proc_in_group = n_proc_start; n_proc_in_group <= n_proc_end; n_proc_in_group *= 2) { 
        // Get the parent group
        MPI_Group parent_group;
        MPI_Comm_group(MPI_COMM_WORLD,&parent_group);
        // Examine my group id
        const int n_group = n_proc/n_proc_in_group;
        assert(n_proc % n_proc_in_group == 0);
        const int group_id = my_rank/n_proc_in_group;
        // Create groups of n_proc_in_group processes
        std::vector<int> ranks;
        ranks.reserve(n_proc_in_group);
        const int offset = n_proc_in_group * group_id;
        for (int i=0; i<n_proc_in_group; i++) ranks[i]=offset+i;
        MPI_Group group;
        MPI_Group_incl(parent_group, n_proc_in_group, &ranks[0], &group);
        // Create a communicator corresponding to this group
        MPI_Comm comm;
        MPI_Comm_create(MPI_COMM_WORLD, group, &comm);
        // Measure a time to create and destroy a new comm
        constexpr int n_trials = 1000;
        double sendbuf = my_rank;
        double recvbuf {0.0};
        MPI_Barrier(MPI_COMM_WORLD);
        double time_offset = get_time_of_day();
        for (int i = 0; i<n_trials; i++) {
            MPI_Allreduce(&sendbuf, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, comm);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        etime[n_proc_in_group] = (get_time_of_day() - time_offset)/n_trials;
        // Free the communicator
        if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
        // Free the group
        MPI_Group_free(&group);
    }
#else
#error The value of the macro `MEASURE_MODE` is invalid.
#endif

    // Output the measured data
    if (my_rank == 0) {
        std::cout << "# n_proc = " << n_proc << std::endl;
        std::cout << "# n_proc_in_group,    etime [s]" << std::endl;
        for (auto itr = etime.begin(); itr != etime.end(); ++itr) 
            std::cout << itr->first << "    " << itr->second << std::endl;
    }

    // Finalize
    MPI_Finalize();
    return 0;
}
