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
#include <algorithm>
// Include header files of MPI
#include "mpi.h"

#define DEBUG

double get_time_of_day(){
    struct timeval t;
    gettimeofday(&t,NULL);
    return (double)(t.tv_sec + t.tv_usec*1.0e-6);
}


class MMRedistributor {
public:
    int n_group;
    std::vector<int> rank_start, rank_end;
    std::vector<int> lm_start, lm_end;
    int n_proc_min_in_group;
    int n_proc_in_my_group;
    int n_proc_for_intragroup_comm;
    int idx_to_my_group;
    int rank_in_my_group;
    int n_lm_treated_by_me;
    MPI_Group group_alltoallv;
    MPI_Group group_gatherv;
    MPI_Comm comm_alltoallv;
    MPI_Comm comm_gatherv;

    MMRedistributor() {};

    void freeCommunicator() {
        if (group_alltoallv != MPI_GROUP_NULL) MPI_Group_free(&group_alltoallv);
        if (group_gatherv   != MPI_GROUP_NULL) MPI_Group_free(&group_gatherv);
        if (comm_alltoallv != MPI_COMM_NULL) MPI_Comm_free(&comm_alltoallv);
        if (comm_gatherv   != MPI_COMM_NULL) MPI_Comm_free(&comm_gatherv);
    }

    void makeCommunicator(MPI_Comm parent_comm,
                          const int n_mm_compo,
                          const bool debug_flag = false) {
        // Get information of parent_comm
        int n_proc, my_rank;
        MPI_Group parent_group;
        MPI_Comm_size(parent_comm, &n_proc);
        MPI_Comm_rank(parent_comm, &my_rank);
        MPI_Comm_group(parent_comm, &parent_group);
        // Divide parent_comm into groups
        n_group = n_mm_compo;
        n_proc_min_in_group = n_proc/n_group;
        int n_rem = n_proc - n_proc_min_in_group * n_group;
        rank_start.resize(n_group);
        rank_end.resize(n_group);
        rank_start[0] = 0;
        for (int i=0; i<n_group; i++) {
            int n_proc_in_group = n_proc_min_in_group;
            if (i < n_rem) n_proc_in_group++;
            rank_end[i] = rank_start[i] + (n_proc_in_group - 1);
            if (i < n_group-1) rank_start[i+1] = rank_end[i] + 1;
        }
        if (debug_flag) {
            if (my_rank == 0) {
                std::cout << "### Information of groups for M2L ###" << std::endl;
                std::cout << "(n_proc = " << n_proc << ")" << std::endl;
                for (int i=0; i<n_group; i++) {
                    std::cout << "i = " << i
                              << " rank_start = " << rank_start[i]
                              << " rank_end = " << rank_end[i]
                              << std::endl;
                }
            }
        }
        // Calculate idx_to_my_group, n_proc_in_my_group, rank_in_my_group
        for (int i=0; i<n_group; i++)
            if ((rank_start[i] <= my_rank) && (my_rank <= rank_end[i])) {
                idx_to_my_group = i;
                n_proc_in_my_group = rank_end[i] - rank_start[i] + 1;
                rank_in_my_group = my_rank - rank_start[i];
            }
        // Make group_alltoallv & comm_alltoallv 
        std::vector<int> ranks;
        int rnk_start = rank_start[idx_to_my_group];
        int rnk_end = rank_end[idx_to_my_group];
        for (int rnk=rnk_start; rnk<=rnk_end; rnk++) ranks.push_back(rnk);
        MPI_Group_incl(parent_group, n_proc_in_my_group, &ranks[0], &group_alltoallv);
        MPI_Comm_create(MPI_COMM_WORLD, group_alltoallv, &comm_alltoallv);
        // Calculate n_proc_for_intragroup_comm
        if (n_proc_min_in_group <= n_group) n_proc_for_intragroup_comm = n_proc_min_in_group;
        else n_proc_for_intragroup_comm = n_group;
        // Make group_gatherv, comm_gatherv
        if (rank_in_my_group < n_proc_for_intragroup_comm) {
            ranks.resize(n_group);
            for (int i=0; i<n_group; i++) ranks[i] = rank_start[i] + rank_in_my_group;
            MPI_Group_incl(parent_group, n_group, &ranks[0], &group_gatherv);
            MPI_Comm_create(MPI_COMM_WORLD, group_gatherv, &comm_gatherv);
        } else {
            group_gatherv = MPI_GROUP_NULL;
            MPI_Comm_create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, &comm_gatherv);
        }
        // Calculate lm_start, lm_end
        const int n_lm_min = n_mm_compo/n_proc_for_intragroup_comm;
        n_rem = n_mm_compo - n_lm_min * n_proc_for_intragroup_comm;
        lm_start.resize(n_proc_for_intragroup_comm);
        lm_end.resize(n_proc_for_intragroup_comm);
        lm_start[0] = 0; 
        for (int i=0; i<n_proc_for_intragroup_comm; i++) {
            int n_lm = n_lm_min;
            if (i < n_rem) n_lm++;
            //if (my_rank == 0) std::cout << "i = " << i << " n_lm = " << n_lm << std::endl;
            lm_end[i] = lm_start[i] + (n_lm - 1);
            if (i < n_proc_for_intragroup_comm-1) lm_start[i+1] = lm_end[i] + 1;
        }
        if (rank_in_my_group < n_proc_for_intragroup_comm) {
            n_lm_treated_by_me = lm_end[rank_in_my_group] - lm_start[rank_in_my_group] + 1;
        } else {
            n_lm_treated_by_me = 0;
        }
        //MPI_Finalize();
        //std::exit(0);
        if (debug_flag) {
            if (my_rank == 0) {
                std::cout << "### Information of division of (l,m) ###" << std::endl;
                std::cout << "(n_mm_compo = " << n_mm_compo << ")" << std::endl;
                for (int i=0; i<n_proc_for_intragroup_comm; i++) {
                    std::cout << "i = " << i
                              << " lm_start = " << lm_start[i]
                              << " lm_end = " << lm_end[i]
                              << std::endl;
                }
            }
        }
    }

    void performAlltoallv(const int n_cell_loc,
                          const int * cell_idx_loc, 
                          const double * mm_loc,
                          int & n_cell_grp,
                          int *& cell_idx_grp,
                          double *& mm_grp) {
        //#### MPI comm. for cell_idx ####
        // Set recvcounts
        int *recvcounts;
        recvcounts = new int[n_proc_in_my_group];
        int itmp = n_cell_loc;
        MPI_Allgather(&itmp, 1, MPI_INT, 
                      recvcounts, 1, MPI_INT, comm_alltoallv);
        // Set recvdispls
        int recvcount_tot;
        int *recvdispls;
        recvdispls = new int[n_proc_in_my_group];
        recvdispls[0] = 0;
        recvcount_tot = recvcounts[0];
        for (int i=1; i<n_proc_in_my_group; i++) {
            recvdispls[i] = recvdispls[i-1] + recvcounts[i-1];
            recvcount_tot += recvcounts[i];
        }
        // Prepare cell_idx_grp
        n_cell_grp = recvcount_tot; 
        if (cell_idx_grp != nullptr) delete [] cell_idx_grp;
        if (n_cell_grp > 0) cell_idx_grp = new int[n_cell_grp];
        else cell_idx_grp = nullptr;
        // Make sendbuf
        int *sendbuf_i;
        if (n_cell_loc > 0) {
            sendbuf_i = new int[n_cell_loc];
            for (int i=0; i<n_cell_loc; i++) sendbuf_i[i] = cell_idx_loc[i];
        } else sendbuf_i = nullptr;
        // Perform MPI comm. 
        MPI_Allgatherv(sendbuf_i,  n_cell_loc, MPI_INT,
                       cell_idx_grp, recvcounts, recvdispls, MPI_INT,
                       comm_alltoallv);
        // Free
        if (sendbuf_i != nullptr) delete [] sendbuf_i;
        //#### MPI comm. for mm ####
        const int n_mm_compo = lm_end.back();
        // Set sendcounts
        int *sendcounts;
        sendcounts = new int[n_proc_in_my_group];
        int sendcount_tot = 0;
        for (int i=0; i<n_proc_in_my_group; i++) {
            if (i < n_proc_for_intragroup_comm) {
                int n_lm = lm_end[i] - lm_start[i] + 1;
                sendcounts[i] = n_lm * n_cell_loc;
            } else sendcounts[i] = 0;
            sendcount_tot += sendcounts[i];
        }
        // Set senddispls
        int *senddispls;
        senddispls = new int[n_proc_in_my_group];
        senddispls[0] = 0;
        for (int i=1; i<n_proc_in_my_group; i++) {
            senddispls[i] = senddispls[i-1] + sendcounts[i-1];
        }
        // Set recvcounts
        MPI_Alltoall(sendcounts,1,MPI_INT,recvcounts,1,MPI_INT,comm_alltoallv);
        // Set recvdispls
        recvcount_tot = recvcounts[0];
        recvdispls[0] = 0;
        for (int i=1; i<n_proc_in_my_group; i++) {
            recvdispls[i] = recvdispls[i-1] + recvcounts[i-1];
            recvcount_tot += recvcounts[i];
        }
        assert(recvcount_tot == n_lm_treated_by_me * n_cell_grp);
        // Prepare mm_grp (recvbuf)
        if (mm_grp != nullptr) delete [] mm_grp;
        if (recvcount_tot > 0) mm_grp = new double[recvcount_tot];
        else mm_grp = nullptr;
        // Make sendbuf and perform MPI_Alltoallv
        double *sendbuf;
        if (sendcount_tot > 0) sendbuf = new double[sendcount_tot];
        else sendbuf = nullptr;
        int idx = 0;
        for (int k=0; k<n_proc_for_intragroup_comm; k++) { // loop for destination
            for (int i=0; i < n_cell_loc; i++) {
                for (int lm=lm_start[k]; lm<=lm_end[k]; lm++) {
                    sendbuf[idx++] = mm_loc[n_mm_compo * i + lm];
                }
            }
        }
        MPI_Alltoallv(sendbuf, sendcounts, senddispls, MPI_DOUBLE,
                      mm_grp, recvcounts, recvdispls, MPI_DOUBLE,
                      comm_alltoallv);
        // Note that the layout of mm_grp is as follows.
        // | aa | bb | cc | dd | ...
        // where each block represents n_lm_treated_by_me MM,
        // and the number of blocks is n_cell_out.
        // 
        // Free
        delete [] sendcounts;
        delete [] senddispls;
        delete [] sendbuf;
        delete [] recvcounts;
        delete [] recvdispls;
    }

    void performGatherv(const int n_cell_grp,
                        const int * cell_idx_grp,
                        const double * mm_grp,
                        int & n_cell_tot,
                        int *& cell_idx_tot,
                        double *& mm_tot) {
        if (comm_gatherv != MPI_COMM_NULL) {
            for (int lm=lm_start[rank_in_my_group]; lm<=lm_end[rank_in_my_group]; lm++) {
                //#### MPI comm. for cell_idx ####
                // Set recvcounts
                int *recvcounts = new int[n_group];
                int itmp = n_cell_grp;
                MPI_Gather(&itmp, 1, MPI_INT, 
                           recvcounts, 1, MPI_INT,
                           lm,comm_gatherv);
                // Set recvdispls
                int recvcount_tot;
                int *recvdispls = new int[n_group];
                if (idx_to_my_group == lm) { // root process
                    recvdispls[0] = 0;
                    recvcount_tot = recvcounts[0];
                    for (int i=1; i<n_group; i++) {
                        recvdispls[i] = recvdispls[i-1] + recvcounts[i-1];
                        recvcount_tot += recvcounts[i];
                    }
                }
                // Prepare cell_idx_tot
                if (idx_to_my_group == lm) { // root process
                    n_cell_tot = recvcount_tot;
                    if (cell_idx_tot != nullptr) delete [] cell_idx_tot;
                    if (recvcount_tot > 0) cell_idx_tot = new int[recvcount_tot];
                    else cell_idx_tot = nullptr;
                } else {
                    n_cell_tot = 0;
                    cell_idx_tot = nullptr;
                }
                // Make sendbuf_i
                int *sendbuf_i;
                if (n_cell_grp > 0) {
                    sendbuf_i = new int[n_cell_grp];
                    for (int i=0; i<n_cell_grp; i++) sendbuf_i[i] = cell_idx_grp[i];
                } else sendbuf_i = nullptr;
                // Perform MPI comm.
                MPI_Gatherv(sendbuf_i, n_cell_grp, MPI_INT,
                            cell_idx_tot, recvcounts, recvdispls, MPI_INT,
                            lm, comm_gatherv);
                // Free
                if (sendbuf_i != nullptr) delete [] sendbuf_i;
                //#### MPI comm. for mm ####
                // Prepare mm_tot
                if (idx_to_my_group == lm) { // root process
                    if (mm_tot != nullptr) delete [] mm_tot;
                    if (recvcount_tot > 0) mm_tot = new double[recvcount_tot];
                    else mm_tot = nullptr;
                } else {
                    mm_tot = nullptr;
                }
                // Prepare sendbuf
                double *sendbuf;
                if (n_cell_grp > 0) {
                    sendbuf = new double[n_cell_grp];
                    int idx {0};
                    for (int i=0; i<n_cell_grp; i++) {
                        const int offset = lm - lm_start[rank_in_my_group];
                        sendbuf[idx++] = mm_grp[n_lm_treated_by_me * i + offset];
                    }
                } else sendbuf = nullptr;
                // Perform MPI comm.
                MPI_Gatherv(sendbuf, n_cell_grp, MPI_DOUBLE,
                            mm_tot, recvcounts, recvdispls, MPI_DOUBLE,
                            lm, comm_gatherv);
                // Free
                delete [] recvcounts;
                delete [] recvdispls;
                if (sendbuf != nullptr) delete [] sendbuf;
            }
        } else {
            n_cell_tot = 0;
            if (cell_idx_tot != nullptr) delete [] cell_idx_tot; 
            cell_idx_tot = nullptr;
            if (mm_tot != nullptr) delete [] mm_tot;
            mm_tot = nullptr;
        }
    }

};

int main(int argc, char *argv[]) {

    // Set the problem
    const int p = 7;
    const int n_mm_compo = (p+1)*(p+1);
    const int NC = 256;
    const int NC3 = NC*NC*NC;
   
    int n_proc, my_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    assert(n_proc > 2*n_mm_compo);
    // This code is valid only under this condition.

    // Make groups and communicators for gathering MM
    MMRedistributor redist;
    redist.makeCommunicator(MPI_COMM_WORLD, n_mm_compo);

    // Make input data
    int n_cell_loc = NC3/n_proc;
    const int diff = NC3 - n_cell_loc * n_proc;
    if (diff > 0) {
        assert(0 < diff && diff <= n_proc);
        if (my_rank < diff) n_cell_loc++;
    }
    //std::cout << "my_rank = " << my_rank << " n_cell_loc = " << n_cell_loc << std::endl;
    int *cidx_loc = new int[n_cell_loc];
    double *mm_loc = new double[n_mm_compo * n_cell_loc];
    memset(cidx_loc, 0, (size_t)(sizeof(int) * n_cell_loc));
    memset(mm_loc, 0, (size_t)(sizeof(double) * n_mm_compo * n_cell_loc));

    double time_offset = get_time_of_day();

    // MPI_Alltoallv in inter-group
    int n_cell_grp;
    int * cidx_grp {nullptr};
    double * mm_grp {nullptr};
    redist.performAlltoallv(n_cell_loc, cidx_loc, mm_loc,
                            n_cell_grp, cidx_grp, mm_grp);

    // Simple check
    //std::cout << "my_rank = " << my_rank << " n_cell_grp = " << n_cell_grp << std::endl;

    // MPI_Gatherv in intra-group
    int n_cell_tot;
    int * cidx_tot {nullptr};
    double * mm_tot {nullptr};
    redist.performGatherv(n_cell_grp, cidx_grp, mm_grp,
                          n_cell_tot, cidx_tot, mm_tot);

    double etime = get_time_of_day() - time_offset;

    // Simple check
    std::cout << "my_rank = " << my_rank
              << " n_cell_tot = " << n_cell_tot
              << " etime = " << etime
              << std::endl;


    // Finalize
    redist.freeCommunicator();
    MPI_Finalize();

    return 0;
}
