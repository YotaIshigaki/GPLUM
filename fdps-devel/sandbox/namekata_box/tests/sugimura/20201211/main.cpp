// Include FDPS header
#include <particle_simulator.hpp>
// Include the standard C++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cassert>
#include <vector>
// Include the user-defined headers
#include "user-defined.hpp"

int main(int argc, char* argv[]) {
    PS::Initialize(argc, argv);
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::S32 my_rank = PS::Comm::getRank();
    constexpr PS::S32 root = 0;
    constexpr PS::S32 checking_rank = 1;

    // Determine # of local particles using a pseudo-random number generator
    const PS::U64 seed = static_cast<PS::U64>(my_rank);
    std::mt19937_64 mt(seed);
    constexpr PS::S32 N_LOC_MAX = 128; // the upper limit of # of local particles (parameter)
    assert(N_LOC_MAX > 0);
    std::uniform_int_distribution<PS::S32> dist(1, N_LOC_MAX);
    PS::S32 n_loc = dist(mt); // set # of local particles
    std::vector<PS::S32> n_loc_list(n_proc);
    PS::Comm::gather(&n_loc, 1, n_loc_list.data(), root);
    PS::Comm::broadcast(n_loc_list.data(), n_proc, root);
    PS::S32 n_glb {0};
    std::vector<PS::S32> displs(n_proc+1);
    displs[0] = 0;
    for (PS::S32 i=0; i<n_proc; i++) {
        n_glb += n_loc_list[i];
        displs[i+1] = displs[i] + n_loc_list[i];
    }
    if (my_rank == checking_rank) {
        std::cout << "n_glb = " << n_glb << std::endl;
        for (PS::S32 i=0; i<n_proc; i++) {
            std::cout << "rank = " << i
                      << ", n_loc = " << n_loc_list[i]
                      << ", displ = " << displs[i]
                      << std::endl; 
        }
    }

    // Compute the matrix elements
    //
    // Instead of the following codes, tree.calcForceAllAndWriteBack()
    // must be performed here.
    //
    constexpr PS::S32 N_NGB_GEN_LB = N_NGB_TYPICAL-5;
    constexpr PS::S32 N_NGB_GEN_UB = N_NGB_TYPICAL+5;
    assert((N_NGB_GEN_LB < N_NGB_MAX) &&
           (N_NGB_GEN_UB < N_NGB_MAX) &&
           (N_NGB_GEN_LB <= N_NGB_GEN_UB));
    std::uniform_int_distribution<PS::S32> dist_ngb(N_NGB_GEN_LB, N_NGB_GEN_UB);
    std::uniform_int_distribution<PS::S32> dist_id(0, n_glb-1);
    std::uniform_real_distribution<PS::F64> dist_val(0, 1.0);
    std::vector<FullParticle> psys(n_loc);
    for (PS::S32 i=0; i<n_loc; i++) {
        const PS::S32 id_i = i + displs[my_rank];
        psys[i].id = id_i;
        const PS::S32 N_ngb_nonzero = dist_ngb(mt);
        for (PS::S32 j=0; j<N_NGB_MAX; j++) {
            if (j < N_ngb_nonzero) {
                PS::S32 id_j;
                do {
                    id_j = dist_id(mt);
                    // Check if id_j overlaps with the existing IDs or not.
                    bool is_overlapped {false};
                    for (PS::S32 jj=0; jj<j; jj++) {
                        if (id_j == psys[i].force.id[jj]) is_overlapped = true;
                    }
                    if (!is_overlapped) break;
                } while(1);
                psys[i].force.id[j] = id_j;
                psys[i].force.val[j] = dist_val(mt);
            } else {
                psys[i].force.id[j] = -1;
                psys[i].force.val[j] = 0;
            }
        }
    }

    // Gather information using MPI comm.
    // [1] Make a send buffer
    std::vector<CommBuffer> sendbuf;
    for (PS::S32 i=0; i<n_loc; i++) {
        for (PS::S32 j=0; j<N_NGB_MAX; j++) {
            if (psys[i].force.val[j] != 0.0) {
                CommBuffer tmp;
                tmp.i = psys[i].id;
                tmp.j = psys[i].force.id[j]; 
                tmp.val = psys[i].force.val[j];
                sendbuf.push_back(tmp);
            }
        }
    }
    // [2] Gather # of buffer elements to RANK 0
    PS::S32 sendcount = sendbuf.size();
    std::vector<PS::S32> recvcounts(n_proc);
    PS::Comm::gather(&sendcount, 1, recvcounts.data(), root);
    // [3] Calculate # of data to be received at RANK 0 and the displacements
    PS::S32 recvcount_tot {0};
    std::vector<PS::S32> recvdispls(n_proc+1);
    if (my_rank == root) {
        recvdispls[0] = 0;
        for (PS::S32 i=0; i<n_proc; i++) {
            recvcount_tot += recvcounts[i];
            recvdispls[i+1] = recvdispls[i] + recvcounts[i];
        }
    }
    PS::Comm::broadcast(&recvcount_tot, 1, root);
    PS::Comm::broadcast(recvcounts.data(), n_proc, root);
    PS::Comm::broadcast(recvdispls.data(), n_proc+1, root);
    if (my_rank == checking_rank) {
        std::cout << "recvcount_tot = " << recvcount_tot << std::endl;
        for (PS::S32 i=0; i<n_proc; i++) {
            std::cout << "rank = " << i
                      << ", recvcount = " << recvcounts[i]
                      << ", recvdispl = " << recvdispls[i]
                      << std::endl;
        }
    }
    // [4] Prepare a receive buffer at RANK 0
    std::vector<CommBuffer> recvbuf(recvcount_tot);
    // [5] Gather the information to RANK 0
    PS::Comm::gatherV(sendbuf.data(), sendcount,
                      recvbuf.data(), recvcounts.data(), recvdispls.data(), root);

    // Output to a file to check
    if (my_rank == root) {
        const std::string file_name = "result.txt";
        std::ofstream ofs;
        ofs.open(file_name.c_str(), std::ios::trunc);
        for (PS::S32 i=0; i<recvcount_tot; i++) {
            ofs << recvbuf[i].i << "   "
                << recvbuf[i].j << "   "
                << recvbuf[i].val
                << std::endl;
        }
        ofs.close();
    }
 
    PS::Finalize();
    return 0;
}
