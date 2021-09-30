// Include the standard C++ headers
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <vector>
#include <unordered_map>
#include <stdlib.h>
#include <algorithm>
// Include header files of MPI
#include "mpi.h"
// Include header files of FFTW
#include <fftw3-mpi.h> 

//#define OUTPUT_DATASIZE

//------------------------------------------------------
// DO NOT CHANGE THE VALUES OF THE FOLLOWING MACROS:
#define IRECV_FIRST_THEN_ISEND                  (0)
#define IRECV_FIRST_THEN_SEND                   (1)
#define IRECV_FIRST_THEN_SEND_RINGLIKE          (2)
#define IRECV_RINGLIKE_FIRST_THEN_SEND_RINGLIKE (3)
#define IRECV_ISEND_MIXED_RINGLIKE              (4)
#define SENDRECV_ONLY_RINGLIKE                  (5)
#define SENDRECV_SEND_RECV_MIXED_RINGLIKE       (6)
#define ALLTOALLV                               (7)
#define IRECV_FIRST_THEN_SEND_SEQUENCER_2D      (8)

#define RINGLIKE                                (0)
#define NESTED_SQUARE                           (1)
#define NESTED_BY_MANHATTAN_DISTANCE            (2)
//------------------------------------------------------

#define COMM_IMPL_TYPE (8)
#define SEQUENCER_TYPE (1)
// the value must be the same as that of
// one of the macros shown above.

double get_time_of_day(){
    struct timeval t;
    gettimeofday(&t,NULL);
    return (double)(t.tv_sec + t.tv_usec*1.0e-6);
}

enum class SequencerType {Ringlike, NestedSquare, NestedByManhattanDistance};

class Sequencer {
private:
    bool is_initialized_;
    int n_proc_x_, n_proc_y_, n_proc_;
    SequencerType type_;
    
    void calcPosition(const int rank, int & x, int & y) {
        assert(is_initialized_);
        y = rank/n_proc_x_;
        x = rank - n_proc_x_ * y;
    }

public:

    Sequencer(): is_initialized_(false) {}

    void initialize(const int n_proc_x,
                    const int n_proc_y,
                    const SequencerType type) {
        assert(n_proc_x > 0 && n_proc_y > 0);
        n_proc_x_ = n_proc_x;
        n_proc_y_ = n_proc_y;
        n_proc_   = n_proc_x_ * n_proc_y_;
        type_ = type;
        is_initialized_ = true;
    }


    void calcRelativePositionOfCommPartner(const int kth, const int my_rank,
                                           int & x, int & y, int & vx, int & vy,
                                           const bool debug = false) {
        // Calculate the position of my_rank
        calcPosition(my_rank, x, y);

        if (type_ == SequencerType::Ringlike) {
            // Calculate the destination rank
            const int dest = (my_rank + kth) % n_proc_;
            int x_dest, y_dest;
            calcPosition(dest, x_dest, y_dest);
            // Calculate the relative position of destination rank (minimum image)
            int vx = x_dest - x;
            int vy = y_dest - y;
            if (vx >  n_proc_x_/2) vx -= n_proc_x_;
            if (vx < -n_proc_x_/2) vx += n_proc_x_;
            if (vy >  n_proc_y_/2) vy -= n_proc_y_;
            if (vy < -n_proc_y_/2) vy += n_proc_y_;

        } else if (type_ == SequencerType::NestedSquare) {
            // Calculate the current level
            int lev {0};
            int k_start {0}, k_end {0};
            for(;;) {
                if (k_start <= kth && kth <= k_end) break;
                lev++;
                k_start = k_end + 1;
                k_end = (2*lev+1) * (2*lev+1) - 1;
            }
            //std::cout << "lev = " << lev << " k_start = " << k_start << " k_end = " << k_end << std::endl;
            // Set the default start point (correspond to k_start)
            vx = lev;
            vy = 0;
            // Set the initial shift of the starting point
            int n_shift_offset {0};
#if 1
#if 0
            // Pattern #1
            n_shift_offset = (my_rank % 4) * lev;
#else
            // Pattern #2
            {
                const int rem_x = x % 2;
                const int rem_y = y % 2;
                if (rem_x == 1 && rem_y == 0) n_shift_offset = 2*lev;
                else if (rem_x == 1 && rem_y == 1) n_shift_offset = 2*(2*lev);
                else if (rem_x == 0 && rem_y == 0) n_shift_offset = 3*(2*lev);
            }
#endif
#endif
            // Calculate the relative position of destination
            int n_shift = kth - k_start + n_shift_offset;
            for(;;) {
                if (vx == lev && vy < lev) {
                    while (n_shift) {
                        if (vy < lev) {
                            vy++; n_shift--;
                        }
                        if (vy == lev) break;
                    }
                } else if (-lev < vx && vy == lev) {
                    while (n_shift) {
                        if (vx > -lev) {
                            vx--; n_shift--;
                        }
                        if (vx == -lev) break;
                    }
                } else if (vx == -lev && -lev < vy) {
                    while (n_shift) {
                        if (vy > -lev) {
                            vy--; n_shift--;
                        }
                        if (vy == -lev) break;
                    }
                } else if (vx < lev && vy == -lev) {
                    while (n_shift) {
                        if (vx < lev) {
                            vx++; n_shift--;
                        }
                        if (vx == lev) break;
                    }
                }
                if (n_shift == 0) break;
            }
        } else if (type_ == SequencerType::NestedByManhattanDistance) {
            // Calculate the current distance
            // Note that the number of cells at Manhattan distance d is
            // given by 4*d for d >= 1 in two dimension. For d = 0,
            // the number of cells is 1. 
            int dis {0};
            int k_start {0}, k_end {0};
            for (;;) {
                if (k_start <= kth && kth <= k_end) break;
                dis++;
                k_start = k_end + 1;
                k_end = k_start + 4*dis - 1;
            }
            if (debug) std::cout << "dis = " << dis << " k_start = " << k_start << " k_end = " << k_end << std::endl;
            // Set the default starting posint (correspond to k_start)
            vx = dis;
            vy = 0;
            if (debug) std::cout << "vx = " << vx << " vy = " << vy << std::endl;
            // Set the initial shift of the starting point
            int n_shift_offset {0};
#if 0
#if 0
            // Pattern #1
            n_shift_offset = (my_rank % 4) * dis;
#else
            // Pattern #2
            {
                const int rem_x = x % 2;
                const int rem_y = y % 2;
                if (rem_x == 1 && rem_y == 0) n_shift_offset = dis;
                else if (rem_x == 1 && rem_y == 1) n_shift_offset = 2*dis;
                else if (rem_x == 0 && rem_y == 1) n_shift_offset = 3*dis;
            }
#endif
#endif
            // Calculate the relative position of destination rank
            int n_shift = kth - k_start + n_shift_offset;
            if (debug) std::cout << "n_shift = " << n_shift << std::endl;
            for (;;) {
                if (debug) std::cout << "vx = " << vx << " vy = " << vy << std::endl;
                // Determine the direction of shift
                if ((1 <= vx && vx <= dis) && vy >= 0) {
                    while (n_shift) {
                        if (vy < dis) {
                            vx--; vy++; n_shift--;
                        }
                        if (vy == dis) break;
                    }
                } else if (vx <= 0 && (1 <= vy && vy <= dis)) {
                    while (n_shift) {
                        if (vx > -dis) {
                            vx--; vy--; n_shift--;
                        }
                        if (vx == -dis) break;
                    }
                } else if ((-dis <= vx && vx <= -1) && vy <= 0) {
                    while (n_shift) {
                        if (vy > -dis) {
                            vx++; vy--; n_shift--;
                        }
                        if (vy == -dis) break;
                    }
                } else if (vx >= 0 && (-dis <= vy && vy <= -1)) {
                    while (n_shift) {
                        if (vx < dis) {
                            vx++; vy++; n_shift--;
                        }
                        if (vx == dis) break;
                    }
                }
                if (n_shift == 0) break;
            }
        }
    }

};

class MMData {
public:
    int cidx;
    double mm;
};

class MMRedistributor {
public:
    int n_group;
    std::vector<int> rank_start, rank_end;
    std::vector<ptrdiff_t> local_0_start, local_0_end; // z
    int n_proc_min_in_group;
    int n_proc_in_my_group;
    int n_proc_in_parent_group;
    int n_proc_for_fft;
    int idx_to_my_group;
    int rank_in_my_group;
    int rank_in_parent_group;
    MPI_Group group_all;
    MPI_Group group_fft;
    MPI_Comm parent_comm;
    MPI_Comm comm_all;
    MPI_Comm comm_fft;

    MMRedistributor() {};

    void freeCommunicator() {
        if (group_all != MPI_GROUP_NULL) MPI_Group_free(&group_all);
        if (group_fft != MPI_GROUP_NULL) MPI_Group_free(&group_fft);
        if (comm_all != MPI_COMM_NULL) MPI_Comm_free(&comm_all);
        if (comm_fft != MPI_COMM_NULL) MPI_Comm_free(&comm_fft);
    }

    void makeCommunicator(MPI_Comm _parent_comm,
                          const int n_mm_compo,
                          const int NC,
                          const bool debug_flag = false) {
        // Set parent_comm
        parent_comm = _parent_comm;
        // Get information of parent_comm
        int n_proc, my_rank;
        MPI_Group parent_group;
        MPI_Comm_size(parent_comm, &n_proc);
        MPI_Comm_rank(parent_comm, &my_rank);
        MPI_Comm_group(parent_comm, &parent_group);
        n_proc_in_parent_group = n_proc;
        rank_in_parent_group = my_rank;
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
        // Make group_all & comm_all
        std::vector<int> ranks;
        int rnk_start = rank_start[idx_to_my_group];
        int rnk_end = rank_end[idx_to_my_group];
        for (int rnk=rnk_start; rnk<=rnk_end; rnk++) ranks.push_back(rnk);
        MPI_Group_incl(parent_group, n_proc_in_my_group, &ranks[0], &group_all);
        MPI_Comm_create(MPI_COMM_WORLD, group_all, &comm_all);
        // Calculate n_proc_for_fft
        if (n_proc_min_in_group <= NC/2) n_proc_for_fft = n_proc_min_in_group;
        else n_proc_for_fft = NC/2;
        // Make group_fft, comm_fft
        if (rank_in_my_group < n_proc_for_fft) {
            ranks.resize(n_proc_for_fft);
            for (int i=0; i<n_proc_for_fft; i++)
                ranks[i] = rank_start[idx_to_my_group] + i;
            MPI_Group_incl(parent_group, n_proc_for_fft, &ranks[0], &group_fft);
            MPI_Comm_create(MPI_COMM_WORLD, group_fft, &comm_fft);
        } else {
            group_fft = MPI_GROUP_NULL;
            MPI_Comm_create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, &comm_fft);
        }
        //MPI_Finalize();
        //std::exit(0);
    }

    void makeDecompInfo(const int NX, 
                        const int NY,
                        const int NZ) {
        // Get MPI_Datype for ptrdiff_t
        bool is_type_commited {false};
        MPI_Datatype type;
        MPI_Type_match_size(MPI_TYPECLASS_INTEGER, sizeof(ptrdiff_t), &type);
        if (type == MPI_DATATYPE_NULL) {
            MPI_Type_contiguous(sizeof(ptrdiff_t), MPI_BYTE, &type);
            MPI_Type_commit(&type);
            is_type_commited = true;
        }
        // Resize local_0_start, local_z_end
        local_0_start.resize(n_proc_for_fft);
        local_0_end.resize(n_proc_for_fft);
        // Calculate z_start, z_end
        if (comm_fft != MPI_COMM_NULL) {
            ptrdiff_t alloc_local, local_n0, _local_0_start;
            ptrdiff_t local_n1, _local_1_start;
            alloc_local = fftw_mpi_local_size_3d_transposed(NZ, NY, NX/2+1, comm_fft,
                                                            &local_n0, &_local_0_start,
                                                            &local_n1, &_local_1_start);
            ptrdiff_t sendbuf[2];
            sendbuf[0] = _local_0_start;
            sendbuf[1] = local_n0;
            ptrdiff_t *recvbuf = new ptrdiff_t[2 * n_proc_for_fft];
            MPI_Gather(sendbuf, 2, type, recvbuf, 2, type, 0, comm_fft);
            if (rank_in_my_group == 0) {
                for (int i=0; i<n_proc_for_fft; i++) {
                    local_0_start[i] = recvbuf[2 * i];
                    local_0_end[i]   = local_0_start[i] + recvbuf[2 * i + 1] - 1;
                }
            }
            delete [] recvbuf;
        }
        MPI_Bcast(&local_0_start[0], n_proc_for_fft, type, 0, comm_all);
        MPI_Bcast(&local_0_end[0], n_proc_for_fft, type, 0, comm_all);
        // Free
        if (is_type_commited) MPI_Type_free(&type);
        // Check
        if (rank_in_parent_group == 0) {
            std::cout << "### Information of decomposition of FFTW ###" << std::endl;
            std::cout << "(NZ = " << NZ << ")" << std::endl;
            for (int i=0; i<n_proc_for_fft; i++) {
                std::cout << "i = " << i
                          << " local_0_start = " << local_0_start[i]
                          << " local_0_end = " << local_0_end[i]
                          << std::endl;
            }
        }
        
    }

    void performComm(const int n_cell_loc,
                     const int * cell_idx_loc, 
                     const double * mm_loc,
                     const int n_mm_compo,
                     const int NX,
                     const int NY,
                     const int NZ,
                     const int * cidx2rnk,
                     int & n_cell_tot,
                     int *& cell_idx_tot,
                     double *& mm_tot,
                     double & etime_comm_only) {
        // Make a list that describes what number of data is sent to each rank
        int *nsend = new int[n_proc_in_parent_group];
        memset(nsend, 0, (size_t)(sizeof(int) * n_proc_in_parent_group));
        for (int n=0; n<n_cell_loc; n++) {
            const int cidx = cell_idx_loc[n];
            const int iz = cidx/(NX*NY);
            int slab_id;
            for (int k=0; k<n_proc_for_fft; k++) {
                if ((local_0_start[k] <= iz) && (iz <= local_0_end[k])) {
                    slab_id = k;
                    break;
                }
            }
            for (int k=0; k<n_group; k++) {
                const int rnk = rank_start[k] + slab_id;
                nsend[rnk]++;
            }
        }
        // Make a list that describes what number of data is sent from other processes
        int *nrecv = new int[n_proc_in_parent_group];
        memset(nrecv, 0, (size_t)(sizeof(int) * n_proc_in_parent_group));
        if (rank_in_my_group < n_proc_for_fft) { // this process receive data
            const int iz_start = local_0_start[rank_in_my_group];
            const int iz_end   = local_0_end[rank_in_my_group];
            for (int iz=iz_start; iz<=iz_end; iz++)
                for (int iy=0; iy<NY; iy++)
                    for (int ix=0; ix<NX; ix++) {
                        const int cidx = ix + NX*(iy + NY*iz);
                        const int rnk = cidx2rnk[cidx];
                        nrecv[rnk]++;
                    }
        }
        // Prepare sendbuf
        int n_dest {0};
        int sendcount_tot {0};
        for (int i=0; i<n_proc_in_parent_group; i++) {
            if (nsend[i] > 0) n_dest++;
            sendcount_tot += nsend[i];
        }
        MMData *sendbuf;
        if (sendcount_tot > 0) {
            sendbuf = new MMData[sendcount_tot];
            //for (int i=0; i<sendrecv_list.size(); i++) {
            //}
        } else sendbuf = nullptr;
        // Prepare recvbuf
        int n_source {0};
        int recvcount_tot {0};
        for (int i=0; i<n_proc_in_parent_group; i++) {
            if (nrecv[i] > 0) n_source++;
            recvcount_tot += nrecv[i];
        }
        MMData *recvbuf;
        if (recvcount_tot > 0) {
            recvbuf = new MMData[recvcount_tot];
        } else recvbuf = nullptr;
        // Get MPI_Datatype for MMData
        MPI_Datatype type;
        MPI_Type_contiguous(sizeof(MMData), MPI_BYTE, &type);
        MPI_Type_commit(&type);
        // Output 
#ifdef OUTPUT_DATASIZE
        std::ostringstream ss;
        ss << "./output/prof-" << std::setfill('0') << std::setw(5) << rank_in_parent_group << ".txt";
        std::string filename = ss.str();
        std::ofstream output_file;
        output_file.open(filename.c_str(), std::ios::trunc);
        for (int i = 0; i < n_proc_in_parent_group; i++) {
            output_file << i << "   " << nsend[i] << "   " << nrecv[i] << std::endl;
        }   
        output_file.close();
        MPI_Barrier(parent_comm);
#endif
#if COMM_IMPL_TYPE == IRECV_FIRST_THEN_SEND_SEQUENCER_2D
        // Make a sequence of ranks
        constexpr int n_proc_x = 32;
        constexpr int n_proc_y = 32;
        assert(n_proc_x * n_proc_y == n_proc_in_parent_group);
        Sequencer seq;
#if SEQUENCER_TYPE == NESTED_SQUARE
        seq.initialize(n_proc_x, n_proc_y, SequencerType::NestedSquare);
#elif SEQUENCER_TYPE == NESTED_BY_MAHATTAN_DISTANCE
        seq.initialize(n_proc_x, n_proc_y, SequencerType::NestedByManhattanDistance);
#else
#error The macro `SEQUENCER_TYPE` has an invalid value.
#endif
        std::vector<bool> is_registered;
        int n_registered {0};
        is_registered.resize(n_proc_in_parent_group);
        for (int i=0; i<n_proc_in_parent_group; i++) is_registered[i]=false;
        std::vector<int> dest_ranks;
        for (int i=0;; i++) {
            int x,y,vx,vy;
            seq.calcRelativePositionOfCommPartner(i, rank_in_parent_group, x, y, vx, vy);
            // Calculate the position of destination rank
            int x_dest = x + vx;
            while (x_dest < 0) { x_dest += n_proc_x; }
            while (x_dest >= n_proc_x) {x_dest -= n_proc_x; }
            int y_dest = y + vy;
            while (y_dest < 0) { y_dest += n_proc_y; }
            while (y_dest >= n_proc_y) { y_dest -= n_proc_y; }
            const int dest = x_dest + y_dest * n_proc_x;
            // Error check
            if (dest < 0 || n_proc_in_parent_group <= dest) {
                std::cout << "i = " << i
                          << " dest = "  << dest
                          << " x_dest = " << x_dest
                          << " y_dest = " << y_dest
                          << " x = " << x
                          << " y = " << y
                          << " vx = " << vx
                          << " vy = " << vy
                          << std::endl;
                assert(1);
            }
            // Registered if the rank indicated by (x,y,vx,vy) is first appearance
            if (is_registered[dest] == false) {
                dest_ranks.push_back(dest);
                is_registered[dest] = true;
                n_registered++;
            }
            if (n_registered == n_proc_in_parent_group) break;
        }
#endif
        // Perform MPI comm. 
        double time_offset = get_time_of_day();
#if COMM_IMPL_TYPE == IRECV_FIRST_THEN_ISEND
        //### Firstly perform MPI_Irecv and then peform MPI_Isend ###
        const int n_req = n_dest + n_source;
        MPI_Request *req;
        if (n_req > 0) req = new MPI_Request[n_req];
        else req = nullptr;
        int idx_req {0};
        // MPI_Irecv
        int recvdisp {0};
        for (int i=0; i<n_proc_in_parent_group; i++) {
            const int recvcount = nrecv[i];
            if (recvcount > 0) {
                const int recvtag = 0;
                MPI_Irecv(&recvbuf[recvdisp], recvcount, type, i, recvtag, parent_comm, &req[idx_req++]);
            }
            recvdisp += recvcount;
        }
        // MPI_Isend
        int senddisp {0};
        for (int i=0; i<n_proc_in_parent_group; i++) {
            const int sendcount = nsend[i];
            if (sendcount > 0) {
                const int sendtag = 0;
                MPI_Isend(&sendbuf[senddisp], sendcount, type, i, sendtag, parent_comm, &req[idx_req++]);
            }
            senddisp += sendcount;
        }
        assert(idx_req == n_req);
        if (n_req > 0) {
            MPI_Status *stat = new MPI_Status[n_req];
            MPI_Waitall(n_req, req, stat);
            delete [] stat;
            delete [] req;
        }
#elif COMM_IMPL_TYPE == IRECV_FIRST_THEN_SEND
        //### Firstly perform MPI_Irecv and then peform MPI_Send ###
        MPI_Request *req;
        if (n_source > 0) req = new MPI_Request[n_source];
        else req = nullptr;
        int idx_req {0};
        // MPI_Irecv
        int recvdisp {0};
        for (int i=0; i<n_proc_in_parent_group; i++) {
            const int recvcount = nrecv[i];
            if (recvcount > 0) {
                const int recvtag = 0;
                MPI_Irecv(&recvbuf[recvdisp], recvcount, type, i, recvtag, parent_comm, &req[idx_req++]);
            }
            recvdisp += recvcount;
        }
        // MPI_Send
        int senddisp {0};
        for (int i=0; i<n_proc_in_parent_group; i++) {
            const int sendcount = nsend[i];
            if (sendcount > 0) {
                const int sendtag = 0;
                MPI_Send(&sendbuf[senddisp], sendcount, type, i, sendtag, parent_comm);
            }
            senddisp += sendcount;
        }
        if (n_source > 0) {
            MPI_Status *stat = new MPI_Status[n_source];
            MPI_Waitall(n_source, req, stat);
            delete [] stat;
            delete [] req;
        }
#elif COMM_IMPL_TYPE == IRECV_FIRST_THEN_SEND_RINGLIKE
        //### Firstly perform MPI_Irecv and then peform MPI_Send ###
        MPI_Request *req;
        if (n_source > 0) req = new MPI_Request[n_source];
        else req = nullptr;
        int idx_req {0};
        // MPI_Irecv
        int recvdisp {0};
        for (int i=0; i<n_proc_in_parent_group; i++) {
            const int recvcount = nrecv[i];
            if (recvcount > 0) {
                const int recvtag = 0;
                MPI_Irecv(&recvbuf[recvdisp], recvcount, type, i, recvtag, parent_comm, &req[idx_req++]);
            }
            recvdisp += recvcount;
        }
        // MPI_Send
        int senddisp {0};
        for (int i=0; i<n_proc_in_parent_group; i++) {
            const int period = n_proc_in_parent_group;
            const int dest   =  (rank_in_parent_group + i) % period;
            const int sendcount = nsend[dest];
            if (sendcount > 0) {
                const int sendtag = 0;
                MPI_Send(&sendbuf[senddisp], sendcount, type, dest, sendtag, parent_comm);
            }
            senddisp += sendcount;
        }
        if (n_source > 0) {
            MPI_Status *stat = new MPI_Status[n_source];
            MPI_Waitall(n_source, req, stat);
            delete [] stat;
            delete [] req;
        }
#elif COMM_IMPL_TYPE == IRECV_FIRST_THEN_SEND_SEQUENCER_2D
        //### Firstly perform MPI_Irecv and then peform MPI_Send ###
        MPI_Request *req;
        if (n_source > 0) req = new MPI_Request[n_source];
        else req = nullptr;
        int idx_req {0};
        // MPI_Irecv
        int recvdisp {0};
        for (int i=0; i<n_proc_in_parent_group; i++) {
            const int recvcount = nrecv[i];
            if (recvcount > 0) {
                const int recvtag = 0;
                MPI_Irecv(&recvbuf[recvdisp], recvcount, type, i, recvtag, parent_comm, &req[idx_req++]);
            }
            recvdisp += recvcount;
        }
        // MPI_Send
        int senddisp {0};
        for (int i=0; i<n_proc_in_parent_group; i++) {
            const int dest   = dest_ranks[i];
            const int sendcount = nsend[dest];
            if (sendcount > 0) {
                const int sendtag = 0;
                MPI_Send(&sendbuf[senddisp], sendcount, type, dest, sendtag, parent_comm);
            }
            senddisp += sendcount;
        }
        if (n_source > 0) {
            MPI_Status *stat = new MPI_Status[n_source];
            MPI_Waitall(n_source, req, stat);
            delete [] stat;
            delete [] req;
        }
#elif COMM_IMPL_TYPE == IRECV_RINGLIKE_FIRST_THEN_SEND_RINGLIKE
        //### Firstly perform MPI_Irecv and then peform MPI_Send ###
        MPI_Request *req;
        if (n_source > 0) req = new MPI_Request[n_source];
        else req = nullptr;
        int idx_req {0};
        // MPI_Irecv
        int recvdisp {0};
        for (int i=0; i<n_proc_in_parent_group; i++) {
            const int period = n_proc_in_parent_group;
            const int source =  (rank_in_parent_group - i + period) % period;
            const int recvcount = nrecv[source];
            if (recvcount > 0) {
                const int recvtag = 0;
                MPI_Irecv(&recvbuf[recvdisp], recvcount, type, source, recvtag, parent_comm, &req[idx_req++]);
            }
            recvdisp += recvcount;
        }
        // MPI_Send irecv_ringlike_first_then_send_ringlike
        int senddisp {0};
        for (int i=0; i<n_proc_in_parent_group; i++) {
            const int period = n_proc_in_parent_group;
            const int dest   =  (rank_in_parent_group + i) % period;
            const int sendcount = nsend[dest];
            if (sendcount > 0) {
                const int sendtag = 0;
                MPI_Send(&sendbuf[senddisp], sendcount, type, dest, sendtag, parent_comm);
            }
            senddisp += sendcount;
        }
        if (n_source > 0) {
            MPI_Status *stat = new MPI_Status[n_source];
            MPI_Waitall(n_source, req, stat);
            delete [] stat;
            delete [] req;
        }
#elif COMM_IMPL_TYPE == IRECV_ISEND_MIXED_RINGLIKE
        //### perform MPI_Isend and MPI_Irecv togather in ringlike loop
        const int n_req = n_dest + n_source;
        MPI_Request *req;
        if (n_req > 0) req = new MPI_Request[n_req];
        else req = nullptr;
        int idx_req {0};
        int senddisp {0};
        int recvdisp {0};
        for (int i=0; i<n_proc_in_parent_group; i++) {
            const int period = n_proc_in_parent_group;
            const int dest   =  (rank_in_parent_group + i) % period;
            const int source =  (rank_in_parent_group - i + period) % period;
            const int sendcount = nsend[dest];
            const int recvcount = nrecv[source]; 
            const int sendtag = 0;
            const int recvtag = 0;
            if (sendcount > 0 && recvcount > 0) {
                MPI_Irecv(&recvbuf[recvdisp], recvcount, type, source, recvtag, parent_comm, &req[idx_req++]);
                MPI_Isend(&sendbuf[senddisp], sendcount, type, dest, sendtag, parent_comm, &req[idx_req++]);
            } else if (sendcount > 0) {
                MPI_Isend(&sendbuf[senddisp], sendcount, type, dest, sendtag, parent_comm, &req[idx_req++]);
            } else if (recvcount > 0) {
                MPI_Irecv(&recvbuf[recvdisp], recvcount, type, source, recvtag, parent_comm, &req[idx_req++]);
            }
            senddisp += sendcount;
            recvdisp += recvcount;
        }
        assert(idx_req == n_req);
        if (n_req > 0) {
            MPI_Status *stat = new MPI_Status[n_req];
            MPI_Waitall(n_req, req, stat);
            delete [] stat;
            delete [] req;
        }
#elif COMM_IMPL_TYPE == SENDRECV_ONLY_RINGLIKE
        // Using MPI_Sendrecv only
        int senddisp {0};
        int recvdisp {0};
        for (int i=0; i<n_proc_in_parent_group; i++) {
            const int period = n_proc_in_parent_group;
            const int dest   =  (rank_in_parent_group + i) % period;
            const int source =  (rank_in_parent_group - i + period) % period;
            const int sendcount = nsend[dest];
            const int recvcount = nrecv[source]; 
            const int sendtag = 0;
            const int recvtag = 0;
            MPI_Status stat;
            MPI_Sendrecv(&sendbuf[senddisp], sendcount, type, dest, sendtag,
                         &recvbuf[recvdisp], recvcount, type, source, recvtag,
                         parent_comm, &stat);
            senddisp += sendcount;
            recvdisp += recvcount;
        }
#elif COMM_IMPL_TYPE == SENDRECV_SEND_RECV_MIXED_RINGLIKE
        // Using MPI_Sendrecv, MPI_Send, MPI_Recv
        int senddisp {0};
        int recvdisp {0};
        for (int i=0; i<n_proc_in_parent_group; i++) {
            const int period = n_proc_in_parent_group;
            const int dest   =  (rank_in_parent_group + i) % period;
            const int source =  (rank_in_parent_group - i + period) % period;
            const int sendcount = nsend[dest];
            const int recvcount = nrecv[source]; 
            const int sendtag = 0;
            const int recvtag = 0;
            MPI_Status stat;
            if (sendcount > 0 && recvcount > 0) {
                MPI_Sendrecv(&sendbuf[senddisp], sendcount, type, dest, sendtag,
                             &recvbuf[recvdisp], recvcount, type, source, recvtag,
                             parent_comm, &stat);
            } else if (sendcount > 0) {
                MPI_Send(&sendbuf[senddisp], sendcount, type, dest, sendtag, parent_comm);
            } else if (recvcount > 0) {
                MPI_Recv(&recvbuf[recvdisp], recvcount, type, source, recvtag, parent_comm, &stat);
            }
            senddisp += sendcount;
            recvdisp += recvcount;
        }
//#elif COMM_IMPL_TYPE == ALLTOALLV
#else
#error the macro `COMM_IMPL_TYPE` is invalid.
#endif
        etime_comm_only = get_time_of_day() - time_offset;

        // Copy from recvbuf to cell_idx_tot & mm_tot

        // Free
        delete [] nsend;
        delete [] nrecv;
        if (sendbuf != nullptr) delete [] sendbuf;
        if (recvbuf != nullptr) delete [] recvbuf;
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
    assert(n_proc >= n_mm_compo);
    fftw_mpi_init();

    // Make groups and communicators for gathering MM
    MMRedistributor redist;
    redist.makeCommunicator(MPI_COMM_WORLD, n_mm_compo, NC);

    // Make information of decomposition of FFTW
    redist.makeDecompInfo(NC, NC, NC);

    // Set the number of local PM cells
    int n_cell_loc = NC3/n_proc;
    const int diff = NC3 - n_cell_loc * n_proc;
    if (diff > 0) {
        assert(0 < diff && diff <= n_proc);
        if (my_rank < diff) n_cell_loc++;
    }
    int *n_cell_tbl = new int[n_proc];
    MPI_Allgather(&n_cell_loc, 1, MPI_INT, n_cell_tbl, 1, MPI_INT, MPI_COMM_WORLD);

    // Set the mapping from cell index to rank number that is responsible for that cell.
    std::vector<int> cidx2rnk(NC3);
    int idx {0};
    for (int n=0; n<n_proc; n++)
        for (int i=0; i<n_cell_tbl[n]; i++)
            cidx2rnk[idx++] = n;

    // Set the input data (cidx_loc, mm_loc)
    int *cidx_loc = new int[n_cell_loc];
    double *mm_loc = new double[n_mm_compo * n_cell_loc];
    int disp {0};
    if (my_rank > 0) for (int n=0; n<my_rank; n++) disp += n_cell_tbl[n];
    idx = disp;
    for (int i=0; i<n_cell_loc; i++) {
        cidx_loc[i] = disp + i;
    }
    memset(mm_loc, 0, (size_t)(sizeof(double) * n_mm_compo * n_cell_loc));

    double time_offset = get_time_of_day();

    // MPI comm.
    int n_cell_tot {0};
    int * cidx_tot {nullptr};
    double * mm_tot {nullptr};
    double etime_comm_only;
    redist.performComm(n_cell_loc, cidx_loc, mm_loc,
                       n_mm_compo, NC, NC, NC, &cidx2rnk[0], 
                       n_cell_tot, cidx_tot, mm_tot,
                       etime_comm_only);

    double etime = get_time_of_day() - time_offset;

    // Simple check
    std::cout << "my_rank = " << my_rank
              << " n_cell_tot = " << n_cell_tot
              << " etime_comm_only = " << etime_comm_only 
              << " etime (total) = " << etime
              << std::endl;

    // Finalize
    redist.freeCommunicator();
    MPI_Finalize();

    return 0;
}
