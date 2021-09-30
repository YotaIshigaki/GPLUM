#pragma once
#include <complex>
#include <unordered_map>
#include "../ps_defs.hpp"
#include "particle_mesh_multipole_defs.hpp"

#define DEBUG_M2L_ENGINE_INITIALIZE
#define DEBUG_M2L_ENGINE_REDISTMM
#define DEBUG_M2L_ENGINE_CONVOLUTION
#define DEBUG_M2L_ENGINE_REDISTLE

namespace ParticleSimulator {
    namespace ParticleMeshMultipole {

        template <class T, int DIM = 4>
        class M2L_Buffer {
        public:
            S32 size_glb_tot_;
            S32 size_glb_[DIM];
            S32 start_glb_[DIM], end_glb_[DIM];
            S32 size_loc_tot_;
            S32 size_loc_[DIM];
            S32 start_loc_[DIM], end_loc_[DIM];
            S32 size_loc_used_tot_;
            S32 size_loc_used_[DIM];
            S32 start_loc_used_[DIM], end_loc_used_[DIM];
            std::vector<T> buf_; // local buffer

            M2L_Buffer() {
                size_glb_tot_ = 0;
                size_loc_tot_ = 0;
                size_loc_used_tot_ = 0;
                for (S32 i=0; i<DIM; i++) {
                    size_glb_[i] = 0;
                    start_glb_[i] = -1;
                    end_glb_[i] = -1;
                    size_loc_[i] = 0;
                    start_loc_[i] = -1;
                    end_loc_[i] = -1;
                    size_loc_used_[i] = 0;
                    start_loc_used_[i] = -1;
                    end_loc_used_[i] = -1;
                }
            }

            void free() {
                std::vector<T>().swap(buf_);
            }

            void clear() {
                for (S32 i=0; i < buf_.size(); i++) buf_[i] = 0;
            }

            void resize() {
                buf_.resize(size_loc_tot_);
            }

            template <class U>
            void copySizeInfoFrom(const M2L_Buffer<U> tmp) {
                size_glb_tot_ = tmp.size_glb_tot_;
                size_loc_tot_ = tmp.size_loc_tot_;
                size_loc_used_tot_ = tmp.size_loc_used_tot_;
                for (S32 i=0; i<DIM; i++) {
                    size_glb_[i] = tmp.size_glb_[i];
                    start_glb_[i] = tmp.start_glb_[i];
                    end_glb_[i] = tmp.end_glb_[i];
                    size_loc_[i] = tmp.size_loc_[i];
                    start_loc_[i] = tmp.start_loc_[i];
                    end_loc_[i] = tmp.end_loc_[i];
                    size_loc_used_[i] = tmp.size_loc_used_[i];
                    start_loc_used_[i] = tmp.start_loc_used_[i];
                    end_loc_used_[i] = tmp.end_loc_used_[i];
                }
            }

            void calcSizeGlb() {
                for (S32 i=0; i<DIM; i++) {
                    size_glb_[i] = end_glb_[i] - start_glb_[i] + 1;
                }
            }

            void calcSizeLoc() {
                for (S32 i=0; i<DIM; i++) {
                    size_loc_[i] = end_loc_[i] - start_loc_[i] + 1;
                }
            }

            void calcSizeLocUsed() {
                for (S32 i=0; i<DIM; i++) {
                    size_loc_used_[i] = end_loc_used_[i] - start_loc_used_[i] + 1;
                }
            }

            void calcSizeTot() {
                size_glb_tot_ = 1;
                size_loc_tot_ = 1;
                size_loc_used_tot_ = 1;
                for (S32 i=0; i<DIM; i++) {
                    size_glb_tot_ *= size_glb_[i];
                    size_loc_tot_ *= size_loc_[i];
                    size_loc_used_tot_ *= size_loc_used_[i];
                }
            }

            void outputSizeInfo(std::ostream & fout=std::cout) const {
                fout << "size_glb_tot_ = " << size_glb_tot_ << std::endl;
                for (S32 i=0; i<DIM; i++) {
                    fout << start_glb_[i] << "   "
                         << end_glb_[i] << "   "
                         << size_glb_[i] << std::endl;
                }
                fout << "size_loc_tot_ = " << size_loc_tot_ << std::endl;
                for (S32 i=0; i<DIM; i++) {
                    fout << start_loc_[i] << "   "
                         << end_loc_[i] << "   " 
                         << size_loc_[i] << std::endl;
                }
                fout << "size_loc_used_tot_ = " << size_loc_used_tot_ << std::endl;
                for (S32 i=0; i<DIM; i++) {
                    fout << start_loc_used_[i] << "   "
                         << end_loc_used_[i] << "   "
                         << size_loc_used_[i] << std::endl;
                }
            }

        };


        class M2L_Engine {
        public:
            using real_t = double;
            using cplx_t = std::complex<real_t>;

            Parameters param_;
            Parameters param_prev_;
            bool first_call_of_fftw_mpi_init_;
            bool first_call_of_fftw_init_threads_;
            bool first_call_of_convolution_;

            // Variables for communication
            S32 mode_;
            S32 n_proc_in_parent_group_;
            S32 rank_in_parent_group_;
            S32 n_group_;
            S32 n_proc_in_my_group_;
            S32 rank_in_my_group_;
            S32 n_proc_min_in_group_;
            std::vector<S32> rank_start_, rank_end_;
            std::vector<S32> lm_start_, lm_end_;
            std::vector<S32> local_0_start_, local_0_end_; // z of (z,y,x)
            std::vector<S32> local_1_start_, local_1_end_; // y of (y,z,x)
            S32 n_proc_for_fft_;
            S32 idx_to_my_group_;
            bool is_mpifft_usable_;
            bool use_mpifft_;
            S32 fft_size_crit_;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Group parent_group_;
            MPI_Comm parent_comm_;
            MPI_Datatype mpi_real_t_;
            MPI_Datatype mpi_cplx_t_;
            MPI_Group group_all_;
            MPI_Group group_fft_;
            MPI_Group group_int_; // inter-group
            MPI_Comm comm_all_;
            MPI_Comm comm_fft_;
            MPI_Comm comm_int_;
#endif

            // Variables for M2L and green function
            M2L_Buffer<real_t> mm_r_, le_r_, gf_r_;
            M2L_Buffer<cplx_t> mm_k_, le_k_, gf_k_;

            // Variables for FFT
            S32 fft_mode_;
            S32 fft_mode_prev_;
            S32 size_fft_[3];
            S32 size_fft_prev_[3];
            fftw_real_t *rbuf_;
            fftw_cplx_t *kbuf_;
            std::vector<fftw_plan> plan_fwd_;
            std::vector<fftw_plan> plan_bkw_;

            // Timing measurement
            TimeProfilePMM time_profile_;

        
            M2L_Engine() {
                first_call_of_fftw_mpi_init_ = true;
                first_call_of_fftw_init_threads_ = true;
                first_call_of_convolution_ = true;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                group_all_ = MPI_GROUP_NULL;
                group_fft_ = MPI_GROUP_NULL;
                group_int_ = MPI_GROUP_NULL;
                comm_all_ = MPI_COMM_NULL;
                comm_fft_ = MPI_COMM_NULL;
                comm_int_ = MPI_COMM_NULL;
#endif
                fft_mode_ = -1;
                fft_mode_prev_ = -1;
                for (S32 i = 0; i < 3; i++) {
                    size_fft_[i] = 0;
                    size_fft_prev_[i] = 0;
                }
                rbuf_ = nullptr;
                kbuf_ = nullptr;
            }
        
            void finalize() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                if (group_all_ != MPI_GROUP_NULL) MPI_Group_free(&group_all_);
                if (group_fft_ != MPI_GROUP_NULL) MPI_Group_free(&group_fft_);
                if (group_int_ != MPI_GROUP_NULL) MPI_Group_free(&group_int_);
                if (comm_all_ != MPI_COMM_NULL) MPI_Comm_free(&comm_all_);
                if (comm_fft_ != MPI_COMM_NULL) MPI_Comm_free(&comm_fft_);
                if (comm_int_ != MPI_COMM_NULL) MPI_Comm_free(&comm_int_);
#endif
                if (rbuf_ != nullptr) fftw_free(rbuf_);
                if (kbuf_ != nullptr) fftw_free(kbuf_);
                for (S32 i = 0; i < plan_fwd_.size(); i++) fftw_destroy_plan(plan_fwd_[i]);
                for (S32 i = 0; i < plan_bkw_.size(); i++) fftw_destroy_plan(plan_bkw_[i]);
            }
       
            void initialize(const Parameters & param,
                            const bool use_mpifft_if_possible,
                            const S32 fft_size_crit,
                            const bool debug_flag = false) {
                F64 time_start = GetWtime();
                // Copy parameters, etc.
                param_prev_ = param_; // save the previous parameters
                param_ = param;
                fft_size_crit_ = fft_size_crit;
                // Initialize 
                if ((param_ != param_prev_) || 
                    (use_mpifft_ != use_mpifft_if_possible)) {
                    const S32 n_mm_compo = (param_.p + 1) * (param_.p + 1);
                    const S32vec n_cell  = param_.n_cell;
                    const S32 bc         = param_.bc;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                    // Set parent_comm
                    parent_comm_ = MPI_COMM_WORLD;
                    // Get information of parent_comm
                    S32 n_proc, my_rank;
                    MPI_Comm_size(parent_comm_, &n_proc);
                    MPI_Comm_rank(parent_comm_, &my_rank);
                    MPI_Comm_group(parent_comm_, &parent_group_);
                    n_proc_in_parent_group_ = n_proc;
                    rank_in_parent_group_ = my_rank;
                    // Get MPI_Datatype for ptrdiff_t, real_t, cplx_t
                    mpi_real_t_ = GetDataType<real_t>();
                    mpi_cplx_t_ = GetDataType<cplx_t>();
                    // Divide parent_comm into groups
                    if (n_proc >= n_mm_compo) {
                        mode_ = 1;
                        n_group_ = n_mm_compo;
                        // Set rank_start_, rank_end_
                        n_proc_min_in_group_ = n_proc/n_group_;
                        S32 n_rem = n_proc - n_proc_min_in_group_ * n_group_;
                        rank_start_.resize(n_group_);
                        rank_end_.resize(n_group_);
                        rank_start_[0] = 0;
                        for (S32 i = 0; i < n_group_; i++) {
                            S32 n_proc_in_group = n_proc_min_in_group_;
                            if (i < n_rem) n_proc_in_group++;
                            rank_end_[i] = rank_start_[i] + (n_proc_in_group - 1);
                            if (i < n_group_ - 1) rank_start_[i+1] = rank_end_[i] + 1;
                        }
                        if (debug_flag) {
                            if (my_rank == 0) {
                                std::cout << "### Information of groups for M2L ###" << std::endl;
                                std::cout << "(n_proc = " << n_proc << ")" << std::endl;
                                for (S32 i = 0; i < n_group_; i++) {
                                    std::cout << "i = " << i
                                              << " rank_start = " << rank_start_[i]
                                              << " rank_end = " << rank_end_[i]
                                              << std::endl;
                                }
                            }
                        }
                        // Calculate idx_to_my_group, n_proc_in_my_group, rank_in_my_group
                        for (S32 i = 0; i < n_group_; i++)
                            if ((rank_start_[i] <= my_rank) && (my_rank <= rank_end_[i])) {
                                idx_to_my_group_ = i;
                                n_proc_in_my_group_ = rank_end_[i] - rank_start_[i] + 1;
                                rank_in_my_group_ = my_rank - rank_start_[i];
                            }
                        // Set lm_start_, lm_end_
                        lm_start_.resize(n_group_);
                        lm_end_.resize(n_group_);
                        for (S32 i = 0; i < n_group_; i++) {
                            lm_start_[i] = i;
                            lm_end_[i] = i;
                            // Each group is responsible for a single (l,m)
                        }
                        // Make group_all & comm_all
                        if (group_all_ != MPI_GROUP_NULL) MPI_Group_free(&group_all_);
                        if (comm_all_ != MPI_COMM_NULL) MPI_Comm_free(&comm_all_);
                        std::vector<S32> ranks;
                        S32 rnk_start = rank_start_[idx_to_my_group_];
                        S32 rnk_end = rank_end_[idx_to_my_group_];
                        for (S32 rnk = rnk_start; rnk <= rnk_end; rnk++) ranks.push_back(rnk);
                        MPI_Group_incl(parent_group_, n_proc_in_my_group_, &ranks[0], &group_all_);
                        MPI_Comm_create(parent_comm_, group_all_, &comm_all_);
                        // Calculate n_proc_for_fft
                        S32 n_proc_max_for_fft;
                        if (bc == BOUNDARY_CONDITION_OPEN) n_proc_max_for_fft = n_cell.z;
                        else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) n_proc_max_for_fft = n_cell.z/2;
                        else assert(false);
                        if (n_proc_min_in_group_ <= n_proc_max_for_fft) n_proc_for_fft_ = n_proc_min_in_group_;
                        else n_proc_for_fft_ = n_proc_max_for_fft;
                        // Make group_fft, comm_fft
                        if (n_proc_for_fft_ > 1) is_mpifft_usable_ = true;
                        else is_mpifft_usable_ = false;
                        if (is_mpifft_usable_) use_mpifft_ = use_mpifft_if_possible; // specified by user.
                        else use_mpifft_ = false;
                        if (use_mpifft_) {
                            if (first_call_of_fftw_mpi_init_) {
                                fftw_mpi_init();
                                first_call_of_fftw_mpi_init_ = false;
                            }
                            if (group_fft_ != MPI_GROUP_NULL) MPI_Group_free(&group_fft_);
                            if (comm_fft_ != MPI_COMM_NULL) MPI_Comm_free(&comm_fft_);
                            if (rank_in_my_group_ < n_proc_for_fft_) {
                                ranks.resize(n_proc_for_fft_);
                                for (S32 i = 0; i < n_proc_for_fft_; i++)
                                    ranks[i] = rank_start_[idx_to_my_group_] + i;
                                MPI_Group_incl(parent_group_, n_proc_for_fft_, &ranks[0], &group_fft_);
                                MPI_Comm_create(parent_comm_, group_fft_, &comm_fft_);
                            } else {
                                group_fft_ = MPI_GROUP_NULL;
                                MPI_Comm_create(parent_comm_, MPI_GROUP_EMPTY, &comm_fft_);
                            }
                            // In this case, we need to calculate local_0_start & local_0_end
                            // using FFTW API.
                            // Resize local_?_start, local_?_end 
                            local_0_start_.resize(n_proc_for_fft_);
                            local_0_end_.resize(n_proc_for_fft_);
                            local_1_start_.resize(n_proc_for_fft_);
                            local_1_end_.resize(n_proc_for_fft_);
                            // Calculate loca_?_start, local_?_end
                            if (comm_fft_ != MPI_COMM_NULL) {
                                ptrdiff_t alloc_local, local_n0, local_0_start;
                                ptrdiff_t local_n1, local_1_start;
                                if (bc == BOUNDARY_CONDITION_OPEN) {
                                    alloc_local = fftw_mpi_local_size_3d_transposed(2*n_cell.z, 2*n_cell.y, n_cell.x+1,
                                                                                    comm_fft_,
                                                                                    &local_n0, &local_0_start,
                                                                                    &local_n1, &local_1_start);
                                } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                                    alloc_local = fftw_mpi_local_size_3d_transposed(n_cell.z, n_cell.y, n_cell.x/2+1,
                                                                                    comm_fft_,
                                                                                    &local_n0, &local_0_start,
                                                                                    &local_n1, &local_1_start);
                                } else {
                                    if (my_rank == 0) 
                                        std::cout << "This boundary condition is not supported yet." << std::endl;
                                    Abort(-1);
                                }
                                S32 sendbuf[4];
                                sendbuf[0] = local_0_start;
                                sendbuf[1] = local_n0;
                                sendbuf[2] = local_1_start;
                                sendbuf[3] = local_n1;
                                S32 * recvbuf = new S32[4 * n_proc_for_fft_];
                                MPI_Gather(sendbuf, 4, MPI_INT,
                                           recvbuf, 4, MPI_INT,
                                           0, comm_fft_);
                                if (rank_in_my_group_ == 0) {
                                    for (S32 i = 0; i < n_proc_for_fft_; i++) {
                                        local_0_start_[i] = recvbuf[4 * i];
                                        local_0_end_[i]   = local_0_start_[i] + recvbuf[4 * i + 1] - 1;
                                        local_1_start_[i] = recvbuf[4 * i + 2];
                                        local_1_end_[i]   = local_1_start_[i] + recvbuf[4 * i + 3] - 1;
                                    }
                                }
                                delete [] recvbuf;
                            }
                            MPI_Bcast(&local_0_start_[0], n_proc_for_fft_, MPI_INT, 0, comm_all_);
                            MPI_Bcast(&local_0_end_[0], n_proc_for_fft_, MPI_INT, 0, comm_all_);
                            MPI_Bcast(&local_1_start_[0], n_proc_for_fft_, MPI_INT, 0, comm_all_);
                            MPI_Bcast(&local_1_end_[0], n_proc_for_fft_, MPI_INT, 0, comm_all_);
                            // Check
                            if (rank_in_parent_group_ == 0) {
                                std::cout << "### Information of decomposition of FFTW ###" << std::endl;
                                std::cout << "(n_cell.z = " << n_cell.z << ")" << std::endl;
                                for (S32 i = 0; i < n_proc_for_fft_; i++) {
                                    std::cout << "i = " << i
                                              << " local_0_start = " << local_0_start_[i]
                                              << " local_0_end = " << local_0_end_[i]
                                              << " local_1_start = " << local_1_start_[i]
                                              << " local_1_end = " << local_1_end_[i]
                                              << std::endl;
                                }
                            }
                        } else {
                            n_proc_for_fft_ = 1; // reset
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                            if (first_call_of_fftw_init_threads_) {
                                fftw_init_threads();
                                first_call_of_fftw_init_threads_ = false;
                            }
#endif
                            group_fft_ = MPI_GROUP_NULL;
                            comm_fft_ = MPI_COMM_NULL;
                            local_0_start_.resize(1);
                            local_0_end_.resize(1);
                            if (bc == BOUNDARY_CONDITION_OPEN) {
                                local_0_start_[0] = 0;
                                local_0_end_[0] = 2 * n_cell.z - 1;
                            } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                                local_0_start_[0] = 0;
                                local_0_end_[0] = n_cell.z - 1;
                            }
                        }
                        // Make group_int, comm_int
                        if (group_int_ != MPI_GROUP_NULL) MPI_Group_free(&group_int_);
                        if (comm_int_ != MPI_COMM_NULL) MPI_Comm_free(&comm_int_);
                        if (rank_in_my_group_ < n_proc_for_fft_) {
                            ranks.resize(n_group_);
                            for (S32 i = 0; i < n_group_; i++)
                                ranks[i] = rank_start_[i] + rank_in_my_group_;
                            MPI_Group_incl(parent_group_, n_group_, &ranks[0], &group_int_);
                            MPI_Comm_create(parent_comm_, group_int_, &comm_int_);
                        } else {
                            group_int_ = MPI_GROUP_NULL;
                            MPI_Comm_create(parent_comm_, MPI_GROUP_EMPTY, &comm_int_);
                        }
                    } else {
                        mode_ = 2;
                        n_group_ = n_proc;
                        n_proc_in_my_group_ = 1; 
                        rank_in_my_group_ = 0; 
                        n_proc_min_in_group_ = 1;
                        n_proc_for_fft_ = 1;
                        idx_to_my_group_ = my_rank;
                        is_mpifft_usable_ = false;
                        use_mpifft_ = false;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                        if (first_call_of_fftw_init_threads_) {
                            fftw_init_threads();
                            first_call_of_fftw_init_threads_ = false;
                        }
#endif
                        // Set rank_start_, rank_end_
                        rank_start_.resize(n_group_);
                        rank_end_.resize(n_group_);
                        for (S32 i = 0; i < n_group_; i++) {
                            rank_start_[i] = i;
                            rank_end_[i] = i;
                        }
                        // Set lm_start_, lm_end_
                        const S32 n_mm_compo_per_proc_min = n_mm_compo/n_proc;
                        S32 n_rem = n_mm_compo - n_mm_compo_per_proc_min * n_proc;
                        lm_start_.resize(n_proc);
                        lm_end_.resize(n_proc);
                        lm_start_[0] = 0;
                        for (S32 i = 0; i < n_proc; i++) {
                            S32 n_mm_compo_per_proc = n_mm_compo_per_proc_min;
                            if (i < n_rem) n_mm_compo_per_proc++;
                            lm_end_[i] = lm_start_[i] + (n_mm_compo_per_proc - 1);
                            if (i < n_proc-1) lm_start_[i+1] = lm_end_[i] + 1;
                        }
                        if (rank_in_parent_group_ == 0) {
                            std::cout << "### Information of decomposition of (l,m) ###" << std::endl;
                            std::cout << "(n_mm_compo = " << n_mm_compo << ")" << std::endl;
                            for (S32 i = 0; i < n_proc; i++) {
                                std::cout << "i = " << i
                                          << " local_0_start = " << lm_start_[i]
                                          << " local_0_end = " << lm_end_[i]
                                          << std::endl;
                            }
                        }
                        // Set local_0_start_, local_0_end_
                        group_fft_ = MPI_GROUP_NULL;
                        comm_fft_ = MPI_COMM_NULL;
                        local_0_start_.resize(n_proc_for_fft_);
                        local_0_end_.resize(n_proc_for_fft_);
                        if (bc == BOUNDARY_CONDITION_OPEN) {
                            local_0_start_[0] = 0;
                            local_0_end_[0] = 2 * n_cell.z - 1;
                        } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                            local_0_start_[0] = 0;
                            local_0_end_[0] = n_cell.z - 1;
                        }
                        // Note that we do not need to set local_1_start_ 
                        // and local_1_end_ because of use_mpifft = false.
                        // (hence, both arrays are not accessed)
                    }
#else // PARTICLE_SIMULATOR_MPI_PARALLEL
                    // Sequential execution/ OpenMP parallelization
                    mode_ = 0;
                    n_proc_in_parent_group_ = 1;
                    rank_in_parent_group_ = 0;
                    n_group_ = 1;
                    n_proc_in_my_group_ = 1;
                    rank_in_my_group_ = 0;
                    n_proc_min_in_group_ = 1;
                    n_proc_for_fft_ = 1;
                    idx_to_my_group_ = 0;
                    is_mpifft_usable_ = false;
                    use_mpifft_ = false;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                    if (first_call_of_fftw_init_threads_) {
                        fftw_init_threads();
                        first_call_of_fftw_init_threads_ = false;
                    }
#endif
                    // Set rank_start_(n_group_), rank_end_(n_group_) [idx_to_my_group_]
                    rank_start_.resize(n_group_);
                    rank_end_.resize(n_group_);
                    rank_start_[0] = 0;
                    rank_end_[0] = 0;
                    // Set lm_start_(n_group_), lm_end_(n_group_) [idx_to_my_group_]
                    lm_start_.resize(n_group_);
                    lm_end_.resize(n_group_);
                    lm_start_[0] = 0;
                    lm_end_[0] = n_mm_compo - 1;
                    // Set local_0_start__, local_0_end_ 
                    local_0_start_.resize(n_proc_for_fft_);
                    local_0_end_.resize(n_proc_for_fft_);
                    if (bc == BOUNDARY_CONDITION_OPEN) {
                        local_0_start_[0] = 0;
                        local_0_end_[0] = 2 * n_cell.z - 1;
                    } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                        local_0_start_[0] = 0;
                        local_0_end_[0] = n_cell.z - 1;
                    }
#endif // PARTICLE_SIMULATOR_MPI_PARALLEL

                    // Set the size information of mm_r, le_r, mm_k, le_k
                    if (bc == BOUNDARY_CONDITION_OPEN) {
                        // Set the size information of mm_r
                        mm_r_.start_glb_[0] = 0;
                        mm_r_.start_glb_[1] = 0;
                        mm_r_.start_glb_[2] = 0;
                        mm_r_.start_glb_[3] = 0;
                        mm_r_.end_glb_[0] = n_mm_compo - 1;
                        mm_r_.end_glb_[1] = 2 * n_cell.x - 1;
                        mm_r_.end_glb_[2] = 2 * n_cell.y - 1;
                        mm_r_.end_glb_[3] = 2 * n_cell.z - 1;
                        mm_r_.calcSizeGlb();
                        if (rank_in_my_group_ < n_proc_for_fft_) {
                            mm_r_.start_loc_[0] = lm_start_[idx_to_my_group_]; 
                            mm_r_.start_loc_[1] = 0;
                            mm_r_.start_loc_[2] = 0;
                            mm_r_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                            mm_r_.end_loc_[0] = lm_end_[idx_to_my_group_];
                            mm_r_.end_loc_[1] = 2 * n_cell.x - 1;
                            mm_r_.end_loc_[2] = 2 * n_cell.y - 1;
                            mm_r_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                            mm_r_.calcSizeLoc();
                            for (S32 i = 0; i < 4; i++) {
                                mm_r_.start_loc_used_[i] = mm_r_.start_loc_[i];
                                mm_r_.end_loc_used_[i] = mm_r_.end_loc_[i];
                            }
                            mm_r_.calcSizeLocUsed();
                        }
                        // Set the size information of le_k
                        if (use_mpifft_) {
                            // In this case, we set the size information 
                            // basend on FFTW_MPI_TRANSPOSED_OUT format.
                            le_k_.start_glb_[0] = 0;
                            le_k_.start_glb_[1] = 0;
                            le_k_.start_glb_[2] = 0;
                            le_k_.start_glb_[3] = 0;
                            le_k_.end_glb_[0] = n_mm_compo - 1;
                            le_k_.end_glb_[1] = 1 + n_cell.x - 1;
                            le_k_.end_glb_[2] = 2 * n_cell.z - 1;
                            le_k_.end_glb_[3] = 2 * n_cell.y - 1;
                            le_k_.calcSizeGlb();
                            if (rank_in_my_group_ < n_proc_for_fft_) {
                                le_k_.start_loc_[0] = lm_start_[idx_to_my_group_]; // only one (l,m)
                                le_k_.start_loc_[1] = 0;
                                le_k_.start_loc_[2] = 0;
                                le_k_.start_loc_[3] = local_1_start_[rank_in_my_group_];
                                le_k_.end_loc_[0] = lm_end_[idx_to_my_group_];
                                le_k_.end_loc_[1] = n_cell.x;
                                le_k_.end_loc_[2] = 2 * n_cell.z - 1;
                                le_k_.end_loc_[3] = local_1_end_[rank_in_my_group_];
                                le_k_.calcSizeLoc();
                                for (S32 i = 0; i < 4; i++) { 
                                    le_k_.start_loc_used_[i] = le_k_.start_loc_[i];
                                    le_k_.end_loc_used_[i] = le_k_.end_loc_[i];
                                }
                                le_k_.calcSizeLocUsed();
                            }
                        } else {
                            // In this case, we set the size information 
                            // in the same way as le_r. 
                            le_k_.start_glb_[0] = 0;
                            le_k_.start_glb_[1] = 0;
                            le_k_.start_glb_[2] = 0;
                            le_k_.start_glb_[3] = 0;
                            le_k_.end_glb_[0] = n_mm_compo - 1;
                            le_k_.end_glb_[1] = 1 + n_cell.x - 1;
                            le_k_.end_glb_[2] = 2 * n_cell.y - 1;
                            le_k_.end_glb_[3] = 2 * n_cell.z - 1;
                            le_k_.calcSizeGlb();
                            if (rank_in_my_group_ < n_proc_for_fft_) {
                                le_k_.start_loc_[0] = lm_start_[idx_to_my_group_]; 
                                le_k_.start_loc_[1] = 0;
                                le_k_.start_loc_[2] = 0;
                                le_k_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                                le_k_.end_loc_[0] = lm_end_[idx_to_my_group_];
                                le_k_.end_loc_[1] = n_cell.x;
                                le_k_.end_loc_[2] = 2 * n_cell.y - 1;
                                le_k_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                                le_k_.calcSizeLoc();
                                for (S32 i = 0; i < 4; i++) {
                                    le_k_.start_loc_used_[i] = le_k_.start_loc_[i];
                                    le_k_.end_loc_used_[i] = le_k_.end_loc_[i];
                                }
                                le_k_.calcSizeLocUsed();
                            }
                        }
                    } else if (bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                        // Set the size information of mm_r
                        mm_r_.start_glb_[0] = 0;
                        mm_r_.start_glb_[1] = 0;
                        mm_r_.start_glb_[2] = 0;
                        mm_r_.start_glb_[3] = 0;
                        mm_r_.end_glb_[0] = n_mm_compo - 1;
                        mm_r_.end_glb_[1] = n_cell.x - 1;
                        mm_r_.end_glb_[2] = n_cell.y - 1;
                        mm_r_.end_glb_[3] = n_cell.z - 1;
                        mm_r_.calcSizeGlb();
                        if (rank_in_my_group_ < n_proc_for_fft_) {
                            mm_r_.start_loc_[0] = lm_start_[idx_to_my_group_];
                            mm_r_.start_loc_[1] = 0;
                            mm_r_.start_loc_[2] = 0;
                            mm_r_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                            mm_r_.end_loc_[0] = lm_end_[idx_to_my_group_];
                            mm_r_.end_loc_[1] = n_cell.x - 1;
                            mm_r_.end_loc_[2] = n_cell.y - 1;
                            mm_r_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                            mm_r_.calcSizeLoc();
                            for (S32 i = 0; i < 4; i++) {
                                mm_r_.start_loc_used_[i] = mm_r_.start_loc_[i];
                                mm_r_.end_loc_used_[i] = mm_r_.end_loc_[i];
                            }
                            mm_r_.calcSizeLocUsed();
                        }
                        // Set the size information of le_k
                        if (use_mpifft_) {
                            // In this case, we set the size information
                            // basend on FFTW_MPI_TRANSPOSED_OUT format.
                            le_k_.start_glb_[0] = 0;
                            le_k_.start_glb_[1] = 0;
                            le_k_.start_glb_[2] = 0;
                            le_k_.start_glb_[3] = 0;
                            le_k_.end_glb_[0] = n_mm_compo - 1;
                            le_k_.end_glb_[1] = 1 + n_cell.x/2 - 1;
                            le_k_.end_glb_[2] = n_cell.z - 1;
                            le_k_.end_glb_[3] = n_cell.y - 1;
                            le_k_.calcSizeGlb();
                            if (rank_in_my_group_ < n_proc_for_fft_) {
                                le_k_.start_loc_[0] = lm_start_[idx_to_my_group_]; // only one (l,m)
                                le_k_.start_loc_[1] = 0;
                                le_k_.start_loc_[2] = 0;
                                le_k_.start_loc_[3] = local_1_start_[rank_in_my_group_];
                                le_k_.end_loc_[0] = lm_end_[idx_to_my_group_];
                                le_k_.end_loc_[1] = n_cell.x/2;
                                le_k_.end_loc_[2] = n_cell.z - 1;
                                le_k_.end_loc_[3] = local_1_end_[rank_in_my_group_];
                                le_k_.calcSizeLoc();
                                for (S32 i = 0; i < 4; i++) {
                                    le_k_.start_loc_used_[i] = le_k_.start_loc_[i];
                                    le_k_.end_loc_used_[i] = le_k_.end_loc_[i];
                                }
                                le_k_.calcSizeLocUsed();
                            }
                        } else {
                            // In this case, we set the size information 
                            // in the same way as le_r.
                            le_k_.start_glb_[0] = 0;
                            le_k_.start_glb_[1] = 0;
                            le_k_.start_glb_[2] = 0;
                            le_k_.start_glb_[3] = 0;
                            le_k_.end_glb_[0] = n_mm_compo - 1;
                            le_k_.end_glb_[1] = 1 + n_cell.x/2 - 1;
                            le_k_.end_glb_[2] = n_cell.y - 1;
                            le_k_.end_glb_[3] = n_cell.z - 1;
                            le_k_.calcSizeGlb();
                            if (rank_in_my_group_ < n_proc_for_fft_) {
                                le_k_.start_loc_[0] = lm_start_[idx_to_my_group_]; 
                                le_k_.start_loc_[1] = 0;
                                le_k_.start_loc_[2] = 0;
                                le_k_.start_loc_[3] = local_0_start_[rank_in_my_group_];
                                le_k_.end_loc_[0] = lm_end_[idx_to_my_group_];
                                le_k_.end_loc_[1] = n_cell.x/2;
                                le_k_.end_loc_[2] = n_cell.y - 1;
                                le_k_.end_loc_[3] = local_0_end_[rank_in_my_group_];
                                le_k_.calcSizeLoc();
                                for (S32 i = 0; i < 4; i++) {
                                    le_k_.start_loc_used_[i] = le_k_.start_loc_[i];
                                    le_k_.end_loc_used_[i] = le_k_.end_loc_[i];
                                }
                                le_k_.calcSizeLocUsed();
                            }
                        }
                    }
                    // Set the size information of le_r_
                    le_r_.copySizeInfoFrom(mm_r_);
                    // Set the size information of mm_k_
                    mm_k_.copySizeInfoFrom(le_k_);
                    mm_k_.start_loc_[0] = 0; // overwrite
                    mm_k_.end_loc_[0] = n_mm_compo - 1; // overwrite
                    mm_k_.calcSizeLoc();
                    // Calculate total sizes
                    mm_r_.calcSizeTot();
                    mm_k_.calcSizeTot();
                    le_k_.calcSizeTot();
                    le_r_.calcSizeTot();
                }
                time_profile_.M2L_initialize += GetWtime() - time_start;
            }

            void setGreenFunction(const F64vec width_cell) {

                if (param_ != param_prev_) {

                    Comm::barrier();
                    F64 time_start = GetWtime();

                    // Extract basic parameters
                    const S32 p = param_.p;
                    const S32 bc = param_.bc;
                    const S32vec n_cell = param_.n_cell;
                    
                    // Free 
                    gf_r_.free();
                    gf_k_.free();

                    // Memory allocation
                    gf_r_.copySizeInfoFrom(mm_r_);
                    gf_k_.copySizeInfoFrom(mm_k_);
                    const S32 LEN2 = (2*p + 1)*(2*p + 1);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                    if (mode_ == 1 && use_mpifft_) {
                        // Set the size information of gf_r[0]
                        gf_r_.start_glb_[0] = 0;
                        gf_r_.end_glb_[0] = LEN2 - 1;
                        gf_r_.calcSizeGlb();
                        if (rank_in_my_group_ < n_proc_for_fft_) {
                            gf_r_.start_loc_[0] = 0; 
                            gf_r_.end_loc_[0] = LEN2 - 1; // all (l,m)
                            gf_r_.calcSizeLoc();
                            gf_r_.start_loc_used_[0] = gf_r_.start_loc_[0];
                            gf_r_.end_loc_used_[0] = gf_r_.end_loc_[0];
                            gf_r_.calcSizeLocUsed();
                        }
                        gf_r_.calcSizeTot();
                        // Set the size of gf_k[0]
                        gf_k_.start_glb_[0] = 0;
                        gf_k_.end_glb_[0] = LEN2 - 1;
                        gf_k_.calcSizeGlb();
                        if (rank_in_my_group_ < n_proc_for_fft_) {
                            gf_k_.start_loc_[0] = 0; 
                            gf_k_.end_loc_[0] = LEN2 - 1; // all (l,m)
                            gf_k_.calcSizeLoc();
                            gf_k_.start_loc_used_[0] = gf_k_.start_loc_[0];
                            gf_k_.end_loc_used_[0] = gf_k_.end_loc_[0];
                            gf_k_.calcSizeLocUsed();
                        }
                        gf_k_.calcSizeTot();
                    } else {
#endif // PARTICLE_SIMULATOR_MPI_PARALLEL
                        // Set the size information of gf_r[0]
                        gf_r_.start_glb_[0] = 0;
                        gf_r_.end_glb_[0] = LEN2 - 1;
                        gf_r_.calcSizeGlb();
                        gf_r_.start_loc_[0] = 0; 
                        gf_r_.end_loc_[0]   = LEN2 - 1; // all (l,m)
                        gf_r_.calcSizeLoc();
                        gf_r_.start_loc_used_[0] = gf_r_.start_loc_[0];
                        gf_r_.end_loc_used_[0] = gf_r_.end_loc_[0];
                        gf_r_.calcSizeLocUsed();
                        gf_r_.calcSizeTot();
                        // Set the size of gf_k[0]
                        gf_k_.start_glb_[0] = 0;
                        gf_k_.end_glb_[0] = LEN2 - 1;
                        gf_k_.calcSizeGlb();
                        gf_k_.start_loc_[0] = 0; 
                        gf_k_.end_loc_[0]   = LEN2 - 1; // all (l,m)
                        gf_k_.calcSizeLoc();
                        gf_k_.start_loc_used_[0] = gf_k_.start_loc_[0];
                        gf_k_.end_loc_used_[0] = gf_k_.end_loc_[0];
                        gf_k_.calcSizeLocUsed();
                        gf_k_.calcSizeTot();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                    }
#endif // PARTICLE_SIMULATOR_MPI_PARALLEL
                    // check
                    //gf_r_.outputSizeInfo();
                    //gf_k_.outputSizeInfo();
                    // allocate 
                    gf_r_.resize();
                    gf_k_.resize();

                    Comm::barrier();
                    time_profile_.M2L_gf_initialize += GetWtime() - time_start;

                    // Calculate gf_r
                    Comm::barrier();
                    time_start = GetWtime();

                    calcGreenFunctionInRealSpace(width_cell);

                    Comm::barrier();
                    time_profile_.M2L_gf_calc_gf_r += GetWtime() - time_start;

#if 0
                    // Check
                    {
                        std::stringstream ss;
                        ss << "gf_r_" << std::setfill('0') << std::setw(5)
                           << rank_in_parent_group_ << ".txt";
                        const std::string filename = ss.str();
                        std::ofstream output_file;
                        output_file.open(filename.c_str(), std::ios::trunc);
                        for (S32 k = gf_r_.start_loc_used_[3]; k <= gf_r_.end_loc_used_[3]; k++)
                            for (S32 j = gf_r_.start_loc_used_[2]; j <= gf_r_.end_loc_used_[2]; j++)
                                for (S32 i = gf_r_.start_loc_used_[1]; i <= gf_r_.end_loc_used_[1]; i++)
                                    for (S32 lm = gf_r_.start_loc_used_[0]; lm <= gf_r_.end_loc_used_[0]; lm++)
                        {
                            const S32 lm_loc = lm - gf_r_.start_loc_[0];
                            const S32 i_loc  = i  - gf_r_.start_loc_[1];
                            const S32 j_loc  = j  - gf_r_.start_loc_[2];
                            const S32 k_loc  = k  - gf_r_.start_loc_[3];
                            const S32 adr = lm_loc
                                          + gf_r_.size_loc_[0] * (i_loc
                                          + gf_r_.size_loc_[1] * (j_loc
                                          + gf_r_.size_loc_[2] * k_loc));
                            const S32 idx = lm
                                          + gf_r_.size_glb_[0] * (i
                                          + gf_r_.size_glb_[1] * (j
                                          + gf_r_.size_glb_[3] * k));
                            output_file << idx << "    " << gf_r_.buf_[adr] << std::endl;
                        }
                        output_file.close();
                    }
                    Finalize();
                    std::exit(0);
#endif


                    // Calculate gf_k using FFT
                    Comm::barrier();
                    time_start = GetWtime();

                    calcGreenFunctionInWavenumberSpace();

                    Comm::barrier();
                    time_profile_.M2L_gf_calc_gf_k += GetWtime() - time_start;

#if 0
                    // Check
                    {
                        std::stringstream ss;
                        ss << "gf_k_" << std::setfill('0') << std::setw(5)
                           << rank_in_parent_group_ << ".txt";
                        const std::string filename = ss.str();
                        std::ofstream output_file;
                        output_file.open(filename.c_str(), std::ios::trunc);
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                        if (use_mpifft_) {
                            for (S32 j = gf_k_.start_loc_used_[3]; j <= gf_k_.end_loc_used_[3]; j++)
                                for (S32 k = gf_k_.start_loc_used_[2]; k <= gf_k_.end_loc_used_[2]; k++)
                                    for (S32 i = gf_k_.start_loc_used_[1]; i <= gf_k_.end_loc_used_[1]; i++)
                                        for (S32 lm = gf_k_.start_loc_used_[0]; lm <= gf_k_.end_loc_used_[0]; lm++)
                            {
                                const S32 lm_loc = lm - gf_k_.start_loc_[0];
                                const S32 i_loc  = i  - gf_k_.start_loc_[1];
                                const S32 k_loc  = k  - gf_k_.start_loc_[2];
                                const S32 j_loc  = j  - gf_k_.start_loc_[3];
                                const S32 adr = lm_loc
                                              + gf_k_.size_loc_[0] * (i_loc
                                              + gf_k_.size_loc_[1] * (k_loc
                                              + gf_k_.size_loc_[2] * j_loc));
                                const S32 idx = lm 
                                              + gf_k_.size_glb_[0] * (i
                                              + gf_k_.size_glb_[1] * (j
                                              + gf_k_.size_glb_[3] * k));
                                // Note that the definition of idx must be the same as 
                                // the case of use_mpifft=false to compare the calculated
                                // result with that of the original PMMM code.
                                output_file << idx << "    " << gf_k_.buf_[adr] << std::endl;
                            }
                        } else {
#endif
                            for (S32 k = gf_k_.start_loc_used_[3]; k <= gf_k_.end_loc_used_[3]; k++)
                                for (S32 j = gf_k_.start_loc_used_[2]; j <= gf_k_.end_loc_used_[2]; j++)
                                    for (S32 i = gf_k_.start_loc_used_[1]; i <= gf_k_.end_loc_used_[1]; i++)
                                        for (S32 lm = gf_k_.start_loc_used_[0]; lm <= gf_k_.end_loc_used_[0]; lm++)
                            {
                                const S32 lm_loc = lm - gf_k_.start_loc_[0];
                                const S32 i_loc  = i  - gf_k_.start_loc_[1];
                                const S32 j_loc  = j  - gf_k_.start_loc_[2];
                                const S32 k_loc  = k  - gf_k_.start_loc_[3];
                                const S32 adr = lm_loc
                                              + gf_k_.size_loc_[0] * (i_loc
                                              + gf_k_.size_loc_[1] * (j_loc
                                              + gf_k_.size_loc_[2] * k_loc));
                                const S32 idx = lm 
                                              + gf_k_.size_glb_[0] * (i
                                              + gf_k_.size_glb_[1] * (j
                                              + gf_k_.size_glb_[2] * k));
                                output_file << idx << "    " << gf_k_.buf_[adr] << std::endl;
                            }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                        }
#endif
                        output_file.close();
                    }
                    Finalize();
                    std::exit(0);
#endif

                    //Finalize();
                    //std::exit(0);

                }

            }

            template <class Cell_t> 
            void redistMM(const MultientranceArray<S32> & rnklst,
                          const int n_cell_loc,
                          const std::vector<Cell_t> & cell_loc) {
                // Extract parameters
                const S32 p = param_.p;
                const S32vec n_cell = param_.n_cell;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                // Redistribute MM data 
                if (mode_ == 1) {
                    using pair_t = std::pair<S32, real_t>;
                    F64 time_start = GetWtime();
                    // Make a send buffer
                    CommBuffer<pair_t> sendbuf;
                    sendbuf.n_comm = n_proc_in_parent_group_;
                    sendbuf.allocCommInfo();
                    for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                        const S32 rnk =  (rank_in_parent_group_ + i) % n_proc_in_parent_group_;
                        sendbuf.ranks[i] = rnk;
                        sendbuf.adr_from_rank[rnk] = i;
                    }
                    sendbuf.clearCounts();
                    sendbuf.count_tot = 0;
                    for (S32 n = 0; n < n_cell_loc; n++) {
                        if (!cell_loc[n].is_mm_defined) continue;
                        const S32 idx = cell_loc[n].idx;
                        const S32 iz = idx/(n_cell.x * n_cell.y);
                        S32 slab_id;
                        for (S32 k = 0; k < n_proc_for_fft_; k++) {
                            if ((local_0_start_[k] <= iz) && (iz <= local_0_end_[k])) {
                                slab_id = k;
                                break;
                            }
                        }
                        for (S32 k = 0; k < n_group_; k++) {
                            const S32 rnk = rank_start_[k] + slab_id;
                            const S32 adr = sendbuf.adr_from_rank[rnk];
                            sendbuf.counts[adr]++;
                            sendbuf.count_tot++;
                        }
                    }
                    if (sendbuf.count_tot > 0) {
                        sendbuf.allocBuffer();
                        sendbuf.calcDispls();
                        sendbuf.clearCounts();
                        for (S32 n = 0; n < n_cell_loc; n++) {
                            if (!cell_loc[n].is_mm_defined) continue;
                            const S32 idx = cell_loc[n].idx;
                            const S32 iz = idx/(n_cell.x * n_cell.y);
                            S32 slab_id;
                            for (S32 k = 0; k < n_proc_for_fft_; k++) {
                                if ((local_0_start_[k] <= iz) && (iz <= local_0_end_[k])) {
                                    slab_id = k;
                                    break;
                                }
                            }
                            for (S32 k = 0; k < n_group_; k++) {
                                const S32 rnk = rank_start_[k] + slab_id;
                                const S32 adr = sendbuf.adr_from_rank[rnk];
                                const S32 adr_buf = sendbuf.displs[adr] 
                                                  + sendbuf.counts[adr];
                                sendbuf.buf[adr_buf].first  = idx;
                                sendbuf.buf[adr_buf].second = cell_loc[n].mm.buf[k];
                                sendbuf.counts[adr]++;
                            }
                        }
                    }
                    // Make a receive buffer
                    CommBuffer<pair_t> recvbuf;
                    recvbuf.n_comm = n_proc_in_parent_group_;
                    recvbuf.allocCommInfo();
                    for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                        recvbuf.ranks[i] = i;
                        recvbuf.adr_from_rank[i] = i;
                    }
                    recvbuf.clearCounts();
                    recvbuf.count_tot = 0;
                    if (rank_in_my_group_ < n_proc_for_fft_) { // this process receive data
                        const S32 iz_start = local_0_start_[rank_in_my_group_];
                        const S32 iz_end   = local_0_end_[rank_in_my_group_];
                        for (S32 iz = iz_start; iz <= iz_end; iz++) {
                            if (iz < n_cell.z) {
                                for (S32 iy = 0; iy < n_cell.y; iy++) {
                                    for (S32 ix = 0; ix < n_cell.x; ix++) {
                                        const S32 idx = ix + n_cell.x * (iy + n_cell.y * iz);
                                        const S32 cnt = rnklst.counts[idx];
                                        if (cnt > 0) {
                                            const S32 adr = rnklst.displs[idx]; 
                                            const S32 rnk = rnklst.data[adr];
                                            recvbuf.counts[rnk]++;
                                            recvbuf.count_tot++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (recvbuf.count_tot > 0) {
                        recvbuf.allocBuffer();
                        recvbuf.calcDispls();
                    }
                    time_profile_.M2L_preproc_mm_r_comm += GetWtime() - time_start;
                    // Perform MPI comm.
                    time_start = GetWtime();
                    performComm(sendbuf, recvbuf);
                    time_profile_.M2L_redist_mm_r += GetWtime() - time_start;
                    // Copy from recvbuf to mm_r
                    time_start = GetWtime();
                    mm_r_.resize();
                    mm_r_.clear();
                    for (S32 i = 0; i < recvbuf.count_tot; i++) {
                        // (1) Convert recvbuf[i].idx to S32vec idx_3d.
                        const S32 idx_1d = recvbuf.buf[i].first; // idx
                        S32vec idx_3d;
                        idx_3d.z = idx_1d / (n_cell.x * n_cell.y);
                        idx_3d.y = (idx_1d - (n_cell.x * n_cell.y) * idx_3d.z) / n_cell.x;
                        idx_3d.x = idx_1d - (n_cell.x * n_cell.y) * idx_3d.z - n_cell.x * idx_3d.y; 
                        // (2) Convert idx_3d to the one-dimensional cell index
                        //     based on the size information of mm_r
                        const S32 i_loc = idx_3d.x - mm_r_.start_loc_[1];
                        const S32 j_loc = idx_3d.y - mm_r_.start_loc_[2];
                        const S32 k_loc = idx_3d.z - mm_r_.start_loc_[3];
                        const S32 idx = i_loc 
                                      + mm_r_.size_loc_[1] * (j_loc
                                      + mm_r_.size_loc_[2] * k_loc);
                        // Note that recvbuf[i].idx is calculated based on
                        // the size of PM mesh, n_cell. On the other hand,
                        // the cell index of mm_r[] is defined based on
                        // its size, which depends on the boundary condition.
                        // This is the reason why we have performed the index
                        // calculation above.
                        assert(0 <= idx && idx < mm_r_.size_loc_tot_);
                        mm_r_.buf_[idx] += recvbuf.buf[i].second;
                    }
                    time_profile_.M2L_postproc_mm_r_comm += GetWtime() - time_start;
                } else if (mode_ == 2) {
                    F64 time_start = GetWtime();
                    // Make a list of # of local PM cells
                    S32 * n_cell_loc_tbl = new S32[n_proc_in_parent_group_];
                    Comm::allGather(&n_cell_loc, 1, n_cell_loc_tbl); 
                    // Make a send buffer
                    CommBuffer<S32> idx_send;
                    CommBuffer<real_t> mm_send;
                    idx_send.n_comm = n_proc_in_parent_group_;
                    mm_send.n_comm  = n_proc_in_parent_group_;
                    idx_send.allocCommInfo();
                    mm_send.allocCommInfo();
                    idx_send.count_tot = 0;
                    mm_send.count_tot = 0;
                    for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                        const S32 n_lm = lm_end_[i] - lm_start_[i] + 1;
                        idx_send.ranks[i] = i;
                        idx_send.counts[i] = n_cell_loc;
                        idx_send.count_tot += n_cell_loc;
                        mm_send.ranks[i] = i;
                        mm_send.counts[i] = n_lm * n_cell_loc;
                        mm_send.count_tot += n_lm * n_cell_loc;
                    }
                    if (idx_send.count_tot > 0) {
                        idx_send.allocBuffer();
                        mm_send.allocBuffer();
                        idx_send.calcDispls();
                        mm_send.calcDispls();
                        idx_send.clearCounts();
                        mm_send.clearCounts();
                        for (S32 n = 0; n < n_cell_loc; n++) {
                            const S32 idx = cell_loc[n].idx;
                            for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                                // set idx_send
                                const S32 adr = idx_send.displs[i]
                                              + idx_send.counts[i];
                                idx_send.buf[adr] = idx;
                                idx_send.counts[i]++;
                                // set mm_send
                                for (S32 lm = lm_start_[i]; lm <= lm_end_[i]; lm++) {
                                    const S32 adr = mm_send.displs[i]
                                                  + mm_send.counts[i];
                                    mm_send.buf[adr] = cell_loc[n].mm.buf[lm];
                                    mm_send.counts[i]++;
                                }
                            }
                        }
                    }
                    // Make a receive buffer
                    CommBuffer<S32> idx_recv;
                    CommBuffer<real_t> mm_recv;
                    idx_recv.n_comm = n_proc_in_parent_group_;
                    mm_recv.n_comm  = n_proc_in_parent_group_;
                    idx_recv.allocCommInfo();
                    mm_recv.allocCommInfo();
                    idx_recv.count_tot = 0;
                    mm_recv.count_tot = 0;
                    const S32 n_lm = lm_end_[rank_in_parent_group_] 
                                   - lm_start_[rank_in_parent_group_]
                                   + 1;
                    for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                        idx_recv.ranks[i] = i;
                        idx_recv.counts[i] = n_cell_loc_tbl[i];
                        idx_recv.count_tot += n_cell_loc_tbl[i];
                        mm_recv.ranks[i] = i;
                        mm_recv.counts[i] = n_lm * n_cell_loc_tbl[i];
                        mm_recv.count_tot += n_lm * n_cell_loc_tbl[i];
                    }
                    if (idx_recv.count_tot > 0) {
                        idx_recv.allocBuffer();
                        mm_recv.allocBuffer();
                        idx_recv.calcDispls();
                        mm_recv.calcDispls();
                    }
                    time_profile_.M2L_preproc_mm_r_comm += GetWtime() - time_start; 
                    // Perform MPI comm. 
                    time_start = GetWtime();
                    performComm(idx_send, idx_recv);
                    performComm(mm_send, mm_recv);
                    time_profile_.M2L_redist_mm_r += GetWtime() - time_start;
                    // Copy from recvbuf to mm_r
                    time_start = GetWtime();
                    mm_r_.resize();
                    mm_r_.clear();
                    for (S32 n = 0; n < idx_recv.n_comm; n++) {
                        const S32 rnk = idx_recv.ranks[n];
                        const S32 disp_for_idx = idx_recv.displs[n];
                        const S32 cnt = idx_recv.counts[n];
                        const S32 disp_for_mm = mm_recv.displs[n];
                        S32 k = 0; 
                        for (S32 i = 0; i < cnt; i++) {
                            // (1) Convert idx_recv.buf[] to S32vec idx_3d.
                            const S32 idx_1d = idx_recv.buf[disp_for_idx + i];
                            S32vec idx_3d;
                            idx_3d.z = idx_1d / (n_cell.x * n_cell.y);
                            idx_3d.y = (idx_1d - (n_cell.x * n_cell.y) * idx_3d.z) / n_cell.x;
                            idx_3d.x = idx_1d - (n_cell.x * n_cell.y) * idx_3d.z - n_cell.x * idx_3d.y; 
                            // (2) Copy
                            const S32 beg = lm_start_[rank_in_parent_group_];
                            const S32 end = lm_end_[rank_in_parent_group_];
                            for (S32 lm = beg; lm <= end; lm++) {
                                // (2-1) convert idx_3d to the one-dimensional cell index
                                //       based on the size information of mm_r
                                const S32 lm_loc = lm       - mm_r_.start_loc_[0];
                                const S32 i_loc  = idx_3d.x - mm_r_.start_loc_[1];
                                const S32 j_loc  = idx_3d.y - mm_r_.start_loc_[2];
                                const S32 k_loc  = idx_3d.z - mm_r_.start_loc_[3];
                                const S32 idx = lm_loc
                                              + mm_r_.size_loc_[0] * (i_loc 
                                              + mm_r_.size_loc_[1] * (j_loc
                                              + mm_r_.size_loc_[2] * k_loc));
                                // Note that idx_recv.buf[i] is calculated based on
                                // the size of PM mesh, n_cell. On the other hand,
                                // the cell index of mm_r[] is defined based on
                                // its size, which depends on the boundary condition.
                                // This is the reason why we have performed the index
                                // calculation above.
                                assert(0 <= idx && idx < mm_r_.size_loc_tot_);
                                // (2-2) copy mm_recv.buf[] to mm_r.buf_[]
                                mm_r_.buf_[idx] = mm_recv.buf[disp_for_mm + k];
                                k++;
                            }
                        }
                    }
                    time_profile_.M2L_postproc_mm_r_comm += GetWtime() - time_start;
                }
#else // PARTICLE_SIMULATOR_MPI_PARALLEL
                // Copy from cell_loc to mm_r_
                F64 time_start = GetWtime();
                mm_r_.resize();
                mm_r_.clear();
                for (S32 i = 0; i < n_cell_loc; i++) {
                    if (!cell_loc[i].is_mm_defined) continue;
                    const S32 idx_1d = cell_loc[i].idx;
                    S32vec idx_3d;
                    idx_3d.z = idx_1d / (n_cell.x * n_cell.y);
                    idx_3d.y = (idx_1d - (n_cell.x * n_cell.y) * idx_3d.z) / n_cell.x;
                    idx_3d.x = idx_1d - (n_cell.x * n_cell.y) * idx_3d.z - n_cell.x * idx_3d.y; 
                    const S32 i_loc = idx_3d.x - mm_r_.start_loc_[1];
                    const S32 j_loc = idx_3d.y - mm_r_.start_loc_[2];
                    const S32 k_loc = idx_3d.z - mm_r_.start_loc_[3];
                    for (S32 lm = 0; lm < cell_loc[i].mm.buf.size(); lm++) {
                        const S32 lm_loc = lm - mm_r_.start_loc_[0];
                        const S32 idx = lm_loc
                                      + mm_r_.size_loc_[0] * (i_loc 
                                      + mm_r_.size_loc_[1] * (j_loc
                                      + mm_r_.size_loc_[2] * k_loc));
                        assert(0 <= idx && idx < mm_r_.size_loc_tot_);
                        mm_r_.buf_[idx] += cell_loc[i].mm.buf[lm];
                    }
                }
                time_profile_.M2L_redist_mm_r += GetWtime() - time_start;
#endif // PARTICLE_SIMULATOR_MPI_PARALLEL

#if 0
                // Check mm_r 
                std::stringstream ss;
                ss << "mm_r_" << std::setfill('0') << std::setw(5)
                   << rank_in_parent_group_ << ".txt";
                const std::string filename = ss.str();
                std::ofstream output_file;
                output_file.open(filename.c_str(), std::ios::trunc);
                for (S32 k = mm_r_.start_loc_used_[3]; k <= mm_r_.end_loc_used_[3]; k++)
                    for (S32 j = mm_r_.start_loc_used_[2]; j <= mm_r_.end_loc_used_[2]; j++)
                        for (S32 i = mm_r_.start_loc_used_[1]; i <= mm_r_.end_loc_used_[1]; i++)
                            for (S32 lm = mm_r_.start_loc_used_[0]; lm <= mm_r_.end_loc_used_[0]; lm++)
                {
                    const S32 idx = lm 
                                  + mm_r_.size_glb_[0] * (i
                                  + mm_r_.size_glb_[1] * (j
                                  + mm_r_.size_glb_[2] * k));
                    const S32 lm_loc = lm - mm_r_.start_loc_[0];
                    const S32 i_loc  = i  - mm_r_.start_loc_[1];
                    const S32 j_loc  = j  - mm_r_.start_loc_[2];
                    const S32 k_loc  = k  - mm_r_.start_loc_[3];
                    const S32 adr = lm_loc
                                  + mm_r_.size_loc_[0] * (i_loc
                                  + mm_r_.size_loc_[1] * (j_loc
                                  + mm_r_.size_loc_[2] * k_loc));
                    if (mm_r_.buf_[adr] != 0) {
                        output_file << idx << "    " << mm_r_.buf_[adr] << std::endl;
                    }
                }
                output_file.close();
                Finalize();
                std::exit(0);
#endif

            }


            void convolution() { 
                // Extract parameters
                const S32 p = param_.p;
                const S32vec n_cell = param_.n_cell;
                const S32 LEN = (p+1)*(p+1);
                // Save the previous FFT configuration
                fft_mode_prev_ = fft_mode_;
                for (S32 i = 0; i < 3; i++) size_fft_prev_[i] = size_fft_[i];
                // Resize buffers
                mm_k_.resize();
                le_k_.resize();
                le_r_.resize();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL               
                if (mode_ == 1 && use_mpifft_) {
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        F64 time_start = GetWtime();
                        // Check if conditions that need to be fulfilled in this IF branch
                        // are actually fulfilled.
                        assert(mm_r_.size_loc_[0] == 1);
                        assert(le_k_.size_loc_[0] == 1);
                        assert(le_r_.size_loc_[0] == 1);
                        assert(mm_r_.start_loc_used_[0] == mm_k_.start_loc_used_[0]);
                        assert(mm_r_.end_loc_used_[0] == mm_k_.end_loc_used_[0]);

                        // Set FFT configuration
                        fft_mode_ = 1;
                        size_fft_[0] = mm_r_.size_glb_[1];
                        size_fft_[1] = mm_r_.size_glb_[2];
                        size_fft_[2] = mm_r_.size_glb_[3];

                        if ((fft_mode_ != fft_mode_prev_) ||
                            (size_fft_[0] != size_fft_prev_[0]) ||
                            (size_fft_[1] != size_fft_prev_[1]) ||
                            (size_fft_[2] != size_fft_prev_[2])) {
                            // Allocate buffer memory
                            ptrdiff_t alloc_local;
                            ptrdiff_t local_n0, local_0_start;
                            ptrdiff_t local_n1, local_1_start;
                            alloc_local = fftw_mpi_local_size_3d_transposed(mm_k_.size_glb_[3],
                                                                            mm_k_.size_glb_[2],
                                                                            mm_k_.size_glb_[1],
                                                                            comm_fft_,
                                                                            &local_n0, &local_0_start,
                                                                            &local_n1, &local_1_start);
                            // Note that an user must indicate the logical size of FFT by
                            // the sizes of a COMPLEX array when using fftw_mpi_local_size_3d_transposed.
                            const ptrdiff_t size_rbuf = 2 * alloc_local;
                            const ptrdiff_t size_kbuf = alloc_local;
                            if (rbuf_ != nullptr) fftw_free(rbuf_);
                            if (kbuf_ != nullptr) fftw_free(kbuf_);
                            rbuf_ = fftw_alloc_real(size_rbuf);
                            kbuf_ = fftw_alloc_complex(size_kbuf);

                            // Destroy plan if needed
                            for (S32 i = 0; i < plan_fwd_.size(); i++) fftw_destroy_plan(plan_fwd_[i]);
                            plan_fwd_.resize(0);
                            for (S32 i = 0; i < plan_bkw_.size(); i++) fftw_destroy_plan(plan_bkw_[i]);
                            plan_bkw_.resize(0);

                            // Create plan
                            fftw_plan tmp;
                            tmp = fftw_mpi_plan_dft_r2c_3d(
                                      mm_r_.size_glb_[3], mm_r_.size_glb_[2], mm_r_.size_glb_[1], 
                                      &rbuf_[0], &kbuf_[0], comm_fft_,
                                      FFTW_PLANNING_RIGOR_FLAG | FFTW_DESTROY_INPUT | FFTW_MPI_TRANSPOSED_OUT);
                            plan_fwd_.push_back(tmp);
                            tmp = fftw_mpi_plan_dft_c2r_3d(
                                      mm_r_.size_glb_[3], mm_r_.size_glb_[2], mm_r_.size_glb_[1], 
                                      &kbuf_[0], &rbuf_[0], comm_fft_,
                                      FFTW_PLANNING_RIGOR_FLAG | FFTW_DESTROY_INPUT | FFTW_MPI_TRANSPOSED_IN);
                            plan_bkw_.push_back(tmp);
                            // Note that fftw_mpi_plan_dft_r2c_3d requires an user to indicate
                            // the logical size of FFT by the size of a REAL array.
                        }

                        // Perform FFT (forward multipole moments)
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for 
#endif
                        for (S32 k = mm_r_.start_loc_used_[3]; k <= mm_r_.end_loc_used_[3]; k++)
                            for (S32 j = mm_r_.start_loc_used_[2]; j <= mm_r_.end_loc_used_[2]; j++)
                                for (S32 i = mm_r_.start_loc_used_[1]; i <= mm_r_.end_loc_used_[1]; i++)
                        {
                            const S32 i_loc = i - mm_r_.start_loc_[1];
                            const S32 j_loc = j - mm_r_.start_loc_[2];
                            const S32 k_loc = k - mm_r_.start_loc_[3];
                            const S32 adr_src  = i_loc
                                               + mm_r_.size_loc_[1] * (j_loc
                                               + mm_r_.size_loc_[2] * k_loc);
                            const S32 i_buf = i - mm_r_.start_loc_used_[1];
                            const S32 j_buf = j - mm_r_.start_loc_used_[2];
                            const S32 k_buf = k - mm_r_.start_loc_used_[3];
                            const S32 adr_dest = i_buf
                                               + (2*(mm_r_.size_glb_[1]/2 + 1)) * (j_buf
                                               + mm_r_.size_glb_[2] * k_buf);
                            // Note that the size of 1st dimension of rbuf is not NX.
                            // (see Section 6.5 of the mamual of FFTW3)
                            rbuf_[adr_dest] = mm_r_.buf_[adr_src];
                        }

                        // CALL FFTW
                        fftw_execute(plan_fwd_[0]);
    
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for 
#endif
                        for (S32 j = mm_k_.start_loc_used_[3]; j <= mm_k_.end_loc_used_[3]; j++)
                            for (S32 k = mm_k_.start_loc_used_[2]; k <= mm_k_.end_loc_used_[2]; k++)
                                for (S32 i = mm_k_.start_loc_used_[1]; i <= mm_k_.end_loc_used_[1]; i++)
                        {
                            const S32 i_buf = i - mm_k_.start_loc_used_[1];
                            const S32 k_buf = k - mm_k_.start_loc_used_[2];
                            const S32 j_buf = j - mm_k_.start_loc_used_[3];
                            const S32 adr_src = i_buf 
                                               + mm_k_.size_glb_[1] * (k_buf
                                               + mm_k_.size_glb_[2] * j_buf);
                            const S32 lm = mm_k_.start_loc_used_[0];
                            const S32 lm_loc = lm - mm_k_.start_loc_[0];
                            const S32 i_loc = i - mm_k_.start_loc_[1];
                            const S32 k_loc = k - mm_k_.start_loc_[2];
                            const S32 j_loc = j - mm_k_.start_loc_[3];
                            const S32 adr_dest = lm_loc
                                               + mm_k_.size_loc_[0] * (i_loc
                                               + mm_k_.size_loc_[1] * (k_loc
                                               + mm_k_.size_loc_[2] * j_loc));
                            // Note that j and k loops are exchanged because we specified
                            // FFTW_MPI_TRANSPOSED_OUT.

#ifdef PARTICLE_SIMULATOR_STD_COMPLEX_NOT_HAVING_SETTER
                            const cplx_t ctmp(kbuf_[adr_src][0],
                                              kbuf_[adr_src][1]);
                            mm_k_.buf_[adr_dest] = ctmp;
                            // [Notes (tag: #64d4cd48)]
                            //     Fujitsu C++ compiler installed in K computer
                            //     does not support std::complex<T>::real(T value)
                            //     and std::complex<T>::imag(T value).
#else
                            mm_k_.buf_[adr_dest].real(kbuf_[adr_src][0]);
                            mm_k_.buf_[adr_dest].imag(kbuf_[adr_src][1]);
#endif
                            // Note that fftw_complex = double[2]
                        }
                        time_profile_.M2L_mm_r_to_mm_k += GetWtime() - time_start;

                        // Collect mm_k[]
                        time_start = GetWtime();
                        // Calculate recvcount
                        const S32 recvcount = mm_k_.size_loc_used_[1]
                                            * mm_k_.size_loc_used_[2]
                                            * mm_k_.size_loc_used_[3];
                        // Set sendbuf (actually recvbuf)
                        cplx_t *recvbuf = new cplx_t[recvcount * n_group_];
                        S32 adr_dest = recvcount * idx_to_my_group_;
                        for (S32 lm = mm_k_.start_loc_used_[0]; lm <= mm_k_.end_loc_used_[0]; lm++) 
                            for (S32 k = mm_k_.start_loc_used_[3]; k <= mm_k_.end_loc_used_[3]; k++)
                                for (S32 j = mm_k_.start_loc_used_[2]; j <= mm_k_.end_loc_used_[2]; j++)
                                    for (S32 i = mm_k_.start_loc_used_[1]; i <= mm_k_.end_loc_used_[1]; i++)
                        {
                            const S32 lm_loc = lm - mm_k_.start_loc_[0];
                            const S32 i_loc  = i  - mm_k_.start_loc_[1];
                            const S32 j_loc  = j  - mm_k_.start_loc_[2];
                            const S32 k_loc  = k  - mm_k_.start_loc_[3];
                            const S32 adr_src = lm_loc
                                              + mm_k_.size_loc_[0] * (i_loc
                                              + mm_k_.size_loc_[1] * (j_loc
                                              + mm_k_.size_loc_[2] * k_loc));
                            recvbuf[adr_dest++] = mm_k_.buf_[adr_src];
                        } // for (lm,k,j,i)
                        // Perform MPI comm.
                        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                                      recvbuf, recvcount, mpi_cplx_t_, comm_int_);
                        // Copy from recvbuf to mm_k
                        adr_dest = 0;
                        for (S32 k = mm_k_.start_loc_used_[3]; k <= mm_k_.end_loc_used_[3]; k++)
                            for (S32 j = mm_k_.start_loc_used_[2]; j <= mm_k_.end_loc_used_[2]; j++)
                                for (S32 i = mm_k_.start_loc_used_[1]; i <= mm_k_.end_loc_used_[1]; i++)
                                    for (S32 lm = mm_k_.start_glb_[0]; lm <= mm_k_.end_glb_[0]; lm++) 
                        {
                            const S32 lm_loc = lm - mm_k_.start_loc_[0];
                            const S32 i_loc  = i  - mm_k_.start_loc_[1];
                            const S32 j_loc  = j  - mm_k_.start_loc_[2];
                            const S32 k_loc  = k  - mm_k_.start_loc_[3];
                            const S32 adr_src = i_loc
                                              + mm_k_.size_loc_[1] * (j_loc
                                              + mm_k_.size_loc_[2] * (k_loc
                                              + mm_k_.size_loc_[3] * lm_loc));
                            mm_k_.buf_[adr_dest++] = recvbuf[adr_src];
                        } // for (k,j,i,lm)
                        // Free
                        delete [] recvbuf;
                        time_profile_.M2L_collect_mm_k += GetWtime() - time_start;

                        // M2L transformation
                        time_start = GetWtime();
                        transform();
                        time_profile_.M2L_transform += GetWtime() - time_start;

                        // Peform FFT (backward local expansion)
                        time_start = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for 
#endif
                        for (S32 j = le_k_.start_loc_used_[3]; j <= le_k_.end_loc_used_[3]; j++)
                            for (S32 k = le_k_.start_loc_used_[2]; k <= le_k_.end_loc_used_[2]; k++)
                                for (S32 i = le_k_.start_loc_used_[1]; i <= le_k_.end_loc_used_[1]; i++)
                        {
                            const S32 i_loc = i - le_k_.start_loc_[1];
                            const S32 k_loc = k - le_k_.start_loc_[2];
                            const S32 j_loc = j - le_k_.start_loc_[3];
                            const S32 adr_src  = i_loc
                                               + le_k_.size_loc_[1] * (k_loc
                                               + le_k_.size_loc_[2] * j_loc);
                            const S32 i_buf = i - le_k_.start_loc_used_[1];
                            const S32 k_buf = k - le_k_.start_loc_used_[2];
                            const S32 j_buf = j - le_k_.start_loc_used_[3];
                            const S32 adr_dest = i_buf 
                                               + le_k_.size_glb_[1] * (k_buf
                                               + le_k_.size_glb_[2] * j_buf);
                            kbuf_[adr_dest][0] = le_k_.buf_[adr_src].real();
                            kbuf_[adr_dest][1] = le_k_.buf_[adr_src].imag();
                        }
        
                        fftw_execute(plan_bkw_[0]);
        
                        const F64 norm = 1.0 / (le_r_.size_glb_[1] * le_r_.size_glb_[2] * le_r_.size_glb_[3]);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for 
#endif
                        for (S32 k = le_r_.start_loc_used_[3]; k <= le_r_.end_loc_used_[3]; k++)
                            for (S32 j = le_r_.start_loc_used_[2]; j <= le_r_.end_loc_used_[2]; j++)
                                for (S32 i = le_r_.start_loc_used_[1]; i <= le_r_.end_loc_used_[1]; i++)
                        {
                            const S32 i_buf = i - le_r_.start_loc_used_[1];
                            const S32 j_buf = j - le_r_.start_loc_used_[2];
                            const S32 k_buf = k - le_r_.start_loc_used_[3];
                            const S32 adr_src  = i_buf
                                               + (2*(le_r_.size_glb_[1]/2 + 1)) * (j_buf
                                               + le_r_.size_glb_[2] * k_buf);
                            const S32 i_loc = i - le_r_.start_loc_[1];
                            const S32 j_loc = j - le_r_.start_loc_[2];
                            const S32 k_loc = k - le_r_.start_loc_[3];
                            const S32 adr_dest = i_loc
                                               + le_r_.size_loc_[1] * (j_loc
                                               + le_r_.size_loc_[2] * k_loc);
                            le_r_.buf_[adr_dest] = norm * rbuf_[adr_src];
                        }
                        time_profile_.M2L_le_k_to_le_r += GetWtime() - time_start;
                    }
                } else {
#endif
                    // In this case, we can use only a single MPI process to calculate
                    // mm_k for a particular set of (l,m).
                    // We use OpenMP parallerization of FFT if available.
                    if (rank_in_my_group_ == 0) {

                        // Check if conditions that need to be fulfilled in this IF branch
                        // are actually fulfilled.
                        assert(mm_r_.start_loc_used_[0] == mm_k_.start_loc_used_[0]);
                        assert(mm_r_.end_loc_used_[0] == mm_k_.end_loc_used_[0]);
                        for (S32 i = 1; i <= 3; i++) {
                            assert(mm_r_.size_loc_[i] == mm_r_.size_glb_[i]);
                            assert(mm_k_.size_loc_[i] == mm_k_.size_glb_[i]);
                        }
    
                        // Calculate mm_k
                        const S32 n_thread = Comm::getNumberOfThread();
                        const S32 fft_size = mm_r_.size_glb_[3] * mm_r_.size_glb_[2] * mm_r_.size_glb_[1];
                        if (fft_size < fft_size_crit_ && n_thread > 1) {
                            // In this case, the size of array is too small to speed up by multithreaded FFT.
                            // Hence, we assign a single FFT to each thread.
                            F64 time_start = GetWtime(); 
                            // Set FFT configuration
                            fft_mode_ = 2;
                            size_fft_[0] = mm_r_.size_glb_[1];
                            size_fft_[1] = mm_r_.size_glb_[2];
                            size_fft_[2] = mm_r_.size_glb_[3];

                            if ((fft_mode_ != fft_mode_prev_) ||
                                (size_fft_[0] != size_fft_prev_[0]) ||
                                (size_fft_[1] != size_fft_prev_[1]) ||
                                (size_fft_[2] != size_fft_prev_[2])) {
 
                                // Memory allocation
                                const S32 size_rbuf = n_thread * fft_size;
                                const S32 size_kbuf = n_thread * fft_size;
                                if (rbuf_ != nullptr) fftw_free(rbuf_);
                                if (kbuf_ != nullptr) fftw_free(kbuf_);
                                rbuf_ = fftw_alloc_real(size_rbuf);
                                kbuf_ = fftw_alloc_complex(size_kbuf);

                                // Destroy plan if needed
                                for (S32 i = 0; i < plan_fwd_.size(); i++) fftw_destroy_plan(plan_fwd_[i]);
                                plan_fwd_.resize(n_thread);
                                for (S32 i = 0; i < plan_bkw_.size(); i++) fftw_destroy_plan(plan_bkw_[i]);
                                plan_bkw_.resize(n_thread);
    
                                // Create plans of FFTW
                                for (S32 i = 0; i < n_thread; i++) {
                                    const S32 offset = fft_size * i;
                                    plan_fwd_[i] = fftw_plan_dft_r2c_3d(mm_r_.size_glb_[3],
                                                                        mm_r_.size_glb_[2],
                                                                        mm_r_.size_glb_[1], 
                                                                        &rbuf_[offset],
                                                                        &kbuf_[offset],
                                                                        FFTW_PLANNING_RIGOR_FLAG | FFTW_DESTROY_INPUT);
                                    plan_bkw_[i] = fftw_plan_dft_c2r_3d(mm_r_.size_glb_[3],
                                                                        mm_r_.size_glb_[2],
                                                                        mm_r_.size_glb_[1], 
                                                                        &kbuf_[offset],
                                                                        &rbuf_[offset],
                                                                        FFTW_PLANNING_RIGOR_FLAG | FFTW_DESTROY_INPUT);
                                }
                            }
    
                            // Perform FFT (forward multipole)
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for 
#endif
                            for (S32 lm = mm_r_.start_loc_used_[0]; lm <= mm_r_.end_loc_used_[0]; lm++){
                                const S32 ith = Comm::getThreadNum();
                                const S32 offset = fft_size * ith;
    
                                S32 lm_loc = lm - mm_r_.start_loc_[0];
                                for (S32 k = mm_r_.start_loc_used_[3]; k <= mm_r_.end_loc_used_[3]; k++)
                                    for (S32 j = mm_r_.start_loc_used_[2]; j <= mm_r_.end_loc_used_[2]; j++)
                                        for (S32 i = mm_r_.start_loc_used_[1]; i <= mm_r_.end_loc_used_[1]; i++)
                                {
                                    const S32 i_loc = i - mm_r_.start_loc_[1];
                                    const S32 j_loc = j - mm_r_.start_loc_[2];
                                    const S32 k_loc = k - mm_r_.start_loc_[3];
                                    const S32 adr_src  = lm_loc
                                                       + mm_r_.size_loc_[0] * (i_loc
                                                       + mm_r_.size_loc_[1] * (j_loc
                                                       + mm_r_.size_loc_[2] * k_loc));
                                    const S32 adr_dest = i 
                                                       + mm_r_.size_glb_[1] * (j
                                                       + mm_r_.size_glb_[2] * k);
                                    rbuf_[adr_dest + offset] = mm_r_.buf_[adr_src];
                                }
                
                                fftw_execute(plan_fwd_[ith]);
                
                                lm_loc = lm - mm_k_.start_loc_[0];
                                for (S32 k = mm_k_.start_loc_used_[3]; k <= mm_k_.end_loc_used_[3]; k++)
                                    for (S32 j = mm_k_.start_loc_used_[2]; j <= mm_k_.end_loc_used_[2]; j++)
                                        for (S32 i = mm_k_.start_loc_used_[1]; i <= mm_k_.end_loc_used_[1]; i++)
                                {
                                    const S32 i_loc = i - mm_k_.start_loc_[1];
                                    const S32 j_loc = j - mm_k_.start_loc_[2];
                                    const S32 k_loc = k - mm_k_.start_loc_[3];
                                    const S32 adr_src  = i 
                                                       + mm_k_.size_glb_[1] * (j
                                                       + mm_k_.size_glb_[2] * k);
                                    const S32 adr_dest = lm_loc
                                                       + mm_k_.size_loc_[0] * (i_loc
                                                       + mm_k_.size_loc_[1] * (j_loc
                                                       + mm_k_.size_loc_[2] * k_loc));
#ifdef PARTICLE_SIMULATOR_STD_COMPLEX_NOT_HAVING_SETTER
                                    const cplx_t ctmp(kbuf_[adr_src + offset][0],
                                                      kbuf_[adr_src + offset][1]);
                                    mm_k_.buf_[adr_dest] = ctmp;
#else
                                    mm_k_.buf_[adr_dest].real(kbuf_[adr_src + offset][0]);
                                    mm_k_.buf_[adr_dest].imag(kbuf_[adr_src + offset][1]);
#endif
                                }
                            }
                            time_profile_.M2L_mm_r_to_mm_k += GetWtime() - time_start;
    
                            // Collect mm_k[]
                            time_start = GetWtime();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                            if (mode_ == 1) {
                                if (comm_int_ != MPI_COMM_NULL) {
                                    const S32 recvcount = mm_k_.size_loc_used_[1]
                                                        * mm_k_.size_loc_used_[2]
                                                        * mm_k_.size_loc_used_[3];
                                    cplx_t *recvbuf = new cplx_t[recvcount * n_group_];
                                    // Set sendbuf (actually recvbuf)
                                    S32 adr_dest = recvcount * idx_to_my_group_;
                                    for (S32 lm = mm_k_.start_loc_used_[0]; lm <= mm_k_.end_loc_used_[0]; lm++) 
                                        for (S32 k = mm_k_.start_loc_used_[3]; k <= mm_k_.end_loc_used_[3]; k++)
                                            for (S32 j = mm_k_.start_loc_used_[2]; j <= mm_k_.end_loc_used_[2]; j++)
                                                for (S32 i = mm_k_.start_loc_used_[1]; i <= mm_k_.end_loc_used_[1]; i++)
                                    {
                                        const S32 lm_loc = lm - mm_k_.start_loc_[0];
                                        const S32 i_loc  = i  - mm_k_.start_loc_[1];
                                        const S32 j_loc  = j  - mm_k_.start_loc_[2];
                                        const S32 k_loc  = k  - mm_k_.start_loc_[3];
                                        const S32 adr_src = lm_loc
                                                          + mm_k_.size_loc_[0] * (i_loc
                                                          + mm_k_.size_loc_[1] * (j_loc
                                                          + mm_k_.size_loc_[2] * k_loc));
                                        recvbuf[adr_dest++] = mm_k_.buf_[adr_src];
                                    } // for (lm,k,j,i)
                                    // Perform MPI comm.
                                    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                                                  recvbuf, recvcount, mpi_cplx_t_, comm_int_);
                                    // Copy from recvbuf to mm_k
                                    adr_dest = 0;
                                    for (S32 k = mm_k_.start_loc_used_[3]; k <= mm_k_.end_loc_used_[3]; k++)
                                        for (S32 j = mm_k_.start_loc_used_[2]; j <= mm_k_.end_loc_used_[2]; j++)
                                            for (S32 i = mm_k_.start_loc_used_[1]; i <= mm_k_.end_loc_used_[1]; i++)
                                                for (S32 lm = mm_k_.start_glb_[0]; lm <= mm_k_.end_glb_[0]; lm++) 
                                    {
                                        const S32 lm_loc = lm - mm_k_.start_loc_[0];
                                        const S32 i_loc  = i  - mm_k_.start_loc_[1];
                                        const S32 j_loc  = j  - mm_k_.start_loc_[2];
                                        const S32 k_loc  = k  - mm_k_.start_loc_[3];
                                        const S32 adr_src = i_loc
                                                          + mm_k_.size_loc_[1] * (j_loc
                                                          + mm_k_.size_loc_[2] * (k_loc
                                                          + mm_k_.size_loc_[3] * lm_loc));
                                        mm_k_.buf_[adr_dest++] = recvbuf[adr_src];
                                    } // for (k,j,i,lm)
                                    // Free
                                    delete [] recvbuf;
                                }
                            } else if (mode_ == 2) {
                                // Calculate recvcounts[]
                                S32 recvcount_tot = 0;
                                S32 *recvcounts = new S32[n_proc_in_parent_group_];
                                for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                                    recvcounts[i] = (lm_end_[i] - lm_start_[i] + 1)
                                                  * mm_k_.size_loc_used_[1]
                                                  * mm_k_.size_loc_used_[2]
                                                  * mm_k_.size_loc_used_[3];
                                    recvcount_tot += recvcounts[i];
                                }
                                // Calculate recvdispls[]
                                S32 *recvdispls = new S32[n_proc_in_parent_group_];
                                recvdispls[0] = 0;
                                for (S32 i = 1; i < n_proc_in_parent_group_; i++) 
                                    recvdispls[i] = recvdispls[i-1] + recvcounts[i-1];
                                // Set sendbuf (actually recvbuf)
                                cplx_t *recvbuf = new cplx_t[recvcount_tot];
                                S32 adr_dest = recvdispls[rank_in_parent_group_];
                                for (S32 lm = mm_k_.start_loc_used_[0]; lm <= mm_k_.end_loc_used_[0]; lm++) 
                                    for (S32 k = mm_k_.start_loc_used_[3]; k <= mm_k_.end_loc_used_[3]; k++)
                                        for (S32 j = mm_k_.start_loc_used_[2]; j <= mm_k_.end_loc_used_[2]; j++)
                                            for (S32 i = mm_k_.start_loc_used_[1]; i <= mm_k_.end_loc_used_[1]; i++)
                                {
                                    const S32 lm_loc = lm - mm_k_.start_loc_[0];
                                    const S32 i_loc  = i  - mm_k_.start_loc_[1];
                                    const S32 j_loc  = j  - mm_k_.start_loc_[2];
                                    const S32 k_loc  = k  - mm_k_.start_loc_[3];
                                    const S32 adr_src = lm_loc
                                                      + mm_k_.size_loc_[0] * (i_loc
                                                      + mm_k_.size_loc_[1] * (j_loc
                                                      + mm_k_.size_loc_[2] * k_loc));
                                    recvbuf[adr_dest++] = mm_k_.buf_[adr_src];
                                } // for (lm,k,j,i)
                                // Perform MPI comm.                            
                                MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                                               recvbuf, recvcounts, recvdispls,
                                               mpi_cplx_t_, parent_comm_);
                                // Copy from recvbuf to mm_k
                                adr_dest = 0;
                                for (S32 k = mm_k_.start_loc_used_[3]; k <= mm_k_.end_loc_used_[3]; k++)
                                    for (S32 j = mm_k_.start_loc_used_[2]; j <= mm_k_.end_loc_used_[2]; j++)
                                        for (S32 i = mm_k_.start_loc_used_[1]; i <= mm_k_.end_loc_used_[1]; i++)
                                            for (S32 lm = mm_k_.start_glb_[0]; lm <= mm_k_.end_glb_[0]; lm++) 
                                {
                                    const S32 lm_loc = lm - mm_k_.start_loc_[0];
                                    const S32 i_loc  = i  - mm_k_.start_loc_[1];
                                    const S32 j_loc  = j  - mm_k_.start_loc_[2];
                                    const S32 k_loc  = k  - mm_k_.start_loc_[3];
                                    const S32 adr_src = i_loc
                                                      + mm_k_.size_loc_[1] * (j_loc
                                                      + mm_k_.size_loc_[2] * (k_loc
                                                      + mm_k_.size_loc_[3] * lm_loc));
                                    mm_k_.buf_[adr_dest++] = recvbuf[adr_src];
                                } // for (k,j,i,lm)
                                // Free
                                delete [] recvcounts;
                                delete [] recvdispls;
                                delete [] recvbuf;
                            } else {
                                assert(false);
                            }
#endif
                            time_profile_.M2L_collect_mm_k += GetWtime() - time_start;

                            // M2L transformation
                            time_start = GetWtime();
                            transform();
                            time_profile_.M2L_transform += GetWtime() - time_start;

                            // Peform FFT (backward local expansion)
                            time_start = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for 
#endif
                            for (S32 lm = le_r_.start_loc_used_[0]; lm <= le_r_.end_loc_used_[0]; lm++) {
                                const S32 ith = Comm::getThreadNum();
                                const S32 offset = fft_size * ith;
    
                                const S32 lm_loc = lm - le_r_.start_loc_[0];
                                for (S32 k = le_k_.start_loc_used_[3]; k <= le_k_.end_loc_used_[3]; k++)
                                    for (S32 j = le_k_.start_loc_used_[2]; j <= le_k_.end_loc_used_[2]; j++)
                                        for (S32 i = le_k_.start_loc_used_[1]; i <= le_k_.end_loc_used_[1]; i++)
                                {
                                    const S32 i_loc = i - le_k_.start_loc_[1];
                                    const S32 j_loc = j - le_k_.start_loc_[2];
                                    const S32 k_loc = k - le_k_.start_loc_[3];
                                    const S32 adr_src  = lm_loc 
                                                       + le_k_.size_loc_[0] * (i_loc
                                                       + le_k_.size_loc_[1] * (j_loc
                                                       + le_k_.size_loc_[2] * k_loc));
                                    const S32 adr_dest = i 
                                                       + le_k_.size_glb_[1] * (j 
                                                       + le_k_.size_glb_[2] * k);
                                    kbuf_[adr_dest + offset][0] = le_k_.buf_[adr_src].real();
                                    kbuf_[adr_dest + offset][1] = le_k_.buf_[adr_src].imag();
                                }
            
                                fftw_execute(plan_bkw_[ith]);
            
                                const F64 norm = 1.0 / (le_r_.size_glb_[1] * le_r_.size_glb_[2] * le_r_.size_glb_[3]);
                                for (S32 k = le_r_.start_loc_used_[3]; k <= le_r_.end_loc_used_[3]; k++)
                                    for (S32 j = le_r_.start_loc_used_[2]; j <= le_r_.end_loc_used_[2]; j++)
                                        for (S32 i = le_r_.start_loc_used_[1]; i <= le_r_.end_loc_used_[1]; i++)
                                {
                                    const S32 i_loc = i - le_r_.start_loc_[1];
                                    const S32 j_loc = j - le_r_.start_loc_[2];
                                    const S32 k_loc = k - le_r_.start_loc_[3];
                                    const S32 adr_src  = i 
                                                       + le_r_.size_glb_[1] * (j
                                                       + le_r_.size_glb_[2] * k);
                                    const S32 adr_dest = lm_loc
                                                       + le_r_.size_loc_[0] * (i_loc
                                                       + le_r_.size_loc_[1] * (j_loc
                                                       + le_r_.size_loc_[2] * k_loc));
                                    le_r_.buf_[adr_dest] = norm * rbuf_[adr_src + offset];
                                }
                            }
                            time_profile_.M2L_le_k_to_le_r += GetWtime() - time_start; 
                        } else {
                            // In this case, we can expect that the FFT calculation become fast by using
                            // a multithreaded FFT.
                            F64 time_start = GetWtime();
                            // Set FFT configuration
                            fft_mode_ = 3;
                            size_fft_[0] = mm_r_.size_glb_[1];
                            size_fft_[1] = mm_r_.size_glb_[2];
                            size_fft_[2] = mm_r_.size_glb_[3];
    
                            if ((fft_mode_ != fft_mode_prev_) ||
                                (size_fft_[0] != size_fft_prev_[0]) ||
                                (size_fft_[1] != size_fft_prev_[1]) ||
                                (size_fft_[2] != size_fft_prev_[2])) {

                                // Memory allocation
                                const S32 size_rbuf = mm_r_.size_glb_[3] * mm_r_.size_glb_[2] * mm_r_.size_glb_[1];
                                const S32 size_kbuf = mm_k_.size_glb_[3] * mm_k_.size_glb_[2] * mm_k_.size_glb_[1];
                                if (rbuf_ != nullptr) fftw_free(rbuf_);
                                if (kbuf_ != nullptr) fftw_free(kbuf_);
                                rbuf_ = fftw_alloc_real(size_rbuf);
                                kbuf_ = fftw_alloc_complex(size_kbuf);
                
                                // Destroy plan if needed
                                for (S32 i = 0; i < plan_fwd_.size(); i++) fftw_destroy_plan(plan_fwd_[i]);
                                plan_fwd_.resize(0);
                                for (S32 i = 0; i < plan_bkw_.size(); i++) fftw_destroy_plan(plan_bkw_[i]);
                                plan_bkw_.resize(0);

                                // Create plan
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                                fftw_plan_with_nthreads(n_thread);
#endif
                                fftw_plan tmp;
                                tmp = fftw_plan_dft_r2c_3d(
                                        mm_r_.size_glb_[3], mm_r_.size_glb_[2], mm_r_.size_glb_[1], 
                                        &rbuf_[0], &kbuf_[0],
                                        FFTW_PLANNING_RIGOR_FLAG | FFTW_DESTROY_INPUT);
                                plan_fwd_.push_back(tmp);
                                tmp = fftw_plan_dft_c2r_3d(
                                        mm_r_.size_glb_[3], mm_r_.size_glb_[2], mm_r_.size_glb_[1], 
                                        &kbuf_[0], &rbuf_[0],
                                        FFTW_PLANNING_RIGOR_FLAG | FFTW_DESTROY_INPUT);
                                plan_bkw_.push_back(tmp);
                            }

                            // Clear
                            //for (S32 i = 0; i < size_rbuf; i++) rbuf[i] = 0;

                            // Perform FFT (forward multipole moments)
                            for (S32 lm = mm_r_.start_loc_used_[0]; lm <= mm_r_.end_loc_used_[0]; lm++) {
                                S32 lm_loc = lm - mm_r_.start_loc_[0];
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for 
#endif
                                for (S32 k = mm_r_.start_loc_used_[3]; k <= mm_r_.end_loc_used_[3]; k++)
                                    for (S32 j = mm_r_.start_loc_used_[2]; j <= mm_r_.end_loc_used_[2]; j++)
                                        for (S32 i = mm_r_.start_loc_used_[1]; i <= mm_r_.end_loc_used_[1]; i++)
                                {
                                    const S32 i_loc = i - mm_r_.start_loc_[1];
                                    const S32 j_loc = j - mm_r_.start_loc_[2];
                                    const S32 k_loc = k - mm_r_.start_loc_[3];
                                    const S32 adr_src  = lm_loc
                                                       + mm_r_.size_loc_[0] * (i_loc
                                                       + mm_r_.size_loc_[1] * (j_loc
                                                       + mm_r_.size_loc_[2] * k_loc));
                                    const S32 adr_dest = i 
                                                       + mm_r_.size_glb_[1] * (j
                                                       + mm_r_.size_glb_[2] * k);
                                    rbuf_[adr_dest] = mm_r_.buf_[adr_src];
                                }

                                // CALL FFTW
                                fftw_execute(plan_fwd_[0]);
    
                                lm_loc = lm - mm_k_.start_loc_[0];
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for 
#endif
                                for (S32 k = mm_k_.start_loc_used_[3]; k <= mm_k_.end_loc_used_[3]; k++)
                                    for (S32 j = mm_k_.start_loc_used_[2]; j <= mm_k_.end_loc_used_[2]; j++)
                                        for (S32 i = mm_k_.start_loc_used_[1]; i <= mm_k_.end_loc_used_[1]; i++)
                                {
                                    const S32 i_loc = i - mm_k_.start_loc_[1];
                                    const S32 j_loc = j - mm_k_.start_loc_[2];
                                    const S32 k_loc = k - mm_k_.start_loc_[3];
                                    const S32 adr_dest = lm_loc
                                                       + mm_k_.size_loc_[0] * (i_loc
                                                       + mm_k_.size_loc_[1] * (j_loc
                                                       + mm_k_.size_loc_[2] * k_loc));
                                    const S32 adr_src  = i 
                                                       + mm_k_.size_glb_[1] * (j 
                                                       + mm_k_.size_glb_[2] * k);
#ifdef PARTICLE_SIMULATOR_STD_COMPLEX_NOT_HAVING_SETTER
                                    const cplx_t ctmp (kbuf_[adr_src][0],
                                                       kbuf_[adr_src][1]);
                                    mm_k_.buf_[adr_dest] = ctmp;
#else
                                    mm_k_.buf_[adr_dest].real(kbuf_[adr_src][0]);
                                    mm_k_.buf_[adr_dest].imag(kbuf_[adr_src][1]);
#endif
                                    // Note that fftw_complex = double[2]
                                }
                            }
                            time_profile_.M2L_mm_r_to_mm_k += GetWtime() - time_start;

#if 0
                            // Check mm_k
                            std::stringstream ss;
                            ss << "mm_k_" << std::setfill('0') << std::setw(5)
                               << rank_in_parent_group_ << ".txt";
                            const std::string filename = ss.str();
                            std::ofstream output_file;
                            output_file.open(filename.c_str(), std::ios::trunc);
                            for (S32 k = mm_k_.start_loc_used_[3]; k <= mm_k_.end_loc_used_[3]; k++)
                                for (S32 j = mm_k_.start_loc_used_[2]; j <= mm_k_.end_loc_used_[2]; j++)
                                    for (S32 i = mm_k_.start_loc_used_[1]; i <= mm_k_.end_loc_used_[1]; i++)
                                        for (S32 lm = mm_k_.start_loc_used_[0]; lm <= mm_k_.end_loc_used_[0]; lm++)
                            {
                                const S32 idx = lm 
                                              + mm_k_.size_glb_[0] * (i
                                              + mm_k_.size_glb_[1] * (j
                                              + mm_k_.size_glb_[2] * k));
                                const S32 lm_loc = lm - mm_k_.start_loc_[0];
                                const S32 i_loc  = i  - mm_k_.start_loc_[1];
                                const S32 j_loc  = j  - mm_k_.start_loc_[2];
                                const S32 k_loc  = k  - mm_k_.start_loc_[3];
                                const S32 adr = lm_loc
                                              + mm_k_.size_loc_[0] * (i_loc
                                              + mm_k_.size_loc_[1] * (j_loc
                                              + mm_k_.size_loc_[2] * k_loc));
                                if (mm_k_.buf_[adr].real() != 0 ||
                                    mm_k_.buf_[adr].imag() != 0) {
                                    output_file << idx << "    " << mm_k_.buf_[adr] << std::endl;
                                }
                            }
                            output_file.close();
                            Finalize();
                            std::exit(0);
#endif


                            // Collect mm_k[]
                            time_start = GetWtime();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                            if (mode_ == 1) {
                                if (comm_int_ != MPI_COMM_NULL) {
                                    const S32 recvcount = mm_k_.size_loc_used_[1]
                                                        * mm_k_.size_loc_used_[2]
                                                        * mm_k_.size_loc_used_[3];
                                    // Set sendbuf (actually recvbuf)
                                    cplx_t *recvbuf = new cplx_t[recvcount * n_group_];
                                    S32 adr_dest = recvcount * idx_to_my_group_;
                                    for (S32 lm = mm_k_.start_loc_used_[0]; lm <= mm_k_.end_loc_used_[0]; lm++) 
                                        for (S32 k = mm_k_.start_loc_used_[3]; k <= mm_k_.end_loc_used_[3]; k++)
                                            for (S32 j = mm_k_.start_loc_used_[2]; j <= mm_k_.end_loc_used_[2]; j++)
                                                for (S32 i = mm_k_.start_loc_used_[1]; i <= mm_k_.end_loc_used_[1]; i++)
                                    {
                                        const S32 lm_loc = lm - mm_k_.start_loc_[0];
                                        const S32 i_loc  = i  - mm_k_.start_loc_[1];
                                        const S32 j_loc  = j  - mm_k_.start_loc_[2];
                                        const S32 k_loc  = k  - mm_k_.start_loc_[3];
                                        const S32 adr_src = lm_loc
                                                          + mm_k_.size_loc_[0] * (i_loc
                                                          + mm_k_.size_loc_[1] * (j_loc
                                                          + mm_k_.size_loc_[2] * k_loc));
                                        recvbuf[adr_dest++] = mm_k_.buf_[adr_src];
                                    } // for (lm,k,j,i)
                                    // Perform MPI comm.
                                    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                                                  recvbuf, recvcount, mpi_cplx_t_, comm_int_);
                                    // Copy from recvbuf to mm_k
                                    adr_dest = 0;
                                    for (S32 k = mm_k_.start_loc_used_[3]; k <= mm_k_.end_loc_used_[3]; k++)
                                        for (S32 j = mm_k_.start_loc_used_[2]; j <= mm_k_.end_loc_used_[2]; j++)
                                            for (S32 i = mm_k_.start_loc_used_[1]; i <= mm_k_.end_loc_used_[1]; i++)
                                                for (S32 lm = mm_k_.start_glb_[0]; lm <= mm_k_.end_glb_[0]; lm++) 
                                    {
                                        const S32 lm_loc = lm - mm_k_.start_loc_[0];
                                        const S32 i_loc  = i  - mm_k_.start_loc_[1];
                                        const S32 j_loc  = j  - mm_k_.start_loc_[2];
                                        const S32 k_loc  = k  - mm_k_.start_loc_[3];
                                        const S32 adr_src = i_loc
                                                          + mm_k_.size_loc_[1] * (j_loc
                                                          + mm_k_.size_loc_[2] * (k_loc
                                                          + mm_k_.size_loc_[3] * lm_loc));
                                        mm_k_.buf_[adr_dest++] = recvbuf[adr_src];
                                    } // for (k,j,i,lm)
                                    // Free
                                    delete [] recvbuf;
                                }
                            } else if (mode_ == 2) {
                                // Calculate recvcounts[]
                                S32 recvcount_tot = 0;
                                S32 *recvcounts = new S32[n_proc_in_parent_group_];
                                for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                                    recvcounts[i] = (lm_end_[i] - lm_start_[i] + 1)
                                                  * mm_k_.size_loc_used_[1]
                                                  * mm_k_.size_loc_used_[2]
                                                  * mm_k_.size_loc_used_[3];
                                    recvcount_tot += recvcounts[i];
                                }
                                // Calculate recvdispls[]
                                S32 *recvdispls = new S32[n_proc_in_parent_group_];
                                recvdispls[0] = 0;
                                for (S32 i = 1; i < n_proc_in_parent_group_; i++) 
                                    recvdispls[i] = recvdispls[i-1] + recvcounts[i-1];
                                // Set sendbuf (actually recvbuf)
                                cplx_t *recvbuf = new cplx_t[recvcount_tot];
                                S32 adr_dest = recvdispls[rank_in_parent_group_];
                                for (S32 lm = mm_k_.start_loc_used_[0]; lm <= mm_k_.end_loc_used_[0]; lm++) 
                                    for (S32 k = mm_k_.start_loc_used_[3]; k <= mm_k_.end_loc_used_[3]; k++)
                                        for (S32 j = mm_k_.start_loc_used_[2]; j <= mm_k_.end_loc_used_[2]; j++)
                                            for (S32 i = mm_k_.start_loc_used_[1]; i <= mm_k_.end_loc_used_[1]; i++)
                                {
                                    const S32 lm_loc = lm - mm_k_.start_loc_[0];
                                    const S32 i_loc  = i  - mm_k_.start_loc_[1];
                                    const S32 j_loc  = j  - mm_k_.start_loc_[2];
                                    const S32 k_loc  = k  - mm_k_.start_loc_[3];
                                    const S32 adr_src = lm_loc
                                                      + mm_k_.size_loc_[0] * (i_loc
                                                      + mm_k_.size_loc_[1] * (j_loc
                                                      + mm_k_.size_loc_[2] * k_loc));
                                    recvbuf[adr_dest++] = mm_k_.buf_[adr_src];
                                } // for (lm,k,j,i)
                                // Perform MPI comm.                            
                                MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                                               recvbuf, recvcounts, recvdispls,
                                               mpi_cplx_t_, parent_comm_);
                                // Copy from recvbuf to mm_k
                                adr_dest = 0;
                                for (S32 k = mm_k_.start_loc_used_[3]; k <= mm_k_.end_loc_used_[3]; k++)
                                    for (S32 j = mm_k_.start_loc_used_[2]; j <= mm_k_.end_loc_used_[2]; j++)
                                        for (S32 i = mm_k_.start_loc_used_[1]; i <= mm_k_.end_loc_used_[1]; i++)
                                            for (S32 lm = mm_k_.start_glb_[0]; lm <= mm_k_.end_glb_[0]; lm++) 
                                {
                                    const S32 lm_loc = lm - mm_k_.start_loc_[0];
                                    const S32 i_loc  = i  - mm_k_.start_loc_[1];
                                    const S32 j_loc  = j  - mm_k_.start_loc_[2];
                                    const S32 k_loc  = k  - mm_k_.start_loc_[3];
                                    const S32 adr_src = i_loc
                                                      + mm_k_.size_loc_[1] * (j_loc
                                                      + mm_k_.size_loc_[2] * (k_loc
                                                      + mm_k_.size_loc_[3] * lm_loc));
                                    mm_k_.buf_[adr_dest++] = recvbuf[adr_src];
                                } // for (k,j,i,lm)
                                // Free
                                delete [] recvcounts;
                                delete [] recvdispls;
                                delete [] recvbuf;
                            } else {
                                assert(false);
                            }
#endif
                            time_profile_.M2L_collect_mm_k += GetWtime() - time_start;

                            // M2L transformation
                            time_start = GetWtime();
                            transform();
                            time_profile_.M2L_transform += GetWtime() - time_start;

                            // Peform FFT (backward local expansion)
                            time_start = GetWtime();
                            for (S32 lm = le_r_.start_loc_used_[0]; lm <= le_r_.end_loc_used_[0]; lm++) {
                                const S32 lm_loc = lm - le_r_.start_loc_[0];
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for 
#endif
                                for (S32 k = le_k_.start_loc_used_[3]; k <= le_k_.end_loc_used_[3]; k++)
                                    for (S32 j = le_k_.start_loc_used_[2]; j <= le_k_.end_loc_used_[2]; j++)
                                        for (S32 i = le_k_.start_loc_used_[1]; i <= le_k_.end_loc_used_[1]; i++)
                                {
                                    const S32 i_loc = i - le_k_.start_loc_[1];
                                    const S32 j_loc = j - le_k_.start_loc_[2];
                                    const S32 k_loc = k - le_k_.start_loc_[3];
                                    const S32 adr_src  = lm_loc 
                                                       + le_k_.size_loc_[0] * (i_loc
                                                       + le_k_.size_loc_[1] * (j_loc
                                                       + le_k_.size_loc_[2] * k_loc));
                                    const S32 adr_dest = i 
                                                       + le_k_.size_glb_[1] * (j 
                                                       + le_k_.size_glb_[2] * k);
                                    kbuf_[adr_dest][0] = le_k_.buf_[adr_src].real();
                                    kbuf_[adr_dest][1] = le_k_.buf_[adr_src].imag();
                                }
            
                                fftw_execute(plan_bkw_[0]);
            
                                const F64 norm = 1.0 / (le_r_.size_glb_[1] * le_r_.size_glb_[2] * le_r_.size_glb_[3]);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for 
#endif
                                for (S32 k = le_r_.start_loc_used_[3]; k <= le_r_.end_loc_used_[3]; k++)
                                    for (S32 j = le_r_.start_loc_used_[2]; j <= le_r_.end_loc_used_[2]; j++)
                                        for (S32 i = le_r_.start_loc_used_[1]; i <= le_r_.end_loc_used_[1]; i++)
                                {
                                    const S32 i_loc = i - le_r_.start_loc_[1];
                                    const S32 j_loc = j - le_r_.start_loc_[2];
                                    const S32 k_loc = k - le_r_.start_loc_[3];
                                    const S32 adr_src  = i 
                                                       + le_r_.size_glb_[1] * (j
                                                       + le_r_.size_glb_[2] * k);
                                    const S32 adr_dest = lm_loc
                                                       + le_r_.size_loc_[0] * (i_loc
                                                       + le_r_.size_loc_[1] * (j_loc
                                                       + le_r_.size_loc_[2] * k_loc));
                                    le_r_.buf_[adr_dest] = norm * rbuf_[adr_src];
                                }
                            }
                            time_profile_.M2L_le_k_to_le_r += GetWtime() - time_start; 
                        }
                    } // END of if (rank_in_my_group_ == 0)
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                }
#endif

#if 0
                // Check le_r
                std::stringstream ss;
                ss << "le_r_" << std::setfill('0') << std::setw(5)
                   << rank_in_parent_group_ << ".txt";
                const std::string filename = ss.str();
                std::ofstream output_file;
                output_file.open(filename.c_str(), std::ios::trunc);
                for (S32 k = le_r_.start_loc_used_[3]; k <= le_r_.end_loc_used_[3]; k++)
                    for (S32 j = le_r_.start_loc_used_[2]; j <= le_r_.end_loc_used_[2]; j++)
                        for (S32 i = le_r_.start_loc_used_[1]; i <= le_r_.end_loc_used_[1]; i++)
                            for (S32 lm = le_r_.start_loc_used_[0]; lm <= le_r_.end_loc_used_[0]; lm++)
                {
                    const S32 idx = lm 
                                  + le_r_.size_glb_[0] * (i
                                  + le_r_.size_glb_[1] * (j
                                  + le_r_.size_glb_[2] * k));
                    const S32 lm_loc = lm - le_r_.start_loc_[0];
                    const S32 i_loc  = i  - le_r_.start_loc_[1];
                    const S32 j_loc  = j  - le_r_.start_loc_[2];
                    const S32 k_loc  = k  - le_r_.start_loc_[3];
                    const S32 adr = lm_loc
                                  + le_r_.size_loc_[0] * (i_loc
                                  + le_r_.size_loc_[1] * (j_loc
                                  + le_r_.size_loc_[2] * k_loc));
                    if (le_r_.buf_[adr] != 0) {
                        output_file << idx << "    " << le_r_.buf_[adr] << std::endl;
                    }
                }
                output_file.close();
                Finalize();
                std::exit(0);
#endif
            }


            template <class Cell_t>
            void redistLE(const MultientranceArray<S32> & rnklst,
                          const S32 n_cell_loc,
                          std::unordered_map<S32, S32> & adr_cell_loc, 
                          std::vector<Cell_t> & cell_loc) {
                // Extract basic parameters
                const S32 n_mm_compo = (param_.p + 1) * (param_.p + 1);
                const S32vec n_cell = param_.n_cell;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                if (mode_ == 1) {
                    using pair_t = std::pair<S32, real_t>;
                    F64 time_start = GetWtime();
                    // Make a send buffer
                    CommBuffer<pair_t> sendbuf;
                    sendbuf.n_comm = n_proc_in_parent_group_;
                    sendbuf.allocCommInfo();
                    for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                        const S32 rnk = (rank_in_parent_group_ + i) % n_proc_in_parent_group_;
                        sendbuf.ranks[i] = rnk;
                        sendbuf.adr_from_rank[rnk] = i;
                    }
                    sendbuf.clearCounts();
                    sendbuf.count_tot = 0;
                    if (rank_in_my_group_ < n_proc_for_fft_) { // this process has LE data
                        const S32 iz_start = local_0_start_[rank_in_my_group_];
                        const S32 iz_end   = local_0_end_[rank_in_my_group_];
                        for (S32 iz = iz_start; iz <= iz_end; iz++) {
                            if (iz < n_cell.z) {
                                for (S32 iy = 0; iy < n_cell.y; iy++) {
                                    for (S32 ix = 0; ix < n_cell.x; ix++) {
                                        const S32 idx = ix + n_cell.x * (iy + n_cell.y * iz);
                                        const S32 count = rnklst.counts[idx];
                                        const S32 displ = rnklst.displs[idx];
                                        for (S32 i = 0; i < count; i++) {
                                            const S32 rnk = rnklst.data[i + displ];
                                            const S32 adr = sendbuf.adr_from_rank[rnk];
                                            sendbuf.counts[adr]++;
                                            sendbuf.count_tot++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (sendbuf.count_tot > 0) {
                        sendbuf.allocBuffer();
                        sendbuf.calcDispls();
                        sendbuf.clearCounts();
                        const S32 iz_start = local_0_start_[rank_in_my_group_];
                        const S32 iz_end   = local_0_end_[rank_in_my_group_];
                        for (S32 iz = iz_start; iz <= iz_end; iz++) {
                            if (iz < n_cell.z) {
                                for (S32 iy = 0; iy < n_cell.y; iy++) {
                                    for (S32 ix = 0; ix < n_cell.x; ix++) {
                                        const S32 idx = ix + n_cell.x * (iy + n_cell.y * iz);
                                        const S32 displ = rnklst.displs[idx];
                                        const S32 count = rnklst.counts[idx];
                                        for (S32 i = 0; i < count; i++) {
                                            const S32 rnk = rnklst.data[i + displ];
                                            const S32 adr = sendbuf.adr_from_rank[rnk];
                                            const S32 i_loc = ix - le_r_.start_loc_[1];
                                            const S32 j_loc = iy - le_r_.start_loc_[2];
                                            const S32 k_loc = iz - le_r_.start_loc_[3];
                                            const S32 adr_src = i_loc
                                                              + le_r_.size_loc_[1] * (j_loc
                                                              + le_r_.size_loc_[2] * k_loc);
                                            const S32 adr_dest = sendbuf.displs[adr] 
                                                               + sendbuf.counts[adr];
                                            sendbuf.buf[adr_dest].first = idx;
                                            sendbuf.buf[adr_dest].second = le_r_.buf_[adr_src];
                                            sendbuf.counts[adr]++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // Make a receive buffer
                    CommBuffer<pair_t> recvbuf;
                    recvbuf.n_comm = n_proc_in_parent_group_;
                    recvbuf.allocCommInfo();
                    for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                        recvbuf.ranks[i] = i;
                        recvbuf.adr_from_rank[i] = i;
                    }
                    recvbuf.clearCounts();
                    recvbuf.count_tot = 0;
                    for (S32 n = 0; n < n_cell_loc; n++) {
                        const S32 idx = cell_loc[n].idx;
                        const S32 iz = idx/(n_cell.x * n_cell.y);
                        int slab_id;
                        for (S32 k = 0; k < n_proc_for_fft_; k++) {
                            if ((local_0_start_[k] <= iz) && (iz <= local_0_end_[k])) {
                                slab_id = k;
                                break;
                            }
                        }
                        for (S32 k = 0; k < n_group_; k++) {
                            const S32 rnk = rank_start_[k] + slab_id;
                            recvbuf.counts[rnk]++;
                            recvbuf.count_tot++; 
                        }
                    }
                    if (recvbuf.count_tot > 0) {
                        recvbuf.allocBuffer();
                        recvbuf.calcDispls();
                    }
                    time_profile_.M2L_preproc_le_r_comm += GetWtime() - time_start;
                    // Perform MPI comm.
                    time_start = GetWtime();
                    performComm(sendbuf, recvbuf);
                    time_profile_.M2L_redist_le_r += GetWtime() - time_start;
                    // Copy from recvbuf to cell_loc
                    time_start = GetWtime();
                    for (S32 n = 0; n < n_cell_loc; n++) cell_loc[n].le.buf.resize(0);
                    for (S32 i = 0; i < recvbuf.count_tot; i++) {
                        const S32 idx = recvbuf.buf[i].first; // idx
                        const S32 adr = adr_cell_loc[idx];
                        // [Note (tag: #cee527db)]
                        //     In C++11, we cannot use operator [] to access 
                        //     const std::unordered_map<,> and have to use member
                        //     function `at` as follows:
                        //
                        //     adr = adr_cell_loc.at(idx);
                        // 
                        //     However, Fujitsu C++ compiler does not support `at`.
                        //     Hence, we had to remove a const qualifier from
                        //     adr_cell_loc and used operator[].
                        const real_t val = recvbuf.buf[i].second; // le
                        cell_loc[adr].le.buf.push_back(val);
                    }
                    time_profile_.M2L_postproc_le_r_comm += GetWtime() - time_start;
#ifdef DEBUG_M2L_ENGINE_REDISTLE
                    for (S32 n = 0; n < n_cell_loc; n++) {
                        const S32 size = cell_loc[n].le.buf.size();
                        if (size > 0 && size != n_mm_compo) {
                            std::cout << "cell_loc[].le is something wrong:" << std::endl;
                            std::cout << " my_rank = " << rank_in_parent_group_ 
                                      << " adr = " << n
                                      << " idx = " << cell_loc[n].idx
                                      << std::endl;
                            assert(false);
                        }
                    }
#endif
                } else if (mode_ == 2) {
                    F64 time_start = GetWtime();
                    // Make a send buffer
                    CommBuffer<S32> idx_send;
                    CommBuffer<real_t> le_send;
                    idx_send.n_comm = n_proc_in_parent_group_;
                    le_send.n_comm  = n_proc_in_parent_group_;
                    idx_send.allocCommInfo();
                    le_send.allocCommInfo();
                    for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                        idx_send.ranks[i] = i;
                        le_send.ranks[i] = i;
                    }
                    idx_send.clearCounts();
                    le_send.clearCounts();
                    idx_send.count_tot = 0;
                    le_send.count_tot = 0;
                    const S32 n_lm = le_r_.size_loc_used_[0];
                    for (S32 iz = le_r_.start_loc_used_[3]; iz <= le_r_.end_loc_used_[3]; iz++) {
                        if (iz < n_cell.z) {
                            for (S32 iy = le_r_.start_loc_used_[2]; iy <= le_r_.end_loc_used_[2]; iy++) {
                                if (iy < n_cell.y) {
                                    for (S32 ix = le_r_.start_loc_used_[1]; ix <= le_r_.end_loc_used_[1]; ix++) {
                                        if (ix < n_cell.x) {
                                            const S32 idx = ix + n_cell.x * (iy + n_cell.y * iz);
                                            const S32 count = rnklst.counts[idx];
                                            const S32 displ = rnklst.displs[idx];
                                            for (S32 i = 0; i < count; i++) {
                                                const S32 rnk = rnklst.data[i + displ];
                                                if (rnk != -1) {
                                                    idx_send.counts[rnk]++;
                                                    idx_send.count_tot++;
                                                    le_send.counts[rnk] += n_lm;
                                                    le_send.count_tot += n_lm;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (idx_send.count_tot > 0) {
                        idx_send.allocBuffer();
                        idx_send.calcDispls();
                        idx_send.clearCounts();
                        le_send.allocBuffer();
                        le_send.calcDispls();
                        le_send.clearCounts();
                        const S32 n_lm = le_r_.size_loc_used_[0];
                        for (S32 iz = le_r_.start_loc_used_[3]; iz <= le_r_.end_loc_used_[3]; iz++) {
                            if (iz < n_cell.z) {
                                for (S32 iy = le_r_.start_loc_used_[2]; iy <= le_r_.end_loc_used_[2]; iy++) {
                                    if (iy < n_cell.y) {
                                        for (S32 ix = le_r_.start_loc_used_[1]; ix <= le_r_.end_loc_used_[1]; ix++) {
                                            if (ix < n_cell.x) {
                                                const S32 idx = ix + n_cell.x * (iy + n_cell.y * iz);
                                                const S32 count = rnklst.counts[idx];
                                                const S32 displ = rnklst.displs[idx];
                                                for (S32 i = 0; i < count; i++) {
                                                    const S32 rnk = rnklst.data[i + displ];
                                                    if (rnk != -1) {
                                                        // set idx_send
                                                        const S32 adr = idx_send.displs[rnk]
                                                                      + idx_send.counts[rnk];
                                                        idx_send.buf[adr] = idx;
                                                        idx_send.counts[rnk]++;
                                                        // set le_send
                                                        for (S32 lm = le_r_.start_loc_used_[0]; lm <= le_r_.end_loc_used_[0]; lm++) {
                                                            const S32 lm_loc = lm - le_r_.start_loc_[0];
                                                            const S32 i_loc  = ix - le_r_.start_loc_[1];
                                                            const S32 j_loc  = iy - le_r_.start_loc_[2];
                                                            const S32 k_loc  = iz - le_r_.start_loc_[3];
                                                            const S32 adr_src = lm_loc
                                                                              + le_r_.size_loc_[0] * (i_loc
                                                                              + le_r_.size_loc_[1] * (j_loc
                                                                              + le_r_.size_loc_[2] * k_loc));
                                                            const S32 adr_dest = le_send.displs[rnk]
                                                                               + le_send.counts[rnk];
                                                            le_send.buf[adr_dest] = le_r_.buf_[adr_src];
                                                            le_send.counts[rnk]++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // Make a receive buffer
                    CommBuffer<S32> idx_recv;
                    CommBuffer<real_t> le_recv;
                    idx_recv.n_comm = n_proc_in_parent_group_;
                    le_recv.n_comm  = n_proc_in_parent_group_;
                    idx_recv.allocCommInfo();
                    le_recv.allocCommInfo();
                    for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                        idx_recv.ranks[i] = i;
                        le_recv.ranks[i] = i;
                    }
                    idx_recv.clearCounts();
                    le_recv.clearCounts();
                    idx_recv.count_tot = 0;
                    le_recv.count_tot = 0;
                    for (S32 n = 0; n < n_cell_loc; n++) {
                        const S32 idx = cell_loc[n].idx;
                        for (S32 i = 0; i < n_proc_in_parent_group_; i++) {
                            const S32 n_lm = lm_end_[i] - lm_start_[i] + 1;
                            idx_recv.counts[i]++;
                            idx_recv.count_tot++;
                            le_recv.counts[i] += n_lm;
                            le_recv.count_tot += n_lm;
                        }
                    }
                    if (idx_recv.count_tot > 0) {
                        idx_recv.allocBuffer();
                        le_recv.allocBuffer();
                        idx_recv.calcDispls();
                        le_recv.calcDispls();
                    }
                    time_profile_.M2L_preproc_le_r_comm += GetWtime() - time_start;
                    // Perform MPI comm.
                    time_start = GetWtime();
                    performComm(idx_send, idx_recv);
                    performComm(le_send, le_recv);
                    time_profile_.M2L_redist_le_r += GetWtime() - time_start;
                    // Copy receive buffers to cell_loc[]
                    time_start = GetWtime();
                    for (S32 n = 0; n < idx_recv.n_comm; n++) {
                        const S32 rnk = idx_recv.ranks[n];
                        const S32 disp_for_idx = idx_recv.displs[n];
                        const S32 cnt = idx_recv.counts[n];
                        const S32 disp_for_le = le_recv.displs[n];
                        S32 k = 0; 
                        for (S32 i = 0; i < cnt; i++) {
                            const S32 idx = idx_recv.buf[disp_for_idx + i];
                            const S32 adr = adr_cell_loc[idx]; // see Note #cee527db.
                            const S32 beg = lm_start_[rnk];
                            const S32 end = lm_end_[rnk];
                            for (S32 lm = beg; lm <= end; lm++) {
                                cell_loc[adr].le.buf[lm] = le_recv.buf[disp_for_le + k];
                                k++;
                            }
                        }
                    }
                    time_profile_.M2L_postproc_le_r_comm += GetWtime() - time_start;
                }
#else // PARTICLE_SIMULATOR_MPI_PARALLEL
                // Copy from le_r_[] to cell_loc[] 
                F64 time_start = GetWtime();
                for (S32 k = le_r_.start_loc_used_[3]; k <= le_r_.end_loc_used_[3]; k++) {
                    if (k < n_cell.z) {
                        for (S32 j = le_r_.start_loc_used_[2]; j <= le_r_.end_loc_used_[2]; j++) {
                            if (j < n_cell.y) {
                                for (S32 i = le_r_.start_loc_used_[1]; i <= le_r_.end_loc_used_[1]; i++) {
                                    if (i < n_cell.x) {
                                        for (S32 lm = le_r_.start_loc_used_[0]; lm <= le_r_.end_loc_used_[0]; lm++) {
                                            const S32 lm_loc = lm - le_r_.start_loc_[0];
                                            const S32 i_loc  = i  - le_r_.start_loc_[1];
                                            const S32 j_loc  = j  - le_r_.start_loc_[2];
                                            const S32 k_loc  = k  - le_r_.start_loc_[3];
                                            const S32 adr_src = lm_loc
                                                              + le_r_.size_loc_[0] * (i_loc
                                                              + le_r_.size_loc_[1] * (j_loc
                                                              + le_r_.size_loc_[2] * k_loc));
                                            const S32 idx = i + n_cell.x * (j + n_cell.y * k);
                                            const S32 adr_dest = adr_cell_loc[idx]; // see Note #cee527db.
                                            cell_loc[adr_dest].le.buf[lm] = le_r_.buf_[adr_src];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                time_profile_.M2L_redist_le_r += GetWtime() - time_start;
#endif // PARTICLE_SIMULATOR_MPI_PARALLEL

#if 0
                // Check le_r
                std::stringstream ss;
                ss << "le_r_" << std::setfill('0') << std::setw(5)
                   << rank_in_parent_group_ << ".txt";
                const std::string filename = ss.str();
                std::ofstream output_file;
                output_file.open(filename.c_str(), std::ios::trunc);
                for (S32 i = 0; i < n_cell_loc; i++) {
                    const S32 idx_1d = cell_loc[i].idx;
                    S32vec idx_3d;
                    idx_3d.z = idx_1d / (n_cell.x * n_cell.y);
                    idx_3d.y = (idx_1d - (n_cell.x * n_cell.y) * idx_3d.z) / n_cell.x;
                    idx_3d.x = idx_1d - (n_cell.x * n_cell.y) * idx_3d.z - n_cell.x * idx_3d.y; 
                    for (S32 lm = 0; lm < cell_loc[i].le.buf.size(); lm++) {
                        const S32 idx = lm 
                                      + n_mm_compo * (idx_3d.x
                                      + n_cell.x * (idx_3d.y
                                      + n_cell.y * idx_3d.z));
                        const real_t val = cell_loc[i].le.buf[lm];
                        if (val != 0) {
                            output_file << idx << "    " << val << std::endl;
                        }
                    }
                }
                output_file.close();
                Finalize();
                std::exit(0);
#endif
            }

            void clearTimeProfile() {
                time_profile_.clear();
            }

            TimeProfilePMM getTimeProfile() const {
                return time_profile_;
            }

        private:

            void calcGreenFunctionInRealSpace(const F64vec cell_length) {
                if (param_.bc == BOUNDARY_CONDITION_OPEN) {
                    calcGreenFunctionInRealSpaceForOpenBC(cell_length);
                } else if (param_.bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                    calcGreenFunctionInRealSpaceForPeriodicBC(cell_length);
                }
            }

            void calcGreenFunctionInRealSpaceForOpenBC(const F64vec cell_length) {
                if (rank_in_my_group_ < n_proc_for_fft_) {
                    // Extract parameters
                    const S32 p = param_.p;
                    const S32 icut = param_.icut;
                    const S32vec n_cell = param_.n_cell;
                    // Local variables
                    const S32 n_thread = Comm::getNumberOfThread();
                    Slm<real_t> * slm = new Slm<real_t>[n_thread];
                    for (S32 i = 0; i < n_thread; i++) slm[i].alloc(2*p);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for (S32 k = gf_r_.start_loc_used_[3]; k <= gf_r_.end_loc_used_[3]; k++)
                        for (S32 j = gf_r_.start_loc_used_[2]; j <= gf_r_.end_loc_used_[2]; j++)
                            for (S32 i = gf_r_.start_loc_used_[1]; i <= gf_r_.end_loc_used_[1]; i++)
                    {
    
                        const S32 ith = Comm::getThreadNum();
                        const S32 i_loc = i - gf_r_.start_loc_[1];
                        const S32 j_loc = j - gf_r_.start_loc_[2];
                        const S32 k_loc = k - gf_r_.start_loc_[3];
    
                        if ((k == gf_r_.size_glb_[3]/2) || 
                            (j == gf_r_.size_glb_[2]/2) || 
                            (i == gf_r_.size_glb_[1]/2)) {
                            for (S32 lm = gf_r_.start_loc_used_[0]; lm <= gf_r_.end_loc_used_[0]; lm++) {
                                const S32 lm_loc = lm - gf_r_.start_loc_[0];
                                const S32 adr = lm_loc
                                              + gf_r_.size_loc_[0] * (i_loc
                                              + gf_r_.size_loc_[1] * (j_loc
                                              + gf_r_.size_loc_[2] * k_loc));
                                gf_r_.buf_[adr] = 0.0;
                            }
                            continue;
                        }
                        const S32 kk = (k > gf_r_.size_glb_[3]/2) ? k - gf_r_.size_glb_[3] : k;
                        const S32 jj = (j > gf_r_.size_glb_[2]/2) ? j - gf_r_.size_glb_[2] : j;
                        const S32 ii = (i > gf_r_.size_glb_[1]/2) ? i - gf_r_.size_glb_[1] : i;
                        // [TODO] 
                        //    In the calculation above, we must also use gf_r_.start_glb_[]
                        //    and gf_r_.end_glb_[] if gf_r_.start_glb_[?] != 0.
                        if (abs(kk) <= icut && 
                            abs(jj) <= icut && 
                            abs(ii) <= icut) {
                            for (S32 lm = gf_r_.start_loc_used_[0]; lm <= gf_r_.end_loc_used_[0]; lm++) {
                                const S32 lm_loc = lm - gf_r_.start_loc_[0];
                                const S32 adr = lm_loc
                                              + gf_r_.size_loc_[0] * (i_loc
                                              + gf_r_.size_loc_[1] * (j_loc
                                              + gf_r_.size_loc_[2] * k_loc));
                                gf_r_.buf_[adr] = 0.0;
                            }
                            continue;
                        }
                        const F64 dx = cell_length.x * F64(ii);
                        const F64 dy = cell_length.y * F64(jj);
                        const F64 dz = cell_length.z * F64(kk);
                        slm[ith].eval_opt(-dx, dy, dz); // eval S_l^{-m}
                        for (S32 lm = gf_r_.start_loc_used_[0]; lm <= gf_r_.end_loc_used_[0]; lm++) {
                            const S32 lm_loc = lm - gf_r_.start_loc_[0];
                            const S32 adr = lm_loc
                                          + gf_r_.size_loc_[0] * (i_loc
                                          + gf_r_.size_loc_[1] * (j_loc
                                          + gf_r_.size_loc_[2] * k_loc));
                            gf_r_.buf_[adr] = slm[ith].buf[lm];
                        }
                    } 
                    // Free memory
                    for (S32 i = 0; i < n_thread; i++) slm[i].freeMem();
                    delete [] slm;
                }
            }
    
            void calcGreenFunctionInRealSpaceForPeriodicBC(const F64vec cell_length) {
                if (rank_in_my_group_ < n_proc_for_fft_) {
                    // Extract parameters
                    const S32 p = param_.p;
                    const S32 icut = param_.icut;
                    const F64vec period_length = param_.pos_unit_cell.getFullLength();
                    const S32vec n_cell = param_.n_cell;
                    const F64 alpha = param_.alpha;
                    const S32 NMAX = param_.nmax;
                    const S32 MMAX = param_.mmax;
                    // Local variables
                    const S32 n_thread = Comm::getNumberOfThread();
                    const F64 bx = period_length.x;
                    const F64 by = period_length.y;
                    const F64 bz = period_length.z;
                    const F64 bx_inv = 1.0/bx;
                    const F64 by_inv = 1.0/by;
                    const F64 bz_inv = 1.0/bz;
                    const F64 bx2_inv = 1.0/(bx*bx);
                    const F64 by2_inv = 1.0/(by*by);
                    const F64 bz2_inv = 1.0/(bz*bz);
                    const F64 vol_inv = 1.0/(bx*by*bz);
                    CutFunc<real_t, cplx_t> *cf;
                    cf = new CutFunc<real_t, cplx_t>[n_thread];
                    Slm<real_t> *slm, *rsum, *ksum;
                    slm  = new Slm<real_t>[n_thread];
                    rsum = new Slm<real_t>[n_thread];
                    ksum = new Slm<real_t>[n_thread];
                    for (S32 i = 0; i < n_thread; i++) {
                        cf[i].init(2*p);
                        slm[i].alloc(2*p);
                        rsum[i].alloc(2*p);
                        ksum[i].alloc(2*p);
                    }
                    {
                        Slm<real_t> tmp;
                        tmp.dry_run(2*p);
                    }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for (S32 k = gf_r_.start_loc_used_[3]; k <= gf_r_.end_loc_used_[3]; k++)
                         for (S32 j = gf_r_.start_loc_used_[2]; j <= gf_r_.end_loc_used_[2]; j++)
                             for (S32 i = gf_r_.start_loc_used_[1]; i <= gf_r_.end_loc_used_[1]; i++)
                    {
    
                        const S32 ith = Comm::getThreadNum();
                        const S32 i_loc = i - gf_r_.start_loc_[1];
                        const S32 j_loc = j - gf_r_.start_loc_[2];
                        const S32 k_loc = k - gf_r_.start_loc_[3];
    
                        // real-space sum
                        rsum[ith].clear();
                        const S32 kk = (k >= gf_r_.size_glb_[3]/2) ? k - gf_r_.size_glb_[3] : k;
                        const S32 jj = (j >= gf_r_.size_glb_[2]/2) ? j - gf_r_.size_glb_[2] : j;
                        const S32 ii = (i >= gf_r_.size_glb_[1]/2) ? i - gf_r_.size_glb_[1] : i;
                        // [TODO] 
                        //    In the calculation above, we must also use gf_r.start_glb_[]
                        //    and gf_r.end_glb_[] if gf_r.start_glb_[?] != 0.
                        S32 nx, ny, nz;
                        for (nz = -NMAX; nz <= NMAX; nz++)
                            for (ny = -NMAX; ny <= NMAX; ny++)
                                for (nx = -NMAX; nx <= NMAX; nx++)
                        {
                            const S32 kkk = kk + nz * gf_r_.size_glb_[3];
                            const S32 jjj = jj + ny * gf_r_.size_glb_[2];
                            const S32 iii = ii + nx * gf_r_.size_glb_[1];
                            if( 0 == (iii|jjj|kkk) ) continue;
                            const F64 dx = cell_length.x * F64(iii);
                            const F64 dy = cell_length.y * F64(jjj);
                            const F64 dz = cell_length.z * F64(kkk);
                            const F64 dr2 = dx*dx + dy*dy + dz*dz;
                            cf[ith].eval_rcut(dr2 * (alpha*alpha));
                            slm[ith].eval_opt(-dx, dy, dz); // eval S_l^{-m}
                            // near cell correction
                            const bool near = (abs(kkk)<=icut && abs(jjj)<=icut && abs(iii)<=icut);
                            if (near) {
                                for(S32 l=0; l<=2*p; l++){
                                    cf[ith].rcut[l] -= 1.0;
                                }
                            }
                            for(S32 l=0; l<=2*p; l++){
                                for(S32 m=0; m<=l; m++){
                                    const cplx_t val = cf[ith].rcut[l] * slm[ith].val_at(l, m);
                                    rsum[ith].accum_at(l, m, val);
                                }
                            }
                        }
           
                        // wave-space sum
                        ksum[ith].clear();
                        S32 mx, my, mz;
                        for (mz = -MMAX; mz <= MMAX; mz++)
                            for (my = -MMAX; my <= MMAX; my++)
                                for (mx = -MMAX; mx <= MMAX; mx++)
                        {
                            if(0 == (mx|my|mz)) continue;
                            const F64 dx = cell_length.x * F64(ii);
                            const F64 dy = cell_length.y * F64(jj);
                            const F64 dz = cell_length.z * F64(kk);
                            //---
                            // The following code is the same as that of the original PMMM code.
                            // But, the use of (i,j,k) is probably WRONG.
                            // const F64 dx = cell_length.x * F64(i);
                            // const F64 dy = cell_length.y * F64(j);
                            // const F64 dz = cell_length.z * F64(k);
                            //---
                            const F64 twopi = 8.0 * std::atan(1.0);
                            const F64 theta = twopi * (dx*mx*bx_inv
                                                      +dy*my*by_inv
                                                      +dz*mz*bz_inv);
                            const cplx_t phase(cos(theta), sin(theta));
                            slm[ith].eval_opt(-mx*bx_inv, my*by_inv, mz*bz_inv);
                            const F64 kn2 = mx*mx*bx2_inv + my*my*by2_inv + mz*mz*bz2_inv;
                            cf[ith].eval_kcut(kn2, alpha);
                            for (S32 l=0; l<=2*p; l++){
                                for (S32 m=0; m<=l; m++){
                                    const cplx_t val = (cf[ith].kcut[l] * phase) * slm[ith].val_at(l, m) * vol_inv;
                                    ksum[ith].accum_at(l, m, val);
                                }
                            }
                        }
                        // store sum
#ifdef PARTICLE_SIMULATOR_PMM_IGNORE_RSPACE
                        rsum[ith].clear();
#endif
#ifdef PARTICLE_SIMULATOR_PMM_IGNORE_KSPACE
                        ksum[ith].clear();
#endif
                        for (S32 lm = gf_r_.start_loc_used_[0]; lm <= gf_r_.end_loc_used_[0]; lm++) {
                            const S32 lm_loc = lm - gf_r_.start_loc_[0];
                            const S32 adr = lm_loc
                                          + gf_r_.size_loc_[0] * (i_loc
                                          + gf_r_.size_loc_[1] * (j_loc
                                          + gf_r_.size_loc_[2] * k_loc));
                            gf_r_.buf_[adr] = rsum[ith].buf[lm] + ksum[ith].buf[lm];
                        }
#ifdef PARTICLE_SIMULATOR_PMM_SHOW_RSUM_AND_KSUM
                        if (0 == (i|j|k)) {
                            rsum.show();
                            ksum.show();
                        };
#endif
                    } // for(k, j, i)
           
                    // Free memory
                    for (S32 i = 0; i < n_thread; i++) {
                        cf[i].freeMem();
                        slm[i].freeMem();
                        rsum[i].freeMem();
                        ksum[i].freeMem();
                    }
                    delete [] cf;
                    delete [] slm;
                    delete [] rsum;
                    delete [] ksum;
                }
            }

            void calcGreenFunctionInWavenumberSpace() {
                // Extract basic parameters
                const S32 p = param_.p;
                const S32vec n_cell = param_.n_cell;
                // Local variables
                fftw_real_t *rbuf = nullptr;
                fftw_cplx_t *kbuf = nullptr;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                if (mode_ == 1 && use_mpifft_) {
                    if (rank_in_my_group_ < n_proc_for_fft_) {
                        // Check if conditions that need to be fulfilled in this IF branch
                        // are actually fulfilled.
                        assert(gf_r_.start_loc_[0] == gf_k_.start_loc_[0]);
                        assert(gf_r_.end_loc_[0] == gf_k_.end_loc_[0]);
                        assert(gf_r_.start_loc_used_[0] == gf_k_.start_loc_used_[0]);
                        assert(gf_r_.end_loc_used_[0] == gf_k_.end_loc_used_[0]);

                        // Allocate buffer memory
                        ptrdiff_t alloc_local;
                        ptrdiff_t local_n0, local_0_start;
                        ptrdiff_t local_n1, local_1_start;
                        alloc_local = fftw_mpi_local_size_3d_transposed(gf_k_.size_glb_[3],
                                                                        gf_k_.size_glb_[2],
                                                                        gf_k_.size_glb_[1],
                                                                        comm_fft_,
                                                                        &local_n0, &local_0_start,
                                                                        &local_n1, &local_1_start);
                        // Note that an user must indicate the logical size of FFT by
                        // the sizes of a COMPLEX array when using fftw_mpi_local_size_3d_transposed.
                        const ptrdiff_t size_rbuf = 2 * alloc_local;
                        const ptrdiff_t size_kbuf = alloc_local;
                        rbuf = fftw_alloc_real(size_rbuf);
                        kbuf = fftw_alloc_complex(size_kbuf);

                        // Create plan
                        fftw_plan plan_fwd = 
                            fftw_mpi_plan_dft_r2c_3d(
                                gf_r_.size_glb_[3], gf_r_.size_glb_[2], gf_r_.size_glb_[1], 
                                &rbuf[0], &kbuf[0], comm_fft_,
                                FFTW_ESTIMATE | FFTW_DESTROY_INPUT | FFTW_MPI_TRANSPOSED_OUT);
                        // Note that fftw_mpi_plan_dft_r2c_3d requires an user to indicate
                        // the logical size of FFT by the size of a REAL array.

                        // Perform FFT
                        for (S32 lm = gf_r_.start_loc_used_[0]; lm <= gf_r_.end_loc_used_[0]; lm++) {
                            const S32 lm_loc = lm - gf_r_.start_loc_[0];
                            for (S32 k = gf_r_.start_loc_used_[3]; k <= gf_r_.end_loc_used_[3]; k++)
                                for (S32 j = gf_r_.start_loc_used_[2]; j <= gf_r_.end_loc_used_[2]; j++)
                                    for (S32 i = gf_r_.start_loc_used_[1]; i <= gf_r_.end_loc_used_[1]; i++)
                            {
                                const S32 i_loc = i - gf_r_.start_loc_[1];
                                const S32 j_loc = j - gf_r_.start_loc_[2];
                                const S32 k_loc = k - gf_r_.start_loc_[3];
                                const S32 adr_src  = lm_loc
                                                   + gf_r_.size_loc_[0] * (i_loc
                                                   + gf_r_.size_loc_[1] * (j_loc
                                                   + gf_r_.size_loc_[2] * k_loc));
                                const S32 i_buf = i - gf_r_.start_loc_used_[1];
                                const S32 j_buf = j - gf_r_.start_loc_used_[2];
                                const S32 k_buf = k - gf_r_.start_loc_used_[3];
                                const S32 adr_dest = i_buf 
                                                   + (2*(gf_r_.size_glb_[1]/2 + 1)) * (j_buf
                                                   + gf_r_.size_glb_[2] * k_buf);
                                // Note that the size of 1st dimension of rbuf is not NX.
                                // (see Section 6.5 of the mamual of FFTW3)
                                rbuf[adr_dest] = gf_r_.buf_[adr_src];
                            }

                            // CALL FFTW
                            fftw_execute(plan_fwd);
   
                            for (S32 j = gf_k_.start_loc_used_[3]; j <= gf_k_.end_loc_used_[3]; j++)
                                for (S32 k = gf_k_.start_loc_used_[2]; k <= gf_k_.end_loc_used_[2]; k++)
                                    for (S32 i = gf_k_.start_loc_used_[1]; i <= gf_k_.end_loc_used_[1]; i++)
                            {
                                const S32 i_buf = i - gf_k_.start_loc_used_[1];
                                const S32 k_buf = k - gf_k_.start_loc_used_[2];
                                const S32 j_buf = j - gf_k_.start_loc_used_[3];
                                const S32 adr_src = i_buf 
                                                  + gf_k_.size_glb_[1] * (k_buf
                                                  + gf_k_.size_glb_[2] * j_buf);
                                const S32 i_loc = i - gf_k_.start_loc_[1];
                                const S32 k_loc = k - gf_k_.start_loc_[2];
                                const S32 j_loc = j - gf_k_.start_loc_[3];
                                const S32 adr_dest = lm_loc
                                                   + gf_k_.size_loc_[0] * (i_loc
                                                   + gf_k_.size_loc_[1] * (k_loc
                                                   + gf_k_.size_loc_[2] * j_loc));
                                // Note that j and k loops are exchanged because we specified
                                // FFTW_MPI_TRANSPOSED_OUT.
                                assert(0 <= adr_src  && adr_src < size_kbuf);
                                assert(0 <= adr_dest && adr_dest < gf_k_.size_loc_tot_);
#ifdef PARTICLE_SIMULATOR_STD_COMPLEX_NOT_HAVING_SETTER
                                const cplx_t ctmp (kbuf[adr_src][0],
                                                   kbuf[adr_src][1]); 
                                gf_k_.buf_[adr_dest] = ctmp;
#else
                                gf_k_.buf_[adr_dest].real(kbuf[adr_src][0]);
                                gf_k_.buf_[adr_dest].imag(kbuf[adr_src][1]);
#endif
                                // Note that fftw_complex = double[2]
                            }
                        }
                        // Destropy the plan
                        fftw_destroy_plan(plan_fwd);

                    }
                } else {
#endif
                    // In this case, we can use only a single MPI process to calculate
                    // gf_k for a particular set of (l,m).
                    // We use OpenMP parallerization of FFT if available.
                    if (rank_in_my_group_ == 0) {

                        // Check if conditions that need to be fulfilled in this IF branch
                        // are actually fulfilled.
                        assert(gf_r_.start_loc_[0] == gf_k_.start_loc_[0]);
                        assert(gf_r_.end_loc_[0] == gf_k_.end_loc_[0]);
                        assert(gf_r_.start_loc_used_[0] == gf_k_.start_loc_used_[0]);
                        assert(gf_r_.end_loc_used_[0] == gf_k_.end_loc_used_[0]);
                        for (S32 i = 1; i <= 3; i++) {
                            assert(gf_r_.size_loc_[i] == gf_r_.size_glb_[i]);
                            assert(gf_k_.size_loc_[i] == gf_k_.size_glb_[i]);
                        }

                        // Calculate gf_k 
                        const S32 n_thread = Comm::getNumberOfThread();
                        const S32 fft_size = gf_r_.size_glb_[3] * gf_r_.size_glb_[2] * gf_r_.size_glb_[1];
                        if (fft_size < fft_size_crit_ && n_thread > 1) {
                            // In this case, the size of array is too small to speed up by multithreaded FFT.
                            // Hence, we assign a single FFT to each thread.

                            // Allocate buffer memory
                            const S32 size_rbuf = n_thread * fft_size;
                            const S32 size_kbuf = n_thread * fft_size ;
                            rbuf = fftw_alloc_real(size_rbuf);
                            kbuf = fftw_alloc_complex(size_kbuf);

                            // Create plan
                            fftw_plan *plan_fwd;
                            plan_fwd = (fftw_plan *) malloc((size_t) sizeof(fftw_plan) * n_thread);
                            for (S32 i=0; i < n_thread; i++) {
                                const S32 offset = fft_size * i;
                                plan_fwd[i] = fftw_plan_dft_r2c_3d(gf_r_.size_glb_[3],
                                                                   gf_r_.size_glb_[2],
                                                                   gf_r_.size_glb_[1], 
                                                                   &rbuf[offset], 
                                                                   &kbuf[offset],
                                                                   FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
                            }

                            // Perform FFT
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for 
#endif
                            for (S32 lm = gf_r_.start_loc_used_[0]; lm <= gf_r_.end_loc_used_[0]; lm++){
                                const S32 ith = Comm::getThreadNum();
                                const S32 offset = fft_size * ith;
    
                                const S32 lm_loc = lm - gf_r_.start_loc_[0];
                                for (S32 k = gf_r_.start_loc_used_[3]; k <= gf_r_.end_loc_used_[3]; k++)
                                    for (S32 j = gf_r_.start_loc_used_[2]; j <= gf_r_.end_loc_used_[2]; j++)
                                        for (S32 i = gf_r_.start_loc_used_[1]; i <= gf_r_.end_loc_used_[1]; i++)
                                {
                                    const S32 i_loc = i - gf_r_.start_loc_[1];
                                    const S32 j_loc = j - gf_r_.start_loc_[2];
                                    const S32 k_loc = k - gf_r_.start_loc_[3];
                                    const S32 adr_src  = lm_loc
                                                       + gf_r_.size_loc_[0] * (i_loc
                                                       + gf_r_.size_loc_[1] * (j_loc
                                                       + gf_r_.size_loc_[2] * k_loc));
                                    const S32 adr_dest = i 
                                                       + gf_r_.size_glb_[1] * (j
                                                       + gf_r_.size_glb_[2] * k);
                                    rbuf[adr_dest + offset] = gf_r_.buf_[adr_src];
                                }
    
                                // CALL FFTW
                                fftw_execute(plan_fwd[ith]);
        
                                for (S32 k = gf_k_.start_loc_used_[3]; k <= gf_k_.end_loc_used_[3]; k++)
                                    for (S32 j = gf_k_.start_loc_used_[2]; j <= gf_k_.end_loc_used_[2]; j++)
                                        for (S32 i = gf_k_.start_loc_used_[1]; i <= gf_k_.end_loc_used_[1]; i++)
                                {
                                    const S32 i_loc = i - gf_k_.start_loc_[1];
                                    const S32 j_loc = j - gf_k_.start_loc_[2];
                                    const S32 k_loc = k - gf_k_.start_loc_[3];
                                    const S32 adr_dest = lm_loc
                                                       + gf_k_.size_loc_[0] * (i_loc
                                                       + gf_k_.size_loc_[1] * (j_loc
                                                       + gf_k_.size_loc_[2] * k_loc));
                                    const S32 adr_src  = i 
                                                       + gf_k_.size_glb_[1] * (j 
                                                       + gf_k_.size_glb_[2] * k);
#ifdef PARTICLE_SIMULATOR_STD_COMPLEX_NOT_HAVING_SETTER
                                    const cplx_t ctmp (kbuf[adr_src + offset][0],
                                                       kbuf[adr_src + offset][1]);
                                    gf_k_.buf_[adr_dest] = ctmp;
#else
                                    gf_k_.buf_[adr_dest].real(kbuf[adr_src + offset][0]);
                                    gf_k_.buf_[adr_dest].imag(kbuf[adr_src + offset][1]);
#endif
                                    // Note that fftw_complex = double[2]
                                }
                            }
    
                            // Destroy the plan
                            for (S32 i = 0; i < n_thread; i++) fftw_destroy_plan(plan_fwd[i]);
                            delete [] plan_fwd;
    
                        } else {
                            // In this case, we can expect that the FFT calculation become fast by using
                            // a multithreaded FFT.
    
                            // Allocate buffer memory
                            const S32 size_rbuf = gf_r_.size_glb_[3] * gf_r_.size_glb_[2] * gf_r_.size_glb_[1];
                            const S32 size_kbuf = gf_k_.size_glb_[3] * gf_k_.size_glb_[2] * gf_k_.size_glb_[1];
                            rbuf = fftw_alloc_real(size_rbuf);
                            kbuf = fftw_alloc_complex(size_kbuf);
    
                            // Create plan
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                            fftw_plan_with_nthreads(n_thread);
#endif
                            fftw_plan plan_fwd = 
                                fftw_plan_dft_r2c_3d(
                                    gf_r_.size_glb_[3], gf_r_.size_glb_[2], gf_r_.size_glb_[1], 
                                    &rbuf[0], &kbuf[0],
                                    FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    
                            // Perform FFT
                            for (S32 lm = gf_r_.start_loc_used_[0]; lm <= gf_r_.end_loc_used_[0]; lm++) {
                                const S32 lm_loc = lm - gf_r_.start_loc_[0];
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for 
#endif
                                for (S32 k = gf_r_.start_loc_used_[3]; k <= gf_r_.end_loc_used_[3]; k++)
                                    for (S32 j = gf_r_.start_loc_used_[2]; j <= gf_r_.end_loc_used_[2]; j++)
                                        for (S32 i = gf_r_.start_loc_used_[1]; i <= gf_r_.end_loc_used_[1]; i++)
                                {
                                    const S32 i_loc = i - gf_r_.start_loc_[1];
                                    const S32 j_loc = j - gf_r_.start_loc_[2];
                                    const S32 k_loc = k - gf_r_.start_loc_[3];
                                    const S32 adr_src  = lm_loc
                                                       + gf_r_.size_loc_[0] * (i_loc
                                                       + gf_r_.size_loc_[1] * (j_loc
                                                       + gf_r_.size_loc_[2] * k_loc));
                                    const S32 adr_dest = i 
                                                       + gf_r_.size_glb_[1] * (j
                                                       + gf_r_.size_glb_[2] * k);
                                    rbuf[adr_dest] = gf_r_.buf_[adr_src];
                                }
    
                                // CALL FFTW
                                fftw_execute(plan_fwd);
        
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for 
#endif
                                for (S32 k = gf_k_.start_loc_used_[3]; k <= gf_k_.end_loc_used_[3]; k++)
                                    for (S32 j = gf_k_.start_loc_used_[2]; j <= gf_k_.end_loc_used_[2]; j++)
                                        for (S32 i = gf_k_.start_loc_used_[1]; i <= gf_k_.end_loc_used_[1]; i++)
                                {
                                    const S32 i_loc = i - gf_k_.start_loc_[1];
                                    const S32 j_loc = j - gf_k_.start_loc_[2];
                                    const S32 k_loc = k - gf_k_.start_loc_[3];
                                    const S32 adr_dest = lm_loc
                                                       + gf_k_.size_loc_[0] * (i_loc
                                                       + gf_k_.size_loc_[1] * (j_loc
                                                       + gf_k_.size_loc_[2] * k_loc));
                                    const S32 adr_src  = i 
                                                       + gf_k_.size_glb_[1] * (j 
                                                       + gf_k_.size_glb_[2] * k);
#ifdef PARTICLE_SIMULATOR_STD_COMPLEX_NOT_HAVING_SETTER
                                    const cplx_t ctmp (kbuf[adr_src][0],
                                                       kbuf[adr_src][1]);
                                    gf_k_.buf_[adr_dest] = ctmp;
#else
                                    gf_k_.buf_[adr_dest].real(kbuf[adr_src][0]);
                                    gf_k_.buf_[adr_dest].imag(kbuf[adr_src][1]);
#endif
                                    // Note that fftw_complex = double[2]
                                }
                            }
    
                            // Destroy the plan
                            fftw_destroy_plan(plan_fwd);
                        } // END of if (rank_in_my_group_ == 0)
                    } 
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                }
#endif
                // Free
                if (rbuf != nullptr) fftw_free(rbuf);
                if (kbuf != nullptr) fftw_free(kbuf);
            }

        
            void transform() {
                // Extract parameters
                const S32 p = param_.p;
                const S32 n_thread = Comm::getNumberOfThread();

                // Local buffers
                Slm<cplx_t> * slm_tmp;
                MultipoleMoment<cplx_t> * mm_tmp;
                LocalExpansion<cplx_t>  * le_tmp;
                slm_tmp = new Slm<cplx_t>[n_thread];
                mm_tmp  = new MultipoleMoment<cplx_t>[n_thread];
                le_tmp  = new LocalExpansion<cplx_t>[n_thread];
                for (S32 i = 0; i < n_thread; i++) {
                    slm_tmp[i].alloc(2*p);
                    mm_tmp[i].alloc(p);
                    le_tmp[i].alloc(p);
                }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                if (mode_ == 1 && use_mpifft_) {
                    // In this case, the data of ?_k[] is stored in the format
                    // of FFTW_MPI_TRANSPOSED_OUT.
                    if (rank_in_my_group_ < n_proc_for_fft_) {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                        for (S32 j = le_k_.start_loc_used_[3]; j <= le_k_.end_loc_used_[3]; j++) {
                            const S32 ith = Comm::getThreadNum();
                            for (S32 k = le_k_.start_loc_used_[2]; k <= le_k_.end_loc_used_[2]; k++) {
                                for (S32 i = le_k_.start_loc_used_[1]; i <= le_k_.end_loc_used_[1]; i++) {
                                    const S32 i_loc = i - le_k_.start_loc_[1];
                                    const S32 k_loc = k - le_k_.start_loc_[2];
                                    const S32 j_loc = j - le_k_.start_loc_[3];
                                    // Copy from mm_k to mm_tmp
                                    for (S32 lm = mm_k_.start_glb_[0]; lm <= mm_k_.end_glb_[0]; lm++) {
                                        const S32 adr_src = lm
                                                          + mm_k_.size_loc_[0] * (i_loc
                                                          + mm_k_.size_loc_[1] * (k_loc
                                                          + mm_k_.size_loc_[2] * j_loc));
                                        mm_tmp[ith].buf[lm] = mm_k_.buf_[adr_src];
                                    }
                                    // Copy from gf_k to slm_tmp
                                    for (S32 lm = gf_k_.start_glb_[0]; lm <= gf_k_.end_glb_[0]; lm++) {
                                        const S32 adr_src = lm
                                                          + gf_k_.size_loc_[0] * (i_loc
                                                          + gf_k_.size_loc_[1] * (k_loc
                                                          + gf_k_.size_loc_[2] * j_loc));
                                        slm_tmp[ith].buf[lm] = gf_k_.buf_[adr_src];
                                    }
                                    // do M2L
                                    //slm_tmp[ith].transform_M2L(mm_tmp[ith], le_tmp[ith], false);
                                    // Copy le_tmp to le_k[]
                                    for (S32 lm = le_k_.start_loc_used_[0]; lm <= le_k_.end_loc_used_[0]; lm++) {
                                        // Find (ldst, mdst)
                                        S32 ldst = 0, mdst;
                                        do {
                                            if ((ldst + 1)*(ldst + 1) > lm) break;
                                            else ldst++;
                                        } while (1);
                                        mdst = lm - ldst * (ldst + 1);
                                        if (mdst < 0) mdst *= -1;
                                        // Compute le_k at (ldst, mdst)
                                        slm_tmp[ith].transform_M2L(mm_tmp[ith], ldst, mdst, le_tmp[ith], false);
                                        // Copy
                                        const S32 lm_loc = lm - le_k_.start_loc_[0];
                                        const S32 adr_dest = lm_loc
                                                           + le_k_.size_loc_[0] * (i_loc
                                                           + le_k_.size_loc_[1] * (k_loc
                                                           + le_k_.size_loc_[2] * j_loc));
                                        le_k_.buf_[adr_dest] = le_tmp[ith].buf[lm];
                                    }
                                }
                            }
                        }
                    }
                } else {
#endif // PARTICLE_SIMULATOR_MPI_PARALLEL
                    // In this case, the data of ?_k[] is stored in the same
                    // way as ?_r[].
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for (S32 k = le_k_.start_loc_used_[3]; k <= le_k_.end_loc_used_[3]; k++) {
                        const S32 ith = Comm::getThreadNum();
                        for (S32 j = le_k_.start_loc_used_[2]; j <= le_k_.end_loc_used_[2]; j++) {
                            for (S32 i = le_k_.start_loc_used_[1]; i <= le_k_.end_loc_used_[1]; i++) {
                                const S32 i_loc = i - le_k_.start_loc_[1];
                                const S32 j_loc = j - le_k_.start_loc_[2];
                                const S32 k_loc = k - le_k_.start_loc_[3];
                                // Copy from mm_k to mm_tmp
                                for (S32 lm = mm_k_.start_glb_[0]; lm <= mm_k_.end_glb_[0]; lm++) {
                                    const S32 adr_src = lm
                                                      + mm_k_.size_loc_[0] * (i_loc
                                                      + mm_k_.size_loc_[1] * (j_loc
                                                      + mm_k_.size_loc_[2] * k_loc));
                                    mm_tmp[ith].buf[lm] = mm_k_.buf_[adr_src];
                                }
                                // Copy from gf_k to slm_tmp
                                for (S32 lm = gf_k_.start_glb_[0]; lm <= gf_k_.end_glb_[0]; lm++) {
                                    const S32 adr_src = lm
                                                      + gf_k_.size_loc_[0] * (i_loc
                                                      + gf_k_.size_loc_[1] * (j_loc
                                                      + gf_k_.size_loc_[2] * k_loc));
                                    slm_tmp[ith].buf[lm] = gf_k_.buf_[adr_src];
                                }
                                // do M2L
                                //slm_tmp[ith].transform_M2L(mm_tmp[ith], le_tmp[ith], false);
                                // Copy le_tmp to le_k[]
                                for (S32 lm = le_k_.start_loc_used_[0]; lm <= le_k_.end_loc_used_[0]; lm++) {
                                    // Find (ldst, mdst)
                                    S32 ldst = 0, mdst;
                                    do {
                                        if ((ldst + 1)*(ldst + 1) > lm) break;
                                        else ldst++;
                                    } while (1);
                                    mdst = lm - ldst * (ldst + 1);
                                    if (mdst < 0) mdst *= -1;
                                    slm_tmp[ith].transform_M2L(mm_tmp[ith], ldst, mdst, le_tmp[ith], false);
                                    // Copy
                                    const S32 lm_loc = lm - le_k_.start_loc_[0];
                                    const S32 adr_dest = lm_loc
                                                       + le_k_.size_loc_[0] * (i_loc
                                                       + le_k_.size_loc_[1] * (j_loc
                                                       + le_k_.size_loc_[2] * k_loc));
                                    le_k_.buf_[adr_dest] = le_tmp[ith].buf[lm];
                                }
                            }
                        }
                    }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                }
#endif
                // Free memory
                for (S32 i = 0; i < n_thread; i++) {
                    slm_tmp[i].freeMem();
                    mm_tmp[i].freeMem();
                    le_tmp[i].freeMem();
                }
                delete [] slm_tmp;
                delete [] mm_tmp;
                delete [] le_tmp;
            }
        
        };
    } // END of namespace of ParticleMeshMultipole
} // END of namespace of ParticleSimulator
