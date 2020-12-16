#pragma once
#include <cassert>
#include <map>
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include <fftw3-mpi.h>
#else
#include <fftw3.h>
#endif

namespace ParticleSimulator {
    namespace ParticleMeshMultipole {
        // Data type
#if defined(PARTICLE_SIMULATOR_PMM_FFTW_FLOAT)
        using fftw_real_t = float;
        using fftw_cplx_t = fftwf_complex;
#define fftw_plan            fftwf_plan
#define fftw_plan_dft_r2c_3d fftwf_plan_dft_r2c_3d
#define fftw_plan_dft_c2r_3d fftwf_plan_dft_c2r_3d 
#define fftw_execute         fftwf_execute
#define fftw_destroy_plan    fftwf_destroy_plan
#elif defined(PARTICLE_SIMULATOR_PMM_FFTW_LONG_DOUBLE)
        using fftw_real_t = long double;
        using fftw_cplx_t = fftwl_complex;
#define fftw_plan            fftwl_plan
#define fftw_plan_dft_r2c_3d fftwl_plan_dft_r2c_3d
#define fftw_plan_dft_c2r_3d fftwl_plan_dft_c2r_3d 
#define fftw_execute         fftwl_execute
#define fftw_destroy_plan    fftwl_destroy_plan
#else
        using fftw_real_t = double;
        using fftw_cplx_t = fftw_complex;
#endif

        // Auxiliary variables to manage call of FFTW functions
        bool first_call_of_fftw_mpi_init {true};
        bool first_call_of_fftw_init_threads {true};


        class Parameters {
        public:
            S32 p; // the maximum order of multipole expansion
            S32 icut; // the minimum cell separation 
            S32 bc; // the boundary condition 
            F64ort pos_unit_cell; // the size of unit cell (real image of the computational box).
            S32vec n_cell; // the number of particle mesh cells
            S32 nmax; // this parameter determines how far virtual image boxes
                      // are to be considered in the calculation of the Ewald
                      // summation in the real space.
            S32 mmax; // this parameter determines how large wavenumber vectors
                      // are to be considered in the calculation of the Ewald
                      // summation in the wavenumber space.
            F64 alpha; // A parameter of Ewald summation. 

            Parameters() {
                p = FDPS_DFLT_VAL_PMMM_P;
                icut = FDPS_DFLT_VAL_PMMM_ICUT;
                bc = BOUNDARY_CONDITION_PERIODIC_XYZ;
                pos_unit_cell.init();
                n_cell = S32vec(FDPS_DFLT_VAL_PMMM_N_CELL_1D,
                                FDPS_DFLT_VAL_PMMM_N_CELL_1D,
                                FDPS_DFLT_VAL_PMMM_N_CELL_1D);
                nmax = FDPS_DFLT_VAL_PMMM_NMAX;
                mmax = FDPS_DFLT_VAL_PMMM_MMAX;
                alpha = FDPS_DFLT_VAL_PMMM_ALPHA;
            }

            const Parameters & operator = (const Parameters & rhs) {
                p             = rhs.p;
                icut          = rhs.icut;
                bc            = rhs.bc;
                pos_unit_cell = rhs.pos_unit_cell; 
                n_cell        = rhs.n_cell;
                nmax          = rhs.nmax;
                mmax          = rhs.mmax;
                alpha         = rhs.alpha;
                return (*this);
            }

            bool operator != (const Parameters & rhs) const {
                return ((p != rhs.p) ||
                        (icut != rhs.icut) ||
                        (bc != rhs.bc) ||
                        (pos_unit_cell.low_ != rhs.pos_unit_cell.low_) ||
                        (pos_unit_cell.high_ != rhs.pos_unit_cell.high_) ||
                        (n_cell != rhs.n_cell) ||
                        (nmax != rhs.nmax) ||
                        (mmax != rhs.mmax) ||
                        (alpha != rhs.alpha));
            }

            bool operator == (const Parameters & rhs) const {
                return ((p == rhs.p) &&
                        (icut == rhs.icut) && 
                        (bc == rhs.bc) &&
                        (pos_unit_cell.low_ == rhs.pos_unit_cell.low_) &&
                        (pos_unit_cell.high_ == rhs.pos_unit_cell.high_) &&
                        (n_cell == rhs.n_cell) &&
                        (nmax == rhs.nmax) &&
                        (mmax == rhs.mmax) &&
                        (alpha == rhs.alpha));
            }

            S32 sanity_check(std::ostream & fout=std::cout) const {
                S32 ret = 0;
                if (p <= 0) {
                    ret = 1;
                    PARTICLE_SIMULATOR_PRINT_ERROR("The order of multipole expansion must be > 0");
                    fout << "p = " << p << std::endl;
                }
                if (icut < 1) {
                    ret = 1;
                    PARTICLE_SIMULATOR_PRINT_ERROR("The cell separation must be >= 1");
                    fout << "icut = " << icut << std::endl;
                }
                if ((bc != BOUNDARY_CONDITION_OPEN) && (bc != BOUNDARY_CONDITION_PERIODIC_XYZ)) {
                    ret = 1;
                    PARTICLE_SIMULATOR_PRINT_ERROR("This boundary condition is not supported");
                    fout << "bc = " << bc << std::endl;
                }
                if (!pos_unit_cell.isValid()) {
                    ret = 1;
                    PARTICLE_SIMULATOR_PRINT_ERROR("The size of unit cell is invalid");
                    fout << "pos_unit_cell = " << pos_unit_cell << std::endl;
                }
                if ((n_cell.x <= 0) || (n_cell.y <= 0) || (n_cell.z <= 0)) {
                    ret = 1;
                    PARTICLE_SIMULATOR_PRINT_ERROR("The number of the cells in each dimension must be > 0");
                    fout << "n_cell = " << n_cell << std::endl;
                }
                if ((nmax <= 0) || (mmax <= 0) || (alpha <= 0.0)) {
                    ret = 1;
                    PARTICLE_SIMULATOR_PRINT_ERROR("The parameters for Ewald summation method are invalid");
                    fout << "nmax  = " << nmax << std::endl;
                    fout << "mmax  = " << mmax << std::endl;
                    fout << "alpha = " << alpha << std::endl;
                }
                return ret;
            }

            void dump(std::ostream & fout=std::cout) const {
                fout << "------------------------------------------------" << std::endl;
                fout << "Parameters of Particle Mesh Multipole method"     << std::endl;
                fout << "    p      = " << p << std::endl;
                fout << "    icut   = " << icut << std::endl;
                fout << "    bc     = " << bc << std::endl;
                fout << "    n_cell = " << n_cell << std::endl;
                fout << "    nmax   = " << nmax << std::endl;
                fout << "    mmax   = " << mmax << std::endl;
                fout << "    alpha  = " << alpha << std::endl;
                fout << "------------------------------------------------" << std::endl;
            }

        };

        template <class T, int DIM = 4>
        class Buffer {
        public:
            S32 size_glb_tot_;
            S32 size_glb_[DIM];
            S32 start_glb_[DIM], end_glb_[DIM];
            S32 size_loc_tot_;
            S32 size_loc_[DIM];
            S32 start_loc_[DIM], end_loc_[DIM];
            S32 capacity_;
            T *buf_; // local buffer

            Buffer() {
                size_glb_tot_ = 0;
                size_loc_tot_ = 0;
                for (S32 i = 0; i < DIM; i++) {
                    size_glb_[i] = 0;
                    start_glb_[i] = -1;
                    end_glb_[i] = -1;
                    size_loc_[i] = 0;
                    start_loc_[i] = -1;
                    end_loc_[i] = -1;
                }
                capacity_ = 0;
                buf_ = nullptr;
            }

            ~Buffer() {
                free();
            }

            void free() {
                if (buf_ != nullptr) {
                    delete [] buf_;
                    buf_ = nullptr;
                    capacity_ = 0;
                }
            }

            void clear() {
                for (S32 i = 0; i < capacity_; i++) buf_[i] = 0;
            }

            void resize() {
                if (size_loc_tot_ > 0) {
                    if (capacity_ == 0) {
                        assert(buf_ == nullptr);
                        buf_ = new T [size_loc_tot_];
                        capacity_ = size_loc_tot_;
                    } else {
                        assert(capacity_ > 0);
                        if (capacity_ < size_loc_tot_) {
                            assert(buf_ != nullptr);
                            delete [] buf_;
                            buf_ = new T [size_loc_tot_];
                            capacity_ = size_loc_tot_;
                        }
                    }
                }
            }

            void moveFrom(Buffer & src) {
                // The following conditions must be fulfilled.
                assert((buf_ == nullptr) && (src.buf_ != nullptr));
                assert(size_loc_tot_ <= src.capacity_);

                capacity_ = src.capacity_;
                buf_      = src.buf_;
                src.capacity_ = 0;
                src.buf_ = nullptr;
            }

            void copySizeInfoFrom(const Buffer & tmp) {
                size_glb_tot_ = tmp.size_glb_tot_;
                size_loc_tot_ = tmp.size_loc_tot_;
                for (S32 i = 0; i < DIM; i++) {
                    size_glb_[i] = tmp.size_glb_[i];
                    start_glb_[i] = tmp.start_glb_[i];
                    end_glb_[i] = tmp.end_glb_[i];
                    size_loc_[i] = tmp.size_loc_[i];
                    start_loc_[i] = tmp.start_loc_[i];
                    end_loc_[i] = tmp.end_loc_[i];
                }
            }

            void copySizeInfoGlbOnlyFrom(const Buffer & tmp) {
                size_glb_tot_ = tmp.size_glb_tot_;
                for (S32 i = 0; i < DIM; i++) {
                    size_glb_[i] = tmp.size_glb_[i];
                    start_glb_[i] = tmp.start_glb_[i];
                    end_glb_[i] = tmp.end_glb_[i];
                }
            }

            void calcSizeGlb() {
                for (S32 i = 0; i < DIM; i++) {
                    size_glb_[i] = end_glb_[i] - start_glb_[i] + 1;
                }
            }

            void calcSizeLoc() {
                for (S32 i = 0; i < DIM; i++) {
                    size_loc_[i] = end_loc_[i] - start_loc_[i] + 1;
                }
            }

            void calcSizeTot() {
                size_glb_tot_ = 1;
                size_loc_tot_ = 1;
                for (S32 i = 0; i < DIM; i++) {
                    size_glb_tot_ *= size_glb_[i];
                    size_loc_tot_ *= size_loc_[i];
                }
            }

            void outputSizeInfo(std::ostream & fout=std::cout) const {
                fout << "size_glb_tot_ = " << size_glb_tot_ << std::endl;
                for (S32 i = 0; i < DIM; i++) {
                    fout << start_glb_[i] << "   "
                         << end_glb_[i] << "   "
                         << size_glb_[i] << std::endl;
                }
                fout << "size_loc_tot_ = " << size_loc_tot_ << std::endl;
                for (S32 i = 0; i < DIM; i++) {
                    fout << start_loc_[i] << "   "
                         << end_loc_[i] << "   " 
                         << size_loc_[i] << std::endl;
                }
            }

            void writeBufferToFile(const std::string & file_name_prefix,
                                   const S32 file_number,
                                   const bool skip_zero = false,
                                   const bool transposed_layout = false) const {
                // Set file name
                std::stringstream ss;
                ss << file_name_prefix << std::setfill('0') << std::setw(5) << file_number << ".txt";
                const std::string filename = ss.str();
                std::ofstream output_file;
                output_file.open(filename.c_str(), std::ios::trunc);
                if (transposed_layout) {
                    for (S32 j = start_loc_[3]; j <= end_loc_[3]; j++)
                    for (S32 k = start_loc_[2]; k <= end_loc_[2]; k++)
                    for (S32 i = start_loc_[1]; i <= end_loc_[1]; i++)
                    for (S32 lm = start_loc_[0]; lm <= end_loc_[0]; lm++)
                    { // for(k,j,i,lm)
                        const S32 idx = lm 
                                      + size_glb_[0] * (i
                                      + size_glb_[1] * (k
                                      + size_glb_[2] * j));
                        const S32 lm_loc = lm - start_loc_[0];
                        const S32 i_loc  = i  - start_loc_[1];
                        const S32 k_loc  = k  - start_loc_[2];
                        const S32 j_loc  = j  - start_loc_[3];
                        const S32 adr = lm_loc
                                      + size_loc_[0] * (i_loc
                                      + size_loc_[1] * (k_loc
                                      + size_loc_[2] * j_loc));
                        if (skip_zero && buf_[adr] == T(0)) continue;
                        output_file << idx << "    " <<buf_[adr] << std::endl;
                    } // for(k,j,i,lm)
                } else {
                    for (S32 k = start_loc_[3]; k <= end_loc_[3]; k++)
                    for (S32 j = start_loc_[2]; j <= end_loc_[2]; j++)
                    for (S32 i = start_loc_[1]; i <= end_loc_[1]; i++)
                    for (S32 lm = start_loc_[0]; lm <= end_loc_[0]; lm++)
                    { // for(k,j,i,lm)
                        const S32 idx = lm 
                                      + size_glb_[0] * (i
                                      + size_glb_[1] * (j
                                      + size_glb_[2] * k));
                        const S32 lm_loc = lm - start_loc_[0];
                        const S32 i_loc  = i  - start_loc_[1];
                        const S32 j_loc  = j  - start_loc_[2];
                        const S32 k_loc  = k  - start_loc_[3];
                        const S32 adr = lm_loc
                                      + size_loc_[0] * (i_loc
                                      + size_loc_[1] * (j_loc
                                      + size_loc_[2] * k_loc));
                        if (skip_zero && buf_[adr] == T(0)) continue;
                        output_file << idx << "    " <<buf_[adr] << std::endl;
                    } // for(k,j,i,lm)
                }
                output_file.close();
            }

        };

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        template <class T>
        class CommBuffer {
        public:
            S32 n_comm;
            S32 count_tot;
            S32 *ranks, *counts, *displs;
            std::unordered_map<S32, S32> adr_from_rank;
            T *buf;

            CommBuffer(): 
                n_comm {0}, count_tot {0}, ranks {nullptr},
                counts {nullptr}, displs {nullptr}, buf {nullptr} {}

            ~CommBuffer() {
                if (ranks  != nullptr) delete [] ranks;
                if (counts != nullptr) delete [] counts;
                if (displs != nullptr) delete [] displs;
                if (buf    != nullptr) delete [] buf;
            }

            void allocCommInfo() {
                if (n_comm > 0) {
                    ranks  = new S32[n_comm];
                    counts = new S32[n_comm];
                    displs = new S32[n_comm];
                }
            }

            void allocBuffer() {
                if (count_tot > 0) buf = new T[count_tot];
            }

            void calcDispls() {
                if (n_comm > 0) {
                    displs[0] = 0;
                    for (S32 i = 1; i < n_comm; i++) 
                        displs[i] = displs[i-1] + counts[i-1];
                }
            }

            void clearCounts() {
                if (n_comm > 0) for (S32 i = 0; i < n_comm; i++) counts[i] = 0;
            }

            void dumpCommInfo(std::ostream & fout=std::cout) const {
                fout << "n_comm = " << n_comm << std::endl;
                fout << "count_tot = " << count_tot << std::endl;
                for (S32 i = 0; i < n_comm; i++) {
                    fout << "i = " << i
                         << "    rank = " << ranks[i]
                         << "    count = " << counts[i]
                         << "    displ = " << displs[i]
                         << std::endl;
                }
            }

        };


        template <class send_t, class recv_t>
        void performComm(const CommBuffer<send_t> & send,
                         CommBuffer<recv_t> & recv,
                         const MPI_Comm comm = MPI_COMM_WORLD,
                         const bool debug_flag = false) {
            // Set up req[]
            S32 n_req_used = 0;
            MPI_Request *req;
            if (recv.n_comm > 0) req = new MPI_Request[recv.n_comm];
            else req = nullptr;
            // MPI_Irecv
            MPI_Datatype mpi_recv_t = GetDataType<recv_t>();
            for (S32 i = 0; i < recv.n_comm; i++) {
                const S32 rnk = recv.ranks[i];
                const S32 count = recv.counts[i];
                const S32 disp = recv.displs[i];
                const S32 tag = 0;
                if (count > 0) {
                    MPI_Irecv(&recv.buf[disp], count, mpi_recv_t,
                              rnk, tag, comm, &req[n_req_used++]);
                }
            }
            // MPI_Send
            MPI_Datatype mpi_send_t = GetDataType<send_t>();
            for (S32 i = 0; i < send.n_comm; i++) {
                const S32 rnk = send.ranks[i];
                const S32 count = send.counts[i];
                const S32 disp = send.displs[i];
                const S32 tag = 0;
                if (count > 0) {
                    MPI_Send(&send.buf[disp], count, mpi_send_t,
                             rnk, tag, comm);
                }
            }
            if (n_req_used > 0) {
                MPI_Status *stat = new MPI_Status[n_req_used];
                MPI_Waitall(n_req_used, req, stat);
                delete [] stat;
            }
            if (recv.n_comm > 0) delete [] req;
        }
#endif // PARTICLE_SIMULATOR_MPI_PARALLEL

        
        class TimeProfilePMM {
        public:
            F64 calc_force_all_and_write_back;

            // To check breakdown
            F64 initialize;
            F64 set_cell;
            F64 set_ip_info_to_cell;
            F64 calc_msum_and_quad0;
            F64 calc_multipole_moment;
            F64 collect_multipole_moment;
            F64 GFC_initialize; // GFC := Green Funtion Calculator
            F64 GFC_calc_gf_r;
            F64 GFC_redist_gf_r;
            F64 GFC_calc_gf_k;
            F64 GFC_redist_gf_k;
            F64 M2L_initialize;
            F64 M2L_preproc_mm_r_comm; // redistMM
            F64 M2L_redist_mm_r; // redistMM
            F64 M2L_postproc_mm_r_comm; // redistMM
            F64 M2L_mm_r_to_mm_k;
            F64 M2L_gather_mm_k_trans;
            F64 M2L_transform;
            F64 M2L_scatter_le_k_trans;
            F64 M2L_le_k_to_le_r;
            F64 M2L_preproc_le_r_comm; // redstLE
            F64 M2L_redist_le_r; // redistLE
            F64 M2L_postproc_le_r_comm; // redistLE
            F64 L2P;

            TimeProfilePMM() {
                clear();
            }

            TimeProfilePMM operator + (const TimeProfilePMM & rhs) const {
                TimeProfilePMM ret;
                ret.calc_force_all_and_write_back = this->calc_force_all_and_write_back
                                                  + rhs.calc_force_all_and_write_back;

                ret.initialize = this->initialize 
                               + rhs.initialize;
                ret.set_cell = this->set_cell 
                             + rhs.set_cell;
                ret.set_ip_info_to_cell = this->set_ip_info_to_cell 
                                        + rhs.set_ip_info_to_cell;
                ret.calc_msum_and_quad0 = this->calc_msum_and_quad0 
                                        + rhs.calc_msum_and_quad0;
                ret.calc_multipole_moment = this->calc_multipole_moment 
                                          + rhs.calc_multipole_moment;
                ret.collect_multipole_moment = this->collect_multipole_moment 
                                             + rhs.collect_multipole_moment;
                ret.GFC_initialize = this->GFC_initialize
                                   + rhs.GFC_initialize;
                ret.GFC_calc_gf_r = this->GFC_calc_gf_r 
                                  + rhs.GFC_calc_gf_r;
                ret.GFC_redist_gf_r = this->GFC_redist_gf_r 
                                    + rhs.GFC_redist_gf_r;
                ret.GFC_calc_gf_k = this->GFC_calc_gf_k 
                                  + rhs.GFC_calc_gf_k;
                ret.GFC_redist_gf_k = this->GFC_redist_gf_k
                                    + rhs.GFC_redist_gf_k;
                ret.M2L_initialize = this->M2L_initialize 
                                   + rhs.M2L_initialize;
                ret.M2L_preproc_mm_r_comm = this->M2L_preproc_mm_r_comm 
                                          + rhs.M2L_preproc_mm_r_comm;
                ret.M2L_redist_mm_r = this->M2L_redist_mm_r 
                                    + rhs.M2L_redist_mm_r;
                ret.M2L_postproc_mm_r_comm = this->M2L_postproc_mm_r_comm 
                                           + rhs.M2L_postproc_mm_r_comm;
                ret.M2L_mm_r_to_mm_k = this->M2L_mm_r_to_mm_k 
                                     + rhs.M2L_mm_r_to_mm_k;
                ret.M2L_gather_mm_k_trans = this->M2L_gather_mm_k_trans
                                          + rhs.M2L_gather_mm_k_trans;
                ret.M2L_transform = this->M2L_transform 
                                  + rhs.M2L_transform;
                ret.M2L_scatter_le_k_trans = this->M2L_scatter_le_k_trans
                                           + rhs.M2L_scatter_le_k_trans;
                ret.M2L_le_k_to_le_r = this->M2L_le_k_to_le_r 
                                     + rhs.M2L_le_k_to_le_r;
                ret.M2L_preproc_le_r_comm = this->M2L_preproc_le_r_comm 
                                          + rhs.M2L_preproc_le_r_comm;
                ret.M2L_redist_le_r = this->M2L_redist_le_r 
                                    + rhs.M2L_redist_le_r;
                ret.M2L_postproc_le_r_comm = this->M2L_postproc_le_r_comm 
                                           + rhs.M2L_postproc_le_r_comm;
                ret.L2P = this->L2P 
                        + rhs.L2P;

                return ret;
            }

            TimeProfilePMM operator / (const F64 s) const {
                TimeProfilePMM ret;
                ret.calc_force_all_and_write_back = calc_force_all_and_write_back / s;

                ret.initialize = initialize / s;
                ret.set_cell = set_cell / s ;
                ret.set_ip_info_to_cell = set_ip_info_to_cell / s;
                ret.calc_msum_and_quad0 = calc_msum_and_quad0 / s;
                ret.calc_multipole_moment = calc_multipole_moment / s;
                ret.collect_multipole_moment = collect_multipole_moment / s;
                ret.GFC_initialize = GFC_initialize / s ;
                ret.GFC_calc_gf_r = GFC_calc_gf_r / s;
                ret.GFC_redist_gf_r = GFC_redist_gf_r / s;
                ret.GFC_calc_gf_k = GFC_calc_gf_k / s;
                ret.GFC_redist_gf_k = GFC_redist_gf_k / s;
                ret.M2L_initialize = M2L_initialize / s;
                ret.M2L_preproc_mm_r_comm = M2L_preproc_mm_r_comm / s;
                ret.M2L_redist_mm_r = M2L_redist_mm_r / s;
                ret.M2L_postproc_mm_r_comm = M2L_postproc_mm_r_comm / s;
                ret.M2L_mm_r_to_mm_k = M2L_mm_r_to_mm_k / s;
                ret.M2L_gather_mm_k_trans = M2L_gather_mm_k_trans / s;
                ret.M2L_transform = M2L_transform / s;
                ret.M2L_scatter_le_k_trans = M2L_scatter_le_k_trans / s;
                ret.M2L_le_k_to_le_r = M2L_le_k_to_le_r / s;
                ret.M2L_preproc_le_r_comm = M2L_preproc_le_r_comm / s;
                ret.M2L_redist_le_r = M2L_redist_le_r / s;
                ret.M2L_postproc_le_r_comm = M2L_postproc_le_r_comm / s;
                ret.L2P = L2P / s;

                return ret;
            }

            const TimeProfilePMM & operator /= (const F64 s) {
                (*this) = (*this) / s;
                return (*this);
            }


            void clear() {
                calc_force_all_and_write_back = 0.0;

                initialize = 0.0;
                set_cell = 0.0;
                set_ip_info_to_cell = 0.0;
                calc_msum_and_quad0 = 0.0;
                calc_multipole_moment = 0.0;
                collect_multipole_moment = 0.0;
                GFC_initialize = 0.0;;
                GFC_calc_gf_r = 0.0;;
                GFC_redist_gf_r = 0.0;;
                GFC_calc_gf_k = 0.0;;
                GFC_redist_gf_k = 0.0;;
                M2L_initialize = 0.0;
                M2L_preproc_mm_r_comm = 0.0;
                M2L_redist_mm_r = 0.0;
                M2L_postproc_mm_r_comm = 0.0;
                M2L_mm_r_to_mm_k = 0.0;
                M2L_gather_mm_k_trans = 0.0;
                M2L_transform = 0.0;
                M2L_scatter_le_k_trans = 0.0;
                M2L_le_k_to_le_r = 0.0;
                M2L_preproc_le_r_comm = 0.0;
                M2L_redist_le_r = 0.0;
                M2L_postproc_le_r_comm = 0.0;
                L2P = 0.0;
            }

            void dump(std::ostream & fout=std::cout) const {
                fout << "calc_force_all_and_write_back = " 
                     << calc_force_all_and_write_back << std::endl;

                fout << "initialize = " << initialize << std::endl;
                fout << "set_cell = " << set_cell << std::endl;
                fout << "set_ip_info_to_cell = " << set_ip_info_to_cell << std::endl;
                fout << "calc_msum_and_quad0 = " << calc_msum_and_quad0 << std::endl;
                fout << "calc_multipole_moment = " << calc_multipole_moment << std::endl;
                fout << "collect_multipole_moment = " << collect_multipole_moment << std::endl;
                fout << "GFC_initialize = " << GFC_initialize << std::endl;
                fout << "GFC_calc_gf_r = " << GFC_calc_gf_r << std::endl;
                fout << "GFC_redist_gf_r = " << GFC_redist_gf_r << std::endl;
                fout << "GFC_calc_gf_k = " << GFC_calc_gf_k << std::endl;
                fout << "GFC_redist_gf_k = " << GFC_redist_gf_k << std::endl;
                fout << "M2L_initialize = " << M2L_initialize << std::endl;
                fout << "M2L_preproc_mm_r_comm = " << M2L_preproc_mm_r_comm << std::endl;
                fout << "M2L_redist_mm_r = " << M2L_redist_mm_r << std::endl;
                fout << "M2L_postproc_mm_r_comm = " << M2L_postproc_mm_r_comm << std::endl;
                fout << "M2L_mm_r_to_mm_k = " << M2L_mm_r_to_mm_k << std::endl;
                fout << "M2L_gather_mm_k_trans = " << M2L_gather_mm_k_trans << std::endl;
                fout << "M2L_transform = " << M2L_transform << std::endl;
                fout << "M2L_scatter_le_k_trans = " << M2L_scatter_le_k_trans << std::endl;
                fout << "M2L_le_k_to_le_r = " << M2L_le_k_to_le_r << std::endl;
                fout << "M2L_preproc_le_r_comm = " << M2L_preproc_le_r_comm << std::endl;
                fout << "M2L_redist_le_r = " << M2L_redist_le_r << std::endl;
                fout << "M2L_postproc_le_r_comm = " << M2L_postproc_le_r_comm << std::endl;
                fout << "L2P = " << L2P << std::endl;
            }
        };


        std::map<S32, S32> primeFactorization(const S32 val) {
            assert(val > 1);
            std::map<S32, S32> ret;
            S32 tmp = val;
            S32 prime_number = 2;
            do {
                // Repeat a division until tmp cannot be divided
                // by the current prime number
                do {
                    if (tmp % prime_number == 0) {
                       tmp /= prime_number;
                       ret[prime_number]++;
                    } else break;
                } while (1);
                // Update prime number
                if (tmp > 1) {
                    prime_number++;
                    if (prime_number == tmp) {
                        ret[prime_number]++;
                        break;
                    }
                } else break;
            } while (1);
            return ret;
        }

        S32 getMaximumFactorLessThanOrEqualTo(const S32 val, const S32 UL) {
            // Perform prime factorization
            std::map<S32, S32> mp = primeFactorization(val);
            // Copy information of prime factorization
            std::vector<S32> numbers, exponents;
            for (auto itr = mp.begin(); itr != mp.end(); ++itr) {
                numbers.push_back(itr->first);
                exponents.push_back(itr->second);
            }
            // Calculate # of combinations
            S32 n_combo = 1;
            for (S32 i = 0; i < exponents.size(); i++)
                n_combo *= (exponents[i] + 1); // +1 comes from the case of 0.
            // Calculate strides of combination space
            std::vector<S32> strides;
            strides.resize(exponents.size());
            strides[0] = 1;
            for (S32 i = 1; i < exponents.size(); i++) {
                strides[i] = strides[i-1] * (exponents[i-1] + 1);
            }
            // Find the maximum factor
            S32 fact_max = -1;
            for (S32 n = 0; n < n_combo; n++) { // n is combination ID
                std::vector<S32> powers;
                powers.resize(exponents.size());
                S32 tmp = n;
                for (S32 i = powers.size() - 1; i >= 0; i--) {
                    powers[i] = tmp / strides[i];
                    tmp -= strides[i] * powers[i];
                }
                assert(tmp == 0);
                // Calculate a factor corresponding to powers
                S32 fact = 1;
                for (S32 i = 0; i < numbers.size(); i++) {
                    S32 fact_tmp = 1;
                    for (S32 k = 0; k < powers[i]; k++) fact_tmp *= numbers[i];
                    fact *= fact_tmp;
                }
                // Compare
                if (fact == UL) {
                    return fact;
                } else if (fact < UL) {
                    if (fact > fact_max) fact_max = fact;
                }
            }
            return fact_max;
        }
    } // END of namespace of ParticleMeshMultipole
} // END of namespace of ParticleSimulator
