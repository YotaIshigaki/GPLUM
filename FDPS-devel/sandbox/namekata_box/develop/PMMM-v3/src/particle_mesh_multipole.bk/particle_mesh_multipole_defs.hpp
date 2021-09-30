#pragma once
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include <fftw3-mpi.h>
#else
#include <fftw3.h>
#endif

namespace ParticleSimulator {
    namespace ParticleMeshMultipole {

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

#if defined(PARTICLE_SIMULATOR_PMM_M2L_USE_FFTW_MEASURE)
#define FFTW_PLANNING_RIGOR_FLAG  FFTW_MEASURE
#elif defined(PARTICLE_SIMULATOR_PMM_M2L_USE_FFTW_PATIENT)
#define FFTW_PLANNING_RIGOR_FLAG  FFTW_PATIENT
#elif defined(PARTICLE_SIMULATOR_PMM_M2L_USE_FFTW_EXHAUSTIVE)
#define FFTW_PLANNING_RIGOR_FLAG  FFTW_EXHAUSTIVE
#else
#define FFTW_PLANNING_RIGOR_FLAG  FFTW_ESTIMATE
#endif

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
            F64 M2L_initialize;
            F64 M2L_gf_initialize;
            F64 M2L_gf_calc_gf_r;
            F64 M2L_gf_calc_gf_k;
            F64 M2L_preproc_mm_r_comm; // redistMM
            F64 M2L_redist_mm_r; // redistMM
            F64 M2L_postproc_mm_r_comm; // redistMM
            F64 M2L_mm_r_to_mm_k;
            F64 M2L_collect_mm_k;
            F64 M2L_transform;
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
                ret.M2L_initialize = this->M2L_initialize 
                                   + rhs.M2L_initialize;
                ret.M2L_gf_initialize = this->M2L_gf_initialize 
                                      + rhs.M2L_gf_initialize;
                ret.M2L_gf_calc_gf_r = this->M2L_gf_calc_gf_r 
                                     + rhs.M2L_gf_calc_gf_r;
                ret.M2L_gf_calc_gf_k = this->M2L_gf_calc_gf_k 
                                     + rhs.M2L_gf_calc_gf_k;
                ret.M2L_preproc_mm_r_comm = this->M2L_preproc_mm_r_comm 
                                          + rhs.M2L_preproc_mm_r_comm;
                ret.M2L_redist_mm_r = this->M2L_redist_mm_r 
                                    + rhs.M2L_redist_mm_r;
                ret.M2L_postproc_mm_r_comm = this->M2L_postproc_mm_r_comm 
                                           + rhs.M2L_postproc_mm_r_comm;
                ret.M2L_mm_r_to_mm_k = this->M2L_mm_r_to_mm_k 
                                     + rhs.M2L_mm_r_to_mm_k;
                ret.M2L_collect_mm_k = this->M2L_collect_mm_k 
                                     + rhs.M2L_collect_mm_k;
                ret.M2L_transform = this->M2L_transform 
                                  + rhs.M2L_transform;
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
                ret.M2L_initialize = M2L_initialize / s;
                ret.M2L_gf_initialize = M2L_gf_initialize / s ;
                ret.M2L_gf_calc_gf_r = M2L_gf_calc_gf_r / s;
                ret.M2L_gf_calc_gf_k = M2L_gf_calc_gf_k / s;
                ret.M2L_preproc_mm_r_comm = M2L_preproc_mm_r_comm / s;
                ret.M2L_redist_mm_r = M2L_redist_mm_r / s;
                ret.M2L_postproc_mm_r_comm = M2L_postproc_mm_r_comm / s;
                ret.M2L_mm_r_to_mm_k = M2L_mm_r_to_mm_k / s;
                ret.M2L_collect_mm_k = M2L_collect_mm_k / s;
                ret.M2L_transform = M2L_transform / s;
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
                M2L_initialize = 0.0;
                M2L_gf_initialize = 0.0;;
                M2L_gf_calc_gf_r = 0.0;;
                M2L_gf_calc_gf_k = 0.0;;
                M2L_preproc_mm_r_comm = 0.0;
                M2L_redist_mm_r = 0.0;
                M2L_postproc_mm_r_comm = 0.0;
                M2L_mm_r_to_mm_k = 0.0;
                M2L_collect_mm_k = 0.0;
                M2L_transform = 0.0;
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
                fout << "M2L_initialize = " << M2L_initialize << std::endl;
                fout << "M2L_gf_initialize = " << M2L_gf_initialize << std::endl;
                fout << "M2L_gf_calc_gf_r = " << M2L_gf_calc_gf_r << std::endl;
                fout << "M2L_gf_calc_gf_k = " << M2L_gf_calc_gf_k << std::endl;
                fout << "M2L_preproc_mm_r_comm = " << M2L_preproc_mm_r_comm << std::endl;
                fout << "M2L_redist_mm_r = " << M2L_redist_mm_r << std::endl;
                fout << "M2L_postproc_mm_r_comm = " << M2L_postproc_mm_r_comm << std::endl;
                fout << "M2L_mm_r_to_mm_k = " << M2L_mm_r_to_mm_k << std::endl;
                fout << "M2L_collect_mm_k = " << M2L_collect_mm_k << std::endl;
                fout << "M2L_transform = " << M2L_transform << std::endl;
                fout << "M2L_le_k_to_le_r = " << M2L_le_k_to_le_r << std::endl;
                fout << "M2L_preproc_le_r_comm = " << M2L_preproc_le_r_comm << std::endl;
                fout << "M2L_redist_le_r = " << M2L_redist_le_r << std::endl;
                fout << "M2L_postproc_le_r_comm = " << M2L_postproc_le_r_comm << std::endl;
                fout << "L2P = " << L2P << std::endl;
            }
        };


    } // END of namespace of ParticleMeshMultipole
} // END of namespace of ParticleSimulator
