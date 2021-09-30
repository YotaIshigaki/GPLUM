#pragma once
#include <complex>
#include "../ps_defs.hpp"
#include "particle_mesh_multipole_utils.hpp"
#include "cell.hpp"
#include "multidimensional_array.hpp"
#include "green_function.hpp"
#include "convolution.hpp"

namespace ParticleSimulator {
    namespace ParticleMeshMultipole {

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

        template <class Tforce, class Tepi, class Tepj, class Tspj>
        class ParticleMeshMultipole {
        public: 
            typedef double real_t;
            typedef std::complex<real_t> cplx_t;
            typedef Cell_FMM<real_t, cplx_t, Tforce, Tepi, Tepj, Tspj> Cell_t;

        private:
            bool first_call_by_initialize;
           
            // Parameters
            Parameters param_;
            Parameters param_prev_;

            // Particle information
            S32 n_loc_tot_;
            F64 msum_;
 
            // Cell information
            S32 n_cell_tot_;
            F64vec width_cell_;
            MultidimensionalArray<Cell_t,3> cell_;

            // Green function
            GreenFunction<real_t, cplx_t> gf_;

        private:

            template <class Tmomloc, class Tmomglb>
            void setParam(const TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj> & tree,
                          const DomainInfo & dinfo) {
                param_.icut = tree.getICut();
                param_.n_cell = tree.getNCell();
                param_.pos_unit_cell = tree.getPosUnitCell();
                param_.bc = dinfo.getBoundaryCondition();
                // Check consistency
                const F64ort pos_root_domain = dinfo.getPosRootDomain();
                if ((param_.bc == BOUNDARY_CONDITION_PERIODIC_XYZ) && 
                     ((param_.pos_unit_cell.low_  != pos_root_domain.low_) ||
                      (param_.pos_unit_cell.high_ != pos_root_domain.high_))) {
                    PARTICLE_SIMULATOR_PRINT_ERROR("param_ is not consistent with a given DomainInfo.")
                    Abort(-1);
                }

            }

            void setCell() {
                Comm::barrier();
                F64 time_start = GetWtime(); 

                const S32 nx = param_.n_cell.x;
                const S32 ny = param_.n_cell.y;
                const S32 nz = param_.n_cell.z;
                n_cell_tot_ = nx * ny * nz;
                width_cell_ = GetWidthOfParticleMeshCell(param_.pos_unit_cell,
                                                         param_.n_cell);
                S32 sizes[3] = {nx, ny, nz};
                cell_.initialize(sizes);
                for (S32 k=0; k<nz; k++){
                    for (S32 j=0; j<ny; j++){
                        for (S32 i=0; i<nx; i++){
                            cell_(k,j,i).init(param_.p);
                            const S32vec idx = S32vec(i,j,k);
                            const F64vec pos = GetCenterOfParticleMeshCell(param_.pos_unit_cell,
                                                                           width_cell_,
                                                                           idx);
                            cell_(k,j,i).setPos(pos);
                        }
                    }
                }

                Comm::barrier();
                F64 time_end = GetWtime(); 
                if (Comm::getRank() == 0)
                    std::cout << "setCell: " << time_end - time_start << " [s]" << std::endl;
            }

            void setGreenFunction() {

                if (!(param_ == param_prev_)) {
                    Comm::barrier();
                    F64 time_start = GetWtime();

                    gf_.free();

                    gf_.init(param_.p, 
                             param_.bc,
                             param_.n_cell.x,
                             param_.n_cell.y,
                             param_.n_cell.z);

                    Comm::barrier();
                    F64 time_end = GetWtime();
                    if (Comm::getRank() == 0)
                        std::cout << "setGreenFunction (init): " << time_end - time_start << " [s]" << std::endl;
                    time_start = time_end;

                    gf_.set(param_.icut,
                            width_cell_,
                            param_.alpha,
                            param_.nmax,
                            param_.mmax);

                    Comm::barrier();
                    time_end = GetWtime();
                    if (Comm::getRank() == 0)
                        std::cout << "setGreenFunction (set): " << time_end - time_start << " [s]" << std::endl;
                    time_start = time_end;

                    gf_.doFFT();

                    Comm::barrier();
                    time_end = GetWtime();
                    if (Comm::getRank() == 0)
                        std::cout << "setGreenFunction (doFFT): " << time_end - time_start << " [s]" << std::endl;
                }

            }

            template <class Tmomloc, class Tmomglb>
            void setParticleToCell(TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj> & tree,
                                   const DomainInfo & dinfo) {
                Comm::barrier();
                F64 time_start = GetWtime();

                n_loc_tot_ = tree.getNumberOfEpiSorted();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                for (S32 k=0; k<param_.n_cell.z; k++) {
                    for (S32 j=0; j<param_.n_cell.y; j++) {
                        for (S32 i=0; i<param_.n_cell.x; i++) {
                            S32vec idx = S32vec(i, j, k);
                            CellInfo<Tforce, Tepi, Tepj, Tspj> info;
                            if (tree.getParticleMeshCellInfo(idx, info) == 0) {
                                cell_(k,j,i).setEpjToCell(info.epj_is_shared, info.n_epj, info.epj_first);
                                cell_(k,j,i).setSpjToCell(info.n_spj, info.spj_first);
                                cell_(k,j,i).setEpiToCell(info.n_epi, info.epi_first);
                                cell_(k,j,i).setForceToCell(info.force_first);
                            }
                        }
                    }
                }

                Comm::barrier();
                F64 time_end = GetWtime();
                if (Comm::getRank() == 0)
                    std::cout << "setParticleToCell: " << time_end - time_start << " [s]" << std::endl;
               
            }

            void calcTotalCharge(const DomainInfo & dinfo) {
                Comm::barrier();
                F64 time_start = GetWtime();

                const F64ort pos_my_domain = dinfo.getPosDomain(Comm::getRank());
                F64 msum_tmp = 0.0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for reduction (+:msum_tmp)
#endif
                for (S32 i=0; i < n_cell_tot_; i++) {
                    msum_tmp += cell_[i].getTotalCharge(pos_my_domain);
#if 0
                    if (Comm::getRank() == 1)
                        std::cout << "i = " << i << " msum = " << cell_[i].getTotalCharge(pos_my_domain) << std::endl;
#endif
#if 1
                    if (Comm::getRank() == 1 && i == 4) {
                        std::cout << "msum[0] = " << cell_[i].getTotalCharge(pos_my_domain) << std::endl;
                        std::cout << "n_epj_grp = " << cell_[i].n_epj.size() << std::endl;
                        for (S32 n=0; n<cell_[i].n_epj.size(); n++) {
                            for (S32 k=0; k<cell_[i].n_epj[n]; k++) {
                                std::stringstream ss;
                                ss << std::setfill('0') << std::setw(6) << cell_[i].epj_first[n][k].getId();
                                std::cout << "id = " << ss.str()
                                          << " qj = " << cell_[i].epj_first[n][k].getCharge() << std::endl;
                            }
                        }
                        std::cout << "n_spj_grp = " << cell_[i].n_spj.size() << std::endl;
                        for (S32 n=0; n<cell_[i].n_spj.size(); n++) {
                            for (S32 k=0; k<cell_[i].n_spj[n]; k++) {
                                std::cout << "(n,k) = " << n << ", " << k 
                                          << " qj = " << cell_[i].spj_first[n][k].getCharge() << std::endl;
                            }
                        }
                    }
#endif
                }
                msum_ = Comm::getSum(msum_tmp);

                Comm::barrier();
                F64 time_end = GetWtime();
                if (Comm::getRank() == 0) {
                    std::cout << "calcTotalCharge: " << time_end - time_start << " [s]" << std::endl;
                    std::cout << "msum_ = " << msum_ << std::endl;
                }
            }

            void calcMultipoleMoments(const DomainInfo & dinfo) {
                Comm::barrier();
                F64 time_start = GetWtime();

                const F64ort pos_my_domain = dinfo.getPosDomain(Comm::getRank());
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                for(int i=0; i<n_cell_tot_; i++){
                    cell_[i].do_P2M(pos_my_domain);
                }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                const S32 LEN = (param_.p + 1)*(param_.p + 1);
                real_t *sendbuf, *recvbuf;
                sendbuf = new real_t[n_cell_tot_ * LEN];
                recvbuf = new real_t[n_cell_tot_ * LEN];
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                for (S32 i=0; i < n_cell_tot_; i++) {
                    for (S32 lm = 0; lm < LEN; lm++) {
                        const S32 k = lm + LEN * i;
                        sendbuf[k] = cell_[i].mm.buf[lm];
                    }
                }
                const S32 count = n_cell_tot_ * LEN;
                MPI_Allreduce(sendbuf, recvbuf, count, MPI_DOUBLE,
                              MPI_SUM, MPI_COMM_WORLD);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                for (S32 i = 0; i < n_cell_tot_; i++) {
                    for (S32 lm = 0; lm < LEN; lm++) {
                        const S32 k = lm + LEN * i;
                        cell_[i].mm.buf[lm] = recvbuf[k];
                    }
                }
                // Release memory
                delete [] sendbuf;
                delete [] recvbuf;
#endif // PARTICLE_SIMULATOR_MPI_PARALLEL

                Comm::barrier();
                F64 time_end = GetWtime();
                if (Comm::getRank() == 0)
                    std::cout << "calcMultipoleMoments: " << time_end - time_start << " [s]" << std::endl;
            }

            void calcForceParticleMesh(const DomainInfo & dinfo, const bool clear_force = true) {
                Comm::barrier();
                F64 time_start = GetWtime();

                M2LConvolution(gf_, cell_);

                Comm::barrier();
                F64 time_end = GetWtime();
                if (Comm::getRank() == 0) 
                    std::cout << "calcForceParticleMesh (M2L): "
                              << time_end - time_start << " [s]" << std::endl;
                time_start = time_end;

                // Dipole correction
                if (param_.bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {

                    const F64ort pos_my_domain = dinfo.getPosDomain(Comm::getRank());
                    F64 dipole_x, dipole_y, dipole_z;
                    dipole_x = dipole_y = dipole_z = 0.0;
                    F64 quad0 = 0.0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for reduction(+:dipole_x), reduction(+:dipole_y), reduction(+:dipole_z), reduction(+:quad0)
#endif
                    for (S32 i = 0; i < n_cell_tot_; i++){
                        dipole_x += cell_[i].mm.buf[3];
                        dipole_y += cell_[i].mm.buf[1];
                        dipole_z += cell_[i].mm.buf[2];
                        quad0 += cell_[i].dispersion(pos_my_domain);
                    }
                    F64vec dipole = F64vec(dipole_x, dipole_y, dipole_z);
                    quad0  = Comm::getSum(quad0);
                    const F64 pi = 4.0 * atan(1.0);
                    dipole *= (4./3.) * pi;
                    std::cout << "dipole = " << dipole << std::endl;
                    std::cout << "quad0  = " << quad0 << std::endl;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for (S32 i = 0; i < n_cell_tot_; i++){
                        cell_[i].le.buf[3] += 2.0 * dipole.x;
                        cell_[i].le.buf[1] -= 2.0 * dipole.y;
                        cell_[i].le.buf[2] += 1.0 * dipole.z;
                        cell_[i].le.buf[0] += ((2./3.) * pi) * quad0;
                        // self energy correction
                        cell_[i].le.buf[0] -= 
                            param_.alpha * (2.0/sqrt(pi)) * cell_[i].mm.buf[0];
                    }
                }

                Comm::barrier();
                time_end = GetWtime();
                if (Comm::getRank() == 0)
                    std::cout << "calcForceParticleMesh (corr.): "
                              << time_end - time_start << " [s]" << std::endl;
                time_start = time_end;

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif          
                for (S32 i = 0; i < n_cell_tot_; i++){
                    cell_[i].do_L2P(clear_force);
                    if (param_.bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                        cell_[i].do_L2P_corr(msum_, param_.alpha);
                    }
                }

                Comm::barrier();
                time_end = GetWtime();
                if (Comm::getRank() == 0)
                    std::cout << "calcForceParticleMesh (L2P): "
                              << time_end - time_start << " [s]" << std::endl;

            }


        public:

            ParticleMeshMultipole() {
                 first_call_by_initialize = true;
                 msum_ = 0.0;
            }
            ~ParticleMeshMultipole() {}

            void initialize(const S32 _p) {
                assert(first_call_by_initialize);
                first_call_by_initialize = false;
                param_.p = _p;
            }

            void reinitialize(const S32 & _p) {
                param_prev_ = param_;
                param_.p = _p;
                // Reset cell information
                cell_.freeMem();
            }

            template <class Tmomloc, class Tmomglb>
            void calcForceAll(TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj> & tree,
                              const DomainInfo & dinfo,
                              const bool clear_force = true) {

                setParam(tree, dinfo);

                setCell();

                setGreenFunction();

                setParticleToCell(tree, dinfo);

                calcTotalCharge(dinfo);

                calcMultipoleMoments(dinfo);

                calcForceParticleMesh(dinfo, clear_force);

            }

            template <class Tmomloc, class Tmomglb, class Tpsys>
            void calcForceAllAndWriteBack(TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj> & tree, 
                                          const DomainInfo & dinfo,
                                          Tpsys & psys,
                                          const bool clear_force = true) {
                calcForceAll(tree, dinfo, clear_force);
                tree.copyForceOriginalOrder();
                for(S32 i=0; i<n_loc_tot_; i++)
                    psys[i].copyFromForcePMM(tree.getForce(i));
            }

        };
    } // END of namespace of ParticleMeshMultipole
} // END of namespace of ParticleSimulator
