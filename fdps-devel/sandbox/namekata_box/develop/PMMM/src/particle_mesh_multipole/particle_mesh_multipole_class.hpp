#pragma once
#include <complex>
#include "../ps_defs.hpp"
#include "particle_mesh_multipole_utils.hpp"
#include "particle.hpp"
#include "cell.hpp"
#include "multidimensional_array.hpp"
#include "green_function.hpp"
#include "convolution.hpp"

namespace ParticleSimulator {
    namespace ParticleMeshMultipole {

        class ParticleMeshMultipole {
        public: 
            typedef double real_t;
            typedef std::complex<real_t> cplx_t;
            typedef Cell_FMM<real_t, cplx_t> Cell_t;

        private:
            bool first_call_by_initialize;
           
            // Parameters
            ParticleMeshMultipoleParameters param_;
            ParticleMeshMultipoleParameters param_prev_;

            // Particle information
            S32 n_loc_tot_;
            F64 msum_;
            Particle *ptcl_;
 
            // Cell information
            S32 n_cell_tot_;
            F64vec width_cell_;
            MultidimensionalArray<Cell_t,3> cell_;

            // Green function
            GreenFunction<real_t, cplx_t> gf_;

        private:

            void setCell() {
                Comm::barrier();
                F64 time_start = GetWtime(); 

                const S32 nx = param_.n_cell.x;
                const S32 ny = param_.n_cell.y;
                const S32 nz = param_.n_cell.z;
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

            void setParticleToCell() {
                Comm::barrier();
                F64 time_start = GetWtime();

                for (S32 i=0; i<n_loc_tot_; i++) {
                    const F64vec pos = ptcl_[i].getPos();
                    const S32vec idx = GetCellIDMeasuredInUnitCell(param_.pos_unit_cell,
                                                                   width_cell_,
                                                                   pos);
                    cell_(idx.z, idx.y, idx.x).setParticleToCell(&ptcl_[i]);
                }

                Comm::barrier();
                F64 time_end = GetWtime();
                if (Comm::getRank() == 0)
                    std::cout << "setParticleToCell: " << time_end - time_start << " [s]" << std::endl;
            }

            void calcTotalCharge() {
                Comm::barrier();
                F64 time_start = GetWtime();

                F64 msum_tmp = 0.0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for reduction (+:msum_tmp)
#endif
                for (S32 i=0; i<n_loc_tot_; i++) {
                    msum_tmp += ptcl_[i].getCharge();
                }
                msum_ = Comm::getSum(msum_tmp);

                Comm::barrier();
                F64 time_end = GetWtime();
                if (Comm::getRank() == 0) {
                    std::cout << "calcTotalCharge: " << time_end - time_start << " [s]" << std::endl;
                    std::cout << "msum_ = " << msum_ << std::endl;
                }
            }

            void calcMultipoleMoments() {
                Comm::barrier();
                F64 time_start = GetWtime();

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                for(int i=0; i<n_cell_tot_; i++){
                    cell_[i].do_P2M();
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

            void calcForceParticleMesh() {
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
                        quad0 += cell_[i].dispersion();
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
                    cell_[i].do_L2P();
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
                 ptcl_ = NULL;
            }
            ~ParticleMeshMultipole() {}

            void initialize(const ParticleMeshMultipoleParameters & _param) {
                assert(first_call_by_initialize);
                first_call_by_initialize = false;
                param_ = _param;

                n_cell_tot_ = param_.n_cell.x
                            * param_.n_cell.y
                            * param_.n_cell.z;
                width_cell_ = GetWidthOfParticleMeshCell(param_.pos_unit_cell,
                                                         param_.n_cell);
            }

            void reinitialize(const ParticleMeshMultipoleParameters & _param) {
                param_prev_ = param_;
                param_ = _param;
                // Reset cell information
                n_cell_tot_ = param_.n_cell.x
                            * param_.n_cell.y
                            * param_.n_cell.z;
                width_cell_ = GetWidthOfParticleMeshCell(param_.pos_unit_cell,
                                                         param_.n_cell);
                cell_.freeMem();
            }

            void setDomainInfo(const DomainInfo & dinfo) {
                // Currently, this function just checks if param_ passed by tree (PP part)
                // is consistent with the given DomainInfo.
                const S32 bc = dinfo.getBoundaryCondition();
                const F64ort pos_root_domain = dinfo.getPosRootDomain();
                if ((param_.bc != bc) || 
                    (bc == BOUNDARY_CONDITION_PERIODIC_XYZ && 
                     ((param_.pos_unit_cell.low_  != pos_root_domain.low_) ||
                      (param_.pos_unit_cell.high_ != pos_root_domain.high_)))) {
                    PARTICLE_SIMULATOR_PRINT_ERROR("param_ is not consistent with a given DomainInfo.")
                    Abort(-1);
                }
            }

            template <class Tpsys>
            void setParticleSystem(const Tpsys & psys,
                                   const bool clear=true) {
                 if (clear) {
                     if (ptcl_ != NULL) delete [] ptcl_;
                     
                     n_loc_tot_ = psys.getNumberOfParticleLocal();
                     ptcl_ = new Particle[n_loc_tot_];
                     for (S32 i = 0; i < n_loc_tot_; i++) {
                         ptcl_[i].pos    = psys[i].getPos();
                         ptcl_[i].charge = psys[i].getCharge();
                     }
                 } else {
                     Particle * ptcl_prev = new Particle[n_loc_tot_];
                     memcpy((void *)ptcl_prev, (void *)ptcl_, (size_t)sizeof(Particle)*n_loc_tot_);
                     delete [] ptcl_;

                     const S32 n_psys = psys.getNumberOfParticleLocal();
                     ptcl_ = new Particle[n_loc_tot_ + n_psys];
                     memcpy((void *)ptcl_, (void *)ptcl_prev, (size_t)sizeof(Particle)*n_loc_tot_);
                     for (S32 i = 0; i < n_psys; i++) {
                         const S32 ii = i + n_loc_tot_;
                         ptcl_[ii].pos  = psys[i].getPos();
                         ptcl_[ii].charge = psys[i].getCharge();
                     }
                     n_loc_tot_ += n_psys;

                     delete [] ptcl_prev;
                 }
            }


            void calcForceAll() {

                setCell();

                setGreenFunction();

                setParticleToCell();

                calcTotalCharge();

                calcMultipoleMoments();

                calcForceParticleMesh();

            }

            template <class Tpsys>
            void calcForceAllAndWriteBack(Tpsys & psys,
                                          const DomainInfo & dinfo) {
                setDomainInfo(dinfo);
                setParticleSystem(psys);
                calcForceAll();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                for (S32 i = 0; i < n_loc_tot_; i++) {
                    const F64vec acc = ptcl_[i].getAcc();
                    const F64 pot = ptcl_[i].getPot();
                    psys[i].copyFromForcePMM(acc, pot);
                }
            }

        };
    } // END of namespace of ParticleMeshMultipole
} // END of namespace of ParticleSimulator
