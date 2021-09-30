#pragma once
#include <cmath>
#include <math.h>
#include <complex>
#include "ps_defs.hpp"

namespace ParticleSimulator {
    namespace Ewald {
        using real_t = double;
        using cplx_t = std::complex<real_t>;
        
        static inline cplx_t expi(real_t theta){ // exp(i\theta) 
            return cplx_t(std::cos(theta), std::sin(theta));
        }
        
        class Waves {
        public:
            S32 MMAX_;
            cplx_t *wx_;
            cplx_t *wy_;
            cplx_t *wz_;

            Waves(): wx_(nullptr), wy_(nullptr), wz_(nullptr) {}
            ~Waves() {
                if (wx_ != nullptr) delete [] wx_;
                if (wy_ != nullptr) delete [] wy_;
                if (wz_ != nullptr) delete [] wz_;
            }
           
            void initialize(const S32 MMAX,
                            const F64vec & period_length,
                            const F64vec & pos) {
                MMAX_ = MMAX;
                wx_ = new cplx_t [2*MMAX+1];
                wy_ = new cplx_t [2*MMAX+1];
                wz_ = new cplx_t [2*MMAX+1];

                const real_t twopi = 8.0 * std::atan(1.0); // 2\pi
                const cplx_t ex = expi(twopi * pos.x / period_length.x);
                const cplx_t ey = expi(twopi * pos.y / period_length.y);
                const cplx_t ez = expi(twopi * pos.z / period_length.z);
                
                cplx_t emx(1.0, 0.0);
                cplx_t emy(1.0, 0.0);
                cplx_t emz(1.0, 0.0);
                
                for(S32 m=0; m<=MMAX_; m++){
                    wx_[MMAX_ + m] = emx;
                    wx_[MMAX_ - m] = std::conj(emx);
                    emx *= ex;
                    wy_[MMAX_ + m] = emy;
                    wy_[MMAX_ - m] = std::conj(emy);
                    emy *= ey;
                    wz_[MMAX_ + m] = emz;
                    wz_[MMAX_ - m] = std::conj(emz);
                    emz *= ez;
                }
            }
        };
        
        class EwaldSum {
        public:
            using value_type = cplx_t;
            bool is_initialized_;
            S32 MMAX_;
            S32 NX_, NY_, NZ_;
            F64vec period_length_;
            cplx_t * sum_;

            EwaldSum(): is_initialized_(false), sum_(nullptr) {}
            ~EwaldSum() {
                if (sum_ != nullptr) delete [] sum_;
            }

            void initialize(const S32 MMAX,
                            const F64vec & period_length) {
                if (is_initialized_) {
                    if (sum_ != nullptr) delete [] sum_;
                }
                MMAX_ = MMAX;
                NX_ = 2*MMAX+1;
                NY_ = 2*MMAX+1;
                NZ_ = 2*MMAX+1;
                period_length_ = period_length;
                sum_ = new cplx_t[NX_*NY_*NZ_];
                is_initialized_ = true;
            }
            
            S32 get_index(const S32 ix, const S32 iy, const S32 iz) const {
                return ix + iy * NX_ + iz * NX_ * NY_;
            }

            S32 size() const {
                return NX_ * NY_ * NZ_;
            }

            cplx_t * getPointer(const S32 i=0) const { return sum_+i; }

            void clear(){
                assert(is_initialized_);
                for (S32 i=0; i<NX_ * NY_ * NZ_; i++){
                    sum_[i] = cplx_t(0.0, 0.0);
                }
            }
       
            void assign(const F64 charge, const F64vec & pos){
                assert(is_initialized_);
                Waves w;
                w.initialize(MMAX_, period_length_, pos);
                for(S32 iz=0; iz<NZ_; iz++){
                    for(S32 iy=0;  iy<NY_; iy++){
                        for(S32 ix=0;  ix<NX_; ix++){
                            const S32 idx = get_index(ix,iy,iz);
                            sum_[idx] += ((charge * w.wz_[iz]) * w.wy_[iy]) * w.wx_[ix];
                        }
                    }
                }
            }

            void reduce() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                const S32 count = size();
                real_t * sendbuf = new real_t[2 * count];
                real_t * recvbuf = new real_t[2 * count];
                // Reduce the real part
                for (S32 i=0; i<count; i++) {
                    sendbuf[2 * i]     = sum_[i].real();
                    sendbuf[2 * i + 1] = sum_[i].imag();
                }
                MPI_Allreduce(sendbuf, recvbuf, 2*count,
                              GetDataType<real_t>(), MPI_SUM, MPI_COMM_WORLD);
                for (S32 i=0; i<count; i++) {
#ifdef PARTICLE_SIMULATOR_STD_COMPLEX_NOT_HAVING_SETTER 
                    const cplx_t ctmp (recvbuf[2 * i], recvbuf[2 * i + 1]);
                    sum_[i] = ctmp;
                    // See Note #64d4cd48 in particle_mesh_multipole/M2L_engine.hpp.
#else
                    sum_[i].real(recvbuf[2 * i]);
                    sum_[i].imag(recvbuf[2 * i + 1]);
#endif
                }
                // Free
                delete [] sendbuf;
                delete [] recvbuf;
#endif
            }
        
            void filter(const F64 alpha){
                assert(is_initialized_);
                const F64 pi = 4.0 * std::atan(1.0);
                const F64 Lx = period_length_.x;
                const F64 Ly = period_length_.y;
                const F64 Lz = period_length_.z;
                const F64 coef = 1.0/(4.0*alpha*alpha);
                for(S32 iz=0; iz<NZ_; iz++){
                    for(S32 iy=0; iy<NY_; iy++){
                        for(S32 ix=0; ix<NX_; ix++){
                            const S32 mx = ix - MMAX_;
                            const S32 my = iy - MMAX_;
                            const S32 mz = iz - MMAX_;
                            const F64 g2 = (4.0*pi*pi) * (mx*mx/(Lx*Lx) 
                                                         +my*my/(Ly*Ly)
                                                         +mz*mz/(Lz*Lz));
                            if (g2 == 0.0) continue;
                            const F64 factor = std::exp(-coef * g2) / g2;
                            const S32 idx = get_index(ix,iy,iz);
                            sum_[idx] *= factor;
                        }
                    }
                }
                const S32 idx = get_index(MMAX_, MMAX_, MMAX_);
                sum_[idx] = 0.0;
            }
        
            void phi_and_grad(const F64vec pos,
                              F64 & phi,
                              F64vec & acc) const {
                assert(is_initialized_);
                const F64 pi = 4.0 * std::atan(1.0);
                const F64 Lx = period_length_.x;
                const F64 Ly = period_length_.y;
                const F64 Lz = period_length_.z;
                const F64 coef = 4.0*pi/(Lx*Ly*Lz);
                Waves w;
                w.initialize(MMAX_, period_length_, pos);
                F64 phi_tmp {0.0};
                F64vec acc_tmp(0.0);
                for (S32 iz=0; iz<NZ_; iz++){
                    for (S32 iy=0; iy<NY_; iy++){
                        for (S32 ix=0; ix<NX_; ix++){
                            const S32 mx = ix - MMAX_;
                            const S32 my = iy - MMAX_;
                            const S32 mz = iz - MMAX_;
                            const F64vec gvec = (2.0*pi) * F64vec(mx/Lx, my/Ly, mz/Lz);
                            const cplx_t cs = (w.wz_[iz] * w.wy_[iy]) * w.wx_[ix];
                            const F64 c = cs.real();
                            const F64 s = cs.imag();
                            const F64 csum = sum_[get_index(ix,iy,iz)].real();
                            const F64 ssum = sum_[get_index(ix,iy,iz)].imag();
                            
                            phi_tmp += c*csum + s*ssum;
                            acc_tmp += (s*csum - c*ssum) * gvec;
                        }
                    }
                }
                phi += coef * phi_tmp;
                acc += - coef * acc_tmp;
            }
        };
        
        template <class Tfp>
        void eval_k_space(const F64ort & pos_unit_cell,
                          const F64 alpha,
                          const S32 MMAX,
                          ParticleSystem<Tfp> & psys)
        {
            const F64vec period_length = pos_unit_cell.getFullLength();
            EwaldSum sum;
            sum.initialize(MMAX, period_length); 
            sum.clear();
            for (S32 i=0; i<psys.getNumberOfParticleLocal(); i++){
                const F64 charge = psys[i].getCharge();
                const F64vec pos = psys[i].getPos();
                sum.assign(charge, pos);
            }
            sum.reduce();
            sum.filter(alpha);
            for (S32 i=0; i<psys.getNumberOfParticleLocal(); i++){
                const F64vec pos = psys[i].getPos();
                F64 phi {0.0};
                F64vec acc(0.0);
                sum.phi_and_grad(pos, phi, acc);
                psys[i].accumulateForceEwald(acc, phi);
            }
        }
        
        inline void cutfunc(const F64 r,
                            const F64 alpha,
                            const F64 rsq,
                            const F64 asq,
                            F64 & pcut,
                            F64 & fcut) {
            const F64 pi = 4.0 * std::atan(1.0);
            const F64 c = 2.0 / std::sqrt(pi);
#ifdef PARTICLE_SIMULATOR_USE_ERFC_IN_MATH_H
            const F64 tmp = erfc(alpha * r);
#else
            const F64 tmp = std::erfc(alpha * r);
#endif
            pcut = tmp;
            fcut = tmp + c * (alpha * r) * std::exp(-asq * rsq);
        }

        inline F64 minimum_image(const F64 dr,
                                 const F64 period_length) { // for 1D
            if (dr >  0.5 * period_length) return (dr - period_length);
            if (dr < -0.5 * period_length) return (dr + period_length);
            return dr;
        }

        inline F64vec minimum_image(const F64vec dr,
                                    const F64vec period_length) { // for 3D
            return F64vec(minimum_image(dr.x, period_length.x),
                          minimum_image(dr.y, period_length.y),
                          minimum_image(dr.z, period_length.z));
        }
        
        template <class Tfp>
        void eval_r_space(const F64ort & pos_unit_cell,
                          const F64 alpha,
                          const S32 NMIR,
                          const F64 msum,
                          ParticleSystem<Tfp> & psys)
        {
            // Local class
            class Epj {
            public:
                F64 charge;
                F64vec pos;
                F64 getCharge() const {return charge;}
                F64vec getPos() const {return pos;}
            };
            const S32 n_proc = Comm::getNumberOfProc();
            const S32 my_rank = Comm::getRank();
            const F64 pi = 4.0 * std::atan(1.0);
            const F64 ainv2 = pi / (alpha * alpha);
            const F64vec period_length = pos_unit_cell.getFullLength();
            // Collect number of particles of other processes.
            const S32 n_loc = psys.getNumberOfParticleLocal();
            ReallocatableArray<S32> n_ptcls;
            n_ptcls.reserveAtLeast(n_proc);
            for (S32 i=0; i<n_proc; i++) n_ptcls[i] = 0;
            Comm::allGather(&n_loc, 1, n_ptcls.getPointer());
            // Calculate force in the real-space
            ReallocatableArray<Epj> p_send, p_recv;
            p_send.reserveEmptyAreaAtLeast(n_loc);
            for (S32 i=0; i<n_loc; i++) {
                p_send[i].charge = psys[i].getCharge();
                p_send[i].pos    = psys[i].getPos();
            }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            for (S32 jj = 0; jj < n_proc; jj++) {
                // Get particle information of process j
                const S32 rank_recv = (my_rank + jj) % n_proc;
                const S32 rank_send = (n_proc + my_rank - jj) % n_proc;
                const S32 n_jp = n_ptcls[rank_recv];
                p_recv.reserveEmptyAreaAtLeast(n_jp);
                MPI_Request req_send, req_recv;
                MPI_Status stat_send, stat_recv;
                MPI_Irecv(p_recv.getPointer(), n_jp, GetDataType<Epj>(),
                          rank_recv, 0, MPI_COMM_WORLD, &req_recv);
                MPI_Isend(p_send.getPointer(), n_loc, GetDataType<Epj>(),
                          rank_send, 0, MPI_COMM_WORLD, &req_send);
                MPI_Wait(&req_send, &stat_send);
                MPI_Wait(&req_recv, &stat_recv);
#else
            const S32 jj = n_proc - 1;
            const S32 n_jp = n_loc;
            p_recv.reserveEmptyAreaAtLeast(n_jp);
            for (S32 i=0; i < n_loc; i++) {
                p_recv[i].charge = p_send[i].charge;
                p_recv[i].pos    = p_send[i].pos;
            }
#endif
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                for(S32 i=0; i<n_loc; i++){
                    F64 phi {0.0};
                    F64vec grad(0.0);
                    for(S32 mz=-NMIR; mz<=NMIR; mz++){
                        for(S32 my=-NMIR; my<=NMIR; my++){
                            for(S32 mx=-NMIR; mx<=NMIR; mx++){
                                const F64vec shift = F64vec(period_length.x * mx,
                                                            period_length.y * my,
                                                            period_length.z * mz);
                                F64 ptmp {0.0};
                                F64vec gtmp(0.0);
                                for(S32 j=0; j<n_jp; j++){
                                    F64vec dr = p_recv[j].getPos() - p_send[i].getPos();
                                    dr = minimum_image(dr, period_length) + shift;
                                    
                                    const F64 r2 = dr * dr;
                                    if (0.0 == r2) continue;
                                    const F64 r  = sqrt(r2);
                                    const F64 rinv = 1.0 / r;
                                    const F64 qri  = p_recv[j].getCharge() * rinv;
                                    const F64 qri3 = qri * (rinv * rinv);
                                    
                                    F64 pcut, fcut;
                                    cutfunc(r, alpha, r2, alpha*alpha, pcut, fcut);
                                    
                                    ptmp += pcut * qri;
                                    gtmp += (fcut * qri3) * dr;
                                } // for(j)
                                phi  += ptmp;
                                grad += gtmp;
                            }
                        }
                    }
                    // self energy
                    if (jj == (n_proc - 1)) {
                        phi -= alpha * (2.0/std::sqrt(pi)) * p_send[i].getCharge();
                        phi -= msum * ainv2; // charged system correction
                    }
                    psys[i].accumulateForceEwald(grad, phi);
                } // for(i)
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            } // for(j_rank)
#endif
            // Free
            n_ptcls.freeMem(1);
            p_send.freeMem(1);
            p_recv.freeMem(1);
        }

        class Ewald {
        public:
            bool is_initialized_;
            F64ort pos_unit_cell_;
            F64 alpha_;
            S32 NMIR_;
            S32 MMAX_;

            Ewald(): is_initialized_(false) {}

            void initialize(const F64ort & pos_unit_cell,
                            const F64 alpha = 2.4,
                            const S32 NMIR = 3,
                            const S32 MMAX = 5) {
                pos_unit_cell_ = pos_unit_cell;
                alpha_ = alpha;
                NMIR_ = NMIR;
                MMAX_ = MMAX;
                is_initialized_ = true;
                // Output parameters
                if (Comm::getRank() == 0) {
                    std::cout << "-------------------------------" << std::endl;
                    std::cout << "Parameters in the Ewald method:" << std::endl;
                    std::cout << "    alpha_ = " << alpha_         << std::endl;
                    std::cout << "    NMIR_  = " << NMIR_          << std::endl;
                    std::cout << "    MMAX_  = " << MMAX_          << std::endl;
                    std::cout << "-------------------------------" << std::endl;
                }
            }
        
            template <class Tfp>    
            void calcForceAllAndWriteBack(ParticleSystem<Tfp> & psys,
                                          const bool clear_force = true) {
                assert(is_initialized_);
                // Clear force
                if (clear_force) 
                    for (S32 i = 0; i < psys.getNumberOfParticleLocal(); i++)
                        psys[i].clearEwald();
                eval_k_space(pos_unit_cell_, alpha_, MMAX_, psys);
                F64 msum_loc = 0.0;
                for (S32 i = 0; i < psys.getNumberOfParticleLocal(); i++)
                    msum_loc += psys[i].getCharge();
                F64 msum = Comm::getSum(msum_loc);
                eval_r_space(pos_unit_cell_, alpha_, NMIR_, msum, psys);
            }
        };

    }; // END of namespace of Ewald
}; // END of namespace of ParticleSimulator
