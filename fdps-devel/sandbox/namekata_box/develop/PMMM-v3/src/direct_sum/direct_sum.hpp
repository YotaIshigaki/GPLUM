#pragma once
#include <cmath>
#include "ps_defs.hpp"

namespace ParticleSimulator {
    namespace DirectSum {

        class DirectSum {
        public:
            bool is_initialized_;

            DirectSum(): is_initialized_(false) {}

            void initialize() {
                is_initialized_ = true;
            }

            template <class Tfp>
            void calcForceAllAndWriteBack(ParticleSystem<Tfp> & psys,
                                          const bool clear_force = true) {
                assert(is_initialized_);
                class Epj {
                public:
                    F64 charge;
                    F64vec pos;
                    F64 getCharge() const {return charge;}
                    F64vec getPos() const {return pos;}
                };
                const S32 n_proc = Comm::getNumberOfProc();
                const S32 my_rank = Comm::getRank();
                // Collect number of particles of other processes
                const S32 n_loc = psys.getNumberOfParticleLocal();
                ReallocatableArray<S32> n_ptcls;
                n_ptcls.reserveAtLeast(n_proc);
                for (S32 i=0; i<n_proc; i++) n_ptcls[i] = 0;
                Comm::allGather(&n_loc, 1, n_ptcls.getPointer());
                // Clear force
                if (clear_force)
                    for (S32 i = 0; i < n_loc; i++)
                        psys[i].clearDirect();
                // Calculate force
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
                        for(S32 j=0; j<n_jp; j++){
                            const F64 qj = p_recv[j].getCharge();
                            const F64vec dr = p_recv[j].getPos() - p_send[i].getPos();
                            const F64 r2 = dr * dr;
                            if (0.0 == r2) continue;
                            const F64 r  = sqrt(r2);
                            const F64 rinv = 1.0 / r;
                            const F64 qri  = qj * rinv;
                            const F64 qri3 = qri * (rinv * rinv);
                            phi += qri;
                            grad += qri3 * dr;
                            // Note that grad is \nabla\phi instead of -\nabla\phi.
                        } // for(j)
                        psys[i].accumulateForceDirect(grad, phi);
                    } // for(i)
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                } // for(j_rank)
#endif
                // Free
                n_ptcls.freeMem(1);
                p_send.freeMem(1);
                p_recv.freeMem(1);
            }

        };

    }; // END of namespace of DirectSum
}; // END of namespace of ParticleSimulator
