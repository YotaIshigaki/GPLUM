#pragma once
#include "../ps_defs.hpp"

namespace ParticleSimulator {
    namespace ParticleMeshMultipole {

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
                         CommBuffer<recv_t> & recv) {
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
                              rnk, tag, MPI_COMM_WORLD, &req[n_req_used++]);
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
                             rnk, tag, MPI_COMM_WORLD);
                }
            }
            if (n_req_used > 0) {
                MPI_Status *stat = new MPI_Status[n_req_used];
                MPI_Waitall(n_req_used, req, stat);
                delete [] stat;
                delete [] req;
            }
        }


#endif // PARTICLE_SIMULATOR_MPI_PARALLEL
    } // END of namespace of ParticleSimulator
} // END of namespace of ParticleSimulator

