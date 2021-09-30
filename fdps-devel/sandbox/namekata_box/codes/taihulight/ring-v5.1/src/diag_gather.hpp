#pragma once

#include "ps_defs.hpp"

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include<mpi.h>
#endif

namespace ParticleSimulator{

    template<class Tspj>
    struct DiagGather {
        MPI_Comm comm_rad;
        MPI_Comm comm_phi;
        MPI_Comm comm_short;
        MPI_Comm comm_long;

        int rank_glb, size_glb;
        int rank_rad, size_rad;
        int rank_phi, size_phi;
        int rank_short, size_short;
        int rank_long,  size_long;

        int icount;
        int *counts;
        int *displs;

        DiagGather (int nrad, int nphi, const MPI_Comm &comm_glb) {
            int color, key;

            MPI_Comm_rank(comm_glb, &rank_glb);
            MPI_Comm_size(comm_glb, &size_glb);
            assert(size_glb == nrad*nphi);

            rank_rad = rank_glb % nrad;
            rank_phi = rank_glb / nrad;

            MPI_Comm_split(comm_glb, color=rank_phi, key=rank_rad, &comm_rad);
            MPI_Comm_split(comm_glb, color=rank_rad, key=rank_phi, &comm_phi);

            MPI_Comm_size(comm_rad, &size_rad);
            MPI_Comm_size(comm_phi, &size_phi);

            assert(nrad == size_rad);
            assert(nphi == size_phi);

            color = rank_phi / nrad;
            key   = rank_phi % nrad;
            if (nphi % nrad) {
                int ncolor = nphi / nrad;
                if (rank_phi >= nrad * ncolor) {
                    color--;
                    key += nrad;
                }
            }
            MPI_Comm_split(comm_phi, color, key, &comm_short);
            MPI_Comm_size(comm_short, &size_short);
            MPI_Comm_rank(comm_short, &rank_short);

            color = rank_phi % nrad;
            key   = rank_phi / nrad;

            MPI_Comm_split(comm_phi, color, key, &comm_long);
            MPI_Comm_size(comm_long, &size_long);
            MPI_Comm_rank(comm_long, &rank_long);

            counts = new int[nrad];
            displs = new int[nrad+1];

            icount = size_long;
            int root = rank_rad;
            MPI_Bcast(&icount, 1, MPI_INT, root, comm_short);
            MPI_Allgather(
                  &icount, 1, MPI_INT,
                  counts,  1, MPI_INT,
                  comm_rad);
            displs[0] = 0;
            for(int i=0; i<nrad; i++){
                displs[i+1] = displs[i] + counts[i];
            }
        }
    
        ~DiagGather() {
            delete [] counts;
            delete [] displs;
            MPI_Comm_free(&comm_rad);
            MPI_Comm_free(&comm_phi);
            MPI_Comm_free(&comm_short);
            MPI_Comm_free(&comm_long);
        }
    
        void show_comm_info(){
            printf("glb:(%2d, %2d), rad:(%2d, %2d), phi:(%2d, %2d), short:(%2d, %2d), long:(%2d, %2d)\n",
                  rank_glb, size_glb,
                  rank_rad, size_rad,
                  rank_phi, size_phi,
                  rank_short, size_short,
                  rank_long, size_long);
        }

        void gather_cm_fast(const Tspj &spj, Tspj spj_super_domain[]) {
            // Memory allocation of local buffers
            ReallocatableArray<Tspj> spj_sub_domain;
            ReallocatableArray<Tspj> work;
            spj_sub_domain.resizeNoInitialize(size_rad);
            work.resizeNoInitialize(size_phi);


            // Compute my_spj_super_domain
            MPI_Allgather(&spj, 1, GetDataType<Tspj>(),
                          spj_sub_domain.getPointer(), 1, GetDataType<Tspj>(),
                          comm_rad);


            int root = rank_phi % size_rad;
            if (rank_rad == root) {
                Tspj my_spj_super_domain;
                my_spj_super_domain.clear();
                for (S32 i=0; i<size_rad; i++){
                    my_spj_super_domain.mass += spj_sub_domain[i].mass;
                    my_spj_super_domain.pos += spj_sub_domain[i].mass * spj_sub_domain[i].pos;
#ifdef PHI_R_TREE
                    my_spj_super_domain.pos_phi += spj_sub_domain[i].mass * spj_sub_domain[i].pos_phi;
                    my_spj_super_domain.pos_r += spj_sub_domain[i].mass * spj_sub_domain[i].pos_r;
#endif
                }
                my_spj_super_domain.pos     /= my_spj_super_domain.mass;
#ifdef PHI_R_TREE
                my_spj_super_domain.pos_phi /= my_spj_super_domain.mass;
                my_spj_super_domain.pos_r   /= my_spj_super_domain.mass;
#endif
                MPI_Allgather(
                      &my_spj_super_domain, 1, GetDataType<Tspj>(),
                      work.getPointer(), 1, GetDataType<Tspj>(),
                      comm_long);
            }
    
            root = rank_rad;
            MPI_Bcast(work.getPointer(), icount, GetDataType<Tspj>(), root, comm_short);
    
            MPI_Allgatherv(
                  work.getPointer(), icount, GetDataType<Tspj>(),
                  spj_super_domain, counts, displs, GetDataType<Tspj>(),
                  comm_rad);

            // Reorder 
            for (S32 i=0; i<size_phi; i++) work[i] = spj_super_domain[i];
            for (S32 i=0; i<size_rad; i++) {
                const S32 nj   = counts[i];
                const S32 joff = displs[i];
                for (S32 j=0; j<nj; j++){
                    spj_super_domain[i + size_rad*j] = work[joff + j];
                }
            }

        }

    };

}
