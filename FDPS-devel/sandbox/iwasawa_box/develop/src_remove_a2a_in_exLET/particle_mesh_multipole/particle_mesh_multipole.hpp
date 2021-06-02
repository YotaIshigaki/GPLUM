#pragma once
// Include header file(s) of FDPS core part
#include "../ps_defs.hpp"
// Include header files of FFTW library
#include <fftw3.h>
// Define macros (for test)
//#define EWALD_FILE
//#define NON_CHARGE_NEUTRAL
// Include header files for PMMM
#include "multidimensional_array.hpp"
#include "fmm.hpp"
#include "particle.hpp"
#include "cell.hpp"
#include "cutoff.hpp"
#include "ewald.hpp"
#include "particle_mesh_utils.hpp"
#include "particle_particle_utils.hpp"

static void print_err(
        const std::vector<double> &err,
        const char * name,
        const int icut,
        const int p)
{
    static char fname[256];
    sprintf(fname, "%s.c%dp%d.dat", name, icut, p);
    FILE *fp = fopen(fname, "w");
    assert(fp);
    const int len = err.size();
    for(int i=0; i<len; i++){
        fprintf(fp, "%e %e\n", double(i)/len, err[i]);
    }
    fclose(fp);
}

namespace ParticleSimulator {

    template <int PFMM, int ICUT>
    class ParticleMeshMultipole {
    // < Definitions of template arguments >
    //     PFMM := the maximum order of multipole expansion
    //     ICUT := the minimu cell separation
    public: 
        typedef Cell_FMM<PFMM> Cell_t;

    private:
        bool first_call_by_initialize;

        // Domain information
        F64ort * pos_domain_;
        F64ort pos_root_domain_;
        S32 boundary_condition_;
        bool periodic_axis_[DIMENSION_LIMIT];
 
        // Particle information
        Particle *ptcl_;
        Particle *ptcl_recv_;
        S32 n_loc_tot_;  // the number of local particles
        S32 n_recv_tot_; // the number of particles which are sent to this process
        S32 n_glb_tot_; // n_loc_tot_ + (# of ptcls in interacting nearby cells)
        F64 msum_; // total charge 

        // Cell information
        MultidimensionalArray<Cell_t,3> cell_;
        F64ort pos_root_cell_;
        S32 NC_, NC3_;
        F64 clen_;

    public:

        ParticleMeshMultipole() {
            first_call_by_initialize = true;            
            // Initialize domain information
            pos_domain_ = NULL;
            periodic_axis_[0] = periodic_axis_[1] = false;
            pos_root_domain_.low_.x  = - LARGE_FLOAT;
            pos_root_domain_.high_.x =   LARGE_FLOAT;
            pos_root_domain_.low_.y  = - LARGE_FLOAT;
            pos_root_domain_.high_.y =   LARGE_FLOAT;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            periodic_axis_[2] = false;
            pos_root_domain_.low_.z  = - LARGE_FLOAT;
            pos_root_domain_.high_.z =   LARGE_FLOAT;
#endif
            boundary_condition_ = BOUNDARY_CONDITION_OPEN;
            // Initialize particle information 
            ptcl_ = NULL;
            // Initialize cell information
            pos_root_cell_.low_.x  = - LARGE_FLOAT;
            pos_root_cell_.high_.x =   LARGE_FLOAT;
            pos_root_cell_.low_.y  = - LARGE_FLOAT;
            pos_root_cell_.high_.y =   LARGE_FLOAT;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            pos_root_cell_.low_.z  = - LARGE_FLOAT;
            pos_root_cell_.high_.z =   LARGE_FLOAT;
#endif
        }

        ~ParticleMeshMultipole() {
            // Finalize domain information
            if (pos_domain_ != NULL) delete [] pos_domain_;
            // Finalize particle information
            if (ptcl_ != NULL) delete [] ptcl_;
            // Finalize cell information
        }

        void initialize() {
            assert(first_call_by_initialize);
            first_call_by_initialize = false;
            // Initialize domain information
            pos_domain_ = new F64ort[Comm::getNumberOfProc()];
            // Initialize particle information
            // Initialize cell information 
        }

        template <class Tdinfo>
        void setDomainInfoPMMM(const Tdinfo & dinfo) {
            for (S32 i = 0; i < Comm::getNumberOfProc(); i++)
                pos_domain_[i] = dinfo.getPosDomain(i);
            pos_root_domain_    = dinfo.getPosRootDomain();
            boundary_condition_ = dinfo.getBoundaryCondition();
            dinfo.getPeriodicAxis(periodic_axis_);
            // Currently, there are the following limitations:
            // (1) Only open boundary condition or 
            //     periodic boundary condition for all directions 
            //     are supported.
            if ((boundary_condition_ != BOUNDARY_CONDITION_OPEN) &&
                (boundary_condition_ != BOUNDARY_CONDITION_PERIODIC_XYZ)) {
                std::cerr << "Boundary condition " << boundary_condition_ 
                          << "is not supported." << std::endl;
                Abort(-1);
            }
            // (2) Basic root cell must be [0,1)^{3} when the boundary
            //     condition is periodic for all directions.
            if (boundary_condition_ == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                if ((pos_root_domain_.low_.x != 0.0) ||
                    (pos_root_domain_.low_.y != 0.0) ||
                    (pos_root_domain_.low_.z != 0.0) ||
                    (pos_root_domain_.high_.x != 1.0) ||
                    (pos_root_domain_.high_.y != 1.0) ||
                    (pos_root_domain_.high_.z != 1.0)) {
                    std::cerr << "Domain other than [0,1)^3 is not supported." << std::endl;
                    std::cerr << pos_root_domain_ << std::endl;
                    Abort(-1);
                }
            } else {
                // When the open boundary, we overwrite pos_root_domain_
                // because it is initialized by LARGE_FLOAT even if an user
                // of FDPS call dinfo.setPosRootDomain() in his/her code.
                pos_root_domain_.low_  = F64vec(0.0, 0.0, 0.0);
                pos_root_domain_.high_ = F64vec(1.0, 1.0, 1.0);
                // [TODO] This limitation MUST be removed in future.
                //    If we adopt a cube which contains all the particles as a root cell,
                //    we can remove this limination.
            }
            // Correct pos_domain_ so that pos_domain_[i] is contained in pos_root_domain_
            for (S32 i = 0; i < Comm::getNumberOfProc(); i++) {
                if (pos_domain_[i].low_.x < pos_root_domain_.low_.x)
                    pos_domain_[i].low_.x = pos_root_domain_.low_.x;
                if (pos_domain_[i].low_.y < pos_root_domain_.low_.y)
                    pos_domain_[i].low_.y = pos_root_domain_.low_.y;
                if (pos_domain_[i].low_.z < pos_root_domain_.low_.z)
                    pos_domain_[i].low_.z = pos_root_domain_.low_.z;

                if (pos_domain_[i].high_.x > pos_root_domain_.high_.x)
                    pos_domain_[i].high_.x = pos_root_domain_.high_.x;
                if (pos_domain_[i].high_.y > pos_root_domain_.high_.y)
                    pos_domain_[i].high_.y = pos_root_domain_.high_.y;
                if (pos_domain_[i].high_.z > pos_root_domain_.high_.z)
                    pos_domain_[i].high_.z = pos_root_domain_.high_.z;
            }
        }

        template <class Tpsys>
        void setParticlePMMM(const Tpsys & psys,
                             const bool clear=true) {
            if (clear) {
                if (ptcl_ != NULL) delete [] ptcl_;
                n_loc_tot_ = psys.getNumberOfParticleLocal();
                ptcl_ = new Particle[n_loc_tot_];
                for (S32 i = 0; i < n_loc_tot_; i++) {
                    ptcl_[i].mass  = psys[i].getChargePMMM();
                    ptcl_[i].pos.x = psys[i].getPos().x;
                    ptcl_[i].pos.y = psys[i].getPos().y;
                    ptcl_[i].pos.z = psys[i].getPos().z;
                }
            } else {
                Particle *ptcl_prev = new Particle[n_loc_tot_];
                for (S32 i = 0; i < n_loc_tot_; i++)
                    ptcl_prev[i] = ptcl_[i];
                delete [] ptcl_;

                S32 n_app = psys.getNumberOfParticleLocal();
                ptcl_ = new Particle[n_loc_tot_ + n_app];
                for (S32 i = 0; i < n_loc_tot_; i++)
                    ptcl_[i] = ptcl_prev[i];
                for (S32 i = 0; i < n_app; i++) {
                    const S32 ii = i + n_loc_tot_;
                    ptcl_[ii].mass  = psys[i].getChargePMMM();
                    ptcl_[ii].pos.x = psys[i].getPos().x;
                    ptcl_[ii].pos.y = psys[i].getPos().y;
                    ptcl_[ii].pos.z = psys[i].getPos().z;
                }
                n_loc_tot_ += n_app;
 
                delete [] ptcl_prev;
            }
        }

        void setCell() {
            // Set pos_root_cell_
            if (boundary_condition_ == BOUNDARY_CONDITION_OPEN) {
#if 0
                for (S32 i = 0; i < n_loc_tot_; i++) {
                    F64vec pos;
                    pos.x = ptcl_[i].pos.x;
                    pos.y = ptcl_[i].pos.y;
                    pos.z = ptcl_[i].pos.z;
                    if (pos.x < pos_root_cell_.low_.x) pos_root_cell_.low_.x = pos.x;
                    if (pos.y < pos_root_cell_.low_.y) pos_root_cell_.low_.y = pos.y;
                    if (pos.z < pos_root_cell_.low_.z) pos_root_cell_.low_.z = pos.z;
                    if (pos.x > pos_root_cell_.high_.x) pos_root_cell_.high_.x = pos.x;
                    if (pos.y > pos_root_cell_.high_.y) pos_root_cell_.high_.y = pos.y;
                    if (pos.z > pos_root_cell_.high_.z) pos_root_cell_.high_.z = pos.z;
                }
#else
                pos_root_cell_ = pos_root_domain_;
#endif
            } else if (boundary_condition_ == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                pos_root_cell_ = pos_root_domain_;
            } else {
                if (Comm::getRank() == 0) {
                    std::cerr << "Boundary condition " << boundary_condition_ 
                              << " is not supported." << std::endl;
                }
                Abort(-1);
            }
            // Calculate # of cells
#if 0
            const S64 n_ptcl_ = Comm::getSum(n_loc_tot_);
            const S32 Kest = n_ptcl_ / (PFMM*PFMM);
            NC_ = 2;
            for (;;) {
                if (NC_*NC_*NC_ > Kest) break;
                NC_ += 2;
            }
#else
            NC_ = 8; // for test
#endif
            // Set NC3_
            NC3_ = NC_ * NC_ * NC_;
            // Set cell_
            cell_.initialize(NC_, NC_, NC_);
            clen_ = 1.0 / NC_;
            for(S32 k=0; k<NC_; k++){
                for(S32 j=0; j<NC_; j++){
                    for(S32 i=0; i<NC_; i++){
                        cell_(k,j,i).set(S32vec(i,j,k), clen_);
                    }
                }
            }
        }

        void assignParticleToCellLocal() {

            msum_ = 0.0;
            for (S32 i=0; i<n_loc_tot_; i++){
                msum_ += ptcl_[i].mass;
                const S32vec idx = cell_nearest(ptcl_[i].pos, clen_);
                assert(0 <= idx.x && idx.x < NC_);
                assert(0 <= idx.y && idx.y < NC_);
                assert(0 <= idx.z && idx.z < NC_);
                cell_(idx.z,idx.y,idx.x).plist.push_back(&ptcl_[i]);
            }
        
            for(S32 k=0; k<NC_; k++)
                for(S32 j=0; j<NC_; j++)
                    for(S32 i=0; i<NC_; i++)
                        cell_(k,j,i).sanity_check();

        }

        void exchangeParticle() {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            Comm::barrier();
            if (Comm::getRank() == 0) std::cout << "PP part: exchange particle (start)" << std::endl;

            // Local variables and initialize them
            const S32 myrank = Comm::getRank(); 
            const S32 nproc = Comm::getNumberOfProc();
            S32 nsend[nproc];
            S32 nsend_disp[nproc];
            S32 nrecv[nproc];
            S32 nrecv_disp[nproc]; 
            for (S32 i = 0; i < nproc; i++) {
                nsend[i] = nsend_disp[i] = nrecv[i] = nrecv_disp[i] = 0;
            }
        
            // Compute idx_domain[], idx_req[]
            S32vec idx_domain_low[nproc];
            S32vec idx_domain_high[nproc];
            S32vec idx_req_low[nproc];
            S32vec idx_req_high[nproc];
            for (S32 i = 0; i < nproc; i++) {
                F64vec pos;
                // idx_domain[]
                // (1) low_
                pos = pos_domain_[i].low_;
                idx_domain_low[i] = cell_nearest(pos, clen_);
                if (idx_domain_low[i].x < 0) idx_domain_low[i].x = 0;
                if (idx_domain_low[i].y < 0) idx_domain_low[i].y = 0;
                if (idx_domain_low[i].z < 0) idx_domain_low[i].z = 0;
                if (idx_domain_low[i].x >= NC_) idx_domain_low[i].x = NC_-1;
                if (idx_domain_low[i].y >= NC_) idx_domain_low[i].y = NC_-1;
                if (idx_domain_low[i].z >= NC_) idx_domain_low[i].z = NC_-1;
                // (2) high_
                pos = pos_domain_[i].high_;
                idx_domain_high[i] = cell_nearest(pos, clen_);
                if (idx_domain_high[i].x < 0) idx_domain_high[i].x = 0;
                if (idx_domain_high[i].y < 0) idx_domain_high[i].y = 0;
                if (idx_domain_high[i].z < 0) idx_domain_high[i].z = 0;
                if (idx_domain_high[i].x >= NC_) idx_domain_high[i].x = NC_-1;
                if (idx_domain_high[i].y >= NC_) idx_domain_high[i].y = NC_-1;
                if (idx_domain_high[i].z >= NC_) idx_domain_high[i].z = NC_-1;
                // idx_req[]
                // (1) low_
                idx_req_low[i]  = idx_domain_low[i];
                idx_req_low[i].x -= ICUT;
                idx_req_low[i].y -= ICUT;
                idx_req_low[i].z -= ICUT;
                if (idx_req_low[i].x < 0) idx_req_low[i].x = 0;
                if (idx_req_low[i].y < 0) idx_req_low[i].y = 0;
                if (idx_req_low[i].z < 0) idx_req_low[i].z = 0;
                if (idx_req_low[i].x >= NC_) idx_req_low[i].x = NC_-1;
                if (idx_req_low[i].y >= NC_) idx_req_low[i].y = NC_-1;
                if (idx_req_low[i].z >= NC_) idx_req_low[i].z = NC_-1;
                // (2) high_
                idx_req_high[i] = idx_domain_high[i];
                idx_req_high[i].x += ICUT;
                idx_req_high[i].y += ICUT;
                idx_req_high[i].z += ICUT;
                if (idx_req_high[i].x < 0) idx_req_high[i].x = 0;
                if (idx_req_high[i].y < 0) idx_req_high[i].y = 0;
                if (idx_req_high[i].z < 0) idx_req_high[i].z = 0;
                if (idx_req_high[i].x >= NC_) idx_req_high[i].x = NC_-1;
                if (idx_req_high[i].y >= NC_) idx_req_high[i].y = NC_-1;
                if (idx_req_high[i].z >= NC_) idx_req_high[i].z = NC_-1;
            }
            // check
#if 0
            std::cout << "pos_domain: " << pos_domain_[myrank]
                      << " (rank = " << myrank << ")" << std::endl;
            printf("idx_domain: (%d, %d, %d) (%d, %d, %d) (rank = %d)\n",
                   idx_domain_low[myrank].x, idx_domain_low[myrank].y, idx_domain_low[myrank].z,
                   idx_domain_high[myrank].x, idx_domain_high[myrank].y, idx_domain_high[myrank].z,
                   myrank);
            printf("idx_req: (%d, %d, %d) (%d, %d, %d) (rank = %d)\n",
                   idx_req_low[myrank].x, idx_req_low[myrank].y, idx_req_low[myrank].z,
                   idx_req_high[myrank].x, idx_req_high[myrank].y, idx_req_high[myrank].z,
                   myrank);
            //Finalize();
            //std::exit(0);
#endif

            // Compute nsend[], nsend_disp[]
            const S32vec idx_my_domain_low  = idx_domain_low[myrank];
            const S32vec idx_my_domain_high = idx_domain_high[myrank];
            for (S32 irank = 0; irank < nproc; irank++) {
                nsend[irank] = 0;
                if (irank == myrank) continue;
                // Calculate the number of particles which are sent
                const S32 kstart = idx_my_domain_low.z;
                const S32 kend   = idx_my_domain_high.z;
                const S32 jstart = idx_my_domain_low.y;
                const S32 jend   = idx_my_domain_high.y;
                const S32 istart = idx_my_domain_low.x;
                const S32 iend   = idx_my_domain_high.x;
                for (S32 k = kstart; k <= kend; k++) {
                    if ((k < idx_req_low[irank].z) || 
                        (k > idx_req_high[irank].z)) continue;
                    for (S32 j = jstart; j <= jend; j++) {
                        if ((j < idx_req_low[irank].y) ||
                            (j > idx_req_high[irank].z)) continue;
                        for (S32 i = istart; i <= iend; i++) {
                            if ((i < idx_req_low[irank].x) ||
                                (i > idx_req_high[irank].x)) continue;
                            nsend[irank] += cell_(k,j,i).plist.size();
                        }
                    }
                }
#if 1
                // check
                std::cout << "nsend (" << myrank << " -> " << irank << "): "
                          << nsend[irank] << std::endl;
#endif
            }
            nsend_disp[0] = 0;
            if (nproc > 1) 
                for (S32 i = 1; i < nproc; i++)
                    nsend_disp[i] = nsend_disp[i-1] + nsend[i-1];


            // Compute nrecv[], nrecv_disp[]
            Comm::allToAll(nsend, 1, nrecv);
            nrecv_disp[0] = 0;
            if (nproc > 1)
                for (S32 i = 1; i < nproc; i++)
                    nrecv_disp[i] = nrecv_disp[i-1] + nrecv[i-1];
            // check 
            for (S32 irank = 0; irank < nproc; irank++) {
                if (irank == myrank) continue;
                std::cout << "nrecv (" << myrank << " <- " << irank << "): "
                          << nrecv[irank] << std::endl;
            }

           
            // Make a buffer for send
            S32 nsend_tot = 0;
            for (S32 i = 0; i < nproc; i++) nsend_tot += nsend[i];
            Particle *ptcl_send = new Particle[nsend_tot];
            S32 disp = 0;
            for (S32 irank = 0; irank < nproc; irank++) {
                if (irank == myrank) continue;
                // Calculate the number of particles which are sent
                const S32 kstart = idx_my_domain_low.z;
                const S32 kend   = idx_my_domain_high.z;
                const S32 jstart = idx_my_domain_low.y;
                const S32 jend   = idx_my_domain_high.y;
                const S32 istart = idx_my_domain_low.x;
                const S32 iend   = idx_my_domain_high.x;
                for (S32 k = kstart; k <= kend; k++) {
                    if ((k < idx_req_low[irank].z) || 
                        (k > idx_req_high[irank].z)) continue;
                    for (S32 j = jstart; j <= jend; j++) {
                        if ((j < idx_req_low[irank].y) ||
                            (j > idx_req_high[irank].z)) continue;
                        for (S32 i = istart; i <= iend; i++) {
                            if ((i < idx_req_low[irank].x) ||
                                (i > idx_req_high[irank].x)) continue;
                            const S32 np = cell_(k,j,i).plist.size();
                            for (S32 n = 0; n < np; n++) {
                               const Particle * ptcl_tmp = cell_(k,j,i).plist[n];
                               ptcl_send[disp] = *ptcl_tmp;
                               disp++;
                            }
                        }
                    }
                }
            }

            // Setup a buffer for recieve
            if (ptcl_recv_ != NULL) delete [] ptcl_recv_;
            n_recv_tot_ = 0;
            for (S32 i = 0; i < nproc; i++) n_recv_tot_ += nrecv[i];
            ptcl_recv_ = new Particle[n_recv_tot_];


            // Peform MPI communication
            Comm::allToAllV(ptcl_send, nsend, nsend_disp, 
                            ptcl_recv_, nrecv, nrecv_disp);


            // Register ptcl_recv_ to cell_
            for (S32 i = 0; i < n_recv_tot_; i++) {
                const S32vec idx = cell_nearest(ptcl_recv_[i].pos, clen_);
                assert(0 <= idx.x && idx.x < NC_);
                assert(0 <= idx.y && idx.y < NC_);
                assert(0 <= idx.z && idx.z < NC_);
                cell_(idx.z,idx.y,idx.x).plist.push_back(&ptcl_recv_[i]);
            } 

            // Release memory
            delete [] ptcl_send;

            Comm::barrier();
            if (Comm::getRank() == 0) std::cout << "PP part: exchange particle (end)" << std::endl;
            //Finalize();
            //std::exit(0);
#endif
        }

        void calcForceParticleMesh() {
            Comm::barrier();
            if (Comm::getRank() == 0) std::cout << "PM part: green function (start)" << std::endl;

            const int NMAX=3;
            const int MMAX=5;
            const double alpha = 2.4;
            GreenFunction<PFMM> gf;
            gf.init(boundary_condition_, NC_, NC_, NC_);
            gf.set(ICUT, clen_, alpha, NMAX, MMAX);
            gf.doFFT();

            Comm::barrier();
            if (Comm::getRank() == 0) {
                std::cout << "PM part: green function (end)" << std::endl;
                std::cout << "PM part: multipole moments (start)" << std::endl;
            }
        
            for(int i=0; i<NC3_; i++){
                cell_[i].do_P2M();
            }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            const int n_mm = (PFMM+1)*(PFMM+1);
            MultipoleMoment<PFMM> *mm_loc, *mm_tot;
            mm_loc = new MultipoleMoment<PFMM>[NC3_];
            mm_tot = new MultipoleMoment<PFMM>[NC3_];
            for (int i = 0; i < NC3_; i++)
                for (int lm = 0; lm < n_mm; lm++)
                    mm_loc[i].buf[lm] = cell_[i].mm.buf[lm];
            const int count = NC3_ * n_mm;
            MPI_Allreduce(mm_loc, mm_tot, count, MPI_DOUBLE,
                          MPI_SUM, MPI_COMM_WORLD);
            for (int i = 0; i < NC3_; i++)
                for (int lm = 0; lm < n_mm; lm++)
                    cell_[i].mm.buf[lm] = mm_tot[i].buf[lm];
#endif
            Comm::barrier();
            if (Comm::getRank() == 0) {
                std::cout << "PM part: multipole moments (end)" << std::endl;
                std::cout << "PM part: M2L convolution (start)" << std::endl;
            }

            M2LConvolution<PFMM,3> (gf, cell_, boundary_condition_);

            Comm::barrier();
            if (Comm::getRank() == 0) {
                std::cout << "PM part: M2L convolution (end)" << std::endl;
                std::cout << "PM part: L2P conversion (start)" << std::endl;
            }

            // Dipole correction
            if (boundary_condition_ == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                Comm::barrier();
                if (Comm::getRank() == 0) {
                    std::cout << "--- applying corrections" << std::endl;
                }

                F64vec dipole(0.0);
                F64 quad0 = 0.0;
                for (S32 i = 0; i < NC3_; i++){
                    dipole.x += cell_[i].mm.buf[3];
                    dipole.y += cell_[i].mm.buf[1];
                    dipole.z += cell_[i].mm.buf[2];
                    quad0 += cell_[i].dispersion();
                }
                const F64 pi = 4.0 * atan(1.0);
                dipole *= (4./3.) * pi;
                printf("quad : %e\n", quad0);
                for (S32 i = 0; i < NC3_; i++){
                    cell_[i].le.buf[3] += 2.0 * dipole.x;
                    cell_[i].le.buf[1] -= 2.0 * dipole.y;
                    cell_[i].le.buf[2] += 1.0 * dipole.z;
                    cell_[i].le.buf[0] += ((2./3.) * pi) * quad0;
                    // self energy correction
                    cell_[i].le.buf[0] -= 
                        alpha * (2.0/sqrt(pi)) * cell_[i].mm.buf[0];
                }
            }
       
            for (S32 i = 0; i < NC3_; i++){
                cell_[i].do_L2P();
                if (boundary_condition_ == BOUNDARY_CONDITION_PERIODIC_XYZ)
                    cell_[i].do_L2P_corr(msum_, alpha);
            }

            Comm::barrier();
            if (Comm::getRank() == 0) {
                std::cout << "PM part: L2P conversion (end)" << std::endl;
            }

        }

        void calcForceSum() {
            const int myrank = Comm::getRank();
            // Compute force sums
            PS::F64vec fpp_loc(0.0), fpm_loc(0.0);
            for(int i=0; i<n_loc_tot_; i++){
                fpp_loc.x += ptcl_[i].mass * ptcl_[i].acc_direct.x;
                fpp_loc.y += ptcl_[i].mass * ptcl_[i].acc_direct.y;
                fpp_loc.z += ptcl_[i].mass * ptcl_[i].acc_direct.z;
                fpm_loc.x += ptcl_[i].mass * ptcl_[i].acc_app.x;
                fpm_loc.y += ptcl_[i].mass * ptcl_[i].acc_app.y;
                fpm_loc.z += ptcl_[i].mass * ptcl_[i].acc_app.z;
            }
            PS::F64vec fpp_tot, fpm_tot;
            fpp_tot = Comm::getSum(fpp_loc);
            fpm_tot = Comm::getSum(fpm_loc);
            if (myrank == 0) {
                printf("PP ftot : (%e, %e, %e)\n", fpp_tot.x, fpp_tot.y, fpp_tot.z);
                printf("PM ftot : (%e, %e, %e)\n", fpm_tot.x, fpm_tot.y, fpm_tot.z);
            }
        
        }

        void checkForce() {
            // Compute the exact force
            if (boundary_condition_ == BOUNDARY_CONDITION_OPEN) {
                // Direct sum (only valid for serial exec.)
#pragma omp parallel for
                for (S32 i=0; i<n_loc_tot_; i++){
                    Particle &pi = ptcl_[i];
                    for (S32 j=0; j<n_loc_tot_; j++){
                        const Particle pj = ptcl_[j];
                        if(j == i) continue;
                        const F64vec dr  = pj.pos - pi.pos;
                        const F64 r2  = dr*dr;
                        const F64 ri2 = 1.0 / r2;
                        const F64 ri  = sqrt(ri2);
                        const F64 ri3 = ri * ri2;
                        pi.phi_direct += pj.mass * ri;
                        pi.acc_direct += (pj.mass * ri3) * dr;
                    }
                }
            } else if (boundary_condition_ == BOUNDARY_CONDITION_PERIODIC_XYZ) {
               // Direct Ewald (only valid for serial exec.)
#ifdef EWALD_FILE
                FILE *fp = fopen("ewald.dat", "r");
                assert(fp);
                int n;
                fscanf(fp, "%d", &n);
                assert(n_loc_tot_ == n);
                for(int i=0; i<n_loc_tot_; i++){
                    fscanf(fp, "%lf %lf %lf %lf",
                            &ptcl_[i].phi_direct,
                            &ptcl_[i].acc_direct.x,
                            &ptcl_[i].acc_direct.y,
                            &ptcl_[i].acc_direct.z);
                }
                fclose(fp);
#else
                const double alpha_ewald = 2.4;
                eval_k_space<5>(n_loc_tot_, alpha_ewald, ptcl_);
                eval_r_space<3>(n_loc_tot_, alpha_ewald, msum_, ptcl_);
#endif
            }

            // Compute and output the sum of potential
            double en_app_loc=0.0, en_dir_loc=0.0;
            for(int i=0; i<n_loc_tot_; i++){
                en_app_loc += 0.5 * ptcl_[i].mass * ptcl_[i].phi_app;
                en_dir_loc += 0.5 * ptcl_[i].mass * ptcl_[i].phi_direct;
            }
            double en_app_tot = Comm::getSum(en_app_loc);
            double en_dir_tot = Comm::getSum(en_dir_loc);
            printf("energy : %24.16e, %24.16e\n", en_app_tot, en_dir_tot);

            // Output the relative differences of potential and force
            if (Comm::getRank() == 0) {
                std::vector<double> err(n_loc_tot_);
                for(int i=0; i<n_loc_tot_; i++) err[i] = ptcl_[i].adiff_rel();
                std::sort(err.begin(), err.end());
                print_err(err, "adiffr", ICUT, PFMM);
        
                for(int i=0; i<n_loc_tot_; i++) err[i] = ptcl_[i].pdiff_rel();
                std::sort(err.begin(), err.end());
                print_err(err, "pdiffr", ICUT, PFMM);
            }

        }


        void calcForceAll() {

            setCell();

            assignParticleToCellLocal();

            calcForceParticleMesh();

            exchangeParticle();

            calcForceParticleParticle<PFMM, ICUT, 3>(NC_, NC_, NC_, cell_,
                                                     boundary_condition_);

            calcForceSum(); // for debug
      
            // Compute the sum of PP and PM forces 
            for(int i=0; i<n_loc_tot_; i++){
                ptcl_[i].move_accp();
            }

            checkForce(); // for debug
        
            Finalize();
            std::exit(0);
       
        }

        template <class Tpsys, class Tdinfo>
        void calcForceAllAndWriteBack(Tpsys & psys,
                                      const Tdinfo & dinfo) {
            setDomainInfoPMMM(dinfo);
            setParticlePMMM(psys);
            calcForceAll();
            for (S32 i = 0; i < n_loc_tot_; i++) {
                F64vec acc;
                acc.x = ptcl_[i].acc_app.x;
                acc.y = ptcl_[i].acc_app.y;
                acc.z = ptcl_[i].acc_app.z;
                F64 pot = ptcl_[i].phi_app;
                psys[i].copyFromForcePMMM(acc, pot);
            }
        }

    };

} // END of namespace of ParticleSimulator
#include "pmmm_impl.hpp"
