#pragma once
#include <complex>
#include "../ps_defs.hpp"
#include "particle_mesh_multipole_defs.hpp"
#include "cell.hpp"
#include "M2L_engine.hpp"

namespace ParticleSimulator {
    namespace ParticleMeshMultipole {

        template <class Tforce, class Tepi>
        class ParticleMeshMultipole {
        public: 
            typedef double real_t;
            typedef std::complex<real_t> cplx_t;
            typedef Cell_FMM<real_t, cplx_t, Tepi, Tforce> Cell_t;

        private:
            bool first_call_by_initialize_;
           
            // Parameters
            Parameters param_;
            S32 p_spj2mm_;
            U32 fftw_planning_rigor_flag_;
            bool use_mpifft_if_possible_;
            S32 fft_size_crit_;

            // Particle information
            S32 n_loc_tot_;
            F64 msum_;
            F64 quad0_;
 
            // Cell information
            S32 n_cell_loc_;
            S32 n_cell_tot_;
            F64vec width_cell_;
            std::vector<Cell_t> cell_loc_;
            std::unordered_map<S32, S32> adr_cell_loc_from_cell_index_;

            // M2L calculation
            M2L_Engine m2l_engine_;

            // Timing measurement
            TimeProfilePMM time_profile_;

        private:

            template <class Tepj, class Tmomloc, class Tmomglb, class Tspj>
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

            template <class Tepj, class Tmomloc, class Tmomglb, class Tspj>
            void setCell(const TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj> & tree) {
                Comm::barrier();
                F64 time_start = GetWtime(); 

                const S32 nx = param_.n_cell.x;
                const S32 ny = param_.n_cell.y;
                const S32 nz = param_.n_cell.z;
                n_cell_tot_ = nx * ny * nz;
                width_cell_ = GetWidthOfParticleMeshCell(param_.pos_unit_cell,
                                                         param_.n_cell);

                const std::vector<S32> list_pm_cell_loc = tree.getListOfLocalParticleMeshCellIndex();
                n_cell_loc_ = list_pm_cell_loc.size();
                //if (Comm::getRank() == 3) std::cout << "n_cell_loc_ = " << n_cell_loc_ << std::endl;
                cell_loc_.resize(n_cell_loc_);
                adr_cell_loc_from_cell_index_.clear();
                for (S32 i=0; i<n_cell_loc_; i++) {
                    const S32 idx = list_pm_cell_loc[i];
                    //if (Comm::getRank() == 3) std::cout << "[#] i = " << i << " idx = " << idx << std::endl;
                    S32vec idx_3d;
                    idx_3d.z = idx / (nx * ny);
                    idx_3d.y = (idx - (nx * ny) * idx_3d.z) / nx;
                    idx_3d.x = idx - (nx * ny) * idx_3d.z - nx * idx_3d.y;
                    const F64vec pos = GetCenterOfParticleMeshCell(param_.pos_unit_cell,
                                                                   width_cell_,
                                                                   idx_3d);
                    cell_loc_[i].clear();
                    cell_loc_[i].init(param_.p);
                    cell_loc_[i].setIdx(idx);
                    cell_loc_[i].setPos(pos);
                    adr_cell_loc_from_cell_index_[idx] = i;
                }

                Comm::barrier();
                time_profile_.set_cell += GetWtime() - time_start;
            }

            template <class Tepj, class Tmomloc, class Tmomglb, class Tspj>
            void setIParticleInfoToCell(TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj> & tree,
                                        const DomainInfo & dinfo) {
                Comm::barrier();
                F64 time_start = GetWtime();

                n_loc_tot_ = tree.getNumberOfEpiSorted();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                for (S32 i=0; i<n_cell_loc_; i++) {
                    const S32 idx = cell_loc_[i].idx;
                    IParticleInfo<Tepi, Tforce> info;
                    if (tree.getIParticleInfoOfParticleMeshCell(idx, info) == 0) {
                        cell_loc_[i].setIParticleInfoToCell(info.n_epi, info.epi_first, info.force_first);
                        // for debug
                        //if (Comm::getRank() == 3 && (idx == 288 || idx == 226)) {
                        //    std::cout << "########" << std::endl;
                        //    std::cout << "i = " << i << " idx = " << idx << std::endl;
                        //    const S32 nx = param_.n_cell.x; 
                        //    const S32 ny = param_.n_cell.y; 
                        //    const S32 nz = param_.n_cell.z;
                        //    const F64ort pos_unit_cell = param_.pos_unit_cell;
                        //    S32vec idx_3d;
                        //    idx_3d.z = idx / (nx * ny);
                        //    idx_3d.y = (idx - (nx * ny) * idx_3d.z) / nx;
                        //    idx_3d.x = idx - (nx * ny) * idx_3d.z - nx * idx_3d.y;
                        //    F64ort box = GetBoxOfParticleMeshCell(idx_3d, pos_unit_cell.low_, width_cell_);
                        //    std::cout << "(ix,iy,iz) = " << idx_3d.x << "   " << idx_3d.y << "   " << idx_3d.z << std::endl;
                        //    std::cout << "box = " << box << std::endl;
                        //    std::cout << "--------" << std::endl;
                        //    for (S32 k=0; k < info.n_epi; k++) {
                        //        const S64 id = info.epi_first[k].getId();
                        //        const F64vec pos = info.epi_first[k].getPos();
                        //        std::cout << "k = " << k
                        //                  << " id = " << id
                        //                  << " pos = " << pos
                        //                  << " T/F = " << box.contains(pos)
                        //                  << std::endl;
                        //    }
                        //}
                    }
                }

                Comm::barrier();
                time_profile_.set_ip_info_to_cell += GetWtime() - time_start;
               
            }

            template <class Tepj, class Tmomloc, class Tmomglb, class Tspj>
            void calcTotalChargeAndDispersion(TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj> & tree,
                                              const DomainInfo & dinfo) {
                Comm::barrier();
                F64 time_start = GetWtime();

                const S32 nx = param_.n_cell.x;
                const S32 ny = param_.n_cell.y;
                const S32 nz = param_.n_cell.z;
                const F64ort pos_my_domain = dinfo.getPosDomain(Comm::getRank());
                F64 msum_loc = 0.0;
                F64 quad0_loc = 0.0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for reduction (+:msum_loc), reduction(+:quad0_loc)
#endif
                for(int i=0; i<n_cell_loc_; i++){
                    const F64vec center = cell_loc_[i].center;
                    const S32 idx = cell_loc_[i].idx;
                    F64 msum, quad0;
                    tree.calcTotalChargeAndDispersionOfParticleMeshCell(idx, center, pos_my_domain, msum, quad0);
                    msum_loc += msum;
                    quad0_loc += quad0;
                }
                msum_ = Comm::getSum(msum_loc);
                quad0_ = Comm::getSum(quad0_loc);

                Comm::barrier();
                time_profile_.calc_msum_and_quad0 += GetWtime() - time_start;
#if 0
                if (Comm::getRank() == 0) {
                    std::cout << "msum_ = " << msum_ << std::endl;
                    std::cout << "quad0_ = " << quad0_ << std::endl;
                }
#endif
            }

            template <class Tepj, class Tmomloc, class Tmomglb, class Tspj>
            void calcMultipoleMoment(TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj> & tree,
                                     const DomainInfo & dinfo) {
                Comm::barrier();
                F64 time_start = GetWtime();

                const S32 n_proc = Comm::getNumberOfProc();
                const S32 my_rank = Comm::getRank();
                const S32 n_thread = Comm::getNumberOfThread();
                MultipoleMoment<real_t, cplx_t> * mm = new MultipoleMoment<real_t, cplx_t>[n_thread];
                for (S32 ith = 0; ith < n_thread; ith++) mm[ith].alloc(param_.p);

                const S32 nx = param_.n_cell.x;
                const S32 ny = param_.n_cell.y;
                const S32 nz = param_.n_cell.z;
                const F64ort pos_my_domain = dinfo.getPosDomain(Comm::getRank());

                // Calculate multipole moments
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for 
#endif
                for (S32 i = 0; i < n_cell_loc_; i++){
                    const S32 ith = Comm::getThreadNum();
                    const S32 idx = cell_loc_[i].idx;
                    const F64vec center = cell_loc_[i].center;
                    mm[ith].clear();
                    if (tree.template calcMultipoleMomentOfParticleMeshCell<real_t, cplx_t>(idx, center, pos_my_domain, mm[ith], p_spj2mm_) == 0) {
                        cell_loc_[i].setMultipoleMoment(mm[ith]);
                    }
                }
                
                for (S32 ith=0; ith<n_thread; ith++) mm[ith].freeMem();
                delete [] mm;

                Comm::barrier();
                time_profile_.calc_multipole_moment += GetWtime() - time_start;
            }

            template <class Tepj, class Tmomloc, class Tmomglb, class Tspj>
            void collectMultipoleMoment(TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj> & tree,
                                        const DomainInfo & dinfo) {
                static S64 n_called = 0;
                n_called++;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                Comm::barrier();
                F64 time_start = GetWtime();

                // Calculate convenient variables
                using pair_t = std::pair<S32, S32>;
                const S32 n_proc = Comm::getNumberOfProc();
                const S32 my_rank = Comm::getRank();
                const S32 LEN = (param_.p + 1) * (param_.p + 1);
                const S32 nx = param_.n_cell.x;
                const S32 ny = param_.n_cell.y;
                const S32 nz = param_.n_cell.z;

                // Collect multipole moments necessary to calculate self-energy correction
                MultientranceArray<S32> rnklst = tree.getRankListFromParticleMeshCellIndex();
                std::vector<pair_t> rnk_idx_pairs_send; // send rank & send cell index
                std::vector<pair_t> rnk_idx_pairs_recv; // recv rank & recv cell index
                for (S32 i=0; i<n_cell_loc_; i++) {
                    const S32 idx = cell_loc_[i].idx;
                    if (cell_loc_[i].is_mm_defined) {
                        const S32 cnt = rnklst.counts[idx];
                        if (cnt > 1) {
                            const S32 disp = rnklst.displs[idx];
                            for (S32 k=1; k<cnt; k++) { // skip k=0 because it is myself.
                                const S32 rnk = rnklst.data[disp + k]; 
                                pair_t tmp;
                                tmp.first = rnk;
                                tmp.second = idx;
                                rnk_idx_pairs_send.push_back(tmp);
                            }
                        }
                    } else {
                        const S32 disp = rnklst.displs[idx];
                        const S32 rnk = rnklst.data[disp];
                        pair_t tmp;
                        tmp.first = rnk;
                        tmp.second = idx;
                        rnk_idx_pairs_recv.push_back(tmp);
                    }
                }
                // Prepare sendbuf
                CommBuffer<S32> idx_send;
                CommBuffer<real_t> mm_send;
                std::unordered_map<S32, S32> cnt_mp;
                for (S32 i=0; i<rnk_idx_pairs_send.size(); i++) {
                    const S32 rnk = rnk_idx_pairs_send[i].first;
                    cnt_mp[rnk]++;
                }
                idx_send.n_comm = cnt_mp.size();
                mm_send.n_comm  = cnt_mp.size();
                idx_send.allocCommInfo();
                mm_send.allocCommInfo();
                S32 adr {0};
                idx_send.count_tot = 0;
                mm_send.count_tot = 0;
                for(auto itr = cnt_mp.begin(); itr != cnt_mp.end(); ++itr) {
                    const S32 rnk = itr->first; // key: rnk
                    const S32 nc = itr->second; // val: # of cells
                    idx_send.ranks[adr] = rnk;
                    idx_send.counts[adr] = nc;
                    idx_send.count_tot += nc;
                    mm_send.ranks[adr] = rnk;
                    mm_send.counts[adr] = LEN * nc;
                    mm_send.count_tot += LEN * nc;
                    adr++;
                }
                idx_send.calcDispls();
                mm_send.calcDispls();
                if (idx_send.count_tot > 0) {
                    idx_send.allocBuffer();
                    mm_send.allocBuffer();
                    idx_send.clearCounts();
                    mm_send.clearCounts();
                    for (S32 i = 0; i < rnk_idx_pairs_send.size(); i++) {
                        const S32 rnk = rnk_idx_pairs_send[i].first;
                        const S32 idx = rnk_idx_pairs_send[i].second;
                        S32 adr_buf {0};
                        for (S32 k=0; k<idx_send.n_comm; k++) {
                            if (rnk == idx_send.ranks[k]) {
                                adr_buf = k;
                            }
                        }
                        {
                            const S32 offset = idx_send.displs[adr_buf] 
                                             + idx_send.counts[adr_buf];
                            idx_send.buf[offset] = idx;
                            idx_send.counts[adr_buf]++;
                        }
                        {
                            S32 offset = mm_send.displs[adr_buf]
                                       + mm_send.counts[adr_buf];
                            S32 adr = adr_cell_loc_from_cell_index_[idx];
                            for (S32 lm=0; lm<LEN; lm++) {
                                mm_send.buf[offset + lm] = cell_loc_[adr].mm.buf[lm];
                            }
                            mm_send.counts[adr_buf] += LEN;
                        }
                    }
                }
                // Prepare recvbuf
                CommBuffer<S32> idx_recv;
                CommBuffer<real_t> mm_recv;
                cnt_mp.clear();
                for (S32 i=0; i<rnk_idx_pairs_recv.size(); i++) {
                    const S32 rnk = rnk_idx_pairs_recv[i].first;
                    cnt_mp[rnk]++;
                }
                idx_recv.n_comm = cnt_mp.size();
                mm_recv.n_comm  = cnt_mp.size();
                idx_recv.allocCommInfo();
                mm_recv.allocCommInfo();
                adr = 0;
                idx_recv.count_tot = 0;
                mm_recv.count_tot = 0;
                for(auto itr = cnt_mp.begin(); itr != cnt_mp.end(); ++itr) {
                    const S32 rnk = itr->first; // key: rnk
                    const S32 nc = itr->second; // val: # of cells
                    idx_recv.ranks[adr] = rnk;
                    idx_recv.counts[adr] = nc;
                    idx_recv.count_tot += nc;
                    mm_recv.ranks[adr] = rnk;
                    mm_recv.counts[adr] = LEN * nc;
                    mm_recv.count_tot += LEN * nc;
                    adr++;
                }
                idx_recv.calcDispls();
                mm_recv.calcDispls();
                idx_recv.allocBuffer();
                mm_recv.allocBuffer();
#if 0
                // Check
                if (n_called == 2) {
                    std::stringstream ss;
                    ss << "mm_comm_info_" << std::setfill('0') << std::setw(5) << my_rank << ".txt";
                    const std::string file_name = ss.str();
                    std::ofstream output_file;
                    output_file.open(file_name.c_str(), std::ios::trunc);
                    output_file << "< idx_send >" << std::endl;
                    idx_send.dumpCommInfo(output_file);
                    output_file << std::endl;
                    output_file << "< idx_recv >" << std::endl;
                    idx_recv.dumpCommInfo(output_file);
                    output_file << std::endl;
                    output_file << "< mm_send >" << std::endl;
                    mm_send.dumpCommInfo(output_file);
                    output_file << std::endl;
                    output_file << "< mm_recv >" << std::endl;
                    mm_recv.dumpCommInfo(output_file);
                    output_file << std::endl;
                    output_file.close();
                    Finalize();
                    std::exit(0);
                }
#endif
                // Perform MPI comm.
                performComm(idx_send, idx_recv);
                performComm(mm_send, mm_recv);
                // Copy the contents of recieve buffers to cell_loc_[].mm.buf[]
                if (idx_recv.n_comm > 0) {
                    for (S32 i = 0; i < idx_recv.count_tot; i++) {
                        const S32 idx = idx_recv.buf[i];
                        const S32 adr = adr_cell_loc_from_cell_index_[idx];
                        const S32 offset = LEN * i;
                        for (S32 lm=0; lm<LEN; lm++) {
                            cell_loc_[adr].mm.buf[lm] = mm_recv.buf[offset + lm];
                        }
                    }
                }

                Comm::barrier();
                F64 time_end = GetWtime();
                time_profile_.collect_multipole_moment += GetWtime() - time_start;
#endif // PARTICLE_SIMULATOR_MPI_PARALLEL
            }

            template <class Tepj, class Tmomloc, class Tmomglb, class Tspj>
            void calcForceParticleMesh(TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj> & tree,
                                       const DomainInfo & dinfo, const bool clear_force = true) {
                // Get data describing which ranks are relavant to each PM cell 
                MultientranceArray<S32> rnklst_from_pm_cell_idx = tree.getRankListFromParticleMeshCellIndex();

                // Determine the way of parallization
                m2l_engine_.initialize(param_, 
                                       fftw_planning_rigor_flag_,
                                       use_mpifft_if_possible_,
                                       fft_size_crit_);

                // Set green function
                m2l_engine_.setGreenFunction();

                // Collect information of MM necessary to perform local task
                m2l_engine_.redistMM(rnklst_from_pm_cell_idx,
                                     n_cell_loc_, cell_loc_);

                // Perform M2L transform
                m2l_engine_.convolution();

                // Collect information of LE necessary to evaluate the potential and
                // its gradient at the positions of local particles.
                m2l_engine_.redistLE(rnklst_from_pm_cell_idx,
                                     n_cell_loc_, adr_cell_loc_from_cell_index_, cell_loc_);
                
                Comm::barrier();
                F64 time_start = GetWtime();

                // Dipole correction
                if (param_.bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {

                    F64 dipole_x, dipole_y, dipole_z;
                    dipole_x = dipole_y = dipole_z = 0.0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for reduction(+:dipole_x), reduction(+:dipole_y), reduction(+:dipole_z)
#endif
                    for (S32 i = 0; i < n_cell_loc_; i++){
                        if (cell_loc_[i].is_mm_defined) {
                            dipole_x += cell_loc_[i].mm.buf[3];
                            dipole_y += cell_loc_[i].mm.buf[1];
                            dipole_z += cell_loc_[i].mm.buf[2];
                        }
                    }
                    dipole_x = Comm::getSum(dipole_x);
                    dipole_y = Comm::getSum(dipole_y);
                    dipole_z = Comm::getSum(dipole_z);
                    F64vec dipole = F64vec(dipole_x, dipole_y, dipole_z);
                    const F64 pi = 4.0 * atan(1.0);
                    dipole *= (4./3.) * pi;
#if 0
                    std::cout << "msum_  = " << msum_ << std::endl;
                    std::cout << "quad0_ = " << quad0_ << std::endl;
                    std::cout << "dipole = " << dipole << std::endl;
#endif
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for (S32 i = 0; i < n_cell_loc_; i++){
                        cell_loc_[i].le.buf[3] += 2.0 * dipole.x;
                        cell_loc_[i].le.buf[1] -= 2.0 * dipole.y;
                        cell_loc_[i].le.buf[2] += 1.0 * dipole.z;
                        cell_loc_[i].le.buf[0] += ((2./3.) * pi) * quad0_;
                        // self energy correction
                        cell_loc_[i].le.buf[0] -= 
                            param_.alpha * (2.0/sqrt(pi)) * cell_loc_[i].mm.buf[0];
                    }
                }

#if 0
                // Check le_r after dipole correction
                static S32 n_called {0};
                n_called++;
                if (n_called == 2) {
                    std::stringstream ss;
                    ss << "le_r_" << std::setfill('0') << std::setw(5)
                       << Comm::getRank() << ".txt";
                    const std::string filename = ss.str();
                    std::ofstream output_file;
                    output_file.open(filename.c_str(), std::ios::trunc);
                    for (S32 i = 0; i < n_cell_loc_; i++) {
                        const S32 idx_1d = cell_loc_[i].idx;
                        const S32 nx = param_.n_cell.x;
                        const S32 ny = param_.n_cell.y;
                        const S32 nz = param_.n_cell.z;
                        S32vec idx_3d;
                        idx_3d.z = idx_1d / (nx * ny);
                        idx_3d.y = (idx_1d - (nx * ny) * idx_3d.z) / nx;
                        idx_3d.x = idx_1d - (nx * ny) * idx_3d.z - nx * idx_3d.y;
                        const S32 LEN = (param_.p + 1) * (param_.p + 1);
                        for (S32 lm = 0; lm < LEN; lm++) {
                            const S32 idx = lm
                                          + LEN * (idx_3d.x
                                          + nx * (idx_3d.y
                                          + ny * idx_3d.z));
                            const real_t val = cell_loc_[i].le.buf[lm];
                            output_file << idx << "    " << val << std::endl;
                        }
                    }
                    output_file.close();
                    Finalize();
                    std::exit(0);
                }
#endif

                Comm::barrier();
                time_start = GetWtime();

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for 
#endif          
                for (S32 i = 0; i < n_cell_loc_; i++){
                    cell_loc_[i].do_L2P(clear_force);
                    if (param_.bc == BOUNDARY_CONDITION_PERIODIC_XYZ) {
                        cell_loc_[i].do_L2P_corr(msum_, param_.alpha);
                    }
                }

                Comm::barrier();
                time_profile_.L2P += GetWtime() - time_start;
                //Finalize();
                //std::exit(0);
            }


        public:

            ParticleMeshMultipole() {
                 first_call_by_initialize_ = true;
                 msum_ = 0.0;
                 quad0_ = 0.0;
            }
            ~ParticleMeshMultipole() {}

            void initialize(const S32 p, 
                            const S32 p_spj2mm,
                            const U32 fftw_planning_rigor_flag = FFTW_MEASURE,
                            const bool use_mpifft_if_possible = true,
                            const S32 fft_size_crit = 10000) {
                // Check error
                assert(first_call_by_initialize_);
                assert(p >= 0);
                assert(p_spj2mm >= 0);
                if ((fftw_planning_rigor_flag != FFTW_ESTIMATE) &&
                    (fftw_planning_rigor_flag != FFTW_MEASURE)  &&
                    (fftw_planning_rigor_flag != FFTW_PATIENT)  &&
                    (fftw_planning_rigor_flag != FFTW_EXHAUSTIVE)) {
                    if (Comm::getRank() == 0) {
                        std::stringstream err_msg;
                        err_msg << "The value of argument fftw_planning_rigor_flag must be\n"
                                << "one of the following:\n"
                                << "--------------------------\n"
                                << "    FFTW_ESTIMATE   (" << FFTW_ESTIMATE   << ")\n"
                                << "    FFTW_MEASURE    (" << FFTW_MEASURE    << ")\n"
                                << "    FFTW_PATIENT    (" << FFTW_PATIENT    << ")\n"
                                << "    FFTW_EXHAUSTIVE (" << FFTW_EXHAUSTIVE << ")\n"
                                << "--------------------------\n"
                                << "On the other hand, the value you specified for\n"
                                << "fftw_planning_rigor_flag was " << fftw_planning_rigor_flag << ".\n"
                                << "The following are possible causes:\n"
                                << "(1) You just passed an invalid value.\n"
                                << "(2) The value is valid but its position in the argument list\n"
                                << "    is wrong.";
                        std::cout << err_msg.str() << std::endl;
                    }
                    Abort(-1);
                }
                // Make tables for Rlm class
                {
                    Rlm<real_t> tmp;
                    tmp.make_table(p); // for P2M & L2P
                }
                // [Note (tag: #144eddb3)]
                //   (1) Tables required by Slm class are created when needed.
                //       (see transform() in M2L_Engine.hpp & 
                //        calcGreenFunctionInRealSpace() in green_function.hpp)
                //   (2) In the current implementation, we do not need to care about p_spj2mm
                //       because MomentMultipoleGeometricCenterPMMM in tree.hpp is implemented
                //       using fmm_org.hpp, where p is a template parameter and hence necessary
                //       tables are made at the compilation time.

                // Store parameters
                first_call_by_initialize_ = false;
                param_.p = p;
                p_spj2mm_ = p_spj2mm;
                fftw_planning_rigor_flag_ = fftw_planning_rigor_flag;
                use_mpifft_if_possible_ = use_mpifft_if_possible;
                fft_size_crit_ = fft_size_crit;
            }

            template <class Tepj, class Tmomloc, class Tmomglb, class Tspj>
            void calcForceAll(TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj> & tree,
                              const DomainInfo & dinfo,
                              const bool clear_force = true) {
                //Comm::barrier(); if (Comm::getRank() == 0) std::cout << "OK1 at " << __FILE__ << std::endl;
                setParam(tree, dinfo);

                //Comm::barrier(); if (Comm::getRank() == 0) std::cout << "OK2 at " << __FILE__ << std::endl;
                setCell(tree);

                //Comm::barrier(); if (Comm::getRank() == 0) std::cout << "OK3 at " << __FILE__ << std::endl;
                setIParticleInfoToCell(tree, dinfo);

                //Comm::barrier(); if (Comm::getRank() == 0) std::cout << "OK4 at " << __FILE__ << std::endl;
                calcTotalChargeAndDispersion(tree, dinfo);

                //Comm::barrier(); if (Comm::getRank() == 0) std::cout << "OK5 at " << __FILE__ << std::endl;
                calcMultipoleMoment(tree, dinfo);

                //Comm::barrier(); if (Comm::getRank() == 0) std::cout << "OK6 at " << __FILE__ << std::endl;
                collectMultipoleMoment(tree, dinfo);

                //Comm::barrier(); if (Comm::getRank() == 0) std::cout << "OK7 at " << __FILE__ << std::endl;
                calcForceParticleMesh(tree, dinfo, clear_force);

            }

            template <class Tepj, class Tmomloc, class Tmomglb, class Tspj, class Tpsys>
            void calcForceAllAndWriteBack(TreeForForce<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj> & tree, 
                                          const DomainInfo & dinfo,
                                          Tpsys & psys,
                                          const bool clear_force = true) {
                Comm::barrier();
                F64 time_start = GetWtime();

                calcForceAll(tree, dinfo, clear_force);
                tree.copyForceOriginalOrder();
                for(S32 i=0; i<n_loc_tot_; i++)
                    psys[i].copyFromForcePMM(tree.getForce(i));

                Comm::barrier();
                time_profile_.calc_force_all_and_write_back += GetWtime() - time_start;
            }

            void clearTimeProfile() {
                time_profile_.clear();
                m2l_engine_.clearTimeProfile();
            }

            TimeProfilePMM getTimeProfile() const {
                return time_profile_ + m2l_engine_.getTimeProfile();
            }

        };
    } // END of namespace of ParticleMeshMultipole
} // END of namespace of ParticleSimulator
