////////////////////////////////////////////////
/// implementaion of methods of TreeForForce ///

#include"tree_walk.hpp"

namespace ParticleSimulator{

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    std::vector<S32> TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getListOfParticleMeshCellIndex(const S32 i) const {
        assert(i < Comm::getNumberOfProc());
        return list_pm_cell_for_mm_[i];
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    std::vector<S32> TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getListOfSharedParticleMeshCellIndex() const {
        return list_pm_cell_shared_for_mm_;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    std::vector<S32> TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getListOfLocalParticleMeshCellIndex() const {
        return list_pm_cell_loc_;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    MultientranceArray<S32> TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getRankListFromParticleMeshCellIndex() const {
        return rnklst_from_pm_cell_idx_;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getIParticleInfoOfParticleMeshCell(const S32 idx, IParticleInfo<Tepi, Tforce> & info) {
        auto itr = ip_info_pm_cell_.find(idx);
        if ( itr != ip_info_pm_cell_.end()) {
            info = ip_info_pm_cell_[idx];
            // [Note (tag: #fe161b03)]
            //     In order to add a const qulifier to this function,
            //     we have to use std::unordered_map::at instead of operator [].
            //     However, Fujitsu C++ compiler does not support `at`.
            //     Hence, we gave up on making this function a const function.
            //     (see also Note #cee527db in particle_mesh_multipole/M2L_engine.hpp)
            return 0;
        } else {
            return -1;
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template <class real_t, class cplx_t>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMultipoleMomentOfParticleMeshCell(const S32 idx, 
                                          const F64vec & center,
                                          const F64ort & pos_my_domain,
                                          MultipoleMoment<real_t, cplx_t> & mm,
                                          const S32 p_spj2mm) const {
        S32 ret = -1;
        const std::vector<S32> list_my_pm_cell = list_pm_cell_for_mm_[Comm::getRank()];
        auto itr_outer = std::find(list_my_pm_cell.begin(),
                                   list_my_pm_cell.end(),
                                   idx);
        if ( itr_outer != list_my_pm_cell.end() ) {
            auto itr_inner = adr_tc_glb_from_pm_cell_idx_.find(idx);
            if ( itr_inner != adr_tc_glb_from_pm_cell_idx_.end() ) {
                const S32 adr_tc = itr_inner->second;
                const S32 adr_tree_sp_first = spj_sorted_.size() - tc_glb_.size();
                F64 r_crit_sq;
                if (theta_ > 0.0) {
                    r_crit_sq = (length_ * length_) / (theta_ * theta_);
                } else {
                    r_crit_sq = -1.0;
                }
                S32vec idx_3d;
                idx_3d.z = idx / (n_cell_.x * n_cell_.y);
                idx_3d.y = (idx - (n_cell_.x * n_cell_.y) * idx_3d.z) / n_cell_.x;
                idx_3d.x = idx - (n_cell_.x * n_cell_.y) * idx_3d.z - n_cell_.x * idx_3d.y;
                CalcMultipoleMomentOfParticleMeshCell
                    <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, Tspj,
                     real_t, cplx_t>
                    (idx_3d, adr_tc, tc_glb_, tp_glb_, epj_sorted_, spj_sorted_,
                     n_leaf_limit_, lev_leaf_limit_, adr_tree_sp_first, r_crit_sq,
                     pos_unit_cell_.low_, icut_,
                     level_pm_cell_, width_pm_cell_, center, pos_my_domain,
                     mm, p_spj2mm);
                ret = 0;
            }
        }
        return ret;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcTotalChargeAndDispersionOfParticleMeshCell(const S32 idx, 
                                                   const F64vec & center,
                                                   const F64ort & pos_my_domain,
                                                   F64 & msum,
                                                   F64 & quad0) const {
        msum = quad0 = 0;
        auto itr = adr_tc_glb_from_pm_cell_idx_.find(idx);
        if ( itr != adr_tc_glb_from_pm_cell_idx_.end() ) {
            const S32 adr_tc = itr->second;
            CalcTotalChargeAndDispersionOfParticleMeshCell
                <TreeCell<Tmomglb>, TreeParticle, Tepj>
                (adr_tc, tc_glb_, tp_glb_, epj_sorted_,
                 center, pos_my_domain, msum, quad0);
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    decomposeParticleMeshCellImpl2(const DomainInfo & dinfo) {
        using pair_t = std::pair<S32, S32>; // (idx,np) or (rnk, np)
        const S32 my_rank = Comm::getRank(); 
        const S32 n_proc  = Comm::getNumberOfProc();
        const S32 n_cell_tot = n_cell_.x * n_cell_.y * n_cell_.z;

        // Initialize 
        if (first_call_decomposeParticleMeshCell_) {
            list_pm_cell_for_mm_ = new std::vector<S32>[n_proc];
            first_call_decomposeParticleMeshCell_ = false;
        } else {
            for (S32 i=0; i<n_proc; i++)
                list_pm_cell_for_mm_[i].clear();
        }

        // Make a list of pairs of a cell index and number of EPJs in that cell
        // and a list of cell indices of PM cells in which local EPIs exist 
        // (list_pm_cell_loc_). 
        ReallocatableArray<pair_t> * idx_np_pairs;
        idx_np_pairs = new ReallocatableArray<pair_t>[n_proc];
        list_pm_cell_loc_.resize(0); // reset
        for (const auto & e: ip_info_pm_cell_) {
            pair_t tmp;
            tmp.first = e.first; // idx
            tmp.second = e.second.n_epi; // np
            idx_np_pairs[my_rank].push_back(tmp);
            list_pm_cell_loc_.push_back(e.first);
        }
        const S32 sendcount = idx_np_pairs[my_rank].size();
        S32 * recvcounts = new S32[n_proc];
        Comm::allGather(&sendcount, 1, recvcounts);
        S32 recvcount_tot = 0;
        for (S32 i=0; i<n_proc; i++) recvcount_tot += recvcounts[i];
        ReallocatableArray<pair_t> recvbuf;
        recvbuf.resizeNoInitialize(recvcount_tot);
        S32 * recvdispls = new S32[n_proc];
        recvdispls[0] = 0;
        for (S32 i=1; i<n_proc; i++)
            recvdispls[i] = recvdispls[i-1] + recvcounts[i-1];
        Comm::allGatherV(idx_np_pairs[my_rank].getPointer(), sendcount,
                         recvbuf.getPointer(), recvcounts, recvdispls);
        S32 offset = 0;
        for (S32 i=0; i < n_proc; i++) {
            if (i != my_rank) {
                for (S32 k=0; k < recvcounts[i]; k++) 
                    idx_np_pairs[i].push_back(recvbuf[offset + k]);
            }
            offset += recvcounts[i];
        }
        recvbuf.freeMem(1);
        delete [] recvcounts;
        delete [] recvdispls;

        // Set rnklst_from_pm_cell_idx_ and relevant variables 
        rnklst_from_pm_cell_idx_.counts.resize(n_cell_tot);
        for (S32 i = 0; i < n_cell_tot; i++) rnklst_from_pm_cell_idx_.counts[i] = 0;
        for (S32 i = 0; i < n_proc; i++) {
            for (S32 k = 0; k < idx_np_pairs[i].size(); k++) {
                const S32 idx = idx_np_pairs[i][k].first; // idx
                rnklst_from_pm_cell_idx_.counts[idx]++;
            }
        }
        rnklst_from_pm_cell_idx_.displs.resize(n_cell_tot);
        rnklst_from_pm_cell_idx_.displs[0] = 0;
        for (S32 i = 1; i < n_cell_tot; i++)
            rnklst_from_pm_cell_idx_.displs[i] = rnklst_from_pm_cell_idx_.displs[i-1]
                                               + rnklst_from_pm_cell_idx_.counts[i-1];
        for (S32 i = 0; i < n_cell_tot; i++) rnklst_from_pm_cell_idx_.counts[i] = 0;
        rnklst_from_pm_cell_idx_.data.resize(recvcount_tot);
        for (S32 i = 0; i < n_proc; i++) {
            for (S32 k = 0; k < idx_np_pairs[i].size(); k++) {
                const S32 idx = idx_np_pairs[i][k].first; // idx
                const S32 adr = rnklst_from_pm_cell_idx_.displs[idx]
                              + rnklst_from_pm_cell_idx_.counts[idx];
                rnklst_from_pm_cell_idx_.data[adr] = i;
                rnklst_from_pm_cell_idx_.counts[idx]++;
            }
        }
#if 0
        // Check
        if (Comm::getRank() == 0) {
            const std::string filename = "rnklst_from_pm_cell_idx_.txt";
            std::ofstream output_file;
            output_file.open(filename.c_str(), std::ios::trunc);
            for (S32 idx=0; idx < n_cell_tot; idx++) {
                const S32 cnt  = rnklst_from_pm_cell_idx_.counts[idx];
                const S32 disp = rnklst_from_pm_cell_idx_.displs[idx];
                output_file << "idx = " << idx
                            << " cnt = " << cnt
                            << " disp = " << disp
                            << std::endl;
                for (S32 i=0; i < cnt; i++) {
                    const S32 rnk = rnklst_from_pm_cell_idx_.data[i + disp];
                    output_file << "    " << rnk << std::endl;
                }
            }
            output_file.close();
        }
        Finalize();
        std::exit(0);
#endif

 
        // Set rnk_np_pairs to calculate list_pm_cell_for_mm_ and etc.
        ReallocatableArray<pair_t> * rnk_np_pairs;
        rnk_np_pairs = new ReallocatableArray<pair_t>[n_cell_tot];
        for (S32 i = 0; i < n_proc; i++) {
            for (S32 k = 0; k < idx_np_pairs[i].size(); k++) {
                const S32 idx = idx_np_pairs[i][k].first; // idx
                const S32 np = idx_np_pairs[i][k].second; // np
                pair_t tmp;
                tmp.first = i; // rnk
                tmp.second = np; // np
                rnk_np_pairs[idx].push_back(tmp);
            }
        }

        // Free memory of unused variables
        for (S32 i=0; i<n_proc; i++) idx_np_pairs[i].freeMem(1);
        delete [] idx_np_pairs;

        // Calculate a list of indice of PM cells for which this process
        // is responsible for the calculation of multipole moment.
#if !defined(PARTICLE_SIMULATOR_USE_EPJ_ONLY_TO_EVAL_MM_IN_PMM)
        // In this case, we use SPJ of a process that is responsible for PM cell
        // of interest to evaluate multipole moments of that PM cell.
        // The process responsible for a PM cell is the process that the number
        // of local particles contained in that PM cell is the largest.
        for (S32 idx = 0; idx < n_cell_tot; idx++) {
            // Find a process having the most particles in a cell
            S32 rnk_max = -1;
            S32 np_max = 0;
            for (S32 i = 0; i < rnk_np_pairs[idx].size(); i++) {
                const S32 rnk = rnk_np_pairs[idx][i].first; // rnk 
                const S32 np = rnk_np_pairs[idx][i].second; // np
                if (np > np_max) {
                    np_max = np; rnk_max = rnk;
                }
            }
            if (rnk_max != -1) {
                list_pm_cell_for_mm_[rnk_max].push_back(idx);
            }
            // Change the order of ranks in rnklst_from_pm_cell_idx_
            const S32 cnt = rnklst_from_pm_cell_idx_.counts[idx];
            if (cnt > 1) {
                // Examine the location of rnk_max in rnklst_from_pm_cell_idx_
                const S32 disp = rnklst_from_pm_cell_idx_.displs[idx];
                S32 adr = disp;
                for (S32 i=0; i<cnt; i++) {
                    const S32 rnk = rnklst_from_pm_cell_idx_.data[adr];
                    if (rnk == rnk_max) break;
                    else adr++;
                }
                // Swap data 
                if (adr != disp) {
                    const S32 rnk = rnklst_from_pm_cell_idx_.data[disp];
                    rnklst_from_pm_cell_idx_.data[disp] = rnk_max;
                    rnklst_from_pm_cell_idx_.data[adr]  = rnk;
                    //std::cout << "swap: " << adr << " <-> " << disp << std::endl;
                }
            }
        }
#else
        assert(false); // not yet determined 
        // In this case, we do not use SPJ to evaluate multipole moments in
        // the PMM module. Instead, we use EPJ only to evaluate them.
        // Because epj_sorted_ consists of both local particles and particles
        // owned by other processes, we need to exclude other processes' 
        // particles in the calculation of multipole moments.
        list_pm_cell_shared_for_mm_.clear();
        for (S32 idx = 0; idx < n_cell_tot; idx++) {
            // Find a process having the most particles in a cell
            const S32 n_rnk = rnk_np_pairs[idx].size();
            S32 rnk_max = -1;
            S32 np_max = 0;
            for (S32 i = 0; i < n_rnk; i++) {
                const S32 rnk = rnk_np_pairs[idx][i].first; // rnk
                const S32 np = rnk_np_pairs[idx][i].second; // np
                list_pm_cell_for_mm_[rnk].push_back(idx);
                if (np > np_max) {
                    np_max = np; rnk_max = rnk;
                }
            }
            if (n_rnk > 1) {
                list_pm_cell_shared_for_mm_.push_back(idx);
            }
            // Change the order of ranks in rnklst_from_pm_cell_idx_
            const S32 cnt = rnklst_from_pm_cell_idx_.counts[idx];
            if (cnt > 1) {
                // Examine the location of rnk_max in rnklst_from_pm_cell_idx_
                const S32 disp = rnklst_from_pm_cell_idx_.displs[idx];
                S32 adr = disp;
                for (S32 i=0; i<cnt; i++) {
                    const S32 rnk = rnklst_from_pm_cell_idx_.data[adr];
                    if (rnk == rnk_max) break;
                    else adr++;
                }
                // Swap data 
                if (adr != disp) {
                    const S32 rnk = rnklst_from_pm_cell_idx_.data[disp];
                    rnklst_from_pm_cell_idx_.data[disp] = rnk_max;
                    rnklst_from_pm_cell_idx_.data[adr]  = rnk;
                }
            }
        }
#endif

        // Free
        delete [] rnk_np_pairs;

#if 0
        // Check
        std::ostringstream ss;
        ss << "list_pm_cell_" 
           << std::setfill('0') << std::setw(5) << my_rank << ".txt";
        std::string filename = ss.str();
        std::ofstream output_file;
        output_file.open(filename.c_str(), std::ios::trunc);
        for (S32 i=0; i < list_pm_cell_for_mm_[my_rank].size(); i++) {
            const S32 idx = list_pm_cell_for_mm_[my_rank][i];
            output_file << idx << std::endl;
        }
        output_file.close();
#endif
#if 0
        if (my_rank == 0) std::cout << "decomposeParticleMeshCellImpl2 ended" << std::endl;
        Finalize();
        std::exit(0);
#endif
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    Tepj * TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getEpjFromId(const S64 id, const Tepj * epj_tmp){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
        {
            if(map_id_to_epj_.empty()){
                S64 n_epj = epj_sorted_.size();
                for(S32 i=0; i<n_epj; i++){
                    if(GetMSB(tp_glb_[i].adr_ptcl_) == 1) continue;
                    Tepj * epj_tmp = epj_sorted_.getPointer(i);
                    S64 id_tmp = epj_tmp->getId();
                    map_id_to_epj_.insert( std::pair<S64, Tepj*>(id_tmp, epj_tmp) );
                }
            }
        }
        Tepj * epj = NULL;
        typename MyMap::iterator it = map_id_to_epj_.find(id);
        if( it != map_id_to_epj_.end() ) epj = it->second;
        return epj;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    size_t TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>
    ::getMemSizeUsed() const {
        size_t tmp = 0;
        for(int i=0; i<Comm::getNumberOfThread(); i++){
            tmp =
                epj_for_force_[i].getMemSize() + spj_for_force_[i].getMemSize()
                + epjr_send_buf_[i].getMemSize() + epjr_send_buf_for_scatter_[i].getMemSize() 
                + epjr_recv_1st_sorted_[i].getMemSize();
        }
        return tmp 
            + tp_glb_.getMemSize()
            + tc_loc_.getMemSize()
            + tc_glb_.getMemSize()
            + epi_sorted_.getMemSize() + epi_org_.getMemSize()
            + epj_sorted_.getMemSize() + epj_org_.getMemSize()
            + spj_sorted_.getMemSize() + spj_org_.getMemSize()
            + ipg_.getMemSize()
            + epj_send_.getMemSize()
            + spj_send_.getMemSize()
            + force_sorted_.getMemSize() + force_org_.getMemSize()
            + epjr_sorted_.getMemSize() + epjr_send_.getMemSize()
            + epjr_recv_.getMemSize() + epjr_recv_1st_buf_.getMemSize()
            + epjr_recv_2nd_buf_.getMemSize();
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    size_t TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>
    ::getUsedMemorySize() const {
        return getMemSizeUsed();
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    initializeRemaining() {
        lev_group_limit_ = lev_leaf_limit_  = 0;
        map_id_to_epj_.clear();
        lev_max_loc_ = lev_max_glb_ = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        const S64 n_proc = Comm::getNumberOfProc();
        n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = n_walk_local_ = 0;
        n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
        wtime_exlet_comm_ = wtime_exlet_a2a_ = wtime_exlet_a2av_ = 0.0;
        wtime_walk_LET_1st_ = wtime_walk_LET_2nd_ = 0.0;
        ni_ave_ = nj_ave_ = 0;
        Comm::barrier();

        /*
        epi_sorted_.setAllocMode(1); // msortLT --- final
        spj_sorted_.setAllocMode(1); //  --- final
        epi_org_.setAllocMode(1); // setPtclLT ---
        spj_org_.setAllocMode(1); // insted of it, use spj_recv
        epj_org_.setAllocMode(1);
        force_sorted_.setAllocMode(1); // -- final
        epj_send_.setAllocMode(1);
        spj_send_.setAllocMode(1);
        
        epj_for_force_    = new ReallocatableArray<Tepj>[n_thread];
        spj_for_force_    = new ReallocatableArray<Tspj>[n_thread];
        for(S32 i=0; i<n_thread; i++){
            epj_for_force_[i].setAllocMode(1);
            spj_for_force_[i].setAllocMode(1);
        }
        */

#ifdef PARTICLE_SIMULATOR_USE_MEMORY_POOL
        epi_sorted_.setAllocMode(1); // msortLT --- final
        spj_sorted_.setAllocMode(1); //  --- final
        epi_org_.setAllocMode(1); // setPtclLT ---
        spj_org_.setAllocMode(1); // insted of it, use spj_recv
        epj_org_.setAllocMode(1);
        force_sorted_.setAllocMode(1); // -- final
        epj_send_.setAllocMode(1);
        spj_send_.setAllocMode(1);
        epj_for_force_    = new ReallocatableArray<Tepj>[n_thread];
        spj_for_force_    = new ReallocatableArray<Tspj>[n_thread];
        for(S32 i=0; i<n_thread; i++){
            epj_for_force_[i].setAllocMode(1);
            spj_for_force_[i].setAllocMode(1);
        }
#else
        epi_org_.setAllocMode(1);
        epj_send_.setAllocMode(1);
        spj_send_.setAllocMode(1);
        epj_for_force_    = new ReallocatableArray<Tepj>[n_thread];
        spj_for_force_    = new ReallocatableArray<Tspj>[n_thread];
        for(S32 i=0; i<n_thread; i++){
            epj_for_force_[i].setAllocMode(1);
            spj_for_force_[i].setAllocMode(1);
        }
#endif

        
        adr_epj_for_force_ = new ReallocatableArray<S32>[n_thread];
        adr_spj_for_force_ = new ReallocatableArray<S32>[n_thread];
        adr_ipg_for_force_ = new ReallocatableArray<S32>[n_thread];
        
        epjr_send_buf_ = new ReallocatableArray<EPJWithR>[n_thread];
        epjr_send_buf_for_scatter_ = new ReallocatableArray<EPJWithR>[n_thread];
        epjr_recv_1st_sorted_ = new ReallocatableArray<EPJWithR>[n_thread];
        epj_neighbor_ = new ReallocatableArray<Tepj>[n_thread];
        n_cell_open_ = new CountT[n_thread];
        Comm::barrier();
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        req_send_ = new MPI_Request[n_proc];
        req_recv_ = new MPI_Request[n_proc];
        status_   = new MPI_Status[n_proc];
#endif
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        if(Comm::getRank() == 0){
            std::cerr<<"used mem size for tree= "<<this->getMemSizeUsed()*1e-9<<"[GB]"<<std::endl;
        }
#endif
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    initialize(const U64 n_glb_tot,
               const F64 theta,
               const U32 n_leaf_limit,
               const U32 n_group_limit) {
        bool err = false;
        if (typeid(TSM) == typeid(SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE)) {
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The argument specification you are using is incorrect.");
        }
        if(is_initialized_ == true){
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("Do not initialize the tree twice");
        }
        if( typeid(TSM) == typeid(SEARCH_MODE_LONG) || 
            typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF) ){
            if(theta < 0.0){
                err = true;
                PARTICLE_SIMULATOR_PRINT_ERROR("The opening criterion of the tree must be >= 0.0");
                std::cout<<"theta = " << theta << std::endl;
            }
        }
        if(n_leaf_limit <= 0){
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The limit number of the particles in the leaf cell must be > 0");
            std::cout<<"n_leaf_limit = " << n_leaf_limit << std::endl;
        }
        if(n_group_limit < n_leaf_limit){
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The limit number of particles in ip graoups msut be >= that in leaf cells");
            std::cout<<"n_group_limit = " << n_group_limit << std::endl;
            std::cout<<"n_leaf_limit = " << n_leaf_limit << std::endl;
        }
        if (err) {
            std::cerr<<"SEARCH_MODE: " << typeid(TSM).name() << std::endl;
            std::cerr<<"Force: " << typeid(Tforce).name() << std::endl;
            std::cerr<<"EPI: " << typeid(Tepi).name() << std::endl;
            std::cerr<<"EPJ: " << typeid(Tepj).name() << std::endl;
            std::cerr<<"SPJ: " << typeid(Tspj).name() << std::endl;
            Abort(-1);
        }
        is_initialized_ = true;
        n_glb_tot_ = n_glb_tot;
        theta_ = theta;
        n_leaf_limit_ = n_leaf_limit;
        n_group_limit_ = n_group_limit;
        initializeRemaining();
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    initialize(const U64 n_glb_tot,
               const S32vec n_cell,
               const S32 icut,
               const F64 theta,
               const U32 n_leaf_limit,
               const U32 n_group_limit) {
        bool err = false;
        if (typeid(TSM) != typeid(SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE)) {
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The argument specification you are using is incorrect.");
        }
        if (is_initialized_ == true){
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("Do not initialize the tree twice");
        }
        if (icut < 1) {
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The cell separation must be >= 1");
            std::cout << "icut = " << icut << std::endl;
        }
        if ((n_cell.x <= 0) || (n_cell.y <= 0) || (n_cell.z <= 0)) {
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The number of the cells in each dimension must be > 0");
            std::cout << "n_cell = " << n_cell << std::endl;
        }
        if (theta < 0.0){
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The opening criterion of the tree must be >= 0.0");
            std::cout << "theta = " << theta << std::endl;
        }
        if (n_leaf_limit <= 0){
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The limit number of the particles in the leaf cell must be > 0");
            std::cout << "n_leaf_limit = " << n_leaf_limit << std::endl;
        }
        if (n_group_limit < n_leaf_limit){
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The limit number of particles in ip groups must be >= that in leaf cells");
            std::cout << "n_leaf_limit  = " << n_leaf_limit << std::endl;
            std::cout << "n_group_limit = " << n_group_limit << std::endl;
        }
        if (err) {
            std::cerr<<"SEARCH_MODE: " << typeid(TSM).name() << std::endl;
            std::cerr<<"Force: " << typeid(Tforce).name() << std::endl;
            std::cerr<<"EPI: " << typeid(Tepi).name() << std::endl;
            std::cerr<<"EPJ: " << typeid(Tepj).name() << std::endl;
            std::cerr<<"SPJ: " << typeid(Tspj).name() << std::endl;
            Abort(-1);
        }
        is_initialized_ = true;
        n_glb_tot_ = n_glb_tot;
        icut_ = icut;
        n_cell_ = n_cell;
        theta_ = theta;
        n_leaf_limit_ = n_leaf_limit;
        n_group_limit_ = n_group_limit; 

        initializeRemaining();
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tpsys>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setParticleLocalTreeImpl(const Tpsys & psys,
                             const bool clear){
        const F64 time_offset = GetWtime();
        const S32 nloc = psys.getNumberOfParticleLocal();
        if(clear){ n_loc_tot_ = 0;}
        //        const S32 offset = 0;
        const S32 offset = n_loc_tot_;
        n_loc_tot_ += nloc;
        epj_org_.resizeNoInitialize(n_loc_tot_);
        epi_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<nloc; i++){
                epi_org_[i].copyFromFP( psys[i] );
                epj_org_[i].copyFromFP( psys[i] );
            }
        }
        else{
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<nloc; i++){
                epi_org_[i+offset].copyFromFP( psys[i] );
                epj_org_[i+offset].copyFromFP( psys[i] );
            }
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"nloc="<<nloc<<std::endl;
        std::cout<<"n_loc_tot_="<<n_loc_tot_<<std::endl;
#endif
        time_profile_.set_particle_local_tree += GetWtime() - time_offset;
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setRootCell(const DomainInfo & dinfo){
        const F64 time_offset = GetWtime();
        calcCenterAndLengthOfRootCell(typename TSM::search_type(), dinfo);
        // set pos_unit_cell_ for PMM
        if(dinfo.getBoundaryCondition() == BOUNDARY_CONDITION_OPEN){
            pos_unit_cell_ = pos_root_cell_;
            //pos_unit_cell_ = F64ort(F64vec(0.0), F64vec(1.0)); // for test
        }
        else{
            pos_unit_cell_ = dinfo.getPosRootDomain();
        }
        repositionPosRootCellToEmbedParticleMesh(typename TSM::search_type());
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
            PARTICLE_SIMULATOR_PRINT_LINE_INFO();
            std::cout<<"length_="<<length_<<" center_="<<center_<<std::endl;
            std::cout<<"pos_root_cell_="<<pos_root_cell_<<std::endl;
#endif
        time_profile_.set_root_cell += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setRootCell(const F64 l, const F64vec & c){
        const F64 time_offset = GetWtime();
        center_ = c;
        length_ = l;
        pos_root_cell_.low_ = center_ - F64vec(length_*0.5);
        pos_root_cell_.high_ = center_ + F64vec(length_*0.5);
        time_profile_.set_root_cell += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    mortonSortLocalTreeOnly(const bool reuse){
        const F64 wtime_offset = GetWtime();
        epi_sorted_.resizeNoInitialize(n_loc_tot_);
        epj_sorted_.resizeNoInitialize(n_loc_tot_);
        adr_org_from_adr_sorted_loc_.resizeNoInitialize(n_loc_tot_);
        tp_glb_.resizeNoInitialize(n_loc_tot_);
        if(!reuse){
            //tp_loc_.resizeNoInitialize(n_loc_tot_);
            //tp_glb_.resizeNoInitialize(n_loc_tot_);
            ReallocatableArray<TreeParticle> tp_buf(n_loc_tot_, n_loc_tot_, 1);
            MortonKey<typename TSM::key_type>::initialize(length_ * 0.5,
                                                          center_,
                                                          pos_unit_cell_.low_,
                                                          idx_unit_cell_,
                                                          level_pm_cell_,
                                                          width_pm_cell_);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                //tp_loc_[i].setFromEP(epj_org_[i], i);
                tp_glb_[i].template setFromEP<Tepj, typename TSM::key_type>(epj_org_[i], i);
            }
#ifdef USE_STD_SORT
            //std::sort(tp_loc_.getPointer(), tp_loc_.getPointer()+n_loc_tot_, 
            //          [](const TreeParticle & l, const TreeParticle & r )
            //          ->bool{return l.getKey() < r.getKey();} );
            std::sort(tp_glb_.getPointer(), tp_glb_.getPointer()+n_loc_tot_, 
                      [](const TreeParticle & l, const TreeParticle & r )
                      ->bool{return l.getKey() < r.getKey();} );            

#else
            //rs_.lsdSort(tp_loc_.getPointer(), tp_buf.getPointer(), 0, n_loc_tot_-1);
            rs_.lsdSort(tp_glb_.getPointer(), tp_buf.getPointer(), 0, n_loc_tot_-1);
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                //const S32 adr = tp_loc_[i].adr_ptcl_;
                const S32 adr = tp_glb_[i].adr_ptcl_;
                adr_org_from_adr_sorted_loc_[i] = adr;
            }
            tp_buf.freeMem(1);
        } // end of if() no reuse
        
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_loc_tot_; i++){
            const S32 adr = adr_org_from_adr_sorted_loc_[i];
            epi_sorted_[i] = epi_org_[adr];
            epj_sorted_[i] = epj_org_[adr];
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epi_sorted_.size()="<<epi_sorted_.size()<<" epj_sorted_.size()="<<epj_sorted_.size()<<std::endl;
#endif
        time_profile_.morton_sort_local_tree += GetWtime() - wtime_offset;
        time_profile_.make_local_tree += GetWtime() - wtime_offset;

        /*
        if(Comm::getRank()==0){
            std::cerr<<"epi_sorted_.size()= "<<epi_sorted_.size()
                     <<" epi_org_.size()= "<<epi_org_.size()
                     <<" epj_sorted_.size()= "<<epj_sorted_.size()
                     <<" epj_org_.size()= "<<epj_org_.size()
                     <<std::endl;
            for(S32 i=0; i<epi_sorted_.size(); i += epi_sorted_.size()/10){
                std::cerr<<"i= "<<i
                         <<" epi_sorted_[i].id= "<<epi_sorted_[i].id
                         <<" epi_sorted_[i].pos= "<<epi_sorted_[i].pos
                         <<std::endl;
            }
            for(S32 i=0; i<epj_sorted_.size(); i += epj_sorted_.size()/10){
                std::cerr<<"i= "<<i
                         <<"epj_sorted_[i].id= "<<epj_sorted_[i].id
                         <<" epj_sorted_[i].pos= "<<epj_sorted_[i].pos
                         <<std::endl;
            }
        }
        */

#if 0
        // for debug
        std::string filename;
        std::ostringstream ss;
        std::ofstream output_file;
        ss << "epi_sorted_" << std::setfill('0')  << std::setw(5) << Comm::getRank() << ".txt";
        filename = ss.str();
        output_file.open(filename.c_str(), std::ios::trunc);
        output_file.setf(std::ios_base::scientific,
                         std::ios_base::floatfield);
        output_file << std::setprecision(15) << std::showpos;
        output_file << "# x, y, z" << std::endl;
        for(S32 i=0; i<epi_sorted_.size(); i++){
            output_file << epi_sorted_[i].pos << std::endl;
        }
        output_file.close();
        Finalize(); std::exit(0);
#endif
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    linkCellLocalTreeOnly(){
        const F64 time_offset = GetWtime();
        LinkCell<TreeCell< Tmomloc >, typename TSM::key_type>
            (tc_loc_,  adr_tc_level_partition_loc_,
             tp_glb_.getPointer(), lev_max_loc_, n_loc_tot_,
             n_leaf_limit_, lev_leaf_limit_);
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"tc_loc_.size()="<<tc_loc_.size()<<std::endl;
        std::cout<<"lev_max_loc_="<<lev_max_loc_<<std::endl;
#endif
        time_profile_.link_cell_local_tree += GetWtime() - time_offset;
        time_profile_.make_local_tree += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnly(){
        F64 time_offset = GetWtime();
        calcMomentLocalTreeOnlyImpl(typename TSM::search_type());
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif
        time_profile_.calc_moment_local_tree += GetWtime() - time_offset;
        time_profile_.make_local_tree_tot = time_profile_.calc_moment_local_tree + time_profile_.make_local_tree;
    }


    // FOR P^3T
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLongScatter){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_,
                   n_leaf_limit_, lev_leaf_limit_);
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLongSymmetry){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_,
                   n_leaf_limit_, lev_leaf_limit_);
    }
    
    // FOR P^3T
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLongCutoffScatter){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_,
                   n_leaf_limit_, lev_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLong){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_,
                   n_leaf_limit_, lev_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLongCutoff){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_,
                   n_leaf_limit_, lev_leaf_limit_);
    }

    // For PM^3
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLongParticleMeshMultipole){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_,
                   n_leaf_limit_, lev_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchShortScatter){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_,
                   n_leaf_limit_, lev_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchShortGather){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epi_sorted_.getPointer(), lev_max_loc_,
                   n_leaf_limit_, lev_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchShortSymmetry){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epi_sorted_.getPointer(), lev_max_loc_,
                   n_leaf_limit_, lev_leaf_limit_);
    }

    ///////////////////////////////
    /// morton sort global tree ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    mortonSortGlobalTreeOnly(const bool reuse){
        F64 time_offset = GetWtime();
        //ReallocatableArray< TreeParticle > tp_buf;
        if(!map_id_to_epj_.empty()){
            map_id_to_epj_.clear();
        }
        assert(map_id_to_epj_.empty());
        tp_glb_.resizeNoInitialize(n_glb_tot_);
        const S32 n_ep_tot = epj_org_.size();
        epj_sorted_.resizeNoInitialize(n_ep_tot);
        adr_org_from_adr_sorted_glb_.resizeNoInitialize(n_glb_tot_);
        if(!reuse){
            ReallocatableArray<TreeParticle> tp_buf(n_glb_tot_, n_glb_tot_, 1);
            //tp_buf.resizeNoInitialize(n_glb_tot_);
#ifdef USE_STD_SORT
            std::sort(tp_glb_.getPointer(), tp_glb_.getPointer()+n_glb_tot_, 
                      [](const TreeParticle & l, const TreeParticle & r )
                      ->bool{return l.getKey() < r.getKey();} );
#else
            rs_.lsdSort(tp_glb_.getPointer(), tp_buf.getPointer(), 0, n_glb_tot_-1);
#endif

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel 
#endif
            for(S32 i=0; i<n_glb_tot_; i++){
                adr_org_from_adr_sorted_glb_[i] = tp_glb_[i].adr_ptcl_;
            }
            tp_buf.freeMem(1);
        }
        if( typeid(typename TSM::force_type) == typeid(TagForceLong)){
            //const S32 n_sp_tot = spj_recv_.size();
            const S32 n_sp_tot = spj_org_.size();
            assert(n_ep_tot+n_sp_tot == n_glb_tot_);
            
            spj_sorted_.resizeNoInitialize(n_sp_tot);

#if 0
            // thread parallelized, but not fast.
            const S32 n_thread = Comm::getNumberOfThread();
            static ReallocatableArray<U32> n_cnt_ep;
            static ReallocatableArray<U32> n_cnt_sp;
            n_cnt_ep.resizeNoInitialize(n_thread);
            n_cnt_sp.resizeNoInitialize(n_thread);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel 
#endif
            {
                S32 id_thread = Comm::getThreadNum();
                n_cnt_ep[id_thread] = n_cnt_sp[id_thread] = 0;
                U32 id_tp_head = (n_glb_tot_/n_thread) * id_thread + std::min(id_thread, (S32)n_glb_tot_%n_thread);
                U32 id_tp_tail = (n_glb_tot_/n_thread) * (id_thread+1) + std::min( (id_thread+1), (S32)n_glb_tot_%n_thread);
                for(U32 i=id_tp_head; i<id_tp_tail; i++){
                    const U32 adr = adr_org_from_adr_sorted_glb_[i];
                    //const U32 adr = tp_glb_[i].adr_ptcl_;
                    if( GetMSB(adr) == 0){
                        n_cnt_ep[id_thread]++;
                    }
                    else{
                        n_cnt_sp[id_thread]++;
                    }
                }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp barrier
#endif
                U32 id_ep_head = 0;
                U32 id_sp_head = 0;
                for(U32 i=0; i<id_thread; i++){
                    id_ep_head += n_cnt_ep[i];
                    id_sp_head += n_cnt_sp[i];
                }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp barrier
#endif
                n_cnt_ep[id_thread] = n_cnt_sp[id_thread] = 0;
                for(U32 i=id_tp_head; i<id_tp_tail; i++){
                    const U32 adr = adr_org_from_adr_sorted_glb_[i];
                    //const U32 adr = tp_glb_[i].adr_ptcl_;
                    if( GetMSB(adr) == 0){
                        const U32 id_dst = id_ep_head + n_cnt_ep[id_thread];
                        epj_sorted_[id_dst] = epj_org_[adr];
                        tp_glb_[i].adr_ptcl_ = id_dst;
                        n_cnt_ep[id_thread]++;
                    }
                    else{
                        const U32 id_dst = id_sp_head + n_cnt_sp[id_thread];
                        spj_sorted_[id_dst] = spj_org_[ClearMSB(adr)];
                        tp_glb_[i].adr_ptcl_ = SetMSB(id_dst);
                        n_cnt_sp[id_thread]++;
                    }
                }
            }
#else
            U32 n_cnt_ep = 0;
            U32 n_cnt_sp = 0;
            for(S32 i=0; i<n_glb_tot_; i++){
                const U32 adr = adr_org_from_adr_sorted_glb_[i];
                if( GetMSB(adr) == 0){
                    epj_sorted_[n_cnt_ep] = epj_org_[adr];
                    tp_glb_[i].adr_ptcl_ = n_cnt_ep;
                    n_cnt_ep++;
                }
                else{
                    spj_sorted_[n_cnt_sp] = spj_org_[ClearMSB(adr)];
                    tp_glb_[i].adr_ptcl_ = SetMSB(n_cnt_sp);
                    n_cnt_sp++;
                }
            }
#endif
        }
        else{
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_glb_tot_; i++){
                const U32 adr = adr_org_from_adr_sorted_glb_[i];
                epj_sorted_[i] = epj_org_[adr];
                tp_glb_[i].adr_ptcl_ = i;
            }
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"tp_glb_.size()="<<tp_glb_.size()<<std::endl;
        std::cout<<"epj_sorted_.size()="<<epj_sorted_.size()<<" spj_sorted_.size()="<<spj_sorted_.size()<<std::endl;
#endif
        time_profile_.morton_sort_global_tree += GetWtime() - time_offset;
        time_profile_.make_global_tree += GetWtime() - time_offset;
    }

    /////////////////////////////
    /// link cell global tree ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    linkCellGlobalTreeOnly(){
        const F64 time_offset = GetWtime();
        LinkCell<TreeCell< Tmomloc >, typename TSM::key_type>
            (tc_glb_, adr_tc_level_partition_glb_,
             tp_glb_.getPointer(), lev_max_glb_, n_glb_tot_,
             n_leaf_limit_, lev_leaf_limit_);
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"tc_glb_.size()="<<tc_glb_.size()<<std::endl;
        std::cout<<"lev_max_glb_="<<lev_max_glb_<<std::endl;
#endif
        time_profile_.link_cell_global_tree += GetWtime() - time_offset;
        time_profile_.make_global_tree += GetWtime() - time_offset;
    }

    
    //////////////////////////
    // CALC MOMENT GLOBAL TREE
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnly(){
        const F64 time_offset = GetWtime();
        calcMomentGlobalTreeOnlyImpl(typename TSM::search_type());
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif
        time_profile_.calc_moment_global_tree += GetWtime() - time_offset;
        time_profile_.make_global_tree_tot = time_profile_.calc_moment_global_tree + time_profile_.make_global_tree;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchLong){
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_, lev_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchLongScatter){
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_, lev_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchLongSymmetry){
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_, lev_leaf_limit_);
    }    

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchLongCutoff){
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_, lev_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchLongParticleMeshMultipole){
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_, lev_leaf_limit_);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchShortScatter){
        CalcMoment(adr_tc_level_partition_glb_, tc_glb_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_glb_,
                   n_leaf_limit_, lev_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchShortGather){
        CalcMoment(adr_tc_level_partition_glb_, tc_glb_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_glb_,
                   n_leaf_limit_, lev_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchShortSymmetry){
        CalcMoment(adr_tc_level_partition_glb_, tc_glb_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_glb_,
                   n_leaf_limit_, lev_leaf_limit_);
    }
    
    ////////////////////
    /// MAKE IPGROUP ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeIPGroup(){
        const F64 time_offset = GetWtime();
        ipg_.clearSize();
        makeIPGroupImpl(typename TSM::force_type());
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"ipg_.size()="<<ipg_.size()<<std::endl;
#endif
        time_profile_.calc_force += GetWtime() - time_offset;
        time_profile_.calc_force__make_ipgroup += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeIPGroupImpl(TagForceLong){
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"n_group_limit_="<<n_group_limit_<<std::endl;
#endif

#ifdef PARTICLE_SIMULATOR_GLB_TREE_CELL_AS_IPG_BOX
        MakeIPGroupLongGLBTreeCellAsIPGBox(ipg_, tc_loc_, tc_glb_, epj_sorted_, 0, 0,
                                           n_group_limit_, lev_group_limit_,
                                           n_leaf_limit_, lev_leaf_limit_);
#else //PARTICLE_SIMULATOR_GLB_TREE_CELL_AS_IPG_BOX
#if 1
        MakeIPGroupUseGLBTreeLong(ipg_, tc_loc_, tc_glb_, epi_sorted_, 0, 0,
                                  n_group_limit_, lev_group_limit_,
                                  n_leaf_limit_, lev_leaf_limit_);
#else
        MakeIPGroupLong(ipg_, tc_loc_, epi_sorted_, 0,
                        n_group_limit_, lev_group_limit_);
#endif

        const S32 n_ipg = ipg_.size();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif //PARTICLE_SIMULATOR_THREAD_PARALLEL
        for(S32 i=0; i<n_ipg; i++){
            const S32 n = ipg_[i].n_ptcl_;
            const S32 adr = ipg_[i].adr_ptcl_;
            ipg_[i].vertex_in_ = GetMinBoxSingleThread(epi_sorted_.getPointer(adr), n);
        }
#endif //PARTICLE_SIMULATOR_GLB_TREE_CELL_AS_IPG_BOX
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeIPGroupImpl(TagForceShort){
        MakeIPGroupShort(ipg_, tc_loc_, epi_sorted_, 0,
                         n_group_limit_, lev_group_limit_);
    }
    

    /////////////////////////////
    /// MAKE INTERACTION LIST ///
    // pj is epj or spj
    template<class Tpj>
    void CopyPjForForceST(const ReallocatableArray<S32> & adr_pj,
                          const ReallocatableArray<Tpj> & pj_sorted,
                          ReallocatableArray<Tpj> & pj_for_force){
        const S32 n_pj = adr_pj.size();
        pj_for_force.resizeNoInitialize(n_pj);
        for(S32 i=0; i<n_pj; i++){
            const S32 adr_pj_src = adr_pj[i];
            pj_for_force[i] = pj_sorted[adr_pj_src];
        }
    }
    template<class Tpj>
    void CopyPjForForceST(const ReallocatableArray<S32> & adr_pj,
                          const ReallocatableArray<Tpj> & pj_sorted,
                          const S32 n_head,
                          const S32 n_tail,
                          ReallocatableArray<Tpj> & pj_for_force){
        pj_for_force.resizeNoInitialize(n_tail);
        for(S32 i=n_head; i<n_tail; i++){
            const S32 adr_pj_src = adr_pj[i];
            pj_for_force[i] = pj_sorted[adr_pj_src];
        }
    }
    
    template<class T>
    struct TraitsForCutoff{
        typedef TagWithoutCutoff type_cutoff;
    };
    template<>
    struct TraitsForCutoff<TagSearchLongCutoff>{
        typedef TagWithCutoff type_cutoff;
    };
    template<>
    struct TraitsForCutoff<TagSearchLongParticleMeshMultipole>{
        typedef TagWithCutoffByCell type_cutoff;
    };
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListLongForZeroTheta(TagWithoutCutoff, const S32 adr_ipg){
        const S32 ith = Comm::getThreadNum();
        const S32 n_tmp = tc_glb_[0].n_ptcl_;
        S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
        epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
        for(S32 ip=0; ip<n_tmp; ip++){
            if( GetMSB(tp_glb_[adr_ptcl_tmp].adr_ptcl_) == 0){
                epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp++]);
            }
            else {
                adr_ptcl_tmp++;
            }
        }
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListLongForZeroTheta(TagWithCutoff, const S32 adr_ipg){
        const S32 ith = Comm::getThreadNum();
        if( !tc_glb_[0].isLeaf(n_leaf_limit_,lev_leaf_limit_) ){
            const F64 r_cut_sq  = epj_sorted_[0].getRSearch() * epj_sorted_[0].getRSearch();
            const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_in_;
            const F64ort cell_box = pos_root_cell_;
            MakeInteractionListLongCutoffEPForZeroTheta
                (tc_glb_, tc_glb_[0].adr_tc_, tp_glb_, 
                 epj_sorted_, epj_for_force_[ith],
                 cell_box, pos_target_box, r_cut_sq,
                 n_leaf_limit_, lev_leaf_limit_); 
        }
        else{
            const S32 n_tmp = tc_glb_[0].n_ptcl_;
            S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
            epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
            for(S32 ip=0; ip<n_tmp; ip++){
                if( GetMSB(tp_glb_[adr_ptcl_tmp].adr_ptcl_) == 0){
                    epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp++]);
                }
                else {
                    adr_ptcl_tmp++;
                }
            }
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListLongForZeroTheta(TagWithCutoffByCell, const S32 adr_ipg){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = GetMinBoxRoundedUpToParticleMeshCell((ipg_[adr_ipg]).vertex_in_,
                                                                           pos_unit_cell_.low_,
                                                                           width_pm_cell_, icut_);
        if( !tc_glb_[0].isLeaf(n_leaf_limit_, lev_leaf_limit_) ){
            MakeInteractionListLongCutoffByCellEPForZeroTheta
                (tc_glb_, tc_glb_[0].adr_tc_, tp_glb_, 
                 epj_sorted_, epj_for_force_[ith],
                 pos_target_box, n_leaf_limit_, lev_leaf_limit_); 
        }
        else{
            const S32 n_tmp = tc_glb_[0].n_ptcl_;
            S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
            epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
            for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
                if( GetMSB(tp_glb_[adr_ptcl_tmp].adr_ptcl_) == 0){
                    if (pos_target_box.contains(epj_sorted_[adr_ptcl_tmp].getPos())) {
                        epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp]);
                    }
                }
            }
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionList(const S32 adr_ipg, const bool clear){
        makeInteractionListImpl(typename TSM::force_type(), adr_ipg, clear);
    }

    static void CalcNumShift(const F64ort root_domain,
                             const F64ort my_outer_boundary,
                             const F64vec shift,
                             S32 & n_shift){
        F64ort root_domain_tmp = root_domain;
        root_domain_tmp.high_ += shift;
        root_domain_tmp.low_ += shift;
        //std::cerr<<"root_domain_tmp= "<<root_domain_tmp
        //         <<" my_outer_boundary= "<<my_outer_boundary
        //         <<std::endl;
        while(my_outer_boundary.overlapped(root_domain_tmp)){
            root_domain_tmp.high_ += shift;
            root_domain_tmp.low_ += shift;
            n_shift++;
        }        
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tep2>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcCenterAndLengthOfRootCellImpl(const Tep2 ep[],
                                      const DomainInfo & dinfo){
        F64ort box_loc;
        box_loc.initNegativeVolume();
        if (typeid(TSM) == typeid(SEARCH_MODE_LONG)) {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
            {
                F64ort box_loc_tmp;
                box_loc_tmp.init();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for nowait
#endif
                for(S32 ip=0; ip<n_loc_tot_; ip++){
                    box_loc_tmp.merge(ep[ip].getPos());
                }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
                {
                    box_loc.merge(box_loc_tmp);
                }
            }
        } else {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
            {
                F64ort box_loc_tmp;
                box_loc_tmp.init();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for nowait
#endif
                for(S32 ip=0; ip<n_loc_tot_; ip++){
                    box_loc_tmp.merge(ep[ip].getPos(), GetMyRSearch(ep[ip])*1.000001);
                }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp critical
#endif
                {
                    box_loc.merge(box_loc_tmp);
                }
            }
        }
        F64vec xlow_loc  = box_loc.low_;
        F64vec xhigh_loc = box_loc.high_;
        F64vec xlow_glb  = Comm::getMinValue(xlow_loc);
        F64vec xhigh_glb = Comm::getMaxValue(xhigh_loc);
        const F64ort my_outer_boundary(xlow_glb, xhigh_glb);
        const F64ort root_domain = dinfo.getPosRootDomain();
        const F64vec shift = root_domain.high_ - root_domain.low_;
        bool pa[DIMENSION_LIMIT];
        dinfo.getPeriodicAxis(pa);
        S32 num_shift_p[DIMENSION_LIMIT];
        S32 num_shift_n[DIMENSION_LIMIT];
        S32 num_shift_max[DIMENSION_LIMIT];
        //std::cerr<<"my_outer_boundary= "<<my_outer_boundary<<std::endl;
        for(S32 cid=0; cid<DIMENSION_LIMIT; cid++){
            if (pa[cid]) {
                num_shift_p[cid] = num_shift_n[cid] = 0;
                F64vec shift_tmp(0.0);
                shift_tmp[cid] = shift[cid];
                CalcNumShift(root_domain, my_outer_boundary, shift_tmp, num_shift_p[cid]);
                CalcNumShift(root_domain, my_outer_boundary, -shift_tmp, num_shift_n[cid]);
                //std::cerr<<"num_shift_p[cid]= "<<num_shift_p[cid]
                //         <<" num_shift_n[cid]= "<<num_shift_n[cid]
                //         <<std::endl;
                num_shift_max[cid] = std::max(num_shift_p[cid], num_shift_n[cid]);
            } else {
                num_shift_max[cid] = 0;
            }
        }
        length_ = 0.0;
        for(S32 cid=0; cid<DIMENSION_LIMIT; cid++){
            if (pa[cid]) {
                F64 length_tmp = (2*num_shift_max[cid]+1)*shift[cid];
                if(length_tmp > length_) length_ = length_tmp;
                center_[cid] = root_domain.getCenter()[cid];
            } else {
                F64 length_tmp = my_outer_boundary.getFullLength()[cid];
                if (length_tmp > length_) length_ = length_tmp;
                center_[cid] = my_outer_boundary.getCenter()[cid];
            }
        }
        length_ *= 1.000001;
        pos_root_cell_.low_ = center_ - F64vec(length_*0.5);
        pos_root_cell_.high_ = center_ + F64vec(length_*0.5);
        //std::cerr<<"pos_root_cell_= "<<pos_root_cell_<<std::endl;
    }

    inline static S32 getMinExpOfTwoGE(const S32 in) { // GE := equal to or greater than
        assert(in > 0);
        S32 tmp = 1, lev = 0;
        while (tmp < in) {
            tmp *= 2;
            lev++;
        }
        return lev;
    }

    inline static S32 getMinExpOfTwoGT(const S32 in) { // GT := greater than
        assert(in > 0);
        S32 tmp = 1, lev = 0;
        while (tmp <= in) {
            tmp *= 2;
            lev++;
        }
        return lev;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    repositionPosRootCellToEmbedParticleMeshImpl(){
        // Compute cell widths
        width_pm_cell_ = GetWidthOfParticleMeshCell(pos_unit_cell_, n_cell_);

        // Compute tree level 
        S32vec n_cell_margin, n_cell_tot_reqmin, level;
        n_cell_margin.x = (icut_/n_cell_.x + 1)*n_cell_.x;
        n_cell_margin.y = (icut_/n_cell_.y + 1)*n_cell_.y;
        n_cell_margin.z = (icut_/n_cell_.z + 1)*n_cell_.z;
        n_cell_tot_reqmin.x = n_cell_.x + 2 * n_cell_margin.x;
        n_cell_tot_reqmin.y = n_cell_.y + 2 * n_cell_margin.y;
        n_cell_tot_reqmin.z = n_cell_.z + 2 * n_cell_margin.z;
        // [Notes]
        // (1) In order to get a correct PM cell ID in the exchange LET,
        //     we must require that a shifted domain is fully contained
        //     in pos_root_cell. 
        // (2) The following is an old implementation:
        //     n_cell_tot_reqmin.x = n_cell_.x + 2*icut_;
        //     n_cell_tot_reqmin.y = n_cell_.y + 2*icut_;
        //     n_cell_tot_reqmin.z = n_cell_.z + 2*icut_;
        while(1) {
            F64vec length;
            length.x = width_pm_cell_.x * n_cell_tot_reqmin.x;
            length.y = width_pm_cell_.y * n_cell_tot_reqmin.y;
            length.z = width_pm_cell_.z * n_cell_tot_reqmin.z;
            bool is_large_enough = true;
            if (length.x < pos_root_cell_.getFullLength().x) {
                n_cell_tot_reqmin.x++; is_large_enough=false;
            }
            if (length.y < pos_root_cell_.getFullLength().y) {
                n_cell_tot_reqmin.y++; is_large_enough=false;
            }
            if (length.z < pos_root_cell_.getFullLength().z) {
                n_cell_tot_reqmin.z++; is_large_enough=false;
            }
            if (is_large_enough) break;
        }
        level.x = getMinExpOfTwoGT(n_cell_tot_reqmin.x);
        level.y = getMinExpOfTwoGT(n_cell_tot_reqmin.y);
        level.z = getMinExpOfTwoGT(n_cell_tot_reqmin.z);
        level_pm_cell_ = std::max(std::max(level.x, level.y), level.z);
        lev_leaf_limit_ = lev_group_limit_ = level_pm_cell_;

        S32vec n_cell_tot;
        n_cell_tot.x = 1<<level_pm_cell_; // 2^{level_pm_cell_}
        n_cell_tot.y = 1<<level_pm_cell_;
        n_cell_tot.z = 1<<level_pm_cell_;
        

        // Compute idx_unit_cell_, pos_root_cell_
        idx_unit_cell_.x = (n_cell_tot.x - n_cell_tot_reqmin.x)/2 + n_cell_margin.x;
        idx_unit_cell_.y = (n_cell_tot.y - n_cell_tot_reqmin.y)/2 + n_cell_margin.y;
        idx_unit_cell_.z = (n_cell_tot.z - n_cell_tot_reqmin.z)/2 + n_cell_margin.z;
        assert(idx_unit_cell_.x > (icut_ - 1));
        assert(idx_unit_cell_.y > (icut_ - 1));
        assert(idx_unit_cell_.z > (icut_ - 1));
        pos_root_cell_.low_.x = pos_unit_cell_.low_.x - width_pm_cell_.x * idx_unit_cell_.x;
        pos_root_cell_.low_.y = pos_unit_cell_.low_.y - width_pm_cell_.y * idx_unit_cell_.y;
        pos_root_cell_.low_.z = pos_unit_cell_.low_.z - width_pm_cell_.z * idx_unit_cell_.z;
        pos_root_cell_.high_.x = pos_root_cell_.low_.x + width_pm_cell_.x * n_cell_tot.x;
        pos_root_cell_.high_.y = pos_root_cell_.low_.y + width_pm_cell_.y * n_cell_tot.y;
        pos_root_cell_.high_.z = pos_root_cell_.low_.z + width_pm_cell_.z * n_cell_tot.z;

        // Recalculate length_ & center_
        const F64vec lenvec = pos_root_cell_.getFullLength();
        length_ = std::max(std::max(lenvec.x, lenvec.y), lenvec.z);
        center_ = pos_root_cell_.getCenter();

#if 0
        std::cout << "pos_unit_cell_ = " << pos_unit_cell_ << std::endl;
        std::cout << "idx_unit_cell_ = " << idx_unit_cell_ << std::endl;
        std::cout << "(n_cell_tot      = " << n_cell_tot << ")" << std::endl;
        std::cout << "level_pm_cell_   = " << level_pm_cell_   << std::endl;
        std::cout << "width_pm_cell_   = " << width_pm_cell_   << std::endl;
        std::cout << "pos_root_cell_   = " << pos_root_cell_   << std::endl;
        std::cout << "pos_root_cell_.getFullLength() = " << pos_root_cell_.getFullLength() << std::endl;
        std::cout << "length_ " << length_ << std::endl;
        std::cout << "center_ " << center_ << std::endl;
#endif
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    copyForceOriginalOrder(){
        epi_org_.freeMem(1);
        force_org_.resizeNoInitialize(n_loc_tot_);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_loc_tot_; i++){
            const S32 adr = adr_org_from_adr_sorted_loc_[i];
            force_org_[adr] = force_sorted_[i];
        }
    }
    
    // return forces in original order
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceDirect(Tfunc_ep_ep pfunc_ep_ep,
                    Tforce force[],
                    const DomainInfo & dinfo,
                    const bool clear){
        if(clear){
            for(S32 i=0; i<n_loc_tot_; i++) force[i].clear();
        }
        Tepj * epj_tmp;
        S32 n_epj_tmp;
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        AllGatherParticle(epj_tmp, n_epj_tmp, epj_org_.getPointer(), n_loc_tot_, dinfo.getPosRootDomain(), pos_root_cell_, pa);
        pfunc_ep_ep(epi_org_.getPointer(), n_loc_tot_, epj_tmp, n_epj_tmp, force);
        delete [] epj_tmp;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceDirectAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                const DomainInfo & dinfo,
                                const bool clear){
        force_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
            for(S32 i=0; i<n_loc_tot_; i++)force_org_[i].clear();
        }
        Tepj * epj_tmp;
        S32 n_epj_tmp;
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        AllGatherParticle(epj_tmp, n_epj_tmp, epj_org_.getPointer(), n_loc_tot_, dinfo.getPosRootDomain().getFullLength(), pos_root_cell_, pa);
        pfunc_ep_ep(epi_org_.getPointer(), n_loc_tot_, epj_tmp, n_epj_tmp, force_org_.getPointer());
        delete [] epj_tmp;
    }

#if 1
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){
        return getNeighborListOneParticleImpl(typename TSM::neighbor_search_type(), ptcl, epj);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagNeighborSearchScatter, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        SearchNeighborListOneParticleScatter(pos_target,    tc_glb_.getPointer(),
                                             tp_glb_.getPointer(), adr, 
                                             epj_sorted_,   epj_neighbor_[id_thread],
                                             n_leaf_limit_, lev_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagNeighborSearchGather, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        F64 r_search_sq = ptcl.getRSearch() * ptcl.getRSearch();
        const S32 adr = 0;
        SearchNeighborListOneParticleGather(pos_target,
                                            r_search_sq,
                                            tc_glb_.getPointer(),
                                            tp_glb_.getPointer(), adr, 
                                            epj_sorted_, epj_neighbor_[id_thread],
                                            n_leaf_limit_, lev_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagNeighborSearchSymmetry, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        F64 r_search_sq = ptcl.getRSearch() * ptcl.getRSearch();
        const S32 adr = 0;
        SearchNeighborListOneParticleSymmetry(pos_target,
                                              r_search_sq,
                                              tc_glb_.getPointer(),
                                              tp_glb_.getPointer(), adr, 
                                              epj_sorted_, epj_neighbor_[id_thread],
                                              n_leaf_limit_, lev_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }
    #ifdef PARTICLE_SIMULATOR_CHECK_SEARCH_MODE
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagNeighborSearchNo, const Tptcl & ptcl, Tepj * & epj){
        return -1;
        // std::cerr<<"not implemented"<<std::endl;
    }
    #endif
#else
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){
        return getNeighborListOneParticleImpl(typename TSM::search_type(), ptcl, epj);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagSearchShortScatter, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        SearchNeighborListOneParticleScatter(pos_target,    tc_glb_.getPointer(),
                                             tp_glb_.getPointer(), adr, 
                                             epj_sorted_,   epj_neighbor_[id_thread],
                                             n_leaf_limit_, lev_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagSearchLongScatter, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        SearchNeighborListOneParticleScatter(pos_target, tc_glb_.getPointer(),
                                             tp_glb_.getPointer(), adr, 
                                             epj_sorted_, epj_neighbor_[id_thread],
                                             n_leaf_limit_, lev_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagSearchShortSymmetry, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        F64 r_search_sq = ptcl.getRSearch() * ptcl.getRSearch();
        const S32 adr = 0;
        SearchNeighborListOneParticleSymmetry(pos_target,
                                              r_search_sq,
                                              tc_glb_.getPointer(),
                                              tp_glb_.getPointer(), adr, 
                                              epj_sorted_, epj_neighbor_[id_thread],
                                              n_leaf_limit_, lev_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagSearchLongSymmetry, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        const F64 r_search_sq = ptcl.getRSearch() * ptcl.getRSearch();
        SearchNeighborListOneParticleSymmetry(pos_target,
                                              r_search_sq,
                                              tc_glb_.getPointer(),
                                              tp_glb_.getPointer(), adr,
                                              epj_sorted_,   epj_neighbor_[id_thread],
                                              n_leaf_limit_, lev_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagSearchShortGather, const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        F64 r_search_sq = ptcl.getRSearch() * ptcl.getRSearch();
        const S32 adr = 0;
        SearchNeighborListOneParticleGather(pos_target,
                                            r_search_sq,
                                            tc_glb_.getPointer(),
                                            tp_glb_.getPointer(), adr, 
                                            epj_sorted_, epj_neighbor_[id_thread],
                                            n_leaf_limit_, lev_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }
#endif

    

    
    // 2016 02/05
    /*
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        SearchNeighborListOneParticleScatter(pos_target,  tc_glb_.getPointer(),
                                             tp_glb_.getPointer(), adr, 
                                             epj_sorted_, epj_neighbor_[id_thread],
                                             n_leaf_limit_, lev_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }
    */
    /*
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){

        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        SearchNeighborListOneParticleGather(pos_target,  tc_glb_.getPointer(),
                                            tp_glb_.getPointer(), adr, 
                                            epj_sorted_, epj_neighbor_[id_thread],
                                            n_leaf_limit_, lev_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;

    }
    */
    /*
    template<class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<SEARCH_MODE_SYMMETRY, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        epj_neighbor_[id_thread].clearSize();
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        SearchNeighborListOneParticleSymmetry(pos_target,  tc_glb_.getPointer(),
                                              tp_glb_.getPointer(), adr, 
                                              epj_sorted_, epj_neighbor_[id_thread],
                                              n_leaf_limit_, lev_leaf_limit_);
        S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }    
    */


    
#if 0
    // under construction
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneIPGroupImpl(TagSearchLongScatter, const Tptcl & ptcl, S32 & nnp){
        const F64vec pos_target = ptcl.getPos();
        const S32 adr = 0;
        const S32 size_old = epj_neighbor_.size();
        SearchNeighborListOneIPGroupScatter(pos_target,    tc_glb_.getPointer(),       adr, 
                                            epj_sorted_,   epj_neighbor_, 
                                            n_leaf_limit_, lev_leaf_limit_);
        nnp = epj_neighbor_.size() - size_old;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneIPGroup(const S32 iipg, S32 & nip, 
                              const Tepi * & epi, S32 & nnp, Tepj * & epj){
        nip = ipg_[iipg].n_ptcl_;
        const S32 adr = ipg_[iipg].adr_ptcl_;
        epi = epi_sorted_[adr].getPointer();
        const S32 head = epj_neighbor_.size();
        getNeighborListOneIPGroupImpl(typename TSM::search_type(), ptcl, nnp);
        epj = epj_neighbor_.getPointer(head);
    }
#endif

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    clearSizeOfArray(){
        tp_glb_.clearSize();
        tc_loc_.clearSize();
        tc_glb_.clearSize();
        epi_sorted_.clearSize();
        epi_org_.clearSize();
        epj_sorted_.clearSize();
        epj_org_.clearSize();
        spj_sorted_.clearSize();
        spj_org_.clearSize();
        ipg_.clearSize();
        epj_send_.clearSize();
        spj_send_.clearSize();
        force_sorted_.clearSize();
        force_org_.clearSize();
        epjr_sorted_.clearSize();
        epjr_send_.clearSize();
        epjr_recv_.clearSize();
        epjr_recv_1st_buf_.clearSize();
        epjr_recv_2nd_buf_.clearSize();
        const S32 n_thread = Comm::getNumberOfThread();
        for(S32 i=0; i<n_thread; i++){
            epj_for_force_[i].clearSize();
            spj_for_force_[i].clearSize();
            epjr_send_buf_[i].clearSize();
            epjr_send_buf_for_scatter_[i].clearSize();
            epjr_recv_1st_sorted_[i].clearSize();
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    F64ort TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getOuterBoundaryOfLocalTreeImpl(TagSearchLongSymmetry){
        return tc_loc_[0].mom_.vertex_out_;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    F64ort TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getInnerBoundaryOfLocalTreeImpl(TagSearchLongSymmetry){
        return tc_loc_[0].mom_.vertex_in_;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    F64ort TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getOuterBoundaryOfLocalTreeImpl(TagSearchShortSymmetry){
        return tc_loc_[0].mom_.vertex_out_;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    F64ort TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getInnerBoundaryOfLocalTreeImpl(TagSearchShortSymmetry){
        return tc_loc_[0].mom_.vertex_in_;
    }
    
    ////////////////////
    // for reuse
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseList(const DomainInfo & dinfo, const bool flag_reuse){
        if(typeid(TSM) == typeid(SEARCH_MODE_LONG)
           && dinfo.getBoundaryCondition() != BOUNDARY_CONDITION_OPEN){
            PARTICLE_SIMULATOR_PRINT_ERROR("The forces w/o cutoff can be evaluated only under the open boundary condition");
            Abort(-1);
        }
        if(!flag_reuse){ comm_table_.clearSize(); }
        //comm_table_.clear();
        exchangeLocalEssentialTreeReuseListImpl(typename TSM::search_type(), dinfo, flag_reuse);
        time_profile_.exchange_LET_tot = time_profile_.make_LET_1st
            + time_profile_.exchange_LET_1st
            + time_profile_.make_LET_2nd
            + time_profile_.exchange_LET_2nd;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListLong(const DomainInfo & dinfo,
                                            const bool flag_reuse){
        F64 time_offset = GetWtime();

        F64 r_crit_sq = LARGE_FLOAT;
        if (theta_ > 0) r_crit_sq = (length_ * length_) / (theta_ * theta_);
        else r_crit_sq = - 1.0; // negative value is used to represent theta_ = 0
        if(!flag_reuse){
            FindScatterParticle<TSM, TreeCell<Tmomloc>, TreeParticle,
                                Tepj, Tspj, WALK_MODE_NORMAL>
                (tc_loc_, tp_glb_,
                 epj_sorted_,
                 comm_table_.n_ep_send_,   comm_table_.adr_ep_send_, 
                 dinfo, n_leaf_limit_, lev_leaf_limit_, 
                 comm_table_.n_sp_send_,   comm_table_.adr_sp_send_, 
                 comm_table_.shift_per_image_,
                 comm_table_.n_image_per_proc_,
                 comm_table_.n_ep_per_image_,
                 comm_table_.n_sp_per_image_,
                 r_crit_sq,
                 icut_, level_pm_cell_);
            ExchangeNumber(comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                           comm_table_.n_sp_send_, comm_table_.n_sp_recv_);
            const S32 n_proc = Comm::getNumberOfProc();
            comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
            comm_table_.n_sp_send_tot_ = comm_table_.adr_sp_send_.size();
            comm_table_.n_ep_recv_tot_ = comm_table_.n_sp_recv_tot_ = 0;
            S32 n_ep_recv_tot_tmp = 0;
            S32 n_sp_recv_tot_tmp = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for reduction(+:n_ep_recv_tot_tmp), reduction(+:n_sp_recv_tot_tmp)
#endif
            for(S32 i=0; i<n_proc; i++){
                n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[i];
                n_sp_recv_tot_tmp += comm_table_.n_sp_recv_[i];
            }
            comm_table_.n_ep_recv_tot_ += n_ep_recv_tot_tmp;
            comm_table_.n_sp_recv_tot_ += n_sp_recv_tot_tmp;
        }

        time_profile_.make_LET_1st += GetWtime() - time_offset;
        time_offset = GetWtime();

        ExchangeLet<TSM, Tepj, Tspj>(epj_sorted_, comm_table_.n_ep_send_,
                                     comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_,
                                     comm_table_.adr_ep_send_,
                                     //epj_recv_,
                                     epj_org_, n_loc_tot_,
                                     tc_loc_, comm_table_.n_sp_send_,
                                     comm_table_.n_sp_recv_, comm_table_.n_sp_per_image_,
                                     comm_table_.adr_sp_send_,
                                     //spj_recv_,
                                     spj_org_,
                                     comm_table_.shift_per_image_,
                                     comm_table_.n_image_per_proc_);

        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchLong,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        exchangeLocalEssentialTreeReuseListLong(dinfo, flag_reuse);
    }    
    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchLongCutoff,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        exchangeLocalEssentialTreeReuseListLong(dinfo, flag_reuse);
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchLongScatter,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        exchangeLocalEssentialTreeReuseListLong(dinfo, flag_reuse);
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchLongSymmetry,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        exchangeLocalEssentialTreeReuseListLong(dinfo, flag_reuse);
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchLongParticleMeshMultipole,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        exchangeLocalEssentialTreeReuseListLong(dinfo, flag_reuse);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchShortScatter,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        F64 time_offset = GetWtime();

        if(!flag_reuse){
            FindScatterParticle<TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj, WALK_MODE_NORMAL>
                (tc_loc_, tp_glb_, epj_sorted_,
                 comm_table_.n_ep_send_,  comm_table_.adr_ep_send_,
                 dinfo, n_leaf_limit_, lev_leaf_limit_,
                 comm_table_.shift_per_image_,
                 comm_table_.n_image_per_proc_,
                 comm_table_.n_ep_per_image_, 0);

            ExchangeNumber(comm_table_.n_ep_send_, comm_table_.n_ep_recv_);

            const S32 n_proc = Comm::getNumberOfProc();
            comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
            comm_table_.n_ep_recv_tot_ = 0;
            S32 n_ep_recv_tot_tmp = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for reduction(+:n_ep_recv_tot_tmp)
#endif
            for(S32 i=0; i<n_proc; i++){
                n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[i];
            }
            comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
        }
        
        time_profile_.make_LET_1st += GetWtime() - time_offset;
        time_offset = GetWtime();

        ExchangeLet(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                    comm_table_.n_ep_per_image_,
                    comm_table_.adr_ep_send_,
                    //epj_recv_,
                    epj_org_, n_loc_tot_, 
                    comm_table_.shift_per_image_,
                    comm_table_.n_image_per_proc_);

        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    }

    ////////////////
    // SYMMETRY
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchShortSymmetry,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        F64 time_offset = GetWtime();

        const S32 n_proc = Comm::getNumberOfProc();
        static ReallocatableArray<S32> n_ep_send_per_proc_1st;
        static ReallocatableArray<S32> n_ep_send_per_proc_2nd;
        static ReallocatableArray<S32> n_ep_recv_per_proc_1st;
        static ReallocatableArray<S32> n_ep_recv_per_proc_2nd;
        static ReallocatableArray<S32> adr_ep_send_1st;
        static ReallocatableArray<S32> adr_ep_send_2nd;
        static ReallocatableArray<F64vec> shift_per_image_1st;
        static ReallocatableArray<F64vec> shift_per_image_2nd;
        static ReallocatableArray<S32> n_image_per_proc_1st;
        static ReallocatableArray<S32> n_image_per_proc_2nd;
        static ReallocatableArray<S32> n_ep_send_per_image_1st;
        static ReallocatableArray<S32> n_ep_send_per_image_2nd;
        static ReallocatableArray<Tepj> ep_recv_1st;
        
        if(!flag_reuse){
            ////////////
            // 1st STEP (send j particles)
            FindScatterParticle<TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj, WALK_MODE_NORMAL>
                (tc_loc_, tp_glb_, epj_sorted_,
                 n_ep_send_per_proc_1st,  adr_ep_send_1st,
                 dinfo, n_leaf_limit_, lev_leaf_limit_,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 n_ep_send_per_image_1st, 0);
            ExchangeNumber(n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st);
            ExchangeParticle(epj_sorted_, n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st,
                             n_ep_send_per_image_1st,
                             adr_ep_send_1st, ep_recv_1st,
                             shift_per_image_1st,
                             n_image_per_proc_1st);
            ////////////
            // 2nd STEP (send j particles)
            FindExchangeParticleDoubleWalk<TreeCell<Tmomloc>, TreeParticle, Tepj>
                (ep_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo,
                 n_leaf_limit_, lev_leaf_limit_,
                 n_ep_send_per_proc_2nd, n_ep_send_per_image_2nd,
                 n_image_per_proc_2nd,
                 adr_ep_send_2nd, shift_per_image_2nd,
                 epj_sorted_,
                 center_, length_);
            /////////////////////
            // 3rd STEP (exchange # of particles again)
            ExchangeNumber(n_ep_send_per_proc_2nd, n_ep_recv_per_proc_2nd);

            ReallocatableArray<S32> n_disp_ep_send_per_proc;
            n_disp_ep_send_per_proc.resizeNoInitialize(n_proc+1);
            n_disp_ep_send_per_proc[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_ep_send_per_proc[i+1] = n_disp_ep_send_per_proc[i]
                    + n_ep_send_per_proc_1st[i]
                    + n_ep_send_per_proc_2nd[i];
            }

            ReallocatableArray<S32> n_disp_image_per_proc;
            ReallocatableArray<S32> n_disp_image_per_proc_1st;
            ReallocatableArray<S32> n_disp_image_per_proc_2nd;
            n_disp_image_per_proc.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc_1st.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc_2nd.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc[0] = n_disp_image_per_proc_1st[0] = n_disp_image_per_proc_2nd[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_image_per_proc[i+1]     = n_disp_image_per_proc[i]     + n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
                n_disp_image_per_proc_1st[i+1] = n_disp_image_per_proc_1st[i] + n_image_per_proc_1st[i];
                n_disp_image_per_proc_2nd[i+1] = n_disp_image_per_proc_2nd[i] + n_image_per_proc_2nd[i];
            }

            const S32 n_image_tot_1st = shift_per_image_1st.size();
            const S32 n_image_tot_2nd = shift_per_image_2nd.size();
            ReallocatableArray<S32> n_disp_ep_send_per_image_1st;
            ReallocatableArray<S32> n_disp_ep_send_per_image_2nd;
            n_disp_ep_send_per_image_1st.resizeNoInitialize(n_image_tot_1st+1);
            n_disp_ep_send_per_image_2nd.resizeNoInitialize(n_image_tot_2nd+1);
            n_disp_ep_send_per_image_1st[0] = n_disp_ep_send_per_image_2nd[0] = 0;
            for(S32 i=0; i<n_image_tot_1st; i++){
                n_disp_ep_send_per_image_1st[i+1] = n_disp_ep_send_per_image_1st[i] + n_ep_send_per_image_1st[i];
            }
            for(S32 i=0; i<n_image_tot_2nd; i++){
                n_disp_ep_send_per_image_2nd[i+1] = n_disp_ep_send_per_image_2nd[i] + n_ep_send_per_image_2nd[i];
            }
            
            const S32 n_image_tot = n_disp_image_per_proc_1st[n_proc]  + n_disp_image_per_proc_2nd[n_proc];
            comm_table_.shift_per_image_.resizeNoInitialize(n_image_tot);
            comm_table_.n_ep_per_image_.resizeNoInitialize(n_image_tot);
            const S32 n_send_tot = n_disp_ep_send_per_proc[n_proc];
            comm_table_.adr_ep_send_.resizeNoInitialize(n_send_tot);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_proc; i++){
                S32 n_ep_cnt = 0;
                S32 n_image_cnt = 0;
                const S32 ep_head    = n_disp_ep_send_per_proc[i];
                const S32 image_head = n_disp_image_per_proc[i];
                const S32 image_head_1st = n_disp_image_per_proc_1st[i];
                const S32 image_end_1st  = n_disp_image_per_proc_1st[i+1];
                for(S32 j=image_head_1st; j<image_end_1st; j++, n_image_cnt++){
                    const S32 ep_head_1st = n_disp_ep_send_per_image_1st[j];
                    const S32 ep_end_1st  = n_disp_ep_send_per_image_1st[j+1];
                    comm_table_.shift_per_image_[image_head+n_image_cnt] = shift_per_image_1st[j];
                    comm_table_.n_ep_per_image_[image_head+n_image_cnt]  = n_ep_send_per_image_1st[j];
                    for(S32 k=ep_head_1st; k<ep_end_1st; k++, n_ep_cnt++){
                        comm_table_.adr_ep_send_[ep_head+n_ep_cnt] = adr_ep_send_1st[k];
                    }
                }
                const S32 image_head_2nd = n_disp_image_per_proc_2nd[i];
                const S32 image_end_2nd  = n_disp_image_per_proc_2nd[i+1];
                for(S32 j=image_head_2nd; j<image_end_2nd; j++, n_image_cnt++){
                    const S32 ep_head_2nd = n_disp_ep_send_per_image_2nd[j];
                    const S32 ep_end_2nd  = n_disp_ep_send_per_image_2nd[j+1];
                    comm_table_.shift_per_image_[image_head+n_image_cnt] = shift_per_image_2nd[j];
                    comm_table_.n_ep_per_image_[image_head+n_image_cnt]  = n_ep_send_per_image_2nd[j];
                    for(S32 k=ep_head_2nd; k<ep_end_2nd; k++, n_ep_cnt++){
                        comm_table_.adr_ep_send_[ep_head+n_ep_cnt] = adr_ep_send_2nd[k];
                    }
                }
            }
            comm_table_.n_ep_send_.resizeNoInitialize(n_proc);
            comm_table_.n_ep_recv_.resizeNoInitialize(n_proc);
            comm_table_.n_image_per_proc_.resizeNoInitialize(n_proc);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_proc; i++){
                comm_table_.n_ep_send_[i] = n_ep_send_per_proc_1st[i] + n_ep_send_per_proc_2nd[i];
                comm_table_.n_ep_recv_[i] = n_ep_recv_per_proc_1st[i] + n_ep_recv_per_proc_2nd[i];
                comm_table_.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            //for(S32 i=0; i<n_image_tot; i++){
            for(S32 i=0; i<n_proc; i++){
                comm_table_.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
            }
            comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
            S32 n_ep_recv_tot_tmp = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for reduction(+:n_ep_recv_tot_tmp)
#endif
            for(S32 i=0; i<n_proc; i++){
                n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[i];
            }
            comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
        } // end of reuse

        time_profile_.make_LET_1st += GetWtime() - time_offset;
        time_offset = GetWtime();

        ExchangeLet(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                    comm_table_.n_ep_per_image_,
                    comm_table_.adr_ep_send_,
                    //epj_recv_,
                    epj_org_, n_loc_tot_,
                    comm_table_.shift_per_image_, comm_table_.n_image_per_proc_);

        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    }


    ////////////////
    // GATHER MODE
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeReuseListImpl(TagSearchShortGather,
                                            const DomainInfo & dinfo,
                                            const bool flag_reuse){
        F64 time_offset = GetWtime();
    
        const S32 n_proc = Comm::getNumberOfProc();
        static ReallocatableArray<S32> n_ep_send_per_proc_1st;
        static ReallocatableArray<S32> n_ep_send_per_proc_2nd;
        static ReallocatableArray<S32> n_ep_recv_per_proc_1st;
        static ReallocatableArray<S32> n_ep_recv_per_proc_2nd;
        static ReallocatableArray<S32> adr_ep_send_1st;
        static ReallocatableArray<S32> adr_ep_send_2nd;
        static ReallocatableArray<F64vec> shift_per_image_1st;
        static ReallocatableArray<F64vec> shift_per_image_2nd;
        static ReallocatableArray<S32> n_image_per_proc_1st;
        static ReallocatableArray<S32> n_image_per_proc_2nd;
        static ReallocatableArray<S32> n_ep_send_per_image_1st;
        static ReallocatableArray<S32> n_ep_send_per_image_2nd;
        static ReallocatableArray<EssentialParticleBase> ep_recv_1st;
        static ReallocatableArray<EssentialParticleBase> epi_base_sorted;
        const S32 n_epi_sorted = epi_sorted_.size();
        epi_base_sorted.resizeNoInitialize(n_epi_sorted);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_epi_sorted; i++){
            epi_base_sorted[i].pos      = epi_sorted_[i].getPos();
            epi_base_sorted[i].r_search = epi_sorted_[i].getRSearch();
        }
        if(!flag_reuse){
            ////////////
            // 1st STEP (send epi_base)
            FindScatterParticle<TreeCell<Tmomloc>, TreeParticle, EssentialParticleBase, Tspj, WALK_MODE_NORMAL>
                (tc_loc_, tp_glb_, epi_base_sorted,
                 n_ep_send_per_proc_1st,  adr_ep_send_1st,
                 dinfo, n_leaf_limit_, lev_leaf_limit_,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 n_ep_send_per_image_1st, 0);

            //std::cerr<<"check b"<<std::endl;
            
            ExchangeNumber(n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st);

            //std::cerr<<"check c"<<std::endl;
            
            // exchange epi_base particles
            ExchangeParticle(epi_base_sorted, n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st,
                             n_ep_send_per_image_1st,
                             adr_ep_send_1st, ep_recv_1st,
                             shift_per_image_1st,
                             n_image_per_proc_1st);

            //std::cerr<<"check d"<<std::endl;
            
            ////////////
            // 2nd STEP (find j particle)
            FindExchangeParticleDoubleWalk<TreeCell<Tmomloc>, TreeParticle, EssentialParticleBase>
                (ep_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo,
                 n_leaf_limit_, lev_leaf_limit_, 
                 n_ep_send_per_proc_2nd, n_ep_send_per_image_2nd,
                 n_image_per_proc_2nd,
                 adr_ep_send_2nd, shift_per_image_2nd,
                 epi_base_sorted,
                 center_, length_);

            //std::cerr<<"check e"<<std::endl;
            
            /////////////////////
            // 3rd STEP (exchange # of particles again)
            ExchangeNumber(n_ep_send_per_proc_2nd, n_ep_recv_per_proc_2nd);

            //std::cerr<<"check f"<<std::endl;
            
            ReallocatableArray<S32> n_disp_ep_send_per_proc;
            n_disp_ep_send_per_proc.resizeNoInitialize(n_proc+1);
            n_disp_ep_send_per_proc[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_ep_send_per_proc[i+1] = n_disp_ep_send_per_proc[i]
                    + n_ep_send_per_proc_1st[i]
                    + n_ep_send_per_proc_2nd[i];
            }

            //std::cerr<<"check g"<<std::endl;
            
            ReallocatableArray<S32> n_disp_image_per_proc;
            ReallocatableArray<S32> n_disp_image_per_proc_1st;
            ReallocatableArray<S32> n_disp_image_per_proc_2nd;
            n_disp_image_per_proc.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc_1st.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc_2nd.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc[0] = n_disp_image_per_proc_1st[0] = n_disp_image_per_proc_2nd[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_image_per_proc[i+1]     = n_disp_image_per_proc[i]     + n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
                n_disp_image_per_proc_1st[i+1] = n_disp_image_per_proc_1st[i] + n_image_per_proc_1st[i];
                n_disp_image_per_proc_2nd[i+1] = n_disp_image_per_proc_2nd[i] + n_image_per_proc_2nd[i];
            }

            //std::cerr<<"check h"<<std::endl;
            
            const S32 n_image_tot_1st = shift_per_image_1st.size();
            const S32 n_image_tot_2nd = shift_per_image_2nd.size();
            ReallocatableArray<S32> n_disp_ep_send_per_image_1st;
            ReallocatableArray<S32> n_disp_ep_send_per_image_2nd;
            n_disp_ep_send_per_image_1st.resizeNoInitialize(n_image_tot_1st+1);
            n_disp_ep_send_per_image_2nd.resizeNoInitialize(n_image_tot_2nd+1);
            n_disp_ep_send_per_image_1st[0] = n_disp_ep_send_per_image_2nd[0] = 0;
            for(S32 i=0; i<n_image_tot_1st; i++){
                n_disp_ep_send_per_image_1st[i+1] = n_disp_ep_send_per_image_1st[i] + n_ep_send_per_image_1st[i];
            }
            for(S32 i=0; i<n_image_tot_2nd; i++){
                n_disp_ep_send_per_image_2nd[i+1] = n_disp_ep_send_per_image_2nd[i] + n_ep_send_per_image_2nd[i];
            }

            //std::cerr<<"check i"<<std::endl;
            
            const S32 n_image_tot = n_disp_image_per_proc_1st[n_proc]  + n_disp_image_per_proc_2nd[n_proc];
            comm_table_.shift_per_image_.resizeNoInitialize(n_image_tot);
            comm_table_.n_ep_per_image_.resizeNoInitialize(n_image_tot);
            const S32 n_send_tot = n_disp_ep_send_per_proc[n_proc];
            comm_table_.adr_ep_send_.resizeNoInitialize(n_send_tot);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_proc; i++){
                S32 n_ep_cnt = 0;
                S32 n_image_cnt = 0;
                const S32 ep_head    = n_disp_ep_send_per_proc[i];
                const S32 image_head = n_disp_image_per_proc[i];
                const S32 image_head_1st = n_disp_image_per_proc_1st[i];
                const S32 image_end_1st  = n_disp_image_per_proc_1st[i+1];
                for(S32 j=image_head_1st; j<image_end_1st; j++, n_image_cnt++){
                    const S32 ep_head_1st = n_disp_ep_send_per_image_1st[j];
                    const S32 ep_end_1st  = n_disp_ep_send_per_image_1st[j+1];
                    comm_table_.shift_per_image_[image_head+n_image_cnt] = shift_per_image_1st[j];
                    comm_table_.n_ep_per_image_[image_head+n_image_cnt]  = n_ep_send_per_image_1st[j];
                    for(S32 k=ep_head_1st; k<ep_end_1st; k++, n_ep_cnt++){
                        comm_table_.adr_ep_send_[ep_head+n_ep_cnt] = adr_ep_send_1st[k];
                    }
                }
                const S32 image_head_2nd = n_disp_image_per_proc_2nd[i];
                const S32 image_end_2nd  = n_disp_image_per_proc_2nd[i+1];
                for(S32 j=image_head_2nd; j<image_end_2nd; j++, n_image_cnt++){
                    const S32 ep_head_2nd = n_disp_ep_send_per_image_2nd[j];
                    const S32 ep_end_2nd  = n_disp_ep_send_per_image_2nd[j+1];
                    comm_table_.shift_per_image_[image_head+n_image_cnt] = shift_per_image_2nd[j];
                    comm_table_.n_ep_per_image_[image_head+n_image_cnt]  = n_ep_send_per_image_2nd[j];
                    for(S32 k=ep_head_2nd; k<ep_end_2nd; k++, n_ep_cnt++){
                        comm_table_.adr_ep_send_[ep_head+n_ep_cnt] = adr_ep_send_2nd[k];
                    }
                }
            }

            //std::cerr<<"check j"<<std::endl;
            
            comm_table_.n_ep_send_.resizeNoInitialize(n_proc);
            comm_table_.n_ep_recv_.resizeNoInitialize(n_proc);
            comm_table_.n_image_per_proc_.resizeNoInitialize(n_proc);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_proc; i++){
                comm_table_.n_ep_send_[i] = n_ep_send_per_proc_1st[i] + n_ep_send_per_proc_2nd[i];
                comm_table_.n_ep_recv_[i] = n_ep_recv_per_proc_1st[i] + n_ep_recv_per_proc_2nd[i];
                comm_table_.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
            }

            //std::cerr<<"check k"<<std::endl;


#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            //for(S32 i=0; i<n_image_tot; i++){ // org, NEED TO CHECK
            for(S32 i=0; i<n_proc; i++){
                comm_table_.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
            }
            comm_table_.n_ep_send_tot_ = comm_table_.adr_ep_send_.size();
            S32 n_ep_recv_tot_tmp = 0;


            
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for reduction(+:n_ep_recv_tot_tmp)
#endif
            for(S32 i=0; i<n_proc; i++){
                n_ep_recv_tot_tmp += comm_table_.n_ep_recv_[i];
            }

            comm_table_.n_ep_recv_tot_ = n_ep_recv_tot_tmp;
            //std::cerr<<"check m"<<std::endl;
        } // end of reuse flag
        //std::cerr<<"check n"<<std::endl;
        time_profile_.make_LET_1st += GetWtime() - time_offset;
        time_offset = GetWtime();

        ExchangeLet(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                    comm_table_.n_ep_per_image_,
                    comm_table_.adr_ep_send_,
                    //epj_recv_,
                    epj_org_, n_loc_tot_,
                    comm_table_.shift_per_image_, comm_table_.n_image_per_proc_);

        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
        //std::cerr<<"check z"<<std::endl;
    }
    
    template<class TSM>
    F64ort GetIpgBoxForInteractionList(const IPGroup<TSM> & ipg){
        return ipg.vertex_;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListIndexShort(){
        static bool first = true;
        static ReallocatableArray<S32> * adr_epj_tmp;
        static ReallocatableArray<S32> * adr_ipg_tmp;
        static ReallocatableArray<S32> * n_disp_epj_tmp;
        const S32 n_thread = Comm::getNumberOfThread();
        if(first){
            adr_epj_tmp = new ReallocatableArray<S32>[n_thread];
            adr_ipg_tmp = new ReallocatableArray<S32>[n_thread];
            n_disp_epj_tmp = new ReallocatableArray<S32>[n_thread];
            first = false;
        }
        const S32 n_ipg = ipg_.size();
        const S32 adr_tc = 0;
        ReallocatableArray<Tspj> spj_dummy;
        ReallocatableArray<S32> adr_spj;
        interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
        interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_ep_cum_prev = 0;
            adr_epj_tmp[ith].clearSize();
            adr_ipg_tmp[ith].clearSize();
            n_disp_epj_tmp[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 i=0; i<n_ipg; i++){
                adr_ipg_tmp[ith].push_back(i);
                n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
                const F64 r_crit_sq = 999.9;
                TargetBox<TSM> target_box;
                target_box.set(ipg_[i]);
                MakeListUsingTreeRecursiveTop
                    <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj,
                     WALK_MODE_NORMAL, TagChopLeafTrue, TagChopNonleafFalse>
                    (tc_glb_,  adr_tc, tp_glb_,
                     epj_sorted_, adr_epj_tmp[ith],
                     target_box,
                     r_crit_sq,
                     n_leaf_limit_,
                     lev_leaf_limit_,
                     F64vec(0.0), 0);

                interaction_list_.n_ep_[i] = adr_epj_tmp[ith].size() - n_ep_cum_prev;
                n_ep_cum_prev = adr_epj_tmp[ith].size();
            }
            n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
        }
        interaction_list_.n_disp_ep_[0] = 0;
        for(S32 i=0; i<n_ipg; i++){
            interaction_list_.n_disp_ep_[i+1] = interaction_list_.n_disp_ep_[i] + interaction_list_.n_ep_[i];
        }
        interaction_list_.adr_ep_.resizeNoInitialize( interaction_list_.n_disp_ep_[n_ipg] );

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_thread; i++){
            for(S32 j=0; j<adr_ipg_tmp[i].size(); j++){
                const S32 adr_ipg = adr_ipg_tmp[i][j];
                S32 adr_ep = interaction_list_.n_disp_ep_[adr_ipg];
                const S32 k_h = n_disp_epj_tmp[i][j];
                const S32 k_e = n_disp_epj_tmp[i][j+1];
                for(S32 k=k_h; k<k_e; k++, adr_ep++){
                    interaction_list_.adr_ep_[adr_ep] = adr_epj_tmp[i][k];
                }
            }
        }
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListIndexLong(){
        static bool first = true;
        static ReallocatableArray<S32> * adr_epj_tmp;
        static ReallocatableArray<S32> * adr_spj_tmp;
        static ReallocatableArray<S32> * adr_ipg_tmp;
        static ReallocatableArray<S32> * n_disp_epj_tmp;
        static ReallocatableArray<S32> * n_disp_spj_tmp;
        const S32 n_thread = Comm::getNumberOfThread();
        if(first){
            adr_epj_tmp = new ReallocatableArray<S32>[n_thread];
            adr_spj_tmp = new ReallocatableArray<S32>[n_thread];
            adr_ipg_tmp = new ReallocatableArray<S32>[n_thread];
            n_disp_epj_tmp = new ReallocatableArray<S32>[n_thread];
            n_disp_spj_tmp = new ReallocatableArray<S32>[n_thread];
            first = false;
        }
        const S32 n_ipg = ipg_.size();
        const S32 adr_tc = 0;
        const S32 adr_tree_sp_first = spj_sorted_.size() - tc_glb_.size();
        F64 r_crit_sq = LARGE_FLOAT;
        if (theta_ > 0.0) r_crit_sq = (length_ * length_) / (theta_ * theta_);
        else r_crit_sq = - 1.0; // negative value is used to represent theta_ = 0
        interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
        interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
        interaction_list_.n_sp_.resizeNoInitialize(n_ipg);
        interaction_list_.n_disp_sp_.resizeNoInitialize(n_ipg+1);
        
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_ep_cum_prev = 0;
            S32 n_sp_cum_prev = 0;
            adr_epj_tmp[ith].clearSize();
            adr_spj_tmp[ith].clearSize();
            adr_ipg_tmp[ith].clearSize();
            n_disp_epj_tmp[ith].clearSize();
            n_disp_spj_tmp[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 i=0; i<n_ipg; i++){
                //F64 m_tmp = 0.0;
                adr_ipg_tmp[ith].push_back(i);
                n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
                n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
                TargetBox<TSM> target_box;
                target_box.set(ipg_[i], icut_);
                if (typeid(TSM) != typeid(SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE)) {
                    MakeListUsingTreeRecursiveTop
                        <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, Tspj,
                         WALK_MODE_NORMAL, TagChopLeafFalse, TagCopyInfoCloseWithTpAdrptcl,
                         TagChopNonleafFalse>
                        (tc_glb_,  adr_tc, tp_glb_,
                         epj_sorted_, adr_epj_tmp[ith],
                         spj_sorted_, adr_spj_tmp[ith],
                         target_box,
                         r_crit_sq, n_leaf_limit_, lev_leaf_limit_,
                         adr_tree_sp_first, F64vec(0.0), 0);
                } else {
                    MakeListUsingTreeRecursiveTop
                        <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, Tspj,
                         WALK_MODE_NORMAL, TagChopLeafTrue, TagCopyInfoCloseWithTpAdrptcl,
                         TagChopNonleafTrue>
                        (tc_glb_,  adr_tc, tp_glb_,
                         epj_sorted_, adr_epj_tmp[ith],
                         spj_sorted_, adr_spj_tmp[ith],
                         target_box,
                         r_crit_sq, n_leaf_limit_, lev_leaf_limit_,
                         adr_tree_sp_first, F64vec(0.0), level_pm_cell_);
                }
                interaction_list_.n_ep_[i] = adr_epj_tmp[ith].size() - n_ep_cum_prev;
                interaction_list_.n_sp_[i] = adr_spj_tmp[ith].size() - n_sp_cum_prev;
                n_ep_cum_prev = adr_epj_tmp[ith].size();
                n_sp_cum_prev = adr_spj_tmp[ith].size();
            }
            n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
            n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
        } // end of OMP

        interaction_list_.n_disp_ep_[0] = 0;
        interaction_list_.n_disp_sp_[0] = 0;
        for(S32 i=0; i<n_ipg; i++){
            interaction_list_.n_disp_ep_[i+1] = interaction_list_.n_disp_ep_[i] + interaction_list_.n_ep_[i];
            interaction_list_.n_disp_sp_[i+1] = interaction_list_.n_disp_sp_[i] + interaction_list_.n_sp_[i];
        }
        interaction_list_.adr_ep_.resizeNoInitialize( interaction_list_.n_disp_ep_[n_ipg] );
        interaction_list_.adr_sp_.resizeNoInitialize( interaction_list_.n_disp_sp_[n_ipg] );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_thread; i++){
            for(S32 j=0; j<adr_ipg_tmp[i].size(); j++){
                const S32 adr_ipg = adr_ipg_tmp[i][j];
                S32 adr_ep = interaction_list_.n_disp_ep_[adr_ipg];
                const S32 k_ep_h = n_disp_epj_tmp[i][j];
                const S32 k_ep_e = n_disp_epj_tmp[i][j+1];
                for(S32 k=k_ep_h; k<k_ep_e; k++, adr_ep++){
                    interaction_list_.adr_ep_[adr_ep] = adr_epj_tmp[i][k];
                }
                S32 adr_sp = interaction_list_.n_disp_sp_[adr_ipg];
                const S32 k_sp_h = n_disp_spj_tmp[i][j];
                const S32 k_sp_e = n_disp_spj_tmp[i][j+1];
                for(S32 k=k_sp_h; k<k_sp_e; k++, adr_sp++){
                    interaction_list_.adr_sp_[adr_sp] = adr_spj_tmp[i][k];
                }
            }
        }
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::    
    freeMem(){
        tp_glb_.freeMem();
        tc_loc_.freeMem();
        tc_glb_.freeMem();
        epi_sorted_.freeMem();
        epi_org_.freeMem();
        epj_sorted_.freeMem();
        epj_org_.freeMem();
        spj_sorted_.freeMem();
        spj_org_.freeMem();
        ipg_.freeMem();
        epj_send_.freeMem();
        spj_send_.freeMem();
        force_sorted_.freeMem();
        force_org_.freeMem();
        epjr_sorted_.freeMem();
        epjr_send_.freeMem();
        epjr_recv_.freeMem();
        epjr_recv_1st_buf_.freeMem();
        epjr_recv_2nd_buf_.freeMem();
        const S32 n_thread = Comm::getNumberOfThread();
        for(S32 i=0; i<n_thread; i++){
            epj_for_force_[i].freeMem();
            spj_for_force_[i].freeMem();
            epjr_send_buf_[i].freeMem();
            epjr_send_buf_for_scatter_[i].freeMem();
            epjr_recv_1st_sorted_[i].freeMem();
        }
    }
    
}

#include"tree_for_force_impl_force.hpp"
