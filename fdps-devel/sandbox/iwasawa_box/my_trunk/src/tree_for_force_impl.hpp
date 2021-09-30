////////////////////////////////////////////////
/// implementaion of methods of TreeForForce ///

#include"tree_walk.hpp"

namespace ParticleSimulator{
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    Tepj * TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getEpjFromId(const S64 id, const Tepj * epj_tmp){
#pragma omp critical
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
            + tp_buf_.getMemSize() + tp_loc_.getMemSize() + tp_glb_.getMemSize()
            + tc_loc_.getMemSize() + tc_glb_.getMemSize()
            + epi_sorted_.getMemSize() + epi_org_.getMemSize()
            + epj_sorted_.getMemSize() + epj_org_.getMemSize()
            + spj_sorted_.getMemSize() + spj_org_.getMemSize()
            + ipg_.getMemSize()
            + epj_send_.getMemSize() + epj_recv_.getMemSize()
            + spj_send_.getMemSize() + spj_recv_.getMemSize()
            + force_.getMemSize() + force_buf_.getMemSize()
            + epjr_sorted_.getMemSize() + epjr_send_.getMemSize()
            + epjr_recv_.getMemSize() + epjr_recv_1st_buf_.getMemSize()
            + epjr_recv_2nd_buf_.getMemSize();
    }
    

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    initialize(const U64 n_glb_tot,
               const F64 theta,
               const U32 n_leaf_limit,
               const U32 n_group_limit){
        if(is_initialized_ == true){
            PARTICLE_SIMULATOR_PRINT_ERROR("Do not initialize the tree twice");
            std::cerr<<"SEARCH_MODE: "<<typeid(TSM).name()<<std::endl;
            std::cerr<<"Force: "<<typeid(Tforce).name()<<std::endl;
            std::cerr<<"EPI: "<<typeid(Tepi).name()<<std::endl;
            std::cerr<<"EPJ: "<<typeid(Tepj).name()<<std::endl;
            std::cerr<<"SPJ: "<<typeid(Tspj).name()<<std::endl;
            Abort(-1);
        }
        comm_table_.initialize();
        map_id_to_epj_.clear();
        is_initialized_ = true;
        n_glb_tot_ = n_glb_tot;
        theta_ = theta;
        n_leaf_limit_ = n_leaf_limit;
        n_group_limit_ = n_group_limit;
        lev_max_loc_ = lev_max_glb_ = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        const S64 n_proc = Comm::getNumberOfProc();
        n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = n_walk_local_ = 0;
        //n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
        ni_ave_ = nj_ave_ = 0;
        bool err = false;
        if(n_leaf_limit_ <= 0){
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The limit number of the particles in the leaf cell must be > 0");
            std::cout<<"n_leaf_limit_= "<<n_leaf_limit_<<std::endl;
        }
        if(n_group_limit_ < n_leaf_limit_){
            err = true;
            PARTICLE_SIMULATOR_PRINT_ERROR("The limit number of particles in ip graoups msut be >= that in leaf cells");
            std::cout<<"n_group_limit_= "<<n_group_limit_<<std::endl;
            std::cout<<"n_leaf_limit_= "<<n_leaf_limit_<<std::endl;
        }
        if( typeid(TSM) == typeid(SEARCH_MODE_LONG) || 
            typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF) ){
            if(theta_ < 0.0){
                err = true;
                PARTICLE_SIMULATOR_PRINT_ERROR("The opening criterion of the tree must be >= 0.0");
                std::cout<<"theta_= "<<theta_<<std::endl;
            }
        }
        if(err){
            std::cout<<"SEARCH_MODE: "<<typeid(TSM).name()<<std::endl;
            std::cout<<"Force: "<<typeid(Tforce).name()<<std::endl;
            std::cout<<"EPI: "<<typeid(Tepi).name()<<std::endl;
            std::cout<<"SPJ: "<<typeid(Tspj).name()<<std::endl;
            ParticleSimulator::Abort(-1);
        }
        epj_for_force_    = new ReallocatableArray<Tepj>[n_thread];
        spj_for_force_    = new ReallocatableArray<Tspj>[n_thread];
        //id_epj_for_force_ = new ReallocatableArray<S32>[n_thread];
        //id_spj_for_force_ = new ReallocatableArray<S32>[n_thread];

        epjr_send_buf_ = new ReallocatableArray<EPJWithR>[n_thread];
        epjr_send_buf_for_scatter_ = new ReallocatableArray<EPJWithR>[n_thread];
        epjr_recv_1st_sorted_ = new ReallocatableArray<EPJWithR>[n_thread];
        epj_neighbor_ = new ReallocatableArray<Tepj>[n_thread];
        n_cell_open_ = new CountT[n_thread];

        Comm::barrier();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        if(Comm::getRank() == 0){
            std::cerr<<"used mem size for tree(0)="<<this->getMemSizeUsed()*1e-9<<"[GB]"<<std::endl;
        }
#endif	
        if( typeid(TSM) == typeid(SEARCH_MODE_LONG) || 
            typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF) || 
            typeid(TSM) == typeid(SEARCH_MODE_LONG_SCATTER) ||
            typeid(TSM) == typeid(SEARCH_MODE_LONG_SYMMETRY)){
            n_sp_send_disp_ = new S32[n_proc+1];
            n_sp_recv_disp_ = new S32[n_proc+1];
            n_ep_sp_send_ = new S32[n_proc * 2];
            n_ep_sp_recv_ = new S32[n_proc * 2];
        }
        else{
            n_ep_sp_send_ = n_ep_sp_recv_ = NULL;
        }
        n_ep_send_disp_ = new S32[n_proc+1];
        n_ep_recv_disp_ = new S32[n_proc+1];
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        req_send_ = new MPI_Request[n_proc];
        req_recv_ = new MPI_Request[n_proc];
	status_   = new MPI_Status[n_proc];
#endif
    }


    //////////////////////////////////////
    // FUNCTIONS CALLED IN calcForceXXX //
    ///////////////
    // SET PARTICLE FOR LOCAL TREE
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tpsys>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setParticleLocalTree(const Tpsys & psys,
                         const bool clear){
        const F64 time_offset = GetWtime();
        const S32 nloc = psys.getNumberOfParticleLocal();
        if(clear){ n_loc_tot_ = 0;}
        const S32 offset = n_loc_tot_;
        n_loc_tot_ += nloc;
        epi_org_.resizeNoInitialize(n_loc_tot_);
        epj_org_.resizeNoInitialize(n_loc_tot_);
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

    ////////////////////////////
    // SET LET TO GLOBAL TREE //
#if 1
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setLocalEssentialTreeToGlobalTree(const bool flag_reuse){
        F64 time_offset = GetWtime();
        epj_org_.swap(epj_sorted_); // epj_sorted_ -> epj_org_
        //epj_sorted_.freeMem();
        this->n_glb_tot_ = n_loc_tot_ + epj_recv_.size() + spj_recv_.size();
        setLocalEssentialTreeToGlobalTreeImpl(typename TSM::force_type(), flag_reuse);
        //this->n_glb_tot_ = tp_glb_.size();
        time_profile_.set_particle_global_tree += GetWtime() - time_offset;
    }
#else
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setLocalEssentialTreeToGlobalTree(const bool flag_reuse){
        F64 time_offset = GetWtime();
        epj_org_.swap(epj_sorted_); // epj_sorted_ -> epj_org_
        //epj_sorted_.freeMem();
        if(!flag_reuse){
            tp_glb_.swap(tp_loc_);
            for(S32 i=0; i<n_loc_tot_; i++) tp_glb_[i].adr_ptcl_ = i;
            //tp_loc_.freeMem();
        }
        setLocalEssentialTreeToGlobalTreeImpl(typename TSM::force_type(), flag_reuse);
        this->n_glb_tot_ = tp_glb_.size();
        time_profile_.set_particle_global_tree += GetWtime() - time_offset;
    }
#endif

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setLocalEssentialTreeToGlobalTreeImpl(TagForceLong,
                                           const bool flag_reuse){
        SetLocalEssentialTreeToGlobalTreeImpl(epj_recv_, spj_recv_, n_loc_tot_,
                                              epj_org_,  spj_org_, tp_glb_,
                                              flag_reuse);
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setLocalEssentialTreeToGlobalTreeImpl(TagForceShort,
                                           const bool flag_reuse){
        SetLocalEssentialTreeToGlobalTreeImpl(epj_recv_, n_loc_tot_,
                                              epj_org_,  tp_glb_,
                                              flag_reuse);
    }
    
    /////////////////////////////
    // SET ROOT CELL OF LOCAL TREE
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setRootCell(const DomainInfo & dinfo){
        const F64 time_offset = GetWtime();
        if(dinfo.getBoundaryCondition() == BOUNDARY_CONDITION_OPEN){
            calcCenterAndLengthOfRootCellOpenImpl(typename TSM::search_type());
        }
        else{
            calcCenterAndLengthOfRootCellPeriodicImpl(typename TSM::search_type());
        }
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

    /////////////////////////////////////////////////////////
    // SORT PARTICLES IN LOCAL TREE (tp_) IN MORTON ORDER
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    mortonSortLocalTreeOnly(const bool reuse,
                            const bool keep_fp_order){
        const F64 wtime_offset = GetWtime();
        F64 wtime_offset_in = 0.0;
        if( reuse && !keep_fp_order ){
            epi_sorted_.swap(epi_org_); // epi_org_ -> epi_sorted_
            //epi_org_.freeMem();
#if 1
            epj_sorted_.swap(epj_org_); // epi_org_ -> epi_sorted_
#else
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                const Tepj epj = epj_org_[i];
                epj_sorted_[i] = epj;
            }
#endif
            time_profile_.morton_sort_local_tree += GetWtime() - wtime_offset;
            time_profile_.make_local_tree += GetWtime() - wtime_offset;
            return;
        }
        epi_sorted_.resizeNoInitialize(n_loc_tot_);
        epj_sorted_.resizeNoInitialize(n_loc_tot_);
        adr_org_from_adr_sorted_.resizeNoInitialize(n_loc_tot_);
        if(!reuse){
            tp_loc_.resizeNoInitialize(n_loc_tot_);
            tp_buf_.resizeNoInitialize(n_loc_tot_);
            wtime_offset_in = GetWtime();
            MortonKey::initialize( length_ * 0.5, center_);
            //wtime_offset_in = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                tp_loc_[i].setFromEP(epj_org_[i], i);
            }
            time_profile_.morton_key_local_tree += GetWtime() - wtime_offset_in;
            wtime_offset_in = GetWtime();
#ifdef USE_STD_SORT
            std::sort(tp_loc_.getPointer(), tp_loc_.getPointer()+n_loc_tot_,
                      [](const TreeParticle & l, const TreeParticle & r )
                      ->bool{return l.getKey() < r.getKey();} );
#else
            rs_.lsdSort(tp_loc_.getPointer(), tp_buf_.getPointer(), 0, n_loc_tot_-1);
#endif

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                const S32 adr = tp_loc_[i].adr_ptcl_;
                adr_org_from_adr_sorted_[i] = adr;
            }
            time_profile_.morton_sort_local_tree += GetWtime() - wtime_offset_in;
        }
        wtime_offset_in = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_loc_tot_; i++){
            const S32 adr = adr_org_from_adr_sorted_[i];
            epi_sorted_[i] = epi_org_[adr];
            epj_sorted_[i] = epj_org_[adr];
#ifdef SHRINK_SIZE_OF_SPJ_SORTED
            tp_loc_[i].adr_ptcl_ = i;
#endif
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"tp_loc_.size()="<<tp_loc_.size()<<" tp_buf_.size()="<<tp_buf_.size()<<std::endl;
        std::cout<<"epi_sorted_.size()="<<epi_sorted_.size()<<" epj_sorted_.size()="<<epj_sorted_.size()<<std::endl;
#endif
        time_profile_.morton_sort_local_tree__reorder += GetWtime() - wtime_offset_in;
        time_profile_.morton_sort_local_tree += GetWtime() - wtime_offset_in;
        //time_profile_.morton_sort_local_tree += GetWtime() - wtime_offset;
        time_profile_.make_local_tree += GetWtime() - wtime_offset;
    }


#if 1
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    mortonSortGlobalTreeOnly(const bool reuse){
        F64 time_offset = GetWtime();
        F64 wtime_offset_in = 0.0;
        if(!map_id_to_epj_.empty()){
            map_id_to_epj_.clear();
        }
        assert(map_id_to_epj_.empty());
        const S32 n_ep_add = epj_recv_.size();
        const S32 n_sp_add = spj_recv_.size();
        const S32 n_ep_tot = n_loc_tot_ + n_ep_add;
#ifdef SHRINK_SIZE_OF_SPJ_SORTED
        epj_sorted_.resizeNoInitialize(n_ep_tot);
#else
        epj_sorted_.resizeNoInitialize(n_glb_tot_);
#endif
        if(!reuse){
            tp_glb_.resizeNoInitialize( n_glb_tot_ );
            tp_glb_.swap(tp_loc_);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                //tp_glb_[i].key_ = tp_loc_[i].key_;
                tp_glb_[i].adr_ptcl_ = i;
            }

            wtime_offset_in = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
            for(S32 i=0; i<n_ep_add; i++){
                const S32 id_src = n_loc_tot_ + i;
                tp_glb_[id_src].setFromEP(epj_org_[id_src], id_src);
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif	    
            for(S32 i=0; i<n_sp_add; i++){
                const S32 id_dst = n_ep_tot + i;
#ifdef SHRINK_SIZE_OF_SPJ_SORTED
                tp_glb_[id_dst].setFromSP(spj_recv_[i], i);
#else
                tp_glb_[id_dst].setFromSP(spj_recv_[i], id_dst);
#endif
            }
            time_profile_.morton_key_global_tree += GetWtime() - wtime_offset_in;
            wtime_offset_in = GetWtime();
            tp_buf_.resizeNoInitialize(n_glb_tot_);
#ifdef USE_STD_SORT
            std::sort(tp_glb_.getPointer(), tp_glb_.getPointer()+n_glb_tot_, 
                      [](const TreeParticle & l, const TreeParticle & r )
                      ->bool{return l.getKey() < r.getKey();} );
#else
            rs_.lsdSort(tp_glb_.getPointer(), tp_buf_.getPointer(), 0, n_glb_tot_-1);
#endif
            time_profile_.morton_sort_global_tree += GetWtime() - wtime_offset_in;
        }
        wtime_offset_in = GetWtime();
        if( typeid(typename TSM::force_type) == typeid(TagForceLong)){
#ifdef SHRINK_SIZE_OF_SPJ_SORTED
            const S32 n_sp_tot = spj_recv_.size();
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
                    const U32 adr = tp_glb_[i].adr_ptcl_;
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
                    const U32 adr = tp_glb_[i].adr_ptcl_;
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
            static ReallocatableArray<U32> adr_ep;
            static ReallocatableArray<U32> adr_sp;
            adr_ep.resizeNoInitialize(n_ep_tot);
            adr_sp.resizeNoInitialize(n_sp_add);
            if(!reuse){
                U32 n_cnt_ep = 0;
                U32 n_cnt_sp = 0;
                for(S32 i=0; i<n_glb_tot_; i++){
                    const U32 adr = tp_glb_[i].adr_ptcl_;
                    if( GetMSB(adr) == 0){
                        epj_sorted_[n_cnt_ep] = epj_org_[adr];
                        tp_glb_[i].adr_ptcl_ = n_cnt_ep;
                        adr_ep[n_cnt_ep] = adr;
                        n_cnt_ep++;
                    }
                    else{
                        spj_sorted_[n_cnt_sp] = spj_org_[ClearMSB(adr)];
                        tp_glb_[i].adr_ptcl_ = SetMSB(n_cnt_sp);
                        adr_sp[n_cnt_sp] = ClearMSB(adr);
                        n_cnt_sp++;
                    }
                }
            }
            else{
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                for(S32 i=0; i<n_ep_tot; i++){
                    S32 adr = adr_ep[i];
                    epj_sorted_[i] = epj_org_[adr];
                }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                for(S32 i=0; i<n_sp_add; i++){
                    S32 adr = adr_sp[i];
                    spj_sorted_[i] = spj_org_[adr];
                }
            }
#endif
#else // SHRINK_SIZE_OF_SPJ_SORTED
            spj_sorted_.resizeNoInitialize( spj_org_.size() );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_glb_tot_; i++){
                const U32 adr = tp_glb_[i].adr_ptcl_;
                if( GetMSB(adr) == 0){
                    epj_sorted_[i] = epj_org_[adr];
                }
                else{
                    spj_sorted_[i] = spj_org_[ClearMSB(adr)];
                }
            }
#endif
        }
        else{
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_glb_tot_; i++){
                const U32 adr = tp_glb_[i].adr_ptcl_;
                epj_sorted_[i] = epj_org_[adr];
            }
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"tp_glb_.size()="<<tp_glb_.size()<<" tp_buf_.size()="<<tp_buf_.size()<<std::endl;
        std::cout<<"epj_sorted_.size()="<<epj_sorted_.size()<<" spj_sorted_.size()="<<spj_sorted_.size()<<std::endl;
#endif
        time_profile_.morton_sort_global_tree__reorder += GetWtime() - wtime_offset_in;
        time_profile_.morton_sort_global_tree += GetWtime() - wtime_offset_in;
        //time_profile_.morton_sort_global_tree += GetWtime() - time_offset;
        time_profile_.make_global_tree += GetWtime() - time_offset;
    }
#else
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    mortonSortGlobalTreeOnly(const bool reuse){
        F64 time_offset = GetWtime();
        if(!map_id_to_epj_.empty()){
            map_id_to_epj_.clear();
        }
        assert(map_id_to_epj_.empty());
#ifdef SHRINK_SIZE_OF_SPJ_SORTED
        const S32 n_ep_tot = n_loc_tot_ + epj_recv_.size();
        epj_sorted_.resizeNoInitialize(n_ep_tot);
#else
        epj_sorted_.resizeNoInitialize(n_glb_tot_);
#endif
        if(!reuse){
            tp_buf_.resizeNoInitialize(n_glb_tot_);
#ifdef USE_STD_SORT
            std::sort(tp_glb_.getPointer(), tp_glb_.getPointer()+n_glb_tot_, 
                      [](const TreeParticle & l, const TreeParticle & r )
                      ->bool{return l.getKey() < r.getKey();} );
#else
            rs_.lsdSort(tp_glb_.getPointer(), tp_buf_.getPointer(), 0, n_glb_tot_-1);
#endif
        }
        if( typeid(typename TSM::force_type) == typeid(TagForceLong)){
#ifdef SHRINK_SIZE_OF_SPJ_SORTED
            const S32 n_sp_tot = spj_recv_.size();
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
                    const U32 adr = tp_glb_[i].adr_ptcl_;
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
                    const U32 adr = tp_glb_[i].adr_ptcl_;
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
                const U32 adr = tp_glb_[i].adr_ptcl_;
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
#else // SHRINK_SIZE_OF_SPJ_SORTED
            spj_sorted_.resizeNoInitialize( spj_org_.size() );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_glb_tot_; i++){
                const U32 adr = tp_glb_[i].adr_ptcl_;
                if( GetMSB(adr) == 0){
                    epj_sorted_[i] = epj_org_[adr];
                }
                else{
                    spj_sorted_[i] = spj_org_[ClearMSB(adr)];
                }
            }
#endif
        }
        else{
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_glb_tot_; i++){
                const U32 adr = tp_glb_[i].adr_ptcl_;
                epj_sorted_[i] = epj_org_[adr];
            }
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"tp_glb_.size()="<<tp_glb_.size()<<" tp_buf_.size()="<<tp_buf_.size()<<std::endl;
        std::cout<<"epj_sorted_.size()="<<epj_sorted_.size()<<" spj_sorted_.size()="<<spj_sorted_.size()<<std::endl;
#endif
        time_profile_.morton_sort_global_tree += GetWtime() - time_offset;
        time_profile_.make_global_tree += GetWtime() - time_offset;
    }
#endif
    
    ///////////////////////////////////
    // LINK TREE CELLS USING "tp_loc_"
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    linkCellLocalTreeOnly(const bool reuse, const bool keep_fp_order){
        const F64 time_offset = GetWtime();
        if(reuse) return;
        LinkCell(tc_loc_,  adr_tc_level_partition_loc_, tp_loc_.getPointer(),
                 lev_max_loc_, n_loc_tot_, n_leaf_limit_);
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"tc_loc_.size()="<<tc_loc_.size()<<std::endl;
        std::cout<<"lev_max_loc_="<<lev_max_loc_<<std::endl;
#endif
        time_profile_.link_cell_local_tree += GetWtime() - time_offset;
        time_profile_.make_local_tree += GetWtime() - time_offset;
    }
    
    /////////////////////////////////
    // CALCULATE MOMENT OF LOCAL TREE 
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnly(){
        F64 wtime_offset = GetWtime();
        calcMomentLocalTreeOnlyImpl(typename TSM::search_type());
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif
        time_profile_.calc_moment_local_tree += GetWtime() - wtime_offset;
        time_profile_.make_local_tree_tot += GetWtime() - wtime_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLong){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLongScatter){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLongSymmetry){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLongCutoffScatter){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchLongCutoff){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchShortScatter){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchShortGather){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epi_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagSearchShortSymmetry){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                   epi_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }

    /////////////////////////////
    /// link cell global tree ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    linkCellGlobalTreeOnly(){
        const F64 time_offset = GetWtime();
        LinkCell(tc_glb_, adr_tc_level_partition_glb_,
                 tp_glb_.getPointer(), lev_max_glb_, n_glb_tot_,
                 n_leaf_limit_);
        /*
#ifdef SHRINK_SIZE_OF_SPJ_SORTED
        if( typeid(typename TSM::force_type) == typeid(TagForceLong)){
            U32 n_cnt_ep = 0;
            U32 n_cnt_sp = 0;
            for(S32 i=0; i<n_glb_tot_; i++){
                const U32 adr = tp_glb_[i].adr_ptcl_;
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
        }
#endif
        */
        
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
#if 1
        calcMomentGlobalTreeOnlyImpl(typename TSM::force_type());
#else
        calcMomentGlobalTreeOnlyImpl(typename TSM::search_type());
#endif
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
#endif
        time_profile_.calc_moment_global_tree += GetWtime() - time_offset;
        time_profile_.make_global_tree_tot = time_profile_.calc_moment_global_tree + time_profile_.make_global_tree;
    }

#if 1
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagForceLong){
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_);
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagForceShort){
        CalcMoment(adr_tc_level_partition_glb_, tc_glb_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_);
    }
#else
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchLong){
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchLongScatter){
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchLongSymmetry){
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_);
    }    

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchLongCutoff){
        CalcMomentLongGlobalTree
            (adr_tc_level_partition_glb_,  tc_glb_.getPointer(),
             tp_glb_.getPointer(),     epj_sorted_.getPointer(),
             spj_sorted_.getPointer(), lev_max_glb_,
             n_leaf_limit_);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchShortScatter){
        CalcMoment(adr_tc_level_partition_glb_, tc_glb_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchShortGather){
        CalcMoment(adr_tc_level_partition_glb_, tc_glb_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentGlobalTreeOnlyImpl(TagSearchShortSymmetry){
        CalcMoment(adr_tc_level_partition_glb_, tc_glb_.getPointer(),
                   epj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_);
    }
#endif
    
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
        time_profile_.calc_force__make_ipgroup += GetWtime() - time_offset;
        time_profile_.make_ipg   += GetWtime() - time_offset;
        //time_profile_.calc_force += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeIPGroupImpl(TagForceLong){
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"n_group_limit_="<<n_group_limit_<<std::endl;
#endif
        /*
        std::cerr<<"tc_loc_.size()= "<<tc_loc_.size()
                 <<" tc_glb_.size()= "<<tc_glb_.size()
                 <<std::endl;
        */
#if 1
        MakeIPGroupUseGLBTreeLong(ipg_, tc_loc_, tc_glb_, epi_sorted_, 0, 0, n_group_limit_, n_leaf_limit_); // NEW
#else
        MakeIPGroupLong(ipg_, tc_loc_, epi_sorted_, 0, n_group_limit_);
#endif
	const S32 n_ipg = ipg_.size();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif	//PARTICLE_SIMULATOR_THREAD_PARALLEL
	for(S32 i=0; i<n_ipg; i++){
	    const S32 n = ipg_[i].n_ptcl_;
	    const S32 adr = ipg_[i].adr_ptcl_;
	    ipg_[i].vertex_ = GetMinBoxSingleThread(epi_sorted_.getPointer(adr), n);
	}
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeIPGroupImpl(TagForceShort){
        MakeIPGroupShort(ipg_, tc_loc_, epi_sorted_, 0, n_group_limit_);
    }
    
    //////////////////////
    // ADD MOMENT AS SP //
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    addMomentAsSpGlobal(){
        F64 wtime_offset = GetWtime();
        /*
        std::cerr<<"A) epj_sorted_.size()= "
                 <<epj_sorted_.size()
                 <<" spj_sorted_.size()= "
                 <<spj_sorted_.size()
                 <<std::endl;
        */
        AddMomentAsSpImpl(typename TSM::force_type(), 
                          tc_glb_, 
                          spj_sorted_.size(),
                          spj_sorted_); 
        /*
        std::cerr<<"B) epj_sorted_.size()= "
                 <<epj_sorted_.size()
                 <<" spj_sorted_.size()= "
                 <<spj_sorted_.size()
                 <<std::endl;
        */
        time_profile_.add_moment_as_sp_global += GetWtime() - wtime_offset;
    }

    //////////////////////////
    //// WRITE BACK FORCE ////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tsys>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    writeBackForce(Tsys & sys,
                   const INTERACTION_LIST_MODE list_mode,
                   const bool keep_fp_order){
        const F64 time_offset = GetWtime();
        bool flag_reuse = false;
        if(list_mode == REUSE_LIST) flag_reuse = true;
        if( (keep_fp_order==false) && (flag_reuse==false) ){
            static Tsys sys_tmp;
            sys_tmp.createParticle(n_loc_tot_);
            sys_tmp.setNumberOfParticleLocal(n_loc_tot_);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                const S32 adr_org = adr_org_from_adr_sorted_[i];
                sys_tmp[i] = sys[adr_org];
            }
            sys.swapPtcl(sys_tmp);
        }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        //for(S32 i=0; i<n_loc_tot_; i++) sys[i].copyFromForce(force_org_[i]);
        for(S32 i=0; i<n_loc_tot_; i++) sys[i].copyFromForce(force_[i]);
        time_profile_.write_back += GetWtime() - time_offset;
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tep2>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcCenterAndLengthOfRootCellOpenNoMargenImpl(const Tep2 ep[]){
        const F64ort min_box  = GetMinBox(ep, n_loc_tot_);
        center_ = min_box.getCenter();
        const F64 tmp0 = (min_box.high_ - center_).getMax();
        const F64 tmp1 = (center_ - min_box.low_).getMax();
        length_ = std::max(tmp0, tmp1) * 2.0 * 1.001;
        pos_root_cell_.low_ = center_ - F64vec(length_*0.5);
        pos_root_cell_.high_ = center_ + F64vec(length_*0.5);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tep2>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcCenterAndLengthOfRootCellOpenWithMargenImpl(const Tep2 ep[]){
	const F64ort min_box  = GetMinBoxWithMargen(ep, n_loc_tot_);
	center_ = min_box.getCenter();
	const F64 tmp0 = (min_box.high_ - center_).getMax();
	const F64 tmp1 = (center_ - min_box.low_).getMax();
	length_ = std::max(tmp0, tmp1) * 2.0 * 1.001;
	pos_root_cell_.low_ = center_ - F64vec(length_*0.5);
	pos_root_cell_.high_ = center_ + F64vec(length_*0.5);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tep2>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcCenterAndLengthOfRootCellPeriodicImpl2(const Tep2 ep[]){
        //F64 rsearch_max_loc = std::numeric_limits<F64>::max() * -0.25;
	F64 rsearch_max_loc = -LARGE_FLOAT;
        F64ort box_loc;
        box_loc.initNegativeVolume();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	
#pragma omp parallel
#endif	
        {
            //F64 rsearch_max_loc_tmp = std::numeric_limits<F64>::max() * -0.25;
	    F64 rsearch_max_loc_tmp = -LARGE_FLOAT;
            F64ort box_loc_tmp;
            box_loc_tmp.init();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp for nowait
#endif	    
            for(S32 ip=0; ip<n_loc_tot_; ip++){
                rsearch_max_loc_tmp = (rsearch_max_loc_tmp > ep[ip].getRSearch()) ? rsearch_max_loc_tmp : ep[ip].getRSearch();
                box_loc_tmp.merge(ep[ip].getPos());
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp critical
#endif	    
            {
                rsearch_max_loc = rsearch_max_loc > rsearch_max_loc_tmp ? rsearch_max_loc : rsearch_max_loc_tmp;
                box_loc.merge(box_loc_tmp);
            }
        }
        F64 rsearch_max_glb = Comm::getMaxValue(rsearch_max_loc);
        F64vec xlow_loc = box_loc.low_;
        F64vec xhigh_loc = box_loc.high_;
        F64vec xlow_glb = Comm::getMinValue(xlow_loc);
        F64vec xhigh_glb = Comm::getMaxValue(xhigh_loc);

        xlow_glb -= F64vec(rsearch_max_glb);
        xhigh_glb += F64vec(rsearch_max_glb);
        center_ = (xhigh_glb + xlow_glb) * 0.5;
        F64 tmp0 = (xlow_glb - center_).applyEach(Abs<F64>()).getMax();
        F64 tmp1 = (xhigh_glb - center_).applyEach(Abs<F64>()).getMax();
        length_ = std::max(tmp0, tmp1) * 2.0 * 2.0;
        pos_root_cell_.low_ = center_ - F64vec(length_*0.5);
        pos_root_cell_.high_ = center_ + F64vec(length_*0.5);
    }

    /////////////////
    // CALC FORCE ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceOnly(Tfunc_ep_ep pfunc_ep_ep,
                  const S32 adr_ipg,
                  const bool clear){
        const S32 offset = ipg_[adr_ipg].adr_ptcl_;
        const S32 n_epi = ipg_[adr_ipg].n_ptcl_;
        const S32 ith = Comm::getThreadNum();
        const S32 n_epj = epj_for_force_[ith].size();
        const S32 n_tail = offset + n_epi;
        if(clear){
            for(S32 i=offset; i<n_tail; i++) force_[i].clear();
        }
        pfunc_ep_ep(epi_sorted_.getPointer(offset),     n_epi,
                    epj_for_force_[ith].getPointer(),   n_epj,
                    force_.getPointer(offset));
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep, class Tfunc_ep_sp>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceOnly(Tfunc_ep_ep pfunc_ep_ep,
                  Tfunc_ep_sp pfunc_ep_sp,
                  const S32 adr_ipg,
                  const bool clear){
        const S32 offset = ipg_[adr_ipg].adr_ptcl_;
        const S32 n_epi = ipg_[adr_ipg].n_ptcl_;
        const S32 ith = Comm::getThreadNum();
        const S32 n_epj = epj_for_force_[ith].size();
        const S32 n_spj = spj_for_force_[ith].size();
        const S32 n_tail = offset + n_epi;
        if(clear){
            for(S32 i=offset; i<n_tail; i++) force_[i].clear();
        }
        pfunc_ep_ep(epi_sorted_.getPointer(offset),     n_epi,
                    epj_for_force_[ith].getPointer(),   n_epj,
                    force_.getPointer(offset));
        pfunc_ep_sp(epi_sorted_.getPointer(offset),     n_epi,
                    spj_for_force_[ith].getPointer(),   n_spj,
                    force_.getPointer(offset));
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    copyForceOriginalOrder(const bool keep_fp_order){
        if(keep_fp_order){
            force_buf_.resizeNoInitialize(n_loc_tot_);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                const S32 adr = adr_org_from_adr_sorted_[i];
                force_buf_[adr] = force_[i];
            }
            force_.swap(force_buf_);
            //force_buf_.freeMem();
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceWalkOnly(Tfunc_ep_ep pfunc_ep_ep,
                      const bool clear){
        F64 time_offset = GetWtime();
        force_.resizeNoInitialize(n_loc_tot_);
        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        const S64 n_ipg = ipg_.size();
        n_walk_local_ += n_ipg;
        S64 ni_tmp = 0;
        S64 nj_tmp = 0;
        S64 n_interaction_ep_ep_tmp = 0;
        for(S32 i=0; i<Comm::getNumberOfThread(); i++) n_cell_open_[i] = 0;
        //PROFILE::Start(profile.calc_force);
        if(n_ipg > 0){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : ni_tmp, nj_tmp, n_interaction_ep_ep_tmp)
#endif	    
            for(S32 i=0; i<n_ipg; i++){
                makeInteractionList(i);
                ni_tmp += ipg_[i].n_ptcl_;
                nj_tmp += epj_for_force_[Comm::getThreadNum()].size();
                n_interaction_ep_ep_tmp += ipg_[i].n_ptcl_ * epj_for_force_[Comm::getThreadNum()].size();
                //calcForceOnly( pfunc_ep_ep, i, clear);
            }
            ni_ave_ = ni_tmp / n_ipg;
            nj_ave_ = nj_tmp / n_ipg;
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
            for(S32 i=1; i<Comm::getNumberOfThread(); i++) n_cell_open_[0] += n_cell_open_[i];
        }
        else{
            ni_ave_ = nj_ave_ = 0;
            n_walk_local_ = n_interaction_ep_ep_local_ = 0;
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"ipg_.size()="<<ipg_.size()<<std::endl;
#endif
        time_profile_.calc_force += GetWtime() - time_offset;
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
        //force_org_.resizeNoInitialize(n_loc_tot_);
        force_.resizeNoInitialize(n_loc_tot_);
        if(clear){
            //for(S32 i=0; i<n_loc_tot_; i++)force_org_[i].clear();
            for(S32 i=0; i<n_loc_tot_; i++)force_[i].clear();
        }
        Tepj * epj_tmp;
        S32 n_epj_tmp;
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        AllGatherParticle(epj_tmp, n_epj_tmp, epj_org_.getPointer(), n_loc_tot_, dinfo.getPosRootDomain().getFullLength(), pos_root_cell_, pa);
        //pfunc_ep_ep(epi_org_.getPointer(), n_loc_tot_, epj_tmp, n_epj_tmp, force_org_.getPointer());
        pfunc_ep_ep(epi_org_.getPointer(), n_loc_tot_, epj_tmp, n_epj_tmp, force_.getPointer());
        delete [] epj_tmp;
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
                                             epj_sorted_,   epj_neighbor_[id_thread], n_leaf_limit_);
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
					    epj_sorted_,   epj_neighbor_[id_thread], n_leaf_limit_);
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
					      epj_sorted_,   epj_neighbor_[id_thread],
                                              n_leaf_limit_);
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
	SearchNeighborListOneParticleScatter(pos_target,    tc_glb_.getPointer(),
                                             tp_glb_.getPointer(), adr, 
                                             epj_sorted_,   epj_neighbor_[id_thread], n_leaf_limit_);
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
                                              n_leaf_limit_);
	S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){
	return getNeighborListOneParticleImpl(typename TSM::search_type(), ptcl, epj);
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    freeMem(){
	tp_buf_.freeMem();
	tp_loc_.freeMem();
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
	epj_recv_.freeMem();
	spj_send_.freeMem();
	spj_recv_.freeMem();
        force_.freeMem();
	force_buf_.freeMem();
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

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    reallocMem(){
	tp_buf_.reallocMem();
	tp_loc_.reallocMem();
	tp_glb_.reallocMem();
	tc_loc_.reallocMem();
	tc_glb_.reallocMem();
	epi_sorted_.reallocMem();
	epi_org_.reallocMem();
	epj_sorted_.reallocMem();
	epj_org_.reallocMem();
	spj_sorted_.reallocMem();
	spj_org_.reallocMem();
	ipg_.reallocMem();
	epj_send_.reallocMem();
	epj_recv_.reallocMem();
	spj_send_.reallocMem();
	spj_recv_.reallocMem();
	force_.reallocMem();
	force_buf_.reallocMem();
	//force_sorted_.reallocMem();
	//force_org_.reallocMem();
	epjr_sorted_.reallocMem();
	epjr_send_.reallocMem();
	epjr_recv_.reallocMem();
	epjr_recv_1st_buf_.reallocMem();
	epjr_recv_2nd_buf_.reallocMem();
        const S32 n_thread = Comm::getNumberOfThread();
	for(S32 i=0; i<n_thread; i++){
	    epj_for_force_[i].reallocMem();
	    spj_for_force_[i].reallocMem();
	    epjr_send_buf_[i].reallocMem();
	    epjr_send_buf_for_scatter_[i].reallocMem();
	    epjr_recv_1st_sorted_[i].reallocMem();
	}
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    clearSizeOfArray(){
	tp_buf_.clearSize();
	tp_loc_.clearSize();
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
	epj_recv_.clearSize();
	spj_send_.clearSize();
	spj_recv_.clearSize();
	force_.clearSize();
	force_buf_.clearSize();
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
        const S32 adr_tree_sp_first = 0;
        const F64 r_crit_sq = 999.9;
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
                TargetBox<TSM> target_box;
                GetTargetBox<TSM>(ipg_[i], target_box);
                MakeListUsingTreeRecursiveTop
                    <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, Tspj,
                     TagWalkModeNormal, TagChopLeafTrue>
                    (tc_glb_,  adr_tc, tp_glb_,
                     epj_sorted_, adr_epj_tmp[ith],
                     spj_dummy,   adr_spj,
                     //pos_target_box,
                     target_box,
                     r_crit_sq, n_leaf_limit_,
                     adr_tree_sp_first, F64vec(0.0));
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
        const S32 adr_tree_sp_first = spj_org_.size();
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
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
                adr_ipg_tmp[ith].push_back(i);
                n_disp_epj_tmp[ith].push_back(n_ep_cum_prev);
                n_disp_spj_tmp[ith].push_back(n_sp_cum_prev);
                TargetBox<TSM> target_box;
                GetTargetBox<TSM>(ipg_[i], target_box);
                MakeListUsingTreeRecursiveTop
                    <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, Tspj,
                     TagWalkModeNormal, TagChopLeafTrue>
                    (tc_glb_,  adr_tc, tp_glb_,
                     epj_sorted_, adr_epj_tmp[ith],
                     spj_sorted_, adr_spj_tmp[ith],
                     target_box,
                     r_crit_sq, n_leaf_limit_,
                     adr_tree_sp_first, F64vec(0.0));
#if 0
                if(Comm::getRank()==0){
                    std::cerr<<"n_ep_cum_prev= "<<n_ep_cum_prev
                             <<" n_sp_cum_prev= "<<n_sp_cum_prev
                             <<std::endl;
                    F64 mass_tmp = 0.0;
                    for(S32 j=n_ep_cum_prev; j<adr_epj_tmp[ith].size(); j++){
                        mass_tmp += epj_sorted_[adr_epj_tmp[ith][j]].getCharge();
                    }
                    for(S32 j=n_sp_cum_prev; j<adr_spj_tmp[ith].size(); j++){
                        mass_tmp += spj_sorted_[adr_spj_tmp[ith][j]].getCharge();
                    }
                    std::cerr<<"mass_tmp= "<<mass_tmp<<std::endl;
                }
#endif
                
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
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalk(Tfunc_ep_ep pfunc_ep_ep,
                    const bool clear){
        F64 time_offset = GetWtime();
        force_.resizeNoInitialize(n_loc_tot_);
        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                //force_sorted_[i].clear();
                force_[i].clear();
            }
        }
        S64 n_interaction_ep_ep_tmp = 0;
        const S64 n_ipg = ipg_.size();
        if(n_ipg > 0){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : n_interaction_ep_ep_tmp)
#endif
            for(S32 i=0; i<n_ipg; i++){
                const S32 ith = Comm::getThreadNum();
                const S32 n_epi = ipg_[i].n_ptcl_;
                const S32 adr_epi_head = ipg_[i].adr_ptcl_;
                const S32 n_epj = interaction_list_.n_ep_[i];
                const S32 adr_epj_head = interaction_list_.n_disp_ep_[i];
                const S32 adr_epj_end  = interaction_list_.n_disp_ep_[i+1];
                n_interaction_ep_ep_tmp += ipg_[i].n_ptcl_ * n_epj;
                epj_for_force_[ith].resizeNoInitialize(n_epj);
                S32 n_cnt = 0;
                for(S32 j=adr_epj_head; j<adr_epj_end; j++, n_cnt++){
                    const S32 adr_epj = interaction_list_.adr_ep_[j];
                    epj_for_force_[ith][n_cnt] = epj_sorted_[adr_epj];
                }
                /*
                pfunc_ep_ep(epi_sorted_.getPointer(adr_epi_head),     n_epi,
                            epj_for_force_[ith].getPointer(),   n_epj,
                            force_sorted_.getPointer(adr_epi_head));
                */
                pfunc_ep_ep(epi_sorted_.getPointer(adr_epi_head),     n_epi,
                            epj_for_force_[ith].getPointer(),   n_epj,
                            force_.getPointer(adr_epi_head));
            }
        }
        n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
        time_profile_.calc_force += GetWtime() - time_offset;
    }

    
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep, class Tfunc_ep_sp>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalk(Tfunc_ep_ep pfunc_ep_ep,
                    Tfunc_ep_sp pfunc_ep_sp,
                    const bool clear){
        F64 time_offset = GetWtime();
        force_.resizeNoInitialize(n_loc_tot_);
        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                //force_sorted_[i].clear();
                force_[i].clear();
            }
        }
        S64 n_interaction_ep_ep_tmp = 0;
        S64 n_interaction_ep_sp_tmp = 0;
        const S64 n_ipg = ipg_.size();
        if(n_ipg > 0){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : n_interaction_ep_ep_tmp, n_interaction_ep_sp_tmp)
#endif
            for(S32 i=0; i<n_ipg; i++){
                const S32 ith = Comm::getThreadNum();
                const S32 n_epi = ipg_[i].n_ptcl_;
                const S32 adr_epi_head = ipg_[i].adr_ptcl_;
                
                const S32 n_epj = interaction_list_.n_ep_[i];
                const S32 adr_epj_head = interaction_list_.n_disp_ep_[i];
                const S32 adr_epj_end  = interaction_list_.n_disp_ep_[i+1];

                const S32 n_spj = interaction_list_.n_sp_[i];
                const S32 adr_spj_head = interaction_list_.n_disp_sp_[i];
                const S32 adr_spj_end  = interaction_list_.n_disp_sp_[i+1];
                n_interaction_ep_ep_tmp += ipg_[i].n_ptcl_ * n_epj;
                n_interaction_ep_sp_tmp += ipg_[i].n_ptcl_ * n_spj;
                epj_for_force_[ith].resizeNoInitialize(n_epj);
                spj_for_force_[ith].resizeNoInitialize(n_spj);
                S32 n_ep_cnt = 0;
                //F64 mass_tmp = 0.0;
                for(S32 j=adr_epj_head; j<adr_epj_end; j++, n_ep_cnt++){
                    const S32 adr_epj = interaction_list_.adr_ep_[j];
                    epj_for_force_[ith][n_ep_cnt] = epj_sorted_[adr_epj];
                    //mass_tmp += epj_sorted_[adr_epj].mass;
                }
                /*
                pfunc_ep_ep(epi_sorted_.getPointer(adr_epi_head),     n_epi,
                            epj_for_force_[ith].getPointer(),   n_epj,
                            force_sorted_.getPointer(adr_epi_head));
                */
                pfunc_ep_ep(epi_sorted_.getPointer(adr_epi_head),     n_epi,
                            epj_for_force_[ith].getPointer(),   n_epj,
                            force_.getPointer(adr_epi_head));
                S32 n_sp_cnt = 0;
                for(S32 j=adr_spj_head; j<adr_spj_end; j++, n_sp_cnt++){
                    const S32 adr_spj = interaction_list_.adr_sp_[j];
                    spj_for_force_[ith][n_sp_cnt] = spj_sorted_[adr_spj];
                    //mass_tmp += spj_sorted_[adr_spj].mass;
                }
                //std::cerr<<"mass_tmp= "<<mass_tmp<<std::endl;
                /*
                pfunc_ep_sp(epi_sorted_.getPointer(adr_epi_head),     n_epi,
                            spj_for_force_[ith].getPointer(),   n_spj,
                            force_sorted_.getPointer(adr_epi_head));
                */
                pfunc_ep_sp(epi_sorted_.getPointer(adr_epi_head),     n_epi,
                            spj_for_force_[ith].getPointer(),   n_spj,
                            force_.getPointer(adr_epi_head));
            }
        }
        n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
        n_interaction_ep_sp_local_ += n_interaction_ep_sp_tmp;
        time_profile_.calc_force += GetWtime() - time_offset;
    }
}
#include"tree_for_force_impl_ex_let.hpp"
#include"tree_for_force_impl_force.hpp"
