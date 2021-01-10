////////////////////////////////////////////////
/// implementaion of methods of TreeForForce ///

#ifdef SUNWAY
extern "C"{
    #include <athread.h>
    #include <cpe_func.h>
    void SLAVE_FUN(GetMinIpgBox)(void *);
}
#endif

namespace ParticleSimulator{
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tep2, class Tep3>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>
    ::scatterEP(S32 n_send[],
                S32 n_send_disp[],
                S32 n_recv[],
                S32 n_recv_disp[],
                ReallocatableArray<Tep2> & ep_send,  // send buffer
                ReallocatableArray<Tep2> & ep_recv,  // recv buffer
                ReallocatableArray<Tep2> * ep_send_buf,  // send buffer
                const ReallocatableArray<Tep3> & ep_org, // original
                const DomainInfo & dinfo){
        F64 time_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
        const F64ort outer_boundary_of_my_tree = tc_loc_[0].mom_.vertex_out_;

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"outer_boundary_of_my_tree="<<outer_boundary_of_my_tree<<std::endl;
#endif
        wtime_walk_LET_1st_ = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif	
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_proc_cum = 0;
            bool pa[DIMENSION];
            dinfo.getPeriodicAxis(pa);
            ep_send_buf[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    		    
#pragma omp for schedule(dynamic, 4)
#endif	    
            for(S32 ib=0; ib<n_proc; ib++){
                n_send[ib] = n_recv[ib] = 0;
                shift_image_domain_[ith].clearSize();
                if(dinfo.getBoundaryCondition() == BOUNDARY_CONDITION_OPEN){
                    shift_image_domain_[ith].push_back( (F64vec)(0.0));
                }
                else{
                    CalcNumberAndShiftOfImageDomain
                        (shift_image_domain_[ith], 
                         dinfo.getPosRootDomain().getFullLength(),
                         outer_boundary_of_my_tree,
                         dinfo.getPosDomain(ib),
                         pa);
                }
                S32 n_ep_offset = ep_send_buf[ith].size();
                const S32 n_image = shift_image_domain_[ith].size();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
                PARTICLE_SIMULATOR_PRINT_LINE_INFO();
                std::cout<<"n_image="<<n_image<<std::endl;
                std::cout<<"dinfo.getPosDomain(ib)="<<dinfo.getPosDomain(ib)<<std::endl;
#endif

                for(S32 ii=0; ii<n_image; ii++){
                    if(my_rank == ib && ii == 0) continue; // skip self image
                    F64ort pos_domain = dinfo.getPosDomain(ib).shift(shift_image_domain_[ith][ii]);
                    const S32 adr_tc_tmp = tc_loc_[0].adr_tc_;
#if 0
                    const S32 adr_root_cell = 0;
                    MakeListUsingOuterBoundaryIteration
                        (tc_loc_.getPointer(),  adr_root_cell,
                         ep_org.getPointer(),   ep_send_buf[ith],
                         pos_domain,            n_leaf_limit_,
                         -shift_image_domain_[ith][ii]);
#else
                    if( !tc_loc_[0].isLeaf(n_leaf_limit_) ){
                        /*
                          MakeListUsingOuterBoundary
                          (tc_loc_.getPointer(),  adr_tc_tmp,
                          ep_org.getPointer(),   ep_send_buf[ith],
                          pos_domain,            n_leaf_limit_,
                          -shift_image_domain_[ith][ii]);
                        */
                        // add M.I. 2016/03/23
                        MakeListUsingOuterBoundaryWithRootCellCheck
                            (tc_loc_.getPointer(),  adr_tc_tmp,
                             ep_org.getPointer(),   ep_send_buf[ith],
                             pos_domain,            n_leaf_limit_,
                             pos_root_cell_,
                             -shift_image_domain_[ith][ii]);
                    }
                    else{
                        const S32 n_tmp = tc_loc_[0].n_ptcl_;
                        S32 adr_tmp = tc_loc_[0].adr_ptcl_;
                        for(S32 iii=0; iii<n_tmp; iii++, adr_tmp++){
                            const F64vec pos_tmp = ep_org[adr_tmp].getPos();
                            const F64 size_tmp = ep_org[adr_tmp].getRSearch();
                            const F64 dis_sq_tmp = pos_domain.getDistanceMinSQ(pos_tmp);
#ifdef ORIGINAL_SCATTER_MODE
                            if(dis_sq_tmp > size_tmp*size_tmp) continue;
#endif
                            //if( pos_root_cell_.notOverlapped(pos_tmp-shift_image_domain_[ith][ii]) ) continue;  // add M.I. 2016/03/12
			    if( pos_root_cell_.notContained(pos_tmp-shift_image_domain_[ith][ii]) ) continue;  // add M.I. 2016/03/12
                            ep_send_buf[ith].increaseSize();
                            ep_send_buf[ith].back() = ep_org[adr_tmp];
                            const F64vec pos_new = ep_send_buf[ith].back().getPos() - shift_image_domain_[ith][ii];
                            ep_send_buf[ith].back().setPos(pos_new);
                        }
                    }
#endif
                }
                n_send[ib] = ep_send_buf[ith].size() - n_ep_offset;
                id_proc_send_[ith][n_proc_cum++] = ib;  // new
            } // omp for
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif	    
            {
                n_send_disp[0] = 0;
                for(S32 i=0; i<n_proc; i++){
                    n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                }
                ep_send.resizeNoInitialize(n_send_disp[n_proc]);
            }
            S32 n_ep_cnt = 0;
            for(S32 ib=0; ib<n_proc_cum; ib++){
                const S32 id = id_proc_send_[ith][ib];
                const S32 adr_ep_tmp = n_send_disp[id];
                const S32 n_ep_tmp = n_send[id];
                for(int ip=0; ip<n_ep_tmp; ip++){
                    ep_send[adr_ep_tmp+ip] = ep_send_buf[ith][n_ep_cnt++];
                }
            }
        } // omp parallel scope
        wtime_walk_LET_1st_ = GetWtime() - wtime_walk_LET_1st_;
        time_profile_.make_LET_1st += GetWtime() - time_offset;
        //time_profile_.exchange_LET_tot += GetWtime() - time_offset;
        time_offset = GetWtime();

        Tcomm_scatterEP_tmp_ = GetWtime();
        Comm::allToAll(n_send, 1, n_recv); // TEST
        n_recv_disp[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
        }

        ep_recv.resizeNoInitialize( n_recv_disp[n_proc] );

        Comm::allToAllV(ep_send.getPointer(), n_send, n_send_disp,
                        ep_recv.getPointer(), n_recv, n_recv_disp);

        Tcomm_scatterEP_tmp_ = GetWtime() - Tcomm_scatterEP_tmp_;

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"ep_send.size()="<<ep_send.size()<<std::endl;
        std::cout<<"ep_recv.size()="<<ep_recv.size()<<std::endl;
#endif
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
        //time_profile_.exchange_LET_tot += GetWtime() - time_offset;
    }





    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    size_t TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>
    ::getMemSizeUsed() const {
        size_t tmp = 0;
        for(int i=0; i<Comm::getNumberOfThread(); i++){
            tmp = id_ep_send_buf_[i].getMemSize() + id_sp_send_buf_[i].getMemSize()
                + epj_for_force_[i].getMemSize() + spj_for_force_[i].getMemSize()
                + ep_send_buf_for_scatter_[i].getMemSize() + shift_image_domain_[i].getMemSize()
                + epjr_send_buf_[i].getMemSize() + epjr_send_buf_for_scatter_[i].getMemSize() 
                + epjr_recv_1st_sorted_[i].getMemSize() + epj_send_buf_[i].getMemSize() 
                + id_ptcl_send_[i].getMemSize() + shift_image_box_[i].getMemSize()
                + ip_disp_[i].getMemSize() + tp_scatter_[i].getMemSize()
                //+ tc_recv_1st_[i].getMemSize()
                + epj_recv_1st_sorted_[i].getMemSize();
        }
        return tmp
#ifndef USE_MEMORY_POOL
            + tp_buf_.getMemSize()
            + tp_buf_2_.getMemSize()
#endif
#ifndef REMOVE_TP_LOC
            + tp_loc_.getMemSize()
#endif
            + tp_glb_.getMemSize()
            + tc_loc_.getMemSize() + tc_glb_.getMemSize()
            //+ epi_sorted_.getMemSize() + epi_org_.getMemSize()
            + epi_sorted_.getMemSize()
            + epj_sorted_.getMemSize() + epj_org_.getMemSize()
            + spj_sorted_.getMemSize() + spj_org_.getMemSize()
            + ipg_.getMemSize()
            + epj_send_.getMemSize() + epj_recv_.getMemSize()
            + spj_send_.getMemSize() + spj_recv_.getMemSize()
            + force_sorted_.getMemSize() + force_org_.getMemSize()
            + epjr_sorted_.getMemSize() + epjr_send_.getMemSize()
            + epjr_recv_.getMemSize() + epjr_recv_1st_buf_.getMemSize()
            + epjr_recv_2nd_buf_.getMemSize()
            + epj_recv_1st_buf_.getMemSize() + epj_recv_2nd_buf_.getMemSize();
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

        is_initialized_ = true;
        n_glb_tot_ = n_glb_tot;
        theta_ = theta;
        n_leaf_limit_ = n_leaf_limit;
        n_group_limit_ = n_group_limit;
        //lev_max_ = 0;
        lev_max_loc_ = lev_max_glb_ = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        const S64 n_proc = Comm::getNumberOfProc();
        const S64 np_ave = (n_glb_tot_ / n_proc);
        n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = n_walk_local_ = 0;
        n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
        wtime_exlet_comm_ = wtime_exlet_a2a_ = wtime_exlet_a2av_ = 0.0;
        wtime_walk_LET_1st_ = wtime_walk_LET_2nd_ = 0.0;
        ni_ave_ = nj_ave_ = 0;
        Comm::barrier();
        if(Comm::getRank() == 0){
            std::cerr<<"np_ave="<<np_ave<<std::endl;
        }

        const F64 np_one_dim = pow( ((F64)np_ave)*1.0001, 1.0/DIMENSION) + 4;
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
        n_surface_for_comm_ = (4*np_one_dim)*6;
#else
        n_surface_for_comm_ = (S64)(((F64)(6)*np_one_dim*np_one_dim+(F64)(8)*np_one_dim)*(F64)(6));
#endif

#ifndef REMOVE_EPI_SORTED
        epi_sorted_.reserve( np_ave+(np_ave+3)/3 + 1024 );
#endif


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
                //PARTICLE_SIMULATOR_PRINT_ERROR("theta must be >= 0.0");
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

        id_ep_send_buf_   = new ReallocatableArray<S32>[n_thread];
        id_sp_send_buf_   = new ReallocatableArray<S32>[n_thread];
        epj_for_force_    = new ReallocatableArray<Tepj>[n_thread];
        spj_for_force_    = new ReallocatableArray<Tspj>[n_thread];
        adr_epj_for_force_ = new ReallocatableArray<S32>[n_thread];
        adr_spj_for_force_ = new ReallocatableArray<S32>[n_thread];

        id_epj_recorder_for_force_ = new ReallocatableArray<S32>[n_thread];
        id_spj_recorder_for_force_ = new ReallocatableArray<S32>[n_thread];
        n_epi_recorder_for_force_ = new ReallocatableArray<S32>[n_thread];
        n_disp_epi_recorder_for_force_ = new ReallocatableArray<S32>[n_thread];
        n_epj_recorder_for_force_ = new ReallocatableArray<S32>[n_thread];
        n_spj_recorder_for_force_ = new ReallocatableArray<S32>[n_thread];
        n_disp_epj_recorder_for_force_ = new ReallocatableArray<S32>[n_thread];
        n_disp_spj_recorder_for_force_ = new ReallocatableArray<S32>[n_thread];

#ifdef SUNWAY
        // The followings are used in the copy from epj_sorted_loc_ to epj_sorted_
        adr_epj_loc_group_head_ = new ReallocatableArray<S32>[64];
        adr_epj_glb_group_head_ = new ReallocatableArray<S32>[64];
        group_size_epj_loc_     = new ReallocatableArray<S32>[64];
        group_size_epj_glb_     = new ReallocatableArray<S32>[64];
#endif
        
        id_proc_send_ = new S32*[n_thread];
        rank_proc_send_ =  new ReallocatableArray<S32>[n_thread]; // new version of id_proc_send
        shift_image_domain_ = new ReallocatableArray<F64vec>[n_thread];
        id_ptcl_send_ = new ReallocatableArray<S32>[n_thread];
        ip_disp_ = new ReallocatableArray<S32>[n_thread];
        tp_scatter_ = new ReallocatableArray< TreeParticle >[n_thread];

        //tc_recv_1st_ = new ReallocatableArray< TreeCell< Tmomloc > >[n_thread];

        shift_image_box_ = new ReallocatableArray<F64vec>[n_thread];// note
        ep_send_buf_for_scatter_ = new ReallocatableArray<Tepj>[n_thread];
        epj_send_buf_ = new ReallocatableArray<Tepj>[n_thread]; // note
        epj_recv_1st_sorted_ = new ReallocatableArray< Tepj >[n_thread];
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
        if( (FORCE_TYPE)(typename TSM::force_type()).force_type == FORCE_TYPE_LONG){
            if(theta_ > 0.0){
                const S64 n_tmp = epi_sorted_.capacity() + (S64)((F64)2000 * pow((0.5 / theta_), DIMENSION));
                const S64 n_new = std::min( std::min(n_tmp, n_glb_tot_+(S64)(100)), (S64)(20000) );
#ifndef USE_MEMORY_POOL
                epj_org_.reserve( n_new + 10000);
                spj_org_.reserve( epj_org_.capacity() );
#endif
                epj_sorted_.reserve( n_new + 10000);
                spj_sorted_.reserve( epj_sorted_.capacity() );
                epj_sorted_loc_.reserve( n_new + 10000);
#ifndef USE_MEMORY_POOL
                spj_sorted_loc_.reserve( epj_sorted_loc_.capacity() );
#endif

            }
            else if(theta == 0.0) {
#ifndef USE_MEMORY_POOL
                epj_org_.reserve(n_glb_tot_ + 100);
                spj_org_.reserve(n_glb_tot_ + 100);
#endif
                epj_sorted_.reserve(n_glb_tot_ + 100);
                spj_sorted_.reserve(n_glb_tot_ + 100);
                epj_sorted_loc_.reserve(n_glb_tot_ + 100);
#ifndef USE_MEMORY_POOL
                spj_sorted_loc_.reserve(n_glb_tot_ + 100);
#endif
            }
            for(S32 i=0; i<n_thread; i++){
                id_ep_send_buf_[i].reserve(100);
                id_sp_send_buf_[i].reserve(100);
            }
            spj_send_.reserve(n_proc+1000);
            spj_recv_.reserve(n_proc+1000);
            n_sp_send_ = new S32[n_proc];
            n_sp_recv_ = new S32[n_proc];
            n_sp_send_disp_ = new S32[n_proc+1];
            n_sp_recv_disp_ = new S32[n_proc+1];
            for(S32 i=0; i<n_thread; i++){
                epj_for_force_[i].reserve(10000);
                spj_for_force_[i].reserve(10000);
            }
            n_ep_sp_send_ = new S32[n_proc * 2];
            n_ep_sp_recv_ = new S32[n_proc * 2];
        }
        else{
            // FOR SHORT MODE
#ifndef USE_MEMORY_POOL
            epj_org_.reserve( epi_sorted_.capacity() + n_surface_for_comm_ );
#endif
            epj_sorted_.reserve( epj_org_.capacity() );
            for(S32 i=0; i<n_thread; i++){
                id_ep_send_buf_[i].reserve(n_surface_for_comm_ * 2 / n_thread);
            }
            for(S32 i=0; i<n_thread; i++){
                epj_for_force_[i].reserve(n_surface_for_comm_ * 2 / n_thread);
                spj_for_force_[i].reserve(1);
            }
            n_ep_sp_send_ = n_ep_sp_recv_ = NULL;
        }
        if(typeid(TSM) == typeid(SEARCH_MODE_LONG_SCATTER)){
            for(S32 ith=0; ith<n_thread; ith++){
                epj_neighbor_[ith].reserve(10);
                //epj_neighbor_[ith].reserve(epj_org_.capacity() / n_thread * 2);
            }
        }
        Comm::barrier();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT		
        if(Comm::getRank() == 0){
            std::cerr<<"used mem size for tree(1)="<<this->getMemSizeUsed()*1e-9<<"[GB]"<<std::endl;
        }
#endif
        //tp_loc_.reserve( epi_org_.capacity() );
#ifndef REMOVE_TP_LOC
        tp_loc_.reserve( epi_sorted_.capacity() );
#endif
        tp_glb_.reserve( epj_org_.capacity() );
#ifndef USE_MEMORY_POOL
        tp_buf_.reserve( epj_org_.capacity() );
#endif
#ifndef USE_MEMORY_POOL        
    #ifdef REMOVE_TP_LOC
        tc_loc_.reserve( 10 );
    #else 
        tc_loc_.reserve( tp_loc_.capacity() / n_leaf_limit_ * N_CHILDREN );
    #endif
        tc_glb_.reserve( tp_glb_.capacity() / n_leaf_limit_ * N_CHILDREN );
#endif // USE_MEMORY_POOL
        //ipg_.reserve( std::min(epi_org_.capacity()/n_group_limit_*4, epi_org_.capacity()) );
        ipg_.reserve( std::min(epi_sorted_.capacity()/n_group_limit_*4, epi_sorted_.capacity()) );

        Comm::barrier();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT	
        if(Comm::getRank() == 0){
            std::cerr<<"used mem size for tree(2)="<<this->getMemSizeUsed()*1e-9<<"[GB]"<<std::endl;
        }
#endif
        for(S32 i=0; i<n_thread; i++){
            id_proc_send_[i] = new S32[n_proc];
            rank_proc_send_[i].reserve(n_proc);
        }
        //epj_send_.reserve(n_surface_for_comm_);
        //epj_recv_.reserve(n_surface_for_comm_);
        epj_send_.reserve(10);
        epj_recv_.reserve(10);

        
        n_ep_send_ = new S32[n_proc];
        n_ep_recv_ = new S32[n_proc];
        n_ep_send_disp_ = new S32[n_proc+1];
        n_ep_recv_disp_ = new S32[n_proc+1];

        //force_org_.reserve(epi_sorted_.capacity());
        //force_sorted_.reserve(epi_sorted_.capacity());        

        Comm::barrier();	
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT	
        if(Comm::getRank() == 0){
            std::cerr<<"used mem size for tree(3)="<<this->getMemSizeUsed()*1e-9<<"[GB]"<<std::endl;
        }
#endif
        // new variables for commnuication of LET
        // for scatterEP

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        req_send_ = new MPI::Request[n_proc];
        req_recv_ = new MPI::Request[n_proc];
#endif

        if( typeid(TSM) == typeid(SEARCH_MODE_SCATTER) 
            || typeid(TSM) == typeid(SEARCH_MODE_SYMMETRY) 
            || typeid(TSM) == typeid(SEARCH_MODE_GATHER)){
            for(S32 i=0; i<n_thread; i++){
                shift_image_domain_[i].reserve(5*5*5);
            }
        }

        if( typeid(TSM) == typeid(SEARCH_MODE_SYMMETRY)
            || typeid(TSM) == typeid(SEARCH_MODE_GATHER) ){
            n_epj_recv_1st_ = new S32[n_proc];
            n_epj_recv_disp_1st_ = new S32[n_proc+1];
            n_epj_recv_2nd_ = new S32[n_proc];
            n_epj_recv_disp_2nd_ = new S32[n_proc+1];
            id_proc_src_ = new S32[n_proc];
            id_proc_dest_ = new S32[n_proc];
            adr_tc_level_partition_recv_1st_ = new S32*[n_thread];
            for(S32 i=0; i<n_thread; i++){
                id_ptcl_send_[i].reserve(n_surface_for_comm_);
                shift_image_box_[i].reserve(5*5*5);
                ip_disp_[i].reserve(3*3*3);
                tp_scatter_[i].reserve(n_surface_for_comm_);
                //tc_recv_1st_[i].reserve(n_surface_for_comm_*4);
                adr_tc_level_partition_recv_1st_[i] = new S32[TREE_LEVEL_LIMIT+2];
            }
        }

        if( typeid(TSM) == typeid(SEARCH_MODE_SCATTER) 
            || typeid(TSM) == typeid(SEARCH_MODE_SYMMETRY) ){
            for(S32 i=0; i<n_thread; i++){
                ep_send_buf_for_scatter_[i].reserve(n_surface_for_comm_);
            }
        }

        // for symmetry
        if( typeid(TSM) == typeid(SEARCH_MODE_SYMMETRY) ){
            // LET 2nd
            epj_recv_1st_buf_.reserve(n_surface_for_comm_);
            epj_recv_2nd_buf_.reserve(n_surface_for_comm_);
            for(S32 i=0; i<n_thread; i++){
                epj_send_buf_[i].reserve(n_surface_for_comm_);
                epj_recv_1st_sorted_[i].reserve(n_surface_for_comm_*4);
            }
        }

        if( typeid(TSM) == typeid(SEARCH_MODE_GATHER) ){
            epjr_sorted_.reserve(epj_org_.capacity());
            epjr_send_.reserve(n_surface_for_comm_);
            epjr_recv_.reserve(n_surface_for_comm_);
            epjr_recv_1st_buf_.reserve(n_surface_for_comm_);
            epjr_recv_2nd_buf_.reserve(n_surface_for_comm_);
            for(S32 i=0; i<n_thread; i++){
                epjr_send_buf_[i].reserve(n_surface_for_comm_);
                epjr_send_buf_for_scatter_[i].reserve(n_surface_for_comm_);
                epjr_recv_1st_sorted_[i].reserve(n_surface_for_comm_);
            }
        }
        Comm::barrier();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT	
        if(Comm::getRank() == 0){
            std::cerr<<"used mem size for tree(3)="<<this->getMemSizeUsed()*1e-9<<"[GB]"<<std::endl;
        }
#endif
    }
    
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
        //epi_org_.resizeNoInitialize(n_loc_tot_);
        epj_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for
#endif	    
            for(S32 i=0; i<nloc; i++){
                //epi_org_[i].copyFromFP( psys[i] );
                epj_org_[i].copyFromFP( psys[i] );
            }
        }
        else{
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for
#endif	    
            for(S32 i=0; i<nloc; i++){
                //epi_org_[i+offset].copyFromFP( psys[i] );
                epj_org_[i+offset].copyFromFP( psys[i] );
            }
        }
        //std::cout<<"step d"<<std::endl;
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

    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    linkCellLocalTreeOnly(){
        const F64 time_offset = GetWtime();
        //#ifdef SUNWAY
#if 1
#ifdef REMOVE_TP_LOC
        LinkCellSunWay(tc_loc_,  adr_tc_level_partition_loc_, tp_glb_.getPointer(), lev_max_loc_, n_loc_tot_, n_leaf_limit_);
#else
        LinkCellSunWay(tc_loc_,  adr_tc_level_partition_loc_, tp_loc_.getPointer(), lev_max_loc_, n_loc_tot_, n_leaf_limit_);
#endif
#else
        LinkCell(tc_loc_,  adr_tc_level_partition_loc_, tp_loc_.getPointer(), lev_max_loc_, n_loc_tot_, n_leaf_limit_);
#endif
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
        calcMomentLocalTreeOnlyImpl(typename TSM::ex_let_stage_type());
        time_profile_.calc_moment_local_tree += GetWtime() - time_offset;
        time_profile_.make_local_tree_tot = time_profile_.calc_moment_local_tree + time_profile_.make_local_tree;
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagExLetOneStage){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(), epj_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcMomentLocalTreeOnlyImpl(TagExLetTwoStage){
        CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(), epi_sorted_.getPointer(), lev_max_loc_, n_leaf_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTree(const DomainInfo & dinfo){
        if(typeid(TSM) == typeid(SEARCH_MODE_LONG)
           && dinfo.getBoundaryCondition() != BOUNDARY_CONDITION_OPEN){
            PARTICLE_SIMULATOR_PRINT_ERROR("The forces w/o cutoff can be evaluated only under the open boundary condition");
            Abort(-1);
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        if(tc_loc_[0].n_ptcl_ <= 0){
            std::cout<<"The number of particles of this process is 0."<<std::endl;
            std::cout<<"tc_loc_[0].n_ptcl_="<<tc_loc_[0].n_ptcl_<<std::endl;
        }
#endif
        exchangeLocalEssentialTreeImpl(typename TSM::search_type(), dinfo);
        time_profile_.exchange_LET_tot = time_profile_.make_LET_1st + time_profile_.exchange_LET_1st + time_profile_.make_LET_2nd + time_profile_.exchange_LET_2nd;
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epj_send_.size()="<<epj_send_.size()<<" spj_send_.size()="<<spj_send_.size()<<std::endl;
        std::cout<<"epj_recv_.size()="<<epj_recv_.size()<<" spj_recv_.size()="<<spj_recv_.size()<<std::endl;
#endif
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchLong, const DomainInfo & dinfo){
        F64 time_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	
#pragma omp parallel
#endif	
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_proc_cum = 0;
            id_ep_send_buf_[ith].reserve(1000);
            id_sp_send_buf_[ith].reserve(1000);
            id_ep_send_buf_[ith].resizeNoInitialize(0);
            id_sp_send_buf_[ith].resizeNoInitialize(0);
            const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
            const S32 adr_tc_tmp = tc_loc_[0].adr_tc_;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp for schedule(dynamic, 4)
#endif	    
            for(S32 ib=0; ib<n_proc; ib++){
                n_ep_send_[ib] = n_sp_send_[ib] = n_ep_recv_[ib] = n_sp_recv_[ib] = 0;
                if(my_rank == ib) continue;
                const F64ort pos_target_domain = dinfo.getPosDomain(ib);
                const S32 n_ep_cum = id_ep_send_buf_[ith].size();
                const S32 n_sp_cum = id_sp_send_buf_[ith].size();
                if(!tc_loc_[0].isLeaf(n_leaf_limit_)){
                    SearchSendParticleLong<TreeCell<Tmomloc>, Tepj>
                        (tc_loc_,               adr_tc_tmp,    epj_sorted_,
                         id_ep_send_buf_[ith],  id_sp_send_buf_[ith],
                         pos_target_domain,  r_crit_sq,   n_leaf_limit_);
                }
                else{
                    F64vec pos_tmp = tc_loc_[0].mom_.getPos();
                    if( (pos_target_domain.getDistanceMinSQ(pos_tmp) <= r_crit_sq * 4.0) ){
                        const S32 n_tmp = tc_loc_[0].n_ptcl_;
                        S32 adr_tmp = tc_loc_[0].adr_ptcl_;
                        for(S32 ip=0; ip<n_tmp; ip++){
                            id_ep_send_buf_[ith].push_back(adr_tmp++);
                        }
                    }
                    else{
                        id_sp_send_buf_[ith].push_back(0); // set root
                    }
                }

                n_ep_send_[ib] = id_ep_send_buf_[ith].size() - n_ep_cum;
                n_sp_send_[ib] = id_sp_send_buf_[ith].size() - n_sp_cum;
                id_proc_send_[ith][n_proc_cum++] = ib;
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp single
#endif	    
            {
                n_ep_send_disp_[0] = n_sp_send_disp_[0] = 0;
                for(S32 ib=0; ib<n_proc; ib++){
                    n_ep_send_disp_[ib+1] = n_ep_send_disp_[ib] + n_ep_send_[ib];
                    n_sp_send_disp_[ib+1] = n_sp_send_disp_[ib] + n_sp_send_[ib];
                }
                epj_send_.resizeNoInitialize( n_ep_send_disp_[n_proc] );
                spj_send_.resizeNoInitialize( n_sp_send_disp_[n_proc] );
            }
            S32 n_ep_cnt = 0;
            S32 n_sp_cnt = 0;
            for(S32 ib=0; ib<n_proc_cum; ib++){
                const S32 id = id_proc_send_[ith][ib];
                S32 adr_ep_tmp = n_ep_send_disp_[id];
                const S32 n_ep_tmp = n_ep_send_[id];
                for(int ip=0; ip<n_ep_tmp; ip++){
                    const S32 id_ep = id_ep_send_buf_[ith][n_ep_cnt++];
                    epj_send_[adr_ep_tmp++] = epj_sorted_[id_ep];
                }
                S32 adr_sp_tmp = n_sp_send_disp_[id];
                const S32 n_sp_tmp = n_sp_send_[id];
                for(int ip=0; ip<n_sp_tmp; ip++){
                    const S32 id_sp = id_sp_send_buf_[ith][n_sp_cnt++];
                    spj_send_[adr_sp_tmp++].copyFromMoment(tc_loc_[id_sp].mom_);
                }
            }
        } // omp parallel scope

        time_profile_.make_LET_1st += GetWtime() - time_offset;
        time_offset = GetWtime();

        wtime_exlet_comm_ = GetWtime();

        for(S32 i=0; i<n_proc; i++){
            n_ep_sp_send_[2*i] = n_ep_send_[i];
            n_ep_sp_send_[2*i+1] = n_sp_send_[i];
        }
        wtime_exlet_a2a_ = GetWtime();
#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<S32, 2> comm_a2a_2d;
        comm_a2a_2d.execute(n_ep_sp_send_, 2, n_ep_sp_recv_);
#else
        Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_); // TEST
#endif //FAST_ALL_TO_ALL_FOR_K


        wtime_exlet_a2a_ = GetWtime() - wtime_exlet_a2a_;

        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_[i] = n_ep_sp_recv_[2*i];
            n_sp_recv_[i] = n_ep_sp_recv_[2*i+1];
        }
        n_ep_recv_disp_[0] = n_sp_recv_disp_[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_disp_[i+1] = n_ep_recv_disp_[i] + n_ep_recv_[i];
            n_sp_recv_disp_[i+1] = n_sp_recv_disp_[i] + n_sp_recv_[i];
        }
        epj_recv_.resizeNoInitialize( n_ep_recv_disp_[n_proc] );
        spj_recv_.resizeNoInitialize( n_sp_recv_disp_[n_proc] );

        wtime_exlet_a2av_ = GetWtime();
#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<Tepj, 2> comm_a2a_epj_2d;
        static CommForAllToAll<Tspj, 2> comm_a2a_spj_2d;
        comm_a2a_epj_2d.executeV(epj_send_, epj_recv_, n_ep_send_, n_ep_recv_);
        comm_a2a_spj_2d.executeV(spj_send_, spj_recv_, n_sp_send_, n_sp_recv_);
#else
        Comm::allToAllV(epj_send_.getPointer(), n_ep_send_, n_ep_send_disp_,
                        epj_recv_.getPointer(), n_ep_recv_, n_ep_recv_disp_);
        Comm::allToAllV(spj_send_.getPointer(), n_sp_send_, n_sp_send_disp_,
                        spj_recv_.getPointer(), n_sp_recv_, n_sp_recv_disp_);
#endif //FAST_ALL_TO_ALL_FOR_K
        wtime_exlet_a2av_ = GetWtime() - wtime_exlet_a2av_;
        wtime_exlet_comm_ = GetWtime() - wtime_exlet_comm_;
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;
        n_let_ep_send_1st_ += (CountT)epj_send_.size();
        n_let_ep_recv_1st_ += (CountT)epj_recv_.size();
        n_let_sp_send_1st_ += (CountT)spj_send_.size();
        n_let_sp_recv_1st_ += (CountT)spj_recv_.size();
    }



    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchLongCutoff, const DomainInfo & dinfo){
        F64 time_offset = GetWtime();
        //std::cout<<"SEARCH_MODE_LONG_CUTOFF"<<std::endl;
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
        static bool first = true;
        static ReallocatableArray<Tepj> * epj_send_buf;
        static ReallocatableArray<Tspj> * spj_send_buf;
        if(first){
            const S32 n_thread = Comm::getNumberOfThread();
            epj_send_buf = new ReallocatableArray<Tepj>[n_thread];
            spj_send_buf = new ReallocatableArray<Tspj>[n_thread];
            for(S32 i=0; i<n_thread; i++){
                epj_send_buf[i].reserve(1000);
                spj_send_buf[i].reserve(1000);
            }
            first = false;
        }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	
#pragma omp parallel
#endif	
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_proc_cum = 0;
            epj_send_buf[ith].clearSize();
            spj_send_buf[ith].clearSize();
            S32 n_epj_cum = 0;
            S32 n_spj_cum = 0;
            const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
            const S32 adr_tc_tmp = tc_loc_[0].adr_tc_;
            const F64 r_cut_sq  = epj_sorted_[0].getRSearch() * epj_sorted_[0].getRSearch();
            ReallocatableArray<F64vec> shift_image_domain(5*5*5);
            bool periodic_axis[DIMENSION];
            dinfo.getPeriodicAxis(periodic_axis);
            const F64ort outer_boundary_of_my_tree = tc_loc_[0].mom_.vertex_out_;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp for schedule(dynamic, 4)
#endif	    
            for(S32 i=0; i<n_proc; i++){
                n_ep_send_[i] = n_sp_send_[i] = n_ep_recv_[i] = n_sp_recv_[i] = 0;
                shift_image_domain.clearSize();
                if(dinfo.getBoundaryCondition() == BOUNDARY_CONDITION_OPEN){
                    shift_image_domain.push_back( F64vec(0.0) );
                }
                else{
                    CalcNumberAndShiftOfImageDomain
                        (shift_image_domain, dinfo.getPosRootDomain().getFullLength(),
                         outer_boundary_of_my_tree, dinfo.getPosDomain(i), periodic_axis);
                }
                S32 n_image = shift_image_domain.size();
                for(S32 j = 0; j < n_image; j++) {
                    if(my_rank == i && j == 0) continue;
                    const F64ort pos_target_domain = dinfo.getPosDomain(i).shift(shift_image_domain[j]);
                    const F64ort cell_box = pos_root_cell_;
                    if( !tc_loc_[0].isLeaf(n_leaf_limit_) ){
/*
                        SearchSendParticleLongCutoff<TreeCell<Tmomloc>, Tepj, Tspj>
                            (tc_loc_,       adr_tc_tmp,
                             epj_sorted_,   epj_send_buf[ith],
                             spj_send_buf[ith],   cell_box,
                             pos_target_domain,   r_crit_sq,
                             r_cut_sq,            n_leaf_limit_,
                             - shift_image_domain[j]);
*/

                        SearchSendParticleLongCutoffWithRootCellCheck
                            <TreeCell<Tmomloc>, Tepj, Tspj>
                            (tc_loc_,       adr_tc_tmp,
                             epj_sorted_,   epj_send_buf[ith],
                             spj_send_buf[ith],   cell_box,
                             pos_target_domain,   r_crit_sq,
                             r_cut_sq,            n_leaf_limit_,
                             pos_root_cell_,
                             - shift_image_domain[j]);

                    }
                    else{
                        const F64 dis_sq_cut = pos_target_domain.getDistanceMinSQ(cell_box);
                        if(dis_sq_cut <= r_cut_sq){
                            const F64vec pos_tmp = tc_loc_[0].mom_.getPos();
                            const F64 dis_sq_crit = pos_target_domain.getDistanceMinSQ(pos_tmp);
                            if(dis_sq_crit <= r_crit_sq*4.0){
                                const S32 n_tmp = tc_loc_[0].n_ptcl_;
                                S32 adr_ptcl_tmp = tc_loc_[0].adr_ptcl_;
                                for(S32 ip=0; ip<n_tmp; ip++){
                                    //const F64vec pos_new = epj_send_buf[ith].back().getPos() - shift_image_domain[j];
                                    //if( pos_root_cell_.notOverlapped( epj_sorted_[adr_ptcl_tmp].getPos()-shift_image_domain[j] ) ) continue;  // added by M.I. 2016/03/12
				    if( pos_root_cell_.notContained( epj_sorted_[adr_ptcl_tmp].getPos()-shift_image_domain[j] ) ) continue;  // added by M.I. 2016/03/12
                                    epj_send_buf[ith].push_back(epj_sorted_[adr_ptcl_tmp++]);
                                    const F64vec pos_new = epj_send_buf[ith].back().getPos() - shift_image_domain[j];
                                    epj_send_buf[ith].back().setPos(pos_new);
                                }
                            }
                            else{
                                //if( pos_root_cell_.notOverlapped(tc_loc_[0].mom_.getPos()-shift_image_domain[j]) ) continue; // added by M.I. 2016/03/12
				if( pos_root_cell_.notContained(tc_loc_[0].mom_.getPos()-shift_image_domain[j]) ) continue; // added by M.I. 2016/03/12
                                spj_send_buf[ith].increaseSize();
                                spj_send_buf[ith].back().copyFromMoment(tc_loc_[0].mom_);
                                const F64vec pos_new = spj_send_buf[ith].back().getPos() - shift_image_domain[j];
                                spj_send_buf[ith].back().setPos(pos_new);
                            }
                        }
                    }
                }
                n_ep_send_[i] = epj_send_buf[ith].size() - n_epj_cum;
                n_sp_send_[i] = spj_send_buf[ith].size() - n_spj_cum;
                n_epj_cum = epj_send_buf[ith].size();
                n_spj_cum = spj_send_buf[ith].size();
                id_proc_send_[ith][n_proc_cum++] = i;
            } // end of for loop
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif	    
            {
                n_ep_send_disp_[0] = n_sp_send_disp_[0] = 0;
                for(S32 ib=0; ib<n_proc; ib++){
                    n_ep_send_disp_[ib+1] = n_ep_send_disp_[ib] + n_ep_send_[ib];
                    n_sp_send_disp_[ib+1] = n_sp_send_disp_[ib] + n_sp_send_[ib];
                }
                epj_send_.resizeNoInitialize( n_ep_send_disp_[n_proc] );
                spj_send_.resizeNoInitialize( n_sp_send_disp_[n_proc] );
            }
            S32 n_ep_cnt = 0;
            S32 n_sp_cnt = 0;
            for(S32 ib=0; ib<n_proc_cum; ib++){
                const S32 id = id_proc_send_[ith][ib];
                S32 adr_ep_tmp = n_ep_send_disp_[id];
                const S32 n_ep_tmp = n_ep_send_[id];
                for(int ip=0; ip<n_ep_tmp; ip++){
                    epj_send_[adr_ep_tmp++] = epj_send_buf[ith][n_ep_cnt++];
                }
                S32 adr_sp_tmp = n_sp_send_disp_[id];
                const S32 n_sp_tmp = n_sp_send_[id];
                for(int ip=0; ip<n_sp_tmp; ip++){
                    spj_send_[adr_sp_tmp++] = spj_send_buf[ith][n_sp_cnt++];
                }
            }
        } // omp parallel scope
        time_profile_.make_LET_1st += GetWtime() - time_offset;
        time_offset = GetWtime();

        Comm::allToAll(n_ep_send_, 1, n_ep_recv_); // TEST
        Comm::allToAll(n_sp_send_, 1, n_sp_recv_); // TEST
        n_ep_recv_disp_[0] = n_sp_recv_disp_[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_disp_[i+1] = n_ep_recv_disp_[i] + n_ep_recv_[i];
            n_sp_recv_disp_[i+1] = n_sp_recv_disp_[i] + n_sp_recv_[i];
        }
        epj_recv_.resizeNoInitialize( n_ep_recv_disp_[n_proc] );
        spj_recv_.resizeNoInitialize( n_sp_recv_disp_[n_proc] );
        Comm::allToAllV(epj_send_.getPointer(), n_ep_send_, n_ep_send_disp_,
                        epj_recv_.getPointer(), n_ep_recv_, n_ep_recv_disp_); // TEST
        Comm::allToAllV(spj_send_.getPointer(), n_sp_send_, n_sp_send_disp_, 
                        spj_recv_.getPointer(), n_sp_recv_, n_sp_recv_disp_); // TEST

        time_profile_.exchange_LET_1st += GetWtime() - time_offset;

        n_let_ep_send_1st_ += (CountT)epj_send_.size();
        n_let_ep_recv_1st_ += (CountT)epj_recv_.size();
        n_let_sp_send_1st_ += (CountT)spj_send_.size();
        n_let_sp_recv_1st_ += (CountT)spj_recv_.size();
        //time_profile_.exchange_LET_tot += time_profile_.make_LET_1st + time_profile_.exchange_LET_1st;
    }



    // FOR P^3T
    // original version
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchLongScatter, const DomainInfo & dinfo){
        F64 time_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	
#pragma omp parallel
#endif	
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_proc_cum = 0;
            id_ep_send_buf_[ith].reserve(1000);
            id_sp_send_buf_[ith].reserve(1000);
            id_ep_send_buf_[ith].resizeNoInitialize(0);
            id_sp_send_buf_[ith].resizeNoInitialize(0);
//#ifdef DEBUG_1023
#ifdef PARTICLE_SIMULATOR_EXCHANGE_LET_ALL
	    /// NEW MUST REMOVED
	    const F64 r_crit_sq = 1e10;
#else //PARTICLE_SIMULATOR_EXCHANGE_LET_ALL
            const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
#endif //PARTICLE_SIMULATOR_EXCHANGE_LET_ALL
            const S32 adr_tc_tmp = tc_loc_[0].adr_tc_;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 ib=0; ib<n_proc; ib++){
                n_ep_send_[ib] = n_sp_send_[ib] = n_ep_recv_[ib] = n_sp_recv_[ib] = 0;
                if(my_rank == ib) continue;
                const F64ort pos_target_domain = dinfo.getPosDomain(ib);
                const S32 n_ep_cum = id_ep_send_buf_[ith].size();
                const S32 n_sp_cum = id_sp_send_buf_[ith].size();
                if(!tc_loc_[0].isLeaf(n_leaf_limit_)){
                    SearchSendParticleLongScatter<TreeCell<Tmomloc>, Tepj>
                        (tc_loc_,               adr_tc_tmp,    epj_sorted_,
                         id_ep_send_buf_[ith],  id_sp_send_buf_[ith],
                         pos_target_domain,  r_crit_sq,   n_leaf_limit_);
                }
                else{
                    F64vec pos_tmp = tc_loc_[0].mom_.getPos();
                    if( (pos_target_domain.getDistanceMinSQ(pos_tmp) <= r_crit_sq * 4.0) ){
                        const S32 n_tmp = tc_loc_[0].n_ptcl_;
                        S32 adr_tmp = tc_loc_[0].adr_ptcl_;
                        for(S32 ip=0; ip<n_tmp; ip++){
                            id_ep_send_buf_[ith].push_back(adr_tmp++);
                        }
                    }
                    else{
                        id_sp_send_buf_[ith].push_back(0); // set root
                    }
                }
                n_ep_send_[ib] = id_ep_send_buf_[ith].size() - n_ep_cum;
                n_sp_send_[ib] = id_sp_send_buf_[ith].size() - n_sp_cum;
                id_proc_send_[ith][n_proc_cum++] = ib;
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp single
#endif
            {
                n_ep_send_disp_[0] = n_sp_send_disp_[0] = 0;
                for(S32 ib=0; ib<n_proc; ib++){
                    n_ep_send_disp_[ib+1] = n_ep_send_disp_[ib] + n_ep_send_[ib];
                    n_sp_send_disp_[ib+1] = n_sp_send_disp_[ib] + n_sp_send_[ib];
                }
                epj_send_.resizeNoInitialize( n_ep_send_disp_[n_proc] );
                spj_send_.resizeNoInitialize( n_sp_send_disp_[n_proc] );
            }
            S32 n_ep_cnt = 0;
            S32 n_sp_cnt = 0;
            for(S32 ib=0; ib<n_proc_cum; ib++){
                const S32 id = id_proc_send_[ith][ib];
                S32 adr_ep_tmp = n_ep_send_disp_[id];
                const S32 n_ep_tmp = n_ep_send_[id];
                for(int ip=0; ip<n_ep_tmp; ip++){
                    const S32 id_ep = id_ep_send_buf_[ith][n_ep_cnt++];
                    epj_send_[adr_ep_tmp++] = epj_sorted_[id_ep];
                }
                S32 adr_sp_tmp = n_sp_send_disp_[id];
                const S32 n_sp_tmp = n_sp_send_[id];
                for(int ip=0; ip<n_sp_tmp; ip++){
                    const S32 id_sp = id_sp_send_buf_[ith][n_sp_cnt++];
                    spj_send_[adr_sp_tmp++].copyFromMoment(tc_loc_[id_sp].mom_);
                }
            }
        } // omp parallel scope

        time_profile_.make_LET_1st += GetWtime() - time_offset;
        time_offset = GetWtime();

        for(S32 i=0; i<n_proc; i++){
            n_ep_sp_send_[2*i] = n_ep_send_[i];
            n_ep_sp_send_[2*i+1] = n_sp_send_[i];
        }

        F64 wtime_offset_tmp = GetWtime();
#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<S32, 2> comm_a2a_2d;
        comm_a2a_2d.execute(n_ep_sp_send_, 2, n_ep_sp_recv_);
#else
        Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_); // TEST
#endif //FAST_ALL_TO_ALL_FOR_K
        time_profile_.exchange_LET_1st__a2a_n += GetWtime() - wtime_offset_tmp;

        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_[i] = n_ep_sp_recv_[2*i];
            n_sp_recv_[i] = n_ep_sp_recv_[2*i+1];
        }
        n_ep_recv_disp_[0] = n_sp_recv_disp_[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_disp_[i+1] = n_ep_recv_disp_[i] + n_ep_recv_[i];
            n_sp_recv_disp_[i+1] = n_sp_recv_disp_[i] + n_sp_recv_[i];
        }
        epj_recv_.resizeNoInitialize( n_ep_recv_disp_[n_proc] );
        spj_recv_.resizeNoInitialize( n_sp_recv_disp_[n_proc] );


#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<Tepj, 2> comm_a2a_epj_2d;
        static CommForAllToAll<Tspj, 2> comm_a2a_spj_2d;

        wtime_offset_tmp = GetWtime();
        comm_a2a_epj_2d.executeV(epj_send_, epj_recv_, n_ep_send_, n_ep_recv_);
        time_profile_.exchange_LET_1st__a2a_ep += GetWtime() - wtime_offset_tmp;

        wtime_offset_tmp = GetWtime();
        comm_a2a_spj_2d.executeV(spj_send_, spj_recv_, n_sp_send_, n_sp_recv_);
        time_profile_.exchange_LET_1st__a2a_sp += GetWtime() - wtime_offset_tmp;
#else
        wtime_offset_tmp = GetWtime();
        Comm::allToAllV(epj_send_.getPointer(), n_ep_send_, n_ep_send_disp_,
                        epj_recv_.getPointer(), n_ep_recv_, n_ep_recv_disp_);
        time_profile_.exchange_LET_1st__a2a_ep += GetWtime() - wtime_offset_tmp;

        wtime_offset_tmp = GetWtime();
        Comm::allToAllV(spj_send_.getPointer(), n_sp_send_, n_sp_send_disp_,
                        spj_recv_.getPointer(), n_sp_recv_, n_sp_recv_disp_);
        time_profile_.exchange_LET_1st__a2a_sp += GetWtime() - wtime_offset_tmp;
#endif //FAST_ALL_TO_ALL_FOR_K
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;


        n_let_ep_send_1st_ += (CountT)epj_send_.size();
        n_let_ep_recv_1st_ += (CountT)epj_recv_.size();
        n_let_sp_send_1st_ += (CountT)spj_send_.size();
        n_let_sp_recv_1st_ += (CountT)spj_recv_.size();
        //time_profile_.exchange_LET_tot += time_profile_.make_LET_1st + time_profile_.exchange_LET_1st;
    }









#if 0
    // FOR P^3T + PM
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchLongCutoffScatter, const DomainInfo & dinfo){
        const F64 time_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            S32 n_proc_cum = 0;
            id_ep_send_buf_[ith].reserve(1000);
            id_sp_send_buf_[ith].reserve(1000);
            shift_image_domain_[ith].reserve(125);
            id_ep_send_buf_[ith].resizeNoInitialize(0);
            id_sp_send_buf_[ith].resizeNoInitialize(0);
            shift_image_domain_[ith].resizeNoInitialize(0);
            const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
            const S32 adr_tc_tmp = tc_loc_[0].adr_tc_;
            bool pa[DIMENSION];
            dinfo.getPeriodicAxis(pa);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 ib=0; ib<n_proc; ib++){
                n_ep_send_[ib] = n_sp_send_[ib] = n_ep_recv_[ib] = n_sp_recv_[ib] = 0;
                shift_image_domain_.clearSize();
                if(dinfo.getBoundaryCondition() == BOUNDARY_CONDITION_OPEN){
                    shift_image_domain_.push_back( F64vec(0.0) );
                }
                else{
                    CalcNumberAndShiftOfImageDomain
                        (shift_image_domain_, dinfo.getPosRootDomain().getFullLength(),
                         outer_boundary_of_my_tree, dinfo.getPosDomain(i), periodic_axis);
                }
                S32 n_image = shift_image_domain.size();
                for(S32 j = 0; j < n_image; j++) {
                    if(my_rank == i && j == 0) continue;
                    const F64ort pos_target_domain =dinfo.getPosDomain(i).shift(shift_image_domain[j]);
                    const F64ort cell_box = pos_root_cell_;
                    if( !tc_loc_[0].isLeaf(n_leaf_limit_) ){
                        SearchSendParticleLongCutoffScatter();
                    }
                    else{}
                }
                n_ep_send_[i] = epj_send_buf[ith].size() - n_epj_cum;
                n_sp_send_[i] = spj_send_buf[ith].size() - n_spj_cum;
                n_epj_cum = epj_send_buf[ith].size();
                n_spj_cum = spj_send_buf[ith].size();
                id_proc_send_[ith][n_proc_cum++] = i;
            } // end of for loop
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp single
#endif
            {
                n_ep_send_disp_[0] = n_sp_send_disp_[0] = 0;
                for(S32 ib=0; ib<n_proc; ib++){
                    n_ep_send_disp_[ib+1] = n_ep_send_disp_[ib] + n_ep_send_[ib];
                    n_sp_send_disp_[ib+1] = n_sp_send_disp_[ib] + n_sp_send_[ib];
                }
                epj_send_.resizeNoInitialize( n_ep_send_disp_[n_proc] );
                spj_send_.resizeNoInitialize( n_sp_send_disp_[n_proc] );
            }
            S32 n_ep_cnt = 0;
            S32 n_sp_cnt = 0;
            for(S32 ib=0; ib<n_proc_cum; ib++){
                const S32 id = id_proc_send_[ith][ib];
                S32 adr_ep_tmp = n_ep_send_disp_[id];
                const S32 n_ep_tmp = n_ep_send_[id];
                for(int ip=0; ip<n_ep_tmp; ip++){
                    const S32 id_ep = id_ep_send_buf_[ith][n_ep_cnt++];
                    epj_send_[adr_ep_tmp++] = epj_sorted_[id_ep];
                }
                S32 adr_sp_tmp = n_sp_send_disp_[id];
                const S32 n_sp_tmp = n_sp_send_[id];
                for(int ip=0; ip<n_sp_tmp; ip++){
                    const S32 id_sp = id_sp_send_buf_[ith][n_sp_cnt++];
                    spj_send_[adr_sp_tmp++].copyFromMoment(tc_loc_[id_sp].mom_);
                }
            }
        } // omp parallel scope

        time_profile_.make_LET_1st += GetWtime() - time_offset;
        time_offset = GetWtime();

        for(S32 i=0; i<n_proc; i++){
            n_ep_sp_send_[2*i] = n_ep_send_[i];
            n_ep_sp_send_[2*i+1] = n_sp_send_[i];
        }
#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<S32, 2> comm_a2a_2d;
        comm_a2a_2d.execute(n_ep_sp_send_, 2, n_ep_sp_recv_);
#else
        Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_); // TEST
#endif //FAST_ALL_TO_ALL_FOR_K
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_[i] = n_ep_sp_recv_[2*i];
            n_sp_recv_[i] = n_ep_sp_recv_[2*i+1];
        }
        n_ep_recv_disp_[0] = n_sp_recv_disp_[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_disp_[i+1] = n_ep_recv_disp_[i] + n_ep_recv_[i];
            n_sp_recv_disp_[i+1] = n_sp_recv_disp_[i] + n_sp_recv_[i];
        }
        epj_recv_.resizeNoInitialize( n_ep_recv_disp_[n_proc] );
        spj_recv_.resizeNoInitialize( n_sp_recv_disp_[n_proc] );
#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<Tepj, 2> comm_a2a_epj_2d;
        static CommForAllToAll<Tspj, 2> comm_a2a_spj_2d;
        comm_a2a_epj_2d.executeV(epj_send_, epj_recv_, n_ep_send_, n_ep_recv_);
        comm_a2a_spj_2d.executeV(spj_send_, spj_recv_, n_sp_send_, n_sp_recv_);
#else
        Comm::allToAllV(epj_send_.getPointer(), n_ep_send_, n_ep_send_disp_,
                        epj_recv_.getPointer(), n_ep_recv_, n_ep_recv_disp_);
        Comm::allToAllV(spj_send_.getPointer(), n_sp_send_, n_sp_send_disp_,
                        spj_recv_.getPointer(), n_sp_recv_, n_sp_recv_disp_);
#endif //FAST_ALL_TO_ALL_FOR_K
        time_profile_.exchange_LET_1st += GetWtime() - time_offset;

        n_let_ep_send_1st_ += (CountT)epj_send_.size();
        n_let_ep_recv_1st_ += (CountT)epj_recv_.size();
        n_let_sp_send_1st_ += (CountT)spj_send_.size();
        n_let_sp_recv_1st_ += (CountT)spj_recv_.size();
        //time_profile_.exchange_LET_tot += time_profile_.make_LET_1st + time_profile_.exchange_LET_1st;
    }
#endif



    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchShortScatter, const DomainInfo & dinfo){
        scatterEP(n_ep_send_,  n_ep_send_disp_,
                  n_ep_recv_,  n_ep_recv_disp_, 
                  epj_send_,   epj_recv_,
                  ep_send_buf_for_scatter_,
                  epj_sorted_, dinfo);	
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"epj_send_.size()="<<epj_send_.size()<<" spj_send_.size()="<<spj_send_.size()<<std::endl;
        std::cout<<"epj_recv_.size()="<<epj_recv_.size()<<" spj_recv_.size()="<<spj_recv_.size()<<std::endl;
#endif
        n_let_ep_send_1st_ += (CountT)epj_send_.size();
        n_let_ep_recv_1st_ += (CountT)epj_recv_.size();
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchShortGather, const DomainInfo & dinfo){
        exchangeLocalEssentialTreeGatherImpl(typename HasRSearch<Tepj>::type(), dinfo);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeGatherImpl(TagRSearch, const DomainInfo & dinfo){
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        std::cout<<"EPJ has RSearch"<<std::endl;
#endif
        exchangeLocalEssentialTreeImpl(TagSearchShortSymmetry(), dinfo);
    }

    /////////////////////////////
    /// link cell global tree ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    linkCellGlobalTreeOnly(){
        const F64 time_offset = GetWtime();
        //#ifdef SUNWAY
#if 1
        LinkCellSunWay(tc_glb_, adr_tc_level_partition_glb_, tp_glb_.getPointer(), lev_max_glb_, n_glb_tot_, n_leaf_limit_);
#else
        LinkCell(tc_glb_, adr_tc_level_partition_glb_, tp_glb_.getPointer(), lev_max_glb_, n_glb_tot_, n_leaf_limit_);
        //LinkCell(tc_glb_, adr_tc_level_partition_, tp_glb_.getPointer(), lev_max_, n_glb_tot_, n_leaf_limit_);
#endif
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
        calcMomentGlobalTreeOnlyImpl(typename TSM::force_type());
        time_profile_.calc_moment_global_tree += GetWtime() - time_offset;
        time_profile_.make_global_tree_tot = time_profile_.calc_moment_global_tree + time_profile_.make_global_tree;
    }
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
        CalcMoment(adr_tc_level_partition_glb_, tc_glb_.getPointer(), epj_sorted_.getPointer(), lev_max_glb_, n_leaf_limit_);
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
        
#ifdef USE_MEMORY_POOL
        MakeIPGroupUseGLBTreeLong(ipg_, etc_loc_, tc_glb_, epi_sorted_, 0, 0, n_group_limit_, n_leaf_limit_); // NEW
#else // USE_MEMORY_POOL
    #if 1
        MakeIPGroupUseGLBTreeLong(ipg_, tc_loc_, tc_glb_, epi_sorted_, 0, 0, n_group_limit_, n_leaf_limit_); // NEW
    #else
        MakeIPGroupLong(ipg_, tc_loc_, epi_sorted_, 0, n_group_limit_);
    #endif
#endif        
	const S32 n_ipg = ipg_.size();
    #ifdef SUNWAY
        unsigned long args[3];
        args[0] = (unsigned long) n_ipg;
        args[1] = (unsigned long) ipg_.getPointer();
        args[2] = (unsigned long) epi_sorted_.getPointer(); 
        /*
#ifdef REMOVE_VERTEX
        args[2] = (unsigned long) epi_sorted_.getPointer();
#else
        args[2] = (unsigned long) epj_sorted_loc_.getPointer();
#endif
        */
        __real_athread_spawn((void *)slave_GetMinIpgBox, args); // inside this, x,y,z->r,phi,z to get ipg box
        athread_join();
    #else // SUNWAY
	for(S32 i=0; i<n_ipg; i++){
	    const S32 n = ipg_[i].n_ptcl_;
	    const S32 adr = ipg_[i].adr_ptcl_;
        #ifdef PHI_R_TREE
            //ipg_[i].vertex_ = GetMinBoxSingleThreadForTree(epj_sorted_loc_.getPointer(adr), n);
            ipg_[i].vertex_ = GetMinBoxSingleThreadForTree(epi_sorted_.getPointer(adr), n);
        #else // PHI_R_TREE
            ipg_[i].vertex_ = GetMinBoxSingleThread(epi_sorted_.getPointer(adr), n);
        #endif // PHI_R_TREE
        }
    #endif // SUNWAY
    #ifdef DEBUG_PRINT_MAKE_IPG
        Comm::barrier();
        if(Comm::getRank() == 0){
            for(S32 i=0; i<n_ipg; i++){
                std::cerr<<"ipg_["<<i<<"].vertex_"<<ipg_[i].vertex_<<std::endl;
            }
        }
        for(S32 i=0; i<n_ipg; i++){
            const S32 n = ipg_[i].n_ptcl_;
            const S32 adr = ipg_[i].adr_ptcl_;
            for(S32 j=adr; j<adr+n; j++){
                const F64 pos_x   = epi_sorted_[j].pos.x;
                const F64 pos_y   = epi_sorted_[j].pos.y;
                const F64 pos_z   = epi_sorted_[j].pos.z;
                F64 pos_phi = atan2(pos_y, pos_x);
                if(pos_phi < 0.0) pos_phi += 8.0 * atan(1.0);
                F64 pos_r   = sqrt(pos_x*pos_x + pos_y*pos_y);
                F64vec pos_tmp = F64vec(pos_phi, pos_r, pos_z);
                if( pos_tmp.x<(ipg_[i].vertex_.low_.x-1e-10) ||  pos_tmp.x>(ipg_[i].vertex_.high_.x+1e-10)
                    || pos_tmp.y<(ipg_[i].vertex_.low_.y-1e-10) ||  pos_tmp.y>(ipg_[i].vertex_.high_.y+1e-10)
                    || pos_tmp.z<(ipg_[i].vertex_.low_.z-1e-10) ||  pos_tmp.z>(ipg_[i].vertex_.high_.z+1e-10) ){
                    std::cerr<<"i= "<<i
                             <<" j= "<<j
                             <<" ipg_[i].vertex_= "<<ipg_[i].vertex_
                             <<" pos_tmp= "<<pos_tmp
                             <<std::endl;
                    assert(0);
                }
            }
        }
        Comm::barrier();
    #endif
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeIPGroupImpl(TagForceShort){
        MakeIPGroupShort(ipg_, tc_loc_, epi_sorted_, 0, n_group_limit_);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionList(const S32 adr_ipg, const bool clear){
        makeInteractionListImpl(typename TSM::search_type(), adr_ipg, clear);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchLong, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_;
        if(clear){
            epj_for_force_[ith].clearSize();
            spj_for_force_[ith].clearSize();
        }

        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
        if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
            MakeInteractionListLongEPSP
                (tc_glb_, tc_glb_[0].adr_tc_, 
                 tp_glb_, epj_sorted_, 
                 epj_for_force_[ith],
                 spj_sorted_, spj_for_force_[ith],
                 pos_target_box, r_crit_sq, n_leaf_limit_);
        }
        else{
            const F64vec pos_tmp = tc_glb_[0].mom_.getPos();
            if( pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq*4.0 ){
                const S32 n_tmp = tc_glb_[0].n_ptcl_;
                S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                spj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                for(S32 ip=0; ip<n_tmp; ip++){
                    if( GetMSB(tp_glb_[adr_ptcl_tmp].adr_ptcl_) == 0){
                        epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp++]);
                    }
                    else{
                        spj_for_force_[ith].pushBackNoCheck(spj_sorted_[adr_ptcl_tmp++]);
                    }
                }
            }
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchLongCutoff, 
                            const S32 adr_ipg, 
                            const bool clear,
                            const bool make_id_list){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_;
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
        const F64 r_cut_sq  = epj_sorted_[0].getRSearch() * epj_sorted_[0].getRSearch();
        if(clear){
            epj_for_force_[ith].clearSize();
            spj_for_force_[ith].clearSize();
        }
        const F64ort cell_box = pos_root_cell_;
        if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
            MakeInteractionListLongCutoffEPSP
                (tc_glb_, tc_glb_[0].adr_tc_, tp_glb_, 
                 epj_sorted_, epj_for_force_[ith],
                 spj_sorted_, spj_for_force_[ith],
                 cell_box,
                 pos_target_box, r_crit_sq, r_cut_sq, n_leaf_limit_); 
        }
        else{
            if(pos_target_box.getDistanceMinSQ(cell_box) <= r_cut_sq){
                const F64vec pos_tmp = tc_glb_[0].mom_.getPos();
                if( pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq * 4.0){
                    const S32 n_tmp = tc_glb_[0].n_ptcl_;
                    S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                    epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                    spj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                    for(S32 ip=0; ip<n_tmp; ip++){
                        if( GetMSB(tp_glb_[adr_ptcl_tmp].adr_ptcl_) == 0){
                            epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp++]);
                        }
                        else{
                            spj_for_force_[ith].pushBackNoCheck(spj_sorted_[adr_ptcl_tmp++]);
                        }
                    }
                }
                else{
                    spj_for_force_[ith].increaseSize();
                    spj_for_force_[ith].back().copyFromMoment(tc_glb_[0].mom_);
                }
            }
        }
    }


#if 1
    // FOR P^3T
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchLongScatter, 
                            const S32 adr_ipg, 
                            const bool clear,
                            const bool make_id_list){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_;
#ifdef PARTICLE_SIMULATOR_INTERACTION_LIST_ALL
        const F64 r_crit_sq = 1e10;
#else //PARTICLE_SIMULATOR_INTERACTION_LIST_ALL
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
#endif //PARTICLE_SIMULATOR_INTERACTION_LIST_ALL
	if(clear){	
	    epj_for_force_[ith].clearSize();
	    spj_for_force_[ith].clearSize();
	}
        if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
            MakeInteractionListLongScatterEPSP
                (tc_glb_, tc_glb_[0].adr_tc_, 
                 tp_glb_, epj_sorted_, 
                 epj_for_force_[ith],
                 spj_sorted_, spj_for_force_[ith],
                 pos_target_box, r_crit_sq, n_leaf_limit_);
        }
        else{
            const F64vec pos_tmp = tc_glb_[0].mom_.getPos();
            if( pos_target_box.getDistanceMinSQ(pos_tmp) <= r_crit_sq*4.0 ){
                const S32 n_tmp = tc_glb_[0].n_ptcl_;
                S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                spj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                for(S32 ip=0; ip<n_tmp; ip++){
                    if( GetMSB(tp_glb_[adr_ptcl_tmp].adr_ptcl_) == 0){
                        epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp++]);
                    }
                    else{
                        spj_for_force_[ith].pushBackNoCheck(spj_sorted_[adr_ptcl_tmp++]);
                    }
                }
            }
        }
    }
#endif


    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchShortScatter, 
                            const S32 adr_ipg, 
                            const bool clear,
                            const bool make_id_list){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_;
	if(clear){	
	    epj_for_force_[ith].clearSize();
	}
#if 0
        const S32 adr_root_cell = 0;
        n_cell_open_[ith] += MakeListUsingOuterBoundaryIteration
            (tc_glb_.getPointer(),     adr_root_cell,
             epj_sorted_.getPointer(), epj_for_force_[ith], 
             pos_target_box,           n_leaf_limit_);
#else
        if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
            MakeListUsingOuterBoundary
                (tc_glb_.getPointer(),     tc_glb_[0].adr_tc_,
                 epj_sorted_.getPointer(), epj_for_force_[ith], 
                 pos_target_box,   n_leaf_limit_);
        }
        else{
            //if( pos_target_box.overlapped( tc_glb_[0].mom_.getVertexOut()) ){
	    if( pos_target_box.contained( tc_glb_[0].mom_.getVertexOut()) ){
                S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                const S32 n_tmp = tc_glb_[0].n_ptcl_;
                //epj_for_force_[ith].reserve( epj_for_force_[ith].size() + n_tmp );
                epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
                    const F64vec pos_tmp = epj_sorted_[adr_ptcl_tmp].getPos();
                    const F64 size_tmp = epj_sorted_[adr_ptcl_tmp].getRSearch();
                    const F64 dis_sq_tmp = pos_target_box.getDistanceMinSQ(pos_tmp);
                    if(dis_sq_tmp > size_tmp*size_tmp) continue;
                    //epj_for_force_[ith].increaseSize();
                    //epj_for_force_[ith].back() = epj_sorted_[adr_ptcl_tmp];
                    //epj_for_force_[ith].push_back(epj_sorted_[adr_ptcl_tmp]);
                    epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp]);
                    const F64vec pos_new = epj_for_force_[ith].back().getPos();
                    epj_for_force_[ith].back().setPos(pos_new);
                }
            }
        }
#endif
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchShortGather,
                            const S32 adr_ipg, 
                            const bool clear,
                            const bool make_id_list){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_;
	if(clear){
	    epj_for_force_[ith].clearSize();
	}
#if 1
        const S32 adr_root_cell = 0;
        n_cell_open_[ith] += MakeListUsingInnerBoundaryIteration
            (tc_glb_.getPointer(),     adr_root_cell,
             epj_sorted_.getPointer(), epj_for_force_[ith],
             pos_target_box,           n_leaf_limit_);
#else
        if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
            MakeListUsingInnerBoundary
                (tc_glb_.getPointer(),     tc_glb_[0].adr_tc_,
                 epj_sorted_.getPointer(), epj_for_force_[ith],
                 pos_target_box,                 n_leaf_limit_);
        }
        else{
            //if( pos_target_box.overlapped( tc_glb_[0].mom_.getVertexIn()) ){
	    if( pos_target_box.contained( tc_glb_[0].mom_.getVertexIn()) ){
                S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                const S32 n_tmp = tc_glb_[0].n_ptcl_;
                //epj_for_force_[ith].reserve( epj_for_force_[ith].size() + n_tmp );
                epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
                    const F64vec pos_tmp = epj_sorted_[adr_ptcl_tmp].getPos();
                    //if( pos_target_box.overlapped( pos_tmp) ){
		    if( pos_target_box.contained( pos_tmp) ){
                        //epj_for_force_[ith].push_back(epj_sorted_[adr_ptcl_tmp]);
                        epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp]);
                        const F64vec pos_new = epj_for_force_[ith].back().getPos();
                        epj_for_force_[ith].back().setPos(pos_new);
                    }
                }
            }
        }
#endif
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchShortSymmetry, 
                            const S32 adr_ipg, 
                            const bool clear,
                            const bool make_id_list){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box_out = (ipg_[adr_ipg]).vertex_;
        const F64ort pos_target_box_in = (ipg_[adr_ipg]).vertex_in;
	if(clear){
	    epj_for_force_[ith].clearSize();
	}
#if 1
        const S32 adr_root_cell = 0;
        n_cell_open_[ith] += MakeListUsingOuterBoundaryAndInnerBoundaryIteration
            (tc_glb_.getPointer(),     adr_root_cell,
             epj_sorted_.getPointer(), epj_for_force_[ith],
             pos_target_box_out,       pos_target_box_in, 
             n_leaf_limit_);
#else
        if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
            MakeListUsingOuterBoundaryAndInnerBoundary
                (tc_glb_.getPointer(),     tc_glb_[0].adr_tc_,
                 epj_sorted_.getPointer(), epj_for_force_[ith],
                 pos_target_box_out,       pos_target_box_in,
                 n_leaf_limit_);
        }
        else{
	    /*
            if( pos_target_box_out.overlapped(tc_glb_[0].mom_.getVertexIn()) 
                || pos_target_box_in.overlapped(tc_glb_[9].mom_.getVertexOut()) ){
	    */
            if( pos_target_box_out.contained(tc_glb_[0].mom_.getVertexIn()) 
                || pos_target_box_in.contained(tc_glb_[9].mom_.getVertexOut()) ){
                S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                const S32 n_tmp = tc_glb_[0].n_ptcl_;
                epj_for_force_[ith].reserveAtLeast( n_tmp );
                for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
                    const F64vec pos_tmp = epj_sorted_[adr_ptcl_tmp].getPos();
                    const F64 size_tmp = epj_sorted_[adr_ptcl_tmp].getRSearch();
                    const F64 dis_sq_tmp = pos_target_box_in.getDistanceMinSQ(pos_tmp);
                    //if( pos_target_box_out.notOverlapped(pos_tmp) && dis_sq_tmp > size_tmp*size_tmp) continue;
		    if( pos_target_box_out.notContained(pos_tmp) && dis_sq_tmp > size_tmp*size_tmp) continue;
                    epj_for_force_[ith].pushBackNoCheck( epj_sorted_[adr_ptcl_tmp] );
                    const F64vec pos_new = epj_for_force_[ith].back().getPos();
                    epj_for_force_[ith].back().setPos(pos_new);
                }
            }
        }
#endif
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
            for(S32 i=offset; i<n_tail; i++) force_sorted_[i].clear();
        }
        pfunc_ep_ep(epi_sorted_.getPointer(offset),     n_epi,
                    epj_for_force_[ith].getPointer(),   n_epj,
                    force_sorted_.getPointer(offset));
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
            for(S32 i=offset; i<n_tail; i++) force_sorted_[i].clear();
        }
        pfunc_ep_ep(epi_sorted_.getPointer(offset),     n_epi,
                    epj_for_force_[ith].getPointer(),   n_epj,
                    force_sorted_.getPointer(offset));
        pfunc_ep_sp(epi_sorted_.getPointer(offset),     n_epi,
                    spj_for_force_[ith].getPointer(),   n_spj,
                    force_sorted_.getPointer(offset));
    }


















    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    copyForceOriginalOrder(){
        const F64 wtime_offset = GetWtime();
        force_org_.resizeNoInitialize(n_loc_tot_);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	
#pragma omp parallel for
#endif	
        for(S32 i=0; i<n_loc_tot_; i++){
#ifdef REMOVE_TP_LOC
            const S32 adr = ClearMSB(tp_glb_[i].adr_ptcl_);
#else
            const S32 adr = ClearMSB(tp_loc_[i].adr_ptcl_);
#endif
            force_org_[adr] = force_sorted_[i];
        }
        time_profile_.copy_force_original_order += GetWtime() - wtime_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForce(Tfunc_ep_ep pfunc_ep_ep,
              const bool clear){
        F64 time_offset = GetWtime();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        force_org_.resizeNoInitialize(n_loc_tot_);
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
                calcForceOnly( pfunc_ep_ep, i, clear);
            }
            ni_ave_ = ni_tmp / n_ipg;
            nj_ave_ = nj_tmp / n_ipg;
            n_interaction_ep_ep_ += n_interaction_ep_ep_tmp;
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
            for(S32 i=1; i<Comm::getNumberOfThread(); i++) n_cell_open_[0] += n_cell_open_[i];
        }
        else{
            ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = 0;
            n_walk_local_ = n_interaction_ep_ep_local_ = 0;
        }
        //PROFILE::Stop(profile.calc_force);
        copyForceOriginalOrder();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"ipg_.size()="<<ipg_.size()<<std::endl;
        std::cout<<"force_sorted_.size()="<<force_sorted_.size()<<std::endl;
        std::cout<<"force_org_.size()="<<force_org_.size()<<std::endl;
#endif
        time_profile_.calc_force += GetWtime() - time_offset;
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceWalkOnly(Tfunc_ep_ep pfunc_ep_ep,
                      const bool clear){
        F64 time_offset = GetWtime();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        force_org_.resizeNoInitialize(n_loc_tot_);
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
            n_interaction_ep_ep_ += n_interaction_ep_ep_tmp;
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
            for(S32 i=1; i<Comm::getNumberOfThread(); i++) n_cell_open_[0] += n_cell_open_[i];
        }
        else{
            ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = 0;
            n_walk_local_ = n_interaction_ep_ep_local_ = 0;
        }
        //PROFILE::Stop(profile.calc_force);
        copyForceOriginalOrder();
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"ipg_.size()="<<ipg_.size()<<std::endl;
        std::cout<<"force_sorted_.size()="<<force_sorted_.size()<<std::endl;
        std::cout<<"force_org_.size()="<<force_org_.size()<<std::endl;
#endif
        time_profile_.calc_force += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep, class Tfunc_ep_sp>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForce(Tfunc_ep_ep pfunc_ep_ep,
              Tfunc_ep_sp pfunc_ep_sp,
              const bool clear){
        const F64 time_offset = GetWtime();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        force_org_.resizeNoInitialize(n_loc_tot_);
        const S64 n_ipg = ipg_.size();
        n_walk_local_ += n_ipg;
        S64 ni_tmp = 0;
        S64 nj_tmp = 0;
        S64 n_interaction_ep_ep_tmp = 0;
        S64 n_interaction_ep_sp_tmp = 0;
        for(S32 i=0; i<Comm::getNumberOfThread(); i++) n_cell_open_[i] = 0;
        //PROFILE::Start(profile.calc_force);
        if(n_ipg > 0){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : ni_tmp, nj_tmp, n_interaction_ep_ep_tmp, n_interaction_ep_sp_tmp)
#endif	    
            for(S32 i=0; i<n_ipg; i++){
                makeInteractionList(i);
                ni_tmp += ipg_[i].n_ptcl_;
                nj_tmp += epj_for_force_[Comm::getThreadNum()].size();
                nj_tmp += spj_for_force_[Comm::getThreadNum()].size();
                n_interaction_ep_ep_tmp += ipg_[i].n_ptcl_ * epj_for_force_[Comm::getThreadNum()].size();
                n_interaction_ep_sp_tmp += ipg_[i].n_ptcl_ * spj_for_force_[Comm::getThreadNum()].size();
                calcForceOnly( pfunc_ep_ep, pfunc_ep_sp, i, clear);
            }
            ni_ave_ = ni_tmp / n_ipg;
            nj_ave_ = nj_tmp / n_ipg;
            n_interaction_ep_ep_ = n_interaction_ep_ep_tmp;
            n_interaction_ep_sp_ = n_interaction_ep_sp_tmp;
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
            n_interaction_ep_sp_local_ += n_interaction_ep_sp_tmp;
            for(S32 i=1; i<Comm::getNumberOfThread(); i++) n_cell_open_[0] += n_cell_open_[i];
        }
        else{
            ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
        }
        //PROFILE::Stop(profile.calc_force);
        copyForceOriginalOrder();

#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"ipg_.size()="<<ipg_.size()<<std::endl;
        std::cout<<"force_sorted_.size()="<<force_sorted_.size()<<std::endl;
        std::cout<<"force_org_.size()="<<force_org_.size()<<std::endl;
#endif
        time_profile_.calc_force += GetWtime() - time_offset;
    }





    //////////////
    /// for multi walk
#if 0
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                            Tfunc_retrieve pfunc_retrieve,
                            const S32 tag_max,
                            const S32 n_walk_limit,
                            const bool clear){
        if(tag_max <= 0){
            PARTICLE_SIMULATOR_PRINT_ERROR("tag_max is illegal. In currente version, tag_max must be 1");
            Abort(-1);
        }
        S32 ret = 0;
        const F64 time_offset = GetWtime();
        ret = calcForceMultiWalkIndexImpl(typename TSM::force_type(),
                                          pfunc_dispatch,
                                          pfunc_retrieve,
                                          tag_max,
                                          n_walk_limit,
                                          clear);
        time_profile_.calc_force += GetWtime() - time_offset;
        return ret;
    }

#if 0    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkIndexImpl(TagForceLong,
                                Tfunc_dispatch pfunc_dispatch,
                                Tfunc_retrieve pfunc_retrieve,
                                const S32 tag_max,
                                const S32 n_walk_limit,
                                const bool clear){
        const F64 offset_core = GetWtime();
        static bool first = true;
        S32 ret = 0;
        S32 tag = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        force_org_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        n_walk_local_ += n_ipg;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        static std::vector<S32> walk_grp;
        static std::vector<S32> walk_grp_disp;
        walk_grp.resize(n_loop_max); // group of walk (walk_grp[i] < n_walk_limit)
        walk_grp_disp.resize(n_loop_max+1);
        static S32  * iw2ith;
        static S32  * iw2cnt;
        static S32  * n_epi_array;
        static S32  * n_epi_array_prev;
        static S32  * n_epj_array;
        static S32  * n_spj_array;
        static S32  ** n_epj_disp_thread; // [n_thread][n_walk]
        static S32  ** n_spj_disp_thread;// [n_thread][n_walk]
        static Tepi ** epi_array; // array of pointer *[n_walk]
        static S32  ** id_epj_array; // array of pointer *[n_walk]
        static S32  ** id_spj_array; // array of pointer *[n_walk]
        static Tforce ** force_array; // array of pointer *[n_walk]
        static Tforce ** force_prev_array; // array of pointer *[n_walk]
        static S32  * cnt_thread;
        static S32  * n_ep_cum_thread;
        static S32  * n_sp_cum_thread;
        static S64 * n_interaction_ep_ep_array;
        static S64 * n_interaction_ep_sp_array;
        if(first){
            iw2ith = new S32[n_walk_limit];
            iw2cnt = new S32[n_walk_limit];
            n_epi_array = new S32[n_walk_limit];
            n_epi_array_prev = new S32[n_walk_limit];
            n_epj_array = new S32[n_walk_limit];
            n_spj_array = new S32[n_walk_limit];
            n_epj_disp_thread = new S32*[n_thread];
            n_spj_disp_thread = new S32*[n_thread];
            epi_array        = new Tepi*[n_walk_limit];
            id_epj_array        = new S32*[n_walk_limit];
            id_spj_array        = new S32*[n_walk_limit];
            force_array      = new Tforce*[n_walk_limit];
            force_prev_array = new Tforce*[n_walk_limit];
            for(int i=0; i<n_thread; i++){
                n_epj_disp_thread[i] = new S32[n_walk_limit];
                n_spj_disp_thread[i] = new S32[n_walk_limit];
            }
            cnt_thread = new S32[n_thread];
            n_ep_cum_thread = new S32[n_thread];
            n_sp_cum_thread = new S32[n_thread];
            n_interaction_ep_ep_array = new S64[n_thread];
            n_interaction_ep_sp_array = new S64[n_thread];
            first = false;
        }
        walk_grp_disp[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            walk_grp[wg] = n_ipg / n_loop_max + 1;
            walk_grp_disp[wg+1] = walk_grp_disp[wg] + walk_grp[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
            walk_grp[wg] = n_ipg / n_loop_max;
            walk_grp_disp[wg+1] = walk_grp_disp[wg] + walk_grp[wg];
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_array[i] = n_interaction_ep_sp_array[i] = 0;
        }
        bool first_loop = true;
        S32 n_walk_prev = 0;
        if(n_ipg > 0){
            n_epj_array[0] = epj_sorted_.size();
            n_spj_array[0] = spj_sorted_.size();
            S32 tmp = pfunc_dispatch(tag, -1, 
                                     (const Tepi **)epi_array, n_epi_array, 
                                     (const Tepj **)epj_sorted_.getPointer(), (const S32**)id_epj_array, n_epj_array, 
                                     (const Tspj **)spj_sorted_.getPointer(), (const S32**)id_spj_array, n_spj_array);
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = walk_grp[wg];
                const S32 walk_grp_head = walk_grp_disp[wg];
                for(int i=0; i<n_thread; i++){
                    n_ep_cum_thread[i] = n_sp_cum_thread[i] = cnt_thread[i] = 0;
                    n_epj_disp_thread[i][0] = 0;
                    n_spj_disp_thread[i][0] = 0;
                    epj_for_force_[i].clearSize();
                    spj_for_force_[i].clearSize();
                }
                const F64 offset_calc_force__core__walk_tree = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4) 
#endif
                for(int iw=0; iw<n_walk; iw++){
                    const S32 id_ipg = walk_grp_head + iw;
                    const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_; 
                    const S32 ith = Comm::getThreadNum();
                    n_epi_array[iw] = ipg_[id_ipg].n_ptcl_;
                    epi_array[iw]   = epi_sorted_.getPointer(first_adr_ip);
                    force_array[iw] = force_sorted_.getPointer(first_adr_ip);
                    makeInteractionListIndexImpl(typename TSM::search_type(), id_ipg, false);
                    n_epj_array[iw] = epj_for_force_[ith].size() - n_ep_cum_thread[ith];
                    n_spj_array[iw] = spj_for_force_[ith].size() - n_sp_cum_thread[ith];
                    n_ep_cum_thread[ith] = epj_for_force_[ith].size();
                    n_sp_cum_thread[ith] = spj_for_force_[ith].size();
                    n_epj_disp_thread[ith][cnt_thread[ith]+1] = n_ep_cum_thread[ith];
                    n_spj_disp_thread[ith][cnt_thread[ith]+1] = n_sp_cum_thread[ith];
                    n_interaction_ep_ep_array[ith] += ((S64)n_epj_array[iw]*(S64)n_epi_array[iw]);
                    n_interaction_ep_sp_array[ith] += ((S64)n_spj_array[iw]*(S64)n_epi_array[iw]);
                    iw2ith[iw] = ith;
                    iw2cnt[iw] = cnt_thread[ith];
                    cnt_thread[ith]++;
                }
                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_array_prev, force_prev_array);
                } // retrieve
                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    S32 n_sp_head = n_spj_disp_thread[ith][cnt];
                    id_epj_array[iw] = id_epj_for_force_[ith].getPointer(n_ep_head);
                    id_spj_array[iw] = id_spj_for_force_[ith].getPointer(n_sp_head);
                }
                Tepj ** epj_array_tmp;
                Tspj ** spj_array_tmp;
                ret += pfunc_dispatch(tag, n_walk,   (const Tepi**)epi_array, n_epi_array, 
                                      (const Tepj **)epj_array_tmp, (const S32**)id_epj_array, n_epj_array, 
                                      (const Tspj **)spj_array_tmp, (const S32**)id_spj_array, n_spj_array);
                first_loop = false;

                for(int iw=0; iw<n_walk; iw++){
                    n_epi_array_prev[iw] = n_epi_array[iw];
                    force_prev_array[iw] = force_array[iw];
                }
                n_walk_prev = n_walk;
            } // end of walk group loop
            ret += pfunc_retrieve(tag, n_walk_prev, n_epi_array_prev, force_prev_array);
        } // if(n_ipg > 0)
        else{
            ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_array[i];
            n_interaction_ep_sp_local_ += n_interaction_ep_sp_array[i];
        }
        time_profile_.calc_force__core += GetWtime() - offset_core;
        const F64 offset_copy_original_order = GetWtime();
        copyForceOriginalOrder();
        time_profile_.calc_force__copy_original_order += GetWtime() - offset_copy_original_order;
        return ret;
    }
#endif
    
    // for short
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkIndexImpl(TagForceShort,
                                Tfunc_dispatch pfunc_dispatch,
                                Tfunc_retrieve pfunc_retrieve,
                                const S32 tag_max,
                                const S32 n_walk_limit,
                                const bool clear){
        /*
        const F64 offset_core = GetWtime();
        static bool first = true;
        S32 ret = 0;
        S32 tag = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        force_org_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        n_walk_local_ += n_ipg;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        static std::vector<S32> walk_grp;
        static std::vector<S32> walk_grp_disp;
        walk_grp.resize(n_loop_max); // group of walk (walk_grp[i] < n_walk_limit)
        walk_grp_disp.resize(n_loop_max+1);
        static S32  * iw2ith;
        static S32  * iw2cnt;
        static S32  * n_epi_array;
        static S32  * n_epi_array_prev;
        static S32  * n_epj_array;
        static S32 **  n_epj_disp_thread; // [n_thread][n_walk]
        static Tepi ** epi_array; // array of pointer *[n_walk]
        static S32  ** id_epj_array; // array of pointer *[n_walk]
        static Tforce ** force_array; // array of pointer *[n_walk]
        static Tforce ** force_prev_array; // array of pointer *[n_walk]
        static S32  * cnt_thread;
        static S32  * n_ep_cum_thread;
        static S64 * n_interaction_ep_ep_array;
        if(first){
            iw2ith = new S32[n_walk_limit];
            iw2cnt = new S32[n_walk_limit];
            n_epi_array = new S32[n_walk_limit];
            n_epi_array_prev = new S32[n_walk_limit];
            n_epj_array = new S32[n_walk_limit];
            n_epj_disp_thread = new S32*[n_thread];
            epi_array        = new Tepi*[n_walk_limit];
            id_epj_array        = new Tepj*[n_walk_limit];
            force_array      = new Tforce*[n_walk_limit];
            force_prev_array = new Tforce*[n_walk_limit];
            for(int i=0; i<n_thread; i++){
                n_epj_disp_thread[i] = new S32[n_walk_limit];
            }
            cnt_thread = new S32[n_thread];
            n_ep_cum_thread = new S32[n_thread];
            n_interaction_ep_ep_array = new S64[n_thread];
            first = false;
        }
        walk_grp_disp[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            walk_grp[wg] = n_ipg / n_loop_max + 1;
            walk_grp_disp[wg+1] = walk_grp_disp[wg] + walk_grp[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
            walk_grp[wg] = n_ipg / n_loop_max;
            walk_grp_disp[wg+1] = walk_grp_disp[wg] + walk_grp[wg];
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_array[i] = 0;
        }
        bool first_loop = true;
        S32 n_walk_prev = 0;
        if(n_ipg > 0){
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = walk_grp[wg];
                const S32 walk_grp_head = walk_grp_disp[wg];
                for(int i=0; i<n_thread; i++){
                    n_ep_cum_thread[i] = cnt_thread[i] = 0;
                    n_epj_disp_thread[i][0] = 0;
                    epj_for_force_[i].clearSize();
                }
                const F64 offset_calc_force__core__walk_tree = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4) 
#endif
                for(int iw=0; iw<n_walk; iw++){
                    const S32 id_ipg = walk_grp_head + iw;
                    const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_; 
                    const S32 ith = Comm::getThreadNum();
                    n_epi_array[iw] = ipg_[id_ipg].n_ptcl_;
                    epi_array[iw]   = epi_sorted_.getPointer(first_adr_ip);
                    force_array[iw] = force_sorted_.getPointer(first_adr_ip);
                    makeInteractionList(id_ipg, false);
                    n_epj_array[iw] = epj_for_force_[ith].size() - n_ep_cum_thread[ith];
                    n_ep_cum_thread[ith] = epj_for_force_[ith].size();
                    n_epj_disp_thread[ith][cnt_thread[ith]+1] = n_ep_cum_thread[ith];
                    n_interaction_ep_ep_array[ith] += (S64)n_epj_array[iw]*(S64)n_epi_array[iw];
                    iw2ith[iw] = ith;
                    iw2cnt[iw] = cnt_thread[ith];
                    cnt_thread[ith]++;
                }
                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_array_prev, force_prev_array);
v                } // retrieve
                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    id_epj_array[iw] = epj_for_force_[ith].getPointer(n_ep_head);
                }
                ret += pfunc_dispatch(tag, n_walk, (const Tepi**)epi_array, n_epi_array, (const S32**)id_epj_array, n_epj_array);
                first_loop = false;

                for(int iw=0; iw<n_walk; iw++){
                    n_epi_array_prev[iw] = n_epi_array[iw];
                    force_prev_array[iw] = force_array[iw];
                }
                n_walk_prev = n_walk;
            } // end of walk group loop
            ret += pfunc_retrieve(tag, n_walk_prev, n_epi_array_prev, force_prev_array);
        } // if(n_ipg > 0)
        else{
            ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_array[i];
        }
        time_profile_.calc_force__core += GetWtime() - offset_core;
        const F64 offset_copy_original_order = GetWtime();
        copyForceOriginalOrder();
        time_profile_.calc_force__copy_original_order += GetWtime() - offset_copy_original_order;
        return ret;
        */
    }
#endif //#if 0


    // no index version
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalk(Tfunc_dispatch pfunc_dispatch,
                       Tfunc_retrieve pfunc_retrieve,
                       const S32 tag_max,
                       const S32 n_walk_limit,
                       const bool clear){
        if(tag_max <= 0){
            PARTICLE_SIMULATOR_PRINT_ERROR("tag_max is illegal. In currente version, tag_max must be 1");
            Abort(-1);
        }
        S32 ret = 0;
        const F64 time_offset = GetWtime();
        ret = calcForceMultiWalkImpl(typename TSM::force_type(),
                                     pfunc_dispatch,
                                     pfunc_retrieve,
                                     tag_max,
                                     n_walk_limit,
                                     clear);
        time_profile_.calc_force += GetWtime() - time_offset;
        return ret;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkImpl(TagForceLong,
                           Tfunc_dispatch pfunc_dispatch,
                           Tfunc_retrieve pfunc_retrieve,
                           const S32 tag_max,
                           const S32 n_walk_limit,
                           const bool clear){
        const F64 offset_core = GetWtime();
        static bool first = true;
        S32 ret = 0;
        S32 tag = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        force_org_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        n_walk_local_ += n_ipg;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        static std::vector<S32> walk_grp;
        static std::vector<S32> walk_grp_disp;
        walk_grp.resize(n_loop_max); // group of walk (walk_grp[i] < n_walk_limit)
        walk_grp_disp.resize(n_loop_max+1);
        static S32  * iw2ith;
        static S32  * iw2cnt;
        static S32  * n_epi_array;
        static S32  * n_epi_array_prev;
        static S32  * n_epj_array;
        static S32  * n_spj_array;
        static S32 ** n_epj_disp_thread; // [n_thread][n_walk]
        static S32 ** n_spj_disp_thread;// [n_thread][n_walk]
        static Tepi ** epi_array; // array of pointer *[n_walk]
        static Tepj ** epj_array; // array of pointer *[n_walk]
        static Tspj ** spj_array; // array of pointer *[n_walk]
        static Tforce ** force_array; // array of pointer *[n_walk]
        static Tforce ** force_prev_array; // array of pointer *[n_walk]
        static S32  * cnt_thread;
        static S32  * n_ep_cum_thread;
        static S32  * n_sp_cum_thread;
        static S64 * n_interaction_ep_ep_array;
        static S64 * n_interaction_ep_sp_array;
        if(first){
            iw2ith = new S32[n_walk_limit];
            iw2cnt = new S32[n_walk_limit];
            n_epi_array = new S32[n_walk_limit];
            n_epi_array_prev = new S32[n_walk_limit];
            n_epj_array = new S32[n_walk_limit];
            n_spj_array = new S32[n_walk_limit];
            n_epj_disp_thread = new S32*[n_thread];
            n_spj_disp_thread = new S32*[n_thread];
            epi_array        = new Tepi*[n_walk_limit];
            epj_array        = new Tepj*[n_walk_limit];
            spj_array        = new Tspj*[n_walk_limit];
            force_array      = new Tforce*[n_walk_limit];
            force_prev_array = new Tforce*[n_walk_limit];
            for(int i=0; i<n_thread; i++){
                n_epj_disp_thread[i] = new S32[n_walk_limit];
                n_spj_disp_thread[i] = new S32[n_walk_limit];
            }
            cnt_thread = new S32[n_thread];
            n_ep_cum_thread = new S32[n_thread];
            n_sp_cum_thread = new S32[n_thread];
            n_interaction_ep_ep_array = new S64[n_thread];
            n_interaction_ep_sp_array = new S64[n_thread];
            first = false;
        }
        walk_grp_disp[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            walk_grp[wg] = n_ipg / n_loop_max + 1;
            walk_grp_disp[wg+1] = walk_grp_disp[wg] + walk_grp[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
            walk_grp[wg] = n_ipg / n_loop_max;
            walk_grp_disp[wg+1] = walk_grp_disp[wg] + walk_grp[wg];
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_array[i] = n_interaction_ep_sp_array[i] = 0;
        }
        bool first_loop = true;
        S32 n_walk_prev = 0;
        if(n_ipg > 0){
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = walk_grp[wg];
                const S32 walk_grp_head = walk_grp_disp[wg];
                for(int i=0; i<n_thread; i++){
                    n_ep_cum_thread[i] = n_sp_cum_thread[i] = cnt_thread[i] = 0;
                    n_epj_disp_thread[i][0] = 0;
                    n_spj_disp_thread[i][0] = 0;
                    epj_for_force_[i].clearSize();
                    spj_for_force_[i].clearSize();
                }
                const F64 offset_calc_force__core__walk_tree = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4) 
#endif
                for(int iw=0; iw<n_walk; iw++){
                    const S32 id_ipg = walk_grp_head + iw;
                    const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_; 
                    const S32 ith = Comm::getThreadNum();
                    n_epi_array[iw] = ipg_[id_ipg].n_ptcl_;
                    epi_array[iw]   = epi_sorted_.getPointer(first_adr_ip);
                    force_array[iw] = force_sorted_.getPointer(first_adr_ip);
                    makeInteractionList(id_ipg, false);
                    n_epj_array[iw] = epj_for_force_[ith].size() - n_ep_cum_thread[ith];
                    n_spj_array[iw] = spj_for_force_[ith].size() - n_sp_cum_thread[ith];
                    n_ep_cum_thread[ith] = epj_for_force_[ith].size();
                    n_sp_cum_thread[ith] = spj_for_force_[ith].size();
                    n_epj_disp_thread[ith][cnt_thread[ith]+1] = n_ep_cum_thread[ith];
                    n_spj_disp_thread[ith][cnt_thread[ith]+1] = n_sp_cum_thread[ith];
                    n_interaction_ep_ep_array[ith] += ((S64)n_epj_array[iw]*(S64)n_epi_array[iw]);
                    n_interaction_ep_sp_array[ith] += ((S64)n_spj_array[iw]*(S64)n_epi_array[iw]);
                    iw2ith[iw] = ith;
                    iw2cnt[iw] = cnt_thread[ith];
                    cnt_thread[ith]++;
                }
                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_array_prev, force_prev_array);
                } // retrieve
                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    S32 n_sp_head = n_spj_disp_thread[ith][cnt];
                    epj_array[iw] = epj_for_force_[ith].getPointer(n_ep_head);
                    spj_array[iw] = spj_for_force_[ith].getPointer(n_sp_head);
                }
                ret += pfunc_dispatch(tag, n_walk, (const Tepi**)epi_array, n_epi_array, (const Tepj**)epj_array, n_epj_array, (const Tspj**)spj_array, n_spj_array);
                first_loop = false;

                for(int iw=0; iw<n_walk; iw++){
                    n_epi_array_prev[iw] = n_epi_array[iw];
                    force_prev_array[iw] = force_array[iw];
                }
                n_walk_prev = n_walk;
            } // end of walk group loop
            ret += pfunc_retrieve(tag, n_walk_prev, n_epi_array_prev, force_prev_array);
        } // if(n_ipg > 0)
        else{
            ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_array[i];
            n_interaction_ep_sp_local_ += n_interaction_ep_sp_array[i];
        }
        time_profile_.calc_force__core += GetWtime() - offset_core;
        const F64 offset_copy_original_order = GetWtime();
        copyForceOriginalOrder();
        time_profile_.calc_force__copy_original_order += GetWtime() - offset_copy_original_order;
        return ret;
    }

    // for short
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkImpl(TagForceShort,
                           Tfunc_dispatch pfunc_dispatch,
                           Tfunc_retrieve pfunc_retrieve,
                           const S32 tag_max,
                           const S32 n_walk_limit,
                           const bool clear){
        const F64 offset_core = GetWtime();
        static bool first = true;
        S32 ret = 0;
        S32 tag = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        force_org_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        n_walk_local_ += n_ipg;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        static std::vector<S32> walk_grp;
        static std::vector<S32> walk_grp_disp;
        walk_grp.resize(n_loop_max); // group of walk (walk_grp[i] < n_walk_limit)
        walk_grp_disp.resize(n_loop_max+1);
        static S32  * iw2ith;
        static S32  * iw2cnt;
        static S32  * n_epi_array;
        static S32  * n_epi_array_prev;
        static S32  * n_epj_array;
        static S32 ** n_epj_disp_thread; // [n_thread][n_walk]
        static Tepi ** epi_array; // array of pointer *[n_walk]
        static Tepj ** epj_array; // array of pointer *[n_walk]
        static Tforce ** force_array; // array of pointer *[n_walk]
        static Tforce ** force_prev_array; // array of pointer *[n_walk]
        static S32  * cnt_thread;
        static S32  * n_ep_cum_thread;
        static S64 * n_interaction_ep_ep_array;
        if(first){
            iw2ith = new S32[n_walk_limit];
            iw2cnt = new S32[n_walk_limit];
            n_epi_array = new S32[n_walk_limit];
            n_epi_array_prev = new S32[n_walk_limit];
            n_epj_array = new S32[n_walk_limit];
            n_epj_disp_thread = new S32*[n_thread];
            epi_array        = new Tepi*[n_walk_limit];
            epj_array        = new Tepj*[n_walk_limit];
            force_array      = new Tforce*[n_walk_limit];
            force_prev_array = new Tforce*[n_walk_limit];
            for(int i=0; i<n_thread; i++){
                n_epj_disp_thread[i] = new S32[n_walk_limit];
            }
            cnt_thread = new S32[n_thread];
            n_ep_cum_thread = new S32[n_thread];
            n_interaction_ep_ep_array = new S64[n_thread];
            first = false;
        }
        walk_grp_disp[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            walk_grp[wg] = n_ipg / n_loop_max + 1;
            walk_grp_disp[wg+1] = walk_grp_disp[wg] + walk_grp[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
            walk_grp[wg] = n_ipg / n_loop_max;
            walk_grp_disp[wg+1] = walk_grp_disp[wg] + walk_grp[wg];
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_array[i] = 0;
        }
        bool first_loop = true;
        S32 n_walk_prev = 0;
        if(n_ipg > 0){
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = walk_grp[wg];
                const S32 walk_grp_head = walk_grp_disp[wg];
                for(int i=0; i<n_thread; i++){
                    n_ep_cum_thread[i] = cnt_thread[i] = 0;
                    n_epj_disp_thread[i][0] = 0;
                    epj_for_force_[i].clearSize();
                }
                const F64 offset_calc_force__core__walk_tree = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4) 
#endif
                for(int iw=0; iw<n_walk; iw++){
                    const S32 id_ipg = walk_grp_head + iw;
                    const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_; 
                    const S32 ith = Comm::getThreadNum();
                    n_epi_array[iw] = ipg_[id_ipg].n_ptcl_;
                    epi_array[iw]   = epi_sorted_.getPointer(first_adr_ip);
                    force_array[iw] = force_sorted_.getPointer(first_adr_ip);
                    makeInteractionList(id_ipg, false);
                    n_epj_array[iw] = epj_for_force_[ith].size() - n_ep_cum_thread[ith];
                    n_ep_cum_thread[ith] = epj_for_force_[ith].size();
                    n_epj_disp_thread[ith][cnt_thread[ith]+1] = n_ep_cum_thread[ith];
                    n_interaction_ep_ep_array[ith] += (S64)n_epj_array[iw]*(S64)n_epi_array[iw];
                    iw2ith[iw] = ith;
                    iw2cnt[iw] = cnt_thread[ith];
                    cnt_thread[ith]++;
                }
                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_array_prev, force_prev_array);
                } // retrieve
                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    epj_array[iw] = epj_for_force_[ith].getPointer(n_ep_head);
                }
                ret += pfunc_dispatch(tag, n_walk, (const Tepi**)epi_array, n_epi_array, (const Tepj**)epj_array, n_epj_array);
                first_loop = false;

                for(int iw=0; iw<n_walk; iw++){
                    n_epi_array_prev[iw] = n_epi_array[iw];
                    force_prev_array[iw] = force_array[iw];
                }
                n_walk_prev = n_walk;
            } // end of walk group loop
            ret += pfunc_retrieve(tag, n_walk_prev, n_epi_array_prev, force_prev_array);
        } // if(n_ipg > 0)
        else{
            ni_ave_ = nj_ave_ = n_interaction_ep_ep_ = n_interaction_ep_sp_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_array[i];
        }
        time_profile_.calc_force__core += GetWtime() - offset_core;
        const F64 offset_copy_original_order = GetWtime();
        copyForceOriginalOrder();
        time_profile_.calc_force__copy_original_order += GetWtime() - offset_copy_original_order;
        return ret;
    }


#if 0    
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
        AllGatherParticle(epj_tmp, n_epj_tmp, epj_org_.getPointer(), n_loc_tot_, dinfo.getPosRootDomain().getFullLength(), pos_root_cell_, pa);
        pfunc_ep_ep(epi_org_.getPointer(), n_loc_tot_, epj_tmp, n_epj_tmp, force);
        delete [] epj_tmp;
    }
#endif

#if 0
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
#endif
    
#if 0
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticleImpl(TagSearchLongScatter, const Tptcl & ptcl, S32 & nnp){
        const F64vec pos_target = ptcl.getPos();
        const S32 id_thread = Comm::getThreadNum();

        const S32 adr = 0;
        const S32 size_old = epj_neighbor_[id_thread].size();
/*
        SearchNeighborListOneParticleScatter(pos_target,    tc_glb_.getPointer(),       adr, 
                                             epj_sorted_,   epj_neighbor_, n_leaf_limit_);
*/
/*
        SearchNeighborListOneParticleScatter(pos_target,    tc_glb_.getPointer(),       
                                             tp_glb_.getPointer(), adr, 
                                             epj_sorted_,   epj_neighbor_[id_thread], n_leaf_limit_);
*/
	bool error = false;
	SearchNeighborListOneParticleScatter(pos_target,    tc_glb_.getPointer(),       
                                             tp_glb_.getPointer(), adr, 
                                             epj_sorted_,   epj_neighbor_[id_thread], n_leaf_limit_,
					     error);
	if(error){ nnp = -1; }
	else{
	    nnp = epj_neighbor_[id_thread].size() - size_old;
	}
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tptcl>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){
        const S32 id_thread = Comm::getThreadNum();
        const S32 head = epj_neighbor_[id_thread].size();
	//std::cerr<<"head="<<head<<std::endl;
        S32 nnp = 0;
        getNeighborListOneParticleImpl(typename TSM::search_type(), ptcl, nnp);
        epj = epj_neighbor_[id_thread].getPointer(head);
	if(nnp == -1){
	    epj_neighbor_[id_thread].clearSize();
	}
        return nnp;
    }
#else
    
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
					      epj_sorted_,   epj_neighbor_[id_thread], n_leaf_limit_);
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
    getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj){
	return getNeighborListOneParticleImpl(typename TSM::search_type(), ptcl, epj);
    }
    
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
	SearchNeighborListOneParticleScatter(pos_target,    tc_glb_.getPointer(),
                                             tp_glb_.getPointer(), adr, 
                                             epj_sorted_,   epj_neighbor_[id_thread], n_leaf_limit_);
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
	SearchNeighborListOneParticleGather(pos_target,    tc_glb_.getPointer(),
                                             tp_glb_.getPointer(), adr, 
                                             epj_sorted_,   epj_neighbor_[id_thread], n_leaf_limit_);
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
	SearchNeighborListOneParticleSymmetry(pos_target,    tc_glb_.getPointer(),
					      tp_glb_.getPointer(), adr, 
					      epj_sorted_,   epj_neighbor_[id_thread], n_leaf_limit_);
	S32 nnp = epj_neighbor_[id_thread].size();
        epj = epj_neighbor_[id_thread].getPointer();
        return nnp;
    }    
    */
#endif

    
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
                                            epj_sorted_,   epj_neighbor_, n_leaf_limit_);
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
    freeMem(){
#ifndef USE_MEMORY_POOL
	tp_buf_.freeMem();
#endif
#ifndef REMOVE_TP_LOC
	tp_loc_.freeMem();
#endif
	tp_glb_.freeMem();
	tc_loc_.freeMem();
	tc_glb_.freeMem();
	epi_sorted_.freeMem();
	//epi_org_.freeMem();
	epj_sorted_.freeMem();
	epj_org_.freeMem();  
	spj_sorted_.freeMem();
	spj_org_.freeMem();
	ipg_.freeMem();
	epj_send_.freeMem();
	epj_recv_.freeMem();
	spj_send_.freeMem();
	spj_recv_.freeMem();
	force_sorted_.freeMem();
	force_org_.freeMem();
	epjr_sorted_.freeMem();
	epjr_send_.freeMem();
	epjr_recv_.freeMem();
	epjr_recv_1st_buf_.freeMem();
	epjr_recv_2nd_buf_.freeMem();
	epj_recv_1st_buf_.freeMem();
	epj_recv_2nd_buf_.freeMem();
        const S32 n_thread = Comm::getNumberOfThread();
	for(S32 i=0; i<n_thread; i++){
	    id_ep_send_buf_[i].freeMem();
	    id_sp_send_buf_[i].freeMem();
	    epj_for_force_[i].freeMem();
	    spj_for_force_[i].freeMem();
	    ep_send_buf_for_scatter_[i].freeMem();
	    shift_image_domain_[i].freeMem();
	    epjr_send_buf_[i].freeMem();
	    epjr_send_buf_for_scatter_[i].freeMem();
	    epjr_recv_1st_sorted_[i].freeMem();
	    epj_send_buf_[i].freeMem();
	    id_ptcl_send_[i].freeMem();
	    shift_image_box_[i].freeMem();
	    ip_disp_[i].freeMem();
	    tp_scatter_[i].freeMem();
	    //tc_recv_1st_[i].freeMem();
	    epj_recv_1st_sorted_[i].freeMem();
	}
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    reallocMem(){
#ifndef USE_MEMORY_POOL
	tp_buf_.reallocMem();
#endif
#ifndef REMOVE_TP_LOC
	tp_loc_.reallocMem();
#endif
	tp_glb_.reallocMem();
	tc_loc_.reallocMem();
	//tc_glb_.reallocMem();
	epi_sorted_.reallocMem();
	//epi_org_.reallocMem();
	epj_sorted_.reallocMem();
	epj_org_.reallocMem();
	spj_sorted_.reallocMem();
	//spj_org_.reallocMem();
	ipg_.reallocMem();
	epj_send_.reallocMem();
	epj_recv_.reallocMem();
	spj_send_.reallocMem();
	spj_recv_.reallocMem();
	force_sorted_.reallocMem();
	force_org_.reallocMem();
	epjr_sorted_.reallocMem();
	epjr_send_.reallocMem();
	epjr_recv_.reallocMem();
	epjr_recv_1st_buf_.reallocMem();
	epjr_recv_2nd_buf_.reallocMem();
	epj_recv_1st_buf_.reallocMem();
	epj_recv_2nd_buf_.reallocMem();
        const S32 n_thread = Comm::getNumberOfThread();
	for(S32 i=0; i<n_thread; i++){
	    id_ep_send_buf_[i].reallocMem();
	    id_sp_send_buf_[i].reallocMem();
	    epj_for_force_[i].reallocMem();
	    spj_for_force_[i].reallocMem();
	    ep_send_buf_for_scatter_[i].reallocMem();
	    shift_image_domain_[i].reallocMem();
	    epjr_send_buf_[i].reallocMem();
	    epjr_send_buf_for_scatter_[i].reallocMem();
	    epjr_recv_1st_sorted_[i].reallocMem();
	    epj_send_buf_[i].reallocMem();
	    id_ptcl_send_[i].reallocMem();
	    shift_image_box_[i].reallocMem();
	    ip_disp_[i].reallocMem();
	    tp_scatter_[i].reallocMem();
	    //tc_recv_1st_[i].reallocMem();
	    epj_recv_1st_sorted_[i].reallocMem();
	}
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    clearSizeOfArray(){
        comm_table_.clearSize();
#ifndef USE_MEMORY_POOL
        tp_buf_.clearSize();
#endif
#ifndef REMOVE_TP_LOC
        tp_loc_.clearSize();
#endif
        tp_glb_.clearSize();
        tc_loc_.clearSize();
        tc_glb_.clearSize();
        epi_sorted_.clearSize();
        //epi_org_.clearSize();
        epj_sorted_.clearSize();
        epj_org_.clearSize();
        spj_sorted_.clearSize();
        spj_org_.clearSize();
        ipg_.clearSize();
        epj_send_.clearSize();
        epj_recv_.clearSize();
        spj_send_.clearSize();
        spj_recv_.clearSize();
        force_sorted_.clearSize();
        force_org_.clearSize();
        epjr_sorted_.clearSize();
        epjr_send_.clearSize();
        epjr_recv_.clearSize();
        epjr_recv_1st_buf_.clearSize();
        epjr_recv_2nd_buf_.clearSize();
        epj_recv_1st_buf_.clearSize();
        epj_recv_2nd_buf_.clearSize();
        const S32 n_thread = Comm::getNumberOfThread();
        for(S32 i=0; i<n_thread; i++){
            id_ep_send_buf_[i].clearSize();
            id_sp_send_buf_[i].clearSize();
            epj_for_force_[i].clearSize();
            spj_for_force_[i].clearSize();
            ep_send_buf_for_scatter_[i].clearSize();
            shift_image_domain_[i].clearSize();
            epjr_send_buf_[i].clearSize();
            epjr_send_buf_for_scatter_[i].clearSize();
            epjr_recv_1st_sorted_[i].clearSize();
            epj_send_buf_[i].clearSize();
            id_ptcl_send_[i].clearSize();
            shift_image_box_[i].clearSize();
            ip_disp_[i].clearSize();
            tp_scatter_[i].clearSize();
            //tc_recv_1st_[i].clearSize();
            epj_recv_1st_sorted_[i].clearSize();
        }
    }

}