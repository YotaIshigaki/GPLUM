#pragma once
#include <unistd.h>

#ifdef SUNWAY
extern "C"{
    #include <athread.h>
    #include <cpe_func.h>
    void SLAVE_FUN(CopyIndirect)(void *);
    void SLAVE_FUN(CopyIndirectInverse)(void *);
    //void SLAVE_FUN(CopyIndirectInverse2)(void *);
    void SLAVE_FUN(CopyIndirect2)(void *);
    void SLAVE_FUN(CopyDirect)(void *);
    void SLAVE_FUN(CopyStride)(void *);
    void SLAVE_FUN(GenMortonKey)(void *);
}
#endif

#ifdef SUNWAY_PREFETCH
enum{
    L1ROFF = 64,
    L2ROFF = 256,
    L1WOFF = 64,
    L2WOFF = 256,
};
#endif

namespace ParticleSimulator{
    F64 wtime_make_key_local_tree = 0.0;
    F64 wtime_sort_local_tree = 0.0;
    F64 wtime_prefixsum_recorder = 0.0;
    F64 wtime_dispatch = 0.0;
    F64 wtime_retrieve = 0.0;
    F64 wtime_calc_force = 0.0;
    F64 wtime_t0, wtime_t1, wtime_t2, wtime_t3, wtime_t4,
        wtime_t5, wtime_t6, wtime_t7, wtime_t8;

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve, class Tfp>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>
    ::calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe
    (Tfunc_dispatch pfunc_dispatch,
     Tfunc_retrieve pfunc_retrieve,
     const S32 tag_max,
     ParticleSystem<Tfp> & psys,
     DomainInfo & dinfo,
     const S32 n_walk_limit,
     const bool clear_force,
     const bool reuse)
    {
        const F64 wtime_offset_0 = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        S32 ret = 0;
        if(!reuse){
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0) std::cerr<<"OK0 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            setParticleLocalTreeImpl(psys, epi_org_, epj_org_, true); // new

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK1 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            setRootCell(dinfo); // original

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK2 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            //mortonSortLocalTreeOnly(); // original
            mortonSortLocalTreeOnlyImpl(epi_org_, epj_org_, epi_sorted_, epj_sorted_loc_); // new

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK3 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            mortonSortFP<Tfp>(psys);    // new

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK4 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            linkCellLocalTreeOnly();   // original


            F64 wtime_offset = GetWtime();
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK5 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            //calcMomentLocalTreeOnly(); // original
            CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(), epj_sorted_loc_.getPointer(), lev_max_loc_, n_leaf_limit_); //new
            time_profile_.calc_moment_local_tree += GetWtime() - wtime_offset;

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK6 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            addMomentAsSpLocalTreeImpl(typename TSM::force_type()); // original
            wtime_calc_force += MPI::Wtime() - wtime_offset_0;
            
            // MUST MODIFY TO USE ALLGATHER of monopole
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK7 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            exchangeLocalEssentialTree3(dinfo); // original
            wtime_offset = MPI::Wtime(); // used for wtime_calc_force

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK8 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            setLocalEssentialTreeToGlobalTree2(); // original

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK9 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            mortonSortGlobalTreeOnly3(); // new

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK10 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            linkCellGlobalTreeOnly(); // original

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK11 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            calcMomentGlobalTreeOnly(); // original

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK12 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            addMomentAsSpGlobalTreeImpl(typename TSM::force_type()); // original

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK13 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            makeIPGroup(); // original

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK14 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            makeAllInteractionListId3(true); // new

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK15 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            wtime_calc_force += MPI::Wtime() - wtime_offset;
            ret += calcForceUsingIdListMultiWalkIndex3(pfunc_dispatch, pfunc_retrieve, n_walk_limit, mw_info_, clear_force, reuse, true); // new

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK16 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
        }
        else{
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK17 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            setParticleLocalTreeImpl(psys, epi_sorted_, epj_sorted_loc_, true); // new
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK18 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            F64 wtime_offset = GetWtime();
            for(S32 i=0; i<n_loc_tot_; i++){
                const S32 adr = adr_epj_loc2glb_[i];
                epj_sorted_[adr] = epj_sorted_loc_[i];
            }
            time_profile_.morton_sort_global_tree += GetWtime() - wtime_offset;
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK19 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            wtime_offset = GetWtime();
            CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(), epj_sorted_loc_.getPointer(), lev_max_loc_, n_leaf_limit_);
            time_profile_.calc_moment_local_tree += GetWtime() - wtime_offset;
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK20 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            addMomentAsSpLocalTreeImpl(typename TSM::force_type()); //original
            wtime_calc_force += MPI::Wtime() - wtime_offset_0;
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK21 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            exchangeLocalEssentialTree3(dinfo, true); //original
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK22 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            const S32 n_epj_add = comm_table_.n_disp_ep_recv_[n_proc];
            const S32 n_spj_add = comm_table_.n_disp_sp_recv_[n_proc];
            wtime_offset = GetWtime(); // used both for time_profile_.morton_sort_global_tree and wtime_calc_force 
            for(S32 i=0; i<n_epj_add; i++){
                const S32 adr = adr_epj_buf2glb_[i];
                epj_sorted_[adr] = comm_table_.ep_recv_[i];
            }
            for(S32 i=0; i<n_spj_add; i++){
                const S32 adr = adr_spj_buf2glb_[i];
                spj_sorted_[adr] = comm_table_.sp_recv_[i];
            }
            time_profile_.morton_sort_global_tree += GetWtime() - wtime_offset;
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK23 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            calcMomentGlobalTreeOnly(); //original
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK24 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            addMomentAsSpGlobalTreeImpl(typename TSM::force_type()); //original
            wtime_calc_force += MPI::Wtime() - wtime_offset;
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK25 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            ret += calcForceUsingIdListMultiWalkIndex3(pfunc_dispatch, pfunc_retrieve, n_walk_limit, mw_info_, clear_force, reuse, true); // new
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK26 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
        }
        //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"CHECK 17"<<std::endl;
        F64 wtime_offset_1 = GetWtime();
        for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_sorted_[i]);
        //Comm::barrier();if(Comm::getRank()==0)std::cerr<<"CHECK 18"<<std::endl;
        time_profile_.calc_force__copy_original_order += GetWtime() - wtime_offset_1;
        time_profile_.calc_force_all += GetWtime() - wtime_offset_0;
        return ret;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfp>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    mortonSortFP(ParticleSystem<Tfp> & sys){
        const F64 time_offset = GetWtime();
#ifdef SUNWAY
//#if 0
        Tfp * sys_buf = new Tfp[n_loc_tot_];
        unsigned long arg[7];
        arg[0] = (unsigned long)(n_loc_tot_);
        arg[1] = (unsigned long)(&sys[0]);
        arg[2] = (unsigned long)(&sys_buf[0]);
        arg[3] = (unsigned long)(sizeof(sys[0]));
        __real_athread_spawn((void*)slave_CopyDirect, arg);
        athread_join();

    #if 1
        arg[1] = (unsigned long)(&sys_buf[0]);
        arg[2] = (unsigned long)(&sys[0]);
        arg[3] = (unsigned long)(sizeof(sys[0]));
        arg[4] = (unsigned long)(adr_ptcl_of_tp_loc_.getPointer());
        __real_athread_spawn((void*)slave_CopyIndirectInverse, arg);
        athread_join();
    #else
        arg[1] = (unsigned long)(&sys_buf[0]);
        arg[2] = (unsigned long)(&sys[0]);
        arg[3] = (unsigned long)(sizeof(sys[0]));
        arg[4] = (unsigned long)(((int*)tp_loc_.getPointer())+2);
        arg[5] = (unsigned long)(sizeof(U32)); // type of address
        arg[6] = (unsigned long)(sizeof(tp_loc_[0])-sizeof(U32)); // stride
        __real_athread_spawn((void*)slave_CopyIndirectInverse2, arg);
        athread_join();        
    #endif
        delete [] sys_buf;
#else
        Tfp * sys_buf = new Tfp[n_loc_tot_];
        std::memcpy(sys_buf, &sys[0], sizeof(Tfp)*n_loc_tot_);
        for(S32 i=0; i<n_loc_tot_; i++){
            const S32 adr = tp_loc_[i].adr_ptcl_;
            sys[i] = sys_buf[adr];
        }
        delete [] sys_buf;
#endif        
        time_profile_.morton_sort_FP += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeAllInteractionListId3(const bool clear){
        const F64 time_offset = GetWtime();
        const S64 n_ipg = ipg_.size();
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
        const S32 sp_tree_offset = spj_sorted_.size() - tc_glb_.size();
        if(n_ipg > 0 && clear){
            id_recorder_for_interaction_.resizeNoInitialize(n_ipg);
            for(S32 ith=0; ith<Comm::getNumberOfThread(); ith++){
                id_epj_recorder_for_force_[ith].clearSize();
                id_spj_recorder_for_force_[ith].clearSize();
                n_epj_recorder_for_force_[ith].clearSize();
                n_spj_recorder_for_force_[ith].clearSize();
                n_disp_epj_recorder_for_force_[ith].clearSize();
                n_disp_spj_recorder_for_force_[ith].clearSize();
                n_disp_epj_recorder_for_force_[ith].push_back(0);
                n_disp_spj_recorder_for_force_[ith].push_back(0);
            }
        }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4)
#endif	    
        for(S32 ig=0; ig<n_ipg; ig++){
            const S32 ith = Comm::getThreadNum();
            const S32 adr_epj_prev = id_epj_recorder_for_force_[ith].size();
            const S32 adr_spj_prev = id_spj_recorder_for_force_[ith].size();
            const F64ort pos_target_box = (ipg_[ig]).vertex_;
            ReallocatableArray<Tepj> ep_tmp;
            ReallocatableArray<Tspj> sp_tmp;
            MakeListUsingTreeRecursive<TSM, MAKE_LIST_MODE_INTERACTION, LIST_CONTENT_ID, TreeCell<Tmomglb>, TreeParticle, Tepj, Tspj>
                (tc_glb_, 0, tp_glb_, 
                 epj_sorted_, ep_tmp, id_epj_recorder_for_force_[ith], 
                 spj_sorted_, sp_tmp, id_spj_recorder_for_force_[ith], 
                 pos_target_box, r_crit_sq, n_leaf_limit_, sp_tree_offset);

            n_disp_epj_recorder_for_force_[ith].push_back( id_epj_recorder_for_force_[ith].size() );
            n_disp_spj_recorder_for_force_[ith].push_back( id_spj_recorder_for_force_[ith].size() );
            const S32 n_epj = id_epj_recorder_for_force_[ith].size() - adr_epj_prev;
            const S32 n_spj = id_spj_recorder_for_force_[ith].size() - adr_spj_prev;
            n_epj_recorder_for_force_[ith].push_back(n_epj);
            n_spj_recorder_for_force_[ith].push_back(n_spj);

            F64 mass = 0.0;
            for(S32 i=adr_epj_prev; i<id_epj_recorder_for_force_[ith].size(); i++){
                S32 adr = id_epj_recorder_for_force_[ith][i];
                //std::cout<<"i= "<<i<<" adr0= "<<adr<<" epj_sorted_[adr].mass= "<<epj_sorted_[adr].mass<<std::endl;
                mass += epj_sorted_[adr].mass;
            }
            for(S32 i=adr_spj_prev; i<id_spj_recorder_for_force_[ith].size(); i++){
                S32 adr = id_spj_recorder_for_force_[ith][i];
                //std::cout<<"i= "<<i<<" adr0= "<<adr<<" spj_sorted_[adr].mass= "<<spj_sorted_[adr].mass<<std::endl;                
                mass += spj_sorted_[adr].mass;
            }
            // FOR DEBUG
            //if(Comm::getRank()%100==0)  std::cerr<<"mass= "<<mass<<std::endl;
        }
        time_profile_.make_all_interaction_list_id += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceUsingIdListMultiWalkIndex3
    (Tfunc_dispatch pfunc_dispatch,
     Tfunc_retrieve pfunc_retrieve,
     const S32 n_walk_limit,
     MultiWalkInfo<Tepi, Tepj, Tspj, Tforce> & mw_info,
     const bool clear_force,
     const bool reuse_list,
     const bool flag_retrieve)
    {
        F64 wtime_offset_out = GetWtime();
        F64 wtime_offset_in = GetWtime();
        S32 ret = 0;
        const S32 n_ipg = ipg_.size();
        n_epi_recorder_for_force_[0].resizeNoInitialize(n_ipg);
        n_disp_epi_recorder_for_force_[0].resizeNoInitialize(n_ipg+1);
        n_disp_epi_recorder_for_force_[0][0] = 0;
        for(S32 i=0; i<n_ipg; i++){
            n_epi_recorder_for_force_[0][i] = ipg_[i].n_ptcl_;
            n_disp_epi_recorder_for_force_[0][i+1] = n_disp_epi_recorder_for_force_[0][i] + n_epi_recorder_for_force_[0][i];
        }
        const S32 ni_tot = n_disp_epi_recorder_for_force_[0][n_ipg];
        force_sorted_.resizeNoInitialize(ni_tot);
        if(clear_force){
            for(S32 ip=0; ip<ni_tot; ip++){
                force_sorted_[ip].clear();
            }
        }
        wtime_prefixsum_recorder = GetWtime() - wtime_offset_in;
        wtime_offset_in = GetWtime();
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
        Comm::barrier();
        if(Comm::getRank()==0)std::cerr<<"OK1 @calcForceUsingIdListMultiWalkIndex3"<<std::endl;
#endif
        pfunc_dispatch(n_ipg,
                       epi_sorted_.getPointer(),
                       n_epi_recorder_for_force_[0].getPointer(),
                       n_disp_epi_recorder_for_force_[0].getPointer(),
                       id_epj_recorder_for_force_[0].getPointer(),
                       n_epj_recorder_for_force_[0].getPointer(),
                       n_disp_epj_recorder_for_force_[0].getPointer(),
                       id_spj_recorder_for_force_[0].getPointer(),
                       n_spj_recorder_for_force_[0].getPointer(),
                       n_disp_spj_recorder_for_force_[0].getPointer(),
                       epj_sorted_.getPointer(),
                       spj_sorted_.getPointer());
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
        Comm::barrier();
        if(Comm::getRank()==0)std::cerr<<"OK2 @calcForceUsingIdListMultiWalkIndex3"<<std::endl;
#endif
        wtime_dispatch = GetWtime() - wtime_offset_in;
        wtime_offset_in = GetWtime();
        if(flag_retrieve){
            pfunc_retrieve(ni_tot, force_sorted_.getPointer());
        }
        wtime_retrieve = GetWtime() - wtime_offset_in;
        const S32 n_epj_tot = n_disp_epj_recorder_for_force_[0][n_ipg];
        const S32 n_spj_tot = n_disp_spj_recorder_for_force_[0][n_ipg];        
        n_epi_loc_ave_ = ni_tot / n_ipg;
        n_epj_loc_ave_ = n_epj_tot / n_ipg;
        n_spj_loc_ave_ = n_spj_tot / n_ipg;
        time_profile_.calc_force = GetWtime() - wtime_offset_out;
        return ret;
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tpsys>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setParticleLocalTreeImpl(const Tpsys & psys,
                             ReallocatableArray<Tepi> & epi,
                             ReallocatableArray<Tepj> & epj,
                             const bool clear){
        const F64 time_offset = GetWtime();
        const S32 nloc = psys.getNumberOfParticleLocal();
        if(clear){ n_loc_tot_ = 0;}
        const S32 offset = n_loc_tot_;
        n_loc_tot_ += nloc;
        epi.resizeNoInitialize(n_loc_tot_);
        epj.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for
#endif	    
            for(S32 i=0; i<nloc; i++){
                epi[i].copyFromFP( psys[i] );
                epj[i].copyFromFP( psys[i] );
            }
        }
        else{
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for
#endif	    
            for(S32 i=0; i<nloc; i++){
                epi[i+offset].copyFromFP( psys[i] );
                epj[i+offset].copyFromFP( psys[i] );
            }
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"nloc="<<nloc<<std::endl;
        std::cout<<"n_loc_tot_="<<n_loc_tot_<<std::endl;
#endif
        time_profile_.set_particle_local_tree += GetWtime() - time_offset;
    }


    ///////////////////////////////
    /// morton sort global tree ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    mortonSortGlobalTreeOnly3(){
        F64 time_offset = GetWtime();
        tp_glb_.resizeNoInitialize(n_glb_tot_);
        tp_buf_.resizeNoInitialize(n_glb_tot_);
#ifdef USE_STD_SORT
        /*
        std::sort(tp_glb_.getPointer(), tp_glb_.getPointer()+n_glb_tot_, 
		  [](const TreeParticle & l, const TreeParticle & r )
		  ->bool{return l.getKey() < r.getKey();} );
        */
        std::sort(tp_glb_.getPointer(), tp_glb_.getPointer()+n_glb_tot_, 
                  LessOPKEY());
#else //USE_STD_SORT
        rs_.lsdSort(tp_glb_.getPointer(), tp_buf_.getPointer(), 0, n_glb_tot_-1);
#endif //USE_STD_SORT
        epj_sorted_.resizeNoInitialize( n_glb_tot_ );
        //std::cerr<<"check 1"<<std::endl;
        if( typeid(TSM) == typeid(SEARCH_MODE_LONG)
            || typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF) 
            || typeid(TSM) == typeid(SEARCH_MODE_LONG_SCATTER) 
            || typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF_SCATTER) ){
            // Long mode
            spj_sorted_.resizeNoInitialize( spj_org_.size() );
            adr_epj_org2glb_.resizeNoInitialize(n_loc_tot_); // new line
            const S32 n_proc = Comm::getNumberOfProc();
            const S32 n_epj_add = comm_table_.n_disp_ep_recv_[n_proc];
            const S32 n_spj_add = comm_table_.n_disp_sp_recv_[n_proc];
            const S32 offset_spj = n_loc_tot_ + n_epj_add;
            //if(Comm::getRank() == 0) std::cerr<<"offset_spj= "<<offset_spj<<std::endl;
            adr_epj_buf2glb_.resizeNoInitialize(n_epj_add);
            adr_spj_buf2glb_.resizeNoInitialize(n_spj_add);

//#ifdef SUNWAY
#if 0
            // under construction
            adr_epj_loc2glb_.resizeNoInitialize(n_loc_tot_);
            unsigned long arg[6];
            arg[0] = (unsigned long)(n_glb_tot_);
            arg[1] = (unsigned long)(n_loc_tot_);
            arg[2] = (unsigned long)(epj_org_.getPointer());
            arg[3] = (unsigned long)(sizeof(epj_org_[0]));
            arg[4] = (unsigned long)(spj_org_.getPointer());
            arg[5] = (unsigned long)(sizeof(spj_org_[0]));
            arg[6] = (unsigned long)(((int*)tp_glb_.getPointer())+2);
            arg[7] = (unsigned long)(sizeof(tp_glb_[0]));
            arg[8] = (unsigned long)(epj_sorted_.getPointer());
#else // SUNWAY
    #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
    #pragma omp parallel for
    #endif	    
            for(S32 i=0; i<n_glb_tot_; i++){
                const U32 adr = tp_glb_[i].adr_ptcl_;
                if( GetMSB(adr) == 0){
                    epj_sorted_[i] = epj_org_[adr];
                    // new lines
                    if( adr < n_loc_tot_){
                        adr_epj_org2glb_[adr] = i;
                    }
                    else{
                        adr_epj_buf2glb_[adr-n_loc_tot_] = i; 
                    }
                    // new lines
                }
                else{
                    const U32 adr_new = ClearMSB(adr);
                    spj_sorted_[i] = spj_org_[adr_new];
                    adr_spj_buf2glb_[adr_new-offset_spj] = i; // new line
                }
            }
            // new lines
            adr_epj_loc2glb_.resizeNoInitialize(n_loc_tot_);
            for(S32 i=0; i<n_loc_tot_; i++){
                const U32 adr_org = tp_loc_[i].adr_ptcl_;
                const U32 adr_glb = adr_epj_org2glb_[adr_org];
                adr_epj_loc2glb_[i] = adr_glb;
            }
            // new lines
#endif // SUNWAY
        }
        else{
            // short mode
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif	    
            for(S32 i=0; i<n_glb_tot_; i++){
                const U32 adr = tp_glb_[i].adr_ptcl_;
                epj_sorted_[i] = epj_org_[adr];
            }
        }
        //std::cerr<<"check 2"<<std::endl;
        time_profile_.morton_sort_global_tree += GetWtime() - time_offset;
        time_profile_.make_global_tree += GetWtime() - time_offset;
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    mortonSortLocalTreeOnlyImpl(const ReallocatableArray<Tepi> & _epi_org,
                                const ReallocatableArray<Tepj> & _epj_org,
                                ReallocatableArray<Tepi> & _epi_sorted,
                                ReallocatableArray<Tepj> & _epj_sorted){
        const F64 time_offset = GetWtime();
        tp_loc_.resizeNoInitialize(n_loc_tot_);
        _epi_sorted.resizeNoInitialize(n_loc_tot_);
        _epj_sorted.resizeNoInitialize(n_loc_tot_);
        MortonKey::initialize( length_ * 0.5, center_);
        F64 wtime_offset_in = GetWtime();
#ifdef SUNWAY
        //try {
        S32 myrank = Comm::getRank();
        S32 nprocs = Comm::getNumberOfProc();
#ifdef GEN_MORTON_KEY_SEQ_EXEC_FOR_DEBUG
        // [*** Note ***]
        //    This part is used to check CPEs of the nodes.
        if (myrank == 0) std::cout << "Enter slave_GenMortonKey()" << std::endl;
        char hostname[256];
        gethostname(hostname,256);
        for (S32 irank=0; irank<nprocs; irank++) {
            if (myrank == 0) 
                std::cout << "checking ... RANK " << irank
                          << " (host: " << hostname << ")" << std::endl;
            if (irank == myrank) {
#endif
                unsigned long args[6];
                args[0] = (unsigned long)(n_loc_tot_);
                args[1] = (unsigned long)(&length_);
                args[2] = (unsigned long)(&center_);
                args[3] = (unsigned long) tp_loc_.getPointer();
                args[4] = (unsigned long) _epj_org.getPointer();
                args[5] = (unsigned long) myrank;
                __real_athread_spawn((void*)slave_GenMortonKey, args);
                athread_join();
#ifdef GEN_MORTON_KEY_SEQ_EXEC_FOR_DEBUG
            }
            if (myrank == 0) 
              std::cout << "RANK " << irank << "(host: " << hostname << ") passed!!" << std::endl;
            Comm::barrier();
        }
        Comm::barrier();
        if (myrank == 0) std::cout << "Leave slave_GenMortonKey()" << std::endl;
#endif
        //}
        //catch (...) {
        //    std::ofstream flog;
        //    std::stringstream fnum;
        //    std::string fname;
        //    fnum << std::setfill('0') << std::setw(5) << Comm::getRank();
        //    fname = "catch" + fnum.str() + ".txt";
        //    flog.open(fname,std::ios::trunc);

        //    char hostname[256];
        //    gethostname(hostname,256);
        //    std::cout << "Catched an error in mortonSortLocalTreeOnlyImpl." << std::endl;
        //    std::cout << "myrank = " << Comm::getRank() << std::endl;
        //    std::cout << "hostname = " << hostname << std::endl;
        //    flog << "Catched an error in mortonSortLocalTreeOnlyImpl." << std::endl;
        //    flog << "myrank = " << Comm::getRank() << std::endl;
        //    flog << "hostname = " << hostname << std::endl;

        //    flog.close();
        //}
#else //SUNWAY
  #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	
  #pragma omp parallel for
  #endif
        for(S32 i=0; i<n_loc_tot_; i++){
            tp_loc_[i].setFromEP(_epj_org[i], i);
        }
#endif //SUNWAY
        wtime_make_key_local_tree = GetWtime() - wtime_offset_in;
        
        wtime_offset_in = GetWtime();
#ifdef USE_STD_SORT
        /*
	std::sort(tp_loc_.getPointer(), tp_loc_.getPointer()+n_loc_tot_, 
		  [](const TreeParticle & l, const TreeParticle & r )
		  ->bool{return l.getKey() < r.getKey();} );
        */
	std::sort(tp_loc_.getPointer(), tp_loc_.getPointer()+n_loc_tot_, 
                  LessOPKEY());
#else //USE_STD_SORT
        tp_buf_.resizeNoInitialize(n_loc_tot_);
        rs_.lsdSort(tp_loc_.getPointer(), tp_buf_.getPointer(), 0, n_loc_tot_-1);
#endif //USE_STD_SORT
        wtime_sort_local_tree = GetWtime() - wtime_offset_in;
        
#ifdef SUNWAY
        // TODO
        adr_ptcl_of_tp_loc_.resizeNoInitialize(n_loc_tot_);
        for(S32 i=0; i<n_loc_tot_; i++){
            adr_ptcl_of_tp_loc_[i] = tp_loc_[i].adr_ptcl_;
        }        
        unsigned long arg[7];
        /*
        arg[0] = (unsigned long)(n_loc_tot_);
        arg[1] = (unsigned long)(((int*)(tp_loc_.getPointer()))+2);
        arg[2] = (unsigned long)(adr_tp_loc_org2sorted_.getPointer());
        arg[3] = (unsigned long)(sizeof(adr_tp_loc_org2sorted_[0]));
        arg[4] = (unsigned long)(sizeof(U64) + sizeof(int));
        __real_athread_spawn((void*)slave_CopyStride, arg);
        athread_join();
        */
        /*
        std::cerr<<"adr_tp_loc_org2sorted_[0]= "<<adr_tp_loc_org2sorted_[0]
                 <<" tp_loc_[0].adr_ptcl_= "<<tp_loc_[0].adr_ptcl_<<std::endl;
        for(S32 i=0; i<n_loc_tot_; i++){
            assert(adr_tp_loc_org2sorted_[i]==tp_loc_[i].adr_ptcl_);
        }
        */
        /*
        for(S32 i=0; i<n_loc_tot_; i++){
            _epi_sorted[i].id = -i-1;
            _epj_sorted[i].id = -i-1;
        }
        */
    #if 1
        arg[0] = (unsigned long)(n_loc_tot_);
        arg[1] = (unsigned long)(_epj_org.getPointer());
        arg[2] = (unsigned long)(_epj_sorted.getPointer());
        arg[3] = (unsigned long)(sizeof(_epj_org[0]));
        arg[4] = (unsigned long)(adr_ptcl_of_tp_loc_.getPointer());
        __real_athread_spawn((void*)slave_CopyIndirectInverse, arg);
        athread_join();

        arg[0] = (unsigned long)(n_loc_tot_);
        arg[1] = (unsigned long)(_epi_org.getPointer());
        arg[2] = (unsigned long)(_epi_sorted.getPointer());
        arg[3] = (unsigned long)(sizeof(_epi_org[0]));
        arg[4] = (unsigned long)(adr_ptcl_of_tp_loc_.getPointer());
        __real_athread_spawn((void*)slave_CopyIndirectInverse, arg);
        athread_join();
    #elif 0
        //unsigned long arg[7];
        arg[0] = (unsigned long)(n_loc_tot_);
        arg[1] = (unsigned long)(_epj_org.getPointer());
        arg[2] = (unsigned long)(_epj_sorted.getPointer());
        arg[3] = (unsigned long)(sizeof(_epj_org[0]));
        arg[4] = (unsigned long)(((int*)tp_loc_.getPointer())+2);
        arg[5] = (unsigned long)(sizeof(U32)); // type of address
        arg[6] = (unsigned long)(sizeof(tp_loc_[0])-sizeof(U32)); // stride 
        __real_athread_spawn((void*)slave_CopyIndirectInverse2, arg);
        athread_join();

        arg[0] = (unsigned long)(n_loc_tot_);
        arg[1] = (unsigned long)(_epi_org.getPointer());
        arg[2] = (unsigned long)(_epi_sorted.getPointer());
        arg[3] = (unsigned long)(sizeof(_epi_org[0]));
        arg[4] = (unsigned long)(((int*)tp_loc_.getPointer())+2);
        arg[5] = (unsigned long)(sizeof(U32)); // type of address
        arg[6] = (unsigned long)(sizeof(tp_loc_[0])-sizeof(U32)); // stride  
        //arg[5] = (unsigned long)(sizeof(tp_loc_[0])); // byte
        __real_athread_spawn((void*)slave_CopyIndirectInverse2, arg);
        athread_join();
    #else
        ReallocatableArray<S32> adr_tmp;
        adr_tmp.resizeNoInitialize(n_loc_tot_);
        for(S32 i=0; i<n_loc_tot_; i++){
            adr_tmp[i] = tp_loc_[i].adr_ptcl_;
        }
        arg[0] = (unsigned long)(n_loc_tot_);
        arg[1] = (unsigned long)(_epj_org.getPointer());
        arg[2] = (unsigned long)(_epj_sorted.getPointer());
        arg[3] = (unsigned long)(sizeof(_epj_org[0]));
        arg[4] = (unsigned long)(adr_tmp.getPointer());
        __real_athread_spawn((void*)slave_CopyIndirectInverse, arg);
        athread_join();
        for(S32 i=0; i<n_loc_tot_; i++){
            const S32 adr = tp_loc_[i].adr_ptcl_;
            _epi_sorted[i] = _epi_org[adr];
            //_epj_sorted[i] = _epj_org[adr];
        }
    #endif
        /*
        for(S32 i=0; i<n_loc_tot_; i++){
            const S32 adr = tp_loc_[i].adr_ptcl_;
            assert(_epi_sorted[i].id == _epi_org[adr].id);
            assert(_epj_sorted[i].id == _epj_org[adr].id);
        }
        */
#else //SUNWAY
  #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	
  #pragma omp parallel for
  #endif	
        for(S32 i=0; i<n_loc_tot_; i++){
            const S32 adr = tp_loc_[i].adr_ptcl_;
            _epi_sorted[i] = _epi_org[adr];
            _epj_sorted[i] = _epj_org[adr];
        }
#endif //SUNWAY
        time_profile_.morton_sort_local_tree += GetWtime() - time_offset;
        time_profile_.make_local_tree += GetWtime() - time_offset;
        //std::cerr<<"end of func"<<std::endl;        
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTree3(const DomainInfo & dinfo,
                                const bool reuse_list){
        const F64 time_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
        tree_outer_pos_.resizeNoInitialize(n_proc);
        if(!reuse_list){
            comm_table_.n_ep_send_[my_rank] = comm_table_.n_sp_send_[my_rank]
                = comm_table_.n_ep_recv_[my_rank] = comm_table_.n_sp_recv_[my_rank] = 0;
            F64ort outer_pos = tc_loc_[0].mom_.getVertexOut();
            Comm::allGather(&outer_pos, 1, tree_outer_pos_.getPointer());
            //Comm::allGather(&(tc_loc_[0].mom_.getVertexOut()), 1, tree_outer_pos_.getPointer());
        }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif	
        {
            const S32 ith = Comm::getThreadNum();
            F64 time_offset0 = GetWtime();
            if(!reuse_list){
                const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
                const F64 longet_len_of_outer_box = tree_outer_pos_[my_rank].getFullLength().getMax();
                const F64 r_crit_sq_2 = (longet_len_of_outer_box * longet_len_of_outer_box) / (theta_ * theta_);
                rank_proc_send_[ith].resizeNoInitialize(0);
                id_ep_send_buf_[ith].resizeNoInitialize(0);
                id_sp_send_buf_[ith].resizeNoInitialize(0);
                ReallocatableArray<Tepj> ep_list_dummy;
                ReallocatableArray<Tspj> sp_list_dummy, sp_list_dummy2;
                S32 n_ep_cum = 0;
                S32 n_sp_cum = 0;
                //id_sp_send_buf_[ith].reserveEmptyAreaAtLeast(n_proc);//reserve enough size
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
                for(S32 ib=0; ib<n_proc; ib++){
                    rank_proc_send_[ith].push_back(ib);
                    n_ep_send_[ib] = n_sp_send_[ib] = n_ep_recv_[ib] = n_sp_recv_[ib] = 0;
                    if(my_rank == ib) continue;
                    /*
                    if(my_rank==8 && ib == 12){
                        std::cerr<<"r(8->12)= "<<tree_outer_pos_[ib].getDistanceMinSQ(tree_outer_pos_[my_rank])
                                 <<std::endl;
                    }
                    */
                    if(tree_outer_pos_[ib].getDistanceMinSQ(tree_outer_pos_[my_rank]) > r_crit_sq_2){
                        //if(0){
                        id_sp_send_buf_[ith].push_back(0); // 0 means spj of root.
                        //id_sp_send_buf_[ith].pushBackNoCheck(0); // 0 means spj of root.
                    }
                    else{
                        const F64ort pos_target_domain = dinfo.getPosDomain(ib);
                        MakeListUsingTreeRecursiveTop<TSM, MAKE_LIST_MODE_LET, LIST_CONTENT_ID, TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj>
                            (tc_loc_, 0, tp_loc_, epj_sorted_loc_, ep_list_dummy,  id_ep_send_buf_[ith], sp_list_dummy, sp_list_dummy2,
                             id_sp_send_buf_[ith], pos_target_domain, r_crit_sq, n_leaf_limit_);
                    }
                    comm_table_.n_ep_send_[ib] = id_ep_send_buf_[ith].size() - n_ep_cum;
                    comm_table_.n_sp_send_[ib] = id_sp_send_buf_[ith].size() - n_sp_cum;
                    n_ep_cum = id_ep_send_buf_[ith].size();
                    n_sp_cum = id_sp_send_buf_[ith].size();
                }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif
                {
                    comm_table_.setDispSend();
                    comm_table_.setSizeSendBuf();
                }
            } // end of if(!reuse_list)
            time_profile_.make_LET_1st += GetWtime() - time_offset0;
            S32 n_ep_cnt = 0;
            S32 n_sp_cnt = 0;
            time_offset0 = GetWtime();
            for(S32 ib=0; ib<rank_proc_send_[ith].size(); ib++){
                const S32 rank        = rank_proc_send_[ith][ib];
                const S32 adr_ep_head = comm_table_.n_disp_ep_send_[rank];
                const S32 adr_ep_tail = comm_table_.n_disp_ep_send_[rank+1];
                const S32 n_ep_tmp    = n_ep_send_[rank];
                for(int ip=adr_ep_head; ip<adr_ep_tail; ip++, n_ep_cnt++){
                    comm_table_.ep_send_[ip] = epj_sorted_loc_[ id_ep_send_buf_[ith][n_ep_cnt] ];
                }
                const S32 adr_sp_head = comm_table_.n_disp_sp_send_[rank];
                const S32 adr_sp_tail = comm_table_.n_disp_sp_send_[rank+1];
                const S32 n_sp_tmp    = n_sp_send_[rank];
                for(int ip=adr_sp_head; ip<adr_sp_tail; ip++, n_sp_cnt++){
                    comm_table_.sp_send_[ip] = spj_sorted_loc_[ id_sp_send_buf_[ith][n_sp_cnt] ];
                }
            }
            time_profile_.make_LET_2nd += GetWtime() - time_offset0;
        }
        /*
        if(Comm::getRank()==8){
            for(int i=0; i<n_proc; i++){
                std::cerr<<"Comm::getRank()= "<<Comm::getRank()
                         <<" i= "<<i
                         <<" comm_table_.n_ep_send_[i]= "<<comm_table_.n_ep_send_[i]
                         <<" comm_table_.n_sp_send_[i]= "<<comm_table_.n_sp_send_[i]
                         <<std::endl;
            }
        }
        */
#if 1
        // new
        comm_table_.setSPTop(spj_sorted_loc_[0]);
        //comm_table_.exchangeLetSunWay(tc_loc_[0].mom_.getVertexOut(), reuse_list);
        //comm_table_.exchangeLetSunWay(tree_outer_pos_, theta_, reuse_list);
        comm_table_.exchangeLetSunWay(tree_outer_pos_, theta_, reuse_list, &dinfo);
#else
        // original
        comm_table_.exchangeLet(reuse_list);
#endif
        time_profile_.exchange_LET_tot += GetWtime() - time_offset;
    }
    
#if 0
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTree3(const DomainInfo & dinfo,
                                const bool reuse_list){
        const F64 time_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
        tree_outer_pos_.resizeNoInitialize(n_proc);
        if(!reuse_list){
            comm_table_.n_ep_send_[my_rank] = comm_table_.n_sp_send_[my_rank]
                = comm_table_.n_ep_recv_[my_rank] = comm_table_.n_sp_recv_[my_rank] = 0;
            F64ort outer_pos = tc_loc_[0].mom_.getVertexOut();
            Comm::allGather(&outer_pos, 1, tree_outer_pos_.getPointer());
            //Comm::allGather(&(tc_loc_[0].mom_.getVertexOut()), 1, tree_outer_pos_.getPointer());
        }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif	
        {
            const S32 ith = Comm::getThreadNum();
            F64 time_offset0 = GetWtime();
            if(!reuse_list){
                const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
                const F64 longet_len_of_outer_box = tree_outer_pos_[my_rank].getFullLength().getMax();
                const F64 r_crit_sq_2 = (longet_len_of_outer_box * longet_len_of_outer_box) / (theta_ * theta_);
                rank_proc_send_[ith].resizeNoInitialize(0);
                id_ep_send_buf_[ith].resizeNoInitialize(0);
                id_sp_send_buf_[ith].resizeNoInitialize(0);
                ReallocatableArray<Tepj> ep_list_dummy;
                ReallocatableArray<Tspj> sp_list_dummy, sp_list_dummy2;
                S32 n_ep_cum = 0;
                S32 n_sp_cum = 0;
                //id_sp_send_buf_[ith].reserveEmptyAreaAtLeast(n_proc);//reserve enough size
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
                for(S32 ib=0; ib<n_proc; ib++){
                    rank_proc_send_[ith].push_back(ib);
                    n_ep_send_[ib] = n_sp_send_[ib] = n_ep_recv_[ib] = n_sp_recv_[ib] = 0;
                    if(my_rank == ib) continue;
                    /*
                    if(my_rank==8 && ib == 12){
                        std::cerr<<"r(8->12)= "<<tree_outer_pos_[ib].getDistanceMinSQ(tree_outer_pos_[my_rank])
                                 <<std::endl;
                    }
                    */
                    if(tree_outer_pos_[ib].getDistanceMinSQ(tree_outer_pos_[my_rank]) > r_crit_sq_2){
                        //if(0){
                        id_sp_send_buf_[ith].push_back(0); // 0 means spj of root.
                        //id_sp_send_buf_[ith].pushBackNoCheck(0); // 0 means spj of root.
                    }
                    else{
                        const F64ort pos_target_domain = dinfo.getPosDomain(ib);
                        MakeListUsingTreeRecursiveTop<TSM, MAKE_LIST_MODE_LET, LIST_CONTENT_ID, TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj>
                            (tc_loc_, 0, tp_loc_, epj_sorted_loc_, ep_list_dummy,  id_ep_send_buf_[ith], sp_list_dummy, sp_list_dummy2,
                             id_sp_send_buf_[ith], pos_target_domain, r_crit_sq, n_leaf_limit_);
                    }
                    comm_table_.n_ep_send_[ib] = id_ep_send_buf_[ith].size() - n_ep_cum;
                    comm_table_.n_sp_send_[ib] = id_sp_send_buf_[ith].size() - n_sp_cum;
                    n_ep_cum = id_ep_send_buf_[ith].size();
                    n_sp_cum = id_sp_send_buf_[ith].size();
                }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp single
#endif
                {
                    comm_table_.setDispSend();
                    comm_table_.setSizeSendBuf();
                }
            } // end of if(!reuse_list)
            time_profile_.make_LET_1st += GetWtime() - time_offset0;
            S32 n_ep_cnt = 0;
            S32 n_sp_cnt = 0;
            time_offset0 = GetWtime();
            for(S32 ib=0; ib<rank_proc_send_[ith].size(); ib++){
                const S32 rank        = rank_proc_send_[ith][ib];
                const S32 adr_ep_head = comm_table_.n_disp_ep_send_[rank];
                const S32 adr_ep_tail = comm_table_.n_disp_ep_send_[rank+1];
                const S32 n_ep_tmp    = n_ep_send_[rank];
                for(int ip=adr_ep_head; ip<adr_ep_tail; ip++, n_ep_cnt++){
                    comm_table_.ep_send_[ip] = epj_sorted_loc_[ id_ep_send_buf_[ith][n_ep_cnt] ];
                }
                const S32 adr_sp_head = comm_table_.n_disp_sp_send_[rank];
                const S32 adr_sp_tail = comm_table_.n_disp_sp_send_[rank+1];
                const S32 n_sp_tmp    = n_sp_send_[rank];
                for(int ip=adr_sp_head; ip<adr_sp_tail; ip++, n_sp_cnt++){
                    comm_table_.sp_send_[ip] = spj_sorted_loc_[ id_sp_send_buf_[ith][n_sp_cnt] ];
                }
            }
            time_profile_.make_LET_2nd += GetWtime() - time_offset0;
        }
        /*
        if(Comm::getRank()==8){
            for(int i=0; i<n_proc; i++){
                std::cerr<<"Comm::getRank()= "<<Comm::getRank()
                         <<" i= "<<i
                         <<" comm_table_.n_ep_send_[i]= "<<comm_table_.n_ep_send_[i]
                         <<" comm_table_.n_sp_send_[i]= "<<comm_table_.n_sp_send_[i]
                         <<std::endl;
            }
        }
        */
#if 1
        // new
        comm_table_.setSPTop(spj_sorted_loc_[0]);
        //comm_table_.exchangeLetSunWay(tc_loc_[0].mom_.getVertexOut(), reuse_list);
        comm_table_.exchangeLetSunWay(tree_outer_pos_, theta_, reuse_list);
#else
        // original
        comm_table_.exchangeLet(reuse_list);
#endif
        time_profile_.exchange_LET_tot += GetWtime() - time_offset;
    }
#endif

    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve, class Tfp>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>
    ::calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge
    (Tfunc_dispatch pfunc_dispatch,
     Tfunc_retrieve pfunc_retrieve,
     const S32 tag_max,
     ParticleSystem<Tfp> & psys,
     DomainInfo & dinfo,
     const S32 n_walk_limit,
     const bool clear_force,
     const bool reuse)
    {
        const F64 wtime_offset_0 = GetWtime();
        bool flag_retrieve = false;
        const S32 n_proc = Comm::getNumberOfProc();
        S32 ret = 0;
        if(!reuse){
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK0 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            setParticleLocalTreeImpl(psys, epi_org_, epj_org_, true); // new
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK1 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            setRootCell(dinfo); // original
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK2 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            mortonSortLocalTreeOnlyImpl(epi_org_, epj_org_, epi_sorted_, epj_sorted_loc_); // new
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK3 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            mortonSortFP<Tfp>(psys);    // new
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK4 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            linkCellLocalTreeOnly();   // original
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK5 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            F64 wtime_offset = GetWtime();
            CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(), epj_sorted_loc_.getPointer(), lev_max_loc_, n_leaf_limit_); //new
            time_profile_.calc_moment_local_tree += GetWtime() - wtime_offset;
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK6 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            addMomentAsSpLocalTreeImpl(typename TSM::force_type()); // original
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK7 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            exchangeLocalEssentialTree3(dinfo); // original
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK8 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            setLocalEssentialTreeToGlobalTree2(); // original
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK9 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            mortonSortGlobalTreeOnly3(); // new
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK10 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            linkCellGlobalTreeOnly(); // original
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK11 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            calcMomentGlobalTreeOnly(); // original
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK12 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            addMomentAsSpGlobalTreeImpl(typename TSM::force_type()); // original
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK13 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            makeIPGroup(); // original
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK14 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            makeAllInteractionListId3(true); // new
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK15 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            ret += calcForceUsingIdListMultiWalkIndex3(pfunc_dispatch, pfunc_retrieve, n_walk_limit, mw_info_, clear_force, reuse, flag_retrieve); // new
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK16 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
        }
        else{
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK17 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            F64 wtime_offset = GetWtime();
#ifdef SUNWAY
            epj_sorted_loc_.resizeNoInitialize(n_loc_tot_);
            unsigned long arg[5];
            arg[0] = (unsigned long)(n_loc_tot_);
            arg[1] = (unsigned long)(epi_sorted_.getPointer());
            arg[2] = (unsigned long)(epj_sorted_loc_.getPointer());
            arg[3] = (unsigned long)(sizeof(epj_sorted_loc_[0]));
            __real_athread_spawn((void*)slave_CopyDirect, arg);
            athread_join();
            //std::cerr<<"END Direct Copy"<<std::endl;
            //std::cerr<<"epj_sorted_loc_[0].id(B)= "<<epj_sorted_loc_[0].id<<std::endl;
            //std::cerr<<"epj_sorted_loc_[1].id(B)= "<<epj_sorted_loc_[1].id<<std::endl;
            //std::cerr<<"epj_sorted_loc_[2].id(B)= "<<epj_sorted_loc_[2].id<<std::endl;            
            //for(S32 i=0; i<n_loc_tot_; i++) assert(epj_sorted_loc_[i].id == epi_sorted_[i].id); // for debug

            /*
            // prefetch version
            for(S32 i=0; i<n_loc_tot_; i+=4){
                asm volatile ("fetchd     %0(%1)" : : "i"(L1ROFF), "r"(epi_sorted_.getPointer()+i));
                asm volatile ("fetchd_e   %0(%1)" : : "i"(L2ROFF), "r"(epi_sorted_.getPointer()+i));
                asm volatile ("fetchd_w   %0(%1)" : : "i"(L1ROFF), "r"(epj_sorted_loc_.getPointer()+i));
                asm volatile ("fetchd_we  %0(%1)" : : "i"(L2ROFF), "r"(epj_sorted_loc_.getPointer()+i));
                epj_sorted_loc_[i+0].copyFromEPI(epi_sorted_[i+0]);
                epj_sorted_loc_[i+1].copyFromEPI(epi_sorted_[i+1]);
                epj_sorted_loc_[i+2].copyFromEPI(epi_sorted_[i+2]);
                epj_sorted_loc_[i+3].copyFromEPI(epi_sorted_[i+3]);
            }
            */
#else
            for(S32 i=0; i<n_loc_tot_; i++){
                epj_sorted_loc_[i].copyFromEPI(epi_sorted_[i]);
            }
#endif
            time_profile_.set_particle_local_tree+= GetWtime() - wtime_offset;
            wtime_t0 = GetWtime() - wtime_offset;
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK18 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif

            /*
#ifndef SUNWAY
            wtime_offset = GetWtime();
            for(S32 i=0; i<n_loc_tot_; i++){
                const S32 adr = adr_epj_loc2glb_[i];
                epj_sorted_[adr] = epj_sorted_loc_[i];
            }
            time_profile_.morton_sort_global_tree += GetWtime() - wtime_offset;
#endif
            */
            
            wtime_t1 = GetWtime() - wtime_offset;
            
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK19 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            wtime_offset = GetWtime();
            CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(), epj_sorted_loc_.getPointer(), lev_max_loc_, n_leaf_limit_);
            time_profile_.calc_moment_local_tree += GetWtime() - wtime_offset;
            wtime_t2 = GetWtime() - wtime_offset;
            
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK20 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            wtime_offset = GetWtime();
            addMomentAsSpLocalTreeImpl(typename TSM::force_type()); //original
            wtime_t3 = GetWtime() - wtime_offset;

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK21 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            wtime_offset = GetWtime();
            exchangeLocalEssentialTree3(dinfo, true); //original
            wtime_t4 = GetWtime() - wtime_offset;

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK22 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            wtime_offset = GetWtime();
            const S32 n_epj_add = comm_table_.n_disp_ep_recv_[n_proc];
            const S32 n_spj_add = comm_table_.n_disp_sp_recv_[n_proc];
#ifdef SUNWAY
//#if 0
            //for(S32 i=0; i<n_loc_tot_; i++)  epj_sorted_[i].id = -i+1; // for debug
            arg[0] = (unsigned long)(n_loc_tot_);
            arg[1] = (unsigned long)(epj_sorted_loc_.getPointer());
            arg[2] = (unsigned long)(epj_sorted_.getPointer());
            arg[3] = (unsigned long)(sizeof(epj_sorted_loc_[0]));
            arg[4] = (unsigned long)(adr_epj_loc2glb_.getPointer());
            __real_athread_spawn((void*)slave_CopyIndirect, arg);
            athread_join();
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK23 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            /*
            for(S32 i=0; i<n_loc_tot_; i++){
                const S32 adr = adr_epj_loc2glb_[i];
                assert(epj_sorted_[adr].id == epj_sorted_loc_[i].id);
            }
            */
            arg[0] = (unsigned long)(n_epj_add);
            arg[1] = (unsigned long)(comm_table_.ep_recv_.getPointer());
            arg[2] = (unsigned long)(epj_sorted_.getPointer());
            arg[3] = (unsigned long)(sizeof(epj_sorted_[0]));
            arg[4] = (unsigned long)(adr_epj_buf2glb_.getPointer());
            __real_athread_spawn((void*)slave_CopyIndirect, arg);
            athread_join();
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK24 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            /*
            for(S32 i=0; i<n_epj_add; i++){
                const S32 adr = adr_epj_buf2glb_[i];
                assert(epj_sorted_[adr].id==comm_table_.ep_recv_[i].id);
            }
            */
            arg[0] = (unsigned long)(n_spj_add);
            arg[1] = (unsigned long)(comm_table_.sp_recv_.getPointer());
            arg[2] = (unsigned long)(spj_sorted_.getPointer());
            arg[3] = (unsigned long)(sizeof(spj_sorted_[0]));
            arg[4] = (unsigned long)(adr_spj_buf2glb_.getPointer());
            __real_athread_spawn((void*)slave_CopyIndirect, arg);
            athread_join();
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK25 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            /*
            for(S32 i=0; i<n_spj_add; i++){
                const S32 adr = adr_spj_buf2glb_[i];
                assert(spj_sorted_[adr].mass == comm_table_.sp_recv_[i].mass);
            } 
            */
#else
            for(S32 i=0; i<n_loc_tot_; i++){
                const S32 adr = adr_epj_loc2glb_[i];
                epj_sorted_[adr] = epj_sorted_loc_[i];
            } 
            for(S32 i=0; i<n_epj_add; i++){
                const S32 adr = adr_epj_buf2glb_[i];
                epj_sorted_[adr] = comm_table_.ep_recv_[i];
            }
            for(S32 i=0; i<n_spj_add; i++){
                const S32 adr = adr_spj_buf2glb_[i];
                spj_sorted_[adr] = comm_table_.sp_recv_[i];
            }
#endif
            time_profile_.morton_sort_global_tree += GetWtime() - wtime_offset;
            
            wtime_t5 = GetWtime() - wtime_offset;

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK26 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            wtime_offset = GetWtime();
            calcMomentGlobalTreeOnly(); //original
            wtime_t6 = GetWtime() - wtime_offset;

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK27 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            wtime_offset = GetWtime();
            addMomentAsSpGlobalTreeImpl(typename TSM::force_type()); //original
            wtime_t7 = GetWtime() - wtime_offset;

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK28 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            wtime_offset = GetWtime();
            ret += calcForceUsingIdListMultiWalkIndex3(pfunc_dispatch, pfunc_retrieve, n_walk_limit, mw_info_, clear_force, reuse, flag_retrieve); // new
            wtime_t8 = GetWtime() - wtime_offset;
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK29 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
        }
        time_profile_.calc_force_all += GetWtime() - wtime_offset_0;
        return ret;
    }
}
