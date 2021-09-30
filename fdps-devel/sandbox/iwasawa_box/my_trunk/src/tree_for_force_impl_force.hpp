#pragma once

#include"tree_walk.hpp"

namespace ParticleSimulator{
    ///////////////////////////////////////////////////
    //
    // FUNCTIONS OF WALK+FORCE WITH DOUBLE BUFFERING 
    //
    ///////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////
    //////////// Walk+Force, Kernel:Ptcl, List:Index ////////////
    // NO MULTI WALK VERSION
    // TYPE: SHORT
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForce(Tfunc_ep_ep pfunc_ep_ep,
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
	F64 offset_walk_tree,offset_dispatch;
        if(n_ipg > 0){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : ni_tmp, nj_tmp, n_interaction_ep_ep_tmp)
#endif
            for(S32 i=0; i<n_ipg; i++){
  	        offset_walk_tree = GetWtime();
	        makeInteractionList(i);
		time_profile_.calc_force__core__walk_tree += GetWtime() - offset_walk_tree;
                ni_tmp += ipg_[i].n_ptcl_;
                nj_tmp += epj_for_force_[Comm::getThreadNum()].size();
                n_interaction_ep_ep_tmp += ipg_[i].n_ptcl_ * epj_for_force_[Comm::getThreadNum()].size();
		offset_dispatch = GetWtime();
                calcForceOnly( pfunc_ep_ep, i, clear);
		time_profile_.calc_force__core__dispatch += GetWtime() - offset_dispatch;
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
        //std::cout<<"force_sorted_.size()="<<force_sorted_.size()<<std::endl;
        //std::cout<<"force_org_.size()="<<force_org_.size()<<std::endl;
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
        force_.resizeNoInitialize(n_loc_tot_);
        const S64 n_ipg = ipg_.size();
        n_walk_local_ += n_ipg;
        S64 ni_tmp = 0;
        S64 nj_tmp = 0;
        S64 n_interaction_ep_ep_tmp = 0;
        S64 n_interaction_ep_sp_tmp = 0;
        for(S32 i=0; i<Comm::getNumberOfThread(); i++) n_cell_open_[i] = 0;
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
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_tmp;
            n_interaction_ep_sp_local_ += n_interaction_ep_sp_tmp;
            for(S32 i=1; i<Comm::getNumberOfThread(); i++) n_cell_open_[0] += n_cell_open_[i];
        }
        else{
            ni_ave_ = nj_ave_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
        }
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
        PARTICLE_SIMULATOR_PRINT_LINE_INFO();
        std::cout<<"ipg_.size()="<<ipg_.size()<<std::endl;
#endif
        time_profile_.calc_force += GetWtime() - time_offset;
    }

    // FUNCTIONS CALLED IN calcForce()
    // ORIGINAL VERSION
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionList(const S32 adr_ipg, const bool clear){
        makeInteractionListImpl(typename TSM::search_type(), adr_ipg, clear);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchLong,
                            const S32 adr_ipg,
                            const bool clear){
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
#if 0
        F64 mass_tmp = 0.0;
        for(S32 jp=0; jp<epj_for_force_[ith].size(); jp++)  mass_tmp += epj_for_force_[ith][jp].mass;
        for(S32 jp=0; jp<spj_for_force_[ith].size(); jp++)  mass_tmp += spj_for_force_[ith][jp].mass;
        assert(mass_tmp == 1.0);
        if(Comm::getRank()==0){
            std::cout<<"makeInteractionListImpl) mass_tmp= "<<mass_tmp
                     <<" n_epj= "<<epj_for_force_[ith].size()
                     <<" n_spj= "<<spj_for_force_[ith].size()
                     <<std::endl;
            std::cout<<"adr_ipg= "<<adr_ipg
                     <<" ipg_[adr_ipg].n_ptcl_= "<<ipg_[adr_ipg].n_ptcl_
                     <<" vertex= "<<(ipg_[adr_ipg]).vertex_
                     <<std::endl;
        }

#endif        
    }



    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchLongCutoff, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_;
        F64 r_crit_sq = LARGE_FLOAT;
        if(theta_ > 0.0){
            r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
        }
        const F64 r_cut_sq  = epj_sorted_[0].getRSearch() * epj_sorted_[0].getRSearch();
        //std::cerr<<"r_crit= "<<sqrt(r_crit_sq)<<std::endl;
	if(clear){
	    epj_for_force_[ith].clearSize();
	    spj_for_force_[ith].clearSize();
	}
        const F64ort cell_box = pos_root_cell_;
        if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
            //std::cerr<<"check a"<<std::endl;
            MakeInteractionListLongCutoffEPSP
                (tc_glb_, tc_glb_[0].adr_tc_, tp_glb_, 
                 epj_sorted_, epj_for_force_[ith],
                 spj_sorted_, spj_for_force_[ith],
                 cell_box,
                 pos_target_box, r_crit_sq, r_cut_sq, n_leaf_limit_); 
        }
        else{
            //std::cerr<<"check b"<<std::endl;
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


    // FOR P^3T
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchLongScatter, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_;
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
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
            else{
                spj_for_force_[ith].increaseSize();
                spj_for_force_[ith].back().copyFromMoment(tc_glb_[0].mom_);
            }
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchLongSymmetry, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box_in = ipg_[adr_ipg].vertex_;
        const F64ort pos_target_box_out = ipg_[adr_ipg].vertex_out_;
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_) * 0.25;
	if(clear){	
	    epj_for_force_[ith].clearSize();
	    spj_for_force_[ith].clearSize();
	}
        if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
            MakeInteractionListLongSymmetryEPSP
                (tc_glb_, tc_glb_[0].adr_tc_, 
                 tp_glb_, epj_sorted_, 
                 epj_for_force_[ith],
                 spj_sorted_, spj_for_force_[ith],
                 pos_target_box_in, pos_target_box_out, r_crit_sq, n_leaf_limit_);
        }
        else{
            const F64vec pos_tmp = tc_glb_[0].mom_.getPos();
            if( pos_target_box_in.getDistanceMinSQ(pos_tmp) <= r_crit_sq*4.0 ){
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

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchShortScatter, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_;
	if(clear){	
	    epj_for_force_[ith].clearSize();
	}
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
                    epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp]);
                    const F64vec pos_new = epj_for_force_[ith].back().getPos();
                    epj_for_force_[ith].back().setPos(pos_new);
                }
            }
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchShortGather, const S32 adr_ipg, const bool clear){
        const S32 ith = Comm::getThreadNum();
        const F64ort pos_target_box = (ipg_[adr_ipg]).vertex_;
	if(clear){
	    epj_for_force_[ith].clearSize();
	}
        if( !tc_glb_[0].isLeaf(n_leaf_limit_) ){
            MakeListUsingInnerBoundary
                (tc_glb_.getPointer(),     tc_glb_[0].adr_tc_,
                 epj_sorted_.getPointer(), epj_for_force_[ith],
                 pos_target_box,                 n_leaf_limit_);
        }
        else{
            if( pos_target_box.overlapped( tc_glb_[0].mom_.getVertexIn()) ){
                S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                const S32 n_tmp = tc_glb_[0].n_ptcl_;
                epj_for_force_[ith].reserveEmptyAreaAtLeast( n_tmp );
                for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
                    const F64vec pos_tmp = epj_sorted_[adr_ptcl_tmp].getPos();
                    if( pos_target_box.overlapped( pos_tmp) ){
                        epj_for_force_[ith].pushBackNoCheck(epj_sorted_[adr_ptcl_tmp]);
                        const F64vec pos_new = epj_for_force_[ith].back().getPos();
                        epj_for_force_[ith].back().setPos(pos_new);
                    }
                }
            }
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeInteractionListImpl(TagSearchShortSymmetry, const S32 adr_ipg, const bool clear){
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
            if( pos_target_box_out.contained(tc_glb_[0].mom_.getVertexIn()) 
                || pos_target_box_in.contained(tc_glb_[9].mom_.getVertexOut()) ){
                S32 adr_ptcl_tmp = tc_glb_[0].adr_ptcl_;
                const S32 n_tmp = tc_glb_[0].n_ptcl_;
                epj_for_force_[ith].reserveAtLeast( n_tmp );
                for(S32 ip=0; ip<n_tmp; ip++, adr_ptcl_tmp++){
                    const F64vec pos_tmp = epj_sorted_[adr_ptcl_tmp].getPos();
                    const F64 size_tmp = epj_sorted_[adr_ptcl_tmp].getRSearch();
                    const F64 dis_sq_tmp = pos_target_box_in.getDistanceMinSQ(pos_tmp);
		    if( pos_target_box_out.notContained(pos_tmp) && dis_sq_tmp > size_tmp*size_tmp) continue;
                    epj_for_force_[ith].pushBackNoCheck( epj_sorted_[adr_ptcl_tmp] );
                    const F64vec pos_new = epj_for_force_[ith].back().getPos();
                    epj_for_force_[ith].back().setPos(pos_new);
                }
            }
        }
#endif
    }


    
    
    //////////////////////////////////////////////////////////////
    //////////// Walk+Force, Kernel:Index, List:Index ////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve, class Tsys>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkIndex2(Tfunc_dispatch pfunc_dispatch,
                             Tfunc_retrieve pfunc_retrieve,
                             Tsys & sys,
                             const S32 tag_max,
                             const S32 n_walk_limit,
                             const bool flag_keep_list,
                             const bool clear){
        const F64 wtime_offset = GetWtime();
        if(tag_max <= 0){
            PARTICLE_SIMULATOR_PRINT_ERROR("tag_max is illegal. In currente version, tag_max must be 1");
            Abort(-1);
        }
        force_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                force_[i].clear();
            }
        }
        S32 ret = 0;
        ret = calcForceMultiWalkIndex2Impl(typename TSM::force_type(),
                                           pfunc_dispatch,
                                           pfunc_retrieve,
                                           sys,
                                           tag_max,
                                           n_walk_limit,
                                           flag_keep_list,
                                           clear);
        time_profile_.calc_force += GetWtime() - wtime_offset;
        return ret;
    }

    //////////// Walk+Force, Kernel:Index, List:Index, Force:Long //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
        template<class Tfunc_dispatch, class Tfunc_retrieve, class Tsys>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkIndex2Impl(TagForceLong,
                                 Tfunc_dispatch pfunc_dispatch,
                                 Tfunc_retrieve pfunc_retrieve,
                                 Tsys & sys,
                                 const S32 tag_max,
                                 const S32 n_walk_limit,
                                 const bool flag_keep_list,
                                 const bool clear){
        const F64 offset_core = GetWtime();
        // send all epj and spj
        Tepi ** epi_dummy = NULL;
        S32 * n_epi_dummy = NULL;
        S32 ** id_epj_dummy = NULL;
        S32 *  n_epj_dummy = NULL;
        S32 ** id_spj_dummy = NULL;
        S32 *  n_spj_dummy = NULL;
        pfunc_dispatch(0, 0, (const Tepi**)epi_dummy, n_epi_dummy,
                       (const S32**)id_epj_dummy, n_epj_dummy,
                       (const S32**)id_spj_dummy, n_spj_dummy,
                       epj_sorted_.getPointer(), epj_sorted_.size(),
                       spj_sorted_.getPointer(), spj_sorted_.size(),
                       true);
        static bool first = true;
        S32 ret = 0;
        S32 tag = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        force_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        if(n_ipg <= 0) return 0;
        const S32 n_ipg_amari = (n_ipg > 0) ? n_ipg%n_walk_limit : 0;
        n_walk_local_ += n_ipg;
        if(flag_keep_list){
            interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
            interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
            interaction_list_.adr_ep_.clearSize();
            interaction_list_.n_sp_.resizeNoInitialize(n_ipg);
            interaction_list_.n_disp_sp_.resizeNoInitialize(n_ipg+1);
            interaction_list_.adr_sp_.clearSize();
        }
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg_amari==0 ? 0 : 1);

        static std::vector<S32> n_walk_ar;
        n_walk_ar.resize(n_loop_max); // group of walk (walk_grp[i] < n_walk_limit)
        static std::vector<S32> n_disp_walk_ar;
        n_disp_walk_ar.resize(n_loop_max+1);

        static S32  * iw2ith;
        static S32  * iw2cnt;
        static S32 ** n_epj_disp_thread; // [n_thread][n_walk+1]
        static S32 ** n_spj_disp_thread;// [n_thread][n_walk+1]
        static Tforce ** force_array; // array of pointer *[n_walk]
        static Tforce ** force_prev_array; // array of pointer *[n_walk]
        static S32  * cnt_thread;
        static S32  * n_ep_cum_thread;
        static S32  * n_sp_cum_thread;
        static S64 * n_interaction_ep_ep_ar;
        static S64 * n_interaction_ep_sp_ar;

        static ReallocatableArray<S32> * adr_epj_tmp;
        static ReallocatableArray<S32> * adr_spj_tmp;
        static ReallocatableArray<S32> * adr_ipg_tmp;

        static ReallocatableArray<S32> n_epi_ar;
        n_epi_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32> n_epi_ar_prev;
        n_epi_ar_prev.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tepi*> epi_ar; // array of pointer *[n_walk]
        epi_ar.resizeNoInitialize(n_walk_limit);

        static ReallocatableArray<S32> n_epj_ar;
        n_epj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32*> id_epj_ar;
        id_epj_ar.resizeNoInitialize(n_walk_limit);

        static ReallocatableArray<S32> n_spj_ar;
        n_spj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32*> id_spj_ar;
        id_spj_ar.resizeNoInitialize(n_walk_limit);

        if(first){
            iw2ith = new S32[n_walk_limit];
            iw2cnt = new S32[n_walk_limit];
            n_epj_disp_thread = new S32*[n_thread];
            n_spj_disp_thread = new S32*[n_thread];
            force_array      = new Tforce*[n_walk_limit];
            force_prev_array = new Tforce*[n_walk_limit];
            for(int i=0; i<n_thread; i++){
                n_epj_disp_thread[i] = new S32[n_walk_limit+1];
                n_spj_disp_thread[i] = new S32[n_walk_limit+1];
            }
            cnt_thread = new S32[n_thread];
            n_ep_cum_thread = new S32[n_thread];
            n_sp_cum_thread = new S32[n_thread];
            n_interaction_ep_ep_ar = new S64[n_thread];
            n_interaction_ep_sp_ar = new S64[n_thread];

            adr_epj_tmp = new ReallocatableArray<S32>[n_thread];
            adr_spj_tmp = new ReallocatableArray<S32>[n_thread];
            adr_ipg_tmp = new ReallocatableArray<S32>[n_thread];
            first = false;
        }
        n_disp_walk_ar[0] = 0;
        const S32 wg_tmp = n_ipg > 0 ? n_ipg%n_loop_max : 0;
        for(int wg=0; wg<wg_tmp; wg++){
            n_walk_ar[wg]        = n_ipg / n_loop_max + 1;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int wg=wg_tmp; wg<n_loop_max; wg++){
            n_walk_ar[wg]        = n_ipg / n_loop_max;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_ar[i] = n_interaction_ep_sp_ar[i] = 0;
        }
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
        const S32 adr_tree_sp_first = spj_org_.size();
        bool first_loop = true;
        S32 n_walk_prev = 0;
        CountT n_epj_for_force_tmp = 0;
        CountT n_spj_for_force_tmp = 0;
        if(n_ipg > 0){
            S32 ip_head = 0, ip_tail=0, ip_head_prev=0, ip_tail_prev=0;
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = n_walk_ar[wg];
                const S32 walk_grp_head = n_disp_walk_ar[wg];
                ip_head = ipg_[walk_grp_head].adr_ptcl_; 
                ip_tail = ipg_[walk_grp_head+n_walk-1].adr_ptcl_ + ipg_[walk_grp_head+n_walk-1].n_ptcl_;
                F64 offset_calc_force__core__walk_tree = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel reduction(+ : n_epj_for_force_tmp, n_spj_for_force_tmp)
#endif
                {
                    const S32 ith = Comm::getThreadNum();
                    n_ep_cum_thread[ith] = n_sp_cum_thread[ith] = cnt_thread[ith] = 0;
                    n_epj_disp_thread[ith][0] = 0;
                    n_spj_disp_thread[ith][0] = 0;
                    adr_epj_tmp[ith].clearSize();
                    adr_spj_tmp[ith].clearSize();
                    adr_ipg_tmp[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        adr_ipg_tmp[ith].push_back(id_ipg);
                        const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_; 
                        const S32 ith = Comm::getThreadNum();
                        n_epi_ar[iw] = ipg_[id_ipg].n_ptcl_;
                        epi_ar[iw]   = epi_sorted_.getPointer(first_adr_ip);
                        //force_array[iw] = force_sorted_.getPointer(first_adr_ip);
                        force_array[iw] = force_.getPointer(first_adr_ip);

                        TargetBox<TSM> target_box;
                        GetTargetBox<TSM>(ipg_[id_ipg], target_box);
                        S32 adr_tc = 0;
                        MakeListUsingTreeRecursiveTop
                            <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, 
                             Tspj, TagWalkModeNormal, TagChopLeafTrue>
                            (tc_glb_,  adr_tc, tp_glb_,
                             epj_sorted_, adr_epj_tmp[ith],
                             spj_sorted_, adr_spj_tmp[ith],
                             target_box,
                             r_crit_sq, n_leaf_limit_,
                             adr_tree_sp_first, F64vec(0.0));

                        n_epj_ar[iw] = adr_epj_tmp[ith].size() - n_ep_cum_thread[ith];
                        n_spj_ar[iw] = adr_spj_tmp[ith].size() - n_sp_cum_thread[ith];

                        n_ep_cum_thread[ith] = adr_epj_tmp[ith].size();
                        n_sp_cum_thread[ith] = adr_spj_tmp[ith].size();

                        n_epj_disp_thread[ith][cnt_thread[ith]+1] = n_ep_cum_thread[ith];
                        n_spj_disp_thread[ith][cnt_thread[ith]+1] = n_sp_cum_thread[ith];

                        n_interaction_ep_ep_ar[ith] += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[iw]);
                        n_interaction_ep_sp_ar[ith] += ((S64)n_spj_ar[iw]*(S64)n_epi_ar[iw]);

                        n_epj_for_force_tmp += (CountT)n_epj_ar[iw];
                        n_spj_for_force_tmp += (CountT)n_spj_ar[iw];

                        iw2ith[iw] = ith;
                        iw2cnt[iw] = cnt_thread[ith];
                        cnt_thread[ith]++;
                    } // end of OMP for
                } // end of OMP parallel scope
                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;
                if(flag_keep_list){
                    interaction_list_.n_disp_ep_[0] = interaction_list_.n_disp_sp_[0] = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_ep_[id_ipg] = n_epj_ar[iw];
                        interaction_list_.n_sp_[id_ipg] = n_spj_ar[iw];
                    }
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_disp_ep_[id_ipg+1] = interaction_list_.n_disp_ep_[id_ipg] + interaction_list_.n_ep_[id_ipg];
                        interaction_list_.n_disp_sp_[id_ipg+1] = interaction_list_.n_disp_sp_[id_ipg] + interaction_list_.n_sp_[id_ipg];
                    }
                    interaction_list_.adr_ep_.resizeNoInitialize( interaction_list_.n_disp_ep_[walk_grp_head+n_walk] );
                    interaction_list_.adr_sp_.resizeNoInitialize( interaction_list_.n_disp_sp_[walk_grp_head+n_walk] );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(S32 i=0; i<n_thread; i++){
                        for(S32 j=0; j<adr_ipg_tmp[i].size(); j++){
                            const S32 adr_ipg = adr_ipg_tmp[i][j];
                            S32 adr_ep = interaction_list_.n_disp_ep_[adr_ipg];
                            const S32 k_ep_h = n_epj_disp_thread[i][j];
                            const S32 k_ep_e = n_epj_disp_thread[i][j+1];
                            for(S32 k=k_ep_h; k<k_ep_e; k++, adr_ep++){
                                interaction_list_.adr_ep_[adr_ep] = adr_epj_tmp[i][k];
                            }
                            S32 adr_sp = interaction_list_.n_disp_sp_[adr_ipg];
                            const S32 k_sp_h = n_spj_disp_thread[i][j];
                            const S32 k_sp_e = n_spj_disp_thread[i][j+1];
                            for(S32 k=k_sp_h; k<k_sp_e; k++, adr_sp++){
                                interaction_list_.adr_sp_[adr_sp] = adr_spj_tmp[i][k];
                            }
                        }
                    }
                } // end of flag_keep_list

                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    S32 n_sp_head = n_spj_disp_thread[ith][cnt];
                    id_epj_ar[iw] = adr_epj_tmp[ith].getPointer(n_ep_head);
                    id_spj_ar[iw] = adr_spj_tmp[ith].getPointer(n_sp_head);
                }
                //time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;

                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar_prev.getPointer(), force_prev_array);
                } // retrieve                
                pfunc_dispatch(0, n_walk, 
                               (const Tepi**)epi_ar.getPointer(),   n_epi_ar.getPointer(),
                               (const S32**)id_epj_ar.getPointer(), n_epj_ar.getPointer(),
                               (const S32**)id_spj_ar.getPointer(), n_spj_ar.getPointer(),
                               epj_sorted_.getPointer(), epj_sorted_.size(),
                               spj_sorted_.getPointer(), spj_sorted_.size(),
                               false);
                if(!first_loop){
                    F64 wtime_offset_wb_int_cp = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for
#endif
                    for(S32 ip=ip_head_prev; ip<ip_tail_prev; ip++){
                        sys[ip].copyFromForce(force_[ip]);
                        sys[ip].integrate();
                        epi_org_[ip].copyFromFP( sys[ip] );
                        epj_org_[ip].copyFromFP( sys[ip] );
                    }
                    time_profile_.calc_force__wb_int_cp += GetWtime() - wtime_offset_wb_int_cp;
                }
                for(int iw=0; iw<n_walk; iw++){
                    n_epi_ar_prev[iw] = n_epi_ar[iw];
                    force_prev_array[iw] = force_array[iw];
                }
                n_walk_prev = n_walk;
                ip_head_prev = ip_head;
                ip_tail_prev = ip_tail;
                first_loop = false;
            } // end of walk group loop
            ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar_prev.getPointer(), force_prev_array);

            F64 wtime_offset_wb_int_cp = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for
#endif
            for(S32 ip=ip_head_prev; ip<ip_tail_prev; ip++){
                sys[ip].copyFromForce(force_[ip]);
                sys[ip].integrate();
                epi_org_[ip].copyFromFP( sys[ip] );
                epj_org_[ip].copyFromFP( sys[ip] );
            }
            time_profile_.calc_force__wb_int_cp += GetWtime() - wtime_offset_wb_int_cp;
        } // if(n_ipg > 0)
        else{
            ni_ave_ = nj_ave_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
            n_epj_for_force_ = n_spj_for_force_ = 0;
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_ar[i];
            n_interaction_ep_sp_local_ += n_interaction_ep_sp_ar[i];
        }
        n_epj_for_force_ += n_epj_for_force_tmp;
        n_spj_for_force_ += n_spj_for_force_tmp;
        time_profile_.calc_force__core += GetWtime() - offset_core;
        return ret;
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                            Tfunc_retrieve pfunc_retrieve,
                            const S32 tag_max,
                            const S32 n_walk_limit,
                            const bool flag_keep_list,
                            const bool clear){
        const F64 wtime_offset = GetWtime();
        if(tag_max <= 0){
            PARTICLE_SIMULATOR_PRINT_ERROR("tag_max is illegal. In currente version, tag_max must be 1");
            Abort(-1);
        }
        force_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                force_[i].clear();
            }
        }
        S32 ret = 0;
        ret = calcForceMultiWalkIndexImpl(typename TSM::force_type(),
                                          pfunc_dispatch,
                                          pfunc_retrieve,
                                          tag_max,
                                          n_walk_limit,
                                          flag_keep_list,
                                          clear);
        time_profile_.calc_force += GetWtime() - wtime_offset;
        return ret;
    }

    //////////// Walk+Force, Kernel:Index, List:Index, Force:Long //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkIndexImpl(TagForceLong,
                                Tfunc_dispatch pfunc_dispatch,
                                Tfunc_retrieve pfunc_retrieve,
                                const S32 tag_max,
                                const S32 n_walk_limit,
                                const bool flag_keep_list,
                                const bool clear){
        const F64 offset_core = GetWtime();
        // send all epj and spj
        Tepi ** epi_dummy = NULL;
        S32 * n_epi_dummy = NULL;
        S32 ** id_epj_dummy = NULL;
        S32 *  n_epj_dummy = NULL;
        S32 ** id_spj_dummy = NULL;
        S32 *  n_spj_dummy = NULL;
        pfunc_dispatch(0, 0, (const Tepi**)epi_dummy, n_epi_dummy,
                       (const S32**)id_epj_dummy, n_epj_dummy,
                       (const S32**)id_spj_dummy, n_spj_dummy,
                       epj_sorted_.getPointer(), epj_sorted_.size(),
                       spj_sorted_.getPointer(), spj_sorted_.size(),
                       true);
        static bool first = true;
        S32 ret = 0;
        S32 tag = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        force_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        if(n_ipg <= 0) return 0;
        const S32 n_ipg_amari = (n_ipg > 0) ? n_ipg%n_walk_limit : 0;
        n_walk_local_ += n_ipg;
        if(flag_keep_list){
            interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
            interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
            interaction_list_.adr_ep_.clearSize();
            interaction_list_.n_sp_.resizeNoInitialize(n_ipg);
            interaction_list_.n_disp_sp_.resizeNoInitialize(n_ipg+1);
            interaction_list_.adr_sp_.clearSize();
        }
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg_amari==0 ? 0 : 1);

        static std::vector<S32> n_walk_ar;
        n_walk_ar.resize(n_loop_max); // group of walk (walk_grp[i] < n_walk_limit)
        static std::vector<S32> n_disp_walk_ar;
        n_disp_walk_ar.resize(n_loop_max+1);

        static S32  * iw2ith;
        static S32  * iw2cnt;
        static S32 ** n_epj_disp_thread; // [n_thread][n_walk+1]
        static S32 ** n_spj_disp_thread;// [n_thread][n_walk+1]
        static Tforce ** force_array; // array of pointer *[n_walk]
        static Tforce ** force_prev_array; // array of pointer *[n_walk]
        static S32  * cnt_thread;
        static S32  * n_ep_cum_thread;
        static S32  * n_sp_cum_thread;
        static S64 * n_interaction_ep_ep_ar;
        static S64 * n_interaction_ep_sp_ar;

        static ReallocatableArray<S32> * adr_epj_tmp;
        static ReallocatableArray<S32> * adr_spj_tmp;
        static ReallocatableArray<S32> * adr_ipg_tmp;

        static ReallocatableArray<S32> n_epi_ar;
        n_epi_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32> n_epi_ar_prev;
        n_epi_ar_prev.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tepi*> epi_ar; // array of pointer *[n_walk]
        epi_ar.resizeNoInitialize(n_walk_limit);

        static ReallocatableArray<S32> n_epj_ar;
        n_epj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32*> id_epj_ar;
        id_epj_ar.resizeNoInitialize(n_walk_limit);

        static ReallocatableArray<S32> n_spj_ar;
        n_spj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32*> id_spj_ar;
        id_spj_ar.resizeNoInitialize(n_walk_limit);

        if(first){
            iw2ith = new S32[n_walk_limit];
            iw2cnt = new S32[n_walk_limit];
            n_epj_disp_thread = new S32*[n_thread];
            n_spj_disp_thread = new S32*[n_thread];
            force_array      = new Tforce*[n_walk_limit];
            force_prev_array = new Tforce*[n_walk_limit];
            for(int i=0; i<n_thread; i++){
                n_epj_disp_thread[i] = new S32[n_walk_limit+1];
                n_spj_disp_thread[i] = new S32[n_walk_limit+1];
            }
            cnt_thread = new S32[n_thread];
            n_ep_cum_thread = new S32[n_thread];
            n_sp_cum_thread = new S32[n_thread];
            n_interaction_ep_ep_ar = new S64[n_thread];
            n_interaction_ep_sp_ar = new S64[n_thread];

            adr_epj_tmp = new ReallocatableArray<S32>[n_thread];
            adr_spj_tmp = new ReallocatableArray<S32>[n_thread];
            adr_ipg_tmp = new ReallocatableArray<S32>[n_thread];
            first = false;
        }
        n_disp_walk_ar[0] = 0;
        const S32 wg_tmp = n_ipg > 0 ? n_ipg%n_loop_max : 0;
        for(int wg=0; wg<wg_tmp; wg++){
            n_walk_ar[wg]        = n_ipg / n_loop_max + 1;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int wg=wg_tmp; wg<n_loop_max; wg++){
            n_walk_ar[wg]        = n_ipg / n_loop_max;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_ar[i] = n_interaction_ep_sp_ar[i] = 0;
        }
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
        const S32 adr_tree_sp_first = spj_org_.size();
        bool first_loop = true;
        S32 n_walk_prev = 0;
        CountT n_epj_for_force_tmp = 0;
        CountT n_spj_for_force_tmp = 0;
        //std::cerr<<"n_loop_max= "<<n_loop_max<<std::endl;
        if(n_ipg > 0){
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = n_walk_ar[wg];
                const S32 walk_grp_head = n_disp_walk_ar[wg];
                F64 offset_calc_force__core__walk_tree = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel reduction(+ : n_epj_for_force_tmp, n_spj_for_force_tmp)
#endif
                {
                    const S32 ith = Comm::getThreadNum();
                    n_ep_cum_thread[ith] = n_sp_cum_thread[ith] = cnt_thread[ith] = 0;
                    n_epj_disp_thread[ith][0] = 0;
                    n_spj_disp_thread[ith][0] = 0;
                    adr_epj_tmp[ith].clearSize();
                    adr_spj_tmp[ith].clearSize();
                    adr_ipg_tmp[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        adr_ipg_tmp[ith].push_back(id_ipg);
                        const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_; 
                        const S32 ith = Comm::getThreadNum();
                        n_epi_ar[iw] = ipg_[id_ipg].n_ptcl_;
                        epi_ar[iw]   = epi_sorted_.getPointer(first_adr_ip);
                        //force_array[iw] = force_sorted_.getPointer(first_adr_ip);
                        force_array[iw] = force_.getPointer(first_adr_ip);

                        TargetBox<TSM> target_box;
                        GetTargetBox<TSM>(ipg_[id_ipg], target_box);
                        S32 adr_tc = 0;
                        MakeListUsingTreeRecursiveTop
                            <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, 
                             Tspj, TagWalkModeNormal, TagChopLeafTrue>
                            (tc_glb_,  adr_tc, tp_glb_,
                             epj_sorted_, adr_epj_tmp[ith],
                             spj_sorted_, adr_spj_tmp[ith],
                             target_box,
                             r_crit_sq, n_leaf_limit_,
                             adr_tree_sp_first, F64vec(0.0));

                        n_epj_ar[iw] = adr_epj_tmp[ith].size() - n_ep_cum_thread[ith];
                        n_spj_ar[iw] = adr_spj_tmp[ith].size() - n_sp_cum_thread[ith];

                        n_ep_cum_thread[ith] = adr_epj_tmp[ith].size();
                        n_sp_cum_thread[ith] = adr_spj_tmp[ith].size();

                        n_epj_disp_thread[ith][cnt_thread[ith]+1] = n_ep_cum_thread[ith];
                        n_spj_disp_thread[ith][cnt_thread[ith]+1] = n_sp_cum_thread[ith];

                        n_interaction_ep_ep_ar[ith] += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[iw]);
                        n_interaction_ep_sp_ar[ith] += ((S64)n_spj_ar[iw]*(S64)n_epi_ar[iw]);

                        n_epj_for_force_tmp += (CountT)n_epj_ar[iw];
                        n_spj_for_force_tmp += (CountT)n_spj_ar[iw];

                        iw2ith[iw] = ith;
                        iw2cnt[iw] = cnt_thread[ith];
                        cnt_thread[ith]++;
                    } // end of OMP for
                } // end of OMP parallel scope
                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;
                if(flag_keep_list){
                    interaction_list_.n_disp_ep_[0] = interaction_list_.n_disp_sp_[0] = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_ep_[id_ipg] = n_epj_ar[iw];
                        interaction_list_.n_sp_[id_ipg] = n_spj_ar[iw];
                    }
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_disp_ep_[id_ipg+1] = interaction_list_.n_disp_ep_[id_ipg] + interaction_list_.n_ep_[id_ipg];
                        interaction_list_.n_disp_sp_[id_ipg+1] = interaction_list_.n_disp_sp_[id_ipg] + interaction_list_.n_sp_[id_ipg];
                    }
                    interaction_list_.adr_ep_.resizeNoInitialize( interaction_list_.n_disp_ep_[walk_grp_head+n_walk] );
                    interaction_list_.adr_sp_.resizeNoInitialize( interaction_list_.n_disp_sp_[walk_grp_head+n_walk] );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(S32 i=0; i<n_thread; i++){
                        for(S32 j=0; j<adr_ipg_tmp[i].size(); j++){
                            const S32 adr_ipg = adr_ipg_tmp[i][j];
                            S32 adr_ep = interaction_list_.n_disp_ep_[adr_ipg];
                            const S32 k_ep_h = n_epj_disp_thread[i][j];
                            const S32 k_ep_e = n_epj_disp_thread[i][j+1];
                            for(S32 k=k_ep_h; k<k_ep_e; k++, adr_ep++){
                                interaction_list_.adr_ep_[adr_ep] = adr_epj_tmp[i][k];
                            }
                            S32 adr_sp = interaction_list_.n_disp_sp_[adr_ipg];
                            const S32 k_sp_h = n_spj_disp_thread[i][j];
                            const S32 k_sp_e = n_spj_disp_thread[i][j+1];
                            for(S32 k=k_sp_h; k<k_sp_e; k++, adr_sp++){
                                interaction_list_.adr_sp_[adr_sp] = adr_spj_tmp[i][k];
                            }
                        }
                    }
                } // end of flag_keep_list

                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    S32 n_sp_head = n_spj_disp_thread[ith][cnt];
                    id_epj_ar[iw] = adr_epj_tmp[ith].getPointer(n_ep_head);
                    id_spj_ar[iw] = adr_spj_tmp[ith].getPointer(n_sp_head);
                }
                //time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;

                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar_prev.getPointer(), force_prev_array);
                } // retrieve                
                pfunc_dispatch(0, n_walk, 
                               (const Tepi**)epi_ar.getPointer(),   n_epi_ar.getPointer(),
                               (const S32**)id_epj_ar.getPointer(), n_epj_ar.getPointer(),
                               (const S32**)id_spj_ar.getPointer(), n_spj_ar.getPointer(),
                               epj_sorted_.getPointer(), epj_sorted_.size(),
                               spj_sorted_.getPointer(), spj_sorted_.size(),
                               false);
                first_loop = false;
                for(int iw=0; iw<n_walk; iw++){
                    n_epi_ar_prev[iw] = n_epi_ar[iw];
                    force_prev_array[iw] = force_array[iw];
                }
                n_walk_prev = n_walk;
            } // end of walk group loop
            ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar_prev.getPointer(), force_prev_array);
        } // if(n_ipg > 0)
        else{
            ni_ave_ = nj_ave_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
            n_epj_for_force_ = n_spj_for_force_ = 0;
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_ar[i];
            n_interaction_ep_sp_local_ += n_interaction_ep_sp_ar[i];
        }
        n_epj_for_force_ += n_epj_for_force_tmp;
        n_spj_for_force_ += n_spj_for_force_tmp;
        time_profile_.calc_force__core += GetWtime() - offset_core;
        return ret;
    }

    //////////// Walk+Force, Kernel:Index, List:Index, Force:Short //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkIndexImpl(TagForceShort,
                                Tfunc_dispatch pfunc_dispatch,
                                Tfunc_retrieve pfunc_retrieve,
                                const S32 tag_max,
                                const S32 n_walk_limit,
                                const bool flag_keep_list,
                                const bool clear){
        const F64 offset_core = GetWtime();
        Tepi ** epi_dummy = NULL;
        S32 * n_epi_dummy = NULL;
        S32 ** id_epj_dummy = NULL;
        S32 *  n_epj_dummy = NULL;
        pfunc_dispatch(0, 0, (const Tepi**)epi_dummy, n_epi_dummy,
                       (const S32**)id_epj_dummy, n_epj_dummy,
                       epj_sorted_.getPointer(), epj_sorted_.size(),
                       true);
        static bool first = true;
        S32 ret = 0;
        S32 tag = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        force_.resizeNoInitialize(n_loc_tot_);
        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        if(n_ipg <= 0) return 0;

        const S32 n_ipg_amari = (n_ipg > 0) ? n_ipg%n_walk_limit : 0;
        n_walk_local_ += n_ipg;

        if(flag_keep_list){
            interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
            interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
            interaction_list_.adr_ep_.clearSize();
        }

        //const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg_amari==0 ? 0 : 1);
        static std::vector<S32> n_walk_ar;
        n_walk_ar.resize(n_loop_max); // group of walk (walk_grp[i] < n_walk_limit)
        static std::vector<S32> n_disp_walk_ar;
        n_disp_walk_ar.resize(n_loop_max+1);

        static S32  * iw2ith;
        static S32  * iw2cnt;
        static S32 ** n_epj_disp_thread; // [n_thread][n_walk]
        static Tforce ** force_array; // array of pointer *[n_walk]
        static Tforce ** force_prev_array; // array of pointer *[n_walk]
        static S32  * cnt_thread;
        static S32  * n_ep_cum_thread;
        static S64 * n_interaction_ep_ep_ar;

        static ReallocatableArray<S32> * adr_epj_tmp;
        static ReallocatableArray<S32> * adr_ipg_tmp;

        static ReallocatableArray<S32> n_epi_ar;
        n_epi_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32> n_epi_ar_prev;
        n_epi_ar_prev.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tepi*> epi_ar; // array of pointer *[n_walk]
        epi_ar.resizeNoInitialize(n_walk_limit);

        static ReallocatableArray<S32> n_epj_ar;
        n_epj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32*> id_epj_ar;
        id_epj_ar.resizeNoInitialize(n_walk_limit);

        if(first){
            iw2ith = new S32[n_walk_limit];
            iw2cnt = new S32[n_walk_limit];
            n_epj_disp_thread = new S32*[n_thread];
            force_array      = new Tforce*[n_walk_limit];
            force_prev_array = new Tforce*[n_walk_limit];
            for(int i=0; i<n_thread; i++){
                n_epj_disp_thread[i] = new S32[n_walk_limit+1];
            }
            cnt_thread = new S32[n_thread];
            n_ep_cum_thread = new S32[n_thread];
            n_interaction_ep_ep_ar = new S64[n_thread];

            adr_epj_tmp = new ReallocatableArray<S32>[n_thread];
            adr_ipg_tmp = new ReallocatableArray<S32>[n_thread];

            first = false;
        }
        n_disp_walk_ar[0] = 0;
        const S32 wg_tmp = n_ipg > 0 ? n_ipg%n_loop_max : 0;
        //for(int wg=0; wg<n_ipg%n_loop_max; wg++){
        for(int wg=0; wg<wg_tmp; wg++){
            n_walk_ar[wg]        = n_ipg / n_loop_max + 1;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        //for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
        for(int wg=wg_tmp; wg<n_loop_max; wg++){
            n_walk_ar[wg]        = n_ipg / n_loop_max;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_ar[i] = 0;
        }
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
        bool first_loop = true;
        S32 n_walk_prev = 0;
        if(n_ipg > 0){
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = n_walk_ar[wg];
                const S32 walk_grp_head = n_disp_walk_ar[wg];
                const F64 offset_calc_force__core__walk_tree = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
                {
                    const S32 ith = Comm::getThreadNum();
                    n_ep_cum_thread[ith] = cnt_thread[ith] = 0;
                    n_epj_disp_thread[ith][0] = 0;
                    adr_epj_tmp[ith].clearSize();
                    adr_ipg_tmp[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4) 
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        adr_ipg_tmp[ith].push_back(id_ipg);
                        const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_; 
                        const S32 ith = Comm::getThreadNum();
                        n_epi_ar[iw] = ipg_[id_ipg].n_ptcl_;
                        epi_ar[iw]   = epi_sorted_.getPointer(first_adr_ip);
                        //force_array[iw] = force_sorted_.getPointer(first_adr_ip);
                        force_array[iw] = force_.getPointer(first_adr_ip);

                        TargetBox<TSM> target_box;
                        GetTargetBox<TSM>(ipg_[id_ipg], target_box);
                        S32 adr_tc = 0;
                        MakeListUsingTreeRecursiveTop
                            <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj,
                             TagWalkModeNormal, TagChopLeafTrue>
                            (tc_glb_,  adr_tc, tp_glb_,
                             epj_sorted_, adr_epj_tmp[ith],
                             target_box,
                             r_crit_sq, n_leaf_limit_,
                             F64vec(0.0));

                        /*
                        F64 mass_tmp = 0.0;
                        for(S32 i=0; i<adr_epj_tmp[ith].size(); i++){
                            mass_tmp += epj_sorted_[ adr_epj_tmp[ith][i] ].mass;
                        }
                        assert(fmod(mass_tmp, 1.0)==0.0);
                        */

                        n_epj_ar[iw] = adr_epj_tmp[ith].size() - n_ep_cum_thread[ith];

                        n_ep_cum_thread[ith] = adr_epj_tmp[ith].size();

                        n_epj_disp_thread[ith][cnt_thread[ith]+1] = n_ep_cum_thread[ith];

                        n_interaction_ep_ep_ar[ith] += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[iw]);

                        iw2ith[iw] = ith;
                        iw2cnt[iw] = cnt_thread[ith];
                        cnt_thread[ith]++;
                    } // end of OMP for
                } // end of OMP parallel scope

                if(flag_keep_list){
                    interaction_list_.n_disp_ep_[0] = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_ep_[id_ipg] = n_epj_ar[iw];
                    }
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_disp_ep_[id_ipg+1] = interaction_list_.n_disp_ep_[id_ipg] + interaction_list_.n_ep_[id_ipg];
                    }
                    interaction_list_.adr_ep_.resizeNoInitialize( interaction_list_.n_disp_ep_[walk_grp_head+n_walk] );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(S32 i=0; i<n_thread; i++){
                        for(S32 j=0; j<adr_ipg_tmp[i].size(); j++){
                            const S32 adr_ipg = adr_ipg_tmp[i][j];
                            S32 adr_ep = interaction_list_.n_disp_ep_[adr_ipg];
                            const S32 k_ep_h = n_epj_disp_thread[i][j];
                            const S32 k_ep_e = n_epj_disp_thread[i][j+1];
                            for(S32 k=k_ep_h; k<k_ep_e; k++, adr_ep++){
                                interaction_list_.adr_ep_[adr_ep] = adr_epj_tmp[i][k];
                            }
                        }
                    }
                }
                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;
#if 1
                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    id_epj_ar[iw] = adr_epj_tmp[ith].getPointer(n_ep_head);
                }
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar_prev.getPointer(), force_prev_array);
                } // retrieve
#else
                // original
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar_prev.getPointer(), force_prev_array);
                } // retrieve
                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    id_epj_ar[iw] = adr_epj_tmp[ith].getPointer(n_ep_head);
                }
#endif
                pfunc_dispatch(0, n_walk, 
                               (const Tepi**)epi_ar.getPointer(),   n_epi_ar.getPointer(),
                               (const S32**)id_epj_ar.getPointer(), n_epj_ar.getPointer(),
                               epj_sorted_.getPointer(), epj_sorted_.size(),
                               false);

                first_loop = false;
                for(int iw=0; iw<n_walk; iw++){
                    n_epi_ar_prev[iw] = n_epi_ar[iw];
                    force_prev_array[iw] = force_array[iw];
                }
                n_walk_prev = n_walk;
            } // end of walk group loop
            ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar_prev.getPointer(), force_prev_array);
        } // if(n_ipg > 0)
        else{
            ni_ave_ = nj_ave_ = 0;
            n_interaction_ep_ep_local_ = 0;
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_ar[i];
        }
        time_profile_.calc_force__core += GetWtime() - offset_core;
        return ret;
    }

    //////////////////////////////////////////////////////////////
    //////////// Walk+Force, Kernel:Ptcl, List:Ptcl //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve, class Tsys>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalk2(Tfunc_dispatch pfunc_dispatch,
                        Tfunc_retrieve pfunc_retrieve,
                        Tsys & sys,
                        const S32 tag_max,
                        const S32 n_walk_limit,
                        const bool clear){
        if(tag_max <= 0){
            PARTICLE_SIMULATOR_PRINT_ERROR("tag_max is illegal. In currente version, tag_max must be 1");
            Abort(-1);
        }
        force_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                force_[i].clear();
            }
        }
        S32 ret = 0;
        const F64 wtime_offset = GetWtime();
        ret = calcForceMultiWalk2Impl(typename TSM::force_type(),
                                     pfunc_dispatch,
                                     pfunc_retrieve,
                                      sys,
                                     tag_max,
                                     n_walk_limit,
                                     clear);
        time_profile_.calc_force += GetWtime() - wtime_offset;
        return ret;
    }

    //////////// Walk+Force, Kernel:Ptcl, List:Ptcl, Force:Long//////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve, class Tsys>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalk2Impl(TagForceLong,
                            Tfunc_dispatch pfunc_dispatch,
                            Tfunc_retrieve pfunc_retrieve,
                            Tsys & sys,
                            const S32 tag_max,
                            const S32 n_walk_limit,
                            const bool clear){
        const F64 offset_core = GetWtime();
        static bool first = true;
        S32 ret = 0;
        S32 tag = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        //const S32 n_ipg_amari = (n_ipg > 0) ? n_ipg%n_walk_limit : 0;
        if(n_ipg <= 0) return 0;
        n_walk_local_ += n_ipg;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        //const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg_amari==0 ? 0 : 1);
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
                n_epj_disp_thread[i] = new S32[n_walk_limit+1];
                n_spj_disp_thread[i] = new S32[n_walk_limit+1];
            }
            cnt_thread = new S32[n_thread];
            n_ep_cum_thread = new S32[n_thread];
            n_sp_cum_thread = new S32[n_thread];
            n_interaction_ep_ep_array = new S64[n_thread];
            n_interaction_ep_sp_array = new S64[n_thread];
            first = false;
        }
        walk_grp_disp[0] = 0;
        //const S32 wg_tmp = n_ipg > 0 ? n_ipg%n_loop_max : 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
        //for(int wg=0; wg<wg_tmp; wg++){
            walk_grp[wg] = n_ipg / n_loop_max + 1;
            walk_grp_disp[wg+1] = walk_grp_disp[wg] + walk_grp[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
        //for(int wg=wg_tmp; wg<n_loop_max; wg++){
            walk_grp[wg] = n_ipg / n_loop_max;
            walk_grp_disp[wg+1] = walk_grp_disp[wg] + walk_grp[wg];
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_array[i] = n_interaction_ep_sp_array[i] = 0;
        }
        bool first_loop = true;
        S32 n_walk_prev = 0;
        CountT n_epj_for_force_tmp = 0;
        CountT n_spj_for_force_tmp = 0;
        if(n_ipg > 0){
            S32 ip_head = 0, ip_tail=0, ip_head_prev=0, ip_tail_prev=0;
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = walk_grp[wg];
                const S32 walk_grp_head = walk_grp_disp[wg];
                ip_head = ipg_[walk_grp_head].adr_ptcl_; 
                ip_tail = ipg_[walk_grp_head+n_walk-1].adr_ptcl_ + ipg_[walk_grp_head+n_walk-1].n_ptcl_;
                for(int i=0; i<n_thread; i++){
                    n_ep_cum_thread[i] = n_sp_cum_thread[i] = cnt_thread[i] = 0;
                    n_epj_disp_thread[i][0] = 0;
                    n_spj_disp_thread[i][0] = 0;
                    epj_for_force_[i].clearSize();
                    spj_for_force_[i].clearSize();
                }
                const F64 offset_calc_force__core__walk_tree = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : n_epj_for_force_tmp, n_spj_for_force_tmp)
#endif
                for(int iw=0; iw<n_walk; iw++){
                    const S32 id_ipg = walk_grp_head + iw;
                    const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_; 
                    const S32 ith = Comm::getThreadNum();
                    n_epi_array[iw] = ipg_[id_ipg].n_ptcl_;
                    epi_array[iw]   = epi_sorted_.getPointer(first_adr_ip);
                    //force_array[iw] = force_sorted_.getPointer(first_adr_ip);
                    force_array[iw] = force_.getPointer(first_adr_ip);
                    makeInteractionList(id_ipg, false);
                    n_epj_array[iw] = epj_for_force_[ith].size() - n_ep_cum_thread[ith];
                    n_spj_array[iw] = spj_for_force_[ith].size() - n_sp_cum_thread[ith];
                    n_ep_cum_thread[ith] = epj_for_force_[ith].size();
                    n_sp_cum_thread[ith] = spj_for_force_[ith].size();
                    n_epj_disp_thread[ith][cnt_thread[ith]+1] = n_ep_cum_thread[ith];
                    n_spj_disp_thread[ith][cnt_thread[ith]+1] = n_sp_cum_thread[ith];
                    n_interaction_ep_ep_array[ith] += ((S64)n_epj_array[iw]*(S64)n_epi_array[iw]);
                    n_interaction_ep_sp_array[ith] += ((S64)n_spj_array[iw]*(S64)n_epi_array[iw]);
                    n_epj_for_force_tmp += (CountT)n_epj_array[iw];
                    n_spj_for_force_tmp += (CountT)n_spj_array[iw];
                    iw2ith[iw] = ith;
                    iw2cnt[iw] = cnt_thread[ith];
                    cnt_thread[ith]++;
                }
                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;

                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    S32 n_sp_head = n_spj_disp_thread[ith][cnt];
                    epj_array[iw] = epj_for_force_[ith].getPointer(n_ep_head);
                    spj_array[iw] = spj_for_force_[ith].getPointer(n_sp_head);
                }
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_array_prev, force_prev_array);
                } // retrieve
                ret += pfunc_dispatch(tag, n_walk, (const Tepi**)epi_array, n_epi_array, (const Tepj**)epj_array, n_epj_array, (const Tspj**)spj_array, n_spj_array);


                if(!first_loop){
                    F64 wtime_offset_wb_int_cp = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for
#endif
                    for(S32 ip=ip_head_prev; ip<ip_tail_prev; ip++){
                        sys[ip].copyFromForce(force_[ip]);
                        sys[ip].integrate();
                        //epi_org_[ip].copyFromFP( sys[ip] );
                        //epj_org_[ip].copyFromFP( sys[ip] );
                    }
                    time_profile_.calc_force__wb_int_cp += GetWtime() - wtime_offset_wb_int_cp;
                }

                for(int iw=0; iw<n_walk; iw++){
                    n_epi_array_prev[iw] = n_epi_array[iw];
                    force_prev_array[iw] = force_array[iw];
                }
                n_walk_prev = n_walk;
                ip_head_prev = ip_head;
                ip_tail_prev = ip_tail;
                first_loop = false;
            } // end of walk group loop
            ret += pfunc_retrieve(tag, n_walk_prev, n_epi_array_prev, force_prev_array);

            F64 wtime_offset_wb_int_cp = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for
#endif
            for(S32 ip=ip_head_prev; ip<ip_tail_prev; ip++){
                sys[ip].copyFromForce(force_[ip]);
                sys[ip].integrate();
                //epi_org_[ip].copyFromFP( sys[ip] );
                //epj_org_[ip].copyFromFP( sys[ip] );
            }
            time_profile_.calc_force__wb_int_cp += GetWtime() - wtime_offset_wb_int_cp;
        } // if(n_ipg > 0)
        else{
            ni_ave_ = nj_ave_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
            n_epj_for_force_ = n_spj_for_force_ = 0;
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_array[i];
            n_interaction_ep_sp_local_ += n_interaction_ep_sp_array[i];
        }
        n_epj_for_force_ += n_epj_for_force_tmp;
        n_spj_for_force_ += n_spj_for_force_tmp;
        time_profile_.calc_force__core += GetWtime() - offset_core;
        return ret;
    }


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
        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        force_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                //force_sorted_[i].clear();
                //force_org_[i].clear();
                force_[i].clear();
            }
        }
        S32 ret = 0;
        const F64 wtime_offset = GetWtime();
        ret = calcForceMultiWalkImpl(typename TSM::force_type(),
                                     pfunc_dispatch,
                                     pfunc_retrieve,
                                     tag_max,
                                     n_walk_limit,
                                     clear);
        time_profile_.calc_force += GetWtime() - wtime_offset;
        return ret;
    }

    //////////// Walk+Force, Kernel:Ptcl, List:Ptcl, Force:Long//////////////
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
        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        //const S32 n_ipg_amari = (n_ipg > 0) ? n_ipg%n_walk_limit : 0;
        if(n_ipg <= 0) return 0;
        n_walk_local_ += n_ipg;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        //const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg_amari==0 ? 0 : 1);
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
                n_epj_disp_thread[i] = new S32[n_walk_limit+1];
                n_spj_disp_thread[i] = new S32[n_walk_limit+1];
            }
            cnt_thread = new S32[n_thread];
            n_ep_cum_thread = new S32[n_thread];
            n_sp_cum_thread = new S32[n_thread];
            n_interaction_ep_ep_array = new S64[n_thread];
            n_interaction_ep_sp_array = new S64[n_thread];
            first = false;
        }
        walk_grp_disp[0] = 0;
        //const S32 wg_tmp = n_ipg > 0 ? n_ipg%n_loop_max : 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
        //for(int wg=0; wg<wg_tmp; wg++){
            walk_grp[wg] = n_ipg / n_loop_max + 1;
            walk_grp_disp[wg+1] = walk_grp_disp[wg] + walk_grp[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
        //for(int wg=wg_tmp; wg<n_loop_max; wg++){
            walk_grp[wg] = n_ipg / n_loop_max;
            walk_grp_disp[wg+1] = walk_grp_disp[wg] + walk_grp[wg];
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_array[i] = n_interaction_ep_sp_array[i] = 0;
        }
        bool first_loop = true;
        S32 n_walk_prev = 0;
        CountT n_epj_for_force_tmp = 0;
        CountT n_spj_for_force_tmp = 0;
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
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : n_epj_for_force_tmp, n_spj_for_force_tmp)
#endif
                for(int iw=0; iw<n_walk; iw++){
                    const S32 id_ipg = walk_grp_head + iw;
                    const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_; 
                    const S32 ith = Comm::getThreadNum();
                    n_epi_array[iw] = ipg_[id_ipg].n_ptcl_;
                    epi_array[iw]   = epi_sorted_.getPointer(first_adr_ip);
                    //force_array[iw] = force_sorted_.getPointer(first_adr_ip);
                    force_array[iw] = force_.getPointer(first_adr_ip);
                    makeInteractionList(id_ipg, false);
                    n_epj_array[iw] = epj_for_force_[ith].size() - n_ep_cum_thread[ith];
                    n_spj_array[iw] = spj_for_force_[ith].size() - n_sp_cum_thread[ith];
                    n_ep_cum_thread[ith] = epj_for_force_[ith].size();
                    n_sp_cum_thread[ith] = spj_for_force_[ith].size();
                    n_epj_disp_thread[ith][cnt_thread[ith]+1] = n_ep_cum_thread[ith];
                    n_spj_disp_thread[ith][cnt_thread[ith]+1] = n_sp_cum_thread[ith];
                    n_interaction_ep_ep_array[ith] += ((S64)n_epj_array[iw]*(S64)n_epi_array[iw]);
                    n_interaction_ep_sp_array[ith] += ((S64)n_spj_array[iw]*(S64)n_epi_array[iw]);
                    n_epj_for_force_tmp += (CountT)n_epj_array[iw];
                    n_spj_for_force_tmp += (CountT)n_spj_array[iw];
                    iw2ith[iw] = ith;
                    iw2cnt[iw] = cnt_thread[ith];
                    cnt_thread[ith]++;
                }
                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;
#if 1
                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    S32 n_sp_head = n_spj_disp_thread[ith][cnt];
                    epj_array[iw] = epj_for_force_[ith].getPointer(n_ep_head);
                    spj_array[iw] = spj_for_force_[ith].getPointer(n_sp_head);
                }
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_array_prev, force_prev_array);
                } // retrieve
#else
                //original
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
#endif
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
            ni_ave_ = nj_ave_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
            n_epj_for_force_ = n_spj_for_force_ = 0;
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_array[i];
            n_interaction_ep_sp_local_ += n_interaction_ep_sp_array[i];
        }
        n_epj_for_force_ += n_epj_for_force_tmp;
        n_spj_for_force_ += n_spj_for_force_tmp;
        time_profile_.calc_force__core += GetWtime() - offset_core;
        return ret;
    }

    //////////// Walk+Force, Kernel:Ptcl, List:Ptcl, Force:Short//////////////
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
        //std::cerr<<"rank(A)= "<<Comm::getRank()<<std::endl;
        const F64 offset_core = GetWtime();
        static bool first = true;
        S32 ret = 0;
        S32 tag = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        if(n_ipg <= 0) return 0;
        //const S32 n_ipg_amari = (n_ipg > 0) ? n_ipg%n_walk_limit : 0;
        n_walk_local_ += n_ipg;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        //const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg_amari==0 ? 0 : 1);
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
                n_epj_disp_thread[i] = new S32[n_walk_limit+1];
            }
            cnt_thread = new S32[n_thread];
            n_ep_cum_thread = new S32[n_thread];
            n_interaction_ep_ep_array = new S64[n_thread];
            first = false;
        }
        walk_grp_disp[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            walk_grp[wg] = (n_ipg/n_loop_max) + 1;
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
                    //force_array[iw] = force_sorted_.getPointer(first_adr_ip);
                    force_array[iw] = force_.getPointer(first_adr_ip);
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
#if 1
                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    epj_array[iw] = epj_for_force_[ith].getPointer(n_ep_head);
                }
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_array_prev, force_prev_array);
                } // retrieve
#else
                // original
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_array_prev, force_prev_array);
                } // retrieve
                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    epj_array[iw] = epj_for_force_[ith].getPointer(n_ep_head);
                }
#endif
                ret += pfunc_dispatch(tag, n_walk, (const Tepi**)epi_array, n_epi_array, (const Tepj**)epj_array, n_epj_array);
                first_loop = false;
                for(int iw=0; iw<n_walk; iw++){
                    n_epi_array_prev[iw] = n_epi_array[iw];
                    force_prev_array[iw] = force_array[iw];
                }
                n_walk_prev = n_walk;
                //std::cerr<<"rank(E)= "<<Comm::getRank()<<std::endl;
            } // end of walk group loop
            ret += pfunc_retrieve(tag, n_walk_prev, n_epi_array_prev, force_prev_array);
        } // if(n_ipg > 0)
        else{
            ni_ave_ = nj_ave_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
        }

        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_array[i];
        }
        time_profile_.calc_force__core += GetWtime() - offset_core;
        return ret;
    }

    //////////////////////////////////////////////////
    //////////// Walk+Force, Kernel:Ptcl, List:Index //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkIndexNew(Tfunc_dispatch pfunc_dispatch,
                               Tfunc_retrieve pfunc_retrieve,
                               const S32 tag_max,
                               const S32 n_walk_limit,
                               const bool flag_keep_list,
                               const bool clear){
        if(tag_max <= 0){
            PARTICLE_SIMULATOR_PRINT_ERROR("tag_max is illegal. In currente version, tag_max must be 1");
            Abort(-1);
        }
        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        force_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                //force_sorted_[i].clear();
                //force_org_[i].clear();
                force_[i].clear();
            }
        }
        S32 ret = 0;
        const F64 time_offset = GetWtime();
        ret = calcForceMultiWalkIndexNewImpl(typename TSM::force_type(),
                                             pfunc_dispatch,
                                             pfunc_retrieve,
                                             tag_max,
                                             n_walk_limit,
                                             flag_keep_list,
                                             clear);
        time_profile_.calc_force += GetWtime() - time_offset;
        return ret;
    }

    //////////// Walk+Force, Kernel:Ptcl, List:Index Force:Long //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkIndexNewImpl(TagForceLong,
                                   Tfunc_dispatch pfunc_dispatch,
                                   Tfunc_retrieve pfunc_retrieve,
                                   const S32 tag_max,
                                   const S32 n_walk_limit,
                                   const bool flag_keep_list,
                                   const bool clear){
        const F64 offset_core = GetWtime();
        static bool first = true;
        S32 ret = 0;
        S32 tag = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        if(n_ipg <= 0) return 0;
        n_walk_local_ += n_ipg;

        if(flag_keep_list){
            interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
            interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
            interaction_list_.adr_ep_.clearSize();
            interaction_list_.n_sp_.resizeNoInitialize(n_ipg);
            interaction_list_.n_disp_sp_.resizeNoInitialize(n_ipg+1);
            interaction_list_.adr_sp_.clearSize();
        }

        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        static std::vector<S32> n_walk_ar;
        n_walk_ar.resize(n_loop_max); // group of walk (walk_grp[i] < n_walk_limit)
        static std::vector<S32> n_disp_walk_ar;
        n_disp_walk_ar.resize(n_loop_max+1);

        static S32  * iw2ith;
        static S32  * iw2cnt;
        static S32 ** n_epj_disp_thread; // [n_thread][n_walk]
        static S32 ** n_spj_disp_thread;// [n_thread][n_walk]
        static Tforce ** force_array; // array of pointer *[n_walk]
        static Tforce ** force_prev_array; // array of pointer *[n_walk]
        static S32  * cnt_thread;
        static S32  * n_ep_cum_thread;
        static S32  * n_sp_cum_thread;
        static S64 * n_interaction_ep_ep_ar;
        static S64 * n_interaction_ep_sp_ar;

        static ReallocatableArray<S32> * adr_epj_tmp;
        static ReallocatableArray<S32> * adr_spj_tmp;
        static ReallocatableArray<S32> * adr_ipg_tmp;

        static ReallocatableArray<S32> n_epi_ar;
        n_epi_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32> n_epi_ar_prev;
        n_epi_ar_prev.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tepi*> epi_ar; // array of pointer *[n_walk]
        epi_ar.resizeNoInitialize(n_walk_limit);

        static ReallocatableArray<S32> n_epj_ar;
        n_epj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tepj*> epj_ar;
        epj_ar.resizeNoInitialize(n_walk_limit);

        static ReallocatableArray<S32> n_spj_ar;
        n_spj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tspj*> spj_ar;
        spj_ar.resizeNoInitialize(n_walk_limit);

        if(first){
            iw2ith = new S32[n_walk_limit];
            iw2cnt = new S32[n_walk_limit];
            n_epj_disp_thread = new S32*[n_thread];
            n_spj_disp_thread = new S32*[n_thread];
            force_array      = new Tforce*[n_walk_limit];
            force_prev_array = new Tforce*[n_walk_limit];
            for(int i=0; i<n_thread; i++){
                n_epj_disp_thread[i] = new S32[n_walk_limit+1];
                n_spj_disp_thread[i] = new S32[n_walk_limit+1];
            }
            cnt_thread = new S32[n_thread];
            n_ep_cum_thread = new S32[n_thread];
            n_sp_cum_thread = new S32[n_thread];
            n_interaction_ep_ep_ar = new S64[n_thread];
            n_interaction_ep_sp_ar = new S64[n_thread];

            adr_epj_tmp = new ReallocatableArray<S32>[n_thread];
            adr_spj_tmp = new ReallocatableArray<S32>[n_thread];
            adr_ipg_tmp = new ReallocatableArray<S32>[n_thread];

            first = false;
        }
        n_disp_walk_ar[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            n_walk_ar[wg]        = n_ipg / n_loop_max + 1;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
            n_walk_ar[wg]        = n_ipg / n_loop_max;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_ar[i] = n_interaction_ep_sp_ar[i] = 0;
        }
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
        const S32 adr_tree_sp_first = spj_org_.size();
        bool first_loop = true;
        S32 n_walk_prev = 0;
        CountT n_epj_for_force_tmp = 0;
        CountT n_spj_for_force_tmp = 0;
        if(n_ipg > 0){
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = n_walk_ar[wg];
                const S32 walk_grp_head = n_disp_walk_ar[wg];
                const F64 offset_calc_force__core__walk_tree = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel reduction(+ : n_epj_for_force_tmp, n_spj_for_force_tmp)
#endif
                {
                    const S32 ith = Comm::getThreadNum();
                    n_ep_cum_thread[ith] = n_sp_cum_thread[ith] = cnt_thread[ith] = 0;
                    n_epj_disp_thread[ith][0] = 0;
                    n_spj_disp_thread[ith][0] = 0;
                    adr_epj_tmp[ith].clearSize();
                    adr_spj_tmp[ith].clearSize();
                    adr_ipg_tmp[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4) 
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        adr_ipg_tmp[ith].push_back(id_ipg);
                        const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_; 
                        const S32 ith = Comm::getThreadNum();
                        n_epi_ar[iw] = ipg_[id_ipg].n_ptcl_;
                        epi_ar[iw]   = epi_sorted_.getPointer(first_adr_ip);
                        //force_array[iw] = force_sorted_.getPointer(first_adr_ip);
                        force_array[iw] = force_.getPointer(first_adr_ip);

                        TargetBox<TSM> target_box;
                        GetTargetBox<TSM>(ipg_[id_ipg], target_box);
                        S32 adr_tc = 0;
                        MakeListUsingTreeRecursiveTop
                            <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, 
                             Tspj, TagWalkModeNormal, TagChopLeafTrue>
                            (tc_glb_,  adr_tc, tp_glb_,
                             epj_sorted_, adr_epj_tmp[ith],
                             spj_sorted_, adr_spj_tmp[ith],
                             target_box,
                             r_crit_sq, n_leaf_limit_,
                             adr_tree_sp_first, F64vec(0.0));

                        /*
                        F64 mass_tmp = 0.0;
                        for(S32 i=0; i<adr_epj_tmp[ith].size(); i++){
                            mass_tmp += epj_sorted_[ adr_epj_tmp[ith][i] ].mass;
                        }
                        for(S32 i=0; i<adr_spj_tmp[ith].size(); i++){
                            mass_tmp += spj_sorted_[ adr_spj_tmp[ith][i] ].mass;
                        }
                        assert(fmod(mass_tmp, 1.0)==0.0);
                        */
                        n_epj_ar[iw] = adr_epj_tmp[ith].size() - n_ep_cum_thread[ith];
                        n_spj_ar[iw] = adr_spj_tmp[ith].size() - n_sp_cum_thread[ith];

                        n_epj_for_force_tmp += (CountT)n_epj_ar[iw];
                        n_spj_for_force_tmp += (CountT)n_spj_ar[iw];

                        n_ep_cum_thread[ith] = adr_epj_tmp[ith].size();
                        n_sp_cum_thread[ith] = adr_spj_tmp[ith].size();

                        n_epj_disp_thread[ith][cnt_thread[ith]+1] = n_ep_cum_thread[ith];
                        n_spj_disp_thread[ith][cnt_thread[ith]+1] = n_sp_cum_thread[ith];

                        n_interaction_ep_ep_ar[ith] += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[iw]);
                        n_interaction_ep_sp_ar[ith] += ((S64)n_spj_ar[iw]*(S64)n_epi_ar[iw]);

                        iw2ith[iw] = ith;
                        iw2cnt[iw] = cnt_thread[ith];
                        cnt_thread[ith]++;

                    } // end of OMP for
                } // end of OMP parallel scope

                if(flag_keep_list){
                    interaction_list_.n_disp_ep_[0] = interaction_list_.n_disp_sp_[0] = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_ep_[id_ipg] = n_epj_ar[iw];
                        interaction_list_.n_sp_[id_ipg] = n_spj_ar[iw];
                    }
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_disp_ep_[id_ipg+1] = interaction_list_.n_disp_ep_[id_ipg] + interaction_list_.n_ep_[id_ipg];
                        interaction_list_.n_disp_sp_[id_ipg+1] = interaction_list_.n_disp_sp_[id_ipg] + interaction_list_.n_sp_[id_ipg];
                    }

                    interaction_list_.adr_ep_.resizeNoInitialize( interaction_list_.n_disp_ep_[walk_grp_head+n_walk] );
                    interaction_list_.adr_sp_.resizeNoInitialize( interaction_list_.n_disp_sp_[walk_grp_head+n_walk] );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(S32 i=0; i<n_thread; i++){
                        for(S32 j=0; j<adr_ipg_tmp[i].size(); j++){
                            const S32 adr_ipg = adr_ipg_tmp[i][j];
                            S32 adr_ep = interaction_list_.n_disp_ep_[adr_ipg];
                            const S32 k_ep_h = n_epj_disp_thread[i][j];
                            const S32 k_ep_e = n_epj_disp_thread[i][j+1];
                            for(S32 k=k_ep_h; k<k_ep_e; k++, adr_ep++){
                                interaction_list_.adr_ep_[adr_ep] = adr_epj_tmp[i][k];
                            }
                            S32 adr_sp = interaction_list_.n_disp_sp_[adr_ipg];
                            const S32 k_sp_h = n_spj_disp_thread[i][j];
                            const S32 k_sp_e = n_spj_disp_thread[i][j+1];
                            for(S32 k=k_sp_h; k<k_sp_e; k++, adr_sp++){
                                interaction_list_.adr_sp_[adr_sp] = adr_spj_tmp[i][k];
                            }
                        }
                    }
                }

                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar_prev.getPointer(), force_prev_array);
                } // retrieve
#if 1
                epj_for_force_[0].clearSize();
                spj_for_force_[0].clearSize();
                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    S32 n_sp_head = n_spj_disp_thread[ith][cnt];
                    S32 * id_epj = adr_epj_tmp[ith].getPointer(n_ep_head);
                    for(S32 jp=0; jp<n_epj_ar[iw]; jp++){
                        epj_for_force_[0].push_back( epj_sorted_[id_epj[jp]] );
                    }
                    S32 * id_spj = adr_spj_tmp[ith].getPointer(n_sp_head);
                    for(S32 jp=0; jp<n_spj_ar[iw]; jp++){
                        spj_for_force_[0].push_back( spj_sorted_[id_spj[jp]] );
                    }
                }
                S64 n_epj_cnt = 0;
                S64 n_spj_cnt = 0;
                for(S32 iw=0; iw<n_walk; iw++){
                    epj_ar[iw] = epj_for_force_[0].getPointer(n_epj_cnt);
                    spj_ar[iw] = spj_for_force_[0].getPointer(n_spj_cnt);
                    n_epj_cnt += n_epj_ar[iw];
                    n_spj_cnt += n_spj_ar[iw];
                }
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar_prev.getPointer(), force_prev_array);
                } // retrieve
#else
                // original
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar_prev.getPointer(), force_prev_array);
                } // retrieve
                epj_for_force_[0].clearSize();
                spj_for_force_[0].clearSize();
                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    S32 n_sp_head = n_spj_disp_thread[ith][cnt];
                    S32 * id_epj = adr_epj_tmp[ith].getPointer(n_ep_head);
                    for(S32 jp=0; jp<n_epj_ar[iw]; jp++){
                        epj_for_force_[0].push_back( epj_sorted_[id_epj[jp]] );
                    }
                    S32 * id_spj = adr_spj_tmp[ith].getPointer(n_sp_head);
                    for(S32 jp=0; jp<n_spj_ar[iw]; jp++){
                        spj_for_force_[0].push_back( spj_sorted_[id_spj[jp]] );
                    }
                }
                S64 n_epj_cnt = 0;
                S64 n_spj_cnt = 0;
                for(S32 iw=0; iw<n_walk; iw++){
                    epj_ar[iw] = epj_for_force_[0].getPointer(n_epj_cnt);
                    spj_ar[iw] = spj_for_force_[0].getPointer(n_spj_cnt);
                    n_epj_cnt += n_epj_ar[iw];
                    n_spj_cnt += n_spj_ar[iw];
                }
#endif
                ret += pfunc_dispatch(tag, n_walk, 
                                      (const Tepi**)epi_ar.getPointer(), n_epi_ar.getPointer(), 
                                      (const Tepj**)epj_ar.getPointer(), n_epj_ar.getPointer(), 
                                      (const Tspj**)spj_ar.getPointer(), n_spj_ar.getPointer());
                first_loop = false;
                for(int iw=0; iw<n_walk; iw++){
                    n_epi_ar_prev[iw] = n_epi_ar[iw];
                    force_prev_array[iw] = force_array[iw];
                }
                n_walk_prev = n_walk;
            } // end of walk group loop
            ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar_prev.getPointer(), force_prev_array);
        } // if(n_ipg > 0)
        else{
            ni_ave_ = nj_ave_ = 0;
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = 0;
            n_epj_for_force_ = n_spj_for_force_ = 0;
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_ar[i];
            n_interaction_ep_sp_local_ += n_interaction_ep_sp_ar[i];
        }
        n_epj_for_force_ += n_epj_for_force_tmp;
        n_spj_for_force_ += n_spj_for_force_tmp;
        /*
        std::cerr<<"n_epj_for_force_= "<<n_epj_for_force_
                 <<" n_spj_for_force_= "<<n_spj_for_force_
                 <<std::endl;
        */
        time_profile_.calc_force__core += GetWtime() - offset_core;
        return ret;
    }
    

    //////////// Walk+Force, Kernel:Ptcl, List:Index Force:Short //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceMultiWalkIndexNewImpl(TagForceShort,
                                   Tfunc_dispatch pfunc_dispatch,
                                   Tfunc_retrieve pfunc_retrieve,
                                   const S32 tag_max,
                                   const S32 n_walk_limit,
                                   const bool flag_keep_list,
                                   const bool clear){
        const F64 offset_core = GetWtime();
        static bool first = true;
        S32 ret = 0;
        S32 tag = 0;
        const S32 n_thread = Comm::getNumberOfThread();
        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        const S32 n_ipg = ipg_.size();
        if(n_ipg <= 0) return 0;
        n_walk_local_ += n_ipg;

        if(flag_keep_list){
            interaction_list_.n_ep_.resizeNoInitialize(n_ipg);
            interaction_list_.n_disp_ep_.resizeNoInitialize(n_ipg+1);
            interaction_list_.adr_ep_.clearSize();
        }

        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        static std::vector<S32> n_walk_ar;
        n_walk_ar.resize(n_loop_max); // group of walk (walk_grp[i] < n_walk_limit)
        static std::vector<S32> n_disp_walk_ar;
        n_disp_walk_ar.resize(n_loop_max+1);

        static S32  * iw2ith;
        static S32  * iw2cnt;
        static S32 ** n_epj_disp_thread; // [n_thread][n_walk]
        static Tforce ** force_array; // array of pointer *[n_walk]
        static Tforce ** force_prev_array; // array of pointer *[n_walk]
        static S32  * cnt_thread;
        static S32  * n_ep_cum_thread;
        static S64 * n_interaction_ep_ep_ar;

        static ReallocatableArray<S32> * adr_epj_tmp;
        static ReallocatableArray<S32> * adr_ipg_tmp;

        static ReallocatableArray<S32> n_epi_ar;
        n_epi_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32> n_epi_ar_prev;
        n_epi_ar_prev.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tepi*> epi_ar; // array of pointer *[n_walk]
        epi_ar.resizeNoInitialize(n_walk_limit);

        static ReallocatableArray<S32> n_epj_ar;
        n_epj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tepj*> epj_ar;
        epj_ar.resizeNoInitialize(n_walk_limit);

        if(first){
            iw2ith = new S32[n_walk_limit];
            iw2cnt = new S32[n_walk_limit];
            n_epj_disp_thread = new S32*[n_thread];
            force_array      = new Tforce*[n_walk_limit];
            force_prev_array = new Tforce*[n_walk_limit];
            for(int i=0; i<n_thread; i++){
                n_epj_disp_thread[i] = new S32[n_walk_limit+1];
            }
            cnt_thread = new S32[n_thread];
            n_ep_cum_thread = new S32[n_thread];
            n_interaction_ep_ep_ar = new S64[n_thread];

            adr_epj_tmp = new ReallocatableArray<S32>[n_thread];
            adr_ipg_tmp = new ReallocatableArray<S32>[n_thread];

            first = false;
        }
        n_disp_walk_ar[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            n_walk_ar[wg]        = n_ipg / n_loop_max + 1;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
            n_walk_ar[wg]        = n_ipg / n_loop_max;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_ar[i] = 0;
        }
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
        bool first_loop = true;
        S32 n_walk_prev = 0;

        if(n_ipg > 0){
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = n_walk_ar[wg];
                const S32 walk_grp_head = n_disp_walk_ar[wg];
                const F64 offset_calc_force__core__walk_tree = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
                {
                    const S32 ith = Comm::getThreadNum();
                    n_ep_cum_thread[ith] = cnt_thread[ith] = 0;
                    n_epj_disp_thread[ith][0] = 0;
                    adr_epj_tmp[ith].clearSize();
                    adr_ipg_tmp[ith].clearSize();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4) 
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        adr_ipg_tmp[ith].push_back(id_ipg);
                        const S32 first_adr_ip = ipg_[id_ipg].adr_ptcl_; 
                        const S32 ith = Comm::getThreadNum();
                        n_epi_ar[iw] = ipg_[id_ipg].n_ptcl_;
                        epi_ar[iw]   = epi_sorted_.getPointer(first_adr_ip);
                        //force_array[iw] = force_sorted_.getPointer(first_adr_ip);
                        force_array[iw] = force_.getPointer(first_adr_ip);

                        TargetBox<TSM> target_box;
                        GetTargetBox<TSM>(ipg_[id_ipg], target_box);
                        S32 adr_tc = 0;
                        MakeListUsingTreeRecursiveTop
                            <TSM, TreeCell<Tmomglb>, TreeParticle, Tepj,
                             TagWalkModeNormal, TagChopLeafTrue>
                            (tc_glb_,  adr_tc, tp_glb_,
                             epj_sorted_, adr_epj_tmp[ith],
                             target_box,
                             r_crit_sq, n_leaf_limit_,
                             F64vec(0.0));

                        /*
                        F64 mass_tmp = 0.0;
                        for(S32 i=0; i<adr_epj_tmp[ith].size(); i++){
                            mass_tmp += epj_sorted_[ adr_epj_tmp[ith][i] ].mass;
                        }
                        assert(fmod(mass_tmp, 1.0)==0.0);
                        */

                        n_epj_ar[iw] = adr_epj_tmp[ith].size() - n_ep_cum_thread[ith];

                        n_ep_cum_thread[ith] = adr_epj_tmp[ith].size();

                        n_epj_disp_thread[ith][cnt_thread[ith]+1] = n_ep_cum_thread[ith];

                        n_interaction_ep_ep_ar[ith] += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[iw]);

                        iw2ith[iw] = ith;
                        iw2cnt[iw] = cnt_thread[ith];
                        cnt_thread[ith]++;
                    } // end of OMP for
                } // end of OMP parallel scope

                if(flag_keep_list){
                    interaction_list_.n_disp_ep_[0] = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_ep_[id_ipg] = n_epj_ar[iw];
                    }
                    for(int iw=0; iw<n_walk; iw++){
                        const S32 id_ipg = walk_grp_head + iw;
                        interaction_list_.n_disp_ep_[id_ipg+1] = interaction_list_.n_disp_ep_[id_ipg] + interaction_list_.n_ep_[id_ipg];
                    }
                    interaction_list_.adr_ep_.resizeNoInitialize( interaction_list_.n_disp_ep_[walk_grp_head+n_walk] );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(S32 i=0; i<n_thread; i++){
                        for(S32 j=0; j<adr_ipg_tmp[i].size(); j++){
                            const S32 adr_ipg = adr_ipg_tmp[i][j];
                            S32 adr_ep = interaction_list_.n_disp_ep_[adr_ipg];
                            const S32 k_ep_h = n_epj_disp_thread[i][j];
                            const S32 k_ep_e = n_epj_disp_thread[i][j+1];
                            for(S32 k=k_ep_h; k<k_ep_e; k++, adr_ep++){
                                interaction_list_.adr_ep_[adr_ep] = adr_epj_tmp[i][k];
                            }
                        }
                    }
                }

                time_profile_.calc_force__core__walk_tree += GetWtime() - offset_calc_force__core__walk_tree;



#if 1

                epj_for_force_[0].clearSize();
                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    S32 * id_epj = adr_epj_tmp[ith].getPointer(n_ep_head);
                    for(S32 jp=0; jp<n_epj_ar[iw]; jp++){
                        epj_for_force_[0].push_back( epj_sorted_[id_epj[jp]] );
                    }
                }
                S64 n_epj_cnt = 0;
                for(S32 iw=0; iw<n_walk; iw++){
                    epj_ar[iw] = epj_for_force_[0].getPointer(n_epj_cnt);
                    n_epj_cnt += n_epj_ar[iw];
                }
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar_prev.getPointer(), force_prev_array);
                } // retrieve
#else
                // original                
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar_prev.getPointer(), force_prev_array);
                } // retrieve
                epj_for_force_[0].clearSize();
                for(S32 iw=0; iw<n_walk; iw++){
                    S32 ith = iw2ith[iw];
                    S32 cnt = iw2cnt[iw];
                    S32 n_ep_head = n_epj_disp_thread[ith][cnt];
                    S32 * id_epj = adr_epj_tmp[ith].getPointer(n_ep_head);
                    for(S32 jp=0; jp<n_epj_ar[iw]; jp++){
                        epj_for_force_[0].push_back( epj_sorted_[id_epj[jp]] );
                    }
                }
                S64 n_epj_cnt = 0;
                for(S32 iw=0; iw<n_walk; iw++){
                    epj_ar[iw] = epj_for_force_[0].getPointer(n_epj_cnt);
                    n_epj_cnt += n_epj_ar[iw];
                }
#endif
                
                ret += pfunc_dispatch(tag, n_walk, 
                                      (const Tepi**)epi_ar.getPointer(), n_epi_ar.getPointer(), 
                                      (const Tepj**)epj_ar.getPointer(), n_epj_ar.getPointer());
                first_loop = false;
                for(int iw=0; iw<n_walk; iw++){
                    n_epi_ar_prev[iw] = n_epi_ar[iw];
                    force_prev_array[iw] = force_array[iw];
                }
                n_walk_prev = n_walk;
            } // end of walk group loop
            ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar_prev.getPointer(), force_prev_array);
        } // if(n_ipg > 0)
        else{
            ni_ave_ = nj_ave_ = 0;
            n_interaction_ep_ep_local_ = 0;
        }
        for(int i=0; i<n_thread; i++){
            n_interaction_ep_ep_local_ += n_interaction_ep_ep_ar[i];
        }
        time_profile_.calc_force__core += GetWtime() - offset_core;
        return ret;
    }
    //////////// Kernel:Ptcl, List:Index //////////////
    //////////////////////////////////////////////////


    ///////////////////////////////////////////////////
    //
    // FUNCTIONS OF FORCE WITHOUT WALK
    // (MUST BE USED AFTER WALK)
    // 
    ///////////////////////////////////////////////////

    ////////////////////////////////////////////////////
    //////////// Force Only, Kernel:Index, List:Index //////////////

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve, class Tsys>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalkForMultiWalk2(Tfunc_dispatch pfunc_dispatch,
                                 Tfunc_retrieve pfunc_retrieve,
                                 Tsys & sys,
                                 const S32 n_walk_limit,
                                 const bool clear){
        F64 wtime_offset = GetWtime();
        calcForceNoWalkForMultiWalk2Impl(typename TSM::force_type(),
                                         pfunc_dispatch,
                                         pfunc_retrieve,
                                         sys,
                                         n_walk_limit,
                                         clear);
        const S64 n_ipg = ipg_.size();
        n_walk_local_ += n_ipg;
        n_interaction_ep_ep_local_ = interaction_list_.n_disp_ep_[n_ipg];
        if( typeid(typename TSM::force_type) == typeid(TagForceLong)){
            n_interaction_ep_sp_local_ = interaction_list_.n_disp_sp_[n_ipg];
        }
        time_profile_.calc_force += GetWtime() - wtime_offset;
    }

    //////////// Force Only, Kernel:Index, List:Index, Force:Long //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve, class Tsys>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalkForMultiWalk2Impl(TagForceLong,
                                     Tfunc_dispatch pfunc_dispatch,
                                     Tfunc_retrieve pfunc_retrieve,
                                     Tsys & sys,
                                     const S32 n_walk_limit,
                                     const bool clear){
        F64 time_offset = GetWtime();
        S32 ret = 0;
        S32 tag = 0;
        static ReallocatableArray<Tepi*> epi_ar;
        epi_ar.resizeNoInitialize(n_walk_limit);
#if 1
        static ReallocatableArray<S32>  n_epi_ar[2];
        n_epi_ar[0].resizeNoInitialize(n_walk_limit);
        n_epi_ar[1].resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tforce*> force_ar[2];
        force_ar[0].resizeNoInitialize(n_walk_limit);
        force_ar[1].resizeNoInitialize(n_walk_limit);
#else
        static ReallocatableArray<S32>  n_epi_ar;
        n_epi_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tforce*> force_ar;
        force_ar.resizeNoInitialize(n_walk_limit);
#endif
        static ReallocatableArray<S32> n_epj_ar;
        n_epj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32*> id_epj_ar;
        id_epj_ar.resizeNoInitialize(n_walk_limit);

        static ReallocatableArray<S32> n_spj_ar;
        n_spj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32*> id_spj_ar;
        id_spj_ar.resizeNoInitialize(n_walk_limit);

        force_.resizeNoInitialize(n_loc_tot_);
        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                //force_sorted_[i].clear();
                //force_org_[i].clear();
                force_[i].clear();
            }
        }
        Tepi ** epi_dummy = NULL;
        S32 * n_epi_dummy = NULL;
        S32 ** id_epj_dummy = NULL;
        S32 *  n_epj_dummy = NULL;
        S32 ** id_spj_dummy = NULL;
        S32 *  n_spj_dummy = NULL;
        pfunc_dispatch(0, 0, (const Tepi**)epi_dummy, n_epi_dummy,
                       (const S32**)id_epj_dummy, n_epj_dummy,
                       (const S32**)id_spj_dummy, n_spj_dummy,
                       epj_sorted_.getPointer(), epj_sorted_.size(),
                       spj_sorted_.getPointer(), spj_sorted_.size(),
                       true);
        const S64 n_ipg = ipg_.size();
        if(n_ipg <= 0) return;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        static std::vector<S32> n_walk_ar;
        n_walk_ar.resize(n_loop_max); // group of walk (n_walk_ar[i] < n_walk_limit)
        static std::vector<S32> n_disp_walk_ar;
        n_disp_walk_ar.resize(n_loop_max+1);
        n_disp_walk_ar[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max + 1;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        bool first_loop = true;
        //F64 time_offset = GetWtime();
        if(n_ipg > 0){
#if 1
            S32 n_walk=-1, n_walk_prev=-1, lane_0=-1, lane_1=-1;
            S32 ip_head = 0, ip_tail=0, ip_head_prev=0, ip_tail_prev=0;
            for(int wg=0; wg<n_loop_max; wg++){
                n_walk = n_walk_ar[wg];
                const S32 n_walk_head = n_disp_walk_ar[wg];
                lane_0 = wg % 2;
                lane_1 = (wg+1) % 2;
                ip_head = ipg_[n_walk_head].adr_ptcl_;
                ip_tail = ipg_[n_walk_head+n_walk-1].adr_ptcl_ + ipg_[n_walk_head+n_walk-1].n_ptcl_;
                for(S32 iw=0; iw<n_walk; iw++){
                    const S32 id_walk = n_walk_head + iw;
                    const S32 first_adr_ip = ipg_[id_walk].adr_ptcl_; 
                    epi_ar[iw] = epi_sorted_.getPointer(first_adr_ip);
                    n_epi_ar[lane_0][iw]  = ipg_[id_walk].n_ptcl_;
                    //force_ar[lane_0][iw] = force_sorted_.getPointer(first_adr_ip);
                    force_ar[lane_0][iw] = force_.getPointer(first_adr_ip);

                    n_epj_ar[iw]  = interaction_list_.n_ep_[id_walk];
                    id_epj_ar[iw] = interaction_list_.adr_ep_.getPointer(interaction_list_.n_disp_ep_[id_walk]);

                    n_spj_ar[iw]  = interaction_list_.n_sp_[id_walk];
                    id_spj_ar[iw] = interaction_list_.adr_sp_.getPointer(interaction_list_.n_disp_sp_[id_walk]);

                    n_interaction_ep_ep_local_ += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[lane_0][iw]);
                    n_interaction_ep_sp_local_ += ((S64)n_spj_ar[iw]*(S64)n_epi_ar[lane_0][iw]);
                }
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar[lane_1].getPointer(), force_ar[lane_1].getPointer());
                }
                pfunc_dispatch(0, n_walk, 
                               (const Tepi**)epi_ar.getPointer(),   n_epi_ar[lane_0].getPointer(),
                               (const S32**)id_epj_ar.getPointer(), n_epj_ar.getPointer(),
                               (const S32**)id_spj_ar.getPointer(), n_spj_ar.getPointer(),
                               epj_sorted_.getPointer(), epj_sorted_.size(),
                               spj_sorted_.getPointer(), spj_sorted_.size(),
                               false);
                if(!first_loop){
                    F64 wtime_offset_wb_int_cp = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
                    for(S32 ip=ip_head_prev; ip<ip_tail_prev; ip++){
                        sys[ip].copyFromForce(force_[ip]);
                        sys[ip].integrate();
                        epi_org_[ip].copyFromFP( sys[ip] );
                        epj_org_[ip].copyFromFP( sys[ip] );
                    }
                    time_profile_.calc_force__wb_int_cp += GetWtime() - wtime_offset_wb_int_cp;
                }
                n_walk_prev = n_walk;
                ip_head_prev = ip_head;
                ip_tail_prev = ip_tail;
                first_loop = false;
            }
            ret += pfunc_retrieve(tag, n_walk, n_epi_ar[lane_0].getPointer(), force_ar[lane_0].getPointer());
            F64 wtime_offset_wb_int_cp = GetWtime();
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for
#endif
            for(S32 ip=ip_head_prev; ip<ip_tail_prev; ip++){
                sys[ip].copyFromForce(force_[ip]);
                sys[ip].integrate();
                epi_org_[ip].copyFromFP( sys[ip] );
                epj_org_[ip].copyFromFP( sys[ip] );
            }
            time_profile_.calc_force__wb_int_cp += GetWtime() - wtime_offset_wb_int_cp;
#else
            //original            
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = n_walk_ar[wg];
                const S32 n_walk_head = n_disp_walk_ar[wg];
                for(S32 iw=0; iw<n_walk; iw++){
                    const S32 id_walk = n_walk_head + iw;
                    const S32 first_adr_ip = ipg_[id_walk].adr_ptcl_; 
                    n_epi_ar[iw]  = ipg_[id_walk].n_ptcl_;
                    epi_ar[iw] = epi_sorted_.getPointer(first_adr_ip);
                    //force_ar[iw] = force_sorted_.getPointer(first_adr_ip);
                    force_ar[iw] = force_.getPointer(first_adr_ip);

                    n_epj_ar[iw]  = interaction_list_.n_ep_[id_walk];
                    id_epj_ar[iw] = interaction_list_.adr_ep_.getPointer(interaction_list_.n_disp_ep_[id_walk]);

                    n_spj_ar[iw]  = interaction_list_.n_sp_[id_walk];
                    id_spj_ar[iw] = interaction_list_.adr_sp_.getPointer(interaction_list_.n_disp_sp_[id_walk]);

                    //n_interaction_ep_ep_local_ += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[iw]);
                    //n_interaction_ep_sp_local_ += ((S64)n_spj_ar[iw]*(S64)n_epi_ar[iw]);
                }
                pfunc_dispatch(0, n_walk, 
                               (const Tepi**)epi_ar.getPointer(),   n_epi_ar.getPointer(),
                               (const S32**)id_epj_ar.getPointer(), n_epj_ar.getPointer(),
                               (const S32**)id_spj_ar.getPointer(), n_spj_ar.getPointer(),
                               epj_sorted_.getPointer(), epj_sorted_.size(),
                               spj_sorted_.getPointer(), spj_sorted_.size(),
                               false);
                ret += pfunc_retrieve(tag, n_walk, n_epi_ar.getPointer(), force_ar.getPointer());
            }
#endif
        }
	time_profile_.calc_force__core += GetWtime() - time_offset;
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalkForMultiWalk(Tfunc_dispatch pfunc_dispatch,
                                Tfunc_retrieve pfunc_retrieve,
                                const S32 n_walk_limit,
                                const bool clear){
        F64 wtime_offset = GetWtime();
        calcForceNoWalkForMultiWalkImpl(typename TSM::force_type(),
                                        pfunc_dispatch,
                                        pfunc_retrieve,
                                        n_walk_limit,
                                        clear);
        const S64 n_ipg = ipg_.size();
        n_walk_local_ += n_ipg;
        n_interaction_ep_ep_local_ = interaction_list_.n_disp_ep_[n_ipg];
        if( typeid(typename TSM::force_type) == typeid(TagForceLong)){
            n_interaction_ep_sp_local_ = interaction_list_.n_disp_sp_[n_ipg];
        }
        time_profile_.calc_force += GetWtime() - wtime_offset;
    }

    //////////// Force Only, Kernel:Index, List:Index, Force:Long //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalkForMultiWalkImpl(TagForceLong,
                                    Tfunc_dispatch pfunc_dispatch,
                                    Tfunc_retrieve pfunc_retrieve,
                                    const S32 n_walk_limit,
                                    const bool clear){
        F64 time_offset = GetWtime();
        S32 ret = 0;
        S32 tag = 0;
        static ReallocatableArray<Tepi*> epi_ar;
        epi_ar.resizeNoInitialize(n_walk_limit);
#if 1
        static ReallocatableArray<S32>  n_epi_ar[2];
        n_epi_ar[0].resizeNoInitialize(n_walk_limit);
        n_epi_ar[1].resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tforce*> force_ar[2];
        force_ar[0].resizeNoInitialize(n_walk_limit);
        force_ar[1].resizeNoInitialize(n_walk_limit);
#else
        static ReallocatableArray<S32>  n_epi_ar;
        n_epi_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tforce*> force_ar;
        force_ar.resizeNoInitialize(n_walk_limit);
#endif
        static ReallocatableArray<S32> n_epj_ar;
        n_epj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32*> id_epj_ar;
        id_epj_ar.resizeNoInitialize(n_walk_limit);

        static ReallocatableArray<S32> n_spj_ar;
        n_spj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32*> id_spj_ar;
        id_spj_ar.resizeNoInitialize(n_walk_limit);

        force_.resizeNoInitialize(n_loc_tot_);
        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                //force_sorted_[i].clear();
                //force_org_[i].clear();
                force_[i].clear();
            }
        }
        Tepi ** epi_dummy = NULL;
        S32 * n_epi_dummy = NULL;
        S32 ** id_epj_dummy = NULL;
        S32 *  n_epj_dummy = NULL;
        S32 ** id_spj_dummy = NULL;
        S32 *  n_spj_dummy = NULL;
        pfunc_dispatch(0, 0, (const Tepi**)epi_dummy, n_epi_dummy,
                       (const S32**)id_epj_dummy, n_epj_dummy,
                       (const S32**)id_spj_dummy, n_spj_dummy,
                       epj_sorted_.getPointer(), epj_sorted_.size(),
                       spj_sorted_.getPointer(), spj_sorted_.size(),
                       true);
        const S64 n_ipg = ipg_.size();
        if(n_ipg <= 0) return;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        static std::vector<S32> n_walk_ar;
        n_walk_ar.resize(n_loop_max); // group of walk (n_walk_ar[i] < n_walk_limit)
        static std::vector<S32> n_disp_walk_ar;
        n_disp_walk_ar.resize(n_loop_max+1);
        n_disp_walk_ar[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max + 1;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        bool first_loop = true;
        //F64 time_offset = GetWtime();
        if(n_ipg > 0){
#if 1
            S32 n_walk=-1, n_walk_prev=-1, lane_0=-1, lane_1=-1;
            for(int wg=0; wg<n_loop_max; wg++){
                n_walk = n_walk_ar[wg];
                const S32 n_walk_head = n_disp_walk_ar[wg];
                lane_0 = wg % 2;
                lane_1 = (wg+1) % 2;
                for(S32 iw=0; iw<n_walk; iw++){
                    const S32 id_walk = n_walk_head + iw;
                    const S32 first_adr_ip = ipg_[id_walk].adr_ptcl_; 
                    epi_ar[iw] = epi_sorted_.getPointer(first_adr_ip);
                    n_epi_ar[lane_0][iw]  = ipg_[id_walk].n_ptcl_;
                    //force_ar[lane_0][iw] = force_sorted_.getPointer(first_adr_ip);
                    force_ar[lane_0][iw] = force_.getPointer(first_adr_ip);

                    n_epj_ar[iw]  = interaction_list_.n_ep_[id_walk];
                    id_epj_ar[iw] = interaction_list_.adr_ep_.getPointer(interaction_list_.n_disp_ep_[id_walk]);

                    n_spj_ar[iw]  = interaction_list_.n_sp_[id_walk];
                    id_spj_ar[iw] = interaction_list_.adr_sp_.getPointer(interaction_list_.n_disp_sp_[id_walk]);

                    n_interaction_ep_ep_local_ += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[lane_0][iw]);
                    n_interaction_ep_sp_local_ += ((S64)n_spj_ar[iw]*(S64)n_epi_ar[lane_0][iw]);
                }
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar[lane_1].getPointer(), force_ar[lane_1].getPointer());
                }
                pfunc_dispatch(0, n_walk, 
                               (const Tepi**)epi_ar.getPointer(),   n_epi_ar[lane_0].getPointer(),
                               (const S32**)id_epj_ar.getPointer(), n_epj_ar.getPointer(),
                               (const S32**)id_spj_ar.getPointer(), n_spj_ar.getPointer(),
                               epj_sorted_.getPointer(), epj_sorted_.size(),
                               spj_sorted_.getPointer(), spj_sorted_.size(),
                               false);
                n_walk_prev = n_walk;
                first_loop = false;
            }
            ret += pfunc_retrieve(tag, n_walk, n_epi_ar[lane_0].getPointer(), force_ar[lane_0].getPointer());
#else
            //original            
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = n_walk_ar[wg];
                const S32 n_walk_head = n_disp_walk_ar[wg];
                for(S32 iw=0; iw<n_walk; iw++){
                    const S32 id_walk = n_walk_head + iw;
                    const S32 first_adr_ip = ipg_[id_walk].adr_ptcl_; 
                    n_epi_ar[iw]  = ipg_[id_walk].n_ptcl_;
                    epi_ar[iw] = epi_sorted_.getPointer(first_adr_ip);
                    //force_ar[iw] = force_sorted_.getPointer(first_adr_ip);
                    force_ar[iw] = force_.getPointer(first_adr_ip);

                    n_epj_ar[iw]  = interaction_list_.n_ep_[id_walk];
                    id_epj_ar[iw] = interaction_list_.adr_ep_.getPointer(interaction_list_.n_disp_ep_[id_walk]);

                    n_spj_ar[iw]  = interaction_list_.n_sp_[id_walk];
                    id_spj_ar[iw] = interaction_list_.adr_sp_.getPointer(interaction_list_.n_disp_sp_[id_walk]);

                    //n_interaction_ep_ep_local_ += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[iw]);
                    //n_interaction_ep_sp_local_ += ((S64)n_spj_ar[iw]*(S64)n_epi_ar[iw]);
                }
                pfunc_dispatch(0, n_walk, 
                               (const Tepi**)epi_ar.getPointer(),   n_epi_ar.getPointer(),
                               (const S32**)id_epj_ar.getPointer(), n_epj_ar.getPointer(),
                               (const S32**)id_spj_ar.getPointer(), n_spj_ar.getPointer(),
                               epj_sorted_.getPointer(), epj_sorted_.size(),
                               spj_sorted_.getPointer(), spj_sorted_.size(),
                               false);
                ret += pfunc_retrieve(tag, n_walk, n_epi_ar.getPointer(), force_ar.getPointer());
            }
#endif
        }
	time_profile_.calc_force__core += GetWtime() - time_offset;
    }

    //////////// Force Only, Kernel:Index, List:Index, Force:Short //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalkForMultiWalkImpl(TagForceShort,
                                    Tfunc_dispatch pfunc_dispatch,
                                    Tfunc_retrieve pfunc_retrieve,
                                    const S32 n_walk_limit,
                                    const bool clear){
        F64 time_offset = GetWtime();
        S32 ret = 0;
        S32 tag = 0;
        static ReallocatableArray<Tepi*> epi_ar;
        epi_ar.resizeNoInitialize(n_walk_limit);
#if 1
        // overlape version
        static ReallocatableArray<S32>  n_epi_ar[2];
        n_epi_ar[0].resizeNoInitialize(n_walk_limit);
        n_epi_ar[1].resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tforce*> force_ar[2];
        force_ar[0].resizeNoInitialize(n_walk_limit);
        force_ar[1].resizeNoInitialize(n_walk_limit);
#else
        static ReallocatableArray<S32>  n_epi_ar;
        n_epi_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tforce*> force_ar;
        force_ar.resizeNoInitialize(n_walk_limit);
#endif
        static ReallocatableArray<S32> n_epj_ar;
        n_epj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<S32*> id_epj_ar;
        id_epj_ar.resizeNoInitialize(n_walk_limit);

        force_.resizeNoInitialize(n_loc_tot_);
        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                //force_sorted_[i].clear();
                //force_org_[i].clear();
                force_[i].clear();
            }
        }
        Tepi ** epi_dummy = NULL;
        S32 * n_epi_dummy = NULL;
        S32 ** id_epj_dummy = NULL;
        S32 *  n_epj_dummy = NULL;

        pfunc_dispatch(0, 0, (const Tepi**)epi_dummy, n_epi_dummy,
                       (const S32**)id_epj_dummy, n_epj_dummy,
                       epj_sorted_.getPointer(), epj_sorted_.size(),
                       true);
        const S64 n_ipg = ipg_.size();
        if(n_ipg <= 0) return;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        static std::vector<S32> n_walk_ar;
        n_walk_ar.resize(n_loop_max); // group of walk (n_walk_ar[i] < n_walk_limit)
        static std::vector<S32> n_disp_walk_ar;
        n_disp_walk_ar.resize(n_loop_max+1);
        n_disp_walk_ar[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max + 1;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        bool first_loop = true;
        //F64 time_offset = GetWtime();
        if(n_ipg > 0){
#if 1
            // overlape version
            S32 n_walk=-1, n_walk_prev=-1, lane_0=-1, lane_1=-1;
            for(int wg=0; wg<n_loop_max; wg++){
                n_walk = n_walk_ar[wg];
                const S32 n_walk_head = n_disp_walk_ar[wg];
                lane_0 = wg % 2;
                lane_1 = (wg+1) % 2;
                for(S32 iw=0; iw<n_walk; iw++){
                    const S32 id_walk = n_walk_head + iw;
                    const S32 first_adr_ip = ipg_[id_walk].adr_ptcl_;
                    epi_ar[iw] = epi_sorted_.getPointer(first_adr_ip);
                    n_epi_ar[lane_0][iw]  = ipg_[id_walk].n_ptcl_;
                    //force_ar[lane_0][iw] = force_sorted_.getPointer(first_adr_ip);
                    force_ar[lane_0][iw] = force_.getPointer(first_adr_ip);

                    n_epj_ar[iw]  = interaction_list_.n_ep_[id_walk];
                    id_epj_ar[iw] = interaction_list_.adr_ep_.getPointer(interaction_list_.n_disp_ep_[id_walk]);

                    n_interaction_ep_ep_local_ += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[lane_0][iw]);
                }
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar[lane_1].getPointer(), force_ar[lane_1].getPointer());
                }
                pfunc_dispatch(0, n_walk, 
                               (const Tepi**)epi_ar.getPointer(),   n_epi_ar[lane_0].getPointer(),
                               (const S32**)id_epj_ar.getPointer(), n_epj_ar.getPointer(),
                               epj_sorted_.getPointer(), epj_sorted_.size(),
                               false);
                n_walk_prev = n_walk;
                first_loop = false;
            }
            ret += pfunc_retrieve(tag, n_walk, n_epi_ar[lane_0].getPointer(), force_ar[lane_0].getPointer());
#else
            // original
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = n_walk_ar[wg];
                const S32 n_walk_head = n_disp_walk_ar[wg];
                for(S32 iw=0; iw<n_walk; iw++){
                    const S32 id_walk = n_walk_head + iw;
                    const S32 first_adr_ip = ipg_[id_walk].adr_ptcl_; 
                    n_epi_ar[iw]  = ipg_[id_walk].n_ptcl_;
                    epi_ar[iw] = epi_sorted_.getPointer(first_adr_ip);
                    //force_ar[iw] = force_sorted_.getPointer(first_adr_ip);
                    force_ar[iw] = force_.getPointer(first_adr_ip);

                    n_epj_ar[iw]  = interaction_list_.n_ep_[id_walk];
                    id_epj_ar[iw] = interaction_list_.adr_ep_.getPointer(interaction_list_.n_disp_ep_[id_walk]);

                    n_interaction_ep_ep_local_ += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[iw]);
                }
                pfunc_dispatch(0, n_walk, 
                               (const Tepi**)epi_ar.getPointer(),   n_epi_ar.getPointer(),
                               (const S32**)id_epj_ar.getPointer(), n_epj_ar.getPointer(),
                               epj_sorted_.getPointer(), epj_sorted_.size(),
                               false);
                ret += pfunc_retrieve(tag, n_walk, n_epi_ar.getPointer(), force_ar.getPointer());
            }
#endif
        }
        time_profile_.calc_force__core += GetWtime() - time_offset;
    }

    ///////////////////////////////////////////////////////////////
    //////////// Force Only, Kernel:Ptcl, List:Index //////////////

    //////////// Force Only, Kernel:Ptcl, List:Index, Force:Long //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalkForMultiWalkNewImpl(TagForceLong,
                                       Tfunc_dispatch pfunc_dispatch,
                                       Tfunc_retrieve pfunc_retrieve,
                                       const S32 n_walk_limit,
                                       const bool clear){

        //F64 time_offset = GetWtime();
        S32 ret = 0;
        S32 tag = 0;
        static ReallocatableArray<Tepi*> epi_ar;
        epi_ar.resizeNoInitialize(n_walk_limit);
#if 1
        static ReallocatableArray<S32>  n_epi_ar[2];
        n_epi_ar[0].resizeNoInitialize(n_walk_limit);
        n_epi_ar[1].resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tforce*> force_ar[2];
        force_ar[0].resizeNoInitialize(n_walk_limit);
        force_ar[1].resizeNoInitialize(n_walk_limit);
#else
        static ReallocatableArray<S32>  n_epi_ar;
        n_epi_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tforce*> force_ar;
        force_ar.resizeNoInitialize(n_walk_limit);
#endif
        
        static ReallocatableArray<S32> n_epj_ar;
        n_epj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tepj*> epj_ar;
        epj_ar.resizeNoInitialize(n_walk_limit);

        static ReallocatableArray<S32> n_spj_ar;
        n_spj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tspj*> spj_ar;
        spj_ar.resizeNoInitialize(n_walk_limit);

        force_.resizeNoInitialize(n_loc_tot_);
        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                //force_sorted_[i].clear();
                //force_org_[i].clear();
                force_[i].clear();
            }
        }
        const S64 n_ipg = ipg_.size();
        if(n_ipg <= 0) return;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        static std::vector<S32> n_walk_ar;
        n_walk_ar.resize(n_loop_max); // group of walk (n_walk_ar[i] < n_walk_limit)
        static std::vector<S32> n_disp_walk_ar;
        n_disp_walk_ar.resize(n_loop_max+1);
        n_disp_walk_ar[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max + 1;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        bool first_loop = true;
        if(n_ipg > 0){
#if 1
            S32 n_walk=-1, n_walk_prev=-1, lane_0=-1, lane_1=-1;
            for(int wg=0; wg<n_loop_max; wg++){
                n_walk = n_walk_ar[wg];
                n_walk_prev = (wg>0) ? n_walk_ar[wg-1] : 0;
                const S32 n_walk_head = n_disp_walk_ar[wg];
                lane_0 = wg % 2;
                lane_1 = (wg+1) % 2;
                for(S32 iw=0; iw<n_walk; iw++){
                    const S32 id_walk = n_walk_head + iw;
                    const S32 first_adr_ip = ipg_[id_walk].adr_ptcl_;
                    epi_ar[iw] = epi_sorted_.getPointer(first_adr_ip);
                    n_epi_ar[lane_0][iw]  = ipg_[id_walk].n_ptcl_;
                    //force_ar[lane_0][iw] = force_sorted_.getPointer(first_adr_ip);
                    force_ar[lane_0][iw] = force_.getPointer(first_adr_ip);
                    n_epj_ar[iw]  = interaction_list_.n_ep_[id_walk];
                    n_spj_ar[iw]  = interaction_list_.n_sp_[id_walk];
                    n_interaction_ep_ep_local_ += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[lane_0][iw]);
                    n_interaction_ep_sp_local_ += ((S64)n_spj_ar[iw]*(S64)n_epi_ar[lane_0][iw]);
                }

                const S64 n_ep_head = interaction_list_.n_disp_ep_[n_walk_head];
                const S64 n_ep_end  = interaction_list_.n_disp_ep_[n_walk_head+n_walk];
                const S64 n_epj_tot = n_ep_end - n_ep_head;
                epj_for_force_[0].resizeNoInitialize(n_epj_tot);
                for(S32 jp=0; jp<n_epj_tot; jp++){
                    epj_for_force_[0][jp] = epj_sorted_[ interaction_list_.adr_ep_[jp+n_ep_head] ];
                }
                const S64 n_sp_head = interaction_list_.n_disp_sp_[n_walk_head];
                const S64 n_sp_end  = interaction_list_.n_disp_sp_[n_walk_head+n_walk];
                const S64 n_spj_tot = n_sp_end - n_sp_head;
                spj_for_force_[0].resizeNoInitialize(n_spj_tot);
                for(S32 jp=0; jp<n_spj_tot; jp++){
                    spj_for_force_[0][jp] = spj_sorted_[ interaction_list_.adr_sp_[jp+n_sp_head] ];
                }
                S64 n_epj_cnt = 0;
                S64 n_spj_cnt = 0;
                epj_ar.resizeNoInitialize(n_walk);
                spj_ar.resizeNoInitialize(n_walk);
                for(S32 iw=0; iw<n_walk; iw++){
                    epj_ar[iw] = epj_for_force_[0].getPointer(n_epj_cnt);
                    spj_ar[iw] = spj_for_force_[0].getPointer(n_spj_cnt);
                    n_epj_cnt += n_epj_ar[iw];
                    n_spj_cnt += n_spj_ar[iw];
                }
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar[lane_1].getPointer(), force_ar[lane_1].getPointer());
                }
                ret += pfunc_dispatch(tag, n_walk, 
                                      (const Tepi**)epi_ar.getPointer(), n_epi_ar[lane_0].getPointer(), 
                                      (const Tepj**)epj_ar.getPointer(), n_epj_ar.getPointer(), 
                                      (const Tspj**)spj_ar.getPointer(), n_spj_ar.getPointer());
                first_loop = false;
            }
            ret += pfunc_retrieve(tag, n_walk, n_epi_ar[lane_0].getPointer(), force_ar[lane_0].getPointer());
#else
            // original
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = n_walk_ar[wg];
                const S32 n_walk_head = n_disp_walk_ar[wg];
                for(S32 iw=0; iw<n_walk; iw++){
                    const S32 id_walk = n_walk_head + iw;
                    const S32 first_adr_ip = ipg_[id_walk].adr_ptcl_; 
                    n_epi_ar[iw]  = ipg_[id_walk].n_ptcl_;
                    epi_ar[iw] = epi_sorted_.getPointer(first_adr_ip);
                    //force_ar[iw] = force_sorted_.getPointer(first_adr_ip);
                    force_ar[iw] = force_.getPointer(first_adr_ip);

                    n_epj_ar[iw]  = interaction_list_.n_ep_[id_walk];
                    n_spj_ar[iw]  = interaction_list_.n_sp_[id_walk];

                    n_interaction_ep_ep_local_ += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[iw]);
                    n_interaction_ep_sp_local_ += ((S64)n_spj_ar[iw]*(S64)n_epi_ar[iw]);
                }

                const S64 n_ep_head = interaction_list_.n_disp_ep_[n_walk_head];
                const S64 n_ep_end  = interaction_list_.n_disp_ep_[n_walk_head+n_walk];
                const S64 n_epj_tot = n_ep_end - n_ep_head;
                epj_for_force_[0].resizeNoInitialize(n_epj_tot);
                for(S32 jp=0; jp<n_epj_tot; jp++){
                    epj_for_force_[0][jp] = epj_sorted_[ interaction_list_.adr_ep_[jp+n_ep_head] ];
                }
                const S64 n_sp_head = interaction_list_.n_disp_sp_[n_walk_head];
                const S64 n_sp_end  = interaction_list_.n_disp_sp_[n_walk_head+n_walk];
                const S64 n_spj_tot = n_sp_end - n_sp_head;
                spj_for_force_[0].resizeNoInitialize(n_spj_tot);
                for(S32 jp=0; jp<n_spj_tot; jp++){
                    spj_for_force_[0][jp] = spj_sorted_[ interaction_list_.adr_sp_[jp+n_sp_head] ];
                }
                S64 n_epj_cnt = 0;
                S64 n_spj_cnt = 0;
                epj_ar.resizeNoInitialize(n_walk);
                spj_ar.resizeNoInitialize(n_walk);
                for(S32 iw=0; iw<n_walk; iw++){
                    epj_ar[iw] = epj_for_force_[0].getPointer(n_epj_cnt);
                    spj_ar[iw] = spj_for_force_[0].getPointer(n_spj_cnt);
                    n_epj_cnt += n_epj_ar[iw];
                    n_spj_cnt += n_spj_ar[iw];
                }                
                ret += pfunc_dispatch(tag, n_walk, 
                                      (const Tepi**)epi_ar.getPointer(), n_epi_ar.getPointer(), 
                                      (const Tepj**)epj_ar.getPointer(), n_epj_ar.getPointer(), 
                                      (const Tspj**)spj_ar.getPointer(), n_spj_ar.getPointer());

                ret += pfunc_retrieve(tag, n_walk, n_epi_ar.getPointer(), force_ar.getPointer());
            }
#endif
        }
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalkForMultiWalkNew(Tfunc_dispatch pfunc_dispatch,
                                   Tfunc_retrieve pfunc_retrieve,
                                   const S32 n_walk_limit,
                                   const bool clear){
        F64 wtime_offset = GetWtime();
        calcForceNoWalkForMultiWalkNewImpl(typename TSM::force_type(),
                                           pfunc_dispatch,
                                           pfunc_retrieve,
                                           n_walk_limit,
                                           clear);
        const S32 n_ipg = ipg_.size();
        n_walk_local_ += n_ipg;
        n_interaction_ep_ep_local_ = interaction_list_.n_disp_ep_[n_ipg];
        if( typeid(typename TSM::force_type) == typeid(TagForceLong)){
            n_interaction_ep_sp_local_ = interaction_list_.n_disp_sp_[n_ipg];
        }
        time_profile_.calc_force += GetWtime() - wtime_offset;
    }

    //////////// Force Only, Kernel:Ptcl, List:Index, Force:Short //////////////
    template<class TSM, class Tforce, class Tepi, class Tepj,
	     class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceNoWalkForMultiWalkNewImpl(TagForceShort,
                                       Tfunc_dispatch pfunc_dispatch,
                                       Tfunc_retrieve pfunc_retrieve,
                                       const S32 n_walk_limit,
                                       const bool clear){
        F64 time_offset = GetWtime();
        S32 ret = 0;
        S32 tag = 0;
        static ReallocatableArray<Tepi*> epi_ar;
        epi_ar.resizeNoInitialize(n_walk_limit);
        
#if 1
        static ReallocatableArray<S32>  n_epi_ar[2];
        n_epi_ar[0].resizeNoInitialize(n_walk_limit);
        n_epi_ar[1].resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tforce*> force_ar[2];
        force_ar[0].resizeNoInitialize(n_walk_limit);
        force_ar[1].resizeNoInitialize(n_walk_limit);
#else
        static ReallocatableArray<S32>  n_epi_ar;
        n_epi_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tforce*> force_ar;
        force_ar.resizeNoInitialize(n_walk_limit);
#endif
        
        static ReallocatableArray<S32> n_epj_ar;
        n_epj_ar.resizeNoInitialize(n_walk_limit);
        static ReallocatableArray<Tepj*> epj_ar;
        epj_ar.resizeNoInitialize(n_walk_limit);

        //force_sorted_.resizeNoInitialize(n_loc_tot_);
        //force_org_.resizeNoInitialize(n_loc_tot_);
        force_.resizeNoInitialize(n_loc_tot_);
        if(clear){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_loc_tot_; i++){
                //force_sorted_[i].clear();
                //force_org_[i].clear();
                force_[i].clear();
            }
        }
        const S64 n_ipg = ipg_.size();
        if(n_ipg <= 0) return;
        const S32 n_loop_max = n_ipg/n_walk_limit + (n_ipg%n_walk_limit==0 ? 0 : 1);
        static std::vector<S32> n_walk_ar;
        n_walk_ar.resize(n_loop_max); // group of walk (n_walk_ar[i] < n_walk_limit)
        static std::vector<S32> n_disp_walk_ar;
        n_disp_walk_ar.resize(n_loop_max+1);
        n_disp_walk_ar[0] = 0;
        for(int wg=0; wg<n_ipg%n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max + 1;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        for(int wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
            n_walk_ar[wg] = n_ipg / n_loop_max;
            n_disp_walk_ar[wg+1] = n_disp_walk_ar[wg] + n_walk_ar[wg];
        }
        bool first_loop = true;
        if(n_ipg > 0){
#if 1
            S32 n_walk=-1, n_walk_prev=-1, lane_0=-1, lane_1=-1;
            for(int wg=0; wg<n_loop_max; wg++){
                n_walk = n_walk_ar[wg];
                n_walk_prev = (wg>0) ? n_walk_ar[wg-1] : 0;
                const S32 n_walk_head = n_disp_walk_ar[wg];
                lane_0 = wg % 2;
                lane_1 = (wg+1) % 2;
                for(S32 iw=0; iw<n_walk; iw++){
                    const S32 id_walk = n_walk_head + iw;
                    const S32 first_adr_ip = ipg_[id_walk].adr_ptcl_;
                    epi_ar[iw] = epi_sorted_.getPointer(first_adr_ip);
                    n_epi_ar[lane_0][iw]  = ipg_[id_walk].n_ptcl_;
                    //force_ar[lane_0][iw] = force_sorted_.getPointer(first_adr_ip);
                    force_ar[lane_0][iw] = force_.getPointer(first_adr_ip);
                    n_epj_ar[iw]  = interaction_list_.n_ep_[id_walk];
                    n_interaction_ep_ep_local_ += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[lane_0][iw]);
                }
                const S64 n_ep_head = interaction_list_.n_disp_ep_[n_walk_head];
                const S64 n_ep_end  = interaction_list_.n_disp_ep_[n_walk_head+n_walk];
                const S64 n_epj_tot = n_ep_end - n_ep_head;
                epj_for_force_[0].resizeNoInitialize(n_epj_tot);
                for(S32 jp=0; jp<n_epj_tot; jp++){
                    epj_for_force_[0][jp] = epj_sorted_[ interaction_list_.adr_ep_[jp+n_ep_head] ];
                }
                S64 n_epj_cnt = 0;
                epj_ar.resizeNoInitialize(n_walk);
                for(S32 iw=0; iw<n_walk; iw++){
                    epj_ar[iw] = epj_for_force_[0].getPointer(n_epj_cnt);
                    n_epj_cnt += n_epj_ar[iw];
                }
                if(!first_loop){
                    ret += pfunc_retrieve(tag, n_walk_prev, n_epi_ar[lane_1].getPointer(), force_ar[lane_1].getPointer());
                }
                ret += pfunc_dispatch(tag, n_walk, 
                                      (const Tepi**)epi_ar.getPointer(), n_epi_ar[lane_0].getPointer(), 
                                      (const Tepj**)epj_ar.getPointer(), n_epj_ar.getPointer());
                first_loop = false;
            }
            ret += pfunc_retrieve(tag, n_walk, n_epi_ar[lane_0].getPointer(), force_ar[lane_0].getPointer());
#else
            // original
            for(int wg=0; wg<n_loop_max; wg++){
                const S32 n_walk = n_walk_ar[wg];
                const S32 n_walk_head = n_disp_walk_ar[wg];
                for(S32 iw=0; iw<n_walk; iw++){
                    const S32 id_walk = n_walk_head + iw;
                    const S32 first_adr_ip = ipg_[id_walk].adr_ptcl_; 
                    n_epi_ar[iw]  = ipg_[id_walk].n_ptcl_;
                    epi_ar[iw] = epi_sorted_.getPointer(first_adr_ip);
                    //force_ar[iw] = force_sorted_.getPointer(first_adr_ip);
                    force_ar[iw] = force_.getPointer(first_adr_ip);

                    n_epj_ar[iw]  = interaction_list_.n_ep_[id_walk];

                    n_interaction_ep_ep_local_ += ((S64)n_epj_ar[iw]*(S64)n_epi_ar[iw]);
                }
                const S64 n_ep_head = interaction_list_.n_disp_ep_[n_walk_head];
                const S64 n_ep_end  = interaction_list_.n_disp_ep_[n_walk_head+n_walk];
                const S64 n_epj_tot = n_ep_end - n_ep_head;
                epj_for_force_[0].resizeNoInitialize(n_epj_tot);
                for(S32 jp=0; jp<n_epj_tot; jp++){
                    epj_for_force_[0][jp] = epj_sorted_[ interaction_list_.adr_ep_[jp+n_ep_head] ];
                }
                S64 n_epj_cnt = 0;
                epj_ar.resizeNoInitialize(n_walk);
                for(S32 iw=0; iw<n_walk; iw++){
                    epj_ar[iw] = epj_for_force_[0].getPointer(n_epj_cnt);
                    n_epj_cnt += n_epj_ar[iw];
                }
                ret += pfunc_dispatch(tag, n_walk, 
                                      (const Tepi**)epi_ar.getPointer(), n_epi_ar.getPointer(), 
                                      (const Tepj**)epj_ar.getPointer(), n_epj_ar.getPointer()); 
                ret += pfunc_retrieve(tag, n_walk, n_epi_ar.getPointer(), force_ar.getPointer());
            }
#endif
        }
	time_profile_.calc_force__core += GetWtime() - time_offset;
    }
}
