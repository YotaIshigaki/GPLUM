namespace ParticleSimulator{


    inline F64ort GetCorrespondingTreeCell(const F64ort & box){
        return MortonKey::getCorrespondingTreeCell(box);
        /*
        const U64 low  = MortonKey::getKey(box.low_);
        const U64 high = MortonKey::getKey(box.high_);
        const U64 tmp = high^low;
        S32 lev = 0;
        S32 id = 0;
        for(S32 i=TREE_LEVEL_LIMIT-1; i>=0; i--){
            id = (tmp >> i*3) & 0x7;
            if(id != 0){
                lev = TREE_LEVEL_LIMIT-i-1;
                break;
            }
        }
        return lev;
        */
    }
    
    ///////////////////////////////
    // small functions for dispatch
    template<class Ttc>
    inline F64ort GetOuterBoundaryOfMyTree(TagSearchLong,
                                           const ReallocatableArray<Ttc> & tc_first){
        return F64ort(-1234.5, -9876.5);
    }
    template<class Ttc>
    inline F64ort GetOuterBoundaryOfMyTree(TagSearchLongCutoff,
                                           const ReallocatableArray<Ttc> & tc_first){
        return tc_first[0].mom_.getVertexOut();
    }
    template<class Ttc>
    inline F64ort GetOuterBoundaryOfMyTree(TagSearchLongScatter,
                                           const ReallocatableArray<Ttc> & tc_first){
        return tc_first[0].mom_.getVertexOut();
    }
    template<class Ttc>
    inline F64ort GetOuterBoundaryOfMyTree(TagSearchLongSymmetry,
                                           const ReallocatableArray<Ttc> & tc_first){
        return tc_first[0].mom_.getVertexOut();
    }
    template<class Ttc>
    inline F64ort GetOuterBoundaryOfMyTree(TagSearchShortGather,
                                           const ReallocatableArray<Ttc> & tc_first){
        return tc_first[0].mom_.getVertexOut();
    }
    template<class Ttc>
    inline F64ort GetOuterBoundaryOfMyTree(TagSearchShortSymmetry,
                                           const ReallocatableArray<Ttc> & tc_first){
        return tc_first[0].mom_.getVertexOut();
    }
    template<class Ttc>
    inline F64ort GetOuterBoundaryOfMyTree(TagSearchShortScatter,
                                           const ReallocatableArray<Ttc> & tc_first){
        return tc_first[0].mom_.getVertexOut();
    }

    inline void ExchangeOuterBoundary(TagSearchLong,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){
#ifdef USE_ALLGATHER_EXLET
        Comm::allGather(my_outer, 1, outer.getPointer());
#endif
        /* do nothing */
    }
    inline void ExchangeOuterBoundary(TagSearchLongCutoff,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){
#ifdef USE_ALLGATHER_EXLET
        Comm::allGather(my_outer, 1, outer.getPointer());
#endif
        /* do nothing */
    }
    inline void ExchangeOuterBoundary(TagSearchLongScatter,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){
#ifdef USE_ALLGATHER_EXLET
        Comm::allGather(my_outer, 1, outer.getPointer());
#endif
        /* do nothing */
    }
    inline void ExchangeOuterBoundary(TagSearchLongSymmetry,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){
        Comm::allGather(my_outer, 1, outer.getPointer());
    }
    inline void ExchangeOuterBoundary(TagSearchShortScatter,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){
#ifdef USE_ALLGATHER_EXLET
        Comm::allGather(my_outer, 1, outer.getPointer());
#endif
        /* do nothing */ }
    inline void ExchangeOuterBoundary(TagSearchShortGather,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){
#ifdef USE_ALLGATHER_EXLET
        Comm::allGather(my_outer, 1, outer.getPointer());
#endif
        /* do nothing */ }
    inline void ExchangeOuterBoundary(TagSearchShortSymmetry,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){
#ifdef USE_ALLGATHER_EXLET
        Comm::allGather(my_outer, 1, outer.getPointer());
#endif
        /* do nothing */ }


    

    inline void ExchangeBoundary(TagSearchLong,
                                 F64ort * my_inner,
                                 F64ort * my_outer,
                                 ReallocatableArray<F64ort> & inner,
                                 ReallocatableArray<F64ort> & outer){
#ifdef USE_ALLGATHER_EXLET
        Comm::allGather(my_outer, 1, outer.getPointer());
#endif
        /* do nothing */
    }
    inline void ExchangeBoundary(TagSearchLongCutoff,
                                 F64ort * my_inner,
                                 F64ort * my_outer,
                                 ReallocatableArray<F64ort> & inner,
                                 ReallocatableArray<F64ort> & outer){
#ifdef USE_ALLGATHER_EXLET
        Comm::allGather(my_outer, 1, outer.getPointer());
#endif
        /* do nothing */
    }
    inline void ExchangeBoundary(TagSearchLongScatter,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){
#ifdef USE_ALLGATHER_EXLET
        Comm::allGather(my_outer, 1, outer.getPointer());
#endif
        /* do nothing */
    }
    inline void ExchangeBoundary(TagSearchLongSymmetry,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){
        Comm::allGather(my_outer, 1, outer.getPointer());
    }
    inline void ExchangeBoundary(TagSearchShortScatter,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){
#ifdef USE_ALLGATHER_EXLET
        Comm::allGather(my_outer, 1, outer.getPointer());
#endif
        /* do nothing */ }
    inline void ExchangeBoundary(TagSearchShortGather,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){
#ifdef USE_ALLGATHER_EXLET
        Comm::allGather(my_outer, 1, outer.getPointer());
#endif
        /* do nothing */ }
    inline void ExchangeBoundary(TagSearchShortSymmetry,
                                      F64ort * my_outer,
                                      ReallocatableArray<F64ort> & outer){
#ifdef USE_ALLGATHER_EXLET
        Comm::allGather(my_outer, 1, outer.getPointer());
#endif
        /* do nothing */ }

    

    ////////////////////
    // for exchange LET
    // only for open boundary
    template<class Tptcl, class Ttc>
    inline void CopyPtclFromTreeToSendBuf(TagSearchBoundaryConditionOpenOnly,
                                          ReallocatableArray<Tptcl> & ptcl_send,
                                          const ReallocatableArray<Ttc> & tc,
                                          const ReallocatableArray<S32> & adr_ptcl_send,
                                          const S32 n_ptcl,
                                          const S32 n_ptcl_offset,
                                          const F64vec & shift){
        for(S32 j=0; j<n_ptcl; j++){
            S32 adr = adr_ptcl_send[n_ptcl_offset+j];
            ptcl_send[n_ptcl_offset+j].copyFromMoment(tc[adr].mom_);
        }
    }

    // for periodic boundary (not only for periodic but also open boudnary)
    template<class Tptcl, class Ttc>
    inline void CopyPtclFromTreeToSendBuf(TagSearchBoundaryConditionOpenPeriodic,
                                          ReallocatableArray<Tptcl> & ptcl_send,
                                          const ReallocatableArray<Ttc> & tc,
                                          const ReallocatableArray<S32> & adr_ptcl_send,
                                          const S32 n_ptcl,
                                          const S32 n_ptcl_offset,
                                          const F64vec & shift){
        for(S32 j=0; j<n_ptcl; j++){
            S32 adr = adr_ptcl_send[n_ptcl_offset+j];
            ptcl_send[n_ptcl_offset+j].copyFromMoment(tc[adr].mom_);
            const F64vec pos_new = ptcl_send[n_ptcl_offset+j].getPos() - shift;
            ptcl_send[n_ptcl_offset+j].setPos(pos_new);
        }
    }
    template<class Tptcl>
    inline void CopyPtclToSendBuf(TagSearchBoundaryConditionOpenOnly,
                                  ReallocatableArray<Tptcl> & ptcl_send,
                                  const ReallocatableArray<Tptcl> & ptcl,
                                  const ReallocatableArray<S32> & adr_ptcl_send,
                                  const S32 n_ptcl,
                                  const S32 n_ptcl_offset,
                                  const F64vec & shift){
        for(S32 j=0; j<n_ptcl; j++){
            S32 adr = adr_ptcl_send[n_ptcl_offset+j];
            ptcl_send[n_ptcl_offset+j] = ptcl[adr];
        }
    }
    template<class Tptcl>
    inline void CopyPtclToSendBuf(TagSearchBoundaryConditionOpenPeriodic,
                                  ReallocatableArray<Tptcl> & ptcl_send,
                                  const ReallocatableArray<Tptcl> & ptcl,
                                  const ReallocatableArray<S32> & adr_ptcl_send,
                                  const S32 n_ptcl,
                                  const S32 n_ptcl_offset,
                                  const F64vec & shift){
        for(S32 j=0; j<n_ptcl; j++){
            S32 adr = adr_ptcl_send[n_ptcl_offset+j];
            ptcl_send[n_ptcl_offset+j] = ptcl[adr];
            ptcl_send[n_ptcl_offset+j].setPos(ptcl[adr].getPos() - shift);
        }
    }
    // for exchange LET
    ////////////////////
    
    // small functions for dispatch
    ///////////////////////////////


    ////////////////
    // for long mode
#ifdef USE_ALLGATHER_EXLET
    template<class TSM, class Ttc, class Ttp, class Tep, class Tsp, class Twalkmode>
    inline void FindScatterParticle(const ReallocatableArray<Ttc> & tc_first,
                                    const ReallocatableArray<Ttp> & tp_first,
                                    const ReallocatableArray<Tep> & ep_first,
                                    ReallocatableArray<S32> & n_ep_send,
                                    ReallocatableArray<S32> & adr_ep_send,
                                    const DomainInfo & dinfo,
                                    const S32 n_leaf_limit,
                                    ReallocatableArray<S32> & n_sp_send,
                                    ReallocatableArray<S32> & adr_sp_send,
                                    ReallocatableArray<F64vec> & shift_per_image,
                                    ReallocatableArray<S32> & n_image_per_proc,
                                    ReallocatableArray<S32> & n_ep_per_image,
                                    ReallocatableArray<S32> & n_sp_per_image,
                                    ReallocatableArray<S32> & rank_send,
                                    ReallocatableArray<S32> & rank_recv,
                                    ReallocatableArray<Tsp> & top_moment,
                                    const F64ort & inner_boundary,
                                    const F64 theta,
                                    const F64 r_crit_sq,
                                    ReallocatableArray<S32> & rank_recv_a2a){
        const auto n_thread = Comm::getNumberOfThread();
        const auto n_proc = Comm::getNumberOfProc();
        const auto my_rank = Comm::getRank();
        static auto first = true;
        static ReallocatableArray<S32> * rank_tmp;
        static ReallocatableArray<F64vec> * shift_per_image_tmp;
        static ReallocatableArray<S32> * adr_ep_send_tmp;
        static ReallocatableArray<S32> * n_ep_per_image_tmp;
        static ReallocatableArray<S32> * adr_sp_send_tmp;
        static ReallocatableArray<S32> * n_sp_per_image_tmp;
        static ReallocatableArray<F64ort> outer_boundary_of_tree;
        static ReallocatableArray<F64ort> minimum_tree_boundary_of_local_tree;
        static ReallocatableArray<S32> * rank_send_tmp;
        static ReallocatableArray<S32> * rank_recv_tmp;
        static ReallocatableArray<S32> * rank_recv_a2a_tmp;
        if(first){
            rank_tmp = new ReallocatableArray<S32>[n_thread];
            shift_per_image_tmp = new ReallocatableArray<F64vec>[n_thread];
            adr_ep_send_tmp = new ReallocatableArray<S32>[n_thread];
            n_ep_per_image_tmp = new ReallocatableArray<S32>[n_thread];
            adr_sp_send_tmp = new ReallocatableArray<S32>[n_thread];
            n_sp_per_image_tmp = new ReallocatableArray<S32>[n_thread];
            rank_send_tmp = new ReallocatableArray<S32>[n_thread];
            rank_recv_tmp = new ReallocatableArray<S32>[n_thread];
            rank_recv_a2a_tmp = new ReallocatableArray<S32>[n_thread];
            first = false;
        }
        outer_boundary_of_tree.resizeNoInitialize(n_proc);
        minimum_tree_boundary_of_local_tree.resizeNoInitialize(n_proc);
        n_ep_send.resizeNoInitialize(n_proc);
        n_sp_send.resizeNoInitialize(n_proc);
        n_image_per_proc.resizeNoInitialize(n_proc);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_proc; i++){
            n_ep_send[i] = n_sp_send[i] = n_image_per_proc[i] = 0;
        }
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        F64ort pos_root_domain = dinfo.getPosRootDomain();
        F64vec len_peri = pos_root_domain.getFullLength();
        for(S32 i=0; i<DIMENSION; i++){
            if(pa[i]==false) len_peri[i] = 0.0;
        }
        auto outer_boundary_of_my_tree = GetOuterBoundaryOfMyTree(typename TSM::search_type(), tc_first);
        auto minimum_tree_boundary_of_my_local_tree = GetCorrespondingTreeCell(inner_boundary);
        assert( minimum_tree_boundary_of_my_local_tree.contained(inner_boundary) );
        if(typeid(TSM) == typeid(SEARCH_MODE_LONG)){
#if 0
            outer_boundary_of_my_tree = minimum_tree_boundary_of_my_local_tree;
#else
            outer_boundary_of_my_tree = inner_boundary;
#endif
            //std::cout<<"my_rank= "<<my_rank<<std::endl;
            //std::cout<<"inner_boundary= "<<inner_boundary<<std::endl;
            //std::cout<<"GetCorrespondingTreeCell(inner_boundary)= "
            //         <<GetCorrespondingTreeCell(inner_boundary)
            //         <<std::endl;
#if 0
            outer_boundary_of_my_tree.low_ = F64vec(-999999.9);
            outer_boundary_of_my_tree.high_ = F64vec(999999.9);
#endif
        }

        ExchangeOuterBoundary(typename TSM::search_type(), &outer_boundary_of_my_tree, outer_boundary_of_tree);
        Comm::allGather(&minimum_tree_boundary_of_my_local_tree, 1, minimum_tree_boundary_of_local_tree.getPointer());

        //if(my_rank==0){
        //    for(S32 i=0; i<n_proc; i++){
        //        std::cout<<"i= "<<i
        //                 <<" outer_boundary_of_tree= "
        //                 <<outer_boundary_of_tree[i]
        //                 <<std::endl;
        //    }
        //}
        
        F64 r_crit_send = minimum_tree_boundary_of_local_tree[my_rank].getFullLength().getMax();
        if(theta > 0.0){
            r_crit_send = r_crit_send*r_crit_send / (theta*theta);
        }
        else{
            r_crit_send = -1.0;
        }

        /*
        if(my_rank==23){
            std::cout<<"minimum_tree_boundary_of_local_tree[23]= "<<minimum_tree_boundary_of_local_tree[23]<<std::endl;
            std::cout<<"tc_first[0].mom_.pos= "<<tc_first[0].mom_.pos<<std::endl;
            for(S32 i=8; i<16; i++){
                std::cout<<"i= "<<i<<" tc_first[i].mom_.pos= "<<tc_first[i].mom_.pos<<std::endl;
            }
        }
        Comm::barrier();
        */
        
PS_OMP_PARALLEL
        {
            const S32 ith = Comm::getThreadNum();
            const S32 n_thread = Comm::getNumberOfThread();
            const S32 head = (n_proc/n_thread)*ith + std::min(n_proc%n_thread, ith);
            const S32 end  = (n_proc/n_thread)*(ith+1) + std::min(n_proc%n_thread, (ith+1));
            rank_send_tmp[ith].clearSize();
            rank_recv_tmp[ith].clearSize();
            rank_recv_a2a_tmp[ith].clearSize();
            for(S32 ib=head; ib<end; ib++){
                const S32 id0 = std::min(ib, my_rank);
                const S32 id1 = std::max(ib, my_rank);
                bool flag_recv_a2a = true;

                /*
                if(my_rank==23 && ib==40){
                    std::cout<<"outer_boundary_of_tree[23]= "<<outer_boundary_of_tree[23]<<std::endl;
                    std::cout<<"outer_boundary_of_tree[40]= "<<outer_boundary_of_tree[40]<<std::endl;
                    std::cout<<"(GetDistanceMinSq(outer_boundary_of_tree[id0], outer_boundary_of_tree[id1], len_peri)= "
                             <<(GetDistanceMinSq(outer_boundary_of_tree[id0], outer_boundary_of_tree[id1], len_peri))
                             <<std::endl;
                }
                */
                
                if( (GetDistanceMinSq(outer_boundary_of_tree[id0], outer_boundary_of_tree[id1], len_peri) <= 0.0)
                    || (GetDistanceMinSq(dinfo.getPosDomain(id0), outer_boundary_of_tree[id1], len_peri) <= 0.0)
                    || (GetDistanceMinSq(outer_boundary_of_tree[id0], dinfo.getPosDomain(id1), len_peri) <= 0.0) ){
                    //|| (GetDistanceMinSq(dinfo.getPosDomain(id0), dinfo.getPosDomain(id1), len_peri) <= 0.0)){
                    
                    //if(Comm::getRank()==0){
                    //    std::cerr<<"id0= "<<id0<<" outer_boundary_of_tree[id0]= "<<outer_boundary_of_tree[id0]<<std::endl;
                    //    std::cerr<<"id1= "<<id1<<" outer_boundary_of_tree[id1]= "<<outer_boundary_of_tree[id1]<<std::endl;
                    //}
                    
                    rank_send_tmp[ith].push_back(ib);
                    rank_recv_tmp[ith].push_back(ib);
                    flag_recv_a2a = false;
                }
                else{

                    /*
                    if(my_rank==23 && ib==40){
                        std::cout<<"r_crit_send= "<<r_crit_send<<std::endl;
                        std::cout<<"minimum_tree_boundary_of_local_tree[ib]= "<<minimum_tree_boundary_of_local_tree[ib]<<std::endl;
                        std::cout<<"top_moment[my_rank].getPos()= "<<top_moment[my_rank].getPos()<<std::endl;
                        std::cout<<"(GetDistanceMinSq(minimum_tree_boundary_of_local_tree[ib], top_moment[my_rank].getPos(), len_peri)= "
                                 <<(GetDistanceMinSq(minimum_tree_boundary_of_local_tree[ib], top_moment[my_rank].getPos(), len_peri))
                                 <<std::endl;
                    }
                    */
                    
                    //if(GetDistanceMinSq(outer_boundary_of_tree[ib], top_moment[my_rank].getPos(), len_peri) <= r_crit_send){
                    /*
                    if( (GetDistanceMinSq(minimum_tree_boundary_of_local_tree[ib], top_moment[my_rank].getPos(), len_peri) <= r_crit_send)
                        || (GetDistanceMinSq(dinfo.getPosDomain(ib), top_moment[my_rank].getPos(), len_peri) <= r_crit_send) 
                        || (GetDistanceMinSq(outer_boundary_of_tree[ib], top_moment[my_rank].getPos(), len_peri) <= r_crit_send)
                        || (GetDistanceMinSq(dinfo.getPosDomain(ib), minimum_tree_boundary_of_local_tree[my_rank], len_peri) <= 0.0)){
                    */
                    if(    (GetDistanceMinSq(dinfo.getPosDomain(ib), top_moment[my_rank].getPos(), len_peri) <= r_crit_send) 
                        || (GetDistanceMinSq(outer_boundary_of_tree[ib], top_moment[my_rank].getPos(), len_peri) <= r_crit_send)
                        || (GetDistanceMinSq(dinfo.getPosDomain(ib), minimum_tree_boundary_of_local_tree[my_rank], len_peri) <= 0.0)){
                    
                        rank_send_tmp[ith].push_back(ib);
                    }

                    //F64 r_crit_recv = outer_boundary_of_tree[ib].getFullLength().getMax();
                    F64 r_crit_recv = minimum_tree_boundary_of_local_tree[ib].getFullLength().getMax();
                    if(theta > 0.0){
                        r_crit_recv = r_crit_recv*r_crit_recv / (theta*theta);
                    }
                    else{
                        r_crit_recv = -1.0;
                    }
                    //if(GetDistanceMinSq(outer_boundary_of_tree[my_rank], top_moment[ib].getPos(), len_peri) <= r_crit_recv){
                    /*
                    if( (GetDistanceMinSq(minimum_tree_boundary_of_local_tree[my_rank], top_moment[ib].getPos(), len_peri) <= r_crit_recv)
                        || (GetDistanceMinSq(dinfo.getPosDomain(my_rank), top_moment[ib].getPos(), len_peri) <= r_crit_recv)
                        || (GetDistanceMinSq(outer_boundary_of_tree[my_rank], top_moment[ib].getPos(), len_peri) <= r_crit_recv)
                        || (GetDistanceMinSq(dinfo.getPosDomain(my_rank), minimum_tree_boundary_of_local_tree[ib], len_peri) <= 0.0)){
                    */
                    if(    (GetDistanceMinSq(dinfo.getPosDomain(my_rank), top_moment[ib].getPos(), len_peri) <= r_crit_recv)
                        || (GetDistanceMinSq(outer_boundary_of_tree[my_rank], top_moment[ib].getPos(), len_peri) <= r_crit_recv)
                        || (GetDistanceMinSq(dinfo.getPosDomain(my_rank), minimum_tree_boundary_of_local_tree[ib], len_peri) <= 0.0)){                    
                        rank_recv_tmp[ith].push_back(ib);
                        flag_recv_a2a = false;
                    }
                }
                if(flag_recv_a2a){
                    rank_recv_a2a_tmp[ith].push_back(ib);
                }
            }
        } // end of OMP scope

        Comm::barrier();
        //exit(1);
        
        PackData(&rank_send_tmp[0], Comm::getNumberOfThread(), rank_send);
        PackData(&rank_recv_tmp[0], Comm::getNumberOfThread(), rank_recv);
        PackData(&rank_recv_a2a_tmp[0], Comm::getNumberOfThread(), rank_recv_a2a);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_proc; i++){
            n_image_per_proc[i] = 0;
        }
        S32 adr_tc = 0;
        S32 adr_tree_sp_first = 0;
PS_OMP_PARALLEL
        {
            S32 ith = Comm::getThreadNum();
            rank_tmp[ith].clearSize();
            shift_per_image_tmp[ith].clearSize();
            adr_ep_send_tmp[ith].clearSize();
            n_ep_per_image_tmp[ith].clearSize();
            adr_sp_send_tmp[ith].clearSize();
            n_sp_per_image_tmp[ith].clearSize();
            ReallocatableArray<Tsp> sp_first;
PS_OMP(omp for schedule(dynamic, 4))
            for(S32 ib=0; ib<rank_send.size(); ib++){
                const S32 rank = rank_send[ib];
                const S32 n_image_tmp_prev = shift_per_image_tmp[ith].size();
                rank_tmp[ith].push_back(rank);
                CalcNumberAndShiftOfImageDomain
                    (shift_per_image_tmp[ith],  dinfo.getPosRootDomain().getFullLength(),
                     outer_boundary_of_my_tree, dinfo.getPosDomain(rank), pa, false);
                const S32 n_image_tmp = shift_per_image_tmp[ith].size();
                n_image_per_proc[rank] = n_image_tmp - n_image_tmp_prev;
                S32 n_ep_prev = adr_ep_send_tmp[ith].size();
                S32 n_sp_prev = adr_sp_send_tmp[ith].size();
                for(S32 j=n_image_tmp_prev; j<n_image_tmp; j++){
                    S32 n_ep_prev_2 = adr_ep_send_tmp[ith].size();
                    S32 n_sp_prev_2 = adr_sp_send_tmp[ith].size();
                    if(my_rank==rank && j==n_image_tmp_prev){
                        n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
                        n_sp_per_image_tmp[ith].push_back(adr_sp_send_tmp[ith].size() - n_sp_prev_2);
                        continue;
                    }
                    F64ort pos_target_domain = dinfo.getPosDomain(rank).shift(shift_per_image_tmp[ith][j]);
                    TargetBox<TSM> target_box;
                    SetTargetBoxExLet(target_box, pos_target_domain, outer_boundary_of_tree[rank].shift(shift_per_image_tmp[ith][j]));
                    MakeListUsingTreeRecursiveTop
                        <TSM, Ttc, TreeParticle, Tep, Tsp, Twalkmode, TagChopLeafFalse, TagCopyInfoCloseNoSp>
                        (tc_first, adr_tc, tp_first,
                         ep_first, adr_ep_send_tmp[ith],
                         sp_first, adr_sp_send_tmp[ith],
                         target_box,
                         r_crit_sq, n_leaf_limit,
                         adr_tree_sp_first, F64vec(0.0));
                    n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
                    n_sp_per_image_tmp[ith].push_back(adr_sp_send_tmp[ith].size() - n_sp_prev_2);
                }

                n_ep_send[rank] = adr_ep_send_tmp[ith].size() - n_ep_prev;
                n_sp_send[rank] = adr_sp_send_tmp[ith].size() - n_sp_prev;
            }
        } // end of OMP scope
        static ReallocatableArray<S32> n_disp_image_per_proc;
        n_disp_image_per_proc.resizeNoInitialize(n_proc+1);
        n_disp_image_per_proc[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_image_per_proc[i+1] = n_disp_image_per_proc[i] + n_image_per_proc[i];
        }
        const S32 n_image_tot = n_disp_image_per_proc[n_proc];
        shift_per_image.resizeNoInitialize( n_image_tot );
        n_ep_per_image.resizeNoInitialize( n_image_tot);
        n_sp_per_image.resizeNoInitialize( n_image_tot);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_thread; i++){
            S32 n_cnt = 0;
            for(S32 j=0; j<rank_tmp[i].size(); j++){
                S32 rank = rank_tmp[i][j];
                S32 offset = n_disp_image_per_proc[rank];
                for(S32 k=0; k<n_image_per_proc[rank]; k++){
                    shift_per_image[offset+k] = shift_per_image_tmp[i][n_cnt];
                    n_ep_per_image[offset+k]  = n_ep_per_image_tmp[i][n_cnt];
                    n_sp_per_image[offset+k]  = n_sp_per_image_tmp[i][n_cnt];
                    n_cnt++;
                }
            }
        }
        static ReallocatableArray<S32> n_disp_ep_per_image;
        n_disp_ep_per_image.resizeNoInitialize(n_image_tot+1);
        static ReallocatableArray<S32> n_disp_sp_per_image;
        n_disp_sp_per_image.resizeNoInitialize(n_image_tot+1);
        n_disp_ep_per_image[0] = 0;
        n_disp_sp_per_image[0] = 0;
        for(S32 i=0; i<n_image_tot; i++){
            n_disp_ep_per_image[i+1] = n_disp_ep_per_image[i] + n_ep_per_image[i];
            n_disp_sp_per_image[i+1] = n_disp_sp_per_image[i] + n_sp_per_image[i];
        }
        const S32 n_ep_send_tot = n_disp_ep_per_image[ n_image_tot ];
        const S32 n_sp_send_tot = n_disp_sp_per_image[ n_image_tot ];
        adr_ep_send.resizeNoInitialize( n_ep_send_tot );
        adr_sp_send.resizeNoInitialize( n_sp_send_tot );
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_thread; i++){
            S32 n_cnt_ep = 0;
            S32 n_cnt_sp = 0;
            for(S32 j=0; j<rank_tmp[i].size(); j++){
                S32 rank = rank_tmp[i][j];
                const S32 adr_image_head = n_disp_image_per_proc[rank];
                const S32 adr_image_end = n_disp_image_per_proc[rank+1];
                for(S32 k=adr_image_head; k<adr_image_end; k++){
                    const S32 adr_ep_head = n_disp_ep_per_image[k];
                    const S32 adr_ep_end = n_disp_ep_per_image[k+1];
                    for(S32 l=adr_ep_head; l<adr_ep_end; l++){
                        adr_ep_send[l] = adr_ep_send_tmp[i][n_cnt_ep++];
                    }
                    const S32 adr_sp_head = n_disp_sp_per_image[k];
                    const S32 adr_sp_end  = n_disp_sp_per_image[k+1];
                    for(S32 l=adr_sp_head; l<adr_sp_end; l++){
                        adr_sp_send[l] = adr_sp_send_tmp[i][n_cnt_sp++];
                    }
                }
            }
        }
    }
#else    
    template<class TSM, class Ttc, class Ttp, class Tep, class Tsp, class Twalkmode>
    inline void FindScatterParticle(const ReallocatableArray<Ttc> & tc_first,
                                    const ReallocatableArray<Ttp> & tp_first,
                                    const ReallocatableArray<Tep> & ep_first,
                                    ReallocatableArray<S32> & n_ep_send,
                                    ReallocatableArray<S32> & adr_ep_send,
                                    const DomainInfo & dinfo,
                                    const S32 n_leaf_limit,
                                    ReallocatableArray<S32> & n_sp_send,
                                    ReallocatableArray<S32> & adr_sp_send,
                                    ReallocatableArray<F64vec> & shift_per_image,
                                    ReallocatableArray<S32> & n_image_per_proc,
                                    ReallocatableArray<S32> & n_ep_per_image,
                                    ReallocatableArray<S32> & n_sp_per_image,
                                    const F64 r_crit_sq){
        const S32 n_thread = Comm::getNumberOfThread();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();        
        static bool first = true;
        static ReallocatableArray<S32> * rank_tmp;
        static ReallocatableArray<F64vec> * shift_per_image_tmp;
        static ReallocatableArray<S32> * adr_ep_send_tmp;
        static ReallocatableArray<S32> * n_ep_per_image_tmp;
        static ReallocatableArray<S32> * adr_sp_send_tmp;
        static ReallocatableArray<S32> * n_sp_per_image_tmp;
        static ReallocatableArray<F64ort> outer_boundary_of_tree;
        if(first){
            rank_tmp = new ReallocatableArray<S32>[n_thread];
            shift_per_image_tmp = new ReallocatableArray<F64vec>[n_thread];
            adr_ep_send_tmp = new ReallocatableArray<S32>[n_thread];
            n_ep_per_image_tmp = new ReallocatableArray<S32>[n_thread];
            adr_sp_send_tmp = new ReallocatableArray<S32>[n_thread];
            n_sp_per_image_tmp = new ReallocatableArray<S32>[n_thread];
            first = false;
        }
        outer_boundary_of_tree.resizeNoInitialize(n_proc);
        n_ep_send.resizeNoInitialize(n_proc);
        n_sp_send.resizeNoInitialize(n_proc);
        n_image_per_proc.resizeNoInitialize(n_proc);
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        F64ort pos_root_domain = dinfo.getPosRootDomain();
        F64vec len_peri = pos_root_domain.getFullLength();
        for(S32 i=0; i<DIMENSION; i++){
            if(pa[i]==false) len_peri[i] = 0.0;
        }
        F64ort outer_boundary_of_my_tree = GetOuterBoundaryOfMyTree(typename TSM::search_type(), tc_first);
        S32 adr_tc = 0;
        S32 adr_tree_sp_first = 0;
        ExchangeOuterBoundary(typename TSM::search_type(), &outer_boundary_of_my_tree, outer_boundary_of_tree);
        // for long symmetry
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            S32 ith = Comm::getThreadNum();
            rank_tmp[ith].clearSize();
            shift_per_image_tmp[ith].clearSize();
            adr_ep_send_tmp[ith].clearSize();
            n_ep_per_image_tmp[ith].clearSize();
            adr_sp_send_tmp[ith].clearSize();
            n_sp_per_image_tmp[ith].clearSize();
            ReallocatableArray<Tsp> sp_first;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4)
#endif
            for(S32 i=0; i<n_proc; i++){
                const S32 n_image_tmp_prev = shift_per_image_tmp[ith].size();
                rank_tmp[ith].push_back(i);
                CalcNumberAndShiftOfImageDomain
                    (shift_per_image_tmp[ith],  dinfo.getPosRootDomain().getFullLength(),
                     outer_boundary_of_my_tree, dinfo.getPosDomain(i), pa, false);
                const S32 n_image_tmp = shift_per_image_tmp[ith].size();
                n_image_per_proc[i] = n_image_tmp - n_image_tmp_prev;
                S32 n_ep_prev = adr_ep_send_tmp[ith].size();
                S32 n_sp_prev = adr_sp_send_tmp[ith].size();
                for(S32 j=n_image_tmp_prev; j<n_image_tmp; j++){
                    S32 n_ep_prev_2 = adr_ep_send_tmp[ith].size();
                    S32 n_sp_prev_2 = adr_sp_send_tmp[ith].size();
                    if(my_rank==i && j==n_image_tmp_prev){
                        n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2); // is 0
                        n_sp_per_image_tmp[ith].push_back(adr_sp_send_tmp[ith].size() - n_sp_prev_2); // is 0
                        continue;
                    }
                    F64ort pos_target_domain = dinfo.getPosDomain(i).shift(shift_per_image_tmp[ith][j]);
                    TargetBox<TSM> target_box;
                    SetTargetBoxExLet(target_box, pos_target_domain, outer_boundary_of_tree[i].shift(shift_per_image_tmp[ith][j]));
                    MakeListUsingTreeRecursiveTop
                        <TSM, Ttc, TreeParticle, Tep, Tsp, Twalkmode, TagChopLeafFalse, TagCopyInfoCloseNoSp>
                        (tc_first, adr_tc, tp_first,
                         ep_first, adr_ep_send_tmp[ith],
                         sp_first, adr_sp_send_tmp[ith],
                         target_box,
                         r_crit_sq, n_leaf_limit,
                         adr_tree_sp_first, F64vec(0.0));
                    n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
                    n_sp_per_image_tmp[ith].push_back(adr_sp_send_tmp[ith].size() - n_sp_prev_2);
                }
                n_ep_send[i] = adr_ep_send_tmp[ith].size() - n_ep_prev;
                n_sp_send[i] = adr_sp_send_tmp[ith].size() - n_sp_prev;
            }
        } // end of OMP scope

        ReallocatableArray<S32> n_disp_image_per_proc;
        n_disp_image_per_proc.resizeNoInitialize(n_proc+1);
        n_disp_image_per_proc[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_image_per_proc[i+1] = n_disp_image_per_proc[i] + n_image_per_proc[i];
        }
        const S32 n_image_tot = n_disp_image_per_proc[n_proc];
        shift_per_image.resizeNoInitialize( n_image_tot );
        n_ep_per_image.resizeNoInitialize( n_image_tot);
        n_sp_per_image.resizeNoInitialize( n_image_tot);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_thread; i++){
            S32 n_cnt = 0;
            for(S32 j=0; j<rank_tmp[i].size(); j++){
                S32 rank = rank_tmp[i][j];
                S32 offset = n_disp_image_per_proc[rank];
                for(S32 k=0; k<n_image_per_proc[rank]; k++){
                    shift_per_image[offset+k] = shift_per_image_tmp[i][n_cnt];
                    n_ep_per_image[offset+k] = n_ep_per_image_tmp[i][n_cnt];
                    n_sp_per_image[offset+k] = n_sp_per_image_tmp[i][n_cnt];
                    n_cnt++;
                }
            }
        }
        ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot+1, n_image_tot+1, 1);
        ReallocatableArray<S32> n_disp_sp_per_image(n_image_tot+1, n_image_tot+1, 1);
        n_disp_ep_per_image[0] = 0;
        n_disp_sp_per_image[0] = 0;
        for(S32 i=0; i<n_image_tot; i++){
            n_disp_ep_per_image[i+1] = n_disp_ep_per_image[i] + n_ep_per_image[i];
            n_disp_sp_per_image[i+1] = n_disp_sp_per_image[i] + n_sp_per_image[i];
        }
        const S32 n_ep_send_tot = n_disp_ep_per_image[ n_image_tot ];
        const S32 n_sp_send_tot = n_disp_sp_per_image[ n_image_tot ];
        adr_ep_send.resizeNoInitialize( n_ep_send_tot );
        adr_sp_send.resizeNoInitialize( n_sp_send_tot );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_thread; i++){
            S32 n_cnt_ep = 0;
            S32 n_cnt_sp = 0;
            for(S32 j=0; j<rank_tmp[i].size(); j++){
                S32 rank = rank_tmp[i][j];
                const S32 adr_image_head = n_disp_image_per_proc[rank];
                const S32 adr_image_end = n_disp_image_per_proc[rank+1];
                for(S32 k=adr_image_head; k<adr_image_end; k++){
                    const S32 adr_ep_head = n_disp_ep_per_image[k];
                    const S32 adr_ep_end = n_disp_ep_per_image[k+1];
                    for(S32 l=adr_ep_head; l<adr_ep_end; l++){
                        adr_ep_send[l] = adr_ep_send_tmp[i][n_cnt_ep++];
                    }
                    const S32 adr_sp_head = n_disp_sp_per_image[k];
                    const S32 adr_sp_end  = n_disp_sp_per_image[k+1];
                    for(S32 l=adr_sp_head; l<adr_sp_end; l++){
                        adr_sp_send[l] = adr_sp_send_tmp[i][n_cnt_sp++];
                    }
                }
            }
        }
    }
#endif

    /////////////////
    // for short mode
#ifdef USE_ALLGATHER_EXLET
    template<class Ttc, class Ttp, class Tep, class Tsp, class Twalkmode>
    inline void FindScatterParticle
    (const ReallocatableArray<Ttc> & tc_first,
     const ReallocatableArray<Ttp> & tp_first,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<S32> & n_ep_send,
     ReallocatableArray<S32> & adr_ep_send,
     const DomainInfo & dinfo,
     const S32 n_leaf_limit,
     ReallocatableArray<F64vec> & shift_per_image,
     ReallocatableArray<S32> & n_image_per_proc,
     ReallocatableArray<S32> & n_ep_per_image,
     ReallocatableArray<S32> & rank_send,
     ReallocatableArray<S32> & rank_recv){
        const auto n_thread = Comm::getNumberOfThread();
        const auto n_proc = Comm::getNumberOfProc();
        const auto my_rank = Comm::getRank();
        rank_send.resizeNoInitialize(n_proc);
        rank_recv.resizeNoInitialize(n_proc);
        n_ep_send.resizeNoInitialize(n_proc);
        n_image_per_proc.resizeNoInitialize(n_proc);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_proc; i++){
            rank_send[i] = rank_recv[i] = n_ep_send[i] = n_image_per_proc[i] = 0;
        }
        static auto first = true;
        static ReallocatableArray<S32> * rank_tmp;
        static ReallocatableArray<F64vec> * shift_per_image_tmp;
        static ReallocatableArray<S32> * adr_ep_send_tmp;
        static ReallocatableArray<S32> * n_ep_per_image_tmp;
        static ReallocatableArray<S32> * rank_send_tmp;
        static ReallocatableArray<S32> * rank_recv_tmp;
        static ReallocatableArray<F64ort> outer_boundary_of_tree;
        if(first){
            rank_tmp = new ReallocatableArray<S32>[n_thread];
            shift_per_image_tmp = new ReallocatableArray<F64vec>[n_thread];
            adr_ep_send_tmp = new ReallocatableArray<S32>[n_thread];
            n_ep_per_image_tmp = new ReallocatableArray<S32>[n_thread];
            rank_send_tmp = new ReallocatableArray<S32>[n_thread];
            rank_recv_tmp = new ReallocatableArray<S32>[n_thread];
            first = false;
        }
        outer_boundary_of_tree.resizeNoInitialize(n_proc);
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        const auto pos_root_domain = dinfo.getPosRootDomain();
        auto len_peri = pos_root_domain.getFullLength();
        for(S32 i=0; i<DIMENSION; i++){
            if(pa[i]==false) len_peri[i] = 0.0;
        }
        const auto outer_boundary_of_my_tree = tc_first[0].mom_.vertex_out_;
        Comm::allGather(&outer_boundary_of_my_tree, 1, outer_boundary_of_tree.getPointer());
PS_OMP_PARALLEL
        {
            const auto ith = Comm::getThreadNum();
            const auto n_thread = Comm::getNumberOfThread();
            const auto head = (n_proc/n_thread)*ith + std::min(n_proc%n_thread, ith);
            const auto end  = (n_proc/n_thread)*(ith+1) + std::min(n_proc%n_thread, (ith+1));
            rank_send_tmp[ith].clearSize();
            rank_recv_tmp[ith].clearSize();
            for(S32 ib=head; ib<end; ib++){
                const S32 id0 = std::min(ib, my_rank);
                const S32 id1 = std::max(ib, my_rank);
                if(GetDistanceMinSq(outer_boundary_of_tree[id0], outer_boundary_of_tree[id1], len_peri) <= 0.0){
                    rank_send_tmp[ith].push_back(ib);
                    rank_recv_tmp[ith].push_back(ib);
                }
            }
        } // end of OMP scope
        PackData(&rank_send_tmp[0], Comm::getNumberOfThread(), rank_send);
        PackData(&rank_recv_tmp[0], Comm::getNumberOfThread(), rank_recv);

PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_proc; i++){
            n_image_per_proc[i] = 0;
        }
        S32 adr_tc = 0;
        S32 adr_tree_sp_first = 0;
PS_OMP_PARALLEL
        {
            auto ith = Comm::getThreadNum();
            rank_tmp[ith].clearSize();
            shift_per_image_tmp[ith].clearSize();
            adr_ep_send_tmp[ith].clearSize();
            n_ep_per_image_tmp[ith].clearSize();
            ReallocatableArray<S32> adr_sp_send_tmp;
            ReallocatableArray<Tsp> sp_first;
            F64 r_crit_sq = 1.0; // dummy
PS_OMP(omp for schedule(dynamic, 4))
            for(S32 ib=0; ib<rank_send.size(); ib++){
                const auto rank = rank_send[ib];
                const auto n_image_tmp_prev = shift_per_image_tmp[ith].size();
                rank_tmp[ith].push_back(rank);
                CalcNumberAndShiftOfImageDomain
                    (shift_per_image_tmp[ith], dinfo.getPosRootDomain().getFullLength(),
                     outer_boundary_of_my_tree, dinfo.getPosDomain(rank), pa, false);
                const auto n_image_tmp = shift_per_image_tmp[ith].size();
                n_image_per_proc[rank] = n_image_tmp - n_image_tmp_prev;
                const auto n_ep_prev = adr_ep_send_tmp[ith].size();
                for(S32 j=n_image_tmp_prev; j<n_image_tmp; j++){
                    const auto n_ep_prev_2 = adr_ep_send_tmp[ith].size();
                    if(my_rank==rank && j==n_image_tmp_prev){
                        n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
                        continue;
                    }
                    const auto pos_target_domain = dinfo.getPosDomain(rank).shift(shift_per_image_tmp[ith][j]);
                    TargetBox<SEARCH_MODE_SCATTER> target_box;
                    target_box.vertex_in_ = pos_target_domain;
                    MakeListUsingTreeRecursiveTop
                        <SEARCH_MODE_SCATTER, Ttc, TreeParticle, Tep, Tsp, Twalkmode, TagChopLeafFalse, TagCopyInfoCloseNoSp>
                        (tc_first, adr_tc, tp_first,
                         ep_first, adr_ep_send_tmp[ith],
                         sp_first, adr_sp_send_tmp,
                         target_box,
                         r_crit_sq, n_leaf_limit,
                         adr_tree_sp_first, F64vec(0.0));
                    n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
                }
                n_ep_send[rank] = adr_ep_send_tmp[ith].size() - n_ep_prev;
            }
        } // end of OMP scope
        ReallocatableArray<S32> n_disp_image_per_proc(n_proc+1, n_proc+1, 1);
        n_disp_image_per_proc[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_image_per_proc[i+1] = n_disp_image_per_proc[i] + n_image_per_proc[i];
        }
        const auto n_image_tot = n_disp_image_per_proc[n_proc];
        shift_per_image.resizeNoInitialize( n_image_tot );
        n_ep_per_image.resizeNoInitialize( n_image_tot);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_thread; i++){
            auto n_cnt = 0;
            for(S32 j=0; j<rank_tmp[i].size(); j++){
                auto rank = rank_tmp[i][j];
                auto offset = n_disp_image_per_proc[rank];
                for(S32 k=0; k<n_image_per_proc[rank]; k++){
                    shift_per_image[offset+k] = shift_per_image_tmp[i][n_cnt];
                    n_ep_per_image[offset+k] = n_ep_per_image_tmp[i][n_cnt];
                    n_cnt++;
                }
            }
        }
        ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot+1, n_image_tot+1, 1);
        n_disp_ep_per_image[0] = 0;
        for(S32 i=0; i<n_image_tot; i++){
            n_disp_ep_per_image[i+1] = n_disp_ep_per_image[i] + n_ep_per_image[i];
        }
        const auto n_ep_send_tot = n_disp_ep_per_image[ n_image_tot ];
        adr_ep_send.resizeNoInitialize( n_ep_send_tot );
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_thread; i++){
            auto n_cnt_ep = 0;
            for(S32 j=0; j<rank_tmp[i].size(); j++){
                auto rank = rank_tmp[i][j];
                const auto adr_image_head = n_disp_image_per_proc[rank];
                const auto adr_image_end = n_disp_image_per_proc[rank+1];
                for(S32 k=adr_image_head; k<adr_image_end; k++){
                    const auto adr_ep_head = n_disp_ep_per_image[k];
                    const auto adr_ep_end = n_disp_ep_per_image[k+1];
                    for(S32 l=adr_ep_head; l<adr_ep_end; l++){
                        adr_ep_send[l] = adr_ep_send_tmp[i][n_cnt_ep++];
                    }
                }
            }
        }
    }
#else
    template<class Ttc, class Ttp, class Tep, class Tsp, class Twalkmode>
    inline void FindScatterParticle
    (const ReallocatableArray<Ttc> & tc_first,
     const ReallocatableArray<Ttp> & tp_first,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<S32> & n_ep_send,
     ReallocatableArray<S32> & adr_ep_send,
     const DomainInfo & dinfo,
     const S32 n_leaf_limit,
     ReallocatableArray<F64vec> & shift_per_image,
     ReallocatableArray<S32> & n_image_per_proc,
     ReallocatableArray<S32> & n_ep_per_image){
        static bool first = true;
        static const S32 n_thread = Comm::getNumberOfThread();
        static ReallocatableArray<S32> * rank_tmp;
        static ReallocatableArray<F64vec> * shift_per_image_tmp;
        static ReallocatableArray<S32> * adr_ep_send_tmp;
        static ReallocatableArray<S32> * n_ep_per_image_tmp;
        if(first){
            rank_tmp = new ReallocatableArray<S32>[n_thread];
            shift_per_image_tmp = new ReallocatableArray<F64vec>[n_thread];
            adr_ep_send_tmp = new ReallocatableArray<S32>[n_thread];
            n_ep_per_image_tmp = new ReallocatableArray<S32>[n_thread];
            first = false;
        }
        const auto n_proc = Comm::getNumberOfProc();
        const auto my_rank = Comm::getRank();
        n_ep_send.resizeNoInitialize(n_proc);
        n_image_per_proc.resizeNoInitialize(n_proc);
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        auto pos_root_domain = dinfo.getPosRootDomain();
        auto len_peri = pos_root_domain.getFullLength();
        for(S32 i=0; i<DIMENSION; i++){
            if(pa[i]==false) len_peri[i] = 0.0;
        }
        const auto outer_boundary_of_my_tree = tc_first[0].mom_.vertex_out_;
        auto adr_tc = 0;
        auto adr_tree_sp_first = 0;
PS_OMP_PARALLEL
        {
            auto ith = Comm::getThreadNum();
            rank_tmp[ith].clearSize();
            shift_per_image_tmp[ith].clearSize();
            adr_ep_send_tmp[ith].clearSize();
            n_ep_per_image_tmp[ith].clearSize();
            ReallocatableArray<S32> adr_sp_send_tmp;
            ReallocatableArray<Tsp> sp_first;
            F64 r_crit_sq = 1.0; // dummy
PS_OMP(omp for schedule(dynamic, 4))
            for(S32 i=0; i<n_proc; i++){
                const auto n_image_tmp_prev = shift_per_image_tmp[ith].size();
                rank_tmp[ith].push_back(i);
                CalcNumberAndShiftOfImageDomain
                    (shift_per_image_tmp[ith], dinfo.getPosRootDomain().getFullLength(),
                     outer_boundary_of_my_tree, dinfo.getPosDomain(i), pa, false);
                const auto n_image_tmp = shift_per_image_tmp[ith].size();
                n_image_per_proc[i] = n_image_tmp - n_image_tmp_prev;
                auto n_ep_prev = adr_ep_send_tmp[ith].size();
                for(S32 j=n_image_tmp_prev; j<n_image_tmp; j++){
                    auto n_ep_prev_2 = adr_ep_send_tmp[ith].size();
                    if(my_rank==i && j==n_image_tmp_prev){
                        n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2); // is 0
                        continue;
                    }
                    auto pos_target_domain = dinfo.getPosDomain(i).shift(shift_per_image_tmp[ith][j]);
                    TargetBox<SEARCH_MODE_SCATTER> target_box;
                    target_box.vertex_in_ = pos_target_domain;
                    MakeListUsingTreeRecursiveTop
                        <SEARCH_MODE_SCATTER, Ttc, TreeParticle, Tep, Tsp, Twalkmode, TagChopLeafFalse, TagCopyInfoCloseNoSp>
                        (tc_first, adr_tc, tp_first,
                         ep_first, adr_ep_send_tmp[ith],
                         sp_first, adr_sp_send_tmp,
                         target_box,
                         r_crit_sq, n_leaf_limit,
                         adr_tree_sp_first, F64vec(0.0));
                    n_ep_per_image_tmp[ith].push_back(adr_ep_send_tmp[ith].size() - n_ep_prev_2);
                }
                n_ep_send[i] = adr_ep_send_tmp[ith].size() - n_ep_prev;
            }
        } // end of OMP scope
        ReallocatableArray<S32> n_disp_image_per_proc(n_proc+1, n_proc+1, 1);
        n_disp_image_per_proc[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_image_per_proc[i+1] = n_disp_image_per_proc[i] + n_image_per_proc[i];
        }
        const auto n_image_tot = n_disp_image_per_proc[n_proc];
        shift_per_image.resizeNoInitialize( n_image_tot );
        n_ep_per_image.resizeNoInitialize( n_image_tot);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_thread; i++){
            S32 n_cnt = 0;
            for(S32 j=0; j<rank_tmp[i].size(); j++){
                auto rank = rank_tmp[i][j];
                auto offset = n_disp_image_per_proc[rank];
                for(S32 k=0; k<n_image_per_proc[rank]; k++){
                    shift_per_image[offset+k] = shift_per_image_tmp[i][n_cnt];
                    n_ep_per_image[offset+k] = n_ep_per_image_tmp[i][n_cnt];
                    n_cnt++;
                }
            }
        }
        ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot+1, n_image_tot+1, 1);
        n_disp_ep_per_image[0] = 0;
        for(S32 i=0; i<n_image_tot; i++){
            n_disp_ep_per_image[i+1] = n_disp_ep_per_image[i] + n_ep_per_image[i];
        }
        const auto n_ep_send_tot = n_disp_ep_per_image[ n_image_tot ];
        adr_ep_send.resizeNoInitialize( n_ep_send_tot );
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_thread; i++){
            auto n_cnt_ep = 0;
            for(S32 j=0; j<rank_tmp[i].size(); j++){
                auto rank = rank_tmp[i][j];
                const auto adr_image_head = n_disp_image_per_proc[rank];
                const auto adr_image_end = n_disp_image_per_proc[rank+1];
                for(S32 k=adr_image_head; k<adr_image_end; k++){
                    const auto adr_ep_head = n_disp_ep_per_image[k];
                    const auto adr_ep_end = n_disp_ep_per_image[k+1];
                    for(S32 l=adr_ep_head; l<adr_ep_end; l++){
                        adr_ep_send[l] = adr_ep_send_tmp[i][n_cnt_ep++];
                    }
                }
            }
        }
    }
#endif
    
    ///////////////////
    // exchange # of LET
    ////// FOR LONG SEARCH
    inline void ExchangeNumberLong(const ReallocatableArray<S32> & n_ep_send,
                                   ReallocatableArray<S32> & n_ep_recv,
                                   const ReallocatableArray<S32> & n_sp_send,
                                   ReallocatableArray<S32> & n_sp_recv,
                                   const ReallocatableArray<S32> & rank_send,
                                   const ReallocatableArray<S32> & rank_recv){
        ReallocatableArray<S32> n_ep_sp_send(rank_send.size()*2, rank_send.size()*2, 1);
        ReallocatableArray<S32> n_ep_sp_recv(rank_recv.size()*2, rank_recv.size()*2, 1);
PS_OMP_PARALLEL_FOR
        for(int i=0; i<rank_send.size(); i++){
            const S32 rank = rank_send[i];
            n_ep_sp_send[i*2]   = n_ep_send[rank];
            n_ep_sp_send[i*2+1] = n_sp_send[rank];
        }
        Comm::isendrecv(n_ep_sp_send.getPointer(), rank_send.getPointer(), 2, rank_send.size(),
                        n_ep_sp_recv.getPointer(), rank_recv.getPointer(), 2, rank_recv.size());
        const S32 n_proc = Comm::getNumberOfProc();
        n_ep_recv.resizeNoInitialize(n_proc);
        n_sp_recv.resizeNoInitialize(n_proc);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv[i] = n_sp_recv[i] = 0;
        }
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<rank_recv.size(); i++){
            const S32 rank = rank_recv[i];
            n_ep_recv[rank] = n_ep_sp_recv[i*2];
            n_sp_recv[rank] = n_ep_sp_recv[i*2+1];
        }
    }
    inline void ExchangeNumberLong(const ReallocatableArray<S32> & n_ep_send,
                                   ReallocatableArray<S32> & n_ep_recv,
                                   const ReallocatableArray<S32> & n_sp_send,
                                   ReallocatableArray<S32> & n_sp_recv){
        const S32 n_proc = Comm::getNumberOfProc();
        static ReallocatableArray<S32> n_ep_sp_send(n_proc*2);
        static ReallocatableArray<S32> n_ep_sp_recv(n_proc*2);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_proc; i++){
            n_ep_sp_send[i*2]   = n_ep_send[i];
            n_ep_sp_send[i*2+1] = n_sp_send[i];
        }
#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<S32, 2> comm_a2a_2d;
        comm_a2a_2d.execute(n_ep_sp_send, 2, n_ep_sp_recv);
#else
        Comm::allToAll(n_ep_sp_send.getPointer(), 2, n_ep_sp_recv.getPointer());
#endif
        n_ep_recv.resizeNoInitialize(n_proc);
        n_sp_recv.resizeNoInitialize(n_proc);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv[i] = n_ep_sp_recv[i*2];
            n_sp_recv[i] = n_ep_sp_recv[i*2+1];
        }
    }

    ///////////////////////
    ////// FOR SHORT SEARCH
    inline void ExchangeNumberShort
    (const ReallocatableArray<S32> & n_ep_send,
     ReallocatableArray<S32> & n_ep_recv,
     const ReallocatableArray<S32> & rank_send,
     const ReallocatableArray<S32> & rank_recv){
        ReallocatableArray<S32> n_ep_send_tmp(rank_send.size(), rank_send.size(), 1);
        ReallocatableArray<S32> n_ep_recv_tmp(rank_recv.size(), rank_recv.size(), 1);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<rank_send.size(); i++){
            const auto rank = rank_send[i];
            n_ep_send_tmp[i] = n_ep_send[rank];
        }
        Comm::isendrecv(n_ep_send_tmp.getPointer(), rank_send.getPointer(), 1, rank_send.size(),
                        n_ep_recv_tmp.getPointer(), rank_recv.getPointer(), 1, rank_recv.size());
        const auto n_proc = Comm::getNumberOfProc();
        n_ep_recv.resizeNoInitialize(n_proc);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv[i] = 0;
        }
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<rank_recv.size(); i++){
            const auto rank = rank_recv[i];
            n_ep_recv[rank] = n_ep_recv_tmp[i];
        }
    }
    
    inline void ExchangeNumberShort
    (const ReallocatableArray<S32> & n_ep_send,
     ReallocatableArray<S32> & n_ep_recv){
        const S32 n_proc = Comm::getNumberOfProc();
        n_ep_recv.resizeNoInitialize(n_proc);
#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<S32, 2> comm_a2a_2d;
        comm_a2a_2d.execute(n_ep_send, 1, n_ep_recv);
#else
        Comm::allToAll(n_ep_send.getPointer(), 1, n_ep_recv.getPointer());
#endif
    }
    // exchange # of LET
    //////////////

    //////////////////
    // EXCHANGE LET //
    //////////////////
    // FOR LONG SEARCH    
#ifdef USE_ALLGATHER_EXLET
    template<class TSM, class Tep, class Tsp, class Ttc>
    inline void ExchangeLet(const ReallocatableArray<Tep> & ep,
                            const ReallocatableArray<S32> & n_ep_send,
                            const ReallocatableArray<S32> & n_ep_recv,
                            const ReallocatableArray<S32> & n_ep_per_image,
                            const ReallocatableArray<S32> & adr_ep_send,
                            ReallocatableArray<Tep> & ep_org,
                            const S32 n_ep_offset, // ep_org[n_ep_offset] = n_ep_recv[0]
                            const ReallocatableArray<Ttc> & tc,
                            const ReallocatableArray<S32> & n_sp_send,
                            const ReallocatableArray<S32> & n_sp_recv,
                            const ReallocatableArray<S32> & n_sp_per_image,
                            const ReallocatableArray<S32> & adr_sp_send,
                            ReallocatableArray<Tsp> & sp_org,
                            const ReallocatableArray<F64vec> & shift_image_domain,
                            const ReallocatableArray<S32> & n_image_per_proc,
                            const ReallocatableArray<S32> & rank_send,
                            const ReallocatableArray<S32> & rank_recv){
        ReallocatableArray<S32> n_disp_ep_send(rank_send.size()+1, rank_send.size()+1, 1);
        ReallocatableArray<S32> n_disp_sp_send(rank_send.size()+1, rank_send.size()+1, 1);
        ReallocatableArray<S32> n_disp_ep_recv(rank_recv.size()+1, rank_recv.size()+1, 1);
        ReallocatableArray<S32> n_disp_sp_recv(rank_recv.size()+1, rank_recv.size()+1, 1);
        ReallocatableArray<S32> n_ep_send_tmp(rank_send.size(), rank_send.size(), 1);
        ReallocatableArray<S32> n_sp_send_tmp(rank_send.size(), rank_send.size(), 1);
        ReallocatableArray<S32> n_ep_recv_tmp(rank_recv.size(), rank_recv.size(), 1);
        ReallocatableArray<S32> n_sp_recv_tmp(rank_recv.size(), rank_recv.size(), 1);
        
        n_disp_ep_send[0] = n_disp_sp_send[0] = n_disp_ep_recv[0] = n_disp_sp_recv[0] = 0;
        for(int i=0; i<rank_send.size(); i++){
            const S32 rank = rank_send[i];
            n_disp_ep_send[i+1] = n_ep_send[rank] + n_disp_ep_send[i];
            n_disp_sp_send[i+1] = n_sp_send[rank] + n_disp_sp_send[i];
            n_ep_send_tmp[i] = n_ep_send[rank];
            n_sp_send_tmp[i] = n_sp_send[rank];
        }
        for(int i=0; i<rank_recv.size(); i++){
            const S32 rank = rank_recv[i];
            n_disp_ep_recv[i+1] = n_ep_recv[rank] + n_disp_ep_recv[i];
            n_disp_sp_recv[i+1] = n_sp_recv[rank] + n_disp_sp_recv[i];
            n_ep_recv_tmp[i] = n_ep_recv[rank];
            n_sp_recv_tmp[i] = n_sp_recv[rank];
        }
        const S32 n_image_tot = n_ep_per_image.size();
        ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot+1, n_image_tot+1, 1);
        ReallocatableArray<S32> n_disp_sp_per_image(n_image_tot+1, n_image_tot+1, 1);
        n_disp_ep_per_image[0] = n_disp_sp_per_image[0] = 0;
        for(S32 i=0; i<n_image_tot; i++){
            n_disp_ep_per_image[i+1] = n_disp_ep_per_image[i] +  n_ep_per_image[i];
            n_disp_sp_per_image[i+1] = n_disp_sp_per_image[i] +  n_sp_per_image[i];
        }
        ReallocatableArray<Tep> ep_send(n_disp_ep_send[rank_send.size()], n_disp_ep_send[rank_send.size()], 1);
        ReallocatableArray<Tsp> sp_send(n_disp_sp_send[rank_send.size()], n_disp_sp_send[rank_send.size()], 1);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4)
#endif
        for(S32 i=0; i<n_image_tot; i++){
            F64vec shift = shift_image_domain[i];
            S32 n_ep = n_ep_per_image[i];
            S32 n_ep_offset = n_disp_ep_per_image[i];
            CopyPtclToSendBuf(typename TSM::search_boundary_type(), ep_send, ep, adr_ep_send, n_ep, n_ep_offset, shift);
            S32 n_sp = n_sp_per_image[i];
            S32 n_sp_offset = n_disp_sp_per_image[i];
            CopyPtclFromTreeToSendBuf(typename TSM::search_boundary_type(), sp_send, tc, adr_sp_send, n_sp, n_sp_offset, shift);
        }
        ep_org.resizeNoInitialize(n_disp_ep_recv[rank_recv.size()] + n_ep_offset);
        sp_org.resizeNoInitialize(n_disp_sp_recv[rank_recv.size()] );
        Comm::isendrecvV(ep_send.getPointer(), rank_send.getPointer(), n_ep_send_tmp.getPointer(), n_disp_ep_send.getPointer(), rank_send.size(),
                         sp_send.getPointer(), rank_send.getPointer(), n_sp_send_tmp.getPointer(), n_disp_sp_send.getPointer(), rank_send.size(), 
                         ep_org.getPointer(n_ep_offset), rank_recv.getPointer(), n_ep_recv_tmp.getPointer(), n_disp_ep_recv.getPointer(), rank_recv.size(),
                         sp_org.getPointer(),  rank_recv.getPointer(), n_sp_recv_tmp.getPointer(), n_disp_sp_recv.getPointer(), rank_recv.size());
    }
#else
    template<class TSM, class Tep, class Tsp, class Ttc>
    inline void ExchangeLet
    (const ReallocatableArray<Tep> & ep,
     const ReallocatableArray<S32> & n_ep_send,
     const ReallocatableArray<S32> & n_ep_recv,
     const ReallocatableArray<S32> & n_ep_per_image,
     const ReallocatableArray<S32> & adr_ep_send,
     ReallocatableArray<Tep> & ep_org,
     const S32 n_ep_offset,
     const ReallocatableArray<Ttc> & tc,
     const ReallocatableArray<S32> & n_sp_send,
     const ReallocatableArray<S32> & n_sp_recv,
     const ReallocatableArray<S32> & n_sp_per_image,
     const ReallocatableArray<S32> & adr_sp_send,
     ReallocatableArray<Tsp> & sp_org,
     const ReallocatableArray<F64vec> & shift_image_domain,
     const ReallocatableArray<S32> & n_image_per_proc){
        const S32 n_proc = Comm::getNumberOfProc();
        //const S32 my_rank = Comm::getRank();
        ReallocatableArray<S32> n_disp_ep_send(n_proc+1, n_proc+1, 1);
        ReallocatableArray<S32> n_disp_sp_send(n_proc+1, n_proc+1, 1);
        ReallocatableArray<S32> n_disp_ep_recv(n_proc+1, n_proc+1, 1);
        ReallocatableArray<S32> n_disp_sp_recv(n_proc+1, n_proc+1, 1);
        n_disp_ep_send[0] = n_disp_sp_send[0] = n_disp_ep_recv[0] = n_disp_sp_recv[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_ep_send[i+1] = n_ep_send[i] + n_disp_ep_send[i];
            n_disp_sp_send[i+1] = n_sp_send[i] + n_disp_sp_send[i];
            n_disp_ep_recv[i+1] = n_ep_recv[i] + n_disp_ep_recv[i];
            n_disp_sp_recv[i+1] = n_sp_recv[i] + n_disp_sp_recv[i];
        }
        const S32 n_image_tot = n_ep_per_image.size();
        //if(my_rank==0) std::cerr<<"n_image_tot= "<<n_image_tot<<std::endl;
        ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot+1, n_image_tot+1, 1);
        ReallocatableArray<S32> n_disp_sp_per_image(n_image_tot+1, n_image_tot+1, 1);
        n_disp_ep_per_image[0] = n_disp_sp_per_image[0] = 0;
        for(S32 i=0; i<n_image_tot; i++){
            n_disp_ep_per_image[i+1] = n_disp_ep_per_image[i] +  n_ep_per_image[i];
            n_disp_sp_per_image[i+1] = n_disp_sp_per_image[i] +  n_sp_per_image[i];
        }
        ReallocatableArray<Tep> ep_send( n_disp_ep_send[n_proc], n_disp_ep_send[n_proc], 1 );
        ReallocatableArray<Tsp> sp_send( n_disp_sp_send[n_proc], n_disp_sp_send[n_proc], 1 );

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4)
#endif
        for(S32 i=0; i<n_image_tot; i++){
            //F64 cm_mass = 0.0;
            //F64vec cm_pos = 0.0;
            F64vec shift = shift_image_domain[i];

            S32 n_ep = n_ep_per_image[i];
            S32 n_ep_offset = n_disp_ep_per_image[i];
            CopyPtclToSendBuf(typename TSM::search_boundary_type(), ep_send, ep, adr_ep_send, n_ep,
                              n_ep_offset, shift);
            //for(S32 i2=n_ep_offset; i2<n_ep_offset+n_ep; i2++){
            //    cm_mass += ep_send[i2].mass;
            //    cm_pos += ep_send[i2].mass * ep_send[i2].pos;
            //}
            
            S32 n_sp = n_sp_per_image[i];
            S32 n_sp_offset = n_disp_sp_per_image[i];
            CopyPtclFromTreeToSendBuf(typename TSM::search_boundary_type(), sp_send, tc, adr_sp_send, n_sp, n_sp_offset, shift);
            /*
            for(S32 i2=n_sp_offset; i2<n_sp_offset+n_sp; i2++){
                cm_mass += sp_send[i2].mass;
                cm_pos += sp_send[i2].mass * sp_send[i2].pos;
            }
            if(Comm::getRank()==0){
                std::cerr<<"i= "<<i
                         <<" cm_mass= "<<cm_mass
                         <<" cm_pos= "<<cm_pos / cm_mass
                         <<std::endl;
            }
            */

        }
        //exit(1);
        ep_org.resizeNoInitialize( n_disp_ep_recv[n_proc] + n_ep_offset);
        sp_org.resizeNoInitialize( n_disp_sp_recv[n_proc] );

#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<Tep, 2> comm_a2a_epj_2d;
        comm_a2a_epj_2d.executeV(ep_send, ep_org, n_ep_send.getPointer(), n_ep_recv.getPointer(), 0, n_ep_offset);
        static CommForAllToAll<Tsp, 2> comm_a2a_spj_2d;
        comm_a2a_spj_2d.executeV(sp_send, sp_org, n_sp_send.getPointer(), n_sp_recv.getPointer(), 0, 0);
#else
        Comm::allToAllV(ep_send.getPointer(), n_ep_send.getPointer(), n_disp_ep_send.getPointer(),
                        ep_org.getPointer(n_ep_offset), n_ep_recv.getPointer(), n_disp_ep_recv.getPointer());
        Comm::allToAllV(sp_send.getPointer(), n_sp_send.getPointer(), n_disp_sp_send.getPointer(),
                        sp_org.getPointer(), n_sp_recv.getPointer(), n_disp_sp_recv.getPointer());
#endif
    }
#endif

    //////////////////
    // FOR SHORT SEARCH
#ifdef USE_ALLGATHER_EXLET
    template<typename TSM, typename Tep>
    inline void ExchangeLet
    (const ReallocatableArray<Tep> & ep,
     const ReallocatableArray<S32> & n_ep_send,
     const ReallocatableArray<S32> & n_ep_recv,
     const ReallocatableArray<S32> & n_ep_per_image,
     const ReallocatableArray<S32> & adr_ep_send,
     ReallocatableArray<Tep> & ep_org,
     const S32 n_ep_offset,
     const ReallocatableArray<F64vec> & shift_image_domain,
     const ReallocatableArray<S32> & n_image_per_proc,
     const ReallocatableArray<S32> & rank_send,
     const ReallocatableArray<S32> & rank_recv){
        ReallocatableArray<S32> n_ep_send_tmp(rank_send.size(), rank_send.size(), 1);
        ReallocatableArray<S32> n_ep_recv_tmp(rank_recv.size(), rank_recv.size(), 1);
        ReallocatableArray<S32> n_disp_ep_send(rank_send.size()+1, rank_send.size()+1, 1);
        ReallocatableArray<S32> n_disp_ep_recv(rank_recv.size()+1, rank_recv.size()+1, 1);
        n_disp_ep_send[0] = n_disp_ep_recv[0] = 0;
        for(int i=0; i<rank_send.size(); i++){
            const S32 rank = rank_send[i];
            n_disp_ep_send[i+1] = n_ep_send[rank] + n_disp_ep_send[i];
            n_ep_send_tmp[i] = n_ep_send[rank];
        }
        for(int i=0; i<rank_recv.size(); i++){
            const S32 rank = rank_recv[i];
            n_disp_ep_recv[i+1] = n_ep_recv[rank] + n_disp_ep_recv[i];
            n_ep_recv_tmp[i] = n_ep_recv[rank];
        }
        const S32 n_image_tot = n_ep_per_image.size();
        ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot+1, n_image_tot+1, 1);
        n_disp_ep_per_image[0] = 0;
        for(S32 i=0; i<n_image_tot; i++){
            n_disp_ep_per_image[i+1] = n_disp_ep_per_image[i] +  n_ep_per_image[i];
        }
        ReallocatableArray<Tep> ep_send(n_disp_ep_send[rank_send.size()], n_disp_ep_send[rank_send.size()], 1);
PS_OMP(omp parallel for schedule(dynamic, 4))
        for(S32 i=0; i<n_image_tot; i++){
            F64vec shift = shift_image_domain[i];
            S32 n_ep = n_ep_per_image[i];
            S32 n_ep_offset = n_disp_ep_per_image[i];
            CopyPtclToSendBuf(typename TSM::search_boundary_type(), ep_send, ep, adr_ep_send, n_ep, n_ep_offset, shift);
        }
        ep_org.resizeNoInitialize(n_disp_ep_recv[rank_recv.size()] + n_ep_offset);
        Comm::isendrecvV(ep_send.getPointer(), rank_send.getPointer(), n_ep_send_tmp.getPointer(), n_disp_ep_send.getPointer(), rank_send.size(),
                         ep_org.getPointer(n_ep_offset), rank_recv.getPointer(), n_ep_recv_tmp.getPointer(), n_disp_ep_recv.getPointer(), rank_recv.size());
    }
#else    
    template<typename TSM, typename Tep>
    inline void ExchangeLet
    (const ReallocatableArray<Tep> & ep,
     const ReallocatableArray<S32> & n_ep_send,
     const ReallocatableArray<S32> & n_ep_recv,
     const ReallocatableArray<S32> & n_ep_per_image,
     const ReallocatableArray<S32> & adr_ep_send,
     ReallocatableArray<Tep> & ep_org,
     const S32 n_ep_offset, 
     const ReallocatableArray<F64vec> & shift_image_domain,
     const ReallocatableArray<S32> & n_image_per_proc){
        const S32 n_proc = Comm::getNumberOfProc();
        ReallocatableArray<S32> n_disp_ep_send(n_proc+1, n_proc+1, 1);
        ReallocatableArray<S32> n_disp_ep_recv(n_proc+1, n_proc+1, 1);
        n_disp_ep_send[0] = n_disp_ep_recv[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_ep_send[i+1] = n_ep_send[i] + n_disp_ep_send[i];
            n_disp_ep_recv[i+1] = n_ep_recv[i] + n_disp_ep_recv[i];
        }
        const S32 n_image_tot = n_ep_per_image.size();
        ReallocatableArray<S32> n_disp_ep_per_image(n_image_tot+1, n_image_tot+1, 1);
        n_disp_ep_per_image[0] = 0;
        for(S32 i=0; i<n_image_tot; i++){
            n_disp_ep_per_image[i+1] = n_disp_ep_per_image[i] +  n_ep_per_image[i];
        }
        ReallocatableArray<Tep> ep_send( n_disp_ep_send[n_proc], n_disp_ep_send[n_proc], 1 );
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4)
#endif
        for(S32 i=0; i<n_image_tot; i++){
            F64vec shift = shift_image_domain[i];
            S32 n_ep = n_ep_per_image[i];
            S32 n_ep_offset = n_disp_ep_per_image[i];
            CopyPtclToSendBuf(typename TSM::search_boundary_type(), ep_send, ep, adr_ep_send, n_ep, n_ep_offset, shift);
        }
        ep_org.resizeNoInitialize( n_disp_ep_recv[n_proc] + n_ep_offset);
#ifdef FAST_ALL_TO_ALL_FOR_K
        static CommForAllToAll<Tep, 2> comm_a2a_epj_2d;
        comm_a2a_epj_2d.executeV(ep_send, ep_org, n_ep_send.getPointer(), n_ep_recv.getPointer(), 0, n_ep_offset);
#else
        Comm::allToAllV(ep_send.getPointer(), n_ep_send.getPointer(), n_disp_ep_send.getPointer(),
                        ep_org.getPointer(n_ep_offset), n_ep_recv.getPointer(), n_disp_ep_recv.getPointer());
#endif
    }
#endif
    // exchange LET
    //////////////

    template<class Ttc, class Ttp, class Tepj>
    inline void FindExchangeParticleDoubleWalk(const ReallocatableArray<Tepj> & epj_A, // received particles
                                               const ReallocatableArray<Ttc> & tc_first_B,
                                               const ReallocatableArray<S32> & n_epj_src_per_proc,
                                               const ReallocatableArray<S32> & n_image_send_per_proc_irnai, // not needed
                                               const DomainInfo & dinfo,
                                               const S32 n_leaf_limit_B,
                                               ReallocatableArray<S32> & n_epj_send_per_proc,
                                               ReallocatableArray<S32> & n_epj_send_per_image,
                                               ReallocatableArray<S32> & n_image_send_per_proc,
                                               ReallocatableArray<S32> & adr_ep_send,
                                               ReallocatableArray<F64vec> & shift_per_image,
                                               const ReallocatableArray<Tepj> & epj_B, // assigned
                                               const F64vec & center_tree,
                                               const F64    & full_len_tree
                                               ){
        const S32 n_proc = Comm::getNumberOfProc();
        //const S32 my_rank = Comm::getRank();
        /*
        if(my_rank==0){
            std::cerr<<"epj_A.size()= "<<epj_A.size()
                     <<" tc_first_B.size()= "<<tc_first_B.size()
                     <<std::endl;
            for(S32 i=0; i<n_proc; i++){
                std::cerr<<"n_epj_src_per_proc[i]= "<<n_epj_src_per_proc[i]<<std::endl;
            }
        }
        */
        const S32 n_thread = Comm::getNumberOfThread();
        const F64ort pos_root_domain = dinfo.getPosRootDomain();
        const F64vec len_root_domain = pos_root_domain.getFullLength();
        n_image_send_per_proc.resizeNoInitialize(n_proc);
        n_epj_send_per_proc.resizeNoInitialize(n_proc);
        ReallocatableArray<S32> n_disp_epj_per_proc;
        n_disp_epj_per_proc.resizeNoInitialize(n_proc+1);
        n_disp_epj_per_proc[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_epj_per_proc[i+1] = n_disp_epj_per_proc[i] + n_epj_src_per_proc[i];
        }
        for(S32 i=0; i<n_proc; i++){
            n_epj_send_per_proc[i] = 0;
        }
        
        static ReallocatableArray<Tepj> * epj_sorted_tmp;
        static ReallocatableArray<F64vec> * shift_per_image_tmp;
        static ReallocatableArray<Ttp> * tp_tmp;
        static ReallocatableArray<Ttc> * tc_tmp;
        static ReallocatableArray<S32> * adr_tc_level_partition_tmp;
        static ReallocatableArray<S32> * adr_ptcl_send_tmp;
        static ReallocatableArray<S32> * rank_dst_tmp;
        static ReallocatableArray<S32> * adr_epj_src_per_image_tmp;
        static ReallocatableArray<S32> * n_epj_src_per_image_tmp;
        static ReallocatableArray<S32> * n_epj_send_per_image_tmp;
        static bool first = true;
        if(first){
            epj_sorted_tmp = new ReallocatableArray<Tepj>[n_thread];
            shift_per_image_tmp = new ReallocatableArray<F64vec>[n_thread];
            tp_tmp = new ReallocatableArray<Ttp>[n_thread];
            tc_tmp = new ReallocatableArray<Ttc>[n_thread];
            adr_tc_level_partition_tmp = new ReallocatableArray<S32>[n_thread];
            adr_ptcl_send_tmp = new ReallocatableArray<S32>[n_thread];
            for(S32 i=0; i<n_thread; i++){
                adr_tc_level_partition_tmp[i].resizeNoInitialize(TREE_LEVEL_LIMIT+2);
            }
            rank_dst_tmp = new ReallocatableArray<S32>[n_thread];
            adr_epj_src_per_image_tmp = new ReallocatableArray<S32>[n_thread];
            n_epj_src_per_image_tmp = new ReallocatableArray<S32>[n_thread];
            n_epj_send_per_image_tmp = new ReallocatableArray<S32>[n_thread];
            first = false;
        }
        ReallocatableArray<S32> rank_src(n_proc, 0, 1);
        for(S32 i=0; i<n_proc; i++){
            n_image_send_per_proc[i] = 0;
            if(n_epj_src_per_proc[i] > 0) rank_src.push_back(i);
        }
        ReallocatableArray<F64> len_peri(DIMENSION, DIMENSION, 1);
        bool pa[DIMENSION];
        dinfo.getPeriodicAxis(pa);
        for(S32 i=0; i<DIMENSION; i++){
            if(pa[i]==false) len_peri[i] = 0.0;
        }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
            const S32 ith = Comm::getThreadNum();
            epj_sorted_tmp[ith].clearSize();
            shift_per_image_tmp[ith].clearSize();
            adr_ptcl_send_tmp[ith].clearSize();
            rank_dst_tmp[ith].clearSize();
            tc_tmp[ith].clearSize();
            tp_tmp[ith].clearSize();
            adr_epj_src_per_image_tmp[ith].clearSize();
            n_epj_src_per_image_tmp[ith].clearSize();
            n_epj_send_per_image_tmp[ith].clearSize();
            //S32 n_ep_send_cum_old = 0;
            //bool first_loop = true;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 1)
            //#pragma omp for schedule(static)
#endif
            for(S32 ib=0; ib<rank_src.size(); ib++){
                const S32 rank_tmp = rank_src[ib];
                rank_dst_tmp[ith].push_back(rank_tmp);
                /*
                if(Comm::getRank() == 0){
                    std::cerr<<"rank_tmp= "<<rank_tmp
                             <<" n_epj_per_proc[rank_tmp]= "
                             <<n_epj_per_proc[rank_tmp]
                             <<std::endl;
                }
                */
                if( n_epj_src_per_proc[rank_tmp] <= 0) continue;
                const S32 adr_ptcl_head = n_disp_epj_per_proc[rank_tmp];
                const S32 adr_ptcl_end = n_disp_epj_per_proc[rank_tmp+1];
                S32vec id_image_old = -9999;
                S32 n_image = 0;
                for(S32 ip=adr_ptcl_head; ip<adr_ptcl_end; ip++){
                    const F64vec pos_target = epj_A[ip].getPos();
                    const S32vec id_image_new = CalcIDOfImageDomain(pos_root_domain, pos_target, pa);
                    if(id_image_old != id_image_new){
                        adr_epj_src_per_image_tmp[ith].push_back(ip);
                        n_epj_src_per_image_tmp[ith].push_back(0);
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
                        shift_per_image_tmp[ith].push_back(F64vec(id_image_new.x*len_root_domain.x, id_image_new.y*len_root_domain.y));
#else
                        shift_per_image_tmp[ith].push_back(F64vec(id_image_new.x*len_root_domain.x, id_image_new.y*len_root_domain.y, id_image_new.z*len_root_domain.z));
#endif
                        id_image_old = id_image_new;
                        n_image++;
                    }
                    n_epj_src_per_image_tmp[ith].back()++;
                }
                /*
                if(Comm::getRank()==0 && rank_tmp==1){
                    std::cerr<<"adr_ptcl_head= "<<adr_ptcl_head
                             <<" adr_ptcl_end= "<<adr_ptcl_end
                             <<std::endl;
                }
                */
                /*
                if(Comm::getRank()==0){
                    std::cerr<<"rank_tmp= "<<rank_tmp
                             <<" n_image= "<<n_image
                             <<" n_disp_epj_per_image_tmp[ith].size()= "<<n_disp_epj_per_image_tmp[ith].size()
                             <<std::endl;
                    for(S32 ip=adr_ptcl_head; ip<adr_ptcl_end; ip++){
                        std::cerr<<"ip, pos= "<<ip<<", "<<epj_A[ip].getPos()<<std::endl;
                    }
                }
                */
                n_image_send_per_proc[rank_tmp] = n_image;
                const S32 adr_image_end  = n_epj_src_per_image_tmp[ith].size();
                const S32 adr_image_head = adr_image_end - n_image;
                const S32 n_epj_send_cum_prev = adr_ptcl_send_tmp[ith].size();
                for(S32 ii=adr_image_head; ii<adr_image_end; ii++){
                    const F64vec shift = shift_per_image_tmp[ith][ii];
                    /*
                    if(Comm::getRank()==0){
                        std::cerr<<"rank_tmp= "<<rank_tmp
                                 <<" ii= "<<ii
                                 <<" shift= "<<shift
                                 <<" n_epj_src_per_image_tmp[ith][ii]= "<<n_epj_src_per_image_tmp[ith][ii]
                                 <<std::endl;
                    }
                    */
                    const F64ort pos_domain = dinfo.getPosDomain(rank_tmp).shift(shift);
                    //const F64ort pos_domain = dinfo.getPosDomain(rank_tmp).shift(-shift);
                    /*
                    if(Comm::getRank() == 0 && rank_tmp==1){
                        std::cerr<<"rank_tmp= "<<rank_tmp
                                 <<" shift= "<<shift
                                 <<" pos_domain= "<<pos_domain
                                 <<" n_disp_epj_per_image_tmp[ith][ii]= "
                                 <<n_disp_epj_per_image_tmp[ith][ii]
                                 <<" n_disp_epj_per_image_tmp[ith][ii+1]= "
                                 <<n_disp_epj_per_image_tmp[ith][ii+1]
                                 <<std::endl;
                    }
                    */
                    ///////////////
                    // MAKE TREE A
                    S32 n_cnt = 0;
                    tp_tmp[ith].resizeNoInitialize(n_epj_src_per_image_tmp[ith][ii]);
                    for(S32 ip=adr_epj_src_per_image_tmp[ith][ii]; ip<adr_epj_src_per_image_tmp[ith][ii]+n_epj_src_per_image_tmp[ith][ii]; ip++, n_cnt++){
                        tp_tmp[ith][n_cnt].setFromEP(epj_A[ip], ip);
                        //if(Comm::getRank()==0) std::cerr<<"epj_A[ip].pos= "<<epj_A[ip].pos<<std::endl;
                    }
                    std::sort(tp_tmp[ith].getPointer(), tp_tmp[ith].getPointer()+n_cnt, LessOPKEY());
                    epj_sorted_tmp[ith].resizeNoInitialize(n_cnt);
                    for(S32 ip=0; ip<n_cnt; ip++){
                        const S32 adr = tp_tmp[ith][ip].adr_ptcl_;
                        epj_sorted_tmp[ith][ip] = epj_A[adr];
                        tp_tmp[ith][ip].adr_ptcl_ = ip;
                        /*
                        if(Comm::getRank()==0 && rank_tmp == 1){
                            std::cout<<"ip= "<<ip<<" epj_A[adr].pos= "<<epj_A[adr].pos<<std::endl;
                        }
                        */
                    }
                    /*
                    if(Comm::getRank() == 0){
                        for(S32 ip=0; ip<n_cnt; ip++){
                            std::cerr<<"epj_sorted_tmp[ith][ip].pos= "<<epj_sorted_tmp[ith][ip].pos<<std::endl;
                        }
                    }
                    */
                    const S32 n_leaf_limit_A = 1;
                    S32 lev_max_A = 0;
                    LinkCellST(tc_tmp[ith], adr_tc_level_partition_tmp[ith].getPointer(),
                               tp_tmp[ith].getPointer(), lev_max_A, n_cnt, n_leaf_limit_A);
                    CalcMomentST(adr_tc_level_partition_tmp[ith].getPointer(),
                                 tc_tmp[ith].getPointer(), 
                                 epj_sorted_tmp[ith].getPointer(), lev_max_A, n_leaf_limit_A);
                    /*
                    if(Comm::getRank() == 0 && rank_tmp == 1){
                        std::cerr<<"tc_tmp[ith][0].mom_.getVertexOut()= "<<tc_tmp[ith][0].mom_.getVertexOut()
                                 <<"tc_tmp[ith][0].mom_.getVertexIn()= "<<tc_tmp[ith][0].mom_.getVertexIn()
                                 <<" epj_sorted_tmp[ith][0].pos= "<<epj_sorted_tmp[ith][0].pos
                                 <<std::endl;
                    }
                    */
                    /*
                    if(Comm::getRank() == 0){
                        S32 err = 0;
                        tc_tmp[ith].getPointer()->checkTree(epj_sorted_tmp[ith].getPointer(),
                                                            tc_tmp[ith].getPointer(),
                                                            center_tree, full_len_tree*0.5,
                                                            n_leaf_limit_A, 1e-4,
                                                            err);
                        tc_tmp[ith].getPointer()->dumpTree(epj_sorted_tmp[ith].getPointer(),
                                                           tc_tmp[ith].getPointer(),
                                                           center_tree, full_len_tree*0.5,
                                                           n_leaf_limit_A);
                    }
                    */
                    const S32 n_epj_send_per_image_prev = adr_ptcl_send_tmp[ith].size();
                    MakeListDoubleWalkTop(tc_tmp[ith], tc_first_B,
                                          epj_A,    pos_domain,
                                          n_leaf_limit_A, n_leaf_limit_B,
                                          adr_ptcl_send_tmp[ith]);
                    n_epj_send_per_image_tmp[ith].push_back(adr_ptcl_send_tmp[ith].size() - n_epj_send_per_image_prev);

                    if(Comm::getRank() == 0 && rank_tmp == 1){
                        /*
                        std::cout<<"rank_tmp= "<<rank_tmp
                                 <<" shift= "<<shift
                                 <<" pos_domain= "<<pos_domain
                                 <<" tc_tmp[ith].size()= "<<tc_tmp[ith].size()
                                 <<" tc_first_B.size()= "<<tc_first_B.size()
                                 <<" n_epj_send_per_image_tmp[ith].back()= "<<n_epj_send_per_image_tmp[ith].back()
                                 <<" adr_ptcl_send_tmp[ith].size()= "<<adr_ptcl_send_tmp[ith].size()
                                 <<std::endl;
                        */
                        /*
                        for(S32 iii=n_epj_send_per_image_prev; iii<adr_ptcl_send_tmp[ith].size(); iii++){
                            std::cerr<<"adr_ptcl_send_tmp[ith][iii]= "<<adr_ptcl_send_tmp[ith][iii]
                                     <<" epj_B[adr_ptcl_send_tmp[ith][iii]].id= "<<epj_B[adr_ptcl_send_tmp[ith][iii]].id
                                     <<" pos= "<<epj_B[adr_ptcl_send_tmp[ith][iii]].pos
                                     <<std::endl;
                        }
                        */
                    }

                    /*
                    const S32 n_ep_per_image = adr_ptcl_send_tmp[ith].size();
                    n_ep_per_image_tmp[ith].push_back(n_ep_per_image);
                    for(S32 jp=0; jp<n_ep_per_image; jp++){
                        const S32 adr = adr_ptcl_send_tmp[ith][jp];
                        //const F64vec pos_j = epj[adr].getPos();
                        //if( pos_root_cell_.notOverlapped(pos_j-shift) ) continue;
                    }
                    */
                } // end of for ii=0 to n_image
                n_epj_send_per_proc[rank_tmp] = adr_ptcl_send_tmp[ith].size() - n_epj_send_cum_prev;
            } // end of OMP for
        } // end of OMP scope
        /*
        if(Comm::getRank()==0){
            for(S32 i=0; i<n_proc; i++){
                std::cerr<<"n_image_send_per_proc[i]= "<<n_image_send_per_proc[i]
                         <<" n_epj_send_per_proc[i]= "<<n_epj_send_per_proc[i]
                         <<std::endl;
            }
        }
        */
        /*
        if(Comm::getRank()==0){
            for(S32 i=0; i<n_thread; i++){
                std::cerr<<"rank_send_tmp[i].size()= "<<rank_send_tmp[i].size()<<std::endl;
                for(S32 j=0; j<rank_send_tmp[i].size(); j++){
                    std::cerr<<"rank_send_tmp[i][j]= "<<rank_send_tmp[i][j]<<std::endl;
                }
            }
        }
        */

        ReallocatableArray<S32> n_disp_image_send_per_proc;
        n_disp_image_send_per_proc.resizeNoInitialize(n_proc+1);
        n_disp_image_send_per_proc[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_image_send_per_proc[i+1] = n_disp_image_send_per_proc[i] + n_image_send_per_proc[i];
        }
        /*
        if(Comm::getRank()==0){
            for(S32 i=0; i<n_proc; i++){
                std::cerr<<"n_image_send_per_proc[i]= "<<n_image_send_per_proc[i]
                         <<" n_disp_image_send_per_proc[i]= "<<n_disp_image_send_per_proc[i]
                         <<std::endl;
            }
        }
        */
        
        /*
        S32 n_image_send_tot = 0;
        for(S32 i=0; i<n_thread; i++){
            n_image_send_tot += shift_per_image_tmp[i].size();
        }
        if(Comm::getRank()==0){
            std::cerr<<"n_image_send_tot= "<<n_image_send_tot<<std::endl;
        }
        */

        const S32 n_image_send_tot = n_disp_image_send_per_proc[n_proc];
        n_epj_send_per_image.resizeNoInitialize(n_image_send_tot);
        shift_per_image.resizeNoInitialize(n_image_send_tot);
        for(S32 i=0; i<n_thread; i++){
            S32 n_cnt_image = 0;
            //S32 n_cnt_ep = 0;
            for(S32 j=0; j<rank_dst_tmp[i].size(); j++){
                const S32 rank = rank_dst_tmp[i][j];
                const S32 adr_image_head = n_disp_image_send_per_proc[rank];
                const S32 adr_image_end  = n_disp_image_send_per_proc[rank+1];
                /*
                if(Comm::getRank()==0){
                    std::cerr<<"rank= "<<rank
                             <<" adr_image_head= "<<adr_image_head
                             <<" adr_image_end= "<<adr_image_end
                             <<std::endl;
                }
                */
                for(S32 k=adr_image_head; k<adr_image_end; k++, n_cnt_image++){
                    n_epj_send_per_image[k] = n_epj_send_per_image_tmp[i][n_cnt_image];
                    shift_per_image[k] = shift_per_image_tmp[i][n_cnt_image];
                    /*
                    if(Comm::getRank()==0){
                        std::cerr<<"k= "<<k
                                 <<" n_epj_send_per_image[k]= "<<n_epj_send_per_image[k]
                                 <<std::endl;
                    }
                    */
                }
            }
        }
        ReallocatableArray<S32> n_disp_epj_send_per_image;
        n_disp_epj_send_per_image.resizeNoInitialize(n_image_send_tot+1);
        n_disp_epj_send_per_image[0] = 0;
        for(S32 i=0; i<n_image_send_tot; i++){
            n_disp_epj_send_per_image[i+1] = n_disp_epj_send_per_image[i] + n_epj_send_per_image[i];
        }
        S32 n_epj_send_tot = 0;
        for(S32 i=0; i<n_proc; i++) n_epj_send_tot += n_epj_send_per_proc[i];
        adr_ep_send.resizeNoInitialize(n_epj_send_tot);
        
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
        for(S32 i=0; i<n_thread; i++){
            S32 n_cnt_image = 0;
            S32 n_cnt_ep = 0;
            for(S32 j=0; j<rank_dst_tmp[i].size(); j++){
                const S32 rank = rank_dst_tmp[i][j];
                const S32 adr_image_head = n_disp_image_send_per_proc[rank];
                const S32 adr_image_end  = n_disp_image_send_per_proc[rank+1];
                for(S32 k=adr_image_head; k<adr_image_end; k++, n_cnt_image++){
                    n_epj_send_per_image[k] = n_epj_send_per_image_tmp[i][n_cnt_image];
                    shift_per_image[k] = shift_per_image_tmp[i][n_cnt_image];
                    const S32 adr_epj_head = n_disp_epj_send_per_image[k];
                    const S32 adr_epj_end  = n_disp_epj_send_per_image[k+1];
                    for(S32 l=adr_epj_head; l<adr_epj_end; l++, n_cnt_ep++){
                        adr_ep_send[l] = adr_ptcl_send_tmp[i][n_cnt_ep];
                        /*
                        if(Comm::getRank()==0){
                            std::cerr<<"l= "<<l
                                     <<" adr_ep_send[l]= "<<adr_ep_send[l]
                                     <<std::endl;
                        }
                        */
                    }
                }
            }
        }
        /*
        if(Comm::getRank()==0){
            for(S32 i=0; i<n_proc; i++){
                std::cout<<"rank= "<<i<<std::endl;
                const S32 adr_image_head = n_disp_image_send_per_proc[i];
                const S32 adr_image_end  = n_disp_image_send_per_proc[i+1];
                for(S32 j=adr_image_head; j<adr_image_end; j++){
                    std::cout<<"image= "<<j
                             <<" shift_per_image[j]= "<<shift_per_image[j]
                             <<std::endl;
                    const S32 adr_epj_head = n_disp_epj_send_per_image[j];
                    const S32 adr_epj_end  = n_disp_epj_send_per_image[j+1];
                    for(S32 k=adr_epj_head; k<adr_epj_end; k++){
                        std::cout<<"k= "<<k
                                 <<" adr_ep_send[k]= "<<adr_ep_send[k]
                                 <<std::endl;
                    }
                }
            }
        }
        */
        //Comm::barrier();
        //exit(1);        
    }
    
    //////////////////
    // add moment as sp
    template<class Ttreecell, class Tspj>
    inline void AddMomentAsSpImpl(TagForceLong,
                                  const ReallocatableArray<Ttreecell> & _tc,
                                  const S32 offset,
                                  ReallocatableArray<Tspj> & _spj){
        //if(Comm::getRank()==0){
        //    std::cerr<<"offset= "<<offset
        //             <<" _tc.size()= "<<_tc.size()
        //             <<" _spj.size()= "<<_spj.size()
        //             <<std::endl;
        //}
        _spj.resizeNoInitialize(offset+_tc.size());
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<_tc.size(); i++){
            _spj[offset+i].copyFromMoment(_tc[i].mom_);
        }
    }    
    template<class Ttreecell, class Tspj>
    inline void AddMomentAsSpImpl(TagForceShort,
                                  const ReallocatableArray<Ttreecell> & _tc,
                                  const S32 offset,
                                  ReallocatableArray<Tspj> & _spj){
        // do nothing
    }


    
    // for long force
    template<class Ttp, class Tepj, class Tspj>
    inline void SetLocalEssentialTreeToGlobalTreeImpl
    (const ReallocatableArray<Tepj> & epj_org,
     const ReallocatableArray<Tspj> & spj_org,
     const S32 n_loc,
     ReallocatableArray<Ttp> & tp_glb,
     const bool flag_reuse = false){
        const S32 n_loc_ep = epj_org.size();
        const S32 n_loc_ep_sp = n_loc_ep + spj_org.size(); 
        tp_glb.resizeNoInitialize( n_loc_ep_sp );
        if(!flag_reuse){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
            {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
                for(S32 i=n_loc; i<n_loc_ep; i++){
                    tp_glb[i].setFromEP(epj_org[i], i);
                }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
                for(S32 i=n_loc_ep; i<n_loc_ep_sp; i++){
                    const S32 i_src = i-n_loc_ep;
                    tp_glb[i].setFromSP(spj_org[i_src], i_src);
                }
            }
        }
    }

    template<class Ttp, class Tepj>
    inline void SetLocalEssentialTreeToGlobalTreeImpl
    (
     const ReallocatableArray<Tepj> & epj_org,
     const S32 n_loc,
     ReallocatableArray<Ttp> & tp_glb,
     const bool flag_reuse = false){
        const S32 n_loc_ep = epj_org.size();
        tp_glb.resizeNoInitialize( n_loc_ep );
        if(!flag_reuse){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
            {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
                for(S32 i=n_loc; i<n_loc_ep; i++){
                    tp_glb[i].setFromEP(epj_org[i], i);
                }
            }
        }
    }

    template<class Ttc>
    inline void SetOuterBoxGlobalTreeForLongCutoffRecursive(Ttc tc[],
                                                            const S32 n_leaf_limit,
                                                            const F64 r_cut,
                                                            const S32 adr_tc,
                                                            const F64 tc_hlen,
                                                            const F64vec tc_cen){
        F64 child_hlen = tc_hlen*0.5;
        for(S32 i=0; i<N_CHILDREN; i++){
            F64vec child_cen = tc_cen + SHIFT_CENTER[i]*tc_hlen;
            tc[adr_tc+i].mom_.vertex_out_.high_ = child_cen + (child_hlen+r_cut);
            tc[adr_tc+i].mom_.vertex_out_.low_  = child_cen - (child_hlen+r_cut);
            S32 child_adr_tc = tc[adr_tc+i].adr_tc_;
            if(tc[adr_tc+i].n_ptcl_ <= 0) continue;
            else if(tc[adr_tc+i].isLeaf(n_leaf_limit)) continue;
            else{
                SetOuterBoxGlobalTreeForLongCutoffRecursive(tc, n_leaf_limit, r_cut, child_adr_tc,
                                                            child_hlen, child_cen);
            }
        }
    }
    
    template<class Ttc, class Tepj>
    inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchLongCutoff,
                                                      Ttc tc[],
                                                      Tepj epj[],
                                                      const S32 n_leaf_limit,
                                                      const F64 tc_hlen,
                                                      const F64vec tc_cen){
        F64 r_cut = epj[0].getRSearch();
        tc[0].mom_.vertex_out_.high_ = tc_cen + (tc_hlen+r_cut);
        tc[0].mom_.vertex_out_.low_  = tc_cen - (tc_hlen+r_cut);
        for(S32 i=1; i<N_CHILDREN; i++) tc[i].mom_.vertex_out_.init();
        if( tc[0].n_ptcl_ < 0 || tc[0].isLeaf(n_leaf_limit) ) return;
        S32 adr_tc = N_CHILDREN;
        SetOuterBoxGlobalTreeForLongCutoffRecursive(tc, n_leaf_limit, r_cut, adr_tc,
                                                    tc_hlen, tc_cen);
    }
    template<class Ttc, class Tepj>
    inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchLong,
                                                      Ttc tc[],
                                                      Tepj epj[],
                                                      const S32 n_leaf_limit,
                                                      const F64 tc_hlen,
                                                      const F64vec tc_cen){
        // do nothing
    }
    template<class Ttc, class Tepj>
    inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchLongScatter,
                                               Ttc tc[],
                                               Tepj epj[],
                                               const S32 n_leaf_limit,
                                               const F64 tc_hlen,
                                               const F64vec tc_cen){
        // do nothing
    }
    template<class Ttc, class Tepj>
    inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchLongSymmetry,
                                               Ttc tc[],
                                               Tepj epj[],
                                               const S32 n_leaf_limit,
                                               const F64 tc_hlen,
                                               const F64vec tc_cen){
        // do nothing
    }
    template<class Ttc, class Tepj>
    inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchShortGather,
                                                      Ttc tc[],
                                                      Tepj epj[],
                                                      const S32 n_leaf_limit,
                                                      const F64 tc_hlen,
                                                      const F64vec tc_cen){
        // do nothing
    }
    template<class Ttc, class Tepj>
    inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchShortScatter,
                                                      Ttc tc[],
                                                      Tepj epj[],
                                                      const S32 n_leaf_limit,
                                                      const F64 tc_hlen,
                                                      const F64vec tc_cen){
        // do nothing
    }
    template<class Ttc, class Tepj>
    inline void SetOuterBoxGlobalTreeForLongCutoffTop(TagSearchShortSymmetry,
                                                      Ttc tc[],
                                                      Tepj epj[],
                                                      const S32 n_leaf_limit,
                                                      const F64 tc_hlen,
                                                      const F64vec tc_cen){
        // do nothing
    }

    // it may works only if tp has one component.
    template<class Tsys, class Ttp, class Tepi, class Tepj>
    inline void CopyFpToEpSortedLocalTree(const Tsys & sys,
                                   const ReallocatableArray<Ttp> & tp,
                                   ReallocatableArray<Tepi> &epi_sorted,
                                   ReallocatableArray<Tepj> &epj_sorted){
        const S32 n_loc = sys.getNumberOfParticleLocal();
        for(S32 i=0; i<n_loc; i++){
            const S32 adr = tp[i].adr_ptcl_;
            epi_sorted[i].copyFromFP(sys[adr]);
            epj_sorted[i].copyFromFP(sys[adr]);
        }
    }


    template<class Tcomm>
    void MakeCommTableFor2StepCommuniction(Tcomm & comm_table,
                                           const ReallocatableArray<S32> & n_ep_send_per_proc_1st,
                                           const ReallocatableArray<S32> & n_image_per_proc_1st,
                                           const ReallocatableArray<F64vec> & shift_per_image_1st,
                                           const ReallocatableArray<S32> & n_ep_send_per_image_1st,
                                           const ReallocatableArray<S32> & adr_ep_send_1st,
                                           const ReallocatableArray<S32> & n_ep_recv_per_proc_1st,
                                           const ReallocatableArray<S32> & n_ep_send_per_proc_2nd,
                                           const ReallocatableArray<S32> & n_image_per_proc_2nd,
                                           const ReallocatableArray<F64vec> & shift_per_image_2nd,
                                           const ReallocatableArray<S32> & n_ep_send_per_image_2nd,
                                           const ReallocatableArray<S32> & adr_ep_send_2nd,
                                           const ReallocatableArray<S32> & n_ep_recv_per_proc_2nd){
        //const auto my_rank = Comm::getRank();
        const auto n_proc = Comm::getNumberOfProc();
        ReallocatableArray<S32> n_disp_ep_send_per_proc(n_proc+1, n_proc+1, 1);
        n_disp_ep_send_per_proc[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_ep_send_per_proc[i+1] = n_disp_ep_send_per_proc[i] + n_ep_send_per_proc_1st[i] + n_ep_send_per_proc_2nd[i];
            //if(my_rank==0){
            //    std::cerr<<"n_disp_ep_send_per_proc[i+1]= "<<n_disp_ep_send_per_proc[i+1]<<std::endl;
            //}
        }
        ReallocatableArray<S32> n_disp_image_per_proc(n_proc+1, n_proc+1, 1);
        ReallocatableArray<S32> n_disp_image_per_proc_1st(n_proc+1, n_proc+1, 1);
        ReallocatableArray<S32> n_disp_image_per_proc_2nd(n_proc+1, n_proc+1, 1);
        n_disp_image_per_proc[0] = n_disp_image_per_proc_1st[0] = n_disp_image_per_proc_2nd[0] = 0;
        for(S32 i=0; i<n_proc; i++){
            n_disp_image_per_proc[i+1]     = n_disp_image_per_proc[i]     + n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
            n_disp_image_per_proc_1st[i+1] = n_disp_image_per_proc_1st[i] + n_image_per_proc_1st[i];
            n_disp_image_per_proc_2nd[i+1] = n_disp_image_per_proc_2nd[i] + n_image_per_proc_2nd[i];
        }

        const S32 n_image_tot_1st = shift_per_image_1st.size();
        const S32 n_image_tot_2nd = shift_per_image_2nd.size();
        ReallocatableArray<S32> n_disp_ep_send_per_image_1st(n_image_tot_1st+1, n_image_tot_1st+1, 1);
        ReallocatableArray<S32> n_disp_ep_send_per_image_2nd(n_image_tot_2nd+1, n_image_tot_2nd+1, 1);
        n_disp_ep_send_per_image_1st[0] = n_disp_ep_send_per_image_2nd[0] = 0;
        for(S32 i=0; i<n_image_tot_1st; i++){
            n_disp_ep_send_per_image_1st[i+1] = n_disp_ep_send_per_image_1st[i] + n_ep_send_per_image_1st[i];
        }
        for(S32 i=0; i<n_image_tot_2nd; i++){
            n_disp_ep_send_per_image_2nd[i+1] = n_disp_ep_send_per_image_2nd[i] + n_ep_send_per_image_2nd[i];
        }
        const S32 n_image_tot = n_disp_image_per_proc_1st[n_proc]  + n_disp_image_per_proc_2nd[n_proc];
        comm_table.shift_per_image_.resizeNoInitialize(n_image_tot);
        comm_table.n_ep_per_image_.resizeNoInitialize(n_image_tot);
        const S32 n_send_tot = n_disp_ep_send_per_proc[n_proc];
        comm_table.adr_ep_send_.resizeNoInitialize(n_send_tot);
PS_OMP_PARALLEL_FOR
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
                comm_table.shift_per_image_[image_head+n_image_cnt] = shift_per_image_1st[j];
                comm_table.n_ep_per_image_[image_head+n_image_cnt]  = n_ep_send_per_image_1st[j];
                for(S32 k=ep_head_1st; k<ep_end_1st; k++, n_ep_cnt++){
                    comm_table.adr_ep_send_[ep_head+n_ep_cnt] = adr_ep_send_1st[k];
                }
            }
            const S32 image_head_2nd = n_disp_image_per_proc_2nd[i];
            const S32 image_end_2nd  = n_disp_image_per_proc_2nd[i+1];
            for(S32 j=image_head_2nd; j<image_end_2nd; j++, n_image_cnt++){
                const S32 ep_head_2nd = n_disp_ep_send_per_image_2nd[j];
                const S32 ep_end_2nd  = n_disp_ep_send_per_image_2nd[j+1];
                comm_table.shift_per_image_[image_head+n_image_cnt] = shift_per_image_2nd[j];
                comm_table.n_ep_per_image_[image_head+n_image_cnt]  = n_ep_send_per_image_2nd[j];
                for(S32 k=ep_head_2nd; k<ep_end_2nd; k++, n_ep_cnt++){
                    comm_table.adr_ep_send_[ep_head+n_ep_cnt] = adr_ep_send_2nd[k];
                }
            }
        }
        comm_table.n_ep_send_.resizeNoInitialize(n_proc);
        comm_table.n_ep_recv_.resizeNoInitialize(n_proc);
        comm_table.n_image_per_proc_.resizeNoInitialize(n_proc);
PS_OMP_PARALLEL_FOR
        for(S32 i=0; i<n_proc; i++){
            comm_table.n_ep_send_[i] = n_ep_send_per_proc_1st[i] + n_ep_send_per_proc_2nd[i];
            comm_table.n_ep_recv_[i] = n_ep_recv_per_proc_1st[i] + n_ep_recv_per_proc_2nd[i];
            comm_table.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
        }
PS_OMP_PARALLEL_FOR 
        for(S32 i=0; i<n_proc; i++){
            comm_table.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
        }
        comm_table.n_ep_send_tot_ = comm_table.adr_ep_send_.size();
        S32 n_ep_recv_tot_tmp = 0;
PS_OMP(omp parallel for reduction(+:n_ep_recv_tot_tmp))
        for(S32 i=0; i<n_proc; i++){
            n_ep_recv_tot_tmp += comm_table.n_ep_recv_[i];
        }
        comm_table.n_ep_recv_tot_ = n_ep_recv_tot_tmp;        
    }
}
