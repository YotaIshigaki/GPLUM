#pragma once

#include<particle_simulator.hpp>

namespace ParticleSimulator{

    // list mode
    struct TagMakeInteractionList{};
    struct TagMakeLetList{};
    struct MAKE_LIST_MODE_LET{
        typedef TagMakeLetList list_mode_type;
    };
    struct MAKE_LIST_MODE_INTERACTION{
        typedef TagMakeInteractionList list_mode_type;
    };

    // make id list
    struct TagListContentId{};
    struct TagListContentPtcl{};
    struct TagListContentPtclId{};
    struct LIST_CONTENT_ID{
        typedef TagListContentId   list_content_type;
    };
    struct LIST_CONTENT_PTCL{
        typedef TagListContentPtcl   list_content_type;
    };
    struct LIST_CONTENT_PTCL_ID{
        typedef TagListContentPtclId list_content_type;
    };

    template<class Ttc>
    U32 GetOpenBit(TagSearchLong,
                    const ReallocatableArray<Ttc> & tc_first,
                    const S32 adr_tc,
                    const F64ort & pos_target_box,
                    const F64 r_crit_sq){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            open_bit |= ( (pos_target_box.getDistanceMinSQ(pos) <= r_crit_sq) << i);
        }
        return open_bit;
    }
    template<class Ttc>
    U32 GetOpenBit(TagSearchLongScatter,
                    const ReallocatableArray<Ttc> & tc_first,
                    const S32 adr_tc,
                    const F64ort & pos_target_box,
                    const F64 r_crit_sq){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            open_bit |= ( (pos_target_box.getDistanceMinSQ(pos) <= r_crit_sq) 
                          || pos_target_box.contained(tc_first[adr_tc+i].mom_.getVertexOut()) ) << i;
        }
        return open_bit;
    }
    template<class Ttc>
    U32 GetOpenBit(TagSearchLongSymmetryOneStage,
                   const ReallocatableArray<Ttc> & tc_first,
                   const S32 adr_tc,
                   const F64ort & pos_target_box,
                   const F64 r_crit_sq){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            open_bit |= ( (pos_target_box.getDistanceMinSQ(pos) <= r_crit_sq) 
                          || pos_target_box.contained(tc_first[adr_tc+i].mom_.getVertexOut()) ) << i;
        }
        return open_bit;
    }
    template<class Ttc>
    U32 GetOpenBit(TagSearchLongSymmetryOneStage,
                    const ReallocatableArray<Ttc> & tc_first,
                    const S32 adr_tc,
                    const F64ort & pos_target_box_in,
                    const F64ort & pos_target_box_out,
                    const F64 r_crit_sq){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            open_bit |= ( (pos_target_box_in.getDistanceMinSQ(pos) <= r_crit_sq) 
                          || pos_target_box_in.contained(tc_first[adr_tc+i].mom_.getVertexOut()) 
                          || pos_target_box_out.contained(tc_first[adr_tc+i].mom_.getVertexIn()) ) << i;
        }
        return open_bit;
    }
    template<class Ttc>
    U32 GetOpenBit(TagSearchShortScatter,
                   const ReallocatableArray<Ttc> & tc_first,
                   const S32 adr_tc,
                   const F64ort & pos_target_box,
                   const F64 r_crit_sq){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bit |= ( pos_target_box.contained( tc_first[adr_tc+i].mom_.getVertexOut() ) << i);
        }
        return open_bit;
    }


    template<class Ttc, class Tsp>
    void CopyInfoDistant(TagForceLong,
                         TagListContentId,
                         const ReallocatableArray<Ttc> & tc_first,
                         const S32 adr_tc,
                         const S32 adr_sp,
                         ReallocatableArray<Tsp> & sp_list,
                         ReallocatableArray<S32> & id_sp_list){
        id_sp_list.push_back(adr_sp);
    }
    /*
    template<class Ttp, class Tep, class Tsp>
    void CopyInfoClose(TagSearchLong,
                       TagMakeLetList,
                       TagListContentId,
                       const ReallocatableArray<Ttp> & tp_first,
                       const S32 adr_ptcl,
                       const S32 n_ptcl,
                       const ReallocatableArray<Tep> & ep_first,
                       ReallocatableArray<Tep> & ep_list,
                       ReallocatableArray<S32> & id_ep_list,
                       const ReallocatableArray<Tsp> & sp_first,
                       ReallocatableArray<Tsp> & sp_list,
                       ReallocatableArray<S32> & id_sp_list){
        S32 cnt_adr_ptcl = adr_ptcl;
        id_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            id_ep_list.pushBackNoCheck(cnt_adr_ptcl);
        }
    }

    template<class Ttp, class Tep, class Tsp>
    void CopyInfoClose(TagSearchLong,
                       TagMakeInteractionList,
                       TagListContentId,
                       const ReallocatableArray<Ttp> & tp_first,
                       const S32 adr_ptcl,
                       const S32 n_ptcl,
                       const ReallocatableArray<Tep> & ep_first,
                       ReallocatableArray<Tep> & ep_list,
                       ReallocatableArray<S32> & id_ep_list,
                       const ReallocatableArray<Tsp> & sp_first,
                       ReallocatableArray<Tsp> & sp_list,
                       ReallocatableArray<S32> & id_sp_list){
        S32 cnt_adr_ptcl = adr_ptcl;
        id_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        id_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                id_ep_list.pushBackNoCheck(cnt_adr_ptcl);
            }
            else{
                id_sp_list.pushBackNoCheck(cnt_adr_ptcl);
            }
        }
    }
    */


    template<class Ttp, class Tep, class Tsp>
    void CopyInfoClose(TagForceLong,
                       TagMakeLetList,
                       TagListContentId,
                       const ReallocatableArray<Ttp> & tp_first,
                       const S32 adr_ptcl,
                       const S32 n_ptcl,
                       const ReallocatableArray<Tep> & ep_first,
                       ReallocatableArray<Tep> & ep_list,
                       ReallocatableArray<S32> & id_ep_list,
                       const ReallocatableArray<Tsp> & sp_first,
                       ReallocatableArray<Tsp> & sp_list,
                       ReallocatableArray<S32> & id_sp_list){
        S32 cnt_adr_ptcl = adr_ptcl;
        id_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            id_ep_list.pushBackNoCheck(cnt_adr_ptcl);
        }
    }

    template<class Ttp, class Tep, class Tsp>
    void CopyInfoClose(TagForceLong,
                       TagMakeInteractionList,
                       TagListContentId,
                       const ReallocatableArray<Ttp> & tp_first,
                       const S32 adr_ptcl,
                       const S32 n_ptcl,
                       const ReallocatableArray<Tep> & ep_first,
                       ReallocatableArray<Tep> & ep_list,
                       ReallocatableArray<S32> & id_ep_list,
                       const ReallocatableArray<Tsp> & sp_first,
                       ReallocatableArray<Tsp> & sp_list,
                       ReallocatableArray<S32> & id_sp_list){
        S32 cnt_adr_ptcl = adr_ptcl;
        id_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        id_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                id_ep_list.pushBackNoCheck(cnt_adr_ptcl);
            }
            else{
                id_sp_list.pushBackNoCheck(cnt_adr_ptcl);
            }
        }
    }

    template<class TSM, class TMLM, class TCM, class Ttc, class Ttp, class Tep, class Tsp>
    void MakeListUsingTreeRecursive(const ReallocatableArray<Ttc> & tc_first,
                                    const S32 adr_tc,
                                    const ReallocatableArray<Ttp> & tp_first,
                                    const ReallocatableArray<Tep> & ep_first,
                                    ReallocatableArray<Tep> & ep_list,
                                    ReallocatableArray<S32> & id_ep_list,
                                    const ReallocatableArray<Tsp> & sp_first,
                                    ReallocatableArray<Tsp> & sp_list,
                                    ReallocatableArray<S32> & id_sp_list,
                                    const F64ort & pos_target_box,
                                    const F64 r_crit_sq,
                                    const S32 n_leaf_limit,
                                    const S32 adr_tree_sp_first=0){ // adress of first sp coming from the (global) tree.

        const Ttc * tc_cur = tc_first.getPointer(adr_tc);
        const F64vec pos_tc = tc_cur->mom_.getPos();
        const S32 n_ptcl = tc_cur->n_ptcl_;
        const S32 adr_ptcl_child = tc_cur->adr_ptcl_;
        const S32 adr_tc_child = tc_cur->adr_tc_;
        if( !(tc_cur->isLeaf(n_leaf_limit)) ){ // not leaf
            U32 open_bit = GetOpenBit(typename TSM::search_type(), tc_first, adr_tc_child, pos_target_box, r_crit_sq*0.25);
            for(S32 i=0; i<N_CHILDREN; i++){
                if( tc_first[adr_tc_child+i].n_ptcl_ <= 0) continue;
                else if( (open_bit>>i) & 0x1 ){ // close
                    MakeListUsingTreeRecursive<TSM, TMLM, TCM, Ttc, Ttp, Tep, Tsp>(tc_first, adr_tc_child+i, tp_first, ep_first, ep_list, id_ep_list, sp_first, sp_list, id_sp_list, pos_target_box, r_crit_sq*0.25, n_leaf_limit, adr_tree_sp_first);
                }
                else{ // far
                    CopyInfoDistant(typename TSM::force_type(), typename TCM::list_content_type(), tc_first, adr_tc_child+i, adr_tc_child+i+adr_tree_sp_first, sp_list, id_sp_list);
                }
            }
        }
        else{ //leaf
            CopyInfoClose(typename TSM::force_type(), typename TMLM::list_mode_type(), typename TCM::list_content_type(), tp_first, adr_ptcl_child, n_ptcl, ep_first, ep_list, id_ep_list, sp_first, sp_list, id_sp_list);
        }
    }



    template<class TSM, class TMLM, class TCM, class Ttc, class Ttp, class Tep, class Tsp>
    void MakeListUsingTreeRecursive(const ReallocatableArray<Ttc> & tc_first,
                                    const S32 adr_tc,
                                    const ReallocatableArray<Ttp> & tp_first,
                                    const ReallocatableArray<Tep> & ep_first,
                                    ReallocatableArray<Tep> & ep_list,
                                    ReallocatableArray<S32> & id_ep_list,
                                    const ReallocatableArray<Tsp> & sp_first,
                                    ReallocatableArray<Tsp> & sp_list,
                                    ReallocatableArray<S32> & id_sp_list,
                                    const F64ort & pos_target_box_in,
                                    const F64ort & pos_target_box_out,
                                    const F64 r_crit_sq,
                                    const S32 n_leaf_limit,
                                    const S32 adr_tree_sp_first=0){ // adress of first sp coming from the (global) tree.

        const Ttc * tc_cur = tc_first.getPointer(adr_tc);
        const F64vec pos_tc = tc_cur->mom_.getPos();
        const S32 n_ptcl = tc_cur->n_ptcl_;
        const S32 adr_ptcl_child = tc_cur->adr_ptcl_;
        const S32 adr_tc_child = tc_cur->adr_tc_;
        if( !(tc_cur->isLeaf(n_leaf_limit)) ){ // not leaf
            U32 open_bit = GetOpenBit(typename TSM::search_type(), tc_first, adr_tc_child, pos_target_box_in, pos_target_box_out, r_crit_sq*0.25);
            for(S32 i=0; i<N_CHILDREN; i++){
                if( tc_first[adr_tc_child+i].n_ptcl_ <= 0) continue;
                else if( (open_bit>>i) & 0x1 ){ // close
                    MakeListUsingTreeRecursive<TSM, TMLM, TCM, Ttc, Ttp, Tep, Tsp>(tc_first, adr_tc_child+i, tp_first, ep_first, ep_list, id_ep_list, sp_first, sp_list, id_sp_list, 
                                                                                   pos_target_box_in, pos_target_box_out, r_crit_sq*0.25, n_leaf_limit, adr_tree_sp_first);
                }
                else{ // far
                    CopyInfoDistant(typename TSM::force_type(), typename TCM::list_content_type(), tc_first, adr_tc_child+i, adr_tc_child+i+adr_tree_sp_first, sp_list, id_sp_list);
                }
            }
        }
        else{ //leaf
            CopyInfoClose(typename TSM::force_type(), typename TMLM::list_mode_type(), typename TCM::list_content_type(), tp_first, adr_ptcl_child, n_ptcl, ep_first, ep_list, id_ep_list, sp_first, sp_list, id_sp_list);
        }
    }



}
