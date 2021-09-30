#pragma once

#include<particle_simulator.hpp>

namespace ParticleSimulator{
    
    void CopyInfoDistant(TagForceLong,
                         const S32 adr_sp,
                         ReallocatableArray<S32> & adr_sp_list){
        adr_sp_list.push_back(adr_sp);
    }

    template<class Ttp, class Tep, class Tsp>
    void CopyInfoClose(TagForceLong,
                       const ReallocatableArray<Ttp> & tp_first,
                       const S32 adr_ptcl,
                       const S32 n_ptcl,
                       const ReallocatableArray<Tep> & ep_first,
                       ReallocatableArray<S32> & adr_ep_list,
                       const ReallocatableArray<Tsp> & sp_first,
                       ReallocatableArray<S32> & adr_sp_list){
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        adr_sp_list.reserveEmptyAreaAtLeast(n_ptcl);

        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
#ifdef REDUCE_MEMORY
            U32 adr = tp_first[cnt_adr_ptcl].adr_ptcl_;
            if( GetMSB(adr) == 0){
                adr_ep_list.pushBackNoCheck(adr);
            }
            else{
                adr_sp_list.pushBackNoCheck(ClearMSB(adr));
            }
#else
            if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
            }
            else{
                adr_sp_list.pushBackNoCheck(cnt_adr_ptcl);
            }
#endif
        }
    }
    
    template<class Ttc>
    U32 GetOpenBit(TagSearchLongScatter,
                   const ReallocatableArray<Ttc> & tc_first,
                   const S32 adr_tc,
                   const F64ort & pos_target_box,
                   const F64 r_crit_sq,
                   const F64 len_peri_x){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            if(tc_first[adr_tc+i].n_ptcl_ <= 0) continue;
            const F64vec pos = tc_first[adr_tc+i].mom_.getPosForTree();
            const F64 dis0 = GetDistanceMinSqPeriodicX(pos_target_box, pos, len_peri_x);
            const F64 dis1 = GetDistanceMinSqPeriodicX(pos_target_box, tc_first[adr_tc+i].mom_.getVertexOut(), len_peri_x);
            open_bit |= ( (dis0 <= r_crit_sq) || (dis1 <= 0.0) ) << i;
        }
        return open_bit;
    }

    template<class TSM, class Ttc, class Ttp, class Tep, class Tsp>
    void MakeListUsingTreeRecursive(const ReallocatableArray<Ttc> & tc_first,
                                    const S32 adr_tc,
                                    const ReallocatableArray<Ttp> & tp_first,
                                    const ReallocatableArray<Tep> & ep_first,
                                    ReallocatableArray<S32> & adr_ep_list,
                                    const ReallocatableArray<Tsp> & sp_first,
                                    ReallocatableArray<S32> & adr_sp_list,
                                    const F64ort & pos_target_box,
                                    const F64 r_crit_sq,
                                    const S32 n_leaf_limit,
                                    const S32 adr_tree_sp_first, // adress of first sp coming from the (global) tree.
                                    const F64 len_peri_x){ 
        const Ttc * tc_cur = tc_first.getPointer(adr_tc);
        //const F64vec pos_tc = tc_cur->mom_.getPosForTree();
        const S32 n_ptcl = tc_cur->n_ptcl_;
        const S32 adr_ptcl_child = tc_cur->adr_ptcl_;
        const S32 adr_tc_child = tc_cur->adr_tc_;
        if( !(tc_cur->isLeaf(n_leaf_limit)) ){ // not leaf
            U32 open_bit = GetOpenBit(typename TSM::search_type(), tc_first, adr_tc_child, pos_target_box, r_crit_sq*0.25, len_peri_x);
            for(S32 i=0; i<N_CHILDREN; i++){
                if( tc_first[adr_tc_child+i].n_ptcl_ <= 0) continue;
                else if( (open_bit>>i) & 0x1 ){ // close
                    MakeListUsingTreeRecursive<TSM, Ttc, Ttp, Tep, Tsp>
                        (tc_first, adr_tc_child+i, tp_first, ep_first, adr_ep_list, sp_first, adr_sp_list,
                         pos_target_box, r_crit_sq*0.25, n_leaf_limit, adr_tree_sp_first, len_peri_x);
                }
                else{ // far
                    CopyInfoDistant(typename TSM::force_type(), adr_tc_child+adr_tree_sp_first+i, adr_sp_list);
                }
            }
        }
        else{ //leaf
            CopyInfoClose(typename TSM::force_type(), tp_first, adr_ptcl_child, n_ptcl, ep_first, adr_ep_list, sp_first, adr_sp_list);
        }
    }
    
    // for LET
    template<class Ttc>
    bool IsOpen(TagSearchLongScatter,
                const ReallocatableArray<Ttc> & tc_first,
                const S32 adr_tc,
                const F64ort & pos_target_box,
                const F64 r_crit_sq,
                const F64 len_peri_x){
        static const F64 PI = 4.0*atan(1.0);
        const F64vec pos = tc_first[adr_tc].mom_.getPosForTree();
        const F64 dis0 = GetDistanceMinSqPeriodicX(pos_target_box, pos, len_peri_x);
        const F64 dis1 = GetDistanceMinSqPeriodicX(pos_target_box, tc_first[adr_tc].mom_.getVertexOut(), len_peri_x);
        bool ret = ( (dis0 <= r_crit_sq) || (dis1 <= 0.0) );
        return ret;
    }

    template<class TSM, class Ttc, class Ttp, class Tep, class Tsp>
    void MakeListUsingTreeRecursiveTop(const ReallocatableArray<Ttc> & tc_first,
                                       const S32 adr_tc,
                                       const ReallocatableArray<Ttp> & tp_first,
                                       const ReallocatableArray<Tep> & ep_first,
                                       ReallocatableArray<S32> & adr_ep_list,
                                       const ReallocatableArray<Tsp> & sp_first,
                                       ReallocatableArray<S32> & adr_sp_list,
                                       const F64ort & pos_target_box,
                                       const F64 r_crit_sq,
                                       const S32 n_leaf_limit,
                                       const S32 adr_tree_sp_first, // adress of first sp coming from the (global) tree.
                                       const F64 len_peri_x){
        if( !IsOpen(typename TSM::search_type(), tc_first, adr_tc, pos_target_box, r_crit_sq, len_peri_x) ){
            S32 n_ptcl = tc_first[adr_tc].n_ptcl_;
            if(n_ptcl > 0){
                CopyInfoDistant(typename TSM::force_type(), adr_tc, adr_sp_list);
            }
            return;
        }
        MakeListUsingTreeRecursive<TSM, Ttc, Ttp, Tep, Tsp>
            (tc_first, adr_tc, tp_first, ep_first, adr_ep_list, sp_first, adr_sp_list,
             pos_target_box, r_crit_sq, n_leaf_limit, adr_tree_sp_first, len_peri_x);
    }
}
