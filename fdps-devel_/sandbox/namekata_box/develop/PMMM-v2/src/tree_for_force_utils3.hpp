// In this file, we implement utility functions used for TreeForForce
// for Particle Mesh Multipole.
#include <unordered_map>
#include <key.hpp>

namespace ParticleSimulator {

    //////
    struct TagChopCellTrue {}; 
    struct TagChopCellFalse {};

    //////
    struct TagCopyIPInfo {}; // copy epi 
    struct TagCopyJPInfo {}; // copy epj

    ///////////
    // Copy functions
    template <class Tkey, class Ttp, class Tep, class Tsp, class Tcellinfo>
    inline void CopyPtclInfo(TagChopCellFalse,
                             TagCopyIPInfo,
                             const S32 n_ptcl,
                             const U32 adr_ptcl, // = tc[].adr_ptcl_
                             const ReallocatableArray<Ttp> & tp,
                             const ReallocatableArray<Tep> & ep,
                             const ReallocatableArray<Tsp> & sp,
                             std::unordered_map<S32vec, Tcellinfo, HashFunctor<S32vec> > & cell_info,
                             const std::vector<S32vec> & list_my_pm_cell) {
        if (n_ptcl > 0) {
            assert(GetMSB(adr_ptcl) == 0);
            const U64 mkey = tp[adr_ptcl].key_;
            const S32vec idx = MortonKey<Tkey>::getPMCellID(mkey);
            auto itr = cell_info.find(idx);   
            if ( itr != cell_info.end()) {
                cell_info[idx].n_epi = n_ptcl;
                cell_info[idx].adr_epi_first = adr_ptcl;
                cell_info[idx].epi_first = ep.getPointer(adr_ptcl);
            } else {
                Tcellinfo tmp;
                tmp.n_epi = n_ptcl;
                tmp.adr_epi_first = adr_ptcl;
                tmp.epi_first = ep.getPointer(adr_ptcl);
                cell_info[idx] = tmp;
            }
        }
    }

    template <class Tkey, class Ttp, class Tep, class Tsp, class Tcellinfo>
    inline void CopyPtclInfo(TagChopCellTrue,
                             TagCopyIPInfo,
                             const S32 n_ptcl,
                             const U32 adr_ptcl, // = tc.adr_ptcl_
                             const ReallocatableArray<Ttp> & tp,
                             const ReallocatableArray<Tep> & ep,
                             const ReallocatableArray<Tsp> & sp,
                             std::unordered_map<S32vec, Tcellinfo, HashFunctor<S32vec> > & cell_info,
                             const std::vector<S32vec> & list_my_pm_cell) {
        if (n_ptcl > 0) {
            assert(GetMSB(adr_ptcl) == 0);
            const U64 mkey = tp[adr_ptcl].key_;
            const S32vec idx = MortonKey<Tkey>::getPMCellID(mkey);
            auto ret = std::find(list_my_pm_cell.begin(),
                                 list_my_pm_cell.end(),
                                 idx);
            if (ret != list_my_pm_cell.end()) {
                auto itr = cell_info.find(idx);   
                if ( itr != cell_info.end()) {
                    cell_info[idx].n_epi = n_ptcl;
                    cell_info[idx].adr_epi_first = adr_ptcl;
                    cell_info[idx].epi_first = ep.getPointer(adr_ptcl);
                } else {
                    Tcellinfo tmp;
                    tmp.n_epi = n_ptcl;
                    tmp.adr_epi_first = adr_ptcl;
                    tmp.epi_first = ep.getPointer(adr_ptcl);
                    cell_info[idx] = tmp;
                }
            }
        }
    }

    template <class Tkey, class Ttp, class Tep, class Tsp, class Tcellinfo>
    inline void CopyPtclInfo(TagChopCellFalse,
                             TagCopyJPInfo,
                             const S32 n_ptcl,
                             const U32 adr_tp, // = tc[].adr_ptcl_
                             const ReallocatableArray<Ttp> & tp,
                             const ReallocatableArray<Tep> & ep,
                             const ReallocatableArray<Tsp> & sp,
                             std::unordered_map<S32vec, Tcellinfo, HashFunctor<S32vec> > & cell_info,
                             const std::vector<S32vec> & list_my_pm_cell) {
        if (n_ptcl > 0) {
            S32 n_epj = 0;
            S32 adr_tp_for_epj_first = -1;
            S32 adr_epj_first = -1;
            S32 n_spj = 0;
            S32 adr_tp_for_spj_first = -1;
            S32 adr_spj_first = -1;
            for (S32 i=0; i<n_ptcl; i++) {
                const U32 adr_ptcl = tp[adr_tp+i].adr_ptcl_;
                if (GetMSB(adr_ptcl) == 0) { // epj
                    n_epj++;
                    if (adr_epj_first == -1) {
                        adr_tp_for_epj_first = adr_tp+i;
                        adr_epj_first = adr_ptcl;
                    }
                } else { // spj
                    n_spj++;
                    if (adr_spj_first == -1) {
                        adr_tp_for_spj_first = adr_tp+i;
                        adr_spj_first = ClearMSB(adr_ptcl);
                    }
                }
            }
            if (n_epj > 0) {
                const U64 mkey = tp[adr_tp_for_epj_first].key_;
                const S32vec idx = MortonKey<Tkey>::getPMCellID(mkey);
                auto itr = cell_info.find(idx);   
                if ( itr != cell_info.end()) {
                    cell_info[idx].n_epj.push_back(n_epj);
                    cell_info[idx].epj_first.push_back(ep.getPointer(adr_epj_first));
                    if (n_spj > 0) {
                        cell_info[idx].n_spj.push_back(n_spj);
                        cell_info[idx].spj_first.push_back(sp.getPointer(adr_spj_first));
                    }
                } else {
                    Tcellinfo tmp;
                    tmp.n_epj.push_back(n_epj);
                    tmp.epj_first.push_back(ep.getPointer(adr_epj_first));
                    if (n_spj > 0) {
                        tmp.n_spj.push_back(n_spj);
                        tmp.spj_first.push_back(sp.getPointer(adr_spj_first));
                    }
                    cell_info[idx] = tmp;
                }
            }
        }
    }

    template <class Tkey, class Ttp, class Tep, class Tsp, class Tcellinfo>
    inline void CopyPtclInfo(TagChopCellTrue,
                             TagCopyJPInfo,
                             const S32 n_ptcl,
                             const U32 adr_tp, // = tc[].adr_ptcl_
                             const ReallocatableArray<Ttp> & tp,
                             const ReallocatableArray<Tep> & ep,
                             const ReallocatableArray<Tsp> & sp,
                             std::unordered_map<S32vec, Tcellinfo, HashFunctor<S32vec> > & cell_info,
                             const std::vector<S32vec> & list_my_pm_cell) {
        if (n_ptcl > 0) {
            S32 n_epj = 0;
            S32 adr_tp_for_epj_first = -1;
            S32 adr_epj_first = -1;
            S32 n_spj = 0;
            S32 adr_tp_for_spj_first = -1;
            S32 adr_spj_first = -1;
            for (S32 i=0; i<n_ptcl; i++) {
                const U32 adr_ptcl = tp[adr_tp+i].adr_ptcl_;
                if (GetMSB(adr_ptcl) == 0) { // epj
                    n_epj++;
                    if (adr_epj_first == -1) {
                        adr_tp_for_epj_first = adr_tp+i;
                        adr_epj_first = adr_ptcl;
                    }
                } else { // spj
                    n_spj++;
                    if (adr_spj_first == -1) { // spj
                        adr_tp_for_spj_first = adr_tp+i;
                        adr_spj_first = ClearMSB(adr_ptcl);
                    }
                }
            }
            if (n_epj > 0) {
                const U64 mkey = tp[adr_tp_for_epj_first].key_;
                const S32vec idx = MortonKey<Tkey>::getPMCellID(mkey);
                auto ret = std::find(list_my_pm_cell.begin(),
                                     list_my_pm_cell.end(),
                                     idx);
                if (ret != list_my_pm_cell.end()) {
                    auto itr = cell_info.find(idx);   
                    if ( itr != cell_info.end()) {
                        cell_info[idx].n_epj.push_back(n_epj);
                        cell_info[idx].epj_first.push_back(ep.getPointer(adr_epj_first));
                        if (n_spj > 0) {
                            cell_info[idx].n_spj.push_back(n_spj);
                            cell_info[idx].spj_first.push_back(sp.getPointer(adr_spj_first));
                        }
                    } else {
                        Tcellinfo tmp;
                        tmp.n_epj.push_back(n_epj);
                        tmp.epj_first.push_back(ep.getPointer(adr_epj_first));
                        if (n_spj > 0) {
                            tmp.n_spj.push_back(n_spj);
                            tmp.spj_first.push_back(sp.getPointer(adr_spj_first));
                        }
                        cell_info[idx] = tmp;
                    }
                }
            }
        }
    }

    template<class Ttc>
    inline U32 GetOpenBit(const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const ReallocatableArray<F64ort> & target_box,
                          const F64 r_crit_sq){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            bool too_close = false;
            for (S32 k=0; k<target_box.size(); k++) {
                const F64 dis = target_box[k].getDistanceMinSq(pos);
                too_close = too_close || (dis <= r_crit_sq);
            }
            open_bit |= (too_close << i);
            open_bit |= (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN); // if true, it should be checked
        }
        return open_bit;
    }

    template<class TSM, class Ttc, class Ttp, class Tep, class Tsp>
    inline void MakeList
    (const ReallocatableArray<Ttc> & tc_first,
     const S32 adr_tc,
     const ReallocatableArray<Ttp> & tp_first,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<S32> & adr_ep_list,
     const ReallocatableArray<Tsp> & sp_first,
     ReallocatableArray<S32> & adr_sp_list,
     const ReallocatableArray<F64ort> & target_box,
     const F64 r_crit_sq,
     const S32 n_leaf_limit,
     const S32 lev_leaf_limit, 
     const S32 adr_tree_sp_first,
     const bool debug_flag){
        const Ttc * tc_cur = tc_first.getPointer(adr_tc);
        const S32 n_ptcl = tc_cur->n_ptcl_;
        const S32 adr_ptcl = tc_cur->adr_ptcl_;
        const S32 adr_tc_child = tc_cur->adr_tc_;
        if( !(tc_cur->isLeaf(n_leaf_limit, lev_leaf_limit)) ){ // not leaf
            U32 open_bit = GetOpenBit(tc_first, adr_tc_child,
                                      target_box, r_crit_sq*0.25);
            for(S32 i=0; i<N_CHILDREN; i++){
                if( !((open_bit>>(i+N_CHILDREN)) & 0x1) ) continue;
                else if( (open_bit>>i) & 0x1 ){ // close
                    MakeList<TSM, Ttc, Ttp, Tep, Tsp>
                        (tc_first, adr_tc_child+i, tp_first,
                         ep_first, adr_ep_list, sp_first, adr_sp_list,
                         target_box, r_crit_sq*0.25, n_leaf_limit, lev_leaf_limit,
                         adr_tree_sp_first, debug_flag);
                }
                else{ // far
                    const S32 adr_sp = adr_tree_sp_first + adr_tc_child + i;
                    adr_sp_list.push_back(adr_sp);
                    //if (debug_flag) {
                    //    std::cout << "[a] adr_sp = " << adr_sp
                    //              << " qj = " << sp_first[adr_sp].getCharge() << std::endl;
                    //}
                }
            }
        }
        else{ //leaf
            S32 cnt_adr_ptcl = adr_ptcl;
            adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
            adr_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
            for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
                if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                    const S32 adr_ep = tp_first[cnt_adr_ptcl].adr_ptcl_;
                    adr_ep_list.pushBackNoCheck(adr_ep);
                }
                else{
                    const S32 adr_sp = ClearMSB(tp_first[cnt_adr_ptcl].adr_ptcl_);
                    adr_sp_list.pushBackNoCheck(adr_sp);
                    //if (debug_flag) {
                    //    std::cout << "[b] adr_sp = " << adr_sp
                    //              << " qj = " << sp_first[adr_sp].getCharge() << std::endl;
                    //}
                }
            }
        }
    }

    template <class Ttc>
    inline bool IsOpen(const ReallocatableArray<Ttc> & tc,
                       const S32 adr_tc,
                       const ReallocatableArray<F64ort> & target_box,
                       const F64 r_crit_sq) {
        bool too_close = false;
        const F64vec pos = tc[adr_tc].mom_.getPos();
        for (S32 i=0; i<target_box.size(); i++) {
            const F64 dis = target_box[i].getDistanceMinSq(pos);
            too_close = too_close || (dis <= r_crit_sq);
        }
        return too_close;
                
    }

    template <class TSM, class Ttc, class Ttp, class Tep, class Tsp, class Tcellinfo,
              class Tchopmode, class Tcopymode>
    inline void CalcParticleMeshCellInfo
    (S32 adr_tc_level_partition[],
     const ReallocatableArray<Ttc> & tc,
     const ReallocatableArray<Ttp> & tp,
     const ReallocatableArray<Tep> & ep,
     const ReallocatableArray<Tsp> & sp,
     const S32 n_leaf_limit,
     const S32 lev_leaf_limit,
     const S32 adr_tree_sp_first,
     const F64 r_crit_sq_on_root_cell,  
     const F64vec & pos_unit_cell,
     const S32 icut,
     const S32 lev_pm_cell,
     const F64vec & width_pm_cell,
     const std::vector<S32vec> & list_my_pm_cell,
     std::unordered_map<S32vec, Tcellinfo, HashFunctor<S32vec> > & cell_info) {
        const S32 head = adr_tc_level_partition[lev_pm_cell];
        const S32 next = adr_tc_level_partition[lev_pm_cell+1];
        for (S32 adr_tc=head; adr_tc<next; adr_tc++) {
            const Ttc * tc_tmp = tc.getPointer(adr_tc);
            const S32 n_ptcl = tc_tmp->n_ptcl_;
            const U32 adr_ptcl = tc_tmp->adr_ptcl_;
            if (n_ptcl == 0) continue;
            else { // PM cell having particles
                if (typeid(Tcopymode) == typeid(TagCopyIPInfo)) {
                    CopyPtclInfo<typename TSM::key_type, Ttp, Tep, Tsp, Tcellinfo>
                        (Tchopmode(), Tcopymode(),
                         n_ptcl, adr_ptcl, tp, ep, sp,
                         cell_info, list_my_pm_cell);
                } else {
#if 0
                    CopyPtclInfo<typename TSM::key_type, Ttp, Tep, Tsp, Tcellinfo>
                        (Tchopmode(), Tcopymode(),
                         n_ptcl, adr_ptcl, tp, ep, sp,
                         cell_info, list_my_pm_cell);
#else
                    F64 r_crit_sq = r_crit_sq_on_root_cell;
                    for (S32 lev=1; lev<=lev_pm_cell; lev++) r_crit_sq *= 0.25;
                    if (r_crit_sq >= 0.0) {
                        // This case corresponds to \theta > 0.
                        // In this case, we use local SPJs to calculate MM if local SPJs
                        // satisfies the opening angle criterion to the PM cells separated
                        // by distance specified by icut.

                        // [1] Check if this tree cell contains epj
                        S32 adr_tp_for_epj_first = -1;
                        S32 adr_epj_first = -1;
                        for (S32 i=0; i<n_ptcl; i++) {
                            const U32 adr_ptcl_tmp = tp[adr_ptcl+i].adr_ptcl_;
                            if (GetMSB(adr_ptcl_tmp) == 0) { // epj
                                adr_tp_for_epj_first = adr_ptcl + i;
                                adr_epj_first = adr_ptcl;
                            }
                        }
                        if (adr_epj_first == -1) continue;

                        // [2] Get the PM cell index
                        const U64 mkey = tp[adr_tp_for_epj_first].key_;
                        const S32vec idx = MortonKey<typename TSM::key_type>::getPMCellID(mkey);
                        if (typeid(Tchopmode) == typeid(TagChopCellTrue)) {
                            auto itr = std::find(list_my_pm_cell.begin(),
                                                 list_my_pm_cell.end(),
                                                 idx);
                            if (itr == list_my_pm_cell.end()) continue;
                        }

                        // [3] Make lists of array index of ep and sp
                        ReallocatableArray<S32> adr_ep_list;
                        ReallocatableArray<S32> adr_sp_list;

                        ReallocatableArray<F64ort> target_box;
                        ReallocatableArray<F64vec> shift_vec;
                        shift_vec.push_back(F64vec( width_pm_cell.x * icut, 0, 0)); // +x
                        shift_vec.push_back(F64vec(-width_pm_cell.x * icut, 0, 0)); // -x
                        shift_vec.push_back(F64vec(0,  width_pm_cell.y * icut, 0)); // +y
                        shift_vec.push_back(F64vec(0, -width_pm_cell.y * icut, 0)); // -y
                        shift_vec.push_back(F64vec(0, 0,  width_pm_cell.z * icut)); // +z
                        shift_vec.push_back(F64vec(0, 0, -width_pm_cell.z * icut)); // -z
                        F64ort box = PMM::GetBoxOfParticleMeshCell(idx, pos_unit_cell, width_pm_cell);
                        for (S32 k=0; k<shift_vec.size(); k++) 
                            target_box.push_back(box.shift(shift_vec[k]));
                        if (IsOpen(tc, adr_tc, target_box, r_crit_sq)) {
                            bool debug_flag = false;
                            if (Comm::getRank() == 0 && idx.x == 0 && idx.y == 0 && idx.z == 0)
                                debug_flag = true;
                            MakeList<TSM, Ttc, Ttp, Tep, Tsp>
                                (tc, adr_tc, tp, ep, adr_ep_list, sp, adr_sp_list, target_box,
                                 r_crit_sq, n_leaf_limit, lev_leaf_limit, adr_tree_sp_first, debug_flag);
                        } else {
                            // In this case, we use SPJ corresponding to this tree cell
                            // for MM calculation.
                            const S32 adr_sp = adr_tree_sp_first + adr_tc;
                            adr_sp_list.push_back(adr_sp);
                        }

                        // [4] Copy to CellInfo
                        Tcellinfo tmp;
                        S32 adr_ep_prev = -1;
                        S32 n_ep = 0;
                        const S32 n_ep_tot = adr_ep_list.size();
                        for (S32 i=0; i<n_ep_tot; i++) {
                            const S32 adr_ep = adr_ep_list[i];
                            if ((adr_ep_prev == -1) ||
                                ((adr_ep_prev != -1) && (adr_ep == (adr_ep_prev+1)))) n_ep++;
                            else {
                                tmp.n_epj.push_back(n_ep);
                                tmp.epj_first.push_back(ep.getPointer(adr_ep_prev - n_ep + 1));
                                n_ep = 1; // reset
                            }
                            if (i == n_ep_tot - 1) { // last element
                                tmp.n_epj.push_back(n_ep);
                                tmp.epj_first.push_back(ep.getPointer(adr_ep - n_ep + 1));
                            }
                            adr_ep_prev = adr_ep;
                        }
                        S32 adr_sp_prev = -1;
                        S32 n_sp = 0;
                        const S32 n_sp_tot = adr_sp_list.size();
                        for (S32 i=0; i<n_sp_tot; i++) {
                            const S32 adr_sp = adr_sp_list[i];
                            if ((adr_sp_prev == -1) ||
                                ((adr_sp_prev != -1) && (adr_sp == (adr_sp_prev+1)))) n_sp++;
                            else {
                                tmp.n_spj.push_back(n_sp);
                                tmp.spj_first.push_back(sp.getPointer(adr_sp_prev - n_sp + 1));
                                n_sp = 1; // reset
                            }
                            if (i == n_sp_tot - 1) { // last elment
                                tmp.n_spj.push_back(n_sp);
                                tmp.spj_first.push_back(sp.getPointer(adr_sp - n_sp + 1));
                            }

                            adr_sp_prev = adr_sp;
                        }
                        auto itr = cell_info.find(idx);
                        if (itr != cell_info.end()) {
                            cell_info[idx].n_epj = tmp.n_epj;
                            cell_info[idx].epj_first = tmp.epj_first;
                            cell_info[idx].n_spj = tmp.n_spj;
                            cell_info[idx].spj_first = tmp.spj_first;
                        } else {
                            cell_info[idx] = tmp;
                        }

                        // [5] Free memory
                        adr_ep_list.freeMem();
                        adr_sp_list.freeMem();
                    } else {
                        // This case corresponds to \theta = 0.
                        // In this case, we record the information of all EPJs and SPJs
                        // received from other processes in CellInfo. 
                        CopyPtclInfo<typename TSM::key_type, Ttp, Tep, Tsp, Tcellinfo>
                            (Tchopmode(), Tcopymode(),
                             n_ptcl, adr_ptcl, tp, ep, sp,
                             cell_info, list_my_pm_cell);
                    }
#endif
                }
            }
        }
    }

} // END of namespace of ParticleSimulator
