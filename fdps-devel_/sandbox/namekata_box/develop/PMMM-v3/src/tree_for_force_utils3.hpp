// In this file, we implement utility functions used for TreeForForce
// for Particle Mesh Multipole.
#include <unordered_map>
#include <key.hpp>

namespace ParticleSimulator {

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
     const S32 adr_tree_sp_first){
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
                         adr_tree_sp_first);
                }
                else{ // far
                    const S32 adr_sp = adr_tree_sp_first + adr_tc_child + i;
                    adr_sp_list.push_back(adr_sp);
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

    ///////////
    // CalcMultipoleMomentOfParticleMeshCell
    ///////////
    template <class TSM, class Ttc, class Ttp, class Tep, class Tsp,
              class real_t, class cplx_t>
    inline void CalcMultipoleMomentOfParticleMeshCell
    (const S32vec & idx,
     const S32 adr_tc,
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
     const F64vec & center,
     const F64ort & pos_my_domain,
     MultipoleMoment<real_t, cplx_t> & mm,
     const S32 p_spj2mm) {
        assert(p_spj2mm >= 0);
        const Ttc * tc_tmp = tc.getPointer(adr_tc);
        const S32 n_ptcl = tc_tmp->n_ptcl_;
        const U32 adr_ptcl = tc_tmp->adr_ptcl_;
        if (n_ptcl > 0) { // PM cell having particles
            ReallocatableArray<S32> adr_ep_list;
            ReallocatableArray<S32> adr_sp_list;

            // Make lists of array index of ep and sp
#if !defined(PARTICLE_SIMULATOR_USE_EPJ_ONLY_TO_EVAL_MM_IN_PMM)
            F64 r_crit_sq = r_crit_sq_on_root_cell;
            for (S32 lev=1; lev<=lev_pm_cell; lev++) r_crit_sq *= 0.25;
            if (r_crit_sq >= 0.0) {
                // This case corresponds to \theta > 0.
                // In this case, we use local SPJs to calculate MM if local SPJs
                // satisfies the opening angle criterion to the PM cells separated
                // by distance specified by icut.

                const F64vec pos = tc_tmp->mom_.getPos();
                const F64ort my_box = GetBoxOfParticleMeshCell(idx, pos_unit_cell, width_pm_cell);
                ReallocatableArray<F64ort> target_box;
                ReallocatableArray<F64vec> shift_vec;
                shift_vec.push_back(F64vec( width_pm_cell.x * (icut+1), 0, 0)); // +x
                shift_vec.push_back(F64vec(-width_pm_cell.x * (icut+1), 0, 0)); // -x
                shift_vec.push_back(F64vec(0,  width_pm_cell.y * (icut+1), 0)); // +y
                shift_vec.push_back(F64vec(0, -width_pm_cell.y * (icut+1), 0)); // -y
                shift_vec.push_back(F64vec(0, 0,  width_pm_cell.z * (icut+1))); // +z
                shift_vec.push_back(F64vec(0, 0, -width_pm_cell.z * (icut+1))); // -z
                for (S32 k=0; k<shift_vec.size(); k++) 
                    target_box.push_back(my_box.shift(shift_vec[k]));
                if (IsOpen(tc, adr_tc, target_box, r_crit_sq)) {
                    MakeList<TSM, Ttc, Ttp, Tep, Tsp>
                        (tc, adr_tc, tp, ep, adr_ep_list, sp, adr_sp_list, target_box,
                         r_crit_sq, n_leaf_limit, lev_leaf_limit, adr_tree_sp_first);
                } else {
                    // In this case, we use SPJ corresponding to this tree cell
                    // for MM calculation.
                    const S32 adr_sp = adr_tree_sp_first + adr_tc;
                    adr_sp_list.push_back(adr_sp);
                }
            } else {
                // This case corresponds to \theta = 0.
                // In this case, we use only EPJs to calculate multipole moments.
                S32 cnt_adr_ptcl = adr_ptcl;
                adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
                adr_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
                for (S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++) {
                    if( GetMSB(tp[cnt_adr_ptcl].adr_ptcl_) == 0){
                        const S32 adr_ep = tp[cnt_adr_ptcl].adr_ptcl_;
                        adr_ep_list.pushBackNoCheck(adr_ep);
                    }
                    else{
                        const S32 adr_sp = ClearMSB(tp[cnt_adr_ptcl].adr_ptcl_);
                        adr_sp_list.pushBackNoCheck(adr_sp);
                    }
                }
            }
#else
            // In this case, we use only ``local" EPJs to calculate multipole moments.
            S32 cnt_adr_ptcl = adr_ptcl;
            adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
            for (S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++) {
                if( GetMSB(tp[cnt_adr_ptcl].adr_ptcl_) == 0){
                    const S32 adr_ep = tp[cnt_adr_ptcl].adr_ptcl_;
                    if (pos_my_domain.contains(ep[adr_ep].getPos())) {
                        adr_ep_list.pushBackNoCheck(adr_ep);
                    }
                }
            }
#endif

            // Calculate multipole moments
            F64 msum {0.0};
            for (S32 i=0; i < adr_ep_list.size(); i++){
                const S32 adr = adr_ep_list[i];
                const F64 charge = ep[adr].getCharge();
                const F64vec pos = ep[adr].getPos();
                msum += charge;
                mm.assign_particle(center, pos, charge);
            }
            MultipoleMoment<real_t, cplx_t> mm_spj;
            mm_spj.alloc(mm.p);
            for (S32 i=0; i < adr_sp_list.size(); i++){
                const S32 adr = adr_sp_list[i];
                const F64 charge = sp[adr].getCharge();
                const F64vec pos = sp[adr].getPos();
                msum += charge;
                mm_spj.clear();
                // Set monopole (l=0)
                mm_spj.buf[0] = charge;
                // Set dipole (l=1)
                if (p_spj2mm >= 1) {
                    const F64vec dipole = GetMyDipole(sp[adr]);
                    mm_spj.buf[1] = 0.5 * dipole.y;
                    mm_spj.buf[2] = - dipole.z;
                    mm_spj.buf[3] = 0.5 * dipole.x;
                }
                // Set quadrupole (l=2)
                if (p_spj2mm >= 2) {
                    const F64mat quad = GetMyQuadrupole(sp[adr]);
                    mm_spj.buf[4] = 0.25 * quad.xy;
                    mm_spj.buf[5] = 0.5  * quad.yz;
                    mm_spj.buf[6] = 0.25 * (2.0 * quad.zz - quad.xx - quad.yy);
                    mm_spj.buf[7] = -0.5 * quad.xz;
                    mm_spj.buf[8] = 0.125 * (quad.xx - quad.yy);
                }
#ifdef PARTICLE_SIMULATOR_PMM_EXPERIMENTAL_FEATURE
                // Set multipole 
                constexpr S32 p_spj = GetMyMultipoleOrder<Tsp>();
                if (p_spj != -1) {
                    constexpr S32 buflen = MultipoleMoment0<p_spj>::length;
                    F64 *buf = new F64[buflen];
                    GetMyMultipole(sp[adr], buflen, buf);
                    const S32 mmlen = mm_spj.buf.size();
                    for (S32 k=0; k<std::min(buflen, mmlen); k++) {
                        mm_spj.buf[k] = buf[k];
                    }
                    delete [] buf;
                }
#endif

                // [Notes]
                // (1) For the data format of MultipoleMoment class,
                //     see the descriptions of Appendix A.6 in Nitadori (2014)
                //     [arXiv:1409.5981v2].
                // (2) In the above, we assume that `dipole` and `quad` are
                //     calculated by the TreeForForce class as
                //     dipole = \sum_{k} q_{k}(r_{k}-r_{g.c.}) and
                //     quad = \sum_{k} q_{k}(r_{k}-r_{g.c.})_{i}(r_{k}-r_{g.c.})_{j},
                //     where
                //     q_{k} is the electric charge of particle k,
                //     r_{k} is the position vector of particle k,
                //     r_{g.c.} is the geometric center of a tree cell,
                //     and i, j take one of x,y,z.
                // M2M transformation
                mm.assign_from_MM(mm_spj, center, pos);
            }
            mm_spj.freeMem();

            // Free memory
            adr_ep_list.freeMem();
            adr_sp_list.freeMem();
        }
    }

    ///////////
    // CalcTotalChargeAndDispersionOfParticleMeshCell
    ///////////
    template <class Ttc, class Ttp, class Tep>
    inline void CalcTotalChargeAndDispersionOfParticleMeshCell
    (const S32 adr_tc,
     const ReallocatableArray<Ttc> & tc,
     const ReallocatableArray<Ttp> & tp,
     const ReallocatableArray<Tep> & ep,
     const F64vec & center,
     const F64ort & pos_my_domain,
     F64 & msum,
     F64 & quad0) {
         const Ttc * tc_tmp = tc.getPointer(adr_tc);
         const S32 n_ptcl = tc_tmp->n_ptcl_;
         const U32 adr_ptcl = tc_tmp->adr_ptcl_;
         if (n_ptcl > 0) { // PM cell having particles
             ReallocatableArray<S32> adr_ep_list;
             adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);

             // Make a list of array index of ep 
             S32 cnt_adr_ptcl = adr_ptcl;
             for (S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++) {
                 if( GetMSB(tp[cnt_adr_ptcl].adr_ptcl_) == 0){ // EPJ
                     const S32 adr_ep = tp[cnt_adr_ptcl].adr_ptcl_;
                     adr_ep_list.pushBackNoCheck(adr_ep);
                 }
             }

             // Calculate total charge and dispersion using the EPJs
             // this process is in charge of.
             msum = quad0 = 0;
             for (S32 i=0; i < adr_ep_list.size(); i++){
                 const S32 adr = adr_ep_list[i];
                 const F64 charge = ep[adr].getCharge();
                 const F64vec pos = ep[adr].getPos();
                 if (pos_my_domain.contains(pos)) {
                     const F64vec dr = pos - center;
                     msum += charge;
                     quad0 += charge * (dr * dr);
                 }
             }

             // Free memory
             adr_ep_list.freeMem();
        }
    }

    ///////////
    // CalcJParticleInfoOfParticleMeshCell
    ///////////
    template <class TSM, class Ttc, class Ttp, class Tep, class Tsp>
    inline void CalcJParticleInfoOfParticleMeshCell
    (S32 adr_tc_level_partition[],
     const ReallocatableArray<Ttc> & tc,
     const ReallocatableArray<Ttp> & tp,
     const ReallocatableArray<Tep> & ep,
     const ReallocatableArray<Tsp> & sp,
     const S32vec & n_cell,
     const S32 lev_pm_cell,
     std::unordered_map<S32, S32> & adr_tc_glb_from_pm_cell_idx) {
        const S32 head = adr_tc_level_partition[lev_pm_cell];
        const S32 next = adr_tc_level_partition[lev_pm_cell+1];
        for (S32 adr_tc=head; adr_tc<next; adr_tc++) {
            const Ttc * tc_tmp = tc.getPointer(adr_tc);
            const S32 n_ptcl = tc_tmp->n_ptcl_;
            const U32 adr_tp = tc_tmp->adr_ptcl_;
            if (n_ptcl == 0) continue;
            else { // PM cell having particles
                const U64 mkey = tp[adr_tp].key_;
                const S32vec idx = MortonKey<typename TSM::key_type>::getPMCellID(mkey);
                if (idx.x < 0 || n_cell.x <= idx.x ||
                    idx.y < 0 || n_cell.y <= idx.y ||
                    idx.z < 0 || n_cell.z <= idx.z) continue;
                const S32 idx_1d = idx.x + n_cell.x * (idx.y + n_cell.y * idx.z);
                adr_tc_glb_from_pm_cell_idx[idx_1d] = adr_tc;
#if 0
                // for debug
                if (idx_1d == 0) {
                    std::cout << "rank = " << Comm::getRank()
                              << " mkey = " << mkey 
                              << " idx.x = " << idx.x
                              << " idx.y = " << idx.y
                              << " idx.z = " << idx.z
                              << std::endl;
                    for (S32 i=0; i<n_ptcl; i++) {
                        const U32 adr_ptcl = tp[adr_tp + i].adr_ptcl_; 
                        if (GetMSB(adr_ptcl) == 0) { // ep
                            std::cout << "EPJ: (" << Comm::getRank()
                                      << ") i = " << i << " "
                                      << ep[adr_ptcl].pos.x << " "
                                      << ep[adr_ptcl].pos.y << " "
                                      << ep[adr_ptcl].pos.z << " "
                                      << ep[adr_ptcl].mass << std::endl;
                        } else { // sp
                            const U32 adr = ClearMSB(adr_ptcl);
                            std::cout << "SPJ: (" << Comm::getRank()
                                      << ") i = " << i << " "
                                      << sp[adr].getPos().x << " "
                                      << sp[adr].getPos().y << " "
                                      << sp[adr].getPos().z << " "
                                      << sp[adr].getCharge() << std::endl;
                        }
                    }
                }
#endif
            }
        }
    }


    ///////////
    // CalcEpiInfoOfParticleMeshCell
    ///////////
    template <class TSM, class Ttc, class Ttp, class Tep, class Tptclinfo>
    inline void CalcEpiInfoOfParticleMeshCell
    (S32 adr_tc_level_partition[],
     const ReallocatableArray<Ttc> & tc,
     const ReallocatableArray<Ttp> & tp,
     const ReallocatableArray<Tep> & ep,
     const S32vec & n_cell,
     const S32 lev_pm_cell,
     std::unordered_map<S32, Tptclinfo> & ip_info) {
        const S32 head = adr_tc_level_partition[lev_pm_cell];
        const S32 next = adr_tc_level_partition[lev_pm_cell+1];
        for (S32 adr_tc=head; adr_tc<next; adr_tc++) {
            const Ttc * tc_tmp = tc.getPointer(adr_tc);
            const S32 n_ptcl = tc_tmp->n_ptcl_;
            const U32 adr_ptcl = tc_tmp->adr_ptcl_;
            if (n_ptcl == 0) continue;
            else { // PM cell having particles
                assert(GetMSB(adr_ptcl) == 0);
                const U64 mkey = tp[adr_ptcl].key_;
                const S32vec idx = MortonKey<typename TSM::key_type>::getPMCellID(mkey);
                if (idx.x < 0 || n_cell.x <= idx.x ||
                    idx.y < 0 || n_cell.y <= idx.y ||
                    idx.z < 0 || n_cell.z <= idx.z) continue;
                const S32 idx_1d = idx.x + n_cell.x * (idx.y + n_cell.y * idx.z);
                Tptclinfo tmp;
                tmp.n_epi = n_ptcl;
                tmp.adr_epi_first = adr_ptcl;
                tmp.epi_first = ep.getPointer(adr_ptcl);
                ip_info[idx_1d] = tmp;
            }
        }
    }

} // END of namespace of ParticleSimulator
