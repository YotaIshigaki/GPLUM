#pragma once

#include<particle_simulator.hpp>

namespace ParticleSimulator{

#if 0
    ////////////
    // debug
    bool f_debug_GetOpenBit {false};
    // debug
    ////////////
#endif

    ////////////
    // walk mode
    struct TagWalkModeNormal{};
    struct TagWalkModeNearestImage{};
    struct WALK_MODE_NORMAL{
        typedef TagWalkModeNormal walk_type;
    };
    struct WALK_MODE_NEAREST_IMAGE{
        typedef TagWalkModeNearestImage walk_type;
    };
    // walk mode
    ////////////

    ////////////
    // walk mode
    struct TagChopLeafTrue{};
    struct TagChopLeafFalse{};
    struct TagChopNonleafTrue{};
    struct TagChopNonleafFalse{};
    // walk mode
    ////////////

    ////////////
    // copy info close mode
    struct TagCopyInfoCloseNormal{};
    struct TagCopyInfoCloseNoSp{};
    struct TagCopyInfoCloseWithTpAdrptcl{};
    // copy info close mode
    ////////////
    
    ////////////
    // targe box class (for long)
    template<class TSM>
    struct TargetBox{
        F64ort vertex_;
    };
    template<>
    struct TargetBox<SEARCH_MODE_LONG>{
        F64ort vertex_;
        void set(const F64ort & vertex_in,
                 const F64ort & vertex_out,
                 const S32 icut) { // used in FindScatterParticle
            vertex_ = vertex_in;
        }
        template<class Tipg>
        void set(const Tipg & ipg,
                 const S32 icut){ // used in makeInteractionListImpl
            vertex_ = ipg.vertex_in_;
        }
        unsigned int contains(const F64vec & pos) const {
            return 1; 
        }
        void dump(std::ostream & fout=std::cerr){
            fout<<"vertex_= "<<vertex_<<std::endl;
        }
    };

    template<>
    struct TargetBox<SEARCH_MODE_LONG_SYMMETRY>{
        F64ort vertex_in_;
        F64ort vertex_out_;
        void set(const F64ort & vertex_in,
                 const F64ort & vertex_out,
                 const S32 icut) { // used in FindScatterParticle
            vertex_in_ = vertex_in;
            vertex_out_ = vertex_out;
        }
        template<class Tipg>
        void set(const Tipg & ipg,
                 const S32 icut){ // used in makeInteractionListImpl
            vertex_in_  = ipg.vertex_in_;
            vertex_out_ = ipg.vertex_out_;
        }
        unsigned int contains(const F64vec pos) const {
            return vertex_out_.contains(pos); 
        }
    };

    template<>
    struct TargetBox<SEARCH_MODE_LONG_SCATTER>{
        F64ort vertex_;
        void set(const F64ort & vertex_in,
                 const F64ort & vertex_out,
                 const S32 icut) { // used in FindScatterParticle
            vertex_ = vertex_in;
        }
        template<class Tipg>
        void set(const Tipg & ipg,
                 const S32 icut){ // used in makeInteractionListImpl
            vertex_ = ipg.vertex_in_;
        }
        unsigned int contains(const F64vec & pos) const {
            return 1; 
        }
    };
    
    template<>
    struct TargetBox<SEARCH_MODE_LONG_CUTOFF>{
        F64ort vertex_;
        void set(const F64ort & vertex_in,
                 const F64ort & vertex_out,
                 const S32 icut) { // used in FindScatterParticle
            vertex_ = vertex_in;
        }
        template<class Tipg>
        void set(const Tipg & ipg,
                 const S32 icut){ // used in makeInteractionListImpl
            vertex_ = ipg.vertex_in_;
        }
        unsigned int contains(const F64vec & pos) const {
            return 1; 
        }
        void dump(std::ostream & fout=std::cerr){
            fout<<"vertex_= "<<vertex_<<std::endl;
        }
    };

    template<>
    struct TargetBox<SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>{
        F64ort vertex_in_;
        S32ort vertex_out_;
        void set(const F64ort & vertex_in,
                 const F64ort & vertex_out,
                 const S32 icut) { // used in FindScatterParticle
            vertex_in_  = vertex_in;
            vertex_out_.low_  = MortonKey<TagKeyMeshBased>::getPMCellID(vertex_in.low_)  - S32vec(icut);
            vertex_out_.high_ = MortonKey<TagKeyMeshBased>::getPMCellID(vertex_in.high_) + S32vec(icut);
        }
        template<class Tipg>
        void set(const Tipg & ipg,
                 const S32 icut){ // used in makeInteractionListImpl() & makeInteractionListIndexLong()
            vertex_in_ = ipg.vertex_in_;
            vertex_out_.low_  = MortonKey<TagKeyMeshBased>::getPMCellID(ipg.vertex_in_.low_)  - S32vec(icut);
            vertex_out_.high_ = MortonKey<TagKeyMeshBased>::getPMCellID(ipg.vertex_in_.high_) + S32vec(icut);
        }
        unsigned int contains(const F64vec & pos) const {
            const S32vec pos_tmp = MortonKey<TagKeyMeshBased>::getPMCellID(pos);
            return vertex_out_.contains(pos_tmp); 
        }
        unsigned int contains(const U64 mkey) const {
            // This is faster than the above `contains`. We can use this function
            // only for EPJ and SPJ having a Morton key.
            const S32vec pos_tmp = MortonKey<TagKeyMeshBased>::getPMCellID(mkey);
            return vertex_out_.contains(pos_tmp); 
        }
        void dump(std::ostream & fout=std::cerr){
            fout<<"vertex_in_= "<<vertex_in_<<std::endl;
            fout<<"vertex_out_= "<<vertex_out_<<std::endl;
        }
    };
    // targe box class (for long)
    ////////////
    
    ////////////
    // targe box class (for short)
    template<>
    struct TargetBox<SEARCH_MODE_SCATTER>{
        F64ort vertex_in_;
        template<class Tipg>
        void set(const Tipg & ipg){ // used in makeInteractionListImpl
            vertex_in_ = ipg.vertex_in_;
        }
        template<class Tep>
        inline bool isInEpList(const ReallocatableArray<Tep> & ep_first,
                               const S32 adr) const {
            const F64 dis_sq = vertex_in_.getDistanceMinSQ(ep_first[adr].getPos());
            const F64 r_crit_sq = ep_first[adr].getRSearch() * ep_first[adr].getRSearch();
            return dis_sq <= r_crit_sq;
        }
    };
    
    template<>
    struct TargetBox<SEARCH_MODE_SYMMETRY>{
        F64ort vertex_out_;
        F64ort vertex_in_;
        template<class Tipg>
        void set(const Tipg & ipg){ // used in makeInteractionListImpl
            vertex_out_ = ipg.vertex_out_;
            vertex_in_ = ipg.vertex_in_;
        }
        template<class Tep>
        inline bool isInEpList(const ReallocatableArray<Tep> & ep_first,
                               const S32 adr) const {
            const F64 dis_sq = vertex_in_.getDistanceMinSQ(ep_first[adr].getPos());
            const F64 r_crit_sq = ep_first[adr].getRSearch() * ep_first[adr].getRSearch();
            return (dis_sq <= r_crit_sq) || (vertex_out_.overlapped(ep_first[adr].getPos()));
        }
    };
    template<>
    struct TargetBox<SEARCH_MODE_GATHER>{
        F64ort vertex_out_;
        template<class Tipg>
        void set(const Tipg & ipg){ // used in makeInteractionListImpl
            vertex_out_ = ipg.vertex_out_;
        }
        template<class Tep>
        inline bool isInEpList(const ReallocatableArray<Tep> & ep_first,
                               const S32 adr) const {
            return vertex_out_.overlapped(ep_first[adr].getPos());
        }
    };
    // targe box class (for short)
    ////////////

    
    ///////////
    // IS OPEN
    template<class TSM, class Ttc>
    inline bool IsOpen(TagSearchLong,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> & target_box,
                       const F64 r_crit_sq,
                       const F64vec & len_peri,
                       TagWalkModeNormal){
        return (tc_first[adr_tc].n_ptcl_ > 0);
    }
    
    template<class TSM, class Ttc>
    inline bool IsOpen(TagSearchLongScatter,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> & target_box,
                       const F64 r_crit_sq,
                       const F64vec & len_peri,
                       TagWalkModeNormal){
        return (tc_first[adr_tc].n_ptcl_ > 0);
    }

    template<class TSM, class Ttc>
    inline bool IsOpen(TagSearchLongSymmetry,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> & target_box,
                       const F64 r_crit_sq,
                       const F64vec & len_peri,
                       TagWalkModeNormal){
        return (tc_first[adr_tc].n_ptcl_ > 0);
    }
    
    template<class TSM, class Ttc>
    inline bool IsOpen(TagSearchLongCutoff,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> & target_box,
                       const F64 r_crit_sq,
                       const F64vec & len_peri,
                       TagWalkModeNormal){
        return ( (target_box.vertex_.overlapped(tc_first[adr_tc].mom_.getVertexOut()))
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
    }

    template<class TSM, class Ttc> 
    inline bool IsOpen(TagSearchLongParticleMeshMultipole,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> & target_box,
                       const F64 r_crit_sq,
                       const F64vec & len_peri,
                       TagWalkModeNormal){
        const S32 n_ptcl = tc_first[adr_tc].n_ptcl_;
        if (n_ptcl > 0) {
            const F64ort vertex_in = tc_first[adr_tc].mom_.getVertexIn();
            S32ort box;
            box.low_ = MortonKey<TagKeyMeshBased>::getPMCellID(vertex_in.low_);
            box.high_ = MortonKey<TagKeyMeshBased>::getPMCellID(vertex_in.high_);
            return (target_box.vertex_out_.overlapped(box));
        } else {
            return false;
        }
        // [Notes]
        // The above implementation is almost the same as the implementation
        // shown below, in which vertex_out_ is class F64ort. However,
        // the following implementation is not secure. 
        // return ( (target_box.vertex_out_.overlapped(tc_first[adr_tc].mom_.getVertexIn()))
        //          &&  (tc_first[adr_tc].n_ptcl_ > 0) );
    }
    
    template<class TSM, class Ttc>
    inline bool IsOpen(TagSearchShortScatter,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> & target_box,
                       const F64 r_crit_sq,
                       const F64vec & len_peri,
                       TagWalkModeNormal){
        return ( (target_box.vertex_in_.overlapped(tc_first[adr_tc].mom_.getVertexOut()))
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
    }
    
    template<class TSM, class Ttc>
    inline bool IsOpen(TagSearchShortSymmetry,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> & target_box,
                       const F64 r_crit_sq,
                       const F64vec & len_peri,
                       TagWalkModeNormal){
        return ( ( (target_box.vertex_in_.overlapped(tc_first[adr_tc].mom_.getVertexOut()))
                 || (target_box.vertex_out_.overlapped(tc_first[adr_tc].mom_.getVertexIn())) )
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
    }
    template<class TSM, class Ttc>
    inline bool IsOpen(TagSearchShortGather,
                       const ReallocatableArray<Ttc> & tc_first,
                       const S32 adr_tc,
                       const TargetBox<TSM> & target_box,
                       const F64 r_crit_sq,
                       const F64vec & len_peri,
                       TagWalkModeNormal){
        return ( target_box.vertex_out_.overlapped(tc_first[adr_tc].mom_.getVertexIn())
                 &&  (tc_first[adr_tc].n_ptcl_ > 0) );
    }    
    // IS OPEN
    ///////////
    
    ///////////
    // COPY INFO DISTANT
    template <class TSM, class Tsp>
    inline void CopyInfoDistant(TagForceShort,
                                TagChopNonleafTrue,
                                const S32 adr_sp,
                                const ReallocatableArray<Tsp> & sp_first,
                                ReallocatableArray<S32> & adr_sp_list,
                                const TargetBox<TSM> & target_box){
        // do nothing
    }
    template <class TSM, class Tsp>
    inline void CopyInfoDistant(TagForceShort,
                                TagChopNonleafFalse,
                                const S32 adr_sp,
                                const ReallocatableArray<Tsp> & sp_first,
                                ReallocatableArray<S32> & adr_sp_list,
                                const TargetBox<TSM> & target_box){
        // do nothing
    }
    template <class TSM, class Tsp>
    inline void CopyInfoDistant(TagForceLong,
                                TagChopNonleafTrue,
                                const S32 adr_sp,
                                const ReallocatableArray<Tsp> & sp_first,
                                ReallocatableArray<S32> & adr_sp_list,
                                const TargetBox<TSM> & target_box){
        if (target_box.contains(sp_first[adr_sp].getPos())) {
            adr_sp_list.push_back(adr_sp);
        }
    }
    template <class TSM, class Tsp>
    inline void CopyInfoDistant(TagForceLong,
                                TagChopNonleafFalse,
                                const S32 adr_sp,
                                const ReallocatableArray<Tsp> & sp_first,
                                ReallocatableArray<S32> & adr_sp_list,
                                const TargetBox<TSM> & target_box){
        adr_sp_list.push_back(adr_sp);
    }
    // COPY INFO DISTANT
    ///////////

    ///////////
    // COPY INFO CLOSE
    template<class TSM, class Ttp, class Tep, class Tsp,
             typename std::enable_if<
                 !std::is_same<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>::value
             >::type* = nullptr
     >
    inline void CopyInfoClose(TagForceLong,
                              TagChopLeafTrue,
                              TagCopyInfoCloseNoSp,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box){
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            assert(GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_)==0);
            if (target_box.contains(ep_first[cnt_adr_ptcl].getPos())) {
                adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
            }
        }
    }
    template<class TSM, class Ttp, class Tep, class Tsp,
             typename std::enable_if<
                 std::is_same<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>::value
             >::type* = nullptr
    >
    inline void CopyInfoClose(TagForceLong,
                              TagChopLeafTrue,
                              TagCopyInfoCloseNoSp,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box){
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            assert(GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_)==0);
            const U64 key = tp_first[cnt_adr_ptcl].key_;
            if (target_box.contains(key)) {
                adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
            }
        }
    }
    template<class TSM, class Ttp, class Tep, class Tsp>
    inline void CopyInfoClose(TagForceLong,
                              TagChopLeafFalse,
                              TagCopyInfoCloseNoSp,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box){
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            assert(GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_)==0);
            adr_ep_list.pushBackNoCheck(cnt_adr_ptcl);
            // for debug
            //if (Comm::getRank() == 1) {
            //    const S64 id = ep_first[cnt_adr_ptcl].id; 
            //    if ( id == 1855 || id == 3609 || id == 3649) {
            //        std::cout << "[CopyInfoClose] id = " << id << std::endl;
            //    }
            //}
        }
    }

    template<class TSM, class Ttp, class Tep, class Tsp,
             typename std::enable_if<
                 !std::is_same<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>::value
             >::type* = nullptr
    >
    inline void CopyInfoClose(TagForceLong,
                              TagChopLeafTrue,
                              TagCopyInfoCloseWithTpAdrptcl,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box){
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        adr_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                const S32 adr_ep = tp_first[cnt_adr_ptcl].adr_ptcl_;
                if (target_box.contains(ep_first[adr_ep].getPos())) {
                    adr_ep_list.pushBackNoCheck(adr_ep);
                }
            }
            else{
                const S32 adr_sp = ClearMSB(tp_first[cnt_adr_ptcl].adr_ptcl_);
                if (target_box.contains(sp_first[adr_sp].getPos())) {
                    adr_sp_list.pushBackNoCheck(adr_sp);
                }
            }
        }
    }
    template<class TSM, class Ttp, class Tep, class Tsp,
             typename std::enable_if<
                 std::is_same<TSM, SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE>::value
             >::type* = nullptr
    >
    inline void CopyInfoClose(TagForceLong,
                              TagChopLeafTrue,
                              TagCopyInfoCloseWithTpAdrptcl,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box){
        S32 cnt_adr_ptcl = adr_ptcl;
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        adr_sp_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
            if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                const S32 adr_ep = tp_first[cnt_adr_ptcl].adr_ptcl_;
                const U64 key = tp_first[cnt_adr_ptcl].key_;
                if (target_box.contains(key)) {
                    adr_ep_list.pushBackNoCheck(adr_ep);
                }
            }
            else{
                const S32 adr_sp = ClearMSB(tp_first[cnt_adr_ptcl].adr_ptcl_);
                const U64 key = tp_first[cnt_adr_ptcl].key_;
                if (target_box.contains(key)) {
                    adr_sp_list.pushBackNoCheck(adr_sp);
                }
            }
        }
    }
    template<class TSM, class Ttp, class Tep, class Tsp>
    inline void CopyInfoClose(TagForceLong,
                              TagChopLeafFalse,
                              TagCopyInfoCloseWithTpAdrptcl,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box){
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

    template<class TSM, class Ttp, class Tep, class Tsp>
    inline void CopyInfoClose(TagForceShort,
                              TagChopLeafTrue,
                              TagCopyInfoCloseNoSp,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box){
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++){
            const S32 adr = adr_ptcl + ip;
            if( target_box.isInEpList(ep_first, adr) ){
                adr_ep_list.pushBackNoCheck(adr);
            }
        }
    }
    template<class TSM, class Ttp, class Tep, class Tsp>
    inline void CopyInfoClose(TagForceShort,
                              TagChopLeafFalse,
                              TagCopyInfoCloseNoSp,
                              const ReallocatableArray<Ttp> & tp_first,
                              const S32 adr_ptcl,
                              const S32 n_ptcl,
                              const ReallocatableArray<Tep> & ep_first,
                              ReallocatableArray<S32> & adr_ep_list,
                              const ReallocatableArray<Tsp> & sp_first,
                              ReallocatableArray<S32> & adr_sp_list,
                              const TargetBox<TSM> & target_box){
        adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
        for(S32 ip=0; ip<n_ptcl; ip++){
            const S32 adr = adr_ptcl + ip;
            adr_ep_list.pushBackNoCheck(adr);
        }
    }
    // COPY INFO CLOSE
    ///////////
    
    /////////////////
    // GET OPEN BITS
    template<class TSM, class Ttc>
    inline U32 GetOpenBit(TagSearchLong,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64 r_crit_sq,
                          const F64vec & len_peri,
                          const S32 lev_crit,
                          TagWalkModeNormal){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            const F64 dis = target_box.vertex_.getDistanceMinSq(pos);
            open_bit |= (dis <= r_crit_sq) << i;
            open_bit |= (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN); // if true, it should be checked
        }
        return open_bit;
    }

    template<class TSM, class Ttc>
    inline U32 GetOpenBit(TagSearchLongScatter,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64 r_crit_sq,
                          const F64vec & len_peri,
                          const S32 lev_crit,
                          TagWalkModeNormal){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            //const F64 dis0 = target_box.vertex_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexOut());
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            const F64 dis1 = target_box.vertex_.getDistanceMinSq(pos);
            open_bit |= ( target_box.vertex_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut()) || dis1 <= r_crit_sq) << i;
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0)
                          << (i + N_CHILDREN) ); // if true, it should be checked
        }
        return open_bit;
    }

    template<class TSM, class Ttc>
    inline U32 GetOpenBit(TagSearchLongSymmetry,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64 r_crit_sq,
                          const F64vec & len_peri,
                          const S32 lev_crit,
                          TagWalkModeNormal){
        U32 open_bit = 0;
#if 1
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            const F64 dis1 = target_box.vertex_in_.getDistanceMinSq(pos);
            open_bit |= ( (target_box.vertex_in_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut())
                           || target_box.vertex_out_.overlapped(tc_first[adr_tc+i].mom_.getVertexIn()) )
                          || dis1 <= r_crit_sq) << i;
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0)
                          << (i + N_CHILDREN) ); // if true, it should be checked
            /*
            const F64 dis0 = std::min( target_box.vertex_in_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexOut()),
                                       target_box.vertex_out_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexIn()) );
            if(Comm::getRank()==0){
                if((target_box.vertex_in_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut()))){
                    std::cerr<<"A) dis0= "<<target_box.vertex_in_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexOut())<<std::endl;
                }
                if((target_box.vertex_out_.overlapped(tc_first[adr_tc+i].mom_.getVertexIn()))){
                    std::cerr<<"B) dis1= "<<target_box.vertex_out_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexIn())<<std::endl;
                }
            }
            */
#if 0
            if(Comm::getRank()==0){
                if(target_box.vertex_in_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexOut()) <= 0.0){
                    
                    if(target_box.vertex_in_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut()) == 0){
                        std::cout<<"A) n_ptcl= "<<tc_first[adr_tc+i].n_ptcl_
                                 <<" adr_tc+i= "<<adr_tc+i
                                 <<" tc_first[adr_tc+i].mom_.getVertexOut()= "<<tc_first[adr_tc+i].mom_.getVertexOut()
                                 <<" tc_first[adr_tc+i].mom_.getVertexIn()= "<<tc_first[adr_tc+i].mom_.getVertexIn()
                                 <<std::endl;
                    }
                        /*
                    std::cerr<<"A) over= "<<target_box.vertex_in_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut())
                        //<<" dis= "<<target_box.vertex_in_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexOut())
                        //<<" target_box.vertex_in_= "<<target_box.vertex_in_
                        //<<" tc_first[adr_tc+i].mom_.getVertexOut()= "<<tc_first[adr_tc+i].mom_.getVertexOut()
                             <<" n_ptcl= "<<tc_first[adr_tc+i].n_ptcl_
                             <<std::endl;
                        */
                }
                /*
                if(target_box.vertex_out_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexIn()) <= 0.0){
                    if(target_box.vertex_out_.overlapped(tc_first[adr_tc+i].mom_.getVertexIn()) == 0){
                        std::cerr<<"B) n_ptcl= "<<tc_first[adr_tc+i].n_ptcl_
                                 <<" tc_first[adr_tc+i].mom_.getVertexOut()= "<<tc_first[adr_tc+i].mom_.getVertexOut()
                                 <<std::endl;
                    }
                }
                */
            }
#endif
        }

        /*
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64 dis0 = std::min( target_box.vertex_in_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexOut()),
                                       target_box.vertex_out_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexIn()) );
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            const F64 dis1 = target_box.vertex_in_.getDistanceMinSq(pos);
            open_bit |= ( dis0 <= 0.0 || dis1 <= r_crit_sq) << i;
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0)
                          << (i + N_CHILDREN) ); // if true, it should be checked
        }
        */
        /*
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            const F64 dis = target_box.vertex_in_.getDistanceMinSq(pos);
            open_bit |= (dis <= r_crit_sq) << i;
            open_bit |= (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN); // if true, it should be checked
        }
        */
#else
        for(S32 i=0; i<N_CHILDREN; i++){
            const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
            const F64 dis = target_box.vertex_in_.getDistanceMinSq(pos);
            open_bit |= (dis <= r_crit_sq) << i;
            open_bit |= (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN); // if true, it should be checked
        }
#endif
        return open_bit;
    }
    
    template<class TSM, class Ttc>
    inline U32 GetOpenBit(TagSearchLongCutoff,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64 r_crit_sq,
                          const F64vec & len_peri,
                          const S32 lev_crit,
                          TagWalkModeNormal){
        U32 open_bit = 0;
        if (r_crit_sq >= 0.0) {
            // theta_ > 0 case
            for(S32 i=0; i<N_CHILDREN; i++){
                const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
                const F64 dis = target_box.vertex_.getDistanceMinSq(pos);
                open_bit |= (dis <= r_crit_sq) << i;
                open_bit |= ( ( (target_box.vertex_.overlapped( tc_first[adr_tc+i].mom_.getVertexOut()))
                                && (tc_first[adr_tc+i].n_ptcl_ > 0) )
                              << (i + N_CHILDREN) ); // if true, it should be checked
            }
        }
        else {
            // theta_ = 0 case
            for(S32 i=0; i<N_CHILDREN; i++){
                open_bit |= 1 << i;
                open_bit |= ( ( (target_box.vertex_.overlapped( tc_first[adr_tc+i].mom_.getVertexOut()))
                                && (tc_first[adr_tc+i].n_ptcl_ > 0) )
                              << (i + N_CHILDREN) ); // if true, it should be checked
            }
        }
        return open_bit;
    }

    template<class TSM, class Ttc>
    inline U32 GetOpenBit(TagSearchLongParticleMeshMultipole,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64 r_crit_sq,
                          const F64vec & len_peri,
                          const S32 lev_crit,
                          TagWalkModeNormal){
        U32 open_bit = 0;
        if (r_crit_sq >= 0.0) {
            // theta_ > 0 case
            for(S32 i=0; i<N_CHILDREN; i++){
                const S32 n_ptcl = tc_first[adr_tc+i].n_ptcl_;
                const F64vec pos = tc_first[adr_tc+i].mom_.getPos();
                const F64 dis = target_box.vertex_in_.getDistanceMinSq(pos);
                const F64ort tc_vertex_in = tc_first[adr_tc+i].mom_.getVertexIn();
                if (n_ptcl > 0) {
                    S32ort box;
                    box.low_  = MortonKey<TagKeyMeshBased>::getPMCellID(tc_vertex_in.low_);
                    box.high_ = MortonKey<TagKeyMeshBased>::getPMCellID(tc_vertex_in.high_);
                    open_bit |= ( (dis <= r_crit_sq) ||
                                  (target_box.vertex_in_.overlapped(tc_vertex_in)) ||
                                  (target_box.vertex_out_.doesNotContain(box)) ||
                                  (tc_first[adr_tc+i].level_ < lev_crit)) << i;
                    open_bit |= ( target_box.vertex_out_.overlapped( box )
                                  << (i + N_CHILDREN) ); // if true, it should be checked
#if 0
                    // for debug
                    if (f_debug_GetOpenBit && i == 0) {
                        std::cout << "target_box.vertex_in_.low_  = " << target_box.vertex_in_.low_ << std::endl;
                        std::cout << "target_box.vertex_in_.high_ = " << target_box.vertex_in_.high_ << std::endl;
                        std::cout << "tc_first[].vertex_in.low_   = " << tc_vertex_in.low_ << std::endl;
                        std::cout << "tc_first[].vertex_in.high_  = " << tc_vertex_in.high_ << std::endl;
                        std::cout << "tc_first[].level_ = " << tc_first[adr_tc+i].level_ << std::endl;
                        std::cout << "lev_crit          = " << lev_crit << std::endl;
                        std::cout << "tc_first[].pos = " << pos << std::endl; 
                        std::cout << "dis = " << dis << std::endl; 
                        std::cout << "r_crit_sq = " << r_crit_sq << std::endl;
                        std::cout << "box.low_  = " << box.low_ << std::endl; 
                        std::cout << "box.high_ = " << box.high_ << std::endl;
                        std::cout << "target_box.vertex_out_.low_  = " << target_box.vertex_out_.low_ << std::endl;
                        std::cout << "target_box.vertex_out_.high_ = " << target_box.vertex_out_.high_ << std::endl;
                        std::cout << "Cond. (1) = " << (dis <= r_crit_sq) << std::endl;
                        std::cout << "Cond. (2) = " << (target_box.vertex_in_.overlapped(tc_vertex_in)) << std::endl;
                        std::cout << "Cond. (3) = " << (target_box.vertex_out_.doesNotContain(box)) << std::endl;
                        std::cout << "Cond. (4) = " << (tc_first[adr_tc+i].level_ < lev_crit) << std::endl;
                    }
                    // for debug
#endif
                }
                // [Notes]
                // The above implementation is almost the same as the implementation
                // shown below, in which vertex_out_ is class F64ort. However,
                // the following implementation is not secure. 
                // open_bit |= ( (dis <= r_crit_sq) ||
                //               (target_box.vertex_out_.doesNotContain(tc_vertex_in)) ||
                //               (tc_first[adr_tc+i].level_ < lev_crit)) << i;
                // open_bit |= ( ( (target_box.vertex_out_.overlapped( tc_vertex_in ))
                //                 && (tc_first[adr_tc+i].n_ptcl_ > 0) )
                //               << (i + N_CHILDREN) ); // if true, it should be checked
            }
        }
        else {
            // theta_ = 0 case
            for(S32 i=0; i<N_CHILDREN; i++){
                const S32 n_ptcl = tc_first[adr_tc+i].n_ptcl_;
                if (n_ptcl > 0) {
                    const F64ort tc_vertex_in = tc_first[adr_tc+i].mom_.getVertexIn();
                    S32ort box;
                    box.low_ = MortonKey<TagKeyMeshBased>::getPMCellID(tc_vertex_in.low_);
                    box.high_ = MortonKey<TagKeyMeshBased>::getPMCellID(tc_vertex_in.high_);
                    // Note that `box` cannot be calculated correctly when n_ptcl_ = 0
                    // since tc_vertex_in is not defined for n_ptcl_ = 0.
                    open_bit |= 1 << i;
                    open_bit |= ( target_box.vertex_out_.overlapped( box ) 
                                  << (i + N_CHILDREN) ); // if true, it should be checked
                }
                // [Notes]
                // The above implementation is almost the same as the implementation
                // shown below, in which vertex_out_ is class F64ort. However,
                // the following implementation is not secure. 
                // open_bit |= 1 << i;
                // open_bit |= ( ( (target_box.vertex_out_.overlapped( tc_first[adr_tc+i].mom_.getVertexIn()))
                //                 && (tc_first[adr_tc+i].n_ptcl_ > 0) )
                //               << (i + N_CHILDREN) ); // if true, it should be checked
            }
        }
        return open_bit;
    }
    
    template<class TSM, class Ttc>
    inline U32 GetOpenBit(TagSearchShortScatter,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          const TargetBox<TSM> & target_box,
                          const F64 r_crit_sq,
                          const F64vec & len_peri,
                          const S32 lev_crit,
                          TagWalkModeNormal){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
#if 0
            const F64 dis = target_box.vertex_in_.getDistanceMinSq(tc_first[adr_tc+i].mom_.getVertexOut());
            open_bit |= ( dis <= 0.0 ) << i;
#else
            open_bit |= ( target_box.vertex_in_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut()) ) << i;
#endif
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN) ); // if true, it should be checked
        }
        return open_bit;
    }

    template<class TSM, class Ttc>
    inline U32 GetOpenBit(TagSearchShortSymmetry,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          //const F64ort & pos_target_box,
                          const TargetBox<TSM> & target_box,
                          const F64 r_crit_sq,
                          const F64vec & len_peri,
                          const S32 lev_crit,
                          TagWalkModeNormal){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bit |= ( target_box.vertex_in_.overlapped(tc_first[adr_tc+i].mom_.getVertexOut())
                          || target_box.vertex_out_.overlapped(tc_first[adr_tc+i].mom_.getVertexIn()) ) << i;            
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0)
                          << (i + N_CHILDREN) ); // if true, it should be checked
        }
        return open_bit;
    }

    template<class TSM, class Ttc>
    inline U32 GetOpenBit(TagSearchShortGather,
                          const ReallocatableArray<Ttc> & tc_first,
                          const S32 adr_tc,
                          //const F64ort & pos_target_box,
                          const TargetBox<TSM> & target_box,
                          const F64 r_crit_sq,
                          const F64vec & len_peri,
                          const S32 lev_crit,
                          TagWalkModeNormal){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bit |= ( target_box.vertex_out_.overlapped(tc_first[adr_tc+i].mom_.getVertexIn()) ) << i;
            open_bit |= ( (tc_first[adr_tc+i].n_ptcl_ > 0) << (i + N_CHILDREN) ); // if true, it should be checked
        }
        return open_bit;
    }
    // GET OPEN BITS
    /////////////////

    // for both long mode and short mode
    template<class TSM, class Ttc, class Ttp, class Tep, class Tsp,
             class Twalkmode, class Tchopleaf, class Tcopyinfoclose, class Tchopnonleaf>
    inline void MakeListUsingTreeRecursive
    (const ReallocatableArray<Ttc> & tc_first,
     const S32 adr_tc,
     const ReallocatableArray<Ttp> & tp_first,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<S32> & adr_ep_list,
     const ReallocatableArray<Tsp> & sp_first,
     ReallocatableArray<S32> & adr_sp_list,
     //const F64ort & pos_target_box,
     const TargetBox<TSM> & target_box,
     const F64 r_crit_sq,
     const S32 n_leaf_limit,
     const S32 lev_leaf_limit, 
     const S32 adr_tree_sp_first, // adress of first sp coming from the (global) tree.
     const F64vec & len_peri,
     const S32 lev_crit){
        const Ttc * tc_cur = tc_first.getPointer(adr_tc);
        const S32 n_ptcl = tc_cur->n_ptcl_;
        const S32 adr_ptcl_child = tc_cur->adr_ptcl_;
        const S32 adr_tc_child = tc_cur->adr_tc_;
        if( !(tc_cur->isLeaf(n_leaf_limit, lev_leaf_limit)) ){ // not leaf
            U32 open_bit = GetOpenBit(typename TSM::search_type(), tc_first, adr_tc_child,
                                      target_box, r_crit_sq*0.25, len_peri, lev_crit,
                                      typename Twalkmode::walk_type());
            for(S32 i=0; i<N_CHILDREN; i++){
                if( !((open_bit>>(i+N_CHILDREN)) & 0x1) ) continue;
                else if( (open_bit>>i) & 0x1 ){ // close
                    MakeListUsingTreeRecursive<TSM, Ttc, Ttp, Tep, Tsp,
                                               Twalkmode, Tchopleaf, Tcopyinfoclose, Tchopnonleaf>
                        (tc_first, adr_tc_child+i, tp_first, ep_first, adr_ep_list, sp_first, adr_sp_list,
                         target_box, r_crit_sq*0.25, n_leaf_limit, lev_leaf_limit, adr_tree_sp_first,
                         len_peri, lev_crit);
                }
                else{ // far
                    CopyInfoDistant(typename TSM::force_type(), Tchopnonleaf(),
                                    adr_tc_child+adr_tree_sp_first+i, sp_first, adr_sp_list,
                                    target_box);
#if 0
                    // for debug
                    if (Comm::getRank() == 0 && adr_tree_sp_first > 0) {
                        const F64 charge = tc_first[adr_tc_child+i].mom_.getCharge();
                        const F64vec pos = tc_first[adr_tc_child+i].mom_.getPos();
                        if ((0.0007238193309 < charge) && (charge < 0.000723819331)) {
                            std::cout << "------" << std::endl;
                            std::cout << "adr_tree_sp_first = " << adr_tree_sp_first << std::endl;
                            std::cout << "adr_tc = " << adr_tc << std::endl;
                            std::cout << "adr_tc_child = " << adr_tc_child << std::endl;
                            std::cout << "i = " << i << std::endl;
                            std::cout << "n_ptcl = " << n_ptcl << std::endl;
                            std::cout << "charge = " << charge << std::endl;
                            std::cout << "pos    = " << pos << std::endl;
                            std::cout << "------" << std::endl;
#if 1
                            // Check the content of a tree cell
                            S32 cnt_adr_ptcl = adr_ptcl_child;
                            for(S32 ip=0; ip<n_ptcl; ip++, cnt_adr_ptcl++){
                                if( GetMSB(tp_first[cnt_adr_ptcl].adr_ptcl_) == 0){
                                    const S32 adr_ep = tp_first[cnt_adr_ptcl].adr_ptcl_;
                                    const U64 key = tp_first[cnt_adr_ptcl].key_;
                                    const S64 id = GetMyId(ep_first[adr_ep]);
                                    const F64 charge = ep_first[adr_ep].getCharge();
                                    const F64vec pos = ep_first[adr_ep].getPos();
                                    const S32vec idx = MortonKey<TagKeyMeshBased>::getPMCellID(key);
                                    std::cout << "ip = " << ip
                                              << " id = " << id
                                              << " idx = " << idx
                                              << " pos = " << pos
                                              << " charge = " << charge
                                              << std::endl;
                                } else{
                                    const S32 adr_sp = ClearMSB(tp_first[cnt_adr_ptcl].adr_ptcl_); 
                                    const U64 key = tp_first[cnt_adr_ptcl].key_;                   
                                    const S64 id = -1;
                                    const F64 charge = sp_first[adr_sp].getCharge();
                                    const F64vec pos = sp_first[adr_sp].getPos();
                                    const S32vec idx = MortonKey<TagKeyMeshBased>::getPMCellID(key);
                                    std::cout << "ip = " << ip
                                              << " id = " << id
                                              << " idx = " << idx
                                              << " pos = " << pos
                                              << " charge = " << charge
                                              << std::endl;
                                }                                                                  
                            }                                                                      
                            std::cout << "------" << std::endl;
#endif
                        }
                    }
                    // for debug
#endif
                }
            }
        }
        else{ //leaf
#if 0
            CopyInfoClose(typename TSM::force_type(), tp_first, adr_ptcl_child, n_ptcl,
                          ep_first, adr_ep_list, sp_first, adr_sp_list);
#elif 0
            adr_ep_list.reserveEmptyAreaAtLeast(n_ptcl);
            for(S32 ip=0; ip<n_ptcl; ip++){
                F64 dis_sq = target_box.vertex_in_.getDistanceMinSQ(ep_first[adr_ptcl_child+ip].pos);
                if(dis_sq > ep_first[adr_ptcl_child+ip].getRSearch()*ep_first[adr_ptcl_child+ip].getRSearch()) continue;
                adr_ep_list.pushBackNoCheck(adr_ptcl_child+ip);
            }
#else
            CopyInfoClose(typename TSM::force_type(), Tchopleaf(), Tcopyinfoclose(), tp_first, adr_ptcl_child, n_ptcl,
                          ep_first, adr_ep_list, sp_first, adr_sp_list, target_box);
#endif
        }
    }

    // for long mode
    template<class TSM, class Ttc, class Ttp, class Tep, class Tsp,
             class Twalkmode, class Tchopleaf, class Tcopyinfoclose, class Tchopnonleaf>
    inline void MakeListUsingTreeRecursiveTop
    (const ReallocatableArray<Ttc> & tc_first,
     const S32 adr_tc,
     const ReallocatableArray<Ttp> & tp_first,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<S32> & adr_ep_list,
     const ReallocatableArray<Tsp> & sp_first,
     ReallocatableArray<S32> & adr_sp_list,
     const TargetBox<TSM> & target_box,
     const F64 r_crit_sq,
     const S32 n_leaf_limit,
     const S32 lev_leaf_limit, 
     const S32 adr_tree_sp_first, // adress of first sp coming from the (global) tree.
     const F64vec & len_peri,
     const S32 lev_crit){
        if ((r_crit_sq >= 0.0) ||
            (typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF)) ||
            (typeid(TSM) == typeid(SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE))) {
            // theta_ > 0 case or PS::SEARCH_MODE_LONG_CUTOFF
            if( IsOpen(typename TSM::search_type(), tc_first, adr_tc, target_box,
                       r_crit_sq, len_peri, typename Twalkmode::walk_type()) ){
                MakeListUsingTreeRecursive<TSM, Ttc, Ttp, Tep, Tsp,
                                           Twalkmode, Tchopleaf, Tcopyinfoclose, Tchopnonleaf>
                    (tc_first, adr_tc, tp_first, ep_first, adr_ep_list, sp_first,
                     adr_sp_list, target_box, r_crit_sq, n_leaf_limit, lev_leaf_limit,
                     adr_tree_sp_first, len_peri, lev_crit);
            } 
        }
        else {
            // theta_ = 0 case with the modes other than PS::SEARCH_MODE_LONG_CUTOFF
            // PS::SEARCH_MODE_LONG_PARTICLE_MESH_MULTIPOLE
            const S32 n_ptcl = tc_first[0].n_ptcl_;
            S32 adr_ptcl = tc_first[0].adr_ptcl_;
            CopyInfoClose(typename TSM::force_type(), Tchopleaf(), Tcopyinfoclose(),
                          tp_first, adr_ptcl, n_ptcl,
                          ep_first, adr_ep_list, sp_first, adr_sp_list, target_box);
        }
    }

    // for short mode
    template<class TSM, class Ttc, class Ttp, class Tep,
             class Twalkmode, class Tchopleaf, class Tchopnonleaf>
    inline void MakeListUsingTreeRecursiveTop
    (const ReallocatableArray<Ttc> & tc_first,
     const S32 adr_tc,
     const ReallocatableArray<Ttp> & tp_first,
     const ReallocatableArray<Tep> & ep_first,
     ReallocatableArray<S32> & adr_ep_list,
     const TargetBox<TSM> & target_box,
     const F64 r_crit_sq,
     const S32 n_leaf_limit,
     const S32 lev_leaf_limit, 
     const F64vec & len_peri,
     const S32 lev_crit){
        S32 adr_tree_sp_first = 0;
        ReallocatableArray<SuperParticleBase> sp_first;
        ReallocatableArray<S32> adr_sp_list;
        MakeListUsingTreeRecursive<TSM, Ttc, Ttp, Tep, SuperParticleBase,
                                   Twalkmode, Tchopleaf, TagCopyInfoCloseNoSp, Tchopnonleaf>
            (tc_first, adr_tc, tp_first, ep_first, adr_ep_list, sp_first,
             adr_sp_list, target_box, r_crit_sq, n_leaf_limit, lev_leaf_limit,
             adr_tree_sp_first, len_peri, lev_crit);
        //adr_sp_list.freeMem();
        //sp_first.freeMem();
    }


    
    ////////////////////
    // FOR DOUBLE WALK
    //////////////////////////
    // for short symmetry or gather mode
    template<class Ttc>
    inline void IsOverlapped(const Ttc * tc_first,
                             const S32 adr_tc,
                             const F64ort & pos_box,
                             const S32 n_leaf_limit,
                             const S32 lev_leaf_limit,
                             U32 & is_overlapped){
        if(is_overlapped) return;
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bit |= (tc_first[adr_tc+i].mom_.getVertexOut().overlapped(pos_box) << i);
        }
        for(S32 i=0; i<N_CHILDREN; i++){
            if(is_overlapped) return;
            const S32 adr_tc_child = adr_tc + i;
            const Ttc * tc_child = tc_first + adr_tc_child;
            const S32 n_child = tc_child->n_ptcl_;
            if(n_child == 0) continue;
            if( (open_bit>>i) & 0x1){
                if( !(tc_child->isLeaf(n_leaf_limit,lev_leaf_limit)) ){
                    // not leaf
                    IsOverlapped(tc_first, tc_first[adr_tc_child].adr_tc_,
                                 pos_box, n_leaf_limit, lev_leaf_limit,
                                 is_overlapped);
                }
                else{
                    // leaf and overlapped
                    is_overlapped = 1;
                    return;
                }
            }
        }
    }

    template<class Ttc>
    inline U32 IsOverlappedTop(const ReallocatableArray<Ttc> & tc_first,
                               const F64ort & pos_box,
                               const S32 n_leaf_limit,
                               const S32 lev_leaf_limit){
        U32 is_overlapped = 0;
        S32 adr_tc = N_CHILDREN;
        if( !tc_first[0].isLeaf(n_leaf_limit,lev_leaf_limit) ){
            IsOverlapped(tc_first.getPointer(), adr_tc, pos_box,
                         n_leaf_limit, lev_leaf_limit, is_overlapped);
        }
        else{
            is_overlapped = tc_first[0].mom_.getVertexOut().overlapped(pos_box);
        }
        return is_overlapped;
    }
    
    template<class Ttca, class Ttcb>
    inline U32 GetOpenBitDoubleWalk(const ReallocatableArray<Ttca> & tc_first_A,
                                    const ReallocatableArray<Ttcb> & tc_first_B,
                                    const S32 adr_tc_B,
                                    const F64ort & pos_domain,
                                    const S32 n_leaf_limit_A,
                                    const S32 lev_leaf_limit_A){
        U32 open_bit = 0;
        for(S32 i=0; i<N_CHILDREN; i++){
            open_bit |= (tc_first_B[adr_tc_B+i].mom_.getVertexOut().overlapped(pos_domain)
                         || IsOverlappedTop(tc_first_A, tc_first_B[adr_tc_B+i].mom_.getVertexIn(), n_leaf_limit_A)) << i;
            open_bit |= (tc_first_B[adr_tc_B+i].n_ptcl_ > 0) << (i + N_CHILDREN); // if false, just skip
        }
        return open_bit;
    }
    
    template<class Ttca, class Ttcb, class Tepb>
    inline void MakeListDoubleWalk(const ReallocatableArray<Ttca> & tc_first_A,
                                   const ReallocatableArray<Ttcb> & tc_first_B,
                                   const S32 adr_tc_B,
                                   const ReallocatableArray<Tepb> & ep_first_B,
                                   const F64ort & pos_domain,
                                   const S32 n_leaf_limit_A,
                                   const S32 lev_leaf_limit_A,
                                   const S32 n_leaf_limit_B,
                                   const S32 lev_leaf_limit_B,
                                   ReallocatableArray<S32> & adr_ptcl_send){
        U32 open_bit = GetOpenBitDoubleWalk(tc_first_A, tc_first_B, adr_tc_B, pos_domain,
                                            n_leaf_limit_A, lev_leaf_limit_A);
        for(S32 i=0; i<N_CHILDREN; i++){
            if( !((open_bit>>(i+N_CHILDREN)) & 0x1) ) continue;
            else if( (open_bit>>i) & 0x1){
                const S32 adr_tc_child = adr_tc_B + i;
                const Ttcb * tc_child = tc_first_B.getPointer() + adr_tc_child;
                const S32 n_child = tc_child->n_ptcl_;
                if( !(tc_child->isLeaf(n_leaf_limit_B,lev_leaf_limit_B)) ){
                    MakeListDoubleWalk(tc_first_A, tc_first_B,
                                       tc_first_B[adr_tc_child].adr_tc_,
                                       ep_first_B, pos_domain,
                                       n_leaf_limit_A, lev_leaf_limit_A,
                                       n_leaf_limit_B, lev_leaf_limit_B,
                                       adr_ptcl_send);
                }
                else{
                    if( pos_domain.overlapped(tc_first_B[adr_tc_child].mom_.getVertexOut()) ){
                        //if(Comm::getRank()==0) std::cerr<<"tc_first_B[adr_tc_child].mom_.getVertexOut()= "<<tc_first_B[adr_tc_child].mom_.getVertexOut()<<std::endl;
                        continue;
                    }
                    else{
                        /*
                        if(Comm::getRank()==0){
                            std::cerr<<"pos_domain= "<<pos_domain
                                     <<" tc_first_B[adr_tc_child].mom_.getVertexOut()= "
                                     <<tc_first_B[adr_tc_child].mom_.getVertexOut()
                                     <<" tc_first_B[adr_tc_child].mom_.getVertexIn()= "
                                     <<tc_first_B[adr_tc_child].mom_.getVertexIn()
                                     <<" overlapped= "
                                     <<IsOverlappedTop(tc_first_A, tc_first_B[adr_tc_child].mom_.getVertexIn(),
                                                       n_leaf_limit_A, lev_leaf_limit_A)
                                     <<std::endl;
                        }
                        */                       
                        S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                        for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                            /*
                            if(Comm::getRank()==0){
                                std::cerr<<"adr_ptcl_tmp= "<<adr_ptcl_tmp<<std::endl;
                            }
                            */
                            adr_ptcl_send.push_back(adr_ptcl_tmp);
                        }
                    }
                    /*
                    S32 adr_ptcl_tmp = tc_child->adr_ptcl_;
                    for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                        if( pos_domain.overlapped(tc_first_B[adr_tc_child].mom_.getVertexOut()) ) continue;
                        adr_ptcl_send.push_back(adr_ptcl_tmp);
                    }
                    */
                }
            }
        }
    }

    
    template<class Ttca, class Ttcb, class Tepb>
    inline void MakeListDoubleWalkTop
    (const ReallocatableArray<Ttca> & tc_first_A, // tree of reciving ptcls
     const ReallocatableArray<Ttcb> & tc_first_B, // tree of assigned ptcls
     const ReallocatableArray<Tepb> & ep_first_B,
     const F64ort & pos_domain,
     const S32 n_leaf_limit_A,
     const S32 lev_leaf_limit_A,
     const S32 n_leaf_limit_B,
     const S32 lev_leaf_limit_B,
     ReallocatableArray<S32> & adr_ptcl_send){
        const S32 adr_tc_B = N_CHILDREN;
        if( !tc_first_B[0].isLeaf(n_leaf_limit_B,lev_leaf_limit_B) ){
            MakeListDoubleWalk(tc_first_A, tc_first_B, adr_tc_B,
                               ep_first_B, pos_domain,
                               n_leaf_limit_A, lev_leaf_limit_A,
                               n_leaf_limit_B, lev_leaf_limit_B,
                               adr_ptcl_send);
        }
        else{
#if 0
            S32 adr_ptcl_tmp = tc_first_B[0].adr_ptcl_;
            const S32 n_child = tc_first_B[0].n_ptcl_;
#else
            // original
            if( tc_first_B[0].mom_.getVertexOut().overlapped(pos_domain) ){
                // already sent
                return;
            }
            else{
                S32 adr_ptcl_tmp = tc_first_B[0].adr_ptcl_;
                const S32 n_child = tc_first_B[0].n_ptcl_;
                // not sent yet (might contain distant ptcls)
                // this is sufficient condition
                for(S32 ip=0; ip<n_child; ip++, adr_ptcl_tmp++){
                    adr_ptcl_send.push_back(adr_ptcl_tmp);
                }
            }
#endif
        }
    }
}
