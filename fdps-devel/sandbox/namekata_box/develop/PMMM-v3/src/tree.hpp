#pragma once
#include<key.hpp>

namespace ParticleSimulator{

    class TreeParticle{
    public:
        U64 key_;
        U32 adr_ptcl_; // U32 because its MSB is used for distinguish if it is EP or SP
        TreeParticle() : key_(SetMSB( U64(0) )), adr_ptcl_(SetMSB( U32(0) )) {}

        template<class Tfp, class Tkey>
        void setFromFP(const Tfp & fp, const U32 adr){
            //key_ = MortonKey<DIMENSION>::getKey( fp.getPos() );
            key_ = MortonKey<Tkey>::getKey( fp.getPos() );
            adr_ptcl_ = adr;
        }
        template<class Tep, class Tkey>
        void setFromEP(const Tep & ep, const U32 adr){
            //key_ = MortonKey<DIMENSION>::getKey( ep.getPos() );
            key_ = MortonKey<Tkey>::getKey( ep.getPos() );
            adr_ptcl_ = adr;
        }
        template<class Tsp, class Tkey>
        void setFromSP(const Tsp & sp, const U32 adr){
            //key_ = MortonKey<DIMENSION>::getKey( sp.getPos() );
            key_ = MortonKey<Tkey>::getKey( sp.getPos() );
            adr_ptcl_ = SetMSB(adr);
        }

        U64 getKey() const {
            return key_;
        }
        // for DEBUG
        void setKey(const U64 key){
            key_ = key;
        }
        void dump(std::ostream & fout=std::cout) const {
            fout<<std::oct<<"key_="<<key_<<std::dec<<std::endl;
            fout<<"adr_ptcl_="<<adr_ptcl_<<std::endl;
        }
    };

    template<class Tmom>
    class TreeCell{
    public:
        S32 n_ptcl_;
        U32 adr_tc_;
        U32 adr_ptcl_;
        S32 level_;
        Tmom mom_;
        TreeCell(){
            n_ptcl_ = adr_tc_ = adr_ptcl_ = level_ = 0;
            mom_.init();
        }
        void clear(){
            n_ptcl_ = adr_tc_ = adr_ptcl_ = level_ = 0;
            mom_.init();
        }
        void clearMoment(){
            mom_.init();
        }
        bool isLeaf(const S32 n_leaf_limit,
                    const S32 lev_leaf_limit) const {
            return (n_ptcl_ == 0) ||
                   (n_ptcl_ <= n_leaf_limit && level_ >= lev_leaf_limit) ||
                   (level_ == TREE_LEVEL_LIMIT);
        }
        
        // for DEBUG
        void dump(std::ostream & fout=std::cout) const {
            fout<<"n_ptcl_="<<n_ptcl_<<std::endl;
            fout<<"adr_tc_="<<adr_tc_<<std::endl;
            fout<<"adr_ptcl_="<<adr_ptcl_<<std::endl;
            mom_.dump(fout);
        }
        
        template<class Tep>
        void dumpTree(const Tep * const first_ep_ptr,
                      const TreeCell * const first_tc_ptr,
                      const F32vec & center,
                      const F32 half_length,
                      const S32 n_leaf_limit,
                      const S32 lev_leaf_limit,
                      std::ostream & fout=std::cout) const {
            fout<<std::endl;
            fout<<"cell info"<<std::endl;
            dump(fout);
            if( !(this->isLeaf(n_leaf_limit, lev_leaf_limit)) ){
                fout<<"this cell is not a leaf"<<std::endl;
                fout<<"half_length="<<half_length<<std::endl;
                fout<<"center="<<center<<std::endl;
                fout<<"level="<<this->level_<<std::endl;
                const TreeCell * child = first_tc_ptr + adr_tc_;
                for(S32 ic=0; ic<N_CHILDREN; ic++){
                    if((child + ic)->n_ptcl_ <= 0) continue;
                    const Tep * ptcl = first_ep_ptr + adr_ptcl_;
                    for(S32 ip=0; ip<n_ptcl_; ip++, ptcl++){
                        if(!IsInBox(ptcl->getPos(), center, half_length)){
                            fout<<"out of box(Cell)"<<std::endl;
                            fout<<"ptcl->getPos()="<<ptcl->getPos()<<std::endl;
                            fout<<"center="<<center<<std::endl;
                            fout<<"half_length="<<half_length<<std::endl;
                        }
                        else{
                            fout<<"in box(Cell)"<<std::endl;
                            fout<<"ptcl->getPos()="<<ptcl->getPos()<<std::endl;
                            fout<<"center="<<center<<std::endl;
                            fout<<"half_length="<<half_length<<std::endl;
                        }
                    }
                    fout<<"octid="<<ic<<std::endl;
                    (child + ic)->dumpTree
                        (first_ep_ptr, first_tc_ptr, 
                         center+SHIFT_CENTER[ic]*half_length, half_length*0.5, 
                         n_leaf_limit, fout);
                }
            }
            else{
                fout<<"this cell is a leaf"<<std::endl;
                fout<<"half_length="<<half_length<<std::endl;
                fout<<"center="<<center<<std::endl;
                fout<<"level="<<this->level_<<std::endl;
                const Tep * ptcl = first_ep_ptr + adr_ptcl_;
                for(S32 ip=0; ip<n_ptcl_; ip++, ptcl++){
                    if(!IsInBox(ptcl->getPos(), center, half_length)){
                        fout<<"out of box(LeafCell)"<<std::endl;
                        fout<<"ptcl->getPos()="<<ptcl->getPos()<<std::endl;
                        fout<<"center="<<center<<std::endl;
                        fout<<"half_length="<<half_length<<std::endl;
                    }
                    else{
                        fout<<"in box(LeafCell)"<<std::endl;
                        fout<<"ptcl->getPos()="<<ptcl->getPos()<<std::endl;
                        fout<<"center="<<center<<std::endl;
                        fout<<"half_length="<<half_length<<std::endl;
                    }
                    //ptcl->dump(fout);
                    fout<<std::endl;
                }
                fout<<std::endl;
            }
        }
        
        template<class Tep>
        void checkTree(const Tep * const first_ep_ptr,
                       const TreeCell * const first_tc_ptr,
                       const F32vec & center,
                       const F32 half_length,
                       const S32 n_leaf_limit,
                       const S32 lev_leaf_limit,
                       const F32 tolerance,
                       S32 & err,
                       std::ostream & fout=std::cout) const {
            if( !(this->isLeaf(n_leaf_limit,lev_leaf_limit)) ){
                //std::cerr<<"adr_tc_="<<adr_tc_<<std::endl;
                const TreeCell * child = first_tc_ptr + adr_tc_;
                for(S32 ic=0; ic<N_CHILDREN; ic++){
                    //std::cerr<<"(child + ic)->n_ptcl_="<<(child + ic)->n_ptcl_<<std::endl;
                    if((child + ic)->n_ptcl_ <= 0) continue;
                    const Tep * ptcl = first_ep_ptr + adr_ptcl_;
                    for(S32 ip=0; ip<n_ptcl_; ip++, ptcl++){
                        if(!IsInBox(ptcl->getPos(), center, half_length, tolerance)){
                            fout<<"out of box(Cell)"<<std::endl;
                            fout<<"ptcl->getPos()="<<ptcl->getPos()<<std::endl;
                            fout<<"center="<<center<<std::endl;
                            fout<<"half_length="<<half_length<<std::endl;
                            fout<<"(center+half_length)-ptcl->getPos()="<<(center+half_length)-ptcl->getPos()<<std::endl;
                            fout<<"(center-half_length)-ptcl->getPos()="<<(center-half_length)-ptcl->getPos()<<std::endl;
                            err++;
                        }
                    }
                    //std::cerr<<"ic="<<ic<<std::endl;
                    //std::cerr<<"SHIFT_CENTER[ic]="<<SHIFT_CENTER[ic]<<std::endl;
                    //std::cerr<<"center+SHIFT_CENTER[ic]*half_length="<<center+SHIFT_CENTER[ic]*half_length<<std::endl;
                    (child + ic)->checkTree
                        (first_ep_ptr, first_tc_ptr,
                         center+SHIFT_CENTER[ic]*half_length, half_length*0.5,
                         n_leaf_limit, tolerance, err, fout);
                }
            }
            else{
                const Tep * ptcl = first_ep_ptr + adr_ptcl_;
                for(S32 ip=0; ip<n_ptcl_; ip++, ptcl++){
                    if(!IsInBox(ptcl->getPos(), center, half_length, tolerance)){
                        fout<<"out of box(Leaf)"<<std::endl;
                        fout<<"center="<<center<<std::endl;
                        fout<<"half_length="<<half_length<<std::endl;
                        fout<<"ptcl->getPos()="<<ptcl->getPos()<<std::endl;
                        fout<<"(center+half_length)-ptcl->getPos()="<<(center+half_length)-ptcl->getPos()<<std::endl;
                        fout<<"(center-half_length)-ptcl->getPos()="<<(center-half_length)-ptcl->getPos()<<std::endl;
                        err++;
                    }
                }
            }
        }
        
        template<class Tep, class Tsp>
        void checkTreeLongGlobalTree(const Tep * const first_ep_ptr,
                                     const Tsp * const first_sp_ptr,
                                     const TreeParticle * const first_tp_ptr,
                                     const TreeCell * const first_tc_ptr,
                                     const F32vec & center,
                                     const F32 half_length,
                                     const S32 n_leaf_limit,
                                     const S32 lev_leaf_limit,
                                     const F32 tolerance,
                                     S32 & err,
                                     std::ostream & fout=std::cout) const {
            if( !(this->isLeaf(n_leaf_limit,lev_leaf_limit)) ){
                const TreeCell * child = first_tc_ptr + adr_tc_;
                for(S32 ic=0; ic<N_CHILDREN; ic++){
                    if((child + ic)->n_ptcl_ <= 0) continue;
                    const TreeParticle * tp_tmp = first_tp_ptr + adr_ptcl_;
                    for(S32 ip=0; ip<n_ptcl_; ip++, tp_tmp++){
                        F32vec pos_tmp;
                        /*
                          const U32 adr = tp_tmp->adr_ptcl_;
                          if( GetMSB(adr) ) pos_tmp = first_sp_ptr[ClearMSB(adr)].getPos();
                          else pos_tmp = first_ep_ptr[adr].getPos();
                        */
                        if(GetMSB(tp_tmp->adr_ptcl_)) pos_tmp = first_sp_ptr[adr_ptcl_+ip].getPos();
                        else pos_tmp = first_ep_ptr[adr_ptcl_+ip].getPos();
                        //if(Comm::getRank() == 0){
                        if(!IsInBox(pos_tmp, center, half_length, tolerance)){
                            fout<<"out of box(Cell)"<<std::endl;
                            fout<<"pos_tmp="<<pos_tmp<<std::endl;
                            fout<<"center="<<center<<std::endl;
                            fout<<"half_length="<<half_length<<std::endl;
                            //fout<<"adr="<<adr<<std::endl;
                            fout<<"adr_ptcl_+ip="<<adr_ptcl_+ip<<std::endl;
                            fout<<"(center+half_length)-pos_tmp="<<(center+half_length)-pos_tmp<<std::endl;
                            fout<<"(center-half_length)-pos_tmp="<<(center-half_length)-pos_tmp<<std::endl;
                            err++;
                            //}
                        }
                    }
                    /*
                    (child + ic)->checkTreeLongGlobalTree < Tep, Tsp >
                        (first_ep_ptr, first_sp_ptr, first_tp_ptr, first_tc_ptr,
                         center+SHIFT_CENTER[ic]*half_length, half_length*0.5,
                         n_leaf_limit, tolerance, err, fout);
                    */
                    (child + ic)->checkTreeLongGlobalTree
                        (first_ep_ptr, first_sp_ptr, first_tp_ptr, first_tc_ptr,
                         center+SHIFT_CENTER[ic]*half_length, half_length*0.5,
                         n_leaf_limit, tolerance, err, fout);
                }
            }
            else{
                const TreeParticle * tp_tmp = first_tp_ptr + adr_ptcl_;
                for(S32 ip=0; ip<n_ptcl_; ip++, tp_tmp++){
                    F32vec pos_tmp;
                    if(GetMSB(tp_tmp->adr_ptcl_)) pos_tmp = first_sp_ptr[adr_ptcl_+ip].getPos();
                    else pos_tmp = first_ep_ptr[adr_ptcl_+ip].getPos();
                    //if(Comm::getRank() == 0){
                    if(!IsInBox(pos_tmp, center, half_length, tolerance)){
                        fout<<"out of box(Leaf)"<<std::endl;
                        fout<<"center="<<center<<std::endl;
                        fout<<"half_length="<<half_length<<std::endl;
                        //fout<<"adr="<<adr<<std::endl;
                        fout<<"adr_ptcl_+ip="<<adr_ptcl_+ip<<std::endl;
                        fout<<"pos_tmp="<<pos_tmp<<std::endl;
                        fout<<"(center+half_length)-pos_tmp="<<(center+half_length)-pos_tmp<<std::endl;
                        fout<<"(center-half_length)-pos_tmp="<<(center-half_length)-pos_tmp<<std::endl;
                        err++;
                    }
                    // }
                }
            }
        }
    };

#if 1
    template<class Tipg>
    class IPGroup{};
    
    template<>
    class IPGroup<TagIpgLongNormal>{
    public:
        S32 n_ptcl_;
        S32 adr_ptcl_;
        F64ort vertex_in_;
        template<class Ttc>
        void copyFromTC(const Ttc & tc){
            n_ptcl_ = tc.n_ptcl_;
            adr_ptcl_ = tc.adr_ptcl_;
        }
    };

    template<>
    class IPGroup<TagIpgIn>{
    public:
        S32 n_ptcl_;
        S32 adr_ptcl_;
        F64ort vertex_in_;
        template<class Ttc>
        void copyFromTC(const Ttc & tc){
            n_ptcl_ = tc.n_ptcl_;
            adr_ptcl_ = tc.adr_ptcl_;
            vertex_in_ = tc.mom_.vertex_in_;
        }
    };
    
    template<>
    class IPGroup<TagIpgInAndOut>{
    public:
        S32 n_ptcl_;
        S32 adr_ptcl_;
        F64ort vertex_in_;
        F64ort vertex_out_;
        template<class Ttc>
        void copyFromTC(const Ttc & tc){
            n_ptcl_ = tc.n_ptcl_;
            adr_ptcl_ = tc.adr_ptcl_;
            vertex_in_ = tc.mom_.vertex_in_;
            vertex_out_ = tc.mom_.vertex_out_;
        }
    };
    
    template<>
    class IPGroup<TagIpgOut>{
    public:
        S32 n_ptcl_;
        S32 adr_ptcl_;
        F64ort vertex_out_;
        template<class Ttc>
        void copyFromTC(const Ttc & tc){
            n_ptcl_ = tc.n_ptcl_;
            adr_ptcl_ = tc.adr_ptcl_;
            vertex_out_ = tc.mom_.vertex_out_;
        }
    };
#else
    template<class TSM>
    class IPGroup{
    public:
        S32 n_ptcl_;
        S32 adr_ptcl_;
        F64ort vertex_;
        template<class Ttc> 
        void copyFromTC(const Ttc & tc){
            CopyFromTCDummy<TSM, Ttc>()(this, tc);
        }
        //
        template<class TSM2, class Ttc2, class Tdummy=void>
        struct CopyFromTCDummy{
            void operator () (IPGroup * ip, const Ttc2 & tc){
                // for SCATTER, FIXED
                ip->n_ptcl_ = tc.n_ptcl_;
                ip->adr_ptcl_ = tc.adr_ptcl_;
                ip->vertex_ = tc.mom_.vertex_in_;
            }
        };
        template<class Ttc2, class Tdummy>
        struct CopyFromTCDummy<SEARCH_MODE_GATHER, Ttc2, Tdummy>{
            void operator () (IPGroup * ip, const Ttc2 & tc){
                // for GAHTER
                ip->n_ptcl_ = tc.n_ptcl_;
                ip->adr_ptcl_ = tc.adr_ptcl_;
                ip->vertex_ = tc.mom_.vertex_out_;
            }
        };

        template<class Ttc2, class Tdummy>
        struct CopyFromTCDummy<SEARCH_MODE_LONG, Ttc2, Tdummy>{
            void operator () (IPGroup * ip, const Ttc2 & tc){
                ip->n_ptcl_ = tc.n_ptcl_;
                ip->adr_ptcl_ = tc.adr_ptcl_;
            }
        };
        template<class Ttc2, class Tdummy>
        struct CopyFromTCDummy<SEARCH_MODE_LONG_CUTOFF, Ttc2, Tdummy>{
            void operator () (IPGroup * ip, const Ttc2 & tc){
                ip->n_ptcl_ = tc.n_ptcl_;
                ip->adr_ptcl_ = tc.adr_ptcl_;
            }
        };
        // for DEBUG
        void dump(std::ostream & fout = std::cout){
            fout<<"n_ptcl_="<<n_ptcl_<<std::endl;
            fout<<"adr_ptcl_="<<adr_ptcl_<<std::endl;
            fout<<"vertex_="<<vertex_<<std::endl;
        }
    };

    template<>
    class IPGroup<SEARCH_MODE_SYMMETRY>{
    public:
        S32 n_ptcl_;
        S32 adr_ptcl_;
        F64ort vertex_;
        F64ort vertex_in_;
        template<class Ttc> 
        void copyFromTC(const Ttc & tc){
            // for SYMMETRY
            n_ptcl_ = tc.n_ptcl_;
            adr_ptcl_ = tc.adr_ptcl_;
            vertex_ = tc.mom_.vertex_out_;
            vertex_in_ = tc.mom_.vertex_in_;
        }
        // for DEBUG
        void dump(std::ostream & fout = std::cout){
            fout<<"n_ptcl_="<<n_ptcl_<<std::endl;
            fout<<"adr_ptcl_="<<adr_ptcl_<<std::endl;
            fout<<"vertex_ing_="<<vertex_in_<<std::endl;
        }
    };

    template<>
    class IPGroup<SEARCH_MODE_LONG_SYMMETRY>{
    public:
        S32 n_ptcl_;
        S32 adr_ptcl_;
        F64ort vertex_;
        F64ort vertex_out_;
        template<class Ttc> 
        void copyFromTC(const Ttc & tc){
            // for SYMMETRY
            n_ptcl_ = tc.n_ptcl_;
            adr_ptcl_ = tc.adr_ptcl_;
            vertex_ = tc.mom_.vertex_in_;
            vertex_out_ = tc.mom_.vertex_out_;
        }
        // for DEBUG
        void dump(std::ostream & fout = std::cout){
            fout<<"n_ptcl_="<<n_ptcl_<<std::endl;
            fout<<"adr_ptcl_="<<adr_ptcl_<<std::endl;
            fout<<"vertex_="<<vertex_<<std::endl;
            fout<<"vertex_out_="<<vertex_out_<<std::endl;
        }
    };    
#endif
    
    ////////////////////////
    // ESSENTIAL PARTICLE //
    class EPXROnly{
        F64 r_search;
        F64vec pos;
    public:
        F64 getRSearch() const { return r_search;}
        F64vec getPos() const { return pos;}
        void setPos(const F64vec & pos_in) { pos = pos_in;}
        template<class Tep>
        void copyFromEP(const Tep & ep){
            r_search = ep.getRSearch();
            pos = ep.getPos();
        }
        template<class Tptcl>
        const EPXROnly & operator = (const Tptcl & ptcl){
            r_search = ptcl.getRSearch();
            pos = ptcl.getPos();
            return (*this);
        }
    };


    //////////////
    /// MOMENT ///
    //////////////
    //  We implement the following classes in this order:
    //     MomentMonopole
    //     MomentMonopoleInOnly
    //     MomentMonopoleOutOnly
    //     MomentMonopoleInAndOut
    //     MomentMonopoleCutoffScatter
    //     MomentQuadrupole
    //     MomentQuadrupoleInOnly
    //     MomentQuadrupoleOutOnly
    //     MomentQuadrupoleInAndOut
    //     MomentMonopoleGeometricCenter
    //     MomentMonopoleGeometricCenterInOnly
    //     MomentMonopoleGeometricCenterOutOnly
    //     MomentMonopoleGeometricCenterInAndOut
    //     MomentDipoleGeometricCenter
    //     MomentDipoleGeometricCenterInOnly
    //     MomentDipoleGeometricCenterOutOnly
    //     MomentDipoleGeometricCenterInAndOut
    //     MomentQuadrupoleGeometricCenter
    //     MomentQuadrupoleGeometricCenterInOnly
    //     MomentQuadrupoleGeometricCenterOutOnly
    //     MomentQuadrupoleGeometricCenterInAndOut
    //     MomentMultipoleGeometricCenterPMMM
    // After that, there are several ``typedef" to give another names.
    // [TODO]
    //     Maybe, the implementation becomes short
    //     if we use the inheritance feature in C++.

    class MomentMonopole {
    public:
        FSP mass;
        FSPvec pos;
        MomentMonopole(){
            mass = 0.0;
            pos = 0.0;
        }
        MomentMonopole(const FSP m, const FSPvec & p){
            mass = m;
            pos = p;
        }
        void init(){
            mass = 0.0;
            pos = 0.0;
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return mass;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            mass += epj.getCharge();
            pos += epj.getCharge() * epj.getPos();
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){}
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentMonopole & mom){
            mass += mom.mass;
            pos += mom.mass * mom.pos;
        }
        void accumulate2(const MomentMonopole & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
        }
    };

    class MomentMonopoleInOnly{
    public:
        FSP mass;
        FSPvec pos;
        F64ort vertex_in_;
        MomentMonopoleInOnly(){
            mass = 0.0;
            pos = 0.0;
            vertex_in_.init();
        }
        MomentMonopoleInOnly(const FSP m, const FSPvec & p){ 
            mass = m;
            pos = p;
        }
        MomentMonopoleInOnly(const FSP m, const FSPvec & p, const F64ort & in){ 
            mass = m;
            pos = p;
            vertex_in_ = in;
        }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            mass = 0.0;
            pos = 0.0;
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return mass;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->mass += epj.getCharge();
            this->pos += epj.getCharge() * epj.getPos();
            (this->vertex_in_).merge(epj.getPos());
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){}
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentMonopoleInOnly & mom){
            this->mass += mom.mass;
            this->pos += mom.mass * mom.pos;
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        void accumulate2(const MomentMonopoleInOnly & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };

    class MomentMonopoleOutOnly{
    public:
        FSP mass;
        FSPvec pos;
        F64ort vertex_out_;
        MomentMonopoleOutOnly(){
            mass = 0.0;
            pos = 0.0;
            vertex_out_.init();
        }
        MomentMonopoleOutOnly(const FSP m, const FSPvec & p){ 
            mass = m;
            pos = p;
        }
        MomentMonopoleOutOnly(const FSP m, const FSPvec & p, const F64ort & out) {
            mass = m;
            pos = p;
            vertex_out_ = out;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        void init(){
            mass = 0.0;
            pos = 0.0;
            vertex_out_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return mass;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->mass += epj.getCharge();
            this->pos += epj.getCharge() * epj.getPos();
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){}
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentMonopoleOutOnly & mom){
            this->mass += mom.mass;
            this->pos += mom.mass * mom.pos;
            (this->vertex_out_).merge(mom.vertex_out_);
        }
        void accumulate2(const MomentMonopoleOutOnly & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
        }
    };

    class MomentMonopoleInAndOut{
    public:
        FSP mass;
        FSPvec pos;
        F64ort vertex_out_;
        F64ort vertex_in_;
        MomentMonopoleInAndOut(){
            mass = 0.0;
            pos = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        MomentMonopoleInAndOut(const FSP m, const FSPvec & p){ 
            mass = m;
            pos = p;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            mass = 0.0;
            pos = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return mass;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->mass += epj.getCharge();
            this->pos += epj.getCharge() * epj.getPos();
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
            (this->vertex_in_).merge(epj.getPos());
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){}
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentMonopoleInAndOut & mom){
            this->mass += mom.mass;
            this->pos += mom.mass * mom.pos;
            (this->vertex_out_).merge(mom.vertex_out_);
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        void accumulate2(const MomentMonopoleInAndOut & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };

    class MomentMonopoleCutoffScatter{// for P^3T + PM
    public:
        FSP mass;
        FSPvec pos;
        F64ort vertex_out_; // cutoff
        F64ort vertex_out2_; // search ep
        F64ort vertex_in_;
        MomentMonopoleCutoffScatter(){
            mass = 0.0;
            pos = 0.0;
            vertex_out_.init();
            vertex_out2_.init();
            vertex_in_.init();
        }
        MomentMonopoleCutoffScatter(const FSP m, const FSPvec & p){ 
            mass = m;
            pos = p;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        F64ort getVertexOut2() const { return vertex_out2_; }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            mass = 0.0;
            pos = 0.0;
            vertex_out_.init();
            vertex_out2_.init();
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return mass;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->mass += epj.getCharge();
            this->pos += epj.getCharge() * epj.getPos();
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
            (this->vertex_out2_).merge(epj.getPos(), epj.getRSearch2());
            (this->vertex_in_).merge(epj.getPos());
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){}
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentMonopoleCutoffScatter & mom){
            this->mass += mom.mass;
            this->pos += mom.mass * mom.pos;
            (this->vertex_out_).merge(mom.vertex_out_);
            (this->vertex_out2_).merge(mom.vertex_out2_);
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        void accumulate2(const MomentMonopoleCutoffScatter & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };


    class MomentQuadrupole{
    public:
        FSPvec pos;
        FSP mass;
        FSPmat quad;
        void init(){
            pos = 0.0;
            mass = 0.0;
            quad = 0.0;
        }
        MomentQuadrupole(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
        }
        MomentQuadrupole(const FSP m, const FSPvec & p, const FSPmat & q){
            mass = m;
            pos = p;
            quad = q;
        }
        FSPvec getPos() const {
            return pos;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            mass += epj.getCharge();
            pos += epj.getCharge() * epj.getPos();
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = epj.getCharge();
            F64vec ptmp = epj.getPos() - this->pos;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quad.xx += cx * ptmp.x;
            this->quad.yy += cy * ptmp.y;
            this->quad.zz += cz * ptmp.z;
            this->quad.xy += cx * ptmp.y;
            this->quad.xz += cx * ptmp.z;
            this->quad.yz += cy * ptmp.z;
#else
            // under construction
#endif
        }
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentQuadrupole & mom){
            mass += mom.mass;
            pos += mom.mass * mom.pos;
        }
        void accumulate2(const MomentQuadrupole & mom){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 mtmp = mom.mass;
            F64vec ptmp = mom.pos - this->pos;
            F64 cx = mtmp * ptmp.x;
            F64 cy = mtmp * ptmp.y;
            F64 cz = mtmp * ptmp.z;
            this->quad.xx += cx * ptmp.x + mom.quad.xx;
            this->quad.yy += cy * ptmp.y + mom.quad.yy;
            this->quad.zz += cz * ptmp.z + mom.quad.zz;
            this->quad.xy += cx * ptmp.y + mom.quad.xy;
            this->quad.xz += cx * ptmp.z + mom.quad.xz;
            this->quad.yz += cy * ptmp.z + mom.quad.yz;
#else
            // under construction
#endif
        }
    };

    class MomentQuadrupoleInOnly{
    public:
        FSP mass;
        FSPvec pos;
        FSPmat quad;
        F64ort vertex_in_;
        MomentQuadrupoleInOnly(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
            vertex_in_.init();
        }
        MomentQuadrupoleInOnly(const FSP m, const FSPvec & p, const FSPmat & q){ 
            mass = m;
            pos = p;
            quad = q;
        }
        MomentQuadrupoleInOnly(const FSP m, const FSPvec & p, const FSPmat & q, const F64ort & in){ 
            mass = m;
            pos = p;
            quad = q;
            vertex_in_ = in;
        }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return mass;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->mass += epj.getCharge();
            this->pos += epj.getCharge() * epj.getPos();
            (this->vertex_in_).merge(epj.getPos());
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = epj.getCharge();
            F64vec ptmp = epj.getPos() - this->pos;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quad.xx += cx * ptmp.x;
            this->quad.yy += cy * ptmp.y;
            this->quad.zz += cz * ptmp.z;
            this->quad.xy += cx * ptmp.y;
            this->quad.xz += cx * ptmp.z;
            this->quad.yz += cy * ptmp.z;
#else
            // under construction
#endif
        }
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentQuadrupoleInOnly & mom){
            this->mass += mom.mass;
            this->pos += mom.mass * mom.pos;
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        void accumulate2(const MomentQuadrupoleInOnly & mom){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 mtmp = mom.mass;
            F64vec ptmp = mom.pos - this->pos;
            F64 cx = mtmp * ptmp.x;
            F64 cy = mtmp * ptmp.y;
            F64 cz = mtmp * ptmp.z;
            this->quad.xx += cx * ptmp.x + mom.quad.xx;
            this->quad.yy += cy * ptmp.y + mom.quad.yy;
            this->quad.zz += cz * ptmp.z + mom.quad.zz;
            this->quad.xy += cx * ptmp.y + mom.quad.xy;
            this->quad.xz += cx * ptmp.z + mom.quad.xz;
            this->quad.yz += cy * ptmp.z + mom.quad.yz;
#else
            // under construction
#endif
        }
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };

    class MomentQuadrupoleOutOnly{
    public:
        FSP mass;
        FSPvec pos;
        FSPmat quad;
        F64ort vertex_out_;
        MomentQuadrupoleOutOnly(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
            vertex_out_.init();
        }
        MomentQuadrupoleOutOnly(const FSP m, const FSPvec & p, const FSPmat & q){ 
            mass = m;
            pos = p;
            quad = q;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        void init(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
            vertex_out_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return mass;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->mass += epj.getCharge();
            this->pos += epj.getCharge() * epj.getPos();
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = epj.getCharge();
            F64vec ptmp = epj.getPos() - this->pos;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quad.xx += cx * ptmp.x;
            this->quad.yy += cy * ptmp.y;
            this->quad.zz += cz * ptmp.z;
            this->quad.xy += cx * ptmp.y;
            this->quad.xz += cx * ptmp.z;
            this->quad.yz += cy * ptmp.z;
#else
            // under construction
#endif
        }
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentQuadrupoleOutOnly & mom){
            this->mass += mom.mass;
            this->pos += mom.mass * mom.pos;
            (this->vertex_out_).merge(mom.vertex_out_);
        }
        void accumulate2(const MomentQuadrupoleOutOnly & mom){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 mtmp = mom.mass;
            F64vec ptmp = mom.pos - this->pos;
            F64 cx = mtmp * ptmp.x;
            F64 cy = mtmp * ptmp.y;
            F64 cz = mtmp * ptmp.z;
            this->quad.xx += cx * ptmp.x + mom.quad.xx;
            this->quad.yy += cy * ptmp.y + mom.quad.yy;
            this->quad.zz += cz * ptmp.z + mom.quad.zz;
            this->quad.xy += cx * ptmp.y + mom.quad.xy;
            this->quad.xz += cx * ptmp.z + mom.quad.xz;
            this->quad.yz += cy * ptmp.z + mom.quad.yz;
#else
            // under construction
#endif
        }
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
        }
    };

    class MomentQuadrupoleInAndOut{
    public:
        FSP mass;
        FSPvec pos;
        FSPmat quad;
        F64ort vertex_out_;
        F64ort vertex_in_;
        MomentQuadrupoleInAndOut(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        MomentQuadrupoleInAndOut(const FSP m, const FSPvec & p, const FSPmat & q){ 
            mass = m;
            pos = p;
            quad = q;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return mass;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->mass += epj.getCharge();
            this->pos += epj.getCharge() * epj.getPos();
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
            (this->vertex_in_).merge(epj.getPos());
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = epj.getCharge();
            F64vec ptmp = epj.getPos() - this->pos;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quad.xx += cx * ptmp.x;
            this->quad.yy += cy * ptmp.y;
            this->quad.zz += cz * ptmp.z;
            this->quad.xy += cx * ptmp.y;
            this->quad.xz += cx * ptmp.z;
            this->quad.yz += cy * ptmp.z;
#else
            // under construction
#endif
        }
        void set(){
            pos = pos / mass;
        }
        void accumulate(const MomentQuadrupoleInAndOut & mom){
            this->mass += mom.mass;
            this->pos += mom.mass * mom.pos;
            (this->vertex_out_).merge(mom.vertex_out_);
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        void accumulate2(const MomentQuadrupoleInAndOut & mom){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 mtmp = mom.mass;
            F64vec ptmp = mom.pos - this->pos;
            F64 cx = mtmp * ptmp.x;
            F64 cy = mtmp * ptmp.y;
            F64 cz = mtmp * ptmp.z;
            this->quad.xx += cx * ptmp.x + mom.quad.xx;
            this->quad.yy += cy * ptmp.y + mom.quad.yy;
            this->quad.zz += cz * ptmp.z + mom.quad.zz;
            this->quad.xy += cx * ptmp.y + mom.quad.xy;
            this->quad.xz += cx * ptmp.z + mom.quad.xz;
            this->quad.yz += cy * ptmp.z + mom.quad.yz;
#else
            // under construction
#endif
        }
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };

    class MomentMonopoleGeometricCenter{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        MomentMonopoleGeometricCenter(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
        }
        MomentMonopoleGeometricCenter(const FSP c, const FSPvec & p, const SSP n){
            n_ptcl = n;
            charge = c;
            pos = p;
        }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            charge += epj.getCharge();
            pos += epj.getPos();
            n_ptcl++;
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){}
        void set(){
            pos = pos / n_ptcl;
        }
        void accumulate(const MomentMonopoleGeometricCenter & mom){
            charge += mom.charge;
            pos += mom.n_ptcl * mom.pos;
            n_ptcl += mom.n_ptcl;
        }
        void accumulate2(const MomentMonopoleGeometricCenter & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
        }
    };

    class MomentMonopoleGeometricCenterInOnly{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        F64ort vertex_in_;
        MomentMonopoleGeometricCenterInOnly(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            vertex_in_.init();
        }
        MomentMonopoleGeometricCenterInOnly(const FSP c, const FSPvec & p, const SSP n){
            n_ptcl = n;
            charge = c;
            pos = p;
        }
        MomentMonopoleGeometricCenterInOnly(const FSP c, const FSPvec & p, const SSP n, const F64ort & in){
            n_ptcl = n;
            charge = c;
            pos = p;
            vertex_in_ = in;
        }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->charge += epj.getCharge();
            this->pos += epj.getPos();
            (this->n_ptcl)++;
            (this->vertex_in_).merge(epj.getPos());
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){}
        void set(){
            pos = pos / n_ptcl;
        }
        void accumulate(const MomentMonopoleGeometricCenterInOnly & mom){
            this->charge += mom.charge;
            this->pos += mom.n_ptcl * mom.pos;
            this->n_ptcl += mom.n_ptcl;
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        void accumulate2(const MomentMonopoleGeometricCenterInOnly & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };

    class MomentMonopoleGeometricCenterOutOnly{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        F64ort vertex_out_;
        MomentMonopoleGeometricCenterOutOnly(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            vertex_out_.init();
        }
        MomentMonopoleGeometricCenterOutOnly(const FSP c, const FSPvec & p, const SSP n){
            n_ptcl = n;
            charge = c;
            pos = p;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            vertex_out_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->charge += epj.getCharge();
            this->pos += epj.getPos();
            (this->n_ptcl)++;
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){}
        void set(){
            pos = pos / n_ptcl;
        }
        void accumulate(const MomentMonopoleGeometricCenterOutOnly & mom){
            this->charge += mom.charge;
            this->pos += mom.n_ptcl * mom.pos;
            this->n_ptcl += mom.n_ptcl;
            (this->vertex_out_).merge(mom.vertex_out_);
        }
        void accumulate2(const MomentMonopoleGeometricCenterOutOnly & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
        }
    };

    class MomentMonopoleGeometricCenterInAndOut{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        F64ort vertex_out_;
        F64ort vertex_in_;
        MomentMonopoleGeometricCenterInAndOut(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        MomentMonopoleGeometricCenterInAndOut(const FSP c, const FSPvec & p, const SSP n){
            n_ptcl = n;
            charge = c;
            pos = p;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->charge += epj.getCharge();
            this->pos += epj.getPos();
            (this->n_ptcl)++;
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
            (this->vertex_in_).merge(epj.getPos());
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){}
        void set(){
            pos = pos / n_ptcl;
        }
        void accumulate(const MomentMonopoleGeometricCenterInAndOut & mom){
            this->charge += mom.charge;
            this->pos += mom.n_ptcl * mom.pos;
            this->n_ptcl += mom.n_ptcl;
            (this->vertex_out_).merge(mom.vertex_out_);
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        void accumulate2(const MomentMonopoleGeometricCenterInAndOut & mom){}
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };

    class MomentDipoleGeometricCenter{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        MomentDipoleGeometricCenter(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
        }
        MomentDipoleGeometricCenter(const FSP c, const FSPvec & p, 
                                    const SSP n, const FSPvec & di){
            charge = c;
            pos = p;
            n_ptcl = n;
            dipole = di;
        }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = dipole = 0.0;
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            charge += epj.getCharge();
            pos += epj.getPos();
            n_ptcl++;
        }
        void accumulate(const MomentDipoleGeometricCenter & mom){
            charge += mom.charge;
            pos += mom.n_ptcl * mom.pos;
            n_ptcl += mom.n_ptcl;
        }
        void set(){
            pos = pos / n_ptcl;
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
            this->dipole += epj.getCharge() * (epj.getPos() - this->pos);
        }
        void accumulate2(const MomentDipoleGeometricCenter & mom){
            //dipole += mom.charge * (mom.pos - this->pos);
            this->dipole += mom.charge * (mom.pos - this->pos) + mom.dipole;
        }
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
            fout<<"dipole="<<dipole<<std::endl;
        }
    };

    class MomentDipoleGeometricCenterInOnly{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        F64ort vertex_in_;
        MomentDipoleGeometricCenterInOnly(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
            vertex_in_.init();
        }
        MomentDipoleGeometricCenterInOnly(const FSP c, const FSPvec & p, 
                                          const SSP n, const FSPvec & di){
            charge = c;
            pos = p;
            n_ptcl = n;
            dipole = di;
        }
        MomentDipoleGeometricCenterInOnly(const FSP c, const FSPvec & p, 
                                          const SSP n, const FSPvec & di,
                                          const F64ort & in){
            charge = c;
            pos = p;
            n_ptcl = n;
            dipole = di;
            vertex_in_ = in;
        }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = dipole = 0.0;
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->charge += epj.getCharge();
            this->pos += epj.getPos();
            (this->n_ptcl)++;
            (this->vertex_in_).merge(epj.getPos());
        }
        void accumulate(const MomentDipoleGeometricCenterInOnly & mom){
            this->charge += mom.charge;
            this->pos += mom.n_ptcl * mom.pos;
            this->n_ptcl += mom.n_ptcl;
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        void set(){
            pos = pos / n_ptcl;
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
            this->dipole += epj.getCharge() * (epj.getPos() - this->pos);
        }
        void accumulate2(const MomentDipoleGeometricCenterInOnly & mom){
            //dipole += mom.charge * (mom.pos - this->pos);
            this->dipole += mom.charge * (mom.pos - this->pos) + mom.dipole;
        }
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
            fout<<"dipole="<<dipole<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };

    class MomentDipoleGeometricCenterOutOnly{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        F64ort vertex_out_;
        MomentDipoleGeometricCenterOutOnly(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
            vertex_out_.init();
        }
        MomentDipoleGeometricCenterOutOnly(const FSP c, const FSPvec & p, 
                                           const SSP n, const FSPvec & di){
            charge = c;
            pos = p;
            n_ptcl = n;
            dipole = di;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = dipole = 0.0;
            vertex_out_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->charge += epj.getCharge();
            this->pos += epj.getPos();
            (this->n_ptcl)++;
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
        }
        void accumulate(const MomentDipoleGeometricCenterOutOnly & mom){
            this->charge += mom.charge;
            this->pos += mom.n_ptcl * mom.pos;
            this->n_ptcl += mom.n_ptcl;
            (this->vertex_out_).merge(mom.vertex_out_);
        }
        void set(){
            pos = pos / n_ptcl;
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
            this->dipole += epj.getCharge() * (epj.getPos() - this->pos);
        }
        void accumulate2(const MomentDipoleGeometricCenterOutOnly & mom){
            //dipole += mom.charge * (mom.pos - this->pos);
            this->dipole += mom.charge * (mom.pos - this->pos) + mom.dipole;
        }
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
            fout<<"dipole="<<dipole<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
        }
    };

    class MomentDipoleGeometricCenterInAndOut{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        F64ort vertex_out_;
        F64ort vertex_in_;
        MomentDipoleGeometricCenterInAndOut(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        MomentDipoleGeometricCenterInAndOut(const FSP c, const FSPvec & p, 
                                            const SSP n, const FSPvec & di){
            charge = c;
            pos = p;
            n_ptcl = n;
            dipole = di;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = dipole = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->charge += epj.getCharge();
            this->pos += epj.getPos();
            (this->n_ptcl)++;
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
            (this->vertex_in_).merge(epj.getPos());
        }
        void accumulate(const MomentDipoleGeometricCenterInAndOut & mom){
            this->charge += mom.charge;
            this->pos += mom.n_ptcl * mom.pos;
            this->n_ptcl += mom.n_ptcl;
            (this->vertex_out_).merge(mom.vertex_out_);
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        void set(){
            pos = pos / n_ptcl;
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
            this->dipole += epj.getCharge() * (epj.getPos() - this->pos);
        }
        void accumulate2(const MomentDipoleGeometricCenterInAndOut & mom){
            //dipole += mom.charge * (mom.pos - this->pos);
            this->dipole += mom.charge * (mom.pos - this->pos) + mom.dipole;
        }
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
            fout<<"dipole="<<dipole<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };

    class MomentQuadrupoleGeometricCenter{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        FSPmat quadrupole;
        MomentQuadrupoleGeometricCenter(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
            quadrupole = 0.0;
        }
        MomentQuadrupoleGeometricCenter(const FSP c, const FSPvec & p, 
                                        const SSP n, const FSPvec & di,
                                        const FSPmat & q){
            charge = c;
            pos = p;
            n_ptcl = n;
            dipole = di;
            quadrupole = q;
        }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = dipole = 0.0;
            quadrupole = 0.0;
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            charge += epj.getCharge();
            pos += epj.getPos();
            n_ptcl++;
        }
        void set(){
            pos = pos / n_ptcl;
        }
        void accumulate(const MomentQuadrupoleGeometricCenter & mom){
            charge += mom.charge;
            pos += mom.n_ptcl * mom.pos;
            n_ptcl += mom.n_ptcl;
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = epj.getCharge();
            F64vec ptmp = epj.getPos() - this->pos;
            this->dipole += ctmp * ptmp;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quadrupole.xx += cx * ptmp.x;
            this->quadrupole.yy += cy * ptmp.y;
            this->quadrupole.zz += cz * ptmp.z;
            this->quadrupole.xy += cx * ptmp.y;
            this->quadrupole.xz += cx * ptmp.z;
            this->quadrupole.yz += cy * ptmp.z;
#else
	    // underconstruction
#endif
        }
        void accumulate2(const MomentQuadrupoleGeometricCenter & mom){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = mom.charge;
            F64vec ptmp = mom.pos - this->pos;
            this->dipole += ctmp * ptmp + mom.dipole;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quadrupole.xx += cx * ptmp.x + mom.quadrupole.xx;
            this->quadrupole.yy += cy * ptmp.y + mom.quadrupole.yy;
            this->quadrupole.zz += cz * ptmp.z + mom.quadrupole.zz;
            this->quadrupole.xy += cx * ptmp.y + mom.quadrupole.xy;
            this->quadrupole.xz += cx * ptmp.z + mom.quadrupole.xz;
            this->quadrupole.yz += cy * ptmp.z + mom.quadrupole.yz;
#else
	    // underconstruction
#endif
        }
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
            fout<<"dipole="<<dipole<<std::endl;
            fout<<"quadrupole="<<quadrupole<<std::endl;
        }
    };

    class MomentQuadrupoleGeometricCenterInOnly{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        FSPmat quadrupole;
        F64ort vertex_in_;
        MomentQuadrupoleGeometricCenterInOnly(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
            quadrupole = 0.0;
            vertex_in_.init();
        }
        MomentQuadrupoleGeometricCenterInOnly(const FSP c, const FSPvec & p, 
                                              const SSP n, const FSPvec & di,
                                              const FSPmat & q){
            charge = c;
            pos = p;
            n_ptcl = n;
            dipole = di;
            quadrupole = q;
        }
        MomentQuadrupoleGeometricCenterInOnly(const FSP c, const FSPvec & p, 
                                              const SSP n, const FSPvec & di,
                                              const FSPmat & q, const F64ort & in){
            charge = c;
            pos = p;
            n_ptcl = n;
            dipole = di;
            quadrupole = q;
            vertex_in_ = in;
        }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = dipole = 0.0;
            quadrupole = 0.0;
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->charge += epj.getCharge();
            this->pos += epj.getPos();
            this->n_ptcl++;
            (this->vertex_in_).merge(epj.getPos());
        }
        void set(){
            pos = pos / n_ptcl;
        }
        void accumulate(const MomentQuadrupoleGeometricCenterInOnly & mom){
            this->charge += mom.charge;
            this->pos += mom.n_ptcl * mom.pos;
            this->n_ptcl += mom.n_ptcl;
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = epj.getCharge();
            F64vec ptmp = epj.getPos() - this->pos;
            this->dipole += ctmp * ptmp;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quadrupole.xx += cx * ptmp.x;
            this->quadrupole.yy += cy * ptmp.y;
            this->quadrupole.zz += cz * ptmp.z;
            this->quadrupole.xy += cx * ptmp.y;
            this->quadrupole.xz += cx * ptmp.z;
            this->quadrupole.yz += cy * ptmp.z;
#else
            // underconstruction
#endif
        }
        void accumulate2(const MomentQuadrupoleGeometricCenterInOnly & mom){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = mom.charge;
            F64vec ptmp = mom.pos - this->pos;
            this->dipole += ctmp * ptmp + mom.dipole;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quadrupole.xx += cx * ptmp.x + mom.quadrupole.xx;
            this->quadrupole.yy += cy * ptmp.y + mom.quadrupole.yy;
            this->quadrupole.zz += cz * ptmp.z + mom.quadrupole.zz;
            this->quadrupole.xy += cx * ptmp.y + mom.quadrupole.xy;
            this->quadrupole.xz += cx * ptmp.z + mom.quadrupole.xz;
            this->quadrupole.yz += cy * ptmp.z + mom.quadrupole.yz;
#else
            // underconstruction
#endif
        }
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
            fout<<"dipole="<<dipole<<std::endl;
            fout<<"quadrupole="<<quadrupole<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };

    class MomentQuadrupoleGeometricCenterOutOnly{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        FSPmat quadrupole;
        F64ort vertex_out_;
        MomentQuadrupoleGeometricCenterOutOnly(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
            quadrupole = 0.0;
            vertex_out_.init();
        }
        MomentQuadrupoleGeometricCenterOutOnly(const FSP c, const FSPvec & p, 
                                               const SSP n, const FSPvec & di,
                                               const FSPmat & q){
            charge = c;
            pos = p;
            n_ptcl = n;
            dipole = di;
            quadrupole = q;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = dipole = 0.0;
            quadrupole = 0.0;
            vertex_out_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->charge += epj.getCharge();
            this->pos += epj.getPos();
            this->n_ptcl++;
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
        }
        void set(){
            pos = pos / n_ptcl;
        }
        void accumulate(const MomentQuadrupoleGeometricCenterOutOnly & mom){
            this->charge += mom.charge;
            this->pos += mom.n_ptcl * mom.pos;
            this->n_ptcl += mom.n_ptcl;
            (this->vertex_out_).merge(mom.vertex_out_);
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = epj.getCharge();
            F64vec ptmp = epj.getPos() - this->pos;
            this->dipole += ctmp * ptmp;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quadrupole.xx += cx * ptmp.x;
            this->quadrupole.yy += cy * ptmp.y;
            this->quadrupole.zz += cz * ptmp.z;
            this->quadrupole.xy += cx * ptmp.y;
            this->quadrupole.xz += cx * ptmp.z;
            this->quadrupole.yz += cy * ptmp.z;
#else
           // underconstruction
#endif
        }
        void accumulate2(const MomentQuadrupoleGeometricCenterOutOnly & mom){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = mom.charge;
            F64vec ptmp = mom.pos - this->pos;
            this->dipole += ctmp * ptmp + mom.dipole;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quadrupole.xx += cx * ptmp.x + mom.quadrupole.xx;
            this->quadrupole.yy += cy * ptmp.y + mom.quadrupole.yy;
            this->quadrupole.zz += cz * ptmp.z + mom.quadrupole.zz;
            this->quadrupole.xy += cx * ptmp.y + mom.quadrupole.xy;
            this->quadrupole.xz += cx * ptmp.z + mom.quadrupole.xz;
            this->quadrupole.yz += cy * ptmp.z + mom.quadrupole.yz;
#else
            // underconstruction
#endif
        }
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
            fout<<"dipole="<<dipole<<std::endl;
            fout<<"quadrupole="<<quadrupole<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
        }
    };

    class MomentQuadrupoleGeometricCenterInAndOut{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        FSPmat quadrupole;
        F64ort vertex_out_;
        F64ort vertex_in_;
        MomentQuadrupoleGeometricCenterInAndOut(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
            quadrupole = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        MomentQuadrupoleGeometricCenterInAndOut(const FSP c, const FSPvec & p, 
                                                const SSP n, const FSPvec & di,
                                                const FSPmat & q){
            charge = c;
            pos = p;
            n_ptcl = n;
            dipole = di;
            quadrupole = q;
        }
        F64ort getVertexOut() const { return vertex_out_; }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = dipole = 0.0;
            quadrupole = 0.0;
            vertex_out_.init();
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->charge += epj.getCharge();
            this->pos += epj.getPos();
            this->n_ptcl++;
            (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
            (this->vertex_in_).merge(epj.getPos());
        }
        void set(){
            pos = pos / n_ptcl;
        }
        void accumulate(const MomentQuadrupoleGeometricCenterInAndOut & mom){
            this->charge += mom.charge;
            this->pos += mom.n_ptcl * mom.pos;
            this->n_ptcl += mom.n_ptcl;
            (this->vertex_out_).merge(mom.vertex_out_);
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = epj.getCharge();
            F64vec ptmp = epj.getPos() - this->pos;
            this->dipole += ctmp * ptmp;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quadrupole.xx += cx * ptmp.x;
            this->quadrupole.yy += cy * ptmp.y;
            this->quadrupole.zz += cz * ptmp.z;
            this->quadrupole.xy += cx * ptmp.y;
            this->quadrupole.xz += cx * ptmp.z;
            this->quadrupole.yz += cy * ptmp.z;
#else
            // underconstruction
#endif
        }
        void accumulate2(const MomentQuadrupoleGeometricCenterInAndOut & mom){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            F64 ctmp = mom.charge;
            F64vec ptmp = mom.pos - this->pos;
            this->dipole += ctmp * ptmp + mom.dipole;
            F64 cx = ctmp * ptmp.x;
            F64 cy = ctmp * ptmp.y;
            F64 cz = ctmp * ptmp.z;
            this->quadrupole.xx += cx * ptmp.x + mom.quadrupole.xx;
            this->quadrupole.yy += cy * ptmp.y + mom.quadrupole.yy;
            this->quadrupole.zz += cz * ptmp.z + mom.quadrupole.zz;
            this->quadrupole.xy += cx * ptmp.y + mom.quadrupole.xy;
            this->quadrupole.xz += cx * ptmp.z + mom.quadrupole.xz;
            this->quadrupole.yz += cy * ptmp.z + mom.quadrupole.yz;
#else
            // underconstruction
#endif
        }
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
            fout<<"dipole="<<dipole<<std::endl;
            fout<<"quadrupole="<<quadrupole<<std::endl;
            fout<<"vertex_out.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out.high_="<<vertex_out_.high_<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };

#ifdef PARTICLE_SIMULATOR_PMM_EXPERIMENTAL_FEATURE
    template <int p>
    class MomentMultipoleGeometricCenterPMMM{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        MultipoleMoment0<p> mm;
        F64ort vertex_in_;
        MomentMultipoleGeometricCenterPMMM(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            mm.clear();
            vertex_in_.init();
        }
        template <int psrc>
        MomentMultipoleGeometricCenterPMMM(const FSP charge,
                                           const FSPvec & pos,
                                           const SSP n_ptcl,
                                           const MultipoleMoment0<psrc> & mm){
            assert(p == psrc);
            this->charge = charge;
            this->pos = pos;
            this->n_ptcl = n_ptcl;
            this->mm = mm;
        }
        template <int psrc>
        MomentMultipoleGeometricCenterPMMM(const FSP charge,
                                           const FSPvec & pos, 
                                           const SSP n_ptcl,
                                           const MultipoleMoment0<psrc> & mm,
                                           const F64ort & vertex_in){
            assert(p == psrc);
            this->charge = charge;
            this->pos = pos;
            this->n_ptcl = n_ptcl;
            this->mm = mm;
            this->vertex_in_ = vertex_in;
        }
        F64ort getVertexIn() const { return vertex_in_; }
        void init(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            mm.clear();
            vertex_in_.init();
        }
        FSPvec getPos() const {
            return pos;
        }
        FSP getCharge() const {
            return charge;
        }
        template<class Tepj>
        void accumulateAtLeaf(const Tepj & epj){
            this->charge += epj.getCharge();
            this->pos += epj.getPos();
            this->n_ptcl++;
            (this->vertex_in_).merge(epj.getPos());
        }
        void set(){
            pos = pos / n_ptcl;
        }
        void accumulate(const MomentMultipoleGeometricCenterPMMM<p> & mom){
            this->charge += mom.charge;
            this->pos += mom.n_ptcl * mom.pos;
            this->n_ptcl += mom.n_ptcl;
            (this->vertex_in_).merge(mom.vertex_in_);
        }
        template<class Tepj>
        void accumulateAtLeaf2(const Tepj & epj){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            const F64 charge_ptcl = epj.getCharge();
            const F64vec pos_ptcl = epj.getPos();
            mm.assign_particle(this->pos, pos_ptcl, charge_ptcl);
#else
            // underconstruction
#endif
        }
        template <int psrc>
        void accumulate2(const MomentMultipoleGeometricCenterPMMM<psrc> & mom){
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            const F64vec pos_src = mom.pos;
            mm.assign_from_MM(mom.mm, this->pos, pos_src);
#else
            // underconstruction
#endif
        }
        // for DEBUG 
        void dump(std::ostream & fout = std::cout) const {
            fout<<"charge="<<charge<<std::endl;
            fout<<"pos="<<pos<<std::endl;
            fout<<"n_ptcl="<<n_ptcl<<std::endl;
            fout<<"vertex_in.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in.high_="<<vertex_in_.high_<<std::endl;
        }
    };
#endif // PARTICLE_SIMULATOR_PMM_EXPERIMENTAL_FEATURE
    
    // for P^3M 
    typedef MomentMonopoleOutOnly    MomentMonopoleCutoff;
    typedef MomentMonopoleCutoff     MomentMonopolePeriodic;
    // [TODO] For what purpose is MomentMonopolePeriodic introduced ?

    // for P^3T
    typedef MomentMonopoleInAndOut MomentMonopoleScatter;
    typedef MomentMonopoleInAndOut MomentMonopoleSymmetry;
    typedef MomentQuadrupoleInAndOut MomentQuadrupoleScatter;
    typedef MomentQuadrupoleInAndOut MomentQuadrupoleSymmetry;
    // Note that there is a special class for P^T + PM.

    // for PM^3
    typedef MomentMonopoleInOnly                   MomentMonopolePMMM;
    typedef MomentQuadrupoleInOnly                 MomentQuadrupolePMMM;
    typedef MomentMonopoleGeometricCenterInOnly    MomentMonopoleGeometricCenterPMMM;
    typedef MomentDipoleGeometricCenterInOnly      MomentDipoleGeometricCenterPMMM;
    typedef MomentQuadrupoleGeometricCenterInOnly  MomentQuadrupoleGeometricCenterPMMM;
    
    class MomentSearchInAndOut{
    public:
        F64ort vertex_out_;
        F64ort vertex_in_;
        void init(){
            vertex_out_.init();
            vertex_in_.init();
        }
        F64ort getVertexOut() const { return vertex_out_; }
        F64ort getVertexIn() const { return vertex_in_; }
        template<class Tep>
        void accumulateAtLeaf(const Tep & ep){
            (this->vertex_out_).merge(ep.getPos(), ep.getRSearch());
            (this->vertex_in_).merge(ep.getPos());
        }
        template<class Tep>
        void accumulateAtLeaf(Tep & ep){
            (this->vertex_out_).merge(ep.getPos(), ep.getRSearch());
            (this->vertex_in_).merge(ep.getPos());
        }
        void set(){}
        void accumulate(const MomentSearchInAndOut & _mom){
            (this->vertex_out_).merge(_mom.vertex_out_);
            (this->vertex_in_).merge(_mom.vertex_in_);
        }
        template<class Tep>
        void accumulateAtLeaf2(Tep & ep){}
        void accumulate2(const MomentSearchInAndOut & _mom){}
        void dump(std::ostream & fout=std::cout) const {
            fout<<"vertex_out_.low_="<<vertex_out_.low_<<std::endl;
            fout<<"vertex_out_.high_="<<vertex_out_.high_<<std::endl;
            fout<<"vertex_in_.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in_.high_="<<vertex_in_.high_<<std::endl;
        }
    };

    class MomentSearchInOnly{
    public:
        //F32ort vertex_in_;
        F64ort vertex_in_;
        void init(){
            vertex_in_.init();
        }
        //F32ort getVertexIn() const { return vertex_in_; }
        F64ort getVertexIn() const { return vertex_in_; }
        template<class Tep>
        void accumulateAtLeaf(const Tep & ep){
            (this->vertex_in_).merge(ep.getPos());
        }
        void set(){}
        void accumulate(const MomentSearchInOnly & _mom){
            (this->vertex_in_).merge(_mom.vertex_in_);
        }
        template<class Tep>
        void accumulateAtLeaf2(const Tep & ep){}
        void accumulate2(const MomentSearchInOnly & _mom){}
        void dump(std::ostream & fout=std::cout) const {
            fout<<"vertex_in_.low_="<<vertex_in_.low_<<std::endl;
            fout<<"vertex_in_.high_="<<vertex_in_.high_<<std::endl;
        }
    };

    ///////////
    /// SPJ ///
    ///////////
    //  We implement the following classes in this order:
    //     SPJMonopole
    //     SPJMonopoleInOnly
    //     SPJMonopoleOutOnly
    //     SPJMonopoleInAndOut
    //     SPJMonopolePMMM
    //     SPJMonopoleCutoffScatter
    //     SPJQuadrupole
    //     SPJQuadrupoleInOnly
    //     SPJQuadrupoleOutOnly
    //     SPJQuadrupoleInAndOut
    //     SPJQuadrupolePMMM
    //     SPJMonopoleGeometricCenter
    //     SPJMonopoleGeometricCenterInOnly
    //     SPJMonopoleGeometricCenterOutOnly
    //     SPJMonopoleGeometricCenterInAndOut
    //     SPJMonopoleGeometricCenterPMMM
    //     SPJDipoleGeometricCenter
    //     SPJDipoleGeometricCenterInOnly
    //     SPJDipoleGeometricCenterOutOnly
    //     SPJDipoleGeometricCenterInAndOut
    //     SPJDipoleGeometricCenterPMMM
    //     SPJQuadrupoleGeometricCenter
    //     SPJQuadrupoleGeometricCenterInOnly
    //     SPJQuadrupoleGeometricCenterOutOnly
    //     SPJQuadrupoleGeometricCenterInAndOut
    //     SPJQuadrupoleGeometricCenterPMMM
    //     SPJMultipoleGeometricCenterPMMM
    // After that, there are several ``typedef" to give another names.
    // The only difference between classes is the definition of member
    // function of convertToMoment().
    // [TODO]
    //     Again, the implementation becomes short
    //     if we use the inheritance feature in C++.
    class SPJMonopole{
    public:
        FSP mass;
        FSPvec pos;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            FSP mass = mom.mass;
            FSPvec pos = mom.pos;
            this->mass = mass;
            this->pos = pos;
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return mass;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopole convertToMoment() const {
            return MomentMonopole(mass, pos);
        }
    };

    class SPJMonopoleInOnly{
    public:
        FSP mass;
        FSPvec pos;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            FSP mass = mom.mass;
            FSPvec pos = mom.pos;
            this->mass = mass;
            this->pos = pos;
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return mass;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopoleInOnly convertToMoment() const {
            return MomentMonopoleInOnly(mass, pos);
        }
    };

    class SPJMonopoleOutOnly{
    public:
        FSP mass;
        FSPvec pos;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            FSP mass = mom.mass;
            FSPvec pos = mom.pos;
            this->mass = mass;
            this->pos = pos;
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return mass;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopoleOutOnly convertToMoment() const {
            return MomentMonopoleOutOnly(mass, pos);
        }
    };

    class SPJMonopoleInAndOut{
    public:
        FSP mass;
        FSPvec pos;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            FSP mass = mom.mass;
            FSPvec pos = mom.pos;
            this->mass = mass;
            this->pos = pos;
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return mass;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopoleInAndOut convertToMoment() const {
            return MomentMonopoleInAndOut(mass, pos);
        }
    };

    class SPJMonopolePMMM{ // for PMMM
    public:
        FSP mass;
        FSPvec pos;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            FSP mass = mom.mass;
            FSPvec pos = mom.pos;
            this->mass = mass;
            this->pos = pos;
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return mass;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopoleInOnly convertToMoment() const {
            return MomentMonopoleInOnly(mass, pos, F64ort(pos, pos));
        }
    };

    class SPJMonopoleCutoffScatter{// for P^3T + PM
    public:
        FSP mass;
        FSPvec pos;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            mass = mom.mass;
            pos = mom.pos;
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return mass;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopoleCutoffScatter convertToMoment() const {
            return MomentMonopoleCutoffScatter(mass, pos);
        }
    };

    class SPJQuadrupole{
    public:
        FSP mass;
        FSPvec pos;
        FSPmat quad;
        FSP getCharge() const {
            return mass;
        }
        FSPmat getQuadrupole() const {
            return quad;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        void copyFromMoment(const MomentQuadrupole & mom){
            mass = mom.mass;
            pos = mom.pos;
            quad = mom.quad;
        }
        MomentQuadrupole convertToMoment() const {
            return MomentQuadrupole(mass, pos, quad);
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
        }
    };

    class SPJQuadrupoleInOnly{
    public:
        FSP mass;
        FSPvec pos;
        FSPmat quad;
        FSP getCharge() const {
            return mass;
        }
        FSPmat getQuadrupole() const {
            return quad;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            FSP mass = mom.mass;
            FSPvec pos = mom.pos;
            FSPmat quad = mom.quad;
            this->mass = mass;
            this->pos = pos;
            this->quad = quad;
        }
        MomentQuadrupoleInOnly convertToMoment() const {
            return MomentQuadrupoleInOnly(mass, pos, quad);
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
        }
        void dump(std::ostream & fout=std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
        }
    };

    class SPJQuadrupoleOutOnly{
    public:
        FSP mass;
        FSPvec pos;
        FSPmat quad;
        FSP getCharge() const {
            return mass;
        }
        FSPmat getQuadrupole() const {
            return quad;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            FSP mass = mom.mass;
            FSPvec pos = mom.pos;
            FSPmat quad = mom.quad;
            this->mass = mass;
            this->pos = pos;
            this->quad = quad;
        }
        MomentQuadrupoleOutOnly convertToMoment() const {
            return MomentQuadrupoleOutOnly(mass, pos, quad);
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
        }
        void dump(std::ostream & fout=std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
        }
    };

    class SPJQuadrupoleInAndOut{
    public:
        FSP mass;
        FSPvec pos;
        FSPmat quad;
        FSP getCharge() const {
            return mass;
        }
        FSPmat getQuadrupole() const {
            return quad;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            FSP mass = mom.mass;
            FSPvec pos = mom.pos;
            FSPmat quad = mom.quad;
            this->mass = mass;
            this->pos = pos;
            this->quad = quad;
        }
        MomentQuadrupoleInAndOut convertToMoment() const {
            return MomentQuadrupoleInAndOut(mass, pos, quad);
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
        }
        void dump(std::ostream & fout=std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
        }
    };

    class SPJQuadrupolePMMM{
    public:
        FSP mass;
        FSPvec pos;
        FSPmat quad;
        FSP getCharge() const {
            return mass;
        }
        FSPmat getQuadrupole() const {
            return quad;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            FSP mass = mom.mass;
            FSPvec pos = mom.pos;
            FSPmat quad = mom.quad;
            this->mass = mass;
            this->pos = pos;
            this->quad = quad;
        }
        MomentQuadrupoleInOnly convertToMoment() const {
            return MomentQuadrupoleInOnly(mass, pos, quad, F64ort(pos, pos));
        }
        void clear(){
            mass = 0.0;
            pos = 0.0;
            quad = 0.0;
        }
        void dump(std::ostream & fout=std::cout) const {
            fout<<"mass="<<mass<<std::endl;
            fout<<"pos="<<pos<<std::endl;
        }
    };

    class SPJMonopoleGeometricCenter{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopoleGeometricCenter convertToMoment() const {
            return MomentMonopoleGeometricCenter(charge, pos, n_ptcl);
        }
    };

    class SPJMonopoleGeometricCenterInOnly{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopoleGeometricCenterInOnly convertToMoment() const {
            return MomentMonopoleGeometricCenterInOnly(charge, pos, n_ptcl);
        }
    };

    class SPJMonopoleGeometricCenterOutOnly{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopoleGeometricCenterOutOnly convertToMoment() const {
            return MomentMonopoleGeometricCenterOutOnly(charge, pos, n_ptcl);
        }
    };

    class SPJMonopoleGeometricCenterInAndOut{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopoleGeometricCenterInAndOut convertToMoment() const {
            return MomentMonopoleGeometricCenterInAndOut(charge, pos, n_ptcl);
        }
    };

    class SPJMonopoleGeometricCenterPMMM{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMonopoleGeometricCenterInOnly convertToMoment() const {
            return MomentMonopoleGeometricCenterInOnly(charge, pos, n_ptcl, F64ort(pos, pos));
        }
    };

    class SPJDipoleGeometricCenter{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
            dipole = mom.dipole;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getDipole() const {
            return dipole;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentDipoleGeometricCenter convertToMoment() const {
            return MomentDipoleGeometricCenter(charge, pos, n_ptcl, dipole);
        }
    };

    class SPJDipoleGeometricCenterInOnly{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
            dipole = mom.dipole;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getDipole() const {
            return dipole;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentDipoleGeometricCenterInOnly convertToMoment() const {
            return MomentDipoleGeometricCenterInOnly(charge, pos, n_ptcl, dipole);
        }
    };

    class SPJDipoleGeometricCenterOutOnly{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
            dipole = mom.dipole;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getDipole() const {
            return dipole;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentDipoleGeometricCenterOutOnly convertToMoment() const {
            return MomentDipoleGeometricCenterOutOnly(charge, pos, n_ptcl, dipole);
        }
    };

    class SPJDipoleGeometricCenterInAndOut{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
            dipole = mom.dipole;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getDipole() const {
            return dipole;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentDipoleGeometricCenterInAndOut convertToMoment() const {
            return MomentDipoleGeometricCenterInAndOut(charge, pos, n_ptcl, dipole);
        }
    };

    class SPJDipoleGeometricCenterPMMM{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
            dipole = mom.dipole;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getDipole() const {
            return dipole;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentDipoleGeometricCenterInOnly convertToMoment() const {
            return MomentDipoleGeometricCenterInOnly(charge, pos, n_ptcl, dipole, F64ort(pos, pos));
        }
    };

    class SPJQuadrupoleGeometricCenter{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        FSPmat quadrupole;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
            dipole = mom.dipole;
            quadrupole = mom.quadrupole;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
            quadrupole = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getDipole() const {
            return dipole;
        }
        FSPmat getQuadrupole() const {
            return quadrupole;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentQuadrupoleGeometricCenter convertToMoment() const {
            return MomentQuadrupoleGeometricCenter(charge, pos, n_ptcl, dipole, quadrupole);
        }
    };

    class SPJQuadrupoleGeometricCenterInOnly{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        FSPmat quadrupole;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
            dipole = mom.dipole;
            quadrupole = mom.quadrupole;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
            quadrupole = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getDipole() const {
            return dipole;
        }
        FSPmat getQuadrupole() const {
            return quadrupole;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentQuadrupoleGeometricCenterInOnly convertToMoment() const {
            return MomentQuadrupoleGeometricCenterInOnly(charge, pos, n_ptcl, dipole, quadrupole);
        }
    };

    class SPJQuadrupoleGeometricCenterOutOnly{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        FSPmat quadrupole;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
            dipole = mom.dipole;
            quadrupole = mom.quadrupole;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
            quadrupole = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getDipole() const {
            return dipole;
        }
        FSPmat getQuadrupole() const {
            return quadrupole;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentQuadrupoleGeometricCenterOutOnly convertToMoment() const {
            return MomentQuadrupoleGeometricCenterOutOnly(charge, pos, n_ptcl, dipole, quadrupole);
        }
    };

    class SPJQuadrupoleGeometricCenterInAndOut{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        FSPmat quadrupole;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
            dipole = mom.dipole;
            quadrupole = mom.quadrupole;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
            quadrupole = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getDipole() const {
            return dipole;
        }
        FSPmat getQuadrupole() const {
            return quadrupole;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentQuadrupoleGeometricCenterInAndOut convertToMoment() const {
            return MomentQuadrupoleGeometricCenterInAndOut(charge, pos, n_ptcl, dipole, quadrupole);
        }
    };

    class SPJQuadrupoleGeometricCenterPMMM{
    public:
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        FSPvec dipole;
        FSPmat quadrupole;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
            dipole = mom.dipole;
            quadrupole = mom.quadrupole;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            dipole = 0.0;
            quadrupole = 0.0;
        }
        FSP getCharge() const {
            return charge;
        }
        FSPvec getDipole() const {
            return dipole;
        }
        FSPmat getQuadrupole() const {
            return quadrupole;
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentQuadrupoleGeometricCenterInOnly convertToMoment() const {
            return MomentQuadrupoleGeometricCenterInOnly(charge, pos, n_ptcl, dipole, quadrupole, F64ort(pos, pos));
        }
    };

#ifdef PARTICLE_SIMULATOR_PMM_EXPERIMENTAL_FEATURE
    template <int p>
    class SPJMultipoleGeometricCenterPMMM{
    public:
        enum {
            order = p,
        };
        SSP n_ptcl;
        FSP charge;
        FSPvec pos;
        MultipoleMoment0<p> mm;
        template<class Tmom>
        void copyFromMoment(const Tmom & mom){
            assert(p == mom.mm.order);
            n_ptcl = mom.n_ptcl;
            charge = mom.charge;
            pos = mom.pos;
            mm = mom.mm;
        }
        void clear(){
            n_ptcl = 0;
            charge = 0.0;
            pos = 0.0;
            mm.clear();
        }
        FSP getCharge() const {
            return charge;
        }
        constexpr static S32 getMultipoleOrder() { return order; }
        template <class Tbuf>
        void getMultipole(S32 buflen, Tbuf *buf) const {
            for (S32 i=0; i < buflen; i++) buf[i] = mm.buf[i];
        }
        FSPvec getPos() const {
            return pos;
        }
        void setPos(const FSPvec & pos_new) {
            pos = pos_new;
        }
        MomentMultipoleGeometricCenterPMMM<p> convertToMoment() const {
            return MomentMultipoleGeometricCenterPMMM<p>(charge, pos, n_ptcl, mm, F64ort(pos, pos));
        }
    };
#endif // PARTICLE_SIMULATOR_PMM_EXPERIMENTAL_FEATURE

    // for P^3M
    typedef SPJMonopoleOutOnly SPJMonopoleCutoff;
    typedef SPJMonopoleCutoff  SPJMonopolePeriodic;
    // [TODO] For what purpose is SPJMonopolePeriodic introduced ?

    // for P^3T
    typedef SPJMonopoleInAndOut SPJMonopoleScatter;
    typedef SPJMonopoleInAndOut SPJMonopoleSymmetry;
    typedef SPJQuadrupoleInAndOut SPJQuadrupoleScatter;
    typedef SPJQuadrupoleInAndOut SPJQuadrupoleSymmetry;

    class SuperParticleBase{
    public:
        void clear(){}
    };

    class EssentialParticleBase{
    public:
        F64vec pos;
        F64    r_search;
        F64vec getPos() const {
            return pos;
        }
        F64 getRSearch() const {
            return r_search;
        }
        void setPos(const F64vec & _pos){
            pos = _pos;
        }
        void clear(){}
    };

}
