#pragma once
#include <unistd.h>
#include <cmath>


namespace ParticleSimulator{
    F64 TREE_Z_COORD_OFFSET = 0.0;
}

namespace ParticleSimulator{

    template<class Tepi, class Tipg, class Tepj, class Tspj>
    void DumpInteractionList(const ReallocatableArray<Tepi> & epi,
                             const ReallocatableArray<Tipg> & ipg,
                             const ReallocatableArray<Tepj> & epj,
                             const ReallocatableArray<int> & adr_epj,
                             const ReallocatableArray<int> & n_epj,
                             const ReallocatableArray<int> & n_disp_epj,
                             const ReallocatableArray<Tspj> & spj,
                             const ReallocatableArray<int> & adr_spj,
                             const ReallocatableArray<int> & n_spj,
                             const ReallocatableArray<int> & n_disp_spj,
                             const int n_show = 20){
        F64 cm_mass = 0.0;
        F64vec cm_pos = 0.0;
        std::cerr<<"ipg.size()= "<<ipg.size()<<std::endl;
        std::cerr<<"n_epj.size()= "<<n_epj.size()<<std::endl;
        std::cerr<<"n_spj.size()= "<<n_spj.size()<<std::endl;
        for(S32 i=0; i<ipg.size(); i += (ipg.size()/n_show)+1){
            std::cerr<<"i= "<<i
                     <<" ipg[i].n_ptcl_= "<<ipg[i].n_ptcl_
                     <<" n_epj[i]= "<<n_epj[i]
                     <<" n_spj[i]= "<<n_spj[i]
                     <<std::endl;
            for(S32 j=ipg[i].adr_ptcl_; j<ipg[i+1].adr_ptcl_; j++){
                std::cerr<<"epi[i].pos= "<<epi[j].pos<<std::endl;
            }
            for(S32 j=n_disp_epj[i]; j<n_disp_epj[i+1]; j++){
                S32 adr = adr_epj[j];
                std::cerr<<"adr= "<<adr
                         <<" epj[adr].mass= "<<epj[adr].mass
                         <<" pos= "<<epj[adr].pos
                         <<std::endl;
                cm_mass += epj[adr].mass;
                cm_pos  += epj[adr].mass*epj[adr].pos;
            }
            for(S32 j=n_disp_spj[i]; j<n_disp_spj[i+1]; j++){
                S32 adr = adr_spj[j];
                std::cerr<<"adr= "<<adr
                         <<" spj[adr].mass= "<<spj[adr].mass
                         <<" pos= "<<spj[adr].pos
                         <<std::endl;
                cm_mass += spj[adr].mass;
                cm_pos  += spj[adr].mass*spj[adr].pos;
            }
        }
        std::cerr<<"cm_mass= "<<cm_mass
                 <<" cm_pos= "<<cm_pos / cm_mass
                 <<std::endl;
    }
    
    template<class Ttc>
    void DumpTreeCellInfo(const ReallocatableArray<Ttc> & tc, const int i){
        std::cerr<<"tc["<<i<<"].n_ptcl_= "<<tc[i].n_ptcl_<<std::endl;
        std::cerr<<"tc["<<i<<"].adr_ptcl_= "<<tc[i].adr_ptcl_<<std::endl;
        std::cerr<<"tc["<<i<<"].adr_tc_= "<<tc[i].adr_tc_<<std::endl;
        std::cerr<<"tc["<<i<<"].level_= "<<tc[i].level_<<std::endl;
        std::cerr<<"tc["<<i<<"].mom_.mass= "<<tc[i].mom_.mass<<std::endl;
        std::cerr<<"tc["<<i<<"].mom_.pos(cartesian)= "<<tc[i].mom_.pos<<std::endl;
#ifdef PHI_R_TREE
        std::cerr<<"tc["<<i<<"].mom_.pos(cylindrical)= "<<tc[i].mom_.pos_phi
                 <<"   "<<tc[i].mom_.pos_r<<"   "<<tc[i].mom_.pos.z<<std::endl;
#endif        
        std::cerr<<"tc["<<i<<"].mom_.vertex_out= "<<tc[i].mom_.vertex_out_<<std::endl;
        std::cerr<<"tc["<<i<<"].mom_.vertex_in= "<<tc[i].mom_.vertex_in_<<std::endl;
    }

    template<>
    void DumpTreeCellInfo(const ReallocatableArray<etcLM> & tc, const int i){
        std::cerr<<"tc["<<i<<"].n_ptcl_= "<<tc[i].n_ptcl_<<std::endl;
        std::cerr<<"tc["<<i<<"].adr_ptcl_= "<<tc[i].adr_ptcl_<<std::endl;
        std::cerr<<"tc["<<i<<"].adr_tc_= "<<tc[i].adr_tc_<<std::endl;
        std::cerr<<"tc["<<i<<"].level_= "<<tc[i].level_<<std::endl;
        std::cerr<<"tc["<<i<<"].mass= "<<tc[i].mass<<std::endl;
        std::cerr<<"tc["<<i<<"].pos= "<<tc[i].pos.x<<"   "<<tc[i].pos.y<<"   "<<tc[i].pos.z<<std::endl;
    }
    
    template<class Ttc>
    void TreeCellCheck(const ReallocatableArray<Ttc> & tc, const std::string & mes, const int check_rank=0){
        if(Comm::getRank() == check_rank){
            std::cerr<<mes<<std::endl;
            std::cerr<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
            std::cerr<<"tc.size()= "<<tc.size()<<std::endl;
            S32 n_cnt_tc_no_empty = 0;
            for(S32 i=0; i<tc.size(); i++){
                if(tc[i].n_ptcl_ > 0){
                    n_cnt_tc_no_empty++;
                    DumpTreeCellInfo(tc, i);
                }
            }
            std::cerr<<"n_cnt_tc_no_empty= "<<n_cnt_tc_no_empty<<std::endl;
        }
    }

    template<class Tptcl>
    void PtclCheck(const ReallocatableArray<Tptcl> & ptcl, const std::string & mes, const int n_show=100, const int check_rank=0){
        if(Comm::getRank() == check_rank){
            std::cerr<<mes<<std::endl;
            std::cerr<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
            std::cerr<<"ptcl.size()= "<<ptcl.size()<<std::endl;
            for(S32 i=0; i<ptcl.size(); i+=(ptcl.size()/n_show)+1){
                std::cerr<<"ptcl["<<i<<"].mass= "<<ptcl[i].mass
                         <<" pos= "<<ptcl[i].pos
                         <<std::endl;
            }
        }
    }

    template<class Tep, class Tsp, class Ttp>
    void PtclCheck2(const ReallocatableArray<Tep> & epj,
                    const ReallocatableArray<Tsp> & spj,
                    const ReallocatableArray<Ttp> & tp,
                    const std::string & mes,
                    const int n_show=100,
                    const int check_rank=0){
        if(Comm::getRank() == check_rank){
            std::cerr<<mes<<std::endl;
            std::cerr<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
            std::cerr<<"epj.size()= "<<epj.size()<<std::endl;
            std::cerr<<"spj.size()= "<<spj.size()<<std::endl;
            std::cerr<<"tp.size()= "<<tp.size()<<std::endl;
            for(S32 i=0; i<tp.size(); i+=(tp.size()/n_show)+1){
                U32 adr = tp[i].adr_ptcl_;
                if(GetMSB(adr)==0){
#ifndef REDUCE_MEMORY
                    adr = i;
#endif
                    std::cerr<<"i= "<<i
                             <<" epj["<<adr<<"].mass= "<<epj[adr].mass
                             <<" pos= "<<epj[adr].pos
                             <<std::endl;
                }
                else{
#ifdef REDUCE_MEMORY
                    adr = ClearMSB(adr);
#else
                    adr = i;
#endif
                    std::cerr<<"i= "<<i
                             <<" spj["<<adr<<"].mass= "<<spj[adr].mass
                             <<" pos= "<<spj[adr].pos
                             <<std::endl;
                }
            }
        }
    }
    
    class WtimeAbs{
    public:
        static F64 zero_point_;
        static F64 ex_ptcl_;
        static F64 set_ptcl_;
        static F64 set_root_cell_;
        static F64 morton_sort_lt_;
        static F64 morton_sort_fp_;
        static F64 link_cell_lt_;
        static F64 calc_moment_lt_;
        static F64 add_moment_lt_;
        static F64 exchange_let_;
        static F64 set_let_;
        static F64 morton_sort_gt_;
        static F64 link_cell_gt_;
        static F64 calc_moment_gt_;
        static F64 add_moment_gt_;
        static F64 make_ipg_;
        static F64 make_list_;
        static F64 balance_walk_;
        static F64 calc_force_;
        void dump(std::ofstream & fout){
            fout<<ex_ptcl_-zero_point_<<"   "
                <<set_ptcl_-zero_point_<<"   "
                <<set_root_cell_-zero_point_<<"   "
                <<morton_sort_lt_-zero_point_<<"   "
                <<std::endl;
        }
    };
    F64 WtimeAbs::zero_point_;
    F64 WtimeAbs::ex_ptcl_;
    F64 WtimeAbs::set_ptcl_;
    F64 WtimeAbs::set_root_cell_;
    F64 WtimeAbs::morton_sort_lt_;
    F64 WtimeAbs::morton_sort_fp_;
    F64 WtimeAbs::link_cell_lt_;
    F64 WtimeAbs::calc_moment_lt_;
    F64 WtimeAbs::add_moment_lt_;
    F64 WtimeAbs::exchange_let_;
    F64 WtimeAbs::set_let_;
    F64 WtimeAbs::morton_sort_gt_;
    F64 WtimeAbs::link_cell_gt_;
    F64 WtimeAbs::calc_moment_gt_;
    F64 WtimeAbs::add_moment_gt_;
    F64 WtimeAbs::make_ipg_;
    F64 WtimeAbs::make_list_;
    F64 WtimeAbs::balance_walk_;
    F64 WtimeAbs::calc_force_;
}



#ifdef SUNWAY
extern "C"{
    #include <athread.h>
    #include <cpe_func.h>
    void SLAVE_FUN(Preproc1_CopyEPJLocToEPJGlb)(void *);
    void SLAVE_FUN(Preproc2_CopyEPJLocToEPJGlb)(void *);
    void SLAVE_FUN(CopyEPJLocToEPJGlb)(void *);
    
    //void SLAVE_FUN(CopyEPJLocToEPJ)(void *);
    void SLAVE_FUN(CopyIndirect)(void *);
    void SLAVE_FUN(CopyIndirectInverse)(void *);
    //void SLAVE_FUN(CopyIndirectInverse2)(void *);
    void SLAVE_FUN(CopyIndirect2)(void *);
    void SLAVE_FUN(CopyDirect)(void *);
    void SLAVE_FUN(CopyStride)(void *);
    void SLAVE_FUN(GenMortonKey)(void *);
    void SLAVE_FUN(MakeListUsingTree)(void *);
    void SLAVE_FUN(balance_nwalk)(void *);
    void SLAVE_FUN(CopyEPISortedToEPJSortedLoc)(void *);
    void SLAVE_FUN(CopyTCToETC)(void *);
    void SLAVE_FUN(CopyEpjOrgToEpiSortedEpjSorted)(void *);
    void SLAVE_FUN(CopyFPToEPJ)(void *);
}
#endif

//#ifdef BW_CHECK_SAFETY
long cpe_cost(const int ki, const int kj, const int ks, const int nsat, const long EPJ_COST, const long SPJ_COST, const long SAT_COST, const long PLANET_COST, const long CONST_COST) {
  return ((kj+3)/4     * EPJ_COST + 
          (ks+3)/4     * SPJ_COST + 
          (nsat+3)/4   * SAT_COST + 
          PLANET_COST) * (ki+3)/4*4 + CONST_COST;
}
//#endif

#ifdef SUNWAY_PREFETCH
enum{
    L1ROFF = 64,
    L2ROFF = 256,
    L1WOFF = 64,
    L2WOFF = 256,
};
#endif

extern int N_SATELLITE;
extern Force_Kernel_Pars cpe_pars;
extern void samplesort(
		const int n,
		void *ptr_inp,
		void *ptr_work,
		void *ptr_out);

namespace ParticleSimulator{
    F64 wtime_make_key_local_tree = 0.0;
    F64 wtime_sort_local_tree = 0.0;
    F64 wtime_prefixsum_recorder = 0.0;
    F64 wtime_dispatch = 0.0;
    F64 wtime_retrieve = 0.0;
    F64 wtime_calc_force = 0.0;
    F64 wtime_t0, wtime_t1, wtime_t2, wtime_t3, wtime_t4,
        wtime_t5, wtime_t6, wtime_t7, wtime_t8;

    F64 wtime_ex_let_allgather_root_ex = 0.0;
    F64 wtime_ex_let_sd_2_allgather_bd = 0.0;
    F64 wtime_ex_let_sd_2_ex_n_d_d = 0.0;
    F64 wtime_ex_let_sd_2_ex_n_d_sd_1d_ring = 0.0;
    F64 wtime_ex_let_sd_2_ex_n_d_sd_sub = 0.0;
    F64 wtime_ex_let_sd_2_allgather_sd = 0.0;
    F64 wtime_ex_let_sd_2_ep_sp_among_sd = 0.0;
    F64 wtime_ex_let_sd_2_ep_sp_in_sd = 0.0;
    F64 wtime_ex_let_sd_2_make_list = 0.0;
    
    S32 n_walk_cpe[64],n_disp_walk[64];
    S32* adr_n_walk;
    
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
            WtimeAbs::set_ptcl_ = GetWtime();
#if 1
            setParticleLocalTreeImpl(psys, epj_org_, true); // new
#else
            setParticleLocalTreeImpl(psys, epi_org_, epj_org_, true); // new
#endif
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            if(Comm::getRank() == 0){
                std::cerr<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
            }
#endif
            
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK1 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            WtimeAbs::set_root_cell_ = GetWtime();
#ifdef PHI_R_TREE
            static const F64 PI = 4.0*atan(1.0);
            F64 len_root_tree = 2.0*PI*1.0001;
            //setRootCell(len_root_tree, F64vec(len_root_tree*0.5, 0.0, 0.0) ); // original
            //setRootCell(len_root_tree, F64vec(len_root_tree*0.5, 0.0, 0.0-TREE_Z_COORD_OFFSET) ); // original
            setRootCell(len_root_tree, F64vec(len_root_tree*0.5, 0.0, TREE_Z_COORD_OFFSET) ); // original
#else
            setRootCell(dinfo); // original
#endif

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            if(Comm::getRank() < 2){
                std::cerr<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
                std::cerr<<"length_= "<<length_<<std::endl;
                std::cerr<<"center_= "<<center_<<std::endl;
            }
            //exit(1);
#endif

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK2 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            WtimeAbs::morton_sort_lt_ = GetWtime();
    #if 1
            mortonSortLocalTreeOnlyImpl(epj_org_, epi_sorted_, epj_sorted_loc_); // new
    #else
            mortonSortLocalTreeOnlyImpl(epi_org_, epj_org_, epi_sorted_, epj_sorted_loc_); // new
    #endif

    #if 0
            U64 key_prev = 0;
            for(S32 i=0; i<epj_sorted_loc_.size(); i++){
                F64 MY_PHI_2 = 2.0*4.0*atan(1.0);
                F64 pos_x = epj_sorted_loc_[i].pos.x;
                F64 pos_y = epj_sorted_loc_[i].pos.y;
                F64 pos_z = epj_sorted_loc_[i].pos.z;
                F64 pos_phi = atan2(pos_y, pos_x);
                if(pos_phi < 0.0) pos_phi += MY_PHI_2;
                F64 pos_r = sqrt(pos_x*pos_x + pos_y*pos_y);
                U64 key = MortonKey::getKey( F64vec(pos_phi, pos_r, pos_z) );
                if(i>0 && Comm::getRank()==0 && i<10){
                    std::cerr<<"i= "<<i
                             <<" key_prev= "<<key_prev
                             <<" key= "<<key
                             <<" tp_loc_[i-1].key= "<<tp_loc_[i-1].key_
                             <<" tp_loc_[i].key= "<<tp_loc_[i].key_
                             <<std::endl;
                    std::cerr<<"pos_phi= "<<pos_phi<<" pos_r= "<<pos_r<<std::endl;
                }
                
                if(i>0 && Comm::getRank()==0){
                    if(key_prev > key){
                        std::cerr<<"i= "<<i
                                 <<" key_prev= "<<key_prev
                                 <<" key= "<<key
                                 <<" tp_loc_[i-1].key= "<<tp_loc_[i-1].key_
                                 <<" tp_loc_[i].key= "<<tp_loc_[i].key_
                                 <<std::endl;
                        std::cerr<<"pos_phi= "<<pos_phi<<" pos_r= "<<pos_r<<std::endl;
                    }
                    assert(key_prev <= key);
                }
                key_prev = key;
            }
            //exit(1);
    #endif
            
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK3 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
            //Finalize();
            //std::exit(0);
#endif
#ifdef DEBUG_PRING_CALC_FORCE_UNSAFE
            if(Comm::getRank() == 0){
                std::cerr<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
                for(S32 i=0; i<100; i++){
                    std::cerr<<std::oct<<"tp_loc_[i].key_= "<<tp_loc_[i].key_<<std::endl;
                    std::cerr<<std::dec;
                    std::cerr<<"epi_sorted_[i].pos= "<<epi_sorted_[i].pos<<std::endl;
                    //std::cerr<<"epj_sorted_loc_[i].pos= "<<epj_sorted_loc_[i].pos_phi<<"   "<<epj_sorted_loc_[i].pos_r<<"   "<<epj_sorted_loc_[i].pos<<std::endl;
                    std::cerr<<"epj_sorted_loc_[i].pos= "<<epj_sorted_loc_[i].pos<<std::endl;
                    std::cerr<<std::endl;
                }
            }
            //exit(1);
#endif
            WtimeAbs::morton_sort_fp_ = GetWtime();
        #ifdef DEBUG_SORT_LOCAL  // for consistency
            mortonSortFP<Tfp>(psys);    // new

            #if defined(DEBUG_PRING_CALC_FORCE_UNSAFE)
            if(Comm::getRank() == 0){
                std::cerr<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
                for(S32 i=0; i<10; i++){
                    std::cerr<<std::oct<<"tp_loc_[i].key_= "<<tp_loc_[i].key_<<std::endl;
                    std::cerr<<std::dec;
                    std::cerr<<"psys[i].pos= "<<psys[i].pos<<std::endl;
                    std::cerr<<"epi_sorted_[i].pos= "<<epi_sorted_[i].pos<<std::endl;
                    std::cerr<<"epj_sorted_loc_[i].pos= "<<epj_sorted_loc_[i].pos_phi<<"   "<<epj_sorted_loc_[i].pos_r<<"   "<<epj_sorted_loc_[i].pos<<std::endl;
                    std::cerr<<std::endl;
                }
            }
            //exit(1);
            #endif            
        #endif

            
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK4 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            WtimeAbs::link_cell_lt_ = GetWtime();
            linkCellLocalTreeOnly();   // original

            F64 wtime_offset = GetWtime();
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK5 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
            //Finalize();
            //std::exit(0);
#endif

#if defined(DEBUG_PRING_CALC_FORCE_UNSAFE)
            if(Comm::getRank() == 0){
                std::cerr<<"(Before) Comm::getRank()= "<<Comm::getRank()<<std::endl;
                std::cerr<<"tc_loc_[0].n_ptcl_= "<<tc_loc_[0].n_ptcl_<<std::endl;
                std::cerr<<"tc_loc_[0].mom_.mass= "<<tc_loc_[0].mom_.mass<<std::endl;
                std::cerr<<"tc_loc_[0].mom_.pos(cartesian)= "<<tc_loc_[0].mom_.pos<<std::endl;
    #ifdef PHI_R_TREE
                std::cerr<<"tc_loc_[0].mom_.pos(cylindrical)= "<<tc_loc_[0].mom_.pos_phi
                         <<"   "<<tc_loc_[0].mom_.pos_r<<"   "<<tc_loc_[0].mom_.pos.z<<std::endl;
    #endif                
                for(S32 i=8; i<8+8; i++){
                    if(tc_loc_[i].n_ptcl_ > 0){
                        std::cerr<<"tc_loc_["<<i<<"].n_ptcl_= "<<tc_loc_[i].n_ptcl_<<std::endl;
                        std::cerr<<"tc_loc_["<<i<<"].mom_.mass= "<<tc_loc_[i].mom_.mass<<std::endl;
                        std::cerr<<"tc_loc_["<<i<<"].mom_.pos(cartesian)= "<<tc_loc_[i].mom_.pos<<std::endl;
    #ifdef PHI_R_TREE
                        std::cerr<<"tc_loc_["<<i<<"].mom_.pos(cylindrical)= "<<tc_loc_[i].mom_.pos_phi
                                 <<"   "<<tc_loc_[i].mom_.pos_r<<"   "<<tc_loc_[i].mom_.pos.z<<std::endl;
    #endif
                    }
                }
            }
#endif
            
            WtimeAbs::calc_moment_lt_ = GetWtime();
            //calcMomentLocalTreeOnly(); // original
            CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(),
                       epj_sorted_loc_.getPointer(), lev_max_loc_, n_leaf_limit_); //new
            time_profile_.calc_moment_local_tree += GetWtime() - wtime_offset;

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK6 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            
#if defined(DEBUG_PRING_CALC_FORCE_UNSAFE)
            if(Comm::getRank() == 0){
                std::cerr<<"(After) Comm::getRank()= "<<Comm::getRank()<<std::endl;
                std::cerr<<"tc_loc_[0].n_ptcl_= "<<tc_loc_[0].n_ptcl_<<std::endl;
                std::cerr<<"tc_loc_[0].mom_.mass= "<<tc_loc_[0].mom_.mass<<std::endl;
                std::cerr<<"tc_loc_[0].mom_.pos(cartesian)= "<<tc_loc_[0].mom_.pos<<std::endl;
    #ifdef PHI_R_TREE
                std::cerr<<"tc_loc_[0].mom_.pos(cylindrical)= "<<tc_loc_[0].mom_.pos_phi
                         <<"   "<<tc_loc_[0].mom_.pos_r<<"   "<<tc_loc_[0].mom_.pos.z<<std::endl;
    #endif
                for(S32 i=8; i<8+8; i++){
                    if(tc_loc_[i].n_ptcl_ > 0){
                        std::cerr<<"tc_loc_["<<i<<"].n_ptcl_= "<<tc_loc_[i].n_ptcl_<<std::endl;
                        std::cerr<<"tc_loc_["<<i<<"].mom_.mass= "<<tc_loc_[i].mom_.mass<<std::endl;
                        std::cerr<<"tc_loc_["<<i<<"].mom_.pos(cartesian)= "<<tc_loc_[i].mom_.pos<<std::endl;
    #ifdef PHI_R_TREE
                        std::cerr<<"tc_loc_["<<i<<"].mom_.pos(cylindrical)= "<<tc_loc_[i].mom_.pos_phi
                                 <<"   "<<tc_loc_[i].mom_.pos_r<<"   "<<tc_loc_[i].mom_.pos.z<<std::endl;
    #endif
                    }
                }
            }
#endif
            
            

            WtimeAbs::add_moment_lt_ = GetWtime();
            addMomentAsSpLocalTreeImpl(typename TSM::force_type()); // original
            wtime_calc_force += MPI::Wtime() - wtime_offset_0;
            
            // MUST MODIFY TO USE ALLGATHER of monopole
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK7 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif

            WtimeAbs::exchange_let_ = GetWtime();
#ifdef USE_SUPER_DOMAIN
            //exchangeLocalEssentialTreeSuperDomain(dinfo, false);
            exchangeLocalEssentialTreeSuperDomain2(dinfo, false);
#else
            exchangeLocalEssentialTree3(dinfo); // original
#endif
            wtime_offset = MPI::Wtime(); // used for wtime_calc_force
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0){
                const S32 n_proc = Comm::getNumberOfProc();
                std::cerr<<"OK8 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
            }
#endif
#if defined(DEBUG_PRING_CALC_FORCE_UNSAFE)
            Comm::barrier();
            if(Comm::getRank()==0){
                std::cerr<<"n_ep_recv= "
                         <<comm_table_.ep_recv_.size()
                         <<std::endl;
                std::cerr<<"n_sp_recv= "
                         <<comm_table_.sp_recv_.size()
                         <<std::endl;
                for(S32 i=0; i<comm_table_.ep_recv_.size(); i += comm_table_.ep_recv_.size()/100 + 1){
                    std::cerr<<"i= "<<i
                             <<"ep_recv[i].mass= "<<comm_table_.ep_recv_[i].mass
                             <<" pos= "<<comm_table_.ep_recv_[i].pos
                             <<std::endl;
                }
                for(S32 i=0; i<comm_table_.sp_recv_.size(); i += comm_table_.sp_recv_.size()/100 + 1){
                    std::cerr<<"i= "<<i
                             <<"sp_recv[i].mass= "<<comm_table_.sp_recv_[i].mass
                             <<" pos= "<<comm_table_.sp_recv_[i].pos
                             <<std::endl;
                }
            }
            Comm::barrier();
            //exit(1);
#endif
            WtimeAbs::set_let_ = GetWtime();
            setLocalEssentialTreeToGlobalTree2(); // original

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0){
                std::cerr<<"OK9 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
            }
#endif
#if defined(DEBUG_PRING_CALC_FORCE_UNSAFE)
            Comm::barrier();
            if(Comm::getRank()==0){
                std::cerr<<"epj_org_.size()= "<<epj_org_.size()<<std::endl;
                std::cerr<<"spj_org_.size()= "<<spj_org_.size()<<std::endl;
                std::cerr<<"tp_glb_.size()= "<<tp_glb_.size()<<std::endl;
                for(S32 itmp=0; itmp<10; itmp++){
                    std::cerr<<"itmp= "<<itmp<<std::endl;
                    std::cerr<<"epj_org_[itmp].mass= "<<epj_org_[itmp].mass
                             <<" epj_org_[itmp].pos= "<<epj_org_[itmp].pos
                             <<std::endl;
    #ifdef REDUCE_MEMORY
                    S32 itmp2 = itmp;
    #else
                    S32 itmp2 = epj_org_.size()+itmp;
    #endif
                    std::cerr<<"spj_org_[itmp2].mass= "<<spj_org_[itmp2].mass
                             <<" spj_org_[itmp2].pos= "<<spj_org_[itmp2].pos
                             <<std::endl;
                    std::cerr<<"tp_glb_[itmp].key_= "<<tp_glb_[itmp].key_
                             <<" tp_glb_[itmp].adr_ptcl_= "<<tp_glb_[itmp].adr_ptcl_
                             <<std::endl;
                    S32 itmp3 = tp_glb_.size()-itmp-1;
                    std::cerr<<"tp_glb_[itmp3].key_= "<<tp_glb_[itmp3].key_
                             <<" tp_glb_[itmp3].adr_ptcl_= "<<tp_glb_[itmp3].adr_ptcl_
                             <<std::endl;
                }
                F64 mass_tmp = 0.0;
                F64vec cm_pos_tmp = 0.0;

                for(S32 itmp=0; itmp<epj_org_.size(); itmp++){
                    mass_tmp += epj_org_[itmp].mass;
                    cm_pos_tmp += epj_org_[itmp].mass*epj_org_[itmp].pos;
                }
    #ifdef REDUCE_MEMORY
                for(S32 itmp=0; itmp<spj_org_.size(); itmp++){
                    mass_tmp += spj_org_[itmp].mass;
                    cm_pos_tmp += spj_org_[itmp].mass*spj_org_[itmp].pos;
                }
    #else //REDUCE_MEMORY
                for(S32 itmp=epj_org_.size(); itmp<spj_org_.size(); itmp++){
                    mass_tmp += spj_org_[itmp].mass;
                    cm_pos_tmp += spj_org_[itmp].mass*spj_org_[itmp].pos;
                }
    #endif //REDUCE_MEMORY
                std::cerr<<"mass_tmp= "<<mass_tmp<<" cm_pos_tmp= "<<cm_pos_tmp/mass_tmp<<std::endl;
                F64 mass_tmp_2 = 0.0;
                F64vec cm_pos_tmp_2 = 0.0;
                for(S32 itmp=0; itmp<tp_glb_.size(); itmp++){

                    if( GetMSB(tp_glb_[itmp].adr_ptcl_) == 0){
    #ifdef REDUCE_MEMORY
                        S32 itmp2 = tp_glb_[itmp].adr_ptcl_;
    #else
                        S32 itmp2 = itmp;
    #endif
                        mass_tmp_2   += epj_org_[itmp2].mass;
                        cm_pos_tmp_2 += epj_org_[itmp2].mass*epj_org_[itmp2].pos;
                    }
                    else{
    #ifdef REDUCE_MEMORY
                        S32 itmp2 = ClearMSB(tp_glb_[itmp].adr_ptcl_);
    #else
                        S32 itmp2 = itmp;
    #endif
                        mass_tmp_2   += spj_org_[itmp2].mass;
                        cm_pos_tmp_2 += spj_org_[itmp2].mass*spj_org_[itmp2].pos;
                    }
                }
                std::cerr<<"mass_tmp_2= "<<mass_tmp_2<<" cm_pos_tmp_2= "<<cm_pos_tmp_2/mass_tmp_2<<std::endl;
                std::cerr<<std::endl;
            }
            Comm::barrier();
            //exit(1);
#endif
            WtimeAbs::morton_sort_gt_ = GetWtime();
            mortonSortGlobalTreeOnly3(); // new
            
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0){
                std::cerr<<"OK10 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
            }
#endif
#if defined(DEBUG_PRING_CALC_FORCE_UNSAFE)
            Comm::barrier();
            if(Comm::getRank()==0){
                std::cerr<<"epj_sorted_.size()= "<<epj_sorted_.size()<<std::endl;
                std::cerr<<"spj_sorted_.size()= "<<spj_sorted_.size()<<std::endl;
                std::cerr<<"tp_glb_.size()= "<<tp_glb_.size()<<std::endl;
                //for(S32 itmp=0; itmp<tp_glb_.size(); itmp += (tp_glb_.size()/2000)+1){
                for(S32 itmp=0; itmp<tp_glb_.size(); itmp += (tp_glb_.size()/100)+1){
                    std::cerr<<"itmp= "<<itmp<<std::endl;
                    U32 itmp2 = tp_glb_[itmp].adr_ptcl_;
                    if( GetMSB(itmp2) == 0){
    #ifdef REDUCE_MEMORY
                        itmp2 = itmp2;
    #else
                        itmp2 = itmp;
    #endif
                        std::cerr<<"itmp2= "<<itmp2
                                 <<" epj_sorted_[itmp2].mass= "<<epj_sorted_[itmp2].mass
                                 <<" epj_sorted_[itmp2].pos= "<<epj_sorted_[itmp2].pos
                                 <<std::endl;
                    }
                    else{
    #ifdef REDUCE_MEMORY                        
                        itmp2 = ClearMSB(itmp2);
    #else
                        itmp2 = itmp;
    #endif
                        std::cerr<<"itmp2= "<<itmp2
                                 <<" spj_sorted_[itmp2].mass= "<<spj_sorted_[itmp2].mass
                                 <<" spj_sorted_[itmp2].pos= "<<spj_sorted_[itmp2].pos
                                 <<std::endl;
                    }
                }
                F64 mass_tmp = 0.0;
                F64vec cm_pos_tmp = 0.0;
                for(S32 itmp=0; itmp<tp_glb_.size(); itmp++){
                    U32 itmp2 = tp_glb_[itmp].adr_ptcl_;
                    if( GetMSB(itmp2) == 0){
    #ifdef REDUCE_MEMORY
                        itmp2 = itmp2;
    #else
                        itmp2 = itmp;
    #endif
                        mass_tmp   += epj_sorted_[itmp2].mass;
                        cm_pos_tmp += epj_sorted_[itmp2].mass*epj_sorted_[itmp2].pos;
                    }
                    else{
    #ifdef REDUCE_MEMORY                        
                        itmp2 = ClearMSB(itmp2);
    #else
                        itmp2 = itmp;
    #endif
                        mass_tmp   += spj_sorted_[itmp2].mass;
                        cm_pos_tmp += spj_sorted_[itmp2].mass*spj_sorted_[itmp2].pos;
                    }
                }
                std::cerr<<"mass_tmp= "<<mass_tmp<<" cm_pos_tmp= "<<cm_pos_tmp/mass_tmp<<std::endl;
                std::cerr<<std::endl;
            }
            Comm::barrier();
            //exit(1);
#endif
            WtimeAbs::link_cell_gt_ = GetWtime();
            linkCellGlobalTreeOnly(); // original

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0){
                std::cerr<<"OK11 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
            }
            //Comm::barrier();
            //exit(1);
#endif
            WtimeAbs::calc_moment_gt_ = GetWtime();
            calcMomentGlobalTreeOnly(); // original

#if defined(DEBUG_PRING_CALC_FORCE_UNSAFE)
            Comm::barrier();
            if(Comm::getRank() == 0){
                //if(1){
                std::cerr<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
                std::cerr<<"tc_glb_.size()= "<<tc_glb_.size()<<std::endl;
                std::cerr<<"tc_glb_[0].n_ptcl_= "<<tc_glb_[0].n_ptcl_<<std::endl;
                std::cerr<<"tc_glb_[0].adr_ptcl_= "<<tc_glb_[0].adr_ptcl_<<std::endl;
                std::cerr<<"tc_glb_[0].adr_tc_= "<<tc_glb_[0].adr_tc_<<std::endl;
                std::cerr<<"tc_glb_[0].mom_.vertex_out= "<<tc_glb_[0].mom_.vertex_out_<<std::endl;
                std::cerr<<"tc_glb_[0].mom_.vertex_in= "<<tc_glb_[0].mom_.vertex_in_<<std::endl;
                std::cerr<<"tc_glb_[0].mom_.mass= "<<tc_glb_[0].mom_.mass<<std::endl;
                std::cerr<<"tc_glb_[0].mom_.pos(cartesian)= "<<tc_glb_[0].mom_.pos<<std::endl;
    #ifdef PHI_R_TREE
                std::cerr<<"tc_glb_[0].mom_.pos(cylindrical)= "<<tc_glb_[0].mom_.pos_phi
                         <<"   "<<tc_glb_[0].mom_.pos_r<<"   "<<tc_glb_[0].mom_.pos.z<<std::endl;
    #endif
                S32 n_cnt_tc_no_empty = 0;
                for(S32 i=0; i<tc_glb_.size(); i++){
                    if(tc_glb_[i].n_ptcl_ > 0){
                        std::cerr<<"tc_glb_["<<i<<"].n_ptcl_= "<<tc_glb_[i].n_ptcl_<<std::endl;
                        std::cerr<<"tc_glb_["<<i<<"].adr_ptcl_= "<<tc_glb_[i].adr_ptcl_<<std::endl;
                        std::cerr<<"tc_glb_["<<i<<"].adr_tc_= "<<tc_glb_[i].adr_tc_<<std::endl;
                        std::cerr<<"tc_glb_["<<i<<"].mom_.vertex_out= "<<tc_glb_[i].mom_.vertex_out_<<std::endl;
                        std::cerr<<"tc_glb_["<<i<<"].mom_.vertex_in= "<<tc_glb_[i].mom_.vertex_in_<<std::endl;
                        std::cerr<<"tc_glb_["<<i<<"].mom_.mass= "<<tc_glb_[i].mom_.mass<<std::endl;
                        std::cerr<<"tc_glb_["<<i<<"].mom_.pos(cartesian)= "<<tc_glb_[i].mom_.pos<<std::endl;
    #ifdef PHI_R_TREE
                        std::cerr<<"tc_glb_["<<i<<"].mom_.pos(cylindrical)= "<<tc_glb_[i].mom_.pos_phi
                                 <<"   "<<tc_glb_[i].mom_.pos_r<<"   "<<tc_glb_[i].mom_.pos.z<<std::endl;
    #endif
                        n_cnt_tc_no_empty++;
                        if(i==3703){
                            for(S32 j=0; j<tc_glb_[i].n_ptcl_; j++){
                                U32 adr_ptcl = tc_glb_[i].adr_ptcl_+j;
                                U32 adr_ptcl_2 = tp_glb_[adr_ptcl].adr_ptcl_;
                                U64 key = tp_glb_[adr_ptcl].key_;
                                std::cerr<<"adr_ptcl= "<<adr_ptcl<<std::endl;;
                                std::cerr<<"adr_ptcl_2= "<<adr_ptcl_2<<std::endl;;
                                std::cerr<<"key= "<<key<<std::endl;;
                                if(GetMSB(adr_ptcl_2)==0){
                                    std::cerr<<"adr_ptcl_2= "<<adr_ptcl_2
                                             <<" epj_sorted_[adr_ptcl_2].mass= "<<epj_sorted_[adr_ptcl_2].mass
                                             <<" pos= "<<epj_sorted_[adr_ptcl].pos
                                             <<std::endl;
                                }
                                else{
                                    std::cerr<<"adr_ptcl_2= "<<ClearMSB(adr_ptcl_2)
                                             <<" spj_sorted_[adr_ptcl_2].mass= "<<spj_sorted_[ClearMSB(adr_ptcl_2)].mass
                                             <<" pos= "<<spj_sorted_[ClearMSB(adr_ptcl_2)].pos
                                             <<std::endl;
                                }
                            }
                        }
                    }
                }
                std::cerr<<"n_cnt_tc_no_empty= "<<n_cnt_tc_no_empty<<std::endl;
            }
            Comm::barrier();
            //exit(1);
            S32 n_spj_sorted_prev = spj_sorted_.size(); // for debuging
#endif

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK12 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            WtimeAbs::add_moment_gt_ = GetWtime();
            addMomentAsSpGlobalTreeImpl(typename TSM::force_type()); // original

#if defined(DEBUG_PRING_CALC_FORCE_UNSAFE)
            Comm::barrier();
            if(Comm::getRank() == 0){
                std::cerr<<"spj_sorted_.size()= "<<spj_sorted_.size()<<std::endl;
                std::cerr<<"n_spj_sorted_prev= "<<n_spj_sorted_prev<<std::endl;
                for(S32 itmp=0; itmp<16; itmp++){
                    std::cerr<<"itmp= "<<itmp<<std::endl;
                    std::cerr<<"spj_sorted_[itmp].mass= "<<spj_sorted_[itmp].mass
                             <<" spj_sorted_[itmp].pos= "<<spj_sorted_[itmp].pos
                             <<std::endl;
                    std::cerr<<"spj_sorted_[n_spj_sorted_prev+itmp].mass= "<<spj_sorted_[n_spj_sorted_prev+itmp].mass
                             <<" spj_sorted_[n_spj_sorted_prev+itmp].pos= "<<spj_sorted_[n_spj_sorted_prev+itmp].pos
                             <<std::endl;
                }
            }
            Comm::barrier();
            //exit(1);
#endif
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0){
                std::cerr<<"OK13 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
            }
            Comm::barrier();
#endif
            WtimeAbs::make_ipg_ = GetWtime();

            makeIPGroup();

#ifdef DEBUG_PRING_CALC_FORCE_UNSAFE
            Comm::barrier();
            if(Comm::getRank() == 0){
                std::cerr<<"OK13.1 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
                for(S32 ig=0; ig<ipg_.size(); ig += (ipg_.size()/100)+1){
                    std::cerr<<"ig= "<<ig
                             <<" ipg_[ig].n_ptcl_= "<<ipg_[ig].n_ptcl_
                             <<" ipg_[ig].adr_ptcl_= "<<ipg_[ig].adr_ptcl_
                             <<" ipg_[ig].vertex_= "<<ipg_[ig].vertex_
                             <<std::endl;
                }
            }
            Comm::barrier();
#endif
            
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK14 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            WtimeAbs::make_list_ = GetWtime();
            makeAllInteractionListId3(true); // new

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK14.1 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            //exit(1);

#if 0
            //[DEBUG]
            //Comm::barrier();
            //if (Comm::getRank() == 0) {
            //    std::cout << "Checking the elapsed time for std::sort in MPE..." << std::endl;
            //    //S32 size = 1048576;
            //    //S32 size = 2*1048576;
            //    S32 size = 4*1048576;
            //    std::vector<S32> sample_data;
            //    MTTS mt;
            //    mt.init_genrand(0);
            //    for (S32 i=0; i<size; i++) {
            //        sample_data.push_back(mt.genrand_int31());
            //    }
            //    const double start_time = MPI_Wtime();
            //    std::sort(sample_data.begin(), sample_data.end());
            //    const double end_time = MPI_Wtime();
            //    std::cout << "std::sort completed!" << std::endl;
            //    std::cout << "time (std::sort) = " << end_time - start_time << "[s] " 
            //              << "(size = " << size << ")" << std::endl;
            //}
            athread_halt();
            Finalize();
            std::exit(0);
#endif
            
            //#if defined(DEBUG_PRING_CALC_FORCE_UNSAFE) && defined(PHI_R_TREE)
            WtimeAbs::balance_walk_ = GetWtime();
#if defined(DEBUG_PRING_CALC_FORCE_UNSAFE)
            Comm::barrier();
            if(Comm::getRank() == 0){
                for(S32 ig=0; ig<ipg_.size(); ig++){
                    F64 m_tmp = 0.0;
                    F64vec cm_pos_tmp  = 0.0;
                    F64vec cm_pos_tmp2 = 0.0;
                    S32 i_ep_head = n_disp_epj_recorder_for_force_[0][ig];
                    S32 i_ep_tail = i_ep_head+n_epj_recorder_for_force_[0][ig];
                    S32 i_sp_head = n_disp_spj_recorder_for_force_[0][ig];
                    S32 i_sp_tail = i_sp_head+n_spj_recorder_for_force_[0][ig];

                    for(S32 ip=i_ep_head; ip<i_ep_tail; ip++){
                        S32 adr = id_epj_recorder_for_force_[0][ip];
                        m_tmp += epj_sorted_[adr].mass;
                        cm_pos_tmp += epj_sorted_[adr].mass*epj_sorted_[adr].pos;
        #ifdef PHI_R_TREE
                        F64 pos_phi = atan2(epj_sorted_[adr].pos.y, epj_sorted_[adr].pos.x);
                        if(pos_phi < 0.0) pos_phi += 8.0*atan(1.0);
                        F64 pos_r = sqrt(epj_sorted_[adr].pos.x*epj_sorted_[adr].pos.x +
                                         epj_sorted_[adr].pos.y*epj_sorted_[adr].pos.y);
        #else
                        cm_pos_tmp2.x += epj_sorted_[adr].mass*epj_sorted_[adr].pos.x;
                        cm_pos_tmp2.y += epj_sorted_[adr].mass*epj_sorted_[adr].pos.y;
        #endif
                        cm_pos_tmp2.z += epj_sorted_[adr].mass*epj_sorted_[adr].pos.z;
                    }
                    for(S32 ip=i_sp_head; ip<i_sp_tail; ip++){
                        S32 adr = id_spj_recorder_for_force_[0][ip];
                        m_tmp += spj_sorted_[adr].mass;
                        cm_pos_tmp += spj_sorted_[adr].mass*spj_sorted_[adr].pos;
        #ifdef PHI_R_TREE
                        cm_pos_tmp2.x += spj_sorted_[adr].mass*spj_sorted_[adr].pos_phi;
                        cm_pos_tmp2.y += spj_sorted_[adr].mass*spj_sorted_[adr].pos_r;
        #else
                        cm_pos_tmp2.x += spj_sorted_[adr].mass*spj_sorted_[adr].pos.x;
                        cm_pos_tmp2.y += spj_sorted_[adr].mass*spj_sorted_[adr].pos.y;
        #endif
                        cm_pos_tmp2.z += spj_sorted_[adr].mass*spj_sorted_[adr].pos.z;
                    }
                    cm_pos_tmp /= m_tmp;
                    cm_pos_tmp2 /= m_tmp;
                    
                    if( (ig % ((ipg_.size())/100+1))==0 ){
                        std::cerr<<"ig= "<<ig
                                 <<" ipg_[ig].n_ptcl_= "<<ipg_[ig].n_ptcl_
                                 <<" ipg_[ig].adr_ptcl_= "<<ipg_[ig].adr_ptcl_
                                 <<" ipg_[ig].vertex_= "<<ipg_[ig].vertex_
                                 <<std::endl;
                        std::cerr<<"n_epj= "<<n_epj_recorder_for_force_[0][ig]
                                 <<" i_ep_head= "<<i_ep_head
                                 <<" i_ep_tail= "<<i_ep_tail
                                 <<std::endl;
                        std::cerr<<"n_spj= "<<n_spj_recorder_for_force_[0][ig]
                                 <<" i_sp_head= "<<i_sp_head
                                 <<" i_sp_tail= "<<i_sp_tail
                                 <<std::endl;
                        std::cerr<<" m_tmp= "<<m_tmp
                                 <<" cm_pos_tmp= "<<cm_pos_tmp
                                 <<" cm_pos_tmp2= "<<cm_pos_tmp2
                                 <<std::endl;
                    }
                }
            }
            Comm::barrier();
            //exit(1);
    #endif

            if(adr_n_walk!=NULL) delete [] adr_n_walk;
            adr_n_walk = new S32[ipg_.size()];
            n_epi_recorder_for_force_[0].resizeNoInitialize(ipg_.size());
            for(S32 i=0; i<ipg_.size(); i++){
              n_epi_recorder_for_force_[0][i] = ipg_[i].n_ptcl_;
            }

            cpe_pars.n_walk_cpe  = n_walk_cpe;
            cpe_pars.n_disp_walk = n_disp_walk;
            cpe_pars.adr_n_walk  = adr_n_walk;
    #ifdef BALANCE_NWALK
            cpe_pars.n_sat       = N_SATELLITE;
            cpe_pars.n_epi       = n_epi_recorder_for_force_[0].getPointer();
            cpe_pars.n_epj       = n_epj_recorder_for_force_[0].getPointer();
            cpe_pars.n_spj       = n_spj_recorder_for_force_[0].getPointer();
            cpe_pars.n_walk      = ipg_.size();

        #ifdef BALANCE_NWALK_OUT
            std::ofstream flog;
            std::stringstream fnum;
            std::string fname;
            fnum << std::setfill('0') << std::setw(5) << Comm::getRank();
            fname = "catch" + fnum.str() + ".txt";
            //flog.open(fname.c_str(),std::ios::trunc);
            flog.open(fname.c_str(),std::ios::app);
            flog<<"rank="<<Comm::getRank()<<" n="<<ipg_.size()
                     <<" nepi="<<n_epi_recorder_for_force_[0].size()
                     <<" nepj="<<n_epj_recorder_for_force_[0].size()
                     <<" nspj="<<n_spj_recorder_for_force_[0].size()
                     <<std::endl;
        #endif //BALANCE_NWALK_OUT
         
        #ifdef BW_CHECK_SAFETY  
            Comm::barrier();
            const double start_time = MPI_Wtime(); 
        #endif
            __real_athread_create(0,(void*)slave_balance_nwalk,(void*)&cpe_pars);
            athread_wait(0);
        #ifdef BW_CHECK_SAFETY
            Comm::barrier();
            const double end_time = MPI_Wtime();
            const int ipg_size_loc = ipg_.size();
            int ipg_size_max = 0;
            MPI_Allreduce(&ipg_size_loc, &ipg_size_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            if (Comm::getRank() == 0) {
                std::cout << "wtime (balance_nwalk) = " << end_time - start_time << " [s]"
                          << "(ipg_.size()[max.] = " << ipg_size_max << ")" << std::endl;
            }
        #endif
           
            
        #ifdef BALANCE_NWALK_OUT
            for(int i=0; i<64; i++){
                //flog<<"i= "<<i<<"  "<<adr_n_walk[i]<<std::endl;
                flog<<"i= "<<i<<"  "<<n_walk_cpe[i]<<std::endl;
            }
            flog.close();
        #endif
            
        #ifdef BW_CHECK_SAFETY
            long ncost[64]={0},nwalk_cpe=0,ni_cpe[64]={0},nj_cpe[64]={0},ns_cpe[64]={0},sum=0,ni_sum=0, nj_sum=0, ns_sum=0;
            for (int i=0; i<64; i++) {
              nwalk_cpe += n_walk_cpe[i];
              for (int j=0; j<n_walk_cpe[i]; j++) {
                int iw = adr_n_walk[n_disp_walk[i]+j];
                ni_cpe[i] += cpe_pars.n_epi[iw];
                nj_cpe[i] += cpe_pars.n_epj[iw];
                ns_cpe[i] += cpe_pars.n_spj[iw];
                ncost[i] += cpe_cost(cpe_pars.n_epi[iw],cpe_pars.n_epj[iw],cpe_pars.n_spj[iw],cpe_pars.n_sat,54,30,100,32,237404);
              }
              sum += ncost[i];
              ni_sum += ni_cpe[i];
              nj_sum += nj_cpe[i];
              ns_sum += ns_cpe[i];
            }
            long sum_mpe=0,ni_mpe=0,nj_mpe=0,ns_mpe=0;
            for(int k=0;k<ipg_.size();k++) {
              sum_mpe += cpe_cost(cpe_pars.n_epi[k], cpe_pars.n_epj[k], cpe_pars.n_spj[k], cpe_pars.n_sat, 54,30,100,32,237404);
              ni_mpe += cpe_pars.n_epi[k];
              nj_mpe += cpe_pars.n_epj[k];
              ns_mpe += cpe_pars.n_spj[k];
            }
            //assert(sum==sum_mpe);
            if(sum!=sum_mpe) {
              if(Comm::getRank()==0) {
              std::cerr<<"Ncost sum not match! CPE = "<<sum<<" MPE = "<<sum_mpe<<std::endl;
              std::cerr<<"Nw sum CPE="<<nwalk_cpe<<" MPE="<<ipg_.size()<<std::endl;
              std::cerr<<"NI sum CPE="<<ni_sum<<" MPE="<<ni_mpe<<std::endl;
              std::cerr<<"NJ sum CPE="<<nj_sum<<" MPE="<<nj_mpe<<std::endl;
              std::cerr<<"NS sum CPE="<<ns_sum<<" MPE="<<ns_mpe<<std::endl;

              fprintf(stderr,"Balance walk check on MPE:\n");
              fprintf(stderr,"CID         total");
              for(int k=0;k<64;k++) fprintf(stderr,"%12d",k);
              fprintf(stderr,"\n");

              int nwtot=0;
              for (int i=0; i<64; i++) nwtot += n_walk_cpe[i];
              fprintf(stderr,"NW_CPE%9d",nwtot);
              for(int k=0;k<64;k++) fprintf(stderr,"%12d",n_walk_cpe[k]);
              fprintf(stderr,"\nNW_DISP          ");
              for(int k=0;k<64;k++) fprintf(stderr,"%12d",n_disp_walk[k]);
              fprintf(stderr,"\nNW_ADR_FIRST     ");
              for(int k=0;k<64;k++) fprintf(stderr,"%12d",adr_n_walk[n_disp_walk[k]]);
              fprintf(stderr,"\n");
              fprintf(stderr,"SUM_I%12ld",ni_sum);
              for(int k=0;k<64;k++) fprintf(stderr,"%12d",ni_cpe[k]);
              fprintf(stderr,"\n");
              fprintf(stderr,"SUM%14ld",sum);
              for(int k=0;k<64;k++) fprintf(stderr,"%12d",ncost[k]);
              fprintf(stderr,"\n");
              }
//#ifdef BALANCE_NWALK_OUT
//              flog.open(fname.c_str(),std::ios::trunc);
//              for(int i=0; i<64; i++) flog<<"["<<i<<"]="<<ncost[i]<<" ";
//              flog.close();
//#endif
            }
        #endif //BW_CHECK_SAFETY
    #else //BALANCE_NWALK
        #ifdef BALANCE_NWALK_ON_MPE
            balanceNwalk();
        #else //BALANCE_NWALK_ON_MPE
            int nw_per_cpe = (ipg_.size()+63)/64;
            int nw_cutoff = ipg_.size()%64;
            if (nw_cutoff==0) nw_cutoff = 64;
            for(int i=0; i<64; i++) n_walk_cpe[i]  = (i<nw_cutoff)?nw_per_cpe:nw_per_cpe-1;
            n_disp_walk[0] = 0;
            for(int i=1; i<64; i++) n_disp_walk[i] = n_disp_walk[i-1] + n_walk_cpe[i-1];
            for(int i=0; i<ipg_.size(); i++) adr_n_walk[i] = i;
        #endif //BALANCE_NWALK_ON_MPE
    #endif //BALANCE_NWALK

#if 0
            //[DEBUG]
            Comm::barrier();
            if (Comm::getRank() == 0) {
                std::cout << "outputing adr_n_walk..." << std::endl;
                std::ofstream ofs;
                std::stringstream filenum;
                std::string filename;
                filenum << std::setfill('0') << std::setw(5) << Comm::getRank();
                filename = "adr_n_walk" + filenum.str() + ".txt";
                ofs.open(filename.c_str(),std::ios::trunc);
                for (S32 cpe_id=0; cpe_id<64; cpe_id++) {
                    S32 i_start = n_disp_walk[cpe_id];
                    S32 i_end   = i_start + n_walk_cpe[cpe_id];
                    S64 ncost_sum = 0;
                    ofs << "#####################" << std::endl;
                    ofs << "# CPE ID = " << cpe_id << std::endl;
                    ofs << "#####################" << std::endl;
                    for (S32 i=i_start; i<i_end; i++) {
                        S32 ig = adr_n_walk[i];
                        S32 ncost = cpe_cost(n_epi_recorder_for_force_[0][ig],
                                             n_epj_recorder_for_force_[0][ig],
                                             n_spj_recorder_for_force_[0][ig],
                                             N_SATELLITE,
                                             54,30,100,32,237404);
                        ofs << ig << " " << ncost << std::endl;
                        ncost_sum += ncost;
                    }
                    ofs << "(ncost_sum = " << ncost_sum << ")" << std::endl;
                }
                ofs.close(); 
                std::cout << "output of adr_n_walk is completed!" << std::endl;
            }
            athread_halt();
            Finalize();
            std::exit(0);
#endif

            /*
            //* For a test run
            athread_halt();
            Finalize();
            std::exit(0);
            */
            
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK15 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            WtimeAbs::calc_force_ = GetWtime();
            wtime_calc_force += MPI::Wtime() - wtime_offset;
            ret += calcForceUsingIdListMultiWalkIndex3(pfunc_dispatch, pfunc_retrieve, n_walk_limit, mw_info_, clear_force, reuse, true); // new

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK16 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
            //if (Comm::getRank() == 0) {
            //    std::cout << "etc_loc_.size()                             = " << etc_loc_.size() << std::endl;
            //    std::cout << "etc_loc_.getMemSize()                       = " << etc_loc_.getMemSize()*1.0e-9 << " [GB]" << std::endl;
            //    std::cout << "adr_tc_level_partition_loc_[lev_max_loc_+1] = " << adr_tc_level_partition_loc_[lev_max_loc_+1] << std::endl;
            //    std::cout << "Required mem. size (etc_loc_)               = " << sizeof(etcLM) * adr_tc_level_partition_loc_[lev_max_loc_+1] * 1.0e-9 << std::endl;
            //    std::cout << "etc_glb_.size()                             = " << etc_glb_.size() << std::endl;
            //    std::cout << "etc_glb_.getMemSize()                       = " << etc_glb_.getMemSize()*1.0e-9 << " [GB]" << std::endl;
            //    std::cout << "adr_tc_level_partition_glb_[lev_max_glb_+1] = " << adr_tc_level_partition_glb_[lev_max_glb_+1] << std::endl;
            //    std::cout << "Required mem. size (etc_glb_)               = " << sizeof(etcLM) * adr_tc_level_partition_glb_[lev_max_glb_+1] * 1.0e-9 << std::endl;
            //}
            //psys.dumpMemSizeUsed(std::cout);
            //dumpMemSizeUsed(std::cout);
            //athread_halt();
            //Finalize();
            //std::exit(0);
#endif
            // not needed, but we can measure the performance at n_loop=0
            etc_loc_.resizeNoInitialize(adr_tc_level_partition_loc_[lev_max_loc_+1]);
            etc_glb_.resizeNoInitialize(adr_tc_level_partition_glb_[lev_max_glb_+1]);
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK17 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif

#if 0
            /*
            athread_halt();
            Finalize();
            std::exit(0);
            */
#endif
        }
        else{
            assert(false);
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
#ifdef PHI_R_TREE
        const F64 len_peri_x = 2.0 * 4.0 * atan(1.0);
#endif
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

#ifdef DEBUG_PRINT_MAKE_LIST_ID3
        Comm::barrier();
        if(Comm::getRank()==0) {
            std::cerr<<"OK01 @makeAllInteractionListId3"<<std::endl;
            std::cerr<<"n_ipg= "<<n_ipg
                     <<" r_crit_sq= "<<r_crit_sq
                     <<" sp_tree_offset= "<<sp_tree_offset
                     <<" len_peri_x= "<<len_peri_x
                     <<std::endl;
        }
        Comm::barrier();
#endif //DEBUG_PRINT_MAKE_LIST_ID3

        //#ifdef SUNWAY
#if 1
        //***** CPE versions w/ stack *****
        //* Resize memories of *recorder_for_force_* (currently, we use *tmp* for test)
        enum {
            N_CPE = 64,
            N_JP_PER_GROUP = 2048,
        };
        id_epj_recorder_for_force_[0].reserve(N_JP_PER_GROUP * n_ipg);
        id_spj_recorder_for_force_[0].reserve(N_JP_PER_GROUP * n_ipg);
        n_disp_epj_recorder_for_force_[0].reserve(n_ipg);
        n_disp_spj_recorder_for_force_[0].reserve(n_ipg);
        n_epj_recorder_for_force_[0].reserve(n_ipg);
        n_spj_recorder_for_force_[0].reserve(n_ipg);
//#ifdef USE_SWAP_IN_MAKING_INTERACTION_LIST
//        //* Prepare swap regions for id lists in the main memory
//        enum {
//            SIZE_EPJ_SWAP = 32768,
//            SIZE_SPJ_SWAP = 8192,
//        };
//        int size_epj_swap[N_CPE];
//        int size_spj_swap[N_CPE];
//        ReallocatableArray<int> id_epj_swap[N_CPE];
//        ReallocatableArray<int> id_spj_swap[N_CPE];
//        unsigned long adr_id_epj_swap[N_CPE];
//        unsigned long adr_id_spj_swap[N_CPE];
//        int f_realloc_epj_swap[N_CPE];
//        int f_realloc_spj_swap[N_CPE];
//        for (S32 cpe_id=0; cpe_id<N_CPE; cpe_id++) {
//            id_epj_swap[cpe_id].reserve(SIZE_EPJ_SWAP);
//            id_spj_swap[cpe_id].reserve(SIZE_SPJ_SWAP);
//        }
//#endif

        const F64 wtime_start = MPI::Wtime();
        //* Make interaction lists
        S32 n_epj_shortfall, n_spj_shortfall;
        S32 ca_id_epj, ca_id_spj;
        for(;;) {
            //** Reset
            n_epj_shortfall = 0;
            n_spj_shortfall = 0;
            ca_id_epj = id_epj_recorder_for_force_[0].capacity();
            ca_id_spj = id_spj_recorder_for_force_[0].capacity();
//#ifdef USE_SWAP_IN_MAKING_INTERACTION_LIST
//            for (S32 cpe_id=0; cpe_id<N_CPE; cpe_id++) {
//               size_epj_swap[cpe_id] = id_epj_swap[cpe_id].capacity();
//               size_spj_swap[cpe_id] = id_spj_swap[cpe_id].capacity();
//               adr_id_epj_swap[cpe_id] = (unsigned long) id_epj_swap[cpe_id].getPointer();
//               adr_id_spj_swap[cpe_id] = (unsigned long) id_spj_swap[cpe_id].getPointer();
//               f_realloc_ep_swap[cpe_id] = 0;
//               f_realloc_sp_swap[cpe_id] = 0;
//            }
//#endif

            //** Set the function arguments
            unsigned long args[32];
            args[0]  = (unsigned long) Comm::getRank();
            args[1]  = (unsigned long) n_ipg; 
            args[2]  = (unsigned long) ipg_.getPointer();
            args[3]  = (unsigned long) tc_glb_.getPointer();
            args[4]  = (unsigned long) tp_glb_.getPointer();
            args[5]  = (unsigned long) id_epj_recorder_for_force_[0].getPointer();
            args[6]  = (unsigned long) id_spj_recorder_for_force_[0].getPointer();
            args[7]  = (unsigned long) n_disp_epj_recorder_for_force_[0].getPointer();
            args[8]  = (unsigned long) n_disp_spj_recorder_for_force_[0].getPointer();
            args[9]  = (unsigned long) n_epj_recorder_for_force_[0].getPointer();
            args[10] = (unsigned long) n_spj_recorder_for_force_[0].getPointer();
            args[11] = (unsigned long) &r_crit_sq;
            args[12] = (unsigned long) n_leaf_limit_;
            args[13] = (unsigned long) sp_tree_offset;
        #ifdef PHI_R_TREE
            args[14] = (unsigned long) &len_peri_x;
        #endif
            args[15] = (unsigned long) ca_id_epj;
            args[16] = (unsigned long) ca_id_spj;
            args[17] = (unsigned long) &n_epj_shortfall;
            args[18] = (unsigned long) &n_spj_shortfall;
//#ifdef USE_SWAP_IN_MAKING_INTERACTION_LIST
//            //-(information about the swap)
//            args[19] = (unsigned long) &size_epj_swap[0];
//            args[20] = (unsigned long) &size_spj_swap[0];
//            args[21] = (unsigned long) &adr_id_epj_swap[0];
//            args[22] = (unsigned long) &adr_id_spj_swap[0];
//            args[23] = (unsigned long) &f_realloc_epj_swap[0];
//            args[24] = (unsigned long) &f_realloc_spj_swap[0];
//#endif

            //** Make interaction lists on CPEs
            __real_athread_spawn((void*)slave_MakeListUsingTree, args);
            athread_join();

//#ifdef USE_SWAP_IN_MAKING_INTERACTION_LIST
//            //** Resize the swap regions if needed
//            S32 f_realloc = 0;
//            for (S32 cpe_id=0; cpe_id<N_CPE; cpe_id++) {
//                if (f_realloc_epj_swap[cpe_id] > 0) {
//                   S32 cap = id_epj_swap[cpe_id].capacity();
//                   id_epj_swap[cpe_id].reserve(cap + SIZE_EPJ_SWAP);
//                   f_realloc++; 
//                }
//                if (f_realloc_spj_swap[cpe_id] > 0) {
//                   S32 cap = id_spj_swap[cpe_id].capacity();
//                   id_spj_swap[cpe_id].reserve(cap + SIZE_SPJ_SWAP);
//                   f_realloc++;
//                }
//            }
//            if (f_realloc > 0) continue;
//#endif

            //** Check the termination condition
            if ((n_epj_shortfall == 0) && (n_spj_shortfall == 0)) break;

            //* Extend id_*_recorder_for_force_[0] if needed
            id_epj_recorder_for_force_[0].reserve(ca_id_epj + n_epj_shortfall);
            id_spj_recorder_for_force_[0].reserve(ca_id_spj + n_spj_shortfall);
        }
        const F64 wtime_end = MPI::Wtime();
        //if (Comm::getRank() == 0) std::cout << "wtime (new) = " << wtime_end - wtime_start << std::endl;
#else //SUNWAY
    #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4)
    #endif
        for(S32 ig=0; ig<n_ipg; ig++){
            const S32 ith = Comm::getThreadNum();
            const S32 adr_epj_prev = id_epj_recorder_for_force_[ith].size();
            const S32 adr_spj_prev = id_spj_recorder_for_force_[ith].size();
            const F64ort pos_target_box = (ipg_[ig]).vertex_;
    #ifdef DEBUG_PRINT_MAKE_LIST_ID3
            if(Comm::getRank()==0) {
                std::cerr<<"OK02 @makeAllInteractionListId3"<<std::endl;
                std::cerr<<"ig= "<<ig<<std::endl;
                std::cerr<<"pos_target_box= "<<pos_target_box<<std::endl;
            }
    #endif //DEBUG_PRINT_MAKE_LIST_ID3

            ReallocatableArray<Tepj> ep_tmp;
            ReallocatableArray<Tspj> sp_tmp;
    #ifdef PHI_R_TREE
            MakeListUsingTreeRecursive<TSM, TreeCell<Tmomglb>, TreeParticle, Tepj, Tspj>
                (tc_glb_, 0, tp_glb_, 
                 epj_sorted_, id_epj_recorder_for_force_[ith], 
                 spj_sorted_, id_spj_recorder_for_force_[ith], 
                 pos_target_box, r_crit_sq, n_leaf_limit_, sp_tree_offset, len_peri_x);
    #else
            MakeListUsingTreeRecursive<TSM, MAKE_LIST_MODE_INTERACTION, LIST_CONTENT_ID, TreeCell<Tmomglb>, TreeParticle, Tepj, Tspj>
                (tc_glb_, 0, tp_glb_, 
                 epj_sorted_, ep_tmp, id_epj_recorder_for_force_[ith], 
                 spj_sorted_, sp_tmp, id_spj_recorder_for_force_[ith], 
                 pos_target_box, r_crit_sq, n_leaf_limit_, sp_tree_offset);
    #endif

    #ifdef DEBUG_PRINT_MAKE_LIST_ID3
            if(Comm::getRank()==0) {
                std::cerr<<"OK03 @makeAllInteractionListId3"<<std::endl;
                std::cerr<<"id_epj_recorder_for_force_[ith].size()= "<<id_epj_recorder_for_force_[ith].size()
                         <<" id_spj_recorder_for_force_[ith].size()= "<<id_spj_recorder_for_force_[ith].size()
                         <<std::endl;
            }
    #endif //DEBUG_PRINT_MAKE_LIST_ID3            
            n_disp_epj_recorder_for_force_[ith].push_back( id_epj_recorder_for_force_[ith].size() );
            n_disp_spj_recorder_for_force_[ith].push_back( id_spj_recorder_for_force_[ith].size() );

            const S32 n_epj = id_epj_recorder_for_force_[ith].size() - adr_epj_prev;
            const S32 n_spj = id_spj_recorder_for_force_[ith].size() - adr_spj_prev;
            n_epj_recorder_for_force_[ith].push_back(n_epj);
            n_spj_recorder_for_force_[ith].push_back(n_spj);

    #ifdef DEBUG_PRINT_MAKE_LIST_ID3
            if(Comm::getRank()==0) {
                std::cerr<<"OK04 @makeAllInteractionListId3"<<std::endl;
            }
            F64 cm_mass = 0.0;
            F64vec cm_pos = 0.0;
            for(S32 i=adr_epj_prev; i<id_epj_recorder_for_force_[ith].size(); i++){
                S32 adr = id_epj_recorder_for_force_[ith][i];
                cm_mass += epj_sorted_[adr].mass;
                cm_pos  += epj_sorted_[adr].mass*epj_sorted_[adr].pos;
            }
            if(Comm::getRank()==0) {
                std::cerr<<"OK05 @makeAllInteractionListId3"<<std::endl;
                std::cerr<<"cm_mass= "<<cm_mass<<std::endl;
            }
            for(S32 i=adr_spj_prev; i<id_spj_recorder_for_force_[ith].size(); i++){
                S32 adr = id_spj_recorder_for_force_[ith][i];
                cm_mass += spj_sorted_[adr].mass;
                cm_pos  += spj_sorted_[adr].mass*spj_sorted_[adr].pos;
            }
            if(Comm::getRank()==0) {
                std::cerr<<"OK06 @makeAllInteractionListId3"<<std::endl;
                std::cerr<<"cm_mass= "<<cm_mass<<std::endl;
            }
            // FOR DEBUG
            cm_pos /= cm_mass;
            if(Comm::getRank() == 0){
                std::cerr<<"OK07 @makeAllInteractionListId3"<<std::endl;
                std::cerr<<"n_epj= "<<n_epj<<" n_spj= "<<n_spj<<std::endl;
                std::cerr<<"ipg_[ig].n_ptcl_= "<<ipg_[ig].n_ptcl_<<std::endl;
                std::cerr<<"ipg_[ig].vertex_= "<<ipg_[ig].vertex_<<std::endl;
                std::cerr<<"cm_mass= "<<cm_mass
                         <<" cm_pos= "<<cm_pos
                         <<std::endl;
            }
    #endif
        }
#endif
        //exit(1);
        
#ifdef DEBUG_PRINT_MAKE_LIST_ID3
        Comm::barrier(); if(Comm::getRank()==0) {std::cerr<<"OK08 @makeAllInteractionListId3"<<std::endl;}
#endif //DEBUG_PRINT_MAKE_LIST_ID3
        time_profile_.make_all_interaction_list_id += GetWtime() - time_offset;

#if 0
        //* Check
        if (Comm::getRank() == 0) {
            std::ofstream ofs;
            std::string fname;
            //** id_epj_tmp
            fname = "./id_epj_old.dat";
            //fname = "./id_epj_new.dat";
            ofs.open(fname.c_str(), std::ios::trunc);
            for (S32 ig=0; ig<n_ipg; ig++) {
                const S32 n_epj = n_epj_recorder_for_force_[0][ig];
                const S32 n_disp_epj = n_disp_epj_recorder_for_force_[0][ig];
                for (S32 j=0; j<n_epj; j++) {
                    const S32 jp = n_disp_epj + j;
                    ofs << id_epj_recorder_for_force_[0][jp] << std::endl;
                }
            }
            ofs.close();
            //** id_spj_tmp
            fname = "./id_spj_old.dat";
            //fname = "./id_spj_new.dat";
            ofs.open(fname.c_str(), std::ios::trunc);
            for (S32 ig=0; ig<n_ipg; ig++) {
                const S32 n_spj = n_spj_recorder_for_force_[0][ig];
                const S32 n_disp_spj = n_disp_spj_recorder_for_force_[0][ig];
                for (S32 j=0; j<n_spj; j++) {
                    const S32 jp = n_disp_spj + j;
                    ofs << id_spj_recorder_for_force_[0][jp] << std::endl;
                }
            }
            ofs.close();
        }
#endif
        /*
        // For a test run
        MPI::COMM_WORLD.Barrier();
        if (Comm::getRank() == 0) std::cout << "makeAllInteractionListId3() ended." << std::endl;
        athread_halt();
        Finalize();
        std::exit(0);
        */
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    balanceNwalk(){
        Comm::barrier();
        const double start_time = MPI_Wtime();

        // Compute the cost of each walk
        std::vector<CostInfo> cost_info;
        cost_info.reserve(ipg_.size());
        for (int i=0; i<ipg_.size(); i++) {
            CostInfo tmp;
            tmp.gid_   = i;
            tmp.ncost_ = cpe_cost(n_epi_recorder_for_force_[0][i],
                                  n_epj_recorder_for_force_[0][i],
                                  n_spj_recorder_for_force_[0][i],
                                  N_SATELLITE,
                                  54,30,100,32,237404);
            cost_info.push_back(tmp);
        }
#if 0
        // [Check]
        if (Comm::getRank() == 0) {
            std::cout << "outputing cost_info before sorting...." << std::endl;
            std::ofstream ofs;
            std::stringstream filenum;
            std::string filename;
            filenum << std::setfill('0') << std::setw(5) << Comm::getRank();
            filename = "./cost_info_org" + filenum.str() + ".txt";
            std::cout << "ipg_.size()      = " << ipg_.size() << std::endl;
            std::cout << "cost_info.size() = " << cost_info.size() << std::endl;
            ofs.open(filename.c_str(), std::ios::trunc);
            for (S32 i=0; i<cost_info.size(); i++) {
                ofs << cost_info[i].gid_ << " " 
                    << cost_info[i].ncost_ << std::endl;
            }
            ofs.close();
            std::cout << "output of cost_info_org is completed!" << std::endl;
        }
#endif
   
        // Sort
        std::sort(cost_info.begin(), cost_info.end(), CostInfo());
#if 0
        // [Check]
        if (Comm::getRank() == 0) {
            std::cout << "outputing cost_info after sorting...." << std::endl;
            std::ofstream ofs;
            std::stringstream filenum;
            std::string filename;
            filenum << std::setfill('0') << std::setw(5) << Comm::getRank();
            filename = "./cost_info_sorted" + filenum.str() + ".txt";
            std::cout << "ipg_.size()      = " << ipg_.size() << std::endl;
            std::cout << "cost_info.size() = " << cost_info.size() << std::endl;
            ofs.open(filename.c_str(), std::ios::trunc);
            for (S32 i=0; i<cost_info.size(); i++) {
                ofs << cost_info[i].gid_ << " " 
                    << cost_info[i].ncost_ << std::endl;
            }
            ofs.close();
            std::cout << "output of cost_info_sorted is completed!" << std::endl;
        }
#endif

        // Assign interaction lists to each CPE
        if (ipg_.size() > NUMBER_OF_CPE) {
            // Initialize cost_info.cpe_id_ 
            for (S32 i=0; i<NUMBER_OF_CPE; i++) 
                cost_info[i].cpe_id_ = i;
            // Make a binary tree in which each tree cell has the minimum value
            // of cost among the child tree cells. Here, we assume that the 
            // the number of CPEs is the power of 2.
            //
            // [Structure of binary tree]
            //  - The total number of tree cells are 126.
            //  - The tree cells are structured as follows:
            //    2 + 4 + 8 + 16 + 32 + 64 
            //  in the width-first order.
            //  - To access the child cells, 
            //  - The leaf cells has the sum of cost for each CPE.
            U64 csum_tree[2*NUMBER_OF_CPE]; // two elements are unused.
            const U32 lev_max = 6; // 2^{6} = 64;
            for (S32 i=0; i<NUMBER_OF_CPE; i++) {
                S32 offset = 2 * (1 << (lev_max-1) - 1); // 2*(2^{lev_max-1}-1)
                csum_tree[offset+i] = cost_info[i].ncost_;
            }
#if 0
            for (S32 lev=lev_max-1; lev>0; lev--) {
                S32 offset = 2 * (1 << (lev-1) - 1); 
                S32 nterms = 1 << lev;
                for (S32 i=0; i<nterms; i += 2) {
                   S32 left = offset + i;
                   S32 right = offset + i + 1;
                }
            }
#endif
            // Compute cost_info.cpe_id & ncost_sum[]
            for (S32 i=NUMBER_OF_CPE; i<ipg_.size(); i++) {
                // Locate the index whose ncost_sum is minimum
#if 0
                cpe_id = find_
#else
                // In this case, we use sequential search to find the minimum.
                S32 offset = 2 * (1 << (lev_max-1) - 1);
                S32 cpe_id = 0;
                U64 csum_ref = csum_tree[offset];
                for (S32 k=1; k<(1<<lev_max); k++) {
                    U64 csum = csum_tree[offset+k];
                    if (csum < csum_ref) {
                        cpe_id = k;
                        csum_ref = csum;
                    }
                }
#endif
                // Set cost_info.cpe_id_
                cost_info[i].cpe_id_ = cpe_id;
                // Update the csum tree
                offset = 2 * (1 << (lev_max-1) - 1);
                csum_tree[offset+cpe_id] += cost_info[i].ncost_;
            }
            // Count # of walks processed in each CPE
            for (S32 cpe_id=0; cpe_id<NUMBER_OF_CPE; cpe_id++)
                n_walk_cpe[cpe_id] = 0;
            for (S32 i=0; i<ipg_.size(); i++) {
                S32 cpe_id = cost_info[i].cpe_id_;
                n_walk_cpe[cpe_id]++;
            }
            // Compute n_disp_walk[]
            n_disp_walk[0] = 0;
            for (S32 cpe_id=1; cpe_id<NUMBER_OF_CPE; cpe_id++)
                n_disp_walk[cpe_id] = n_disp_walk[cpe_id-1] + n_walk_cpe[cpe_id-1];
            // Reset n_walk_cpe[]
            for (S32 cpe_id=0; cpe_id<NUMBER_OF_CPE; cpe_id++)
                n_walk_cpe[cpe_id] = 0;
            // Compute adr_n_walk[]
            for (S32 i=0; i<ipg_.size(); i++) {
                S32 cpe_id = cost_info[i].cpe_id_;
                S32 adr = n_disp_walk[cpe_id] + n_walk_cpe[cpe_id];
                adr_n_walk[adr] = cost_info[i].gid_;
                n_walk_cpe[cpe_id]++;
            }
        } else {
            // In this case, we use ipg_.size() CPEs only
            for (S32 i=0; i<ipg_.size(); i++) n_walk_cpe[i]=1;
            n_disp_walk[0] = 0;
            for (S32 i=1; i<NUMBER_OF_CPE; i++)
                n_disp_walk[i] = n_disp_walk[i-1] + n_walk_cpe[i-1];
            for (S32 i=0; i<ipg_.size(); i++) adr_n_walk[i] = i;
        }
#if 0
        // [Check]
        if (Comm::getRank() == 0) {
            std::cout << "outputing adr_n_walk ...." << std::endl;
            std::ofstream ofs;
            std::stringstream filenum;
            std::string filename;
            filenum << std::setfill('0') << std::setw(5) << Comm::getRank();
            filename = "./adr_n_walk" + filenum.str() + ".txt";
            std::cout << "ipg_.size()      = " << ipg_.size() << std::endl;
            ofs.open(filename.c_str(), std::ios::trunc);
            for (S32 cpe_id=0; cpe_id<NUMBER_OF_CPE; cpe_id++) {
               S32 i_start = n_disp_walk[cpe_id];
               S32 i_end   = i_start + n_walk_cpe[cpe_id];
               U64 csum = 0;
               for (S32 i=i_start; i<i_end; i++) {
                   S32 ig = adr_n_walk[i];
                   U64 ncost = cpe_cost(n_epi_recorder_for_force_[0][ig],
                                        n_epj_recorder_for_force_[0][ig],
                                        n_spj_recorder_for_force_[0][ig],
                                        N_SATELLITE,
                                        54,30,100,32,237404);
                   ofs << ig << " " << ncost << std::endl;
                   csum += ncost;
               }
               ofs << "(ncost_sum = " << csum << ")" << std::endl;
            }
            ofs.close();
            std::cout << "output of adr_n_walk is completed!" << std::endl;
        }
#endif

#if 0
        // [Debug]
        Comm::barrier();
        const double end_time = MPI_Wtime();
        if (Comm::getRank() == 0)
            std::cout << "wtime(balanceNwalk) = " 
                      << end_time - start_time << std::endl;
        athread_halt();
        Finalize();
        std::exit(0);
#endif
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
        S64 n_epj_tot = 0;
        S64 n_spj_tot = 0;
        n_epi_recorder_for_force_[0].resizeNoInitialize(n_ipg);
        n_disp_epi_recorder_for_force_[0].resizeNoInitialize(n_ipg+1);
        n_disp_epi_recorder_for_force_[0][0] = 0;
        for(S32 i=0; i<n_ipg; i++){
            n_epi_recorder_for_force_[0][i] = ipg_[i].n_ptcl_;
            n_disp_epi_recorder_for_force_[0][i+1] = n_disp_epi_recorder_for_force_[0][i] + n_epi_recorder_for_force_[0][i];
            n_interaction_ep_ep_local_ += (S64)n_epi_recorder_for_force_[0][i] * (S64)n_epj_recorder_for_force_[0][i];
            n_interaction_ep_sp_local_ += (S64)n_epi_recorder_for_force_[0][i] * (S64)n_spj_recorder_for_force_[0][i];
            n_epj_tot += (S64)n_epj_recorder_for_force_[0][i];
            n_spj_tot += (S64)n_spj_recorder_for_force_[0][i];
        }
        const S32 ni_tot = n_disp_epi_recorder_for_force_[0][n_ipg];
        n_epi_loc_ave_ = (F64)ni_tot / n_ipg;
        n_epj_loc_ave_ = (F64)n_epj_tot / n_ipg;
        n_spj_loc_ave_ = (F64)n_spj_tot / n_ipg;        
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
        /*
        const S32 n_epj_tot = n_disp_epj_recorder_for_force_[0][n_ipg];
        const S32 n_spj_tot = n_disp_spj_recorder_for_force_[0][n_ipg];
        n_epi_loc_ave_ = ni_tot / n_ipg;
        n_epj_loc_ave_ = n_epj_tot / n_ipg;
        n_spj_loc_ave_ = n_spj_tot / n_ipg;
        */
        time_profile_.calc_force = GetWtime() - wtime_offset_out;
        return ret;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tpsys>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setParticleLocalTreeImpl(const Tpsys & psys,
                             ReallocatableArray<Tepj> & epj,
                             const bool clear){
        const F64 time_offset = GetWtime();
        const S32 nloc = psys.getNumberOfParticleLocal();
        if(clear){ n_loc_tot_ = 0;}
        const S32 offset = n_loc_tot_;
        n_loc_tot_ += nloc;
        epj.resizeNoInitialize(n_loc_tot_);
#if 0
        //[DEBUG:start]
        Comm::barrier();
        if (Comm::getRank() == 0) std::cout << "Start point of setParticleLocalTreeImpl()" << std::endl;
        S32 num_err = 0;
        for (S32 i=0; i<nloc; i++) {
            F64vec pos = psys[i].getPos();
            if (fabs(pos.z) < 1.0e-300) num_err++;
        }
        if (num_err /= 0)
            std::cout << "num_err = " << num_err << " "
                      << "(RANK: " << Comm::getRank() << ")" << std::endl;
        Finalize();
        std::exit(0);
        //[DEBUG:end]
#endif
        if(clear){
#ifdef SUNWAY
            unsigned long args[3];
            args[0] = (unsigned long)nloc;
            args[1] = (unsigned long)(&psys[0]);
            args[2] = (unsigned long)epj.getPointer();
            __real_athread_spawn((void *)slave_CopyFPToEPJ, args);
            athread_join();
#else //SUNWAY
            for(S32 i=0; i<nloc; i++){
                epj[i].pos  = psys[i].pos;
                epj[i].mass = psys[i].mass;
                epj[i].vel  = psys[i].vel;
                epj[i].id   = psys[i].id;
            }
#endif //SUNWAY            
        }
        else{
            for(S32 i=0; i<nloc; i++){
                epj[i+offset].copyFromFP( psys[i] );
            }
        }
#if 0
        //[DEBUG:start]
        Comm::barrier();
        if (Comm::getRank() == 0) std::cout << "Checking epj in setParticleLocalTreeImpl()" << std::endl;
        S32 num_err = 0;
        for (S32 i=0; i<nloc; i++) {
            F64vec pos = epj[i].getPos();
            if (fabs(pos.z) < 1.0e-300) num_err++;
        }
        if (num_err /= 0)
            std::cout << "num_err = " << num_err << " "
                      << "(RANK: " << Comm::getRank() << ")" << std::endl;
        Finalize();
        std::exit(0);
        //[DEBUG:end]
#endif

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

#ifdef USE_STD_SORT
        std::sort(tp_glb_.getPointer(), tp_glb_.getPointer()+n_glb_tot_, 
                  LessOPKEY());
#elif defined(USE_CPE_SORT) 
        tp_buf_.resizeNoInitialize(n_glb_tot_);
        //ReallocatableArray< TreeParticle > tp_work_;
        //tp_work_.resizeNoInitialize(n_glb_tot_);
        //        std::cerr<<"rank "<<Comm::getRank()<<" glb sort started "<<std::endl;

//        std::memcpy(tp_buf_.getPointer(),tp_glb_.getPointer(),sizeof(TreeParticle)*n_glb_tot_);
//        std::sort(tp_buf_.getPointer(), tp_buf_.getPointer()+n_glb_tot_, 
//        //        std::sort(tp_glb_.getPointer(), tp_glb_.getPointer()+n_glb_tot_, 
//                  LessOPKEY());

        samplesort(n_glb_tot_,tp_glb_.getPointer(),tp_glb_.getPointer(),tp_buf_.getPointer());
        //        std::memcpy(tp_glb_.getPointer(),tp_buf_.getPointer(),sizeof(TreeParticle)*n_glb_tot_);

//        std::cerr<<"before\n";
//        tp_glb_.print("glb");
//        tp_g_buf_.print("buf");
//        std::swap(tp_glb_,tp_g_buf_);
//        std::cerr<<"after\n";
//        tp_glb_.print("glb");
//        tp_g_buf_.print("buf");
        
        //        std::cerr<<"rank "<<Comm::getRank()<<" swap n="<<n_glb_tot_<<std::endl;
        //        fflush(stderr);
//    	std::memcpy(&tp_work_,&tp_glb_, sizeof(tp_work_));
//    	std::memcpy(&tp_glb_, &tp_buf_, sizeof(tp_work_));
//    	std::memcpy(&tp_buf_, &tp_work_,sizeof(tp_work_));
//        std::cerr<<"tp_work_.size"<<tp_work_.size()<<std::endl;
//        std::cerr<<"tp_glb_.size "<<tp_glb_.size()<<std::endl;
//        std::cerr<<"tp_buf_.size "<<tp_buf_.size()<<std::endl;
        //        std::memcpy(tp_glb_.getPointer(),tp_buf_.getPointer(),sizeof(TreeParticle)*n_glb_tot_);
        // std::swap(tp_glb_,tp_buf_);
        //        std::cerr<<"rank "<<Comm::getRank()<<" glb sort finished "<<std::endl;
        //        fflush(stderr);
#else //USE_STD_SORT
        tp_buf_.resizeNoInitialize(n_glb_tot_);
        rs_.lsdSort(tp_glb_.getPointer(), tp_buf_.getPointer(), 0, n_glb_tot_-1);
#endif //USE_STD_SORT
#ifdef REDUCE_MEMORY
        epj_sorted_.resizeNoInitialize( epj_org_.size() );
#else
        epj_sorted_.resizeNoInitialize( n_glb_tot_ );
#endif
        //std::cerr<<"check 1"<<std::endl;
        if( typeid(TSM) == typeid(SEARCH_MODE_LONG)
            || typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF) 
            || typeid(TSM) == typeid(SEARCH_MODE_LONG_SCATTER) 
            || typeid(TSM) == typeid(SEARCH_MODE_LONG_CUTOFF_SCATTER) ){
#ifdef DEBUG_PRINT_SORT_GLB_3
            Comm::barrier();
            if(Comm::getRank()==0){
                std::cerr<<"tp_glb_.size()= "<<tp_glb_.size()<<std::endl;
                for(S32 i=0; i<tp_glb_.size(); i++){
                    std::cerr<<"i= "<<i
                             <<" tp_glb_[i].key= "<<tp_glb_[i].key_
                             <<" adr_ptcl= "<<tp_glb_[i].adr_ptcl_
                             <<std::endl;
                    if(GetMSB(tp_glb_[i].adr_ptcl_)==0){
                        std::cerr<<"adr_ptcl_2= "<<tp_glb_[i].adr_ptcl_
                                 <<" epj_org.mass= "<<epj_org_[tp_glb_[i].adr_ptcl_].mass
                                 <<" pos= "<<epj_org_[tp_glb_[i].adr_ptcl_].pos
                                 <<std::endl;
                    }
                    else{
                        std::cerr<<"adr_ptcl_2= "<<ClearMSB(tp_glb_[i].adr_ptcl_)
                                 <<" spj_org.mass= "<<spj_org_[ClearMSB(tp_glb_[i].adr_ptcl_)].mass
                                 <<" pos= "<<spj_org_[ClearMSB(tp_glb_[i].adr_ptcl_)].pos
                                 <<std::endl;
                    }
                }
            }
            Comm::barrier();
            //exit(1);
#endif //DEBUG_PRINT_SORT_GLB_3
            // Long mode
            spj_sorted_.resizeNoInitialize( spj_org_.size() );
            adr_epj_org2glb_.resizeNoInitialize(n_loc_tot_); // new line
            const S32 n_proc = Comm::getNumberOfProc();
#ifdef USE_SUPER_DOMAIN
            const S32 n_epj_add = comm_table_.ep_recv_.size();
            const S32 n_spj_add = comm_table_.sp_recv_.size();
#else
            const S32 n_epj_add = comm_table_.n_disp_ep_recv_[n_proc];
            const S32 n_spj_add = comm_table_.n_disp_sp_recv_[n_proc];
#endif
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
    #ifdef REDUCE_MEMORY
            S32 n_cnt_ep = 0;
            S32 n_cnt_sp = 0;
            for(S32 i=0; i<n_glb_tot_; i++){
                const U32 adr = tp_glb_[i].adr_ptcl_;
                if( GetMSB(adr) == 0){
                    epj_sorted_[n_cnt_ep] = epj_org_[adr];
                    tp_glb_[i].adr_ptcl_ = n_cnt_ep;
                    // new lines
                    if( adr < n_loc_tot_){
                        // particles in own process
                        adr_epj_org2glb_[adr] = n_cnt_ep;
                    }
                    else{
                        // particles from other process
                        adr_epj_buf2glb_[adr-n_loc_tot_] = n_cnt_ep; 
                    }
                    n_cnt_ep++;
                    // new lines
                }
                else{
                    const U32 adr_new = ClearMSB(adr);
                    spj_sorted_[n_cnt_sp] = spj_org_[adr_new];
                    adr_spj_buf2glb_[adr_new] = n_cnt_sp; // new line
                    tp_glb_[i].adr_ptcl_ = SetMSB((U32)n_cnt_sp);
                    n_cnt_sp++;
                }
            }
            // new lines
            adr_epj_loc2glb_.resizeNoInitialize(n_loc_tot_);
            for(S32 i=0; i<n_loc_tot_; i++){
                const U32 adr_org = tp_loc_[i].adr_ptcl_;
                const U32 adr_glb = adr_epj_org2glb_[adr_org];
                adr_epj_loc2glb_[i] = adr_glb;
            }
            // The following is added by DN
            adr_epj_glb2loc_.resizeNoInitialize(epj_sorted_.size());
            for (S32 i=0; i<epj_sorted_.size(); i++) adr_epj_glb2loc_[i] = -1;
            for (S32 i=0; i<adr_epj_loc2glb_.size(); i++) {
                const S32 adr = adr_epj_loc2glb_[i];
                adr_epj_glb2loc_[adr] = i;
            }
    #else // REDUCE_MEMORY
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
            // The following is added by DN
            adr_epj_glb2loc_.resizeNoInitialize(epj_sorted_.size());
            for (S32 i=0; i<epj_sorted_.size(); i++) adr_epj_glb2loc_[i] = -1;
            for (S32 i=0; i<adr_epj_loc2glb_.size(); i++) {
                const S32 adr = adr_epj_loc2glb_[i];
                adr_epj_glb2loc_[adr] = i;
            }
    #endif // REDUCE_MEMORY

    #ifdef DEBUG_PRINT_SORT_GLB_3
            Comm::barrier();
            if(Comm::getRank()==0){
                std::cerr<<"tp_glb_.size()= "<<tp_glb_.size()<<std::endl;
                for(S32 i=0; i<tp_glb_.size(); i++){
                    std::cerr<<"i= "<<i
                             <<" tp_glb_[i].key= "<<tp_glb_[i].key_
                             <<" adr_ptcl= "<<tp_glb_[i].adr_ptcl_
                             <<std::endl;
                    if(GetMSB(tp_glb_[i].adr_ptcl_)==0){
                        std::cerr<<"adr_ptcl_2= "<<tp_glb_[i].adr_ptcl_
                                 <<" epj_sorted.mass= "<<epj_sorted_[tp_glb_[i].adr_ptcl_].mass
                                 <<" pos= "<<epj_sorted_[tp_glb_[i].adr_ptcl_].pos
                                 <<std::endl;
                    }
                    else{
                        std::cerr<<"adr_ptcl_2= "<<ClearMSB(tp_glb_[i].adr_ptcl_)
                                 <<" spj_sorted.mass= "<<spj_sorted_[ClearMSB(tp_glb_[i].adr_ptcl_)].mass
                                 <<" pos= "<<spj_sorted_[ClearMSB(tp_glb_[i].adr_ptcl_)].pos
                                 <<std::endl;
                    }
                }
            }
            Comm::barrier();
            //exit(1);
    #endif //DEBUG_PRINT_SORT_GLB_3
            
#if VERSION_MSORT_GLB_TREE_REUSE == 0
#elif VERSION_MSORT_GLB_TREE_REUSE == 1
            unsigned long args[8];
            //* Count the number of groups 
            int n_groups[NUMBER_OF_CPE];
            args[0] = (unsigned long) Comm::getRank();
            args[1] = (unsigned long) n_loc_tot_;
            args[2] = (unsigned long) adr_epj_loc2glb_.getPointer();
            args[3] = (unsigned long) n_groups;
            __real_athread_spawn((void *)slave_Preproc1_CopyEPJLocToEPJGlb, args);
            athread_join();
            //* Reallocate memory if needed
            unsigned long adr_adr_epj_loc_group_head_[NUMBER_OF_CPE];
            unsigned long adr_adr_epj_glb_group_head_[NUMBER_OF_CPE];
            unsigned long adr_group_size_epj_loc_[NUMBER_OF_CPE];
            for (S32 cpe_id=0; cpe_id<NUMBER_OF_CPE; cpe_id++) {
                adr_epj_loc_group_head_[cpe_id].resizeNoInitialize(n_groups[cpe_id]);
                adr_epj_glb_group_head_[cpe_id].resizeNoInitialize(n_groups[cpe_id]);
                group_size_epj_loc_[cpe_id].resizeNoInitialize(n_groups[cpe_id]);
                adr_adr_epj_loc_group_head_[cpe_id] = (unsigned long) adr_epj_loc_group_head_[cpe_id].getPointer();
                adr_adr_epj_glb_group_head_[cpe_id] = (unsigned long) adr_epj_glb_group_head_[cpe_id].getPointer();
                adr_group_size_epj_loc_[cpe_id]     = (unsigned long) group_size_epj_loc_[cpe_id].getPointer();
            }
            //* Set group information
            args[0] = (unsigned long) Comm::getRank();
            args[1] = (unsigned long) n_loc_tot_;
            args[2] = (unsigned long) adr_epj_loc2glb_.getPointer();
            args[3] = (unsigned long) adr_adr_epj_loc_group_head_;
            args[4] = (unsigned long) adr_adr_epj_glb_group_head_;
            args[5] = (unsigned long) adr_group_size_epj_loc_;
            __real_athread_spawn((void *)slave_Preproc2_CopyEPJLocToEPJGlb, args);
            athread_join();
#elif VERSION_MSORT_GLB_TREE_REUSE == 2
            unsigned long args[8];
            //* Count the number of groups 
            int n_groups[NUMBER_OF_CPE];
            args[0] = (unsigned long) Comm::getRank();
            args[1] = (unsigned long) n_glb_tot_;
            args[2] = (unsigned long) adr_epj_glb2loc_.getPointer();
            args[3] = (unsigned long) n_groups;
            __real_athread_spawn((void *)slave_Preproc1_CopyEPJLocToEPJGlb, args);
            athread_join();
            //* Reallocate memory if needed
            unsigned long adr_adr_epj_glb_group_head_[NUMBER_OF_CPE];
            unsigned long adr_adr_epj_loc_group_head_[NUMBER_OF_CPE];
            unsigned long adr_group_size_epj_glb_[NUMBER_OF_CPE];
            for (S32 cpe_id=0; cpe_id<NUMBER_OF_CPE; cpe_id++) {
                adr_epj_glb_group_head_[cpe_id].resizeNoInitialize(n_groups[cpe_id]);
                adr_epj_loc_group_head_[cpe_id].resizeNoInitialize(n_groups[cpe_id]);
                group_size_epj_glb_[cpe_id].resizeNoInitialize(n_groups[cpe_id]);
                adr_adr_epj_glb_group_head_[cpe_id] = (unsigned long) adr_epj_glb_group_head_[cpe_id].getPointer();
                adr_adr_epj_loc_group_head_[cpe_id] = (unsigned long) adr_epj_loc_group_head_[cpe_id].getPointer();
                adr_group_size_epj_glb_[cpe_id]     = (unsigned long) group_size_epj_glb_[cpe_id].getPointer();
            }
            //* Set group information
            args[0] = (unsigned long) Comm::getRank();
            args[1] = (unsigned long) n_glb_tot_;
            args[2] = (unsigned long) adr_epj_glb2loc_.getPointer();
            args[3] = (unsigned long) adr_adr_epj_glb_group_head_;
            args[4] = (unsigned long) adr_adr_epj_loc_group_head_;
            args[5] = (unsigned long) adr_group_size_epj_glb_;
            __real_athread_spawn((void *)slave_Preproc2_CopyEPJLocToEPJGlb, args);
            athread_join();
#if 0
            // [DEBUG]
            S32 n_group_max_loc = 0;
            for (S32 cpe_id=0; cpe_id<NUMBER_OF_CPE; cpe_id++)
                if (n_groups[cpe_id] > n_group_max_loc)
                    n_group_max_loc = n_groups[cpe_id];
            S32 n_group_max = 0;
            MPI_Allreduce(&n_group_max_loc, &n_group_max, 1, MPI_INT,
                          MPI_MAX, MPI_COMM_WORLD);
            if (Comm::getRank() == 0) std::cout << "n_group_max = " << n_group_max << std::endl;
            for (S32 cpe_id=0; cpe_id<NUMBER_OF_CPE; cpe_id++)
                if (n_groups[cpe_id] == n_group_max)
                    std::cout << "Rank " << Comm::getRank() << ", "
                              << "CPE " << cpe_id << " has n_group_max !" << std::endl;
            if (Comm::getRank() == 0) {
                std::ofstream ofs;
                std::stringstream filenum;
                std::string filename;
                filenum << std::setfill('0') << std::setw(5) << Comm::getRank();
                //filename = "group_info_old" + filenum.str()+ ".txt";
                filename = "group_info_new" + filenum.str()+ ".txt";
                ofs.open(filename.c_str(), std::ios::trunc);
                for (S32 cpe_id=0; cpe_id<NUMBER_OF_CPE; cpe_id++) {
                    ofs << "###################################################" << std::endl;
                    ofs << "CPE ID = " << cpe_id << ", "
                        << "n_group = " << n_groups[cpe_id] 
                        << std::endl;
                    ofs << "###################################################" << std::endl;
                    for (S32 ig=0; ig<n_groups[cpe_id]; ig++) {
                        ofs << adr_epj_glb_group_head_[cpe_id][ig] << " "
                            << adr_epj_loc_group_head_[cpe_id][ig] << " "
                            << group_size_epj_glb_[cpe_id][ig] << " "
                            << std::endl;
                    }
                }
                ofs.close();
            }
            //* For a test run
            Comm::barrier();
            if (Comm::getRank() == 0) std::cout << "Preprocess of CopyEPJLocToEPJGlb() completed!" << std::endl;
            athread_halt();
            Finalize();
            std::exit(0);
#endif
#else // VERSION_MSORT_GLB_TREE_REUSE
#error The value of the macro `VERSION_MSORT_GLB_TREE_REUSE` is incorrect.
#endif // VERSION_MSORT_GLB_TREE_REUSE

            
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
    mortonSortLocalTreeOnlyImpl(const ReallocatableArray<Tepj> & _epj_org,
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
            if (myrank == irank) 
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
            if (myrank == irank) 
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

#if 0
    Comm::barrier();
    const S32 target_rank=2;
    if (Comm::getRank() == target_rank) {
        std::cout << "outputing tp_loc_[] in mortonSortLocalTreeOnlyImpl()." << std::endl;
        std::ofstream ofs;
        std::stringstream filenum;
        std::string filename;
        filenum << std::setfill('0') << std::setw(5) << Comm::getRank();
        filename = "./tp_loc_org" + filenum.str() + ".dat";
        ofs.open(filename.c_str(), std::ios::binary|std::ios::trunc);
        S32 tp_num = tp_loc_.size();
        std::cout << "n_loc_tot_         = " << n_loc_tot_ << std::endl;
        std::cout << "tp_loc_.size()     = " << tp_num << std::endl;
        std::cout << "sizeof(tp_loc_[0]) = " << sizeof(tp_loc_[0]) << std::endl;
        ofs.write((char *)&tp_num, sizeof(S32));
        ofs.write((char *)&tp_loc_[0], tp_num * sizeof(tp_loc_[0]));
        ofs.close();
        std::cout << "output of tp_loc_[] is completed!" << std::endl;
    }
    athread_halt();
    Finalize();
    std::exit(0);
#endif
        
        wtime_offset_in = GetWtime();
    #ifdef USE_STD_SORT
	std::sort(tp_loc_.getPointer(), tp_loc_.getPointer()+n_loc_tot_, 
                  LessOPKEY());
#elif defined(USE_CPE_SORT) 
        tp_buf_.resizeNoInitialize(n_loc_tot_);
        //        std::cerr<<"rank "<<Comm::getRank()<<" loc sort started "<<std::endl;

//        std::memcpy(tp_buf_.getPointer(),tp_loc_.getPointer(),sizeof(TreeParticle)*n_loc_tot_);
//        std::sort(tp_buf_.getPointer(), tp_buf_.getPointer()+n_loc_tot_, 
//        std::sort(tp_loc_.getPointer(), tp_loc_.getPointer()+n_loc_tot_, 
//                  MoreOPKEY());
//        {
//    	std::ofstream flog;
//    	std::stringstream fnum;
//    	std::string fname;
//    	fnum << std::setfill('0') << std::setw(5) << Comm::getRank();
//    	fname = "catch" + fnum.str() + ".txt";
//    	flog.open(fname.c_str(),std::ios::trunc);
//        flog<<n_loc_tot_<<" ";
//    	for(int i=0; i<n_loc_tot_; i++) flog<<tp_loc_[i].key_<<" "<<tp_loc_[i].adr_ptcl_<<"\n";
//        flog.close();
//        }
//    	abort();

//        for (int k=0; k<4; k++ ) {
//            if(k==Comm::getRank()) {
//              printf("rank %d\n",k);
        samplesort(n_loc_tot_,tp_loc_.getPointer(),tp_loc_.getPointer(),tp_buf_.getPointer());
//            }
//            Comm::barrier();
//          }

//        std::memcpy(tp_loc_.getPointer(),tp_l_buf_.getPointer(),sizeof(TreeParticle)*n_loc_tot_);
        //        std::swap(tp_loc_,tp_buf_);
//        ReallocatableArray< TreeParticle > tp_temp;
//    	std::memcpy(&tp_temp, &tp_loc_, sizeof(tp_temp));
//    	std::memcpy(&tp_loc_, &tp_buf_, sizeof(tp_temp));
//    	std::memcpy(&tp_buf_, &tp_temp, sizeof(tp_temp));
        
        //std::cerr<<"rank "<<Comm::getRank()<<" loc sort finished "<<std::endl;
    #else //USE_STD_SORT
        tp_buf_.resizeNoInitialize(n_loc_tot_);
        rs_.lsdSort(tp_loc_.getPointer(), tp_buf_.getPointer(), 0, n_loc_tot_-1);
    #endif //USE_STD_SORT
        wtime_sort_local_tree = GetWtime() - wtime_offset_in;
        
    #ifdef SUNWAY
        // TODO
        unsigned long arg[5];
        arg[0] = (unsigned long)(n_loc_tot_);
        arg[1] = (unsigned long)(_epj_org.getPointer());
        arg[2] = (unsigned long)(_epi_sorted.getPointer());
        arg[3] = (unsigned long)(_epj_sorted.getPointer());
        arg[4] = (unsigned long)(tp_loc_.getPointer());
        __real_athread_spawn((void*)slave_CopyEpjOrgToEpiSortedEpjSorted, arg);
        athread_join();
        
        #ifdef DEBUG_SORT_LOCAL
        adr_ptcl_of_tp_loc_.resizeNoInitialize(n_loc_tot_);
        for(S32 i=0; i<n_loc_tot_; i++){
            const S32 adr = tp_loc_[i].adr_ptcl_;
            adr_ptcl_of_tp_loc_[i] = adr;
            assert(_epj_sorted[i].id    == _epj_org[adr].id);
            assert(_epj_sorted[i].mass  == _epj_org[adr].mass);
            assert(_epj_sorted[i].pos.x == _epj_org[adr].pos.x);
            assert(_epj_sorted[i].pos.y == _epj_org[adr].pos.y);
            assert(_epj_sorted[i].pos.z == _epj_org[adr].pos.z);
            assert(_epj_sorted[i].vel.x == _epj_org[adr].vel.x);
            assert(_epj_sorted[i].vel.y == _epj_org[adr].vel.y);
            assert(_epj_sorted[i].vel.z == _epj_org[adr].vel.z);
            assert(_epj_sorted[i].pos_phi == _epj_org[adr].pos_phi);
            assert(_epj_sorted[i].pos_r == _epj_org[adr].pos_r);
            
            assert(_epi_sorted[i].id    == _epj_org[adr].id);
            assert(_epi_sorted[i].mass  == _epj_org[adr].mass);
            assert(_epi_sorted[i].pos.x == _epj_org[adr].pos.x);
            assert(_epi_sorted[i].pos.y == _epj_org[adr].pos.y);
            assert(_epi_sorted[i].pos.z == _epj_org[adr].pos.z);
            assert(_epi_sorted[i].vel.x == _epj_org[adr].vel.x);
            assert(_epi_sorted[i].vel.y == _epj_org[adr].vel.y);
            assert(_epi_sorted[i].vel.z == _epj_org[adr].vel.z);

        }
        //exit(1);
        #endif

    #else //SUNWAY
        #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	
  #pragma omp parallel for
        #endif	
        for(S32 i=0; i<n_loc_tot_; i++){
            const S32 adr = tp_loc_[i].adr_ptcl_;
            _epj_sorted[i] = _epj_org[adr];
        }
    #endif //SUNWAY

#if 0
    Comm::barrier();
    if (Comm::getRank() == 0) std::cout << "Checking _epj_sorted in mortonSortLocalTreeOnlyImpl()." << std::endl;
    S32 num_err=0;
    for (S32 i=0; i<n_loc_tot_; i++) {
    	F64vec pos = _epj_sorted[i].getPos();
        if (fabs(pos.x) < 1.0e-300) num_err++;
    }
    if (num_err > 0) 
        std::cout << "num_err = " << num_err << " "
                  << "(RANK: " << Comm::getRank() << ")"
                  << std::endl;
    Finalize();
    std::exit(0);
#endif

    #ifdef REDUCE_MEMORY
        S32 n_cnt_ep = 0;
        S32 n_cnt_sp = 0;
        for(S32 i=0; i<n_loc_tot_; i++){
            const U32 adr = tp_loc_[i].adr_ptcl_;
            if( GetMSB(adr) == 0){
                tp_loc_[i].adr_ptcl_ = n_cnt_ep;
                n_cnt_ep++;
            }
            else{
                const U32 adr_new = ClearMSB(adr);
                tp_loc_[i].adr_ptcl_ = SetMSB((U32)n_cnt_sp);
                n_cnt_sp++;
            }
        }
    #endif
        
        time_profile_.morton_sort_local_tree += GetWtime() - time_offset;
        time_profile_.make_local_tree += GetWtime() - time_offset;
        //std::cerr<<"end of func"<<std::endl;        
    }


#ifndef USE_SUPER_DOMAIN
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTree3(const DomainInfo & dinfo,
                                const bool reuse_list){
        const F64 time_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
    #ifdef PHI_R_TREE
        const F64 PI = 4.0*atan(1.0);
    #endif
        tree_outer_pos_.resizeNoInitialize(n_proc);
        if(!reuse_list){
            comm_table_.n_ep_send_[my_rank] = comm_table_.n_sp_send_[my_rank]
                = comm_table_.n_ep_recv_[my_rank] = comm_table_.n_sp_recv_[my_rank] = 0;
            const F64 time_offset_tmp = GetWtime();
            F64ort outer_pos = tc_loc_[0].mom_.getVertexOut();
            Comm::allGather(&outer_pos, 1, tree_outer_pos_.getPointer());
            //Comm::allGather(&(tc_loc_[0].mom_.getVertexOut()), 1, tree_outer_pos_.getPointer());
            wtime_ex_let_allgather_root_ex += GetWtime() - time_offset_tmp;
    #ifdef DEBUG_PRINT_EX_LET_3
            if(Comm::getRank() == 0){
                for(S32 ib=0; ib<n_proc; ib++){
                    std::cerr<<"ib= "<<ib
                             <<" tree_outer_pos_[ib]= "<<tree_outer_pos_[ib]
                             <<" dinfo.getPosDomain(ib)= "<<dinfo.getPosDomain(ib)
                             <<std::endl;
                }
            }
    #endif
        }
    #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
    #endif
        {
            const S32 ith = Comm::getThreadNum();
            F64 time_offset0 = GetWtime();
            if(!reuse_list){
    #ifdef SEND_ALL_LET
                const F64 r_crit_sq = (length_ * length_) / (0.001 * 0.001);
    #else
                const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
    #endif
                const F64 longet_len_of_outer_box = tree_outer_pos_[my_rank].getFullLength().getMax();

    #ifdef SEND_ALL_LET
                const F64 r_crit_sq_2 = (longet_len_of_outer_box * longet_len_of_outer_box) / (0.001 * 0.001);
    #else
                const F64 r_crit_sq_2 = (longet_len_of_outer_box * longet_len_of_outer_box) / (theta_ * theta_);
    #endif
                rank_proc_send_[ith].resizeNoInitialize(0);
                id_ep_send_buf_[ith].resizeNoInitialize(0);
                id_sp_send_buf_[ith].resizeNoInitialize(0);
    #ifndef PHI_R_TREE
                ReallocatableArray<Tepj> ep_list_dummy;
    #endif
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
    #ifdef PHI_R_TREE
                    if(GetDistanceMinSqPeriodicX(tree_outer_pos_[ib], tree_outer_pos_[my_rank], 2.0*PI) > r_crit_sq_2){
                        id_sp_send_buf_[ith].push_back(0); // 0 means spj of root.
                    }
                    else{
                        const F64ort pos_target_domain = dinfo.getPosDomain(ib);
#ifdef DEBUG_PRINT_EX_LET_3
                        if(my_rank==0){
                            std::cerr<<"ib= "<<ib
                                     <<" dinfo.getPosDomain(ib)= "<<dinfo.getPosDomain(ib)
                                     <<" dinfo.getPosDomain(my_rank)= "<<dinfo.getPosDomain(my_rank)
                                     <<" tree_outer_pos_[ib]= "<<tree_outer_pos_[ib]
                                     <<" tree_outer_pos_[my_rank]= "<<tree_outer_pos_[my_rank]
                                     <<std::endl;
                        }
#endif
                        MakeListUsingTreeRecursiveTop<TSM, TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj>
                            (tc_loc_, 0, tp_loc_, epj_sorted_loc_, id_ep_send_buf_[ith],
                             sp_list_dummy, id_sp_send_buf_[ith], pos_target_domain,
                             r_crit_sq, n_leaf_limit_, 0, 2.0*PI);
                    }
    #else
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
    #endif
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
        } // end of OMP parallel scope

    #ifdef DEBUG_PRINT_EX_LET_3
        if(Comm::getRank() == 1){
            std::cerr<<"func: "<<__FUNCTION__
                     <<" line: "<<__LINE__
                     <<std::endl;
            for(S32 i=0; i<n_proc; i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.n_ep_send_[i]= "<<comm_table_.n_ep_send_[i]
                         <<" comm_table_.n_sp_send_[i]= "<<comm_table_.n_sp_send_[i]
                         <<std::endl;
            }
        }
        //exit(1);
    #endif

    #if 1
        // new
        comm_table_.setSPTop(spj_sorted_loc_[0]);
        comm_table_.exchangeLetSunWay(tree_outer_pos_, theta_, reuse_list, &dinfo);
    #else
        // original
        comm_table_.exchangeLet(reuse_list);
    #endif
        time_profile_.exchange_LET_tot += GetWtime() - time_offset;
    #ifdef DEBUG_PRINT_EX_LET_3
        /*
        PS::F64 cm_mass_tmp = 0.0;
        PS::F64vec cm_pos_tmp = 0.0;
        for(S32 i=0; i<comm_table_.ep_recv_.size(); i++){
            cm_mass_tmp += comm_table_.ep_recv_[i].mass;
            cm_pos_tmp  += comm_table_.ep_recv_[i].mass*comm_table_.ep_recv_[i].pos;
        }
        */
        if(Comm::getRank() == 1){
            for(S32 i=0; i<n_proc; i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.n_ep_send_[i]= "<<comm_table_.n_ep_send_[i]
                         <<" comm_table_.n_sp_send_[i]= "<<comm_table_.n_sp_send_[i]
                         <<" comm_table_.n_ep_recv_[i]= "<<comm_table_.n_ep_recv_[i]
                         <<" comm_table_.n_sp_recv_[i]= "<<comm_table_.n_sp_recv_[i]
                         <<std::endl;
            }
        }
        //exit(1);
    #endif
    }
#endif //(not) USE_SUPER_DOMAIN


#ifdef USE_SUPER_DOMAIN
        template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeSuperDomain(const DomainInfo & dinfo,
                                          const bool reuse_list){
        const F64 time_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 n_proc_1d = dinfo.getNDomain(0);
        const S32 n_proc_sub = n_proc / n_proc_1d;
        const S32 my_rank = Comm::getRank();
        const S32 my_rank_1d = dinfo.getRank1d(0);
        const S32 my_rank_sub = dinfo.getRankSub(0);
        const F64 PI = 4.0*atan(1.0);
        static ReallocatableArray<F64ort> outer_pos_super_domain;
        static ReallocatableArray<F64ort> outer_pos_sub_domain;
        MPI_Comm comm_sub = dinfo.getCommSub(0);
        MPI_Comm comm_1d = dinfo.getComm1d(0);
        outer_pos_super_domain.resizeNoInitialize(n_proc_1d);
        outer_pos_sub_domain.resizeNoInitialize(n_proc_sub);
        //////////////////////
        // initialize
        if(!reuse_list){
            ///////////////////////////
            // exchange outer boundary of super domain
            F64ort my_outer_pos_sub_domain = tc_loc_[0].mom_.getVertexOut();
            outer_pos_sub_domain.resizeNoInitialize(n_proc_sub);
            MPI_Allgather(&my_outer_pos_sub_domain, 1, GetDataType<F64ort>(),
                          outer_pos_sub_domain.getPointer(), 1, GetDataType<F64ort>(),
                          comm_sub);
            F64ort my_outer_pos_super_domain;
            for(S32 i=0; i<n_proc_sub; i++) my_outer_pos_super_domain.merge(outer_pos_sub_domain[i]);
            MPI_Allgather(&my_outer_pos_super_domain, 1, GetDataType<F64ort>(),
                          outer_pos_super_domain.getPointer(), 1, GetDataType<F64ort>(),
                          comm_1d);
            // exchange outer boundary of super domain
            ///////////////////////////
    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN
            /*
            Comm::barrier();
            if(my_rank == 0){
                for(S32 i=0; i<n_proc_1d; i++){
                    std::cerr<<"outer_pos_super_domain[i]= "<<outer_pos_super_domain[i]<<std::endl;
                }
            }
            //Comm::barrier();
            //exit(1);
            */
    #endif            
            const F64 my_longest_len_super_domain = outer_pos_super_domain[my_rank_1d].getFullLength().getMax();
            const F64 r_crit_sq_send_super_domain = (my_longest_len_super_domain * my_longest_len_super_domain) / (theta_ * theta_);
            const F64 r_crit_sq_tree_walk = (length_ * length_) / (theta_ * theta_);
            const F64 len_peri_x = 2.0*PI;
            comm_table_.clearSize();
            ReallocatableArray<Tspj> sp_list_dummy;
            S32 n_ep_cum = 0;
            S32 n_sp_cum = 0;
            S32 n_ep_cum_2 = 0;
            S32 n_sp_cum_2 = 0;
            F64 time_offset0 = GetWtime();
            for(S32 ib=0; ib<n_proc_1d; ib++){
                // send ptcl
                if(GetDistanceMinSqPeriodicX(outer_pos_super_domain[ib], my_outer_pos_super_domain, len_peri_x) > r_crit_sq_send_super_domain){
                    comm_table_.rank_sd_sd_send_.push_back(ib);
                }
                else{
                    if(my_rank_1d == ib){
                        for(S32 jb=0; jb<n_proc_sub; jb++){
                            const S32 rank_target = ib*n_proc_sub + jb;
                            if(my_rank == rank_target) continue;
                            const F64ort pos_target_domain = dinfo.getPosDomain(rank_target);
                            MakeListUsingTreeRecursiveTop<TSM, TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj>
                                (tc_loc_, 0, tp_loc_,
                                 epj_sorted_loc_, comm_table_.adr_ep_domain_send_,
                                 sp_list_dummy,   comm_table_.adr_sp_domain_send_,
                                 pos_target_domain, r_crit_sq_tree_walk, n_leaf_limit_, 0, len_peri_x);
                            comm_table_.rank_ptcl_domain_send_.push_back(rank_target);
                            comm_table_.n_ep_domain_send_.push_back(comm_table_.adr_ep_domain_send_.size() - n_ep_cum);
                            comm_table_.n_sp_domain_send_.push_back(comm_table_.adr_sp_domain_send_.size() - n_sp_cum);
                            n_ep_cum = comm_table_.adr_ep_domain_send_.size();
                            n_sp_cum = comm_table_.adr_sp_domain_send_.size();
                        }
                    }
                    else{
                        const F64ort pos_target_domain = outer_pos_super_domain[ib];
                        MakeListUsingTreeRecursiveTop<TSM, TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj>
                            (tc_loc_, 0, tp_loc_,
                             epj_sorted_loc_, comm_table_.adr_ep_sd_send_,
                             sp_list_dummy,   comm_table_.adr_sp_sd_send_,
                             pos_target_domain, r_crit_sq_tree_walk, n_leaf_limit_, 0, len_peri_x);
                            comm_table_.rank_ptcl_sd_send_.push_back(ib);
                            comm_table_.n_ep_sd_send_.push_back(comm_table_.adr_ep_sd_send_.size() - n_ep_cum_2);
                            comm_table_.n_sp_sd_send_.push_back(comm_table_.adr_sp_sd_send_.size() - n_sp_cum_2);
                            n_ep_cum_2 = comm_table_.adr_ep_sd_send_.size();
                            n_sp_cum_2 = comm_table_.adr_sp_sd_send_.size();
                    }
                }
                // recv ptcl
                const F64 target_longest_len_super_domain = outer_pos_super_domain[ib].getFullLength().getMax();
                const F64 r_crit_sq_recv_super_domain = (target_longest_len_super_domain * target_longest_len_super_domain) / (theta_ * theta_);
                if(GetDistanceMinSqPeriodicX(my_outer_pos_super_domain, outer_pos_super_domain[ib], len_peri_x) > r_crit_sq_recv_super_domain){
                    comm_table_.rank_sd_sd_recv_.push_back(ib);
                }
                else{
                    if(my_rank_1d == ib){
                        for(S32 jb=0; jb<n_proc_sub; jb++){
                            const S32 rank_target = ib*n_proc_sub + jb;
                            if(my_rank == rank_target) continue;
                            comm_table_.rank_ptcl_domain_recv_.push_back(rank_target);
                        }
                    }
                    else{
                        comm_table_.rank_ptcl_sd_recv_.push_back(ib);
                    }
                }
            }
            time_profile_.make_LET_1st += GetWtime() - time_offset0;
    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN
            Comm::barrier();
            if(my_rank==0){
                for(S32 i=0; i<comm_table_.rank_sd_sd_send_.size(); i++){
                    std::cerr<<"comm_table_.rank_sd_sd_send_[i]= "<<comm_table_.rank_sd_sd_send_[i]<<std::endl;
                }
                for(S32 i=0; i<comm_table_.rank_ptcl_domain_send_.size(); i++){
                    std::cerr<<"comm_table_.rank_ptcl_domain_send_[i]= "<<comm_table_.rank_ptcl_domain_send_[i]
                             <<" comm_table_.n_ep_domain_send_[i]= "<<comm_table_.n_ep_domain_send_[i]
                             <<" comm_table_.n_sp_domain_send_[i]= "<<comm_table_.n_sp_domain_send_[i]
                             <<std::endl;
                }
                for(S32 i=0; i<comm_table_.rank_ptcl_sd_send_.size(); i++){
                    std::cerr<<"comm_table_.rank_ptcl_sd_send_[i]= "<<comm_table_.rank_ptcl_sd_send_[i]
                             <<" comm_table_.n_ep_sd_send_[i]= "<<comm_table_.n_ep_sd_send_[i]
                             <<" comm_table_.n_sp_sd_send_[i]= "<<comm_table_.n_sp_sd_send_[i]
                             <<std::endl;
                }
            }
            //Comm::barrier();
            //exit(1);
    #endif
            
            ///////////////////
            // exchange # of ep and sp (domain -> domain) [in the same super domain]
            comm_table_.n_ep_domain_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size());
            comm_table_.n_sp_domain_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size());
            comm_table_.req_ep_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size());
            comm_table_.req_sp_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size());
            comm_table_.stat_ep_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size());
            comm_table_.stat_sp_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size());
            comm_table_.req_ep_send_.resizeNoInitialize(comm_table_.rank_ptcl_domain_send_.size());
            comm_table_.req_sp_send_.resizeNoInitialize(comm_table_.rank_ptcl_domain_send_.size());
            comm_table_.stat_ep_send_.resizeNoInitialize(comm_table_.rank_ptcl_domain_send_.size());
            comm_table_.stat_sp_send_.resizeNoInitialize(comm_table_.rank_ptcl_domain_send_.size());
            for(S32 i=0; i<comm_table_.rank_ptcl_domain_send_.size(); i++){
                MPI_Isend( comm_table_.n_ep_domain_send_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_domain_send_[i], 0, MPI_COMM_WORLD, comm_table_.req_ep_send_.getPointer(i));
                MPI_Isend( comm_table_.n_sp_domain_send_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_domain_send_[i], 1, MPI_COMM_WORLD, comm_table_.req_sp_send_.getPointer(i));
            }
            for(S32 i=0; i<comm_table_.rank_ptcl_domain_recv_.size(); i++){
                MPI_Irecv( comm_table_.n_ep_domain_recv_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_domain_recv_[i], 0, MPI_COMM_WORLD, comm_table_.req_ep_recv_.getPointer(i));
                MPI_Irecv( comm_table_.n_sp_domain_recv_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_domain_recv_[i], 1, MPI_COMM_WORLD, comm_table_.req_sp_recv_.getPointer(i));
            }            
            MPI_Waitall(comm_table_.rank_ptcl_domain_recv_.size(), comm_table_.req_ep_recv_.getPointer(), comm_table_.stat_ep_recv_.getPointer());
            MPI_Waitall(comm_table_.rank_ptcl_domain_send_.size(), comm_table_.req_ep_send_.getPointer(), comm_table_.stat_ep_send_.getPointer());
            MPI_Waitall(comm_table_.rank_ptcl_domain_recv_.size(), comm_table_.req_sp_recv_.getPointer(), comm_table_.stat_sp_recv_.getPointer());
            MPI_Waitall(comm_table_.rank_ptcl_domain_send_.size(), comm_table_.req_sp_send_.getPointer(), comm_table_.stat_sp_send_.getPointer());

            comm_table_.n_disp_ep_domain_send_.resizeNoInitialize(comm_table_.rank_ptcl_domain_send_.size()+1);
            comm_table_.n_disp_sp_domain_send_.resizeNoInitialize(comm_table_.rank_ptcl_domain_send_.size()+1);
            comm_table_.n_disp_ep_domain_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size()+1);
            comm_table_.n_disp_sp_domain_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size()+1);
            comm_table_.n_disp_ep_domain_send_[0] = comm_table_.n_disp_sp_domain_send_[0] = comm_table_.n_disp_ep_domain_recv_[0] = comm_table_.n_disp_sp_domain_recv_[0] = 0;
            for(S32 i=0; i<comm_table_.rank_ptcl_domain_send_.size(); i++){
                comm_table_.n_disp_ep_domain_send_[i+1] = comm_table_.n_disp_ep_domain_send_[i] + comm_table_.n_ep_domain_send_[i];
                comm_table_.n_disp_sp_domain_send_[i+1] = comm_table_.n_disp_sp_domain_send_[i] + comm_table_.n_sp_domain_send_[i];
            }
            for(S32 i=0; i<comm_table_.rank_ptcl_domain_recv_.size(); i++){
                comm_table_.n_disp_ep_domain_recv_[i+1] = comm_table_.n_disp_ep_domain_recv_[i] + comm_table_.n_ep_domain_recv_[i];
                comm_table_.n_disp_sp_domain_recv_[i+1] = comm_table_.n_disp_sp_domain_recv_[i] + comm_table_.n_sp_domain_recv_[i];
            }
            // exchange # of ep and sp (domain -> domain)
            ///////////////////

            ///////////////////
            // exchange # of ep and sp (domain -> super domain)
            comm_table_.n_ep_sd_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size());
            comm_table_.n_sp_sd_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size());
            comm_table_.req_ep_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size());
            comm_table_.req_sp_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size());
            comm_table_.stat_ep_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size());
            comm_table_.stat_sp_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size());
            comm_table_.req_ep_send_.resizeNoInitialize(comm_table_.rank_ptcl_sd_send_.size());
            comm_table_.req_sp_send_.resizeNoInitialize(comm_table_.rank_ptcl_sd_send_.size());
            comm_table_.stat_ep_send_.resizeNoInitialize(comm_table_.rank_ptcl_sd_send_.size());
            comm_table_.stat_sp_send_.resizeNoInitialize(comm_table_.rank_ptcl_sd_send_.size());            
            for(S32 i=0; i<comm_table_.rank_ptcl_sd_send_.size(); i++){
                MPI_Isend( comm_table_.n_ep_sd_send_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_sd_send_[i], 0, comm_1d, comm_table_.req_ep_send_.getPointer(i));
                MPI_Isend( comm_table_.n_sp_sd_send_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_sd_send_[i], 1, comm_1d, comm_table_.req_sp_send_.getPointer(i));
            }
            for(S32 i=0; i<comm_table_.rank_ptcl_sd_recv_.size(); i++){
                MPI_Irecv( comm_table_.n_ep_sd_recv_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_sd_recv_[i], 0, comm_1d, comm_table_.req_ep_recv_.getPointer(i));
                MPI_Irecv( comm_table_.n_sp_sd_recv_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_sd_recv_[i], 1, comm_1d, comm_table_.req_sp_recv_.getPointer(i));
            }            
            MPI_Waitall(comm_table_.rank_ptcl_sd_recv_.size(), comm_table_.req_ep_recv_.getPointer(), comm_table_.stat_ep_recv_.getPointer());
            MPI_Waitall(comm_table_.rank_ptcl_sd_send_.size(), comm_table_.req_ep_send_.getPointer(), comm_table_.stat_ep_send_.getPointer());
            MPI_Waitall(comm_table_.rank_ptcl_sd_recv_.size(), comm_table_.req_sp_recv_.getPointer(), comm_table_.stat_sp_recv_.getPointer());
            MPI_Waitall(comm_table_.rank_ptcl_sd_send_.size(), comm_table_.req_sp_send_.getPointer(), comm_table_.stat_sp_send_.getPointer());

            comm_table_.n_disp_ep_sd_send_.resizeNoInitialize(comm_table_.rank_ptcl_sd_send_.size()+1);
            comm_table_.n_disp_sp_sd_send_.resizeNoInitialize(comm_table_.rank_ptcl_sd_send_.size()+1);
            comm_table_.n_disp_ep_sd_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size()+1);
            comm_table_.n_disp_sp_sd_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size()+1);
            comm_table_.n_disp_ep_sd_send_[0] = comm_table_.n_disp_sp_sd_send_[0] = comm_table_.n_disp_ep_sd_recv_[0] = comm_table_.n_disp_sp_sd_recv_[0] = 0;
            for(S32 i=0; i<comm_table_.rank_ptcl_sd_send_.size(); i++){
                comm_table_.n_disp_ep_sd_send_[i+1] = comm_table_.n_disp_ep_sd_send_[i] + comm_table_.n_ep_sd_send_[i];
                comm_table_.n_disp_sp_sd_send_[i+1] = comm_table_.n_disp_sp_sd_send_[i] + comm_table_.n_sp_sd_send_[i];
            }
            for(S32 i=0; i<comm_table_.rank_ptcl_sd_recv_.size(); i++){
                comm_table_.n_disp_ep_sd_recv_[i+1] = comm_table_.n_disp_ep_sd_recv_[i] + comm_table_.n_ep_sd_recv_[i];
                comm_table_.n_disp_sp_sd_recv_[i+1] = comm_table_.n_disp_sp_sd_recv_[i] + comm_table_.n_sp_sd_recv_[i];
            }
            // exchange # of ep and sp (domain -> super domain)
            ///////////////////

            ///////////////////
            // exchange # of ep and sp (domain -> super domain) in super domain
            comm_table_.n_ep_domain_recv_2_.resizeNoInitialize(n_proc_sub);
            comm_table_.n_sp_domain_recv_2_.resizeNoInitialize(n_proc_sub);
            comm_table_.n_disp_ep_domain_recv_2_.resizeNoInitialize(n_proc_sub+1);
            comm_table_.n_disp_sp_domain_recv_2_.resizeNoInitialize(n_proc_sub+1);
            comm_table_.n_ep_domain_send_2_ = comm_table_.n_disp_ep_sd_recv_[ comm_table_.rank_ptcl_sd_recv_.size() ];
            MPI_Allgather( &(comm_table_.n_ep_domain_send_2_), 1, GetDataType<S32>(),
                           comm_table_.n_ep_domain_recv_2_.getPointer(), 1, GetDataType<S32>(),
                           comm_sub);
            comm_table_.n_sp_domain_send_2_ = comm_table_.n_disp_sp_sd_recv_[ comm_table_.rank_ptcl_sd_recv_.size() ];
            MPI_Allgather( &(comm_table_.n_sp_domain_send_2_), 1, GetDataType<S32>(),
                           comm_table_.n_sp_domain_recv_2_.getPointer(), 1, GetDataType<S32>(),
                           comm_sub);
            comm_table_.n_disp_ep_domain_recv_2_[0] = comm_table_.n_disp_sp_domain_recv_2_[0] = 0;
            comm_table_.n_ep_domain_recv_2_[my_rank_sub] = comm_table_.n_sp_domain_recv_2_[my_rank_sub] = 0;
            for(S32 i=0; i<n_proc_sub; i++){
                comm_table_.n_disp_ep_domain_recv_2_[i+1] = comm_table_.n_disp_ep_domain_recv_2_[i] + comm_table_.n_ep_domain_recv_2_[i];
                comm_table_.n_disp_sp_domain_recv_2_[i+1] = comm_table_.n_disp_sp_domain_recv_2_[i] + comm_table_.n_sp_domain_recv_2_[i];
            }
            // exchange # of ep and sp (domain -> super domain) in super domain
            ///////////////////
#ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN
            if(my_rank == 0){
                std::cerr<<std::endl;
                comm_table_.dumpNumber();
                std::cerr<<std::endl;
            }
#endif
            //time_profile_.make_LET_2nd += GetWtime() - time_offset0;            
        } // end of reuse
        // initialize
        //////////////////////

        ///////////////////////////
        // allgather SPJ of super domain
        static ReallocatableArray<Tspj> spj_sub_domain;
        spj_sub_domain.resizeNoInitialize(n_proc_sub);
        MPI_Allgather(&spj_sorted_loc_[0], 1, GetDataType<Tspj>(),
                      spj_sub_domain.getPointer(), 1, GetDataType<Tspj>(),
                      comm_sub);
#ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN
        if(my_rank == 0){
            std::cerr<<std::endl;
            for(S32 i=0; i<n_proc_sub; i++){
                std::cerr<<"i="<<i<<std::endl;
                spj_sub_domain[i].dump();
            }
            std::cerr<<std::endl;
        }
#endif
        Tspj my_spj_super_domain;
        my_spj_super_domain.clear();
        for(S32 i=0; i<n_proc_sub; i++){
            my_spj_super_domain.mass += spj_sub_domain[i].mass;
            my_spj_super_domain.pos += spj_sub_domain[i].mass * spj_sub_domain[i].pos;
            my_spj_super_domain.pos_phi += spj_sub_domain[i].mass * spj_sub_domain[i].pos_phi;
            my_spj_super_domain.pos_r += spj_sub_domain[i].mass * spj_sub_domain[i].pos_r;
        }
        my_spj_super_domain.pos     /= my_spj_super_domain.mass;
        my_spj_super_domain.pos_phi /= my_spj_super_domain.mass;
        my_spj_super_domain.pos_r   /= my_spj_super_domain.mass;
#ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"super domain"<<std::endl;
            my_spj_super_domain.dump();
        }
#endif
        static ReallocatableArray<Tspj> spj_super_domain;
        spj_super_domain.resizeNoInitialize(n_proc_1d);
        MPI_Allgather(&my_spj_super_domain, 1, GetDataType<Tspj>(),
                      spj_super_domain.getPointer(), 1, GetDataType<Tspj>(),
                      comm_1d);
#ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"spj super domain"<<std::endl;
            for(S32 i=0; i<n_proc_1d; i++){
                std::cerr<<"i= "<<i<<" n_proc_1d= "<<n_proc_1d<<std::endl;
                spj_super_domain[i].dump();
            }
        }
#endif
        comm_table_.sp_recv_.reserve(comm_table_.rank_sd_sd_recv_.size());
        comm_table_.sp_recv_.clearSize();
        for(S32 i=0; i<comm_table_.rank_sd_sd_recv_.size(); i++){
            S32 rank = comm_table_.rank_sd_sd_recv_[i];
            Tspj spj = spj_super_domain[rank];
            comm_table_.sp_recv_.pushBackNoCheck(spj);
        }
#ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN
        Comm::barrier();
        if(my_rank == 0){
            for(S32 i=0; i<comm_table_.rank_sd_sd_recv_.size(); i++){
                std::cerr<<"comm_table_.rank_sd_sd_recv_[i]= "<<comm_table_.rank_sd_sd_recv_[i]<<std::endl;
                std::cerr<<"comm_table_.sp_recv_[i].mass= "<<comm_table_.sp_recv_[i].mass<<std::endl;
                std::cerr<<"comm_table_.sp_recv_[i].pos_phi= "<<comm_table_.sp_recv_[i].pos_phi<<std::endl;
                std::cerr<<"comm_table_.sp_recv_[i].pos_r= "<<comm_table_.sp_recv_[i].pos_r<<std::endl;
                std::cerr<<std::endl;
            }
        }
        //exit(1);
#endif
        S32 adr_sp_recv_offset = comm_table_.sp_recv_.size();
        // allgather SPJ of super domain
        ///////////////////////////

        ////////////////////////////////////
        // exchange EPJ and SPJ of super domain among super domains
        S32 n_proc_ep_send = 0;
        S32 n_proc_sp_send = 0;
        S32 n_proc_ep_recv = 0;
        S32 n_proc_sp_recv = 0;
        comm_table_.ep_send_.resizeNoInitialize( comm_table_.n_disp_ep_sd_send_[comm_table_.rank_ptcl_sd_send_.size()] );
        comm_table_.sp_send_.resizeNoInitialize( comm_table_.n_disp_sp_sd_send_[comm_table_.rank_ptcl_sd_send_.size()] );
        comm_table_.ep_recv_.resizeNoInitialize( comm_table_.n_disp_ep_sd_recv_[comm_table_.rank_ptcl_sd_recv_.size()] );
        comm_table_.sp_recv_.resizeNoInitialize( comm_table_.n_disp_sp_sd_recv_[comm_table_.rank_ptcl_sd_recv_.size()] + adr_sp_recv_offset);
        for(S32 i=0; i<comm_table_.rank_ptcl_sd_send_.size(); i++){
            S32 rank_send = comm_table_.rank_ptcl_sd_send_[i];
            if( comm_table_.n_ep_sd_send_[i] > 0){
                S32 n_ep = comm_table_.n_ep_sd_send_[i];
                S32 adr_ep_head = comm_table_.n_disp_ep_sd_send_[i];
                S32 adr_ep_tail = comm_table_.n_disp_ep_sd_send_[i+1];
                for(S32 ip=adr_ep_head; ip<adr_ep_tail; ip++){
                    comm_table_.ep_send_[ip] = epj_sorted_loc_[ comm_table_.adr_ep_sd_send_[ip] ];
                }
                MPI_Isend( comm_table_.ep_send_.getPointer(adr_ep_head), n_ep, GetDataType<Tepj>(), rank_send, 0, comm_1d, comm_table_.req_ep_send_.getPointer(n_proc_ep_send));
                n_proc_ep_send++;
            }
            if( comm_table_.n_sp_sd_send_[i] > 0){
                S32 n_sp = comm_table_.n_sp_sd_send_[i];
                S32 adr_sp_head = comm_table_.n_disp_sp_sd_send_[i];
                S32 adr_sp_tail = comm_table_.n_disp_sp_sd_send_[i+1];
                for(S32 ip=adr_sp_head; ip<adr_sp_tail; ip++){
                    comm_table_.sp_send_[ip] = spj_sorted_loc_[ comm_table_.adr_sp_sd_send_[ip] ];
                }
                MPI_Isend( comm_table_.sp_send_.getPointer(adr_sp_head), n_sp, GetDataType<Tspj>(), rank_send, 0, comm_1d, comm_table_.req_sp_send_.getPointer(n_proc_sp_send));
                n_proc_sp_send++;
            }            
        }
        for(S32 i=0; i<comm_table_.rank_ptcl_sd_recv_.size(); i++){
            S32 rank_recv = comm_table_.rank_ptcl_sd_recv_[i];
            if( comm_table_.n_ep_sd_recv_[i] > 0){
                S32 n_ep = comm_table_.n_ep_sd_recv_[i];
                S32 adr_ep_head = comm_table_.n_disp_ep_sd_recv_[i];
                MPI_Irecv( comm_table_.ep_recv_.getPointer(adr_ep_head), n_ep, GetDataType<Tepj>(), rank_recv, 0, comm_1d, comm_table_.req_ep_recv_.getPointer(n_proc_ep_recv));
                n_proc_ep_recv++;
            }
            if( comm_table_.n_sp_sd_recv_[i] > 0){
                S32 n_sp = comm_table_.n_sp_sd_recv_[i];
                S32 adr_sp_head = comm_table_.n_disp_sp_sd_recv_[i] + adr_sp_recv_offset;
                MPI_Irecv( comm_table_.sp_recv_.getPointer(adr_sp_head), n_sp, GetDataType<Tspj>(), rank_recv, 0, comm_1d, comm_table_.req_sp_recv_.getPointer(n_proc_sp_recv));
                n_proc_sp_recv++;
            }
        }
        MPI_Waitall(n_proc_ep_recv, comm_table_.req_ep_recv_.getPointer(), comm_table_.stat_ep_recv_.getPointer());
        MPI_Waitall(n_proc_ep_send, comm_table_.req_ep_send_.getPointer(), comm_table_.stat_ep_send_.getPointer());
        MPI_Waitall(n_proc_sp_recv, comm_table_.req_sp_recv_.getPointer(), comm_table_.stat_sp_recv_.getPointer());
        MPI_Waitall(n_proc_sp_send, comm_table_.req_sp_send_.getPointer(), comm_table_.stat_sp_send_.getPointer());

        S32 adr_sp_recv_offset_2 = comm_table_.n_disp_sp_sd_recv_[comm_table_.rank_ptcl_sd_recv_.size()] + adr_sp_recv_offset;
        S32 adr_ep_recv_offset   = comm_table_.n_disp_ep_sd_recv_[comm_table_.rank_ptcl_sd_recv_.size()];

    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN
/*
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"OK 01 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            std::cerr<<"comm_table_.ep_recv_.size()= "<<comm_table_.ep_recv_.size()
                     <<" comm_table_.sp_recv_.size()= "<<comm_table_.sp_recv_.size()
                     <<std::endl;
            for(S32 i=0; i<adr_ep_recv_offset; i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.ep_recv_[i].mass= "<<comm_table_.ep_recv_[i].mass
                         <<" comm_table_.ep_recv_[i].pos= "<<comm_table_.ep_recv_[i].pos
                         <<std::endl;
            }
            for(S32 i=0; i<adr_sp_recv_offset_2; i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.sp_recv_[i].mass= "<<comm_table_.sp_recv_[i].mass
                         <<" comm_table_.sp_recv_[i].pos= "<<comm_table_.sp_recv_[i].pos
                         <<std::endl;
            }
        }
        Comm::barrier();
        //exit(1);
        */
    #endif

        // in subdomain
        comm_table_.ep_recv_.resizeNoInitialize(adr_ep_recv_offset   + comm_table_.n_disp_ep_domain_recv_2_[n_proc_sub]);
        comm_table_.sp_recv_.resizeNoInitialize(adr_sp_recv_offset_2 + comm_table_.n_disp_sp_domain_recv_2_[n_proc_sub]);

    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN
     /*
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"OK 02 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            std::cerr<<"comm_table_.ep_recv_.size()= "<<comm_table_.ep_recv_.size()
                     <<" comm_table_.sp_recv_.size()= "<<comm_table_.sp_recv_.size()
                     <<std::endl;
            std::cerr<<"comm_table_.n_ep_domain_send_2_= "<<comm_table_.n_ep_domain_send_2_
                     <<" comm_table_.n_sp_domain_send_2_= "<<comm_table_.n_sp_domain_send_2_
                     <<std::endl;
            std::cerr<<"adr_sp_recv_offset= "<<adr_sp_recv_offset<<std::endl;
            for(S32 i=0; i<comm_table_.ep_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.ep_recv_[i].mass= "<<comm_table_.ep_recv_[i].mass
                         <<" comm_table_.ep_recv_[i].pos= "<<comm_table_.ep_recv_[i].pos
                         <<std::endl;
            }
            for(S32 i=0; i<comm_table_.sp_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.sp_recv_[i].mass= "<<comm_table_.sp_recv_[i].mass
                         <<" comm_table_.sp_recv_[i].pos= "<<comm_table_.sp_recv_[i].pos
                         <<std::endl;
            }
        }
        Comm::barrier();
        //exit(1);
        */
    #endif

        n_proc_ep_send = n_proc_sp_send = n_proc_ep_recv = n_proc_sp_recv = 0;
        for(S32 i=0; i<n_proc_sub; i++){
            if(i==my_rank_sub) continue;
            S32 n_ep = comm_table_.n_ep_domain_send_2_;
            S32 n_sp = comm_table_.n_sp_domain_send_2_;
            MPI_Isend( comm_table_.ep_recv_.getPointer(), n_ep, GetDataType<Tepj>(), i, 0, comm_sub, comm_table_.req_ep_send_.getPointer(n_proc_ep_send));
            MPI_Isend( comm_table_.sp_recv_.getPointer(adr_sp_recv_offset), n_sp, GetDataType<Tspj>(), i, 1, comm_sub, comm_table_.req_sp_send_.getPointer(n_proc_sp_send));
            n_proc_ep_send++;
            n_proc_sp_send++;
        }
        for(S32 i=0; i<n_proc_sub; i++){
            if(i==my_rank_sub) continue;
            S32 n_ep = comm_table_.n_ep_domain_recv_2_[i];
            S32 n_sp = comm_table_.n_sp_domain_recv_2_[i];
            MPI_Irecv( comm_table_.ep_recv_.getPointer(comm_table_.n_disp_ep_domain_recv_2_[i]+adr_ep_recv_offset), n_ep, GetDataType<Tepj>(), i, 0, comm_sub, comm_table_.req_ep_recv_.getPointer(n_proc_ep_recv));
            MPI_Irecv( comm_table_.sp_recv_.getPointer(comm_table_.n_disp_sp_domain_recv_2_[i]+adr_sp_recv_offset_2), n_sp, GetDataType<Tspj>(), i, 1, comm_sub, comm_table_.req_sp_recv_.getPointer(n_proc_sp_recv));
            n_proc_ep_recv++;
            n_proc_sp_recv++;
        }
        MPI_Waitall(n_proc_ep_recv, comm_table_.req_ep_recv_.getPointer(), comm_table_.stat_ep_recv_.getPointer());
        MPI_Waitall(n_proc_ep_send, comm_table_.req_ep_send_.getPointer(), comm_table_.stat_ep_send_.getPointer());
        MPI_Waitall(n_proc_sp_recv, comm_table_.req_sp_recv_.getPointer(), comm_table_.stat_sp_recv_.getPointer());
        MPI_Waitall(n_proc_sp_send, comm_table_.req_sp_send_.getPointer(), comm_table_.stat_sp_send_.getPointer());
        
    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN
        /*
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"OK 03 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            std::cerr<<"comm_table_.ep_recv_.size()= "<<comm_table_.ep_recv_.size()
                     <<" comm_table_.sp_recv_.size()= "<<comm_table_.sp_recv_.size()
                     <<std::endl;
            for(S32 i=0; i<comm_table_.ep_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.ep_recv_[i].mass= "<<comm_table_.ep_recv_[i].mass
                         <<" comm_table_.ep_recv_[i].pos= "<<comm_table_.ep_recv_[i].pos
                         <<std::endl;
            }
            for(S32 i=0; i<comm_table_.sp_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.sp_recv_[i].mass= "<<comm_table_.sp_recv_[i].mass
                         <<" comm_table_.sp_recv_[i].pos= "<<comm_table_.sp_recv_[i].pos
                         <<std::endl;
            }
        }
        Comm::barrier();
        //exit(1);
        */
    #endif

        S32 adr_ep_recv_offset_2 = comm_table_.n_disp_ep_domain_recv_2_[n_proc_sub] + adr_ep_recv_offset;
        S32 adr_sp_recv_offset_3 = comm_table_.n_disp_sp_domain_recv_2_[n_proc_sub] + adr_sp_recv_offset_2;
        
    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN
     /*
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"OK 04 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            std::cerr<<"comm_table_.ep_recv_.size()= "<<comm_table_.ep_recv_.size()
                     <<" comm_table_.sp_recv_.size()= "<<comm_table_.sp_recv_.size()
                     <<std::endl;
            std::cerr<<"adr_ep_recv_offset= "<<adr_ep_recv_offset
                     <<" adr_ep_recv_offset_2= "<<adr_ep_recv_offset_2
                     <<std::endl;
            std::cerr<<"adr_sp_recv_offset= "<<adr_sp_recv_offset
                     <<" adr_sp_recv_offset_2= "<<adr_sp_recv_offset_2
                     <<" adr_sp_recv_offset_3= "<<adr_sp_recv_offset_3
                     <<std::endl;
            for(S32 i=0; i<comm_table_.ep_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.ep_recv_[i].mass= "<<comm_table_.ep_recv_[i].mass
                         <<" comm_table_.ep_recv_[i].pos= "<<comm_table_.ep_recv_[i].pos
                         <<std::endl;
            }
            for(S32 i=0; i<comm_table_.sp_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.sp_recv_[i].mass= "<<comm_table_.sp_recv_[i].mass
                         <<" comm_table_.sp_recv_[i].pos= "<<comm_table_.sp_recv_[i].pos
                         <<std::endl;
            }
        }
        Comm::barrier();
        //exit(1);
        */
    #endif

        // exchange EPJ and SPJ of super domain
        ////////////////////////////////////

        ////////////////////////////////////////
        // exchange EPJ and SPJ in super domain
        n_proc_ep_send = n_proc_sp_send = n_proc_ep_recv = n_proc_sp_recv = 0;
        comm_table_.ep_send_.resizeNoInitialize(comm_table_.n_disp_ep_domain_send_[comm_table_.rank_ptcl_domain_send_.size()]);
        comm_table_.sp_send_.resizeNoInitialize(comm_table_.n_disp_sp_domain_send_[comm_table_.rank_ptcl_domain_send_.size()]);
        comm_table_.ep_recv_.resizeNoInitialize(comm_table_.n_disp_ep_domain_recv_[comm_table_.rank_ptcl_domain_recv_.size()] + adr_ep_recv_offset_2);
        comm_table_.sp_recv_.resizeNoInitialize(comm_table_.n_disp_sp_domain_recv_[comm_table_.rank_ptcl_domain_recv_.size()] + adr_sp_recv_offset_3);

    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN
     /*
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"OK 05 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            std::cerr<<"comm_table_.ep_send_.size()= "<<comm_table_.ep_send_.size()
                     <<" comm_table_.sp_send_.size()= "<<comm_table_.sp_send_.size()
                     <<std::endl;
            std::cerr<<"comm_table_.ep_recv_.size()= "<<comm_table_.ep_recv_.size()
                     <<" comm_table_.sp_recv_.size()= "<<comm_table_.sp_recv_.size()
                     <<std::endl;

            std::cerr<<"comm_table_.rank_ptcl_domain_send_.size()= "<<comm_table_.rank_ptcl_domain_send_.size()
                     <<" comm_table_.rank_ptcl_domain_recv_.size()= "<<comm_table_.rank_ptcl_domain_recv_.size()
                     <<std::endl;
            for(S32 i=0; i<comm_table_.rank_ptcl_domain_send_.size(); i++){
                std::cerr<<"comm_table_.rank_ptcl_domain_send_[i]= "<<comm_table_.rank_ptcl_domain_send_[i]
                         <<" comm_table_.n_ep_domain_send_[i]= "<<comm_table_.n_ep_domain_send_[i]
                         <<" comm_table_.n_sp_domain_send_[i]= "<<comm_table_.n_sp_domain_send_[i]
                         <<std::endl;
            }
            for(S32 i=0; i<comm_table_.rank_ptcl_domain_recv_.size(); i++){
                std::cerr<<"comm_table_.rank_ptcl_domain_recv_[i]= "<<comm_table_.rank_ptcl_domain_recv_[i]
                         <<" comm_table_.n_ep_domain_recv_[i]= "<<comm_table_.n_ep_domain_recv_[i]
                         <<" comm_table_.n_sp_domain_recv_[i]= "<<comm_table_.n_sp_domain_recv_[i]
                         <<std::endl;
            }
        }
        Comm::barrier();
        //exit(1);
        */
    #endif

        n_proc_ep_send = n_proc_sp_send = n_proc_ep_recv = n_proc_sp_recv = 0;
        for(S32 i=0; i<comm_table_.rank_ptcl_domain_send_.size(); i++){
            S32 rank_send = comm_table_.rank_ptcl_domain_send_[i];
            if( comm_table_.n_ep_domain_send_[i] > 0){
                S32 n_ep = comm_table_.n_ep_domain_send_[i];
                S32 adr_ep_head = comm_table_.n_disp_ep_domain_send_[i];
                S32 adr_ep_tail = comm_table_.n_disp_ep_domain_send_[i+1];
                for(S32 ip=adr_ep_head; ip<adr_ep_tail; ip++){
                    //comm_table_.ep_send_[ip] = epj_sorted_loc_[ id_ep_send_buf_[0][ip] ];
                    comm_table_.ep_send_[ip] = epj_sorted_loc_[ comm_table_.adr_ep_domain_send_[ip] ];
                }
                MPI_Isend( comm_table_.ep_send_.getPointer(adr_ep_head), n_ep, GetDataType<Tepj>(), rank_send, 0, MPI_COMM_WORLD, comm_table_.req_ep_send_.getPointer(n_proc_ep_send));
                n_proc_ep_send++;
            }
            if( comm_table_.n_sp_domain_send_[i] > 0){
                S32 n_sp = comm_table_.n_sp_domain_send_[i];
                S32 adr_sp_head = comm_table_.n_disp_sp_domain_send_[i];
                S32 adr_sp_tail = comm_table_.n_disp_sp_domain_send_[i+1];
                for(S32 ip=adr_sp_head; ip<adr_sp_tail; ip++){
                    //comm_table_.sp_send_[ip] = spj_sorted_loc_[ id_sp_send_buf_[0][ip] ];
                    comm_table_.sp_send_[ip] = spj_sorted_loc_[ comm_table_.adr_sp_domain_send_[ip] ];
                }                
                MPI_Isend( comm_table_.sp_send_.getPointer(adr_sp_head), n_sp, GetDataType<Tspj>(), rank_send, 1, MPI_COMM_WORLD, comm_table_.req_sp_send_.getPointer(n_proc_sp_send));
                n_proc_sp_send++;
            }
        }

    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN
     /*
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"OK 06 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            std::cerr<<"comm_table_.ep_send_.size()= "<<comm_table_.ep_send_.size()
                     <<" comm_table_.sp_send_.size()= "<<comm_table_.sp_send_.size()
                     <<std::endl;
            for(S32 i=0; i<comm_table_.ep_send_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.ep_send_[i].mass= "<<comm_table_.ep_send_[i].mass
                         <<" comm_table_.ep_send_[i].pos= "<<comm_table_.ep_send_[i].pos
                         <<std::endl;
            }
            for(S32 i=0; i<comm_table_.sp_send_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.sp_send_[i].mass= "<<comm_table_.sp_send_[i].mass
                         <<" comm_table_.sp_send_[i].pos= "<<comm_table_.sp_send_[i].pos
                         <<std::endl;
            }
        }
        Comm::barrier();
        //exit(1);
        */
    #endif

        for(S32 i=0; i<comm_table_.rank_ptcl_domain_recv_.size(); i++){
            S32 rank_recv = comm_table_.rank_ptcl_domain_recv_[i];
            if( comm_table_.n_ep_domain_recv_[i] > 0){
                S32 n_ep = comm_table_.n_ep_domain_recv_[i];
                S32 adr_ep_head = comm_table_.n_disp_ep_domain_recv_[i] + adr_ep_recv_offset_2;
                MPI_Irecv( comm_table_.ep_recv_.getPointer(adr_ep_head), n_ep, GetDataType<Tepj>(), rank_recv, 0, MPI_COMM_WORLD, comm_table_.req_ep_recv_.getPointer(n_proc_ep_recv));
                n_proc_ep_recv++;
            }
            if( comm_table_.n_sp_domain_recv_[i] > 0){
                S32 n_sp = comm_table_.n_sp_domain_recv_[i];
                S32 adr_sp_head = comm_table_.n_disp_sp_domain_recv_[i] + adr_sp_recv_offset_3;
                MPI_Irecv( comm_table_.sp_recv_.getPointer(adr_sp_head), n_sp, GetDataType<Tspj>(), rank_recv, 1, MPI_COMM_WORLD, comm_table_.req_sp_recv_.getPointer(n_proc_sp_recv));
                n_proc_sp_recv++;
            }
        }
        MPI_Waitall(n_proc_ep_recv, comm_table_.req_ep_recv_.getPointer(), comm_table_.stat_ep_recv_.getPointer());
        MPI_Waitall(n_proc_ep_send, comm_table_.req_ep_send_.getPointer(), comm_table_.stat_ep_send_.getPointer());
        MPI_Waitall(n_proc_sp_recv, comm_table_.req_sp_recv_.getPointer(), comm_table_.stat_sp_recv_.getPointer());
        MPI_Waitall(n_proc_sp_send, comm_table_.req_sp_send_.getPointer(), comm_table_.stat_sp_send_.getPointer());
        // exchange EPJ and SPJ of normal domain
        ////////////////////////////////////////

    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN
/*
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"OK 07 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            std::cerr<<"comm_table_.ep_recv_.size()= "<<comm_table_.ep_recv_.size()
                     <<" comm_table_.sp_recv_.size()= "<<comm_table_.sp_recv_.size()
                     <<std::endl;
            for(S32 i=0; i<comm_table_.ep_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.ep_recv_[i].mass= "<<comm_table_.ep_recv_[i].mass
                         <<" comm_table_.ep_recv_[i].pos= "<<comm_table_.ep_recv_[i].pos
                         <<std::endl;
            }
            for(S32 i=0; i<comm_table_.sp_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.sp_recv_[i].mass= "<<comm_table_.sp_recv_[i].mass
                         <<" comm_table_.sp_recv_[i].pos= "<<comm_table_.sp_recv_[i].pos
                         <<std::endl;
            }
        }
        Comm::barrier();
        //exit(1);
        */
    #endif

    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
        Comm::barrier(); 
        if(Comm::getRank() == 0){ 
            std::cerr<<"OK 08 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            F64 cm_mass = 0.0;
            F64vec cm_pos = 0.0;
            std::cerr<<"comm_table_.ep_recv_.size()= "<<comm_table_.ep_recv_.size()
                     <<" comm_table_.sp_recv_.size()= "<<comm_table_.sp_recv_.size()
                     <<std::endl;
            for(S32 i=0; i<n_loc_tot_; i++){
                cm_mass += epj_sorted_loc_[i].mass;
                cm_pos += epj_sorted_loc_[i].mass * epj_sorted_loc_[i].pos;
            }
            std::cerr<<"A) cm_mass= "<<cm_mass
                     <<" cm_pos= "<<cm_pos / cm_mass
                     <<std::endl;            
            for(S32 i=0; i<comm_table_.ep_recv_.size(); i++){
                //std::cerr<<"i= "<<i<<" mass= "<<comm_table_.ep_recv_[i].mass<<std::endl;
                cm_mass += comm_table_.ep_recv_[i].mass;
                cm_pos += comm_table_.ep_recv_[i].mass * comm_table_.ep_recv_[i].pos;
            }
            std::cerr<<"B) cm_mass= "<<cm_mass
                     <<" cm_pos= "<<cm_pos / cm_mass
                     <<std::endl;
            for(S32 i=0; i<comm_table_.sp_recv_.size(); i++){

                cm_mass += comm_table_.sp_recv_[i].mass;
                cm_pos += comm_table_.sp_recv_[i].mass * comm_table_.sp_recv_[i].pos;
            }
            cm_pos /= cm_mass;
            std::cerr<<"C) cm_mass= "<<cm_mass
                     <<" cm_mass= "<<cm_pos
                     <<std::endl;
        }
        //exit(1);
    #endif

        time_profile_.exchange_LET_tot += GetWtime() - time_offset;
        
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeSuperDomain2(const DomainInfo & dinfo,
                                           const bool reuse_list){
        const F64 time_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 n_proc_1d = dinfo.getNDomain(0);
        const S32 n_proc_sub = n_proc / n_proc_1d;
        const S32 my_rank = Comm::getRank();
        const S32 my_rank_1d = dinfo.getRank1d(0);
        const S32 my_rank_sub = dinfo.getRankSub(0);
        const F64 PI = 4.0*atan(1.0);
        static ReallocatableArray<F64ort> outer_pos_super_domain;
        static ReallocatableArray<F64ort> outer_pos_sub_domain;
        MPI_Comm comm_sub = dinfo.getCommSub(0);
        MPI_Comm comm_1d = dinfo.getComm1d(0);
        outer_pos_super_domain.resizeNoInitialize(n_proc_1d);
        outer_pos_sub_domain.resizeNoInitialize(n_proc_sub);
        F64 wtime_offset_in;
        #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
        Comm::barrier();
        if (Comm::getRank() == 0) std::cerr << "OK1 @exchangeLocalEssentialTreeSuperDomain2" << std::endl;
        #endif
        //////////////////////
        // initialize
        if(!reuse_list){
            #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
            Comm::barrier();
            if (Comm::getRank() == 0) std::cerr << "(a) OK2 @exchangeLocalEssentialTreeSuperDomain2" << std::endl;
            #endif
            ///////////////////////////
            // exchange outer boundary of super domain
            wtime_offset_in = GetWtime();
            F64ort my_outer_pos_sub_domain = tc_loc_[0].mom_.getVertexOut();
            outer_pos_sub_domain.resizeNoInitialize(n_proc_sub);
            MPI_Allgather(&my_outer_pos_sub_domain, 1, GetDataType<F64ort>(),
                          outer_pos_sub_domain.getPointer(), 1, GetDataType<F64ort>(),
                          comm_sub);
            #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
            Comm::barrier();
            if (Comm::getRank() == 0) std::cerr << "(a) OK3 @exchangeLocalEssentialTreeSuperDomain2" << std::endl;
            #endif
            F64ort my_outer_pos_super_domain;
            for(S32 i=0; i<n_proc_sub; i++) my_outer_pos_super_domain.merge(outer_pos_sub_domain[i]);
            #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
            Comm::barrier();
            if (Comm::getRank() == 0) std::cerr << "(a) OK4 @exchangeLocalEssentialTreeSuperDomain2" << std::endl;
            #endif
            #if 1
            if(my_rank_sub == 0){
                MPI_Allgather(&my_outer_pos_super_domain, 1, GetDataType<F64ort>(),
                              outer_pos_super_domain.getPointer(), 1, GetDataType<F64ort>(),
                              comm_1d);
            }
            MPI_Bcast(outer_pos_super_domain.getPointer(), n_proc_1d, GetDataType<F64ort>(), 0, comm_sub);
            #else
            MPI_Allgather(&my_outer_pos_super_domain, 1, GetDataType<F64ort>(),
                          outer_pos_super_domain.getPointer(), 1, GetDataType<F64ort>(),
                          comm_1d);
            #endif
            #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
            Comm::barrier();
            if (Comm::getRank() == 0) std::cerr << "(a) OK5 @exchangeLocalEssentialTreeSuperDomain2" << std::endl;
            #endif
            wtime_ex_let_sd_2_allgather_bd += GetWtime() - wtime_offset_in ;
            // exchange outer boundary of super domain
            ///////////////////////////
            #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
            Comm::barrier();
            if(my_rank == 0){
                for(S32 i=0; i<n_proc_1d; i++){
                    std::cerr<<"outer_pos_super_domain[i]= "<<outer_pos_super_domain[i]<<std::endl;
                }
            }
            Comm::barrier();
            #endif
            const F64 my_longest_len_super_domain = outer_pos_super_domain[my_rank_1d].getFullLength().getMax();
            const F64 r_crit_sq_send_super_domain = (my_longest_len_super_domain * my_longest_len_super_domain) / (theta_ * theta_);
            const F64 r_crit_sq_tree_walk = (length_ * length_) / (theta_ * theta_);
            const F64 len_peri_x = 2.0*PI;
            comm_table_.clearSize();
            ReallocatableArray<Tspj> sp_list_dummy;
            S32 n_ep_cum = 0;
            S32 n_sp_cum = 0;
            S32 n_ep_cum_2 = 0;
            S32 n_sp_cum_2 = 0;
            wtime_offset_in = GetWtime();
            for(S32 ib=0; ib<n_proc_1d; ib++){
                // send ptcl
                if(GetDistanceMinSqPeriodicX(outer_pos_super_domain[ib], my_outer_pos_super_domain, len_peri_x) > r_crit_sq_send_super_domain){
                    comm_table_.rank_sd_sd_send_.push_back(ib);
                }
                else{
                    if(my_rank_1d == ib){
                        for(S32 jb=0; jb<n_proc_sub; jb++){
                            const S32 rank_target = ib*n_proc_sub + jb;
                            if(my_rank == rank_target) continue;
                            const F64ort pos_target_domain = dinfo.getPosDomain(rank_target);
                            MakeListUsingTreeRecursiveTop<TSM, TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj>
                                (tc_loc_, 0, tp_loc_,
                                 epj_sorted_loc_, comm_table_.adr_ep_domain_send_,
                                 sp_list_dummy,   comm_table_.adr_sp_domain_send_,
                                 pos_target_domain, r_crit_sq_tree_walk, n_leaf_limit_, 0, len_peri_x);
                            comm_table_.rank_ptcl_domain_send_.push_back(rank_target);
                            comm_table_.n_ep_domain_send_.push_back(comm_table_.adr_ep_domain_send_.size() - n_ep_cum);
                            comm_table_.n_sp_domain_send_.push_back(comm_table_.adr_sp_domain_send_.size() - n_sp_cum);
                            n_ep_cum = comm_table_.adr_ep_domain_send_.size();
                            n_sp_cum = comm_table_.adr_sp_domain_send_.size();

                        }
                    }
                    else{
                        const F64ort pos_target_domain = outer_pos_super_domain[ib];
                        MakeListUsingTreeRecursiveTop<TSM, TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj>
                            (tc_loc_, 0, tp_loc_,
                             epj_sorted_loc_, comm_table_.adr_ep_sd_send_,
                             sp_list_dummy,   comm_table_.adr_sp_sd_send_,
                             pos_target_domain, r_crit_sq_tree_walk, n_leaf_limit_, 0, len_peri_x);
                            comm_table_.rank_ptcl_sd_send_.push_back(ib);
                            comm_table_.n_ep_sd_send_.push_back(comm_table_.adr_ep_sd_send_.size() - n_ep_cum_2);
                            comm_table_.n_sp_sd_send_.push_back(comm_table_.adr_sp_sd_send_.size() - n_sp_cum_2);
                            n_ep_cum_2 = comm_table_.adr_ep_sd_send_.size();
                            n_sp_cum_2 = comm_table_.adr_sp_sd_send_.size();
                    }
                }
                // recv ptcl
                const F64 target_longest_len_super_domain = outer_pos_super_domain[ib].getFullLength().getMax();
                const F64 r_crit_sq_recv_super_domain = (target_longest_len_super_domain * target_longest_len_super_domain) / (theta_ * theta_);
                if(GetDistanceMinSqPeriodicX(my_outer_pos_super_domain, outer_pos_super_domain[ib], len_peri_x) > r_crit_sq_recv_super_domain){
                    comm_table_.rank_sd_sd_recv_.push_back(ib);
                }
                else{
                    if(my_rank_1d == ib){
                        for(S32 jb=0; jb<n_proc_sub; jb++){
                            const S32 rank_target = ib*n_proc_sub + jb;
                            if(my_rank == rank_target) continue;
                            comm_table_.rank_ptcl_domain_recv_.push_back(rank_target);
                        }
                    }
                    else{
                        comm_table_.rank_ptcl_sd_recv_.push_back(ib);
                    }
                }
            }
            wtime_ex_let_sd_2_make_list += GetWtime() - wtime_offset_in;
            #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
            Comm::barrier();
            if(my_rank==0){
                for(S32 i=0; i<comm_table_.rank_sd_sd_send_.size(); i++){
                    std::cerr<<"comm_table_.rank_sd_sd_send_[i]= "<<comm_table_.rank_sd_sd_send_[i]<<std::endl;
                }
                for(S32 i=0; i<comm_table_.rank_ptcl_domain_send_.size(); i++){
                    std::cerr<<"comm_table_.rank_ptcl_domain_send_[i]= "<<comm_table_.rank_ptcl_domain_send_[i]
                             <<" comm_table_.n_ep_domain_send_[i]= "<<comm_table_.n_ep_domain_send_[i]
                             <<" comm_table_.n_sp_domain_send_[i]= "<<comm_table_.n_sp_domain_send_[i]
                             <<std::endl;
                }
                for(S32 i=0; i<comm_table_.rank_ptcl_sd_send_.size(); i++){
                    std::cerr<<"comm_table_.rank_ptcl_sd_send_[i]= "<<comm_table_.rank_ptcl_sd_send_[i]
                             <<" comm_table_.n_ep_sd_send_[i]= "<<comm_table_.n_ep_sd_send_[i]
                             <<" comm_table_.n_sp_sd_send_[i]= "<<comm_table_.n_sp_sd_send_[i]
                             <<std::endl;
                }
            }
            Comm::barrier();
            #endif
            
            ///////////////////
            // exchange # of ep and sp (domain -> domain) [in the same super domain]
            wtime_offset_in = GetWtime();
            comm_table_.n_ep_domain_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size());
            comm_table_.n_sp_domain_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size());
            comm_table_.req_ep_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size());
            comm_table_.req_sp_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size());
            comm_table_.stat_ep_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size());
            comm_table_.stat_sp_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size());
            comm_table_.req_ep_send_.resizeNoInitialize(comm_table_.rank_ptcl_domain_send_.size());
            comm_table_.req_sp_send_.resizeNoInitialize(comm_table_.rank_ptcl_domain_send_.size());
            comm_table_.stat_ep_send_.resizeNoInitialize(comm_table_.rank_ptcl_domain_send_.size());
            comm_table_.stat_sp_send_.resizeNoInitialize(comm_table_.rank_ptcl_domain_send_.size());
            for(S32 i=0; i<comm_table_.rank_ptcl_domain_send_.size(); i++){
                MPI_Isend( comm_table_.n_ep_domain_send_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_domain_send_[i], 0, MPI_COMM_WORLD, comm_table_.req_ep_send_.getPointer(i));
                MPI_Isend( comm_table_.n_sp_domain_send_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_domain_send_[i], 1, MPI_COMM_WORLD, comm_table_.req_sp_send_.getPointer(i));
            }
            for(S32 i=0; i<comm_table_.rank_ptcl_domain_recv_.size(); i++){
                MPI_Irecv( comm_table_.n_ep_domain_recv_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_domain_recv_[i], 0, MPI_COMM_WORLD, comm_table_.req_ep_recv_.getPointer(i));
                MPI_Irecv( comm_table_.n_sp_domain_recv_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_domain_recv_[i], 1, MPI_COMM_WORLD, comm_table_.req_sp_recv_.getPointer(i));
            }

            MPI_Waitall(comm_table_.rank_ptcl_domain_recv_.size(), comm_table_.req_ep_recv_.getPointer(), comm_table_.stat_ep_recv_.getPointer());
            MPI_Waitall(comm_table_.rank_ptcl_domain_send_.size(), comm_table_.req_ep_send_.getPointer(), comm_table_.stat_ep_send_.getPointer());
            MPI_Waitall(comm_table_.rank_ptcl_domain_recv_.size(), comm_table_.req_sp_recv_.getPointer(), comm_table_.stat_sp_recv_.getPointer());
            MPI_Waitall(comm_table_.rank_ptcl_domain_send_.size(), comm_table_.req_sp_send_.getPointer(), comm_table_.stat_sp_send_.getPointer());

            comm_table_.n_disp_ep_domain_send_.resizeNoInitialize(comm_table_.rank_ptcl_domain_send_.size()+1);
            comm_table_.n_disp_sp_domain_send_.resizeNoInitialize(comm_table_.rank_ptcl_domain_send_.size()+1);
            comm_table_.n_disp_ep_domain_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size()+1);
            comm_table_.n_disp_sp_domain_recv_.resizeNoInitialize(comm_table_.rank_ptcl_domain_recv_.size()+1);
            comm_table_.n_disp_ep_domain_send_[0] = comm_table_.n_disp_sp_domain_send_[0] = comm_table_.n_disp_ep_domain_recv_[0] = comm_table_.n_disp_sp_domain_recv_[0] = 0;
            for(S32 i=0; i<comm_table_.rank_ptcl_domain_send_.size(); i++){
                comm_table_.n_disp_ep_domain_send_[i+1] = comm_table_.n_disp_ep_domain_send_[i] + comm_table_.n_ep_domain_send_[i];
                comm_table_.n_disp_sp_domain_send_[i+1] = comm_table_.n_disp_sp_domain_send_[i] + comm_table_.n_sp_domain_send_[i];
            }
            for(S32 i=0; i<comm_table_.rank_ptcl_domain_recv_.size(); i++){
                comm_table_.n_disp_ep_domain_recv_[i+1] = comm_table_.n_disp_ep_domain_recv_[i] + comm_table_.n_ep_domain_recv_[i];
                comm_table_.n_disp_sp_domain_recv_[i+1] = comm_table_.n_disp_sp_domain_recv_[i] + comm_table_.n_sp_domain_recv_[i];
            }

            n_ep_send_tot_loc += comm_table_.n_disp_ep_domain_send_[comm_table_.rank_ptcl_domain_send_.size()];
            n_sp_send_tot_loc += comm_table_.n_disp_sp_domain_send_[comm_table_.rank_ptcl_domain_send_.size()];
            
            wtime_ex_let_sd_2_ex_n_d_d += GetWtime() - wtime_offset_in;
            // exchange # of ep and sp (domain -> domain)
            ///////////////////

            ///////////////////
            // exchange # of ep and sp (domain -> super domain)
            wtime_offset_in = GetWtime();
            comm_table_.n_ep_sd_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size());
            comm_table_.n_sp_sd_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size());
            comm_table_.req_ep_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size());
            comm_table_.req_sp_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size());
            comm_table_.stat_ep_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size());
            comm_table_.stat_sp_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size());
            comm_table_.req_ep_send_.resizeNoInitialize(comm_table_.rank_ptcl_sd_send_.size());
            comm_table_.req_sp_send_.resizeNoInitialize(comm_table_.rank_ptcl_sd_send_.size());
            comm_table_.stat_ep_send_.resizeNoInitialize(comm_table_.rank_ptcl_sd_send_.size());
            comm_table_.stat_sp_send_.resizeNoInitialize(comm_table_.rank_ptcl_sd_send_.size());            
            for(S32 i=0; i<comm_table_.rank_ptcl_sd_send_.size(); i++){
                MPI_Isend( comm_table_.n_ep_sd_send_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_sd_send_[i], 0, comm_1d, comm_table_.req_ep_send_.getPointer(i));
                MPI_Isend( comm_table_.n_sp_sd_send_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_sd_send_[i], 1, comm_1d, comm_table_.req_sp_send_.getPointer(i));
            }
            for(S32 i=0; i<comm_table_.rank_ptcl_sd_recv_.size(); i++){
                MPI_Irecv( comm_table_.n_ep_sd_recv_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_sd_recv_[i], 0, comm_1d, comm_table_.req_ep_recv_.getPointer(i));
                MPI_Irecv( comm_table_.n_sp_sd_recv_.getPointer(i), 1, GetDataType<S32>(), comm_table_.rank_ptcl_sd_recv_[i], 1, comm_1d, comm_table_.req_sp_recv_.getPointer(i));
            }            
            MPI_Waitall(comm_table_.rank_ptcl_sd_recv_.size(), comm_table_.req_ep_recv_.getPointer(), comm_table_.stat_ep_recv_.getPointer());
            MPI_Waitall(comm_table_.rank_ptcl_sd_send_.size(), comm_table_.req_ep_send_.getPointer(), comm_table_.stat_ep_send_.getPointer());
            MPI_Waitall(comm_table_.rank_ptcl_sd_recv_.size(), comm_table_.req_sp_recv_.getPointer(), comm_table_.stat_sp_recv_.getPointer());
            MPI_Waitall(comm_table_.rank_ptcl_sd_send_.size(), comm_table_.req_sp_send_.getPointer(), comm_table_.stat_sp_send_.getPointer());

            comm_table_.n_disp_ep_sd_send_.resizeNoInitialize(comm_table_.rank_ptcl_sd_send_.size()+1);
            comm_table_.n_disp_sp_sd_send_.resizeNoInitialize(comm_table_.rank_ptcl_sd_send_.size()+1);
            comm_table_.n_disp_ep_sd_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size()+1);
            comm_table_.n_disp_sp_sd_recv_.resizeNoInitialize(comm_table_.rank_ptcl_sd_recv_.size()+1);
            comm_table_.n_disp_ep_sd_send_[0] = comm_table_.n_disp_sp_sd_send_[0] = comm_table_.n_disp_ep_sd_recv_[0] = comm_table_.n_disp_sp_sd_recv_[0] = 0;
            for(S32 i=0; i<comm_table_.rank_ptcl_sd_send_.size(); i++){
                comm_table_.n_disp_ep_sd_send_[i+1] = comm_table_.n_disp_ep_sd_send_[i] + comm_table_.n_ep_sd_send_[i];
                comm_table_.n_disp_sp_sd_send_[i+1] = comm_table_.n_disp_sp_sd_send_[i] + comm_table_.n_sp_sd_send_[i];
            }
            for(S32 i=0; i<comm_table_.rank_ptcl_sd_recv_.size(); i++){
                comm_table_.n_disp_ep_sd_recv_[i+1] = comm_table_.n_disp_ep_sd_recv_[i] + comm_table_.n_ep_sd_recv_[i];
                comm_table_.n_disp_sp_sd_recv_[i+1] = comm_table_.n_disp_sp_sd_recv_[i] + comm_table_.n_sp_sd_recv_[i];
            }

            n_ep_send_tot_loc += comm_table_.n_disp_ep_sd_send_[comm_table_.rank_ptcl_sd_send_.size()];
            n_sp_send_tot_loc += comm_table_.n_disp_sp_sd_send_[comm_table_.rank_ptcl_sd_send_.size()];
            
            wtime_ex_let_sd_2_ex_n_d_sd_1d_ring += GetWtime() - wtime_offset_in;
            // exchange # of ep and sp (domain -> super domain)
            ///////////////////

            ///////////////////
            // exchange # of ep and sp (domain -> super domain) in super domain
            wtime_offset_in = GetWtime();
            comm_table_.n_ep_domain_recv_2_.resizeNoInitialize(n_proc_sub);
            comm_table_.n_sp_domain_recv_2_.resizeNoInitialize(n_proc_sub);
            comm_table_.n_disp_ep_domain_recv_2_.resizeNoInitialize(n_proc_sub+1);
            comm_table_.n_disp_sp_domain_recv_2_.resizeNoInitialize(n_proc_sub+1);
            comm_table_.n_ep_domain_send_2_ = comm_table_.n_disp_ep_sd_recv_[ comm_table_.rank_ptcl_sd_recv_.size() ];
            MPI_Allgather( &(comm_table_.n_ep_domain_send_2_), 1, GetDataType<S32>(),
                           comm_table_.n_ep_domain_recv_2_.getPointer(), 1, GetDataType<S32>(),
                           comm_sub);
            comm_table_.n_sp_domain_send_2_ = comm_table_.n_disp_sp_sd_recv_[ comm_table_.rank_ptcl_sd_recv_.size() ];
            MPI_Allgather( &(comm_table_.n_sp_domain_send_2_), 1, GetDataType<S32>(),
                           comm_table_.n_sp_domain_recv_2_.getPointer(), 1, GetDataType<S32>(),
                           comm_sub);
            comm_table_.n_disp_ep_domain_recv_2_[0] = comm_table_.n_disp_sp_domain_recv_2_[0] = 0;
            comm_table_.n_ep_domain_recv_2_[my_rank_sub] = comm_table_.n_sp_domain_recv_2_[my_rank_sub] = 0;
            for(S32 i=0; i<n_proc_sub; i++){
                comm_table_.n_disp_ep_domain_recv_2_[i+1] = comm_table_.n_disp_ep_domain_recv_2_[i] + comm_table_.n_ep_domain_recv_2_[i];
                comm_table_.n_disp_sp_domain_recv_2_[i+1] = comm_table_.n_disp_sp_domain_recv_2_[i] + comm_table_.n_sp_domain_recv_2_[i];
            }

            n_ep_send_tot_loc += comm_table_.n_disp_ep_domain_recv_2_[n_proc_sub];
            n_sp_send_tot_loc += comm_table_.n_disp_sp_domain_recv_2_[n_proc_sub];
            
            wtime_ex_let_sd_2_ex_n_d_sd_sub += GetWtime() - wtime_offset_in;
            // exchange # of ep and sp (domain -> super domain) in super domain
            ///////////////////
            #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
            if(my_rank == 0){
                std::cerr<<std::endl;
                comm_table_.dumpNumber();
                std::cerr<<std::endl;
            }
            #endif
        } // end of reuse
        // initialize
        //////////////////////
        #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
        Comm::barrier();
        if (Comm::getRank() == 0) std::cerr << "OK2 @exchangeLocalEssentialTreeSuperDomain2" << std::endl;
        #endif

        ///////////////////////////
        // allgather SPJ of super domain
        wtime_offset_in = GetWtime();
        static ReallocatableArray<Tspj> spj_sub_domain;
        spj_sub_domain.resizeNoInitialize(n_proc_sub);
        #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
        Comm::barrier();
        if (Comm::getRank() == 0) std::cerr << "OK3 @exchangeLocalEssentialTreeSuperDomain2" << std::endl;
        #endif
        MPI_Allgather(&spj_sorted_loc_[0], 1, GetDataType<Tspj>(),
                      spj_sub_domain.getPointer(), 1, GetDataType<Tspj>(),
                      comm_sub);
        #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
        Comm::barrier();
        if (Comm::getRank() == 0) std::cerr << "OK4 @exchangeLocalEssentialTreeSuperDomain2" << std::endl;
        if(my_rank == 0){
            std::cerr<<std::endl;
            for(S32 i=0; i<n_proc_sub; i++){
                std::cerr<<"i="<<i<<std::endl;
                spj_sub_domain[i].dump();
            }
            std::cerr<<std::endl;
        }
        #endif
        Tspj my_spj_super_domain;
        my_spj_super_domain.clear();
        for(S32 i=0; i<n_proc_sub; i++){
            my_spj_super_domain.mass += spj_sub_domain[i].mass;
            my_spj_super_domain.pos += spj_sub_domain[i].mass * spj_sub_domain[i].pos;
            my_spj_super_domain.pos_phi += spj_sub_domain[i].mass * spj_sub_domain[i].pos_phi;
            my_spj_super_domain.pos_r += spj_sub_domain[i].mass * spj_sub_domain[i].pos_r;
        }
        my_spj_super_domain.pos     /= my_spj_super_domain.mass;
        my_spj_super_domain.pos_phi /= my_spj_super_domain.mass;
        my_spj_super_domain.pos_r   /= my_spj_super_domain.mass;
        #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
        Comm::barrier();
        if (Comm::getRank() == 0) std::cerr << "OK5 @exchangeLocalEssentialTreeSuperDomain2" << std::endl;
        if(my_rank == 0){
            std::cerr<<"super domain"<<std::endl;
            my_spj_super_domain.dump();
        }
        #endif
        static ReallocatableArray<Tspj> spj_super_domain;
        spj_super_domain.resizeNoInitialize(n_proc_1d);
        #if 1
        if(my_rank_sub == 0){
            MPI_Allgather(&my_spj_super_domain, 1, GetDataType<Tspj>(),
                          spj_super_domain.getPointer(), 1, GetDataType<Tspj>(),
                          comm_1d);
        }
        MPI_Bcast(spj_super_domain.getPointer(), n_proc_1d, GetDataType<Tspj>(), 0, comm_sub);
        #else
        MPI_Allgather(&my_spj_super_domain, 1, GetDataType<Tspj>(),
                      spj_super_domain.getPointer(), 1, GetDataType<Tspj>(),
                      comm_1d);
        #endif
        #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
        Comm::barrier();
        if (Comm::getRank() == 0) std::cerr << "OK6 @exchangeLocalEssentialTreeSuperDomain2" << std::endl;
        if(my_rank == 0){
            std::cerr<<"spj super domain"<<std::endl;
            for(S32 i=0; i<n_proc_1d; i++){
                std::cerr<<"i= "<<i<<" n_proc_1d= "<<n_proc_1d<<std::endl;
                spj_super_domain[i].dump();
            }
        }
        #endif
        comm_table_.sp_recv_.reserve(comm_table_.rank_sd_sd_recv_.size());
        comm_table_.sp_recv_.clearSize();
        for(S32 i=0; i<comm_table_.rank_sd_sd_recv_.size(); i++){
            S32 rank = comm_table_.rank_sd_sd_recv_[i];
            Tspj spj = spj_super_domain[rank];
            comm_table_.sp_recv_.pushBackNoCheck(spj);
        }
        #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
        Comm::barrier();
        if (Comm::getRank() == 0) std::cerr << "OK7 @exchangeLocalEssentialTreeSuperDomain2" << std::endl;
        if(my_rank == 0){
            for(S32 i=0; i<comm_table_.rank_sd_sd_recv_.size(); i++){
                std::cerr<<"comm_table_.rank_sd_sd_recv_[i]= "<<comm_table_.rank_sd_sd_recv_[i]<<std::endl;
                std::cerr<<"comm_table_.sp_recv_[i].mass= "<<comm_table_.sp_recv_[i].mass<<std::endl;
                std::cerr<<"comm_table_.sp_recv_[i].pos_phi= "<<comm_table_.sp_recv_[i].pos_phi<<std::endl;
                std::cerr<<"comm_table_.sp_recv_[i].pos_r= "<<comm_table_.sp_recv_[i].pos_r<<std::endl;
                std::cerr<<std::endl;
            }
        }
        //exit(1);
        #endif
        wtime_ex_let_sd_2_allgather_sd += GetWtime() - wtime_offset_in;
        S32 adr_sp_recv_offset = comm_table_.sp_recv_.size();
        // allgather SPJ of super domain
        ///////////////////////////

        ////////////////////////////////////
        // exchange EPJ and SPJ of super domain among super domains
        wtime_offset_in = GetWtime();
        S32 n_proc_ep_send = 0;
        S32 n_proc_sp_send = 0;
        S32 n_proc_ep_recv = 0;
        S32 n_proc_sp_recv = 0;
        comm_table_.ep_send_.resizeNoInitialize( comm_table_.n_disp_ep_sd_send_[comm_table_.rank_ptcl_sd_send_.size()] );
        comm_table_.sp_send_.resizeNoInitialize( comm_table_.n_disp_sp_sd_send_[comm_table_.rank_ptcl_sd_send_.size()] );
        comm_table_.ep_recv_.resizeNoInitialize( comm_table_.n_disp_ep_sd_recv_[comm_table_.rank_ptcl_sd_recv_.size()] );
        comm_table_.sp_recv_.resizeNoInitialize( comm_table_.n_disp_sp_sd_recv_[comm_table_.rank_ptcl_sd_recv_.size()] + adr_sp_recv_offset);
        for(S32 i=0; i<comm_table_.rank_ptcl_sd_send_.size(); i++){
            S32 rank_send = comm_table_.rank_ptcl_sd_send_[i];
            if( comm_table_.n_ep_sd_send_[i] > 0){
                S32 n_ep = comm_table_.n_ep_sd_send_[i];
                S32 adr_ep_head = comm_table_.n_disp_ep_sd_send_[i];
                S32 adr_ep_tail = comm_table_.n_disp_ep_sd_send_[i+1];
                for(S32 ip=adr_ep_head; ip<adr_ep_tail; ip++){
                    comm_table_.ep_send_[ip] = epj_sorted_loc_[ comm_table_.adr_ep_sd_send_[ip] ];
                }
                MPI_Isend( comm_table_.ep_send_.getPointer(adr_ep_head), n_ep, GetDataType<Tepj>(), rank_send, 0, comm_1d, comm_table_.req_ep_send_.getPointer(n_proc_ep_send));
                n_proc_ep_send++;
            }
            if( comm_table_.n_sp_sd_send_[i] > 0){
                S32 n_sp = comm_table_.n_sp_sd_send_[i];
                S32 adr_sp_head = comm_table_.n_disp_sp_sd_send_[i];
                S32 adr_sp_tail = comm_table_.n_disp_sp_sd_send_[i+1];
                for(S32 ip=adr_sp_head; ip<adr_sp_tail; ip++){
                    comm_table_.sp_send_[ip] = spj_sorted_loc_[ comm_table_.adr_sp_sd_send_[ip] ];
                }
                MPI_Isend( comm_table_.sp_send_.getPointer(adr_sp_head), n_sp, GetDataType<Tspj>(), rank_send, 0, comm_1d, comm_table_.req_sp_send_.getPointer(n_proc_sp_send));
                n_proc_sp_send++;
            }            
        }
        for(S32 i=0; i<comm_table_.rank_ptcl_sd_recv_.size(); i++){
            S32 rank_recv = comm_table_.rank_ptcl_sd_recv_[i];
            if( comm_table_.n_ep_sd_recv_[i] > 0){
                S32 n_ep = comm_table_.n_ep_sd_recv_[i];
                S32 adr_ep_head = comm_table_.n_disp_ep_sd_recv_[i];
                MPI_Irecv( comm_table_.ep_recv_.getPointer(adr_ep_head), n_ep, GetDataType<Tepj>(), rank_recv, 0, comm_1d, comm_table_.req_ep_recv_.getPointer(n_proc_ep_recv));
                n_proc_ep_recv++;
            }
            if( comm_table_.n_sp_sd_recv_[i] > 0){
                S32 n_sp = comm_table_.n_sp_sd_recv_[i];
                S32 adr_sp_head = comm_table_.n_disp_sp_sd_recv_[i] + adr_sp_recv_offset;
                MPI_Irecv( comm_table_.sp_recv_.getPointer(adr_sp_head), n_sp, GetDataType<Tspj>(), rank_recv, 0, comm_1d, comm_table_.req_sp_recv_.getPointer(n_proc_sp_recv));
                n_proc_sp_recv++;
            }
        }
        MPI_Waitall(n_proc_ep_recv, comm_table_.req_ep_recv_.getPointer(), comm_table_.stat_ep_recv_.getPointer());
        MPI_Waitall(n_proc_ep_send, comm_table_.req_ep_send_.getPointer(), comm_table_.stat_ep_send_.getPointer());
        MPI_Waitall(n_proc_sp_recv, comm_table_.req_sp_recv_.getPointer(), comm_table_.stat_sp_recv_.getPointer());
        MPI_Waitall(n_proc_sp_send, comm_table_.req_sp_send_.getPointer(), comm_table_.stat_sp_send_.getPointer());

        S32 adr_sp_recv_offset_2 = comm_table_.n_disp_sp_sd_recv_[comm_table_.rank_ptcl_sd_recv_.size()] + adr_sp_recv_offset;
        S32 adr_ep_recv_offset   = comm_table_.n_disp_ep_sd_recv_[comm_table_.rank_ptcl_sd_recv_.size()];

    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
        /*
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"OK 01 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            std::cerr<<"comm_table_.ep_recv_.size()= "<<comm_table_.ep_recv_.size()
                     <<" comm_table_.sp_recv_.size()= "<<comm_table_.sp_recv_.size()
                     <<std::endl;
            for(S32 i=0; i<adr_ep_recv_offset; i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.ep_recv_[i].mass= "<<comm_table_.ep_recv_[i].mass
                         <<" comm_table_.ep_recv_[i].pos= "<<comm_table_.ep_recv_[i].pos
                         <<std::endl;
            }
            for(S32 i=0; i<adr_sp_recv_offset_2; i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.sp_recv_[i].mass= "<<comm_table_.sp_recv_[i].mass
                         <<" comm_table_.sp_recv_[i].pos= "<<comm_table_.sp_recv_[i].pos
                         <<std::endl;
            }
        }
        Comm::barrier();
        //exit(1);
        */
    #endif

        // in subdomain
        comm_table_.ep_recv_.resizeNoInitialize(adr_ep_recv_offset   + comm_table_.n_disp_ep_domain_recv_2_[n_proc_sub]);
        comm_table_.sp_recv_.resizeNoInitialize(adr_sp_recv_offset_2 + comm_table_.n_disp_sp_domain_recv_2_[n_proc_sub]);

    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
     /*
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"OK 02 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            std::cerr<<"comm_table_.ep_recv_.size()= "<<comm_table_.ep_recv_.size()
                     <<" comm_table_.sp_recv_.size()= "<<comm_table_.sp_recv_.size()
                     <<std::endl;
            std::cerr<<"comm_table_.n_ep_domain_send_2_= "<<comm_table_.n_ep_domain_send_2_
                     <<" comm_table_.n_sp_domain_send_2_= "<<comm_table_.n_sp_domain_send_2_
                     <<std::endl;
            std::cerr<<"adr_sp_recv_offset= "<<adr_sp_recv_offset<<std::endl;
            for(S32 i=0; i<comm_table_.ep_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.ep_recv_[i].mass= "<<comm_table_.ep_recv_[i].mass
                         <<" comm_table_.ep_recv_[i].pos= "<<comm_table_.ep_recv_[i].pos
                         <<std::endl;
            }
            for(S32 i=0; i<comm_table_.sp_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.sp_recv_[i].mass= "<<comm_table_.sp_recv_[i].mass
                         <<" comm_table_.sp_recv_[i].pos= "<<comm_table_.sp_recv_[i].pos
                         <<std::endl;
            }
        }
        Comm::barrier();
        //exit(1);
        */
    #endif

        n_proc_ep_send = n_proc_sp_send = n_proc_ep_recv = n_proc_sp_recv = 0;
        for(S32 i=0; i<n_proc_sub; i++){
            if(i==my_rank_sub) continue;
            S32 n_ep = comm_table_.n_ep_domain_send_2_;
            S32 n_sp = comm_table_.n_sp_domain_send_2_;
            MPI_Isend( comm_table_.ep_recv_.getPointer(), n_ep, GetDataType<Tepj>(), i, 0, comm_sub, comm_table_.req_ep_send_.getPointer(n_proc_ep_send));
            MPI_Isend( comm_table_.sp_recv_.getPointer(adr_sp_recv_offset), n_sp, GetDataType<Tspj>(), i, 1, comm_sub, comm_table_.req_sp_send_.getPointer(n_proc_sp_send));
            n_proc_ep_send++;
            n_proc_sp_send++;
        }
        for(S32 i=0; i<n_proc_sub; i++){
            if(i==my_rank_sub) continue;
            S32 n_ep = comm_table_.n_ep_domain_recv_2_[i];
            S32 n_sp = comm_table_.n_sp_domain_recv_2_[i];
            MPI_Irecv( comm_table_.ep_recv_.getPointer(comm_table_.n_disp_ep_domain_recv_2_[i]+adr_ep_recv_offset), n_ep, GetDataType<Tepj>(), i, 0, comm_sub, comm_table_.req_ep_recv_.getPointer(n_proc_ep_recv));
            MPI_Irecv( comm_table_.sp_recv_.getPointer(comm_table_.n_disp_sp_domain_recv_2_[i]+adr_sp_recv_offset_2), n_sp, GetDataType<Tspj>(), i, 1, comm_sub, comm_table_.req_sp_recv_.getPointer(n_proc_sp_recv));
            n_proc_ep_recv++;
            n_proc_sp_recv++;
        }
        MPI_Waitall(n_proc_ep_recv, comm_table_.req_ep_recv_.getPointer(), comm_table_.stat_ep_recv_.getPointer());
        MPI_Waitall(n_proc_ep_send, comm_table_.req_ep_send_.getPointer(), comm_table_.stat_ep_send_.getPointer());
        MPI_Waitall(n_proc_sp_recv, comm_table_.req_sp_recv_.getPointer(), comm_table_.stat_sp_recv_.getPointer());
        MPI_Waitall(n_proc_sp_send, comm_table_.req_sp_send_.getPointer(), comm_table_.stat_sp_send_.getPointer());
        
    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
        /*
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"OK 03 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            std::cerr<<"comm_table_.ep_recv_.size()= "<<comm_table_.ep_recv_.size()
                     <<" comm_table_.sp_recv_.size()= "<<comm_table_.sp_recv_.size()
                     <<std::endl;
            for(S32 i=0; i<comm_table_.ep_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.ep_recv_[i].mass= "<<comm_table_.ep_recv_[i].mass
                         <<" comm_table_.ep_recv_[i].pos= "<<comm_table_.ep_recv_[i].pos
                         <<std::endl;
            }
            for(S32 i=0; i<comm_table_.sp_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.sp_recv_[i].mass= "<<comm_table_.sp_recv_[i].mass
                         <<" comm_table_.sp_recv_[i].pos= "<<comm_table_.sp_recv_[i].pos
                         <<std::endl;
            }
        }
        Comm::barrier();
        //exit(1);
        */
    #endif

        S32 adr_ep_recv_offset_2 = comm_table_.n_disp_ep_domain_recv_2_[n_proc_sub] + adr_ep_recv_offset;
        S32 adr_sp_recv_offset_3 = comm_table_.n_disp_sp_domain_recv_2_[n_proc_sub] + adr_sp_recv_offset_2;
        
    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
     /*
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"OK 04 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            std::cerr<<"comm_table_.ep_recv_.size()= "<<comm_table_.ep_recv_.size()
                     <<" comm_table_.sp_recv_.size()= "<<comm_table_.sp_recv_.size()
                     <<std::endl;
            std::cerr<<"adr_ep_recv_offset= "<<adr_ep_recv_offset
                     <<" adr_ep_recv_offset_2= "<<adr_ep_recv_offset_2
                     <<std::endl;
            std::cerr<<"adr_sp_recv_offset= "<<adr_sp_recv_offset
                     <<" adr_sp_recv_offset_2= "<<adr_sp_recv_offset_2
                     <<" adr_sp_recv_offset_3= "<<adr_sp_recv_offset_3
                     <<std::endl;
            for(S32 i=0; i<comm_table_.ep_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.ep_recv_[i].mass= "<<comm_table_.ep_recv_[i].mass
                         <<" comm_table_.ep_recv_[i].pos= "<<comm_table_.ep_recv_[i].pos
                         <<std::endl;
            }
            for(S32 i=0; i<comm_table_.sp_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.sp_recv_[i].mass= "<<comm_table_.sp_recv_[i].mass
                         <<" comm_table_.sp_recv_[i].pos= "<<comm_table_.sp_recv_[i].pos
                         <<std::endl;
            }
        }
        Comm::barrier();
        //exit(1);
        */
    #endif
        wtime_ex_let_sd_2_ep_sp_among_sd += GetWtime() - wtime_offset_in;
        // exchange EPJ and SPJ of super domain
        ////////////////////////////////////

        ////////////////////////////////////////
        // exchange EPJ and SPJ in super domain
        wtime_offset_in = GetWtime();
        n_proc_ep_send = n_proc_sp_send = n_proc_ep_recv = n_proc_sp_recv = 0;
        comm_table_.ep_send_.resizeNoInitialize(comm_table_.n_disp_ep_domain_send_[comm_table_.rank_ptcl_domain_send_.size()]);
        comm_table_.sp_send_.resizeNoInitialize(comm_table_.n_disp_sp_domain_send_[comm_table_.rank_ptcl_domain_send_.size()]);
        comm_table_.ep_recv_.resizeNoInitialize(comm_table_.n_disp_ep_domain_recv_[comm_table_.rank_ptcl_domain_recv_.size()] + adr_ep_recv_offset_2);
        comm_table_.sp_recv_.resizeNoInitialize(comm_table_.n_disp_sp_domain_recv_[comm_table_.rank_ptcl_domain_recv_.size()] + adr_sp_recv_offset_3);

    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
     /*
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"OK 05 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            std::cerr<<"comm_table_.ep_send_.size()= "<<comm_table_.ep_send_.size()
                     <<" comm_table_.sp_send_.size()= "<<comm_table_.sp_send_.size()
                     <<std::endl;
            std::cerr<<"comm_table_.ep_recv_.size()= "<<comm_table_.ep_recv_.size()
                     <<" comm_table_.sp_recv_.size()= "<<comm_table_.sp_recv_.size()
                     <<std::endl;

            std::cerr<<"comm_table_.rank_ptcl_domain_send_.size()= "<<comm_table_.rank_ptcl_domain_send_.size()
                     <<" comm_table_.rank_ptcl_domain_recv_.size()= "<<comm_table_.rank_ptcl_domain_recv_.size()
                     <<std::endl;
            for(S32 i=0; i<comm_table_.rank_ptcl_domain_send_.size(); i++){
                std::cerr<<"comm_table_.rank_ptcl_domain_send_[i]= "<<comm_table_.rank_ptcl_domain_send_[i]
                         <<" comm_table_.n_ep_domain_send_[i]= "<<comm_table_.n_ep_domain_send_[i]
                         <<" comm_table_.n_sp_domain_send_[i]= "<<comm_table_.n_sp_domain_send_[i]
                         <<std::endl;
            }
            for(S32 i=0; i<comm_table_.rank_ptcl_domain_recv_.size(); i++){
                std::cerr<<"comm_table_.rank_ptcl_domain_recv_[i]= "<<comm_table_.rank_ptcl_domain_recv_[i]
                         <<" comm_table_.n_ep_domain_recv_[i]= "<<comm_table_.n_ep_domain_recv_[i]
                         <<" comm_table_.n_sp_domain_recv_[i]= "<<comm_table_.n_sp_domain_recv_[i]
                         <<std::endl;
            }
        }
        Comm::barrier();
        //exit(1);
        */
    #endif

        n_proc_ep_send = n_proc_sp_send = n_proc_ep_recv = n_proc_sp_recv = 0;
        for(S32 i=0; i<comm_table_.rank_ptcl_domain_send_.size(); i++){
            S32 rank_send = comm_table_.rank_ptcl_domain_send_[i];
            if( comm_table_.n_ep_domain_send_[i] > 0){
                S32 n_ep = comm_table_.n_ep_domain_send_[i];
                S32 adr_ep_head = comm_table_.n_disp_ep_domain_send_[i];
                S32 adr_ep_tail = comm_table_.n_disp_ep_domain_send_[i+1];
                for(S32 ip=adr_ep_head; ip<adr_ep_tail; ip++){
                    //comm_table_.ep_send_[ip] = epj_sorted_loc_[ id_ep_send_buf_[0][ip] ];
                    comm_table_.ep_send_[ip] = epj_sorted_loc_[ comm_table_.adr_ep_domain_send_[ip] ];
                }
                MPI_Isend( comm_table_.ep_send_.getPointer(adr_ep_head), n_ep, GetDataType<Tepj>(), rank_send, 0, MPI_COMM_WORLD, comm_table_.req_ep_send_.getPointer(n_proc_ep_send));
                n_proc_ep_send++;
            }
            if( comm_table_.n_sp_domain_send_[i] > 0){
                S32 n_sp = comm_table_.n_sp_domain_send_[i];
                S32 adr_sp_head = comm_table_.n_disp_sp_domain_send_[i];
                S32 adr_sp_tail = comm_table_.n_disp_sp_domain_send_[i+1];
                for(S32 ip=adr_sp_head; ip<adr_sp_tail; ip++){
                    //comm_table_.sp_send_[ip] = spj_sorted_loc_[ id_sp_send_buf_[0][ip] ];
                    comm_table_.sp_send_[ip] = spj_sorted_loc_[ comm_table_.adr_sp_domain_send_[ip] ];
                }                
                MPI_Isend( comm_table_.sp_send_.getPointer(adr_sp_head), n_sp, GetDataType<Tspj>(), rank_send, 1, MPI_COMM_WORLD, comm_table_.req_sp_send_.getPointer(n_proc_sp_send));
                n_proc_sp_send++;
            }
        }

    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
     /*
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"OK 06 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            std::cerr<<"comm_table_.ep_send_.size()= "<<comm_table_.ep_send_.size()
                     <<" comm_table_.sp_send_.size()= "<<comm_table_.sp_send_.size()
                     <<std::endl;
            for(S32 i=0; i<comm_table_.ep_send_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.ep_send_[i].mass= "<<comm_table_.ep_send_[i].mass
                         <<" comm_table_.ep_send_[i].pos= "<<comm_table_.ep_send_[i].pos
                         <<std::endl;
            }
            for(S32 i=0; i<comm_table_.sp_send_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.sp_send_[i].mass= "<<comm_table_.sp_send_[i].mass
                         <<" comm_table_.sp_send_[i].pos= "<<comm_table_.sp_send_[i].pos
                         <<std::endl;
            }
        }
        Comm::barrier();
        //exit(1);
        */
    #endif

        for(S32 i=0; i<comm_table_.rank_ptcl_domain_recv_.size(); i++){
            S32 rank_recv = comm_table_.rank_ptcl_domain_recv_[i];
            if( comm_table_.n_ep_domain_recv_[i] > 0){
                S32 n_ep = comm_table_.n_ep_domain_recv_[i];
                S32 adr_ep_head = comm_table_.n_disp_ep_domain_recv_[i] + adr_ep_recv_offset_2;
                MPI_Irecv( comm_table_.ep_recv_.getPointer(adr_ep_head), n_ep, GetDataType<Tepj>(), rank_recv, 0, MPI_COMM_WORLD, comm_table_.req_ep_recv_.getPointer(n_proc_ep_recv));
                n_proc_ep_recv++;
            }
            if( comm_table_.n_sp_domain_recv_[i] > 0){
                S32 n_sp = comm_table_.n_sp_domain_recv_[i];
                S32 adr_sp_head = comm_table_.n_disp_sp_domain_recv_[i] + adr_sp_recv_offset_3;
                MPI_Irecv( comm_table_.sp_recv_.getPointer(adr_sp_head), n_sp, GetDataType<Tspj>(), rank_recv, 1, MPI_COMM_WORLD, comm_table_.req_sp_recv_.getPointer(n_proc_sp_recv));
                n_proc_sp_recv++;
            }
        }
        MPI_Waitall(n_proc_ep_recv, comm_table_.req_ep_recv_.getPointer(), comm_table_.stat_ep_recv_.getPointer());
        MPI_Waitall(n_proc_ep_send, comm_table_.req_ep_send_.getPointer(), comm_table_.stat_ep_send_.getPointer());
        MPI_Waitall(n_proc_sp_recv, comm_table_.req_sp_recv_.getPointer(), comm_table_.stat_sp_recv_.getPointer());
        MPI_Waitall(n_proc_sp_send, comm_table_.req_sp_send_.getPointer(), comm_table_.stat_sp_send_.getPointer());
        wtime_ex_let_sd_2_ep_sp_in_sd += GetWtime() - wtime_offset_in;
        // exchange EPJ and SPJ of normal domain
        ////////////////////////////////////////

    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
/*
        Comm::barrier();
        if(my_rank == 0){
            std::cerr<<"OK 07 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            std::cerr<<"comm_table_.ep_recv_.size()= "<<comm_table_.ep_recv_.size()
                     <<" comm_table_.sp_recv_.size()= "<<comm_table_.sp_recv_.size()
                     <<std::endl;
            for(S32 i=0; i<comm_table_.ep_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.ep_recv_[i].mass= "<<comm_table_.ep_recv_[i].mass
                         <<" comm_table_.ep_recv_[i].pos= "<<comm_table_.ep_recv_[i].pos
                         <<std::endl;
            }
            for(S32 i=0; i<comm_table_.sp_recv_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" comm_table_.sp_recv_[i].mass= "<<comm_table_.sp_recv_[i].mass
                         <<" comm_table_.sp_recv_[i].pos= "<<comm_table_.sp_recv_[i].pos
                         <<std::endl;
            }
        }
        Comm::barrier();
        //exit(1);
        */
    #endif

    #ifdef DEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
        Comm::barrier(); 
        if(Comm::getRank() == 0){ 
            std::cerr<<"OK 08 @exchangeLocalEssentialTreeSuperDomain"<<std::endl;
            F64 cm_mass = 0.0;
            F64vec cm_pos = 0.0;
            std::cerr<<"comm_table_.ep_recv_.size()= "<<comm_table_.ep_recv_.size()
                     <<" comm_table_.sp_recv_.size()= "<<comm_table_.sp_recv_.size()
                     <<std::endl;
            for(S32 i=0; i<n_loc_tot_; i++){
                cm_mass += epj_sorted_loc_[i].mass;
                cm_pos += epj_sorted_loc_[i].mass * epj_sorted_loc_[i].pos;
            }
            std::cerr<<"A) cm_mass= "<<cm_mass
                     <<" cm_pos= "<<cm_pos / cm_mass
                     <<std::endl;
            for(S32 i=0; i<comm_table_.ep_recv_.size(); i++){
                //std::cerr<<"i= "<<i<<" mass= "<<comm_table_.ep_recv_[i].mass<<std::endl;
                cm_mass += comm_table_.ep_recv_[i].mass;
                cm_pos += comm_table_.ep_recv_[i].mass * comm_table_.ep_recv_[i].pos;
            }
            std::cerr<<"B) cm_mass= "<<cm_mass
                     <<" cm_pos= "<<cm_pos / cm_mass
                     <<std::endl;
            for(S32 i=0; i<comm_table_.sp_recv_.size(); i++){

                cm_mass += comm_table_.sp_recv_[i].mass;
                cm_pos += comm_table_.sp_recv_[i].mass * comm_table_.sp_recv_[i].pos;
            }
            cm_pos /= cm_mass;
            std::cerr<<"C) cm_mass= "<<cm_mass
                     <<" cm_mass= "<<cm_pos
                     <<std::endl;
        }
        //exit(1);
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
        Comm::barrier();
        const F64 wtime_offset_0 = GetWtime();
        bool flag_retrieve = false;
        const S32 n_proc = Comm::getNumberOfProc();
        S32 ret = 0;
        if(!reuse){
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK0 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif

#if defined(DEBUG_PRING_CALC_FORCE_LOOP_MERGE)
            F64 cm_mass_loc_tmp_00 = 0.0;
            F64vec cm_pos_loc_tmp_00 = 0.0;
            for(S32 i=0; i<psys.getNumberOfParticleLocal(); i++){
                cm_mass_loc_tmp_00 += psys[i].mass;
                #ifdef PHI_R_TREE
                cm_pos_loc_tmp_00.x  += psys[i].mass*psys[i].pos.y*cos(psys[i].pos.x);
                cm_pos_loc_tmp_00.y  += psys[i].mass*psys[i].pos.y*sin(psys[i].pos.x);
                cm_pos_loc_tmp_00.z  += psys[i].mass*psys[i].pos.z;
                #else
                cm_pos_loc_tmp_00  += psys[i].mass*psys[i].pos;
                #endif
            }
            F64vec cm_pos_glb_tmp_00 = Comm::getSum(cm_pos_loc_tmp_00);
            F64 cm_mass_glb_tmp_00 = Comm::getSum(cm_mass_loc_tmp_00);
            cm_pos_loc_tmp_00 /= cm_mass_loc_tmp_00;
            cm_pos_glb_tmp_00 /= cm_mass_glb_tmp_00;
            if(Comm::getRank() == 0){
                std::cerr<<"rank= "<<Comm::getRank()
                         <<" psys.getNumberOfParticleLocal()= "<<psys.getNumberOfParticleLocal()
                         <<" cm_mass_glb_tmp_00= "<<cm_mass_glb_tmp_00
                         <<" cm_pos_glb_tmp_00= "<<cm_pos_glb_tmp_00
                         <<std::endl;
            }
            
            //exit(1);
#endif
            #if 1
            setParticleLocalTreeImpl(psys, epj_org_, true); // new
            #else
            setParticleLocalTreeImpl(psys, epi_org_, epj_org_, true); // new
            #endif
#if defined(DEBUG_PRING_CALC_FORCE_LOOP_MERGE)
            if(Comm::getRank() == 0){
                std::cerr<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
                for(S32 i=0; i<10; i++){
                    std::cerr<<"psys[i].pos= "<<psys[i].pos<<std::endl;
                    std::cerr<<"epj_org_[i].pos= "<<epj_org_[i].pos<<std::endl;
                    std::cerr<<std::endl;
                }
            }
            //exit(1);
#endif
            
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK1 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif


#ifdef PHI_R_TREE
            static const F64 PI = 4.0*atan(1.0);
            F64 len_root_tree = 2.0*PI*1.0001;
            //setRootCell(len_root_tree, F64vec(len_root_tree*0.5, 0.0, 0.0-TREE_Z_COORD_OFFSET) ); // original
            setRootCell(len_root_tree, F64vec(len_root_tree*0.5, 0.0, TREE_Z_COORD_OFFSET) ); // original
            /*
            if(Comm::getRank() < 2){
                std::cerr<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
                std::cerr<<"length_= "<<length_<<std::endl;
                std::cerr<<"center_= "<<center_<<std::endl;
            }
            exit(1);
            */
#else
            setRootCell(dinfo); // original
#endif

            
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK2 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            #if 1
            mortonSortLocalTreeOnlyImpl(epj_org_, epi_sorted_, epj_sorted_loc_); // new
            #else
            mortonSortLocalTreeOnlyImpl(epi_org_, epj_org_, epi_sorted_, epj_sorted_loc_); // new
            #endif
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK3 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif

    #ifdef DEBUG_SORT_LOCAL // for consistency
            mortonSortFP<Tfp>(psys);    // new
    #endif
            
#if defined(DEBUG_PRING_CALC_FORCE_LOOP_MERGE)
            F64 cm_mass_loc_tmp_0 = 0.0;
            F64vec cm_pos_loc_tmp_0 = 0.0;
            for(S32 i=0; i<epj_sorted_loc_.size(); i++){
                cm_mass_loc_tmp_0 += epj_sorted_loc_[i].mass;
                cm_pos_loc_tmp_0  += epj_sorted_loc_[i].mass*epj_sorted_loc_[i].pos;
            }
            F64vec cm_pos_glb_tmp_0 = Comm::getSum(cm_pos_loc_tmp_0);
            F64 cm_mass_glb_tmp_0 = Comm::getSum(cm_mass_loc_tmp_0);
            cm_pos_loc_tmp_0 /= cm_mass_loc_tmp_0;
            cm_pos_glb_tmp_0 /= cm_mass_glb_tmp_0;
            if(Comm::getRank() == 0){
                std::cerr<<"rank= "<<Comm::getRank()
                         <<"epj_sorted_loc_.size()= "<<epj_sorted_loc_.size()
                         <<" cm_mass_glb_tmp_0= "<<cm_mass_glb_tmp_0
                         <<" cm_pos_glb_tmp_0= "<<cm_pos_glb_tmp_0
                         <<" tc_loc_[0].mom_.pos= "<<tc_loc_[0].mom_.pos
                         <<std::endl;
            }
            F64 cm_mass_loc_tmp_1 = 0.0;
            F64vec cm_pos_loc_tmp_1 = 0.0;
            for(S32 i=0; i<psys.getNumberOfParticleLocal(); i++){
                cm_mass_loc_tmp_1 += psys[i].mass;
                #ifdef PHI_R_TREE
                cm_pos_loc_tmp_1.x  += psys[i].mass*psys[i].pos.y*cos(psys[i].pos.x);
                cm_pos_loc_tmp_1.y  += psys[i].mass*psys[i].pos.y*sin(psys[i].pos.x);
                cm_pos_loc_tmp_1.z  += psys[i].mass*psys[i].pos.z;
                #else
                cm_pos_loc_tmp_1  += psys[i].mass*psys[i].pos;
                #endif
            }
            F64vec cm_pos_glb_tmp_1 = Comm::getSum(cm_pos_loc_tmp_1);
            F64 cm_mass_glb_tmp_1 = Comm::getSum(cm_mass_loc_tmp_1);
            cm_pos_loc_tmp_1 /= cm_mass_loc_tmp_1;
            cm_pos_glb_tmp_1 /= cm_mass_glb_tmp_1;
            if(Comm::getRank() == 0){
                std::cerr<<"rank= "<<Comm::getRank()
                         <<" psys.getNumberOfParticleLocal()= "<<psys.getNumberOfParticleLocal()
                         <<" cm_mass_glb_tmp_1= "<<cm_mass_glb_tmp_1
                         <<" cm_pos_glb_tmp_1= "<<cm_pos_glb_tmp_1
                         <<" tc_loc_[0].mom_.pos= "<<tc_loc_[0].mom_.pos
                         <<std::endl;
            }
            
            //exit(1);
#endif
            
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

#if defined(DEBUG_PRING_CALC_FORCE_LOOP_MERGE)
            F64 cm_mass_tmp = 0.0;
            F64vec cm_pos_tmp = 0.0;
            for(S32 i=0; i<epi_sorted_.size(); i++){
                cm_mass_tmp += epi_sorted_[i].mass;
                cm_pos_tmp  += epi_sorted_[i].mass*epi_sorted_[i].pos;
            }
            F64vec cm_pos_glb_tmp2 = Comm::getSum(cm_pos_tmp);
            F64 cm_mass_glb_tmp2 = Comm::getSum(cm_mass_tmp);
            cm_pos_tmp /= cm_mass_tmp;
            cm_pos_glb_tmp2 /= cm_mass_glb_tmp2;
            if(Comm::getRank() % 3 == 0){
                std::cerr<<"rank= "<<Comm::getRank()
                         <<" cm_pos_glb_tmp2= "<<cm_pos_glb_tmp2
                         <<std::endl;
            }
            if( fabs((cm_pos_tmp.x-tc_loc_[0].mom_.pos.x)/tc_loc_[0].mom_.pos.x) > 1e-8
                || fabs((cm_pos_tmp.y-tc_loc_[0].mom_.pos.y)/tc_loc_[0].mom_.pos.y) > 1e-8){
                std::cerr<<"cm_mass_tmp= "<<cm_mass_tmp
                         <<" tc_loc_[0].mom_.mass= "<<tc_loc_[0].mom_.mass
                         <<" cm_pos_tmp= "<<cm_pos_tmp
                         <<" tc_loc_[0].mom_.pos= "<<tc_loc_[0].mom_.pos
    #ifdef PHI_R_TREE
                         <<" tc_loc_[0].mom_.pos_phi= "<<tc_loc_[0].mom_.pos_phi
                         <<" tc_loc_[0].mom_.pos_r= "<<tc_loc_[0].mom_.pos_r
    #endif
                         <<std::endl;
            }
            if(Comm::getRank() == 2){
                std::cerr<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
                std::cerr<<"cm_mass_tmp= "<<cm_mass_tmp<<std::endl;
                std::cerr<<"cm_pos_tmp= "<<cm_pos_tmp<<std::endl;
                std::cerr<<"tc_loc_[0].mom_.mass= "<<tc_loc_[0].mom_.mass<<std::endl;
                std::cerr<<"tc_loc_[0].mom_.pos.x= "<<tc_loc_[0].mom_.pos.x<<std::endl;
                std::cerr<<"tc_loc_[0].mom_.pos.y= "<<tc_loc_[0].mom_.pos.y<<std::endl;
                std::cerr<<"tc_loc_[0].mom_.pos.z= "<<tc_loc_[0].mom_.pos.z<<std::endl;
    #ifdef PHI_R_TREE
                std::cerr<<"tc_loc_[0].mom_.pos_phi= "<<tc_loc_[0].mom_.pos_phi<<std::endl;
                std::cerr<<"tc_loc_[0].mom_.pos_r= "<<tc_loc_[0].mom_.pos_r<<std::endl;
    #endif
                for(S32 i=8; i<8+8; i++){
                    std::cerr<<"tc_loc_["<<i<<"].mom_.mass= "<<tc_loc_[i].mom_.mass<<std::endl;
                    std::cerr<<"tc_loc_["<<i<<"].mom_.pos(cartesian)= "<<tc_loc_[i].mom_.pos<<std::endl;
    #ifdef PHI_R_TREE
                    std::cerr<<"tc_loc_["<<i<<"].mom_.pos(cylindrical)= "<<tc_loc_[i].mom_.pos_phi
                             <<"   "<<tc_loc_[i].mom_.pos_r<<"   "<<tc_loc_[i].mom_.pos.z<<std::endl;
    #endif
                }
            }
            //exit(1);
#endif

            
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK6 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            addMomentAsSpLocalTreeImpl(typename TSM::force_type()); // original
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK7 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            
#ifdef USE_SUPER_DOMAIN
            //exchangeLocalEssentialTreeSuperDomain(dinfo); // original
            exchangeLocalEssentialTreeSuperDomain2(dinfo); // original
#else
            exchangeLocalEssentialTree3(dinfo); // original
#endif
#if defined(DEBUG_PRING_CALC_FORCE_LOOP_MERGE)
            F64 cm_mass_glb_tmp = cm_mass_tmp;
            F64vec cm_pos_glb_tmp = cm_mass_tmp*cm_pos_tmp;
            for(S32 i=0; i<comm_table_.ep_recv_.size(); i++){
                cm_mass_glb_tmp += comm_table_.ep_recv_[i].mass;
                cm_pos_glb_tmp += comm_table_.ep_recv_[i].mass*comm_table_.ep_recv_[i].pos;
            }
            for(S32 i=0; i<comm_table_.sp_recv_.size(); i++){
                cm_mass_glb_tmp += comm_table_.sp_recv_[i].mass;
                cm_pos_glb_tmp += comm_table_.sp_recv_[i].mass*comm_table_.sp_recv_[i].pos;
            }
            cm_pos_glb_tmp /= cm_mass_glb_tmp;
            if(Comm::getRank()==0){
                std::cerr<<"cm_mass_glb_tmp= "<<cm_mass_glb_tmp
                         <<" cm_pos_glb_tmp= "<<cm_pos_glb_tmp
                         <<std::endl;
            }
            //exit(1);
#endif
            
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


#if defined(DEBUG_PRING_CALC_FORCE_LOOP_MERGE)
            if(Comm::getRank() == 0){
                std::cerr<<"Comm::getRank()= "<<Comm::getRank()<<std::endl;
                std::cerr<<"tc_glb_[0].mom_.mass= "<<tc_glb_[0].mom_.mass<<std::endl;
                std::cerr<<"tc_glb_[0].mom_.pos.x= "<<tc_glb_[0].mom_.pos.x<<std::endl;
                std::cerr<<"tc_glb_[0].mom_.pos.y= "<<tc_glb_[0].mom_.pos.y<<std::endl;
                std::cerr<<"tc_glb_[0].mom_.pos.z= "<<tc_glb_[0].mom_.pos.z<<std::endl;
    #ifdef PHI_R_TREE                
                std::cerr<<"tc_glb_[0].mom_.pos_phi= "<<tc_glb_[0].mom_.pos_phi<<std::endl;
                std::cerr<<"tc_glb_[0].mom_.pos_r= "<<tc_glb_[0].mom_.pos_r<<std::endl;
    #endif
                for(S32 i=8; i<8+8; i++){
                    std::cerr<<"tc_glb_["<<i<<"].mom_.mass= "<<tc_glb_[i].mom_.mass<<std::endl;
                    std::cerr<<"tc_glb_["<<i<<"].mom_.pos(cartesian)= "<<tc_glb_[i].mom_.pos<<std::endl;
    #ifdef PHI_R_TREE
                    std::cerr<<"tc_glb_["<<i<<"].mom_.pos(cylindrical)= "<<tc_glb_[i].mom_.pos_phi
                             <<"   "<<tc_glb_[i].mom_.pos_r<<"   "<<tc_glb_[i].mom_.pos.z<<std::endl;
    #endif
                }
            }
            //exit(1);
#endif
            
            
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
            if(Comm::getRank()==0)std::cerr<<"OK14.1 @calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe"<<std::endl;
#endif
            
            if(adr_n_walk!=NULL) delete [] adr_n_walk;
            adr_n_walk = new S32[ipg_.size()];
            n_epi_recorder_for_force_[0].resizeNoInitialize(ipg_.size());
            for(S32 i=0; i<ipg_.size(); i++){
                n_epi_recorder_for_force_[0][i] = ipg_[i].n_ptcl_;
            }

            cpe_pars.n_walk_cpe  = n_walk_cpe;
            cpe_pars.n_disp_walk = n_disp_walk;
            cpe_pars.adr_n_walk  = adr_n_walk;
    #ifdef BALANCE_NWALK
            cpe_pars.n_sat       = N_SATELLITE;
            cpe_pars.n_epi       = n_epi_recorder_for_force_[0].getPointer();
            cpe_pars.n_epj       = n_epj_recorder_for_force_[0].getPointer();
            cpe_pars.n_spj       = n_spj_recorder_for_force_[0].getPointer();
            cpe_pars.n_walk      = ipg_.size();

        #ifdef BALANCE_NWALK_OUT
            std::ofstream flog;
            std::stringstream fnum;
            std::string fname;
            fnum << std::setfill('0') << std::setw(5) << Comm::getRank();
            fname = "catch" + fnum.str() + ".txt";
            //flog.open(fname.c_str(),std::ios::trunc);
            flog.open(fname.c_str(), std::ios::app);
            flog<<"rank="<<Comm::getRank()<<" n="<<ipg_.size()
                <<" nepi="<<n_epi_recorder_for_force_[0].size()
                <<" nepj="<<n_epj_recorder_for_force_[0].size()
                <<" nspj="<<n_spj_recorder_for_force_[0].size()
                <<std::endl;
        #endif //BALANCE_NWALK_OUT

            __real_athread_create(0,(void*)slave_balance_nwalk,(void*)&cpe_pars);
            athread_wait(0);

        #ifdef BALANCE_NWALK_OUT
            for(int i=0; i<64; i++){
                //flog<<"i= "<<i<<"  "<<adr_n_walk[i]<<std::endl;
                flog<<"i= "<<i<<"  "<<n_walk_cpe[i]<<std::endl;
            }
            flog.close();
        #endif
            
        #ifdef BW_CHECK_SAFETY
            long ncost[64]={0},nwalk_cpe=0,ni_cpe[64]={0},nj_cpe[64]={0},ns_cpe[64]={0},sum=0,ni_sum=0, nj_sum=0, ns_sum=0;
            for (int i=0; i<64; i++) {
              nwalk_cpe += n_walk_cpe[i];
              for (int j=0; j<n_walk_cpe[i]; j++) {
                int iw = adr_n_walk[n_disp_walk[i]+j];
                ni_cpe[i] += cpe_pars.n_epi[iw];
                nj_cpe[i] += cpe_pars.n_epj[iw];
                ns_cpe[i] += cpe_pars.n_spj[iw];
                ncost[i] += cpe_cost(cpe_pars.n_epi[iw],cpe_pars.n_epj[iw],cpe_pars.n_spj[iw],cpe_pars.n_sat,54,30,100,32,237404);
              }
              sum += ncost[i];
              ni_sum += ni_cpe[i];
              nj_sum += nj_cpe[i];
              ns_sum += ns_cpe[i];
            }
            long sum_mpe=0,ni_mpe=0,nj_mpe=0,ns_mpe=0;
            for(int k=0;k<ipg_.size();k++) {
              sum_mpe += cpe_cost(cpe_pars.n_epi[k], cpe_pars.n_epj[k], cpe_pars.n_spj[k], cpe_pars.n_sat, 54,30,100,32,237404);
              ni_mpe += cpe_pars.n_epi[k];
              nj_mpe += cpe_pars.n_epj[k];
              ns_mpe += cpe_pars.n_spj[k];
            }
            //assert(sum==sum_mpe);
            if(sum!=sum_mpe) {
              if(Comm::getRank()==0) {
              std::cerr<<"Ncost sum not match! CPE = "<<sum<<" MPE = "<<sum_mpe<<std::endl;
              std::cerr<<"Nw sum CPE="<<nwalk_cpe<<" MPE="<<ipg_.size()<<std::endl;
              std::cerr<<"NI sum CPE="<<ni_sum<<" MPE="<<ni_mpe<<std::endl;
              std::cerr<<"NJ sum CPE="<<nj_sum<<" MPE="<<nj_mpe<<std::endl;
              std::cerr<<"NS sum CPE="<<ns_sum<<" MPE="<<ns_mpe<<std::endl;

              fprintf(stderr,"Balance walk check on MPE:\n");
              fprintf(stderr,"CID         total");
              for(int k=0;k<64;k++) fprintf(stderr,"%12d",k);
              fprintf(stderr,"\n");

              int nwtot=0;
              for (int i=0; i<64; i++) nwtot += n_walk_cpe[i];
              fprintf(stderr,"NW_CPE%9d",nwtot);
              for(int k=0;k<64;k++) fprintf(stderr,"%12d",n_walk_cpe[k]);
              fprintf(stderr,"\nNW_DISP          ");
              for(int k=0;k<64;k++) fprintf(stderr,"%12d",n_disp_walk[k]);
              fprintf(stderr,"\nNW_ADR_FIRST     ");
              for(int k=0;k<64;k++) fprintf(stderr,"%12d",adr_n_walk[n_disp_walk[k]]);
              fprintf(stderr,"\n");
              fprintf(stderr,"SUM_I%12ld",ni_sum);
              for(int k=0;k<64;k++) fprintf(stderr,"%12d",ni_cpe[k]);
              fprintf(stderr,"\n");
              fprintf(stderr,"SUM%14ld",sum);
              for(int k=0;k<64;k++) fprintf(stderr,"%12d",ncost[k]);
              fprintf(stderr,"\n");
              }
//#ifdef BALANCE_NWALK_OUT
//              flog.open(fname.c_str(),std::ios::trunc);
//              for(int i=0; i<64; i++) flog<<"["<<i<<"]="<<ncost[i]<<" ";
//              flog.close();
//#endif
            }
        #endif //BW_CHECK_SAFETY
    #else //BALANCE_NWALK
        #ifdef BALANCE_NWALK_ON_MPE
            balanceNwalk();
        #else
            int nw_per_cpe = (ipg_.size()+63)/64;
            int nw_cutoff = ipg_.size()%64;
            if (nw_cutoff==0) nw_cutoff = 64;
            for(int i=0; i<64; i++) n_walk_cpe[i]  = (i<nw_cutoff)?nw_per_cpe:nw_per_cpe-1;
            n_disp_walk[0] = 0;
            for(int i=1; i<64; i++) n_disp_walk[i] = n_disp_walk[i-1] + n_walk_cpe[i-1];
            for(int i=0; i<ipg_.size(); i++) adr_n_walk[i] = i;
        #endif
    #endif //BALANCE_NWALK


#if defined(DEBUG_PRING_CALC_FORCE_LOOP_MERGE)
            Comm::barrier();
            if(Comm::getRank()==0){
                std::cerr<<"XXXXX"<<std::endl;
                DumpInteractionList(epi_sorted_, ipg_,
                                    epj_sorted_, id_epj_recorder_for_force_[0], n_epj_recorder_for_force_[0],
                                    n_disp_epj_recorder_for_force_[0],
                                    spj_sorted_, id_spj_recorder_for_force_[0], n_spj_recorder_for_force_[0],
                                    n_disp_spj_recorder_for_force_[0], 20);
            }
            Comm::barrier();
#endif
            
            
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK15 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            ret += calcForceUsingIdListMultiWalkIndex3(pfunc_dispatch, pfunc_retrieve, n_walk_limit, mw_info_, clear_force, reuse, flag_retrieve); // new
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK16 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif

#ifdef CALC_MOMENT_LEAN
            // TODO: move to CPE
            etc_loc_.resizeNoInitialize(adr_tc_level_partition_loc_[lev_max_loc_+1]);
            etc_glb_.resizeNoInitialize(adr_tc_level_partition_glb_[lev_max_glb_+1]);
            #ifdef SUNWAY
            unsigned long arg[3];
            arg[0] = (unsigned long)(adr_tc_level_partition_loc_[lev_max_loc_+1]);
            arg[1] = (unsigned long)(tc_loc_.getPointer());
            arg[2] = (unsigned long)(etc_loc_.getPointer());
            __real_athread_spawn((void*)slave_CopyTCToETC, arg);
            athread_join();
            
            arg[0] = (unsigned long)(adr_tc_level_partition_glb_[lev_max_glb_+1]);
            arg[1] = (unsigned long)(tc_glb_.getPointer());
            arg[2] = (unsigned long)(etc_glb_.getPointer());
            __real_athread_spawn((void*)slave_CopyTCToETC, arg);
            athread_join();
            
            #else
            for(S32 i=0; i<adr_tc_level_partition_loc_[lev_max_loc_+1]; i++){
                etc_loc_[i].adr_tc_ = tc_loc_[i].adr_tc_;
                etc_loc_[i].adr_ptcl_ = tc_loc_[i].adr_ptcl_;
                etc_loc_[i].n_ptcl_ = tc_loc_[i].n_ptcl_;
                etc_loc_[i].level_ = tc_loc_[i].level_;
            }
            #endif
#endif
#if defined(DEBUG_PRING_CALC_FORCE_LOOP_MERGE)
            Comm::barrier();
            TreeCellCheck(tc_loc_,  "tc_loc_,  first_loop");
            Comm::barrier();
            TreeCellCheck(etc_loc_, "etc_loc_, first_loop, mass and pos are not crrect becaus these values are update later");
            Comm::barrier();
            //exit(1);
#endif
        }
        else{
            //////////////////////////////////
            /////////// REUSE LIST /////////// 
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK17 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            F64 wtime_offset = GetWtime();
#ifdef SUNWAY
            unsigned long arg[5];
    #ifdef PHI_R_TREE
            epj_sorted_loc_.resizeNoInitialize(n_loc_tot_);
            arg[0] = (unsigned long)(n_loc_tot_);
            arg[1] = (unsigned long)(epi_sorted_.getPointer());
            arg[2] = (unsigned long)(epj_sorted_loc_.getPointer());
            __real_athread_spawn((void*)slave_CopyEPISortedToEPJSortedLoc, arg);
            athread_join();
    #else
            epj_sorted_loc_.resizeNoInitialize(n_loc_tot_);
            arg[0] = (unsigned long)(n_loc_tot_);
            arg[1] = (unsigned long)(epi_sorted_.getPointer());
            arg[2] = (unsigned long)(epj_sorted_loc_.getPointer());
            arg[3] = (unsigned long)(sizeof(epj_sorted_loc_[0]));
            __real_athread_spawn((void*)slave_CopyDirect, arg);
            athread_join();


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
    #endif
#else
            for(S32 i=0; i<n_loc_tot_; i++){
                epj_sorted_loc_[i].copyFromEPI(epi_sorted_[i]);
            }
#endif
            time_profile_.set_particle_local_tree+= GetWtime() - wtime_offset;
            wtime_t0 = GetWtime() - wtime_offset;
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            for(S32 i=0; i<n_loc_tot_; i++) assert(epj_sorted_loc_[i].id == epi_sorted_[i].id); // for debug
            if(Comm::getRank()==0)std::cerr<<"OK18 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
            
#if defined(DEBUG_PRING_CALC_FORCE_LOOP_MERGE)
            Comm::barrier();
            PtclCheck(epi_sorted_, "epi_sorted_, reuse", 100, 0);
            PtclCheck(epj_sorted_loc_, "epj_sorted_loc_, reuse", 100, 0);
            Comm::barrier();
            //exit(1);
#endif
            
            wtime_t1 = GetWtime() - wtime_offset;
            
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK19 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;

#endif
            wtime_offset = GetWtime();

#ifdef CALC_MOMENT_LEAN
            //ReallocatableArray<etcLM> etc_loc;
            //etc_loc.resizeNoInitialize(tc_loc_.size());
            //etc_loc.resizeNoInitialize(adr_tc_level_partition_loc_[lev_max_loc_+1]);
            CalcMomentLean(adr_tc_level_partition_loc_, etc_loc_.getPointer(), epj_sorted_loc_.getPointer(), lev_max_loc_, n_leaf_limit_);

            //exit(1);
#else
            CalcMoment(adr_tc_level_partition_loc_, tc_loc_.getPointer(), epj_sorted_loc_.getPointer(), lev_max_loc_, n_leaf_limit_);
#endif
            time_profile_.calc_moment_local_tree += GetWtime() - wtime_offset;
            wtime_t2 = GetWtime() - wtime_offset;
            
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK20 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
#if defined(DEBUG_PRING_CALC_FORCE_LOOP_MERGE)
            Comm::barrier();
            TreeCellCheck(etc_loc_, "etc_loc_, after calc moment");
            Comm::barrier();
            //exit(1);
#endif
            
            wtime_offset = GetWtime();
#ifdef CALC_MOMENT_LEAN
            addMomentAsSpLocalTreeLeanImpl(typename TSM::force_type()); //original
#else
            addMomentAsSpLocalTreeImpl(typename TSM::force_type()); //original
#endif
            wtime_t3 = GetWtime() - wtime_offset;

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK21 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
            //athread_halt();
            //Finalize();
            //std::exit(0);
#endif
            wtime_offset = GetWtime();
#ifdef USE_SUPER_DOMAIN
            //exchangeLocalEssentialTreeSuperDomain(dinfo, true);
            exchangeLocalEssentialTreeSuperDomain2(dinfo, true);
#else
            exchangeLocalEssentialTree3(dinfo, true); //original
#endif

#if defined(DEBUG_PRING_CALC_FORCE_LOOP_MERGE)
            Comm::barrier();
            PtclCheck(comm_table_.ep_recv_, "ep_recv, reuse", 100, 0);
            PtclCheck(comm_table_.sp_recv_, "sp_recv, reuse", 100, 0);
            Comm::barrier();
            //exit(1);
#endif
            
            wtime_t4 = GetWtime() - wtime_offset;

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK22 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
            //athread_halt();
            //Finalize();
            //std::exit(0);
#endif
            wtime_offset = GetWtime();
            
#ifdef USE_SUPER_DOMAIN
            const S32 n_epj_add = comm_table_.ep_recv_.size();
            const S32 n_spj_add = comm_table_.sp_recv_.size();
#else
            const S32 n_epj_add = comm_table_.n_disp_ep_recv_[n_proc];
            const S32 n_spj_add = comm_table_.n_disp_sp_recv_[n_proc];
#endif
            
#ifdef SUNWAY

    #if 1
        #if VERSION_MSORT_GLB_TREE_REUSE == 0
            arg[0] = (unsigned long)(epj_sorted_.size());
            arg[1] = (unsigned long)(epj_sorted_loc_.getPointer());
            arg[2] = (unsigned long)(epj_sorted_.getPointer());
            arg[3] = (unsigned long)(adr_epj_glb2loc_.getPointer());
            __real_athread_spawn((void*)slave_CopyEPJLocToEPJGlb, arg);
            athread_join();
        #elif VERSION_MSORT_GLB_TREE_REUSE == 1
            int n_groups[NUMBER_OF_CPE];
            unsigned long adr_adr_epj_loc_group_head_[NUMBER_OF_CPE];
            unsigned long adr_adr_epj_glb_group_head_[NUMBER_OF_CPE];
            unsigned long adr_group_size_epj_loc_[NUMBER_OF_CPE];
            for (S32 i=0; i<NUMBER_OF_CPE; i++) {
                n_groups[i] = adr_epj_loc_group_head_[i].size();
                adr_adr_epj_loc_group_head_[i] = (unsigned long) adr_epj_loc_group_head_[i].getPointer();
                adr_adr_epj_glb_group_head_[i] = (unsigned long) adr_epj_glb_group_head_[i].getPointer();
                adr_group_size_epj_loc_[i]     = (unsigned long) group_size_epj_loc_[i].getPointer();
            }
            arg[0] = Comm::getRank();
            arg[1] = (unsigned long)(epj_sorted_loc_.size());
            arg[2] = (unsigned long)(epj_sorted_loc_.getPointer());
            arg[3] = (unsigned long)(epj_sorted_.getPointer());
            arg[4] = (unsigned long)(n_groups);
            arg[5] = (unsigned long)(adr_adr_epj_loc_group_head_);
            arg[6] = (unsigned long)(adr_adr_epj_glb_group_head_);
            arg[7] = (unsigned long)(adr_group_size_epj_loc_);
            __real_athread_spawn((void*)slave_CopyEPJLocToEPJGlb, arg);
            athread_join();
        #elif VERSION_MSORT_GLB_TREE_REUSE == 2
            int n_groups[NUMBER_OF_CPE];
            unsigned long adr_adr_epj_glb_group_head_[NUMBER_OF_CPE];
            unsigned long adr_adr_epj_loc_group_head_[NUMBER_OF_CPE];
            unsigned long adr_group_size_epj_glb_[NUMBER_OF_CPE];
            for (S32 i=0; i<NUMBER_OF_CPE; i++) {
                n_groups[i] = adr_epj_glb_group_head_[i].size();
                adr_adr_epj_glb_group_head_[i] = (unsigned long) adr_epj_glb_group_head_[i].getPointer();
                adr_adr_epj_loc_group_head_[i] = (unsigned long) adr_epj_loc_group_head_[i].getPointer();
                adr_group_size_epj_glb_[i]     = (unsigned long) group_size_epj_glb_[i].getPointer();
            }
            arg[0] = Comm::getRank();
            arg[1] = (unsigned long)(epj_sorted_.size());
            arg[2] = (unsigned long)(epj_sorted_loc_.getPointer());
            arg[3] = (unsigned long)(epj_sorted_.getPointer());
            arg[4] = (unsigned long)(n_groups);
            arg[5] = (unsigned long)(adr_adr_epj_glb_group_head_);
            arg[6] = (unsigned long)(adr_adr_epj_loc_group_head_);
            arg[7] = (unsigned long)(adr_group_size_epj_glb_);
            //if (Comm::getRank() == 0) {//debug
            __real_athread_spawn((void*)slave_CopyEPJLocToEPJGlb, arg);
            athread_join();
            //}// debug
            #if 0
            //* For test run
            Comm::barrier();
            if (Comm::getRank() == 0) std::cout << "... copy ends. ..." << std::endl;
            athread_halt(); 
            Finalize();
            std::exit(0); 
            #endif
        #else // VERSION_MSORT_GLB_TREE_REUSE
            #error The value of the macro `VERSION_MSORT_GLB_TREE_REUSE` is incorrect.
        #endif // VERSION_MSORT_GLB_TREE_REUSE
            
            /*
            arg[0] = (unsigned long)(epj_sorted_.size());
            arg[1] = (unsigned long)(epj_sorted_loc_.getPointer());
            arg[2] = (unsigned long)(epj_sorted_.getPointer());
            arg[3] = (unsigned long)(adr_epj_glb2loc_.getPointer());
            __real_athread_spawn((void*)slave_CopyEPJLocToEPJ, arg);
            athread_join();
            */
    #else //#if 1

            //--- Old implementation [start] -----
            //for(S32 i=0; i<n_loc_tot_; i++)  epj_sorted_[i].id = -i+1; // for debug
            arg[0] = (unsigned long)(n_loc_tot_);
            arg[1] = (unsigned long)(epj_sorted_loc_.getPointer());
            arg[2] = (unsigned long)(epj_sorted_.getPointer());
            arg[3] = (unsigned long)(sizeof(epj_sorted_loc_[0]));
            arg[4] = (unsigned long)(adr_epj_loc2glb_.getPointer());
            __real_athread_spawn((void*)slave_CopyIndirect, arg);
            athread_join();
            //--- Old implementation [end] -----
    #endif //#if 1

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
#else // SUNWAY
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
#endif // SUNWAY
            time_profile_.morton_sort_global_tree += GetWtime() - wtime_offset;
            
            wtime_t5 = GetWtime() - wtime_offset;

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK26 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif

#if defined(DEBUG_PRING_CALC_FORCE_LOOP_MERGE)
            Comm::barrier();
            PtclCheck2(epj_sorted_, spj_sorted_, tp_glb_, "before moment glb, reuse");
            Comm::barrier();
            //exit(1);
#endif
            
            wtime_offset = GetWtime();
            
#ifdef CALC_MOMENT_LEAN
            CalcMomentLongGlobalTreeLean
                (adr_tc_level_partition_glb_,  etc_glb_.getPointer(), 
                 tp_glb_.getPointer(),         epj_sorted_.getPointer(),
                 spj_sorted_.getPointer(),     lev_max_glb_,
                 n_leaf_limit_);
            time_profile_.calc_moment_global_tree += GetWtime() - wtime_offset;
            //exit(1);
#else
            //calcMomentGlobalTreeOnly(); //original
            CalcMomentLongGlobalTree
                (adr_tc_level_partition_glb_,  tc_glb_.getPointer(), 
                 tp_glb_.getPointer(),         epj_sorted_.getPointer(),
                 spj_sorted_.getPointer(),     lev_max_glb_,
                 n_leaf_limit_);
            time_profile_.calc_moment_global_tree += GetWtime() - wtime_offset;
            //exit(1);
#endif
            wtime_t6 = GetWtime() - wtime_offset;

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK27 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif

#if defined(DEBUG_PRING_CALC_FORCE_LOOP_MERGE)
            Comm::barrier();
            TreeCellCheck(etc_glb_, "etc_glb_, after calc moment");
            Comm::barrier();
            //exit(1);
#endif
            
            wtime_offset = GetWtime();
#ifdef CALC_MOMENT_LEAN            
            addMomentAsSpGlobalTreeLeanImpl(typename TSM::force_type()); //original
#else
            addMomentAsSpGlobalTreeImpl(typename TSM::force_type()); //original
#endif
            wtime_t7 = GetWtime() - wtime_offset;

#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK28 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif

#if defined(DEBUG_PRING_CALC_FORCE_LOOP_MERGE)
            Comm::barrier();
            if(Comm::getRank()==0){
                std::cerr<<"YYYYY"<<std::endl;
                DumpInteractionList(epi_sorted_, ipg_,
                                    epj_sorted_, id_epj_recorder_for_force_[0], n_epj_recorder_for_force_[0],
                                    n_disp_epj_recorder_for_force_[0],
                                    spj_sorted_, id_spj_recorder_for_force_[0], n_spj_recorder_for_force_[0],
                                    n_disp_spj_recorder_for_force_[0], 20);
            }
            Comm::barrier();
#endif
            
            wtime_offset = GetWtime();
            ret += calcForceUsingIdListMultiWalkIndex3(pfunc_dispatch, pfunc_retrieve, n_walk_limit, mw_info_, clear_force, reuse, flag_retrieve); // new
            wtime_t8 = GetWtime() - wtime_offset;
#ifdef PARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0)std::cerr<<"OK29 @calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge"<<std::endl;
#endif
        }
        //NAN_TEST(epi_sorted_,"[nan-1]");
        //NAN_TEST(epj_sorted_,"[nan-1]");

#if defined(DEBUG_PRING_CALC_FORCE_LOOP_MERGE)
        Comm::barrier();
        PtclCheck(epi_sorted_, "epi_sorted_, after integration", 100000);
        Comm::barrier();
        if(reuse){ exit(1);}
#endif
        
        time_profile_.calc_force_all += GetWtime() - wtime_offset_0;
        return ret;
    }
}
