#pragma once

//#define UNORDERED_SET

#ifdef UNORDERED_SET 
#include<unordered_set>
#else
#include<set>
#endif

#include<sort.hpp>
#include<tree.hpp>
#include<comm_table.hpp>
#include<multi_walk.hpp>
#include<tree_for_force_utils.hpp>

extern "C"{
    #include "cpe_func.h"
}

namespace ParticleSimulator{
    ///////////////////////////////////////
    /// TREE FOR FORCE CLASS DEFINITION ///
    template<
        class TSM, // search mode
        class Tforce, // USER def
        class Tepi, // USER def
        class Tepj, // USER def
        class Tmomloc, // PS or USER def
        class Tmomglb, // PS or USER def
        class Tspj // PS or USER def
        >
    class TreeForForce{

    public:
        //F64 length_; // length of a side of the root cell
        //F64vec center_; // new member (not used)
        //F64ort pos_root_cell_;
        CommTable<Tepj, Tspj> comm_table_; // for reuse the list
        ReallocatableArray< Tepi > epi_sorted_;
        ReallocatableArray< Tepj > epj_sorted_loc_;
    private:

        //CommTable<Tepj, Tspj> comm_table_; // for reuse the list

        TimeProfile time_profile_;

        CountT n_interaction_ep_ep_local_, n_interaction_ep_sp_local_, n_walk_local_;
        CountT n_let_ep_send_1st_, n_let_ep_recv_1st_, n_let_sp_send_1st_, n_let_sp_recv_1st_,
            n_let_ep_send_2nd_, n_let_ep_recv_2nd_;
        CountT * n_cell_open_;
        CountT n_proc_send_exchange_LET_1st__icomm_sp_, n_proc_recv_exchange_LET_1st__icomm_sp_, 
            n_proc_send_exchange_LET_1st__icomm_ep_, n_proc_recv_exchange_LET_1st__icomm_ep_; 

        F64 wtime_exlet_comm_;
        F64 wtime_exlet_a2a_;
        F64 wtime_exlet_a2av_;
        F64 Tcomm_scatterEP_tmp_;
        F64 wtime_walk_LET_1st_, wtime_walk_LET_2nd_;

        bool is_initialized_;

        S64 n_interaction_ep_ep_;
        S64 n_interaction_ep_sp_;
        S32 ni_ave_;
        S32 nj_ave_;
        S32 n_epi_loc_ave_;
        S32 n_epj_loc_ave_;
        S32 n_spj_loc_ave_;
        
        RadixSort<U64, 8> rs_;
        S32 n_loc_tot_; // # of all kinds of particles in local process
        S64 n_glb_tot_; // # of all kinds of particles in all processes
        S32 n_leaf_limit_;
        S32 n_group_limit_;

        //S32 adr_tc_level_partition_[TREE_LEVEL_LIMIT+2];

        //S32 lev_max_;
        F64 theta_;
        
        F64 length_; // length of a side of the root cell
        F64vec center_; // new member (not used)
        F64ort pos_root_cell_;


        
#ifdef USE_MEMORY_POOL
        TempArray< Tepj >epj_org_;
        TempArray< TreeParticle > tp_buf_;
        TempArray< TreeParticle > tp_buf_2_;
        TempArray< Tspj > spj_org_;
        TempArray< TreeCell< Tmomloc > > tc_loc_;
        TempArray< TreeCell< Tmomglb > > tc_glb_;
    #ifdef REMOVE_TP_LOC
        ReallocatableArray< TreeParticle > tp_glb_;
        TempArray< U64 > tp_loc_adr_ptcl_;
        //TempArray< U32 > tp_loc_adr_ptcl_;
    #else // REMOVE_TP_LOC
        ReallocatableArray< TreeParticle > tp_loc_, tp_glb_;
    #endif // REMOVE_TP_LOC
        //ReallocatableArray< TreeCell< Tmomloc > > tc_loc_;
        //ReallocatableArray< TreeCell< Tmomglb > > tc_glb_;
        ReallocatableArray< etcLM > etc_loc_;
        ReallocatableArray< etcLM > etc_glb_;
        ReallocatableArray< Tepj > epj_sorted_;
        ReallocatableArray< Tspj > spj_sorted_, spj_sorted_loc_; 
        
#else
#ifdef REMOVE_TP_LOC
        ReallocatableArray< TreeParticle > tp_buf_, tp_glb_;
#else
        ReallocatableArray< TreeParticle > tp_buf_, tp_loc_, tp_glb_;
#endif
        ReallocatableArray< TreeCell< Tmomloc > > tc_loc_;
        ReallocatableArray< TreeCell< Tmomglb > > tc_glb_;
        ReallocatableArray< Tepj > epj_sorted_, epj_org_;
        ReallocatableArray< Tspj > spj_sorted_, spj_org_;
        ReallocatableArray< Tspj > spj_sorted_loc_; // for reuse of interaction list
        ReallocatableArray<etcLM> etc_loc_;
        ReallocatableArray<etcLM> etc_glb_;
#endif
        ReallocatableArray< IPGroup<TSM> > ipg_;
        ReallocatableArray<Tepj> epj_send_;
        ReallocatableArray<Tepj> epj_recv_;
        ReallocatableArray<Tspj> spj_send_;
        ReallocatableArray<Tspj> spj_recv_;

        S32 * n_ep_send_; // * n_proc
        S32 * n_sp_send_; // * n_proc
        S32 * n_ep_send_disp_; // * n_proc+1
        S32 * n_sp_send_disp_; // * n_proc+1
        S32 * n_ep_recv_; // * n_proc
        S32 * n_sp_recv_; // * n_proc
        S32 * n_ep_recv_disp_; // * n_proc+1
        S32 * n_sp_recv_disp_; // * n_proc+1

        S32 * n_ep_sp_send_; // 2 * n_proc: even id is # of EP and odd id is # of SP
        S32 * n_ep_sp_recv_; // 2 * n_proc: even id is # of EP and odd id is # of SP

        
        ReallocatableArray<S32> * id_ep_send_buf_;
        ReallocatableArray<S32> * id_sp_send_buf_;

        S32 ** id_proc_send_; // id_proc_send_[n_thread][n_proc]
        ReallocatableArray<S32> * rank_proc_send_; // new version of id_proc_send. have to replace

        ReallocatableArray<Tforce> force_sorted_;
        ReallocatableArray<Tforce> force_org_;

        ReallocatableArray<Tepj> * epj_for_force_;
        ReallocatableArray<Tspj> * spj_for_force_;

        ReallocatableArray<S32> * adr_epj_for_force_;
        ReallocatableArray<S32> * adr_spj_for_force_;

        ReallocatableArray<S32> * id_epj_recorder_for_force_;
        ReallocatableArray<S32> * id_spj_recorder_for_force_;
    public:
        ReallocatableArray<S32> * n_epi_recorder_for_force_;
        ReallocatableArray<S32> * n_disp_epi_recorder_for_force_;
        ReallocatableArray<S32> * n_epj_recorder_for_force_;
        ReallocatableArray<S32> * n_spj_recorder_for_force_;
        ReallocatableArray<S32> * n_disp_epj_recorder_for_force_;
        ReallocatableArray<S32> * n_disp_spj_recorder_for_force_;
    private:
        ReallocatableArray<S32>   adr_epj_loc2glb_;

        //---- the followings are added by DN
        ReallocatableArray<S32>   adr_epj_glb2loc_;
        ReallocatableArray<S32> * adr_epj_loc_group_head_;
        ReallocatableArray<S32> * adr_epj_glb_group_head_;
        ReallocatableArray<S32> * group_size_epj_loc_;
        ReallocatableArray<S32> * group_size_epj_glb_;
        //----
        
        //ReallocatableArray<S32>   adr_epj_glb2loc_; // added by DN
        ReallocatableArray<S32>   adr_epj_buf2glb_;
        ReallocatableArray<S32>   adr_spj_buf2glb_;
        ReallocatableArray<S32>   adr_epj_org2glb_;
        //ReallocatableArray<U32>   adr_tp_loc_org2sorted_;
        ReallocatableArray<U32>   adr_ptcl_of_tp_loc_;
        ReallocatableArray<IdRecorderForInteraction> id_recorder_for_interaction_;
    public:
        S32 adr_tc_level_partition_loc_[TREE_LEVEL_LIMIT+2];
        S32 adr_tc_level_partition_glb_[TREE_LEVEL_LIMIT+2];
        S32 lev_max_loc_;
        S32 lev_max_glb_;
    private:
        ReallocatableArray<F64ort> tree_outer_pos_;
        
        S32 n_surface_for_comm_;

        // new variables for commnuication of LET
        // for scatterEP
        ReallocatableArray<Tepj> * ep_send_buf_for_scatter_;
        ReallocatableArray<F64vec> * shift_image_domain_;


        //PROFILE::Profile profile;

        // for gather mode
        class EPJWithR{
            Tepj epj_;
            F64 r_search_;
        public:
            Tepj getEPJ() const { return epj_; }
            F64vec getPos() const { return epj_.getPos(); }
            F64 getCharge() const { return epj_.getCharge(); }
            F64 getRSearch() const { return r_search_; }
            void copyFromEPJ(const Tepj & epj){ epj_ = epj; }
            void setRSearch(const F64 r_search){ r_search_ = r_search; }
            void setPos(const F64vec & pos){ epj_.setPos(pos);}
        };
        ReallocatableArray<EPJWithR> epjr_sorted_; // cirectly copied from EPJ + RSearch
        ReallocatableArray<EPJWithR> epjr_send_;
        ReallocatableArray<EPJWithR> epjr_recv_;
        ReallocatableArray<EPJWithR> epjr_recv_1st_buf_;
        ReallocatableArray<EPJWithR> epjr_recv_2nd_buf_;
        ReallocatableArray<EPJWithR> * epjr_send_buf_; // for 1st communication
        ReallocatableArray<EPJWithR> * epjr_send_buf_for_scatter_;
        ReallocatableArray<EPJWithR> * epjr_recv_1st_sorted_;

        MultiWalkInfo<Tepi, Tepj, Tspj, Tforce> mw_info_;

        // for symmeteric mode
        ReallocatableArray<Tepj> epj_recv_1st_buf_;
        ReallocatableArray<Tepj> epj_recv_2nd_buf_;
        ReallocatableArray<Tepj> * epj_send_buf_; // for 1st communication
        S32 * n_epj_recv_1st_;
        S32 * n_epj_recv_disp_1st_;
        S32 * n_epj_recv_2nd_;
        S32 * n_epj_recv_disp_2nd_;
        S32 * id_proc_src_;
        S32 * id_proc_dest_;
        ReallocatableArray<S32> * id_ptcl_send_;
        ReallocatableArray<F64vec> * shift_image_box_;
        ReallocatableArray<S32> * ip_disp_;
        // new val
        ReallocatableArray< TreeParticle > * tp_scatter_;
        //ReallocatableArray< TreeCell< Tmomloc > > * tc_recv_1st_;
        ReallocatableArray<Tepj> * epj_recv_1st_sorted_;
        S32 ** adr_tc_level_partition_recv_1st_;

        //F64 periodic_length_[DIMENSTION];


        
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI::Request * req_send_;
        MPI::Request * req_recv_;
#endif

        template<class Tep2, class Tep3>
        inline void scatterEP(S32 n_send[],
                              S32 n_send_disp[],
                              S32 n_recv[],
                              S32 n_recv_disp[],
                              ReallocatableArray<Tep2> & ep_send,  // send buffer
                              ReallocatableArray<Tep2> & ep_recv,  // recv buffer
                              ReallocatableArray<Tep2> * ep_end_buf,  // send buffer
                              const ReallocatableArray<Tep3> & ep_org, // original
                              const DomainInfo & dinfo);

        template<class Tep2, class Tep3>
        inline void scatterEPForGather(S32 n_send[],
                                       S32 n_send_disp[],
                                       S32 n_recv[],
                                       S32 n_recv_disp[],
                                       ReallocatableArray<Tep2> & ep_send,  // send buffer
                                       ReallocatableArray<Tep2> & ep_recv,  // recv buffer
                                       const ReallocatableArray<Tep3> & ep_org, // original
                                       const DomainInfo & dinfo);

        void calcMomentLocalTreeOnlyImpl(TagExLetOneStage);
        void calcMomentLocalTreeOnlyImpl(TagExLetTwoStage);
        /*
        void calcMomentLocalTreeOnlyImpl(TagSearchLong);
        void calcMomentLocalTreeOnlyImpl(TagSearchLongCutoff);
        void calcMomentLocalTreeOnlyImpl(TagSearchLongScatter);
        void calcMomentLocalTreeOnlyImpl(TagSearchLongCutoffScatter);
        void calcMomentLocalTreeOnlyImpl(TagSearchShortScatter);
        void calcMomentLocalTreeOnlyImpl(TagSearchShortGather);
        void calcMomentLocalTreeOnlyImpl(TagSearchShortSymmetry);
        */

        void exchangeLocalEssentialTreeImpl(TagSearchLong, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchLongCutoff, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchLongScatter, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchLongCutoffScatter, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchShortScatter, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchShortGather, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeImpl(TagSearchShortSymmetry, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeGatherImpl(TagRSearch, const DomainInfo & dinfo);
        void exchangeLocalEssentialTreeGatherImpl(TagNoRSearch, const DomainInfo & dinfo);

        void setLocalEssentialTreeToGlobalTreeImpl(TagForceShort);
        void setLocalEssentialTreeToGlobalTreeImpl(TagForceLong);
        void setLocalEssentialTreeToGlobalTreeImpl3(TagForceShort);
        void setLocalEssentialTreeToGlobalTreeImpl3(TagForceLong);

        // probably, calcMomentGlobalTreeOnlyImpl is classified depending on force_type.
        void calcMomentGlobalTreeOnlyImpl(TagForceLong);
        void calcMomentGlobalTreeOnlyImpl(TagForceShort);
        //void calcMomentGlobalTreeOnlyImpl(TagSearchLong);
        //void calcMomentGlobalTreeOnlyImpl(TagSearchLongCutoff);
        //void calcMomentGlobalTreeOnlyImpl(TagSearchLongScatter);
        //void calcMomentGlobalTreeOnlyImpl(TagSearchShortScatter);
        //void calcMomentGlobalTreeOnlyImpl(TagSearchShortGather);
        //void calcMomentGlobalTreeOnlyImpl(TagSearchShortSymmetry);

        void addMomentAsSpLocalTreeImpl(TagForceLong);
        void addMomentAsSpLocalTreeImpl(TagForceShort);
        void addMomentAsSpGlobalTreeImpl(TagForceLong);
        void addMomentAsSpGlobalTreeImpl(TagForceShort);
        
        void addMomentAsSpLocalTreeLeanImpl(TagForceLong);
        void addMomentAsSpGlobalTreeLeanImpl(TagForceLong);
        
        void makeIPGroupImpl(TagForceLong);
        void makeIPGroupImpl(TagForceShort);

        void makeInteractionListImpl(TagSearchLong,          const S32 adr_ipg, const bool clear);
        void makeInteractionListImpl(TagSearchLongCutoff,    const S32 adr_ipg, const bool clear, const bool make_id_list=false);
        void makeInteractionListImpl(TagSearchLongScatter,   const S32 adr_ipg, const bool clear, const bool make_id_list=false);
        void makeInteractionListImpl(TagSearchShortScatter,  const S32 adr_ipg, const bool clear, const bool make_id_list=false);
        void makeInteractionListImpl(TagSearchShortGather,   const S32 adr_ipg, const bool clear, const bool make_id_list=false);
        void makeInteractionListImpl(TagSearchShortSymmetry, const S32 adr_ipg, const bool clear, const bool make_id_list=false);
        void makeInteractionListImpl2(TagSearchLong,         const S32 adr_ipg, const bool clear, const bool make_id_list=false);

        void makeInteractionListId(const S32 adr_ipg, const S32 ith, const bool clear=true);

        void makeAllInteractionListId(const bool clear=true);
        void makeAllInteractionListId3(const bool clear=true);        

        //void makeInteractionListIndexImpl(TagSearchLong, const S32 adr_ipg, const bool clear);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkImpl(TagForceLong,
                                   Tfunc_dispatch pfunc_dispatch,
                                   Tfunc_retrieve pfunc_retrieve,
                                   const S32 tag_max,
                                   const S32 n_walk_limit,
                                   const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkImpl(TagForceShort,
                                   Tfunc_dispatch pfunc_dispatch,
                                   Tfunc_retrieve pfunc_retrieve,
                                   const S32 tag_max,
                                   const S32 n_walk_limit,
                                   const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkIndexImpl(TagForceLong,
                                        Tfunc_dispatch pfunc_dispatch,
                                        Tfunc_retrieve pfunc_retrieve,
                                        const S32 tag_max,
                                        const S32 n_walk_limit,
                                        const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkIndexImpl(TagForceShort,
                                        Tfunc_dispatch pfunc_dispatch,
                                        Tfunc_retrieve pfunc_retrieve,
                                        const S32 tag_max,
                                        const S32 n_walk_limit,
                                        const bool clear=true);
	

        void checkMakeGlobalTreeImpl(TagForceLong, S32 & err, const F64vec & center, const F64 tolerance, std::ostream & fout);
        void checkMakeGlobalTreeImpl(TagForceShort, S32 & err, const F64vec & center, const F64 tolerance, std::ostream & fout);

        void checkCalcMomentLocalTreeImpl(TagSearchLong, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentLocalTreeImpl(TagSearchLongCutoff, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentLocalTreeImpl(TagSearchShortScatter, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentLocalTreeImpl(TagSearchShortGather, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentLocalTreeImpl(TagSearchShortSymmetry, const F64 tolerance, std::ostream & fout);
        
        void checkCalcMomentGlobalTreeImpl(TagSearchLong, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentGlobalTreeImpl(TagSearchLongCutoff, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentGlobalTreeImpl(TagSearchShortScatter, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentGlobalTreeImpl(TagSearchShortGather, const F64 tolerance, std::ostream & fout);
        void checkCalcMomentGlobalTreeImpl(TagSearchShortSymmetry, const F64 tolerance, std::ostream & fout);


        ///////////////////////
        ///// for open boundary
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchLongSymmetryOneStage){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epj_org_.getPointer());
        }
        // for P^3T
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchLongScatter){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epj_org_.getPointer());
        }
        // for P^3T
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchLongCutoffScatter){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epj_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchShortScatter){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epj_org_.getPointer());
        }
        /*
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchShortGather){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epi_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchShortSymmetry){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epi_org_.getPointer());
        }
        */
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchLongCutoff){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epj_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchLong){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epj_org_.getPointer());
        }
        template<class Tep2>
        void calcCenterAndLengthOfRootCellOpenNoMargenImpl(const Tep2 ep[]);
        template<class Tep2>
        void calcCenterAndLengthOfRootCellOpenWithMargenImpl(const Tep2 ep[]);

        /////////////
        //// PERIODIC
        // for P^3T
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchLongSymmetryOneStage){
            calcCenterAndLengthOfRootCellPeriodicImpl2(epj_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchLongCutoffScatter){
            calcCenterAndLengthOfRootCellPeriodicImpl2(epj_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchShortScatter){
            calcCenterAndLengthOfRootCellPeriodicImpl2(epj_org_.getPointer());
        }
        /*
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchShortGather){
            calcCenterAndLengthOfRootCellPeriodicImpl2(epi_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchShortSymmetry){
            calcCenterAndLengthOfRootCellPeriodicImpl2(epi_org_.getPointer());
        }
        */
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchLongCutoff){
            calcCenterAndLengthOfRootCellPeriodicImpl2(epj_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchLong){}
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchLongScatter){}

        template<class Tep2>
        void calcCenterAndLengthOfRootCellPeriodicImpl2(const Tep2 ep[]);

        void checkMortonSortGlobalTreeOnlyImpl(TagForceLong, std::ostream & fout);
        void checkMortonSortGlobalTreeOnlyImpl(TagForceShort, std::ostream & fout);

        void checkMakeInteractionListImpl(TagSearchLong,
                                          const DomainInfo & dinfo,
                                          const S32 adr_ipg,
                                          const S32 ith,
                                          const F64 tolerance,
                                          std::ostream & fout);
        void checkMakeInteractionListImpl(TagSearchLongCutoff,
					  const DomainInfo & dinfo,
					  const S32 adr_ipg,
					  const S32 ith,
					  const F64 tolerance,
					  std::ostream & fout);
        void checkMakeInteractionListImpl(TagSearchShortScatter,
                                          const DomainInfo & dinfo,
                                          const S32 adr_ipg,
                                          const S32 ith,
                                          const F64 tolerance,
                                          std::ostream & fout);
        void checkMakeInteractionListImpl(TagSearchShortGather,
                                          const DomainInfo & dinfo,
                                          const S32 adr_ipg,
                                          const S32 ith,
                                          const F64 tolerance,
                                          std::ostream & fout);
        void checkMakeInteractionListImpl(TagSearchShortSymmetry,
                                          const DomainInfo & dinfo,
                                          const S32 adr_ipg,
                                          const S32 ith,
                                          const F64 tolerance,
                                          std::ostream & fout);
        void checkExchangeLocalEssentialTreeImpl(TagForceLong,
                                                 const DomainInfo & dinfo,
                                                 const F64 tolerance,
                                                 std::ostream & fout);
        void checkExchangeLocalEssentialTreeForLongImpl(TagSearchLong,
                                                        const DomainInfo & dinfo,
                                                        const F64 tolerance, 
                                                        std::ostream & fout);
        void checkExchangeLocalEssentialTreeForLongImpl(TagSearchLongCutoff,
                                                        const DomainInfo & dinfo,
                                                        const F64 tolerance, 
                                                        std::ostream & fout);
        void checkExchangeLocalEssentialTreeImpl(TagForceShort,
                                                 const DomainInfo & dinfo,
                                                 const F64 tolerance,
                                                 std::ostream & fout);
        template<class Tep2>
        void checkExchangeLocalEssentialTreeForShortImpl
        (TagSearchShortScatter,
         const DomainInfo & dinfo,
         const Tep2 ep_tmp[],
         const S32 jp_head,
         const S32 jp_tail,
         const S32 rank_target,
         const S32 n_image_per_proc,
         const ReallocatableArray<F64vec> & shift_image_domain,
         S32 & n_recv_per_proc,
         ReallocatableArray<F64vec> & pos_direct);
        
        template<class Tep2>
        void checkExchangeLocalEssentialTreeForShortImpl
        (TagSearchShortGather,
         const DomainInfo & dinfo,
         const Tep2 ep_tmp[],
         const S32 jp_head,
         const S32 jp_tail,
         const S32 rank_target,
         const S32 n_image_per_proc,
         const ReallocatableArray<F64vec> & shift_image_domain,
         S32 & n_recv_per_proc,
         ReallocatableArray<F64vec> & pos_direct);
	
        template<class Tep2>
        void checkExchangeLocalEssentialTreeForShortImpl
        (TagSearchShortSymmetry,
         const DomainInfo & dinfo,
         const Tep2 ep_tmp[],
         const S32 jp_head,
         const S32 jp_tail,
         const S32 rank_target,
         const S32 n_image_per_proc,
         const ReallocatableArray<F64vec> & shift_image_domain,
         S32 & n_recv_per_proc,
         ReallocatableArray<F64vec> & pos_direct);
	
    public:
// new 
        void setPrefixOfProfile(const char * str){
            //profile.setPrefix(str);
        }

        // for neighbour search
        ReallocatableArray<Tepj> * epj_neighbor_;
        TimeProfile getTimeProfile() const {return time_profile_;}
        void clearTimeProfile(){time_profile_.clear();}
        CountT getNumberOfWalkLocal() const { return n_walk_local_; }
        CountT getNumberOfInteractionEPEPLocal() const { return n_interaction_ep_ep_local_; }
        CountT getNumberOfInteractionEPSPLocal() const { return n_interaction_ep_sp_local_; }
        CountT getNumberOfWalkGlobal() const { return Comm::getSum(n_walk_local_); }
        CountT getNumberOfInteractionEPEPGlobal() const { return Comm::getSum(n_interaction_ep_ep_local_); }
        CountT getNumberOfInteractionEPSPGlobal() const { return Comm::getSum(n_interaction_ep_sp_local_); }

        CountT getNumberOfLETEPSend1stLocal() const {return n_let_ep_send_1st_;}
        CountT getNumberOfLETEPRecv1stLocal() const {return n_let_ep_recv_1st_;}
        CountT getNumberOfLETSPSend1stLocal() const {return n_let_sp_send_1st_;}
        CountT getNumberOfLETSPRecv1stLocal() const {return n_let_sp_recv_1st_;}
        CountT getNumberOfLETEPSend2ndLocal() const {return n_let_ep_send_2nd_;}
        CountT getNumberOfLETEPRecv2ndLocal() const {return n_let_ep_recv_2nd_;}

        CountT getNumberOfLETEPSend1stGlobal() const {return Comm::getSum(n_let_ep_send_1st_);}
        CountT getNumberOfLETEPRecv1stGlobal() const {return Comm::getSum(n_let_ep_recv_1st_);}
        CountT getNumberOfLETSPSend1stGlobal() const {return Comm::getSum(n_let_sp_send_1st_);}
        CountT getNumberOfLETSPRecv1stGlobal() const {return Comm::getSum(n_let_sp_recv_1st_);}
        CountT getNumberOfLETEPSend2ndGlobal() const {return Comm::getSum(n_let_ep_send_2nd_);}
        CountT getNumberOfLETEPRecv2ndGlobal() const {return Comm::getSum(n_let_ep_recv_2nd_);}

        CountT getNumberOfCellOpenLocal() const {return n_cell_open_[0]; }
        CountT getNumberOfCellOpenGlobal() const {return Comm::getSum(n_cell_open_[0]); }
        CountT getNumberOfCellGlobal() const {return tc_glb_.size(); }

        CountT getNumberOfProcSendLET1stICommSP() const {return n_proc_send_exchange_LET_1st__icomm_sp_;}
        CountT getNumberOfProcRecvLET1stICommSP() const {return n_proc_recv_exchange_LET_1st__icomm_sp_;}
        CountT getNumberOfProcSendLET1stICommEP() const {return n_proc_send_exchange_LET_1st__icomm_ep_;}
        CountT getNumberOfProcRecvLET1stICommEP() const {return n_proc_recv_exchange_LET_1st__icomm_ep_;}

        CountT getNumberOfEPIAverageLocal() const {return n_epi_loc_ave_;}
        CountT getNumberOfEPJAverageLocal() const {return n_epj_loc_ave_;}
        CountT getNumberOfSPJAverageLocal() const {return n_spj_loc_ave_;}
        CountT getNumberOfEPIAverageGlobal() const {return Comm::getSum(n_epi_loc_ave_)/Comm::getNumberOfProc();}
        CountT getNumberOfEPJAverageGlobal() const {return Comm::getSum(n_epj_loc_ave_)/Comm::getNumberOfProc();}
        CountT getNumberOfSPJAverageGlobal() const {return Comm::getSum(n_spj_loc_ave_)/Comm::getNumberOfProc();}
        
        CountT getNumberOfIPG() const {return ipg_.size();}

        void clearNumberOfInteraction(){
            n_interaction_ep_ep_local_ = n_interaction_ep_sp_local_ = n_walk_local_ = 0;
        }
        void clearCounterAll(){
            n_let_ep_send_1st_ = n_let_ep_recv_1st_ = n_let_sp_send_1st_
                = n_let_sp_recv_1st_  = n_let_ep_send_2nd_ = n_let_ep_recv_2nd_ =  0;
            clearNumberOfInteraction();
            clearTimeProfile();
            //time_profile_.clear();
        }

        TreeForForce() : is_initialized_(false){}

        size_t getMemSizeUsed()const;
	
        void setNInteractionEPEP(const S64 n_ep_ep){
            n_interaction_ep_ep_ = n_ep_ep;
        }
        void setNInteractionEPSP(const S64 n_ep_sp){
            n_interaction_ep_sp_ = n_ep_sp;
        }

        S64 getNInteractionEPEP() const { return n_interaction_ep_ep_;}
        S64 getNInteractionEPSP() const { return n_interaction_ep_sp_;}


        void initialize(const U64 n_glb_tot,
                        const F64 theta=0.7,
                        const U32 n_leaf_limit=8,
                        const U32 n_group_limit=64);

        void reallocMem();
        void freeMem();
        void clearSizeOfArray();

        template<class Tpsys>
        void setParticleLocalTree(const Tpsys & psys, const bool clear=true);
        template<class Tpsys>
        void setParticleLocalTreeImpl(const Tpsys & psys,
                                      ReallocatableArray<Tepi> & epi,
                                      ReallocatableArray<Tepj> & epj,
                                      const bool clear);

#ifdef USE_MEMORY_POOL
        template<class Tpsys, class Tepjarray>
        void setParticleLocalTreeImpl(const Tpsys & psys,
                                      Tepjarray & epj,
                                      const bool clear);
#else
        template<class Tpsys>
        void setParticleLocalTreeImpl(const Tpsys & psys,
                                      ReallocatableArray<Tepj> & epj,
                                      const bool clear);        
#endif
        void setRootCell(const DomainInfo & dinfo);
        void setRootCell(const F64 l, const F64vec & c=F64vec(0.0));
        template<class Ttree>  void copyRootCell(const Ttree & tree);
        void mortonSortLocalTreeOnly();
        /*
        void mortonSortLocalTreeOnlyImpl(const ReallocatableArray<Tepi> & _epi_org,
                                         const ReallocatableArray<Tepj> & _epj_org,
                                         ReallocatableArray<Tepi> & _epi_sorted,
                                         ReallocatableArray<Tepj> & _epj_sorted);
        */
#ifdef USE_MEMORY_POOL
        template<class Tepjarray0, class Tepiarray0, class Tepjarray1>
        void mortonSortLocalTreeOnlyImpl(const Tepjarray0 & _epj_org,
                                         Tepiarray0 & _epi_sorted,
                                         Tepjarray1 & _epj_sorted);
#else
        void mortonSortLocalTreeOnlyImpl(const ReallocatableArray<Tepj> & _epj_org,
                                         ReallocatableArray<Tepi> & _epi_sorted,
                                         ReallocatableArray<Tepj> & _epj_sorted);
#endif
        
        void linkCellLocalTreeOnly();
        void linkCellGlobalTreeOnly();
        void calcMomentLocalTreeOnly();
        void calcMomentGlobalTreeOnly();
        //void addMomentAsSp(); // for multiwalk (send index)
        void makeIPGroup();
        void exchangeLocalEssentialTree(const DomainInfo & dinfo);
        void exchangeLocalEssentialTree2(const DomainInfo & dinfo, const bool reuse_list=false);
        void exchangeLocalEssentialTree3(const DomainInfo & dinfo, const bool reuse_list=false);
#ifdef USE_SUPER_DOMAIN
        void exchangeLocalEssentialTreeSuperDomain(const DomainInfo & dinfo, const bool reuse_list=false);
    #ifdef STATIC_DD
        void exchangeLocalEssentialTreeSuperDomain2(const DomainInfo & dinfo, const F64 r_search, const bool reuse_list=false);
    #else
        void exchangeLocalEssentialTreeSuperDomain2(const DomainInfo & dinfo, const bool reuse_list=false);
    #endif
#endif
        void setLocalEssentialTreeToGlobalTree();
        void setLocalEssentialTreeToGlobalTree2();
        void mortonSortGlobalTreeOnly();
        void mortonSortGlobalTreeOnly3();
        
        //S32 getNumberOfIPG() const { return ipg_.size();}
        void makeInteractionList(const S32 adr_ipg, const bool clear=true);
        void balanceNwalk(); // added by D.N.
        template<class Tfunc_ep_ep>
        void calcForceOnly(Tfunc_ep_ep pfunc_ep_ep,
                           const S32 adr_ipg,
                           const bool clear = true);
        template<class Tfunc_ep_ep, class Tfunc_ep_sp>
        void calcForceOnly(Tfunc_ep_ep pfunc_ep_ep,
                           Tfunc_ep_sp pfunc_ep_sp,
                           const S32 adr_ipg,
                           const bool clear = true);
        void copyForceOriginalOrder();
        template<class Tfunc_ep_ep>
        void calcForce(Tfunc_ep_ep pfunc_ep_ep,
                       const bool clear=true);
	
        template<class Tfunc_ep_ep>
        void calcForceWalkOnly(Tfunc_ep_ep pfunc_ep_ep,
                               const bool clear=true);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalk(Tfunc_dispatch pfunc_dispatch,
                               Tfunc_retrieve pfunc_retrieve,
                               const S32 tag_max,
                               const S32 n_walk_limit,
                               const bool clear=true);
	
        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                                    Tfunc_retrieve pfunc_retrieve,
                                    const S32 tag_max,
                                    const S32 n_walk_limit,
                                    const bool clear=true);

        template<class Tfunc_ep_ep, class Tfunc_ep_sp>
        void calcForce(Tfunc_ep_ep pfunc_ep_ep,
                       Tfunc_ep_sp pfunc_ep_sp,
                       const bool clear=true);

        template<class Tfunc_ep_ep, class Tfunc_ep_sp>
        void calcForce2(Tfunc_ep_ep pfunc_ep_ep,
                        Tfunc_ep_sp pfunc_ep_sp,
                        const bool clear,
                        const bool make_id_list=true);


        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                   Tpsys & psys,
                                   const bool clear=true){
            calcForce(pfunc_ep_ep, clear);
            const F64 time_offset = GetWtime();
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            time_profile_.calc_force += GetWtime() - time_offset;
        }

        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                   Tfunc_ep_sp pfunc_ep_sp,
                                   Tpsys & psys,
                                   const bool clear=true){
            calcForce(pfunc_ep_ep, pfunc_ep_sp, clear);
            const F64 time_offset = GetWtime();
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            time_profile_.calc_force += GetWtime() - time_offset;
        }



        Tforce getForce(const S32 i) const { return force_org_[i]; }

        ///////////////////////
        /// CHECK FUNCTIONS ///
        void checkMortonSortLocalTreeOnly(std::ostream & fout = std::cout);
        void checkMortonSortGlobalTreeOnly(std::ostream & fout = std::cout);
        void checkMakeLocalTree(const F64 tolerance = 1e-6, std::ostream & fout = std::cout);
        void checkMakeGlobalTree(const F64 tolerance = 1e-6, std::ostream & fout = std::cout);
        void checkCalcMomentLocalTree(const F64 tolerance = 1e-5, std::ostream & fout = std::cout);
        void checkCalcMomentGlobalTree(const F64 tolerance = 1e-5, std::ostream & fout = std::cout);
        void checkExchangeLocalEssentialTree(const DomainInfo & dinfo, 
                                             const F64 tolerance = 1e-5, 
                                             std::ostream & fout = std::cout);
        void checkMakeIPGroup(const F64 tolerance = 1e-5, std::ostream & fout = std::cout);
        void checkMakeInteractionList(const DomainInfo & dinfo,
                                      const S32 adr_ipg = 0, 
                                      const S32 ith = 0, 
                                      const F64 tolerance = 1e-5, 
                                      std::ostream & fout = std::cout){
            checkMakeInteractionListImpl(TSM::search_type(),  dinfo, adr_ipg, ith, tolerance, fout);
        }
        template<class Tfunc_ep_ep, class Tfunc_compare>
        void checkForce(Tfunc_ep_ep pfunc_ep_ep,
                        Tfunc_compare func_compare,
                        const DomainInfo & dinfo,
                        std::ostream & fout=std::cout);

        ////////////
        // for neighbour search APIs
        template<class Tptcl>
        S32 getNeighborListOneParticle(const Tptcl & ptcl, Tepj * & epj);
        template<class Tptcl>
        S32 getNeighborListOneParticleImpl(TagSearchShortScatter, const Tptcl & ptcl, Tepj * & epj);
        template<class Tptcl>
        S32 getNeighborListOneParticleImpl(TagSearchShortGather, const Tptcl & ptcl, Tepj * & epj);
        template<class Tptcl>
        S32 getNeighborListOneParticleImpl(TagSearchShortSymmetry, const Tptcl & ptcl, Tepj * & epj);
        template<class Tptcl>
        S32 getNeighborListOneParticleImpl(TagSearchLongScatter, const Tptcl & ptcl, Tepj * & epj);

        ////////////////////////////
        /// HIGH LEVEL FUNCTIONS ///
        //////////////////
        // FOR LONG FORCE
        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAll(Tfunc_ep_ep pfunc_ep_ep,
                          Tpsys & psys,
                          DomainInfo & dinfo,
                          const bool clear_force=true){
            setParticleLocalTree(psys);
            setRootCell(dinfo);
            mortonSortLocalTreeOnly();
            linkCellLocalTreeOnly();
            calcMomentLocalTreeOnly();
            exchangeLocalEssentialTree(dinfo);
            setLocalEssentialTreeToGlobalTree();
            mortonSortGlobalTreeOnly();
            linkCellGlobalTreeOnly();
            calcMomentGlobalTreeOnly();
            makeIPGroup();
            calcForce(pfunc_ep_ep, clear_force);
        }

        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllWalkOnly(Tfunc_ep_ep pfunc_ep_ep, 
                                  Tpsys & psys,
                                  DomainInfo & dinfo,
                                  const bool clear_force=true){
            setParticleLocalTree(psys);
            setRootCell(dinfo);
            mortonSortLocalTreeOnly();
            linkCellLocalTreeOnly();
            calcMomentLocalTreeOnly();
            exchangeLocalEssentialTree(dinfo);
            setLocalEssentialTreeToGlobalTree();
            mortonSortGlobalTreeOnly();
            linkCellGlobalTreeOnly();
            calcMomentGlobalTreeOnly();
            makeIPGroup();
            calcForceWalkOnly(pfunc_ep_ep, clear_force);
        }

        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllWithCheck(Tfunc_ep_ep pfunc_ep_ep, 
                                   Tpsys & psys,
                                   DomainInfo & dinfo,
                                   const bool clear_force=true){
            setParticleLocalTree(psys);
            setRootCell(dinfo);
            mortonSortLocalTreeOnly();
            checkMortonSortLocalTreeOnly(); // check morton sort
            linkCellLocalTreeOnly();
            checkMakeLocalTree(); // check link cell
            calcMomentLocalTreeOnly();
            checkCalcMomentLocalTree(); // check calc moment
            exchangeLocalEssentialTree(dinfo);
            checkExchangeLocalEssentialTree(dinfo); // check ex let
            setLocalEssentialTreeToGlobalTree();
            mortonSortGlobalTreeOnly();
            checkMortonSortGlobalTreeOnly(); // check morton sort 
            linkCellGlobalTreeOnly();
            checkMakeGlobalTree(); // check link cell
            calcMomentGlobalTreeOnly();
            checkCalcMomentGlobalTree(); // check calc moment 
            makeIPGroup();
            checkMakeIPGroup(); // check  make ipg
            calcForce(pfunc_ep_ep, clear_force);
        }

        template<class Tfunc_ep_ep>
        void calcForceMakingTree(Tfunc_ep_ep pfunc_ep_ep, 
                                 DomainInfo & dinfo,
                                 const bool clear_force=true){
            setRootCell(dinfo);
            mortonSortLocalTreeOnly();
            linkCellLocalTreeOnly();
            calcMomentLocalTreeOnly();
            exchangeLocalEssentialTree(dinfo);
            setLocalEssentialTreeToGlobalTree();
            mortonSortGlobalTreeOnly();
            linkCellGlobalTreeOnly();
            calcMomentGlobalTreeOnly();
            makeIPGroup();
            calcForce(pfunc_ep_ep, clear_force);
        }

        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                      Tpsys & psys,
                                      DomainInfo & dinfo,
                                      const bool clear_force = true){
            calcForceAll(pfunc_ep_ep, psys, dinfo, clear_force); 
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
        }


        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllWalkOnlyAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, 
                                              Tpsys & psys,
                                              DomainInfo & dinfo,
                                              const bool clear_force = true){
            calcForceAllWalkOnly(pfunc_ep_ep, psys, dinfo, clear_force); 
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
        }


        template<class Tfunc_ep_ep, class Tpsys>
        void calcForceAllAndWriteBackWithCheck(Tfunc_ep_ep pfunc_ep_ep, 
                                               Tpsys & psys,
                                               DomainInfo & dinfo,
                                               const bool clear_force = true){
            calcForceAllWithCheck(pfunc_ep_ep, psys, dinfo, clear_force);
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
        }



        //////////////////
        // FOR LONG FORCE
        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAll(Tfunc_ep_ep pfunc_ep_ep, 
                          Tfunc_ep_sp pfunc_ep_sp,  
                          Tpsys & psys,
                          DomainInfo & dinfo,
                          const bool clear_force=true){
            setParticleLocalTree(psys);
            setRootCell(dinfo);
            mortonSortLocalTreeOnly();
            linkCellLocalTreeOnly();
            calcMomentLocalTreeOnly();
            exchangeLocalEssentialTree(dinfo);
            setLocalEssentialTreeToGlobalTree();
            mortonSortGlobalTreeOnly();
            linkCellGlobalTreeOnly();
            calcMomentGlobalTreeOnly();
            makeIPGroup();
            calcForce(pfunc_ep_ep, pfunc_ep_sp, clear_force);
        }

        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAllWithCheck(Tfunc_ep_ep pfunc_ep_ep, 
                                   Tfunc_ep_sp pfunc_ep_sp,  
                                   Tpsys & psys,
                                   DomainInfo & dinfo,
                                   const bool clear_force,
                                   std::ostream & fout){

            fout<<"setParticleLocalTree"<<std::endl;
            setParticleLocalTree(psys);
            fout<<"setRootCell"<<std::endl;
            setRootCell(dinfo);
            fout<<"mortonSortLocalTreeOnly"<<std::endl;
            mortonSortLocalTreeOnly();

            checkMortonSortLocalTreeOnly(fout); // check morton sort
            fout<<"linkCellLocalTreeOnly"<<std::endl;
            linkCellLocalTreeOnly();
            checkMakeLocalTree(1e-6, fout); // check link cell
            fout<<"calcMomentLocalTreeOnly"<<std::endl;
            calcMomentLocalTreeOnly();
            checkCalcMomentLocalTree(1e-6, fout); // check calc moment
            fout<<"exchangeLocalEssentialTree"<<std::endl;
            exchangeLocalEssentialTree(dinfo);
            checkExchangeLocalEssentialTree(dinfo); // check ex let
            fout<<"setLocalEssentialTreeToGlobalTree"<<std::endl;
            setLocalEssentialTreeToGlobalTree();
            fout<<"mortonSortGlobalTreeOnly"<<std::endl;
            mortonSortGlobalTreeOnly();
            checkMortonSortGlobalTreeOnly(fout); // check morton sort 
            fout<<"linkCellGlobalTreeOnly"<<std::endl;
            linkCellGlobalTreeOnly();
            checkMakeGlobalTree(1e-6, fout); // check link cell
            fout<<"calcMomentGlobalTreeOnly"<<std::endl;
            calcMomentGlobalTreeOnly();
            checkCalcMomentGlobalTree(1e-6, fout); // check calc moment 
            fout<<"makeIPGroup"<<std::endl;
            makeIPGroup();
            checkMakeIPGroup(1e-6, fout); // check  make ipg
            fout<<"calcForce"<<std::endl;
            calcForce(pfunc_ep_ep, pfunc_ep_sp, clear_force);
        }

        template<class Tfunc_ep_ep, class Tfunc_ep_sp>
        void calcForceMakingTree(Tfunc_ep_ep pfunc_ep_ep, 
                                 Tfunc_ep_sp pfunc_ep_sp,  
                                 DomainInfo & dinfo,
                                 const bool clear_force=true){
            setRootCell(dinfo);
            mortonSortLocalTreeOnly();
            linkCellLocalTreeOnly();
            calcMomentLocalTreeOnly();
            exchangeLocalEssentialTree(dinfo);
            setLocalEssentialTreeToGlobalTree();
            mortonSortGlobalTreeOnly();
            linkCellGlobalTreeOnly();
            calcMomentGlobalTreeOnly();
            makeIPGroup();
            calcForce(pfunc_ep_ep, pfunc_ep_sp, clear_force);
        }

        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAllAndWriteBack(Tfunc_ep_ep pfunc_ep_ep, 
                                      Tfunc_ep_sp pfunc_ep_sp,  
                                      Tpsys & psys,
                                      DomainInfo & dinfo,
                                      const bool clear_force=true){

            clearSizeOfArray();
            calcForceAll(pfunc_ep_ep, pfunc_ep_sp, psys, dinfo, clear_force);
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
        }





        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAllAndWriteBackWithCheck(Tfunc_ep_ep pfunc_ep_ep,
                                               Tfunc_ep_sp pfunc_ep_sp,  
                                               Tpsys & psys,
                                               DomainInfo & dinfo,
                                               const bool clear_force,
                                               std::ostream & fout=std::cout){
            calcForceAllWithCheck(pfunc_ep_ep, pfunc_ep_sp, psys, dinfo, clear_force, fout);
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
        }


        template<class Tfunc_dispatch, class Tfunc_retrieve, class Tpsys>
        S32 calcForceAllMultiWalk(Tfunc_dispatch pfunc_dispatch,
                                  Tfunc_retrieve pfunc_retrieve,
                                  const S32 tag_max,
                                  Tpsys & psys,
                                  DomainInfo & dinfo,
                                  const S32 n_walk_limit,
                                  const bool clear=true){
            S32 ret = 0;
            setParticleLocalTree(psys);
            setRootCell(dinfo);
            mortonSortLocalTreeOnly();
            linkCellLocalTreeOnly();
            calcMomentLocalTreeOnly();
            exchangeLocalEssentialTree(dinfo);
            setLocalEssentialTreeToGlobalTree();
            mortonSortGlobalTreeOnly();
            linkCellGlobalTreeOnly();
            calcMomentGlobalTreeOnly();
            makeIPGroup();
            ret = calcForceMultiWalk(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, clear);
            return ret;
        }

        template<class Tfunc_dispatch, class Tfunc_retrieve, class Tpsys>
        S32 calcForceAllAndWriteBackMultiWalk(Tfunc_dispatch pfunc_dispatch,
                                              Tfunc_retrieve pfunc_retrieve,
                                              const S32 tag_max,
                                              Tpsys & psys,
                                              DomainInfo & dinfo,
                                              const S32 n_walk_limit,
                                              const bool clear=true){
            S32 ret = 0;
            ret = calcForceAllMultiWalk(pfunc_dispatch, pfunc_retrieve,
                                        tag_max, psys, dinfo, n_walk_limit, clear);
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            return ret;
        }


        ///////
        // new 
        template<class Tfunc_dispatch, class Tfunc_retrieve, class Tpsys>
        S32 calcForceAllMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                                       Tfunc_retrieve pfunc_retrieve,
                                       const S32 tag_max,
                                       Tpsys & psys,
                                       DomainInfo & dinfo,
                                       const S32 n_walk_limit,
                                       const bool clear=true){
            S32 ret = 0;
            setParticleLocalTree(psys);
            setRootCell(dinfo);
            mortonSortLocalTreeOnly();
            linkCellLocalTreeOnly();
            calcMomentLocalTreeOnly();
            exchangeLocalEssentialTree(dinfo);
            setLocalEssentialTreeToGlobalTree();
            mortonSortGlobalTreeOnly();
            linkCellGlobalTreeOnly();
            calcMomentGlobalTreeOnly();
            addMomentAsSpImpl(typename TSM::force_type());
            makeIPGroup();
            ret = calcForceMultiWalkIndex(pfunc_dispatch, pfunc_retrieve, tag_max, n_walk_limit, clear);
            return ret;
        }
        template<class Tfunc_dispatch, class Tfunc_retrieve, class Tpsys>
        S32 calcForceAllAndWriteBackMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                                                   Tfunc_retrieve pfunc_retrieve,
                                                   const S32 tag_max,
                                                   Tpsys & psys,
                                                   DomainInfo & dinfo,
                                                   const S32 n_walk_limit,
                                                   const bool clear=true){
            S32 ret = 0;
            ret = calcForceAllMultiWalkIndex(pfunc_dispatch, pfunc_retrieve,
                                             tag_max, psys, dinfo, n_walk_limit, clear);
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            return ret;
        }
	

        template<class Tfunc_ep_ep>
        void calcForceDirect(Tfunc_ep_ep pfunc_ep_ep,
                             Tforce force[],
                             const DomainInfo & dinfo,
                             const bool clear=true);

        template<class Tfunc_ep_ep>
        void calcForceDirectAndWriteBack(Tfunc_ep_ep pfunc_ep_ep,
                                         const DomainInfo & dinfo,
                                         const bool clear=true);

        void dump(std::ostream & fout){
            fout<<"n_loc_tot_="<<n_loc_tot_<<std::endl;
            fout<<"n_glb_tot_="<<n_glb_tot_<<std::endl;
            fout<<"length_="<<length_<<std::endl;
            fout<<"center_="<<center_<<std::endl;
            fout<<"pos_root_cell_="<<pos_root_cell_<<std::endl;
        }

        ////////////////////////////////
        // for reuse of interaction list
        void exchangeLocalEssentialTreeUsingCommTable();

        template<class Tfunc_ep_ep, class Tfunc_ep_sp>
        void calcForceUsingIdList(Tfunc_ep_ep pfunc_ep_ep, 
                                  Tfunc_ep_sp pfunc_ep_sp,
                                  const bool clear_force);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceUsingIdListMultiWalk(Tfunc_dispatch pfunc_dispatch,
                                           Tfunc_retrieve pfunc_retrieve,
                                           const S32 n_walk_limit,
                                           MultiWalkInfo<Tepi, Tepj, Tspj, Tforce> & mw_info,
                                           const bool clear_force);

        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceUsingIdListMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                                               Tfunc_retrieve pfunc_retrieve,
                                               const S32 n_walk_limit,
                                               MultiWalkInfo<Tepi, Tepj, Tspj, Tforce> & mw_info,
                                               const bool clear_force);


        template<class Tfunc_dispatch, class Tfunc_retrieve>
        S32 calcForceUsingIdListMultiWalkIndex3
        (Tfunc_dispatch pfunc_dispatch,
         Tfunc_retrieve pfunc_retrieve,
         const S32 n_walk_limit,
         MultiWalkInfo<Tepi, Tepj, Tspj, Tforce> & mw_info,
         const bool clear_force,
         const bool reuse_list = true,
         const bool flag_retrieve=true);

        template<class Tfp>
        void mortonSortFP(ParticleSystem<Tfp> & sys);
        
        template<class Tfunc_dispatch, class Tfunc_retrieve, class Tfp>
        S32 calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe
        (Tfunc_dispatch pfunc_dispatch,
         Tfunc_retrieve pfunc_retrieve,
         const S32 tag_max,
         ParticleSystem<Tfp> & psys,
         DomainInfo & dinfo,
         const S32 n_walk_limit,
         const bool clear_force,
         const bool reuse);


        template<class Tfunc_dispatch, class Tfunc_retrieve, class Tfp>
        S32 calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge
        (Tfunc_dispatch pfunc_dispatch,
         Tfunc_retrieve pfunc_retrieve,
         const S32 tag_max,
         ParticleSystem<Tfp> & psys,
         DomainInfo & dinfo,
         const S32 n_walk_limit,
         const bool clear_force,
         const bool reuse);

        void dumpMemSizeUsed(std::ostream & fout){
            const S32 thread_num = 1;
            const S32 cpe_num = 64;
            F64 norm_fac=1e-9;
            std::string str_unit=" [GB]";
            U64 size_of_cpe_buffers = 0;
            for (S32 cpe_id=0; cpe_id<cpe_num; cpe_id++)
                size_of_cpe_buffers += (adr_epj_loc_group_head_[cpe_id].getMemSize()
                                       +adr_epj_glb_group_head_[cpe_id].getMemSize()
                                       +group_size_epj_loc_[cpe_id].getMemSize()
                                       +group_size_epj_glb_[cpe_id].getMemSize());
            F64 sum = (double) (tp_buf_.getMemSize() 
                                //+ tp_loc_.getMemSize()
                              + tp_glb_.getMemSize()
                              + tc_loc_.getMemSize()
                              + tc_glb_.getMemSize()
                                //+ tc_recv_1st_[0].getMemSize()
                              + epi_sorted_.getMemSize()
                              + epj_sorted_.getMemSize()
                              + epj_org_.getMemSize()
                              + spj_sorted_.getMemSize()
                              + spj_org_.getMemSize()
                              + epj_sorted_loc_.getMemSize()
                              + spj_sorted_loc_.getMemSize()
                              + ipg_.getMemSize()
                              + epj_send_.getMemSize()
                              + epj_recv_.getMemSize()
                              + spj_send_.getMemSize()
                              + spj_recv_.getMemSize()
                              + id_ep_send_buf_[0].getMemSize()
                              + id_sp_send_buf_[0].getMemSize()
                              + rank_proc_send_[0].getMemSize()
                              + force_sorted_.getMemSize()
                              + force_org_.getMemSize()
                              + epj_for_force_[0].getMemSize()
                              + spj_for_force_[0].getMemSize()
                              + adr_epj_for_force_[0].getMemSize()
                              + adr_spj_for_force_[0].getMemSize()
                              + id_epj_recorder_for_force_[0].getMemSize()
                              + id_spj_recorder_for_force_[0].getMemSize()
                              + n_epi_recorder_for_force_[0].getMemSize()
                              + n_disp_epi_recorder_for_force_[0].getMemSize()
                              + n_epj_recorder_for_force_[0].getMemSize()
                              + n_spj_recorder_for_force_[0].getMemSize()
                              + n_disp_epj_recorder_for_force_[0].getMemSize()
                              + n_disp_spj_recorder_for_force_[0].getMemSize()
                              + adr_epj_loc2glb_.getMemSize()
                              + adr_epj_glb2loc_.getMemSize()
                              + size_of_cpe_buffers
                              + adr_epj_buf2glb_.getMemSize()
                              + adr_spj_buf2glb_.getMemSize()
                              + adr_epj_org2glb_.getMemSize()
                              + adr_ptcl_of_tp_loc_.getMemSize()
                              + id_recorder_for_interaction_.getMemSize()
                              + tree_outer_pos_.getMemSize()
                              + ep_send_buf_for_scatter_[0].getMemSize()
                              + shift_image_domain_[0].getMemSize()
                              + epjr_sorted_.getMemSize()
                              + epjr_send_.getMemSize()
                              + epjr_recv_.getMemSize()
                              + epjr_recv_1st_buf_.getMemSize()
                              + epjr_recv_2nd_buf_.getMemSize()
                              + epjr_send_buf_[0].getMemSize()
                              + epjr_send_buf_for_scatter_[0].getMemSize()
                              + epjr_recv_1st_sorted_[0].getMemSize()
                              + epj_recv_1st_buf_.getMemSize()
                              + epj_recv_2nd_buf_.getMemSize()
                              + epj_send_buf_[0].getMemSize()
                              + id_ptcl_send_[0].getMemSize()
                              + shift_image_box_[0].getMemSize()
                              + ip_disp_[0].getMemSize()
                              + tp_scatter_[0].getMemSize()
                              + epj_recv_1st_sorted_[0].getMemSize()
                              + etc_loc_.getMemSize()
                              + etc_glb_.getMemSize())
                              * norm_fac;
            F64 sum_max;
            MPI::COMM_WORLD.Allreduce(&sum,&sum_max,1,MPI::DOUBLE,MPI::MAX);

            if (Comm::getRank() == 0) {
                fout<<"*** TreeForForce class ***" << std::endl;
                fout<<"tp_buf_                           = " << tp_buf_.getMemSize()*norm_fac                           << str_unit << std::endl;
                fout<<"tp_buf_2_                         = " << tp_buf_2_.getMemSize()*norm_fac                           << str_unit << std::endl;
                //fout<<"tp_loc_                           = " << tp_loc_.getMemSize()*norm_fac                           << str_unit << std::endl;
                fout<<"tp_glb_                           = " << tp_glb_.getMemSize()*norm_fac                           << str_unit << std::endl;
                fout<<"tc_loc_                           = " << tc_loc_.getMemSize()*norm_fac                           << str_unit << std::endl; 
                fout<<"tc_glb_                           = " << tc_glb_.getMemSize()*norm_fac                           << str_unit << std::endl;
                //fout<<"tc_recv_1st_                      = " << tc_recv_1st_[0].getMemSize()*norm_fac                      << str_unit << std::endl;

                fout<<"epi_sorted_                       = " << epi_sorted_.getMemSize()*norm_fac                       << str_unit << std::endl;
                fout<<"epj_sorted_                       = " << epj_sorted_.getMemSize()*norm_fac                       << str_unit << std::endl;
                fout<<"epj_org_                          = " << epj_org_.getMemSize()*norm_fac                          << str_unit << std::endl;
                fout<<"spj_sorted_                       = " << spj_sorted_.getMemSize()*norm_fac                       << str_unit << std::endl;
                fout<<"spj_org_                          = " << spj_org_.getMemSize()*norm_fac                          << str_unit << std::endl;
                fout<<"epj_sorted_loc_                   = " << epj_sorted_loc_.getMemSize()*norm_fac                   << str_unit << std::endl;
                fout<<"spj_sorted_loc_                   = " << spj_sorted_loc_.getMemSize()*norm_fac                   << str_unit <<std::endl;
                fout<<"ipg_                              = " << ipg_.getMemSize()*norm_fac                              << str_unit <<std::endl;

                fout<<"epj_send_                         = " << epj_send_.getMemSize()*norm_fac                         << str_unit <<std::endl;
                fout<<"epj_recv_                         = " << epj_recv_.getMemSize()*norm_fac                         << str_unit <<std::endl;
                fout<<"spj_send_                         = " << spj_send_.getMemSize()*norm_fac                         << str_unit <<std::endl;
                fout<<"spj_recv_                         = " << spj_recv_.getMemSize()*norm_fac                         << str_unit <<std::endl;

                fout<<"id_ep_send_buf_[0]                = " << id_ep_send_buf_[0].getMemSize()*norm_fac                << str_unit <<std::endl;
                fout<<"id_sp_send_buf_[0]                = " << id_sp_send_buf_[0].getMemSize()*norm_fac                << str_unit <<std::endl;
                fout<<"rank_proc_send_[0]                = " << rank_proc_send_[0].getMemSize()*norm_fac                << str_unit <<std::endl;

                fout<<"force_sorted_                     = " << force_sorted_.getMemSize()*norm_fac                     << str_unit <<std::endl;
                fout<<"force_org_                        = " << force_org_.getMemSize()*norm_fac                        << str_unit << std::endl;
    
                fout<<"epj_for_force_[0]                 = " << epj_for_force_[0].getMemSize()*norm_fac                 << str_unit <<std::endl;
                fout<<"spj_for_force_[0]                 = " << spj_for_force_[0].getMemSize()*norm_fac                 << str_unit <<std::endl;

                fout<<"adr_epj_for_force_[0]             = " << adr_epj_for_force_[0].getMemSize()*norm_fac             << str_unit << std::endl;
                fout<<"adr_spj_for_force_[0]             = " << adr_spj_for_force_[0].getMemSize()*norm_fac             << str_unit << std::endl;

                fout<<"id_epj_recorder_for_force_[0]     = " << id_epj_recorder_for_force_[0].getMemSize()*norm_fac     << str_unit << std::endl;
                fout<<"id_spj_recorder_for_force_[0]     = " << id_spj_recorder_for_force_[0].getMemSize()*norm_fac     << str_unit << std::endl;
                fout<<"n_epi_recorder_for_force_[0]      = " << n_epi_recorder_for_force_[0].getMemSize()*norm_fac      << str_unit << std::endl;
                fout<<"n_disp_epi_recorder_for_force_[0] = " << n_disp_epi_recorder_for_force_[0].getMemSize()*norm_fac << str_unit << std::endl;
                fout<<"n_epj_recorder_for_force_[0]      = " << n_epj_recorder_for_force_[0].getMemSize()*norm_fac      << str_unit << std::endl;
                fout<<"n_spj_recorder_for_force_[0]      = " << n_spj_recorder_for_force_[0].getMemSize()*norm_fac      << str_unit << std::endl;
                fout<<"n_disp_epj_recorder_for_force_[0] = " << n_disp_epj_recorder_for_force_[0].getMemSize()*norm_fac << str_unit << std::endl;
                fout<<"n_disp_spj_recorder_for_force_[0] = " << n_disp_spj_recorder_for_force_[0].getMemSize()*norm_fac << str_unit << std::endl;

                fout<<"adr_epj_loc2glb_                  = " << adr_epj_loc2glb_.getMemSize()*norm_fac                  << str_unit << std::endl;

                fout<<"adr_epj_glb2loc_                  = " << adr_epj_glb2loc_.getMemSize()*norm_fac                  << str_unit << std::endl;
                fout<<"sum(adr_epj_loc_group_head_       = " << size_of_cpe_buffers*norm_fac                            << str_unit << std::endl;
                fout<<"    adr_epj_glb_group_head_         "                                                                        << std::endl;
                fout<<"    group_size_epj_loc_             "                                                                        << std::endl;
                fout<<"    group_size_epj_glb_)            "                                                                        << std::endl;

                fout<<"adr_epj_buf2glb_                  = " << adr_epj_buf2glb_.getMemSize()*norm_fac                  << str_unit << std::endl;
                fout<<"adr_spj_buf2glb_                  = " << adr_spj_buf2glb_.getMemSize()*norm_fac                  << str_unit << std::endl;
                fout<<"adr_epj_org2glb_                  = " << adr_epj_org2glb_.getMemSize()*norm_fac                  << str_unit << std::endl;
                fout<<"adr_ptcl_of_tp_loc_               = " << adr_ptcl_of_tp_loc_.getMemSize()*norm_fac               << str_unit << std::endl;
                fout<<"id_recorder_for_interaction_      = " << id_recorder_for_interaction_.getMemSize()*norm_fac      << str_unit << std::endl;
                
                fout<<"tree_outer_pos_                   = " << tree_outer_pos_.getMemSize()*norm_fac                   << str_unit << std::endl;
                fout<<"ep_send_buf_for_scatter_[0]       = " << ep_send_buf_for_scatter_[0].getMemSize()*norm_fac       << str_unit << std::endl;
                fout<<"shift_image_domain_[0]            = " << shift_image_domain_[0].getMemSize()*norm_fac            << str_unit << std::endl;

                fout<<"epjr_sorted_                      = " << epjr_sorted_.getMemSize()*norm_fac                      << str_unit << std::endl;
                fout<<"epjr_send_                        = " << epjr_send_.getMemSize()*norm_fac                        << str_unit << std::endl;
                fout<<"epjr_recv_                        = " << epjr_recv_.getMemSize()*norm_fac                        << str_unit << std::endl;
                fout<<"epjr_recv_1st_buf_                = " << epjr_recv_1st_buf_.getMemSize()*norm_fac                << str_unit << std::endl;
                fout<<"epjr_recv_2nd_buf_                = " << epjr_recv_2nd_buf_.getMemSize()*norm_fac                << str_unit << std::endl;
                fout<<"epjr_send_buf_[0]                 = " << epjr_send_buf_[0].getMemSize()*norm_fac                 << str_unit << std::endl;
                fout<<"epjr_send_buf_for_scatter_[0]     = " << epjr_send_buf_for_scatter_[0].getMemSize()*norm_fac     << str_unit << std::endl;
                fout<<"epjr_recv_1st_sorted_[0]          = " << epjr_recv_1st_sorted_[0].getMemSize()*norm_fac          << str_unit << std::endl;

                fout<<"epj_recv_1st_buf_                 = " << epj_recv_1st_buf_.getMemSize()*norm_fac                 << str_unit << std::endl;
                fout<<"epj_recv_2nd_buf_                 = " << epj_recv_2nd_buf_.getMemSize()*norm_fac                 << str_unit << std::endl;
                fout<<"epj_send_buf_[0]                  = " << epj_send_buf_[0].getMemSize()*norm_fac                  << str_unit << std::endl;
                fout<<"id_ptcl_send_[0]                  = " << id_ptcl_send_[0].getMemSize()*norm_fac                  << str_unit << std::endl;
                fout<<"shift_image_box_[0]               = " << shift_image_box_[0].getMemSize()*norm_fac               << str_unit << std::endl;
                fout<<"ip_disp_[0]                       = " << ip_disp_[0].getMemSize()*norm_fac                       << str_unit << std::endl;
                fout<<"tp_scatter_[0]                    = " << tp_scatter_[0].getMemSize()*norm_fac                    << str_unit << std::endl;
                fout<<"epj_recv_1st_sorted_[0]           = " << epj_recv_1st_sorted_[0].getMemSize()*norm_fac           << str_unit << std::endl;

                fout<<"etc_loc_                          = " << etc_loc_.getMemSize()*norm_fac                          << str_unit << std::endl;
                fout<<"etc_glb_                          = " << etc_glb_.getMemSize()*norm_fac                          << str_unit << std::endl;

                fout<<"sum                               = " << sum                                                     << str_unit << std::endl;
                fout<<"sum (tree,max.)                   = " << sum_max                                                 << str_unit << std::endl;
                fout<<"sum-epi_sorted_-epj_sorted_loc_   = " <<sum-(epi_sorted_.getMemSize()+epj_sorted_loc_.getMemSize())*norm_fac<< str_unit << std::endl;
                comm_table_.dumpMemSizeUsed(fout);
            }
        }

        
#ifdef REMOVE_VERTEX
        F64vec getCenterOfTree(){
            return center_;
        }
        F64 getHalfLengthOfTree(){
            return 0.5*length_;
        }
#endif
        
    };

    template<class Tforce, class Tepi, class Tepj, class Tmom=void, class Tsp=void>
    class TreeForForceLong{
    public:
        typedef TreeForForce
        <SEARCH_MODE_LONG,
         Tforce, Tepi, Tepj,
         Tmom, Tmom, Tsp> Normal;

        typedef TreeForForce
        <SEARCH_MODE_LONG_CUTOFF,
         Tforce, Tepi, Tepj,
         Tmom, Tmom, Tsp> WithCutoff;

        typedef TreeForForce
        <SEARCH_MODE_LONG_SCATTER,
         Tforce, Tepi, Tepj,
         Tmom, Tmom, Tsp> WithScatterSearch; // for P^3T

        typedef TreeForForce
        <SEARCH_MODE_LONG_CUTOFF_SCATTER,
         Tforce, Tepi, Tepj,
         Tmom, Tmom, Tsp> WithCutoffScatterSearch; // for P^3T
    };

    template<class Tforce, class Tepi, class Tepj>
    class TreeForForceLong<Tforce, Tepi, Tepj, void, void>{
    public:
        // for P^3T
        typedef TreeForForce
        <SEARCH_MODE_LONG_SCATTER,
         Tforce, Tepi, Tepj,
         MomentMonopoleScatter,
         MomentMonopoleScatter,
         SPJMonopoleScatter> MonopoleWithScatterSearch;

        typedef TreeForForce
        <SEARCH_MODE_LONG_SYMMETRY_ONE_STAGE,
         Tforce, Tepi, Tepj,
         MomentMonopoleInAndOut,
         MomentMonopoleInAndOut,
         SPJMonopoleInAndOut> MonopoleWithSymmetrySearch2;

        typedef TreeForForce
        <SEARCH_MODE_LONG_SCATTER,
         Tforce, Tepi, Tepj,
         MomentQuadrupoleScatter,
         MomentQuadrupoleScatter,
         SPJQuadrupoleScatter> QuadrupoleWithScatterSearch;

        // for P^3T + PM
        typedef TreeForForce
        <SEARCH_MODE_LONG_CUTOFF_SCATTER,
         Tforce, Tepi, Tepj,
         MomentMonopoleCutoffScatter,
         MomentMonopoleCutoffScatter,
         SPJMonopoleCutoffScatter> MonopoleWithCutoffScatterSearch;

        typedef TreeForForce
        <SEARCH_MODE_LONG,
         Tforce, Tepi, Tepj,
         MomentMonopole,
         MomentMonopole,
         SPJMonopole> Monopole;
	
        typedef TreeForForce
        <SEARCH_MODE_LONG_CUTOFF,
         Tforce, Tepi, Tepj,
         MomentMonopoleCutoff,
         MomentMonopoleCutoff,
         SPJMonopoleCutoff> MonopoleWithCutoff;
	
        typedef TreeForForce
        <SEARCH_MODE_LONG,
         Tforce, Tepi, Tepj,
         MomentQuadrupole,
         MomentQuadrupole,
         SPJQuadrupole> Quadrupole;

        typedef TreeForForce
        <SEARCH_MODE_LONG,
         Tforce, Tepi, Tepj,
         MomentMonopoleGeometricCenter,
         MomentMonopoleGeometricCenter,
         SPJMonopoleGeometricCenter> MonopoleGeometricCenter;

        typedef TreeForForce
        <SEARCH_MODE_LONG,
         Tforce, Tepi, Tepj,
         MomentDipoleGeometricCenter,
         MomentDipoleGeometricCenter,
         SPJDipoleGeometricCenter> DipoleGeometricCenter;

        typedef TreeForForce
        <SEARCH_MODE_LONG,
         Tforce, Tepi, Tepj,
         MomentQuadrupoleGeometricCenter,
         MomentQuadrupoleGeometricCenter,
         SPJQuadrupoleGeometricCenter> QuadrupoleGeometricCenter;
    };

    template<class Tforce, class Tepi, class Tepj>
    class TreeForForceShort{
    public:

        typedef TreeForForce 
        <SEARCH_MODE_SYMMETRY,
         Tforce, Tepi, Tepj,
         MomentSearchInAndOut,
         MomentSearchInAndOut,
         SuperParticleBase> Symmetry;

        typedef TreeForForce 
        <SEARCH_MODE_GATHER,
         Tforce, Tepi, Tepj,
         MomentSearchInAndOut,
         MomentSearchInOnly,
         SuperParticleBase> Gather;

        // send_tree: out
        // recv_tree: in
        // loc_tree: in
        // glb_tree: out
        typedef TreeForForce 
        <SEARCH_MODE_SCATTER,
         Tforce, Tepi, Tepj,
         MomentSearchInAndOut,
         MomentSearchInAndOut,
         SuperParticleBase> Scatter;
    };
}
#include"tree_for_force_impl.hpp"
#include"tree_for_force_impl2.hpp"
#include"tree_for_force_impl3.hpp"
#include"tree_for_force_impl_exchange_let.hpp"
#include"tree_for_force_check_impl.hpp"

