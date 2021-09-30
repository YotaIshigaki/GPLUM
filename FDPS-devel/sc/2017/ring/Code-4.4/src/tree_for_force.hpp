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
        S32 n_epi_ave_;
        S32 n_epj_ave_;
        S32 n_spj_ave_;
        
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

        ReallocatableArray< TreeParticle > tp_buf_, tp_loc_, tp_glb_;
        ReallocatableArray< TreeCell< Tmomloc > > tc_loc_;
        ReallocatableArray< TreeCell< Tmomglb > > tc_glb_;
        ReallocatableArray< Tepi > epi_sorted_, epi_org_;
        ReallocatableArray< Tepj > epj_sorted_, epj_org_;
        ReallocatableArray< Tspj > spj_sorted_, spj_org_;
        ReallocatableArray< Tepj > epj_sorted_loc_; // for reuse of interaction list
        ReallocatableArray< Tspj > spj_sorted_loc_; // for reuse of interaction list
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
        ReallocatableArray<S32> * n_epi_recorder_for_force_;
        ReallocatableArray<S32> * n_disp_epi_recorder_for_force_;
        ReallocatableArray<S32> * n_epj_recorder_for_force_;
        ReallocatableArray<S32> * n_spj_recorder_for_force_;
        ReallocatableArray<S32> * n_disp_epj_recorder_for_force_;
        ReallocatableArray<S32> * n_disp_spj_recorder_for_force_;
        ReallocatableArray<S32>   adr_epj_loc2glb_;
        ReallocatableArray<S32>   adr_epj_buf2glb_;
        ReallocatableArray<S32>   adr_spj_buf2glb_;
        ReallocatableArray<S32>   adr_epj_org2glb_;
        //ReallocatableArray<U32>   adr_tp_loc_org2sorted_;
        ReallocatableArray<U32>   adr_ptcl_of_tp_loc_;
        ReallocatableArray<IdRecorderForInteraction> id_recorder_for_interaction_;
        S32 adr_tc_level_partition_loc_[TREE_LEVEL_LIMIT+2];
        S32 adr_tc_level_partition_glb_[TREE_LEVEL_LIMIT+2];
        S32 lev_max_loc_;
        S32 lev_max_glb_;
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
        ReallocatableArray< TreeCell< Tmomloc > > * tc_recv_1st_;
        ReallocatableArray<Tepj> * epj_recv_1st_sorted_;
        S32 ** adr_tc_level_partition_recv_1st_;
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
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchShortGather){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epi_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellOpenImpl(TagSearchShortSymmetry){
            calcCenterAndLengthOfRootCellOpenNoMargenImpl(epi_org_.getPointer());
        }
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
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchShortGather){
            calcCenterAndLengthOfRootCellPeriodicImpl2(epi_org_.getPointer());
        }
        void calcCenterAndLengthOfRootCellPeriodicImpl(TagSearchShortSymmetry){
            calcCenterAndLengthOfRootCellPeriodicImpl2(epi_org_.getPointer());
        }
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

        CountT getNumberOfEPIAverage() const {return n_epi_ave_;}
        CountT getNumberOfEPJAverage() const {return n_epj_ave_;}
        CountT getNumberOfSPJAverage() const {return n_spj_ave_;}
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
        
        void setRootCell(const DomainInfo & dinfo);
        void setRootCell(const F64 l, const F64vec & c=F64vec(0.0));
        template<class Ttree>  void copyRootCell(const Ttree & tree);
        void mortonSortLocalTreeOnly();
        void mortonSortLocalTreeOnlyImpl(const ReallocatableArray<Tepi> & _epi_org,
                                         const ReallocatableArray<Tepj> & _epj_org,
                                         ReallocatableArray<Tepi> & _epi_sorted,
                                         ReallocatableArray<Tepj> & _epj_sorted);
        void linkCellLocalTreeOnly();
        void linkCellGlobalTreeOnly();
        void calcMomentLocalTreeOnly();
        void calcMomentGlobalTreeOnly();
        //void addMomentAsSp(); // for multiwalk (send index)
        void makeIPGroup();
        void exchangeLocalEssentialTree(const DomainInfo & dinfo);
        void exchangeLocalEssentialTree2(const DomainInfo & dinfo, const bool reuse_list=false);
        void exchangeLocalEssentialTree3(const DomainInfo & dinfo, const bool reuse_list=false);
        void setLocalEssentialTreeToGlobalTree();
        void setLocalEssentialTreeToGlobalTree2();
        void mortonSortGlobalTreeOnly();
        void mortonSortGlobalTreeOnly3();
        
        //S32 getNumberOfIPG() const { return ipg_.size();}
        void makeInteractionList(const S32 adr_ipg, const bool clear=true);
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



        template<class Tfunc_ep_ep, class Tfunc_ep_sp, class Tpsys>
        void calcForceAllAndWriteBackReuseList(Tfunc_ep_ep pfunc_ep_ep, 
                                               Tfunc_ep_sp pfunc_ep_sp,  
                                               Tpsys & psys,
                                               DomainInfo & dinfo,
                                               const bool clear_force=true,
                                               const bool reuse = false){
            const F64 time_offset = GetWtime();
            if(!reuse){
                clearSizeOfArray();
                setParticleLocalTree(psys);
                setRootCell(dinfo);
                mortonSortLocalTreeOnly();
                linkCellLocalTreeOnly();
                calcMomentLocalTreeOnly();
                addMomentAsSpLocalTreeImpl(typename TSM::force_type());
                epj_sorted_loc_.resizeNoInitialize(epj_sorted_.size());
                for(S32 i=0; i<epj_sorted_.size(); i++){
                    epj_sorted_loc_[i] = epj_sorted_[i];
                }
                exchangeLocalEssentialTree2(dinfo);
                setLocalEssentialTreeToGlobalTree2();
                mortonSortGlobalTreeOnly();
                linkCellGlobalTreeOnly();
                calcMomentGlobalTreeOnly();
                addMomentAsSpGlobalTreeImpl(typename TSM::force_type());
                makeIPGroup();
                makeAllInteractionListId(true);
                calcForceUsingIdList(pfunc_ep_ep, pfunc_ep_sp, clear_force); // not yet
            }
            else{
                for(S32 i=0; i<spj_org_.size(); i++) spj_org_[i].clear();
                for(S32 i=0; i<spj_sorted_.size(); i++) spj_sorted_[i].clear();
                for(S32 i=0; i<spj_sorted_loc_.size(); i++) spj_sorted_loc_[i].clear();
                for(S32 i=0; i<Comm::getNumberOfThread(); i++){
                    for(S32 j=0; j<spj_for_force_[i].size(); j++) spj_for_force_[i][j].clear();
                }
                setParticleLocalTree(psys);
                F64 wtime_offset = GetWtime();
                for(S32 i=0; i<n_loc_tot_; i++){
                    const S32 adr = tp_loc_[i].adr_ptcl_;
                    epi_sorted_[i] = epi_org_[adr];
                    epj_sorted_[i] = epj_org_[adr];
                }
                time_profile_.morton_sort_local_tree += GetWtime() - wtime_offset;
                calcMomentLocalTreeOnly();
                addMomentAsSpLocalTreeImpl(typename TSM::force_type());
                epj_sorted_loc_.resizeNoInitialize(epj_sorted_.size());
                for(S32 i=0; i<epj_sorted_.size(); i++){
                    epj_sorted_loc_[i] = epj_sorted_[i];
                }
                exchangeLocalEssentialTree2(dinfo, true);
                const S32 n_proc = Comm::getNumberOfProc();
                const S32 n_ep_add = this->comm_table_.n_disp_ep_recv_[n_proc];
                const S32 offset_ep = this->n_loc_tot_;
                wtime_offset = GetWtime();
                for(S32 i=0; i<n_ep_add; i++){
                    this->epj_org_[offset_ep+i] = this->comm_table_.ep_recv_[i];
                }
                time_profile_.set_LET_ep_to_global_tree += GetWtime() - wtime_offset;                
                const S32 n_sp_add = this->comm_table_.n_disp_sp_recv_[n_proc];
                const S32 offset_sp = this->n_loc_tot_ + n_ep_add;
                wtime_offset = GetWtime();                
                for(S32 i=0; i<n_sp_add; i++){
                    this->spj_org_[offset_sp+i] = this->comm_table_.sp_recv_[i];
                }
                time_profile_.set_LET_sp_to_global_tree += GetWtime() - wtime_offset;
                epj_sorted_.resizeNoInitialize( n_glb_tot_ );
                spj_sorted_.resizeNoInitialize( spj_org_.size() ); // it is needed because spj from global tree is push_back into spj_sorted.
                wtime_offset = GetWtime();                
                for(S32 i=0; i<n_glb_tot_; i++){
                    const U32 adr = tp_glb_[i].adr_ptcl_;
                    if( GetMSB(adr) == 0){
                        epj_sorted_[i] = epj_org_[adr];
                    }
                    else{
                        spj_sorted_[i] = spj_org_[ClearMSB(adr)];
                    }
                }
                time_profile_.morton_sort_global_tree += GetWtime() - wtime_offset;                
                calcMomentGlobalTreeOnly();
                addMomentAsSpGlobalTreeImpl(typename TSM::force_type());
                calcForceUsingIdList(pfunc_ep_ep, pfunc_ep_sp, clear_force);
            }
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            time_profile_.calc_force_all += GetWtime() - time_offset;
        }

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


        template<class Tfunc_dispatch, class Tfunc_retrieve, class Tpsys>
        S32 calcForceAllAndWriteBackReuseListMultiWalkImpl(Tfunc_dispatch pfunc_dispatch,
                                                           Tfunc_retrieve pfunc_retrieve,
                                                           const S32 tag_max,
                                                           Tpsys & psys,
                                                           DomainInfo & dinfo,
                                                           const S32 n_walk_limit,
                                                           const bool clear_force,
                                                           const bool reuse,
                                                           const int  mode=1){
            const F64 wtime_offset_0 = GetWtime();
            S32 ret = 0;
            if(!reuse){
                clearSizeOfArray();
                setParticleLocalTree(psys);
                setRootCell(dinfo);
                mortonSortLocalTreeOnly();
                linkCellLocalTreeOnly();
                calcMomentLocalTreeOnly();
                addMomentAsSpLocalTreeImpl(typename TSM::force_type());
                epj_sorted_loc_.resizeNoInitialize(epj_sorted_.size());
                for(S32 i=0; i<epj_sorted_.size(); i++){
                    epj_sorted_loc_[i] = epj_sorted_[i];
                }
                exchangeLocalEssentialTree2(dinfo);
                setLocalEssentialTreeToGlobalTree2();
                mortonSortGlobalTreeOnly();
                linkCellGlobalTreeOnly();
                calcMomentGlobalTreeOnly();
                addMomentAsSpGlobalTreeImpl(typename TSM::force_type());
                makeIPGroup();
                makeAllInteractionListId(true);
                //ret += calcForceUsingIdListMultiWalk(pfunc_dispatch, pfunc_retrieve, n_walk_limit, mw_info_, clear_force);                
                if(mode==0){
                    //ret += calcForceUsingIdListMultiWalk(pfunc_dispatch, pfunc_retrieve, n_walk_limit, mw_info_, clear_force); // only differenc
                }
                else if(mode==1){
                    // send id only
                    ret += calcForceUsingIdListMultiWalkIndex(pfunc_dispatch, pfunc_retrieve, n_walk_limit, mw_info_, clear_force); // only differenc                    
                }
            }
            else{
                /*
                // probably not needed
                for(S32 i=0; i<spj_org_.size(); i++) spj_org_[i].clear();
                for(S32 i=0; i<spj_sorted_.size(); i++) spj_sorted_[i].clear();
                for(S32 i=0; i<spj_sorted_loc_.size(); i++) spj_sorted_loc_[i].clear();
                for(S32 i=0; i<Comm::getNumberOfThread(); i++){
                    for(S32 j=0; j<spj_for_force_[i].size(); j++) spj_for_force_[i][j].clear();
                }
                */
                setParticleLocalTree(psys);
                epj_sorted_loc_.resizeNoInitialize(n_loc_tot_);
                F64 wtime_offset = GetWtime();
#ifdef MY_UNROLL
                epi_sorted_.reserveEmptyAreaAtLeast(8);
                epj_sorted_.reserveEmptyAreaAtLeast(8);
                epi_org_.reserveEmptyAreaAtLeast(8);
                epj_org_.reserveEmptyAreaAtLeast(8);
                epj_sorted_loc_.reserveEmptyAreaAtLeast(8);
                tp_loc_.reserveEmptyAreaAtLeast(8);
                //for(S32 i=0; i<n_loc_tot_/8+1; i++){
                for(S32 i=0; i<n_loc_tot_/8; i++){
                    const S32 adr0 = tp_loc_[i+0].adr_ptcl_;
                    const S32 adr1 = tp_loc_[i+1].adr_ptcl_;
                    const S32 adr2 = tp_loc_[i+2].adr_ptcl_;
                    const S32 adr3 = tp_loc_[i+3].adr_ptcl_;
                    const S32 adr4 = tp_loc_[i+4].adr_ptcl_;
                    const S32 adr5 = tp_loc_[i+5].adr_ptcl_;
                    const S32 adr6 = tp_loc_[i+6].adr_ptcl_;
                    const S32 adr7 = tp_loc_[i+7].adr_ptcl_;
                    epi_sorted_[i+0] = epi_org_[adr0];
                    epj_sorted_[i+0] = epj_sorted_loc_[i+0] = epj_org_[adr0];
                    epi_sorted_[i+1] = epi_org_[adr1];
                    epj_sorted_[i+1] = epj_sorted_loc_[i+1] = epj_org_[adr1];
                    epi_sorted_[i+2] = epi_org_[adr2];
                    epj_sorted_[i+2] = epj_sorted_loc_[i+2] = epj_org_[adr2];
                    epi_sorted_[i+3] = epi_org_[adr3];
                    epj_sorted_[i+3] = epj_sorted_loc_[i+3] = epj_org_[adr3];
                    epi_sorted_[i+4] = epi_org_[adr4];
                    epj_sorted_[i+4] = epj_sorted_loc_[i+4] = epj_org_[adr4];
                    epi_sorted_[i+5] = epi_org_[adr5];
                    epj_sorted_[i+5] = epj_sorted_loc_[i+5] = epj_org_[adr5];
                    epi_sorted_[i+6] = epi_org_[adr6];
                    epj_sorted_[i+6] = epj_sorted_loc_[i+6] = epj_org_[adr6];
                    epi_sorted_[i+7] = epi_org_[adr7];
                    epj_sorted_[i+7] = epj_sorted_loc_[i+7] = epj_org_[adr7];
                }
                for(S32 i=(n_loc_tot_/8)*8; i<n_loc_tot_; i++){
                    const S32 adr = tp_loc_[i].adr_ptcl_;
                    epi_sorted_[i] = epi_org_[adr];
                    epj_sorted_[i] = epj_sorted_loc_[i] = epj_org_[adr];
                }
#else
                for(S32 i=0; i<n_loc_tot_; i++){
                    const S32 adr = tp_loc_[i].adr_ptcl_;
                    epi_sorted_[i] = epi_org_[adr];
                    //epj_sorted_[i] = epj_org_[adr];
                    epj_sorted_[i] = epj_sorted_loc_[i] = epj_org_[adr];
                }
#endif
                time_profile_.morton_sort_local_tree += GetWtime() - wtime_offset;
                //std::cerr<<"time_profile_.morton_sort_local_tree= "<<time_profile_.morton_sort_local_tree<<std::endl;
                calcMomentLocalTreeOnly();
                addMomentAsSpLocalTreeImpl(typename TSM::force_type());
                /*
                epj_sorted_loc_.resizeNoInitialize(epj_sorted_.size());
                wtime_offset = GetWtime();   
                for(S32 i=0; i<epj_sorted_.size(); i++){
                    epj_sorted_loc_[i] = epj_sorted_[i];
                }
                time_profile_.morton_sort_local_tree += GetWtime() - wtime_offset;
                */
                exchangeLocalEssentialTree2(dinfo, true);
                const S32 n_proc = Comm::getNumberOfProc();
                const S32 n_ep_add = this->comm_table_.n_disp_ep_recv_[n_proc];
                const S32 offset_ep = this->n_loc_tot_;
                wtime_offset = GetWtime();
                for(S32 i=0; i<n_ep_add; i++){
                    this->epj_org_[offset_ep+i] = this->comm_table_.ep_recv_[i];
                }
                time_profile_.set_LET_ep_to_global_tree += GetWtime() - wtime_offset;
                const S32 n_sp_add = this->comm_table_.n_disp_sp_recv_[n_proc];
                const S32 offset_sp = this->n_loc_tot_ + n_ep_add;
                wtime_offset = GetWtime();                
                for(S32 i=0; i<n_sp_add; i++){
                    this->spj_org_[offset_sp+i] = this->comm_table_.sp_recv_[i];
                }
                time_profile_.set_LET_sp_to_global_tree += GetWtime() - wtime_offset;
                epj_sorted_.resizeNoInitialize( n_glb_tot_ );
                spj_sorted_.resizeNoInitialize( spj_org_.size() ); // it is needed because spj from global tree is push_back into spj_sorted.
                wtime_offset = GetWtime();
                //#ifdef MY_UNROLL
#ifdef MY_UNROLL
                //tp_glb_.reserveEmptyAreaAtLeast(8);
                //epj_sorted_.reserveEmptyAreaAtLeast(8);
                //epj_org_.reserveEmptyAreaAtLeast(8);                
                for(S32 i=0; i<n_glb_tot_/8; i++){
                    const U32 adr_0 = tp_glb_[i*8+0].adr_ptcl_;
                    const U32 adr_1 = tp_glb_[i*8+1].adr_ptcl_;
                    const U32 adr_2 = tp_glb_[i*8+2].adr_ptcl_;
                    const U32 adr_3 = tp_glb_[i*8+3].adr_ptcl_;
                    const U32 adr_4 = tp_glb_[i*8+4].adr_ptcl_;
                    const U32 adr_5 = tp_glb_[i*8+5].adr_ptcl_;
                    const U32 adr_6 = tp_glb_[i*8+6].adr_ptcl_;
                    const U32 adr_7 = tp_glb_[i*8+7].adr_ptcl_;
                    
                    if( GetMSB(adr_0) == 0) epj_sorted_[i*8+0] = epj_org_[adr_0];
                    else spj_sorted_[i*8+0] = spj_org_[ClearMSB(adr_0)];

                    if( GetMSB(adr_1) == 0) epj_sorted_[i*8+1] = epj_org_[adr_1];
                    else spj_sorted_[i*8+1] = spj_org_[ClearMSB(adr_1)];

                    if( GetMSB(adr_2) == 0) epj_sorted_[i*8+2] = epj_org_[adr_2];
                    else spj_sorted_[i*8+2] = spj_org_[ClearMSB(adr_2)];
                    
                    if( GetMSB(adr_3) == 0) epj_sorted_[i*8+3] = epj_org_[adr_3];
                    else spj_sorted_[i*8+3] = spj_org_[ClearMSB(adr_3)];

                    if( GetMSB(adr_4) == 0) epj_sorted_[i*8+4] = epj_org_[adr_4];
                    else spj_sorted_[i*8+4] = spj_org_[ClearMSB(adr_4)];

                    if( GetMSB(adr_5) == 0) epj_sorted_[i*8+5] = epj_org_[adr_5];
                    else spj_sorted_[i*8+5] = spj_org_[ClearMSB(adr_5)];

                    if( GetMSB(adr_6) == 0) epj_sorted_[i*8+6] = epj_org_[adr_6];
                    else spj_sorted_[i*8+6] = spj_org_[ClearMSB(adr_6)];

                    if( GetMSB(adr_7) == 0) epj_sorted_[i*8+7] = epj_org_[adr_7];
                    else spj_sorted_[i*8+7] = spj_org_[ClearMSB(adr_7)];
                }
                for(S32 i=(n_glb_tot_/8)*8; i<n_glb_tot_; i++){
                    const U32 adr = tp_glb_[i].adr_ptcl_;
                    if( GetMSB(adr) == 0) epj_sorted_[i] = epj_org_[adr];
                    else spj_sorted_[i] = spj_org_[ClearMSB(adr)];
                } 
#else
                for(S32 i=0; i<n_glb_tot_; i++){
                    const U32 adr = tp_glb_[i].adr_ptcl_;
                    if( GetMSB(adr) == 0){
                        epj_sorted_[i] = epj_org_[adr];
                    }
                    else{
                        spj_sorted_[i] = spj_org_[ClearMSB(adr)];
                    }
                }
#endif
                time_profile_.morton_sort_global_tree += GetWtime() - wtime_offset;
                calcMomentGlobalTreeOnly();
                addMomentAsSpGlobalTreeImpl(typename TSM::force_type());
                if(mode == 0){
                    //ret += calcForceUsingIdListMultiWalk(pfunc_dispatch, pfunc_retrieve, n_walk_limit, mw_info_, clear_force); // only differenc
                }
                else if(mode==1){
                    ret += calcForceUsingIdListMultiWalkIndex(pfunc_dispatch, pfunc_retrieve, n_walk_limit, mw_info_, clear_force); // only differenc                    
                }
            }
            F64 wtime_offset_1 = GetWtime();            
            for(S32 i=0; i<n_loc_tot_; i++) psys[i].copyFromForce(force_org_[i]);
            time_profile_.calc_force__copy_original_order += GetWtime() - wtime_offset_1;
            time_profile_.calc_force_all += GetWtime() - wtime_offset_0;
            return ret;
        }










        



        
        

        template<class Tfunc_dispatch, class Tfunc_retrieve, class Tpsys>
        S32 calcForceAllAndWriteBackReuseListMultiWalk(Tfunc_dispatch pfunc_dispatch,
                                                       Tfunc_retrieve pfunc_retrieve,
                                                       const S32 tag_max,
                                                       Tpsys & psys,
                                                       DomainInfo & dinfo,
                                                       const S32 n_walk_limit,
                                                       const bool clear_force=true,
                                                       const bool reuse = false){
            const S32 mode = 0;
            const S32 ret = calcForceAllAndWriteBackReuseListMultiWalkImpl<Tfunc_dispatch, Tfunc_retrieve, Tpsys>
                (pfunc_dispatch, pfunc_retrieve, tag_max, psys,
                 dinfo, n_walk_limit, clear_force, reuse,  mode);
            return ret;
        }

        template<class Tfunc_dispatch, class Tfunc_retrieve, class Tpsys>
        S32 calcForceAllAndWriteBackReuseListMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                                                            Tfunc_retrieve pfunc_retrieve,
                                                            const S32 tag_max,
                                                            Tpsys & psys,
                                                            DomainInfo & dinfo,
                                                            const S32 n_walk_limit,
                                                            const bool clear_force=true,
                                                            const bool reuse = false){
            const S32 mode = 1;
            const S32 ret = calcForceAllAndWriteBackReuseListMultiWalkImpl<Tfunc_dispatch, Tfunc_retrieve, Tpsys>
                (pfunc_dispatch, pfunc_retrieve, tag_max, psys,
                 dinfo, n_walk_limit, clear_force, reuse,  mode);
            return ret;
        }        



        void dumpMemSizeUsed(std::ostream & fout){
            //F64 sum_loc,sum; 
            //sum_loc = (double)
            //        (tp_buf_.getMemSize()+tp_loc_.getMemSize()+tp_glb_.getMemSize()
            //        +epi_sorted_.getMemSize()+epi_org_.getMemSize()+epj_sorted_.getMemSize()
            //        +epj_org_.getMemSize()+spj_sorted_.getMemSize()+spj_org_.getMemSize()
            //        +epj_sorted_loc_.getMemSize()+spj_sorted_loc_.getMemSize()+ipg_.getMemSize()
            //        +epj_send_.getMemSize()+epj_recv_.getMemSize()+spj_send_.getMemSize()
            //        +spj_recv_.getMemSize()+id_ep_send_buf_[0].getMemSize()+id_sp_send_buf_[0].getMemSize()
            //        +rank_proc_send_[0].getMemSize()+force_org_.getMemSize()+epj_for_force_[0].getMemSize()
            //        +spj_for_force_[0].getMemSize()) / 1e9;
            //if (Comm::getRank() == 0) {
            //    fout<<"tp_buf_.getMemSize()= "<<tp_buf_.getMemSize()<<std::endl;
            //    fout<<"tp_loc_.getMemSize()= "<<tp_loc_.getMemSize()<<std::endl;
            //    fout<<"tp_glb_.getMemSize()= "<<tp_glb_.getMemSize()<<std::endl;
            //    fout<<"epi_sorted_.getMemSize()= "<<epi_sorted_.getMemSize()<<std::endl;
            //    fout<<"epi_org_.getMemSize()= "<<epi_org_.getMemSize()<<std::endl;            
            //    fout<<"epj_sorted_.getMemSize()= "<<epj_sorted_.getMemSize()<<std::endl;
            //    fout<<"epj_org_.getMemSize()= "<<epj_org_.getMemSize()<<std::endl;
            //    fout<<"spj_sorted_.getMemSize()= "<<spj_sorted_.getMemSize()<<std::endl;
            //    fout<<"spj_org_.getMemSize()= "<<spj_org_.getMemSize()<<std::endl;
            //    fout<<"epj_sorted_loc_.getMemSize()= "<<epj_sorted_loc_.getMemSize()<<std::endl;
            //    fout<<"spj_sorted_loc_.getMemSize()= "<<spj_sorted_loc_.getMemSize()<<std::endl;
            //    fout<<"ipg_.getMemSize()= "<<ipg_.getMemSize()<<std::endl;
            //    fout<<"epj_send_.getMemSize()= "<<epj_send_.getMemSize()<<std::endl;
            //    fout<<"epj_recv_.getMemSize()= "<<epj_recv_.getMemSize()<<std::endl;
            //    fout<<"spj_send_.getMemSize()= "<<spj_send_.getMemSize()<<std::endl;
            //    fout<<"spj_recv_.getMemSize()= "<<spj_recv_.getMemSize()<<std::endl;
            //    fout<<"id_ep_send_buf_[0].getMemSize()= "<<id_ep_send_buf_[0].getMemSize()<<std::endl;
            //    fout<<"id_sp_send_buf_[0].getMemSize()= "<<id_sp_send_buf_[0].getMemSize()<<std::endl;
            //    fout<<"rank_proc_send_[0].getMemSize()= "<<rank_proc_send_[0].getMemSize()<<std::endl;
            //    fout<<"force_sorted_.getMemSize()= "<<force_org_.getMemSize()<<std::endl;
            //    fout<<"epj_for_force_.getMemSize()= "<<epj_for_force_[0].getMemSize()<<std::endl;
            //    fout<<"spj_for_force_.getMemSize()= "<<spj_for_force_[0].getMemSize()<<std::endl;
            //    fout<<"sum= "<< sum_loc << " [GB]" << std::endl;
            //}
            //MPI::COMM_WORLD.Allreduce(&sum_loc,&sum,1,MPI::DOUBLE,MPI::MAX);
            //if (Comm::getRank() == 0) {
            //    fout<<"sum(tree,max.) = " << sum << " [GB]" << std::endl;
            //}

            if (Comm::getRank() == 0) {
            fout<<"tp_buf_.getMemSize()= "<<tp_buf_.getMemSize()<<std::endl;
            fout<<"tp_loc_.getMemSize()= "<<tp_loc_.getMemSize()<<std::endl;
            fout<<"tp_glb_.getMemSize()= "<<tp_glb_.getMemSize()<<std::endl;
            fout<<"epi_sorted_.getMemSize()= "<<epi_sorted_.getMemSize()<<std::endl;
            fout<<"epi_org_.getMemSize()= "<<epi_org_.getMemSize()<<std::endl;            
            fout<<"epj_sorted_.getMemSize()= "<<epj_sorted_.getMemSize()<<std::endl;
            fout<<"epj_org_.getMemSize()= "<<epj_org_.getMemSize()<<std::endl;
            fout<<"spj_sorted_.getMemSize()= "<<spj_sorted_.getMemSize()<<std::endl;
            fout<<"spj_org_.getMemSize()= "<<spj_org_.getMemSize()<<std::endl;
            fout<<"epj_sorted_loc_.getMemSize()= "<<epj_sorted_loc_.getMemSize()<<std::endl;
            fout<<"spj_sorted_loc_.getMemSize()= "<<spj_sorted_loc_.getMemSize()<<std::endl;
            fout<<"ipg_.getMemSize()= "<<ipg_.getMemSize()<<std::endl;
            fout<<"epj_send_.getMemSize()= "<<epj_send_.getMemSize()<<std::endl;
            fout<<"epj_recv_.getMemSize()= "<<epj_recv_.getMemSize()<<std::endl;
            fout<<"spj_send_.getMemSize()= "<<spj_send_.getMemSize()<<std::endl;
            fout<<"spj_recv_.getMemSize()= "<<spj_recv_.getMemSize()<<std::endl;
            fout<<"id_ep_send_buf_[0].getMemSize()= "<<id_ep_send_buf_[0].getMemSize()<<std::endl;
            fout<<"id_sp_send_buf_[0].getMemSize()= "<<id_sp_send_buf_[0].getMemSize()<<std::endl;
            fout<<"rank_proc_send_[0].getMemSize()= "<<rank_proc_send_[0].getMemSize()<<std::endl;
            fout<<"force_sorted_.getMemSize()= "<<force_org_.getMemSize()<<std::endl;
            fout<<"epj_for_force_.getMemSize()= "<<epj_for_force_[0].getMemSize()<<std::endl;
            fout<<"spj_for_force_.getMemSize()= "<<spj_for_force_[0].getMemSize()<<std::endl;


            fout<<"id_epj_recorder_for_force_[0].getMemSize()= "<<id_epj_recorder_for_force_[0].getMemSize()<<std::endl;
            fout<<"id_spj_recorder_for_force_[0].getMemSize()= "<<id_spj_recorder_for_force_[0].getMemSize()<<std::endl;
            fout<<"n_epi_recorder_for_force_[0].getMemSize()= "<<n_epi_recorder_for_force_[0].getMemSize()<<std::endl;
            fout<<"n_disp_epi_recorder_for_force_[0].getMemSize()= "<<n_disp_epi_recorder_for_force_[0].getMemSize()<<std::endl;
            fout<<"n_epj_recorder_for_force_[0].getMemSize()= "<<n_epj_recorder_for_force_[0].getMemSize()<<std::endl;
            fout<<"n_spj_recorder_for_force_[0].getMemSize()= "<<n_spj_recorder_for_force_[0].getMemSize()<<std::endl;
            fout<<"n_disp_epj_recorder_for_force_[0].getMemSize()= "<<n_disp_epj_recorder_for_force_[0].getMemSize()<<std::endl;
            fout<<"n_disp_spj_recorder_for_force_[0].getMemSize()= "<<n_disp_spj_recorder_for_force_[0].getMemSize()<<std::endl;
            fout<<"adr_epj_loc2glb_.getMemSize()= "<<adr_epj_loc2glb_.getMemSize()<<std::endl;
            fout<<"adr_epj_buf2glb_.getMemSize()= "<<adr_epj_buf2glb_.getMemSize()<<std::endl;
            fout<<"adr_spj_buf2glb_.getMemSize()= "<<adr_spj_buf2glb_.getMemSize()<<std::endl;;
            fout<<"adr_epj_org2glb_.getMemSize()= "<<adr_epj_org2glb_.getMemSize()<<std::endl;
            fout<<"adr_ptcl_of_tp_loc_.getMemSize()= "<<adr_ptcl_of_tp_loc_.getMemSize()<<std::endl;
            fout<<"id_recorder_for_interaction_.getMemSize()= "<<id_recorder_for_interaction_.getMemSize()<<std::endl;
            
            fout<<"sum= "<<
                (double)
                (tp_buf_.getMemSize()+tp_loc_.getMemSize()+tp_glb_.getMemSize()
                +epi_sorted_.getMemSize()+epi_org_.getMemSize()+epj_sorted_.getMemSize()
                +epj_org_.getMemSize()+spj_sorted_.getMemSize()+spj_org_.getMemSize()
                +epj_sorted_loc_.getMemSize()+spj_sorted_loc_.getMemSize()+ipg_.getMemSize()
                +epj_send_.getMemSize()+epj_recv_.getMemSize()+spj_send_.getMemSize()
                +spj_recv_.getMemSize()+id_ep_send_buf_[0].getMemSize()+id_sp_send_buf_[0].getMemSize()
                +rank_proc_send_[0].getMemSize()+force_org_.getMemSize()+epj_for_force_[0].getMemSize()
                 +spj_for_force_[0].getMemSize()
                 +id_epj_recorder_for_force_[0].getMemSize()
                 +id_spj_recorder_for_force_[0].getMemSize()
                 +n_epi_recorder_for_force_[0].getMemSize()
                 +n_disp_epi_recorder_for_force_[0].getMemSize()
                 +n_epj_recorder_for_force_[0].getMemSize()
                 +n_spj_recorder_for_force_[0].getMemSize()
                 +n_disp_epj_recorder_for_force_[0].getMemSize()
                 +n_disp_spj_recorder_for_force_[0].getMemSize()
                 +adr_epj_loc2glb_.getMemSize()
                 +adr_epj_buf2glb_.getMemSize()
                 +adr_spj_buf2glb_.getMemSize()
                 +adr_epj_org2glb_.getMemSize()
                 +adr_ptcl_of_tp_loc_.getMemSize()
                 +id_recorder_for_interaction_.getMemSize()) / 1e9
                <<" [GB]"
                <<std::endl;
            comm_table_.dumpMemSizeUsed(fout);
            }
        }

        

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

