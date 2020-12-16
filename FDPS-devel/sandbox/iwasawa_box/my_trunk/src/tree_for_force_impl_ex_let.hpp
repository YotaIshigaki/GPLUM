#pragma once

#include"tree_walk.hpp"

namespace ParticleSimulator{
    


    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeLong(const DomainInfo & dinfo,
                                   const bool flag_reuse){
        F64 wtime_offset;
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
        if(!flag_reuse){
            wtime_offset = GetWtime();
            //comm_table_.clearSize();
            FindScatterParticle<TSM, TreeCell<Tmomloc>, TreeParticle,
                                Tepj, Tspj, TagWalkModeNormal>
                (tc_loc_, tp_loc_,
                 epj_sorted_, 
                 comm_table_.n_ep_send_,   comm_table_.adr_ep_send_, 
                 dinfo,          n_leaf_limit_,
                 comm_table_.n_sp_send_,   comm_table_.adr_sp_send_, 
                 comm_table_.shift_per_image_,
                 comm_table_.n_image_per_proc_,
                 comm_table_.n_ep_per_image_,
                 comm_table_.n_sp_per_image_,
                 r_crit_sq);
            time_profile_.make_LET_1st += GetWtime() - wtime_offset;
            wtime_offset = GetWtime();
            ExchangeNumber(comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                           comm_table_.n_sp_send_, comm_table_.n_sp_recv_);
            time_profile_.exchange_LET_1st += GetWtime() - wtime_offset;
        }
        wtime_offset = GetWtime();
        ExchangeLet<TSM, Tepj, Tspj>(epj_sorted_, comm_table_.n_ep_send_,
                                     comm_table_.n_ep_recv_, comm_table_.n_ep_per_image_,
                                     comm_table_.adr_ep_send_, epj_recv_,
                                     tc_loc_, comm_table_.n_sp_send_,
                                     comm_table_.n_sp_recv_, comm_table_.n_sp_per_image_,
                                     comm_table_.adr_sp_send_, spj_recv_,
                                     comm_table_.shift_per_image_,
                                     comm_table_.n_image_per_proc_);
        time_profile_.exchange_LET_1st += GetWtime() - wtime_offset;
        SumOfArray(comm_table_.n_ep_send_.getPointer(), comm_table_.n_ep_send_.size(), n_let_ep_send_1st_);
        SumOfArray(comm_table_.n_sp_send_.getPointer(), comm_table_.n_sp_send_.size(), n_let_sp_send_1st_);
        SumOfArray(comm_table_.n_ep_recv_.getPointer(), comm_table_.n_ep_recv_.size(), n_let_ep_recv_1st_);
        SumOfArray(comm_table_.n_sp_recv_.getPointer(), comm_table_.n_sp_recv_.size(), n_let_sp_recv_1st_);
        #if 0
        F64 mass_tmp = 0.0;
        for(S32 i=0; i<epi_sorted_.size(); i++) mass_tmp += epi_sorted_[i].mass;
        for(S32 i=0; i<epj_recv_.size(); i++) mass_tmp += epj_recv_[i].mass;
        for(S32 i=0; i<spj_recv_.size(); i++) mass_tmp += spj_recv_[i].mass;
        std::cerr<<"exLET) mass_tmp= "<<mass_tmp<<std::endl;
        #endif
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTree(const DomainInfo & dinfo,
                               const bool flag_reuse){
        if(typeid(TSM) == typeid(SEARCH_MODE_LONG)
           && dinfo.getBoundaryCondition() != BOUNDARY_CONDITION_OPEN){
            PARTICLE_SIMULATOR_PRINT_ERROR("The forces w/o cutoff can be evaluated only under the open boundary condition");
            Abort(-1);
        }
        if(!flag_reuse){ comm_table_.clearSize(); }
        //comm_table_.clear();
        exchangeLocalEssentialTreeImpl(typename TSM::search_type(), dinfo,flag_reuse);
        /*
        time_profile_.exchange_LET_tot = time_profile_.make_LET_1st
            + time_profile_.exchange_LET_1st
            + time_profile_.make_LET_2nd
            + time_profile_.exchange_LET_2nd;
        */
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchLong,
                                   const DomainInfo & dinfo,
                                   const bool flag_reuse){
        //F64 time_offset = GetWtime();
        exchangeLocalEssentialTreeLong(dinfo, flag_reuse);
	//time_profile_.exchange_LET_1st += GetWtime() - time_offset;
    }    
    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchLongCutoff,
                                   const DomainInfo & dinfo,
                                   const bool flag_reuse){
        exchangeLocalEssentialTreeLong(dinfo, flag_reuse);
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchLongScatter,
                                   const DomainInfo & dinfo,
                                   const bool flag_reuse){
        exchangeLocalEssentialTreeLong(dinfo, flag_reuse);
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
    class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchLongSymmetry,
                                   const DomainInfo & dinfo,
                                   const bool flag_reuse){
        exchangeLocalEssentialTreeLong(dinfo, flag_reuse);
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchShortScatter,
                                   const DomainInfo & dinfo,
                                   const bool flag_reuse){
        F64 wtime_offset;
        if(!flag_reuse){
            wtime_offset = GetWtime();
            FindScatterParticle<TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj,
                                TagWalkModeNormal>
                (tc_loc_, tp_loc_, epj_sorted_,
                 comm_table_.n_ep_send_,  comm_table_.adr_ep_send_,
                 dinfo,          n_leaf_limit_,
                 comm_table_.shift_per_image_,
                 comm_table_.n_image_per_proc_,
                 comm_table_.n_ep_per_image_);
            time_profile_.make_LET_1st += GetWtime() - wtime_offset;
            wtime_offset = GetWtime();
            ExchangeNumber(comm_table_.n_ep_send_, comm_table_.n_ep_recv_);
            time_profile_.exchange_LET_1st += GetWtime() - wtime_offset;
        }
        wtime_offset = GetWtime();
        ExchangeLet(epj_sorted_,
                    comm_table_.n_ep_send_,
                    comm_table_.n_ep_recv_,
                    comm_table_.n_ep_per_image_,
                    comm_table_.adr_ep_send_, epj_recv_,
                    comm_table_.shift_per_image_,
                    comm_table_.n_image_per_proc_);
        time_profile_.exchange_LET_1st += GetWtime() - wtime_offset;
        SumOfArray(comm_table_.n_ep_send_.getPointer(), comm_table_.n_ep_send_.size(), n_let_ep_send_1st_);
        SumOfArray(comm_table_.n_ep_recv_.getPointer(), comm_table_.n_ep_recv_.size(), n_let_ep_recv_1st_);
    }


    ////////////////
    // SYMMETRY
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchShortSymmetry,
                                   const DomainInfo & dinfo,
                                   const bool flag_reuse){
        const S32 n_proc = Comm::getNumberOfProc();
        static ReallocatableArray<S32> n_ep_send_per_proc_1st;
        static ReallocatableArray<S32> n_ep_send_per_proc_2nd;
        static ReallocatableArray<S32> n_ep_recv_per_proc_1st;
        static ReallocatableArray<S32> n_ep_recv_per_proc_2nd;
        static ReallocatableArray<S32> adr_ep_send_1st;
        static ReallocatableArray<S32> adr_ep_send_2nd;
        static ReallocatableArray<F64vec> shift_per_image_1st;
        static ReallocatableArray<F64vec> shift_per_image_2nd;
        static ReallocatableArray<S32> n_image_per_proc_1st;
        static ReallocatableArray<S32> n_image_per_proc_2nd;
        static ReallocatableArray<S32> n_ep_send_per_image_1st;
        static ReallocatableArray<S32> n_ep_send_per_image_2nd;
        static ReallocatableArray<Tepj> ep_recv_1st;
        
        if(!flag_reuse){
            ////////////
            // 1st STEP (send j particles)
            FindScatterParticle<TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj,
                                TagWalkModeNormal>
                (tc_loc_, tp_loc_, epj_sorted_,
                 n_ep_send_per_proc_1st,  adr_ep_send_1st,
                 dinfo,          n_leaf_limit_,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 n_ep_send_per_image_1st);
            ExchangeNumber(n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st);
            ExchangeParticle(epj_sorted_, n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st,
                             n_ep_send_per_image_1st,
                             adr_ep_send_1st, ep_recv_1st,
                             shift_per_image_1st,
                             n_image_per_proc_1st);
            ////////////
            // 2nd STEP (send j particles)
            FindExchangeParticleDoubleWalk<TreeCell<Tmomloc>, TreeParticle, Tepj>
                (ep_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo,
                 n_leaf_limit_,
                 n_ep_send_per_proc_2nd, n_ep_send_per_image_2nd,
                 n_image_per_proc_2nd,
                 adr_ep_send_2nd, shift_per_image_2nd,
                 epj_sorted_,
                 center_, length_);
            /////////////////////
            // 3rd STEP (exchange # of particles again)
            ExchangeNumber(n_ep_send_per_proc_2nd, n_ep_recv_per_proc_2nd);

            ReallocatableArray<S32> n_disp_ep_send_per_proc;
            n_disp_ep_send_per_proc.resizeNoInitialize(n_proc+1);
            n_disp_ep_send_per_proc[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_ep_send_per_proc[i+1] = n_disp_ep_send_per_proc[i]
                    + n_ep_send_per_proc_1st[i]
                    + n_ep_send_per_proc_2nd[i];
            }

            ReallocatableArray<S32> n_disp_image_per_proc;
            ReallocatableArray<S32> n_disp_image_per_proc_1st;
            ReallocatableArray<S32> n_disp_image_per_proc_2nd;
            n_disp_image_per_proc.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc_1st.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc_2nd.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc[0] = n_disp_image_per_proc_1st[0] = n_disp_image_per_proc_2nd[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_image_per_proc[i+1]     = n_disp_image_per_proc[i]     + n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
                n_disp_image_per_proc_1st[i+1] = n_disp_image_per_proc_1st[i] + n_image_per_proc_1st[i];
                n_disp_image_per_proc_2nd[i+1] = n_disp_image_per_proc_2nd[i] + n_image_per_proc_2nd[i];
            }

            const S32 n_image_tot_1st = shift_per_image_1st.size();
            const S32 n_image_tot_2nd = shift_per_image_2nd.size();
            ReallocatableArray<S32> n_disp_ep_send_per_image_1st;
            ReallocatableArray<S32> n_disp_ep_send_per_image_2nd;
            n_disp_ep_send_per_image_1st.resizeNoInitialize(n_image_tot_1st+1);
            n_disp_ep_send_per_image_2nd.resizeNoInitialize(n_image_tot_2nd+1);
            n_disp_ep_send_per_image_1st[0] = n_disp_ep_send_per_image_2nd[0] = 0;
            for(S32 i=0; i<n_image_tot_1st; i++){
                n_disp_ep_send_per_image_1st[i+1] = n_disp_ep_send_per_image_1st[i] + n_ep_send_per_image_1st[i];
            }
            for(S32 i=0; i<n_image_tot_2nd; i++){
                n_disp_ep_send_per_image_2nd[i+1] = n_disp_ep_send_per_image_2nd[i] + n_ep_send_per_image_2nd[i];
            }
            
            const S32 n_image_tot = n_disp_image_per_proc_1st[n_proc]  + n_disp_image_per_proc_2nd[n_proc];
            comm_table_.shift_per_image_.resizeNoInitialize(n_image_tot);
            comm_table_.n_ep_per_image_.resizeNoInitialize(n_image_tot);
            const S32 n_send_tot = n_disp_ep_send_per_proc[n_proc];
            comm_table_.adr_ep_send_.resizeNoInitialize(n_send_tot);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
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
                    comm_table_.shift_per_image_[image_head+n_image_cnt] = shift_per_image_1st[j];
                    comm_table_.n_ep_per_image_[image_head+n_image_cnt]  = n_ep_send_per_image_1st[j];
                    for(S32 k=ep_head_1st; k<ep_end_1st; k++, n_ep_cnt++){
                        comm_table_.adr_ep_send_[ep_head+n_ep_cnt] = adr_ep_send_1st[k];
                    }
                }
                const S32 image_head_2nd = n_disp_image_per_proc_2nd[i];
                const S32 image_end_2nd  = n_disp_image_per_proc_2nd[i+1];
                for(S32 j=image_head_2nd; j<image_end_2nd; j++, n_image_cnt++){
                    const S32 ep_head_2nd = n_disp_ep_send_per_image_2nd[j];
                    const S32 ep_end_2nd  = n_disp_ep_send_per_image_2nd[j+1];
                    comm_table_.shift_per_image_[image_head+n_image_cnt] = shift_per_image_2nd[j];
                    comm_table_.n_ep_per_image_[image_head+n_image_cnt]  = n_ep_send_per_image_2nd[j];
                    for(S32 k=ep_head_2nd; k<ep_end_2nd; k++, n_ep_cnt++){
                        comm_table_.adr_ep_send_[ep_head+n_ep_cnt] = adr_ep_send_2nd[k];
                    }
                }
            }
            comm_table_.n_ep_send_.resizeNoInitialize(n_proc);
            comm_table_.n_ep_recv_.resizeNoInitialize(n_proc);
            comm_table_.n_image_per_proc_.resizeNoInitialize(n_proc);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_proc; i++){
                comm_table_.n_ep_send_[i] = n_ep_send_per_proc_1st[i] + n_ep_send_per_proc_2nd[i];
                comm_table_.n_ep_recv_[i] = n_ep_recv_per_proc_1st[i] + n_ep_recv_per_proc_2nd[i];
                comm_table_.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_image_tot; i++){
                comm_table_.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
            }
        } // end of reuse
        ExchangeLet(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                    comm_table_.n_ep_per_image_,
                    comm_table_.adr_ep_send_, epj_recv_,
                    comm_table_.shift_per_image_, comm_table_.n_image_per_proc_);
    }


    ////////////////
    // GATHER MODE
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTreeImpl(TagSearchShortGather,
                                   const DomainInfo & dinfo,
                                   const bool flag_reuse){
        const S32 n_proc = Comm::getNumberOfProc();
        static ReallocatableArray<S32> n_ep_send_per_proc_1st;
        static ReallocatableArray<S32> n_ep_send_per_proc_2nd;
        static ReallocatableArray<S32> n_ep_recv_per_proc_1st;
        static ReallocatableArray<S32> n_ep_recv_per_proc_2nd;
        static ReallocatableArray<S32> adr_ep_send_1st;
        static ReallocatableArray<S32> adr_ep_send_2nd;
        static ReallocatableArray<F64vec> shift_per_image_1st;
        static ReallocatableArray<F64vec> shift_per_image_2nd;
        static ReallocatableArray<S32> n_image_per_proc_1st;
        static ReallocatableArray<S32> n_image_per_proc_2nd;
        static ReallocatableArray<S32> n_ep_send_per_image_1st;
        static ReallocatableArray<S32> n_ep_send_per_image_2nd;
        //static ReallocatableArray<Tepi> ep_recv_1st;
        static ReallocatableArray<EssentialParticleBase> ep_recv_1st;
        static ReallocatableArray<EssentialParticleBase> epi_base_sorted;
        const S32 n_epi_sorted = epi_sorted_.size();
        epi_base_sorted.resizeNoInitialize(n_epi_sorted);
        for(S32 i=0; i<n_epi_sorted; i++){
            epi_base_sorted[i].pos      = epi_sorted_[i].getPos();
            epi_base_sorted[i].r_search = epi_sorted_[i].getRSearch();
        }
        
        if(!flag_reuse){
            ////////////
            // 1st STEP (send epi_base)
            FindScatterParticle<TreeCell<Tmomloc>, TreeParticle, EssentialParticleBase, Tspj, TagWalkModeNormal>
                (tc_loc_, tp_loc_, epi_base_sorted,
                 n_ep_send_per_proc_1st,  adr_ep_send_1st,
                 dinfo,          n_leaf_limit_,
                 shift_per_image_1st,
                 n_image_per_proc_1st,
                 n_ep_send_per_image_1st);
            ExchangeNumber(n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st);
            ExchangeParticle(epi_base_sorted, n_ep_send_per_proc_1st, n_ep_recv_per_proc_1st,
                             n_ep_send_per_image_1st,
                             adr_ep_send_1st, ep_recv_1st,
                             shift_per_image_1st,
                             n_image_per_proc_1st);
            ////////////
            // 2nd STEP (find j particle)
            FindExchangeParticleDoubleWalk<TreeCell<Tmomloc>, TreeParticle, EssentialParticleBase>
                (ep_recv_1st, tc_loc_, n_ep_recv_per_proc_1st, n_image_per_proc_1st, dinfo,
                 n_leaf_limit_,
                 n_ep_send_per_proc_2nd, n_ep_send_per_image_2nd,
                 n_image_per_proc_2nd,
                 adr_ep_send_2nd, shift_per_image_2nd,
                 epi_base_sorted,
                 center_, length_);
            
            /////////////////////
            // 3rd STEP (exchange # of particles again)
            ExchangeNumber(n_ep_send_per_proc_2nd, n_ep_recv_per_proc_2nd);

            ReallocatableArray<S32> n_disp_ep_send_per_proc;
            n_disp_ep_send_per_proc.resizeNoInitialize(n_proc+1);
            n_disp_ep_send_per_proc[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_ep_send_per_proc[i+1] = n_disp_ep_send_per_proc[i]
                    + n_ep_send_per_proc_1st[i]
                    + n_ep_send_per_proc_2nd[i];
            }

            ReallocatableArray<S32> n_disp_image_per_proc;
            ReallocatableArray<S32> n_disp_image_per_proc_1st;
            ReallocatableArray<S32> n_disp_image_per_proc_2nd;
            n_disp_image_per_proc.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc_1st.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc_2nd.resizeNoInitialize(n_proc+1);
            n_disp_image_per_proc[0] = n_disp_image_per_proc_1st[0] = n_disp_image_per_proc_2nd[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_image_per_proc[i+1]     = n_disp_image_per_proc[i]     + n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
                n_disp_image_per_proc_1st[i+1] = n_disp_image_per_proc_1st[i] + n_image_per_proc_1st[i];
                n_disp_image_per_proc_2nd[i+1] = n_disp_image_per_proc_2nd[i] + n_image_per_proc_2nd[i];
            }

            const S32 n_image_tot_1st = shift_per_image_1st.size();
            const S32 n_image_tot_2nd = shift_per_image_2nd.size();
            ReallocatableArray<S32> n_disp_ep_send_per_image_1st;
            ReallocatableArray<S32> n_disp_ep_send_per_image_2nd;
            n_disp_ep_send_per_image_1st.resizeNoInitialize(n_image_tot_1st+1);
            n_disp_ep_send_per_image_2nd.resizeNoInitialize(n_image_tot_2nd+1);
            n_disp_ep_send_per_image_1st[0] = n_disp_ep_send_per_image_2nd[0] = 0;
            for(S32 i=0; i<n_image_tot_1st; i++){
                n_disp_ep_send_per_image_1st[i+1] = n_disp_ep_send_per_image_1st[i] + n_ep_send_per_image_1st[i];
            }
            for(S32 i=0; i<n_image_tot_2nd; i++){
                n_disp_ep_send_per_image_2nd[i+1] = n_disp_ep_send_per_image_2nd[i] + n_ep_send_per_image_2nd[i];
            }
            
            const S32 n_image_tot = n_disp_image_per_proc_1st[n_proc]  + n_disp_image_per_proc_2nd[n_proc];
            comm_table_.shift_per_image_.resizeNoInitialize(n_image_tot);
            comm_table_.n_ep_per_image_.resizeNoInitialize(n_image_tot);
            const S32 n_send_tot = n_disp_ep_send_per_proc[n_proc];
            comm_table_.adr_ep_send_.resizeNoInitialize(n_send_tot);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
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
                    comm_table_.shift_per_image_[image_head+n_image_cnt] = shift_per_image_1st[j];
                    comm_table_.n_ep_per_image_[image_head+n_image_cnt]  = n_ep_send_per_image_1st[j];
                    for(S32 k=ep_head_1st; k<ep_end_1st; k++, n_ep_cnt++){
                        comm_table_.adr_ep_send_[ep_head+n_ep_cnt] = adr_ep_send_1st[k];
                    }
                }
                const S32 image_head_2nd = n_disp_image_per_proc_2nd[i];
                const S32 image_end_2nd  = n_disp_image_per_proc_2nd[i+1];
                for(S32 j=image_head_2nd; j<image_end_2nd; j++, n_image_cnt++){
                    const S32 ep_head_2nd = n_disp_ep_send_per_image_2nd[j];
                    const S32 ep_end_2nd  = n_disp_ep_send_per_image_2nd[j+1];
                    comm_table_.shift_per_image_[image_head+n_image_cnt] = shift_per_image_2nd[j];
                    comm_table_.n_ep_per_image_[image_head+n_image_cnt]  = n_ep_send_per_image_2nd[j];
                    for(S32 k=ep_head_2nd; k<ep_end_2nd; k++, n_ep_cnt++){
                        comm_table_.adr_ep_send_[ep_head+n_ep_cnt] = adr_ep_send_2nd[k];
                    }
                }
            }
            comm_table_.n_ep_send_.resizeNoInitialize(n_proc);
            comm_table_.n_ep_recv_.resizeNoInitialize(n_proc);
            comm_table_.n_image_per_proc_.resizeNoInitialize(n_proc);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_proc; i++){
                comm_table_.n_ep_send_[i] = n_ep_send_per_proc_1st[i] + n_ep_send_per_proc_2nd[i];
                comm_table_.n_ep_recv_[i] = n_ep_recv_per_proc_1st[i] + n_ep_recv_per_proc_2nd[i];
                comm_table_.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
            for(S32 i=0; i<n_image_tot; i++){
                comm_table_.n_image_per_proc_[i] = n_image_per_proc_1st[i] + n_image_per_proc_2nd[i];
            }
        } // end of reuse flag
        ExchangeLet(epj_sorted_, comm_table_.n_ep_send_, comm_table_.n_ep_recv_,
                    comm_table_.n_ep_per_image_,
                    comm_table_.adr_ep_send_, epj_recv_,
                    comm_table_.shift_per_image_, comm_table_.n_image_per_proc_);
    }

    
}
