#include<multi_walk.hpp>

namespace ParticleSimulator{
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setLocalEssentialTreeToGlobalTree2(){
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 offset = this->n_loc_tot_;
        const S32 n_ep_add = this->comm_table_.n_disp_ep_recv_[n_proc];
        S32 n_sp_add = 0;
        if( (FORCE_TYPE)(typename TSM::force_type()).force_type == FORCE_TYPE_LONG){
            n_sp_add = this->comm_table_.n_disp_sp_recv_[n_proc];
        }
        const S32 offset2 = this->n_loc_tot_ + n_ep_add;
        this->n_glb_tot_ = this->n_loc_tot_ + n_ep_add + n_sp_add;
        this->tp_glb_.resizeNoInitialize( this->n_glb_tot_ );
        this->epj_org_.resizeNoInitialize( offset2 );
        if( (FORCE_TYPE)(typename TSM::force_type()).force_type == FORCE_TYPE_LONG){
            this->spj_org_.resizeNoInitialize( this->n_glb_tot_ );
        }
        else{
            this->spj_org_.resizeNoInitialize( 0 );
        }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	
#pragma omp parallel 
#endif	
        {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp for
#endif	    
            for(S32 i=0; i<offset; i++){
                this->tp_glb_[i] = this->tp_loc_[i]; // NOTE: need to keep tp_loc_[]?
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp for
#endif	    
            for(S32 i=0; i<n_ep_add; i++){
                this->epj_org_[offset+i] = this->comm_table_.ep_recv_[i];
                this->tp_glb_[offset+i].setFromEP(this->comm_table_.ep_recv_[i], offset+i);
            }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp for
#endif	    
            for(S32 i=0; i<n_sp_add; i++){
                this->spj_org_[offset2+i] = this->comm_table_.sp_recv_[i];
                this->tp_glb_[offset2+i].setFromSP(this->comm_table_.sp_recv_[i], offset2+i);
            }
        }
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep, class Tfunc_ep_sp>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceUsingIdList(Tfunc_ep_ep pfunc_ep_ep,
                         Tfunc_ep_sp pfunc_ep_sp,
                         const bool clear_force){
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        force_org_.resizeNoInitialize(n_loc_tot_);
        const S64 n_ipg = ipg_.size();
        if(n_ipg > 0){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for schedule(dynamic, 4) 
#endif	    
            for(S32 ig=0; ig<n_ipg; ig++){
                F64 mass = 0.0;
                const S32 ith     = id_recorder_for_interaction_[ig].ith_;
                const S32 n_epj   = id_recorder_for_interaction_[ig].n_epj_;
                const S32 adr_epj = id_recorder_for_interaction_[ig].adr_epj_;
                epj_for_force_[ith].resizeNoInitialize(n_epj);
                for(S32 jp=0; jp<n_epj; jp++){
                    const S32 adr_j = id_epj_recorder_for_force_[ith][adr_epj+jp];
                    epj_for_force_[ith][jp] = epj_sorted_[adr_j];
                    mass += epj_for_force_[ith][jp].getCharge();
                }
                const S32 n_spj   = id_recorder_for_interaction_[ig].n_spj_;
                const S32 adr_spj = id_recorder_for_interaction_[ig].adr_spj_;
                spj_for_force_[ith].resizeNoInitialize(n_spj);
                for(S32 jp=0; jp<n_spj; jp++){
                    const S32 adr_j = id_spj_recorder_for_force_[ith][adr_spj+jp];
                    spj_for_force_[ith][jp] = spj_sorted_[adr_j];
                    mass += spj_for_force_[ith][jp].getCharge();
                }
                const S32 n_epi   = ipg_[ig].n_ptcl_;
                const S32 adr_epi = ipg_[ig].adr_ptcl_;
                if(clear_force){
                    for(S32 ip=0; ip<n_epi; ip++){
                        force_sorted_[adr_epi+ip].clear();
                    }
                }
                pfunc_ep_ep(epi_sorted_.getPointer(adr_epi),    n_epi,
                            epj_for_force_[ith].getPointer(),   n_epj,
                            force_sorted_.getPointer(adr_epi));
                pfunc_ep_sp(epi_sorted_.getPointer(adr_epi),    n_epi,
                            spj_for_force_[ith].getPointer(),   n_spj,
                            force_sorted_.getPointer(adr_epi));
            }
        }
        copyForceOriginalOrder();
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceUsingIdListMultiWalk(Tfunc_dispatch pfunc_dispatch,
                                  Tfunc_retrieve pfunc_retrieve,
                                  const S32 n_walk_limit,
                                  MultiWalkInfo<Tepi, Tepj, Tspj, Tforce> & mw_info,
                                  const bool clear_force){
        std::vector< std::pair<S32, std::pair<S32, S32> > > iw2ith_adr;
        iw2ith_adr.reserve(n_walk_limit);
        S32 ret = 0;
        S32 tag = 0;
        mw_info.initialize(n_walk_limit);
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        force_org_.resizeNoInitialize(n_loc_tot_);
        if(clear_force){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for
#endif
            for(S32 ip=0; ip<n_loc_tot_; ip++){
                force_org_[ip].clear();
                force_sorted_[ip].clear();
            }
        }
        S32 n_walk_prev;
        const S64 n_ipg = ipg_.size();
        const S32 n_walk_grp = mw_info.getNWalkGrp(n_ipg);
        F64 pot0 = 0.0; // for debug
        //std::cerr<<"n_ipg= "<<n_ipg<<" n_walk_limit= "<<n_walk_limit<<" n_walk_grp= "<<n_walk_grp<<std::endl;
        for(S32 iwg=0; iwg<n_walk_grp; iwg++){
            const S32 n_walk = n_ipg/n_walk_grp + ( (iwg < n_ipg%n_walk_grp) ? 1 : 0 );
            const S32 walk_grp_head = (n_ipg/n_walk_grp)*iwg + ( (iwg < n_ipg%n_walk_grp) ? iwg : n_ipg%n_walk_grp);
            //std::cerr<<"n_walk= "<<n_walk<<" walk_grp_head= "<<walk_grp_head<<std::endl;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel
#endif
            {
                const S32 my_thread_id = Comm::getThreadNum();
                epj_for_force_[my_thread_id].clearSize();
                spj_for_force_[my_thread_id].clearSize();
                S32 cnt_epj = 0;
                S32 cnt_spj = 0;
                //F64 mass0 = 0.0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp for schedule(dynamic, 4) 
#endif	    
                for(S32 iw=0; iw<n_walk; iw++){
                    //mass0 = 0.0;
                    iw2ith_adr[iw].first = my_thread_id;
                    const S32 id_ipg = walk_grp_head + iw;
                    // copy epi
                    mw_info.n_epi_[iw] = ipg_[id_ipg].n_ptcl_;
                    const S32 adr_epi  = ipg_[id_ipg].adr_ptcl_;
                    mw_info.epi_[iw]   = epi_sorted_.getPointer(adr_epi);
                    mw_info.force_[iw] = force_sorted_.getPointer(adr_epi);
                    // copy epj
                    const S32 ith      = id_recorder_for_interaction_[id_ipg].ith_;
                    mw_info.n_epj_[iw] = id_recorder_for_interaction_[id_ipg].n_epj_;
                    const S32 adr_epj  = id_recorder_for_interaction_[id_ipg].adr_epj_;
                    //std::cerr<<"adr_epj= "<<adr_epj<<std::endl;
                    iw2ith_adr[iw].second.first = cnt_epj;
                    cnt_epj += mw_info.n_epj_[iw];
                    for(S32 jp=0; jp<mw_info.n_epj_[iw]; jp++){
                        const S32 adr_j = id_epj_recorder_for_force_[ith][adr_epj+jp];
                        epj_for_force_[my_thread_id].push_back( epj_sorted_[adr_j] );
                        //std::cout<<"epj_sorted_[adr_j].pos= "<<epj_sorted_[adr_j].pos<<std::endl;
                        //mass0 += epj_sorted_[adr_j].mass;
                    }
                    // copy spj
                    mw_info.n_spj_[iw] = id_recorder_for_interaction_[id_ipg].n_spj_;
                    const S32 adr_spj  = id_recorder_for_interaction_[id_ipg].adr_spj_;
                    iw2ith_adr[iw].second.second = cnt_spj;
                    cnt_spj += mw_info.n_spj_[iw];
                    for(S32 jp=0; jp<mw_info.n_spj_[iw]; jp++){
                        const S32 adr_j = id_spj_recorder_for_force_[ith][adr_spj+jp];
                        spj_for_force_[my_thread_id].push_back( spj_sorted_[adr_j] );
                        //mass0 += spj_sorted_[adr_j].mass;
                    }
                }
                //std::cout<<"mass0= "<<mass0<<std::endl;
            } // end of prallel section
            if(iwg!=0){
                ret += pfunc_retrieve(tag, n_walk_prev, mw_info.n_epi_prev_, (Tforce**)mw_info.force_prev_);
            }
            //F64 mass1 = 0.0;
            for(S32 iw=0; iw<n_walk; iw++){
                //mass1 = 0.0;
                S32 ith = iw2ith_adr[iw].first;
                mw_info.epj_[iw] = epj_for_force_[ith].getPointer(iw2ith_adr[iw].second.first);
                mw_info.spj_[iw] = spj_for_force_[ith].getPointer(iw2ith_adr[iw].second.second);
                for(S32 i=0; i<mw_info.n_epj_[iw]; i++){
                    //mass1 += (mw_info.epj_[iw]+i)->mass;
                    //std::cout<<"(mw_info.epj_[iw]+i)->pos= "<<(mw_info.epj_[iw]+i)->pos<<std::endl;
                }
                for(S32 i=0; i<mw_info.n_spj_[iw]; i++){
                    //mass1 += (mw_info.spj_[iw]+i)->mass;
                    //std::cout<<"(mw_info.spj_[iw]+i)->pos= "<<(mw_info.spj_[iw]+i)->pos<<std::endl;
                }
            }

            F64 eps2 = Tepi::eps*Tepi::eps;
            //std::cerr<<"eps2= "<<eps2<<std::endl;
            for(S32 iw=0; iw<n_walk; iw++){
                for(S32 ip=0; ip<mw_info.n_epi_[iw]; ip++){
                    //F64 pot_tmp0 = 0.0;
                    for(S32 jp=0; jp<mw_info.n_epj_[iw]; jp++){
                        F64vec rij = mw_info.epi_[iw][ip].pos - mw_info.epj_[iw][jp].pos;
                        //pot_tmp0 -= mw_info.epi_[iw][ip].mass*mw_info.epj_[iw][jp].mass / sqrt(rij*rij+eps2);
                        //pot_tmp0 -= mw_info.epj_[iw][jp].mass / sqrt(rij*rij+eps2);
                    }
                    for(S32 jp=0; jp<mw_info.n_spj_[iw]; jp++){
                        F64vec rij = mw_info.epi_[iw][ip].pos - mw_info.spj_[iw][jp].pos;
                        //pot_tmp0 -= mw_info.epi_[iw][ip].mass*mw_info.spj_[iw][jp].mass / sqrt(rij*rij+eps2);
                        //pot_tmp0 -= mw_info.spj_[iw][jp].mass / sqrt(rij*rij+eps2);
                    }
                    //pot_tmp0 += mw_info.epi_[iw][ip].mass*mw_info.epi_[iw][ip].mass / sqrt(eps2);
                    //std::cout<<"pot_tmp0= "<<pot_tmp0<<std::endl;
                    //pot0 += pot_tmp0;
                }
            }
            ret += pfunc_dispatch(tag, n_walk,
                                  (const Tepi**)mw_info.epi_, mw_info.n_epi_,
                                  (const Tepj**)mw_info.epj_, mw_info.n_epj_,
                                  (const Tspj**)mw_info.spj_, mw_info.n_spj_);
            for(S32 iw=0; iw<n_walk; iw++){
                mw_info.n_epi_prev_[iw] = mw_info.n_epi_[iw];
                mw_info.force_prev_[iw] = mw_info.force_[iw];
            }
            n_walk_prev = n_walk;
        }
        ret += pfunc_retrieve(tag, n_walk_prev, mw_info.n_epi_prev_, (Tforce**)mw_info_.force_prev_);
        copyForceOriginalOrder();
        return ret;
    }



    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_dispatch, class Tfunc_retrieve>
    S32 TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceUsingIdListMultiWalkIndex(Tfunc_dispatch pfunc_dispatch,
                                       Tfunc_retrieve pfunc_retrieve,
                                       const S32 n_walk_limit,
                                       MultiWalkInfo<Tepi, Tepj, Tspj, Tforce> & mw_info,
                                       const bool clear_force){
        std::vector< std::pair<S32, std::pair<S32, S32> > > iw2ith_adr;
        iw2ith_adr.reserve(n_walk_limit);
        S32 ret = 0;
        S32 tag = 0;
        mw_info.initialize(n_walk_limit);
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        force_org_.resizeNoInitialize(n_loc_tot_);
        if(clear_force){
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for
#endif
            for(S32 ip=0; ip<n_loc_tot_; ip++){
                force_org_[ip].clear();
                force_sorted_[ip].clear();
            }
        }
        pfunc_dispatch(tag, 0,
                       (const Tepi**)mw_info.epi_,    mw_info.n_epi_,
                       (const S32**)mw_info.adr_epj_, mw_info.n_epj_,
                       (const S32**)mw_info.adr_spj_, mw_info.n_spj_,
                       epj_sorted_.getPointer(), epj_sorted_.size(),
                       spj_sorted_.getPointer(), spj_sorted_.size(),
                       true);
        
        S32 n_walk_prev;
        const S64 n_ipg = ipg_.size();
        const S32 n_walk_grp = mw_info.getNWalkGrp(n_ipg);
        F64 pot0 = 0.0; // for debug
        //std::cerr<<"n_ipg= "<<n_ipg<<" n_walk_limit= "<<n_walk_limit<<" n_walk_grp= "<<n_walk_grp<<std::endl;
        for(S32 iwg=0; iwg<n_walk_grp; iwg++){
            const S32 n_walk = n_ipg/n_walk_grp + ( (iwg < n_ipg%n_walk_grp) ? 1 : 0 );
            const S32 walk_grp_head = (n_ipg/n_walk_grp)*iwg + ( (iwg < n_ipg%n_walk_grp) ? iwg : n_ipg%n_walk_grp);
            //std::cerr<<"n_walk= "<<n_walk<<" walk_grp_head= "<<walk_grp_head<<std::endl;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel
#endif
            {
                const S32 my_thread_id = Comm::getThreadNum();
                adr_epj_for_force_[my_thread_id].clearSize();
                adr_spj_for_force_[my_thread_id].clearSize();
                S32 cnt_epj = 0;
                S32 cnt_spj = 0;
                F64 mass0 = 0.0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4) 
#endif	    
                for(S32 iw=0; iw<n_walk; iw++){
                    mass0 = 0.0;
                    iw2ith_adr[iw].first = my_thread_id;
                    const S32 id_ipg = walk_grp_head + iw;
                    // copy epi
                    mw_info.n_epi_[iw] = ipg_[id_ipg].n_ptcl_;
                    const S32 adr_epi  = ipg_[id_ipg].adr_ptcl_;
                    mw_info.epi_[iw]   = epi_sorted_.getPointer(adr_epi);
                    mw_info.force_[iw] = force_sorted_.getPointer(adr_epi);
                    // copy epj
                    const S32 ith      = id_recorder_for_interaction_[id_ipg].ith_;
                    mw_info.n_epj_[iw] = id_recorder_for_interaction_[id_ipg].n_epj_;
                    const S32 adr_epj  = id_recorder_for_interaction_[id_ipg].adr_epj_;
                    iw2ith_adr[iw].second.first = cnt_epj;
                    cnt_epj += mw_info.n_epj_[iw];
                    for(S32 jp=0; jp<mw_info.n_epj_[iw]; jp++){
                        const S32 adr_j = id_epj_recorder_for_force_[ith][adr_epj+jp];
                        adr_epj_for_force_[my_thread_id].push_back( adr_j );
                        //epj_for_force_[my_thread_id].push_back( epj_sorted_[adr_j] );
                        mass0 += epj_sorted_[adr_j].mass;
                    }
                    // copy spj
                    mw_info.n_spj_[iw] = id_recorder_for_interaction_[id_ipg].n_spj_;
                    const S32 adr_spj  = id_recorder_for_interaction_[id_ipg].adr_spj_;
                    iw2ith_adr[iw].second.second = cnt_spj;
                    cnt_spj += mw_info.n_spj_[iw];
                    for(S32 jp=0; jp<mw_info.n_spj_[iw]; jp++){
                        const S32 adr_j = id_spj_recorder_for_force_[ith][adr_spj+jp];
                        adr_spj_for_force_[my_thread_id].push_back( adr_j );
                        //spj_for_force_[my_thread_id].push_back( spj_sorted_[adr_j] );
                        mass0 += spj_sorted_[adr_j].mass;
                    }
                }
                //std::cout<<"mass0= "<<mass0<<std::endl;
            } // end of prallel section
            if(iwg!=0){
                ret += pfunc_retrieve(tag, n_walk_prev, mw_info.n_epi_prev_, (Tforce**)mw_info.force_prev_);
            }
            for(S32 iw=0; iw<n_walk; iw++){
                S32 ith = iw2ith_adr[iw].first;
                //mw_info.epj_[iw] = epj_for_force_[ith].getPointer(iw2ith_adr[iw].second.first);
                //mw_info.spj_[iw] = spj_for_force_[ith].getPointer(iw2ith_adr[iw].second.second);
                mw_info.adr_epj_[iw] = adr_epj_for_force_[ith].getPointer(iw2ith_adr[iw].second.first);
                mw_info.adr_spj_[iw] = adr_spj_for_force_[ith].getPointer(iw2ith_adr[iw].second.second);                
            }
            pfunc_dispatch(tag, n_walk,
                           (const Tepi**)mw_info.epi_, mw_info.n_epi_,
                           (const S32**)mw_info.adr_epj_, mw_info.n_epj_,
                           (const S32**)mw_info.adr_spj_, mw_info.n_spj_,
                           epj_sorted_.getPointer(), epj_sorted_.size(),
                           spj_sorted_.getPointer(), spj_sorted_.size(),
                           false);
            for(S32 iw=0; iw<n_walk; iw++){
                mw_info.n_epi_prev_[iw] = mw_info.n_epi_[iw];
                mw_info.force_prev_[iw] = mw_info.force_[iw];
            }
            n_walk_prev = n_walk;
        }
        ret += pfunc_retrieve(tag, n_walk_prev, mw_info.n_epi_prev_, (Tforce**)mw_info_.force_prev_);
        copyForceOriginalOrder();
        return ret;
    }


    /////////////////////////
    /// add moment as spj ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    addMomentAsSpLocalTreeImpl(TagForceLong){
        spj_sorted_loc_.resizeNoInitialize(tc_loc_.size());
        for(S32 i=0; i<tc_loc_.size(); i++){
            spj_sorted_loc_[i].copyFromMoment(tc_loc_[i].mom_);
        }
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    addMomentAsSpLocalTreeImpl(TagForceShort){ /* do nothing */ }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    addMomentAsSpGlobalTreeImpl(TagForceLong){
        S32 n_spj_prev = spj_sorted_.size();
        //std::cerr<<"epi_sorted_.size()= "<<epi_sorted_.size()<<std::endl;
        //std::cerr<<"epj_sorted_.size()= "<<epj_sorted_.size()<<std::endl;
        //std::cerr<<"spj_sorted_.size()= "<<spj_sorted_.size()<<std::endl;
        //std::cerr<<"spj_org_.size()= "<<spj_org_.size()<<std::endl;
        spj_sorted_.resizeNoInitialize(spj_sorted_.size()+tc_glb_.size());
        for(S32 i=0; i<tc_glb_.size(); i++){
            spj_sorted_[n_spj_prev+i].copyFromMoment(tc_glb_[i].mom_);
        }
    }
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    addMomentAsSpGlobalTreeImpl(TagForceShort){ /* do nothing */ }



    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeAllInteractionListId(const bool clear){
        const S64 n_ipg = ipg_.size();
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
        const S32 sp_tree_offset = spj_sorted_.size() - tc_glb_.size();
        if(n_ipg > 0 && clear){
            id_recorder_for_interaction_.resizeNoInitialize(n_ipg);
            for(S32 ith=0; ith<Comm::getNumberOfThread(); ith++){
                id_epj_recorder_for_force_[ith].clearSize();
                id_spj_recorder_for_force_[ith].clearSize();
            }
        }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for schedule(dynamic, 4)
#endif	    
        for(S32 ig=0; ig<n_ipg; ig++){
            F64 mass = 0.0;
            const S32 ith = Comm::getThreadNum();
            const S32 adr_epj_prev = id_epj_recorder_for_force_[ith].size();
            const S32 adr_spj_prev = id_spj_recorder_for_force_[ith].size();
            const F64ort pos_target_box = (ipg_[ig]).vertex_;
            ReallocatableArray<Tepj> ep_tmp;
            ReallocatableArray<Tspj> sp_tmp;
            MakeListUsingTreeRecursive<TSM, MAKE_LIST_MODE_INTERACTION, LIST_CONTENT_ID, TreeCell<Tmomglb>, TreeParticle, Tepj, Tspj>
                (tc_glb_, 0, tp_glb_, 
                 epj_sorted_, ep_tmp, id_epj_recorder_for_force_[ith], 
                 spj_sorted_, sp_tmp, id_spj_recorder_for_force_[ith], 
                 pos_target_box, r_crit_sq, n_leaf_limit_, sp_tree_offset);
            const S32 n_epj = id_epj_recorder_for_force_[ith].size() - adr_epj_prev;
            const S32 n_spj = id_spj_recorder_for_force_[ith].size() - adr_spj_prev;
            /*
            for(S32 jp=adr_epj_prev; jp<id_epj_recorder_for_force_[ith].size(); jp++){
                //if(Comm::getRank() == 0) std::cerr<<"id_epj_recorder_for_force_[ith][jp]= "<<id_epj_recorder_for_force_[ith][jp]<<std::endl;
                mass += epj_sorted_[ id_epj_recorder_for_force_[ith][jp] ].getCharge();
            }
            for(S32 jp=adr_spj_prev; jp<id_spj_recorder_for_force_[ith].size(); jp++){
                //if(Comm::getRank() == 0) std::cerr<<"id_spj_recorder_for_force_[ith][jp]= "<<id_spj_recorder_for_force_[ith][jp]<<std::endl;
                mass += spj_sorted_[ id_spj_recorder_for_force_[ith][jp] ].getCharge();
            }
            //if(Comm::getRank() == 0)std::cerr<<"mass= "<<mass<<std::endl;
            */
            id_recorder_for_interaction_[ig].set(ig, ith, n_epj, adr_epj_prev, n_spj, adr_spj_prev);
        }
    }



    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    exchangeLocalEssentialTree2(const DomainInfo & dinfo, const bool reuse_list){
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 my_rank = Comm::getRank();
        if(!reuse_list){
            comm_table_.n_ep_send_[my_rank] = comm_table_.n_sp_send_[my_rank] = comm_table_.n_ep_recv_[my_rank] = comm_table_.n_sp_recv_[my_rank] = 0;
        }
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif	
        {
            const S32 ith = Comm::getThreadNum();
            if(!reuse_list){
                const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
                rank_proc_send_[ith].resizeNoInitialize(0);
                id_ep_send_buf_[ith].resizeNoInitialize(0);
                id_sp_send_buf_[ith].resizeNoInitialize(0);
                ReallocatableArray<Tepj> ep_list_dummy;
                ReallocatableArray<Tspj> sp_list_dummy, sp_list_dummy2;
                S32 n_ep_cum = 0;
                S32 n_sp_cum = 0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp for schedule(dynamic, 4)
#endif	    
                for(S32 ib=0; ib<n_proc; ib++){
                    rank_proc_send_[ith].push_back(ib);
                    n_ep_send_[ib] = n_sp_send_[ib] = n_ep_recv_[ib] = n_sp_recv_[ib] = 0;
                    if(my_rank == ib) continue;
                    const F64ort pos_target_domain = dinfo.getPosDomain(ib);
                    MakeListUsingTreeRecursive<TSM, MAKE_LIST_MODE_LET, LIST_CONTENT_ID, TreeCell<Tmomloc>, TreeParticle, Tepj, Tspj>
                        (tc_loc_, 0, tp_loc_, epj_sorted_, ep_list_dummy,  id_ep_send_buf_[ith], sp_list_dummy, sp_list_dummy2, id_sp_send_buf_[ith], pos_target_domain, r_crit_sq, n_leaf_limit_);
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
            S32 n_ep_cnt = 0;
            S32 n_sp_cnt = 0;
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
        }
        comm_table_.exchangeLet(reuse_list);
        /*
        if(Comm::getRank() == 0){
            for(S32 i=0; i<comm_table_.ep_send_.size(); i++){
                std::cout<<"i= "<<i<<" ep_send_[i].pos= "<<comm_table_.ep_send_[i].pos<<" mass= "<<comm_table_.ep_send_[i].mass<<std::endl;
            }
            for(S32 i=0; i<comm_table_.sp_send_.size(); i++){
                std::cout<<"i= "<<i<<" sp_send_[i].pos= "<<comm_table_.ep_send_[i].pos<<" mass= "<<comm_table_.sp_send_[i].mass<<std::endl;
            }
            for(S32 i=0; i<comm_table_.ep_recv_.size(); i++){
                std::cout<<"i= "<<i<<" ep_recv_[i].pos= "<<comm_table_.ep_recv_[i].pos<<" mass= "<<comm_table_.ep_recv_[i].mass<<std::endl;
            }
            for(S32 i=0; i<comm_table_.sp_recv_.size(); i++){
                std::cout<<"i= "<<i<<" sp_recv_[i].pos= "<<comm_table_.ep_recv_[i].pos<<" mass= "<<comm_table_.sp_recv_[i].mass<<std::endl;
            }
        }
        */
    }
}
