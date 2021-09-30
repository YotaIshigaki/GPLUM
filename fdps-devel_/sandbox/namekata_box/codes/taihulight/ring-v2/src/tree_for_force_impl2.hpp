#include<multi_walk.hpp>

#ifdef SUNWAY
extern "C"{
    #include<athread.h>
    #include<cpe_func.h>
    void SLAVE_FUN(CopyIndirect)(void *);
    void SLAVE_FUN(CopyDirect)(void *);
    void SLAVE_FUN(CopyStride)(void *);
    void SLAVE_FUN(CopyTCMomToSPJ)(void *);
    void SLAVE_FUN(CopyETCMomToSPJ)(void *);
    void SLAVE_FUN(SetEpjTpFromBufferCpe)(void *);
    void SLAVE_FUN(SetSpjTpFromBufferCpe)(void *);
}
#endif

namespace ParticleSimulator{
    template<class Ttparray, class Tepjarray0, class Tepjarray1>
    void SetEpjTpFromBuffer(const S32 n,
                            const S32 offset,
                            Ttparray & tp_glb,
                            Tepjarray0 & epj_org,
                            Tepjarray1 & epj_buf,
                            const F64vec & center,
                            const F64 & hlen){
#ifdef SUNWAY
        unsigned long arg[7];
        arg[0] = (unsigned long)(n);
        arg[1] = (unsigned long)(offset);
        arg[2] = (unsigned long)(tp_glb.getPointer());
        arg[3] = (unsigned long)(epj_org.getPointer());
        arg[4] = (unsigned long)(epj_buf.getPointer());
        arg[5] = (unsigned long)(&center);
        arg[6] = (unsigned long)(&hlen);
        __real_athread_spawn((void*)slave_SetEpjTpFromBufferCpe,
                             arg);
        athread_join();
        for(S32 i=0; i<n; i++){
            //epj_org[offset+i] = epj_buf[i];
            tp_glb[offset+i].setFromEP(epj_buf[i], offset+i);
        }
#else
    #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
    #endif
        for(S32 i=0; i<n; i++){
            epj_org[offset+i] = epj_buf[i];
            tp_glb[offset+i].setFromEP(epj_buf[i], offset+i);
        }        
#endif
    }
    
    template<class Ttparray, class Tspjarray0, class Tspjarray1>
    void SetSpjTpFromBuffer(const S32 n,
                            const S32 offset,
                            Ttparray & tp_glb,
                            Tspjarray0 & spj_org,
                            Tspjarray1 & spj_buf,
                            const F64vec & center,
                            const F64 & hlen){
#ifdef SUNWAY
        unsigned long arg[7];
        arg[0] = (unsigned long)(n);
        arg[1] = (unsigned long)(offset);
        arg[2] = (unsigned long)(tp_glb.getPointer());
        arg[3] = (unsigned long)(spj_org.getPointer());
        arg[4] = (unsigned long)(spj_buf.getPointer());
        arg[5] = (unsigned long)(&center);
        arg[6] = (unsigned long)(&hlen);
        __real_athread_spawn((void*)slave_SetSpjTpFromBufferCpe,
                             arg);
        athread_join();
#else
    #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
    #endif
    #ifdef REDUCE_MEMORY
        for(S32 i=0; i<n; i++){
            spj_org[i] = spj_buf[i];
            tp_glb[offset+i].setFromSP(spj_buf[i], i);
        }
    #else // REDUCE_MEMORY
        for(S32 i=0; i<n; i++){
            spj_org[offset+i] = spj_buf[i];
            tp_glb[offset+i].setFromSP(spj_buf[i], offset+i);
        }
    #endif
#endif
    }
}

namespace ParticleSimulator{
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    setLocalEssentialTreeToGlobalTree2(){
        const F64 time_offset = GetWtime();
        const S32 n_proc = Comm::getNumberOfProc();
        const S32 offset = this->n_loc_tot_;
#ifdef DEBUG_PRINT_SETLET_2
        Comm::barrier();
        if(Comm::getRank() == 0){
            std::cerr<<"this->comm_table_.ep_recv_.size()= "<<this->comm_table_.ep_recv_.size()<<std::endl;
            std::cerr<<"this->comm_table_.sp_recv_.size()= "<<this->comm_table_.sp_recv_.size()<<std::endl;
        }
#endif // DEBUG_PRINT_SETLET_2
#ifdef USE_SUPER_DOMAIN
        const S32 n_ep_add = this->comm_table_.ep_recv_.size();
        S32 n_sp_add = 0;
        if( (FORCE_TYPE)(typename TSM::force_type()).force_type == FORCE_TYPE_LONG){
            n_sp_add = this->comm_table_.sp_recv_.size();
        }
#else // USE_SUPER_DOMAIN
        const S32 n_ep_add = this->comm_table_.n_disp_ep_recv_[n_proc];
        S32 n_sp_add = 0;
        if( (FORCE_TYPE)(typename TSM::force_type()).force_type == FORCE_TYPE_LONG){
            n_sp_add = this->comm_table_.n_disp_sp_recv_[n_proc];
        }
#endif // USE_SUPER_DOMAIN
        
#ifdef DEBUG_PRINT_SETLET_2
        Comm::barrier();
        if(Comm::getRank() == 0){
            std::cerr<<"n_ep_add= "<<n_ep_add<<std::endl;
            std::cerr<<"n_sp_add= "<<n_sp_add<<std::endl;
        }
#endif // DEBUG_PRINT_SETLET_2
        const S32 offset2 = this->n_loc_tot_ + n_ep_add;
        this->n_glb_tot_ = this->n_loc_tot_ + n_ep_add + n_sp_add;
        try{
            this->tp_glb_.resizeNoInitialize( this->n_glb_tot_ );
        }
        catch(std::bad_alloc e){
            std::cerr<<"check 1: "<<e.what()<<std::endl;
            std::cerr<<"n_glb_tot_= "<<n_glb_tot_<<std::endl;
            std::cerr<<"rank= "<<Comm::getRank()<<std::endl;
            exit(1);
        }
        try{
            this->epj_org_.resizeNoInitialize( offset2 );
        }
        catch(std::bad_alloc e){
            std::cerr<<"check 2: "<<e.what()<<std::endl;
            std::cerr<<"rank= "<<Comm::getRank()<<std::endl;
            std::cerr<<"this->n_loc_tot_= "<<this->n_loc_tot_<<std::endl;
            std::cerr<<"n_ep_add= "<<n_ep_add<<std::endl;
            std::cerr<<"n_sp_add= "<<n_sp_add<<std::endl;
            std::cerr<<"offset2= "<<offset2<<std::endl;
            exit(1);
        }
        try{
            if( (FORCE_TYPE)(typename TSM::force_type()).force_type == FORCE_TYPE_LONG){
#ifdef REDUCE_MEMORY
                this->spj_org_.resizeNoInitialize( n_sp_add );
#else // REDUCE_MEMORY
                this->spj_org_.resizeNoInitialize( this->n_glb_tot_ );
#endif // REDUCE_MEMORY
            }
            else{
                this->spj_org_.resizeNoInitialize( 0 );
            }
        }
        catch(std::bad_alloc e){
            std::cerr<<"check 3: "<<e.what()<<std::endl;
            std::cerr<<"rank= "<<Comm::getRank()<<std::endl;
            std::cerr<<"n_glb_tot_= "<<n_glb_tot_<<std::endl;
            exit(1);
        }
#ifdef DEBUG_PRINT_SETLET_2
        Comm::barrier();
        if(Comm::getRank() == 0)  std::cerr<<"n_glb_tot_="<<this->n_glb_tot_<<std::endl;
#endif // DEBUG_PRINT_SETLET_2
        
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel
#endif
        {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
            for(S32 i=0; i<offset; i++){
#ifndef REMOVE_TP_LOC
                this->tp_glb_[i] = this->tp_loc_[i]; // NOTE: need to keep tp_loc_[]?
#endif
#ifdef REDUCE_MEMORY
                epj_org_[i] = epj_sorted_loc_[i];
#endif
            }
#ifdef DEBUG_PRINT_SETLET_2
            Comm::barrier();
            if(Comm::getRank() == 0)  std::cerr<<"check 1"<<std::endl;
#endif
            SetEpjTpFromBuffer(n_ep_add, offset, tp_glb_, epj_org_,
                               comm_table_.ep_recv_,
                               this->getCenterOfTree(),
                               this->getHalfLengthOfTree());
            /*
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
            for(S32 i=0; i<n_ep_add; i++){
                this->epj_org_[offset+i] = this->comm_table_.ep_recv_[i];
                this->tp_glb_[offset+i].setFromEP(this->comm_table_.ep_recv_[i], offset+i);
            }
            */


            
#ifdef DEBUG_PRINT_SETLET_2
            Comm::barrier(); if(Comm::getRank() == 0)  std::cerr<<"check 2"<<std::endl;
#endif

            SetSpjTpFromBuffer(n_sp_add, offset2, tp_glb_,
                               spj_org_,
                               comm_table_.sp_recv_,
                               this->getCenterOfTree(),
                               this->getHalfLengthOfTree());
            
            /*            
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for
#endif
            for(S32 i=0; i<n_sp_add; i++){
#ifdef REDUCE_MEMORY
                this->spj_org_[i] = this->comm_table_.sp_recv_[i];
                this->tp_glb_[offset2+i].setFromSP(this->comm_table_.sp_recv_[i], i);
#else
                this->spj_org_[offset2+i] = this->comm_table_.sp_recv_[i];
                this->tp_glb_[offset2+i].setFromSP(this->comm_table_.sp_recv_[i], offset2+i);
#endif
            }
            */
        }


        
#ifdef DEBUG_PRINT_SETLET_2
        Comm::barrier();
        if(Comm::getRank() == 0){
            std::cerr<<"OK03 @DEBUG_PRINT_SETLET_2"<<std::endl;
            for(S32 i=0; i<tp_glb_.size(); i++){
                std::cerr<<"i= "<<i
                         <<" tp_glb_[i].key_= "<<tp_glb_[i].key_
                         <<" adr_ptcl= "<<tp_glb_[i].adr_ptcl_
                         <<std::endl;
                if(GetMSB(tp_glb_[i].adr_ptcl_) == 0){
                    std::cerr<<"adr_ptcl_2= "<<tp_glb_[i].adr_ptcl_
                             <<" epj.mass= "<<epj_org_[tp_glb_[i].adr_ptcl_].mass
                             <<" epj.pos= "<<epj_org_[tp_glb_[i].adr_ptcl_].pos
                             <<std::endl;
                }
                else{
                    std::cerr<<"adr_ptcl_2= "<<ClearMSB(tp_glb_[i].adr_ptcl_)
                             <<" spj.mass= "<<epj_org_[ClearMSB(tp_glb_[i].adr_ptcl_)].mass
                             <<" spj.pos= "<<epj_org_[ClearMSB(tp_glb_[i].adr_ptcl_)].pos
                             <<std::endl;
                }
            }
        }
        Comm::barrier();
#endif
        time_profile_.set_particle_global_tree += GetWtime() - time_offset;
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    template<class Tfunc_ep_ep, class Tfunc_ep_sp>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    calcForceUsingIdList(Tfunc_ep_ep pfunc_ep_ep,
                         Tfunc_ep_sp pfunc_ep_sp,
                         const bool clear_force){
        F64 time_offset = GetWtime();
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
        time_profile_.calc_force += GetWtime() - time_offset;
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
        F64 time_offset = GetWtime();
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
                    epj_for_force_[my_thread_id].reserveEmptyAreaAtLeast(mw_info.n_epj_[iw]);
                    for(S32 jp=0; jp<mw_info.n_epj_[iw]; jp++){
                        const S32 adr_j = id_epj_recorder_for_force_[ith][adr_epj+jp];
                        epj_for_force_[my_thread_id].pushBackNoCheck( epj_sorted_[adr_j] );
                        //epj_for_force_[my_thread_id].push_back( epj_sorted_[adr_j] );
                        //std::cout<<"epj_sorted_[adr_j].pos= "<<epj_sorted_[adr_j].pos<<std::endl;
                        //mass0 += epj_sorted_[adr_j].mass;
                    }
                    // copy spj
                    mw_info.n_spj_[iw] = id_recorder_for_interaction_[id_ipg].n_spj_;
                    const S32 adr_spj  = id_recorder_for_interaction_[id_ipg].adr_spj_;
                    iw2ith_adr[iw].second.second = cnt_spj;
                    cnt_spj += mw_info.n_spj_[iw];
                    spj_for_force_[my_thread_id].reserveEmptyAreaAtLeast(mw_info.n_spj_[iw]);                    
                    for(S32 jp=0; jp<mw_info.n_spj_[iw]; jp++){
                        const S32 adr_j = id_spj_recorder_for_force_[ith][adr_spj+jp];
                        spj_for_force_[my_thread_id].pushBackNocheck( spj_sorted_[adr_j] );
                        //spj_for_force_[my_thread_id].push_back( spj_sorted_[adr_j] );
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
        time_profile_.calc_force += GetWtime() - time_offset;        
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
        F64 time_offset = GetWtime();        
        std::vector< std::pair<S32, std::pair<S32, S32> > > iw2ith_adr;
        iw2ith_adr.reserve(n_walk_limit);
        S32 ret = 0;
        S32 tag = 0;
        mw_info.initialize(n_walk_limit);
        force_sorted_.resizeNoInitialize(n_loc_tot_);
        force_org_.resizeNoInitialize(n_loc_tot_);
        CountT n_walk_tot = 0;
        CountT n_epi_tot = 0;
        CountT n_epj_tot = 0;
        CountT n_spj_tot = 0;
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
        F64 wtime_offset = GetWtime();
        for(S32 iwg=0; iwg<n_walk_grp; iwg++){
            const S32 n_walk = n_ipg/n_walk_grp + ( (iwg < n_ipg%n_walk_grp) ? 1 : 0 );
            const S32 walk_grp_head = (n_ipg/n_walk_grp)*iwg + ( (iwg < n_ipg%n_walk_grp) ? iwg : n_ipg%n_walk_grp);
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel
#endif
            {
                const S32 my_thread_id = Comm::getThreadNum();
                adr_epj_for_force_[my_thread_id].clearSize();
                adr_spj_for_force_[my_thread_id].clearSize();
                S32 cnt_epj = 0;
                S32 cnt_spj = 0;
                //F64 mass0 = 0.0;
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp for schedule(dynamic, 4) reduction(+: n_epi_tot, n_epj_tot, e_spj_tot)
#endif	    
                for(S32 iw=0; iw<n_walk; iw++){
                    //std::cerr<<"iw= "<<iw<<" n_walk= "<<n_walk<<std::endl;
                    //mass0 = 0.0;
                    iw2ith_adr[iw].first = my_thread_id;
                    const S32 id_ipg = walk_grp_head + iw;
                    // copy epi
                    mw_info.n_epi_[iw] = ipg_[id_ipg].n_ptcl_;
                    n_epi_tot += ipg_[id_ipg].n_ptcl_;
                    const S32 adr_epi  = ipg_[id_ipg].adr_ptcl_;
                    mw_info.epi_[iw]   = epi_sorted_.getPointer(adr_epi);
                    mw_info.force_[iw] = force_sorted_.getPointer(adr_epi);
                    
                    // copy epj
                    const S32 ith      = id_recorder_for_interaction_[id_ipg].ith_;
                    mw_info.n_epj_[iw] = id_recorder_for_interaction_[id_ipg].n_epj_;
                    const S32 adr_epj  = id_recorder_for_interaction_[id_ipg].adr_epj_;
                    iw2ith_adr[iw].second.first = cnt_epj;
                    cnt_epj += mw_info.n_epj_[iw];
#ifdef MY_UNROLL
                    const S32 offset_epj = adr_epj_for_force_[my_thread_id].size();
                    adr_epj_for_force_[my_thread_id].resizeNoInitialize(offset_epj+mw_info.n_epj_[iw]);
                    adr_epj_for_force_[my_thread_id].reserveEmptyAreaAtLeast(8);
                    id_epj_recorder_for_force_[ith].reserveEmptyAreaAtLeast(8);
                    const S32 epj_loop_tail = mw_info.n_epj_[iw]/8+1;
                    for(S32 jp=0; jp<epj_loop_tail; jp++){
                        const S32 loc_src = adr_epj + 8*jp;
                        const S32 adr_j_0 = id_epj_recorder_for_force_[ith][loc_src+0];
                        const S32 adr_j_1 = id_epj_recorder_for_force_[ith][loc_src+1];
                        const S32 adr_j_2 = id_epj_recorder_for_force_[ith][loc_src+2];
                        const S32 adr_j_3 = id_epj_recorder_for_force_[ith][loc_src+3];
                        const S32 adr_j_4 = id_epj_recorder_for_force_[ith][loc_src+4];
                        const S32 adr_j_5 = id_epj_recorder_for_force_[ith][loc_src+5];
                        const S32 adr_j_6 = id_epj_recorder_for_force_[ith][loc_src+6];
                        const S32 adr_j_7 = id_epj_recorder_for_force_[ith][loc_src+7];
                        
                        const S32 loc_dst = offset_epj + 8*jp;
                        adr_epj_for_force_[my_thread_id][loc_dst+0] = adr_j_0;
                        adr_epj_for_force_[my_thread_id][loc_dst+1] = adr_j_1;
                        adr_epj_for_force_[my_thread_id][loc_dst+2] = adr_j_2;
                        adr_epj_for_force_[my_thread_id][loc_dst+3] = adr_j_3;
                        adr_epj_for_force_[my_thread_id][loc_dst+4] = adr_j_4;
                        adr_epj_for_force_[my_thread_id][loc_dst+5] = adr_j_5;
                        adr_epj_for_force_[my_thread_id][loc_dst+6] = adr_j_6;
                        adr_epj_for_force_[my_thread_id][loc_dst+7] = adr_j_7;
                    }
                    //std::cerr<<" check b"<<std::endl;
#else
                    adr_epj_for_force_[my_thread_id].reserveEmptyAreaAtLeast(mw_info.n_epj_[iw]);
                    for(S32 jp=0; jp<mw_info.n_epj_[iw]; jp++){
                        const S32 adr_j = id_epj_recorder_for_force_[ith][adr_epj+jp];
                        adr_epj_for_force_[my_thread_id].pushBackNoCheck( adr_j );
                        //adr_epj_for_force_[my_thread_id].push_back( adr_j );
                        //mass0 += epj_sorted_[adr_j].mass;
                    }
#endif
                    // copy spj
                    mw_info.n_spj_[iw] = id_recorder_for_interaction_[id_ipg].n_spj_;
                    const S32 adr_spj  = id_recorder_for_interaction_[id_ipg].adr_spj_;
                    iw2ith_adr[iw].second.second = cnt_spj;
                    cnt_spj += mw_info.n_spj_[iw];
                    //std::cerr<<" check c"<<std::endl;
#ifdef MY_UNROLL
                    const S32 offset_spj = adr_spj_for_force_[my_thread_id].size();
                    adr_spj_for_force_[my_thread_id].resizeNoInitialize(offset_spj+mw_info.n_spj_[iw]);
                    adr_spj_for_force_[my_thread_id].reserveEmptyAreaAtLeast(8);
                    id_spj_recorder_for_force_[ith].reserveEmptyAreaAtLeast(8);
                    const S32 spj_loop_tail = mw_info.n_spj_[iw]/8+1;
                    for(S32 jp=0; jp<mw_info.n_spj_[iw]/8+1; jp++){
                        const S32 loc_src = adr_spj + 8*jp;
                        const S32 adr_j_0 = id_spj_recorder_for_force_[ith][loc_src+0];
                        const S32 adr_j_1 = id_spj_recorder_for_force_[ith][loc_src+1];
                        const S32 adr_j_2 = id_spj_recorder_for_force_[ith][loc_src+2];
                        const S32 adr_j_3 = id_spj_recorder_for_force_[ith][loc_src+3];
                        const S32 adr_j_4 = id_spj_recorder_for_force_[ith][loc_src+4];
                        const S32 adr_j_5 = id_spj_recorder_for_force_[ith][loc_src+5];
                        const S32 adr_j_6 = id_spj_recorder_for_force_[ith][loc_src+6];
                        const S32 adr_j_7 = id_spj_recorder_for_force_[ith][loc_src+7];
                        const S32 loc_dst = offset_spj + 8*jp;
                        adr_spj_for_force_[my_thread_id][loc_dst+0] = adr_j_0;
                        adr_spj_for_force_[my_thread_id][loc_dst+1] = adr_j_1;
                        adr_spj_for_force_[my_thread_id][loc_dst+2] = adr_j_2;
                        adr_spj_for_force_[my_thread_id][loc_dst+3] = adr_j_3;
                        adr_spj_for_force_[my_thread_id][loc_dst+4] = adr_j_4;
                        adr_spj_for_force_[my_thread_id][loc_dst+5] = adr_j_5;
                        adr_spj_for_force_[my_thread_id][loc_dst+6] = adr_j_6;
                        adr_spj_for_force_[my_thread_id][loc_dst+7] = adr_j_7;
                    }
#else
                    adr_spj_for_force_[my_thread_id].reserveEmptyAreaAtLeast(mw_info.n_spj_[iw]);
                    for(S32 jp=0; jp<mw_info.n_spj_[iw]; jp++){
                        const S32 adr_j = id_spj_recorder_for_force_[ith][adr_spj+jp];
                        adr_spj_for_force_[my_thread_id].pushBackNoCheck( adr_j );
                        //adr_spj_for_force_[my_thread_id].push_back( adr_j );
                        //mass0 += spj_sorted_[adr_j].mass;
                    }
#endif
                    //std::cerr<<" check d"<<std::endl;
                }
                n_epj_tot += cnt_epj;
                n_spj_tot += cnt_spj;
            } // end of prallel section
            time_profile_.set_adr_epj_spj_for_CPE += GetWtime() - wtime_offset;
            //std::cerr<<"check a"<<std::endl;
            if(iwg!=0){
                ret += pfunc_retrieve(tag, n_walk_prev, mw_info.n_epi_prev_, (Tforce**)mw_info.force_prev_);
            }
            //std::cerr<<"check b"<<std::endl;
            for(S32 iw=0; iw<n_walk; iw++){
                S32 ith = iw2ith_adr[iw].first;
                mw_info.adr_epj_[iw] = adr_epj_for_force_[ith].getPointer(iw2ith_adr[iw].second.first);
                mw_info.adr_spj_[iw] = adr_spj_for_force_[ith].getPointer(iw2ith_adr[iw].second.second);                
            }
            //std::cerr<<"check c"<<std::endl;
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
        n_epi_loc_ave_ = n_epi_tot / n_ipg;
        n_epj_loc_ave_ = n_epj_tot / n_ipg;
        n_spj_loc_ave_ = n_spj_tot / n_ipg;
        time_profile_.calc_force += GetWtime() - time_offset;
        return ret;
    }

    /////////////////////////
    /// add moment as spj ///
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    addMomentAsSpLocalTreeImpl(TagForceLong){
        const F64 time_offset = GetWtime();
        spj_sorted_loc_.resizeNoInitialize(tc_loc_.size());
#ifdef SUNWAY
        // TODO
        /*
        unsigned long arg[5];
        arg[0] = (unsigned long)(tc_loc_.size());
        arg[1] = (unsigned long)((int*)(tc_loc_.getPointer())+4);
        arg[2] = (unsigned long)(spj_sorted_loc_.getPointer());
        arg[3] = (unsigned long)(sizeof(spj_sorted_loc_[0]));
        arg[4] = (unsigned long)(4*4+sizeof(F64ort)*2); // in byte
        __real_athread_spawn((void*)slave_CopyStride, arg);
        athread_join();
        */
        //F64 starttime;
        //starttime = MPI::Wtime();
        //for(S32 i=0; i<tc_loc_.size(); i++){
        //    spj_sorted_loc_[i].copyFromMoment(tc_loc_[i].mom_);
        //}
        //if (Comm::getRank() == 0) std::cout << "etime(old) = " << MPI::Wtime()-starttime << std::endl;
        //std::ofstream output_file;
        //output_file.open("result_old.dat");
        //   for (S32 i=0; i<tc_loc_.size();i++) {
        //      output_file << spj_sorted_loc_[i].mass  << "  " 
        //                  << spj_sorted_loc_[i].pos.x << "  "
        //                  << spj_sorted_loc_[i].pos.y << "  "
        //                  << spj_sorted_loc_[i].pos.z << std::endl;
        //   }
        //output_file.close();
        //------- 
        //starttime = MPI::Wtime();
        unsigned long args[4];
        args[0] = (unsigned long)(tc_loc_.size());
        args[1] = (unsigned long)0;
        args[2] = (unsigned long)(tc_loc_.getPointer());
        args[3] = (unsigned long)(spj_sorted_loc_.getPointer());
        __real_athread_spawn((void*)slave_CopyTCMomToSPJ, args);
        athread_join();
        //if (Comm::getRank() == 0) std::cout << "etime(new) = " << MPI::Wtime()-starttime << std::endl;
        //output_file.open("result_new.dat");
        //   for (S32 i=0; i<tc_loc_.size();i++) {
        //      output_file << spj_sorted_loc_[i].mass  << "  " 
        //                  << spj_sorted_loc_[i].pos.x << "  "
        //                  << spj_sorted_loc_[i].pos.y << "  "
        //                  << spj_sorted_loc_[i].pos.z << std::endl;
        //   }
        //output_file.close();
        //* For test run
        //if (Comm::getRank() == 0) std::cout << "End:addMomentAsSpLocalTreeImpl" << std::endl;
        //athread_halt();
        //Finalize();
        //std::exit(0);
#else // SUNWAY
        for(S32 i=0; i<tc_loc_.size(); i++){
            spj_sorted_loc_[i].copyFromMoment(tc_loc_[i].mom_);
        }
#endif // SUNWAY
        time_profile_.add_moment_as_sp_local_tree += GetWtime() - time_offset;
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    addMomentAsSpLocalTreeLeanImpl(TagForceLong){
        const F64 time_offset = GetWtime();
        spj_sorted_loc_.resizeNoInitialize(tc_loc_.size());
#ifdef SUNWAY
        unsigned long args[4];
        args[0] = (unsigned long)(etc_loc_.size());
        args[1] = (unsigned long)0;
        args[2] = (unsigned long)(etc_loc_.getPointer());
        args[3] = (unsigned long)(spj_sorted_loc_.getPointer());
        __real_athread_spawn((void*)slave_CopyETCMomToSPJ, args);
        athread_join();
#else // SUNWAY
        for(S32 i=0; i<tc_loc_.size(); i++){
            spj_sorted_loc_[i].copyFromMoment(tc_loc_[i].mom_);
        }
#endif // SUNWAY
        time_profile_.add_moment_as_sp_local_tree += GetWtime() - time_offset;
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    addMomentAsSpLocalTreeImpl(TagForceShort){
        //std::cerr<<"sort"<<std::endl;
        /* do nothing */
    }

    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    addMomentAsSpGlobalTreeImpl(TagForceLong){
        const F64 time_offset = GetWtime();
#ifdef REDUCE_MEMORY
        //S32 n_spj_prev = spj_sorted_.size();
        S32 n_spj_prev = comm_table_.sp_recv_.size();
#else
        S32 n_spj_prev = n_glb_tot_; // new line
#endif
        //spj_sorted_.resizeNoInitialize(spj_sorted_.size()+tc_glb_.size()); // probably bug
        spj_sorted_.resizeNoInitialize(n_spj_prev + tc_glb_.size());
        //#ifdef SUNWAY
#if 1
        // TO DO
        /*
        unsigned long arg[5];
        arg[0] = (unsigned long)(tc_glb_.size());
        arg[1] = (unsigned long)((int*)(tc_glb_.getPointer())+4);
        arg[2] = (unsigned long)(spj_sorted_.getPointer()+n_spj_prev);
        arg[3] = (unsigned long)(sizeof(spj_sorted_[0]));
        arg[4] = (unsigned long)(4*4+sizeof(F64ort)*2); // in byte
        __real_athread_spawn((void*)slave_CopyStride, arg);
        athread_join();
        */
        /*
        for(S32 i=0; i<tc_glb_.size(); i++){
            assert( spj_sorted_[n_spj_prev+i].mass == tc_glb_[i].mom_.mass );
        }
        */
        //F64 starttime;
        //starttime = MPI::Wtime();
        //for(S32 i=0; i<tc_glb_.size(); i++){
        //    spj_sorted_[n_spj_prev+i].copyFromMoment(tc_glb_[i].mom_);
        //}        
        //if (Comm::getRank() == 0) std::cout << "etime(old) = " << MPI::Wtime()-starttime << std::endl;
        //std::ofstream output_file;
        //output_file.open("result_old_glb.dat");
        //   for (S32 i=0; i<tc_glb_.size();i++) {
        //      output_file << spj_sorted_[n_spj_prev+i].mass  << "  " 
        //                  << spj_sorted_[n_spj_prev+i].pos.x << "  "
        //                  << spj_sorted_[n_spj_prev+i].pos.y << "  "
        //                  << spj_sorted_[n_spj_prev+i].pos.z << std::endl;
        //   }
        //output_file.close();
        //------- 
        //starttime = MPI::Wtime();
        unsigned long args[4];
        args[0] = (unsigned long)(tc_glb_.size());
        args[1] = (unsigned long)(n_spj_prev);
        args[2] = (unsigned long)(tc_glb_.getPointer());
        args[3] = (unsigned long)(spj_sorted_.getPointer());
        __real_athread_spawn((void*)slave_CopyTCMomToSPJ, args);
        athread_join();
        //if (Comm::getRank() == 0) std::cout << "etime(new) = " << MPI::Wtime()-starttime << std::endl;
        //output_file.open("result_new_glb.dat");
        //   for (S32 i=0; i<tc_glb_.size();i++) {
        //      output_file << spj_sorted_[n_spj_prev+i].mass  << "  " 
        //                  << spj_sorted_[n_spj_prev+i].pos.x << "  "
        //                  << spj_sorted_[n_spj_prev+i].pos.y << "  "
        //                  << spj_sorted_[n_spj_prev+i].pos.z << std::endl;
        //   }
        //output_file.close();
        //* For test run
        //if (Comm::getRank() == 0) std::cout << "End:addMomentAsSpGlobalTreeImpl" << std::endl;
        //athread_halt();
        //Finalize();
        //std::exit(0);
#else // #if 1
        for(S32 i=0; i<tc_glb_.size(); i++){
            spj_sorted_[n_spj_prev+i].copyFromMoment(tc_glb_[i].mom_);
        }
#endif // #if 1
        time_profile_.add_moment_as_sp_global_tree += GetWtime() - time_offset;
    }


    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    addMomentAsSpGlobalTreeLeanImpl(TagForceLong){
        const F64 time_offset = GetWtime();
#ifdef REDUCE_MEMORY
        //S32 n_spj_prev = spj_sorted_.size();
        S32 n_spj_prev = comm_table_.sp_recv_.size();
#else
        S32 n_spj_prev = n_glb_tot_;
#endif
        spj_sorted_.resizeNoInitialize(n_spj_prev + tc_glb_.size());
#ifdef SUNWAY
        unsigned long args[4];
        args[0] = (unsigned long)(etc_glb_.size());
        args[1] = (unsigned long)(n_spj_prev);
        args[2] = (unsigned long)(etc_glb_.getPointer());
        args[3] = (unsigned long)(spj_sorted_.getPointer());
        __real_athread_spawn((void*)slave_CopyETCMomToSPJ, args);
        athread_join();
#else
        for(S32 i=0; i<tc_glb_.size(); i++){
            spj_sorted_[n_spj_prev+i].copyFromMoment(tc_glb_[i].mom_);
        }
#endif
        time_profile_.add_moment_as_sp_global_tree += GetWtime() - time_offset;        
    }
    
    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    addMomentAsSpGlobalTreeImpl(TagForceShort){ /* do nothing */ }



    template<class TSM, class Tforce, class Tepi, class Tepj,
             class Tmomloc, class Tmomglb, class Tspj>
    void TreeForForce<TSM, Tforce, Tepi, Tepj, Tmomloc, Tmomglb, Tspj>::
    makeAllInteractionListId(const bool clear){
        const F64 time_offset = GetWtime();        
        const S64 n_ipg = ipg_.size();
        const F64 r_crit_sq = (length_ * length_) / (theta_ * theta_);
#ifdef REDUCE_MEMORY
        const S32 sp_tree_offset =  comm_table_.sp_recv_.size();
#else
        const S32 sp_tree_offset = spj_sorted_.size() - tc_glb_.size();
#endif
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
            id_recorder_for_interaction_[ig].set(ig, ith, n_epj, adr_epj_prev, n_spj, adr_spj_prev);
        }
        time_profile_.make_all_interaction_list_id += GetWtime() - time_offset;
    }
}

