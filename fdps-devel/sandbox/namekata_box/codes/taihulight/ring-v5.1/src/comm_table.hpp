#pragma once

#include<particle_simulator.hpp>

namespace ParticleSimulator{
    F64 wtime_alltoallv_ep;
    F64 wtime_alltoallv_sp;
    F64 wtime_alltoallv_sp_allgather;
    F64 wtime_alltoallv_sp_isendrecv;
    F64 wtime_alltoallv_sp_isendrecv_1st;
    F64 wtime_alltoallv_sp_isendrecv_2nd;
    CountT n_ep_send_tot_loc;
    CountT n_sp_send_tot_loc;
    CountT n_sp_send_tot_isendrecv_loc;

    F64 wtime_ex_let_exchange_number;
    
#ifdef USE_SUPER_DOMAIN
    template<class Tepj, class Tspj>
    class CommTable{
    public:
        ReallocatableArray<S32> rank_sd_sd_send_;
        ReallocatableArray<S32> rank_sd_sd_recv_;
        ReallocatableArray<S32> rank_ptcl_sd_send_;
        ReallocatableArray<S32> rank_ptcl_sd_recv_;
        ReallocatableArray<S32> rank_ptcl_domain_send_;
        ReallocatableArray<S32> rank_ptcl_domain_recv_;

        ReallocatableArray<S32> adr_ep_sd_send_;
        ReallocatableArray<S32> adr_sp_sd_send_;

        ReallocatableArray<Tepj> ep_recv_;
        ReallocatableArray<Tspj> sp_recv_;

        ReallocatableArray<S32> n_ep_domain_send_;
        ReallocatableArray<S32> n_sp_domain_send_;
        ReallocatableArray<S32> n_ep_sd_send_;
        ReallocatableArray<S32> n_sp_sd_send_;

        ReallocatableArray<S32> n_ep_domain_recv_;
        ReallocatableArray<S32> n_sp_domain_recv_;
        ReallocatableArray<S32> n_ep_sd_recv_;
        ReallocatableArray<S32> n_sp_sd_recv_;

        ReallocatableArray<S32> n_disp_ep_domain_send_;
        ReallocatableArray<S32> n_disp_sp_domain_send_;
        ReallocatableArray<S32> n_disp_ep_sd_send_;
        ReallocatableArray<S32> n_disp_sp_sd_send_;

        ReallocatableArray<S32> n_disp_ep_domain_recv_;
        ReallocatableArray<S32> n_disp_sp_domain_recv_;
        ReallocatableArray<S32> n_disp_ep_sd_recv_;
        ReallocatableArray<S32> n_disp_sp_sd_recv_;


        S32 n_ep_domain_send_2_;
        S32 n_sp_domain_send_2_;
        ReallocatableArray<S32> n_ep_domain_recv_2_;
        ReallocatableArray<S32> n_sp_domain_recv_2_;
        ReallocatableArray<S32> n_disp_ep_domain_recv_2_;
        ReallocatableArray<S32> n_disp_sp_domain_recv_2_;
        

        ReallocatableArray<MPI_Request> req_ep_send_;
        ReallocatableArray<MPI_Request> req_ep_recv_;
        ReallocatableArray<MPI_Request> req_sp_send_;
        ReallocatableArray<MPI_Request> req_sp_recv_;
        ReallocatableArray<MPI_Status> stat_ep_send_;
        ReallocatableArray<MPI_Status> stat_ep_recv_;
        ReallocatableArray<MPI_Status> stat_sp_send_;
        ReallocatableArray<MPI_Status> stat_sp_recv_;

        ReallocatableArray<Tepj> ep_send_;
        ReallocatableArray<Tspj> sp_send_;

        ReallocatableArray<S32> adr_ep_domain_send_;
        ReallocatableArray<S32> adr_sp_domain_send_;
        

    public:
        void clearSize(){
            adr_ep_domain_send_.clearSize();
            adr_sp_domain_send_.clearSize();
            
            rank_sd_sd_send_.clearSize();
            rank_sd_sd_recv_.clearSize();
            rank_ptcl_sd_send_.clearSize();
            rank_ptcl_sd_recv_.clearSize();
            rank_ptcl_domain_send_.clearSize();
            rank_ptcl_domain_recv_.clearSize();

            adr_ep_sd_send_.clearSize();
            adr_sp_sd_send_.clearSize();

            ep_recv_.clearSize();
            sp_recv_.clearSize();

            n_ep_domain_send_.clearSize();
            n_sp_domain_send_.clearSize();
            n_ep_sd_send_.clearSize();
            n_sp_sd_send_.clearSize();

            n_ep_domain_recv_.clearSize();
            n_sp_domain_recv_.clearSize();
            n_ep_sd_recv_.clearSize();
            n_sp_sd_recv_.clearSize();

            n_disp_ep_domain_send_.clearSize();
            n_disp_sp_domain_send_.clearSize();
            n_disp_ep_sd_send_.clearSize();
            n_disp_sp_sd_send_.clearSize();

            n_disp_ep_domain_recv_.clearSize();
            n_disp_sp_domain_recv_.clearSize();
            n_disp_ep_sd_recv_.clearSize();
            n_disp_sp_sd_recv_.clearSize();            

            req_ep_send_.clearSize();
            req_ep_recv_.clearSize();
            req_sp_send_.clearSize();
            req_sp_recv_.clearSize();
            stat_ep_send_.clearSize();
            stat_ep_recv_.clearSize();
            stat_sp_send_.clearSize();
            stat_sp_recv_.clearSize();

            n_ep_domain_recv_2_.clearSize();
            n_sp_domain_recv_2_.clearSize();
            n_disp_ep_domain_recv_2_.clearSize();
            n_disp_sp_domain_recv_2_.clearSize();

            ep_send_.clearSize();
            sp_send_.clearSize();

        }

        void dumpNumber(std::ostream & fout=std::cerr){
            fout<<"sd-sd"<<std::endl;
            for(S32 i=0; i<rank_sd_sd_send_.size(); i++){
                fout<<"i= "<<i<<" rank_sd_sd_send_[i]= "<<rank_sd_sd_send_[i]<<std::endl;
            }
            for(S32 i=0; i<rank_sd_sd_recv_.size(); i++){
                fout<<"i= "<<i<<" rank_sd_sd_recv_[i]= "<<rank_sd_sd_recv_[i]<<std::endl;
            }
            fout<<std::endl;
            fout<<"ptcl-domain (direct)"<<std::endl;
            for(S32 i=0; i<rank_ptcl_domain_send_.size(); i++){
                fout<<"i= "<<i
                    <<" rank_ptcl_domain_send_[i]= "<<rank_ptcl_domain_send_[i]
                    <<" n_ep_domain_send_[i]= "<<n_ep_domain_send_[i]
                    <<" n_sp_domain_send_[i]= "<<n_sp_domain_send_[i]
                    <<std::endl;
            }
            for(S32 i=0; i<rank_ptcl_domain_recv_.size(); i++){
                fout<<"i= "<<i
                    <<" rank_ptcl_domain_recv_[i]= "<<rank_ptcl_domain_recv_[i]
                    <<" n_ep_domain_recv_[i]= "<<n_ep_domain_recv_[i]
                    <<" n_sp_domain_recv_[i]= "<<n_sp_domain_recv_[i]
                    <<std::endl;
            }

            fout<<std::endl;
            fout<<"ptcl-sd"<<std::endl;
            for(S32 i=0; i<rank_ptcl_sd_send_.size(); i++){
                fout<<"i= "<<i
                    <<" rank_ptcl_sd_send_[i]= "<<rank_ptcl_sd_send_[i]
                    <<" n_ep_sd_send_[i]= "<<n_ep_sd_send_[i]
                    <<" n_sp_sd_send_[i]= "<<n_sp_sd_send_[i]
                    <<std::endl;
            }
            for(S32 i=0; i<rank_ptcl_sd_recv_.size(); i++){
                fout<<"i= "<<i
                    <<" rank_ptcl_sd_recv_[i]= "<<rank_ptcl_sd_recv_[i]
                    <<" n_ep_sd_recv_[i]= "<<n_ep_sd_recv_[i]
                    <<" n_sp_sd_recv_[i]= "<<n_sp_sd_recv_[i]
                    <<std::endl;
            }

            fout<<std::endl;
            fout<<"ptcl-domain (indirect)"<<std::endl;
            fout<<"n_ep_domain_send_2_= "<<n_ep_domain_send_2_
                <<" n_sp_domain_send_2_= "<<n_sp_domain_send_2_
                <<std::endl;
            S32 n_proc_sub = n_ep_domain_recv_2_.size();
            for(S32 i=0; i<n_proc_sub; i++){
                fout<<"i= "<<i
                    <<" n_ep_domain_recv_2_[i]= "<<n_ep_domain_recv_2_[i]
                    <<" n_sp_domain_recv_2_[i]= "<<n_sp_domain_recv_2_[i]
                    <<std::endl;
            }
            fout<<"dumpNumber() completed."<<std::endl;
        }

        void initialize(){
            /*
            const S32 n_proc = Comm::getNumberOfProc();
            n_ep_send_.reserve(n_proc);
            n_sp_send_.reserve(n_proc);
            n_disp_ep_send_.reserve(n_proc+1);
            n_disp_sp_send_.reserve(n_proc+1);
            n_ep_recv_.reserve(n_proc);
            n_sp_recv_.reserve(n_proc);
            n_disp_ep_recv_.reserve(n_proc+1);
            n_disp_sp_recv_.reserve(n_proc+1);
            */
            /*
            rank_sp_isend_.reserve(n_proc);
            rank_sp_irecv_.reserve(n_proc);
            rank_ep_send_.reserve(n_proc);
            rank_ep_recv_.reserve(n_proc);
            tree_outer_pos_.reserve(n_proc);
            */
            /*
            n_sp_send_suprer_domain_.reserve(n_proc);
            n_sp_recv_suprer_domain_.reserve(n_proc);
            rank_1d_send_super_domain_.reserve(n_proc);
            rank_1d_recv_super_domain_.reserve(n_proc);
            rank_send_.reserve(n_proc);
            rank_recv_.reserve(n_proc);            

            req_ep_send_.reserve(n_proc);
            req_sp_send_.reserve(n_proc);
            req_ep_recv_.reserve(n_proc);
            req_sp_recv_.reserve(n_proc);
            stat_ep_send_.reserve(n_proc);
            stat_sp_send_.reserve(n_proc);
            stat_ep_recv_.reserve(n_proc);
            stat_sp_recv_.reserve(n_proc);
            */
            
            //n_ep_sp_send_.reserve(n_proc*2);
            //n_ep_sp_recv_.reserve(n_proc*2);
            //n_ep_sp_send_   = new S32[n_proc*2];
            //n_ep_sp_recv_   = new S32[n_proc*2];
            //rank_sp_top_recv_.reserve(n_proc);
            //sp_top_recv_ = new Tspj[n_proc];
        }
        /*

        */
        /*
        void setDispSend(){
            const S32 n_proc = Comm::getNumberOfProc();
            n_disp_ep_send_[0] = n_disp_sp_send_[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_ep_send_[i+1] = n_disp_ep_send_[i] + n_ep_send_[i];
                n_disp_sp_send_[i+1] = n_disp_sp_send_[i] + n_sp_send_[i];
            }
        }
        */
        /*
        void setDispRecv(){
            const S32 n_proc = Comm::getNumberOfProc();
            n_disp_ep_recv_[0] = n_disp_sp_recv_[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_ep_recv_[i+1] = n_disp_ep_recv_[i] + n_ep_recv_[i];
                n_disp_sp_recv_[i+1] = n_disp_sp_recv_[i] + n_sp_recv_[i];
            }
        }
        */
        /*
        void setSizeSendBuf(){
            const S32 n_proc = Comm::getNumberOfProc();
            ep_send_.resizeNoInitialize( n_disp_ep_send_[n_proc] );
            sp_send_.resizeNoInitialize( n_disp_sp_send_[n_proc] );
        }
        */
        /*
        void setSizeRecvBuf(){
            const S32 n_proc = Comm::getNumberOfProc();
            ep_recv_.resizeNoInitialize( n_disp_ep_recv_[n_proc] );
            sp_recv_.resizeNoInitialize( n_disp_sp_recv_[n_proc] );
        }
        */
#if 0        
        
        void exchangeLet(const bool reuse = false){
            if(!reuse){
                const S32 n_proc = Comm::getNumberOfProc();
                for(S32 i=0; i<n_proc; i++){
                    n_ep_sp_send_[2*i] = n_ep_send_[i];
                    n_ep_sp_send_[2*i+1] = n_sp_send_[i];
                }
                Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_);
                for(S32 i=0; i<n_proc; i++){
                    n_ep_recv_[i] = n_ep_sp_recv_[2*i];
                    n_sp_recv_[i] = n_ep_sp_recv_[2*i+1];
                }
                setDispRecv();
                setSizeRecvBuf();
            }
            Comm::allToAllV(ep_send_.getPointer(), n_ep_send_, n_disp_ep_send_,
                            ep_recv_.getPointer(), n_ep_recv_, n_disp_ep_recv_);
            Comm::allToAllV(sp_send_.getPointer(), n_sp_send_, n_disp_sp_send_,
                            sp_recv_.getPointer(), n_sp_recv_, n_disp_sp_recv_);
        }

        void setSPTop(const Tspj & _sp_top){
            sp_top_ = _sp_top;
        }

        void allgahterTreeOuterPos(const F64ort & _tree_outer_pos){
            Comm::allGather(&_tree_outer_pos, 1, tree_outer_pos_.getPointer());
        }

        void exchangeLetSunWay(const ReallocatableArray<F64ort> & _outer_pos_super_domain,
                               const F64 _theta,
                               const bool reuse = false,
                               const DomainInfo * dinfo = NULL){
            const S32 n_proc = Comm::getNumberOfProc();
            const S32 n_proc_1d = dinfo.getNDomain(0);
            const S32 n_proc_sub = n_proc / n_proc_1d;
            const S32 my_rank = Comm::getRank();
            MPI_Comm comm_1d = dinfo.getComm1d(0);
            MPI_Comm comm_sub = dinfo.getCommSub(0);
            if(!reuse){
                // new line
                rank_sp_isend_.clearSize();
                rank_sp_irecv_.clearSize();
                rank_ep_send_.clearSize();
                rank_ep_recv_.clearSize();
                rank_sp_top_recv_.clearSize();
                S32 n_cnt_send = 0;
                S32 n_cnt_recv = 0;
                Comm::barrier();
    #ifdef PHI_R_TREE
                const F64 length_peri_x = dinfo->getPosRootDomain().getFullLength()[0];
                //if(my_rank == 0) std::cerr<<"length_peri_x= "<<length_peri_x<<std::endl;
    #endif //PHI_R_TREE
                
    #ifdef SUNWAY_DEVELOP
                //const F64 longest_len_of_outer_box_send = _tree_outer_pos[my_rank].getFullLength().getMax();
                const F64 longest_len_of_outer_box_send = _outer_pos_super_domain[my_rank].getFullLength().getMax();
                const F64 r_crit_sq_send = (longest_len_of_outer_box_send * longest_len_of_outer_box_send) / (_theta * _theta);
                for(S32 i=0; i<n_proc_1d; i++){
                    //if(my_rank == i) continue;
                    //const F64 longest_len_of_outer_box_recv = _tree_outer_pos[i].getFullLength().getMax();
                    const F64 longest_len_of_outer_box_recv = _outer_pos_super_domain[i].getFullLength().getMax();
                    const F64 r_crit_sq_recv = (longest_len_of_outer_box_recv * longest_len_of_outer_box_recv) / (_theta * _theta);
        #ifdef PHI_R_TREE
                    if(GetDistanceMinSqPeriodicX(_outer_pos_super_domain[i], _outer_pos_super_domain[my_rank_1d], length_peri_x) <= r_crit_sq_send){
        #else
                    if(_tree_outer_pos[i].getDistanceMinSQ(_tree_outer_pos[my_rank]) <= r_crit_sq_send){
        #endif
                        for(S32 j=0; j<n_proc_sub; j++){
                            // to be consistent with make LET function in tree_for_force_impl3.hpp
                            S32 target_rank = i*n_proc_sub + j;
                            MPI_Isend(n_sp_send_+target_rank, 1, GetDataType<S32>(), target_rank, 0, MPI_COMM_WORLD, &req_sp_send_[n_cnt_send]);
                            MPI_Isend(n_ep_send_+target_rank, 1, GetDataType<S32>(), target_rank, 0, MPI_COMM_WORLD, &req_ep_send_[n_cnt_send]);
                            n_cnt_send++;
                        }
                    }
                    else{
                        assert(n_sp_send_[i]==1);
                        assert(n_ep_send_[i]==0);
                    }
        #ifdef PHI_R_TREE
                    //if(GetDistanceMinSqPeriodicX(_tree_outer_pos[my_rank], _tree_outer_pos[i], length_peri_x) <= r_crit_sq_recv){
                    if(GetDistanceMinSqPeriodicX(_outer_pos_super_domain[my_rank_1d], _outer_pos_super_domain[i], length_peri_x) <= r_crit_sq_recv){
        #else
                    if(_tree_outer_pos[my_rank].getDistanceMinSQ(_tree_outer_pos[i]) <= r_crit_sq_recv){
        #endif
                        req_sp_recv_[n_cnt_recv] = MPI::COMM_WORLD.Irecv(n_sp_recv_+i, 1, GetDataType<S32>(), i, 0);
                        req_ep_recv_[n_cnt_recv] = MPI::COMM_WORLD.Irecv(n_ep_recv_+i, 1, GetDataType<S32>(), i, 1);
                        n_cnt_recv++;
                    }
                    else{
                        n_sp_recv_[i] = 1;
                        n_ep_recv_[i] = 0;
                        //rank_sp_top_recv_.pushBackNoCheck(i);
                    }
                }
                MPI::Request::Waitall(n_cnt_send, req_ep_send_.getPointer());
                MPI::Request::Waitall(n_cnt_send, req_sp_send_.getPointer());
                MPI::Request::Waitall(n_cnt_recv, req_ep_recv_.getPointer());
                MPI::Request::Waitall(n_cnt_recv, req_sp_recv_.getPointer());

                for(S32 i=0; i<n_proc; i++){
                    if(n_sp_send_[i] > 1){
                        rank_sp_isend_.pushBackNoCheck(i);
                    }
                    if(n_ep_send_[i] > 0) rank_ep_send_.pushBackNoCheck(i);
                }
                for(S32 i=0; i<n_proc; i++){
                    if(n_sp_recv_[i] > 1) rank_sp_irecv_.pushBackNoCheck(i);
                    if(n_ep_recv_[i] > 0) rank_ep_recv_.pushBackNoCheck(i);
                }

                // debug
        #ifdef DEBUG_EXLET_SUNWAY
                ReallocatableArray<S32> rank_sp_isend_tmp;
                ReallocatableArray<S32> rank_ep_send_tmp;
                ReallocatableArray<S32> rank_sp_irecv_tmp;
                ReallocatableArray<S32> rank_ep_recv_tmp;
                for(S32 i=0; i<n_proc; i++){
                    n_ep_sp_send_[2*i]   = n_ep_send_[i];
                    n_ep_sp_send_[2*i+1] = n_sp_send_[i];
                    if(n_sp_send_[i] > 1){
                        rank_sp_isend_tmp.push_back(i);
                    }
                    if(n_ep_send_[i] > 0) rank_ep_send_tmp.push_back(i);
                }
                Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_);
                for(S32 i=0; i<n_proc; i++){
                    if(n_ep_sp_recv_[2*i+1] > 1) rank_sp_irecv_tmp.push_back(i);
                    if(n_ep_sp_recv_[2*i] > 0) rank_ep_recv_tmp.push_back(i);
                }
                assert(rank_sp_irecv_tmp.size() == rank_sp_irecv_.size());
                assert(rank_ep_recv_tmp.size()  == rank_ep_recv_.size());
                assert(rank_sp_isend_tmp.size() == rank_sp_isend_.size());
                assert(rank_ep_send_tmp.size()  == rank_ep_send_.size());
        #endif
                // debug
    #else //SUNWAY_DEVELOP
                for(S32 i=0; i<n_proc; i++){
                    n_ep_sp_send_[2*i]   = n_ep_send_[i];
                    n_ep_sp_send_[2*i+1] = n_sp_send_[i];
                    if(n_sp_send_[i] > 1){
                        rank_sp_isend_.pushBackNoCheck(i);
                    }
                    if(n_ep_send_[i] > 0) rank_ep_send_.pushBackNoCheck(i);
                }
                Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_);
                for(S32 i=0; i<n_proc; i++){
                    n_ep_recv_[i] = n_ep_sp_recv_[2*i];
                    n_sp_recv_[i] = n_ep_sp_recv_[2*i+1];
                    if(n_sp_recv_[i] > 1) rank_sp_irecv_.pushBackNoCheck(i);
                    if(n_ep_recv_[i] > 0) rank_ep_recv_.pushBackNoCheck(i);
                }
    #endif //SUNWAY_DEVELOP
                setDispRecv();
                setSizeRecvBuf();
                
            }
            /////////
            // set EP
            n_ep_send_tot_loc = n_disp_ep_send_[n_proc];
            F64 wtime_offset_ep = GetWtime() ;
    #ifndef SUNWAY_DEVELOP
            //#if 0
            Comm::allToAllV(ep_send_.getPointer(), n_ep_send_, n_disp_ep_send_,
                            ep_recv_.getPointer(), n_ep_recv_, n_disp_ep_recv_);

    #else //SUNWAY_DEVELOP
            const S32 n_proc_ep_send = rank_ep_send_.size();
            const S32 n_proc_ep_recv = rank_ep_recv_.size();
            n_proc_ep_send_ = n_proc_ep_send;
            n_proc_ep_recv_ = n_proc_ep_recv;
            req_ep_send_.resizeNoInitialize(n_proc_ep_send);
            req_ep_recv_.resizeNoInitialize(n_proc_ep_recv);
            //std::cerr<<"my_rank(D)= "<<my_rank<<std::endl;
            for(S32 i=0; i<n_proc_ep_send; i++){
                const S32 rank_send   = rank_ep_send_[i];
                const S32 n_ep_send      = n_ep_send_[rank_send];
                const S32 n_disp_ep_send = n_disp_ep_send_[rank_send];
                const S32 tag = 0;
                Tepj * ep_send = ep_send_.getPointer(n_disp_ep_send);
                req_ep_send_[i] = MPI::COMM_WORLD.Isend
                    (ep_send, n_ep_send, GetDataType<Tepj>(), rank_send, tag);
            }
            for(S32 i=0; i<n_proc_ep_recv; i++){
                const S32 rank_recv   = rank_ep_recv_[i];
                const S32 n_ep_recv      = n_ep_recv_[rank_recv];
                const S32 n_disp_ep_recv = n_disp_ep_recv_[rank_recv];
                const S32 tag = 0;
                Tepj * ep_recv        = ep_recv_.getPointer(n_disp_ep_recv);
                req_ep_recv_[i] = MPI::COMM_WORLD.Irecv
                    (ep_recv, n_ep_recv, GetDataType<Tepj>(), rank_recv, tag);
            }
            MPI::Request::Waitall(n_proc_ep_send, req_ep_send_.getPointer());
            MPI::Request::Waitall(n_proc_ep_recv, req_ep_recv_.getPointer());
    #endif //SUNWAY_DEVELOP
            wtime_alltoallv_ep = GetWtime() - wtime_offset_ep;
            //band_width_ep = n_ep_send_tot_loc*sizeof(Tepj) / wtime_alltoallv_ep;
            // set EP
            /////////

            
            //std::cerr<<"my_rank(E)= "<<my_rank<<std::endl;
            
            /////////
            // set SP
            F64 wtime_offset_sp = GetWtime();
    #ifndef SUNWAY_DEVELOP
            //#if 1
            Comm::allToAllV(sp_send_.getPointer(), n_sp_send_, n_disp_sp_send_,
                            sp_recv_.getPointer(), n_sp_recv_, n_disp_sp_recv_);
    #else //SUNWAY_DEVELOP
            Comm::allGather(&sp_top_, 1, sp_top_recv_);
            wtime_alltoallv_sp_allgather += GetWtime() - wtime_offset_sp;
            
            const S32 n_proc_sp_send = rank_sp_isend_.size();
            const S32 n_proc_sp_recv = rank_sp_irecv_.size();

            n_proc_sp_isend_ = n_proc_sp_send;
            n_proc_sp_irecv_ = n_proc_sp_recv;
            req_sp_send_.resizeNoInitialize(n_proc_sp_send);
            req_sp_recv_.resizeNoInitialize(n_proc_sp_recv);
            n_sp_send_tot_loc = n_disp_sp_send_[n_proc];
            F64 wtime_offset_isendrecv = GetWtime();
            for(S32 i=0; i<n_proc_sp_send; i++){
                const S32 rank_send   = rank_sp_isend_[i];
                const S32 n_sp_send      = n_sp_send_[rank_send];
                n_sp_send_tot_isendrecv_loc += n_sp_send;
                const S32 n_disp_sp_send = n_disp_sp_send_[rank_send];
                const S32 tag = 1;
                Tspj * sp_send = sp_send_.getPointer(n_disp_sp_send);
                req_sp_send_[i] = MPI::COMM_WORLD.Isend
                    (sp_send, n_sp_send, GetDataType<Tspj>(), rank_send, tag);
            }
            for(S32 i=0; i<n_proc_sp_recv; i++){
                const S32 rank_recv   = rank_sp_irecv_[i];
                const S32 n_sp_recv      = n_sp_recv_[rank_recv];
                const S32 n_disp_sp_recv = n_disp_sp_recv_[rank_recv];
                const S32 tag = 1;
                Tspj * sp_recv        = sp_recv_.getPointer(n_disp_sp_recv);
                req_sp_recv_[i] = MPI::COMM_WORLD.Irecv
                    (sp_recv, n_sp_recv, GetDataType<Tspj>(), rank_recv, tag);
            }
            MPI::Request::Waitall(n_proc_sp_send, req_sp_send_.getPointer());
            MPI::Request::Waitall(n_proc_sp_recv, req_sp_recv_.getPointer());

            wtime_alltoallv_sp_isendrecv += GetWtime() - wtime_offset_isendrecv;
            //band_width_sp = n_sp_send_tot_loc*sizeof(Tepj) / wtime_alltoallv_ep;
            
            for(S32 i=0; i<n_proc; i++){
                if(n_sp_recv_[i] == 1){
                    sp_recv_[n_disp_sp_recv_[i]] = sp_top_recv_[i];
                }
            }

            // for debug
        #ifdef DEBUG_EXLET_SUNWAY
            ReallocatableArray<Tspj> sp_recv_tmp;
            sp_recv_tmp.resizeNoInitialize(n_disp_sp_recv_[n_proc]);
            Comm::allToAllV(sp_send_.getPointer(), n_sp_send_, n_disp_sp_send_,
                            sp_recv_tmp.getPointer(), n_sp_recv_, n_disp_sp_recv_);
            for(S32 i=0; i<n_disp_sp_recv_[n_proc]; i++){
                if(sp_recv_tmp[i].mass != sp_recv_[i].mass){
                    std::cerr<<"sp_recv_tmp[i].mass= "<<sp_recv_tmp[i].mass
                             <<" sp_recv_[i].mass= "<<sp_recv_[i].mass
                             <<std::endl;
                }
                assert(sp_recv_tmp[i].mass == sp_recv_[i].mass);
            }
        #endif //DEBUG_EXLET_SUNWAY

    #endif //SUNWAY_DEVELOP
            wtime_alltoallv_sp += GetWtime() - wtime_offset_sp;
            // set SP
            /////////
            //std::cerr<<"my_rank(F)= "<<my_rank<<std::endl;
        }


#endif

        void dumpMemSizeUsed(std::ostream & fout){
            F64 norm_fac=1e-9;
            std::string str_unit=" [GB]";
            F64 sum = (double)(rank_sd_sd_send_.getMemSize()
                             + rank_sd_sd_send_.getMemSize()
                             + rank_sd_sd_recv_.getMemSize()
                             + rank_ptcl_sd_send_.getMemSize()
                             + rank_ptcl_sd_recv_.getMemSize()
                             + rank_ptcl_domain_send_.getMemSize()
                             + rank_ptcl_domain_recv_.getMemSize()
                             + adr_ep_sd_send_.getMemSize()
                             + adr_sp_sd_send_.getMemSize()
                             + ep_recv_.getMemSize()
                             + sp_recv_.getMemSize()
                             + n_ep_domain_send_.getMemSize()
                             + n_sp_domain_send_.getMemSize()
                             + n_ep_sd_send_.getMemSize()
                             + n_sp_sd_send_.getMemSize()
                             + n_ep_domain_recv_.getMemSize()
                             + n_sp_domain_recv_.getMemSize()
                             + n_ep_sd_recv_.getMemSize()
                             + n_sp_sd_recv_.getMemSize()
                             + n_disp_ep_domain_send_.getMemSize()
                             + n_disp_sp_domain_send_.getMemSize()
                             + n_disp_ep_sd_send_.getMemSize()
                             + n_disp_sp_sd_send_.getMemSize()
                             + n_disp_ep_domain_recv_.getMemSize()
                             + n_disp_sp_domain_recv_.getMemSize()
                             + n_disp_ep_sd_recv_.getMemSize()
                             + n_disp_sp_sd_recv_.getMemSize()
                             + n_ep_domain_recv_2_.getMemSize()
                             + n_sp_domain_recv_2_.getMemSize()
                             + n_disp_ep_domain_recv_2_.getMemSize()
                             + n_disp_sp_domain_recv_2_.getMemSize()
                             + req_ep_send_.getMemSize()
                             + req_ep_recv_.getMemSize()
                             + req_sp_send_.getMemSize()
                             + req_sp_recv_.getMemSize()
                             + stat_ep_send_.getMemSize()
                             + stat_ep_recv_.getMemSize()
                             + stat_sp_send_.getMemSize()
                             + stat_sp_recv_.getMemSize()
                             + ep_send_.getMemSize()
                             + sp_send_.getMemSize()
                             + adr_ep_domain_send_.getMemSize()
                             + adr_sp_domain_send_.getMemSize())
                             * norm_fac;
        
            fout<<"*** CommTable class ***" << std::endl;
            fout<<"rank_sd_sd_send_         = " << rank_sd_sd_send_.getMemSize()*norm_fac        << str_unit << std::endl;
            fout<<"rank_sd_sd_send_         = " << rank_sd_sd_send_.getMemSize()*norm_fac        << str_unit << std::endl;
            fout<<"rank_sd_sd_recv_         = " << rank_sd_sd_recv_.getMemSize()*norm_fac        << str_unit << std::endl;
            fout<<"rank_ptcl_sd_send_       = " << rank_ptcl_sd_send_.getMemSize()*norm_fac      << str_unit << std::endl;
            fout<<"rank_ptcl_sd_recv_       = " << rank_ptcl_sd_recv_.getMemSize()*norm_fac      << str_unit << std::endl;
            fout<<"rank_ptcl_domain_send_   = " << rank_ptcl_domain_send_.getMemSize()*norm_fac  << str_unit << std::endl;
            fout<<"rank_ptcl_domain_recv_   = " << rank_ptcl_domain_recv_.getMemSize()*norm_fac  << str_unit << std::endl;

            fout<<"adr_ep_sd_send_          = " << adr_ep_sd_send_.getMemSize()*norm_fac         << str_unit << std::endl;
            fout<<"adr_sp_sd_send_          = " << adr_sp_sd_send_.getMemSize()*norm_fac         << str_unit << std::endl;

            fout<<"ep_recv_                 = " << ep_recv_.getMemSize()*norm_fac                << str_unit << std::endl;
            fout<<"sp_recv_                 = " << sp_recv_.getMemSize()*norm_fac                << str_unit << std::endl;

            fout<<"n_ep_domain_send_        = " << n_ep_domain_send_.getMemSize()*norm_fac       << str_unit << std::endl;
            fout<<"n_sp_domain_send_        = " << n_sp_domain_send_.getMemSize()*norm_fac       << str_unit << std::endl;
            fout<<"n_ep_sd_send_            = " << n_ep_sd_send_.getMemSize()*norm_fac           << str_unit << std::endl;
            fout<<"n_sp_sd_send_            = " << n_sp_sd_send_.getMemSize()*norm_fac           << str_unit << std::endl;

            fout<<"n_ep_domain_recv_        = " << n_ep_domain_recv_.getMemSize()*norm_fac       << str_unit << std::endl;
            fout<<"n_sp_domain_recv_        = " << n_sp_domain_recv_.getMemSize()*norm_fac       << str_unit << std::endl;
            fout<<"n_ep_sd_recv_            = " << n_ep_sd_recv_.getMemSize()*norm_fac           << str_unit << std::endl;
            fout<<"n_sp_sd_recv_            = " << n_sp_sd_recv_.getMemSize()*norm_fac           << str_unit << std::endl;

            fout<<"n_disp_ep_domain_send_   = " << n_disp_ep_domain_send_.getMemSize()*norm_fac  << str_unit << std::endl;
            fout<<"n_disp_sp_domain_send_   = " << n_disp_sp_domain_send_.getMemSize()*norm_fac  << str_unit << std::endl;
            fout<<"n_disp_ep_sd_send_       = " << n_disp_ep_sd_send_.getMemSize()*norm_fac      << str_unit << std::endl;
            fout<<"n_disp_sp_sd_send_       = " << n_disp_sp_sd_send_.getMemSize()*norm_fac      << str_unit << std::endl;

            fout<<"n_disp_ep_domain_recv_   = " << n_disp_ep_domain_recv_.getMemSize()*norm_fac  << str_unit << std::endl;
            fout<<"n_disp_sp_domain_recv_   = " << n_disp_sp_domain_recv_.getMemSize()*norm_fac  << str_unit << std::endl;
            fout<<"n_disp_ep_sd_recv_       = " << n_disp_ep_sd_recv_.getMemSize()*norm_fac      << str_unit << std::endl;
            fout<<"n_disp_sp_sd_recv_       = " << n_disp_sp_sd_recv_.getMemSize()*norm_fac      << str_unit << std::endl;

            fout<<"n_ep_domain_recv_2_      = "<< n_ep_domain_recv_2_.getMemSize()*norm_fac      << str_unit << std::endl;
            fout<<"n_sp_domain_recv_2_      = "<< n_sp_domain_recv_2_.getMemSize()*norm_fac      << str_unit << std::endl;
            fout<<"n_disp_ep_domain_recv_2_ = "<< n_disp_ep_domain_recv_2_.getMemSize()*norm_fac << str_unit << std::endl;
            fout<<"n_disp_sp_domain_recv_2_ = "<< n_disp_sp_domain_recv_2_.getMemSize()*norm_fac << str_unit << std::endl;

            fout<<"req_ep_send_             = "<< req_ep_send_.getMemSize()*norm_fac             << str_unit << std::endl;
            fout<<"req_ep_recv_             = "<< req_ep_recv_.getMemSize()*norm_fac             << str_unit << std::endl;
            fout<<"req_sp_send_             = "<< req_sp_send_.getMemSize()*norm_fac             << str_unit << std::endl;
            fout<<"req_sp_recv_             = "<< req_sp_recv_.getMemSize()*norm_fac             << str_unit << std::endl;
            fout<<"stat_ep_send_            = "<< stat_ep_send_.getMemSize()*norm_fac            << str_unit << std::endl;
            fout<<"stat_ep_recv_            = "<< stat_ep_recv_.getMemSize()*norm_fac            << str_unit << std::endl;
            fout<<"stat_sp_send_            = "<< stat_sp_send_.getMemSize()*norm_fac            << str_unit << std::endl;
            fout<<"stat_sp_recv_            = "<< stat_sp_recv_.getMemSize()*norm_fac            << str_unit << std::endl;

            fout<<"ep_send_                 = "<< ep_send_.getMemSize()*norm_fac                 << str_unit << std::endl;
            fout<<"sp_send_                 = "<< sp_send_.getMemSize()*norm_fac                 << str_unit << std::endl;

            fout<<"adr_ep_domain_send_      = "<< adr_ep_domain_send_.getMemSize()*norm_fac      << str_unit << std::endl;
            fout<<"adr_sp_domain_send_      = "<< adr_sp_domain_send_.getMemSize()*norm_fac      << str_unit << std::endl;

            fout<<"sum                      = "<< sum << str_unit << std::endl;
        }
    };
    
#else // USE_SUPER_DOMAIN

    template<class Tepj, class Tspj>
    class CommTable{
    public:
        S32 * n_ep_send_;
        S32 * n_sp_send_;
        S32 * n_ep_recv_;
        S32 * n_sp_recv_;
        S32 * n_disp_ep_send_;
        S32 * n_disp_sp_send_;
        S32 * n_disp_ep_recv_;
        S32 * n_disp_sp_recv_;
        ReallocatableArray<Tepj> ep_send_;
        ReallocatableArray<Tspj> sp_send_;
        ReallocatableArray<Tepj> ep_recv_;
        ReallocatableArray<Tspj> sp_recv_;
        S32 * n_ep_sp_send_;
        S32 * n_ep_sp_recv_;
        ReallocatableArray<S32>  rank_sp_isend_;
        ReallocatableArray<S32>  rank_sp_irecv_;
        ReallocatableArray<S32>  rank_sp_top_recv_;
        ReallocatableArray<S32>  rank_ep_send_;
        ReallocatableArray<S32>  rank_ep_recv_;
        Tspj * sp_top_recv_;
        Tspj   sp_top_;
        ReallocatableArray<MPI::Request> req_ep_send_;
        ReallocatableArray<MPI::Request> req_ep_recv_;
        ReallocatableArray<MPI::Request> req_sp_send_;
        ReallocatableArray<MPI::Request> req_sp_recv_;
        CountT n_proc_ep_send_;
        CountT n_proc_ep_recv_;
        CountT n_proc_sp_isend_;
        CountT n_proc_sp_irecv_;
        ReallocatableArray<F64ort> tree_outer_pos_;
        
    public:
        void initialize(){
            const S32 n_proc = Comm::getNumberOfProc();
            n_ep_send_ = new S32[n_proc];
            n_sp_send_ = new S32[n_proc];
            n_ep_recv_ = new S32[n_proc];
            n_sp_recv_ = new S32[n_proc];
            n_disp_ep_send_ = new S32[n_proc+1];
            n_disp_sp_send_ = new S32[n_proc+1];
            n_disp_ep_recv_ = new S32[n_proc+1];
            n_disp_sp_recv_ = new S32[n_proc+1];
            n_ep_sp_send_   = new S32[n_proc*2];
            n_ep_sp_recv_   = new S32[n_proc*2];
            rank_sp_isend_.reserve(n_proc);
            rank_sp_irecv_.reserve(n_proc);
            rank_sp_top_recv_.reserve(n_proc);
            rank_ep_send_.reserve(n_proc);
            rank_ep_recv_.reserve(n_proc);
            sp_top_recv_ = new Tspj[n_proc];
            tree_outer_pos_.reserve(n_proc);
            req_ep_send_.reserve(n_proc);
            req_sp_send_.reserve(n_proc);
            req_ep_recv_.reserve(n_proc);
            req_sp_recv_.reserve(n_proc);
        }

        void clearSize(){
            rank_sp_isend_.clearSize();
            rank_sp_irecv_.clearSize();
            ep_send_.clearSize();
            sp_send_.clearSize();
            ep_recv_.clearSize();
            sp_recv_.clearSize();
        }

        void setDispSend(){
            const S32 n_proc = Comm::getNumberOfProc();
            n_disp_ep_send_[0] = n_disp_sp_send_[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_ep_send_[i+1] = n_disp_ep_send_[i] + n_ep_send_[i];
                n_disp_sp_send_[i+1] = n_disp_sp_send_[i] + n_sp_send_[i];
            }
        }

        void setDispRecv(){
            const S32 n_proc = Comm::getNumberOfProc();
            n_disp_ep_recv_[0] = n_disp_sp_recv_[0] = 0;
            for(S32 i=0; i<n_proc; i++){
                n_disp_ep_recv_[i+1] = n_disp_ep_recv_[i] + n_ep_recv_[i];
                n_disp_sp_recv_[i+1] = n_disp_sp_recv_[i] + n_sp_recv_[i];
            }
        }

        void setSizeSendBuf(){
            const S32 n_proc = Comm::getNumberOfProc();
            ep_send_.resizeNoInitialize( n_disp_ep_send_[n_proc] );
            sp_send_.resizeNoInitialize( n_disp_sp_send_[n_proc] );
        }

        void setSizeRecvBuf(){
            const S32 n_proc = Comm::getNumberOfProc();
            ep_recv_.resizeNoInitialize( n_disp_ep_recv_[n_proc] );
            sp_recv_.resizeNoInitialize( n_disp_sp_recv_[n_proc] );
        }

        void exchangeLet(const bool reuse = false){
            if(!reuse){
                const S32 n_proc = Comm::getNumberOfProc();
                for(S32 i=0; i<n_proc; i++){
                    n_ep_sp_send_[2*i] = n_ep_send_[i];
                    n_ep_sp_send_[2*i+1] = n_sp_send_[i];
                }
                Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_);
                for(S32 i=0; i<n_proc; i++){
                    n_ep_recv_[i] = n_ep_sp_recv_[2*i];
                    n_sp_recv_[i] = n_ep_sp_recv_[2*i+1];
                }
                setDispRecv();
                setSizeRecvBuf();
            }
            Comm::allToAllV(ep_send_.getPointer(), n_ep_send_, n_disp_ep_send_,
                            ep_recv_.getPointer(), n_ep_recv_, n_disp_ep_recv_);
            Comm::allToAllV(sp_send_.getPointer(), n_sp_send_, n_disp_sp_send_,
                            sp_recv_.getPointer(), n_sp_recv_, n_disp_sp_recv_);
        }

        void setSPTop(const Tspj & _sp_top){
            sp_top_ = _sp_top;
        }

        void allgahterTreeOuterPos(const F64ort & _tree_outer_pos){
            Comm::allGather(&_tree_outer_pos, 1, tree_outer_pos_.getPointer());
        }


#if 1
        void exchangeLetSunWay(const ReallocatableArray<F64ort> & _tree_outer_pos,
                               const F64 _theta,
                               const bool reuse = false,
                               const DomainInfo * dinfo = NULL){
            const S32 n_proc = Comm::getNumberOfProc();
            const S32 my_rank = Comm::getRank();

            if(!reuse){
                // new line
                //allgahterTreeOuterPos(_tree_outer_pos);
                rank_sp_isend_.clearSize();
                rank_sp_irecv_.clearSize();
                rank_ep_send_.clearSize();
                rank_ep_recv_.clearSize();
                rank_sp_top_recv_.clearSize();
                S32 n_cnt_send = 0;
                S32 n_cnt_recv = 0;
                //Comm::barrier();
                
    #ifdef PHI_R_TREE
                const F64 length_peri_x = dinfo->getPosRootDomain().getFullLength()[0];
                //if(my_rank == 0) std::cerr<<"length_peri_x= "<<length_peri_x<<std::endl;
    #endif //PHI_R_TREE
                
    #ifdef SUNWAY_DEVELOP
    //#if 1
                const F64 longest_len_of_outer_box_send = _tree_outer_pos[my_rank].getFullLength().getMax();
                const F64 r_crit_sq_send = (longest_len_of_outer_box_send * longest_len_of_outer_box_send) / (_theta * _theta);

                F64 wtime_offset_0 = GetWtime();
                for(S32 i=0; i<n_proc; i++){
                    if(my_rank == i) continue;
                    const F64 longest_len_of_outer_box_recv = _tree_outer_pos[i].getFullLength().getMax();
                    const F64 r_crit_sq_recv = (longest_len_of_outer_box_recv * longest_len_of_outer_box_recv) / (_theta * _theta);
        #ifdef PHI_R_TREE
                    if(GetDistanceMinSqPeriodicX(_tree_outer_pos[i], _tree_outer_pos[my_rank], length_peri_x) <= r_crit_sq_send){
        #else // PHI_R_TREE
                    if(_tree_outer_pos[i].getDistanceMinSQ(_tree_outer_pos[my_rank]) <= r_crit_sq_send){
        #endif // PHI_R_TREE
                        // to be consistent with make LET function in tree_for_force_impl3.hpp
                        //rank_sp_isend_.pushBackNoCheck(i);
                        //rank_ep_send_.pushBackNoCheck(i);
                        req_sp_send_[n_cnt_send] = MPI::COMM_WORLD.Isend(n_sp_send_+i, 1, GetDataType<S32>(), i, 0);
                        req_ep_send_[n_cnt_send] = MPI::COMM_WORLD.Isend(n_ep_send_+i, 1, GetDataType<S32>(), i, 1);
                        n_cnt_send++;
                    }
                    else{
                        assert(n_sp_send_[i]==1);
                        assert(n_ep_send_[i]==0);
                    }
        #ifdef PHI_R_TREE
                    //if(GetDistanceMinSqPeriodicX(_tree_outer_pos[i], _tree_outer_pos[my_rank], length_peri_x) <= r_crit_sq_send){
                    if(GetDistanceMinSqPeriodicX(_tree_outer_pos[my_rank], _tree_outer_pos[i], length_peri_x) <= r_crit_sq_recv){
        #else  // PHI_R_TREE
                    if(_tree_outer_pos[my_rank].getDistanceMinSQ(_tree_outer_pos[i]) <= r_crit_sq_recv){
        #endif // PHI_R_TREE
                        req_sp_recv_[n_cnt_recv] = MPI::COMM_WORLD.Irecv(n_sp_recv_+i, 1, GetDataType<S32>(), i, 0);
                        req_ep_recv_[n_cnt_recv] = MPI::COMM_WORLD.Irecv(n_ep_recv_+i, 1, GetDataType<S32>(), i, 1);
                        n_cnt_recv++;
                    }
                    else{
                        n_sp_recv_[i] = 1;
                        n_ep_recv_[i] = 0;
                    }
                }
                MPI::Request::Waitall(n_cnt_send, req_ep_send_.getPointer());
                MPI::Request::Waitall(n_cnt_send, req_sp_send_.getPointer());
                MPI::Request::Waitall(n_cnt_recv, req_ep_recv_.getPointer());
                MPI::Request::Waitall(n_cnt_recv, req_sp_recv_.getPointer());
                wtime_ex_let_exchange_number += GetWtime() - wtime_offset_0;
                //std::cerr<<"my_rank= "<<my_rank<<std::endl;
                /*
                for(S32 i=0; i<n_proc; i++){
                    if(n_sp_send_[i] > 1){
                        rank_sp_isend_.pushBackNoCheck(i);
                    }
                    if(n_ep_send_[i] > 0) rank_ep_send_.pushBackNoCheck(i);
                }
                for(S32 i=0; i<n_proc; i++){
                    if(n_sp_recv_[i] > 1) rank_sp_irecv_.pushBackNoCheck(i);
                    if(n_ep_recv_[i] > 0) rank_ep_recv_.pushBackNoCheck(i);
                }
                */
                for(S32 i=0; i<n_proc; i++){
                    if(n_sp_send_[i] > 1 || (n_sp_send_[i] == 1 && n_ep_send_[i] > 0) ){
                        rank_sp_isend_.pushBackNoCheck(i);
                    }
                    if(n_ep_send_[i] > 0) rank_ep_send_.pushBackNoCheck(i);
                }
                for(S32 i=0; i<n_proc; i++){
                    if(n_sp_recv_[i] > 1 || (n_sp_recv_[i] == 1 && n_ep_recv_[i] > 0) ){
                        rank_sp_irecv_.pushBackNoCheck(i);
                    }
                    if(n_ep_recv_[i] > 0) rank_ep_recv_.pushBackNoCheck(i);
                }
                
                // debug
        #ifdef DEBUG_EXLET_SUNWAY
                ReallocatableArray<S32> rank_sp_isend_tmp;
                ReallocatableArray<S32> rank_ep_send_tmp;
                ReallocatableArray<S32> rank_sp_irecv_tmp;
                ReallocatableArray<S32> rank_ep_recv_tmp;
                for(S32 i=0; i<n_proc; i++){
                    n_ep_sp_send_[2*i]   = n_ep_send_[i];
                    n_ep_sp_send_[2*i+1] = n_sp_send_[i];
                    //if(n_sp_send_[i] > 1){
                    if(n_sp_send_[i] > 1 || (n_sp_send_[i] == 1 && n_ep_send_[i] > 0) ){
                        rank_sp_isend_tmp.push_back(i);
                    }
                    if(n_ep_send_[i] > 0) rank_ep_send_tmp.push_back(i);
                }
                Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_);
                for(S32 i=0; i<n_proc; i++){
                    if(n_ep_sp_recv_[2*i+1] > 1 || (n_ep_sp_recv_[2*i+1] == 1 && n_ep_sp_recv_[2*i] > 0) ){
                        rank_sp_irecv_tmp.push_back(i);
                    }
                    //if(n_ep_sp_recv_[2*i+1] > 1) rank_sp_irecv_tmp.push_back(i);
                    if(n_ep_sp_recv_[2*i] > 0) rank_ep_recv_tmp.push_back(i);
                }
                assert(rank_sp_irecv_tmp.size() == rank_sp_irecv_.size());
                assert(rank_ep_recv_tmp.size()  == rank_ep_recv_.size());
                assert(rank_sp_isend_tmp.size() == rank_sp_isend_.size());
                assert(rank_ep_send_tmp.size()  == rank_ep_send_.size());
        #endif
                // debug
    #else
                for(S32 i=0; i<n_proc; i++){
                    n_ep_sp_send_[2*i]   = n_ep_send_[i];
                    n_ep_sp_send_[2*i+1] = n_sp_send_[i];
                    if(n_sp_send_[i] > 1){
                        rank_sp_isend_.pushBackNoCheck(i);
                    }
                    if(n_ep_send_[i] > 0) rank_ep_send_.pushBackNoCheck(i);
                }
                Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_);
                for(S32 i=0; i<n_proc; i++){
                    n_ep_recv_[i] = n_ep_sp_recv_[2*i];
                    n_sp_recv_[i] = n_ep_sp_recv_[2*i+1];
                    if(n_sp_recv_[i] > 1) rank_sp_irecv_.pushBackNoCheck(i);
                    if(n_ep_recv_[i] > 0) rank_ep_recv_.pushBackNoCheck(i);
                }
    #endif
                setDispRecv();
                setSizeRecvBuf();
            }
            /*
            Comm::allGather(&sp_top_, 1, sp_top_recv_);
            const S32 n_recv_sp_top = rank_sp_top_recv_.size();
            for(S32 i=0; i<n_recv_sp_top; i++){
                const S32 target_rank = rank_sp_top_recv_[i];
                sp_recv_[n_disp_sp_recv_[target_rank]] = sp_top_recv_[target_rank];
            }
            */
            /////////
            // set EP
            n_ep_send_tot_loc = n_disp_ep_send_[n_proc];
            F64 wtime_offset_ep = GetWtime() ;
    #ifndef SUNWAY_DEVELOP
            //#if 0
            Comm::allToAllV(ep_send_.getPointer(), n_ep_send_, n_disp_ep_send_,
                            ep_recv_.getPointer(), n_ep_recv_, n_disp_ep_recv_);

    #else
            const S32 n_proc_ep_send = rank_ep_send_.size();
            const S32 n_proc_ep_recv = rank_ep_recv_.size();
            n_proc_ep_send_ = n_proc_ep_send;
            n_proc_ep_recv_ = n_proc_ep_recv;
            req_ep_send_.resizeNoInitialize(n_proc_ep_send);
            req_ep_recv_.resizeNoInitialize(n_proc_ep_recv);
            //std::cerr<<"my_rank(D)= "<<my_rank<<std::endl;
            for(S32 i=0; i<n_proc_ep_send; i++){
                const S32 rank_send   = rank_ep_send_[i];
                const S32 n_ep_send      = n_ep_send_[rank_send];
                const S32 n_disp_ep_send = n_disp_ep_send_[rank_send];
                const S32 tag = 0;
                Tepj * ep_send = ep_send_.getPointer(n_disp_ep_send);
                req_ep_send_[i] = MPI::COMM_WORLD.Isend
                    (ep_send, n_ep_send, GetDataType<Tepj>(), rank_send, tag);
            }
            for(S32 i=0; i<n_proc_ep_recv; i++){
                const S32 rank_recv   = rank_ep_recv_[i];
                const S32 n_ep_recv      = n_ep_recv_[rank_recv];
                const S32 n_disp_ep_recv = n_disp_ep_recv_[rank_recv];
                const S32 tag = 0;
                Tepj * ep_recv        = ep_recv_.getPointer(n_disp_ep_recv);
                req_ep_recv_[i] = MPI::COMM_WORLD.Irecv
                    (ep_recv, n_ep_recv, GetDataType<Tepj>(), rank_recv, tag);
            }
            MPI::Request::Waitall(n_proc_ep_send, req_ep_send_.getPointer());
            MPI::Request::Waitall(n_proc_ep_recv, req_ep_recv_.getPointer());
    #endif
            wtime_alltoallv_ep = GetWtime() - wtime_offset_ep;
            //band_width_ep = n_ep_send_tot_loc*sizeof(Tepj) / wtime_alltoallv_ep;
            // set EP
            /////////

            /////////
            // set SP
            F64 wtime_offset_sp = GetWtime();
    #ifndef SUNWAY_DEVELOP
            //#if 1
            Comm::allToAllV(sp_send_.getPointer(), n_sp_send_, n_disp_sp_send_,
                            sp_recv_.getPointer(), n_sp_recv_, n_disp_sp_recv_);
    #else //SUNWAY_DEVELOP
            Comm::allGather(&sp_top_, 1, sp_top_recv_);
            wtime_alltoallv_sp_allgather += GetWtime() - wtime_offset_sp;
            
            const S32 n_proc_sp_send = rank_sp_isend_.size();
            const S32 n_proc_sp_recv = rank_sp_irecv_.size();

            n_proc_sp_isend_ = n_proc_sp_send;
            n_proc_sp_irecv_ = n_proc_sp_recv;
            req_sp_send_.resizeNoInitialize(n_proc_sp_send);
            req_sp_recv_.resizeNoInitialize(n_proc_sp_recv);
            n_sp_send_tot_loc = n_disp_sp_send_[n_proc];
            F64 wtime_offset_isendrecv = GetWtime();

            for(S32 i=0; i<n_proc_sp_send; i++){
                const S32 rank_send   = rank_sp_isend_[i];
                const S32 n_sp_send      = n_sp_send_[rank_send];
                n_sp_send_tot_isendrecv_loc += n_sp_send;
                const S32 n_disp_sp_send = n_disp_sp_send_[rank_send];
                const S32 tag = 1;
                Tspj * sp_send = sp_send_.getPointer(n_disp_sp_send);
                req_sp_send_[i] = MPI::COMM_WORLD.Isend
                    (sp_send, n_sp_send, GetDataType<Tspj>(), rank_send, tag);
            }
            for(S32 i=0; i<n_proc_sp_recv; i++){
                const S32 rank_recv   = rank_sp_irecv_[i];
                const S32 n_sp_recv      = n_sp_recv_[rank_recv];
                const S32 n_disp_sp_recv = n_disp_sp_recv_[rank_recv];
                const S32 tag = 1;
                Tspj * sp_recv        = sp_recv_.getPointer(n_disp_sp_recv);
                req_sp_recv_[i] = MPI::COMM_WORLD.Irecv
                    (sp_recv, n_sp_recv, GetDataType<Tspj>(), rank_recv, tag);
            }
            MPI::Request::Waitall(n_proc_sp_send, req_sp_send_.getPointer());
            MPI::Request::Waitall(n_proc_sp_recv, req_sp_recv_.getPointer());

            wtime_alltoallv_sp_isendrecv += GetWtime() - wtime_offset_isendrecv;
            //band_width_sp = n_sp_send_tot_loc*sizeof(Tepj) / wtime_alltoallv_ep;
            
            for(S32 i=0; i<n_proc; i++){
                if(n_sp_recv_[i] == 1 && n_ep_recv_[i] == 0){ // bug fix 
                    sp_recv_[n_disp_sp_recv_[i]] = sp_top_recv_[i];
                }
            }

            // for debug
         #ifdef DEBUG_EXLET_SUNWAY
            ReallocatableArray<Tspj> sp_recv_tmp;
            sp_recv_tmp.resizeNoInitialize(n_disp_sp_recv_[n_proc]);
            Comm::allToAllV(sp_send_.getPointer(), n_sp_send_, n_disp_sp_send_,
                            sp_recv_tmp.getPointer(), n_sp_recv_, n_disp_sp_recv_);
            for(S32 i=0; i<n_disp_sp_recv_[n_proc]; i++){
                if(sp_recv_tmp[i].mass != sp_recv_[i].mass){
                    std::cerr<<"my_rank= "<<my_rank
                             <<" i= "<<i
                             <<" n_disp_sp_recv_[n_proc]= "<<n_disp_sp_recv_[n_proc]
                             <<std::endl;
                    std::cerr<<" sp_recv_tmp[i].mass= "<<sp_recv_tmp[i].mass
                             <<" sp_recv_[i].mass= "<<sp_recv_[i].mass
                             <<" sp_recv_tmp[i].pos= "<<sp_recv_tmp[i].pos
                             <<" sp_recv_[i].pos= "<<sp_recv_[i].pos
                             <<std::endl;
                    std::cerr<<" sp_recv_tmp[i+1].mass= "<<sp_recv_tmp[i+1].mass
                             <<" sp_recv_[i+1].mass= "<<sp_recv_[i+1].mass
                             <<" sp_recv_tmp[i+1].pos= "<<sp_recv_tmp[i+1].pos
                             <<" sp_recv_[i+1].pos= "<<sp_recv_[i+1].pos
                             <<std::endl;
                    std::cerr<<" sp_recv_tmp[i-1].mass= "<<sp_recv_tmp[i-1].mass
                             <<" sp_recv_[i-1].mass= "<<sp_recv_[i-1].mass
                             <<" sp_recv_tmp[i-1].pos= "<<sp_recv_tmp[i-1].pos
                             <<" sp_recv_[i-1].pos= "<<sp_recv_[i-1].pos
                             <<std::endl;
                    for(S32 j=0; j<n_proc; j++){
                        std::cerr<<"j= "<<j
                                 <<" n_sp_recv_[j]= "<<n_sp_recv_[j]
                                 <<" n_ep_recv_[j]= "<<n_ep_recv_[j]
                                 <<std::endl;
                    }
                }
                assert(sp_recv_tmp[i].mass == sp_recv_[i].mass);
            }
         #endif

    #endif //SUNWAY_DEVELOP
            wtime_alltoallv_sp += GetWtime() - wtime_offset_sp;
            // set SP
            /////////
            //std::cerr<<"my_rank(F)= "<<my_rank<<std::endl;
        }
                
#elif 0
        void exchangeLetSunWay(const ReallocatableArray<F64ort> & _tree_outer_pos,
                               const F64 _theta,
                               const bool reuse = false){
            const S32 n_proc = Comm::getNumberOfProc();
            const S32 my_rank = Comm::getRank();
            if(!reuse){
                // new line
                //allgahterTreeOuterPos(_tree_outer_pos);
                rank_sp_isend_.clearSize();
                rank_sp_irecv_.clearSize();
                rank_ep_send_.clearSize();
                rank_ep_recv_.clearSize();
                rank_sp_top_recv_.clearSize();
                S32 n_cnt_send = 0;
                S32 n_cnt_recv = 0;
                Comm::barrier();
                
    #ifdef SUNWAY_DEVELOP
                const F64 longest_len_of_outer_box_send = _tree_outer_pos[my_rank].getFullLength().getMax();
                const F64 r_crit_sq_send = (longest_len_of_outer_box_send * longest_len_of_outer_box_send) / (_theta * _theta);
                for(S32 i=0; i<n_proc; i++){
                    if(my_rank == i) continue;
                    const F64 longest_len_of_outer_box_recv = _tree_outer_pos[i].getFullLength().getMax();
                    const F64 r_crit_sq_recv = (longest_len_of_outer_box_recv * longest_len_of_outer_box_recv) / (_theta * _theta);
                    if(_tree_outer_pos[i].getDistanceMinSQ(_tree_outer_pos[my_rank]) <= r_crit_sq_send){
                        // to be consistent with make LET function in tree_for_force_impl3.hpp
                        rank_sp_isend_.pushBackNoCheck(i);
                        rank_ep_send_.pushBackNoCheck(i);
                        req_sp_send_[n_cnt_send] = MPI::COMM_WORLD.Isend(n_sp_send_+i, 1, GetDataType<S32>(), i, 0);
                        req_ep_send_[n_cnt_send] = MPI::COMM_WORLD.Isend(n_ep_send_+i, 1, GetDataType<S32>(), i, 1);
                        n_cnt_send++;
                    }
                    else{
                        assert(n_sp_send_[i]==1);
                        assert(n_ep_send_[i]==0);
                    }
                    if(_tree_outer_pos[my_rank].getDistanceMinSQ(_tree_outer_pos[i]) <= r_crit_sq_recv){
                        rank_sp_irecv_.pushBackNoCheck(i);
                        rank_ep_recv_.pushBackNoCheck(i);
                        req_sp_recv_[n_cnt_recv] = MPI::COMM_WORLD.Irecv(n_sp_recv_+i, 1, GetDataType<S32>(), i, 0);
                        req_ep_recv_[n_cnt_recv] = MPI::COMM_WORLD.Irecv(n_ep_recv_+i, 1, GetDataType<S32>(), i, 1);
                        n_cnt_recv++;
                        //if(i==0) std::cout<<" my_rank= "<<my_rank<<std::endl;
                    }
                    else{
                        n_sp_recv_[i] = 1;
                        n_ep_recv_[i] = 0;
                        rank_sp_top_recv_.pushBackNoCheck(i);
                    }
                }
                MPI::Request::Waitall(n_cnt_send, req_ep_send_.getPointer());
                MPI::Request::Waitall(n_cnt_send, req_sp_send_.getPointer());
                MPI::Request::Waitall(n_cnt_recv, req_ep_recv_.getPointer());
                MPI::Request::Waitall(n_cnt_recv, req_sp_recv_.getPointer());


                // debug
                for(S32 i=0; i<n_proc; i++){
                    n_ep_sp_send_[2*i]   = n_ep_send_[i];
                    n_ep_sp_send_[2*i+1] = n_sp_send_[i];
                    if(n_sp_send_[i] > 1){
                        rank_sp_isend_.pushBackNoCheck(i);
                    }
                    if(n_ep_send_[i] > 0) rank_ep_send_.pushBackNoCheck(i);
                }
                Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_);
                for(S32 i=0; i<n_proc; i++){
                    assert(n_ep_recv_[i] == n_ep_sp_recv_[2*i]);
                    assert(n_sp_recv_[i] == n_ep_sp_recv_[2*i+1]);
                }
                //std::cerr<<"my_rank(D)= "<<my_rank<<std::endl;
                
                //std::cerr<<"my_rank(C)= "<<my_rank<<std::endl;
                /*
                if(my_rank==0){
                    for(S32 i=0; i<n_proc; i++){
                        std::cerr<<"i= "<<i
                                 <<" n_ep_recv_[i]= "<<n_ep_recv_[i]
                                 <<" n_disp_ep_recv_[i]= "<<n_disp_ep_recv_[i]
                                 <<" n_sp_recv_[i]= "<<n_sp_recv_[i]
                                 <<" n_disp_sp_recv_[i]= "<<n_disp_sp_recv_[i]
                                 <<std::endl;
                    }
                }
                */

    #else
                for(S32 i=0; i<n_proc; i++){
                    n_ep_sp_send_[2*i]   = n_ep_send_[i];
                    n_ep_sp_send_[2*i+1] = n_sp_send_[i];
                    if(n_sp_send_[i] > 1){
                        rank_sp_isend_.pushBackNoCheck(i);
                    }
                    if(n_ep_send_[i] > 0) rank_ep_send_.pushBackNoCheck(i);
                }
                Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_);
                for(S32 i=0; i<n_proc; i++){
                    n_ep_recv_[i] = n_ep_sp_recv_[2*i];
                    n_sp_recv_[i] = n_ep_sp_recv_[2*i+1];
                    if(n_sp_recv_[i] > 1) rank_sp_irecv_.pushBackNoCheck(i);
                    if(n_ep_recv_[i] > 0) rank_ep_recv_.pushBackNoCheck(i);
                }
    #endif
                setDispRecv();
                setSizeRecvBuf();
                
            }
            //std::cerr<<"CHECK C"<<std::endl;            
            /*
            Comm::allGather(&sp_top_, 1, sp_top_recv_);
            const S32 n_recv_sp_top = rank_sp_top_recv_.size();
            for(S32 i=0; i<n_recv_sp_top; i++){
                const S32 target_rank = rank_sp_top_recv_[i];
                sp_recv_[n_disp_sp_recv_[target_rank]] = sp_top_recv_[target_rank];
            }
            */
            /////////////////////////////////////////////////////////////
            ///////////////////////// BELOW TODO ////////////////////////
            /////////
            // set EP
            F64 wtime_offset_ep = GetWtime() ;
    #ifndef SUNWAY_DEVELOP
            Comm::allToAllV(ep_send_.getPointer(), n_ep_send_, n_disp_ep_send_,
                            ep_recv_.getPointer(), n_ep_recv_, n_disp_ep_recv_);

    #else
            const S32 n_proc_ep_send = rank_ep_send_.size();
            const S32 n_proc_ep_recv = rank_ep_recv_.size();
            n_proc_ep_send_ = n_proc_ep_send;
            n_proc_ep_recv_ = n_proc_ep_recv;
            req_ep_send_.resizeNoInitialize(n_proc_ep_send);
            req_ep_recv_.resizeNoInitialize(n_proc_ep_recv);
            for(S32 i=0; i<n_proc_ep_send; i++){
                const S32 rank_send   = rank_ep_send_[i];
                const S32 n_ep_send      = n_ep_send_[rank_send];
                const S32 n_disp_ep_send = n_disp_ep_send_[rank_send];
                const S32 tag = 0;
                Tepj * ep_send = ep_send_.getPointer(n_disp_ep_send);
                req_ep_send_[i] = MPI::COMM_WORLD.Isend
                    (ep_send, n_ep_send, GetDataType<Tepj>(), rank_send, tag);
            }
            for(S32 i=0; i<n_proc_ep_recv; i++){
                const S32 rank_recv   = rank_ep_recv_[i];
                const S32 n_ep_recv      = n_ep_recv_[rank_recv];
                const S32 n_disp_ep_recv = n_disp_ep_recv_[rank_recv];
                const S32 tag = 0;
                Tepj * ep_recv        = ep_recv_.getPointer(n_disp_ep_recv);
                req_ep_recv_[i] = MPI::COMM_WORLD.Irecv
                    (ep_recv, n_ep_recv, GetDataType<Tepj>(), rank_recv, tag);
            }
            MPI::Request::Waitall(n_proc_ep_send, req_ep_send_.getPointer());
            MPI::Request::Waitall(n_proc_ep_recv, req_ep_recv_.getPointer());
    #endif
            wtime_alltoallv_ep += GetWtime() - wtime_offset_ep;
            // set EP
            /////////
            
            //std::cerr<<"my_rank(E)= "<<my_rank<<std::endl;
            
            /////////
            // set SP
            F64 wtime_offset_sp = GetWtime() ;
    #ifndef SUNWAY_DEVELOP            
            Comm::allToAllV(sp_send_.getPointer(), n_sp_send_, n_disp_sp_send_,
                            sp_recv_.getPointer(), n_sp_recv_, n_disp_sp_recv_);
    #else
            Comm::allGather(&sp_top_, 1, sp_top_recv_);
            const S32 n_proc_sp_send = rank_sp_isend_.size();
            const S32 n_proc_sp_recv = rank_sp_irecv_.size();

            n_proc_sp_isend_ = n_proc_sp_send;
            n_proc_sp_irecv_ = n_proc_sp_recv;
            req_sp_send_.resizeNoInitialize(n_proc_sp_send);
            req_sp_recv_.resizeNoInitialize(n_proc_sp_recv);
            for(S32 i=0; i<n_proc_sp_send; i++){
                const S32 rank_send   = rank_sp_isend_[i];
                const S32 n_sp_send      = n_sp_send_[rank_send];
                const S32 n_disp_sp_send = n_disp_sp_send_[rank_send];
                const S32 tag = 1;
                Tspj * sp_send = sp_send_.getPointer(n_disp_sp_send);
                req_sp_send_[i] = MPI::COMM_WORLD.Isend
                    (sp_send, n_sp_send, GetDataType<Tspj>(), rank_send, tag);
            }
            for(S32 i=0; i<n_proc_sp_recv; i++){
                const S32 rank_recv   = rank_sp_irecv_[i];
                const S32 n_sp_recv      = n_sp_recv_[rank_recv];
                const S32 n_disp_sp_recv = n_disp_sp_recv_[rank_recv];
                const S32 tag = 1;
                Tspj * sp_recv        = sp_recv_.getPointer(n_disp_sp_recv);
                req_sp_recv_[i] = MPI::COMM_WORLD.Irecv
                    (sp_recv, n_sp_recv, GetDataType<Tspj>(), rank_recv, tag);
            }
            MPI::Request::Waitall(n_proc_sp_send, req_sp_send_.getPointer());
            MPI::Request::Waitall(n_proc_sp_recv, req_sp_recv_.getPointer());

            for(S32 i=0; i<n_proc; i++){
                //if(n_sp_recv_[i] == 1 && n_ep_recv_[i] == 0){
                if(n_sp_recv_[i] == 1){
                    sp_recv_[n_disp_sp_recv_[i]] = sp_top_recv_[i];
                }
            }
    #endif
            wtime_alltoallv_sp += GetWtime() - wtime_offset_sp;
            // set SP
            /////////
            //std::cerr<<"my_rank(F)= "<<my_rank<<std::endl;
        }

#else
        void exchangeLetSunWay(const F64ort & _tree_outer_pos,
                               const bool reuse = false){
            const S32 n_proc = Comm::getNumberOfProc();
            if(!reuse){
                rank_sp_isend_.clearSize();
                rank_sp_irecv_.clearSize();
                rank_ep_send_.clearSize();
                rank_ep_recv_.clearSize();
                for(S32 i=0; i<n_proc; i++){
                    n_ep_sp_send_[2*i]   = n_ep_send_[i];
                    n_ep_sp_send_[2*i+1] = n_sp_send_[i];
                    if(n_sp_send_[i] > 1){
                        rank_sp_isend_.pushBackNoCheck(i);
                    }
                    if(n_ep_send_[i] > 0) rank_ep_send_.pushBackNoCheck(i);
                }
                Comm::allToAll(n_ep_sp_send_, 2, n_ep_sp_recv_);
                for(S32 i=0; i<n_proc; i++){
                    n_ep_recv_[i] = n_ep_sp_recv_[2*i];
                    n_sp_recv_[i] = n_ep_sp_recv_[2*i+1];
                    if(n_sp_recv_[i] > 1) rank_sp_irecv_.pushBackNoCheck(i);
                    if(n_ep_recv_[i] > 0) rank_ep_recv_.pushBackNoCheck(i);
                }
                setDispRecv();
                setSizeRecvBuf();
                // new line
                allgahterTreeOuterPos(_tree_outer_pos);
            }

            /////////
            // set EP
            F64 wtime_offset_ep = GetWtime() ;
    #if 0
            Comm::allToAllV(ep_send_.getPointer(), n_ep_send_, n_disp_ep_send_,
                            ep_recv_.getPointer(), n_ep_recv_, n_disp_ep_recv_);

    #else
            const S32 n_proc_ep_send = rank_ep_send_.size();
            const S32 n_proc_ep_recv = rank_ep_recv_.size();
            n_proc_ep_send_ = n_proc_ep_send;
            n_proc_ep_recv_ = n_proc_ep_recv;
            req_ep_send_.resizeNoInitialize(n_proc_ep_send);
            req_ep_recv_.resizeNoInitialize(n_proc_ep_recv);
            for(S32 i=0; i<n_proc_ep_send; i++){
                const S32 rank_send   = rank_ep_send_[i];
                const S32 n_ep_send      = n_ep_send_[rank_send];
                const S32 n_disp_ep_send = n_disp_ep_send_[rank_send];
                const S32 tag = 0;
                Tepj * ep_send = ep_send_.getPointer(n_disp_ep_send);
                req_ep_send_[i] = MPI::COMM_WORLD.Isend
                    (ep_send, n_ep_send, GetDataType<Tepj>(), rank_send, tag);
            }
            for(S32 i=0; i<n_proc_ep_recv; i++){
                const S32 rank_recv   = rank_ep_recv_[i];
                const S32 n_ep_recv      = n_ep_recv_[rank_recv];
                const S32 n_disp_ep_recv = n_disp_ep_recv_[rank_recv];
                const S32 tag = 0;
                Tepj * ep_recv        = ep_recv_.getPointer(n_disp_ep_recv);
                req_ep_recv_[i] = MPI::COMM_WORLD.Irecv
                    (ep_recv, n_ep_recv, GetDataType<Tepj>(), rank_recv, tag);
            }
            MPI::Request::Waitall(n_proc_ep_send, req_ep_send_.getPointer());
            MPI::Request::Waitall(n_proc_ep_recv, req_ep_recv_.getPointer());            
    #endif
            wtime_alltoallv_ep += GetWtime() - wtime_offset_ep;
            // set EP
            /////////

            /////////
            // set SP
            F64 wtime_offset_sp = GetWtime() ;
    #if 0            
            Comm::allToAllV(sp_send_.getPointer(), n_sp_send_, n_disp_sp_send_,
                            sp_recv_.getPointer(), n_sp_recv_, n_disp_sp_recv_);
    #else
            Comm::allGather(&sp_top_, 1, sp_top_recv_);
            const S32 n_proc_sp_send = rank_sp_isend_.size();
            const S32 n_proc_sp_recv = rank_sp_irecv_.size();

            n_proc_sp_isend_ = n_proc_sp_send;
            n_proc_sp_irecv_ = n_proc_sp_recv;
            req_sp_send_.resizeNoInitialize(n_proc_sp_send);
            req_sp_recv_.resizeNoInitialize(n_proc_sp_recv);
            for(S32 i=0; i<n_proc_sp_send; i++){
                const S32 rank_send   = rank_sp_isend_[i];
                const S32 n_sp_send      = n_sp_send_[rank_send];
                const S32 n_disp_sp_send = n_disp_sp_send_[rank_send];
                const S32 tag = 1;
                Tspj * sp_send = sp_send_.getPointer(n_disp_sp_send);
                req_sp_send_[i] = MPI::COMM_WORLD.Isend
                    (sp_send, n_sp_send, GetDataType<Tspj>(), rank_send, tag);
            }
            for(S32 i=0; i<n_proc_sp_recv; i++){
                const S32 rank_recv   = rank_sp_irecv_[i];
                const S32 n_sp_recv      = n_sp_recv_[rank_recv];
                const S32 n_disp_sp_recv = n_disp_sp_recv_[rank_recv];
                const S32 tag = 1;
                Tspj * sp_recv        = sp_recv_.getPointer(n_disp_sp_recv);
                req_sp_recv_[i] = MPI::COMM_WORLD.Irecv
                    (sp_recv, n_sp_recv, GetDataType<Tspj>(), rank_recv, tag);
            }
            MPI::Request::Waitall(n_proc_sp_send, req_sp_send_.getPointer());
            MPI::Request::Waitall(n_proc_sp_recv, req_sp_recv_.getPointer());

            for(S32 i=0; i<n_proc; i++){
                //if(n_sp_recv_[i] == 1 && n_ep_recv_[i] == 0){
                if(n_sp_recv_[i] == 1){
                    sp_recv_[n_disp_sp_recv_[i]] = sp_top_recv_[i];
                }
            }
    #endif
            wtime_alltoallv_sp += GetWtime() - wtime_offset_sp;
            // set SP
            /////////
        }
#endif
        /*
        void dumpMemSizeUsed(std::ostream & fout){
            fout<<"ep_send_.getMemSize()= "<<ep_send_.getMemSize()<<std::endl;
            fout<<"sp_send_.getMemSize()= "<<sp_send_.getMemSize()<<std::endl;
            fout<<"ep_recv_.getMemSize()= "<<ep_recv_.getMemSize()<<std::endl;
            fout<<"sp_recv_.getMemSize()= "<<sp_recv_.getMemSize()<<std::endl;
            fout<<"rank_sp_isend_.getMemSize()= "<<rank_sp_isend_.getMemSize()<<std::endl;
            fout<<"rank_sp_irecv_.getMemSize()= "<<rank_sp_irecv_.getMemSize()<<std::endl;
            fout<<"rank_sp_top_recv_.getMemSize()= "<<rank_sp_top_recv_.getMemSize()<<std::endl;
            fout<<"rank_ep_send_.getMemSize()= "<<rank_ep_send_.getMemSize()<<std::endl;
            fout<<"rank_ep_recv_.getMemSize()= "<<rank_ep_recv_.getMemSize()<<std::endl;
            fout<<"sum= "<<
                (double)(ep_send_.getMemSize()
                         +sp_send_.getMemSize()
                         +ep_recv_.getMemSize()
                         +sp_recv_.getMemSize()
                         +rank_sp_isend_.getMemSize()
                         +rank_sp_irecv_.getMemSize()
                         +rank_sp_top_recv_.getMemSize()
                         +rank_ep_send_.getMemSize()
                         +rank_ep_recv_.getMemSize()) / 1e9
                <<" [GB]"
                <<std::endl;
        }
        */
    };
#endif
}


