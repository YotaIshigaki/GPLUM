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
                Comm::barrier();
                
#ifdef SUNWAY_DEVELOP
//#if 1
                const F64 longest_len_of_outer_box_send = _tree_outer_pos[my_rank].getFullLength().getMax();
                const F64 r_crit_sq_send = (longest_len_of_outer_box_send * longest_len_of_outer_box_send) / (_theta * _theta);
                for(S32 i=0; i<n_proc; i++){
                    if(my_rank == i) continue;
                    const F64 longest_len_of_outer_box_recv = _tree_outer_pos[i].getFullLength().getMax();
                    const F64 r_crit_sq_recv = (longest_len_of_outer_box_recv * longest_len_of_outer_box_recv) / (_theta * _theta);
                    if(_tree_outer_pos[i].getDistanceMinSQ(_tree_outer_pos[my_rank]) <= r_crit_sq_send){
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
                    if(_tree_outer_pos[my_rank].getDistanceMinSQ(_tree_outer_pos[i]) <= r_crit_sq_recv){
                        //rank_sp_irecv_.pushBackNoCheck(i);
                        //rank_ep_recv_.pushBackNoCheck(i);
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
                /*
                // debug
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
            //std::cerr<<"my_rank(C)= "<<my_rank<<std::endl;
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

            //std::cerr<<"my_rank(E)= "<<my_rank<<std::endl;
            
            /////////
            // set SP
            F64 wtime_offset_sp = GetWtime();
#ifndef SUNWAY_DEVELOP
            //#if 1
            Comm::allToAllV(sp_send_.getPointer(), n_sp_send_, n_disp_sp_send_,
                            sp_recv_.getPointer(), n_sp_recv_, n_disp_sp_recv_);
#else
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
    #if 0
            const S32 n_size_crit = 512 / sizeof(Tspj);
            //const S32 n_size_crit = 10;
            const S32 tag = 0;
            S32 n_proc_sp_send_tmp = 0;
            S32 n_proc_sp_recv_tmp = 0;
            for(S32 j=0; j<n_proc_sp_send; j++){
                const S32 rank_send = rank_sp_isend_[j];
                const S32 n_sp_send = n_sp_send_[rank_send];
                if(n_sp_send < n_size_crit) continue;
                n_sp_send_tot_isendrecv_loc += n_sp_send;
                const S32 n_disp_sp_send = n_disp_sp_send_[rank_send];
                const S32 tag = 1;
                Tspj * sp_send = sp_send_.getPointer(n_disp_sp_send);
                req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend(sp_send, n_sp_send, GetDataType<Tspj>(), rank_send, tag);
            }
            for(S32 j=0; j<n_proc_sp_recv; j++){
                const S32 rank_recv = rank_sp_irecv_[j];
                const S32 n_sp_recv = n_sp_recv_[rank_recv];
                if(n_sp_recv < n_size_crit) continue;
                const S32 n_disp_sp_recv = n_disp_sp_recv_[rank_recv];
                const S32 tag = 1;
                Tspj * sp_recv = sp_recv_.getPointer(n_disp_sp_recv);
                req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv(sp_recv, n_sp_recv, GetDataType<Tspj>(), rank_recv, tag);
            }
            MPI::Request::Waitall(n_proc_sp_send_tmp, req_sp_send_.getPointer());
            MPI::Request::Waitall(n_proc_sp_recv_tmp, req_sp_recv_.getPointer());
            wtime_alltoallv_sp_isendrecv_1st += GetWtime() - wtime_offset_isendrecv;
            F64 wtime_offset_isendrecv_2nd = GetWtime();
            n_proc_sp_send_tmp = n_proc_sp_recv_tmp = 0;
            for(S32 j=0; j<n_proc_sp_send; j++){
                const S32 rank_send = rank_sp_isend_[j];
                const S32 n_sp_send = n_sp_send_[rank_send];
                if(n_sp_send >= n_size_crit) continue;
                n_sp_send_tot_isendrecv_loc += n_sp_send;
                const S32 n_disp_sp_send = n_disp_sp_send_[rank_send];
                const S32 tag = 1;
                Tspj * sp_send = sp_send_.getPointer(n_disp_sp_send);
                req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend(sp_send, n_sp_send, GetDataType<Tspj>(), rank_send, tag);
            }
            for(S32 j=0; j<n_proc_sp_recv; j++){
                const S32 rank_recv = rank_sp_irecv_[j];
                const S32 n_sp_recv = n_sp_recv_[rank_recv];
                if(n_sp_recv >= n_size_crit) continue;
                const S32 n_disp_sp_recv = n_disp_sp_recv_[rank_recv];
                const S32 tag = 1;
                Tspj * sp_recv = sp_recv_.getPointer(n_disp_sp_recv);
                req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv(sp_recv, n_sp_recv, GetDataType<Tspj>(), rank_recv, tag);
            }
            MPI::Request::Waitall(n_proc_sp_send_tmp, req_sp_send_.getPointer());
            MPI::Request::Waitall(n_proc_sp_recv_tmp, req_sp_recv_.getPointer());
            wtime_alltoallv_sp_isendrecv_2nd += GetWtime() - wtime_offset_isendrecv_2nd;
    #elif 0
            assert(dinfo != NULL);
            S32 law = 4;
            for(S32 i=0; i<law; i++){
                //Comm::barrier(); if(my_rank == 0) std::cerr<<"i(0)= "<<i<<std::endl;
                S32 n_proc_sp_send_tmp = 0;
                S32 n_proc_sp_recv_tmp = 0;
                S32 tag = 0;
                if(my_rank%law == 0){
                    //if(my_rank == 0) std::cerr<<"(0-i)%law= "<<(0-i)%law<<std::endl;
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        if(rank_send%law==(0+law+i)%law){
                            //if(Comm::getRank()==0) std::cerr<<"rank_send= "<<rank_send<<" n_sp_send_[rank_send]= "<<n_sp_send_[rank_send]<<std::endl;
                            n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        if(rank_recv%law==(0+law-i)%law){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);
                        }
                    }
                }
                if(my_rank%law == 1){
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        if(rank_send%law==(1+law+i)%law){
                            n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        if(rank_recv%law==(1+law-i)%law){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);
                        }
                    }
                }
                if(my_rank%law == 2){
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        if(rank_send%law==(2+law+i)%law){
                            n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        if(rank_recv%law==(2+law-i)%law){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);
                        }
                    }
                }
                if(my_rank%law == 3){
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        if(rank_send%law==(3+law+i)%law){
                            n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        if(rank_recv%law==(3+law-i)%law){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);
                        }
                    }                    
                }
                MPI::Request::Waitall(n_proc_sp_send_tmp, req_sp_send_.getPointer());
                MPI::Request::Waitall(n_proc_sp_recv_tmp, req_sp_recv_.getPointer());
                //MPI::COMM_WORLD.Barrier();
            }

            /*
            S32 law = 8;
            for(S32 i=0; i<law; i++){
                //Comm::barrier(); if(my_rank == 0) std::cerr<<"i(0)= "<<i<<std::endl;
                S32 n_proc_sp_send_tmp = 0;
                S32 n_proc_sp_recv_tmp = 0;
                S32 tag = 0;
                if(my_rank%law == 0){
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        if(rank_send%law==(0+law+i)%law){
                            n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        if(rank_recv%law==(0+law-i)%law){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);
                        }
                    }
                }
                if(my_rank%law == 1){
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                        if(rank_send%law==(1+law+i)%law){
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        if(rank_recv%law==(1+law-i)%law){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);
                        }
                    }
                }
                if(my_rank%law == 2){
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        if(rank_send%law==(2+law+i)%law){
                            n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        if(rank_recv%law==(2+law-i)%law){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);
                        }
                    }
                }
                if(my_rank%law == 3){
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        if(rank_send%law==(3+law+i)%law){
                            n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        if(rank_recv%law==(3+law-i)%law){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);
                        }
                    }                    
                }
                if(my_rank%law == 4){
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        if(rank_send%law==(4+law+i)%law){
                            n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        if(rank_recv%law==(4+law-i)%law){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);
                        }
                    }
                }
                if(my_rank%law == 5){
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        if(rank_send%law==(5+law+i)%law){
                            n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        if(rank_recv%law==(5+law-i)%law){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);
                        }
                    }
                }
                if(my_rank%law == 6){
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        if(rank_send%law==(6+law+i)%law){
                            n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        if(rank_recv%law==(6+law-i)%law){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);
                        }
                    }
                }
                if(my_rank%law == 7){
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        if(rank_send%law==(7+law+i)%law){
                            n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        if(rank_recv%law==(7+law-i)%law){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);
                        }
                    }                    
                }
                MPI::Request::Waitall(n_proc_sp_send_tmp, req_sp_send_.getPointer());
                MPI::Request::Waitall(n_proc_sp_recv_tmp, req_sp_recv_.getPointer());
            }            
            */
    #elif 0
            assert(dinfo != NULL);
            for(S32 i=0; i<4; i++){
                S32 n_proc_sp_send_tmp = 0;
                S32 n_proc_sp_recv_tmp = 0;
                const S32 my_rank_x = dinfo->getRank1d(0);
                const S32 my_rank_y = dinfo->getRank1d(1);
                const S32 n_proc_y = dinfo->getNDomain(1);;
                const S32 tag = i+1;
                if(i==0){
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        const S32 rank_send_x = rank_send / n_proc_y;
                        const S32 rank_send_y = rank_send % n_proc_y;
                        if(rank_send_x >= my_rank_x && rank_send_y > my_rank_y){
                            n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        const S32 rank_recv_x = rank_recv / n_proc_y;
                        const S32 rank_recv_y = rank_recv % n_proc_y;
                        if(rank_recv_x <= my_rank_x && rank_recv_y < my_rank_y){
                            //if(rank_recv_x >= my_rank_x && rank_recv_y > my_rank_y){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);                            
                        }
                    }
                }
                else if(i==1){
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        const S32 rank_send_x = rank_send / n_proc_y;
                        const S32 rank_send_y = rank_send % n_proc_y;
                        if(rank_send_x < my_rank_x && rank_send_y >= my_rank_y){
                            n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        const S32 rank_recv_x = rank_recv / n_proc_y;
                        const S32 rank_recv_y = rank_recv % n_proc_y;
                        if(rank_recv_x > my_rank_x && rank_recv_y <= my_rank_y){
                            //if(rank_recv_x < my_rank_x && rank_recv_y >= my_rank_y){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);                            
                        }
                    }
                }
                else if(i==2){
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        const S32 rank_send_x = rank_send / n_proc_y;
                        const S32 rank_send_y = rank_send % n_proc_y;
                        if(rank_send_x <= my_rank_x && rank_send_y < my_rank_y){
                            n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        const S32 rank_recv_x = rank_recv / n_proc_y;
                        const S32 rank_recv_y = rank_recv % n_proc_y;
                        if(rank_recv_x >= my_rank_x && rank_recv_y > my_rank_y){
                            //if(rank_recv_x <= my_rank_x && rank_recv_y < my_rank_y){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);                            
                        }
                    }
                }
                else if(i==3){
                    for(S32 j=0; j<n_proc_sp_send; j++){
                        const S32 rank_send   = rank_sp_isend_[j];
                        const S32 rank_send_x = rank_send / n_proc_y;
                        const S32 rank_send_y = rank_send % n_proc_y;
                        if(rank_send_x > my_rank_x && rank_send_y <= my_rank_y){
                            n_sp_send_tot_isendrecv_loc += n_sp_send_[rank_send];
                            req_sp_send_[n_proc_sp_send_tmp++] = MPI::COMM_WORLD.Isend
                                (sp_send_.getPointer(n_disp_sp_send_[rank_send]), n_sp_send_[rank_send], GetDataType<Tspj>(), rank_send, tag);
                        }
                    }
                    for(S32 j=0; j<n_proc_sp_recv; j++){
                        const S32 rank_recv   = rank_sp_irecv_[j];
                        const S32 rank_recv_x = rank_recv / n_proc_y;
                        const S32 rank_recv_y = rank_recv % n_proc_y;
                        if(rank_recv_x < my_rank_x && rank_recv_y >= my_rank_y){
                            //if(rank_recv_x > my_rank_x && rank_recv_y <= my_rank_y){
                            req_sp_recv_[n_proc_sp_recv_tmp++] = MPI::COMM_WORLD.Irecv
                                (sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]), n_sp_recv_[rank_recv], GetDataType<Tspj>(), rank_recv, tag);                            
                        }
                    }
                }
                MPI::Request::Waitall(n_proc_sp_send_tmp, req_sp_send_.getPointer());
                MPI::Request::Waitall(n_proc_sp_recv_tmp, req_sp_recv_.getPointer());
            }
    #elif 0
            for(int i=1; i<n_proc; i++){
                const S32 rank_send = (my_rank+n_proc+i) % n_proc;
                const S32 rank_recv = (my_rank+n_proc-i) % n_proc;
                const S32 n_sp_send = (n_sp_send_[rank_send] > 1) ? n_sp_send_[rank_send] : 0;
                const S32 n_sp_recv = (n_sp_recv_[rank_recv] > 1) ? n_sp_recv_[rank_recv] : 0;
                n_sp_send_tot_isendrecv_loc += n_sp_send;
                if(my_rank % (i*2) < i){
                    if(n_sp_send > 0){
                        Tspj * sp_send = sp_send_.getPointer(n_disp_sp_send_[rank_send]);
                        MPI::COMM_WORLD.Send(sp_send, n_sp_send, GetDataType<Tspj>(), rank_send, 0);
                    }
                    if(n_sp_recv > 0){
                        Tspj * sp_recv = sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]);
                        MPI::COMM_WORLD.Recv(sp_recv, n_sp_recv, GetDataType<Tspj>(), rank_recv, 0);
                    }
                }
                else{
                    if(n_sp_recv > 0){
                        Tspj * sp_recv = sp_recv_.getPointer(n_disp_sp_recv_[rank_recv]);
                        MPI::COMM_WORLD.Recv(sp_recv, n_sp_recv, GetDataType<Tspj>(), rank_recv, 0);
                    }
                    if(n_sp_send > 0){
                        Tspj * sp_send = sp_send_.getPointer(n_disp_sp_send_[rank_send]);
                        MPI::COMM_WORLD.Send(sp_send, n_sp_send, GetDataType<Tspj>(), rank_send, 0);
                    }
                }
            }
    #else
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
    #endif
            wtime_alltoallv_sp_isendrecv += GetWtime() - wtime_offset_isendrecv;
            //band_width_sp = n_sp_send_tot_loc*sizeof(Tepj) / wtime_alltoallv_ep;
            
            for(S32 i=0; i<n_proc; i++){
                if(n_sp_recv_[i] == 1){
                    sp_recv_[n_disp_sp_recv_[i]] = sp_top_recv_[i];
                }
            }
            /*
            // for debug
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
            */
#endif
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
        
    };
}


