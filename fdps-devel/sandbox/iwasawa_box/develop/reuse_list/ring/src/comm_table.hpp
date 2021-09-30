#pragma once

#include<particle_simulator.hpp>

namespace ParticleSimulator{
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
        ReallocatableArray<S32> id_ep_send_;
        ReallocatableArray<S32> id_sp_send_;
        ReallocatableArray<Tepj> ep_send_;
        ReallocatableArray<Tspj> sp_send_;
        ReallocatableArray<Tepj> ep_recv_;
        ReallocatableArray<Tspj> sp_recv_;
        S32 * n_ep_sp_send_;
        S32 * n_ep_sp_recv_;
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
        }

        void clearSize(){
            id_ep_send_.clearSize();
            id_sp_send_.clearSize();
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
            /*
            if(Comm::getRank() == 0){
                for(S32 i=0; i<10; i++){
                    std::cerr<<"ep_send_[i].pos= "<<ep_send_[i].pos<<std::endl;
                }
            }
            */
            Comm::allToAllV(ep_send_.getPointer(), n_ep_send_, n_disp_ep_send_,
                            ep_recv_.getPointer(), n_ep_recv_, n_disp_ep_recv_);
            Comm::allToAllV(sp_send_.getPointer(), n_sp_send_, n_disp_sp_send_,
                            sp_recv_.getPointer(), n_sp_recv_, n_disp_sp_recv_);
        }
    };
}
