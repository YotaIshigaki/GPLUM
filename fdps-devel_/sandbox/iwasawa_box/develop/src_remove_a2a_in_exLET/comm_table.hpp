#pragma once

#include<particle_simulator.hpp>

namespace ParticleSimulator{

    class CommTable{
    public:
        S32 n_ep_send_tot_, n_ep_recv_tot_, n_sp_send_tot_, n_sp_recv_tot_;
        ReallocatableArray<S32> n_ep_send_;
        ReallocatableArray<S32> n_ep_recv_;
        ReallocatableArray<S32> adr_ep_send_;
        
        ReallocatableArray<S32> n_sp_send_;
        ReallocatableArray<S32> n_sp_recv_;
        ReallocatableArray<S32> adr_sp_send_;

        ReallocatableArray<F64vec> shift_per_image_;
        ReallocatableArray<S32> n_image_per_proc_;

        ReallocatableArray<S32> n_ep_per_image_;
        ReallocatableArray<S32> n_sp_per_image_;

        ReallocatableArray<S32> rank_send_;
        ReallocatableArray<S32> rank_recv_;

        ReallocatableArray<S32> rank_recv_allgather_;
        ReallocatableArray<F64ort> pos_domain_allgather_;
        
    public:

        void setPosDomainAllgather(const DomainInfo & dinfo,
                                   const F64ort & pos_root_cell){
            pos_domain_allgather_.resizeNoInitialize(rank_recv_allgather_.size());
            for(S32 i=0; i<pos_domain_allgather_.size(); i++){
                pos_domain_allgather_[i] = dinfo.getPosDomain(pos_root_cell, rank_recv_allgather_[i]);
            }
        }
        
        void initialize(){
            const S32 n_proc = Comm::getNumberOfProc();
            n_ep_send_.resizeNoInitialize(n_proc);
            n_ep_recv_.resizeNoInitialize(n_proc);
            n_sp_send_.resizeNoInitialize(n_proc);
            n_sp_recv_.resizeNoInitialize(n_proc);
            n_ep_send_tot_ = n_ep_recv_tot_ = n_sp_send_tot_ = n_sp_recv_tot_ = 0;
        }

        void clearSize(){
            adr_ep_send_.clearSize();
            adr_sp_send_.clearSize();
            shift_per_image_.clearSize();
            n_ep_per_image_.clearSize();
            n_sp_per_image_.clearSize();
        }

        void clearCounter(){
            //n_ep_send_tot_ = n_ep_recv_tot_ = n_sp_send_tot_ = n_sp_recv_tot_ = 0;
        }
        
        void setCommRank(){
            const S32 n_proc = Comm::getNumberOfProc();
            for(S32 i=0; i<n_proc; i++){
                if(n_ep_send_[i] > 0 || n_sp_send_[i] > 0){
                    rank_send_.push_back(i);
                }
                if(n_ep_recv_[i] > 0 || n_sp_recv_[i] > 0){
                    rank_recv_.push_back(i);
                }
            }
        }

        void setCommRankShort(){
            const S32 n_proc = Comm::getNumberOfProc();
            for(S32 i=0; i<n_proc; i++){
                if(n_ep_send_[i] > 0){
                    rank_send_.push_back(i);
                }
                if(n_ep_recv_[i] > 0){
                    rank_recv_.push_back(i);
                }
            }
        }
    };
}
