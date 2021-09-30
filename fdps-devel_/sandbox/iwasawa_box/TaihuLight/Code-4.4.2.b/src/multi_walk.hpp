#pragma once

namespace ParticleSimulator{
    template<class Tepi, class Tepj, class Tspj, class Tforce>
    class MultiWalkInfo{
    public:
        bool initialized_;
        S32 n_thread_;
        S32 n_walk_limit_;
        S32 * n_epi_;
        S32 * n_epi_prev_;
        S32 * n_epj_;
        S32 * n_spj_;
        Tepi ** epi_; // array of pointer
        Tepj ** epj_; // array of pointer
        Tspj ** spj_; // array of pointer
        Tforce ** force_;
        Tforce ** force_prev_;
        S32  ** adr_epj_; // array of pointer
        S32  ** adr_spj_; // array of pointer
        
        MultiWalkInfo(): initialized_(false){}
        void initialize(const S32 _n_walk_limit){
            if(initialized_) return;
            n_thread_ = Comm::getNumberOfThread();
            n_walk_limit_ = _n_walk_limit;
            n_epi_ = new S32[n_walk_limit_];
            n_epi_prev_ = new S32[n_walk_limit_];
            n_epj_ = new S32[n_walk_limit_];
            n_spj_ = new S32[n_walk_limit_];
            epi_        = new Tepi*[n_walk_limit_];
            epj_        = new Tepj*[n_walk_limit_];
            spj_        = new Tspj*[n_walk_limit_];
            force_      = new Tforce*[n_walk_limit_];
            force_prev_ = new Tforce*[n_walk_limit_];
            adr_epj_    = new S32*[n_walk_limit_];
            adr_spj_    = new S32*[n_walk_limit_];            
            initialized_ = true;
        }
        S32 getNWalkGrp(const S32 _n_ipg){
            return _n_ipg/n_walk_limit_ + (_n_ipg%n_walk_limit_==0 ? 0 : 1);
        }
    };
}
