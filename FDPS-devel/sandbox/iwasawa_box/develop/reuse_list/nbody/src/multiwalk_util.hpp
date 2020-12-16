#include<particle_simulator.hpp>

namespace ParticleSimulator{

    class MultiWalkInfo{
    private:
        S32 n_thread_;
        S32 n_walk_limit_;
    public:
        std::vector<S32> walk_grp_;
        std::vector<S32> walk_grp_disp_;
        S32  * iw2ith_;
        S32  * iw2cnt_;
        S32  * n_epi_array_;
        S32  * n_epi_array_prev_;
        S32  * n_epj_array_;
        S32  * n_spj_array_;
        S32 ** n_epj_disp_thread_; // [n_thread][n_walk]
        S32 ** n_spj_disp_thread_;// [n_thread][n_walk]
        S32  * cnt_thread_;
        S32  * n_ep_cum_thread_;
        S32  * n_sp_cum_thread_;
        S64  * n_interaction_ep_ep_array_;
        S64  * n_interaction_ep_sp_array_;

        void initialize(const S32 _n_walk_limit, const S32 _n_thread){
            n_walk_limit_ = _n_walk_limit;
            n_thread_ = _n_thread;
            iw2ith_ = new S32[n_walk_limit_];
            iw2cnt_ = new S32[n_walk_limit_];
            n_epi_array_ = new S32[n_walk_limit_];
            n_epi_array_prev_ = new S32[n_walk_limit_];
            n_epj_array_ = new S32[n_walk_limit_];
            n_spj_array_ = new S32[n_walk_limit_];
            n_epj_disp_thread_ = new S32*[n_thread_];
            n_spj_disp_thread_ = new S32*[n_thread_];
            for(int i=0; i<n_thread_; i++){
                n_epj_disp_thread_[i] = new S32[n_walk_limit_];
                n_spj_disp_thread_[i] = new S32[n_walk_limit_];
            }
            cnt_thread_ = new S32[n_thread_];
            n_ep_cum_thread_ = new S32[n_thread_];
            n_sp_cum_thread_ = new S32[n_thread_];
            n_interaction_ep_ep_array_ = new S64[n_thread_];
            n_interaction_ep_sp_array_ = new S64[n_thread_];
        }

        void setWalkGrp(const S32 n_ipg, const S32 n_loop_max){
            walk_grp_.resize(n_loop_max); // group of walk (walk_grp[i] < n_walk_limit)
            walk_grp_disp_.resize(n_loop_max+1);
            walk_grp_disp_[0] = 0;
            for(S32 wg=0; wg<n_ipg%n_loop_max; wg++){
                walk_grp_[wg] = n_ipg / n_loop_max + 1;
                walk_grp_disp[wg+1] = walk_grp_disp[wg] + walk_grp[wg];
            }
            for(S32 wg=n_ipg%n_loop_max; wg<n_loop_max; wg++){
                walk_grp[wg] = n_ipg / n_loop_max;
                walk_grp_disp[wg+1] = walk_grp_disp[wg] + walk_grp[wg];
            }
            for(S32 i=0; i<n_thread; i++){
                n_interaction_ep_ep_array[i] = n_interaction_ep_sp_array[i] = 0;
            }
        }
    };

}

