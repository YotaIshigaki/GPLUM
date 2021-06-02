#pragma once

#include<sstream>

#ifdef SUNWAY
extern "C"{
    #include <athread.h>
    #include "cpe_func.h"
    extern void SLAVE_FUN(RotateAndShiftZ)(void*);
}
#endif

extern std::ofstream fout_calc_force_1st_min;
extern std::ofstream fout_calc_force_1st_max;
extern std::ofstream fout_ex_ptcl_min;
extern std::ofstream fout_ex_ptcl_max;

extern std::string s_ex_ptcl_max;
extern std::string s_ex_ptcl_min;
extern std::string s_calc_force_1st_max;
extern std::string s_calc_force_1st_min;

//double TIME_SYSTEM = 0.0; // GLOBAL
double time_sys = 0.0; // GLOBAL

double wtime_calc_force = 0.0;
double wtime_calc_force_d_d = 0.0;
double wtime_calc_force_d_s = 0.0;
double wtime_calc_force_p_ds = 0.0;

struct DataForOutput{
    PS::S32 rank;
    PS::F64 time_sys;
    PS::S32 n_loc;
    PS::S32 n_ipg;
    PS::CountT n_epi_ave_loc;
    PS::CountT n_epj_ave_loc;
    PS::CountT n_spj_ave_loc;
    void dump(std::ofstream & fout){
        fout<<"rank= "<<rank<<std::endl;
        fout<<"time_sys= "<<time_sys<<std::endl;
        fout<<"n_loc= "<<n_loc<<std::endl;
        fout<<"n_ipg= "<<n_ipg<<std::endl;
        fout<<"n_epi_ave_loc= "<<n_epi_ave_loc<<std::endl;
        fout<<"n_epj_ave_loc= "<<n_epj_ave_loc<<std::endl;
        fout<<"n_spj_ave_loc= "<<n_spj_ave_loc<<std::endl;
    }
};

void AllreduceMaxLoc(const long val_in,
                     const int rank_in,
                     long & val_out,
                     int & rank_out){
    struct{
        long val;
        int rank;
    } loc, glb;
    loc.val  = val_in;
    loc.rank = rank_in;
    MPI::COMM_WORLD.Allreduce(&loc, &glb, 1, MPI::LONG_INT, MPI::MAXLOC);
    val_out = glb.val;
    rank_out = glb.rank;
}

void AllreduceMinLoc(const long val_in,
                     const int rank_in,
                     long & val_out,
                     int & rank_out){
    struct{
        long val;
        int rank;
    } loc, glb;
    loc.val  = val_in;
    loc.rank = rank_in;
    MPI::COMM_WORLD.Allreduce(&loc, &glb, 1, MPI::LONG_INT, MPI::MINLOC);
    val_out = glb.val;
    rank_out = glb.rank;
}


void AllreduceMaxLoc(const double val_in, const int rank_in,
                     double & val_out, int & rank_out){
    struct{
        double val;
        int rank;
    } loc, glb;
    loc.val  = val_in;
    loc.rank = rank_in;
    MPI::COMM_WORLD.Allreduce(&loc, &glb, 1, MPI::DOUBLE_INT, MPI::MAXLOC);
    val_out = glb.val;
    rank_out = glb.rank;
}

void AllreduceMinLoc(const double val_in,
                     const int rank_in,
                     double & val_out,
                     int & rank_out){
    struct{
        double val;
        int rank;
    } loc, glb;
    loc.val  = val_in;
    loc.rank = rank_in;
    MPI::COMM_WORLD.Allreduce(&loc, &glb, 1, MPI::DOUBLE_INT, MPI::MINLOC);
    val_out = glb.val;
    rank_out = glb.rank;
}

void DumpProfileParticleSystem(const PS::TimeProfile & prof,
                               std::ofstream & fout){
    fout<<"exchange_particle= "<<prof.exchange_particle<<std::endl;
    fout<<"    find_particle= "<<prof.exchange_particle__find_particle<<std::endl;
    fout<<"    communication= "<<prof.exchange_particle__exchange_particle<<std::endl;
}

void DumpProfileDomain(const PS::TimeProfile & prof,
                       std::ofstream & fout){
    fout<<"decompose_domain= "<<prof.decompose_domain<<std::endl;
}

void DumpProfileTree(const PS::TimeProfile & prof,
                     std::ofstream & fout){
    fout<<"calc_force_all= "<<prof.calc_force_all<<std::endl;
    fout<<"    set_particle_local_tree= "<<prof.set_particle_local_tree<<std::endl;
    fout<<"    morton_sort_local_tree= "<<prof.morton_sort_local_tree<<std::endl;
    //fout<<"        wtime_make_key_local_tree= "<<PS::wtime_make_key_local_tree<<std::endl;
    //fout<<"        wtime_sort_local_tree= "<<PS::wtime_sort_local_tree<<std::endl;
    fout<<"    morton_sort_FP= "<<prof.morton_sort_FP<<std::endl;
    fout<<"    link_cell_local_tree= "<<prof.link_cell_local_tree<<std::endl;
    fout<<"    calc_moment_local_tree= "<<prof.calc_moment_local_tree<<std::endl;
    fout<<"    add_moment_as_sp_local_tree= "<<prof.add_moment_as_sp_local_tree<<std::endl;
    fout<<"    exchange_LET_tot= "<<prof.exchange_LET_tot<<std::endl;
    //fout<<"        wtime_ex_let_allgather_root_ex= "<<PS::wtime_ex_let_allgather_root_ex<<std::endl;
    fout<<"        wtime_ex_let_sd_2_allgather_bd= "<<PS::wtime_ex_let_sd_2_allgather_bd<<std::endl;
    fout<<"        wtime_ex_let_sd_2_ex_n_d_d= "<<PS::wtime_ex_let_sd_2_ex_n_d_d<<std::endl;
    fout<<"        wtime_ex_let_sd_2_ex_n_d_sd_1d_ring= "<<PS::wtime_ex_let_sd_2_ex_n_d_sd_1d_ring<<std::endl;
    fout<<"        wtime_ex_let_sd_2_ex_n_d_sd_sub= "<<PS::wtime_ex_let_sd_2_ex_n_d_sd_sub<<std::endl;
    fout<<"        wtime_ex_let_sd_2_allgather_sd= "<<PS::wtime_ex_let_sd_2_allgather_sd<<std::endl;
    fout<<"        wtime_ex_let_sd_2_ep_sp_among_sd= "<<PS::wtime_ex_let_sd_2_ep_sp_among_sd<<std::endl;
    fout<<"        wtime_ex_let_sd_2_ep_sp_in_sd= "<<PS::wtime_ex_let_sd_2_ep_sp_in_sd<<std::endl;
    fout<<"        wtime_ex_let_sd_2_make_list= "<<PS::wtime_ex_let_sd_2_make_list<<std::endl;
    /*
    fout<<"        ex_let_allgather_root_ex= "<<PS::wtime_ex_let_allgather_root_ex<<std::endl;
    fout<<"        ex_#_of_ptcl= "<<PS::wtime_ex_let_exchange_number<<std::endl;
    fout<<"        make_LET= "<<prof.make_LET_1st<<std::endl;
    fout<<"        set_ptcl_send_buf= "<<prof.make_LET_2nd<<std::endl;
    fout<<"        wtime_alltoallv_ep= "<<PS::wtime_alltoallv_ep<<std::endl;
    fout<<"        wtime_alltoallv_sp= "<<PS::wtime_alltoallv_sp<<std::endl;
    fout<<"            wtime_alltoallv_sp_allgather= "<<PS::wtime_alltoallv_sp_allgather<<std::endl;
    fout<<"            wtime_alltoallv_sp_isendrecv= "<<PS::wtime_alltoallv_sp_isendrecv<<std::endl;
    //fout<<"                wtime_alltoallv_sp_isendrecv_1st= "<<PS::wtime_alltoallv_sp_isendrecv_1st<<std::endl;
    //fout<<"                wtime_alltoallv_sp_isendrecv_2nd= "<<PS::wtime_alltoallv_sp_isendrecv_2nd<<std::endl;
    */
    fout<<"    copy_tc_to_etc_loc= "<<prof.copy_tc_to_etc_loc<<std::endl;
    fout<<"    set_particle_global_tree= "<<prof.set_particle_global_tree<<std::endl;
    fout<<"    morton_sort_global_tree= "<<prof.morton_sort_global_tree<<std::endl;
    fout<<"    link_cell_global_tree= "<<prof.link_cell_global_tree<<std::endl;
    fout<<"    calc_moment_global_tree= "<<prof.calc_moment_global_tree<<std::endl;
    fout<<"    add_moment_as_sp_global_tree= "<<prof.add_moment_as_sp_global_tree<<std::endl;
    fout<<"    make_ipgroup= "<<prof.calc_force__make_ipgroup<<std::endl;
    fout<<"    make_all_interaction_list_id= "<<prof.make_all_interaction_list_id<<std::endl;
    fout<<"    copy_tc_to_etc_glb= "<<prof.copy_tc_to_etc_glb<<std::endl;
    fout<<"    balance_nwalk= "<<prof.balance_nwalk<<std::endl;
        
    fout<<"    calc_force= "<<prof.calc_force<<std::endl;
    fout<<"        wtime_dispatch= "<<PS::wtime_dispatch<<std::endl;
    //fout<<"            wtime_kick= "<<wtime_kick<<std::endl;
    //fout<<"            wtime_drift= "<<wtime_drift<<std::endl;
    //fout<<"            wtime_calc_energy= "<<wtime_calc_energy<<std::endl;
    fout<<"        wtime_prefixsum_recorder= "<<PS::wtime_prefixsum_recorder<<std::endl;
    //fout<<"        wtime_retrieve= "<<PS::wtime_retrieve<<std::endl;
    fout<<"  calc_force_all-calc_force= "<<prof.calc_force_all-prof.calc_force<<std::endl;
    //fout<<"    copy_back= "<<prof.calc_force__copy_original_order<<std::endl;
    fout<<"  sum= "
        <<prof.set_particle_local_tree
        + prof.morton_sort_local_tree
        + prof.morton_sort_FP
        + prof.link_cell_local_tree
        + prof.calc_moment_local_tree
        + prof.add_moment_as_sp_local_tree
        + prof.exchange_LET_tot
        + prof.copy_tc_to_etc_loc
        + prof.set_particle_global_tree
        + prof.morton_sort_global_tree
        + prof.link_cell_global_tree
        + prof.calc_moment_global_tree
        + prof.add_moment_as_sp_global_tree
        + prof.calc_force__make_ipgroup
        + prof.make_all_interaction_list_id
        + prof.copy_tc_to_etc_glb
        + prof.balance_nwalk
        + prof.calc_force
        + prof.calc_force__copy_original_order
        <<std::endl;
    /*
    fout<<"calc_on_MPE= "
        <<prof.set_particle_local_tree
        + prof.morton_sort_local_tree
        + prof.morton_sort_FP
        + prof.link_cell_local_tree
        + prof.add_moment_as_sp_local_tree
        + prof.exchange_LET_tot
        + prof.set_particle_global_tree
        + prof.morton_sort_global_tree
        + prof.link_cell_global_tree
        + prof.add_moment_as_sp_global_tree
        + prof.calc_force__make_ipgroup
        + prof.make_all_interaction_list_id
        + prof.calc_force__copy_original_order
        <<std::endl;
    */
    fout<<std::endl;

}


template<class Tsys>
void ClearCounterAll(Tsys & system,
                     PS::DomainInfo & dinfo){
    TREE_POINTER->clearTimeProfile();
    TREE_POINTER->clearCounterAll();
    system.clearTimeProfile();
    system.clearCounterAll();
    dinfo.clearTimeProfile();
    wtime_calc_force = wtime_calc_force_d_d
        = wtime_calc_force_d_s = wtime_calc_force_p_ds = 0.0;
    wtime_copy_all_epj = wtime_copy_all_spj = wtime_calc_epj
        = wtime_calc_spj = wtime_calc_pj = 0.0;
    wtime_kick = wtime_drift = wtime_calc_energy = 0.0;
    PS::wtime_prefixsum_recorder = PS::wtime_dispatch = PS::wtime_retrieve = 0.0;
    PS::wtime_alltoallv_sp = PS::wtime_alltoallv_ep = 0.0;
    PS::wtime_alltoallv_sp_allgather = PS::wtime_alltoallv_sp_isendrecv = 0.0;
    PS::n_ep_send_tot_loc = PS::n_sp_send_tot_loc = PS::n_sp_send_tot_isendrecv_loc = 0;
    PS::wtime_alltoallv_sp_isendrecv_1st = PS::wtime_alltoallv_sp_isendrecv_2nd = 0.0;
    PS::wtime_ex_let_exchange_number = PS::wtime_ex_let_allgather_root_ex = 0.0;
    PS::wtime_ex_let_allgather_root_ex         = PS::wtime_ex_let_sd_2_allgather_bd
        = PS::wtime_ex_let_sd_2_ex_n_d_d       = PS::wtime_ex_let_sd_2_ex_n_d_sd_1d_ring
        = PS::wtime_ex_let_sd_2_ex_n_d_sd_sub  = PS::wtime_ex_let_sd_2_allgather_sd
        = PS::wtime_ex_let_sd_2_ep_sp_among_sd = PS::wtime_ex_let_sd_2_ep_sp_in_sd
        = PS::wtime_ex_let_sd_2_make_list = 0.0;
}

template<class Tsys>
void DumpTimingProfile(Tsys & system,
                       PS::DomainInfo & dinfo,
                       std::ofstream & fout_wtime){
    PS::TimeProfile prof_tree    = TREE_POINTER->getTimeProfile();
    PS::TimeProfile prof_system  = system.getTimeProfile();
    PS::TimeProfile prof_domain  = dinfo.getTimeProfile();

    PS::CountT n_fp_send_tot_glb = system.getNumberOfParticleSendGlobal();
        
    PS::CountT n_ep_send_tot_glb = PS::Comm::getSum(PS::n_ep_send_tot_loc);
    PS::CountT n_sp_send_tot_glb = PS::Comm::getSum(PS::n_sp_send_tot_loc);
    PS::CountT n_sp_send_tot_isendrecv_glb = PS::Comm::getSum(PS::n_sp_send_tot_isendrecv_loc);
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();

    PS::CountT n_epi_ave_loc = TREE_POINTER->getNumberOfEPIAverageLocal();
    PS::CountT n_epj_ave_loc = TREE_POINTER->getNumberOfEPJAverageLocal();
    PS::CountT n_spj_ave_loc = TREE_POINTER->getNumberOfSPJAverageLocal();
        
    PS::CountT n_epi_ave_glb = TREE_POINTER->getNumberOfEPIAverageGlobal();
    PS::CountT n_epj_ave_glb = TREE_POINTER->getNumberOfEPJAverageGlobal();
    PS::CountT n_spj_ave_glb = TREE_POINTER->getNumberOfSPJAverageGlobal();

    PS::S32 my_rank = PS::Comm::getRank();
    PS::CountT n_epi_max, n_epj_max, n_spj_max;
    PS::S32 rank_epi_max, rank_epj_max, rank_spj_max;
    AllreduceMaxLoc(n_epi_ave_loc, my_rank, n_epi_max, rank_epi_max);
    AllreduceMaxLoc(n_epj_ave_loc, my_rank, n_epj_max, rank_epj_max);
    AllreduceMaxLoc(n_spj_ave_loc, my_rank, n_spj_max, rank_spj_max);

    PS::CountT n_epi_min, n_epj_min, n_spj_min;
    PS::S32 rank_epi_min, rank_epj_min, rank_spj_min;
    AllreduceMinLoc(n_epi_ave_loc, my_rank, n_epi_min, rank_epi_min);
    AllreduceMinLoc(n_epj_ave_loc, my_rank, n_epj_min, rank_epj_min);
    AllreduceMinLoc(n_spj_ave_loc, my_rank, n_spj_min, rank_spj_min);

    PS::S64 n_glb = system.getNumberOfParticleGlobal();
    PS::CountT n_loc = system.getNumberOfParticleLocal();
    PS::CountT n_loc_min, n_loc_max;
    PS::S32 rank_n_loc_min, rank_n_loc_max;
    AllreduceMinLoc(n_loc, my_rank, n_loc_min, rank_n_loc_min);
    AllreduceMaxLoc(n_loc, my_rank, n_loc_max, rank_n_loc_max);
        
    PS::F64 wtime_ex_ptcl_max, wtime_ex_ptcl_min;
    PS::S32 rank_ex_ptcl_max, rank_ex_ptcl_min;
    AllreduceMinLoc(prof_system.exchange_particle, my_rank, wtime_ex_ptcl_min, rank_ex_ptcl_min);
    AllreduceMaxLoc(prof_system.exchange_particle, my_rank, wtime_ex_ptcl_max, rank_ex_ptcl_max);

    PS::F64 wtime_calc_force_1st_max, wtime_calc_force_1st_min;
    PS::S32 rank_calc_force_1st_max, rank_calc_force_1st_min;
    AllreduceMinLoc(PS::wtime_calc_force, my_rank, wtime_calc_force_1st_min, rank_calc_force_1st_min);
    AllreduceMaxLoc(PS::wtime_calc_force, my_rank, wtime_calc_force_1st_max, rank_calc_force_1st_max);

    long int disp_walk_ex_ptcl_max, disp_walk_ex_ptcl_min;
    PS::S32 rank_disp_walk_ex_ptcl_max, rank_disp_walk_ex_ptcl_min;
    AllreduceMinLoc(PS::DISP_WALK_EX_PTCL, my_rank, disp_walk_ex_ptcl_min, rank_disp_walk_ex_ptcl_min);
    AllreduceMaxLoc(PS::DISP_WALK_EX_PTCL, my_rank, disp_walk_ex_ptcl_max, rank_disp_walk_ex_ptcl_max);
        
    PS::S32 n_disp_walk_ex_ptcl_max = PS::N_LOC_ORG;// number of particle
    PS::S32 n_disp_walk_ex_ptcl_min = PS::N_LOC_ORG;
    MPI::COMM_WORLD.Bcast(&n_disp_walk_ex_ptcl_max, 1, MPI_INT, rank_disp_walk_ex_ptcl_max);
    MPI::COMM_WORLD.Bcast(&n_disp_walk_ex_ptcl_min, 1, MPI_INT, rank_disp_walk_ex_ptcl_min);

    PS::S32 n_not_move_max = PS::N_NOT_MOVE;
    PS::S32 n_not_move_min = PS::N_NOT_MOVE;
    MPI::COMM_WORLD.Bcast(&n_not_move_max, 1, MPI_INT, rank_disp_walk_ex_ptcl_max);
    MPI::COMM_WORLD.Bcast(&n_not_move_min, 1, MPI_INT, rank_disp_walk_ex_ptcl_min);

    PS::F64 hit_ratio_max = PS::HIT_RATIO;
    PS::F64 hit_ratio_min = PS::HIT_RATIO;
    MPI::COMM_WORLD.Bcast(&PS::HIT_RATIO, 1, MPI_DOUBLE, rank_disp_walk_ex_ptcl_max);
    MPI::COMM_WORLD.Bcast(&PS::HIT_RATIO, 1, MPI_DOUBLE, rank_disp_walk_ex_ptcl_min);
    
    DataForOutput data_for_output;
    data_for_output.rank = my_rank;
    data_for_output.time_sys = time_sys;
    data_for_output.n_loc = n_loc;
    data_for_output.n_ipg = TREE_POINTER->getNumberOfIPG();
    data_for_output.n_epi_ave_loc = n_epi_ave_loc;
    data_for_output.n_epj_ave_loc = n_epj_ave_loc;
    data_for_output.n_spj_ave_loc = n_spj_ave_loc;
    
    if(PS::Comm::getRank() == rank_calc_force_1st_min){
        fout_calc_force_1st_min.open(s_calc_force_1st_min.c_str(), std::fstream::app);
        data_for_output.dump(fout_calc_force_1st_min);
        DumpProfileDomain(prof_domain, fout_calc_force_1st_min);
        DumpProfileParticleSystem(prof_system, fout_calc_force_1st_min);
        DumpProfileTree(prof_tree, fout_calc_force_1st_min);
        fout_calc_force_1st_min<<std::endl;
        fout_calc_force_1st_min.close();
    }
    if(PS::Comm::getRank() == rank_calc_force_1st_max){
        fout_calc_force_1st_max.open(s_calc_force_1st_max.c_str(), std::fstream::app);
        data_for_output.dump(fout_calc_force_1st_max);
        DumpProfileDomain(prof_domain, fout_calc_force_1st_max);
        DumpProfileParticleSystem(prof_system, fout_calc_force_1st_max);
        DumpProfileTree(prof_tree, fout_calc_force_1st_max);
        fout_calc_force_1st_max<<std::endl;
        fout_calc_force_1st_max.close();
    }
    
    if(PS::Comm::getRank() == rank_ex_ptcl_min){
        fout_ex_ptcl_min.open(s_ex_ptcl_min.c_str(), std::fstream::app);
        data_for_output.dump(fout_ex_ptcl_min);
        DumpProfileDomain(prof_domain, fout_ex_ptcl_min);
        DumpProfileParticleSystem(prof_system, fout_ex_ptcl_min);
        DumpProfileTree(prof_tree, fout_ex_ptcl_min);
        fout_ex_ptcl_min<<std::endl;
        fout_ex_ptcl_min.close();
    }
    if(PS::Comm::getRank() == rank_ex_ptcl_max){
        fout_ex_ptcl_max.open(s_ex_ptcl_max.c_str(), std::fstream::app);
        data_for_output.dump(fout_ex_ptcl_max);
        DumpProfileDomain(prof_domain, fout_ex_ptcl_max);
        DumpProfileParticleSystem(prof_system, fout_ex_ptcl_max);
        DumpProfileTree(prof_tree, fout_ex_ptcl_max);
        fout_ex_ptcl_max<<std::endl;
        fout_ex_ptcl_max.close();
    }

    PS::S64 n_op_ep_ep = 47;
    PS::S64 n_op_ep_sp = 31;
    //PS::S64 n_op_st_ptcl = 62;
    PS::S64 n_op_st_ptcl = 54;
    PS::S64 n_int_ep_ep_glb = TREE_POINTER->getNumberOfInteractionEPEPGlobal();
    PS::S64 n_int_ep_sp_glb = TREE_POINTER->getNumberOfInteractionEPSPGlobal();
    PS::S64 n_int_st_ptcl_glb = N_SATELLITE*system.getNumberOfParticleGlobal();
    PS::F64 n_op_tot = (PS::F64)n_int_ep_ep_glb*(PS::F64)n_op_ep_ep
        + (PS::F64)n_int_ep_sp_glb*(PS::F64)n_op_ep_sp
        + (PS::F64)n_int_st_ptcl_glb*(PS::F64)n_op_st_ptcl;

    PS::S64 n_op_tot_s64 = (PS::S64)n_int_ep_ep_glb*(PS::S64)n_op_ep_ep
        + (PS::S64)n_int_ep_sp_glb*(PS::S64)n_op_ep_sp
        + (PS::S64)n_int_st_ptcl_glb*(PS::S64)n_op_st_ptcl;
        
    PS::S64 n_tree_cell_loc_loc = TREE_POINTER->adr_tc_level_partition_loc_[TREE_POINTER->lev_max_loc_+1];
    PS::S64 n_tree_cell_glb_loc = TREE_POINTER->adr_tc_level_partition_glb_[TREE_POINTER->lev_max_glb_+1];
    PS::S64 n_tree_cell_loc_glb = PS::Comm::getSum(n_tree_cell_loc_loc);
    PS::S64 n_tree_cell_glb_glb = PS::Comm::getSum(n_tree_cell_glb_loc);
        
    if(PS::Comm::getRank() == 0){
        fout_wtime<<std::setprecision(15);
        fout_wtime<<"rank= "<<PS::Comm::getRank()<<std::endl;
        fout_wtime<<"time_sys= "<<time_sys<<std::endl;
        fout_wtime<<"n_glb= "<<n_glb<<std::endl;
        fout_wtime<<"n_loc= "<<n_loc<<std::endl;
        fout_wtime<<"n_proc= "<<n_proc<<std::endl;
        fout_wtime<<"<n_epi>(local)= "<<n_epi_ave_loc<<std::endl;
        fout_wtime<<"<n_epj>(local)= "<<n_epj_ave_loc<<std::endl;
        fout_wtime<<"<n_spj>(local)= "<<n_spj_ave_loc<<std::endl;
        fout_wtime<<"<n_epi>(global)= "<<n_epi_ave_glb<<std::endl;
        fout_wtime<<"<n_epj>(global)= "<<n_epj_ave_glb<<std::endl;
        fout_wtime<<"<n_spj>(global)= "<<n_spj_ave_glb<<std::endl;
        fout_wtime<<"n_tree_cell_loc= "<<n_tree_cell_loc_loc
                  <<" n_tree_cell_glb= "<<n_tree_cell_glb_loc
                  <<std::endl;
        fout_wtime<<"n_tree_cell_loc(ave)= "<<(PS::F64)n_tree_cell_loc_glb/n_proc
                  <<" n_tree_cell_glb(ave)= "<<(PS::F64)n_tree_cell_glb_glb/n_proc
                  <<std::endl;
        fout_wtime<<"n_op_tot= "<<n_op_tot<<std::endl;
        fout_wtime<<"n_op_tot_s64= "<<n_op_tot_s64<<std::endl;
            
        fout_wtime<<"n_loc_max= "<<n_loc_max<<" rank_n_loc_max= "<<rank_n_loc_max<<std::endl;
        fout_wtime<<"n_loc_min= "<<n_loc_min<<" rank_n_loc_min= "<<rank_n_loc_min<<std::endl;
            
        fout_wtime<<"n_epi_max= "<<n_epi_max<<" rank_epi_max= "<<rank_epi_max<<std::endl;
        fout_wtime<<"n_epi_min= "<<n_epi_min<<" rank_epi_min= "<<rank_epi_min<<std::endl;
        fout_wtime<<"n_epj_max= "<<n_epj_max<<" rank_epj_max= "<<rank_epj_max<<std::endl;
        fout_wtime<<"n_epj_min= "<<n_epj_min<<" rank_epj_min= "<<rank_epj_min<<std::endl;
        fout_wtime<<"n_spj_max= "<<n_spj_max<<" rank_spj_max= "<<rank_spj_max<<std::endl;
        fout_wtime<<"n_spj_min= "<<n_spj_min<<" rank_spj_min= "<<rank_spj_min<<std::endl;
            
        fout_wtime<<"n_ipg= "<<TREE_POINTER->getNumberOfIPG()<<std::endl;
        fout_wtime<<"n_fp_send_tot_glb= "<<n_fp_send_tot_glb
                  <<" n_fp_send_tot_glb/n_proc= "<<n_fp_send_tot_glb/n_proc
                  <<std::endl;            
        fout_wtime<<"n_ep_send_tot_glb= "<<n_ep_send_tot_glb
                  <<" n_ep_send_tot_glb/n_proc= "<<n_ep_send_tot_glb/n_proc
                  <<std::endl;
        fout_wtime<<"n_sp_send_tot_glb= "<<n_sp_send_tot_glb
                  <<" n_sp_send_tot_glb/n_proc= "<<n_sp_send_tot_glb/n_proc
                  <<std::endl;
        fout_wtime<<"n_sp_send_tot_isendrecv_glb= "<<n_sp_send_tot_isendrecv_glb
                  <<" n_sp_send_tot_isendrecv_glb/n_proc= "<<n_sp_send_tot_isendrecv_glb/n_proc
                  <<std::endl;
        /*
          fout_wtime<<"band_width(ep, per proc)= "
          <<(double)(n_ep_send_tot_glb*sizeof(Epj)) / PS::wtime_alltoallv_ep / 1e9 / n_proc
          <<" [GB/s]"
          <<std::endl;
          fout_wtime<<"band_width(sp, per proc)= "
          <<(double)(n_sp_send_tot_glb*sizeof(SPJMonopoleScatterSW)) / PS::wtime_alltoallv_sp / 1e9 / n_proc
          <<" [GB/s]"
          <<std::endl;
          fout_wtime<<"band_width(sp(isend+irecv), per proc)= "
          <<(double)(n_sp_send_tot_isendrecv_glb*sizeof(SPJMonopoleScatterSW)) / PS::wtime_alltoallv_sp_isendrecv / 1e9 / n_proc
          <<" [GB/s]"
          <<std::endl;            
        */
        fout_wtime<<"speed= "<<n_op_tot / (prof_tree.calc_force_all+prof_system.exchange_particle) * 1e-12 <<" [Tflops]"<<std::endl;
        std::cerr<<"speed= "<<n_op_tot / (prof_tree.calc_force_all+prof_system.exchange_particle) * 1e-12 <<" [Tflops]"<<std::endl;
        DumpProfileDomain(prof_domain, fout_wtime);
        DumpProfileParticleSystem(prof_system, fout_wtime);
        DumpProfileTree(prof_tree, fout_wtime);
    }
    ClearCounterAll(system, dinfo);
}



template<class Tpsys, class Tsat, class Tpla>
void CalcForceFromPlanet(Tpsys & system,
                         Tsat & satellite,
                         const Tpla & planet){
    const PS::F64 eps2 = Epi::eps * Epi::eps;
    const PS::S32 n = system.getNumberOfParticleLocal();
    const PS::S32 n_sat = satellite.size();
    const PS::F64 mj = planet.mass;
    const PS::F64vec rj = planet.pos;
#pragma omp parallel
    for(PS::S32 i=0; i<n; i++) {
        const PS::F64vec rij = system[i].pos - rj;
        //const PS::F64 r_sq = rij * rij;
        const PS::F64 r_sq = rij * rij + eps2;
        const PS::F64 r_inv = 1.0 / sqrt(r_sq);
        const PS::F64 pij = mj * r_inv;
        const PS::F64 mri3 = pij * r_inv * r_inv;
        system[i].acc -= mri3 * rij;
        system[i].pot -= pij * 2.0; // factor two is just a trick to calculate energy
    }
#pragma omp parallel
    for(PS::S32 i=0; i<n_sat; i++) {
        const PS::F64vec rij = satellite[i].pos - rj;
        //const PS::F64 r_sq = rij * rij;
        const PS::F64 r_sq = rij * rij + eps2;        
        const PS::F64 r_inv = 1.0 / sqrt(r_sq);
        const PS::F64 pij = mj * r_inv;
        const PS::F64 mri3 = pij * r_inv * r_inv;
        satellite[i].acc -= mri3 * rij;
        satellite[i].pot -= pij * 2.0; // factor two is just a trick to calculate energy
    }
}

template<class Tpsys, class Tptcl, class Tforce>
void CalcForceBetweenDustAndSatellite(const Tpsys & system,
                                      Tptcl & satellite,
                                      Tforce force_dust[]){
    const PS::S32 n = system.getNumberOfParticleLocal();
    const PS::S32 n_sat = satellite.size();
    static ForceSat * force_sat_loc;
    static ForceSat * force_sat_glb;
    static NgbSat   * ngb_sat_loc;
    static NgbSat   * ngb_sat_glb;
    static PS::F64  * r_ngb_sq_sat;
    static RNgbSqRank * r_ngb_sq_rank_in;
    static RNgbSqRank * r_ngb_sq_rank_out;
    static bool first = true;
    if(first){
        force_sat_loc = new ForceSat[n_sat];
        force_sat_glb = new ForceSat[n_sat];
        ngb_sat_loc   = new NgbSat[n_sat];
        ngb_sat_glb   = new NgbSat[n_sat];
        r_ngb_sq_sat  = new PS::F64[n_sat];
        r_ngb_sq_rank_in  = new RNgbSqRank[n_sat];
        r_ngb_sq_rank_out = new RNgbSqRank[n_sat];
    }
    for(int i=0; i<n_sat; i++){
        r_ngb_sq_sat[i] = PS::LARGE_FLOAT;
        force_sat_loc[i].clear();
        force_sat_glb[i].clear();
        ngb_sat_loc[i].clear();
        ngb_sat_glb[i].clear();
    }
#pragma omp parallel
    for(PS::S32 i=0; i<n; i++) {
        const PS::F64vec ri = system[i].pos;
        //const PS::F64 r_coll_i = system[i].r_coll;
        const PS::F64 r_coll_i = Epj::r_coll;        
        force_dust[i].clear();
        force_dust[i].r_ngb_sq = PS::LARGE_FLOAT;
        for(PS::S32 j=0; j<n_sat; j++){
            const PS::F64 r_coll_sq = (satellite[j].r_coll+r_coll_i)*(satellite[j].r_coll+r_coll_i);
            const PS::F64 mj = satellite[j].mass;
            const PS::F64vec rij = ri - satellite[j].pos;
            const PS::F64 r_sq = rij * rij;
            const PS::F64 r_inv = 1.0 / sqrt(r_sq);
            const PS::F64 pij = mj * r_inv;
            const PS::F64 mri3 = pij * r_inv * r_inv;
            const PS::F64vec aij = mri3 * rij;
            force_dust[i].acc -= aij;
            force_dust[i].pot -= pij;
            force_sat_loc[j].acc += aij;
            force_sat_loc[j].pot -= pij;
            if(r_coll_sq < rij*rij){
                force_dust[i].n_coll++;
                force_sat_loc[j].n_coll += 1.0+1e-10;
                if(rij*rij < force_dust[i].r_ngb_sq){
                    // for dust
                    force_dust[i].r_ngb_sq = rij*rij;
                    //force_dust[i].ngb.mass = satellite[j].mass;
                    //force_dust[i].ngb.pos = satellite[j].pos;
                    //force_dust[i].ngb.vel = satellite[j].vel;
                }
                if(rij*rij < r_ngb_sq_sat[j]){
                    r_ngb_sq_sat[j] = rij*rij;
                    ngb_sat_loc[j].mass = system[i].mass;
                    ngb_sat_loc[j].pos = system[i].pos;
                    ngb_sat_loc[j].vel = system[i].vel;
                }
            }
        }
    }
    MPI_Allreduce( (double*)force_sat_loc,  (double*)force_sat_glb, 5*n_sat, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(PS::S32 i=0; i<n_sat; i++){
        r_ngb_sq_rank_in[i].r_ngb_sq = r_ngb_sq_sat[i];
        r_ngb_sq_rank_in[i].rank = PS::Comm::getRank();
    }
    for(PS::S32 i=0; i<n_sat; i++){
        MPI_Allreduce(r_ngb_sq_rank_in+i, r_ngb_sq_rank_out+i, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
        MPI_Bcast(ngb_sat_loc+i, 1, PS::GetDataType<NgbSat>(), r_ngb_sq_rank_out[i].rank, MPI_COMM_WORLD);
    }
    for(PS::S32 i=0; i<n_sat; i++){
        satellite[i].acc = force_sat_glb[i].acc;
        satellite[i].pot = force_sat_glb[i].pot;
    }
}

template<class Tptcl>
void kick(Tptcl ptcl[],
          const PS::S32 n,
          const PS::F64 dt) {
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++) {
        ptcl[i].vel  += ptcl[i].acc * dt;
    }
}

template<class Tptcl, class Tforce>
void kick(Tptcl ptcl[],
          const Tforce force[],
          const PS::S32 n,
          const PS::F64 dt) {
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++) {
        ptcl[i].vel  += force[i].acc * dt;
    }
}

template<class Tpsys>
void kick(Tpsys & system,
          const PS::F64 dt) {
    const PS::S32 n = system.getNumberOfParticleLocal();
    kick(&system[0], n, dt);
}


template<class Tptcl>
void drift(Tptcl ptcl[],
           const PS::S32 n,
           const PS::F64 dt) {
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++) {
        ptcl[i].pos  += ptcl[i].vel * dt;
    }
}

template<class Tpsys>
void drift(Tpsys & system,
           const PS::F64 dt) {
    const PS::S32 n = system.getNumberOfParticleLocal();
    drift(&system[0], n, dt);
}


void printHelp() {
    std::cerr<<"o: dir name of output (default: ./result)"<<std::endl;
    std::cerr<<"t: theta (default: 0.5)"<<std::endl;
    std::cerr<<"T: time_end (default: 10.0)"<<std::endl;
    std::cerr<<"s: time_step (default: 1.0 / 128.0)"<<std::endl;
    std::cerr<<"d: dt_diag (default: 1.0 / 8.0)"<<std::endl;
    std::cerr<<"D: dt_snap (default: 1.0)"<<std::endl;
    std::cerr<<"l: n_leaf_limit (default: 8)"<<std::endl;
    std::cerr<<"n: n_group_limit (default: 64)"<<std::endl;
    std::cerr<<"N: n_tot (default: 1024)"<<std::endl;
    std::cerr<<"h: help"<<std::endl;
}

void makeOutputDirectory(char * dir_name) {
    struct stat st;
    if(stat(dir_name, &st) != 0) {
        PS::S32 ret = -1;
        if(PS::Comm::getRank() == 0)
            ret = mkdir(dir_name, 0777);
        PS::Comm::broadcast(&ret, 1);
        if(ret == 0) {
            if(PS::Comm::getRank() == 0)
                fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
        } else {
            fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
            PS::Abort();
        }
    }
}


template<class Ttree, class Tfdust, class Tdis, class Tret, class Tfp,
         class Tsat, class Tpla>
void CalcForce(Ttree  & tree,
               Tfdust & force_dust_from_satellite,
               Tdis   func_dis,
               Tret   func_ret,
               PS::ParticleSystem<Tfp> & system,
               Tsat & satellite,
               const Tpla   & planet,
               PS::DomainInfo & dinfo,
               const PS::S32 n_walk_limit,
               const bool reuse_flag){
#if defined(SUNWAY) && defined(SUNWAY_FORCE_KERNEL)
        // In this case, we use Long's force kernel.
        //cpe_pars.energy_flag = 0;
        cpe_pars.energy_flag = 0;
        cpe_pars.first_flag = 1; // Half-kick
        cpe_pars.last_flag = 0; // Full-drift
#endif
    tree.calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe
        (func_dis, func_ret, 0,
         system, dinfo, n_walk_limit,
         false,  reuse_flag);
    
#ifdef FORCE_SUM_CHECK
    PS::F64vec  force_sum_glb = 0.0;
    force_sum_glb.x = PS::Comm::getSum(FORCE_SUM_LOC.x);
    force_sum_glb.y = PS::Comm::getSum(FORCE_SUM_LOC.y);
    force_sum_glb.z = PS::Comm::getSum(FORCE_SUM_LOC.z);
    if(PS::Comm::getRank() ==0 ){
        std::cerr<<"my_rank= "<<PS::Comm::getRank()
                 <<" force_sum_glb (on CPEs)= "<<force_sum_glb
                 <<std::endl;
    }
#endif
    
}


template<class Ttree, class Tfdust, class Tdis, class Tret, class Tfp,
         class Tsat, class Tpla>
void CalcForceLoopMerge(Ttree  & tree,
                        Tfdust & force_dust_from_satellite,
                        Tdis   func_dis,
                        Tret   func_ret,
                        PS::ParticleSystem<Tfp> & system,
                        Tsat & satellite,
                        const Tpla   & planet,
                        PS::DomainInfo & dinfo,
                        const PS::S32 n_walk_limit,
                        const PS::S32 n_loop_merge,
                        const bool initial=false){
    bool reuse_flag = false;
    if(initial){
        // initial step
#if defined(SUNWAY) && defined(SUNWAY_FORCE_KERNEL)
        // In this case, we use Long's force kernel.
        //cpe_pars.energy_flag = 0;
        cpe_pars.energy_flag = 0;
        cpe_pars.first_flag = 1; // Half-kick
        cpe_pars.last_flag = 0; // Full-drift
#endif

        ClearCounterAll(system, dinfo);
        
        tree.calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge
            (func_dis, func_ret, 0,
             system, dinfo, n_walk_limit,
             false,  reuse_flag);
        /*
        tree.calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe
            (func_dis, func_ret, 0,
             system, dinfo, n_walk_limit,
             false,  reuse_flag);
        */

        DumpTimingProfile(system, dinfo, fout_wtime);
        
#ifdef FORCE_SUM_CHECK
        PS::F64vec  force_sum_glb = 0.0;
        force_sum_glb.x = PS::Comm::getSum(FORCE_SUM_LOC.x);
        force_sum_glb.y = PS::Comm::getSum(FORCE_SUM_LOC.y);
        force_sum_glb.z = PS::Comm::getSum(FORCE_SUM_LOC.z);
        if(PS::Comm::getRank() ==0 ){
            std::cerr<<"my_rank= "<<PS::Comm::getRank()
                     <<" force_sum_glb (on CPEs)= "<<force_sum_glb
                     <<std::endl;
        }
#endif
        return;
    }
    // loop merge
    DISPATCH_MODE = 0;
    for(PS::S32 i=0; i<n_loop_merge; i++){
        if(i==n_loop_merge-1)  DISPATCH_MODE = 1;
        LapTimer::measure(i, LapTimer::BEGIN_OF_STEP);
        LapTimer::step_counter()=i;
#if defined(SUNWAY) && defined(SUNWAY_FORCE_KERNEL)
        // In this case, we use Long's force kernel.
        // We assume that func_dis is DispatchKernelSPKickDrift()
        cpe_pars.energy_flag = 0;
        cpe_pars.first_flag = 0; // Full-Kick
        cpe_pars.last_flag = 0; // Full-drift
#endif
        tree.calcForceAllAndWriteBackReuseListMultiWalkIndexLoopMerge
            (func_dis, func_ret, 0,
             system, dinfo, n_walk_limit,
             false,  reuse_flag);
#ifdef FORCE_SUM_CHECK
        PS::F64vec  force_sum_glb = 0.0;
        force_sum_glb.x = PS::Comm::getSum(FORCE_SUM_LOC.x);
        force_sum_glb.y = PS::Comm::getSum(FORCE_SUM_LOC.y);
        force_sum_glb.z = PS::Comm::getSum(FORCE_SUM_LOC.z);
        if(PS::Comm::getRank() == 0){
            std::cerr<<"my_rank= "<<PS::Comm::getRank()
                     <<" force_sum_glb (on CPEs)= "<<force_sum_glb
                     <<std::endl;
        }
#endif
        //nan_test((EPI*)EPI_POINTER, system.getNumberOfParticleLocal(), "nan-1");
        
        LapTimer::measure(i, LapTimer::END_OF_STEP);
        reuse_flag = true;
        DISPATCH_MODE = 0;


        #if 1
        //////////////
        // TIMING DUMP
        DumpTimingProfile(system, dinfo, fout_wtime);
        // TIMING DUMP
        //////////////
        #endif // #if 0
        //PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "timing dump ends (CalcForceLoopMerge)" << std::endl;
    }
    LapTimer::summary();
}



template<class T>
std::string ToS(const T & t){
    std::stringstream os;
    os<<std::setfill('0')<<std::setw(5)<<t;
    return os.str();
}

template<class Tsys>
void Rotate(Tsys & sys, const PS::S32 n, const PS::F64 theta){
    if(n <= 0) return;
    PS::F64 cth = cos(theta);
    PS::F64 sth = sin(theta);
    //std::cerr<<"sys[0].pos= "<<sys[0].pos<<std::endl;
    for(PS::S32 i=0; i<n; i++){
        const PS::F64 x_new  = cth*sys[i].pos.x - sth*sys[i].pos.y;
        const PS::F64 y_new  = sth*sys[i].pos.x + cth*sys[i].pos.y;
        const PS::F64 vx_new = cth*sys[i].vel.x - sth*sys[i].vel.y;
        const PS::F64 vy_new = sth*sys[i].vel.x + cth*sys[i].vel.y;
        sys[i].pos.x = x_new;
        sys[i].pos.y = y_new;
        sys[i].vel.x = vx_new;
        sys[i].vel.y = vy_new;
    }
    //std::cerr<<"sys[0].pos= "<<sys[0].pos<<std::endl;
    std::cerr<<std::endl;
}


inline void CalculateBoundaryOfDomain(const PS::S32 &np,
				      const PS::F64vec pos_sample[],
				      const PS::S32 cid,
				      const PS::S32 &istart,
				      const PS::S32 &iend,
				      const PS::F64ort pos_root_domain,
				      PS::F64 & xlow,
				      PS::F64 & xhigh) {
    if(istart == 0) {
	xlow  = pos_root_domain.low_[cid];
    } 
    else {
	xlow  = 0.5 * (pos_sample[istart-1][cid] + pos_sample[istart][cid]);
    }
    if(iend == np - 1) {
	xhigh = pos_root_domain.high_[cid];
    } 
    else {
	xhigh = 0.5 * (pos_sample[iend][cid] + pos_sample[iend+1][cid]);
    }
}

template<class Tsys>
inline void DomainDecision(PS::DomainInfo & dinfo,
			   const Tsys & system,
                           const PS::S32 n_smp_ave){
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    assert(n_proc%4 == 0);
    if(n_proc%4 == 0){
	int n_loc_tmp = system.getNumberOfParticleLocal();
        PS::S32 freq_smp = (n_loc_tmp < n_smp_ave) ? 1 : n_loc_tmp / n_smp_ave;
	int * n_recv_tmp = new int[n_proc];
	int * n_recv_disp_tmp = new int[n_proc+1];
	PS::F64vec * pos_loc_tmp = new PS::F64vec[n_loc_tmp];
        PS::F64 counter = 0.0;
	//for(int i=0; i<n_loc_tmp; i++){
        PS::S32 n_loc = 0;
        PS::S32 n_cnt = 0;
        for(n_loc =0; n_cnt<n_loc_tmp; n_loc++, n_cnt+=freq_smp){
	    pos_loc_tmp[n_loc].z = system[n_loc].pos.z;
	    if(system[n_loc].pos.x * system[n_loc].pos.y > 0.0){
		pos_loc_tmp[n_loc].x = fabs(system[n_loc].pos.x);
		pos_loc_tmp[n_loc].y = fabs(system[n_loc].pos.y);
	    }
	    else{
		pos_loc_tmp[n_loc].x = fabs(system[n_loc].pos.y);
		pos_loc_tmp[n_loc].y = fabs(system[n_loc].pos.x);
	    }
	}
	PS::Comm::allGather(&n_loc, 1, n_recv_tmp);
	n_recv_disp_tmp[0] = 0;
	for(int i=0; i<n_proc; i++){
	    n_recv_disp_tmp[i+1] = n_recv_disp_tmp[i] + n_recv_tmp[i];
	}
	int n_glb_tmp = n_recv_disp_tmp[n_proc];
	PS::F64vec * pos_glb_tmp = new PS::F64vec[n_glb_tmp];
	PS::Comm::allGatherV(pos_loc_tmp, n_loc, pos_glb_tmp, n_recv_tmp, n_recv_disp_tmp);
	PS::S32 n_proc_quat = n_proc / 4;
	fout_wtime<<"n_proc, n_proc_quat="<<n_proc<<" "<<n_proc_quat<<std::endl;
	PS::S32 nx_quat = sqrt((PS::F64)n_proc_quat-0.000001)+1;
	while( n_proc_quat % nx_quat != 0) nx_quat++;
	PS::S32 ny_quat = n_proc_quat / nx_quat;
	PS::S32 nx = nx_quat*2;
	PS::S32 ny = ny_quat*2;
	PS::S32 nz = 1;
	if(PS::Comm::getRank() == 0){
            std::cout<<"nx_quat, ny_quat, nx, ny, nz= "
                     <<nx_quat<<" "<<ny_quat<<" "<<nx<<" "<<ny<<" "<<nz<<std::endl;            
            /*
	    fout_wtime<<"nx_quat, ny_quat, nx, ny, nz= "
		      <<nx_quat<<" "<<ny_quat<<" "<<nx<<" "<<ny<<" "<<nz<<std::endl;
            */
	}
	dinfo.setNumberOfDomainMultiDimension(nx, ny, nz);
	PS::F64ort * pos_domain_tmp = new PS::F64ort[n_proc];
	if(PS::Comm::getRank() == 0){
	    fout_wtime<<"n_glb_tmp="<<n_glb_tmp<<std::endl;
	    PS::S32 * istart = new PS::S32[n_proc_quat];
	    PS::S32 * iend   = new PS::S32[n_proc_quat];
	    PS::F64ort pos_root_domain_tmp = PS::F64ort( PS::F64vec(0.0, 0.0, -PS::LARGE_FLOAT), PS::F64vec(PS::LARGE_FLOAT, PS::LARGE_FLOAT, PS::LARGE_FLOAT));
            std::sort(pos_glb_tmp, pos_glb_tmp+n_glb_tmp, PS::LessOPX());
	    fout_wtime<<"pos_glb_tmp[n_glb_tmp-1].x="<<pos_glb_tmp[n_glb_tmp-1].x<<std::endl;
	    for(PS::S32 i = 0; i < n_proc_quat; i++) {
		istart[i] = ((PS::S64)(i) * (PS::S64)(n_glb_tmp)) / (PS::S64)(n_proc_quat);
		if(i > 0) iend[i-1] = istart[i] - 1;
	    }
	    iend[n_proc_quat-1] = n_glb_tmp - 1;
	    for(PS::S32 ix = 0; ix<nx_quat; ix++) {
		PS::S32 ix0 =  ix      * ny_quat * nz;
		PS::S32 ix1 = (ix + 1) * ny_quat * nz;
		PS::F64 x0, x1;
		CalculateBoundaryOfDomain(n_glb_tmp, pos_glb_tmp, 0, istart[ix0], iend[ix1-1], pos_root_domain_tmp, x0, x1);
		for(PS::S32 i=0; i<ny_quat*nz; i++) {
		    PS::S32 offset = (nx_quat+ix)*ny;
		    pos_domain_tmp[offset+i].low_[0]  = x0;
		    pos_domain_tmp[offset+i].high_[0] = x1;
		    pos_domain_tmp[offset+ny_quat+i].low_[0]  = x0;
		    pos_domain_tmp[offset+ny_quat+i].high_[0] = x1;

		    offset = (nx_quat-1-ix)*ny;
		    pos_domain_tmp[offset+i].low_[0]  = -x1;
		    pos_domain_tmp[offset+i].high_[0] = -x0;
		    pos_domain_tmp[offset+ny_quat+i].low_[0]  = -x1;
		    pos_domain_tmp[offset+ny_quat+i].high_[0] = -x0;
		}
	    }

	    for(PS::S32 ix = 0; ix<nx_quat; ix++) {
		PS::S32 ix0 =  ix      * ny_quat * nz;
		PS::S32 ix1 = (ix + 1) * ny_quat * nz;
                std::sort(pos_glb_tmp+istart[ix0], pos_glb_tmp+iend[ix1-1]+1, PS::LessOPY());
		PS::S32 n_tmp_y = iend[ix1-1] - istart[ix0] + 1;
		for(PS::S32 iy = 0; iy<ny_quat; iy++) {
		    PS::S32 iy0 = ix0 +  iy      * nz;
		    PS::S32 iy1 = ix0 + (iy + 1) * nz;
		    PS::F64 y0, y1;
		    CalculateBoundaryOfDomain(n_tmp_y, pos_glb_tmp+istart[ix0], 1, istart[iy0]-istart[ix0], iend[iy1-1]-istart[ix0], pos_root_domain_tmp, y0, y1);
		    for(PS::S32 i=0; i<nz; i++) {
			PS::S32 offset = (nx_quat+ix)*ny + ny_quat + iy;
			pos_domain_tmp[offset+i].low_[1]  = y0;
			pos_domain_tmp[offset+i].high_[1] = y1;
			offset = (nx_quat-ix-1)*ny + ny_quat + iy;
			pos_domain_tmp[offset+i].low_[1]  = y0;
			pos_domain_tmp[offset+i].high_[1] = y1;
			offset = (nx_quat+ix)*ny + ny_quat-iy-1;
			pos_domain_tmp[offset+i].low_[1]  = -y1;
			pos_domain_tmp[offset+i].high_[1] = -y0;
			offset = (nx_quat-ix-1)*ny + ny_quat-iy-1;
			pos_domain_tmp[offset+i].low_[1]  = -y1;
			pos_domain_tmp[offset+i].high_[1] = -y0;
		    }
		}
	    }
	    delete [] istart;
	    delete [] iend;
	}
	PS::Comm::broadcast(pos_domain_tmp, n_proc);
	for(PS::S32 i=0; i<n_proc; i++){
	    pos_domain_tmp[i].low_.z  = -PS::LARGE_FLOAT;
	    pos_domain_tmp[i].high_.z =  PS::LARGE_FLOAT;
	    dinfo.setPosDomain(i, pos_domain_tmp[i]);
	}
	delete [] n_recv_tmp;
	delete [] n_recv_disp_tmp;
	delete [] pos_loc_tmp;
	delete [] pos_glb_tmp;
	delete [] pos_domain_tmp;
    }
}

#ifdef PHI_R_TREE

    template<class Tepi, class Tfp>
        void CopyEpiToFpWithCoordinateTransform(Tfp * fp, const Tepi * epi, const PS::S32 n_loc){
        static const PS::F64 PI = 4.0*atan(1.0);
        for(PS::S32 i=0; i<n_loc; i++){
            const PS::F64 x = epi[i].pos.x;
            const PS::F64 y = epi[i].pos.y;
            const PS::F64 z = epi[i].pos.z;
            //if(x==0.0 || y==0.0 || z==0.0)std::cerr<<"x= "<<x<<" y= "<<y<<" z= "<<z<<std::endl;
            PS::F64 phi = atan2(y, x);
            if(phi < 0.0) phi += 2.0*PI;
            else if(phi >= 2.0*PI) phi -= 2.0*PI;
            const PS::F64 r = sqrt(x*x+y*y);
            fp[i].pos  = PS::F64vec(phi, r, z);
            fp[i].vel  = epi[i].vel;
            fp[i].id   = epi[i].id;
            fp[i].mass = epi[i].mass;
        }
    }

    template<class Tfp>
    void CoordinateChangeCylindricalToCartesian(Tfp * fp){
        const PS::F64 phi = fp->pos.x;
        const PS::F64 r = fp->pos.y;
        fp->pos.x = r*cos(phi);
        fp->pos.y = r*sin(phi);
    }
    template<class Tfp>
    void CoordinateChangeCartesianToCylindrical(Tfp * fp){
        static const PS::F64 PI = 4.0*atan(1.0);
        const PS::F64 x = fp->pos.x;
        const PS::F64 y = fp->pos.y;
        const PS::F64 z = fp->pos.z;
        PS::F64 phi = atan2(y, x);
        phi = (phi < 0.0) ? 2.0*PI+phi : phi;
        const PS::F64 r = sqrt(x*x+y*y);
        fp->pos  = PS::F64vec(phi, r, z);
    }
#endif

template<class Tptcl>
void CalcCM(const Tptcl ptcl[], const PS::S32 n,
            PS::F64 & cm_mass_glb, PS::F64vec & cm_pos_glb){
    PS::F64 cm_mass_loc = 0.0;
    PS::F64vec cm_pos_loc = 0.0;
    for(PS::S32 i=0; i<n; i++){
        cm_mass_loc += ptcl[i].mass;
        cm_pos_loc += ptcl[i].mass * ptcl[i].pos;
    }
    cm_mass_glb = PS::Comm::getSum(cm_mass_loc);
    cm_pos_glb = PS::Comm::getSum(cm_pos_loc);
    cm_pos_glb /= cm_mass_glb;
}

template<class Tptcl>
void CalcCMCyl(const Tptcl ptcl[], const PS::S32 n,
               PS::F64 & cm_mass_glb, PS::F64vec & cm_pos_glb){
    PS::F64 cm_mass_loc = 0.0;
    PS::F64vec cm_pos_loc = 0.0;
    for(PS::S32 i=0; i<n; i++){
        cm_mass_loc += ptcl[i].mass;
        cm_pos_loc.x += ptcl[i].mass * ptcl[i].pos.y*cos(ptcl[i].pos.x);
        cm_pos_loc.y += ptcl[i].mass * ptcl[i].pos.y*sin(ptcl[i].pos.x);
        cm_pos_loc.z += ptcl[i].mass * ptcl[i].pos.z;
    }
    cm_mass_glb = PS::Comm::getSum(cm_mass_loc);
    cm_pos_glb = PS::Comm::getSum(cm_pos_loc);
    cm_pos_glb /= cm_mass_glb;
}

template<class Tsys>
void WriteSnp(const Tsys & system_grav,
              PS::S32 & snp_id,
              const char * dir_name,
              const PS::S32 n_out,
              const std::string comment=""){
    PS::S32 n_loc = system_grav.getNumberOfParticleLocal();
    std::ofstream fout_snp;
    std::string s_snp;
    s_snp = dir_name + ("/snp_" + ToS(snp_id) + "_" + ToS(PS::Comm::getRank()) + ".dat");
    std::cerr<<"s_snp: "<<s_snp<<std::endl;
    fout_snp.open(s_snp.c_str());
    if(!fout_snp.is_open()) exit(1);
    fout_snp<<std::setprecision(15);
    //PS::Comm::barrier();
    fout_snp<<"# n_loc= "<<n_loc<<" N_SATELLITE= "<<N_SATELLITE<<comment<<std::endl;
    PS::S32 interval = 1;
    if(n_out > 0 && n_loc > n_out){
        interval = n_loc / n_out;
    }
    for(PS::S32 i=0; i<N_SATELLITE; i++){
        fout_snp<<i<<"   "<<SATELLITE[i].id<<"   "<<SATELLITE[i].pos<<"   "<<SATELLITE[i].vel<<std::endl;
    }
    for(PS::S32 i=0; i<n_loc; i+=interval){
        fout_snp<<i<<"   "<<system_grav[i].id<<"   "<<system_grav[i].pos<<"   "<<system_grav[i].vel<<std::endl;
    }

    fout_snp.close();
    snp_id++;
    //PS::Comm::barrier();
}

template<class Tpi, class Tpj>
void CalcPotGrav(const Tpi & pi,
                 const Tpj & pj,
                 const PS::F64 r_nat,
                 PS::F64 & eng){
    PS::F64vec rij = pi.pos - pj.pos;
    PS::F64 r_sq = rij*rij;
    PS::F64 r = sqrt(r_sq);
#if 1
    eng -= pi.mass*pj.mass/r;
#else
    PS::F64 r_nat_inv = 1.0 / r_nat;
    PS::F64 r_nat_sq = r_nat * r_nat;
    if(r <= r_nat){
        //eng -= 0.5 * pi.mass * pj.mass * r_nat_inv * r_nat_inv * r_nat_inv * (r-r_nat) * (r-r_nat);
        eng -= (-0.5*pi.mass*pj.mass*r_nat_inv*r_nat_inv*r_nat_inv*(r_sq-r_nat_sq)) + (pi.mass*pj.mass*r_nat_inv);
    }
    else{
        //eng -= pi.mass * pj.mass * (1.0/r - 1.0/r_nat); // 1.0/r_nat is an offset to make contegious potential
        eng -= pi.mass*pj.mass/r;
    }
#endif
}

template<class Tpi, class Tpj>
void CalcPotSprg(const Tpi & pi,
                 const Tpj & pj,
                 const PS::F64 r_nat, // natural lentgh
                 const PS::F64 kappa, // k / m_red
                 PS::F64 & eng){
    PS::F64vec rij = pi.pos - pj.pos;
    PS::F64 r_sq = rij*rij;
    PS::F64 r = sqrt(r_sq);
    PS::F64 dl = r - r_nat;
    if(dl >= 0.0) return;
    PS::F64 m_red = (pi.mass*pj.mass) / (pi.mass+pj.mass);
    PS::F64 k = kappa * m_red;
    eng += 0.5*k*dl*dl;
}

template<class Tptcl, class Tsat>
void CalcEnergy(const PS::ParticleSystem<Tptcl> & ptcl,
                const Tsat * sat,
                const PS::S32 n_sat,
                const PS::F64 r_coll,
                const PS::F64 kappa,
                PS::F64 & eng_tot_glb,
                PS::F64 & eng_kin_glb,
                PS::F64 & eng_pot_grav_glb,
                PS::F64 & eng_pot_sprg_glb){
    eng_tot_glb = eng_kin_glb = eng_pot_grav_glb = eng_pot_sprg_glb = 0.0;
    Tptcl * ptcl_glb = NULL;
    PS::S64 n_ptcl_glb = 0;
    ptcl.allGatherParticle(ptcl_glb, n_ptcl_glb);
    assert(n_ptcl_glb == ptcl.getNumberOfParticleGlobal());

    PS::F64 eng_kin_loc = 0.0;
    PS::F64 eng_pot_grav_loc = 0.0;
    PS::F64 eng_pot_sprg_loc = 0.0;
    Tptcl planet;
    planet.mass = 1.0;
    planet.pos = 0.0;
    const PS::S32 n_ptcl_loc = ptcl.getNumberOfParticleLocal();

    for(PS::S32 i=0; i<n_ptcl_loc; i++){
        PS::F64 e_tmp_0 = 0.0;
        PS::F64 e_sprg_tmp_0 = 0.0;
        for(PS::S32 j=0; j<n_ptcl_glb; j++){
            if(ptcl[i].id == ptcl_glb[j].id) continue;
            CalcPotGrav(ptcl[i], ptcl_glb[j], r_coll*2.0, e_tmp_0);
            CalcPotSprg(ptcl[i], ptcl_glb[j], r_coll*2.0, kappa, e_sprg_tmp_0);
        }
        e_tmp_0 *= 0.5;
        e_sprg_tmp_0 *= 0.5;
        PS::F64 e_tmp_1 = 0.0;
        PS::F64 e_sprg_tmp_1 = 0.0;
        for(PS::S32 j=0; j<n_sat; j++){
            CalcPotGrav(ptcl[i], sat[j], r_coll*2.0, e_tmp_1);
            CalcPotSprg(ptcl[i], sat[j], r_coll*2.0, kappa, e_sprg_tmp_1);
        }
        eng_pot_grav_loc += (e_tmp_0 + e_tmp_1);
        eng_pot_sprg_loc += (e_sprg_tmp_0 + e_sprg_tmp_1);
    }
    
    if(PS::Comm::getRank() == 0){
        PS::F64 e_tmp_2 = 0.0;
        PS::F64 e_sprg_tmp_2 = 0.0;
        for(PS::S32 i=0; i<n_sat; i++){
            for(PS::S32 j=i+1; j<n_sat; j++){
                CalcPotGrav(sat[i], sat[j], r_coll*2.0, e_tmp_2);
                CalcPotSprg(sat[i], sat[j], r_coll*2.0, kappa, e_sprg_tmp_2);
            }
        }
        eng_pot_grav_loc += e_tmp_2;
        eng_pot_sprg_loc += e_sprg_tmp_2;
    }
    
    if(PS::Comm::getRank() == 0){
        PS::F64 e_tmp_3 = 0.0;
        for(PS::S32 i=0; i<n_sat; i++){
            eng_kin_loc += 0.5 * sat[i].mass * sat[i].vel * sat[i].vel;
            CalcPotGrav(sat[i], planet, r_coll*2.0, e_tmp_3);
        }
        eng_pot_grav_loc += e_tmp_3;
    }

    PS::F64 e_tmp_4 = 0.0;
    for(PS::S32 i=0; i<n_ptcl_loc; i++){
        eng_kin_loc += 0.5 * ptcl[i].mass * ptcl[i].vel * ptcl[i].vel;
        CalcPotGrav(ptcl[i], planet, r_coll*2.0, eng_pot_grav_loc);
    }
    eng_pot_grav_loc += e_tmp_4;
    
    eng_pot_grav_glb = PS::Comm::getSum(eng_pot_grav_loc);
    eng_pot_sprg_glb = PS::Comm::getSum(eng_pot_sprg_loc);
    eng_kin_glb = PS::Comm::getSum(eng_kin_loc);
    eng_tot_glb = eng_pot_grav_glb + eng_pot_sprg_glb + eng_kin_glb;

    delete [] ptcl_glb;
}


template<class Tpi, class Tpj, class Tfi>
void CalcForceImpl(const Tpi & pi,
                   const Tpj & pj,
                   const PS::F64 r_nat,
                   Tfi & fi){
    PS::F64vec rij = pi.pos -pj.pos;
    PS::F64    r_sq = rij*rij;
    PS::F64    r_inv = 1.0 / sqrt(r_sq);
    PS::F64    pot  = pj.mass * r_inv;
    fi.acc -= pot * r_inv * r_inv * rij;
}

template<class Tptcl, class Tsat, class Tforce>
void CalcForceDirect(const PS::ParticleSystem<Tptcl> & ptcl,
                     const Tsat * sat,
                     const PS::S32 n_sat,
                     const PS::F64 r_coll,
                     const PS::F64 kappa,
                     Tforce force[]){
    Tptcl * ptcl_glb = NULL;
    PS::S64 n_ptcl_glb = 0;
    ptcl.allGatherParticle(ptcl_glb, n_ptcl_glb);
    assert(n_ptcl_glb == ptcl.getNumberOfParticleGlobal());
        
    Tptcl planet;
    planet.mass = 1.0;
    planet.pos = 0.0;
    const PS::S32 n_ptcl_loc = ptcl.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_ptcl_loc; i++){
        force[i].clear();
        CalcForceImpl(ptcl[i], planet, r_coll*2.0, force[i]);
        for(PS::S32 j=0; j<n_ptcl_glb; j++){
            if(ptcl[i].id == ptcl_glb[j].id) continue;
            CalcForceImpl(ptcl[i], ptcl_glb[j], r_coll*2.0, force[i]);
            //CalcPotSprg(ptcl[i], ptcl_glb[j], r_coll*2.0, kappa, eng_pot_sprg_loc);
        }
        for(PS::S32 j=0; j<n_sat; j++){
            CalcForceImpl(ptcl[i], sat[j], r_coll*2.0, force[i]);
            //CalcPotSprg(ptcl[i], sat[j], r_coll*2.0, kappa, eng_pot_sprg_loc);
        }
    }
    for(PS::S32 i=0; i<n_sat; i++){
        CalcForceImpl(sat[i], planet, r_coll*2.0, force[i+n_ptcl_loc]);
        for(PS::S32 j=0; j<n_ptcl_glb; j++){
            CalcForceImpl(sat[i], ptcl_glb[j], r_coll*2.0, force[i+n_ptcl_loc]);
            //CalcPotSprg(sat[i], ptcl_glb[j], r_coll*2.0, kappa, eng_pot_sprg_loc);
        }
        for(PS::S32 j=i+1; j<n_sat; j++){
            CalcForceImpl(sat[i], sat[j], r_coll*2.0, force[i+n_ptcl_loc]);
            //CalcPotSprg(sat[i], sat[j], r_coll*2.0, kappa, eng_pot_sprg_loc);
        }
    }
    
    delete [] ptcl_glb;
}


void DivideNProc(PS::S32 & nx,
                 PS::S32 & ny,
                 PS::S32 & nz,
                 const PS::S32 n_proc,
                 const PS::F64 delta_ax){
#ifdef PHI_R_TREE
    PS::F64 PI = 4.0*atan(1.0);
    /*
    PS::S32 ny = 2;
    PS::S32 nz = 1;
    PS::S32 nx = n_proc / ny;
    */
    nz = 1;
    ny = 1;
    nx = n_proc / ny;
    double dx = 2.0*PI   / nx;
    double dy = delta_ax / ny;
    double ratio = (dx < dy) ? dx / dy : dy / dx;
    PS::S32 ny_tmp = ny;
    PS::S32 nx_tmp = nx;
    double dx_tmp = dx;
    double dy_tmp = dy;
    double ratio_tmp = ratio;
    do{
        ny = ny_tmp;
        nx = nx_tmp;
        ratio = ratio_tmp;
        ny_tmp += 1;
        while( n_proc % ny_tmp != 0) ny_tmp++;
        nx_tmp = n_proc / ny_tmp;
        dx_tmp = 2.0*PI   / nx_tmp;
        dy_tmp = delta_ax / ny_tmp;
        ratio_tmp = (dx_tmp < dy_tmp) ? dx_tmp / dy_tmp : dy_tmp / dx_tmp;
    }while( fabs(ratio_tmp-1.0) < fabs(ratio-1.0));

    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0){
        std::cerr<<"nx= "<<nx<<" ny= "<<ny<<" nz= "<<nz<<std::endl;
        fout_wtime<<"nx= "<<nx<<" ny= "<<ny<<" nz= "<<nz<<std::endl;
    }

#elif 1
    PS::S32 nx = sqrt((PS::F64)n_proc-0.000001)+1;
    while( n_proc % nx != 0) nx++;
    PS::S32 ny = n_proc / nx;
    PS::S32 nz = 1;
#else
    PS::S32 nx = 400, ny = 396, nz = 1; 
#endif
}

template<class Tptcl>
PS::F64 CalcZCoordAbsMax(const Tptcl ptcl[],
                         const PS::S32 n){
    // return also negative value
    PS::F64 z_coord_offset_loc = 0.0;
    for(PS::S32 i=0; i<n; i++){
        if( fabs(z_coord_offset_loc) < fabs(ptcl[i].pos.z) ) z_coord_offset_loc = ptcl[i].pos.z;
    }
    PS::F64 z_coord_max_glb = PS::Comm::getMaxValue(z_coord_offset_loc);
    PS::F64 z_coord_min_glb = PS::Comm::getMinValue(z_coord_offset_loc);
    PS::F64 ret = ( fabs(z_coord_max_glb) > fabs(z_coord_min_glb) ) ? z_coord_max_glb : z_coord_min_glb;
    return ret;
}

template<class Tptcl>
PS::F64 CalcZCoordShift(const Tptcl ptcl[],
                        const PS::S32 n,
                        const PS::F64 z_offset){
    PS::F64 z = CalcZCoordAbsMax(ptcl, n);
    if(z > 0) z += z_offset;
    else z -= z_offset;
    return z;
}

template<class Tptcl>
void RotateAndShiftPtcl(const Tptcl * ptcl,
                        const PS::S32 n,
                        const PS::F64 theta_rot,
                        const PS::F64 z_offset,
                        PS::F64 & z_coord_shift){
    const PS::F64 cth = std::cos(theta_rot);
    const PS::F64 sth = std::sin(theta_rot);
    PS::F64 z_shift_ptcl_ar[64];

    unsigned long args[5];
    //** Rotate the particles
    args[0] = (unsigned long) n;
    args[1] = (unsigned long) ptcl;
    args[2] = (unsigned long) &cth;
    args[3] = (unsigned long) &sth;
    args[4] = (unsigned long) z_shift_ptcl_ar;
    __real_athread_spawn((void*)slave_RotateAndShiftZ, args);
    athread_join();
    PS::F64 z_shift_sat_ar[64];
    //** Rotate the satellites
    args[0] = (unsigned long) N_SATELLITE;
    args[1] = (unsigned long) SATELLITE.getPointer();
    args[2] = (unsigned long) &cth;
    args[3] = (unsigned long) &sth;
    args[4] = (unsigned long) z_shift_sat_ar;
    __real_athread_spawn((void*)slave_RotateAndShiftZ, args);
    athread_join();
    PS::F64 z_coord_offset_loc = 0.0;
    for(PS::S32 i=0; i<64; i++){
        if(fabs(z_shift_ptcl_ar[i]) > fabs(z_coord_offset_loc) ){
            z_coord_offset_loc = z_shift_ptcl_ar[i];
        }
        if(fabs(z_shift_sat_ar[i]) > fabs(z_coord_offset_loc) ){
            z_coord_offset_loc = z_shift_sat_ar[i];
        }
    }
    PS::F64 z_shift_max_tmp = PS::Comm::getMaxValue(z_coord_offset_loc);
    PS::F64 z_shift_min_tmp = PS::Comm::getMinValue(z_coord_offset_loc);
    z_coord_shift = ( fabs(z_shift_max_tmp) > fabs(z_shift_min_tmp) ) ? z_shift_max_tmp : z_shift_min_tmp;
    if(z_coord_shift > 0.0) z_coord_shift += z_offset;
    else z_coord_shift -= z_offset;
}

template<class Tptcl>
void SetPtclId(Tptcl ptcl[],
               const PS::S32 n_loc){
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    PS::S32 * n_loc_array = new PS::S32[n_proc];
    PS::S32 n_send = n_loc;
    MPI_Allgather(&n_send, 1, MPI_INT, n_loc_array, 1, MPI_INT, MPI_COMM_WORLD);
    PS::S64 id_offset = 0;
    PS::S32 my_rank = PS::Comm::getRank();
    for(PS::S32 i=0; i<my_rank; i++){
        id_offset += (PS::S64)(n_loc_array[i]);
    }
    for(PS::S32 i=0; i<n_loc; i++){
        ptcl[i].id = id_offset + (PS::S64)(i);
    }
    delete [] n_loc_array;
}

template<class Tptcl>
PS::F64vec CalcForceSumDirect(PS::ParticleSystem<Tptcl>  & sys){
    Tptcl * jp = NULL;
    PS::S64 nj = 0;
    PS::Comm::barrier();
    std::cerr<<"before allgather ptcls"<<std::endl;
    sys.allGatherParticle(jp, nj);
    PS::Comm::barrier();
    std::cerr<<"after allgather ptcls"<<std::endl;
    PS::F64 cm_mass = 0.0;
    PS::F64vec cm_pos = 0.0;
    for(PS::S32 i=0; i<nj; i++){
        cm_mass += jp[i].mass;
        cm_pos  += jp[i].mass*jp[i].pos;
    }
    cm_pos /= cm_mass;
    std::cerr<<"my_rank= "<<PS::Comm::getRank()
             <<" cm_mass= "<<cm_mass
             <<" cm_pos= "<<cm_pos
             <<std::endl;
    PS::S64 ni = sys.getNumberOfParticleLocal();
    PS::F64vec force_sum_loc = 0.0;

    const PS::S32 n_seg = 8;
    PS::F64 seg_acc[n_seg+1];
    seg_acc[0] = 1e-30;
    seg_acc[n_seg] = 1e30;
    seg_acc[1] = 1e-10;
    for(PS::S32 i=2; i<n_seg; i++){
        seg_acc[i] = seg_acc[i-1]*10.0;
    }
    for(PS::S64 seg=1; seg<n_seg+1; seg++){    
        for(PS::S64 i=0; i<ni; i++){
            if(PS::Comm::getRank()==0){
                if(i % ((ni/100)+1) == 0){
                    std::cerr<<"i= "<<i<<" ni= "<<ni<<" nj= "<<nj<<std::endl;
                }
            }
            for(PS::S64 j=0; j<nj; j++){
                PS::F64vec acc_tmp = 0.0;
                if(sys[i].id==jp[j].id) continue;
                PS::CalcForce2Body(sys[i].pos, jp[j].pos, jp[j].mass,
                                   acc_tmp);
                if(sqrt(acc_tmp*acc_tmp)>=seg_acc[seg-1] && sqrt(acc_tmp*acc_tmp)<seg_acc[seg]){
                    force_sum_loc += acc_tmp;
                }
            }
            if(PS::Comm::getRank()==8){
                std::cerr<<" seg= "<<seg
                         <<" seg_acc= "<<seg_acc[seg]
                         <<" force_sum_loc= "<<force_sum_loc<<std::endl;
            }
        }
    }
    std::cerr<<"my_rank= "<<PS::Comm::getRank()
             <<" force_sum_loc= "<<force_sum_loc
             <<std::endl;
    PS::F64vec force_sum_glb = PS::Comm::getSum(force_sum_loc);
    delete [] jp;
    return force_sum_loc;
}
