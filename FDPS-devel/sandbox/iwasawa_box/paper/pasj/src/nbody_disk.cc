#include<iostream>
#include<fstream>
#include<unistd.h>
#include<particle_simulator.hpp>
#include"class.hpp"

#ifdef USEPHANTOMGRAPE_X86
#include"phantomquad_x86.hpp"
#endif
#ifdef USEPHANTOMGRAPE_K
#include"phantomquad.hpp"
#endif

void DumpTimeProfile(const PS::TimeProfile & tp, const PS::TimeProfile & tp_max, const PS::S32 rank_max[], std::ostream & fout){
    PS::S32 id = 0;
    fout<<"collect_sample_particle= "<<tp.collect_sample_particle<<", max= "<<tp_max.collect_sample_particle<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"decompose_domain= "<<tp.decompose_domain<<", max= "<<tp_max.decompose_domain<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"exchange_particle= "<<tp.exchange_particle<<", max= "<<tp_max.exchange_particle<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"set_particle_local_tree= "<<tp.set_particle_local_tree<<", max= "<<tp_max.set_particle_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"set_root_cell= "<<tp.set_root_cell<<", max= "<<tp_max.set_root_cell<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_local_tree= "<<tp.make_local_tree<<", max= "<<tp_max.make_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"calc_moment_local_tree= "<<tp.calc_moment_local_tree<<", max= "<<tp_max.calc_moment_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_LET_1st= "<<tp.make_LET_1st<<", max= "<<tp_max.make_LET_1st<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"exchange_LET_1st= "<<tp.exchange_LET_1st<<", max= "<<tp_max.exchange_LET_1st<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_LET_2nd= "<<tp.make_LET_2nd<<", max= "<<tp_max.make_LET_2nd<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"exchange_LET_2nd= "<<tp.exchange_LET_2nd<<", max= "<<tp_max.exchange_LET_2nd<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"set_particle_global_tree= "<<tp.set_particle_global_tree<<", max= "<<tp_max.set_particle_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_global_tree= "<<tp.make_global_tree<<", max= "<<tp_max.make_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"calc_moment_global_tree= "<<tp.calc_moment_global_tree<<", max= "<<tp_max.calc_moment_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"calc_force= "<<tp.calc_force<<", max= "<<tp_max.calc_force<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<std::endl;
} 

void DumpTimeProfile0(const PS::TimeProfile & tp, const PS::TimeProfile & tp_max, const PS::S32 rank_max[], std::ostream & fout, const PS::S64 n_smp_int=1){
    PS::S32 id = 0;
    fout<<"collect_sample_particle= "<<tp.collect_sample_particle / n_smp_int<<", max= "<<tp_max.collect_sample_particle / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"decompose_domain= "<<tp.decompose_domain / n_smp_int<<", max= "<<tp_max.decompose_domain / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;


    fout<<"tp.decompose_domain__gather_particle= "<<tp.decompose_domain__gather_particle / n_smp_int<<", max= "<<tp_max.decompose_domain__gather_particle / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"tp.decompose_domain__sort_particle_1st= "<<tp.decompose_domain__sort_particle_1st / n_smp_int<<", max= "<<tp_max.decompose_domain__sort_particle_1st / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"tp.decompose_domain__sort_particle_2nd= "<<tp.decompose_domain__sort_particle_2nd / n_smp_int<<", max= "<<tp_max.decompose_domain__sort_particle_2nd / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"tp.decompose_domain__sort_particle_3rd= "<<tp.decompose_domain__sort_particle_3rd / n_smp_int<<", max= "<<tp_max.decompose_domain__sort_particle_3rd / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;


    fout<<"tp.decompose_domain__setup= "<<tp.decompose_domain__setup / n_smp_int<<", max= "<<tp_max.decompose_domain__setup / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"tp.decompose_domain__determine_coord_1st= "<<tp.decompose_domain__determine_coord_1st / n_smp_int<<", max= "<<tp_max.decompose_domain__determine_coord_1st / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"tp.decompose_domain__migrae_particle_1st= "<<tp.decompose_domain__migrae_particle_1st / n_smp_int<<", max= "<<tp_max.decompose_domain__migrae_particle_1st / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"tp.decompose_domain__determine_coord_2nd= "<<tp.decompose_domain__determine_coord_2nd / n_smp_int<<", max= "<<tp_max.decompose_domain__determine_coord_2nd / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"tp.decompose_domain__determine_coord_3rd= "<<tp.decompose_domain__determine_coord_3rd / n_smp_int<<", max= "<<tp_max.decompose_domain__determine_coord_3rd / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"tp.decompose_domain__exchange_pos_domain= "<<tp.decompose_domain__exchange_pos_domain / n_smp_int<<", max= "<<tp_max.decompose_domain__exchange_pos_domain / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;

    fout<<std::endl;
} 

void DumpTimeProfile1(const PS::TimeProfile & tp, const PS::TimeProfile & tp_max, const PS::S32 rank_max[], std::ostream & fout, const PS::S64 n_smp_int=1){
    //PS::S32 id = 6;
    PS::S32 id = 12;
    fout<<"exchange_particle= "<<tp.exchange_particle / n_smp_int<<", max= "<<tp_max.exchange_particle / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"exchange_particle__find_particle= "<<tp.exchange_particle__find_particle / n_smp_int
        <<", max= "<<tp_max.exchange_particle__find_particle / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"exchange_particle__exchange_particle= "<<tp.exchange_particle__exchange_particle / n_smp_int
        <<", max= "<<tp_max.exchange_particle__exchange_particle / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<std::endl;
} 

void DumpTimeProfile2(const PS::TimeProfile & tp, const PS::TimeProfile & tp_max, const PS::S32 rank_max[], std::ostream & fout, const PS::S64 n_smp_int=1){
    //PS::S32 id = 7;
    //PS::S32 id = 13;
    PS::S32 id = 15;
    fout<<"set_particle_local_tree= "<<tp.set_particle_local_tree / n_smp_int<<", max= "<<tp_max.set_particle_local_tree / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"set_root_cell= "<<tp.set_root_cell / n_smp_int<<", max= "<<tp_max.set_root_cell / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_local_tree= "<<tp.make_local_tree / n_smp_int<<", max= "<<tp_max.make_local_tree / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"calc_moment_local_tree= "<<tp.calc_moment_local_tree / n_smp_int<<", max= "<<tp_max.calc_moment_local_tree / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_LET_1st= "<<tp.make_LET_1st / n_smp_int<<", max= "<<tp_max.make_LET_1st / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"exchange_LET_1st= "<<tp.exchange_LET_1st / n_smp_int<<", max= "<<tp_max.exchange_LET_1st / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_LET_2nd= "<<tp.make_LET_2nd / n_smp_int<<", max= "<<tp_max.make_LET_2nd / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"exchange_LET_2nd= "<<tp.exchange_LET_2nd / n_smp_int<<", max= "<<tp_max.exchange_LET_2nd / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"set_particle_global_tree= "<<tp.set_particle_global_tree / n_smp_int<<", max= "<<tp_max.set_particle_global_tree / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"make_global_tree= "<<tp.make_global_tree / n_smp_int<<", max= "<<tp_max.make_global_tree / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"calc_moment_global_tree= "<<tp.calc_moment_global_tree / n_smp_int<<", max= "<<tp_max.calc_moment_global_tree / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"calc_force= "<<tp.calc_force / n_smp_int<<", max= "<<tp_max.calc_force / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;

    fout<<"  make_local_tree_tot= "<<tp.make_local_tree_tot / n_smp_int<<", max= "<<tp_max.make_local_tree_tot / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"  make_global_tree_tot= "<<tp.make_global_tree_tot / n_smp_int<<", max= "<<tp_max.make_global_tree_tot / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"  exchange_LET_tot= "<<tp.exchange_LET_tot / n_smp_int<<", max= "<<tp_max.exchange_LET_tot / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;

    fout<<"  morton_sort_local_tree= "<<tp.morton_sort_local_tree / n_smp_int<<", max= "<<tp_max.morton_sort_local_tree / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"  link_cell_local_tree= "<<tp.link_cell_local_tree / n_smp_int<<", max= "<<tp_max.link_cell_local_tree / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"  morton_sort_global_tree= "<<tp.morton_sort_global_tree / n_smp_int<<", max= "<<tp_max.morton_sort_global_tree / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"  link_cell_global_tree= "<<tp.link_cell_global_tree / n_smp_int<<", max= "<<tp_max.link_cell_global_tree / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;

    fout<<"  exchange_LET_1st__a2a_n= "<<tp.exchange_LET_1st__a2a_n / n_smp_int<<", max= "<<tp_max.exchange_LET_1st__a2a_n / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"  exchange_LET_1st__a2a_ep= "<<tp.exchange_LET_1st__a2a_ep / n_smp_int<<", max= "<<tp_max.exchange_LET_1st__a2a_ep / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"  exchange_LET_1st__icomm_ep= "<<tp.exchange_LET_1st__icomm_ep / n_smp_int<<", max= "<<tp_max.exchange_LET_1st__icomm_ep / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"  exchange_LET_1st__a2a_sp= "<<tp.exchange_LET_1st__a2a_sp / n_smp_int<<", max= "<<tp_max.exchange_LET_1st__a2a_sp / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;
    fout<<"  exchange_LET_1st__icomm_sp= "<<tp.exchange_LET_1st__icomm_sp / n_smp_int<<", max= "<<tp_max.exchange_LET_1st__icomm_sp / n_smp_int<<", rank= "<<rank_max[id++]<<std::endl;

    fout<<std::endl;
} 


void getTimeProfileMax(const PS::TimeProfile & tp, const PS::S32 rank, PS::TimeProfile & tp_max, PS::S32 rank_max[]){
    //F64 decompose_domain__sort_particle_1st;
    //  F64 decompose_domain__sort_particle_2nd;
    //  F64 decompose_domain__sort_particle_3rd;
    //  F64 decompose_domain__gather_particle;
    PS::S32 id = 0;
    PS::Comm::getMaxValue(tp.collect_sample_particle, rank, tp_max.collect_sample_particle, rank_max[id++]);
    PS::Comm::getMaxValue(tp.decompose_domain, rank, tp_max.decompose_domain, rank_max[id++]);

    PS::Comm::getMaxValue(tp.decompose_domain__gather_particle, rank, tp_max.decompose_domain__gather_particle, rank_max[id++]);
    PS::Comm::getMaxValue(tp.decompose_domain__sort_particle_1st, rank, tp_max.decompose_domain__sort_particle_1st, rank_max[id++]);
    PS::Comm::getMaxValue(tp.decompose_domain__sort_particle_2nd, rank, tp_max.decompose_domain__sort_particle_2nd, rank_max[id++]);
    PS::Comm::getMaxValue(tp.decompose_domain__sort_particle_3rd, rank, tp_max.decompose_domain__sort_particle_3rd, rank_max[id++]);

    PS::Comm::getMaxValue(tp.decompose_domain__setup, rank, tp_max.decompose_domain__setup, rank_max[id++]);
    PS::Comm::getMaxValue(tp.decompose_domain__determine_coord_1st, rank, tp_max.decompose_domain__determine_coord_1st, rank_max[id++]);
    PS::Comm::getMaxValue(tp.decompose_domain__migrae_particle_1st, rank, tp_max.decompose_domain__migrae_particle_1st, rank_max[id++]);
    PS::Comm::getMaxValue(tp.decompose_domain__determine_coord_2nd, rank, tp_max.decompose_domain__determine_coord_2nd, rank_max[id++]);
    PS::Comm::getMaxValue(tp.decompose_domain__determine_coord_3rd, rank, tp_max.decompose_domain__determine_coord_3rd, rank_max[id++]);
    PS::Comm::getMaxValue(tp.decompose_domain__exchange_pos_domain, rank, tp_max.decompose_domain__exchange_pos_domain, rank_max[id++]);

    PS::Comm::getMaxValue(tp.exchange_particle, rank, tp_max.exchange_particle, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_particle__find_particle, rank, tp_max.exchange_particle__find_particle, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_particle__exchange_particle, rank, tp_max.exchange_particle__exchange_particle, rank_max[id++]);

    PS::Comm::getMaxValue(tp.set_particle_local_tree, rank, tp_max.set_particle_local_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.set_root_cell, rank, tp_max.set_root_cell, rank_max[id++]);
    PS::Comm::getMaxValue(tp.make_local_tree, rank, tp_max.make_local_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.calc_moment_local_tree, rank, tp_max.calc_moment_local_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.make_LET_1st, rank, tp_max.make_LET_1st, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_LET_1st, rank, tp_max.exchange_LET_1st, rank_max[id++]);
    PS::Comm::getMaxValue(tp.make_LET_2nd, rank, tp_max.make_LET_2nd, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_LET_2nd, rank, tp_max.exchange_LET_2nd, rank_max[id++]);
    PS::Comm::getMaxValue(tp.set_particle_global_tree, rank, tp_max.set_particle_global_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.make_global_tree, rank, tp_max.make_global_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.calc_moment_global_tree, rank, tp_max.calc_moment_global_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.calc_force, rank, tp_max.calc_force, rank_max[id++]);

    PS::Comm::getMaxValue(tp.make_local_tree_tot, rank, tp_max.make_local_tree_tot, rank_max[id++]);
    PS::Comm::getMaxValue(tp.make_global_tree_tot, rank, tp_max.make_global_tree_tot, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_LET_tot, rank, tp_max.exchange_LET_tot, rank_max[id++]);

    PS::Comm::getMaxValue(tp.morton_sort_local_tree, rank, tp_max.morton_sort_local_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.link_cell_local_tree, rank, tp_max.link_cell_local_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.morton_sort_global_tree, rank, tp_max.morton_sort_global_tree, rank_max[id++]);
    PS::Comm::getMaxValue(tp.link_cell_global_tree, rank, tp_max.link_cell_global_tree, rank_max[id++]);

    PS::Comm::getMaxValue(tp.exchange_LET_1st__a2a_n, rank, tp_max.exchange_LET_1st__a2a_n, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_LET_1st__a2a_ep, rank, tp_max.exchange_LET_1st__a2a_ep, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_LET_1st__icomm_ep, rank, tp_max.exchange_LET_1st__icomm_ep, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_LET_1st__a2a_sp, rank, tp_max.exchange_LET_1st__a2a_sp, rank_max[id++]);
    PS::Comm::getMaxValue(tp.exchange_LET_1st__icomm_sp, rank, tp_max.exchange_LET_1st__icomm_sp, rank_max[id++]);

    //if(PS::Comm::getRank() == 0){ std::cout<<"tp_max.calc_force="<<tp_max.calc_force<<std::endl; }
}

template<class Tsys, class Ttree>
void DumpTimeProfile(const PS::DomainInfo & dinfo,
		     const Tsys & system,
		     const Ttree & tree,
		     std::ostream & fout,
		     const PS::S64 n_smp_int=1){
    PS::TimeProfile tp_dinfo = dinfo.getTimeProfile();
    PS::TimeProfile tp_system = system.getTimeProfile();
    PS::TimeProfile tp_tree = tree.getTimeProfile();
    PS::TimeProfile tp_dinfo_max, tp_system_max, tp_tree_max;
    const PS::S32 n_profile = 100;
    PS::S32 rank_dinfo_max[n_profile], rank_system_max[n_profile], rank_tree_max[n_profile]; 
    getTimeProfileMax(tp_dinfo, PS::Comm::getRank(), tp_dinfo_max, rank_dinfo_max);
    getTimeProfileMax(tp_system, PS::Comm::getRank(), tp_system_max, rank_system_max);
    getTimeProfileMax(tp_tree, PS::Comm::getRank(), tp_tree_max, rank_tree_max);

    PS::CountT n_int_ep_ep = tree.getNumberOfInteractionEPEPGlobal();
    PS::CountT n_int_ep_sp = tree.getNumberOfInteractionEPSPGlobal();
    /*
    PS::CountT n_op_tot = n_int_ep_ep_grav * n_op_ep_ep_grav + n_int_ep_sp_grav * n_op_ep_sp_grav;
    fout_tcal<<"time_sys= "<<time_sys<<" n_loop= "<<n_loop<<std::endl;
    fout_tcal<<"speed= "<<(PS::F64)(n_op_tot)/(wtime_tot)*1e-12<<"[Tflops]"<<std::endl;
    fout_tcal<<"efficiency= "<<(PS::F64)(n_op_tot)/(wtime_tot)/(flops_per_proc*PS::Comm::getNumberOfProc())<<std::endl;
    fout_tcal<<"wtime_tot= "<<wtime_tot<<std::endl;
    fout_tcal<<"n_op_tot= "<<n_op_tot<<std::endl;
    */
    PS::CountT n_let_send_loc = tree.getNumberOfLETEPSend1stLocal() + tree.getNumberOfLETSPSend1stLocal() + tree.getNumberOfLETEPSend2ndLocal();
    PS::CountT n_let_recv_loc = tree.getNumberOfLETEPRecv1stLocal() + tree.getNumberOfLETSPRecv1stLocal() + tree.getNumberOfLETEPRecv2ndLocal();
    PS::CountT n_let_tot_glb = tree.getNumberOfLETEPRecv1stGlobal() + tree.getNumberOfLETSPRecv1stGlobal() + tree.getNumberOfLETEPRecv2ndGlobal();
    fout<<"system.getNumberOfParticleGlobal()= "<<system.getNumberOfParticleGlobal()
        <<" system.getNumberOfParticleLocal()= "<<system.getNumberOfParticleLocal()
        <<" system.getNumberOfParticleGlobal()+n_let_tot_glb= "<<system.getNumberOfParticleGlobal() + (n_let_tot_glb / n_smp_int)<<std::endl;
    fout<<"n_let_send_loc= "<<n_let_send_loc / n_smp_int
        <<" n_let_recv_loc= "<<n_let_recv_loc / n_smp_int
        <<" n_let_ep_send_1st_local= "<<tree.getNumberOfLETEPSend1stLocal() / n_smp_int
        <<" n_let_ep_recv_1st_local= "<<tree.getNumberOfLETEPRecv1stLocal() / n_smp_int
        <<" n_let_sp_send_1st_local= "<<tree.getNumberOfLETSPSend1stLocal() / n_smp_int
        <<" n_let_sp_recv_1st_local= "<<tree.getNumberOfLETSPRecv1stLocal() / n_smp_int
        <<" n_let_ep_send_2nd_local= "<<tree.getNumberOfLETEPSend2ndLocal() / n_smp_int
        <<" n_let_ep_recv_2nd_local= "<<tree.getNumberOfLETEPRecv2ndLocal() / n_smp_int<<std::endl;
    fout<<"n_let_tot_glb= "<<n_let_tot_glb / n_smp_int
        <<" n_let_ep_send_1st_global= "<<tree.getNumberOfLETEPSend1stGlobal() / n_smp_int
        <<" n_let_ep_recv_1st_global= "<<tree.getNumberOfLETEPRecv1stGlobal() / n_smp_int
        <<" n_let_sp_send_1st_global= "<<tree.getNumberOfLETSPSend1stGlobal() / n_smp_int
        <<" n_let_sp_recv_1st_global= "<<tree.getNumberOfLETSPRecv1stGlobal() / n_smp_int
        <<" n_let_ep_send_2nd_global= "<<tree.getNumberOfLETEPSend2ndGlobal() / n_smp_int
        <<" n_let_ep_recv_2nd_global= "<<tree.getNumberOfLETEPRecv2ndGlobal() / n_smp_int<<std::endl;
    fout<<"n_ptcl_send_local= "<<system.getNumberOfParticleSendLocal() / n_smp_int
	<<" n_ptcl_recv_local= "<<system.getNumberOfParticleRecvLocal() / n_smp_int
	<<" n_ptcl_send_golbal= "<<system.getNumberOfParticleSendGlobal() / n_smp_int
	<<" n_ptcl_recv_golbal= "<<system.getNumberOfParticleRecvGlobal() / n_smp_int<<std::endl;
    //fout_tcal<<"n_cell_open_local= "<<tree_profile.getNumberOfCellOpenLocal()
    //	     <<" n_cell_open_global= "<<tree_profile.getNumberOfCellOpenGlobal()<<std::endl;


    PS::F64 n_proc_let_send_icomm_ep = (PS::F64)tree.getNumberOfProcSendLET1stICommEP();
    PS::F64 n_proc_let_recv_icomm_ep = (PS::F64)tree.getNumberOfProcRecvLET1stICommEP();	    
    PS::F64 n_proc_let_send_icomm_sp = (PS::F64)tree.getNumberOfProcSendLET1stICommSP();
    PS::F64 n_proc_let_recv_icomm_sp = (PS::F64)tree.getNumberOfProcRecvLET1stICommSP();
    PS::F64 n_proc_let_send_icomm_ep_max = 0;
    PS::F64 n_proc_let_recv_icomm_ep_max = 0;
    PS::F64 n_proc_let_send_icomm_sp_max = 0;
    PS::F64 n_proc_let_recv_icomm_sp_max = 0;
    PS::S32 rank_max_s_ep = 0; 
    PS::S32 rank_max_r_ep = 0; 
    PS::S32 rank_max_s_sp = 0; 
    PS::S32 rank_max_r_sp = 0; 
    PS::Comm::getMaxValue(n_proc_let_send_icomm_ep, PS::Comm::getRank(), n_proc_let_send_icomm_ep_max, rank_max_s_ep);
    PS::Comm::getMaxValue(n_proc_let_recv_icomm_ep, PS::Comm::getRank(), n_proc_let_recv_icomm_ep_max, rank_max_r_ep);
    PS::Comm::getMaxValue(n_proc_let_send_icomm_sp, PS::Comm::getRank(), n_proc_let_send_icomm_sp_max, rank_max_s_sp);
    PS::Comm::getMaxValue(n_proc_let_recv_icomm_sp, PS::Comm::getRank(), n_proc_let_recv_icomm_sp_max, rank_max_r_sp);

    fout<<"n_proc_let_send_icomm_ep= "<<n_proc_let_send_icomm_ep / n_smp_int
        <<" max= "<<n_proc_let_send_icomm_ep_max / n_smp_int
        <<" rank= "<<rank_max_s_ep<<std::endl;
    fout<<"n_proc_let_recv_icomm_ep= "<<n_proc_let_recv_icomm_ep / n_smp_int
        <<" max= "<<n_proc_let_recv_icomm_ep_max / n_smp_int
        <<" rank= "<<rank_max_r_ep<<std::endl;
    fout<<"n_proc_let_send_icomm_sp= "<<n_proc_let_send_icomm_sp / n_smp_int
        <<" max= "<<n_proc_let_send_icomm_sp_max / n_smp_int
        <<" rank= "<<rank_max_s_sp<<std::endl;
    fout<<"n_proc_let_recv_icomm_sp= "<<n_proc_let_recv_icomm_sp / n_smp_int
        <<" max= "<<n_proc_let_recv_icomm_sp_max / n_smp_int
        <<" rank= "<<rank_max_r_sp<<std::endl;
    //<<" n_proc_let_recv_icomm_ep= "<<n_proc_let_recv_icomm_ep / n_smp_int
    //  <<" n_proc_let_send_icomm_sp= "<<n_proc_let_send_icomm_sp / n_smp_int
    //  <<" n_proc_let_recv_icomm_sp= "<<n_proc_let_recv_icomm_sp / n_smp_int<<std::endl;

    fout<<std::endl;
    DumpTimeProfile0(tp_dinfo, tp_dinfo_max, rank_dinfo_max, fout, n_smp_int);
    DumpTimeProfile1(tp_system, tp_system_max, rank_system_max, fout, n_smp_int);
    fout<<"n_int_ep_ep= "<<n_int_ep_ep / n_smp_int<<" n_int_ep_sp= "<<n_int_ep_sp / n_smp_int<<std::endl;
    fout<<"ni_ave= "<<(PS::F64)(system.getNumberOfParticleGlobal()) / (tree.getNumberOfWalkGlobal()/n_smp_int)
	     <<" nj_ave(EP-EP)= "<<(PS::F64)(n_int_ep_ep)/n_smp_int/system.getNumberOfParticleGlobal()
	     <<" nj_ave(EP-SP)= "<<(PS::F64)(n_int_ep_sp)/n_smp_int/system.getNumberOfParticleGlobal()<<std::endl;
    DumpTimeProfile2(tp_tree, tp_tree_max, rank_tree_max, fout, n_smp_int);
    
    fout<<std::endl;
    fout<<std::endl;
    fout<<std::endl;

}


class FileHeader{
public:
    PS::S64 n_tot_glb;
    PS::S64 n_disk_glb;
    PS::S64 n_bulge_glb;
    PS::S64 n_dark_glb;
    PS::S64 n_tot_loc;
    PS::F64 time;
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%e\n", time);
        fprintf(fp, "%lld\n", n_tot_glb);
        fprintf(fp, "%lld\n", n_disk_glb);
        fprintf(fp, "%lld\n", n_bulge_glb);
        fprintf(fp, "%lld\n", n_dark_glb);
        fprintf(fp, "%lld\n", n_tot_loc);
    }
};


template<class Tpsys>
void WriteNemoAscii(const Tpsys & psys,
                    const PS::F32 time_sys,
                    const PS::S32 snp_id,
                    const char * dir_name){
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    PS::S64 n_glb = 0;
    FPGrav * fp;
    PS::AllGatherParticle(fp, n_glb, &psys[0], n_loc);
    if(PS::Comm::getRank () == 0){
        const PS::S32 STRINGSIZE = 1024;
        char sout[STRINGSIZE];
        sprintf(sout,"%s/snap%5d.dat", dir_name, snp_id);
        for(int i=0;i<STRINGSIZE;i++)if(sout[i]==' ')sout[i]='0';
        std::ofstream foutput;
        foutput.open(sout);
        foutput<<std::setprecision(15);
        foutput<<n_glb<<std::endl;
        foutput<<"3"<<std::endl;
        foutput<<time_sys<<std::endl;
        for(PS::S64 i=0; i<n_glb; i++) foutput<<fp[i].mass<<std::endl;
        for(PS::S64 i=0; i<n_glb; i++) foutput<<fp[i].pos<<std::endl;
        for(PS::S64 i=0; i<n_glb; i++) foutput<<fp[i].vel<<std::endl;
        foutput.close();
    }
    delete [] fp;
}

template<class Tpsys>
void Kick(Tpsys & system,
          const PS::F64 dt){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(int i=0; i<n; i++){
        system[i].vel  += system[i].acc * dt;
    }
}

template<class Tpsys>
void Drift(Tpsys & system,
           const PS::F64 dt){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(int i=0; i<n; i++){
        system[i].pos  += system[i].vel * dt;
    }
}

template<class Tpsys>
void CalcEnergy(const Tpsys & system,
                PS::F64 & etot,
                PS::F64 & ekin,
                PS::F64 & epot,
                const bool clear=true){
    if(clear){
        etot = ekin = epot = 0.0;
    }
    PS::F64 etot_loc = 0.0;
    PS::F64 ekin_loc = 0.0;
    PS::F64 epot_loc = 0.0;
    const PS::S32 nbody = system.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<nbody; i++){
        ekin_loc += system[i].mass * system[i].vel * system[i].vel;
        epot_loc += system[i].mass * system[i].pot;
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    etot_loc = ekin_loc + epot_loc;
    etot = PS::Comm::getSum(etot_loc);
    epot = PS::Comm::getSum(epot_loc);
    ekin = PS::Comm::getSum(ekin_loc);
}

void ReadSnapShot(PS::ParticleSystem<FPGrav> & sys, 
                  const PS::S64 & n_factor,
                  PS::F64 & time_sys,
                  PS::S64 & n_tot_glb,
                  PS::S64 & n_tot_loc,
                  PS::S64 & n_disk_glb,
                  PS::S64 & n_bulge_glb,
                  PS::S64 & n_dark_glb,
                  const char * input_dir_name, 
                  const PS::S64 n_snp_init,
                  const PS::S64 n_snp_init_limit=1024){

    std::cout<<"n_snp_init="<<n_snp_init<<std::endl;

    const PS::S32 n_proc_glb = PS::Comm::getNumberOfProc();
    const PS::S32 my_rank_glb = PS::Comm::getRank();
    assert( n_proc_glb >= n_snp_init);
    assert(n_proc_glb % n_snp_init == 0);
    MPI_Comm comm_sub;
    int color = my_rank_glb % n_snp_init;

    std::cout<<"color"<<std::endl;

    MPI_Comm_split(MPI_COMM_WORLD, color, my_rank_glb, &comm_sub);
    int my_rank_sub, n_proc_sub;
    MPI_Comm_rank(comm_sub, &my_rank_sub);
    MPI_Comm_size(comm_sub, &n_proc_sub);
    std::cout<<"n_proc_sub="<<n_proc_sub<<" my_rank_sub="<<my_rank_sub<<std::endl;


    char sinput[1024];
    std::ifstream fin;

    PS::S64 n_tot_glb_snp, n_disk_glb_snp, n_bulge_glb_snp, n_dark_glb_snp;
    PS::S64 n_tot_loc_org = 0;
    if( PS::Comm::getRank() < n_snp_init ){
        PS::S64 XXX;
        sprintf(sinput, "%s/disk_%d_%d.dat", input_dir_name, n_snp_init_limit, PS::Comm::getRank());
        std::cout<<"sinput: "<<sinput<<std::endl;
        fin.open(sinput);
        fin>>n_tot_glb_snp;
        fin>>n_disk_glb_snp;
        fin>>n_bulge_glb_snp;
        fin>>n_dark_glb_snp;
        fin>>XXX;
        n_tot_loc_org += XXX;
        fin>>XXX;
        n_tot_loc_org += XXX;
        fin>>XXX;
        n_tot_loc_org += XXX;
        fin>>time_sys;
    }
    else{
        n_tot_loc_org = 0;
    }
    PS::Comm::broadcast(&n_tot_glb_snp, 1, 0);
    PS::Comm::broadcast(&n_disk_glb_snp, 1, 0);
    PS::Comm::broadcast(&n_bulge_glb_snp, 1, 0);
    PS::Comm::broadcast(&n_dark_glb_snp, 1, 0);
    PS::Comm::broadcast(&time_sys, 1, 0);
    MPI_Bcast(&n_tot_loc_org, 1, PS::GetDataType<PS::S64>(), 0, comm_sub);
    //PS::S64 n_tot_glb_org = PS::Comm::getSum(n_tot_loc_org); // at first evaluate n_tot_glb_org befor broadcast n_tot_loc_org
    //n_tot_glb = n_tot_glb_org * n_factor;
    //std::cerr<<"time_sys="<<time_sys<<std::endl;
    //std::cout<<"n_tot_loc_org="<<n_tot_loc_org<<std::endl;


    PS::S64 * id = new PS::S64[n_tot_loc_org];
    PS::F64 * mass = new PS::F64[n_tot_loc_org];
    PS::F64vec * pos = new PS::F64vec[n_tot_loc_org];
    PS::F64vec * vel = new PS::F64vec[n_tot_loc_org];
    if( PS::Comm::getRank() < n_snp_init ){
        for(long i=0; i<n_tot_loc_org; i++){
            fin >> id[i] >> mass[i] >> pos[i] >> vel[i];
        }
        fin.close();
    }
    MPI_Bcast(id, n_tot_loc_org, PS::GetDataType<PS::S64>(), 0, comm_sub);
    MPI_Bcast(mass, n_tot_loc_org, PS::GetDataType<PS::F64>(), 0, comm_sub);
    MPI_Bcast(pos, n_tot_loc_org, PS::GetDataType<PS::F64vec>(), 0, comm_sub);
    MPI_Bcast(vel, n_tot_loc_org, PS::GetDataType<PS::F64vec>(), 0, comm_sub);

    //std::cout<<"finish of coping of snap shot. each proc has n_tot_loc_org ptcls"<<std::endl;

    n_tot_loc = (n_tot_loc_org * n_factor) / n_proc_sub;
    PS::S64 id_head = n_tot_loc * my_rank_sub;
    if(my_rank_sub < (n_tot_loc_org * n_factor) % n_proc_sub ){
        n_tot_loc++;
        id_head += my_rank_sub;
    }
    else{
        id_head += (n_tot_loc_org * n_factor) % n_proc_sub;
    }
    //std::cout<<"my_rank_sub="<<my_rank_sub<<" n_proc_sub="<<n_proc_sub<<std::endl;
    //std::cout<<"id_head="<<id_head<<" n_tot_loc="<<n_tot_loc<<std::endl;

    sys.initialize();
    sys.createParticle(n_tot_loc*3+100);
    sys.setNumberOfParticleLocal(n_tot_loc);
    n_tot_glb = PS::Comm::getSum(n_tot_loc);
    PS::F64 m_factor = (PS::F64)n_tot_glb_snp / (PS::F64)n_tot_glb;
    std::cout<<"m_factor="<<m_factor<<std::endl;
    PS::MTTS mt;
    mt.init_genrand(my_rank_glb);
    PS::S64 cnt = 0;
    for(PS::S64 i=id_head; i<id_head+n_tot_loc; i++){
        PS::S64 i_tmp = i % n_tot_loc_org;
        PS::F64 theta = 2.0 * M_PI * mt.genrand_real2();
        PS::F64 cth = cos(theta);
        PS::F64 sth = sin(theta);
        //PS::F64 cth = 1;
        //PS::F64 sth = 0;
        sys[cnt].id = id[i_tmp];
        sys[cnt].mass = mass[i_tmp] * m_factor;
        sys[cnt].pos.x = cth * pos[i_tmp].x - sth * pos[i_tmp].y;
        sys[cnt].pos.y = sth * pos[i_tmp].x + cth * pos[i_tmp].y;
        sys[cnt].vel.x = cth * vel[i_tmp].x - sth * vel[i_tmp].y;
        sys[cnt].vel.y = sth * vel[i_tmp].x + cth * vel[i_tmp].y;
        sys[cnt].pos.z = pos[i_tmp].z;
        sys[cnt].vel.z = vel[i_tmp].z;
        cnt++;
    }

    assert(cnt == n_tot_loc);

    PS::S64 n_disk_loc = 0;
    PS::S64 n_bulge_loc = 0;
    PS::S64 n_dark_loc = 0;
    for(PS::S64 i=0; i<n_tot_loc; i++){
        if( sys[i].id < n_disk_glb_snp) n_disk_loc++;
        else if( sys[i].id < n_disk_glb_snp + n_bulge_glb_snp) n_bulge_loc++;
        else n_dark_loc++;
    }

    n_disk_glb = PS::Comm::getSum(n_disk_loc);
    n_bulge_glb = PS::Comm::getSum(n_bulge_loc);
    n_dark_glb = PS::Comm::getSum(n_dark_loc);

    std::cout<<"n_disk_glb="<<n_disk_glb<<std::endl;
    std::cout<<"n_bulge_glb="<<n_bulge_glb<<std::endl;
    std::cout<<"n_dark_glb="<<n_dark_glb<<std::endl;

    assert( n_disk_glb+n_bulge_glb+n_dark_glb == n_tot_glb);

    PS::S64 * n_disk_array = new PS::S64[n_proc_glb];
    PS::S64 * n_bulge_array = new PS::S64[n_proc_glb];
    PS::S64 * n_dark_array = new PS::S64[n_proc_glb];

    PS::Comm::allGather(&n_disk_loc, 1, n_disk_array);
    PS::Comm::allGather(&n_bulge_loc, 1, n_bulge_array);
    PS::Comm::allGather(&n_dark_loc, 1, n_dark_array);

    PS::S64 * n_disk_array_disp = new PS::S64[n_proc_glb+1];
    PS::S64 * n_bulge_array_disp = new PS::S64[n_proc_glb+1];
    PS::S64 * n_dark_array_disp = new PS::S64[n_proc_glb+1];

    n_disk_array_disp[0] = n_bulge_array_disp[0] = n_dark_array_disp[0] = 0;
    for(PS::S64 i=0; i<n_proc_glb; i++){
        n_disk_array_disp[i+1] = n_disk_array_disp[i] + n_disk_array[i];
        n_bulge_array_disp[i+1] = n_bulge_array_disp[i] + n_bulge_array[i];
        n_dark_array_disp[i+1] = n_dark_array_disp[i] + n_dark_array[i];
    }

    std::cout<<"n_disk_array_disp[n_proc_glb]="<<n_disk_array_disp[n_proc_glb]<<std::endl;
    std::cout<<"n_bulge_array_disp[n_proc_glb]="<<n_bulge_array_disp[n_proc_glb]<<std::endl;
    std::cout<<"n_dark_array_disp[n_proc_glb]="<<n_dark_array_disp[n_proc_glb]<<std::endl;

    PS::S64 id_disk = n_disk_array_disp[my_rank_glb];
    PS::S64 id_bulge = n_disk_array_disp[n_proc_glb] + n_bulge_array_disp[my_rank_glb];
    PS::S64 id_dark = n_disk_array_disp[n_proc_glb] + n_bulge_array_disp[n_proc_glb] + n_dark_array_disp[my_rank_glb];

    std::cout<<"id_disk="<<id_disk<<std::endl;
    std::cout<<"id_bulge="<<id_bulge<<std::endl;
    std::cout<<"id_dark="<<id_dark<<std::endl;

    for(PS::S64 i=0; i<n_tot_loc; i++){
        if( sys[i].id < n_disk_glb_snp){
            sys[i].id = id_disk;
            id_disk++;
        }
        else if( sys[i].id < n_disk_glb_snp + n_bulge_glb_snp){
            sys[i].id = id_bulge;
            id_bulge++;
        }
        else{
            sys[i].id = id_dark;
            id_dark++;
        }
    }

    std::cout<<"id_disk2="<<id_disk<<std::endl;
    std::cout<<"id_bulge2="<<id_bulge<<std::endl;
    std::cout<<"id_dark2="<<id_dark<<std::endl;

    PS::S64 cnt_disk = 0;
    PS::S64 cnt_bulge = 0;
    PS::S64 cnt_dark = 0;
    for(PS::S64 i=0; i<n_tot_loc; i++){
        if( sys[i].id < n_disk_glb) cnt_disk++;
        else if( sys[i].id < n_disk_glb + n_bulge_glb) cnt_bulge++;
        else cnt_dark++;
    }

    std::cout<<"cnt_disk="<<cnt_disk<<std::endl;
    std::cout<<"cnt_bulge="<<cnt_bulge<<std::endl;
    std::cout<<"cnt_dark="<<cnt_dark<<std::endl;

    assert(cnt_disk == n_disk_loc);
    assert(cnt_bulge == n_bulge_loc);
    assert(cnt_dark == n_dark_loc);

    PS::F64 mass_CM_loc = 0.0;
    PS::F64vec pos_CM_loc = 0.0;
    PS::F64vec vel_CM_loc = 0.0;

    for(PS::S64 i=0; i<n_tot_loc; i++){
        mass_CM_loc += sys[i].mass;
        pos_CM_loc += sys[i].mass * sys[i].pos;
        vel_CM_loc += sys[i].mass * sys[i].vel;
    }

    std::cerr<<"mass_CM_loc="<<mass_CM_loc<<std::endl;
    std::cerr<<"pos_CM_loc="<<pos_CM_loc<<std::endl;
    std::cerr<<"vel_CM_loc="<<vel_CM_loc<<std::endl;

    PS::F64 mass_CM_glb = PS::Comm::getSum(mass_CM_loc);
    PS::F64vec pos_CM_glb = PS::Comm::getSum(pos_CM_loc);
    PS::F64vec vel_CM_glb = PS::Comm::getSum(vel_CM_loc);

    pos_CM_glb /= mass_CM_glb;
    vel_CM_glb /= mass_CM_glb;

    std::cerr<<"mass_CM_glb="<<mass_CM_glb<<std::endl;
    std::cerr<<"pos_CM_glb="<<pos_CM_glb<<std::endl;
    std::cerr<<"vel_CM_glb="<<vel_CM_glb<<std::endl;

    for(PS::S64 i=0; i<n_tot_loc; i++){
        sys[i].pos -= pos_CM_glb;
        sys[i].vel -= vel_CM_glb;
    }

#if 0
    mass_CM_loc = 0.0;
    pos_CM_loc = 0.0;
    vel_CM_loc = 0.0;
    for(PS::S64 i=0; i<n_tot_loc; i++){
        mass_CM_loc += sys[i].mass;
        pos_CM_loc += sys[i].mass * sys[i].pos;
        vel_CM_loc += sys[i].mass * sys[i].vel;
    }
    mass_CM_glb = PS::Comm::getSum(mass_CM_loc);
    pos_CM_glb = PS::Comm::getSum(pos_CM_loc);
    vel_CM_glb = PS::Comm::getSum(vel_CM_loc);
    pos_CM_glb /= mass_CM_glb;
    vel_CM_glb /= mass_CM_glb;
    std::cerr<<"pos_CM_glb="<<pos_CM_glb<<std::endl;
    std::cerr<<"vel_CM_glb="<<vel_CM_glb<<std::endl;
#endif

    std::cout<<"before delete"<<std::endl;

    delete [] n_disk_array;
    delete [] n_bulge_array;
    delete [] n_dark_array;
    delete [] n_disk_array_disp;
    delete [] n_bulge_array_disp;
    delete [] n_dark_array_disp;

    delete [] id;
    delete [] mass;
    delete [] pos;
    delete [] vel;

    std::cout<<"after delete"<<std::endl;

}

int main(int argc, char *argv[]){
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

    PS::Initialize(argc, argv);
    PS::F64 wtime_all = PS::GetWtime();


    const PS::F64 flops_per_proc = 1.28e11; // for K computer (8thread)
    const PS::F64 n_op_ep_ep_grav = 29.0;
    const PS::F64 n_op_ep_sp_grav = 64.0;

    PS::F64 Tbegin = PS::GetWtime();
    PS::F32 theta = 0.5;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::S32 n_smp_ave = 30;
    PS::F32 time_end = 10.0;
    PS::S64 factor = 1;
    PS::S64 n_snp_init = 1;
    //PS::S64 n_ave = 10000000;
    PS::S64 reduce_factor = 1;
    char input_dir_name[1024];
    char output_dir_name[1024];
    int c;
    while((c=getopt(argc,argv,"i:o:t:T:n:N:f:s:l:x:h")) != -1){
        switch(c){
        case 'i':
            sprintf(input_dir_name, optarg);
            break;
        case 'o':
            sprintf(output_dir_name, optarg);
            break;
        case 't':
            theta = atof(optarg);
            //std::cerr<<"theta="<<theta<<std::endl;
            break;
        case 'T':
            time_end = atof(optarg);
            //std::cerr<<"time_end="<<time_end<<std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            //std::cerr<<"n_group_limit="<<n_group_limit<<std::endl;
            break;
        case 'N':
            n_snp_init = atoi(optarg);
            //std::cerr<<"n_group_limit="<<n_group_limit<<std::endl;
            break;
        case 'f':
            factor = atol(optarg);
            //std::cerr<<"factor="<<factor<<std::endl;
            break;
        case 's':
            n_smp_ave = atoi(optarg);
            //std::cerr<<"n_smp_ave="<<n_smp_ave<<std::endl;
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            //std::cerr<<"n_leaf_limit="<<n_leaf_limit<<std::endl;
            break;
        case 'x':
            //n_ave = atol(optarg);
            reduce_factor = atol(optarg);
            break;
        case 'h':
            std::cerr<<"i: input file name (nemo ascii)"<<std::endl;
            std::cerr<<"o: dir name of output"<<std::endl;
            std::cerr<<"t: theta (dafult: 0.5)"<<std::endl;
            std::cerr<<"T: time_end (dafult: 10.0)"<<std::endl;
            std::cerr<<"n: n_group_limit (dafult: 64.0)"<<std::endl;
            std::cerr<<"N: n_snp_init (# of initial snap shot used. dafult: 1)"<<std::endl;
            std::cerr<<"f: n_factor (must >= 1)"<<std::endl;
            std::cerr<<"s: n_smp_ave (dafult: 30)"<<std::endl;
            std::cerr<<"l: n_leaf_limit (dafult: 8)"<<std::endl;
            //std::cerr<<"x: n_ave (dafult: 10000000)"<<std::endl;
            std::cerr<<"x: reduce_factor (integer) (dafult: 1)"<<std::endl;
            return 0;
        }
    }
    std::ofstream fout_eng;
    std::ofstream fout_tcal;
    std::ofstream fout_log;



    if(1){
        char sout_de[1024];
        char sout_tcal[1024];
        char sout_log[1024];
        sprintf(sout_de, "%s/t-de_%05d_%05d.dat", output_dir_name, PS::Comm::getNumberOfProc(), PS::Comm::getRank() );
        sprintf(sout_tcal, "%s/t-tcal_%05d_%05d.dat", output_dir_name, PS::Comm::getNumberOfProc(), PS::Comm::getRank());
        sprintf(sout_log, "%s/log_%05d_%05d.dat", output_dir_name, PS::Comm::getNumberOfProc(), PS::Comm::getRank());
        if(PS::Comm::getRank() == 0){
            std::cerr<<sout_de<<std::endl;
            std::cerr<<sout_tcal<<std::endl;
            std::cerr<<sout_log<<std::endl;
        }
        fout_eng.open(sout_de);
        fout_tcal.open(sout_tcal);
        fout_log.open(sout_log);
    }
    fout_log<<"Comm::getNumberOfProc()="<<PS::Comm::getNumberOfProc()<<std::endl;
    fout_log<<"Comm::getNumberOfThread()="<<PS::Comm::getNumberOfThread()<<std::endl;
    fout_log<<"output_dir_name="<<output_dir_name<<std::endl;
    fout_log<<"input_dir_name="<<input_dir_name<<std::endl;
    fout_log<<"n_snp_init="<<n_snp_init<<std::endl;
    fout_log<<"theta="<<theta<<std::endl;
    fout_log<<"time_end="<<time_end<<std::endl;
    fout_log<<"n_group_limit="<<n_group_limit<<std::endl;
    fout_log<<"n_smp_ave="<<n_smp_ave<<std::endl;
    fout_log<<"n_leaf_limit="<<n_leaf_limit<<std::endl;

    PS::ParticleSystem<FPGrav> system_grav;

    PS::F64 time_sys = 0.0;
    PS::S64 n_tot_glb = 0;
    PS::S64 n_tot_loc = 0;

    PS::S64 n_disk_glb = 0;
    PS::S64 n_bulge_glb = 0;
    PS::S64 n_dark_glb = 0;
    PS::S64 n_snp_init_limit = 1024;
    ReadSnapShot(system_grav, factor, time_sys, n_tot_glb, n_tot_loc, n_disk_glb, n_bulge_glb, n_dark_glb, input_dir_name, n_snp_init, n_snp_init_limit);
    system_grav.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave);

    fout_log<<"finish set particle: time="<<PS::GetWtime() - Tbegin<<std::endl;

    fout_log<<"n_tot_glb="<<n_tot_glb<<std::endl;
    fout_log<<"n_tot_loc="<<n_tot_loc<<std::endl;
    fout_log<<"n_disk_glb="<<n_disk_glb<<std::endl;
    fout_log<<"n_bulge_glb="<<n_bulge_glb<<std::endl;
    fout_log<<"n_dark_glb="<<n_dark_glb<<std::endl;

/*
    if(n_ave < n_tot_loc){
	PS::S64 n_int = n_tot_loc / n_ave;
	PS::S64 n_loc_new = 0;
	for(PS::S64 i=0; i<n_tot_loc; i += n_int){
	    system_grav[n_loc_new++] = system_grav[i];
	}
	if(PS::Comm::getRank() == 0){
	    std::cout<<"n_loc_new="<<n_loc_new<<" n_tot_loc="<<n_tot_loc<<std::endl;
	}
	n_tot_loc = n_loc_new;
    }
    system_grav.setNumberOfParticleLocal(n_tot_loc);
*/
	PS::S64 n_loc_new = 0;
    if(reduce_factor > 1){
        for(PS::S64 i=0; i<n_tot_loc; i += reduce_factor ){
            system_grav[n_loc_new++] = system_grav[i];
        }
        n_tot_loc = n_loc_new;
    }
    system_grav.setNumberOfParticleLocal(n_tot_loc);

    const PS::F32 coef_ema = 0.2;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    fout_log<<"finish dinof.initialize="<<PS::GetWtime() - Tbegin<<std::endl;

    PS::MT::initialize();
    dinfo.collectSampleParticle(system_grav, true);
    fout_log<<"finish dinof.collect="<<PS::GetWtime() - Tbegin<<std::endl;
    dinfo.decomposeDomain();
    fout_log<<"dinfo.getPosDomain(PS::Comm::getRank())="<<dinfo.getPosDomain(PS::Comm::getRank())<<std::endl;
    fout_log<<"finish dinof.decompose="<<PS::GetWtime() - Tbegin<<std::endl;
    PS::F64ort * pos_domain_org = new PS::F64ort[PS::Comm::getNumberOfProc()];
    for(int i=0; i<PS::Comm::getNumberOfProc(); i++){
        pos_domain_org[i] = dinfo.getPosDomain(i);
    }


    system_grav.exchangeParticle(dinfo);
    //system_grav.exchangeParticleSafely(dinfo);
    //system_grav.exchangeParticleSafeMode(dinfo);

    fout_log<<"finish dinof.exchangeparticle="<<PS::GetWtime() - Tbegin<<std::endl;

    PS::S64 n_grav_loc = system_grav.getNumberOfParticleLocal();
    const PS::S64 n_grav_glb = system_grav.getNumberOfParticleGlobal();

    fout_log<<"check 1"<<std::endl;

    PS::TreeForForceLong<ForceGrav, EPIGrav, EPJGrav>::Quadrupole tree_grav;
    //tree_grav.initialize(n_grav_glb, 0.7, n_leaf_limit, n_group_limit);
    tree_grav.initialize(n_grav_glb, 0.7, 8, 256);

    fout_log<<"system_grav.getNumberOfParticleLocal()="<<system_grav.getNumberOfParticleLocal()<<std::endl;
    fout_log<<"system_grav.getNumberOfParticleGlobal()="<<system_grav.getNumberOfParticleGlobal()<<std::endl;

    for(PS::S64 i=0; i<system_grav.getNumberOfParticleLocal(); i++){
	if( (PS::F64)n_grav_glb > 2e4 ){
	    system_grav[i].r_search = 1.0 * pow((PS::F64)n_grav_glb/1e4, -1.0/3.0);
	}
	else{
	    system_grav[i].r_search = 1.0;
	}
    }

    fout_log<<"check 3"<<std::endl;

    PS::TreeForForceShort<ForceRSearch, EPIRSearch, EPJRSearch>::Gather tree_search;
    tree_search.initialize(n_grav_glb, theta, n_leaf_limit, n_group_limit);

    fout_log<<"check 4"<<std::endl;

    PS::S64 cnt = 0;
    for(bool repeat = true; repeat==true;){
        std::cout<<"cnt="<<cnt<<std::endl;
        cnt++;
        bool repeat_loc = false;
        repeat = false;
        tree_search.calcForceAll(CalcRSearch(), system_grav, dinfo, false);
        for(int i=0; i<system_grav.getNumberOfParticleLocal() ; i++){
            if(tree_search.getForce(i).r_search != 0.0){
                repeat_loc = true;
                system_grav[i].r_search = tree_search.getForce(i).r_search * 1.5;
            }
            else{
                system_grav[i].r_search = tree_search.getForce(i).r_search;
            }
        }
        repeat = PS::Comm::synchronizeConditionalBranchOR(repeat_loc);
    }
    fout_log<<"check 5"<<std::endl;
    for(int i=0; i<system_grav.getNumberOfParticleLocal(); i++){
        system_grav[i].r_search = tree_search.getForce(i).r_cut;
        //if(PS::Comm::getRank() == 0) system_grav[i].dump(std::cout);
    }
    tree_search.freeMem();

    fout_log<<"check 6"<<std::endl;

    PS::F64vec * pos_org = new PS::F64vec[system_grav.getNumberOfParticleLocal()];
    PS::F64vec * vel_org = new PS::F64vec[system_grav.getNumberOfParticleLocal()];
    PS::F64 * r_search_org = new PS::F64[system_grav.getNumberOfParticleLocal()];
    PS::S64 n_loc_org = system_grav.getNumberOfParticleLocal();
    for(int i=0; i<system_grav.getNumberOfParticleLocal(); i++){
	pos_org[i] = system_grav[i].pos;
	vel_org[i] = system_grav[i].vel;
	r_search_org[i] = system_grav[i].r_search;
    }
    fout_log<<"check 7"<<std::endl;

    const PS::F32 dt = 1.0/128.0; // 7.66e4 yr
    PS::S64 n_loop_head = 20;
    PS::S64 n_smp_int = 10;
    PS::S64 n_loop_max = n_loop_head + n_smp_int;
    PS::F64 wtime_offset;
    PS::F64 wtime_tot;
    PS::S64 n_loop = 0;


#ifndef SIMPLE_VERSION
    ///////////
    // ShortScatter
    system_grav.setNumberOfParticleLocal(n_loc_org);
    for(int i=0; i<system_grav.getNumberOfParticleLocal(); i++){
        system_grav[i].pos = pos_org[i];
        system_grav[i].vel = vel_org[i];
        system_grav[i].r_search = r_search_org[i];
    }
    PS::MT::initialize();
    //dinfo.collectSampleParticle(system_grav, true);
    //dinfo.decomposeDomain();
    for(int i=0; i<PS::Comm::getNumberOfProc(); i++){
        dinfo.setPosDomain(i, pos_domain_org[i]);
    }
    fout_tcal<<"dinfo.getPosDomain(PS::Comm::getRank())= "
             <<dinfo.getPosDomain(PS::Comm::getRank())<<std::endl;
    system_grav.exchangeParticle(dinfo);
    PS::TreeForForceShort<ForceGrav, EPIGrav, EPJGrav>::Scatter tree_scatter;
    tree_scatter.initialize(n_grav_glb, theta, n_leaf_limit, n_group_limit);
    tree_grav.calcForceAllAndWriteBack(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo);
    n_loop = 0;

    while(n_loop < n_loop_max){
        Kick(system_grav, dt*0.5);
        PS::Comm::barrier();
        if(n_loop == n_loop_head){
            wtime_tot = 0.0;
            dinfo.clearTimeProfile();
            system_grav.clearCounterAll();
            tree_scatter.clearCounterAll();
            wtime_offset = PS::GetWtime();
        }
        Drift(system_grav, dt);
        if( n_loop < 4){
            dinfo.collectSampleParticle(system_grav);
            dinfo.decomposeDomain();
        }
        else{
            dinfo.collectSampleParticle(system_grav, true);
            dinfo.decomposeDomainMultiStep();
        }
        system_grav.exchangeParticle(dinfo);
        tree_scatter.calcForceAll(CalcForceDummy<EPIGrav, EPJGrav, ForceGrav>(), system_grav, dinfo, true);
        tree_grav.calcForceAllAndWriteBack(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo, true);
        Kick(system_grav, dt*0.5);
        n_loop++;
    }
    fout_tcal<<"**********************"<<std::endl;
    fout_tcal<<"MODE= SHORT_SCATTER"<<std::endl;
    DumpTimeProfile(dinfo, system_grav, tree_scatter, fout_tcal, n_smp_int);
    tree_scatter.freeMem();
#endif //SIMPLE_VERSION

#ifndef SIMPLE_VERSION
    ///////////
    // ShortGather
    system_grav.setNumberOfParticleLocal(n_loc_org);
    for(int i=0; i<system_grav.getNumberOfParticleLocal(); i++){
        system_grav[i].pos = pos_org[i];
        system_grav[i].vel = vel_org[i];
        system_grav[i].r_search = r_search_org[i];
    }
    PS::MT::initialize();
    //dinfo.collectSampleParticle(system_grav, true);
    //dinfo.decomposeDomain();
    for(int i=0; i<PS::Comm::getNumberOfProc(); i++){
        dinfo.setPosDomain(i, pos_domain_org[i]);
    }
    fout_tcal<<"dinfo.getPosDomain(PS::Comm::getRank())= "
             <<dinfo.getPosDomain(PS::Comm::getRank())<<std::endl;
    system_grav.exchangeParticle(dinfo);
    PS::TreeForForceShort<ForceGrav, EPIGrav, EPJGrav>::Gather tree_gather;
    tree_gather.initialize(n_grav_glb, theta, n_leaf_limit, n_group_limit);
    tree_grav.calcForceAllAndWriteBack(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo);
    n_loop = 0;
    while(n_loop < n_loop_max){
        Kick(system_grav, dt*0.5);
        PS::Comm::barrier();
	if(n_loop == n_loop_head){
	    wtime_tot = 0.0;
	    dinfo.clearTimeProfile();
	    system_grav.clearCounterAll();
	    tree_gather.clearCounterAll();
	    wtime_offset = PS::GetWtime();
	}
        Drift(system_grav, dt);
        if( n_loop < 4){
            dinfo.collectSampleParticle(system_grav);
            dinfo.decomposeDomain();
        }
	else{
	    dinfo.collectSampleParticle(system_grav, true);
            dinfo.decomposeDomainMultiStep();
	}
        system_grav.exchangeParticle(dinfo);
        tree_gather.calcForceAll(CalcForceDummy<EPIGrav, EPJGrav, ForceGrav>(), system_grav, dinfo, true);
        tree_grav.calcForceAllAndWriteBack(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo, true);
        Kick(system_grav, dt*0.5);
        n_loop++;
    }
    fout_tcal<<"**********************"<<std::endl;
    fout_tcal<<"MODE= SHORT_GATHER"<<std::endl;
    DumpTimeProfile(dinfo, system_grav, tree_gather, fout_tcal, n_smp_int);
    tree_gather.freeMem();
#endif //SIMPLE_VERSION

#ifndef SIMPLE_VERSION
    ///////////
    // ShortSymmetry
    system_grav.setNumberOfParticleLocal(n_loc_org);
    for(int i=0; i<system_grav.getNumberOfParticleLocal(); i++){
        system_grav[i].pos = pos_org[i];
        system_grav[i].vel = vel_org[i];
        system_grav[i].r_search = r_search_org[i];
    }
    PS::MT::initialize();
    //dinfo.collectSampleParticle(system_grav, true);
    //dinfo.decomposeDomain();
    for(int i=0; i<PS::Comm::getNumberOfProc(); i++){
        dinfo.setPosDomain(i, pos_domain_org[i]);
    }
    fout_tcal<<"dinfo.getPosDomain(PS::Comm::getRank())= "
             <<dinfo.getPosDomain(PS::Comm::getRank())<<std::endl;
    system_grav.exchangeParticle(dinfo);
    PS::TreeForForceShort<ForceGrav, EPIGrav, EPJGrav>::Symmetry tree_symmetry;
    tree_symmetry.initialize(n_grav_glb, theta, n_leaf_limit, n_group_limit);
    tree_grav.calcForceAllAndWriteBack(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo);
    n_loop = 0;
    while(n_loop < n_loop_max){
        Kick(system_grav, dt*0.5);
        PS::Comm::barrier();
        if(n_loop == n_loop_head){
            wtime_tot = 0.0;
            dinfo.clearTimeProfile();
            system_grav.clearCounterAll();
            tree_symmetry.clearCounterAll();
            wtime_offset = PS::GetWtime();
        }
        Drift(system_grav, dt);
        if( n_loop < 4){
            dinfo.collectSampleParticle(system_grav);
            dinfo.decomposeDomain();
        }
        else{
            dinfo.collectSampleParticle(system_grav, true);
            dinfo.decomposeDomainMultiStep();
        }
        system_grav.exchangeParticle(dinfo);
        tree_symmetry.calcForceAll(CalcForceDummy<EPIGrav, EPJGrav, ForceGrav>(), system_grav, dinfo, true);
        tree_grav.calcForceAllAndWriteBack(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo, true);
        Kick(system_grav, dt*0.5);
	n_loop++;
    }
    fout_tcal<<"**********************"<<std::endl;
    fout_tcal<<"MODE= SHORT_SYMMETRY"<<std::endl;
    DumpTimeProfile(dinfo, system_grav, tree_symmetry, fout_tcal, n_smp_int);
    tree_symmetry.freeMem();
#endif //SIMPLE_VERSION


    ///////////
    // Long (mono) no cutoff
    system_grav.setNumberOfParticleLocal(n_loc_org);
    for(int i=0; i<system_grav.getNumberOfParticleLocal(); i++){
        system_grav[i].pos = pos_org[i];
        system_grav[i].vel = vel_org[i];
        system_grav[i].r_search = r_search_org[i];
    }
    PS::MT::initialize();
    //dinfo.collectSampleParticle(system_grav, true);
    //dinfo.decomposeDomain();
    for(int i=0; i<PS::Comm::getNumberOfProc(); i++){
        dinfo.setPosDomain(i, pos_domain_org[i]);
    }
    fout_tcal<<"dinfo.getPosDomain(PS::Comm::getRank())= "
             <<dinfo.getPosDomain(PS::Comm::getRank())<<std::endl;
    system_grav.exchangeParticle(dinfo);
    PS::TreeForForceLong<ForceGrav, EPIGrav, EPJGrav>::Monopole tree_long_mono;
    tree_long_mono.initialize(n_grav_glb, theta, n_leaf_limit, n_group_limit);
    tree_grav.calcForceAllAndWriteBack(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo);
    n_loop = 0;
    while(n_loop < n_loop_max){
        Kick(system_grav, dt*0.5);
        PS::Comm::barrier();
        if(n_loop == n_loop_head){
            wtime_tot = 0.0;
            dinfo.clearTimeProfile();
            system_grav.clearCounterAll();
            tree_long_mono.clearCounterAll();
            wtime_offset = PS::GetWtime();
        }
        Drift(system_grav, dt);
        if( n_loop < 4){
            dinfo.collectSampleParticle(system_grav);
            dinfo.decomposeDomain();
        }
        else{
            dinfo.collectSampleParticle(system_grav, true);
            dinfo.decomposeDomainMultiStep();
        }
        system_grav.exchangeParticle(dinfo);
        tree_long_mono.calcForceAll(CalcForceDummy<EPIGrav, EPJGrav, ForceGrav>(), CalcForceDummy<EPIGrav, PS::SPJMonopole, ForceGrav>(), system_grav, dinfo, true);
        tree_grav.calcForceAllAndWriteBack(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo, true);
        Kick(system_grav, dt*0.5);
        n_loop++;
    }
    fout_tcal<<"**********************"<<std::endl;
    fout_tcal<<"MODE= LONG_NO_CUTOFF_MONO"<<std::endl;
    DumpTimeProfile(dinfo, system_grav, tree_long_mono, fout_tcal, n_smp_int);
    tree_long_mono.freeMem();


    ///////////
    // Long (quad) no cutoff
    system_grav.setNumberOfParticleLocal(n_loc_org);
    for(int i=0; i<system_grav.getNumberOfParticleLocal(); i++){
        system_grav[i].pos = pos_org[i];
        system_grav[i].vel = vel_org[i];
        system_grav[i].r_search = r_search_org[i];
    }
    PS::MT::initialize();
    //dinfo.collectSampleParticle(system_grav, true);
    //dinfo.decomposeDomain();
    for(int i=0; i<PS::Comm::getNumberOfProc(); i++){
        dinfo.setPosDomain(i, pos_domain_org[i]);
    }
    fout_tcal<<"dinfo.getPosDomain(PS::Comm::getRank())= "
             <<dinfo.getPosDomain(PS::Comm::getRank())<<std::endl;
    system_grav.exchangeParticle(dinfo);
    PS::TreeForForceLong<ForceGrav, EPIGrav, EPJGrav>::Quadrupole tree_long_quad;
    tree_long_quad.initialize(n_grav_glb, theta, n_leaf_limit, n_group_limit);
    tree_grav.calcForceAllAndWriteBack(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo);
    n_loop = 0;
    while(n_loop < n_loop_max){
        Kick(system_grav, dt*0.5);
        PS::Comm::barrier();
        if(n_loop == n_loop_head){
            wtime_tot = 0.0;
            dinfo.clearTimeProfile();
            system_grav.clearCounterAll();
            tree_long_quad.clearCounterAll();
            wtime_offset = PS::GetWtime();
        }
        Drift(system_grav, dt);
        if( n_loop < 4){
            dinfo.collectSampleParticle(system_grav);
            dinfo.decomposeDomain();
        }
        else{
            dinfo.collectSampleParticle(system_grav, true);
            dinfo.decomposeDomainMultiStep();
        }
        system_grav.exchangeParticle(dinfo);
        tree_long_quad.calcForceAll(CalcForceDummy<EPIGrav, EPJGrav, ForceGrav>(), CalcForceDummy<EPIGrav, PS::SPJQuadrupole, ForceGrav>(), system_grav, dinfo, true);
        tree_grav.calcForceAllAndWriteBack(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo, true);
        Kick(system_grav, dt*0.5);
        n_loop++;
    }
    fout_tcal<<"**********************"<<std::endl;
    fout_tcal<<"MODE= LONG_NO_CUTOFF_QUAD"<<std::endl;
    DumpTimeProfile(dinfo, system_grav, tree_long_quad, fout_tcal, n_smp_int);
    tree_long_quad.freeMem();


    ///////////
    // Long Cutoff
    system_grav.setNumberOfParticleLocal(n_loc_org);
    for(int i=0; i<system_grav.getNumberOfParticleLocal(); i++){
        system_grav[i].pos = pos_org[i];
        system_grav[i].vel = vel_org[i];
        system_grav[i].r_search = r_search_org[i];
    }
    PS::F64 r_min_loc = 9999999.9;
    PS::F64 r_max_loc = -9999999.9;
    for(int i=0; i<system_grav.getNumberOfParticleLocal(); i++){
        for(int k=0; k<3; k++){
            if( r_min_loc > system_grav[i].pos[k]){
                r_min_loc = system_grav[i].pos[k];
            }
            if( r_max_loc < system_grav[i].pos[k]){
                r_max_loc = system_grav[i].pos[k];
            }
        }
    }
    const PS::F64 r_max_glb = PS::Comm::getMaxValue(r_max_loc);
    const PS::F64 r_min_glb = PS::Comm::getMinValue(r_min_loc);
    const PS::F64 r_cut = 3.0 * ((r_max_glb - r_min_glb) / pow(system_grav.getNumberOfParticleGlobal(),  1.0/3.0) / 2.0);
    for(int i=0; i<system_grav.getNumberOfParticleLocal(); i++){
        system_grav[i].r_search = r_cut;
    }
    PS::MT::initialize();
    //dinfo.collectSampleParticle(system_grav, true);
    //dinfo.decomposeDomain();
    for(int i=0; i<PS::Comm::getNumberOfProc(); i++){
        dinfo.setPosDomain(i, pos_domain_org[i]);
    }
    fout_tcal<<"dinfo.getPosDomain(PS::Comm::getRank())= "
             <<dinfo.getPosDomain(PS::Comm::getRank())<<std::endl;
    system_grav.exchangeParticle(dinfo);
    PS::TreeForForceLong<ForceGrav, EPIGrav, EPJGrav>::MonopoleWithCutoff tree_long_cutoff;
    tree_long_cutoff.initialize(n_grav_glb, theta, n_leaf_limit, n_group_limit);
    tree_grav.calcForceAllAndWriteBack(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo);
    n_loop = 0;
    while(n_loop < n_loop_max){
        Kick(system_grav, dt*0.5);
        PS::Comm::barrier();
        if(n_loop == n_loop_head){
            wtime_tot = 0.0;
            dinfo.clearTimeProfile();
            system_grav.clearCounterAll();
            tree_long_cutoff.clearCounterAll();
            wtime_offset = PS::GetWtime();
        }
        Drift(system_grav, dt);
        if( n_loop < 4){
            dinfo.collectSampleParticle(system_grav);
            dinfo.decomposeDomain();
        }
        else{
            dinfo.collectSampleParticle(system_grav, true);
            dinfo.decomposeDomainMultiStep();
        }
        system_grav.exchangeParticle(dinfo);
        tree_long_cutoff.calcForceAll(CalcForceDummy<EPIGrav, EPJGrav, ForceGrav>(), CalcForceDummy<EPIGrav, PS::SPJMonopoleCutoff, ForceGrav>(), system_grav, dinfo, true);
        tree_grav.calcForceAllAndWriteBack(CalcForceEpEp(), CalcForceSpEp(), system_grav, dinfo, true);
        Kick(system_grav, dt*0.5);
        n_loop++;
    }
    fout_tcal<<"**********************"<<std::endl;
    fout_tcal<<"MODE= LONG_CUTOFF"<<std::endl;
    DumpTimeProfile(dinfo, system_grav, tree_long_cutoff, fout_tcal, n_smp_int);
    tree_long_cutoff.freeMem();


    wtime_all = PS::GetWtime() - wtime_all;
    fout_tcal<<"wtime_all= "<<wtime_all<<std::endl;

    PS::Finalize();
    return 0;
}
