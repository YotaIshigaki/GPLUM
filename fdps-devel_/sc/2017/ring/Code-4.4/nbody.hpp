#pragma once

#include<sstream>

double wtime_calc_force = 0.0;
double wtime_calc_force_d_d = 0.0;
double wtime_calc_force_d_s = 0.0;
double wtime_calc_force_p_ds = 0.0;

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

void DumpProfileParticleSystem(const PS::TimeProfile & prof,
                               std::ofstream & fout){
    fout<<"exchange_particle= "<<prof.exchange_particle<<std::endl;
    fout<<"    find_particle= "<<prof.exchange_particle__find_particle<<std::endl;
    fout<<"    communication= "<<prof.exchange_particle__exchange_particle<<std::endl;
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
    fout<<"        wtime_alltoallv_ep= "<<PS::wtime_alltoallv_ep<<std::endl;
    fout<<"        wtime_alltoallv_sp= "<<PS::wtime_alltoallv_sp<<std::endl;
    fout<<"            wtime_alltoallv_sp_allgather= "<<PS::wtime_alltoallv_sp_allgather<<std::endl;
    fout<<"            wtime_alltoallv_sp_isendrecv= "<<PS::wtime_alltoallv_sp_isendrecv<<std::endl;
    fout<<"                wtime_alltoallv_sp_isendrecv_1st= "<<PS::wtime_alltoallv_sp_isendrecv_1st<<std::endl;
    fout<<"                wtime_alltoallv_sp_isendrecv_2nd= "<<PS::wtime_alltoallv_sp_isendrecv_2nd<<std::endl;
    fout<<"    set_particle_global_tree= "<<prof.set_particle_global_tree<<std::endl;
    fout<<"    morton_sort_global_tree= "<<prof.morton_sort_global_tree<<std::endl;
    fout<<"    link_cell_global_tree= "<<prof.link_cell_global_tree<<std::endl;
    fout<<"    calc_moment_global_tree= "<<prof.calc_moment_global_tree<<std::endl;
    fout<<"    add_moment_as_sp_global_tree= "<<prof.add_moment_as_sp_global_tree<<std::endl;
    fout<<"    make_ipgroup= "<<prof.calc_force__make_ipgroup<<std::endl;
    fout<<"    make_all_interaction_list_id= "<<prof.make_all_interaction_list_id<<std::endl;
    fout<<"    calc_force= "<<prof.calc_force<<std::endl;
    fout<<"        wtime_dispatch= "<<PS::wtime_dispatch<<std::endl;
    fout<<"            wtime_kick= "<<wtime_kick<<std::endl;
    fout<<"            wtime_drift= "<<wtime_drift<<std::endl;
    fout<<"            wtime_calc_energy= "<<wtime_calc_energy<<std::endl;
    fout<<"        wtime_prefixsum_recorder= "<<PS::wtime_prefixsum_recorder<<std::endl;
    fout<<"        wtime_retrieve= "<<PS::wtime_retrieve<<std::endl;
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
        + prof.set_particle_global_tree
        + prof.morton_sort_global_tree
        + prof.link_cell_global_tree
        + prof.calc_moment_global_tree
        + prof.add_moment_as_sp_global_tree
        + prof.calc_force__make_ipgroup
        + prof.make_all_interaction_list_id
        + prof.calc_force
        + prof.calc_force__copy_original_order
        <<std::endl;
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
    fout<<std::endl;

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
    const PS::F64 time_offset = PS::GetWtime();
    const PS::F64 time_offset_d_d = PS::GetWtime();
    PS::S32 n_loc = system.getNumberOfParticleLocal();
#if defined(SUNWAY) && defined(SUNWAY_FORCE_KERNEL)
        // In this case, we use Long's force kernel.
        cpe_pars.energy_flag = 0;
        cpe_pars.first_flag = 1; // Half-kick
        cpe_pars.last_flag = 0; // Full-drift
#endif
    tree.calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe
        (func_dis, func_ret, 0,
         system, dinfo, n_walk_limit,
         false,  reuse_flag);
    wtime_calc_force_d_d += PS::GetWtime() - time_offset_d_d;
    wtime_calc_force += PS::GetWtime() - time_offset;
    //PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "OK(CalcForce)" << std::endl;
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
                        const PS::S32 n_loop_merge){
    bool reuse_flag = false;
    DISPATCH_MODE = 0;

    for(PS::S32 i=0; i<n_loop_merge; i++){
        if(i==n_loop_merge-1)  DISPATCH_MODE = 1;
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
        reuse_flag = true;
        DISPATCH_MODE = 0;
        //////////////
        // TIMING DUMP
        PS::TimeProfile prof_tree    = TREE_POINTER->getTimeProfile();
        PS::TimeProfile prof_system  = system.getTimeProfile();
        PS::CountT n_ep_send_tot_glb = PS::Comm::getSum(PS::n_ep_send_tot_loc);
        PS::CountT n_sp_send_tot_glb = PS::Comm::getSum(PS::n_sp_send_tot_loc);
        PS::CountT n_sp_send_tot_isendrecv_glb = PS::Comm::getSum(PS::n_sp_send_tot_isendrecv_loc);
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        if(PS::Comm::getRank() == 0){
            fout_debug<<"n_proc= "<<n_proc<<std::endl;
            fout_debug<<"<n_epi>= "<<TREE_POINTER->getNumberOfEPIAverage()<<std::endl;
            fout_debug<<"<n_epj>= "<<TREE_POINTER->getNumberOfEPJAverage()<<std::endl;
            fout_debug<<"<n_spj>= "<<TREE_POINTER->getNumberOfSPJAverage()<<std::endl;
            fout_debug<<"n_ipg= "<<TREE_POINTER->getNumberOfIPG()<<std::endl;
            fout_debug<<"n_ep_send_tot_glb= "<<n_ep_send_tot_glb
                      <<" n_ep_send_tot_glb/n_proc= "<<n_ep_send_tot_glb/n_proc
                      <<std::endl;
            fout_debug<<"n_sp_send_tot_glb= "<<n_sp_send_tot_glb
                      <<" n_sp_send_tot_glb/n_proc= "<<n_sp_send_tot_glb/n_proc
                      <<std::endl;
            fout_debug<<"n_sp_send_tot_isendrecv_glb= "<<n_sp_send_tot_isendrecv_glb
                      <<" n_sp_send_tot_isendrecv_glb/n_proc= "<<n_sp_send_tot_isendrecv_glb/n_proc
                      <<std::endl;
            fout_debug<<"band_width(ep, per proc)= "
                      <<(double)(n_ep_send_tot_glb*sizeof(Epj)) / PS::wtime_alltoallv_ep / 1e9 / n_proc
                      <<" [GB/s]"
                      <<std::endl;
            fout_debug<<"band_width(sp, per proc)= "
                      <<(double)(n_sp_send_tot_glb*sizeof(SPJMonopoleScatterSW)) / PS::wtime_alltoallv_sp / 1e9 / n_proc
                      <<" [GB/s]"
                      <<std::endl;
            fout_debug<<"band_width(sp(isend+irecv), per proc)= "
                      <<(double)(n_sp_send_tot_isendrecv_glb*sizeof(SPJMonopoleScatterSW)) / PS::wtime_alltoallv_sp_isendrecv / 1e9 / n_proc
                      <<" [GB/s]"
                      <<std::endl;            
            DumpProfileParticleSystem(prof_system, fout_debug);
            DumpProfileTree(prof_tree, fout_debug);
        }
        TREE_POINTER->clearTimeProfile();
        system.clearTimeProfile();
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
        // TIMING DUMP
        //////////////
        //PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "timing dump ends (CalcForceLoopMerge)" << std::endl;

        //////////////
        PS::Comm::barrier(); 
        fout_debug << "n_loc = " 
                   << PS::Comm::getRank() << " "
                   << system.getNumberOfParticleLocal() 
                   << std::endl;
        PS::Comm::barrier(); 
        //////////////

    }

    // For test run
    //fout_debug.close();
    //athread_halt();
    //PS::Finalize();
    //std::exit(0);
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
	fout_debug<<"n_proc, n_proc_quat="<<n_proc<<" "<<n_proc_quat<<std::endl;
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
	    fout_debug<<"nx_quat, ny_quat, nx, ny, nz= "
		      <<nx_quat<<" "<<ny_quat<<" "<<nx<<" "<<ny<<" "<<nz<<std::endl;
            */
	}
	dinfo.setNumberOfDomainMultiDimension(nx, ny, nz);
	PS::F64ort * pos_domain_tmp = new PS::F64ort[n_proc];
	if(PS::Comm::getRank() == 0){
	    fout_debug<<"n_glb_tmp="<<n_glb_tmp<<std::endl;
	    PS::S32 * istart = new PS::S32[n_proc_quat];
	    PS::S32 * iend   = new PS::S32[n_proc_quat];
	    PS::F64ort pos_root_domain_tmp = PS::F64ort( PS::F64vec(0.0, 0.0, -PS::LARGE_FLOAT), PS::F64vec(PS::LARGE_FLOAT, PS::LARGE_FLOAT, PS::LARGE_FLOAT));
            std::sort(pos_glb_tmp, pos_glb_tmp+n_glb_tmp, PS::LessOPX());
	    fout_debug<<"pos_glb_tmp[n_glb_tmp-1].x="<<pos_glb_tmp[n_glb_tmp-1].x<<std::endl;
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
