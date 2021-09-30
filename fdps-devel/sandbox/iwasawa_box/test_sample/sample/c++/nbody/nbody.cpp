#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>
#ifdef ENABLE_PHANTOM_GRAPE_X86
#include <gp5util.h>
#endif
#ifdef ENABLE_GPU_CUDA
#include"force_gpu_cuda.hpp"
#endif
#include "user-defined.hpp"

/*
// for debug
PS::S32 DispatchKernelIndex(const PS::S32 tag,
                            const PS::S32 n_walk,
                            const FPGrav ** epi,
                            const PS::S32 *  n_epi,
                            const PS::S32 ** id_epj,
                            const PS::S32 *  n_epj,
                            const PS::S32 ** id_spj,
                            const PS::S32 *  n_spj,
                            const FPGrav * epj,
                            const PS::S32 n_epj_tot,
                            const PS::SPJMonopole * spj,
                            const PS::S32 n_spj_tot,
                            const bool send_flag){
    if(send_flag == true){
        return 0;
    }
    else{
        return 0;
    }
}

PS::S32 RetrieveKernelIndex(const PS::S32 tag,
                            const PS::S32 n_walk,
                            const PS::S32 * ni,
                            FPGrav      ** force){
    return 0;
}
*/

void MakePlummerModel(const double mass_glb,
                      const long long int n_glb,
                      const long long int n_loc,
                      double *& mass,
                      PS::F64vec *& pos,
                      PS::F64vec *& vel,
                      const double eng = -0.25,
                      const int seed = 0){

    assert(eng < 0.0);
    static const double PI = atan(1.0) * 4.0;
    const double r_cutoff = 22.8 / (-3.0 * PI * mass_glb * mass_glb / (64.0 * -0.25)); // 22.8 is cutoff in Nbody units
    mass = new double[n_loc];
    pos = new PS::F64vec[n_loc];
    vel = new PS::F64vec[n_loc];

    PS::MTTS mt;
    mt.init_genrand( PS::Comm::getRank() );
    for(int i=0; i<n_loc; i++){
        mass[i] = mass_glb / n_glb;
        double r_tmp = 9999.9;
        while(r_tmp > r_cutoff){ 
            double m_tmp = mt.genrand_res53();
            r_tmp = 1.0 / sqrt( pow(m_tmp, (-2.0/3.0)) - 1.0);
        }
        double phi = 2.0 * PI * mt.genrand_res53();
        double cth = 2.0 * (mt.genrand_real2() - 0.5);
        double sth = sqrt(1.0 - cth*cth);
        pos[i][0] = r_tmp * sth * cos(phi);
        pos[i][1] = r_tmp * sth * sin(phi);
        pos[i][2] = r_tmp * cth;
        while(1){
            const double v_max = 0.1;
            const double v_try = mt.genrand_res53();
            const double v_crit = v_max * mt.genrand_res53();
            if(v_crit < v_try * v_try * pow( (1.0 - v_try * v_try), 3.5) ){
                const double ve = sqrt(2.0) * pow( (r_tmp*r_tmp + 1.0), -0.25 );
                phi = 2.0 * PI * mt.genrand_res53();
                cth = 2.0 * (mt.genrand_res53() - 0.5);
                sth = sqrt(1.0 - cth*cth);
                vel[i][0] = ve * v_try * sth * cos(phi);
                vel[i][1] = ve * v_try * sth * sin(phi);
                vel[i][2] = ve * v_try * cth;
                break;
            }
        }
    }

    PS::F64vec cm_pos = 0.0;
    PS::F64vec cm_vel = 0.0;
    double  cm_mass = 0.0;
    for(size_t i=0; i<n_loc; i++){
        cm_pos += mass[i] * pos[i];
        cm_vel += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
    for(size_t i=0; i<n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }

    const double r_scale = -3.0 * PI * mass_glb * mass_glb / (64.0 * eng);
    const double coef = 1.0 / sqrt(r_scale);
    for(size_t i=0; i<n_loc; i++){
        pos[i] *= r_scale;
        vel[i] *= coef;
    }

    double r_max_sq = -1.0;
    for(int i=0; i<n_loc; i++){
        if(r_max_sq < pos[i] * pos[i]){
            r_max_sq = pos[i] * pos[i];
        }
    }
    std::cout<<"r_max= "<<sqrt(r_max_sq)<<std::endl;
}


template<class Tpsys>
void SetParticlesPlummer(Tpsys & psys,
                         const PS::S64 n_glb,
                         PS::S32 & n_loc,  
                         PS::F32 & t_sys=0.0){

    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    /*
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
    */
    n_loc = n_glb;
    psys.setNumberOfParticleLocal(n_loc);

    PS::F64 * mass;
    PS::F64vec * pos;
    PS::F64vec * vel;
    t_sys = 0.0;

    const PS::F64 m_tot = 1.0;
    const PS::F64 eng = -0.25;
    MakePlummerModel(m_tot, n_glb, n_loc, mass, pos, vel, eng);
    //PS::S64 i_h = n_glb/n_proc*my_rank;
    //if( n_glb % n_proc  > my_rank) i_h += my_rank;
    //else i_h += n_glb % n_proc;
    for(size_t i=0; i<n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        //psys[i].id = i_h + i;
        psys[i].id = i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}

void makeColdUniformSphere(const PS::F64 mass_glb,
                           const PS::S64 n_glb,
                           const PS::S64 n_loc,
                           PS::F64 *& mass,
                           PS::F64vec *& pos,
                           PS::F64vec *& vel,
                           const PS::F64 eng = -0.25,
                           const PS::S32 seed = 0) {
    
    assert(eng < 0.0);
    {
        PS::MTTS mt;
        mt.init_genrand(0);
        for(PS::S32 i = 0; i < n_loc; i++){
            mass[i] = mass_glb / n_glb;
            const PS::F64 radius = 3.0;
            do {
                pos[i][0] = (2. * mt.genrand_res53() - 1.) * radius;
                pos[i][1] = (2. * mt.genrand_res53() - 1.) * radius;
                pos[i][2] = (2. * mt.genrand_res53() - 1.) * radius;
            }while(pos[i] * pos[i] >= radius * radius);
            vel[i][0] = 0.0;
            vel[i][1] = 0.0;
            vel[i][2] = 0.0;
        }
    }

    PS::F64vec cm_pos  = 0.0;
    PS::F64vec cm_vel  = 0.0;
    PS::F64    cm_mass = 0.0;
    for(PS::S32 i = 0; i < n_loc; i++){
        cm_pos  += mass[i] * pos[i];
        cm_vel  += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
    for(PS::S32 i = 0; i < n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }
}

template<class Tpsys>
void setParticlesColdUniformSphere(Tpsys & psys,
                                   const PS::S32 n_glb,
                                   PS::S32 & n_loc) {

    n_loc = n_glb; 
    psys.setNumberOfParticleLocal(n_loc);

    PS::F64    * mass = new PS::F64[n_loc];
    PS::F64vec * pos  = new PS::F64vec[n_loc];
    PS::F64vec * vel  = new PS::F64vec[n_loc];
    const PS::F64 m_tot = 1.0;
    const PS::F64 eng   = -0.25;
    makeColdUniformSphere(m_tot, n_glb, n_loc, mass, pos, vel, eng);
    for(PS::S32 i = 0; i < n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos  = pos[i];
        psys[i].vel  = vel[i];
        psys[i].id   = i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}

template<class Tpsys>
void kick(Tpsys & system,
          const PS::F64 dt) {
    PS::S32 n = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < n; i++) {
        system[i].vel  += system[i].acc * dt;
    }
}

template<class Tpsys>
void drift(Tpsys & system,
           const PS::F64 dt) {
    PS::S32 n = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < n; i++) {
        system[i].pos  += system[i].vel * dt;
    }
}

template<class Tpsys>
void calcEnergy(const Tpsys & system,
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
    for(PS::S32 i = 0; i < nbody; i++){
        ekin_loc += system[i].mass * system[i].vel * system[i].vel;
        epot_loc += system[i].mass * (system[i].pot + system[i].mass / FPGrav::eps);
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    etot_loc  = ekin_loc + epot_loc;
    etot = PS::Comm::getSum(etot_loc);
    epot = PS::Comm::getSum(epot_loc);
    ekin = PS::Comm::getSum(ekin_loc);    
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
    PS::S32 ret;
    if (PS::Comm::getRank() == 0) {
        if (stat(dir_name, &st) != 0) {
            ret = mkdir(dir_name, 0777);
        } else {
            ret = 0; // the directory named dir_name already exists.
        }
    } 
    PS::Comm::broadcast(&ret, 1);
    if (ret == 0) {
        if (PS::Comm::getRank() == 0)
            fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
    } else {
        if (PS::Comm::getRank() == 0)
            fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
        PS::Abort();
    }
}

PS::F64 FPGrav::eps = 1.0/32.0;
PS::F64 WTIME_KERNEL = 0.0;
PS::F64 WTIME_H2D = 0.0;
PS::F64 WTIME_D2H = 0.0;
PS::F64 WTIME_SEND_EPJ = 0.0;
int main(int argc, char *argv[]) {
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

    PS::Initialize(argc, argv);
    const PS::S32 tag_max = 1;
    const PS::S32 n_walk_limit = 200;
    PS::F32 theta = 0.5;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    //PS::F32 time_end = 10.0;
    PS::F32 time_end = 1.0 / 64.0;
    PS::F32 dt = 1.0 / 128.0;
    PS::F32 dt_diag = 1.0 / 8.0;
    PS::F32 dt_snap = 1.0;
    char dir_name[1024];
    PS::S64 n_tot = 1024;
    PS::S32 c;
    sprintf(dir_name,"./result");
    opterr = 0;
    while((c=getopt(argc,argv,"i:o:d:D:t:T:l:n:N:hs:")) != -1){
        switch(c){
        case 'o':
            sprintf(dir_name,optarg);
            break;
        case 't':
            theta = atof(optarg);
            std::cerr << "theta =" << theta << std::endl;
            break;
        case 'T':
            time_end = atof(optarg);
            std::cerr << "time_end = " << time_end << std::endl;
            break;
        case 's':
            dt = atof(optarg);
            std::cerr << "time_step = " << dt << std::endl;
            break;
        case 'd':
            dt_diag = atof(optarg);
            std::cerr << "dt_diag = " << dt_diag << std::endl;
            break;
        case 'D':
            dt_snap = atof(optarg);
            std::cerr << "dt_snap = " << dt_snap << std::endl;
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            std::cerr << "n_leaf_limit = " << n_leaf_limit << std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            std::cerr << "n_group_limit = " << n_group_limit << std::endl;
            break;
        case 'N':
            n_tot = atoi(optarg);
            std::cerr << "n_tot = " << n_tot << std::endl;
            break;
        case 'h':
            if(PS::Comm::getRank() == 0) {
                printHelp();
            }
            PS::Finalize();
            return 0;
        default:
            if(PS::Comm::getRank() == 0) {
                std::cerr<<"No such option! Available options are here."<<std::endl;
                printHelp();
            }
            PS::Abort();
        }
    }

    makeOutputDirectory(dir_name);

    std::ofstream fout_eng;

    if(PS::Comm::getRank() == 0) {
        char sout_de[1024];
	sprintf(sout_de, "%s/t-de.dat", dir_name);
        fout_eng.open(sout_de);
        fprintf(stdout, "This is a sample program of N-body simulation on FDPS!\n");
        fprintf(stdout, "Number of processes: %d\n", PS::Comm::getNumberOfProc());
        fprintf(stdout, "Number of threads per process: %d\n", PS::Comm::getNumberOfThread());
    }

    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
    PS::S32 n_loc    = 0;
    PS::F32 time_sys = 0.0;
    if(PS::Comm::getRank() == 0) {
        //setParticlesColdUniformSphere(system_grav, n_tot, n_loc);
        SetParticlesPlummer(system_grav, n_tot, n_loc, time_sys);
    } else {
        system_grav.setNumberOfParticleLocal(n_loc);
    }

    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.decomposeDomainAll(system_grav);
    system_grav.exchangeParticle(dinfo);
    n_loc = system_grav.getNumberOfParticleLocal();
    
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_open();
    g5_set_eps_to_all(FPGrav::eps);
#endif
    
    PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole tree_grav;
    tree_grav.initialize(n_tot, theta, n_leaf_limit, n_group_limit);
#ifdef REUSE_LIST_MODE
    #ifdef MULTI_WALK
    //const PS::S32 n_walk_limit = 200;
    //const PS::S32 tag_max = 1;
    tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
                                                RetrieveKernel,
                                                tag_max,
                                                system_grav,
                                                dinfo,
                                                n_walk_limit,
                                                true,
                                                PS::MAKE_LIST_FOR_REUSE);
    #elif MULTI_WALK_INDEX
    //const PS::S32 n_walk_limit = 200;
    //const PS::S32 tag_max = 1;
    tree_grav.calcForceAllAndWriteBackMultiWalkIndex(DispatchKernelIndex,
                                                     RetrieveKernel,
                                                     tag_max,
                                                     system_grav,
                                                     dinfo,
                                                     n_walk_limit,
                                                     true,
                                                     PS::MAKE_LIST_FOR_REUSE);
    #else //MULTI_WALK
    tree_grav.calcForceAllAndWriteBack(CalcGravity<FPGrav>,
                                       CalcGravity<PS::SPJMonopole>,
                                       system_grav,
                                       dinfo,
                                       true,
                                       PS::MAKE_LIST_FOR_REUSE);
    #endif //MULTI_WALK
    /*
    PS::Comm::barrier();
    if(PS::Comm::getRank()==0) std::cerr<<"OK"<<std::endl;
    PS::Comm::barrier();
    exit(1);
    */
#else //REUSE_LIST_MODE
    
    #ifdef MULTI_WALK
    //const PS::S32 n_walk_limit = 200;
    //const PS::S32 tag_max = 1;
    tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
                                                RetrieveKernel,
                                                tag_max,
                                                system_grav,
                                                dinfo,
                                                n_walk_limit);
    #elif MULTI_WALK_INDEX
    //const PS::S32 n_walk_limit = 200;
    //const PS::S32 tag_max = 1;
    tree_grav.calcForceAllAndWriteBackMultiWalkIndex(DispatchKernelIndex,
                                                     RetrieveKernel,
                                                     tag_max,
                                                     system_grav,
                                                     dinfo,
                                                     n_walk_limit,
                                                     true,
                                                     PS::MAKE_LIST);

    #else //MULTI_WALK
    tree_grav.calcForceAllAndWriteBack(CalcGravity<FPGrav>,
                                       CalcGravity<PS::SPJMonopole>,
                                       system_grav,
                                       dinfo);
    #endif //MULTI_WALK
#endif //REUSE_LIST_MODE
    PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
    calcEnergy(system_grav, Etot0, Ekin0, Epot0);
    if(PS::Comm::getRank()==0){
        std::cerr<<"Epot0= "<<Epot0
                 <<" Ekin0= "<<Ekin0
                 <<" Etot0= "<<Etot0
                 <<std::endl;
    }
    //PS::Comm::barrier();
    //exit(1);

#if 1
    PS::F64 time_diag = 0.0;
    PS::F64 time_snap = 0.0;
    PS::S64 n_loop = 0;
    PS::S32 id_snap = 0;
    PS::F64 wtime_loop = 0.0;

    while(time_sys < time_end){
        if( (time_sys >= time_snap) || ( (time_sys + dt) - time_snap ) > (time_snap - time_sys) ){
            char filename[256];
            sprintf(filename, "%s/%04d.dat", dir_name, id_snap++);
            FileHeader header;
            header.time   = time_sys;
            header.n_body = system_grav.getNumberOfParticleGlobal();
            system_grav.writeParticleAscii(filename, header);
            time_snap += dt_snap;
        }

        calcEnergy(system_grav, Etot1, Ekin1, Epot1);
        PS::TimeProfile tp = tree_grav.getTimeProfile();
        PS::CountT n_int_ep_ep = tree_grav.getNumberOfInteractionEPEPGlobal();
        PS::CountT n_int_ep_sp = tree_grav.getNumberOfInteractionEPSPGlobal();
        if(PS::Comm::getRank() == 0){
            //if( (time_sys >= time_diag) || ( (time_sys + dt) - time_diag ) > (time_diag - time_sys) ){
            if(1){
                tp.dump(std::cout);
                std::cout<<"tp.calc_force__core= "<<tp.calc_force__core
                         <<" tp.calc_force__copy_original_order= "<<tp.calc_force__copy_original_order
                         <<" tp.add_moment_as_sp_local= "<<tp.add_moment_as_sp_local
                         <<" tp.add_moment_as_sp_global= "<<tp.add_moment_as_sp_global
                         <<std::endl;
                std::cout<<"WTIME_KERNEL= "<<WTIME_KERNEL
                         <<" WTIME_H2D= "<<WTIME_H2D
                         <<" WTIME_D2H= "<<WTIME_D2H
                         <<" WTIME_SEND_EPJ= "<<WTIME_SEND_EPJ
                         <<std::endl;
                std::cout<<"n_int_ep_ep= "<<n_int_ep_ep
                         <<" n_int_ep_sp= "<<n_int_ep_sp
                         <<std::endl;
                std::cout<<"speed(kernel only)= "<<double((n_int_ep_ep+n_int_ep_sp)*38.0) / WTIME_KERNEL * 1e-9 <<" [Gflops], speed(total)= "
                         <<double((n_int_ep_ep+n_int_ep_sp)*38.0) / wtime_loop * 1e-9<<" [Gflops]"
                         <<std::endl;

                fout_eng << time_sys << "   " << (Etot1 - Etot0) / Etot0 << std::endl;
                fprintf(stdout, "time: %10.7f energy error: %+e wtime_loop: %+e\n",
                        time_sys, (Etot1 - Etot0) / Etot0, wtime_loop);
                time_diag += dt_diag;
            }
        }
        WTIME_KERNEL = WTIME_H2D = WTIME_D2H = WTIME_SEND_EPJ = 0.0;
        tree_grav.clearTimeProfile();
        tree_grav.clearCounterAll();

        kick(system_grav, dt * 0.5);
        
        time_sys += dt;
        drift(system_grav, dt);

#ifdef REUSE_LIST_MODE
        /*
        if(PS::Comm::getRank()==0){
            std::cerr<<"Epot1= "<<Epot1
                     <<" Ekin1= "<<Ekin1
                     <<" Etot1= "<<Etot1
                     <<std::endl;
            std::cerr<<"A) system_grav[[0].pos= "<<system_grav[0].pos
                     <<" pot= "<<system_grav[0].pot
                     <<std::endl;
        }
        PS::Comm::barrier();
        */
    #ifdef MULTI_WALK
        if(n_loop % 4 == 0){
            dinfo.decomposeDomainAll(system_grav);
            system_grav.exchangeParticle(dinfo);
            tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
                                                        RetrieveKernel,
                                                        tag_max,
                                                        system_grav,
                                                        dinfo,
                                                        n_walk_limit,
                                                        true,
                                                        PS::MAKE_LIST_FOR_REUSE);
        }
        else{
            tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
                                                        RetrieveKernel,
                                                        tag_max,
                                                        system_grav,
                                                        dinfo,
                                                        n_walk_limit,
                                                        true,
                                                        PS::REUSE_LIST);
        }
    #elif MULTI_WALK_INDEX
        double wtime_0 = PS::GetWtime();
        if(n_loop % 4 == 0){
            dinfo.decomposeDomainAll(system_grav);
            system_grav.exchangeParticle(dinfo);
            tree_grav.calcForceAllAndWriteBackMultiWalkIndex(DispatchKernelIndex,
                                                             RetrieveKernel,
                                                             tag_max,
                                                             system_grav,
                                                             dinfo,
                                                             n_walk_limit,
                                                             true,
                                                             PS::MAKE_LIST_FOR_REUSE);
        }
        else{
            tree_grav.calcForceAllAndWriteBackMultiWalkIndex(DispatchKernelIndex,
                                                             RetrieveKernel,
                                                             tag_max,
                                                             system_grav,
                                                             dinfo,
                                                             n_walk_limit,
                                                             true,
                                                             PS::REUSE_LIST);
        }
        wtime_loop = PS::GetWtime()-wtime_0;
        //std::cerr<<"wtime= "<<wtime_loop<<std::endl;
    #else //MULTI_WALK
        if(n_loop % 4 == 0){
        //if(1){
            dinfo.decomposeDomainAll(system_grav);
            system_grav.exchangeParticle(dinfo);
            tree_grav.calcForceAllAndWriteBack(CalcGravity<FPGrav>,
                                               CalcGravity<PS::SPJMonopole>,
                                               system_grav,
                                               dinfo,
                                               true,
                                               PS::MAKE_LIST_FOR_REUSE);
        }
        else{
            tree_grav.calcForceAllAndWriteBack(CalcGravity<FPGrav>,
                                               CalcGravity<PS::SPJMonopole>,
                                               system_grav,
                                               dinfo,
                                               true,
                                               PS::REUSE_LIST);
        }
    #endif //MULTI_WALK
        
#else //REUSE_LIST_MODE
        if(n_loop % 4 == 0){
            dinfo.decomposeDomainAll(system_grav);
        }
        system_grav.exchangeParticle(dinfo);
    #ifdef MULTI_WALK
        e;
        tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
                                                    RetrieveKernel,
                                                    tag_max,
                                                    system_grav,
                                                    dinfo,
                                                    n_walk_limit,
                                                    true);

    #elif MULTI_WALK_INDEX
    tree_grav.calcForceAllAndWriteBackMultiWalkIndex(DispatchKernelIndex,
                                                     RetrieveKernel,
                                                     tag_max,
                                                     system_grav,
                                                     dinfo,
                                                     n_walk_limit,
                                                     true,
                                                     PS::MAKE_LIST);

    #else //MULTI_WALK
        tree_grav.calcForceAllAndWriteBack(CalcGravity<FPGrav>,
                                           CalcGravity<PS::SPJMonopole>,
                                           system_grav,
                                           dinfo);
    #endif //MULTI_WALK
#endif // REUSE_LIST
        
        kick(system_grav, dt * 0.5);
        
        n_loop++;
    }
    
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_close();
#endif

#endif

    PS::Finalize();
    return 0;
}
