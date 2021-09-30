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
#include"nbody.hpp"

const bool FLAG_KEEP_FP_ORDER = false;
//const bool FLAG_KEEP_FP_ORDER = true;
//int n_loop_reuse = 1;
int n_loop_reuse = 8;

PS::F64 WTIME_KICK = 0.0;
PS::F64 WTIME_DRIFT = 0.0;
PS::F64 WTIME_KICK_DRIFT = 0.0;
PS::F64 FPGrav0::eps = 1.0/32.0;
PS::F64 FPGrav::eps = 1.0/32.0;

PS::F64 FPGrav0::dt_tree = 0.0;

PS::F64 WTIME_KERNEL = 0.0;
PS::F64 WTIME_SEND_IP = 0.0;
PS::F64 WTIME_COPY_IP = 0.0;
PS::F64 WTIME_SEND_ID = 0.0;
PS::F64 WTIME_COPY_ID = 0.0;
PS::F64 WTIME_SEND_JP = 0.0;
PS::F64 WTIME_COPY_JP = 0.0;
PS::F64 WTIME_RECV_FORCE = 0.0;
PS::F64 WTIME_COPY_FORCE = 0.0;
PS::F64 WTIME_COPY_IP_ID = 0.0;
PS::F64 WTIME_COPY_IP_JP = 0.0;
PS::F64 WTIME_H2D_IP     = 0.0;
PS::F64 WTIME_H2D_LIST   = 0.0;
PS::F64 WTIME_D2H_FORCE   = 0.0;
PS::F64 WTIME_H2D_ALL_PTCL = 0.0;
/*
PS::F64 WTIME_H2D = 0.0;
PS::F64 WTIME_D2H = 0.0;
PS::F64 WTIME_SEND_EPJ = 0.0;
PS::F64 WTIME_COPY_FOR_SEND_EPJ = 0.0;
PS::F64 WTIME_COPY_FOR_H2D = 0.0;
PS::F64 WTIME_COPY_FOR_D2H = 0.0;
PS::F64 WTIME_SEND_IP = 0.0;
PS::F64 WTIME_COPY_IP = 0.0;
PS::F64 WTIME_COMM_ID = 0.0;
PS::F64 WTIME_COPY_ID = 0.0;
*/

PS::S32 N_STEP = 0;
PS::F64 DT_TREE_HOST;
PS::S32 CONSTRUCTION_STEP = 1;

int main(int argc, char *argv[]) {
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

    PS::Initialize(argc, argv);

    //std::cerr<<"MPI_Wtick()= "<<MPI_Wtick()<<std::endl;

    const PS::S32 tag_max = 1;
    //const PS::S32 n_walk_limit = 200;
    const PS::S32 n_walk_limit = 1000;
    //const PS::S32 n_walk_limit = 2000;
    //const PS::S32 n_walk_limit = 4000;

    PS::F32 theta = 0.5;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::F32 dt = 1.0 / 128.0;
    //PS::F32 dt = 1.0 / 1024.0;
    PS::F32 dt_diag = 1.0 / 8.0;
    PS::F32 dt_snap = 1.0;
    //PS::F32 time_end = 10.0;
    //PS::F32 time_end = 1.0;
    //PS::F32 time_end = dt*32;
    //PS::F32 time_end = dt*16;
    //PS::F32 time_end = dt*64;
    PS::F32 time_end = dt*64;
    char dir_name[1024];
    PS::S64 n_tot = 1024;
    PS::S32 c;
    sprintf(dir_name,"./result");
    opterr = 0;
    while((c=getopt(argc,argv,"r:i:o:d:D:t:T:l:n:N:hs:")) != -1){
        switch(c){
        case 'r':
            n_loop_reuse = atoi(optarg);
            std::cerr << "n_loop_reuse =" << n_loop_reuse << std::endl;
            break;
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
    DT_TREE_HOST = dt;
    FPGrav0::dt_tree = dt;
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

    /////////// start new lines //////////
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_open();
    g5_set_eps_to_all(FPGrav::eps);
#endif

    PS::ParticleSystem<FPGrav0> system_grav0;
    system_grav0.initialize();
    PS::S32 n_smp = 200;
    system_grav0.setAverageTargetNumberOfSampleParticlePerProcess(n_smp);
    PS::S32 n_loc    = 0;
    PS::F32 time_sys = 0.0;
#if 0
    if(PS::Comm::getRank() == 0) {
        SetParticlesPlummer(system_grav0, n_tot, n_loc);
    } else {
        system_grav0.setNumberOfParticleLocal(n_loc);
    }
#else
    if(PS::Comm::getRank() == 0) {
        SetParticlesColdUniformSphere(system_grav0, n_tot, n_loc);
    } else {
        system_grav0.setNumberOfParticleLocal(n_loc);
    }
#endif
    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.decomposeDomainAll(system_grav0);
    system_grav0.exchangeParticle(dinfo);
    n_loc = system_grav0.getNumberOfParticleLocal();
    PS::TreeForForceLong<ForceGrav0, EPIGrav0, EPJGrav0>::Monopole tree_grav0;
    tree_grav0.initialize(n_tot, theta, n_leaf_limit, n_group_limit);
    N_STEP = 0;
    tree_grav0.calcForceAllAndWriteBack(CalcGravity<EPIGrav0, EPJGrav0, ForceGrav0>,
                                        CalcGravity<EPIGrav0, PS::SPJMonopole, ForceGrav0>,
                                        system_grav0,
                                        dinfo,
                                        true,
                                        PS::MAKE_LIST,
                                        FLAG_KEEP_FP_ORDER);
    n_loc    = system_grav0.getNumberOfParticleLocal();
    //PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
    //calcEnergy(system_grav0, Etot0, Ekin0, Epot0);
    Energy E0, E1;
    calcEnergy(system_grav0, E0.tot, E0.kin, E0.pot);
    if(PS::Comm::getRank()==0){
        std::cerr<<"E0.tot= "<<E0.tot
                 <<" E0.kin= "<<E0.kin
                 <<" E0.pot= "<<E0.pot
                 <<std::endl;
    }
    kick(system_grav0, dt * 0.5);
    drift(system_grav0, dt);
    tree_grav0.freeMem();
    /////////// end new lines //////////

    PS::ParticleSystem<FP> system_grav;
    system_grav.initialize();
    system_grav.setAverageTargetNumberOfSampleParticlePerProcess(n_smp);
    system_grav.setNumberOfParticleLocal(n_loc);
    for(PS::S32 i=0; i<n_loc; i++){
        system_grav[i].mass = system_grav0[i].mass;
        system_grav[i].pos  = system_grav0[i].pos;
        system_grav[i].vel  = system_grav0[i].vel;
    }
    system_grav0.freeMem();
    PS::TreeForForceLong<Force, EPI, EPJ>::Monopole tree_grav;
    tree_grav.initialize(n_tot, theta, n_leaf_limit, n_group_limit);


    PS::F64 time_diag = 0.0;
    PS::F64 time_snap = 0.0;
    PS::S64 n_loop = 0;
    PS::S32 id_snap = 0;
    PS::F64 wtime_loop = 0.0;
    while(time_sys < time_end){
        PS::INTERACTION_LIST_MODE int_mode = PS::MAKE_LIST_FOR_REUSE;
        PS::Comm::barrier();
        double wtime_0 = PS::GetWtime();
#ifdef REUSE_LIST_MODE
        if(n_loop % n_loop_reuse == 0){
            dinfo.decomposeDomainAll(system_grav);
            system_grav.exchangeParticle(dinfo);
            int_mode = PS::MAKE_LIST_FOR_REUSE;
            CONSTRUCTION_STEP = 1;
            N_STEP = 0;
        }
        else{
            int_mode = PS::REUSE_LIST;
            CONSTRUCTION_STEP = 0;
            N_STEP++;
        }
        double wtime_1 = PS::GetWtime();
    #ifdef MULTI_WALK
        tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
                                                     RetrieveKernel,
                                                     tag_max,
                                                     system_grav,
                                                     dinfo,
                                                     n_walk_limit,
                                                     true,
                                                     int_mode,
                                                     FLAG_KEEP_FP_ORDER);
    #elif defined(MULTI_WALK_INDEX)
        tree_grav.calcForceAllAndWriteBackMultiWalkIndex(DispatchKernelIndex,
                                                         RetrieveKernel,
                                                         tag_max,
                                                         system_grav,
                                                         dinfo,
                                                         n_walk_limit,
                                                         true,
                                                         int_mode,
                                                         FLAG_KEEP_FP_ORDER);
    #elif defined(MULTI_WALK_INDEX_2)
        tree_grav.calcForceAllAndWriteBackMultiWalkIndex(DispatchKernelIndex2,
                                                         RetrieveKernel2,
                                                         tag_max,
                                                         system_grav,
                                                         dinfo,
                                                         n_walk_limit,
                                                         true,
                                                         int_mode,
                                                         FLAG_KEEP_FP_ORDER);
    #elif defined(MULTI_WALK_INDEX_STREAM)
        tree_grav.calcForceAllAndWriteBackMultiWalkIndex2(DispatchKernelIndexStream2,
                                                          RetrieveKernelStream,
                                                          tag_max,
                                                          system_grav,
                                                          dinfo,
                                                          n_walk_limit,
                                                          true,
                                                          int_mode,
                                                          FLAG_KEEP_FP_ORDER);
    #else //MULTI_WALK
        tree_grav.calcForceAllAndWriteBack2(CalcGravity<EPI, EPJ, Force>,
                                            CalcGravity<EPI, PS::SPJMonopole, Force>,
                                            system_grav,
                                            dinfo,
                                            true,
                                            int_mode,
                                            FLAG_KEEP_FP_ORDER);
    #endif //MULTI_WALK
#else //REUSE_LIST_MODE
        if(n_loop % n_loop_reuse == 0){
            dinfo.decomposeDomainAll(system_grav);
        }
        system_grav.exchangeParticle(dinfo);
        double wtime_1 = PS::GetWtime();
    #ifdef MULTI_WALK
        tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
                                                     RetrieveKernel,
                                                    tag_max,
                                                    system_grav,
                                                    dinfo,
                                                    n_walk_limit,
                                                    true,
                                                    PS::MAKE_LIST,
                                                    FLAG_KEEP_FP_ORDER);

  #elif defined(MULTI_WALK_STREAM)
        tree_grav.calcForceAllAndWriteBackMultiWalk2(DispatchKernelStream,
                                                     RetrieveKernelStream,
                                                     tag_max,
                                                     system_grav,
                                                     dinfo,
                                                     n_walk_limit,
                                                     true,
                                                     PS::MAKE_LIST,
                                                     FLAG_KEEP_FP_ORDER);


    #elif defined(MULTI_WALK_INDEX)
        tree_grav.calcForceAllAndWriteBackMultiWalkIndex(DispatchKernelIndex,
                                                         RetrieveKernel,
                                                         tag_max,
                                                         system_grav,
                                                         dinfo,
                                                         n_walk_limit,
                                                         true,
                                                         PS::MAKE_LIST,
                                                         FLAG_KEEP_FP_ORDER);

    #elif defined(MULTI_WALK_INDEX_STREAM)
        tree_grav.calcForceAllAndWriteBackMultiWalkIndex2(DispatchKernelIndexStream2,
                                                         RetrieveKernelStream,
                                                         tag_max,
                                                         system_grav,
                                                         dinfo,
                                                         n_walk_limit,
                                                         true,
                                                         PS::MAKE_LIST,
                                                         FLAG_KEEP_FP_ORDER);

    #else //MULTI_WALK
        tree_grav.calcForceAllAndWriteBack2(CalcGravity<EPI, EPJ, Force>,
                                           CalcGravity<EPI, PS::SPJMonopole, Force>,
                                           system_grav,
                                           dinfo,
                                           true,
                                           PS::MAKE_LIST,
                                           FLAG_KEEP_FP_ORDER);
    #endif //MULTI_WALK
#endif // REUSE_LIST
        if(n_loop % n_loop_reuse == 0){
            std::cout<<"BIG STEP START"<<std::endl;
        }
#if 0
    #ifdef CHECK_ENERGY
        kick(system_grav, dt * 0.5);
        calcEnergy(system_grav, E1.tot, E1.kin, E1.pot);
        std::cerr<<"E1.tot= "<<E1.tot<<" E1.kin= "<<E1.kin<<" E1.pot= "<<E1.pot
                 <<" rel_err= "<<(E1.tot - E0.tot) / E0.tot
                 <<std::endl;
        kick(system_grav, dt * 0.5);
        drift(system_grav, dt);
    #else
        kick_drift(system_grav, dt);
    #endif
#endif
        time_sys += dt;
        PS::Comm::barrier();
        wtime_loop = PS::GetWtime()-wtime_0;
        DumpTimeProfile(system_grav, dinfo, tree_grav, wtime_loop, n_loop);
        if(n_loop==0){
            tree_grav.dumpIpg(fout_eng);
        }
        n_loop++;
    }
    
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_close();
#endif

    PS::Finalize();
    return 0;
}
