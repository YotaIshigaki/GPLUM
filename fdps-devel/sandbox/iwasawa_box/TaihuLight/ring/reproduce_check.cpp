#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/time.h>
#include<cmath>
#include<particle_simulator.hpp>
#include"debug_tools.hpp"
#include"user-defined.hpp"
#include"force_sunway_impl.hpp"
#include"ic.hpp"
//#include"collision.hpp"
#include"nbody.hpp"
#include"get_pos_domain.hpp"

#ifdef SUNWAY
extern "C"{
    #include <athread.h>
    #include "cpe_func.h"
    extern void SLAVE_FUN(Rotate)(void*);
    extern void SLAVE_FUN(RotateCyl)(void*);
    extern void SLAVE_FUN(CopyEpiToFpWithCoordinateTransformCPE)(void*);
    extern void SLAVE_FUN(CopyEPIToFP)(void*);
    extern void SLAVE_FUN(RotateAndShiftZ)(void*);
}
#endif


static double wtime(){
   struct timeval tv;
   gettimeofday(&tv, NULL);
   return 1.e-6*tv.tv_usec +tv.tv_sec;
}

char HOST_NAME[170000][64];

MPI_Comm MY_MPI_COMM_SUB;
MPI_Comm MY_MPI_COMM_1D;
PS::S32 MY_RANK_SUB;
PS::S32 N_LOOP_GLB;

PS::F64 MASS_DUST_TOTAL = 0.0;
PS::F64vec CM_POS_DUST_TOTAL = 0.0;
PS::F64 Epi::eps      = 1.0/32.0;
PS::F64 Epj::r_coll   = 0.5;
PS::F64 Epj::r_search = 0.5;
PS::F64 Planet::r_coll = 0.5; // identical to Satellite::r_coll; tetative value; modified later.
//PS::F64 DT_TREE = 1.0 / 64.0;
//PS::F64 DT_TREE = 1.0 / 128.0;
PS::F64 DT_TREE = 1.0 / 256.0;
//PS::F64 DT_TREE = 1.0 / 512.0;
//PS::F64 DT_TREE = 1.0 / 2048.0;

Planet PLANET(1.0, PS::F64vec(0.0), PS::F64vec(0.0));

#ifdef DEBUG_CALC_FORCE
//PS::S32 N_SATELLITE = 0;
PS::S32 N_SATELLITE = 128;
#else
//PS::S32 N_SATELLITE = 0;
//PS::S32 N_SATELLITE = 1000;
//PS::S32 N_SATELLITE = 64;
PS::S32 N_SATELLITE = 128;
#endif

PS::ReallocatableArray<Satellite> SATELLITE;
PS::ReallocatableArray<SatelliteForce> SATELLITE_FORCE;
PS::ReallocatableArray<SatelliteForce> SATELLITE_FORCE_LOC;

std::ofstream fout_debug;
std::ofstream fout_wtime;
std::ofstream fout_host_name;

std::ofstream fout_ex_ptcl_max;
std::ofstream fout_ex_ptcl_min;
std::ofstream fout_calc_force_1st_max;
std::ofstream fout_calc_force_1st_min;

std::string s_ex_ptcl_max;
std::string s_ex_ptcl_min;
std::string s_calc_force_1st_max;
std::string s_calc_force_1st_min;


int main(int argc, char *argv[]) {

    double wtime_start_of_main = wtime();
    double wtime_main_0;
    unsigned long args[6];
    
    PS::Initialize(argc, argv);
    PS::Comm::barrier();
    LapTimer::initialize();
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef USE_MEMORY_POOL
    size_t size_memrory_pool = 1200000000;
    PS::MemoryPool::initialize(size_memrory_pool);
#endif
    
    LapTimer::initialize();

    if (PS::Comm::getRank() == 0)
       std::cerr << "etime(MPI startup) = " << wtime() - wtime_start_of_main << std::endl;

#ifdef SUNWAY
    wtime_main_0 = wtime();
    athread_init();
    if (PS::Comm::getRank() == 0)
       std::cerr << "etime(athread_init) = " << wtime() - wtime_main_0 << std::endl;
#endif
    
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    
    PS::S32 snp_id = 0;
    const PS::F64 PI = 4.0 * atan(1.0);
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::S32 my_rank = PS::Comm::getRank();

    MY_RANK_MPI = my_rank;
    
#ifdef COLLISION_TEST
    const PS::S32 n_loop_merge = 1;
#else
    //const PS::S32 n_loop_merge = 1;
    //const PS::S32 n_loop_merge = 8;
    //const PS::S32 n_loop_merge = 16;
    //const PS::S32 n_loop_merge = 32;
    const PS::S32 n_loop_merge = 64;
    //const PS::S32 n_loop_merge = 4;
#endif
    //const PS::S32 n_smp_ave = 128;
    //const PS::S32 n_smp_ave = 200;
    const PS::S32 n_smp_ave = 512;
    //const PS::S32 n_smp_ave = 2048;
    //const PS::S32 n_smp_ave = 8192;
    
    PS::F32 theta = 0.5;
    //PS::S32 n_leaf_limit = 8;
    PS::S32 n_leaf_limit = 16;
    PS::S32 n_group_limit = 64;
    //PS::F32 dt = 1.0 / 64.0;
    //PS::F32 time_end = DT_TREE*32;
    //PS::F32 time_end = DT_TREE*64;
    //PS::F32 time_end = DT_TREE*128;
    //PS::F32 time_end = DT_TREE*320; // performance measurement
    //PS::F32 time_end = DT_TREE*512;
    //PS::F32 time_end = DT_TREE*1024;
    //PS::F32 time_end = DT_TREE*2048; // 2.5 Kepler periods
    //PS::F32 time_end = DT_TREE*4096;
    //PS::F32 time_end = DT_TREE*8192; // 10 Kepler periods
    //PS::F32 time_end = DT_TREE*32768; // 40 Kepler periods
    PS::F32 time_end = 1.0; // for performance measurement
    //PS::F32 time_end = 1.5;
    //PS::F32 time_end = 2.0;
    PS::F32 dt_diag = 1.0 / 8.0;
    PS::F32 dt_snap = 1.0;
    char dir_name[1024];
    //PS::S64 n_glb = 1024, n_loc;
    PS::S64 n_glb = 1024;
    PS::S64 n_loc = 0;
    PS::F64 ax_cen = 1.0;
    //PS::F64 delta_ax = 0.1;
    //PS::F64 delta_ax = 0.01;
    PS::F64 delta_ax = 1e-3;
    //PS::F64 delta_ax = 1e-4;
    //PS::F64 delta_ax = 1e-5;
    PS::F64 ax_in  = ax_cen - 0.5*delta_ax;
    PS::F64 ax_out = ax_cen + 0.5*delta_ax;
    PS::F64 ecc_rms = 0.0;
    PS::F64 inc_rms = 16.0;
#ifdef ENERGY_CHECK
    PS::F64 tau = 1.0e-6;
#else
    //PS::F64 tau = 4.0;
    PS::F64 tau = 1.0;
    //PS::F64 tau = 1.0e-6;
#endif
    PS::F64 ratio_r_phy_r_hill = 1.0;
    //PS::F64 ratio_r_phy_r_hill = 1000.0;
    PS::S32 c;
    PS::F64ort pos_domain;
    PS::F64 theta_rot;
    PS::F64 starttime,endtime;

    sprintf(dir_name,"./result");
    opterr = 0;
    while((c=getopt(argc,argv,"i:o:d:D:t:T:l:n:N:hs:")) != -1){
        switch(c){
        case 'o':
            sprintf(dir_name,optarg);
            break;
        case 't':
            theta = atof(optarg);
            //std::cerr << "theta =" << theta << std::endl;
            break;
        case 'T':
            time_end = atof(optarg);
            //std::cerr << "time_end = " << time_end << std::endl;
            break;
        case 's':
            DT_TREE = atof(optarg);
            //std::cerr << "time_step = " << DT_TREE << std::endl;
            break;
        case 'd':
            dt_diag = atof(optarg);
            //std::cerr << "dt_diag = " << dt_diag << std::endl;
            break;
        case 'D':
            dt_snap = atof(optarg);
            //std::cerr << "dt_snap = " << dt_snap << std::endl;
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            //std::cerr << "n_leaf_limit = " << n_leaf_limit << std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            //std::cerr << "n_group_limit = " << n_group_limit << std::endl;
            break;
        case 'N':
            //n_glb = atol(optarg);
            n_glb = atol(optarg) * PS::Comm::getNumberOfProc();
            //std::cerr << "n_glb = " << n_glb << std::endl;
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

    //* Make a directory
    wtime_main_0 = wtime();
    //makeOutputDirectory(dir_name);
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0)
       std::cerr << "etime(mkdir) = " << wtime() - wtime_main_0 << std::endl;

    //* Open fout_debug
    wtime_main_0 = wtime();
    std::string s_debug = dir_name + ("/debug_" + ToS(PS::Comm::getRank()) + ".dat");
    if (PS::Comm::getRank() == 0) {
        std::cerr<<"s_debug: "<<s_debug<<std::endl;
        fout_debug.open(s_debug.c_str());
        if(!fout_debug.is_open()) exit(1);
        std::cerr << "etime(fopen) = " << wtime() - wtime_main_0 << std::endl;
    }
    std::string s_wtime = dir_name + ("/wtime_" + ToS(PS::Comm::getRank()) + ".dat");
    if (PS::Comm::getRank() == 0) {
        std::cerr<<"s_wtime: "<<s_wtime<<std::endl;
        fout_wtime.open(s_wtime.c_str());
        if(!fout_wtime.is_open()) exit(1);
        std::cerr << "etime(fopen) = " << wtime() - wtime_main_0 << std::endl;
    }
    

    std::string s_host_name = dir_name + ("/host_name_" + ToS(PS::Comm::getRank()) + ".dat");
    if (PS::Comm::getRank() == 0) {
        std::cerr<<"s_host_name: "<<s_host_name<<std::endl;
        fout_host_name.open(s_host_name.c_str());
        if(!fout_host_name.is_open()) exit(1);
        std::cerr << "etime(fopen) = " << wtime() - wtime_main_0 << std::endl;
    }

    char host_name_tmp[64];
    gethostname(host_name_tmp, 64);
    PS::Comm::barrier();
    MPI_Allgather(host_name_tmp, 64, MPI_CHAR,
                  HOST_NAME[0], 64, MPI_CHAR, MPI_COMM_WORLD);
    /*
    if (my_rank == 0) {
        for (PS::S32 i_rank=0; i_rank<n_proc; i_rank++) {
            fout_host_name << "RANK= " << i_rank
                           << " host= " <<HOST_NAME[i_rank]
                           << std::endl;
        }
        fout_host_name.flush();
    }
    */
    
    s_ex_ptcl_max = std::string(dir_name) + "/ex_ptcl_max.dat";
    s_ex_ptcl_min = std::string(dir_name) + "/ex_ptcl_min.dat";
    s_calc_force_1st_max = std::string(dir_name) + "/calc_force_force_1st_max.dat";
    s_calc_force_1st_min = std::string(dir_name) + "/calc_force_force_1st_min.dat";
    
    //* Open fout_eng
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

#ifdef PHI_R_TREE
    PS::S32 nx, ny, nz;
    #if 0
    //ny = 4;
    ny = 1;
    nz = 1;
    nx = n_proc / ny;
    #else
    DivideNProc(nx, ny, nz, n_proc, delta_ax);
    #endif
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

    PS::F64 r_phy = sqrt(tau*(ax_out*ax_out - ax_in*ax_in) / n_glb);
    PS::F64 r_hill = r_phy / ratio_r_phy_r_hill;
    #ifdef ENERGY_CHECK
    PS::F64 mass_dust = (r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*3.0*PLANET.mass*0.5 * 1e14;
    #else
    PS::F64 mass_dust = (r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*3.0*PLANET.mass*0.5;
    #endif
    
    time_sys = 0.0;// MOVE TO GLOBAL VALUABE
    PS::F64 dens = mass_dust * n_glb / (PI*(ax_out*ax_out - ax_in*ax_in));

    wtime_main_0 = wtime();
    
    #ifdef PHI_R_TREE
    GetPosDomainCyl(delta_ax, pos_domain, nx, ny);
    SetParticleKeplerDiskCyl2(system_grav, n_glb, ax_in, ax_out, dens, pos_domain, r_phy, true);
    #else //PHI_R_TREE
    GetPosDomain3(delta_ax,pos_domain);
    SetParticleKeplerDisk3(system_grav, n_glb, ax_in, ax_out, dens, pos_domain);
    #endif //PHI_R_TREE
    n_loc = system_grav.getNumberOfParticleLocal();
    n_glb = system_grav.getNumberOfParticleGlobal();
    SetPtclId(&system_grav[0], n_loc);
    Epj::r_search = r_hill * 3.0;
    Epj::r_coll   = r_phy;
    Epi::eps      = 0.0;
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cerr << "etime(IC) = " << wtime() - wtime_main_0 << std::endl;
    //////////
    // set offset of tree center in z direction
    PS::TREE_Z_COORD_OFFSET = 0.0;
    #ifdef SHIFT_CENTER_Z
    //CalcZCoordShift(&system_grav[0], n_loc, PS::TREE_Z_COORD_OFFSET);
    PS::TREE_Z_COORD_OFFSET = CalcZCoordShift(&system_grav[0], n_loc, r_phy*0.1);
    #endif
    PS::Comm::barrier();
    if(PS::Comm::getRank()==0){
        std::cerr<<"PS::TREE_Z_COORD_OFFSET= "<<PS::TREE_Z_COORD_OFFSET
                 <<" r_phy= "<<r_phy
                 <<std::endl;;
    }
    PS::Comm::barrier();
    // set offset of tree center in z direction
    //////////
    
    ChangePtclToSatellite(system_grav, SATELLITE, N_SATELLITE);
    SATELLITE_FORCE.resizeNoInitialize(64*N_SATELLITE); // We need large memory to use Long's force kenrel
    SATELLITE_FORCE_LOC.resizeNoInitialize(64*N_SATELLITE); // We need large memory to use Long's force kenrel
    for(PS::S32 ip=0; ip<SATELLITE_FORCE.size(); ip++) SATELLITE_FORCE[ip].clear();
    for(PS::S32 ip=0; ip<SATELLITE_FORCE_LOC.size(); ip++) SATELLITE_FORCE_LOC[ip].clear();
    n_loc = system_grav.getNumberOfParticleLocal();
    n_glb = system_grav.getNumberOfParticleGlobal();
    MASS_DUST_TOTAL = mass_dust*n_glb;
    PS::F64vec cm_pos_dust_loc = 0.0;
    for(PS::S32 i=0; i<n_loc; i++){
        cm_pos_dust_loc += system_grav[i].mass*system_grav[i].pos;
    }
    CM_POS_DUST_TOTAL.x = PS::Comm::getSum(cm_pos_dust_loc.x);
    CM_POS_DUST_TOTAL.y = PS::Comm::getSum(cm_pos_dust_loc.y);
    CM_POS_DUST_TOTAL.z = PS::Comm::getSum(cm_pos_dust_loc.z);
    CM_POS_DUST_TOTAL /= MASS_DUST_TOTAL;
    if(PS::Comm::getRank() == 0){
        std::cerr<<"after make satellite"<<std::endl;
        std::cerr<<"dens= "<<dens<<" total mass= "<<mass_dust*n_glb<<" mass_dust= "<<mass_dust
                 <<" CM_POS_DUST_TOTAL= "<<CM_POS_DUST_TOTAL
                 <<std::endl;
        if(N_SATELLITE > 0) std::cerr<<"SATELLITE[0].r_coll= "<<SATELLITE[0].r_coll<<std::endl;
        std::cerr<<"r_hill= "<<r_hill<<" mass_dust_total= "<<mass_dust * n_glb<<std::endl;
        std::cerr<<"Epj::r_search= "<<Epj::r_search<<std::endl;
    }
    NAN_TEST(system_grav,"[nan-0]");
    
#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cerr << "OK1 @main"  << std::endl;
#endif
    system_grav.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave);
#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cerr << "OK2 @main"  << std::endl;
#endif
    
    const PS::F32 coef_ema = 0.3;
    //const PS::F32 coef_ema = 1.0;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.delta_ax_ = delta_ax;
#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cerr << "OK3 @main"  << std::endl;
#endif
    dinfo.setNumberOfDomainMultiDimension(nx, ny, nz);
    if(PS::Comm::getRank() == 0){
        std::cerr<<"dinfo.getNDomain(0)= "<<dinfo.getNDomain(0)
                 <<" dinfo.getNDomain(1)= "<<dinfo.getNDomain(1)
                 <<" dinfo.getNDomain(2)= "<<dinfo.getNDomain(2)
                 <<std::endl;
    }
    n_loc = system_grav.getNumberOfParticleLocal();
    
#ifdef PHI_R_TREE
    ///////////////////////////////
    ///////// set phi and r in FP
    wtime_main_0 = wtime();
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_X);
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cerr << "etime(boundary) = " << wtime() - wtime_main_0 << std::endl;
    wtime_main_0 = wtime();
    dinfo.setPosRootDomain(PS::F64vec(0.0, -PI, -PI), PS::F64vec(2.0*PI, PI, PI)); // phi, r, z
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cerr << "etime(rootdomain) = " << wtime() - wtime_main_0 << std::endl;
    
    wtime_main_0 = wtime();
#ifdef STATIC_DD
    dinfo.setPosDomainCyl();
#else
    dinfo.collectSampleParticle(system_grav, true, 1.0);
    dinfo.decomposeDomainMultiStep3(true, false); // 1st: sample sort, 2nd split y-coord
#endif
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cerr << "etime(DD) = " << wtime() - wtime_main_0 << std::endl;
    wtime_main_0 = wtime();
    
#ifdef USE_SUPER_DOMAIN
    MY_MPI_COMM_SUB = dinfo.getCommSub(0);
    MY_MPI_COMM_1D  = dinfo.getComm1d(0); 
    MY_RANK_SUB     = dinfo.getRankSub(0);
#endif
    
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cerr << "etime(multistep2) = " << wtime() - wtime_main_0 << std::endl;
    PS::Comm::barrier();
    wtime_main_0 = wtime();
    
    n_loc = system_grav.getNumberOfParticleLocal();

    system_grav.exchangeParticle8v3(dinfo);
    system_grav.freeCommunicationBuffer();

    ANGULAR_MOMENTUM_TEST0(system_grav);
    PS::Comm::barrier();
    n_loc = system_grav.getNumberOfParticleLocal();
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cerr << "etime(exch_ptcl) = " << wtime() - wtime_main_0 << std::endl;

    /*
    PS::Comm::barrier();
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64 phi = atan2(system_grav[i].pos.y, system_grav[i].pos.x);
        if(phi < 0.0) phi += 2.0*4.0*atan(1.0);
        if(phi >= dinfo.getPosDomain(my_rank).high_.x || phi < dinfo.getPosDomain(my_rank).low_.x){
            std::cerr<<"rank="<<my_rank
                     <<" phi= "<<phi
                     <<" pos= "<<system_grav[i].pos
                     <<" getPosDomain(my_rank).high_.x= "<<dinfo.getPosDomain(my_rank).high_.x
                     <<" getPosDomain(my_rank).low_.x= "<<dinfo.getPosDomain(my_rank).low_.x
                     <<std::endl;
        }
        assert(phi < dinfo.getPosDomain(my_rank).high_.x && phi >= dinfo.getPosDomain(my_rank).low_.x);
    }
    */
    
    #ifdef USE_SUPER_DOMAIN
    wtime_main_0 = wtime();
    
    #ifdef STATIC_DD
    dinfo.setPosSuperDomainCyl(delta_ax);
    #else // STATIC_DD
    dinfo.initializeSuperDomain();
    #endif // STATIC_DD

    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cerr << "etime(superdomain) = " << wtime() - wtime_main_0 << std::endl;
    //exit(1);
    #endif // USE_SUPER_DOMAIN
#else //PHI_R_TREE
    dinfo.initializePosDomain(pos_domain);
#endif //PHI_R_TREER_PHI_TREE

#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cerr << "OK4 @main"  << std::endl;
#endif
    

    
#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cerr << "OK6 @main"  << std::endl;
#endif
    //checkLoadBalanceOfNumberOfParticle(system_grav);

    //=========================================================
    // Setup constant parameters used in Long's force kernel
    //=========================================================
#if defined(SUNWAY) && defined(SUNWAY_FORCE_KERNEL)
    //const PS::F64 Tdur = 5.0/32.0;
    //const PS::F64 Tdur = 0.0025*2.0*PI;
    const PS::F64 Tdur = DT_TREE*64;
    const PS::F64 ln_e_refl = std::log(0.5);
    //const PS::F64 ln_e_refl = std::log(0.1);
    // Tdur      := the oscillation period, which must be much less than
    //              the Kepler period (=2\pi in our unit) and be larger
    //              than DT_TREE to resolve oscillation motions.
    // ln_e_relf := the natural logarithmic of the reflection coefficient
    cpe_pars.dt          = DT_TREE;
    cpe_pars.r_ep        = Epj::r_coll;
    cpe_pars.r_sat       = Planet::r_coll;
#if 1
    cpe_pars.kappa       = std::pow((2.0*PI/Tdur),2.0); // k^{'}
    cpe_pars.eta         = 4.0*PI/(Tdur*std::sqrt(1.0+std::pow((PI/ln_e_refl),2.0))); // \eta^{'}
    //cpe_pars.kappa       = 1.0e-16; // k^{'}
    //cpe_pars.eta         = 1.0e-16; // \eta^{'}
#else
    // debug
    cpe_pars.kappa       = 1.0e-16; // k^{'}
    cpe_pars.eta         = 1.0e-16; // \eta^{'}
#endif
    cpe_pars.m_particle  = system_grav[0].mass; // the mass of a particle
    if (N_SATELLITE > 0) {
       cpe_pars.m_satelite = SATELLITE[0].mass; // the mass of a satellite
    } else {
       cpe_pars.m_satelite = 0.0; 
    }
    cpe_pars.m_planet    = PLANET.mass;
#endif

    //* Output the calculation parameters
    if (PS::Comm::getRank() == 0) {
        std::cerr << "DT_TREE      = " << DT_TREE << std::endl;
#if defined(SUNWAY) && defined(SUNWAY_FORCE_KERNEL)
        std::cerr << "Tdur         = " << Tdur << std::endl;
        std::cerr << "ln_e_refl    = " << ln_e_refl << std::endl;
#endif
        std::cerr << "r_ep         = " << cpe_pars.r_ep << std::endl;
        std::cerr << "r_sat        = " << cpe_pars.r_sat << std::endl;
        std::cerr << "k_dash       = " << cpe_pars.kappa << std::endl;
        std::cerr << "eta_dash     = " << cpe_pars.eta << std::endl;
        std::cerr << "m_planet     = " << cpe_pars.m_planet << std::endl;
        std::cerr << "n_loop_merge = " << n_loop_merge << std::endl;
        std::cerr << "n_smp_ave    = " << n_smp_ave << std::endl;
        std::cerr << "theta        = " << theta << std::endl;
    }


    ///////////////
    /// chose TREE TYPE
    typedef Epi EPI;
    typedef Epj EPJ;
    typedef Force FORCE;
    typedef SPJMonopoleScatterSW SPJ;
    TreeType tree_grav;
    tree_grav.initialize(n_glb, theta, n_leaf_limit, n_group_limit);

    n_loc = system_grav.getNumberOfParticleLocal();
    n_glb = system_grav.getNumberOfParticleGlobal();

    PS::ReallocatableArray<Satellite> SATELLITE_prev;
    PS::ReallocatableArray<SatelliteForce> SATELLITE_FORCE_prev;
    PS::ReallocatableArray<SatelliteForce> SATELLITE_FORCE_LOC_prev;

    SATELLITE_prev.resizeNoInitialize(N_SATELLITE);
    SATELLITE_FORCE_prev.resizeNoInitialize(64*N_SATELLITE);
    SATELLITE_FORCE_LOC_prev.resizeNoInitialize(64*N_SATELLITE);
    for(PS::S32 ip=0; ip<SATELLITE_FORCE_prev.size(); ip++) SATELLITE_FORCE_prev[ip].clear();
    for(PS::S32 ip=0; ip<SATELLITE_FORCE_LOC_prev.size(); ip++) SATELLITE_FORCE_LOC_prev[ip].clear();
    
    PS::ParticleSystem<FPGrav> system_grav_prev;
    system_grav_prev.initialize();
    system_grav_prev.setNumberOfParticleLocal(n_loc);

    PS::ReallocatableArray<Force> force_dust_from_satellite;
    
    /*
    PS::ReallocatableArray<PS::F64vec> * ang;
    ang = new PS::ReallocatableArray<PS::F64vec>[2];
    ang[0].resizeNoInitialize(n_loc);
    ang[1].resizeNoInitialize(n_loc);
    for(PS::S32 i=0; i<n_loc; i++){
        ang[0][i] = 0.0;
        ang[1][i] = 0.0;
    }
    */
    PS::F64vec ang[2];
    ang[0] = 0.0;
    ang[1] = 0.0;
    for(PS::S32 loop=0; loop<2; loop++){
        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"CHECK A"<<std::endl;
        PS::Comm::barrier();
        for(PS::S32 i=0; i<N_SATELLITE; i++){
            SATELLITE_prev[i] = SATELLITE[i];
        }

        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"CHECK B"<<std::endl;
        PS::Comm::barrier();
        
        for(PS::S32 i=0; i<64*N_SATELLITE; i++){
            SATELLITE_FORCE_prev[i] = SATELLITE_FORCE[i];
            SATELLITE_FORCE_LOC_prev[i] = SATELLITE_FORCE_LOC[i];
        }
        
        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"CHECK C"<<std::endl;
        PS::Comm::barrier();
        
        for(PS::S32 i=0; i<n_loc; i++){
            system_grav_prev[i] = system_grav[i];
        }

        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"CHECK D"<<std::endl;
        PS::Comm::barrier();
        
        TREE_POINTER = & tree_grav;
        tree_grav.epi_sorted_.copyPointer(&system_grav[0]);
        tree_grav.epi_sorted_.setSize(system_grav.getNumberOfParticleLocal());
        tree_grav.epi_sorted_.setCapacity(system_grav.getNumberOfParticleLocal());
        tree_grav.epj_sorted_loc_.copyPointer(&system_grav[0]);
        tree_grav.epj_sorted_loc_.setSize(system_grav.getNumberOfParticleLocal());
        tree_grav.epj_sorted_loc_.setCapacity(system_grav.getNumberOfParticleLocal());

        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"CHECK E"<<std::endl;
        PS::Comm::barrier();
        
        CalcForceLoopMerge(tree_grav,
                           force_dust_from_satellite,
                           DispatchKernelWithSP<EPI, EPJ, SPJ>,
                           RetrieveKernel<FORCE>,
                           system_grav,
                           SATELLITE,
                           PLANET,
                           dinfo,
                           N_WALK_LIMIT,
                           false,
                           true);

        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"CHECK F"<<std::endl;
        PS::Comm::barrier();
        
        for(PS::S32 i=0; i<n_loc; i++){
            ang[loop] += system_grav[i].mass*(system_grav[i].pos^system_grav[i].vel);
        }
        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"ang[loop]= "<<ang[loop]<<std::endl;
        PS::Comm::barrier();
        
        for(PS::S32 i=0; i<n_loc; i++){
            system_grav[i] = system_grav_prev[i];
        }
        for(PS::S32 i=0; i<N_SATELLITE; i++){
            SATELLITE[i] = SATELLITE_prev[i];
        }
        for(PS::S32 i=0; i<64*N_SATELLITE; i++){
            SATELLITE_FORCE[i] = SATELLITE_FORCE_prev[i];
            SATELLITE_FORCE_LOC[i] = SATELLITE_FORCE_LOC_prev[i];
        }
    }
    PS::Comm::barrier();
    if(PS::Comm::getRank()==0) std::cerr<<"end of the loop "<<std::endl;
    PS::Comm::barrier();
    
    PS::S32 err_loc = 0;
    if(ang[0].x != ang[1].x) err_loc = 1;
    if(ang[0].y != ang[1].y) err_loc = 1;
    if(ang[0].z != ang[1].z) err_loc = 1;

    int err_glb = PS::Comm::getSum(err_loc);
    if(err_glb == 0){
        if(PS::Comm::getRank()==0){
            std::cerr<<"PASS"<<std::endl;
        }
        athread_halt();
        PS::Finalize();
        return 0;
    }
    else{
        if(PS::Comm::getRank()==0){
            std::cerr<<"FAIL"<<std::endl;
        }
    }
    
    PS::TempArray<PS::S32> err_array;
    err_array.free();
    err_array.resizeNoInitialize(n_proc);
    MPI_Allgather(&err_loc, 1, MPI_INT, err_array.getPointer(), 1, MPI_INT, MPI_COMM_WORLD);

    for (PS::S32 irank=0; irank<n_proc; irank++) {
        if(PS::Comm::getRank()==0){        
            if(irank % 100==0){
                if(PS::Comm::getRank()==0) std::cerr<<"irank= "<<irank<<std::endl;
            }
        }
        if(err_array[irank] == 1){
            if(irank==PS::Comm::getRank()){
                std::cerr << "RANK= " << irank
                          << " HOST= " << HOST_NAME[irank] << std::endl;
                std::cerr<<"ang[0]= "<<ang[0]
                         <<" ang[1]= "<<ang[1]
                         <<std::endl;
            }
            PS::Comm::barrier();
        }
    }
    if (PS::Comm::getRank() == 0) fout_wtime.close();
    
#ifdef SUNWAY
    athread_halt();
#endif
    
    PS::Finalize();

    return 0;
}

