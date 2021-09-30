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

//PS::S32 N_SATELLITE = 0;
//PS::S32 N_SATELLITE = 1000;
//PS::S32 N_SATELLITE = 64;
PS::S32 N_SATELLITE = 128;

PS::ReallocatableArray<Satellite> SATELLITE;
PS::ReallocatableArray<SatelliteForce> SATELLITE_FORCE;
PS::ReallocatableArray<SatelliteForce> SATELLITE_FORCE_LOC;

std::ofstream fout_debug;
std::ofstream fout_wtime;

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
    size_t size_memrory_pool = 1500000000;
    PS::MemoryPool::initialize(size_memrory_pool);
#endif
    
    LapTimer::initialize();

    if (PS::Comm::getRank() == 0)
       std::cout << "etime(MPI startup) = " << wtime() - wtime_start_of_main << std::endl;

#ifdef SUNWAY
    wtime_main_0 = wtime();
    athread_init();
    if (PS::Comm::getRank() == 0)
       std::cout << "etime(athread_init) = " << wtime() - wtime_main_0 << std::endl;
#endif
    
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

    PS::S32 snp_id = 0;
    const PS::F64 PI = 4.0 * atan(1.0);
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::S32 my_rank = PS::Comm::getRank();

    MY_RANK_MPI = my_rank;
    
    //const PS::S32 n_loop_merge = 16;
    //const PS::S32 n_loop_merge = 32;
    const PS::S32 n_loop_merge = 64;
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
    PS::F32 time_end = DT_TREE*128;
    //PS::F32 time_end = DT_TREE*320; // performance measurement
    //PS::F32 time_end = DT_TREE*512;
    //PS::F32 time_end = DT_TREE*1024;
    //PS::F32 time_end = DT_TREE*2048; // 2.5 Kepler periods
    //PS::F32 time_end = DT_TREE*4096;
    //PS::F32 time_end = DT_TREE*8192; // 10 Kepler periods
    //PS::F32 time_end = DT_TREE*32768; // 40 Kepler periods
    //PS::F32 time_end = 1.0; // for performance measurement
    //PS::F32 time_end = 1.5;
    //PS::F32 time_end = 2.0;
    PS::F32 dt_diag = 1.0 / 8.0;
    PS::F32 dt_snap = 1.0;
    char dir_name[1024];
    PS::S64 n_glb = 1024;
    PS::S64 n_loc = 0;
    PS::F64 ax_cen = 1.0;
    //PS::F64 delta_ax = 0.1;
    //PS::F64 delta_ax = 0.01;
    PS::F64 delta_ax = 1e-3;
    PS::F64 ax_in  = ax_cen - 0.5*delta_ax;
    PS::F64 ax_out = ax_cen + 0.5*delta_ax;
    PS::F64 ecc_rms = 0.0;
    PS::F64 inc_rms = 16.0;
    PS::F64 tau = 1.0;
    PS::F64 ratio_r_phy_r_hill = 1.0;
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

    //* Open fout_debug
    wtime_main_0 = wtime();
    std::string s_debug = dir_name + ("/debug_" + ToS(PS::Comm::getRank()) + ".dat");
    if (PS::Comm::getRank() == 0) {
        std::cerr<<"s_debug: "<<s_debug<<std::endl;
        fout_debug.open(s_debug.c_str());
        if(!fout_debug.is_open()) exit(1);
        std::cout << "etime(fopen) = " << wtime() - wtime_main_0 << std::endl;
    }
    std::string s_wtime = dir_name + ("/wtime_" + ToS(PS::Comm::getRank()) + ".dat");
    if (PS::Comm::getRank() == 0) {
        std::cerr<<"s_wtime: "<<s_wtime<<std::endl;
        fout_wtime.open(s_wtime.c_str());
        if(!fout_wtime.is_open()) exit(1);
        std::cout << "etime(fopen) = " << wtime() - wtime_main_0 << std::endl;
    }
    

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

    //===================================
    //   Create ParticleSystem object
    //===================================
    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();

    //=================================
    //   Determine the process grid
    //=================================
#ifdef PHI_R_TREE
    PS::S32 nx, ny, nz;
    DivideNProc(nx, ny, nz, n_proc, delta_ax);
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

    //================================================== 
    //  Set parameters for planetary ring simulation 
    //================================================== 
    PS::F64 r_phy = sqrt(tau*(ax_out*ax_out - ax_in*ax_in) / n_glb);
    PS::F64 r_hill = r_phy / ratio_r_phy_r_hill;
    PS::F64 mass_dust = (r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*3.0*PLANET.mass*0.5;
    
    time_sys = 0.0;// MOVE TO GLOBAL VALUABE
    PS::F64 dens = mass_dust * n_glb / (PI*(ax_out*ax_out - ax_in*ax_in));

    wtime_main_0 = wtime();
    
    //================================================== 
    //  Make the initial condition
    //================================================== 
#ifdef PHI_R_TREE
    GetPosDomainCyl(delta_ax, pos_domain, nx, ny);
    SetParticleKeplerDiskCyl2(system_grav, n_glb, ax_in, ax_out, dens, pos_domain, true);
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
    if (PS::Comm::getRank() == 0) std::cout << "etime(IC) = " << wtime() - wtime_main_0 << std::endl;

    //* Set offset of tree center in z direction
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
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "OK1 @main"  << std::endl;
#endif
    system_grav.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave);
#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "OK2 @main"  << std::endl;
#endif
    
    //================================
    //   Create DomainInfo object
    //================================
    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.delta_ax_ = delta_ax;
#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "OK3 @main"  << std::endl;
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
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(boundary) = " << wtime() - wtime_main_0 << std::endl;
    wtime_main_0 = wtime();
    dinfo.setPosRootDomain(PS::F64vec(0.0, -PI, -PI), PS::F64vec(2.0*PI, PI, PI)); // phi, r, z
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(rootdomain) = " << wtime() - wtime_main_0 << std::endl;
    
    wtime_main_0 = wtime();
#ifdef STATIC_DD
    dinfo.setPosDomainCyl();
#else
    dinfo.collectSampleParticle(system_grav, true, 1.0);
    dinfo.decomposeDomainMultiStep3(true, false); // 1st: sample sort, 2nd split y-coord
#endif
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(DD) = " << wtime() - wtime_main_0 << std::endl;
    wtime_main_0 = wtime();
    
#ifdef USE_SUPER_DOMAIN
    MY_MPI_COMM_SUB = dinfo.getCommSub(0);
    MY_MPI_COMM_1D  = dinfo.getComm1d(0); 
    MY_RANK_SUB     = dinfo.getRankSub(0);
#endif
    
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(multistep2) = " << wtime() - wtime_main_0 << std::endl;
    PS::Comm::barrier();
    wtime_main_0 = wtime();
    

    //==================================
    //   Exchange Particle
    //==================================
    //system_grav.exchangeParticleSw(dinfo, false);
    //system_grav.exchangeParticle8(dinfo);
    //system_grav.exchangeParticle8v2(dinfo);
    system_grav.exchangeParticle8v3(dinfo);
    system_grav.freeCommunicationBuffer();

    //==================================
    //   Check particle positions
    //==================================
    n_loc = system_grav.getNumberOfParticleLocal();
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) {
        std::cout << "Check if particles are in my domain...[after exchangeParticle()]" << std::endl;
    }
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64 phi = atan2(system_grav[i].pos.y, system_grav[i].pos.x);
        if(phi < 0.0) phi += 2.0*4.0*atan(1.0);
        if(phi >= dinfo.getPosDomain(my_rank).high_.x || phi < dinfo.getPosDomain(my_rank).low_.x){
            std::cerr << "rank = "<< my_rank << " " 
                      << "i = " << i << " "
                      << "n_loc = " << n_loc << " "
                      << "phi = " << phi << " "
                      << "id = " << system_grav[i].id << " "
                      << "pos = " << system_grav[i].pos << std::endl;
        }
        assert(phi < dinfo.getPosDomain(my_rank).high_.x && phi >= dinfo.getPosDomain(my_rank).low_.x);
    }
    //athread_halt();
    //PS::Finalize();
    //std::exit(0);


    ANGULAR_MOMENTUM_TEST0(system_grav);
    PS::Comm::barrier();
#ifdef DUMP_TARGET_PTCL
    PS::DumpTargetPtcl(system_grav,
                       PS::CHECK_ID,
                       PS::N_CHECK_ID);
#endif
#ifdef DUMP_MEMORY_SIZE
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0){
        std::cerr<<"after ex ptcls (1st)"<<std::endl;
    }
    system_grav.dumpMemSizeUsed(std::cerr);
    PS::Comm::barrier();
#endif
    
    n_loc = system_grav.getNumberOfParticleLocal();
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(exch_ptcl) = " << wtime() - wtime_main_0 << std::endl;
    
#ifdef USE_SUPER_DOMAIN
    wtime_main_0 = wtime();
#ifdef STATIC_DD
    dinfo.setPosSuperDomainCyl(delta_ax);
#else // STATIC_DD
    dinfo.initializeSuperDomain();
#endif // STATIC_DD
#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0){
        for(PS::S32 i=0; i<dinfo.getNDomain(0); i++){
            fout_debug<<"i= "
                      <<" dinfo.getPosSuperDomain(i)= "
                      <<dinfo.getPosSuperDomain(i)
                      <<std::endl;
        }
    }
    PS::Comm::barrier();
#endif // DEBUG_PRINT_MAIN
    
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(superdomain) = " << wtime() - wtime_main_0 << std::endl;
    //exit(1);
#endif // USE_SUPER_DOMAIN
#else //PHI_R_TREE
    dinfo.initializePosDomain(pos_domain);
#endif //PHI_R_TREER_PHI_TREE


    
#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "OK6 @main"  << std::endl;
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
    cpe_pars.kappa       = std::pow((2.0*PI/Tdur),2.0); // k^{'}
    cpe_pars.eta         = 4.0*PI/(Tdur*std::sqrt(1.0+std::pow((PI/ln_e_refl),2.0))); // \eta^{'}
    //cpe_pars.kappa       = 1.0e-16; // k^{'}
    //cpe_pars.eta         = 1.0e-16; // \eta^{'}
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
        std::cout << "DT_TREE      = " << DT_TREE << std::endl;
#if defined(SUNWAY) && defined(SUNWAY_FORCE_KERNEL)
        std::cout << "Tdur         = " << Tdur << std::endl;
        std::cout << "ln_e_refl    = " << ln_e_refl << std::endl;
#endif
        std::cout << "r_ep         = " << cpe_pars.r_ep << std::endl;
        std::cout << "r_sat        = " << cpe_pars.r_sat << std::endl;
        std::cout << "k_dash       = " << cpe_pars.kappa << std::endl;
        std::cout << "eta_dash     = " << cpe_pars.eta << std::endl;
        std::cout << "m_planet     = " << cpe_pars.m_planet << std::endl;
        std::cout << "n_loop_merge = " << n_loop_merge << std::endl;
        std::cout << "n_smp_ave    = " << n_smp_ave << std::endl;
        std::cout << "theta        = " << theta << std::endl;
    }


    //=============================================
    //   Create and setup TreeForForce object
    //=============================================
    typedef Epi EPI;
    typedef Epj EPJ;
    typedef Force FORCE;
    typedef SPJMonopoleScatterSW SPJ;
    TreeType tree_grav;
    TREE_POINTER = & tree_grav;
    tree_grav.epi_sorted_.copyPointer(&system_grav[0]);
    tree_grav.epi_sorted_.setSize(system_grav.getNumberOfParticleLocal());
    tree_grav.epi_sorted_.setCapacity(system_grav.getNumberOfParticleLocal());
#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier();if(PS::Comm::getRank()==0)std::cerr<<"OK7 @main"<<std::endl;
#endif
    PS::ReallocatableArray<Force> force_dust_from_satellite;
    tree_grav.initialize(n_glb, theta, n_leaf_limit, n_group_limit);
    PS::Comm::barrier();
    
#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier();if(PS::Comm::getRank()==0)std::cerr<<"OK8 @main"<<std::endl;
#endif

#ifdef WRITE_SNP
    if(my_rank==0) WriteSnp(system_grav, snp_id, dir_name, 10000);
    PS::Comm::barrier();
#endif

#ifdef FORCE_SUM_CHECK_DIRECT
    PS::F64vec force_sum_direct_glb = CalcForceSumDirect(system_grav);
    PS::Comm::barrier();
    if(PS::Comm::getRank()==0){
        std::cerr<<"force_sum_direct_glb= "
                 <<force_sum_direct_glb
                 <<std::endl;
    }
    PS::Comm::barrier();
#endif

#ifdef FORCE_CHECK
    n_loc = system_grav.getNumberOfParticleLocal();
    SatelliteForce * fi_di_org = new SatelliteForce[n_loc+N_SATELLITE]; // NOTE: not only for satellite force
    CalcForceDirect(system_grav, SATELLITE.getPointer(), N_SATELLITE, Epj::r_coll, cpe_pars.kappa, fi_di_org);

    std::vector< std::pair<PS::S64, PS::S32> > id_adr(n_loc);
    for(PS::S32 i=0; i<n_loc; i++){
        id_adr[i] = std::make_pair(system_grav[i].id, i);
    }
    std::sort(id_adr.begin(), id_adr.end());
    SatelliteForce * fi_di_srt = new SatelliteForce[n_loc];
    for(PS::S32 i=0; i<n_loc; i++){
        const PS::S32 adr = id_adr[i].second;
        fi_di_srt[i] = fi_di_org[adr];
    }

    cpe_pars.energy_flag = 0;
    cpe_pars.first_flag = 0; // Full-kick
    cpe_pars.last_flag = 0;  // Full-drift
    tree_grav.calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe
        (DispatchKernelWithSP<EPI, EPJ, SPJ>,
         RetrieveKernel<FORCE>, 0,
         system_grav, dinfo, N_WALK_LIMIT,
         false,  false);
    n_loc = system_grav.getNumberOfParticleLocal();
    
#ifndef REMOVE_EPI_SORTED
    for(PS::S32 i=0; i<n_loc; i++){
        system_grav[i].pos = (((EPI*)(EPI_POINTER))+i)->pos;
        system_grav[i].vel = (((EPI*)(EPI_POINTER))+i)->vel;
        system_grav[i].mass = (((EPI*)(EPI_POINTER))+i)->mass; // by M.I.
        system_grav[i].id = (((EPI*)(EPI_POINTER))+i)->id; // by M.I.
    }
#endif
    
    for(PS::S32 i=0; i<n_loc; i++){
        id_adr[i] = std::make_pair(system_grav[i].id, i);
    }
    std::sort(id_adr.begin(), id_adr.end());
    SatelliteForce * fi_tr_srt = new SatelliteForce[n_loc];
    for(PS::S32 i=0; i<n_loc; i++){
        const PS::S32 adr = id_adr[i].second;
        fi_tr_srt[i].acc = system_grav[adr].vel;
    }
    
    if(PS::Comm::getRank()==0){
        for(PS::S32 i=0; i<n_loc; i++){
            std::cerr<<"i= "
                     <<" id= "<<id_adr[i].first
                     <<" fi_di= "<<fi_di_srt[i].acc
                     <<" fi_tr= "<<fi_tr_srt[i].acc
                     <<" err= "<<(fi_tr_srt[i].acc.x - fi_di_srt[i].acc.x) / fabs(fi_di_srt[i].acc.x)
                     <<" "<<(fi_tr_srt[i].acc.y - fi_di_srt[i].acc.y) / fabs(fi_di_srt[i].acc.y)
                     <<" "<<(fi_tr_srt[i].acc.z - fi_di_srt[i].acc.z) / fabs(fi_di_srt[i].acc.z)
                     <<std::endl;
        }
    }
    
    PS::Finalize();
    return 0;
#endif


#ifdef REMOVE_EPI_SORTED
        tree_grav.epi_sorted_.copyPointer(&system_grav[0]);
        tree_grav.epi_sorted_.setSize(system_grav.getNumberOfParticleLocal());
        tree_grav.epi_sorted_.setCapacity(system_grav.getNumberOfParticleLocal());
#endif
#ifdef REMOVE_EPJ_SORTED_LOC
        tree_grav.epj_sorted_loc_.copyPointer(&system_grav[0]);
        tree_grav.epj_sorted_loc_.setSize(system_grav.getNumberOfParticleLocal());
        tree_grav.epj_sorted_loc_.setCapacity(system_grav.getNumberOfParticleLocal());
#endif

        /*
    CalcForce(tree_grav,
              force_dust_from_satellite,
              DispatchKernelWithSP<EPI, EPJ, SPJ>,
              RetrieveKernel<FORCE>,
              system_grav,
              SATELLITE,
              PLANET,
              dinfo,
              N_WALK_LIMIT,
              false);
        */

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

#ifdef DUMP_MEMORY_SIZE
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0){
        std::cerr<<"after force calc (1st)"<<std::endl;
    }
    tree_grav.dumpMemSizeUsed(std::cerr);
    //tree_grav.comm_table_.dumpMemSizeUsed(std::cerr);
    PS::Comm::barrier();
#endif

    //exit(1);
    
    n_loc = system_grav.getNumberOfParticleLocal();
    n_glb = system_grav.getNumberOfParticleGlobal();
    nan_test((EPI*)EPI_POINTER, system_grav.getNumberOfParticleLocal(), "nan-0");
     
#ifdef DEBUG_PRINT_MAIN
     PS::Comm::barrier(); if(PS::Comm::getRank()==0)std::cerr<<"OK9 @main"<<std::endl;
     CalcCM((EPI*)EPI_POINTER, n_loc, MASS_DUST_TOTAL, CM_POS_DUST_TOTAL);
     if(PS::Comm::getRank()==0){
         std::cerr<<" MASS_DUST_TOTAL= "<<MASS_DUST_TOTAL
                  <<" CM_POS_DUST_TOTAL= "<<CM_POS_DUST_TOTAL
                  <<std::endl;         
     }
#endif
     //exit(1);

#ifndef REMOVE_EPI_SORTED
     for(PS::S32 i=0; i<n_loc; i++){
         system_grav[i].pos = (((EPI*)(EPI_POINTER))+i)->pos;
         system_grav[i].vel = (((EPI*)(EPI_POINTER))+i)->vel;
         system_grav[i].mass = (((EPI*)(EPI_POINTER))+i)->mass; // by M.I.
         system_grav[i].id = (((EPI*)(EPI_POINTER))+i)->id; // by M.I.
     }
#endif
     
#ifdef CHECK_AT_AFTER_FIRST_LOOP
    PS::Comm::barrier();
    if(PS::Comm::getRank()==0){
        std::cerr<<"end of first loop"<<std::endl;
        std::cerr<<"n_loc= "<<n_loc<<std::endl;
        for(PS::S32 i=0; i<n_loc; i+=n_loc/1000){
            std::cerr<<"i= "<<i
                     <<" sys[i].id= "<<system_grav[i].id
                     <<" mass= "<<system_grav[i].mass
                     <<" pos= "<<system_grav[i].pos
                     <<" vel= "<<system_grav[i].vel
                     <<std::endl;
        }
    }
    exit(1);
#endif
     PS::Comm::barrier();if(PS::Comm::getRank()==0)std::cerr<<"OK10 @main"<<std::endl;
     
    //system_grav.dumpMemSizeUsed(fout_wtime);
    //tree_grav.dumpMemSizeUsed(fout_wtime);

#if defined(ROTATE_PTCLS) && defined(SHIFT_CENTER_Z)
    starttime = PS::GetWtime();
    RotateAndShiftPtcl(&system_grav[0],
                       system_grav.getNumberOfParticleLocal(),
                       -DT_TREE,
                       0.1*r_phy,
                       PS::TREE_Z_COORD_OFFSET);
    endtime = PS::GetWtime();
    if (PS::Comm::getRank() == 0)
        fout_wtime << "wtime_calc_shift_z = " << endtime - starttime << std::endl;
    #ifdef DEBUG_SHIFT_CENTER_Z
    if (PS::Comm::getRank() == 0)
        std::cerr<<"A) PS::TREE_Z_COORD_OFFSET(0)= "<<PS::TREE_Z_COORD_OFFSET<<std::endl;
    PS::TREE_Z_COORD_OFFSET = 0.0;
    PS::TREE_Z_COORD_OFFSET = CalcZCoordShift(&system_grav[0], system_grav.getNumberOfParticleLocal(), r_phy*0.1);
    if (PS::Comm::getRank() == 0)
        std::cerr<<"A) PS::TREE_Z_COORD_OFFSET(1)= "<<PS::TREE_Z_COORD_OFFSET<<std::endl;
    #endif
#else
    #ifdef ROTATE_PTCLS
     //* Rotate the coordinate of particles and satellites
     #ifdef SUNWAY
     starttime = PS::GetWtime();
     theta_rot = - DT_TREE;
     const PS::F64 cth = std::cos(theta_rot);
     const PS::F64 sth = std::sin(theta_rot);
     //unsigned long args[4];
     //** Rotate the particles
     args[0] = (unsigned long) system_grav.getNumberOfParticleLocal();
     args[1] = (unsigned long) &system_grav[0];
     args[2] = (unsigned long) &cth;
     args[3] = (unsigned long) &sth;
     __real_athread_spawn((void*)slave_Rotate, args);
     athread_join();
     //** Rotate the satellites
     args[0] = (unsigned long) N_SATELLITE;
     args[1] = (unsigned long) SATELLITE.getPointer();
     args[2] = (unsigned long) &cth;
     args[3] = (unsigned long) &sth;
     __real_athread_spawn((void*)slave_Rotate, args);
     athread_join();
     endtime = PS::GetWtime();
     if (PS::Comm::getRank() == 0)
         fout_wtime << "wtime_rotate = " << endtime - starttime << std::endl;
    #else // SUNWAY
     theta_rot = - DT_TREE;
     Rotate(system_grav,n_loc,theta_rot);
     Rotate(SATELLITE,N_SATELLITE,theta_rot);
     endtime = PS::GetWtime();
     if (PS::Comm::getRank() == 0)
         fout_wtime << "wtime_rotate = " << endtime - starttime << std::endl;
    #endif //SUNWAY
    #endif
    #ifdef SHIFT_CENTER_Z
     starttime = PS::GetWtime();
     z_coord_offset_loc = 0.0;
     for(PS::S32 i=0; i<n_loc; i++){
         if( z_coord_offset_loc < fabs(system_grav[i].pos.z) ) z_coord_offset_loc = fabs(system_grav[i].pos.z);
     }
     PS::TREE_Z_COORD_OFFSET = PS::Comm::getMaxValue(z_coord_offset_loc) + 1.0e-6;
     endtime = PS::GetWtime();
     if (PS::Comm::getRank() == 0)
         fout_wtime << "wtime_calc_shift_z = " << endtime - starttime << std::endl;
    #endif
#endif

     
     //* Update time_sys
     time_sys += DT_TREE;

#ifdef WRITE_SNP
        if(my_rank==0) WriteSnp(system_grav, snp_id, dir_name, 10000);
        PS::Comm::barrier();
        //PS::Finalize();
        //return 0;
#endif
     

    //* Compute the total anugular momentum of the system at t=dt.
    ANGULAR_MOMENTUM_TEST0(system_grav);
    PS::Comm::barrier();

    
    //* Compute energy    
    PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
#if !defined(SUNWAY) || !defined(SUNWAY_FORCE_KERNEL)
    //calcEnergy(system_grav, SATELLITE, SATELLITE_FORCE, Etot0, Ekin0, Epot0);
    CalcEnergy(&system_grav[0], FORCE_ARRAY,
               system_grav.getNumberOfParticleLocal(),
               SATELLITE.getPointer(), SATELLITE_FORCE.getPointer(), N_SATELLITE,
               Etot0, Ekin0, Epot0);
    if(PS::Comm::getRank() == 0){
        std::cerr<<"Ekin0= "<<Ekin0<<" Epot0= "<<Epot0<<" Etot0= "<<Etot0<<std::endl;
    }

    kick(system_grav, DT_TREE*0.5);
    kick(SATELLITE.getPointer(), SATELLITE_FORCE.getPointer(), SATELLITE.size(), DT_TREE*0.5);
    drift(system_grav, DT_TREE);
    drift(SATELLITE.getPointer(), SATELLITE.size(), DT_TREE);
#endif

#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier();if(PS::Comm::getRank()==0)std::cerr<<"OK11 @main"<<std::endl;
#endif
    PS::F64 time_diag = 0.0;
    PS::F64 time_snap = 0.0;
    PS::S64 n_loop = 0;
    PS::S32 id_snap = 0;
    PS::ReallocatableArray<FPGrav> ptcl_old;
    PS::ReallocatableArray<FPGrav> ptcl_new;
    //* Clear time profiles
    TREE_POINTER->clearTimeProfile();
    TREE_POINTER->clearCounterAll();
    system_grav.clearTimeProfile();
    system_grav.clearCounterAll();
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

    //=========================================================
    //   Main loop
    //=========================================================
    while(time_sys < time_end){
    //for(int loop=0; loop<1; loop++){
        PS::Comm::barrier();
        if (PS::Comm::getRank() == 0) {
             std::cout << "#######################" << std::endl;
             std::cout << "  n_loop = " << n_loop << std::endl;
             std::cout << "#######################" << std::endl;
             fout_wtime << "#######################" << std::endl;
             fout_wtime << "  n_loop = " << n_loop << std::endl;
             fout_wtime << "#######################" << std::endl;
             fout_debug << "#######################" << std::endl;
             fout_debug << "  n_loop = " << n_loop << std::endl;
             fout_debug << "#######################" << std::endl;
        }
        N_LOOP_GLB = n_loop;
        PS::Comm::barrier();
        PS::F64 wtime_offset = PS::GetWtime();

        PS::Comm::barrier();
        NAN_TEST(system_grav,"[nan-a]");

        //#####
        //dinfo.decomposeDomainAll(system_grav); // test
        //------
        //dinfo.collectSampleParticle(system_grav, true, 1.0);
        //dinfo.decomposeDomain2();
        //#####
        
        PS::Comm::barrier();
        NAN_TEST(system_grav,"[nan-b]");
#ifdef DEBUG_PRINT_MAIN
        //PS::Comm::barrier();
        //if(PS::Comm::getRank()==0)std::cerr<<"OK12-"<<n_loop<<" @main"<<std::endl;
#endif
        n_loc = system_grav.getNumberOfParticleLocal();

#ifdef DEBUG_PRINT_MAIN
        PS::Comm::barrier();
        CalcCM((EPI*)EPI_POINTER, n_loc, MASS_DUST_TOTAL, CM_POS_DUST_TOTAL);
        if(PS::Comm::getRank()==0){
            std::cerr<<" (EPI) MASS_DUST_TOTAL= "<<MASS_DUST_TOTAL
                     <<" CM_POS_DUST_TOTAL= "<<CM_POS_DUST_TOTAL
                     <<std::endl;
        }
        CalcCMCyl(&system_grav[0], n_loc, MASS_DUST_TOTAL, CM_POS_DUST_TOTAL);
        if(PS::Comm::getRank()==0){
            std::cerr<<" (FP.1) MASS_DUST_TOTAL= "<<MASS_DUST_TOTAL
                     <<" CM_POS_DUST_TOTAL= "<<CM_POS_DUST_TOTAL
                     <<std::endl;
            std::cerr<<" if rotate ptcl is on EPI != FP.1"<<std::endl;
        }
#endif
        
        PS::Comm::barrier();
        PS::WtimeAbs::zero_point_ = MPI::Wtime();

#ifdef STATIC_DD
        //system_grav.exchangeParticleSw(dinfo, false);
        //system_grav.exchangeParticle8(dinfo);
        //system_grav.exchangeParticle8v2(dinfo);
        system_grav.exchangeParticle8v3(dinfo);
        system_grav.freeCommunicationBuffer();
#else // STATIC_DD
        if( (n_loop+1) % 50 ==0){
            dinfo.decomposeDomainMultiStep3(true, false); // 1st: sample sort, 2nd split y-coord
        }
    #ifdef FAST_EX_PTCL
        if(n_loop <= 3){
            //system_grav.exchangeParticle8v2(dinfo);
            system_grav.exchangeParticle8v3(dinfo);
        }
        else{
            //system_grav.exchangeParticle8v2(dinfo);
            system_grav.exchangeParticle8v3(dinfo);
        }
    #else // FAST_EX_PTCL
        system_grav.exchangeParticleSw(dinfo, false);
    #endif // FAST_EX_PTCL
#endif // STATIC_DD

#ifdef DUMP_TARGET_PTCL
    PS::DumpTargetPtcl(system_grav,
                       PS::CHECK_ID,
                       PS::N_CHECK_ID);
#endif
        
#ifdef DUMP_MEMORY_SIZE
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0){
        std::cerr<<"after ex ptcls"<<std::endl;
    }
    system_grav.dumpMemSizeUsed(std::cerr);
    PS::Comm::barrier();
#endif


    
#ifdef REMOVE_EPI_SORTED
        tree_grav.epi_sorted_.copyPointer(&system_grav[0]);
        tree_grav.epi_sorted_.setSize(system_grav.getNumberOfParticleLocal());
        tree_grav.epi_sorted_.setCapacity(system_grav.getNumberOfParticleLocal());
#endif
#ifdef REMOVE_EPJ_SORTED_LOC
        tree_grav.epj_sorted_loc_.copyPointer(&system_grav[0]);
        tree_grav.epj_sorted_loc_.setSize(system_grav.getNumberOfParticleLocal());
        tree_grav.epj_sorted_loc_.setCapacity(system_grav.getNumberOfParticleLocal());
#endif
        n_loc = system_grav.getNumberOfParticleLocal();

#ifdef DEBUG_PRINT_MAIN
        CalcCMCyl(&system_grav[0], n_loc, MASS_DUST_TOTAL, CM_POS_DUST_TOTAL);
        if(PS::Comm::getRank()==0){
            std::cerr<<" (FP.2) MASS_DUST_TOTAL= "<<MASS_DUST_TOTAL
                     <<" CM_POS_DUST_TOTAL= "<<CM_POS_DUST_TOTAL
                     <<std::endl;         
        }
#endif
        NAN_TEST(system_grav,"[nan-c]")

        
#ifdef DEBUG_PRINT_MAIN
        PS::Comm::barrier(); if(PS::Comm::getRank()==0)std::cerr<<"OK15-"<<n_loop<<" @main"<<std::endl;
#endif
        
        CalcForceLoopMerge(tree_grav,
                           force_dust_from_satellite,
                           DispatchKernelWithSPKickDrift<EPI, EPJ, SPJ>,
                           RetrieveKernel<FORCE>,
                           system_grav,
                           SATELLITE,
                           PLANET,
                           dinfo,
                           N_WALK_LIMIT,
                           n_loop_merge);

#ifdef DUMP_MEMORY_SIZE
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0){
        std::cerr<<"after force calc"<<std::endl;
    }
    tree_grav.dumpMemSizeUsed(std::cerr);
    //tree_grav.comm_table_.dumpMemSizeUsed(std::cerr);
    PS::Comm::barrier();
#endif
        
        //PS::Comm::barrier();
        //exit(1);

        //exit(1);
            
        //system_grav.dumpMemSizeUsed(fout_wtime);
        //tree_grav.dumpMemSizeUsed(fout_wtime);
        //PS::Comm::barrier();
        //fout_wtime.close();
        //athread_halt();
        //PS::Finalize();
        //return 0;

        NAN_TEST(system_grav,"[nan-d]");
        nan_test((EPI*)EPI_POINTER, system_grav.getNumberOfParticleLocal(), "nan-e");
        if (PS::Comm::getRank() == 0){
            for(PS::S32 i=0; i<5; i++){
                std::cerr<<"A) system_grav[i].pos= "<<system_grav[i].pos
                         <<" id= "<<system_grav[i].id
                         <<std::endl;
            }
        }
        starttime = PS::GetWtime();
        n_loc = system_grav.getNumberOfParticleLocal();

#ifndef REMOVE_EPI_SORTED
        
    #ifdef SUNWAY
        //#if 0
        unsigned long args_copy[3];
        args_copy[0] = (unsigned long) n_loc;
        args_copy[1] = (unsigned long) EPI_POINTER;
        args_copy[2] = (unsigned long) &system_grav[0];
        __real_athread_spawn((void*)slave_CopyEPIToFP, args_copy);
        athread_join();
    #else
        for(PS::S32 i=0; i<n_loc; i++){
            system_grav[i].pos = (((EPI*)(EPI_POINTER))+i)->pos;
            system_grav[i].mass = (((EPI*)(EPI_POINTER))+i)->mass;
            system_grav[i].vel = (((EPI*)(EPI_POINTER))+i)->vel;
            system_grav[i].id = (((EPI*)(EPI_POINTER))+i)->id;
        }
    #endif
#endif //REMOVE_EPI_SORTED
        
        endtime = PS::GetWtime();


        
        if (PS::Comm::getRank() == 0){
            for(PS::S32 i=0; i<5; i++){
                std::cerr<<"B) system_grav[i].pos= "<<system_grav[i].pos
                         <<" id= "<<system_grav[i].id
                         <<std::endl;
            }
        }
        
        
        if (PS::Comm::getRank() == 0)
            fout_wtime << "wtime_copy_epi_to_psys = " << endtime - starttime << std::endl;


#if defined(ROTATE_PTCLS) && defined(SHIFT_CENTER_Z)
     
    starttime = PS::GetWtime();
    RotateAndShiftPtcl(&system_grav[0],
                       system_grav.getNumberOfParticleLocal(),
                       -DT_TREE * n_loop_merge,
                       0.1*r_phy,
                       PS::TREE_Z_COORD_OFFSET);

    endtime = PS::GetWtime();
    if (PS::Comm::getRank() == 0)
        fout_wtime << "wtime_calc_shift_z = " << endtime - starttime << std::endl;
    #ifdef DEBUG_SHIFT_CENTER_Z
    if (PS::Comm::getRank() == 0)
        std::cerr<<"B) PS::TREE_Z_COORD_OFFSET(0)= "<<PS::TREE_Z_COORD_OFFSET<<std::endl;
    PS::TREE_Z_COORD_OFFSET = 0.0;
    PS::TREE_Z_COORD_OFFSET = CalcZCoordShift(&system_grav[0], system_grav.getNumberOfParticleLocal(), r_phy*0.1);
        if (PS::Comm::getRank() == 0)
        std::cerr<<"B) PS::TREE_Z_COORD_OFFSET(1)= "<<PS::TREE_Z_COORD_OFFSET<<std::endl;
    #endif


#else
    #ifdef ROTATE_PTCLS
        //* Rotate the coordinate of particles and satellites
        #ifdef SUNWAY
        starttime = PS::GetWtime();
        #ifdef DEBUG_0709
        theta_rot = - DT_TREE * n_loop_merge*10;
        #else
        theta_rot = - DT_TREE * n_loop_merge;
        #endif
        const PS::F64 cth = std::cos(theta_rot);
        const PS::F64 sth = std::sin(theta_rot);
        unsigned long args[4];
        //** Rotate the particles
        args[0] = (unsigned long) system_grav.getNumberOfParticleLocal();
        args[1] = (unsigned long) &system_grav[0];
        args[2] = (unsigned long) &cth;
        args[3] = (unsigned long) &sth;
        __real_athread_spawn((void*)slave_Rotate, args);
        athread_join();
        //** Rotate the satellites
        args[0] = (unsigned long) N_SATELLITE;
        args[1] = (unsigned long) SATELLITE.getPointer();
        args[2] = (unsigned long) &cth;
        args[3] = (unsigned long) &sth;
        __real_athread_spawn((void*)slave_Rotate, args);
        athread_join();
        endtime = PS::GetWtime();
        if (PS::Comm::getRank() == 0)
            fout_wtime << "wtime_rotate = " << endtime - starttime << std::endl;
        #else //SUNWAY
        theta_rot  = - DT_TREE * n_loop_merge;
        Rotate(system_grav,n_loc,theta_rot);
        Rotate(SATELLITE,N_SATELLITE,theta_rot);
        endtime = PS::GetWtime();
        if (PS::Comm::getRank() == 0)
            fout_wtime << "wtime_rotate = " << endtime - starttime << std::endl;
        #endif //SUNWAY
    #endif
    #ifdef SHIFT_CENTER_Z
        starttime = PS::GetWtime();
        z_coord_offset_loc = 0.0;
        for(PS::S32 i=0; i<n_loc; i++){
            if( z_coord_offset_loc < fabs(system_grav[i].pos.z) ) z_coord_offset_loc = fabs(system_grav[i].pos.z);
        }
        PS::TREE_Z_COORD_OFFSET = PS::Comm::getMaxValue(z_coord_offset_loc) + 1.0e-6;
        endtime = PS::GetWtime();
        if (PS::Comm::getRank() == 0)
            fout_wtime << "wtime_calc_shift_z = " << endtime - starttime << std::endl;
    #endif
#endif
        
#ifdef WRITE_SNP
        if(my_rank==0) WriteSnp(system_grav, snp_id, dir_name, 10000);
        PS::Comm::barrier();
        //PS::Finalize();
        //return 0;
#endif

        if (PS::Comm::getRank() == 0){
            for(PS::S32 i=0; i<5; i++){
                std::cerr<<"C) system_grav[i].pos= "<<system_grav[i].pos
                         <<" id= "<<system_grav[i].id
                         <<std::endl;
            }
        }
        
        //* Update time_sys
        time_sys += n_loop_merge*DT_TREE;
        
        //* Some tests 
        NAN_TEST(system_grav,"[nan-e]");
        ANGULAR_MOMENTUM_TEST(system_grav);


        PS::Comm::barrier();
        PS::F64 wtime_per_loop = PS::GetWtime() - wtime_offset;
        if(PS::Comm::getRank() == 0){            
            //std::cout<<"wtime_per_loop= "<<wtime_per_loop
            //         <<" wtime_per_loop/n_loop_merge= "<<wtime_per_loop/n_loop_merge
            //         <<std::endl;
            fout_wtime<<"wtime_per_loop= "<<wtime_per_loop
                      <<" wtime_per_loop/n_loop_merge= "<<wtime_per_loop/n_loop_merge
                      <<std::endl;
#if !defined(SUNWAY) || !defined(SUNWAY_FORCE_KERNEL)
            std::cout<<"(ENERGY_TOT-Etot0)/Etot0= "<<(ENERGY_TOT-Etot0)/Etot0<<std::endl;
#endif
        }
        n_loop++;
    }
    
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_close();
#endif

    if (PS::Comm::getRank() == 0) fout_wtime.close();
    
#ifdef SUNWAY
    athread_halt();
#endif

#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0){
        std::cerr<<"&system_grav[0]= "<<&system_grav[0]
                 <<" tree_grav.epi_sorted_.getPointer()= "<<tree_grav.epi_sorted_.getPointer()
                 <<" tree_grav.epj_sorted_loc_.getPointer()= "<<tree_grav.epj_sorted_loc_.getPointer()
                 <<std::endl;
    }
    PS::Comm::barrier();
#endif
    
    PS::Finalize();

    return 0;
}

