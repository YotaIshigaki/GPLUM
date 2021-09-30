#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/time.h>
#include<cmath>
#include<particle_simulator.hpp>
#ifdef ENABLE_PHANTOM_GRAPE_X86
#include <gp5util.h>
#endif
#ifdef ENABLE_GPU_CUDA
#define MULTI_WALK
#include"force_gpu_cuda.hpp"
#endif
#include "user-defined.hpp"
#include"force_sunway_impl.hpp"
#include"ic.hpp"
#include"collision.hpp"
#include"nbody.hpp"
#include"get_pos_domain.hpp"

#ifdef SUNWAY
extern "C"{
   #include <athread.h>
   #include "cpe_func.h"
   extern void SLAVE_FUN(Rotate)(void*);
}
#endif

static double wtime(){
   struct timeval tv;
   gettimeofday(&tv, NULL);
   return 1.e-6*tv.tv_usec +tv.tv_sec;
}

PS::F64 Epi::eps      = 1.0/32.0;
PS::F64 Epj::r_coll   = 0.5;
PS::F64 Epj::r_search = 0.5;
PS::F64 Planet::r_coll = 0.5; // identical to Satellite::r_coll; tetative value; modified later.
//PS::F64 DT_TREE = 1.0 / 64.0;
PS::F64 DT_TREE = 1.0 / 128.0;

Planet PLANET(1.0, PS::F64vec(0.0), PS::F64vec(0.0));
//PS::S32 N_SATELLITE = 0;
//PS::S32 N_SATELLITE = 1000;
//PS::S32 N_SATELLITE = 64;
PS::S32 N_SATELLITE = 128;
PS::ReallocatableArray<Satellite> SATELLITE;
PS::ReallocatableArray<SatelliteForce> SATELLITE_FORCE;
PS::ReallocatableArray<SatelliteForce> SATELLITE_FORCE_LOC;

std::ofstream fout_debug;

int main(int argc, char *argv[]) {
    double start_time = wtime();
    PS::Initialize(argc, argv);
    if (PS::Comm::getRank() == 0)
       std::cout << "etime(init) = " << wtime() - start_time << std::endl;
#ifdef SUNWAY    
    athread_init();
#endif
    if (PS::Comm::getRank() == 0)
       std::cout << "etime(1) = " << wtime() - start_time << std::endl;
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

    const PS::F64 PI = 4.0 * atan(1.0);
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::S32 my_rank = PS::Comm::getRank();
#ifdef COLLISION_TEST
    const PS::S32 n_loop_merge = 1;
#else
    const PS::S32 n_loop_merge = 32;
#endif
    //const PS::S32 n_smp_ave = 30;
    const PS::S32 n_smp_ave = 128;
    //const PS::S32 n_smp_ave = 512;
    //const PS::S32 n_smp_ave = 2048;
    //const PS::S32 n_smp_ave = 8192;
    
    PS::F32 theta = 0.5;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    //PS::F32 dt = 1.0 / 64.0;
    //PS::F32 time_end = DT_TREE*64;
    //PS::F32 time_end = DT_TREE*128;
    PS::F32 time_end = DT_TREE*320; // performance measurement
    //PS::F32 time_end = DT_TREE*512;
    //PS::F32 time_end = DT_TREE*2048; // 2.5 Kepler periods
    //PS::F32 time_end = DT_TREE*8192; // 10 Kepler periods
    //PS::F32 time_end = DT_TREE*32768; // 40 Kepler periods
    PS::F32 dt_diag = 1.0 / 8.0;
    PS::F32 dt_snap = 1.0;
    char dir_name[1024];
    PS::S64 n_glb = 1024;
    PS::F64 ax_cen = 1.0;
    //PS::F64 delta_ax = 0.01;
    PS::F64 delta_ax = 0.1;
    PS::F64 ax_in  = ax_cen - 0.5*delta_ax;
    PS::F64 ax_out = ax_cen + 0.5*delta_ax;
    PS::F64 ecc_rms = 0.0;
    PS::F64 inc_rms = 16.0;
    PS::F64 tau = 1.0;
    //PS::F64 tau = 1.0e-6;
    PS::F64 ratio_r_phy_r_hill = 1.0;
    PS::S32 c;
    PS::F64ort pos_domain;
    PS::F64 starttime,endtime;
#ifdef CHECK_NAN
    PS::S32 num_nan_loc,num_nan;
#endif

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
    makeOutputDirectory(dir_name);
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0)
       std::cout << "etime(2) = " << wtime() - start_time << std::endl;

    std::string s_debug = dir_name + ("/debug_" + ToS(PS::Comm::getRank()) + ".dat");
    fout_debug.open(s_debug.c_str());

    std::ofstream fout_eng;
    if(PS::Comm::getRank() == 0) {
        char sout_de[1024];
        sprintf(sout_de, "%s/t-de.dat", dir_name);
        fout_eng.open(sout_de);
        fprintf(stdout, "This is a sample program of N-body simulation on FDPS!\n");
        fprintf(stdout, "Number of processes: %d\n", PS::Comm::getNumberOfProc());
        fprintf(stdout, "Number of threads per process: %d\n", PS::Comm::getNumberOfThread());
    }
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0)
       std::cout << "etime(3) = " << wtime() - start_time << std::endl;


    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
#ifdef COLLISION_TEST
    //* Collision test: collision of two particles
    //* Choose test_mode
    PS::S32 test_mode = 0; // head-on collision test; no planet gravity
    //PS::S32 test_mode = 1; // two particles orbit around the planet
    //* Set up the particles
    n_glb = 2; // DO NOT CHANGE
    const PS::F64 tau_v = 1.0; // `v` means virtual
    const PS::S64 n_glb_v = (PS::S64) 1.0e11; // `v` means virtual
    const PS::F64 ratio_r_phy_r_hill_v = 1.0; // `v` means virtual
    PS::F64 r_phy = sqrt(tau_v*(ax_out*ax_out - ax_in*ax_in) / n_glb_v);
    PS::F64 r_hill = r_phy / ratio_r_phy_r_hill_v;
    PS::F64 mass_dust = (r_hill/ax_cen)*(r_hill/ax_cen)*(r_hill/ax_cen)*3.0*PLANET.mass*0.5;
    //PS::F64 dz = 3.0 * r_hill;
    PS::F64 dz = 20.0 * r_hill;
    //* Output to check
    if(PS::Comm::getRank() == 0){
        std::cout << "tau_v                = " << tau_v << std::endl;
        std::cout << "n_glb_v              = " << n_glb_v << std::endl;
        std::cout << "ratio_r_phy_r_hill_v = " << ratio_r_phy_r_hill_v << std::endl;
        std::cout << "r_phy                = " << r_phy << std::endl;
        std::cout << "r_hill               = " << r_hill << std::endl;
        std::cout << "dz                   = " << dz << std::endl;
        std::cout << "ax_cen               = " << ax_cen << std::endl;
        std::cout << "mass_dust            = " << mass_dust << std::endl;
    }
    CollisionTest_IC(system_grav,n_glb,mass_dust,ax_cen,dz,test_mode);
    Epj::r_search = r_hill * 10.0;
    Epj::r_coll   = r_phy;
    Epi::eps      = r_phy;

    //* Set up the satellites as mass-less particles 
    PS::MTTS mt;
    mt.init_genrand(0); 
    const PS::F64 large_factor = 1.0e16;
    SATELLITE.resizeNoInitialize(N_SATELLITE);
    for (PS::S32 i=0; i<N_SATELLITE; i++) {
        SATELLITE[i].mass  = 0.0;
        SATELLITE[i].pos.x = 2.0 * (mt.genrand_res53() - 0.5) * large_factor;
        SATELLITE[i].pos.y = 2.0 * (mt.genrand_res53() - 0.5) * large_factor;
        SATELLITE[i].pos.z = 2.0 * (mt.genrand_res53() - 0.5) * large_factor;
        SATELLITE[i].vel   = getVel(SATELLITE[i].pos);
        SATELLITE[i].id    = i + n_glb;
        // Note that the positions of the satellites are set to
        // large values. 
    }
    Planet::r_coll = Epj::r_coll; 
    SATELLITE_FORCE.resizeNoInitialize(64*N_SATELLITE); // We need large memory to use Long's force kenrel
    SATELLITE_FORCE_LOC.resizeNoInitialize(64*N_SATELLITE); // We need large memory to use Long's force kenrel
    for(PS::S32 i=0; i<SATELLITE_FORCE.size(); i++) SATELLITE_FORCE[i].clear();
    for(PS::S32 i=0; i<SATELLITE_FORCE_LOC.size(); i++) SATELLITE_FORCE_LOC[i].clear();

    //* Set up the planet
    if (test_mode == 0) {
        PLANET.mass = 0.0;
    }

    //* Compute # of particles and initialize time_sys
    PS::S32 n_loc = system_grav.getNumberOfParticleLocal();
    n_glb = system_grav.getNumberOfParticleGlobal();
    PS::F32 time_sys = 0.0;

    //* Set pos_domain (used for dinfo.initializePosDomain)
    GetPosDomain3(delta_ax,pos_domain); 

    //* For test run
    //std::cout << "n_loc = " << n_loc << " (myrank = " << PS::Comm::getRank() << ")" << std::endl;
    //fout_debug.close();
    //athread_halt();
    //PS::Finalize();
    //return 0;
#else
#if 1
    // planetary ring
    //tau  = n*PI*r_phy^2/(PI*(ax_out^2-ax_in^2))
    PS::F64 r_phy = sqrt(tau*(ax_out*ax_out - ax_in*ax_in) / n_glb);

    PS::F64 r_hill = r_phy / ratio_r_phy_r_hill;
    // r_hill = ax*qbrt(2*mass_dust / (3*mass_pla));
    PS::F64 mass_dust = (r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*3.0*PLANET.mass*0.5;
    PS::F32 time_sys = 0.0;
    PS::F64 power = 0.0;
    PS::S32 seed = 0;
    PS::F64 dens = mass_dust * n_glb / (PI*(ax_out*ax_out - ax_in*ax_in));
    //#######
    //---- Choice #1 -----
    //SetParticleKeplerDisk(system_grav, n_glb, ax_in, ax_out, ecc_rms, inc_rms, dens, PLANET.mass, 0.0, 1.0, power, seed);
    //---- Choice #2 -----
    //PS::F64 dz = 3.0 * r_hill;
    //PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(4) = " << wtime() - start_time << std::endl;
    //GetPosDomain3(delta_ax,pos_domain);
    //PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(5) = " << wtime() - start_time << std::endl;
    //SetParticleKeplerDisk2(system_grav, n_glb, ax_in, ax_out, dens, dz, pos_domain);
    //PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(6) = " << wtime() - start_time << std::endl;
    //---- Choice #3 -----
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(4) = " << wtime() - start_time << std::endl;
    GetPosDomain3(delta_ax,pos_domain);
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(5) = " << wtime() - start_time << std::endl;
    SetParticleKeplerDisk3(system_grav, n_glb, ax_in, ax_out, dens, pos_domain);
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(6) = " << wtime() - start_time << std::endl;
    //#######
    PS::S32 n_loc    = system_grav.getNumberOfParticleLocal();
    Epj::r_search = r_hill * 10.0;
    Epj::r_coll   = r_phy;
    Epi::eps      = r_phy;
    //Epi::eps      = 0.0;
    
    //if(PS::Comm::getRank() == 0){std::cout<<"system_grav[0].pos= "<<system_grav[0].pos<<std::endl;}
    ChangePtclToSatellite(system_grav, SATELLITE, N_SATELLITE);
    SATELLITE_FORCE.resizeNoInitialize(64*N_SATELLITE); // We need large memory to use Long's force kenrel
    SATELLITE_FORCE_LOC.resizeNoInitialize(64*N_SATELLITE); // We need large memory to use Long's force kenrel
    for(PS::S32 ip=0; ip<SATELLITE_FORCE.size(); ip++) SATELLITE_FORCE[ip].clear();
    for(PS::S32 ip=0; ip<SATELLITE_FORCE_LOC.size(); ip++) SATELLITE_FORCE_LOC[ip].clear();

    if(PS::Comm::getRank() == 0){
        std::cerr<<"dens= "<<dens<<" total mass= "<<mass_dust*n_glb<<" mass_dust= "<<mass_dust<<std::endl;
        if(N_SATELLITE > 0) std::cerr<<"SATELLITE[0].r_coll= "<<SATELLITE[0].r_coll<<std::endl;
        std::cerr<<"r_hill= "<<r_hill<<" mass_dust_total= "<<mass_dust * n_glb<<std::endl;
    }
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0)
       std::cout << "etime(7) = " << wtime() - start_time << std::endl;
    // remove particles overlapping with the satellite
    //RemovePartcile(system_grav, satellite);
    n_loc    = system_grav.getNumberOfParticleLocal();
    n_glb    = system_grav.getNumberOfParticleGlobal();

    //* For test run
    std::cout << "n_loc = " << n_loc << " (myrank = " << PS::Comm::getRank() << ")" << std::endl;
    //fout_debug.close();
    //athread_halt();
    //PS::Finalize();
    //return 0;
#elif 0 // for plummer model
    PS::S32 n_loc_tmp = 0;
    PS::F32 time_sys = 0.0;
    PS::ReallocatableArray<Satellite> satellite; // for compile
    SetParticlePlummer(system_grav, n_glb, n_loc_tmp, time_sys);
    PS::S32 n_loc    = system_grav.getNumberOfParticleLocal();
    //Epj::r_coll = 0.001;
    //Epj::r_search = Epj::r_coll * 10.0;
    Epj::r_coll = 0.0;
    Epj::r_search = 0.0;
    //Epj::r_search = 1.0;
    //Epj::r_search = 10.0;
    planet.mass = 0.0;

#else // for uniform sphere
    Epj::r_coll   = 1.0e-4;
    Epj::r_search = 10.0 * Epj::r_coll;

    PS::S32 n_loc_tmp = 0;
    PS::F32 time_sys = 0.0;
    SetParticleUniformSphere(system_grav, n_glb, n_loc_tmp, time_sys);

    ChangePtclToSatellite(system_grav, SATELLITE, N_SATELLITE);
    SATELLITE_FORCE.resizeNoInitialize(64*N_SATELLITE); // We need large memory to use Long's force kenrel
    SATELLITE_FORCE_LOC.resizeNoInitialize(64*N_SATELLITE); // We need large memory to use Long's force kenrel
    for(PS::S32 ip=0; ip<SATELLITE_FORCE.size(); ip++) SATELLITE_FORCE[ip].clear();
    for(PS::S32 ip=0; ip<SATELLITE_FORCE_LOC.size(); ip++) SATELLITE_FORCE_LOC[ip].clear();

    PS::S32 n_loc = system_grav.getNumberOfParticleLocal();
    n_glb = system_grav.getNumberOfParticleGlobal();

    //if (PS::Comm::getRank() == 0) {
    //    std::ofstream IC_file;
    //    std::string IC_filename = dir_name;
    //    IC_filename += "/IC.dat";
    //    IC_file.open(IC_filename.c_str());
    //    for (PS::S32 i=0; i<n_loc; i++) {
    //       IC_file << system_grav[i].pos.x << "  " 
    //               << system_grav[i].pos.y << "  "
    //               << system_grav[i].pos.z << std::endl;
    //    }
    //    IC_file.close();
    //    std::cout << "IC.dat is output." << std::endl;
    //}

    //* For test run
    //fout_debug.close();
    //athread_halt();
    //PS::Finalize();
    //return 0;
#endif
#endif // End of COLLISION_TEST

#ifdef CHECK_NAN
    num_nan_loc = 0;
    for (PS::S32 i=0; i<system_grav.getNumberOfParticleLocal(); i++) {
         if ( (std::isfinite(system_grav[i].pos.x) != true) ||
              (std::isfinite(system_grav[i].pos.y) != true) ||
              (std::isfinite(system_grav[i].pos.z) != true)) {
             num_nan_loc++;
         }
    }
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0)
        std::cout << "[nan-0] " << num_nan << std::endl;
#endif

    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "OK1 @main"  << std::endl;

    system_grav.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave);

    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "OK2 @main"  << std::endl;

    const PS::F32 coef_ema = 0.3;
    //const PS::F32 coef_ema = 0.5;
    //const PS::F32 coef_ema = 1.0;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);

    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "OK3 @main"  << std::endl;

    double wtime0 = PS::GetWtime();
#if 1
    //DomainDecision(dinfo, system_grav, n_smp_ave);

    PS::S32 nx = sqrt((PS::F64)n_proc-0.000001)+1;
    while( n_proc % nx != 0) nx++;
    PS::S32 ny = n_proc / nx;
    PS::S32 nz = 1;
    dinfo.setNumberOfDomainMultiDimension(nx, ny, nz);

    dinfo.initializePosDomain(pos_domain);

    //dinfo.collectSampleParticle(system_grav, true, 1.0);
    //dinfo.decomposeDomain2();
#else
    PS::S32 nx = sqrt((PS::F64)n_proc-0.000001)+1;
    while( n_proc % nx != 0) nx++;
    PS::S32 ny = n_proc / nx;
    PS::S32 nz = 1;
    dinfo.setNumberOfDomainMultiDimension(nx, ny, nz);
    PS::Comm::barrier();
    //if(PS::Comm::getRank()==0)std::cerr<<"CHECK B"<<std::endl;
    dinfo.decomposeDomainAll(system_grav);
#endif
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "OK4 @main"  << std::endl;
    double wtime_domain = PS::GetWtime() - wtime0;
    if(PS::Comm::getRank()==0)std::cerr<<"wtime_domain= "<<wtime_domain<<std::endl;


#if 0
    // now under construction
    // for efficient exchange particle
    F64ort my_domain = dinfo.getPosDomain(my_rank);
#endif
    
    double wtime1 = PS::GetWtime();
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "OK5 @main"  << std::endl;
    
#if 1
    //system_grav.exchangeParticle(dinfo);
#else
    system_grav.exchangeParticle2(dinfo);
#endif
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "OK6 @main"  << std::endl;
    system_grav.freeCommunicationBuffer();
    PS::Comm::barrier();
    //system_grav.dumpMemSizeUsed(fout_debug);

    //if(PS::Comm::getRank()==0)std::cerr<<"CHECK D"<<std::endl;
    double wtime_ex_ptcl = PS::GetWtime() - wtime1;
    //if(PS::Comm::getRank()==0)std::cerr<<"wtime_ex_ptcl= "<<wtime_ex_ptcl<<std::endl;
    
    n_loc = system_grav.getNumberOfParticleLocal();
    PS::Comm::barrier();
    //if(PS::Comm::getRank()==0)std::cerr<<"CHECK E"<<std::endl;
    //for(PS::S32 i=0; i<n_loc; i++) system_grav[i].rank_org = PS::Comm::getRank();
    //PS::Comm::barrier();if(PS::Comm::getRank()==0)std::cerr<<"CHECK F"<<std::endl;

    PS::S32 numPtclLocalMin,numPtclLocalMax;
    PS::S32 rankMin,rankMax;
    MPI::COMM_WORLD.Allreduce(&n_loc,&numPtclLocalMin,1,MPI::INT,MPI::MIN);
    MPI::COMM_WORLD.Allreduce(&n_loc,&numPtclLocalMax,1,MPI::INT,MPI::MAX);
    if (PS::Comm::getRank() == 0) {
       std::cout << "numPtclLocalMin = " << numPtclLocalMin << std::endl; 
       std::cout << "numPtclLocalMax = " << numPtclLocalMax << std::endl;   
    }
    if (n_loc == numPtclLocalMin)
        std::cout << "rankMin = " << PS::Comm::getRank() << std::endl;
    if (n_loc == numPtclLocalMax)
        std::cout << "rankMax = " << PS::Comm::getRank() << std::endl;

    //=========================================================
    // Setup constant parameters used in Long's force kernel
    //=========================================================
#if defined(SUNWAY) && defined(SUNWAY_FORCE_KERNEL)
    const PS::F64 Tdur = 5.0/32.0; 
    const PS::F64 ln_e_refl = std::log(0.5); 
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
        std::cout << "Tdur         = " << Tdur << std::endl;
        std::cout << "ln_e_refl    = " << ln_e_refl << std::endl;
        std::cout << "r_ep         = " << cpe_pars.r_ep << std::endl;
        std::cout << "r_sat        = " << cpe_pars.r_sat << std::endl;
        std::cout << "k_dash       = " << cpe_pars.kappa << std::endl;
        std::cout << "eta_dash     = " << cpe_pars.eta << std::endl;
        std::cout << "m_planet     = " << cpe_pars.m_planet << std::endl;
        std::cout << "n_loop_merge = " << n_loop_merge << std::endl;
        std::cout << "n_smp_ave    = " << n_smp_ave << std::endl;
        std::cout << "theta        = " << theta << std::endl;
    }


#if 1

    ///////////////
    /// chose TREE TYPE
    typedef Epi EPI;
    typedef Epj EPJ;
    typedef Force FORCE;
    //typedef PS::SPJMonopoleScatter SPJ;
    typedef SPJMonopoleScatterSW SPJ;
    //PS::TreeForForceLong<FORCE, EPI, EPJ>::MonopoleWithScatterSearch tree_grav;
    TreeType tree_grav;
    TREE_POINTER = & tree_grav;
    PS::Comm::barrier();if(PS::Comm::getRank()==0)std::cerr<<"OK7 @main"<<std::endl;
    PS::ReallocatableArray<Force> force_dust_from_satellite;
    tree_grav.initialize(n_glb, theta, n_leaf_limit, n_group_limit);
    PS::Comm::barrier();
    //system_grav.dumpMemSizeUsed(fout_debug);
    //tree_grav.dumpMemSizeUsed(fout_debug);

    PS::Comm::barrier();if(PS::Comm::getRank()==0)std::cerr<<"OK8 @main"<<std::endl;
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
    n_loc = system_grav.getNumberOfParticleLocal();
    n_glb = system_grav.getNumberOfParticleGlobal();
     for(PS::S32 i=0; i<n_loc; i++){
         system_grav[i].pos = (((EPI*)(EPI_POINTER))+i)->pos;
         system_grav[i].vel = (((EPI*)(EPI_POINTER))+i)->vel;
     }
     time_sys += DT_TREE;
    //system_grav.dumpMemSizeUsed(fout_debug);
    //tree_grav.dumpMemSizeUsed(fout_debug);
    //fout_debug.close();
    //athread_halt();
    //PS::Finalize();
    //return 0;

#if 0
    //* Compute the total angular momentum of the particles + the satellites.
    PS::F64vec L0_loc,L0_tot;
    L0_loc = 0.0;
    for (PS::S32 i=0; i<n_loc; i++) {
        L0_loc += system_grav[i].mass * (system_grav[i].pos ^ system_grav[i].vel);
    }
    for (PS::S32 i=0; i<N_SATELLITE; i++) {
        L0_loc += SATELLITE[i].mass * (SATELLITE[i].pos ^ SATELLITE[i].vel);
    }
    MPI::COMM_WORLD.Allreduce(&L0_loc.x,&L0_tot.x,1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&L0_loc.y,&L0_tot.y,1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&L0_loc.z,&L0_tot.z,1,MPI::DOUBLE,MPI::SUM);
    if (PS::Comm::getRank() == 0) {
       std::cout << "L0_tot.x = " << L0_tot.x << std::endl;
       std::cout << "L0_tot.y = " << L0_tot.y << std::endl;
       std::cout << "L0_tot.z = " << L0_tot.z << std::endl;
    }
#endif

    
#if 1 
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
    PS::Comm::barrier();if(PS::Comm::getRank()==0)std::cerr<<"OK11 @main"<<std::endl;

#if 1
    PS::F64 time_diag = 0.0;
    PS::F64 time_snap = 0.0;
    PS::S64 n_loop = 0;
    PS::S32 id_snap = 0;
    while(time_sys < time_end){
        if (PS::Comm::getRank() == 0) {
             //std::cout << "#######################" << std::endl;
             //std::cout << "  n_loop = " << n_loop << std::endl;
             //std::cout << "#######################" << std::endl;
             fout_debug << "#######################" << std::endl;
             fout_debug << "  n_loop = " << n_loop << std::endl;
             fout_debug << "#######################" << std::endl;
        }
        PS::Comm::barrier();
        PS::F64 wtime_offset = PS::GetWtime();


#ifdef CHECK_NAN
        num_nan_loc = 0;
        for (PS::S32 i=0; i<system_grav.getNumberOfParticleLocal(); i++) {
             if ( (std::isfinite(system_grav[i].pos.x) != true) ||
                  (std::isfinite(system_grav[i].pos.y) != true) ||
                  (std::isfinite(system_grav[i].pos.z) != true)) {
                 num_nan_loc++;
             }
        }
        MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
        if (PS::Comm::getRank() == 0)
            std::cout << "[nan-a] " << num_nan << " @n_loop = " << n_loop << std::endl;
#endif

        //dinfo.decomposeDomainAll(system_grav); // test
        //------
        //dinfo.collectSampleParticle(system_grav, true, 1.0);
        //dinfo.decomposeDomain2();
        //------

#ifdef CHECK_NAN
        num_nan_loc = 0;
        for (PS::S32 i=0; i<system_grav.getNumberOfParticleLocal(); i++) {
             if ( (std::isfinite(system_grav[i].pos.x) != true) ||
                  (std::isfinite(system_grav[i].pos.y) != true) ||
                  (std::isfinite(system_grav[i].pos.z) != true)) {
                 num_nan_loc++;
             }
        }
        MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
        if (PS::Comm::getRank() == 0)
            std::cout << "[nan-b] " << num_nan << " @n_loop = " << n_loop << std::endl;
#endif

        //PS::Comm::barrier();
        //if(PS::Comm::getRank()==0)std::cerr<<"OK12-"<<n_loop<<" @main"<<std::endl;

        //system_grav.exchangeParticle(dinfo);
        system_grav.exchangeParticle2(dinfo);
        //system_grav.exchangeParticle3(dinfo);

        //PS::Comm::barrier();
        //if(PS::Comm::getRank()==0)std::cerr<<"OK13-"<<n_loop<<" @main"<<std::endl;

#ifdef CHECK_NAN
        num_nan_loc = 0;
        for (PS::S32 i=0; i<system_grav.getNumberOfParticleLocal(); i++) {
             if ( (std::isfinite(system_grav[i].pos.x) != true) ||
                  (std::isfinite(system_grav[i].pos.y) != true) ||
                  (std::isfinite(system_grav[i].pos.z) != true)) {
                 num_nan_loc++;
             }
        }
        MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
        if (PS::Comm::getRank() == 0)
            std::cout << "[nan-c] " << num_nan << " @n_loop = " << n_loop << std::endl;
#endif

        n_loc = system_grav.getNumberOfParticleLocal();
#if 0
        MPI::COMM_WORLD.Allreduce(&n_loc,&numPtclLocalMin,1,MPI::INT,MPI::MIN);
        MPI::COMM_WORLD.Allreduce(&n_loc,&numPtclLocalMax,1,MPI::INT,MPI::MAX);
        if (PS::Comm::getRank() == 0) {
           std::cout << "numPtclLocalMin = " << numPtclLocalMin << std::endl;   
           std::cout << "numPtclLocalMax = " << numPtclLocalMax << std::endl;   
        }
        if (n_loc == numPtclLocalMin)
            std::cout << "rankMin = " << PS::Comm::getRank() << std::endl;
        if (n_loc == numPtclLocalMax)
            std::cout << "rankMax = " << PS::Comm::getRank() << std::endl;
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
        //system_grav.dumpMemSizeUsed(fout_debug);
        //tree_grav.dumpMemSizeUsed(fout_debug);
        //PS::Comm::barrier();
        //fout_debug.close();
        //athread_halt();
        //PS::Finalize();
        //return 0;

#ifdef CHECK_NAN
        num_nan_loc = 0;
        for (PS::S32 i=0; i<system_grav.getNumberOfParticleLocal(); i++) {
             if ( (std::isfinite(system_grav[i].pos.x) != true) ||
                  (std::isfinite(system_grav[i].pos.y) != true) ||
                  (std::isfinite(system_grav[i].pos.z) != true)) {
                 num_nan_loc++;
             }
        }
        MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
        if (PS::Comm::getRank() == 0)
            std::cout << "[nan-d] " << num_nan << " @n_loop = " << n_loop << std::endl;
#endif

        starttime = PS::GetWtime();
        for(PS::S32 i=0; i<n_loc; i++){
            system_grav[i].pos = (((EPI*)(EPI_POINTER))+i)->pos;
            system_grav[i].vel = (((EPI*)(EPI_POINTER))+i)->vel;
        }
        endtime = PS::GetWtime();
        if (PS::Comm::getRank() == 0)
            fout_debug << "wtime_copy_epi_to_psys = " << endtime - starttime << std::endl;

#ifdef ROTATE_PTCLS
        //* Rotate the coordinate of particles and satellites
#ifdef SUNWAY
        //std::ofstream fp;
        //fp.open("before.dat",std::ios::trunc);
        //for (PS::S32 i=0; i<n_loc; i++) {
        //   fp << system_grav[i].pos.x << "  "
        //      << system_grav[i].pos.y << "  "
        //      << system_grav[i].pos.z << "  "
        //      << std::endl;
        //}
        //fp.close();

        //const PS::F64 theta = - DT_TREE * n_loop_merge;
        //Rotate(system_grav,n_loc,theta);
        //Rotate(SATELLITE,N_SATELLITE,theta);
        //
        //fp.open("after_org.dat",std::ios::trunc);
        //for (PS::S32 i=0; i<n_loc; i++) {
        //   fp << system_grav[i].pos.x << "  "
        //      << system_grav[i].pos.y << "  "
        //      << system_grav[i].pos.z << "  "
        //      << std::endl;
        //}
        //fp.close();

        starttime = PS::GetWtime();
        PS::F64 theta = - DT_TREE * n_loop_merge;
        if (n_loop == 0) theta = - DT_TREE * (n_loop_merge+1);
        const PS::F64 cth = std::cos(theta);
        const PS::F64 sth = std::sin(theta);
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
            fout_debug << "wtime_rotate = " << endtime - starttime << std::endl;

        //fp.open("after_new.dat",std::ios::trunc);
        //for (PS::S32 i=0; i<n_loc; i++) {
        //   fp << system_grav[i].pos.x << "  "
        //      << system_grav[i].pos.y << "  "
        //      << system_grav[i].pos.z << "  "
        //      << std::endl;
        //}
        //fp.close();
#else
        PS::F64 theta  = - DT_TREE * n_loop_merge;
        if (n_loop == 0) theta = - DT_TREE * (n_loop_merge+1);
        Rotate(system_grav,n_loc,theta);
        Rotate(SATELLITE,N_SATELLITE,theta);
        endtime = PS::GetWtime();
        if (PS::Comm::getRank() == 0)
            fout_debug << "wtime_rotate = " << endtime - starttime << std::endl;
#endif
#endif
        
        //* Update time_sys
        time_sys += n_loop_merge*DT_TREE;

#ifdef CHECK_NAN
        num_nan_loc = 0;
        for (PS::S32 i=0; i<system_grav.getNumberOfParticleLocal(); i++) {
             if ( (std::isfinite(system_grav[i].pos.x) != true) ||
                  (std::isfinite(system_grav[i].pos.y) != true) ||
                  (std::isfinite(system_grav[i].pos.z) != true)) {
                 num_nan_loc++;
             }
        }
        MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
        if (PS::Comm::getRank() == 0)
            std::cout << "[nan-e] " << num_nan << std::endl;
#endif

#if 0
        //* rmin,rmax
        PS::F64 r2min_loc,r2min,r2max_loc,r2max;
        r2min_loc =  1.0e300;
        r2max_loc = -1.0e300;
        for (PS::S32 i=0; i<system_grav.getNumberOfParticleLocal(); i++) {
             PS::F64vec pos = system_grav[i].pos;
             PS::F64 r2 = pos * pos;
             if (r2 < r2min_loc) r2min_loc = r2;
             if (r2 > r2max_loc) r2max_loc = r2;
            
        }
        MPI::COMM_WORLD.Allreduce(&r2min_loc,&r2min,1,MPI::DOUBLE,MPI::MIN);
        MPI::COMM_WORLD.Allreduce(&r2max_loc,&r2max,1,MPI::DOUBLE,MPI::MAX);
        //std::cout << "r2min_loc = " << r2min_loc << "  " << PS::Comm::getRank() << std::endl;
        //std::cout << "r2max_loc = " << r2max_loc << "  " << PS::Comm::getRank() << std::endl;
        if (PS::Comm::getRank() == 0) {
            std::cout << "r2min = " << r2min << std::endl;
            std::cout << "r2max = " << r2max << std::endl;
        }
#endif

#if 0
        //* Compute the total angular momentum of the particles + the satellites.
        PS::F64vec L_loc,L_tot;
        L_loc = 0.0;
        for (PS::S32 i=0; i<n_loc; i++) {
            L_loc += system_grav[i].mass * (system_grav[i].pos ^ system_grav[i].vel);
        }
        for (PS::S32 i=0; i<N_SATELLITE; i++) {
            L_loc += SATELLITE[i].mass * (SATELLITE[i].pos ^ SATELLITE[i].vel);
        }
        MPI::COMM_WORLD.Allreduce(&L_loc.x,&L_tot.x,1,MPI::DOUBLE,MPI::SUM);
        MPI::COMM_WORLD.Allreduce(&L_loc.y,&L_tot.y,1,MPI::DOUBLE,MPI::SUM);
        MPI::COMM_WORLD.Allreduce(&L_loc.z,&L_tot.z,1,MPI::DOUBLE,MPI::SUM);
        if (PS::Comm::getRank() == 0) {
           std::cout << "L_tot.x = " << L_tot.x << " " 
                     << "rel.err. = " << std::abs((L_tot.x-L0_tot.x)/L0_tot.x) << std::endl;
           std::cout << "L_tot.y = " << L_tot.y << " " 
                     << "rel.err. = " << std::abs((L_tot.y-L0_tot.y)/L0_tot.y) << std::endl;
           std::cout << "L_tot.z = " << L_tot.z << " "
                     << "rel.err. = " << std::abs((L_tot.z-L0_tot.z)/L0_tot.z) << std::endl;
        }
#endif

#if 0
        //* Output ptcl data (binary)
        if (((n_loop+1)%1 == 0)) {
           std::ofstream ofs;
           std::stringstream fnum1,fnum2;
           fnum1 << ToS(id_snap);
           fnum2 << ToS(PS::Comm::getRank());
           std::string fname;
           fname = "./output/ptcl" + fnum1.str() + "-" + fnum2.str() + ".dat";
           ofs.open(fname.c_str(),std::ios::binary|std::ios::trunc);
           int numPtcl = system_grav.getNumberOfParticleLocal();
           ofs.write((char *) &numPtcl, sizeof(int));
           for (PS::S32 i=0; i<numPtcl; i++) {
               ofs.write((char *) &system_grav[i].pos, 3*sizeof(double));
           }
           ofs.close();
           id_snap++;
        }
#endif

#ifdef COLLISION_TEST
        //* Output ptcl data
        std::ofstream ofs;
        std::stringstream fnum;
        fnum << ToS(n_loop);
        std::string fname;
        fname = "./output/ptcl" + fnum.str() + ".txt";
        ofs.open(fname.c_str(),std::ios::trunc);
        for (PS::S32 i=0; i<system_grav.getNumberOfParticleLocal(); i++) {
            ofs << system_grav[i].pos.x << "  "
                << system_grav[i].pos.y << "  "
                << system_grav[i].pos.z << "  "
                << system_grav[i].vel.x << "  "
                << system_grav[i].vel.y << "  "
                << system_grav[i].vel.z << "  "
                << std::endl;
        }
        PS::F64 dl = (system_grav[1].pos.x - system_grav[0].pos.x)/(2.0*r_phy);
        PS::F64 vx_sum = system_grav[0].vel.x + system_grav[1].vel.x;
        ofs << "# (dl,mom.tot.) = " << dl << " " << vx_sum << std::endl;
        ofs.close();
#endif

        PS::Comm::barrier();
        PS::F64 wtime_per_loop = PS::GetWtime() - wtime_offset;
        PS::CountT n_proc_ep_send_tot  = PS::Comm::getSum(tree_grav.comm_table_.n_proc_ep_send_);
        PS::CountT n_proc_ep_recv_tot  = PS::Comm::getSum(tree_grav.comm_table_.n_proc_ep_recv_);
        PS::CountT n_proc_sp_isend_tot = PS::Comm::getSum(tree_grav.comm_table_.n_proc_sp_isend_);
        PS::CountT n_proc_sp_irecv_tot = PS::Comm::getSum(tree_grav.comm_table_.n_proc_sp_irecv_);
        if(PS::Comm::getRank() == 0){
            fout_debug<<"dinfo.pos_domain[my_rank].low_= "<<dinfo.getPosDomain(my_rank).low_
                      <<" dinfo.pos_domain[my_rank].high_= "<<dinfo.getPosDomain(my_rank).high_
                      <<std::endl;
            fout_debug<<"n_proc_ep_send_= "<<tree_grav.comm_table_.n_proc_ep_send_
                      <<" n_proc_ep_recv_= "<<tree_grav.comm_table_.n_proc_ep_recv_
                      <<std::endl;
            fout_debug<<"n_proc_sp_isend_= "<<tree_grav.comm_table_.n_proc_sp_isend_
                      <<" n_proc_sp_irecv_= "<<tree_grav.comm_table_.n_proc_sp_irecv_
                      <<std::endl;
            fout_debug<<"n_proc_ep_send_tot= "<<n_proc_ep_send_tot
                      <<" n_proc_ep_send_tot/n_proc= "<<((PS::F64)n_proc_ep_send_tot)/n_proc
                      <<std::endl;
            fout_debug<<"n_proc_sp_isend_tot= "<<n_proc_sp_isend_tot
                      <<" n_proc_sp_isend_tot/n_proc= "<<((PS::F64)n_proc_sp_isend_tot)/n_proc
                      <<std::endl;
            fout_debug<<"n_proc_sp_irecv_tot= "<<n_proc_sp_irecv_tot
                      <<" n_proc_sp_irecv_tot/n_proc= "<<((PS::F64)n_proc_sp_irecv_tot)/n_proc
                      <<std::endl;
        }
        if(PS::Comm::getRank() == 0){            
            //std::cout<<"wtime_per_loop= "<<wtime_per_loop
            //         <<" wtime_per_loop/n_loop_merge= "<<wtime_per_loop/n_loop_merge
            //         <<std::endl;
            fout_debug<<"wtime_per_loop= "<<wtime_per_loop
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

#endif

#endif

#endif
    fout_debug.close();
#ifdef SUNWAY
    athread_halt();
#endif
    PS::Finalize();
    return 0;
}

