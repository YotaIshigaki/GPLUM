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

PS::F64 MASS_DUST_TOTAL = 0.0;
PS::F64vec CM_POS_DUST_TOTAL = 0.0;
PS::F64 Epi::eps      = 1.0/32.0;
PS::F64 Epj::r_coll   = 0.5;
PS::F64 Epj::r_search = 0.5;
PS::F64 Planet::r_coll = 0.5; // identical to Satellite::r_coll; tetative value; modified later.
//PS::F64 DT_TREE = 1.0 / 64.0;
//PS::F64 DT_TREE = 1.0 / 128.0;
PS::F64 DT_TREE = 1.0 / 256.0;
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
#ifdef COLLISION_TEST
    const PS::S32 n_loop_merge = 1;
#else
    //const PS::S32 n_loop_merge = 1;
    const PS::S32 n_loop_merge = 3;
    //const PS::S32 n_loop_merge = 8;
    //const PS::S32 n_loop_merge = 16;
    //const PS::S32 n_loop_merge = 32;
    //const PS::S32 n_loop_merge = 64;
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
    PS::F32 time_end = DT_TREE*320; // performance measurement
    //PS::F32 time_end = DT_TREE*512;
    //PS::F32 time_end = DT_TREE*1024;
    //PS::F32 time_end = DT_TREE*2048; // 2.5 Kepler periods
    //PS::F32 time_end = DT_TREE*4096;
    //PS::F32 time_end = DT_TREE*8192; // 10 Kepler periods
    //PS::F32 time_end = DT_TREE*32768; // 40 Kepler periods
    PS::F32 dt_diag = 1.0 / 8.0;
    PS::F32 dt_snap = 1.0;
    char dir_name[1024];
    //PS::S64 n_glb = 1024, n_loc;
    PS::S64 n_glb = 1024;
    PS::S64 n_loc = 0;
    PS::F64 ax_cen = 1.0;
    //PS::F64 delta_ax = 0.1;
    PS::F64 delta_ax = 0.01;
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
       std::cout << "etime(mkdir) = " << wtime() - wtime_main_0 << std::endl;

    //* Open fout_debug
    wtime_main_0 = wtime();
    std::string s_debug = dir_name + ("/debug_" + ToS(PS::Comm::getRank()) + ".dat");
    if (PS::Comm::getRank() == 0) {
        std::cerr<<"s_debug: "<<s_debug<<std::endl;
        fout_debug.open(s_debug.c_str());
        if(!fout_debug.is_open()) exit(1);
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

    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();

#ifdef PHI_R_TREE
    /*
    PS::S32 ny = 2;
    PS::S32 nz = 1;
    PS::S32 nx = n_proc / ny;
    */
    /*
    PS::S32 ny = delta_ax / sqrt((2.0*PI*delta_ax)/((PS::F64)n_proc));
    while( n_proc % ny != 0) ny++;
    if(ny <= 0) ny = 1;
    PS::S32 nx = n_proc / ny;
    PS::S32 nz = 1;
    */

    PS::S32 ny = 1;
    PS::S32 nx = n_proc / ny;
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
    PS::S32 nz = 1;
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0){
        std::cerr<<"nx= "<<nx<<" ny= "<<ny<<" nz= "<<nz<<std::endl;
        fout_debug<<"nx= "<<nx<<" ny= "<<ny<<" nz= "<<nz<<std::endl;
    }

#elif 1
    PS::S32 nx = sqrt((PS::F64)n_proc-0.000001)+1;
    while( n_proc % nx != 0) nx++;
    PS::S32 ny = n_proc / nx;
    PS::S32 nz = 1;
#else
    PS::S32 nx = 400, ny = 396, nz = 1; 
#endif
    
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
    if (test_mode == 0) PLANET.mass = 0.0;

    //* Compute # of particles and initialize time_sys
    n_loc = system_grav.getNumberOfParticleLocal();
    n_glb = system_grav.getNumberOfParticleGlobal();
    //PS::F32 time_sys = 0.0; // MOVE TO GLOBAL VALUABE
    time_sys = 0.0; // MOVE TO GLOBAL VALUABE
    
    //* Set pos_domain (used for dinfo.initializePosDomain)
    GetPosDomain3(delta_ax,pos_domain); 

    //* For test run
    //std::cout << "n_loc = " << n_loc << " (myrank = " << PS::Comm::getRank() << ")" << std::endl;
    //fout_debug.close();
    //athread_halt();
    //PS::Finalize();
    //return 0;
#else // NOT COLLISION_TEST
    #if 1
    // planetary ring
    //tau  = n*PI*r_phy^2/(PI*(ax_out^2-ax_in^2))
    PS::F64 r_phy = sqrt(tau*(ax_out*ax_out - ax_in*ax_in) / n_glb);

    PS::F64 r_hill = r_phy / ratio_r_phy_r_hill;
    // r_hill = ax*qbrt(2*mass_dust / (3*mass_pla));

    #ifdef ENERGY_CHECK
    //PS::F64 mass_dust = (r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*3.0*PLANET.mass*0.5 * 1e9;
    PS::F64 mass_dust = (r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*3.0*PLANET.mass*0.5 * 1e14;
    #else
    PS::F64 mass_dust = (r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*3.0*PLANET.mass*0.5;
    //PS::F64 mass_dust = (r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*(r_hill/((ax_out+ax_in)*0.5))*3.0*PLANET.mass*0.5 * 1e9;
    #endif
    
    time_sys = 0.0;// MOVE TO GLOBAL VALUABE
    PS::F64 power = 0.0;
    PS::S32 seed = 0;
    PS::F64 dens = mass_dust * n_glb / (PI*(ax_out*ax_out - ax_in*ax_in));

    wtime_main_0 = wtime();

    #ifdef PHI_R_TREE
    //#if 0
    GetPosDomainCyl(delta_ax, pos_domain, nx, ny);
    SetParticleKeplerDiskCyl2(system_grav, n_glb, ax_in, ax_out, dens, pos_domain, true);
    #else //PHI_R_TREE
        #if 1
    GetPosDomain3(delta_ax,pos_domain);
        #else
    GetPosDomain3(delta_ax,pos_domain,400,396);
        #endif
    SetParticleKeplerDisk3(system_grav, n_glb, ax_in, ax_out, dens, pos_domain);
    #endif //PHI_R_TREE
    n_loc = system_grav.getNumberOfParticleLocal();
    n_glb = system_grav.getNumberOfParticleGlobal();

    //////////
    // set offset of tree center in z direction
    #ifdef SHIFT_CENTER_Z
    PS::F64 z_coord_offset_loc = 0.0;
    for(PS::S32 i=0; i<n_loc; i++){
        if( z_coord_offset_loc < fabs(system_grav[i].pos.z) ) z_coord_offset_loc = fabs(system_grav[i].pos.z);
    }
        #if 1
    PS::TREE_Z_COORD_OFFSET = PS::Comm::getMaxValue(z_coord_offset_loc) + 1.0e-6;
    if(PS::Comm::getRank() == 0) std::cerr<<"PS::TREE_Z_COORD_OFFSET= "<<PS::TREE_Z_COORD_OFFSET<<std::endl;
        #else
    PS::TREE_Z_COORD_OFFSET = 1.0e-8;
    for(PS::S32 i=0; i<n_loc; i++){
        if(system_grav[i].pos.z < 0.0) system_grav[i].pos.z = -1.0e-10;
        if(system_grav[i].pos.z > 0.0) system_grav[i].pos.z = 1.0e-10;
    }
        #endif
    #endif
    //PS::TREE_Z_COORD_OFFSET = PS::Comm::getMaxValue(z_coord_offset_loc) * 2.1;
    //if(PS::Comm::getRank() == 0) std::cerr<<"PS::TREE_Z_COORD_OFFSET= "<<PS::TREE_Z_COORD_OFFSET<<std::endl;
    //PS::TREE_Z_COORD_OFFSET = 0.0;
    // set offset of tree center in z direction
    //////////
    
        #if defined(PHI_R_TREE) || !defined(PHI_R_TREE)
    PS::S32 * n_loc_array = new PS::S32[n_proc];
    MPI_Allgather(&n_loc, 1, MPI_INT, n_loc_array, 1, MPI_INT, MPI_COMM_WORLD);
    PS::S32 id_offset = 0;
    for(PS::S32 i=0; i<my_rank; i++){
        id_offset += n_loc_array[i];
    }
    for(PS::S32 i=0; i<n_loc; i++){
        system_grav[i].id = id_offset + i;
    }
        #endif
    
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cout << "etime(IC) = " << wtime() - wtime_main_0 << std::endl;

    //Epj::r_search = r_hill * 10.0;
    Epj::r_search = r_hill * 3.0;
    //Epj::r_search = r_hill;
    Epj::r_coll   = r_phy;
    //Epi::eps      = r_phy;
    Epi::eps      = 0.0;
    
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

    
    //* Remove particles overlapping with the satellite
    //RemovePartcile(system_grav, satellite);
    //n_loc = system_grav.getNumberOfParticleLocal();

    //* For test run
    //std::cout << "n_loc = " << n_loc << " (myrank = " << PS::Comm::getRank() << ")" << std::endl;
    //fout_debug.close();
    //athread_halt();
    //PS::Finalize();
    //return 0;
    
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
    #endif // for uniform sphere
#endif // End of COLLISION_TEST

    NAN_TEST(system_grav,"[nan-0]");
#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "OK1 @main"  << std::endl;
#endif
    system_grav.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave);
#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "OK2 @main"  << std::endl;
#endif
    const PS::F32 coef_ema = 0.3;
    //const PS::F32 coef_ema = 0.5;
    //const PS::F32 coef_ema = 1.0;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "OK3 @main"  << std::endl;
#endif
    /*
#ifdef PHI_R_TREE
    PS::S32 ny = 2;
    PS::S32 nx = n_proc / ny;
    PS::S32 nz = 1;    
#elif 1
    PS::S32 nx = sqrt((PS::F64)n_proc-0.000001)+1;
    while( n_proc % nx != 0) nx++;
    PS::S32 ny = n_proc / nx;
    PS::S32 nz = 1;
#else
    PS::S32 nx = 400, ny = 396, nz = 1; 
#endif
    */
    dinfo.setNumberOfDomainMultiDimension(nx, ny, nz);
    if(PS::Comm::getRank() == 0){
        std::cerr<<"dinfo.getNDomain(0)= "<<dinfo.getNDomain(0)
                 <<" dinfo.getNDomain(1)= "<<dinfo.getNDomain(1)
                 <<" dinfo.getNDomain(2)= "<<dinfo.getNDomain(2)
                 <<std::endl;
    }
    n_loc = system_grav.getNumberOfParticleLocal();
#if 0
    for(PS::S32 i=0; i<n_loc; i++){
        system_grav[i].r_init = sqrt(system_grav[i].pos*system_grav[i].pos);
    }
#endif
    
#ifdef PHI_R_TREE
    ///////////////////////////////
    ///////// set phi and r in FP
    #if 0
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64 phi = atan2(system_grav[i].pos.y, system_grav[i].pos.x);
        if(phi < 0.0) phi += 2.0*PI;
        PS::F64 r = sqrt(system_grav[i].pos.x*system_grav[i].pos.x + system_grav[i].pos.y*system_grav[i].pos.y);
        system_grav[i].pos.x = phi;
        system_grav[i].pos.y = r;
    }
    #endif
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_X);
    dinfo.setPosRootDomain(PS::F64vec(0.0, -PI, -PI), PS::F64vec(2.0*PI, PI, PI)); // phi, r, z
    //system_grav.adjustPositionIntoRootDomain(dinfo);
    wtime_main_0 = wtime();
    dinfo.collectSampleParticle(system_grav, true, 1.0);
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(collect) = " << wtime() - wtime_main_0 << std::endl;
    wtime_main_0 = wtime();
    //dinfo.decomposeDomainMultiStep2(true, true, 1.0);
    //dinfo.decomposeDomainMultiStep2(true, false); // 1st: sample sort, 2nd split y-coord
    dinfo.decomposeDomainMultiStep3(true, false); // 1st: sample sort, 2nd split y-coord

#ifdef USE_SUPER_DOMAIN    
    MY_MPI_COMM_SUB = dinfo.getCommSub(0);
    MY_MPI_COMM_1D  = dinfo.getComm1d(0); 
    MY_RANK_SUB     = dinfo.getRankSub(0);
#endif
    
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(multistep2) = " << wtime() - wtime_main_0 << std::endl;
    PS::Comm::barrier();
    /*
    if (PS::Comm::getRank() == 0){
        for(PS::S32 i=0; i<n_proc; i++){
            std::cerr<<"i= "<<i<<" PosDomain(i)= "<<dinfo.getPosDomain(i)<<std::endl;
        }
    }
    */
    //exit(1);
    
    //system_grav.exchangeParticle4(dinfo);
    //system_grav.exchangeParticle5(dinfo);
    wtime_main_0 = wtime();

    n_loc = system_grav.getNumberOfParticleLocal();
    #if defined(DEBUG_EX_PTCL7)
    for(PS::S32 i=0; i<n_loc; i += (n_loc/1000) ){
        PS::F64 x_new = dinfo.getPosDomain(my_rank).low_.x-1e-6;
        if(x_new >= 0.0) system_grav[i].pos.x = x_new;
    }
    for(PS::S32 i=50; i<n_loc; i += (n_loc/1000) ){
        PS::F64 x_new = dinfo.getPosDomain(my_rank).high_.x+1e-6;
        if(x_new >= 0.0) system_grav[i].pos.x = x_new;
    }
    #endif //DEBUG_EX_PTCL7
    #if defined(DEBUG_EX_PTCL8)
    PS::S32 n_cnt = 0;
    for(PS::S32 i=(n_loc/1000); i<n_loc; i += (n_loc/1000), n_cnt++ ){
        double dx = 1e-5;
        PS::F64 phi_new = dinfo.getPosDomain(my_rank).high_.x+dx;
        if(n_cnt%2==1) phi_new = dinfo.getPosDomain(my_rank).low_.x-dx;
        PS::F64 r_new = 1.0;
        system_grav[i].pos.x = r_new*cos(phi_new);
        system_grav[i].pos.y = r_new*sin(phi_new);
    }
    #endif
    //system_grav.exchangeParticle6(dinfo);
    //system_grav.exchangeParticle7(dinfo);
    system_grav.exchangeParticle8(dinfo);
    
    n_loc = system_grav.getNumberOfParticleLocal();
    PS::Comm::barrier();
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64 phi = atan2(system_grav[i].pos.y, system_grav[i].pos.x);
        if(phi < 0.0) phi += 2.0*4.0*atan(1.0);
        if(phi >= dinfo.getPosDomain(my_rank).high_.x || phi < dinfo.getPosDomain(my_rank).low_.x){
            std::cerr<<"rank="<<my_rank<<" phi= "<<phi<<" pos= "<<system_grav[i].pos<<std::endl;
        }
        assert(phi < dinfo.getPosDomain(my_rank).high_.x && phi >= dinfo.getPosDomain(my_rank).low_.x);
    }
    //std::cerr<<"rank= "<<my_rank<<" n_loc= "<<n_loc<<std::endl;
    //exit(1);
    
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(exch_ptcl) = " << wtime() - wtime_main_0 << std::endl;
    n_loc = system_grav.getNumberOfParticleLocal();
    #ifdef USE_SUPER_DOMAIN
    wtime_main_0 = wtime();
    dinfo.initializeSuperDomain();
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "etime(superdomain) = " << wtime() - wtime_main_0 << std::endl;
    //exit(1);
    #endif

    // for debug
    /*
    if(PS::Comm::getRank() == 0){
        for(PS::S32 i=0; i<n_proc; i++){
            std::cerr<<"i= "<<i<<" dinfo.getPosDomain(i)= "<<dinfo.getPosDomain(i)<<std::endl;
        }
    }
    athread_halt();
    PS::Finalize();
    return 0;    
    */
#else //PHI_R_TREE
    dinfo.initializePosDomain(pos_domain);
#endif //PHI_R_TREER_PHI_TREE
#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier(); if (PS::Comm::getRank() == 0) std::cout << "OK4 @main"  << std::endl;
#endif
    system_grav.freeCommunicationBuffer();
    //system_grav.dumpMemSizeUsed(fout_debug);
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
#ifdef DEBUG_PRINT_MAIN
    PS::Comm::barrier();if(PS::Comm::getRank()==0)std::cerr<<"OK7 @main"<<std::endl;
#endif
    PS::ReallocatableArray<Force> force_dust_from_satellite;
    tree_grav.initialize(n_glb, theta, n_leaf_limit, n_group_limit);
    PS::Comm::barrier();
    //system_grav.dumpMemSizeUsed(fout_debug);
    //tree_grav.dumpMemSizeUsed(fout_debug);
    
#ifdef CHECK_EX_PTCL
    n_loc = system_grav.getNumberOfParticleLocal();

    for(PS::S32 ip=0; ip<n_loc; ip++){
        PS::F64 r_inv = 1.0 / sqrt(system_grav[ip].pos*system_grav[ip].pos);
        PS::F64vec acc = -r_inv * r_inv * r_inv * system_grav[ip].pos;
        system_grav[ip].vel += acc*DT_TREE*0.5;
        system_grav[ip].pos += system_grav[ip].vel*DT_TREE;
        system_grav[ip].vel += acc*DT_TREE*0.5;
    }
    //system_grav.exchangeParticle4(dinfo);
    //system_grav.exchangeParticle5(dinfo);
    //system_grav.exchangeParticle6(dinfo);
    //system_grav.exchangeParticle7(dinfo);
    system_grav.exchangeParticle8(dinfo);
    
    n_loc = system_grav.getNumberOfParticleLocal();
    n_glb = system_grav.getNumberOfParticleGlobal();

    //std::cerr<<"rank= "<<my_rank<<" n_loc= "<<n_loc<<std::endl;
    
    PS::TimeProfile prof_tree    = TREE_POINTER->getTimeProfile();
    PS::TimeProfile prof_system  = system_grav.getTimeProfile();
    PS::TimeProfile prof_domain  = dinfo.getTimeProfile();
    PS::F64 wtime_ex_ptcl_max, wtime_ex_ptcl_min;
    PS::S32 rank_ex_ptcl_max, rank_ex_ptcl_min;
    AllreduceMinLoc(prof_system.exchange_particle, my_rank, wtime_ex_ptcl_min, rank_ex_ptcl_min);
    AllreduceMaxLoc(prof_system.exchange_particle, my_rank, wtime_ex_ptcl_max, rank_ex_ptcl_max);

    DataForOutput data_for_output;
    data_for_output.rank = my_rank;
    data_for_output.time_sys = -1.0;
    data_for_output.n_loc = n_loc;
    data_for_output.n_ipg = 0;
    data_for_output.n_epi_ave_loc = 0;
    data_for_output.n_epj_ave_loc = 0;
    data_for_output.n_spj_ave_loc = 0;
    
    if(PS::Comm::getRank() == rank_ex_ptcl_min){
        fout_ex_ptcl_min.open(s_ex_ptcl_min.c_str(), std::fstream::app);
        data_for_output.dump(fout_ex_ptcl_min);
        DumpProfileParticleSystem(prof_system, fout_ex_ptcl_min);
        //DumpProfileTree(prof_tree, fout_ex_ptcl_min);
        fout_ex_ptcl_min<<std::endl;
        fout_ex_ptcl_min.close();
    }
    if(PS::Comm::getRank() == rank_ex_ptcl_max){
        fout_ex_ptcl_max.open(s_ex_ptcl_max.c_str(), std::fstream::app);
        data_for_output.dump(fout_ex_ptcl_max);
        DumpProfileParticleSystem(prof_system, fout_ex_ptcl_max);
        //DumpProfileTree(prof_tree, fout_ex_ptcl_max);
        fout_ex_ptcl_max<<std::endl;
        fout_ex_ptcl_max.close();
    }
#endif

    
#ifdef INITIAL_LOAD_BALANCE
    // dummy loop for load balance.
    cpe_pars.dt = 0.0;
    for(PS::S32 n=0; n<20; n++){
    //for(PS::S32 n=0; n<1; n++){
        PS::F64 wtime_beg = PS::GetWtime();
        PS::wtime_calc_force = 0.0;
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

        PS::F64 wtime_calc_force = PS::GetWtime() - wtime_beg;
        PS::F64 wtime_calc_force_min = PS::Comm::getMinValue(PS::wtime_calc_force);
        PS::F64 wtime_calc_force_max = PS::Comm::getMaxValue(PS::wtime_calc_force);
        if (PS::Comm::getRank() == 0) {
            std::cout << "PS::wtime_calc_force = " << PS::wtime_calc_force << std::endl;
            std::cout << "wtime_calc_force_min = " << wtime_calc_force_min << std::endl;
            std::cout << "wtime_calc_force_max = " << wtime_calc_force_max << std::endl;
        }

        n_loc = system_grav.getNumberOfParticleLocal();
        n_glb = system_grav.getNumberOfParticleGlobal();
        /*
        for(PS::S32 i=0; i<n_loc; i++){
            system_grav[i].pos = (((EPI*)(EPI_POINTER))+i)->pos;
            system_grav[i].vel = (((EPI*)(EPI_POINTER))+i)->vel;
        }
        */
        
        //if (n <= 3) {
        ////dinfo.collectSampleParticle(system_grav, true, (PS::F64)(PS::Comm::getRank()+1)); // std::bad_alloc
        //static PS::MTTS mt;
        //static bool first = true;
        //if (first) {
        //    mt.init_genrand( PS::Comm::getRank() );
        //    first = false;
        //}
        //PS::F64 wgh = mt.genrand_res53();
        //dinfo.collectSampleParticle(system_grav, true, wgh);
        //} else {
#ifdef DEBUG_PRINT_MAIN
        PS::Comm::barrier();if(PS::Comm::getRank()==0)std::cerr<<"OK A @main"<<std::endl;
#endif
        dinfo.collectSampleParticle(system_grav, true, PS::wtime_calc_force);
#ifdef DEBUG_PRINT_MAIN
        PS::Comm::barrier();if(PS::Comm::getRank()==0)std::cerr<<"OK B @main"<<std::endl;
#endif
        //}
        //dinfo.decomposeDomainMultiStep2(false, false);
        dinfo.decomposeDomainMultiStep3(false, false);
#ifdef DEBUG_PRINT_MAIN
        PS::Comm::barrier();if(PS::Comm::getRank()==0)std::cerr<<"OK C @main"<<std::endl;
#endif
        //system_grav.exchangeParticle2(dinfo);
        //system_grav.exchangeParticle4(dinfo);
        //system_grav.exchangeParticle5(dinfo);
        //system_grav.exchangeParticle6(dinfo);
        //system_grav.exchangeParticle7(dinfo);
        system_grav.exchangeParticle8(dinfo);
        
#ifdef DEBUG_PRINT_MAIN
        PS::Comm::barrier();if(PS::Comm::getRank()==0)std::cerr<<"OK D @main"<<std::endl;
#endif
        
        //checkLoadBalanceOfNumberOfParticle(system_grav);
        exit(1);
    }
    
    cpe_pars.dt = DT_TREE;
    //* For test run
    if (PS::Comm::getRank() == 0) fout_debug.close();
    athread_halt();
    PS::Finalize();
    return 0;
#endif

    /*
#ifdef DEBUG_CALC_FORCE
    PS::S32 id_target = 500;
    for(PS::S32 i=0; i<n_loc; i++){
        if(system_grav[i].id == id_target){
            std::cerr<<"system_grav[i].id= "<<system_grav[i].id
                     <<" system_grav[i].pos= "<<system_grav[i].pos<<std::endl;
        }
    }
#endif
    */
#ifdef DEBUG_PRINT_MAIN    
    PS::Comm::barrier();if(PS::Comm::getRank()==0)std::cerr<<"OK8 @main"<<std::endl;
#endif
#if 0
    for(PS::S32 i=0; i<n_loc; i++){
        if(system_grav[i].id == 10000){
            std::cerr<<"A) system_grav[i].id= "<<system_grav[i].id<<" pos= "<<system_grav[i].pos<<std::endl;
        }
        #ifdef PHI_R_TREE
        PS::F64 r_now = sqrt(system_grav[i].pos.y*system_grav[i].pos.y + system_grav[i].pos.z*system_grav[i].pos.z);
        #else
        PS::F64 r_now = sqrt(system_grav[i].pos*system_grav[i].pos);
        #endif
        if( (fabs(r_now - system_grav[i].r_init)/system_grav[i].r_init) >= 0.01){
            std::cerr<<"A"<<std::endl;
            std::cerr<<"system_grav[i].id= "<<system_grav[i].id<<" pos= "<<system_grav[i].pos<<std::endl;
            std::cerr<<"r_now= "<<r_now<<" r_init= "<<system_grav[i].r_init<<std::endl;
        }
    }
#endif

#ifdef WRITE_SNP
    if(my_rank==0) WriteSnp(system_grav, snp_id, dir_name, 10000);
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
    /*
    if(PS::Comm::getRank()==0){
        std::cerr<<"DIRECT FORCE"<<std::endl;
        for(PS::S32 i=0; i<N_SATELLITE; i++){
            std::cerr<<"i= "<<i<<" id= "<<SATELLITE[i].id
                     <<" fi[i].acc= "<<fi[i+n_loc].acc
                     <<std::endl;
        }
        for(PS::S32 i=0; i<n_loc; i++){
            std::cerr<<"i= "<<i<<" id= "<<system_grav[i].id
                     <<" fi[i].acc= "<<fi[i].acc
                     <<std::endl;
        }
    }
    */

    cpe_pars.energy_flag = 0;
    cpe_pars.first_flag = 0; // Full-kick
    cpe_pars.last_flag = 0;  // Full-drift
    tree_grav.calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe
        (DispatchKernelWithSP<EPI, EPJ, SPJ>,
         RetrieveKernel<FORCE>, 0,
         system_grav, dinfo, N_WALK_LIMIT,
         false,  false);
    n_loc = system_grav.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_loc; i++){
        system_grav[i].pos = (((EPI*)(EPI_POINTER))+i)->pos;
        system_grav[i].vel = (((EPI*)(EPI_POINTER))+i)->vel;
        system_grav[i].mass = (((EPI*)(EPI_POINTER))+i)->mass; // by M.I.
        system_grav[i].id = (((EPI*)(EPI_POINTER))+i)->id; // by M.I.
    }

    for(PS::S32 i=0; i<n_loc; i++){
        id_adr[i] = std::make_pair(system_grav[i].id, i);
    }
    std::sort(id_adr.begin(), id_adr.end());
    SatelliteForce * fi_tr_srt = new SatelliteForce[n_loc];
    for(PS::S32 i=0; i<n_loc; i++){
        const PS::S32 adr = id_adr[i].second;
        fi_tr_srt[i].acc = system_grav[adr].vel;
    }
    
    /*
    if(PS::Comm::getRank()==0){
        std::cerr<<"TREE FORCE"<<std::endl;
        for(PS::S32 i=0; i<N_SATELLITE; i++){
            std::cerr<<"i= "<<i<<" id= "<<SATELLITE[i].id
                     <<" fi[i].acc= "<<SATELLITE[i].vel
                     <<std::endl;
        }
        for(PS::S32 i=0; i<n_loc; i++){
            std::cerr<<"i= "<<i<<" id= "<<system_grav[i].id
                     <<" acc= "<<system_grav[i].vel
                     <<std::endl;
        }
    }
    */
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

#ifdef INTEG_CHECK
    n_loc = system_grav.getNumberOfParticleLocal();
    SatelliteForce * fi_di_org = new SatelliteForce[n_loc+N_SATELLITE]; // NOTE: not only for satellite force
    CalcForceDirect(system_grav, SATELLITE.getPointer(), N_SATELLITE, Epj::r_coll, cpe_pars.kappa, fi_di_org);
    
    std::vector< std::pair<PS::S64, PS::S32> > id_adr(n_loc);
    for(PS::S32 i=0; i<n_loc; i++){
        id_adr[i] = std::make_pair(system_grav[i].id, i);
    }
    std::sort(id_adr.begin(), id_adr.end());
    PS::F64vec * vel_di_srt = new PS::F64vec[n_loc+N_SATELLITE];
    PS::F64vec * pos_di_srt = new PS::F64vec[n_loc+N_SATELLITE];
    for(PS::S32 i=0; i<n_loc; i++){
        const PS::S32 adr = id_adr[i].second;
        vel_di_srt[i] = system_grav[adr].vel + (fi_di_org[adr].acc * DT_TREE); // full kick
        pos_di_srt[i] = system_grav[adr].pos + (vel_di_srt[i] * DT_TREE); // full drift
    }
    for(PS::S32 i=0; i<N_SATELLITE; i++){
        vel_di_srt[i+n_loc] = SATELLITE[i].vel + (fi_di_org[i+n_loc].acc * DT_TREE); // full kick
        pos_di_srt[i+n_loc] = SATELLITE[i].pos + (vel_di_srt[i+n_loc] * DT_TREE); // full drift
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
    for(PS::S32 i=0; i<n_loc; i++){
        system_grav[i].pos = (((EPI*)(EPI_POINTER))+i)->pos;
        system_grav[i].vel = (((EPI*)(EPI_POINTER))+i)->vel;
        system_grav[i].mass = (((EPI*)(EPI_POINTER))+i)->mass; // by M.I.
        system_grav[i].id = (((EPI*)(EPI_POINTER))+i)->id; // by M.I.
    }

    for(PS::S32 i=0; i<n_loc; i++){
        id_adr[i] = std::make_pair(system_grav[i].id, i);
    }
    std::sort(id_adr.begin(), id_adr.end());

    if(PS::Comm::getRank()==2){
        for(PS::S32 i=0; i<N_SATELLITE; i++){
            std::cerr<<"i= "
                     <<" id= "<<SATELLITE[i].id
                     <<" pos_di= "<<pos_di_srt[i+n_loc]
                     <<" pos_tr= "<<SATELLITE[i].pos
                     <<" vel_di= "<<vel_di_srt[i+n_loc]
                     <<" vel_tr= "<<SATELLITE[i].vel
                     <<" pos_err= "<<SATELLITE[i].pos - pos_di_srt[i+n_loc]
                     <<" vel_err= "<<SATELLITE[i].vel - vel_di_srt[i+n_loc]
                     <<std::endl;
        }
        /*
        for(PS::S32 i=0; i<n_loc; i++){
            const PS::S32 adr = id_adr[i].second;
            std::cerr<<"i= "
                     <<" id= "<<id_adr[i].first
                     <<" pos_di= "<<pos_di_srt[i]
                     <<" pos_tr= "<<system_grav[adr].pos
                     <<" vel_di= "<<vel_di_srt[i]
                     <<" vel_tr= "<<system_grav[adr].vel
                     <<" pos_err= "<<system_grav[adr].pos - pos_di_srt[i]
                     <<" vel_err= "<<system_grav[adr].vel - vel_di_srt[i]
                     <<std::endl;
        }
        */
    }
    
    PS::Finalize();
    return 0;
#endif
    
#ifdef ENERGY_CHECK
     PS::F64 eng_tot_init, eng_kin_init, eng_pot_grav_init, eng_pot_sprg_init; 
     CalcEnergy(system_grav, SATELLITE.getPointer(), N_SATELLITE, Epj::r_coll, cpe_pars.kappa,
                eng_tot_init, eng_kin_init, eng_pot_grav_init, eng_pot_sprg_init);
     if(PS::Comm::getRank()==0){
         std::cerr<<"eng_tot_init= "<<eng_tot_init
                  <<"eng_tot_init(tmp)= "<<eng_kin_init+eng_pot_grav_init
                  <<" eng_kin_init= "<<eng_kin_init
                  <<" eng_pot_grav_init= "<<eng_pot_grav_init
                  <<" eng_pot_sprg_init= "<<eng_pot_sprg_init
                  <<std::endl;
     }

     //if(my_rank==0) WriteSnp(system_grav, snp_id, dir_name, 10000);
     PS::Comm::barrier();
     
#if 0
     cpe_pars.energy_flag = 0;
     cpe_pars.first_flag = 0; // Full-kick
     cpe_pars.last_flag = 0;  // Full-drift
     tree_grav.calcForceAllAndWriteBackReuseListMultiWalkIndexUnsafe
         (DispatchKernelWithSP<EPI, EPJ, SPJ>,
          RetrieveKernel<FORCE>, 0,
          system_grav, dinfo, N_WALK_LIMIT,
          false,  false);
     n_loc = system_grav.getNumberOfParticleLocal();
     for(PS::S32 i=0; i<n_loc; i++){
         system_grav[i].pos = (((EPI*)(EPI_POINTER))+i)->pos;
         system_grav[i].vel = (((EPI*)(EPI_POINTER))+i)->vel;
         system_grav[i].mass = (((EPI*)(EPI_POINTER))+i)->mass; // by M.I.
         system_grav[i].id = (((EPI*)(EPI_POINTER))+i)->id; // by M.I.
     }
#else
    n_loc = system_grav.getNumberOfParticleLocal();
    SatelliteForce * fi_di_org = new SatelliteForce[n_loc+N_SATELLITE]; // NOTE: not only for satellite force
    CalcForceDirect(system_grav, SATELLITE.getPointer(), N_SATELLITE, Epj::r_coll, cpe_pars.kappa, fi_di_org);
    
    std::vector< std::pair<PS::S64, PS::S32> > id_adr(n_loc);
    for(PS::S32 i=0; i<n_loc; i++){
        id_adr[i] = std::make_pair(system_grav[i].id, i);
    }
    std::sort(id_adr.begin(), id_adr.end());
    PS::F64vec * vel_di_srt = new PS::F64vec[n_loc+N_SATELLITE];
    PS::F64vec * pos_di_srt = new PS::F64vec[n_loc+N_SATELLITE];
    for(PS::S32 i=0; i<n_loc; i++){
        const PS::S32 adr = id_adr[i].second;
        vel_di_srt[i] = system_grav[adr].vel + (fi_di_org[adr].acc * DT_TREE); // full kick
        //pos_di_srt[i] = system_grav[adr].pos + (vel_di_srt[i] * DT_TREE); // full drift
        pos_di_srt[i] = system_grav[adr].pos + 0.5*(vel_di_srt[i]+system_grav[adr].vel)*DT_TREE; // full drift
    }
    for(PS::S32 i=0; i<n_loc; i++){
        const PS::S32 adr = id_adr[i].second;
        system_grav[adr].vel = vel_di_srt[i];
        system_grav[adr].pos = pos_di_srt[i];
    }

    for(PS::S32 i=0; i<N_SATELLITE; i++){
        vel_di_srt[i+n_loc] = SATELLITE[i].vel + (fi_di_org[i+n_loc].acc * DT_TREE); // full kick
        pos_di_srt[i+n_loc] = SATELLITE[i].pos + 0.5*(vel_di_srt[i+n_loc]+SATELLITE[i].vel)*DT_TREE; // full drift
    }
    for(PS::S32 i=0; i<N_SATELLITE; i++){
        SATELLITE[i].vel = vel_di_srt[i+n_loc];
        SATELLITE[i].pos = pos_di_srt[i+n_loc];
    }    
    
#endif
     
     PS::F64 eng_tot_now, eng_kin_now, eng_pot_grav_now, eng_pot_sprg_now; 
     CalcEnergy(system_grav, SATELLITE.getPointer(), N_SATELLITE, Epj::r_coll, cpe_pars.kappa,
                eng_tot_now, eng_kin_now, eng_pot_grav_now, eng_pot_sprg_now);
     if(PS::Comm::getRank()==0){
         std::cerr<<"eng_tot_now= "<<eng_tot_now
                  <<"eng_tot_now(tmp)= "<<eng_kin_now+eng_pot_grav_now
                  <<" eng_kin_now= "<<eng_kin_now
                  <<" eng_pot_grav_now= "<<eng_pot_grav_now
                  <<" eng_pot_sprg_now= "<<eng_pot_sprg_now
                  <<std::endl;
     }

     //if(my_rank==0) WriteSnp(system_grav, snp_id, dir_name, 10000);
     PS::Comm::barrier();
     
     PS::Finalize();
     return 0;
     
#endif //ENERGY_CHECK

#ifdef ENERGY_CHECK_2
     PS::F64 eng_tot_init, eng_kin_init, eng_pot_grav_init, eng_pot_sprg_init; 
     CalcEnergy(system_grav, SATELLITE.getPointer(), N_SATELLITE, Epj::r_coll, cpe_pars.kappa,
                eng_tot_init, eng_kin_init, eng_pot_grav_init, eng_pot_sprg_init);
     if(PS::Comm::getRank()==0){
         std::cerr<<"eng_tot_init= "<<eng_tot_init
                  <<"eng_tot_init(tmp)= "<<eng_kin_init+eng_pot_grav_init
                  <<" eng_kin_init= "<<eng_kin_init
                  <<" eng_pot_grav_init= "<<eng_pot_grav_init
                  <<" eng_pot_sprg_init= "<<eng_pot_sprg_init
                  <<std::endl;
     }
     PS::Comm::barrier();
#endif //ENERGY_CHECK_2
     
    
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
     nan_test((EPI*)EPI_POINTER, system_grav.getNumberOfParticleLocal(), "nan-0");



     
#ifdef ENERGY_CHECK
     PS::F64 e_tot_glb_init, e_kin_glb_init, e_pot_glb_init; 
     CalcEnergy(system_grav, SATELLITE, e_tot_glb_init, e_kin_glb_init, e_pot_glb_init);
     if(PS::Comm::getRank()==0){
         std::cerr<<"e_tot_glb_init= "<<e_tot_glb_init
                  <<"e_kin_glb_init= "<<e_kin_glb_init
                  <<"e_pot_glb_init= "<<e_pot_glb_init
                  <<std::endl;
     }
#endif //ENERGY_CHECK
     
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

     /*     
#ifdef PHI_R_TREE
     #ifdef USE_MY_ATAN2
     double wtime_atan2 = PS::GetWtime();
     args[0] = (unsigned long) n_loc;
     args[1] = (unsigned long) (EPI*)EPI_POINTER;
     args[2] = (unsigned long) &system_grav[0];
     __real_athread_spawn((void*)slave_CopyEpiToFpWithCoordinateTransformCPE, args);
     athread_join();
     wtime_atan2 = PS::GetWtime() - wtime_atan2;
     if(PS::Comm::getRank() == 0){
         std::cerr<<"wtime_atan2= "<<wtime_atan2<<std::endl;
     }
     for(PS::S32 i=0; i<n_loc; i++){
         double x_tmp = ((EPI*)EPI_POINTER)[i].pos.x;
         double y_tmp = ((EPI*)EPI_POINTER)[i].pos.y;
         double phi_tmp = atan2(y_tmp, x_tmp);
         if(phi_tmp < 0.0) phi_tmp += 2.0*4.0*atan(1.0);
         else if(phi_tmp >= 2.0*4.0*atan(1.0)) phi_tmp -= 2.0*4.0*atan(1.0);
         assert( fabs(system_grav[i].pos.x - phi_tmp) < 1e-14);
     }
     #else
     double wtime_atan2 = PS::GetWtime();
     CopyEpiToFpWithCoordinateTransform(&system_grav[0], (EPI*)EPI_POINTER, n_loc);
     wtime_atan2 = PS::GetWtime() - wtime_atan2;
     if(PS::Comm::getRank() == 0){
         std::cerr<<"wtime_atan2= "<<wtime_atan2<<std::endl;
     }
     exit(1);
     #endif
     for(PS::S32 i=0; i<n_loc; i++){
         if(system_grav[i].pos.x<0.0 || system_grav[i].pos.x>=2.0*PI){
             std::cerr<<"system_grav[i].pos= "<<system_grav[i].pos<<std::endl;
         }
         assert(system_grav[i].pos.x >= 0.0 && system_grav[i].pos.x < 2.0*PI);
     }
#else
     for(PS::S32 i=0; i<n_loc; i++){
         system_grav[i].pos = (((EPI*)(EPI_POINTER))+i)->pos;
         system_grav[i].vel = (((EPI*)(EPI_POINTER))+i)->vel;
     }
#endif
     */

     for(PS::S32 i=0; i<n_loc; i++){
         system_grav[i].pos = (((EPI*)(EPI_POINTER))+i)->pos;
         system_grav[i].vel = (((EPI*)(EPI_POINTER))+i)->vel;
         system_grav[i].mass = (((EPI*)(EPI_POINTER))+i)->mass; // by M.I.
         system_grav[i].id = (((EPI*)(EPI_POINTER))+i)->id; // by M.I.
     }

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
     
     
#if 0
    for(PS::S32 i=0; i<n_loc; i++){
        if(system_grav[i].id == 10000){
            std::cerr<<"B) system_grav[i].id= "<<system_grav[i].id<<" pos= "<<system_grav[i].pos<<std::endl;
        }
        #ifdef PHI_R_TREE
        PS::F64 r_now = sqrt(system_grav[i].pos.y*system_grav[i].pos.y + system_grav[i].pos.z*system_grav[i].pos.z);
        #else
        PS::F64 r_now = sqrt(system_grav[i].pos*system_grav[i].pos);
        #endif
        if( (fabs(r_now - system_grav[i].r_init)/system_grav[i].r_init) >= 0.01){
            std::cerr<<"B"<<std::endl;
            std::cerr<<"system_grav[i].id= "<<system_grav[i].id<<" pos= "<<system_grav[i].pos<<std::endl;
            std::cerr<<"r_now= "<<r_now<<" r_init= "<<system_grav[i].r_init<<std::endl;
        }        
    }
#endif
     PS::Comm::barrier();if(PS::Comm::getRank()==0)std::cerr<<"OK10 @main"<<std::endl;
     
#ifdef DEBUG_CALC_FORCE
    if(my_rank == 0) std::cerr<<"CHECK CALC FORCE"<<std::endl;
    MPI_Allgather(&n_loc, 1, MPI_INT, n_loc_array, 1, MPI_INT, MPI_COMM_WORLD);
    PS::S32 * n_disp_loc_array = new PS::S32[n_proc+1];
    n_disp_loc_array[0] = 0;
    for(PS::S32 i=0; i<n_proc; i++){
        n_disp_loc_array[i+1] = n_loc_array[i] + n_disp_loc_array[i];
    }
    assert(n_glb == n_disp_loc_array[n_proc]);
    FPGrav * ptcl_glb = new FPGrav[ n_glb ];
    MPI_Allgatherv(&system_grav[0], n_loc, PS::GetDataType<FPGrav>(),
                   ptcl_glb, n_loc_array, n_disp_loc_array, PS::GetDataType<FPGrav>(),
                   MPI_COMM_WORLD);
    /*
    for(PS::S32 i=0; i<n_loc; i++){
        if(ptcl_glb[i].id == id_target){
            PS::F64vec rij_sun = (((EPI*)EPI_POINTER)+i)->pos;
            PS::F64 inv_r_sun = 1.0 / sqrt(rij_sun*rij_sun);
            PS::F64vec acc_sun = - 1.0 * inv_r_sun * inv_r_sun * inv_r_sun * rij_sun;
            std::cerr<<"ptcl_glb[i].id= "<<ptcl_glb[i].id
                     <<" ptcl_glb[i].pos= "<<ptcl_glb[i].pos
                     <<std::endl;
            std::cerr<<"system_grav[i].id= "<<system_grav[i].id
                //<<" FORCE_ARRAY[i].acc= "<<FORCE_ARRAY[i].acc
                //<<" acc(net)= "<<FORCE_ARRAY[i].acc - acc_sun
                //<<" rij_sun= "<<rij_sun
                     <<" system_grav[i].pos= "<<system_grav[i].pos<<std::endl;
        }
    }
    */
    
    #ifdef PHI_R_TREE
    for(PS::S32 i=0; i<n_glb; i++){
        CoordinateChangeCylindricalToCartesian(&ptcl_glb[i]);
    }
    for(PS::S32 i=0; i<n_loc; i++){
        CoordinateChangeCylindricalToCartesian(&system_grav[i]);
    }
    #endif
    PS::F64vec * acc_loc = new PS::F64vec[n_loc];
    for(PS::S32 i=0; i<n_loc; i++){
        acc_loc[i] = 0.0;
        PS::F64 eps_sq = EPI::eps * EPI::eps;
        for(PS::S32 j=0; j<n_glb; j++){
            if(ptcl_glb[j].id == system_grav[i].id) continue;
            PS::F64vec rij = system_grav[i].pos - ptcl_glb[j].pos;
            //PS::F64 r_inv = 1.0 / sqrt(rij*rij);
            PS::F64 r_inv = 1.0 / sqrt(rij*rij+eps_sq);
            PS::F64 m_r = ptcl_glb[j].mass * r_inv;
            acc_loc[i] -= m_r * r_inv * r_inv * rij;
        }
    }
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64vec r_pla = system_grav[i].pos;
        PS::F64 r_inv = 1.0 / sqrt(r_pla*r_pla);
        PS::F64 m_r = PLANET.mass * r_inv;
        PS::F64vec acc_pla = -m_r * r_inv * r_inv * r_pla;
        PS::F64vec da = FORCE_ARRAY[i].acc - acc_loc[i];
        if( ( (da*da)/(acc_loc[i]*acc_loc[i]) ) >= 1.0){
            std::cerr<<"system_grav[i].id= "<<system_grav[i].id
                     <<" system_grav[i].pos= "<<system_grav[i].pos
                     <<" FORCE_ARRAY[i].acc= "<<FORCE_ARRAY[i].acc
                     <<" acc_loc[i]= "<<acc_loc[i]
                     <<std::endl;
        }
    }
    PS::S64 n_interaction_ep_ep_loc = 0.0;
    PS::S64 n_interaction_ep_sp_loc = 0.0;
    PS::S64 n_ipg = tree_grav.getNumberOfIPG();
    for(PS::S32 i=0; i<n_ipg; i++){
        n_interaction_ep_ep_loc += (PS::S64)(tree_grav.n_epi_recorder_for_force_[0][i])*(PS::S64)(tree_grav.n_epj_recorder_for_force_[0][i]);
        n_interaction_ep_sp_loc += (PS::S64)(tree_grav.n_epi_recorder_for_force_[0][i])*(PS::S64)(tree_grav.n_spj_recorder_for_force_[0][i]);
    }
    PS::S64 n_interaction_ep_ep_glb = PS::Comm::getSum(n_interaction_ep_ep_loc);
    PS::S64 n_interaction_ep_sp_glb = PS::Comm::getSum(n_interaction_ep_sp_loc);
    PS::F64 da_loc = 0.0;
    PS::F64 * acc_err_loc = new PS::F64[n_loc];
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64vec da = FORCE_ARRAY[i].acc - acc_loc[i];
        acc_err_loc[i] = sqrt( (da*da) / (acc_loc[i]*acc_loc[i]) );
        da_loc += (da*da) / (acc_loc[i]*acc_loc[i]);
    }
    PS::F64 da_glb = PS::Comm::getSum(da_loc);
    if(PS::Comm::getRank() == 0){
        std::cerr<<"theta= "<<theta<<" da= "<<sqrt(da_glb / n_glb)<<" n_glb= "<<n_glb
                 <<" dinfo.getNDomain(0)= "<<dinfo.getNDomain(0)<<" dinfo.getNDomain(1)= "<<dinfo.getNDomain(1)
                 <<" n_ep_ave= "<<(PS::F64)n_interaction_ep_ep_glb/n_glb
                 <<" n_sp_ave= "<<(PS::F64)n_interaction_ep_sp_glb/n_glb
                 <<std::endl;
    }
    PS::F64 * acc_err_glb = new PS::F64[n_glb];
    MPI_Allgatherv(acc_err_loc, n_loc, PS::GetDataType<PS::F64>(),
                   acc_err_glb, n_loc_array, n_disp_loc_array, PS::GetDataType<PS::F64>(),
                   MPI_COMM_WORLD);
    std::sort(acc_err_glb, acc_err_glb+n_glb);
    if(PS::Comm::getRank() == 0){
        for(PS::S32 i=0; i<n_glb; i++){
            fout_debug<<acc_err_glb[i]<<"   "<<(PS::F64)(i+1)/n_glb<<std::endl;
        }
        fout_debug.close();
    }
    exit(1);
#endif
    //exit(1);

    
    //system_grav.dumpMemSizeUsed(fout_debug);
    //tree_grav.dumpMemSizeUsed(fout_debug);

#if defined(ROTATE_PTCLS) && defined(SHIFT_CENTER_Z)
    #ifdef DEBUG_SHIFT_CENTER_Z
    z_coord_offset_loc = 0.0;
    for(PS::S32 i=0; i<n_loc; i++){
        if( z_coord_offset_loc < fabs(system_grav[i].pos.z) ) z_coord_offset_loc = fabs(system_grav[i].pos.z);
    }
    PS::TREE_Z_COORD_OFFSET = PS::Comm::getMaxValue(z_coord_offset_loc);
    if (PS::Comm::getRank() == 0)
        std::cerr<<"PS::TREE_Z_COORD_OFFSET(0)= "<<PS::TREE_Z_COORD_OFFSET<<std::endl;
    #endif
     
    starttime = PS::GetWtime();
    theta_rot = - DT_TREE;
    const PS::F64 cth = std::cos(theta_rot);
    const PS::F64 sth = std::sin(theta_rot);
    PS::F64 z_shift_ptcl_ar[64];
    PS::F64 z_shift_sat_ar[64];
    //unsigned long args[4];
    //** Rotate the particles
    args[0] = (unsigned long) system_grav.getNumberOfParticleLocal();
    args[1] = (unsigned long) &system_grav[0];
    args[2] = (unsigned long) &cth;
    args[3] = (unsigned long) &sth;
    args[4] = (unsigned long) z_shift_ptcl_ar;
    __real_athread_spawn((void*)slave_RotateAndShiftZ, args);
    athread_join();
    //** Rotate the satellites
    args[0] = (unsigned long) N_SATELLITE;
    args[1] = (unsigned long) SATELLITE.getPointer();
    args[2] = (unsigned long) &cth;
    args[3] = (unsigned long) &sth;
    args[4] = (unsigned long) z_shift_sat_ar;
    __real_athread_spawn((void*)slave_RotateAndShiftZ, args);
    athread_join();
    z_coord_offset_loc = 0.0;
    for(PS::S32 i=0; i<64; i++){
        if(fabs(z_shift_ptcl_ar[i]) > fabs(z_coord_offset_loc) ){
            z_coord_offset_loc = z_shift_ptcl_ar[i];
        }
        if(fabs(z_shift_sat_ar[i]) > fabs(z_coord_offset_loc) ){
            z_coord_offset_loc = z_shift_sat_ar[i];
        }
    }
    //z_coord_offset_loc = (fabs(z_shift_ptcl_loc) > fabs(z_shift_sat_loc)) ? z_shift_ptcl_loc : z_shift_sat_loc;
    PS::F64 z_shift_max_tmp = PS::Comm::getMaxValue(z_coord_offset_loc);
    PS::F64 z_shift_min_tmp = PS::Comm::getMinValue(z_coord_offset_loc);
    PS::TREE_Z_COORD_OFFSET = ( fabs(z_shift_max_tmp) > fabs(z_shift_min_tmp) ) ? z_shift_max_tmp : z_shift_min_tmp;
    #ifdef DEBUG_SHIFT_CENTER_Z
    if (PS::Comm::getRank() == 0)
        std::cerr<<"PS::TREE_Z_COORD_OFFSET(1)= "<<PS::TREE_Z_COORD_OFFSET<<std::endl;
    #endif
    if(PS::TREE_Z_COORD_OFFSET > 0.0){PS::TREE_Z_COORD_OFFSET +=  Epj::r_coll*0.1;}
    else{PS::TREE_Z_COORD_OFFSET -=  Epj::r_coll*0.1;}
    endtime = PS::GetWtime();
    if (PS::Comm::getRank() == 0)
        fout_debug << "wtime_calc_shift_z = " << endtime - starttime << std::endl;     
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
         fout_debug << "wtime_rotate = " << endtime - starttime << std::endl;
    #else // SUNWAY
     theta_rot = - DT_TREE;
     Rotate(system_grav,n_loc,theta_rot);
     Rotate(SATELLITE,N_SATELLITE,theta_rot);
     endtime = PS::GetWtime();
     if (PS::Comm::getRank() == 0)
         fout_debug << "wtime_rotate = " << endtime - starttime << std::endl;
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
         fout_debug << "wtime_calc_shift_z = " << endtime - starttime << std::endl;
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
     
    //fout_debug.close();
    //athread_halt();
    //PS::Finalize();
    //return 0;

    //* Compute the total anugular momentum of the system at t=dt.
    ANGULAR_MOMENTUM_TEST0(system_grav);  

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
             fout_debug << "#######################" << std::endl;
             fout_debug << "  n_loop = " << n_loop << std::endl;
             fout_debug << "#######################" << std::endl;
        }
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
            
        //system_grav.exchangeParticle2(dinfo);
        //system_grav.exchangeParticle4(dinfo);
        //system_grav.exchangeParticle5(dinfo);
        //system_grav.exchangeParticle6(dinfo);
        //system_grav.exchangeParticle7(dinfo);
        
#ifdef FAST_EX_PTCL
        if(n_loop <= 3){
            system_grav.exchangeParticle8(dinfo, false);
        }
        else{
            system_grav.exchangeParticle8(dinfo, true);
        }
#else
        dinfo.decomposeDomainMultiStep3(true, false); // 1st: sample sort, 2nd split y-coord
        system_grav.exchangeParticle8(dinfo, false);
#endif
        n_loc = system_grav.getNumberOfParticleLocal();
#ifdef DEBUG_PRINT_MAIN
        std::cerr<<"rank= "<<my_rank<<" n_loc= "<<n_loc<<std::endl;
#endif

        /*
    #ifdef SHIFT_CENTER_Z
        starttime = PS::GetWtime();
        z_coord_offset_loc = 0.0;
        for(PS::S32 i=0; i<n_loc; i++){
            if( z_coord_offset_loc < fabs(system_grav[i].pos.z) ) z_coord_offset_loc = fabs(system_grav[i].pos.z);
        }
        PS::TREE_Z_COORD_OFFSET = PS::Comm::getMaxValue(z_coord_offset_loc) + 1.0e-6;
        endtime = PS::GetWtime();
        if (PS::Comm::getRank() == 0)
            fout_debug << "wtime_calc_shift_z = " << endtime - starttime << std::endl;
    #endif
        */
        
#ifdef DEBUG_PRINT_MAIN
        CalcCMCyl(&system_grav[0], n_loc, MASS_DUST_TOTAL, CM_POS_DUST_TOTAL);
        if(PS::Comm::getRank()==0){
            std::cerr<<" (FP.2) MASS_DUST_TOTAL= "<<MASS_DUST_TOTAL
                     <<" CM_POS_DUST_TOTAL= "<<CM_POS_DUST_TOTAL
                     <<std::endl;         
        }
#endif
        
#ifdef DEBUG_PRINT_MAIN
        //PS::Comm::barrier();
        //if(PS::Comm::getRank()==0)std::cerr<<"OK13-"<<n_loop<<" @main"<<std::endl;
#endif
        NAN_TEST(system_grav,"[nan-c]")

#if 0
        checkLoadBalanceOfNumberOfParticle(system_grav);
#endif
        

#if 0        
        ptcl_old.resizeNoInitialize(n_loc);
        for(PS::S32 i=0; i<n_loc; i++){
            ptcl_old[i] = system_grav[i];
        }
        std::ofstream fout_before;
        std::string s_before = std::string(dir_name) + ("/before_" + ToS(PS::Comm::getRank()) + ".dat");
        fout_before.open(s_before.c_str());
        for(PS::S32 i=0; i<n_loc; i++){
            fout_before<<ptcl_old[i].pos<<std::endl;
        }
        fout_before.close();
#endif
        
#ifdef DEBUG_PRINT_MAIN
        PS::Comm::barrier(); if(PS::Comm::getRank()==0)std::cerr<<"OK14-"<<n_loop<<" @main"<<std::endl;
#endif
#if 0        
        for(PS::S32 i=0; i<n_loc; i++){
            if(system_grav[i].id == 10000){
                std::cerr<<"C) system_grav[i].id= "<<system_grav[i].id<<" pos= "<<system_grav[i].pos<<std::endl;
            }
     #ifdef PHI_R_TREE
            PS::F64 r_now = sqrt(system_grav[i].pos.y*system_grav[i].pos.y + system_grav[i].pos.z*system_grav[i].pos.z);
     #else
            PS::F64 r_now = sqrt(system_grav[i].pos*system_grav[i].pos);
     #endif
            if( (fabs(r_now - system_grav[i].r_init)/system_grav[i].r_init) >= 0.01){
                std::cerr<<"C"<<std::endl;
                std::cerr<<"system_grav[i].id= "<<system_grav[i].id<<" pos= "<<system_grav[i].pos<<std::endl;
                std::cerr<<"r_now= "<<r_now<<" r_init= "<<system_grav[i].r_init<<std::endl;
            }
            system_grav[i].r_init = r_now;
        }
#endif
        
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

        //PS::Comm::barrier();
        //exit(1);

#ifdef ENERGY_CHECK_2
     PS::F64 eng_tot_now, eng_kin_now, eng_pot_grav_now, eng_pot_sprg_now; 
     CalcEnergy(system_grav, SATELLITE.getPointer(), N_SATELLITE, Epj::r_coll, cpe_pars.kappa,
                eng_tot_now, eng_kin_now, eng_pot_grav_now, eng_pot_sprg_now);
     if(PS::Comm::getRank()==0){
         std::cerr<<"eng_tot_now= "<<eng_tot_now
                  <<"eng_tot_now(tmp)= "<<eng_kin_now+eng_pot_grav_now
                  <<" eng_kin_now= "<<eng_kin_now
                  <<" eng_pot_grav_now= "<<eng_pot_grav_now
                  <<" eng_pot_sprg_now= "<<eng_pot_sprg_now
                  <<std::endl;
     }
     PS::Comm::barrier();
#endif //ENERGY_CHECK
        
        //exit(1);
            
        //system_grav.dumpMemSizeUsed(fout_debug);
        //tree_grav.dumpMemSizeUsed(fout_debug);
        //PS::Comm::barrier();
        //fout_debug.close();
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
        endtime = PS::GetWtime();


        
        if (PS::Comm::getRank() == 0){
            for(PS::S32 i=0; i<5; i++){
                std::cerr<<"B) system_grav[i].pos= "<<system_grav[i].pos
                         <<" id= "<<system_grav[i].id
                         <<std::endl;
            }
        }
        
        
        if (PS::Comm::getRank() == 0)
            fout_debug << "wtime_copy_epi_to_psys = " << endtime - starttime << std::endl;


#if defined(ROTATE_PTCLS) && defined(SHIFT_CENTER_Z)
    #ifdef DEBUG_SHIFT_CENTER_Z
    z_coord_offset_loc = 0.0;
    for(PS::S32 i=0; i<n_loc; i++){
        if( z_coord_offset_loc < fabs(system_grav[i].pos.z) ) z_coord_offset_loc = fabs(system_grav[i].pos.z);
    }
    PS::TREE_Z_COORD_OFFSET = PS::Comm::getMaxValue(z_coord_offset_loc);
    if (PS::Comm::getRank() == 0)
        std::cerr<<"PS::TREE_Z_COORD_OFFSET(0)= "<<PS::TREE_Z_COORD_OFFSET<<std::endl;
    #endif
     
    starttime = PS::GetWtime();
    theta_rot = - DT_TREE * n_loop_merge;
    const PS::F64 cth = std::cos(theta_rot);
    const PS::F64 sth = std::sin(theta_rot);
    PS::F64 z_shift_ptcl_ar[64];
    PS::F64 z_shift_sat_ar[64];
    //unsigned long args[4];
    //** Rotate the particles
    args[0] = (unsigned long) system_grav.getNumberOfParticleLocal();
    args[1] = (unsigned long) &system_grav[0];
    args[2] = (unsigned long) &cth;
    args[3] = (unsigned long) &sth;
    args[4] = (unsigned long) z_shift_ptcl_ar;
    __real_athread_spawn((void*)slave_RotateAndShiftZ, args);
    athread_join();
    //** Rotate the satellites
    args[0] = (unsigned long) N_SATELLITE;
    args[1] = (unsigned long) SATELLITE.getPointer();
    args[2] = (unsigned long) &cth;
    args[3] = (unsigned long) &sth;
    args[4] = (unsigned long) z_shift_sat_ar;
    __real_athread_spawn((void*)slave_RotateAndShiftZ, args);
    athread_join();
    z_coord_offset_loc = 0.0;
    for(PS::S32 i=0; i<64; i++){
        if(fabs(z_shift_ptcl_ar[i]) > fabs(z_coord_offset_loc) ){
            z_coord_offset_loc = z_shift_ptcl_ar[i];
        }
        if(fabs(z_shift_sat_ar[i]) > fabs(z_coord_offset_loc) ){
            z_coord_offset_loc = z_shift_sat_ar[i];
        }
    }
    //z_coord_offset_loc = (fabs(z_shift_ptcl_loc) > fabs(z_shift_sat_loc)) ? z_shift_ptcl_loc : z_shift_sat_loc;
    PS::F64 z_shift_max_tmp = PS::Comm::getMaxValue(z_coord_offset_loc);
    PS::F64 z_shift_min_tmp = PS::Comm::getMinValue(z_coord_offset_loc);
    PS::TREE_Z_COORD_OFFSET = ( fabs(z_shift_max_tmp) > fabs(z_shift_min_tmp) ) ? z_shift_max_tmp : z_shift_min_tmp;
    #ifdef DEBUG_SHIFT_CENTER_Z
    if (PS::Comm::getRank() == 0)
        std::cerr<<"PS::TREE_Z_COORD_OFFSET(1)= "<<PS::TREE_Z_COORD_OFFSET<<std::endl;
    #endif
    if(PS::TREE_Z_COORD_OFFSET > 0.0){PS::TREE_Z_COORD_OFFSET +=  Epj::r_coll*0.1;}
    else{PS::TREE_Z_COORD_OFFSET -=  Epj::r_coll*0.1;}
    endtime = PS::GetWtime();
    if (PS::Comm::getRank() == 0)
        fout_debug << "wtime_calc_shift_z = " << endtime - starttime << std::endl;     
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
            fout_debug << "wtime_rotate = " << endtime - starttime << std::endl;
        #else //SUNWAY
        theta_rot  = - DT_TREE * n_loop_merge;
        Rotate(system_grav,n_loc,theta_rot);
        Rotate(SATELLITE,N_SATELLITE,theta_rot);
        endtime = PS::GetWtime();
        if (PS::Comm::getRank() == 0)
            fout_debug << "wtime_rotate = " << endtime - starttime << std::endl;
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
            fout_debug << "wtime_calc_shift_z = " << endtime - starttime << std::endl;
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
        
#if 0
        ptcl_new.resizeNoInitialize(n_loc);
        for(PS::S32 i=0; i<n_loc; i++){
            ptcl_new[i] = system_grav[i];
        }
        std::ofstream fout_after;
        std::string s_after = std::string(dir_name) + "/after_" + ToS(PS::Comm::getRank()) + ".dat";
        fout_after.open(s_after.c_str());
        for(PS::S32 i=0; i<n_loc; i++){
            fout_after<<ptcl_new[i].pos<<std::endl;
        }
        fout_after.close();
#endif
        
        //* Update time_sys
        time_sys += n_loop_merge*DT_TREE;

        //std::cerr<<"time_sys= "<<time_sys<<std::endl;
        
        NAN_TEST(system_grav,"[nan-e]");

#if 0
        checkRingWidth(system_grav);
#endif

        ANGULAR_MOMENTUM_TEST(system_grav);

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
        /*
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
        */
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

        //* Check memory usage
        system_grav.dumpMemSizeUsed(fout_debug);
        tree_grav.dumpMemSizeUsed(fout_debug);
        //break;
    }
    
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_close();
#endif

    if (PS::Comm::getRank() == 0) fout_debug.close();


    
#ifdef SUNWAY
    athread_halt();
#endif
    PS::Finalize();

    return 0;
}

