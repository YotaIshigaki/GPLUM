// Include the standard C++ headers
#include <cmath>
#include <math.h>
#include <cfloat>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <time.h>
// Include the header file of FDPS
#include <particle_simulator.hpp>
#include <particle_mesh_multipole.hpp>
// Include the header file of Phantom-GRAPE library
#ifdef ENABLE_PHANTOM_GRAPE_X86
#include <gp5util.h>
#endif
// Include user-defined headers
#include "user_defined.hpp"

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

void setupIC(PS::ParticleSystem<FP_nbody> & psys,
             PS::DomainInfo & dinfo) {
    // Parameters
    //const int NP=128;
    //const int NP=256;
    //const int NP=512;
    const int NP=32768;
    const long seed=19810614;
    const bool flag_charge_neutral {true};

    if (PS::Comm::getRank() == 0) {
        psys.setNumberOfParticleLocal(NP);
        // Initialize pseudo-random generator
        srand48(seed);
        // Set psys except for mass
        for(PS::S32 i = 0; i < NP; i++){
            psys[i].id  = i;
            psys[i].pos = PS::F64vec(drand48(), drand48(), drand48());
            psys[i].vel = 0.0;
            psys[i].acc = 0.0;
            psys[i].pot = 0.0;
        }
        // Set mass 
        PS::F64 msum = 0.0;
        for (PS::S32 i = 0; i < NP; i++){
            msum += (psys[i].mass  = drand48() * (1.0/NP));
        }
        // Adjust mass so that the system is charge neutral.
        if (flag_charge_neutral) {
            for (PS::S32 i = 0; i < NP; i++){
                psys[i].mass -= msum / NP;
            }
        }
        // Output 
        std::string filename = "IC.txt";
        std::ofstream output_file;
        output_file.open(filename.c_str(), std::ios::trunc);
        for (PS::S32 i = 0; i < NP; i++) 
            output_file << psys[i].pos << std::endl;
        output_file.close();
    } else {
        psys.setNumberOfParticleLocal(0);
    }

    // Set domain information
    //dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0),
                           PS::F64vec(1.0, 1.0, 1.0));
}

void calcForceSum(PS::ParticleSystem<FP_nbody> & psys) {
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    PS::F64vec fpp_loc(0.0);
    PS::F64vec fpm_loc(0.0);
    for (PS::S32 i = 0; i < n_loc; i++) {
        fpp_loc += psys[i].mass * psys[i].acc;
        fpm_loc += psys[i].mass * psys[i].acc_pm;
    }
    PS::F64vec fpp_tot;
    PS::F64vec fpm_tot;
    fpp_tot = PS::Comm::getSum(fpp_loc);
    fpm_tot = PS::Comm::getSum(fpm_loc);
    if (PS::Comm::getRank() == 0) {
        std::cout << "PP ftot : " << fpp_tot << std::endl;
        std::cout << "PM ftot : " << fpm_tot << std::endl;
    }
}

int main(int argc, char *argv[]) {
    // Configure stdout & stderr
    std::cout << std::setprecision(15);
    std::cerr << std::setprecision(15);

    //* Initialize FDPS
    PS::Initialize(argc, argv);

    // Make a directory
    char dir_name[1024];
    sprintf(dir_name,"./result");
    makeOutputDirectory(dir_name);

    // Display # of MPI processes and threads
    PS::S32 nprocs = PS::Comm::getNumberOfProc();
    PS::S32 nthrds = PS::Comm::getNumberOfThread();
    if (PS::Comm::getRank() == 0) {
        std::cout << "===========================================" << std::endl
                  << " This is a sample code for PMMM with FDPS!"  << std::endl
                  << " # of processes is " << nprocs               << std::endl
                  << " # of thread is    " << nthrds               << std::endl
                  << "===========================================" << std::endl;
    }

    // Make an instance of ParticleSystem and initialize it
    PS::ParticleSystem<FP_nbody> psys;
    psys.initialize();

    // Make an instance of DomainInfo and initialize it
    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);

    // Make an initial condition and initialize the particle system
    setupIC(psys, dinfo);

    // Perform domain decomposition and particle exchange
    dinfo.decomposeDomainAll(psys);
    psys.exchangeParticle(dinfo);

    // Make an instance of TreeForForce and initialize it
    const PS::S32 n_tot = 3 * psys.getNumberOfParticleLocal();
    const PS::S32vec n_cell = PS::S32vec(8, 8, 8);
    const PS::S32 icut = 2;
    const PS::F32 theta = 0.0;
    const PS::U32 n_leaf_limit = 8;
    const PS::U32 n_group_limit = 64;
    PS::TreeForForceLong<FP_nbody, FP_nbody, FP_nbody>::MonopoleGeometricCenterParticleMeshMultipole tree;
    tree.initialize(n_tot, n_cell, icut, theta, n_leaf_limit, n_group_limit);

    // Initialize Phantom-GRAPE library
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_open();
    g5_set_eps_to_all(FP_nbody::eps);
#endif

    // Perform force calculation (PP part)
    tree.calcForceAllAndWriteBack(CalcForce<FP_nbody>,
                                  CalcForce<PS::SPJMonopoleGeometricCenterPMMM>,
                                  psys, dinfo);

    // Make an instance of ParticleMeshMultipole class and initialize it
    const PS::S32 p = 5;
    PS::PMM::ParticleMeshMultipole<FP_nbody, FP_nbody, FP_nbody, PS::SPJMonopoleGeometricCenterPMMM> pm;
    pm.initialize(p);

    // Perform force calculation (PM part)
    pm.calcForceAllAndWriteBack(tree,dinfo,psys);

    // Output the result of PM part
    std::ostringstream ss;
    ss << "PM" << std::setfill('0') << std::setw(5) << PS::Comm::getRank() << ".txt";
    std::string filename = ss.str();
    std::ofstream output_file;
    output_file.open(filename.c_str(), std::ios::trunc);
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        output_file << psys[i].pos.x << "   "
                    << psys[i].pos.y << "   "
                    << psys[i].pos.z << "   "
                    << psys[i].acc_pm.x << "   "
                    << psys[i].acc_pm.y << "   "
                    << psys[i].acc_pm.z << "   "
                    << psys[i].pot_pm << std::endl;
    }
    output_file.close();

    // Calculate the sum of acceleration
    calcForceSum(psys);


    // Main loop for time integration 


    // Finalize Phantom-GRAPE library
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_close();
#endif
    // Finalize FDPS
    PS::Finalize();
    return 0;
}
