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
#include <direct_sum.hpp>
#include <ewald.hpp>
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

void setupIC(const PS::S32 NP,
             const bool charge_neutral,
             const PS::BOUNDARY_CONDITION bc,
             PS::ParticleSystem<FP_nbody> & psys,
             PS::DomainInfo & dinfo) {
    // Initialize pseudo-random generator
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::S32 my_rank = PS::Comm::getRank();
    const long seed=19810614 + my_rank;
    srand48(seed);
    const PS::S64 NP_tot = NP * n_proc;
    // Set psys except for mass
    psys.setNumberOfParticleLocal(NP);
    for(PS::S32 i = 0; i < NP; i++){
        psys[i].id  = i;
        psys[i].pos = PS::F64vec(drand48(), drand48(), drand48());
        psys[i].vel = 0.0;
    }
    // Set mass 
    PS::F64 msum = 0.0;
    for (PS::S32 i = 0; i < NP; i++){
        msum += (psys[i].mass  = drand48() * (1.0/NP_tot));
    }
    msum = PS::Comm::getSum(msum);
    // Adjust mass so that the system is charge neutral.
    if (charge_neutral) {
        for (PS::S32 i = 0; i < NP; i++){
            psys[i].mass -= msum / NP_tot;
        }
    }
    // Check total charge
    msum = 0.0;
    for (PS::S32 i = 0; i < NP; i++) {
        msum += psys[i].mass;
    }
    msum = PS::Comm::getSum(msum);
    if (my_rank == 0) std::cout << "msum (@IC) = " << msum << std::endl;

    // Set domain information
    dinfo.setBoundaryCondition(bc);
    dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0),
                           PS::F64vec(1.0, 1.0, 1.0));
}


int main(int argc, char *argv[]) {
    // Configure stdout & stderr
    std::cout << std::setprecision(15);
    std::cerr << std::setprecision(15);

    //* Initialize FDPS
    PS::Initialize(argc, argv);

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

    // Make a directory
    char dir_name[1024];
    sprintf(dir_name,"./result");
    makeOutputDirectory(dir_name);

    // Make an instance of ParticleSystem and initialize it
    PS::ParticleSystem<FP_nbody> psys;
    psys.initialize();

    // Make an instance of DomainInfo and initialize it
    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);

    // Make an initial condition 
    // NP := number of particles per process (different from main-dev.cpp!)
    const PS::S32 NP=4096;
    //const PS::S32 NP=1048576;
    const bool charge_neutral=true;
    //const bool charge_neutral=false;
    //const PS::BOUNDARY_CONDITION bc=PS::BOUNDARY_CONDITION_OPEN;
    const PS::BOUNDARY_CONDITION bc=PS::BOUNDARY_CONDITION_PERIODIC_XYZ;

    setupIC(NP, charge_neutral, bc, psys, dinfo);

    // Perform domain decomposition and particle exchange
    dinfo.decomposeDomainAll(psys);
    psys.exchangeParticle(dinfo);
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    const PS::S32 n_glb = psys.getNumberOfParticleGlobal();

    // Initialize Phantom-GRAPE library
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_open();
    g5_set_eps_to_all(FP_nbody::eps);
#endif

    // Make an instance of TreeForForce, initialize it, and perform force calculation of PP part
    const PS::S32 n_tot = 3 * n_loc;
    const PS::S32vec n_cell = PS::S32vec(32, 32, 32);
    const PS::S32 icut = 2;
    const PS::F32 theta = 0.0;
    //const PS::F32 theta = 0.5;
    const PS::U32 n_leaf_limit = 8;
    const PS::U32 n_group_limit = 64;
    PS::TreeForForceLong<FP_nbody, FP_nbody, FP_nbody>::QuadrupoleGeometricCenterParticleMeshMultipole tree;
    tree.initialize(n_tot, n_cell, icut, theta, n_leaf_limit, n_group_limit);
    tree.calcForceAllAndWriteBack(CalcForceEpEp<FP_nbody>,
                                  CalcForceEpSp<PS::SPJQuadrupoleGeometricCenterPMMM>,
                                  psys, dinfo);

    // Make an instance of ParticleMeshMultipole class, initialize it, and force calculation of PM part
    //const PS::S32 p = 4;
    //const PS::S32 p = 5;
    //const PS::S32 p = 6;
    const PS::S32 p = 7;
    //const PS::S32 p = 8;
    //const PS::S32 p = 9;
    //const PS::S32 p = 10;
    //const PS::S32 p = 11;
    //const PS::S32 p = 12;
    const PS::S32 p_spj2mm = 5;
    const bool use_mpifft = true;
    //const bool use_mpifft = false;
    PS::PMM::ParticleMeshMultipole<FP_nbody, FP_nbody> pm;
    pm.initialize(p, p_spj2mm, FFTW_MEASURE, use_mpifft);
    pm.calcForceAllAndWriteBack(tree,dinfo,psys); // to make green function
    if (PS::Comm::getRank() == 0) {
        std::cout << "To check time of the calculation of Green function:" << std::endl;
        PS::PMM::TimeProfilePMM time_prof = pm.getTimeProfile();
        time_prof.dump();
        std::cout << "-----" << std::endl;
    }
    pm.clearTimeProfile();

    // Measure performance
    const PS::S32 n_trials = 1;
    //const PS::S32 n_trials = 32;
    PS::F64 t_start, t_end;
    PS::F64 etime_pp = 0.0, etime_pm = 0.0;
    for (PS::S32 n = 0; n < n_trials; n++) {
        // Measure PP part
        PS::Comm::barrier();
        t_start = PS::GetWtime();
        tree.calcForceAllAndWriteBack(CalcForceEpEp<FP_nbody>,
                                      CalcForceEpSp<PS::SPJQuadrupoleGeometricCenterPMMM>,
                                      psys, dinfo);
        PS::Comm::barrier();
        t_end = PS::GetWtime();
        etime_pp += t_end - t_start;
        // Measure PM part
        PS::Comm::barrier();
        t_start = PS::GetWtime();
        pm.calcForceAllAndWriteBack(tree,dinfo,psys);
        PS::Comm::barrier();
        t_end = PS::GetWtime();
        etime_pm += t_end - t_start;
    }
    etime_pp /= n_trials;
    etime_pm /= n_trials;
    PS::PMM::TimeProfilePMM time_prof = pm.getTimeProfile();
    time_prof /= n_trials;

    // Output the measured performce 
    if (PS::Comm::getRank() == 0) {
        std::cout << "etime (PP) = " << etime_pp << std::endl;
        std::cout << "etime (PM) = " << etime_pm << std::endl;
        time_prof.dump();
    }


    // Finalize Phantom-GRAPE library
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_close();
#endif
    // Finalize FDPS
    PS::Finalize();
    return 0;
}
