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
    const int NP=128;
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
    } else {
        psys.setNumberOfParticleLocal(0);
    }

    // Set domain information
    //dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0),
                           PS::F64vec(1.0, 1.0, 1.0));
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

    // Make an instance of PMMM object
    //const PS::S32 PFMM=7;
    const PS::S32 PFMM=5;
    const PS::S32 ICUT=2;
    PS::ParticleMeshMultipole<PFMM,ICUT> pmmm;
    pmmm.initialize();
    pmmm.calcForceAllAndWriteBack(psys,dinfo);

    // Initialize Phantom-GRAPE library
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_open();
    g5_set_eps_to_all(FP_nbody::eps);
#endif

    // Perform force calculation


    // Main loop for time integration 


    // Finalize Phantom-GRAPE library
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_close();
#endif
    // Finalize FDPS
    PS::Finalize();
    return 0;
}
