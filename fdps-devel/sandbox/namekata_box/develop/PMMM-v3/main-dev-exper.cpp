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

    if (PS::Comm::getRank() == 0) {
        psys.setNumberOfParticleLocal(NP);
        // Initialize pseudo-random generator
        const long seed=19810614;
        srand48(seed);
        // Set psys except for mass
        for(PS::S32 i = 0; i < NP; i++){
            psys[i].id  = i;
            psys[i].pos = PS::F64vec(drand48(), drand48(), drand48());
            psys[i].vel = 0.0;
        }
        // Set mass 
        PS::F64 msum = 0.0;
        for (PS::S32 i = 0; i < NP; i++){
            msum += (psys[i].mass  = drand48() * (1.0/NP));
        }
        // Adjust mass so that the system is charge neutral.
        if (charge_neutral) {
            for (PS::S32 i = 0; i < NP; i++){
                psys[i].mass -= msum / NP;
            }
        }
        // Output 
        std::string filename = "IC.txt";
        std::ofstream output_file;
        output_file.open(filename.c_str(), std::ios::trunc);
        msum = 0.0;
        for (PS::S32 i = 0; i < NP; i++) {
            output_file << psys[i].pos.x << " "
                        << psys[i].pos.y << " "
                        << psys[i].pos.z << " "
                        << psys[i].mass <<  std::endl;
            msum += psys[i].mass;
        }
        output_file.close();
        std::cout << "msum (@IC) = " << msum << std::endl;
    } else {
        psys.setNumberOfParticleLocal(0);
    }

    // Set domain information
    dinfo.setBoundaryCondition(bc);
    dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0),
                           PS::F64vec(1.0, 1.0, 1.0));
}

void calcForceSum(PS::ParticleSystem<FP_nbody> & psys) {
    // Output the result of PP part
    {
        std::ostringstream ss;
        ss << "PP" << std::setfill('0') << std::setw(5) << PS::Comm::getRank() << ".txt";
        std::string filename = ss.str();
        std::ofstream output_file;
        output_file.open(filename.c_str(), std::ios::trunc);
        for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
#if 0
            output_file << psys[i].pos.x << "   "
                        << psys[i].pos.y << "   "
                        << psys[i].pos.z << "   "
                        << psys[i].mtot << std::endl;
#else
            output_file << psys[i].pos.x << "   "
                        << psys[i].pos.y << "   "
                        << psys[i].pos.z << "   "
                        << psys[i].acc.x << "   "
                        << psys[i].acc.y << "   "
                        << psys[i].acc.z << "   "
                        << psys[i].pot << std::endl;
#endif
        }
        output_file.close();
    }
    // Output the result of PM part
    {
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
    }
    // Compute the sum of PP and PM parts
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    PS::F64vec fpp_loc(0.0);
    PS::F64vec fpm_loc(0.0);
    for (PS::S32 i = 0; i < n_loc; i++) {
        fpp_loc += psys[i].mass * psys[i].acc;
        fpm_loc += psys[i].mass * psys[i].acc_pm;
        psys[i].pot += psys[i].pot_pm;
        psys[i].acc += psys[i].acc_pm;
        psys[i].pot_pm = 0.0;
        psys[i].acc_pm = 0.0;
    }
    PS::F64vec fpp_tot;
    PS::F64vec fpm_tot;
    fpp_tot = PS::Comm::getSum(fpp_loc);
    fpm_tot = PS::Comm::getSum(fpm_loc);
    if (PS::Comm::getRank() == 0) {
        std::cout << "PP ftot : " << fpp_tot << std::endl;
        std::cout << "PM ftot : " << fpm_tot << std::endl;
    }
    // Output the sum
    {
        std::ostringstream ss;
        ss << "total" << std::setfill('0') << std::setw(5) << PS::Comm::getRank() << ".txt";
        std::string filename = ss.str();
        std::ofstream output_file;
        output_file.open(filename.c_str(), std::ios::trunc);
        for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
            output_file << psys[i].pos.x << "   "
                        << psys[i].pos.y << "   "
                        << psys[i].pos.z << "   "
                        << psys[i].acc.x << "   "
                        << psys[i].acc.y << "   "
                        << psys[i].acc.z << "   "
                        << psys[i].pot << std::endl;
        }
        output_file.close();
    }
}

void calcForceDirect(PS::ParticleSystem<FP_nbody> & psys) { 
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cout << "the direct sum method started." << std::endl;
    PS::Comm::barrier();

    PS::DirectSum::DirectSum ds;
    ds.initialize();
    ds.calcForceAllAndWriteBack(psys);

    // Output 
    std::ostringstream ss;
    ss << "direct" << std::setfill('0') << std::setw(5) << PS::Comm::getRank() << ".txt";
    std::string filename = ss.str();
    std::ofstream output_file;
    output_file.open(filename.c_str(), std::ios::trunc);
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        output_file << psys[i].pos.x << "   "
                    << psys[i].pos.y << "   "
                    << psys[i].pos.z << "   "
                    << psys[i].acc_exact.x << "   "
                    << psys[i].acc_exact.y << "   "
                    << psys[i].acc_exact.z << "   "
                    << psys[i].pot_exact << std::endl;
    }
    output_file.close();

    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cout << "the direct sum method ended." << std::endl;
    PS::Comm::barrier();
}

void calcForceEwald(PS::ParticleSystem<FP_nbody> & psys) {
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cout << "the Ewald method started." << std::endl;
    PS::Comm::barrier();

    const PS::F64ort pos_unit_cell = PS::F64ort(PS::F64vec(0.0), PS::F64vec(1.0));
    PS::Ewald::Ewald ewald;
    ewald.initialize(pos_unit_cell);
    ewald.calcForceAllAndWriteBack(psys);

    // Output 
    std::ostringstream ss;
    ss << "Ewald" << std::setfill('0') << std::setw(5) << PS::Comm::getRank() << ".txt";
    std::string filename = ss.str();
    std::ofstream output_file;
    output_file.open(filename.c_str(), std::ios::trunc);
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        output_file << psys[i].pos.x << "   "
                    << psys[i].pos.y << "   "
                    << psys[i].pos.z << "   "
                    << psys[i].acc_exact.x << "   "
                    << psys[i].acc_exact.y << "   "
                    << psys[i].acc_exact.z << "   "
                    << psys[i].pot_exact << std::endl;
    }
    output_file.close();

    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cout << "the Ewald method ended." << std::endl;
    PS::Comm::barrier();
}

void calcExactSolution(PS::ParticleSystem<FP_nbody> & psys,
                       const PS::BOUNDARY_CONDITION bc) {
    if (bc == PS::BOUNDARY_CONDITION_OPEN) {
        calcForceDirect(psys);
    } else if (bc == PS::BOUNDARY_CONDITION_PERIODIC_XYZ) {
        calcForceEwald(psys);
    }
}

void perfErrorAnalysis(PS::ParticleSystem<FP_nbody> & psys) {
    // Local class
    class ErrorInfo {
    public:
        PS::F64 p_l1, g_l1;
        // p := potential, g := gradient of potential,
        // l1 := L1 relative error.
    };
    const PS::S32 n_loc = psys.getNumberOfParticleLocal(); 
    const PS::S32 n_glb = psys.getNumberOfParticleGlobal();
    std::vector<ErrorInfo> err(n_loc);
    PS::F64 err_p_rms_loc {0.0};
    PS::F64 err_g_rms_loc {0.0};
    PS::F64 errmax_p {0.0};
    PS::F64 errmax_g {0.0};
    bool errmax_p_is_defined {false};
    bool errmax_g_is_defined {false};
    PS::S32 i_max_p, i_max_g;
    for (PS::S32 i=0; i<n_loc; i++) {
        // potential
        err[i].p_l1 = std::fabs(psys[i].pot - psys[i].pot_exact)
                    / std::fabs(psys[i].pot_exact);
        if (err[i].p_l1 > errmax_p) {
            errmax_p_is_defined = true;
            errmax_p = err[i].p_l1;
            i_max_p = i;
        }
        err_p_rms_loc += err[i].p_l1 * err[i].p_l1;
        // gradient
        const PS::F64vec afid  = psys[i].acc_exact;
        const PS::F64vec adiff = psys[i].acc - psys[i].acc_exact;
        err[i].g_l1 = std::fabs(std::sqrt(adiff * adiff))
                    / std::fabs(std::sqrt(afid  * afid));
        if (err[i].g_l1 > errmax_g) {
            errmax_g_is_defined = true;
            errmax_g = err[i].g_l1;
            i_max_g = i;
        }
        err_g_rms_loc += err[i].g_l1 * err[i].g_l1;
    }
    const PS::F64 err_p_rms = std::sqrt(PS::Comm::getSum(err_p_rms_loc)/n_glb);
    const PS::F64 err_g_rms = std::sqrt(PS::Comm::getSum(err_g_rms_loc)/n_glb);
    if (PS::Comm::getRank() == 0) {
        std::cout << "RMS relative error (potential) = " << err_p_rms << std::endl;
        std::cout << "RMS relative error (gradient)  = " << err_g_rms << std::endl;
    }
    if (errmax_p) {
        std::cout << "my_rank = " << PS::Comm::getRank()
                  << " i_max_p = " << i_max_p 
                  << " pot  = " << psys[i_max_p].pot
                  << " pot_exact = " << psys[i_max_p].pot_exact
                  << std::endl;
    }
    if (errmax_g) {
        std::cout << "my_rank = " << PS::Comm::getRank()
                  << " i_max_g = " << i_max_g
                  << " acc  = " << psys[i_max_g].acc
                  << " acc_exact = " << psys[i_max_g].acc_exact
                  << std::endl;
    }

    // Output 
    {
        std::ostringstream ss;
        ss << "err_p_" << std::setfill('0') << std::setw(5) << PS::Comm::getRank() << ".txt";
        std::string filename = ss.str();
        std::ofstream output_file;
        output_file.open(filename.c_str(), std::ios::trunc);
        for (PS::S32 i = 0; i < n_loc; i++) {
            output_file << err[i].p_l1 << std::endl;
        }
        output_file.close();
    }
    {
        std::ostringstream ss;
        ss << "err_g_" << std::setfill('0') << std::setw(5) << PS::Comm::getRank() << ".txt";
        std::string filename = ss.str();
        std::ofstream output_file;
        output_file.open(filename.c_str(), std::ios::trunc);
        for (PS::S32 i = 0; i < n_loc; i++) {
            output_file << err[i].g_l1 << std::endl;
        }
        output_file.close();
    }
}

int main(int argc, char *argv[]) {
    // Configure stdout & stderr
    std::cout << std::setprecision(15);
    std::cerr << std::setprecision(15);

    // Set parameters
    //const PS::S32 NP=128;
    //const PS::S32 NP=256;
    //const PS::S32 NP=512;
    //const PS::S32 NP=1024;
    const PS::S32 NP=4096;
    //const PS::S32 NP=8192;
    //const PS::S32 NP=16384;
    //const PS::S32 NP=32768;
    //const PS::S32 NP=65536;
    const bool charge_neutral=true;
    //const bool charge_neutral=false;
    //const PS::BOUNDARY_CONDITION bc=PS::BOUNDARY_CONDITION_OPEN;
    const PS::BOUNDARY_CONDITION bc=PS::BOUNDARY_CONDITION_PERIODIC_XYZ;

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

    // Make instances of ParticleSystem and initialize them
    PS::ParticleSystem<FP_nbody> psys0, psys;
    psys0.initialize();
    psys.initialize();

    // Make an instance of DomainInfo and initialize it
    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);

    // Make an initial condition and initialize the particle system
    setupIC(NP, charge_neutral, bc, psys, dinfo);

    // Perform domain decomposition and particle exchange
    dinfo.decomposeDomainAll(psys);
    psys.exchangeParticle(dinfo);

    // Copy psys to psys0 (psys0 is used to store the result of the case of \theta = 0)
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    const PS::S32 n_glb = psys.getNumberOfParticleGlobal();
    psys0.setNumberOfParticleLocal(n_loc);
    for (PS::S32 i = 0; i<n_loc; i++) {
        psys0[i].id  = psys[i].id;
        psys0[i].mass = psys[i].mass;
        psys0[i].pos = psys[i].pos;
        psys0[i].vel = psys[i].vel;
    }

    // Make instances of TreeForForce and initialize it
    const PS::S32 n_tot = 3 * psys.getNumberOfParticleLocal();
    const PS::S32vec n_cell = PS::S32vec(8, 8, 8);
    //const PS::S32vec n_cell = PS::S32vec(16, 16, 16);
    //const PS::S32vec n_cell = PS::S32vec(32, 32, 32);
    const PS::S32 icut = 2;
    const PS::F32 theta = 0.0;
    //const PS::F32 theta = 0.5;
    const PS::U32 n_leaf_limit = 8;
    const PS::U32 n_group_limit = 64;
    constexpr PS::S32 p_spj = 7;
    PS::TreeForForceLong<FP_nbody, FP_nbody, FP_nbody>::MultipoleGeometricCenterParticleMeshMultipole<p_spj> tree0, tree;
    tree0.initialize(n_tot, n_cell, icut, 0.0, n_leaf_limit, n_group_limit);
    tree.initialize(n_tot, n_cell, icut, theta, n_leaf_limit, n_group_limit);

    // Initialize Phantom-GRAPE library
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_open();
    g5_set_eps_to_all(FP_nbody::eps);
#endif

    // Perform force calculation (PP part, \theta = 0)
    PS::F64 t_start, t_end;
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cout << "PP part (\\theta = 0) started."  << std::endl;
    PS::Comm::barrier();
    t_start = PS::GetWtime();
    tree0.calcForceAllAndWriteBack(CalcForceEpEp<FP_nbody>,
                                   CalcForceEpSp<p_spj>,
                                   psys0, dinfo);
    PS::Comm::barrier();
    t_end = PS::GetWtime();
    if (PS::Comm::getRank() == 0) {
        std::cout << "PP part (\\theta = 0) ended."  << std::endl;
        std::cout << "etime (PP, \\theta = 0) = " << t_end - t_start << std::endl;
    }
    PS::Comm::barrier();

    // Perform force calculation (PP part, \theta > 0)
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cout << "PP part (\\theta > 0) started." << std::endl;
    PS::Comm::barrier();
    t_start = PS::GetWtime();
    tree.calcForceAllAndWriteBack(CalcForceEpEp<FP_nbody>,
                                  CalcForceEpSp<p_spj>,
                                  psys, dinfo);
    PS::Comm::barrier();
    t_end = PS::GetWtime();
    if (PS::Comm::getRank() == 0) {
        std::cout << "PP part (\\theta > 0) ended." << std::endl;
        std::cout << "etime (PP, \\theta > 0) = " << t_end - t_start << std::endl;
    }
    PS::Comm::barrier();

    //PS::Finalize();
    //std::exit(0);

    //====================================
    // Check the correctness of PP part
    //====================================
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cout << "Checking the PP part...." << std::endl;
    PS::Comm::barrier();
    // (1) mtot
    const PS::F64 eps4mtot = 1.0e-10;
    for (PS::S32 i = 0; i < n_loc; i++) {
        const PS::F64 diff = psys[i].mtot - psys0[i].mtot;
        const PS::F64 avg  = 0.5 * (psys[i].mtot + psys0[i].mtot);
        if (avg != 0.0) {
            const PS::F64 reldiff = std::abs(diff/avg);
            if (reldiff > eps4mtot) {
                std::cout << "rank = " << PS::Comm::getRank()
                          << " i = " << i 
                          << " id = " << psys[i].id
                          << " reldiff = " << reldiff
                          << " mtot0 = " << psys0[i].mtot
                          << " mtot = " << psys[i].mtot
                          << std::endl;
            }
        }
    }
    // (2) pot & acc
    PS::F64 err_pp_pot_loc {0.0};
    PS::F64 err_pp_grad_loc {0.0};
    bool f_pp_pot {false}, f_pp_grad {false};
    PS::S32 i_pp_pot, i_pp_grad;
    PS::S64 id_pp_pot, id_pp_grad;
    PS::F64 errmax_pp_pot {0.0};
    PS::F64 errmax_pp_grad {0.0};
    for (PS::S32 i = 0; i < n_loc; i++) {
        PS::F64 d2 = (psys[i].pot - psys0[i].pot)
                    *(psys[i].pot - psys0[i].pot);
        PS::F64 n2 = psys0[i].pot * psys0[i].pot;
        if (n2 > 0.0) {
            const PS::F64 err = d2/n2;
            err_pp_pot_loc += err;
            if (err > errmax_pp_pot) {
                f_pp_pot = true;
                errmax_pp_pot = err;
                i_pp_pot = i;
                id_pp_pot = psys[i].id;
            }
        }
        d2 = (psys[i].acc - psys0[i].acc)
            *(psys[i].acc - psys0[i].acc);
        n2 = psys0[i].acc * psys0[i].acc;
        if (n2 > 0.0) {
            const PS::F64 err = d2/n2;
            err_pp_grad_loc += err;
            if (err > errmax_pp_grad) {
                f_pp_grad = true;
                errmax_pp_grad = err;
                i_pp_grad = i;
                id_pp_grad = psys[i].id;
            }
        }
    }
    const PS::F64 err_pp_pot  = std::sqrt(PS::Comm::getSum(err_pp_pot_loc)/n_glb);
    const PS::F64 err_pp_grad = std::sqrt(PS::Comm::getSum(err_pp_grad_loc)/n_glb);
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) {
        std::cout << "RMS relative error (PP, pot., \\theta=0 & >0)  = " << err_pp_pot << std::endl;
        std::cout << "RMS relative error (PP, grad., \\theta=0 & >0) = " << err_pp_grad << std::endl;
    }
    if (f_pp_pot) {
        std::cout << "my_rank = " << PS::Comm::getRank()
                  << " i_pp_pot = " << i_pp_pot 
                  << " id_pp_pot = " << id_pp_pot 
                  << " pot = " << psys[i_pp_pot].pot
                  << " pot0 = " << psys0[i_pp_pot].pot 
                  << std::endl;
    }
    if (f_pp_grad) {
        std::cout << "my_rank = " << PS::Comm::getRank()
                  << " i_pp_grad = " << i_pp_grad 
                  << " id_pp_grad = " << id_pp_grad 
                  << " acc = " << psys[i_pp_grad].acc
                  << " acc0 = " << psys0[i_pp_grad].acc
                  << std::endl;
    }
    PS::Comm::barrier();

    // Make instances of ParticleMeshMultipole class and initialize them
    const PS::S32 p = 2;
    //const PS::S32 p = 3;
    //const PS::S32 p = 4;
    //const PS::S32 p = 5;
    //const PS::S32 p = 6;
    //const PS::S32 p = 7;
    //const PS::S32 p = 8;
    //const PS::S32 p = 9;
    //const PS::S32 p = 10;
    //const PS::S32 p = 11;
    //const PS::S32 p = 12;
    const PS::S32 p_spj2mm = 5;
    PS::PMM::ParticleMeshMultipole<FP_nbody, FP_nbody> pm0, pm;
    pm0.initialize(p, p_spj2mm);
    pm.initialize(p, p_spj2mm);

    // Perform force calculation (PM part, \theta = 0)
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cout << "PM part (\\theta = 0) started." << std::endl;
    PS::Comm::barrier();
    t_start = PS::GetWtime();

    pm0.calcForceAllAndWriteBack(tree0,dinfo,psys0);

    PS::Comm::barrier();
    t_end = PS::GetWtime();
    if (PS::Comm::getRank() == 0) {
        std::cout << "PM part (\\theta = 0) ended." << std::endl;
        std::cout << "etime (PM, \\theta = 0) = " << t_end - t_start << std::endl;
    }
    PS::Comm::barrier();

    // Perform force calculation (PM part, \theta > 0)
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cout << "PM part (\\theta > 0) started." << std::endl;
    PS::Comm::barrier();
    t_start = PS::GetWtime();

    pm.calcForceAllAndWriteBack(tree,dinfo,psys);

    PS::Comm::barrier();
    t_end = PS::GetWtime();
    if (PS::Comm::getRank() == 0) {
        std::cout << "PM part (\\theta > 0) ended." << std::endl;
        std::cout << "etime (PM, \\theta > 0) = " << t_end - t_start << std::endl;
    }
    PS::Comm::barrier();

    //====================================
    // Check the correctness of PM part
    //====================================
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) std::cout << "Checking the PM part...." << std::endl;
    PS::Comm::barrier();
    // (2) pot & acc
    PS::F64 err_pm_pot_loc {0.0};
    PS::F64 err_pm_grad_loc {0.0};
    bool f_pm_pot {false}, f_pm_grad {false};
    PS::S32 i_pm_pot, i_pm_grad;
    PS::S64 id_pm_pot, id_pm_grad;
    PS::F64 errmax_pm_pot {0.0};
    PS::F64 errmax_pm_grad {0.0};
    for (PS::S32 i = 0; i < n_loc; i++) {
        PS::F64 d2 = (psys[i].pot_pm - psys0[i].pot_pm)
                    *(psys[i].pot_pm - psys0[i].pot_pm);
        PS::F64 n2 = psys0[i].pot_pm * psys0[i].pot_pm;
        if (n2 > 0.0) {
            const PS::F64 err = d2/n2;
            err_pm_pot_loc += err;
            if (err > errmax_pm_pot) {
                f_pm_pot = true;
                errmax_pm_pot = err;
                i_pm_pot = i;
                id_pm_pot = psys[i].id;
            }
        }
        d2 = (psys[i].acc_pm - psys0[i].acc_pm)
            *(psys[i].acc_pm - psys0[i].acc_pm);
        n2 = psys0[i].acc_pm * psys0[i].acc_pm;
        if (n2 > 0.0) {
            const PS::F64 err = d2/n2;
            err_pm_grad_loc += err;
            if (err > errmax_pm_grad) {
                f_pm_grad = true;
                errmax_pm_grad = err;
                i_pm_grad = i;
                id_pm_grad = psys[i].id;
            }
        }
    }
    const PS::F64 err_pm_pot  = std::sqrt(PS::Comm::getSum(err_pm_pot_loc)/n_glb);
    const PS::F64 err_pm_grad = std::sqrt(PS::Comm::getSum(err_pm_grad_loc)/n_glb);
    PS::Comm::barrier();
    if (PS::Comm::getRank() == 0) {
        std::cout << "RMS relative error (PM, pot., \\theta=0 & >0)  = " << err_pm_pot << std::endl;
        std::cout << "RMS relative error (PM, grad., \\theta=0 & >0) = " << err_pm_grad << std::endl;
    }
    if (f_pm_pot) {
        std::cout << "my_rank = " << PS::Comm::getRank() 
                  << " i_pm_pot = " << i_pm_pot 
                  << " id_pm_pot = " << id_pm_pot
                  << " pot = " << psys[i_pm_pot].pot_pm
                  << " pot0 = " << psys0[i_pm_pot].pot_pm 
                  << std::endl;
    }
    if (f_pm_grad) {
        std::cout << "my_rank = " << PS::Comm::getRank()
                  << " i_pm_grad = " << i_pm_grad
                  << " id_pm_grad = " << id_pm_grad
                  << " acc = " << psys[i_pm_grad].acc_pm
                  << " acc0 = " << psys0[i_pm_grad].acc_pm
                  << std::endl;
    }
    PS::Comm::barrier();

    // Calculate the sum of acceleration
    calcForceSum(psys);

    // Calculate exact solution
    calcExactSolution(psys, bc);

    // Compare the calculated force and the exact solution
    perfErrorAnalysis(psys);

    // Main loop for time integration 


    // Finalize Phantom-GRAPE library
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_close();
#endif
    // Finalize FDPS
    PS::Finalize();
    return 0;
}
