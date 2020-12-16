// Include the standard C++ headers
#include <cmath>
#include <math.h>
#include <cfloat>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <time.h>
// Include the header file of FDPS
#include <particle_simulator.hpp>
// Include the header file of Phantom-GRAPE library
#if defined(ENABLE_PHANTOM_GRAPE_X86)
#include <gp5util.h>
#endif
// Include user-defined headers
#include "macro_defs.hpp"
#include "mathematical_constants.h"
#include "physical_constants.h"
#include "user_defined.hpp"
#include "ic.hpp"
#include "leapfrog.hpp"
#include "prediction.hpp"

// DO NOT CHANGE THE VALUES OF THE MACROS BELOW:
#define SHARED_VARIABLE_TIMESTEP       (1)
#define INDIVIDUAL_TWO_STAGED_TIMESTEP (2)

#define USE_PREDICTED_DENSITY_ETC

#define CODE_MODE (2)


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

void calcDensity(PS::ParticleSystem<FP_sph> & psys,
                 PS::DomainInfo & dinfo, 
                 PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather & tree) {
#if defined(ENABLE_VARIABLE_SMOOTHING_LENGTH)
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    const PS::S64 n_glb = psys.getNumberOfParticleGlobal();
    // Determine the density and the smoothing length so that Eq.(6) in Springel (2005) 
    // holds within a specified accuracy.
    SCF_smth = 1.25;
    PS::S32 iter = 0;
    for (;;) {
        iter++;
        // Compute density, etc.
        tree.calcForceAllAndWriteBack(CalcDensity(), psys, dinfo);
        // Check convergence
        PS::S32 n_compl_loc = 0;
        for (PS::S32 i = 0; i < n_loc; i++) {
            if (psys[i].flag == 1) n_compl_loc++;
        }
        const PS::S64 n_compl = PS::Comm::getSum(n_compl_loc);
        if (n_compl == n_glb) break;
    }
    // Reset SCF_smth
    SCF_smth = 1.0;
#else
    SCF_smth = 1.0;
    tree.calcForceAllAndWriteBack(CalcDensity(), psys, dinfo);
#endif
}

void setEntropy(PS::ParticleSystem<FP_sph>& psys){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        psys[i].setEntropy();
    }
}

void setPressure(PS::ParticleSystem<FP_sph>& psys){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        psys[i].setPressure();
    }
}

void setupIC(PS::ParticleSystem<FP_sph>& psys,
             PS::DomainInfo& dinfo,
             PS::F64 & time_dump,
             PS::F64 & dt_dump,
             PS::F64 & time_end) {
    // Make an initial condition at MPI rank 0.
    PS::BOUNDARY_CONDITION bc;
    PS::F64ort pos_root_domain;
    if (PS::Comm::getRank() == 0) {
#if (INITIAL_CONDITION == 0)
        MakeGlassIC(psys, bc, pos_root_domain,
                    time_dump, dt_dump, time_end);
#elif (INITIAL_CONDITION == 1)
        EvrardTestIC(psys, bc, pos_root_domain,
                     time_dump, dt_dump, time_end, 1);
#else
#error Invalid IC number is specified.
#endif

        // Check the distribution of particles
        const std::string filename = "IC.txt";
        std::ofstream output_file;
        output_file.open(filename.c_str(),std::ios::trunc);
        output_file.setf(std::ios_base::scientific,
                         std::ios_base::floatfield);
        output_file << std::setprecision(15) << std::showpos;
        for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
           output_file << psys[i].pos.x << " "
                       << psys[i].pos.y << " "
                       << psys[i].pos.z << " "
                       << std::endl;
        }
        output_file.close();
    }
    else {
        psys.setNumberOfParticleLocal(0);
    }

    // Broadcast 
    PS::Comm::broadcast(&bc, 1, 0);
    PS::Comm::broadcast(&pos_root_domain, 1, 0);
    PS::Comm::broadcast(&eps_grav, 1, 0);
    PS::Comm::broadcast(&dt_dump, 1, 0);
    PS::Comm::broadcast(&time_dump, 1, 0);
    PS::Comm::broadcast(&time_end, 1, 0);
    PS::Comm::broadcast(&dt_max, 1, 0);

    // Set the boundary condition and the size of the computational domain if needed.
    dinfo.setBoundaryCondition(bc);
    if (bc != PS::BOUNDARY_CONDITION_OPEN) {
        dinfo.setPosRootDomain(pos_root_domain.low_,
                               pos_root_domain.high_);
    }

    // Compute the average mass of SPH particles
    PS::F64 m_sum_loc = 0.0; 
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++)
        m_sum_loc += psys[i].getCharge();
    mass_avg = PS::Comm::getSum(m_sum_loc) / psys.getNumberOfParticleGlobal();

    // Output information to stdout
    if (PS::Comm::getRank() == 0)
        std::cout << "setupIC() is completed." << std::endl;
    //PS::Finalize();
    //std::exit(0);
}


PS::F64 getTimeStep(const PS::ParticleSystem<FP_sph>& psys) {
   PS::F64 dt = DBL_MAX; 
   if (dt_max > 0.0) dt = dt_max;

   // Timescale for SPH system
   for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
#if defined(ENABLE_GRAVITY_INTERACT)
      const PS::F64 acc = std::sqrt((psys[i].acc_grav + psys[i].acc_hydro)
                                   *(psys[i].acc_grav + psys[i].acc_hydro));
      if (acc > 0.0)
          dt = std::min(dt, CFL_dyn * std::sqrt(eps_grav / acc));
#endif
#if defined(ENABLE_HYDRO_INTERACT)
      dt = std::min(dt, psys[i].dt);
#endif
   }
   return PS::Comm::getMinValue(dt);
}

void checkConservativeVariables(const PS::ParticleSystem<FP_sph>& psys) {
    PS::F64    ekin_loc = 0.0;
    PS::F64    epot_loc = 0.0;
    PS::F64    eth_loc  = 0.0; 
    PS::F64vec mom_loc  = 0.0; 
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        ekin_loc += 0.5 * psys[i].mass * psys[i].vel * psys[i].vel;
        epot_loc += 0.5 * psys[i].mass * (psys[i].pot_grav + psys[i].mass / eps_grav);
        eth_loc  += psys[i].mass * psys[i].eng;
        mom_loc  += psys[i].mass * psys[i].vel;
    }
    PS::F64 ekin    = PS::Comm::getSum(ekin_loc);
    PS::F64 epot    = PS::Comm::getSum(epot_loc);
    PS::F64 eth     = PS::Comm::getSum(eth_loc);
    PS::F64vec mom  = PS::Comm::getSum(mom_loc);

    static bool is_initialized = false;
    static PS::F64 emech_ini, etot_ini;
    if (is_initialized == false) {
        emech_ini = ekin + epot;
        etot_ini  = ekin + epot + eth;
        is_initialized = true;
    }

    if (PS::Comm::getRank() == 0){
        const PS::F64 emech = ekin + epot;
        const PS::F64 etot  = ekin + epot + eth;
        const PS::F64 relerr_mech = std::fabs((emech - emech_ini)/emech_ini);
        const PS::F64 relerr_tot  = std::fabs((etot  - etot_ini)/etot_ini);
        std::cout << "-------------------------" << std::endl;
        std::cout << "E_kin  = " << ekin  << std::endl;
        std::cout << "E_pot  = " << epot  << std::endl;
        std::cout << "E_th   = " << eth   << std::endl;
        std::cout << "E_mech = " << emech << " (" << relerr_mech << ")" << std::endl;
        std::cout << "E_tot  = " << etot  << " (" << relerr_tot  << ")" << std::endl;
        std::cout << "Mom (x) = " << mom.x << std::endl;
        std::cout << "Mom (y) = " << mom.y << std::endl;
        std::cout << "Mom (z) = " << mom.z << std::endl;
        std::cout << "-------------------------" << std::endl;
    }
}

void setLevel(PS::ParticleSystem<FP_sph>& psys) {
    const PS::F64 r_crit = 0.4;
    const PS::F64 r2_crit = r_crit * r_crit;
    PS::S32 n_ptcl_at_Lv0 = 0;
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        const PS::F64 r2 = psys[i].pos * psys[i].pos;
        if (r2 > r2_crit) {
            psys[i].level = 0;
            n_ptcl_at_Lv0++;
        }
        else psys[i].level = 1;
    }
    n_ptcl_at_Lv0 = PS::Comm::getSum(n_ptcl_at_Lv0);
    if (PS::Comm::getRank() == 0) {
        std::cout << "n_ptcl_at_Lv0 = " << n_ptcl_at_Lv0 << std::endl;
    }
}

void checkDensityFluctuation(PS::ParticleSystem<FP_sph>& psys) {
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    const PS::S64 n_glb = psys.getNumberOfParticleGlobal();
    // Compute the average of density
    PS::F64 tmp = 0.0;
    for (PS::S32 i = 0; i < n_loc; i++) tmp += psys[i].dens;
    const PS::F64 dens_avg = PS::Comm::getSum(tmp) / n_glb;
    // Compute the dispersion of density
    tmp = 0.0;
    for (PS::S32 i = 0; i < n_loc; i++) 
        tmp += std::pow(psys[i].dens - dens_avg, 2.0);
    const PS::F64 dens_disp = std::sqrt(PS::Comm::getSum(tmp)/n_glb);
    // Output the status 
    const PS::F64 fluc_str = dens_disp / dens_avg;
    if (PS::Comm::getRank() == 0) {
        std::cout << "---------------------------" << std::endl;
        std::cout << "avg.       = " << dens_avg << std::endl;
        std::cout << "disp.      = " << dens_disp << std::endl;
        std::cout << "disp./avg. = " << fluc_str << std::endl;
    }
    // Output data and end the simulation if the fluctuation is small
    const PS::F64 eps = 3.05e-3;
    if (fluc_str < eps) {
        char filename[256];
        sprintf(filename, "result/glass_data.txt");
        psys.writeParticleBinary(filename, &FP_sph::writeBinaryPos);
        if (PS::Comm::getRank() == 0) {
            std::cout << "A glass-like distribution is obtained." << std::endl;
            std::cout << "The particle position data is output as file " 
                      << filename << std::endl;
        }
        PS::Finalize();
        std::exit(0);
    }

}


int main(int argc, char* argv[]){
    // Configure stdout & stderr
    std::cout << std::setprecision(15);
    std::cerr << std::setprecision(15);

    // Initialize FDPS
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
                  << " This is a sample program of "               << std::endl
                  << " Nbody + SPH particle system on FDPS!"       << std::endl
                  << " # of processes is " << nprocs               << std::endl
                  << " # of thread is    " << nthrds               << std::endl
                  << "===========================================" << std::endl;
    }

    // Make an instance of ParticleSystem and initialize it
    PS::ParticleSystem<FP_sph> psys;
    psys.initialize();

    // Make an instance of DomainInfo and initialize it
    PS::DomainInfo dinfo;
    dinfo.initialize();

    // Define local variables
    PS::F64 dt, time_dump, dt_dump, time_end;

    // Make an initial condition and initialize the particle systems
    setupIC(psys, dinfo, time_dump, dt_dump, time_end);

    // Perform domain decomposition 
    dinfo.decomposeDomainAll(psys);
    
    // Perform particle exchange
    psys.exchangeParticle(dinfo);

    // Make tree structures
    const PS::S64 numPtcl = psys.getNumberOfParticleLocal();

    const PS::F32 theta_grav = 0.5;
    PS::TreeForForceLong<Force_grav, EP_grav, EP_grav>::Monopole tree_grav;
    tree_grav.initialize(3 * numPtcl, theta_grav);

    PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather tree_dens;
    tree_dens.initialize(3 * numPtcl);

    PS::TreeForForceShort<Force_hydro, EP_hydro, EP_hydro>::Symmetry tree_hydro;
    tree_hydro.initialize(3 * numPtcl);

    PS::TreeForForceShort<Force_rate, EP_hydro, EP_hydro>::Gather tree_rate;
    tree_rate.initialize(3 * numPtcl);

#if defined(ENABLE_PHANTOM_GRAPE_X86)
    g5_open();
    g5_set_eps_to_all(eps_grav);
#endif

    // Peform force calculations 
    //- Gravity calculations
#if defined(ENABLE_GRAVITY_INTERACT)
    tree_grav.calcForceAllAndWriteBack(CalcGravity<EP_grav>,
                                       CalcGravity<PS::SPJMonopole>,
                                       psys, dinfo);
#endif

    //- SPH calculations
#if defined(ENABLE_HYDRO_INTERACT)
    calcDensity(psys, dinfo, tree_dens);
    setEntropy(psys);
    setPressure(psys);
    tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), psys, dinfo);
#endif

    // Get timestep
    dt = getTimeStep(psys);

    // Calculate energies 
    checkConservativeVariables(psys);

    // Main loop for time integration
#if CODE_MODE == SHARED_VARIABLE_TIMESTEP
    //#########################
    //    Normal mode 
    //#########################
    PS::S32 nstep = 1;
    for (PS::F64 time = 0.0; time < time_end; time += dt, nstep++){
        if (PS::Comm::getRank() == 0) {
            std::cout << "nstep = " << nstep 
                      << " dt = " << dt 
                      << " time = " << time 
                      << " time_end = " << time_end
                      << std::endl;
        }

        // Leap frog: Initial Kick & Full Drift
        InitialKick(psys, dt);
        FullDrift(psys, dt);
        if (dinfo.getBoundaryCondition() != PS::BOUNDARY_CONDITION_OPEN) {
            psys.adjustPositionIntoRootDomain(dinfo);
        }

        // Leap frog: Predict
        Predict(psys, dt);

        // Perform domain decomposition again
        dinfo.decomposeDomainAll(psys);

        // Exchange the particles between the (MPI) processes
        psys.exchangeParticle(dinfo);

        // Peform force calculations
        PS::F64 t_start; 
        //- Gravity calculations
        PS::Comm::barrier(); t_start = PS::GetWtime();
#if defined(ENABLE_GRAVITY_INTERACT)
        tree_grav.calcForceAllAndWriteBack(CalcGravity<EP_grav>,
                                           CalcGravity<PS::SPJMonopole>,
                                           psys, dinfo);
#endif
        PS::Comm::barrier();
        if (PS::Comm::getRank() == 0) std::cout << "t_grav = " << (PS::GetWtime() - t_start) << std::endl;
        //- SPH calculations
        PS::Comm::barrier(); t_start = PS::GetWtime();
#if defined(ENABLE_HYDRO_INTERACT)
        calcDensity(psys, dinfo, tree_dens);
        setPressure(psys);
        tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), psys, dinfo);
#endif
        PS::Comm::barrier();
        if (PS::Comm::getRank() == 0) std::cout << "t_hydro = " << (PS::GetWtime() - t_start) << std::endl;

        // Get a new timestep
        dt = getTimeStep(psys);

        // Leap frog: Final Kick
        FinalKick(psys, dt);

        // Output result files
        if (time > time_dump){
           FileHeader header;
           header.time      = time;
           header.numPtcl   = psys.getNumberOfParticleGlobal();
           char filename[256];
           static int ndump = 1;
           sprintf(filename, "result/sph%05d.txt", ndump);
           psys.writeParticleAscii(filename, header);
           if (PS::Comm::getRank() == 0){
              std::cout << "============================================" << std::endl;
              std::cout << "output " << filename << " at time = " << time << std::endl;
              std::cout << "============================================" << std::endl;
           }
           time_dump += dt_dump;
           ndump++;
        }

        // Calculate energies
        checkConservativeVariables(psys);

        // Check the amplitude of density fluctuation
#if defined(CHECK_DENSITY_FLUCTUATION)
        if (nstep % 100 == 0) checkDensityFluctuation(psys);
#endif

    }
#elif CODE_MODE == INDIVIDUAL_TWO_STAGED_TIMESTEP
    //#########################
    //    Special mode 
    //#########################
    // Set level and save basic quantities
#if defined(USE_PREDICTED_DENSITY_ETC)
    tree_rate.calcForceAllAndWriteBack(CalcRate(), psys, dinfo);
#endif
    setLevel(psys);
    saveDataForPrediction(psys);

    constexpr PS::S32 n_level = 2;
    constexpr PS::S32 ratio = 10;
    PS::F64 level2dt[n_level];
    level2dt[1] = 0.0003;
    level2dt[0] = ratio * level2dt[1];
    dt = level2dt[1]; // reset
    PS::F64 each_dt[n_level] = {0.0, 0.0};
    PS::S32 nstep = 1;
    for (PS::F64 time = 0.0; time < time_end; time += dt, nstep++){
        if (PS::Comm::getRank() == 0) {
            std::cout << "nstep = " << nstep 
                      << " dt = " << dt 
                      << " time = " << time 
                      << " time_end = " << time_end
                      << std::endl;
        }

        // Set dt & sync
        bool sync = false; 
        each_dt[1]  = level2dt[1];
        if (nstep % ratio > 1) {
            each_dt[0] += level2dt[1];
        } else if (nstep % ratio == 1) {
            each_dt[0] = level2dt[1];
        } else {
            each_dt[0] = level2dt[0];
            sync = true;
        }
        
        //if (PS::Comm::getRank() == 0) {
        //    std::cout << "dt[0] = " << each_dt[0] 
        //              << " dt[1] = " << each_dt[1]
        //              << " sync = " << sync << std::endl;
        //}
        //if (nstep == 200) break;

#if 1 // for debug
        // Leap frog: Initial Kick & Full Drift
        InitialKick(psys, each_dt);
        FullDrift(psys, each_dt);
        if (dinfo.getBoundaryCondition() != PS::BOUNDARY_CONDITION_OPEN) {
            psys.adjustPositionIntoRootDomain(dinfo);
        }

        // Leap frog: Predict
        Predict(psys, each_dt);

        // Perform domain decomposition again
        dinfo.decomposeDomainAll(psys);

        // Exchange the particles between the (MPI) processes
        psys.exchangeParticle(dinfo);

        // Peform force calculations
        PS::F64 t_start; 
        //- Gravity calculations
        PS::Comm::barrier(); t_start = PS::GetWtime();
#if defined(ENABLE_GRAVITY_INTERACT)
        tree_grav.calcForceAllAndWriteBack(CalcGravity<EP_grav>,
                                           CalcGravity<PS::SPJMonopole>,
                                           psys, dinfo);
#endif
        PS::Comm::barrier();
        if (PS::Comm::getRank() == 0) std::cout << "t_grav = " << (PS::GetWtime() - t_start) << std::endl;
        //- SPH calculations
        PS::Comm::barrier(); t_start = PS::GetWtime();
#if defined(ENABLE_HYDRO_INTERACT)
        calcDensity(psys, dinfo, tree_dens);
#if defined(USE_PREDICTED_DENSITY_ETC)
        if (sync == false) {
            predictDensityEtc(psys, each_dt);
        }
#endif
        setPressure(psys);
        tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), psys, dinfo);
#endif
        PS::Comm::barrier();
        if (PS::Comm::getRank() == 0) std::cout << "t_hydro = " << (PS::GetWtime() - t_start) << std::endl;

        // Leap frog: Final Kick
        FinalKick(psys, each_dt);

        // if this step is synchronization step, reset level.
        if (sync) {
            // Calculate necessary quantities to predict
            // dens,smth,gradh,divv,rotv of level 0 particles
            // and store them to *_save.
#if defined(USE_PREDICTED_DENSITY_ETC)
            tree_rate.calcForceAllAndWriteBack(CalcRate(), psys, dinfo);
#endif
            setLevel(psys);
            saveDataForPrediction(psys);
        }

        // Output result files
        if (time > time_dump && sync){
           FileHeader header;
           header.time      = time;
           header.numPtcl   = psys.getNumberOfParticleGlobal();
           char filename[256];
           static int ndump = 1;
           sprintf(filename, "result/sph%05d.txt", ndump);
           psys.writeParticleAscii(filename, header);
           if (PS::Comm::getRank() == 0){
              std::cout << "============================================" << std::endl;
              std::cout << "output " << filename << " at time = " << time << std::endl;
              std::cout << "============================================" << std::endl;
           }
           time_dump += dt_dump;
           ndump++;
        }

        // Calculate energies
        checkConservativeVariables(psys);

        // Check the amplitude of density fluctuation
#if defined(CHECK_DENSITY_FLUCTUATION)
        if (nstep % 100 == 0) checkDensityFluctuation(psys);
#endif
#endif // for debug

    }
#else
#error The macro `CODE_MODE` has an invalid value.
#endif

#if defined(ENABLE_PHANTOM_GRAPE_X86)
    g5_close();
#endif
    // Finalize FDPS
    PS::Finalize();
    return 0;
}

