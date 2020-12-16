// Include the C++ headers
#include "common.h"
// Include the header file of FDPS
#include <particle_simulator.hpp>
// Include the header file of Phantom-GRAPE library
#if defined(ENABLE_PHANTOM_GRAPE_X86)
#include <gp5util.h>
#endif
// Include user-defined headers
#include "macro_defs.h"
#include "debug_utilities.h"
#include "mathematical_constants.h"
#include "physical_constants.h"
#include "run_parameters.h"
#include "user_defined.h"
#include "ic.h"
#include "leapfrog.h"
#include "cooling_heating.h"
#include "star_formation.h"
#include "stellar_feedback.h"
#include "io.h"
#include "hydrodynamics.h"
#include "timestep.h"

void setCircularVelocity(PS::ParticleSystem<FP_gas>& psys){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        const PS::F64vec acc = psys[i].acc_grav + psys[i].acc_hydro;
        const PS::F64 r = std::sqrt(psys[i].pos * psys[i].pos);
        const PS::F64 R = std::sqrt(psys[i].pos.x * psys[i].pos.x 
                                   +psys[i].pos.y * psys[i].pos.y);
        const PS::F64 phi = std::atan2(psys[i].pos.y, psys[i].pos.x);
        const PS::F64 theta = std::atan2(R, psys[i].pos.z);
        PS::F64vec base_vect_r;
        base_vect_r.x = std::sin(theta) * std::cos(phi);
        base_vect_r.y = std::sin(theta) * std::sin(phi);
        base_vect_r.z = std::cos(theta);
        const PS::F64 vel_circ = std::sqrt(r * std::abs(acc * base_vect_r));
        psys[i].vel.x = - vel_circ * std::sin(phi);
        psys[i].vel.y =   vel_circ * std::cos(phi);
        psys[i].vel.z = 0.0;
    }
#if 0
   // [for debug] 
   std::string filename;
   std::ostringstream ss;
   std::ofstream output_file;
   // Output the velocity field
   ss << "velc_fld" << std::setfill('0') << std::setw(5) << PS::Comm::getRank() << ".txt";
   filename = ss.str();
   output_file.open(filename.c_str(),std::ios::trunc);
   output_file.setf(std::ios_base::scientific,
                    std::ios_base::floatfield);
   output_file << std::setprecision(15) << std::showpos;
   for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
       output_file << psys[i].pos.x << " " 
                   << psys[i].pos.y << " "
                   << psys[i].vel.x << " " 
                   << psys[i].vel.y << " "
                   << std::endl;
   }
   output_file.close();
   // Output the rotation curve
   ss.str(""); // clear
   ss << "rot_curve" << std::setfill('0') << std::setw(5) << PS::Comm::getRank() << ".txt";
   filename = ss.str();
   output_file.open(filename.c_str(),std::ios::trunc);
   output_file.setf(std::ios_base::scientific,
                    std::ios_base::floatfield);
   output_file << std::setprecision(15) << std::showpos;
   for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
       const PS::F64 r = std::sqrt(psys[i].pos * psys[i].pos);
       const PS::F64 R = std::sqrt(psys[i].pos.x * psys[i].pos.x
                                  +psys[i].pos.y * psys[i].pos.y);
       const PS::F64 v = std::sqrt(psys[i].vel * psys[i].vel);
       const PS::F64 T = 2.0 * math_const::pi * r / v;
       output_file << R << " " << v << " " << T << std::endl;
   }
   output_file.close();
   //PS::Finalize();
   //std::exit(0);
#endif
}

void printHelp() {
    std::cout << "-r NUMBER" << std::endl;
    std::cout << "       Execute in the restart run mode. NUMBER is file number used for restart." << std::endl;
    std::cout << "       If NUMBER <= 0, execute in the initial run mode." << std::endl;
    std::cout << "-h     Show this help" << std::endl;
}

int main(int argc, char* argv[]){
    // Configure stdout & stderr
    std::cout << std::setprecision(15);
    std::cerr << std::setprecision(15);

    // Initialize FDPS and display parallelization information
    PS::Initialize(argc, argv);
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::S32 my_rank = PS::Comm::getRank();
    const PS::S32 n_thrd = PS::Comm::getNumberOfThread();
    if (PS::Comm::getRank() == 0) {
        std::cout << "===========================================" << std::endl
                  << " This is ASURA-FDPS code."                   << std::endl
                  << " # of processes is " << n_proc               << std::endl
                  << " # of thread is    " << n_thrd               << std::endl
                  << "===========================================" << std::endl;
    }

    // Initialize CELib
    CELibInit();

    // Process command-line arguments
    PS::S32 rst_file_num = -1;
    PS::S32 c;
    opterr = 0; // A global variable defined in unistd.h
    while((c = getopt(argc,argv,"r:h")) != -1){
        switch(c){
        case 'r':
            rst_file_num = atof(optarg);
            if (PS::Comm::getRank() == 0) {
                std::cout << "Restart file number = " << rst_file_num << std::endl;
            }
            break;
        case 'h':
            if (PS::Comm::getRank() == 0) {
                printHelp();
            }
            PS::Finalize();
            return 0;
        default:
            if (PS::Comm::getRank() == 0) {
                std::cout << "No such option! Available options are here." << std::endl;
                printHelp();
            }
            PS::Abort();
        }
    }
    
    // Set the run mode and initialize some run parameters
    if (rst_file_num <= 0) run_param::setRunStat(run_param::RunStatus::InitialRun);
    else run_param::setRunStat(run_param::RunStatus::RestartRun);
    run_param::init(static_cast<std::size_t>(my_rank));

    // Make a directory
    char dir_name[1024];
    sprintf(dir_name,"./result");
    makeOutputDirectory(dir_name);

    // Open a file to record debug information 
    dbg_utils::fopen();

    // Make instances of ParticleSystem and initialize them
    PS::ParticleSystem<FP_nbody> psys_nbody;
    PS::ParticleSystem<FP_star> psys_star;
    PS::ParticleSystem<FP_gas> psys_gas;
    psys_nbody.initialize();
    psys_star.initialize();
    psys_gas.initialize();

    // Make an instance of DomainInfo and initialize it
    PS::DomainInfo dinfo;
    dinfo.initialize();

    // Prepare the initial data
    switch (run_param::run_stat) {
        //======================
        //   Initial Run case
        //======================
        case run_param::RunStatus::InitialRun:
            // Make an initial condition at MPI rank 0.
            if (PS::Comm::getRank() == 0) {
                GalaxyIC(psys_nbody, psys_star, psys_gas);
            }
            else {
                psys_nbody.setNumberOfParticleLocal(0);
                psys_star.setNumberOfParticleLocal(0);
                psys_gas.setNumberOfParticleLocal(0);
            }

            // Broadcast run parameters
            run_param::setup();

            // Set the boundary condition and the size of the computational domain if needed.
            dinfo.setBoundaryCondition(run_param::basic::bc);
            if (run_param::basic::bc != PS::BOUNDARY_CONDITION_OPEN) {
                dinfo.setPosRootDomain(run_param::basic::pos_root_domain.low_,
                                       run_param::basic::pos_root_domain.high_);
            }

            // Perform domain decomposition 
            dinfo.collectSampleParticle(psys_nbody);
            dinfo.collectSampleParticle(psys_star,false);
            dinfo.collectSampleParticle(psys_gas,false);
            dinfo.decomposeDomain();

            // Perform particle exchange
            psys_nbody.exchangeParticle(dinfo);
            psys_star.exchangeParticle(dinfo);
            psys_gas.exchangeParticle(dinfo);

            break;
        //======================
        //   Restart Run case
        //======================
        case run_param::RunStatus::RestartRun:
            // Read files of run-parameters for restart
            {
                const std::string str_dir_name = dir_name;
                std::stringstream ss;
                ss << "RST_fn" << std::setw(5) << std::setfill('0') << rst_file_num;
                const std::string file_id = ss.str();
                run_param::setup(str_dir_name, file_id);
            }
            // Read a file of particle data for restart
            {
                char basename[256], format[256];
                sprintf(format,"%%s_np%%05d_r%%05d.dat");
                sprintf(basename, "result/nbody_RST_fn%05d",rst_file_num);
                psys_nbody.readParticleBinary(static_cast<const char * const>(basename),
                                              static_cast<const char * const>(format),
                                              &FP_nbody::readRestartData);
                //--- test(start) ---
                //TestHeader header_nbody;
                //psys_nbody.readParticleAscii(static_cast<const char * const>(basename),
                //                             static_cast<const char * const>(format),
                //                             header_nbody,
                //                             &FP_nbody::readAscii);
                //--- test(end) ---
                sprintf(basename, "result/star_RST_fn%05d",rst_file_num);
                psys_star.readParticleBinary(static_cast<const char * const>(basename),
                                             static_cast<const char * const>(format),
                                             &FP_star::readRestartData);
                sprintf(basename, "result/gas_RST_fn%05d",rst_file_num);
                psys_gas.readParticleBinary(static_cast<const char * const>(basename),
                                            static_cast<const char * const>(format),
                                            &FP_gas::readRestartData);
            }
            // Read a file of FDPS internal data for restart
            PS::F64ort pos_my_domain;
            {
                std::stringstream ss;
                ss << "result/fdps_int_RST_fn"
                   << std::setw(5) << std::setfill('0') << run_param::io::ndump_rst
                   << "_r"
                   << std::setw(5) << std::setfill('0') << my_rank
                   << ".dat";
                const std::string filename = ss.str();
                std::ifstream ifs;
                ifs.open(filename.c_str(), std::ios::in | std::ios::binary);
                ifs.read((char *)&pos_my_domain, sizeof(pos_my_domain));
                ifs.close();
            }

            // Set the boundary condition and the size of the computational domain if needed.
            dinfo.setBoundaryCondition(run_param::basic::bc);
            if (run_param::basic::bc != PS::BOUNDARY_CONDITION_OPEN) {
                dinfo.setPosRootDomain(run_param::basic::pos_root_domain.low_,
                                       run_param::basic::pos_root_domain.high_);
            }

            // Set domain decomposition information
            dinfo.decomposeDomainAll(const_cast<const PS::F64ort &>(pos_my_domain));

            // check
            //for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            //    if (psys_gas[i].id == 2606203)
            //        psys_gas[i].dump("[aft]", i, __func__, dbg_utils::fout);
            //}
            //PS::Finalize(); std::exit(0);

            break;
    }

    // Make tree structures
    const PS::S64 numPtclSPH  = std::max(psys_gas.getNumberOfParticleLocal(),1);
    const PS::S64 numPtclStar = std::max(psys_star.getNumberOfParticleLocal(), 1);
    const PS::S64 numPtclAll = psys_nbody.getNumberOfParticleLocal() + numPtclStar + numPtclSPH;

    const PS::F32 theta_grav = 0.5;
    PS::TreeForForceLong<Force_grav, EP_grav, EP_grav>::Monopole tree_grav;
    tree_grav.initialize(3 * numPtclAll, theta_grav);

    PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather tree_dens;
    tree_dens.initialize(3 * numPtclSPH);

    PS::TreeForForceShort<Force_hydro, EP_hydro, EP_hydro>::Symmetry tree_hydro;
    tree_hydro.initialize(3 * numPtclSPH);

    // Initialize the Phantom-GRAPE library and set the gravitational softening
#if defined(ENABLE_PHANTOM_GRAPE_X86)
    g5_open();
    g5_set_eps_to_all(run_param::sim::eps_grav);
#else
    EP_grav::eps = run_param::sim::eps_grav;
#endif

    // Peform force calculations 
    //- Gravity calculations
    tree_grav.setParticleLocalTree(psys_nbody);
    tree_grav.setParticleLocalTree(psys_star,false);
    tree_grav.setParticleLocalTree(psys_gas,false);
    tree_grav.calcForceMakingTree(CalcGravity<EP_grav>(),
                                  CalcGravity<PS::SPJMonopole>(),
                                  dinfo);
    for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
        psys_nbody[i].copyFromForce(tree_grav.getForce(i));
    }
    PS::S32 offset = psys_nbody.getNumberOfParticleLocal();
    for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
        psys_star[i].copyFromForce(tree_grav.getForce(i+offset));
    }
    offset += psys_star.getNumberOfParticleLocal();
    for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
        psys_gas[i].copyFromForce(tree_grav.getForce(i+offset));
    }

    //- SPH calculations
    if (run_param::isInitialRun()) {
        calcDensity(psys_gas, psys_star, dinfo, tree_dens, 0.0, 0.0);
    }
    // [Important]
    //    Note that we do not need to perform the densify calculation when
    //    RestartRun mode because the restart file contains density.
#if defined(USE_ENTROPY)
    setEntropy(psys_gas);
#endif
    setPressure(psys_gas);
    tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), psys_gas, dinfo);

    // Set the initial velocity of gas particles
    if (run_param::isInitialRun()) {
        setCircularVelocity(psys_gas);
    }

    // Get timestep
    PS::F64 dt = getTimeStep(psys_nbody, psys_star, psys_gas);

    // Calculate energies 
    run_param::calcConservedQuantities(psys_nbody, psys_star, psys_gas);

    // Main loop for time integration
    const PS::F64 time_start_main_loop = PS::GetWtime();
    const PS::S32 log_intvl = 10;
    for (;
         run_param::basic::time < run_param::basic::time_end;
         run_param::basic::time += dt, run_param::basic::nstep++) {
        dbg_utils::nstep = run_param::basic::nstep;
        //dbg_utils::fout << "nstep = " << run_param::basic::nstep << std::endl;
        if (run_param::basic::nstep % log_intvl == 0 && PS::Comm::getRank() == 0) {
        //if (PS::Comm::getRank() == 0) {
            std::cout << "nstep = " << run_param::basic::nstep 
                      << " dt = " << dt 
                      << " time = " << run_param::basic::time 
                      << " time_end = " << run_param::basic::time_end
                      << std::endl;
        }

        // Leap frog: Initial Kick & Full Drift
        InitialKick(psys_nbody, dt);
        InitialKick(psys_star, dt);
        InitialKick(psys_gas, dt);
        FullDrift(psys_nbody, dt);
        FullDrift(psys_star, dt);
        FullDrift(psys_gas, dt);
        if (dinfo.getBoundaryCondition() != PS::BOUNDARY_CONDITION_OPEN) {
            psys_nbody.adjustPositionIntoRootDomain(dinfo);
            psys_star.adjustPositionIntoRootDomain(dinfo);
            psys_gas.adjustPositionIntoRootDomain(dinfo);
        }

        // Leap frog: Predict
        Predict(psys_gas, dt);

        // Perform domain decomposition again
        dinfo.collectSampleParticle(psys_nbody);
        dinfo.collectSampleParticle(psys_star,false);
        dinfo.collectSampleParticle(psys_gas,false);
        dinfo.decomposeDomain();

        // Exchange the particles between the (MPI) processes
        psys_nbody.exchangeParticle(dinfo);
        psys_star.exchangeParticle(dinfo);
        psys_gas.exchangeParticle(dinfo);

        // Update array index so that we can access FP using a neighbor list of EPJ
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) psys_gas[i].idx = i;

        // Peform force calculations
        PS::F64 t_start; 
        //- Gravity calculations
        PS::Comm::barrier(); t_start = PS::GetWtime();
        tree_grav.setParticleLocalTree(psys_nbody);
        tree_grav.setParticleLocalTree(psys_star,false);
        tree_grav.setParticleLocalTree(psys_gas,false);
        tree_grav.calcForceMakingTree(CalcGravity<EP_grav>(),
                                      CalcGravity<PS::SPJMonopole>(),
                                      dinfo);
        for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
            psys_nbody[i].copyFromForce(tree_grav.getForce(i));
        }
        PS::S32 offset = psys_nbody.getNumberOfParticleLocal();
        for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
            psys_star[i].copyFromForce(tree_grav.getForce(i+offset));
        }
        offset += psys_star.getNumberOfParticleLocal();
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            psys_gas[i].copyFromForce(tree_grav.getForce(i+offset));
        }
        PS::Comm::barrier();
        if (run_param::basic::nstep % log_intvl == 0 && PS::Comm::getRank() == 0)
            std::cout << "t_grav (1) = " << (PS::GetWtime() - t_start) << std::endl;
        //- SPH calculations
        PS::Comm::barrier(); t_start = PS::GetWtime();
        calcDensity(psys_gas, psys_star, dinfo, tree_dens,
                    run_param::basic::time, dt);
        setPressure(psys_gas);
        tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), psys_gas, dinfo);
        PS::Comm::barrier();
        if (run_param::basic::nstep % log_intvl == 0 && PS::Comm::getRank() == 0)
            std::cout << "t_hydro (1) = " << (PS::GetWtime() - t_start) << std::endl;

        // Leap frog: Final Kick
        FinalKick(psys_nbody, dt);
        FinalKick(psys_star, dt);
        FinalKick(psys_gas, dt);

        // Syncronize thermodynamic quantities
#if defined(ENABLE_STELLAR_FEEDBACK) || \
    defined(ENABLE_COOLING_HEATING) || \
    defined(ENABLE_STAR_FORMATION)
        setPressure(psys_gas);
#endif

        // SF feedback calculation
#if defined(ENABLE_STELLAR_FEEDBACK)
        const bool flag_fb = StellarFeedback(psys_star, psys_gas, tree_dens,
                                             run_param::basic::time, dt);
        if (flag_fb) {
           calcDensity(psys_gas, psys_star, dinfo, tree_dens, 0.0, 0.0);
#if defined(USE_ENTROPY)
           setEntropy(psys_gas);
#endif
           setPressure(psys_gas);
        }
#endif

        // Cooling & Heating calculations
#if defined(ENABLE_COOLING_HEATING)
        CoolingHeating(psys_gas, dt);
#endif

        // Star formation
#if defined(ENABLE_STAR_FORMATION)
        const bool flag_sf = StarFormation(psys_gas, psys_star,
                                           run_param::basic::time, dt);
        if (flag_sf) {
           calcDensity(psys_gas, psys_star, dinfo, tree_dens, 0.0, 0.0);
#if defined(USE_ENTROPY)
           setEntropy(psys_gas);
#endif
           setPressure(psys_gas);
        }
#endif

        // Output data for restart
        if (run_param::basic::time > run_param::io::time_dump_rst){
        //if (run_param::basic::nstep == 10){
            // Output run-parameters files for restart
            {
                // (i) Update some variables temporarily because the restart
                //     simulation starts from the next step.
                const PS::S64 nstep = run_param::basic::nstep;
                const PS::F64 time = run_param::basic::time;
                const PS::S32 ndump = run_param::io::ndump;
                const PS::F64 time_dump = run_param::io::time_dump;
                const PS::S32 ndump_rst = run_param::io::ndump_rst;
                const PS::F64 time_dump_rst = run_param::io::time_dump_rst;
                run_param::basic::nstep++;
                run_param::basic::time += dt;
                if (run_param::basic::time > run_param::io::time_dump) {
                    // In this case, normal output will be done in this step.
                    run_param::io::ndump++;
                    run_param::io::time_dump += run_param::io::dt_dump; 
                }
                run_param::io::ndump_rst++;
                run_param::io::time_dump_rst += run_param::io::dt_dump_rst;
                // (ii) Output
                const std::string str_dir_name = dir_name;
                std::stringstream ss;
                ss << "RST_fn" << std::setw(5) << std::setfill('0') << ndump_rst;
                const std::string file_id = ss.str();
                run_param::writeFile(dir_name, file_id);
                // (iii) Restore to the current values
                run_param::basic::nstep = nstep;
                run_param::basic::time = time;
                run_param::io::ndump = ndump;
                run_param::io::time_dump = time_dump;
                run_param::io::ndump_rst = ndump_rst;
                run_param::io::time_dump_rst = time_dump_rst;
            }
            // Output particle data for restart
            {
                char basename[256], format[256];
                sprintf(format,"%%s_np%%05d_r%%05d.dat");
                sprintf(basename, "result/nbody_RST_fn%05d", run_param::io::ndump_rst);
                psys_nbody.writeParticleBinary(static_cast<const char * const>(basename),
                                               static_cast<const char * const>(format),
                                               &FP_nbody::writeRestartData);
                //--- test(start) ---
                //TestHeader header_nbody;
                //header_nbody.time = run_param::basic::time;
                //psys_nbody.writeParticleAscii(static_cast<const char * const>(basename),
                //                              static_cast<const char * const>(format),
                //                              header_nbody,
                //                              &FP_nbody::writeAscii);
                //--- test(end) ---
                sprintf(basename, "result/star_RST_sn%05d", run_param::io::ndump_rst);
                psys_star.writeParticleBinary(static_cast<const char * const>(basename),
                                              static_cast<const char * const>(format),
                                              &FP_star::writeRestartData);
                sprintf(basename, "result/gas_RST_sn%05d", run_param::io::ndump_rst);
                psys_gas.writeParticleBinary(static_cast<const char * const>(basename),
                                             static_cast<const char * const>(format),
                                             &FP_gas::writeRestartData);
            }
            // Output FDPS internal data for restart
            {
                const PS::F64ort pos_my_domain = dinfo.getPosDomain(my_rank);
                std::stringstream ss;
                ss << "result/fdps_int_RST_fn"
                   << std::setw(5) << std::setfill('0') << run_param::io::ndump_rst
                   << "_r"
                   << std::setw(5) << std::setfill('0') << my_rank
                   << ".dat";
                const std::string filename = ss.str();
                std::ofstream ofs;
                ofs.open(filename.c_str(), std::ios::trunc | std::ios::binary);
                ofs.write((char *)&pos_my_domain, sizeof(pos_my_domain));
                ofs.close();
            }
            // Notification
            if (PS::Comm::getRank() == 0){
               std::cout << "============================================" << std::endl;
               std::cout << "output restart files at time = " << run_param::basic::time << std::endl;
               std::cout << "============================================" << std::endl;
            }
            run_param::io::time_dump_rst += run_param::io::dt_dump_rst;
            run_param::io::ndump_rst++;
        }
        // Output normal output files
        if (run_param::basic::time > run_param::io::time_dump){
            // Output run-parameters files
            {
                // (i) Update some variables temporarily because the restart
                //     simulation starts from the next step.
                const PS::S64 nstep = run_param::basic::nstep;
                const PS::F64 time = run_param::basic::time;
                const PS::S32 ndump = run_param::io::ndump;
                const PS::F64 time_dump = run_param::io::time_dump;
                const PS::S32 ndump_rst = run_param::io::ndump_rst;
                const PS::F64 time_dump_rst = run_param::io::time_dump_rst;
                run_param::basic::nstep++;
                run_param::basic::time += dt;
                run_param::io::ndump++;
                run_param::io::time_dump += run_param::io::dt_dump; 
                run_param::io::ndump_rst++;
                run_param::io::time_dump_rst += run_param::io::dt_dump_rst;
                // (ii) Output
                const std::string str_dir_name = dir_name;
                std::stringstream ss;
                ss << "fn" << std::setw(5) << std::setfill('0') << ndump_rst;
                const std::string file_id = ss.str();
                run_param::writeFile(dir_name, file_id);
                // (iii) Restore to the current values
                run_param::basic::nstep = nstep;
                run_param::basic::time = time;
                run_param::io::ndump = ndump;
                run_param::io::time_dump = time_dump;
                run_param::io::ndump_rst = ndump_rst;
                run_param::io::time_dump_rst = time_dump_rst;
            }
            // Output particle data
            {
                FileHeader header_nbody, header_star, header_gas;
                header_nbody.time    = run_param::basic::time;
                header_nbody.numPtcl = psys_nbody.getNumberOfParticleGlobal();
                header_star.time     = run_param::basic::time;
                header_star.numPtcl  = psys_star.getNumberOfParticleGlobal();
                header_gas.time      = run_param::basic::time;
                header_gas.numPtcl   = psys_gas.getNumberOfParticleGlobal();
                char filename[256];
                sprintf(filename, "result/nbody%05d.txt", run_param::io::ndump);
                psys_nbody.writeParticleAscii(filename, header_nbody);
                sprintf(filename, "result/star%05d.txt", run_param::io::ndump);
                psys_star.writeParticleAscii(filename, header_star);
                sprintf(filename, "result/gas%05d.txt", run_param::io::ndump);
                psys_gas.writeParticleAscii(filename, header_gas);
            }
            // Output FDPS internal data
            {
                const PS::F64ort pos_my_domain = dinfo.getPosDomain(my_rank);
                std::stringstream ss;
                ss << "result/fdps_int_fn"
                   << std::setw(5) << std::setfill('0') << run_param::io::ndump_rst
                   << "_r"
                   << std::setw(5) << std::setfill('0') << my_rank
                   << ".dat";
                const std::string filename = ss.str();
                std::ofstream ofs;
                ofs.open(filename.c_str(), std::ios::trunc | std::ios::binary);
                ofs.write((char *)&pos_my_domain, sizeof(pos_my_domain));
                ofs.close();
            }
            // Notification
            if (PS::Comm::getRank() == 0){
               std::cout << "============================================" << std::endl;
               std::cout << "output files at time = " << run_param::basic::time << std::endl;
               std::cout << "============================================" << std::endl;
            }
            run_param::io::time_dump += run_param::io::dt_dump;
            run_param::io::ndump++;
        }

        // Peform force calculations
        //- Gravity calculations
        PS::Comm::barrier(); t_start = PS::GetWtime();
        tree_grav.setParticleLocalTree(psys_nbody);
        tree_grav.setParticleLocalTree(psys_star,false);
        tree_grav.setParticleLocalTree(psys_gas,false);
        tree_grav.calcForceMakingTree(CalcGravity<EP_grav>(),
                                      CalcGravity<PS::SPJMonopole>(),
                                      dinfo);
        for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
            psys_nbody[i].copyFromForce(tree_grav.getForce(i));
        }
        offset = psys_nbody.getNumberOfParticleLocal();
        for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
            psys_star[i].copyFromForce(tree_grav.getForce(i+offset));
        }
        offset += psys_star.getNumberOfParticleLocal();
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            psys_gas[i].copyFromForce(tree_grav.getForce(i+offset));
        }
        PS::Comm::barrier();
        if (run_param::basic::nstep % log_intvl == 0 && PS::Comm::getRank() == 0)
            std::cout << "t_grav (2) = " << (PS::GetWtime() - t_start) << std::endl;
        //- SPH calculations (density is already calculated)
        PS::Comm::barrier(); t_start = PS::GetWtime();
        tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), psys_gas, dinfo);
        PS::Comm::barrier();
        if (run_param::basic::nstep % log_intvl == 0 && PS::Comm::getRank() == 0)
            std::cout << "t_hydro (2) = " << (PS::GetWtime() - t_start) << std::endl;

        // Get a new timestep
        dt = getTimeStep(psys_nbody, psys_star, psys_gas);

        // Calculate energies
        if (run_param::basic::nstep % log_intvl == 0) 
            run_param::calcConservedQuantities(psys_nbody, psys_star, psys_gas);

        // For debug
        if (run_param::basic::nstep == 1) break;
        //if (run_param::basic::nstep == 2) break;
        //if (run_param::basic::nstep == 10) break;
        //if (run_param::basic::nstep == 11) break;
        //if (run_param::basic::nstep == 100) break;
    }
    const PS::F64 time_end_main_loop = PS::GetWtime();

    // Close the file to record debug information
    dbg_utils::fout << "elapsed time = " << time_end_main_loop - time_start_main_loop << " [s]" << std::endl;
    dbg_utils::fclose();

#if defined(ENABLE_PHANTOM_GRAPE_X86)
    g5_close();
#endif
    // Finalize FDPS
    PS::Finalize();
    return 0;
}

