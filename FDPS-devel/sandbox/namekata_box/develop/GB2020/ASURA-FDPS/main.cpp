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
#include "timing_utilities.h"
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
#include "SPH_kernel.h"

#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_NORMAL_RUN
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

void setCircularVelocity(PS::ParticleSystem<FP_star>& psys){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        const PS::F64vec acc = psys[i].acc;
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

void CoordinateTransform(PS::ParticleSystem<FP_dm> & psys_dm,
                         PS::ParticleSystem<FP_star> & psys_star,
                         PS::ParticleSystem<FP_gas> & psys_gas) {
    const PS::F64 dz = - 0.1;
    const PS::F64 angle = 0.25 * math_const::pi;
    for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
        const PS::F64vec pos = psys_dm[i].pos;
        const PS::F64vec vel = psys_dm[i].vel;
        psys_dm[i].pos.y =   pos.y * std::cos(angle) - pos.z * std::sin(angle);
        psys_dm[i].pos.z = - pos.y * std::sin(angle) + pos.z * std::cos(angle) - dz;
        psys_dm[i].vel.y =   vel.y * std::cos(angle) - vel.z * std::sin(angle);
        psys_dm[i].vel.z = - vel.y * std::sin(angle) + vel.z * std::cos(angle);
    }
    for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
        const PS::F64vec pos = psys_star[i].pos;
        const PS::F64vec vel = psys_star[i].vel;
        psys_star[i].pos.y =   pos.y * std::cos(angle) - pos.z * std::sin(angle);
        psys_star[i].pos.z = - pos.y * std::sin(angle) + pos.z * std::cos(angle) - dz;
        psys_star[i].vel.y =   vel.y * std::cos(angle) - vel.z * std::sin(angle);
        psys_star[i].vel.z = - vel.y * std::sin(angle) + vel.z * std::cos(angle);
    }
    for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
        const PS::F64vec pos = psys_gas[i].pos;
        const PS::F64vec vel = psys_gas[i].vel;
        psys_gas[i].pos.y =   pos.y * std::cos(angle) - pos.z * std::sin(angle);
        psys_gas[i].pos.z = - pos.y * std::sin(angle) + pos.z * std::cos(angle) - dz;
        psys_gas[i].vel.y =   vel.y * std::cos(angle) - vel.z * std::sin(angle);
        psys_gas[i].vel.z = - vel.y * std::sin(angle) + vel.z * std::cos(angle);
    }
    if (PS::Comm::getRank() == 0) {
        std::cout << "A coordinate transfrom is performed." << std::endl;
    }
}

#endif

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

    // Check the profile of SPH kernels
    //outputSPHKernelProfile();
    //PS::Finalize(); std::exit(0);

    // Make instances of ParticleSystem and initialize them
    PS::ParticleSystem<FP_dm> psys_dm;
    PS::ParticleSystem<FP_star> psys_star;
    PS::ParticleSystem<FP_gas> psys_gas;
    psys_dm.initialize();
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
            // Make an initial condition 
            // (Run parameter must be set here)
#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_NORMAL_RUN
#if 0
            // Make an initial condition at MPI rank 0.
            if (PS::Comm::getRank() == 0) {
                GalaxyIC(psys_dm, psys_star, psys_gas);
            } else {
                psys_dm.setNumberOfParticleLocal(0);
                psys_star.setNumberOfParticleLocal(0);
                psys_gas.setNumberOfParticleLocal(0);
            }
#else
            // Make an initial condition at each rank
            GordonBellIC(psys_dm, psys_star, psys_gas, dinfo);
#endif
#else
            // Make an initial condition at MPI rank 0.
            if (PS::Comm::getRank() == 0) {
#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_GLASS_DATA_GENERATION_MODE
                GlassDataGenerationModeIC(psys_dm, psys_star, psys_gas);
#elif ASURA_FDPS_EXEC_MODE == ASURA_FDPS_PARTICLE_COLLISION_TEST
                ParticleCollisionTestIC(psys_dm, psys_star, psys_gas);
#elif ASURA_FDPS_EXEC_MODE == ASURA_FDPS_SHOCK_TUBE_TEST
                ShockTubeTestIC(psys_dm, psys_star, psys_gas);
#elif ASURA_FDPS_EXEC_MODE == ASURA_FDPS_SURFACE_TENSION_TEST
                SurfaceTensionTestIC(psys_dm, psys_star, psys_gas);
#elif ASURA_FDPS_EXEC_MODE == ASURA_FDPS_POINT_EXPLOSION_TEST
                PointExplosionTestIC(psys_dm, psys_star, psys_gas);
#else
#error The value of macro ASURA_FDPS_EXEC_MODE is wrong.
#endif
            } else {
                psys_dm.setNumberOfParticleLocal(0);
                psys_star.setNumberOfParticleLocal(0);
                psys_gas.setNumberOfParticleLocal(0);
            }
#endif

            // Check
#if 0
             {
                std::stringstream ss;
                ss << "pos_ini_dm_" << std::setw(5) << std::setfill('0') << PS::Comm::getRank() << ".txt";
                const std::string file_name = ss.str();
                std::ofstream ofs;
                ofs.open(file_name.c_str(), std::ios::trunc);
                for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
                    ofs << psys_dm[i].pos << std::endl;
                }
                ofs.close();
            }
            {
                std::stringstream ss;
                ss << "pos_ini_star_" << std::setw(5) << std::setfill('0') << PS::Comm::getRank() << ".txt";
                const std::string file_name = ss.str();
                std::ofstream ofs;
                ofs.open(file_name.c_str(), std::ios::trunc);
                for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
                    ofs << psys_star[i].pos << std::endl;
                }
                ofs.close();
            }
            {
                std::stringstream ss;
                ss << "pos_ini_gas_" << std::setw(5) << std::setfill('0') << PS::Comm::getRank() << ".txt";
                const std::string file_name = ss.str();
                std::ofstream ofs;
                ofs.open(file_name.c_str(), std::ios::trunc);
                for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
                    ofs << psys_gas[i].pos << std::endl;
                }
                ofs.close();
            }
            PS::Finalize(); std::exit(0);
#endif

            // Broadcast run parameters
            run_param::setup();

            // Set the boundary condition and the size of the computational domain if needed.
            dinfo.setBoundaryCondition(run_param::basic::bc);
            if (run_param::basic::bc != PS::BOUNDARY_CONDITION_OPEN) {
                dinfo.setPosRootDomain(run_param::basic::pos_root_domain.low_,
                                       run_param::basic::pos_root_domain.high_);
            }

            // Perform domain decomposition 
            bool clear_flag {true};
            dinfo.collectSampleParticle(psys_dm, clear_flag);
            if (psys_dm.getNumberOfParticleGlobal() > 0) clear_flag = false;
            dinfo.collectSampleParticle(psys_star, clear_flag);
            if (psys_star.getNumberOfParticleGlobal() > 0) clear_flag = false;
            dinfo.collectSampleParticle(psys_gas, clear_flag);
            dinfo.decomposeDomain();

            // Perform particle exchange
            psys_dm.exchangeParticle(dinfo);
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
                sprintf(basename, "result/dm_RST_fn%05d",rst_file_num);
                psys_dm.readParticleBinary(static_cast<const char * const>(basename),
                                           static_cast<const char * const>(format),
                                           &FP_dm::readRestartData);
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
                   << std::setw(5) << std::setfill('0') << rst_file_num
                   << "_r"
                   << std::setw(5) << std::setfill('0') << my_rank
                   << ".dat";
                const std::string filename = ss.str();
                std::ifstream ifs;
                ifs.open(filename.c_str(), std::ios::in | std::ios::binary);
                if (!ifs) {
                    std::cerr << "Could not open file " << filename << std::endl;
                    assert(false);
                }
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

            // Check if the local particles are contained in pos_my_domain
            for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
                if (pos_my_domain.doesNotContain(psys_dm[i].pos)) {
                    std::cout << "There is a particle that is outside of the local domain!" << std::endl;
                    std::cout << "rank = " << my_rank
                              << " type = DM"
                              << " id = " << psys_dm[i].id
                              << " idx = " << i
                              << std::endl;
                    std::cout << "pos = " << psys_dm[i].pos
                              << " pos_my_domain = " << pos_my_domain
                              << std::endl;
                    assert(false);
                }
            }
            for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
                if (pos_my_domain.doesNotContain(psys_star[i].pos)) {
                    std::cout << "There is a particle that is outside of the local domain!" << std::endl;
                    std::cout << "rank = " << my_rank
                              << " type = Star"
                              << " id = " << psys_star[i].id
                              << " idx = " << i
                              << std::endl;
                    std::cout << "pos = " << psys_star[i].pos
                              << " pos_my_domain = " << pos_my_domain
                              << std::endl;
                    assert(false);
                }
            }
            for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
                if (pos_my_domain.doesNotContain(psys_gas[i].pos)) {
                    std::cout << "There is a particle that is outside of the local domain!" << std::endl;
                    std::cout << "rank = " << my_rank
                              << " type = Gas"
                              << " id = " << psys_gas[i].id
                              << " idx = " << i
                              << std::endl;
                    std::cout << "pos = " << psys_gas[i].pos
                              << " pos_my_domain = " << pos_my_domain
                              << std::endl;
                    assert(false);
                }
            }

            break;
    }

    // Make instances of TreeForForce class
    const PS::S64 numPtclSPH  = std::max(psys_gas.getNumberOfParticleLocal(),1);
    const PS::S64 numPtclStar = std::max(psys_star.getNumberOfParticleLocal(), 1);
    const PS::S64 numPtclAll = psys_dm.getNumberOfParticleLocal() + numPtclStar + numPtclSPH;

#if defined(ASURA_FDPS_ENABLE_GRAVITY)
    PS::TreeForForce<PS::SEARCH_MODE_LONG,
                     Force_grav, EP_grav, EP_grav,
                     Mom_grav, Mom_grav, SPJ_grav> tree_grav;
    tree_grav.initialize(3 * numPtclAll,
                         run_param::grav::soft::theta,
                         run_param::grav::soft::n_leaf_limit,
                         run_param::grav::soft::n_group_limit);
#endif

    PS::TreeForForceShort<Force_knl_sz, EPI_knl_sz, EPJ_knl_sz>::Gather tree_knl_sz;
    tree_knl_sz.initialize(3 * numPtclSPH, 0.0,
                           run_param::sph::n_leaf_limit,
                           run_param::sph::n_group_limit);

    PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather tree_dens;
    tree_dens.initialize(3 * numPtclSPH, 0.0,
                         run_param::sph::n_leaf_limit,
                         run_param::sph::n_group_limit);

    PS::TreeForForceShort<Force_hydro, EP_hydro, EP_hydro>::Symmetry tree_hydro;
    tree_hydro.initialize(3 * numPtclSPH, 0.0,
                          run_param::sph::n_leaf_limit,
                          run_param::sph::n_group_limit);

    // Initialize the Phantom-GRAPE library and set the gravitational softening
#if defined(ENABLE_PHANTOM_GRAPE_X86)
#error Not supported yet.
    //g5_open();
    //g5_set_eps_to_all(run_param::sim::eps_grav);
#endif

    // Peform force calculations
    PS::F64 t_start, t_end;
    time_prof.clear();
#if defined(ASURA_FDPS_ENABLE_GRAVITY)
    //- Gravity calculations
    barrier(); t_start = PS::GetWtime();
    tree_grav.clearTimeProfile();
    tree_grav.setParticleLocalTree(psys_dm);
    tree_grav.setParticleLocalTree(psys_star,false);
    tree_grav.setParticleLocalTree(psys_gas,false);
    tree_grav.calcForceMakingTree(CalcGravity<EP_grav>(),
                                  CalcGravity<SPJ_grav>(),
                                  dinfo);
    for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
        psys_dm[i].copyFromForce(tree_grav.getForce(i));
    }
    PS::S32 offset = psys_dm.getNumberOfParticleLocal();
    for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
        psys_star[i].copyFromForce(tree_grav.getForce(i+offset));
    }
    offset += psys_star.getNumberOfParticleLocal();
    for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
        psys_gas[i].copyFromForce(tree_grav.getForce(i+offset));
    }
    barrier();
    time_prof.calc_gravity_1st = PS::GetWtime() - t_start;
#else
    for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
        psys_dm[i].clearGravitationalForce();
    }
    for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
        psys_star[i].clearGravitationalForce();
    }
    for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
        psys_gas[i].clearGravitationalForce();
    }
#endif
    //- SPH calculations
    if (run_param::isInitialRun()) {
        calcKernelSize(psys_gas, psys_star, dinfo, tree_knl_sz,
                       0.0, 0.0,
                       FunctionUseType::Main);
    }
    // [Important]
    //    Note that we do not need to perform the kernel size calculation
    //    when RestartRun mode because the restart file contains the kernel size.
    calcDensityAndPressure(psys_gas, psys_star, dinfo, tree_dens,
                           0.0, 0.0,
                           FunctionUseType::Main);
    calcSoundSpeed(psys_gas);
    tree_hydro.clearTimeProfile();
    tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), psys_gas, dinfo);
    time_prof.calc_hydro_force_1st = tree_hydro.getTimeProfile().getTotalTime();


#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_NORMAL_RUN
#if 0
    // Set the initial velocity of gas particles
    if (run_param::isInitialRun()) {
        setCircularVelocity(psys_gas);
        // Perform a coordinate transformation
#ifdef ASURA_FDPS_PERFORM_COORD_TRANS
        CoordinateTransform(psys_dm, psys_star, psys_gas);
        // Perform domain decomposition again
        bool clear_flag {true};
        dinfo.collectSampleParticle(psys_dm, clear_flag);
        if (psys_dm.getNumberOfParticleGlobal() > 0) clear_flag = false;
        dinfo.collectSampleParticle(psys_star, clear_flag);
        if (psys_star.getNumberOfParticleGlobal() > 0) clear_flag = false;
        dinfo.collectSampleParticle(psys_gas, clear_flag);
        dinfo.decomposeDomain();
        // Perform particle exchange again
        psys_dm.exchangeParticle(dinfo);
        psys_star.exchangeParticle(dinfo);
        psys_gas.exchangeParticle(dinfo);
        // Perform interaction calculation again
        time_prof.clear();
#if defined(ASURA_FDPS_ENABLE_GRAVITY)
        //- Gravity calculations
        tree_grav.clearTimeProfile();
        tree_grav.setParticleLocalTree(psys_dm);
        tree_grav.setParticleLocalTree(psys_star,false);
        tree_grav.setParticleLocalTree(psys_gas,false);
        tree_grav.calcForceMakingTree(CalcGravity<EP_grav>(),
                                      CalcGravity<SPJ_grav>(),
                                      dinfo);
        for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
            psys_dm[i].copyFromForce(tree_grav.getForce(i));
        }
        PS::S32 offset = psys_dm.getNumberOfParticleLocal();
        for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
            psys_star[i].copyFromForce(tree_grav.getForce(i+offset));
        }
        offset += psys_star.getNumberOfParticleLocal();
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            psys_gas[i].copyFromForce(tree_grav.getForce(i+offset));
        }
        time_prof.calc_gravity_1st = tree_grav.getTimeProfile().getTotalTime();
#else
        for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
            psys_dm[i].clearGravitationalForce();
        }
        for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
            psys_star[i].clearGravitationalForce();
        }
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            psys_gas[i].clearGravitationalForce();
        }
#endif
        //- SPH calculations
        calcKernelSize(psys_gas, psys_star, dinfo, tree_knl_sz,
                       0.0, 0.0,
                       FunctionUseType::Main);
        calcDensityAndPressure(psys_gas, psys_star, dinfo, tree_dens,
                               0.0, 0.0,
                               FunctionUseType::Main);
        calcSoundSpeed(psys_gas);
        tree_hydro.clearTimeProfile();
        tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), psys_gas, dinfo);
        time_prof.calc_hydro_force_1st = tree_hydro.getTimeProfile().getTotalTime();
#endif
    }
#else
    // Set the initial velocity of gas particles
    if (run_param::isInitialRun()) {
        setCircularVelocity(psys_gas);
        setCircularVelocity(psys_star);
    }
#endif
#endif

    // Get timestep
    PS::F64 dt = getTimeStep(psys_dm, psys_star, psys_gas);

    // Calculate energies 
    run_param::calcConservedQuantities(psys_dm, psys_star, psys_gas);

    // Main loop for time integration
    const PS::F64 time_start_main_loop = PS::GetWtime();
    const PS::S32 log_intvl = 10;
    for (;
         run_param::basic::time < run_param::basic::time_end;
         run_param::basic::time += dt, run_param::basic::nstep++) {
        dbg_utils::nstep = run_param::basic::nstep;
        //dbg_utils::fout << "nstep = " << run_param::basic::nstep << std::endl;
        //if (run_param::basic::nstep % log_intvl == 0 && PS::Comm::getRank() == 0) {
        if (PS::Comm::getRank() == 0) {
            std::cout << "nstep = " << run_param::basic::nstep 
                      << " dt = " << dt 
                      << " time = " << run_param::basic::time 
                      << " time_end = " << run_param::basic::time_end
                      << std::endl;
        }
        PS::Comm::barrier();
        const PS::F64 t_start_inside_main_loop = PS::GetWtime();

        // Leap frog: Initial Kick & Full Drift
        InitialKick(psys_dm, dt);
        InitialKick(psys_star, dt);
        InitialKick(psys_gas, dt);
        FullDrift(psys_dm, dt);
        FullDrift(psys_star, dt);
        FullDrift(psys_gas, dt);
        if (dinfo.getBoundaryCondition() != PS::BOUNDARY_CONDITION_OPEN) {
            psys_dm.adjustPositionIntoRootDomain(dinfo);
            psys_star.adjustPositionIntoRootDomain(dinfo);
            psys_gas.adjustPositionIntoRootDomain(dinfo);
        }

        // Leap frog: Predict
        Predict(psys_gas, dt);

        // Perform domain decomposition again
        bool clear_flag {true};
#if defined(ASURA_FDPS_ENABLE_COMPUTATIONAL_TIME_BASED_DD)
        const PS::F32 weight = time_prof.getTotalTime();
#else
        const PS::F32 weight = psys_dm.getNumberOfParticleLocal()
                             + psys_star.getNumberOfParticleLocal()
                             + psys_gas.getNumberOfParticleLocal();;
#endif
        dinfo.collectSampleParticle(psys_dm, clear_flag, weight);
        if (psys_dm.getNumberOfParticleGlobal() > 0) clear_flag = false;
        dinfo.collectSampleParticle(psys_star, clear_flag, weight);
        if (psys_star.getNumberOfParticleGlobal() > 0) clear_flag = false;
        dinfo.collectSampleParticle(psys_gas, clear_flag, weight);
        dinfo.decomposeDomain();

        // Exchange the particles between the (MPI) processes
        psys_dm.exchangeParticle(dinfo);
        psys_star.exchangeParticle(dinfo);
        psys_gas.exchangeParticle(dinfo);

        // Update array index so that we can access FP using a neighbor list of EPJ
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            psys_gas[i].idx = i;
        }

        // Peform force calculations
        time_prof.clear();
#if defined(ASURA_FDPS_ENABLE_GRAVITY)
        //- Gravity calculations
        tree_grav.clearTimeProfile(); // barrier
        tree_grav.setParticleLocalTree(psys_dm);
        tree_grav.setParticleLocalTree(psys_star,false);
        tree_grav.setParticleLocalTree(psys_gas,false);
        tree_grav.calcForceMakingTree(CalcGravity<EP_grav>(),
                                      CalcGravity<SPJ_grav>(),
                                      dinfo);
        for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
            psys_dm[i].copyFromForce(tree_grav.getForce(i));
        }
        PS::S32 offset = psys_dm.getNumberOfParticleLocal();
        for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
            psys_star[i].copyFromForce(tree_grav.getForce(i+offset));
        }
        offset += psys_star.getNumberOfParticleLocal();
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            psys_gas[i].copyFromForce(tree_grav.getForce(i+offset));
        }
        time_prof.calc_gravity_1st = tree_grav.getTimeProfile().getTotalTime();
#else
        for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
            psys_dm[i].clearGravitationalForce();
        }
        for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
            psys_star[i].clearGravitationalForce();
        }
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            psys_gas[i].clearGravitationalForce();
        }
#endif
        // barrier 
        //- SPH calculations
        predictFeedbackRadius(psys_star, dinfo, tree_hydro,
                              run_param::basic::time, dt);
        // barrier
        calcKernelSize(psys_gas, psys_star, dinfo, tree_knl_sz,
                       run_param::basic::time, dt,
                       FunctionUseType::Main);
        // barrier
        calcDensityAndPressure(psys_gas, psys_star, dinfo, tree_dens,
                               run_param::basic::time, dt,
                               FunctionUseType::Main);
        calcSoundSpeed(psys_gas);
        calcArtificialViscosity(psys_gas, dt);
        tree_hydro.clearTimeProfile();
        // barrier
        //if (run_param::basic::nstep == 11345) {
        //    tree_hydro.debug_flag_ = true;
        //} else {
        //    tree_hydro.debug_flag_ = false;
        //}
        tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), psys_gas, dinfo);
        time_prof.calc_hydro_force_1st = tree_hydro.getTimeProfile().getTotalTime();
        //if (run_param::basic::nstep == 11345) {
        //    PS::Finalize();
        //    std::exit(0);
        //}

        // Leap frog: Final Kick
        FinalKick(psys_dm, dt);
        FinalKick(psys_star, dt);
        FinalKick(psys_gas, dt);

        // Syncronize thermodynamic quantities
        // (The specific energy is updated at the last Kick operator.
        //  Hence, the pressure should be recalculated for consistency)
#if defined(ASURA_FDPS_ENABLE_STELLAR_FEEDBACK) || \
    defined(ASURA_FDPS_ENABLE_COOLING_HEATING) || \
    defined(ASURA_FDPS_ENABLE_STAR_FORMATION)
        calcDensityAndPressure(psys_gas, psys_star, dinfo, tree_dens,
                               run_param::basic::time, dt,
                               FunctionUseType::PostProcess);
        calcSoundSpeed(psys_gas);
#endif

        // SF feedback calculation
#if defined(ASURA_FDPS_ENABLE_STELLAR_FEEDBACK)
        const bool flag_fb = StellarFeedback(psys_star, psys_gas, tree_dens,
                                             run_param::basic::time, dt);
        if (flag_fb) {
            calcDensityAndPressure(psys_gas, psys_star, dinfo, tree_dens,
                                   0.0, 0.0,
                                   FunctionUseType::PostProcess);
            calcSoundSpeed(psys_gas);
        }
#endif

        // Cooling & Heating calculations
#if defined(ASURA_FDPS_ENABLE_COOLING_HEATING)
        CoolingHeating(psys_gas, dt);
        calcDensityAndPressure(psys_gas, psys_star, dinfo, tree_dens,
                               0.0, 0.0,
                               FunctionUseType::PostProcess);
        calcSoundSpeed(psys_gas);
#endif

        // Star formation
#if defined(ASURA_FDPS_ENABLE_STAR_FORMATION)
        const bool flag_sf = StarFormation(psys_gas, psys_star,
                                           run_param::basic::time, dt);
        if (flag_sf) {
            calcKernelSize(psys_gas, psys_star, dinfo, tree_knl_sz,
                           0.0, 0.0,
                           FunctionUseType::PostProcess);
            calcDensityAndPressure(psys_gas, psys_star, dinfo, tree_dens,
                                   0.0, 0.0,
                                   FunctionUseType::PostProcess);
            calcSoundSpeed(psys_gas);
        }
#endif

        // [Only for glass data generation mode]
#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_GLASS_DATA_GENERATION_MODE
        checkDensityFluctuation(psys_gas);
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
                if (time > run_param::io::time_dump) {
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
                run_param::writeFile(str_dir_name, file_id);
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
                sprintf(basename, "result/dm_RST_fn%05d", run_param::io::ndump_rst);
                psys_dm.writeParticleBinary(static_cast<const char * const>(basename),
                                            static_cast<const char * const>(format),
                                            &FP_dm::writeRestartData);
                sprintf(basename, "result/star_RST_fn%05d", run_param::io::ndump_rst);
                psys_star.writeParticleBinary(static_cast<const char * const>(basename),
                                              static_cast<const char * const>(format),
                                              &FP_star::writeRestartData);
                sprintf(basename, "result/gas_RST_fn%05d", run_param::io::ndump_rst);
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
                if (!ofs) {
                    std::cout << "Could not open file " << filename << std::endl;
                    assert(false);
                }
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
                FileHeader header_dm, header_star, header_gas;
                header_dm.time    = run_param::basic::time;
                header_dm.numPtcl = psys_dm.getNumberOfParticleGlobal();
                header_star.time     = run_param::basic::time;
                header_star.numPtcl  = psys_star.getNumberOfParticleGlobal();
                header_gas.time      = run_param::basic::time;
                header_gas.numPtcl   = psys_gas.getNumberOfParticleGlobal();
                char filename[256];
                sprintf(filename, "result/dm%05d.txt", run_param::io::ndump);
                psys_dm.writeParticleAscii(filename, header_dm);
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

#if defined(ASURA_FDPS_ENABLE_STELLAR_FEEDBACK) || \
    defined(ASURA_FDPS_ENABLE_COOLING_HEATING) || \
    defined(ASURA_FDPS_ENABLE_STAR_FORMATION)
        // Peform force calculations
        //- Gravity calculations
#if defined(ASURA_FDPS_ENABLE_GRAVITY)
        tree_grav.clearTimeProfile();
        tree_grav.setParticleLocalTree(psys_dm);
        tree_grav.setParticleLocalTree(psys_star,false);
        tree_grav.setParticleLocalTree(psys_gas,false);
        tree_grav.calcForceMakingTree(CalcGravity<EP_grav>(),
                                      CalcGravity<SPJ_grav>(),
                                      dinfo);
        for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
            psys_dm[i].copyFromForce(tree_grav.getForce(i));
        }
        offset = psys_dm.getNumberOfParticleLocal();
        for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
            psys_star[i].copyFromForce(tree_grav.getForce(i+offset));
        }
        offset += psys_star.getNumberOfParticleLocal();
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            psys_gas[i].copyFromForce(tree_grav.getForce(i+offset));
        }
        time_prof.calc_gravity_2nd += tree_grav.getTimeProfile().getTotalTime(); 
#else
        for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
            psys_dm[i].clearGravitationalForce();
        }
        for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
            psys_star[i].clearGravitationalForce();
        }
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            psys_gas[i].clearGravitationalForce();
        }
#endif
        //- SPH calculations (density is already calculated)
        tree_hydro.clearTimeProfile();
        tree_hydro.calcForceAllAndWriteBack(CalcHydroForce(), psys_gas, dinfo);
        time_prof.calc_hydro_force_2nd += tree_hydro.getTimeProfile().getTotalTime();
#endif

        // Output time_prof
        PS::Comm::barrier();
        //if (run_param::basic::nstep % log_intvl == 0 && PS::Comm::getRank() == 0) {
        if (PS::Comm::getRank() == 0) {
            std::cout << "---------------------------------------" << std::endl;
            std::cout << "Result of the time measurement:" << std::endl;
            std::cout << "main loop (w/ barrier) = " << PS::GetWtime() - t_start_inside_main_loop << std::endl;
            time_prof.dump();
            std::cout << "---------------------------------------" << std::endl;
        }
        if (run_param::basic::nstep == 11450) {
            dbg_utils::fout << "---------------------------------------" << std::endl;
            dbg_utils::fout << "Result of the time measurement:" << std::endl;
            dbg_utils::fout << "main loop (w/ barrier) = " << PS::GetWtime() - t_start_inside_main_loop << std::endl;
            time_prof.dump(dbg_utils::fout);
            dbg_utils::fout << "---------------------------------------" << std::endl;
        }

        // Get a new timestep
        dt = getTimeStep(psys_dm, psys_star, psys_gas);

        // Calculate energies
        if (run_param::basic::nstep % log_intvl == 0) 
            run_param::calcConservedQuantities(psys_dm, psys_star, psys_gas);

        // For debug
        //break;
        //if (run_param::basic::nstep == 1) break;
        //if (run_param::basic::nstep == 10) break;
        //if (run_param::basic::nstep == 100) break;
        //if (run_param::basic::nstep == 1000) break;
        //if (run_param::basic::nstep == 10000) break;
        if (run_param::basic::nstep == 11450) break;
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

