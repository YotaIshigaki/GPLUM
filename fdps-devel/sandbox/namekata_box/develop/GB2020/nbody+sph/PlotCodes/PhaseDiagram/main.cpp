/* C++ headers */
#include "common.h"
/* FDPS header */
#include <particle_simulator.hpp>
/* DISLIN header */
#include "discpp.h"
/* User-defined headers */
#include "macro_defs.h"
#include "mathematical_constants.h"
#include "physical_constants.h"
#include "run_parameters.h"
#include "user_defined.h"
#include "leapfrog.h"
#include "cooling_heating.h"
#include "star_formation.h"
#include "stellar_feedback.h"
#include "io.h"
#include "hydrodynamics.h"
#include "timestep.h"

int main(int argc, char *argv[]) {
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

    // Initialize run parameters
    run_param::setRunStat(run_param::RunStatus::RestartRun);
    run_param::init(static_cast<std::size_t>(my_rank));

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

    // Make figures
    //const std::string dir_name = "../../work-hokusai-bw/cooling_wo_fuv/result";
    //const std::string dir_name = "../../work-hokusai-bw/cooling_w_fuv/result";
    //const std::string dir_name = "../../work-hokusai-bw/cooling_wo_fuv+sf/result";
    const std::string dir_name = "../../work-hokusai-bw/cooling_wo_fuv+sf+fb/result";
    const PS::S32 file_num_bgn = 80;
    const PS::S32 file_num_end = 80;
    for (PS::S32 file_num = file_num_bgn; file_num <= file_num_end; file_num++) {
        // Notification
        if (my_rank == 0) {
            std::cout << "Aanalyzing file number " << file_num << std::endl;
        }

        // Read simulation data
        // (i) Read files of run-parameters 
        {   
            std::stringstream ss;
            ss << "fn" << std::setw(5) << std::setfill('0') << file_num;
            const std::string file_id = ss.str();
            run_param::setup(dir_name, file_id, false);
        }
        // (ii) Read data of gas particle
        {
            FileHeader header;
            char filename[256];
            sprintf(filename, "%s/gas%05d.txt", dir_name.c_str(), file_num);
            psys_gas.readParticleAscii(filename, header);
            std::cout << "t     = " << header.time << std::endl;
            std::cout << "N_gas = " << header.numPtcl << std::endl;
        }

        // Perform domain decomposition
        dinfo.decomposeDomainAll(psys_gas);

        // Perform particle exchange
        psys_gas.exchangeParticle(dinfo);

        // SPH calculations
#if defined(USE_ENTROPY)
        setEntropy(psys_gas);
        // Necessary because partile data do not contain the data of entropy
        // and the thermodynamic variables are calculated via entropy.
#endif
        setPressure(psys_gas);

        // Calculate nH, T data
        const PS::S32 n_loc = psys_gas.getNumberOfParticleLocal();
        const PS::S32 n_glb = psys_gas.getNumberOfParticleGlobal();
        std::vector<PS::F64> nH_loc, T_loc;
        nH_loc.resize(n_loc);
        T_loc.resize(n_loc);
        for (PS::S32 i = 0; i < n_loc; i++) {
            nH_loc[i] = psys_gas[i].dens * run_param::unit::dens * run_param::ism::Xhydrogen / phys_const::Mhydrogen; 
            T_loc[i] = psys_gas[i].getTemperature() * run_param::unit::temp;
        }
        std::vector<PS::F64> nH_glb, T_glb;
        nH_glb.resize(n_glb);
        T_glb.resize(n_glb);
        int n_recv_tot;
        PS::Comm::gatherVAll(nH_loc.data(), n_loc, nH_glb.data(), n_recv_tot);
        PS::Comm::gatherVAll(T_loc.data(), n_loc, T_glb.data(), n_recv_tot);
        
        // Output
        if (my_rank == 0) {
            std::stringstream ss;
            ss << "phase_diagram" << std::setw(5) << std::setfill('0') << file_num << ".txt";
            const std::string filename = ss.str();
            std::ofstream ofs;
            ofs.open(filename.c_str(), std::ios::trunc);
            for (PS::S32 i = 0; i < n_glb; i++) {
                ofs << nH_glb[i] << "    " << T_glb[i] << std::endl;
            }
            ofs.close();
        }

    }

    // Finalize FDPS
    PS::Finalize();

    return 0;
}
