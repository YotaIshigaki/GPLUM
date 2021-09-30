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
// plot-related headers
#include "plot_utilities.hpp"

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
    const std::string dir_name = "../../work-hokusai-bw/cooling_wo_fuv+sf+fb/result";
    //const PS::S32 file_num_bgn = 2;
    //const PS::S32 file_num_end = 115;
    const PS::S32 file_num_bgn = 100;
    const PS::S32 file_num_end = 100;
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

        // Make figure
        dislin_util::Plotter<PS::F64> plotter;
        const int nrows = 1;
        const int ncols = 2;
        const int size_panel[] = {1000, 1000};
        const int nx_max_poss = 1024; // the highest possible resolution in X
        const int ny_max_poss = 1024; // the highest possible resolution in Y
        const double npx_in_page_unit[2]
            = {static_cast<double>(nx_max_poss) / size_panel[0],
               static_cast<double>(ny_max_poss) / size_panel[1]};
#if 0
        const std::string save_fmt {"PS"};
        std::stringstream ss;
        ss << "fig" << std::setw(5) << std::setfill('0') << file_num << ".ps";
        const std::string filename = ss.str();
#else
        const std::string save_fmt {"PNG"};
        std::stringstream ss;
        ss << "fig" << std::setw(5) << std::setfill('0') << file_num << ".png";
        const std::string filename = ss.str();
#endif
        plotter.setPageLayoutAndOpenPage(nrows, ncols, size_panel,
                                         npx_in_page_unit, save_fmt, filename);

        PS::F64vec origin(0.0, 0.0, 0.0);
        const PS::F64 alpha = 0.0;
        const PS::F64 beta = 0.0;
        const PS::F64 gamma = 0.0;
        PS::F64ort2 bbox;
        bbox.low_.x = -0.1;
        bbox.low_.y = -0.1;
        bbox.high_.x = 0.1;
        bbox.high_.y = 0.1;
        const PS::S32 nx = nx_max_poss; 
        const PS::S32 ny = ny_max_poss; 
        plotter.setMesh(origin, alpha, beta, gamma, bbox, nx, ny);

        for (PS::S32 panel_id = 0; panel_id < 2; panel_id++) {
            dislin_util::ScalarQuantity qty;
            switch (panel_id) {
                case 0:
                    qty = dislin_util::ScalarQuantity::NumberDensity;
                    break;
                case 1:
                    qty = dislin_util::ScalarQuantity::Temperature;
                    break;
                default:
                    assert(false); 
            }
            plotter.setPanelConfig(panel_id, qty);
            plotter.setupPlotData(psys_gas);
            const std::string panel_text = "";
            plotter.drawPanel(panel_id, panel_text);
        }

        plotter.closePage();

    }

    // Finalize FDPS
    PS::Finalize();

    return 0;
}
