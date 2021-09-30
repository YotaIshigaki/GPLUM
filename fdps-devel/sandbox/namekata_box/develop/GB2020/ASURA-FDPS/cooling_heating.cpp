/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "mathematical_constants.h"
#include "physical_constants.h"
#include "debug_utilities.h"
#include "timing_utilities.h"
#include "run_parameters.h"
#include "user_defined.h"
#include "cooling_heating.h"
/* ASRCH header */
extern "C" {
    #include "ASRCH.h"
}

//#define DEBUG_PRINT_CALC_COOLING_HEATING

void CoolingHeating(PS::ParticleSystem<FP_gas> & psys,
                    const PS::F64 dt_normalized) {
    barrier();
    const PS::F64 t_start = PS::GetWtime();
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();

    // Make interpolation table
    ASRCHMakeRedshiftInterpolateData(0.0);

    for (PS::S32 i = 0; i < n_loc; i++) {
        // Set input values for ASRCH
        const double nH = psys[i].dens * run_param::unit::dens * run_param::ism::Xhydrogen / phys_const::Mhydrogen; 
        const double log10nH = std::log10(nH);
        const double T = psys[i].getTemperature() * run_param::unit::temp;
        const double log10T = std::log10(T);
        const double dt = dt_normalized * run_param::unit::time; // [s]
#ifdef ASURA_FDPS_ENABLE_FUV_HEATING
        const double FUVHeating = 1.0e-24 * run_param::ism::epsilon_FUV * run_param::ism::G0_FUV
                                * (run_param::ism::Zmetal/run_param::ism::Zmetal_solar) / nH;
        // [Note]
        //    The unit of FUVHeaiting shoule be [erg cm^3/s],
        //    which will be converted into [erg/g/s] by multiplying
        //    nH^{2}/rho in ASRCH.
#else
        const double FUVHeating = 0.0;
#endif
        
        // Compute effective cooling rate using ASRCH
#ifdef ASRCH_OLD_API
        const double dUdt = ASRCHGetCoolingHeatingValueExactSolverWithFUVHeating(log10nH, log10T, run_param::ism::Zmetal, dt, FUVHeating);
#else
        struct ASRCHStructInput input;
        input.LogDensity = log10nH
        input.LogTemperature = log10T;
        input.Metallicity = run_param::ism::Zmetal;
        input.dt = dt;
        input.FUVHeating = FUVHeating;
        input.SoftXrayScaling = 0.0;
        input.Ntot = 0.0;

        const double dUdt = ASRCHGetCoolingHeatingValueExactSolverWithFUVHeating(input);
#endif

        // Update thermal quantities
        const double Unow = psys[i].eng;
        double Unew = Unow + (dUdt * dt)/run_param::unit::spen; // Unew is normalized by the simulation units

        // Update thermodynamic variables such as pressure
        psys[i].eng = Unew;

        // Apply temperature limits
        psys[i].applyTemperatureLimits();

        // Nan check
#if defined(ENABLE_NAN_CHECK)
        if (std::isnan(psys[i].eng) ||
            std::isinf(psys[i].eng) ||
            (psys[i].eng <= 0.0) ||
            std::isnan(FUVHeating)  ||
            std::isinf(FUVHeating)) {
            psys[i].dump("[nan]", i, __func__, dbg_utils::fout);
            dbg_utils::fout << "---- Input and Output (ASRCH) ----" << std::endl;
            dbg_utils::fout << "nH   = " << nH << " [cm^-3]" << std::endl;
            dbg_utils::fout << "T    = " << T << " [K]" << std::endl;
            dbg_utils::fout << "dt   = " << dt << " [s]" << std::endl;
            dbg_utils::fout << "FUV  = " << FUVHeating << " [erg cm^3/s]" << std::endl;
            dbg_utils::fout << "dUdt = " << dUdt << std::endl;
            dbg_utils::fout << "Unow = " << Unow << std::endl;
            dbg_utils::fout << "Unew = " << Unew << std::endl;
            dbg_utils::fout << "----------------------------------" << std::endl;
#if defined(FORCE_QUIT_IN_NAN_CHECK)
            assert(false);
#endif
        }
#endif

    }

    // Some Tests
#if 0
    if (PS::Comm::getRank() == 0) {
        // Set input values for ASRCH
        const double nH = 299.689; 
        const double log10nH = std::log10(nH);
        const double T = 76319.6;
        const double log10T = std::log10(T);
        const double dt = 1.76372e12; // [s]
        const double FUVHeating = 0.0;
        
        // Compute effective cooling rate using ASRCH
        const double dUdt = ASRCHGetCoolingHeatingValueExactSolverWithFUVHeating(log10nH, log10T, run_param::ism::Zmetal, dt, FUVHeating);

        // Update thermal quantities
        //const double Unow = 0.0834644;
        const double Unow = (phys_const::kBoltz * T)/((run_param::ism::gamma - 1.0) * run_param::ism::mu * phys_const::Mproton) / run_param::unit::spen;
        double Unew = Unow + (dUdt * dt)/run_param::unit::spen; // Unew is normalized by the simulation units
        
        // Check
        std::cout << "---- Input and Output (ASRCH, test) ----" << std::endl;
        std::cout << "nH   = " << nH << " [cm^-3]" << std::endl;
        std::cout << "T    = " << T << " [K]" << std::endl;
        std::cout << "dt   = " << dt << " [s]" << std::endl;
        std::cout << "FUV  = " << FUVHeating << " [erg cm^3/s]" << std::endl;
        std::cout << "dUdt = " << dUdt << std::endl;
        std::cout << "Unow = " << Unow << std::endl;
        std::cout << "Unew = " << Unew << std::endl;
        std::cout << "----------------------------------" << std::endl;
    }
    //CalcEquilibriumState(run_param::ism::Zmetal,
    //                     run_param::ism::epsilon_FUV,
    //                     0.0 * run_param::ism::G0_FUV);
    PS::Finalize();
    std::exit(0);
#endif

    // Update time_prof
    barrier();
    time_prof.cooling += PS::GetWtime() - t_start;
}

void CalcEquilibriumState(const PS::F64 Zmetal,
                          const PS::F64 epsilon_FUV,
                          const PS::F64 G0_FUV) {
    const PS::S32 my_rank = PS::Comm::getRank();

    // Create a grid of number density
    const PS::S32 n = 4096;
    const PS::F64 log10nH_min = - 6.0;
    const PS::F64 log10nH_max = 6.0;
    const PS::F64 dlog10nH = (log10nH_max - log10nH_min)/(PS::F64)(n-1);
    std::vector<PS::F64> log10nH(n);
    for (PS::S32 i = 0; i < n; i++)
        log10nH[i] = log10nH_min + dlog10nH * (PS::F64)i;

    // Open output file at Rank 0
    std::ofstream ofs;
    if (my_rank == 0) ofs.open("equilibrium_state.txt", std::ios::trunc);
    
    // Make interpolation table
    ASRCHMakeRedshiftInterpolateData(0.0);

    // Calculate the equilibrium temperature for each density
    //const PS::S32 itermax = 8192;
    const PS::S32 itermax = 1024;
    //const PS::S32 itermax = 1;
    const PS::F64 Tini = 1.0e4;
    for (PS::S32 i = 0; i < n; i++) {
        // Set fixed input values for ASRCH
        PS::F64 T = Tini;
        const PS::F64 nH = std::pow(10.0, log10nH[i]);
        const PS::F64 dt = 100.0 * phys_const::Gyr / (PS::F64) itermax; 
        const PS::F64 FUVHeating = 1.0e-24 * epsilon_FUV * G0_FUV
                                  * (Zmetal/run_param::ism::Zmetal_solar) / nH;
        for (PS::S32 iter = 0; iter < itermax; iter++) {
            // Set the other input values for ASRCH
            const PS::F64 log10T = std::log10(T);

            // Compute effective cooling rate using ASRCH
#ifdef ASRCH_OLD_API
            const PS::F64 dUdt = ASRCHGetCoolingHeatingValueExactSolverWithFUVHeating(log10nH[i], log10T, Zmetal, dt, FUVHeating);
#else
            struct ASRCHStructInput input;
            input.LogDensity = log10nH[i]
            input.LogTemperature = log10T;
            input.Metallicity = Zmetal;
            input.dt = dt;
            input.FUVHeating = FUVHeating;
            input.SoftXrayScaling = 0.0;
            input.Ntot = 0.0;

            const PS::F64 dUdt = ASRCHGetCoolingHeatingValueExactSolverWithFUVHeating(input);
#endif

            // Update thermal quantities
            const PS::F64 Tprev = T;
            const PS::F64 Unow = (T/run_param::unit::temp) / (run_param::ism::mu * (run_param::ism::gamma - 1.0));
            const PS::F64 Unew = Unow + (dUdt * dt)/run_param::unit::spen;
            T = (run_param::ism::mu * (run_param::ism::gamma - 1.0)) * Unew * run_param::unit::temp;
        }

        // Output file at Rank 0
        const PS::F64 log10T = std::log10(T);
        if (my_rank == 0) ofs << log10nH[i] << "    " << log10T << "    " << std::endl;
        
    }

    // Close the output file at Rank 0
    if (my_rank == 0) ofs.close();

}

