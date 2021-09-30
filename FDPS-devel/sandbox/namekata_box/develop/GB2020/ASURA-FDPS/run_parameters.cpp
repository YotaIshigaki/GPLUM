/* C++ headers */
#include "common.h"
/* FDPS header */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "physical_constants.h"
#include "user_defined.h"
#include "run_parameters.h"
/* CELib header */
#include "CELib.h"

namespace run_parameters {

    RunStatus run_stat;
    bool is_initialized {false};

    namespace pseudo_random_number_generator {
        std::mt19937_64 mt;
    }

    namespace unit_system {
        PS::F64 mass; // unit mass [g]
        PS::F64 leng; // unit length [cm]
        PS::F64 time; // unit time [s]
        PS::F64 velc; // unit velocity [cm/s]
        PS::F64 eng;  // unit energy [erg]
        PS::F64 dens; // unit density [g/cm^3]
        PS::F64 temp; // unit temperature [K]
        PS::F64 spen; // unit specific energy [erg/g]
    }

    namespace basic {
        PS::BOUNDARY_CONDITION bc; // boundary condition
        PS::F64ort pos_root_domain; // size of computational box (in the simulation unit)
        PS::S64 nstep {1}; // number of steps
        PS::F64 time {0.0}; // current time in the simulation (in the simulation unit)
        PS::F64 time_end; // end time of the simulation (in the simulation unit)
    }

    namespace IO {
        PS::S32 ndump {1}; // next file number of particle data
        PS::F64 time_dump; // next time of output of particle data (in the simulation unit)
        PS::F64 dt_dump; // time interval of output of particle data (in the simulation unit)

        PS::S32 ndump_rst {1}; // next file number of particle data for restart
        PS::F64 time_dump_rst; // next time of output of particle data for restart (in the simulation unit)
        PS::F64 dt_dump_rst; // time interval of output of particle data for restart (in the simulation unit)
    }

    namespace conserved_quantities {
        PS::F64 M_dm_ini {0.0}; // initial dark matter mass
        PS::F64 M_star_ini {0.0}; // initial stellar mass
        PS::F64 M_gas_ini {0.0}; // initial gas mass
        PS::F64 M_elms_ini[CELibYield_Number] = {0.0}; // initial elements' mass
        PS::F64 M_tot_ini {0.0}; // initial total mass
        PS::F64 E_kin_ini {0.0}; // initial kinetic energy
        PS::F64 E_pot_ini {0.0}; // initial potential energy
        PS::F64 E_th_ini {0.0}; // initial thermal energy
        PS::F64 E_tot_ini {0.0}; // initial total energy
        PS::F64vec P_ini(0.0); // initial linear momentum
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
        PS::F64 L_ini {0.0};
#else
        PS::F64vec L_ini(0.0); // initial angular momentum
#endif

        PS::F64 M_dm; // dark matter mass
        PS::F64 M_star; // stellar mass
        PS::F64 M_gas; // gas mass
        PS::F64 M_elms[CELibYield_Number]; // elements' mass
        PS::F64 M_tot; // total mass
        PS::F64 E_kin; // kinetic energy
        PS::F64 E_pot; // potential energy
        PS::F64 E_th; // thermal energy
        PS::F64 E_tot; // total energy
        PS::F64vec P; // momentum
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
        PS::F64 L;
#else
        PS::F64vec L; // angular momentum
#endif

        PS::F64 E_lost_via_ch {0.0}; // energy loss by cooling/heating
        PS::F64 E_gain_via_fb {0.0}; // energy gain by feedback
    }

    namespace gravity {
        namespace soft {
            PS::F64 theta {0.5}; // Opening angle criterion
            PS::S32 n_group_limit {64}; // Number of particles in group tree cell
            PS::S32 n_leaf_limit {8}; // Number of particles in leaf tree cell
            PS::F64 eps_dm; // Gravitational softening for DM particles
            PS::F64 eps_star; // Gravitational softening for star particles
            PS::F64 eps_gas; // Gravitational softening for gas particles
            PS::F64 dt_fid; // Fiducial time step for the Soft part of the P3T scheme (dt_max)
            PS::F64 dt_actl; // Actual time step for the Soft part of the P3T scheme
            PS::F64 CFL {0.3}; // the CFL number for dynamics (only used in the case of variable time step)
            /* PMMM parameters will be added here */

        }
        namespace hard {
            PS::F64 eps; // Gravitational softening between star cluster particles
            PS::F64 eta; // Coefficient for time step for the 4th order Hermite integrator
            PS::F64 eta_s; // Coefficient for the initial time step used in the the 4th order Hermite integrator
            PS::F64 dt_limit; // Limit for the time step in the Hard part calculation
        }
    }

    namespace SPH {
        PS::S32 n_group_limit {64};
        PS::S32 n_leaf_limit {8};
        PS::S32 n_jp_limit {4096};
        PS::F64 h2h_next {1.25}; // Conversion factor from the current smoothing length to the next smoothing length estimate
        //PS::F64 h2h_next {1.1}; // Conversion factor from the current smoothing length to the next smoothing length estimate
        //PS::F64 h2h_next {1.01}; // Conversion factor from the current smoothing length to the next smoothing length estimate
        PS::F64 CFL {0.3}; // the CFL number for hydrodynamics
        //PS::F64 CFL {0.2}; // the CFL number for hydrodynamics
        PS::F64 eps {1.0e-4}; // the tolerance parameter for the Lagragian constraint
        //PS::S32 N_ngb {16}; // Number of neighbor particles (4^{2})
        //PS::S32 N_ngb {25}; // Number of neighbor particles (5^{2})
        //PS::S32 N_ngb {32}; // Number of neighbor particles
        //PS::S32 N_ngb {36}; // Number of neighbor particles (6^{2})
        //PS::S32 N_ngb {49}; // Number of neighbor particles (7^{2})
        //PS::S32 N_ngb {64}; // Number of neighbor particles (8^{2}, 4^{3})
        //PS::S32 N_ngb {81}; // Number of neighbor particles (9^{2})
        //PS::S32 N_ngb {100}; // Number of neighbor particles (10^{2})
        //PS::S32 N_ngb {121}; // Number of neighbor particles (11^{2})
        PS::S32 N_ngb {125}; // Number of neighbor particles (5^{3})
        //PS::S32 N_ngb {128}; // Number of neighbor particles
        //PS::S32 N_ngb {216}; // Number of neighbor particles
        PS::S32 N_ngb_mgn {1}; // Margin in number of neighbor particles
        PS::S32 N_ngb_limit {1024}; // Maximum number of neighbor particles (only used if using a differentiable constraint)
        PS::F64 alpha_AV_ini {1.0}; // the initial value of \alpha in the artificial viscosity
        PS::F64 alpha_AV_min {0.1}; // the minimum value of \alpha in the artificial viscosity
        PS::F64 alpha_AV_max {3.0}; // the maximum value of \alpha in the artificial viscosity
        PS::F64 alpha2beta_AV {0.0}; // the conversion factor from \alpha to \beta in the artificial viscosity
        PS::F64 ell_AV {0.05}; // A parameter determing the decay timescale of \alpha in the artificial viscosity
                               // (see Eq(9) in Cullen & Dehnen (2010))
        PS::F64 T_lower_limit {1.0e1}; // The lower temperature limit. The unit is K.
        PS::F64 T_upper_limit {1.0e8}; // The upper temperature limit. The unit is K.
    }

    namespace ISM {
        // Solar mass abundance
        // Note that the following values are those of solar abundance
        // (the recommended value in Asplund et al.(2009)[ARA+A,47,481]).
        const PS::F64 Xhydrogen_solar = 0.7381;
        const PS::F64 Yhelium_solar = 0.2485;
        const PS::F64 Zmetal_solar = 1.0 - Xhydrogen_solar - Yhelium_solar;
        PS::F64 Zcarbon_solar;
        PS::F64 Znitrogen_solar;
        PS::F64 Zoxygen_solar;
        PS::F64 Zneon_solar;
        PS::F64 Zmagnesium_solar;
        PS::F64 Zsilicon_solar;
        PS::F64 Zsulfur_solar;
        PS::F64 Zcalcium_solar;
        PS::F64 Ziron_solar;
        PS::F64 Znickel_solar;
        PS::F64 Zeuropium_solar;

        // Mass abundance of ISM of interest
        PS::F64 Xhydrogen {Xhydrogen_solar};
        PS::F64 Yhelium {Yhelium_solar};
        PS::F64 Zmetal {Zmetal_solar};

        // Gas properties    
        PS::F64 gamma {5.0/3.0}; // specific heat ratio or heat capacity ratio
        //PS::F64 gamma {1.4}; // specific heat ratio or heat capacity ratio
        PS::F64 mu {0.6}; // mean molecular weight relative to the mass of hydrogen
   
        // FUV background
        PS::F64 epsilon_FUV {0.05};
        PS::F64 G0_FUV {1.7};
    }

    namespace star_formation {
        // Star formation criterions
        PS::F64 nH_threshold {1.0e2}; // [cm^{-3}]
        PS::F64 T_threshold {1.0e2}; // [K]
        PS::F64 C_star {0.1};
        PS::S32 f_spawn {3}; // mass ratio of m_{*,spawn} to m_{gas,0}; should be an integer.
    }

    namespace stellar_feedback {
        PS::S32 N_ngb {128}; // number of neighbor gas particles around a FB particle
        PS::S32 N_ngb_mgn {1}; // margin in number of neighbor gas particles around a FB particle
    }


    void setRunStat(const RunStatus val) {
        run_stat = val;
        is_initialized = true;
        if (PS::Comm::getRank() == 0) {
            switch (run_stat) {
            case RunStatus::InitialRun:
                std::cout << "run_param::run_stat is set to InitialRun." << std::endl;
                break;
            case RunStatus::RestartRun:
                std::cout << "run_param::run_stat is set to RestartRun." << std::endl;
                break;
            default:
                std::cout << "run_param::run_stat is set to an invalid value." << std::endl;
                assert(false);
            }
        }
    }

    bool isInitialRun() {
        assert(is_initialized);
        if (run_stat == RunStatus::InitialRun) return true;
        else return false;
    }

    bool isRestartRun() {
        assert(is_initialized);
        if (run_stat == RunStatus::RestartRun) return true;
        else return false;
    }

    void init(const std::size_t seed) {
        // This function initializes variables that cannot be
        // initialized at compile time.
        // section "prng"
        prng::mt.seed(seed);
        // section "ism"
        ism::Zcarbon_solar = std::pow(10.0, 8.43 - 12.0) * ism::Xhydrogen_solar * (phys_const::Mcarbon/phys_const::Mhydrogen);
        ism::Znitrogen_solar = std::pow(10.0, 7.83 - 12.0) * ism::Xhydrogen_solar * (phys_const::Mnitrogen/phys_const::Mhydrogen);
        ism::Zoxygen_solar = std::pow(10.0, 8.69 - 12.0) * ism::Xhydrogen_solar * (phys_const::Moxygen/phys_const::Mhydrogen);
        ism::Zneon_solar = std::pow(10.0, 7.93 - 12.0) * ism::Xhydrogen_solar * (phys_const::Mneon/phys_const::Mhydrogen);
        ism::Zmagnesium_solar = std::pow(10.0, 7.6 - 12.0) * ism::Xhydrogen_solar * (phys_const::Mmagnesium/phys_const::Mhydrogen);
        ism::Zsilicon_solar = std::pow(10.0, 7.51 - 12.0) * ism::Xhydrogen_solar * (phys_const::Msilicon/phys_const::Mhydrogen);
        ism::Zsulfur_solar = std::pow(10.0, 7.12 - 12.0) * ism::Xhydrogen_solar * (phys_const::Msulfur/phys_const::Mhydrogen);
        ism::Zcalcium_solar = std::pow(10.0, 6.34 - 12.0) * ism::Xhydrogen_solar * (phys_const::Mcalcium/phys_const::Mhydrogen);
        ism::Ziron_solar = std::pow(10.0, 7.5 - 12.0) * ism::Xhydrogen_solar * (phys_const::Miron/phys_const::Mhydrogen);
        ism::Znickel_solar = std::pow(10.0, 6.22 - 12.0) * ism::Xhydrogen_solar * (phys_const::Mnickel/phys_const::Mhydrogen);
        ism::Zeuropium_solar = std::pow(10.0, 0.52 - 12.0) * ism::Xhydrogen_solar * (phys_const::Meuropium/phys_const::Mhydrogen);
    }

    void writeRankDependentFile(const std::string dir_name,
                                const std::string file_id) {
        // write the internal state of pseudo-random number generator into a file
        {
            std::stringstream ss;
            ss << dir_name << "/run_param_prng_" << file_id
               << "_r" << std::setw(5) << std::setfill('0') << PS::Comm::getRank()
               << ".txt";
            const std::string filename = ss.str();
            std::ofstream ofs;
            ofs.open(filename.c_str(), std::ios::trunc);
            if (!ofs) {
                std::cout << "Could not open file " << filename << std::endl;
                assert(false);
            }
            ofs << prng::mt;
            ofs.close();
        }
    }

    void readRankDependentFile(const std::string dir_name,
                               const std::string file_id) {
        // read the internal state of pseudo-random number generator from a file
        {
            std::stringstream ss;
            ss << dir_name << "/run_param_prng_" << file_id
               << "_r" << std::setw(5) << std::setfill('0') << PS::Comm::getRank()
               << ".txt";
            const std::string filename = ss.str();
            std::ifstream ifs;
            ifs.open(filename.c_str(), std::ios::in);
            if (!ifs) {
                std::cout << "Could not open file " << filename << std::endl;
                assert(false);
            }
            ifs >> prng::mt;
            ifs.close();
        }
    }

    void writeRankSharedFile(const std::string dir_name,
                             const std::string file_id) {
        std::stringstream ss;
        ss << dir_name << "/run_param_common_" << file_id << ".dat";
        const std::string filename = ss.str();
        std::ofstream ofs;
        ofs.open(filename.c_str(), std::ios::trunc | std::ios::binary);
        if (!ofs) {
            std::cout << "Could not open file " << filename << std::endl;
            assert(false);
        }
        // section "unit"
        ofs.write((char *)&unit::mass, sizeof(unit::mass));
        ofs.write((char *)&unit::leng, sizeof(unit::leng));
        ofs.write((char *)&unit::time, sizeof(unit::time));
        // section "basic"
        ofs.write((char *)&basic::bc, sizeof(basic::bc));
        ofs.write((char *)&basic::pos_root_domain, sizeof(basic::pos_root_domain));
        ofs.write((char *)&basic::nstep, sizeof(basic::nstep));
        ofs.write((char *)&basic::time, sizeof(basic::time));
        ofs.write((char *)&basic::time_end, sizeof(basic::time_end));
        // section "io"
        ofs.write((char *)&io::ndump, sizeof(io::ndump));
        ofs.write((char *)&io::time_dump, sizeof(io::time_dump));
        ofs.write((char *)&io::dt_dump, sizeof(io::dt_dump));
        ofs.write((char *)&io::ndump_rst, sizeof(io::ndump_rst));
        ofs.write((char *)&io::time_dump_rst, sizeof(io::time_dump_rst));
        ofs.write((char *)&io::dt_dump_rst, sizeof(io::dt_dump_rst));
        // section "cq"
        ofs.write((char *)&cq::M_dm_ini, sizeof(cq::M_dm_ini));
        ofs.write((char *)&cq::M_star_ini, sizeof(cq::M_star_ini));
        ofs.write((char *)&cq::M_gas_ini, sizeof(cq::M_gas_ini));
        ofs.write((char *)&cq::M_elms_ini[0], CELibYield_Number * sizeof(cq::M_elms_ini[0]));
        ofs.write((char *)&cq::M_tot_ini, sizeof(cq::M_tot_ini));
        ofs.write((char *)&cq::E_kin_ini, sizeof(cq::E_kin_ini));
        ofs.write((char *)&cq::E_pot_ini, sizeof(cq::E_pot_ini));
        ofs.write((char *)&cq::E_th_ini, sizeof(cq::E_th_ini));
        ofs.write((char *)&cq::E_tot_ini, sizeof(cq::E_tot_ini));
        ofs.write((char *)&cq::P_ini, sizeof(cq::P_ini));
        ofs.write((char *)&cq::L_ini, sizeof(cq::L_ini));
        ofs.write((char *)&cq::E_lost_via_ch, sizeof(cq::E_lost_via_ch));
        ofs.write((char *)&cq::E_gain_via_fb, sizeof(cq::E_gain_via_fb));
        // section "grav"
        ofs.write((char *)&grav::soft::theta, sizeof(grav::soft::theta));
        ofs.write((char *)&grav::soft::n_group_limit, sizeof(grav::soft::n_group_limit));
        ofs.write((char *)&grav::soft::n_leaf_limit, sizeof(grav::soft::n_leaf_limit));
        ofs.write((char *)&grav::soft::eps_dm, sizeof(grav::soft::eps_dm));
        ofs.write((char *)&grav::soft::eps_star, sizeof(grav::soft::eps_star));
        ofs.write((char *)&grav::soft::eps_gas, sizeof(grav::soft::eps_gas));
        ofs.write((char *)&grav::soft::dt_fid, sizeof(grav::soft::dt_fid));
        ofs.write((char *)&grav::soft::CFL, sizeof(grav::soft::CFL));
        ofs.write((char *)&grav::hard::eps, sizeof(grav::hard::eps));
        ofs.write((char *)&grav::hard::eta, sizeof(grav::hard::eta));
        ofs.write((char *)&grav::hard::eta_s, sizeof(grav::hard::eta_s));
        ofs.write((char *)&grav::hard::dt_limit, sizeof(grav::hard::dt_limit));
        // section "sph"
        ofs.write((char *)&sph::n_group_limit, sizeof(sph::n_group_limit));
        ofs.write((char *)&sph::n_leaf_limit, sizeof(sph::n_leaf_limit));
        ofs.write((char *)&sph::n_jp_limit, sizeof(sph::n_jp_limit));
        ofs.write((char *)&sph::h2h_next, sizeof(sph::h2h_next));
        ofs.write((char *)&sph::CFL, sizeof(sph::CFL));
        ofs.write((char *)&sph::eps, sizeof(sph::eps));
        ofs.write((char *)&sph::N_ngb, sizeof(sph::N_ngb));
        ofs.write((char *)&sph::N_ngb_mgn, sizeof(sph::N_ngb_mgn));
        ofs.write((char *)&sph::N_ngb_limit, sizeof(sph::N_ngb_limit));
        ofs.write((char *)&sph::alpha_AV_ini, sizeof(sph::alpha_AV_ini));
        ofs.write((char *)&sph::alpha_AV_min, sizeof(sph::alpha_AV_min));
        ofs.write((char *)&sph::alpha_AV_max, sizeof(sph::alpha_AV_max));
        ofs.write((char *)&sph::alpha2beta_AV, sizeof(sph::alpha2beta_AV));
        ofs.write((char *)&sph::ell_AV, sizeof(sph::ell_AV));
        ofs.write((char *)&sph::T_lower_limit, sizeof(sph::T_lower_limit));
        ofs.write((char *)&sph::T_upper_limit, sizeof(sph::T_upper_limit));
        // section "ism"
        ofs.write((char *)&ism::Xhydrogen, sizeof(ism::Xhydrogen));
        ofs.write((char *)&ism::Yhelium, sizeof(ism::Yhelium));
        ofs.write((char *)&ism::Zmetal, sizeof(ism::Zmetal));
        ofs.write((char *)&ism::gamma, sizeof(ism::gamma));
        ofs.write((char *)&ism::mu, sizeof(ism::mu));
        ofs.write((char *)&ism::epsilon_FUV, sizeof(ism::epsilon_FUV));
        ofs.write((char *)&ism::G0_FUV, sizeof(ism::G0_FUV));
        // section "sf"
        ofs.write((char *)&sf::nH_threshold, sizeof(sf::nH_threshold));
        ofs.write((char *)&sf::T_threshold, sizeof(sf::T_threshold));
        ofs.write((char *)&sf::C_star, sizeof(sf::C_star));
        ofs.write((char *)&sf::f_spawn, sizeof(sf::f_spawn));
        // section "fb"
        ofs.write((char *)&fb::N_ngb, sizeof(fb::N_ngb));
        ofs.write((char *)&fb::N_ngb_mgn, sizeof(fb::N_ngb_mgn));
        ofs.close();
    }

    void readRankSharedFile(const std::string dir_name,
                            const std::string file_id) {
        std::stringstream ss;
        ss << dir_name << "/run_param_common_" << file_id << ".dat";
        const std::string filename = ss.str();
        std::ifstream ifs;
        ifs.open(filename.c_str(), std::ios::in | std::ios::binary);
        if (!ifs) {
            std::cout << "Could not open file " << filename << std::endl;
            assert(false);
        }
        // section "unit"
        ifs.read((char *)&unit::mass, sizeof(unit::mass));
        ifs.read((char *)&unit::leng, sizeof(unit::leng));
        ifs.read((char *)&unit::time, sizeof(unit::time));
        // section "basic"
        ifs.read((char *)&basic::bc, sizeof(basic::bc));
        ifs.read((char *)&basic::pos_root_domain, sizeof(basic::pos_root_domain));
        ifs.read((char *)&basic::nstep, sizeof(basic::nstep));
        ifs.read((char *)&basic::time, sizeof(basic::time));
        ifs.read((char *)&basic::time_end, sizeof(basic::time_end));
        // section "io"
        ifs.read((char *)&io::ndump, sizeof(io::ndump));
        ifs.read((char *)&io::time_dump, sizeof(io::time_dump));
        ifs.read((char *)&io::dt_dump, sizeof(io::dt_dump));
        ifs.read((char *)&io::ndump_rst, sizeof(io::ndump_rst));
        ifs.read((char *)&io::time_dump_rst, sizeof(io::time_dump_rst));
        ifs.read((char *)&io::dt_dump_rst, sizeof(io::dt_dump_rst));
        // section "cq"
        ifs.read((char *)&cq::M_dm_ini, sizeof(cq::M_dm_ini));
        ifs.read((char *)&cq::M_star_ini, sizeof(cq::M_star_ini));
        ifs.read((char *)&cq::M_gas_ini, sizeof(cq::M_gas_ini));
        ifs.read((char *)&cq::M_elms_ini[0], CELibYield_Number * sizeof(cq::M_elms_ini[0]));
        ifs.read((char *)&cq::M_tot_ini, sizeof(cq::M_tot_ini));
        ifs.read((char *)&cq::E_kin_ini, sizeof(cq::E_kin_ini));
        ifs.read((char *)&cq::E_pot_ini, sizeof(cq::E_pot_ini));
        ifs.read((char *)&cq::E_th_ini, sizeof(cq::E_th_ini));
        ifs.read((char *)&cq::E_tot_ini, sizeof(cq::E_tot_ini));
        ifs.read((char *)&cq::P_ini, sizeof(cq::P_ini));
        ifs.read((char *)&cq::L_ini, sizeof(cq::L_ini));
        ifs.read((char *)&cq::E_lost_via_ch, sizeof(cq::E_lost_via_ch));
        ifs.read((char *)&cq::E_gain_via_fb, sizeof(cq::E_gain_via_fb));
        // section "grav"
        ifs.read((char *)&grav::soft::theta, sizeof(grav::soft::theta));
        ifs.read((char *)&grav::soft::n_group_limit, sizeof(grav::soft::n_group_limit));
        ifs.read((char *)&grav::soft::n_leaf_limit, sizeof(grav::soft::n_leaf_limit));
        ifs.read((char *)&grav::soft::eps_dm, sizeof(grav::soft::eps_dm));
        ifs.read((char *)&grav::soft::eps_star, sizeof(grav::soft::eps_star));
        ifs.read((char *)&grav::soft::eps_gas, sizeof(grav::soft::eps_gas));
        ifs.read((char *)&grav::soft::dt_fid, sizeof(grav::soft::dt_fid));
        ifs.read((char *)&grav::soft::CFL, sizeof(grav::soft::CFL));
        ifs.read((char *)&grav::hard::eps, sizeof(grav::hard::eps));
        ifs.read((char *)&grav::hard::eta, sizeof(grav::hard::eta));
        ifs.read((char *)&grav::hard::eta_s, sizeof(grav::hard::eta_s));
        ifs.read((char *)&grav::hard::dt_limit, sizeof(grav::hard::dt_limit));
        // section "sph"
        ifs.read((char *)&sph::n_group_limit, sizeof(sph::n_group_limit));
        ifs.read((char *)&sph::n_leaf_limit, sizeof(sph::n_leaf_limit));
        ifs.read((char *)&sph::n_jp_limit, sizeof(sph::n_jp_limit));
        ifs.read((char *)&sph::h2h_next, sizeof(sph::h2h_next));
        ifs.read((char *)&sph::CFL, sizeof(sph::CFL));
        ifs.read((char *)&sph::eps, sizeof(sph::eps));
        ifs.read((char *)&sph::N_ngb, sizeof(sph::N_ngb));
        ifs.read((char *)&sph::N_ngb_mgn, sizeof(sph::N_ngb_mgn));
        ifs.read((char *)&sph::N_ngb_limit, sizeof(sph::N_ngb_limit));
        ifs.read((char *)&sph::alpha_AV_ini, sizeof(sph::alpha_AV_ini));
        ifs.read((char *)&sph::alpha_AV_min, sizeof(sph::alpha_AV_min));
        ifs.read((char *)&sph::alpha_AV_max, sizeof(sph::alpha_AV_max));
        ifs.read((char *)&sph::alpha2beta_AV, sizeof(sph::alpha2beta_AV));
        ifs.read((char *)&sph::ell_AV, sizeof(sph::ell_AV));
        ifs.read((char *)&sph::T_lower_limit, sizeof(sph::T_lower_limit));
        ifs.read((char *)&sph::T_upper_limit, sizeof(sph::T_upper_limit));
        // section "ism"
        ifs.read((char *)&ism::Xhydrogen, sizeof(ism::Xhydrogen));
        ifs.read((char *)&ism::Yhelium, sizeof(ism::Yhelium));
        ifs.read((char *)&ism::Zmetal, sizeof(ism::Zmetal));
        ifs.read((char *)&ism::gamma, sizeof(ism::gamma));
        ifs.read((char *)&ism::mu, sizeof(ism::mu));
        ifs.read((char *)&ism::epsilon_FUV, sizeof(ism::epsilon_FUV));
        ifs.read((char *)&ism::G0_FUV, sizeof(ism::G0_FUV));
        // section "sf"
        ifs.read((char *)&sf::nH_threshold, sizeof(sf::nH_threshold));
        ifs.read((char *)&sf::T_threshold, sizeof(sf::T_threshold));
        ifs.read((char *)&sf::C_star, sizeof(sf::C_star));
        ifs.read((char *)&sf::f_spawn, sizeof(sf::f_spawn));
        // section "fb"
        ifs.read((char *)&fb::N_ngb, sizeof(fb::N_ngb));
        ifs.read((char *)&fb::N_ngb_mgn, sizeof(fb::N_ngb_mgn));
        ifs.close();
    }

    void bcast(const int root) {
        // section "unit"
        PS::Comm::broadcast(&unit::mass, 1, root);
        PS::Comm::broadcast(&unit::leng, 1, root);
        PS::Comm::broadcast(&unit::time, 1, root);
        unit::velc = unit::leng / unit::time;
        unit::eng  = unit::mass * SQ(unit::velc);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
        unit::dens = unit::mass / SQ(unit::leng);
#else
        unit::dens = unit::mass / CUBE(unit::leng);
#endif
        unit::temp = (unit::eng/phys_const::kBoltz) * (phys_const::Mproton/unit::mass);
        unit::spen = unit::eng / unit::mass;
        // output simulation units
        if (PS::Comm::getRank() == 0) {
            std::cout << "unit::mass = " << unit::mass << " [g]"
                      << " = " << unit::mass/phys_const::Msolar << " [Msolar]"
                      << std::endl;
            std::cout << "unit::leng = " << unit::leng << " [cm]"
                      << " = " << unit::leng/phys_const::kpc << " [kpc]"
                      << std::endl;
            std::cout << "unit::time = " << unit::time << " [s]"
                      << " = " << unit::time/phys_const::Gyr << " [Gyr]"
                      << std::endl;
            std::cout << "unit::velc = " << unit::velc << " [cm/s]"
                      << " = " << unit::velc/phys_const::km << " [km/s]"
                      << std::endl;
            std::cout << "unit::eng = " << unit::eng << " [erg]"
                      << std::endl;
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
            std::cout << "unit::dens = " << unit::dens << " [g/cm^2]"
                      << " = " << unit::dens/phys_const::Mhydrogen << " [cm^-2]"
                      << std::endl;
#else
            std::cout << "unit::dens = " << unit::dens << " [g/cm^3]"
                      << " = " << unit::dens/phys_const::Mhydrogen << " [cm^-3]"
                      << std::endl;
#endif
            std::cout << "unit::temp = " << unit::temp << " [K]"
                      << std::endl;
            std::cout << "unit::spen = " << unit::spen << " [erg/g]"
                      << std::endl;
        }
        // section "basic"
        PS::Comm::broadcast(&basic::bc, 1, root);
        PS::Comm::broadcast(&basic::pos_root_domain, 1, root);
        PS::Comm::broadcast(&basic::nstep, 1, root);
        PS::Comm::broadcast(&basic::time, 1, root);
        PS::Comm::broadcast(&basic::time_end, 1, root);
        // section "io"
        PS::Comm::broadcast(&io::ndump, 1, root);
        PS::Comm::broadcast(&io::time_dump, 1, root);
        PS::Comm::broadcast(&io::dt_dump, 1, root);
        PS::Comm::broadcast(&io::ndump_rst, 1, root);
        PS::Comm::broadcast(&io::time_dump_rst, 1, root);
        PS::Comm::broadcast(&io::dt_dump_rst, 1, root);
        // section "cq"
        PS::Comm::broadcast(&cq::M_dm_ini, 1, root);
        PS::Comm::broadcast(&cq::M_star_ini, 1, root);
        PS::Comm::broadcast(&cq::M_gas_ini, 1, root);
        PS::Comm::broadcast(&cq::M_elms_ini[0], CELibYield_Number, root);
        PS::Comm::broadcast(&cq::M_tot_ini, 1, root);
        PS::Comm::broadcast(&cq::E_kin_ini, 1, root);
        PS::Comm::broadcast(&cq::E_pot_ini, 1, root);
        PS::Comm::broadcast(&cq::E_th_ini, 1, root);
        PS::Comm::broadcast(&cq::E_tot_ini, 1, root);
        PS::Comm::broadcast(&cq::P_ini, 1,root);
        PS::Comm::broadcast(&cq::L_ini, 1, root);
        PS::Comm::broadcast(&cq::E_lost_via_ch, 1, root);
        PS::Comm::broadcast(&cq::E_gain_via_fb, 1, root);
        // section "grav"
        PS::Comm::broadcast(&grav::soft::theta, 1, root);
        PS::Comm::broadcast(&grav::soft::n_group_limit, 1, root);
        PS::Comm::broadcast(&grav::soft::n_leaf_limit, 1, root);
        PS::Comm::broadcast(&grav::soft::eps_dm, 1, root);
        PS::Comm::broadcast(&grav::soft::eps_star, 1, root);
        PS::Comm::broadcast(&grav::soft::eps_gas, 1, root);
        PS::Comm::broadcast(&grav::soft::dt_fid, 1, root);
        PS::Comm::broadcast(&grav::soft::CFL, 1, root);
        PS::Comm::broadcast(&grav::hard::eps, 1, root);
        PS::Comm::broadcast(&grav::hard::eta, 1, root);
        PS::Comm::broadcast(&grav::hard::eta_s, 1, root);
        PS::Comm::broadcast(&grav::hard::dt_limit, 1, root);
        // section "sph"
        PS::Comm::broadcast(&sph::n_group_limit, 1, root);
        PS::Comm::broadcast(&sph::n_leaf_limit, 1, root);
        PS::Comm::broadcast(&sph::n_jp_limit, 1, root);
        PS::Comm::broadcast(&sph::h2h_next, 1, root);
        PS::Comm::broadcast(&sph::CFL, 1, root);
        PS::Comm::broadcast(&sph::eps, 1, root);
        PS::Comm::broadcast(&sph::N_ngb, 1, root);
        PS::Comm::broadcast(&sph::N_ngb_mgn, 1, root);
        PS::Comm::broadcast(&sph::N_ngb_limit, 1, root);
        PS::Comm::broadcast(&sph::alpha_AV_ini, 1, root);
        PS::Comm::broadcast(&sph::alpha_AV_min, 1, root);
        PS::Comm::broadcast(&sph::alpha_AV_max, 1, root);
        PS::Comm::broadcast(&sph::alpha2beta_AV, 1, root);
        PS::Comm::broadcast(&sph::ell_AV, 1, root);
        PS::Comm::broadcast(&sph::T_lower_limit, 1, root);
        PS::Comm::broadcast(&sph::T_upper_limit, 1, root);
        // section "ism"
        PS::Comm::broadcast(&ism::Xhydrogen, 1, root);
        PS::Comm::broadcast(&ism::Yhelium, 1, root);
        PS::Comm::broadcast(&ism::Zmetal, 1, root);
        PS::Comm::broadcast(&ism::gamma, 1, root);
        PS::Comm::broadcast(&ism::mu, 1, root);
        PS::Comm::broadcast(&ism::epsilon_FUV, 1, root);
        PS::Comm::broadcast(&ism::G0_FUV, 1, root);
        // section "sf"
        PS::Comm::broadcast(&sf::nH_threshold, 1, root);
        PS::Comm::broadcast(&sf::T_threshold, 1, root);
        PS::Comm::broadcast(&sf::C_star, 1, root);
        PS::Comm::broadcast(&sf::f_spawn, 1, root);
        // section "fb"
        PS::Comm::broadcast(&fb::N_ngb, 1, root);
        PS::Comm::broadcast(&fb::N_ngb_mgn, 1, root);
    }

    void setup(const std::string dir_name,
               const std::string file_id,
               const bool read_rank_dep_file) {
        const PS::S32 root = 0;
        switch(run_stat) {
            case RunStatus::InitialRun:
                bcast(root);
                break;
            case RunStatus::RestartRun:
                if (read_rank_dep_file) readRankDependentFile(dir_name, file_id);
                if (PS::Comm::getRank() == root) readRankSharedFile(dir_name, file_id);
                bcast(root);
                break;
            default:
                assert(false);
        }
    }

    void writeFile(const std::string dir_name,
                   const std::string file_id) {
        const PS::S32 root = 0;
        writeRankDependentFile(dir_name, file_id);
        if (PS::Comm::getRank() == root) writeRankSharedFile(dir_name, file_id);
    }

    void calcConservedQuantities(const PS::ParticleSystem<FP_dm>& psys_dm,
                                 const PS::ParticleSystem<FP_star>& psys_star,
                                 const PS::ParticleSystem<FP_gas>& psys_gas) {
        PS::F64 M_dm_loc {0.0};
        PS::F64 M_star_loc {0.0};
        PS::F64 M_gas_loc {0.0};
        PS::F64 M_elms_loc[CELibYield_Number] = {0.0};
        PS::F64 E_kin_loc {0.0};
        PS::F64 E_pot_loc {0.0};
        PS::F64 E_th_loc {0.0};
        PS::F64vec P_loc {0.0};
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
        PS::F64 L_loc {0.0};
#else
        PS::F64vec L_loc {0.0};
#endif
        for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
            const PS::F64 mass = psys_dm[i].mass;
            const PS::F64vec pos = psys_dm[i].pos;
            const PS::F64vec vel = psys_dm[i].vel;
            const PS::F64vec mom = mass * vel;
            M_dm_loc += mass;
            E_kin_loc += 0.5 * mass * SQ(vel);
#if defined(ASURA_FDPS_ENABLE_GRAVITY)
            E_pot_loc += 0.5 * mass * (psys_dm[i].pot + mass / run_param::grav::soft::eps_dm);
#endif
            P_loc += mom;
            L_loc += pos ^ mom;
        }
        for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
            const PS::F64 mass = psys_star[i].mass;
            const PS::F64vec pos = psys_star[i].pos;
            const PS::F64vec vel = psys_star[i].vel;
            const PS::F64vec mom = mass * vel;
            M_star_loc += mass;
            E_kin_loc += 0.5 * mass * SQ(vel);
#if defined(ASURA_FDPS_ENABLE_GRAVITY)
            E_pot_loc += 0.5 * mass * (psys_star[i].pot + mass / run_param::grav::soft::eps_star);
#endif
            P_loc += mom;
            L_loc += pos ^ mom;
        }
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            const PS::F64 mass = psys_gas[i].mass;
            const PS::F64vec pos = psys_gas[i].pos;
            const PS::F64vec vel = psys_gas[i].vel;
            const PS::F64vec mom = mass * vel;
            M_gas_loc += mass;
            for (PS::S32 k = 0; k < CELibYield_Number; k++)
                M_elms_loc[k] += mass * psys_gas[i].mabn[k];
            E_kin_loc += 0.5 * mass * SQ(vel);
#if defined(ASURA_FDPS_ENABLE_GRAVITY)
            E_pot_loc += 0.5 * mass * (psys_gas[i].pot_grav + mass / run_param::grav::soft::eps_gas);
#endif
            E_th_loc += mass * psys_gas[i].eng;
            P_loc += mom;
            L_loc += pos ^ mom;
        }
        cq::M_dm = PS::Comm::getSum(M_dm_loc);
        cq::M_star = PS::Comm::getSum(M_star_loc);
        cq::M_gas = PS::Comm::getSum(M_gas_loc);
        for (PS::S32 k = 0; k < CELibYield_Number; k++)
            cq::M_elms[k] = PS::Comm::getSum(M_elms_loc[k]);
        cq::M_tot = cq::M_dm + cq::M_star + cq::M_gas; 
        cq::E_kin = PS::Comm::getSum(E_kin_loc);
        cq::E_pot = PS::Comm::getSum(E_pot_loc);
        cq::E_th = PS::Comm::getSum(E_th_loc);
        cq::E_tot = cq::E_kin + cq::E_pot + cq::E_th;
        cq::P = PS::Comm::getSum(P_loc);
        cq::L = PS::Comm::getSum(L_loc);
    
        static bool first_call = true;
        if (first_call) {
            if (run_param::isInitialRun()) {
                cq::M_dm_ini = cq::M_dm;
                cq::M_star_ini = cq::M_star;
                cq::M_gas_ini = cq::M_gas;
                for (PS::S32 k = 0; k < CELibYield_Number; k++)
                    cq::M_elms_ini[k] = cq::M_elms[k];
                cq::M_tot_ini = cq::M_tot;
                cq::E_kin_ini = cq::E_kin;
                cq::E_pot_ini = cq::E_pot;
                cq::E_th_ini = cq::E_th;
                cq::E_tot_ini = cq::E_tot;
                cq::P_ini = cq::P;
                cq::L_ini = cq::L;
            }
            first_call = false;
        }
    
        if (PS::Comm::getRank() == 0){
            const std::string str_rel_err = "[rel.err.] ";
            const std::string str_abs_err = "[abs.err.] ";
            const PS::F64 M_err = std::fabs((cq::M_tot - cq::M_tot_ini)/cq::M_tot_ini);
            const PS::F64 E_err = std::fabs((cq::E_tot - cq::E_tot_ini)/cq::E_tot_ini);
            std::string P_err_type[3];
            PS::F64vec P_err;
            if (cq::P_ini.x != 0.0) {
                P_err_type[0] = str_rel_err;
                P_err.x = std::fabs((cq::P.x - cq::P_ini.x)/cq::P_ini.x);
            } else {
                P_err_type[0] = str_abs_err;
                P_err.x = std::fabs(cq::P.x);
            }
            if (cq::P_ini.y != 0.0) {
                P_err_type[1] = str_rel_err;
                P_err.y = std::fabs((cq::P.y - cq::P_ini.y)/cq::P_ini.y);
            } else {
                P_err_type[1] = str_abs_err;
                P_err.y = std::fabs(cq::P.y);
            }
#if !defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
            if (cq::P_ini.z != 0.0) {
                P_err_type[2] = str_rel_err;
                P_err.z = std::fabs((cq::P.z - cq::P_ini.z)/cq::P_ini.z);
            } else {
                P_err_type[2] = str_abs_err;
                P_err.z = std::fabs(cq::P.z);
            }
#endif
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
            std::string L_err_type;
            PS::F64 L_err;
            if (cq::L_ini != 0.0) {
                L_err_type = str_rel_err;
                L_err = std::fabs((cq::L - cq::L_ini)/cq::L_ini);
            } else {
                L_err_type = str_abs_err;
                L_err = std::fabs(cq::L);
            }
#else
            std::string L_err_type[3];
            PS::F64vec L_err;
            if (cq::L_ini.x != 0.0) {
                L_err_type[0] = str_rel_err;
                L_err.x = std::fabs((cq::L.x - cq::L_ini.x)/cq::L_ini.x);
            } else {
                L_err_type[0] = str_abs_err;
                L_err.x = std::fabs(cq::L.x);
            }
            if (cq::L_ini.y != 0.0) {
                L_err_type[1] = str_rel_err;
                L_err.y = std::fabs((cq::L.y - cq::L_ini.y)/cq::L_ini.y);
            } else {
                L_err_type[1] = str_abs_err;
                L_err.y = std::fabs(cq::L.y);
            }
            if (cq::L_ini.z != 0.0) {
                L_err_type[2] = str_rel_err;
                L_err.z = std::fabs((cq::L.z - cq::L_ini.z)/cq::L_ini.z);
            } else {
                L_err_type[2] = str_abs_err;
                L_err.z = std::fabs(cq::L.z);
            }
#endif

            std::cout << "------------------------------------------------------------------" << std::endl;
            std::cout << "Mass (dark matter)   = " << cq::M_dm << std::endl;
            std::cout << "Mass (star)          = " << cq::M_star << std::endl;
            std::cout << "Mass (gas)           = " << cq::M_gas << std::endl;
            std::cout << "Mass (total)         = " << cq::M_tot << " (" << M_err << ")" << std::endl;
            std::cout << "Energy (kinetic)     = " << cq::E_kin << std::endl;
            std::cout << "Energy (potential)   = " << cq::E_pot << std::endl;
            std::cout << "Energy (thermal)     = " << cq::E_th  << std::endl;
            std::cout << "Energy (total)       = " << cq::E_tot << " (" << E_err << ")" << std::endl;
            std::cout << "Momentum (x)         = " << cq::P.x << " (" << P_err_type[0] << P_err.x << ")" <<std::endl;
            std::cout << "Momentum (y)         = " << cq::P.y << " (" << P_err_type[1] << P_err.y << ")" << std::endl;
#if !defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
            std::cout << "Momentum (z)         = " << cq::P.z << " (" << P_err_type[2] << P_err.z << ")" << std::endl;
#endif
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
            std::cout << "Angular momentum (z) = " << cq::L   << " (" << L_err_type    << L_err   << ")" << std::endl;
#else
            std::cout << "Angular momentum (x) = " << cq::L.x << " (" << L_err_type[0] << L_err.x << ")" << std::endl;
            std::cout << "Angular momentum (y) = " << cq::L.y << " (" << L_err_type[1] << L_err.y << ")" << std::endl;
            std::cout << "Angular momentum (z) = " << cq::L.z << " (" << L_err_type[2] << L_err.z << ")" << std::endl;
#endif
            std::cout << "------------------------------------------------------------------" << std::endl;
        }
    }
}
