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
        PS::F64vec L_ini(0.0); // initial angular momentum

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
        PS::F64vec L; // angular momentum

        PS::F64 E_lost_via_ch {0.0}; // energy loss by cooling
        PS::F64 E_gain_via_fb {0.0}; // energy gain by feedback
    }

    namespace simulation {
        PS::F64 eps_grav; // gravitational softening
        PS::F64 mass_avg; // average mass of SPH particles
        PS::F64 SCF_smth {1.0}; // scale factor for smoothing length
        PS::F64 dt_max {0.0}; // maximum allowance of timestep
        const PS::F64 CFL_dyn = 0.3; // coefficient used to limit a timestep in terms of dynamics 
        const PS::F64 CFL_hydro = 0.3; // CFL(Courant-Friedrichs-Lewy) number for hydrodynamics
        const PS::F64 alpha_AV = 1.0; // coefficient of artificial viscosity
        const PS::S32 N_neighbor = 50; // number of neighbor particles
        const PS::S32 N_neighbor_FB = 128; // number of neighbor gas particles around a FB particle
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
        const PS::F64 Xhydrogen = Xhydrogen_solar;
        const PS::F64 Yhelium = Yhelium_solar;
        const PS::F64 Zmetal = Zmetal_solar;

        // Gas properties    
        const PS::F64 gamma = 5.0/3.0; // specific heat ratio or heat capacity ratio
        const PS::F64 mu = 0.6; // mean molecular weight relative to the mass of hydrogen
   
        // FUV background
        const PS::F64 epsilon_FUV = 0.05;
        const PS::F64 G0_FUV = 1.7;
    }

    namespace star_formation {
        // Star formation criterions
        const PS::F64 nH_threshold = 1.0e2; // [cm^{-3}]
        const PS::F64 T_threshold = 1.0e2; // [K]
        const PS::F64 C_star = 0.1;
        const PS::S32 f_spawn = 3; // mass ratio of m_{*,spawn} to m_{gas,0}; should be an integer.
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
        // section "sim"
        ofs.write((char *)&sim::eps_grav, sizeof(sim::eps_grav));
        ofs.write((char *)&sim::mass_avg, sizeof(sim::mass_avg));
        ofs.write((char *)&sim::dt_max, sizeof(sim::dt_max));
        ofs.close();
    }

    void readRankSharedFile(const std::string dir_name,
                            const std::string file_id) {
        std::stringstream ss;
        ss << dir_name << "/run_param_common_" << file_id << ".dat";
        const std::string filename = ss.str();
        std::ifstream ifs;
        ifs.open(filename.c_str(), std::ios::in | std::ios::binary);
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
        // section "sim"
        ifs.read((char *)&sim::eps_grav, sizeof(sim::eps_grav));
        ifs.read((char *)&sim::mass_avg, sizeof(sim::mass_avg));
        ifs.read((char *)&sim::dt_max, sizeof(sim::dt_max));
        ifs.close();
    }

    void bcast(const int root) {
        // section "unit"
        PS::Comm::broadcast(&unit::mass, 1, root);
        PS::Comm::broadcast(&unit::leng, 1, root);
        PS::Comm::broadcast(&unit::time, 1, root);
        unit::velc = unit::leng / unit::time;
        unit::eng  = unit::mass * SQ(unit::velc);
        unit::dens = unit::mass / CUBE(unit::leng);
        unit::temp = (unit::eng/phys_const::kBoltz) * (phys_const::Mhydrogen/unit::mass);
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
            std::cout << "unit::dens = " << unit::dens << " [g/cm^3]"
                      << " = " << unit::dens/phys_const::Mhydrogen << " [cm^-3]"
                      << std::endl;
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
        // section "sim"
        PS::Comm::broadcast(&sim::eps_grav, 1, root);
        PS::Comm::broadcast(&sim::mass_avg, 1, root);
        PS::Comm::broadcast(&sim::dt_max, 1, root);
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

    void calcConservedQuantities(const PS::ParticleSystem<FP_nbody>& psys_nbody,
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
        PS::F64vec L_loc {0.0};
        for (PS::S32 i = 0; i < psys_nbody.getNumberOfParticleLocal(); i++) {
            const PS::F64 mass = psys_nbody[i].mass;
            const PS::F64vec pos = psys_nbody[i].pos;
            const PS::F64vec vel = psys_nbody[i].vel;
            const PS::F64vec mom = mass * vel;
            M_dm_loc += mass;
            E_kin_loc += 0.5 * mass * SQ(vel);
            E_pot_loc += 0.5 * mass * (psys_nbody[i].pot + mass / run_param::sim::eps_grav);
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
            E_pot_loc += 0.5 * mass * (psys_star[i].pot + mass / run_param::sim::eps_grav);
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
            E_pot_loc += 0.5 * mass * (psys_gas[i].pot_grav + mass / run_param::sim::eps_grav);
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
            const PS::F64 M_err = std::fabs((cq::M_tot - cq::M_tot_ini)/cq::M_tot_ini);
            const PS::F64 E_err = std::fabs((cq::E_tot - cq::E_tot_ini)/cq::E_tot_ini);
            PS::F64vec P_err;
            P_err.x = std::fabs((cq::P.x - cq::P_ini.x)/cq::P_ini.x);
            P_err.y = std::fabs((cq::P.y - cq::P_ini.y)/cq::P_ini.y);
            P_err.z = std::fabs((cq::P.z - cq::P_ini.z)/cq::P_ini.z);
            PS::F64vec L_err;
            L_err.x = std::fabs((cq::L.x - cq::L_ini.x)/cq::L_ini.x);
            L_err.y = std::fabs((cq::L.y - cq::L_ini.y)/cq::L_ini.y);
            L_err.z = std::fabs((cq::L.z - cq::L_ini.z)/cq::L_ini.z);

            std::cout << "------------------------------------------------------------------" << std::endl;
            std::cout << "Mass (dark matter)   = " << cq::M_dm << std::endl;
            std::cout << "Mass (star)          = " << cq::M_star << std::endl;
            std::cout << "Mass (gas)           = " << cq::M_gas << std::endl;
            std::cout << "Mass (total)         = " << cq::M_tot << " (" << M_err << ")" << std::endl;
            std::cout << "Energy (kinetic)     = " << cq::E_kin << std::endl;
            std::cout << "Energy (potential)   = " << cq::E_pot << std::endl;
            std::cout << "Energy (thermal)     = " << cq::E_th  << std::endl;
            std::cout << "Energy (total)       = " << cq::E_tot << " (" << E_err << ")" << std::endl;
            std::cout << "Momentum (x)         = " << cq::P.x << " (" << P_err.x << ")" <<std::endl;
            std::cout << "Momentum (y)         = " << cq::P.y << " (" << P_err.y << ")" << std::endl;
            std::cout << "Momentum (z)         = " << cq::P.z << " (" << P_err.z << ")" << std::endl;
            std::cout << "Angular momentum (x) = " << cq::L.x << " (" << L_err.x << ")" << std::endl;
            std::cout << "Angular momentum (y) = " << cq::L.y << " (" << L_err.y << ")" << std::endl;
            std::cout << "Angular momentum (z) = " << cq::L.z << " (" << L_err.z << ")" << std::endl;
            std::cout << "------------------------------------------------------------------" << std::endl;
        }
    }
}
