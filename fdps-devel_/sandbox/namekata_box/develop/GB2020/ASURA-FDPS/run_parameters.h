#pragma once
/* C++ header */
#include "common.h"
/* FDPS header */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "user_defined.h"
/* CELib header */
#include "CELib.h"

namespace run_parameters {
    enum class RunStatus {InitialRun, RestartRun};

    extern RunStatus run_stat;

    namespace pseudo_random_number_generator {
        extern std::mt19937_64 mt;
    }
    namespace prng = pseudo_random_number_generator;

    namespace unit_system {
        extern PS::F64 mass;
        extern PS::F64 leng;
        extern PS::F64 time;
        extern PS::F64 velc;
        extern PS::F64 eng;
        extern PS::F64 dens;
        extern PS::F64 temp;
        extern PS::F64 spen;
    }
    namespace unit = unit_system;

    namespace basic {
        extern PS::BOUNDARY_CONDITION bc;
        extern PS::F64ort pos_root_domain;
        extern PS::S64 nstep;
        extern PS::F64 time;
        extern PS::F64 time_end;
    }

    namespace IO {
        extern PS::S32 ndump;
        extern PS::F64 time_dump;
        extern PS::F64 dt_dump;

        extern PS::S32 ndump_rst;
        extern PS::F64 time_dump_rst;
        extern PS::F64 dt_dump_rst;
    }
    namespace io = IO;

    namespace conserved_quantities {
        extern PS::F64 M_dm_ini; 
        extern PS::F64 M_star_ini;
        extern PS::F64 M_gas_ini;
        extern PS::F64 M_elms_ini[CELibYield_Number];
        extern PS::F64 M_tot_ini;
        extern PS::F64 E_kin_ini;
        extern PS::F64 E_pot_ini;
        extern PS::F64 E_th_ini;
        extern PS::F64 E_tot_ini;
        extern PS::F64vec P_ini;
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
        extern PS::F64 L_ini;
#else
        extern PS::F64vec L_ini;
#endif

        extern PS::F64 M_dm;
        extern PS::F64 M_star;
        extern PS::F64 M_gas;
        extern PS::F64 M_elms[CELibYield_Number];
        extern PS::F64 M_tot;
        extern PS::F64 E_kin;
        extern PS::F64 E_pot;
        extern PS::F64 E_th;
        extern PS::F64 E_tot;
        extern PS::F64vec P;
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
        extern PS::F64 L;
#else
        extern PS::F64vec L;
#endif

        extern PS::F64 E_lost_via_ch; 
        extern PS::F64 E_gain_via_fb;
    }
    namespace cq = conserved_quantities;

    namespace gravity {
        namespace soft {
            extern PS::F64 theta;
            extern PS::S32 n_group_limit;
            extern PS::S32 n_leaf_limit;
            extern PS::F64 eps_dm;
            extern PS::F64 eps_star;
            extern PS::F64 eps_gas;
            extern PS::F64 dt_fid;
            extern PS::F64 dt_actl;
            extern PS::F64 CFL;
            /* PMMM parameters will be added here */
        }
        namespace hard {
            extern PS::F64 eps;
            extern PS::F64 eta;
            extern PS::F64 eta_s;
            extern PS::F64 dt_limit;
        }
    }
    namespace grav = gravity;

    namespace SPH {
        extern PS::S32 n_group_limit;
        extern PS::S32 n_leaf_limit;
        extern PS::S32 n_jp_limit;
        extern PS::F64 h2h_next;
        extern PS::F64 CFL;
        extern PS::F64 eps;
        extern PS::S32 N_ngb;
        extern PS::S32 N_ngb_mgn;
        extern PS::S32 N_ngb_limit;
        extern PS::F64 alpha_AV_ini;
        extern PS::F64 alpha_AV_min;
        extern PS::F64 alpha_AV_max;
        extern PS::F64 alpha2beta_AV;
        extern PS::F64 ell_AV;
        extern PS::F64 T_lower_limit;
        extern PS::F64 T_upper_limit;
    }
    namespace sph = SPH;

    namespace ISM {
        // Solar mass abundance
        extern const PS::F64 Xhydrogen_solar;
        extern const PS::F64 Yhelium_solar;
        extern const PS::F64 Zmetal_solar;
        extern PS::F64 Zcarbon_solar;
        extern PS::F64 Znitrogen_solar;
        extern PS::F64 Zoxygen_solar;
        extern PS::F64 Zneon_solar;
        extern PS::F64 Zmagnesium_solar;
        extern PS::F64 Zsilicon_solar;
        extern PS::F64 Zsulfur_solar;
        extern PS::F64 Zcalcium_solar;
        extern PS::F64 Ziron_solar;
        extern PS::F64 Znickel_solar;
        extern PS::F64 Zeuropium_solar;

        // Mass abundance
        extern PS::F64 Xhydrogen;
        extern PS::F64 Yhelium;
        extern PS::F64 Zmetal;

        // Gas properties
        extern PS::F64 gamma;
        extern PS::F64 mu;

        // FUV background
        extern PS::F64 epsilon_FUV;
        extern PS::F64 G0_FUV;
    }
    namespace ism = ISM;

    namespace star_formation {
        extern PS::F64 nH_threshold;
        extern PS::F64 T_threshold;
        extern PS::F64 C_star;
        extern PS::S32 f_spawn;
    }
    namespace sf = star_formation;

    namespace stellar_feedback {
        extern PS::S32 N_ngb;
        extern PS::S32 N_ngb_mgn;
    }
    namespace fb = stellar_feedback;


    void setRunStat(const RunStatus val);
    bool isInitialRun();
    bool isRestartRun();
    void init(const std::size_t seed);
    void writeRankDependentFile(const std::string dir_name,
                                const std::string file_id);
    void readRankDependentFile(const std::string dir_name,
                               const std::string file_id);
    void writeRankSharedFile(const std::string dir_name,
                             const std::string file_id);
    void readRankSharedFile(const std::string dir_name,
                            const std::string file_id);
    void bcast(const int root = 0);
    void setup(const std::string dir_name = "",
               const std::string file_id = "",
               const bool read_rank_dep_file = true);
    void writeFile(const std::string dir_name,
                   const std::string file_id);
    void calcConservedQuantities(const PS::ParticleSystem<FP_dm>& psys_dm,
                                 const PS::ParticleSystem<FP_star>& psys_star,
                                 const PS::ParticleSystem<FP_gas>& psys_gas);

}
namespace run_param = run_parameters;
