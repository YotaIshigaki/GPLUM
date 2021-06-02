/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "debug_utilities.h"
#include "mathematical_constants.h"
#include "physical_constants.h"
#include "run_parameters.h"
#include "user_defined.h"
#include "star_formation.h"
/* CELib header */
#include "CELib.h"

bool StarFormation(PS::ParticleSystem<FP_gas> & psys_gas,
                   PS::ParticleSystem<FP_star> & psys_star,
                   const PS::F64 t, const PS::F64 dt) {
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::S32 my_rank = PS::Comm::getRank();

    // Make an instance of a distribution
    std::uniform_real_distribution<double> dist(0.0, 1.0);

#ifdef ENABLE_SF_CHECK
    // Calculate the total masses of gas and stars before star formation
    PS::F64 M_gas_bef {0.0};
    for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++)
        M_gas_bef += psys_gas[i].mass;
    M_gas_bef = PS::Comm::getSum(M_gas_bef);
    PS::F64 M_star_bef {0.0};
    for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++)
        M_star_bef += psys_star[i].mass;
    M_star_bef = PS::Comm::getSum(M_star_bef);
#endif

    // Process star formation
    const PS::S32 n_loc = psys_gas.getNumberOfParticleLocal();
    PS::S32 flag_loc {0};
    std::vector<PS::S32> idx;
    for (PS::S32 i = 0; i < n_loc; i++) {

        const PS::F64 nH = psys_gas[i].dens * run_param::unit::dens * run_param::ism::Xhydrogen / phys_const::Mhydrogen;
        const PS::F64 T = psys_gas[i].getTemperature() * run_param::unit::temp;
        // Check if star formation criterions are fulfilled or not
        if ((nH > run_param::sf::nH_threshold) &&
            (T  < run_param::sf::T_threshold) &&
            (psys_gas[i].divv < 0.0)) {
            // Calculate the probability of star formation
            const PS::F64 m = psys_gas[i].mass;
            PS::F64 mspawn = psys_gas[i].mass0 / static_cast<double>(run_param::sf::f_spawn);
            const PS::F64 rho = psys_gas[i].dens * run_param::unit::dens;
            const PS::F64 tdyn = 1.0/std::sqrt(4.0 * math_const::pi * phys_const::Ggrav * rho) / run_param::unit::time;
            const PS::F64 p = (m/mspawn)*(1.0 - exp(-run_param::sf::C_star * dt / tdyn));
            // Calculate properties of a newly-born star
            const PS::F64 prnd = dist(run_param::prng::mt);
            if (prnd < p) {
                flag_loc = 1;
                // Update gas particle and check if this gas particle
                // should be removed or not
                psys_gas[i].n_stars++;
                if (psys_gas[i].n_stars == run_param::sf::f_spawn) {
                    idx.push_back(i);
                    mspawn = psys_gas[i].mass; // re-calculate to minimize error
                }
                psys_gas[i].mass -= mspawn;
                // Note that other quantities depending on mass, e.g. dens,
                // are re-calculated in the next density calculation.

                // Set a newly-born star particle
                FP_star ptcl;
                ptcl.pid = psys_gas[i].id;
                ptcl.id  = psys_gas[i].n_stars - 1;
                ptcl.mass0 = mspawn;
                ptcl.mass = mspawn;
                ptcl.pos = psys_gas[i].pos;
                ptcl.vel = psys_gas[i].vel;
                ptcl.acc = psys_gas[i].acc_grav;
                ptcl.pot = psys_gas[i].pot_grav;
                // feedback information
                ptcl.copyAbundanceFrom(psys_gas[i].mabn);
                ptcl.t_form = t;
                ptcl.t_SNII = t
                    + CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = dist(run_param::prng::mt),
                        .InitialMass_in_Msun = mspawn * run_param::unit::mass/phys_const::Msolar,
                        .Metallicity = ptcl.getMetallicity(),
                        .Count = 0,
                        },CELibFeedbackType_SNII)
                        * phys_const::yr / run_param::unit::time;
                ptcl.t_SNIa = t
                    + CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = dist(run_param::prng::mt),
                        .InitialMass_in_Msun = mspawn * run_param::unit::mass/phys_const::Msolar,
                        .Metallicity = ptcl.getMetallicity(),
                        .Count = 0,
                        },CELibFeedbackType_SNIa)
                        * phys_const::yr / run_param::unit::time;
                ptcl.t_AGB = t
                    + CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = dist(run_param::prng::mt),
                        .InitialMass_in_Msun = mspawn * run_param::unit::mass/phys_const::Msolar,
                        .Metallicity = ptcl.getMetallicity(),
                        .Count = 0,
                        },CELibFeedbackType_AGB)
                        * phys_const::yr / run_param::unit::time;
                ptcl.t_NSM = t
                    + CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = dist(run_param::prng::mt),
                        .InitialMass_in_Msun = mspawn * run_param::unit::mass/phys_const::Msolar,
                        .Metallicity = ptcl.getMetallicity(),
                        .Count = 0,
                        },CELibFeedbackType_NSM)
                        * phys_const::yr / run_param::unit::time;
                ptcl.FBrad = psys_gas[i].smth;
                psys_star.addOneParticle(ptcl);
            }
        }
    }

    // Delete gas particles whose mass is zero
    if (idx.size() > 0) {
        psys_gas.removeParticle(idx.data(), idx.size());
    }

    // Calculate the return value
    const PS::S32 flag_glb = PS::Comm::getSum(flag_loc);

#ifdef ENABLE_SF_CHECK
    // Calculate the total masses of gas and stars after star formation
    if (flag_glb > 0) {
        PS::F64 M_gas_aft {0.0};
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++)
            M_gas_aft += psys_gas[i].mass;
        M_gas_aft = PS::Comm::getSum(M_gas_aft);
        PS::F64 M_star_aft {0.0};
        for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++)
            M_star_aft += psys_star[i].mass;
        M_star_aft = PS::Comm::getSum(M_star_aft);
        // Output the checked result
        if (my_rank == 0) {
            std::cout << "---- the result of some checks in StarFormation ----" << std::endl;
            std::cout << "M_gas (bef.)  = " << M_gas_bef << std::endl;
            std::cout << "M_gas (aft.)  = " << M_gas_aft << std::endl;
            std::cout << "dM_gas        = " << M_gas_aft - M_gas_bef << std::endl;
            std::cout << "M_star (bef.) = " << M_star_bef << std::endl;
            std::cout << "M_star (aft.) = " << M_star_aft << std::endl;
            std::cout << "dM_star       = " << M_star_aft - M_star_bef << std::endl;
            std::cout << "M_tot (bef.)  = " << M_gas_bef + M_star_bef << std::endl;
            std::cout << "M_tot (aft.)  = " << M_gas_aft + M_star_aft << std::endl;
        }
    }
#endif

    // Return the return value
    if (flag_glb > 0) return true;
    else return false;

} 
