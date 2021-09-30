/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "physical_constants.h"
#include "run_parameters.h"
#include "user_defined.h"
#include "SPH_kernel.h"
#include "stellar_feedback.h"
/* CELib header */
#include "CELib.h"

// Local class
class FBParticle {
public:
    PS::F64vec pos;
    PS::F64 dens;
    PS::F64 FBrad;
    PS::F64 Energy;
    PS::F64 EjectaMass;
    PS::F64 Elements[CELibYield_Number];

    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64 getRSearch() const {
        return this->FBrad;
    }
};

bool StellarFeedback(PS::ParticleSystem<FP_star> & psys_star,
                     PS::ParticleSystem<FP_gas> & psys_gas,
                     PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather & tree,
                     const PS::F64 time, const PS::F64 dt) {
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::S32 my_rank = PS::Comm::getRank();

    // Make an instance of a distribution
    std::uniform_real_distribution<double> dist(0.0, 1.0);

#ifdef ENABLE_FB_CHECK
    // Calculate the total masses of gas and stars before feedback
    PS::F64 M_gas_bef {0.0};
    for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++)
        M_gas_bef += psys_gas[i].mass;
    M_gas_bef = PS::Comm::getSum(M_gas_bef);
    PS::F64 M_star_bef {0.0};
    for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++)
        M_star_bef += psys_star[i].mass;
    M_star_bef = PS::Comm::getSum(M_star_bef);
#endif

    // Make a list of local feedback particles
    std::vector<FBParticle> ptcl_loc;
    for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
        // Compute Elements of this particles
        PS::F64 Elements[CELibYield_Number];
        for (PS::S32 k = 0; k < CELibYield_Number; k++)
            Elements[k] = psys_star[i].mass * psys_star[i].mabn[k];
   
        // Clear a buffer storing feedback
        struct CELibStructFeedbackOutput fb;
        fb.Energy = 0.0;
        fb.EjectaMass = 0.0;
        fb.RemnantMass = 0.0;
        for (PS::S32 k = 0; k < CELibYield_Number; k++) fb.Elements[k] = 0.0;
        // SNII feedback
        while (psys_star[i].feedbackAsSNII(time,dt)) {
            PS::U64 cnt = psys_star[i].getSNIICount();
            struct CELibStructFeedbackOutput tmp;
            tmp = CELibGetFeedback((struct CELibStructFeedbackInput){
                        .Mass = psys_star[i].mass0,
                        .MassConversionFactor = run_param::unit::mass/phys_const::Msolar,
                        .Metallicity = psys_star[i].getMetallicity(),
                        .Elements = Elements,
                        .Count = cnt,
                        },CELibFeedbackType_SNII);
            fb.Energy     += tmp.Energy;
            fb.EjectaMass += tmp.EjectaMass;
            for (PS::S32 k = 0; k < CELibYield_Number; k++)
                fb.Elements[k] += tmp.Elements[k];
            psys_star[i].setSNIICount(++cnt);
            psys_star[i].t_SNII = psys_star[i].t_form
                    + CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = dist(run_param::prng::mt),
                        .InitialMass_in_Msun = psys_star[i].mass0 * run_param::unit::mass/phys_const::Msolar,
                        .Metallicity = psys_star[i].getMetallicity(),
                        .Count = cnt,
                        },CELibFeedbackType_SNII)
                        * phys_const::yr / run_param::unit::time;
        }
        // SNIa feedback
        while (psys_star[i].feedbackAsSNIa(time,dt)) {
            PS::U64 cnt = psys_star[i].getSNIaCount();
            struct CELibStructFeedbackOutput tmp;
            tmp = CELibGetFeedback((struct CELibStructFeedbackInput){
                        .Mass = psys_star[i].mass0,
                        .MassConversionFactor = run_param::unit::mass/phys_const::Msolar,
                        .Metallicity = psys_star[i].getMetallicity(),
                        .Elements = Elements,
                        .Count = cnt,
                        },CELibFeedbackType_SNIa);
            fb.Energy     += tmp.Energy;
            fb.EjectaMass += tmp.EjectaMass;
            for (PS::S32 k = 0; k < CELibYield_Number; k++)
                fb.Elements[k] += tmp.Elements[k];
            psys_star[i].setSNIaCount(++cnt);
            psys_star[i].t_SNIa = psys_star[i].t_form
                    + CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = dist(run_param::prng::mt),
                        .InitialMass_in_Msun = psys_star[i].mass0 * run_param::unit::mass/phys_const::Msolar,
                        .Metallicity = psys_star[i].getMetallicity(),
                        .Count = cnt,
                        },CELibFeedbackType_SNIa)
                        * phys_const::yr / run_param::unit::time;
        }
        // AGB feedback
        while (psys_star[i].feedbackAsAGB(time,dt)) {
            PS::U64 cnt = psys_star[i].getAGBCount();
            struct CELibStructFeedbackOutput tmp;
            tmp = CELibGetFeedback((struct CELibStructFeedbackInput){
                        .Mass = psys_star[i].mass0,
                        .MassConversionFactor = run_param::unit::mass/phys_const::Msolar,
                        .Metallicity = psys_star[i].getMetallicity(),
                        .Elements = Elements,
                        .Count = cnt,
                        },CELibFeedbackType_AGB);
            fb.Energy     += tmp.Energy;
            fb.EjectaMass += tmp.EjectaMass;
            for (PS::S32 k = 0; k < CELibYield_Number; k++)
                fb.Elements[k] += tmp.Elements[k];
            psys_star[i].setAGBCount(++cnt);
            psys_star[i].t_AGB = psys_star[i].t_form
                    + CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = dist(run_param::prng::mt),
                        .InitialMass_in_Msun = psys_star[i].mass0 * run_param::unit::mass/phys_const::Msolar,
                        .Metallicity = psys_star[i].getMetallicity(),
                        .Count = cnt,
                        },CELibFeedbackType_AGB)
                        * phys_const::yr / run_param::unit::time;
        }
        // NSM feedback
        while (psys_star[i].feedbackAsNSM(time,dt)) {
            PS::U64 cnt = psys_star[i].getNSMCount();
            struct CELibStructFeedbackOutput tmp;
            tmp = CELibGetFeedback((struct CELibStructFeedbackInput){
                        .Mass = psys_star[i].mass0,
                        .MassConversionFactor = run_param::unit::mass/phys_const::Msolar,
                        .Metallicity = psys_star[i].getMetallicity(),
                        .Elements = Elements,
                        .Count = cnt,
                        },CELibFeedbackType_NSM);
            fb.Energy     += tmp.Energy;
            fb.EjectaMass += tmp.EjectaMass;
            for (PS::S32 k = 0; k < CELibYield_Number; k++)
                fb.Elements[k] += tmp.Elements[k];
            psys_star[i].setNSMCount(++cnt);
            psys_star[i].t_NSM = psys_star[i].t_form
                    + CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = dist(run_param::prng::mt),
                        .InitialMass_in_Msun = psys_star[i].mass0 * run_param::unit::mass/phys_const::Msolar,
                        .Metallicity = psys_star[i].getMetallicity(),
                        .Count = cnt,
                        },CELibFeedbackType_NSM)
                        * phys_const::yr / run_param::unit::time;
        }
        // Add to a list if there is a feedback
        if (fb.Energy > 0.0 || fb.EjectaMass > 0.0) {
            // Convert the unit into the simulation unit
            fb.Energy /= run_param::unit::eng;
            fb.EjectaMass *= (phys_const::Msolar/run_param::unit::mass);
            for (PS::S32 k = 0; k < CELibYield_Number; k++)
                fb.Elements[k] *= (phys_const::Msolar/run_param::unit::mass);
            // Update the mass of star particle of interest
            psys_star[i].mass -= fb.EjectaMass;
            // Add to the list
            FBParticle tmp;
            tmp.pos        = psys_star[i].pos;
            tmp.dens       = psys_star[i].dens;
            tmp.FBrad      = psys_star[i].FBrad;
            tmp.Energy     = fb.Energy;
            tmp.EjectaMass = fb.EjectaMass;
            for (PS::S32 k = 0; k < CELibYield_Number; k++)
                tmp.Elements[k] = fb.Elements[k];
            ptcl_loc.push_back(tmp);
        }
    }

    // Broadcast FB particles 
    const PS::S32 sendcount = ptcl_loc.size();
    std::vector<PS::S32> recvcounts, recvdispls;
    recvcounts.resize(n_proc);
    recvdispls.resize(n_proc);
    PS::Comm::allGather(&sendcount, 1, &recvcounts[0]);
    PS::S32 recvcount_tot {0};
    for (PS::S32 i = 0; i < n_proc; i++) recvcount_tot += recvcounts[i];
    if (recvcount_tot == 0) {
        // In this case, no feedback occurs.
        return false;
    }
    recvdispls[0] = 0;
    for (PS::S32 i = 1; i < n_proc; i++) recvdispls[i] = recvdispls[i-1] + recvcounts[i-1];
    std::vector<FBParticle> ptcl_glb;
    ptcl_glb.resize(recvcount_tot);
    PS::Comm::allGatherV(&ptcl_loc[0], sendcount, &ptcl_glb[0], &recvcounts[0], &recvdispls[0]);

    // Notification (for debug)
    if (my_rank == 0) {
        std::cout << "There are " << recvcount_tot << " feedbacks!" << std::endl;
        PS::F64 E_input {0.0};
        for (PS::S32 i = 0; i < ptcl_glb.size(); i++) E_input += ptcl_glb[i].Energy;
        const PS::F64 E_input_cgs = E_input * run_param::unit::eng;
        std::cout << "Input energy is " << E_input_cgs << " [erg]!" << std::endl;
    }

    // Apply stellar feedback
    for (PS::S32 i = 0; i < ptcl_glb.size(); i++) {
        // Extract i-particle info.
        const PS::F64vec xi = ptcl_glb[i].pos;
        const PS::F64 dens = ptcl_glb[i].dens;
        const PS::F64 hi = ptcl_glb[i].FBrad;
        const PS::F64 Energy = ptcl_glb[i].Energy;
        const PS::F64 EjectaMass = ptcl_glb[i].EjectaMass;
        PS::F64 Elements[CELibYield_Number];
        for (PS::S32 k = 0; k < CELibYield_Number; k++)
            Elements[k] = ptcl_glb[i].Elements[k];
        // Distribute energy and mass to neighbor gas particles
        EP_hydro *epj;
        const PS::S32 n_epj = tree.getNeighborListOneParticle(ptcl_glb[i], epj);
        if (n_epj > 0) {
            // Distribute energy and mass to local gas particles
            for (PS::S32 j = 0; j < n_epj; j++) {
                if (epj[j].type == 0 && epj[j].rank == my_rank) {
                    const PS::F64 mj = epj[j].mass;
                    const PS::F64vec dr = xi - epj[j].pos;
                    const PS::F64 rij = std::sqrt(dr * dr);
                    const PS::F64 w = mj * W(rij,hi) / dens;
                    // Save the current values
                    const PS::F64 jj = epj[j].idx;
                    const PS::F64 M_old = psys_gas[jj].mass;
                    const PS::F64 U_old = psys_gas[jj].mass * psys_gas[jj].eng;
                    PS::F64 Elms_old[CELibYield_Number];
                    for (PS::S32 k = 0; k < CELibYield_Number; k++)
                        Elms_old[k] = psys_gas[jj].mass * psys_gas[jj].mabn[k];
                    // Calculate the updated values
                    const PS::F64 M_new = M_old + w * EjectaMass;
                    const PS::F64 U_new = U_old + w * Energy;
                    PS::F64 Elms_new[CELibYield_Number];
                    for (PS::S32 k = 0; k < CELibYield_Number; k++)
                        Elms_new[k] = Elms_old[k] + w * Elements[k];
                    // Update gas particle
                    psys_gas[jj].mass = M_new;
                    psys_gas[jj].eng  = U_new/M_new;
                    for (PS::S32 k = 0; k < CELibYield_Number; k++) 
                       psys_gas[jj].mabn[k] = Elms_new[k] / M_new;
                }
            }
        }
    }

#ifdef ENABLE_FB_CHECK
    // Calculate the total masses of gas and stars after feedback
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
        std::cout << "---- the result of some checks in StellarFeedback ----" << std::endl;
        std::cout << "M_gas (bef.)  = " << M_gas_bef << std::endl;
        std::cout << "M_gas (aft.)  = " << M_gas_aft << std::endl;
        std::cout << "dM_gas        = " << M_gas_aft - M_gas_bef << std::endl;
        std::cout << "M_star (bef.) = " << M_star_bef << std::endl;
        std::cout << "M_star (aft.) = " << M_star_aft << std::endl;
        std::cout << "dM_star       = " << M_star_aft - M_star_bef << std::endl;
        std::cout << "M_tot (bef.)  = " << M_gas_bef + M_star_bef << std::endl;
        std::cout << "M_tot (aft.)  = " << M_gas_aft + M_star_aft << std::endl;
    }
#endif

    return true;
}
