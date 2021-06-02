/* C++ header */
#include "common.h"
/* FDPS header */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "debug_utilities.h"
#include "run_parameters.h"
#include "user_defined.h"
#include "hydrodynamics.h"

void calcDensity(PS::ParticleSystem<FP_gas> & psys_gas,
                 PS::ParticleSystem<FP_star> & psys_star,
                 PS::DomainInfo & dinfo, 
                 PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather & tree,
                 const PS::F64 t, const PS::F64 dt) {
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::S32 my_rank = PS::Comm::getRank();
    // Create an instance of ParticleSystem class for FB particles
    PS::ParticleSystem<FP_star> psys_fb;
    psys_fb.initialize();
    psys_fb.setNumberOfParticleLocal(0);
    std::vector<PS::S32> idx;
    for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
        if (psys_star[i].feedbackAsSNII(t,dt) ||
            psys_star[i].feedbackAsSNIa(t,dt) ||
            psys_star[i].feedbackAsAGB(t,dt)  ||
            psys_star[i].feedbackAsNSM(t,dt)) {
            psys_fb.addOneParticle(psys_star[i]);
            idx.push_back(i);
        }
    }
    PS::S32 n_fb = psys_fb.getNumberOfParticleGlobal();
    // Determine the density and the smoothing length so that Eq.(6) in Springel (2005) 
    // holds within a specified accuracy.
    const PS::S64 n_glb = psys_gas.getNumberOfParticleGlobal()
                        + psys_fb.getNumberOfParticleGlobal();
    run_param::sim::SCF_smth = 1.25;
    PS::S32 iter = 0;
    for (;;) {
        iter++;
        // Compute density, etc.
        tree.setParticleLocalTree(psys_gas);
        tree.setParticleLocalTree(psys_fb,false);
        tree.calcForceMakingTree(CalcDensity(), dinfo);
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            psys_gas[i].copyFromForce(tree.getForce(i));
        }
        PS::S32 offset = psys_gas.getNumberOfParticleLocal();
        for (PS::S32 i = 0; i < psys_fb.getNumberOfParticleLocal(); i++) {
            psys_fb[i].copyFromForce(tree.getForce(i+offset));
        }

        // Nan check
#if defined(ENABLE_NAN_CHECK)
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            if (std::isnan(psys_gas[i].smth) ||
                std::isinf(psys_gas[i].smth) ||
                (psys_gas[i].smth <= 0.0) ||
                std::isnan(psys_gas[i].dens) ||
                std::isinf(psys_gas[i].dens) ||
                (psys_gas[i].dens <= 0.0)) {
                psys_gas[i].dump("[nan]", i, __func__, dbg_utils::fout);
#if defined(FORCE_QUIT_IN_NAN_CHECK)
                std::cout << "iter = " << iter << std::endl;
                assert(false);
#endif
            }
        }
        for (PS::S32 i = 0; i < psys_fb.getNumberOfParticleLocal(); i++) {
            if (std::isnan(psys_fb[i].FBrad) ||
                std::isinf(psys_fb[i].FBrad) ||
                (psys_fb[i].FBrad <= 0.0)    ||
                std::isnan(psys_fb[i].dens)  ||
                std::isinf(psys_fb[i].dens)  ||
                (psys_fb[i].dens <= 0.0)) {
                psys_fb[i].dump("[nan]", i, __func__, dbg_utils::fout);
#if defined(FORCE_QUIT_IN_NAN_CHECK)
                std::cout << "iter = " << iter << std::endl;
                assert(false);
#endif
            }
        }
#endif

        // Check convergence
        PS::S32 n_compl_loc = 0;
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            if (psys_gas[i].flag == 1) n_compl_loc++;
        }
        for (PS::S32 i = 0; i < psys_fb.getNumberOfParticleLocal(); i++) {
            if (psys_fb[i].flag == 1) n_compl_loc++;
        }
        const PS::S64 n_compl = PS::Comm::getSum(n_compl_loc);
        if (n_compl == n_glb) break;
    }
    // Reset SCF_smth
    run_param::sim::SCF_smth = 1.0;
    // Update the feedback radius
    for (PS::S32 i = 0; i < idx.size(); i++) {
        psys_star[idx[i]].dens  = psys_fb[i].dens;
        psys_star[idx[i]].FBrad = psys_fb[i].FBrad;
    }
}

void setEntropy(PS::ParticleSystem<FP_gas>& psys){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        psys[i].setEntropy();
    }
}

void setPressure(PS::ParticleSystem<FP_gas>& psys){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        psys[i].setPressure();
    }
}
