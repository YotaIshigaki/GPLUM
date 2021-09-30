/* C++ header */
#include "common.h"
/* FDPS header */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "debug_utilities.h"
#include "timing_utilities.h"
#include "run_parameters.h"
#include "user_defined.h"
#include "hydrodynamics.h"

void predictFeedbackRadius(PS::ParticleSystem<FP_star> & psys_star,
                           PS::DomainInfo & dinfo,
                           PS::TreeForForceShort<Force_hydro, EP_hydro, EP_hydro>::Symmetry & tree,
                           const PS::F64 t, const PS::F64 dt) {
    // Predict feedback radius using the most nearest local gas particle
    for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
        if (psys_star[i].feedbackAsSNII(t,dt) ||
            psys_star[i].feedbackAsSNIa(t,dt) ||
            psys_star[i].feedbackAsAGB(t,dt)  ||
            psys_star[i].feedbackAsNSM(t,dt)) {
            EP_hydro * epj;
            const PS::S32 N_ngb_found = tree.getNeighborListOneParticle(psys_star[i], 1, epj);
            if (N_ngb_found > 0) psys_star[i].FBrad = epj->h;
#if 0
            // Check [for debug] 
            if (PS::Comm::getRank() == 492) {
                const PS::F64vec pos_i = psys_star[i].pos;
                const PS::F64 h_i = psys_star[i].FBrad;
                std::cout << "pos_i = " << pos_i << std::endl;
                std::cout << "h_i   = " << h_i << std::endl;
                for (PS::S32 j = 0; j < N_ngb_found; j++) {
                    const PS::F64vec pos_j = epj[j].pos;
                    const PS::F64 h_j = epj[j].h;
                    const PS::F64vec dr = pos_i - pos_j;
                    const PS::F64 rij = std::sqrt(dr * dr);
                    std::cout << "j = " << j
                              << " pos_j = " << pos_j
                              << " h_j = " << h_j
                              << " rij/h_j = " << rij/h_j
                              << std::endl;
                }
            }
#endif
        }
    }
}

void calcKernelSize(PS::ParticleSystem<FP_gas> & psys_gas,
                    PS::ParticleSystem<FP_star> & psys_star,
                    PS::DomainInfo & dinfo, 
                    PS::TreeForForceShort<Force_knl_sz, EPI_knl_sz, EPJ_knl_sz>::Gather & tree,
                    const PS::F64 t, const PS::F64 dt,
                    const FunctionUseType use_type) {
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::S32 my_rank = PS::Comm::getRank();
    et_rij_calc = 0.0;
    et_bisec = 0.0;
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
    // Predict the smoothing length 
    if (use_type == FunctionUseType::Main) {
        constexpr PS::F64 coeff_dec_lim = 0.9;
        constexpr PS::F64 coeff_inc_lim = 1.25;
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
#if 0
            const PS::F64 h_prev = psys_gas[i].h;
            const PS::F64 h_pred_tmp = h_prev + psys_gas[i].h_dot_prev * dt;
            psys_gas[i].h = std::min(coeff_inc_lim * h_prev,
                                     std::max(coeff_dec_lim * h_prev,
                                     h_pred_tmp));
            psys_gas[i].h_prev = h_prev;
#else
            psys_gas[i].h_prev = psys_gas[i].h;
            psys_gas[i].h *= run_param::sph::h2h_next;
#endif
        }
    }
    // Determine the smoothing length so that the number of neighbor
    // particles is in a specified range or the effective number of
    // neighbor particles becomes a specified value.
    const PS::S64 n_glb_gas = psys_gas.getNumberOfParticleGlobal();
    const PS::S64 n_glb_fb = psys_fb.getNumberOfParticleGlobal();
    const PS::S64 n_glb = n_glb_gas + n_glb_fb;
    for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) psys_gas[i].flag = 0;
    for (PS::S32 i = 0; i < psys_fb.getNumberOfParticleLocal(); i++) psys_fb[i].flag = 0;
    tree.clearTimeProfile();
    tree.clearNumberOfInteraction();
    PS::S32 iter = 0;
    for (;;) {
        iter++;
        // for debug
        if (PS::Comm::getRank() == 0) std::cout << "iter = " << iter << std::endl;
        // Compute kernel size, etc.
        bool clear_flag {true};
        tree.setParticleLocalTree(psys_gas, clear_flag);
        if (n_glb_gas > 0) clear_flag = false;
        tree.setParticleLocalTree(psys_fb, clear_flag);
        tree.calcForceMakingTree(CalcKernelSize(), dinfo);
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
            if (std::isnan(psys_gas[i].h) ||
                std::isinf(psys_gas[i].h) ||
                (psys_gas[i].h <= 0.0)) {
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
                (psys_fb[i].FBrad <= 0.0) ) {
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
        // TODO: measure between the end of tree and getSum().
    }
    // Calculate the time rate of change of the smoothing length
    if (use_type == FunctionUseType::Main) {
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            psys_gas[i].h_dot_prev = (psys_gas[i].h - psys_gas[i].h_prev) / dt;
        }
    }
    // Update the feedback radius
    for (PS::S32 i = 0; i < idx.size(); i++) {
        psys_star[idx[i]].FBrad = psys_fb[i].FBrad;
    }
    // Update time_prof
    if (use_type == FunctionUseType::Main)
        time_prof.calc_knl_sz_1st += tree.getTimeProfile().getTotalTime();
    else
        time_prof.calc_knl_sz_2nd += tree.getTimeProfile().getTotalTime(); 

    // Check
    if (run_param::basic::nstep == 11450) {
    //if (run_param::basic::nstep == 11329) {
    //if (run_param::basic::nstep == 11350) {
        //dbg_utils::fout << "et_rij_calc = " << et_rij_calc << std::endl;
        //dbg_utils::fout << "et_bisec    = " << et_bisec << std::endl;
        //tree.getTimeProfile().dump(dbg_utils::fout);
        //dbg_utils::fout << tree.getNumberOfEpiSorted() << "   "
        //                << tree.getNumberOfEpjSorted() << std::endl;
        //PS::Finalize(); std::exit(0);
    }
}

void calcDensityAndPressure(PS::ParticleSystem<FP_gas> & psys_gas,
                            PS::ParticleSystem<FP_star> & psys_star,
                            PS::DomainInfo & dinfo, 
                            PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather & tree,
                            const PS::F64 t, const PS::F64 dt,
                            const FunctionUseType use_type) {
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
#if 0
    // Set output flag
    if (run_param::basic::nstep == 11345) {
        output_flag = true;
        tree.debug_flag_ = true;
    } else {
        output_flag = false;
        tree.debug_flag_ = false;
    }
#endif
    // Clear internal logs
    tree.clearTimeProfile();
    tree.clearNumberOfInteraction();
    // Compute density, etc.
    bool clear_flag {true};
    tree.setParticleLocalTree(psys_gas, clear_flag);
    if (psys_gas.getNumberOfParticleGlobal() > 0) clear_flag = false;
    tree.setParticleLocalTree(psys_fb, clear_flag);
    tree.calcForceMakingTree(CalcDensityAndPressure(), dinfo);
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
        if (std::isnan(psys_gas[i].dens) ||
            std::isinf(psys_gas[i].dens) ||
            (psys_gas[i].dens <= 0.0)) {
            psys_gas[i].dump("[nan]", i, __func__, dbg_utils::fout);
#if defined(FORCE_QUIT_IN_NAN_CHECK)
            assert(false);
#endif
        }
    }
    for (PS::S32 i = 0; i < psys_fb.getNumberOfParticleLocal(); i++) {
        if (std::isnan(psys_fb[i].dens)  ||
            std::isinf(psys_fb[i].dens)  ||
            (psys_fb[i].dens <= 0.0)) {
            psys_fb[i].dump("[nan]", i, __func__, dbg_utils::fout);
#if defined(FORCE_QUIT_IN_NAN_CHECK)
            assert(false);
#endif
        }
    }
#endif

    // Update the density at the particle position
    for (PS::S32 i = 0; i < idx.size(); i++) {
        psys_star[idx[i]].dens  = psys_fb[i].dens;
    }

    // Update time_prof
    if (use_type == FunctionUseType::Main)
        time_prof.calc_density_1st += tree.getTimeProfile().getTotalTime();
    else
        time_prof.calc_density_2nd += tree.getTimeProfile().getTotalTime();

#if 0
    // [DEBUG] Output data to analyze LET size and the length of interaction list
    if (run_param::basic::nstep == 11345) {
        std::stringstream ss;
        ss << "fp_pos_h_" << std::setw(5) << std::setfill('0') << PS::Comm::getRank() << ".txt";
        const std::string file_name = ss.str();
        std::ofstream ofs;
        ofs.open(file_name.c_str(), std::ios::trunc);
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            const PS::F64vec pos = psys_gas[i].pos;
            const PS::F64 h = psys_gas[i].h;
            ofs << pos.x << "   " << pos.y << "   " << pos.z  << "   "
                << h << "   " << 1 << std::endl;
        }
        for (PS::S32 i = 0; i < psys_fb.getNumberOfParticleLocal(); i++) {
            const PS::F64vec pos = psys_fb[i].pos;
            const PS::F64 h = psys_fb[i].FBrad;
            ofs << pos.x << "   " << pos.y << "   " << pos.z  << "   "
                << h << "   " << 2 << std::endl;
        }
        ofs.close();
        PS::Finalize(); std::exit(0);
    }
#endif

    // Check
    //if (run_param::basic::nstep == 11450) {
    //    //tree.getTimeProfile().dump(dbg_utils::fout);
    //    dbg_utils::fout << tree.getNumberOfEpiSorted() << "   "
    //                    << tree.getNumberOfEpjSorted() << std::endl;
    //    PS::Finalize(); std::exit(0);
    //}
}

void calcSoundSpeed(PS::ParticleSystem<FP_gas>& psys){
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        psys[i].calcSoundSpeed();
    }
}

void calcArtificialViscosity(PS::ParticleSystem<FP_gas>& psys,
                             const PS::F64 dt) {
#if ASURA_FDPS_ARTIFICIAL_VISCOSITY == ASURA_FDPS_CONSTANT_VISCOSITY
    // do nothing
#elif ASURA_FDPS_ARTIFICIAL_VISCOSITY == ASURA_FDPS_MORRIS_MONAGHAN_1997
    const PS::F64 ell = run_param::sph::ell_AV;
    const PS::F64 alpha_min = run_param::sph::alpha_AV_min;
    const PS::F64 alpha_max = run_param::sph::alpha_AV_max;
    for (PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        PS::F64 alpha = psys[i].alpha;
        const PS::F64 tau = psys[i].h/(2.0 * ell * psys[i].snds);
        const PS::F64 S = std::max(-psys[i].divv, 0.0);
        alpha += ((alpha_min - alpha) / tau + (alpha_max - alpha) * S) * dt;
        psys[i].alpha = std::min(alpha_max, std::max(alpha_min, alpha));
    }
#elif ASURA_FDPS_ARTIFICIAL_VISCOSITY == ASURA_FDPS_CULLEN_DEHNEN_2010
#error The artificial viscosity by Cullen & Dehnen (2010) is not supported yet.
#endif
}

void checkDensityFluctuation(PS::ParticleSystem<FP_gas>& psys) {
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    const PS::S64 n_glb = psys.getNumberOfParticleGlobal();
    // Compute the average of density
    PS::F64 tmp = 0.0;
    for (PS::S32 i = 0; i < n_loc; i++) tmp += psys[i].dens;
    const PS::F64 dens_avg = PS::Comm::getSum(tmp) / n_glb;
    // Compute the dispersion of density
    tmp = 0.0;
    for (PS::S32 i = 0; i < n_loc; i++) 
        tmp += std::pow(psys[i].dens - dens_avg, 2.0);
    const PS::F64 dens_disp = std::sqrt(PS::Comm::getSum(tmp)/n_glb);
    // Output the status 
    const PS::F64 fluc_str = dens_disp / dens_avg;
    if (PS::Comm::getRank() == 0) {
        std::cout << "---------------------------" << std::endl;
        std::cout << "avg.       = " << dens_avg << std::endl;
        std::cout << "disp.      = " << dens_disp << std::endl;
        std::cout << "disp./avg. = " << fluc_str << std::endl;
    }
    // Output data and end the simulation if the fluctuation is small
    constexpr PS::F64 eps = 1.672e-3;
    if (fluc_str < eps) {
        char filename[256];
        sprintf(filename, "result/glass_data.txt");
        psys.writeParticleBinary(filename, &FP_gas::writeGlassData);
        if (PS::Comm::getRank() == 0) {
            std::cout << "A glass-like distribution is obtained." << std::endl;
            std::cout << "The particle position data is output as file " 
                      << filename << std::endl;
        }
        PS::Finalize();
        std::exit(0);
    }
}
