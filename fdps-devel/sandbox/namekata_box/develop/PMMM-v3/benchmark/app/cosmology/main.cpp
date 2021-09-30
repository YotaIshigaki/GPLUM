// Include header files of the C++ standard template library
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>
#include <float.h>
#include <cstdio>
#include <cstdlib>
// Include header files of FDPS
#include <particle_simulator.hpp>
#include <particle_mesh_multipole.hpp>
#include <ewald.hpp>
// Include header file of Phantom-GRAPE library
#ifdef ENABLE_PHANTOM_GRAPE_X86
#include <gp5util.h>
#endif
// Include user-defined headers
#include "mathematical_constants.h"
#include "cosmology.hpp"
#include "run_param.hpp"
#include "user_defined.hpp"
#include "timestep.hpp"
#include "leapfrog.hpp"
#include "io.hpp"

template <class Tpsys>
void setup_random_particle(Tpsys & psys,
                           run_param & this_run,
                           const PS::S32 npart_total)
{
    PS::S32 rank = PS::Comm::getRank();
    PS::S64 npart_local = (rank == 0) ? npart_total : 0;

    psys.setNumberOfParticleLocal(npart_local);

    this_run.npart_total = npart_total;
    this_run.mpi_nproc = PS::Comm::getNumberOfProc();
    this_run.mpi_rank = rank;

    const long seed=19810614;
    srand48(seed);

    for (PS::S32 i=0;i<npart_local;i++) {
        psys[i].id = i;
        psys[i].mass = 3.0*this_run.cosm.omegam/(8.0*math_const::pi*(PS::F32)npart_total);
        psys[i].pos = PS::F64vec(drand48(), drand48(), drand48());
        psys[i].vel = 0.0;
    }
}


template <class Tpsys>
void read_SB_particle(Tpsys &psys, 
                      run_param &this_run, 
                      const char * input_file){
    psys.readParticleBinary(input_file);
    this_run.npart_total = psys.getNumberOfParticleGlobal();
    this_run.npart_local = psys.getNumberOfParticleLocal();
    this_run.mpi_nproc = PS::Comm::getNumberOfProc();
    this_run.mpi_rank = PS::Comm::getRank();
}

template<class Tpsys>
void read_param_file(Tpsys & psys,
                     run_param & this_run,
                     const char * input_file){
    std::cout << "input_file=" << input_file << std::endl;
    FILE * param_fp;
    param_fp = fopen(input_file, "r");
    if (param_fp == NULL) {
        fprintf(stderr, "File %s not found in input_params.\n", input_file);
        PS::Abort();
        std::exit(EXIT_FAILURE);
    }
    int ret = fscanf(param_fp, "%d", &this_run.mode);
    if (ret == EOF){
        fprintf(stderr, "Input error of mode in input_params()\n");
        PS::Abort();
        std::exit(EXIT_FAILURE);
    }
    ret = fscanf(param_fp, "%lf", &FP_grav::eps);
    ret = fscanf(param_fp, "%lf", &this_run.theta);
    ret = fscanf(param_fp, "%f", &this_run.zend);
    if (PS::Comm::getRank() == 0){
        std::cout << "Gravitational softening = " << FP_grav::eps << std::endl;
        std::cout << "this_run.theta = " << this_run.theta << std::endl;
        std::cout << "this_run.zend = " << this_run.zend << std::endl;
    }
    if (this_run.mode == run_param::SANTABARBARA){
        this_run.znow = 63.0;
        this_run.cosm.omegam = 1.0;
        this_run.cosm.omegav = this_run.cosm.omegab = this_run.cosm.omeganu = 0.0;
        FP_grav::H0 = 50.0; // Hubble parameter at present time
        FP_grav::Lbnd = 64.0; // the physical length of the side of a (cubic) computational box
        char snap_name[1024];
        fscanf(param_fp, "%s", snap_name);
        if (PS::Comm::getRank() == 0) std::cout << "snap_name:" << snap_name << std::endl;
        read_SB_particle(psys, this_run, snap_name);
    }
    else if (this_run.mode == run_param::READ_FILE){
        char snap_name[1024];
        fscanf(param_fp, "%s", snap_name);
        if (PS::Comm::getRank() == 0) std::cout << "snap_name:" << snap_name << std::endl;
        ret = fscanf(param_fp, "%f", &this_run.znow);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omegam);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omegav);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omegab);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omeganu);
        ret = fscanf(param_fp, "%lf", &FP_grav::H0);
        ret = fscanf(param_fp, "%lf", &FP_grav::Lbnd);
        read_SB_particle(psys, this_run, snap_name);    
    }
    else if (this_run.mode == run_param::RANDOM){
        ret = fscanf(param_fp, "%f", &this_run.znow);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omegam);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omegav);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omegab);
        ret = fscanf(param_fp, "%f", &this_run.cosm.omeganu);
        ret = fscanf(param_fp, "%lld", &this_run.npart_total);
        if (PS::Comm::getRank() == 0) std::cout << "npart_total = " << this_run.npart_total << std::endl; 
        FP_grav::H0 = 50.0;
        FP_grav::Lbnd = 64.0;
        setup_random_particle(psys, this_run, this_run.npart_total);
    }
    else{
        PS::Abort();
        std::exit(EXIT_FAILURE);
    }    
    this_run.anow = 1.0 / (1.0 + this_run.znow);
    this_run.tnow = this_run.cosm.atotime(this_run.anow);
    this_run.update_expansion(this_run.tnow);
    this_run.input_params(param_fp);
}

void calc_total_force(PS::ParticleSystem<FP_grav> & psys) {
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    for (PS::S32 i=0; i<n_loc; i++) psys[i].calcTotalForce();
}

void check_error(PS::ParticleSystem<FP_grav> & psys) {
    if (PS::Comm::getRank() == 0) {
        std::cout << "check_error() started." << std::endl;
        if (FP_grav::eps != 0.0) {
            std::cout << "Note that eps is NOT 0!" << std::endl;
        }
    }
    // Make an temporal instance of ParticleSystem
    const PS::S32 n_glb = psys.getNumberOfParticleGlobal();
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    PS::ParticleSystem<FP_grav> psys_tmp;
    psys_tmp.setNumberOfParticleLocal(n_loc);
    for (PS::S32 i=0; i<n_loc; i++) {
        psys_tmp[i].id = psys[i].id;
        psys_tmp[i].mass = psys[i].mass;
        psys_tmp[i].pos = psys[i].pos;
        psys_tmp[i].vel = psys[i].vel;
    }
    // Perform interaction calculation with the Ewald method
    const PS::F64ort pos_unit_cell = PS::F64ort(PS::F64vec(0.0), PS::F64vec(1.0));
    const PS::F64 alpha = 2.4;
    //const PS::S32 NMIR = 3;
    const PS::S32 NMIR = 5;
    const PS::S32 MMAX = 5;
    //const PS::S32 MMAX = 7;
    PS::Ewald::Ewald ewald;
    ewald.initialize(pos_unit_cell,alpha,NMIR,MMAX);
    ewald.calcForceAllAndWriteBack(psys_tmp);
    {
        std::stringstream ss;
        ss << "ewald" << std::setw(5) << std::setfill('0') << PS::Comm::getRank() << ".txt";
        const std::string file_name = ss.str();
        std::ofstream output_file;
        output_file.open(file_name.c_str(), std::ios::trunc);
        output_file << std::setprecision(15) << std::scientific;
        for (PS::S32 i = 0; i < n_loc; i++) {
            output_file << psys_tmp[i].id << "   "
                        << psys_tmp[i].acc << "   "
                        << psys_tmp[i].pot << std::endl;
        }
        output_file.close();
    }
    // Compare the result of the Ewald method with that of PMMM
    std::stringstream ss;
    ss << "err_info_" << std::setw(5) << std::setfill('0') << PS::Comm::getRank() << ".txt";
    const std::string filename = ss.str();
    std::ofstream log_file;
    log_file.open(filename.c_str(), std::ios::trunc);
    PS::F64 rms_relerr_p {0.0};
    PS::F64 rms_relerr_a {0.0};
    PS::F64 relerr_p_max {0.0};
    PS::F64 relerr_a_max {0.0};
    for (PS::S32 i=0; i<n_loc; i++) {
        log_file << "i = " << i 
                 << " rank = " << PS::Comm::getRank()
                 << " id = " << psys[i].getId() << std::endl;
        PS::F64 err;
        // Potential
        const PS::F64 pot = psys[i].pot - psys[i].pot_pm;
        const PS::F64 pot_exact = - psys_tmp[i].pot;
        log_file << "pot(pp) = " << psys[i].pot << std::endl;
        log_file << "pot(pm) = " << -psys[i].pot_pm << std::endl;
        log_file << "pot(tot.) = " << pot << std::endl;
        log_file << "pot_exact = " << pot_exact << std::endl;
        assert(pot_exact != 0.0);
        err = std::fabs(pot - pot_exact) / std::fabs(pot_exact);
        log_file << "err = " << err << std::endl;
        rms_relerr_p += (err * err);
        if (err > relerr_p_max) {
            relerr_p_max = err;
        }
        log_file << "rms_relerr_p = " << rms_relerr_p << std::endl;
        // Acceleration
        const PS::F64vec acc = psys[i].acc + psys[i].acc_pm;
        const PS::F64vec acc_exact = psys_tmp[i].acc;
        const PS::F64vec adiff = acc - acc_exact;
        const PS::F64vec anorm = acc_exact;
        log_file << "acc(pp) = " << psys[i].acc << std::endl;
        log_file << "acc(pm) = " << psys[i].acc_pm << std::endl;
        log_file << "acc(tot.) = " << acc << std::endl;
        log_file << "acc_exact = " << acc_exact << std::endl;
        assert(anorm.x * anorm.y * anorm.z != 0.0);
        err = std::sqrt(adiff * adiff)
            / std::sqrt(anorm * anorm);
        log_file << "err = " << err << std::endl;
        rms_relerr_a += (err * err);
        log_file << "rms_relerr_a = " << rms_relerr_a << std::endl;
        if (err > relerr_a_max) {
            relerr_a_max = err;
        }
        log_file << std::endl;
    }
    log_file.close();
    rms_relerr_p = std::sqrt(PS::Comm::getSum(rms_relerr_p)/n_glb);
    rms_relerr_a = std::sqrt(PS::Comm::getSum(rms_relerr_a)/n_glb);
    relerr_p_max = PS::Comm::getMaxValue(relerr_p_max);
    relerr_a_max = PS::Comm::getMaxValue(relerr_a_max);
    if (PS::Comm::getRank() == 0) {
        std::cout << "RMS relative error of potential    = " << rms_relerr_p << std::endl;
        std::cout << "(Maximum relative error = " << relerr_p_max << ")" << std::endl;
        std::cout << "RMS relative error of acceleration = " << rms_relerr_a << std::endl;
        std::cout << "(Maximum relative error = " << relerr_a_max << ")" << std::endl;
    }
    PS::Finalize();
    std::exit(0);
} 

template<class Tpsys>
void calc_energy(const Tpsys & psys,
                 PS::F64 & etot,
                 PS::F64 & ekin,
                 PS::F64 & epot,
                 const bool clear=true){
    if (clear) {
        etot = ekin = epot = 0.0;
    }
    PS::F64 etot_loc = 0.0;
    PS::F64 ekin_loc = 0.0;
    PS::F64 epot_loc = 0.0;
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    if (FP_grav::eps > 0.0) {
        for(PS::S32 i = 0; i < n_loc; i++){
            ekin_loc += psys[i].mass * psys[i].vel * psys[i].vel;
            epot_loc += psys[i].mass * ( psys[i].pot 
                                       + psys[i].mass / FP_grav::eps);
        }
    } else if (FP_grav::eps == 0.0) {
        for(PS::S32 i = 0; i < n_loc; i++){
            ekin_loc += psys[i].mass * psys[i].vel * psys[i].vel;
            epot_loc += psys[i].mass * psys[i].pot;
        }
    } else {
        assert(false);
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    etot_loc = ekin_loc + epot_loc;
    etot = PS::Comm::getSum(etot_loc);
    epot = PS::Comm::getSum(epot_loc);
    ekin = PS::Comm::getSum(ekin_loc);
}

int main(int argc, char *argv[])
{

    // Set output format
    std::cout << std::setprecision(15) << std::scientific;
    std::cerr << std::setprecision(15) << std::scientific;
   
    // Initialize FDPS 
    PS::Initialize(argc, argv);

    // Make an instance of ParticleSystem and initialize it 
    PS::ParticleSystem<FP_grav> psys;
    psys.initialize();

    // Make an instance of DomainInfo and initialize it
    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), 
                           PS::F64vec(1.0, 1.0, 1.0));

    // Read particle data
    run_param this_run;
    read_param_file(psys, this_run, argv[1]);

    // Perform domain decomposition & particle exchange
    dinfo.decomposeDomainAll(psys);
    psys.adjustPositionIntoRootDomain(dinfo);
    psys.exchangeParticle(dinfo);

    // Make an instance of TreeForForce and initialize it
    PS::S32 n_loc = psys.getNumberOfParticleLocal();
    //const PS::S32vec n_cell = PS::S32vec(8, 8, 8);
    //const PS::S32vec n_cell = PS::S32vec(16, 16, 16);
    const PS::S32vec n_cell = PS::S32vec(32, 32, 32);
    const PS::S32 icut = 2;
    const PS::U32 n_leaf_limit = 8;
    //const PS::U32 n_group_limit = 64;
    const PS::U32 n_group_limit = 128;
    //const PS::U32 n_group_limit = 256;
    PS::TreeForForceLong<Force_grav, EPI_grav, EPJ_grav>::QuadrupoleParticleMeshMultipole tree;
    tree.initialize(3*n_loc, n_cell, icut, this_run.theta, n_leaf_limit, n_group_limit);
    //tree.initialize(3*n_loc, n_cell, icut, 0.1, n_leaf_limit, n_group_limit);
    //tree.initialize(3*n_loc, n_cell, icut, 0.0, n_leaf_limit, n_group_limit);

    // Initialize Phantom-GRAPE library
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_open();
    g5_set_eps_to_all(FP_grav::eps);
#endif

    // Peform force calculation (PP part)
    //FP_grav::eps = 0.0; // for check
    tree.calcForceAllAndWriteBack(CalcGravity<EPJ_grav>(),
                                  CalcGravity<PS::SPJQuadrupolePMMM>(),
                                  psys, dinfo);

    // Make an instance of ParticleMeshMultipole class and initialize it
    //const PS::S32 p = 2;
    //const PS::S32 p = 3;
    const PS::S32 p = 4;
    //const PS::S32 p = 5;
    //const PS::S32 p = 7;
    const PS::S32 p_spj2mm = 5;
    PS::PMM::ParticleMeshMultipole<Force_grav, EPI_grav> pm;
    pm.initialize(p, p_spj2mm);

    // Perform force calculation (PM part)
    pm.calcForceAllAndWriteBack(tree, dinfo, psys);

    // Compare the calculated force with the exact solution (check)
    //check_error(psys);

    // Calculate total force
    calc_total_force(psys);

    // Calculate energy
    PS::F64 ekin0, epot0, etot0;
    calc_energy(psys, etot0, ekin0, epot0, true);
    if (PS::Comm::getRank() == 0) {
        std::cout << "-----------------------------------------" << std::endl;
        std::cout << " Energies at the start of the simulation:" << std::endl;
        std::cout << "    ekin0 = " << ekin0 << std::endl;
        std::cout << "    epot0 = " << epot0 << std::endl;
        std::cout << "    etot0 = " << etot0 << std::endl;
        std::cout << "-----------------------------------------" << std::endl;
    }
    //PS::Finalize();
    //std::exit(0);

    // Time integration
    PS::F64 dtime = calc_dtime(psys, this_run);
    PS::F64 dtime_prev, dtime_mid;
    this_run.output_diag(dtime);
    drift(psys, dinfo, 0.5*dtime);
    dinfo.decomposeDomainAll(psys);
    psys.adjustPositionIntoRootDomain(dinfo);    
    psys.exchangeParticle(dinfo);
    this_run.npart_local = psys.getNumberOfParticleLocal();

    // Output for debugging
#if 0
    std::stringstream ss;
    ss << "ic" << std::setw(5) << std::setfill('0') << PS::Comm::getRank() << ".dat";
    const std::string filename = ss.str();
    std::ofstream output_file;
    output_file.open(filename.c_str(), std::ios::trunc|std::ios::binary);
    for(PS::S32 i = 0; i < psys.getNumberOfParticleLocal(); i++) {
        output_file.write((char *)&psys[i].id, sizeof(PS::S64));
        output_file.write((char *)&psys[i].mass, sizeof(PS::F64));
        output_file.write((char *)&psys[i].pos, 3 * sizeof(PS::F64));
    }
    output_file.close();
    PS::Finalize();
    std::exit(0);
#endif
    this_run.step = 0;
    while(this_run.znow > this_run.zend) {
        // Output current time
        if (PS::Comm::getRank() == 0) {
            std::cout <<  "this_run.step=" << this_run.step 
                      << " this_run.znow=" << this_run.znow
                      << " this_run.zend=" << this_run.zend
                      << " this_run.tnow=" << this_run.tnow
                      << " dtime=" << dtime
                      << std::endl;
        }

        // Perform force calculation
        PS::F64 t_start, t_mid, t_end;
        PS::Comm::barrier(); t_start = PS::GetWtime();
        //if (this_run.step == 90) FP_grav::eps = 0.0; // for debug
        tree.calcForceAllAndWriteBack(CalcGravity<EPJ_grav>(),
                                      CalcGravity<PS::SPJQuadrupolePMMM>(),
                                      psys, dinfo);
        PS::Comm::barrier(); t_mid = PS::GetWtime();
        pm.calcForceAllAndWriteBack(tree, dinfo, psys);
        PS::Comm::barrier(); t_end = PS::GetWtime();
        if (PS::Comm::getRank() == 0) {
            std::cout << "t_PP = " << (t_mid - t_start)
                      << " t_PM = " << (t_end - t_mid)
                      << std::endl;
        }
        // --- for debug [start] ---
        //if (this_run.step == 90) {
        //    check_error(psys);
        //}
        // --- for debug [end] ---
        calc_total_force(psys);

        // Update time, etc.
        this_run.tnow += 0.5*dtime;
        this_run.update_expansion(this_run.tnow);
     
#if 0 
        // Check energy errors
        if (this_run.step % 100 == 0) {
            n_loc = psys.getNumberOfParticleLocal();
            // (i) Save vel
            std::vector<PS::F64vec> vel_save;
            vel_save.resize(n_loc);
            for (PS::S32 i=0; i<n_loc; i++) vel_save[i] = psys[i].vel;
            // (ii) K(dt/2)
            kick(psys, 0.5*dtime, this_run);
            // (iii) Calculate energies
            PS::F64 ekin, epot, etot;
            calc_energy(psys, etot, ekin, epot, true);
            if (PS::Comm::getRank() == 0) {
                std::cout << "-------------------------------------------------" << std::endl;
                std::cout << " Energy errors:" << std::endl;
                std::cout << "    Rel.err. (etot) = " << std::abs((etot-etot0)/etot0) << std::endl;
                std::cout << "    Rel.err. (ekin) = " << std::abs((ekin-ekin0)/ekin0) << std::endl;
                std::cout << "    Abs.err. (epot) = " << std::abs((epot-epot0)) << std::endl;
                std::cout << "    etot  = " << etot << std::endl;
                std::cout << "    etot0 = " << etot0 << std::endl;
                std::cout << "    ekin  = " << ekin << std::endl;
                std::cout << "    ekin0 = " << ekin0 << std::endl;
                std::cout << "    epot  = " << epot << std::endl;
                std::cout << "    epot0 = " << epot0 << std::endl;
                std::cout << " Note that epot0 in the exact solution must be 0 " << std::endl;
                std::cout << " because of the periodic B.C." << std::endl;
                std::cout << "-------------------------------------------------" << std::endl;
            }
            // (iv) Restore vel
            for (PS::S32 i=0; i<n_loc; i++) psys[i].vel = vel_save[i];
        }
#endif

        // Kick
        kick(psys, dtime, this_run);

        // Update time again
        this_run.tnow += 0.5*dtime;
        this_run.update_expansion(this_run.tnow);

        // Calculate time step
        dtime_prev = dtime;
        dtime = calc_dtime(psys, this_run);
        dtime_mid = 0.5*(dtime_prev + dtime);

        // Drift
        drift(psys, dinfo, dtime_mid);

        // Perforce domain decomposition and particle exchange
        dinfo.decomposeDomainAll(psys);
        psys.adjustPositionIntoRootDomain(dinfo);
        psys.exchangeParticle(dinfo);

        // Output
        this_run.npart_local = psys.getNumberOfParticleLocal();
        output_data_in_run(psys, this_run);
        this_run.step++;
        this_run.output_diag(dtime);

        // For test
        //if (this_run.step == 5) break;
        if (this_run.step == 100) break;
        //if (this_run.step == 500) break;
    }
    // Finalize Phantom-GRAPE library
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_close();
#endif
    // Finalize FDPS
    PS::Finalize();
    return 0;
}
