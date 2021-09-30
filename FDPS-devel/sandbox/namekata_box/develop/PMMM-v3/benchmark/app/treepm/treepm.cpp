#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>
#include <float.h>
#include <cstdio>
#include <cstdlib>
#include <particle_simulator.hpp>
#include <particle_mesh.hpp>

#include "treepm.hpp"
#include "run_param.hpp"
#include "cosmology.hpp"
#include "mpi.h"

#include "prototype.h"

static PS::F64 EPS_FOR_PP = 1.0;

template <class Tp>
void setup_random_particle(Tp & ptcl,
                           run_param & this_run,
                           const PS::S32 npart_total)
{
    PS::S32 rank = PS::Comm::getRank();
    PS::S64 npart_local = (rank == 0) ? npart_total : 0;

    ptcl.setNumberOfParticleLocal(npart_local);

    this_run.npart_total = npart_total;
    this_run.mpi_nproc = PS::Comm::getNumberOfProc();
    this_run.mpi_rank = rank;

    PS::MT::init_genrand(0);

    for (PS::S32 i=0;i<npart_local;i++) {
        ptcl[i].mass = 3.0*this_run.cosm.omegam/(8.0*AC::PI*(PS::F32)npart_total);
        ptcl[i].pos[0] = PS::MT::genrand_res53();
        ptcl[i].pos[1] = PS::MT::genrand_res53();
        ptcl[i].pos[2] = PS::MT::genrand_res53();
        ptcl[i].vel[0] = 0.0;
        ptcl[i].vel[1] = 0.0;
        ptcl[i].vel[2] = 0.0;
        //ptcl[i].eps = 0.1/pow(npart_total, 1.0/3.0);
        ptcl[i].eps = EPS_FOR_PP;
    }
}


template <class Tp>
void read_SB_particle(Tp &ptcl, 
                      run_param &this_run, 
                      const char * input_file){
    ptcl.readParticleBinary(input_file);
    this_run.npart_total = ptcl.getNumberOfParticleGlobal();
    this_run.npart_local = ptcl.getNumberOfParticleLocal();
    this_run.mpi_nproc = PS::Comm::getNumberOfProc();
    this_run.mpi_rank = PS::Comm::getRank();
    for (PS::S32 i=0;i<this_run.npart_local;i++) {
       //ptcl[i].eps = 0.1/pow(this_run.npart_total, 1.0/3.0);
       ptcl[i].eps = EPS_FOR_PP;
    }
}

template<class Tptcl>
void read_param_file(Tptcl & ptcl,
                     run_param & this_run,
                     const char * input_file){
    std::cerr << "input_file=" << input_file << std::endl;
    FILE * param_fp;
    param_fp = fopen(input_file, "r");
    if (param_fp == NULL) {
       fprintf(stderr,
               "File %s not found in input_params.\n",
               input_file);
       PS::Abort();
       exit(1);
    }
    int ret = fscanf(param_fp, "%d", &this_run.mode);
    if (ret == EOF){
       fprintf(stderr, "Input error of mode in input_params()\n");
       PS::Abort();
       exit(1);
    }
    ret = fscanf(param_fp, "%lf", &EPS_FOR_PP);
    ret = fscanf(param_fp, "%lf", &this_run.theta);
    ret = fscanf(param_fp, "%f", &this_run.zend);
    if (PS::Comm::getRank() == 0){
       std::cerr << "this_run.theta=" << this_run.theta << std::endl;
       std::cerr << "this_run.zend=" << this_run.zend << std::endl;
    }
    if (this_run.mode == run_param::SANTABARBARA){
       this_run.znow = 63.0;
       this_run.cosm.omegam = 1.0;
       this_run.cosm.omegav = this_run.cosm.omegab = this_run.cosm.omeganu = 0.0;
       FPtreepm::H0 = 50.0;
       FPtreepm::Lbnd = 64.0;
       char snap_name[1024];
       fscanf(param_fp, "%s", snap_name);
       if (PS::Comm::getRank() == 0) std::cerr << "snap_name:" << snap_name << std::endl;
       read_SB_particle(ptcl, this_run, snap_name);
    }
    else if (this_run.mode == run_param::READ_FILE){
       char snap_name[1024];
       fscanf(param_fp, "%s", snap_name);
       if (PS::Comm::getRank() == 0) std::cerr << "snap_name:" << snap_name << std::endl;
       ret = fscanf(param_fp, "%f", &this_run.znow);
       ret = fscanf(param_fp, "%f", &this_run.cosm.omegam);
       ret = fscanf(param_fp, "%f", &this_run.cosm.omegav);
       ret = fscanf(param_fp, "%f", &this_run.cosm.omegab);
       ret = fscanf(param_fp, "%f", &this_run.cosm.omeganu);
       ret = fscanf(param_fp, "%lf", &FPtreepm::H0);
       ret = fscanf(param_fp, "%lf", &FPtreepm::Lbnd);
       read_SB_particle(ptcl, this_run, snap_name);    
    }
    else if (this_run.mode == run_param::RANDOM){
       ret = fscanf(param_fp, "%f", &this_run.znow);
       ret = fscanf(param_fp, "%f", &this_run.cosm.omegam);
       ret = fscanf(param_fp, "%f", &this_run.cosm.omegav);
       ret = fscanf(param_fp, "%f", &this_run.cosm.omegab);
       ret = fscanf(param_fp, "%f", &this_run.cosm.omeganu);
       ret = fscanf(param_fp, "%lld", &this_run.npart_total);
       setup_random_particle(ptcl, this_run, this_run.npart_total);    
    }
    else{
       PS::Abort();
       exit(1);
    }    
    this_run.anow = 1.0 / (1.0 + this_run.znow);
    this_run.tnow = this_run.cosm.atotime(this_run.anow);
    this_run.update_expansion(this_run.tnow);
    this_run.input_params(param_fp);
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
    for(PS::S32 i = 0; i < n_loc; i++){
        ekin_loc += psys[i].mass * psys[i].vel * psys[i].vel;
        epot_loc += psys[i].mass * ( psys[i].pot + psys[i].pot_pm 
                                   + psys[i].mass / psys[i].eps);
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    etot_loc = ekin_loc + epot_loc;
    etot = PS::Comm::getSum(etot_loc);
    epot = PS::Comm::getSum(epot_loc);
    ekin = PS::Comm::getSum(ekin_loc);    
}

int main(int argc, char **argv)
{
    // Set output format 
    std::cout << std::setprecision(15) << std::scientific;
    std::cerr << std::setprecision(15) << std::scientific;

    // Initialize FDPS
    PS::Initialize(argc, argv);

    // Make an instance of ParticleSystem class and initialize it
    PS::ParticleSystem<FPtreepm> ptcl;
    ptcl.initialize();

    // Read particle data
    run_param this_run;
    read_param_file(ptcl, this_run, argv[1]);

    // Make an instance of DomainInfo class and initialize it
    PS::DomainInfo domain_info;
    domain_info.initialize();
    domain_info.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    domain_info.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), 
                                 PS::F64vec(1.0, 1.0, 1.0));

    // Domain decomposition and particle exchange
    domain_info.decomposeDomainAll(ptcl);
    ptcl.adjustPositionIntoRootDomain(domain_info);
    ptcl.exchangeParticle(domain_info);
    this_run.npart_local = ptcl.getNumberOfParticleLocal();

    // Make an instance of TreeForForce class and initialize it
    PS::TreeForForceLong<Result_treepm, EPItreepm, EPJtreepm>::MonopoleWithCutoff treepm_tree;
    //treepm_tree.initialize(3*ptcl.getNumberOfParticleLocal(), this_run.theta);
    treepm_tree.initialize(3*ptcl.getNumberOfParticleLocal(), 0.0);

    // Initialize the Phantom-GRAPE library
#ifdef ENABLE_PHANTOM_GRAPE_X86
    //g5_open();
    pg5_gen_s2_force_table(EPS_FOR_PP, 3.0/SIZE_OF_MESH);
#endif

    // Perform interaction calculation
    treepm_tree.calcForceAllAndWriteBack
        (calc_pp_force<EPJtreepm>(),
         calc_pp_force<PS::SPJMonopoleCutoff>(),
         ptcl,
         domain_info);

    PS::PM::ParticleMesh pm;
    pm.calcForceAllAndWriteBack(ptcl, domain_info);
    for (PS::S32 i=0; i<ptcl.getNumberOfParticleLocal(); i++) {
        PS::F32vec pos = ptcl[i].pos;
        ptcl[i].pot_pm = pm.getPotential(pos);
    }

    // Calculate energy
    PS::F64 etot0, ekin0, epot0;
    calc_energy(ptcl, etot0, ekin0, epot0, true);
    if (PS::Comm::getRank() == 0) {
        std::cout << "--------------------------------" << std::endl;
        std::cout << " Result of energy calculation:" << std::endl;
        std::cout << "    ekin0 = " << ekin0 << std::endl;
        std::cout << "    epot0 = " << epot0 << std::endl;
        std::cout << "    etot0 = " << etot0 << std::endl;
        std::cout << "--------------------------------" << std::endl;
    }
    PS::Finalize();
    std::exit(0);

    // Calculate the 1st time step and drift particles by dtime/2
    PS::F64 dtime = calc_dtime(ptcl, this_run);
    PS::F64 dtime_prev, dtime_mid;
    this_run.output_diag(dtime);
    drift_ptcl(ptcl, domain_info, 0.5*dtime);

    // Domain decomposition and particle exchange
    domain_info.decomposeDomainAll(ptcl);
    ptcl.adjustPositionIntoRootDomain(domain_info);    
    ptcl.exchangeParticle(domain_info);

    // Update this_run and output some data
    this_run.npart_local = ptcl.getNumberOfParticleLocal();
    Map2D map_2d;
    PS::F64 res = 2.0;
    map_2d.initialize(ptcl, this_run, res);

    this_run.step = 0;    
    while(this_run.znow > this_run.zend) {
        // Output to stdout
        if (PS::Comm::getRank() == 0) {
            std::cout <<  "this_run.step=" << this_run.step 
                      << " this_run.znow=" << this_run.znow
                      << " this_run.zend=" << this_run.zend
                      << " dtime = " << dtime
                      << std::endl;
        }

        // Perform interaction calculation
        PS::F64 t_start, t_mid, t_end;
        PS::Comm::barrier(); t_start = PS::GetWtime();
        treepm_tree.calcForceAllAndWriteBack
            (calc_pp_force<EPJtreepm>(),
             calc_pp_force<PS::SPJMonopoleCutoff>(),
             ptcl,
             domain_info);
        PS::Comm::barrier(); t_mid = PS::GetWtime();
        pm.calcForceAllAndWriteBack(ptcl, domain_info);
        PS::Comm::barrier(); t_end = PS::GetWtime();
        if (PS::Comm::getRank() == 0) {
            std::cout << "t_PP = " << (t_mid - t_start)
                      << " t_PM = " << (t_end - t_mid)
                      << std::endl;
        }

        // Update time, etc.
        this_run.tnow += 0.5*dtime;
        this_run.update_expansion(this_run.tnow);
       
        // Kick 
        kick_ptcl(ptcl, dtime, this_run);

        // Update time, etc. again
        this_run.tnow += 0.5*dtime;
        this_run.update_expansion(this_run.tnow);

        // Calculate time step
        dtime_prev = dtime;
        dtime = calc_dtime(ptcl, this_run);
        dtime_mid = 0.5*(dtime_prev + dtime);

        // Drift
        drift_ptcl(ptcl, domain_info, dtime_mid);

        // Perform domain decomposition and particle exchange
        domain_info.decomposeDomainAll(ptcl);
        ptcl.adjustPositionIntoRootDomain(domain_info);
        ptcl.exchangeParticle(domain_info);

        // Output
        this_run.npart_local = ptcl.getNumberOfParticleLocal();
        //output_data_in_run(ptcl, this_run);
        output_data_in_run(ptcl, this_run, map_2d);
        this_run.step++;
        this_run.output_diag(dtime);
         
        if (this_run.step == 5) break;
    }

#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_close();
#endif    
    PS::Finalize();
    return(0);
}
