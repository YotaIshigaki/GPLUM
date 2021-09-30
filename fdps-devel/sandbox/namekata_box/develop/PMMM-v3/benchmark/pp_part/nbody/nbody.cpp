#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>
#include "user-defined.hpp"

void makeColdUniformCube(const PS::F64 mass_glb,
                         const PS::S64 n_loc,
                         PS::F64 *& mass,
                         PS::F64vec *& pos,
                         PS::F64vec *& vel) {
    const PS::S32 n_proc  = PS::Comm::getNumberOfProc(); 
    const PS::S32 my_rank = PS::Comm::getRank();
    const PS::S64 n_glb = PS::Comm::getSum(n_loc);
    PS::MTTS mt;
    mt.init_genrand(my_rank);
    for(PS::S32 i = 0; i < n_loc; i++){
        mass[i] = mass_glb / n_glb;
        pos[i][0] = mt.genrand_res53();
        pos[i][1] = mt.genrand_res53();
        pos[i][2] = mt.genrand_res53();
        vel[i][0] = 0.0;
        vel[i][1] = 0.0;
        vel[i][2] = 0.0;
    }

    PS::F64vec cm_pos  = 0.0;
    PS::F64vec cm_vel  = 0.0;
    PS::F64    cm_mass = 0.0;
    for(PS::S32 i = 0; i < n_loc; i++){
        cm_pos  += mass[i] * pos[i];
        cm_vel  += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos = PS::Comm::getSum(cm_pos)/cm_mass;
    cm_vel = PS::Comm::getSum(cm_vel)/cm_mass;
    for(PS::S32 i = 0; i < n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }
}

template<class Tpsys>
void setParticlesColdUniformCube(Tpsys & psys,
                                 const PS::S32 n_loc) {

    psys.setNumberOfParticleLocal(n_loc);
    PS::F64    * mass = new PS::F64[n_loc];
    PS::F64vec * pos  = new PS::F64vec[n_loc];
    PS::F64vec * vel  = new PS::F64vec[n_loc];
    const PS::F64 m_tot = 1.0;
    makeColdUniformCube(m_tot, n_loc, mass, pos, vel);
    for(PS::S32 i = 0; i < n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos  = pos[i];
        psys[i].vel  = vel[i];
        psys[i].id   = i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}



PS::F64 FPGrav::eps = 1.0/32.0;

int main(int argc, char *argv[]) {
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

    PS::Initialize(argc, argv);
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::S32 n_thrd = PS::Comm::getNumberOfThread();

    PS::F32 theta = 0.5;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::F32 time_end = 10.0;
    PS::F32 dt = 1.0 / 128.0;
    PS::F32 dt_diag = 1.0 / 8.0;
    PS::F32 dt_snap = 1.0;
    //PS::S64 n_ptcl_per_proc = 1000000; // 1M
    //PS::S64 n_ptcl_per_proc = 100000; // 100k
    //PS::S64 n_ptcl_per_proc = 10000; // 10k
    PS::S64 n_ptcl_per_proc = 1000; // 1k
    char dir_name[1024];
    sprintf(dir_name,"./result");

    // Make an instance of ParticleSystem class and set initial condition
    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
    setParticlesColdUniformCube(system_grav, n_ptcl_per_proc);

    // Make an instance of DomainInfo class and perform domain decomposition
    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.decomposeDomainAll(system_grav);
    system_grav.exchangeParticle(dinfo);
    PS::S32 n_loc = system_grav.getNumberOfParticleLocal();
   
    // Make an instance of TreeForForce class and perform force calculation 
    PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::QuadrupoleGeometricCenter tree_grav;
    tree_grav.initialize(3*n_loc, theta, n_leaf_limit, n_group_limit);
    tree_grav.calcForceAllAndWriteBack(CalcGravityEpEp,
                                       CalcGravityEpSp,
                                       system_grav,
                                       dinfo);

    // Measure the performance
    PS::F64 etime_start, etime_end;
    PS::F64 etime = 0;
    const PS::S32 n_trials = 64;
    for (PS::S32 i=0; i<n_trials; i++) {
        PS::Comm::barrier();
        etime_start = PS::GetWtime();
        tree_grav.calcForceAllAndWriteBack(CalcGravityEpEp,
                                           CalcGravityEpSp,
                                           system_grav,
                                           dinfo);
        PS::Comm::barrier();
        etime_end = PS::GetWtime();
        etime += (etime_end - etime_start);
    }
    const PS::F64 etime_1step = etime/n_trials;
    if (PS::Comm::getRank() == 0) {
        std::cout << "------------ Summary ------------" << std::endl;
        std::cout << "n_proc          = " << n_proc << std::endl;
        std::cout << "n_thrd          = " << n_thrd << std::endl;
        std::cout << "n_ptcl_per_proc = " << n_ptcl_per_proc << std::endl;
        std::cout << "n_tot           = " << n_ptcl_per_proc * n_proc << std::endl;
        std::cout << "theta           = " << theta << std::endl;
        std::cout << "n_trials        = " << n_trials << std::endl;
        std::cout << "etime           = " << etime << " [s]" << std::endl;
        std::cout << "etime_1step     = " << etime_1step << " [s]" << std::endl;
        std::cout << "---------------------------------" << std::endl;
    }
    

    PS::Finalize();
    return 0;
}
