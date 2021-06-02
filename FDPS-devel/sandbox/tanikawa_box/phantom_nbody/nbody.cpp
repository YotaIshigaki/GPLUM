#include <gp5util.h>

#include <particle_simulator.hpp>
#include "user-defined.hpp"

template <class Tpsys>
void setParticleColdUniformSphere(Tpsys & psys,
                                  const PS::S32 n_glb) {

    PS::S32 rank  = PS::Comm::getRank();
    PS::S32 n_loc = (rank == 0) ? n_glb : 0;
    psys.setNumberOfParticleLocal(n_loc);

    PS::MT::init_genrand(rank);
    for(PS::S32 i = 0; i < n_loc; i++) {
        psys[i].mass = 1.0 / (PS::F32)n_glb;
        const PS::F64 radius = 3.0;
        do {
            psys[i].pos[0] = (2. * PS::MT::genrand_res53() - 1.) * radius;
            psys[i].pos[1] = (2. * PS::MT::genrand_res53() - 1.) * radius;
            psys[i].pos[2] = (2. * PS::MT::genrand_res53() - 1.) * radius;
        }while(psys[i].pos * psys[i].pos >= radius * radius);
        psys[i].vel[0] = 0.0;
        psys[i].vel[1] = 0.0;
        psys[i].vel[2] = 0.0;        
    }
}

template<class Tpsys>
void predict(Tpsys & system,
             const PS::F64 dt) {
    PS::S32 n_loc = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < n_loc; i++) {
        system[i].predict(dt);
    }
}

template<class Tpsys>
void correct(Tpsys & system,
             const PS::F64 dt) {
    PS::S32 n_loc = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < n_loc; i++) {
        system[i].correct(dt);
    }
}

template<class Tpsys>
PS::F64 calcEnergy(const Tpsys & system) {

    PS::F64 etot = 0.0;
    PS::F64 etot_loc = 0.0;
    PS::F64 ekin_loc = 0.0;
    PS::F64 epot_loc = 0.0;

    const PS::S32 n_loc = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < n_loc; i++){
        ekin_loc += system[i].mass * system[i].vel * system[i].vel;
        epot_loc += system[i].mass * (system[i].pot + system[i].mass / FPGrav::eps);
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    etot_loc  = ekin_loc + epot_loc;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    etot = PS::Comm::getSum(etot_loc);
#else
    etot = etot_loc;
#endif

    return etot;
}

int main(int argc, char *argv[]) {
    PS::F32 time  = 0.0;
    PS::F32 tend  = 10.0;
    PS::F32 dtime = 1.0 / 128.0;
    PS::F32 dtout = 1.0 / 8.0;
    PS::S64 ntot  = 1024;
    
    PS::Initialize(argc, argv);

    g5_open();
    g5_set_eps_to_all(FPGrav::eps);
    
    PS::DomainInfo dinfo;
    dinfo.initialize();

    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();

    PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::
        Monopole tree_grav;
    tree_grav.initialize(ntot);

    setParticleColdUniformSphere(system_grav, ntot);

    dinfo.decomposeDomainAll(system_grav);

    system_grav.exchangeParticle(dinfo);
    
    tree_grav.calcForceAllAndWriteBack
        (CalcGravity<FPGrav>(),
         CalcGravity<PS::SPJMonopole>(),
         system_grav,
         dinfo);
    
    PS::F64 etot0 = calcEnergy(system_grav);
    if(PS::Comm::getRank() == 0) {
        fprintf(stderr, "time: %10.7f energy: %+e energy error: %+e\n",
                time, etot0, (etot0 - etot0) / etot0);
    }

    while(time < tend) {

        predict(system_grav, dtime);        
        dinfo.decomposeDomainAll(system_grav);
        system_grav.exchangeParticle(dinfo);        
        tree_grav.calcForceAllAndWriteBack
            (CalcGravity<FPGrav>(),
             CalcGravity<PS::SPJMonopole>(),
             system_grav,
             dinfo);
        correct(system_grav, dtime);
        
        time += dtime;
        PS::F64 etot1 = calcEnergy(system_grav);
        if(fmod(time, dtout) == 0.0 &&
           PS::Comm::getRank() == 0) {
                fprintf(stderr, "time: %10.7f energy: %+e energy error: %+e\n",
                        time, etot1, (etot1 - etot0) / etot0);
        }
                
    }

    system_grav.writeParticleAscii("result.snp");

    g5_close();
        
    PS::Finalize();

    return 0;
}
