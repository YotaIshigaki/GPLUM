#pragma once
#include <particle_simulator.hpp>
#include "run_param.hpp"

template <class Tpsys>
void kick(Tpsys & psys, const PS::F64 dt, run_param & this_run) {
    const PS::S64 n_loc = psys.getNumberOfParticleLocal();
    const PS::F64 om = this_run.cosm.omegam;
    const PS::F64 ov = this_run.cosm.omegav;
    const PS::F64 anow = this_run.cosm.timetoa(this_run.tnow);
    const PS::F64 at = sqrt(1.e0+om*(1.e0/anow-1.e0)+ov*(SQR(anow)-1.e0))/anow;
    const PS::F64 bt = 1.0/CUBE(anow);
    const PS::F64 atdt1 = 1.0+at*dt;
    const PS::F64 vfact = (2.0-atdt1)/atdt1;
    const PS::F64 afact = bt*dt/atdt1;
    for (PS::S64 i = 0; i < n_loc; i++) {
        psys[i].vel = vfact * psys[i].vel + afact*psys[i].acc;
    }
}

template <class Tpsys>
void drift(Tpsys & psys, const PS::DomainInfo & dinfo, const PS::F64 dt) {
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    for (PS::S64 i = 0; i < n_loc; i++) psys[i].pos += psys[i].vel * dt;
    psys.adjustPositionIntoRootDomain(dinfo);
}
