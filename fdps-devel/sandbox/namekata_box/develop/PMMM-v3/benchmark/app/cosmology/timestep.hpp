#pragma once
// Include header files of C++ library
#include <limits>
// Include header files of FDPS
#include <particle_simulator.hpp>
// Include user-defined headers
#include "cosmology.hpp"

template <class Tpsys>
PS::F64 calc_dtime(Tpsys & psys, run_param & this_run) {
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    PS::F64 dtime = std::numeric_limits<PS::F64>::max();
    for (PS::S64 i = 0; i < n_loc; i++) {
        dtime = fmin(dtime, psys[i].calcDtime(this_run));
    }

    if (this_run.noutput >= 1) {
        COSM::real_t zred_next;
        COSM::real_t zred_next_output;
        COSM::real_t time_next_output;

        if (this_run.znow < this_run.output_timing[this_run.output_indx] && 
            this_run.output_indx+1 < this_run.noutput) {
            zred_next_output = this_run.output_timing[this_run.output_indx+1];
        } else {
            zred_next_output = this_run.output_timing[this_run.output_indx];
        }
        zred_next = this_run.cosm.timetoz(this_run.tnow+dtime);

        if (zred_next < zred_next_output){
            time_next_output = this_run.cosm.ztotime(zred_next_output/1.0001);
            dtime = time_next_output - this_run.tnow;
        }
    }
    return PS::Comm::getMinValue(dtime);
}

