#pragma once
/* C++ headers */
#include "common.h"
/* FDPS header */
#include <particle_simulator.hpp>

/* Class for time measurement */
class TimeProfile {
public:
    PS::F64 calc_gravity_1st, calc_gravity_2nd;
    PS::F64 calc_knl_sz_1st, calc_knl_sz_2nd;
    PS::F64 calc_density_1st, calc_density_2nd;
    PS::F64 calc_hydro_force_1st, calc_hydro_force_2nd;
    PS::F64 feedback;
    PS::F64 cooling;
    PS::F64 star_formation;

    TimeProfile();

    void clear();
    PS::F64 getTotalTime() const;
    void dump(std::ostream & fout = std::cout) const;
    TimeProfile operator + (const TimeProfile & rhs) const;
};

void barrier();

extern TimeProfile time_prof;
