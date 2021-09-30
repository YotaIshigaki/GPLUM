/* C++ header */
#include "common.h"
/* FDPS header */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "timing_utilities.h"

/* Class for time measurement */
TimeProfile::TimeProfile() {
    this->clear();
}
void TimeProfile::clear () {
    this->calc_gravity_1st = 0.0;
    this->calc_gravity_2nd = 0.0;
    this->calc_knl_sz_1st = 0.0;
    this->calc_knl_sz_2nd = 0.0;
    this->calc_density_1st = 0.0;
    this->calc_density_2nd = 0.0;
    this->calc_hydro_force_1st = 0.0;
    this->calc_hydro_force_2nd = 0.0;
    this->feedback = 0.0;
    this->cooling = 0.0;
    this->star_formation = 0.0;
}
PS::F64 TimeProfile::getTotalTime() const {
    const PS::F64 sum = this->calc_gravity_1st 
                      + this->calc_gravity_2nd
                      + this->calc_knl_sz_1st
                      + this->calc_knl_sz_2nd
                      + this->calc_density_1st
                      + this->calc_density_2nd
                      + this->calc_hydro_force_1st
                      + this->calc_hydro_force_2nd;
                      + this->feedback
                      + this->cooling
                      + this->star_formation;
    return sum;
}
void TimeProfile::dump(std::ostream & fout) const {
    fout << "total_time= " << this->getTotalTime() << std::endl;
    fout << "  calc_gravity_1st = " << this->calc_gravity_1st << std::endl;
    fout << "  calc_gravity_2nd = " << this->calc_gravity_2nd << std::endl;
    fout << "  calc_knl_sz_1st = " << this->calc_knl_sz_1st << std::endl;
    fout << "  calc_knl_sz_2nd = " << this->calc_knl_sz_2nd << std::endl;
    fout << "  calc_density_1st = " << this->calc_density_1st << std::endl;
    fout << "  calc_density_2nd = " << this->calc_density_2nd << std::endl;
    fout << "  calc_hydro_force_1st = " << this->calc_hydro_force_1st << std::endl;
    fout << "  calc_hydro_force_2nd = " << this->calc_hydro_force_2nd << std::endl;
    fout << "  feedback = " << this->feedback << std::endl;
    fout << "  cooling = " << this->cooling << std::endl;
    fout << "  star_formation = " << this->star_formation << std::endl;
}
TimeProfile TimeProfile::operator + (const TimeProfile & rhs) const {
    TimeProfile ret;
    ret.calc_gravity_1st = this->calc_gravity_1st + rhs.calc_gravity_1st;
    ret.calc_gravity_2nd = this->calc_gravity_2nd + rhs.calc_gravity_2nd;
    ret.calc_knl_sz_1st = this->calc_knl_sz_1st + rhs.calc_knl_sz_1st;
    ret.calc_knl_sz_2nd = this->calc_knl_sz_2nd + rhs.calc_knl_sz_2nd;
    ret.calc_density_1st = this->calc_density_1st + rhs.calc_density_1st;
    ret.calc_density_2nd = this->calc_density_2nd + rhs.calc_density_2nd;
    ret.calc_hydro_force_1st = this->calc_hydro_force_1st + rhs.calc_hydro_force_1st;
    ret.calc_hydro_force_2nd = this->calc_hydro_force_2nd + rhs.calc_hydro_force_2nd;
    ret.feedback = this->feedback + rhs.feedback;
    ret.cooling = this->cooling + rhs.cooling;
    ret.star_formation = this->star_formation + rhs.star_formation;
}

void barrier() {
#ifdef ASURA_FDPS_ENABLE_BARRIER
    PS::Comm::barrier();
#endif
}

TimeProfile time_prof;
