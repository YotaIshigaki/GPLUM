#pragma once
#include "particle_simulator.hpp"

namespace simulation_units {

extern double unit_leng;
extern double unit_mass;
extern double unit_velc;
extern double unit_dens;
extern double unit_time;
extern double unit_engy;
extern double unit_pres;
extern double unit_temp;
extern double unit_rate;
extern double unit_accl;

extern void setup_sim_units();
extern void read_sim_units();

}
