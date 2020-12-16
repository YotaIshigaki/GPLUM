#pragma once
#include <particle_simulator.hpp>
#include "Particle_Class.h"
#include "SPH_Objects_Class.h"
#include "Misc_Class.h"

extern void Initialize(int argc, char* argv[], 
                       SPH_Objects& SPH_objs,
                       Fluctuation_Monitor& f_monitor);
extern void Finalize();
