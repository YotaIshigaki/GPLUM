#pragma once
#include <particle_simulator.hpp>
#include "Particle_Class.h"
#include "Nbody_Objects_Class.h"

extern void Initialize(int argc, char* argv[], 
                       Nbody_Objects& Nbody_objs);
extern void Finalize();
