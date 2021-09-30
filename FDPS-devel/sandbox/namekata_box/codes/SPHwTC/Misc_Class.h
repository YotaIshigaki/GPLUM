#pragma once
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "Particle_Class.h"

//======================================
//* Class Decl.: Fluctuation_Monitor
//======================================
class Fluctuation_Monitor {
   private:
      PS::F64 check_time,check_interval;
   public:
      // Constructors
      Fluctuation_Monitor();
      // Members
      void set_config(PS::F64 start_time, PS::F64 interval);
      void dump_velocity(PS::ParticleSystem<SPH_FP>& system, PS::F64 dt);
      void check(PS::ParticleSystem<SPH_FP>& system,
                 PS::F64 time, PS::F64 stop_time);
};
