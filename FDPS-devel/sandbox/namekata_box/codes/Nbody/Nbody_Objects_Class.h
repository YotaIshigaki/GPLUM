#pragma once
#include <particle_simulator.hpp>
#include "Particle_Class.h"

//================================
//* Class Decl.: Nbody_Objects
//================================
class Nbody_Objects {
   public:
      PS::ParticleSystem<Nbody_FP> system;
      PS::DomainInfo dinfo;
      // Constructors
      // Methods for Nbody simulation
      void init_tree();
      void calc_gravity();
      PS::F64 calc_timestep();
};

