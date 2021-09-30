#pragma once
#include <particle_simulator.hpp>
#include <particle_mesh.hpp>
#include "Particle_Class.h"

//================================
//* Class Decl.: Nbody_Objects
//================================
class Nbody_Objects {
   public:
      PS::ParticleSystem<Nbody_FP> system;
      PS::DomainInfo dinfo;
      PS::TreeForForceLong<Nbody_PP_Results, Nbody_EP, Nbody_EP>::MonopoleWithCutoff pp_tree;
      PS::PM::ParticleMesh pm;
      // Constructors
      // Methods for Nbody simulation
      void init_tree();
      void calc_gravity();
      PS::F64 calc_timestep();
};

