#pragma once
#include <particle_simulator.hpp>
#include "Particle_Class.h"

//================================
//* Class Decl.: SPH_Objects
//================================
class SPH_Objects {
   public:
      PS::ParticleSystem<SPH_FP> system;
      PS::DomainInfo dinfo;
      PS::TreeForForceShort<SPH_Gather_Results, SPH_EP, SPH_EP>::Gather gat_tree;
      PS::TreeForForceShort<SPH_Symmetry_Results, SPH_EP, SPH_EP>::Symmetry sym_tree;
      // Constructors
      // Members for SPH simulation
      void init_tree();
      void update_h();
      void calc_rho();
      void calc_pres();
      void calc_hydroforce(PS::F64 dt);
      PS::F64 calc_timestep();
};

