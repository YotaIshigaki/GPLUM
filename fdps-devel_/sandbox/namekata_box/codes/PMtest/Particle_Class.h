#pragma once
#include <particle_simulator.hpp>

//======================================
//* Class Decl.: Nbody_PP_Results
//======================================
class Nbody_PP_Results
{
   public:
      PS::F64 pot;
      PS::F64vec agrv;
      void clear();
};

//======================================
//* Class Decl.: Nbody_FP 
//======================================
class Nbody_FP
{
   public:
      //---------------------------------------------------------------
      // The definitions of member variables 
      // 
      //  id      : ID of Nbody particles
      //  m       : particle mass
      //  eps     : gravitational softening
      //  rc      : cutoff radius
      //  x       : position of particle
      //  v       : velocity of particle
      //  v_half  : velocity of particle at t=t^{n+1/2}
      //  agrv    : acceleration by gravity
      //  pot     : gravitational potential
      //---------------------------------------------------------------
      PS::S64 id;
	   PS::F64 m;
      PS::F64 eps;
      PS::F64 rc;
	   PS::F64vec x; 
	   PS::F64vec v,v_half; 
	   PS::F64vec agrv; 
      PS::F64 pot;
      // Constructor & destructor
      Nbody_FP();
      ~Nbody_FP();
      // Member functions required by FDPS
	   PS::F64 getCharge() const;
      PS::F64 getChargeParticleMesh() const;
	   PS::F64vec getPos() const;
      PS::F64 getRSearch() const;
	   void setPos(const PS::F64vec& x);
      void copyFromForce(const Nbody_PP_Results& result);
      void copyFromForceParticleMesh(const PS::F64 apm);
};

//=========================================
//* Class Decl.: Nbody_EP 
//=========================================
class Nbody_EP
{
   public:
      PS::S64 id;
	   PS::F64 m;
      PS::F64 eps;
      PS::F64 rc;
	   PS::F64vec x;
      // Constructor & destructor
      Nbody_EP();
      ~Nbody_EP();
      // Member functions required by FDPS
      PS::F64 getCharge() const;
	   PS::F64vec getPos() const;
      PS::F64 getRSearch() const;
	   void setPos(const PS::F64vec& x);
      void copyFromFP(const Nbody_FP& FP);
      // The other member functions
};

//================================
//* Class Decl.: Calc_gravity
//================================
class Calc_gravity{
   public:
      void operator () (const Nbody_EP* const ep_i,
                        const PS::S32 Nip,
                        const Nbody_EP* const ep_j,
                        const PS::S32 Njp,
                        Nbody_PP_Results* const result);
};
