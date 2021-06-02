#pragma once
#include <particle_simulator.hpp>

//=====================================
//* Class Decl.: SPH_Gather_Results
//=====================================
class SPH_Gather_Results {
   public:
      PS::S32 flag;
      PS::F64 rho;
      PS::F64 h;
      PS::F64 vsig;
      PS::F64 divv,Dtdivv,trSS,trWW;
      PS::F64vec rotv;
      PS::F64 fAV,alpha_loc;
      PS::F64 f_grdh;
      void clear();
};

//======================================
//* Class Decl.: SPH_Symmetry_Results
//======================================
class SPH_Symmetry_Results {
   public:
      PS::F64vec ap,avis;
      PS::F64 gamad,gamvis;
      PS::F64 vsig;
      void clear();
};

//======================================
//* Class Decl.: SPH_FP 
//======================================
class SPH_FP
{
   public:
      //---------------------------------------------------------------
      // The definitions of member variables 
      // 
      //  id      : ID of SPH particles
      //  m       : particle mass
      //  h       : SPH kernel length
      //  x       : position of particle
      //  v       : velocity of particle
      //  v_half  : velocity of particle at t=t^{n+1/2}
      //  ap      : acceleration by pressure gradient
      //  avis    : acceleration by viscousity
      //  rho     : mass density
      //  pres    : pressure
      //  u       : specific internal energy
      //  u_half  : specific internal energy at t=t^{n+1/2}
      //  Tgas    : gas temperature
      //  Tgr     : grain temperature
      //  mu      : mean molecular weight
      //  gameff  : effective specific heat ratio
      //  cs      : sound speed
      //  gamad   : adiabatic compressional heating rate
      //  gamvis  : artificial viscosity heating
      //  alpha   : variable \alpha by Morris (1997)
      //  fAV     : the so-called Balsara switch
      //---------------------------------------------------------------
      PS::S64 id;
      PS::S32 flag;
	   PS::F64 m; 
	   PS::F64 h,h_prev;
	   PS::F64vec x; 
	   PS::F64vec v,v_half; 
	   PS::F64vec ap,avis,atot; 
      PS::F64 rho,pres,u,u_half;
      PS::F64 Tgas,Tgr,mu,gameff,cs;
      PS::F64 gamad,gamvis;
      PS::F64 divv,rotv,Dtdivv,trSS,trWW;
      PS::F64 alpha,alpha_loc;
      PS::F64 fAV;
      PS::F64 vsig;
      PS::F64 f_grdh;
      // Constructor & destructor
      SPH_FP();
      ~SPH_FP();
      // Member functions required by FDPS
	   PS::F64 getCharge() const;
      PS::F64 getChargeParticleMesh() const;
	   PS::F64vec getPos() const;
	   PS::F64 getRSearch() const;
	   void setPos(const PS::F64vec& x);
      void copyFromForce(const SPH_Gather_Results& result);
      void copyFromForce(const SPH_Symmetry_Results& result);
      void copyFromForceParticleMesh(const PS::F64 apm);
	   void writeAscii(FILE* fp) const;
	   void readAscii(FILE* fp);
      // The other member functions
      void calcPressure();
};

//=========================================
//* Class Decl.: SPH_EP 
//=========================================
class SPH_EP
{
   public:
      PS::S64 id;
	   PS::F64 m;
	   PS::F64 h,h_prev;
	   PS::F64vec x;
	   PS::F64vec v;
      PS::F64vec ap,avis,atot;
	   PS::F64 rho;
	   PS::F64 pres;
	   PS::F64 cs;
      PS::F64 divv,Dtdivv,trSS,trWW;
      PS::F64 alpha;
      PS::F64 fAV;
      PS::F64 vsig;
      PS::F64 f_grdh;
      // Constructor & destructor
      SPH_EP();
      ~SPH_EP();
      // Member functions required by FDPS
	   PS::F64vec getPos() const;
	   PS::F64 getRSearch() const;
	   void setPos(const PS::F64vec& x);
      void copyFromFP(const SPH_FP& FP);
      // The other member functions
};

//=========================================
//* Class Decl.: SPH_IO
//=========================================
class SPH_IO
{
   public:
      PS::S64 id;
      PS::F64 m;
      PS::F64 h;
      PS::F64vec x;
      PS::F64vec v;
      PS::F64vec ap,avis;
      PS::F64 u,mu,gameff;
      PS::F64 gamad,gamvis;
      PS::F64 alpha;
      // Constructors & destructor
      // Member functions
      void copyFromFP(const SPH_FP& FP);
      void copyToFP(SPH_FP& FP) const;
};


//=============================
//* Class Decl.: Calc_density
//=============================
class Calc_density{
   public:
      void operator () (const SPH_EP* const ep_i,
                        const PS::S32 Nip,
                        const SPH_EP* const ep_j,
                        const PS::S32 Njp,
                        SPH_Gather_Results* const result);
};

//================================
//* Class Decl.: Calc_divv
//================================
class Calc_divv{
   public:
      void operator () (const SPH_EP* const ep_i,
                        const PS::S32 Nip,
                        const SPH_EP* const ep_j,
                        const PS::S32 Njp,
                        SPH_Gather_Results* const result);
};

//================================
//* Class Decl.: Calc_rotv
//================================
class Calc_rotv{
   public:
      void operator () (const SPH_EP* const ep_i,
                        const PS::S32 Nip,
                        const SPH_EP* const ep_j,
                        const PS::S32 Njp,
                        SPH_Gather_Results* const result);
};

//=====================================
//* Class Decl.: Calc_vsig_S05
//=====================================
class Calc_vsig_S05{
   public:
      void operator () (const SPH_EP* const ep_i,
                        const PS::S32 Nip,
                        const SPH_EP* const ep_j,
                        const PS::S32 Njp,
                        SPH_Symmetry_Results* const result);
};

//================================
//* Class Decl.: Calc_vsig_CD10
//================================
class Calc_vsig_CD10{
   public:
      void operator () (const SPH_EP* const ep_i,
                        const PS::S32 Nip,
                        const SPH_EP* const ep_j,
                        const PS::S32 Njp,
                        SPH_Gather_Results* const result);
};

//================================
//* Class Decl.: Calc_RoC_velcfld
//================================
class Calc_RoC_velcfld{
   public:
      void operator () (const SPH_EP* const ep_i,
                        const PS::S32 Nip,
                        const SPH_EP* const ep_j,
                        const PS::S32 Njp,
                        SPH_Gather_Results* const result);
};

//================================
//* Class Decl.: Calc_alpha_loc
//================================
class Calc_alpha_loc{
   public:
      void operator () (const SPH_EP* const ep_i,
                        const PS::S32 Nip,
                        const SPH_EP* const ep_j,
                        const PS::S32 Njp,
                        SPH_Gather_Results* const result);
};


//================================
//* Class Decl.: Calc_fAV
//================================
class Calc_fAV{
   public:
      void operator () (const SPH_EP* const ep_i,
                        const PS::S32 Nip,
                        const SPH_EP* const ep_j,
                        const PS::S32 Njp,
                        SPH_Gather_Results* const result);
};

//================================
//* Class Decl.: Calc_gradh_term
//================================
class Calc_gradh_term{
   public:
      void operator () (const SPH_EP* const ep_i,
                        const PS::S32 Nip,
                        const SPH_EP* const ep_j,
                        const PS::S32 Njp,
                        SPH_Gather_Results* const result);
};

//================================
//* Class Decl.: Calc_hydro_force
//================================
class Calc_hydro_force{
   public:
      void operator () (const SPH_EP* const ep_i,
                        const PS::S32 Nip,
                        const SPH_EP* const ep_j,
                        const PS::S32 Njp,
                        SPH_Symmetry_Results* const result);
};


//=====================================
//* Class Decl.: Calc_kernel_length
//=====================================
class Calc_kernel_length{
   public:
      void operator () (const SPH_EP* const ep_i,
                        const PS::S32 Nip,
                        const SPH_EP* const ep_j,
                        const PS::S32 Njp,
                        SPH_Gather_Results* const result);
};

