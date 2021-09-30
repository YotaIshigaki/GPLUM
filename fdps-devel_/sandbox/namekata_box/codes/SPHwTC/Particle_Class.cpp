/* Standard headers */
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cfloat>
/* Libraries */
//#ifndef __GNUC__
//#  define  __attribute__(x)  /*NOTHING*/
//#endif
#include <Eigen/Dense>
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#include <omp.h>
#endif
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "preprocess_keywords.h"
#include "Math_Functions.hpp"
#include "Kernel_Functions.h"
#include "Particle_Class.h"
#include "Calculation_Conditions.h"

//* Constants
#if (Simulation_Dimension == 1)
static const PS::S32 dimSpace=1;
#elif (Simulation_Dimension == 2)
static const PS::S32 dimSpace=2;
#elif (Simulation_Dimension == 3)
static const PS::S32 dimSpace=3;
#endif

//* Prototype declaration
PS::F64 get_next_h(const PS::S32 nnb, const PS::F64 h_old);

//==================================
//* Class Defs.: SPH_Gather_Results
//==================================
void SPH_Gather_Results::clear() {
   flag      = 0;
   rho       = 0.0e0;
   h         = 0.0e0;
   vsig      = 0.0e0;
   divv      = 0.0e0;
   rotv      = 0.0e0;
   Dtdivv    = 0.0e0;
   trSS      = 0.0e0;
   trWW      = 0.0e0;
   fAV       = 0.0e0;
   alpha_loc = 0.0e0;
   f_grdh    = 0.0e0;
}

//======================================
//* Class Defs.: SPH_Symmetry_Results
//======================================
void SPH_Symmetry_Results::clear() {
   ap     = 0.0e0;
   avis   = 0.0e0;
   gamad  = 0.0e0;
   gamvis = 0.0e0;
   vsig   = 0.0e0;
}

//======================================
//* Class Defs.: SPH_FP 
//======================================
// Constructor & destructor
SPH_FP::SPH_FP() {};
SPH_FP::~SPH_FP() {};
// Member functions required by FDPS
PS::F64 SPH_FP::getCharge() const{
   return m;
}

PS::F64 SPH_FP::getChargeParticleMesh() const{
   return m;
}

PS::F64vec SPH_FP::getPos() const{
   return x;
}

PS::F64 SPH_FP::getRSearch() const{
   return h;
}
void SPH_FP::setPos(const PS::F64vec& x){
   this->x = x;
}
void SPH_FP::copyFromForce(const SPH_Gather_Results& result){}
void SPH_FP::copyFromForce(const SPH_Symmetry_Results& result){}
void SPH_FP::copyFromForceParticleMesh(const PS::F64 apm){}
// The other member functions 

//=========================================
//* Class Defs.: SPH_EP 
//=========================================
// Constructor & destructor
SPH_EP::SPH_EP() {};
SPH_EP::~SPH_EP() {};
// Member functions required by FDPS
PS::F64vec SPH_EP::getPos() const{
   return x;
}
PS::F64 SPH_EP::getRSearch() const{
   return h;
}
void SPH_EP::setPos(const PS::F64vec& x){
   this->x = x;
}
void SPH_EP::copyFromFP(const SPH_FP& FP){
   id     = FP.id;
   m      = FP.m;
   h      = FP.h;
   h_prev = FP.h_prev;
   x      = FP.x;
   v      = FP.v;
   ap     = FP.ap;
   avis   = FP.avis;
   atot   = FP.atot;
   rho    = FP.rho;
   pres   = FP.pres;
   cs     = FP.cs;
   divv   = FP.divv;
   Dtdivv = FP.Dtdivv;
   trSS   = FP.trSS;
   trWW   = FP.trWW;
   alpha  = FP.alpha;
   fAV    = FP.fAV;
   vsig   = FP.vsig;
   f_grdh = FP.f_grdh;
}

//=========================================
//* Class Defs.: SPH_EP 
//=========================================
void SPH_IO::copyFromFP(const SPH_FP& FP) {
   id     = FP.id;
   m      = FP.m;
   h      = FP.h;
   x      = FP.x;
   v      = FP.v;
   ap     = FP.ap;
   avis   = FP.avis;
   u      = FP.u;
   mu     = FP.mu;
   gameff = FP.gameff;
   gamad  = FP.gamad;
   gamvis = FP.gamvis;
   alpha  = FP.alpha;
}

void SPH_IO::copyToFP(SPH_FP& FP) const {
   FP.id     = id;
   FP.m      = m;
   FP.h      = h;
   FP.x      = x;
   FP.v      = v;
   FP.ap     = ap;
   FP.avis   = avis;
   FP.u      = u;
   FP.mu     = mu;
   FP.gameff = gameff;
   FP.gamad  = gamad;
   FP.gamvis = gamvis;
   FP.alpha  = alpha;
}

//=============================
//* Class Defs.: Calc_density
//=============================
void Calc_density::operator () (const SPH_EP* const ep_i,
                                const PS::S32 Nip,
                                const SPH_EP* const ep_j,
                                const PS::S32 Njp,
                                SPH_Gather_Results* const result) {
   using namespace kernel_functions;
   //* Local variables
   PS::F64vec dx;
   PS::F64 rij;

   for (PS::S32 i=0; i<Nip; i++){
      for (PS::S32 j=0; j<Njp; j++){
          dx = ep_j[j].x - ep_i[i].x;
          rij = std::sqrt(dx*dx);
          result[i].rho += ep_j[j].m * W_SPH(rij,ep_i[i].h);
      }
   }
}

//=======================================
//* Class Defs.: Calc_divv
//=======================================
void Calc_divv::operator ()(const SPH_EP* const ep_i,
                            const PS::S32 Nip,
                            const SPH_EP* const ep_j,
                            const PS::S32 Njp,
                            SPH_Gather_Results* const result){
   using namespace kernel_functions;
   //* Local variables
   PS::F64vec dx,dv;
   PS::F64 rij,xdotv,dWi;

   for (PS::S32 i=0; i<Nip; i++) {
      for (PS::S32 j=0; j<Njp; j++) {
         dx = ep_i[i].x - ep_j[j].x;
  			rij = std::sqrt(dx * dx);
         if ((ep_i[i].id == ep_j[j].id) && (rij == 0.0e0)) continue;
			dv = ep_i[i].v - ep_j[j].v;
         xdotv = dx * dv;
         dWi = dWdr_SPH(rij,ep_i[i].h);
         result[i].divv -= (ep_j[j].m * dWi) * xdotv;
		}
      result[i].divv /= ep_i[i].rho;
	}
}

//=======================================
//* Class Defs.: Calc_rotv
//=======================================
void Calc_rotv::operator ()(const SPH_EP* const ep_i,
                            const PS::S32 Nip,
                            const SPH_EP* const ep_j,
                            const PS::S32 Njp,
                            SPH_Gather_Results* const result){
   using namespace kernel_functions;
   //* Local variables
   PS::F64vec dx,dv;
   PS::F64 rij,xdotv,dWi;

   for (PS::S32 i=0; i<Nip; i++) {
      for (PS::S32 j=0; j<Njp; j++) {
         dx = ep_i[i].x - ep_j[j].x;
  			rij = std::sqrt(dx * dx);
         if ((ep_i[i].id == ep_j[j].id) && (rij == 0.0e0)) continue;
			dv = ep_i[i].v - ep_j[j].v;
         xdotv = dx * dv;
         dWi = dWdr_SPH(rij,ep_i[i].h);
         result[i].rotv += (ep_j[j].m * dWi) * (dv ^ dx); 
		}
      result[i].rotv /= ep_i[i].rho;
	}
}

//=====================================
//* Class Defs.: Calc_vsig_S05
//=====================================
void Calc_vsig_S05::operator ()(const SPH_EP* const ep_i,
                                const PS::S32 Nip,
                                const SPH_EP* const ep_j,
                                const PS::S32 Njp,
                                SPH_Symmetry_Results* const result){
   //* Local variables
   PS::F64vec dx,dv;
   PS::F64 rij,xdotv,wij;

   for (PS::S32 i=0; i<Nip; i++) {
      for (PS::S32 j=0; j<Njp; j++) {
         dx = ep_i[i].x - ep_j[j].x;
         rij = std::sqrt(dx * dx);
         if ((ep_i[i].id == ep_j[j].id) && (rij == 0.0e0)) continue;
         dv = ep_i[i].v - ep_j[j].v;
         xdotv = dx * dv;
         wij = 0.0e0;
         if (xdotv < 0.0e0) wij=xdotv/rij;
         result[i].vsig = std::max(result[i].vsig,
                                  ep_i[i].cs + ep_j[j].cs - 3.0e0*wij);
      }
   }

}

//=======================================
//* Class Defs.: Calc_vsig_CD10
//=======================================
void Calc_vsig_CD10::operator ()(const SPH_EP* const ep_i,
                                 const PS::S32 Nip,
                                 const SPH_EP* const ep_j,
                                 const PS::S32 Njp,
                                 SPH_Gather_Results* const result){
   //* Local variables
   PS::F64vec dx,dv;
   PS::F64 rij,xdotv,cs_ij;

   for (PS::S32 i=0; i<Nip; i++) {
      for (PS::S32 j=0; j<Njp; j++) {
         dx = ep_i[i].x - ep_j[j].x;
  			rij = std::sqrt(dx * dx);
         if ((ep_i[i].id == ep_j[j].id) && (rij == 0.0e0)) continue;
			dv = ep_i[i].v - ep_j[j].v;
         xdotv = dx * dv;
         if (xdotv < 0.0e0) {
            xdotv = xdotv/rij;
         } else {
            xdotv = 0.0e0;
         }
         cs_ij = 0.5e0 * (ep_i[i].cs + ep_j[j].cs);
         result[i].vsig = std::max(result[i].vsig,cs_ij-xdotv);
		}
	}
}

//=======================================
//* Class Defs.: Calc_RoC_velcfld
//=======================================
void Calc_RoC_velcfld::operator ()(const SPH_EP* const ep_i,
                                   const PS::S32 Nip,
                                   const SPH_EP* const ep_j,
                                   const PS::S32 Njp,
                                   SPH_Gather_Results* const result){
   using namespace Eigen;
   using namespace kernel_functions;
   //* Local variables
   PS::F64vec dx,dv,da;
   PS::F64 rij,xdotv,Wi,dWi,p,wj;
#if (Simulation_Dimension == 1)
   PS::F64 T,Tinv,D,E;
#elif ((Simulation_Dimension == 2) || (Simulation_Dimension == 3))
   Matrix<PS::F64,dimSpace,dimSpace> T,Tinv,D,E,V,V2,Vtrns,A,S,W,Wtrns;
#endif

   //* Compute divv, Dtdivv, trSS, trWW
   for (PS::S32 i=0; i<Nip; i++) {
#if (Simulation_Dimension == 1)
      T = 0.0e0;
      D = 0.0e0;
      E = 0.0e0;
#elif ((Simulation_Dimension == 2) || (Simulation_Dimension == 3))
      T = T.Zero();
      D = D.Zero();
      E = E.Zero();
#endif
      for (PS::S32 j=0; j<Njp; j++) {
         dx = ep_i[i].x - ep_j[j].x;
         rij = std::sqrt(dx * dx);
         if ((ep_i[i].id == ep_j[j].id) && (rij == 0.0e0)) continue;
			dv = ep_i[i].v    - ep_j[j].v;
         da = ep_i[i].atot - ep_j[j].atot;

         //* Weights
         //-[Common Part]
         p   = dimSpace + 2.0e0;
         //-[Choice 1]
         //Wi  = W_SPH(rij,ep_i[i].h);
         //wj  = ep_j[j].m * Wi / ep_j[j].rho;
         //wj  = (ep_j[j].gameff-1.0e0) * ep_j[j].m * ep_j[j].u 
         //    * Wi / ep_j[j].pres;
         //-[Choice 2]
         dWi = dWdr_SPH(rij,ep_i[i].h);
         wj  = ep_j[j].m * pow(ep_i[i].h,p) * dWi / ep_j[j].rho;
         //wj  = (ep_j[j].gameff-1.0e0) * ep_j[j].m * ep_j[j].u
         //    * pow(ep_i[i].h,p) * dWi / ep_j[j].pres;

#if (Simulation_Dimension == 1)
         T += wj * dx.x * dx.x;
         D += wj * dv.x * dx.x;
         E += wj * da.x * dx.x;
      }
      Tinv = 1.0e0/T;
      result[i].divv   = D * Tinv;
      result[i].Dtdivv = E * Tinv - result[i].divv * result[i].divv;
      result[i].trSS   = 0.0e0;
      result[i].trWW   = 0.0e0;
      // Because shear and vortex do not exit in 1D simulations.
#else
#if (Simulation_Dimension == 2)
        T(0,0) += wj * dx.x * dx.x;
        T(1,1) += wj * dx.y * dx.y;
        T(0,1) += wj * dx.x * dx.y;

        D(0,0) += wj * dv.x * dx.x;
        D(1,1) += wj * dv.y * dx.y;
        D(0,1) += wj * dv.x * dx.y;

        E(0,0) += wj * da.x * dx.x;
        E(1,1) += wj * da.y * dx.y;
        E(0,1) += wj * da.x * dx.y;
     }
     T(1,0) = T(0,1);
     D(1,0) = D(0,1);
     E(1,0) = E(0,1);
#elif (Simulation_Dimension == 3)
         T(0,0) += wj * dx.x * dx.x;
         T(1,1) += wj * dx.y * dx.y;
         T(2,2) += wj * dx.z * dx.z;
         T(0,1) += wj * dx.x * dx.y;
         T(0,2) += wj * dx.x * dx.z;
         T(1,2) += wj * dx.y * dx.z;

         D(0,0) += wj * dv.x * dx.x;
         D(1,1) += wj * dv.y * dx.y;
         D(2,2) += wj * dv.z * dx.z;
         D(0,1) += wj * dv.x * dx.y;
         D(0,2) += wj * dv.x * dx.z;
         D(1,2) += wj * dv.y * dx.z;

         E(0,0) += wj * da.x * dx.x;
         E(1,1) += wj * da.y * dx.y;
         E(2,2) += wj * da.z * dx.z;
         E(0,1) += wj * da.x * dx.y;
         E(0,2) += wj * da.x * dx.z;
         E(1,2) += wj * da.y * dx.z;
      }
      T(1,0) = T(0,1);
      T(2,0) = T(0,2);
      T(2,1) = T(1,2);
      D(1,0) = D(0,1);
      D(2,0) = D(0,2);
      D(2,1) = D(1,2);
      E(1,0) = E(0,1);
      E(2,0) = E(0,2);
      E(2,1) = E(1,2);
#endif
      Tinv = T.inverse();
      V = D * Tinv;
      A = E * Tinv;
      result[i].divv = V.trace();
      V2 = V * V;
      T = A - V2;
      result[i].Dtdivv = T.trace();
      Vtrns = V.transpose();
      S = 0.5e0*(V + Vtrns);
      for (PS::S32 k=0; k<dimSpace; k++) 
         S(k,k) = S(k,k) - result[i].divv/dimSpace;
      T = S * S;
      result[i].trSS = T.trace();
      W = 0.5e0*(V - Vtrns);
      Wtrns = W.transpose();
      T = W * Wtrns;
      result[i].trWW = T.trace();
#endif
   }

}

//=======================================
//* Class Defs.: Calc_alpha_loc
//=======================================
void Calc_alpha_loc::operator ()(const SPH_EP* const ep_i,
                                 const PS::S32 Nip,
                                 const SPH_EP* const ep_j,
                                 const PS::S32 Njp,
                                 SPH_Gather_Results* const result){
   using namespace math_functions;
   using namespace kernel_functions;
   using namespace calculation_conditions;
   //* Local variables
   PS::F64vec dx;
   PS::F64 rij,h,h2,Ri,s,xi,Ai;

   for (PS::S32 i=0; i<Nip; i++) {
      Ri = 0.0e0;
      h  = ep_i[i].h;
      h2 = h*h;
      for (PS::S32 j=0; j<Njp; j++) {
         dx = ep_i[i].x - ep_j[j].x;
         rij = std::sqrt(dx * dx);
         Ri += sign(ep_j[j].divv) * ep_j[j].m * W_SPH(rij,h);
      }
      Ri /= ep_i[i].rho;
      s = 2.0e0 * pow(1.0e0-Ri,4.0e0) * ep_i[i].divv;
      s = s*s;
      if (ep_i[i].trSS > 0.0e0) {
         xi = s/(s+ep_i[i].trSS);
      } else {
         if (ep_i[i].divv < 0.0e0) {
            xi = 1.0e0;
         } else {
            xi = 0.0e0;
         }
      }
      Ai = xi * std::max(-ep_i[i].Dtdivv, 0.0e0);
      result[i].alpha_loc = alpha_SPH_max * h2 * Ai
                         / (ep_i[i].vsig * ep_i[i].vsig + h2 * Ai);
   }

}

//=======================================
//* Class Defs.: Calc_fAV
//=======================================
void Calc_fAV::operator ()(const SPH_EP* const ep_i,
                           const PS::S32 Nip,
                           const SPH_EP* const ep_j,
                           const PS::S32 Njp,
                           SPH_Gather_Results* const result){
   using namespace kernel_functions;
	//* Local variables
	PS::F64vec dx,dv,rotv;
	PS::F64 rij,xdotv,dWi,divv,rotv_abs;
   PS::F64 coef;

	for (PS::S32 i=0; i<Nip; i++) {
      divv = 0.0e0;
      rotv = 0.0e0;
		for (PS::S32 j=0; j<Njp; j++) {
         dx = ep_i[i].x - ep_j[j].x;
  			rij = std::sqrt(dx * dx);
         if ((ep_i[i].id == ep_j[j].id) && (rij == 0.0e0)) continue; 
			dv = ep_i[i].v - ep_j[j].v;
         xdotv = dx * dv;
         dWi  = dWdr_SPH(rij,ep_i[i].h);
         coef = ep_j[j].m * dWi;
         divv -= coef * xdotv;
			rotv += coef * (dv ^ dx);
	   }
      rotv_abs = std::sqrt(rotv * rotv)/ep_i[i].rho;
      divv = std::abs(divv)/ep_i[i].rho;
      result[i].fAV = divv/(divv + rotv_abs + 1.0e-4 * ep_i[i].cs/ep_i[i].h);
	}

}

//================================
//* Class Defs.: Calc_gradh_term
//================================
void Calc_gradh_term::operator ()(const SPH_EP* const ep_i,
                                  const PS::S32 Nip,
                                  const SPH_EP* const ep_j,
                                  const PS::S32 Njp,
                                  SPH_Gather_Results* const result){
   using namespace kernel_functions;
   //* Local variables
   PS::F64 hi,rij;
   PS::F64vec dx;

   for (PS::S32 i=0; i<Nip; i++) {
      hi = ep_i[i].h;
      result[i].f_grdh = ep_i[i].m * dWdh_SPH(0.0e0,hi);
      for (PS::S32 j=0; j<Njp; j++) {
         dx = ep_i[i].x - ep_j[j].x;
         rij = std::sqrt(dx * dx);
         result[i].f_grdh += ep_j[j].m * dWdh_SPH(rij,hi);
       }
       result[i].f_grdh = 1.0e0/(1.0e0+(hi * result[i].f_grdh)/(dimSpace * ep_i[i].rho));
   }

}

//==================================
//* Class Defs.: Calc_hydro_force
//==================================
void Calc_hydro_force::operator ()(const SPH_EP* const ep_i,
                                   const PS::S32 Nip,
                                   const SPH_EP* const ep_j,
                                   const PS::S32 Njp,
                                   SPH_Symmetry_Results* const result){
   using namespace kernel_functions;
   //* Local variables
   PS::F64vec dx,dv;
   PS::F64 rij,xdotv;
   PS::F64 dWi,dWj,dWij;
   PS::F64 povrho2_i,povrho2_j;
   PS::F64 hij,muij,alpha_ij,fAV_ij;
   PS::F64 ap;

   //* Pressure-gradient acceleration + viscous acceleration
   for (PS::S32 i=0; i<Nip; i++) {
      povrho2_i = ep_i[i].pres/(ep_i[i].rho * ep_i[i].rho);
      for (PS::S32 j=0; j<Njp; j++) {
         dx = ep_i[i].x - ep_j[j].x;
         rij = std::sqrt(dx * dx);
         if ((ep_i[i].id == ep_j[j].id) && (rij == 0.0e0)) continue; 
         dv = ep_i[i].v - ep_j[j].v;
         xdotv = dx * dv;
         dWi  = dWdr_SPH(rij,ep_i[i].h);
         dWj  = dWdr_SPH(rij,ep_j[j].h);
         dWij = 0.5e0*(dWi+dWj);

         //--------------------------------------------------------
         //* Pure-Hydro part
         povrho2_j = ep_j[j].pres/(ep_j[j].rho * ep_j[j].rho);
         ap = ep_j[j].m * (ep_i[i].f_grdh * povrho2_i*dWi 
                          +ep_j[j].f_grdh * povrho2_j*dWj);
         result[i].ap -= (ap * dx);
         result[i].gamad += ep_j[j].m * xdotv * dWij;
         //--------------------------------------------------------

         //--------------------------------------------------------
         //* Artificial Viscosity Part
         if (xdotv < 0.0e0) {
            alpha_ij = 0.5e0 * (ep_i[i].alpha + ep_j[j].alpha);
            fAV_ij   = 0.5e0 * (ep_i[i].fAV + ep_j[j].fAV);
            //* [1] Monaghan & Gingold (1993)
            //hij = 0.5e0 * (ep_i[i].hi + ep_j[j].h);
            //muij = hij * xdotv/(rij*rij + 0.01e0*hij*hij)
            //beta_ij = 4.0e0*alpha_ij // LHS corresponds to 2\beta.
            //muij = (-alpha_ij*(ep_i[i].cs+ep_j[j].cs) 
            //        +beta_ij*muij)*muij/(ep_i[i].rho+ep_j[j].rho);
            //muij = ep_j[j].m*muij*dWij;
            //* [2] Monaghan (1997)
            muij = xdotv/rij; // LHS corresponds to w_{ij} in Springel+05.
            muij = -alpha_ij * ((ep_i[i].cs + ep_j[j].cs)-3.0e0*muij)*muij 
                 / (ep_i[i].rho + ep_j[j].rho);
            muij = ep_j[j].m * muij * dWij;
            result[i].avis   -= muij * dx * fAV_ij;
            result[i].gamvis += muij * xdotv * fAV_ij;
         }
         //--------------------------------------------------------
      }
      result[i].gamad  *= povrho2_i * ep_i[i].f_grdh;
      result[i].gamvis *= 0.5e0;
   }
}


//=====================================
//* Class Defs: Calc_kernel_length 
//=====================================
void Calc_kernel_length::operator () (const SPH_EP* const ep_i,
                                      const PS::S32 Nip,
                                      const SPH_EP* const ep_j,
                                      const PS::S32 Njp,
                                      SPH_Gather_Results* const result) {
   using namespace calculation_conditions;
   //* Local variables
   PS::F64vec dx;
   PS::F64 h,h2,hmax,h_L,h_U,dh,dh_prev;
   PS::F64 rij2;
   PS::S32 nnb,numNoChange;
   bool normalTermination,forceQuit;

   for (PS::S32 i=0; i<Nip; i++) {
      //* Load ith particle info.
      h       = ep_i[i].h_prev;
      hmax    = ep_i[i].h;
      h2      = h*h;
      h_L     = 0.0e0;
      h_U     = hmax;
      dh_prev = h_U-h_L;
      //* Determine gather neighbors and smoothing length
      //  simulataneously.
      numNoChange=0;
      for(;;) {
         //* Count # of neighbors
         nnb = 0;
         for(PS::S32 j=0; j<Njp; j++) {
            dx = ep_j[j].x - ep_i[i].x;
            rij2 = dx * dx;
            if ((rij2 < h2) && (ep_i[i].id != ep_j[j].id)) nnb = nnb + 1;
         }
         //* Check termination conditions
         normalTermination = (numPtclNeighbLB <= nnb) &&
                             (nnb <= numPtclNeighbUB);
         forceQuit = ((h == hmax) && (nnb < numPtclNeighbLB)) ||
                      (numNoChange == 4);
         if ((normalTermination == true) || (forceQuit == true)) {
            result[i].h = h;
            if (forceQuit == true) { 
               result[i].flag = 1;  
               //if (PS::Comm::getRank() == 0) {
               //   if (ep_i[i].id == 0) 
               //      std::cout << "nnb = " << nnb << std::endl;
               //}
            }
            break;
         }
         //* Update h_L & h_U
         if (nnb < numPtclNeighbLB) {
            if (h_L < h) h_L = h;
         }
         else if (numPtclNeighbUB < nnb) {
            if (h < h_U) h_U = h;
         }
         dh = h_U-h_L;
         if (dh == dh_prev) {
            numNoChange = numNoChange + 1;
         }
         else {
            dh_prev = dh;
            numNoChange = 0;
         }
         //* Update h
         h = get_next_h(nnb,h);
         if ((h <= h_L) || (h == h_U)) {
            // In this case, we switch to the bisection search.
            // The inclusion of '=' in the if statement is very
            // important to escape a limit cycle.
            h = 0.5e0*(h_L+h_U);
         }
         else if (h_U < h) {
            h = h_U;
         }
         h2 = h*h;
       } // END of iteration loop
   } // END of particle loop
}

/*-------------------------------------------------------------------*/
///////////////////////     F U N C T I O N     ///////////////////////
/////////////////////// < G E T _ N E X T _ H > ///////////////////////
/*-------------------------------------------------------------------*/
PS::F64 get_next_h(const PS::S32 nnb, const PS::F64 h_old) {
   using namespace calculation_conditions;
   //* Local parameters
   const PS::F64 pow=1.0e0/dimSpace;
   //* Local variables
   PS::F64 a,s,sinv,p;
   PS::F64 sm,sp,s0;

   sm= std::pow(((PS::F64)(numPtclNeighb)
               /((PS::F64)(numPtclNeighb)+2.0e0*dimSpace)),pow);
   sp= std::pow(((PS::F64)(numPtclNeighb)
               /((PS::F64)(numPtclNeighb)-2.0e0*dimSpace)),pow);
   s = std::pow(((PS::F64)(numPtclNeighb)
               /std::max((PS::F64)(nnb),1.0e0)),pow);
   //* [1] Simplest update method
   //return h_old*s;
   //* [2] Hernquist & Katz(1989)'s update method
   //ret = 0.5d0*h_old*(1.0d0+s);
   //* [3] Thacker et al.(2000)'s update method
   if (s < 1.0e0) {
      a = 0.2e0*(1.0e0+s*s);
      s0 = (1.0e0-sm)/std::sqrt(std::log(1.5e0));
   }
   else if (s >= 1.0e0) {
      sinv = 1.0e0/s;
      a = 0.2e0*(1.0e0+sinv*sinv*sinv);
      s0 = (sp-1.0e0)/std::sqrt(std::log(1.5e0));
   }
   p = 0.5e0*std::exp(-((s-1.0e0)*(s-1.0e0)/(s0*s0)));
   //return h_old*(1.0e0-a+a*s);
   return h_old*((1.0e0-p)*(1.0e0-a+a*s)+p);
   // see Eq.(10),(11) in Thacker et al.(2000)[MNRAS,319,619]

}
 
