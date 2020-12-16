/* Standard headers */
#include <cmath>
#include <limits>
#include <vector>
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "preprocess_keywords.h"
#include "Particle_Class.h"
#include "SPH_Objects_Class.h"
#include "Calculation_Conditions.h"

//================================
//* Class Defs.: SPH_Objects
//================================
void SPH_Objects::init_tree() {
   //PS::S32 numPtclLocal = system.getNumberOfParticleGlobal();
   //PS::U64 n_glb_tot = 3 * numPtclLocal;
   PS::S32 numPtclLocal = system.getNumberOfParticleLocal();
   PS::U64 n_glb_tot = 3 * numPtclLocal;
   if (PS::Comm::getRank() == 0)
      std::cout << "n_glb_tot = " << n_glb_tot << std::endl;
   gat_tree.initialize(n_glb_tot);
   sym_tree.initialize(n_glb_tot);
}
void SPH_Objects::update_h() {
   using namespace calculation_conditions;
   //* Local variables
   PS::S32 iter=0,flagSumLocal,flagSum;
   PS::S32 numPtclLocal = system.getNumberOfParticleLocal();
   std::vector<PS::F64> extFac;
   SPH_Gather_Results result;

   //* Save original kernel length
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      system[i].h_prev = system[i].h;
      extFac.push_back(extFacNBSearch);
   }
   //* Update h so that # of neighbors is close to numPticlNeighb
   for(;;) {
      iter++;
      //if (PS::Comm::getRank() == 0) 
      //   std::cout << "iter = " << iter << std::endl;
      //* Set an extended kernel length
      for (PS::S32 i=0; i<numPtclLocal; i++)
         system[i].h = extFac[i] * system[i].h_prev;
      //* Update h
      gat_tree.calcForceAll(Calc_kernel_length(), system, dinfo);
      for (PS::S32 i=0; i<numPtclLocal; i++) {
         result = gat_tree.getForce(i);
         system[i].flag = result.flag;
         system[i].h    = result.h;
      }
      //* Check error
      flagSumLocal=0; flagSum=0;
      for (PS::S32 i=0; i<numPtclLocal; i++) flagSumLocal += system[i].flag;
      flagSum = PS::Comm::getSum(flagSumLocal);
      if (flagSum == 0) {
         // In this case, there is no error.
         break;
      } else {
         // In this case, we repeat neighbor search with a larger extFac. 
         if (flagSumLocal > 0) {
            for (PS::S32 i=0; i<numPtclLocal; i++)
               if (system[i].flag == 1)
                  extFac[i] *= 1.1e0; 
         }
         //if (PS::Comm::getRank() == 0) 
         //   std::cout << "flagSum = " << flagSum << std::endl;
      }
   }
}
void SPH_Objects::calc_rho() {
	gat_tree.calcForceAll(Calc_density(), system, dinfo);
   for (PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++) {
      SPH_Gather_Results result = gat_tree.getForce(i);
      system[i].rho = result.rho;
   }

}
void SPH_Objects::calc_pres() {
   for (PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++){
#if (Treatment_EnergyEq != Isothermal)
      system[i].Tgas = system[i].mu * (system[i].gameff - 1.0e0) * system[i].u;
#endif
      system[i].pres = system[i].rho * system[i].Tgas / system[i].mu;
      system[i].cs   = std::sqrt(system[i].gameff * system[i].Tgas / system[i].mu);
   }
}
void SPH_Objects::calc_hydroforce(PS::F64 dt) {
   using namespace calculation_conditions;
   //* Local variables
   PS::S32 numPtclLocal = system.getNumberOfParticleLocal();

   //* [1] Compute Balsara switch
#if (Balsara_Switch == 1)
   gat_tree.calcForceAll(Calc_fAV(), system, dinfo);
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      SPH_Gather_Results result = gat_tree.getForce(i);
      system[i].fAV = result.fAV;
   }
#endif

   //* [2] Update \alpha
#if (Morris_Monaghan_Switch == 1)
   //* Compute divv
   gat_tree.calcForceAll(Calc_divv(), system, dinfo);
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      SPH_Gather_Results result = gat_tree.getForce(i);
      system[i].divv = result.divv;

		//* Update \alpha
      PS::F64 tau = system[i].h / (0.2e0 * system[i].cs);
      PS::F64 s = std::max(-system[i].divv, 0.0e0);
      PS::F64 dadt = - (system[i].alpha - alpha_SPH_min)/tau + s;
      system[i].alpha += dadt * dt;

      //* Adjust alpha so that alpha_SPH_min < alpha < alpha_SPH_max
      if (system[i].alpha < alpha_SPH_min) 
         system[i].alpha = alpha_SPH_min;
		if (system[i].alpha > alpha_SPH_max)
         system[i].alpha = alpha_SPH_max; 
   }
#elif (Cullen_Dehnen_Switch == 1)
   //* Compute v_{sig} and atot
   gat_tree.calcForceAll(Calc_vsig_CD10(), system, dinfo);
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      SPH_Gather_Results result = gat_tree.getForce(i);
      system[i].vsig = result.vsig;
      system[i].atot = system[i].ap 
                     + system[i].avis;
   }

   //* Compute divv,Dtdivv,trSS,trWW
   gat_tree.calcForceAll(Calc_RoC_velcfld(), system, dinfo);
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      SPH_Gather_Results result = gat_tree.getForce(i);
      system[i].divv   = result.divv;
      system[i].Dtdivv = result.Dtdivv;
      system[i].trSS   = result.trSS;
      system[i].trWW   = result.trWW;
   }

   //* Compute \alpha_{loc}
   gat_tree.calcForceAll(Calc_alpha_loc(), system, dinfo);
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      SPH_Gather_Results result = gat_tree.getForce(i);
      system[i].alpha_loc = result.alpha_loc;
   }

   //* Update \alpha
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      if (system[i].alpha < system[i].alpha_loc) {
         system[i].alpha = system[i].alpha_loc;
      } else {
         tau = system[i].h / (2.0e0 * 0.05e0 * system[i].vsig);
         system[i].alpha = system[i].alpha_loc
                         + (system[i].alpha 
                           -system[i].alpha_loc)*exp(-dt/tau);
      }
   }
#endif

   //* [3] Compute grad-h term
   gat_tree.calcForceAll(Calc_gradh_term(), system, dinfo);
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      SPH_Gather_Results result = gat_tree.getForce(i);
      system[i].f_grdh = result.f_grdh;
   }

   //* [4] Compute accelerations
   try {
      sym_tree.calcForceAll(Calc_hydro_force(), system, dinfo);
   } catch (std::bad_alloc& e) {
      if (PS::Comm::getRank() == 0)
         std::cout << e.what() << " in sym_tree.calcForceAll()." << std::endl;
         PS::Finalize();
         std::exit(1);
   }
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      SPH_Symmetry_Results result = sym_tree.getForce(i);
      system[i].ap     = result.ap;
      system[i].avis   = result.avis;
      system[i].gamad  = result.gamad;
      system[i].gamvis = result.gamvis;
   }
}
PS::F64 SPH_Objects::calc_timestep() {
   using namespace calculation_conditions;
   //* Local constants
   const PS::F64 eps=1.0e-100;

   //* Compute signal velocity
   sym_tree.calcForceAll(Calc_vsig_S05(), system, dinfo);
   //* Compute local minimum
   PS::F64 dt = std::numeric_limits<double>::max();
   for (PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++) {
      SPH_Symmetry_Results result = sym_tree.getForce(i);
      system[i].vsig = result.vsig;

      //* [2] Absolute velocity and acceleration
      PS::F64 v = std::sqrt(system[i].v * system[i].v) + system[i].cs;
      PS::F64vec atot = system[i].ap + system[i].avis;
      PS::F64 atot_abs = std::sqrt(atot * atot);

      //* [3] Find local minimum timescale
      PS::F64 dt_ith = std::min(system[i].h/system[i].vsig,
                                v/(atot_abs + eps));
      dt = std::min(dt,dt_ith);
   }
   //* Compute global minium
	return CFL_SPH * PS::Comm::getMinValue(dt);

}

