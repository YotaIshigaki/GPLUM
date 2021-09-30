/* Standard headers */
#include <cmath>
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "Particle_Class.h"
#include "Calculation_Conditions.h"
#include "Misc_Class.h"

//======================================
//* Class Defs.: Fluctuation_Monitor
//======================================
Fluctuation_Monitor::Fluctuation_Monitor() {
   check_time     = 0.0e0;
   check_interval = 0.0e0;
}
void Fluctuation_Monitor::set_config(PS::F64 start_time,
                                     PS::F64 interval) {
   check_time     = start_time;
   check_interval = interval;
}
void Fluctuation_Monitor::dump_velocity(PS::ParticleSystem<SPH_FP>& system,
                                        PS::F64 dt) {
   using namespace calculation_conditions;

   for (PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++) {
      PS::F64 tau = system[i].h / system[i].cs;
      PS::F64 facDump = std::exp(-0.1e0 * (CFL_SPH/0.1e0) * (dt/tau));
      system[i].v *= facDump;
   }

}
void Fluctuation_Monitor::check(PS::ParticleSystem<SPH_FP>& system,
                                PS::F64 time, PS::F64 stop_time){

   //* Error handling
   if ((check_time <= 0.0e0) || (check_interval <= 0.0e0)) {
      std::cout << "Error detected in check_fluctuation()." << std::endl;
      std::exit(1);
   }

   if (time > check_time) {
      //* Local variables
      PS::S32 numPtclLocal = system.getNumberOfParticleLocal();
      PS::S32 numPtcl = system.getNumberOfParticleGlobal();
      //* Compute average density
      PS::F64 rho_avrg = 0.0e0;
      for (PS::S32 i=0; i<numPtclLocal; i++) rho_avrg += system[i].rho;
      rho_avrg = PS::Comm::getSum(rho_avrg);
      rho_avrg /= numPtcl; 
      //* Compute density dispersion
      PS::F64 rho_disp = 0.0e0;
      for (PS::S32 i=0; i<numPtclLocal; i++) {
         rho_disp += std::pow(system[i].rho - rho_avrg,2); 
      }
      rho_disp = PS::Comm::getSum(rho_disp);
      rho_disp /= numPtcl;
      rho_disp = std::sqrt(rho_disp);
      //* Output the result
      if (PS::Comm::getRank() == 0) {
         std::cout << "----- Density fluctuation Monitor -----" << std::endl;
         std::cout << "   time/stop_time = " << time/stop_time << std::endl;
         std::cout << "   rho_avrg       = " << rho_avrg << std::endl;
         std::cout << "   rho_disp       = " << rho_disp << std::endl;
         std::cout << "   disp/avrg      = " << rho_disp/rho_avrg << std::endl;
         std::cout << "---------------------------------------" << std::endl;
      }
      //* Set next check time
      check_time += check_interval;
   }

}
