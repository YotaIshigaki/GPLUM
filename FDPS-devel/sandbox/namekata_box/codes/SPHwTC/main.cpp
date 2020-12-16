/* Standard headers */
#include <cstdio>
#include <iostream>
#include <cmath>
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "preprocess_keywords.h"
#include "Simulation_Units.h"
#include "Calculation_Conditions.h"
#include "SPH_Objects_Class.h"
#include "Misc_Class.h"
#include "Initialize.h"
#include "Time_Integration.h"
#include "IO.h"

/*-------------------------------------------------------------------*/
////////////////////// M A I N   F U N C T I O N //////////////////////
/*-------------------------------------------------------------------*/
int main(int argc, char* argv[]) {
   using namespace simulation_units;
   using namespace calculation_conditions;
   using namespace time_integrator;

   //* Local variables
   PS::S64 nstep;
   SPH_Objects SPH_objs;
   Fluctuation_Monitor f_monitor;

   //* Initialize various qtys.
   Initialize(argc,argv,SPH_objs,f_monitor);

   //* Time Integration (main loop)
	nstep = 0; IO_ctrl.time = 0.0e0;
	for (;;){
      // Update # of steps
      nstep++;

      // Compute global timestep
      IO_ctrl.dt = SPH_objs.calc_timestep();
      if ((nstep == 1) || (nstep % 10 == 0)) {
         std::cout << "--------------------------------" << std::endl;
         std::cout << "nstep = " << nstep << std::endl;
         std::cout << "time  = " << IO_ctrl.time << std::endl;
         std::cout << "dt    = " << IO_ctrl.dt   << std::endl;
      }
      
		//* Leapfrog: KD (1st stage)
      Leapfrog_KickDrift(SPH_objs.system,
                         SPH_objs.dinfo,
                         IO_ctrl.dt);

      //* Compute gravitational accelerations
      //calc_gravity(SPH_objs.system,
      //             SPH_objs.dinfo);

		//* Update density, pressure, etc. and compute accelerations
      SPH_objs.update_h();
      SPH_objs.calc_rho();
		SPH_objs.calc_pres();
		SPH_objs.calc_hydroforce(IO_ctrl.dt);

		//* Leap frog: K (2nd stage)
		Leapfrog_FinalKick(SPH_objs.system,
                         IO_ctrl.dt);

      //* Dump velocity and check density fluctuation
#if (Current_Code_Mode == Make_Glass_Mode)
      f_monitor.dump_velocity(SPH_objs.system,
                              IO_ctrl.dt);
      f_monitor.check(SPH_objs.system,
                      IO_ctrl.time, 
                      IO_ctrl.stop_time);
#endif

		//* Output
      IO_ctrl.write_data(SPH_objs);

      //* Update time
      IO_ctrl.time += IO_ctrl.dt;

      //* Termination conditions
      if (nstep == nstepMax) break;
      if (IO_ctrl.time > IO_ctrl.stop_time) break;
      
	}

   //* Finalize main
   Finalize();
   return 0;
}
