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
#include "Nbody_Objects_Class.h"
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
   Nbody_Objects Nbody_objs;

   //* Initialize various qtys.
   Initialize(argc,argv,Nbody_objs);

   //* Time Integration (main loop)
	nstep = 0; IO_ctrl.time = 0.0e0;
	for (;;){
      // Update # of steps
      nstep++;

      // Compute global timestep
      IO_ctrl.dt = Nbody_objs.calc_timestep();
      if ((nstep == 1) || (nstep % 10 == 0)) {
         std::cout << "--------------------------------" << std::endl;
         std::cout << "nstep = " << nstep << std::endl;
         std::cout << "time  = " << IO_ctrl.time << std::endl;
         std::cout << "dt    = " << IO_ctrl.dt   << std::endl;
      }
      
		//* Leapfrog: KD (1st stage)
      Leapfrog_KickDrift(Nbody_objs.system,
                         Nbody_objs.dinfo,
                         IO_ctrl.dt);

      //* Compute gravitational accelerations
      Nbody_objs.calc_gravity();

		//* Leap frog: K (2nd stage)
		Leapfrog_FinalKick(Nbody_objs.system,
                         IO_ctrl.dt);

		//* Output
      IO_ctrl.write_data(Nbody_objs);

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
