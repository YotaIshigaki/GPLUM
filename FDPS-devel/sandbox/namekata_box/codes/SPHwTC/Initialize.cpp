/* Standard headers */
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
/* FDPS header */
#include <particle_simulator.hpp>
/* User-defined header */
#include "preprocess_keywords.h"
#include "Particle_Class.h"
#include "SPH_Objects_Class.h"
#include "Misc_Class.h"
#include "Simulation_Units.h"
#include "Initial_Conditions.h"
#include "IO.h"
#include "Initialize.h"

/*-------------------------------------------------------------------*/
///////////////////////   S U B R O U T I N E   ///////////////////////
/////////////////////// < I N I T I A L I Z E > ///////////////////////
/*-------------------------------------------------------------------*/
void Initialize(int argc, char* argv[],
                SPH_Objects& SPH_objs,
                Fluctuation_Monitor& f_monitor) {
   using namespace simulation_units;
   //* Local variables

   //* Initialize FDPS and particle systems
   PS::Initialize(argc,argv);
   SPH_objs.system.initialize();
   
   //* Other initializations
   IO_ctrl.initialize();
   setup_sim_units();

#if (Restart_Mode == 0)
   //* Set initial condition
#if (Current_Code_Mode == Make_Glass_Mode)
   PM_test_IC(SPH_objs.system, 
              SPH_objs.dinfo);
   //make_glass_IC(SPH_objs.system,
   //              SPH_objs.dinfo,
   //              f_monitor);
#elif (Current_Code_Mode == Shock_Tube_Test_Mode)
   shock_tube_IC(SPH_objs.system,
                 SPH_objs.dinfo);
#else
#error Invalid setting for `Current_Code_Mode`.
#endif

	//* Make trees and initiale them
	SPH_objs.init_tree();

	//* Compute density/pressure/accelerations
   SPH_objs.update_h();
   SPH_objs.calc_rho();
   SPH_objs.calc_pres();
	SPH_objs.calc_hydroforce(0.0e0);

   //* Output I.C.
   IO_ctrl.write_data(SPH_objs);

#elif (Restart_Mode == 1)
   //* Read restart data
   IO_ctrl.read_data(SPH_objs,0);
#else
#error Invalid setting for `Restart_Mode`.
#endif

   //* Output STDOUT
   if (PS::Comm::getRank() == 0) 
      std::cout << "Initialize() completed!" << std::endl;
}

/*-------------------------------------------------------------------*/
////////////////////////  S U B R O U T I N E  ////////////////////////
////////////////////////  < F I N A L I Z E >  ////////////////////////
/*-------------------------------------------------------------------*/
void Finalize() {

   //* Finalize FDPS
   PS::Finalize();

}
