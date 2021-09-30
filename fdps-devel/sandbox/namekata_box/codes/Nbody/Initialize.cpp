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
#include "Nbody_Objects_Class.h"
#include "Simulation_Units.h"
#include "Initial_Conditions.h"
#include "IO.h"
#include "Initialize.h"

/*-------------------------------------------------------------------*/
///////////////////////   S U B R O U T I N E   ///////////////////////
/////////////////////// < I N I T I A L I Z E > ///////////////////////
/*-------------------------------------------------------------------*/
void Initialize(int argc, char* argv[],
                Nbody_Objects& Nbody_objs) {
   using namespace simulation_units;
   //* Local variables

   //* Initialize FDPS and particle systems
   PS::Initialize(argc,argv);
   Nbody_objs.system.initialize();
   
   //* Other initializations
   IO_ctrl.initialize();
   setup_sim_units();

#if (Restart_Mode == 0)
   //* Set initial condition
   PM_test_IC(Nbody_objs.system, 
              Nbody_objs.dinfo);

	//* Make trees and initiale them
	Nbody_objs.init_tree();

	//* Compute gravity
	Nbody_objs.calc_gravity();

   //* Output I.C.
   IO_ctrl.write_data(Nbody_objs);

#elif (Restart_Mode == 1)
   //* Read restart data
   IO_ctrl.read_data(Nbody_objs,0);
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
