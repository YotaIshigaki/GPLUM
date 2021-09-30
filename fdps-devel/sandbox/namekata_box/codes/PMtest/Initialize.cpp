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
#include "Misc_Class.h"
#include "Initial_Conditions.h"
#include "Initialize.h"

/*-------------------------------------------------------------------*/
///////////////////////   S U B R O U T I N E   ///////////////////////
/////////////////////// < I N I T I A L I Z E > ///////////////////////
/*-------------------------------------------------------------------*/
void Initialize(Nbody_Objects& Nbody_objs, 
                Crystal_Parameters& NaCl_params,
                PS::S32 nstep) {
   //* Local variables
   static bool first_call=true;

   //* Check
   if (PS::Comm::getRank() == 0)
      std::cout << "first_call = " << first_call << std::endl;

   //* Initialize FDPS objects
   if (first_call == true) {
      Nbody_objs.system.initialize();
      Nbody_objs.dinfo.initialize();
   }
   // [Note #1]
   //    This is because system.initialize() is not allowed to
   //    be called twice. 
   
   //* Set initial condition
#if (NaCl_Crystal_Mode == 1)
   NaCl_IC(Nbody_objs.system, 
           Nbody_objs.dinfo,
           NaCl_params);
#elif (NaCl_Crystal_Mode == 0)
   PM_test_IC(Nbody_objs.system, 
              Nbody_objs.dinfo,
              NaCl_params);
#else
#error Invalid value of NaCl_Crystal_Mode.
#endif

	//* Make trees and initiale them
   if (first_call == true) {
	   Nbody_objs.init_tree();
   }
   // [Note #2]
   //    The same reason as Note #1.

	//* Compute gravity
	Nbody_objs.calc_gravity();

   //* Update first_call
   first_call = false;

}

/*-------------------------------------------------------------------*/
////////////////////////  S U B R O U T I N E  ////////////////////////
////////////////////////  < F I N A L I Z E >  ////////////////////////
/*-------------------------------------------------------------------*/
void Finalize() {

   //* Finalize FDPS
   PS::Finalize();

}
