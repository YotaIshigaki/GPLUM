/* Standard headers */
#include <iostream>
#include <cmath>
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "preprocess_keywords.h"
#include "Calculation_Conditions.h"

namespace calculation_conditions {

//******** Numerial Resolution ***************************************
//const unsigned long numPtclPlaned=pow(2,18);

//******** Gravity parameters ****************************************
//const double ThetaCritGrv=0.5e0;

//******** SPH parameters ********************************************
//   numPtclNeighb   : number of neighbor particles
//   numPtclNeighbLB : the lower bound of the number of neighbor particles
//   numPtclNeighhUB : the upper bound of the number of neighbor particles
//   extFacNBSearch  : an extention factor
//   alpha_SPH_min   : the minimum of \alpha parameter of A.V.
//   alpha_SPH_max   : the maximum of \alpha parameter of A.V.
//   CFL_SPH         : Courant-Friedrichs-Lewy (CFL) number for 
//                     hydrodynamic part.
//
//*********************************************************************
//const int numPtclNeighb = 50;
const int numPtclNeighb = 96;
const int numPtclNeighbLB = numPtclNeighb - 1;
const int numPtclNeighbUB = numPtclNeighb + 1;
const double extFacNBSearch = 1.25e0;
const double alpha_SPH_min = 0.1e0;
const double alpha_SPH_max = 1.0e0; 
const double CFL_SPH = 0.1e0;

//******** IO parameters **********************************************
//   nstepMax      : maximum number of steps
//   tstop_CGS     : the termination time of the simulation (in CGS)
//   tdump_ini_CGS : the initial time of output (in CGS)
//   dtdump_CGS    : the time interval of output (in CGS)
//   tdump         : the time of output
//*********************************************************************
//const int nstepMax=1000000;
//const int nstepMax=100000;
//const int nstepMax=10000;
const int nstepMax=3000;
//const int nstepMax=100;
//const int nstepMax=10;
//const int nstepMax=1;

const double tstop_CGS=0.25e0;
const double dtdump_CGS=0.25e0;
const double tdump_ini_CGS=dtdump_CGS;


} // END of calculation_conditions
