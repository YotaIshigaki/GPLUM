#pragma once

namespace calculation_conditions {

//******** Numerial Resolution ***************************************

//******** Gravity parameters ****************************************


//******** SPH parameters ********************************************
extern const int numPtclNeighb;
extern const int numPtclNeighbLB;
extern const int numPtclNeighbUB;
extern const double extFacNBSearch;
extern const double alpha_SPH_min;
extern const double alpha_SPH_max;
extern const double CFL_SPH;

//******** IO parameters **********************************************
extern const int nstepMax;
extern const double tstop_CGS;
extern const double dtdump_CGS;
extern const double tdump_ini_CGS; 

} // END of calculation_conditions
