#include "config.h"
#include "SetUpTestRun.h"

static bool first = true;
// #ifdef TASK_NFW
// extern double M200_NFW;
// #endif
void ForceFromExternalPotentials(void){

#ifdef TASK_KEPLER //{
    KeplerPotential();
#endif // TASK_KEPLER //}

#ifndef GRAVITY_RUN
    return ;
#else  //GRAVITY_RUN

    if(first == true){ // Call initializations of external potentials if necessary.
#ifdef TASK_AGNTORUS
        InitializeAGNTorusExternalPotentials();
#endif //TASK_AGNTORUS
        first = false;
    }
    
#ifdef TASK_MW
    MilkyWayPotentialHaloDisk();
#endif //TASK_MW
#ifdef TASK_NFW
    ExternalPotentialForNFW();
#endif //TASK_MW
#ifdef TASK_AGNTORUS
    AGNForceFromExternalPotentials(); // set init function
#endif //TASK_AGNTORUS

    return ;
#endif //GRAVITY_RUN

}
