#include "config.h"
#if (defined(HAVE_GRAPE7) || defined(HAVE_GRAPE6A) || defined(HAVE_GRAPE5) || defined(HAVE_PHANTOM_GRAPE))
#include <gp5util.h>
#else
#include "GRAPEEmulator.h"
#define HAVE_GRAPE_EMULATION    (ON)
#endif

#define ReleaseGRAPE_EveryTime (ON)
#define LogPrint    (ON)

static bool HoldFlag = OFF;
static double HoldTime;
void HoldGRAPE(void){

    if(HoldFlag == OFF){
#ifdef HAVE_GRAPE_EMULATION
	    g5_open_emu();
#else
	    g5_open();
#endif
        HoldTime = GetElapsedTime();
        HoldFlag = ON;
#ifdef PRINT_LOG_HOLD_RELEASE_GRAPE 
#ifndef HAVE_GRAPE_EMULATION
        fprintf(stderr,"Hold GRAPE!! Host name is %s\n",MPIGetProcessorName());
        fflush(NULL);
#endif
#endif
    }

    return;
}

void ReleaseGRAPE(void){

    if(((ReleaseGRAPE_EveryTime == ON)&&(HoldFlag == ON))
            ||((HoldFlag == ON)&&(GetElapsedTime()-HoldTime > GRAPE_HOLDING_TIME_IN_SEC))){
#ifdef HAVE_GRAPE_EMULATION
        g5_close_emu();
#else
        g5_close();
#endif
        HoldFlag = OFF;
#ifdef PRINT_LOG_HOLD_RELEASE_GRAPE 
#ifndef HAVE_GRAPE_EMULATION
        fprintf(stderr,"Release GRAPE!! %g sec. Host name is %s\n",
                GetElapsedTime()-HoldTime,MPIGetProcessorName());
        fflush(NULL);
#endif
#endif
    }

    return;
}

