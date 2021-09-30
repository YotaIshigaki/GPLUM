#include "config.h"
#include "HydroMisc.h"
#include "IMFParameters.h"
#include "NeighborSearch.h"

#ifdef TASK_AGNTORUS
static void ExtraOperationsForHydroInAGNRun(void);
#endif // TASK_AGNTORUS

void CalcExtraOperationsForHydro(void){

#ifdef TASK_AGNTORUS
    ExtraOperationsForHydroInAGNRun();
#endif // TASK_AGNTORUS

    return ;
}

#ifdef TASK_AGNTORUS
#define DURATOIN_TIME_AGNTORUS (3.e+6) //in yr.
#define SF_RATE_AGNTORUS   (0.1) //in  Msun/yr.
static bool first_agn = true;
static double TotalEventNumber; 
static double ReferenceMass,InvReferenceMass;
static double InputEnergyInAnEvent;
static int EventCounter = 0;
static void InitExtraOperationsForHydroInAGNRun(void){

    // Get specific SN rate.
    double SpecificSNrate = 0.0074;
#if (UseSFModelSpawn)
    ReferenceMass = Phydro[0]->SpawnMass;
#else
    ReferenceMass = Phydro[0]->Mass;
#endif
    InvReferenceMass = 1.0/ReferenceMass;
    //InputEnergyInAnEvent = SNIIEnergyPerNumber*(SpecificSNrate*ReferenceMass)/GetUnitEnergy();
    InputEnergyInAnEvent = SNIIEnergyPerNumber*(SpecificSNrate*(ReferenceMass*Pall.UnitMass/MSUN_CGS))
                            *(SQ(Pall.UnitTime)/(Pall.UnitMass*SQ(Pall.UnitLength))); // ergs. g (cm/s) ^2
    TotalEventNumber = (SF_RATE_AGNTORUS*MSUN_CGS/Pall.UnitMass)*DURATOIN_TIME_AGNTORUS/ReferenceMass;

    fprintf(stderr,"Reference Mass is %g [Msun]\n",ReferenceMass*Pall.UnitMass/MSUN_CGS);
    fprintf(stderr,"InputEnergyInAnEvnt is %g [in simulation unit]\n",InputEnergyInAnEvent);
    fprintf(stderr,"TotalEventNumber is %g\n",TotalEventNumber);

    return;
}

static void ExtraOperationsForHydroInAGNRun(void){

    if(Pall.TCurrent*Pall.UnitTime > DURATOIN_TIME_AGNTORUS*YEAR_CGS)
        return ;

    if(first_agn == true){
        InitExtraOperationsForHydroInAGNRun();
        first_agn = false;
    }

    double EventNumberInThisStep = 
        (SF_RATE_AGNTORUS*MSUN_CGS/Pall.UnitMass)*(Pall.dtnow*Pall.UnitTime/YEAR_CGS)*InvReferenceMass;
    int IntegerEventNumberInThisStep = (int)EventNumberInThisStep;
    double Residual = EventNumberInThisStep-(double)IntegerEventNumberInThisStep;
    fprintf(stderr,"SN in this step is %d + %g\n",IntegerEventNumberInThisStep,Residual);
    // call a dice.
    if(gsl_rng_uniform(RandomGenerator) < Residual)
        IntegerEventNumberInThisStep ++;
    fprintf(stderr,"The number of SNe in this step is %d, Actives %ld\n",
            IntegerEventNumberInThisStep,Pall.NActivesHydro_t);

    int LocalEvents = 0;
#if 1
    double Rmax = 32*PC_CGS/Pall.UnitLength;
    int Actives = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
        double R2 = SQ(PhydroBody(i)->PosP[0])+SQ(PhydroBody(i)->PosP[1]);
        if(R2 < SQ(Rmax)){
            Actives ++;
        }
        }
    }

    //double EventCriterion = (double)IntegerEventNumberInThisStep/(double)Pall.NActivesHydro_t;
    double EventCriterion = (double)IntegerEventNumberInThisStep/(double)Actives;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
        double R2 = SQ(PhydroBody(i)->PosP[0])+SQ(PhydroBody(i)->PosP[1]);
        if(R2 < SQ(Rmax)){
            if(gsl_rng_uniform(RandomGenerator) < EventCriterion){
                Phydro[i]->DQheat += InputEnergyInAnEvent;
                fprintf(stderr,"input energy[%d,%ld] = %g, %g, original U is %g, T = %g\n",
                        i,PhydroBody(i)->GlobalID,Phydro[i]->DQheat,InputEnergyInAnEvent,Phydro[i]->U*PhydroBody(i)->Mass,
                        Pall.ConvertUtoT*Phydro[i]->U);
                LocalEvents ++;
            }
        }
        }
    }
#else
    double EventCriterion = (double)IntegerEventNumberInThisStep/(double)Pall.Nhydro_t;
    for(int i=0;i<Pall.Nhydro;i++){
        if(gsl_rng_uniform(RandomGenerator) < EventCriterion){
            Phydro[i]->DQheat += InputEnergyInAnEvent;
            fprintf(stderr,"input energy[%d] = %g, %g, original U is %g\n",
                    i,Phydro[i]->DQheat,InputEnergyInAnEvent,Phydro[i]->U*PhydroBody(i)->Mass);
            LocalEvents ++;
        }
    }
#endif
    int GlobalEvents;
    MPI_Allreduce(&LocalEvents,&GlobalEvents,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    EventCounter += GlobalEvents;
    fprintf(stderr,"The cumulative and total number of SNe %d and %g\n",EventCounter,TotalEventNumber);

    // if(Pall.TCurrent*Pall.UnitTime/YEAR_CGS > DURATOIN_TIME_AGNTORUS)
        // exit(1);

    return ;
}
#endif // TASK_AGNTORUS
