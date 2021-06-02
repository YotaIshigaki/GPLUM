#include "config.h"
#include "NeighborSearch.h"
#include "KernelFunctions.h"
#include "StellarFeedback.h"
#include "SizeDetermination.h"
#ifdef USE_RADIATION_PRESSURE //{
#include "RadiationPressure.h"
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
#include "StellarWind.h"
#endif // USE_STELLAR_WIND //}
#ifdef USE_MOMENTUM_FEEDBACK //{
#include "CalcEffectiveSurfaceArea.h"
#endif // USE_MOMENTUM_FEEDBACK //}

/*! \file StellarFeedback.c
 * \brief Feedback routines for type II/Ia SNe.
 */

//#define TARGET (413121)

static int StellarFeedbackExportFlagsMaxAllocated = 0;
struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticle; 

static int MaxIterationTimes = 20;
static bool OverMaxIterationTimes = false;
static double LocalKernelMax = 0.e0;

#define USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE
#define __CHECK_SUM__
#define __CHECK_WEIGHT__
// static int Test;

static double SNIINumberPerMass;
static double SNIIEnergy;
static double EnergyConvertToSimulationUnit;
static double InternalEnergyConvertToSimulationUnit;
static double MassConvertToSimulationUnit;
static double MassConvertToMsun;
static double KernelPeak;


static inline double __attribute__((always_inline)) KernelStellarFeedback(const double r, const double InvKerneli){

    double u = r*InvKerneli;
    double coef = M_1_PI*CUBE(InvKerneli);
    if(u<1.e0){
        return (coef*(1.e0 - 1.5*SQ(u) + 0.75*CUBE(u)));
    } else if (u<2.e0){
        return (coef*(0.25*CUBE(2.e0-u)));
    } else {
        return 0.e0;
    }
}

#ifdef USE_CELIB //{
/*!
 * This function initializes the whole functions of CELib.
 */
void InitializeStellarFeedback(void){

    CELibSetRunParameterIMFType(CHEMICALEVOLUTION_IMFTYPE);
    CELibSetRunParameterSNIIRange(CHEMICALEVOLUTION_SNII_UPPERMASS,CHEMICALEVOLUTION_SNII_LOWERMASS);
    CELibSetRunParameterSNIIYieldsTableID(CHEMICALEVOLUTION_SNII_YIELD_TYPE);
    CELibSetRunParameterSNIIHyperNovaFraction(CHEMICALEVOLUTION_SNII_HYPERNOVA_FRACTION);

    CELibSetRunParameterSNIIBinNumber(0);

#if CHEMICALEVOLUTION_SNIa_TYPE == -1
    CELibSetRunParameterSNIaType(1);
#else
    CELibSetRunParameterSNIaType(CHEMICALEVOLUTION_SNIa_TYPE);
#endif

    CELibSetRunParameterSNIaYieldsTableID(CHEMICALEVOLUTION_SNIa_YIELD_TYPE);
    CELibSetRunParameterSNIaYieldsTableModelID(CHEMICALEVOLUTION_SNIa_YIELD_MODEL);

#if CHEMICALEVOLUTION_SNIa_TYPE == 0
    CELibSetRunParameterSNIaRange(CHEMICALEVOLUTION_SNIa_UPPERMASS,CHEMICALEVOLUTION_SNIa_UPPERMASS);
#endif // CHEMICALEVOLUTION_SNIa_TYPE

    CELibSetRunParameterSNIaNassociation(CHEMICALEVOLUTION_SNIa_EVENTNUMBER);

#if CHEMICALEVOLUTION_POPIII_IMF == 0
    CELibSetRunParameterPopIIIIMF(0);
    CELibSetRunParameterPopIIISNe(0);
    CELibSetRunParameterPopIIIAGB(0);
    CELibSetRunParameterPopIIILifeTime(0);
#else
    CELibSetRunParameterPopIIIIMF(1);
    CELibSetRunParameterPopIIISNe(1);
    CELibSetRunParameterPopIIIAGB(1);
    CELibSetRunParameterPopIIILifeTime(1);
#endif


#ifdef USE_CELIB_AGB //{
    CELibSetRunParameterAGBBinNumber(CHEMICALEVOLUTION_AGB_NBIN);
    CELibSetRunParameterAGBBinType(CHEMICALEVOLUTION_AGB_BIN_TYPE);
    CELibSetRunParameterAGBBinTimeInterval(CHEMICALEVOLUTION_AGB_INTERVAL);
#endif // USE_CELIB_AGB //}

#ifdef USE_CELIB_NSM //{
    CELibSetRunParameterNSMDTDPowerLawIndex(CHEMICALEVOLUTION_NSM_DTD_INDEX);
    CELibSetRunParameterNSMDTDOffsetForPower(CHEMICALEVOLUTION_NSM_DTD_OFFSET);
#endif //USE_CELIB_NSM //}


    if(MPIGetMyID() == MPI_ROOT_RANK){
        CELibSetRunParameterTestMode(true);
    } else {
        CELibSetRunParameterTestMode(false);
    }
    
    // CELibSetRunParameterTestMode(false);
    // CELibSetRunParameterIntegrationSteps(10);

    CELibInit();

    if(MPIGetMyID() == MPI_ROOT_RANK){
        CELibShowVersion();
        CELibShowCurrentStatus();
    }

    EnergyConvertToSimulationUnit = GetUnitEnergy();
    InternalEnergyConvertToSimulationUnit = GetUnitSpecificEnergy();

    MassConvertToSimulationUnit = MSUN_CGS/Pall.UnitMass;
    MassConvertToMsun = Pall.UnitMass/MSUN_CGS;

    // SNIINumberPerMass = CELibGetSNIINumberPerMass();
    //SNIIEnergy = CELibGetSNIIEnergy();

#ifdef SNII_PEAK_TEMPERATURE //{
    KernelPeak = KernelStellarFeedback(0,1.0);
#endif // SNII_PEAK_TEMPERATURE //}

    return ;
}

int StellarFeedbackGetIMFType(void){
    return CELibGetRunParameterIMFType();
}

int StellarFeedbackGetSNIaType(void){
    return CELibGetRunParameterSNIaType();
}


static inline int __attribute__((always_inline)) CheckEventTime(const int Index){

#ifdef USE_CELIB_AGB // {
    if(Pstar[Index]->EventTimeSNII < 0.e0) 
        return NONE;

#if 0
    if(Pall.TCurrent > Pstar[Index]->EventTime){
        if(Pstar[Index]->SNIaCount == -1){ // II 
            return CELibFeedbackType_SNII;
        } else {
            return CELibFeedbackType_SNIa;
        }
    }
#else
    if(Pstar[Index]->SNIICount != NONE){ // II 
        if(Pall.TCurrent > Pstar[Index]->EventTimeSNII){
            return CELibFeedbackType_SNII;
        }
    } else {
        if(Pall.TCurrent > Pstar[Index]->EventTimeSNIa){
            return CELibFeedbackType_SNIa;
        }
    }
#endif

    if(Pall.TCurrent > Pstar[Index]->EventTimeAGB){
        return CELibFeedbackType_AGB;
    }

#ifdef USE_CELIB_NSM //{
    if(Pall.TCurrent > Pstar[Index]->EventTimeNSM){
        return CELibFeedbackType_NSM;
    }
#endif //USE_CELIB_NSM //}

    return NONE;
#else // USE_CELIB_AGB //}//{
    if(Pstar[Index]->EventTimeSNII < 0.e0) 
        return NONE;
#if 0
    if(Pall.TCurrent > Pstar[Index]->EventTime){
        if(Pstar[Index]->SNIaCount == -1){ // II 
            return CELibFeedbackType_SNII;
        } else {
            return CELibFeedbackType_SNIa;
        }
    }
#else
    if(Pstar[Index]->SNIICount != NONE){ // II 
        if(Pall.TCurrent > Pstar[Index]->EventTimeSNII){
            return CELibFeedbackType_SNII;
        }
    } else {
        if(Pall.TCurrent > Pstar[Index]->EventTimeSNIa){
            return CELibFeedbackType_SNIa;
        }
    }
#endif 
    return NONE;
#endif // USE_CELIB_AGB //}

}

int CountFeedbackNumber(void){
    int counter = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarBody(i)->Active){
            int CurrentFeedBackType = CheckEventTime(i);
            if(CurrentFeedBackType != NONE){
                counter ++;
            }
        }
    }
    return counter;
}


static int StellarFeedbackNContactedDomains;
static int *StellarFeedbackContactedDomainID;

//#define StellarFeedbackRadiusFactInc   (1.14) // 1.5 ^ (1.3)
#define StellarFeedbackRadiusFactInc   (1.2596) //  2^(0.333)
#define StellarFeedbackRadiusFactDec   (0.79) // 0.75 ^ (1.3)
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
static double SmoothedMassConversionFactor;
#endif //USE_SMOOTHED_NEIGHBOR_NUMBER

struct StructStellarFeedbackExport{
    double    Radius;  // Feedback radius (2xRadius is the actual feedback radius).
    double    Pos[3];  // Position.
    int       Leaf;
#ifdef __CHECK_WEIGHT__
    int       GlobalID;
#endif // __CHECK_WEIGHT__
};

struct StructStellarFeedbackImport{
    int    Leaf;
    int    Nlist; // Nlist.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double SmoothedNumber; // Smoothed mas..
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
    double Density;     // Mass weighted normalization factor.
#ifdef SET_SNII_TEMPERATURE
    double GasMass;     // Local gas mass
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
    double DistanceMin;
    unsigned long int DistanceMinGlobalID;
#endif
#ifdef __CHECK_SUM__ //{
    int CheckSum;
#endif // __CHECK_SUM__ //}
};

struct StructStellarFeedbackLocalInfo{
    int Nlist;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
    double SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER
    double Density;
#ifdef SET_SNII_TEMPERATURE
    double GasMass;// Local gas mass.
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
    double DistanceMin;
    unsigned long int DistanceMinGlobalID;
#endif
#ifdef __CHECK_SUM__ //{
    int CheckSum;
#endif // __CHECK_SUM__ //}
}; 

static int HydroUpdateFlagSize = 0;
static bool *HydroUpdateFlag = NULL;

static inline void __attribute__((always_inline)) StellarFeedbackAllocateContanctedDomainID(void){
    StellarFeedbackContactedDomainID = malloc(sizeof(int)*MPIGetNumProcs());
    return ;
}


static inline bool __attribute__((always_inline)) CheckLocalExternalDomainsContacted(const int MyID, const int ExtID){

    for(int k=0;k<3;k++){
        if((EdgesForHydro[MyID].PosMax[k] < EdgesForHydro[ExtID].PosMin[k])||
           (EdgesForHydro[MyID].PosMin[k] > EdgesForHydro[ExtID].PosMax[k]))  return false;
    }
    return true;

}


/*
 * This function checkes the number of contacted domains by comparing the local
 * domain edge to the external domains. 
 */
static inline void __attribute__((always_inline)) CheckContactedDomain(void){
    int NProcs = MPIGetNumProcs();
    int MyID = MPIGetMyID();
    StellarFeedbackNContactedDomains = 0;
    for(int i=0;i<NProcs-1;i++){
        int NodeID = CommunicationTable[i].SendRank;
        assert(MPIGetMyID() != NodeID);
        if(CheckLocalExternalDomainsContacted(MyID,NodeID)){
            StellarFeedbackContactedDomainID[StellarFeedbackNContactedDomains] = i;
            StellarFeedbackNContactedDomains ++;
        }
    }
    return ;
}


static inline bool __attribute__((always_inline)) OverlapDomainStellarFeedback(double Pos[restrict], const double h, const int NodeID){ 

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForHydro[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

static inline bool __attribute__((always_inline)) CheckInLocalDomain(double Pos[], double Kernel, const int MyID){
    for(int k=0;k<3;k++){
        if(Pos[k]+2.e0*Kernel > EdgesForHydro[MyID].PosMax[k]) return false;
        if(Pos[k]-2.e0*Kernel < EdgesForHydro[MyID].PosMin[k]) return false;
    }
    return true;
}


// #define __EXPORT_ALL__
static inline int __attribute__((always_inline)) CheckStellarFeedbackExportFlagsModified(const int NodeIndex, const int NProcs, bool StellarFeedbackExportFlags[restrict], const int NActives, struct StructActiveStellarFeedbackParticle ActiveStellarFeedbackParticle[restrict]){

    if(NActives == 0) 
        return 0;

    int ExportNodeID = CommunicationTable[NodeIndex].SendRank;
    // Node By Node Comparison
    double BoxCenter[] = {GravityNode[0].Pos[0],GravityNode[0].Pos[1],GravityNode[0].Pos[2]};
    if(!OverlapDomainStellarFeedback(BoxCenter,LocalKernelMax,ExportNodeID)){
        return 0;
    }

    int NExport = 0;
    for(int i=0;i<NActives;i++){
        int Offset = i*NProcs;
        //if(StellarFeedbackExportFlags[i][NProcs-1]){
        if(StellarFeedbackExportFlags[Offset+NProcs-1]){
#ifndef __EXPORT_ALL__
            if(OverlapDomainStellarFeedback(ActiveStellarFeedbackParticle[i].Pos,2.0*ActiveStellarFeedbackParticle[i].Radius,ExportNodeID)){
#endif // __EXPORT_ALL__
                //ExportFlags[i][NodeIndex] = true;
                StellarFeedbackExportFlags[Offset+NodeIndex] = true;
                NExport ++;
#ifndef __EXPORT_ALL__
            }
#endif // __EXPORT_ALL__
        }
    }

    return NExport;
}

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
static int GetSmoothedNumberofNeighbors(double Pos[restrict], const double Kernel, int Neighbors[restrict], double *SmoothedNumber){

    int nlist = 0;
    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    do {
        CurrentNodeID = GetNeighborsSmoothedNumberIterativeApproach(CurrentNodeID,
                Pos,2.e0*Kernel,&nlist,Neighbors,SmoothedNumber);
    }while(CurrentNodeID != RootNodeID);

    return nlist;
}
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}

/// Functions for the momentum feedback.
#ifdef USE_MOMENTUM_FEEDBACK //{

static inline double __attribute__((always_inline)) SNRMetalScalling(const double Z){
#define Zsun (0.01304)

    const double norm = (2.0/(pow(0.01,-0.14)));

    if(Z < 0.01*Zsun) return 2;
    else              return norm*pow(Z/Zsun,-0.14);
}

/*
 * This function returns the terminal momentum in the CGS units.
 */
static inline double __attribute__((always_inline)) GetPt(const double E_in_erg, const double nH, const double Z){
    const double ConvertToCGS = MSUN_CGS*VELOCITY_KMS_CGS;
    const double E_index = 13.0/14.0;
    const double nH_index = -1.0/7.0;
    double fz = SNRMetalScalling(Z);

    //return ConvertToCGS*4.8e5*pow(E_in_erg/1.e+51,E_index)*pow(fmax(nH,1.e-3),nH_index)*fz*sqrt(fz); //Msun*km/s
    //return ConvertToCGS*TERMINAL_MOMENTUM_EFFICIENCY*
        //4.8e5*pow(1.0,E_index)*pow(fmax(nH,1.e-3),nH_index)*fz*sqrt(fz); //Msun*km/s

    return ConvertToCGS*TERMINAL_MOMENTUM_EFFICIENCY*
        4.8e5*pow(E_in_erg/1.e+51,E_index)*pow(fmax(nH,1.e-3),nH_index)*fz*sqrt(fz); //Msun*km/s
    // return ConvertToCGS*TERMINAL_MOMENTUM_EFFICIENCY*
        // 4.8e5*pow(E_in_erg/1.e+51,E_index)*pow(nH,nH_index)*fz*sqrt(fz); //Msun*km/s
}


/*
 * This function returns the momentum corresponding to the input energy.
 */
static inline double __attribute__((always_inline)) GetPej(const double E_in_erg, const double M_in_Msun){
    return sqrt(2.0*M_in_Msun*MSUN_CGS*E_in_erg); // cgs
}

static void MomentumFeedbackKick(void){

    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            if((fabs(Phydro[i]->MomentumFB[0]) > 0.e0)
                    ||(fabs(Phydro[i]->MomentumFB[1]) > 0.e0)
                    ||(fabs(Phydro[i]->MomentumFB[2]) > 0.e0)){

                double p_old[] = {PhydroBody(i)->Mass*PhydroBody(i)->Vel[0],
                                  PhydroBody(i)->Mass*PhydroBody(i)->Vel[1],
                                  PhydroBody(i)->Mass*PhydroBody(i)->Vel[2]};

                // count new mass
                double NewMass = 0.e0;
                for(int k=0;k<CELibYield_Number;k++){
                    NewMass += Phydro[i]->Elements[k];
                }
                double iNewMass = 1.0/NewMass;

                double vold[] = {PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2]};

#if 0
                if(PhydroBody(i)->GlobalID == TARGET){
                        fprintf(stderr,"--Momentum increment \n");
                        fprintf(stderr,"o %g %g %g\n",
                            Phydro[i]->Mass*PhydroBody(i)->Vel[0],
                            Phydro[i]->Mass*PhydroBody(i)->Vel[1],
                            Phydro[i]->Mass*PhydroBody(i)->Vel[2]);

                        fprintf(stderr,"p %g %g %g\n",
                            Phydro[i]->MomentumFB[0],
                            Phydro[i]->MomentumFB[1],
                            Phydro[i]->MomentumFB[2]);
                        fflush(NULL);
                }
#endif

#ifdef USE_MOMENTUM_FEEDBACK_CARRY_OVER //{
                if(0.5*Phydro[i]->Mass*NORM(PhydroBody(i)->Vel)<NORM(Phydro[i]->MomentumFB)){
                    // carry over mode.
                    PhydroBody(i)->Vel[0] = (Phydro[i]->Mass*PhydroBody(i)->Vel[0]+MOMENTUM_FEEDBACK_CARRY_OVER_FRACTION*Phydro[i]->MomentumFB[0])*iNewMass;
                    PhydroBody(i)->Vel[1] = (Phydro[i]->Mass*PhydroBody(i)->Vel[1]+MOMENTUM_FEEDBACK_CARRY_OVER_FRACTION*Phydro[i]->MomentumFB[1])*iNewMass;
                    PhydroBody(i)->Vel[2] = (Phydro[i]->Mass*PhydroBody(i)->Vel[2]+MOMENTUM_FEEDBACK_CARRY_OVER_FRACTION*Phydro[i]->MomentumFB[2])*iNewMass;

                    Phydro[i]->MomentumFB[0] -= MOMENTUM_FEEDBACK_CARRY_OVER_FRACTION*Phydro[i]->MomentumFB[0];
                    Phydro[i]->MomentumFB[1] -= MOMENTUM_FEEDBACK_CARRY_OVER_FRACTION*Phydro[i]->MomentumFB[1];
                    Phydro[i]->MomentumFB[2] -= MOMENTUM_FEEDBACK_CARRY_OVER_FRACTION*Phydro[i]->MomentumFB[2];

                } else {
                    PhydroBody(i)->Vel[0] = (Phydro[i]->Mass*PhydroBody(i)->Vel[0]+Phydro[i]->MomentumFB[0])*iNewMass;
                    PhydroBody(i)->Vel[1] = (Phydro[i]->Mass*PhydroBody(i)->Vel[1]+Phydro[i]->MomentumFB[1])*iNewMass;
                    PhydroBody(i)->Vel[2] = (Phydro[i]->Mass*PhydroBody(i)->Vel[2]+Phydro[i]->MomentumFB[2])*iNewMass;
                    Phydro[i]->MomentumFB[0] = 
                    Phydro[i]->MomentumFB[1] = 
                    Phydro[i]->MomentumFB[2] = 0.e0;
                }
#ifdef USE_MAX_ALPHA_AT_SNREGIONS //{
                Phydro[i]->Alpha = VISCOSITY_MAX_AT_SNREGIONS;
#endif // USE_MAX_ALPHA_AT_SNREGIONS //}

#else // USE_MOMENTUM_FEEDBACK_CARRY_OVER //}//{


                //double iMass = 1.0/Phydro[i]->Mass;
                PhydroBody(i)->Vel[0] = (Phydro[i]->Mass*PhydroBody(i)->Vel[0]+Phydro[i]->MomentumFB[0])*iNewMass;
                PhydroBody(i)->Vel[1] = (Phydro[i]->Mass*PhydroBody(i)->Vel[1]+Phydro[i]->MomentumFB[1])*iNewMass;
                PhydroBody(i)->Vel[2] = (Phydro[i]->Mass*PhydroBody(i)->Vel[2]+Phydro[i]->MomentumFB[2])*iNewMass;

#if 1

                ////////////////////////////
                double v = NORM(PhydroBody(i)->Vel)*(Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;
                double vo = NORM(vold)*(Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;
                //double vh = NORM(Phydro[i]->VelP)*(Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;

                const double vel_unit = (Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;
                if(v > 1000){
                        fprintf(stderr,"--Momentum increment %ld \n",PhydroBody(i)->GlobalID);
                        fprintf(stderr,"o %g %g %g\n",
                            Phydro[i]->Mass*vold[0],
                            Phydro[i]->Mass*vold[1],
                            Phydro[i]->Mass*vold[2]);

                        fprintf(stderr,"p %g %g %g\n",
                            Phydro[i]->MomentumFB[0],
                            Phydro[i]->MomentumFB[1],
                            Phydro[i]->MomentumFB[2]);
                        fflush(NULL);

                    // if(v > 10000){
                        fprintf(stderr,"Vel[%ld] = %g km/s , %g km/s high speed! %s:%d\n",PhydroBody(i)->GlobalID,v,vo,__FUNCTION__,__LINE__);
                        double dmass = NewMass-Phydro[i]->Mass;
                        fprintf(stderr,"Mass %g %g / dmass = %g\n",Phydro[i]->Mass,NewMass,dmass);
                        fprintf(stderr,"Momentum %g %g %g\n",Phydro[i]->MomentumFB[0],Phydro[i]->MomentumFB[1],Phydro[i]->MomentumFB[2]);
                        if(dmass > 0.e0){
                            fprintf(stderr,"v %g %g %g\n",vel_unit*Phydro[i]->MomentumFB[0]/dmass,
                                    vel_unit*Phydro[i]->MomentumFB[1]/dmass,
                                    vel_unit*Phydro[i]->MomentumFB[2]/dmass);
                        }
                        fprintf(stderr,"vold %g %g %g\n",vel_unit*vold[0],vel_unit*vold[1],vel_unit*vold[2]);
                        fflush(NULL);
                    // } else {
                        //fprintf(stderr,"Vel[%ld] = %g km/s from %g km/s \n",PhydroBody(i)->GlobalID,v,vo);
                    // }
                }
                ////////////////////////////

#endif

#if 0
                //if(PhydroBody(i)->GlobalID == 13420)
                if(PhydroBody(i)->GlobalID == 5448){
                    fprintf(stderr,"pVel[%ld] = %g km/s , %g km/s high speed! %s:%d\n",PhydroBody(i)->GlobalID,v,vo,__FUNCTION__,__LINE__);
                    double dmass = NewMass-Phydro[i]->Mass;
                    fprintf(stderr,"Mass %g %g / dmass = %g\n",Phydro[i]->Mass,NewMass,dmass);
                    fprintf(stderr,"Momentum %g %g %g\n",Phydro[i]->MomentumFB[0],Phydro[i]->MomentumFB[1],Phydro[i]->MomentumFB[2]);
                    fprintf(stderr,"v %g %g %g\n",vel_unit*Phydro[i]->MomentumFB[0]/dmass,
                            vel_unit*Phydro[i]->MomentumFB[1]/dmass,
                            vel_unit*Phydro[i]->MomentumFB[2]/dmass);
                    fprintf(stderr,"vold %g %g %g\n",vel_unit*vold[0],vel_unit*vold[1],vel_unit*vold[2]);
                    fflush(NULL);
                }
#endif

                //Phydro[i]->MomentumFB[0] = Phydro[i]->MomentumFB[1] = Phydro[i]->MomentumFB[2] = 0.e0;
#ifdef USE_MAX_ALPHA_AT_SNREGIONS //{

                double p_new[] = {PhydroBody(i)->Mass*PhydroBody(i)->Vel[0],
                                  PhydroBody(i)->Mass*PhydroBody(i)->Vel[1],
                                  PhydroBody(i)->Mass*PhydroBody(i)->Vel[2]};

#if 0
                double f[] = {
                    fabs(p_old[0])>0?p_new[0]/p_old[0]:0.e0,
                    fabs(p_old[1])>0?p_new[1]/p_old[1]:0.e0,
                    fabs(p_old[2])>0?p_new[2]/p_old[2]:0.e0};

                double f_norm = 10.0*NORM(f);
                //Phydro[i]->Alpha = fmax(Phydro[i]->Alpha,f_norm);
                Phydro[i]->Alpha = fmax(Phydro[i]->Alpha,fmin(f_norm,VARIABLE_VISCOSITY_MAX));
#endif
#if 0
                //double f_norm = 10000*(NORM(p_old)>0?DISTANCE(p_new,p_old)/NORM(p_old):0.e0);
                Phydro[i]->Alpha = fmax(Phydro[i]->Alpha,fmin(sqrt(f_norm),VARIABLE_VISCOSITY_MAX));
#endif

#if 0
                // Mach number
                double v_norm = DISTANCE(p_new,p_old)/PhydroBody(i)->Mass;
                double cs = sqrt(Pall.GGm1*Phydro[i]->UPred);
                Phydro[i]->Alpha = fmax(Phydro[i]->Alpha,fmin(10.0*(v_norm+Phydro[i]->Vsig)/cs,VARIABLE_VISCOSITY_MAX));
#endif

                //Phydro[i]->Alpha = VARIABLE_VISCOSITY_MAX;
                Phydro[i]->Alpha = VISCOSITY_MAX_AT_SNREGIONS;
#endif // USE_MAX_ALPHA_AT_SNREGIONS //}


                Phydro[i]->MomentumFB[0] = 
                Phydro[i]->MomentumFB[1] = 
                Phydro[i]->MomentumFB[2] = 0.e0;
#endif // USE_MOMENTUM_FEEDBACK_CARRY_OVER //}
            }
        }
    }

    return ;
}

/*
 * This function evaluates Rsne. The value has been converted into the
 * simulation units.
 * See sec 2.3.2 in Hopkins et al. (2018), Modelling supernova feedback.
 * Note that the post-shock thermal energy outsize Rcool decays \propto (r/Rcool)^{-6/5}.
 */
static inline double __attribute__((always_inline)) GetRcool(const double E_in_erg, const double nH, const double Z){

    const double factor = PC_CGS/Pall.UnitLength;

    const double E_index = 2.0/7.0;
    const double nH_index = -3.0/7.0;

    double fz = SNRMetalScalling(Z);                   
    double Rcool = 28.4*pow(E_in_erg/1.e+51,E_index)*pow(nH,nH_index)*fz*sqrt(fz);  // pc
    return factor*Rcool;
}

#endif // USE_MOMENTUM_FEEDBACK //}


/*
 * This function returns a structure which involves various values (the neighbor
 * number, the smootehd neighbor number, the density, the local gas mass, the
 * minimum distance, and the globalID of the closest particle).
 */
struct StructStellarFeedbackLocalInfo RetrunStellarFeedbackLocalInfo(double Pos[restrict], const double Radius){

    static int Neighbors[MaxNeighborSize];
    struct StructStellarFeedbackLocalInfo TempStellarFeedbackLocalInfo = {
        .Nlist = 0,
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        .SmoothedNumber = 0.0,
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
        .Density = 0.e0,
#ifdef SET_SNII_TEMPERATURE //{
        .GasMass = 0.e0,
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT
        .DistanceMin = 0.e0,
        .DistanceMinGlobalID = NONE,
#endif
#ifdef __CHECK_SUM__ //{
        .CheckSum = 0,
#endif //__CHECK_SUM__ //}
    };

    double InvRadiusi = 1.e0/Radius;
    TempStellarFeedbackLocalInfo.Nlist = 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            GetSmoothedNumberofNeighbors(Pos,Radius,Neighbors,&(TempStellarFeedbackLocalInfo.SmoothedNumber));
#else // USE_SMOOTHED_NEIGHBOR_NUMBER 
            GetNeighborsLimited(Pos,2.0*Radius,Neighbors);
            //GetNeighborsDirect(Pos,2*Radius,Neighbors);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}

    for(int i=0;i<TempStellarFeedbackLocalInfo.Nlist;i++){
        int leaf = Neighbors[i];
        double xij[3];
#ifdef PERIODIC_RUN
        xij[0] = PeriodicDistance(Pos[0],PhydroPosP(leaf)[0],0);
        xij[1] = PeriodicDistance(Pos[1],PhydroPosP(leaf)[1],1);
        xij[2] = PeriodicDistance(Pos[2],PhydroPosP(leaf)[2],2);
#else // PERIODIC_RUN
        xij[0] = Pos[0]-PhydroPosP(leaf)[0];
        xij[1] = Pos[1]-PhydroPosP(leaf)[1];
        xij[2] = Pos[2]-PhydroPosP(leaf)[2];
#endif // PERIODIC_RUN

        double r = NORM(xij);
        double w = KernelStellarFeedback(r,InvRadiusi);

        TempStellarFeedbackLocalInfo.Density += PhydroMass(leaf)*w;
        //assert(PhydroMass(leaf)*w > 0.e0);
#ifdef SET_SNII_TEMPERATURE
        TempStellarFeedbackLocalInfo.GasMass += PhydroMass(leaf);
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
        if(i==0){
            TempStellarFeedbackLocalInfo.DistanceMinGlobalID = PhydroBody(leaf)->GlobalID;
            TempStellarFeedbackLocalInfo.DistanceMin = r;
        } else {
            if (TempStellarFeedbackLocalInfo.DistanceMin > r){
                TempStellarFeedbackLocalInfo.DistanceMinGlobalID = PhydroBody(leaf)->GlobalID;
                TempStellarFeedbackLocalInfo.DistanceMin = r;
            }
        }
#endif
#ifdef __CHECK_SUM__ //{
        TempStellarFeedbackLocalInfo.CheckSum += PhydroBody(leaf)->GlobalID;
#endif // __CHECK_SUM__ //}
    }

    return TempStellarFeedbackLocalInfo;
}


/*
 * This function checks the convergence of the feedback radius, using bisetion
 * method. If the radius is convergered, this function returns a "true" flag.
 * If not, it returns a "false" flag.
 */
static inline bool __attribute__((always_inline)) CheckNeighborNumberAndUpdateFeedbackRadius_i(const int NProcs, bool StellarFeedbackExportFlags_i[restrict], struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticle_i){ 

#ifdef USE_SN_INPUT_PARTICLE_NUMBER //{
    int NBmin = SN_INPUT_PARTICLE_NUMBER-SN_INPUT_PARTICLE_NUMBER_MARGIN;
    int NBmax = SN_INPUT_PARTICLE_NUMBER+SN_INPUT_PARTICLE_NUMBER_MARGIN;
#else // USE_SN_INPUT_PARTICLE_NUMBER
    int NBmin = Pall.Ns-Pall.Npm;
    int NBmax = Pall.Ns+Pall.Npm;
#endif // USE_SN_INPUT_PARTICLE_NUMBER //}

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double Nlist = ActiveStellarFeedbackParticle_i->SmoothedNumber*
        SmoothedMassConversionFactor*CUBE(ActiveStellarFeedbackParticle_i->Radius);
#else // USE_SMOOTHED_NEIGHBOR_NUMBER
    int Nlist = ActiveStellarFeedbackParticle_i->Nlist;
#endif  // USE_SMOOTHED_NEIGHBOR_NUMBER //}

    // Here, the convergence condition is satisfied. 
    if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
        StellarFeedbackExportFlags_i[NProcs-1] = false;
        return true;
    //}else if((Nlist<=(NBmax))&&(ActiveStellarFeedbackParticle_i->Rvalue>0.e0)&&(ActiveStellarFeedbackParticle_i->Lvalue>0.e0)){
    }else if(((0.8*NBmin)<=Nlist)&&(Nlist<=(2.0*NBmax))&&(ActiveStellarFeedbackParticle_i->Rvalue>0.e0)&&(ActiveStellarFeedbackParticle_i->Lvalue>0.e0)){
        if(ActiveStellarFeedbackParticle_i->Rvalue-ActiveStellarFeedbackParticle_i->Lvalue < 1.e-6*ActiveStellarFeedbackParticle_i->Lvalue){
            StellarFeedbackExportFlags_i[NProcs-1] = false;
            return true;
        }
    }
    if(StellarFeedbackExportFlags_i[NProcs-1]){
        if(Nlist<NBmin){
            ActiveStellarFeedbackParticle_i->Lvalue = fmax(ActiveStellarFeedbackParticle_i->Lvalue,ActiveStellarFeedbackParticle_i->Radius);
        } else if(Nlist>NBmax){
            if(ActiveStellarFeedbackParticle_i->Rvalue > 0.e0){
                ActiveStellarFeedbackParticle_i->Rvalue = fmin(ActiveStellarFeedbackParticle_i->Rvalue,ActiveStellarFeedbackParticle_i->Radius);
            }else{
                ActiveStellarFeedbackParticle_i->Rvalue = ActiveStellarFeedbackParticle_i->Radius;
            }
        }

        if((ActiveStellarFeedbackParticle_i->Lvalue>0.e0)&&(ActiveStellarFeedbackParticle_i->Rvalue>0.e0)){
            ActiveStellarFeedbackParticle_i->Radius = cbrt(0.5*(CUBE(ActiveStellarFeedbackParticle_i->Lvalue)+CUBE(ActiveStellarFeedbackParticle_i->Rvalue)));
        }else{
            if((ActiveStellarFeedbackParticle_i->Rvalue == 0.e0)&&(ActiveStellarFeedbackParticle_i->Lvalue > 0.e0)){
                ActiveStellarFeedbackParticle_i->Radius *= StellarFeedbackRadiusFactInc;
            }else if((ActiveStellarFeedbackParticle_i->Rvalue > 0.e0)&&(ActiveStellarFeedbackParticle_i->Lvalue == 0.e0)){
                ActiveStellarFeedbackParticle_i->Radius *= StellarFeedbackRadiusFactDec;
            }
        }
    }
    return false;
}


static inline void __attribute__((always_inline)) UpdateStellarFeedbackRadiusLocalModified(const int Index, const int NActives, 
        struct StructActiveStellarFeedbackParticle ActiveStellarFeedbackParticle[restrict], int Neighbors[restrict], const int MyID, const int NProcs, 
            bool StellarFeedbackExportFlags[restrict]){

    //if(CheckNeighborNumberAndUpdateFeedbackRadius_i(NProcs,
                //StellarFeedbackExportFlags[Index],ActiveStellarFeedbackParticle+Index) == true)
    if(CheckNeighborNumberAndUpdateFeedbackRadius_i(NProcs,
                StellarFeedbackExportFlags+Index*NProcs,ActiveStellarFeedbackParticle+Index) == true)
        return;

    do{
        if(!CheckInLocalDomain(ActiveStellarFeedbackParticle[Index].Pos,ActiveStellarFeedbackParticle[Index].Radius,MyID)) return;
        for(int i=0;i<StellarFeedbackNContactedDomains;i++){
            int NodeID = StellarFeedbackContactedDomainID[i];
            if(OverlapDomainStellarFeedback(ActiveStellarFeedbackParticle[Index].Pos,2.0*ActiveStellarFeedbackParticle[Index].Radius,NodeID)) return;
        }
        struct StructStellarFeedbackLocalInfo TemporalData 
            = RetrunStellarFeedbackLocalInfo(ActiveStellarFeedbackParticle[Index].Pos,ActiveStellarFeedbackParticle[Index].Radius);
        ActiveStellarFeedbackParticle[Index].Nlist = TemporalData.Nlist;
        ActiveStellarFeedbackParticle[Index].Density = TemporalData.Density;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        ActiveStellarFeedbackParticle[Index].SmoothedNumber = TemporalData.SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
        ActiveStellarFeedbackParticle[Index].GasMass = TemporalData.GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
        if(TemporalData.Nlist > 0){
            ActiveStellarFeedbackParticle[Index].DistanceMin = TemporalData.DistanceMin;
            ActiveStellarFeedbackParticle[Index].DistanceMinGlobalID = TemporalData.DistanceMinGlobalID;
        }
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef __CHECK_SUM__ //{
        ActiveStellarFeedbackParticle[Index].CheckSum = TemporalData.CheckSum;
#endif // __CHECK_SUM__ //}
        ActiveStellarFeedbackParticle[Index].IterationCount ++;
    }while(CheckNeighborNumberAndUpdateFeedbackRadius_i(NProcs,
                StellarFeedbackExportFlags+Index*NProcs,ActiveStellarFeedbackParticle+Index) == false);

    return;
}

static inline int __attribute__((always_inline)) CheckNeighborNumberAndUpdateFeedbackRadiusModified(const int NActives, const int NProcs, bool StellarFeedbackExportFlags[restrict], struct StructActiveStellarFeedbackParticle ActiveStellarFeedbackParticle[restrict]){ 

#ifdef USE_SN_INPUT_PARTICLE_NUMBER //{
    int NBmin = SN_INPUT_PARTICLE_NUMBER-SN_INPUT_PARTICLE_NUMBER_MARGIN;
    int NBmax;
    if(OverMaxIterationTimes){
        NBmax = SN_INPUT_PARTICLE_NUMBER+2*SN_INPUT_PARTICLE_NUMBER_MARGIN;
    } else {
        NBmax = SN_INPUT_PARTICLE_NUMBER+SN_INPUT_PARTICLE_NUMBER_MARGIN;
    }
#else // USE_SN_INPUT_PARTICLE_NUMBER
    int NBmin = Pall.Ns-Pall.Npm;
    int NBmax = Pall.Ns+Pall.Npm;
#endif // USE_SN_INPUT_PARTICLE_NUMBER //}

    int NLocalActiveLeaves = 0;
    for(int i=0;i<NActives;i++){
        int Offset = i*NProcs;
        //if(StellarFeedbackExportFlags[i][NProcs-1]){ 
        if(StellarFeedbackExportFlags[Offset+NProcs-1]){ 
#ifdef USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE //{
        if(ActiveStellarFeedbackParticle[i].LocalUpdateFlags == false){ 
#endif //USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE //}
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            double Nlist = ActiveStellarFeedbackParticle[i].SmoothedNumber*
                SmoothedMassConversionFactor*CUBE(ActiveStellarFeedbackParticle[i].Radius);
#else
            int Nlist = ActiveStellarFeedbackParticle[i].Nlist;
#endif //USE_MOOTHED_NEIGHBOR_NUMBER //}

            if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
                StellarFeedbackExportFlags[Offset+NProcs-1] = false;
            }else if(((0.8*NBmin)<=Nlist)&&(Nlist<=(2.0*NBmax))&&(ActiveStellarFeedbackParticle[i].Rvalue>0.e0)&&(ActiveStellarFeedbackParticle[i].Lvalue>0.e0)){
                if(ActiveStellarFeedbackParticle[i].Rvalue-ActiveStellarFeedbackParticle[i].Lvalue < 1.e-6*ActiveStellarFeedbackParticle[i].Lvalue){
                    StellarFeedbackExportFlags[Offset+NProcs-1] = false;
                }
            }
            if(StellarFeedbackExportFlags[Offset+NProcs-1]){
                if(Nlist<NBmin){
                    ActiveStellarFeedbackParticle[i].Lvalue = fmax(ActiveStellarFeedbackParticle[i].Lvalue,ActiveStellarFeedbackParticle[i].Radius);
                } else if(Nlist>NBmax){
                    if(ActiveStellarFeedbackParticle[i].Rvalue > 0.e0){
                        ActiveStellarFeedbackParticle[i].Rvalue = fmin(ActiveStellarFeedbackParticle[i].Rvalue,ActiveStellarFeedbackParticle[i].Radius);
                    }else{
                        ActiveStellarFeedbackParticle[i].Rvalue = ActiveStellarFeedbackParticle[i].Radius;
                    }
                }

                if((ActiveStellarFeedbackParticle[i].Lvalue>0.e0)&&(ActiveStellarFeedbackParticle[i].Rvalue>0.e0)){
                    ActiveStellarFeedbackParticle[i].Radius = cbrt(0.5*(CUBE(ActiveStellarFeedbackParticle[i].Lvalue)+CUBE(ActiveStellarFeedbackParticle[i].Rvalue)));
                }else{
                    if((ActiveStellarFeedbackParticle[i].Rvalue == 0.e0)&&(ActiveStellarFeedbackParticle[i].Lvalue > 0.e0)){
                        ActiveStellarFeedbackParticle[i].Radius *= StellarFeedbackRadiusFactInc;
                    }else if((ActiveStellarFeedbackParticle[i].Rvalue > 0.e0)&&(ActiveStellarFeedbackParticle[i].Lvalue == 0.e0)){
                        ActiveStellarFeedbackParticle[i].Radius *= StellarFeedbackRadiusFactDec;
                    }
                }
                NLocalActiveLeaves ++;
            }
            ActiveStellarFeedbackParticle[i].IterationCount ++;
#ifdef USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE //{
        } else {
            NLocalActiveLeaves ++;
        }
#endif //USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE //}
        }
    }

    return NLocalActiveLeaves;

}

static int MaxIterationTimesForInitialGuess = 10;
double StellarFeedbackRadiusInitialGuess(const int Index){

    int Neighbors[MaxNeighborSize];
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
    double Radius = 2.0*PstarBody(Index)->Eps; 

    if(Pall.Nhydro == 0)
        return Radius;

    double Pos[3] = {PstarBody(Index)->Pos[0],PstarBody(Index)->Pos[1],PstarBody(Index)->Pos[2]};
    int Iteration = 0;
    do{ 
        int Nlist = 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                GetSmoothedNumberofNeighbors(Pos,Radius,Neighbors,&SmoothedNumber);
#else // USE_SMOOTHED_NEIGHBOR_NUMBER 
                GetNeighborsLimited(Pos,2.0*Radius,Neighbors);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
        if(Nlist > 0){
            double Dist2 = DISTANCE2(Pos,NBCache[Neighbors[0]].Pos);
            double RadiusMin = NBCache[Neighbors[0]].Kernel;
            for(int i=1;i<Nlist;i++){
                int leaf = Neighbors[i];
                double CurrentDist2 = DISTANCE2(Pos,NBCache[leaf].Pos);
                if(Dist2 > CurrentDist2){
                    Dist2 = CurrentDist2;
                    RadiusMin = NBCache[leaf].Kernel;
                }
            }
            return RadiusMin;
        }
        Radius *= StellarFeedbackRadiusFactInc; 
#if 0
        if(Nlist > 2*Pall.Ns){
            double Dist2 = DISTANCE2(Pos,NBCache[Neighbors[0]].Pos);
            double RadiusMin = NBCache[Neighbors[0]].Kernel;
            for(int i=1;i<Nlist;i++){
                int leaf = Neighbors[i];
                double CurrentDist2 = DISTANCE2(Pos,NBCache[Neighbors[i]].Pos);
                if(Dist2 > CurrentDist2){
                    Dist2 = CurrentDist2;
                    RadiusMin = Radius;
                }
            }
            return RadiusMin;
        } else if(Nlist == 0){
            return 2.0*Radius;
        } else {
            return NBCache[Neighbors[0]].Kernel;
        }
#endif
        Iteration ++;
    } while(Iteration < MaxIterationTimesForInitialGuess);
    return 2.0*Radius;
}

#define _ResetTiming_ 2
static void ResetKernelSize(const int NActives, const int Niteration, const int NProcs, 
        const int IndexList[restrict], bool StellarFeedbackExportFlags[restrict]){

    if((Niteration+1)%(_ResetTiming_*MaxIterationTimes) == 0){
        for(int i=0;i<NActives;i++){
            int Offset = i*NProcs;
            if(StellarFeedbackExportFlags[Offset+NProcs-1]){ 
                ActiveStellarFeedbackParticle[i].Rvalue = ActiveStellarFeedbackParticle[i].Lvalue = 0.e0;
                ActiveStellarFeedbackParticle[i].Radius = 2.0*PstarBody(IndexList[i])->Eps
                                            *(gsl_rng_uniform(RandomGenerator)+0.5);
            }
        }
    }

}


//static bool **StellarFeedbackExportFlagsLog;
//static bool *StellarFeedbackExportFlagsLog; //! Export data log for stellar feedback.
static bool *StellarFeedbackExportFlags = NULL;

static bool first = true;
static void CalcFeedbackRadius(int NActives, const int IndexList[restrict], const int TypeList[restrict]){

    if(first){
        StellarFeedbackAllocateContanctedDomainID();
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
        SmoothedMassConversionFactor = (4.0*M_PI/3.0)*8.0;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER
        first = false;
    }

    OverMaxIterationTimes = false;

    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();
    int Neighbors[MaxNeighborSize];

    //static int StellarFeedbackExportFlagsMaxAllocated = 0;
    //static bool (*StellarFeedbackExportFlags)[NProcs];

    if(StellarFeedbackExportFlagsMaxAllocated < MAX(NActives,NAdditionUnit)){
        if(StellarFeedbackExportFlagsMaxAllocated > 0){
            free(StellarFeedbackExportFlags);
            free(ActiveStellarFeedbackParticle);
        }
        StellarFeedbackExportFlagsMaxAllocated = (int)(MAX(ForAngelsShare*NActives,NAdditionUnit));
        StellarFeedbackExportFlags = malloc(sizeof(bool)*StellarFeedbackExportFlagsMaxAllocated*NProcs);
        ActiveStellarFeedbackParticle = malloc(sizeof(struct StructActiveStellarFeedbackParticle)*StellarFeedbackExportFlagsMaxAllocated);
        //StellarFeedbackExportFlagsLog = (bool **)StellarFeedbackExportFlags;
        // StellarFeedbackExportFlagsLog = StellarFeedbackExportFlags;
    }

    for(int k=0;k<NActives;k++){
        int Offset = k*NProcs;
        StellarFeedbackExportFlags[Offset+NProcs-1] = ON;
        ActiveStellarFeedbackParticle[k].Index = IndexList[k];
        ActiveStellarFeedbackParticle[k].Pos[0] = PstarPos(IndexList[k])[0];
        ActiveStellarFeedbackParticle[k].Pos[1] = PstarPos(IndexList[k])[1]; 
        ActiveStellarFeedbackParticle[k].Pos[2] = PstarPos(IndexList[k])[2]; 
        //ActiveStellarFeedbackParticle[k].Radius = 2.0*PstarBody(IndexList[k])->Eps;
        ActiveStellarFeedbackParticle[k].Radius = StellarFeedbackRadiusInitialGuess(IndexList[k]);
        // ActiveStellarFeedbackParticle[k].TypeII = ActiveStellarFeedbackParticle[k].TypeIa = false;
        /*
        if(Pstar[IndexList[k]]->SNIaCount == -1){
            ActiveStellarFeedbackParticle[k].TypeII = true;
        } else {
            ActiveStellarFeedbackParticle[k].TypeIa = true;
        }
        */
        ActiveStellarFeedbackParticle[k].Type = TypeList[k];
        if(TypeList[k] == CELibFeedbackType_AGB){
            ActiveStellarFeedbackParticle[k].Count = Pstar[IndexList[k]]->AGBCount;
#ifdef USE_CELIB_NSM //{
        } else if(TypeList[k] == CELibFeedbackType_NSM){
            ActiveStellarFeedbackParticle[k].Count = Pstar[IndexList[k]]->NSMCount;
#endif // USE_CELIB_NSM //}
        } else {
            ActiveStellarFeedbackParticle[k].Count = Pstar[IndexList[k]]->SNIaCount;
        }
        ActiveStellarFeedbackParticle[k].InitialMass = Pstar[IndexList[k]]->InitialMass;
        ActiveStellarFeedbackParticle[k].Metallicity = Pstar[IndexList[k]]->Z;
        ActiveStellarFeedbackParticle[k].Rvalue = ActiveStellarFeedbackParticle[k].Lvalue = 0.e0;
        ActiveStellarFeedbackParticle[k].IterationCount = 0;
    }

    int BitMask = 0x01; 
    int NExportThisTime[NProcs-1];
    int NImportThisTime[NProcs-1];
    int NExportThisTimeNew[NProcs];
    int NImportThisTimeNew[NProcs];

    struct StructStellarFeedbackExport *StellarFeedbackExportSend[NProcs-1];
    struct StructStellarFeedbackExport *StellarFeedbackExportRecv = NULL;
    struct StructStellarFeedbackImport *StellarFeedbackImportSend = NULL;
    struct StructStellarFeedbackImport *StellarFeedbackImportRecv[NProcs-1];
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

#ifdef PRINT_LOG_KERNEL_ITERATION
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"Feedback Iteration");
#endif // PRINT_LOG_KERNEL_ITERATION

    int NActiveLeaves;
    MPI_Allreduce(&NActives,&NActiveLeaves,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    CheckContactedDomain();
    // Node center
    double BoxCenter[] = {GravityNode[0].Pos[0],GravityNode[0].Pos[1],GravityNode[0].Pos[2]};

    int Niteration = 0;
    do{
#ifdef PRINT_LOG_KERNEL_ITERATION
        if(MPIGetMyID()==MPI_ROOT_RANK)
            fprintf(stderr,":[%d] = %d ",Niteration,NActiveLeaves);

#endif // PRINT_LOG_KERNEL_ITERATION
        for(int i=0;i<NActives;i++){ 
            //int Offset = i*NActives;
            int Offset = i*NProcs;
            //if(StellarFeedbackExportFlags[i][NProcs-1]){
            if(StellarFeedbackExportFlags[Offset+NProcs-1]){
                if(Niteration > MaxIterationTimes){
// #ifdef SET_SNII_TEMPERATURE
                    // fprintf(stderr,"%ld %d %g %g %g\n",PstarBody(ActiveStellarFeedbackParticle[i].Index)->GlobalID,
                        // ActiveStellarFeedbackParticle[i].Nlist,ActiveStellarFeedbackParticle[i].GasMass,
                        // ActiveStellarFeedbackParticle[i].Lvalue,ActiveStellarFeedbackParticle[i].Rvalue);
// #else 
                    // fprintf(stderr,"%ld %d | %g %g\n",PstarBody(ActiveStellarFeedbackParticle[i].Index)->GlobalID,
                        // ActiveStellarFeedbackParticle[i].Nlist,
                        // ActiveStellarFeedbackParticle[i].Lvalue,ActiveStellarFeedbackParticle[i].Rvalue);
// #endif 
                    // fflush(NULL);
                }
                ActiveStellarFeedbackParticle[i].Nlist = 0;
                ActiveStellarFeedbackParticle[i].Density = 0.0;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                ActiveStellarFeedbackParticle[i].SmoothedNumber = 0.e0;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
                ActiveStellarFeedbackParticle[i].GasMass = 0.e0;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
                ActiveStellarFeedbackParticle[i].DistanceMin = 2.0*ActiveStellarFeedbackParticle[i].Radius;
#endif // MAXIMUM_ENERGY_INPUT //}
                for(int k=0;k<NProcs-1;k++)
                    StellarFeedbackExportFlags[Offset+k] = false;
                LocalKernelMax = fmax(LocalKernelMax,
                        DISTANCE(BoxCenter,ActiveStellarFeedbackParticle[i].Pos)+2.0*ActiveStellarFeedbackParticle[i].Radius);
#ifdef __CHECK_SUM__ //{
                ActiveStellarFeedbackParticle[i].CheckSum = 0;
#endif // __CHECK_SUM__ //}
            }
        }

        bool ExportFlags[NActives];

        for(int i=0;i<NProcs-1;i++){
            memset(ExportFlags,0,sizeof(bool)*NActives);

            NExportThisTime[i] = CheckStellarFeedbackExportFlagsModified(i,
                    NProcs,StellarFeedbackExportFlags,NActives,ActiveStellarFeedbackParticle);

            CheckSizeofBufferExportSendIndex(NExportThisTime[i],
                    sizeof(struct StructStellarFeedbackExport),i);
            CheckSizeofBufferImportRecvIndex(NExportThisTime[i],
                    sizeof(struct StructStellarFeedbackImport),i);
            StellarFeedbackExportSend[i] = BufferExportSend[i];
            StellarFeedbackImportRecv[i] = BufferImportRecv[i];

            int NExport = 0;
            if(NExportThisTime[i] > 0){
            for(int k=0;k<NActives;k++){
                int Offset = k*NProcs;
                if(StellarFeedbackExportFlags[Offset+NProcs-1]){
                    if(StellarFeedbackExportFlags[Offset+i]&BitMask){ 
                        StellarFeedbackExportSend[i][NExport].Pos[0] = ActiveStellarFeedbackParticle[k].Pos[0];
                        StellarFeedbackExportSend[i][NExport].Pos[1] = ActiveStellarFeedbackParticle[k].Pos[1];
                        StellarFeedbackExportSend[i][NExport].Pos[2] = ActiveStellarFeedbackParticle[k].Pos[2];
                        StellarFeedbackExportSend[i][NExport].Radius = ActiveStellarFeedbackParticle[k].Radius;
                        StellarFeedbackExportSend[i][NExport].Leaf = k;
#ifdef __CHECK_WEIGHT__
                        StellarFeedbackExportSend[i][NExport].GlobalID =
                            PstarBody(ActiveStellarFeedbackParticle[k].Index)->GlobalID;
#endif //__CHECK_WEIGHT__
                        NExport ++;
                    }
                }
            }
            }
            NExportThisTime[i] = NExport;
        }

        int NImportThisTime2[NProcs];
        int NExportThisTime2[NProcs];
        NImportThisTime2[MPIGetMyID()] = 0;
        for(int i=0;i<NProcs-1;i++){
            NExportThisTime2[CommunicationTable[i].SendRank] = NExportThisTime[i];
        }
        MPI_Alltoall(NExportThisTime2,1,MPI_INT,NImportThisTime2,1,MPI_INT,MPI_COMM_WORLD);
        int NImport = 0;
        for(int i=0;i<NProcs-1;i++){
            NImportThisTime[i] = NImportThisTime2[CommunicationTable[i].RecvRank];
            NImport += NImportThisTime[i];
        }
        int NImportAll = NImport;

        CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructStellarFeedbackExport));
        CheckSizeofBufferImportSend(NImportAll,sizeof(struct StructStellarFeedbackImport));
        StellarFeedbackExportRecv = BufferExportRecv;
        StellarFeedbackImportSend = BufferImportSend; 

        NImport = 0;
        int counter_send = 0;
        int counter_recv = 0;

        int SendFlag,RecvFlag;
        for(int i=0;i<NProcs-1;i++){
            if(NExportThisTime[i]>0){
                MPI_Isend(StellarFeedbackExportSend[i],
                    NExportThisTime[i]*sizeof(struct StructStellarFeedbackExport),
                        MPI_BYTE,CommunicationTable[i].SendRank,TAG_STELLARFEEDBACK_EXPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
                MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
                counter_send ++;
            }
            if(NImportThisTime[i]>0){
                MPI_Irecv(StellarFeedbackExportRecv+NImport,
                    NImportThisTime[i]*sizeof(struct StructStellarFeedbackExport),
                        MPI_BYTE,CommunicationTable[i].RecvRank,TAG_STELLARFEEDBACK_EXPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
                MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
                counter_recv ++;
                NImport += NImportThisTime[i];
            }
        }

        for(int i=0;i<NActives;i++){  // Check local
            int Offset = i*NProcs;
            //if(StellarFeedbackExportFlags[i][NProcs-1]){
            if(StellarFeedbackExportFlags[Offset+NProcs-1]){
                struct StructStellarFeedbackLocalInfo TemporalData = 
                    RetrunStellarFeedbackLocalInfo(ActiveStellarFeedbackParticle[i].Pos,ActiveStellarFeedbackParticle[i].Radius);
                ActiveStellarFeedbackParticle[i].Nlist = TemporalData.Nlist;
                ActiveStellarFeedbackParticle[i].Density = TemporalData.Density;

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                ActiveStellarFeedbackParticle[i].SmoothedNumber = TemporalData.SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
                ActiveStellarFeedbackParticle[i].GasMass = TemporalData.GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
                if(TemporalData.Nlist > 0){
                    ActiveStellarFeedbackParticle[i].DistanceMin = TemporalData.DistanceMin;
                    ActiveStellarFeedbackParticle[i].DistanceMinGlobalID = TemporalData.DistanceMinGlobalID;
                }
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef __CHECK_SUM__ //{
                ActiveStellarFeedbackParticle[i].CheckSum = TemporalData.CheckSum;
                // dprintlmpi(TemporalData.CheckSum);
#endif //__CHECK_SUM__ //}

#ifdef USE_STELLARFEEDBACK_RADIUS_LOCAL_UPDATE //{
                /// Insert Local Update Routine here.
                ActiveStellarFeedbackParticle[i].LocalUpdateFlags = false;
                int IsLocal = 0;
                for(int k=0;k<NProcs-1;k++){
                    //if(StellarFeedbackExportFlags[i][k]&BitMask){
                    if(StellarFeedbackExportFlags[Offset+k]&BitMask){
                        IsLocal ++;
                    }
                }
                if(IsLocal == 0){
                    ActiveStellarFeedbackParticle[i].LocalUpdateFlags = true;
                    //UpdateStellarFeedbackRadiusLocal(i,ActiveStellarFeedbackParticle,Neighbors,
                            //MyID,NProcs,StellarFeedbackExportFlags);
                    UpdateStellarFeedbackRadiusLocalModified(i,NActives,ActiveStellarFeedbackParticle,Neighbors,
                            MyID,NProcs,StellarFeedbackExportFlags);
                }
#endif // USE_KERNEL_LOCAL_UPDATE //}
            }
        }

        MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
        MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);

        for(int i=0;i<NImportAll;i++){ // For imported data.
            struct StructStellarFeedbackLocalInfo TemporalData = 
                RetrunStellarFeedbackLocalInfo(StellarFeedbackExportRecv[i].Pos,StellarFeedbackExportRecv[i].Radius);

            StellarFeedbackImportSend[i].Nlist = TemporalData.Nlist;
            StellarFeedbackImportSend[i].Density = TemporalData.Density;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            StellarFeedbackImportSend[i].SmoothedNumber = TemporalData.SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
            StellarFeedbackImportSend[i].GasMass = TemporalData.GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
            StellarFeedbackImportSend[i].DistanceMin = TemporalData.DistanceMin;
            StellarFeedbackImportSend[i].DistanceMinGlobalID = TemporalData.DistanceMinGlobalID;
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef __CHECK_SUM__ //{
            StellarFeedbackImportSend[i].CheckSum = TemporalData.CheckSum;
#endif //__CHECK_SUM__ //}
            StellarFeedbackImportSend[i].Leaf = StellarFeedbackExportRecv[i].Leaf;
        }

        NImportAll = 0;
        int NImportAllNew = 0;
        for(int i=0;i<NProcs-1;i++){
            NImportThisTimeNew[i] = 0;
            for(int k=0;k<NImportThisTime[i];k++){
                if(StellarFeedbackImportSend[NImportAll].Nlist > 0){
                    StellarFeedbackImportSend[NImportAllNew] = StellarFeedbackImportSend[NImportAll];
                    NImportThisTimeNew[i] ++;
                    NImportAllNew ++;
                }
                NImportAll ++;
            }
        }

        int NImportThisTimeNew2[NProcs];
        int NExportThisTimeNew2[NProcs];
        NImportThisTimeNew2[MPIGetMyID()] = 0;
        for(int i=0;i<NProcs-1;i++){
            NImportThisTimeNew2[CommunicationTable[i].SendRank] = NImportThisTimeNew[i];
        }
        MPI_Alltoall(NImportThisTimeNew2,1,MPI_INT,NExportThisTimeNew2,1,MPI_INT,MPI_COMM_WORLD);
        for(int i=0;i<NProcs-1;i++){
            NExportThisTimeNew[i] = NExportThisTimeNew2[CommunicationTable[i].RecvRank];
        }

        NImport = 0;
        counter_send = counter_recv = 0;
        for(int i=0;i<NProcs-1;i++){
            if(NImportThisTimeNew[i]>0){
                MPI_Isend(StellarFeedbackImportSend+NImport,
                    NImportThisTimeNew[i]*sizeof(struct StructStellarFeedbackImport),
                        MPI_BYTE,CommunicationTable[i].SendRank,TAG_STELLARFEEDBACK_IMPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
                MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
                counter_send ++;
            }
            if(NExportThisTimeNew[i]>0){
                MPI_Irecv(StellarFeedbackImportRecv[i],
                    NExportThisTimeNew[i]*sizeof(struct StructStellarFeedbackImport),
                        MPI_BYTE,CommunicationTable[i].RecvRank,TAG_STELLARFEEDBACK_IMPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
                MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
                counter_recv ++;
            }
            NImport += NImportThisTimeNew[i];
        }
        MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
        MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);

        for(int i=0;i<NProcs-1;i++){
            for(int k=0;k<NExportThisTimeNew[i];k++){ 
                int leaf = StellarFeedbackImportRecv[i][k].Leaf;
                ActiveStellarFeedbackParticle[leaf].Nlist += StellarFeedbackImportRecv[i][k].Nlist;
                ActiveStellarFeedbackParticle[leaf].Density += StellarFeedbackImportRecv[i][k].Density;
                
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                ActiveStellarFeedbackParticle[leaf].SmoothedNumber += StellarFeedbackImportRecv[i][k].SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
                ActiveStellarFeedbackParticle[leaf].GasMass += StellarFeedbackImportRecv[i][k].GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
                if(ActiveStellarFeedbackParticle[leaf].DistanceMin > StellarFeedbackImportRecv[i][k].DistanceMin){
                    ActiveStellarFeedbackParticle[leaf].DistanceMin = StellarFeedbackImportRecv[i][k].DistanceMin;
                    ActiveStellarFeedbackParticle[leaf].DistanceMinGlobalID = StellarFeedbackImportRecv[i][k].DistanceMinGlobalID;
                }
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef __CHECK_SUM__ //{
                ActiveStellarFeedbackParticle[leaf].CheckSum += StellarFeedbackImportRecv[i][k].CheckSum;
#endif //__CHECK_SUM__ //}
            }
        }

        // if(Niteration > MaxIterationTimes)
            // break;
        assert(Niteration < 1000);
        if(Niteration > MaxIterationTimes)
            OverMaxIterationTimes = true;
        ResetKernelSize(NActives,Niteration,NProcs,IndexList,StellarFeedbackExportFlags);

        int NLocalActiveLeaves = CheckNeighborNumberAndUpdateFeedbackRadiusModified(NActives,
                NProcs,StellarFeedbackExportFlags,ActiveStellarFeedbackParticle);
        MPI_Allreduce(&NLocalActiveLeaves,&NActiveLeaves,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

        Niteration ++;
    } while (0<NActiveLeaves);

#ifdef EVALUATE_KERNEL_BY_ITERATION //{
#ifdef PRINT_LOG_KERNEL_ITERATION //{
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"\n");
#else // PRINT_LOG_KERNEL_ITERATION 
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"%d iterations for kernel determination.\n",Niteration);
#endif // PRINT_LOG_KERNEL_ITERATION //}
#endif // EVALUATE_KERNEL_BY_ITERATION //}

    return;
}


struct StructEnergyHeavyElements{
    double Pos[3];
    double Radius;
    int    Type;
    // bool TypeII;
    // bool TypeIa;
    double Density;    // The mass weighted nomalization factor.
    int k_hydro_localmin;
#ifdef SET_SNII_TEMPERATURE
    double GasMass;// Local gas mass.
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
    double DistanceMin;
    unsigned long int DistanceMinGlobalID;
#endif
    struct CELibStructFeedbackOutput  CELibData;
#ifdef USE_RADIATION_PRESSURE //{
    double MomentumRP;
    double Metallicity;
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_MOMENTUM_FEEDBACK //{
    double Pt;
    double Pej;
    double NumberDensity;
    double SA;
    double WeightCorrection[EffSAVecSize];
    double Rcool;
    double Vel[3];
#endif // USE_MOMENTUM_FEEDBACK //}
    double WeightSum;    // The kernel normalization factor.
#ifdef  __CHECK_WEIGHT__
#ifdef USE_STAR_TIMESTEP_LIMITER //{ 
    int k;
#endif // USE_STAR_TIMESTEP_LIMITER //}
    int Leaf;
    int GlobalID;
#endif // __CHECK_WEIGHT__
};// *EnergyHeavyElements;

static inline void __attribute__((always_inline)) IncrementElements(double DestElements[restrict], double SrcElements[restrict], const double Weight){

    DestElements[CELibYield_H]  += SrcElements[CELibYield_H] *Weight;
    DestElements[CELibYield_He] += SrcElements[CELibYield_He]*Weight;
    DestElements[CELibYield_C]  += SrcElements[CELibYield_C] *Weight;
    DestElements[CELibYield_N]  += SrcElements[CELibYield_N] *Weight;
    DestElements[CELibYield_O]  += SrcElements[CELibYield_O] *Weight;
    DestElements[CELibYield_Ne] += SrcElements[CELibYield_Ne]*Weight;
    DestElements[CELibYield_Mg] += SrcElements[CELibYield_Mg]*Weight;
    DestElements[CELibYield_Si] += SrcElements[CELibYield_Si]*Weight;
    DestElements[CELibYield_S]  += SrcElements[CELibYield_S] *Weight;
    DestElements[CELibYield_Ca] += SrcElements[CELibYield_Ca]*Weight;
    DestElements[CELibYield_Fe] += SrcElements[CELibYield_Fe]*Weight;
    DestElements[CELibYield_Ni] += SrcElements[CELibYield_Ni]*Weight;
    DestElements[CELibYield_Eu] += SrcElements[CELibYield_Eu]*Weight;

    return;
}

long int GID;
long int GIDHydro;
int tmpLeaf;
static void CalcWeights(const double r, const double xij[restrict], const double InvKerneli, const double dwi, const double NumberDensityi, const double WeightCorrection[restrict], const double wk_norm, const double pw_norm, const int leaf, double *wk, double pvec[]){

    double InvKernelj = 1.0/Phydro[leaf]->KernelPred;
    double dwj = dSPHKernel(r,InvKernelj);

#ifdef USE_GRAD_N //{
    double NumberDensityj = Phydro[leaf]->NumberDensity;
#else // USE_GRAD_N //}//{
    double NumberDensityj = Phydro[leaf]->RhoPred/Phydro[leaf]->Mass; //1/Vi
#error 000
#endif // USE_GRAD_N //}
    double surface_area = fabs(dwi*r/SQ(NumberDensityi) + dwj*r/SQ(NumberDensityj));
    double wk_tmp = 0.5*(1.0-1.0/sqrt(1.0+surface_area/(M_PI*SQ(r))));

    double wk_vec[EffSAVecSize] = {0.e0};
    if(xij[0]>0){ 
        wk_vec[0] = wk_tmp*xij[0]/r;
        wk_vec[1] = 0.0;
    }else{
        wk_vec[0] = 0.0;
        wk_vec[1] = wk_tmp*xij[0]/r;
    }
    if(xij[1]>0){ 
        wk_vec[2] = wk_tmp*xij[1]/r;
        wk_vec[3] = 0.0;
    }else{
        wk_vec[2] = 0.0;
        wk_vec[3] = wk_tmp*xij[1]/r;
    }
    if(xij[2]>0){ 
        wk_vec[4] = wk_tmp*xij[2]/r;
        wk_vec[5] = 0.0;
    } else {
        wk_vec[4] = 0.0;
        wk_vec[5] = wk_tmp*xij[2]/r;
    }

    wk_tmp *= wk_norm;
    *wk = wk_tmp;

    double pnorm_tmp = 0.e0;
    double pvec_tmp[3] = {0.e0};
    for(int k=0;k<3;k++){
        double q = 0.e0;
        int i1 = 2*k;
        int i2 = i1+1;
        double q_i1 = fabs(WeightCorrection[i1]);
        double q_i2 = fabs(WeightCorrection[i2]);
        if((q_i1>0.e0)&&(q_i2>0.e0)){
            double rr = q_i2/q_i1;
            double rr2 = SQ(rr);
            //if(wk_vec[i1] != 0.e0){
            if(fpclassify(fabs(wk_vec[i1])) != FP_ZERO){
                q += wk_norm*wk_vec[i1]*sqrt(0.5*(1.e0+rr2));
            }else{
                q += wk_norm*wk_vec[i2]*sqrt(0.5*(1.e0+1.e0/rr2));
            }
        }else{
            q += wk_norm*(wk_vec[i1]+wk_vec[i2]);
        }
        pvec_tmp[k] = -q;
        pnorm_tmp += SQ(pvec_tmp[k]);
    }
    double pnorm = sqrt(pnorm_tmp);

    *wk = pnorm*pw_norm;
    for(int l=0;l<3;l++){
        //pvec[l] = pvec_tmp[l]/pnorm;
        pvec[l] = pvec_tmp[l]*pw_norm;
    }
    if(*wk > 1.0){
        *wk = 0.e0;
        for(int l=0;l<3;l++){
            pvec[l] = 0.e0;
        }
    }

    if(*wk > 1.0){
        fprintf(stderr,"Too large[%02d] %g *** // Leaf = %d %ld\n",MPIGetMyID(),*wk,tmpLeaf,PhydroBody(tmpLeaf)->GlobalID);
        fprintf(stderr,"## wk, wk_norm, pnorm, pw_norm = %g %g %g %g\n",*wk,wk_norm,pnorm,pw_norm);
        fprintf(stderr,"## wk_tmp = %g, surface = %g\n",
                0.5*(1.0-1.0/sqrt(1.0+surface_area/(M_PI*SQ(r)))),
                fabs(dwi*r/SQ(NumberDensityi) + dwj*r/SQ(NumberDensityj)));
        fprintf(stderr,"## dwi,r,NumberDensity %g %g %g  = %g / %g %g %g = %g\n",
                dwi,r,NumberDensityi,
                dwi*r/SQ(NumberDensityi),
                dwj,r,NumberDensityj,
                dwj*r/SQ(NumberDensityj));
        fprintf(stderr,"## dwi = %g r,InvKerneli = %g %g // %g\n",
                    dSPHKernel(r,InvKerneli),r,InvKerneli,r*InvKerneli);
        fprintf(stderr,"## dwj = %g r,InvKernelj = %g %g // %g\n",
                    dSPHKernel(r,InvKernelj),r,InvKernelj,r*InvKernelj);
        fprintf(stderr,"## wk_t = %g %g %g\n",*wk,pnorm,pw_norm);
        fprintf(stderr,"## 2hi/r,2hj/r = %g %g\n",2.0/InvKerneli/r,2.0/InvKernelj/r);
        fflush(NULL);
        assert(1==2);
    }
    if(*wk < 0.0){
        fprintf(stderr,"Negative %g \n",*wk);
    }
    return ;
}

static void Update_k_hydro_localmin(const int leaf, const int k_hydro_localmin){

    if(Phydro[leaf]->Active)
        return ;

    double EraMinimum = Pall.EraLocal+0.1*Pall.dtnow;
    //if(Phydro[leaf]->k_hydro - k_hydro_localmin > MAX_K_LOCAL_FB){
        //double dt_localmin = Pall.dtmin*exp2((double)(k_hydro_localmin+MAX_K_LOCAL_FB));
    if(Phydro[leaf]->k_hydro_localmin > k_hydro_localmin+MAX_K_LOCAL_FB){
    if(Phydro[leaf]->k_hydro > k_hydro_localmin+MAX_K_LOCAL_FB){
        Phydro[leaf]->k_hydro_localmin = k_hydro_localmin+MAX_K_LOCAL_FB;
        double dt_localmin = Pall.dtmin*exp2((double)(k_hydro_localmin+MAX_K_LOCAL_FB));
        int step = (int)(Pall.EraLocal/dt_localmin);
        double NextUpdateEra;
        do{
            NextUpdateEra = step*dt_localmin;
            step ++;
        }while(NextUpdateEra < EraMinimum);
        Phydro[leaf]->NextUpdateEra = fmin(Phydro[leaf]->NextUpdateEra,NextUpdateEra);
    }
    }
    return ;
}


/*
 * This function distributes the released energy and heavy elements to the
 * surrouding ISM.
 * The definition of the weight is as follows:
 *     Weight = Mass*w/Density.
 * This weight function guarantees that \sum Weight = 1, except a rounding
 * error.
 */
static void DistributeEnergyHeavyElements(const int NExplosion, struct StructEnergyHeavyElements EnergyHeavyElements[restrict], const bool Local){

    int Offset = ReturnCalcSizeOffset(CS_TypeSN);
    int Neighbors[MaxNeighborSize];
    for(int i=0;i<NExplosion;i++){
        double InvRadiusi = 1.e0/EnergyHeavyElements[i].Radius;
        double InvDensityi = 1.e0/EnergyHeavyElements[i].Density;
        double InvWeightSum = 1.e0/EnergyHeavyElements[i].WeightSum;
#ifdef USE_NEIGHBOR_LIST //{
        int Nlist;
        if(Local){
            struct StructGetLocalNeighborList NB = GetLocalNeighrborList(i+Offset);
            if(NB.Nlist>0){
                Nlist = NB.Nlist;
                for(int k=0;k<Nlist;k++){
                    Neighbors[k] = NB.Neighbors[k];
                }
            } else {
                Nlist = GetNeighborsLimited(EnergyHeavyElements[i].Pos,2.e0*EnergyHeavyElements[i].Radius,Neighbors);
            }
        } else {
            Nlist = GetNeighborsLimited(EnergyHeavyElements[i].Pos,2.e0*EnergyHeavyElements[i].Radius,Neighbors);
        }
#else // USE_NEIGHBOR_LIST //}//{

#ifdef USE_MOMENTUM_FEEDBACK //{
        int Nlist;
        
        if((EnergyHeavyElements[i].Type == StellarFeedbackType_SNII)
         ||(EnergyHeavyElements[i].Type == StellarFeedbackType_SNIa)
#ifdef USE_STELLAR_WIND //{
         ||(EnergyHeavyElements[i].Type == StellarFeedbackType_SW)
#endif // USE_STELLAR_WIND //}
           ){
            // Nlist = GetNeighborsPairsLimited(EnergyHeavyElements[i].Pos,2.e0*EnergyHeavyElements[i].Radius,Neighbors);
            Nlist = GetNeighborsPairsLimitedNBCacheIndex(EnergyHeavyElements[i].Pos,2.e0*EnergyHeavyElements[i].Radius,Neighbors);
        } else {
            // Nlist = GetNeighborsLimited(EnergyHeavyElements[i].Pos,2.e0*EnergyHeavyElements[i].Radius,Neighbors);
            Nlist = GetNeighborsLimitedNBCacheIndex(EnergyHeavyElements[i].Pos,2.e0*EnergyHeavyElements[i].Radius,Neighbors);
        }
#else // USE_MOMENTUM_FEEDBACK //}//{
        int Nlist = GetNeighborsLimited(EnergyHeavyElements[i].Pos,2.e0*EnergyHeavyElements[i].Radius,Neighbors);
#endif // // USE_MOMENTUM_FEEDBACK //}

#endif // USE_NEIGHBOR_LIST //}

#ifdef USE_MOMENTUM_FEEDBACK //{
        double NumberDensityi = EnergyHeavyElements[i].NumberDensity;
        double InvKerneli = 1.0/EnergyHeavyElements[i].Radius;
        double wk_norm = 1.0/(EnergyHeavyElements[i].SA);
        double pw_norm = 1.0/(EnergyHeavyElements[i].WeightCorrection[6]);

        if((EnergyHeavyElements[i].Type == StellarFeedbackType_SNII)
         ||(EnergyHeavyElements[i].Type == StellarFeedbackType_SNIa)
#ifdef USE_STELLAR_WIND //{
         ||(EnergyHeavyElements[i].Type == StellarFeedbackType_SW)
#endif // USE_STELLAR_WIND //}
           ){
            if(isinf(wk_norm)||isinf(pw_norm)){
                if(Local){
                    fprintf(stderr,"wk_norm/L %g / %g\n",wk_norm,EnergyHeavyElements[i].SA);
                    fprintf(stderr,"pw_norm/L %g / %g\n",pw_norm,EnergyHeavyElements[i].WeightCorrection[6]);
                } else {
                    fprintf(stderr,"wk_norm/E %g / %g\n",wk_norm,EnergyHeavyElements[i].SA);
                    fprintf(stderr,"pw_norm/E %g / %g\n",pw_norm,EnergyHeavyElements[i].WeightCorrection[6]);
                }
                fflush(NULL);
            }
        }



#endif // USE_MOMENTUM_FEEDBACK //}

        // dprintlmpi(Nlist);
        double wwww = 0.e0;
        for(int k=0;k<Nlist;k++){
            int index = Neighbors[k];
            int leaf = NBCache[index].Leaf;
            double xij[3],Posj[3];

            Posj[0] = NBCache[index].Pos[0];
            Posj[1] = NBCache[index].Pos[1];
            Posj[2] = NBCache[index].Pos[2];
#ifdef PERIODIC_RUN 
            xij[0] = PeriodicDistance(EnergyHeavyElements[i].Pos[0],Posj[0],0);
            xij[1] = PeriodicDistance(EnergyHeavyElements[i].Pos[1],Posj[1],1);
            xij[2] = PeriodicDistance(EnergyHeavyElements[i].Pos[2],Posj[2],2);
#else // PERIODIC_RUN 
            xij[0] = EnergyHeavyElements[i].Pos[0]-Posj[0];
            xij[1] = EnergyHeavyElements[i].Pos[1]-Posj[1];
            xij[2] = EnergyHeavyElements[i].Pos[2]-Posj[2];
#endif // PERIODIC_RUN

            double r = NORM(xij); 
            //double w = KernelStellarFeedback(r,InvRadiusi);
            double w = SPHKernel(r,InvRadiusi);

            double Weight = PhydroMass(leaf)*w*InvDensityi;
            // double m = PhydroMass(leaf);

            if(EnergyHeavyElements[i].Type == StellarFeedbackType_SNII){
#ifdef USE_MOMENTUM_FEEDBACK //{
                if(EnergyHeavyElements[i].CELibData.EjectaMass>0.e0){
                    double InvR = 1.0/r;

                    double wk,pvec[3];
                    double dwi = dSPHKernel(r,InvKerneli);
                    GID = EnergyHeavyElements[i].GlobalID;
                    GIDHydro = PhydroBody(leaf)->GlobalID;
                    tmpLeaf = leaf;
                    CalcWeights(r,xij,InvKerneli,dwi,NumberDensityi,EnergyHeavyElements[i].WeightCorrection,wk_norm,pw_norm,leaf,&wk,pvec);

                    wwww += wk;

                    //if(wk==0.e0) continue;
                    if(fpclassify(fabs(wk))==FP_ZERO) continue;
                    if(wk>1.e0){
                        fprintf(stderr,"Too large weight [%02d]  : %g %d\n",MPIGetMyID(),wk,Local?0:1);
                        fprintf(stderr,"GID %ld\n",EnergyHeavyElements[i].GlobalID);
                        fprintf(stderr,"//r = %g, xij = %g %g %g\n",r,xij[0],xij[1],xij[2]);
                        fprintf(stderr,"//InvKerneli, dwi, NumberDensityi = %g %g %g\n",InvKerneli,dwi,NumberDensityi);
                        fprintf(stderr,"//wk_norm = %g pw_nowm = %g wk = %g, pvec = %g %g %g\n",
                                wk_norm,pw_norm,wk,pvec[0],pvec[1],pvec[2]);
                        fflush(NULL);
                    }

                    double dm = wk*EnergyHeavyElements[i].CELibData.EjectaMass; // Simulation units

                    double E = SNE_EFFICIENCY*EnergyHeavyElements[i].CELibData.Energy; // Energy in simulation units
                    double E_in_erg = E/EnergyConvertToSimulationUnit; // Energy in erg
                    
                    double p_egy_j = sqrt(wk*2.0*(Phydro[leaf]->Mass+dm)*Pall.UnitMass*E_in_erg);
                    double Pt = GetPt(E_in_erg,
                            Pall.ConvertNumberDensityToCGS*Phydro[leaf]->RhoPred,Phydro[leaf]->Z);

#if 0
                    if(PhydroBody(leaf)->GlobalID == TARGET){
                     fprintf(stderr,"pre: %g / %g %g / %g %g / %g %g\n",wk,
                             dm,dm*Pall.UnitMass/MSUN_CGS,E,E_in_erg,p_egy_j,Pt);
                    }
#endif

                    // fprintf(stderr,"pre: %g %g / %g %g / %g %g / %g %g\n",wk,pnorm,
                            // dm,dm*Pall.UnitMass/MSUN_CGS,E,E_in_erg,p_egy_j,Pt);

                    const double VelFactor = Pall.UnitLength/Pall.UnitTime; // convert cgs units
		            double dmom[3] = {0.e0}; 
                    if(wk*Pt<p_egy_j){ 
                        // momentum conserving 
                        for(int l=0;l<3;l++){
                            //dmom[l] = pvec[l]*wk*Pt; 
                            dmom[l] = pvec[l]*Pt; 
                            dmom[l] += dm*Pall.UnitMass*EnergyHeavyElements[i].Vel[l]*VelFactor;
                        }
                        double dE_therm_j = wk*E; 
                        if(r>EnergyHeavyElements[i].Rcool){
                            dE_therm_j *= pow(EnergyHeavyElements[i].Rcool/r,6.5); 
                        }
                        Phydro[leaf]->DQheat += dE_therm_j;

#if 0
                    if(PhydroBody(leaf)->GlobalID == TARGET){
                         fprintf(stderr,"mom: %g %g, %g\n",
                                NORM(dmom),NORM2(dmom)/(2.0*(Phydro[leaf]->Mass+dm)*Pall.UnitMass),
                               dE_therm_j*EnergyConvertToSimulationUnit);
                    }
#endif
                        // fprintf(stderr,"mom: %g %g, %g\n",
                                // NORM(dmom),NORM2(dmom)/(2.0*(Phydro[leaf]->Mass+dm)*Pall.UnitMass),
                                // dE_therm_j*EnergyConvertToSimulationUnit);
                    }else{ 
                        // energy conserving
                        for(int l=0;l<3;l++){
                            //dmom[l] = pvec[l]*p_egy_j;
                            dmom[l] = pvec[l]*EnergyHeavyElements[i].WeightCorrection[6]*p_egy_j;
                            dmom[l] += dm*Pall.UnitMass*EnergyHeavyElements[i].Vel[l]*VelFactor;
                        }
#if 0
                    if(PhydroBody(leaf)->GlobalID == TARGET){
                         fprintf(stderr,"ene: %g %g\n",
                                 NORM(dmom),
                                 NORM2(dmom)/(2.0*(Phydro[leaf]->Mass+dm)*Pall.UnitMass));
                    }
#endif
                        // fprintf(stderr,"ene: %g %g\n",
                                // NORM(dmom),
                                // NORM2(dmom)/(2.0*(Phydro[leaf]->Mass+dm)*Pall.UnitMass));
                    }

                    const double MomentumConvertToSimulationUnit = 1.0/((Pall.UnitMass*Pall.UnitLength)/Pall.UnitTime);
#if 0
                    if(PhydroBody(leaf)->GlobalID == TARGET){
                        const double vel_unit = Pall.UnitLength/Pall.UnitTime/VELOCITY_KMS_CGS;
                        if(wk*Pt<p_egy_j){
                            fprintf(stderr,"Momentum conserving mode\n");
                            fprintf(stderr,"dm %g wk %g E %g, Rcool %g\n",dm,wk,E,EnergyHeavyElements[i].Rcool);
                            fprintf(stderr,"pvec %g %g %g, Pt %g\n",pvec[0],pvec[1],pvec[2],Pt);
                            fprintf(stderr,"vel %g %g %g\n",vel_unit*EnergyHeavyElements[i].Vel[0],
                                    vel_unit*EnergyHeavyElements[i].Vel[1],
                                    vel_unit*EnergyHeavyElements[i].Vel[2]);
                        } else {
                            fprintf(stderr,"Energy conserving mode\n");
                            fprintf(stderr,"wk %g M %g dm %g E %g\n",wk,Phydro[leaf]->Mass,dm,E_in_erg);
                            fprintf(stderr,"p_egy_j %g wk6 %g\n",p_egy_j,EnergyHeavyElements[i].WeightCorrection[6]);
                            fprintf(stderr,"pvec %g %g %g, Pt %g\n",pvec[0],pvec[1],pvec[2],Pt);
                            fprintf(stderr,"vel %g %g %g\n",vel_unit*EnergyHeavyElements[i].Vel[0],
                                    vel_unit*EnergyHeavyElements[i].Vel[1],
                                    vel_unit*EnergyHeavyElements[i].Vel[2]);
                            fprintf(stderr,"wk_norm %g pw_norm %g\n",wk_norm,pw_norm);
                        }
                        fprintf(stderr,"--Momentum \n");
                        fprintf(stderr,"p %g %g %g\n",
                            dmom[0]*MomentumConvertToSimulationUnit,
                            dmom[1]*MomentumConvertToSimulationUnit,
                            dmom[2]*MomentumConvertToSimulationUnit);

                        fprintf(stderr,"t %g %g %g\n",
                            Phydro[leaf]->MomentumFB[0] + dmom[0]*MomentumConvertToSimulationUnit,
                            Phydro[leaf]->MomentumFB[1] + dmom[1]*MomentumConvertToSimulationUnit,
                            Phydro[leaf]->MomentumFB[2] + dmom[2]*MomentumConvertToSimulationUnit);
                        fflush(NULL);
                    }
#endif

                    if(isnan(dmom[0])||isnan(dmom[1])||isnan(dmom[2])){
                        const double vel_unit = Pall.UnitLength/Pall.UnitTime/VELOCITY_KMS_CGS;
                        fprintf(stderr,"Momentum not a number(SNII)\n");
                        if(wk*Pt<p_egy_j){
                            fprintf(stderr,"Momentum conserving mode\n");
                            fprintf(stderr,"dm %g wk %g E %g, Rcool %g\n",dm,wk,E,EnergyHeavyElements[i].Rcool);
                            fprintf(stderr,"pvec %g %g %g, Pt %g\n",pvec[0],pvec[1],pvec[2],Pt);
                            fprintf(stderr,"vel %g %g %g\n",vel_unit*EnergyHeavyElements[i].Vel[0],
                                    vel_unit*EnergyHeavyElements[i].Vel[1],
                                    vel_unit*EnergyHeavyElements[i].Vel[2]);
                        } else {
                            fprintf(stderr,"Energy conserving mode\n");
                            fprintf(stderr,"wk %g M %g dm %g E %g\n",wk,Phydro[leaf]->Mass,dm,E_in_erg);
                            fprintf(stderr,"p_egy_j %g wk6 %g\n",p_egy_j,EnergyHeavyElements[i].WeightCorrection[6]);
                            fprintf(stderr,"pvec %g %g %g, Pt %g\n",pvec[0],pvec[1],pvec[2],Pt);
                            fprintf(stderr,"vel %g %g %g\n",vel_unit*EnergyHeavyElements[i].Vel[0],
                                    vel_unit*EnergyHeavyElements[i].Vel[1],
                                    vel_unit*EnergyHeavyElements[i].Vel[2]);
                            fprintf(stderr,"wk_norm %g pw_norm %g\n",wk_norm,pw_norm);
                        }
                        fflush(NULL);
                    } else {
                        Phydro[leaf]->MomentumFB[0] += dmom[0]*MomentumConvertToSimulationUnit;
                        Phydro[leaf]->MomentumFB[1] += dmom[1]*MomentumConvertToSimulationUnit;
                        Phydro[leaf]->MomentumFB[2] += dmom[2]*MomentumConvertToSimulationUnit;
                    }

                    IncrementElements(Phydro[leaf]->Elements,EnergyHeavyElements[i].CELibData.Elements,wk);

#ifdef USE_STAR_TIMESTEP_LIMITER //{
                    //Update_k_hydro_localmin(leaf,EnergyHeavyElements[i].k_hydro_localmin);
                    //dprintlmpi(EnergyHeavyElements[i].k);
                    Update_k_hydro_localmin(leaf,EnergyHeavyElements[i].k);
                    //Phydro[leaf]->k_hydro_localmin = MIN(EnergyHeavyElements[i].k+1,Phydro[leaf]->k_hydro_localmin);
#endif // USE_STAR_TIMESTEP_LIMITER //}
                }
#else // USE_MOMENTUM_FEEDBACK //}//{

#ifdef MAXIMUM_ENERGY_INPUT //{
                if(PhydroBody(leaf)->GlobalID == EnergyHeavyElements[i].DistanceMinGlobalID)
                    Phydro[leaf]->DQheat += EnergyHeavyElements[i].CELibData.Energy;
#else 
                Phydro[leaf]->DQheat += EnergyHeavyElements[i].CELibData.Energy*Weight;
#endif // MAXIMUM_ENERGY_INPUT //}

                IncrementElements(Phydro[leaf]->Elements,EnergyHeavyElements[i].CELibData.Elements,Weight);
#endif // USE_MOMENTUM_FEEDBACK //}


            } else if(EnergyHeavyElements[i].Type == StellarFeedbackType_SNIa){

#ifdef USE_MOMENTUM_FEEDBACK //{
                if(EnergyHeavyElements[i].CELibData.EjectaMass>0.e0){
                    double InvR = 1.0/r;

                    double wk,pvec[3];
                    double dwi = dSPHKernel(r,InvKerneli);
                    CalcWeights(r,xij,InvKerneli,dwi,NumberDensityi,EnergyHeavyElements[i].WeightCorrection,wk_norm,pw_norm,leaf,&wk,pvec);
                    //if(wk==0.e0) continue;
                    if(fpclassify(fabs(wk))==FP_ZERO) continue;

                    double dm = wk*EnergyHeavyElements[i].CELibData.EjectaMass; // Simulation units

                    double E = EnergyHeavyElements[i].CELibData.Energy; // Energy in simulation unit
                    double E_in_erg = E/EnergyConvertToSimulationUnit; // Energy in erg.
                    
                    double p_egy_j = sqrt(wk*2.0*(Phydro[leaf]->Mass+dm)*Pall.UnitMass*E_in_erg);
                    double Pt = GetPt(E_in_erg,
                            Pall.ConvertNumberDensityToCGS*Phydro[leaf]->RhoPred,Phydro[leaf]->Z);

                    const double VelFactor = Pall.UnitLength/Pall.UnitTime; // convert cgs units
		            double dmom[3] = {0.e0}; 
                    if(wk*Pt<p_egy_j){ 
                        // momentum conserving 
                        for(int l=0;l<3;l++){
                            //dmom[l] = pvec[l]*wk*Pt; 
                            dmom[l] = pvec[l]*Pt; 
                            dmom[l] += dm*Pall.UnitMass*EnergyHeavyElements[i].Vel[l]*VelFactor;
                        }
                        double dE_therm_j = wk*E; 
                        if (r>EnergyHeavyElements[i].Rcool){
                            dE_therm_j *= pow(EnergyHeavyElements[i].Rcool/r,6.5); 
                        }
                        Phydro[leaf]->DQheat += dE_therm_j;
                    }else{ 
                        // energy conserving
                        for(int l=0;l<3;l++){
                            //dmom[l] = pvec[l]*p_egy_j;
                            dmom[l] = pvec[l]*EnergyHeavyElements[i].WeightCorrection[6]*p_egy_j;
                            dmom[l] += dm*Pall.UnitMass*EnergyHeavyElements[i].Vel[l]*VelFactor;
                        }
                    }

                    const double MomentumConvertToSimulationUnit = 1.0/((Pall.UnitMass*Pall.UnitLength)/Pall.UnitTime);

#if 0
                    if(PhydroBody(leaf)->GlobalID == TARGET){
                        const double vel_unit = Pall.UnitLength/Pall.UnitTime/VELOCITY_KMS_CGS;
                        fprintf(stderr,"--TypeIa \n");
                        if(wk*Pt<p_egy_j){
                            fprintf(stderr,"Momentum conserving mode\n");
                            fprintf(stderr,"dm %g wk %g E %g, Rcool %g\n",dm,wk,E,EnergyHeavyElements[i].Rcool);
                            fprintf(stderr,"pvec %g %g %g, Pt %g\n",pvec[0],pvec[1],pvec[2],Pt);
                            fprintf(stderr,"vel %g %g %g\n",vel_unit*EnergyHeavyElements[i].Vel[0],
                                    vel_unit*EnergyHeavyElements[i].Vel[1],
                                    vel_unit*EnergyHeavyElements[i].Vel[2]);
                        } else {
                            fprintf(stderr,"Energy conserving mode\n");
                            fprintf(stderr,"wk %g M %g dm %g E %g\n",wk,Phydro[leaf]->Mass,dm,E_in_erg);
                            fprintf(stderr,"p_egy_j %g wk6 %g\n",p_egy_j,EnergyHeavyElements[i].WeightCorrection[6]);
                            fprintf(stderr,"pvec %g %g %g, Pt %g\n",pvec[0],pvec[1],pvec[2],Pt);
                            fprintf(stderr,"vel %g %g %g\n",vel_unit*EnergyHeavyElements[i].Vel[0],
                                    vel_unit*EnergyHeavyElements[i].Vel[1],
                                    vel_unit*EnergyHeavyElements[i].Vel[2]);
                            fprintf(stderr,"wk_norm %g pw_norm %g\n",wk_norm,pw_norm);
                        }
                        fprintf(stderr,"--Momentum \n");
                        fprintf(stderr,"p %g %g %g\n",
                            dmom[0]*MomentumConvertToSimulationUnit,
                            dmom[1]*MomentumConvertToSimulationUnit,
                            dmom[2]*MomentumConvertToSimulationUnit);

                        fprintf(stderr,"t %g %g %g\n",
                            Phydro[leaf]->MomentumFB[0] + dmom[0]*MomentumConvertToSimulationUnit,
                            Phydro[leaf]->MomentumFB[1] + dmom[1]*MomentumConvertToSimulationUnit,
                            Phydro[leaf]->MomentumFB[2] + dmom[2]*MomentumConvertToSimulationUnit);
                        fflush(NULL);
                    }
#endif 
                    if(isnan(dmom[0])||isnan(dmom[1])||isnan(dmom[2])){
                        const double vel_unit = Pall.UnitLength/Pall.UnitTime/VELOCITY_KMS_CGS;
                        fprintf(stderr,"Momentum not a number(SNIa)\n");
                        if(wk*Pt<p_egy_j){
                            fprintf(stderr,"Momentum conserving mode\n");
                            fprintf(stderr,"dm %g wk %g E %g, Rcool %g\n",dm,wk,E,EnergyHeavyElements[i].Rcool);
                            fprintf(stderr,"pvec %g %g %g, Pt %g\n",pvec[0],pvec[1],pvec[2],Pt);
                            fprintf(stderr,"vel %g %g %g\n",vel_unit*EnergyHeavyElements[i].Vel[0],
                                    vel_unit*EnergyHeavyElements[i].Vel[1],
                                    vel_unit*EnergyHeavyElements[i].Vel[2]);
                        } else {
                            fprintf(stderr,"Energy conserving mode\n");
                            fprintf(stderr,"wk %g M %g dm %g E %g\n",wk,Phydro[leaf]->Mass,dm,E_in_erg);
                            fprintf(stderr,"p_egy_j %g wk6 %g\n",p_egy_j,EnergyHeavyElements[i].WeightCorrection[6]);
                            fprintf(stderr,"pvec %g %g %g, Pt %g\n",pvec[0],pvec[1],pvec[2],Pt);
                            fprintf(stderr,"vel %g %g %g\n",vel_unit*EnergyHeavyElements[i].Vel[0],
                                    vel_unit*EnergyHeavyElements[i].Vel[1],
                                    vel_unit*EnergyHeavyElements[i].Vel[2]);
                            fprintf(stderr,"wk_norm %g pw_norm %g\n",wk_norm,pw_norm);
                        }
                        fflush(NULL);
                    } else {
                        Phydro[leaf]->MomentumFB[0] += dmom[0]*MomentumConvertToSimulationUnit;
                        Phydro[leaf]->MomentumFB[1] += dmom[1]*MomentumConvertToSimulationUnit;
                        Phydro[leaf]->MomentumFB[2] += dmom[2]*MomentumConvertToSimulationUnit;
                    }

                    IncrementElements(Phydro[leaf]->Elements,EnergyHeavyElements[i].CELibData.Elements,wk);
                    //Update_k_hydro_localmin(leaf,EnergyHeavyElements[i].k);
                }
#else // USE_MOMENTUM_FEEDBACK //}//{


#ifdef MAXIMUM_ENERGY_INPUT //{
                if(PhydroBody(leaf)->GlobalID == EnergyHeavyElements[i].DistanceMinGlobalID)
                    Phydro[leaf]->DQheat += EnergyHeavyElements[i].CELibData.Energy;
#else
                Phydro[leaf]->DQheat += EnergyHeavyElements[i].CELibData.Energy*Weight;
#endif // MAXIMUM_ENERGY_INPUT  //}


                IncrementElements(Phydro[leaf]->Elements,EnergyHeavyElements[i].CELibData.Elements,Weight);
#endif // USE_MOMENTUM_FEEDBACK //}
            } else if(EnergyHeavyElements[i].Type == StellarFeedbackType_AGB){
                IncrementElements(Phydro[leaf]->Elements,EnergyHeavyElements[i].CELibData.Elements,Weight);
            }
#ifdef USE_CELIB_NSM //{
            else if(EnergyHeavyElements[i].Type == StellarFeedbackType_NSM){
                IncrementElements(Phydro[leaf]->Elements,EnergyHeavyElements[i].CELibData.Elements,Weight);
            }
#endif // USE_CELIB_NSM //}
#ifdef USE_RADIATION_PRESSURE //{
            else if(EnergyHeavyElements[i].Type == StellarFeedbackType_RP){
                // Compute radiation feedback
                double _xij[] = {Phydro[leaf]->PosP[0]-EnergyHeavyElements[i].Pos[0],
                                 Phydro[leaf]->PosP[1]-EnergyHeavyElements[i].Pos[1],
                                 Phydro[leaf]->PosP[2]-EnergyHeavyElements[i].Pos[2]};

                double InvR = 1.0/NORM(_xij);
                double Factor = EnergyHeavyElements[i].MomentumRP*w*InvWeightSum;
                Phydro[leaf]->MomentumRP[0] += Factor*_xij[0]*InvR;
                Phydro[leaf]->MomentumRP[1] += Factor*_xij[1]*InvR;
                Phydro[leaf]->MomentumRP[2] += Factor*_xij[2]*InvR;
#if 1
                fprintf(stderr,"! MomentumRP %d/%d : %g %g %g/ %g %g %g/ %g\n",leaf,Nlist,
                        Phydro[leaf]->MomentumRP[0],Phydro[leaf]->MomentumRP[1],Phydro[leaf]->MomentumRP[2],
                        w,InvWeightSum,EnergyHeavyElements[i].WeightSum,
                        EnergyHeavyElements[i].Radius);
#endif
            }
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
            else if(EnergyHeavyElements[i].Type == StellarFeedbackType_SW){
#ifdef USE_MOMENTUM_FEEDBACK //{
                double wk,pvec[3];
                double dwi = dSPHKernel(r,InvKerneli);
                CalcWeights(r,xij,InvKerneli,dwi,NumberDensityi,EnergyHeavyElements[i].WeightCorrection,wk_norm,pw_norm,leaf,&wk,pvec);
                if(wk==0.e0) continue;

                double E = EnergyHeavyElements[i].CELibData.Energy; // Energy in simulation units
                double E_in_erg = E/EnergyConvertToSimulationUnit; // Energy in erg
                    
                double p_egy_j = sqrt(wk*2.0*(Phydro[leaf]->Mass)*Pall.UnitMass*E_in_erg);
                double Pt = GetPt(E_in_erg,
                        Pall.ConvertNumberDensityToCGS*Phydro[leaf]->RhoPred,Phydro[leaf]->Z);

                const double VelFactor = Pall.UnitLength/Pall.UnitTime; // convert cgs units
                double dmom[3] = {0.e0}; 
                if(wk*Pt<p_egy_j){ 
                    // momentum conserving 
                    for(int l=0;l<3;l++){
                        dmom[l] = pvec[l]*wk*Pt; 
                    }
                    double dE_therm_j = wk*E; 
                    if(r>EnergyHeavyElements[i].Rcool){
                        dE_therm_j *= pow(EnergyHeavyElements[i].Rcool/r,6.5); 
                    }
                    Phydro[leaf]->DQheat += dE_therm_j;
                }else{ 
                    // energy conserving
                    for(int l=0;l<3;l++){
                        dmom[l] = pvec[l]*p_egy_j;
                    }
                }
                const double MomentumConvertToSimulationUnit = 1.0/((Pall.UnitMass*Pall.UnitLength)/Pall.UnitTime);
                
#if 0
                    if(PhydroBody(leaf)->GlobalID == TARGET){
                        const double vel_unit = Pall.UnitLength/Pall.UnitTime/VELOCITY_KMS_CGS;
                        fprintf(stderr,"--SW \n");
                        fprintf(stderr,"--Momentum \n");
                        fprintf(stderr,"p %g %g %g\n",
                            dmom[0]*MomentumConvertToSimulationUnit,
                            dmom[1]*MomentumConvertToSimulationUnit,
                            dmom[2]*MomentumConvertToSimulationUnit);

                        fprintf(stderr,"t %g %g %g\n",
                            Phydro[leaf]->MomentumFB[0] + dmom[0]*MomentumConvertToSimulationUnit,
                            Phydro[leaf]->MomentumFB[1] + dmom[1]*MomentumConvertToSimulationUnit,
                            Phydro[leaf]->MomentumFB[2] + dmom[2]*MomentumConvertToSimulationUnit);
                        fflush(NULL);
                    }
#endif

                Phydro[leaf]->MomentumFB[0] += dmom[0]*MomentumConvertToSimulationUnit;
                Phydro[leaf]->MomentumFB[1] += dmom[1]*MomentumConvertToSimulationUnit;
                Phydro[leaf]->MomentumFB[2] += dmom[2]*MomentumConvertToSimulationUnit;

                //Update_k_hydro_localmin(leaf,EnergyHeavyElements[i].k);
#else // USE_MOMENTUM_FEEDBACK //}//{

#ifdef MAXIMUM_ENERGY_INPUT //{
                if(PhydroBody(leaf)->GlobalID == EnergyHeavyElements[i].DistanceMinGlobalID)
                    Phydro[leaf]->DQheat += EnergyHeavyElements[i].CELibData.Energy;
#else // MAXIMUM_ENERGY_INPUT //}//{
                Phydro[leaf]->DQheat += EnergyHeavyElements[i].CELibData.Energy*Weight;
#endif // MAXIMUM_ENERGY_INPUT  //}
#endif // USE_MOMENTUM_FEEDBACK //}

            }
#endif // USE_STELLAR_WIND //}

            // turn on Hydro Update Flag
            HydroUpdateFlag[leaf] = true;
        }
        // gprintlmpi(wwww);
    }

    return ;
}


static inline void __attribute__((always_inline)) ConvertSNIIEjectaintoUnitMass(double Elements[restrict]){

    for(int i=0;i<CELibYield_Number;i++){
        Elements[i] *= MSUN_CGS/Pall.UnitMass;
    }
    return ;
}


static inline void __attribute__((always_inline)) ConvertSNIaEjectaintoUnitMass(double Elements[restrict]){

    for(int i=0;i<CELibYield_Number;i++){
        Elements[i] *= MSUN_CGS/Pall.UnitMass;
    }
    return ;
}

static inline void __attribute__((always_inline)) ConvertEjectaUnits(struct CELibStructFeedbackOutput *CELibData, const int type){

    if((type == CELibFeedbackType_SNII)||(type == CELibFeedbackType_SNIa)){
        CELibData->Energy *= EnergyConvertToSimulationUnit;
    } else {
        CELibData->Energy = 0.e0;
    }
    double m = 0.e0;
    for(int i=0;i<CELibYield_Number;i++){
        CELibData->Elements[i] *= MassConvertToSimulationUnit;
        m += CELibData->Elements[i];
    }
    CELibData->EjectaMass *= MassConvertToSimulationUnit;
    CELibData->RemnantMass *= MassConvertToSimulationUnit;


    if(fabs(m-CELibData->EjectaMass) > 0.01*fabs(m)){
        fprintf(stderr,"%d:%d %g %g %d\n",MPIGetMyID(),__LINE__,m,CELibData->EjectaMass,type);
        fflush(NULL);
        assert(fabs(m-CELibData->EjectaMass) > 0.01*m);
    }


    return ;
}

static inline void __attribute__((always_inline)) UpdateStarParticleMass(double Mass, const int Index){

    double Fraction = (Pstar[Index]->Mass-Mass)/Pstar[Index]->Mass;
    for(int i=0;i<CELibYield_Number;i++){
        Pstar[Index]->Elements[i] *= Fraction;
    }
    double MassExcludeH = 0.e0;
    for(int i=1;i<CELibYield_Number;i++){
        MassExcludeH += Pstar[Index]->Elements[i];
    }

    Pstar[Index]->Mass -= Mass;
    PstarBody(Index)->Mass = Pstar[Index]->Mass;
    Pstar[Index]->Elements[CELibYield_H] = Pstar[Index]->Mass - MassExcludeH;

    return ;
}

static inline double __attribute__((always_inline)) GetSNIITemperatureBase(struct StructActiveStellarFeedbackParticle AP_i){

#if SET_SNII_TEMPERATURE_VARIABLE_MODEL==0 //{
    double log_density = log10(Pall.ConvertNumberDensityToCGS*AP_i.Density);
    const double TmaxLog = log10(5.e+8);
    const double TminLog = log10(5.e+6);
    if(log_density > +2){
        return 5.e+8;
    } else if (log_density < -2){
        return 5.e+6;
    } else {
        double grad = (TmaxLog-TminLog)/4;
        return pow(10.0,grad*(log_density+2.0)+TminLog);
    } 
#elif SET_SNII_TEMPERATURE_VARIABLE_MODEL==1 //}//{
    const double TminLog = log10(5.e+6);
    double range = 2.0*gsl_rng_uniform(RandomGenerator);
    return pow(10.0,TminLog+range);
#else // SET_SNII_TEMPERATURE_VARIABLE_MODEL==0 //}
#error Wrong value in SET_SNII_TEMPERATURE_VARIABLE_MODEL
    return 0.e0;
#endif // SET_SNII_TEMPERATURE_VARIABLE_MODEL==0 //}
}


/*
 * This function calculates the amounts of the release energy and heavy elements
 * of star particles, and then it scatters them to the surrounding ISM.
 */
static void ReleaseEnergyHeavyElements(const int NExplosion, const int StellarFeedBackList[restrict]){

    const int NProcs = MPIGetNumProcs();
    MPI_Status  mpi_status;

    struct StructEnergyHeavyElements EnergyHeavyElements[NExplosion+1];

    // Allocate HydroUpdateFlag
    if(Pall.Nhydro > HydroUpdateFlagSize){
        HydroUpdateFlagSize = ForAngelsShare*MAX(Pall.Nhydro,NAdditionUnit);
        free(HydroUpdateFlag);
        HydroUpdateFlag = malloc(sizeof(bool)*HydroUpdateFlagSize); 
    }

    for(int i=0;i<HydroUpdateFlagSize;i++)
        HydroUpdateFlag[i] = false;

    int Counter[StellarFeedbackType_Number] = {0};

#warning Something wrong regarding StellarFeedBackList
    int StellarFeedBackList_back[NExplosion];
    for(int i=0;i<NExplosion;i++){
        StellarFeedBackList_back[i] = StellarFeedBackList[i];
    }

    // Pack release energy/mass/elements.
    for(int i=0;i<NExplosion;i++){
        EnergyHeavyElements[i].Pos[0] = ActiveStellarFeedbackParticle[i].Pos[0];
        EnergyHeavyElements[i].Pos[1] = ActiveStellarFeedbackParticle[i].Pos[1];
        EnergyHeavyElements[i].Pos[2] = ActiveStellarFeedbackParticle[i].Pos[2];
        EnergyHeavyElements[i].Radius = ActiveStellarFeedbackParticle[i].Radius;
        EnergyHeavyElements[i].Density = ActiveStellarFeedbackParticle[i].Density;
        EnergyHeavyElements[i].k_hydro_localmin = ActiveStellarFeedbackParticle[i].k_hydro_localmin;
        EnergyHeavyElements[i].Type = ActiveStellarFeedbackParticle[i].Type;
#ifdef __CHECK_WEIGHT__
        EnergyHeavyElements[i].Leaf = ActiveStellarFeedbackParticle[i].Index;
        //EnergyHeavyElements[i].GlobalID = PstarBody(ActiveStellarFeedbackParticle[i].Index)->GlobalID;
        EnergyHeavyElements[i].GlobalID = PstarBody(ActiveStellarFeedbackParticle[i].Index)->GlobalID+
                                          Pstar[ActiveStellarFeedbackParticle[i].Index]->NthChildren;
#endif // __CHECK_WEIGHT__

#ifdef USE_STAR_TIMESTEP_LIMITER //{
        EnergyHeavyElements[i].k = ActiveStellarFeedbackParticle[i].k;
#endif // USE_STAR_TIMESTEP_LIMITER //}

#warning Something wrong regarding StellarFeedBackList
        //int leaf = StellarFeedBackList[i];
        int leaf = StellarFeedBackList_back[i];
        if(ActiveStellarFeedbackParticle[i].Type == StellarFeedbackType_SNII){ // type II
            EnergyHeavyElements[i].CELibData 
                = CELibGetFeedback((struct CELibStructFeedbackInput){
                        .Mass = ActiveStellarFeedbackParticle[i].InitialMass,
                        .MassConversionFactor = MassConvertToMsun,
                        .Metallicity = ActiveStellarFeedbackParticle[i].Metallicity,
                        .Elements = Pstar[leaf]->Elements,
                        .Count = Pstar[leaf]->SNIICount,
#ifdef USE_CELIB_SNII_INDIVIDUAL //{
                        .Mode = CELibSNIIRateModelID_Individual,
#else // USE_CELIB_SNII_INDIVIDUAL //}//{
                        .Mode = CELibSNIIRateModelID_Once,
#endif // USE_CELIB_SNII_INDIVIDUAL //}
                        },CELibFeedbackType_SNII);

            ConvertEjectaUnits(&(EnergyHeavyElements[i].CELibData),CELibFeedbackType_SNII);

            UpdateStarParticleMass(EnergyHeavyElements[i].CELibData.EjectaMass,leaf);

#ifdef USE_MOMENTUM_FEEDBACK //{
            double E_in_erg = EnergyHeavyElements[i].CELibData.Energy/EnergyConvertToSimulationUnit;
            double M_in_Msun = EnergyHeavyElements[i].CELibData.EjectaMass*Pall.UnitMass/MSUN_CGS;
            EnergyHeavyElements[i].Pej = GetPej(E_in_erg,M_in_Msun);
            EnergyHeavyElements[i].Pt = GetPt(E_in_erg,
                    ActiveStellarFeedbackParticle[i].nH_ave,ActiveStellarFeedbackParticle[i].Z_ave);

            const double MomentumConvertToSimulationUnit = 1.0/((Pall.UnitMass*Pall.UnitLength)/Pall.UnitTime);
            EnergyHeavyElements[i].Pej *= MomentumConvertToSimulationUnit;
            EnergyHeavyElements[i].Pt  *= MomentumConvertToSimulationUnit;

            EnergyHeavyElements[i].NumberDensity = ActiveStellarFeedbackParticle[i].NumberDensity;
            EnergyHeavyElements[i].SA = ActiveStellarFeedbackParticle[i].SA;
            if(isinf(1.0/ActiveStellarFeedbackParticle[i].SA)){
                fprintf(stderr,"INF [%02d] %d %d %g\n",
                        MPIGetMyID(),__LINE__,i,ActiveStellarFeedbackParticle[i].SA);
            }

            for(int l=0;l<EffSAVecSize;l++){
                EnergyHeavyElements[i].WeightCorrection[l] 
                    = ActiveStellarFeedbackParticle[i].WeightCorrection[l]; 
            }

            EnergyHeavyElements[i].Rcool = GetRcool(E_in_erg,
                    ActiveStellarFeedbackParticle[i].nH_ave,ActiveStellarFeedbackParticle[i].Z_ave);

            EnergyHeavyElements[i].Vel[0] = PstarBody(leaf)->Vel[0];
            EnergyHeavyElements[i].Vel[1] = PstarBody(leaf)->Vel[1];
            EnergyHeavyElements[i].Vel[2] = PstarBody(leaf)->Vel[2];
#endif //USE_MOMENTUM_FEEDBACK //}


#ifdef SET_SNII_TEMPERATURE //{
            EnergyHeavyElements[i].GasMass = ActiveStellarFeedbackParticle[i].GasMass;
#endif //SET_SNII_TEMPERATURE //}

#ifdef MAXIMUM_ENERGY_INPUT //{
            EnergyHeavyElements[i].DistanceMin = ActiveStellarFeedbackParticle[i].DistanceMin;
            EnergyHeavyElements[i].DistanceMinGlobalID = ActiveStellarFeedbackParticle[i].DistanceMinGlobalID;
#endif // MAXIMUM_ENERGY_INPUT //}


#if defined(PRESERVE_SNII_EVENTRATE) //{
            double MassInMsun = Pstar[leaf]->InitialMass*Pall.UnitMass/MSUN_CGS;
            if(Pstar[leaf]->TypeIIProb){
#error SNIINumberPerMass is not evaluated.
                double prob = SNIINumberPerMass*MassInMsun;
                if(prob >= 1.0){ // Do nothing!
                    EnergyHeavyElements[i].CELibData.Energy = SNIIEnergy*MassInMsun;
                } else { // Put energy of one SN.
                    EnergyHeavyElements[i].CELibData.Energy = SNIIEnergy*EnergyConvertToSimulationUnit;
                }
            } else {
                EnergyHeavyElements[i].CELibData.Energy = 0.e0;
            }
#elif defined(SET_SNII_TEMPERATURE) //} PRESERVE_SNII_EVENTRATE //{
            double MassInMsun = ActiveStellarFeedbackParticle[i].InitialMass*MassConvertToMsun;
            double Usn_in_sim_unit = EnergyHeavyElements[i].CELibData.Energy
                        /(EnergyHeavyElements[i].GasMass);
            double Tsn = Pall.ConvertUtoT*Usn_in_sim_unit;

#ifdef USE_SET_SNII_TEMPERATURE_VARIABLE //{
            double TsnBase = GetSNIITemperatureBase(ActiveStellarFeedbackParticle[i]);
            fprintf(stderr,"TsnBase : %g %g\n",Pall.ConvertNumberDensityToCGS*ActiveStellarFeedbackParticle[i].Density,TsnBase);
#else // USE_SET_SNII_TEMPERATURE_VARIABLE //}//{
            double TsnBase = SNII_TEMPERATURE;
#endif // USE_SET_SNII_TEMPERATURE_VARIABLE //}

#ifdef SNII_PEAK_TEMPERATURE //{
            double prob = Tsn/(KernelPeak*TsnBase);
#else //SNII_PEAK_TEMPERATURE //} //{
            double prob = Tsn/TsnBase;
#endif //SNII_PEAK_TEMPERATURE //}


            if(prob >= 1.0){ // Do nothing.
                EnergyHeavyElements[i].CELibData.Energy = EnergyConvertToSimulationUnit*SNIINumberPerMass*SNIIEnergy*MassInMsun;
            }else{
                if(prob > gsl_rng_uniform(RandomGenerator)){
#ifdef SNII_PEAK_TEMPERATURE //{
                    EnergyHeavyElements[i].CELibData.Energy = 
                        Pall.ConvertTtoU*KernelPeak*TsnBase*EnergyHeavyElements[i].GasMass;
#else //SNII_PEAK_TEMPERATURE //} //{
                    EnergyHeavyElements[i].CELibData.Energy = 
                        Pall.ConvertTtoU*TsnBase*EnergyHeavyElements[i].GasMass;
#endif //SNII_PEAK_TEMPERATURE //}
                }else{
                    EnergyHeavyElements[i].CELibData.Energy = 0.e0;
                }
            }
            fprintf(stderr,"p = %g, Tsn = %g [k], E = %g, %g [erg], Texp = %g [K]\n",
                    prob,Tsn,EnergyHeavyElements[i].CELibData.Energy,
                    EnergyHeavyElements[i].CELibData.Energy/EnergyConvertToSimulationUnit,
                    Pall.ConvertUtoT*
                    EnergyHeavyElements[i].CELibData.Energy/EnergyHeavyElements[i].GasMass
                    );
            fflush(NULL);
#endif // SET_SNII_TEMPERATURE //}


            EnergyHeavyElements[i].Type = StellarFeedbackType_SNII;
            Counter[StellarFeedbackType_SNII] ++;
        } else if (ActiveStellarFeedbackParticle[i].Type == StellarFeedbackType_SNIa){ // type Ia
            EnergyHeavyElements[i].CELibData 
                = CELibGetFeedback((struct CELibStructFeedbackInput){
                        .Mass = ActiveStellarFeedbackParticle[i].InitialMass,
                        .MassConversionFactor = MassConvertToMsun,
                        .Metallicity = ActiveStellarFeedbackParticle[i].Metallicity,
                        },CELibFeedbackType_SNIa);
                    
            ConvertEjectaUnits(&(EnergyHeavyElements[i].CELibData),CELibFeedbackType_SNIa);

#ifdef USE_MOMENTUM_FEEDBACK //{
            double E_in_erg = EnergyHeavyElements[i].CELibData.Energy/EnergyConvertToSimulationUnit;
            double M_in_Msun = EnergyHeavyElements[i].CELibData.EjectaMass/MassConvertToSimulationUnit;
            EnergyHeavyElements[i].Pej = GetPej(E_in_erg,M_in_Msun);
            EnergyHeavyElements[i].Pt = GetPt(E_in_erg,
                    ActiveStellarFeedbackParticle[i].nH_ave,ActiveStellarFeedbackParticle[i].Z_ave);
            const double MomentumConvertToSimulationUnit = 1.0/((Pall.UnitMass*Pall.UnitLength)/Pall.UnitTime);
            EnergyHeavyElements[i].Pej *= MomentumConvertToSimulationUnit;
            EnergyHeavyElements[i].Pt  *= MomentumConvertToSimulationUnit;

            EnergyHeavyElements[i].NumberDensity = ActiveStellarFeedbackParticle[i].NumberDensity;
            EnergyHeavyElements[i].SA = ActiveStellarFeedbackParticle[i].SA;
            for(int l=0;l<EffSAVecSize;l++){
                EnergyHeavyElements[i].WeightCorrection[l] 
                    = ActiveStellarFeedbackParticle[i].WeightCorrection[l]; 
            }

            EnergyHeavyElements[i].Rcool = GetRcool(E_in_erg,
                    ActiveStellarFeedbackParticle[i].nH_ave,ActiveStellarFeedbackParticle[i].Z_ave);

            EnergyHeavyElements[i].Vel[0] = PstarBody(leaf)->Vel[0];
            EnergyHeavyElements[i].Vel[1] = PstarBody(leaf)->Vel[1];
            EnergyHeavyElements[i].Vel[2] = PstarBody(leaf)->Vel[2];
#endif //USE_MOMENTUM_FEEDBACK //}

            if((PstarBody(leaf)->Mass > EnergyHeavyElements[i].CELibData.EjectaMass)&&(EnergyHeavyElements[i].CELibData.EjectaMass>0.e0)){
                UpdateStarParticleMass(EnergyHeavyElements[i].CELibData.EjectaMass,leaf);
            } else {
                EnergyHeavyElements[i].CELibData.EjectaMass = 0.e0;
                EnergyHeavyElements[i].CELibData.RemnantMass = Pstar[leaf]->Mass;
                for(int k=0;k<CELibYield_Number;k++){
                    EnergyHeavyElements[i].CELibData.Elements[k] = 0.e0;
                }
                fprintf(stderr,"\t//SNIa ejecta mass is larger than the stellar mass.\n");
                fprintf(stderr,"\t//This might be because the particle mass is too small or SNIaNassociation is too large.\n");
                fprintf(stderr,"\t//InitMass = %g Msun, CurrentMass = %g, Ejecta mass = %g Msun, SNIaNassociation = %d\n",
                        Pstar[leaf]->InitialMass*Pall.UnitMass/MSUN_CGS,
                        PstarBody(leaf)->Mass*Pall.UnitMass/MSUN_CGS,
                        EnergyHeavyElements[i].CELibData.EjectaMass,
                        CHEMICALEVOLUTION_SNIa_EVENTNUMBER);
                fprintf(stderr,"\t//Metallicity = %g\n",ActiveStellarFeedbackParticle[i].Metallicity);

            }

#ifdef MAXIMUM_ENERGY_INPUT //{
            EnergyHeavyElements[i].DistanceMin = ActiveStellarFeedbackParticle[i].DistanceMin;
            EnergyHeavyElements[i].DistanceMinGlobalID = ActiveStellarFeedbackParticle[i].DistanceMinGlobalID;
#endif // MAXIMUM_ENERGY_INPUT //}

            EnergyHeavyElements[i].Type = StellarFeedbackType_SNIa;
            Counter[StellarFeedbackType_SNIa] ++;
        } 
#ifdef USE_CELIB_AGB
        else if (ActiveStellarFeedbackParticle[i].Type == StellarFeedbackType_AGB){ // AGB
            EnergyHeavyElements[i].CELibData
                = CELibGetFeedback((struct CELibStructFeedbackInput){
                        .Mass = ActiveStellarFeedbackParticle[i].InitialMass,
                        .MassConversionFactor = MassConvertToMsun,
                        .Metallicity = ActiveStellarFeedbackParticle[i].Metallicity,
                        .Count = ActiveStellarFeedbackParticle[i].Count, 
                        .Elements = Pstar[leaf]->Elements,
                        },CELibFeedbackType_AGB);

            // Convert into the unit mass.
            ConvertEjectaUnits(&(EnergyHeavyElements[i].CELibData),CELibFeedbackType_AGB);
            
            if(EnergyHeavyElements[i].CELibData.EjectaMass>0.e0){
                if(PstarBody(leaf)->Mass > EnergyHeavyElements[i].CELibData.EjectaMass){
                    UpdateStarParticleMass(EnergyHeavyElements[i].CELibData.EjectaMass,leaf);
                } else {
                    double ej = EnergyHeavyElements[i].CELibData.EjectaMass;
                    EnergyHeavyElements[i].CELibData.EjectaMass = 0.e0;
                    EnergyHeavyElements[i].CELibData.RemnantMass = Pstar[leaf]->Mass;
                    for(int k=0;k<CELibYield_Number;k++){
                        EnergyHeavyElements[i].CELibData.Elements[k] = 0.e0;
                    }
                    fprintf(stderr,"\t//AGB ejecta mass is larger than the stellar mass.\n");
                    fprintf(stderr,"\t//This might be because the particle mass is too small or SNIaNassociation is too large.\n");
                    fprintf(stderr,"\t//InitMass = %g Msun, Current Mass = %g, Original ejecta mass = %g Msun, SNIaNassociation = %d\n",
                            Pstar[leaf]->InitialMass*Pall.UnitMass/MSUN_CGS,
                            PstarBody(leaf)->Mass*Pall.UnitMass/MSUN_CGS,
                            EnergyHeavyElements[i].CELibData.EjectaMass*Pall.UnitMass/MSUN_CGS,
                            CHEMICALEVOLUTION_SNIa_EVENTNUMBER);
                    fprintf(stderr,"\t//Metallicity = %g\n",ActiveStellarFeedbackParticle[i].Metallicity);

                }
            }

            EnergyHeavyElements[i].Type = StellarFeedbackType_AGB;
            Counter[StellarFeedbackType_AGB] ++;
        }
#endif // USE_CELIB_AGB
#ifdef USE_CELIB_NSM
        else if (ActiveStellarFeedbackParticle[i].Type == CELibFeedbackType_NSM){ // NSM
            EnergyHeavyElements[i].CELibData
                = CELibGetFeedback((struct CELibStructFeedbackInput){
                        .Mass = ActiveStellarFeedbackParticle[i].InitialMass,
                        .MassConversionFactor = MassConvertToMsun,
                        },CELibFeedbackType_NSM);

            // Convert into the unit mass.
            ConvertEjectaUnits(&(EnergyHeavyElements[i].CELibData),CELibFeedbackType_NSM);

            if((PstarBody(leaf)->Mass > EnergyHeavyElements[i].CELibData.EjectaMass)&&(EnergyHeavyElements[i].CELibData.EjectaMass>0.e0)){
                UpdateStarParticleMass(EnergyHeavyElements[i].CELibData.EjectaMass,leaf);
            } else {
                EnergyHeavyElements[i].CELibData.EjectaMass = 0.e0;
                EnergyHeavyElements[i].CELibData.RemnantMass = Pstar[leaf]->Mass;
                for(int k=0;k<CELibYield_Number;k++){
                    EnergyHeavyElements[i].CELibData.Elements[k] = 0.e0;
                }
                fprintf(stderr,"\t//NSM ejecta mass is larger than the stellar mass.\n");
                fprintf(stderr,"\t//This might be because the particle mass is too small or SNIaNassociation is too large.\n");
                fprintf(stderr,"\t//InitMass = %g Msun, Ejecta mass = %g Msun, SNIaNassociation = %d\n",
                        PstarBody(leaf)->Mass,EnergyHeavyElements[i].CELibData.EjectaMass,
                        CHEMICALEVOLUTION_SNIa_EVENTNUMBER);

            }

            EnergyHeavyElements[i].Type = StellarFeedbackType_NSM;
            Counter[StellarFeedbackType_NSM] ++;
        }
#endif // USE_CELIB_NSM
#ifdef USE_RADIATION_PRESSURE //{
        else if (ActiveStellarFeedbackParticle[i].Type == StellarFeedbackType_RP){ // RP
            UpdateBolometricLuminosity(leaf);

            const double Factor = Pall.UnitTime/(Pall.UnitMass*Pall.UnitLength);
            // Total released momentum in simulation units.
            EnergyHeavyElements[i].MomentumRP 
                = Factor*GetRadiationPressure(leaf)*PstarBody(leaf)->dt*Pall.UnitTime;

            EnergyHeavyElements[i].WeightSum = ActiveStellarFeedbackParticle[i].WeightSum;
            EnergyHeavyElements[i].Metallicity = ActiveStellarFeedbackParticle[i].Metallicity;

#if 1
            fprintf(stderr,"MomentumRP : %g / %g %g %g\n",
                    EnergyHeavyElements[i].MomentumRP,
                    EnergyHeavyElements[i].WeightSum,
                    EnergyHeavyElements[i].Radius,
                    EnergyHeavyElements[i].Metallicity);
#endif

            EnergyHeavyElements[i].Type = StellarFeedbackType_RP;
            Counter[StellarFeedbackType_RP] ++;
        }
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
        else if (ActiveStellarFeedbackParticle[i].Type == StellarFeedbackType_SW){ // SW
            UpdateStellarWindEnergy(leaf);
            // dprintlmpi(leaf);fflush(NULL);

            EnergyHeavyElements[i].CELibData.Energy = 
                Pstar[leaf]->StellarWindEnergy*EnergyConvertToSimulationUnit;
            // gprintlmpi(EnergyHeavyElements[i].CELibData.Energy);fflush(NULL);


#ifdef USE_MOMENTUM_FEEDBACK //{
            double E_in_erg = EnergyHeavyElements[i].CELibData.Energy/EnergyConvertToSimulationUnit;
            EnergyHeavyElements[i].Pej = 0.e0;
            EnergyHeavyElements[i].Pt = GetPt(E_in_erg,
                    ActiveStellarFeedbackParticle[i].nH_ave,ActiveStellarFeedbackParticle[i].Z_ave);
            const double MomentumConvertToSimulationUnit = 1.0/((Pall.UnitMass*Pall.UnitLength)/Pall.UnitTime);
            EnergyHeavyElements[i].Pej *= MomentumConvertToSimulationUnit;
            EnergyHeavyElements[i].Pt  *= MomentumConvertToSimulationUnit;

            EnergyHeavyElements[i].NumberDensity = ActiveStellarFeedbackParticle[i].NumberDensity;
            EnergyHeavyElements[i].SA = ActiveStellarFeedbackParticle[i].SA;
            for(int l=0;l<EffSAVecSize;l++){
                EnergyHeavyElements[i].WeightCorrection[l] 
                    = ActiveStellarFeedbackParticle[i].WeightCorrection[l]; 
            }

            EnergyHeavyElements[i].Rcool = GetRcool(E_in_erg,
                    ActiveStellarFeedbackParticle[i].nH_ave,ActiveStellarFeedbackParticle[i].Z_ave);

            EnergyHeavyElements[i].Vel[0] = PstarBody(leaf)->Vel[0];
            EnergyHeavyElements[i].Vel[1] = PstarBody(leaf)->Vel[1];
            EnergyHeavyElements[i].Vel[2] = PstarBody(leaf)->Vel[2];
#endif //USE_MOMENTUM_FEEDBACK //}


#ifdef MAXIMUM_ENERGY_INPUT //{
            EnergyHeavyElements[i].DistanceMin = ActiveStellarFeedbackParticle[i].DistanceMin;
            EnergyHeavyElements[i].DistanceMinGlobalID = ActiveStellarFeedbackParticle[i].DistanceMinGlobalID;
#endif // MAXIMUM_ENERGY_INPUT //}

            EnergyHeavyElements[i].WeightSum = ActiveStellarFeedbackParticle[i].WeightSum;

            EnergyHeavyElements[i].Type = StellarFeedbackType_SW;
            Counter[StellarFeedbackType_SW] ++;
            // dprintlmpi(Counter[StellarFeedbackType_SW]);fflush(NULL);
        }
#endif // USE_STELLAR_WIND //}
    }

    // Export
    struct StructEnergyHeavyElements *EnergyHeavyElementsExportSend[NProcs-1];
    struct StructEnergyHeavyElements *EnergyHeavyElementsExportRecv = NULL;
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];


    int NExportThisTime[NProcs-1];
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime[i] = 0;
        for(int k=0;k<NExplosion;k++){
            int Offset = k*NProcs;
            //if((StellarFeedbackExportFlags[Offset+i])
             //||(CheckExporFlagGatherScatter(k,CommunicationTable[i].SendRank))){
             if(CheckExporFlagGatherScatter(k,CommunicationTable[i].SendRank)){
                NExportThisTime[i] ++;
            }
        }

        CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructEnergyHeavyElements),i);
        EnergyHeavyElementsExportSend[i] = BufferExportSend[i];

        NExportThisTime[i] = 0; 
        for(int k=0;k<NExplosion;k++){
            int Offset = k*NProcs;
            //if((StellarFeedbackExportFlags[Offset+i])
             if(CheckExporFlagGatherScatter(k,CommunicationTable[i].SendRank)){
                EnergyHeavyElementsExportSend[i][NExportThisTime[i]] = EnergyHeavyElements[k];
                NExportThisTime[i] ++;
            }
        }
    }

    int NImport = 0;
    int NImportThisTime[NProcs-1];
    int NImportThisTime2[NProcs];
    int NExportThisTime2[NProcs];
    NImportThisTime2[MPIGetMyID()] = 0;
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime2[CommunicationTable[i].SendRank] = NExportThisTime[i];
    }
    MPI_Alltoall(NExportThisTime2,1,MPI_INT,NImportThisTime2,1,MPI_INT,MPI_COMM_WORLD);
    for(int i=0;i<NProcs-1;i++){
        NImportThisTime[i] = NImportThisTime2[CommunicationTable[i].RecvRank];
        NImport += NImportThisTime[i];
    }

    CheckSizeofBufferExportRecv(NImport,sizeof(struct StructEnergyHeavyElements));
    EnergyHeavyElementsExportRecv = BufferExportRecv;

    NImport = 0;
    int counter_send = 0;
    int counter_recv = 0;
    int SendFlag,RecvFlag;
    for(int i=0;i<NProcs-1;i++){
        if(NExportThisTime[i]>0){
            MPI_Isend(EnergyHeavyElementsExportSend[i],
                NExportThisTime[i]*sizeof(struct StructEnergyHeavyElements),
                MPI_BYTE,CommunicationTable[i].SendRank,TAG_STELLARFEEDBACK_EXPORT,
                    MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
        if(NImportThisTime[i]>0){
            MPI_Irecv(EnergyHeavyElementsExportRecv+NImport,
                NImportThisTime[i]*sizeof(struct StructEnergyHeavyElements),
                MPI_BYTE,CommunicationTable[i].RecvRank,TAG_STELLARFEEDBACK_EXPORT,
                    MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
            NImport += NImportThisTime[i];
        }
    }

    // Test = 0;
    // Distribution for local ISM particles
    DistributeEnergyHeavyElements(NExplosion,EnergyHeavyElements,true);

    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);

    // Test = 1;
    // Distribution for external ISM particles.
    DistributeEnergyHeavyElements(NImport,EnergyHeavyElementsExportRecv,false);

    return ;
}

static void StellarFeedbackEndProcedure(const int NExplosion, const int IndexList[restrict], const int TypeList[restrict]){

    const double FactorMass = Pall.UnitMass/MSUN_CGS;
    const double FactorTime = YEAR_CGS/Pall.UnitTime;

    // Turn on the TypeII flags.
    for(int i=0;i<NExplosion;i++){
        double NextEventTime;
        int leaf = IndexList[i];
        int Type; 
        double Metallicity = Pstar[leaf]->Z;
        double Mass_in_Msun = Pstar[leaf]->InitialMass*MassConvertToMsun;
        int Count;
#ifdef USE_STAR_TIMESTEP_LIMITER //{
        Pstar[leaf]->dt_fb = 0.e0;
#endif // USE_STAR_TIMESTEP_LIMITER //}

#ifdef USE_CELIB_SNII_INDIVIDUAL //{
        if(TypeList[i] == StellarFeedbackType_SNII){ // II/Ia
            Pstar[leaf]->SNIICount ++;
            Pstar[leaf]->EventTimeSNII = Pstar[leaf]->FormationTime
                +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                    .R = gsl_rng_uniform(RandomGenerator),
                    .InitialMass_in_Msun = FactorMass*Pstar[leaf]->InitialMass,
                    .Metallicity = Pstar[leaf]->Z,
                    .Count = Pstar[leaf]->SNIICount,
                    .Mode = CELibSNIIRateModelID_Individual,
                    .Nassociation = CHEMICALEVOLUTION_SNII_EVENTNUMBER,
                    },CELibFeedbackType_SNII)
                    *FactorTime;


#ifdef USE_FEEDBACK_TIMESTEP //{
            Pstar[leaf]->dt_fb = fmax(CHEMICALEVOLUTION_SNII_TIMEINTERVAL_MIN*FactorTime,
                    2.0*
                    (CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 1.e0,
                        .InitialMass_in_Msun = FactorMass*Pstar[leaf]->InitialMass,
                        .Metallicity = Pstar[leaf]->Z,
                        .Count = Pstar[leaf]->SNIICount,
                        .Mode = CELibSNIIRateModelID_Individual,
                        .Nassociation = CHEMICALEVOLUTION_SNII_EVENTNUMBER,
                        },CELibFeedbackType_SNII)-
                    CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = 0.e0,
                        .InitialMass_in_Msun = FactorMass*Pstar[leaf]->InitialMass,
                        .Metallicity = Pstar[leaf]->Z,
                        .Count = Pstar[leaf]->SNIICount,
                        .Nassociation = CHEMICALEVOLUTION_SNII_EVENTNUMBER,
                        .Mode = CELibSNIIRateModelID_Individual,
                        },CELibFeedbackType_SNII))
                        *FactorTime);
#endif // USE_FEEDBACK_TIMESTEP //}

    
            if(Pstar[leaf]->EventTimeSNII > Pall.TEnd){
                Pstar[leaf]->SNIICount = NONE;
                Pstar[leaf]->EventTimeSNII = 10*Pall.TEnd;
                Pstar[leaf]->SNIaCount = 0;
                Pstar[leaf]->EventTimeSNIa = Pstar[leaf]->FormationTime
                        +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                            .R = gsl_rng_uniform(RandomGenerator),
                            .InitialMass_in_Msun = FactorMass*Pstar[leaf]->InitialMass,
                            .Metallicity = Pstar[leaf]->Z,
                            .Count = Pstar[leaf]->SNIaCount,
                            },CELibFeedbackType_SNIa)
                            *FactorTime;
#ifdef USE_FEEDBACK_TIMESTEP //{
                Pstar[leaf]->dt_fb = 0.e0;
#endif // USE_FEEDBACK_TIMESTEP //}
            } 
        } else if(TypeList[i] == StellarFeedbackType_SNIa){ // Ia
            // fprintf(stderr,"@@@@@ increment? %d %d\n",Pstar[leaf]->SNIICount,Pstar[leaf]->SNIaCount);
            Pstar[leaf]->SNIaCount ++;
            Pstar[leaf]->EventTimeSNIa = Pstar[leaf]->FormationTime
                    +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = gsl_rng_uniform(RandomGenerator),
                        .InitialMass_in_Msun = FactorMass*Pstar[leaf]->InitialMass,
                        .Metallicity = Pstar[leaf]->Z,
                        .Count = Pstar[leaf]->SNIaCount,
                        },CELibFeedbackType_SNIa)
                        *FactorTime;
        }
#else // USE_CELIB_SNII_INDIVIDUAL //}//{
        if((TypeList[i] == StellarFeedbackType_SNII)||(TypeList[i] == StellarFeedbackType_SNIa)){ // II/Ia
            if(TypeList[i] == StellarFeedbackType_SNII){
                Pstar[leaf]->TypeII = true;
                Pstar[leaf]->SNIICount = NONE;
                Pstar[leaf]->SNIaCount = -1;
            } else if (TypeList[i] == StellarFeedbackType_SNIa){
                Pstar[leaf]->TypeIa = false;
            }

            Pstar[leaf]->SNIaCount ++;
            Pstar[leaf]->EventTimeSNIa = Pstar[leaf]->FormationTime
                    +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = gsl_rng_uniform(RandomGenerator),
                        .InitialMass_in_Msun = FactorMass*Pstar[leaf]->InitialMass,
                        .Metallicity = Pstar[leaf]->Z,
                        .Count = Pstar[leaf]->SNIaCount,
                        },CELibFeedbackType_SNIa)
                        *FactorTime;
        } 
#endif // USE_CELIB_SNII_INDIVIDUAL //}
        
        else if(TypeList[i] == StellarFeedbackType_AGB){ // AGB
            Pstar[leaf]->AGBCount ++;
            Pstar[leaf]->EventTimeAGB = Pstar[leaf]->FormationTime
                    +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = gsl_rng_uniform(RandomGenerator),
                        .InitialMass_in_Msun = FactorMass*Pstar[leaf]->InitialMass,
                        .Metallicity = Pstar[leaf]->Z,
                        .Count = Pstar[leaf]->AGBCount,
                        },CELibFeedbackType_AGB)
                        *FactorTime;
#ifdef USE_CELIB_NSM //{
        } else if(TypeList[i] == StellarFeedbackType_NSM){ // NSM
            Pstar[leaf]->NSMCount ++;
            Pstar[leaf]->EventTimeNSM = Pstar[leaf]->FormationTime
                    +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                        .R = gsl_rng_uniform(RandomGenerator),
                        .InitialMass_in_Msun = FactorMass*Pstar[leaf]->InitialMass,
                        .Metallicity = Pstar[leaf]->Z,
                        .Count = Pstar[leaf]->NSMCount,
                        },CELibFeedbackType_NSM)
                        *FactorTime;
#endif // USE_CELIB_NSM //}
        }
    }


    // Increment SPH particle mass.
    for(int i=0;i<Pall.Nhydro;i++){
        if(HydroUpdateFlag[i]){
            double Mass = 0.e0;
            for(int k=0;k<CELibYield_Number;k++){
                Mass += Phydro[i]->Elements[k];
            }
            if(Mass < 0.e0){
                fprintf(stderr,"GID = %ld\n",PhydroBody(i)->GlobalID);
                fprintf(stderr,"Mass = %g -> %g\n",Phydro[i]->Mass,Mass);
                for(int k=0;k<CELibYield_Number;k++){
                    fprintf(stderr,"Elements %g\n",Phydro[i]->Elements[k]);
                }
            }

            Phydro[i]->Mass = PhydroBody(i)->Mass = Mass;

            assert(Phydro[i]->Mass > 0.e0);

            double MassLightElements = 
                Phydro[i]->Elements[CELibYield_H]+Phydro[i]->Elements[CELibYield_He];
            Phydro[i]->Z = (Mass-MassLightElements)/Mass;
        }
    }

    return ;
}

static void AllocateActiveStellarFeedbackParticle(const int NActives){

    int NProcs = MPIGetNumProcs();
    if(StellarFeedbackExportFlagsMaxAllocated < MAX(NActives,NAdditionUnit)){
        StellarFeedbackExportFlagsMaxAllocated = (int)(MAX(ForAngelsShare*NActives,NAdditionUnit));
        StellarFeedbackExportFlags = realloc(StellarFeedbackExportFlags,
                sizeof(bool)*StellarFeedbackExportFlagsMaxAllocated*NProcs);
        ActiveStellarFeedbackParticle = realloc(ActiveStellarFeedbackParticle,
                sizeof(struct StructActiveStellarFeedbackParticle)*StellarFeedbackExportFlagsMaxAllocated);
    }

    return ;
}


void StellarFeedback(void){

#ifdef EVALUATE_SIZES_ALL_TOGETHER //{

    if(ReturnCalcSizeElementNumber(CS_TypeSN,true)
#ifdef USE_RADIATION_PRESSURE //{
            +ReturnCalcSizeElementNumber(CS_TypeRP,true) 
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
            +ReturnCalcSizeElementNumber(CS_TypeSW,true) 
#endif // USE_STELLAR_WIND //}
            == 0){
        return ;
    }

    int NExplosionSN = ReturnCalcSizeElementNumber(CS_TypeSN,false);
#ifdef USE_RADIATION_PRESSURE //{
    int NExplosionRP = ReturnCalcSizeElementNumber(CS_TypeRP,false);
#else // USE_RADIATION_PRESSURE //}//{
    int NExplosionRP = 0;
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
    int NExplosionSW = ReturnCalcSizeElementNumber(CS_TypeSW,false);
#else // USE_STELLAR_WIND //}//{
    int NExplosionSW = 0;
#endif // USE_STELLAR_WIND //}
    int NExplosion = NExplosionSN + NExplosionRP + NExplosionSW;

#if 0
/////////////////////////////////////////////////////////


#ifdef USE_RADIATION_PRESSURE //{
    if(ReturnCalcSizeElementNumber(CS_TypeSN,true)+ReturnCalcSizeElementNumber(CS_TypeRP,true) == 0){
        return ;
    }
    int NExplosionSN = ReturnCalcSizeElementNumber(CS_TypeSN,false);
    int NExplosionRP = ReturnCalcSizeElementNumber(CS_TypeRP,false);
    int NExplosion = NExplosionSN + NExplosionRP;
#if 1
    if(NExplosion > 0){
        dprintlmpi(NExplosionSN);
        dprintlmpi(NExplosionRP);
        dprintlmpi(NExplosion);
    }
#endif
#else // USE_RADIATION_PRESSURE //}//{
    int NExplosion = ReturnCalcSizeElementNumber(CS_TypeSN,false);
    if(ReturnCalcSizeElementNumber(CS_TypeSN,true) == 0){
        return ;
    }
#endif  // USE_RADIATION_PRESSURE //}

/////////////////////////////////////////////////////////
#endif


    int StellarFeedBackList[NExplosion+1];
    int FeedBackType[NExplosion+1];

    AllocateActiveStellarFeedbackParticle(NExplosion);
    CalcSizeSetSNInfo(ActiveStellarFeedbackParticle);
#ifdef USE_RADIATION_PRESSURE //{
    CalcSizeSetRPInfo(ActiveStellarFeedbackParticle+NExplosionSN);
#endif  // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
    CalcSizeSetSWInfo(ActiveStellarFeedbackParticle+NExplosionSN+NExplosionRP);
#endif  // USE_STELLAR_WIND //}
    for(int i=0;i<NExplosion;i++){
        StellarFeedBackList[i] = ActiveStellarFeedbackParticle[i].Index;
        FeedBackType[i] = ActiveStellarFeedbackParticle[i].Type;
    }


#ifdef USE_MOMENTUM_FEEDBACK //{
    if(ReturnCalcSizeElementNumber(CS_TypeSN,true)
#ifdef USE_STELLAR_WIND //{
            +ReturnCalcSizeElementNumber(CS_TypeSW,true)
#endif // USE_STELLAR_WIND //}
            > 0){
        CalcEffectiveSurfaceAreaPrev(NExplosion,ActiveStellarFeedbackParticle);
        CalcEffectiveSurfaceAreaSum(NExplosion,ActiveStellarFeedbackParticle);
        CalcEffectiveSurfaceAreaVec(NExplosion,ActiveStellarFeedbackParticle);
#ifdef TASK_TEST_MOMENTUMFEEDBACK //{
        CalcEffectiveSurfaceAreaVecCheck(NExplosion,ActiveStellarFeedbackParticle);
#endif //TASK_TEST_MOMENTUMFEEDBACK //}
        // fflush(NULL);
        // MPI_Barrier(MPI_COMM_WORLD); assert(1==2);
    }

#endif // USE_MOMENTUM_FEEDBACK //}

#else // EVALUATE_SIZES_ALL_TOGETHER //}//{

#error This part will be discarded. 
///////////////////////////////////////////////////////////////////

    // pick up star particles which explode during this time-step.
    int NExplosion = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarBody(i)->Active){
            if(CheckEventTime(i) != NONE){
                NExplosion ++;
            }
        }
    }
    // fflush(NULL);

    int GlobalNExplosion = 0;
    MPI_Allreduce(&NExplosion,&GlobalNExplosion,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if(GlobalNExplosion == 0){
        return ;
    }
    if(MPIGetMyID()==MPI_ROOT_RANK)
        dprintlmpi(GlobalNExplosion);

    int StellarFeedBackList[NExplosion+1];
    int FeedBackType[NExplosion+1];
    int counter = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarBody(i)->Active){
            int CurrentFeedBackType = CheckEventTime(i);
            if(CurrentFeedBackType != NONE){
                StellarFeedBackList[counter] = i;
                FeedBackType[counter] = CurrentFeedBackType;
                counter ++;
            }
        }
    }

    // fflush(NULL);

    // Feedback radius.
    CalcFeedbackRadius(NExplosion,StellarFeedBackList,FeedBackType);

#error This part will be discarded. 
///////////////////////////////////////////////////////////////////

#endif // EVALUATE_SIZES_ALL_TOGETHER //}

    // Release energy and heavy elements to the surrouding ISM.
    ReleaseEnergyHeavyElements(NExplosion,StellarFeedBackList);

#ifdef USE_MOMENTUM_FEEDBACK //{
    MomentumFeedbackKick();
#endif // USE_MOMENTUM_FEEDBACK //}

    // End procedures
    StellarFeedbackEndProcedure(NExplosion,StellarFeedBackList,FeedBackType);

#ifdef USE_RADIATION_PRESSURE //{
    RadiationPressureKick();
#endif // USE_RADIATION_PRESSURE //}

    return ;
}

void StellarFeedbackGatherScatterESA(void){

    if(ReturnCalcSizeElementNumber(CS_TypeSN,true)
#ifdef USE_RADIATION_PRESSURE //{
            +ReturnCalcSizeElementNumber(CS_TypeRP,true) 
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
            +ReturnCalcSizeElementNumber(CS_TypeSW,true) 
#endif // USE_STELLAR_WIND //}
            == 0){
        return ;
    }


    int NExplosionSN = ReturnCalcSizeElementNumber(CS_TypeSN,false);
#ifdef USE_RADIATION_PRESSURE //{
    int NExplosionRP = ReturnCalcSizeElementNumber(CS_TypeRP,false);
#else // USE_RADIATION_PRESSURE //}//{
    int NExplosionRP = 0;
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
    int NExplosionSW = ReturnCalcSizeElementNumber(CS_TypeSW,false);
#else // USE_STELLAR_WIND //}//{
    int NExplosionSW = 0;
#endif // USE_STELLAR_WIND //}
    int NExplosion = NExplosionSN + NExplosionRP + NExplosionSW;

    int StellarFeedBackList[NExplosion+1];
    int FeedBackType[NExplosion+1];

    AllocateActiveStellarFeedbackParticle(NExplosion);
    CalcSizeSetSNInfo(ActiveStellarFeedbackParticle);
#ifdef USE_RADIATION_PRESSURE //{
    CalcSizeSetRPInfo(ActiveStellarFeedbackParticle+NExplosionSN);
#endif  // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
    CalcSizeSetSWInfo(ActiveStellarFeedbackParticle+NExplosionSN+NExplosionRP);
#endif  // USE_STELLAR_WIND //}
    for(int i=0;i<NExplosion;i++){
        StellarFeedBackList[i] = ActiveStellarFeedbackParticle[i].Index;
        FeedBackType[i] = ActiveStellarFeedbackParticle[i].Type;
    }


#ifdef USE_MOMENTUM_FEEDBACK //{
    if(ReturnCalcSizeElementNumber(CS_TypeSN,true)
#ifdef USE_STELLAR_WIND //{
            +ReturnCalcSizeElementNumber(CS_TypeSW,true)
#endif // USE_STELLAR_WIND //}
            > 0){

        CalcEffectiveSurfaceAreaPrev(NExplosion,ActiveStellarFeedbackParticle);
        CalcEffectiveSurfaceAreaSum(NExplosion,ActiveStellarFeedbackParticle);
        CalcEffectiveSurfaceAreaVec(NExplosion,ActiveStellarFeedbackParticle);
#if 0
        CalcEffectiveSurfaceAreaVecCheck(NExplosion,ActiveStellarFeedbackParticle);
        fflush(NULL);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(1);
#endif
    }
#endif // USE_MOMENTUM_FEEDBACK //}


    // Release energy and heavy elements to the surrouding ISM.
    ReleaseEnergyHeavyElements(NExplosion,StellarFeedBackList);

#ifdef USE_MOMENTUM_FEEDBACK //{
    MomentumFeedbackKick();
#endif // USE_MOMENTUM_FEEDBACK //}

    // End procedures
    StellarFeedbackEndProcedure(NExplosion,StellarFeedBackList,FeedBackType);

#ifdef USE_RADIATION_PRESSURE //{
    RadiationPressureKick();
#endif // USE_RADIATION_PRESSURE //}

    return ;
}

#ifdef TASK_TEST_STELLARFEEDBACK //{
void StellarFeedbackTestCalcFeedbackRadius(int NActives, const int IndexList[restrict], const int TypeList[restrict], double Radius[restrict], int Nlist[restrict], int CheckSum[restrict]){
    CalcFeedbackRadius(NActives,IndexList,TypeList);
    for(int i=0;i<NActives;i++){
        Radius[i] = ActiveStellarFeedbackParticle[i].Radius;
        Nlist[i] = ActiveStellarFeedbackParticle[i].Nlist;
#ifdef __CHECK_SUM__
        CheckSum[i] = ActiveStellarFeedbackParticle[i].CheckSum;
#endif // __CHECK_SUM__
        gprintlmpi(ActiveStellarFeedbackParticle[i].Density);
    }

    return ;
}

void StellarFeedbackTestReleaseEnergyHeavyElements(const int NExplosion, const int IndexList[restrict], const int TypeList[restrict]){
    ReleaseEnergyHeavyElements(NExplosion,IndexList);
    StellarFeedbackEndProcedure(NExplosion,IndexList,TypeList);
    return ;
}
#endif // TASK_TEST_STELLARFEEDBACK //}


#endif // USE_CELIB //}

