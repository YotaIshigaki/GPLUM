#include "config.h"
#include "PlantHydroTree.h"
#include "NeighborSearch.h"
#include "KernelFunctions.h"
#include "HIIregion.h"
#include "StellarFeedback.h"
#include "SizeDetermination.h"
#ifdef USE_RADIATION_PRESSURE //{
#include "RadiationPressure.h"
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
#include "StellarWind.h"
#endif // USE_STELLAR_WIND //}

static int Niteration_for_debug = 0;

/*! \file SizeDetermination.c
 * \brief The kernel sizes of SPH particles, the sizes of HII regions, and
 * feedback radii are determined in this routine. 
 */

#define UPDATE_SIZE_LOCAL
//#define ADD_PERTURBATION 
#define USE_DEBUG_MODE

// #define RadiusFactInc   (1.14) // 1.5 ^ (1/3)
// #define RadiusFactDec   (0.79) // 0.75 ^ (1/3)
#define RadiusFactInc   (1.2599) // 2^(1/3)
#define RadiusFactDec   (0.79372) // 0.5^(1/3)
#define RadiusFactInc_First   (3) // 1.5 ^ (1.3)
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
static double SmoothedMassConversionFactor;
#endif //USE_SMOOTHED_NEIGHBOR_NUMBER

static int MaxIterationTimes = 20;
static bool OverMaxIterationTimes = false;
static double LocalKernelMax = 0.e0;

#ifdef USE_NEIGHBOR_LIST //{
static int *LocalNeighborList;
static double *LNLK;
static int LocalNeighborListStride = 0;
#endif // USE_NEIGHBOR_LIST //}
#ifdef USE_NEIGHBOR_LIST_AND_SORT //{
static int NnbSort;
#endif // USE_NEIGHBOR_LIST_AND_SORT //}


/// For HII region
#define __EXPORT_TIMESTEP__
#define __PHOTON_COUNT_BASE__
/// For HII region

/// For SN region
#define __CHECK_SUM__
#define __CHECK_WEIGHT__
/// For SN region


////////////////////////////////////////////////////////////////////////////////////////////////

struct StructCSExport{
    int       Type;
    double    Radius;  // Radius
    double    Pos[3];  // Position.
    int       Leaf;
    bool      ExtraIterationFlag;
#ifdef USE_DEBUG_MODE
     unsigned long int GlobalID;
#endif // USE_DEBUG_MODE

    // HII region
    // double    Radius;  // Radius.
#ifdef __EXPORT_TIMESTEP__
    int       k_star;
#endif // __EXPORT_TIMESTEP__
};


#if 0
// This part can be update with union
struct StructCSImport{
    int       Type;
    double    SmoothedNumber;    // Mass.
    int       Nlist;   // Nlist.
    int       Leaf;
    bool      ExtraIterationFlag;
#ifdef USE_DEBUG_MODE
    unsigned long int GlobalID;
#endif // USE_DEBUG_MODE
    // HII region
#ifdef __PHOTON_COUNT_BASE__ //{
    double    PhotonCount;    // Total absorbed photon number per sec.
    double    PhotonCountDistanceMin;
#else // __PHOTON_COUNT_BASE__ //} //{
    double    Mass;    // Mass.
    double    MassDistanceMin;
#endif // __PHOTON_COUNT_BASE__ //}
#ifdef SET_SNII_TEMPERATURE //{
    double    GasMass;
#endif //SET_SNII_TEMPERATURE //}
    double    DistanceMin;
    double    Density;
#ifdef USE_SINK_PARTICLE //{
    double    PotentialMin;
    double    MassTotal;
    double    VCOM[3];
#endif // USE_SINK_PARTICLE //}
#ifdef USE_RADIATION_PRESSURE //{
    double    MetalMass;
#endif // USE_RADIATION_PRESSURE //{
    double    WeightSum;
#ifdef __CHECK_SUM__ //{
    int       CheckSum;
#endif //__CHECK_SUM__ //}
};
#else

struct StructCSImportHydroParticle {
    int k_hydro_localmin;    
#ifdef USE_SINK_PARTICLE //{
    double    PotentialMin;
    double    DensityMax;
    double    MassTotal;
    double    VCOM[3];
#ifdef USE_SINK_NOLOCALSINK_CONDITION //{
    bool      NoLocalSink;
#endif // USE_SINK_NOLOCALSINK_CONDITION //}
#endif // USE_SINK_PARTICLE //}
    bool      ExtraIterationFlag;
};

struct StructCSImportHIIParticle {
#ifdef __PHOTON_COUNT_BASE__ //{
    double    PhotonCount;    // Total absorbed photon number per sec.
    double    PhotonCountDistanceMin;
#else // __PHOTON_COUNT_BASE__ //} //{
    double    Mass;    // Mass.
    double    MassDistanceMin;
#endif // __PHOTON_COUNT_BASE__ //}
    double    DistanceMin;
};

struct StructCSImportSNParticle {
    double    Density;
#ifdef SET_SNII_TEMPERATURE //{
    double    GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
    double    DistanceMin;
    unsigned long int  DistanceMinGlobalID;
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef USE_MOMENTUM_FEEDBACK //{  
    double nH_ave;
    double Z_ave;
    double NumberDensity;
#endif // USE_MOMENTUM_FEEDBACK //}
    double    WeightSum;
#ifdef __CHECK_SUM__ //{
    int       CheckSum;
#endif //__CHECK_SUM__ //}
};

#ifdef USE_RADIATION_PRESSURE //{
struct StructCSImportRPParticle {
    double    MetalMass;
    double    GasMass;
    double    WeightSum;
};
#endif // USE_RADIATION_PRESSURE //{

#ifdef USE_STELLAR_WIND //{
struct StructCSImportSWParticle {
#ifdef MAXIMUM_ENERGY_INPUT //{
    double    DistanceMin;
    unsigned long int  DistanceMinGlobalID;
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef USE_MOMENTUM_FEEDBACK //{  
    double nH_ave;
    double Z_ave;
    double NumberDensity;
#endif // USE_MOMENTUM_FEEDBACK //}
    double    WeightSum;
};
#endif // USE_STELLAR_WIND //}

struct StructCSImport{
    int       Type;
    int       k_hydro_localmin;
    double    SmoothedNumber;    // Mass.
    int       Nlist;   // Nlist.
    int       Leaf;
#ifdef USE_DEBUG_MODE
    unsigned long int GlobalID;
#endif // USE_DEBUG_MODE
    union {
        struct StructCSImportHydroParticle Hydro; 
        struct StructCSImportHIIParticle   HII; 
        struct StructCSImportSNParticle    SN; 
#ifdef USE_RADIATION_PRESSURE //{
        struct StructCSImportRPParticle    RP; 
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
        struct StructCSImportSWParticle    SW; 
#endif // USE_STELLAR_WIND //}
    } Body;
};

#endif


static int CS_NContactedDomains;
static int *CS_ContactedDomainID;

static inline void CS_AllocateContactedDomainID(void){
    CS_ContactedDomainID = malloc(sizeof(int)*MPIGetNumProcs());
    return ;
}

static inline bool CS_CheckLocalExternalDomainsContacted(const int MyDomainID, const int ExtDomainID){

    for(int k=0;k<3;k++){
        if((EdgesForHydro[MyDomainID].PosMax[k] < EdgesForHydro[ExtDomainID].PosMin[k])||
           (EdgesForHydro[MyDomainID].PosMin[k] > EdgesForHydro[ExtDomainID].PosMax[k]))  return false;
    }
    return true;

}

/*
 * This function checkes that how many contacted domains around the local
 * domain. Comparing the local domain edge to the external domain, 
 */
static inline void CS_CheckContactedDomain(void){
    const int NProcs = MPIGetNumProcs();
    const int MyID = MPIGetMyID();
    CS_NContactedDomains = 0;
    for(int i=0;i<NProcs-1;i++){
        int DomainID = CommunicationTable[i].SendRank;
        if(CS_CheckLocalExternalDomainsContacted(MyID,DomainID)){
            CS_ContactedDomainID[CS_NContactedDomains] = i;
            CS_NContactedDomains ++;
        }
    }
    return ;
}


struct StructCSActiveHydroParticle{
    double Rho;
#ifdef USE_SINK_PARTICLE //{
    double PotentialMin;
    double DensityMax;
    double MassTotal;
    double VCOM[3];
#ifdef USE_SINK_NOLOCALSINK_CONDITION //{
    bool   NoLocalSink;
#endif // USE_SINK_NOLOCALSINK_CONDITION //}
#endif // USE_SINK_PARTICLE //}
    int CacheIndex;
    int ExtraIteration;
}; 

struct StructCSActiveHIIParticle{
#ifdef __PHOTON_COUNT_BASE__ //{
    double PhotonCount;
#else //__PHOTON_COUNT_BASE__ //} //{
    double Mass;
#endif // __PHOTON_COUNT_BASE__ //}
    double LyAphoton;  // The number of Lyman continum photons [s^{-1}].
    // double Radius;

    double DistanceMin;
#ifdef __PHOTON_COUNT_BASE__ //{
    double PhotonCountDistanceMin;
#else //__PHOTON_COUNT_BASE__ //} //{
    double MassDistanceMin;
#endif // __PHOTON_COUNT_BASE__ //}
    bool HIIRegion;
};

struct StructCSActiveSNParticle{
    double Density;    // The mass weighted nomalization factor.
#ifdef SET_SNII_TEMPERATURE
    double GasMass;// Local gas mass.
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
    double DistanceMin;
    unsigned long int DistanceMinGlobalID;
#endif
#ifdef USE_MOMENTUM_FEEDBACK //{
    double nH_ave;
    double Z_ave;
    double NumberDensity;
#endif // USE_MOMENTUM_FEEDBACK //}
#ifdef __CHECK_SUM__ //{
    int CheckSum;
#endif // __CHECK_SUM__ //}
    int Type; // 1=TypeII, 2=TypeIa, 3=AGB
    int Count;
    double InitialMass;
    double Metallicity;
    double WeightSum;
    int IterationCount;
};

#ifdef USE_RADIATION_PRESSURE //{
struct StructCSActiveRPParticle{
    double GasMass;     // Total gas mass
    double MetalMass;   // Total metal mass sum_i m_i*Z_i
    double WeightSum;   // sum_i wi;
    int IterationCount; // can be removed 
};
#endif // USE_RADIATION_PRESSURE //}

struct StructCSActiveSWParticle{
#ifdef MAXIMUM_ENERGY_INPUT
    double DistanceMin;
    unsigned long int DistanceMinGlobalID;
#endif
#ifdef USE_MOMENTUM_FEEDBACK //{
    double nH_ave;
    double Z_ave;
    double NumberDensity;
#endif // USE_MOMENTUM_FEEDBACK //}
    double InitialMass;
    double Metallicity;
    double WeightSum;
};

static int CSExportFlagsMaxAllocated = 0;
struct StructActiveParticle{
    int Type;
    int Index; // NBCache[Index].
    int Nlist; // Number of neighbors.
    bool LocalUpdateFlags; // Flag for the local update.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double SmoothedNumber; // Smoothed mas..
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
    double Pos[3];
    double Radius; // 2xRadius is the support scale.
    double Rvalue;
    double Lvalue;
    int k_hydro_localmin;
    int k;
    int Counter;
#ifdef USE_DEBUG_MODE //{
    unsigned long long int GlobalID;
#endif // USE_DEBUG_MODE //}
    union {
        struct StructCSActiveHydroParticle Hydro; 
        struct StructCSActiveHIIParticle   HII; 
        struct StructCSActiveSNParticle    SN; 
#ifdef USE_RADIATION_PRESSURE //{
        struct StructCSActiveRPParticle    RP; 
#endif // USE_RADIATION_PRESSURE //{
#ifdef USE_STELLAR_WIND //{
        struct StructCSActiveSWParticle    SW; 
#endif // USE_STELLAR_WIND //{
    } Body;
} *ActiveParticle;

static int HydroEntry(struct StructActiveParticle AP[], bool CalcSizeExportFlags[], const int NProcs){
    
    int NActives = 0;
#if 1
    int RootNodeID = 0; 
    int NumberofLeaves = HydroNode[RootNodeID].NumberofLeaves;
    int header = HydroNode[RootNodeID].Leaves;
    for(int i=0;i<NumberofLeaves;i++){
        int leaf = header + i; 
        if(HydroRoot.Leaves[leaf] < 0) continue;
        if(NBCache[leaf].Active){
            int Index = NBCache[leaf].Leaf;
            CalcSizeExportFlags[NActives*NProcs+NProcs-1] = true;
            AP[NActives].Index = Index;
            AP[NActives].Nlist = Phydro[Index]->Nlist;
            AP[NActives].Type = CS_TypeHydro;
            AP[NActives].Radius = NBCache[leaf].Kernel;
            /*
            if(NBCache[leaf].Kernel == 0){
                fprintf(stderr,"Special alart! Kernel size = 0\n");
                dlprintlmpi(PhydroBody(Index)->GlobalID);
                dprintlmpi(leaf);
                gprintlmpi(AP[NActives].Radius);
                gprintlmpi(Phydro[Index]->Kernel);
                gprintlmpi(Phydro[Index]->KernelPred);
                fflush(NULL);
            }
            */
            AP[NActives].Rvalue = AP[NActives].Lvalue = 0.e0;
            AP[NActives].Pos[0] = NBCache[leaf].Pos[0];
            AP[NActives].Pos[1] = NBCache[leaf].Pos[1];
            AP[NActives].Pos[2] = NBCache[leaf].Pos[2];
            AP[NActives].k = Phydro[i]->k_hydro;
            // AP[NActives].Pos[0] = PhydroBody(Index)->PosP[0];
            // AP[NActives].Pos[1] = PhydroBody(Index)->PosP[1];
            // AP[NActives].Pos[2] = PhydroBody(Index)->PosP[2];
            AP[NActives].Body.Hydro.CacheIndex = leaf;
            AP[NActives].Counter = 0;
#ifdef USE_DEBUG_MODE //{
            AP[NActives].GlobalID = PhydroBody(Index)->GlobalID;
#endif // USE_DEBUG_MODE //}
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
            AP[NActives].Body.Hydro.SmoothedNumber = Phydro[NBCache[leaf].Leaf]->SmoothedNumber;
#endif //ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
            NActives ++;
        }
    }
#else 
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            CalcSizeExportFlags[NActives*NProcs+NProcs-1] = true;
            AP[NActives].Index = i;
            AP[NActives].Nlist = Phydro[i]->Nlist;
            AP[NActives].Type = CS_TypeHydro;
            AP[NActives].Radius = Phydro[i]->Kernel;
            AP[NActives].Rvalue = AP[NActives].Lvalue = 0.e0;
            AP[NActives].Pos[0] = PhydroBody(i)->Pos[0];
            AP[NActives].Pos[1] = PhydroBody(i)->Pos[1];
            AP[NActives].Pos[2] = PhydroBody(i)->Pos[2];
            AP[NActives].Body.Hydro.CacheIndex = i;
#error

            NActives ++;
        }
    }
#endif
    return NActives;
}

static int HIIregionEntry(struct StructActiveParticle AP[], bool CalcSizeExportFlags[], const int NProcs){

    int NActiveYoungStars = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(Pstar[i]->HIIflag){
            int Offset = NActiveYoungStars*NProcs;
            CalcSizeExportFlags[Offset+NProcs-1] = true;

            AP[NActiveYoungStars].Index = i;
            AP[NActiveYoungStars].Nlist = 0;
            AP[NActiveYoungStars].Type = CS_TypeHII;
            AP[NActiveYoungStars].Rvalue = AP[NActiveYoungStars].Lvalue = 0.e0;
            AP[NActiveYoungStars].Pos[0] = PstarBody(i)->PosP[0];
            AP[NActiveYoungStars].Pos[1] = PstarBody(i)->PosP[1];
            AP[NActiveYoungStars].Pos[2] = PstarBody(i)->PosP[2];
            AP[NActiveYoungStars].k = PstarBody(i)->k;
            // Set the initial Radius.
            if((Pstar[i]->StromgrenRadius>0.e0)&&(Pstar[i]->StromgrenRadius*Pall.UnitLength < 50*PC_CGS)){ // Reuse old one.
                AP[NActiveYoungStars].Radius = Pstar[i]->StromgrenRadius;
            } else{ // Set an initial guess.
                AP[NActiveYoungStars].Radius = 0.25*PstarBody(i)->Eps;
            }

            // fprintf(stderr,"-- %d %d %g %d:%s\n",NActiveYoungStars,i,AP[NActiveYoungStars].Radius,__LINE__,__FUNCTION__);
            // fflush(NULL);

#ifdef __PHOTON_COUNT_BASE__
            AP[NActiveYoungStars].Body.HII.PhotonCount = 0.e0;
#else //__PHOTON_COUNT_BASE__
            AP[NActiveYoungStars].Body.HII.Mass = 0.e0;
#endif // __PHOTON_COUNT_BASE__

#ifdef PRESERVE_SNII_EVENTRATE
            double MassInMsun = Pstar[i]->InitialMass*Pall.UnitMass/MSUN_CGS;
            double prob = SNIINumber*MassInMsun;
            if(prob >= 1.0){
#   ifdef USE_ASRFLX_NLY // {
                double Age = Pall.UnitTime*(Pall.TCurrent - Pstar[i]->FormationTime)/YEAR_CGS;
                AP[NActiveYoungStars].Body.HII.LyAphoton = (PstarBody(i)->Mass*Pall.UnitMass/MSUN_CGS)*
                    ASRFLXGetNLY(Age,Pstar[i]->Z);
#   else //USE_ASRFLX_NLY //}//{
                AP[NActiveYoungStars].Body.HII.LyAphoton = (PstarBody(i)->Mass*Pall.UnitMass/MSUN_CGS)*
                    ReturnNumberofLyPhoton(i,0);
#   endif // USE_ASRFLX_NLY // }
            } else {
                if(Pstar[i]->TypeIIProb){
#   ifdef USE_ASRFLX_NLY // {
                    double Age = Pall.UnitTime*(Pall.TCurrent - Pstar[i]->FormationTime)/YEAR_CGS;
                    AP[NActiveYoungStars].Body.HII.LyAphoton = (1.0/SNIINumber)*
                        ASRFLXGetNLY(Age,Pstar[i]->Z);
#   else //USE_ASRFLX_NLY //}//{
                    AP[NActiveYoungStars].Body.HII.LyAphoton = (1.0/SNIINumber)*ReturnNumberofLyPhoton(i,0);
#   endif // USE_ASRFLX_NLY // }
                } 
            }
#else
#   ifdef USE_ASRFLX_NLY // {
            double Age = Pall.UnitTime*(Pall.TCurrent - Pstar[i]->FormationTime)/YEAR_CGS;
            AP[NActiveYoungStars].Body.HII.LyAphoton = (PstarBody(i)->Mass*Pall.UnitMass/MSUN_CGS)*
                ASRFLXGetNLY(Age,Pstar[i]->Z);
                //ASRFLXGetNLY(Age,0);
            // fprintf(stderr,"HII %g %g\n",AP[NActiveYoungStars].Body.HII.LyAphoton,
                    // (PstarBody(i)->Mass*Pall.UnitMass/MSUN_CGS)*ReturnNumberofLyPhoton(i,0));
#   else //USE_ASRFLX_NLY //}//{
            AP[NActiveYoungStars].Body.HII.LyAphoton = (PstarBody(i)->Mass*Pall.UnitMass/MSUN_CGS)*
                ReturnNumberofLyPhoton(i,0);
#   endif // USE_ASRFLX_NLY // }
#endif
            AP[NActiveYoungStars].Body.HII.HIIRegion = true;

            AP[NActiveYoungStars].Counter = 0;
#ifdef USE_DEBUG_MODE //{
            AP[NActiveYoungStars].GlobalID = PstarBody(i)->GlobalID;
#endif // USE_DEBUG_MODE //}

            NActiveYoungStars ++;
        }
    }
    return NActiveYoungStars;
}

static inline int __attribute__((always_inline)) CheckEventTime(const int Index){

#ifdef USE_CELIB //{
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
#endif // USE_CELIB_NSM //}


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
#endif // USE_CELIB //}
}

static int StellarFeedbackEntry(struct StructActiveParticle AP[], bool CalcSizeExportFlags[], const int NProcs){

    int NExplosoin = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarBody(i)->Active){
            int CurrentFeedbackType = CheckEventTime(i);
            if(CurrentFeedbackType != NONE){
                int Offset = NExplosoin*NProcs;
                CalcSizeExportFlags[Offset+NProcs-1] = true;
                AP[NExplosoin].Index = i;
                AP[NExplosoin].Nlist = 0;
                AP[NExplosoin].Type = CS_TypeSN;
                AP[NExplosoin].Rvalue = AP[NExplosoin].Lvalue = 0.e0;
                AP[NExplosoin].Pos[0] = PstarPos(i)[0];
                AP[NExplosoin].Pos[1] = PstarPos(i)[1]; 
                AP[NExplosoin].Pos[2] = PstarPos(i)[2]; 
                AP[NExplosoin].k = PstarBody(i)->k;
#ifdef CELib //}
                AP[NExplosoin].Radius = StellarFeedbackRadiusInitialGuess(i);
#else // CELib //}//{
                AP[NExplosoin].Radius = 2.0*PstarBody(i)->Eps; 
#endif // CELib //}
                AP[NExplosoin].Body.SN.Type = CurrentFeedbackType;
#ifdef USE_CELIB //{
                if(CurrentFeedbackType == CELibFeedbackType_AGB){
                    AP[NExplosoin].Body.SN.Count = Pstar[i]->AGBCount;
#ifdef USE_CELIB_NSM //{
                } else if(CurrentFeedbackType == CELibFeedbackType_NSM){
                    AP[NExplosoin].Body.SN.Count = Pstar[i]->NSMCount;
#endif // USE_CELIB_NSM //}
                } else {                        
                    AP[NExplosoin].Body.SN.Count = Pstar[i]->SNIaCount;
                }
#endif // USE_CELIB //}
                AP[NExplosoin].Body.SN.InitialMass = Pstar[i]->InitialMass;
                AP[NExplosoin].Body.SN.Metallicity = Pstar[i]->Z;
                AP[NExplosoin].Body.SN.IterationCount = 0;
                AP[NExplosoin].Body.SN.Type = CurrentFeedbackType;
                AP[NExplosoin].Counter = 0;
#ifdef USE_DEBUG_MODE //{
                AP[NExplosoin].GlobalID = PstarBody(i)->GlobalID;
#endif // USE_DEBUG_MODE //}
                NExplosoin ++;
            }
        }
    }

    return NExplosoin;
}

#ifdef USE_RADIATION_PRESSURE //{

static int MaxIterationTimesForInitialGuessRP = 10;
double RadiationPressureRadiusInitialGuess(const int Index){

    int Neighbors[MaxNeighborSize];
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
    if(Pstar[Index]->RadiusRP > 0)
        return Pstar[Index]->RadiusRP;

    // if(Pall.Nhydro == 0)
        //return PstarBody(Index)->Eps;
    //double Radius = Pstar[Index]->RadiusRP; //return PstarBody(Index)->Eps;
    double Radius = PstarBody(Index)->Eps;

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
        Radius *= RadiusFactInc; 
        Iteration ++;
    } while(Iteration < MaxIterationTimesForInitialGuessRP);
    return 2.0*Radius;
}

static int RadiationPressureEntry(struct StructActiveParticle AP[], bool CalcSizeExportFlags[], const int NProcs){

    int NRadiationPressure = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarBody(i)->Active){
            //if(Pall.TCurrent-Pstar[i]->FormationTime < GetRadiationPressureSourceLifeTime()){
            if(CheckRadiationPressureSource(i) == true){
                int Offset = NRadiationPressure*NProcs;
                CalcSizeExportFlags[Offset+NProcs-1] = true;
                AP[NRadiationPressure].Index = i;
                AP[NRadiationPressure].Nlist = 0;
                AP[NRadiationPressure].Type = CS_TypeRP;
                AP[NRadiationPressure].Rvalue = AP[NExplosoin].Lvalue = 0.e0;
                AP[NRadiationPressure].Pos[0] = PstarPos(i)[0];
                AP[NRadiationPressure].Pos[1] = PstarPos(i)[1]; 
                AP[NRadiationPressure].Pos[2] = PstarPos(i)[2]; 
                AP[NRadiationPressure].k = PstarBody(i)->k;

                AP[NRadiationPressure].Radius = RadiationPressureRadiusInitialGuess(i);
                AP[NRadiationPressure].Counter = 0;
#ifdef USE_DEBUG_MODE //{
                AP[NRadiationPressure].GlobalID = PstarBody(i)->GlobalID;
#endif // USE_DEBUG_MODE //}
                NRadiationPressure ++;
            }
        }
    }

    return NRadiationPressure;
}

#endif // USE_RADIATION_PRESSURE //}

#ifdef USE_STELLAR_WIND //{

static int MaxIterationTimesForInitialGuessSW = 10;
double StellarWindRadiusInitialGuess(const int Index){

    int Neighbors[MaxNeighborSize];
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
    if(Pstar[Index]->RadiusSW > 0)
        return Pstar[Index]->RadiusSW;

    double Radius = PstarBody(Index)->Eps;

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
        Radius *= RadiusFactInc; 
        Iteration ++;
    } while(Iteration < MaxIterationTimesForInitialGuessSW);
    return 2.0*Radius;
}

static int StellarWindEntry(struct StructActiveParticle AP[], bool CalcSizeExportFlags[], const int NProcs){

    int NStellarWind = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarBody(i)->Active){
            if(CheckStellarWindSource(i) == true){
                int Offset = NStellarWind*NProcs;
                CalcSizeExportFlags[Offset+NProcs-1] = true;
                AP[NStellarWind].Index = i;
                AP[NStellarWind].Nlist = 0;
                AP[NStellarWind].Type = CS_TypeSW;
                AP[NStellarWind].Rvalue = AP[NStellarWind].Lvalue = 0.e0;
                AP[NStellarWind].Pos[0] = PstarPos(i)[0];
                AP[NStellarWind].Pos[1] = PstarPos(i)[1]; 
                AP[NStellarWind].Pos[2] = PstarPos(i)[2]; 
                AP[NStellarWind].k = PstarBody(i)->k;
                AP[NStellarWind].Radius = StellarWindRadiusInitialGuess(i);
                AP[NStellarWind].Counter = 0;
#ifdef USE_DEBUG_MODE //{
                AP[NStellarWind].GlobalID = PstarBody(i)->GlobalID;
#endif // USE_DEBUG_MODE //}
                NStellarWind ++;
            }
        }
    }

    return NStellarWind;
}

#endif // USE_STELLAR_WIND //}

static inline bool __attribute__((always_inline)) CS_OverlapDomain(double Pos[restrict], const double h, const int NodeID){ 

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForHydro[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

/*
 * This function return number of particles which should export to the node ID of [Index].
 */
static inline int __attribute__((always_inline)) CS_CheckExportFlags(const int DomainID, const int NProcs, bool CSExportFlags[restrict], const int NActives, struct StructActiveParticle AP[restrict]){

    if((Pall.Nhydro+Pall.Nstars) == 0)
        return 0;

    int ExportDomainID = CommunicationTable[DomainID].SendRank;

    double BoxCenter[] = {HydroNode[0].Pos[0],HydroNode[0].Pos[1],HydroNode[0].Pos[2]}; // need check
    if(!CS_OverlapDomain(BoxCenter,LocalKernelMax,ExportDomainID)){
        return 0;
    }

    int NExport = 0;
    for(int i=0;i<NActives;i++){
        int Offset = i*NProcs;
        if(CSExportFlags[Offset+NProcs-1]){
            if(CS_OverlapDomain(AP[i].Pos,2.0*AP[i].Radius,ExportDomainID)){
                CSExportFlags[Offset+DomainID] = true;
                NExport ++;
            }
        }
    }

    return NExport;
}

struct StructCSHydroLocalInfo {
    int Nlist;
    int k_hydro_localmin;
#ifdef USE_SINK_PARTICLE //{
    double PotentialMin;
    double DensityMax;
    double MassTotal;
    double VCOM[3];
#ifdef USE_SINK_NOLOCALSINK_CONDITION //{
    bool   NoLocalSink;
#endif // USE_SINK_NOLOCALSINK_CONDITION //}
#endif // USE_SINK_PARTICLE //}
};

#ifdef USE_SINK_NOLOCALSINK_CONDITION //{
static bool CheckNoLocalSinkCondition(const double h, const double Pos[restrict]){

    //double h2 = SQ(2*h);

    for(int i=0;i<Pall.Nsink;i++){
        double h2 = SQ(2*h + Psink[i]->AccretionRadius);
        if(h2 > DISTANCE2(Pos,PsinkBody(i)->PosP)){
            return false;
        }
    }

    return true;
}
#endif // USE_SINK_NOLOCALSINK_CONDITION //}

static struct StructCSHydroLocalInfo CS_GetNumberofNeighbors(double Pos[restrict], const double Kernel, int Neighbors[restrict]){

    struct StructCSHydroLocalInfo TempHydroLocalInfo = {
        .Nlist = 0,
        .k_hydro_localmin = MAXIMUM_TIME_HIERARCHY,
#ifdef USE_SINK_PARTICLE //{
        .PotentialMin = 0.e0,
        .DensityMax = 0.e0,
        .MassTotal = 0.e0,
        .VCOM[0] = 0.e0,
        .VCOM[1] = 0.e0,
        .VCOM[2] = 0.e0,
#ifdef USE_SINK_NOLOCALSINK_CONDITION //{
        .NoLocalSink = true,
#endif // USE_SINK_NOLOCALSINK_CONDITION //}
#endif // USE_SINK_PARTICLE //}
    };

    int Iteration = 0;
    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    do {
        CurrentNodeID = GetNeighborsIterativeApproach(CurrentNodeID,Pos,2.e0*Kernel,&(TempHydroLocalInfo.Nlist),Neighbors);

#ifdef USE_SINK_PARTICLE //{
        if((TempHydroLocalInfo.Nlist>0)&&(Iteration==0)){
#ifdef USE_SMOOTHED_POTENTIAL //{
            TempHydroLocalInfo.PotentialMin = Phydro[Neighbors[0]]->Pot;
#else // USE_SMOOTHED_POTENTIAL //}//{
            TempHydroLocalInfo.PotentialMin = PhydroBody(Neighbors[0])->Pot/PhydroBody(Neighbors[0])->Mass;
#endif // USE_SMOOTHED_POTENTIAL //}
            TempHydroLocalInfo.DensityMax = Phydro[Neighbors[0]]->RhoPred;
        }

        for(int i=0;i<TempHydroLocalInfo.Nlist;i++){
            int index = Neighbors[i];

            TempHydroLocalInfo.k_hydro_localmin = MIN(TempHydroLocalInfo.k_hydro_localmin,Phydro[index]->k_hydro);

            TempHydroLocalInfo.MassTotal += PhydroBody(index)->Mass;
            TempHydroLocalInfo.VCOM[0] += PhydroBody(index)->Mass*PhydroBody(index)->Vel[0];
            TempHydroLocalInfo.VCOM[1] += PhydroBody(index)->Mass*PhydroBody(index)->Vel[1];
            TempHydroLocalInfo.VCOM[2] += PhydroBody(index)->Mass*PhydroBody(index)->Vel[2];

#ifdef USE_SMOOTHED_POTENTIAL //{
            TempHydroLocalInfo.PotentialMin = fmin(TempHydroLocalInfo.PotentialMin,Phydro[index]->Pot);
#else // USE_SMOOTHED_POTENTIAL //}//{
            TempHydroLocalInfo.PotentialMin = fmin(TempHydroLocalInfo.PotentialMin,PhydroBody(index)->Pot/PhydroBody(index)->Mass);
#endif // USE_SMOOTHED_POTENTIAL //}
            TempHydroLocalInfo.DensityMax = fmax(TempHydroLocalInfo.DensityMax,Phydro[index]->RhoPred);
        }
#ifdef USE_SINK_NOLOCALSINK_CONDITION //{
        if(!CheckNoLocalSinkCondition(Kernel,Pos)){
            TempHydroLocalInfo.NoLocalSink = false;
        }
#endif // USE_SINK_NOLOCALSINK_CONDITION //}
#endif // USE_SINK_PARTICLE //}
        Iteration ++;
    }while(CurrentNodeID != RootNodeID);

    return TempHydroLocalInfo;
}

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
static int CS_GetSmoothedNumberofNeighbors(double Pos[restrict], const double Kernel, int Neighbors[restrict], double *SmoothedNumber){

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

struct StructCSHIILocalInfo {
    int Nlist;
    int k_hydro_localmin;
#ifdef __PHOTON_COUNT_BASE__ //{
    double PhotonCount;
    double PhotonCountDistanceMin;
#else //__PHOTON_COUNT_BASE__ //} //{
    double Mass;
    double MassDistanceMin;
#endif // __PHOTON_COUNT_BASE__ //}
    double DistanceMin;
};

/*
 * This function counts (1) the total number of photons which can be absorbed in
 * the neighbor particles, or (2) the total mass of the neighbor particles.
 */
static struct StructCSHIILocalInfo CS_ReturnHIILocalInfo(double Pos[restrict], const double Radius, int Neighbors[restrict]){

    struct StructCSHIILocalInfo TempHIILocalInfo = {
        .Nlist = 0,
        .k_hydro_localmin = MAXIMUM_TIME_HIERARCHY,
        .DistanceMin = 0.0,
        //.DistMin = 10000,
#ifdef __PHOTON_COUNT_BASE__ //{
        .PhotonCount = 0.e0,
        .PhotonCountDistanceMin = 0.0,
#else //__PHOTON_COUNT_BASE__ //} //{
        .Mass = 0.e0,
        .MassDistanceMin = 0.0,
#endif // __PHOTON_COUNT_BASE__ //}
    };

    int Iteration = 0;
    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    do {
        CurrentNodeID = GetNeighborsIterativeApproach(CurrentNodeID,Pos,
                2.e0*Radius,&(TempHIILocalInfo.Nlist),Neighbors);


        if((TempHIILocalInfo.Nlist>0)&&(Iteration==0)){
            TempHIILocalInfo.DistanceMin = DISTANCE2(PhydroBody(Neighbors[0])->PosP,Pos);
        }

        int Nlist = 0;
        for(int i=0;i<TempHIILocalInfo.Nlist;i++){
            int index = Neighbors[i];
            if(Phydro[index]->UPred*Pall.ConvertUtoT > 1.e+4) continue;
            if(Phydro[index]->RhoPred*Pall.ConvertNumberDensityToCGS < 1.e-1) continue;
            TempHIILocalInfo.k_hydro_localmin = fmin(TempHIILocalInfo.k_hydro_localmin,Phydro[index]->k_hydro);
            Nlist ++;

#ifdef __PHOTON_COUNT_BASE__ //{ // fiducial model.
            const double InvProtonMass = 1.0/PROTON_MASS_CGS;
            double X = (1.0-HeliumAbandance-Phydro[index]->Z);
            double n_H = (X*Phydro[index]->RhoPred*Pall.UnitMass*InvProtonMass)/CUBE(Pall.UnitLength);
            double n_e = n_H;
            double r3 = 3*PhydroBody(index)->Mass/(4*M_PI*Phydro[index]->RhoPred);
            double CurrentPhotonCount = (4.0*M_PI/3.0)*n_H*n_e*ReturnRecombinationConstantAlpha()*r3*CUBE(Pall.UnitLength);

            TempHIILocalInfo.PhotonCount += CurrentPhotonCount;
#else //__PHOTON_COUNT_BASE__ //} //{
            double CurrentMass = (1.0-HeliumAbandance-Phydro[index]->Z)*PhydroBody(index)->Mass;
            TempHIILocalInfo.Mass += CurrentMass;
#endif // __PHOTON_COUNT_BASE__ //}

            // Evaluate minimum distance.
            double TmpDistanceMin = DISTANCE2(PhydroBody(index)->PosP,Pos);
            if(TempHIILocalInfo.DistanceMin >= TmpDistanceMin){
                TempHIILocalInfo.DistanceMin = TmpDistanceMin;
#ifdef __PHOTON_COUNT_BASE__ //{
                TempHIILocalInfo.PhotonCountDistanceMin = CurrentPhotonCount;
                // gprintlmpi(CurrentPhotonCount);
#else //__PHOTON_COUNT_BASE__ //} //{
                TempHIILocalInfo.MassDistanceMin = CurrentMass;
#endif // __PHOTON_COUNT_BASE__ //}
            }
        }
        TempHIILocalInfo.Nlist = Nlist;
        Iteration ++;
    }while(CurrentNodeID != RootNodeID);

    return TempHIILocalInfo;
}

struct StructCSStellarFeedbackLocalInfo{
    int Nlist;
    int k_hydro_localmin;
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
#ifdef USE_MOMENTUM_FEEDBACK //{
    double nH_ave;
    double Z_ave;
    double NumberDensity;
#endif // USE_MOMENTUM_FEEDBACK //}
    double WeightSum;
#ifdef __CHECK_SUM__ //{
    int CheckSum;
#endif // __CHECK_SUM__ //}
}; 


/*
 * This function returns a structure which involves various values (the neighbor
 * number, the smootehd neighbor number, the density, the local gas mass, the
 * minimum distance, and the globalID of the closest particle).
 */
struct StructCSStellarFeedbackLocalInfo CS_RetrunStellarFeedbackLocalInfo(double Pos[restrict], const double Radius, int Neighbors[restrict]){

    //static int Neighbors[MaxNeighborSize];
    struct StructCSStellarFeedbackLocalInfo TempStellarFeedbackLocalInfo = {
        .Nlist = 0,
        .k_hydro_localmin = MAXIMUM_TIME_HIERARCHY,
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
#ifdef USE_MOMENTUM_FEEDBACK //{
        .nH_ave = 0.e0,
        .Z_ave = 0.e0,
        .NumberDensity = 0.e0,
#endif // USE_MOMENTUM_FEEDBACK //}
        .WeightSum = 0.e0,
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

        TempStellarFeedbackLocalInfo.k_hydro_localmin = MIN(TempStellarFeedbackLocalInfo.k_hydro_localmin,Phydro[leaf]->k_hydro);

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
        //double w = KernelStellarFeedback(r,InvRadiusi);
        double w = SPHKernel(r,InvRadiusi);

        TempStellarFeedbackLocalInfo.WeightSum += w;
        TempStellarFeedbackLocalInfo.Density += Phydro[leaf]->Mass*w;
        //assert(PhydroMass(leaf)*w > 0.e0);
#ifdef SET_SNII_TEMPERATURE
        TempStellarFeedbackLocalInfo.GasMass += Phydro[leaf]->Mass;
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
#ifdef USE_MOMENTUM_FEEDBACK //{
        TempStellarFeedbackLocalInfo.nH_ave += Pall.ConvertNumberDensityToCGS*Phydro[leaf]->RhoPred;
        TempStellarFeedbackLocalInfo.Z_ave += Phydro[leaf]->Z;
        TempStellarFeedbackLocalInfo.NumberDensity += w;
#endif // USE_MOMENTUM_FEEDBACK //}
#ifdef __CHECK_SUM__ //{
        TempStellarFeedbackLocalInfo.CheckSum += PhydroBody(leaf)->GlobalID;
#endif // __CHECK_SUM__ //}
    }

    return TempStellarFeedbackLocalInfo;
}

#ifdef USE_RADIATION_PRESSURE //{

struct StructCSRadiationPressureLocalInfo{
    int Nlist;
    int k_hydro_localmin;
    double GasMass;// Local gas mass.
    double MetalMass;// Local metal mass.
    double WeightSum;
}; 


/*
 * This function returns a structure which involves various values (the neighbor
 * number, the local gas and metal mass, and the sum of kernel).
 */
struct StructCSRadiationPressureLocalInfo CS_RetrunRadiationPressureLocalInfo(double Pos[restrict], const double Radius, int Neighbors[restrict]){

    struct StructCSRadiationPressureLocalInfo TempRadiationPressureLocalInfo = {
        .Nlist = 0,
        .k_hydro_localmin = MAXIMUM_TIME_HIERARCHY,
        .GasMass = 0.e0,
        .MetalMass = 0.e0,
        .WeightSum = 0.e0,
    };

    double InvRadiusi = 1.e0/Radius;
    TempRadiationPressureLocalInfo.Nlist = 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            GetSmoothedNumberofNeighbors(Pos,Radius,Neighbors,&(TempStellarFeedbackLocalInfo.SmoothedNumber));
#else // USE_SMOOTHED_NEIGHBOR_NUMBER 
            GetNeighborsLimited(Pos,2.0*Radius,Neighbors);
            //GetNeighborsDirect(Pos,2*Radius,Neighbors);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}

    for(int i=0;i<TempRadiationPressureLocalInfo.Nlist;i++){
        int leaf = Neighbors[i];

        TempRadiationPressureLocalInfo.k_hydro_localmin = MIN(TempRadiationPressureLocalInfo.k_hydro_localmin,Phydro[leaf]->k_hydro);

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
        double w = SPHKernel(r,InvRadiusi);

        TempRadiationPressureLocalInfo.WeightSum += w;
        TempRadiationPressureLocalInfo.GasMass   += Phydro[leaf]->Mass;
        TempRadiationPressureLocalInfo.MetalMass += Phydro[leaf]->Mass*Phydro[leaf]->Z;
    }

    return TempRadiationPressureLocalInfo;
}

#endif // USE_RADIATION_PRESSURE //}

#ifdef USE_STELLAR_WIND //{
struct StructCSStellarWindLocalInfo{
    int Nlist;
    int k_hydro_localmin;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
    double SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER
#ifdef MAXIMUM_ENERGY_INPUT
    double DistanceMin;
    unsigned long int DistanceMinGlobalID;
#endif
#ifdef USE_MOMENTUM_FEEDBACK //{
    double nH_ave;
    double Z_ave;
    double NumberDensity;
#endif // USE_MOMENTUM_FEEDBACK //}
    double WeightSum;
}; 


/*
 * This function returns a structure which involves various values.
 */
struct StructCSStellarWindLocalInfo CS_RetrunStellarWindLocalInfo(double Pos[restrict], const double Radius, int Neighbors[restrict]){

    struct StructCSStellarWindLocalInfo TempStellarWindLocalInfo = {
        .Nlist = 0,
        .k_hydro_localmin = MAXIMUM_TIME_HIERARCHY,
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        .SmoothedNumber = 0.0,
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef MAXIMUM_ENERGY_INPUT
        .DistanceMin = 0.e0,
        .DistanceMinGlobalID = NONE,
#endif
#ifdef USE_MOMENTUM_FEEDBACK //{
        .nH_ave = 0.e0,
        .Z_ave = 0.e0,
        .NumberDensity = 0.e0,
#endif // USE_MOMENTUM_FEEDBACK //}
        .WeightSum = 0.e0,
    };

    double InvRadiusi = 1.e0/Radius;
    TempStellarWindLocalInfo.Nlist = 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            GetSmoothedNumberofNeighbors(Pos,Radius,Neighbors,&(TempStellarWindLocalInfo.SmoothedNumber));
#else // USE_SMOOTHED_NEIGHBOR_NUMBER 
            GetNeighborsLimited(Pos,2.0*Radius,Neighbors);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}

    for(int i=0;i<TempStellarWindLocalInfo.Nlist;i++){
        int leaf = Neighbors[i];
    
        TempStellarWindLocalInfo.k_hydro_localmin = fmin(TempStellarWindLocalInfo.k_hydro_localmin,Phydro[leaf]->k_hydro);

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
        double w = SPHKernel(r,InvRadiusi);

        TempStellarWindLocalInfo.WeightSum += w;
#ifdef MAXIMUM_ENERGY_INPUT
        if(i==0){
            TempStellarWindLocalInfo.DistanceMinGlobalID = PhydroBody(leaf)->GlobalID;
            TempStellarWindLocalInfo.DistanceMin = r;
        } else {
            if (TempStellarWindLocalInfo.DistanceMin > r){
                TempStellarWindLocalInfo.DistanceMinGlobalID = PhydroBody(leaf)->GlobalID;
                TempStellarWindLocalInfo.DistanceMin = r;
            }
        }
#endif
#ifdef USE_MOMENTUM_FEEDBACK //{
        TempStellarWindLocalInfo.nH_ave += Pall.ConvertNumberDensityToCGS*Phydro[leaf]->RhoPred;
        TempStellarWindLocalInfo.Z_ave += Phydro[leaf]->Z;
        TempStellarWindLocalInfo.NumberDensity += w;
#endif // USE_MOMENTUM_FEEDBACK //}
    }

    return TempStellarWindLocalInfo;
}
#endif // USE_STELLAR_WIND //}

#ifdef USE_NEIGHBOR_LIST_AND_SORT //{

struct StructSortNeighborList{
    int Index;
    double R2;
};

static int SortNeighrboListCmpR2(const void *x, const void *y){
    const struct StructSortNeighborList *pointer1 = x;
    const struct StructSortNeighborList *pointer2 = y;
    if( pointer1->R2 > pointer2->R2)
        return 1;
    else if( pointer1->R2 < pointer2->R2)
        return -1;
    else
        return 0;
}


static inline void SortNeighborListAndSetKernelSize(struct StructActiveParticle *AP_i, int Neighbors[]){
 
    if((NnbSort > AP_i->Nlist)&&(AP_i->Nlist>LocalNeighborListStride)){
        //set
        struct StructSortNeighborList SortNeighborList[AP_i->Nlist];
        for(int i=0;i<AP_i->Nlist;i++){
            SortNeighborList[i].Index = Neighbors[i];
            SortNeighborList[i].R2 = DISTANCE2(AP_i->Pos,Phydro[Neighbors[i]]->PosP);
        }

        // sort
        qsort(SortNeighborList,AP_i->Nlist,sizeof(struct StructSortNeighborList),
                 (int(*)(const void*, const void*))SortNeighrboListCmpR2);

        // fprintf(stderr,"Sort : %d %d %g %g\n",AP_i->Nlist,Pall.Ns,
                // AP_i->Radius,0.5*sqrt(SortNeighborList[AP_i->Nlist-1].R2));
        // pickup
        AP_i->Nlist = Pall.Ns;
        AP_i->Radius = 0.5*sqrt(SortNeighborList[AP_i->Nlist-1].R2);

        for(int i=0;i<AP_i->Nlist;i++){
            Neighbors[i] = SortNeighborList[i].Index;
        }
    }
 
    return ;
}
#endif // USE_NEIGHBOR_LIST_AND_SORT //}

#ifdef USE_NEIGHBOR_LIST //{
static inline void CopyNeighborList(struct StructActiveParticle *AP_i, const int Neighbors[], const int Offset){

    if(AP_i->Nlist < LocalNeighborListStride-1){
        LocalNeighborList[Offset+LocalNeighborListStride-1] = AP_i->Nlist;
        LNLK[Offset/LocalNeighborListStride] = AP_i->Radius;
        for(int k=0;k<AP_i->Nlist;k++){
            LocalNeighborList[Offset+k] = Neighbors[k];
        }

    }
    return ;
}

struct StructGetLocalNeighborList GetLocalNeighrborList(const int Index){

    struct StructGetLocalNeighborList NB;
    NB.Nlist = LocalNeighborList[LocalNeighborListStride*(Index+1)-1];
    NB.Neighbors = LocalNeighborList+LocalNeighborListStride*Index;
    NB.Kernel = LNLK[Index];
    
    return NB;
}

#endif // USE_NEIGHBOR_LIST //}



static inline void __attribute__((always_inline)) CS_OverwriteNeighborInfo(
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER  //{
        const double Nlist, 
#else 
        const int Nlist, 
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER  //}
        const int leaf){

    int Index = NBCache[leaf].Leaf;
    Phydro[Index]->Nlist = Nlist;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
    Phydro[Index]->SmoothedNumber = Nlist;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER
    Phydro[Index]->Kernel = 
    Phydro[Index]->KernelPred = NBCache[leaf].Kernel;

    return ;
}

static inline bool __attribute__((always_inline)) CS_CheckNeighborNumberAndUpdateHydroRadius_i(const int NProcs, bool CSExportFlags_i[restrict], struct StructActiveParticle *AP_i){ 
    // Index can be removed.

    const int NBmin = Pall.Ns-Pall.Npm;
    const int NBmax = Pall.Ns+Pall.Npm;
    const int NBminLimit = 0.8*NBmin;
    const int NBmaxLimit = 2.0*NBmax;
#ifdef USE_MAXIMUM_KERNEL_SIZE
#ifdef MAXIMUM_KERNEL_SIZE
    const double MaximumKernelSize = Pall.AdaptiveSofteningFactor*MAXIMUM_KERNEL_SIZE*KPC_CGS/Pall.UnitLength;
#else
#error Set MAXIMUM_KERNEL_SIZE
#endif
#endif

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double Nlist = AP_i->SmoothedNumber*
        SmoothedMassConversionFactor*CUBE(AP_i->Radius);
#else // USE_SMOOTHED_NEIGHBOR_NUMBER
    int Nlist = AP_i->Nlist;
#endif  // USE_SMOOTHED_NEIGHBOR_NUMBER //}

    if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
        CSExportFlags_i[NProcs-1] = false;
        // OverwriteNeighborInfo(Nlist,leaf);
        return true;
    }else if((Nlist>=(NBminLimit))&&(Nlist<=(NBmaxLimit))
            &&((AP_i->Rvalue*AP_i->Lvalue)>0.e0)){
            //&&(AP_i->Rvalue>0.e0)&&(AP_i->Lvalue>0.e0)){
        if(AP_i->Rvalue-AP_i->Lvalue < 1.e-3*AP_i->Lvalue){
            CSExportFlags_i[NProcs-1] = false;
            // OverwriteNeighborInfo(Nlist,leaf);
            return true;
        }
    }
#ifdef USE_MINIMUM_KERNEL_SIZE
    else if (((NBmax)<Nlist)&&(AP_i->Radius<=0.5*PhydroBody(AP_i->Index)->Eps*Pall.AdaptiveSofteningFactor)){
        AP_i->Radius = 0.5*PhydroBody(AP_i->Index)->Eps*Pall.AdaptiveSofteningFactor;
        CSExportFlags_i[NProcs-1] = false;
        // OverwriteNeighborInfo(Nlist,leaf);
        return true;
    }
#endif
#ifdef USE_MAXIMUM_KERNEL_SIZE
    else if (((NBmin)>Nlist)&&(AP_i->Radius>MaximumKernelSize)){
        AP_i->Radius = MaximumKernelSize;
        CSExportFlags_i[NProcs-1] = false;
        //OverwriteNeighborInfo(Nlist,leaf);
        return true;
    }
#endif
    if(CSExportFlags_i[NProcs-1]){
        if(Nlist<NBmin){
            AP_i->Lvalue = fmax(AP_i->Lvalue,AP_i->Radius);
        } else if(Nlist>NBmax){
            if(AP_i->Rvalue > 0.e0){
                AP_i->Rvalue = fmin(AP_i->Rvalue,AP_i->Radius);
            }else{
                AP_i->Rvalue = AP_i->Radius;
            }
        }

        if((AP_i->Lvalue*AP_i->Rvalue)>0.e0){
            AP_i->Radius = cbrt(0.5*(CUBE(AP_i->Lvalue)+CUBE(AP_i->Rvalue)));
        }else{
            if((AP_i->Rvalue == 0.e0)&&(AP_i->Lvalue > 0.e0)){
                AP_i->Radius *= RadiusFactInc;
            }else if((AP_i->Rvalue > 0.e0)&&(AP_i->Lvalue == 0.e0)){
                AP_i->Radius *= RadiusFactDec;
            }
        }
    }
    
    return false;
}

// anker
static inline bool __attribute__((always_inline)) CS_CheckLocalMassAndUpdateHIIRadius_i(const int NProcs, bool CSExportFlags_i[NProcs], struct StructActiveParticle *AP_i){ 

#define ConvergenceFactor  (1.e-2)

    const double Rmax = 50*PC_CGS/Pall.UnitLength;

#ifdef __PHOTON_COUNT_BASE__
    double Qest = AP_i->Body.HII.PhotonCount;
#else // __PHOTON_COUNT_BASE__
    double fact = (3.0/(4.0*M_PI))*(ReturnAlpha()/(SQ(PROTON_MASS_CGS)));
    double Qest = fact*(SQ((AP_i->Body.HII.Mass)*Pall.UnitMass)/(CUBE((2.0*(AP_i->Radius))*Pall.UnitLength)));
#endif // __PHOTON_COUNT_BASE__
    double Qmin = (1.0-ConvergenceFactor)*AP_i->Body.HII.LyAphoton;
    double Qmax = (1.0+ConvergenceFactor)*AP_i->Body.HII.LyAphoton;

    if(Niteration_for_debug > 100)
        fprintf(stderr,"Qest/// %d : %g %g %g // %g %g %g\n",AP_i->Nlist,Qmin,Qest,Qmax,
                AP_i->Lvalue,AP_i->Radius,AP_i->Rvalue);

#ifdef USE_HIIREGION_ITERATION_LIMIT
    if(HIIREGION_ITERATION_LIMIT_COUNT < AP_i->Counter){
        CSExportFlags_i[NProcs-1] = false;
        if(Qest > Qmax){
            AP_i->Radius = 0.e0;
        } 
        return true;
    } else 
#endif
    if(((Qmin)<=Qest)&&(Qest<=(Qmax))){
        CSExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__
        AP_i->Body.HII.PhotonCount = 1.0;
#else // __PHOTON_COUNT_BASE__
        AP_i->Body.HII.Mass = 1.0;
#endif // __PHOTON_COUNT_BASE__
        return true;
    } else if((AP_i->Nlist < 5)&&(Qest < Qmax)){
        CSExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__
        AP_i->Body.HII.PhotonCount = 6.0;
#else // __PHOTON_COUNT_BASE__
        AP_i->Body.HII.Mass = 6.0;
#endif // __PHOTON_COUNT_BASE__
        AP_i->Radius *= cbrt(Qest/Qmax);
        return true;
    } else if((AP_i->Radius > Rmax)&&(Qest < Qmax)){
        CSExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__
        AP_i->Body.HII.PhotonCount = 2.0;
#else // __PHOTON_COUNT_BASE__
        AP_i->Body.HII.Mass = 2.0;
#endif // __PHOTON_COUNT_BASE__
        return true;
    } else if(AP_i->Lvalue > Rmax){
        CSExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__ //{
        AP_i->Body.HII.PhotonCount = 3.0;
#else // __PHOTON_COUNT_BASE__ // { //}
        AP_i->Body.HII.Mass = 3.0;
#endif // __PHOTON_COUNT_BASE__ //}
        AP_i->Radius = AP_i->Lvalue;
        return true;
    } else if((AP_i->Nlist == 1)&&(Qest < Qmax)){
        CSExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__
        AP_i->Body.HII.PhotonCount = 7.0;
#else // __PHOTON_COUNT_BASE__
        AP_i->Body.HII.Mass = 7.0;
#endif // __PHOTON_COUNT_BASE__
        AP_i->Body.HII.HIIRegion = false;
        return true;
    } else if((AP_i->Nlist == 0)&&(AP_i->Radius > Rmax)){
        CSExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__ //{
        AP_i->Body.HII.PhotonCount = 4.0;
#else // __PHOTON_COUNT_BASE__ // { //}
        AP_i->Body.HII.Mass = 4.0;
#endif // __PHOTON_COUNT_BASE__ //}
        AP_i->Body.HII.HIIRegion = false;
        return true;
    } else if((Qest<=Qmax)&&((AP_i->Rvalue*AP_i->Lvalue)>0.e0)){
        if(AP_i->Rvalue-AP_i->Lvalue < 1.e-3*AP_i->Rvalue){
            CSExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__ //{
            AP_i->Body.HII.PhotonCount = 5.0;
#else // __PHOTON_COUNT_BASE__ // { //}
            AP_i->Body.HII.Mass = 5.0;
#endif // __PHOTON_COUNT_BASE__ //}
            if(Qest == 0.0){
                AP_i->Body.HII.HIIRegion = false;
            }
            return true;
        }
    } 
    else if((Qest>Qmax)&&(AP_i->Nlist == 1)){
        CSExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__
        AP_i->Body.HII.PhotonCount = 6.0;
#else // __PHOTON_COUNT_BASE__
        AP_i->Body.HII.Mass = 6.0;
#endif // __PHOTON_COUNT_BASE__
        AP_i->Radius *= cbrt(Qmax/Qest);
        return true;
    }
#ifdef __PHOTON_COUNT_BASE__ //{
        else if((AP_i->Body.HII.PhotonCountDistanceMin > AP_i->Body.HII.LyAphoton))
#else // __PHOTON_COUNT_BASE__ // { //}
#error USE PHOTON COUNT BASE MODE
#endif // __PHOTON_COUNT_BASE__ //}
        {
            CSExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__ //{
            AP_i->Body.HII.PhotonCount = 4.0;
#else // __PHOTON_COUNT_BASE__ // { //}
            AP_i->Body.HII.Mass = 4.0;
#endif // __PHOTON_COUNT_BASE__ //}
            AP_i->Body.HII.HIIRegion = false;
            return true;
    }


    if(CSExportFlags_i[NProcs-1]){
        if(Qest<Qmin){
            AP_i->Lvalue = fmax(AP_i->Lvalue,AP_i->Radius);
        } else if(Qest>Qmax){
            if((AP_i->Rvalue*AP_i->Lvalue > 0.e0)
                    &&(AP_i->Rvalue-AP_i->Lvalue < 1.e-3*AP_i->Rvalue)){
                AP_i->Lvalue *= 0.5;
            }else{
                if(AP_i->Rvalue > 0.e0){
                    AP_i->Rvalue = fmin(AP_i->Rvalue,AP_i->Radius);
                }else{
                    AP_i->Rvalue = AP_i->Radius;
                }
            }
        }

        if((AP_i->Lvalue*AP_i->Rvalue)>0.e0){
            AP_i->Radius = cbrt(0.5*(CUBE(AP_i->Lvalue)+CUBE(AP_i->Rvalue)));
        }else{
            if((AP_i->Rvalue == 0.e0)&&(AP_i->Lvalue > 0.e0)){
                //AP_i->Radius *= RadiusFactInc_First;
                AP_i->Radius *= RadiusFactInc;
            }else if((AP_i->Rvalue > 0.e0)&&(AP_i->Lvalue == 0.e0)){
                AP_i->Radius *= RadiusFactDec;
            }
        }
    }
    if(Niteration_for_debug > 100){
        fprintf(stderr,"NewValue /// %d %g %g %g\n",AP_i->Nlist,
                AP_i->Lvalue,AP_i->Radius,AP_i->Rvalue);
        fflush(NULL);
    }

    AP_i->Counter ++;

    return false;
}

//#define StellarFeedbackRadiusFactInc   (1.14) // 1.5 ^ (1.3)
#define StellarFeedbackRadiusFactInc   (1.2596) //  2^(0.333)
#define StellarFeedbackRadiusFactDec   (0.79) // 0.75 ^ (1.3)
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
static double SmoothedMassConversionFactor;
#endif //USE_SMOOTHED_NEIGHBOR_NUMBER
/*
 * This function checks the convergence of the feedback radius, using bisetion
 * method. If the radius is convergered, this function returns a "true" flag.
 * If not, it returns a "false" flag.
 */
static inline bool __attribute__((always_inline)) CS_CheckNeighborNumberAndUpdateStellarFeedbackRadius_i(const int NProcs, bool CSExportFlags_i[restrict], struct StructActiveParticle *AP_i){ 

#ifdef USE_SN_INPUT_PARTICLE_NUMBER //{
    const int NBmin = SN_INPUT_PARTICLE_NUMBER-SN_INPUT_PARTICLE_NUMBER_MARGIN;
    const int NBmax = SN_INPUT_PARTICLE_NUMBER+SN_INPUT_PARTICLE_NUMBER_MARGIN;
#else // USE_SN_INPUT_PARTICLE_NUMBER
    const int NBmin = Pall.Ns-Pall.Npm;
    const int NBmax = Pall.Ns+Pall.Npm;
#endif // USE_SN_INPUT_PARTICLE_NUMBER //}

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double Nlist = AP_i->SmoothedNumber*
        SmoothedMassConversionFactor*CUBE(AP_i->Radius);
#else // USE_SMOOTHED_NEIGHBOR_NUMBER
    int Nlist = AP_i->Nlist;
#endif  // USE_SMOOTHED_NEIGHBOR_NUMBER //}

#ifdef USE_SNE_ITERATION_LIMIT
    if(SNE_ITERATION_LIMIT_COUNT < AP_i->Counter){
        CSExportFlags_i[NProcs-1] = false;
        return true;
    } else 
#endif
    // Here, the convergence condition is satisfied. 
    if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
        CSExportFlags_i[NProcs-1] = false;
        return true;
    }else if(((0.8*NBmin)<=Nlist)&&(Nlist<=(2.0*NBmax))&&(AP_i->Rvalue>0.e0)&&(AP_i->Lvalue>0.e0)){
        if(AP_i->Rvalue-AP_i->Lvalue < 1.e-6*AP_i->Lvalue){
            CSExportFlags_i[NProcs-1] = false;
            return true;
        }
    }
    if(CSExportFlags_i[NProcs-1]){
        if(Nlist<NBmin){
            AP_i->Lvalue = fmax(AP_i->Lvalue,AP_i->Radius);
        } else if(Nlist>NBmax){
            if(AP_i->Rvalue > 0.e0){
                AP_i->Rvalue = fmin(AP_i->Rvalue,AP_i->Radius);
            }else{
                AP_i->Rvalue = AP_i->Radius;
            }
        }

        if((AP_i->Lvalue>0.e0)&&(AP_i->Rvalue>0.e0)){
            AP_i->Radius = cbrt(0.5*(CUBE(AP_i->Lvalue)+CUBE(AP_i->Rvalue)));
        }else{
            if((AP_i->Rvalue == 0.e0)&&(AP_i->Lvalue > 0.e0)){
                AP_i->Radius *= StellarFeedbackRadiusFactInc;
            }else if((AP_i->Rvalue > 0.e0)&&(AP_i->Lvalue == 0.e0)){
                AP_i->Radius *= StellarFeedbackRadiusFactDec;
            }
        }
    }

    AP_i->Counter ++;

    return false;
}

#ifdef USE_RADIATION_PRESSURE //{
static inline bool __attribute__((always_inline)) CS_CheckNeighborNumberAndUpdateRadiationPressureRadius_i(const int NProcs, bool CSExportFlags_i[restrict], struct StructActiveParticle *AP_i){ 

    const int NBmin = RADIATION_PRESSURE_INPUT_PARTICLE_NUMBER-RADIATION_PRESSURE_INPUT_PARTICLE_NUMBER_MERGIN;
    const int NBmax = RADIATION_PRESSURE_INPUT_PARTICLE_NUMBER+RADIATION_PRESSURE_INPUT_PARTICLE_NUMBER_MERGIN;

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double Nlist = AP_i->SmoothedNumber*
        SmoothedMassConversionFactor*CUBE(AP_i->Radius);
#else // USE_SMOOTHED_NEIGHBOR_NUMBER
    int Nlist = AP_i->Nlist;
#endif  // USE_SMOOTHED_NEIGHBOR_NUMBER //}

    if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
        CSExportFlags_i[NProcs-1] = false;
        // OverwriteNeighborInfo(Nlist,leaf);
        return true;
    }else if((Nlist>=(0.8*NBmin))&&(Nlist<=(2*NBmax))
            &&(AP_i->Rvalue>0.e0)&&(AP_i->Lvalue>0.e0)){
        if(AP_i->Rvalue-AP_i->Lvalue < 1.e-6*AP_i->Lvalue){
            CSExportFlags_i[NProcs-1] = false;
            return true;
        }
    }
    if(CSExportFlags_i[NProcs-1]){
        if(Nlist<NBmin){
            AP_i->Lvalue = fmax(AP_i->Lvalue,AP_i->Radius);
        } else if(Nlist>NBmax){
            if(AP_i->Rvalue > 0.e0){
                AP_i->Rvalue = fmin(AP_i->Rvalue,AP_i->Radius);
            }else{
                AP_i->Rvalue = AP_i->Radius;
            }
        }

        if((AP_i->Lvalue>0.e0)&&(AP_i->Rvalue>0.e0)){
            AP_i->Radius = cbrt(0.5*(CUBE(AP_i->Lvalue)+CUBE(AP_i->Rvalue)));
        }else{
            if((AP_i->Rvalue == 0.e0)&&(AP_i->Lvalue > 0.e0)){
                AP_i->Radius *= RadiusFactInc;
            }else if((AP_i->Rvalue > 0.e0)&&(AP_i->Lvalue == 0.e0)){
                AP_i->Radius *= RadiusFactDec;
            }
        }
    }
    
    return false;
}
#endif // USE_RADIATION_PRESSURE //}

#ifdef USE_STELLAR_WIND //{
#define StellarWindRadiusFactInc (StellarFeedbackRadiusFactInc) //  2^(0.333)
#define StellarWindRadiusFactDec (StellarFeedbackRadiusFactDec) // 0.75 ^ (1.3)
static inline bool __attribute__((always_inline)) CS_CheckNeighborNumberAndUpdateStellarWindRadius_i(const int NProcs, bool CSExportFlags_i[restrict], struct StructActiveParticle *AP_i){ 

#ifdef USE_SN_INPUT_PARTICLE_NUMBER //{
    const int NBmin = SN_INPUT_PARTICLE_NUMBER-SN_INPUT_PARTICLE_NUMBER_MARGIN;
    const int NBmax = SN_INPUT_PARTICLE_NUMBER+SN_INPUT_PARTICLE_NUMBER_MARGIN;
#else // USE_SN_INPUT_PARTICLE_NUMBER
    const int NBmin = Pall.Ns-Pall.Npm;
    const int NBmax = Pall.Ns+Pall.Npm;
#endif // USE_SN_INPUT_PARTICLE_NUMBER //}

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double Nlist = AP_i->SmoothedNumber*
        SmoothedMassConversionFactor*CUBE(AP_i->Radius);
#else // USE_SMOOTHED_NEIGHBOR_NUMBER
    int Nlist = AP_i->Nlist;
#endif  // USE_SMOOTHED_NEIGHBOR_NUMBER //}

    // Here, the convergence condition is satisfied. 
    if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
        CSExportFlags_i[NProcs-1] = false;
        return true;
    }else if(((0.8*NBmin)<=Nlist)&&(Nlist<=(2.0*NBmax))&&(AP_i->Rvalue>0.e0)&&(AP_i->Lvalue>0.e0)){
        if(AP_i->Rvalue-AP_i->Lvalue < 1.e-6*AP_i->Lvalue){
            CSExportFlags_i[NProcs-1] = false;
            return true;
        }
    }
    if(CSExportFlags_i[NProcs-1]){
        if(Nlist<NBmin){
            AP_i->Lvalue = fmax(AP_i->Lvalue,AP_i->Radius);
        } else if(Nlist>NBmax){
            if(AP_i->Rvalue > 0.e0){
                AP_i->Rvalue = fmin(AP_i->Rvalue,AP_i->Radius);
            }else{
                AP_i->Rvalue = AP_i->Radius;
            }
        }

        if((AP_i->Lvalue>0.e0)&&(AP_i->Rvalue>0.e0)){
            AP_i->Radius = cbrt(0.5*(CUBE(AP_i->Lvalue)+CUBE(AP_i->Rvalue)));
        }else{
            if((AP_i->Rvalue == 0.e0)&&(AP_i->Lvalue > 0.e0)){
                AP_i->Radius *= StellarWindRadiusFactInc;
            }else if((AP_i->Rvalue > 0.e0)&&(AP_i->Lvalue == 0.e0)){
                AP_i->Radius *= StellarWindRadiusFactDec;
            }
        }
    }
    return false;
}
#endif // USE_STELLAR_WIND //}

static inline bool __attribute__((always_inline)) CS_CheckConvergence_i(const int NProcs, bool CSExportFlags_i[restrict], struct StructActiveParticle *AP_i, const int Type){ 
    if(AP_i->Type == CS_TypeHydro){
        return CS_CheckNeighborNumberAndUpdateHydroRadius_i(NProcs,CSExportFlags_i,AP_i); 
    } else if(AP_i->Type == CS_TypeHII){
        return CS_CheckLocalMassAndUpdateHIIRadius_i(NProcs,CSExportFlags_i,AP_i); 
    } else if(AP_i->Type == CS_TypeSN){
        return CS_CheckNeighborNumberAndUpdateStellarFeedbackRadius_i(NProcs,CSExportFlags_i,AP_i); 
    }
#ifdef USE_RADIATION_PRESSURE //{
    else if(AP_i->Type == CS_TypeRP){
        return CS_CheckNeighborNumberAndUpdateRadiationPressureRadius_i(NProcs,CSExportFlags_i,AP_i); 
    }
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
    else if(AP_i->Type == CS_TypeSW){
        return CS_CheckNeighborNumberAndUpdateStellarWindRadius_i(NProcs,CSExportFlags_i,AP_i); 
    }
#endif // USE_STELLAR_WIND //}
    return false;
}


static void CS_UpdateRadiusAndOthers_i(struct StructActiveParticle *AP_i, int Neighbors[restrict]){

    if(AP_i->Type == CS_TypeHydro){
#if 0
        AP_i->Nlist = 
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            CS_GetSmoothedNumberofNeighbors(AP_i->Pos,AP_i->Radius,
                    Neighbors,&(AP_i->SmoothedNumber));
#else // USE_SMOOTHED_NEIGHBOR_NUMBER 
            CS_GetNumberofNeighbors(AP_i->Pos,AP_i->Radius,Neighbors);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
        // if(AP_i->GlobalID == 6){
            // fprintf(stderr," -- %d %g -- \n",AP_i->Nlist,AP_i->Radius);
        // }
#endif 
        struct StructCSHydroLocalInfo TemporalData = 
            CS_GetNumberofNeighbors(AP_i->Pos,AP_i->Radius,Neighbors);
        AP_i->Nlist = TemporalData.Nlist;
        AP_i->k_hydro_localmin = TemporalData.k_hydro_localmin;

#ifdef USE_SINK_PARTICLE //{
        AP_i->Body.Hydro.PotentialMin = TemporalData.PotentialMin;
        AP_i->Body.Hydro.DensityMax = TemporalData.DensityMax;
        AP_i->Body.Hydro.MassTotal = TemporalData.MassTotal;
        AP_i->Body.Hydro.VCOM[0] = TemporalData.VCOM[0];
        AP_i->Body.Hydro.VCOM[1] = TemporalData.VCOM[1];
        AP_i->Body.Hydro.VCOM[2] = TemporalData.VCOM[2];
#ifdef USE_SINK_NOLOCALSINK_CONDITION //{
        AP_i->Body.Hydro.NoLocalSink = TemporalData.NoLocalSink;
#endif // USE_SINK_NOLOCALSINK_CONDITION //}
#endif // USE_SINK_PARTICLE //}

    } else if(AP_i->Type == CS_TypeHII){
        struct StructCSHIILocalInfo TemporalData =
            CS_ReturnHIILocalInfo(AP_i->Pos,AP_i->Radius,Neighbors);

        AP_i->Nlist = TemporalData.Nlist;
        AP_i->k_hydro_localmin = TemporalData.k_hydro_localmin;
#ifdef __PHOTON_COUNT_BASE__
        AP_i->Body.HII.PhotonCount = TemporalData.PhotonCount;
        AP_i->Body.HII.PhotonCountDistanceMin = TemporalData.PhotonCountDistanceMin;
#else //__PHOTON_COUNT_BASE__
        AP_i->Body.HII.Mass = TemporalData.Mass;
        AP_i->Body.HII.MassDistanceMin = TemporalData.MassDistanceMin;
#endif // __PHOTON_COUNT_BASE__
        AP_i->Body.HII.DistanceMin = TemporalData.DistanceMin;
    } else if(AP_i->Type == CS_TypeSN){
        struct StructCSStellarFeedbackLocalInfo TemporalData = 
            CS_RetrunStellarFeedbackLocalInfo(AP_i->Pos,AP_i->Radius,Neighbors);

        AP_i->Nlist = TemporalData.Nlist;
        AP_i->k_hydro_localmin = TemporalData.k_hydro_localmin;
        AP_i->Body.SN.Density = TemporalData.Density;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        AP_i->Body.SN.SmoothedNumber = TemporalData.SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
        AP_i->Body.SN.GasMass = TemporalData.GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
        if(TemporalData.Nlist > 0){
            AP_i->Body.SN.DistanceMin = TemporalData.DistanceMin;
            AP_i->Body.SN.DistanceMinGlobalID = TemporalData.DistanceMinGlobalID;
        }
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef USE_MOMENTUM_FEEDBACK //{
        AP_i->Body.SN.nH_ave = TemporalData.nH_ave;
        AP_i->Body.SN.Z_ave = TemporalData.Z_ave;
        AP_i->Body.SN.NumberDensity = TemporalData.NumberDensity;
#endif // USE_MOMENTUM_FEEDBACK //}
        AP_i->Body.SN.WeightSum = TemporalData.WeightSum;
#ifdef __CHECK_SUM__ //{
        AP_i->Body.SN.CheckSum = TemporalData.CheckSum;
#endif //__CHECK_SUM__ //}
    }
#ifdef USE_RADIATION_PRESSURE //{
    else if(AP_i->Type == CS_TypeRP){
        struct StructCSRadiationPressureLocalInfo TemporalData = 
            CS_RetrunRadiationPressureLocalInfo(AP_i->Pos,AP_i->Radius,Neighbors);

        AP_i->Nlist = TemporalData.Nlist;
        AP_i->k_hydro_localmin = TemporalData.k_hydro_localmin;
        AP_i->Body.RP.GasMass = TemporalData.GasMass;
        AP_i->Body.RP.MetalMass = TemporalData.MetalMass;
        AP_i->Body.RP.WeightSum = TemporalData.WeightSum;
    }
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
    else if(AP_i->Type == CS_TypeSW){
        struct StructCSStellarWindLocalInfo TemporalData = 
            CS_RetrunStellarWindLocalInfo(AP_i->Pos,AP_i->Radius,Neighbors);

        AP_i->Nlist = TemporalData.Nlist;
        AP_i->k_hydro_localmin = TemporalData.k_hydro_localmin;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        AP_i->Body.SW.SmoothedNumber = TemporalData.SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef MAXIMUM_ENERGY_INPUT //{
        if(TemporalData.Nlist > 0){
            AP_i->Body.SW.DistanceMin = TemporalData.DistanceMin;
            AP_i->Body.SW.DistanceMinGlobalID = TemporalData.DistanceMinGlobalID;
        }
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef USE_MOMENTUM_FEEDBACK //{
        AP_i->Body.SW.nH_ave = TemporalData.nH_ave;
        AP_i->Body.SW.Z_ave = TemporalData.Z_ave;
        AP_i->Body.SW.NumberDensity = TemporalData.NumberDensity;
#endif // USE_MOMENTUM_FEEDBACK //}
        AP_i->Body.SW.WeightSum = TemporalData.WeightSum;
    }
#endif // USE_STELLAR_WIND //}
    return ;
}


static void CS_UpdateRadiusAndOthersImported_i(struct StructCSImport *CSImportSend_i, struct StructCSExport *CSExportRecv_i, int Neighbors[restrict]){

    if(CSExportRecv_i->Type == CS_TypeHydro){
        struct StructCSHydroLocalInfo TemporalData = 
            CS_GetNumberofNeighbors(CSExportRecv_i->Pos,CSExportRecv_i->Radius,Neighbors);
        CSImportSend_i->Nlist = TemporalData.Nlist;
        CSImportSend_i->k_hydro_localmin = TemporalData.k_hydro_localmin;

#ifdef USE_SINK_PARTICLE //{
        CSImportSend_i->Body.Hydro.PotentialMin = TemporalData.PotentialMin;
        CSImportSend_i->Body.Hydro.DensityMax = TemporalData.DensityMax;
        CSImportSend_i->Body.Hydro.MassTotal = TemporalData.MassTotal;
        CSImportSend_i->Body.Hydro.VCOM[0] = TemporalData.VCOM[0];
        CSImportSend_i->Body.Hydro.VCOM[1] = TemporalData.VCOM[1];
        CSImportSend_i->Body.Hydro.VCOM[2] = TemporalData.VCOM[2];
#ifdef USE_SINK_NOLOCALSINK_CONDITION //{
        CSImportSend_i->Body.Hydro.NoLocalSink = TemporalData.NoLocalSink;
#endif // USE_SINK_NOLOCALSINK_CONDITION //}
#endif // USE_SINK_PARTICLE //}

    } else if(CSExportRecv_i->Type == CS_TypeHII){
        struct StructCSHIILocalInfo TemporalData =
            CS_ReturnHIILocalInfo(CSExportRecv_i->Pos,CSExportRecv_i->Radius,Neighbors);

        CSImportSend_i->Nlist = TemporalData.Nlist;
        CSImportSend_i->k_hydro_localmin = TemporalData.k_hydro_localmin;
#ifdef __PHOTON_COUNT_BASE__
        CSImportSend_i->Body.HII.PhotonCount = TemporalData.PhotonCount;
        CSImportSend_i->Body.HII.PhotonCountDistanceMin = TemporalData.PhotonCountDistanceMin;
#else //__PHOTON_COUNT_BASE__
        CSImportSend_i->Body.HII.Mass = TemporalData.Mass;
        CSImportSend_i->Body.HII.MassDistMin = TemporalData.MassDistanceMin;
#endif // __PHOTON_COUNT_BASE__
        CSImportSend_i->Body.HII.DistanceMin = TemporalData.DistanceMin;

    } else if(CSExportRecv_i->Type == CS_TypeSN){
        struct StructCSStellarFeedbackLocalInfo TemporalData = 
            CS_RetrunStellarFeedbackLocalInfo(CSExportRecv_i->Pos,CSExportRecv_i->Radius,Neighbors);

        CSImportSend_i->Nlist = TemporalData.Nlist;
        CSImportSend_i->k_hydro_localmin = TemporalData.k_hydro_localmin;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        CSImportSend_i->SmoothedNumber = TemporalData.SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
        CSImportSend_i->Body.SN.Density = TemporalData.Density;
#ifdef SET_SNII_TEMPERATURE //{
        CSImportSend_i->Body.SN.GasMass = TemporalData.GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
        if(TemporalData.Nlist > 0){
            CSImportSend_i->Body.SN.DistanceMin = TemporalData.DistanceMin;
            CSImportSend_i->Body.SN.DistanceMinGlobalID = TemporalData.DistanceMinGlobalID;
        }
#endif // MAXIMUM_ENERGY_INPUT //}

#ifdef USE_MOMENTUM_FEEDBACK //{  
        CSImportSend_i->Body.SN.nH_ave = TemporalData.nH_ave;
        CSImportSend_i->Body.SN.Z_ave = TemporalData.Z_ave;
        CSImportSend_i->Body.SN.NumberDensity = TemporalData.NumberDensity;
#endif // USE_MOMENTUM_FEEDBACK //}

        CSImportSend_i->Body.SN.WeightSum = TemporalData.WeightSum;
#ifdef __CHECK_SUM__ //{
        CSImportSend_i->Body.SN.CheckSum = TemporalData.CheckSum;
#endif //__CHECK_SUM__ //}
    }
#ifdef USE_RADIATION_PRESSURE //{
    else if(CSExportRecv_i->Type == CS_TypeRP){
        struct StructCSRadiationPressureLocalInfo TemporalData = 
            CS_RetrunRadiationPressureLocalInfo(CSExportRecv_i->Pos,CSExportRecv_i->Radius,Neighbors);

        CSImportSend_i->Nlist = TemporalData.Nlist;
        CSImportSend_i->k_hydro_localmin = TemporalData.k_hydro_localmin;
        CSImportSend_i->Body.RP.GasMass = TemporalData.GasMass;
        CSImportSend_i->Body.RP.MetalMass = TemporalData.MetalMass;
        CSImportSend_i->Body.RP.WeightSum = TemporalData.WeightSum;
    }
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
    else if(CSExportRecv_i->Type == CS_TypeSW){
        struct StructCSStellarWindLocalInfo TemporalData = 
            CS_RetrunStellarWindLocalInfo(CSExportRecv_i->Pos,CSExportRecv_i->Radius,Neighbors);

        CSImportSend_i->Nlist = TemporalData.Nlist;
        CSImportSend_i->k_hydro_localmin = TemporalData.k_hydro_localmin;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        CSImportSend_i->SmoothedNumber = TemporalData.SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef MAXIMUM_ENERGY_INPUT //{
        if(TemporalData.Nlist > 0){
            CSImportSend_i->Body.SW.DistanceMin = TemporalData.DistanceMin;
            CSImportSend_i->Body.SW.DistanceMinGlobalID = TemporalData.DistanceMinGlobalID;
        }
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef USE_MOMENTUM_FEEDBACK //{  
        CSImportSend_i->Body.SW.nH_ave = TemporalData.nH_ave;
        CSImportSend_i->Body.SW.Z_ave = TemporalData.Z_ave;
        CSImportSend_i->Body.SW.NumberDensity = TemporalData.NumberDensity;
#endif // USE_MOMENTUM_FEEDBACK //}
        CSImportSend_i->Body.SW.WeightSum = TemporalData.WeightSum;
    }
#endif // USE_STELLAR_WIND //}

    CSImportSend_i->Leaf = CSExportRecv_i->Leaf;
    CSImportSend_i->Type = CSExportRecv_i->Type;
#ifdef USE_DEBUG_MODE //{
    CSImportSend_i->GlobalID = CSExportRecv_i->GlobalID;
#endif // USE_DEBUG_MODE //}
    return ;
}


static inline bool  __attribute__((always_inline)) CS_CheckInLocalDomain(double Pos[], double Extent /* Radiusx2 */, const int MyID){
    for(int k=0;k<3;k++){
        if(Pos[k]+Extent > EdgesForHydro[MyID].PosMax[k]) return false;
        if(Pos[k]-Extent < EdgesForHydro[MyID].PosMin[k]) return false;
    }
    return true;
}

static inline bool __attribute__((always_inline)) CS_OverlapDomainKernel(double Pos[restrict], const double h /* 2*Radius */, const int NodeID){ 

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForHydro[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

static inline void __attribute__((always_inline)) CS_UpdateRadiusLocal(struct StructActiveParticle *AP_i, const int MyID, const int NProcs, bool CSExportFlags_i[restrict], const int Offset){

    if(CS_CheckConvergence_i(NProcs,CSExportFlags_i,AP_i,AP_i->Type) == true)
        return;

    int Neighbors[MaxNeighborSize];
    int counter = 0;
    do{
        if(!CS_CheckInLocalDomain(AP_i->Pos,2.0*AP_i->Radius,MyID)) return;
        for(int i=0;i<CS_NContactedDomains;i++){
            int NodeID = CS_ContactedDomainID[i];
            if(CS_OverlapDomainKernel(AP_i->Pos,2.0*AP_i->Radius,NodeID)) return;
        }

        CS_UpdateRadiusAndOthers_i(AP_i,Neighbors);

#ifdef USE_NEIGHBOR_LIST_AND_SORT //{
        SortNeighborListAndSetKernelSize(AP_i,Neighbors);
#endif // USE_NEIGHBOR_LIST_AND_SORT //}
#ifdef USE_NEIGHBOR_LIST //{
        CopyNeighborList(AP_i,Neighbors,Offset);
#endif // USE_NEIGHBOR_LIST //}

        counter ++;
        if(counter > 10) return ;
    }while(CS_CheckConvergence_i(NProcs,CSExportFlags_i,AP_i,AP_i->Type) == false);

    return;
}


static inline int __attribute__((always_inline)) CS_CheckConvergence(const int NActives, const int NProcs, bool CSExportFlags[restrict], struct StructActiveParticle AP[restrict]){ 


    int NLocalActiveLeaves = 0;
    for(int i=0;i<NActives;i++){
        int Offset = i*NProcs;
        if(CSExportFlags[Offset+NProcs-1]){ 
#ifdef UPDATE_SIZE_LOCAL //{
        if(AP[i].LocalUpdateFlags == false){ 
#endif // UPDATE_SIZE_LOCAL //} 
            if(CS_CheckConvergence_i(NProcs,CSExportFlags+Offset,AP+i,AP[i].Type) == false){
                NLocalActiveLeaves ++;
            }
#ifdef UPDATE_SIZE_LOCAL //{
        } else {
            NLocalActiveLeaves ++;
        }
#endif // UPDATE_SIZE_LOCAL //}
        }
    }

    return NLocalActiveLeaves;
}

#ifdef ADD_PERTURBATION //{
#define _ResetTiming_ 1
static void ResetLRvalues(const int NActives, const int Niteration, const int NProcs, bool CSExportFlags[restrict], struct StructActiveParticle AP[restrict]){ 

    if((Niteration+1)%(_ResetTiming_*MaxIterationTimes) == 0){
#if 0
        for(int i=0;i<NActives;i++){
            int Offset = i*NProcs;
            if(CSExportFlags[Offset+NProcs-1]){ 
                AP[i].Rvalue = AP[i].Lvalue = 0.e0;
            }
        }
#else 
        for(int i=0;i<NActives;i++){
            int Offset = i*NProcs;
            if(CSExportFlags[Offset+NProcs-1]){ 
                AP[i].Lvalue *= 0.5;
                AP[i].Rvalue *= 2.0; 
            }
        }
#endif
    }
    return ;
}
#endif // ADD_PERTURBATION //}

static void CalcMomentInv(double B[3][3]){

#if (DIMENSION == 2)
    B[2][2] = 1.0;
#endif // DIMENSION == 2

#if (DIMENSION == 1)
    B[1][1] = B[2][2] = 1.0;
#endif // DIMENSION == 2

    double C[3][3];
    double det = B[0][0]*B[1][1]*B[2][2] + B[1][0]*B[2][1]*B[0][2] + B[2][0]*B[0][1]*B[1][2]
                -B[0][0]*B[2][1]*B[1][2] - B[2][0]*B[1][1]*B[0][2] - B[1][0]*B[0][1]*B[2][2];
    det += 1.e-16;

    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            C[i][j] = B[i][j];
        }
    }

    B[0][0] = (C[1][1]*C[2][2]-C[1][2]*C[2][1])/det;
    B[0][1] = (C[0][2]*C[2][1]-C[0][1]*C[2][2])/det;
    B[0][2] = (C[0][1]*C[1][2]-C[0][2]*C[1][1])/det;
    B[1][0] = (C[1][2]*C[2][0]-C[1][0]*C[2][2])/det;
    B[1][1] = (C[0][0]*C[2][2]-C[0][2]*C[2][0])/det;
    B[1][2] = (C[0][2]*C[1][0]-C[0][0]*C[1][2])/det;
    B[2][0] = (C[1][0]*C[2][1]-C[1][1]*C[2][0])/det;
    B[2][1] = (C[0][1]*C[2][0]-C[0][0]*C[2][1])/det;
    B[2][2] = (C[0][0]*C[1][1]-C[0][1]*C[1][0])/det;

    return ;
}

extern struct StructActiveStellarFeedbackParticle *ActiveStellarFeedbackParticle;

static void RestoreAllData(const int NActives, struct StructActiveParticle AP[restrict]){

    for(int i=0;i<NActives;i++){
        if(AP[i].Type == CS_TypeHydro){
            int Index = AP[i].Index;
            Phydro[Index]->Kernel = Phydro[Index]->KernelPred = AP[i].Radius;
            int CacheIndex = AP[i].Body.Hydro.CacheIndex;
            NBCache[CacheIndex].Kernel = AP[i].Radius;
        } else if(AP[i].Type == CS_TypeHII){
            int Index = AP[i].Index;
#ifdef __PHOTON_COUNT_BASE__
            Pstar[Index]->Density = AP[i].Body.HII.PhotonCount;
#else // __PHOTON_COUNT_BASE__
            Pstar[Index]->Density = AP[i].Body.HII.Mass;
#endif // __PHOTON_COUNT_BASE__
            if(AP[i].Body.HII.HIIRegion){
                Pstar[Index]->StromgrenRadius = 2.0*AP[i].Radius;
            } else {
                Pstar[Index]->StromgrenRadius = NONE;
            }
        } else if(AP[i].Type == CS_TypeSN){
        }
#ifdef USE_RADIATION_PRESSURE //{
        else if(AP[i].Type == CS_TypeRP){
            int Index = AP[i].Index;
            Pstar[Index]->RadiusRP = AP[i].Radius;
        }
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
        else if(AP[i].Type == CS_TypeSW){
            int Index = AP[i].Index;
            Pstar[Index]->RadiusSW = AP[i].Radius;
        }
#endif // USE_STELLAR_WIND //}
    }

    return ;
}

int CSEntryNumbers[CS_TypeNumber];
int CSEntryOffset[CS_TypeNumber];
int GlobalCSEntryNumbers[CS_TypeNumber];
int ReturnCalcSizeElementNumber(const int Type, const bool Global){
    if(Global == true){
        return GlobalCSEntryNumbers[Type];
    } else {
        return CSEntryNumbers[Type];
    }
}

int ReturnCalcSizeOffset(const int Type){
    return CSEntryOffset[Type];
}


void CalcSizeGetHydroInfo_i(const int Index, double *PotentialMin, double *DensityMax, double VCOM[], bool *NoLocalSink){

#ifdef USE_SINK_PARTICLE //{
    *PotentialMin = ActiveParticle[Index].Body.Hydro.PotentialMin;
    *DensityMax = ActiveParticle[Index].Body.Hydro.DensityMax;
    
    VCOM[0] = ActiveParticle[Index].Body.Hydro.VCOM[0];
    VCOM[1] = ActiveParticle[Index].Body.Hydro.VCOM[1];
    VCOM[2] = ActiveParticle[Index].Body.Hydro.VCOM[2];
#ifdef USE_SINK_NOLOCALSINK_CONDITION //{
    *NoLocalSink = ActiveParticle[Index].Body.Hydro.NoLocalSink;
#else // USE_SINK_NOLOCALSINK_CONDITION //}//{
    *NoLocalSink = true;
#endif // USE_SINK_NOLOCALSINK_CONDITION //}
#endif // USE_SINK_PARTICLE //}

    return ;
}

void CalcSizeSetSNInfo(struct StructActiveStellarFeedbackParticle ActiveStellarFeedbackParticle[]){

    int Offset = CSEntryOffset[CS_TypeSN];
    for(int i=0;i<CSEntryNumbers[CS_TypeSN];i++){
        int Index = Offset+i;
        ActiveStellarFeedbackParticle[i].Index = ActiveParticle[Index].Index;

        ActiveStellarFeedbackParticle[i].Pos[0] = ActiveParticle[Index].Pos[0];
        ActiveStellarFeedbackParticle[i].Pos[1] = ActiveParticle[Index].Pos[1];
        ActiveStellarFeedbackParticle[i].Pos[2] = ActiveParticle[Index].Pos[2];
        ActiveStellarFeedbackParticle[i].Radius = ActiveParticle[Index].Radius;
        ActiveStellarFeedbackParticle[i].WeightSum = ActiveParticle[Index].Body.SN.WeightSum;

        ActiveStellarFeedbackParticle[i].Type = ActiveParticle[Index].Body.SN.Type;
        ActiveStellarFeedbackParticle[i].Count = ActiveParticle[Index].Body.SN.Count;
        ActiveStellarFeedbackParticle[i].InitialMass = ActiveParticle[Index].Body.SN.InitialMass;
        ActiveStellarFeedbackParticle[i].Metallicity = ActiveParticle[Index].Body.SN.Metallicity;
        ActiveStellarFeedbackParticle[i].IterationCount = 0;

        ActiveStellarFeedbackParticle[i].Nlist = ActiveParticle[Index].Nlist;
        ActiveStellarFeedbackParticle[i].Density = ActiveParticle[Index].Body.SN.Density;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        ActiveStellarFeedbackParticle[i].SmoothedNumber = ActiveParticle[Index].SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
        ActiveStellarFeedbackParticle[i].GasMass = ActiveParticle[Index].Body.SN.GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
        ActiveStellarFeedbackParticle[i].DistanceMin = ActiveParticle[Index].Body.SN.DistanceMin;
        ActiveStellarFeedbackParticle[i].DistanceMinGlobalID = ActiveParticle[Index].Body.SN.DistanceMinGlobalID;
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef USE_MOMENTUM_FEEDBACK //{
        ActiveStellarFeedbackParticle[i].nH_ave = ActiveParticle[Index].Body.SN.nH_ave;
        ActiveStellarFeedbackParticle[i].Z_ave = ActiveParticle[Index].Body.SN.Z_ave;
        ActiveStellarFeedbackParticle[i].NumberDensity = ActiveParticle[Index].Body.SN.NumberDensity;
#endif // USE_MOMENTUM_FEEDBACK //}
#ifdef __CHECK_SUM__ //{
        ActiveStellarFeedbackParticle[i].CheckSum = ActiveParticle[Index].Body.SN.CheckSum;
#endif //__CHECK_SUM__ //}
#ifdef USE_STAR_TIMESTEP_LIMITER //{
        ActiveStellarFeedbackParticle[i].k = ActiveParticle[Index].k;
#endif // USE_STAR_TIMESTEP_LIMITER //}
    }

    return ;
}

#ifdef USE_RADIATION_PRESSURE //{
void CalcSizeSetRPInfo(struct StructActiveStellarFeedbackParticle ActiveStellarFeedbackParticle[]){

    int Offset = CSEntryOffset[CS_TypeRP];
    for(int i=0;i<CSEntryNumbers[CS_TypeRP];i++){
        int Index = Offset+i;
        ActiveStellarFeedbackParticle[i].Index = ActiveParticle[Index].Index;

        ActiveStellarFeedbackParticle[i].Pos[0] = ActiveParticle[Index].Pos[0];
        ActiveStellarFeedbackParticle[i].Pos[1] = ActiveParticle[Index].Pos[1];
        ActiveStellarFeedbackParticle[i].Pos[2] = ActiveParticle[Index].Pos[2];
        ActiveStellarFeedbackParticle[i].Radius = ActiveParticle[Index].Radius;
        ActiveStellarFeedbackParticle[i].WeightSum = ActiveParticle[Index].Body.RP.WeightSum;

        ActiveStellarFeedbackParticle[i].Type = StellarFeedbackType_RP;
        ActiveStellarFeedbackParticle[i].Metallicity = 
            ActiveParticle[Index].Body.RP.MetalMass/ActiveParticle[Index].Body.RP.GasMass;
        ActiveStellarFeedbackParticle[i].IterationCount = 0;

        ActiveStellarFeedbackParticle[i].Nlist = ActiveParticle[Index].Nlist;
#if 0
        fprintf(stderr,"## %d %g %g %g\n",ActiveStellarFeedbackParticle[i].Nlist,
                ActiveStellarFeedbackParticle[i].Radius,
                ActiveStellarFeedbackParticle[i].WeightSum,
                ActiveStellarFeedbackParticle[i].Metallicity);
#endif
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        ActiveStellarFeedbackParticle[i].SmoothedNumber = ActiveParticle[Index].SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef USE_STAR_TIMESTEP_LIMITER //{
        ActiveStellarFeedbackParticle[i].k = ActiveParticle[Index].k;
#endif // USE_STAR_TIMESTEP_LIMITER //}
    }

    return ;
}
#endif // USE_RADIATION_PRESSURE //}

#ifdef USE_STELLAR_WIND //{
void CalcSizeSetSWInfo(struct StructActiveStellarFeedbackParticle ActiveStellarFeedbackParticle[]){

    int Offset = CSEntryOffset[CS_TypeSW];
    for(int i=0;i<CSEntryNumbers[CS_TypeSW];i++){
        int Index = Offset+i;
        ActiveStellarFeedbackParticle[i].Index = ActiveParticle[Index].Index;

        ActiveStellarFeedbackParticle[i].Pos[0] = ActiveParticle[Index].Pos[0];
        ActiveStellarFeedbackParticle[i].Pos[1] = ActiveParticle[Index].Pos[1];
        ActiveStellarFeedbackParticle[i].Pos[2] = ActiveParticle[Index].Pos[2];
        ActiveStellarFeedbackParticle[i].Radius = ActiveParticle[Index].Radius;
        ActiveStellarFeedbackParticle[i].WeightSum = ActiveParticle[Index].Body.SW.WeightSum;

        ActiveStellarFeedbackParticle[i].Nlist = ActiveParticle[Index].Nlist;
        ActiveStellarFeedbackParticle[i].Type = StellarFeedbackType_SW;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
        ActiveStellarFeedbackParticle[i].SmoothedNumber = ActiveParticle[Index].SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef MAXIMUM_ENERGY_INPUT //{
        ActiveStellarFeedbackParticle[i].DistanceMin = ActiveParticle[Index].Body.SW.DistanceMin;
        ActiveStellarFeedbackParticle[i].DistanceMinGlobalID = ActiveParticle[Index].Body.SW.DistanceMinGlobalID;
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef USE_MOMENTUM_FEEDBACK //{
        ActiveStellarFeedbackParticle[i].nH_ave = ActiveParticle[Index].Body.SW.nH_ave;
        ActiveStellarFeedbackParticle[i].Z_ave = ActiveParticle[Index].Body.SW.Z_ave;
        ActiveStellarFeedbackParticle[i].NumberDensity = ActiveParticle[Index].Body.SW.NumberDensity;
#endif // USE_MOMENTUM_FEEDBACK //}
#ifdef USE_STAR_TIMESTEP_LIMITER //{
        ActiveStellarFeedbackParticle[i].k = ActiveParticle[Index].k;
#endif // USE_STAR_TIMESTEP_LIMITER //}
    }

    return ;
}
#endif // USE_STELLAR_WIND //}

static void FinalProcedure(const int NActives, struct StructActiveParticle AP[restrict]){

    for(int i=0;i<NActives;i++){
        if(AP[i].Type == CS_TypeHydro){
#ifdef USE_SINK_PARTICLE //{
            double InvMass = 1.0/AP[i].Body.Hydro.MassTotal;
            AP[i].Body.Hydro.VCOM[0] *= InvMass;
            AP[i].Body.Hydro.VCOM[1] *= InvMass;
            AP[i].Body.Hydro.VCOM[2] *= InvMass;

            int index = AP[i].Index;
#ifdef USE_SMOOTHED_POTENTIAL //{
            double Pot = Phydro[index]->Pot;
#else // USE_SMOOTHED_POTENTIAL //}//{
            double Pot = PhydroBody(index)->Pot/PhydroBody(index)->Mass;
#endif // USE_SMOOTHED_POTENTIAL //}
            if(Pot == AP[i].Body.Hydro.PotentialMin){
                AP[i].Body.Hydro.PotentialMin = 1;
            } else {
                AP[i].Body.Hydro.PotentialMin = 0;
            }
            if(Phydro[index]->RhoPred == AP[i].Body.Hydro.DensityMax){
                AP[i].Body.Hydro.DensityMax = 1;
            } else {
                AP[i].Body.Hydro.DensityMax = 0;
            }
#endif // USE_SINK_PARTICLE //}
        } else if(AP[i].Type == CS_TypeHII){
        } else if(AP[i].Type == CS_TypeSN){
        }
#ifdef USE_RADIATION_PRESSURE //{
        else if(AP[i].Type == CS_TypeRP){
        }
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
        else if(AP[i].Type == CS_TypeSW){
        }
#endif // USE_STELLAR_WIND //}
    }

    return ;
}


static inline void TooManyIterationReportHydro(const int Index, const int Niteration){

    if((Niteration > 100)&&(Niteration%10) == 0){
        fprintf(stderr,"-- %lld | %d | %g %g %g | %1.15g %1.15g <> %1.15g \n",ActiveParticle[Index].GlobalID,
                ActiveParticle[Index].Nlist,
                ActiveParticle[Index].Pos[0],ActiveParticle[Index].Pos[1],ActiveParticle[Index].Pos[2],
                ActiveParticle[Index].Radius,ActiveParticle[Index].Rvalue,ActiveParticle[Index].Lvalue);

        if(ActiveParticle[Index].Nlist == 0){
            int _index = ActiveParticle[Index].Index;
            StructureReportPhydro(_index);
            fflush(NULL);
        }

        fflush(NULL);
        if(isnan(ActiveParticle[Index].Pos[0])){
            int _index = ActiveParticle[Index].Index;
            fprintf(stderr,"-- %lld %d detect nan hydro | Rho %g T %g h %g\n",ActiveParticle[Index].GlobalID,_index,
                    Pall.ConvertNumberDensityToCGS*Phydro[_index]->RhoPred,
                    Pall.ConvertUtoT*Phydro[_index]->UPred,
                    Phydro[_index]->Kernel);
            StructureReportPhydro(_index);
            fflush(NULL);
            assert(1==2);
        }
        if(isinf(ActiveParticle[Index].Pos[0])){
            int _index = ActiveParticle[Index].Index;
            fprintf(stderr,"-- %lld %d detect nan hydro | Rho %g T %g h %g\n",ActiveParticle[Index].GlobalID,_index,
                    Pall.ConvertNumberDensityToCGS*Phydro[_index]->RhoPred,
                    Pall.ConvertUtoT*Phydro[_index]->UPred,
                    Phydro[_index]->Kernel);
            StructureReportPhydro(_index);
            fflush(NULL);
            assert(1==2);
        }
    }

    return ;
}

static inline void TooManyIterationReportHII(const int Index, const int Niteration){

    if((Niteration > 20)&&(Niteration%10) == 0){
        fprintf(stderr,"Too many iteration(HII): %d %d %1.15g %1.15g | %1.15g\n",Index,
                ActiveParticle[Index].Nlist,
                ActiveParticle[Index].Lvalue*Pall.UnitLength/PC_CGS,
                ActiveParticle[Index].Rvalue*Pall.UnitLength/PC_CGS,
                ActiveParticle[Index].Radius*Pall.UnitLength/PC_CGS);
#if 0
        if(ActiveParticle[Index].Rvalue-ActiveParticle[Index].Lvalue < 1.e-3*ActiveParticle[Index].Rvalue){
            ActiveParticle[Index].Lvalue *= 0.5;
            ActiveParticle[Index].Rvalue *= 2;
        }
#endif
#if 0
        fprintf(stderr,"Too many iteration(HII): %d %d %1.15g %1.15g | %1.15g\n",Index,
                ActiveParticle[Index].Nlist,
                ActiveParticle[Index].Lvalue*Pall.UnitLength/PC_CGS,
                ActiveParticle[Index].Rvalue*Pall.UnitLength/PC_CGS,
                ActiveParticle[Index].Radius*Pall.UnitLength/PC_CGS);
#endif
    }

#if 0
    else if(Niteration > 20){
        fprintf(stderr,"Too many iteration(HII): %d %1.15g %1.15g | %1.15g\n",
                ActiveParticle[Index].Nlist,
                ActiveParticle[Index].Lvalue*Pall.UnitLength/PC_CGS,
                ActiveParticle[Index].Rvalue*Pall.UnitLength/PC_CGS,
                ActiveParticle[Index].Radius*Pall.UnitLength/PC_CGS);
    }
#endif
    return ;
}


static inline void DataFullCheckHydro(const int Index){

    if(isnan(ActiveParticle[Index].Pos[0])){
        int _index = ActiveParticle[Index].Index;
        fprintf(stderr,"-- %lld | %d | %g %g %g | %g %g <> %g \n",ActiveParticle[Index].GlobalID,
                ActiveParticle[Index].Nlist,
                ActiveParticle[Index].Pos[0],ActiveParticle[Index].Pos[1],ActiveParticle[Index].Pos[2],
                ActiveParticle[Index].Radius,ActiveParticle[Index].Rvalue,ActiveParticle[Index].Lvalue);
        fprintf(stderr,"-- %lld %d detect nan hydro | Rho %g T %g h %g\n",ActiveParticle[Index].GlobalID,_index,
                Pall.ConvertNumberDensityToCGS*Phydro[_index]->RhoPred,
                Pall.ConvertUtoT*Phydro[_index]->UPred,
                Phydro[_index]->Kernel);
        fflush(NULL);
        assert(1==2);
    }
    if(isinf(ActiveParticle[Index].Pos[0])){
        int _index = ActiveParticle[Index].Index;
        fprintf(stderr,"-- %lld | %d | %g %g %g | %g %g <> %g \n",ActiveParticle[Index].GlobalID,
                ActiveParticle[Index].Nlist,
                ActiveParticle[Index].Pos[0],ActiveParticle[Index].Pos[1],ActiveParticle[Index].Pos[2],
                ActiveParticle[Index].Radius,ActiveParticle[Index].Rvalue,ActiveParticle[Index].Lvalue);
        fprintf(stderr,"-- %lld %d detect nan hydro | Rho %g T %g h %g\n",ActiveParticle[Index].GlobalID,_index,
                Pall.ConvertNumberDensityToCGS*Phydro[_index]->RhoPred,
                Pall.ConvertUtoT*Phydro[_index]->UPred,
                Phydro[_index]->Kernel);
        fflush(NULL);
        assert(1==2);
    }
   

    return ;
}

static inline void DataFullCheckHII(const int Index){

    if(isnan(ActiveParticle[Index].Pos[0])){
        int _index = ActiveParticle[Index].Index;
        fprintf(stderr,"-- %lld | %d | %g %g %g | %g %g <> %g \n",ActiveParticle[Index].GlobalID,
                ActiveParticle[Index].Nlist,
                ActiveParticle[Index].Pos[0],ActiveParticle[Index].Pos[1],ActiveParticle[Index].Pos[2],
                ActiveParticle[Index].Radius,ActiveParticle[Index].Rvalue,ActiveParticle[Index].Lvalue);
        fprintf(stderr,"-- %lld %d detect nan HII | Formation time %g\n",ActiveParticle[Index].GlobalID,_index,
                Pstar[_index]->FormationTime);
        fflush(NULL);
        assert(1==2);
    }
    if(isinf(ActiveParticle[Index].Pos[0])){
        int _index = ActiveParticle[Index].Index;
        fprintf(stderr,"-- %lld | %d | %g %g %g | %g %g <> %g \n",ActiveParticle[Index].GlobalID,
                ActiveParticle[Index].Nlist,
                ActiveParticle[Index].Pos[0],ActiveParticle[Index].Pos[1],ActiveParticle[Index].Pos[2],
                ActiveParticle[Index].Radius,ActiveParticle[Index].Rvalue,ActiveParticle[Index].Lvalue);
        fprintf(stderr,"-- %lld %d detect nan HII | Formation time %g\n",ActiveParticle[Index].GlobalID,_index,
                Pstar[_index]->FormationTime);
        fflush(NULL);
        assert(1==2);
    }

    return ;
}

static inline void DataFullCheckSNe(const int Index){

    if(isnan(ActiveParticle[Index].Pos[0])){
        int _index = ActiveParticle[Index].Index;
        fprintf(stderr,"-- %lld | %d | %g %g %g | %g %g <> %g \n",ActiveParticle[Index].GlobalID,
                ActiveParticle[Index].Nlist,
                ActiveParticle[Index].Pos[0],ActiveParticle[Index].Pos[1],ActiveParticle[Index].Pos[2],
                ActiveParticle[Index].Radius,ActiveParticle[Index].Rvalue,ActiveParticle[Index].Lvalue);
        fprintf(stderr,"-- %lld %d detect nan SNe | Formation time %g %d\n",ActiveParticle[Index].GlobalID,_index,
                Pstar[_index]->FormationTime, ActiveParticle[Index].Body.SN.Type);
        fflush(NULL);
        assert(1==2);
    }
    if(isinf(ActiveParticle[Index].Pos[0])){
        int _index = ActiveParticle[Index].Index;
        fprintf(stderr,"-- %lld | %d | %g %g %g | %g %g <> %g \n",ActiveParticle[Index].GlobalID,
                ActiveParticle[Index].Nlist,
                ActiveParticle[Index].Pos[0],ActiveParticle[Index].Pos[1],ActiveParticle[Index].Pos[2],
                ActiveParticle[Index].Radius,ActiveParticle[Index].Rvalue,ActiveParticle[Index].Lvalue);
        fprintf(stderr,"-- %lld %d detect nan SNe | Formation time %g %d\n",ActiveParticle[Index].GlobalID,_index,
                Pstar[_index]->FormationTime,ActiveParticle[Index].Body.SN.Type);
        fflush(NULL);
        assert(1==2);
    }

    return ;
}

static inline int ExchangeParticleNumbersForExport(int NExportThisTime[], int NImportThisTime[]){

    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();

    int NImportThisTime2[NProcs];
    int NExportThisTime2[NProcs];
    NImportThisTime2[MyID] = 0;
    NExportThisTime2[MyID] = 0;
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime2[CommunicationTable[i].SendRank] = NExportThisTime[i];
    }
    MPI_Alltoall(NExportThisTime2,1,MPI_INT,NImportThisTime2,1,MPI_INT,MPI_COMM_WORLD);
    int NImport = 0;
    for(int i=0;i<NProcs-1;i++){
        NImportThisTime[i] = NImportThisTime2[CommunicationTable[i].RecvRank];
        NImport += NImportThisTime[i];
    }
    return NImport;
}


static bool first = true;
void CalcSize(void){

    double TimingResultThisRoutine = GetElapsedTime();

#ifdef EVALUATE_KERNEL_BY_ITERATION
    if(first){
        CS_AllocateContactedDomainID();
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
        SmoothedMassConversionFactor = (4.0*M_PI/3.0)*8.0;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER
#ifdef USE_NEIGHBOR_LIST //{
        LocalNeighborListStride = Pall.Ns+Pall.Npm+1;
#endif // USE_NEIGHBOR_LIST //}
#ifdef USE_NEIGHBOR_LIST_AND_SORT //{
        NnbSort = NEIGHBOR_LIST_SORT_LIMIT_FACTOR*(Pall.Ns+Pall.Npm);
#endif // USE_NEIGHBOR_LIST_AND_SORT //}
        first = false;
    }

    OverMaxIterationTimes = false;

#ifdef USE_BARYON_COMM //}
    if(MPI_BARYON_COMM_WORLD == MPI_COMM_NULL){
        return ;
    }
    int MyID = MPIGetBaryonMyID();
    int NProcs = MPIGetBaryonNumProcs();
#else // USE_BARYON_COMM //}//{
    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();
#endif // USE_BARYON_COMM //}
    int Neighbors[MaxNeighborSize];

    static bool *CSExportFlags = NULL;


    int NLocalActiveYoungStars = CheckHIIflags(0);
#ifdef USE_CELIB //{
    int NLocalActiveFeedback = CountFeedbackNumber();
#else // USE_CELIB //}//{
    int NLocalActiveFeedback = 0;
#endif //USE_CELIB //}

#ifdef USE_RADIATION_PRESSURE //{
    int NLocalActiveRadiationPressure = CountRadiationPressureSourceNumber();
#else // USE_RADIATION_PRESSURE //}//{
    int NLocalActiveRadiationPressure = 0;
#endif // USE_RADIATION_PRESSURE //}

#ifdef USE_STELLAR_WIND //{
    int NLocalActiveStellarWind = CountStellarWindSourceNumber();
#else // USE_STELLAR_WIND //}//{
    int NLocalActiveStellarWind = 0;
#endif // USE_STELLAR_WIND //}

    int MaxEntry = MAX(Pall.Nhydro+(NLocalActiveYoungStars+NLocalActiveFeedback
                +NLocalActiveRadiationPressure+NLocalActiveStellarWind),NAdditionUnit);

    if(CSExportFlagsMaxAllocated < MaxEntry){
        if(CSExportFlagsMaxAllocated > 0){
            free(CSExportFlags);
            free(ActiveParticle);
#ifdef USE_NEIGHBOR_LIST //{
            free(LocalNeighborList);
            free(LNLK);
#endif // USE_NEIGHBOR_LIST //}
        }
        CSExportFlagsMaxAllocated = MAX(ForAngelsShare*MaxEntry,NAdditionUnit);
        CSExportFlags = malloc(sizeof(bool)*CSExportFlagsMaxAllocated*NProcs);
        ActiveParticle = malloc(sizeof(struct StructActiveParticle)*CSExportFlagsMaxAllocated);
#ifdef USE_NEIGHBOR_LIST //{
        LocalNeighborList = malloc(sizeof(int)*CSExportFlagsMaxAllocated*LocalNeighborListStride);
        LNLK = malloc(sizeof(double)*CSExportFlagsMaxAllocated);
#endif // USE_NEIGHBOR_LIST //}
    }

    for(int i=0;i<CS_TypeNumber;i++){
        CSEntryNumbers[i] = 0; 
    }

    int NActives = 0;
    CSEntryNumbers[CS_TypeHydro] = HydroEntry(ActiveParticle,CSExportFlags,NProcs);
    CSEntryOffset[CS_TypeHydro] = NActives;
    NActives += CSEntryNumbers[CS_TypeHydro];

    CSEntryNumbers[CS_TypeSN] = StellarFeedbackEntry(ActiveParticle+NActives,CSExportFlags+NActives*NProcs,NProcs);
    CSEntryOffset[CS_TypeSN] = NActives;
    NActives += CSEntryNumbers[CS_TypeSN];

    CSEntryNumbers[CS_TypeHII] += HIIregionEntry(ActiveParticle+NActives,CSExportFlags+NActives*NProcs,NProcs);
    CSEntryOffset[CS_TypeHII] = NActives;
    NActives += CSEntryNumbers[CS_TypeHII];

#ifdef USE_RADIATION_PRESSURE //{
    CSEntryNumbers[CS_TypeRP] += RadiationPressureEntry(ActiveParticle+NActives,CSExportFlags+NActives*NProcs,NProcs);
    CSEntryOffset[CS_TypeRP] = NActives;
    NActives += CSEntryNumbers[CS_TypeRP];
#endif // USE_RADIATION_PRESSURE //}

#ifdef USE_STELLAR_WIND //{
    CSEntryNumbers[CS_TypeSW] += StellarWindEntry(ActiveParticle+NActives,CSExportFlags+NActives*NProcs,NProcs);
    CSEntryOffset[CS_TypeSW] = NActives;
    NActives += CSEntryNumbers[CS_TypeSW];
#endif // USE_STELLAR_WIND //}

#ifdef USE_BARYON_COMM //}
    MPI_Allreduce(CSEntryNumbers,GlobalCSEntryNumbers,CS_TypeNumber,MPI_INT,MPI_SUM,MPI_BARYON_COMM_WORLD);
#else // USE_BARYON_COMM //}//{
    MPI_Allreduce(CSEntryNumbers,GlobalCSEntryNumbers,CS_TypeNumber,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif // USE_BARYON_COMM //}


    const int BitMask = 0x01; 
    int NExportThisTime[NProcs-1];
    int NImportThisTime[NProcs-1];
    int NExportThisTimeNew[NProcs];
    int NImportThisTimeNew[NProcs];

    struct StructCSExport *CSExportSend[NProcs-1];
    struct StructCSExport *CSExportRecv = NULL;
    struct StructCSImport *CSImportSend = NULL;
    struct StructCSImport *CSImportRecv[NProcs-1];
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

#ifdef PRINT_LOG_KERNEL_ITERATION
    if(MyID==MPI_ROOT_RANK)
        fprintf(stderr,"Kernel Iteration");
#endif // PRINT_LOG_KERNEL_ITERATION

    int NActiveLeaves;
#ifdef USE_BARYON_COMM //}
    MPI_Allreduce(&NActives,&NActiveLeaves,1,MPI_INT,MPI_SUM,MPI_BARYON_COMM_WORLD);
#else // USE_BARYON_COMM //}//{
    MPI_Allreduce(&NActives,&NActiveLeaves,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif // USE_BARYON_COMM //}

    ///////////////////////////////////////////////////////////////
    // Baryon domain & hydro domain 
    // 
    CS_CheckContactedDomain();
    double BoxCenter[] = {HydroNode[0].Pos[0],HydroNode[0].Pos[1],HydroNode[0].Pos[2]}; // need check
    //double BoxCenter[] = {GravityNode[0].Pos[0],GravityNode[0].Pos[1],GravityNode[0].Pos[2]};

    int Niteration = 0;
    Niteration_for_debug = Niteration;
    do{
#ifdef PRINT_LOG_KERNEL_ITERATION
        if(MPIGetMyID()==MPI_ROOT_RANK)
            fprintf(stderr,":[%d] = %d ",Niteration,NActiveLeaves);
#endif // PRINT_LOG_KERNEL_ITERATION

        // Initialization part
        LocalKernelMax = 0.e0;  
        int Checker[CS_TypeNumber] = {0};
        int CheckerNlist[CS_TypeNumber] = {0};
        for(int i=0;i<NActives;i++){  // Clear Entries
            int Offset = i*NProcs;

            if(CSExportFlags[Offset+NProcs-1]){
#ifdef USE_NEIGHBOR_LIST //{
                int LNLOffset = i*LocalNeighborListStride;
                LocalNeighborList[LNLOffset+LocalNeighborListStride-1] = NONE;
#endif // USE_NEIGHBOR_LIST //}

                ActiveParticle[i].k_hydro_localmin = MAXIMUM_TIME_HIERARCHY;

                if(ActiveParticle[i].Type == CS_TypeHydro){
                    CheckerNlist[CS_TypeHydro] += ActiveParticle[i].Nlist;
                } else if(ActiveParticle[i].Type == CS_TypeHII){
                    CheckerNlist[CS_TypeHII] += ActiveParticle[i].Nlist;
                } else if(ActiveParticle[i].Type == CS_TypeSN){
                    CheckerNlist[CS_TypeSN] += ActiveParticle[i].Nlist;
                } 
#ifdef USE_RADIATION_PRESSURE //{
                else if(ActiveParticle[i].Type == CS_TypeRP){
                    CheckerNlist[CS_TypeRP] += ActiveParticle[i].Nlist;
                }
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
                else if(ActiveParticle[i].Type == CS_TypeSW){
                    CheckerNlist[CS_TypeSW] += ActiveParticle[i].Nlist;
                }
#endif // USE_STELLAR_WIND //}
 
                int oldNlist = ActiveParticle[i].Nlist;
                ActiveParticle[i].Nlist = 0;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                ActiveParticle[i].SmoothedNumber = 0.e0;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
                for(int k=0;k<NProcs-1;k++)
                    CSExportFlags[Offset+k] = false;
                LocalKernelMax = fmax(LocalKernelMax,DISTANCE(BoxCenter,ActiveParticle[i].Pos)+2.0*ActiveParticle[i].Radius);
                if(ActiveParticle[i].Type == CS_TypeHydro){
                    //DataFullCheckHydro(i);
                    TooManyIterationReportHydro(i,Niteration);

#ifdef USE_SINK_PARTICLE //{
                    ActiveParticle[i].Body.Hydro.PotentialMin = 0.e0;
                    ActiveParticle[i].Body.Hydro.DensityMax = 0.e0;
                    ActiveParticle[i].Body.Hydro.MassTotal = 0.e0;
                    ActiveParticle[i].Body.Hydro.VCOM[0] = 
                    ActiveParticle[i].Body.Hydro.VCOM[1] = 
                    ActiveParticle[i].Body.Hydro.VCOM[2] = 0.e0;
#ifdef USE_SINK_NOLOCALSINK_CONDITION //{
                    ActiveParticle[i].Body.Hydro.NoLocalSink = true;
#endif // USE_SINK_NOLOCALSINK_CONDITION //}
#endif // USE_SINK_PARTICLE //}
                    Checker[CS_TypeHydro] ++;
                } else if(ActiveParticle[i].Type == CS_TypeHII){
                    //DataFullCheckHII(i);
                    ActiveParticle[i].Nlist = oldNlist;
                    TooManyIterationReportHII(i,Niteration);
                    ActiveParticle[i].Nlist = 0;

                    ActiveParticle[i].Body.HII.DistanceMin = 0.e0;
#ifdef __PHOTON_COUNT_BASE__
                    ActiveParticle[i].Body.HII.PhotonCount = ActiveParticle[i].Body.HII.PhotonCountDistanceMin = 0.e0;
#else //__PHOTON_COUNT_BASE__
                    ActiveParticle[i].Body.HII.Mass = ActiveParticle[i].Body.HII.MassDistanceMin = 0.e0;
#endif // __PHOTON_COUNT_BASE__
                    Checker[CS_TypeHII] ++;
                } else if(ActiveParticle[i].Type == CS_TypeSN){
                    //DataFullCheckSNe(i);
                    ActiveParticle[i].Body.SN.Density = 0.0;
                    ActiveParticle[i].Body.SN.WeightSum = 0.0;
#ifdef SET_SNII_TEMPERATURE //{
                    ActiveParticle[i].Body.SN.GasMass = 0.e0;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
                    ActiveParticle[i].Body.SN.DistanceMin = 2.0*ActiveParticle[i].Radius;
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef USE_MOMENTUM_FEEDBACK //{
                    ActiveParticle[i].Body.SN.nH_ave = 0.e0;
                    ActiveParticle[i].Body.SN.Z_ave  = 0.e0;
                    ActiveParticle[i].Body.SN.NumberDensity = 0.e0;
#endif // USE_MOMENTUM_FEEDBACK //}
                    Checker[CS_TypeSN] ++;
                }
#ifdef USE_RADIATION_PRESSURE //{
                else if(ActiveParticle[i].Type == CS_TypeRP){
                    ActiveParticle[i].Body.RP.GasMass = 0.0;
                    ActiveParticle[i].Body.RP.MetalMass = 0.0;
                    ActiveParticle[i].Body.RP.WeightSum = 0.0;
                    Checker[CS_TypeRP] ++;
                }
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
                else if(ActiveParticle[i].Type == CS_TypeSW){
                    ActiveParticle[i].Body.SW.WeightSum = 0.0;
#ifdef MAXIMUM_ENERGY_INPUT //{
                    ActiveParticle[i].Body.SW.DistanceMin = 2.0*ActiveParticle[i].Radius;
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef USE_MOMENTUM_FEEDBACK //{
                    ActiveParticle[i].Body.SW.nH_ave = 0.e0;
                    ActiveParticle[i].Body.SW.Z_ave  = 0.e0;
                    ActiveParticle[i].Body.SW.NumberDensity = 0.e0;
#endif // USE_MOMENTUM_FEEDBACK //}
                    Checker[CS_TypeSW] ++;
                }
#endif // USE_STELLAR_WIND //}
            }
        }

        if((Niteration > 100)&&(Niteration%10) == 0){
            MPI_Allreduce(MPI_IN_PLACE,Checker,CS_TypeNumber,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE,CheckerNlist,CS_TypeNumber,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

            if(MPIGetMyID() == MPI_ROOT_RANK){
                fprintf(stderr,"Iteration %d: %d %d %d\n",Niteration,Checker[0],Checker[1],Checker[2]);
                fprintf(stderr,"Nlist %d: %d %d %d | %d %d %d \n",Niteration,
                        CheckerNlist[0],CheckerNlist[1],CheckerNlist[2],
                        CheckerNlist[0]/MAX(Checker[0],1),
                        CheckerNlist[1]/MAX(Checker[1],1),
                        CheckerNlist[2]/MAX(Checker[2],1));
                fflush(NULL);
            }
        }


        int NExportMaxThisTime = 0;
        for(int i=0;i<NProcs-1;i++){
            NExportThisTime[i] = CS_CheckExportFlags(i,NProcs,CSExportFlags,NActives,ActiveParticle);

            CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructCSExport),i);
            CheckSizeofBufferImportRecvIndex(NExportThisTime[i],sizeof(struct StructCSImport),i);
            CSExportSend[i] = BufferExportSend[i];
            CSImportRecv[i] = BufferImportRecv[i];

            int NExport = 0;
            if(NExportThisTime[i] > 0){
                for(int k=0;k<NActives;k++){
                    int Offset = k*NProcs;
                    if(CSExportFlags[Offset+NProcs-1]){ 
                        if(CSExportFlags[Offset+i]&BitMask){ 
                            CSExportSend[i][NExport].Type = ActiveParticle[k].Type;
                            CSExportSend[i][NExport].Pos[0] = ActiveParticle[k].Pos[0];
                            CSExportSend[i][NExport].Pos[1] = ActiveParticle[k].Pos[1];
                            CSExportSend[i][NExport].Pos[2] = ActiveParticle[k].Pos[2];
                            CSExportSend[i][NExport].Radius = ActiveParticle[k].Radius;
                            CSExportSend[i][NExport].Leaf = k;
#ifdef USE_DEBUG_MODE //{
                            CSExportSend[i][NExport].GlobalID = ActiveParticle[k].GlobalID;
#endif // USE_DEBUG_MODE //}
                            NExport ++;
                        }
                    }
                }
            }
            NExportThisTime[i] = NExport;
            NExportMaxThisTime = MAX(NExportMaxThisTime,NExport);
        }

        //int NImportAll = ExchangeParticleNumbersForExport(NExportThisTime,NImportThisTime);
#if 1
        int NImportThisTime2[NProcs];
        int NExportThisTime2[NProcs];
        NImportThisTime2[MyID] = 0;
        NExportThisTime2[MyID] = 0;
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
#endif


        CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructCSExport));
        CheckSizeofBufferImportSend(NImportAll,sizeof(struct StructCSImport));
        CSExportRecv = BufferExportRecv;
        CSImportSend = BufferImportSend; 

        NImport = 0;

        int counter_send = 0;
        int counter_recv = 0;

        int SendFlag,RecvFlag;
        for(int i=0;i<NProcs-1;i++){
            if(NExportThisTime[i]>0){
                MPI_Isend(CSExportSend[i],
                    NExportThisTime[i]*sizeof(struct StructCSExport),
                        MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_KERNEL_EXPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
                MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
                counter_send ++;
            }
            if(NImportThisTime[i]>0){
                MPI_Irecv(CSExportRecv+NImport,
                    NImportThisTime[i]*sizeof(struct StructCSExport),
                        MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_KERNEL_EXPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
                MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
                counter_recv ++;
                NImport += NImportThisTime[i];
            }
        }


        double TimeNBS = GetElapsedTime();
        for(int i=0;i<NActives;i++){  // Check local
            int Offset = i*NProcs;
            if(CSExportFlags[Offset+NProcs-1]){
                CS_UpdateRadiusAndOthers_i(ActiveParticle+i,Neighbors);

#ifdef UPDATE_SIZE_LOCAL //{
                /// Insert Local Update Routine here.
                ActiveParticle[i].LocalUpdateFlags = false;
                int IsLocal = 0;
                for(int k=0;k<NProcs-1;k++){
                    if(CSExportFlags[Offset+k]&BitMask){
                        IsLocal ++;
                    }
                }

                if(IsLocal == 0){
#   ifdef USE_NEIGHBOR_LIST_AND_SORT //{
                    SortNeighborListAndSetKernelSize(ActiveParticle+i,Neighbors);
#   endif // USE_NEIGHBOR_LIST_AND_SORT //}
#   ifdef USE_NEIGHBOR_LIST //{
                    CopyNeighborList(ActiveParticle+i,Neighbors,i*LocalNeighborListStride);
#   endif // USE_NEIGHBOR_LIST //}
                    ActiveParticle[i].LocalUpdateFlags = true;
                    CS_UpdateRadiusLocal(ActiveParticle+i,MyID,NProcs,CSExportFlags+Offset,
#   ifdef USE_NEIGHBOR_LIST //{
                            i*LocalNeighborListStride);
#   else // USE_NEIGHBOR_LIST //}//{
                            0);
#   endif // USE_NEIGHBOR_LIST //}//{
                }
#endif // UPDATE_SIZE_LOCAL // }
            }
        }
        TimingResults.HydroKernelNeighborSearchThisStep += GetElapsedTime()-TimeNBS;

        double TimeComm = GetElapsedTime();
        MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
        MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);
        TimingResults.HydroKernelCommThisStep += GetElapsedTime()-TimeComm;


        TimeNBS = GetElapsedTime();
        for(int i=0;i<NImportAll;i++){ // Imported
            CS_UpdateRadiusAndOthersImported_i(CSImportSend+i,CSExportRecv+i,Neighbors);
        }
        TimingResults.HydroKernelNeighborSearchThisStep += GetElapsedTime()-TimeNBS;

#if 1
        NImportAll = 0;
        int NImportAllNew = 0;
        for(int i=0;i<NProcs-1;i++){
            NImportThisTimeNew[i] = 0;
            for(int k=0;k<NImportThisTime[i];k++){
                if(CSImportSend[NImportAll].Nlist > 0){
                    CSImportSend[NImportAllNew] = CSImportSend[NImportAll];
                    NImportThisTimeNew[i] ++;
                    NImportAllNew ++;
                }
                NImportAll ++;
            }
        }

#if 0
static inline int ExchangeParticleNumbersForImport(int NExportThisTimeNew[restrict], int NImportThisTimeNew[restrict]){

    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();

    int NImportThisTimeNew2[NProcs];
    int NExportThisTimeNew2[NProcs];
    NExportThisTimeNew2[MyID] = 0;
    NImportThisTimeNew2[MyID] = 0;
    for(int i=0;i<NProcs-1;i++){
        NImportThisTimeNew2[CommunicationTable[i].SendRank] = NImportThisTimeNew[i];
    }
    MPI_Alltoall(NImportThisTimeNew2,1,MPI_INT,NExportThisTimeNew2,1,MPI_INT,MPI_COMM_WORLD);
    for(int i=0;i<NProcs-1;i++){
        NExportThisTimeNew[i] = NExportThisTimeNew2[CommunicationTable[i].RecvRank];
    }

    return ;
}
#endif

#if 1
        int NImportThisTimeNew2[NProcs];
        int NExportThisTimeNew2[NProcs];
        NExportThisTimeNew2[MyID] = 0;
        NImportThisTimeNew2[MyID] = 0;
        for(int i=0;i<NProcs-1;i++){
            // NImportThisTimeNew2[CommunicationTable[i].SendRank] = NImportThisTimeNew[i];
             NImportThisTimeNew2[CommunicationTable[i].RecvRank] = NImportThisTimeNew[i];
        }
        MPI_Alltoall(NImportThisTimeNew2,1,MPI_INT,NExportThisTimeNew2,1,MPI_INT,MPI_COMM_WORLD);
        for(int i=0;i<NProcs-1;i++){
            //NExportThisTimeNew[i] = NExportThisTimeNew2[CommunicationTable[i].RecvRank];
            NExportThisTimeNew[i] = NExportThisTimeNew2[CommunicationTable[i].SendRank];
        }
#endif

#else
        for(int i=0;i<NProcs-1;i++){
            NExportThisTimeNew[i] = NExportThisTime[i];
            NImportThisTimeNew[i] = NImportThisTime[i];
        }
#endif


        NImport = 0;
        counter_send = counter_recv = 0;
        for(int i=0;i<NProcs-1;i++){
            if(NImportThisTimeNew[i]>0){
                MPI_Isend(CSImportSend+NImport,
                    NImportThisTimeNew[i]*sizeof(struct StructCSImport),
                        //MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_KERNEL_IMPORT+i,
                        MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_KERNEL_IMPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
                MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
                counter_send ++;
            }
            if(NExportThisTimeNew[i]>0){
                MPI_Irecv(CSImportRecv[i],
                    NExportThisTimeNew[i]*sizeof(struct StructCSImport),
                        //MPI_BYTE,CommunicationTable[i].RecvRank,TAG_SPH_KERNEL_IMPORT+i,
                        MPI_BYTE,CommunicationTable[i].SendRank,TAG_SPH_KERNEL_IMPORT+i,
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
                int leaf = CSImportRecv[i][k].Leaf;
                ActiveParticle[leaf].k_hydro_localmin = fmin(ActiveParticle[leaf].k_hydro_localmin,CSImportRecv[i][k].k_hydro_localmin);
                if(CSImportRecv[i][k].Type == CS_TypeHydro){
                    ActiveParticle[leaf].Nlist += CSImportRecv[i][k].Nlist;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                    ActiveParticle[leaf].SmoothedNumber += CSImportRecv[i][k].SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}

#ifdef USE_SINK_PARTICLE //{
                    ActiveParticle[leaf].Body.Hydro.PotentialMin = fmin(ActiveParticle[leaf].Body.Hydro.PotentialMin,CSImportRecv[i][k].Body.Hydro.PotentialMin);
                    ActiveParticle[leaf].Body.Hydro.DensityMax = fmax(ActiveParticle[leaf].Body.Hydro.DensityMax,CSImportRecv[i][k].Body.Hydro.DensityMax);
                    ActiveParticle[leaf].Body.Hydro.MassTotal += CSImportRecv[i][k].Body.Hydro.MassTotal;
                    ActiveParticle[leaf].Body.Hydro.VCOM[0] += CSImportRecv[i][k].Body.Hydro.VCOM[0];
                    ActiveParticle[leaf].Body.Hydro.VCOM[1] += CSImportRecv[i][k].Body.Hydro.VCOM[1];
                    ActiveParticle[leaf].Body.Hydro.VCOM[2] += CSImportRecv[i][k].Body.Hydro.VCOM[2];
#ifdef USE_SINK_NOLOCALSINK_CONDITION //{
                    if(!CSImportRecv[i][k].Body.Hydro.NoLocalSink){
                        ActiveParticle[leaf].Body.Hydro.NoLocalSink = false;
                    }
#endif // USE_SINK_NOLOCALSINK_CONDITION //}
#endif // USE_SINK_PARTICLE //}
                    
                } else if(CSImportRecv[i][k].Type == CS_TypeHII){
                    ActiveParticle[leaf].Nlist += CSImportRecv[i][k].Nlist;
#ifdef __PHOTON_COUNT_BASE__ //{
                    ActiveParticle[leaf].Body.HII.PhotonCount += CSImportRecv[i][k].Body.HII.PhotonCount;
#else // __PHOTON_COUNT_BASE__ //}//{
                    ActiveParticle[leaf].Body.HII.Mass += CSImportRecv[i][k].Body.HII.Mass;
#endif //__PHOTON_COUNT_BASE__ //}
                    if((ActiveParticle[leaf].Body.HII.DistanceMin > CSImportRecv[i][k].Body.HII.DistanceMin)&&(CSImportRecv[i][k].Nlist > 0)){
                        ActiveParticle[leaf].Body.HII.DistanceMin = CSImportRecv[i][k].Body.HII.DistanceMin;
#ifdef __PHOTON_COUNT_BASE__ //{
                        ActiveParticle[leaf].Body.HII.PhotonCountDistanceMin = CSImportRecv[i][k].Body.HII.PhotonCountDistanceMin;
#else // __PHOTON_COUNT_BASE__ //}//{
                        ActiveParticle[leaf].Body.HII.MassDistanceMin = CSImportRecv[i][k].MassDistanceMin;
#endif //__PHOTON_COUNT_BASE__ //}
                    }
                } else if(CSImportRecv[i][k].Type == CS_TypeSN){
                    ActiveParticle[leaf].Nlist += CSImportRecv[i][k].Nlist;
                    ActiveParticle[leaf].Body.SN.Density += CSImportRecv[i][k].Body.SN.Density;
                
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                    ActiveParticle[leaf].Body.SN.SmoothedNumber += CSImportRecv[i][k].SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef SET_SNII_TEMPERATURE //{
                    ActiveParticle[leaf].Body.SN.GasMass += CSImportRecv[i][k].Body.SN.GasMass;
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
                    if(ActiveParticle[leaf].Body.SN.DistanceMin > CSImportRecv[i][k].Body.SN.DistanceMin){
                        ActiveParticle[leaf].Body.SN.DistanceMin = CSImportRecv[i][k].Body.SN.DistanceMin;
                        ActiveParticle[leaf].Body.SN.DistanceMinGlobalID = CSImportRecv[i][k].Body.SN.DistanceMinGlobalID;
                    }
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef USE_MOMENTUM_FEEDBACK //{
                    ActiveParticle[leaf].Body.SN.nH_ave += CSImportRecv[i][k].Body.SN.nH_ave;
                    ActiveParticle[leaf].Body.SN.Z_ave += CSImportRecv[i][k].Body.SN.Z_ave;
                    ActiveParticle[leaf].Body.SN.nH_ave /= fmax((double)ActiveParticle[leaf].Nlist,1.0);
                    ActiveParticle[leaf].Body.SN.Z_ave /= fmax((double)ActiveParticle[leaf].Nlist,1.0);
                    ActiveParticle[leaf].Body.SN.NumberDensity += CSImportRecv[i][k].Body.SN.NumberDensity;
#endif //USE_MOMENTUM_FEEDBACK //}
                    ActiveParticle[leaf].Body.SN.WeightSum += CSImportRecv[i][k].Body.SN.WeightSum;
#ifdef __CHECK_SUM__ //{
                    ActiveParticle[leaf].Body.SN.CheckSum += CSImportRecv[i][k].Body.SN.CheckSum;
#endif //__CHECK_SUM__ //}
                }
#ifdef USE_RADIATION_PRESSURE //{
                else if(CSImportRecv[i][k].Type == CS_TypeRP){
                    ActiveParticle[leaf].Nlist += CSImportRecv[i][k].Nlist;
                    ActiveParticle[leaf].Body.RP.GasMass += CSImportRecv[i][k].Body.RP.GasMass;
                    ActiveParticle[leaf].Body.RP.MetalMass += CSImportRecv[i][k].Body.RP.MetalMass;
                    ActiveParticle[leaf].Body.RP.WeightSum += CSImportRecv[i][k].Body.RP.WeightSum;
                }
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
                else if(CSImportRecv[i][k].Type == CS_TypeSW){
                    ActiveParticle[leaf].Nlist += CSImportRecv[i][k].Nlist;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
                    ActiveParticle[leaf].Body.SW.SmoothedNumber += CSImportRecv[i][k].SmoothedNumber;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef MAXIMUM_ENERGY_INPUT //{
                    if(ActiveParticle[leaf].Body.SW.DistanceMin > CSImportRecv[i][k].Body.SW.DistanceMin){
                        ActiveParticle[leaf].Body.SW.DistanceMin = CSImportRecv[i][k].Body.SW.DistanceMin;
                        ActiveParticle[leaf].Body.SW.DistanceMinGlobalID = CSImportRecv[i][k].Body.SW.DistanceMinGlobalID;
                    }
#endif // MAXIMUM_ENERGY_INPUT //}
#ifdef USE_MOMENTUM_FEEDBACK //{
                    ActiveParticle[leaf].Body.SW.nH_ave += CSImportRecv[i][k].Body.SW.nH_ave;
                    ActiveParticle[leaf].Body.SW.Z_ave += CSImportRecv[i][k].Body.SW.Z_ave;
                    ActiveParticle[leaf].Body.SW.nH_ave /= fmax((double)ActiveParticle[leaf].Nlist,1.0);
                    ActiveParticle[leaf].Body.SW.Z_ave /= fmax((double)ActiveParticle[leaf].Nlist,1.0);
                    ActiveParticle[leaf].Body.SW.NumberDensity += CSImportRecv[i][k].Body.SW.NumberDensity;
#endif //USE_MOMENTUM_FEEDBACK //}
                    ActiveParticle[leaf].Body.SW.WeightSum += CSImportRecv[i][k].Body.SW.WeightSum;
                }
#endif // USE_STELLAR_WIND //}
            }
        }

#ifdef ADD_PERTURBATION //{
        // assert(Niteration < 1000);
        if(Niteration > 10*MaxIterationTimes) break;
        if(Niteration > MaxIterationTimes)
            OverMaxIterationTimes = true;
        ResetLRvalues(NActives,Niteration,NProcs,CSExportFlags,ActiveParticle);
#endif // ADD_PERTURBATION //}

        int NLocalActiveLeaves = CS_CheckConvergence(NActives,NProcs,CSExportFlags,ActiveParticle); 
        MPI_Allreduce(&NLocalActiveLeaves,&NActiveLeaves,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

        Niteration ++;
        Niteration_for_debug = Niteration;
#ifndef ADD_PERTURBATION //{
        if(Niteration > 10*MaxIterationTimes)
            break;
#endif // ADD_PERTURBATION //}
    } while (0<NActiveLeaves);

#else // EVALUATE_KERNEL_BY_ITERATION
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            Phydro[i]->Kernel = KERNEL_FACTOR*
                pow(Phydro[i]->Mass/Phydro[i]->Rho,1.0/((double)DIMENSION));
            Phydro[i]->Rho = 0.e0;
        }
    }
#endif // EVALUATE_KERNEL_BY_ITERATION

    // Restore all data for Phydro/Pstar
    RestoreAllData(NActives,ActiveParticle);

    // Final procedure for VCOM
    FinalProcedure(NActives,ActiveParticle);


    PlantHydroTreeKernelMaxUpdate();
#ifdef EVALUATE_KERNEL_BY_ITERATION
#ifdef PRINT_LOG_KERNEL_ITERATION
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"\n");
#else // PRINT_LOG_KERNEL_ITERATION
    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"%d iterations for kernel determination.\n",Niteration);
#endif // PRINT_LOG_KERNEL_ITERATION
#endif // EVALUATE_KERNEL_BY_ITERATION

    TimingResults.HydroKernelThisStep += GetElapsedTime()-TimingResultThisRoutine;

    return;
}
