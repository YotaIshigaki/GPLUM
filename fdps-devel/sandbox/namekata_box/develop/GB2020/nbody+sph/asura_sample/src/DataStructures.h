/* 
 * This file contains all of data structures in the simulation code.
 * Structure StructPbody includes global information of all type particles.
 * Structure StructPhydro includes information about concerning for gas
 * particles.  Structure StructPStar includes information about concerning
 * for star particles.  StructPbody and StructPgas and StructPstar are
 * linked by bi-direction links.
 */
#pragma once 

#include <CELib.h>

typedef struct StructPbody_tag  StructPbody,*StructPbodyptr;
typedef struct StructPhydro_tag StructPhydro,*StructPhydroptr;
typedef struct StructPstar_tag  StructPstar,*StructPstarptr;
typedef struct StructPsink_tag  StructPsink,*StructPsinkptr;

extern struct StructPall{
    /* General Parameters */
    /* For Unit */
    double UnitMass;
    double UnitLength;
    double UnitTime;
    double DegreeOfFreedom;
    double HeliumWeight;
    double MeanMolecularWeight;
    double ConvertCoolingRate;
    double ConvertTtoU;
    double ConvertUtoT;
    double ConvertDensityToCGS;       // Get density in CGS unit.
    double ConvertNumberDensityToCGS; // To multiple this factor to the density in the simulation unit,
                                       // you can get the Number density in CGS unit.
    double ConvertErgsPreSec;
    double hubble;      // h x 100 km/s/Mpc
    /* For Unit */ 
    double GravConst;   // Gravitational Constant.
    /* General Parameters */

    /* Cosmological Parameters */ // This simulation code assumes flat universes.
    double OmegaB;      // OmegaB = rho_b/rho_crit. 
    double OmegaCDM;    // OmegaCDM = rho_cdm/rho_crit.
    double OmegaM;      // OmegaM = OmegaB + OmegaCDM.
    double OmegaL;      // OmegaL = 1.0 - OmegaM. 
    double Hubble;      // The Hubble paramter at z=0. hubble x 100 kms/s/Mpc(in simulation unit)
    double TCMB;        // The CMB temprature at now.
    double InitialRedshift;  // The initial redshift of the simulation.
    double FrozenRedshift;   // Softening length is comovingly evolved until z > FrozenRedsfhit. 
    double Lbox[3];        // For Periodic condition : a size of the simulation box.
                           // The volume should start from 0 for every direction.
    double Lboxh[3];       // For Periodic condition : Lbox/2
    double BoxCenter[3];   // For Periodic condition : Center of the simluation box.

    // Parameters depending redshift
    double HubbleZ;      // The Hubble paramter at the given redshift.
    double RHOcZ;        // Critical density at the given redshift. 
    double OverDensityZ; // OverDensity at the given redshift. // In SCDM, this value is fixed at 178.
    /* Cosmological Parameters */
    
    /* Simulation Parameters */
    double AdaptiveSofteningFactor;   // Adaptive Softening Factor. If you not use this, you must set 1.0.
    /* Simulation Parameters */

    /* Global Time Step */
    double    TStart;       // Simulation start time.
    double    TEnd;         // Simulation end time.
    double    TCurrent;     // Current time. 
    double    Redshift;     // Current Redshift.
    double    Era;       // dtmin*2^kmax + TCurrent. //if t>Era, BuildHierarchicalTimeStep()
    double    EraLocal;  // Tstep*dtmin.
    double    EraLocalEnd;    // dtmin*2^kmax. 
    double    EraStart;  // EraStart = TCurrent at the first step of HTS.
    double    dtnow;        // Current dt.
    double    dtmin;        // minimum time step within current block.
    double    dtmax;        // maximum time step within current block.
    unsigned int    TStep;
    unsigned int    TStepTotal;
    int       kmax;
#ifndef USE_VARIABLE_TIMESTEP
    double    dt_const; // dt for constant run.
#endif
    //unsigned long int TActiveParticlesAll; // total updated number of the simultion. 
    unsigned long int TActiveParticles;     // The total updated number of particles. 
    unsigned long int TActiveParticlesDM;   // The total updated number of dark matter particles. 
    unsigned long int TActiveParticlesHydro;// The total updated number of hydro particles. 
    unsigned long int TActiveParticlesStar; // The total updated number of star particles. 
    unsigned long int TActiveParticlesSink; // The total updated number of sink particles. 
    /* Global Time Step */

    /* Number of particles in current node */
    unsigned long int    NActivesAll;   // The number of active particles.
    unsigned long int    NActives;      // The number of active particles.
    unsigned long int    NActivesDM;    // The number of active "DM" particles.
    unsigned long int    NActivesHydro; // The number of active "Hydro" particles.
    unsigned long int    NActivesStars; // The number of active "Star" particles.
    unsigned long int    NActivesSink;  // The number of active "Sink" particles.

    unsigned long int    NActivesAll_t;   // The total number of active particles.
    unsigned long int    NActives_t;      // The total number of active particles.
    unsigned long int    NActivesDM_t;    // The total number of active "DM" particles.
    unsigned long int    NActivesHydro_t; // The total number of active "Hydro" particles.
    unsigned long int    NActivesStars_t; // The total number of active "Star" particles.
    unsigned long int    NActivesSink_t;  // The total number of active "Star" particles.

    unsigned long int    Ntotal;
    ///// unsigned long int    NDMBoundary;
    unsigned long int    NDM;
    unsigned long int    Nhydro;
    unsigned long int    Nstars;
    unsigned long int    Nsink;
    /* Number of particles in current node */

    /* Number of particles in all nodes */
    unsigned long int    Ntotal_t;
    ///// unsigned long int    NDMBoundary_t;
    unsigned long int    NDM_t;
    unsigned long int    Nhydro_t;
    unsigned long int    Nstars_t;
    unsigned long int    Nsink_t;
    /* Number of particles in all nodes */

    /* For SPH Particles */
    unsigned int Ns;
    unsigned int Npm; 
    double Gamma;   // Gamma : Polytropic index.
    double Gm1;     // Gamma - 1
    double GGm1;    // Gamma (Gamma-1)    
    double CS;      // The sound speed for the isothermal run.

    double HydroAlpha;   // Artificial Viscosty parameter, alpha 
    double HydroBeta;    // Artificial Viscosty parameter, beta 
    double HydroEta2;    // Artificial Viscosty parameter, eta^2
    double CoolingEnergyLoss; // The total amount of the energy loss by a radiative cooling.
    double ViscousAlphaMin;  // The minimum value of the alpha.
    double ViscousAlphaMax;  // The maximum value of the alpha.
    double ViscousS;         // The source factor for variable alpha model.
    double ViscousL;         // The damping scale for variable alpha model.
    /* For SPH Particles */

    /* For Sink Paricles */
    double SinkThresholdDensity;
    /* For Sink Paricles */

    /* For Run Status */
    int RunStatus;
    int NumProcs;
    int OutPutFileNumber;
    double OutPutInterval;
    double DampInterval;
    char BaseFileName[MaxCharactersInLine];
    char ASCIIFileName[MaxCharactersInLine];
    char RestartFileName[MaxCharactersInLine];
    /* For Run Status */

} Pall;

extern struct StructTimingResults{
    /* For Count CPU Time */
    double Total;

    double Decomposition;

    double Gravity;
    double GravityTree;
    double GravityComm;

    double Hydro;
    double HydroNeighborSearch;
    double HydroTree;
    double HydroComm;

    double HydroKernel;
    double HydroKernelNeighborSearch;
    double HydroKernelComm;

    double HydroDensity;
    double HydroDensityNeighborSearch;
    double HydroDensityComm;

    double HydroAcc;
    double HydroAccNeighborSearch;
    double HydroAccComm;

    double Cooling;

    double Feedback;
    double FeedbackNeighborSearch;
    double FeedbackComm;

    double Sink;
    double SinkNeighborSearch;
    double SinkTree;
    double SinkFormation;
    double SinkAccretion;
    double SinkComm;

    double Integral;

    double SortStructures;
    double KeyGeneration;
    double TimeStep;

    double GravImbarance;
    double HydroImbarance;
    double Starformation; 

    // This Step
    double TotalThisStep;

    double DecompositionThisStep;

    double GravityThisStep;
    double GravityTreeThisStep;
    double GravityCommThisStep;


    double HydroThisStep;
    double HydroNeighborSearchThisStep;
    double HydroTreeThisStep;
    double HydroCommThisStep;

    double HydroKernelThisStep;
    double HydroKernelNeighborSearchThisStep;
    double HydroKernelCommThisStep;

    double HydroDensityThisStep;
    double HydroDensityNeighborSearchThisStep;
    double HydroDensityCommThisStep;

    double HydroAccThisStep;
    double HydroAccNeighborSearchThisStep;
    double HydroAccCommThisStep;

    double CoolingThisStep;

    double FeedbackThisStep;
    double FeedbackNeighborSearchThisStep;
    double FeedbackCommThisStep;

    double SinkThisStep;
    double SinkNeighborSearchThisStep;
    double SinkTreeThisStep;
    double SinkFormationThisStep;
    double SinkAccretionThisStep;
    double SinkCommThisStep;

    double IntegralThisStep;

    double SortStructuresThisStep;
    double KeyGenerationThisStep;
    double TimeStepThisStep;

    double GravImbaranceThisStep;
    double HydroImbaranceThisStep;
    double StarformationThisStep;
    /* For Count CPU Time */
} TimingResults;


struct StructPbody_tag{
    StructPbodyptr Next; // <TMP> <LEAN>
    void    (*Baryon);   // <TMP> <LEAN>

    unsigned long int    GlobalID; // Global, unique ID for each particles.
    unsigned long long int OrderingKey; // <TMP> <LEAN> Key for ordering.
    unsigned long int    NextLeaf;   // <TMP> <LEAN> Next Leaf for Tree (Grav). This indicates NextLeaf in Pbody[NextLeaf].

    bool    Use;    // <TMP> <LEAN> If this structure is not use, the flag = false. 
    bool    Active; // <LEAN> This is used for the individual time step.

    int InteractionList; // <TMP> <LEAN>

    /* Physical quantities */
    double    Pos[3];       // <MIX> The position of the particle.
    double    PosP[3];      // <TMP> <LEAN> The predictor of the particle position.
    double    Vel[3];       // <MIX> The velocity of the particle.
    double    Velh[3];      // <MIX> <LEAN> The "half-step slided" velocity of the particle. For leapfrog integration.
    double    Acc[3];       // The acceleration of the particle. 
    double    AccOld[3];    // <LEAN> The acceleration of the particle. 
    double    Pot;          // The potential energy.
    double    Mass;         // <MIX> The mass of the particle.
    double    Eps;          // <MIX> The softening length.
    /* Physical quantities */

    /* Type of Particle */
    short    Type;
    /* Type of Particle */
    /* Time Step */
    short     k;            // <TMP> <LEAN>
    double    dt;           // <LEAN>
    double    EraLocal;  // <TMP> <LEAN>
    /* Time Step */
};


#define PbodyHydro(_index)              ((StructPhydroptr)(Pbody[_index]->Baryon))
#define PbodyHydroVelP(_index)          (((StructPhydroptr)(Pbody[_index]->Baryon))->VelP)
#define PbodyHydroKernel(_index)        (((StructPhydroptr)(Pbody[_index]->Baryon))->Kernel)
#define PbodyHydroKernelPred(_index)    (((StructPhydroptr)(Pbody[_index]->Baryon))->KernelPred)
#define PbodyHydroRho(_index)           (((StructPhydroptr)(Pbody[_index]->Baryon))->Rho)
#define PbodyHydroRhoPred(_index)       (((StructPhydroptr)(Pbody[_index]->Baryon))->RhoPred)
#define PbodyHydroAcc(_index)           (((StructPhydroptr)(Pbody[_index]->Baryon))->HydroAcc)
#define PbodyHydroU(_index)             (((StructPhydroptr)(Pbody[_index]->Baryon))->U)
#define PbodyHydroUPred(_index)         (((StructPhydroptr)(Pbody[_index]->Baryon))->UPred)
#define PbodyHydroDu(_index)            (((StructPhydroptr)(Pbody[_index]->Baryon))->Du)
#define PbodyHydroDuCooling(_index)     (((StructPhydroptr)(Pbody[_index]->Baryon))->DuCooling)
#define PbodyHydroDuPrev(_index)        (((StructPhydroptr)(Pbody[_index]->Baryon))->DuPrev)
#define PbodyHydroDiv(_index)           (((StructPhydroptr)(Pbody[_index]->Baryon))->Div)
#define PbodyHydroRot(_index)           (((StructPhydroptr)(Pbody[_index]->Baryon))->Rot)
#define PbodyHydroVsig(_index)          (((StructPhydroptr)(Pbody[_index]->Baryon))->Vsig)
#define PbodyHydroZ(_index)             (((StructPhydroptr)(Pbody[_index]->Baryon))->Z)
#define PbodyHydroZII(_index)           (((StructPhydroptr)(Pbody[_index]->Baryon))->ZII)
#define PbodyHydroZIa(_index)           (((StructPhydroptr)(Pbody[_index]->Baryon))->ZIa)

#define PbodyStar(_index)       ((StructPstarptr)(Pbody[_index]->Baryon))
#define PbodyStarZ(_index)      (((StructPstarptr)(Pbody[_index]->Baryon))->Z)
#define PbodyStarZII(_index)    (((StructPstarptr)(Pbody[_index]->Baryon))->ZII)
#define PbodyStarZIa(_index)    (((StructPstarptr)(Pbody[_index]->Baryon))->ZIa)
#define PbodyStarTypeII(_index) (((StructPstarptr)(Pbody[_index]->Baryon))->TypeII)
#define PbodyStarTypeIa(_index) (((StructPstarptr)(Pbody[_index]->Baryon))->TypeIa)

#define PbodySink(_index)       ((StructPsinkptr)(Pbody[_index]->Baryon))
#define PbodySinkAM(_index)     (((StructPsinkptr)(Pbody[_index]->Baryon))->AM)
#define PbodySinkFormationTime(_index)      (((StructPsinkptr)(Pbody[_index]->Baryon))->FormationTime)
#define PbodySinkZ(_index)      (((StructPsinkptr)(Pbody[_index]->Baryon))->Z)
#define PbodySinkZII(_index)    (((StructPsinkptr)(Pbody[_index]->Baryon))->ZII)
#define PbodySinkZIa(_index)    (((StructPsinkptr)(Pbody[_index]->Baryon))->ZIa)


struct StructPhydro_tag{
    StructPhydroptr Next;   // <TMP> <LEAN>
    StructPbodyptr  Body;   // <TMP> <LEAN>
    //StructNodeptr   HostNode;   //The pointer of the lowest node that this particle hosts.

    bool Use;                // <TMP> <LEAN>
    unsigned int ExportFlag; // <TMP> <LEAN> if you use NProcs>32, increase this bit field.
    unsigned int Nlist;      // <LEAN>
    unsigned long int NextLeaf; // <TMP> <LEAN> Next Leaf for Tree (Neighbor Search).
                                // This flag indicates NextLeaf in Phydro[NextLeaf].
    bool    Active;          // <TMP> <LEAN> if "OFF", the particle skips the cooling routine.
    bool    CoolingFlag;     // <TMP> <LEAN> if "OFF", the particle skips the cooling routine.
    //bool    HIIFlag;         // <TMP> <LEAN>

    /* ID in Leaves. */
    int Leaf; // <LEAN> The ID to the HydroTree leaves.
    /* ID in Leaves. */

    /* Time step */
    bool      GravityKickFlag; // <LEAN> The kick flag of gravitataional acceleration.
    int       k_hydro;         // <LEAN>
    double    dt_hydro;        // <LEAN>
    double    EraLocal_hydro;  // <LEAN> The local time in a current block.
//#ifdef HYDRO_TIMESTEP_LIMITER
    bool      HydroTimeStepLimiterFlag; // <TMP> <LEAN>
    int       k_hydro_localmin;  // <LEAN>
    int       k_hydro_localmin_old;  // <LEAN> 
    double    NextUpdateEra;     // <LEAN> The next update timing in Era. This variable is used when the HydroTimeStepLimiterFlag is true.
    double    dt_hydro_localmin; // <LEAN> 
//#endif
    /* Time step */

    double    Mass;         // <TMP> <LEAN>
    double    PosP[3];      // <TMP> <LEAN>
    double    VelP[3];      // <TMP> <LEAN>

    double    Rho;
    double    RhoPred;    // <TMP> <LEAN> Predictor of density.
    double    EnergyDensity; // Energy density used in DISPH
    double    EnergyDensityPred; // <TMP> <LEAN> Predictor of the energy density used in DISPH.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER
    double    SmoothedNumber;    // <TMP> <LEAN> Smoothed Number.
#endif
    double    Kernel;     // Size of kernel. 0<- rij/Kernel ->2
    double    KernelPred; // <TMP> <LEAN> Predictor of kernel. 0<- rij/Kernel ->2
// #ifdef USE_GRAD_H //{
    double    Gradh; // <LEAN> Grad-h term
// #ifdef USE_GRAD_N //{
    double    NumberDensity; // <LEAN> Number density
    double    GradN;         // <LEAN> Grad-N term
    double    fij;           // <LEAN> Part of fij. fij = 1-fij/(Uj)
// #endif // USE_GRAD_N //}
// #endif // USE_GRAD_H //}
    double    Div;        // <LEAN> 
    double    Rot[3];     // <LEAN> 
    double    F;          // <LEAN> 
    double    HydroAcc[3];
    double    U;          // The internal energy of gas.
    double    UPred;      // <TMP> <LEAN>
    double    Du;         // Du/Dt.
    double    DuPrev;     // <LEAN> Du/Dt.
    double    DuCooling;  // <LEAN> The amount of the energy loss by the radiative cooling.

    // Variables for SPSPH
    double    PseudoDensity;          // The fundamental quantity of SPSPH
    double    PseudoDensityPred;      // <TMP> The predictor of y
    double    Zw;         // The weight of SPHSPH
    double    ZwPred;     // <TMP> The predictor of the weight used in SPSPH
    double    DZw;        // The time derivative of the weight used in SPSPH
    double    Ddif;       // The diffusion parameter for SPSPH


    double    DQheat;     // <LEAN> The amount of the heating energy by SNe, which injects during dt.
    double    dMass;      // <LEAN> The additional mass by SNe.
    double    Z;          // Metallicity of gas.
    double    ZII;        // The mass of metal by SNII.
    double    ZIa;        // The mass of metal by SNIa.
    double    dZII;       // The mass of metal by SNII.
    double    dZIa;       // The mass of metal by SNIa.
    double    Vsig;       // <LEAN> Signal Velocity  or Maximum Velocity difference.
    double    Alpha;      //
    double    DAlpha;     // <LEAN>

    // Metal diffusion
    double    DZdiff;
    double    Sxy;
    double    Sxz;
    double    Syx;
    double    Syz;
    double    Szx;
    double    Szy;

#if VISCOSITY_TYPE == 1 //{
    double    Bxx;
    double    Bxy;
    double    Bxz;
    double    Byx;
    double    Byy;
    double    Byz;
    double    Bzx;
    double    Bzy;
    double    Bzz;
#endif // VISCOSITY_TYPE == 1 //}

    //double    HIILastContactTime;

    // H2 formation
#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    double  G0;
    double  fH2;
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#if (UseSFModelSpawn)
    short   SpawnTimes; // <LEAN>
    double  SpawnMass;  // <LEAN>
#endif
#ifdef USE_CELIB
    double    Elements[CELibYield_Number]; // <MIX> H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe, Ni, Eu
#endif // USE_CELIB
//#ifdef USE_FUVFEEDBACK //{
    double GradRho;               // <LEAN>
    double G0thin;                // <LEAN>
    double G0thinLocal;           // <LEAN>
    double G0extLocal;            // <LEAN>
    double G0thinNextUpdateTime;  // <TMP> <LEAN>
//#endif // USE_FUVFEEDBACK //}

#ifdef USE_RADIATION_PRESSURE //{
    double MomentumRP[3];               // momentum injection by radiation
#endif //USE_RADIATION_PRESSURE //{

#ifdef USE_MOMENTUM_FEEDBACK //{
    double MomentumFB[3];               // momentum injection by SNe
#endif //USE_MOMENTUM_FEEDBACK //{

// #ifdef USE_MULTIPHASE_MODEL //{
    double Uhot;            // Hot gas temperature
    double Rhohot;          // Hot gas density
    double Mhot;            // Hot gas mass
    double MultiphaseEndTime; // 10xt_cross
    bool   MultiphaseFlag;  // If true, the gas particle is evolved under the multiphase mode.
// #endif //USE_MULTIPHASE_MODEL //}

#ifdef USE_PARTICLE_SPLITTING //{
    bool   SplitFlag;
    int    SplitGeneration; // Set particle generation.
    int    SplitNthChildren;        // Set particle generation.
    unsigned long int    AncestorGlobalID; // Ancestor's Global ID.
#endif // USE_PARTICLE_SPLITTING //}
#ifdef USE_SMOOTHED_POTENTIAL //{
    double Pot;
#endif // USE_SMOOTHED_POTENTIAL //}

#ifdef USE_PARTICLE_TAG
    int Tag;
#endif
};




#define PhydroBody(_index)   (Phydro[_index]->Body)
#define PhydroPos(_index)    (Phydro[_index]->Body->Pos)
#define PhydroPosP(_index)   (Phydro[_index]->PosP)
#define PhydroVel(_index)    (Phydro[_index]->Body->Vel)
#define PhydroAcc(_index)    (Phydro[_index]->Body->Acc)
#define PhydroMass(_index)   (Phydro[_index]->Mass)
#define PhydroActive(_index) (Phydro[_index]->Active)

struct StructPstar_tag{
    StructPstarptr  Next;   // <TMP> <LEAN>
    StructPbodyptr  Body;   // <TMP> <LEAN>

    bool    Use;            // <TMP> <LEAN>
    bool    TypeII;         // fb flag
    bool    TypeIa;         // fb flag
    bool    HIIflag;        // <TMP> <LEAN>
#if defined(PRESERVE_SNII_EVENTRATE)
    bool    TypeIIProb;     // fb flag
#endif

    short   IMFTYPE;
    short   NthChildren;
    unsigned long int   ParentGlobalID; // Global ID of parent particles.

    double  InitialMass;    // <MIX> <LEAN> Initial Mass 
    double  Mass;           // <MIX> Current Mass
    double  FormationTime;  // 
    double  Z;              // Metallicity of gas.
    double  ZII;            // The mass of metal by TypeII.
    double  ZIa;            // The mass of metal by TypeIa.
    double  TempForm;       // <LEAN> Temperature(t=FormationTime)
    double  RhoForm;        // <LEAN> Density(t=FormationTime) // UnitMass/UnitLength^3
    double  StromgrenRadius;// <LEAN> 
    double  Density;        // <TMP> <LEAN> 
#ifdef USE_CELIB //{
    int       SNIICount;  // Use this counter for SNeII
    int       SNIaCount;  // Use this counter for SNeIa
    double    EventTimeSNII;    // Next event time of SNeII
    double    EventTimeSNIa;    // Next event time of SNeIa
// #ifdef USE_CELIB_AGB //{
    int       AGBCount;         // Number count for AGB mass loss.
    double    EventTimeAGB; // Next event time of AGB mass loss.
// #endif // USE_CELIB_AGB //}
#ifdef USE_CELIB_NSM //{
    int       NSMCount;         // Number of NSM event.
    double    EventTimeNSM; // Next event timing of NSM.
#endif // USE_CELIB_NSM //}
    double    Elements[CELibYield_Number]; // <MIX> H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe, Ni, Eu
#endif // USE_CELIB //}

//#ifdef USE_FUVFEEDBACK //{
    double LFUV;
//#endif //USE_FUVFEEDBACK //}

#ifdef USE_RADIATION_PRESSURE //{
    double BolometricLuminosity;
    double RadiusRP;
#endif //USE_RADIATION_PRESSURE //{
#ifdef USE_STELLAR_WIND //{
    double StellarWindEnergy;
    double RadiusSW;
#endif //USE_STELLAR_WIND //{
            
#ifdef USE_FEEDBACK_TIMESTEP //{
    double dt_fb; // <TMP> <LEAN>
#endif // USE_FEEDBACK_TIMESTEP //}

#ifdef USE_MOMENTUM_FEEDBACK //{
    double nH_ave; // Neighbor averaged number density
    double Z_ave;  // Neighbor averaged metallicity
#endif // USE_MOMENTUM_FEEDBACK //}

#ifdef USE_PARTICLE_SPLITTING //{
    int    SplitGeneration; // Set particle generation.
    int    SplitNthChildren;        // Set particle generation.
    unsigned long int   GrandParentGlobalID; // Grand Parent Global ID.
#endif // USE_PARTICLE_SPLITTING //}
};



#define PstarBody(_index) (Pstar[_index]->Body)
#define PstarPos(_index)  (Pstar[_index]->Body->Pos)
#define PstarPosP(_index) (Pstar[_index]->Body->PosP)
#define PstarVel(_index)  (Pstar[_index]->Body->Vel)
#define PstarAcc(_index)  (Pstar[_index]->Body->Acc)
#define PstarMass(_index) (Pstar[_index]->Body->Mass)
#define PstarActive(_index) (Pstar[_index]->Body->Active)


struct StructPsink_tag{
    StructPsinkptr  Next;   // <TMP> <LEAN> 
    StructPbodyptr  Body;   // <TMP> <LEAN> 

    bool    Use;            // <TMP> <LEAN> 
    unsigned long int   ParentGlobalID; // Global ID of parent particles.

    double  PosP[3];        // <TMP> <LEAN> 
    double  VelP[3];        // <TMP> <LEAN> 
    double  dt_localmin;    // <LEAN> The local minimum time-steps.

    int     NumberofAbsorbedParticles;    // The number of absorbed gas particles.
    double  AccretionRadius;// The sink-hydro accretion radius.
    double  MergingDistance;// The sink-sink merging distance.
    double  AM[3];          // The values of angular momentum.
    double  FormationTime;  // The formation epoch.
    double  Z;              // The metallicity of the sink particle.
    double  ZII;            // The mass of metal due to TypeII SNe.
    double  ZIa;            // The mass of metal due to TypeIa SNe.
    double  AccretionMass;     // The accreted mass onto this sink particle.
    double  AccretionMassGas;  // The accretion/absorbtion gas mass.
    double  AccretionMassStar; // The accretion/absorbtion star mass.
    double  AccretionMassToBH; // The accretion mass onto the BH.

#ifdef USE_CELIB
    double  Elements[CELibYield_Number]; // <MIX> H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe, Ni
#endif // USE_CELIB
#ifdef USE_PARTICLE_SPLITTING //{
    int    SplitGeneration; // Set particle generation.
    int    SplitNthChildren;        // Set particle generation.
    unsigned long int   GrandParentGlobalID; // Grand Parent Global ID.
#endif // USE_PARTICLE_SPLITTING //}
};

#define PsinkBody(_index) (Psink[_index]->Body)
#define PsinkPos(_index)  (Psink[_index]->Body->Pos)
#define PsinkPosP(_index) (Psink[_index]->Body->PosP)
#define PsinkVel(_index)  (Psink[_index]->Body->Vel)
#define PsinkAcc(_index)  (Psink[_index]->Body->Acc)
#define PsinkMass(_index) (Psink[_index]->Body->Mass)
#define PsinkActive(_index) (Psink[_index]->Body->Active)

/** Element Pointers **/
extern StructPbodyptr   PbodyElements;
extern StructPhydroptr  PhydroElements;
extern StructPstarptr   PstarElements;
extern StructPsinkptr   PsinkElements;
extern unsigned int PbodyElementsSize;
extern unsigned int PhydroElementsSize;
extern unsigned int PstarElementsSize;
extern unsigned int PsinkElementsSize;
/** Element Pointers **/

/** Direct Access Pointers **/
extern StructPbodyptr   *Pbody;
extern StructPhydroptr  *Phydro;
extern StructPstarptr   *Pstar;
extern StructPsinkptr   *Psink;
extern unsigned int PbodySize;
extern unsigned int PhydroSize;
extern unsigned int PstarSize;
extern unsigned int PsinkSize;
/** Direct Access Pointers **/


/* Communication Buffer */
extern int *NumberofBufferExportSendAllocated;
extern int NumberofBufferExportRecvAllocated;
extern int NumberofBufferImportSendAllocated;
extern int *NumberofBufferImportRecvAllocated;
extern void **BufferExportSend;
extern void *BufferExportRecv;
extern void *BufferImportSend;
extern void **BufferImportRecv;
/* Communication Buffer */

/* Hydro Interaction Buffer */
extern int NumberofBufferHydroInteractionFlagsAllocated;
extern bool *BufferHydroInteractionFlags;
/* Hydro INteraction  Buffer */

/* Structures for Domain Decomposition */
extern struct StructInfoBiSection{
    int Left;
    int Right;
    double Pos; // center of the node 
    int Axis; // Axis for bisection 
} *InfoBiSection;
/* Structures for Domain Decomposition */

/* Structures for Communication */
extern struct StructCommunicationTable{
    int Level;
    int SendRank;
    int SendSize;
    int RecvRank;
    int RecvSize;
#ifdef USE_BARYON_COMM //{
    int HydroSendRank;
    int HydroRecvRank;
    int ActiveHydroSendRank;
    int ActiveHydroRecvRank;
    int BaryonSendRank;
    int BaryonRecvRank;
    int ActiveBaryonSendRank;
    int ActiveBaryonRecvRank;
#endif // USE_BARYON_COMM //}
}   *CommunicationTable;
/* Structures for Communication */

//container for particle data exchange

/* Structures for Domain Decomposition */

/* Definition of structures for tree */
#define	TreeNsub		(8)  // Three dimension
#define	TreeDiagonal	(0.8661e0) // Diagonal of a cube with the width = 2

typedef struct	StructGravityRoot_tag  StructGravityRoot,*StructGravityRootptr;
struct StructGravityRoot_tag {
	int	NumberofLeaves;            // Current number of used nodes.
	int	NumberofAllocatedLeaves;   // Maximum allocated memory for *Leaves.
	int	*Leaves;                   // Leaves, which is morton-ordered local particles index list.
	int	NumberofNodes;             // Whole allocated nodes for this tree.
	int	NumberofAllocatedNodes;    // Whole allocated nodes for this tree.
    int NumberofNodeCreationLimit; // Maximum number of particles for the lowest node.
	int	NumberofLeavesInGroup;     // n_g for tree + grape.

    short BaseNumberofChildren;    // The maximum number of children is derived from the cube of this number.

    short CurrentMaxLevel;         // The maximum level of this tree. 
    short MaxLevel;                // Allowing maximum level of this tree. 

    double PosMax[3];              // Particle distribution edges, maximum.
    double PosMin[3];              // Particle distribution edges, minimum.
    double Width;                  // We here assume the root is a cube.
    double WidthFactor[TreeMaxNodeLevel]; // Width * WidthFactor[LevelofNode] = Width of the Node.
	double OpeningAngle;
};

extern struct StructGravityNode{
	int Next;
	int Parent;
	int Children;
	int Sister;
	//int Traverse;

    unsigned long long int OrderingKey; // Key for ordering.

	short	Level;            // Level of node in belonging tree structure. // this can be short.
	int Leaves;               // This member points a Pbody or Phydro in the top of the node.
	short NumberofChildren;   // Number of child-node in this node.
	int NumberofLeaves;       // Number of particles in this node.
	int	NumberofActiveLeaves; // Number of `active' particles in this node. 

	double	Pos[3];           // Center of this node.
	double	COM[3];           // Center of mass of this node.
	double	Mass;             // Mass of this node. 
    double  DistanceMax;      // The maximum distance between particles and the center of this node.
#ifdef USE_SYMMETRIZED_SOFTENING
    double  Eps2;             // Mass weighted (softening length)^2.
    double  EpsMax;           // The maximum gravitational softening length.
    double  EpsMin;           // The minimum gravitational softening length.
#endif
} *GravityNode;

extern struct StructGravityCache{
    double Pos[3];
    double Mass;
    double Eps;
    int  Leaf;
    bool Active;
} *GravityCache;
extern struct StructGravityAccPotCache{
    double Acc[3];
    double Pot;
    int  InteractionList;
} *GravityAccPotCache;
extern StructGravityRoot GravityRoot;

typedef struct	StructHydroRoot_tag  StructHydroRoot,*StructHydroRootptr;

struct StructHydroRoot_tag {
    int LifeGauge;                    // Life gauge of this tree. Check also ``HYDROTREE_LIFEGAUGE''.
	int	NumberofLeaves;               // Current number of used nodes.
	int	NumberofAllocatedLeaves;      // Maximum allocated memory for *Leaves.
	int	*Leaves;                      // Leaves, which is morton-ordered local particles index list.
	int	NumberofNodes;                // Whole allocated nodes for this tree.
	int	NumberofAllocatedNodes;       // Whole allocated nodes for this tree.
    int NumberofNodeCreationLimit;    // Maximum number of particles for the lowest node.

    short BaseNumberofChildren;       // The maximum number of children is derived from the cube of this number.

    short CurrentMaxLevel;            // The maximum level of this tree. 
    short MaxLevel;                   // Allowing maximum level of this tree. 

    double PosMax[3];                 // Maximum courner of particle extention. 
                                      // Not the root box.
    double PosMin[3];                 // Minimum courner of particle extention. 
                                      // Not the root box.
    double Width;                     // We here assume that the root is a cube.
    double WidthFactor[TreeMaxNodeLevel]; // Width * WidthFactor[LevelofNode] = Width of the Node.
};

extern struct StructHydroNode{
	int Next;
	int Parent;
	int Children;
	int Sister;
	//int Traverse; // This index is necessary for the update of the node information. 

    unsigned long long int OrderingKey; // Key for ordering.

	short Level;              // Level of node in belonging tree structure. // this can be short.
	short NumberofChildren;   // Number of child-node in this node.
	int Leaves;               // This member points a Pbody or Phydro in the top of the node.
	int NumberofLeaves;       // Number of particles in this node.
	int	NumberofActiveLeaves; // Number of `active' particles in this node. 

	double	Pos[3];           // Center of this node.
	double	KernelMax;        // Maximum length of kernels.
    double  DistanceMax;      // Maximum distance between particles and the center of this node.
} *HydroNode,*HydroNodeImported,*StellarNode,*SinkNode;

typedef struct StructNBCache_tag  StructNBCache;
struct StructNBCache_tag{
    double Pos[3];            // Position of this leaf.
    double Kernel;            // Kernel size of this leaf.
    int  Leaf;                // ID to refer the `Phydro' array.
    bool Active;              // Active leaf check flag.
};
extern StructHydroRoot HydroRoot;
extern StructNBCache *NBCache;

extern StructHydroRoot HydroRootImported;
extern StructNBCache *NBImportedCache;

extern StructHydroRoot StellarRoot;
extern StructNBCache *StellarNBCache;

extern StructHydroRoot SinkRoot;
extern StructNBCache *SinkNBCache;


typedef struct	StructYoungStarRoot_tag  StructYoungStarRoot,*StructYoungStarRootptr;
struct StructYoungStarRoot_tag {
	int	NumberofLeaves;            // Current number of used nodes.
	int	NumberofAllocatedLeaves;   // Maximum allocated memory for *Leaves.
	int	*Leaves;                   // Leaves, which is morton-ordered local particles index list.
	int	NumberofNodes;             // Whole allocated nodes for this tree.
	int	NumberofAllocatedNodes;    // Whole allocated nodes for this tree.
    int NumberofNodeCreationLimit; // Maximum number of particles for the lowest node.
	int	NumberofLeavesInGroup;     // n_g for tree + grape.

    short BaseNumberofChildren;    // The maximum number of children is derived from the cube of this number.

    short CurrentMaxLevel;         // The maximum level of this tree. 
    short MaxLevel;                // Allowing maximum level of this tree. 

    double PosMax[3];              // Particle distribution edges, maximum.
    double PosMin[3];              // Particle distribution edges, minimum.
    double Width;                  // We here assume the root is a cube.
    double WidthFactor[TreeMaxNodeLevel]; // Width * WidthFactor[LevelofNode] = Width of the Node.
	double OpeningAngle;
};

extern struct StructYoungStarNode{
	int Next;
	int Parent;
	int Children;
	int Sister;

	short	Level;            // Level of node in belonging tree structure. // this can be short.
	int Leaves;               // This member points a Pbody or Phydro in the top of the node.
	short NumberofChildren;   // Number of child-node in this node.
	int NumberofLeaves;       // Number of particles in this node.
	int	NumberofActiveLeaves; // Number of `active' particles in this node. 

	double	Pos[3];           // Center of this node.
	double	COM[3];           // Center of mass of this node.
	double	LFUV;             // Mass of this node. 
    double  DistanceMax;      // The maximum distance between particles and the center of this node.
} *YoungStarNode;

extern struct StructYoungStarCache{
    double Pos[3];
    double LFUV;
    double Eps;
    int  Leaf;
    bool Active;
} *YoungStarCache;

extern struct StructYoungStarResultCache{
    double Acc[3];
    double Pot;
    int  InteractionList;
} *YoungStarResultCache;
extern StructYoungStarRoot YoungStarRoot;
/* Definition of structures for tree */

/* 
 * This structure holds the edges of the particle disributions.
 */
extern struct StructEdges{ 
    double PosMax[3];
    double PosMin[3];
} *EdgesForHydro,
  *EdgesForActiveHydro,
  *EdgesForGravity,
  *EdgesForStars,
  *EdgesForSink,
  *EdgesForYoungStars,
  *EdgesForActiveSink;

/* Friend-of-Friend */
extern int NFOFGroups;
extern int FOFSize;
extern int FOFCatalogSize;

extern struct StructFOF{
    int Next;
    int Head;
    int Tail;
} *FOF;

extern struct StructFOFCatalog{
    //int ID;
    int Number;
    int Head;
    double Mass;
    double Pos[3];
    double Vel[3];
} *FOFCatalog;
/* Friend-of-Friend */

/* Log file pointers */
extern FILE *FpEnergy;
extern FILE *FpMomentum;
extern FILE *FpAMomentum;
extern FILE *FpElapsedTime;
extern FILE *FpStepLog;
extern char FnameEnergy[MaxCharactersInLine];
extern char FnameMomentum[MaxCharactersInLine];
extern char FnameAMomentum[MaxCharactersInLine];
extern char FnameElapsedTime[MaxCharactersInLine];
extern char FnameStepLog[MaxCharactersInLine];
/* Log file pointers */

/* Random Generator */
extern gsl_rng *RandomGenerator;
/* Random Generator */

/* MPI Communicators */
#ifdef USE_BARYON_COMM //{
extern MPI_Group hydro_group;
extern MPI_Group activehydro_group;
extern MPI_Comm MPI_HYDRO_COMM_WORLD;
extern MPI_Comm MPI_ACTIVEHYDRO_COMM_WORLD;
extern MPI_Group baryon_group;
extern MPI_Comm MPI_BARYON_COMM_WORLD;
extern MPI_Group activebaryon_group;
extern MPI_Comm MPI_ACTIVEBARYON_COMM_WORLD;
#endif // USE_BARYON_COMM //}
/* MPI Communicators */

