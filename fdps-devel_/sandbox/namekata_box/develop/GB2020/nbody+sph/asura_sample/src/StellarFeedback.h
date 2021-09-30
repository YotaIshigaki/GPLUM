#pragma once 

#define EffSAVecSize (7)

#define __CHECK_SUM__
#define __CHECK_WEIGHT__


enum{
    StellarFeedbackType_SNII = CELibFeedbackType_SNII,
    StellarFeedbackType_SNIa = CELibFeedbackType_SNIa,
    StellarFeedbackType_AGB  = CELibFeedbackType_AGB,
    StellarFeedbackType_NSM  = CELibFeedbackType_NSM,
#ifdef USE_RADIATION_PRESSURE //{
    StellarFeedbackType_RP,
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
    StellarFeedbackType_SW,
#endif // USE_STELLAR_WIND //}
    StellarFeedbackType_Number,
};

#if 1

struct StructActiveStellarFeedbackParticle {
    int Index; // 
    int Nlist; // Number of neighbors.
    int k_hydro_localmin;
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double SmoothedNumber; // Smoothed mas..
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
    double Pos[3];
    double Radius;  // Feedback radius is 2*Radius.
    double Density; // The mass weighted nomalization factor.
#ifdef SET_SNII_TEMPERATURE
    double GasMass; // Local gas mass.
#endif //SET_SNII_TEMPERATURE
#ifdef MAXIMUM_ENERGY_INPUT
    double DistanceMin;
    unsigned long int DistanceMinGlobalID;
#endif
#ifdef USE_MOMENTUM_FEEDBACK //{
    double nH_ave;
    double Z_ave;
    double NumberDensity;
    double SA;
    double WeightCorrection[EffSAVecSize];
    double CheckWeight;
#endif // USE_MOMENTUM_FEEDBACK //}
    double WeightSum;
    bool LocalUpdateFlags; // Flag for the local update.
    int Type; // Follow enum StellarFeedbackType_*
    int Count;
    double InitialMass;
    double Metallicity;
    double Rvalue;
    double Lvalue;
    int IterationCount;
#ifdef __CHECK_SUM__ //{
    int CheckSum;
#endif // __CHECK_SUM__ //}
#ifdef USE_STAR_TIMESTEP_LIMITER //{
    int k;
#endif // USE_STAR_TIMESTEP_LIMITER //]
};
#else

struct StructActiveStellarFeedbackSNParticle {
    double Density; // The mass weighted nomalization factor.
#ifdef SET_SNII_TEMPERATURE //{
    double GasMass; // Local gas mass.
#endif //SET_SNII_TEMPERATURE //}
#ifdef MAXIMUM_ENERGY_INPUT //{
    double DistanceMin;
    unsigned long int DistanceMinGlobalID;
#endif // MAXIMUM_ENERGY_INPUT //}
    double InitialMass;
    double Metallicity;
    int Count;
    double WeightSum;
#ifdef __CHECK_SUM__ //{
    int CheckSum;
#endif // __CHECK_SUM__ //}

};

#ifdef USE_RADIATION_PRESSURE //{
struct StructActiveStellarFeedbackRPParticle {
    // double GasMass;
    // double MetalMass;
    double Metallicity;
    double WeightSum;
};
#endif // USE_RADIATION_PRESSURE //}

struct StructActiveStellarFeedbackParticle {
    int Index; // 
    int Nlist; // Number of neighbors.
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
    double SmoothedNumber; // Smoothed mas..
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
    double Pos[3];
    double Radius;  // Feedback radius is 2*Radius.
    union{
        struct StructActiveStellarFeedbackSNParticle SN;
#ifdef USE_RADIATION_PRESSURE //{
        struct StructActiveStellarFeedbackRPParticle RP;
#endif // USE_RADIATION_PRESSURE //}
    } Body;
    bool LocalUpdateFlags; // Flag for the local update.
    int Type; // Follow enum StellarFeedbackType_*
    double Rvalue;
    double Lvalue;
    int IterationCount;
};
#endif

void InitializeStellarFeedback(void);
int StellarFeedbackGetIMFType(void);
int StellarFeedbackGetIaType(void);
int CountFeedbackNumber(void);
double StellarFeedbackRadiusInitialGuess(const int Index);
double StellarFeedbackGetNextExplosionTime(const int mode, const double Metallicity, const double InitialMass, const int Count);
struct StructStellarFeedbackLocalInfo RetrunStellarFeedbackLocalInfo(double Pos[restrict], const double Radius);
void StellarFeedback(void);
void StellarFeedbackGatherScatterESA(void);

#ifdef TASK_TEST_STELLARFEEDBACK //{
void StellarFeedbackTestCalcFeedbackRadius(int NActives, const int IndexList[restrict], const int TypeList[restrict], double Radius[restrict], int Nlist[restrict], int CheckSum[restrict]);
void StellarFeedbackTestReleaseEnergyHeavyElements(const int NExplosion, const int IndexList[restrict], const int TypeList[restrict]);
#endif // TASK_TEST_STELLARFEEDBACK //}

