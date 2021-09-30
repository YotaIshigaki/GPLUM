#include "config.h"

static double RP_lifetime;
static double dt_RP_max;
//static double (*RPfunction)(const int); // Evaluate the amount of radiation pressure from a star particle

#define Zsun (0.0134)

/*
 * This function returns the released momentum of radation in cgs units.
 * Functional form is based on Cerevino et al. 2014
 */
static inline double __attribute__((always_inline)) RPfunctionCerevino(const int index){

#ifdef USE_RADIATION_PRESSURE //{
    return (RADIATION_PRESSURE_ETA1 + RADIATION_PRESSURE_ETA2*RADIATION_PRESSURE_TAUIR)*
        Pstar[index]->BolometricLuminosity/LIGHT_C_CGS;
#endif // USE_RADIATION_PRESSURE //}
}

/*
 * This function returns the released momentum of radation in cgs units.
 * Functional form is based on Okamoto et al. 2014
 */
static inline double __attribute__((always_inline)) RPfunctionOkamoto(const int index){

#ifdef  USE_RADIATION_PRESSURE //{
    return (RADIATION_PRESSURE_ETA1 + RADIATION_PRESSURE_TAU0*(Pstar[index]->Z/Zsun))*
        Pstar[index]->BolometricLuminosity/LIGHT_C_CGS;
#endif // USE_RADIATION_PRESSURE //}
}

double GetRadiationPressure(int index){

#if RADIATION_PRESSURE_MODEL == 0 //{
    return RPfunctionCerevino(index);
#elif RADIATION_PRESSURE_MODEL == 1 //}//{
    return RPfunctionOkamoto(index);
#else  // RADIATION_PRESSURE_MODEL //}//{
    return 0;
#endif // RADIATION_PRESSURE_MODEL //}

}

void UpdateBolometricLuminosity(const int index){

    const double MassConversionFactor = Pall.UnitMass/MSUN_CGS;
    const double TimeConversionFactor = Pall.UnitTime/YEAR_CGS;
#ifdef  USE_RADIATION_PRESSURE //{
    Pstar[index]->BolometricLuminosity = 
        MassConversionFactor*Pstar[index]->InitialMass*
        ASRFLXGetBol((Pall.TCurrent-Pstar[index]->FormationTime)*TimeConversionFactor,
                Pstar[index]->Z); 
#endif // USE_RADIATION_PRESSURE //}

    return ;
}


void InitRadiationPressure(void){

    RP_lifetime = RADIATION_PRESSURE_SOURCE_LIFETIME/Pall.UnitTime;
    dt_RP_max = RADIATION_PRESSURE_MAX_TIMESTEP/Pall.UnitTime;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Init RaidationPressure\n");
        fprintf(stderr," RP_lifetime = %g\n",RP_lifetime);
        fprintf(stderr," dt_RP_max = %g\n",dt_RP_max);
    }
    return;
}

double GetRadiationPressureSourceLifeTime(void){
    return RP_lifetime;
}

double GetRadiationPressureMaxTimestep(void){
    return dt_RP_max;
}

bool CheckRadiationPressureSource(const int index){
    if((Pstar[index]->TypeII == false)
            &&(Pall.TCurrent-Pstar[index]->FormationTime < GetRadiationPressureSourceLifeTime())){
        return true;
    } else {
        return false;
    }
}

bool CheckRadiationPressureSourceThroughPbody(const int index){
    if((PbodyStar(index)->TypeII == false)
            &&(Pall.TCurrent-PbodyStar(index)->FormationTime < GetRadiationPressureSourceLifeTime())){
        return true;
    } else {
        return false;
    }
}

int CountRadiationPressureSourceNumber(void){
    int counter = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarBody(i)->Active){
            //if(Pall.TCurrent-Pstar[i]->FormationTime < GetRadiationPressureSourceLifeTime()){
            if(CheckRadiationPressureSource(i)){
                counter ++;
            }
        }
    }
    return counter;
}

void RadiationPressureKick(void){

#ifdef  USE_RADIATION_PRESSURE //{
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            double iMass = 1.0/Phydro[i]->Mass;
            PhydroBody(i)->Vel[0] += Phydro[i]->MomentumRP[0]*iMass;
            PhydroBody(i)->Vel[1] += Phydro[i]->MomentumRP[1]*iMass;
            PhydroBody(i)->Vel[2] += Phydro[i]->MomentumRP[2]*iMass;
            Phydro[i]->MomentumRP[0] = Phydro[i]->MomentumRP[1] = Phydro[i]->MomentumRP[2] = 0.e0;
        }
    }
#endif // USE_RADIATION_PRESSURE //}

    return ;
}
