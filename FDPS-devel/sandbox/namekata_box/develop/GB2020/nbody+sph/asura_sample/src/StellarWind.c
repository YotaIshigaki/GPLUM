#include "config.h"
#include "StellarWindEnergy.h"

static double SW_lifetime;
static double dt_SW_max;

#define SWScaleFactor (0.002)

/*
 * This function returns an index which covers the Age;
 * StellarWindTime[index] < Age < StellarWindTime[index+1];
 */
static inline int __attribute__((always_inline)) StellarWindFindTargetAgeTableIndex(const double Age){

    if(Age < StellarWindTime[0]){
        return 0;
    } else if(Age >= StellarWindTime[StellarWindDataSize-2]){
        return StellarWindDataSize-2;
    } else {
        int LValue = 0;
        int RValue = StellarWindDataSize-1;
        int Pivot = (LValue+RValue)>>1;
        do {
            if(Age < StellarWindTime[Pivot]){
                RValue = Pivot;
            } else {
                LValue = Pivot;
            }

            Pivot = (LValue+RValue)>>1;
            if(RValue-LValue == 1)
                return LValue;

        } while((RValue-LValue) != 0);

    }

}

static inline int __attribute__((always_inline)) StellarWindFindTargetMetalTableIndex(const double Z){

    const int Size = sizeof(StellarWindZ)/sizeof(double);
    for(int i=0;i<Size-1;i++){
        if(Z < StellarWindZ[i]){
            //fprintf(stderr,"%g <  %g\n",Z,StellarWindZ[i]);
            return MAX(i-1,0);
        }
    }
    return Size-2;
}


static inline double __attribute__((always_inline)) LinearInterpolation(const double x, const double x1, const double x0, const double y1, const double y0){
    double grad = (y1-y0)/(x1-x0);
    return grad*(x-x0)+y0;
}

static double UpdateStellarWindEnergybyTableInterpolation(const int index){ // in ergs/Msun

    double Z = Pstar[index]->Z;
    const double Factor = Pall.UnitTime/YEAR_CGS;
    double Age1 = fmax(0.0,(Pall.TCurrent-Pstar[index]->FormationTime)*Factor);
    double Age0 = fmax(0.0,(Pall.TCurrent-Pstar[index]->FormationTime-PstarBody(index)->dt)*Factor);

    int TableID1 = StellarWindFindTargetAgeTableIndex(Age1);
    int TableID0 = StellarWindFindTargetAgeTableIndex(Age0);

    double SWEnergy = 0.e0;
    if(Z<=StellarWindZ[0]){
        SWEnergy = LinearInterpolation(Age1,
                     StellarWindTime[TableID1+1],StellarWindTime[TableID1],
                     StellarWindEnergy[0][TableID1+1],StellarWindEnergy[0][TableID1])
                  -LinearInterpolation(Age0,
                     StellarWindTime[TableID0+1],StellarWindTime[TableID0],
                     StellarWindEnergy[0][TableID0+1],StellarWindEnergy[0][TableID0]);
    } else if(Z>=StellarWindZ[7]){
        SWEnergy = LinearInterpolation(Age1,
                     StellarWindTime[TableID1+1],StellarWindTime[TableID1],
                     StellarWindEnergy[7][TableID1+1],StellarWindEnergy[7][TableID1])
                  -LinearInterpolation(Age0,
                     StellarWindTime[TableID0+1],StellarWindTime[TableID0],
                     StellarWindEnergy[7][TableID0+1],StellarWindEnergy[7][TableID0]);
    } else { // interpolation in Z
        int MetalTableID = StellarWindFindTargetMetalTableIndex(Z);

        double SWEnergy0 = LinearInterpolation(Age1,
                            StellarWindTime[TableID1+1],StellarWindTime[TableID1],
                            StellarWindEnergy[MetalTableID][TableID1+1],StellarWindEnergy[MetalTableID][TableID1])
                          -LinearInterpolation(Age0,
                            StellarWindTime[TableID0+1],StellarWindTime[TableID0],
                            StellarWindEnergy[MetalTableID][TableID0+1],StellarWindEnergy[MetalTableID][TableID0]);

        double SWEnergy1 = LinearInterpolation(Age1,
                            StellarWindTime[TableID1+1],StellarWindTime[TableID1],
                            StellarWindEnergy[MetalTableID+1][TableID1+1],StellarWindEnergy[MetalTableID+1][TableID1])
                          -LinearInterpolation(Age0,
                            StellarWindTime[TableID0+1],StellarWindTime[TableID0],
                            StellarWindEnergy[MetalTableID+1][TableID0+1],StellarWindEnergy[MetalTableID+1][TableID0]);

        SWEnergy = LinearInterpolation(Z,
                    StellarWindZ[MetalTableID+1],StellarWindZ[MetalTableID],
                    SWEnergy1,SWEnergy0);
    }

    return SWEnergy*1.e51; 
}

void UpdateStellarWindEnergy(const int index){ // in simulation units

    const double MassConversionFactor = Pall.UnitMass/MSUN_CGS;
    const double TimeConversionFactor = Pall.UnitTime/YEAR_CGS;
#ifdef  USE_STELLAR_WIND //{

    Pstar[index]->StellarWindEnergy = 
#if STELLAR_WIND_TYPE == 0
        SWScaleFactor*MassConversionFactor*Pstar[index]->InitialMass*
        ASRFLXGetBol((Pall.TCurrent-Pstar[index]->FormationTime)*TimeConversionFactor,
                Pstar[index]->Z)*PstarBody(index)->dt*Pall.UnitTime; 
        // Note: ASRFLXGetBol returns the specific bolometric luminosity of an SSP particle.
        // The unit is erg/s/Msun.
        // To convert the stellar wind energy, we need to multiple dt.
        // Since StellarWindEnergy holds energy in units of erg, we need to
        // multiple dt*Pall.UnitTime.
        // The scale SWScaleFactor comes from figure 1 in Agertz+(2013)
#if 0
    fprintf(stderr,"SW: %d Mass = %g [Msun], Lb = %g [erg/s], dt %g [yr], SWE %g [erg]\n",index,
        MassConversionFactor*Pstar[index]->InitialMass,
        MassConversionFactor*Pstar[index]->InitialMass*
        ASRFLXGetBol((Pall.TCurrent-Pstar[index]->FormationTime)*TimeConversionFactor,
                Pstar[index]->Z),
        PstarBody(index)->dt*Pall.UnitTime/YEAR_CGS,
        SWScaleFactor*MassConversionFactor*Pstar[index]->InitialMass*
        ASRFLXGetBol((Pall.TCurrent-Pstar[index]->FormationTime)*TimeConversionFactor,
                Pstar[index]->Z)*PstarBody(index)->dt*Pall.UnitTime);
#endif

#elif STELLAR_WIND_TYPE == 1 
         MassConversionFactor*Pstar[index]->InitialMass*
            UpdateStellarWindEnergybyTableInterpolation(index);
#else 
#error Incrrect STELLAR_WIND_TYPE 
#endif
#endif // USE_STELLAR_WIND //}

    return ;
}


static double CheckStellarWindEnergybyTableInterpolation(const double Z, const double Age_in_year){ // in ergs/Msun

    int TableID = StellarWindFindTargetAgeTableIndex(Age_in_year);

    double SWEnergy = 0.e0;
    if(Z<=StellarWindZ[0]){
        SWEnergy = LinearInterpolation(Age_in_year,
                     StellarWindTime[TableID+1],StellarWindTime[TableID],
                     StellarWindEnergy[0][TableID+1],StellarWindEnergy[0][TableID]);
    } else if(Z>=StellarWindZ[7]){
        SWEnergy = LinearInterpolation(Age_in_year,
                     StellarWindTime[TableID+1],StellarWindTime[TableID],
                     StellarWindEnergy[7][TableID+1],StellarWindEnergy[7][TableID]);
    } else { // interpolation in Z
        int MetalTableID = StellarWindFindTargetMetalTableIndex(Z);

        double SWEnergy0 = LinearInterpolation(Age_in_year,
                            StellarWindTime[TableID+1],StellarWindTime[TableID],
                            StellarWindEnergy[MetalTableID][TableID+1],StellarWindEnergy[MetalTableID][TableID]);

        double SWEnergy1 = LinearInterpolation(Age_in_year,
                            StellarWindTime[TableID+1],StellarWindTime[TableID],
                            StellarWindEnergy[MetalTableID+1][TableID+1],StellarWindEnergy[MetalTableID+1][TableID]);

        SWEnergy = LinearInterpolation(Z,
                    StellarWindZ[MetalTableID+1],StellarWindZ[MetalTableID],
                    SWEnergy1,SWEnergy0);
    }

    return SWEnergy; 
}

static void WriteStellarWindDataTest(void){

    const double Z[] = {0,0.00002,0.0002,0.002,0.01,0.03};

    // Random
    double Size = sizeof(Z)/sizeof(double);
    for(int k=0;k<Size;k++){
        FILE *fp;
        char fname[MaxCharactersInLine];
        Snprintf(fname,"./Log/StellarWindTest.%03d",k);
        FileOpen(fp,fname,"w");
        fprintf(fp,"#Z=%g\n",Z[k]);

        int TestSize = 10000; 
        double Tmin = 0;   //year
        double Tmax = 3e7; //year
        for(int i=0;i<TestSize;i++){

            double Age = Tmax*gsl_rng_uniform(RandomGenerator);
            double e = CheckStellarWindEnergybyTableInterpolation(Z[k],Age);
            fprintf(fp,"%g %g\n",Age,e);
        }

        fclose(fp);
    }

    return ;
}


static void WriteStellarWindData(void){

    FILE *fp;
    FileOpen(fp,"./Log/StellarWind.dat","w");

    int size = sizeof(StellarWindTime)/sizeof(double);
    dprintlmpi(size);

    for(int i=0;i<size;i++){
        fprintf(fp,"%g %g %g %g %g %g %g %g %g\n",
                StellarWindTime[i],
                StellarWindEnergy[0][i], //Z=0
                StellarWindEnergy[1][i], //Z=0.00001
                StellarWindEnergy[2][i], //Z=0.0001
                StellarWindEnergy[3][i], //Z=0.001
                StellarWindEnergy[4][i], //Z=0.004
                StellarWindEnergy[5][i], //Z=0.008
                StellarWindEnergy[6][i], //Z=0.02
                StellarWindEnergy[7][i]  //Z=0.05
                );
    }

    fclose(fp);

    return ;
}

void InitStellarWind(void){

    SW_lifetime = STELLAR_WIND_SOURCE_LIFETIME/Pall.UnitTime;
    dt_SW_max = STELLAR_WIND_MAX_TIMESTEP/Pall.UnitTime;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Init StellarWind\n");
        fprintf(stderr," SW_lifetime = %g\n",SW_lifetime);
        fprintf(stderr," dt_SW_max = %g\n",dt_SW_max);

        // WriteStellarWindData
        WriteStellarWindData();
        WriteStellarWindDataTest();
    }


    return;
}

double GetStellarWindSourceLifeTime(void){
    return SW_lifetime;
}

double GetStellarWindMaxTimestep(void){
    return dt_SW_max;
}

bool CheckStellarWindSource(const int index){

    if((Pstar[index]->TypeII == false)
            &&(Pall.TCurrent-Pstar[index]->FormationTime < GetStellarWindSourceLifeTime())){
        return true;
    } else {
        return false;
    }
}

bool CheckStellarWindSourceThroughPbody(const int index){
    if((PbodyStar(index)->TypeII == false)
            &&(Pall.TCurrent-PbodyStar(index)->FormationTime < GetStellarWindSourceLifeTime())){
        return true;
    } else {
        return false;
    }
}

int CountStellarWindSourceNumber(void){

    int counter = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarBody(i)->Active){
            if(CheckStellarWindSource(i)){
                counter ++;
            }
        }
    }
    return counter;
}
