#include "config.h"
#include "IMFParameters.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "LifeTime.h"

static gsl_interp_accel *AccAccessAge[LIFETIME_Z];
static gsl_spline *SplineAccessAge[LIFETIME_Z];
static gsl_interp_accel *AccAccessMass[LIFETIME_Z];
static gsl_spline *SplineAccessMass[LIFETIME_Z];
static double LogLifeTimeArray[LIFETIME_Z][LIFETIME_M];

void InitializeLifeTimeLookUpTable(void){

    int i,j;
    double MassArray[LIFETIME_M];

    for(i=0;i<LIFETIME_Z;i++){
        AccAccessMass[i] = gsl_interp_accel_alloc();
        //SplineAccessMass[i] = gsl_spline_alloc(gsl_interp_cspline,LIFETIME_M);
        SplineAccessMass[i] = gsl_spline_alloc(gsl_interp_linear,LIFETIME_M);
        gsl_spline_init(SplineAccessMass[i], LifeTimeMass, LifeTimeZMF[i], LIFETIME_M);
        //gsl_spline_init(SplineAccessAge[i], LifeTimeMass, LifeTimeZMF[i],LIFETIME_M);
    }


    for(i=0;i<LIFETIME_M;i++){
        MassArray[i] = LifeTimeMass[LIFETIME_M-1-i];
        for(j=0;j<LIFETIME_Z;j++){
            LogLifeTimeArray[j][i] = log10(LifeTimeZMF[j][LIFETIME_M-1-i]);
        }
    }

    for(i=0;i<LIFETIME_Z;i++){
        AccAccessAge[i] = gsl_interp_accel_alloc();
        //SplineAccessAge[i] = gsl_spline_alloc(gsl_interp_cspline,LIFETIME_M);
        SplineAccessAge[i] = gsl_spline_alloc(gsl_interp_linear,LIFETIME_M);
        //gsl_spline_init(SplineAccessMass[i], MassArray, LogLifeTimeArray[i], LIFETIME_M);
        gsl_spline_init(SplineAccessAge[i],LogLifeTimeArray[i],MassArray,LIFETIME_M);
    }

    return;
}

/*  
 * This function returns "Just" Dying Star's Mass at the given Age.
 * If the given Age overs LookUp Table, this function returns "-1".  
 * "Age" must be in the unit of "YEAR".
 */
double GetDyingStarMassFromAge(double Age, const double Metal){

    if(Age < 0.e0){
        fprintf(stderr,"Age is %g\n",Age);
        return -1.0;
    }

    int i;
    //double LogAge = log10(Age)+10.0;
    //const static double ConvetTime = Pall.UnitTime/YEAR_CGS;  
    double LogAge = log10(Age);

    if(Metal <= LifeTimeZ[0]){
        if( (LogAge > LogLifeTimeArray[0][0])
                && (LogAge < LogLifeTimeArray[0][LIFETIME_M-1])){
            //if(get_myid() == 0)
                //fprintf(stderr,"Mass =  %g\n",gsl_spline_eval(SplineAccessAge[0], LogAge, AccAccessAge[0]));
            return (gsl_spline_eval(SplineAccessAge[0], LogAge, AccAccessAge[0]));
        }
    } else if(Metal >= LifeTimeZ[LIFETIME_Z-1]){
        if( (LogAge > LogLifeTimeArray[LIFETIME_Z-1][0])
                && (LogAge < LogLifeTimeArray[LIFETIME_Z-1][LIFETIME_M-1])){
            return gsl_spline_eval(SplineAccessAge[LIFETIME_Z-1], LogAge, AccAccessAge[LIFETIME_Z-1]);
        }
    } else {
        i=0;
        while(LifeTimeZ[i] < Metal) i++;
        if( (LogAge > MAX(LogLifeTimeArray[i-1][0],LogLifeTimeArray[i][0]))
                &&(LogAge < MIN(LogLifeTimeArray[i-1][LIFETIME_M-1],LogLifeTimeArray[i][LIFETIME_M-1])) ){
            return (((LifeTimeZ[i]-Metal)*gsl_spline_eval(SplineAccessAge[i-1], LogAge, AccAccessAge[i-1])
                        +(Metal-LifeTimeZ[i-1])*gsl_spline_eval(SplineAccessAge[i], LogAge, AccAccessAge[i]))/
                        (LifeTimeZ[i]-LifeTimeZ[i-1]));
        }
    }

    return -1.e0;
}

/*  
 * This function returns "Just" Dying Star's Mass at the given Age.
 * If the given Age overs LookUp Table, this function returns "-1".  
 * "Age" must be in the unit of "YEAR".
 * Moreover, this function returns Mmax or Mmin when age is out of the table defined range.
 */
double GetDyingStarMassFromAgeWithClip(const double Age, const double Metal){

    int i;
    //double LogAge = log10(Age)+10.0;
    //const static double ConvetTime = Pall.UnitTime/YEAR_CGS;  
    double LogAge = log10(Age);

    if(Metal <= LifeTimeZ[0]){
        if( (LogAge > LogLifeTimeArray[0][0])
                && (LogAge < LogLifeTimeArray[0][LIFETIME_M-1])){
            //if(get_myid() == 0)
                //fprintf(stderr,"Mass =  %g\n",gsl_spline_eval(SplineAccessAge[0], LogAge, AccAccessAge[0]));
            return (gsl_spline_eval(SplineAccessAge[0], LogAge, AccAccessAge[0]));
        } else if (LogAge < LogLifeTimeArray[0][0]){
            return (gsl_spline_eval(SplineAccessAge[0], LogLifeTimeArray[0][0], AccAccessAge[0]));
        } else if (LogAge > LogLifeTimeArray[0][LIFETIME_M-1]){
            return (gsl_spline_eval(SplineAccessAge[0], LogLifeTimeArray[0][LIFETIME_M-1], AccAccessAge[0]));
        }
    } else if(Metal >= LifeTimeZ[LIFETIME_Z-1]){
        if( (LogAge > LogLifeTimeArray[LIFETIME_Z-1][0])
                && (LogAge < LogLifeTimeArray[LIFETIME_Z-1][LIFETIME_M-1])){
            return gsl_spline_eval(SplineAccessAge[LIFETIME_Z-1], LogAge, AccAccessAge[LIFETIME_Z-1]);
        } else if (LogAge < LogLifeTimeArray[LIFETIME_Z-1][0]){
            return gsl_spline_eval(SplineAccessAge[LIFETIME_Z-1], LogLifeTimeArray[LIFETIME_Z-1][0], AccAccessAge[LIFETIME_Z-1]);
        } else if (LogAge > LogLifeTimeArray[LIFETIME_Z-1][LIFETIME_M-1]){
            return gsl_spline_eval(SplineAccessAge[LIFETIME_Z-1], LogLifeTimeArray[LIFETIME_Z-1][LIFETIME_M-1], AccAccessAge[LIFETIME_Z-1]);
        }
    } else {
        i=0;
        while(LifeTimeZ[i] < Metal) i++;
        if( (LogAge > MAX(LogLifeTimeArray[i-1][0],LogLifeTimeArray[i][0]))
                &&(LogAge < MIN(LogLifeTimeArray[i-1][LIFETIME_M-1],LogLifeTimeArray[i][LIFETIME_M-1])) ){
            return (((LifeTimeZ[i]-Metal)*gsl_spline_eval(SplineAccessAge[i-1], LogAge, AccAccessAge[i-1])
                        +(Metal-LifeTimeZ[i-1])*gsl_spline_eval(SplineAccessAge[i], LogAge, AccAccessAge[i]))/
                        (LifeTimeZ[i]-LifeTimeZ[i-1]));
        } else if (LogAge < MAX(LogLifeTimeArray[i-1][0],LogLifeTimeArray[i][0])){
            return (((LifeTimeZ[i]-Metal)*gsl_spline_eval(SplineAccessAge[i-1], MAX(LogLifeTimeArray[i-1][0],LogLifeTimeArray[i][0]), AccAccessAge[i-1])
                        +(Metal-LifeTimeZ[i-1])*gsl_spline_eval(SplineAccessAge[i], MAX(LogLifeTimeArray[i-1][0],LogLifeTimeArray[i][0]), AccAccessAge[i]))/
                        (LifeTimeZ[i]-LifeTimeZ[i-1]));
        } else if (LogAge > MIN(LogLifeTimeArray[i-1][LIFETIME_M-1],LogLifeTimeArray[i][LIFETIME_M-1])){
            return (((LifeTimeZ[i]-Metal)*gsl_spline_eval(SplineAccessAge[i-1], MIN(LogLifeTimeArray[i-1][LIFETIME_M-1],LogLifeTimeArray[i][LIFETIME_M-1]), AccAccessAge[i-1])
                        +(Metal-LifeTimeZ[i-1])*gsl_spline_eval(SplineAccessAge[i], MIN(LogLifeTimeArray[i-1][LIFETIME_M-1],LogLifeTimeArray[i][LIFETIME_M-1]), AccAccessAge[i]))/
                        (LifeTimeZ[i]-LifeTimeZ[i-1]));
        }
    }

    return -1.e0;
}

/*  
 * The range of the stellar mass, which can obtain "LifeTime" correctly, is 0.6 <-> 120.
 * The returned value is the LifeTime in yr.
 */
double GetLifeTimeFromStellarMass(const double Mass, const double Metal){

    int i;

    if(Metal <= LifeTimeZ[0]){
        return gsl_spline_eval(SplineAccessMass[0], Mass, AccAccessMass[0]);
    } else if(Metal >= LifeTimeZ[LIFETIME_Z-1]){
        return gsl_spline_eval(SplineAccessMass[LIFETIME_Z-1], Mass, AccAccessMass[LIFETIME_Z-1]);
    } else {
        i=0;
        while(LifeTimeZ[i] < Metal) i++;
        return (((LifeTimeZ[i]-Metal)*gsl_spline_eval(SplineAccessMass[i-1], Mass, AccAccessMass[i-1])
                +(Metal-LifeTimeZ[i-1])*gsl_spline_eval(SplineAccessMass[i], Mass, AccAccessMass[i]))/
                        (LifeTimeZ[i]-LifeTimeZ[i-1]));
    }
}

static double IMFNormalizedConstant;
static double InvIMFPower;

/* ! This function computes the nomalization of the IMF (Initial Mass Function) 
 */
void InitializeIMF(void){
    IMFNormalizedConstant = ((-ImfPower+1)/(pow(MIMFMax,-ImfPower+1)-pow(MIMFMin,-ImfPower+1)));
    InvIMFPower = 1.0/(-ImfPower);
    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"IMF normalized constant(%g) = %g\n",-ImfPower,IMFNormalizedConstant);
    return ;
}

double IMFSNIINumberPerMass(void){
    return (IMFNormalizedConstant*
            (pow(MSNIIMax,-ImfPower)-pow(MSNIIMin,-ImfPower))*InvIMFPower);
}

double IMFdSNIINumberPerMass(const double Mass){
    return (IMFNormalizedConstant*(pow(Mass,-ImfPower-1)));
}

// TimeMax and TimeMin are in unit of YEAR. 
double IMFSNIINumberPerMassLimitedTimeRange(const double TimeMax, const double TimeMin, const double Metallicity){
    double MassMax =  GetDyingStarMassFromAgeWithClip(TimeMin,Metallicity);
    double MassMin =  GetDyingStarMassFromAgeWithClip(TimeMax,Metallicity);
    if(MassMax < MSNIIMin) MassMax = MSNIIMin;
    if(MassMin < MSNIIMin) MassMin = MSNIIMin;
    if(MassMax > MSNIIMax) MassMax = MSNIIMax;
    if(MassMin > MSNIIMax) MassMin = MSNIIMax;

    return (IMFNormalizedConstant*
            (pow(MassMax,-ImfPower)-pow(MassMin,-ImfPower))*InvIMFPower);
}

double IMFSNIINumberPerMassLimitedRange(const double MassMax, const double MassMin){
    return (IMFNormalizedConstant*
            (pow(MassMax,-ImfPower)-pow(MassMin,-ImfPower))*InvIMFPower);
}

double IMFSNIIEjectaMassPerMass(void){
    return ( IMFNormalizedConstant*
            ((pow(MSNIIMax,-ImfPower+1)-pow(MSNIIMin,-ImfPower+1))/(-ImfPower+1)
             -1.4*(pow(MSNIIMax,-ImfPower)-pow(MSNIIMin,-ImfPower))*InvIMFPower));
}

double IMFSNIIEjectaMassPerMassLimitedRange(const double MassMax, const double MassMin){
    return ( IMFNormalizedConstant*
            ((pow(MassMax,-ImfPower+1)-pow(MassMin,-ImfPower+1))/(-ImfPower+1)
             -1.4*(pow(MassMax,-ImfPower)-pow(MassMin,-ImfPower))*InvIMFPower));
}

// TimeMax and TimeMin are in unit of YEAR. 
double IMFSNIIEjectaMassPerMassLimitedTimeRange(const double TimeMax, const double TimeMin, const double Metallicity){

    double MassMax =  GetDyingStarMassFromAgeWithClip(TimeMin,Metallicity);
    double MassMin =  GetDyingStarMassFromAgeWithClip(TimeMax,Metallicity);
    return ( IMFNormalizedConstant*
            ((pow(MassMax,-ImfPower+1)-pow(MassMin,-ImfPower+1))/(-ImfPower+1)
             -1.4*(pow(MassMax,-ImfPower)-pow(MassMin,-ImfPower))*InvIMFPower));
}

double IMFSNIIYieldMassPerMass(void){ // Steinmetz Muller 
    return ( IMFNormalizedConstant*
            (0.357*(pow(MSNIIMax,-ImfPower+1)-pow(MSNIIMin,-ImfPower+1))/(-ImfPower+1)
             -2.2*(pow(MSNIIMax,-ImfPower)-pow(MSNIIMin,-ImfPower))*InvIMFPower));
}

double IMFSNIIYieldMassPerMassLimitedRange(const double MassMax, const double MassMin){
    return ( IMFNormalizedConstant*
            (0.357*(pow(MassMax,-ImfPower+1)-pow(MassMin,-ImfPower+1))/(-ImfPower+1)
             -2.2*(pow(MassMax,-ImfPower)-pow(MassMin,-ImfPower))*InvIMFPower));
}

void IMFSNIINumberEjectaYieldPerMassLimitedRange(double MassMax, double MassMin,
        double *SNIINumber, double *SNIIEjectaMass, double *SNIIYieldMass){

    double coef = IMFNormalizedConstant*InvIMFPower;
    double mu = pow(MassMax,-ImfPower);
    double ml = pow(MassMin,-ImfPower);
    double mul = mu-ml;
    double coef1 = IMFNormalizedConstant/(-ImfPower+1.0);
    double mu1 = pow(MassMax,-ImfPower+1);
    double ml1 = pow(MassMin,-ImfPower+1);
    double mul1 = mu1-ml1;


    *SNIINumber = coef*mul; 
    *SNIIEjectaMass = coef1*mul1-coef*1.4*mul;
    *SNIIYieldMass = coef1*0.357*mul1-coef*2.2*mul;
    return;   
}

double IMFLongestLifeTime(void){
    int i;
    double MaxLifeTime =  LifeTimeZMF[0][20];

    for(i=0;i<LIFETIME_Z;i++)
        MaxLifeTime = MAX(MaxLifeTime, LifeTimeZMF[i][20]);
    return MaxLifeTime;
}

double IMFShortestLifeTime(void){
    int i;
    double MinLifeTime =  LifeTimeZMF[0][LIFETIME_M-1];

    for(i=0;i<LIFETIME_Z;i++)
        MinLifeTime = MIN(MinLifeTime, LifeTimeZMF[i][LIFETIME_M-1]);

    return MinLifeTime; 
}
