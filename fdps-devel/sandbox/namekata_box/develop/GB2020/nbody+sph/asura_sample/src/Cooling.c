#include "config.h"
#include "Cooling.h"
#include "Heating.h"
#include "H2.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#ifdef USE_SPAANS1997_COOLING_FUNCTIONS
#include "CoolingTable.h"
static void InitializeCoolingTableSpaans1997(void);
static double CoolingFunction(double T, const int mode);
static void MakeCoolingTable(const int mode);
static void CalcCoolingSpaans1997(void);
#endif // USE_SPAANS1997_COOLING_FUNCTIONS
#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
#include "CoolingTableSpaans2008.h"
static void InitializeCoolingTableSpaans2008(void);
#endif // USE_SPAANS2008_COOLING_FUNCTIONS
#ifdef USE_CLOUDY_COOLING_FUNCTIONS
#include "CloudyCooling.h"
static void CalcCoolingCloudy(void);
#endif // USE_CLOUDY_COOLING_FUNCTIONS

static char PathCoolingHeatingLog[] = "./CoolingHeatingLog"; 
#ifdef USE_INVERSE_COMPTON_COOLING
static double InverseComptonCooling(const double Ti, const double ne, const double redshift);
#endif // USE_INVERSE_COMPTON_COOLING

static double logTmin,logTmax;
static double sqrt11,isqrt11;

static double T_cut;
static double InitialSkipCoolingRoutine;
static double FUVConversionFactor; // [erg/g/s = cm^2/s^3]

void InitializeCoolingTable(void){

#ifdef COOLING_RUN //{
    MakeDir(PathCoolingHeatingLog);
    Pall.ConvertCoolingRate = GetUnitCoolingRate();
    FUVConversionFactor = GetUnitCoolingRate(); // [erg/g/s = cm^2/s^3]

#if defined(USE_ASRCH)
    if(MPIGetMyID() == MPI_ROOT_RANK){
        ASRCHShowInfo();
        ASRCHShowDimension();
#ifndef USE_EXACT_COOLING_SOLVER //{
        fprintf(stderr,"ASRCH always uses the exact cooling solver\n");
#endif //USE_EXACT_COOLING_SOLVER //}
        fflush(NULL);
    }

    ASRCHSetRunParameterReionizationRedshift(ASRCH_REIONIZATION_REDSHIFT);
    ASRCHSetRunParameterVerbose(false);





#elif defined(USE_CLOUDY_COOLING_FUNCTIONS)

#if (defined(USE_SPAANS1997_COOLING_FUNCTIONS)||defined(USE_SPAANS2008_COOLING_FUNCTIONS))
#error Select one of the three cooling functions.
#endif
    if(MPIGetMyID() == MPI_ROOT_RANK){
        ReportCloudyCoolingFunction();
        ReportDimensionofCloudyCoolingFunction();
        fflush(NULL);
    }

#elif defined(USE_SPAANS2008_COOLING_FUNCTIONS) // USE_*_COOLING_FUNCTIONS

#if (defined(USE_SPAANS1997_COOLING_FUNCTIONS)||defined(USE_CLOUDY_COOLING_FUNCTIONS))
#error Select one of the three cooling functions.
#endif
    InitializeCoolingTableSpaans2008();
#elif defined(USE_SPAANS1997_COOLING_FUNCTIONS) // USE_*_COOLING_FUNCTIONS

#if (defined(USE_SPAANS2008_COOLING_FUNCTIONS)||defined(USE_CLOUDY_COOLING_FUNCTIONS))
#error Select one of the three cooling functions.
#endif
    InitializeCoolingTableSpaans1997();

#else // USE_*_COOLING_FUNCTIONS
#error Select one of the three cooling functions.
#endif // USE_*_COOLING_FUNCTIONS


#ifdef COOLING_CUTOFF_TEMPERATURE //{
    T_cut = log10(COOLING_CUTOFF_TEMPERATURE);
#endif // COOLING_CUTOFF_TEMPERATURE //}

#ifdef INITIAL_SKIP_COOLING_ROUTINE //{
    InitialSkipCoolingRoutine = INITIAL_SKIP_COOLING_ROUTINE/Pall.UnitTime;
#else  // INITIAL_SKIP_COOLING_ROUTINE //{//}
    InitialSkipCoolingRoutine = 0.e0;
#endif // INITIAL_SKIP_COOLING_ROUTINE //}

    sqrt11 = sqrt(1.1);
    isqrt11 = 1.e0/sqrt11;


#endif //COOLING_RUN //}
    return;
}


static inline double __attribute__((always_inline)) EvaluateColumnDensity(const int index){

    double X = 1.e0 - HeliumAbandance - Phydro[index]->Z;
    double nH = Pall.ConvertNumberDensityToCGS*Phydro[index]->Rho*X;

    //return pow(nH,2.0/3.0)*1.4e20;
    //return pow(nH,2.0/3.0)*3e19;
#ifdef USE_FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE //{
    return nH*FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE;
#else // USE_FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE //}//{
    // 2 is a calibration factor.
    return nH*(Phydro[index]->Rho/fabs(2*Phydro[index]->GradRho))*Pall.UnitLength;
#endif // USE_FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE //}

}

#ifdef USE_ASRCH //{
static double CoolingSolverASRCH(const int index){

	int iter = 0;
    double dt = Phydro[index]->dt_hydro;
#ifdef SET_SNII_TEMPERATURE_GAS //{
    double R = gsl_rng_uniform(RandomGenerator);
    double USN = Phydro[index]->DQheat/PhydroMass(index);
    double TSN = Pall.ConvertUtoT*USN;
    double Heating = 0.0;
    if(R < TSN/SNII_TEMPERATURE){
	    Heating = Pall.ConvertTtoU*SNII_TEMPERATURE/dt;
    } else {
        Heating = 0;
    }
#else //}//{
	double Heating = Phydro[index]->DQheat/(PhydroMass(index)*dt); // [erg/g/s]
#endif //}

    double U = Phydro[index]->U;
    double Uold = U, Uinit = U, Uupper = U, Ulower = U;
	double Metallicity = Phydro[index]->Z;
    double Rho = Phydro[index]->Rho;
    double X = 1.e0 - HeliumAbandance - Metallicity;
    double nH = Pall.ConvertNumberDensityToCGS*Rho*X;

#ifdef USE_FUVFEEDBACK //{

#define FUV_Epsilon (0.05)
#define FUV_Zsolar (0.0134)
#   ifdef USE_FUVFEEDBACK_CONSTANT_VALUE //{
    //double nH = Pall.ConvertNumberDensityToCGS*Rho;
    double Gamma_Heating = 1.e-24*FUV_Epsilon*FUVFEEDBACK_CONSTANT_VALUE*(Metallicity/FUV_Zsolar); // [erg/s]
    double FUVHeating = FUVConversionFactor*nH*Gamma_Heating; // [erg/g/s]
#   else // USE_FUVFEEDBACK_CONSTANT_VALUE //}//{

    double Ncol = EvaluateColumnDensity(index);

    double Zdep = (Metallicity/FUV_Zsolar);
#       ifdef USE_FUVFEEDBACK_LOCAL_CORRECTION //{
#           ifdef USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //{
    double G = (Phydro[index]->G0thin-Phydro[index]->G0thinLocal)*exp(-Zdep*FUVFEEDBACK_SIGMA1000*Ncol); //
#           else // USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //}//{
    double G = (Phydro[index]->G0thin-Phydro[index]->G0thinLocal)*exp(-FUVFEEDBACK_SIGMA1000*Ncol); //
#           endif // USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //}
    G += Phydro[index]->G0extLocal;
#       else // USE_FUVFEEDBACK_LOCAL_CORRECTION //}//{
#           ifdef USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //{
    double G = Phydro[index]->G0thin*exp(-Zdep*FUVFEEDBACK_SIGMA1000*Ncol); //
#           else // USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //}//{
    double G = Phydro[index]->G0thin*exp(-FUVFEEDBACK_SIGMA1000*Ncol); //
#           endif // USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //}
#       endif // USE_FUVFEEDBACK_LOCAL_CORRECTION //}

    double Gamma_Heating = 1.e-24*FUV_Epsilon*G*Zdep; // [erg/s]
    double FUVHeating = FUVConversionFactor*nH*Gamma_Heating/Rho; // [erg/g/s]
#   endif // USE_FUVFEEDBACK_CONSTANT //}

#else // USE_FUVFEEDBACK //}//{

    double Gamma_Heating = 0.e0; // [erg/s]
    double FUVHeating = 0.e0; 
#endif // USE_FUVFEEDBACK //}

#if (SingleStepCoolingSkip) //{
    if(Phydro[index]->DQheat > 0.e0){
        Phydro[index]->DQheat = 0.e0;
        return Heating;
    } else {
        Heating = 0.0;
    }
#endif // SingleStepCoolingSkip //}

    double dudt_exact = 0.0;
    double LognH = log10(nH);
    double LogT = log10(Pall.ConvertUtoT*U);
    Phydro[index]->DQheat = 0.e0;
    const double factor = CUBE(Pall.UnitTime)/SQ(Pall.UnitLength); // [erg/g/s], [cm^2/s^3]
    return factor*ASRCHGetCoolingHeatingValueExactSolverWithFUVHeating(LognH,LogT,
            Metallicity,dt*Pall.UnitTime,Gamma_Heating/nH);
}

static void CalcCoolingASRCH(void){

    ASRCHMakeRedshiftInterpolateData(Pall.Redshift);

    int counter = 0;

    double EnergyLossThisStep = 0.e0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            Phydro[i]->DuCooling = CoolingSolverASRCH(i);
            double dt = Phydro[i]->dt_hydro;

#define USE_BOTTOM_TEMPERATURE
#ifdef USE_BOTTOM_TEMPERATURE //{
            if(Pall.ConvertUtoT*Phydro[i]->U < 10.0){
                Phydro[i]->DuCooling = Pall.ConvertTtoU*(11.0-Pall.ConvertUtoT*Phydro[i]->U)/dt;
            }
#endif // USE_BOTTOM_TEMPERATURE //}

#define USE_TOP_TEMPERATURE
#ifdef USE_TOP_TEMPERATURE //{
            if(Phydro[i]->U+Phydro[i]->DuCooling*dt > Pall.ConvertTtoU*1.e9){
                Phydro[i]->DuCooling = (Pall.ConvertTtoU*1.e9-Phydro[i]->U)/dt;
            }
#endif // USE_TOP_TEMPERATURE //}

            Phydro[i]->UPred = Phydro[i]->U += Phydro[i]->DuCooling*dt;
            assert(Phydro[i]->U>0.e0);

            EnergyLossThisStep -= PhydroMass(i)*Phydro[i]->DuCooling*dt;
        }
    }
    Pall.CoolingEnergyLoss += EnergyLossThisStep;

    return ;
}
#endif // USE_ASRCH //}

#ifdef USE_MULTIPHASE_MODEL //{
/*
 * This is a function which evaluates uhot.
 * This function is kindly provided by T. Okamoto.
 */
static inline double __attribute__((always_inline)) find_uhot(double t_cool_in_year, double nH_in_cgs, double u_in_simulation_units){ // t_cool (yr), nH (/cc), u (code units)
    const double Tfactor = pow(3.16227766e7/Pall.ConvertUtoT, 1.0/3.0);
    return pow((t_cool_in_year/2.7e7)*nH_in_cgs*u_in_simulation_units, 2.0/3.0) * Tfactor;
}

static inline double __attribute__((always_inline)) GetGammaFUV(const double nH, const double Metallicity, const double Ncol, const double G0thin, const double G0thinLocal, const double G0extLocal){
    
#define FUV_Epsilon (0.05)
#define FUV_Zsolar (0.0134)

#ifdef USE_FUVFEEDBACK //{

#   ifdef USE_FUVFEEDBACK_CONSTANT_VALUE //{
    return 1.e-24*FUV_Epsilon*FUVFEEDBACK_CONSTANT_VALUE*(Metallicity/FUV_Zsolar); // [erg/s]
    // double FUVHeating = FUVConversionFactor*nH*Gamma_Heating; // [erg/g/s]
#   else // USE_FUVFEEDBACK_CONSTANT_VALUE //}//{

    double Zdep = (Metallicity/FUV_Zsolar);
#       ifdef USE_FUVFEEDBACK_LOCAL_CORRECTION //{
#           ifdef USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //{
    double G = (G0thin-G0thinLocal)*exp(-Zdep*FUVFEEDBACK_SIGMA1000*Ncol); //
#           else // USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //}//{
    double G = (G0thin-G0thinLocal)*exp(-FUVFEEDBACK_SIGMA1000*Ncol); //
#           endif // USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //}
    G += G0extLocal;
#       else // USE_FUVFEEDBACK_LOCAL_CORRECTION //}//{
#           ifdef USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //{
    double G = G0thin*exp(-Zdep*FUVFEEDBACK_SIGMA1000*Ncol); //
#           else // USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //}//{
    double G = G0thin*exp(-FUVFEEDBACK_SIGMA1000*Ncol); //
#           endif // USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //}
#       endif // USE_FUVFEEDBACK_LOCAL_CORRECTION //}
    return 1.e-24*FUV_Epsilon*G*Zdep; // [erg/s]
    //double FUVHeating = FUVConversionFactor*nH*Gamma_Heating/Rho; // [erg/g/s]

#   endif // USE_FUVFEEDBACK_CONSTANT //}

#else // USE_FUVFEEDBACK //}//{
    return 0.e0; // [erg/s]
#endif // USE_FUVFEEDBACK //}

#undef FUV_Epsilon 
#undef FUV_Zsolar 

}

#   if defined(USE_ASRCH)

static double MultiphaseCoolingSolverASRCH(const int index){

    const double factor = CUBE(Pall.UnitTime)/SQ(Pall.UnitLength); // [erg/g/s], [cm^2/s^3]
    double dt = Phydro[index]->dt_hydro;
    double Mass = Phydro[index]->Mass;
    double InvMass = 1.0/Mass;

    double U = Phydro[index]->U;
    U += Phydro[index]->DQheat*InvMass;
    Phydro[index]->U = Phydro[index]->UPred = U;

    // Evaluate the sound crossing time
    double cs_effective = sqrt(Pall.GGm1*U);
    double t_cross = 2.0*Phydro[index]->Kernel/cs_effective;
    const double t_cross_factor = Pall.UnitTime/YEAR_CGS;
    double t_cross_in_year = t_cross*t_cross_factor;

    double Rho = Phydro[index]->Rho;
	double Metallicity = Phydro[index]->Z;
    double X = 1.e0 - HeliumAbandance - Metallicity;
    double nH = Pall.ConvertNumberDensityToCGS*Rho*X;

    double Ehot = Phydro[index]->Mhot*Phydro[index]->Uhot + Phydro[index]->DQheat;

    if(U*Mass <= Ehot){
        fprintf(stderr,"Oh! E<Ehot %d %ld %g %g %g | %g %g %g %g\n",__LINE__,
                PhydroBody(index)->GlobalID,
                U*Mass,U,Mass,
                Ehot,
                Phydro[index]->Uhot,
                Phydro[index]->Mhot,
                Phydro[index]->DQheat);
        fflush(NULL);
        assert(U*Mass > Ehot);
    }

    double dudt_hot = 0.e0;
    double Mhot = 0.e0;
    double Uhot = 0.e0;
    double Mhot_prev = 0.e0;
    double Uhot_prev = 0.e0;
    Phydro[index]->MultiphaseFlag = false;
    if(Ehot > 0.e0){
#if 1
        Uhot = Phydro[index]->Uhot;
        Mhot = Phydro[index]->Mhot;
        if(Phydro[index]->DQheat > 0.0){
            Uhot = fmax(find_uhot(MULTIPHASE_MODEL_CRITERION*t_cross_in_year,nH,U),Uhot);
            Phydro[index]->MultiphaseEndTime = Pall.TCurrent+MULTIPHASE_MODEL_MAX_LIFETIME*t_cross;

            if(Uhot > 0){
                Mhot = Ehot/Uhot;
            }

#if 0
            if(Mhot > Mass){
                fprintf(stderr,"Oh! Mhot>Mass -%d- %ld | %g %g| %g %g\n",__LINE__,
                        PhydroBody(index)->GlobalID,
                        Mass,Mhot,
                        U*Mass,Ehot);
                fprintf(stderr,"Oh! Mhot>Mass (Uhot/Thot)  -%d- %g %g | (U/T) %g %g | \n",__LINE__,
                        Uhot,Pall.ConvertUtoT*Uhot,U,Pall.ConvertUtoT*U);

                fprintf(stderr,"Oh! Mhot>Mass Condition | %g < %g || Uh/U %g %g \n",
                        Uhot,Phydro[index]->U*(Uhot*Mhot)/(Phydro[index]->U*Phydro[index]->Mass),
                        Uhot,Phydro[index]->U);
                fflush(NULL);
            }
#endif

            if(Mhot > Mass){
                // This case, Uhot <= U. Thus, the multiphase mode should be
                // turned off.
                // fprintf(stderr,"Oh! Mhot>Mass!!!! -%d- %ld | Mass/Mhot %g %g\n",__LINE__,
                        // PhydroBody(index)->GlobalID,Mass,Mhot);

                Phydro[index]->Mhot = Phydro[index]->Rhohot = 0.e0;
                Uhot = Mhot = 0.e0;
            } else if(Mhot < 0.e0){
                // fprintf(stderr,"Oh! Mhot<0.0!!!! -%d- %ld | Mass/Mhot %g %g\n",__LINE__,
                        // PhydroBody(index)->GlobalID,Mass,Mhot);

                Phydro[index]->Uhot = Phydro[index]->Rhohot = 0.e0;
                Uhot = Mhot = 0.e0;
            } else if(Uhot < U){ // Need care for the case Uhot < U
                // fprintf(stderr,"Oh! Uhot<U!!!! -%d- %ld | U/Uhot %g %g\n",__LINE__,
                        // PhydroBody(index)->GlobalID,U,Uhot);
                Phydro[index]->Uhot = Phydro[index]->Mhot = Phydro[index]->Rhohot = 0.e0;
                Uhot = Mhot = 0.e0;
            }
            Mhot_prev = Mhot;
            Uhot_prev = Uhot;
        }

        if(Uhot > 0.e0){
            Mhot_prev = Mhot;
            Uhot_prev = Uhot;

            double Rhohot = Rho*U/Uhot;
            double nHhot = Pall.ConvertNumberDensityToCGS*Rhohot*X;
            double LognHhot = log10(nHhot);
            double LogThot  = log10(Pall.ConvertUtoT*Uhot);
            dudt_hot = factor*ASRCHGetCoolingHeatingValueExactSolverWithFUVHeating(LognHhot
                    ,LogThot, Metallicity,dt*Pall.UnitTime,0.e0);
            if(isnan(dudt_hot)){
                fprintf(stderr,"Oh! dudt_hot is nan! -%d- %g %g %g %g\n",__LINE__,
                        LognHhot,LogThot,Metallicity,dt*Pall.UnitTime);
                fflush(NULL);
            }
            Uhot += dudt_hot*dt;
            if(Uhot < 0.e0){
                Uhot = 0.e0;
                Mhot = 0.e0;
                dudt_hot = -Uhot_prev/dt;
                Phydro[index]->MultiphaseFlag = false;
                fprintf(stderr,"Use the special condition\n");
            }

            // Integrate Uhot here
            Phydro[index]->Uhot = Uhot;
            Phydro[index]->Mhot = Mhot;
            Phydro[index]->Rhohot = Rhohot;
            Phydro[index]->MultiphaseFlag = true;
        }
#else
        Uhot = find_uhot(MULTIPHASE_MODEL_CRITERION*t_cross_in_year,nH,U);

        if(Uhot > 0 ){
            Mhot = Ehot/Uhot;
        }

        if(Mhot > Mass){
            Phydro[index]->Mhot = Phydro[index]->Rhohot = 0.e0;
        } else if(Mhot < 0.e0){
            //Phydro[index]->Uhot = Phydro[index]->Mhot = Phydro[index]->Rhohot = 0.e0;
            Phydro[index]->Uhot = Phydro[index]->Rhohot = 0.e0;
        } else if(Uhot < U){ // Need care for the case Uhot < U
            Phydro[index]->Uhot = Phydro[index]->Mhot = Phydro[index]->Rhohot = 0.e0;
        } else {
            Mhot_prev = Mhot;
            Uhot_prev = Uhot;

            double Rhohot = Rho*U/Uhot;
            double nHhot = Pall.ConvertNumberDensityToCGS*Rhohot*X;
            double LognHhot = log10(nHhot);
            double LogThot  = log10(Pall.ConvertUtoT*Uhot);
            dudt_hot = factor*ASRCHGetCoolingHeatingValueExactSolverWithFUVHeating(LognHhot
                    ,LogThot, Metallicity,dt*Pall.UnitTime,0.e0);
            if(isnan(dudt_hot)){
                fprintf(stderr,"Oh! dudt_hot is nan! %g %g %g %g\n",
                        LognHhot,LogThot,Metallicity,dt*Pall.UnitTime);
                fflush(NULL);
            }
            Uhot += dudt_hot*dt;
            if(Uhot < 0.e0){
                Uhot = 0.e0;
                Mhot = 0.e0;
                dudt_hot = -Uhot_prev/dt;
                Phydro[index]->MultiphaseFlag = false;
                fprintf(stderr,"Use special condition\n");
            }
            // Integrate Uhot here
            Phydro[index]->Uhot = Uhot;
            Phydro[index]->Mhot = Mhot;
            Phydro[index]->Rhohot = Rhohot;
            Phydro[index]->MultiphaseFlag = true;
        }
#endif
    } else { // For a non-multiphase particle.
        Phydro[index]->Uhot = Phydro[index]->Mhot = Phydro[index]->Rhohot = 0.e0;
    }


    double Mcold = (Mass-Mhot_prev);
    if(Mcold <= 0.e0){
        fprintf(stderr,"Oh! Mcold<=0 -%d- %ld | %g %g %g %d\n",__LINE__,
                PhydroBody(index)->GlobalID,
                Mass,Mcold,Mhot_prev,Phydro[index]->MultiphaseFlag?1:0);
        fflush(NULL);
        assert(Mcold > 0.e0);
    }
    double Ucold = (Mass*U-Mhot_prev*Uhot_prev)/Mcold;
    if(Ucold <= 0.e0){
        fprintf(stderr,"Oh! Ucold<=0 -%d- %ld | %g %g %g\n",__LINE__,
                PhydroBody(index)->GlobalID,
                U,Ucold,Phydro[index]->Uhot);
        fflush(NULL);
        assert(Ucold > 0.e0);
    }
    double Rhocold = Rho*U/Ucold;
    double nHcold = Pall.ConvertNumberDensityToCGS*Rhocold*X;
    double LognHcold = log10(nHcold);
    double LogTcold  = log10(Pall.ConvertUtoT*Ucold);
    double Gamma_Heating_FUV = GetGammaFUV(nHcold,Metallicity,EvaluateColumnDensity(index),
                Phydro[index]->G0thin,Phydro[index]->G0thinLocal,Phydro[index]->G0extLocal);
    
    double dudt_cold = factor*ASRCHGetCoolingHeatingValueExactSolverWithFUVHeating(LognHcold
            ,LogTcold, Metallicity,dt*Pall.UnitTime,Gamma_Heating_FUV/nHcold);

    if(isnan(dudt_cold)){
        fprintf(stderr,"Oh! dudt_cold is nan! -%d- %g %g %g %g %g | %g %g %g %g | %g %g\n",__LINE__,
                LognHcold,LogTcold,Metallicity,dt*Pall.UnitTime,Gamma_Heating_FUV/nHcold,
                nHcold,Ucold,Mass,Mhot_prev,Mass*U,Mhot_prev*Uhot_prev);

        StructureReportPhydro(index);
        fflush(NULL);
    }

    Ucold += dudt_cold*dt;


    Phydro[index]->DQheat = 0.e0;
    double Unew = (Mcold*Ucold+Mhot*Uhot)*InvMass;
    double dudt = (Unew-U)/dt;

    if(Pall.ConvertUtoT*(U+dudt*dt) >= 1.e9){
        fprintf(stderr,"Strange! -%d- %ld %g | %g %g | %g %g %g %g\n",__LINE__,
                PhydroBody(index)->GlobalID,
                Pall.ConvertUtoT*(U+dudt*dt),
                Pall.ConvertUtoT*(U),
                Pall.ConvertUtoT*(dudt*dt),
                Mcold,dudt_cold,Mhot_prev,dudt_hot);
        fflush(NULL);
    }

    if(Phydro[index]->MultiphaseFlag){
        if(Phydro[index]->MultiphaseEndTime < Pall.TCurrent){
            Phydro[index]->Uhot = Phydro[index]->Mhot = Phydro[index]->Rhohot = 0.e0;
            Phydro[index]->MultiphaseFlag = false;
        // } else if(Phydro[index]->Uhot < 2*Phydro[index]->U){
            // Phydro[index]->Uhot = Phydro[index]->Mhot = Phydro[index]->Rhohot = 0.e0;
            // Phydro[index]->MultiphaseFlag = false;
        } else {

            double LognH = log10(nH);
            double LogT  = log10(Pall.ConvertUtoT*(U+dudt*dt));
            const double InvYEAR_CGS = 1.0/YEAR_CGS;
            double t_cool_in_year = ASRCHGetCoolingTime(LognH,LogT,Metallicity)*InvYEAR_CGS;
            if(t_cool_in_year > 10*t_cross_in_year){
                 fprintf(stderr,"End condition t_cool = %g, 10*t_cross = %g\n",
                        t_cool_in_year,10*t_cross_in_year);

                double Rhohot = Rho*U/Uhot_prev;
                double nHhot = Pall.ConvertNumberDensityToCGS*Rhohot*X;
                double Thot  = (Pall.ConvertUtoT*Uhot_prev);

                fprintf(stderr,"Ends -%d-  %g %g %g %g %g %g %g \n",__LINE__,
                        nH,Pall.ConvertUtoT*U,nHhot,Thot,Mass,Mhot,Mhot_prev);

                Phydro[index]->Uhot = Phydro[index]->Mhot = Phydro[index]->Rhohot = 0.e0;
                Phydro[index]->MultiphaseFlag = false;
            }
        }
    }

    return dudt;
}

static void CalcMultiphaseCoolingASRCH(void){

    ASRCHMakeRedshiftInterpolateData(Pall.Redshift);

    double EnergyLossThisStep = 0.e0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){

            Phydro[i]->DuCooling = MultiphaseCoolingSolverASRCH(i);
            double dt = Phydro[i]->dt_hydro;

#define USE_BOTTOM_TEMPERATURE
#ifdef USE_BOTTOM_TEMPERATURE //{
            //if(Pall.ConvertUtoT*Phydro[i]->U < 10.0){
            if(Pall.ConvertUtoT*(Phydro[i]->U+Phydro[i]->DuCooling*dt) < 10.0){
                //Phydro[i]->DuCooling = Pall.ConvertTtoU*(11.0-Pall.ConvertUtoT*Phydro[i]->U)/dt;
#if 0
                fprintf(stderr,"Oh!Cool %ld %g | %g %g %g\n",
                        PhydroBody(i)->GlobalID,
                        Pall.ConvertUtoT*(Phydro[i]->U+Phydro[i]->DuCooling*dt),
                        Pall.ConvertUtoT*(Phydro[i]->U),
                        Pall.ConvertUtoT*(Phydro[i]->DuCooling*dt),dt);
                fflush(NULL);
#endif
                Phydro[i]->DuCooling = (Pall.ConvertTtoU*10.0-Phydro[i]->U)/dt;
            }
#endif // USE_BOTTOM_TEMPERATURE //}

#define USE_TOP_TEMPERATURE
#ifdef USE_TOP_TEMPERATURE //{
            if(Pall.ConvertUtoT*(Phydro[i]->U+Phydro[i]->DuCooling*dt) > 1.e9){
                fprintf(stderr,"Oh!Hot %ld %g | %g %g\n",
                        PhydroBody(i)->GlobalID,
                        Pall.ConvertUtoT*(Phydro[i]->U+Phydro[i]->DuCooling*dt),
                        Pall.ConvertUtoT*(Phydro[i]->U),
                        Pall.ConvertUtoT*(Phydro[i]->DuCooling*dt));
                fflush(NULL);
                Phydro[i]->DuCooling = (Pall.ConvertTtoU*1.e9-Phydro[i]->U)/dt;
            }
#endif // USE_TOP_TEMPERATURE //}

            Phydro[i]->UPred = Phydro[i]->U += Phydro[i]->DuCooling*dt;
            assert(Phydro[i]->U>0.e0);

#if 0
            if(Phydro[i]->MultiphaseFlag){
                fprintf(stderr,"!! Etot, Ehot = %g %g\n",
                        Phydro[i]->U*Phydro[i]->Mass,Phydro[i]->Uhot*Phydro[i]->Mhot);
            }
#endif

            EnergyLossThisStep -= PhydroMass(i)*Phydro[i]->DuCooling*dt;
        }
    }
    Pall.CoolingEnergyLoss += EnergyLossThisStep;

    return ;
}

#   endif // USE_*_COOLING_FUNCTIONS
#endif // USE_MULTIPHASE_MODEL //}//{


void CalcCooling(void){

#ifdef COOLING_RUN //{
    if(Pall.TCurrent < InitialSkipCoolingRoutine)
        return ;

    double TimingResultThisRoutine = GetElapsedTime();

#ifdef USE_MULTIPHASE_MODEL //{
#   if defined(USE_ASRCH)
    CalcMultiphaseCoolingASRCH();
#   else
#   error This mode cannot work without ASRCH
#   endif // USE_*_COOLING_FUNCTIONS
#else // USE_MULTIPHASE_MODEL //}//{
#   if defined(USE_ASRCH)
    CalcCoolingASRCH();
#   elif defined(USE_CLOUDY_COOLING_FUNCTIONS)
    CalcCoolingCloudy();
#   elif defined(USE_SPAANS2008_COOLING_FUNCTIONS)
    CalcCoolingSpaans2008();
#   elif defined(USE_SPAANS1997_COOLING_FUNCTIONS)
    CalcCoolingSpaans1997();
#   endif // USE_*_COOLING_FUNCTIONS
#endif // USE_MULTIPHASE_MODEL //}

    TimingResults.CoolingThisStep = GetElapsedTime()-TimingResultThisRoutine;
#endif // COOLING_RUN //}

    return;
}

#ifdef USE_CLOUDY_COOLING_FUNCTIONS //{
static inline double ReturnCloudyCoolingRate(const double Rho, const double U, const double Metallicity) __attribute__((always_inline));
static inline double ReturnCloudyCoolingRate(const double Rho, const double U, const double Metallicity){

	double X = 1.e0 - HeliumAbandance - Metallicity;
	double LogT = log10(Pall.ConvertUtoT*U);
	double nH = Pall.ConvertNumberDensityToCGS*Rho*X;
	double LognH = log10(nH);


    if(LogT < CLOUDY_COOLINGHEATING_TEMPERATURE_MAX){
        return -Pall.ConvertCoolingRate*
            SQ(nH)*(ReturnCloudyCoolingHeatingValue(LognH,LogT,Metallicity))/Rho;
    } else {
        double Lambda = ReturnCloudyCoolingHeatingValue(LognH,CLOUDY_COOLINGHEATING_TEMPERATURE_MAX,Metallicity);
        // Extraplation
        double Ti = Pall.ConvertUtoT*U;
        const double Tmax = pow(10.0,CLOUDY_COOLINGHEATING_TEMPERATURE_MAX);
        Lambda *= sqrt(Ti/Tmax);
        return -Pall.ConvertCoolingRate*SQ(nH)*Lambda/Rho;
    }

    // return -Pall.ConvertCoolingRate*
        // SQ(nH)*(ReturnCloudyCoolingHeatingValue(LognH,LogT,Metallicity))/Rho;
}

#ifdef USE_FUVFEEDBACK //{
#define REDUCEDFUVHEATINGRATE_TMIN (9e+3)
#define REDUCEDFUVHEATINGRATE_TMAX (1.e+4)
const static double iBase = 1.0/(REDUCEDFUVHEATINGRATE_TMAX-REDUCEDFUVHEATINGRATE_TMIN);
static double ReducedFUVHeatingRate(const double Heating, const double Ti){

    if(Ti < REDUCEDFUVHEATINGRATE_TMIN){
        return Heating;
    }else if(Ti > REDUCEDFUVHEATINGRATE_TMAX){
        return 0.e0;
    } else {
        return Heating*0.5*(1.0-tanh(2.0*(Ti-REDUCEDFUVHEATINGRATE_TMIN)*iBase-1.0));
    }
}
#else  // USE_FUVFEEDBACK //}//{
static double ReducedFUVHeatingRate(const double Heating, const double Ti){
    return 0.e0;
}
#endif // USE_FUVFEEDBACK //}

static double CoolingSolverCloudy(const int index){

	int iter = 0;
    double dt = Phydro[index]->dt_hydro;
#ifdef SET_SNII_TEMPERATURE_GAS //{
    double R = gsl_rng_uniform(RandomGenerator);
    double USN = Phydro[index]->DQheat/PhydroMass(index);
    double TSN = Pall.ConvertUtoT*USN;
    double Heating = 0.0;
    if(R < TSN/SNII_TEMPERATURE){
	    Heating = Pall.ConvertTtoU*SNII_TEMPERATURE/dt;
    } else {
        Heating = 0;
    }
#else //}//{
	double Heating = Phydro[index]->DQheat/(PhydroMass(index)*dt); // [erg/g/s]
#endif //}

    double U = Phydro[index]->U;
    double Uold = U, Uinit = U, Uupper = U, Ulower = U;
	double Metallicity = Phydro[index]->Z;
    double Rho = Phydro[index]->Rho;

#ifdef USE_FUVFEEDBACK //{

#define FUV_Epsilon (0.05)
#define FUV_Zsolar (0.0134)
#   ifdef USE_FUVFEEDBACK_CONSTANT_VALUE //{
    double nH = Pall.ConvertNumberDensityToCGS*Phydro[index]->Rho;
    double Gamma_Heating = 1.e-24*FUV_Epsilon*FUVFEEDBACK_CONSTANT_VALUE*(Metallicity/FUV_Zsolar); // [erg/s]
    double FUVHeating = FUVConversionFactor*nH*Gamma_Heating; // [erg/g/s]
    // fprintf(stderr,"C! %g %g %g\n",FUVFEEDBACK_CONSTANT_VALUE,Gamma_Heating,FUVHeating);
    // fprintf(stderr,"C! %g\n",FUVConversionFactor);
#   else // USE_FUVFEEDBACK_CONSTANT_VALUE //}//{
    double nH = Pall.ConvertNumberDensityToCGS*Phydro[index]->Rho;

#       ifdef USE_FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE //{
    double Ncol = nH*FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE;
#       else // USE_FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE //}//{
    double Ncol = nH*(Phydro[index]->Rho/fabs(2*Phydro[index]->GradRho))*Pall.UnitLength;
#       endif // USE_FUVFEEDBACK_COLUMNDENSITY_CONSTANT_SCALE //}

    double Zdep = (Metallicity/FUV_Zsolar);
#       ifdef USE_FUVFEEDBACK_LOCAL_CORRECTION //{
#           ifdef USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //{
    double G = (Phydro[index]->G0thin-Phydro[index]->G0thinLocal)*exp(-Zdep*FUVFEEDBACK_SIGMA1000*Ncol); //
#           else // USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //}//{
    double G = (Phydro[index]->G0thin-Phydro[index]->G0thinLocal)*exp(-FUVFEEDBACK_SIGMA1000*Ncol); //
#           endif // USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //}
    G += Phydro[index]->G0extLocal;
#       else // USE_FUVFEEDBACK_LOCAL_CORRECTION //}//{
#           ifdef USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //{
    double G = Phydro[index]->G0thin*exp(-Zdep*FUVFEEDBACK_SIGMA1000*Ncol); //
#           else // USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //}//{
    double G = Phydro[index]->G0thin*exp(-FUVFEEDBACK_SIGMA1000*Ncol); //
#           endif // USE_FUVFEEDBACK_METAL_DEPEND_SIGMA1000  //}
#       endif // USE_FUVFEEDBACK_LOCAL_CORRECTION //}

    double Gamma_Heating = 1.e-24*FUV_Epsilon*G*Zdep; // [erg/s]
    double FUVHeating = FUVConversionFactor*nH*Gamma_Heating/Rho; // [erg/g/s]
#   endif // USE_FUVFEEDBACK_CONSTANT //}

#else // USE_FUVFEEDBACK //}//{

    double Gamma_Heating = 0.e0; // [erg/s]
    double FUVHeating = 0.e0; 
#endif // USE_FUVFEEDBACK //}

#if (SingleStepCoolingSkip) //{
    if(Phydro[index]->DQheat > 0.e0){
        Phydro[index]->DQheat = 0.e0;
        return Heating;
    } else {
        Heating = 0.0;
    }
#endif // SingleStepCoolingSkip //}

#if 1
#ifdef USE_EXACT_COOLING_SOLVER //{
    double dudt_exact = 0.0;
    {   
        double X = 1.e0 - HeliumAbandance - Metallicity;
        double nH = Pall.ConvertNumberDensityToCGS*Rho*X;
        double LognH = log10(nH);
        double LogT = log10(Pall.ConvertUtoT*U);
	    // Heating = (Phydro[index]->DQheat/dt)*(Pall.UnitMass*SQ(Pall.UnitLength)/CUBE(Pall.UnitTime)); // [erg/s]
	    Phydro[index]->DQheat = 0.e0;
        const double factor = CUBE(Pall.UnitTime)/SQ(Pall.UnitLength); // [erg/g/s], [cm^2/s^3]
        //return factor*ReturnCloudyCoolingHeatingValueExactSolver(LognH,LogT,
                 //Metallicity,dt*Pall.UnitTime);
        // return factor*ReturnCloudyCoolingHeatingValueExactSolverWithFUVHeating(LognH,LogT,
                // Metallicity,dt*Pall.UnitTime,(Heating+Gamma_Heating)/nH);
        return factor*ReturnCloudyCoolingHeatingValueExactSolverWithFUVHeating(LognH,LogT,
                Metallicity,dt*Pall.UnitTime,Gamma_Heating/nH);

        // dudt_exact = factor*ReturnCloudyCoolingHeatingValueExactSolverWithFUVHeating(LognH,LogT,
                // Metallicity,dt*Pall.UnitTime,(Heating+Gamma_Heating)/nH);
    }
#endif //USE_EXACT_COOLING_SOLVER //}
#endif


    double LambdaNet = ReturnCloudyCoolingRate(Rho,U,Metallicity); // [erg/g/s]
	if(U < Uold+(LambdaNet+Heating+ReducedFUVHeatingRate(FUVHeating,Pall.ConvertUtoT*U))*dt){
		Uupper *= sqrt11;
		while(Uupper < Uold+(ReturnCloudyCoolingRate(Rho,Uupper,Metallicity)+Heating+ReducedFUVHeatingRate(FUVHeating,Pall.ConvertUtoT*Uupper))*dt){
			Uupper *= sqrt11;
            if(Pall.ConvertUtoT*Uupper >= 1.e+9){
                break;
            }
		}
	}

    // int counter = 0;
    // double Ull = Pall.ConvertTtoU*2.7;
    // double U10k = Pall.ConvertTtoU*10;
    if(U > Uold+(LambdaNet+Heating+ReducedFUVHeatingRate(FUVHeating,Pall.ConvertUtoT*U))*dt){ 
		Ulower *= isqrt11;
		while(Ulower > Uold+(ReturnCloudyCoolingRate(Rho,Ulower,Metallicity)+Heating+ReducedFUVHeatingRate(FUVHeating,Pall.ConvertUtoT*Uupper))*dt){
			Ulower *= isqrt11;

            if(Pall.ConvertUtoT*Ulower <= 1.e+1){
                break;
            }
		}
	}


    double DeltaU;
	do{
		U = 0.5*(Uupper + Ulower);
        LambdaNet = ReturnCloudyCoolingRate(Rho,U,Metallicity);

		//if( (U-(Uold+(LambdaNet+Heating))*dt) > 0.e0){
		if( (U-(Uold+(LambdaNet+Heating+ReducedFUVHeatingRate(FUVHeating,Pall.ConvertUtoT*U))*dt)) > 0.e0){
			Uupper = U;
        } else {
			Ulower = U;
        }

		DeltaU = Uupper-Ulower;
		iter ++;

		if(iter >= (MaxIteration-10))
			fprintlmpi(U);
	} while( (fabs(DeltaU/U) > 1.e-6) && (iter < MaxIteration) );


	if( iter >= MaxIteration ){
		fprintf(stderr,"Failed to converge : function %s line %d file %s\n",
                __FUNCTION__,__LINE__,__FILE__);
		fprintf(stderr,"U = %g\n Rho = %g\n dt = %g\n",Phydro[index]->U,Phydro[index]->Rho,dt);
        fflush(NULL);
        MPI_Abort(MPI_COMM_WORLD,-1);
		exit(CoolingConvergenceFailure);
	}

	Phydro[index]->DQheat = 0.e0;

	return ( (U-Uinit)/dt );
}

static void CalcCoolingCloudy(void){

    CloudyDataInterpolationRedshift(Pall.Redshift);

    int counter = 0;

    double EnergyLossThisStep = 0.e0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            Phydro[i]->DuCooling = CoolingSolverCloudy(i);
            double dt = Phydro[i]->dt_hydro;

#define USE_BOTTOM_TEMPERATURE
#ifdef USE_BOTTOM_TEMPERATURE //{
            if(Pall.ConvertUtoT*Phydro[i]->U < 10.0){
                Phydro[i]->DuCooling = Pall.ConvertTtoU*(11.0-Pall.ConvertUtoT*Phydro[i]->U)/dt;
            }
#endif // USE_BOTTOM_TEMPERATURE //}

#if 0
#define USE_TOP_TEMPERATURE
#ifdef USE_TOP_TEMPERATURE //{
            if(Pall.ConvertUtoT*Phydro[i]->U > 1.e9){
                Phydro[i]->DuCooling = Pall.ConvertTtoU*(1.e9-Pall.ConvertUtoT*Phydro[i]->U)/dt;
            }
#endif // USE_TOP_TEMPERATURE //}
#else
#define USE_TOP_TEMPERATURE
#ifdef USE_TOP_TEMPERATURE //{
            if(Phydro[i]->U+Phydro[i]->DuCooling*dt > Pall.ConvertTtoU*1.e9){
                Phydro[i]->DuCooling = (Pall.ConvertTtoU*1.e9-Phydro[i]->U)/dt;
            }
#endif // USE_TOP_TEMPERATURE //}
#endif

            Phydro[i]->UPred = Phydro[i]->U += Phydro[i]->DuCooling*dt;
            assert(Phydro[i]->U>0.e0);

            EnergyLossThisStep -= PhydroMass(i)*Phydro[i]->DuCooling*dt;
        }
    }
    Pall.CoolingEnergyLoss += EnergyLossThisStep;

    return ;
}
#endif // USE_CLOUDY_COOLING_FUNCTIONS //}


#if defined(USE_SPAANS1997_COOLING_FUNCTIONS)
#define CFmin   (1.e1)  // 10^0.6
#define CFmax   (1.e+8) // 10^8
#define NCOOL	20000
#define NCOOLINGCURVE	2
static double Lmd0[NCOOLINGCURVE][NCOOL];
static double CFStep;

enum {
    MetalZero,
    MetalSolar,
    CFNumber,
};

static void InitializeCoolingTableSpaans1997(void){

    logTmin = log10(CFmin);
    logTmax = log10(CFmax);

    MakeCoolingTable(MetalZero);
    MakeCoolingTable(MetalSolar);

    double metal[4] = {0.0,1.e-4,0.02,0.05};
    char fname[MaxCharactersInLine];

    for(int i=0;i<4;i++){
        FILE *fp;
        sprintf(fname,"Cooling.%02d",i);
        FileOpen(fp,fname,"w");
        for(double T=1.e1;T<1.e+8;T*=1.1){
            fprintf(fp,"%e %e\n",log10(T),
                log10(MetalInterpolatedCoolingRate(T,metal[i])));
        }
        fclose(fp);
    }

    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"Conversion Factor of Cooling Rate = %g\n",Pall.ConvertCoolingRate);

    return;

}

/*
 * If the temperature reaches to the upper limiter of the cooling tables, the
 * cooling rate should be proportional to T^0.5.
 */
static void MakeCoolingTable(const int mode){

	int iw;
	double w,t,step;
    //int SizeTable = sizeof(CoolingTable_s0)/(sizeof(double));
	FILE *fp;
    char fname[MaxCharactersInLine];

    double *Source0;
    double *Source1;
    double *Source2;
    double *Source3;
    double *Dest; 

    if(mode == MetalZero){
        Source0 = CoolingTable_s0;
        Source1 = CoolingTable_s1;
        Source2 = CoolingTable_s2;
        Source3 = CoolingTable_s3;
        Dest = Lmd0[0];
        Snprintf(fname,"%s/CoolingTable.z0",PathCoolingHeatingLog);
	    FileOpen(fp,fname,"w");
    }else if (mode == MetalSolar){
        Source0 = CoolingTable_v0;
        Source1 = CoolingTable_v1;
        Source2 = CoolingTable_v2;
        Source3 = CoolingTable_v3;
        Dest = Lmd0[1];
        Snprintf(fname,"%s/CoolingTable.z002",PathCoolingHeatingLog);
	    FileOpen(fp,fname,"w");
    } else {
        MPI_Finalize();
        fprintf(stderr,"Cooling Table Generation Failure : function %s line %d file %s\n",
                __FUNCTION__,__LINE__,__FILE__);
        exit(CoolingTableGenerationFailure);
    }

	t = 8.e0;
	iw = (int)(5.e0*t - 2.e0);

	step = (logTmax-logTmin)/NCOOL;
	CFStep = 1.e0/step; // NCOOL/(log10(CFmax)-log10(CFmin))
                        // When you calculate (log10(T)-log10(CFmin)) * CFstep,
                        // you can obtain the index "i" which is the index 
                        // in order to access the cooling curve by Lmd[][i]

	for(int i=0;i<NCOOL;i++){
		t = 1.e0 + step*i;
		if(t<8.0){
			iw = (int)(5.e0*t - 2.e0);
			w = t - (iw + 2.e0)/5.e0;
			Dest[i] = Source0[iw]+Source1[iw]*w+Source2[iw]*SQ(w)+Source3[iw]*CUBE(w);
		} else {
			iw = (int)(5.e0*8.e0 - 2.e0);
			w = 8.e0 - (iw + 2.e0)/5.e0;
            Dest[i] = Source0[iw]+Source1[iw]*w+Source2[iw]*SQ(w)+Source3[iw]*CUBE(w);
		}
		fprintf(fp,"%e %e\n",pow(10.0,t),Dest[i]);
	}
	fclose(fp);


    step = (8.0-1.0)/NCOOL;
    if(mode == MetalZero){
        Snprintf(fname,"%s/CoolingTable.z0",PathCoolingHeatingLog);
	    FileOpen(fp,fname,"w");
        for(int i=0;i<NCOOL;i++){
            t = 1.e0 + step*i;
            fprintf(fp,"%e %e\n",pow(10.0,t),
                    log10(CoolingFunction(pow(10.0,t),mode))-32.688);
        }
    }else if (mode == MetalSolar){
        Snprintf(fname,"%s/CoolingTable.z002",PathCoolingHeatingLog);
	    FileOpen(fp,fname,"w");
        for(int i=0;i<NCOOL;i++){
            t = 1.e0 + step*i;
            fprintf(fp,"%e %e\n",pow(10.0,t),
                    log10(CoolingFunction(pow(10.0,t),mode))-32.688);
        }
    }
	fclose(fp);

	return;
}

static double CoolingFunction(double T, const int mode){

    int TOP,BOTTOM;
    double logT,L,CENTER;
    logT = log10(T);
#define InvDumpingFactor   (300)
	
#ifdef MOLECULAR_COOLING

#ifdef COOLING_CUTOFF_TEMPERATURE
	if(logT<T_cut){
		L = 0.e0;
#else
	if(logT<1.e0){
		L = 0.e0;
#endif

#else
	if(logT<4.e0){ // exp
        if(logT<4.0){
            L = 0.e0;
        }else{
            L = Lmd0[mode][(int)((4.e0-logTmin)*CFStep)];
            L = pow(10.e0,L) * exp((logT-4.e0)*InvDumpingFactor);
        }
#endif
	} else if(logT<CFmax){
		CENTER = (logT-1.e0)*CFStep;
		BOTTOM = MAX((int)CENTER,0);
		TOP = MIN(BOTTOM+1,NCOOL-1);
		BOTTOM = TOP-1;
		L = ((TOP-CENTER)*Lmd0[mode][BOTTOM]+(CENTER-BOTTOM)*Lmd0[mode][TOP]);
		L = pow(10.e0,L);
	} else {
		L = pow(10.e0,L);
	}
	return (Pall.ConvertCoolingRate*L);
}

double MetalInterpolatedCoolingRate(const double T, double metal){

    double Lrate,logmetal;
    static const double ZeroMetal = 0.02*1.e-5;
    static const double log000 = -6.6990; //log10(0.02*1.e-5)
    static const double log002 = -1.6990; //log10(0.02)
    static const double idSZ = 1.e0/(-1.6990-(-6.6990)); //1/dSZ

    if(metal<ZeroMetal){
        Lrate = CoolingFunction(T,MetalZero);
    } else if(metal<0.02){
        logmetal = log10(metal);
        Lrate = (log002-logmetal)*CoolingFunction(T,MetalZero)+
                            (logmetal-log000)*CoolingFunction(T,MetalSolar);
        Lrate *= idSZ;
    }else{
        Lrate = CoolingFunction(T,MetalSolar);
    }

    return Lrate;
}

static double MetalDependingCoolingRate(const int index, const double ui, const double rhoi, const double g0i, const double Metallicity){

	double X = 1.e0 - HeliumAbandance - Metallicity;
	double Ti = Pall.ConvertUtoT*ui;
	double nH = Pall.ConvertNumberDensityToCGS*rhoi*X;

	double LRad = SQ(nH)*(MetalInterpolatedCoolingRate(Ti,Metallicity));
#ifdef USE_INVERSE_COMPTON_COOLING //{
	double ne = nH * ( 1.e0-0.5*HeliumAbandance )/X;
	double LCc = InverseComptonCooling(Ti,ne,Pall.Redshift);
#else // USE_INVERSE_COMPTON_COOLING
	double LCc = 0.e0;
#endif // USE_INVERSE_COMPTON_COOLING //}

#ifdef USE_FARULTRAVIOLET_HEATING //{
	double HeatingFUV = FUV(Ti,rhoi,g0i,Metallicity);
#else // USE_FARULTRAVIOLET_HEATING
	double HeatingFUV = 0.e0; 
#endif // USE_FARULTRAVIOLET_HEATING //}

    if(Ti < 10.0){
        LRad = LCc = 0.0;
    }
	
	double DuRad = ( (HeatingFUV-LRad-LCc)/rhoi );
	return ( DuRad );
}


// This routine returns only du_cooling for the particle of the ID = index.
static double CoolingSolver(const int index){

	int iter = 0;

    double dt = Phydro[index]->dt_hydro;

    double U = Phydro[index]->U;
    double Uold,Uinit,Uupper,Ulower;
    Uold = Uinit = Uupper = Ulower = U;
	double Metallicity = Phydro[index]->Z;
	double Heating = Phydro[index]->DQheat/(PhydroMass(index)*dt);

#if (SingleStepCoolingSkip)
    if(Phydro[index]->DQheat > 0.e0){
        Phydro[index]->DQheat = 0.e0;
        return Heating;
    }
#endif // SingleStepCoolingSkip

    double Rho = Phydro[index]->Rho;
    //double G0 = Phydro[index]->G0;
    double G0 = 1.0;
    double LambdaNet = MetalDependingCoolingRate(index,U,Rho,G0,Metallicity);
	if(U < Uold+(LambdaNet+Heating)*dt){
		Uupper *= sqrt11;
		while(Uupper < 
                Uold+(MetalDependingCoolingRate(index,Uupper,Rho,G0,Metallicity)+Heating)*dt){
			Uupper *= sqrt11;
		}
	}
    if(U > Uold+(LambdaNet+Heating)*dt){ 
		Ulower *= isqrt11;
		while(Ulower > 
                Uold+(MetalDependingCoolingRate(index,Ulower,Rho,G0,Metallicity)+Heating)*dt){
			Ulower *= isqrt11;
		}
	}

    double DeltaU;
	do{
		U = 0.5*(Uupper + Ulower);
        LambdaNet = MetalDependingCoolingRate(index,U,Rho,G0,Metallicity);

		//if( (U-(Uold+(LambdaNet+Heating))*dt) > 0.e0)
		if( (U-(Uold+(LambdaNet+Heating)*dt)) > 0.e0){
			Uupper = U;
        } else {
			Ulower = U;
        }

		DeltaU = Uupper-Ulower;

		iter ++;

		if(iter >= (MaxIteration-10))
			fprintlmpi(U);
	} while( (fabs(DeltaU/U) > 1.e-6) && (iter < MaxIteration) );

	if( iter >= MaxIteration ){
		fprintf(stderr,"Failed to converge : function %s line %d file %s\n",
                __FUNCTION__,__LINE__,__FILE__);
		fprintf(stderr,"U = %g\n Rho = %g\n dt = %g\n",
			Phydro[index]->U,Phydro[index]->Rho,dt);
        fflush(NULL);
		exit(CoolingConvergenceFailure);
	}

	Phydro[index]->DQheat = 0.e0;

    if(U < 0.e0){
        fprintf(stderr,"In %s, line %d, function %s",__FILE__,__LINE__,__FUNCTION__);
        fprintf(stderr,"U Uinit = %g %g\n",U,Uinit);
        MPI_Finalize();
        exit(CoolingNotFP_NORMAL);
    }

	return ( (U-Uinit)/dt );
}


static void CalcCoolingSpaans1997(void){

    double EnergyLossThisStep = 0.e0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){

            Phydro[i]->DuCooling = CoolingSolver(i);
            Phydro[i]->UPred = Phydro[i]->U += Phydro[i]->DuCooling*Phydro[i]->dt_hydro;
            EnergyLossThisStep -= PhydroMass(i)*Phydro[i]->DuCooling*Phydro[i]->dt_hydro;
        }
    }
    Pall.CoolingEnergyLoss += EnergyLossThisStep;

    return ;
}
#endif // USE_SPAANS1997_COOLING_FUNCTIONS



#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////
//#define SIZE_TABLE_FH2 11
//#define SIZE_TABLE_G0 7
//#define SIZE_TABLE_Z 4
//#define SIZE_TABLE_T 104
//static double TableFH2[SIZE_TABLE_FH2];
//static double TableT[SIZE_TABLE_T];
//static double TableZ[SIZE_TABLE_Z];
//static double TableG0[SIZE_TABLE_G0];
//static double TableLambda[SIZE_TABLE_T][SIZE_TABLE_Z][SIZE_TABLE_G0][SIZE_TABLE_FH2];

static double Zsolar;   // solar metallicity

static void Linear1D(const int n, double x[n], double y[n], double xi, double *yi);
static void Spline1D(const int n, double x[n], double y[n], double xi, double *yi);
static void Spline2D(const int n1, const int n2, 
              double x1[n1], double x2[n2], double y[n1][n2], 
              double x1i, double x2i, double *yi);
static void Spline3D(const int n1, const int n2, const int n3, 
              double x1[n1], double x2[n2], double x3[n3], double y[n1][n2][n3], 
              double x1i, double x2i, double x3i, double *yi);
static void Spline4D(const int n1, const int n2, const int n3, const int n4, 
              double x1[n1], double x2[n2], double x3[n3], double x4[n4], double y[n1][n2][n3][n4], 
              double x1i, double x2i, double x3i, double x4i, double *yi);
static double CoolingSolverSpaans2008(const int index);
static double CoolingRateSpaans2008(const int index, double ui, double rhoi, double Metallicity, double fH2i, double g0i);

static int FlagCoolingTable = OFF;

static void InitializeCoolingTableSpaans2008(void){

    logTmin = log10(1.e+1);
    logTmax = log10(1.e+8);
    Zsolar = 0.02;

    InitializeCoolingRateLookUpTable();

    if(MPIGetMyID() == MPI_ROOT_RANK){
        //double metal[4] = {0.0,0.001,0.02,0.05};
        double metal[] = {0.001*Zsolar,0.01*Zsolar,0.1*Zsolar,0.5*Zsolar,1.0*Zsolar};
        double g0[] = {0.01,1.0,10.0,100.0,1000.0};
        double fH2[] = {0.001,0.01,0.1,0.5,1.0};
        char fname[MaxCharactersInLine];

        for(int i=0;i<5;i++){
            FILE *fp;
            // metallicity dependence
            sprintf(fname,"%s/CoolingSpaans2008.Z%02d",PathCoolingHeatingLog,i);
            FileOpen(fp,fname,"w");
            fprintf(fp,"#(Z,H2,G0)=(%f,%f,%f)\n",metal[i],fH2[2],g0[1]);
            for(double T=1.e1;T<1.e+8;T*=1.1){
                fprintf(fp,"%e %e\n",log10(T),
                    log10(InterpolatedCoolingRateSpaans2008(T,metal[i],g0[1],fH2[2])/Pall.ConvertCoolingRate));
            }
            fclose(fp);
            // FUV dependence
            sprintf(fname,"%s/CoolingSpaans2008.FUV%02d",PathCoolingHeatingLog,i);
            FileOpen(fp,fname,"w");
            fprintf(fp,"#(Z,H2,G0)=(%f,%f,%f)\n",metal[3],fH2[2],g0[i]);
            for(double T=1.e1;T<1.e+8;T*=1.1){
                fprintf(fp,"%e %e\n",log10(T),
                    log10(InterpolatedCoolingRateSpaans2008(T,metal[3],g0[i],fH2[2])/Pall.ConvertCoolingRate));
            }
            fclose(fp);
            // H2 fraction dependence
            sprintf(fname,"%s/CoolingSpaans2008.H2%02d",PathCoolingHeatingLog,i);
            FileOpen(fp,fname,"w");
            fprintf(fp,"#(Z,H2,G0)=(%f,%f,%f)\n",metal[3],fH2[i],g0[1]);
            for(double T=1.e1;T<1.e+8;T*=1.1){
                fprintf(fp,"%e %e\n",log10(T),
                    log10(InterpolatedCoolingRateSpaans2008(T,metal[3],g0[1],fH2[i])/Pall.ConvertCoolingRate));
            }
            fclose(fp);
        }
        fprintf(stderr,"Conversion Factor of Cooling Rate = %g\n",Pall.ConvertCoolingRate);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    return;
}

//double TableLogLambda[SIZE_TABLE_Z][SIZE_TABLE_G0][SIZE_TABLE_FH2][SIZE_TABLE_T];
static gsl_interp_accel *AccAccessLogT[SIZE_TABLE_Z][SIZE_TABLE_G0][SIZE_TABLE_FH2];
static gsl_spline       *SplineAccessLogT[SIZE_TABLE_Z][SIZE_TABLE_G0][SIZE_TABLE_FH2];
void InitializeCoolingRateLookUpTable(void){

    if(FlagCoolingTable == ON)  return;

    for(int im=0;im<SIZE_TABLE_Z  ;im++) TableLogZ  [im] = log10(TableLogZ  [im]);
    for(int ig=0;ig<SIZE_TABLE_G0 ;ig++) TableLogG0 [ig] = log10(TableLogG0 [ig]);
    for(int ih=0;ih<SIZE_TABLE_FH2;ih++) TableLogFH2[ih] = log10(TableLogFH2[ih]); 

    //InitializeLifeTimeLookUpTable();
    for(int im=0;im<SIZE_TABLE_Z;im++){
        for(int ig=0;ig<SIZE_TABLE_G0;ig++){
            for(int ih=0;ih<SIZE_TABLE_FH2;ih++){
                AccAccessLogT[im][ig][ih] = gsl_interp_accel_alloc();
                //SplineAccessLogT[im][ig][ih] = gsl_spline_alloc(gsl_interp_linear,SIZE_TABLE_T);
                SplineAccessLogT[im][ig][ih] = gsl_spline_alloc(gsl_interp_cspline,SIZE_TABLE_T);
                gsl_spline_init(SplineAccessLogT[im][ig][ih], TableLogT, TableLogLambda[im][ig][ih], SIZE_TABLE_T);
            }
        }
    }
    
    FlagCoolingTable = ON;

    return;
}


void CalcCoolingSpaans2008(void){
    
    double EnergyLossThisStep = 0.e0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            Phydro[i]->DuCooling = CoolingSolverSpaans2008(i);
            Phydro[i]->UPred = Phydro[i]->U += Phydro[i]->DuCooling*Phydro[i]->dt_hydro;

            EnergyLossThisStep -= PhydroMass(i)*Phydro[i]->DuCooling*Phydro[i]->dt_hydro;
#ifdef TASK_GALACTIC_CENTER
            Phydro[i]->dZII += PhydroMass(i)*Phydro[i]->DuCooling*Phydro[i]->dt_hydro;
            if(Phydro[i]->DuCooling < 0.0)
                Phydro[i]->dZIa += PhydroMass(i)*Phydro[i]->DuCooling*Phydro[i]->dt_hydro;
#if 1
            double HeatingFUV = FUV(Pall.ConvertUtoT*Phydro[i]->UPred,Phydro[i]->Rho,
                    Phydro[i]->G0,Phydro[i]->Z)/Phydro[i]->Rho;
            Phydro[i]->ZIa += PhydroMass(i)*HeatingFUV*Phydro[i]->dt_hydro;
            Phydro[i]->ZII += PhydroMass(i)*(Phydro[i]->DuCooling-HeatingFUV)*Phydro[i]->dt_hydro;
#endif 
#endif //TASK_GALACTIC_CENTER
        }
    }
    Pall.CoolingEnergyLoss += EnergyLossThisStep;

    return;
}

/*
 * This routine returns only du_cooling for the particle of the ID = index.
 */
static double CoolingSolverSpaans2008(const int index){

	int iter = 0;

    double dt = Phydro[index]->dt_hydro;
    double U = Phydro[index]->U;
    double Uold = U, Uinit = U, Uupper = U, Ulower = U;
    double Rho = Phydro[index]->Rho;
	double Metallicity = Phydro[index]->Z;
	double Heating = Phydro[index]->DQheat/(PhydroMass(index)*dt);

    double G0 = Phydro[index]->G0;
    double fH2 = Phydro[index]->fH2;

#if (SingleStepCoolingSkip)
    if(Phydro[index]->DQheat > 0.e0){
        Phydro[index]->DQheat = 0.e0;
        return Heating;
    }
#endif // SingleStepCoolingSkip

    double LambdaNet = CoolingRateSpaans2008(index,U,Rho,Metallicity,fH2,G0);
	if(U < Uold+(LambdaNet+Heating)*dt){
		Uupper *= sqrt11;
		while(Uupper < 
                Uold+(CoolingRateSpaans2008(index,Uupper,Rho,Metallicity,fH2,G0)+Heating)*dt){
			Uupper *= sqrt11;
		}
	}
    if(U > Uold+(LambdaNet+Heating)*dt){ 
		Ulower *= isqrt11;
		while(Ulower > 
                Uold+(CoolingRateSpaans2008(index,Ulower,Rho,Metallicity,fH2,G0)+Heating)*dt){
			Ulower *= isqrt11;
		}
	}

    double DeltaU;
	do{
		U = 0.5*(Uupper + Ulower);
        LambdaNet = CoolingRateSpaans2008(index,U,Rho,Metallicity,fH2,G0);
		//if( (U-(Uold+(LambdaNet+Heating))*dt) > 0.e0)
		if( (U-(Uold+(LambdaNet+Heating)*dt)) > 0.e0)
			Uupper = U;
		else
			Ulower = U;
		DeltaU = Uupper-Ulower;
		iter ++;
		if(iter >= (MaxIteration-10)){
			double T = U*Pall.ConvertUtoT;
			eprintlmpi(T);
			eprintlmpi(fH2);
			fprintlmpi(G0);
			fprintlmpi(Metallicity);
			eprintlmpi(LambdaNet);
        }
	} while( (fabs(DeltaU/U) > 1.e-6) && (iter < MaxIteration) );

	if( iter >= MaxIteration ){
		fprintf(stderr,"Failed to converge : function %s line %d file %s\n",
                __FUNCTION__,__LINE__,__FILE__);
		fprintf(stderr," U = %g\n T = %g\n Rho = %g\n dt = %g\n",
			Phydro[index]->U,Phydro[index]->U*Pall.ConvertUtoT,Phydro[index]->Rho,dt);
        fflush(NULL);
		exit(CoolingConvergenceFailure);
	}

	Phydro[index]->DQheat = 0.e0;

    if(U < 0.e0){
        fprintf(stderr,"In %s, line %d, function %s",__FILE__,__LINE__,__FUNCTION__);
        fprintf(stderr,"U Uinit = %g %g\n",U,Uinit);
        MPI_Finalize();
        exit(CoolingNotFP_NORMAL);
    }
	return ( (U-Uinit)/dt );
}

static double CoolingRateSpaans2008(const int index, double ui, double rhoi, double Metallicity, double fH2i, double g0i){
	double X = 1.e0 - HeliumAbandance - Metallicity;
	double Ti = Pall.ConvertUtoT*ui;
	double nH = Pall.ConvertNumberDensityToCGS*rhoi*X;
    double LRad;// = SQ(nH)*InterpolatedCoolingRateSpaans2008(Ti,Metallicity,g0i,fH2i);
#define T_UPPER_LIMIT_FOR_SPAANS2008 (1.e+8)
    if(Ti < T_UPPER_LIMIT_FOR_SPAANS2008){ // 1.e+8 is the limit of the cooling function by Spaans (2008).
        LRad = SQ(nH)*InterpolatedCoolingRateSpaans2008(Ti,Metallicity,g0i,fH2i);
    } else {
        double Lambda = InterpolatedCoolingRateSpaans2008(T_UPPER_LIMIT_FOR_SPAANS2008,Metallicity,g0i,fH2i);
        // Extraplation
        Lambda *= sqrt(Ti/T_UPPER_LIMIT_FOR_SPAANS2008);
        LRad = SQ(nH)*Lambda;
    }

#ifdef USE_INVERSE_COMPTON_COOLING //{
	double ne = nH * ( 1.e0-0.5*HeliumAbandance )/X;
	double LCc = InverseComptonCooling(Ti,ne,Pall.Redshift);
#else // USE_INVERSE_COMPTON_COOLING
	double LCc = 0.e0;
#endif //USE_INVERSE_COMPTON_COOLING //}

#ifdef USE_FARULTRAVIOLET_HEATING //{
	double HeatingFUV = FUV(Ti,rhoi,g0i,Metallicity);
#else  // USE_FARULTRAVIOLET_HEATING
	double HeatingFUV = 0.e0;
#endif // USE_FARULTRAVIOLET_HEATING //}
	
	double DuRad = ( (HeatingFUV-LRad-LCc)/rhoi );
	return ( DuRad );
}

double InterpolatedCoolingRateSpaans2008(double Tg, double Zg, double G0g, double fH2g){
    //double GetLifeTimeFromStellarMass(const double Mass, const double Metal);

#ifdef COOLING_CUTOFF_TEMPERATURE
	if(log10(Tg)<T_cut)         return 0.e0;
#else  // COOLING_CUTOFF_TEMPERATURE
	if(log10(Tg)<TableLogT[0])  return 0.e0;
#endif // COOLING_CUTOFF_TEMPERATURE
     
    Zg = log10(MAX(Zg/Zsolar,pow(10,TableLogZ[0])));
    G0g = log10(MAX(G0g,pow(10,TableLogG0[0])));
    fH2g = log10(MAX(fH2g,pow(10,TableLogFH2[0])));

    // linear interpolation
    int im = 0, ig = 0, ih = 0;
    while(TableLogZ  [im] < Zg)        im++;
    while(TableLogG0 [ig] < G0g)       ig++;
    while(TableLogFH2[ih] < fH2g)      ih++;
    int ipm, ipg, iph;
    int inm = 0;
    if(Zg>=TableLogZ[SIZE_TABLE_Z-1])           inm = SIZE_TABLE_Z-1;
    int ing = 0;
    if(G0g>=TableLogG0[SIZE_TABLE_G0-1])        ing = SIZE_TABLE_G0-1;
    int inh = 0;
    if(fH2g>=TableLogFH2[SIZE_TABLE_FH2-1])     inh = SIZE_TABLE_FH2-1;
    double l = 0.e0, lintp = 0.e0;
    double dZinv = 1.0/fabs(TableLogZ  [im]-TableLogZ  [abs(im-1)]);
    double dGinv = 1.0/fabs(TableLogG0 [ig]-TableLogG0 [abs(ig-1)]);
    double dHinv = 1.0/fabs(TableLogFH2[ih]-TableLogFH2[abs(ih-1)]);
    if( (Zg>TableLogZ[0] && Zg<TableLogZ[SIZE_TABLE_Z-1])
            && (G0g>TableLogG0[0] &&  G0g<TableLogG0[SIZE_TABLE_G0-1]) 
            && (fH2g>TableLogFH2[0] && fH2g<TableLogFH2[SIZE_TABLE_FH2-1])
            ){
        for(int imm=0; imm<2; imm++){
            for(int igg=0; igg<2; igg++){
                for(int ihh=0; ihh<2; ihh++){
                    inm = im-1+imm;
                    ing = ig-1+igg;
                    inh = ih-1+ihh;
                    ipm = im - imm;
                    ipg = ig - igg;
                    iph = ih - ihh;
                    l = gsl_spline_eval(SplineAccessLogT[inm][ing][inh],log10(Tg),AccAccessLogT[inm][ing][inh]);
                    lintp += l*fabs((Zg-TableLogZ[ipm])*(G0g-TableLogG0[ipg])*(fH2g-TableLogFH2[iph]));
                }
            }
        }
        return Pall.ConvertCoolingRate*pow(10.0,lintp*dZinv*dGinv*dHinv);
    }
    else if( !(Zg>TableLogZ[0] && Zg<TableLogZ[SIZE_TABLE_Z-1])
            && (G0g>TableLogG0[0] &&  G0g<TableLogG0[SIZE_TABLE_G0-1]) 
            && (fH2g>TableLogFH2[0] && fH2g<TableLogFH2[SIZE_TABLE_FH2-1])
            ){
        for(int igg=0; igg<2; igg++){
            for(int ihh=0; ihh<2; ihh++){
                ing = ig-1+igg;
                inh = ih-1+ihh;
                ipg = ig - igg;
                iph = ih - ihh;
                l = gsl_spline_eval(SplineAccessLogT[inm][ing][inh],log10(Tg),AccAccessLogT[inm][ing][inh]);
                lintp += l*fabs((G0g-TableLogG0[ipg])*(fH2g-TableLogFH2[iph]));
            }
        }
        return Pall.ConvertCoolingRate*pow(10.0,lintp*dGinv*dHinv);
    }
    else if( (Zg>TableLogZ[0] && Zg<TableLogZ[SIZE_TABLE_Z-1])
            && !(G0g>TableLogG0[0] &&  G0g<TableLogG0[SIZE_TABLE_G0-1]) 
            && (fH2g>TableLogFH2[0] && fH2g<TableLogFH2[SIZE_TABLE_FH2-1])
            ){
        for(int imm=0; imm<2; imm++){
            for(int ihh=0; ihh<2; ihh++){
                inm = im-1+imm;
                inh = ih-1+ihh;
                ipm = im - imm;
                iph = ih - ihh;
                l = gsl_spline_eval(SplineAccessLogT[inm][ing][inh],log10(Tg),AccAccessLogT[inm][ing][inh]);
                lintp += l*fabs((Zg-TableLogZ[ipm])*(fH2g-TableLogFH2[iph]));
            }
        }
        return Pall.ConvertCoolingRate*pow(10.0,lintp*dZinv*dHinv);
    }
    else if( (Zg>TableLogZ[0] && Zg<TableLogZ[SIZE_TABLE_Z-1])
            && (G0g>TableLogG0[0] &&  G0g<TableLogG0[SIZE_TABLE_G0-1]) 
            && !(fH2g>TableLogFH2[0] && fH2g<TableLogFH2[SIZE_TABLE_FH2-1])
            ){
        for(int imm=0; imm<2; imm++){
            for(int igg=0; igg<2; igg++){
                inm = im-1+imm;
                ing = ig-1+igg;
                ipm = im - imm;
                ipg = ig - igg;
                l = gsl_spline_eval(SplineAccessLogT[inm][ing][inh],log10(Tg),AccAccessLogT[inm][ing][inh]);
                lintp += l*fabs((Zg-TableLogZ[ipm])*(G0g-TableLogG0[ipg]));
            }
        }
        return Pall.ConvertCoolingRate*pow(10.0,lintp*dZinv*dGinv);
    }
    else if( !(Zg>TableLogZ[0] && Zg<TableLogZ[SIZE_TABLE_Z-1])
            && !(G0g>TableLogG0[0] &&  G0g<TableLogG0[SIZE_TABLE_G0-1]) 
            && (fH2g>TableLogFH2[0] && fH2g<TableLogFH2[SIZE_TABLE_FH2-1])
            ){
        for(int ihh=0; ihh<2; ihh++){
            inh = ih-1+ihh;
            iph = ih - ihh;
            l = gsl_spline_eval(SplineAccessLogT[inm][ing][inh],log10(Tg),AccAccessLogT[inm][ing][inh]);
            lintp += l*fabs(fH2g-TableLogFH2[iph]);
        }
        return Pall.ConvertCoolingRate*pow(10.0,lintp*dHinv);
    }
    else if( !(Zg>TableLogZ[0] && Zg<TableLogZ[SIZE_TABLE_Z-1])
            && (G0g>TableLogG0[0] &&  G0g<TableLogG0[SIZE_TABLE_G0-1]) 
            && !(fH2g>TableLogFH2[0] && fH2g<TableLogFH2[SIZE_TABLE_FH2-1])
            ){
        for(int igg=0; igg<2; igg++){
            ing = ig-1+igg;
            ipg = ig - igg;
            l = gsl_spline_eval(SplineAccessLogT[inm][ing][inh],log10(Tg),AccAccessLogT[inm][ing][inh]);
            lintp += l*fabs(G0g-TableLogG0[ipg]);
        }
        return Pall.ConvertCoolingRate*pow(10.0,lintp*dGinv);
    }
    else if( (Zg>TableLogZ[0] && Zg<TableLogZ[SIZE_TABLE_Z-1])
            && !(G0g>TableLogG0[0] &&  G0g<TableLogG0[SIZE_TABLE_G0-1]) 
            && !(fH2g>TableLogFH2[0] && fH2g<TableLogFH2[SIZE_TABLE_FH2-1])
            ){
        for(int imm=0; imm<2; imm++){
            inm = im-1+imm;
            ipm = im - imm;
            l = gsl_spline_eval(SplineAccessLogT[inm][ing][inh],log10(Tg),AccAccessLogT[inm][ing][inh]);
            lintp += l*fabs(Zg-TableLogZ[ipm]);
        }
        return Pall.ConvertCoolingRate*pow(10.0,lintp*dZinv);
    }
    else{
        inm = 0;
        if(Zg>=TableLogZ[SIZE_TABLE_Z-1])          inm = SIZE_TABLE_Z-1;
        ing = 0;
        if(G0g>=TableLogG0[SIZE_TABLE_G0-1])       ing = SIZE_TABLE_G0-1;
        inh = 0;
        if(fH2g>=TableLogFH2[SIZE_TABLE_FH2-1])    inh = SIZE_TABLE_FH2-1;
        lintp = gsl_spline_eval(SplineAccessLogT[inm][ing][inh],log10(Tg),AccAccessLogT[inm][ing][inh]);
        return Pall.ConvertCoolingRate*pow(10.0,lintp);
    }
}

static void Linear1D(const int n, double x[n], double y[n], double xi, double *yi){
    for(int i=0; i<n-1; i++){
        if(xi>x[i] && xi<x[i+1]){
            *yi = y[i] + (y[i+1] - y[i])*(xi - x[i])/(x[i+1] - x[i]);
            return ;
        }
    }
}

static void Spline1D(const int n, double x[n], double y[n], double xi, double *yi){
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,n);
    gsl_spline_init(spline, x, y, n);
    *yi = gsl_spline_eval(spline, xi, acc);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return ;
}

static void Spline2D(const int n1, const int n2, 
              double x1[n1], double x2[n2], double y[n1][n2], 
              double x1i, double x2i, double *yi){

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,n2);
    double y1d[n1];
    for(int i1=0; i1<n1; i1++){
        gsl_spline_init(spline, x2, y[i1], n2);
        for(int i2=0; i2<n2; i2++)
            y1d[i1] = gsl_spline_eval(spline, x2i, acc);
    }
    gsl_interp_accel_free(acc);
    gsl_spline_free(spline);

    Spline1D(n1,x1,y1d,x1i,yi);
    return ;
}

static void Spline3D(const int n1, const int n2, const int n3, 
              double x1[n1], double x2[n2], double x3[n3], double y[n1][n2][n3], 
              double x1i, double x2i, double x3i, double *yi){
    
    gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
    gsl_spline *spline1 = gsl_spline_alloc(gsl_interp_cspline,n3);
    double y2d[n1][n2];
    for(int i1=0; i1<n1; i1++){
        for(int i2=0; i2<n2; i2++){
            gsl_spline_init(spline1, x3, y[i1][i2], n3);
            for(int i3=0; i3<n3; i3++){
                y2d[i1][i2] = gsl_spline_eval(spline1, x3i, acc1);
            } 
        }
    }
    gsl_interp_accel_free(acc1);
    gsl_spline_free(spline1);

    gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
    gsl_spline *spline2 = gsl_spline_alloc(gsl_interp_cspline,n2);
    double y1d[n1];
    for(int i1=0; i1<n1; i1++){
        gsl_spline_init(spline2, x2, y2d[i1], n2);
        for(int i2=0; i2<n2; i2++)
            y1d[i1] = gsl_spline_eval(spline2, x2i, acc2);
    }
    gsl_interp_accel_free(acc2);
    gsl_spline_free(spline2);

    Spline1D(n1,x1,y1d,x1i,yi);
    return ;
}

static void Spline4D(const int n1, const int n2, const int n3, const int n4, 
              double x1[n1], double x2[n2], double x3[n3], double x4[n4], double y[n1][n2][n3][n4], 
              double x1i, double x2i, double x3i, double x4i, double *yi){
    
    gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
    gsl_spline *spline1 = gsl_spline_alloc(gsl_interp_cspline,n4);
    double y3d[n1][n2][n3];
    for(int i1=0; i1<n1; i1++){
        for(int i2=0; i2<n2; i2++){
            for(int i3=0; i3<n3; i3++){
                gsl_spline_init(spline1, x4, y[i1][i2][i3], n4);
                for(int i4=0; i4<n4; i4++){
                    y3d[i1][i2][i3] = gsl_spline_eval(spline1, x4i, acc1);
                } 
            }
        }
    }
    gsl_interp_accel_free(acc1);
    gsl_spline_free(spline1);
    
    gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
    gsl_spline *spline2 = gsl_spline_alloc(gsl_interp_cspline,n3);
    double y2d[n1][n2];
    for(int i1=0; i1<n1; i1++){
        for(int i2=0; i2<n2; i2++){
            gsl_spline_init(spline2, x3, y3d[i1][i2], n3);
            for(int i3=0; i3<n3; i3++){
                y2d[i1][i2] = gsl_spline_eval(spline2, x3i, acc2);
            } 
        }
    }
    gsl_interp_accel_free(acc2);
    gsl_spline_free(spline2);

    gsl_interp_accel *acc3 = gsl_interp_accel_alloc();
    gsl_spline *spline3 = gsl_spline_alloc(gsl_interp_cspline,n2);
    double y1d[n1];
    for(int i1=0; i1<n1; i1++){
        gsl_spline_init(spline3, x2, y2d[i1], n2);
        for(int i2=0; i2<n2; i2++)
            y1d[i1] = gsl_spline_eval(spline3, x2i, acc3);
    }
    gsl_interp_accel_free(acc3);
    gsl_spline_free(spline3);

    Spline1D(n1,x1,y1d,x1i,yi);
    return ;
}
#endif // USE_SPAANS2008_COOLING_FUNCTIONS

#ifdef USE_INVERSE_COMPTON_COOLING
static double InverseComptonCooling(const double Ti, const double ne, const double redshift){

#define InvDumpingFactorCompton   (0.001)
    double L_Cc = SQ(SQ(1.e0+redshift))*ne*Ti;
	if(Ti<1.e+3){
        L_Cc *= 0.e0;
    }else if(Ti<1.e+4){
        L_Cc *= exp(1.e-2*(Ti-1.e+4));
    }
	return (Pall.ConvertCoolingRate*5.41e-36*L_Cc);
}
#endif // USE_INVERSE_COMPTON_COOLING
