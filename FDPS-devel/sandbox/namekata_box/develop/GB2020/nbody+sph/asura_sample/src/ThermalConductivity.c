#include "config.h"

static double ThermalConductivityUnitConvertFactor = 0.e0;
#define _VERBOSE_MODE_
static void WriteThermalConductivityRate(void);
void InitThermalConductivity(void){
#ifdef USE_THERMAL_CONDUCTIVITY // {
    // ThermalConductivityUnitConvertFactor = (1.0/Pall.UnitMass)*CUBE(Pall.UnitTime)/Pall.UnitLength;
#ifdef _VERBOSE_MODE_
    WriteThermalConductivityRate();
#endif //_VERBOSE_MODE_
#endif // USE_THERMAL_CONDUCTIVITY // }

    return ;
}

/*! This function returns the coefficient of thermal conductivity.  According to
 * Parker 1953, the coefficient K is
 *  2.5 \times 10^3 T^{1/2} [erg/K/s/cm] for T<4.47 \times 10^4 K
 *  1.24\times 10^-6 T^{5/2} [erg/K/s/cm] for T > 4.47 \times 10^4 K.
 * The Spitzer value (Spitzer 1962) is 
 *  2.5 \times 10^5 T_4^{1/2} for 10^{-2} < T_4 < 1,
 *  6.0 \times 10^3 T_4^{5/2} for T_4 > 1,
 * where T_4 = T/1.e4. 
 * For the Reinicke Meyer-ter-Vehn Verification test, this function returns a
 * constant value of 1.e21 [erg/ev/s/cm] which is 11.720 [erg/K/s/cm] since
 * ev/k_B = 11604.505(20) K.
 */
#define _Tcrit_Parker_ (4.47e4)
#define _Tcrit_Spitzer_ (1)
double CalcThermalConductivity(const int Index){

#ifdef USE_THERMAL_CONDUCTIVITY //{
#if (defined(TASK_TEST_1D_THERMAL_CONDUCTIVITY)||defined(TASK_TEST_3D_THERMAL_CONDUCTIVITY)) //{
    return 1;
#elif defined(TASK_TEST_REINICKE_MEYER_TER_VEHN) //}//{
    return ThermalConductivityUnitConvertFactor*11.72;
#else // TASK_TEST_REINICKE_MEYER_TER_VEHN_VERIFICATION //}//{
#if THERMAL_CONDUCTIVITY_COEF_TYPE==0 //{
    double T = Pall.ConvertUtoT*Phydro[Index]->UPred;
    if(T<_Tcrit_Parker_){
        return ThermalConductivityUnitConvertFactor*2.5e3*sqrt(T);
    } else {
        return ThermalConductivityUnitConvertFactor*1.24e-6*SQ(T)*sqrt(T);
    }
#elif THERMAL_CONDUCTIVITY_COEF_TYPE==1 // THERMAL_CONDUCTIVITY_COEF_TYPE //}//{ // Use Spitzer type
    double T_4 = 1.e-4*Pall.ConvertUtoT*Phydro[Index]->UPred;
    if(T_4<_Tcrit_Spitzer_){
        return ThermalConductivityUnitConvertFactor*2.5e5*sqrt(T_4);
    } else {
        return ThermalConductivityUnitConvertFactor*6.0e3*SQ(T_4)*sqrt(T_4);
    }
#endif // _ThermalConductivity_ParkerType_ //}
#endif // TASK_TEST_REINICKE_MEYER_TER_VEHN //}
#else // USE_THERMAL_CONDUCTIVITY //}//{
    return 0.e0;
#endif // USE_THERMAL_CONDUCTIVITY //}
}

double CalcThermalConductivityByTemperature(const double T){

#ifdef USE_THERMAL_CONDUCTIVITY //{
#if (defined(TASK_TEST_1D_THERMAL_CONDUCTIVITY)||defined(TASK_TEST_3D_THERMAL_CONDUCTIVITY)) //{
    return 1;
#elif defined(TASK_TEST_REINICKE_MEYER_TER_VEHN) //}//{
    return ThermalConductivityUnitConvertFactor*11.72;
#else // TASK_TEST_REINICKE_MEYER_TER_VEHN_VERIFICATION //}//{
#if THERMAL_CONDUCTIVITY_COEF_TYPE==0 //{
    if(T<_Tcrit_Parker_){
        return ThermalConductivityUnitConvertFactor*2.5e3*sqrt(T);
    } else {
        return ThermalConductivityUnitConvertFactor*1.24e-6*SQ(T)*sqrt(T);
    }
#elif THERMAL_CONDUCTIVITY_COEF_TYPE==1 // THERMAL_CONDUCTIVITY_COEF_TYPE //}//{ // Use Spitzer type
    double T_4 = 1.e-4*T;
    if(T_4<_Tcrit_Spitzer_){
        return ThermalConductivityUnitConvertFactor*2.5e5*sqrt(T_4);
    } else {
        return ThermalConductivityUnitConvertFactor*6.0e3*SQ(T_4)*sqrt(T_4);
    }
#endif // THERMAL_CONDUCTIVITY_COEF_TYPE //}
#endif // TASK_TEST_REINICKE_MEYER_TER_VEHN //}
#else // USE_THERMAL_CONDUCTIVITY //}//{
    return 0.e0;
#endif // USE_THERMAL_CONDUCTIVITY //}
}

static void WriteThermalConductivityRate(void){

    FILE *fp;
    char fname[MaxCharactersInLine];

    MakeDir("./data");
    sprintf(fname,"./data/ThermalConductivity.dat");

    int Nbin = 100;
    double Tmin = 10;
    double Tmax = 1.e8;

    double logTmin = log10(Tmin);
    double logTmax = log10(Tmax);
    double logdT = (logTmax-logTmin)/Nbin;

    FileOpen(fp,fname,"w");
    for(int i=0;i<Nbin;i++){
        double T = pow(10.0,logTmin+logdT*i);

        double K_Parker;
        if(T<_Tcrit_Parker_){
            K_Parker = 2.5e3*sqrt(T);
        } else {
            K_Parker = 1.24e-6*SQ(T)*sqrt(T);
        }

        double T_4 = 1.e-4*T;
        double K_Spitzer;
        if(T_4<_Tcrit_Spitzer_){
            K_Spitzer = 2.5e5*sqrt(T_4);
        } else {
            K_Spitzer = 6.0e3*SQ(T_4)*sqrt(T_4);
        }

        fprintf(fp,"%g %g %g %g %g\n",T,K_Parker,K_Spitzer,
                ThermalConductivityUnitConvertFactor*K_Parker,
                ThermalConductivityUnitConvertFactor*K_Spitzer);
    }

    fclose(fp);

    return ;
}
