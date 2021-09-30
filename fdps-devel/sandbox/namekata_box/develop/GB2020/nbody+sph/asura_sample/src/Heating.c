#include "config.h"
#include "Cooling.h"
#include "Heating.h"

static char PathCoolingHeatingLog[] = "./CoolingHeatingLog"; 

//#define G_0     (1.0)   //
//#define G_0     (1.7)   // From Draine (1987)
//#define G_0     (1.6e-3)
#define ep      (0.05)  // Efficiency for less than 10000 K.

/*
 * The photoelectronic heating of small grains and PAHs.
 * Non dependence on local metallicity, which means that 
 * this routine assumes the solar abundance for metallicity.
 */

static double MetalIndependentFUV(const double Rho, const double G0, const double Metallicity);
static double MetalDependentFUV(const double Rho, const double G0, const double Metallicity);

static double ConvertFUV = 0.e0;
static double ScaleLengthFUV = 0.e0;

void InitializeFarUltraVioletField(void){

#ifdef USE_FARULTRAVIOLET_HEATING //{

#if (defined(USE_SPAANS1997_COOLING_FUNCTIONS) || defined(USE_SPAANS2008_COOLING_FUNCTIONS))

#ifdef USE_FARULTRAVIOLET_HEATING
    //ConvertFUV = GetUnitEnergyPerUnitTime();
    ConvertFUV = GetUnitCoolingRate();
    //ScaleLengthFUV = ModelParameters.Rd;

    if(MPIGetMyID()==MPI_ROOT_RANK){ // check balance
        MakeDir(PathCoolingHeatingLog);

        fprintf(stderr,"Convert FUV = %g\n",ConvertFUV);
        //fprintf(stderr,"ScaleLengthFUV = %g\n",ScaleLengthFUV);
        
        assert(Pall.ConvertCoolingRate > 0.e0);
        char fname[MaxCharactersInLine];
        FILE *fp;
#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
        double Zsolar = 0.02;
        double metal[] = {0.001*Zsolar,0.01*Zsolar,0.1*Zsolar,0.5*Zsolar,1.0*Zsolar};
        double g0[] = {0.01,1.0,10.0,100.0,1000.0};
        double fH2[] = {0.001,0.01,0.1,0.5,1.0};
        for(int i=0;i<5;i++){
            // Z-dep
            sprintf(fname,"%s/CoolingHeatingBalance.Z%02d.dat",PathCoolingHeatingLog,i);
            FileOpen(fp,fname,"w");
            fprintf(fp,"# (Z,G0,H2)=(%f,%f,%f)\n",metal[i],g0[1],fH2[1]);
            for(double nH=1.e-4;nH<1.e+4;nH*=1.1){
                double Tmin = 1.e1;
                double Rho = nH/Pall.ConvertNumberDensityToCGS;
                double MinValue = fabs(FUV(Tmin,Rho,g0[1],metal[i])
                        - SQ(nH)*InterpolatedCoolingRateSpaans2008(Tmin,metal[i],g0[1],fH2[1]));
                for(double T=1.e1;T<1.e+8;T*=1.1){
                    double v = (FUV(Tmin,Rho,g0[1],metal[i])
                            - SQ(nH)*InterpolatedCoolingRateSpaans2008(T,metal[i],g0[1],fH2[1]));
                    if(fabs(v)<MinValue){
                        MinValue = fmin(MinValue,fabs(v));
                        Tmin = T;
                    }
                    if(v < 0.e0)
                        break;
                }
                fprintf(fp,"%e %e %e %e\n",log10(nH),log10(Tmin),
                        FUV(Tmin,Rho,g0[1],metal[i]),
                        SQ(nH)*InterpolatedCoolingRateSpaans2008(Tmin,metal[i],g0[1],fH2[1]));
            }
            fclose(fp);
            // FUV-dep
            sprintf(fname,"%s/CoolingHeatingBalance.FUV%02d.dat",PathCoolingHeatingLog,i);
            FileOpen(fp,fname,"w");
            fprintf(fp,"# (Z,G0,H2)=(%f,%f,%f)\n",metal[3],g0[i],fH2[1]);
            for(double nH=1.e-4;nH<1.e+4;nH*=1.1){
                double Tmin = 1.e1;
                double Rho = nH/Pall.ConvertNumberDensityToCGS;
                double MinValue = fabs(FUV(Tmin,Rho,g0[i],metal[3])
                        - SQ(nH)*InterpolatedCoolingRateSpaans2008(Tmin,metal[3],g0[i],fH2[1]));
                for(double T=1.e1;T<1.e+8;T*=1.1){
                    double v = (FUV(Tmin,Rho,g0[i],metal[3])
                            - SQ(nH)*InterpolatedCoolingRateSpaans2008(T,metal[3],g0[i],fH2[1]));
                    if(fabs(v)<MinValue){
                        MinValue = fmin(MinValue,fabs(v));
                        Tmin = T;
                    }
                    if(v < 0.e0)
                        break;
                }
                fprintf(fp,"%e %e %e %e\n",log10(nH),log10(Tmin),
                        FUV(Tmin,Rho,g0[i],metal[3]),
                        SQ(nH)*InterpolatedCoolingRateSpaans2008(Tmin,metal[3],g0[i],fH2[3]));
            }
            fclose(fp);
            // H2-dep
            sprintf(fname,"%s/CoolingHeatingBalance.H2%02d.dat",PathCoolingHeatingLog,i);
            FileOpen(fp,fname,"w");
            fprintf(fp,"# (Z,G0,H2)=(%f,%f,%f)\n",metal[3],g0[1],fH2[i]);
            for(double nH=1.e-4;nH<1.e+4;nH*=1.1){
                double Tmin = 1.e1;
                double Rho = nH/Pall.ConvertNumberDensityToCGS;
                double MinValue = fabs(FUV(Tmin,Rho,g0[1],metal[3])
                        - SQ(nH)*InterpolatedCoolingRateSpaans2008(Tmin,metal[3],g0[1],fH2[i]));
                for(double T=1.e1;T<1.e+8;T*=1.1){
                    double v = (FUV(Tmin,Rho,g0[1],metal[3])
                            - SQ(nH)*InterpolatedCoolingRateSpaans2008(T,metal[3],g0[1],fH2[i]));
                    if(fabs(v)<MinValue){
                        MinValue = fmin(MinValue,fabs(v));
                        Tmin = T;
                    }
                    if(v < 0.e0)
                        break;
                }
                fprintf(fp,"%e %e %e %e\n",log10(nH),log10(Tmin),
                        FUV(Tmin,Rho,g0[1],metal[3]),
                        SQ(nH)*InterpolatedCoolingRateSpaans2008(Tmin,metal[3],g0[1],fH2[i]));
            }
            fclose(fp);
        }
#else // USE_SPAANS1997_COOLING_FUNCTIONS
        double G0 = 1.0;
        double metal = 0.02;
        sprintf(fname,"%s/CoolingHeatingBalance.dat",PathCoolingHeatingLog);
        FileOpen(fp,fname,"w");
        for(double nH=1.e-4;nH<1.e+4;nH*=1.1){
            double Tmin = 1.e1;
            double Rho = nH/Pall.ConvertNumberDensityToCGS;
            double fH2 = 0.01; 
            double MinValue = fabs(FUV(Tmin,Rho,G0,metal) - SQ(nH)*MetalInterpolatedCoolingRate(Tmin,metal));
            for(double T=1.e1;T<1.e+8;T*=1.1){
                double v = (FUV(Tmin,Rho,G0,metal) - SQ(nH)*MetalInterpolatedCoolingRate(T,metal));
                if(fabs(v)<MinValue){
                    MinValue = fmin(MinValue,fabs(v));
                    Tmin = T;
                }
                if(v < 0.e0)
                    break;
            }
            fprintf(fp,"%e %e %e %e\n",log10(nH),log10(Tmin),
                    FUV(Tmin,Rho,G0,metal),SQ(nH)*MetalInterpolatedCoolingRate(Tmin,metal));
        }
        fclose(fp);
#endif // USE_SPAANS2008_COOLING_FUNCTIONS
    }
#endif

#endif //(defined(USE_SPAANS1997_COOLING_FUNCTIONS) || defined(USE_SPAANS2008_COOLING_FUNCTIONS))

#endif // USE_FARULTRAVIOLET_HEATING //}
    return;
}

#define Tmin    9000.0
#define Tmax    10000.0
const static double iBase = 1.0/(Tmax-Tmin);

static double CalcFUVEpsilon(const double NumberDensity, const double Temperature, const double ne){

    // How to evaluate the electron number density in various parameters?

    return 0.e0;
}

/*
 * This function evluates the strength of the local far ultra-violet (FUV)
 * field. If the FARULTRAVIOLET_HEATING_CUTOFF flag turns on, the FUV field has
 * the temperature cut off at T=10^4 K. If this flag turns off, the parameter
 * epsilon is evaluated as a function of density, temperature and G0.
 *
 */
double FUV(const double Ti, const double Rho, const double G0, const double Metallicity){

#ifdef USE_FARULTRAVIOLET_HEATING
#ifdef FARULTRAVIOLET_HEATING_CUTOFF
    double HeatingRate = 0.0;
    if(Ti < Tmin){
        HeatingRate += MetalIndependentFUV(Rho,G0,Metallicity);
    }else if(Ti > Tmax){
        HeatingRate += 0.e0;
    } else {
	    //double Ti = Pall.ConvertUtoT*U;
        HeatingRate += MetalIndependentFUV(Rho,G0,Metallicity)*0.5*(1.0-tanh(2.0*(Ti-Tmin)*iBase-1.0));
    }
	return (HeatingRate);
#else 
#define phi_PAH (0.5) // The fiducial value obtained by Wolfire+2003
    double HeatingRate = 0.e0;
    //double epsilon = CalcFUVEpsilon(NumberDensity,Temperature,ne);
#error
#endif
#endif

}

double FUVatR(double R){
    // G0 =1 at R=2Rd;
    return 1.0*exp(-R/ScaleLengthFUV)/exp(-2.0);
}

/*
 * This routine returns the amount of a FUV field.
 */
static double MetalIndependentFUV(const double Rho, const double G0, const double Metallicity){
	double nH = Pall.ConvertNumberDensityToCGS*Rho*(1.e0 - HeliumAbandance - Metallicity);
    return (ConvertFUV*nH*1.e-24*ep*G0);
}

/*
 * This routine returns the amount of a FUV field.
 */
static double MetalDependentFUV(const double Rho, const double G0, const double Metallicity){
	double nH = Pall.ConvertNumberDensityToCGS*Rho*(1.e0 - HeliumAbandance - Metallicity);
    return (ConvertFUV*nH*1.e-24*ep*G0*(50.e0*Metallicity));
}

