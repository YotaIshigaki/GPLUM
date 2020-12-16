#include "config.h"

// These functions are tuned only for the flat, lambda dominated universe.

double CalcCurrentTtoZ(void){

    const static double mTwo_Third = -2.0/3.0;
    const static double Three_Second = 3.0/2.0;

    if(Pall.OmegaM > 0.99){ // The Einstein-de Sitter universe.
        double tau = Three_Second*Pall.TCurrent*Pall.Hubble;

        return (pow(tau,mTwo_Third)-1.0);
    } else { // General case in flat universe.
        //double aml3 = Pall.OmegaM/(1.0-Pall.OmegaM); // a_ml^3
        double aml = cbrt(Pall.OmegaM/(1.0-Pall.OmegaM)); // a_ml
        double tau = Three_Second*Pall.Hubble*sqrt(1.0-Pall.OmegaM)*Pall.TCurrent;

        return (MAX(cbrt(4.0/SQ(exp(tau)-exp(-tau)))/aml-1.0,0.e0));
    }
}

double CalcCurrentZtoT(void){

    const static double Two_Third = 2.0/3.0;
    const static double mThree_Second = -3.0/2.0;

    if(Pall.OmegaM > 0.99){ // The Einstein-de Sitter universe.
        return (Two_Third*pow(1.0+Pall.Redshift,mThree_Second)/Pall.Hubble);
    } else { // General case in flat universe.
        double aml3 = Pall.OmegaM/(1.0-Pall.OmegaM); // a_ml^3 
        double iHSqrt1mOmegaM = 1.0/(Pall.Hubble*sqrt(1.0-Pall.OmegaM));
        double a3_aml3 = 1.0/(aml3*CUBE(1.0+Pall.Redshift));

        return (Two_Third*iHSqrt1mOmegaM*log(sqrt(a3_aml3)+sqrt(1.0+a3_aml3)));
    }
}

double CalcTtoZ(const double TCurrent){

    const static double mTwo_Third = -2.0/3.0;
    const static double Three_Second = 3.0/2.0;

    if(Pall.OmegaM > 0.99){ // The Einstein-de Sitter universe.
        double tau = Three_Second*TCurrent*Pall.Hubble;

        return (pow(tau,mTwo_Third)-1.0);
    } else { // General case in flat universe.
        //double aml3 = Pall.OmegaM/(1.0-Pall.OmegaM); // a_ml^3
        double aml = cbrt(Pall.OmegaM/(1.0-Pall.OmegaM)); // a_ml
        double tau = Three_Second*Pall.Hubble*sqrt(1.0-Pall.OmegaM)*TCurrent;

        return (MAX(cbrt(4.0/SQ(exp(tau)-exp(-tau)))/aml-1.0,0.e0));
    }
}

double CalcZtoT(const double Redshift){

    const static double Two_Third = 2.0/3.0;
    const static double mThree_Second = -3.0/2.0;

    if(Pall.OmegaM > 0.99){ // The Einstein-de Sitter universe.
        return (Two_Third*pow(1.0+Redshift,mThree_Second)/Pall.Hubble);
    } else { // General case in flat universe.
        double aml3 = Pall.OmegaM/(1.0-Pall.OmegaM); // a_ml^3 
        double iHSqrt1mOmegaM = 1.0/(Pall.Hubble*sqrt(1.0-Pall.OmegaM));
        double a3_aml3 = 1.0/(aml3*CUBE(1.0+Redshift));

        return (Two_Third*iHSqrt1mOmegaM*log(sqrt(a3_aml3)+sqrt(1.0+a3_aml3)));
    }
}

double CalcHubbleZ(void){

    return 0.0;
}

double CalcOmegaMZ(void){

    return 0.0;
}

double CalcOmegaLZ(void){

    return 0.0;
}

#if 0
#define Nsize (10)

    Pall.UnitLength = MPC_CGS;
    Pall.UnitTime = 10.0*GIGAYEAR_CGS;
    Pall.UnitMass = 1.e+11*MSUN_CGS;
    double Zs = 100;
    double Ze = 0;
    double Zstep = (Zs-Ze)/Nsize; 
    double HubbleConvertionFactor = GetUnitHubble();
    Pall.OmegaM = 0.3e0;
    Pall.OmegaL = 0.7e0;
    Pall.Hubble = 70.e0;
    //Pall.Hubble = 0.7*HubbleConvertionFactor;
    Pall.Hubble = 0.7*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);
    for(int i=0;i<Nsize;i++){
        Pall.Redshift = Zs - Zstep*(i+1);
        Pall.TCurrent = CalcZtoT();
        fprintf(stderr,"z = %g, t = %g, z = %g\n",Pall.Redshift,Pall.TCurrent,CalcTtoZ());
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(0);
#endif
#if 0
#define Nsize (1000)
    double Zs = 100;
    double Ze = 0;
    double Zstep = (Zs-Ze)/Nsize; 

    /*
    Pall.OmegaM = 0.3e0;
    Pall.OmegaL = 0.7e0;
    Pall.Hubble = 70.e0;
    Pall.hubble = 0.7;
    Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);
    */
    for(int i=0;i<Nsize;i++){
        Pall.Redshift = Zs - Zstep*(i+1);
        Pall.TCurrent = CalcCurrentZtoT();
        /*
        fprintf(stderr,"z = %g, t = %g, %g %g %g %g\n",Pall.Redshift,Pall.TCurrent,
                CalcHubbleParameterZ(),CalcCriticalDensityZ(),
                CalcOverDensityZ(),CalcVirializedDensityZ());
                */
        fprintf(stderr,"%g %g %g %g %g %g\n",Pall.Redshift,Pall.TCurrent,
                CalcHubbleParameterZ(),CalcCriticalDensityZ(),
                CalcOverDensityZ(),CalcVirializedDensityZ());
    }


    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(0);
#endif
