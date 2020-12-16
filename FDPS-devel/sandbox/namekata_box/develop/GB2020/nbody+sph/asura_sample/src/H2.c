#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "Astro.h"
#include "H2.h"
//#include <cpgplot.h>

static double G0 = 1.0;
static double N0 = 1.e+14;  // cm^-2
static double k = 3.0/4.0; // Draine & Bertoldi (1996)

/*
static gsl_interp_accel *AccAccessLogT;
static gsl_spline *SplineAccessLogT;
void InitializeH2fractionLookUpTable(void){
    for(int it=0; it<SIZE_TABLE_T; it++){
        TableLogT[it] = 0.0;
        TableH2Fraction[it] = 0.0;
    }
    AccAccessLogT = gsl_interp_accel_alloc();
    SplineAccessLogT = gsl_spline_alloc(gsl_interp_cspline,SIZE_TABLE_T);
    gsl_spline_init(SplineAccessLogT,TableLogT,TableH2Fraction,SIZE_TABLE_T); 
}

double InterpolatedTemperatureDependentFH2(double Tg, double G0g){

    // interpolation
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,nT);
    gsl_spline_init(spline, tmp, lamb, nT);
    return gsl_spline_eval(spline, t, acc);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}
*/

// f_H2 = f_H2(T,G0) 
double TemperatureDependentFH2(double Tg, double G0g){
    if(Tg<=100.0)
        return 1.0;
    else{
        double a = log10(MAX(G0g,1.0));
        double b = -2.2*(a+1.0);
        return MAX(1.0*pow(Tg/100.0,b),1.e-6);
    }
}

// H2 formation factor
static double Rf(double Tk){
    double TkMax = 500.0;   // [K] Jura (1976)
    double T0 = 100.0;      // [K] Buch & Zhang (1991)
    double mu = 5.0;        // Papadopoulos et al. (2002)
    if(Tk<=TkMax)
        return  3.5e-17*mu*sqrt(Tk/100.0)/SQ(1.0+Tk/T0);
    else
        return  0.e0;
}

// H2 sheilding factor
static double fs(double N2){
    if(N2<=N0)  return  1.0;
    else        return  pow(N2/N0,-k);
}

// H-H2 collisional destruction rate of H2 (Martin et al. 1998)
static double Rate1(double Tk){
    return 1.0;
}

// H2-H2 collisional destruction rate of H2 (Martin et al. 1998)
static double Rate2(double Tk){
    return 1.0;
}

//  
static double K(double Rc){
    /*
    double k = 3.0/4.0; // Draine & Bertoldi (1996)
    double f3k = 3.0 - k;
    double f4k = 4.0 - k;
    double f42k = 4.0 - 2.0*k;
    return 3.0*pow(2.0,-k+3.0)*(2.0*f3k/f4k-1.0)/(f3k*f42k);
    */
    return 0.976;   // for k=3/4
}

// 
double dn2(double n1, double n2, double Tk){
    double kd = 4.0e-11;    // s^-1
    double n = n1 + 2.0*n2; // n(HI) + 2.0*n(H2)
    double Rc = 1.0;
    double dn2dt = Rf(Tk)*SQ(n)
                -(2.0*Rf(Tk)*n + G0*kd*pow(N0/(n2*Rc),k)*K(Rc))*n2
                -Rate1(Tk)*(n-2.0*n2) - Rate2(Tk)*SQ(n2);
    return dn2dt;
}


