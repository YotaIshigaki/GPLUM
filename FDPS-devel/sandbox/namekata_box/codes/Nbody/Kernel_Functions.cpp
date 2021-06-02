/* Standard headers */
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
/* User-defined headers */
#include "preprocess_keywords.h"
#include "Kernel_Functions.h"

namespace kernel_functions {
// In this namespace, we define representative SPH kernels and 
// their derivatives w.r.t. distance r and smoothing length h.
// 
// [References]
//   (1) Monaghan (1985) [Computer Physics Reports,3,71-124].
//   (2) Monaghan (1992) [ARA+A,30,543].
//   (3) Hernquist & Katz (1989) [ApJS,70,419].
//   (4) Thacker et al.(2000) [MNRAS,319,619]

//* Prototype declaration of static functions
static PS::F64 W_M4(const PS::F64 r, const PS::F64 h);
static PS::F64 dWdr_M4(const PS::F64 r, const PS::F64 h);
static PS::F64 dWdh_M4(const PS::F64 r, const PS::F64 h);
static PS::F64 phi_M4(const PS::F64 r, const PS::F64 h);
static PS::F64 dphidr_M4(const PS::F64 r, const PS::F64 h);
static PS::F64 dphidh_M4(const PS::F64 r, const PS::F64 h);
static PS::F64 dWdr_M4_TC92(const PS::F64 r, const PS::F64 h);
static PS::F64 dWdh_M4_TC92(const PS::F64 r, const PS::F64 h);
static PS::F64 W_M5(const PS::F64 r, const PS::F64 h);
static PS::F64 dWdr_M5(const PS::F64 r, const PS::F64 h);
static PS::F64 dWdh_M5(const PS::F64 r, const PS::F64 h);
static PS::F64 phi_M5(const PS::F64 r, const PS::F64 h);
static PS::F64 dphidr_M5(const PS::F64 r, const PS::F64 h);
static PS::F64 dphidh_M5(const PS::F64 r, const PS::F64 h);
static PS::F64 W_M6(const PS::F64 r, const PS::F64 h);
static PS::F64 dWdr_M6(const PS::F64 r, const PS::F64 h);
static PS::F64 dWdh_M6(const PS::F64 r, const PS::F64 h);
static PS::F64 phi_M6(const PS::F64 r, const PS::F64 h);
static PS::F64 dphidr_M6(const PS::F64 r, const PS::F64 h);
static PS::F64 dphidh_M6(const PS::F64 r, const PS::F64 h);
static PS::F64 W_CT(const PS::F64 r, const PS::F64 h);
static PS::F64 dWdr_CT(const PS::F64 r, const PS::F64 h);
static PS::F64 dWdh_CT(const PS::F64 r, const PS::F64 h);
static PS::F64 phi_CT(const PS::F64 r, const PS::F64 h);
static PS::F64 dphidr_CT(const PS::F64 r, const PS::F64 h);
static PS::F64 dphidh_CT(const PS::F64 r, const PS::F64 h);



/*-------------------------------------------------------------------*/
/////////////////////////// F U N C T I O N ///////////////////////////
///////////////////////////  < W _ S P H >  ///////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 W_SPH(const PS::F64 r, const PS::F64 h){

#if ((SPH_Kernel == M4_Cubic_Spline) || (SPH_Kernel == M4_CS_w_TC92))
   return W_M4(r,h);
#elif (SPH_Kernel == M5_Quartic_Spline)
   return W_M5(r,h);
#elif (SPH_Kernel == M6_Quintic_Spline)
   return W_M6(r,h);
#elif (SPH_Kernel == Core_Triangle)
   return W_CT(r,h);
#else
#error Macro-definition `SPH_Kernel` is invalid.
#endif

}

/*-------------------------------------------------------------------*/
/////////////////////////   F U N C T I O N   /////////////////////////
///////////////////////// < D W D R _ S P H > /////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dWdr_SPH(const PS::F64 r, const PS::F64 h){

#if (SPH_Kernel == M4_Cubic_Spline)
   return dWdr_M4(r,h);
#elif (SPH_Kernel == M4_CS_w_TC92)
   return dWdr_M4_TC92(r,h);
#elif (SPH_Kernel == M5_Quartic_Spline)
   return dWdr_M5(r,h);
#elif (SPH_Kernel == M6_Quintic_Spline)
   return dWdr_M6(r,h);
#elif (SPH_Kernel == Core_Triangle)
   return dWdr_CT(r,h);
#else
#error Macro-definition `SPH_Kernel` is invalid.
#endif

}

/*-------------------------------------------------------------------*/
/////////////////////////   F U N C T I O N   /////////////////////////
///////////////////////// < D W D R _ S P H > /////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dWdh_SPH(const PS::F64 r, const PS::F64 h){

#if ((SPH_Kernel == M4_Cubic_Spline) || (SPH_Kernel == M4_CS_w_TC92))
   return dWdh_M4(r,h);
#elif (SPH_Kernel == M5_Quartic_Spline)
   return dWdh_M5(r,h);
#elif (SPH_Kernel == M6_Quintic_Spline)
   return dWdh_M6(r,h)
#elif (SPH_Kernel == Core_Triangle)
   return dWdh_CT(r,h)
#else
#error Macro-definition `SPH_Kernel` is invalid.
#endif

}

/*-------------------------------------------------------------------*/
//////////////////////////  F U N C T I O N  //////////////////////////
////////////////////////// < P H I _ S P H > //////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 phi_SPH(const PS::F64 r, const PS::F64 h){

#if ((SPH_Kernel == M4_Cubic_Spline) || (SPH_Kernel == M4_CS_w_TC92))
   return phi_M4(r,h);
#elif (SPH_Kernel == M5_Quartic_Spline)
   return phi_M5(r,h);
#elif (SPH_Kernel == M6_Quintic_Spline)
   return phi_M6(r,h);
#elif (SPH_Kernel == Core_Triangle)
   return phi_CT(r,h);
#else
#error Macro-definition `SPH_Kernel` is invalid.
#endif

}

/*-------------------------------------------------------------------*/
///////////////////////     F U N C T I O N     ///////////////////////
/////////////////////// < D P H I D R _ S P H > ///////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dphidr_SPH(const PS::F64 r, const PS::F64 h){

#if ((SPH_Kernel == M4_Cubic_Spline) || (SPH_Kernel == M4_CS_w_TC92))
   return dphidr_M4(r,h);
#elif (SPH_Kernel == M5_Quartic_Spline)
   return dphidr_M5(r,h);
#elif (SPH_Kernel == M6_Quintic_Spline)
   return dphidr_M6(r,h);
#elif (SPH_Kernel == Core_Triangle)
   return dphidr_CT(r,h);
#else
#error Macro-definition `SPH_Kernel` is invalid.
#endif

}

/*-------------------------------------------------------------------*/
///////////////////////     F U N C T I O N     ///////////////////////
/////////////////////// < D P H I D H _ S P H > ///////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dphidh_SPH(const PS::F64 r, const PS::F64 h){

#if ((SPH_Kernel == M4_Cubic_Spline) || (SPH_Kernel == M4_CS_w_TC92))
   return dphidh_M4(r,h);
#elif (SPH_Kernel == M5_Quartic_Spline)
   return dphidh_M5(r,h);
#elif (SPH_Kernel == M6_Quintic_Spline)
   return dphidh_M6(r,h);
#elif (SPH_Kernel == Core_Triangle)
   return dphidh_CT(r,h);
#else
#error Macro-definition `SPH_Kernel` is invalid.
#endif

}

/*-------------------------------------------------------------------*/
/////////////////////////// F U N C T I O N ///////////////////////////
///////////////////////////   < W _ M 4 >   ///////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 W_M4(const PS::F64 r, const PS::F64 h){
  // M4 Cubic spline kernel 
  // (see Eq. (4) in Springel (2005)[MNRAS,364,1105])

  //* Local variables
  PS::F64 u,u2,s,cc;

  u = r/h;
#if (Simulation_Dimension == 1)
  cc=4.0e0/(3.0e0*h);
#elif (Simulation_Dimension == 2)
  cc=4.0e1/(7.0e0*M_PI*h*h);
#elif (Simulation_Dimension == 3)
  cc=8.0e0/(M_PI*h*h*h);
#endif
  if (u <= 0.5e0) {
     u2 = u*u;
     return cc*(1.0e0+u2*6.0e0*(-1.0e0+u));
  }
  else if ((0.5e0 < u) && (u <= 1.0e0)) {
     s = 1.0e0-u;
     return cc*2.0e0*s*s*s;
  }
  else {
     return 0.0e0;
  }
         
}

/*-------------------------------------------------------------------*/
/////////////////////////   F U N C T I O N   /////////////////////////
/////////////////////////  < D W D R _ M 4 >  /////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dWdr_M4(const PS::F64 r, const PS::F64 h){
   // This subroutine gives 
   //     \dfrac{\partial W_M4(r,h)}{\partial r}\dfrac{1}{r},
   // which is used to evaluate \nabla W(r,h).

   //* Local variables
   PS::F64 u,cc;

   u=r/h;
#if (Simulation_Dimension == 1)
   cc = -8.0e0/(h*h);
#elif (Simulation_Dimension == 2)
   cc = -2.4e2/(7.0e0*M_PI*h*h*h);
#elif (Simulation_Dimension == 3)
   cc = -48.0e0/(M_PI*h*h*h*h);
#endif
   if (u <= 0.5e0) {
      return cc*u*(2.0e0-3.0e0*u)/(r);
   }
   else if ((0.5e0 < u) && (u < 1.0e0)) {
      return cc*(1.0e0-u)*(1.0e0-u)/(r);
   }
   else {
      return 0.0e0;
   }
} 

/*-------------------------------------------------------------------*/
/////////////////////////   F U N C T I O N   /////////////////////////
/////////////////////////  < D W D H _ M 4 >  /////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dWdh_M4(const PS::F64 r, const PS::F64 h){
   // This subroutine gives 
   //     \dfrac{\partial W(r,h)}{\partial h}.

   //* Local variables
   PS::F64 u,u2,s,cc;

   u=r/h;
#if (Simulation_Dimension == 1)
   cc=-4.0e0/(3.0e0*h*h);
   if (u <= 0.5e0) {
      u2 = u*u;
      return cc*(1.0e0+u2*(-1.8e1+2.4e1*u));
   }
   else if ((0.5e0 < u) && (u < 1.0e0)) {
      s = 1.0e0-u;
      return cc*2.0e0*s*s*(1.0e0-4.0e0*u);
   }
   else {
      return 0.0e0;
   }
#elif (Simulation_Dimension == 2)
   cc=-8.0e1/(7.0e0*M_PI*h*h*h);
   if (u <= 0.5e0) {
      u2 = u*u;
      return cc*(1.0e0+u2*(-12.0e0+15.0e0*u));
   }
   else if ((0.5e0 < u) && (u < 1.0e0)) {
      s = 1.0e0-u;
      return cc*s*s*(2.0e0-5.0e0*u);
   }
   else {
      return 0.0e0;
   }
#elif (Simulation_Dimension == 3)
   cc=-24.0e0/(M_PI*h*h*h*h);
   if (u <= 0.5e0) {
      u2 = u*u;
      return cc*(1.0e0
                +u2*(-10.0e0
                     +12.0e0*u));
   }
   else if ((0.5e0 < u) && (u < 1.0e0)) {
      s = 1.0e0-u;
      return cc*2.0e0*s*s*(1.0e0-2.0e0*u);
   }
   else {
      return 0.0e0;
   }
#endif

}

/*-------------------------------------------------------------------*/
/////////////////////////   F U N C T I O N   /////////////////////////
/////////////////////////   < P H I _ M 4 >   /////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 phi_M4(const PS::F64 r, const PS::F64 h) {
   // This function is equivalent to the equation (A2) described in 
   // Price & Monaghan (2007)[MNRAS,374,1347], but `h` is defined 
   // so that \phi(r)=1/r^{2} at r=h. Therefore, the coefficients 
   // of the polynomials here are different from those in PM07.
   // (one can obtain the following polynomials by the transformation
   //  of {`q` --> 2u, AND, h --> h/2} in PM07.)

   //* Local variables
   PS::F64 u,u2;

   u = r/h;
   if (u <= 0.5e0) {
      u2  = u*u;
      return (-2.8e0
              +u2*(16.0e0/3.0e0
                  +u2*(-9.6e0
                       +6.4e0*u)))/h;
   }
   else if ((0.5e0 < u) && (u < 1.0e0)) {
      u2 = u*u;
      return (1.0e0/(15.0e0*u)
             -3.2e0
             +u2*(32.0e0/3.0e0
                 +u*(-16.0e0
                     +u*(9.6e0
                        -32.0e0*u/15.0e0))))/h;
   }
   else {
      return -1.0e0/r;
   }

}

/*-------------------------------------------------------------------*/
////////////////////////    F U N C T I O N    ////////////////////////
//////////////////////// < D P H I D R _ M 4 > ////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dphidr_M4(const PS::F64 r, const PS::F64 h) {
   // This function is equivalent to the equation (A1) described in 
   // Price & Monaghan (2007)[MNRAS,374,1347].

   //* Local variables
   PS::F64 u,u2,cc;

   u = r/h;
   cc = 32.0e0/(h*h);
   if (u <= 0.5e0) {
      u2 = u*u;
      return cc*u*(1.0e0/3.0e0
               +u2*(-1.2e0+u));
   }
   else if ((0.5e0 < u) && (u < 1.0e0)) {
      u2 = u*u;
      return cc*(-1.0e0/(480.0e0*u2)
                +u*(2.0e0/3.0e0
                   +u*(-1.5e0
                       +u*(1.2e0
                          -u/3.0e0))));
   }
   else {
      return 1.0e0/(r*r);
   }

}

/*-------------------------------------------------------------------*/
////////////////////////    F U N C T I O N    ////////////////////////
//////////////////////// < D P H I D H _ M 4 > ////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dphidh_M4(const PS::F64 r, const PS::F64 h) {
   // This function is equivalent to the equation (A3) described in 
   // Price & Monaghan (2007)[MNRAS,374,1347]. Note that we should 
   // convert L.H.S. of Eq.(A3) as well as R.H.S. of Eq.(A3) because
   // L.H.S. include h, which should be converted to (H/2)!

   //* Local variables
   PS::F64 u,u2,cc;

   u = r/h;
   cc = 2.0e0/(h*h);
   if (u <= 0.5e0) {
      u2 = u*u;
      return cc*(1.4e0
                +u2*(-8.0e0
                     +u2*(24.0e0
                         -19.2e0*u)));
   }
   else if ((0.5e0 < u) && (u < 1.0e0)) {
      u2 = u*u;
      return cc*(1.6e0
                +u2*(-16.0e0
                     +u*(32.0e0
                        +u*(-24.0e0
                            +6.4e0*u))));
   }
   else {
      return 0.0e0;
   }

}

/*-------------------------------------------------------------------*/
/////////////////////       F U N C T I O N       /////////////////////
///////////////////// < D W D R _ M 4 _ T C 9 2 > /////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dWdr_M4_TC92(const PS::F64 r, const PS::F64 h) {
   // This subroutine gives 
   //     \dfrac{\partial W_M4(r,h)}{\partial r}\dfrac{1}{r},
   // which is used to evaluate \nabla W(r,h).
   // But, a modification is applied, if r/h is less than 1/3,
   // according to Thomas & Couchman (1992)[MNRAS,257,11].

   //* Local variables
   PS::F64 u,cc;

   u=r/h;
#if (Simulation_Dimension == 1)
   cc = -8.0e0/(h*h);
#elif (Simulation_Dimension == 2)
   cc = -2.4e2/(7.0e0*M_PI*h*h*h);
#elif (Simulation_Dimension == 3)
   cc = -48.0e0/(M_PI*h*h*h*h);
#endif
   if (u <= 1.0e0/3.0e0) {
      return cc*(1.0e0/3.0e0)/(r);
   }
   else if ((1.0e0/3.0e0 < u) && (u <= 0.5e0)) {
      return cc*u*(2.0e0-3.0e0*u)/(r);
   }
   else if ((0.5e0 < u) && (u < 1.0e0)) {
      return cc*(1.0e0-u)*(1.0e0-u)/(r);
   }
   else {
      return 0.0e0;
   }

}

/*-------------------------------------------------------------------*/
/////////////////////////   F U N C T I O N   /////////////////////////
/////////////////////////  < D W D H _ M 4 >  /////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dWdh_M4_TC92(const PS::F64 r, const PS::F64 h) {
   // This subroutine gives 
   //     \dfrac{\partial W(r,h)}{\partial h}.
   // But, a modification is applied, if r/h is less than 1/3,
   // (similar to Thomas & Couchman (1992)[MNRAS,257,11])

   //* Local variables
   PS::F64 u,u2,s,cc;

   // [Note]
   //   This function does not support multidimention except for 3D.
   u=r/h;
   cc=-24.0e0/(M_PI*h*h*h*h);
   if (u <= 1.0e0/3.0e0) {
      return cc*(1.0e0/3.0e0);
   }
   else if ((1.0e0/3.0e0 < u) && (u <= 0.5e0)) {
      u2 = u*u;
      return cc*(1.0e0
                +u2*(-10.0e0
                     +12.0e0*u));
   }
   else if ((0.5e0 < u) && (u < 1.0e0)) {
      s = 1.0e0-u;
      return cc*2.0e0*s*s*(1.0e0-2.0e0*u);
   }
   else {
      return 0.0e0;
   }

}

/*-------------------------------------------------------------------*/
/////////////////////////// F U N C T I O N ///////////////////////////
///////////////////////////   < W _ M 5 >   ///////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 W_M5(const PS::F64 r, const PS::F64 h) {
   // M5 Quartic spline kernel

   //* Local parameters
   const PS::F64 a=1.0e0/5.0e0;
   const PS::F64 b=3.0e0/5.0e0;
   //* Local variables
   PS::F64 u,u2,s,cc;

   u=r/h;
#if (Simulation_Dimension == 1)
   cc=(pow(2.5e0,5.0e0))/(24.0e0*h);
#elif (Simulation_Dimension == 2)
   cc=(pow(2.5e0,6.0e0))*9.6e1/(1.119e3*M_PI*h*h);
#elif (Simulation_Dimension == 3)
   cc=(pow(2.5e0,7.0e0))/(20.0e0*M_PI*h*h*h);
#endif
   if (u <= a) {
      u2 = u*u;
      return cc*(0.368e0+u2*(-2.4e0+6.0e0*u2));
   }
   else if ((a < u) && (u <= b)) {
      return cc*(0.352e0
                +u*(0.32e0
                   +u*(-4.8e0
                       +u*(8.0e0
                          -4.0e0*u))));
   }
   else if ((b < u) && (u <= 1.0e0)) {
      s = 1.0e0-u;
      return cc*s*s*s*s;
   }
   else {
      return 0.0e0;
   }

}

/*-------------------------------------------------------------------*/
/////////////////////////   F U N C T I O N   /////////////////////////
/////////////////////////  < D W D R _ M 5 >  /////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dWdr_M5(const PS::F64 r, const PS::F64 h) {
   // M5 Quartic spline kernel

   //* Local parameters
   const PS::F64 a=1.0e0/5.0e0;
   const PS::F64 b=3.0e0/5.0e0;
   //* Local variables
   PS::F64 u,cc;

   u=r/h;
#if (Simulation_Dimension == 1)
   cc=-(pow(2.5e0,5.0e0))/(6.0e0*h*h);
#elif (Simulation_Dimension == 2)
   cc=-(pow(2.5e0,6.0e0))*3.84e2/(1.119e3*M_PI*h*h*h);
#elif (Simulation_Dimension == 3)
   cc=-(pow(2.5e0,7.0e0))/(5.0e0*M_PI*h*h*h*h);
#endif
   if (u <= a) {
      return cc*u*(1.2e0-6.0e0*u*u)/(r);
   }
   else if ((a < u) && (u <= b)) {
      return cc*(-0.08e0
                 +u*(2.4e0
                    +u*(-6.0e0
                        +4.0e0*u)))/(r);
   }
   else if ((b < u) && (u <= 1.0e0)) {
      u = 1.0e0-u;
      return cc*u*u*u/(r);
   }
   else {
      return 0.0e0;
   }

}

/*-------------------------------------------------------------------*/
//////////////////////////  F U N C T I O N  //////////////////////////
////////////////////////// < D W D H _ M 5 > //////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dWdh_M5(const PS::F64 r, const PS::F64 h) {
   // M5 Quartic spline kernel

   //* Local parameters
   PS::F64 a=1.0e0/5.0e0;
   PS::F64 b=3.0e0/5.0e0;
   //* Local variables
   PS::F64 u,u2,s,cc;

   u=r/h;
#if (Simulation_Dimension == 1)
   cc=-5.0e0*(pow(2.5e0,5.0e0))/(24.0e0*h*h);
   if (u <= a) {
      u2 = u*u;
      return cc*(7.36e-2+u2*6.0e0*(-0.24e0+u2));
   }
   else if ((a < u) && (u <= b)) {
      return cc*(7.04e-2+u*(0.128e0+u*(-2.88e0+u*4.0e0*(1.6e0-u))));
   }
   else if ((b < u) && (u <= 1.0e0)) {
      s = 1.0e0-u;
      return cc*(0.2e0-u)*s*s*s;
   }
   else {
      return 0.0e0;
   }
#elif (Simulation_Dimension == 2)
   cc=-5.76e2*(pow(2.5e0,6.0e0))/(1.199e3*M_PI*h*h*h);
   if (u <= a) {
      u2 = u*u;
      return cc*(4.6e1/3.75e2+u2*(-1.6e0+6.0e0*u2));
   }
   else if ((a < u) && (u <= b)) {
      return cc*(4.4e1/3.75e2
                +u*(0.16e0
                   +u*(-3.2e0
                      +u*4.0e0*(5.0e0/3.0e0-u))));
   }
   else if ((b < u) && (u <= 1.0e0)) {
      return cc*(1.0e0/3.0e0
                +u*(-2.0e0
                    +u*(4.0e0
                        +u*(-1.0e1/3.0e0+u))));
   }
   else {
      return 0.0e0;
   }
#elif (Simulation_Dimension == 3)
   cc=-7.0e0*(pow(2.5e0,7.0e0))/(20.0e0*M_PI*h*h*h*h);
   if (u <= a) {
      u2 = u*u;
      return cc*(138.0e0/875.0e0
                +u2*(-12.0e0/7.0e0
                     +6.0e0*u2));
   }
   else if ((a < u) && (u <= b)) {
      return cc*((4.0e0*(3.3e1
                        +5.0e0*u*(8.0e0
                                 -2.5e1*u*(6.0e0
                                          +u*(-1.2e1
                                              +7.0e0*u)))))/8.75e2);
   }
   else if ((b < u) && (u <= 1.0e0)) {
      s = 1.0e0/7.0e0;
      return cc*(3.0e0*s
                +u*(-1.6e1*s
                    +u*(3.0e1*s+(-2.4e1*s+u)*u)));
   }
   else {
      return 0.0e0;
   }
#endif
   
}

/*-------------------------------------------------------------------*/
/////////////////////////// F U N C T I O N ///////////////////////////
/////////////////////////// < P H I _ M 5 > ///////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 phi_M5(const PS::F64 r, const PS::F64 h) {
   //* Local parameters
   const PS::F64 a=1.0e0/5.0e0;
   const PS::F64 b=3.0e0/5.0e0;
   //* Local variables
   PS::F64 u,u2,cc;

   u=r/h;
   cc=-pow(2.5e0,7.0e0)/(5.0e0*h);
   if (u <= a) {
      u2 = u*u;
      return cc*(1.199e3/4.6875e4
                +u2*(-2.3e1/3.75e2
                     +u2*(0.12e0
                         -u2/7.0e0)));
   }
   else if ((a < u) && (u <= b)) {
      u2 = u*u;
      return cc*(2.0e0/(1.640625e6*u)
                +1.198e3/4.6875e4
                +u2*(-2.2e1/3.75e2
                     +u*(-2.0e0/7.5e1
                        +u*(0.24e0
                            +u*(-4.0e0/1.5e1
                                +2.0e0*u/2.1e1)))));
   }
   else if ((b < u) && (u <= 1.0e0)) {
      u2 = u*u;
      return cc*(-4.37e2/(3.28125e5*u)
                 +1.0e0/3.0e1
                 +u2*(-1.0e0/6.0e0
                      +u*(1.0e0/3.0e0
                         +u*(-0.3e0
                             +u*(2.0e0/1.5e1
                                -u/4.2e1)))));
   }
   else {
      return -1.0e0/r;
   }
   
}

/*-------------------------------------------------------------------*/
////////////////////////    F U N C T I O N    ////////////////////////
//////////////////////// < D P H I D R _ M 5 > ////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dphidr_M5(const PS::F64 r, const PS::F64 h) {
   //* Local parameters
   const PS::F64 a=1.0e0/5.0e0;
   const PS::F64 b=3.0e0/5.0e0;
   //* Local variables
   PS::F64 u,u2,cc;

   u=r/h;
   cc=pow(2.5e0,7.0e0)/(5.0e0*h*h);
   if (u <= a) {
      u2 = u*u;
      return cc*u*(4.6e1/3.75e2
                  +u2*(-0.48e0
                       +6.0e0*u2/7.0e0));
   }
   else if ((a < u) && (u <= b)) {
      u2 = u*u;
      return cc*(2.0e0/(1.640625e6*u2)
                +4.4e1*u/3.75e2
                +u2*(0.08e0
                    +u*(-0.96e0
                        +u*(4.0e0/3.0e0
                           -4.0e0*u/7.0e0))));
   }
   else if ((b < u) && (u <= 1.0e0)) {
      u2 = u*u;
      return cc*(-4.37e2/(3.28125e5*u2)
                 +u/3.0e0
                 +u2*(-1.0e0
                      +u*(1.2e0
                         +u*(-2.0e0/3.0e0
                             +u/7.0e0))));
   }
   else {
      return 1.0e0/(r*r);
   }
}

/*-------------------------------------------------------------------*/
////////////////////////    F U N C T I O N    ////////////////////////
//////////////////////// < D P H I D H _ M 5 > ////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dphidh_M5(const PS::F64 r, const PS::F64 h) {
   //* Local parameters
   const PS::F64 a=1.0e0/5.0e0;
   const PS::F64 b=3.0e0/5.0e0;
   //* Local variables
   PS::F64 u,u2,cc;

   u=r/h;
   cc=pow(2.5e0,7.0e0)/(5.0e0*h*h);
   if (u <= a) {
      u2 = u*u;
      return cc*(1.199e3/4.6875e4
                +u2*(-0.184e0
                     +u2*(0.6e0
                         -u2)));
   }
   else if ((a < u) && (u <= b)) {
      u2 = u*u;
      return cc*(1.198e3/4.6875e4
                +u2*(-0.176e0
                     +u*(-8.0e0/7.5e1
                        +u*(1.2e0
                            +u*(-1.6e0
                                +2.0e0*u/3.0e0)))));
   }
   else if ((b < u) && (u <= 1.0e0)) {
      u2 = u*u;
      return cc*(1.0e0/3.0e1
                +u2*(-0.5e0
                     +u*(4.0e0/3.0e0
                        +u*(-1.5e0
                            +u*(0.8e0
                               -u/6.0e0)))));
   }
   else {
      return 0.0e0;
   }

}

/*-------------------------------------------------------------------*/
/////////////////////////// F U N C T I O N ///////////////////////////
///////////////////////////   < W _ M 6 >   ///////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 W_M6(const PS::F64 r, const PS::F64 h) {
   // M6 Quintic spline kernel
   //* Local parameters
   const PS::F64 a=1.0e0/3.0e0;
   const PS::F64 b=2.0e0/3.0e0;
   //* Local variables
   PS::F64 u,u2,s,cc;

   u=r/h;
#if (Simulation_Dimension == 1)
   cc=pow(3.0e0,5.0e0)/(4.0e1*h);
#elif (Simulation_Dimension == 2)
   cc=pow(3.0e0,7.0e0)*7.0e0/(4.78e2*M_PI*h*h);
#elif (Simulation_Dimension == 3)
   cc=pow(3.0e0,7.0e0)/(4.0e1*M_PI*h*h*h);
#endif
   if (u <= a) {
      u2 = u*u;
      return cc*(2.2e1/8.1e1
                +u2*(-2.0e1/9.0e0
                     +u2*1.0e1*(1.0e0-u)));
   }
   else if ((a < u) && (u <= b)) {
      return cc*(1.7e1/8.1e1
                +u*(2.5e1/2.7e1
                   +u*(-7.0e1/9.0e0
                       +u*(5.0e1/3.0e0
                          +u*5.0e0*(-3.0e0+u)))));
   }
   else if ((b < u) && (u <= 1.0e0)) {
      s = 1.0e0-u;
      u2 = s*s;
      return cc*u2*u2*s;
   }
   else {
      return 0.0e0;
   }

}

/*-------------------------------------------------------------------*/
/////////////////////////   F U N C T I O N   /////////////////////////
/////////////////////////  < D W D R _ M 6 >  /////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dWdr_M6(const PS::F64 r, const PS::F64 h) {
   // M6 Quintic spline kernel
   //* Local parameters
   const PS::F64 a=1.0e0/3.0e0;
   const PS::F64 b=2.0e0/3.0e0;
   //* Local variables
   PS::F64 u,u2,s,cc;

   u=r/h;
#if (Simulation_Dimension == 1)
   cc=-3.0e0**5.0e0/(8.0e0*h*h);
#elif (Simulation_Dimension == 2)
   cc=-pow(3.0e0,7.0e0)*3.5e1/(4.78e2*M_PI*h*h*h);
#elif (Simulation_Dimension == 3)
   cc=-pow(3.0e0,7.0e0)/(8.0e0*M_PI*h*h*h*h);
#endif
   if (u <= a) {
      u2 = u*u;
      return cc*u*(8.0e0/9.0e0
                  +u2*(-8.0e0
                       +1.0e1*u))/(r);
   }
   else if ((a < u) && (u <= b)) {
      return cc*(-5.0e0/2.7e1
                 +u*(2.8e1/9.0e0
                    +u*(-1.0e1
                        +u*(1.2e1
                           -5.0e0*u))))/(r);
   }
   else if ((b < u) && (u <= 1.0e0)) {
      s = 1.0e0-u;
      s = s*s;
      return cc*s*s/(r);
   }
   else {
      return 0.0e0;
   }

}

/*-------------------------------------------------------------------*/
//////////////////////////  F U N C T I O N  //////////////////////////
////////////////////////// < D M D H _ M 6 > //////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dWdh_M6(const PS::F64 r, const PS::F64 h) {
   //* Local parameters
   const PS::F64 a=1.0e0/3.0e0;
   const PS::F64 b=2.0e0/3.0e0;
   //* Local variables
   PS::F64 u,u2,cc;

   u=r/h;
#if (Simulation_Dimension == 1)
   cc=-6.0e0*(pow(3.0e0,5.0e0))/(4.0e1*h*h);
   if (u <= a) {
      u2 = u*u;
      return cc*(1.1e1/2.43e2
                +u2*(-1.0e1/9.0e0
                     +u2*(2.5e1/3.0e0
                         -1.0e1*u)));
   }
   else if ((a < u) && (u <= b)) {
      return cc*(1.7e1/4.86e2
                +u*(2.5e1/8.1e1
                   +u*(-3.5e1/9.0e0
                       +u*(1.0e2/9.0e0
                          +u*5.0e0*(-2.5e0+u)))));
   }
   else if ((b < u) && (u <= 1.0e0)) {
      return cc*(1.0e0/6.0e0
                +u*(-5.0e0/3.0e0
                    +u*(5.0e0
                       +u*(-2.0e1/3.0e0
                          +u*(2.5e1/6.0e0-u)))));
   }
   else {
      return 0.0e0;
   }
#elif (Simulation_Dimension == 2)
   cc=-4.9e1*(pow(3.0e0,7.0e0))/(4.78e2*M_PI*h*h*h);
   if (u <= a) {
      u2 = u*u;
      return cc*(4.4e1/5.67e2
                +u2*(-8.0e1/6.3e1
                     +u2*1.0e1*(6.0e0/7.0e0-u)));
   }
   else if ((a < u) && (u <= b)) {
      return cc*(3.4e1/5.67e2
                +u*(2.5e1/6.3e1
                   +u*(-4.0e1/9.0e0
                       +u*(2.5e2/2.1e1
                          +u*5.0e0*(-1.8e1/7.0e0+u)))));
   }
   else if ((b < u) && (u <= 1.0e0)) {
      return cc*(2.0e0/7.0e0
                +u*(-1.5e1/7.0e0
                    +u*(4.0e1/7.0e0
                       +u*(-5.0e1/7.0e0
                           +u*(3.0e1/7.0e0-u)))));
   }
   else {
      return 0.0e0;
   }
#elif (Simulation_Dimension == 3)
   cc=-pow(3.0e0,7.0e0)/(5.0e0*M_PI*h*h*h*h);
   if (u <= a) {
      u2 = u*u;
      return cc*(1.1e1/1.08e2
                +u2*(-2.5e1/1.8e1
                     +u2*(8.75e0
                         -1.0e1*u)));
   }
   else if ((a < u) && (u <= b)) {
      return cc*(17.0e0/216.0e0
                +u*(25.0e0/54.0e0
                   +u*(-175.0e0/36.0e0
                       +u*(12.5e0
                          +u*(-13.125e0
                              +5.0e0*u)))));
   }
   else if ((b < u) && (u <= 1.0e0)) {
      return cc*(0.375e0
                +u*(-2.5e0
                    +u*(6.25e0
                       +u*(-7.5e0
                           +(4.375e0-u)*u))));
   }
   else {
      return 0.0e0;
   }
#endif
}

/*-------------------------------------------------------------------*/
/////////////////////////// F U N C T I O N ///////////////////////////
/////////////////////////// < P H I _ M 6 > ///////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 phi_M6(const PS::F64 r, const PS::F64 h) {
   //* Local parameters
   const PS::F64 a=1.0e0/3.0e0;
   const PS::F64 b=2.0e0/3.0e0;
   //* Local variables
   PS::F64 u,u2,cc;

   u=r/h;
   cc=-pow(3.0e0,7.0e0)/(10.0e0*h);
   if (u <= a) {
      u2 = u*u;
      return cc*(2.39e2/1.5309e4
                +u2*(-1.1e1/2.43e2
                     +u2*(1.0e0/9.0e0
                         +u2*(-5.0e0/2.1e1
                              +5.0e0*u/2.8e1))));
   }
   else if ((a < u) && (u <= b)) {
      u2 = u*u;
      return cc*(5.0e0/(367416.0e0*u)
                +473.0e0/30618.0e0
                +u2*(-17.0e0/486.0e0
                     +u*(-25.0e0/324.0e0
                         +u*(7.0e0/18.0e0
                            +u*(-5.0e0/9.0e0
                                +u*(5.0e0/14.0e0
                                   -5.0e0*u/56.0e0))))));
   }
   else if ((b < u) && (u <= 1.0e0)) {
      u2 = u*u;
      return cc*(-169.0e0/(122472.0e0*u)
                 +1.0e0/42.0e0
                 +u2*(-1.0e0/6.0e0
                      +u*(5.0e0/12.0e0
                         +u*(-0.5e0
                             +u*(1.0e0/3.0e0
                                +u*(-5.0e0/42.0e0
                                    +u/56.0e0))))));
	}
   else {
      return -1.0e0/r;
   }

}

/*-------------------------------------------------------------------*/
////////////////////////    F U N C T I O N    ////////////////////////
//////////////////////// < D P H I D R _ M 6 > ////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dphidr_M6(const PS::F64 r, const PS::F64 h) {
   //* Local parameters
   const PS::F64 a=1.0e0/3.0e0;
   const PS::F64 b=2.0e0/3.0e0;
   //* Local variables
   PS::F64 u,u2,cc;

   u=r/h;
   cc=pow(3.0e0,7.0e0)/(10.0e0*h*h);
   if (u <= a) {
      u2 = u*u;
      return cc*u*(22.0e0/243.0e0
                  +u2*(-4.0e0/9.0e0
                       +u2*(10.0e0/7.0e0
                           -1.25e0*u)));
   }
   else if ((a < u) && (u <= b)) {
      u2 = u*u;
      return cc*(5.0e0/(367416.0e0*u2)
                +u*(17.0e0/243.0e0
                   +u*(25.0e0/108.0e0
                      +u*(-14.0e0/9.0e0
                          +u*(25.0e0/9.0e0
                             +u*(-15.0e0/7.0e0
                                 +0.625e0*u))))));
   }
   else if ((b < u) && (u <= 1.0e0)) {
      u2 = u*u;
      return cc*(-169.0e0/(122472.0e0*u2)
                 +u*(1.0e0/3.0e0
                    +u*(-1.25e0
                        +u*(2.0e0
                           +u*(-5.0e0/3.0e0
                               +u*(5.0e0/7.0e0
                                  -u/8.0e0))))));
   }
   else {
      return 1.0e0/(r*r);
   }

}

/*-------------------------------------------------------------------*/
////////////////////////    F U N C T I O N    ////////////////////////
//////////////////////// < D P H I D H _ M 6 > ////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dphidh_M6(const PS::F64 r, const PS::F64 h) {
   //* Local parameters
   const PS::F64 a=1.0e0/3.0e0;
   const PS::F64 b=2.0e0/3.0e0;
   //* Local variables
   PS::F64 u,u2,cc;

   u=r/h;
   cc=pow(3.0e0,7.0e0)/(10.0e0*h*h);
   if (u <= a) {
      u2 = u*u;
      return cc*(2.39e2/1.5309e4
                +u2*(-1.1e1/8.1e1
                     +u2*(5.0e0/9.0e0
                         +u2*(-5.0e0/3.0e0
                              +1.0e1*u/7.0e0))));
   }
   else if ((a < u) && (u <= b)) {
      u2 = u*u;
      return cc*(4.73e2/3.0618e4
                +u2*(-1.7e1/1.62e2
                     +u*(-2.5e1/8.1e1
                         +u*(3.5e1/1.8e1
                            +u*(-1.0e1/3.0e0
                                +u*(2.5e0
                                   -5.0e0*u/7.0e0))))));
   }
   else if ((b < u) && (u <= 1.0e0)) {
      u2 = u*u;
      return cc*(1.0e0/4.2e1
                +u2*(-0.5e0
                     +u*(5.0e0/3.0e0
                        +u*(-2.5e0
                            +u*(2.0e0
                               +u*(-5.0e0/6.0e0
                                   +u/7.0e0))))));
   }
   else {
      return 0.0e0;
   }

}

/*-------------------------------------------------------------------*/
//////////////////////////  F U N C T I O N  //////////////////////////
//////////////////////////    < W _ C T >    //////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 W_CT(const PS::F64 r, const PS::F64 h) {
   // The core triangle kernel function described in
   // Read et al. (2010)[MNRAS,405,1503].
   // See Eq.(45) in their paper.

   //* Local parameters
   const PS::F64 a=1.0e0/3.0e0;
   const PS::F64 c0=1.0043895747599452e0;
   const PS::F64 b=1.2222222222222223e0;
   //* Local variables
   PS::F64 u,cc;

   u = r/h;
#if (Simulation_Dimension == 1)
   cc = 2.16e2/(1.7e2*h);
#elif (Simulation_Dimension == 2)
   cc = 1.296e3/(2.3e2*M_PI*h*h);
#elif (Simulation_Dimension == 3)
   cc = 8.0e0/(M_PI*c0*h*h*h);
#endif
   if (u <= a) {
      return cc*(b-2.0e0*u);
   }
   else if ((a < u) && (u <= 0.5e0)) {
      return cc*(1.0e0 - 6.0e0*u*u + 6.0e0*u*u*u);
   }
   else if ((0.5e0 < u) && (u <= 1.0e0)) {
      return cc*2.0e0*(1.0e0-u)*(1.0e0-u)*(1.0e0-u);
   }
   else {
      return 0.0e0;
   }

}

/*-------------------------------------------------------------------*/
/////////////////////////   F U N C T I O N   /////////////////////////
/////////////////////////  < D W D R _ C T >  /////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dWdr_CT(const PS::F64 r, const PS::F64 h) {
   // This subroutine gives
   //     \dfrac{\partial W(r,h)}{\partial r}\dfrac{1}{r}
   // for the CT kernel described in Read et al. (2010)[MNRAS,405,1503].
   //* Local parameters
   const PS::F64 a=1.0e0/3.0e0;
   const PS::F64 c0=1.0043895747599452e0;
   //* Local variables
   PS::F64 u,cc;

   u = r/h;
#if (Simulation_Dimension == 1) 
   cc = -1.296e3/(1.7e2*h*h);
#elif (Simulation_Dimension == 2) 
   cc = -7.776e3/(2.3e2*M_PI*h*h*h);
#elif (Simulation_Dimension == 3) 
   cc = -48.0e0/(M_PI*c0*h*h*h*h);
#endif
   if (u <= a) {
      return cc*a/(r);
   }
   else if ((a < u) && (u <= 0.5e0)) {
      return cc*(2.0e0*u - 3.0e0*u*u)/(r);
   }
   else if ((0.5e0 < u) && (u <= 1.0e0)) {
      return cc*(1.0e0-u)*(1.0e0-u)/(r);
   }
   else {
      return 0.0e0;
   }
   // Note that this definition is the same as the kernel gradient
   // used in Thomas & Couchman (1992)[MNRAS,257,11]

}

/*-------------------------------------------------------------------*/
//////////////////////////  F U N C T I O N  //////////////////////////
////////////////////////// < D W D H _ C T > //////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dWdh_CT(const PS::F64 r, const PS::F64 h) {
   //* Local parameters
   const PS::F64 a=1.0e0/3.0e0;
   const PS::F64 a0=-2.6666666666666665e0;
   const PS::F64 c0=1.0043895747599452e0;
   const PS::F64 b=1.2222222222222223e0;
   // Note that 
   //   (1) a0 is actually -16\alpha + 24\alpha^{2} for \alpha=1/3.
   //* Local variables
   PS::F64 u,u2,s,cc;

   u=r/h;
#if (Simulation_Dimension == 1)
   cc=-2.16e2/(1.7e2*h*h);
   if (u <= a) {
      return cc*(1.1e1/9.0e0-4.0e0*u);
   }
   else if ((a < u) && (u <= 0.5e0)) {
      u2 = u*u;
      return cc*(1.0e0+u2*(-1.8e1+2.4e1*u));
   }
   else if ((0.5e0 < u) && (u < 1.0e0)) {
      s = 1.0e0-u;
      return cc*2.0e0*s*s*(1.0e0-4.0e0*u);
   }
   else {
      return 0.0e0;
   }
#elif (Simulation_Dimension == 2)
   cc=-2.592e3/(2.3e2*M_PI*h*h*h);
   if (u <= a) {
      return cc*(1.1e1/9.0e0-3.0e0*u);
   }
   else if ((a < u) && (u <= 0.5e0)) {
      u2 = u*u;
      return cc*(1.0e0+u2*(-1.2e1+1.5e1*u));
   }
   else if ((0.5e0 < u) && (u < 1.0e0)) {
      s = 1.0e0-u;
      return cc*s*s*(2.0e0-5.0e0*u);
   }
   else {
      return 0.0e0;
   }
#elif (Simulation_Dimension == 3)
   cc=-24.0e0/(M_PI*c0*h*h*h*h);
   if (u <= a) {
      return cc*(a0*u+b);
   }
   else if ((a < u) && (u <= 0.5e0)) {
      return cc*(1.0e0
                -10.0e0*u*u
                +12.0e0*u*u*u);
   }
   else if ((0.5e0 < u) && (u < 1.0e0)) {
      return cc*2.0e0*(1.0e0-u)*(1.0e0-u)*(1.0e0-2.0e0*u);
   }
   else {
      return 0.0e0;
   }
#endif

}

/*-------------------------------------------------------------------*/
/////////////////////////// F U N C T I O N ///////////////////////////
/////////////////////////// < P H I _ C T > ///////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 phi_CT(const PS::F64 r, const PS::F64 h) {
   //* Local parameters
   const PS::F64 a=1.0e0/3.0e0;
   const PS::F64 c0=1.0043895747599452e0;
   //* Local variables
   PS::F64 u,u2,cc;

   u=r/h;
   cc=-32.0e0/(c0*h);
   if (u <= a) {
      u2  = u*u;
      return cc*(115.0e0/1296.0e0
                +u2*(-11.0e0/54.0e0
                     +u/6.0e0));
   }
   else if ((a < u) && (u <= 0.5e0)) {
      u2  = u*u;
      return cc*(1.0e0/(7290.0e0*u)
                +7.0e0/80.0e0
                +u2*(-1.0e0/6.0e0
                     +u2*(3.0e-1
                         -0.2e0*u)));
   }
   else if ((0.5e0 < u) && (u <= 1.0e0)) {
      u2  = u*u;
      return cc*(-227.0e0/(1.1664e5*u)
                +0.1e0
                +u2*(-1.0e0/3.0e0
                     +u*(0.5e0
                         +u*(-3.0e-1
                             +u/1.5e1))));
   }
   else {
      return -1.0e0/r;
   }

}

/*-------------------------------------------------------------------*/
////////////////////////    F U N C T I O N    ////////////////////////
//////////////////////// < D P H I D R _ C T > ////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dphidr_CT(const PS::F64 r, const PS::F64 h) {
   //* Local parameters
   const PS::F64 a=1.0e0/3.0e0;
   const PS::F64 c0=1.0043895747599452e0;
   //* Local variables
   PS::F64 u,u2,u3,cc;

   u=r/h;
   cc=32.0e0/(c0*r*r);
   if (u <= a) {
      return cc*u*u*u*(11.0e0/27.0e0-0.5e0*u);
   }
   else if ((a < u) && (u <= 0.5e0)) {
      u2  = u*u; u3 = u2*u;
      return cc*(1.0e0/7.29e3
                +u3*(1.0e0/3.0e0
                    +u2*(-6.0e0/5.0e0
                         +1.0e0*u)));
   }
   else if ((0.5e0 < u) && (u <= 1.0e0)) {
      u2  = u*u; u3 = u2*u;
      return cc*(-2.27e2/1.1664e5
                 +u3*(2.0e0/3.0e0
                     +u*(-1.5e0
                         +u*(6.0e0/5.0e0
                            -u/3.0e0))));
   }
   else {
      return 1.0e0/(r*r);
   }

}

/*-------------------------------------------------------------------*/
////////////////////////    F U N C T I O N    ////////////////////////
//////////////////////// < D P H I D H _ C T > ////////////////////////
/*-------------------------------------------------------------------*/
PS::F64 dphidh_CT(const PS::F64 r, const PS::F64 h) {
   //* Local parameters
   const PS::F64 a=1.0e0/3.0e0;
   const PS::F64 c0=1.0043895747599452e0;
   //* Local variables
   PS::F64 u,u2,cc;

   u=r/h;
   cc=32.0e0/(c0*h*h);
   if (u <= a) {
      u2 = u*u;
      return cc*(1.15e2/1.296e3
                +u2*(-11.0e0/18.0e0
                     +2.0e0*u/3.0e0));
   }
   else if ((a < u) && (u <= 0.5e0)) {
      u2 = u*u;
      return cc*(7.0e0/8.0e1
                +u2*(-0.5e0
                     +u2*(1.5e0
                         -1.2e0*u)));
   }
   else if ((0.5e0 < u) && (u <= 1.0e0)) {
      u2 = u*u;
      return cc*(0.1e0
                +u2*(-1.0e0
                     +u*(2.0e0
                        +u*(-1.5e0
                            +0.4e0*u))));
   }
   else {
      return 0.0e0;
   }

}

/*-------------------------------------------------------------------*/
////////////////          S U B R O U T I N E          ////////////////
//////////////// < C H E C K _ S P H _ K E R N E L S > ////////////////
/*-------------------------------------------------------------------*/
void check_SPH_kernels(void) {
   // This function checks and outputs the profiles of various 
   // SPH kernel functions. If using MPI, you must call this function
   // by RANK 0 only (otherwise, a number of file IO take place).

   //* Local parameters
   const int ndim=128;
   const std::string dir_name="SPH_kernels";
   //* Local variables
   int i,ierr;
   std::string filename,shell_cmd;
   std::ofstream output_file;
   PS::F64 dr,h;
   PS::F64 *r;
   PS::F64 *W,*dWdr,*dWdh,*phi,*dphidr,*dphidh;

   //* Memory allocation 
   r      = new PS::F64[ndim];
   W      = new PS::F64[ndim];
   dWdr   = new PS::F64[ndim];
   dWdh   = new PS::F64[ndim];
   phi    = new PS::F64[ndim];
   dphidr = new PS::F64[ndim];
   dphidh = new PS::F64[ndim];

   //* Define sampling points
   dr = 3.0e0/(double)(ndim);
   for (i=0;i<ndim;i++) r[i] = 0.5e0*dr + dr*(double)i;

   //* Make directory
   shell_cmd = "mkdir -p " + dir_name;
   system(shell_cmd.c_str());

   //* Output M4 Cubic spline kernel
   filename = dir_name + "/W_M4.txt";
   output_file.open(filename.c_str(),std::ios::trunc);
   h=1.825742e0;
   for (i=0;i<ndim;i++) {
      W[i]       = W_M4(r[i],h);
      dWdr[i]    = r[i]*dWdr_M4(r[i],h);
      dWdh[i]    = dWdh_M4(r[i],h);
      phi[i]     = phi_M4(r[i],h);
      dphidr[i]  = dphidr_M4(r[i],h);
      dphidh[i]  = dphidh_M4(r[i],h);
   }
   output_file.setf(std::ios_base::scientific,
                    std::ios_base::floatfield);
   for (i=0;i<ndim;i++) {
      output_file << std::setprecision(15) 
                  << std::showpos
                  << r[i] << " "
                  << W[i] << " "
                  << dWdr[i] << " "
                  << dWdh[i] << " "
                  << phi[i] << " "
                  << dphidr[i] << " "
                  << dphidh[i] << " "
                  << std::endl;
      // Here, we need " "(single space) to separate each value.
   }
   output_file.close();
   // h is set based on Dehnen & Alley (2012)[MNRAS,425,1068].
   // See H/h values in Table 1 of their paper.

   //* Output M5 Quatic spline kernel
   filename = dir_name + "/W_M5.txt";
   output_file.open(filename.c_str(),std::ios::trunc);
   h=2.018932e0;
   for (i=0;i<ndim;i++) {
      W[i]       = W_M5(r[i],h);
      dWdr[i]    = r[i]*dWdr_M5(r[i],h);
      dWdh[i]    = dWdh_M5(r[i],h);
      phi[i]     = phi_M5(r[i],h);
      dphidr[i]  = dphidr_M5(r[i],h);
      dphidh[i]  = dphidh_M5(r[i],h);
   }
   output_file.setf(std::ios_base::scientific,
                    std::ios_base::floatfield);
   for (i=0;i<ndim;i++) {
      output_file << std::setprecision(16) 
                  << std::showpos
                  << r[i] << " "
                  << W[i] << " "
                  << dWdr[i] << " "
                  << dWdh[i] << " "
                  << phi[i] << " "
                  << dphidr[i] << " "
                  << dphidh[i] << " "
                  << std::endl;
   }
   output_file.close();

   //* Output M6 Quintic spline kernel
   filename = dir_name + "/W_M6.txt";
   output_file.open(filename.c_str(),std::ios::trunc);
   h=2.195775e0;
   for (i=0;i<ndim;i++) {
      W[i]       = W_M6(r[i],h);
      dWdr[i]    = r[i]*dWdr_M6(r[i],h);
      dWdh[i]    = dWdh_M6(r[i],h);
      phi[i]     = phi_M6(r[i],h);
      dphidr[i]  = dphidr_M6(r[i],h);
      dphidh[i]  = dphidh_M6(r[i],h);
   }
   output_file.setf(std::ios_base::scientific,
                    std::ios_base::floatfield);
   for (i=0;i<ndim;i++) {
      output_file << std::setprecision(16) 
                  << std::showpos
                  << r[i] << " "
                  << W[i] << " "
                  << dWdr[i] << " "
                  << dWdh[i] << " "
                  << phi[i] << " "
                  << dphidr[i] << " "
                  << dphidh[i] << " "
                  << std::endl;
   }
   output_file.close();

   //* Output CT kernel
   filename = dir_name + "/W_CT.txt";
   output_file.open(filename.c_str(),std::ios::trunc);
   h=1.825742e0;
   for (i=0;i<ndim;i++) {
      W[i]       = W_CT(r[i],h);
      dWdr[i]    = r[i]*dWdr_CT(r[i],h);
      dWdh[i]    = dWdh_CT(r[i],h);
      phi[i]     = phi_CT(r[i],h);
      dphidr[i]  = dphidr_CT(r[i],h);
      dphidh[i]  = dphidh_CT(r[i],h);
   }
   output_file.setf(std::ios_base::scientific,
                    std::ios_base::floatfield);
   for (i=0;i<ndim;i++) {
      output_file << std::setprecision(16) 
                  << std::showpos
                  << r[i] << " "
                  << W[i] << " "
                  << dWdr[i] << " "
                  << dWdh[i] << " "
                  << phi[i] << " "
                  << dphidr[i] << " "
                  << dphidh[i] << " "
                  << std::endl;
   }
   output_file.close();

   //* Output gnuplot script
   filename = dir_name + "/plot_SPH_kernels.plt";
   output_file.open(filename.c_str(),std::ios::trunc);
   output_file << "set terminal postscript enhanced color eps" << std::endl;
   output_file << "set size square" << std::endl;
   output_file << "# [1] W" << std::endl;
   output_file << "set xrange [0.0:2.5]" << std::endl;
   output_file << "set xlabel 'r'" << std::endl;
   output_file << "set ylabel 'W'" << std::endl;
   output_file << "set output 'W.eps'" << std::endl;
   output_file << "plot 'W_M4.txt' u 1:2 w l lw 2 title 'cubic spline', \\" << std::endl;
   output_file << "     'W_M5.txt' u 1:2 w l lw 2 title 'quartic spline', \\" << std::endl;
   output_file << "     'W_M6.txt' u 1:2 w l lw 2 title 'quintic spline', \\" << std::endl;
   output_file << "     'W_CT.txt' u 1:2 w l lw 2 title 'core triangle'" << std::endl;
   output_file.close();

   //* Release memory
   delete[] r,W,dWdr,dWdh,phi,dphidr,dphidh;

}

} // kernel_functions
