/*************************************************************
 *                                                           *
 *  FMM Library for CPPMD  MD simulation program.            *
 *                                                           *
 *  Created by:  Novichkov Gleb, PhD at ERI, In              *
 *  Date:        2011.01 - 2011.05                           *
 *                                                           *
 *                  ~~~  Math Routines  ~~~                  *
 *                                                           *
 *************************************************************/
#include  <stdio.h>        //  for printf()
#include  <math.h>         //  for sqrt ()
#include  <complex>        //  for complex<double> type
#include  <stdlib.h>       //  for abs()
#include  "coefficients.h"

using namespace std;


extern "C++" {
   //  routines for debugging.  Provided by:  io_routines.cpp
  void check_value (const char * routine_name, const char * var_name, complex <double> value);
  void check_value (const char * routine_name, const char * var_name, double value);
}

//
//  ======== Factorial Routines ========
//

//  ======== Computes ln (Gamma (x)) ========
double ln_gamma (double x)
{
  static double c_coef[6] = {
     76.18009172947146,
    -86.50532032941677,
     24.01409824083091,
    -1.231739572450155,
     0.1208650973866179E-2,
    -0.5395239384953e-5 };

#define PI    3.141592653589793238462643383279502884197169399375
#define coef1 (sqrt(2*PI))

  double x1, x2, y, z;
  double ln_val;

  x1 = x + .5;
  x2 = x1 + 5;  // x + .5 + 5

  y = x1*log(x2) - x2;

  double series = 1.000000000190015;

  z = x;
  for (int i = 0; i<6; i++)
    series += c_coef[i]/(++z);

  ln_val = y + log(coef1*series/x);
  return ln_val;

}

//  ======== Computes Gamma-function Gamma(x) ========
double gamma_function (double x) {
  return exp (ln_gamma(x));
}

//  ======== Computes n! ========
double factorial (int n)
{

  double fact;
  static double tabulated[33]={1, 1, 2, 6, 24, 120, 720, 5040};
  tabulated[28] = 304888344611713860501504000000.0;
  tabulated[29] = 8841761993739701954543616000000.0;
  tabulated[30] = 265252859812191058636308480000000.0;
  tabulated[31] = 8222838654177922817725562880000000.0;
  tabulated[32] = 263130836933693530167218012160000000.0;


  static int top = 7;

  if (n<0) {
    printf ("factorial() :: Negative parameter: n=%d\n", n);
    return -1.0;
  }

  if (n>=33)
    return gamma_function(n+1);


  if (28<=n && n<33)
    return  tabulated[n];


  fact = 1.0;
  int i;
  while (top<n) {
    i = top + 1;
    tabulated[i] = double(i)*tabulated[top];
    top++;
  }
  return tabulated[n];
}


//
//  ======== Power Computation Routines ========
//


//  ======== Computes x^n, n>=0. ========

double my_pow (double x, int n)
{
#define BREAK_EVEN_POW 23           // obtained from benchmarking

  double x_n = 1;

  if (n<0 || n>BREAK_EVEN_POW) {
    x_n = pow(x, n);
  }
  else {
//  TO DO:  replace this by devide/conquere algorithm
    for (int j = 0; j<n; j++)
      x_n *=x;
  }
  return x_n;
}





//  ======== Computes i^n ========
std :: complex<double> pow_I (long m)
{
  // this static cause error Intel icpc -static
  //  static complex<double> pows [4] = {
  complex<double> pows [4] = {
      complex<double> (1,0),
      complex<double> (0,1),
      complex<double> (-1,0),
      complex<double> (0,-1)
    };

  return pows [m & 3];
/*
  // The (original) code below is slower compared to reading from memory.
  // Hardware:  Intel(R) Core(TM)2 Duo CPU     E8500  @ 3.16GHz.

  switch (m & 3){  // % 4 gives signed values ...
    case  0:  return std :: complex<double> (1,0);
    case  1:  return std :: complex<double> (0,1);
    case  2:  return std :: complex<double> (-1,0);
    case  3:  return std :: complex<double> (0,-1);
  }
*/
}

//  ======== Computes (-1)^n ========
double pow_m1(long m)
{
  double pows [2] = {1.0, -1.0};
  return pows [m & 1];

/*
  // The (original) code below is slower compared to reading from memory.
  // Hardware:  Intel(R) Core(TM)2 Duo CPU     E8500  @ 3.16GHz.

  return (m & 1) ? -1.0 : 1.0;  // % 2 gives signed values ...
*/

}


//
//  Coordinate Transformation Routines
//

//  ====== Computes radial coordinates of a point position relative to a centre point
void  transform_to_radial(
        const  double position_cartesian[3],      //  [in]    point position, 3D
        const  double centre_cartesian[3],        //  [in]    centre position, 3D
               double relative_pos_radial[3]      //  [out]   relative position, in radial coordinates
      )
{
#define TWO_PI    (2*3.141592653589793238462643383279502884197169399375)
  double delta[3];
  double rho, phi, theta, sin_th;

  phi = theta = rho = 0;

//  ====== computing dx, dy, dz and distance to the centre (rho)
  for (int i = 0; i<3; i++) {
    delta [i] = position_cartesian[i] - centre_cartesian[i];
    rho += delta [i] * delta [i];
  }
  rho = sqrt (rho);

//  ====== computing angles, if not in the center
  if (rho>0) {
    theta  = acos (delta [2]/rho);            // -- rho*cos (theta) = z;
    phi    = atan2 (delta [1], delta[0]);

    if (phi<0)
      phi += TWO_PI;
  }
//  ====== returning relative position, in radial coordinates
  relative_pos_radial [0] = rho;
  relative_pos_radial [1] = theta;
  relative_pos_radial [2] = phi;
}

void  transform_to_radial_pbc(
        const  double position_cartesian[3],      //  [in]    point position, 3D
        const  double centre_cartesian[3],        //  [in]    centre position, 3D
               double relative_pos_radial[3]      //  [out]   relative position, in radial coordinates
      )
{
#define TWO_PI    (2*3.141592653589793238462643383279502884197169399375)
  double delta[3];
  double rho, phi, theta, sin_th;

  phi = theta = rho = 0;

//  ====== computing dx, dy, dz and distance to the centre (rho)
  for (int i = 0; i<3; i++) {
    delta [i] = position_cartesian[i] - centre_cartesian[i];
    if(delta[i]>0.5){
      printf("shift %d %f - 1.0\n",i,delta[i]);
      delta[i]-=1.0;
    }
    if(delta[i]<-0.5){
      printf("shift %d %f + 1.0\n",i,delta[i]);
      delta[i]+=1.0;
    }
    rho += delta [i] * delta [i];
  }
  rho = sqrt (rho);

//  ====== computing angles, if not in the center
  if (rho>0) {
    theta  = acos (delta [2]/rho);            // -- rho*cos (theta) = z;
    phi    = atan2 (delta [1], delta[0]);

    if (phi<0)
      phi += TWO_PI;
  }
//  ====== returning relative position, in radial coordinates
  relative_pos_radial [0] = rho;
  relative_pos_radial [1] = theta;
  relative_pos_radial [2] = phi;
}


//  ====== Computes coefficients A^m_n ======
double
compute_A (
    long m,       //  top index
    long n        //  bottom index
  )
{
  double  x1 = factorial (n-m);
  double  x2 = factorial (n+m);
  return pow_m1 (n)/sqrt (x1*x2);
}
