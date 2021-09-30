/*************************************************************
 *                                                           *
 *  FMM Library for CPPMD  MD simulation program.            *
 *                                                           *
 *  Created by:  Novichkov Gleb, PhD at ERI, In              *
 *  Date:        2011.01 - 2011.05                           *
 *                                                           *
 *  Multipole Expansion Coefficients Computation Routines.   *
 *                                                           *
 *************************************************************/

#include <stdio.h>        //  for printf()
#include <math.h>         //  for sqrt ()
#include <complex>        //  for complex numbers
#include <stdlib.h>       //  for abs (int)

#include "defs.h"
#include "cell.h"
#include "coefficients.h"


using namespace std;


/*
  Obsolete:  External Data of pre-computed coefficients, now (from 2011.05.26 on) passed through arguments

  extern "C++" Coefficients * coeffs;
  extern "C++" complex <double> * sph_harmonic_storage;

*/

extern "C++" {


//  External Routines


//  Provided by:  math_routines.cpp
  double  pow_m1 (long n);


//  Used by:      convert_MEC2LEC(),
//                translate_MECs(),
//                translate_LECs(),
//                compute_cell_MECs(),
//                LEC_to_potential_force()
//  Provided by:  math_routines.cpp
  void  transform_to_radial(
          const  double position_cartesian[3],      //  [in]    point position, 3D
          const  double centre_cartesian[3],        //  [in]    centre position, 3D
                 double relative_pos_radial[3]      //  [out]   relative radial coordinates
        );
  void  transform_to_radial_pbc(
          const  double position_cartesian[3],      //  [in]    point position, 3D
          const  double centre_cartesian[3],        //  [in]    centre position, 3D
                 double relative_pos_radial[3]      //  [out]   relative radial coordinates
        );

//  Routines for Debugging.
//  Functionality:  checks if the `value` is NAN.  If it is, displays variable name and the name
//                  of the routine that contains it.
//  Provided by:    io_routines.cpp
  void check_value (const char * routine_name, const char * var_name, complex <double> value);
  void check_value (const char * routine_name, const char * var_name, double value);
}



//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Harmonics Computation Routines  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//
//  ~~~~~~~~~~~~~~~~~~~~~~  Spherical Harmonics:  coefficients Y^m_n(theta, phi)/rho^{n+1}  ~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
void    compute_spherical_harmonics (
                  complex  <double>  sph_harmonics [],  //  [out] spherical harmonics (coefficients Y^m_n/rho^{n+1})
                  const     int      expansion_order,   //  [in]  expansion order
                  const     double   coords_sph[3],     //  [in]  {rho, theta, phi}
                  const Coefficients * coeffs           //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
                )
{

  const double & rho   = coords_sph[0];
  const double & theta = coords_sph[1];
  const double & phi   = coords_sph[2];

  complex  <double> H;

  double x      = cos (theta);
  double sqroot = sin (theta); //?sqrt ((1-x)*(1+x));

  double p_mm = 1.0;

  double rho_np1;
  double rho_mp1 = rho;

  for (int m = 0, int_fact=1; m<=expansion_order; m++, int_fact += 2, rho_mp1*=rho) {

    complex <double> expm = exp(complex <double> (0, phi*m ));

    H = coeffs->Y_mn [INDEX(m,m)] * p_mm * expm / rho_mp1;

    sph_harmonics [INDEX(-m,m)] = std::conj (H);
    sph_harmonics [INDEX(m,m)]  = H;

    double p_m_nm1 = 0;                 //  P^m_{n-1}
    double p_m_nm2;                     //  P^m_{n-2}
    double alpha, beta;

    double p_mn = p_mm;
    rho_np1 = rho_mp1;
    for (int n = m+1; n <= expansion_order; n++) {
      rho_np1 *= rho;
      p_m_nm2 = p_m_nm1;
      p_m_nm1 = p_mn;

      //  Updating p_mn
      alpha = x*(((n)<<1) - 1);       //  x*(2*n - 1)
      beta  = double (n + m - 1);
      p_mn = (alpha*p_m_nm1 - beta*p_m_nm2)/(n-m);

      //  Spherical Harmonic coefficients
      H = coeffs->Y_mn [INDEX(m,n)] * p_mn * expm / rho_np1;
      sph_harmonics [INDEX(m,n)]  = H;
      sph_harmonics [INDEX(-m,n)] = std::conj (H);
    }   //  end of n-loop

    //  updating p_mm
    p_mm *= -int_fact*sqroot;
  } //  end of m-loop

} // compute_spherical_harmonics()


//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~  Solid Harmonics:  coefficients Y^m_n(theta, phi)*rho^{n} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
void    compute_solid_harmonics (
                  complex  <double>  sld_harmonics [],  //  [out] solid harmonics (coefficients Y^m_n(\theta, \phi)*rho^n)
                  const     int      expansion_order,   //  [in]  expansion order
                  const     double   coords_sph[3],     //  [in]  {rho, theta, phi}
                  const Coefficients * coeffs           //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
                )
{

  const double & rho   = coords_sph[0];
  const double & theta = coords_sph[1];
  const double & phi   = coords_sph[2];

  complex  <double> H;

  double x      = cos (theta);
  double sqroot = sin (theta); //?sqrt ((1-x)*(1+x));

  double p_mm = 1.0;

  double rho_n;
  double rho_m = 1.0;

  for (int m = 0, int_fact=1; m<=expansion_order; m++, int_fact += 2, rho_m*=rho) {

    complex <double> expm = exp(complex <double> (0, phi*m ));

    H = coeffs->Y_mn [INDEX(m,m)] * p_mm * expm * rho_m;
    sld_harmonics [INDEX(-m,m)] = std::conj (H);
    sld_harmonics [INDEX(m,m)]  = H;

    double p_m_nm1 = 0;                 //  P^m_{n-1}
    double p_m_nm2;                     //  P^m_{n-2}
    double alpha, beta;

    double p_mn = p_mm;
    rho_n = rho_m;
    for (int n = m+1; n <= expansion_order; n++) {
      rho_n *= rho;
      p_m_nm2 = p_m_nm1;
      p_m_nm1 = p_mn;

      //  Updating p_mn
      alpha = x*(((n)<<1) - 1);       // == x*(2*n - 1)
      beta  = double (n + m - 1);
      p_mn = (alpha*p_m_nm1 - beta*p_m_nm2)/(n-m);

      //  Solid Harmonic coefficients
      H = coeffs->Y_mn [INDEX(m,n)] * p_mn * expm * rho_n;
      sld_harmonics [INDEX(m,n)]  = H;
      sld_harmonics [INDEX(-m,n)] = std::conj (H);

    }
    //  updating p_mm
    p_mm *= -int_fact*sqroot;
  }

} //  compute_solid_harmonics()



//
//  Computes solid harmonics multiplied by charge (for MECs computation)
//
void    compute_solid_harmonics_charge (
                  complex  <double>  sld_harmonics [],  //  [in/out] solid harmonics (Y^m_n(\theta, \phi)*rho^n})
                                                        //  NOTE: updated by "+="
                  const     int      expansion_order,   //  [in]  expansion order
                  const     double   coords_sph[3],     //  [in]  {rho, theta, phi}
                  const     double   & charge,          //  [in]  charge
                  const Coefficients * coeffs           //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
                )
{

  const double & rho   = coords_sph[0];
  const double & theta = coords_sph[1];
  const double & phi   = coords_sph[2];

  complex  <double> H;

  double x      = cos (theta);
  double sqroot = sin (theta); //?sqrt ((1-x)*(1+x));

  double p_mm = 1.0;

  double rho_n;
  double rho_m = 1.0;

  double p_m_nm1 = 0;                 //  P^m_{n-1}
  double p_m_nm2;                     //  P^m_{n-2}
  double alpha, beta;
  double p_mn = p_mm;

  // unwinding loop :: m = 0
  {
    sld_harmonics [0] += charge;

    p_m_nm1 = 0;
    p_mn  = p_mm;
    rho_n = rho_m;

    for (int n = 1; n <= expansion_order; n++) {
      rho_n *= rho;
      p_m_nm2 = p_m_nm1;
      p_m_nm1 = p_mn;

      //  Updating p_mn
      alpha = x*(((n)<<1) - 1);       // == x*(2*n - 1)
      beta  = double (n  - 1);
      p_mn = (alpha*p_m_nm1 - beta*p_m_nm2)/n;

      //  Solid Harmonic coefficients
      H = charge * coeffs->Y_mn [INDEX(0,n)] * p_mn * rho_n;
      sld_harmonics [INDEX(0, n)] += H;
    }
    //  updating p_mm
    p_mm  = -sqroot;
    rho_m = rho;
  }

  for (int m = 1, int_fact=3; m<=expansion_order; m++, int_fact += 2, rho_m*=rho) {

    complex <double> expm = exp(complex <double> (0, -phi*m /*orig: phi*m*/));

    H = charge * coeffs->Y_mn [INDEX(m,m)] * p_mm * expm * rho_m;

    //  m>0, so "+=" is OK
    sld_harmonics [INDEX(-m, m)] += std::conj (H);
    sld_harmonics [INDEX( m, m)] += H;

    p_m_nm1 = 0;
    p_mn  = p_mm;
    rho_n = rho_m;

    for (int n = m+1; n <= expansion_order; n++) {
      rho_n *= rho;
      p_m_nm2 = p_m_nm1;
      p_m_nm1 = p_mn;

      //  Updating p_mn
      alpha = x*(((n)<<1) - 1);       // == x*(2*n - 1)
      beta  = double (n + m - 1);
      p_mn = (alpha*p_m_nm1 - beta*p_m_nm2)/(n-m);

      //  Solid Harmonic coefficients
      H = charge * coeffs->Y_mn [INDEX(m,n)] * p_mn * expm * rho_n;

      sld_harmonics [INDEX( m, n)] += H;
      sld_harmonics [INDEX(-m, n)] += std::conj (H);

    }
    //  updating p_mm
    p_mm *= -int_fact*sqroot;
  }

} //  compute_solid_harmonics_charge()


void    compute_solid_harmonics_theta_deriv (
                  complex  <double>  sld_harmonics [],        //  [out] solid harmonics (Y^m_n(\theta, \phi)*rho^{n})
                  complex  <double>  sld_harm_theta_der [],   //  [out] theta derivatives of solid harmonics (d/d\theta [Y^m_n(\theta, \phi)*rho^{n} ])
                  const     int      expansion_order,         //  [in]  expansion order
                  const     double   coords_sph[3],           //  [in]  {rho, theta, phi}
                  const Coefficients * coeffs                 //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
                )
{

  const double & rho   = coords_sph[0];
  const double & theta = coords_sph[1];
  const double & phi   = coords_sph[2];

  complex <double> H, Htd;
  complex <double> LPtd;

  double x      = cos (theta);
  double sqroot = sin (theta); //?sqrt ((1-x)*(1+x));

  double xsqm1 = (x-1)*(x+1);

  double p_mm = 1.0;

  double rho_n;
  double rho_m = 1.0;

  for (int m = 0, int_fact=1; m<=expansion_order; m++, int_fact += 2, rho_m*=rho) {

    complex <double> expm = exp(complex <double> (0, phi*m ));

    //  Solid Harmonics
    H = coeffs->Y_mn [INDEX(m,m)] * p_mm * expm * rho_m;
    sld_harmonics [INDEX(-m,m)] = std::conj (H);
    sld_harmonics [INDEX(m,m)]  = H;

    //  theta derivatives of solid harmonics


    LPtd = -sqroot * x * m*p_mm/xsqm1;      //  Theta derivative of P^m_n(cos(\theta))
    Htd = coeffs->Y_mn [INDEX(m,m)] * LPtd * expm * rho_m;
    sld_harm_theta_der [INDEX(-m,m)] = std::conj (Htd);
    sld_harm_theta_der [INDEX(m,m)]  = Htd;


    double p_m_nm1 = 0;                 //  P^m_{n-1}
    double p_m_nm2;                     //  P^m_{n-2}
    double alpha, beta;

    double p_mn = p_mm;
    rho_n = rho_m;
    for (int n = m+1; n <= expansion_order; n++) {
      rho_n *= rho;
      p_m_nm2 = p_m_nm1;
      p_m_nm1 = p_mn;

      //  Updating p_mn
      alpha = x*(((n)<<1) - 1);       // == x*(2*n - 1)
      beta  = double (n + m - 1);
      p_mn = (alpha*p_m_nm1 - beta*p_m_nm2)/(n-m);

      //  Solid Harmonic coefficients
      H = coeffs->Y_mn [INDEX(m,n)] * p_mn * expm * rho_n;
      sld_harmonics [INDEX(-m,n)] = std::conj (H);
      sld_harmonics [INDEX(m,n)]  = H;

      //  Theta Derivative of Solid Harmonics
      LPtd = -sqroot * (n*x*p_mn - (n+m)*p_m_nm1)/ xsqm1;       //  Theta derivative of P^m_n(cos(\theta))
      Htd = coeffs->Y_mn [INDEX(m,n)] * LPtd * expm * rho_n;
      sld_harm_theta_der [INDEX(-m,n)] = std::conj (Htd);
      sld_harm_theta_der [INDEX(m,n)]  = Htd;

    }
    //  updating p_mm
    p_mm *= -int_fact*sqroot;
  }

} //  compute_solid_harmonics_theta_deriv()




//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END of Harmonics Computation Routines ~~~~~~~~~~~~~~~~~~~~~~~~~
//


//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Translation Routines  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
#include <omp.h>
//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Translation:  MECs -> LECs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
void  convert_MEC2LEC (
        complex <double>  new_LEC [],     //  [out] LEC coefficients
        complex <double>  old_MEC [],     //  [in]  MEC coefficients
        const double   LEC_centre [3],    //  [in]  position of the center local expansion
        const double   MEC_centre [3],    //  [in]  position of the center multipole expansion
        const int &    p,                 //  [in]  expansion order
        const Coefficients * coeffs,      //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
        complex <double> * sph_harmonic_storage   //  [in]  pointer to the pre-allocated buffer to store Harmonic coefficients
  )
{
  double    MEC_centre_sph [3];
  //  double &  rho    = MEC_centre_sph[0];
  //  double &  alpha  = MEC_centre_sph[1];
  //  double &  beta   = MEC_centre_sph[2];

  transform_to_radial (MEC_centre, LEC_centre, MEC_centre_sph);

  complex <double> * H_mn = sph_harmonic_storage;
  // If dynamic allocation is needed:  new complex <double> [(2*p+1)*(2*p+1)];

  compute_spherical_harmonics (H_mn, 2*p, MEC_centre_sph /*{rho, theta, phi}*/, coeffs);

//  THIS 4-NESTED LOOP TAKES ABOUT 95% of all computational time


  complex<double>  L_kj;


/*
  #pragma	omp parallel \
			    shared (new_LEC, old_MEC, H_mn, coeffs) \
				  private (j, k, m, n, sum, L_kj) \
				  default(none)
*/

/*  Empty `#pragma	omp parallel' does not stall */
#if 0
  #pragma	omp parallel
  {
  	
  }
#endif

/*  This codelets causes stall: */
#if 0
  #pragma	omp parallel 
  {
  	#pragma omp for
  	for (int i=0; i<10; i++);
  }
#endif

/*  This codelets causes stall: */
#if 0
  #pragma	omp parallel default(none)
  {
  	#pragma omp for
  	for (int i=0; i<10; i++);
  }
#endif

#if 0
  #pragma	omp parallel default(none)
  {
  	//printf ("OMP:: num_threads = %d\n", omp_get_num_threads ()); 
  	for (int i=0; i<10; i++);
  }
#endif


#if 0
  #pragma	omp parallel default(none)
  {
  	//printf ("OMP:: tid = %d\n", omp_get_thread_num ()); 
  	for (int i=0; i<10; i++);
  	//printf ("OMP:: tid = %d\n", omp_get_thread_num ()); 
  }
#endif



#pragma	omp parallel shared(new_LEC, old_MEC, H_mn, coeffs)
  {

		// int nthreads = omp_get_num_threads();	
		// int tid = omp_get_thread_num();	

		#pragma omp for 
    for (int j = 0; j<=p; j++) {  /* TODO:  correct this:		for (j=tid; j<p; j+=nthreads) {	*/
      for (int k = 0; k<=j; k++) {
			  complex<double>  sum = 0;

        //  ~~~~~~~~~~~~~~~ Loop computing coefficients L^k_j (cf. 5.5 Greengard-Rohlin's paper in Acta Numerica (1997), pp 229-269)
        for (int n = 0; n<=p; n++) {
          for (int m = -n; m<=n; m++) {
            sum += old_MEC [INDEX(m,n)] * coeffs->T_kjmn [SUPERINDEX(k,j,m,n,p)] * H_mn [INDEX(m-k, j+n)];
          } // m-loop
        } // n-loop


        new_LEC [INDEX(-k,j)] = std::conj (sum/*L_kj*/);
        new_LEC [INDEX(k,j)]  = sum; /*L_kj;*/

      }   //  k-loop
    }   //  j-loop


  }  // end of OMP PRAGMA




  /*  If H_mn was dynamically allocated:  delete [] H_mn; */

} //  convert_MEC2LEC()

void  convert_MEC2LEC_pbc (
        complex <double>  new_LEC [],     //  [out] LEC coefficients
        complex <double>  old_MEC [],     //  [in]  MEC coefficients
        const double   LEC_centre [3],    //  [in]  position of the center local expansion
        const double   MEC_centre [3],    //  [in]  position of the center multipole expansion
        const int &    p,                 //  [in]  expansion order
        const Coefficients * coeffs,      //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
        complex <double> * sph_harmonic_storage   //  [in]  pointer to the pre-allocated buffer to store Harmonic coefficients
  )
{
  double    MEC_centre_sph [3];
  //  double &  rho    = MEC_centre_sph[0];
  //  double &  alpha  = MEC_centre_sph[1];
  //  double &  beta   = MEC_centre_sph[2];

  transform_to_radial_pbc (MEC_centre, LEC_centre, MEC_centre_sph);

  complex <double> * H_mn = sph_harmonic_storage;
  // If dynamic allocation is needed:  new complex <double> [(2*p+1)*(2*p+1)];

  compute_spherical_harmonics (H_mn, 2*p, MEC_centre_sph /*{rho, theta, phi}*/, coeffs);

//  THIS 4-NESTED LOOP TAKES ABOUT 95% of all computational time


  complex<double>  L_kj;

#pragma	omp parallel shared(new_LEC, old_MEC, H_mn, coeffs)
  {

		// int nthreads = omp_get_num_threads();	
		// int tid = omp_get_thread_num();	

		#pragma omp for 
    for (int j = 0; j<=p; j++) {  /* TODO:  correct this:		for (j=tid; j<p; j+=nthreads) {	*/
      for (int k = 0; k<=j; k++) {
			  complex<double>  sum = 0;

        //  ~~~~~~~~~~~~~~~ Loop computing coefficients L^k_j (cf. 5.5 Greengard-Rohlin's paper in Acta Numerica (1997), pp 229-269)
        for (int n = 0; n<=p; n++) {
          for (int m = -n; m<=n; m++) {
            sum += old_MEC [INDEX(m,n)] * coeffs->T_kjmn [SUPERINDEX(k,j,m,n,p)] * H_mn [INDEX(m-k, j+n)];
          } // m-loop
        } // n-loop


        new_LEC [INDEX(-k,j)] = std::conj (sum/*L_kj*/);
        new_LEC [INDEX(k,j)]  = sum; /*L_kj;*/

      }   //  k-loop
    }   //  j-loop


  }  // end of OMP PRAGMA




  /*  If H_mn was dynamically allocated:  delete [] H_mn; */

} //  convert_MEC2LEC_pbc()


//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Translation:  MECs -> MECs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
void
translate_MECs (
        complex <double>  new_MEC [],     //  [out]  new multipole expansion coefficients
        complex <double>  old_MEC [],     //  [in]  old multipole expansion coefficients
        const double   new_centre [3],    //  [in]  position of the new center of the expansion
        const double   old_centre [3],    //  [in]  position of the center of the original expansion
        const int &    p,                 //  [in]  expansion order
        const Coefficients * coeffs       //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
  )
{

  double    MEC_centre_sph [3];
  //  double &  rho    = MEC_centre_sph[0];
  //  double &  alpha  = MEC_centre_sph[1];
  //  double &  beta   = MEC_centre_sph[2];

  transform_to_radial (old_centre, new_centre, MEC_centre_sph);

  complex <double> * H_mn = new complex <double> [(2*p+1)*(2*p+1)];

  // complex <double> * H_mn = sph_harmonic_storage;
  // If dynamic allocation is needed:  new complex <double> [(2*p+1)*(2*p+1)];

  compute_solid_harmonics (H_mn, 2*p, MEC_centre_sph /*{rho, theta, phi}*/, coeffs);

  int j, k;
  complex<double>  sum;

  for (j = 0; j<=p; j++) {
    for (k = 0; k<=j; k++) {

      sum = 0;

      for (int n = 0; n<=j; n++){
        for (int m = max(n-j+k,-n); m<=min(j-n+k,n); m++){
        /* Equivalent to  for (int m = -n; m<=n; m++){    if (abs (k-m) > abs (j-n)) continue; ...       //  suggested by the L. Greengard.           ...           }
         */
          sum += old_MEC[INDEX (k-m, j-n)] * coeffs->S_kjmn [SUPERINDEX(k,j,m,n,p)] * H_mn [INDEX(-m, n)];;
        }
      }

      new_MEC [INDEX(-k,j)] = std::conj (sum);
      new_MEC [INDEX(k,j)]  = sum;

    } //  end of k-loop
  } //  end of j-loop
  delete [] H_mn;
}   // end of translate_MECs()

//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Translation:  LECs -> LECs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
void
translate_LECs (
        complex <double>  new_LEC [],     //  [out]  new multipole expansion coefficients.  direct writing, NO "+="
        complex <double>  old_LEC [],     //  [in]  old multipole expansion coefficients
        const double   new_centre [3],    //  [in]  position of the new center of the expansion
        const double   old_centre [3],    //  [in]  position of the center of the original expansion
        const int &    p,                 //  [in]  expansion order
        const Coefficients * coeffs       //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
  )
{

  double    LEC_centre_sph [3];
  //  double &  rho    = LEC_centre_sph[0];
  //  double &  alpha  = LEC_centre_sph[1];
  //  double &  beta   = LEC_centre_sph[2];

  transform_to_radial (old_centre, new_centre, LEC_centre_sph);

  complex <double> * H_mn = new complex <double> [(2*p+1)*(2*p+1)];

  // complex <double> * H_mn = sph_harmonic_storage;
  // If dynamic allocation is needed:  new complex <double> [(2*p+1)*(2*p+1)];

  compute_solid_harmonics (H_mn, 2*p, LEC_centre_sph /*{rho, theta, phi}*/, coeffs);

  int j, k;
  complex<double>  sum, Y;

  for (j = 0; j<=p; j++) {
    for (k = 0; k<=j; k++) {

      sum = 0;

      for (int n = j; n<=p; n++){
//        for (int m = max(n-j+k,-n); m<=min(j-n+k,n); m++){
        for (int m = -n; m<=n; m++){
          if (abs (m-k) > abs (n-j)) continue;
        /* Equivalent to  for (int m = -n; m<=n; m++){    if (abs (k-m) > abs (j-n)) continue; ...       //  suggested by the L. Greengard.           ...           }
         */
          sum += old_LEC[INDEX (m, n)] * coeffs->R_kjmn [SUPERINDEX(k,j,m,n,p)] * H_mn [INDEX(m-k, n-j)];
        }
      }

      new_LEC [INDEX(-k,j)] = std::conj (sum);
      new_LEC [INDEX(k,j)]  = sum;

    } //  end of k-loop
  } //  end of j-loop
  delete [] H_mn;
}   // end of translate_LECs()




//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Evaluation Routines  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Computes MECs of a Cell  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
void
compute_cell_MECs (
          Cell & cell,                      //  [in/out]  Warning: cell MECs are updated by "+=".  Therefore check if clearing beforehand is necessary
          const Coefficients * coeffs       //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
    )
{

  int  p = Cell::get_expansion_order();                   //  expansion order p
  unsigned long  n_charges = cell.particles.size();       //  number of particles (charges) in the cell

  double *  & charges   = cell.charges;                   //  array of charges in the cell
  double *  & positions = cell.positions;                 //  positions of the charges, array of triples (x_i, y_i, z_i)

  complex <double> * & MEC = cell.MEC;                    //  reference to the MECs of the cell

  double coord_sph [3];             //  spherical coordinates ( rho, theta, phi );

  for (unsigned long i = 0; i<n_charges; i++) {
    transform_to_radial (positions + 3*i, cell.centre, coord_sph);
    compute_solid_harmonics_charge (MEC, p, coord_sph, charges[i], coeffs);
  }

}




//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Evaluation: LECs -> Potential/Force  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
void
LEC_to_potential_force (
        double & potential,             //  [in/out] potential (updated by "+=" )
        double field[3],                //  [out] field
        double position [3],            //  [in]  Particle's position
        double center   [3],            //  [in]  Center of the cell
        const complex<double> LEC[],    //  [in]  LECs
        const int & expansion_order,    //  [in]  order of expansion
        const Coefficients * coeffs     //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
      )
{

  double  rel_position_sph [3];

  double &  rho    = rel_position_sph [0];
  double &  theta  = rel_position_sph [1];
  double &  phi    = rel_position_sph [2];

  double  rho_n;

  complex<double>  PHI, PHI_r, PHI_phi, PHI_theta;
  complex<double>  YL_mn, LY_theta_mn, sum_YL_mn, sum_mYL_mn, sum_LY_theta_mn;

  complex<double>  * H_mn   = new complex<double> [Cell::get_num_coeff()];
  complex<double>  * Htd_mn = new complex<double> [Cell::get_num_coeff()];


  transform_to_radial (position, center, rel_position_sph);
  compute_solid_harmonics_theta_deriv (H_mn, Htd_mn, expansion_order, rel_position_sph, coeffs);

  PHI = 0;
  PHI_r = 0;
  PHI_phi = 0;
  PHI_theta = 0;

  rho_n = 1.0;
  for (int n = 0; n <= expansion_order; n++, rho_n*=rho){

    sum_YL_mn  = sum_mYL_mn = sum_LY_theta_mn = 0;

    for (int m = -n; m<=n; m++) {
      YL_mn             = H_mn[ INDEX(m,n) ] * LEC [ INDEX(m,n) ];
      sum_LY_theta_mn  += Htd_mn[INDEX (m, n)]*LEC [ INDEX(m,n) ];;
      sum_YL_mn        += YL_mn;
      sum_mYL_mn       += double(m)*YL_mn;
    }  // end of m-loop

    PHI     += sum_YL_mn;
    PHI_r   += double(n) * sum_YL_mn / rho;
    PHI_phi += sum_mYL_mn;
    PHI_theta += sum_LY_theta_mn;

  } // end of n-loop

  PHI_phi *= std::complex<double> (0, 1);

  potential = real (PHI);

  double  PHI_x, PHI_y, PHI_z;

  double  RPHI_r     = real (PHI_r);
  double  RPHI_theta = real (PHI_theta);
  double  RPHI_phi   = real (PHI_phi);


  PHI_x = cos(phi)*sin(theta)*RPHI_r + cos(phi)*cos(theta)/rho*RPHI_theta - sin(phi)/sin(theta)/rho*RPHI_phi;
  PHI_y = sin(phi)*sin(theta)*RPHI_r + sin(phi)*cos(theta)/rho*RPHI_theta + cos(phi)/sin(theta)/rho*RPHI_phi;
  PHI_z =          cos(theta)*RPHI_r -          sin(theta)/rho*RPHI_theta;

  field [0] = -PHI_x;
  field [1] = -PHI_y;
  field [2] = -PHI_z;

  delete [] H_mn;
  delete [] Htd_mn;
}





