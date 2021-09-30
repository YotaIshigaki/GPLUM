/*************************************************************
 *                                                           *
 *  FMM Library for CPPMD  MD simulation program.            *
 *                                                           *
 *  Created by:  Novichkov Gleb, PhD at ERI, In              *
 *  Date:        2011.01 - 2011.05                           *
 *                                                           *
 *  ~~ Class to store pre-computed coefficients needed for   *
 *     FMM computation.  Source file.  ~~                    *
 *                                                           *
 *************************************************************/

#include "defs.h"
#include "coefficients.h"
#include <cstdlib>              //  for abs (int &)
#include <iostream>

using namespace std;

extern "C++"{

//  provided by expansion_routines.cpp
  double compute_A (
    long m,       //  top index
    long n        //  bottom index
  );

  complex <double> pow_I (long m);

  double pow_m1 (long m);

}

Coefficients :: Coefficients ()
{
  expansion_order = 0;
  A_mn = Y_mn = NULL;
  S_kjmn = T_kjmn = R_kjmn = NULL;
}

Coefficients :: Coefficients (int exp_order) : expansion_order (exp_order) {

  A_mn = new double [(4*expansion_order+1)*(4*expansion_order+1)];
  Y_mn = new double [(4*expansion_order+1)*(4*expansion_order+1)];
  int & p = expansion_order;

  S_kjmn = new complex <double> [(p+1)*(p+1)*(p+1)*(p+1)];
  T_kjmn = new complex <double> [(p+1)*(p+1)*(p+1)*(p+1)];
  R_kjmn = new complex <double> [(p+1)*(p+1)*(p+1)*(p+1)];

  compute_A_mn ();
  compute_Y_mn ();
  compute_STR_kjmn ();

}

void Coefficients :: set_expansion_order (int exp_order) {  expansion_order = exp_order; }

void Coefficients :: compute_A_mn ()
{
  for (int n = 0; n<=2*expansion_order; n++)
  for (int m = -n; m<=n; m++) {
      A_mn [INDEX(m,n)] = compute_A (m,n);
  }
}

void Coefficients :: compute_Y_mn ()
{
  for (int n = 0; n<=2*expansion_order; n++)
  for (int m = -n; m<=n; m++) {

      double part_fact = 1.0;
      for (int i = n - abs(m) + 1; i <= n + abs(m); i++)
        part_fact *= i;
      Y_mn [INDEX(m,n)]= 1.0/sqrt (part_fact);
  }

}

extern "C++" void check_value (const char * routine_name, const char * var_name, std::complex <double> value);
#include <assert.h>
void Coefficients :: compute_STR_kjmn ()
{
  int & p = expansion_order;
  for (int j=0; j<=p; j++)
  for (int k=-j; k<=j; k++)
    for (int n=0; n<=p; n++)
    for (int m=-n; m<=n; m++) {
      double & Amn      = A_mn[ INDEX(m,n) ];
      double & Akj      = A_mn[ INDEX(k,j) ];

      if (abs (k-m) <= abs (j-n)) {
        S_kjmn [SUPERINDEX(k,j,m,n,p)] = pow_I( abs(k) - abs(m) - abs(k-m)) * Amn * A_mn[ INDEX(k-m, j-n) ] / Akj;
        R_kjmn [SUPERINDEX(k,j,m,n,p)] = pow_I( abs(m) - abs(m-k) - abs(k)) * A_mn[ INDEX(m-k, n-j) ]* Akj * pow_m1(n+j) / Amn ;

#ifdef FMM_CHECKVALUES
        check_value (__FILE__, __FUNCTION__, S_kjmn [SUPERINDEX(k,j,m,n,p)]);
        check_value (__FILE__, __FUNCTION__, R_kjmn [SUPERINDEX(k,j,m,n,p)]);
        assert (S_kjmn [SUPERINDEX(k,j,m,n,p)]==S_kjmn [SUPERINDEX(k,j,m,n,p)]);
        assert (R_kjmn [SUPERINDEX(k,j,m,n,p)]==R_kjmn [SUPERINDEX(k,j,m,n,p)]);
#endif
      }
/*      if (abs (m-k) <= j + n) :: This condition is void, as -j<=k<=j, -n<=m<=n, a) m-k<=m+j<=n+j, and b) k-m <= j+m <= j-m <= j+n */
      T_kjmn [SUPERINDEX(k,j,m,n,p)] = pow_I( abs(k-m) - abs(k) - abs(m)) * Amn * Akj * pow_m1(n) / A_mn[ INDEX(m-k, j+n) ];

#ifdef FMM_CHECKVALUES
      check_value (__FILE__, __FUNCTION__, T_kjmn [SUPERINDEX(k,j,m,n,p)]);
      assert (T_kjmn [SUPERINDEX(k,j,m,n,p)]==T_kjmn [SUPERINDEX(k,j,m,n,p)]);
#endif
    }

}


Coefficients :: ~Coefficients()
{
  if (A_mn)   {  delete [] A_mn;    A_mn = NULL;  }
  if (Y_mn)   {  delete [] Y_mn;    Y_mn = NULL;  }
  if (S_kjmn) {  delete [] S_kjmn;  S_kjmn = NULL;  }
  if (T_kjmn) {  delete [] T_kjmn;  T_kjmn = NULL;  }
  if (R_kjmn) {  delete [] R_kjmn;  R_kjmn = NULL;  }
}
