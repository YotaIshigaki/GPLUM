/*************************************************************
 *                                                           *
 *  FMM Library for CPPMD  MD simulation program.            *
 *                                                           *
 *  Created by:  Novichkov Gleb, PhD at ERI, In              *
 *  Date:        2011.01 - 2011.05                           *
 *                                                           *
 *  ~~  Class to store pre-computed coefficients needed for  *
 *      FMM computation.  Header file.  ~~                   *
 *                                                           *
 *************************************************************/
#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H

#include <complex>

class Coefficients {

  int expansion_order;

public:
  double * A_mn;                     //  A_mn = (-1)^n/sqrt {(n-m)!(n+m)!}
  double * Y_mn;                     //  Y_mn = sqrt {(n-abs(m))!/(n+abs(m))!}
  std::complex <double> * S_kjmn;    //  S_kjmn = I^{abs(k)-abs(m)-abs(k-m)}*A^m_n*A^{k-m}_{j-n}/A^k_j
  std::complex <double> * T_kjmn;    //  T_kjmn =
  std::complex <double> * R_kjmn;    //  R_kjmn =
  Coefficients ();
  Coefficients (int exp_order);

  void set_expansion_order (int exp_order);
  ~Coefficients ();


  void compute_A_mn ();
  void compute_Y_mn ();
  void compute_STR_kjmn ();


};
#endif  //  #ifndef COEFFICIENTS_H
