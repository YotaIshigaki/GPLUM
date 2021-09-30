/*************************************************************
 *                                                           *
 *  FMM Library for CPPMD  MD simulation program.            *
 *                                                           *
 *  Created by:  Novichkov Gleb, PhD at ERI, Inc             *
 *               Ohno Yousuke, PhD at RIKEN                  *
 *  Date:        2011.01 - 2011.06                           *
 *                                                           *
 *            ~~~ CommMEC class header file  ~~~             *
 *                                                           *
 *************************************************************/
#ifndef CommMEC_H
#define CommMEC_H

#include "defs.h"                       //  for ucellindex_t
#include <complex>                      //  for MECs

class CommMEC {



public:
  char * data;
  int N;
  int data_size;
  static int num_coeff;


  CommMEC();
  ~CommMEC();

  void  allocate(const int & size);
  void  allocate();
  
  void  empty ();
  
  void  add (const ucellindex_t & cell_idx, const std::complex<double> * MECs);

	ucellindex_t  get_idx (const int & i);

	std::complex<double> * get_MECs (const int & i);

  static void set_num_coeff (const int &);
};
#endif
