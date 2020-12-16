/*************************************************************
 *                                                           *
 *  FMM Library for CPPMD  MD simulation program.            *
 *                                                           *
 *  Created by:  Novichkov Gleb, PhD at ERI, Inc             *
 *               Ohno Yousuke, PhD at RIKEN                  *
 *  Date:        2011.01 - 2011.06                           *
 *                                                           *
 *            ~~~ CommMEC class source file  ~~~             *
 *                                                           *
 *************************************************************/
#include "CommMEC.h"
#include <stdio.h>              // for printf
#include <string.h>             // for memcpy
#include <stdlib.h>             // for realloc()


int CommMEC :: num_coeff = 0;


CommMEC :: CommMEC()
{
  N=0;
  data_size = 0;
  data = NULL;
}

void CommMEC :: add (const ucellindex_t & cell_idx, const std::complex<double> * MECs)
{
  if (num_coeff == 0) {
    printf ("CommMEC :: num_coeff = %d, add () is not performed.\n", num_coeff);
    return;
  }

  int additional_size = sizeof (ucellindex_t) + sizeof (std::complex<double>)*num_coeff;

  char * p = (char *) realloc (data, data_size + additional_size);
  if (p) {
    data=p;

    memcpy (data+data_size,                       &cell_idx, sizeof(ucellindex_t));
    memcpy (data+data_size+sizeof (ucellindex_t), MECs,      sizeof(std::complex<double>)*num_coeff);

    data_size+=additional_size;

    N++;
  } else {
    printf ("CommMEC :: realloc() failed. returned value is %p\n", p);
  }
}

void CommMEC :: allocate(const int & size)
{
  data = (char *) malloc(size);
  data_size = size;
}


void CommMEC :: allocate()
{
	if (N>0) {
	  int size = N*(sizeof (ucellindex_t) + sizeof (std::complex<double>)*num_coeff);
		allocate (size);
	}
	else {
#ifdef FMM_VERBOSE	
    printf ("CommMEC :: N = %d, so allocate () is not performed.\n", N);
#endif
  }
}

void  CommMEC :: empty () 
{
  N = data_size = 0;
  if (data != NULL) { free (data);  data = NULL; }
}


ucellindex_t  CommMEC :: get_idx (const int & i)
{
  int size = sizeof (ucellindex_t) + sizeof (std::complex<double>)*num_coeff;
  char * p = data + size*i;

  ucellindex_t  idx;

  memcpy (&idx, p, sizeof (ucellindex_t));
  return idx;
}

std::complex<double> * CommMEC :: get_MECs (const int & i)
{
  int size = sizeof (ucellindex_t) + sizeof (std::complex<double>)*num_coeff;
  char * p = data + size*i + sizeof (ucellindex_t);
  return (std::complex<double> *) p;
}

CommMEC :: ~CommMEC()
{
  empty();
}


void  CommMEC :: set_num_coeff (const int &_num_coeff) { CommMEC :: num_coeff = _num_coeff; }
