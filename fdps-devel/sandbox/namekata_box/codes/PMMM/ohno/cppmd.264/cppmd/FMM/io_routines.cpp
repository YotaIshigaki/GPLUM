/*************************************************************
 *                                                           *
 *  FMM Library for CPPMD  MD simulation program.            *
 *                                                           *
 *  Created by:  Novichkov Gleb, PhD at ERI, In              *
 *  Date:        2011.01 - 2011.05                           *
 *                                                           *
 *                ~~~  I/O Routines  ~~~                     *
 *                                                           *
 *************************************************************/
#include <stdio.h>        //  for printf()
#include <vector>         //  for vector class
#include <complex>        //  for complex class
#include "cell.h"

using namespace std;

void print_list (vector<ucellindex_t> list) {
  for (int i = 0; i<list.size(); i++)
    printf ("%ld\n",(long) list [i]);
//  printf ("Total %d elements\n", list.size());
}

void print_point (double point[], int dimension) {
  for (int i = 0; i<dimension; i++)
    printf ("%.3f\n", point [i]);
}

void print_point (float point[], int dimension) {
  for (int i = 0; i<dimension; i++)
    printf ("%.3f\n", point [i]);
}
