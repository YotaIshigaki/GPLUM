/*************************************************************
 *                                                           *
 *  FMM Library for CPPMD  MD simulation program.            *
 *                                                           *
 *  Created by:  Novichkov Gleb, PhD at ERI, Inc.            *
 *  Date:        2011.01 - 2011.05                           *
 *                                                           *
 *          ~~~  Cell class.  Source file.  ~~~              *
 *                                                           *
 *************************************************************/

#include "defs.h"                       //  FMM definitions file
#include "cell.h"                       //  Cell class header file

#include <complex>                      //  for MEC, LEC, etc
#include <vector>                       //  for vector class
#include <stdio.h>                      //  for fprintf(), etc

#include <string.h>                     //  for memset()


int Cell :: expansion_order;
int Cell :: num_coeff;

//
//                                          Methods related to Expansion Order, Number of Coefficients
//

void Cell :: set_expansion_order (int n)
{
  expansion_order = n;
  num_coeff = (n+1)*(n+1); /*  After optimization must be: n*(n+1)/2  */;
}

int  Cell :: get_expansion_order () {
  return expansion_order;
}

int Cell :: get_num_coeff () {
  return num_coeff;
}


//
//                                          Memory allocation methods
//
void Cell :: allocate()
{
  MEC = new std::complex<double> [num_coeff];
  LEC = new std::complex<double> [num_coeff];
  for (int i = 0; i<num_coeff; i++) {  MEC[i] = LEC[i] = 0; }
}

void Cell :: deallocate()
{

  destroy_particles_data();

  if (MEC) {
    delete [] MEC;    MEC = NULL;
  }
  if (LEC) {
    delete [] LEC;    LEC = NULL;
  }
}


//
//                                          Constructors/Destructors
//
Cell :: Cell()
{
  centre [0] = centre [1] = centre [2] = 0;

  positions = NULL;
  charges   = NULL;

  if (num_coeff) {
    allocate();
  }
  else {
    MEC = LEC = NULL;
  }
}

Cell ::  ~Cell()
{
  deallocate();
}


//
//                                          In-Cell Particles Related Methods
//
void  Cell :: clear_mecs() {
  if (MEC)
    memset (MEC, 0, num_coeff*sizeof (std::complex<double>));
}

void  Cell :: clear_lecs() {
  if (LEC)
    memset (LEC, 0, num_coeff*sizeof (std::complex<double>));
}

void  Cell :: destroy_particles_data()
{
  if (positions) {
    delete[] positions; positions = NULL;
  }

  if (charges) {
    delete[] charges; charges = NULL;
  }
}


void Cell :: allocate_positions_charges()
{
  unsigned long size = particles.size();
  positions = new double [3*size];
  charges   = new double [size];
}


void Cell :: copy_positions_charges (double all_positions[], double all_charges[])
{
/** Obsolete:

    unsigned long size = particles.size();
    unsigned long  pos_index;
    printf ("**** Cell :: copy_positions_charges() :: size = %ld\n", particles.size() );
 */

  ucellindex_t   pos_index;
  for (ucellindex_t i = 0; i<particles.size(); i++) {
    pos_index = particles [i];
//    printf ("Cell :: copy_positions_charges :: pos_index = %ld\n", (long) pos_index);
    positions[3*i + 0] = all_positions [3*pos_index+0];
    positions[3*i + 1] = all_positions [3*pos_index+1];
    positions[3*i + 2] = all_positions [3*pos_index+2];
    charges [i] = all_charges[pos_index];
  }
}


//
//                                          Methods to set/get cell index
//
void Cell :: set_my_index (ucellindex_t index) {
  my_index = index;
}

ucellindex_t Cell :: get_my_index() {
  return my_index;
}


//
//                                          For Debugging/Cell Info Inspection
//
void Cell :: show_centre()
{
  printf ("my_index = %ld\tcenter = (%.4f , %.4f , %.4f)\n",
    (unsigned long) my_index, centre[0], centre[1], centre[2]);
}

void Cell :: show_positions_charges()
{
  unsigned long size = particles.size();
  printf("*** Cell %ld : %ld particles ***\n", (unsigned long) my_index, size);
  for (unsigned long i = 0; i<size; i++) {
    printf ("%5ld : (%.3f, %.3f, %.3f) ; charge = %.2f\n",
              (unsigned long)particles[i],
              positions[3*i + 0],
              positions[3*i + 1],
              positions[3*i + 2],
              charges [i]);
  }
}
