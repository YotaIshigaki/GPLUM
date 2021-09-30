/*************************************************************
 *                                                           *
 *  FMM Library for CPPMD  MD simulation program.            *
 *                                                           *
 *  Created by:  Novichkov Gleb, PhD at ERI, Inc.            *
 *  Date:        2011.01 - 2011.05                           *
 *                                                           *
 *          ~~~  Cell class.  Header file.  ~~~              *
 *                                                           *
 *************************************************************/

// *********************** CHANGES **********************
//  2011.05.17:
//    - set_num_terms() -> set_expansion_order()  (number of terms = expansion order + 1)
//    - get_num_terms() -> get_expansion_order()
//  2011.05.30:  added the following methods
//    - clear_mecs ()       -- to clear values in MECs to zero
//    - clear_lecs ()       -- to clear values in LECs to zero
//    - clear_particles ()  -- to clear particle information
//


#ifndef CELL_H
#define CELL_H

#include "defs.h"                       //  FMM definitions file
#include <complex>                      //  for MEC, LEC, etc
#include <vector>                       //  for vector class


class Cell {

  static int expansion_order;
  static int num_coeff;
  ucellindex_t   my_index;

public:
//  ===== Cell Information
  double centre [3];

//  ===== Expansion Coefficients
  std::complex <double> * MEC;  // Multipole Expansion Coefficients
  std::complex <double> * LEC;  // Local Expansion Coefficients

//  ===== Particles Information
  std::vector <ucellindex_t> particles;          //  list of particles belonging to the cell (contains their indices)
  double * positions;                       //  list of the positions of these particles (in (x, y, z) triples)
  double * charges;                         //  list of the positions of these particles (in (x, y, z) triples)

//  std::vector <double>  particles_pos;    //  list of the positions of these particles (in (x, y, z) triples)
//  std::vector <double>  charges;          //  list of the positions of these particles (in (x, y, z) triples)

  void allocate_positions_charges();
  void copy_positions_charges (double all_positions[], double all_charges[]);


//  ===== Memory Allocation methods
  void allocate();
  void deallocate();

//  ===== Index matters
  void      set_my_index (ucellindex_t);
  ucellindex_t   get_my_index();


//  ===== Debug methods
  void      show_centre();
  void      show_positions_charges();

//  ===== Service methods
  static void set_expansion_order (int n);
  static  int get_expansion_order ();
  static  int get_num_coeff ();

//  ===== Methods to clear cell data
  void  clear_mecs();
  void  clear_lecs();
  void  destroy_particles_data();

//  ===== Constructors/Destructors
  Cell();
  ~Cell();

};

#endif

