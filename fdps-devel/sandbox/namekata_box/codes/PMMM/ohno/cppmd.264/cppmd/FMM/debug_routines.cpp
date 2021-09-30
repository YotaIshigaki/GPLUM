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

#include "cell.h"
#include <stdio.h>        //  for printf()


enum TREE_ISP {MEC, LEC};

void inspect_tree(Cell ** & cell_tree, int finest_level, TREE_ISP type = MEC)
{

  printf ("Expansion Order: %d\n", Cell::get_expansion_order() );

  std::complex <double> * coeff_ptr;

  for (int l = 0; l<=finest_level; l++) {
    ucellindex_t cells_per_level = 1 << (3*l);

    printf ("Level %d : %ld cells\n", l, (unsigned long)cells_per_level);
    for (ucellindex_t i = 0; i<cells_per_level; i++) {

      if (type == MEC)
        coeff_ptr = cell_tree[l][i].MEC;
      else
        coeff_ptr = cell_tree[l][i].LEC;

      printf ("Cell %6ld : ", (unsigned long) i);

      for (int j = 0; j < Cell :: get_num_coeff(); j++)
        printf ("(%.3e, %.3e)\t", real(coeff_ptr[j]), imag(coeff_ptr[j]));

      printf ("\n");
    }  // for (i)

  }  // for (l)

}

int cmp_subtree(
  Cell ** & cell_tree1,
  Cell ** & cell_tree2,
  ucellindex_t  root_idx,
  int start_level,
  int finest_level,
  TREE_ISP type = MEC)
{

//  printf ("number of terms: %d\n", Cell::get_num_terms());

  std::complex <double> * coeff_ptr1, * coeff_ptr2;
  int retval = 0;

  for (int l = start_level; l<=finest_level; l++) {
    int depth = l - start_level;
    ucellindex_t cells_per_level = 1 << (3*depth);

//    printf ("Level %d : %ld cells\n", l, (unsigned long)cells_per_level);
    for (ucellindex_t i = 0; i<cells_per_level; i++) {

      ucellindex_t cell_idx = root_idx*(1<< (3*depth)) + i;

      if (type == MEC) {
        coeff_ptr1 = cell_tree1[l][cell_idx].MEC;
        coeff_ptr2 = cell_tree2[l][cell_idx].MEC;
      }
      else {
        coeff_ptr1 = cell_tree1[l][cell_idx].LEC;
        coeff_ptr2 = cell_tree2[l][cell_idx].LEC;
      }
//      printf ("Cell %6ld : ", (unsigned long) i);

      for (int j = 0; j < Cell :: get_num_coeff(); j++) {
//        if (coeff_ptr1[j] != coeff_ptr2[j] && l>1) {
				double delta = abs(coeff_ptr1[j] - coeff_ptr2[j]);
				double rel_delta = 0;
				if (abs(coeff_ptr1[j])>0) rel_delta = delta/abs(coeff_ptr1[j]);

        if (rel_delta>1e-8 && l>1) {
          fprintf (stdout, "*** ERROR!!! %s coefficients are different."
                  "  Level %d, cell %ld, coefficient %d :: (%.3e, %.3e) vs (%.3e, %.3e);\t |diff| = %.3e\n",
                  type == MEC ? "MEC" : "LEC", l, (unsigned long)cell_idx, j,
                  real(coeff_ptr1[j]), imag(coeff_ptr1[j]),
                  real(coeff_ptr2[j]), imag(coeff_ptr2[j]),
                  rel_delta
//                  real(coeff_ptr1[j])-real(coeff_ptr2[j]), imag(coeff_ptr1[j])-imag(coeff_ptr2[j])
                  );
          fflush (stdout);
          retval = 1;
        } else {
//          fprintf (stdout, "OK. %s coefficients are the same."
//                  "  Level %d, cell %ld, coefficient %d\n",
//                  type == MEC ? "MEC" : "LEC", l, (unsigned long)cell_idx, j);
//          fflush (stdout);
        }

      } // for (j)
//        printf ("(%.3e, %.3e)\t", real(coeff_ptr[j]), imag(coeff_ptr[j]));
//      printf ("\n");
    }  // for (i)
  }  // for (l)
  return retval;
}


void check_value (const char * routine_name, const char * var_name, std::complex <double> value)
{
  if (isnan (std::real (value)) || isnan (std::imag (value)))
    fprintf (stderr, "%s :: %s = (%.3e, %.3e)\n", routine_name, var_name, std::real (value), std::imag(value));
}

void check_value (const char * routine_name, const char * var_name, double value)
{
  if (isnan (value))
    fprintf (stderr, "%s :: %s = %.3e\n", routine_name, var_name, value);
}
