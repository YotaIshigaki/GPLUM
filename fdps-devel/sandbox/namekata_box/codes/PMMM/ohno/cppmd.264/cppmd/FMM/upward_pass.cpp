/*************************************************************
 *                                                           *
 *  FMM Library for CPPMD  MD simulation program.            *
 *                                                           *
 *  Created by:  Novichkov Gleb, PhD at ERI, In              *
 *  Date:        2011.01 - 2011.05                           *
 *                                                           *
 *              ~~~  Upward Pass Routines  ~~~               *
 *                                                           *
 *************************************************************/

#include <mpi.h>                          //  for MPI routines
#include "cell.h"
#include "coefficients.h"
#include <stdio.h>                        //  for printf(), etc
#include <stdlib.h>                       //  for exit()

using namespace std;

extern "C++" {

//  Used by:  upward_pass(), upward_pass_mpi()
//  Provided by: expansion_routines.cpp
void
compute_cell_MECs (
          Cell & cell,                      //  [in/out]
          const Coefficients * coeffs       //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
    );

//  Used by translate_and_merge_MEC(), upward_pass_mpi()
//  Provided by: expansion_routines.cpp
void
translate_MECs (
        complex <double>  new_MEC [],              //  [out]  new multipole expansion coefficients
        complex <double>  old_MEC [],              //  [in]  old multipole expansion coefficients
        const double   new_centre [3],             //  [in]  position of the new center of the expansion
        const double   old_centre [3],             //  [in]  position of the center of the original expansion
        const int &    p,                          //  [in]  expansion order
        const Coefficients * coeffs                //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
  );


//  Purpose:  Computes the index of the cell that contains the point
//  Returns:  Index of the cell
//  Used by:      sort_particles_to_cells ()
//  Provided by:  cell_routines.cpp
ucellindex_t  compute_point_index (
            const double point[],                  //  [in]  position of the point in R^d
            const int &  level,                    //  [in]  level of refinement
            const int &  dimension                 //  [in]  dimension
          );

//  Purpose:  computes center of the cell
//  Returns:  center of the cell, as array double []
//  Used by:  upward_pass(), upward_pass_mpi()
//  Provided by: cell_routines.cpp
void compute_cell_center (
          double center[],                         //  [out] center of the cell
          const ucellindex_t & cell_index,         //  [in]  index of the cell
          const int & level,                       //  [in]  level of refinement
          const int & dimension = 3                //  [in]  dimension of the system
        );
}



//  =======================================================================
//  1.  Translates MECs of the source cell to the center of the target cell
//  2.  Adds the translated MECs to the MECs of the target cell
//  Used by:  upward_pass ()
//  =======================================================================
void translate_and_merge_MEC (
        Cell & target_cell,                        //  [out] cell, MECs are translated to the center of this cell.  Then, added to MECs of this cell.
  const Cell & source_cell,                        //  [in]  cell, MECs of which are
  const Coefficients * coeffs                      //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
     )
{
  if (source_cell.MEC == NULL) {
    fprintf (stderr,  "merge_MEC() :: memory for Multipole Expansion Coefficients "\
                      "for the source cell has not been allocated.  Aborting.\n");
    exit (-1);
  }

  if (target_cell.MEC == NULL) {
    fprintf (stderr,  "merge_MEC() :: memory for Multipole Expansion Coefficients "\
                      "for the target cell has not been allocated.  Aborting.\n");
    exit (-1);
  }

//  buffer for the new MEC expansion
  int  num_coeff = Cell :: get_num_coeff();

  std::complex <double> * MEC_buffer = new std::complex <double> [num_coeff];

//  translating MECs to the new center
  translate_MECs (
    MEC_buffer,                                    //  [out] new expansion
    source_cell.MEC,                               //  [in] old expansion
    target_cell.centre,                            //  [in] new center
    source_cell.centre,                            //  [in] old center
    Cell::get_expansion_order(),                   //  [in] expansion order p
    coeffs                                         //  [in] reference to the structure with precomputed coefficients
  );

//  merging with MECs of the target cell
  for (int i = 0; i<num_coeff; i++)
    target_cell.MEC[i] += MEC_buffer[i];

  delete  [] MEC_buffer;
}

//  =======================================================================
//  Sorts Particle To Cells, at a given level (in the algorithm, the finest level is used)
//  Used by:  upward_pass ()
//  =======================================================================
void    sort_particles_to_cells (
            double position[],                     //  [in]  positions of the particles (in (x, y, z) triples)
            double charges[],                      //  [in]  charges of the particles
            unsigned long total,                   //  [in]  total number of particles
            int    level,                          //  [in]  level of refinement
            Cell   cell_array []                   //  [in/out]  array of cells
          )
{
  ucellindex_t index;

  for (unsigned long i = 0; i<total; i++) {
    index = compute_point_index (&position[3*i], level, 3);
    cell_array [index].particles.push_back(i);

#if 0 // FMM_DEBUG
    fprintf (stdout, "%5ld : (%.3f, %.3f, %.3f) \t : index = %5ld\n",
      i, position[3*i], position[3*i+1], position[3*i+2], (unsigned long)index);
#endif

#if 0
    cell_array [index].charges.push_back (charges[i]);
    cell_array [index].particles_pos.push_back (position [i+0]);
    cell_array [index].particles_pos.push_back (position [i+1]);
    cell_array [index].particles_pos.push_back (position [i+2]);
#endif
  }
}

//  ======================================================= Sequential ====
//  Upward Pass Routine
//  Used by:  main() or any other FMM method caller
//  =======================================================================
void  upward_pass (
        Cell  ** & cell_tree,                      //  [in/out]  Tree of cells
        double position[],                         //  [in]  Array containing all particles
        double charges [],                         //  [in]  Array containing all charges
        const int & finest_level,                  //  [in]  Finest level of the tree hierarchy
        const unsigned long & total,               //  [in]  Total number of particles
        const Coefficients * coeffs                //  [in]  Pointer to the coefficients structure (with pre-computed coefficients)
      )
{
  ucellindex_t index;
  ucellindex_t child_index;

//  ===== Step 1 : a) Sort particles to cells, b) create MECs for each cell =====

#if FMM_VERBOSE
  fprintf (stdout, "STEP 1\n");
  fprintf (stdout, "LEVEL %d\n", finest_level);
  fprintf (stdout, "Sorting particles to cells, computing MECs.\n");
#endif

  sort_particles_to_cells (position, charges, total, finest_level, cell_tree [finest_level]);

//  @Finest Level:
//    a) compute center of each cell
//    b) create MECs for each cell

  for (ucellindex_t i = 0; i<(1<<(3*finest_level)); i++) {

    Cell & cur_cell = cell_tree[finest_level][i];

    cur_cell.set_my_index (i);
    compute_cell_center (cur_cell.centre, i, finest_level, 3);

//  Added on 2011.05.30 by Gleb
    cur_cell.destroy_particles_data();
//  End of Gleb addition

    cur_cell.allocate_positions_charges();
    cur_cell.copy_positions_charges (position, charges);

//  Added on 2011.05.30 by Gleb
    cur_cell.clear_mecs ();
//  End of Gleb addition

    compute_cell_MECs (cur_cell, coeffs);

  }


//  ===== Step 2 :: Creating MECs of the parent by shifting and merging MECs of the children  =====
#if FMM_VERBOSE
  fprintf (stdout, "STEP 2\n");
  fprintf (stdout, "Creating MECs of the parent by shifting and merging MECs of the children\n");
#endif

  for (int l = finest_level - 1; l>=2; l--) {

#if FMM_VERBOSE
    printf ("LEVEL %d\n", l);
#endif

    for (ucellindex_t i = 0; i < (1 << (3*l)); i++) {   //  VECTORIZABLE
      Cell & parent_cell = cell_tree[l][i];
      parent_cell.set_my_index (i);
      compute_cell_center (parent_cell.centre, i, l);
//  Added on 2011.05.30 by Gleb
      parent_cell.clear_mecs();
//  End of addition
      for (int k = 0; k<8; k++) {
        child_index = (i<<3) + k;
        translate_and_merge_MEC (parent_cell/*out*/, cell_tree[l+1][child_index]/*in*/, coeffs);
      }
    }
  }   //  for (int l = finest_level - 1, ...)

}     //  upward_pass ()



//  ===================================================== MPI Parallel ====
//  Upward Pass Routine which uses MPI
//  Used by:  main()
//  =======================================================================
#define   MPI_UPWARD_PASS     (1<<12)

bool is_receiver (int myrank, int l, int rank_finest_level, int dim) {
  return myrank % (1 << ((rank_finest_level-l)*dim)) == 0 ? 1 : 0;
}


void  upward_pass_mpi (
        Cell  **  & cell_tree,                     //  [in/out]  Tree of cells
        double position[],                         //  [in]  Array containing all particles
        double charges [],                         //  [in]  Array containing all charges
        const int & myrank,                        //  [in]  rank
        const int & n_active,                      //  [in]  number of active ranks
        const int & rank_finest_level,             //  [in]  rank-finest level
        const int & finest_level,                  //  [in]  Finest level of the tree hierarchy
        const unsigned long & total,               //  [in]  Total number of particles
        const Coefficients * coeffs,               //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
        MPI_Comm fmm_world = MPI_COMM_WORLD        //  [in]  Communicator
      )
{
  ucellindex_t index;
  ucellindex_t child_index;


  MPI_Status status;
  MPI_Request request;

//  ===== Step 1 : a) Sort particles to cells, b) create MECs for each cell =====
#if FMM_VERBOSE
  fprintf (stdout, "Rank %d :: LEVEL %d\n", myrank, finest_level);
  fprintf (stdout, "Rank %d :: Sorting particles to cells, computing MECs.\n", myrank);  fflush (stdout);
#endif

// for now, can be done for all ranks.
  sort_particles_to_cells (position, charges, total, finest_level, cell_tree [finest_level]);

//  @Finest Level:
//    a) compute center of each cell
//    b) create MECs for each cell

// ************** FINEST LEVEL PROCESSING. (for now, can be done for all ranks.)
  for (ucellindex_t i = 0; i<(1<<(3*finest_level)); i++) {

    Cell & cur_cell = cell_tree[finest_level][i];

    cur_cell.set_my_index (i);
    compute_cell_center (cur_cell.centre, i, finest_level);



//  Added on 2011.05.30 by Gleb
    cur_cell.destroy_particles_data();
//  End of Gleb addition

    cur_cell.allocate_positions_charges();
    cur_cell.copy_positions_charges (position, charges);

//  Added on 2011.05.30 by Gleb
    cur_cell.clear_mecs ();
//  End of Gleb addition

    compute_cell_MECs (cur_cell, coeffs);

  }


//  ===== Step 2 :: Creating MECs of the parent by shifting and merging MECs of the children  =====

#if FMM_VERBOSE
  fprintf (stdout, "Rank %d :: Creating MECs of the parent by shifting and merging MECs of the children.\n", myrank);
  fprintf (stdout, "Rank %d :: Upward Subtree Traversal.\n", myrank);
  fflush (stdout);
#endif

//  ===== Each Rank Upward Traverses its subtree =====
  for (int l = finest_level - 1; l>=rank_finest_level; l--) {

#if FMM_VERBOSE
    fprintf (stdout, "Rank %d :: LEVEL %d\n", myrank, l); fflush (stdout);
#endif

    int depth = l - rank_finest_level;

    ucellindex_t parent_idx;
    for (ucellindex_t i = 0; i < (1 << (3*depth)); i++) {
      parent_idx = (myrank << (3*depth)) + i;

      Cell & parent_cell = cell_tree[l][parent_idx];
      parent_cell.set_my_index (parent_idx);
      compute_cell_center (parent_cell.centre, parent_idx, l);
//  Added on 2011.05.30 by Gleb
      parent_cell.clear_mecs();
//  End of addition
      for (int k = 0; k<8; k++) {
        child_index = (parent_idx<<3) + k;
        translate_and_merge_MEC (parent_cell /*out*/, cell_tree[l+1][child_index]/*in*/, coeffs);
      }
    } // end of i-loop
  } // end of l-loop


//  ====================================== Parallel-scan-like merging of MECs =====
#if FMM_VERBOSE
  fprintf (stdout, "Rank %d :: Parallel-scan-like merging of MECs\n", myrank); fflush (stdout);
#endif

//  buffers for the MEC expansion receiving/processing
  int  num_coeff = Cell :: get_num_coeff();
  std::complex <double> * received_MEC_buffer = new std::complex <double> [num_coeff];
  std::complex <double> * shifted_MEC_buffer  = new std::complex <double> [num_coeff];

  ucellindex_t cell_idx_send;
  ucellindex_t cell_idx_recv;

  bool  is_to_send = 1;

  for (int l = rank_finest_level - 1; l>=2; l--) {

    int mpi_tag = MPI_UPWARD_PASS + l;

    if (is_receiver (myrank, l, rank_finest_level, 3)) {  // receiver

      int value = 0;
      int depth_to_rf = rank_finest_level - l;

      cell_idx_send  = ((ucellindex_t)myrank) >> ((depth_to_rf - 1)*3);
      cell_idx_recv  = ((ucellindex_t)myrank) >> (depth_to_rf*3);

//    Merging MECs from child to parent by copying.
      Cell & receiving_cell = cell_tree[l][cell_idx_recv];
      Cell & sending_cell   = cell_tree[l+1][cell_idx_send];

      receiving_cell.set_my_index (cell_idx_recv);
      compute_cell_center (receiving_cell.centre, cell_idx_recv, l);

      sending_cell.set_my_index (cell_idx_send);
      compute_cell_center (sending_cell.centre, cell_idx_send, l+1);



//    Translating MECs to the center of receiving (parent) cell
      translate_MECs (
        shifted_MEC_buffer,                        //  [out] new expansion
        sending_cell.MEC,                          //  [in] old expansion
        receiving_cell.centre,                     //  [in] new center
        sending_cell.centre,                       //  [in] old center
        Cell::get_expansion_order(),               //  [in] expansion order
        coeffs                                     //  [in] reference to the structure with precomputed coefficients
      );

//  Added on 2011.05.30 by Gleb
      receiving_cell.clear_mecs();
//  End of addition

//    Merging with MECs of the receiving (parent) cell
      for (int i = 0; i<num_coeff; i++)
        receiving_cell.MEC[i] += shifted_MEC_buffer[i];


//    Receive MECs from all nodes that should send and process it.
      for (int i = 1; i<8; i++) {

        //  calculating from whom we should receive data
        int sender_rank = myrank + (i << ((depth_to_rf-1)*3));

        //  receiving MECs
        MPI_Recv(received_MEC_buffer, num_coeff, MPI::DOUBLE_COMPLEX, sender_rank, mpi_tag, fmm_world, &status);

        //  ~~~~~~~~~~~~~~~~ processing received MECs ~~~~~~~~~~~~~~~~

        //  compute index of the cell that sends MECs. Then, compute its center.
        cell_idx_send = cell_idx_recv<<3 + i;
        int level = l+1;
        Cell & sending_cell   = cell_tree[level][cell_idx_send];
        sending_cell.set_my_index (cell_idx_send);
        compute_cell_center (sending_cell.centre, cell_idx_send, level);

        //  translating MECs to the new center
        translate_MECs (
          shifted_MEC_buffer,                      //  [out]  new expansion
          received_MEC_buffer,                     //  [in]  old expansion
          receiving_cell.centre,                   //  [in]  new center (center of the cell that receives MECs)
          sending_cell.centre,                     //  [in]  old center (center of the sending cell)
          Cell::get_expansion_order(),             //  [in]  expansion order p
          coeffs                                   //  [in]  reference to the structure with precomputed coefficients
        );

        //  merging with MECs of the target cell
        for (int i = 0; i<num_coeff; i++)
          receiving_cell.MEC[i] += shifted_MEC_buffer[i];
      } // end of i-loop

    } else  { // sender
      //  check if the process should send.  It should if it did not send before (during previous iteration)
      if (is_to_send) {

        //  calculating to which rank we should send data
        int r = myrank % (1 << ((rank_finest_level-l)*3));
        int receiver_rank = myrank - r;

        //  computing read/write indices
        int depth = rank_finest_level - l;

        //  cell_idx_r  = myrank >> ((depth - 1)*3);
        cell_idx_send = myrank >> ((depth - 1)*3);

        //  preparing/selecting data to send
        Cell & sending_cell   = cell_tree[l+1][cell_idx_send];
        sending_cell.set_my_index (cell_idx_send);

        //  sending data to the recepient
        MPI_Isend (sending_cell.MEC, num_coeff, MPI::DOUBLE_COMPLEX, receiver_rank, mpi_tag, fmm_world, &request);
        //  Wait before the buffer (sending_cell.MEC) can be used (accessed) again
        MPI_Wait(&request, &status);

        //  setting flag to `false` so we do not send again
        is_to_send = false;
      }   // if (is_to_send)
    }   //  else
  }   //  end of l-loop (levels)

  delete[] shifted_MEC_buffer;
  delete[] received_MEC_buffer;

}     // upward_pass_mpi ()

