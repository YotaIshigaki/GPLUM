/*************************************************************
 *                                                           *
 *  FMM Library for CPPMD  MD simulation program.            *
 *                                                           *
 *  Created by:  Novichkov Gleb, PhD at ERI, Inc             *
 *               Ohno Yousuke, PhD at RIKEN                  *
 *  Date:        2011.01 - 2011.06                           *
 *                                                           *
 *            ~~~  Downward Pass Routines  ~~~               *
 *                                                           *
 *************************************************************/

#include <mpi.h>                          //  for MPI routines

#include "defs.h"                         //  FMM definitions file
#include "cell.h"                         //  Cell class header file
#include "coefficients.h"                 //  Precomputed coefficients class header
#include "CommMEC.h"                      //  Class to store packet (cell_idx, cell_MECs[]) data (used in downward_pass_mpi())

#include <stdio.h>                        //  for printf(), etc
#include <stdlib.h>                       //  for exit()
#include <cstring>

using namespace std;

extern "C++" {

/*************************************************************
 *      Purpose:  Computes the interaction list of a cell
 *      Returns:  Interaction list of the cell
 *                (indices of the cells on the list).
 *      Type:     A major routine
 *      Provided by: cell_routines.cpp
 **/

vector<ucellindex_t> find_interaction_list (
    const ucellindex_t &    cell_index,         //  [in]  cell's index
    const int          &    level,              //  [in]  level of refinement
    const BOUNDARY_CONDITION & 
          boundary_condition =
          DEFAULT_BOUNDARY_CONDITION,           //  [in]  boundary condition
    const int          &  dimension = 3         //  [in]  dimension
  );
vector<ucellindex_t> find_interaction_list_pbc (
						vector<double> & shift,
    const ucellindex_t &    cell_index,         //  [in]  cell's index
    const int          &    level,              //  [in]  level of refinement
    const int          &  dimension = 3         //  [in]  dimension
  );


void find_interaction_list (
    vector<ucellindex_t> & interaction_list,    //  [out] interaction list for the cell `cell_index'
    const ucellindex_t &    cell_index,         //  [in]  cell's index
    const int          &    level,              //  [in]  level of refinement
    const BOUNDARY_CONDITION &
          boundary_condition =
          DEFAULT_BOUNDARY_CONDITION,           //  [in]  boundary condition
    const int          &  dimension = 3         //  [in]  dimension
  );

/*********************************************************************************
 *      Purpose:  Computes the restricted (to a `rank') interaction list of a cell
 *      Returns:  Interaction list of the cell (indices of the cells on the list),
 *                which belong to the rank `rank'
 *      Type:     A major routine
 **/
void  find_restricted_interaction_list(
          vector<ucellindex_t> &  interaction_list,      //  [out] interaction list restricted to `rank'
          const ucellindex_t   &  cell_index,            //  [in]  cell's index
          const int            &  level,                 //  [in]  level of refinement
          const int            &  rank_finest_level,     //  [in]  rank-finest level
          const int            &  rank,                  //  [in]  rank
          const BOUNDARY_CONDITION &
                boundary_condition =
                DEFAULT_BOUNDARY_CONDITION,              //  [in]  boundary condition
          const int &             dimension = 3          //  [in]  dimension
        );
void  find_restricted_interaction_list_pbc(
          vector<ucellindex_t> &  interaction_list,      //  [out] interaction list restricted to `rank'
	  vector<double> & shift,
          const ucellindex_t   &  cell_index,            //  [in]  cell's index
          const int            &  level,                 //  [in]  level of refinement
          const int            &  rank_finest_level,     //  [in]  rank-finest level
          const int            &  rank,                  //  [in]  rank
          const int &             dimension = 3          //  [in]  dimension
        );


void  compute_cell_center (
          double center[],                               //  [out] center of the cell
          const ucellindex_t & cell_index,               //  [in]  index of the cell
          const int & level,                             //  [in]  level of refinement
          const int & dimension = 3                      //  [in]  dimension of the system
        );

//  Translates LECs -> LECs
void  translate_LECs (
          complex <double>  new_LEC [],                  //  [out]  new multipole expansion coefficients
          complex <double>  old_LEC [],                  //  [in]  old multipole expansion coefficients
          const double   new_centre [3],                 //  [in]  position of the new center of the expansion
          const double   old_centre [3],                 //  [in]  position of the center of the original expansion
          const int &    p,                              //  [in]  expansion order
          const Coefficients * coeffs                    //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
        );


//  =====================  Converts MECs -> LECs  ===
//  Output: Multipole Expansion Coefficients (array)
//  =================================================
void  convert_MEC2LEC (
          complex <double>  new_LEC [],                  //  [out] LEC coefficients
          complex <double>  old_MEC [],                  //  [in]  MEC coefficients
          const double   LEC_centre [3],                 //  [in]  position of the center local expansion
          const double   MEC_centre [3],                 //  [in]  position of the center multipole expansion
          const int &    p,                              //  [in]  expansion order
          const Coefficients * coeffs,                   //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
          complex <double> * sph_harmonic_storage        //  [in]  pointer to the pre-allocated buffer to store Harmonic coefficients
        );
void  convert_MEC2LEC_pbc (
          complex <double>  new_LEC [],                  //  [out] LEC coefficients
          complex <double>  old_MEC [],                  //  [in]  MEC coefficients
          const double   LEC_centre [3],                 //  [in]  position of the center local expansion
          const double   MEC_centre [3],                 //  [in]  position of the center multipole expansion
          const int &    p,                              //  [in]  expansion order
          const Coefficients * coeffs,                   //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
          complex <double> * sph_harmonic_storage        //  [in]  pointer to the pre-allocated buffer to store Harmonic coefficients
        );

//  ================ Computes potential from LECs ===
//  Output: Potential Value (complex)
//  =================================================
complex<double>
  evaluate_LEC (
          double position [3],
          double center   [3],
          complex<double> LEC[],
          int num_terms
        );

//  ====== Computes potential and force from LECs ===
//  Output: Potential Value (complex)
//  =================================================
void LEC_to_potential_force (
          double & potential,                            //  [in/out] potential (updated by "+=" )
          double field[3],                               //  [out] field
          double position [3],                           //  [in]  Particle's position
          double center   [3],                           //  [in]  Center of the cell
          const complex<double> LEC[],                   //  [in]  LECs
          const int & expansion_order,                   //  [in]  order of expansion
          const Coefficients * coeffs                    //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
        );

//  ====== Computes potential and force from LECs ===
//  Output: List of nearest neighbors (cell indices)
//  =================================================
vector<ucellindex_t>
  find_nearest_neighbors(
          const ucellindex_t  &   cell_index,            //  [in]  index of the cell
          const int           &   level,                 //  [in]  level of refinement
          const int           &   dimension,             //  [in]  dimension
          const BOUNDARY_CONDITION
                boundary_condition =
                DEFAULT_BOUNDARY_CONDITION               //  [in]  boundary condition
        );

}
//  ======================================== End of External Routines =====



void  compute_potential_force (
          double my_position[3],                         // [in]
          double other_position[3],                      // [in]
          double my_charge,                              // [in]
          double other_charge,                           // [in]
          double & my_potential,                         // [in/out] potential value is updated (+=)
          double my_force[3]                             // [in/out] force value is updated (+=)
        );
void  compute_potential_force_pbc (
          double my_position[3],                         // [in]
          double other_position[3],                      // [in]
          double my_charge,                              // [in]
          double other_charge,                           // [in]
          double & my_potential,                         // [in/out] potential value is updated (+=)
          double my_force[3]                             // [in/out] force value is updated (+=)
        );
void  compute_potential_force_shifted (
          double my_position[3],                         // [in]
          double other_position[3],                      // [in]
          double my_charge,                              // [in]
          double other_charge,                           // [in]
	  double shift[3],                               // [in]
          double & my_potential,                         // [in/out] potential value is updated (+=)
          double my_force[3]                             // [in/out] force value is updated (+=)
        );


//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Downward Pass Routines ~~~~~~~
//

//  ======================================================= Sequential ====
//  Downward Pass Routine
//  Used by:  main()
//  =======================================================================
void  downward_pass (
        double & total_potential,                        //  [out] computed total potential
#ifdef  USE_POTENTIAL_ARRAY
        double  potential [],                            //  [out] computed potential for each particle
#endif
        double  force [],                                //  [out] computed force (in triples (Fx, Fy, Fz).
        Cell  ** &  cell_tree,                           //  [in/out]  Tree of cells
        int     finest_level,                            //  [in]  Finest level of the tree hierarchy
        double  position [],                             //  [in]  Array containing all particles
        double  charges [],                              //  [in]  Array containing all charges
        const Coefficients * coeffs,                     //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
        complex <double>   * sph_harmonic_storage,       //  [in]  pointer to the pre-allocated buffer to store Harmonic coefficients
        const BOUNDARY_CONDITION &
              boundary_condition =
              DEFAULT_BOUNDARY_CONDITION                 //  [in]  Boundary Condition
      )
{
  ucellindex_t child_idx, parent_idx;
  vector<ucellindex_t> interaction_list;

//  NOTE: initialization of LECs at level 1 to zero is not necessary because
//        LECs are already initialized to zero by the cell constructor

// Allocating buffer to store shifted LECs of the parent
  int  num_coeff = Cell :: get_num_coeff();
  complex <double> * LEC_buffer = new complex <double> [num_coeff];


//  ======== Downward Pass ========
  for (int level = 2; level<=finest_level; level++) {

#if FMM_VERBOSE
    fprintf (stdout, "LEVEL %d\n", level);  fflush (stdout);
#endif

//  loop over cells on level `level'.
    for (ucellindex_t i = 0; i < (1 << (3*level)) /*8^level*/; i++) {
      parent_idx = i>>3;

      Cell & cur_cell = cell_tree[level][i];
      Cell & parent_cell = cell_tree[level-1][parent_idx];

//  ===== Step 3a):  shifting parent's LECs to the cell i =====



//  shifting LECs to the center of the current cell and store them in the LEC of the current cell
      translate_LECs (
        cur_cell.LEC,                                    //  [out] new expansion.    direct writing, NO "+="
        parent_cell.LEC,                                 //  [in] old expansion
        cur_cell.centre,                                 //  [in] new center
        parent_cell.centre,                              //  [in] old center
        Cell::get_expansion_order(),                     //  [in] expansion order
        coeffs                                           //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
      );


//  ===== Step 3b):  Computing far field contribution from the cells in the interaction list
      interaction_list = find_interaction_list (i, level, boundary_condition);
//        find_interaction_list (interaction_list /*out*/, i, level, boundary_condition);


      for (unsigned long k = 0; k < interaction_list.size(); k++) {
        ucellindex_t  idx     = interaction_list[k];
        Cell &   il_cell = cell_tree[level][idx];


        // converting MECs (interacting cell) -> LECs (current cell i) (stored to LEC_buffer)
        // TODO:  perhaps, we can do "+=" for the output and use cell LECs instead of buffer.  Maybe it will be faster.
        convert_MEC2LEC (
          LEC_buffer,                                    //  [out] LEC coefficients.    direct writing, NO "+="
          il_cell.MEC,                                   //  [in]  MEC coefficients
          cur_cell.centre,                               //  [in]  position of the new center of the expansion
          il_cell.centre,                                //  [in]  position of the center of the original expansion
          Cell::get_expansion_order(),                   //  [in] expansion order
          coeffs,                                        //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
          sph_harmonic_storage
        );

        // adding the converted LECs to the LECs of the current cell
        for (int ic = 0; ic<Cell::get_num_coeff(); ic++)
          cur_cell.LEC[ic] += LEC_buffer[ic];

      }   //  end of interaction list k-loop

      interaction_list.empty();

    }   //  end of the cells i-loop

  }   //  end of the levels level-loop
  delete [] LEC_buffer;


//  ===== Steps 4 & 5: potential at finest level  =====
  ucellindex_t particle_id;
  vector<ucellindex_t> nearest_cells;

  for (ucellindex_t j = 0; j < (1 << (3*finest_level)) /*8^level*/; j++) {
    Cell &  cur_cell = cell_tree [finest_level][j];
    nearest_cells = find_nearest_neighbors (j, finest_level, 3, boundary_condition/*ZERO*/);

    //  loop over particles in the cell
    for (unsigned long i = 0; i<cur_cell.particles.size(); i++) {
      particle_id = cur_cell.particles[i];

//  ===== Step 4:  evaluating Far Field potential (LECs) at finest level  =====

      double field[3] = {0, 0, 0};

      double poten_tmp = 0;

      LEC_to_potential_force (
              poten_tmp,                                 //  [in/out] potential (updated by "+=" )
              field,                                     //  [out] field
              &cur_cell.positions[3*i],                  //  [in]  Particle's position
              cur_cell.centre,                           //  [in]  Center of the cell
              cur_cell.LEC,                              //  [in]  LECs
              Cell::get_expansion_order(),               //  [in] expansion order
              coeffs                                     //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
            );

      total_potential+=poten_tmp;

#ifdef  USE_POTENTIAL_ARRAY
      potential [particle_id]+=poten_tmp;
#endif

      double & charge = cur_cell.charges[i];
      force [particle_id*3 + 0] = charge*field [0];
      force [particle_id*3 + 1] = charge*field [1];
      force [particle_id*3 + 2] = charge*field [2];


//  ===== Step 5:  Direct interaction with particles in the nearest cells, at finest level  =====

      //  loop over nearest cells (to get their particles)
      ucellindex_t idx;

      double near_pot = 0;
      double my_force[3] = {0., 0., 0.};

      for (int kc = 0; kc < nearest_cells.size(); kc++) {
        idx = nearest_cells[kc];  // index of the nearest cell
        Cell &  near_cell = cell_tree[finest_level][idx];


        //  loop over particles in the current nearest cell
        ucellindex_t near_particle_id;
        for (unsigned long kp = 0; kp < near_cell.particles.size(); kp++) {
          near_particle_id = near_cell.particles [kp];
          if (particle_id != near_particle_id) {
            compute_potential_force (
               &cur_cell.positions [3*i],                // [in]  position of the current particle
               &near_cell.positions [3*kp],              // [in]  position of the other particle
               cur_cell.charges[i],                      // [in]  charge of the current particle
               near_cell.charges [kp],                   // [in]  charge of the other particle
               near_pot,                                 // [in/out]  current particle's potential value is updated (+=)
               my_force                                  // [in/out]  current particle's force is updated (+=)
            );
          }  // if (particle_id != near_particle_id)

        }
      }   //  end of kc-loop over nearest cells

#ifdef  USE_POTENTIAL_ARRAY
      potential [particle_id] += near_pot;
#endif
      total_potential +=near_pot;


      force [3*particle_id + 0] += my_force [0];
      force [3*particle_id + 1] += my_force [1];
      force [3*particle_id + 2] += my_force [2];

    }   //  end of i-loop over particles in the cell
  }   //  end of j-loop over the cells

}  //  downward_pass()


//  ===================================================== MPI Parallel ====
//  Downward Pass Routine which uses MPI
//  Used by:  main()
//  =======================================================================
#define   MPI_DOWNWARD_PASS             (1<<13)
#define   MPI_DOWNWARD_PASS_TOP2RF_IL   (1<<14)
#define   MPI_DOWNWARD_PASS_RF2FINEST   (1<<15)


void  downward_pass_mpi (
        double & total_potential,                        //  [out] computed total potential
#ifdef  USE_POTENTIAL_ARRAY
        double  potential [],                            //  [out] computed potential for each particle
#endif
        double  force [],                                //  [out] computed force (in triples (Fx, Fy, Fz).
        Cell  ** &  cell_tree,                           //  [in/out]  Tree of cells
        int     myrank,                                  //  [in]  Rank
        int     rank_finest_level,                       //  [in]  Rank-finest level
        int     finest_level,                            //  [in]  Finest level of the tree hierarchy
        double  position [],                             //  [in]  Array containing all particles
        double  charges [],                              //  [in]  Array containing all charges
        const Coefficients * coeffs,                     //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
        complex <double> * sph_harmonic_storage,         //  [in]  pointer to the pre-allocated buffer to store Harmonic coefficients
        const BOUNDARY_CONDITION &
              boundary_condition =
              DEFAULT_BOUNDARY_CONDITION,                //  Boundary Condition
        MPI_Comm fmm_world = MPI_COMM_WORLD              //  communicator
      )
{

//  MPI vars
  MPI_Status status;
  MPI_Request request;


  ucellindex_t child_idx, parent_idx;
  vector<ucellindex_t> interaction_list;


//  allocating buffer to store shifted LECs of the parent
  int  num_coeff = Cell :: get_num_coeff();

  std::complex <double> * LEC_buffer = new std::complex <double> [num_coeff];
  std::complex <double> * received_LEC_buffer = new std::complex <double> [num_coeff];
// NOTE: received_MEC_buffer might not be needed in the future
  std::complex <double> * received_MEC_buffer = new std::complex <double> [num_coeff];



//  ======================================================================================= Downward Pass ========



//
//  =========================================================== Downward Pass:  Level 2 -> Rank-finest level ===
//
  bool  is_sender   = false;
  int   start_level = 1;
  int   stop_level  = rank_finest_level;

  vector<double> pbc_shift;

  if (myrank == 0)
    is_sender = true;

  //  printf("Rank:%d  rank_finest_level %d\n",myrank,rank_finest_level); fflush (stdout);
  for (int level = 1; level <= rank_finest_level; level++) {
    // depth from the current cell to stop_level
    int depth   = rank_finest_level - level;
    int mpi_tag = MPI_DOWNWARD_PASS + level;

    if ( myrank % (1 << (depth*3)) == 0) {     //  active, sender or receiver

      //      printf("Rank:%d active\n",myrank);fflush (stdout);

//    computing indexes of the parent cell/child cell (child cell for the sending rank)
      parent_idx = myrank >> ((depth+1)*3);
      child_idx  = myrank >> (depth*3);

//    References to Parent Cell, Child Cell
      Cell & parent_cell = cell_tree [level-1][parent_idx];
      Cell & child_cell  = cell_tree [level][child_idx];

//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sender/Receiver ~~~~~~~~~~~~
      if (is_sender) {    //  Sender
	//	printf("Rank:%d is_sender level %d\n",myrank,level);
//      in case of 'sender' child is processed by the same rank, so storing directly.
        translate_LECs (
          child_cell.LEC,                                //  [out] new expansion.    direct writing, NO "+="
          parent_cell.LEC,                               //  [in] old expansion
          child_cell.centre,                             //  [in] new center
          parent_cell.centre,                            //  [in] old center
          Cell::get_expansion_order(),                   //  [in] expansion order
          coeffs                                         //  [in] pointer to the coefficients structure (with pre-computed coefficients)
        );


//      send LECs of the parent to the children (computing their ranks first)
        for (int i = 1; i<8; i++) {
          int receiver_rank = myrank + (i << (depth*3));   /* recepient rank */
          MPI_Isend (parent_cell.LEC, num_coeff, MPI::DOUBLE_COMPLEX, receiver_rank, mpi_tag, fmm_world , &request);
          //  wait before the buffer (parent_cell.LEC) can be used (accessed) again
          MPI_Wait(&request, &status);
        }

      } else {           // Receiver
	//	printf("Rank:%d ! is_sender level %d\n",myrank,level);
//  Receiving LECs of the parent.  Then, shifting them to the child.

        //  calculating from which rank we should receive data
        int r = myrank % (1 << ((depth + 1)*3));
        int sender_rank = myrank - r;

        //  receiving data from the sender
        MPI_Recv(received_LEC_buffer, num_coeff, MPI::DOUBLE_COMPLEX, sender_rank, mpi_tag, fmm_world , &status);

        //  received data once, we don't receive it again.  only send
        is_sender = true;

        //  processing received data
        translate_LECs (
          child_cell.LEC,                                //  [out] new expansion
          received_LEC_buffer,                           //  [in] old expansion
          child_cell.centre,                             //  [in] new center
          parent_cell.centre,                            //  [in] old center
          Cell::get_expansion_order(),                   //  [in] expansion order
          coeffs                                         //  [in] pointer to the coefficients structure (with pre-computed coefficients)
        );
      }
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF Sender/Receiver ~~~~~~~~~~~~



//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Contribution from Interaction List Cells ~~~~~
      double il_cell_center[3];
      if(boundary_condition==PBC){
	pbc_shift.clear();
	if(level==1){
	  interaction_list.clear();
	}else{
	  interaction_list = find_interaction_list_pbc (pbc_shift, child_idx, level);
	//interaction_list = find_interaction_list (child_idx, level, boundary_condition);
	  printf("Rank %d level %d size of interaction_list %d\n",myrank,level,interaction_list.size());
	if(myrank==0)
	  {
	    for(unsigned long k = 0; k < interaction_list.size(); k++) {
	      ucellindex_t  il_cell_idx     = interaction_list[k];
	      printf("0x%08x 0x%08x %f %f %f\n",child_idx,il_cell_idx,pbc_shift[k*3],pbc_shift[k*3+1],pbc_shift[k*3+2]);
	      //	      printf("0x%08x 0x%08x\n",child_idx,il_cell_idx);
	    }
	  }
	}
      }else{
	interaction_list = find_interaction_list (child_idx, level, boundary_condition);
//      find_interaction_list (interaction_list, child_idx, level, boundary_condition);
      }

      for (unsigned long k = 0; k < interaction_list.size(); k++) {
        ucellindex_t  il_cell_idx     = interaction_list[k];

        int il_cell_rank = RANK(il_cell_idx, depth);
        int tag = MPI_DOWNWARD_PASS_TOP2RF_IL;

	{
	  
	}

        MPI_Sendrecv (child_cell.MEC,      num_coeff,  MPI::DOUBLE_COMPLEX, il_cell_rank, tag,
                      received_MEC_buffer, num_coeff,  MPI::DOUBLE_COMPLEX, il_cell_rank, tag,     fmm_world, &status);

        // setting index and center of the interacting cell
        compute_cell_center (il_cell_center, il_cell_idx, level);

        // converting MECs (interacting cell) -> LECs (current cell i) (stored to LEC_buffer)
        // TODO:  perhaps, we can do "+=" for the output and use cell LECs instead of buffer.  Maybe it will be faster.
	if(boundary_condition==PBC){
	  // add pbc_shift
	  //	  if(myrank==0)
	    {
	      //	      printf("0x%08x 0x%08x %f %f %f\n",child_idx,il_cell_idx,pbc_shift[k*3],pbc_shift[k*3+1],pbc_shift[k*3+2]);
	      //printf("0x%08x 0x%08x\n",child_idx,il_cell_idx);
	  }
	}
        convert_MEC2LEC (
          LEC_buffer,                                    //  [out] LEC coefficients.    direct writing, NO "+="
          received_MEC_buffer,                           //  [in]  MEC coefficients
          child_cell.centre,                             //  [in]  position of the new center of the expansion
          il_cell_center,                                //  [in]  position of the center of the original expansion
          Cell::get_expansion_order(),                   //  [in] expansion order
          coeffs,                                        //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
          sph_harmonic_storage                           //  [in]  pointer to the pre-allocated buffer to store Harmonic coefficients
        );

        // adding the converted LECs to the LECs of the current cell
        for (int ic = 0; ic<Cell::get_num_coeff(); ic++)
          child_cell.LEC[ic] += LEC_buffer[ic];

      }   //  end of interaction list loop
      interaction_list.empty();
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF Contribution from Interaction List Cells' ~~~~~

    }  //  end of 'active, sender or receiver'
    else{
      //      printf("Rank:%d inactive  myrank % (1 << (%d*3)) %d\n",myrank, depth, myrank % (1 << (depth*3)));fflush (stdout);
    }

  }   //  End of level-loop (level=2,.., rank_finest_level)


//
//  =========================================================== Downward Pass:  Rank-Finest level -> Finest level ===
//

  CommMEC :: set_num_coeff (num_coeff);
  int n_active = 1 << (3*rank_finest_level);
  bool * is_rank_receiver = new bool [n_active];

//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Level Loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (int level = rank_finest_level+1; level <= finest_level; level++) {
    int depth = level - rank_finest_level;            // depth from rank_finest_level to current level

#if FMM_VERBOSE
    fprintf (stdout, "RANK %d :: LEVEL %3d  depth = %d.  RF -> FINEST traversal.\n", myrank, level, depth);  fflush (stdout);
#endif

//  Performing 2 tasks:
//  1. Shifting Parent LECs to the children
//  2. Processing cells from the interaction list

    vector <CommMEC>  send_comm_mec (n_active);
    vector <CommMEC>  recv_comm_mec (n_active);

    //  Loop over cells on this level (owned by the rank)
    for (ucellindex_t ic = 0; ic < (1 << (3*depth)); ic++) {

      ucellindex_t cell_idx  = (myrank<<(depth*3)) + ic;
      Cell & cur_cell = cell_tree [level][cell_idx];

      parent_idx = cell_idx >> 3;
      Cell & parent_cell = cell_tree [level-1][parent_idx];
      compute_cell_center (parent_cell.centre, parent_idx, level-1);

//    shift LECs from parent cell to the child
      translate_LECs (
        cur_cell.LEC,                                    //  [out] new expansion.    direct writing, NO "+="
        parent_cell.LEC,                                 //  [in] old expansion
        cur_cell.centre,                                 //  [in] new center
        parent_cell.centre,                              //  [in] old center
        Cell::get_expansion_order(),                     //  [in] expansion order
        coeffs
      );



//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Contribution from Interaction List Cells ~~~~~

      memset (is_rank_receiver, false, n_active*sizeof(bool));
      interaction_list = find_interaction_list (cell_idx, level, boundary_condition);
//      find_interaction_list (interaction_list, cell_idx, level, boundary_condition);
      // Interaction List contribution
      for (unsigned long k = 0; k < interaction_list.size(); k++) {

        ucellindex_t  il_cell_idx  = interaction_list[k];
        int il_cell_rank = RANK(il_cell_idx, depth);


        if (il_cell_rank == myrank) {

          //  ~~~  Converting MECs (interacting cell) -> LECs (current cell ic) (stored to LEC_buffer)  ~~~
          double il_cell_center[3];
          compute_cell_center (il_cell_center, il_cell_idx, level);


          // TODO:  perhaps, we can do "+=" for the output and use cell LECs instead of buffer.  Maybe it will be faster.
          convert_MEC2LEC (
            LEC_buffer,                                  //  [out] LEC coefficients.    direct writing, NO "+="
            cell_tree [level][il_cell_idx].MEC,          //  [in]  MEC coefficients
            cur_cell.centre,                             //  [in]  position of the new center of the expansion
            il_cell_center,                              //  [in]  position of the center of the original expansion
            Cell::get_expansion_order(),                 //  [in]  expansion order
            coeffs,                                      //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
            sph_harmonic_storage                         //  [in]  pointer to the pre-allocated buffer to store Harmonic coefficients
          );


          //  ~~~  Adding the converted LECs to the LECs of the current cell ~~~
          for (int ic = 0; ic<Cell::get_num_coeff(); ic++)
            cur_cell.LEC[ic] += LEC_buffer[ic];

        }  else  {

          is_rank_receiver [il_cell_rank] = true;

        }  // if {} else {}

      } // end of k-loop (interaction list loop)

      interaction_list.empty();

      for (int ir = 0; ir<n_active; ir++)  if (is_rank_receiver [ir]) {
        send_comm_mec[ir].add (cell_idx, cur_cell.MEC);
      }

    } // end of ic-loop (loop over cells on this level)



//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Exchanging buffered (idx{k}, MECs{k}) data ~~~

    //  Exchanging the sizes of packed CommMECs buffers, allocating buffers for receiving data
    for (int irank = 0; irank<n_active; irank++) {
      if (irank == myrank)  continue;

      int tag =  MPI_DOWNWARD_PASS_RF2FINEST + level;
      MPI_Sendrecv (&send_comm_mec[irank].N,   1,   MPI_INT,   irank,   tag,
                    &recv_comm_mec[irank].N,   1,   MPI_INT,   irank,   tag,  fmm_world, &status);
      recv_comm_mec[irank].allocate();


    } // irank-loop

    //  Exchanging buffered ComMEC data between ranks
    for (int irank = 0; irank<n_active; irank++) {
      if (irank == myrank)  continue;
      int tag =  MPI_DOWNWARD_PASS_RF2FINEST + level;
      MPI_Sendrecv (send_comm_mec[irank].data,  send_comm_mec[irank].data_size,  MPI_BYTE,  irank,  tag,
                    recv_comm_mec[irank].data,  recv_comm_mec[irank].data_size,  MPI_BYTE,  irank,  tag,  fmm_world,  &status);
    } // irank-loop





//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Processing received buffered (idx{k}, MECs{k}) data ~~~
    for (int irank = 0; irank<n_active; irank++) {

      ucellindex_t  sender_idx;
      ucellindex_t  receiver_idx;
      vector<ucellindex_t>  sender_rint_list;
      complex <double> * received_MECs;


      for (int isc=0; isc < recv_comm_mec[irank].N; isc++) {
        sender_idx    = recv_comm_mec[irank].get_idx(isc);
        received_MECs = recv_comm_mec[irank].get_MECs(isc);

	if(boundary_condition==PBC){
	  pbc_shift.clear();
        find_restricted_interaction_list_pbc (sender_rint_list /*out*/,
					      pbc_shift,
                                          sender_idx,
                                          level,
                                          rank_finest_level,
                                          myrank);
	}else{
        // find iteraction list of `sender_idx' restricted to the rank `myrank'
        find_restricted_interaction_list (sender_rint_list /*out*/,
                                          sender_idx,
                                          level,
                                          rank_finest_level,
                                          myrank,
                                          boundary_condition);
	}
        for (int ilc = 0; ilc<sender_rint_list.size(); ilc++) {

          //  Shifting MECs of `sender_idx' (that came from rank `irank') to LECs of `receiver_idx' (that belongs to `myrank').
          receiver_idx = sender_rint_list [ilc];
          Cell & receiver_cell = cell_tree[level][receiver_idx];

          double   sender_center[3];
          compute_cell_center (sender_center, sender_idx, level);

	  //	  printf("0x%08x 0x%08x %f %f %f\n",sender_idx,receiver_idx,pbc_shift[ilc*3],pbc_shift[ilc*3+1],pbc_shift[ilc*3+2]);
          // TODO:  perhaps, we can do "+=" for the output and use cell LECs instead of buffer.  Maybe it will be faster.
	  if(boundary_condition==PBC){
	    sender_center[0] += pbc_shift[ilc*3];
	    sender_center[1] += pbc_shift[ilc*3+1];
	    sender_center[2] += pbc_shift[ilc*3+2];
	  }
          convert_MEC2LEC (
            LEC_buffer,                                  //  [out] LEC coefficients.    direct writing, NO "+="
            received_MECs,                               //  [in]  received MEC coefficients
            receiver_cell.centre,                        //  [in]  position of the new center of the expansion
            sender_center,                               //  [in]  position of the center of the original expansion
            Cell::get_expansion_order(),                 //  [in] expansion order
            coeffs,                                      //  [in]  pointer to the coefficients structure (with pre-computed coefficients)
            sph_harmonic_storage                         //  [in]  pointer to the pre-allocated buffer to store Harmonic coefficients
          );

          //  Adding the converted LECs to the LECs of the receiver cell
          for (int ic = 0; ic<Cell::get_num_coeff(); ic++)
            receiver_cell.LEC[ic] += LEC_buffer[ic];

        }  // end of ilc-loop (loop over the cells in the interaction list of (sender_idx) restricted to `myrank'

        sender_rint_list.clear();
      }  // end of isc-loop
    }  // end of irank-loop
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END of Processing received buffered (idx{k}, MECs{k}) data ~~~
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF Contribution from Interaction List Cells' ~~~~~

  }  // end of level-loop (level)
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END of Level Loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//  Releasing buffers as we do not need them anymore
  delete [] is_rank_receiver;
  delete [] LEC_buffer;
  delete [] received_LEC_buffer;
  delete [] received_MEC_buffer;

//
//  ==================================================== Steps 4 & 5: Potential/Force compution at finest level  ===
//
  ucellindex_t particle_id, cell_idx;
  vector<ucellindex_t> nearest_cells;
  int depth = finest_level - rank_finest_level;

  //  loop over cells in this rank
  for (ucellindex_t ic = 0; ic < (1 << (3*depth)); ic++) {
    cell_idx  = (myrank<<(depth*3)) + ic;
    nearest_cells = find_nearest_neighbors (cell_idx, finest_level, 3, boundary_condition);
    Cell &  cell = cell_tree [finest_level][cell_idx];

    //  loop over particles in the cell
    for (unsigned long jp = 0; jp<cell.particles.size(); jp++) {
      particle_id = cell.particles[jp];

//  ===== Step 4:  evaluating Far Field potential (LECs) at finest level  =====

      double field [3] = {0, 0, 0};
      double poten_tmp = 0;

      LEC_to_potential_force (
              poten_tmp,                                 //  [in/out] potential (updated by "+=" )
              field,                                     //  [in/out] field
              &cell.positions[3*jp],                     //  [in]  Particle's position
              cell.centre,                               //  [in]  Center of the cell
              cell.LEC,                                  //  [in]  LECs
              Cell::get_expansion_order(),               //  [in] expansion order
              coeffs
            );

      double & charge = cell.charges[jp];

      total_potential += charge*poten_tmp;

#ifdef  USE_POTENTIAL_ARRAY
      potential [particle_id] += charge*poten_tmp;
#endif

      force [particle_id*3 + 0] = charge*field [0];
      force [particle_id*3 + 1] = charge*field [1];
      force [particle_id*3 + 2] = charge*field [2];


//  ===== Step 5:  Direct interaction with particles in the nearest cells, at finest level  =====

      ucellindex_t idx;
      double near_pot = 0, my_force[3] = {0.,0.,0.};

      //  loop over nearest cells (to get their particles)
      for (int kc = 0; kc < nearest_cells.size(); kc++) {
        idx = nearest_cells[kc];
        Cell &  near_cell = cell_tree [finest_level][idx];

        //  loop over particles in the current nearest cell
        ucellindex_t near_particle_id;
        for (unsigned long kp = 0; kp < near_cell.particles.size(); kp++) {
          near_particle_id = near_cell.particles [kp];
          if (particle_id != near_particle_id)
	    if(boundary_condition==PBC){
            compute_potential_force_pbc (
               &cell.positions [3*jp],                   // [in]  position of the current particle
               &near_cell.positions [3*kp],              // [in]  position of the other particle
               cell.charges[jp],                         // [in]  charge of the current particle
               near_cell.charges [kp],                   // [in]  charge of the other particle
               near_pot,                                 // [in/out]  current particle's potential value is updated (+=)
               my_force                                  // [in/out]  current particle's force is updated (+=)
            );
	    }else{
            compute_potential_force (
               &cell.positions [3*jp],                   // [in]  position of the current particle
               &near_cell.positions [3*kp],              // [in]  position of the other particle
               cell.charges[jp],                         // [in]  charge of the current particle
               near_cell.charges [kp],                   // [in]  charge of the other particle
               near_pot,                                 // [in/out]  current particle's potential value is updated (+=)
               my_force                                  // [in/out]  current particle's force is updated (+=)
            );
	    }
        }  //  end of kp-loop over particles in the nearest cells
      }   //  end of kc-loop over nearest cells


      total_potential += charge*near_pot;

#ifdef  USE_POTENTIAL_ARRAY
      potential [particle_id] += charge*near_pot;
#endif

      force [3*particle_id + 0] += my_force [0];
      force [3*particle_id + 1] += my_force [1];
      force [3*particle_id + 2] += my_force [2];

    }   //  end of j-loop over particles in the cell
  }   //  end of i-loop over the cells
//  ==================================================== END of Steps 4 & 5: Potential/Force compution at finest level  ===

}   //  downward_pass_mpi()


//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Computational Routines-Helpers  ~~~~~~~~~~~~~~~
//

//  computes potential and force between two charges
void compute_potential_force (
          double my_position[3],        // [in]
          double other_position[3],     // [in]
          double my_charge,             // [in]
          double other_charge,          // [in]
          double & my_potential,        // [in/out] potential value is updated (+=)
          double my_force[3]            // [in/out] force value is updated (+=)
        )
{
  double rho_sq, rho;
  double dr[3];
  rho_sq = 0;

  for (int i = 0; i<3; i++) {
    dr [i] = my_position[i] - other_position[i];
    rho_sq += dr[i] * dr[i];
  }

  rho = sqrt (rho_sq);

  if (rho> 0) {
    my_potential += other_charge/rho;
    double rho_rho_sq = rho*rho_sq;
    my_force [0] += my_charge*other_charge*dr[0]/rho_rho_sq;
    my_force [1] += my_charge*other_charge*dr[1]/rho_rho_sq;
    my_force [2] += my_charge*other_charge*dr[2]/rho_rho_sq;
  } else {
    fprintf (stderr, "FMM (downward pass) ::  compute_potential_force() :: rho = %e\n", rho); fflush (stderr);
  }

}

void compute_potential_force_pbc (
          double my_position[3],        // [in]
          double other_position[3],     // [in]
          double my_charge,             // [in]
          double other_charge,          // [in]
          double & my_potential,        // [in/out] potential value is updated (+=)
          double my_force[3]            // [in/out] force value is updated (+=)
        )
{
  double rho_sq, rho;
  double dr[3];
  rho_sq = 0;

  for (int i = 0; i<3; i++) {
    dr [i] = my_position[i] - other_position[i];
    if(dr[i]>0.5)dr[i]-=1.0;
    if(dr[i]<-0.5)dr[i]+=1.0;
    rho_sq += dr[i] * dr[i];
  }

  rho = sqrt (rho_sq);

  if (rho> 0) {
    my_potential += other_charge/rho;
    double rho_rho_sq = rho*rho_sq;
    my_force [0] += my_charge*other_charge*dr[0]/rho_rho_sq;
    my_force [1] += my_charge*other_charge*dr[1]/rho_rho_sq;
    my_force [2] += my_charge*other_charge*dr[2]/rho_rho_sq;
  } else {
    fprintf (stderr, "FMM (downward pass) ::  compute_potential_force() :: rho = %e\n", rho); fflush (stderr);
  }

}

void compute_potential_force_shifted (
          double my_position[3],        // [in]
          double other_position[3],     // [in]
          double my_charge,             // [in]
          double other_charge,          // [in]
	  double shift[3],              // [in]
          double & my_potential,        // [in/out] potential value is updated (+=)
          double my_force[3]            // [in/out] force value is updated (+=)
        )
{
  double rho_sq, rho;
  double dr[3];
  rho_sq = 0;

  for (int i = 0; i<3; i++) {
    dr [i] = my_position[i] - other_position[i] + shift[i];
    rho_sq += dr[i] * dr[i];
  }

  rho = sqrt (rho_sq);

  if (rho> 0) {
    my_potential += other_charge/rho;
    double rho_rho_sq = rho*rho_sq;
    my_force [0] += my_charge*other_charge*dr[0]/rho_rho_sq;
    my_force [1] += my_charge*other_charge*dr[1]/rho_rho_sq;
    my_force [2] += my_charge*other_charge*dr[2]/rho_rho_sq;
  } else {
    fprintf (stderr, "FMM (downward pass) ::  compute_potential_force() :: rho = %e\n", rho); fflush (stderr);
  }

}

//  computes potential between two charges
double  compute_potential (
          double my_position[3],
          double other_position[3],
          double other_charge
        )
{
  double rho, delta, potential;
  rho = 0;
  for (int i = 0; i<3; i++) {
    delta = my_position[i] - other_position[i];
    rho += delta * delta;
  }

  potential = other_charge/sqrt(rho);

  return potential;
}


//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Routines to broadcast/collect data by MPI  ~~~~~
//

//  ===================================================== MPI Parallel ====
//  Routine to gather particles' potentials from all (slave) ranks
//  to master rank.  Currently not used.
//  =======================================================================
void gather_potential(
          double potential[],                     //  [out]  array of particle potentials
          double computed_potential[],            //  [in]  array of potentials computed by FMM
          Cell ** cell_tree,                      //  [in]  tree of cells
          int myrank,                             //  [in]  rank
          int n_active,                           //  [in]  number of active ranks
          int rank_finest_level,                  //  [in]  rank finest level
          int finest_level,                       //  [in]  finest level (bottom level of the tree hierarchy)
          MPI_Comm fmm_world = MPI_COMM_WORLD     //  [in]  communicator
        )
{
  int depth = finest_level - rank_finest_level;

//  MPI variables
  MPI_Status status;
  MPI_Request request;


//  packing/unpacking variables (for (particle_id, potential) packets transmission)
  char * buffer;
  int   buffer_size;
  int   buffer_position;


  // allocating buffer for packing/unpacking
  buffer_size = sizeof(int)+sizeof(double);
  buffer = (char*) malloc (buffer_size);

  //  computational variables
  int     n_particles;
  int     particle_id;
  double  pot;

  if (myrank == 0) {

    // loop over cells
    for (ucellindex_t icell = 0; icell<(1<<(3*depth)); icell++) {
      ucellindex_t cell_idx = icell;
      Cell & cell = cell_tree[finest_level][cell_idx];
      n_particles = cell.particles.size();

      // loop over particles in the cell
      for (int ip = 0; ip < n_particles; ip++) {
        particle_id = cell.particles[ip];
        pot = computed_potential[particle_id];
        potential[particle_id]=pot;
      }
    }


    // loop over ranks
    for (int rk = 1; rk<n_active; rk++) {

      // loop over cells in each rank
      for (ucellindex_t icell = 0; icell<(1<<(3*depth)); icell++) {
        int n_particles;

        //  receive number of particles in the cell
        int tag = 0; // TODO:  replace by meaningful number
        MPI_Recv(&n_particles, 1, MPI_INT, rk, tag, fmm_world , &status);

        //  loop over particles to receive potential data (particle_id, potential)
        for (int ip = 0; ip < n_particles; ip++) {

          //  Receiving potential data (particle_id, potential)
          int particle_tag = 0; // TODO: replace by a meaningful number
          MPI_Recv(buffer, buffer_size, MPI_PACKED, rk, particle_tag, fmm_world , &status);

          //  Unpacking
          buffer_position = 0;
          MPI_Unpack (buffer, buffer_size, &buffer_position, &particle_id, 1, MPI_INT,    fmm_world );
          MPI_Unpack (buffer, buffer_size, &buffer_position, &pot,         1, MPI_DOUBLE, fmm_world );

          //  Processing
          potential[particle_id] = pot;
        }   //  ip-loop over particles;   for (int ip = 0; ...)
      }   //  icell-loop over cells
    }   //  rk-loop over ranks
  }   //  if (rank == 0)

  else {

    // loop over cells of this rank
    for (ucellindex_t icell = 0; icell<(1<<(3*depth)); icell++) {
      ucellindex_t cell_idx = (myrank<<(depth*3)) + icell;
      Cell & cell = cell_tree[finest_level][cell_idx];
      n_particles = cell.particles.size();

      //  loop over particles in the cell
      for (int ip = 0; ip < n_particles; ip++) {
        particle_id = cell.particles[ip];
        pot = computed_potential[particle_id];

      //  Packing
        buffer_position = 0;
        MPI_Pack (&particle_id,   1, MPI_INT,    buffer, buffer_size, &buffer_position, fmm_world );
        MPI_Pack (&pot,           1, MPI_DOUBLE, buffer, buffer_size, &buffer_position, fmm_world );

      //  Sending
        int tag = 0; // TODO:  some meaningful value for the tag
        MPI_Isend (buffer, buffer_size, MPI_PACKED, 0, tag, fmm_world , &request);
      //  Wait before re-using the buffer!!!
        MPI_Wait(&request, &status);

      }  //  ip-loop over particles
    }  //  icell-loop over cells

  }  // end of myrank > 0

  if (buffer) free (buffer);

}


//  ===================================================== MPI Parallel ====
//  Particle data (positions, charges) broadcast routine.
//  Master Rank:  Broadcasts input data
//  Slave Rank:   Allocates position/charges arrays, then receives
//                position/charges data
//  Used by:      main()
//  =======================================================================
void    bcast_particle_data (
          int         myrank,                        //  [in]  rank
          int         n_active,                      //  [in]  number of active ranks
 unsigned long     &  n_particles,                   //  [in/out]  number of particles
          double * &  position,                      //  [in/out]  particle positions, in triples (x, y, z)
          double * &  charges,                       //  [in/out]  charges of the particles
          MPI_Comm fmm_world = MPI_COMM_WORLD        //  [in]  Communicator
        )
{

  MPI_Bcast(&n_particles, 1, MPI_LONG, 0 /* ROOT_RANK */, fmm_world );
#if FMM_VERBOSE
  fprintf (stdout, "rank %d :: n_particles = %ld\n", myrank, n_particles); fflush(stdout);
#endif
//  Slave Ranks :: allocating memory for positions/charges data
  if (myrank != 0) {
    position = new double [3*n_particles];
    charges  = new double [n_particles];
  }

//  Broadcasting particles positions, charges from the master rank
  MPI_Bcast(position, 3*n_particles, MPI_DOUBLE, 0 /*ROOT_RANK*/, fmm_world );
  MPI_Bcast(charges,  n_particles,   MPI_DOUBLE, 0 /*ROOT_RANK*/, fmm_world );

//  Release positions/charges memory for the inactive ranks (ranks which do not participate in computations)
  if (myrank >= n_active) {
    delete [] position;
    delete [] charges;
  }
}
