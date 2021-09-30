/*************************************************************
 *                                                           *
 *  FMM Library for CPPMD  MD simulation program.            *
 *                                                           *
 *  Created by:  Novichkov Gleb, PhD at ERI, Inc.            *
 *  Date:        2011.01 - 2011.05                           *
 *                                                           *
 *   ~~~  Routines to contruct tree hierarchy cells.  ~~~    *
 *                                                           *
 *************************************************************/

#include "defs.h"                       //  FMM definitions file
#include "cell.h"                       //  Cell class header file

#include <stdio.h>                      //  for printf()
#include <string.h>                     //  for strstr()
#include <stdlib.h>                     //  for atoi()
#include <time.h>                       //  for time(), etc.
#include <math.h>                       //  for abs()
#include <vector>                       //  for vector class


#define test_bit(var, pos)  !!(((ucellindex_t)var) & ((ucellindex_t)1<< (pos)))

using namespace std;

inline ucellindex_t ciabs(cellindex_t val) { return val >=0 ? val : -val; }


extern "C++" {

//  output routines (for debug)
  void print_list (vector<ucellindex_t> list);
  void print_point (float point[], int dimension);
}


//
//  ***  Service Routines for Cell Tree Hierarchy  ***
//

/*********************************************************************************
 *      Purpose: computes position of the cell in the i-th direction.
 *      Returns: value in the range 0 .. 2^dimension - 1
 *      Type:    A service routine for the tree construction.
 **/
ucellindex_t project_i (
          const ucellindex_t & cell_index,         //  [in]  index of the cell
          const int & directon,                    //  [in]  the direction
          const int & level,                       //  [in]  level of refinement
          const int & dimension = 3                //  [in]  dimension of the system
        )
{
  ucellindex_t cell_pos_i = 0;
  for (int k = 1; k<=level; k++) {
    cell_pos_i += ucellindex_t (test_bit (cell_index, dimension*(level-k) + directon)) << (level-k);
  }
  return cell_pos_i;
}

/*********************************************************************************
 *      Purpose:  Finds neighbours of cell in the i-th direction
 *      Returns:  A list of positions of the neighbors, in the i-th direction
 *      Type:     A service routine for the tree construction.
 **/
vector <ucellindex_t>
find_neighbours_i (
          const ucellindex_t & cell_index,         //  [in]  index of the cell
          const int & direction,                   //  [in]  the direction
          const int & level,                       //  [in]  level of refinement
          const int & dimension,                   //  [in]  dimension of the system
          const BOUNDARY_CONDITION &
          boundary_condition =
          DEFAULT_BOUNDARY_CONDITION               //  [in]  boundary condition
        )
{

  vector <ucellindex_t> neighbours;

  //  Finding Cell Position in the i-direction
  ucellindex_t cell_pos_i = project_i (cell_index, direction, level, dimension);
  if (boundary_condition == PBC) {

    cellindex_t left_cell_i  = (cellindex_t (cell_pos_i) - 1 ) & ((1<<level) - 1);
    cellindex_t right_cell_i = (cellindex_t (cell_pos_i) + 1 ) & ((1<<level) - 1);

    neighbours.push_back (left_cell_i);
    neighbours.push_back (cell_pos_i);
    neighbours.push_back (right_cell_i);
  }
  else if (boundary_condition == ZERO) {
//  ========= left-most cell
    if (cell_pos_i == 0) {
      neighbours.push_back (cell_pos_i);
      neighbours.push_back (cell_pos_i+1);
    }
    else
//  ========= right-most cell
    if (cell_pos_i == (1<<level) - 1) {
      neighbours.push_back (cell_pos_i);
      neighbours.push_back (cell_pos_i-1);
    }
//  ========= interior cell
    else {
      neighbours.push_back (cell_pos_i-1);
      neighbours.push_back (cell_pos_i);
      neighbours.push_back (cell_pos_i+1);
    }
  }
  return neighbours;
}
vector <ucellindex_t>
find_neighbours_i_pbc (
          const ucellindex_t & cell_index,         //  [in]  index of the cell
          const int & direction,                   //  [in]  the direction
          const int & level,                       //  [in]  level of refinement
          const int & dimension,                   //  [in]  dimension of the system
	  vector<double> & shift          // [out] offset PBC 
        )
{

  vector <ucellindex_t> neighbours;

  //  Finding Cell Position in the i-direction
  ucellindex_t cell_pos_i = project_i (cell_index, direction, level, dimension);
  {

    cellindex_t left_cell_i  = (cellindex_t (cell_pos_i) - 1 ) & ((1<<level) - 1);
    cellindex_t right_cell_i = (cellindex_t (cell_pos_i) + 1 ) & ((1<<level) - 1);
    
    double left = ((cell_pos_i == 0) ? -1.0 : 0.0);
    double right = ((cell_pos_i == (1<<level) - 1) ? 1.0 : 0.0);

    neighbours.push_back (left_cell_i);
    neighbours.push_back (cell_pos_i);
    neighbours.push_back (right_cell_i);
    shift.push_back(left);
    shift.push_back(double(0.0));
    shift.push_back(right);
  }
  return neighbours;
}

/*********************************************************************************
 *      Purpose:  Computes partial indices (PI) of the neighbours from their positions,
 *                in the i-th direction
 *      Returns:  A list of partial indices of the cells,  the i-th direction
 *      Type:     A service routine for the tree construction.
 **/
vector <ucellindex_t>
compute_PI (
    const vector <ucellindex_t> & neighbours_list,   //  [in]  list of neighbours positions, in the i-th direction
    const int & direction,                           //  [in]  the direction
    const int & level,                               //  [in]  level of refinement
    const int & dimension                            //  [in]  the dimension
  )
{
  vector <ucellindex_t> PI_list;

  for (int i = 0; i<neighbours_list.size(); i++) {

    ucellindex_t cell_pos_i = neighbours_list[i];
    ucellindex_t p_i = 0;

    for (int k = 1; k<=level; k++)
      p_i += ucellindex_t (test_bit (cell_pos_i, level - k)) << (dimension*(level-k) + direction);

    PI_list.push_back (p_i);

  }

  return PI_list;
}



/*********************************************************************************
 *      Purpose:  Merges (bitwise) partial indices (PI) of the cells
 *      Returns:  A list of merged partial indices of the cells.
 *                the size of the list is |list1| x |list2|.
 *      Type:     A service routine for the tree construction.
 **/
vector<ucellindex_t>
merge_PI (
      const vector <ucellindex_t> & list1,         //  [in]  first list
      const vector <ucellindex_t> & list2          //  [in]  second list
    )
{
  vector <ucellindex_t> result_list;

  if (list2.size()>0)
    for (int i = 0; i<list1.size(); i++)
    for (int j = 0; j<list2.size(); j++) {
      ucellindex_t p_i = list1[i] + list2[j];
      result_list.push_back(p_i);
    }
  else
    for (int i = 0; i<list1.size(); i++){
      ucellindex_t p_i = list1[i];
      result_list.push_back(p_i);
    }
  return result_list;
}

vector<ucellindex_t>
merge_PI_pbc (
      const vector <ucellindex_t> & list1,         //  [in]  first list
      const vector <ucellindex_t> & list2,         //  [in]  second list
      const int                   & direction,
      const int                   & dimension,
      const vector <double>       & shift1D,
      vector <double>             & shift
    )
{
  vector <ucellindex_t> result_list;
  vector <double> result_shift;

  if (list2.size()>0)
    for (int i = 0; i<list1.size(); i++)
    for (int j = 0; j<list2.size(); j++) {
      ucellindex_t p_i = list1[i] + list2[j];
      result_list.push_back(p_i);
      for(int d=0;d<direction;d++)result_shift.push_back(shift[j*dimension+d]);
      result_shift.push_back(shift1D[i]+shift[j*dimension+direction]);
      for(int d=direction+1;d<dimension;d++)result_shift.push_back(shift[j*dimension+d]);
    }
  else{
    for (int i = 0; i<list1.size(); i++){
      ucellindex_t p_i = list1[i];
      result_list.push_back(p_i);
      for(int d=0;d<direction;d++)result_shift.push_back(double(0.0));
      result_shift.push_back(shift1D[i]);
      for(int d=direction+1;d<dimension;d++)result_shift.push_back(double(0.0));
    }
  }
  shift.swap(result_shift);
  return result_list;
}

/*********************************************************************************
 *      Purpose:  Determines if two cell are well separated
 *      Returns:  1, if two cells are well separate, 0 -- otherwise
 *      Type:     A service routine for the tree construction.
 **/
bool is_well_separated (
        const ucellindex_t &  cell1_index,         //  [in]  index of the first cell
        const ucellindex_t &  cell2_index,         //  [in]  index of the second cell
        const int  &          level,               //  [in]  level of refinement
        const int  &          dimension,           //  [in]  dimension
        const BOUNDARY_CONDITION &
              boundary_condition =
              DEFAULT_BOUNDARY_CONDITION           //  [in]  boundary condition
      )
{
  bool is_near = true;
  ucellindex_t cell1_pos;
  ucellindex_t cell2_pos;
  ucellindex_t distance;

  for (int i = 0; i<dimension; i++) {

    cell1_pos = project_i (cell1_index, i, level, dimension);
    cell2_pos = project_i (cell2_index, i, level, dimension);
    distance = ciabs(cellindex_t(cell1_pos - cell2_pos));
    if (boundary_condition == PBC) {
    ////////// TODO  Level 2  distance is depend replica or original
      if (distance == ((1<<level) - 1))
        distance = 1;

/*
   //
   // if nearest neighbors include of distance 2.  To be checked and tested.
   //
      ucellindex_t n = ucellindex_t (distance >> (level - 1));
      distance  = distance - (n<<(level - 1));
*/
    }
    if (distance>1) {
      is_near = false;
      break;
    }

  }

  return !is_near;
}

bool is_well_separated_pbc (
        const ucellindex_t &  cell1_index,         //  [in]  index of the first cell
        const ucellindex_t &  cell2_index,         //  [in]  index of the second cell
        const int  &          level,               //  [in]  level of refinement
        const int  &          dimension,           //  [in]  dimension
        const double shift[]
      )
{
  bool is_near = true;
  ucellindex_t cell1_pos;
  ucellindex_t cell2_pos;
  ucellindex_t distance;

  for (int i = 0; i<dimension; i++) {

    cell1_pos = project_i (cell1_index, i, level, dimension);
    cell2_pos = project_i (cell2_index, i, level, dimension);
    distance = ciabs(cellindex_t(cell1_pos - cell2_pos + int(shift[i]*(1<<level))));

    if (distance>1) {
      is_near = false;
      break;
    }

  }

  return !is_near;
}

/*********************************************************************************
 *      Purpose:  Computes the local index of the child (inside parent cell)
 *      Returns:  Local index of the child
 *      Type:     A service routine for the tree construction.
 */
ucellindex_t   compute_local_index (
            const ucellindex_t &  index,           //  [in]
            const int &  level,                    //  [in]
            const int &  finest_level,             //  [in]
            const int &  dimension = 3             //  [in]
          )
{
  ucellindex_t   local_index;

  //  divide by 2^(d*(L_max - L)) and take the last d bits
  local_index = index >> (dimension*(finest_level - level)) & ((1<<dimension) - 1);
  return local_index;
}


//
//  ***  Major Routines for Cell Tree Hierarchy  ***
//
/*********************************************************************************
 *      Purpose:  Finds nearest neighbors of the cell
 *      Returns:  A list of neighbors (their indices)
 *      Type:     A major routine
 **/
vector<ucellindex_t>
find_nearest_neighbors(
    const ucellindex_t  &   cell_index,            //  [in]  index of the cell
    const int           &   level,                 //  [in]  level of refinement
    const int           &   dimension,             //  [in]  dimension
    const BOUNDARY_CONDITION
          boundary_condition =
          DEFAULT_BOUNDARY_CONDITION               //  [in]  boundary condition
  )
{

  vector<ucellindex_t> neighbors_list;

  for (int i = 0 ; i<dimension; i++) {

//  ======= find neighbors in the direction i
    vector <ucellindex_t> neighbors_i = find_neighbours_i (cell_index, i, level, dimension, boundary_condition);

//  ======= Partial Index (PI) Computation
//    1. compute PI^i for each neighbour
//    2. put computed value in the temporary list L_tmp
    vector <ucellindex_t> L_tmp = compute_PI (neighbors_i, i, level, dimension);

//  ======= Merge PIs
//  step 3 :: Using list L_{i-1} of PIs obtained in the previous steps, create
//            a new list by adding each element of L_tmp to each element of L_{i-1}
//            Obtained new list has size |L_{i-1}|x|L_tmp|.
    neighbors_list = merge_PI (L_tmp, neighbors_list);
  }

  return neighbors_list;
}

vector<ucellindex_t>
find_nearest_neighbors_pbc(
    const ucellindex_t  &   cell_index,            //  [in]  index of the cell
    const int           &   level,                 //  [in]  level of refinement
    const int           &   dimension,             //  [in]  dimension
    vector<double>      &   shift                  //  [out]
  )
{

  vector<ucellindex_t> neighbors_list;

  for (int i = 0 ; i<dimension; i++) {

//  ======= find neighbors in the direction i
    vector <double> shift1D;
    vector <ucellindex_t> neighbors_i = find_neighbours_i_pbc (cell_index, i, level, dimension, shift1D);

//  ======= Partial Index (PI) Computation
//    1. compute PI^i for each neighbour
//    2. put computed value in the temporary list L_tmp
    vector <ucellindex_t> L_tmp = compute_PI (neighbors_i, i, level, dimension);

//  ======= Merge PIs
//  step 3 :: Using list L_{i-1} of PIs obtained in the previous steps, create
//            a new list by adding each element of L_tmp to each element of L_{i-1}
//            Obtained new list has size |L_{i-1}|x|L_tmp|.
    neighbors_list = merge_PI_pbc (L_tmp, neighbors_list, i, dimension, shift1D, shift);
  }

  return neighbors_list;
}





/*********************************************************************************
 *      Purpose:  Computes the index of the cell that contains the point
 *      Returns:  Index of the cell
 *      Type:     A major routine
 **/
ucellindex_t  compute_point_index (
            const double point[],                  //  [in]  position of the point in R^d
            const int &  level,                    //  [in]  level of refinement
            const int &  dimension                 //  [in]  dimension
          )
{

  ucellindex_t cell_pos_i;
  ucellindex_t index;

  index = 0;
  for (int l = 1; l<=level; l++) {

    int local_index = 0;
    for (int i = 0; i<dimension; i++) {
      double pos = point[i];
      //      if(pos>=1.0)pos-=1.0;   // TODO: replace by correct version.
      //      if(pos<0.0)pos+=1.0;
      cell_pos_i   = ucellindex_t ((1<<l) * pos);
      local_index += (cell_pos_i & 1)<<i;
    }

    index += ucellindex_t(local_index)<<(dimension*(level - l));
  }
  return index;
}

/*********************************************************************************
 *      Purpose: computes center of the cell
 *      Returns: center of the cell, as array double []
 *      Type:    routine necessary for the upward/downward passes
 **/
void compute_cell_center (
          double center[],                         //  [out] center of the cell
          const ucellindex_t & cell_index,         //  [in]  index of the cell
          const int & level,                       //  [in]  level of refinement
          const int & dimension                    //  [in]  dimension of the system
        )
{
  ucellindex_t cell_pos_i, div;

  for (int i = 0; i<dimension; i++) {
    cell_pos_i = project_i (cell_index, i, level, dimension);
    div = 1 << (level + 1);
    center[i] = (double((cell_pos_i<<1) + 1))/double(div);
  }
}


//
//  *************************************************************  Interaction List Compuation Routines  *********
//


/*********************************************************************************
 *      Purpose:  Computes the interaction list of a cell
 *      Returns:  Interaction list of the cell (indices of the cells on the list).
 *      Type:     A major routine
 **/
void find_interaction_list (
    vector<ucellindex_t> &  interaction_list,   //  [out] interaction list for the cell `cell_index'
    const ucellindex_t &    cell_index,         //  [in]  cell's index
    const int          &    level,              //  [in]  level of refinement
    const BOUNDARY_CONDITION &
          boundary_condition =
          DEFAULT_BOUNDARY_CONDITION,           //  [in]  boundary condition
    const int          &    dimension = 3       //  [in]  dimension
  )
{
  vector<ucellindex_t>  nearest_parents;
  ucellindex_t  parent_index;

//  ======= computing index of the cell's parent
  parent_index = cell_index >> dimension;

//  ======= finding nearest neighbors of the parent
  nearest_parents = find_nearest_neighbors (parent_index, level-1, dimension, boundary_condition);


//  ======= finding well separated children of the parent's neighbors
//    for each child of the each neighbour of the parent,
//      1. check, if it near
//      2. if yes, add to the interaction list;
//
  for ( ; nearest_parents.size()>0 ; ) {

    ucellindex_t parent = nearest_parents.back();
    nearest_parents.pop_back();

    for (int i = 0; i<(1<<dimension); i++) {
      ucellindex_t child = (parent<<dimension) + i;
      if (is_well_separated (child, cell_index, level, dimension, boundary_condition))
        interaction_list.push_back (child);
    }

  }
//  ======= done
}

/*********************************************************************************
 *      Purpose:  Computes the interaction list of a cell
 *      Returns:  Interaction list of the cell (indices of the cells on the list).
 *      Type:     A major routine
 **/

vector<ucellindex_t>  find_interaction_list (
    const ucellindex_t &    cell_index,         //  [in]  cell's index
    const int          &    level,              //  [in]  level of refinement
    const BOUNDARY_CONDITION &
          boundary_condition =
          DEFAULT_BOUNDARY_CONDITION,           //  [in]  boundary condition
    const int          &  dimension = 3         //  [in]  dimension
  )
{
  vector<ucellindex_t>  interaction_list;
  vector<ucellindex_t>  nearest_parents;
  ucellindex_t  parent_index;

//  ======= computing index of the cell's parent
  parent_index = cell_index >> dimension;

//  ======= finding nearest neighbors of the parent
  nearest_parents = find_nearest_neighbors (parent_index, level-1, dimension, boundary_condition);


//  ======= finding well separated children of the parent's neighbors
//    for each child of the each neighbour of the parent,
//      1. check, if it near
//      2. if yes, add to the interaction list;
//
  for ( ; nearest_parents.size()>0 ; ) {

    ucellindex_t parent = nearest_parents.back();
    nearest_parents.pop_back();

    for (int i = 0; i<(1<<dimension); i++) {
      ucellindex_t child = (parent<<dimension) + i;
      if (is_well_separated (child, cell_index, level, dimension, boundary_condition))
        interaction_list.push_back (child);
    }

  }
//  ======= done

  return interaction_list;
}

vector<ucellindex_t>  find_interaction_list_pbc (
						 vector<double> & shift,
    const ucellindex_t &    cell_index,         //  [in]  cell's index
    const int          &    level,              //  [in]  level of refinement
    const int          &  dimension = 3         //  [in]  dimension
  )
{
  vector<ucellindex_t>  interaction_list;
  vector<ucellindex_t>  nearest_parents;
  ucellindex_t  parent_index;
  double *shift_tmp;
  
  shift_tmp = new double[dimension];

//  ======= computing index of the cell's parent
  parent_index = cell_index >> dimension;

//  ======= finding nearest neighbors of the parent
  vector<double>  parent_shift;
  nearest_parents = find_nearest_neighbors_pbc (parent_index, level-1, dimension, parent_shift);


//  ======= finding well separated children of the parent's neighbors
//    for each child of the each neighbour of the parent,
//      1. check, if it near
//      2. if yes, add to the interaction list;
//
  for ( ; nearest_parents.size()>0 ; ) {

    ucellindex_t parent = nearest_parents.back();
    nearest_parents.pop_back();
    for(int d=dimension-1;d>=0;d--){
      shift_tmp[d] = parent_shift.back();
      parent_shift.pop_back();
    }
    
    for (int i = 0; i<(1<<dimension); i++) {
      ucellindex_t child = (parent<<dimension) + i;
      if (is_well_separated_pbc (child, cell_index, level, dimension, shift_tmp)){
        interaction_list.push_back (child);
	for(int d=0;d<dimension;d++){
	  shift.push_back(shift_tmp[d]);
	}
      }
    }

  }
//  ======= done

  return interaction_list;
}



/*********************************************************************************
 *      Purpose:  Computes the restricted (to `rank') interaction list of a cell
 *      Returns:  Interaction list of the cell (indices of the cells on the list),
 *                which belong to the rank `rank'
 *      Type:     A major routine
 **/
void find_restricted_interaction_list(
    vector<ucellindex_t> & interaction_list,   //  [out] interaction list restricted to `rank'
    const ucellindex_t   &  cell_index,         //  [in]  cell's index
    const int &             level,              //  [in]  level of refinement
    const int &             rank_finest_level,  //  [in]  rank-finest level
    const int &             rank,               //  [in]  rank
    const BOUNDARY_CONDITION &
          boundary_condition =
          DEFAULT_BOUNDARY_CONDITION,           //  [in]  boundary condition
    const int &             dimension = 3       //  [in]  dimension
  )
{
  vector<ucellindex_t>  nearest_parents;
  ucellindex_t  parent_index;


//  ======= computing index of the cell's parent
  parent_index = cell_index >> dimension;

//  ======= finding nearest neighbors of the parent
  nearest_parents = find_nearest_neighbors (parent_index, level-1, dimension, boundary_condition);


//  `depth' is need to determine if `rank' owns the cell
  int depth = level - rank_finest_level;


//  ======= finding well separated children of the parent's neighbors
//    for each child of the each neighbour of the parent,
//      1. check, if it near
//      2. if yes, add to the interaction list;
//
  for ( ; nearest_parents.size()>0 ; ) {

    ucellindex_t parent = nearest_parents.back();
    nearest_parents.pop_back();

    for (int i = 0; i<(1<<dimension); i++) {
      ucellindex_t child = (parent<<dimension) + i;
      if (RANK(child, depth) == rank)
        if ( is_well_separated (child, cell_index, level, dimension, boundary_condition)) {
          interaction_list.push_back (child);
        }
    }

  }
//  ======= done

}

void find_restricted_interaction_list_pbc(
    vector<ucellindex_t> & interaction_list,   //  [out] interaction list restricted to `rank'
    vector<double>       &   shift,
    const ucellindex_t   &  cell_index,         //  [in]  cell's index
    const int &             level,              //  [in]  level of refinement
    const int &             rank_finest_level,  //  [in]  rank-finest level
    const int &             rank,               //  [in]  rank
    const int &             dimension = 3       //  [in]  dimension
  )
{
  vector<ucellindex_t>  nearest_parents;
  ucellindex_t  parent_index;
  double *shift_tmp;
  
  shift_tmp = new double[dimension];


//  ======= computing index of the cell's parent
  parent_index = cell_index >> dimension;

//  ======= finding nearest neighbors of the parent
  vector<double>  parent_shift;
  nearest_parents = find_nearest_neighbors_pbc (parent_index, level-1, dimension, parent_shift);


//  `depth' is need to determine if `rank' owns the cell
  int depth = level - rank_finest_level;


//  ======= finding well separated children of the parent's neighbors
//    for each child of the each neighbour of the parent,
//      1. check, if it near
//      2. if yes, add to the interaction list;
//
  for ( ; nearest_parents.size()>0 ; ) {

    ucellindex_t parent = nearest_parents.back();
    nearest_parents.pop_back();
    for(int d=dimension-1;d>=0;d--){
      shift_tmp[d] = parent_shift.back();
      parent_shift.pop_back();
    }

    for (int i = 0; i<(1<<dimension); i++) {
      ucellindex_t child = (parent<<dimension) + i;
      if (RANK(child, depth) == rank)
        if ( is_well_separated_pbc (child, cell_index, level, dimension, shift_tmp)) {
          interaction_list.push_back (child);
	  for(int d=0;d<dimension;d++){
	    shift.push_back(shift_tmp[d]);
	  }
        }
    }

  }
//  ======= done

}

/*********************************************************************************
 *      Computes:   a) interaction list of a cell restricted to `rank'
 *                  b) list of ranks != `rank' that own cells of the (full) interaction list
 *      Returns:    a) and b)
 *      Type:       A major routine
 **/
void find_restricted_interaction_list_ranks (
    vector<ucellindex_t> &  interaction_list,   //  [out] interaction list restricted to `rank'
    vector<int>          &  involved_ranks,     //  [out] list of ranks != `rank' that contain interaction list cells
    const ucellindex_t   &  cell_index,         //  [in]  cell's index
    const int &             level,              //  [in]  level of refinement
    const int &             rank_finest_level,  //  [in]  rank-finest level
    const int &             rank,               //  [in]  rank
    const int &             n_ranks,            //  [in]  number of ranks
    const BOUNDARY_CONDITION &
          boundary_condition =
          DEFAULT_BOUNDARY_CONDITION,           //  [in]  boundary condition
    const int &             dimension = 3       //  [in]  dimension
  )
{
//  vector<ucellindex_t>  interaction_list;
  vector<ucellindex_t>  nearest_parents;
  ucellindex_t  parent_index;

  int * involved_ranks_flag = new int [n_ranks];
  memset (involved_ranks_flag, 0, n_ranks*sizeof(int));

//  ======= computing index of the cell's parent
  parent_index = cell_index >> dimension;

//  ======= finding nearest neighbors of the parent
  nearest_parents = find_nearest_neighbors (parent_index, level-1, dimension, boundary_condition);


//  `depth' is need to determine if `rank' owns the cell
  int depth = level - rank_finest_level;
  int cell_owner_rank = RANK (cell_index, depth);

//  ======= finding well separated children of the parent's neighbors
//    for each child of the each neighbour of the parent,
//      1. check, if it near
//      2. if yes, add to the interaction list;
//
  for ( ; nearest_parents.size()>0 ; ) {

    ucellindex_t parent = nearest_parents.back();
    nearest_parents.pop_back();

    for (int i = 0; i<(1<<dimension); i++) {
      ucellindex_t child = (parent<<dimension) + i;

      if ( is_well_separated (child, cell_index, level, dimension, boundary_condition)) {
        int rk = RANK(child, depth);
        if (rk == rank /*cell_owner_rank rank*/)
          interaction_list.push_back (child);
        else
          involved_ranks_flag[rk] = 1;

      }
    }

  }
//  ======= done

  for (int ir = 0; ir<n_ranks; ir++) if (involved_ranks_flag[ir]) {
    involved_ranks.push_back (ir);
  }
  delete [] involved_ranks_flag;
}
