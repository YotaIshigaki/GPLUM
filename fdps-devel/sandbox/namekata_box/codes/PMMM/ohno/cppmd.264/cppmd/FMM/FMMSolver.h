/*************************************************************
 *                                                           *
 *  FMM Library for CPPMD  MD simulation program.            *
 *                                                           *
 *  Created by:  Novichkov Gleb, PhD at ERI, Inc             *
 *               Ohno Yousuke, PhD at RIKEN                  *
 *                                                           *
 *  Date:        2011.01 - 2011.05                           *
 *                                                           *
 *         ~~~  FMM Solver Class. Header File.  ~~~          *
 *                                                           *
 *************************************************************/

#ifndef  FMMSOLVER_H
#define  FMMSOLVER_H
#include <mpi.h>

#include "defs.h"                         //  FMM definitions file
#include "cell.h"                         //  Cell class header file
#include "coefficients.h"                 //  Precomputed coefficients class header

#ifdef ERI

#include "ParticleInfo.h"                 //  CPPMD Header file.  Located in `cppmd/Common'
#else

#include <ParticleInfo.h>                 //  CPPMD Header file.  Located in `cppmd/Common'
#endif



//
//    Sequential FMM Solver
//
class FMMSolver {

protected:
  Cell ** cell_tree;
  int finest_level;           //  Finest Level
  int p;                      //  Order of Multipole Expansion


  unsigned long _total;
  double * _positions;
  double * _charges;

  double scale;
  double fscale;

public:

  Coefficients * coeffs;
  std::complex <double>  * spherical_harmonics;                    //  Storage for the spherical harmonics.
                                                              //  to avoid allocation/deallocation in co

  FMMSolver();
  FMMSolver (int L, int expansion_order);
  virtual ~FMMSolver();

  virtual   void  setParticles (
                      unsigned long  total,                   // [in]
                      double positions     [],                // [in]
                      double charges       []                 // [in]
                  );

  virtual   void  setParticles (
                      ParticleArray & particles               // [in]
                  );


//  virtual   void  setParticles () = 0;                       // Implemented in FMMSolverParallel class


  virtual   void  RunSolver (
#ifdef  USE_POTENTIAL_ARRAY
                      double potential     [],                // [out]
#endif
                      double forces        [],                // [out]
                      double & total_potential,               // [out]
                      BOUNDARY_CONDITION  boudary_condition   // [in]
                      = DEFAULT_BOUNDARY_CONDITION
                  );
/*
  virtual   void  RunSolver (
                      ParticleArray & particles,              // [in]
                      double & total_potential,               // [out]
                      BOUNDARY_CONDITION  boudary_condition   // [in]
                      = DEFAULT_BOUNDARY_CONDITION
                  );
*/
  void destroyParticles();


  unsigned long get_total();

  void set_scale(const double _scale);

  void scale_force(double *force);

protected:
  void fmm_create_tree ();
  void fmm_clear_cell_tree_particles ();
  void fmm_clear_cell_tree_MECs ();
  void fmm_delete_tree ();
};


//
//    Parallel FMM Solver
//
class FMMSolverParallel : public FMMSolver {

  int rank_finest_level;      //  rank-finest level (for MPI)
  int numproc, myrank;        //  number of processes, current rank
  int n_active;               //  number of active (participating in FMM computation) ranks
  MPI_Comm myworld;           //  communicator

public:

  FMMSolverParallel (int L, int expansion_order, MPI_Comm lcomm = MPI_COMM_WORLD);
  virtual ~FMMSolverParallel();

//  Usage of setParticles():
//    master rank:  setParticle (total, position, charges)
//    slave rank:   setParticle ()
//
  virtual   void  setParticles (
                      unsigned long total = 0,                 // [in]
                      double * positions  = NULL,              // [in]
                      double * charges    = NULL               // [in]
                  );

  virtual   void  setParticles (
                      ParticleArray & particles                // [in]
                  );

  template<class PA, class GPA>
  void update_particles (       const PA& particlearray,
                                const std::vector<TypeRange>& typerange,
                                const std::vector<int>& self_longset_index,
                                const GPA& ghost,
                                const std::vector<TypeRange>& ghosttyperange,
                                const std::vector<int>& ghost_longset_index
                              );

  template<class PA, class GPA>
  unsigned long allocate_particles (     const PA& particlearray,
                                         const std::vector<TypeRange>& typerange,
                                         const std::vector<int>& self_longset_index,
                                         const GPA& ghost,
                                         const std::vector<TypeRange>& ghosttyperange,
                                         const std::vector<int>& ghost_longset_index
                                         );


  virtual   void  RunSolver (
#ifdef  USE_POTENTIAL_ARRAY
                      double potential [],                    // [out]
#endif
                      double forces [],                       // [out]
                      double & total_potential,               // [out]
                      BOUNDARY_CONDITION  boudary_condition   // [in]
                      = DEFAULT_BOUNDARY_CONDITION
                  );

  int   get_n_active ();
/*
  virtual   void  RunSolver (
                      ParticleArray & particles,              // [in]
                      double & total_potential,               // [out]
                      BOUNDARY_CONDITION  boudary_condition   // [in]
                      = DEFAULT_BOUNDARY_CONDITION
                  );
*/



};


#endif
