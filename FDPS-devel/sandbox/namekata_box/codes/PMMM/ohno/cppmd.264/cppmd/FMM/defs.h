/*************************************************************
 *                                                           *
 *        FMM Library for CPPMD  MD simulation program.      *
 *                                                           *
 *        Created by:  Novichkov Gleb, PhD at ERI, Inc.      *
 *        Date:        2011.01 - 2011.05                     *
 *                                                           *
 *       ~~~  Common basic data types; definitions ~~~       *
 *                                                           *
 *************************************************************/
#ifndef  FMMDEFS_H
#define  FMMDEFS_H


//
//  Basis data types used in FMM program
//
enum  BOUNDARY_CONDITION {ZERO, PBC};

#ifdef HAVE_INT128

  typedef __uint128_t           uint128;
  typedef __int128_t            int128;
  typedef __uint128_t ucellindex_t;
  typedef __int128_t cellindex_t;
#define MPI_CELLINDEX MPI_DOUBLE_COMPLEX

#else

  typedef unsigned long long ucellindex_t;
  typedef long long cellindex_t;
#define MPI_CELLINDEX MPI_INTEGER8

#endif


//
//  Definitions used in FMM program
//
#define  DEFAULT_BOUNDARY_CONDITION     ZERO

#define  INDEX(m,n)  (n)*(n) + (n) + (m)
#define  SUPERINDEX(k,j,m,n,p)  (INDEX(k,j))*((p)+1)*((p)+1) + INDEX(m,n)

#define  RANK(cell_idx, depth)  (cell_idx) >> (3*(depth))

#ifndef NULL
  #define NULL (void *)0
#endif


#define  FMM_OMP_CHUNK    2

#endif  // FMMDEFS_H
