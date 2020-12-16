#pragma once

//* Constants
#define N_CHILDREN (8)

#ifdef USE_96BIT_KEY
#define TREE_LEVEL_LIMIT (31)
#else
#define TREE_LEVEL_LIMIT (21)
#endif

//* Basic data types
typedef int S32_;
typedef unsigned int U32_;
#ifdef PARTICLE_SIMULATOR_ALL_64BIT_PRECISION
typedef double F32_;
#else
typedef float F32_;
#endif
typedef long long int S64_;
typedef unsigned long long int U64_;
typedef double F64_;

//* Vector types
typedef struct {
   S32_ x,y,z;
} S32vec_;
typedef struct {
   F32_ x,y,z;
} F32vec_;
typedef struct {
   F64_ x,y,z;
} F64vec_;

//* Symmetric matrix types
typedef struct {
   F32_ xx, yy, zz, xy, xz, yz;
} F32mat_;
typedef struct {
   F64_ xx, yy, zz, xy, xz, yz;
} F64mat_;

//* Orthotope types
typedef struct {
   F32vec_ low_;
   F32vec_ high_;
} F32ort_;
typedef struct {
   F64vec_ low_;
   F64vec_ high_;
} F64ort_;

/*
#ifdef USE_96BIT_KEY
typedef struct {
    U64_ up_;
    U32_ lo_;
} KeyT_;
#else
typedef struct {
    U64_ up_;
} KeyT_;
#endif
*/
