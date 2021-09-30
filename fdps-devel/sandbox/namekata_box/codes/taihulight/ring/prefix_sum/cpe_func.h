#pragma once
//================
/* Header files */
//================
/* Standard C headers */
/* Headers that are specific for Sunway-Taihulight */
/* User-defined headers */
#include "ps_defs_cpe.h"


//====================================================
/* Macros and macro functions used in CPE functions */
//====================================================
//** the macros used in copy.c
#define NUMBER_OF_CPE (64)
#define LOCAL_MEMORY_SIZE_LIMIT (40960)
//#define LOCAL_MEMORY_SIZE_LIMIT (32768)
//#define LOCAL_MEMORY_SIZE_LIMIT (16384)
#define ELEMENT_SIZE_LIMIT (2048)
//----------------------------------------------------

//** the macros used in Long's original code kernel_sunway_simd.h
//#define N_WALK_LIMIT 1280
#define NI_LIMIT     128
#define NSAT_LIMIT   132
#define LDM_B_LIMIT  50000
//#define N_CPE_LIMIT  64

//#define NI_LIMIT_MM  262144
//#define N_EPJ_LIMIT	 20000000
//#define N_SPJ_LIMIT	 10000000

#define EPI_SIZE 64
#define EPJ_SIZE 64
#define SPJ_SIZE 32
#define FORCE_SIZE 128
#define FORCEMM_SIZE 32

#define REG_PUTR(var, dst) asm volatile ("putr %0,%1\n"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1\n"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0\n":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0\n":"=r"(var))
#define NOP() asm volatile ("nop\n")
//----------------------------------------------------


//=========================
/* Structure definitions */
//=========================
typedef double REAL;

typedef struct{
    REAL r[4]; //x,y,z,m
    REAL v[4]; //vx,vy,vz,id
} EpiMM;
typedef EpiMM EpjMM;

typedef struct{
	REAL r[4];
} SpjMM;

typedef struct{
  REAL ax[4],ay[4],az[4],pot[4];
} ForceMM4;

typedef struct{
  REAL ax,ay,az,pot;
} ForceMM;

typedef struct{
  EpiMM* epi;
  EpjMM* epj;
  SpjMM* spj;
  EpiMM* sat;
  ForceMM* force_sat;
  int*   n_disp_epi;
  int*   n_disp_epj;
  int*   n_disp_spj;
  int*   n_epi;
  int*   n_epj;
  int*   n_spj;
  int*   adr_epj;
  int*   adr_spj;
  
  //double r_coll; 
  double r_ep;
  double r_sat;
  double kappa; // k^{'}
  double eta; // \eta^{'}
  double m_particle;
  double m_satelite;
  double m_planet;
  //double eps2;
  double dt;
  double ekin;
  double pot;
  int    n_walk;
  int    n_sat;
  int    energy_flag;
  int    first_flag;
  int    last_flag;

#ifdef SUNWAY_FORCE_KERNEL_CHECK
  ForceMM* forcei;
#endif

} Force_Kernel_Pars;


//* EssentialParticleJ (used in CalcMoment, CalcMomentLongGlobalTree, GenMortonKey)
typedef struct {
   F64vec_ pos;
   F64_ mass;
   F64vec_ vel;
   S64_ id;
} epjLM;

//* SuperParticleJ (used in CalcMoment, CalcMomentLongGlobalTree, GenMortonKey)
typedef struct {
   F64vec_ pos;
   F64_ mass;
#ifdef USE_QUADRUPOLE
   F32mat_ quad;
#endif
} spjLM;

//* TreeParticle (the order of vars. must be consistent with TreeParticle class in FDPS)
typedef struct {
   U64_ key_;
   U32_ adr_ptcl_;
   int pad_;
} tpLM;

//* TreeCell (the order of vars must be consistent with TreeCell class in FDPS)
typedef struct {
   S32_ n_ptcl_;
   U32_ adr_tc_;
   U32_ adr_ptcl_;
   S32_ level_;
   F32_ mass;
   F32vec_ pos;
#ifdef USE_QUADRUPOLE
   F32mat_ quad;
#endif
   F64ort_ vertex_out_;
   F64ort_ vertex_in_;
} tcLM;


//==========================
/* Prototype declarations */
//==========================
void CopyIndirect(void *);
void CopyIndirectInverse(void *);
void CopyIndirectInverse2(void *);
//void CopyIndirect2(void *);
void CopyDirect(void *);
void CopyStride(void *);
void CopyStrideByHand(void *);
void GenMortonKey(void *);
void CalcMoment(void *);
void CalcMomentLongGlobalTree(void *);
void ForceKernelSunWay1st(void* pars);
void ForceKernelSunWay2nd(void* pars);

