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
//** the macros used in copy_msort_glb_tree_reuse.c
//#define VERSION_MSORT_GLB_TREE_REUSE (2)
#define VERSION_MSORT_GLB_TREE_REUSE (0)
//----------------------------------------------------

//** the macros used in Long's original code kernel_sunway_simd.h
//#define N_WALK_LIMIT 1280
//#define NI_LIMIT     128
#define NI_LIMIT     256
#define NSAT_LIMIT   132
#define LDM_B_LIMIT  50000
//#define N_CPE_LIMIT  64

//#define NI_LIMIT_MM  262144
//#define N_EPJ_LIMIT	 20000000
//#define N_SPJ_LIMIT	 10000000

#define EPI_SIZE 64
#define EPJ_SIZE 64
#define SPJ_SIZE 32
#define TP_SIZE  16
#define FORCE_SIZE 128
#define FORCEMM_SIZE 32

#define NSPLIT 64
#define NTHREADS 64
#define OVERSAMPKLE 200


#define REG_PUTR(var, dst) asm volatile ("putr %0,%1\n"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1\n"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile ("getr %0\n":"=r"(var))
#define REG_GETC(var) asm volatile ("getc %0\n":"=r"(var))
#define NOP() asm volatile ("nop\n")
//----------------------------------------------------

extern int MY_RANK_MPI;
extern F64vec_ FORCE_SUM_LOC;

//=========================
/* Structure definitions */
//=========================
typedef double REAL;

typedef struct{
    REAL r[4]; //x,y,z,m
    REAL v[4]; //vx,vy,vz,id
} EpiMM;
typedef EpiMM EpjLM;

typedef struct{
    REAL r[4]; //x,y,z,m
    REAL v[4]; //vx,vy,vz,id
#ifdef PHI_R_TREE
    //REAL s[2];
#endif
} EpjMM;

typedef struct{
    REAL r[4];
#ifdef PHI_R_TREE
    REAL s[2];
#endif
} SpjMM;

typedef struct{
	REAL r[4];
} SpjLM;

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
  int*   n_walk_cpe;
  int*   n_disp_walk;
  int*   adr_n_walk;
  
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

//* FullParticle (used in Rotate, ex_ptcl)
typedef struct {
   F64vec_ pos;
   F64_ mass;
   F64vec_ vel;
   S64_ id;
} fpLM;

//* EssentialParticleI (used in setparticle for cylindrical coordinate)
typedef struct {
   F64vec_ pos;
   F64_ mass;
   F64vec_ vel;
   S64_ id;
} epiLM;

//* EssentialParticleJ (used in CalcMoment, CalcMomentLongGlobalTree, GenMortonKey)
typedef struct {
   F64vec_ pos;
   F64_ mass;
   F64vec_ vel;
   S64_ id;
#ifdef PHI_R_TREE
    //F64_ pos_phi;
    //F64_ pos_r;
#endif
} epjLM;

//* SuperParticleJ (used in CalcMoment, CalcMomentLongGlobalTree, GenMortonKey)
typedef struct {
   F64vec_ pos;
   F64_ mass;
#ifdef USE_QUADRUPOLE
   F32mat_ quad;
#endif
#ifdef PHI_R_TREE
    F64_ pos_phi;
    F64_ pos_r;
    //F64_ pad[2];
#endif
} spjLM;

//* TreeParticle (the order of vars. must be consistent with TreeParticle class in FDPS)

#ifdef USE_96BIT_KEY
typedef struct {
   KeyT_ key_;
   U32_ adr_ptcl_;
} tpLM;
#else
typedef struct {
   KeyT_ key_;
   U32_ adr_ptcl_;
   int pad_;
} tpLM;
#endif

/*
typedef struct {
   U64_ key_;
   U32_ adr_ptcl_;
   int pad_;
} tpLM;
*/


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
#ifndef REMOVE_VERTEX
   F64ort_ vertex_out_;
   F64ort_ vertex_in_;
#endif
#ifdef PHI_R_TREE
    F32_ pos_phi;
    F32_ pos_r;
#endif
} tcLM;

//* IPGroup (the order of vars. must be consistent with IPGroup class in FDPS)
typedef struct {
   S32_ n_ptcl_;
   S32_ adr_ptcl_;
   F64ort_ vertex_;
} ipgLM;

//* sample sort 
typedef struct KeyIndex {
	unsigned long key;
	unsigned long index;
} KeyIndex;

typedef struct sort_pars {
    int N;
    void* glb_input;//[N]
    void* glb_work;//[N]
    void* glb_result;//[N]

    unsigned long* glb_parts;//[NSPLIT-1];
    int* glb_task_off;//[NTHREADS+1];
    int* glb_status; //[NTHREADS];
    int* glb_count; //[NTHREADS][NSPLIT];
} sort_pars;

typedef struct{
	float x, y, z;
} floatvec_;

typedef unsigned short U16;

typedef struct {
    U32_ adr_tc_;
    U32_ adr_ptcl_;
    S32_ n_ptcl_;
    S32_ level_;
    F32_ mass;
    F32vec_ pos;
    //float mass;
    //floatvec_ pos;
    /*
#ifdef PHI_R_TREE
    F32_ pos_phi;
    F32_ pos_r;
#endif
    */
    /*
#ifdef __cplusplus
	template <typename Tp>
	void copyFrom(const Tp &src){
		n_ptcl_   = src.n_ptcl_;
		adr_tc_   = src.adr_tc_;
		adr_ptcl_ = src.adr_ptcl_;
		level_    = src.level_;

		mass   =  src.mom_.mass;
		pos.x  =  src.mom_.pos.x;
		pos.y  =  src.mom_.pos.y;
		pos.z  =  src.mom_.pos.z;
	}
	template <typename Tp>
	void copyTo(Tp &dst) const {
		dst.n_ptcl_   = n_ptcl_;
		dst.adr_tc_   = adr_tc_;
		dst.adr_ptcl_ = adr_ptcl_;
		dst.level_    = level_;

		dst.mom_.mass   =  mass;
		dst.mom_.pos.x  =  pos.x;
		dst.mom_.pos.y  =  pos.y;
		dst.mom_.pos.z  =  pos.z;
	}
#endif
    */
} etcLM;


//==========================
/* Prototype declarations */
//==========================
void Preproc1_CopyEPJLocToEPJGlb(void *);
void Preproc2_CopyEPJLocToEPJGlb(void *);
void CopyEPJLocToEPJGlb(void *);

//void CopyEPJLocToEPJ(void *);
void CopyIndirect(void *);
void CopyIndirectInverse(void *);
void CopyIndirectInverse2(void *);
//void CopyIndirect2(void *);
void CopyDirect(void *);
void CopyStride(void *);
void CopyStrideByHand(void *);
void CopyTCMomToSPJ(void *);
void GenMortonKey(void *);
void LinkCell(void *);
void CalcMoment(void *);
void CalcMomentLongGlobalTree(void *);
void MakeListUsingTree(void *);
void Rotate(void *);
void ForceKernelSunWay1st(void* pars);
void ForceKernelSunWay2nd(void* pars);

void balance_nwalk(void* pars);

void count_kernel(void* pars);
void split_kernel(void* pars);
void final_kernel(void* pars);

void CopyEPISortedToEPJSorted(void *);

void CalcMomentLean(void *);
void CalcMomentLongGlobalTreeLean(void *);

void CopyTCToETC(void *);

void CheckIsInDomain(void *);

void CopyTpLocAdrPtclToAdrPtcl(void *);
void CopyEpjOrgToEpiSortedEpjSorted(void *);

void GetMinIpgBox(void *);

void CopyFPToEPJ(void *);
void CopyEPIToFP(void *);

void SetParticleToSendSendBuffer(void *);

void RotateAndShiftZ(void *);

void SetTpAdrPtcl(void *);

void SetTpLocAdrPtclFromTpCpe(void *);

void SetAdrGlb1stCpe(void *);
void SetAdrGlb2ndCpe(void *);
void SetPtclSortedGlbCpe(void *);

void SetEpjTpFromBufferCpe(void *);
void SetSpjTpFromBufferCpe(void *);
