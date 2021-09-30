#pragma once

//#define N_WALK_LIMIT 500
#define N_WALK_LIMIT 1
#define NI_LIMIT     100000
#define NJ_LIMIT     1000000
#define N_CPE_LIMIT  64

//typedef float REAL;
typedef double REAL;

typedef struct{
	REAL x;
	REAL y;
	REAL z;
        int  id;
} EpiMM;

typedef struct{
	REAL x;
	REAL y;
	REAL z;
        int adr;
	REAL mass;
        int id;
} EpjMM;

typedef struct{
	REAL x;
	REAL y;
	REAL z;
	REAL mass;
} SpjMM;

typedef struct{
	REAL ax;
	REAL ay;
	REAL az;
	REAL pot;
        REAL r_ngb_sq;
        int  n_coll;
        int  adr_ngb;
} ForceMM;

extern int     disp_mm[N_WALK_LIMIT+2][3];
extern EpiMM   epi_mm[NI_LIMIT];
extern EpjMM   epj_mm[NJ_LIMIT];
extern SpjMM   spj_mm[NJ_LIMIT];
extern ForceMM force_mm[NI_LIMIT];
extern double  r_coll_sq_mm;

void DispatchKernelSunWay(const int n_walk, const double eps2);
void RetrieveKernelSunWay(const int n_epi);

