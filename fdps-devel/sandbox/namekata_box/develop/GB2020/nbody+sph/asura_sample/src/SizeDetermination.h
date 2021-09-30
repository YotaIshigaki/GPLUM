#pragma once 

enum {
    CS_TypeHydro,
    CS_TypeHII,
    CS_TypeSN,
#ifdef USE_RADIATION_PRESSURE //{
    CS_TypeRP,
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
    CS_TypeSW,
#endif // USE_STELLAR_WIND //}
    CS_TypeNumber,
};


#include "StellarFeedback.h"

void CalcSize(void);
int ReturnCalcSizeElementNumber(const int Type, const bool Global);
int ReturnCalcSizeOffset(const int Type);
//void CalcSizeGetHydroInfo_i(const int Index, double *PotentialMin, double VCOM[]);
//void CalcSizeGetHydroInfo_i(const int Index, double *PotentialMin, double VCOM[], bool *NoLocalSink);
void CalcSizeGetHydroInfo_i(const int Index, double *PotentialMin, double *DensityMax, double VCOM[], bool *NoLocalSink);
void CalcSizeSetSNInfo(struct StructActiveStellarFeedbackParticle ActiveStellarFeedbackParticle[]);
#ifdef USE_RADIATION_PRESSURE //{
void CalcSizeSetRPInfo(struct StructActiveStellarFeedbackParticle ActiveStellarFeedbackParticle[]);
#endif // USE_RADIATION_PRESSURE //}
#ifdef USE_STELLAR_WIND //{
void CalcSizeSetSWInfo(struct StructActiveStellarFeedbackParticle ActiveStellarFeedbackParticle[]);
#endif // USE_STELLAR_WIND //}

#ifdef USE_NEIGHBOR_LIST //{
struct StructGetLocalNeighborList{
    int Nlist;
    int *Neighbors;
    double Kernel;
};

struct StructGetLocalNeighborList GetLocalNeighrborList(const int Index);
#endif // USE_NEIGHBOR_LIST //}
