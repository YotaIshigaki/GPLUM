#ifndef PAIRLISTZERODIPOLE_H
#define PAIRLISTZERODIPOLE_H

#include "PairListInteraction.h"

#ifdef CHECK_ENERGY
extern double zerodipole_energy;
#endif

template<typename PA>
double calc_self_energy_zerodipole(const int *iindex,
				   const PA& particlei,
				   const int npl,
				   const double cutoff2);

void pairlistloopljcf1_zerodipole(const double alpha,
				  const int *jindex,
				  const int *lj,
				  const int npair,
				  const double posi[3],
				  const double chargei,
				  const double (*ljmp)[4],
				  const PosChargeArray& particlej,
				  const double cutoff2,
				  double force[3]);

void pairlistloopljcfe1_zerodipole(const double alpha,
				   const int *jindex,
				   const int *lj,
				   const int npair,
				   const double posi[3],
				   const double chargei,
				   const double (*ljmp)[4],
				   const PosChargeArray& particlej,
				   const double cutoff2,
				   double force[3],
				   double &energy);

void pairlistloopljcfe1_zerodipole(const double alpha,
				   const int *jindex,
				   const int *lj,
				   const int npair,
				   const double posi[3],
				   const double chargei,
				   const double (*ljmp)[4],
				   const PosChargeArray& particlej,
				   const double cutoff2,
				   double force[3],
				   double &energy,
				   double &virial);

template<typename PA, typename GPA>
void pairlistloopljcfe_zerodipole(const int (*jindex)[MAX_PAIR],
				  const int (*lj)[MAX_PAIR],
				  const int *npair,
				  const int *iindex,
				  const PA& particlei,
				  const GPA& particlej,
				  const int npl,
				  const double cutoff2,
				  double (*force)[3],
				  double& energy,
				  double& virial,
				  const OperationSelector operations);

template<typename PA, typename GPA>
void pairlistloopljcfe_zerodipole(const int (*jindex)[MAX_PAIR],
				  const int (*lj)[MAX_PAIR],
				  const int *npair,
				  const int *iindex,
				  const PA& particlei,
				  const GPA& particlej,
				  const int npl,
				  const int *selfnpair,
				  const double cutoff2,
				  double (*force)[3],
				  double& energy,
				  double& virial,
				  const OperationSelector operations);

void pairlistloopljcf1_zerodipole_ljshift(const double alpha,
					   const int *jindex,
					   const int *lj,
					   const int npair,
					   const double posi[3],
					   const double chargei,
					   const double (*ljmp)[4],
					   const PosChargeArray& particlej,
					   const double cutoff2,
					  double force[3]);

void pairlistloopljcfe1_zerodipole_ljshift(const double alpha,
					   const int *jindex,
					   const int *lj,
					   const int npair,
					   const double posi[3],
					   const double chargei,
					   const double (*ljmp)[4],
					   const PosChargeArray& particlej,
					   const double cutoff2,
					   double force[3],
					   double &energy);

void pairlistloopljcfe1_zerodipole_ljshift(const double alpha,
					   const int *jindex,
					   const int *lj,
					   const int npair,
					   const double posi[3],
					   const double chargei,
					   const double (*ljmp)[4],
					   const PosChargeArray& particlej,
					   const double cutoff2,
					   double force[3],
					   double &energy,
					   double &virial);

template<typename PA, typename GPA>
void pairlistloopljcfe_zerodipole_ljshift(const int (*jindex)[MAX_PAIR],
					  const int (*lj)[MAX_PAIR],
					  const int *npair,
					  const int *iindex,
					  const PA& particlei,
					  const GPA& particlej,
					  const int npl,
					  const double cutoff2,
					  double (*force)[3],
					  double& energy,
					  double& virial,
					  const OperationSelector operations);

template<typename PA, typename GPA>
void pairlistloopljcfe_zerodipole_ljshift(const int (*jindex)[MAX_PAIR],
					  const int (*lj)[MAX_PAIR],
					  const int *npair,
					  const int *iindex,
					  const PA& particlei,
					  const GPA& particlej,
					  const int npl,
					  const int *selfnpair,
					  const double cutoff2,
					  double (*force)[3],
					  double& energy,
					  double& virial,
					  const OperationSelector operations);

template<typename PA, typename GPA>
  void pairlistloopljcfe_zerodipole(const int (*jindex)[MAX_PAIR],
				    const int (*lj)[MAX_PAIR],
				    const int *npair,
				    const int *iindex,
				    const PA& particlei,
				    const GPA& particlej,
				    const int npl,
				    const double cutoff2,
				    double (*force)[3],
				    double& energy,
				    double& virial,
				    const OperationSelector operations);

template<typename PA, typename GPA>
  void pairlistloopljcfe_zerodipole(const int (*jindex)[MAX_PAIR],
				    const int (*lj)[MAX_PAIR],
				    const int *npair,
				    const int *iindex,
				    const PA& particlei,
				    const GPA& particlej,
				    const int npl,
				    const int *selfnpair,
				    const double cutoff2,
				    double (*force)[3],
				    double& energy,
				    double& virial,
				    const OperationSelector operations);

template<typename PA, typename GPA>
void pairlistloopljcfe_zerodipole0_ljshift(const int (*jindex)[MAX_PAIR],
					  const int (*lj)[MAX_PAIR],
					  const int *npair,
					  const int *iindex,
					  const PA& particlei,
					  const GPA& particlej,
					  const int npl,
					  const double cutoff2,
					  double (*force)[3],
					  double& energy,
					  double& virial,
					  const OperationSelector operations);

template<typename PA, typename GPA>
void pairlistloopljcfe_zerodipole0_ljshift(const int (*jindex)[MAX_PAIR],
					  const int (*lj)[MAX_PAIR],
					  const int *npair,
					  const int *iindex,
					  const PA& particlei,
					  const GPA& particlej,
					  const int npl,
					  const int *selfnpair,
					  const double cutoff2,
					  double (*force)[3],
					  double& energy,
					  double& virial,
					  const OperationSelector operations);


template<typename PA, typename GPA>
  void pairlistloopljcfe_zerodipole0(const int (*jindex)[MAX_PAIR],
				     const int (*lj)[MAX_PAIR],
				     const int *npair,
				     const int *iindex,
				     const PA& particlei,
				     const GPA& particlej,
				     const int npl,
				     const double cutoff2,
				     double (*force)[3],
				     double& energy,
				     double& virial,
				     const OperationSelector operations);

template<typename PA, typename GPA>
  void pairlistloopljcfe_zerodipole0(const int (*jindex)[MAX_PAIR],
				     const int (*lj)[MAX_PAIR],
				     const int *npair,
				     const int *iindex,
				     const PA& particlei,
				     const GPA& particlej,
				     const int npl,
				     const int *selfnpair,
				     const double cutoff2,
				     double (*force)[3],
				     double& energy,
				     double& virial,
				     const OperationSelector operations);

#endif
