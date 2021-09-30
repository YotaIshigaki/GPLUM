#ifndef PAIRLISTINTERACTION_H
#define PAIRLISTINTERACTION_H

#include "OperationSelector.h"
#include "PairList.h"

#ifdef CHECK_ENERGY
extern double shortpair_energy;
extern double ljpair_energy;
#endif

struct waterexclude{
  int exid[2];

waterexclude() {exid[0]=-1;exid[1]=-1;};
};

void pairlistloopljfe1(const int *jindex,
                        const double (*lj)[4],
                        const int npair,
                        const double posi[3],
                        const double chargei,
                        const double (*ljmp)[4],
                        const PosChargeArray& particlej,
                        const double cutoff2,
                        double force[3],
                        double& energy);

void pairlistloopljcfe1(const int *jindex,
                        const int *lj,
                        const int npair,
                        const double posi[3],
                        const double chargei,
                        const double (*ljmp)[4],
                        const PosChargeArray& particlej,
                        const double cutoff2,
                        double force[3],
                        double& energy);
void pairlistloopljcfe1(const int *jindex,
                        const int *lj,
                        const int npair,
                        const double posi[3],
                        const double chargei,
                        const double (*ljmp)[4],
                        const PosChargeArray& particlej,
                        const double cutoff2,
                        double force[3],
                        double& energy,
			double& virial);

void pairlistloopljcfe1(const int *jindex,
                        const double (*lj)[4],
                        const int npair,
                        const double posi[3],
                        const double chargei,
                        const double (*ljmp)[4],
                        const PosChargeArray& particlej,
                        const double cutoff2,
                        double force[3],
                        double& energy);

void pairlistloopljcfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const PosChargeArray& poschargei,
                       const AtomtypeArray& atomia,
                       const PosChargeArray& poschargej,
                       const int npl,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy);

void pairlistloopljcfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const PosChargeArray& poschargei,
                       const AtomtypeArray& atomia,
                       const PosChargeArray& poschargej,
                       const int npl,
                       const int *selfnpair,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations);

template<typename PA, typename GPA>
void pairlistloopljfe(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljfe(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe_wolf_ljshift(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe_wolf_ljshift(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe_ewald_ljshift(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe_ewald_ljshift(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe1_ewald(const double alpha,
                              const int *jindex,
                              const int *lj,
                              const int npair,
                              const double posi[3],
                              const double chargei,
                              const double (*ljmp)[4],
                              const PosChargeArray& particlej,
                              const double cutoff2,
                              double force[3],
                              double& energy);
template<typename PA, typename GPA>
void pairlistloopljcfe_ewald(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe_ewald(const int (*jindex)[MAX_PAIR],
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

void pairlistloopljcfe1(const int *jindex,
                        const int *lj,
                        const int npair,
                        const double posi[3],
                        const double chargei,
                        const double (*ljmp)[4],
                        const ParticleArray& particlej,
                        const double cutoff2,
                        double force[3],
                        double& energy);

void pairlistloopljcfe1(const double (*pos)[3],
                        const double *charge,
                        const double (*lj)[4],
                        const int npair,
                        const double posi[3],
                        const double chargei,
                        const double cutoff2,
                        double force[3],
                        double& energy);

extern "C" void set_shift_coefficient_(double *cutoff);

extern "C" void pairlistloopljfe1_i_posc_f_(const int *jindex,
                                             const int *lj,
                                             const int *npair,
                                             const double posi[3],
                                             const double *chargei,
                                             const double (*ljmp)[4],
                                             const double *posc,
                                             const double *cutoff2,
                                             double force[3],
                                             double* energy);

extern "C" void pairlistloopljcfe1_i_posc_f_(const int *jindex,
                                             const int *lj,
                                             const int *npair,
                                             const double posi[3],
                                             const double *chargei,
                                             const double (*ljmp)[4],
                                             const double *posc,
                                             const double *cutoff2,
                                             double force[3],
                                             double* energy);

extern "C" void pairlistloopljewaldfe1_i_posc_f_(const double *alpha,
                                             const int *jindex,
                                             const int *lj,
                                             const int *npair,
                                             const double posi[3],
                                             const double *chargei,
                                             const double (*ljmp)[4],
                                             const double *posc,
                                             const double *cutoff2,
                                             double force[3],
                                             double* energy);

extern "C" void pairlistloopljcfe1_i_f_(const int *jindex,
                                        const double (*lj)[4],
                                        const int *npair,
                                        const double posi[3],
                                        const double *chargei,
                                        const double *posc,
                                        const double *cutoff2,
                                        double force[3],
                                        double* energy);

extern "C" void pairlistloopljcfe1_f_(const double (*pos)[3],
                                      const double *charge,
                                      const double (*lj)[4],
                                      const int *npair,
                                      const double posi[3],
                                      const double *chargei,
                                      const double *cutoff2,
                                      double force[3],
                                      double *energy);

void pairlistloopljcfe(const double (*pos)[MAX_PAIR][3],
                       const double (*charge)[MAX_PAIR],
                       const double (*lj)[MAX_PAIR][4],
                       const int *npair,
                       const double (*posi)[3],
                       const double *chargei,
                       const int npl,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy);


template<typename GPA>
void cellpairinteraction1(const double posi[3],
			  const double chargei,
			  const GPA& particlej,
			  const std::vector<TypeRange>& typerangej,
			  const std::vector<int>& shorttarget_index,
			  Force& force,
			  const int icell, const int i,
			  const bool self);
template<typename GPA>
void cellpairinteraction1(const double posi[3],
			  const double chargei,
			  const GPA& particlej,
			  const std::vector<TypeRange>& typerangej,
			  const std::vector<int>& shorttarget_index,
			  Force& force,
			  double& energy,
			  const int icell, const int i,
			  const bool self);

template<typename GPA>
void cellpairsinteraction(const CombinedParticleArray& particlei,
			  const std::vector<TypeRange>& typerangei,
			  const GPA& particlej,
			  const std::vector<TypeRange>& typerangej,
			  const std::vector<std::vector<int> >& shorttarget_index,
			  ForceArray& force,
			  double& energy,
			  const bool self,
			  const OperationSelector operations);


template<typename GPA>
void cellpairinteraction1(const ParticlePosCharge posci,
			  const std::vector<waterexclude> &waterexcludelist,
			  const GPA& particlej,
			  const std::vector<TypeRange>& typerangej,
			  const std::vector<int>& shorttarget_index,
			  Force& force,
			  double& energy,
			  const int icell, const int i,
			  const bool self);

template<typename GPA>
void cellpairsinteraction(const CombinedParticleArray& particlei,
			  const std::vector<TypeRange>& typerangei,
			  const std::vector<waterexclude> &waterexcludelist,
			  const GPA& particlej,
			  const std::vector<TypeRange>& typerangej,
			  const std::vector<std::vector<int> >& shorttarget_index,
			  ForceArray& force,
			  double& energy,
			  const bool self,
			  const OperationSelector operations);


#endif

