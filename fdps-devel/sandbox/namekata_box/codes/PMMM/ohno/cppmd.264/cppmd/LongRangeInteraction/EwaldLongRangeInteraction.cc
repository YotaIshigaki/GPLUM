#include "EwaldInterface.h"
#include "EwaldLongRangeInteraction.h"

EwaldChargeAssign::EwaldChargeAssign(int unitid,
                                     const LongRangeParameter& _param)
  : param(_param),
    ewaldInterface(new EwaldModule::EwaldInterface(unitid, param)) {}

EwaldChargeAssign::EwaldChargeAssign(int unitid,
                                     const LongRangeParameter& _param,
                                     EwaldModule::EwaldInterface* ei)
  : param(_param),
    ewaldInterface(ei) {}

void EwaldChargeAssign::setSide(double side)
{
  ewaldInterface->setSide(side);
}

template<class PA, class GPA>
void EwaldChargeAssign::assign(PA& particlearray,
                               const std::vector<ParticleRange>& selfrange,
                               GPA& ghost,
                               const std::vector<ParticleRange>& ghostrange,
                               GridData& gridcharge,
                               const std::vector<ParticleRange>& self_selfenergy_range,
                               const std::vector<ParticleRange>& ghost_selfenergy_range)
{
  ewaldInterface->chargeAssign(param, particlearray, selfrange, 
                               ghost, ghostrange, gridcharge,
                               self_selfenergy_range,ghost_selfenergy_range);
}
#ifdef OLDPARTICLE
template
void EwaldChargeAssign::assign(ParticleArray& particlearray,
                               const std::vector<ParticleRange>& selfrange,
                               ParticleArray& ghost,
                               const std::vector<ParticleRange>& ghostrange,
                               GridData& gridcharge,
                               const std::vector<ParticleRange>& self_selfenergy_range,
                               const std::vector<ParticleRange>& ghost_selfenergy_range);
#else
template
void EwaldChargeAssign::assign(CombinedParticleArray& particlearray,
                               const std::vector<ParticleRange>& selfrange,
                               GhostParticleArray& ghost,
                               const std::vector<ParticleRange>& ghostrange,
                               GridData& gridcharge,
                               const std::vector<ParticleRange>& self_selfenergy_range,
                               const std::vector<ParticleRange>& ghost_selfenergy_range);
#endif

template<class PA, class GPA>
void EwaldChargeAssign::backinterpolate
                      (PA& particlearray,
                       const std::vector<ParticleRange>& selfrange,
                       GPA& ghost,
                       const std::vector<ParticleRange>& ghostrange,
                       GridData& gridpotential)
{
  ewaldInterface->backInterpolate(param, particlearray, selfrange,
                                  ghost, ghostrange, gridpotential);
}
#ifdef OLDPARTICLE
template
void EwaldChargeAssign::backinterpolate
                      (ParticleArray& particlearray,
                       const std::vector<ParticleRange>& selfrange,
                       ParticleArray& ghost,
                       const std::vector<ParticleRange>& ghostrange,
                       GridData& gridpotential);
#else
template
void EwaldChargeAssign::backinterpolate
                      (CombinedParticleArray& particlearray,
                       const std::vector<ParticleRange>& selfrange,
                       GhostParticleArray& ghost,
                       const std::vector<ParticleRange>& ghostrange,
                       GridData& gridpotential);
#endif

template<class PA, class GPA>
void EwaldChargeAssign::addenergy(PA& particlearray,
                                  const std::vector<ParticleRange>& selfrange,
                                  GPA& ghost,
                                  const std::vector<ParticleRange>& ghostrange,
                                  double& energy)
{
  energy += ewaldInterface->getSelfEnergy();
  energy += ewaldInterface->getDipoleEnergy();
}
#ifdef OLDPARTICLE
template
void EwaldChargeAssign::addenergy(ParticleArray& particlearray,
                                  const std::vector<ParticleRange>& selfrange,
                                  ParticleArray& ghost,
                                  const std::vector<ParticleRange>& ghostrange,
                                  double& energy);
#else
template
void EwaldChargeAssign::addenergy(CombinedParticleArray& particlearray,
                                  const std::vector<ParticleRange>& selfrange,
                                  GhostParticleArray& ghost,
                                  const std::vector<ParticleRange>& ghostrange,
                                  double& energy);
#endif

EwaldPoissonSolver::EwaldPoissonSolver(int unitid,
                                       const LongRangeParameter& _param)
  : param(_param),
    ewaldInterface(new EwaldModule::EwaldInterface(unitid, param)) {}

EwaldPoissonSolver::EwaldPoissonSolver(int unitid,
                                       const LongRangeParameter& _param,
                                       EwaldModule::EwaldInterface* ei)
  : param(_param),
    ewaldInterface(ei) {}

void EwaldPoissonSolver::solvePoisson(GridData& gridcharge,
                                      GridData& gridpotential, double& energy)
{
  energy = ewaldInterface->solvePoisson(param, gridcharge, gridpotential);
}

