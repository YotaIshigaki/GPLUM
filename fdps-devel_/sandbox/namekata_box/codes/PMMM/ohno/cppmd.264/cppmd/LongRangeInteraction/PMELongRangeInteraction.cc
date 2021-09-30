#include "PMEInterface.h"
#include "PMELongRangeInteraction.h"

PMEChargeAssign::PMEChargeAssign(int unitid, const LongRangeParameter& _param) 
  : param(_param), pmeInterface(new PMEModule::PMEInterface(unitid, param)) {
  //  std::cout << "construct PMEChargeAssign pmeInterface " << pmeInterface << std::endl;
}

PMEChargeAssign::PMEChargeAssign(int unitid, const LongRangeParameter& _param,
                                 PMEModule::PMEInterface* pmei) 
  : param(_param), pmeInterface(pmei) {
  //  std::cout << "construct PMEChargeAssign pmeInterface " << pmeInterface << std::endl;
}

PMEChargeAssign::~PMEChargeAssign() {
  /*
  if (pmeInterface) {
    //    std::cout << "delete(pmeInterface);" << std::endl;
    delete(pmeInterface);
  }
  */
}

template<class PA, class GPA>
void PMEChargeAssign::assign(PA& particlearray,
                             const std::vector<ParticleRange>& selfrange,
                             GPA& ghost,
                             const std::vector<ParticleRange>& ghostrange,
                             GridData& gridcharge,
                             const std::vector<ParticleRange>& self_selfenergy_range,
                             const std::vector<ParticleRange>& ghost_selfenergy_range)
{
  /*
  std::cout << "PMEChargeAssign::assign use pmeInterface " << pmeInterface << std::endl;
  std::cout << "size of selfrange " << selfrange.size() << std::endl;
  std::cout << "size of ghostrange " << ghostrange.size() << std::endl;
  std::cout << " assign to gridcharge " << &gridcharge  << std::endl;
  */
  //  pmeInterface->chargeAssign(param, particlearray, selfrange, gridcharge);
  pmeInterface->chargeAssign(param, particlearray, selfrange, 
                             ghost, ghostrange, gridcharge,
                             self_selfenergy_range,ghost_selfenergy_range);
}
#ifdef OLDPARTICLE
template
void PMEChargeAssign::assign(ParticleArray& particlearray,
                             const std::vector<ParticleRange>& selfrange,
                             ParticleArray& ghost,
                             const std::vector<ParticleRange>& ghostrange,
                             GridData& gridcharge,
                             const std::vector<ParticleRange>& self_selfenergy_range,
                             const std::vector<ParticleRange>& ghost_selfenergy_range);
#else
template
void PMEChargeAssign::assign(CombinedParticleArray& particlearray,
                             const std::vector<ParticleRange>& selfrange,
                             GhostParticleArray& ghost,
                             const std::vector<ParticleRange>& ghostrange,
                             GridData& gridcharge,
                             const std::vector<ParticleRange>& self_selfenergy_range,
                             const std::vector<ParticleRange>& ghost_selfenergy_range);
#endif
template<class PA, class GPA>
void PMEChargeAssign::backinterpolate
                      (PA& particlearray,
                       const std::vector<ParticleRange>& selfrange,
                       GPA& ghost,
                       const std::vector<ParticleRange>& ghostrange,
                       GridData& gridpotential)
{
  //  std::cout << "PMEChargeAssign::backinterpolate" << std::endl;
  //  pmeInterface->backInterpolate(param, particlearray, selfrange, gridpotential);
  pmeInterface->backInterpolate(param, particlearray, selfrange, 
                                ghost, ghostrange, gridpotential);
}
#ifdef OLDPARTICLE
template
void PMEChargeAssign::backinterpolate
                      (ParticleArray& particlearray,
                       const std::vector<ParticleRange>& selfrange,
                       ParticleArray& ghost,
                       const std::vector<ParticleRange>& ghostrange,
                       GridData& gridpotential);
#else
template
void PMEChargeAssign::backinterpolate
                      (CombinedParticleArray& particlearray,
                       const std::vector<ParticleRange>& selfrange,
                       GhostParticleArray& ghost,
                       const std::vector<ParticleRange>& ghostrange,
                       GridData& gridpotential);
#endif

template<class PA, class GPA>
void PMEChargeAssign::addenergy(PA& particlearray,
                                const std::vector<ParticleRange>& selfrange,
                                GPA& ghost,
                                const std::vector<ParticleRange>& ghostrange,
                                double& energy)
{
  energy += pmeInterface->getSelfEnergy();
  energy += pmeInterface->getDipoleEnergy();
  //  printf("longenergy %e + %e + %e = %e\n",energy-pmeInterface->getSelfEnergy()-pmeInterface->getDipoleEnergy(),pmeInterface->getSelfEnergy(), pmeInterface->getDipoleEnergy(), energy);
}
#ifdef OLDPARTICLE
template
void PMEChargeAssign::addenergy(ParticleArray& particlearray,
                                const std::vector<ParticleRange>& selfrange,
                                ParticleArray& ghost,
                                const std::vector<ParticleRange>& ghostrange,
                                double& energy);
#else
template
void PMEChargeAssign::addenergy(CombinedParticleArray& particlearray,
                                const std::vector<ParticleRange>& selfrange,
                                GhostParticleArray& ghost,
                                const std::vector<ParticleRange>& ghostrange,
                                double& energy);
#endif

PMEModule::PMEInterface* PMEChargeAssign::get_pmeInterface()
{
  return pmeInterface;
}

PMEPoissonSolver::PMEPoissonSolver(int unitid,
                                   const LongRangeParameter& _param) 
  : param(_param), pmeInterface(new PMEModule::PMEInterface(unitid, param)) {
  //  std::cout << "construct PMEPoissonSolver pmeInterface " << pmeInterface << std::endl;
}

PMEPoissonSolver::PMEPoissonSolver(int unitid,
                                   const LongRangeParameter& _param,
                                   PMEModule::PMEInterface* pmei) 
  : param(_param), pmeInterface(pmei) {
  //  std::cout << "construct PMEPoissonSolver with *pmeInterface " << pmeInterface << std::endl;
}

PMEPoissonSolver::~PMEPoissonSolver() {
  /*
  if (pmeInterface) {
    delete(pmeInterface);
  }
  */
}

void PMEPoissonSolver::solvePoisson(GridData& gridcharge,
                                    GridData& gridpotential, double& energy,
				    double & virial)
{
  //  std::cout << "PMEChargeAssign::solvePoisson" << std::endl;
  energy = pmeInterface->solvePoisson(param, gridcharge, gridpotential,
				      virial);
}

PMEModule::PMEInterface* PMEPoissonSolver::get_pmeInterface()
{
  return pmeInterface;
}

template<>
void PMELongRangeInteraction::initialize()
{
  //  std::cout << "poissonsolver.get_pmeInterface()->initialize()" << std::endl;
  poissonsolver.get_pmeInterface()->initialize();
  //  std::cout << "poissonsolver.get_pmeInterface()->initialize() done" << std::endl;
}
