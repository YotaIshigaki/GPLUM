#include <iostream>
#include "DummyLongRangeInteraction.h"

DummyChargeAssign::DummyChargeAssign(int unitid) 
    : ChargeAssign(unitid)
{
}

DummyChargeAssign::DummyChargeAssign(int unitid, const LongRangeParameter& _param) 
    : ChargeAssign(unitid,_param)
{
  //   std::cout << "construct DummyChargeAssign" << std::endl;
 
}

void DummyChargeAssign::assign(ParticleArray& particlearray, 
                               const std::vector<ParticleRange>& selfrange,
                               ParticleArray& ghost, 
                               const std::vector<ParticleRange>& ghostrange,
                               GridData& gridcharge) 
{
  //  std::cout << "  DummyChargeAssign::assign " << unit_identifier << std::endl;
}

void DummyChargeAssign::backinterpolate(ParticleArray& particlearray, 
                                        const std::vector<ParticleRange>& selfrange,
                                        ParticleArray& ghost, 
                                        const std::vector<ParticleRange>& ghostrange,
                                        GridData& gridpotential)
{
  std::cout << "  DummyChargeAssign::backinterpolate " << unit_identifier << std::endl;
}

DummyPoissonSolver::DummyPoissonSolver(int unitid) 
    : PoissonSolver(unitid) 
{
}

DummyPoissonSolver::DummyPoissonSolver(int unitid, const LongRangeParameter& _param) 
    : PoissonSolver(unitid,_param) 
{
}

void DummyPoissonSolver::solvePoisson(GridData& gridcharge,
                                      GridData& gridpotential, double& energy, double &virial)
{
  std::cout << "  DummyPoissonSolver " << unit_identifier << std::endl;
}
