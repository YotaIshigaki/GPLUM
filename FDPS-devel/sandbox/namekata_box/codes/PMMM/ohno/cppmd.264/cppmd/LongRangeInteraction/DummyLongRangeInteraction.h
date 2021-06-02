#ifndef DUMMYLONGRANGEINTERACTION_H
#define DUMMYLONGRANGEINTERACTION_H
#include "LongRangeInteraction.h"

class DummyChargeAssign : public ChargeAssign {
public:

  DummyChargeAssign() : ChargeAssign(0) {}

  DummyChargeAssign(int unitid);
  DummyChargeAssign(int unitid, const LongRangeParameter& _param);
  void assign(ParticleArray& particlearray, 
              const std::vector<ParticleRange>& selfrange,
              ParticleArray& ghost, 
              const std::vector<ParticleRange>& ghostrange,
              GridData& gridcharge);
  void backinterpolate(ParticleArray& particlearray, 
                       const std::vector<ParticleRange>& selfrange,
                       ParticleArray& ghost, 
                       const std::vector<ParticleRange>& ghostrange,
                       GridData& gridpotential);
};

class DummyPoissonSolver : public PoissonSolver {
public:
  DummyPoissonSolver() : PoissonSolver(0) {}
  DummyPoissonSolver(int unitid);
  DummyPoissonSolver(int unitid, const LongRangeParameter& _param);

  void solvePoisson(GridData& gridcharge,
                    GridData& gridpotential, double& energy, double &virial);
};

typedef LongRangeInteraction<DummyChargeAssign, DummyPoissonSolver> DummyLongRangeInteraction;

#endif
