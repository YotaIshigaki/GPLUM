// -*- mode: C++; -*-
#ifndef REALSPACELONGRANGEINTERACTION_H
#define REALSPACELONGRANGEINTERACTION_H
#include "LongRangeInteraction.h"

class RealSpaceChargeAssign : public ChargeAssign {
public:
  RealSpaceChargeAssign(int unitid, const LongRangeParameter& _param);
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

class RealSpacePoissonSolver : public PoissonSolver {
public:
  RealSpacePoissonSolver(int unitid, const LongRangeParameter& _param);

  void solvePoisson(GridData& gridcharge,
                    GridData& gridpotential, double& energy);
};

typedef LongRangeInteraction<RealSpaceChargeAssign, RealSpacePoissonSolver>
    RealSpaceLongRangeInteraction;

#endif  // REALSPACELONGRANGEINTERACTION_H
