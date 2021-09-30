#ifndef MULTIGRIDLONGRANGEINTERACTION_H
#define MULTIGRIDLONGRANGEINTERACTION_H

#include "LongRangeInteraction.h"
#include "MultigridFwd.h"
#include "Array3D.h"

class MultigridChargeAssign : public ChargeAssign {
public:
  typedef PMEModule::Array3D<double> R3D;

  MultigridChargeAssign(int unitid, const LongRangeParameter& _param);
  ~MultigridChargeAssign();

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
  void addenergy(ParticleArray& particlearray,
                 const std::vector<ParticleRange>& selfrange,
                 ParticleArray& ghost,
                 const std::vector<ParticleRange>& ghostrange,
                 double& energy);
private:
  const LongRangeParameter& param;
  int nodeid;
  MultigridModule::ChargeAssignmentMethod* caMethod;
  GeometryXYZ node_geometry;
  R3D::Array *p_rho;
  R3D::Array *p_phi;
  std::vector<Position> cd;
  std::vector<double> charge;
  std::vector<Force> fc;
  bool selfEnergyCalculated;
  double selfEnergy;
};

class MultigridPoissonSolver : public PoissonSolver {
public:
  typedef PMEModule::Array3D<double> R3D;

  MultigridPoissonSolver(int unitid, const LongRangeParameter& _param);
  ~MultigridPoissonSolver();
  
  void solvePoisson(GridData& gridcharge,
                    GridData& gridpotential, double& energy);
private:
  const LongRangeParameter& param;
  int nodeid;
#ifdef ORIGINAL_MGPOISSON
  MGPoisson* poisson;
#else
  MultigridModule::MGPoisson* poisson;
#endif
  GeometryXYZ node_geometry;
  SpaceVector<int> pos;
  SpaceVector<int> localNumGrid;
  SpaceVector<double> localSide;
  R3D::Array *p_rho;
  R3D::Array *p_phi;
};

typedef LongRangeInteraction<MultigridChargeAssign, MultigridPoissonSolver>
    MultigridLongRangeInteraction;

#endif // MULTIGRIDLONGRANGEINTERACTION_H
