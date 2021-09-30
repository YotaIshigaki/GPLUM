#ifndef PMELONGRANGEINTERACTION_H
#define PMELONGRANGEINTERACTION_H

#include "LongRangeInteraction.h"
#include "PMEInterfaceFwd.h"

class PMEChargeAssign : public ChargeAssign {
public:
  PMEChargeAssign(int unitid, const LongRangeParameter& _param);
  PMEChargeAssign(int unitid, const LongRangeParameter& _param,
                  PMEModule::PMEInterface* pmei);
  ~PMEChargeAssign();

  template<class PA, class GPA>
  void assign(PA& particlearray,
              const std::vector<ParticleRange>& selfrange,
              GPA& ghost,
              const std::vector<ParticleRange>& ghostrange,
              GridData& gridcharge,
              const std::vector<ParticleRange>& self_selfenergy_range,
              const std::vector<ParticleRange>& ghost_selfenergy_range);
  template<class PA, class GPA>
  void backinterpolate(PA& particlearray,
                       const std::vector<ParticleRange>& selfrange,
                       GPA& ghost,
                       const std::vector<ParticleRange>& ghostrange,
                       GridData& gridpotential);
  template<class PA, class GPA>
  void addenergy(PA& particlearray,
                 const std::vector<ParticleRange>& selfrange,
                 GPA& ghost,
                 const std::vector<ParticleRange>& ghostrange,
                 double& energy);
  PMEModule::PMEInterface* get_pmeInterface();
private:
  const LongRangeParameter& param;
  PMEModule::PMEInterface* pmeInterface;
};

class PMEPoissonSolver : public PoissonSolver {
public:
  PMEPoissonSolver(int unitid, const LongRangeParameter& _param);
  PMEPoissonSolver(int unitid, const LongRangeParameter& _param,
                   PMEModule::PMEInterface* pmei);
  ~PMEPoissonSolver();
  
  void solvePoisson(GridData& gridcharge,
                    GridData& gridpotential, double& energy,
		    double &virial);
  PMEModule::PMEInterface* get_pmeInterface();
private:
  const LongRangeParameter& param;
  PMEModule::PMEInterface* pmeInterface;
  
};

typedef LongRangeInteraction<PMEChargeAssign, PMEPoissonSolver, PMEModule::PMEInterface>
    PMELongRangeInteraction;

#endif // PMELONGRANGEINTERACTION_H
