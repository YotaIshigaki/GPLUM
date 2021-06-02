#ifndef EWALDLONGRANGEINTERACTION_H
#define EWALDLONGRANGEINTERACTION_H
#include "LongRangeInteraction.h"
#include "EwaldInterfaceFwd.h"

class EwaldChargeAssign : public ChargeAssign {
public:
  EwaldChargeAssign(int unitid, const LongRangeParameter& _param);
  EwaldChargeAssign(int unitid, const LongRangeParameter& _param,
                    EwaldModule::EwaldInterface* ei);
  ~EwaldChargeAssign() {
    /*
    if (ewaldInterface) {
      delete(ewaldInterface);
    }
    */
  }

  void setSide(double side);

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
private:
  const LongRangeParameter& param;
  EwaldModule::EwaldInterface* ewaldInterface;
};

class EwaldPoissonSolver : public PoissonSolver {
public:
  EwaldPoissonSolver(int unitid, const LongRangeParameter& _param);
  EwaldPoissonSolver(int unitid, const LongRangeParameter& _param,
                     EwaldModule::EwaldInterface* ei);
  ~EwaldPoissonSolver() {
    /*
    if (ewaldInterface) {
      delete(ewaldInterface);
    }
    */
  }

  void solvePoisson(GridData& gridcharge,
                    GridData& gridpotential, double& energy);
private:
  const LongRangeParameter& param;
  EwaldModule::EwaldInterface* ewaldInterface;
};

typedef LongRangeInteraction<EwaldChargeAssign, EwaldPoissonSolver,EwaldModule::EwaldInterface>
    EwaldLongRangeInteraction;
#endif // EWALDLONGRANGEINTERACTION_H
