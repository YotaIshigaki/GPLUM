#ifndef EWALDINTERFACEBASE_H
#define EWALDINTERFACEBASE_H

#include "UseMPI.h"
#include "LongRangeInteraction.h"

namespace EwaldModule {

struct EwaldBaseContext {
  EwaldBaseContext(const LongRangeParameter& _param) : param(_param) {}
  const LongRangeParameter& param;
};

template <typename T, typename TC=EwaldBaseContext>
class EwaldInterfaceBase {
public:
  typedef TC* Context;
  typedef TC BaseContext;

  EwaldInterfaceBase(const LongRangeParameter& param) 
    : method(), cd(), charge(), fc(), selfEnergyCalculated(false),
      isEwaldBaseNode(true) {
    //    std::cout << "construct EwaldInterfaceBase" << std::endl;
  }
  virtual ~EwaldInterfaceBase() {}

  inline void clear_cd_charge()
  {
    cd.clear();
    charge.clear();
  }

  /* (PME/Ewald)Interface -> (PME/Ewald)Method methods */
  template<class PA>
  void chargeAssign(const LongRangeParameter& param,
                    PA& particlearray,
                    const std::vector<ParticleRange>& selfrange,
                    GridData& gridcharge) {
    if (!isEwaldBaseNode) return;
    //        std::cout << "EwaldInterfaceBase::chargeAssign" << std::endl;
    TC context(param);
    clear_cd_charge();
    convert(particlearray, selfrange);
    calculateSelfEnergy();
    chargeAssignMain(&context);
    copyToGrid(gridcharge);
  }
  template<class PA, class GPA>
  void chargeAssign(const LongRangeParameter& param,
                    PA& particlearray,
                    const std::vector<ParticleRange>& selfrange,
                    GPA& ghost,
                    const std::vector<ParticleRange>& ghostrange,
                    GridData& gridcharge,
                    const std::vector<ParticleRange>& self_selfenergy_range,
                    const std::vector<ParticleRange>& ghost_selfenergy_range) {
    if (!isEwaldBaseNode) return;
    TC context(param);
    if (!selfEnergyCalculated) {
      clear_cd_charge();
      convert(particlearray, self_selfenergy_range);
      convert(ghost, ghost_selfenergy_range);
      calculateSelfEnergy();
    }
    clear_cd_charge();
    convert(particlearray, selfrange);
    convert(ghost, ghostrange);
    chargeAssignMain(&context);
  }
  double solvePoisson(const LongRangeParameter& param,
                      const GridData& gridcharge, GridData& gridpotential,
		      double &virial) {
    if (!isEwaldBaseNode) return 0.0;
    TC context(param);
    double potEwald = method->calculateEwaldEnergy(&context,virial);
    return potEwald;
  }
  template<class PA>
  void backInterpolate(const LongRangeParameter& param,
                       PA& particlearray,
                       const std::vector<ParticleRange>& selfrange,
                       const GridData& gridpotential) {
    if (!isEwaldBaseNode) return;
    TC context(param);
    clear_cd_charge();
    convert(particlearray, selfrange);
    copyFromGrid(gridpotential);
    fc.assign(cd.size(), Force());
    backInterpolateMain(&context);
    convertForce(particlearray, selfrange);
  }
  template<class PA, class GPA>
  void backInterpolate(const LongRangeParameter& param,
                       PA& particlearray,
                       const std::vector<ParticleRange>& selfrange,
                       GPA& ghost,
                       const std::vector<ParticleRange>& ghostrange,
                       const GridData& gridpotential) {
    if (!isEwaldBaseNode) return;
    TC context(param);
    fc.assign(cd.size(), Force());
    backInterpolateMain(&context);
    convertForce(particlearray, selfrange, ghost, ghostrange);
  }
  double getSelfEnergy() { return method->getSelfEnergy(); }
  double getDipoleEnergy() { return method->getDipoleEnergy(); }
  
  /* (PME/Ewald)Method -> (PME/Ewald)Interface methods */
  const SpaceVector<double>& getSide(Context pContext) {
    return pContext->param.boxSize;
  }
  double getVolume(Context pContext) {
    const SpaceVector<double>& box = pContext->param.boxSize;
    return box[0]*box[1]*box[2];
  }
  const SpaceVector<double> getWaveVector(Context pContext,
                                          const SpaceVector<int>& kv) {
    const SpaceVector<double>& box = pContext->param.boxSize;
    return SpaceVector<double>(static_cast<double>(kv[0])/box[0],
                               static_cast<double>(kv[1])/box[1],
                               static_cast<double>(kv[2])/box[2]);
  }
  /* end of (PME/Ewald)Method -> (PME/Ewald)Interface methods */
private:
  virtual void chargeAssignMain(Context pContext)=0;
  virtual void backInterpolateMain(Context pContext)=0;

  template<class PA>
  void convert(PA& particlearray,
               const std::vector<ParticleRange>& selfrange) {
    for (std::vector<ParticleRange>::size_type it=0;it < selfrange.size();
         ++it) {
      for (int i = selfrange[it].begin; i < selfrange[it].end; ++i) {
        cd.push_back(getpos(particlearray,i));
        charge.push_back(getcharge(particlearray,i));
      }
    }
  }
  template<class PA>
  void convertForce(PA& particlearray,
                    const std::vector<ParticleRange>& selfrange)
  {
    std::vector<Force>::const_iterator itfc = fc.begin();
    for (std::vector<ParticleRange>::size_type it=0;it < selfrange.size();
         ++it) {
      for (int i = selfrange[it].begin; i < selfrange[it].end; ++i,++itfc) {
           getforce(particlearray,i) += *itfc;
      }
    }
  }
  template<class PA, class GPA>
  void convertForce(PA& particlearray,
                    const std::vector<ParticleRange>& selfrange,
                    GPA& ghost,
                    const std::vector<ParticleRange>& ghostrange)
  {
    std::vector<Force>::const_iterator itfc = fc.begin();
    for (std::vector<ParticleRange>::size_type it=0;it < selfrange.size();
         ++it) {
      for (int i = selfrange[it].begin; i < selfrange[it].end; ++i,++itfc) {
        getforce(particlearray,i) += *itfc;
      }
    }
    for (std::vector<ParticleRange>::size_type it=0;it < ghostrange.size();
         ++it) {
      for (int i = ghostrange[it].begin; i < ghostrange[it].end; ++i,++itfc) {
           getforce(ghost,i) += *itfc;
      }
    }
  }
  virtual void copyFromGrid(const GridData& data)=0;
  virtual void copyToGrid(GridData& data)=0;

protected:
  void calculateSelfEnergy() {
    if (!selfEnergyCalculated) {
      method->calculateSelfEnergy(charge);
      selfEnergyCalculated = true;
    }
  }

  T* method;
  std::vector<Position> cd;
  std::vector<double> charge;
  std::vector<Force> fc;
  bool selfEnergyCalculated;
  bool isEwaldBaseNode;
};
}
#endif
