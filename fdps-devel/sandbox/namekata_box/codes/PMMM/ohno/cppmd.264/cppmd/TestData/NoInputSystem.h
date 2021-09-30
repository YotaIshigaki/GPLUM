#ifndef NOINPUTSYSTEM_H
#define NOINPUTSYSTEM_H

#include <string>
#include "SpaceVector.h"
#include "ParticleInfo.h"
#include "Common.h"

template<class PA>
class NoInputSystem {
public:

  SpaceVector<double> latticeSpacing;
  SpaceVector<int> latticeNum;
  SpaceVector<double> side;

  PA particle;
  std::vector<PotentialModel> potentialmodel;

  int natom;

  virtual ~NoInputSystem(){}

  virtual int setLattice()=0;
  virtual void writePDB()=0;

};

#endif
