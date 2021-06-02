#ifndef ARFCC_H
#define ARFCC_H

#include "NoInputSystem.h"

template<class PA>
class ArFCC : public NoInputSystem<PA> {
public:

  ArFCC();
  ArFCC(int lnum);
  ArFCC(SpaceVector<int> lnum);

  int setLattice();
  void writePDB();
};

#endif
