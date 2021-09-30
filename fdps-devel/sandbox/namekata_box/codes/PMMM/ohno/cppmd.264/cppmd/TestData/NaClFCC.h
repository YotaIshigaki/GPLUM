#ifndef NACLFCC_H
#define NACLFCC_H

#include "NoInputSystem.h"

template<class PA>
class NaClFCC : public NoInputSystem<PA> {
public:

  NaClFCC();
  NaClFCC(int lnum);
  NaClFCC(SpaceVector<int> lnum);

  int setLattice();
  void writePDB();
};

#endif

