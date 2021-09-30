#ifndef WATERLATTICE_H
#define WATERLATTICE_H

#include "NoInputSystem.h"

template<class PA>
class WaterLattice : public NoInputSystem<PA> {
public:

  WaterLattice();
  WaterLattice(int lnum);
  WaterLattice(SpaceVector<int> lnum);

  int setLattice();
  void writePDB();
  
};  

#endif
