#ifndef SMALLBALLCORE_H
#define SMALLBALLCORE_H

#include "CellMethodCore.h"

namespace CellMethodModule {

class SmallBallCore : public CellMethodCore {
public:
  SmallBallCore(CellMethod* method)
    : CellMethodCore(method) {}

  void assignCalcCell(Context pContext);
private:
  class CheckArray {
  public:
    CheckArray(const CellRange& range)
      : base(range.min),
        strides((range.max[1]-range.min[1]+1)*(range.max[2]-range.min[2]+1),
                range.max[2]-range.min[2]+1, 1),
        data((range.max[0]-range.min[0]+1)*strides[0], false) {}
    bool operator()(int x, int y, int z) const {
      return data[getIndex(x, y, z)];
    }
    void set(int x, int y, int z) {
      data[getIndex(x, y, z)] = true;
    }
    void reset(int x, int y, int z) {
      data[getIndex(x, y, z)] = false;
    }
  private:
    int getIndex(int x, int y, int z) const {
      return (x-base.x)*strides.x+(y-base.y)*strides.y+(z-base.z)*strides.z;
    }
    Indexes base;
    Indexes strides;
    std::vector<bool> data;
  };
  struct SortIndex {
    int layer;
    int directionFactor;
    double distance;
    Indexes ipos;
    int id;
    int adoptedID;
    bool isAdopted;
  };
  struct lessSortIndex {
    bool operator()(SortIndex q1, SortIndex q2) const
    {
      return (q1.layer < q2.layer) ||
            ((q1.layer == q2.layer) &&
             (q1.directionFactor < q2.directionFactor)) ||
            ((q1.layer == q2.layer) &&
             (q1.directionFactor == q2.directionFactor) &&
             (q1.distance < q2.distance));
    }
  };
};
}
#endif
