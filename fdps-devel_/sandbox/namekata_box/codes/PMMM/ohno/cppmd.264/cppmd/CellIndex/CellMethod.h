#ifndef CELLMETHOD_H
#define CELLMETHOD_H

#include "CellMethodCore.h"

namespace CellMethodModule {

#ifdef CELL_ORDER_CPPMD0
const int Ix = 2;
const int Iy = 1;
const int Iz = 0;
#else
const int Ix = 0;
const int Iy = 1;
const int Iz = 2;
#endif

struct CellMethodContext {
  typedef SpaceVector<int> Indexes;

  CellMethodContext(const int _celldiv,
                    const SpaceVector<double>& _cellsize,
                    const SpaceVector<double>& _margin,
                    const double _cutoff)
    : celldiv(_celldiv), cellsize(_cellsize),
      margin(_margin), cutoff(_cutoff),
      cutoffCellNum(int(ceil((cutoff+2.0*margin[Ix])/cellsize[Ix])),
                    int(ceil((cutoff+2.0*margin[Iy])/cellsize[Iy])),
                    int(ceil((cutoff+2.0*margin[Iz])/cellsize[Iz]))),
      cellShape(_celldiv) {}

  const int celldiv;
  const SpaceVector<double>& cellsize;
  const SpaceVector<double>& margin;
  const double cutoff;
  const Indexes cutoffCellNum;
  const Indexes cellShape;
};

class CellMethod {
public:
  typedef CellMethodContext* Context;
  typedef CellMethodCore* CorePtr;
  typedef SpaceVector<int> Indexes;
  typedef CellMethodCore::Pair Pair;
  typedef std::map<int, int> TargetMap;

  CellMethod() : core(), assigned(false), targets(), targetmap(), tmppairs() {}
  ~CellMethod() {
    delete(core);
  }

  void setCellIndexType(CellIndexType cellindextype);

  /* CellMethod -> CellMethodCore methods */
  void target(const int celldiv,
              const SpaceVector<int>& ci, 
              const SpaceVector<double>& cellsize,
              const SpaceVector<double>& margin,
              const double cutoff,
              std::vector<int>& send_target_cell, 
              std::vector<int>& recv_target_cell, 
              std::vector<int>& short_target_cell,
              ShiftCellArray& shift_cell,
              std::vector<Pair>& ghost_cell_pairs);
  /* end of CellMethod -> CellMethodCore methods */

  /* CellMethodCore -> CellMethod methods */
  Indexes getCutoffCellNum(Context pContext) { return pContext->cutoffCellNum; }
  bool inCutoffSphere(Context pContext, const Indexes& ipos) {
    Indexes vpos;
    vpos[Ix] = ipos.x;
    vpos[Iy] = ipos.y;
    vpos[Iz] = ipos.z;
    return (min_cell_distance(Indexes(0), vpos,
                              pContext->cellsize,
                              pContext->margin) < pContext->cutoff);
  }
  double calculateDistance(Context pContext, const Indexes& ipos) {
    Position v(pContext->cellsize[Ix]*ipos.x,
               pContext->cellsize[Iy]*ipos.y,
               pContext->cellsize[Iz]*ipos.z);
    return v.norm();
  }
  Indexes getCellShape(Context pContext) { return pContext->cellShape; }
  /* end of CellMethodCore -> CellMethod methods */
private:
  CorePtr core;
  bool assigned;
  std::vector<std::vector<int>*> targets;
  TargetMap targetmap;
  std::vector<Pair> tmppairs;
};
}
#endif
