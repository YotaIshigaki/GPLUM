#ifndef CELLMETHODIMPL_H
#define CELLMETHODIMPL_H

#include "HalfShellCore.h"
#include "SmallBallCore.h"
#include "HalfShellCoreImpl.h"
#include "SmallBallCoreImpl.h"

namespace CellMethodModule {

void CellMethod::setCellIndexType(CellIndexType cellindextype)
{
  delete(core);
  core = 0;
  assigned = false;
  targets.clear();
  targetmap.clear();
  switch (cellindextype) {
  case HalfShell:
    core = new HalfShellCore(this);
    break;
  case FullShell:
    core = new HalfShellCore(this, true);
    break;
  case SmallBall:
    core = new SmallBallCore(this);
    break;
#if 1
  default:
    break;
#endif
  }
}

void CellMethod::target(const int celldiv,
                        const SpaceVector<int>& ci, 
                        const SpaceVector<double>& cellsize,
                        const SpaceVector<double>& margin,
                        const double cutoff,
                        std::vector<int>& send_target_cell, 
                        std::vector<int>& recv_target_cell, 
                        std::vector<int>& short_target_cell,
                        ShiftCellArray& shift_cell,
                        std::vector<Pair>& ghost_cell_pairs)
{
  CellMethodContext context(celldiv, cellsize, margin, cutoff);
  if (!assigned) {
    core->assignCalcCell(&context);
    assigned = true;
  }
  recv_target_cell.clear();
  short_target_cell.clear();
  shift_cell.clear();
  Indexes zero(0);
  std::vector<int> tidvec(core->getCalcVec().size(), -1);
  for (std::vector<Indexes>::size_type ivec=0;
       ivec<core->getCalcVec().size();++ivec) {
    const Indexes& calcv = core->getCalcVec()[ivec];
    Indexes calc;
    calc[Ix] = calcv.x;
    calc[Iy] = calcv.y;
    calc[Iz] = calcv.z;
    if (calc == zero) continue;
    Indexes recv = ci + calc;
    Indexes send = ci - calc;
    Indexes shift(0);
    for (int i=0;i<3;++i) {
      if (recv[i] < 0) {
        recv[i] += celldiv;
        shift[i] = -1;
      }
      if (recv[i] >= celldiv) {
        recv[i] -= celldiv;
        shift[i] = 1;
      }
      if (send[i] < 0) {
        send[i] += celldiv;
      }
      if (send[i] >= celldiv) {
        send[i] -= celldiv;
      }
    }
    int recv_cell = recv.x+recv.y*celldiv+recv.z*celldiv*celldiv;
    tidvec[ivec] = recv_cell;
    recv_target_cell.push_back(recv_cell);
    if(shift.x==1) shift_cell.plusx.push_back(recv_cell);
    if(shift.x==-1) shift_cell.minusx.push_back(recv_cell);
    if(shift.y==1) shift_cell.plusy.push_back(recv_cell);
    if(shift.y==-1) shift_cell.minusy.push_back(recv_cell);
    if(shift.z==1) shift_cell.plusz.push_back(recv_cell);
    if(shift.z==-1) shift_cell.minusz.push_back(recv_cell);
    int send_cell = send.x+send.y*celldiv+send.z*celldiv*celldiv;
    send_target_cell.push_back(send_cell);
  }

  int ci_cell = ci.x+ci.y*celldiv+ci.z*celldiv*celldiv;
  tmppairs.clear();
  for (std::vector<Pair>::const_iterator it = ghost_cell_pairs.begin();
       it != ghost_cell_pairs.end();++it) {
    if (it->first == ci_cell) {
      short_target_cell.push_back(it->second);
    }
    else if (it->second == ci_cell) {
      short_target_cell.push_back(it->first);
    }
    else {
      tmppairs.push_back(*it);
    }
  }
  ghost_cell_pairs.swap(tmppairs);
  for (std::vector<Pair>::const_iterator it = core->getCellPairs().begin();
       it != core->getCellPairs().end();++it) {
    int tid1 = tidvec[it->first];
    int tid2 = tidvec[it->second];
    if (tid1 < 0 && tid2 >= 0) {
      short_target_cell.push_back(tid2);
    }
    else if (tid2 < 0 && tid1 >= 0) {
      short_target_cell.push_back(tid1);
    }
    else if (tid1 >= 0 && tid2 >= 0) {
      TargetMap::iterator ittarget;
      if ((ittarget = targetmap.find(tid1)) != targetmap.end()) {
        targets[ittarget->second]->push_back(tid2);
      }
      else if ((ittarget = targetmap.find(tid2)) != targetmap.end()) {
        targets[ittarget->second]->push_back(tid1);
      }
      else {
        ghost_cell_pairs.push_back(Pair(tid1, tid2));
      }
    }
  }
  targetmap[ci_cell] = targets.size();
  targets.push_back(&short_target_cell);
}
}
#endif
