#ifndef SMALLBALLCOREIMPL_H
#define SMALLBALLCOREIMPL_H

#include "CellMethodCore.h"

namespace CellMethodModule {

void SmallBallCore::assignCalcCell(Context pContext)
{
  CellRange checkRange;
  Indexes cutoffCellNum = method->getCutoffCellNum(pContext);
  checkRange.min = -cutoffCellNum;
  checkRange.max = cutoffCellNum;
  CheckArray checkArray(checkRange);
  std::vector<SortIndex> sortIndexes;
  for (int ix = checkRange.min.x;ix <= checkRange.max.x;++ix) {
    for (int iy = checkRange.min.y;iy <= checkRange.max.y;++iy) {
      for (int iz = checkRange.min.z;iz <= checkRange.max.z;++iz) {
        Indexes ipos(ix, iy, iz);
        if (method->inCutoffSphere(pContext, ipos)) {
          SortIndex si;
          si.layer = std::max(abs(ix), std::max(abs(iy), abs(iz)));
          si.directionFactor = 
              ((ix == si.layer)?1:0) +
              ((iy == si.layer)?2:0) +
              ((iz == si.layer)?4:0) +
              ((ix == -si.layer)?8:0) +
              ((iy == -si.layer)?16:0) +
              ((iz == -si.layer)?32:0);
          si.distance = method->calculateDistance(pContext, ipos);
          si.ipos = ipos;
          si.isAdopted = false;
          sortIndexes.push_back(si);
          checkArray.set(ix, iy, iz);
        }
      }
    }
  }

  std::sort(sortIndexes.begin(), sortIndexes.end(), lessSortIndex());
  int id = 0;
  for (std::vector<SortIndex>::iterator it = sortIndexes.begin();
       it != sortIndexes.end();++it,++id) {
    it->id = id;
  }
  if (sortIndexes.size() == 0) {
    SortIndex si;
    si.layer = 0;
    si.ipos[0] = 0;
    si.isAdopted = false;
    sortIndexes.push_back(si);
  }
  assert(sortIndexes.size() > 0);
  assert(sortIndexes[0].ipos == Indexes(0));
  sortIndexes[0].isAdopted = true;
  checkArray.reset(0, 0, 0);
  std::vector<Pair> pairIDs(1);
  pairIDs[0].first = 0;
  pairIDs[0].second = 0;

  for (std::vector<SortIndex>::iterator it = sortIndexes.begin();
       it != sortIndexes.end();++it) {
    if (!checkArray(it->ipos[0], it->ipos[1], it->ipos[2])) {
      continue;
    }
    bool found = false;
    for (std::vector<SortIndex>::iterator jt = sortIndexes.begin();
         jt != it+1;++jt) {
      if (jt->isAdopted) {
        continue;
      }
      for (std::vector<SortIndex>::iterator kt = sortIndexes.begin();
           kt != jt;++kt) {
        if (kt->isAdopted) {
          Indexes dpos = jt->ipos - kt->ipos;
          if (dpos == it->ipos || -dpos == it->ipos) {
            jt->isAdopted = true;
            found = true;
            break;
          }
        }
      }
      if (found) {
        for (std::vector<SortIndex>::iterator kt = sortIndexes.begin();
             kt != jt;++kt) {
          if (!kt->isAdopted) {
            continue;
          }
          Indexes dpos = jt->ipos - kt->ipos;
          bool overFlag = false;
          for (int j=0;j<3;++j) {
            if (dpos[j] < -cutoffCellNum[j] || dpos[j] > cutoffCellNum[j]) {
              overFlag = true;
            }
            if (2*dpos[j] == -method->getCellShape(pContext)[j]) {
              dpos[j] = -dpos[j];
            }
          }
          if (!overFlag && 
              (checkArray(dpos[0], dpos[1], dpos[2]) ||
               checkArray(-dpos[0], -dpos[1], -dpos[2]))) {
            checkArray.reset(dpos[0], dpos[1], dpos[2]);
            checkArray.reset(-dpos[0], -dpos[1], -dpos[2]);
            Pair pairID;
            if (jt->id <= kt->id) {
              pairID.first = jt->id;
              pairID.second = kt->id;
            }
            else {
              pairID.first = kt->id;
              pairID.second = jt->id;
            }
            pairIDs.push_back(pairID);
          }
        }
        break;
      }
    }
  }

  calcVec.clear();
  int adoptedID = 0;
  for (std::vector<SortIndex>::iterator it = sortIndexes.begin();
       it != sortIndexes.end();++it) {
    if (it->isAdopted) {
      it->adoptedID = adoptedID;
      calcVec.push_back(it->ipos);
      ++adoptedID;
    }
    else {
      it->adoptedID = -1;
    }
  }

  cellRange.min = 0;
  cellRange.max = 0;
  for (std::vector<SortIndex>::iterator it = sortIndexes.begin();
       it != sortIndexes.end();++it) {
    if (it->isAdopted) {
      for (int i=0;i<3;++i) {
        cellRange.min[i] = std::min(it->ipos[i], cellRange.min[i]);
        cellRange.max[i] = std::max(it->ipos[i], cellRange.max[i]);
      }
    }
  }

  cellPairs.clear();
  for (std::vector<Pair>::iterator it = pairIDs.begin();
       it != pairIDs.end();++it) {
    int icell = sortIndexes.at(it->first).adoptedID;
    int jcell = sortIndexes.at(it->second).adoptedID;
    assert(icell >= 0 && jcell >= 0);
    cellPairs.push_back(Pair(icell, jcell));
  }
}

}
#endif
