#ifndef HALFSHELLCOREIMPL_H
#define HALFSHELLCOREIMPL_H

#include "CellMethodCore.h"

namespace CellMethodModule {

void HalfShellCore::assignCalcCell(Context pContext)
{
  CellRange checkRange;
  Indexes cutoffCellNum = method->getCutoffCellNum(pContext);
  checkRange.min = -cutoffCellNum;
  checkRange.max = cutoffCellNum;
  calcVec.clear();
  fullFlags.clear();
  Indexes zero(0);
  calcVec.push_back(zero);
  fullFlags.push_back(0);
  for (int ix = checkRange.min.x;ix <= checkRange.max.x;++ix) {
    for (int iy = checkRange.min.y;iy <= checkRange.max.y;++iy) {
      for (int iz = checkRange.min.z;iz <= checkRange.max.z;++iz) {
        Indexes ipos(ix, iy, iz);
	if (method->inCutoffSphere(pContext, ipos)) {
          if (ix > 0 || (ix == 0 && iy > 0) ||
             (ix == 0 && iy == 0 && iz > 0)) {
            calcVec.push_back(ipos);
            fullFlags.push_back(0);
          }
          else if (fullFlag && ipos != zero) {
            calcVec.push_back(ipos);
            fullFlags.push_back(1);
          }
        }
      }
    }
  }

  cellRange.min = 0;
  cellRange.max = 0;
  cellPairs.clear();
  for (std::vector<Indexes>::size_type i = 0;i < calcVec.size();++i) {
    const Indexes& ipos = calcVec[i];
    for (int j=0;j<3;++j) {
      cellRange.min[j] = std::min(ipos[j], cellRange.min[j]);
      cellRange.max[j] = std::max(ipos[j], cellRange.max[j]);
    }
    if (!fullFlags[i]) {
      cellPairs.push_back(Pair(0, i));
    }
  }
}

}
#endif
