#ifndef CELLMETHODCORE_H
#define CELLMETHODCORE_H

#include "Common.h"
#include "CellMethodFwd.h"

namespace CellMethodModule {
  typedef SpaceVector<int> Indexes;
  typedef struct {
    Indexes min;
    Indexes max;
  } CellRange;

class CellMethodCore {
public:
  typedef CellMethod_::Context Context;
  typedef std::pair<int,int> Pair;
  CellMethodCore(CellMethod* _method)
    : method(_method), calcVec(), cellPairs(), cellRange(), fullFlags() {}

  virtual ~CellMethodCore(){}

  virtual void assignCalcCell(Context pContext) = 0;
  const std::vector<Indexes>& getCalcVec() { return calcVec; }
  const std::vector<Pair>& getCellPairs() { return cellPairs; }
  const CellRange& getCellRange() { return cellRange; }
  
protected:
  CellMethod* method;
  std::vector<Indexes> calcVec;
  std::vector<Pair> cellPairs;
  CellRange cellRange;
  std::vector<bool> fullFlags;
};
}
#endif
