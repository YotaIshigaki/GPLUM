#ifndef HALFSHELLCORE_H
#define HALFSHELLCORE_H

#include "CellMethodCore.h"

namespace CellMethodModule {

class HalfShellCore : public CellMethodCore {
public:
  HalfShellCore(CellMethod* method, bool _fullFlag=false)
    : CellMethodCore(method), fullFlag(_fullFlag) {}

  void assignCalcCell(Context pContext);
private:
  bool fullFlag;
};
}
#endif
