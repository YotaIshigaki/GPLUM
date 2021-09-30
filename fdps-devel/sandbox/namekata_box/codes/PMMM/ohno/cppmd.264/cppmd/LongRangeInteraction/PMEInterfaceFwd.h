#ifndef PMEINTERFACEFWD_H
#define PMEINTERFACEFWD_H

namespace EwaldModule {
class EwaldBaseContext;
}
namespace PMEModule {
typedef EwaldModule::EwaldBaseContext PMEContext;
class PMEInterface;
namespace PMEInterface_ {
  typedef PMEContext* Context;
}
}
#endif
