#ifndef PARALLELCOORDINATOR
#define PARALLELCOORDINATOR

#include "OperationSelector.h"

class ParallelCoordinator {
public:

  virtual ~ParallelCoordinator() {}

  virtual OperationSelector selectOperations(int unit_id, 
                                             const bool withlong,
                                             const bool withbond) 
  {
    OperationSelector ope;
#if 1
    ope.doShortEnergyToHalf = false;        // ???
    ope.doShortLongCommunication = false;   // ???
#endif
    if(unit_id&0x1){
      ope.doIntegration = false;
      ope.doShortrangecalculation = false;
      ope.doLongrangecalculation = true;
      ope.doCovalentBondcaculation = false;
    }else{
      ope.doIntegration = true;
      ope.doShortrangecalculation = true;
      if(withlong){
        ope.doLongrangecalculation = true;
      }else{
        ope.doLongrangecalculation = false;
      }
      if(withbond){
        ope.doCovalentBondcaculation = true;
      }else{
        ope.doCovalentBondcaculation = false;
      }
    }
    ope.doExchangeForce = true;
#ifdef USE_PAIRLIST
    ope.doMakePairList = true;
#else
    ope.doMakePairList = false;
#endif
    ope.doReturnShortForce = true;
    ope.doExchangeForceIndexed = false;
    ope.doEnergycalculation = true;
    ope.doVirialcalculation = false;
    return ope;
  }
};

#endif
