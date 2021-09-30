#ifndef OPERATIONSELECTOR_H
#define OPERATIONSELECTOR_H

struct OperationSelector {
  bool doIntegration;
  bool doShortrangecalculation;
  bool doLongrangecalculation;
  bool doCovalentBondcaculation;
  bool doExchangeForce;
  bool doExchangeForceIndexed;
  bool doShortEnergyToHalf;
  bool doShortLongCommunication;
  bool doReturnShortForce;
  bool doMakePairList;
  bool doEnergycalculation;
  bool doVirialcalculation;
  bool doCorrectTranslation;
};

#endif

