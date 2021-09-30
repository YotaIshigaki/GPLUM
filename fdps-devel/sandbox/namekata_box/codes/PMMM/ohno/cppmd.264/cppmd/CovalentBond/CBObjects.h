#ifndef CBOBJECTS_H
#define CBOBJECTS_H

#include "CBInterfaceFwd.h"

#include <cstdio>

namespace CBModule {

class CBObjects {
public:
  CBObjects(CBInterface* _cbInterface) 
    : cbInterface(_cbInterface),
    scnb(2.0), scee(1.2), scnbInv(1.0/scnb), sceeInv(1.0/scee), cancel_short(true), cancel_short_virial(false) 
    {
      //      printf("CBObjects(CBInterface* _cbInterface)\n");
    }
#ifdef ORG_CB
  template<typename PR>
  void calcCBBond(CBInterface_::Context pContext);
  template<typename PR>
  void calcCBAngle(CBInterface_::Context pContext);
  template<typename PR>
  void calcCBTorsion(CBInterface_::Context pContext);
  template<typename PR>
  void calcCBImproper(CBInterface_::Context pContext);
#endif
  void set_cancel_short(){
    cancel_short = true;
  }
  void unset_cancel_short(){
    cancel_short = false;
  }
  void set_cancel_short_virial(){
    cancel_short_virial = true;
  }
  void unset_cancel_short_virial(){
    cancel_short_virial = false;
  }

  // DEBUG LJC cancel 
  double tmpshort;

private:
#ifdef ORG_CB
  template<typename PR>
  void subtractCBParticleForce(CBInterface_::Context pContext,
                               PR pi, 
                               PR pj, 
                               CBInterface_::ForceLocation fi,
                               CBInterface_::ForceLocation fj,
			       const Position& d);
  template<typename PR>
  void subtractCBParticleForce(CBInterface_::Context pContext,
                               PR pi, 
                               PR pj, 
                               CBInterface_::ForceLocation fi,
                               CBInterface_::ForceLocation fj,
			       const Position& d,
			       double& virial);
  template<typename PR>
  void calcCB14Interaction(CBInterface_::Context pContext,
                           PR pi, 
                           PR pl, 
                           CBInterface_::ForceLocation fi,
                           CBInterface_::ForceLocation fl,
			   const Position& d);

  template <typename PR, typename T, typename ITERATOR> 
  void calcCBTorsionImproper(CBInterface_::Context pContext);
#endif

  CBInterface* cbInterface;
  double scnb,scee,scnbInv,sceeInv;
  bool cancel_short;
  bool cancel_short_virial;
};
}
#endif
