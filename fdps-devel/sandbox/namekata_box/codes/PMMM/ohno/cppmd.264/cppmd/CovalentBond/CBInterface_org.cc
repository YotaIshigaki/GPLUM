#include "CBInterface.h"

class CBInterface {
public:

  template<typename PR>
    void calcBond(const BondArray& bond, double& energy, double& shortenergy, double& shortvirial){
    CBContext context(bond, energy, shortenergy, shortvirial);
    cbObjects.calcCBBond<PR>(&context);
  }

  template<typename PR>
    void calcAngle(const AngleArray& angle, double& energy, double& shortenergy, double& shortvirial){
    CBContext context(angle, energy, shortenergy, shortvirial);
    cbObjects.calcCBAngle<PR>(&context);
  }

  template<typename PR>
    void calcTorsion(const TorsionArray& torsion, double& energy, double& shortenergy, double& shortvirial){
    CBContext context(torsion, energy, shortenergy, shortvirial);
    cbObjects.calcCBTorsion<PR>(&context);
  }

  template<typename PR>
    void calcImproper(const ImproperArray& improper, double& energy, double& shortenergy, double& shortvirial){
    CBContext context(improper, energy, shortenergy, shortvirial);
    cbObjects.calcCBImproper<PR>(&context);
  }
};

}
