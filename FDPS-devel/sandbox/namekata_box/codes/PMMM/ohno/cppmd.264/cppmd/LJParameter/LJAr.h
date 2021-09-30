#include "UnitParameter.h"
#include "LJParameter.h"

static LJMixparameterArray ljmixAr;

class LJAr : public LJParameter {
public:
  double LJCutoffEnergyCorrection;  // It may be onley Ar-Ar interaction
  LJAr () : LJParameter("Ar",
                        UnitParameter::normalizeLength(3.4e-10),
                        UnitParameter::normalizeEnergy_eV(120*UnitParameter::Boltzmann))
  {
    ljmixAr.resize2d(1,1);
    LJMixParameter armix;
    LJParameter mix = LJParameterStuff::mixParameter(*this, *this);
    armix.potFirst = mix.getPotentialFirst();
    armix.potSecond = mix.getPotentialSecond();
    armix.forceFirst = mix.getForceFirst();
    armix.forceSecond = mix.getForceSecond();
    ljmixAr[0][0] = armix;
  }

  double calcLJCutoffEnergyCorrection(double cutoff)
  {
    // It may be onley Ar-Ar interaction
    LJCutoffEnergyCorrection = 8./9.*M_PI*(pow(cutoff,-9) - 3.*pow(cutoff,-3));
    return LJCutoffEnergyCorrection;
  }
};
