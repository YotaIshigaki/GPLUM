#ifndef EWALDBASE_H
#define EWALDBASE_H

#include "SurfaceDipole.h"

namespace EwaldModule {

class EwaldBase {
public:
  EwaldBase(const double _cutoff,
            const double _alpha=0.0, const double _kCutoff=0.0,
	    const int surfaceDipole=0)
    : cutoff(_cutoff), alpha(_alpha), kCutoff(_kCutoff),
      calcSurfaceDipole(surfaceDipole), surface(),
      selfEnergy(0.0), shareValue(1.0),
      erfcCutValue(1e-6),       /* erfc is regarded as 0 below this value (auto cut_off) */
      ewaldExpCutValue(5e-9)    /* exp  is regarded as 0 below this value (auto cut_off) */
  {
    //    std::cout << " EwaldBase alpha " << alpha << "  kCutoff " << kCutoff <<  std::endl;
  }

  virtual ~EwaldBase(){}

  virtual void setSide(double side);
  double getAlpha() const { return alpha;}
  double getKCutoff() const { return kCutoff;}
  int getSurfaceDipole() const { return calcSurfaceDipole;}
  void calculateSelfEnergy(const std::vector<double>& charge) {
    for (std::vector<double>::const_iterator it = charge.begin();
         it != charge.end();++it) {
      selfEnergy += (-(*it)*(*it)*getAlpha()*M_2_SQRTPI*.5);
    }
printf("selfEnergy %e\n",selfEnergy);
  }
  void clearDipoleMoment(double volume) {
    if (calcSurfaceDipole) {
      surface.initialize(volume);
    }
  }
  void calculateDipoleMoment(const std::vector<Position>& cd,
                        const std::vector<double>& charge)
  {
    if (calcSurfaceDipole) {
       surface.addDipoleMoment(cd, charge);
    }
  }
  void calculateDipoleForce(std::vector<Force>& fc, 
                       const std::vector<double>& charge)
  {
    if (calcSurfaceDipole) {
       surface.addDipoleForce(fc, charge);
    }
  }
  SpaceVector<double>& getDipoleMoment() { return surface.dipoleMoment; }
  double getDipoleEnergy() { return shareValue*surface.dipoleEnergy(); }
  double getSelfEnergy() { return selfEnergy; }
  void setShare(double _shareValue) { shareValue = _shareValue; }
  double getShare() { return shareValue; }

private:
  double cutoff;
  double alpha;
  double kCutoff;
  int calcSurfaceDipole;
  SurfaceDipole surface;
  double selfEnergy;
  double shareValue;
  double erfcCutValue;
  double ewaldExpCutValue;
};

void estimate_alpha_kCutoff(double &alpha, double &kCutoff,
                            double side, double cutoff, 
                            double erfcCutValue=1.0e-6,
                            double ewaldExpCutValue=5.0e-9);
double estimate_alpha(double cutoff, double erfcCutValue=1.0e-6);
double estimate_kCutoff(double side, double cutoff, double alpha);
}
#endif
