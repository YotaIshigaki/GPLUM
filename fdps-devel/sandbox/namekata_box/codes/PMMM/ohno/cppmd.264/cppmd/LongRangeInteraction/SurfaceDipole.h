#ifndef SURFACEDIPOLE_H
#define SURFACEDIPOLE_H

#include <vector>
#include "Common.h"

namespace EwaldModule {

class SurfaceDipole {
public:
  void initialize(double volume) {
    dipoleMoment = 0.0;
    dipoleFactor = 2.0*M_PI/(3.0*volume);
  }
  void addDipoleMoment(const std::vector<Position>& cd,
                       const std::vector<double>& charge) {
    for (std::vector<Position>::size_type i=0;i<cd.size();++i) {
      dipoleMoment += charge[i]*cd[i];
    }
  }  
  
  template<class PA>
  void addDipoleMoment(const PA& particlearray,
		       const std::vector<ParticleRange>& particlerange) {

    int pr;
    for(pr=0;pr<particlerange.size();pr++){
      int i;
      int si=particlerange[pr].begin;
      int ei=particlerange[pr].end;
      for(i=si;i<ei;i++){
	dipoleMoment += getcharge(particlearray,i)*getpos(particlearray,i);
      }
    }
  } 
  
  double dipoleEnergy() {
    return dipoleFactor*dipoleMoment.norm2();
  }
  void addDipoleForce(std::vector<Force>& fc,
                      const std::vector<double>& charge) const {
    for (std::vector<Force>::size_type i=0;i<fc.size();++i) {
      fc[i] -= 2.0*dipoleFactor*charge[i]*dipoleMoment;
    }
  }
  
  SpaceVector<double> dipoleMoment;
  double dipoleFactor;
};
}
#endif
