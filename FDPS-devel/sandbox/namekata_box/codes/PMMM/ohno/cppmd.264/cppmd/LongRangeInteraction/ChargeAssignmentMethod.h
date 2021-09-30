#ifndef CHARGEASSIGNMENTMETHOD_H
#define CHARGEASSIGNMENTMETHOD_H

#include <vector>

#include <sys/time.h>

#include "Common.h"
#include "Array3D.h"

namespace MultigridModule {

class ChargeAssignmentMethod{
public:
  typedef ChargeAssignmentMethod* Ptr;
  typedef PMEModule::Array3D<double>::Array GridAry3;

  ChargeAssignmentMethod();
  ChargeAssignmentMethod(double _alpha,
			 double _originX, double _originY, double _originZ) :
    alpha(_alpha),
    mgAlpha (_alpha*sqrt(2.0)),
    mgAlpha2(mgAlpha*mgAlpha),
    mgAlpha3(mgAlpha*mgAlpha*mgAlpha),
    mgAlpha5(mgAlpha*mgAlpha*mgAlpha*mgAlpha*mgAlpha),
    origin(_originX, _originY, _originZ)
  {}
  virtual ~ChargeAssignmentMethod() {}

  virtual void initialize( const SpaceVector<double>& side,
			   const SpaceVector<int>& nodeShape,
			   SpaceVector<int>& nGrid,
			   std::vector<int>& localPosition) = 0;

  virtual void chargeAssignment(const std::vector<Position>& particleCrd,
				const std::vector<double>& particleCharge,
				GridAry3& rho, const long step) = 0;

  virtual void backInterpolation(const std::vector<Position>& particleCrd,
				 const std::vector<double>& particleCharge,
				 std::vector<Force>& force,
				 GridAry3& rho, GridAry3& phi, const long step) = 0;

  virtual double calcReciprocalPotential(GridAry3& rho, GridAry3& phi) = 0;

  double getHx(){return h[0];};
  double getHy(){return h[1];};
  double getHz(){return h[2];};

  int getNx(){return nGrid[0];};
  int getNy(){return nGrid[1];};
  int getNz(){return nGrid[2];};

  // compute self energy
  double calcSelfPotential(const std::vector<double>& particleCharge)
  {
    double selfEnergy = 0.0e0;
    for(std::vector<double>::size_type i=0; i<particleCharge.size(); ++i)
      selfEnergy += particleCharge[i]*particleCharge[i];
    return selfEnergy * (-alpha*M_2_SQRTPI*.5);
  }

  std::vector<int>& getLocalNumGridBegin(){return localNumGridBegin;}
  std::vector<int>& getLocalNumGridEnd(){return localNumGridEnd;}
  std::vector<int>& getLocalNumGrid(){return localNumGrid;}

protected:

  double alpha;
  double mgAlpha, mgAlpha2, mgAlpha3, mgAlpha5;
  SpaceVector<double> origin;

  std::vector<int> localNumGridBegin;
  std::vector<int> localNumGridEnd;
  std::vector<int> localIndex;
  std::vector<int> localNumGrid;

  SpaceVector<double> h;
  SpaceVector<double> hInv;

  SpaceVector<double> side;
  SpaceVector<double> halfSide;

  SpaceVector<int> nGrid;

};
}
#endif
