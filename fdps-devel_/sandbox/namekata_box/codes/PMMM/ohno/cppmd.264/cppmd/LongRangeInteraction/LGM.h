#ifndef LGM_H
#define LGM_H

#include "ChargeAssignmentMethod.h"

namespace MultigridModule {

class LGM : public ChargeAssignmentMethod {
public:
  typedef LGM* Ptr;

  LGM();
  LGM(double _mgBeta, 
      double _originX, double _originY, double _originZ, 
      double _lgmCACutoff, double _lgmBICutoff) : 
    ChargeAssignmentMethod(_mgBeta, _originX, _originY, _originZ ),
    lgmCACutoff(_lgmCACutoff),lgmCACutoff2(_lgmCACutoff*_lgmCACutoff),
    lgmBICutoff(_lgmBICutoff),lgmBICutoff2(_lgmBICutoff*_lgmBICutoff),
    caCoeff(0), biCoeffFc(0), biCoeffPot(0)
  {}
  virtual ~LGM() {}

  void initialize( const SpaceVector<double>& _side,
		   const SpaceVector<int>& nodeShape,
		   SpaceVector<int>& nGrid,
		   std::vector<int>& localIndex );

  void chargeAssignment(const std::vector<Position>& particleCrd,
			const std::vector<double>& particleCharge,
			ChargeAssignmentMethod::GridAry3& rho,
			const long step);

  void backInterpolation(const std::vector<Position>& particleCrd,
			 const std::vector<double>& particleCharge,
			 std::vector<Force>& force,
			 ChargeAssignmentMethod::GridAry3& rho,
			 ChargeAssignmentMethod::GridAry3& phi,
			 const long step);

  double calcReciprocalPotential(ChargeAssignmentMethod::GridAry3& rho,
				 ChargeAssignmentMethod::GridAry3& phi);

  void getGridRange(std::vector<int>& begin,
		    std::vector<int>& end,
		    const Position& coordinate)
  {
    int myIndex;
    for (int idim=0; idim<3; ++idim) {
      myIndex  = int( floor(coordinate[idim]*hInv[idim]) );
      assert(-numGridMargin <= myIndex && myIndex <= nGrid[idim]+numGridMargin );
      begin[idim] = gridRangeBegin[myIndex+numGridMargin][idim];
      end[idim]   = gridRangeEnd[myIndex+numGridMargin][idim];
      assert(0 <= begin[idim]);
      assert(end[idim] <= localNumGrid[idim]);
    }
  }

  void initialize(const std::vector<Position>& particleCrd);

  void clearTables();

private:

  std::vector<double> expTableX;
  std::vector<double> expTableY;
  std::vector<double> expTableZ;

  typedef struct {
    std::vector<double> exp;
    std::vector<double> distance;
    std::vector<double> distance2;
  } GridTables;

  GridTables xGrid;
  GridTables yGrid;
  GridTables zGrid;

  double lgmCACutoff, lgmCACutoff2;
  double lgmBICutoff, lgmBICutoff2;

  int getNodeNum(int gridNum, int idim){
    return int( floor(gridNum/localNumGrid[idim]) );
  }

  void generateGridRange();
  SpaceVector<int> numGridCutoff;
  std::vector<SpaceVector<int> > gridRangeBegin, gridRangeEnd;
  double caCoeff, biCoeffFc, biCoeffPot;
  int numGridMargin;

};
}
#endif
