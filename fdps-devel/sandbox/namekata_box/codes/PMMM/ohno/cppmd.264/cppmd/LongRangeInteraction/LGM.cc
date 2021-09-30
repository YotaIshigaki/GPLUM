#include <iostream>
#include "LGM.h"

using namespace std;
using namespace MultigridModule;

void LGM::initialize( const SpaceVector<double>& _side,
		      const SpaceVector<int>& nodeShape,
		      SpaceVector<int>& _nGrid,
		      std::vector<int>& _localIndex )
{
  double lgmMaxCutoff = max(lgmCACutoff, lgmBICutoff);

  nGrid( _nGrid[0], _nGrid[1], _nGrid[2] ); 

  side ( _side.getX(),
         _side.getY(),
         _side.getZ() );

  halfSide ( _side.getX()*.5,
             _side.getY()*.5,
             _side.getZ()*.5 );

  h ( side[0]/static_cast<double>(nGrid[0]),
      side[1]/static_cast<double>(nGrid[1]),
      side[2]/static_cast<double>(nGrid[2]) );

  hInv ( 1.0e0/h[0], 1.0e0/h[1], 1.0e0/h[2] );

  numGridCutoff ( static_cast<int>( ceil( lgmMaxCutoff/h[0] ) ),
		  static_cast<int>( ceil( lgmMaxCutoff/h[1] ) ),
		  static_cast<int>( ceil( lgmMaxCutoff/h[2] ) ) );

  double margin = side[0];
  numGridMargin = 0;

  for (int idim=0; idim<3; ++idim)
    localIndex.push_back( _localIndex[idim] );

  for (int idim=0; idim<3; ++idim) {

    numGridMargin = max(numGridMargin, static_cast<int>(margin/h[idim])*2 );

    int numGrid = nGrid[idim]/nodeShape[idim];
    int numGridIndexBegin = numGrid * localIndex[idim];
    int numGridIndexEnd = numGrid*( localIndex[idim]+1 );

    localNumGridBegin.push_back( numGridIndexBegin );
    localNumGridEnd.push_back( numGridIndexEnd );
    localNumGrid.push_back( numGridIndexEnd - numGridIndexBegin );

  }

  // generate calculation grid range list
  generateGridRange();

  // set coefficient chargeassignment and backInterpolation
  caCoeff    = mgAlpha3/pow(M_PI,1.5);
  biCoeffFc  = -2.0*mgAlpha5*h[0]*h[1]*h[2]/pow(M_PI,1.5);
  biCoeffPot = h[0]*h[1]*h[2]*.5;

  expTableX.resize(localNumGrid[0]);
  for(int igrid=0; igrid<localNumGrid[0]; ++igrid){
    expTableX[igrid]=0;
  }

  expTableY.resize(localNumGrid[1]);
  for(int jgrid=0; jgrid<localNumGrid[1]; ++jgrid){
    expTableY[jgrid]=0;
  }

  expTableZ.resize(localNumGrid[2]);
  for(int kgrid=0; kgrid<localNumGrid[2]; ++kgrid){
    expTableY[kgrid]=0;
  }

};

void LGM::chargeAssignment(const std::vector<Position>& particleCrd,
			   const std::vector<double>& particleCharge,
			   ChargeAssignmentMethod::GridAry3& rho,
			   const long step)
{

  // number of atoms
  int particleNum = particleCharge.size();

  double norm2;
  double pcx, pcy, pcz;
  double distanceX, distanceY, distanceZ;
  double distanceX2, distanceY2, distanceZ2;

  std::vector<int> begin(3), end(3);

  for(int icharge=0; icharge<particleNum; ++icharge){

    // particle coordinate
    pcx = particleCrd[icharge][0];
    pcy = particleCrd[icharge][1];
    pcz = particleCrd[icharge][2];

    getGridRange(begin, end, particleCrd[icharge]);

    double xgridBegin = origin[0] + 
      h[0]*static_cast<double>( begin[0] + localNumGridBegin[0] - 1) - pcx;
    double ygridBegin = origin[1] + 
      h[1]*static_cast<double>( begin[1] + localNumGridBegin[1] - 1) - pcy;
    double zgridBegin = origin[2] + 
      h[2]*static_cast<double>( begin[2] + localNumGridBegin[2] - 1) - pcz;

    // x-grid
    distanceX = xgridBegin;
    for(int igrid=begin[0]; igrid<end[0]; ++igrid){
      distanceX += h[0];
      expTableX[igrid] = exp(-mgAlpha2*distanceX*distanceX);
    }
    // y-grid
    distanceY = ygridBegin;
    for(int jgrid=begin[1]; jgrid<end[1]; ++jgrid){
      distanceY += h[1];
      expTableY[jgrid] = exp(-mgAlpha2*distanceY*distanceY);
    }
    // z-grid
    distanceZ = zgridBegin;
    for(int kgrid=begin[2]; kgrid<end[2]; ++kgrid){
      distanceZ += h[2];
      expTableZ[kgrid] = exp(-mgAlpha2*distanceZ*distanceZ);
    }

    // charge assignment
    distanceX = xgridBegin;
    for(int igrid=begin[0]; igrid<end[0]; ++igrid){
      distanceX += h[0];
      distanceX2 = distanceX*distanceX;

      distanceY = ygridBegin;
      for(int jgrid=begin[1]; jgrid<end[1]; ++jgrid){
	distanceY += h[1];
	distanceY2 = distanceY*distanceY;

	distanceZ = zgridBegin;
	for(int kgrid=begin[2]; kgrid<end[2]; ++kgrid){
	  distanceZ += h[2];
	  distanceZ2 = distanceZ*distanceZ;

	  norm2 = distanceX2 + distanceY2 + distanceZ2;

	  // cutoff
	  if( norm2 < lgmCACutoff2 ){
	    rho[igrid][jgrid][kgrid] += caCoeff * particleCharge[icharge] *
	      expTableX[igrid] * expTableY[jgrid] * expTableZ[kgrid];
	  }

	}  // loop kgrid
      }  // loop jgrid
    }  // loop igrid
  }  // particle loop
}

void LGM::backInterpolation(const std::vector<Position>& particleCrd,
			    const std::vector<double>& particleCharge,
			    std::vector<Force>& force,
			    ChargeAssignmentMethod::GridAry3& rho,
			    ChargeAssignmentMethod::GridAry3& phi,
			    const long step)
{
  // number of atoms
  int particleNum = particleCharge.size();

  double norm2;
  double pcx, pcy, pcz;
  double fix, fiy, fiz;

  // tmp value
  double cf;
  std::vector<int> begin(3), end(3);

  double distanceX, distanceY, distanceZ;
  double distanceX2, distanceY2, distanceZ2;

  // number of atoms
  for(int icharge=0; icharge<particleNum; ++icharge){

    // particle coordinate
    pcx = particleCrd[icharge][0];
    pcy = particleCrd[icharge][1];
    pcz = particleCrd[icharge][2];

    fix = 0.0; fiy = 0.0; fiz = 0.0;

    getGridRange(begin, end, particleCrd[icharge]);

    double xgridBegin = origin[0] + 
      h[0]*static_cast<double>( begin[0] + localNumGridBegin[0] - 1) - pcx;
    double ygridBegin = origin[1] + 
      h[1]*static_cast<double>( begin[1] + localNumGridBegin[1] - 1) - pcy;
    double zgridBegin = origin[2] + 
      h[2]*static_cast<double>( begin[2] + localNumGridBegin[2] - 1) - pcz;

    // x-grid
    distanceX = xgridBegin;
    for(int igrid=begin[0]; igrid<end[0]; ++igrid){
      distanceX += h[0];
      expTableX[igrid] = exp(-mgAlpha2*distanceX*distanceX);
    }
    // y-grid
    distanceY = ygridBegin;
    for(int jgrid=begin[1]; jgrid<end[1]; ++jgrid){
      distanceY += h[1];
      expTableY[jgrid] = exp(-mgAlpha2*distanceY*distanceY);
    }
    // z-grid
    distanceZ = zgridBegin;
    for(int kgrid=begin[2]; kgrid<end[2]; ++kgrid){
      distanceZ += h[2];
      expTableZ[kgrid] = exp(-mgAlpha2*distanceZ*distanceZ);
    }

    // back interpolation
    distanceX = xgridBegin;
    for(int igrid=begin[0]; igrid<end[0]; ++igrid){
      distanceX += h[0];
      distanceX2 = distanceX*distanceX;

      distanceY = ygridBegin;
      for(int jgrid=begin[1]; jgrid<end[1]; ++jgrid){
	distanceY += h[1];
	distanceY2 = distanceY*distanceY;

	distanceZ = zgridBegin;
	for(int kgrid=begin[2]; kgrid<end[2]; ++kgrid){
	  distanceZ += h[2];
	  distanceZ2 = distanceZ*distanceZ;

	  norm2 = distanceX2 + distanceY2 + distanceZ2;

	  // cutoff
	  if( norm2 < lgmBICutoff2 ){
	    cf = phi[igrid][jgrid][kgrid]*expTableX[igrid]*expTableY[jgrid]*expTableZ[kgrid];
	    fix += distanceX*cf;
	    fiy += distanceY*cf;
	    fiz += distanceZ*cf;
	  }
	} // loop kgrid
      } // loop jgrid
    } // loop igrid

    force[icharge][0] += biCoeffFc*particleCharge[icharge]*fix;
    force[icharge][1] += biCoeffFc*particleCharge[icharge]*fiy;
    force[icharge][2] += biCoeffFc*particleCharge[icharge]*fiz;
    
  } // particle loop

}

double LGM::calcReciprocalPotential(ChargeAssignmentMethod::GridAry3& rho,
				    ChargeAssignmentMethod::GridAry3& phi)
{
  double recPot = 0.0e0;

  for(int i=0; i<localNumGrid[0]; ++i)
    for(int j=0; j<localNumGrid[1]; ++j)
      for(int k=0; k<localNumGrid[2]; ++k)
	recPot += rho[i][j][k]*phi[i][j][k];

  return biCoeffPot*recPot;
}

void LGM::generateGridRange()
{
  std::vector<int> gridNum(3);
  gridNum[0] = nGrid[0] + 2*numGridMargin + 1;
  gridNum[1] = nGrid[1] + 2*numGridMargin + 1;
  gridNum[2] = nGrid[2] + 2*numGridMargin + 1;

  int maxGridNum = 0;
  for(int idim=0; idim<3; ++idim)
    maxGridNum = max(maxGridNum,gridNum[idim]);

  gridRangeBegin.resize( maxGridNum );
  gridRangeEnd.  resize( maxGridNum );

  for(int idim=0; idim<3; ++idim){

    for(int igrid=-numGridMargin; igrid<nGrid[idim]+numGridMargin; ++igrid){

      int begin = igrid - numGridCutoff[idim] + 1;
      int end   = igrid + numGridCutoff[idim] + 1;
      int mgrid   = igrid + numGridMargin;

      gridRangeBegin[mgrid][idim] = max(localNumGridBegin[idim],begin);
      gridRangeEnd  [mgrid][idim] = min(localNumGridEnd[idim],    end);

      gridRangeBegin[mgrid][idim] -= localNumGridBegin[idim];
      gridRangeEnd  [mgrid][idim] -= localNumGridBegin[idim];

    } // igrid
  } // idim
}
