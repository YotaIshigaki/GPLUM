#ifdef CPPMD_ENABLE_OLDPME
#include "PMEInterface.h"
#endif
#include "PMEArray.h"
#include <cmath>

using namespace PMEModule;

PMEArray::PMEArray(const FFT3D::Ptr& _fft, double _alpha,
                   const SpaceVector<std::size_t>& _axis):
    arrayPtr(), virialarrayPtr(), fft(_fft),
    nGrid(fft->getSize(fft->getComplexAxis(0)),
          fft->getSize(fft->getComplexAxis(1)),
          fft->getSize(fft->getComplexAxis(2))),
    alpha(_alpha),
    axis(_axis[fft->getComplexAxis(0)],
         _axis[fft->getComplexAxis(1)],
         _axis[fft->getComplexAxis(2)]),
    useb2(false), b2(),
    side(0.0)
{
  arrayPtr = FFT3D::createReal3D(nGrid[0], nGrid[1], nGrid[2]);
  virialarrayPtr = FFT3D::createReal3D(nGrid[0], nGrid[1], nGrid[2]);
}

PMEArray::PMEArray(const FFT3D::Ptr& _fft, double _alpha,
                   const SpaceVector<std::size_t>& _axis,
                   const std::vector<double>& _b2X,
                   const std::vector<double>& _b2Y,
                   const std::vector<double>& _b2Z):
    arrayPtr(), virialarrayPtr(), fft(_fft),
    nGrid(fft->getSize(fft->getComplexAxis(0)),
          fft->getSize(fft->getComplexAxis(1)),
          fft->getSize(fft->getComplexAxis(2))),
    alpha(_alpha),
    axis(_axis[fft->getComplexAxis(0)],
         _axis[fft->getComplexAxis(1)],
         _axis[fft->getComplexAxis(2)]),
    useb2(true), b2(),
    side(0.0)
{
  std::vector<double>& b2X = b2[fft->getComplexAxis(0)];
  std::vector<double>& b2Y = b2[fft->getComplexAxis(1)];
  std::vector<double>& b2Z = b2[fft->getComplexAxis(2)];
  b2X.insert(b2X.end(), _b2X.begin(), _b2X.end());
  b2Y.insert(b2Y.end(), _b2Y.begin(), _b2Y.end());
  b2Z.insert(b2Z.end(), _b2Z.begin(), _b2Z.end());
  arrayPtr = FFT3D::createReal3D(nGrid[0], nGrid[1], nGrid[2]);
  virialarrayPtr = FFT3D::createReal3D(nGrid[0], nGrid[1], nGrid[2]);
}

#ifdef CPPMD_ENABLE_OLDPME
const FFT3D::Real3D& PMEArray::calculate(PMEInterface* pmeInterface,
                                         PMEInterface::Context pContext)
{
  //  std::cout << "PMEArray::calculate" << std::endl;
  if (checkUpdate(pmeInterface, pContext)) {
    //    std::cout << "PMEArray::calculate call update" << std::endl;
    update(pmeInterface, pContext);
  }
  FFT3D::Real3D& BC = *arrayPtr;
  //  std::cout << "PMEArray::calculate return *" << arrayPtr << std::endl;
  return BC;
}
#endif

const FFT3D::Real3D& PMEArray::getvirialarray()
{
  FFT3D::Real3D& VBC = *virialarrayPtr;
  return VBC;
}

#ifdef CPPMD_ENABLE_OLDPME
bool PMEArray::checkUpdate(PMEInterface* pmeInterface,
                           PMEInterface::Context pContext)
{
  bool updateFlag = false;
  //  std::cout << "PMEArray::checkUpdate" << std::endl;
  for (int i=0;i<SpaceVector<double>::Dim;++i) {
    if (side[i] != pmeInterface->getSide(pContext)[i]) {
      side[i] = pmeInterface->getSide(pContext)[i];
      updateFlag = true;
    }
  }
  return updateFlag;
}

void PMEArray::update(PMEInterface* pmeInterface,
                      PMEInterface::Context pContext)
{
  //std::cout << "PMEArray::update " << std::endl;
  FFT3D::Real3D& BC = *arrayPtr;
  FFT3D::Real3D& VBC = *virialarrayPtr;
  const double piVInv = 1.0/(M_PI*pmeInterface->getVolume(pContext));
  const double pi2_alpha2Inv = M_PI*M_PI/(alpha*alpha);
  const std::vector<double>& b2X = b2[0];
  const std::vector<double>& b2Y = b2[1];
  const std::vector<double>& b2Z = b2[2];
  SpaceVector<int> mv;
  for (int mx=fft->getComplexStart(0);mx<fft->getComplexEnd(0);++mx) {
    mv[axis[0]] = minByAbs(mx, mx-nGrid.getX());
    for (int my=fft->getComplexStart(1);my<fft->getComplexEnd(1);++my) {
      mv[axis[1]] = minByAbs(my, my-nGrid.getY());
      for (int mz=fft->getComplexStart(2);mz<fft->getComplexEnd(2);++mz) {
        mv[axis[2]] = minByAbs(mz, mz-nGrid.getZ());
        double m2 = pmeInterface->getWaveVector(pContext, mv).norm2();
        if (m2 > 0.0) {
          BC[mx][my][mz] = piVInv*exp(-m2*pi2_alpha2Inv)/m2;   //  Ulrich Essmann 1995  (3.9) (4.7)
	  VBC[mx][my][mz] = piVInv*exp(-m2*pi2_alpha2Inv)*pi2_alpha2Inv;   //  Ulrich Essmann 1995  (3.9) (4.7)
          if (useb2) {  /// Smooth PME
            BC[mx][my][mz] *= b2X[mx]*b2Y[my]*b2Z[mz];  //  Ulrich Essmann 1995  (4.7) (4.8)
            VBC[mx][my][mz] *= b2X[mx]*b2Y[my]*b2Z[mz];  //  Ulrich Essmann 1995  (4.7) (4.8)

          }
        }
        else {
          BC[mx][my][mz] = 0.0;
          VBC[mx][my][mz] = 0.0;
        }
      }
    }
  }
}
#endif

const FFT3D::Real3D& PMEArray::calculate(const SpaceVector<double> boxSize)
{
  //  std::cout << "PMEArray::calculate" << std::endl;
  if (checkUpdate(boxSize)) {
    //    std::cout << "PMEArray::calculate call update" << std::endl;
    update(boxSize);
  }
  FFT3D::Real3D& BC = *arrayPtr;
  //  std::cout << "PMEArray::calculate return *" << arrayPtr << std::endl;
  return BC;
}

bool PMEArray::checkUpdate(const SpaceVector<double> boxSize)
{
  bool updateFlag = false;
  //  std::cout << "PMEArray::checkUpdate" << std::endl;
  for (int i=0;i<SpaceVector<double>::Dim;++i) {
    if (side[i] != boxSize[i]) {
      side[i] =boxSize[i];
      updateFlag = true;
    }
  }
  return updateFlag;
}

#define USE_VIRIAL
void PMEArray::update(const SpaceVector<double> boxSize)
{
  //std::cout << "PMEArray::update " << std::endl;
  FFT3D::Real3D& BC = *arrayPtr;
#ifdef USE_VIRIAL
  FFT3D::Real3D& VBC = *virialarrayPtr;
#endif
  const double piVInv = 1.0/(M_PI*boxSize[0]*boxSize[1]*boxSize[2]);
  const double pi2_alpha2Inv = M_PI*M_PI/(alpha*alpha);
  const std::vector<double>& b2X = b2[0];
  const std::vector<double>& b2Y = b2[1];
  const std::vector<double>& b2Z = b2[2];
  SpaceVector<int> mv;
  for (int mx=fft->getComplexStart(0);mx<fft->getComplexEnd(0);++mx) {
    mv[axis[0]] = minByAbs(mx, mx-nGrid.getX());
    for (int my=fft->getComplexStart(1);my<fft->getComplexEnd(1);++my) {
      mv[axis[1]] = minByAbs(my, my-nGrid.getY());
      for (int mz=fft->getComplexStart(2);mz<fft->getComplexEnd(2);++mz) {
        mv[axis[2]] = minByAbs(mz, mz-nGrid.getZ());
        double m2 = (mv[0]/boxSize[0])*(mv[0]/boxSize[0])
	  + (mv[1]/boxSize[1])*(mv[1]/boxSize[1])
	  + (mv[2]/boxSize[2])*(mv[2]/boxSize[2]);
        if (m2 > 0.0) {
          BC[mx][my][mz] = piVInv*exp(-m2*pi2_alpha2Inv)/m2;   //  Ulrich Essmann 1995  (3.9) (4.7)
#ifdef USE_VIRIAL
	  VBC[mx][my][mz] = piVInv*exp(-m2*pi2_alpha2Inv)*pi2_alpha2Inv;   //  Ulrich Essmann 1995  (3.9) (4.7)
#endif
          if (useb2) {  /// Smooth PME
            BC[mx][my][mz] *= b2X[mx]*b2Y[my]*b2Z[mz];  //  Ulrich Essmann 1995  (4.7) (4.8)
#ifdef USE_VIRIAL
	    VBC[mx][my][mz] *= b2X[mx]*b2Y[my]*b2Z[mz];  //  Ulrich Essmann 1995  (4.7) (4.8)
#endif
          }
        }
        else {
          BC[mx][my][mz] = 0.0;
#ifdef USE_VIRIAL
	  VBC[mx][my][mz] = 0.0;
#endif
        }
      }
    }
  }
}

