#ifndef PMEARRAY_H
#define PMEARRAY_H

#include "UseConfig.h"

#ifdef CPPMD_ENABLE_SIMPLE_FFT
#include "FFT.h"
typedef FFTModule::FFT3D FFT3D;
#else
#include "FFT3D.h"
#endif
#include "Common.h"
#ifdef CPPMD_ENABLE_OLDPME
#include "PMEInterfaceFwd.h"
#endif

namespace PMEModule {

class PMEArray {
public:
  typedef PMEArray* Ptr;

  PMEArray(const FFT3D::Ptr& fft, double alpha,
	   const SpaceVector<std::size_t>& axis);
  PMEArray(const FFT3D::Ptr& fft, double alpha,
	   const SpaceVector<std::size_t>& axis,
	   const std::vector<double>& b2X,
	   const std::vector<double>& b2Y,
	   const std::vector<double>& b2Z);
  ~PMEArray() {
    delete(arrayPtr);
  }
  
#ifdef CPPMD_ENABLE_OLDPME
  const FFT3D::Real3D& calculate(PMEInterface* pmeInterface,
                                 PMEInterface_::Context pContext);
#endif
  const FFT3D::Real3D& calculate(const SpaceVector<double> boxSize);
  const FFT3D::Real3D& getvirialarray();
private:
#ifdef CPPMD_ENABLE_OLDPME
  bool checkUpdate(PMEInterface* pmeInterface, PMEInterface_::Context pContext);
  void update(PMEInterface* pmeInterface, PMEInterface_::Context pContext);
#endif
  bool checkUpdate(const SpaceVector<double> boxSize);
  void update(const SpaceVector<double> boxSize);

  static int minByAbs(int a, int b) {
    return (abs(a)<abs(b)?(a):(b));
  }

  FFT3D::Real3D_Ptr arrayPtr;
  FFT3D::Real3D_Ptr virialarrayPtr;
  FFT3D::Ptr fft;
  SpaceVector<std::size_t> nGrid;
  double alpha;
  SpaceVector<std::size_t> axis;
  bool useb2;
  std::vector<double> b2[3];
  SpaceVector<double> side;
};
}
#endif
