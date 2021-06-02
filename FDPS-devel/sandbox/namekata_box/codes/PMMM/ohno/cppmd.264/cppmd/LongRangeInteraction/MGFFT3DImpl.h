#ifndef MGFFT3DIMPL_H
#define MGFFT3DIMPL_H
#include "MGFFT3D.h"

namespace MultigridModule {

void copyFromFFT3D(const FFT3D::Ptr fft, double *p)
{
  const FFT3D::Real3D& r3d = fft->getReal();
  int nz = fft->getSize(2);
  for (int ix=0;ix < fft->getSize(0);++ix) {
    for (int iy=0;iy < fft->getSize(1);++iy,p+=nz) {
      const double* pstart = &r3d[ix][iy][0];
      std::copy(pstart, pstart+nz, p);
    }
  }
}

void copyToFFT3D(FFT3D::Ptr fft, const double *p)
{
  FFT3D::Real3D& r3d = fft->getReal();
  int nz = fft->getSize(2);
  for (int ix=0;ix < fft->getSize(0);++ix) {
    for (int iy=0;iy < fft->getSize(1);++iy,p+=nz) {
      std::copy(p, p+nz, &r3d[ix][iy][0]);
    }
  }
}
}
#endif
