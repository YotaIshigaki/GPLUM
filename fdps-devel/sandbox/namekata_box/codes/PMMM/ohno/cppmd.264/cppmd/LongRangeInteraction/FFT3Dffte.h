#ifndef FFT3DFFTE_H
#define FFT3DFFTE_H
#ifdef USE_FFTE

#include "FFT3D.h"

extern "C" {
    void zfft3d_(std::complex<double>* a, int *nx, int *ny, int *nz, int *iopt);
}

namespace PMEModule {

class FFT3Dffte: public FFT3D {
    public:
	FFT3Dffte(int nx, int ny, int nz);
	~FFT3Dffte() {}
	void forward();
	void inverse();
	void backward();
	static FFT3D::Ptr createFFT3D(int nx, int ny, int nz) {
	    return FFT3D::Ptr(new FFT3Dffte(nx, ny, nz));
	}
	static int getAllowableLength(int n);
	StorageType getComplexStorageType() { return FULL_Z; }
    private:
	double factor;
};

FFT3Dffte::FFT3Dffte(int _nx, int _ny, int _nz):
    FFT3D::FFT3D(_nx, _ny, _nz, 1), factor(nx*ny*nz)
{
    int iopt=0;
    zfft3d_(FFT3D::pc3d->data(), &nz, &ny, &nx, &iopt);
}

void FFT3Dffte::forward()
{
    int iopt=-1;
    Complex3D& c3d = *(FFT3D::pc3d);
    Real3D& r3d = *(FFT3D::pr3d);
    for (int ix=0;ix<nx;++ix) {
	for (int iy=0;iy<ny;++iy) {
	    for (int iz=0;iz<nz;++iz) {
		c3d[ix][iy][iz] = r3d[ix][iy][iz];
	    }
	}
    }
    zfft3d_(c3d.data(), &nz, &ny, &nx, &iopt);
}

void FFT3Dffte::inverse()
{
    int iopt=1;
    Complex3D& c3d = *(FFT3D::pc3d);
    Real3D& r3d = *(FFT3D::pr3d);
    zfft3d_(c3d.data(), &nz, &ny, &nx, &iopt);
    for (int ix=0;ix<nx;++ix) {
	for (int iy=0;iy<ny;++iy) {
	    for (int iz=0;iz<nz;++iz) {
		r3d[ix][iy][iz] = c3d[ix][iy][iz].real();
	    }
	}
    }
}

void FFT3Dffte::backward()
{
    int iopt=1;
    Complex3D& c3d = *(FFT3D::pc3d);
    Real3D& r3d = *(FFT3D::pr3d);
    zfft3d_(c3d.data(), &nz, &ny, &nx, &iopt);
    for (int ix=0;ix<nx;++ix) {
	for (int iy=0;iy<ny;++iy) {
	    for (int iz=0;iz<nz;++iz) {
		r3d[ix][iy][iz] = factor*c3d[ix][iy][iz].real();
	    }
	}
    }
}

int FFT3Dffte::getAllowableLength(int n)
{
    for (;;++n) {
      int m = n;
      while (m % 2 == 0) m /= 2;
      while (m % 3 == 0) m /= 3;
      while (m % 5 == 0) m /= 5;
      if (m == 1) break;
    }
    return n;
}
}
#endif
#endif
