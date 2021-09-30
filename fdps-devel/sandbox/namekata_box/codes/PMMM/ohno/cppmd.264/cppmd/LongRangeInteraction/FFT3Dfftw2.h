#ifndef FFT3DFFTW2_H
#define FFT3DFFTW2_H
#ifdef USE_FFTW2

#include "FFT3D.h"
#include "rfftw.h"

namespace PMEModule {

class FFT3Dfftw2: public FFT3D {
    public:
	FFT3Dfftw2(int nx, int ny, int nz);
	~FFT3Dfftw2();
	void forward();
	void inverse();
	void backward();
	static FFT3D::Ptr createFFT3D(int nx, int ny, int nz) {
	    return FFT3D::Ptr(new FFT3Dfftw2(nx, ny, nz));
	}
    private:
	fftw_complex* top;
	double normalizeFactor;
	rfftwnd_plan forwardPlan;
	rfftwnd_plan inversePlan;
};

FFT3Dfftw2::FFT3Dfftw2(int _nx, int _ny, int _nz):
    FFT3D::FFT3D(_nx, _ny, _nz),
    top((fftw_complex *)FFT3D::pc3d->data()), normalizeFactor(1.0/(nx*ny*nz)),
    forwardPlan(rfftw3d_create_plan(_nx, _ny, _nz, FFTW_REAL_TO_COMPLEX,
	    FFTW_ESTIMATE | FFTW_IN_PLACE)),
    inversePlan(rfftw3d_create_plan(_nx, _ny, _nz, FFTW_COMPLEX_TO_REAL,
	    FFTW_ESTIMATE | FFTW_IN_PLACE))
{
}

FFT3Dfftw2::~FFT3Dfftw2()
{
    rfftwnd_destroy_plan(forwardPlan);
    rfftwnd_destroy_plan(inversePlan);
}

void FFT3Dfftw2::forward()
{
  //  std::cout << "FFT3Dfftw2::forward" << std::endl;
    rfftwnd_one_real_to_complex(forwardPlan, FFT3D::pr3d->data(), 0);
}

void FFT3Dfftw2::inverse()
{
    rfftwnd_one_complex_to_real(inversePlan, top, 0);
    scaling(normalizeFactor);
}
void FFT3Dfftw2::backward()
{
    rfftwnd_one_complex_to_real(inversePlan, top, 0);
}
}
#endif
#endif
