#ifndef FFT3DFFTW3_H
#define FFT3DFFTW3_H
#ifdef USE_FFTW3

#include "FFT3D.h"
#include "fftw3.h"

namespace PMEModule {

class FFT3Dfftw3: public FFT3D {
    public:
	FFT3Dfftw3(int nx, int ny, int nz);
	~FFT3Dfftw3();
	void forward();
	void inverse();
	void backward();
	static FFT3D::Ptr createFFT3D(int nx, int ny, int nz) {
	    return FFT3D::Ptr(new FFT3Dfftw3(nx, ny, nz));
	}
    private:
	fftw_complex* top;
	double normalizeFactor;
	fftw_plan forwardPlan;
	fftw_plan inversePlan;
};

FFT3Dfftw3::FFT3Dfftw3(int _nx, int _ny, int _nz):
    FFT3D::FFT3D(_nx, _ny, _nz),
    top((fftw_complex *)FFT3D::pc3d->data()),
    normalizeFactor(1.0/(nx*ny*nz)),
    forwardPlan(fftw_plan_dft_r2c_3d(_nx, _ny, _nz,
	    FFT3D::pr3d->data(), top, FFTW_ESTIMATE)),
    inversePlan(fftw_plan_dft_c2r_3d(_nx, _ny, _nz,
	    top, FFT3D::pr3d->data(), FFTW_ESTIMATE))
{
}

FFT3Dfftw3::~FFT3Dfftw3()
{
    fftw_destroy_plan(forwardPlan);
    fftw_destroy_plan(inversePlan);
    fftw_cleanup();
}

void FFT3Dfftw3::forward()
{
    fftw_execute(forwardPlan);
}

void FFT3Dfftw3::inverse()
{
    fftw_execute(inversePlan);
    scaling(normalizeFactor);
}

void FFT3Dfftw3::backward()
{
    fftw_execute(inversePlan);
}
}
#endif
#endif
