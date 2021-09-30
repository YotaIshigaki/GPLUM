#ifndef FFT3DMKL_H
#define FFT3DMKL_H
#ifdef USE_MKL

#include "FFT3D.h"
#include "mkl_dfti.h"

namespace PMEModule {
class FFT3Dmkl: public FFT3D {
    public:
	FFT3Dmkl(int nx, int ny, int nz);
	~FFT3Dmkl();
	void forward();
	void inverse();
	void backward();
	static FFT3D::Ptr createFFT3D(int nx, int ny, int nz) {
	    return FFT3D::Ptr(new FFT3Dmkl(nx, ny, nz));
	}
    private:
	DFTI_DESCRIPTOR_HANDLE forwardHandle;
	DFTI_DESCRIPTOR_HANDLE inverseHandle;
	DFTI_DESCRIPTOR_HANDLE backwardHandle;
};

FFT3Dmkl::FFT3Dmkl(int _nx, int _ny, int _nz):
    FFT3D::FFT3D(_nx, _ny, _nz, 2),
    forwardHandle(0),inverseHandle(0),backwardHandle(0)
{
    MKL_LONG length[3];
    MKL_LONG status;
    length[0] = _nx;
    length[1] = _ny;
    length[2] = _nz;
    MKL_LONG strides[4];
    strides[0] = 0;
    strides[1] = _ny*(_nz/2+1);
    strides[2] = _nz/2+1;
    strides[3] = 1;
    DftiCreateDescriptor(&forwardHandle, DFTI_DOUBLE, DFTI_REAL, 3, length);
    DftiSetValue(forwardHandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    DftiSetValue(forwardHandle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DftiSetValue(forwardHandle, DFTI_OUTPUT_STRIDES, strides);
    status = DftiCommitDescriptor(forwardHandle);
    if (status) {
	throw std::runtime_error(DftiErrorMessage(status));
    }

    DftiCreateDescriptor(&inverseHandle, DFTI_DOUBLE, DFTI_REAL, 3, length);
    DftiSetValue(inverseHandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    DftiSetValue(inverseHandle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DftiSetValue(inverseHandle, DFTI_BACKWARD_SCALE, 1.0/(_nx*_ny*_nz));
    DftiSetValue(inverseHandle, DFTI_INPUT_STRIDES, strides);
    status = DftiCommitDescriptor(inverseHandle);
    if (status) {
	throw std::runtime_error(DftiErrorMessage(status));
    }

    DftiCreateDescriptor(&backwardHandle, DFTI_DOUBLE, DFTI_REAL, 3, length);
    DftiSetValue(backwardHandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    DftiSetValue(backwardHandle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DftiSetValue(backwardHandle, DFTI_INPUT_STRIDES, strides);
    status = DftiCommitDescriptor(backwardHandle);
    if (status) {
	throw std::runtime_error(DftiErrorMessage(status));
    }
}

FFT3Dmkl::~FFT3Dmkl()
{
    DftiFreeDescriptor(&forwardHandle);
    DftiFreeDescriptor(&inverseHandle);
    DftiFreeDescriptor(&backwardHandle);
}

void FFT3Dmkl::forward()
{
    DftiComputeForward(forwardHandle, FFT3D::pr3d->data(), FFT3D::pc3d->data());
}

void FFT3Dmkl::inverse()
{
    DftiComputeBackward(inverseHandle, FFT3D::pc3d->data(), FFT3D::pr3d->data());
}

void FFT3Dmkl::backward()
{
    DftiComputeBackward(backwardHandle, FFT3D::pc3d->data(), FFT3D::pr3d->data());
}
}
#endif
#endif
