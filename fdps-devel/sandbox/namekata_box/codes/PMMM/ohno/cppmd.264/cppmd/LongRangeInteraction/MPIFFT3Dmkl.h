#ifndef MPIFFT3DMKL_H
#define MPIFFT3DMKL_H
#ifdef USE_MKL
#ifdef USE_MPI

#include "MPIFFT3D.h"
#include "mkl_cdft.h"

namespace PMEModule {

class MPIFFT3Dmkl: public FFT3D {
    public:
	MPIFFT3Dmkl(int nx, int ny, int nz,
	  int local_nx, int x_start,
	  int local_out_nx, int out_x_start,
	  int local_size);
	~MPIFFT3Dmkl();
	void forward();
	void inverse();
	void backward();
	static MPI_Comm getCommunicator() {
	    if (pComm) return pComm->communicator;
	    return MPI_COMM_WORLD;
	}
	static FFT3D::Ptr createFFT3D(int nx, int ny, int nz);
	StorageType getComplexStorageType() { return FULL_Z; }
	DivideType getRealDivideType() {
            return DIVIDE_X;
        }
    private:
	MPI_Comm fftWorld;
	Complex3D_Ptr pw3d;
	std::vector<std::complex<double> > workvec;
	DFTI_DESCRIPTOR_DM_HANDLE fftHandle;
	double factor;
};

FFT3D::Ptr MPIFFT3Dmkl::createFFT3D(int nx, int ny, int nz)
{
    DFTI_DESCRIPTOR_DM_HANDLE handle;
    MKL_LONG lengths[3] = {nx, ny, nz};
    long status = DftiCreateDescriptorDM(getCommunicator(),
    	&handle, DFTI_DOUBLE, DFTI_COMPLEX, 3, lengths);
    if (status) {
	throw std::runtime_error(DftiErrorMessage(status));
    }
    MKL_LONG local_size;
    status = DftiGetValueDM(handle, CDFT_LOCAL_SIZE, &local_size);
    if (status) {
	throw std::runtime_error(DftiErrorMessage(status));
    }
    MKL_LONG local_nx;
    status = DftiGetValueDM(handle, CDFT_LOCAL_NX, &local_nx);
    if (status) {
	throw std::runtime_error(DftiErrorMessage(status));
    }
    MKL_LONG x_start;
    status = DftiGetValueDM(handle, CDFT_LOCAL_X_START, &x_start);
    if (status) {
	throw std::runtime_error(DftiErrorMessage(status));
    }
    MKL_LONG out_nx;
    status = DftiGetValueDM(handle, CDFT_LOCAL_OUT_NX, &out_nx);
    if (status) {
	throw std::runtime_error(DftiErrorMessage(status));
    }
    MKL_LONG out_x_start;
    status = DftiGetValueDM(handle, CDFT_LOCAL_OUT_X_START,
	    &out_x_start);
    if (status) {
	throw std::runtime_error(DftiErrorMessage(status));
    }
    DftiFreeDescriptorDM(&handle);
    return FFT3D::Ptr(new MPIFFT3Dmkl(nx, ny, nz,
	    local_nx, x_start, out_nx, out_x_start, local_size));
}

MPIFFT3Dmkl::MPIFFT3Dmkl(int _nx, int _ny, int _nz,
			 int local_nx, int x_start,
			 int out_nx, int out_x_start,
			 int local_size):
    FFT3D::FFT3D(out_nx, _ny, _nz, 1, local_size, local_size),
    fftWorld(getCommunicator()),
    pw3d(createComplex3D(pc3d->data(), local_nx, _ny, _nz)),
    workvec(local_size), fftHandle(0), factor(1.0/(_nx*_ny*_nz))
    
{
    size[0] = _nx;
    size[1] = _ny;
    size[2] = _nz;
    Dims rshape = {{local_nx, ny, nz}};
    pr3d->reshape(rshape);
    Dims rbases = {{x_start, 0, 0}};
    pr3d->reindex(rbases);
    pw3d->reindex(rbases);
    Dims cbases = {{out_x_start, 0, 0}};
    pc3d->reindex(cbases);
    setStartEnd();

    MKL_LONG length[3];
    length[0] = _nx;
    length[1] = _ny;
    length[2] = _nz;
    long status;
    status = DftiCreateDescriptorDM(fftWorld, &fftHandle, DFTI_DOUBLE,
				    DFTI_COMPLEX, 3, length);
    if (status) {
	throw std::runtime_error(DftiErrorMessage(status));
    }
    status = DftiSetValueDM(fftHandle, CDFT_WORKSPACE, &workvec[0]);
    if (status) {
	throw std::runtime_error(DftiErrorMessage(status));
    }
    status = DftiCommitDescriptorDM(fftHandle);
    if (status) {
	throw std::runtime_error(DftiErrorMessage(status));
    }
}

MPIFFT3Dmkl::~MPIFFT3Dmkl()
{
    delete(pw3d);
    int finalized;
    MPI_Finalized(&finalized);
    if (finalized) return;
    DftiFreeDescriptorDM(&fftHandle);
}

void MPIFFT3Dmkl::forward()
{
    Real3D& r3d = *pr3d;
    Complex3D& w3d = *pw3d;
    for (int ix=realStart[0];ix<realEnd[0];++ix) {
	for (int iy=0;iy<ny;++iy) {
	    for (int iz=0;iz<nz;++iz) {
		w3d[ix][iy][iz] = r3d[ix][iy][iz];
	    }
	}
    }
    DftiComputeForwardDM(fftHandle, pc3d->data());
}

void MPIFFT3Dmkl::inverse()
{
    DftiComputeBackwardDM(fftHandle, pc3d->data());
    Real3D& r3d = *pr3d;
    Complex3D& w3d = *pw3d;
    for (int ix=realStart[0];ix<realEnd[0];++ix) {
	for (int iy=0;iy<ny;++iy) {
	    for (int iz=0;iz<nz;++iz) {
		r3d[ix][iy][iz] = factor*w3d[ix][iy][iz].real();
	    }
	}
    }
}

void MPIFFT3Dmkl::backward()
{
    DftiComputeBackwardDM(fftHandle, pc3d->data());
    Real3D& r3d = *pr3d;
    Complex3D& w3d = *pw3d;
    for (int ix=realStart[0];ix<realEnd[0];++ix) {
	for (int iy=0;iy<ny;++iy) {
	    for (int iz=0;iz<nz;++iz) {
		r3d[ix][iy][iz] = w3d[ix][iy][iz].real();
	    }
	}
    }
}
}
#endif
#endif
#endif
