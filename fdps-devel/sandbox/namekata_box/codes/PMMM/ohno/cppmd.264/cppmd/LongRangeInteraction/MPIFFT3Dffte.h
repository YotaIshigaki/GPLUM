#ifndef MPIFFT3DFFTE_H
#define MPIFFT3DFFTE_H
#ifdef USE_FFTE
#ifdef USE_MPI

#include "MPIFFT3D.h"

extern "C" {
    void pzfft3d_(std::complex<double>* a, std::complex<double>* b,
	int *nx, int *ny, int *nz, int *icomm, int *npu, int *iopt);
}

namespace PMEModule {

class MPIFFT3Dffte: public FFT3D {
    public:
	MPIFFT3Dffte(int nx, int ny, int nz, int npu, int myrank);
	MPIFFT3Dffte(int nx, int ny, int nz, int npu, int myrank,
		     MPI_Comm comm);
	~MPIFFT3Dffte() {
          delete(pw3d);
          int finalized;
          MPI_Finalized(&finalized);
          if (finalized) return;
          if (hasOwnComm) MPI_Comm_free(&fftWorld);
        }
	void forward();
	void inverse();
	void backward();
	static FFT3D::Ptr createFFT3D(int nx, int ny, int nz) {
	    int npu = getNumberOfProcessor();
	    int myrank = getRank();
	    if (nx % npu == 0 && nz % npu == 0) {
		return FFT3D::Ptr(new MPIFFT3Dffte(nx, ny, nz, npu, myrank));
	    }
	    else {
		MPI_Comm myOwnComm = getOwnCommunicator(nx, ny, nz, npu);
		if (myrank < npu) {
		  myrank = getRank(myOwnComm);
		}
		return FFT3D::Ptr(new MPIFFT3Dffte(nx, ny, nz, npu, myrank,
						   myOwnComm));
	    }
	}
	static MPI_Comm getCommunicator() {
	    if (pComm) return pComm->communicator;
	    return MPI_COMM_WORLD;
	}
	static int getAllowableLength(int n);
	StorageType getComplexStorageType() { return FULL_Z; }
	DivideType getRealDivideType() {
            return DIVIDE_X;
        }
    private:
	static int getNumberOfProcessor() {
	    int npu;
	    MPI_Comm comm = (MPI_Comm)getCommunicator();
	    MPI_Comm_size(comm, &npu);
	    return npu;
	}
	static int getRank() {
	    int rank;
	    MPI_Comm comm = (MPI_Comm)getCommunicator();
	    MPI_Comm_rank(comm, &rank);
	    return rank;
	}
	static int getRank(MPI_Comm comm) {
	    int rank;
	    MPI_Comm_rank(comm, &rank);
	    return rank;
	}
	static MPI_Comm getOwnCommunicator(int nx, int ny, int nz, int &npu);
	void initialize();

	int nax, nay, naz, npu;
	MPI_Comm fftWorld;
	int icomm;
	int myrank;
	std::vector<std::complex<double> > wvec;
	Complex3D_Ptr pw3d;
	bool hasOwnComm;
};

MPIFFT3Dffte::MPIFFT3Dffte(int _nx, int _ny, int _nz, int _npu, int _myrank):
    FFT3D::FFT3D(_nx/_npu, _ny, _nz, 1),
    nax(_nx), nay(_ny), naz(_nz), npu(_npu),
    fftWorld(getCommunicator()), icomm(MPI_Comm_c2f(fftWorld)),
    myrank(_myrank), wvec(_nx*_ny*_nz/npu, 0.0),
    pw3d(createComplex3D(&wvec[0], _nx/npu, _ny, _nz)),
    hasOwnComm(false)
{
    size[0] = nax;
    size[1] = nay;
    size[2] = naz;
    initialize();
}

MPIFFT3Dffte::MPIFFT3Dffte(int _nx, int _ny, int _nz, int _npu, int _myrank,
			   MPI_Comm comm):
    FFT3D::FFT3D(_nx/_npu, _ny, _nz, 1),
    nax(_nx), nay(_ny), naz(_nz), npu(_npu),
    fftWorld(comm), icomm(MPI_Comm_c2f(fftWorld)),
    myrank(_myrank), wvec(_nx*_ny*_nz/npu, 0.0),
    pw3d(createComplex3D(&wvec[0], _nx/npu, _ny, _nz)),
    hasOwnComm(true)
{
    int nodeNum;
    MPI_Comm_size(fftWorld, &nodeNum);
    size[0] = nax;
    size[1] = nay;
    size[2] = naz;
    if (myrank >= npu) {
	setNullStartEnd();
    }
    else {
	initialize();
    }
}

void MPIFFT3Dffte::initialize()
{
    if (nax % npu != 0) {
	throw(std::runtime_error("nx must be a multiple of npu"));
    }
    if (naz % npu != 0) {
	throw(std::runtime_error("nz must be a multiple of npu"));
    }
    Dims rbases = {{(nax/npu)*myrank, 0, 0}};
    pr3d->reindex(rbases);
    pw3d->reindex(rbases);
    Dims cshape = {{nax, nay, naz/npu}};
    pc3d->reshape(cshape);
    Dims cbases = {{0, 0, (naz/npu)*myrank}};
    pc3d->reindex(cbases);
    setStartEnd();

    int iopt=0;
    Complex3D& w3d = *pw3d;
    Complex3D& c3d = *pc3d;
    pzfft3d_(w3d.data(), c3d.data(), &naz, &nay, &nax, &icomm, &npu, &iopt);
}

void MPIFFT3Dffte::forward()
{
    if (myrank >= npu) return;
    int iopt=-2;
    Complex3D& w3d = *pw3d;
    Complex3D& c3d = *pc3d;
    Real3D& r3d = *pr3d;
    int sx = r3d.index_bases()[0];
    int ex = sx+r3d.shape()[0];
    for (int ix=sx;ix<ex;++ix) {
	for (int iy=0;iy<nay;++iy) {
	    for (int iz=0;iz<naz;++iz) {
		w3d[ix][iy][iz] = r3d[ix][iy][iz];
	    }
	}
    }
    pzfft3d_(w3d.data(), c3d.data(), &naz, &nay, &nax, &icomm, &npu, &iopt);
}

void MPIFFT3Dffte::inverse()
{
    if (myrank >= npu) return;
    int iopt=2;
    Complex3D& c3d = *pc3d;
    Complex3D& w3d = *pw3d;
    Real3D& r3d = *pr3d;
    pzfft3d_(c3d.data(), w3d.data(), &naz, &nay, &nax, &icomm, &npu, &iopt);
    int sx = r3d.index_bases()[0];
    int ex = sx+r3d.shape()[0];
    for (int ix=sx;ix<ex;++ix) {
	for (int iy=0;iy<nay;++iy) {
	    for (int iz=0;iz<naz;++iz) {
		r3d[ix][iy][iz] = w3d[ix][iy][iz].real();
	    }
	}
    }
}

void MPIFFT3Dffte::backward()
{
    if (myrank >= npu) return;
    int iopt=2;
    Complex3D& c3d = *pc3d;
    Complex3D& w3d = *pw3d;
    Real3D& r3d = *pr3d;
    pzfft3d_(c3d.data(), w3d.data(), &naz, &nay, &nax, &icomm, &npu, &iopt);
    int sx = r3d.index_bases()[0];
    int ex = sx+r3d.shape()[0];
    double factor = nax*nay*naz;
    for (int ix=sx;ix<ex;++ix) {
	for (int iy=0;iy<nay;++iy) {
	    for (int iz=0;iz<naz;++iz) {
		r3d[ix][iy][iz] = factor*w3d[ix][iy][iz].real();
	    }
	}
    }
}

int MPIFFT3Dffte::getAllowableLength(int n)
{
    int npu = getNumberOfProcessor();
    {
	int m = npu;
	while (m % 2 == 0) m /= 2;
	while (m % 3 == 0) m /= 3;
	while (m % 5 == 0) m /= 5;
	if (m > 1) {
	    throw std::runtime_error("npu must be (2^l)*(3^m)*(5^n)");
	}
    }
    for (;;++n) {
	if (n % npu > 0) continue;
	int m = n;
	while (m % 2 == 0) m /= 2;
	while (m % 3 == 0) m /= 3;
	while (m % 5 == 0) m /= 5;
	if (m == 1) break;
    }
    return n;
}

MPI_Comm MPIFFT3Dffte::getOwnCommunicator(int nx, int ny, int nz, int& npu)
{
    for (int enpu=npu;enpu > 0;--enpu) {
	if (nx % enpu == 0 && nz % enpu == 0) {
	    npu = enpu;
	    int rank;
	    MPI_Comm comm = getCommunicator();
	    MPI_Comm_rank(comm, &rank);
	    int color = (rank < npu)?0:1;
	    MPI_Comm newcomm;
	    MPI_Comm_split(comm, color, rank, &newcomm);
	    return newcomm;
	}
    }
    throw(std::runtime_error("number of processor error"));
}
}
#endif
#endif
#endif
