#include "UseConfig.h"
#ifdef USE_MPI
#ifdef USE_FFTW2
#include "MPIFFT3Dfftw2.h"
#endif
#ifdef USE_FFTE
#include "MPIFFT3Dffte.h"
#endif
#ifdef USE_MKL
#include "MPIFFT3Dmkl.h"
#endif
#endif /* USE_MPI */
#ifdef USE_FFTW2
#include "FFT3Dfftw2.h"
#endif
#ifdef USE_FFTW3
#include "FFT3Dfftw3.h"
#endif
#ifdef USE_FFTE
#include "FFT3Dffte.h"
#endif
#ifdef USE_MKL
#include "FFT3Dmkl.h"
#endif
#include "FFT3D.h"
#include<sstream>

using namespace PMEModule;

Communicator* FFT3D::pComm = 0;

void FFT3D::initialize()
{
    if (initialized) return;
#ifdef USE_MPI
#ifdef USE_FFTW2
    registerLibrary("fftw2", false, MPIFFT3Dfftw2::createFFT3D);
#endif
#ifdef USE_FFTE
    registerLibrary("ffte", false, MPIFFT3Dffte::createFFT3D,
	MPIFFT3Dffte::getAllowableLength);
#endif
#ifdef USE_MKL
    registerLibrary("mkl", false, MPIFFT3Dmkl::createFFT3D);
#endif
#endif /* USE_MPI */
#ifdef USE_FFTW2
    registerLibrary("fftw2", true, FFT3Dfftw2::createFFT3D);
#endif
#ifdef USE_FFTW3
    registerLibrary("fftw3", true, FFT3Dfftw3::createFFT3D);
#endif
#ifdef USE_FFTE
    registerLibrary("ffte", true, FFT3Dffte::createFFT3D,
	FFT3Dffte::getAllowableLength);
#endif
#ifdef USE_MKL
    registerLibrary("mkl", true, FFT3Dmkl::createFFT3D);
#endif
    initialized = true;
}

std::string FFT3D::getLibraryNames()
{
    initialize();
    std::ostringstream os;
    for (std::vector<LibraryNameFunc>::const_iterator it=libraries.begin();
	it != libraries.end();++it) {
	if (os.str().length() > 0) os << "/";
	os << it->libraryName;
    }
    return os.str();
}

FFT3D::Ptr FFT3D::createFFT3D(int nx, int ny, int nz, bool single)
{
    initialize();
    for (std::vector<LibraryNameFunc>::const_iterator it=libraries.begin();
	it != libraries.end();++it) {
        if (single && !it->single) continue;
	return it->func(nx, ny, nz);
    }
    throw std::runtime_error("no FFT libraries");
}

FFT3D::Ptr FFT3D::createFFT3D(const std::string& name,
                              int nx, int ny, int nz, bool single)
{
    initialize();
    for (std::vector<LibraryNameFunc>::const_iterator it=libraries.begin();
	it != libraries.end();++it) {
        if (single && !it->single) continue;
	if (name == it->libraryName) {
	    return it->func(nx, ny, nz);
	}
    }
    throw std::runtime_error("FFT library ("+name+") not support");
}

int FFT3D::getAllowableLength(int n, const std::string& name)
{
  initialize();
  for (std::vector<LibraryNameFunc>::const_iterator it=libraries.begin();
      it != libraries.end();++it) {
	if (name == "" || name == it->libraryName) {
	  if (it->funcA) {
	    return it->funcA(n);
	  }
	  else {
	    return n;
	  }
	}
  }
  if (libraries.size() == 0) {
    throw std::runtime_error("no FFT libraries");
  }
  return 0;
}

bool FFT3D::initialized = false;
std::vector<FFT3D::LibraryNameFunc> FFT3D::libraries;
