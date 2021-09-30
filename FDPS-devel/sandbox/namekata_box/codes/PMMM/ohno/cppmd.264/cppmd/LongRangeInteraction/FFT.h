#ifndef FFT_H
#define FFT_H

#include "UseMPI.h"
#include "Array3D.h"
#ifdef CPPMD_ENABLE_LARGEMODEL
#include <Common.h>
#endif 
#include <complex>
#include <iostream>
#include <vector>


namespace FFTModule {

  struct Communicator {
#ifdef USE_MPI
    MPI_Comm communicator;
#endif
  };


class FFT3D {
 public:
  typedef PMEModule::Array3D<double> R3D;
  typedef PMEModule::Array3D<std::complex<double> > C3D;
  typedef R3D::ArrayRef Real3D;
  typedef R3D::Array Real3DArray;
  typedef R3D::Dims Dims;
  typedef C3D::ArrayRef Complex3D;

  typedef FFT3D* Ptr;
  typedef Real3D* Real3D_Ptr;
  typedef Complex3D* Complex3D_Ptr;


  typedef enum {
    FULL_Z,
    HALF_Z_EVEN,
    HALF_Z_ODD
  } StorageType;
  typedef enum {
    DIVIDE_NONE,
    DIVIDE_X
  } DivideType;  

  FFT3D(int _nx, int _ny, int _nz, int memtype=0,
	int rvecSize=0, int cvecSize=0);

  // FFT3D.h  98 - 101
  virtual ~FFT3D() {
    delete(pc3d);
    delete(pr3d);
  };

  static Ptr createFFT3D(int nx, int ny, int nz, bool single=false);

  // FFT3d.h 109 - 114
  void assign(double value);

  // FFT3d.h 107
  Real3D& getReal() { return *pr3d; } 

  // FFT3d.h 102
  virtual void forward()=0;

  // FFT3d.h 103
  virtual void backward()=0;

  // FFT3D.h 106
  void get_r3d_n(int &rnx, int &rny, int &rnz){rnx=r3d_nx;rny=r3d_ny;rnz=r3d_nz;}

  // FFT3d.h  115 - 122
  virtual StorageType getComplexStorageType() {
    if (nz % 2 == 0) {
      return HALF_Z_EVEN;
    }
    else {
      return HALF_Z_ODD;
    }
  }

  // FFT3D.h 123-125
  virtual DivideType getRealDivideType() {
    return DIVIDE_NONE;
  }

  // FFT3D.h 166-168
  int getRealStart(int i) {
    return realStart[i];
  }

  // FFT3D.h 169-171
  int getRealEnd(int i) {
    return realEnd[i];
  }

  // FFT3D.h 172-174
  int getComplexStart(int i) {
    return complexStart[i];
  }
  // FFT3D.h 175-177
  int getComplexEnd(int i) {
    return complexEnd[i];
  }

  // FFT3D.h 178-180
  int getSize(int i) {
    return size[i];
  }

  // FFT3D.h 181 - 183
  int getComplexAxis(int i) {
    return complexAxis[i];
  }

  // FFT3D.h 165
  static void setCommunicator(Communicator* _pComm) { pComm = _pComm; }

  // FFT3D.h 159-161
  static Real3D_Ptr createReal3D(int nx, int ny, int nz) {
    return Real3D_Ptr(R3D::createArray(nx, ny, nz));
  }

  // FFT3D.h 293 - 295
  static Real3D_Ptr createReal3D(double* data, int nx, int ny, int nz) {
    return Real3D_Ptr(R3D::createArray(data, nx, ny, nz));
  }
  // FFT3D.h 296 - 299
  static Complex3D_Ptr createComplex3D(std::complex<double>* data,
				       int nx, int ny, int nz) {
    return Complex3D_Ptr(C3D::createArray(data, nx, ny, nz));
  }
  // FFT3D.h 300 - 303
  static Complex3D_Ptr createComplex3D(double* data,
				       int nx, int ny, int nz) {
    return createComplex3D((std::complex<double>*)data, nx, ny, nz);
  }

  void setStartEnd();

  double convolution(const FFT3D::Real3D& BC);
  double convolution(const FFT3D::Real3D& BC, 
		     const FFT3D::Real3D& VBC, double &v);

  int r3d_nx;
  int r3d_ny;
  int r3d_nz;
  int nx;
  int ny;
  int nz;


  Real3D_Ptr pr3d;       //! Mesh assigned charge
  Complex3D_Ptr pc3d;    //! Complex reciprocal grid (results of forward), reuse pr3d for fftw2
  std::vector<double> rvec;  //! array for remap to pr3d
  std::vector<std::complex<double> > cvec;
  int realStart[3];
  int realEnd[3];
  int complexStart[3];
  int complexEnd[3];
  int size[3];
  int complexAxis[3];

  int timer_forward;
  int timer_backward;

 protected:
  static Communicator* pComm;
};



}
#endif

