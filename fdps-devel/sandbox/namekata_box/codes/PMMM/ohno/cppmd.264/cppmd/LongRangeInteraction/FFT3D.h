#ifndef FFT3D_H
#define FFT3D_H

#include "Array3D.h"
#ifdef CPPMD_ENABLE_LARGEMODEL
#include <Common.h>
#endif 
#include <complex>
#include <iostream>
#include <vector>

namespace PMEModule {

struct Communicator;

class FFT3D {
 public:
  typedef Array3D<double> R3D;
  typedef Array3D<std::complex<double> > C3D;
  typedef R3D::ArrayRef Real3D;
  typedef R3D::Array Real3DArray;
  typedef R3D::Dims Dims;
  typedef C3D::ArrayRef Complex3D;

  typedef FFT3D* Ptr;
  typedef FFT3D::Real3D* Real3D_Ptr;
  typedef struct {
    std::string libraryName;
    FFT3D::Ptr (*func)(int nx, int ny, int nz);
    int (*funcA)(int n);
    bool single;
  } LibraryNameFunc;
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
        int rvecSize=0, int cvecSize=0):
      nx(_nx), ny(_ny), nz(_nz),
      pr3d(), pc3d(), rvec(), cvec(),
      realStart(), realEnd(),
      complexStart(), complexEnd()
  {
    size[0] = nx;
    size[1] = ny;
    size[2] = nz;
    r3d_nx = nx;
    r3d_ny = ny;
    r3d_nz = nz;
    complexAxis[0] = 0;
    complexAxis[1] = 1;
    complexAxis[2] = 2;
    switch (memtype) {
      case 0: /* in-place Real 3D FFT */
        if (rvecSize == 0) {
          rvecSize = nx*ny*(2*(nz/2+1));
        }
        rvec.assign(rvecSize, 0.0);
	r3d_nz = 2*(nz/2+1);
        pr3d = createReal3D(&rvec[0], nx, ny, 2*(nz/2+1));
        pc3d = createComplex3D(pr3d->data(), nx, ny, nz/2+1);
        break;
      case 1: /* out-of-place Complex 3D FFT */
        if (rvecSize == 0) {
          rvecSize = nx*ny*nz;
        }
        if (cvecSize == 0) {
          cvecSize = nx*ny*nz;
        }
        rvec.assign(rvecSize, 0.0);
        cvec.assign(cvecSize, 0.0);
        pr3d = createReal3D(&rvec[0], nx, ny, nz);
        pc3d = createComplex3D(&cvec[0], nx, ny, nz);
        break;
      case 2: /* out-of-place Real 3D FFT */
        if (rvecSize == 0) {
          rvecSize = nx*ny*nz;
        }
        if (cvecSize == 0) {
          cvecSize = nx*ny*(nz/2+1);
        }
        rvec.assign(rvecSize, 0.0);
        cvec.assign(cvecSize, 0.0);
        pr3d = createReal3D(&rvec[0], nx, ny, nz);
        pc3d = createComplex3D(&cvec[0], nx, ny, nz/2+1);
        break;
      default: 
        throw(std::runtime_error("FFT3D::FFT3D() invalid memtype"));
    }
    setStartEnd();
  }
  virtual ~FFT3D() {
    delete(pc3d);
    delete(pr3d);
  };
  virtual void forward()=0;
  virtual void backward()=0;
  virtual void inverse()=0;
  /* getReal() and getComplex() may be equal */
  void get_r3d_n(int &rnx, int &rny, int &rnz){rnx=r3d_nx;rny=r3d_ny;rnz=r3d_nz;}
  Real3D& getReal() { return *pr3d; } 
  Complex3D& getComplex() { return *pc3d; }
  void assign(double value) {
    Real3D& r3d = *pr3d;
    for (size_t i=0;i<r3d.num_elements();++i) {
      r3d.data()[i] = value;
    }
  }
  virtual StorageType getComplexStorageType() {
    if (nz % 2 == 0) {
      return HALF_Z_EVEN;
    }
    else {
      return HALF_Z_ODD;
    }
  }
  virtual DivideType getRealDivideType() {
    return DIVIDE_NONE;
  }
  void scaling(double factor) {
    Real3D& r3d = *pr3d;
    for (size_t i=0;i<r3d.num_elements();++i) {
      r3d.data()[i] *= factor;
    }
  }
  void printReal(std::ostream& fout=std::cout)
  {
    Real3D& r3d = *pr3d;
    for (int ix=realStart[0];ix<realEnd[0];++ix) {
      for (int iy=realStart[1];iy<realEnd[1];++iy) {
        for (int iz=realStart[2];iz<realEnd[2];++iz) {
          fout << "R(" << ix << "," << iy << "," << iz << ") " 
              << r3d[ix][iy][iz] << std::endl;
        }
      }
    }
  }
  void printComplex(std::ostream& fout=std::cout)
  {
    Complex3D& c3d = *pc3d;
    for (int ix=complexStart[0];ix<complexEnd[0];++ix) {
      for (int iy=complexStart[1];iy<complexEnd[1];++iy) {
        for (int iz=complexStart[2];iz<complexEnd[2];++iz) {
          fout << "C(" << ix << "," << iy << "," << iz << ") "
              << c3d[ix][iy][iz] << std::endl;
        }
      }
    }
  }
  static Ptr createFFT3D(int nx, int ny, int nz, bool single=false);
  static Ptr createFFT3D(const std::string& name,
                         int nx, int ny, int nz, bool single=false);
  static Real3D_Ptr createReal3D(int nx, int ny, int nz) {
    return Real3D_Ptr(R3D::createArray(nx, ny, nz));
  }
  static std::string getLibraryNames();
  static void initialize();
  static int getAllowableLength(int n, const std::string& name="");
  static void setCommunicator(Communicator* _pComm) { pComm = _pComm; }
  int getRealStart(int i) {
    return realStart[i];
  }
  int getRealEnd(int i) {
    return realEnd[i];
  }
  int getComplexStart(int i) {
    return complexStart[i];
  }
  int getComplexEnd(int i) {
    return complexEnd[i];
  }
  int getSize(int i) {
    return size[i];
  }
  int getComplexAxis(int i) {
    return complexAxis[i];
  }
  double convolution(const FFT3D::Real3D& BC)
  {
    FFT3D::Complex3D& FQ = *pc3d;  // FQ is Fourie transform of Q
    forward();
    int nzh = complexEnd[2]-1;
    double convnorm = 0.0;
    switch (getComplexStorageType()) {
      case FFT3D::FULL_Z:
        for (int mx=complexStart[0];mx<complexEnd[0];++mx) {
          for (int my=complexStart[1];my<complexEnd[1];++my) {
            for (int mz=complexStart[2];mz<complexEnd[2]; ++mz) {
              convnorm += BC[mx][my][mz]*norm(FQ[mx][my][mz]);  // Ulrich Essmann 1995 (4.7)
              FQ[mx][my][mz] *= BC[mx][my][mz];
            }
          }
        }
        convnorm *= 0.5;
        break;
      case FFT3D::HALF_Z_EVEN:
        for (int mx=complexStart[0];mx<complexEnd[0];++mx) {
          for (int my=complexStart[1];my<complexEnd[1];++my) {
            convnorm += 0.5*BC[mx][my][0]*norm(FQ[mx][my][0]);
            FQ[mx][my][0] *= BC[mx][my][0];
            convnorm += 0.5*BC[mx][my][nzh]*norm(FQ[mx][my][nzh]);
            FQ[mx][my][nzh] *= BC[mx][my][nzh];
            for (int mz=1;mz<nzh;++mz) {
              convnorm += BC[mx][my][mz]*norm(FQ[mx][my][mz]);
              FQ[mx][my][mz] *= BC[mx][my][mz];
            }
          }
        }
        break;
      case FFT3D::HALF_Z_ODD:
        for (int mx=complexStart[0];mx<complexEnd[0];++mx) {
          for (int my=complexStart[1];my<complexEnd[1];++my) {
            convnorm += 0.5*BC[mx][my][0]*norm(FQ[mx][my][0]);
            FQ[mx][my][0] *= BC[mx][my][0];
            for (int mz=1;mz<complexEnd[2];++mz) {
              convnorm += BC[mx][my][mz]*norm(FQ[mx][my][mz]);
              FQ[mx][my][mz] *= BC[mx][my][mz];
            }
          }
        }
        break;
    }
    backward();
    return convnorm;
  }
  double convolution(const FFT3D::Real3D& BC, 
		     const FFT3D::Real3D& VBC, double &v)
  {
    FFT3D::Complex3D& FQ = *pc3d;  // FQ is Fourie transform of Q
    forward();
    int nzh = complexEnd[2]-1;
    double convnorm = 0.0;
    double virial = 0.0;
    switch (getComplexStorageType()) {
      case FFT3D::FULL_Z:
        for (int mx=complexStart[0];mx<complexEnd[0];++mx) {
          for (int my=complexStart[1];my<complexEnd[1];++my) {
            for (int mz=complexStart[2];mz<complexEnd[2]; ++mz) {
              convnorm += BC[mx][my][mz]*norm(FQ[mx][my][mz]);  // Ulrich Essmann 1995 (4.7)
	      virial += VBC[mx][my][mz]*norm(FQ[mx][my][mz]);
              FQ[mx][my][mz] *= BC[mx][my][mz];
            }
          }
        }
        convnorm *= 0.5;
        virial *= 0.5;
        break;
      case FFT3D::HALF_Z_EVEN:
        for (int mx=complexStart[0];mx<complexEnd[0];++mx) {
          for (int my=complexStart[1];my<complexEnd[1];++my) {
            convnorm += 0.5*BC[mx][my][0]*norm(FQ[mx][my][0]);
	    virial += 0.5*VBC[mx][my][0]*norm(FQ[mx][my][0]);
            FQ[mx][my][0] *= BC[mx][my][0];
            convnorm += 0.5*BC[mx][my][nzh]*norm(FQ[mx][my][nzh]);
            virial += 0.5*VBC[mx][my][nzh]*norm(FQ[mx][my][nzh]);
            FQ[mx][my][nzh] *= BC[mx][my][nzh];
            for (int mz=1;mz<nzh;++mz) {
              convnorm += BC[mx][my][mz]*norm(FQ[mx][my][mz]);
              virial += VBC[mx][my][mz]*norm(FQ[mx][my][mz]);
              FQ[mx][my][mz] *= BC[mx][my][mz];
            }
          }
        }
        break;
      case FFT3D::HALF_Z_ODD:
        for (int mx=complexStart[0];mx<complexEnd[0];++mx) {
          for (int my=complexStart[1];my<complexEnd[1];++my) {
            convnorm += 0.5*BC[mx][my][0]*norm(FQ[mx][my][0]);
            virial += 0.5*VBC[mx][my][0]*norm(FQ[mx][my][0]);
            FQ[mx][my][0] *= BC[mx][my][0];
            for (int mz=1;mz<complexEnd[2];++mz) {
              convnorm += BC[mx][my][mz]*norm(FQ[mx][my][mz]);
              virial += VBC[mx][my][mz]*norm(FQ[mx][my][mz]);
              FQ[mx][my][mz] *= BC[mx][my][mz];
            }
          }
        }
        break;
    }
    backward();
    v = -(convnorm + 2.0*virial);
    return convnorm;
  }
 protected:
  typedef Complex3D* Complex3D_Ptr;

  static Real3D_Ptr createReal3D(double* data, int nx, int ny, int nz) {
    return Real3D_Ptr(R3D::createArray(data, nx, ny, nz));
  }
  static Complex3D_Ptr createComplex3D(std::complex<double>* data,
                                       int nx, int ny, int nz) {
    return Complex3D_Ptr(C3D::createArray(data, nx, ny, nz));
  }
  static Complex3D_Ptr createComplex3D(double* data,
                                       int nx, int ny, int nz) {
    return createComplex3D((std::complex<double>*)data, nx, ny, nz);
  }
  void setStartEnd() {
#ifdef CPPMD_ENABLE_LARGEMODEL
    if(DebugLog::verbose>1){
    std::cout << "setStartEnd()" << std::endl;
#endif
    for (int i=0;i<3;++i) {
      realStart[i] = pr3d->index_bases()[i];
      realEnd[i] = realStart[i]+pr3d->shape()[i];
      complexStart[i] = pc3d->index_bases()[i];
      complexEnd[i] = complexStart[i]+pc3d->shape()[i];
    }
    if (getComplexStorageType() != FULL_Z) {
      realEnd[2] = std::min(realEnd[2], nz);
    }
  }
  void setNullStartEnd() {
    for (int i=0;i<3;++i) {
      realStart[i] = 0;
      realEnd[i] = 0;
      complexStart[i] = 0;
      complexEnd[i] = 0;
    }
  }

  static Communicator* pComm;

  int r3d_nx;
  int r3d_ny;
  int r3d_nz;
  int nx;
  int ny;
  int nz;
  Real3D_Ptr pr3d;
  Complex3D_Ptr pc3d;
  std::vector<double> rvec;
  std::vector<std::complex<double> > cvec;
  int realStart[3];
  int realEnd[3];
  int complexStart[3];
  int complexEnd[3];
  int size[3];
  int complexAxis[3];
 private:
  static void registerLibrary(const std::string& name, bool single,
                              FFT3D::Ptr (*func)(int nx, int ny, int nz),
                              int (*funcA)(int n)=0) {
    LibraryNameFunc lib;
    lib.libraryName = name;
    lib.func = func;
    lib.funcA = funcA;
    lib.single = single;
    libraries.push_back(lib);
  }
  static bool initialized;
  static std::vector<LibraryNameFunc> libraries;
};
}
#endif
