/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#ifndef kernel_h
#define kernel_h
#include "sort.h"
#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

const int  P2 = P * P;                                          //!< P^2
const int  P4 = P2 * P2;                                        //!< P^4

//! Unified CPU/GPU kernel class
class KernelBase : public Sort {
protected:
  vect                 X0;                                      //!< Center of root cell
  real                 R0;                                      //!< Radius of root cell
  C_iter               Ci0;                                     //!< icells.begin()
  C_iter               Cj0;                                     //!< jcells.begin()

  int                  ATOMS;                                   //!< Number of atom types in Van der Waals
  std::vector<real>    RSCALE;                                  //!< Scaling parameter for Van der Waals
  std::vector<real>    GSCALE;                                  //!< Scaling parameter for Van der Waals

  CellMap              sourceBegin;                             //!< Define map for offset of source cells
  CellMap              sourceSize;                              //!< Define map for size of source cells
  CellMap              targetBegin;                             //!< Define map for offset of target cells

  real *factorial;                                              //!< Factorial
  real *prefactor;                                              //!< \f$ \sqrt{ \frac{(n - |m|)!}{(n + |m|)!} } \f$
  real *Anm;                                                    //!< \f$ (-1)^n / \sqrt{ \frac{(n + m)!}{(n - m)!} } \f$
  complex *Cnm;                                                 //!< M2L translation matrix \f$ C_{jn}^{km} \f$
public:
  real NP2P;                                                    //!< Number of P2P kernel calls
  real NM2P;                                                    //!< Number of M2P kernel calls
  real NM2L;                                                    //!< Number of M2L kernel calls

protected:
//! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
  void evalMultipole(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) const {
    const complex I(0.,1.);                                     // Imaginary unit
    real x = std::cos(alpha);                                   // x = cos(alpha)
    real y = std::sin(alpha);                                   // y = sin(alpha)
    real fact = 1;                                              // Initialize 2 * m + 1
    real pn = 1;                                                // Initialize Legendre polynomial Pn
    real rhom = 1;                                              // Initialize rho^m
    for( int m=0; m!=P; ++m ) {                                 // Loop over m in Ynm
      complex eim = std::exp(I * real(m * beta));               //  exp(i * m * beta)
      real p = pn;                                              //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
      int nmn = m * m;                                          //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * prefactor[npn] * eim;               //  rho^m * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
      real p1 = p;                                              //  Pnm-1
      p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;// theta derivative of r^n * Ynm
      rhom *= rho;                                              //  rho^m
      real rhon = rhom;                                         //  rho^n
      for( int n=m+1; n!=P; ++n ) {                             //  Loop over n in Ynm
        int npm = n * n + n + m;                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
        Ynm[npm] = rhon * p * prefactor[npm] * eim;             //   rho^n * Ynm
        Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
        real p2 = p1;                                           //   Pnm-2
        p1 = p;                                                 //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;// theta derivative
        rhon *= rho;                                            //   Update rho^n
      }                                                         //  End loop over n in Ynm
      pn = -pn * fact * y;                                      //  Pn
      fact += 2;                                                //  2 * m + 1
    }                                                           // End loop over m in Ynm
  }

//! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
  void evalLocal(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) const {
    const complex I(0.,1.);                                     // Imaginary unit
    real x = std::cos(alpha);                                   // x = cos(alpha)
    real y = std::sin(alpha);                                   // y = sin(alpha)
    real fact = 1;                                              // Initialize 2 * m + 1
    real pn = 1;                                                // Initialize Legendre polynomial Pn
    real rhom = 1.0 / rho;                                      // Initialize rho^(-m-1)
    for( int m=0; m!=2*P; ++m ) {                               // Loop over m in Ynm
      complex eim = std::exp(I * real(m * beta));               //  exp(i * m * beta)
      real p = pn;                                              //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
      int nmn = m * m;                                          //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * prefactor[npn] * eim;               //  rho^(-m-1) * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
      real p1 = p;                                              //  Pnm-1
      p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;// theta derivative of r^n * Ynm
      rhom /= rho;                                              //  rho^(-m-1)
      real rhon = rhom;                                         //  rho^(-n-1)
      for( int n=m+1; n!=2*P; ++n ) {                           //  Loop over n in Ynm
        int npm = n * n + n + m;                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
        Ynm[npm] = rhon * p * prefactor[npm] * eim;             //   rho^n * Ynm for m > 0
        Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
        real p2 = p1;                                           //   Pnm-2
        p1 = p;                                                 //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;// theta derivative
        rhon /= rho;                                            //   rho^(-n-1)
      }                                                         //  End loop over n in Ynm
      pn = -pn * fact * y;                                      //  Pn
      fact += 2;                                                //  2 * m + 1
    }                                                           // End loop over m in Ynm
  }

public:
//! Constructor
  KernelBase() : X0(0), R0(0), NP2P(0), NM2P(0), NM2L(0) {}
//! Destructor
  ~KernelBase() {}

//! Set center of root cell
  void setX0(vect x0) {X0 = x0;}
//! Set radius of root cell
  void setR0(real r0) {R0 = r0;}

//! Get center of root cell
  vect getX0() const {return X0;}
//! Get radius of root cell
  real getR0() const {return R0;}

//! Precalculate M2L translation matrix
  void preCalculation() {
    const complex I(0.,1.);                                     // Imaginary unit
    factorial = new real  [P];                                  // Factorial
    prefactor = new real  [4*P2];                               // sqrt( (n - |m|)! / (n + |m|)! )
    Anm       = new real  [4*P2];                               // (-1)^n / sqrt( (n + m)! / (n - m)! )
    Cnm       = new complex [P4];                               // M2L translation matrix Cjknm

    factorial[0] = 1;                                           // Initialize factorial
    for( int n=1; n!=P; ++n ) {                                 // Loop to P
      factorial[n] = factorial[n-1] * n;                        //  n!
    }                                                           // End loop to P

    for( int n=0; n!=2*P; ++n ) {                               // Loop over n in Anm
      for( int m=-n; m<=n; ++m ) {                              //  Loop over m in Anm
        int nm = n*n+n+m;                                       //   Index of Anm
        int nabsm = abs(m);                                     //   |m|
        real fnmm = EPS;                                        //   Initialize (n - m)!
        for( int i=1; i<=n-m; ++i ) fnmm *= i;                  //   (n - m)!
        real fnpm = EPS;                                        //   Initialize (n + m)!
        for( int i=1; i<=n+m; ++i ) fnpm *= i;                  //   (n + m)!
        real fnma = 1.0;                                        //   Initialize (n - |m|)!
        for( int i=1; i<=n-nabsm; ++i ) fnma *= i;              //   (n - |m|)!
        real fnpa = 1.0;                                        //   Initialize (n + |m|)!
        for( int i=1; i<=n+nabsm; ++i ) fnpa *= i;              //   (n + |m|)!
        prefactor[nm] = std::sqrt(fnma/fnpa);                   //   sqrt( (n - |m|)! / (n + |m|)! )
        Anm[nm] = ODDEVEN(n)/std::sqrt(fnmm*fnpm);              //   (-1)^n / sqrt( (n + m)! / (n - m)! )
      }                                                         //  End loop over m in Anm
    }                                                           // End loop over n in Anm

    for( int j=0, jk=0, jknm=0; j!=P; ++j ) {                   // Loop over j in Cjknm
      for( int k=-j; k<=j; ++k, ++jk ){                         //  Loop over k in Cjknm
        for( int n=0, nm=0; n!=P; ++n ) {                       //   Loop over n in Cjknm
          for( int m=-n; m<=n; ++m, ++nm, ++jknm ) {            //    Loop over m in Cjknm
            const int jnkm = (j+n)*(j+n)+j+n+m-k;               //     Index C_{j+n}^{m-k}
            Cnm[jknm] = std::pow(I,real(abs(k-m)-abs(k)-abs(m)))//     Cjknm
                      * real(ODDEVEN(j)*Anm[nm]*Anm[jk]/Anm[jnkm]) * EPS;
          }                                                     //    End loop over m in Cjknm
        }                                                       //   End loop over n in Cjknm
      }                                                         //  End loop over in k in Cjknm
    }                                                           // End loop over in j in Cjknm
  }

//! Free temporary allocations
  void postCalculation() {
    delete[] factorial;                                         // Free factorial
    delete[] prefactor;                                         // Free sqrt( (n - |m|)! / (n + |m|)! )
    delete[] Anm;                                               // Free (-1)^n / sqrt( (n + m)! / (n - m)! )
    delete[] Cnm;                                               // Free M2L translation matrix Cjknm
  }

//! Set paramters for Van der Waals
  void setVanDerWaals(int atoms, double *rscale, double *gscale) {
    assert(atoms <= 16);                                        // Change GPU constant memory alloc if needed
    THETA = .1;                                                 // Force opening angle to be small
    ATOMS = atoms;                                              // Set number of atom types
    RSCALE.resize(ATOMS*ATOMS);                                 // Resize rscale vector
    GSCALE.resize(ATOMS*ATOMS);                                 // Resize gscale vector
    for( int i=0; i!=ATOMS*ATOMS; ++i ) {                       // Loop over scale vector
      RSCALE[i] = rscale[i];                                    //  Set rscale vector
      GSCALE[i] = gscale[i];                                    //  Set gscale vector
    }                                                           // End loop over scale vector
  }
};

class Kernel : public KernelBase {
public:
  void initialize();                                            //!< Initialize kernels
#ifndef SPARC_SIMD
  void P2M(C_iter Ci) const;                                    //!< Evaluate P2M kernel on CPU
#else
  void P2M(Cell* Ci) const;                                    //!< Evaluate P2M kernel on CPU
#endif
  void M2M(C_iter Ci, C_iter Cj) const;                         //!< Evaluate M2M kernel on CPU
  void M2L(C_iter Ci, C_iter Cj) const;                         //!< Evaluate M2L kernel on CPU
  void M2P(C_iter Ci, C_iter Cj) const;                         //!< Evaluate M2P kernel on CPU
  void P2P(C_iter Ci, C_iter Cj) const;                         //!< Evaluate P2P kernel on CPU
  void L2L(C_iter Ci, C_iter Cj) const;                         //!< Evaluate L2L kernel on CPU
  void L2P(C_iter Ci) const;                                    //!< Evaluate L2P kernel on CPU
  void finalize();                                              //!< Finalize kernels
};

#include "cpuP2P.h"
#ifdef FMM_CARTESIAN
#include "cpuCartesianLaplace.h"
#else
#include "cpuSphericalLaplace.h"
#endif
#endif
