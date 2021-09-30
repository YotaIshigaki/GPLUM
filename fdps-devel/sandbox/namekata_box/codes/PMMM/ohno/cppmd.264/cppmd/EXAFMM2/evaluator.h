#include "kernel.h"

class Evaluator : public Kernel {
public:
  bool printNow;

public:
  Evaluator() : printNow(true) {}

  double getTime() const {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return double(tv.tv_sec+tv.tv_usec*1e-6);
  }

  void upwardPass() const {
    double tic, toc;

    int rankOffset = 13 * numCells;
    for( int i=0; i<numCells; i++ ) {
      for_m Multipole[i+rankOffset][m] = 0;
      for_l Local[i][l] = 0;
    }

    tic = getTime();
    P2M();
    toc = getTime();
    if( printNow ) printf("P2M     : %lf : %f GFlops\n",toc-tic,148*numBodies/(toc-tic)*1e-9);

    tic = getTime();
    M2M();
    toc = getTime();
    if( printNow ) printf("M2M     : %lf : %f GFlops\n",toc-tic,955*numCells/(toc-tic)*1e-9);
  }

  void downwardPass() {
    double tic, toc;

    tic = getTime();
    M2L();
    toc = getTime();
    if( printNow ) printf("M2L     : %lf : %f GFlops\n",toc-tic,2417*numCells*189/(toc-tic)*1e-9);

    tic = getTime();
    L2L();
    toc = getTime();
    if( printNow ) printf("L2L     : %lf : %f GFlops\n",toc-tic,1902*numCells/(toc-tic)*1e-9);

    tic = getTime();
    L2P();
    toc = getTime();
    if( printNow ) printf("L2P     : %lf : %f GFlops\n",toc-tic,558*numBodies/(toc-tic)*1e-9);

#if DO_P2P
    tic = getTime();
    P2P();
    toc = getTime();
    if( printNow ) printf("P2P     : %lf : %f GFlops\n",toc-tic,22.*numBodies*numBodies*27/(toc-tic)/numLeafs*1e-9);
#endif
  }

  void direct() {
    int prange = numImages == 0 ? 0 : pow(3,numImages - 1);
    real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
    for( int i=0; i<100; i++ ) {
      real Po = 0, Fx = 0, Fy = 0, Fz = 0;
      int jx[3];
      for( jx[2]=-prange; jx[2]<=prange; jx[2]++ ) {
        for( jx[1]=-prange; jx[1]<=prange; jx[1]++ ) {
          for( jx[0]=-prange; jx[0]<=prange; jx[0]++ ) {
            for( int j=0; j<numBodies; j++ ) {
              real dist[3];
              for_3d dist[d] = Jbodies[i][d] - Jbodies[j][d] - jx[d] * 2 * RGlob[d];
              real R2 = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];
              real invR2 = 1.0 / R2;
              if( R2 == 0 ) invR2 = 0;
              real invR = Jbodies[j][3] * sqrt(invR2);
              real invR3 = invR2 * invR;
              Po += invR;
              Fx -= dist[0] * invR3;
              Fy -= dist[1] * invR3;
              Fz -= dist[2] * invR3;
            }
          }
        }
      }
      diff1 += (Ibodies[i][0] - Po) * (Ibodies[i][0] - Po);
      norm1 += Po * Po;
      diff2 += (Ibodies[i][1] - Fx) * (Ibodies[i][1] - Fx)
             + (Ibodies[i][2] - Fy) * (Ibodies[i][2] - Fy)
             + (Ibodies[i][3] - Fz) * (Ibodies[i][3] - Fz);
      norm2 += Fx * Fx + Fy * Fy + Fz * Fz;
    }
    printf("Err Pot : %lf\n",sqrt(diff1/norm1));
    printf("Err Forc: %lf\n",sqrt(diff2/norm2));
  }
};
