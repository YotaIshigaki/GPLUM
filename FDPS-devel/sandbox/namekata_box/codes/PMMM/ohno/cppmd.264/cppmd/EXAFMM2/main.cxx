#if Serial
#include "serialfmm.h"
#else
#include "parallelfmm.h"
#endif

const int N = 1000;
int main() {
  double tic, toc;
#if Serial
  SerialFMM FMM;
#else
  ParallelFMM FMM;
#endif
  srand48(FMM.MPIRANK);
  int numBodies = N + int(drand48() * N);
  int maxBodies = 2 * N;
  FMM.allocate(maxBodies,2,1);
  FMM.numBodies = numBodies;
  if( FMM.printNow ) {
    printf("N       : %d\n",FMM.numBodies);
    printf("Levels  : %d\n",FMM.maxLevel);
    printf("Images  : %d\n",FMM.numImages);
    printf("------------------\n");
  }

  int mpisize = FMM.MPISIZE;
  int maxPartition[3] = {1, 1, 1};
  int dim = 0;
  while( mpisize != 1 ) {
    int ndiv = 2;
    if( (mpisize % 3) == 0 ) ndiv = 3;
    maxPartition[dim] *= ndiv;
    mpisize /= ndiv;
    dim = (dim + 1) % 3;
  }

  tic = FMM.getTime();
  FMM.partitioner(maxPartition,1);
  toc = FMM.getTime();
  if( FMM.printNow ) printf("Part    : %lf\n",toc-tic);

  for( int it=0; it<1; it++ ) {
    int ix[3] = {0, 0, 0};
    srand48(FMM.MPIRANK+it*FMM.MPISIZE);
    FMM.numBodies = N + int(drand48() * N);
    FMM.R0 = .5;
    for_3d FMM.RGlob[d] = FMM.R0 * FMM.numPartition[FMM.maxGlobLevel][d];
    FMM.getGlobIndex(ix,FMM.MPIRANK,FMM.maxGlobLevel);
    for_3d FMM.X0[d] = 2 * FMM.R0 * (ix[d] + .5);
    srand48(FMM.MPIRANK);
    real average = 0;
    for( int i=0; i<FMM.numBodies; i++ ) {
      FMM.Jbodies[i][0] = drand48() + 2 * FMM.R0 * ix[0];
      FMM.Jbodies[i][1] = drand48() + 2 * FMM.R0 * ix[1];
      FMM.Jbodies[i][2] = drand48() + 2 * FMM.R0 * ix[2];
      FMM.Jbodies[i][3] = (drand48() - .5) / FMM.numBodies;
      average += FMM.Jbodies[i][3];
    }
    average /= FMM.numBodies;
    for( int i=0; i<FMM.numBodies; i++ ) {
      FMM.Jbodies[i][3] -= average;
    }
  
    tic = FMM.getTime();
    FMM.sortBodies();
    toc = FMM.getTime();
    if( FMM.printNow ) printf("Sort    : %lf\n",toc-tic);
  
    tic = FMM.getTime();
    FMM.buildTree();
    toc = FMM.getTime();
    if( FMM.printNow ) printf("Tree    : %lf\n",toc-tic);
  
    FMM.upwardPass();
  
#if Serial
#else
    tic = FMM.getTime();
    FMM.P2PSend();
    toc = FMM.getTime();
    if( FMM.printNow ) printf("P2P Send: %lf\n",toc-tic);

    tic = FMM.getTime();
    FMM.P2PRecv();
    toc = FMM.getTime();
    if( FMM.printNow ) printf("P2P Recv: %lf\n",toc-tic);

    tic = FMM.getTime();
    FMM.M2LSend();
    toc = FMM.getTime();
    if( FMM.printNow ) printf("M2L Send: %lf\n",toc-tic);

    tic = FMM.getTime();
    FMM.M2LRecv();
    toc = FMM.getTime();
    if( FMM.printNow ) printf("M2L Recv: %lf\n",toc-tic);

    tic = FMM.getTime();
    FMM.rootGather();
    toc = FMM.getTime();
    if( FMM.printNow ) printf("Gather  : %lf\n",toc-tic);
  
    FMM.globM2M();
  
    FMM.globM2L();
#endif
  
    tic = FMM.getTime();
    FMM.periodicM2L();
    toc = FMM.getTime();
    if( FMM.printNow ) printf("M2L Peri: %lf\n",toc-tic);
  
#if Serial
#else
    tic = FMM.getTime();
    FMM.globL2L();
    toc = FMM.getTime();
    if( FMM.printNow ) printf("L2L Glob: %lf\n",toc-tic);
#endif
  
    FMM.downwardPass();
  
    tic = FMM.getTime();
#if Serial
    FMM.direct();
#else
    FMM.globDirect();
#endif
    toc = FMM.getTime();
    if( FMM.printNow ) printf("Direct  : %lf\n",toc-tic);
  }
  FMM.deallocate();

}
