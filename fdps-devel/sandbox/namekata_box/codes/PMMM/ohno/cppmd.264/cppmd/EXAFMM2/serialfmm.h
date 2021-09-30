#include <mpi.h>
#include "evaluator.h"

class SerialFMM : public Evaluator {
protected:
  int bodiesDispl[26];
  int bodiesCount[26];
  int multipoleDispl[26];
  int multipoleCount[26];
  int leafsDispl[26];
  int leafsCount[26];
  int IX[10][3];
  int gatherLevel;
  MPI_Comm MPI_COMM_LOCAL, MPI_COMM_GLOBAL;

private:
  void checkPartition(int *maxPartition) {
    int partitionSize = 1;
    for_3d partitionSize *= maxPartition[d];
    assert( MPISIZE == partitionSize );
    int checkLevel[3], partition[3];
    for_3d partition[d] = maxPartition[d];
    for( int d=0; d<3; d++ ) {
      int lev = 1;
      while( partition[d] != 1 ) {
        int ndiv = 2;
        if( (partition[d] % 3) == 0 ) ndiv = 3;
        partition[d] /= ndiv;
        lev++;
      }
      checkLevel[d] = lev - 1;
    }
    maxGlobLevel = FMMMAX(FMMMAX(checkLevel[0],checkLevel[1]),checkLevel[2]);
    for_3d numPartition[0][d] = 1;
    for_3d partition[d] = maxPartition[d];
    for( int lev=1; lev<=maxGlobLevel; lev++ ) {
      for( int d=0; d<3; d++ ) {
        int ndiv = 2;
        if( (partition[d] % 3) == 0 ) ndiv = 3;
        if( checkLevel[d] < maxGlobLevel && lev == 1 ) ndiv = 1;
        numPartition[lev][d] = ndiv * numPartition[lev-1][d];
        partition[d] /= ndiv;
      }
    }
  }

  void setSendCounts() {
    int multipoleType[3] = {8 * maxLevel,
                            4 * ((1 << (maxLevel + 1)) - 2),
                            2 * (((1 << (2 * maxLevel + 2)) - 1) / 3 - 1)};
    int leafsType[3] = {1, (1 << maxLevel), (1 << (2 * maxLevel))};
    int bodiesType[3];
    for_3d bodiesType[d] = leafsType[d] * maxBodies / numLeafs * 2;
    int i = 0;
    int ix[3];
    bodiesDispl[0] = multipoleDispl[0] = leafsDispl[0] = 0;
    for( ix[2]=-1; ix[2]<=1; ix[2]++ ) {
      for( ix[1]=-1; ix[1]<=1; ix[1]++ ) {
        for( ix[0]=-1; ix[0]<=1; ix[0]++ ) {
          if( ix[0] != 0 || ix[1] != 0 || ix[2] != 0 ) {
            int zeros = 0;
            for_3d zeros += ix[d] == 0;
            bodiesCount[i] = bodiesType[zeros];
            multipoleCount[i] = multipoleType[zeros];
            leafsCount[i] = leafsType[zeros];
            if( i > 0 ) {
              bodiesDispl[i] = bodiesDispl[i-1] + bodiesCount[i-1];
              multipoleDispl[i] = multipoleDispl[i-1] + multipoleCount[i-1];
              leafsDispl[i] = leafsDispl[i-1] + leafsCount[i-1];
            }
            i++;
          }
        }
      }
    }
    assert( numSendBodies >= bodiesDispl[25] + bodiesCount[25] );
    assert( numSendCells == multipoleDispl[25] + multipoleCount[25] );
    assert( numSendLeafs == leafsDispl[25] + leafsCount[25] );
  }

protected:
  inline void getIndex(int i, int *ix, real diameter) const {
    for_3d ix[d] = int((Jbodies[i][d] + R0 - X0[d]) / diameter);
  }

  void sort(int *key) const {
    int Imax = key[0];
    int Imin = key[0];
    for( int i=0; i<numBodies; i++ ) {
      Imax = FMMMAX(Imax,key[i]);
      Imin = FMMMIN(Imin,key[i]);
    }
    int numBucket = Imax - Imin + 1;
    int *bucket = new int [numBucket];
    for( int i=0; i<numBucket; i++ ) bucket[i] = 0;
    for( int i=0; i<numBodies; i++ ) bucket[key[i]-Imin]++;
    for( int i=1; i<numBucket; i++ ) bucket[i] += bucket[i-1];
    for( int i=numBodies-1; i>=0; --i ) {
      bucket[key[i]-Imin]--;
      int inew = bucket[key[i]-Imin];
      for_4d Ibodies[inew][d] = Jbodies[i][d];
    }
    for( int i=0; i<numBodies; i++ ) {
      for_4d Jbodies[i][d] = Ibodies[i][d];
      for_4d Ibodies[i][d] = 0;
    }
    delete[] bucket;
  }

public:
  void allocate(int N, int L, int I) {
    maxLevel = L;
    maxBodies = N;
    numImages = I;
    numCells = ((1 << 3 * (L + 1)) - 1) / 7;
    numLeafs = 1 << 3 * L;
    numSendCells = 64 * L + 48 * ((1 << (L + 1)) - 2) + 12 * (((1 << (2 * L + 2)) - 1) / 3 - 1);
    numSendLeafs = 8 + 12 * (1 << L) + 6 * (1 << (2 * L));
    numSendBodies = numSendLeafs * maxBodies / numLeafs * 2;
    Ibodies = new real [maxBodies][4]();
    Jbodies = new real [maxBodies+numSendBodies][4]();
    Multipole = new real [27*numCells][MTERM]();
    Local = new real [numCells][LTERM]();
    Leafs = new int [27*numLeafs][2]();
    globMultipole = new real [2*MPISIZE][MTERM]();
    globLocal = new real [10][LTERM]();
    sendJbodies = new float [numSendBodies][4]();
    recvJbodies = new float [numSendBodies][4]();
    sendMultipole = new float [numSendCells][MTERM]();
    recvMultipole = new float [numSendCells][MTERM]();
    sendLeafs = new int [numSendLeafs][2]();
    recvLeafs = new int [numSendLeafs][2]();
  }

  void deallocate() {
    delete[] Ibodies;
    delete[] Jbodies;
    delete[] Multipole;
    delete[] Local;
    delete[] Leafs;
    delete[] globMultipole;
    delete[] globLocal;
    delete[] sendJbodies;
    delete[] recvJbodies;
    delete[] sendMultipole;
    delete[] recvMultipole;
    delete[] sendLeafs;
    delete[] recvLeafs;
  }

  void partitioner(int *maxPartition, int level=0) {
    checkPartition(maxPartition);
    numGlobCells = 0;
    for( int lev=0; lev<=maxGlobLevel; lev++ ) {
      globLevelOffset[lev] = numGlobCells;
      numGlobCells += numPartition[lev][0] * numPartition[lev][1] * numPartition[lev][2];
    }
    getGlobIndex(IX[maxGlobLevel],MPIRANK,maxGlobLevel);
    for( int lev=maxGlobLevel; lev>0; lev-- ) {
      for_3d IX[lev-1][d] = IX[lev][d] * numPartition[lev-1][d] / numPartition[lev][d];
    }
    setSendCounts();
    gatherLevel = level;
    assert(gatherLevel <= maxGlobLevel);
    int ix[3], numChild[3];
    for_3d numChild[d] = numPartition[maxGlobLevel][d] / numPartition[gatherLevel][d];
    for_3d ix[d] = IX[maxGlobLevel][d] % numChild[d];
    int key = ix[0] + (ix[1] + ix[2] * numChild[1]) * numChild[0];
    int color = getGlobKey(IX[gatherLevel],gatherLevel);
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &MPI_COMM_LOCAL);
    MPI_Comm_split(MPI_COMM_WORLD, key, color, &MPI_COMM_GLOBAL);
  }

  void sortBodies() const {
    int *key = new int [numBodies];
    real diameter = 2 * R0 / (1 << maxLevel);
    int ix[3] = {0, 0, 0};
    for( int i=0; i<numBodies; i++ ) {
      getIndex(i,ix,diameter);
      key[i] = getKey(ix,maxLevel);
    }
    sort(key);
    delete[] key;
  }

  void buildTree() const {
    int rankOffset = 13 * numLeafs;
    for( int i=rankOffset; i<numLeafs+rankOffset; i++ ) {
      Leafs[i][0] = Leafs[i][1] = 0;
    }
    real diameter = 2 * R0 / (1 << maxLevel);
    int ix[3] = {0, 0, 0};
    getIndex(0,ix,diameter);
    int ileaf = getKey(ix,maxLevel,false) + rankOffset;
    Leafs[ileaf][0] = 0;
    for( int i=0; i<numBodies; i++ ) {
      getIndex(i,ix,diameter);
      int inew = getKey(ix,maxLevel,false) + rankOffset;
      if( ileaf != inew ) {
        Leafs[ileaf][1] = Leafs[inew][0] = i;
        ileaf = inew;
      }
    }
    Leafs[ileaf][1] = numBodies;
    for( int i=rankOffset; i<numLeafs+rankOffset; i++ ) {
      assert( Leafs[i][1] != Leafs[i][0] );
    }
  }
};
