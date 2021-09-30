#ifndef BIN2HDF_H
#define BIN2HDF_H
#include <mpi.h>
#include "HDFSnapshotFile.h"
#include "Converter.h"

class Bin2HDF {
public:
  char **argv;
  int argc;
  MPI_Comm comm;
  int comm_rank;
  int comm_size;
  std::string args;

  std::string filenamebase;
  int number_of_restorefiles;
  std::string hdfname;
  int with_waterlist;
  int shrink;
  HDFSnapshot::HDFSnapshotFile::SnapshotType snapshottype;

  Bin2HDF(int argc, char* argv[], MPI_Comm comm);

  void outputSingleHDF(Converter<CombinedParticleArray>& converter);
  void dumpParallelHDF(Converter<CombinedParticleArray>& converter);

  void NormalizeString(std::string &s);
  template<typename T>
    void Read(T &t);
  void ReadSnapshotType(HDFSnapshot::HDFSnapshotFile::SnapshotType &t);

  int ParseArguments(int arqgc, char* argv[]);

};
#endif
