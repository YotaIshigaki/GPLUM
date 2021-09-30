#ifndef HDF2RST_H
#define HDF2RST_H
#include <mpi.h>
#include "HDFSnapshotFile.h"
#include "Converter.h"

class HDF2RST {
public:
  char **argv;
  int argc;
  MPI_Comm comm;
  int comm_rank;
  int comm_size;
  std::string args;

  std::string rstname;
  //  int number_of_restorefiles;
  std::string hdfname;
  //  int with_waterlist;
  //  int shrink;
  size_t step;
  HDFSnapshot::HDFSnapshotFile::SnapshotType snapshottype;

  int number_of_restorefiles;
  int number_of_localrestorefiles;
  AtomID total_number_of_particle;
  std::vector<int> restorefileindexes;
  AtomID maxid;
  AtomID local_maxid;

  HDF2RST(int argc, char* argv[], MPI_Comm comm);

  void restoreSingleHDF();
  void restoreParallelDumpHDF();

  void NormalizeString(std::string &s);
  template<typename T>
    void Read(T &t);
  void ReadSnapshotType(HDFSnapshot::HDFSnapshotFile::SnapshotType &t);

  int ParseArguments(int arqgc, char* argv[]);

};
#endif
