#ifndef HDF2BIN_H
#define HDF2BIN_H
#include <mpi.h>
#include "HDFSnapshotFile.h"
#include "Converter.h"

class HDF2Bin {
public:
  char **argv;
  int argc;
  MPI_Comm comm;
  int comm_rank;
  int comm_size;
  std::string args;

  std::string binname;
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

  HDF2Bin(int argc, char* argv[], MPI_Comm comm);

  void dumpBin(const std::string& binname,
	       const CombinedParticleArray& particlearray,
	       const std::vector<TypeRange>& typerangearray,
	       const std::vector<CovalentBondInfo::BondList>& bondlistarray,
	       const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
	       const WaterList& waterlist,
	       const ShakeList& shakelist,
	       const double[15],
	       const int[4],
	       const long t);
  void restoreHDF(const std::string& filename,
		  size_t step,
		  CombinedParticleArray& pa, std::vector<TypeRange>& tr,
		  std::vector<CovalentBondInfo::BondList>& bondlistarray,
		  std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
		  WaterList& waterlist,
		  ShakeList& shakelist,
		  long& t,
		  double od[15],
		  int oi[4]);
  void restoreSingleHDF();
  void restoreParallelDumpHDF();

  void NormalizeString(std::string &s);
  template<typename T>
    void Read(T &t);
  void ReadSnapshotType(HDFSnapshot::HDFSnapshotFile::SnapshotType &t);

  int ParseArguments(int arqgc, char* argv[]);

};
#endif
