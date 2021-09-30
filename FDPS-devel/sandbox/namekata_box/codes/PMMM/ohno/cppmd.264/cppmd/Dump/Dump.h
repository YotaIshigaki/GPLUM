#ifndef DUMP_H
#define DUMP_H

#include <cstdio>
#include <vector>
#include "ParticleInfo.h"
#include "AmberRestrt.h"

#ifdef TTHA_VITRO
#define TTHA_NCHAIN 4
#else
#ifdef TTHA_VIVO
#define TTHA_NCHAIN 8
#endif
#endif


class Dump {
 public:
  amberfile::AmberRestrt amber_crd;
  SpaceVector<double> boxsize;
  int num_particle;
  int num_node;
  int node_id;
#ifdef DEBUG_UNITOUTPUT
  int num_copy;
#endif  // DEBUG_UNITOUTPUT
  int is_restart;
  FILE *fp_;
#ifdef TTHA_NCHAIN
  FILE *fp_ttha[TTHA_NCHAIN];
#endif


  struct mpi_buffer {
    char *begin;
    int *np;
    int *atomid;
    Position *position;
    Atomtype *atomtype;
    Velocity *velocity;

    mpi_buffer()
    {
      begin = NULL;
      np = NULL;
      atomid = NULL;
      position = NULL;
      atomtype = NULL;
      velocity = NULL;
    }
  };
  
  mpi_buffer buffer;
  int buffer_size;

  int exclude_type;
  std::vector<Atomtype> exclude_types;
  std::vector<Atomtype> type_list;

#ifdef DEBUG_UNITOUTPUT
  Dump(const SpaceVector<double>& bs,
       int np, int nn, int me, int exclude, int num_copy = 1, int restart=0);

  Dump(const SpaceVector<double>& bs,
       int np, int nn, int me, int exclude,
       const std::string& filename, int num_copy = 1, int restart=0);
#else  // DEBUG_UNITOUTPUT
  Dump(const SpaceVector<double>& bs,
       int np, int nn, int me, int exclude, int restart=0);

  Dump(const SpaceVector<double>& bs,
       int np, int nn, int me, int exclude,
       const std::string& filename, int restart=0);
#endif  // DEBUG_UNITOUTPUT

  ~Dump();

  void init(const std::string& filename);
  template<typename PA>
  void setAmberCRD(const PA& pa,
                   const std::vector<TypeRange>& typerange);
  void DumpAmberCRD(const ParticleArray& pa,
                    const std::vector<TypeRange>& typerange,
                    const std::string& filename);
  template<typename PA>
  void DumpAmberCRD(const PA& pa,
                    const std::vector<TypeRange>& typerange);

  template<typename PA>
  void GatherDumpAmberCRD(const PA& pa,
                          const std::vector<TypeRange>& typerange);

  void setBoxSize(const SpaceVector<double>& bs);

  void DumpCoordinatesToFile();
  void setExcludeTypes();
};


#include "CovalentBondInfo.h"

enum BinaryDumpMode{
  BDUMP,
  BRESTORE
};

template<class PA>
class BinaryDump {
 public:
  explicit BinaryDump(const int _mode=BDUMP);
  explicit BinaryDump(const std::string& filename, const int _mode=BDUMP);
  ~BinaryDump();
  void dumpmode();
  void restoremode();
  int init(const std::string& filename);
  int close();

  void dump_binary_basic(const PA& particlearray,
			 const std::vector<TypeRange>& typerangearray,
			 const std::vector<CovalentBondInfo::BondList>& bondlistarray,
			 const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
			 const WaterList& waterlist,
			 const ShakeList& shakelist,
			 const long t);
  void dump_binary_optional_double(const double *od, const int num);
  void dump_binary_optional_int(const int *oi, const int num);

  void restore_binary_typerangearray(std::vector<TypeRange>& typerangearray);
  void restore_binary_basic(PA& particlearray,
			    std::vector<TypeRange>& typerangearray,
			    std::vector<CovalentBondInfo::BondList>& bondlistarray,
			    std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
			    WaterList& waterlist,
			    ShakeList& shakelist,
			    long& t);

  void merge_binary_basic_by_atomid(PA& particlearray,
				    CovalentBondInfo::BondList& bondlist,
				    WaterList& waterlist,
				    ShakeList& shakelist,
				    long& t);

  void restore_binary_optional_double(double *od, const int num);
  void restore_binary_optional_int(int *oi, const int num);

  int exact_num_particle();

 private:
  void writeint(const int i);
  void writelong(const long l);
  void writedouble(const double d);
  int readint();
  long readlong();
  double readdouble();
  
  FILE *fp;
  int mode;
  bool native_byteorder;
  long count;
  int num_particle;
};

#endif
