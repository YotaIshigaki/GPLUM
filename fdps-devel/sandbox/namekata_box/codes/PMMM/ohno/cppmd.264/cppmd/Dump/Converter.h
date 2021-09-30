#ifndef CONVERTER_H
#define CONVERTER_H

#include "Dump.h"
#include <iosfwd>
#include <string>

template<class PA>
class Converter {
public:
  MPI_Comm comm;
  int rank;
  int num_rank;

  std::string filenamebase;
  int number_of_restorefiles;
  int number_of_localrestorefiles;
  std::vector< BinaryDump<PA> > restorefiles;
  std::vector<int> restorefileindexes;

  std::vector< PA > particlearrays;
  std::vector< std::vector<TypeRange> > typerangearrays;
  std::vector< std::vector<CovalentBondInfo::BondList> > bondlistarrays;
  std::vector< std::vector<CovalentBondInfo::BondList> > bondlistarray_idxs;
  std::vector< WaterList > waterlists;
  std::vector< ShakeList > shakelists;
  std::vector< long > timesteps;
  std::vector< double > kenergys;
  std::vector< double > penergys;
  std::vector< double > virials;
  std::vector< double > total_penergys;
  std::vector< double > eta_poss;
  std::vector< double > eta_vels;
  std::vector< double > eta_forces;
  std::vector< double > logv_poss;
  std::vector< double > logv_vels;
  std::vector< double > logv_forces;
  std::vector< SpaceVector<double> > boxsizes;
  std::vector< double > volumes;
  std::vector< double > tljcecs;
  std::vector< int > ris;
  std::vector< int > pis;
  std::vector< int > dcis;
  std::vector< int > dris;

  AtomID total_size_of_particlearray;
  AtomID local_number_of_particle;
  AtomID total_number_of_particle;
  AtomID maxid;
  AtomID local_maxid;

  PA particlearray;
  std::vector<TypeRange> typerangearray;
  std::vector<CovalentBondInfo::BondList> bondlistarray;
  WaterList waterlist;
  ShakeList shakelist;

  /* depden Integrator.cc */
  double kenergy;
  double penergy;
  double virial;
  double total_penergy;
  double eta_pos;
  double eta_vel;
  double eta_force;
  double logv_pos;
  double logv_vel;
  double logv_force;
  SpaceVector<double> boxsize;
  double volume;
  double tljcec;
  int ri;
  int pi;
  int dci;
  int dri;


  Converter(const std::string& fbase, int nrsts, MPI_Comm comm);
  ~Converter();

  int open_restorefiles();

  int read_restorefiles();

  void read_optionals();

  void read_optional();

  void read_rawimages();

  int read_typerangearrays();

  int merge_particle();
  int merge_particle_reindex();
  int merge_particle_by_atomid(PA& dpa, 
			       const PA& spa, const std::vector<TypeRange>& typerange);

  int merge_bondlistarray();
  int merge_waterlists();
  int merge_waterlists_atomid();
  int merge_shakelists_atomid();

  int dump_amber_rst(std::string& filename);
  int dump_amber_crd(std::string& filename);

  int gather_size();
  int gather_particle(PA& dpa, const PA& spa);
  int gather_bondlist(std::vector<CovalentBondInfo::BondList>& dbl,
		      const std::vector<CovalentBondInfo::BondList>& sbl);
  int gather_waterlist(WaterList& dwl, const WaterList& swl);
  int gather_shakelist(ShakeList& dsl, const ShakeList& ssl);

};

#endif 
