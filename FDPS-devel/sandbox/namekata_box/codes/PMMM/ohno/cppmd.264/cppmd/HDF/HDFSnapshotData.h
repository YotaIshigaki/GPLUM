#ifndef HDFSNAPSHOTDATA_H
#define HDFSNAPSHOTDATA_H

#include <cstdio>
#include <vector>
#include "CovalentBondInfo.h"
#include "ParticleInfo.h"
#include "HDFSnapshotFile.h"

enum HDFDumpMode{
  HDFDUMP,
  HDFRESTORE
};

class HDFDump {
 public:
  explicit HDFDump(const int _mode=HDFDUMP);
  explicit HDFDump(const std::string& filename, const int _mode=HDFDUMP);
  ~HDFDump();
  void dumpmode();
  void restoremode();
  int init(const std::string& filename);
  int close();

  void setParameter(double deltat, int step,
		    int reduction_counter, int print_counter,
		    int crd_counter, int rst_counter,
		    AtomID size_of_particle);
  void setInfo(int rank, int num_rank, string* version, AtomID size_of_particle, HDFSnapshot::HDFSnapshotFile::SnapshotType snapshot_type);
  void setIntegration(double LJEnergyCorrection,
		      double KineticEnergy,
		      double PotentialEnergy,
		      double TotalEnergy,
		      double Virial,
		      double EtaPosition,
		      double EtaVelocity,
		      double EtaForce,
		      double LogVPosition,
		      double LogVVelocity,
		      double LogVForce,
		      double Volume );
  void setBoxDatatype(SpaceVector<double> boxsize, int type, SpaceVector<double> angle);

  HDFSnapshot::HDFSnapshotFile::SnapshotType getSnapshotType() {
    return snapshot->getSnapshotType();
  }

  void getInfo(int& rank, int& num_rank, string& version, AtomID& size_of_particle,
	       HDFSnapshot::HDFSnapshotFile::SnapshotType& snapshot_type);
  void getParameter(double& deltat, int& step,
		    int& reduction_counter, int& print_counter,
		    int& crd_counter, int& rst_counter,
		    AtomID& size_of_particle);
  void getIntegration(double& LJEnergyCorrection,
		      double& KineticEnergy,
		      double& PotentialEnergy,
		      double& TotalEnergy,
		      double& Virial,
		      double& EtaPosition,
		      double& EtaVelocity,
		      double& EtaForce,
		      double& LogVPosition,
		      double& LogVVelocity,
		      double& LogVForce,
		      double& Volume ); 
  void getBoxDatatype(SpaceVector<double>& boxsize, int& type, SpaceVector<double>& angle);
  template<class PA>
    void getParticle(PA& particlearray,std::vector<TypeRange>& typerangearray);
  void getBondlists(std::vector<CovalentBondInfo::BondList>& bondlistarray);
  void getWaterlist(WaterList& waterlist);
  void getShakelist(ShakeList& shakelist);


  template<class PA>
  void setParticle(const PA& particlearray,
		   const std::vector<TypeRange>& typerangearray);
  void setBondlists(const std::vector<CovalentBondInfo::BondList>& bondlistarray);
  void setWaterlist(const WaterList& waterlist);
  void setShakelist(const ShakeList& shakelist);

  void dump();
  void restore(size_t step=0);

 private:
  HDFSnapshot::HDFSnapshotFile *snapshot;
  int mode;
  long long num_atom;
};

#endif
