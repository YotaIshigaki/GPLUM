#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "MPIParallel.h"
#include "MPIParallelLongRange.h"
#include <vector>
#include <cmath>
#include "Common.h"
#include "ParticleInfo.h"
#include "CovalentBondInfo.h"
#include "CalcForce.h"
#include "Timer.h"
#include "Settle.h"
#include "Config.h"

#include "Dump.h"
#ifdef USE_HDF
#include "HDFSnapshotData.h"
#endif

//! integrate
/*! Input  : Mass, Charge, Atomtype, ID
  I/O    : Position,  Velocuty, Force, total Energy
*/
template<class IPA, class IGPA>
class Integrator {
 public:
  double energy;
  double deltat;
  int unit_identifier;
  int short_id;

  std::vector<TypeRange> typerangearray;
  std::vector<CovalentBondInfo::BondList> bondlistarray;
  std::vector<CovalentBondInfo::BondList> bondlistarray_idx;
  ForceArray force;

  //  ParticleArray ghost;
  std::vector<ParticleRange> targetparticlerange;
  std::vector<TypeRange> targettyperange;
  std::vector<CovalentBondInfo::BondList> targetbond;
  std::vector<int> recvsetid;
  std::map<int,int> recvsetid_to_index;
  ForceArray ghostforce;

  // control interval
  int pi;
  int ri;
  int dci;
  int dri;


  SpaceVector<double> boxsize;
  double volume;
  SpaceVector<double> initial_boxsize;
  SpaceVector<double> inv_initial_boxsize;
  double boxsize_minimum_scale;
  double tljcec;
 
  // constant temerature, pressure
  double ref_temperature;
  // nose-hoover NVT
  double tau;
  double eta_inv_mass;
  double eta_pos;
  double eta_vel;
  double eta_force;
  // andersen-hoover NPT  MARK E. TUCKERMAN, et.al Molecular Physics, 1996, Vol.87, No.5, 1117-1157
  double ref_pressure;
  double tauv;
  double pressure;
  double diff_pressure_volume;
  double logv_mass;
  double logv_inv_mass;
  double logv_pos;
  double logv_vel;
  double logv_force;

  bool reset_ghost_longset_index;

  Settle settle;
  Config::TempCtrl temp_control;

  // Roll Back
  IPA backup_particle;

  std::vector<TypeRange> backup_typerangearray;
  std::vector<CovalentBondInfo::BondList> backup_bondlistarray;
  std::vector<CovalentBondInfo::BondList> backup_bondlistarray_idx;
  WaterList backup_waterlist;
  ShakeList backup_shakelist;
  long backup_time;
  int backup_pi;
  int backup_ri;
  int backup_dci;
  int backup_dri;
  double backup_eta_pos;
  double backup_eta_vel;
  double backup_eta_force;
  double backup_logv_pos;
  double backup_logv_vel;
  double backup_logv_force;
  SpaceVector<double> backup_boxsize;
  double backup_volume;
  double backup_tljcec;

  //! force calculation intervals
  /*!
    Indicate no force calculation time steps after calculated step
    0 : every time step
   */
  int mti_bond;   //!< bond calculation interval 
  int mti_short;  //!< short range calculation interval 
  int mti_long;   //!< long range calculation interval

  bool pre_calcforce;   //! true: calcforce before integration when initial data not include force (ex. Amber Restart).
#ifdef USE_HDF
  HDFDump hdfdump;
  HDFDump hdfrestore;
#endif
#ifdef BINARY_DUMP
  BinaryDump<IPA> binarydump;
  BinaryDump<IPA> binaryrestore;
#endif
  bool restore;

  int timer1;
  int timer_move; 
  int timer_exp;
  int timer_expcalc;
  int timer_calc;
  int timer_exf;
  int timer_rede;
  int timer_settle;
  int timer_shake;
  int timer_crd;
  int timer_rst;

  Integrator();

  Integrator(const int unitid, 
             const int shortid,
             std::vector<TypeRange>& tra,
             std::vector<CovalentBondInfo::BondList>& ba,
             std::vector<CovalentBondInfo::BondList>& ba_idx,
             const SpaceVector<double> bs,
	     const double bms,
             const Config::TempCtrl& tc
             );
  
  virtual ~Integrator() {}

  inline
  void mergerForces(CombinedParticleArray& particlearray);
  inline
  void mergerForces(ParticleArray& particlearray);

  /*
  inline
  void reduce_momentum(SpaceVector<double>& m_sum, 
                       SpaceVector<double>& am_sum,
                       const SpaceVector<double> momentum, 
                       const SpaceVector<double> angular_momentum);
  */

  template<typename PA>
  inline
  void reduce(double& ke, double& energy, Force& sf,
              const double lke, const double le,
              const PA& particlearray);

  inline
  void reduce(double& ke, double& energy,
              const double lke, const double le);
  inline
  void reduce(double& ke, double& energy, double &virial,
              const double lke, const double le, const double lv);
  inline
  void reduce(double& ke, double& energy, double& virial, double &total_mass, Position& center_of_mass,
	      const double lke, const double le, const double lv, const double tm, const Position cm);

  inline
  void dump_detail(const IPA& particlearray,
                   const WaterList& waterlist,
                   const long t);
#ifdef USE_HDF
  int open_hdfdump(const std::string& hdfname);
  void close_hdfdump();
  void dump_hdf_snapshot(const IPA& particlearray,
			 const std::vector<TypeRange>& typerangearray,
			 const std::vector<CovalentBondInfo::BondList>& bondlistarray,
			 const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
			 const WaterList& waterlist,
			 const ShakeList& shakelist,
			 const double kenergy,
			 const double penergy,
			 const double virial,
			 const double total_penergy,
			 const long t,
			 const long num_total,
			 const int num_rank);
  int open_hdfrestore(const std::string& hdfname);
  void close_hdfrestore();
  void restore_hdf_snapshot(IPA& particlearray,
			    std::vector<TypeRange>& typerangearray,
			    std::vector<CovalentBondInfo::BondList>& bondlistarray,
			    std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
			    WaterList& waterlist,
			    ShakeList& shakelist,
			    double& kenergy,
			    double& penergy,
			    double& virial,
			    double& total_penergy,
			    long& t,
			    const long num_total,
			    const int num_rank);
#endif

#ifdef BINARY_DUMP
  inline
  void dump_binary_snapshot(const IPA& particlearray,
			    const std::vector<TypeRange>& typerangearray,
			    const std::vector<CovalentBondInfo::BondList>& bondlistarray,
			    const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
			    const WaterList& waterlist,
			    const ShakeList& shakelist,
			    const double kenergy,
			    const double penergy,
			    const double virial,
			    const double total_penergy,
			    const long t);
  inline
  void restore_binary_snapshot(IPA& particlearray,
			       std::vector<TypeRange>& typerangearray,
			       std::vector<CovalentBondInfo::BondList>& bondlistarray,
			       std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
			       WaterList& waterlist,
			       ShakeList& shakelist,
			       double& kenergy,
			       double& penergy,
			       double& virial,
			       double& total_penergy,
			       long& t);
#endif

  inline
  void backup_for_rollback(const IPA& particlearray,
			   const std::vector<TypeRange>& typerangearray,
			   const std::vector<CovalentBondInfo::BondList>& bondlistarray,
			   const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
			   const WaterList& waterlist,
			   const ShakeList& shakelist,
			   const long t);
  inline
  void rollback_from_backup(IPA& particle,
			    std::vector<TypeRange>& typerangearray,
			    std::vector<CovalentBondInfo::BondList>& bondlistarray,
			    std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
			    WaterList& waterlist,
			    ShakeList& shakelist,
			    long& time);
  inline
  bool check_move_over_and_rollback(const bool not_over_flow,
				    const int over_move,
				    const int max_moveinter,
				    IPA& particlearray,
				    std::vector<TypeRange>& typerangearray,
				    std::vector<CovalentBondInfo::BondList>& bondlistarray,
				    std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
				    WaterList& waterlist,
				    ShakeList& shakelist,
				    long& t,
				    int& moveinter);
  inline
  void velocity_verlet_pos_vel(IPA& particlearray,
			       const std::vector<TypeRange>& typerangearray,
			       const double deltat);

  inline
  void velocity_verlet_npt_pos_vel(IPA& particlearray,
				   const std::vector<TypeRange>& typerangearray,
				   const double deltat);

  inline
  void velocity_verlet_npt_pos_vel_vtcorr(IPA& particlearray,
					  const std::vector<TypeRange>& typerangearray,
					  const double deltat, const double inv_nt);

  inline
  void velocity_verlet_velocity(IPA& particlearray,
				const std::vector<TypeRange>& typerangearray,
				const double deltat);
  
  bool time_integrations(IPA& particlearray, 
                                 std::vector<CovalentBondInfo::BondList>& bondlistarray,
                                 std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
                                 WaterList& waterlist,
                                 ShakeList& shakelist,
                                 std::vector<int>& setid,
                                 std::map<int,int>& setid_to_index,
                                 CalcForce& calcforce, 
                                 double dt, const long tmax, 
                                 OperationSelector operations,
                                 MPICommunicator& communicator,
                                 MPIReceiverForLong& longreceiver,
                                 MPISenderForLong& longsender,
                                 PostProcess& postprocess,
                                 const ShakeList& sl,
                                 int shake_type, int shake_max_iterate, double shake_tolerance,
                                 int reduce_interval, int printinterval,
                                 int dump_crd_inter, int dump_rst_inter, Dump& dump, Dump& dump_rst,
                                 int moveinterval,
                                 int num_total,
                                 int num_freedom,
                                 int& dump_restart,
				 long& last_t);

#ifdef BINARY_DUMP
  int open_binarydump(const std::string& filename);
  int open_binaryrestore(const std::string& filename);
  void close_binarydump();
  void close_binaryrestore();
#endif

 private:
};

#endif
