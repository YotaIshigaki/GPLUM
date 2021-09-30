#ifndef CACULATIONUNIT_H
#define CACULATIONUNIT_H

#include "MPIParallel.h"
#include "MPIParallelLongRange.h"
#include <map>
#include "Common.h"
#include "ParticleInfo.h"
#include "ShortRangeInteraction.h"
#include "ParallelCoordinator.h"
#include "CellIndex.h"
#include "CubicCell.h"
#include "HalfShell.h"
#include "CovalentBondInfo.h"
#include "CovalentBond.h"
#include "LongRangeInteraction.h"
#include "DummyLongRangeInteraction.h"
#include "CalcForce.h"
#include "Integrator.h"
#include "Config.h"

#include "Dump.h"

#include "TestParticle.h"

#include "Log.h"


//! minimam calculation unit for single node
template<class PA, class GPA>
class CalculationUnit {
public:
  int unit_identifier;                //!< unit ID
                                      //! may be same as node_id, rank
  int short_id;                       //!< ID in short communicator

  int num_total;                      //!< number of total particle
  int num_freedom;
  PA particlearray;        //!< partile data of thie node
  WaterList waterlist;
  ShakeList shakelist;
  std::vector<int> setid;
  std::map<int,int> setid_to_index;
  std::vector<TypeRange> typerangearray;
  std::vector<CovalentBondInfo::BondList> bondlistarray;
  std::vector<CovalentBondInfo::BondList> bondlistarray_idx;
  CovalentBondParameterList covalent_bond_parameter_list;
  std::vector< std::vector<int> > targetset;
  std::vector< std::pair<int,int> > ghostpairs;
  SpaceVector<double> boxsize;
  Integrator<PA,GPA> integrator;          //!< integrator of thie node
  CalcForce calcforce;            //!< force calclator
  OperationSelector operations;
  PostProcess postprocess;
  /*
    PostProcess has postprocess_receive for received j-particle (may be periodic shift), 
    merge for postprocess of move particle,
    select_move_particle and move_inside_node for preprocess of move particle.
    name miss match.
  */
  MPICommunicator communicator;    //!< communicator

  MPIReceiverForLong longreceiver;  //!< communicator to receive for long
  MPISenderForLong longsender;      //!< communicator to send for long

  //  CalculationUnit(){}

  CalculationUnit(MPICommunicator _communicator,
                  MPIReceiverForLong _lrecv,
                  MPISenderForLong _lsend,
                  PostProcess _postprocess,
                  ShortRange::CoulombType cltype,
                  const LongRangeParameter longrangeparam,
                  const int unitid,
                  const int shortid,
                  const SpaceVector<double> bs,
		  const double bms,
                  const MPI_Comm long_comm,
                  const Config::TempCtrl& temp_control
                  );

  CalculationUnit(PA& _particlearray, 
                  std::vector<int> sid,
                  std::vector<TypeRange> typerange,
                  WaterList _waterlist,
                  ShakeList _shakelist,
                  std::vector<CovalentBondInfo::BondList> bondlist,
                  std::vector<CovalentBondInfo::BondList> bondlist_idx,
                  const CovalentBondParameterList& cbplist,
                  std::vector< std::vector<int> > targetset,
                  std::vector< std::pair<int,int> > ghostpairs,
                  OperationSelector ope,
                  MPICommunicator& _communicator,
                  MPIReceiverForLong& _lrecv,
                  MPISenderForLong& _lsend,
                  PostProcess _postprocess,
                  const ShortRange::CoulombType cltype,
                  const LongRangeParameter longrangeparam,
                  const int unitid, const int shortid,
                  const SpaceVector<double> bs, 
		  const double bms,
		  const double cutoff, 
#ifdef USE_PAIRLIST
                  const double pmargin,
#endif
                  int num_t, int num_f,
                  const MPI_Comm long_comm,
                  const Config::TempCtrl& temp_control,
                  bool expire_reverse=false
                  );
  
  //  CalculationUnit& operator=(const CalculationUnit& cu);

  //! integrate
  bool startTimeIntegration(double dt, long tmax, 
                            const ShakeList& sl,
                            int shake_type, int shake_max_iterate, double shake_tolerance,
                            int reduce_interval, int printinter,
                            int dump_crd_inter, int dump_rst_inter, Dump& dump, Dump& dump_rst,
                            int moveinter,
                            int& write_restart,
			    long& last_t);

  void initialize(SpaceVector<int> celldiv3d, SpaceVector<int> nodediv3d){
    calcforce.initialize(particlearray.size(), celldiv3d, nodediv3d, operations);
  }

  size_t getParticleArray(ParticleArray &pa);
};




#endif // CACULATIONUNIT_H
