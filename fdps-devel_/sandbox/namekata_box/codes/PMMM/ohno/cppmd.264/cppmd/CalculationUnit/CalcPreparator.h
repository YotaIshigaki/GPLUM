
#ifndef CALCPREPARATOR_H
#define CALCPREPARATOR_H

#include "MPIParallel.h"
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
#include "CalculationUnit.h"
#include "Config.h"

#include "Dump.h"

#include "TestParticle.h"

#include "Log.h"


namespace CalcPreparator_ {


/*
  INPUT for all
  node_id

  INPUT for makeoperations only
  bool withlong
  bool withbond

  INPUT for makeoperations and maketarget
  CellIndexType cellindextype

  INPUT for maketarget only
  int num_node :
  num_particle
  SpaceVector<double> cellmargin :
  
  INPUT for maketarget and constructcommunicator
  max_particle_in_cell



  INPUT for maketarget and constructcalculationunit
  SpaceVector<double> boxsize :
  cutoff

  INPUT-OUTPUT for maketarget
  num_total_set
  celldiv

  INPUT for makeparticle only
  ParticleArray& pa
  std::vector<PotentialModel>& pm
  //  SpaceVector<double> boxsize unused

  INPUT for makebond only
  CovalentBondList bl
  bool excludewaterbond
  
  INPUT for constructcommunicator only
  std::map<int,int> uidtorank

  INPUT for constructcalculationunit only
  ShortRange::CoulombType cltype
*/

void makeoperations(const int node_id,
                    const bool withlong,
                    const bool withbond,
		    const bool withcorrecttrans,
                    const CellIndexType cellindextype,
                    const LongRangeMPIPlan longplan=Combine0D,
                    const int num_long_only=0);
/*
void maketarget(const int node_id,
                int &short_id,
                std::vector<int>& node_id_of_shorts,
                const int num_node,
                const int num_particle,
                int& num_total_set,
                int& celldiv,
                const int max_particle_in_cell,
                const SpaceVector<double> boxsize, 
                const SpaceVector<double> cellmargin,
                const CellIndexType cellindextype,
                const double cutoff,
                const LongRangeMPIPlan longplan=Combine0D,
                const int num_long_only=0);
*/

void makegeometry(const int node_id,
                  const int opts_nodediv3d[],
                  const int opts_celldiv3d_in_node[],
                  const int num_node,
                  const int num_short,
                  int &short_id,
                  std::vector<int>& node_id_of_shorts,
                  int& num_total_set,
                  int& celldiv,
                  SpaceVector<int>& ndiv,
                  SpaceVector<int>& cdiv_in_node,
                  SpaceVector<int>& cdiv);
void maketarget(const int node_id,
                const int short_id,
                const int num_short,
                const int num_particle,
                const int max_particle_in_cell,
                const SpaceVector<double> boxsize, 
                const SpaceVector<double> cellmargin,
                const CellIndexType cellindextype,
                const double cutoff,
		const double boxsize_minimum_scale,
                int& num_total_set);

template<class PA>
void makeparticle(const int short_id,
                  const int total_num_particle,
                  const int total_num_freedom,
#if 1
                  const int num_copy,
                  const int num_copy_local,
#endif
                  const PA& pa, 
                  const WaterList& wl,
                  const ShakeList& sl,
                  const std::vector<PotentialModel>& pm,
                  const double ljcec);

void makebond(const int short_id,
              CovalentBondList bl,
              bool excludewaterbond,
              const ShakeList& sl,
              int shake_type);

//! set long range node geometry
void makelongrangegeometry(const LongRangeMPIPlan longplan,
                           const int num_long_only,
                           int& num_long,
                           GeometryXYZ& longgeometry);

/*! make MPI_COMM for longrange and set ID of long range node
  @attention set longrange communicator to mpi_comm_long
  @param[in] longplan select assign pattern of long range nodes
  @param[in] longgeometry geometry of long range nodes
  @param[in] num_long total number of long range nodes
  @param[in] num_short total number of short range nodes
  @param[in] mpi_comm_all parent MPI communicator
  @param[in] node_id rank at MPI communicator identified by mpi_comm_all
  @param[in] short_id identifier for short range or NO_SHORT_ID if this node is not in short range
  @param[out] long_id identifier for long range or NO_LONG_ID if this node is not in long range
 */
void makelongcomm(const LongRangeMPIPlan longplan,
                  const GeometryXYZ longgeometry,
                  const int num_long,
                  const int num_short,
                  const MPI_Comm mpi_comm_all,
                  const int node_id,
                  const int short_id,
                  int& long_id);

//! make list short cell id for long range node
/*!
  @param[in] longgeometry geometry of long range nodes
  @param[in] fringe fringe of space decomosition
  @param[in] long_id identifier for long range
  @param[out] cellid_list list of short cell identifier
 */
void makelongrangerequiredcell(const GeometryXYZ longgeometry,
                               const double fringe,
                               const int long_id,
                               std::vector<int>& cellid_list);

//! lists cell id to send long range node
void makecellidlistsendtolong(const GeometryXYZ longgeometry,
                              const double fringe,
                              const int num_long,
                              const int short_id,
                              std::vector< std::vector<int> >& send_list);

//! list cell id for selfenergy
void makeselfenergycell_list(const GeometryXYZ longgeometry,
                             const int long_id,
                             const std::vector<int>& cellid_list,
                             std::vector<int>& selfenergycell_list);

//! make MPI_Comm each long node and short nodes that send to long
void make_long_short_comm(const int num_node,
                          const int num_long,
                          const int node_id,
                          const int short_id,
                          const int long_id,
                          const MPI_Comm mpi_comm_all,
                          const std::vector< std::vector<int> >& send_to_long_list,
                          const std::vector<int>& long_reqcell_list,
                          std::vector<MPI_Comm>& mpi_comm_ls_list,
                          std::vector<int>& long_short_id,
                          std::vector<int>& sender_local_id,
                          std::map<int,int>& idinlong_to_longrank,
                          std::vector<int>& receiver_local_id,
                          std::vector< std::vector<int> >& long_recv_set_id_lists,
                          std::vector<int>& sender_global_id_list,
                          std::vector<int>& reciever_global_id_list);

//! set LongRangeParameter
void makelongrangeparameter(double cutoff,
                            double alpha,
                            double kCutoff,
                            int surfaceDipole,
                            SpaceVector<double> boxSize,
                            SpaceVector<double> gridLengths,
                            int order,
                            PMEType pmeType,
                            SpaceVector<int> grid_num,
                            GeometryXYZ node_geometry,
                            MGType multigridType,
                            int multigridIteration,
                            int multigridFFTlevel);

void constructcommunicator(const int short_id,
                           const int max_particle_in_cell,
                           const std::map<int,int>& shortidtorank,
                           const int num_long,
                           const int long_id,
                           const std::vector<MPI_Comm>& mpi_comm_ls_list,
                           const std::vector<int>& long_short_id,
                           const std::vector<int>& sender_local_id,
                           const std::vector< std::vector<int> >& long_recv_set_id_lists,
                           const std::vector<int>& receiver_local_rank,
                           const std::vector< std::vector<int> >& send_to_long_list,
                           const std::map<int,int>& idinlong_to_longrank,
                           const int short_comm_pattern,
                           const int move_comm_pattern);

template<class PA, class GPA>
    CalculationUnit<PA,GPA>* constructcalculationunit(const int node_id,
                                          const int short_id,
                                          const SpaceVector<double> boxsize,
                                          const CellIndexType cellindextype,
                                          const double cutoff,
#ifdef USE_PAIRLIST
                                          const double pmargin,
#endif
                                          const Config::TempCtrl& temp_control,
                                          const ShortRange::CoulombType cltype,
                                          const CovalentBondParameterList& covalent_bond_parameter_list);

//! construct set list of self set for long
void constructselflongset(const int long_id, const int short_id,
                          const std::vector<int>& short_set_list,
                          const std::vector< std::vector<int> >& send_to_long_list,
                          const std::vector<int>& selfenergycell_list,
                          std::vector<int>& longset_index,
                          std::vector<int>& selfenergycell_index);

void constructghostlongset(const int long_id, const int short_id,
                           const std::vector<int>& long_reqcell_list,
                           const std::vector< std::vector<int> >& long_recv_set_id_lists,
                           const std::vector<int>& short_set_list,
                           const std::vector<int>& selfenergycell_list,
                           std::vector<int>& ghostlongset,
                           std::vector<int>& ghost_selfenergycell_list);
#ifdef CPPMD_ENABLE_FMM
 void construct_fmm_target_cell(CalcForce &cf);
#endif  // CPPMD_ENABLE_FMM
#ifdef CPPMD_ENABLE_PMMM
 void construct_pmmm_target_cell(CalcForce &cf);
 void construct_pmmm(CalcForce &cf, int num_pp, int num_mm, double total_charge);
#endif  // CPPMD_ENABLE_PMMM

}

  
#endif // CALCPREPARATOR_H

