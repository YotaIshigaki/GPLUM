#ifndef CALCFORCE_H
#define CALCFORCE_H

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif  // HAVE_CONFIG_H
#include <mpi.h>
#include <map>
#include "Common.h"
#include "ParticleInfo.h"
#include "OperationSelector.h"
#include "ShortRangeInteraction.h"
#include "ShortRangeInteractionSet.h"
#include "CovalentBondInfo.h"
#include "CovalentBond.h"
#include "LongRangeInteraction.h"
#include "DummyLongRangeInteraction.h"
#ifdef USE_MR3
#include "mr3_interface.h"
#endif

#ifdef CPPMD_ENABLE_LONGRANGE
# if defined(CPPMD_ENABLE_FMM)
#  ifdef CPPMD_ENABLE_MR3EXAFMM
#   include "MR3_FMM.h"
#  else
#   include "FMMLongRangeInteraction.h"
#  endif
# elif defined(CPPMD_ENABLE_PME)
#   include "PME.h"
# elif defined(CPPMD_ENABLE_OLDPME)
#   include "PMEInterface.h"
#   include "PMELongRangeInteraction.h"
# elif defined(CPPMD_ENABLE_EWALD)
#   include "EwaldInterface.h"
#   include "EwaldLongRangeInteraction.h"
# elif defined(CPPMD_ENABLE_STMEWALD)
#   include "StMELongRangeInteraction.h"
# elif defined(CPPMD_ENABLE_PMMM)
#   include "PMMMInteraction.h"
# endif
#endif  // CPPMD_ENABLE_LONGRANGE

#ifdef USE_PAIRLIST
#include "PairList.h"
#include "PairListInteraction.h"
#include "PairListZeroDipole.h"
#endif


//! Calculate force and potential energy
// Input : Position, Charge, Atomtype, ID
// Output : Force, total Energy
//template<class CBPARTICLE>
class CalcForce {
public:
  int unit_identifier;
  int short_id;

  std::vector<int> shortset;
  std::vector< std::vector<int> > self_shorttarget_set;
  std::vector< std::vector<int> > self_shorttarget_index;
  std::vector<int> longset;
  std::vector<int> self_longset_index;
  std::vector<int> self_selfenergycell_index;
  ParticleArray ghost;
  std::vector< std::vector<int> > ghost_shorttarget_set;
  std::vector< std::vector<int> > ghost_shorttarget_index;
  std::vector<int> ghostlongset;
  std::vector<int> ghost_selfenergycell_list;
  std::vector<int> ghost_longset_index;
  std::vector<int> ghost_selfenergycell_index;
#ifdef CPPMD_ENABLE_FMM
  std::vector< std::vector<int> > self_fmmtarget_set;
  std::vector< std::vector<int> > self_fmmtarget_index;
  std::vector< std::vector<int> > ghost_fmmtarget_set;
  std::vector< std::vector<int> > ghost_fmmtarget_index;
#endif  // CPPMD_ENABLE_FMM
#ifdef CPPMD_ENABLE_PMMM
  std::vector< std::vector<int> > self_pmmmtarget_set;
  std::vector< std::vector<int> > self_pmmmtarget_index;
  std::vector< std::vector<int> > ghost_pmmmtarget_set;
  std::vector< std::vector<int> > ghost_pmmmtarget_index;
  std::vector< waterexclude > waterexcludelist;
#endif
  CovalentBondList bondlist;
  LJInteractionSet lj;
  EwaldRealAndLJInteractionSet ewald_and_lj;
  CoulombAndLJInteractionSet cl_and_lj;
  LongRangeParameter longrangeparam;
  DummyLongRangeInteraction longrange;

#ifdef CPPMD_ENABLE_LONGRANGE
# if defined(CPPMD_ENABLE_FMM)
#  ifdef CPPMD_ENABLE_MR3EXAFMM
  MR3_FMM mr3fmm;
#  else
  FMMLongRangeInteraction fmmlongrange;
#  endif
# elif defined(CPPMD_ENABLE_PME)
  PMEModuleInterface *pmemoduleinterface;
  PMELongRangeInteraction pmelongrange;
# elif defined(CPPMD_ENABLE_OLDPME)
  PMEModule::PMEInterface *pmeinterface;
  PMELongRangeInteraction pmelongrange;
# elif defined(CPPMD_ENABLE_EWALD)
  EwaldModule::EwaldInterface *ewaldinterface;
  EwaldLongRangeInteraction ewaldlongrange;
# elif defined(CPPMD_ENABLE_STMEWALD)
  StMELongRangeInteraction stmelongrange;
# elif defined(CPPMD_ENABLE_PMMM)
  PMMMLocalInteraction pmmmlocal;
  PMMMLongRangeInteraction pmmmlongrange;
# endif
#endif  // CPPMD_ENABLE_LONGRANGE

#ifdef OLDPARTICLE
  CovalentBond<CBModule::CBInterface::ParticleLocation> covalentbond;
#else
  CovalentBond<CBModule::CBInterface::ParticleIndex> covalentbond;
#endif
  ShortRange::CoulombType cltype;
  std::vector< std::pair<int,int> > ghost_pair_set;
  std::vector< std::pair<int,int> > ghost_pair_index;

  double cutoff;
  double cutoff2;

  ForceArray shortforce;
  ForceArray ghostshortforce;
  double shortenergyself;
  double shortenergy;
  double shortenergyexcluded;
  double longenergy;
  double bondenergy;
  double shortvirial;
  double longvirial;

#ifdef USE_PAIRLIST
  double pmargin2;

# ifdef INDEX_PAIRLIST
  int (*pid)[MAX_PAIR];
  int (*plj)[MAX_PAIR];
  int *npair;
  int *iid;
  double (*fi)[3];
  int npl;

  int (*selfpid)[MAX_PAIR];
  int (*selfplj)[MAX_PAIR];
  int *selfnpair;
  int *selfiid;
  double (*selffi)[3];
  int selfnpl;

# else  // INDEX_PAIRLIST
  int (*pid)[MAX_PAIR];
  double (*pos)[MAX_PAIR][3];
  double (*charge)[MAX_PAIR];
  double (*plj)[MAX_PAIR][4];
  int *npair;
  int *iid;
  double (*posi)[3];
  double *chargei;
  double (*fi)[3];
  int npl;

  int (*selfpid)[MAX_PAIR];
  double (*selfpos)[MAX_PAIR][3];
  double (*selfcharge)[MAX_PAIR];
  double (*selfplj)[MAX_PAIR][4];
  int *selfnpair;
  int *selfiid;
  double (*selfposi)[3];
  double *selfchargei;
  double (*selffi)[3];
  int selfnpl;

# endif  // INDEX_PAIRLIST
#endif  // USE_PAIRLIST

  int timer_short;
  int timer_long;
  int timer_cb;
  int timer_cellpair;
  int timer_long1;
  int timer_long2;
# ifdef CPPMD_ENABLE_PMMM
  int timer_pmmm_pm_top_half;
  int timer_pmmm_pm_comm;
  int timer_pmmm_pm_bottom_half;
  int timer_pmmm_mm_prep;
  int timer_pmmm_mm_calc;
  int timer_pmmm_mm_post;
# endif


  //  CalcForce(){}

  CalcForce(std::vector<int>& sset,
            const MPI_Comm long_comm=MPI_COMM_WORLD);

  CalcForce(std::vector<int> sset, const ShortRange::CoulombType _cltype,
            const LongRangeParameter lrp, 
            const int unitid, const int shortid, const double cf=0.0,
#ifdef USE_PAIRLIST
            const double pm=1.0,
            const int num_list=1,
#endif  // USE_PAIRLIST
            const MPI_Comm long_comm=MPI_COMM_WORLD);

  CalcForce(std::vector<int> sset, std::vector< std::vector<int> > targetset, 
            std::vector< std::pair<int,int> > ghostpairs,
            const ShortRange::CoulombType _cltype,
            const LongRangeParameter lrp,
            const int unitid, const int shortdid, const double cf=0.0,
#ifdef USE_PAIRLIST
            const double pm=1.0,
            const int num_list=1,
#endif  // USE_PAIRLIST
            const MPI_Comm long_comm=MPI_COMM_WORLD,
            bool expire_reverse=false);

  void register_timers();

  void set_cutoff(double cf);
#ifdef USE_PAIRLIST
  void set_pmargin(double pm);
  void allocate_list(size_t num_particle);
#endif  // USE_PAIRLIST

  void clear_shorttarget_set();

  void generate_shorttarget_set(std::vector< std::vector<int> >& targetset,
                                bool expire_reverse=false);
#ifdef CPPMD_ENABLE_FMM
  void clear_fmmtarget_set();
  void generate_fmmtarget_set(std::vector< std::vector<int> >& fmmtargetset,
			      bool expire_reverse=false);
#endif  // CPPMD_ENABLE_FMM
#ifdef CPPMD_ENABLE_PMMM
  void clear_pmmmtarget_set();
  void generate_pmmmtarget_set(std::vector< std::vector<int> >& pmmmtargetset,
			      bool expire_reverse=false);
#endif  // CPPMD_ENABLE_PMMM

  
  void convert_setid_to_index(std::map<int,int>& id_to_index,
                              const std::vector<int>& id,
                              std::vector<int>& index);
  void selfshort_zero();
  void ghostshort_zero();
  void short_zero();

  void long_zero();

  void convert_setid_to_index(std::map<int,int>& id_to_index,
                              const std::vector< std::pair<int,int> >& pairid,
                              std::vector< std::pair<int,int> >& pairindex);

  void set_ghost_longset_index(const std::map<int,int>& id_to_index);

  void initialize(size_t pa_size, SpaceVector<int> celldiv3d, SpaceVector<int> nodediv3d,OperationSelector operations);

  template<class PA>
    double
    calcCorrectionExcludedWater_full_fixed(const PA& particlearray,
					   const WaterList& waterlist,
					   const double cutoff2,
					   const ShortRange::CoulombType clt);

#ifdef CPPMD_ENABLE_PMMM
  template<class PA>
  void calcPMMM_top_half(const PA& particlearray, 
			 const std::vector<TypeRange>& typerangearray,
			 OperationSelector operations);

  void make_waterexcludelist(const size_t num,
			     const WaterList& waterlist,
			     std::vector< waterexclude > &wex);
#endif

  template<class PA, class GPA>
  void calcForce_local(PA& particlearray, 
		       const std::vector<TypeRange>& typerangearray, 
		       const std::vector<CovalentBondInfo::BondList>& bondlistarray,
		       const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
		       const WaterList& waterlist,
		       GPA& ghost, 
		       const std::vector<int>& ghostsetid, 
		       std::map<int,int>& recvsetid_to_index, 
		       const std::vector<TypeRange>& ghosttyperange, 
		       const std::vector<CovalentBondInfo::BondList>& ghostbond,
		       ForceArray& force, 
		       ForceArray& ghostforce, 
		       double& energy,
		       double& virial,
		       OperationSelector operations);
  template<class PA, class GPA>
  void calcLongForce(PA& particlearray, 
			       const std::vector<TypeRange>& typerangearray, 
			       GPA& ghost, 
			       const std::vector<TypeRange>& ghosttyperange);
  template<class PA, class GPA>
  void calcForce(PA& particlearray, 
                 const std::vector<TypeRange>& typerangearray, 
                 const std::vector<CovalentBondInfo::BondList>& bondlistarray,
                 const std::vector<CovalentBondInfo::BondList>& bondlistarray_idx,
		 const WaterList& waterlist,
                 GPA& ghost, 
                 const std::vector<int>& ghostsetid, 
                 std::map<int,int>& recvsetid_to_index, 
                 const std::vector<TypeRange>& ghosttyperange, 
                 const std::vector<CovalentBondInfo::BondList>& ghostbond,
                 ForceArray& force, 
                 ForceArray& ghostforce, 
                 double& energy,
		 double& virial,
                 OperationSelector operations,
		 bool with_local=true);
};
#endif  // CALCFORCE_H
