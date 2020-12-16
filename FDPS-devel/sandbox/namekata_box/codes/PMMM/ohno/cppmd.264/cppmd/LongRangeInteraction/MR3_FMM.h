#ifndef MR3FMM_H
#define MR3FMM_H

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif  // HAVE_CONFIG_H
#include <mpi.h>
#include <map>
#include "Common.h"
#include "ParticleInfo.h"
#include "OperationSelector.h"
#include "LongRangeParameter.h"
//#include "../EXAFMM2/parallelfmm.h"

class ParallelFMM;

class MR3_FMM {
public:
  MR3_FMM(int, const LongRangeParameter& param, MPI_Comm) {
    xmax[0] = param.boxSize.x;
    xmax[1] = param.boxSize.y;
    xmax[2] = param.boxSize.z;
    allocated = false;
  }

  ~MR3_FMM();

  void initialize(size_t size, SpaceVector<int> celldiv3d, SpaceVector<int> nodediv3d);

  //! NOT YET
  template<class PA, class GPA>
    void calcforce(PA& particlearray,
		   std::vector<TypeRange>& typerangearray,
		   std::vector< std::vector<int> >& self_shorttarget_index,
		   GPA& ghost,
		   std::vector<TypeRange>& ghosttyperange, 
		   std::vector< std::vector<int> >& ghost_shorttarget_index,
		   std::vector< std::pair<int,int> >& ghost_pair_index,
		   ForceArray& shortforce, 
		   ForceArray& ghostshortforce,
		   double& shortenergyself,
		   double& shortenergy,
		   const double cutoff2);

  template<class PA, class GPA>
    void calccoulomb(PA& particlearray,
		     const std::vector<TypeRange>& typerangearray,
		     std::vector<int>& self_longset_index,
		      std::vector<int>& self_selfenergycell_index,
		     GPA& ghost,
		     const std::vector<TypeRange>& ghosttyperange, 
		     std::vector<int>& ghost_longset_index,
		     std::vector<int>& ghost_selfenergycell_index,
		     double& longenergy);

  double xmax[3];   //! max of position for compatibility to MR3

  int numnode3d[3];      // 
  int numcell3d[3];
  int numlocalcell3d[3];
  size_t num_cell;
  size_t num_local_cell;
 
  bool allocated;

  ParallelFMM *parallelfmm;
};
#endif
