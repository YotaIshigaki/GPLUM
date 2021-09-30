#ifndef MR3INTERFACE_H
#define MR3INTERFACE_H

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif  // HAVE_CONFIG_H
#include <mpi.h>
#include <map>
#include "Common.h"
#include "ParticleInfo.h"
#include "OperationSelector.h"

class MR3 {
 public:
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
    void calcvdw(PA& particlearray,
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

  double *xi, *qi, *force, *ei;
  int *atypei;

  int ni, nj;
  int nat;
  double *rscale, *gscale;
  int tblno;
  double xmax;
  int periodicflag;

};
#endif
