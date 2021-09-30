#ifndef SHORTRANGEINTERACTIONSET_H
#define SHORTRANGEINTERACTIONSET_H

#include "ShortRangeInteraction.h"


//! Interaction loop with Array of TypeRange
/*! T : Force/Potential class : Coulomb, EwaldReal, LJ, LJCoulomb, LJEwaldReal
  POTENTIALMODEL : PotentialModel (defined in Common.h) : depend on T : CoulombType for Coulomb or EwaldReal, LJType for LJ, LJCoulombType for LJCoulomb or LJEwaldReal
!!!! This version suport only between same type
!!!! TODO : Coulomb or EwaldReal interaction between LJCoulombType and CoulombType,  LJ interaction between LJCoulombType and LJType

i-set\ j-set LJ  LJCoulomb  Coulomb
LJ           LJ  LJ         No
LJCoulomb    LJ  LJCoulomb  Coulomb
Coulomb      No  Coulomb    Coulomb
In this class, only diagonal pair 
 */ 
template<ShortRange::InteractionType T, int POTENTIALMODEL, typename TPI, typename TPJ>
class ShortRangeInteractionSet {
public:

  //! self interaction for each set of i-particle
  void loopself_inset(const std::vector<TPI>& particleseti, 
                      const std::vector<TypeRange>& typerange, 
                      ForceArray& forceseti, 
                      double &energy);

  void loopself_inset(const std::vector<TPI>& particleseti, 
                      const std::vector<TypeRange>& typerange, 
                      ForceArray& forceseti, 
                      double &energy, 
                      const double cutoff2);

  //! self interaction for i-particle sets
  void loopself(const std::vector<TPI>& particleseti, 
                const std::vector<TypeRange>& typerange, 
                const std::vector< std::vector<int> >& targetindex,
                ForceArray& forceseti, 
                double &energy);

  void loopself(const std::vector<TPI>& particlei, 
                const std::vector<TypeRange>& typerange, 
                const std::vector< std::vector<int> >& targetindex,
                ForceArray& forceseti, 
                double &energy, 
                const double cutoff2);

  //! interaction between i-particle sets and j-particle sets
  void loopisetjsetij(const std::vector<TPI>& particleseti, 
                      const std::vector<TypeRange>& typerange, 
                      const std::vector<TPJ>& particlesetj,
                      const std::vector<TypeRange>& ghosttyperange, 
                      const std::vector< std::vector<int> >& ghosttargetindex,
                      ForceArray& forceseti, 
                      ForceArray& forcesetj, 
                      double &energy);

  void loopisetjsetij(const std::vector<TPI>& particleseti, 
                      const std::vector<TypeRange>& typerange, 
                      const std::vector<TPJ>& particlesetj,
                      const std::vector<TypeRange>& ghosttyperange, 
                      const std::vector< std::vector<int> >& ghosttargetindex,
                      ForceArray& forceseti, 
                      ForceArray& forcesetj, 
                      double &energy,
                      const double cutoff2);

  void loopijset(const std::vector<TPI>& particleseti, 
                 const std::vector<TypeRange>& typerange, 
                 const std::vector<TPJ>& particlesetj,
                 const std::vector<TypeRange>& ghosttyperange, 
                 const std::vector< std::vector<int> >& ghosttargetindex,
                 ForceArray& forceseti, 
                 ForceArray& forcesetj, 
                 double &energy);

  void loopijset(const std::vector<TPI>& particlei, 
                 const std::vector<TypeRange>& typerange, 
                 const std::vector<TPJ>& particlej,
                 const std::vector<TypeRange>& ghosttyperange, 
                 const std::vector< std::vector<int> >& ghosttargetindex,
                 ForceArray& forceseti, 
                 ForceArray& forcesetj, 
                 double &energy, 
                 const double cutoff2);

  void allloop(const std::vector<TPI>& particleseti,
               const std::vector<TypeRange>& typerange, 
               const std::vector< std::vector<int> >& targetindex,
               const std::vector<TPJ>& particlesetj,
               const std::vector<TypeRange>& ghosttyperange, 
               const std::vector< std::vector<int> >& ghosttargetindex,
               const std::vector< std::pair<int,int> >& ghostpairindex,
               ForceArray& forceseti, 
               ForceArray& forcesetj, 
               double &energyself, double &energy);

  void allloop(const std::vector<TPI>& particleseti,
               const std::vector<TypeRange>& typerange, 
               const std::vector< std::vector<int> >& targetindex,
               const std::vector<TPJ>& particlesetj,
               const std::vector<TypeRange>& ghosttyperange, 
               const std::vector< std::vector<int> >& ghosttargetindex,
               const std::vector< std::pair<int,int> >& ghostpairindex,
               ForceArray& forceseti, 
               ForceArray& forcesetj, 
               double &energyself, double &energy, const double cutoff2);
};

template<ShortRange::CoulombType CT, typename TPI, typename TPJ>
class AllTypeShortRangeInteractionSet {
public:
  LJ lj;
  LJCoulomb ljc;
  Coulomb cl;
  LJEwaldReal ljer;
  EwaldReal er;
#ifdef RDTSCTIME
  unsigned long long starttime;
  unsigned long long stoptime;
  unsigned long long ljtime, ljctime, cltime;
#endif

  void allloop(const std::vector<TPI>& particleseti,
               const std::vector<TypeRange>& typerange, 
               const std::vector< std::vector<int> >& targetindex,
               const std::vector<TPJ>& particlejset,
               const std::vector<TypeRange>& ghosttyperange, 
               const std::vector< std::vector<int> >& ghosttargetindex,
               const std::vector< std::pair<int,int> >& ghostpairindex,
               ForceArray& forceseti,
               ForceArray& forcesetj,
               double &energyself, double &energy);

  void allloop(const std::vector<TPI>& particleseti,
               const std::vector<TypeRange>& typerange, 
               const std::vector< std::vector<int> >& targetindex,
               const std::vector<TPJ>& particlesetj,
               const std::vector<TypeRange>& ghosttyperange, 
               const std::vector< std::vector<int> >& ghosttargetindex,
               const std::vector< std::pair<int,int> >& ghostpairindex,
               ForceArray& forceseti,
               ForceArray& forcesetj,
               double &energyself, double &energy, 
               const double cutoff2);

};

typedef ShortRangeInteractionSet<ShortRange::LJ, AnyLJPotential, Particle, Particle> LJInteractionSet;


typedef AllTypeShortRangeInteractionSet<ShortRange::OriginalCoulomb, Particle, Particle> CoulombAndLJInteractionSet;
typedef AllTypeShortRangeInteractionSet<ShortRange::ForEwaldReal, Particle, Particle> EwaldRealAndLJInteractionSet;

#endif
