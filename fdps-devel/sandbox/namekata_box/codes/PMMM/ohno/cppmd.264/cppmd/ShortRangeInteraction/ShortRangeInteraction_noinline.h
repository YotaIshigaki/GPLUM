#include "ShortRangeInteraction.h"

template<class T, typename TPI, typename TPJ>
void loopj(T& interaction, TPI& particlei, std::vector<TPJ>& particlej,
           size_t begin, size_t end, Force& forcei, double& energyi)
{
  Force force(0.0,0.0,0.0);
  double energy=0.0;
  for(long j=begin;j<end;j++){
    double energyij = 0.0;
    double dp;
    Position d = particlei.position - particlej[j].position;
    double r2 = d.norm2();
    interaction.Interaction(r2,particlei,particlej[j],energyij,dp);
    Force forceij = d*dp;
    energy += energyij;
    particlej[j].force -= forceij;
    force += forceij;
  }
  forcei += force;
  energyi += energy;
}

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
template<class T, int RANGETYPE, typename TPI, typename TPJ>
class ShortRangeInteractionSet {
public:
  T interaction;

  //! self interaction for each set of i-particle
  void loopself_inset(std::vector<TPI>& particlei, 
                      std::vector<TypeRange>& typerange, 
                      double &energy) {
    for(size_t si=0;si<typerange.size();si++){
      ParticleRange r = getrange<RANGETYPE>(typerange[si]);
      for(long i=r.begin;i<r.end;i++){
        Force forcei(0.0,0.0,0.0);
        double energyi=0.0;
        for(long j=i+1;j<r.end;j++){
          double energyij = 0.0;
          double dp;
          Position d = particlei[i].position - particlei[j].position;
          double r2 = d.norm2();
          interaction.Interaction(r2,particlei[i],particlei[j],energyij,dp);
          Force forceij = d*dp;
          energyi += energyij;
          particlei[j].force -= forceij;
          forcei += forceij;
        }
        particlei[i].force += forcei;
        energy += energyi;
      }
    }
  }
  //! self interaction for i-particle sets
  void loopself(std::vector<TPI>& particlei, 
                std::vector<TypeRange>& typerange, 
                std::vector< std::vector<int> >& targetindex,
                double &energy) {
    loopself_inset(particlei, typerange, energy);
    for(size_t si=0;si<typerange.size();si++){
      ParticleRange r = getrange<RANGETYPE>(typerange[si]);
      for(size_t t=0;t<targetindex[si].size();t++){
        int ti = targetindex[si][t];
        if(ti>si){
          ParticleRange tr = getrange<RANGETYPE>(typerange[ti]);
          for(long i=r.begin;i<r.end;i++){
            Force forcei(0.0,0.0,0.0);
            double energyi=0.0;
            for(long j=tr.begin;j<tr.end;j++){
              double energyij = 0.0;
              double dp;
              Position d = particlei[i].position - particlei[j].position;
              double r2 = d.norm2();
              interaction.Interaction(r2,particlei[i],particlei[j],energyij,dp);
              Force forceij = d*dp;
              energyi += energyij;
              particlei[j].force -= forceij;
              forcei += forceij;
            }
            particlei[i].force += forcei;
            energy += energyi;
          }
        }
      }
    }
  }

  //! interaction between i-particle sets and j-particle sets
  void loopij(std::vector<TPI>& particlei, 
              std::vector<TypeRange>& typerange, 
              std::vector<TPJ>& particlej,
              std::vector<TypeRange>& ghosttyperange, 
              std::vector< std::vector<int> >& ghosttargetindex,
              double &energy) {
   for(size_t si=0;si<typerange.size();si++){
     ParticleRange r = getrange<RANGETYPE>(typerange[si]);
      for(long i=r.begin;i<r.end;i++){
        Force forcei(0.0,0.0,0.0);
        double energyi=0.0;
        for(size_t t=0;t<ghosttargetindex[si].size();t++){
          int ti = ghosttargetindex[si][t];
          ParticleRange tr = getrange<RANGETYPE>(ghosttyperange[ti]);
          for(long j=tr.begin;j<tr.end;j++){
            double energyij = 0.0;
            double dp;
            Position d = particlei[i].position - particlej[j].position;
            double r2 = d.norm2();
            interaction.Interaction(r2,particlei[i],particlej[j],energyij,dp);
            Force forceij = d*dp;
            energyi += energyij;
            particlej[j].force -= forceij;
            forcei += forceij;
          }
        }
        particlei[i].force += forcei;
        energy += energyi;
      }
    }
  }
};

