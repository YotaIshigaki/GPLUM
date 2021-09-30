#include "ShortRangeInteractionSet.h"

// without reaction force, cutoff
template<ShortRange::InteractionType T, typename TPI, typename TPJ>
static void loopj(const TPI& particlei, const std::vector<TPJ>& particlesetj,
                  int begin, int end, 
                  Force& forcei, double& energyi)
{
  double fx=0.0, fy=0.0, fz=0.0;
  double energy=0.0;
  int j;
  //#pragma omp parallel for private(j) reduction(+ : energy,fx,fy,fz)
  for(j=begin;j<end;j++){
    double energyij = 0.0;
    double dp;
    Position d = particlei.position - particlesetj[j].position;
    double r2 = d.norm2();
    Interaction<T>(r2,particlei,particlesetj[j],energyij,dp);
    Force forceij = d*dp;
    energy += energyij;
    fx += forceij.x;
    fy += forceij.y;
    fz += forceij.z;
  }
  forcei.x += fx;
  forcei.y += fy;
  forcei.z += fz;
  energyi += energy;
}

template<ShortRange::InteractionType T, typename TPI, typename TPJ>
static void loopij(const std::vector<TPI>& particleseti, size_t i_begin, size_t i_end, 
                   const std::vector<TPJ>& particlesetj, size_t j_begin, size_t j_end,
                   ForceArray& forceseti, double &energy)
{
  for(size_t i=i_begin;i<i_end;i++){
    Force force(0.0,0.0,0.0);
    double energyi=0.0;
    loopj<T>(particleseti[i], particlesetj, j_begin, j_end,
             force,  energyi);
    forceseti[i] += force;
    energy += energyi;
  }
}


template<ShortRange::InteractionType T, typename TPI, typename TPJ>
static void loopj(const TPI& particlei, const std::vector<TPJ>& particlesetj,
                  int begin, int end, 
                  Force& forcei, ForceArray& forcesetj, double& energyi)
{
  double fx=0.0, fy=0.0, fz=0.0;
  double energy=0.0;
  int j;
  //#pragma omp parallel for private(j) reduction(+ : energy,fx,fy,fz)
  for(j=begin;j<end;j++){
    double energyij = 0.0;
    double dp;
    Position d = particlei.position - particlesetj[j].position;
    double r2 = d.norm2();
    Interaction<T>(r2,particlei,particlesetj[j],energyij,dp);
    Force forceij = d*dp;
    energy += energyij;
    forcesetj[j] -= forceij;
    fx += forceij.x;
    fy += forceij.y;
    fz += forceij.z;
  }
  forcei.x += fx;
  forcei.y += fy;
  forcei.z += fz;
  energyi += energy;
}

template<ShortRange::InteractionType T, typename TPI, typename TPJ>
static void loopij(const std::vector<TPI>& particleseti, size_t i_begin, size_t i_end, 
                   const std::vector<TPJ>& particlesetj, size_t j_begin, size_t j_end,
                   ForceArray& forceseti, ForceArray& forcesetj, double &energy)
{
  for(size_t i=i_begin;i<i_end;i++){
    Force force(0.0,0.0,0.0);
    double energyi=0.0;
    loopj<T>(particleseti[i], particlesetj, j_begin, j_end,
             force, forcesetj, energyi);
    forceseti[i] += force;
    energy += energyi;
  }
}

template<ShortRange::InteractionType T, int POTENTIALMODEL, typename TPI, typename TPJ>
static void loopjset(const TPI& particlei, const std::vector<TPJ>& particlesetj,
                     const std::vector<TypeRange>& ghosttyperange, 
                     const std::vector<int>& ghosttargetindex,
                     Force& forcei, ForceArray& forcesetj, double& energy)
{
  for(size_t t=0;t<ghosttargetindex.size();t++){
    int ti = ghosttargetindex[t];
    int j_begin = getbegin<POTENTIALMODEL>(ghosttyperange[ti]);
    int j_end = getend<POTENTIALMODEL>(ghosttyperange[ti]);
    Force force(0.0,0.0,0.0);
    double energyi=0.0;
    loopj<T>(particlei, particlesetj, j_begin, j_end,
             force, forcesetj, energyi);
    forcei += force;
    energy += energyi;
  }
}

template<ShortRange::InteractionType T, typename TPI, typename TPJ>
static void loopj(const TPI& particlei, const std::vector<TPJ>& particlesetj,
                  int begin, int end, 
                  Force& forcei, double& energyi,
                  const double cutoff2)
{
  double fx=0.0,fy=0.0,fz=0.0;
  double energy=0.0;
  int j;
  //#pragma omp parallel for private(j) reduction(+ : energy,fx,fy,fz)
  for(j=begin;j<end;j++){
    double energyij = 0.0;
    double dp;
    Position d = particlei.position - particlesetj[j].position;
    double r2 = d.norm2();
    if(r2<cutoff2){
      Interaction<T>(r2,particlei,particlesetj[j],energyij,dp);
      Force forceij = d*dp;
      energy += energyij;
      fx += forceij.x;
      fy += forceij.y;
      fz += forceij.z;
    }
  }
  forcei.x += fx;
  forcei.y += fy;
  forcei.z += fz;
  energyi += energy;
}

template<ShortRange::InteractionType T, typename TPI, typename TPJ>
static void loopj(const TPI& particlei, const std::vector<TPJ>& particlesetj,
                  int begin, int end, 
                  Force& forcei, ForceArray& forcesetj, double& energyi,
                  const double cutoff2)
{
  double fx=0.0,fy=0.0,fz=0.0;
  double energy=0.0;
  int j;
  //#pragma omp parallel for private(j) reduction(+ : energy,fx,fy,fz)
  for(j=begin;j<end;j++){
    double energyij = 0.0;
    double dp;
    Position d = particlei.position - particlesetj[j].position;
    double r2 = d.norm2();
    if(r2<cutoff2){
      Interaction<T>(r2,particlei,particlesetj[j],energyij,dp);
      Force forceij = d*dp;
      energy += energyij;
      forcesetj[j] -= forceij;
      fx += forceij.x;
      fy += forceij.y;
      fz += forceij.z;
    }
  }
  forcei.x += fx;
  forcei.y += fy;
  forcei.z += fz;
  energyi += energy;
}

template<ShortRange::InteractionType T, typename TPI, typename TPJ>
static void loopij(const std::vector<TPI>& particleseti, size_t i_begin, size_t i_end, 
                   const std::vector<TPJ>& particlesetj, size_t j_begin, size_t j_end,
                   ForceArray& forceseti, ForceArray& forcesetj, double &energy, 
                   const double cutoff2)
{
  for(size_t i=i_begin;i<i_end;i++){
    Force force(0.0,0.0,0.0);
    double energyi=0.0;
    loopj<T>(particleseti[i], particlesetj, j_begin, j_end,
             force, forcesetj, energyi, cutoff2);
    forceseti[i] += force;
    energy += energyi;
  }
}

template<ShortRange::InteractionType T, int POTENTIALMODELI, int POTENTIALMODELJ, typename TPJ>
static void looppairset(const std::vector<TPJ>& particlesetj,
                        const std::vector<TypeRange>& ghosttyperange, 
	                   const std::vector<std::pair<int, int> >& ghostpairindex,
                    ForceArray& forcesetj,
	                double& energy)
{
  for(std::vector<std::pair<int, int> >::const_iterator
      it = ghostpairindex.begin();it != ghostpairindex.end();++it) {
    int si = it->first;
    int sj = it->second;
    int i_begin = getbegin<POTENTIALMODELI>(ghosttyperange[si]);
    int i_end = getend<POTENTIALMODELI>(ghosttyperange[si]);
    int j_begin = getbegin<POTENTIALMODELJ>(ghosttyperange[sj]);
    int j_end = getend<POTENTIALMODELJ>(ghosttyperange[sj]);
    loopij<T>(particlesetj, i_begin, i_end, particlesetj, j_begin, j_end, forcesetj, forcesetj,energy);
  }
}

template<ShortRange::InteractionType T, int POTENTIALMODELI, int POTENTIALMODELJ, typename TPJ>
static void looppairset(const std::vector<TPJ>& particlesetj,
                        const std::vector<TypeRange>& ghosttyperange, 
	                const std::vector<std::pair<int, int> >& ghostpairindex,
                    ForceArray& forcesetj,
	                double& energy, const double cutoff2)
{
  for(std::vector<std::pair<int, int> >::const_iterator
      it = ghostpairindex.begin();it != ghostpairindex.end();++it) {
    int si = it->first;
    int sj = it->second;
    int i_begin = getbegin<POTENTIALMODELI>(ghosttyperange[si]);
    int i_end = getend<POTENTIALMODELI>(ghosttyperange[si]);
    int j_begin = getbegin<POTENTIALMODELJ>(ghosttyperange[sj]);
    int j_end = getend<POTENTIALMODELJ>(ghosttyperange[sj]);
    loopij<T>(particlesetj, i_begin, i_end, particlesetj, j_begin, j_end,
              forcesetj, forcesetj,
              energy, cutoff2);
  }
}

template<ShortRange::CoulombType CT, typename TPJ>
static void looppairset(const std::vector<TPJ>& particlesetj,
                        const std::vector<TypeRange>& ghosttyperange, 
	                const std::vector<std::pair<int, int> >& ghostpairindex,
                     ForceArray& forcesetj,
	                double& energy)
{
  looppairset<ShortRange::LJ,OnlyLJPotential,AnyLJPotential,TPJ>
    (particlesetj, ghosttyperange, ghostpairindex, forcesetj, energy);
  if(CT==ShortRange::OriginalCoulomb){
    looppairset<ShortRange::LJCoulomb,
                LJCoulombPotential,LJCoulombPotential,TPJ>
      (particlesetj, ghosttyperange, ghostpairindex, forcesetj, energy);
    looppairset<ShortRange::Coulomb,
                LJCoulombPotential,OnlyCoulombPotential,TPJ>
      (particlesetj, ghosttyperange, ghostpairindex, forcesetj, energy);
    looppairset<ShortRange::Coulomb,
                OnlyCoulombPotential,AnyCoulombPotential,TPJ>
      (particlesetj, ghosttyperange, ghostpairindex, forcesetj, energy);
  }else if(CT==ShortRange::ForEwaldReal){
    looppairset<ShortRange::LJEwaldReal,
                LJCoulombPotential,LJCoulombPotential,TPJ>
      (particlesetj, ghosttyperange, ghostpairindex, forcesetj, energy);
    looppairset<ShortRange::EwaldReal,
                LJCoulombPotential,OnlyCoulombPotential,TPJ>
      (particlesetj, ghosttyperange, ghostpairindex, forcesetj, energy);
    looppairset<ShortRange::EwaldReal,
                OnlyCoulombPotential,AnyCoulombPotential,TPJ>
      (particlesetj, ghosttyperange, ghostpairindex, forcesetj, energy);
  }
}

template<ShortRange::CoulombType CT, typename TPJ>
static void looppairset(const std::vector<TPJ>& particlesetj,
                        const std::vector<TypeRange>& ghosttyperange, 
	                const std::vector<std::pair<int, int> >& ghostpairindex,
                     ForceArray& forcesetj,
	                double& energy, const double cutoff2)
{
  looppairset<ShortRange::LJ,OnlyLJPotential,AnyLJPotential,TPJ>
    (particlesetj, ghosttyperange, ghostpairindex, forcesetj, energy, cutoff2);

  // TODO missing LJC - LJ

  if(CT==ShortRange::OriginalCoulomb){
    looppairset<ShortRange::LJCoulomb,
                LJCoulombPotential,LJCoulombPotential,TPJ>
      (particlesetj, ghosttyperange, ghostpairindex, forcesetj, energy, cutoff2);
    looppairset<ShortRange::Coulomb,
                LJCoulombPotential,OnlyCoulombPotential,TPJ>
      (particlesetj, ghosttyperange, ghostpairindex, forcesetj, energy, cutoff2);
    looppairset<ShortRange::Coulomb,
                OnlyCoulombPotential,AnyCoulombPotential,TPJ>
      (particlesetj, ghosttyperange, ghostpairindex, forcesetj, energy, cutoff2);
  }else if(CT==ShortRange::ForEwaldReal){
    looppairset<ShortRange::LJEwaldReal,
                LJCoulombPotential,LJCoulombPotential,TPJ>
      (particlesetj, ghosttyperange, ghostpairindex, forcesetj, energy, cutoff2);
    looppairset<ShortRange::EwaldReal,
                LJCoulombPotential,OnlyCoulombPotential,TPJ>
      (particlesetj, ghosttyperange, ghostpairindex, forcesetj, energy, cutoff2);
    looppairset<ShortRange::EwaldReal,
                OnlyCoulombPotential,AnyCoulombPotential,TPJ>
      (particlesetj, ghosttyperange, ghostpairindex, forcesetj, energy, cutoff2);
  }
}

template<ShortRange::InteractionType T, int POTENTIALMODEL, typename TPI, typename TPJ>
static void loopjset(const TPI& particlei, const std::vector<TPJ>& particlesetj,
                     const std::vector<TypeRange>& ghosttyperange, 
                     const std::vector<int>& ghosttargetindex,
                     Force& forcei, ForceArray& forcesetj, double& energy, 
                     const double cutoff2)
{
  for(size_t t=0;t<ghosttargetindex.size();t++){
    int ti = ghosttargetindex[t];
    int j_begin = getbegin<POTENTIALMODEL>(ghosttyperange[ti]);
    int j_end = getend<POTENTIALMODEL>(ghosttyperange[ti]);
    Force force(0.0,0.0,0.0);
    double energyi=0.0;
    loopj<T>(particlei, particlesetj, j_begin, j_end,
             force, forcesetj, energyi, cutoff2);
    forcei += force;
    energy += energyi;
  }
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
template<ShortRange::InteractionType T, int POTENTIALMODEL, typename TPI, typename TPJ>
void
ShortRangeInteractionSet<T,POTENTIALMODEL,TPI,TPJ>::loopself_inset(
    const std::vector<TPI>& particleseti, 
    const std::vector<TypeRange>& typerange, 
    ForceArray& forceseti, double &energy) 
{
  for(size_t si=0;si<typerange.size();si++){
    ParticleRange r = getrange<POTENTIALMODEL>(typerange[si]);
    for(long i=r.begin;i<r.end;i++){
      Force force(0.0,0.0,0.0);
      double energyi=0.0;
      loopj<T>(particleseti[i], particleseti, i+1, r.end,
               force, forceseti, energyi);
      forceseti[i] += force;
      energy += energyi;
    }
  }
}

template<ShortRange::InteractionType T, int POTENTIALMODEL, typename TPI, typename TPJ>
void
ShortRangeInteractionSet<T,POTENTIALMODEL,TPI,TPJ>::loopself_inset(
    const std::vector<TPI>& particleseti, 
    const std::vector<TypeRange>& typerange, 
    ForceArray& forceseti, double &energy, 
    const double cutoff2) 
{
  for(size_t si=0;si<typerange.size();si++){
    ParticleRange r = getrange<POTENTIALMODEL>(typerange[si]);
    for(long i=r.begin;i<r.end;i++){
      Force force(0.0,0.0,0.0);
      double energyi=0.0;
      loopj<T>(particleseti[i], particleseti, i+1, r.end,
               force, forceseti, energyi, cutoff2);
      forceseti[i] += force;
      energy += energyi;
    }
  }
}

  //! self interaction for i-particle sets
template<ShortRange::InteractionType T, int POTENTIALMODEL, typename TPI, typename TPJ>
void
ShortRangeInteractionSet<T,POTENTIALMODEL,TPI,TPJ>::loopself(
    const std::vector<TPI>& particleseti, 
    const std::vector<TypeRange>& typerange, 
    const std::vector< std::vector<int> >& targetindex,
    ForceArray& forceseti, double &energy) 
{
  loopself_inset(particleseti, typerange, forceseti, energy);
  for(size_t si=0;si<typerange.size();si++){
    ParticleRange r = getrange<POTENTIALMODEL>(typerange[si]);
    for(size_t t=0;t<targetindex[si].size();t++){
      int ti = targetindex[si][t];
        //      if(ti>si)
      {
        ParticleRange tr = getrange<POTENTIALMODEL>(typerange[ti]);
        for(long i=r.begin;i<r.end;i++){
          Force force(0.0,0.0,0.0);
          double energyi=0.0;
          loopj<T>(particleseti[i], particleseti, tr.begin, tr.end,
                   force, forceseti, energyi);
          forceseti[i] += force;
          energy += energyi;
        }
      }
    }
  }
}

template<ShortRange::InteractionType T, int POTENTIALMODEL, typename TPI, typename TPJ>
void
ShortRangeInteractionSet<T,POTENTIALMODEL,TPI,TPJ>::loopself(
    const std::vector<TPI>& particleseti, 
    const std::vector<TypeRange>& typerange, 
    const std::vector< std::vector<int> >& targetindex,
    ForceArray& forceseti, double &energy, const double cutoff2)
{
  loopself_inset(particleseti, typerange, forceseti, energy);
  for(size_t si=0;si<typerange.size();si++){
    ParticleRange r = getrange<POTENTIALMODEL>(typerange[si]);
    for(size_t t=0;t<targetindex[si].size();t++){
      int ti = targetindex[si][t];
        //      if(ti>si)
      {
        ParticleRange tr = getrange<POTENTIALMODEL>(typerange[ti]);
        for(long i=r.begin;i<r.end;i++){
          Force force(0.0,0.0,0.0);
          double energyi=0.0;
          loopj<T>(particleseti[i], particleseti, tr.begin, tr.end,
                   force, forceseti, energyi, cutoff2);
          forceseti[i] += force;
          energy += energyi;
        }
      }
    }
  }
}

  //! interaction between i-particle sets and j-particle sets
template<ShortRange::InteractionType T, int POTENTIALMODEL, typename TPI, typename TPJ>
void
ShortRangeInteractionSet<T,POTENTIALMODEL,TPI,TPJ>::loopisetjsetij(
    const std::vector<TPI>& particleseti, 
    const std::vector<TypeRange>& typerange, 
    const std::vector<TPJ>& particlesetj,
    const std::vector<TypeRange>& ghosttyperange, 
    const std::vector< std::vector<int> >& ghosttargetindex,
    ForceArray& forceseti, ForceArray& forcesetj, double &energy)
{
  for(size_t si=0;si<typerange.size();si++){
    ParticleRange r = getrange<POTENTIALMODEL>(typerange[si]);
    for(size_t t=0;t<ghosttargetindex[si].size();t++){
      int ti = ghosttargetindex[si][t];
      ParticleRange tr = getrange<POTENTIALMODEL>(ghosttyperange[ti]);
      loopij<T>(particleseti, r.begin, r.end,
                particlesetj, tr.begin, tr.end, 
                forceseti, forcesetj, energy);
    }
  }
}

template<ShortRange::InteractionType T, int POTENTIALMODEL, typename TPI, typename TPJ>
void
ShortRangeInteractionSet<T,POTENTIALMODEL,TPI,TPJ>::loopisetjsetij(
    const std::vector<TPI>& particleseti, 
    const std::vector<TypeRange>& typerange, 
    const std::vector<TPJ>& particlesetj,
    const std::vector<TypeRange>& ghosttyperange, 
    const std::vector< std::vector<int> >& ghosttargetindex,
    ForceArray& forceseti, ForceArray& forcesetj, double &energy, 
    const double cutoff2)
{
  for(size_t si=0;si<typerange.size();si++){
    ParticleRange r = getrange<POTENTIALMODEL>(typerange[si]);
    for(size_t t=0;t<ghosttargetindex[si].size();t++){
      int ti = ghosttargetindex[si][t];
      ParticleRange tr = getrange<POTENTIALMODEL>(ghosttyperange[ti]);
      loopij<T>(particleseti, r.begin, r.end,
                particlesetj, tr.begin, tr.end, 
                forceseti, forcesetj, energy,
                cutoff2);
    }
  }
}

template<ShortRange::InteractionType T, int POTENTIALMODEL, typename TPI, typename TPJ>
void
ShortRangeInteractionSet<T,POTENTIALMODEL,TPI,TPJ>::loopijset(
    const std::vector<TPI>& particleseti, 
    const std::vector<TypeRange>& typerange, 
    const std::vector<TPJ>& particlesetj,
    const std::vector<TypeRange>& ghosttyperange, 
    const std::vector< std::vector<int> >& ghosttargetindex,
    ForceArray& forceseti, ForceArray& forcesetj, double &energy)
{
  for(size_t si=0;si<typerange.size();si++){
    ParticleRange r = getrange<POTENTIALMODEL>(typerange[si]);
    for(int i=r.begin;i<r.end;i++){
      loopjset<T,POTENTIALMODEL,TPI,TPJ>(particleseti[i], 
                                         particlesetj, ghosttyperange, ghosttargetindex[si], 
                                         forceseti[i], forcesetj, energy);
    }
  }
}

template<ShortRange::InteractionType T, int POTENTIALMODEL, typename TPI, typename TPJ>
void
ShortRangeInteractionSet<T,POTENTIALMODEL,TPI,TPJ>::loopijset(
    const std::vector<TPI>& particleseti, 
    const std::vector<TypeRange>& typerange, 
    const std::vector<TPJ>& particlesetj,
    const std::vector<TypeRange>& ghosttyperange, 
    const std::vector< std::vector<int> >& ghosttargetindex,
    ForceArray& forceseti, ForceArray& forcesetj, double &energy, 
    const double cutoff2) 
{
  for(size_t si=0;si<typerange.size();si++){
    ParticleRange r = getrange<POTENTIALMODEL>(typerange[si]);
    for(int i=r.begin;i<r.end;i++){
      loopjset<T,POTENTIALMODEL,TPI,TPJ>(particleseti[i], 
                                         particlesetj, ghosttyperange, ghosttargetindex[si], 
                                         forceseti[i], forcesetj, energy, cutoff2);
    }
  }
}

template<ShortRange::InteractionType T, int POTENTIALMODEL, typename TPI, typename TPJ>
void
ShortRangeInteractionSet<T,POTENTIALMODEL,TPI,TPJ>::allloop(
    const std::vector<TPI>& particleseti,
    const std::vector<TypeRange>& typerange, 
    const std::vector< std::vector<int> >& targetindex,
    const std::vector<TPJ>& particlesetj,
    const std::vector<TypeRange>& ghosttyperange, 
    const std::vector< std::vector<int> >& ghosttargetindex,
               const std::vector< std::pair<int,int> >& ghostpairindex,
    ForceArray& forceseti,                                   
    ForceArray& forcesetj,                 
               double &energyself,
    double &energy)
{
  for(long si=0;si<long(typerange.size());si++){
    ParticleRange r = getrange<POTENTIALMODEL>(typerange[si]);
    // i-T   i-T,j-T
    for(int i=r.begin;i<r.end;i++){
      loopj<T>( particleseti[i], particleseti, i+1, r.end,
                forceseti[i], forceseti, energyself);
      for(long t=0;t<long(targetindex[si].size());t++){
        int ti = targetindex[si][t];
        //        if(ti>si)
        {
          ParticleRange tr = getrange<POTENTIALMODEL>(typerange[ti]);
          loopj<T>(particleseti[i], particleseti, tr.begin, tr.end,
                   forceseti[i], forceseti, energyself);
          
        }
      }
      loopjset<T,POTENTIALMODEL,TPI,TPJ>(particleseti[i], 
                                         particlesetj, 
                                         ghosttyperange, 
                                         ghosttargetindex[si], 
                                         forceseti[i], 
                                         forcesetj,
                                         energy);   
    }
  }
  looppairset<T,POTENTIALMODEL,POTENTIALMODEL,TPJ>(particlesetj,
                                                   ghosttyperange,
                                                   ghostpairindex,
                                                   forcesetj,
                                                   energy);
}

template<ShortRange::InteractionType T, int POTENTIALMODEL, typename TPI, typename TPJ>
void
ShortRangeInteractionSet<T,POTENTIALMODEL,TPI,TPJ>::allloop(
    const std::vector<TPI>& particleseti,
    const std::vector<TypeRange>& typerange, 
    const std::vector< std::vector<int> >& targetindex,
    const std::vector<TPJ>& particlesetj,
    const std::vector<TypeRange>& ghosttyperange, 
    const std::vector< std::vector<int> >& ghosttargetindex,
                                  const std::vector< std::pair<int,int> >& ghostpairindex,
    ForceArray& forceseti,                                   
    ForceArray& forcesetj,                 
                                  double &energyself,
                                  double &energy, const double cutoff2)
{
  for(long si=0;si<long(typerange.size());si++){
    ParticleRange r = getrange<POTENTIALMODEL>(typerange[si]);
      // i-T   i-T,j-T
    for(int i=r.begin;i<r.end;i++){
      loopj<T>( particleseti[i], particleseti, i+1, r.end,
                forceseti[i], forceseti, energyself, cutoff2);
      for(long t=0;t<long(targetindex[si].size());t++){
        int ti = targetindex[si][t];
          //      if(ti>si)
        {
          ParticleRange tr = getrange<POTENTIALMODEL>(typerange[ti]);
          loopj<T>(particleseti[i], particleseti, tr.begin, tr.end,
                   forceseti[i], forceseti, energyself, cutoff2);

        }
      }
      loopjset<T,POTENTIALMODEL,TPI,TPJ>(particleseti[i], 
                                         particlesetj, 
                                         ghosttyperange, 
                                         ghosttargetindex[si], 
                                         forceseti[i], forcesetj, energy,
                                         cutoff2);   
    }
  }
  looppairset<T,POTENTIALMODEL,POTENTIALMODEL,TPJ>(particlesetj,
                                                   ghosttyperange,
                                                   ghostpairindex,
                                                   forcesetj,
                                                   energy, cutoff2);
}


/* 
   if(CT==ShortRange::OriginalCoulomb)   if(CT==ShortRange::ForEwaldReal) 
   TLJC ShortRange::LJCoulomb            ShortRange::LJEwaldReal
   TCL ShortRange::Coulomb               ShortRange::EwaldReal
*/
template<ShortRange::InteractionType TLJC, ShortRange::InteractionType TCL, 
        typename TPI, typename TPJ>
void withcoulombloop(const std::vector<TPI>& particleseti,
                     const std::vector<TypeRange>& typerange, 
                     const std::vector<int>& targetindex,
                     ParticleRange ljcr, ParticleRange clr,
                     const std::vector<TPJ>& particlesetj,
                     const std::vector<TypeRange>& ghosttyperange, 
                     const std::vector<int>& ghosttargetindex, 
                     ForceArray& forceseti,
                     ForceArray& forcesetj,
                     double &energyself,
                     double &energy)
{
  for(int i=ljcr.begin;i<ljcr.end;i++){
    // vs self LJC
    loopj<TLJC>(particleseti[i], 
                particleseti, i+1, ljcr.end,
                forceseti[i], forceseti, energyself);
    // vs i-LJC
    for(size_t t=0;t<targetindex.size();t++){
      int ti = targetindex[t];
      //        if(ti>si)
      {
        ParticleRange tr = getrange<LJCoulombPotential>(typerange[ti]);
        loopj<TLJC>(particleseti[i], 
                    particleseti, tr.begin, tr.end,
                    forceseti[i], forceseti, energyself);
      }
    }         
    // vs self CL
    loopj<TCL>(particleseti[i], 
               particleseti, clr.begin, clr.end,
               forceseti[i], forceseti, energyself);
    // vs i-CL
    for(size_t t=0;t<targetindex.size();t++){
      int ti = targetindex[t];
      //        if(ti>si)
      {
        ParticleRange tr = getrange<OnlyCoulombPotential>(typerange[ti]);
        loopj<TCL>(particleseti[i], 
                   particleseti, tr.begin, tr.end,
                   forceseti[i], forceseti, energyself);
      }
    }         
    // vs j-LJ
    loopjset<ShortRange::LJ,OnlyLJPotential,TPI,TPJ>(particleseti[i], 
                                                     particlesetj, 
                                                     ghosttyperange, 
                                                     ghosttargetindex, 
                                                     forceseti[i],
                                                     forcesetj,
                                                     energy);
    // vs j-LJC
    loopjset<TLJC,LJCoulombPotential,TPI,TPJ>(particleseti[i],
                                              particlesetj, 
                                              ghosttyperange, 
                                              ghosttargetindex, 
                                              forceseti[i],
                                              forcesetj,
                                              energy);
    // vs j-CL
    loopjset<TCL,OnlyCoulombPotential,TPI,TPJ>(particleseti[i], 
                                               particlesetj, 
                                               ghosttyperange, 
                                               ghosttargetindex, 
                                               forceseti[i],
                                               forcesetj,
                                               energy);
  }
}

/*
  if(CT==ShortRange::OriginalCoulomb)   if(CT==ShortRange::ForEwaldReal) 
  TCL ShortRange::Coulomb               ShortRange::EwaldReal
*/
template<ShortRange::InteractionType TCL, 
         typename TPI, typename TPJ>
void onlycoulombloop(const std::vector<TPI>& particleseti,
                     const std::vector<TypeRange>& typerange, 
                     const std::vector<int>& targetindex,
                     ParticleRange clr,
                     const std::vector<TPJ>& particlesetj,
                     const std::vector<TypeRange>& ghosttyperange, 
                     const std::vector<int>& ghosttargetindex, 
                     ForceArray& forceseti,
                     ForceArray& forcesetj,
                     double &energyself,
                     double &energy)
{
  for(int i=clr.begin;i<clr.end;i++){
    loopj<TCL>(particleseti[i], particleseti, i+1, clr.end,
               forceseti[i], forceseti, energy);
    for(size_t t=0;t<targetindex.size();t++){
      int ti = targetindex[t];
      //        if(ti>si)
      {
        int j_begin = getbegin<AnyCoulombPotential>(typerange[ti]);
        int j_end = getend<AnyCoulombPotential>(typerange[ti]);
        loopj<TCL>(particleseti[i], particleseti, j_begin, j_end,
                   forceseti[i], forceseti, energyself);
      }
    }
    loopjset<TCL,AnyCoulombPotential,TPI,TPJ>(particleseti[i], 
                                              particlesetj, 
                                              ghosttyperange, 
                                              ghosttargetindex, 
                                              forceseti[i],
                                              forcesetj,
                                              energy);
  }
}

template<ShortRange::CoulombType CT, typename TPI, typename TPJ>
void 
AllTypeShortRangeInteractionSet<CT,TPI,TPJ>::allloop(
    const std::vector<TPI>& particleseti,
    const std::vector<TypeRange>& typerange, 
    const std::vector< std::vector<int> >& targetindex,
    const std::vector<TPJ>& particlesetj,
    const std::vector<TypeRange>& ghosttyperange, 
    const std::vector< std::vector<int> >& ghosttargetindex,
                                                     const std::vector< std::pair<int,int> >& ghostpairindex,
    ForceArray& forceseti,
    ForceArray& forcesetj,
                                                     double &energyself,
                                                     double &energy)
{
  for(long si=0;si<long(typerange.size());si++){
    ParticleRange ljr = getrange<OnlyLJPotential>(typerange[si]);
    ParticleRange ljcr = getrange<LJCoulombPotential>(typerange[si]);
    ParticleRange clr = getrange<OnlyCoulombPotential>(typerange[si]);
#ifdef RDTSCTIME
    starttime = rdtsc();
#endif
    // i-LJ vs  i-LJ,LJC,j-LJ,LJC
    for(int i=ljr.begin;i<ljr.end;i++){
      // vs self LJ, self LJC
      loopj<ShortRange::LJ>(particleseti[i], particleseti, i+1, ljcr.end,
                            forceseti[i], forceseti, energyself);
      // vs i-LJ, i-LJC
      for(size_t t=0;t<targetindex[si].size();t++){
        int ti = targetindex[si][t];
        //        if(ti>si)
        {
          int j_begin = getbegin<AnyLJPotential>(typerange[ti]);
          int j_end = getend<AnyLJPotential>(typerange[ti]);
          loopj<ShortRange::LJ>(particleseti[i], particleseti, j_begin, j_end,
                                forceseti[i], forceseti, energyself);
        }
      }
      // vs j-LJ, j-LJC
      loopjset<ShortRange::LJ,AnyLJPotential,TPI,TPJ>(particleseti[i], 
                                                      particlesetj, 
                                                      ghosttyperange, 
                                                      ghosttargetindex[si], 
                                                      forceseti[i],
                                                      forcesetj,
                                                      energy);
    }

#ifdef RDTSCTIME
    stoptime = rdtsc();
    ljtime += lap_rdtsc(starttime,stoptime);
    starttime = rdtsc();
#endif
    // i-LJC vs i-LJC,CL,j-LJC,CL
    if(CT==ShortRange::OriginalCoulomb){
      withcoulombloop<ShortRange::LJCoulomb,ShortRange::Coulomb>(
          particleseti, typerange, targetindex[si], ljcr, clr,
          particlesetj, ghosttyperange, ghosttargetindex[si],
          forceseti, forcesetj, energyself, energy);
    }else if(CT==ShortRange::ForEwaldReal){
      withcoulombloop<ShortRange::LJEwaldReal,ShortRange::EwaldReal>(
          particleseti, typerange, targetindex[si], ljcr, clr,
          particlesetj, ghosttyperange, ghosttargetindex[si],
          forceseti, forcesetj, energyself, energy);
    }else{
      std::cout << "Unknown CoulombType" << std::endl;
    }
#ifdef RDTSCTIME
    stoptime = rdtsc();
    ljctime += lap_rdtsc(starttime,stoptime);
    starttime = rdtsc();
#endif
    // i-CL  i-CL,j-LJC,CL
    if(CT==ShortRange::OriginalCoulomb){
      onlycoulombloop<ShortRange::Coulomb>(
          particleseti, typerange, targetindex[si], clr,
          particlesetj, ghosttyperange, ghosttargetindex[si],
          forceseti, forcesetj, energyself, energy );
    }else if(CT==ShortRange::ForEwaldReal){
      onlycoulombloop<ShortRange::EwaldReal>(
          particleseti, typerange, targetindex[si], clr,
          particlesetj, ghosttyperange, ghosttargetindex[si],
          forceseti, forcesetj, energyself, energy );
    }
#ifdef RDTSCTIME
    stoptime = rdtsc();
    cltime += lap_rdtsc(starttime,stoptime);
#endif
  }
  looppairset<CT,TPJ>(particlesetj, ghosttyperange, ghostpairindex, forcesetj, energy);
}

////////////// without cutoff


////////////// with cutoff

/* 
   if(CT==ShortRange::OriginalCoulomb)   if(CT==ShortRange::ForEwaldReal) 
   TLJC ShortRange::LJCoulomb            ShortRange::LJEwaldReal
   TCL ShortRange::Coulomb               ShortRange::EwaldReal
*/
template<ShortRange::InteractionType TLJC, ShortRange::InteractionType TCL, 
        typename TPI, typename TPJ>
void withcoulombloop(const std::vector<TPI>& particleseti,
                     const std::vector<TypeRange>& typerange, 
                     const std::vector<int>& targetindex,
                     ParticleRange ljcr, ParticleRange clr,
                     const std::vector<TPJ>& particlesetj,
                     const std::vector<TypeRange>& ghosttyperange, 
                     const std::vector<int>& ghosttargetindex, 
                     ForceArray& forceseti,
                     ForceArray& forcesetj,
                     double &energyself,
                     double &energy,
                     double cutoff2)
{
  for(int i=ljcr.begin;i<ljcr.end;i++){
    // vs self LJC
    loopj<TLJC>(particleseti[i], 
                particleseti, i+1, ljcr.end,
                forceseti[i], forceseti, energyself, cutoff2);
    // vs i-LJC
    for(size_t t=0;t<targetindex.size();t++){
      int ti = targetindex[t];
      //        if(ti>si)
      {
        ParticleRange tr = getrange<LJCoulombPotential>(typerange[ti]);
        loopj<TLJC>(particleseti[i], 
                    particleseti, tr.begin, tr.end,
                    forceseti[i], forceseti, energyself, cutoff2);
      }
    }         
    // vs self CL
    loopj<TCL>(particleseti[i], 
               particleseti, clr.begin, clr.end,
               forceseti[i], forceseti, energyself, cutoff2);
    // vs i-CL
    for(size_t t=0;t<targetindex.size();t++){
      int ti = targetindex[t];
      //        if(ti>si)
      {
        ParticleRange tr = getrange<OnlyCoulombPotential>(typerange[ti]);
        loopj<TCL>(particleseti[i], 
                   particleseti, tr.begin, tr.end,
                   forceseti[i], forceseti, energyself, cutoff2);
      }
    }         
    // vs j-LJ
    loopjset<ShortRange::LJ,OnlyLJPotential,TPI,TPJ>(particleseti[i], 
                                                     particlesetj, 
                                                     ghosttyperange, 
                                                     ghosttargetindex, 
                                                     forceseti[i],
                                                     forcesetj,
                                                     energy, cutoff2);
    // vs j-LJC
    loopjset<TLJC,LJCoulombPotential,TPI,TPJ>(particleseti[i],
                                              particlesetj, 
                                              ghosttyperange, 
                                              ghosttargetindex, 
                                              forceseti[i],
                                              forcesetj,
                                              energy, cutoff2);
    // vs j-CL
    loopjset<TCL,OnlyCoulombPotential,TPI,TPJ>(particleseti[i], 
                                               particlesetj, 
                                               ghosttyperange, 
                                               ghosttargetindex, 
                                               forceseti[i],
                                               forcesetj,
                                               energy, cutoff2);
  }
}

/*
  if(CT==ShortRange::OriginalCoulomb)   if(CT==ShortRange::ForEwaldReal) 
  TCL ShortRange::Coulomb               ShortRange::EwaldReal
*/
template<ShortRange::InteractionType TCL, 
         typename TPI, typename TPJ>
void onlycoulombloop(const std::vector<TPI>& particleseti,
                     const std::vector<TypeRange>& typerange, 
                     const std::vector<int>& targetindex,
                     ParticleRange clr,
                     const std::vector<TPJ>& particlesetj,
                     const std::vector<TypeRange>& ghosttyperange, 
                     const std::vector<int>& ghosttargetindex, 
                     ForceArray& forceseti,
                     ForceArray& forcesetj,
                     double &energyself,
                     double &energy, double cutoff2)
{
  for(int i=clr.begin;i<clr.end;i++){
    loopj<TCL>(particleseti[i], particleseti, i+1, clr.end,
               forceseti[i], forceseti, energyself, cutoff2);
    for(size_t t=0;t<targetindex.size();t++){
      int ti = targetindex[t];
      //        if(ti>si)
      {
        int j_begin = getbegin<AnyCoulombPotential>(typerange[ti]);
        int j_end = getend<AnyCoulombPotential>(typerange[ti]);
        loopj<TCL>(particleseti[i], particleseti, j_begin, j_end,
                   forceseti[i], forceseti, energyself, cutoff2);
      }
    }
    loopjset<TCL,AnyCoulombPotential,TPI,TPJ>(particleseti[i], 
                                              particlesetj, 
                                              ghosttyperange, 
                                              ghosttargetindex, 
                                              forceseti[i],
                                              forcesetj,
                                              energy, cutoff2);
  }
}

template<ShortRange::CoulombType CT, typename TPI, typename TPJ>
void 
AllTypeShortRangeInteractionSet<CT,TPI,TPJ>::allloop(
    const std::vector<TPI>& particleseti,
    const std::vector<TypeRange>& typerange, 
    const std::vector< std::vector<int> >& targetindex,
    const std::vector<TPJ>& particlesetj,
    const std::vector<TypeRange>& ghosttyperange, 
    const std::vector< std::vector<int> >& ghosttargetindex,
    const std::vector< std::pair<int,int> >& ghostpairindex,
    ForceArray& forceseti,
    ForceArray& forcesetj,
    double &energyself,
    double &energy, 
    const double cutoff2)
{
  for(long si=0;si<long(typerange.size());si++){
    ParticleRange ljr = getrange<OnlyLJPotential>(typerange[si]);
    ParticleRange ljcr = getrange<LJCoulombPotential>(typerange[si]);
    ParticleRange clr = getrange<OnlyCoulombPotential>(typerange[si]);
#ifdef RDTSCTIME
    starttime = rdtsc();
#endif
    // i-LJ   i-LJ,LJC,j-LJ,LJC
    for(int i=ljr.begin;i<ljr.end;i++){
      // vs self LJ, self LJC
      loopj<ShortRange::LJ>(particleseti[i], particleseti, i+1, ljcr.end,
                            forceseti[i], forceseti, energyself, cutoff2);
      // vs i-LJ, i-LJC
      for(size_t t=0;t<targetindex[si].size();t++){
        int ti = targetindex[si][t];
        //        if(ti>si)
        {
          int j_begin = getbegin<AnyLJPotential>(typerange[ti]);
          int j_end = getend<AnyLJPotential>(typerange[ti]);
          loopj<ShortRange::LJ>(particleseti[i], particleseti, j_begin, j_end,
                                forceseti[i], forceseti, energyself, cutoff2);
        }
      }
      // vs j-LJ, j-LJC
      loopjset<ShortRange::LJ,AnyLJPotential,TPI,TPJ>(particleseti[i], 
                                                      particlesetj, 
                                                      ghosttyperange, 
                                                      ghosttargetindex[si], 
                                                      forceseti[i],
                                                      forcesetj,
                                                      energy, cutoff2);
    }
#ifdef RDTSCTIME
    stoptime = rdtsc();
    ljtime += lap_rdtsc(starttime,stoptime);
    starttime = rdtsc();
#endif
    // i-LJC i-LJC,CL,j-LJC,CL
    if(CT==ShortRange::OriginalCoulomb){
      withcoulombloop<ShortRange::LJCoulomb,ShortRange::Coulomb>(
          particleseti, typerange, targetindex[si], ljcr, clr,
          particlesetj, ghosttyperange, ghosttargetindex[si],
          forceseti, forcesetj, energyself, energy,
          cutoff2);
    }else if(CT==ShortRange::ForEwaldReal){
      withcoulombloop<ShortRange::LJEwaldReal,ShortRange::EwaldReal>(
          particleseti, typerange, targetindex[si], ljcr, clr,
          particlesetj, ghosttyperange, ghosttargetindex[si],
          forceseti, forcesetj, energyself, energy,
          cutoff2);
    }else{
      std::cout << "Unknown CoulombType" << std::endl;
    }
#ifdef RDTSCTIME
    stoptime = rdtsc();
    ljctime += lap_rdtsc(starttime,stoptime);
    starttime = rdtsc();
#endif
    // i-CL  i-CL,j-LJC,CL
    if(CT==ShortRange::OriginalCoulomb){
      onlycoulombloop<ShortRange::Coulomb>(
          particleseti, typerange, targetindex[si], clr,
          particlesetj, ghosttyperange, ghosttargetindex[si],
          forceseti, forcesetj, energyself, energy,
          cutoff2);
    }else if(CT==ShortRange::ForEwaldReal){
      onlycoulombloop<ShortRange::EwaldReal>(
          particleseti, typerange, targetindex[si], clr,
          particlesetj, ghosttyperange, ghosttargetindex[si],
          forceseti, forcesetj, energyself, energy,
          cutoff2);
    }

#ifdef RDTSCTIME
    stoptime = rdtsc();
    cltime += lap_rdtsc(starttime,stoptime);
#endif
  }
  looppairset<CT,TPJ>(particlesetj, ghosttyperange, ghostpairindex,
                      forcesetj, energy, cutoff2);
}




template class ShortRangeInteractionSet<ShortRange::LJ, AnyLJPotential, Particle, Particle>;
template class AllTypeShortRangeInteractionSet<ShortRange::OriginalCoulomb, Particle, Particle>;
template class AllTypeShortRangeInteractionSet<ShortRange::ForEwaldReal, Particle, Particle>;

