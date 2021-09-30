#ifndef COMMON_H
#define COMMON_H

#ifdef OVERLAP
#include <mpi.h>
#endif
#include <vector>
#include <map>
#include <cstddef>
#include "SpaceVector.h"

namespace DebugLog {
  extern int verbose;
  extern int particle_dump;
}

typedef SpaceVector<double> Position;
typedef SpaceVector<double> Velocity;
typedef SpaceVector<double> Force;
typedef SpaceVector<double> Acceleration;


#define NO_SHORT_ID -1
#define NO_LONG_ID -1


enum CellIndexUpdateMode{
  UpdateFalse,
  UpdateTrue,
  UpdateAuto,
  UpdateAccurate
};

typedef std::vector<int> PairList;

enum PotentialModel{
  NoPotential,
  OnlyLJPotential,
  OnlyCoulombPotential,
  LJCoulombPotential,
  AnyLJPotential,
  AnyCoulombPotential,
  AnyPotential,
  UnkonwnPotential,
};

/*
enum RangeType{
  AllType,
  LJType,
  LJCoulombType,
  LJAllType,
  CoulombType,
  CoulombAllType,
  UnkonwnType,
};
*/

//! index range of particle subset
/*! begin >= end means empty subset
 */
struct ParticleRange{
  int begin;
  int end;
  //  ParticleRange() : begin(0), end(0) {}

  void shift(int offset)
  {
    begin += offset;
    end += offset;
  }
  int number_of_particle()
  {
    return end-begin;
  }
};

//! index ranges of particle subsets divided by force type
/*!
  TypeRange typerange;

  [typerange.being, typerange.end)                     : all particle
  [typerange.lj.begin, typerange.lj.end)               : lj particle
  [typerange.ljcoulomb.begin, typerange.ljcoulomb.end) : ljcoulomb particle
  [typerange.coulomb.begin, typerange.coulomb.end)     : coulomb particle
 */
struct TypeRange : public ParticleRange {
  ParticleRange lj;
  ParticleRange ljcoulomb;
  ParticleRange coulomb;

  void shift(int offset)
  {
    /*
    lj.begin += offset;
    lj.end += offset;
    ljcoulomb.begin += offset;
    ljcoulomb.end += offset;
    coulomb.begin += offset;
    coulomb.end += offset;
    all.begin += offset;
    all.end += offset;
    */
    ParticleRange::shift(offset);
    lj.shift(offset);
    ljcoulomb.shift(offset);
    coulomb.shift(offset);
  }

  PotentialModel particle_index_to_type(int i)
  {
    PotentialModel t;
    if((i>=lj.begin)&&(i<lj.end)){
      t = OnlyLJPotential;
    }else if((i>=ljcoulomb.begin)&&(i<ljcoulomb.end)){
      t = LJCoulombPotential;
    }else if((i>=coulomb.begin)&&(i<coulomb.end)){
      t = OnlyCoulombPotential;
    }else{
      t = UnkonwnPotential;
      printf("index %d not in this range\n", i);
    }
    return t;
  }
};

template<int T>
inline ParticleRange getrange(const TypeRange& tr)
{
  if(T==OnlyCoulombPotential){
    return tr.coulomb;
  }else if(T==LJCoulombPotential){
    return tr.ljcoulomb;
  }else if(T==OnlyLJPotential){
    return tr.lj;
  }else if(T==AnyLJPotential){
    TypeRange alj;
    alj.begin = tr.lj.begin;
    alj.end = tr.ljcoulomb.end;
    return alj;
  }else{
    return tr;
  }
}

/* 
   getbeginhave<T>,  getendhave<T>
   return begin/end of particle range which interact specified by T
   for LJ         lj and ljcoulomb
   for LJCoulomb  ljcoulomb
   for Coulomb    ljcoulomb and coulomb
   range must be continuos
    lj.end = ljcoulomb.begin, ljcoulomb.end = coulomb.begin

   !!! use these function, ljcoulomb caluculate separately
     type      i-range                    j-range
     LJ        getbegin/endhave<LJ>       getbegin/endhave<LJ>
     Coulomb   getbegin/endhave<Coulomb>  getbegin/endhave<Coulomb>

   !!! To reduce calculation using ljcoulomb
      type      i-range              j-range
     LJ         getrange<LJ>         getbegin/endhave<LJ>
                getrange<LJCoulomb>  getrange<LJ>
     LJCoulomb  getrange<LJCoulomb>  getrange<LJCoulomb>
     Coulomb    getrange<LJCoulomb>  getrange<Coulomb>
                getrange<Coulomb>    getbegin/endhave<Coulomb> 
 */
//! return begin of range which particles has Type T force/potential
template<int T>
inline int getbegin(const TypeRange& tr)
{
  if(T==AnyCoulombPotential){
    return tr.ljcoulomb.begin;
  }else if(T==OnlyCoulombPotential){
    return tr.coulomb.begin;
  }else if(T==LJCoulombPotential){
    return tr.ljcoulomb.begin;
  }else if(T==AnyLJPotential){
    return tr.lj.begin;
  }else if(T==OnlyLJPotential){
    return tr.lj.begin;
  }else{
    return tr.begin;
  }
}
//! return begin of range which particles has Type T force/potential
template<int T>
inline int getend(const TypeRange& tr)
{
  if(T==AnyCoulombPotential){
    return tr.coulomb.end;
  }else if(T==OnlyCoulombPotential){
    return tr.coulomb.end;
  }else if(T==LJCoulombPotential){
    return tr.ljcoulomb.end;
  }else if(T==AnyLJPotential){
    return tr.ljcoulomb.end;
  }else if(T==OnlyLJPotential){
    return tr.lj.end;
  }else{
    return tr.end;
  }
}

inline
int number_of_particle(std::vector<TypeRange> targettyperange)
{
  int num=0;
  for(std::vector<TypeRange>::size_type i=0;i<targettyperange.size();i++){
    num+=targettyperange[i].number_of_particle();
  }
  return num;
}

inline
void make_idset_to_index(std::vector< std::vector<int> >& idset,
                         std::map<int,int>& idset_to_index)
{
  for(int i=0;i<int(idset.size());i++){
    for(int si=0;si<int(idset[i].size());si++){
      idset_to_index.insert(std::pair<int,int>(idset[i][si],i));
    }
  }
}

inline
void make_id_to_index(std::vector<int> setid,
                      std::map<int,int>& setid_to_index)
{
  for(size_t i=0;i<setid.size();i++){
    setid_to_index.insert(std::pair<int,int>(setid[i],i));
  }
}


// Periodic Condition

//! Shift position move over boundary to opsite side
template<class T>
inline int periodic_shift(T& p, T size)
{
  if(p<T(0)){
    p+=size;
    return +1;
  }else if(p>=size){
    p-=size;
    return -1;
  }
  return 0;
}

template<class T>
inline void periodic_shift(T& p, T size, int& flag)
{
  flag = 0;
  if(p<T(0)){
    p+=size;
    flag = -1;
  }else if(p>=size){
    p-=size;
    flag = 1;
  }
}

//! Shift relative position to (-size/2,+size/2]
/*!
  make |p| <= size/2
 */
template<class T>
inline void periodic_relative_shift(T& p, T size)
{
  if(p<=-size/T(2)){
    p+=size;
  }else if(p>size/T(2)){
    p-=size;
  }
}

//! check move out boundary
template<class T>
inline bool out_of_bounds(SpaceVector<T>p, SpaceVector<T>b)
{
  if( (p.x<T(0)) || (p.x>b.x)
      || (p.y<T(0)) || (p.y>b.y)
      || (p.z<T(0)) || (p.z>b.z) ){
    return true;
  }
  return false;
}

//! Required Shift type
enum ReqShift{
  PlusX,
  MinusX,
  PlusY,
  MinusY,
  PlusZ,
  MinusZ
};

//! Array of cells required periodic shift for target cell
struct ShiftCellArray{
  std::vector<int> plusx;
  std::vector<int> minusx;
  std::vector<int> plusy;
  std::vector<int> minusy;
  std::vector<int> plusz;
  std::vector<int> minusz;

  void clear(){
    plusx.clear();
    minusx.clear();
    plusy.clear();
    minusy.clear();
    plusz.clear();
    minusz.clear();
  }
};


inline unsigned long long rdtsc() 
{
#ifdef RDTSCTIME
  unsigned long long ret;
#ifdef __x86_64
  //  int counter=0;
  unsigned  low, high;
  __asm__ __volatile__ ("cpuid" : : : "eax", "ebx", "ecx", "edx");
  //  __asm__ __volatile__ ("rdtsc" : "=A" (ret));
  //  __asm__ __volatile__ ("rdpmc" : "=a" (low), "=d" (high) : "c" (counter));
  __asm__ __volatile__ ("rdtsc" : "=a" (low), "=d" (high));
  ret = (unsigned long long)(low)|((unsigned long long)(high)<<32);
#endif
  return ret;
#else
  return 0ULL;
#endif
}

inline unsigned long long lap_rdtsc(unsigned long long start,
                                    unsigned long long stop) 
{
  if(stop<=start){
    return 0xFFFFFFFFFFFFFFFFLL - start + stop;
  }else{
    return stop-start;
  }
}

#ifdef OVERLAP
namespace MPICHECKER {
  extern int count;
  extern MPI_Request* reqs;
  extern int index;
  extern int flag;
  extern MPI_Status* status;
  extern bool done;
  void mpi_checker();
  }
#endif

#ifdef FXPA
void hpc_start(int);
void hpc_stop(int);
#endif

#endif
