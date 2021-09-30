#ifndef PARTICLEINFO_H
#define PARTICLEINFO_H
#include <vector>
#include <cstdio>
#include <cstdlib>

#include "Common.h"

#define DEBUG 0

#define ATOMTYPE_WO 14
#define ATOMTYPE_WH 11   // H:0 HO:1 HS:2 HC:3 H1:4 H2:5 H3:6 HP:7 HA:8 H4:9 H5:10 WH:11
#define MAX_HBOND 4      // max no of H bond per one heavy atom

typedef int Atomtype;
typedef int AtomID;

//! all particle data in one struct
struct Particle{
  Position position;
  Velocity velocity; 
  Force force;
  double charge;
  double mass;
  double inv_mass;
  Atomtype atomtype;
  AtomID atomid;
  /*
  Particle() : position(0.0,0.0,0.0), 
               velocity(0.0,0.0,0.0), 
               force(0.0,0.0,0.0), 
               charge(0.0), mass(0.0), atomtype(0)
  {}
  Particle(const Particle *p) : position(p->position),
                                velocity(p->velocity),
                                force(p->force),
                                charge(p->charge),
                                mass(p->mass),
                                atomtype(p->atomtype)
  {}
  */
};


//! particle data, Array of Struct form
typedef std::vector<Particle> ParticleArray;

typedef struct {
  std::vector<Particle>::size_type begin;
  std::vector<Particle>::size_type end;
}Particleset;

//! particle double read data for force calculation
struct ParticlePosCharge{
  Position position;
  double charge;
};
typedef std::vector<ParticlePosCharge> PosChargeArray;

//! particle int read data for force calculation
typedef std::vector<Atomtype> AtomtypeArray;

//! particle read data for force calculation
struct ParticlePosChargeAT{
  Position position;
  double charge;
  Atomtype atomtype;
};
typedef std::vector<ParticlePosChargeAT> PosChargeATArray;

//! particle doubel write data for force calculation
typedef std::vector<Force> ForceArray;

//! particle int transfer data for MPIComm
typedef std::vector<AtomID> AtomIDArray;

//! particle used integration
struct ParticleParameter2{
  Velocity velocity; 
  double mass;
  double inv_mass;
};
typedef std::vector<ParticleParameter2> ParticleParameter2Array;

struct GhostParticleArray{
  PosChargeArray poscharge;
  AtomtypeArray atomtype;
  ForceArray force;
  AtomIDArray atomid;

  GhostParticleArray()
  : poscharge(), 
    atomtype(),
    force(), atomid()
  {
    poscharge.clear();
    atomtype.clear();
    force.clear();
    atomid.clear();
  }
  GhostParticleArray(const ParticleArray& pa)
  : poscharge(),
    atomtype(),
    force(), atomid()
  {
    size_t num = pa.size();
    poscharge.resize(num);
    atomtype.resize(num);
    force.resize(num);
    atomid.resize(num);
    for(size_t i=0;i<num;i++){
      poscharge[i].position = pa[i].position;
      /* poscharge[i].position.x = pa[i].position.x; */
      /* poscharge[i].position.y = pa[i].position.y; */
      /* poscharge[i].position.z = pa[i].position.z; */
      poscharge[i].charge = pa[i].charge;
      atomtype[i] = pa[i].atomtype;
      force[i] = pa[i].force;
      /* force[i].x = pa[i].force.x; */
      /* force[i].y = pa[i].force.y; */
      /* force[i].z = pa[i].force.z; */
      atomid[i] = pa[i].atomid;
    }
  }
  void clear()
  {
    poscharge.clear();
    atomtype.clear();
    force.clear();
    atomid.clear();
  }
  void resize(size_t num)
  {
    poscharge.resize(num);
    atomtype.resize(num);
    force.resize(num);
    atomid.resize(num);
  }
  void push_back(const Particle p)
  {
    ParticlePosCharge pc;
    pc.position = p.position;
    /* pc.position.x = p.position.x; */
    /* pc.position.y = p.position.y; */
    /* pc.position.z = p.position.z; */
    pc.charge = p.charge;
    atomtype.push_back(p.atomtype);
    poscharge.push_back(pc);
    force.push_back(p.force);
    atomid.push_back(p.atomid);
  }
  size_t size() const
  {
    return poscharge.size();
  }
};

struct CombinedParticleArray{
  typedef PosChargeArray::size_type size_type;

  PosChargeArray poscharge;
  AtomtypeArray atomtype;
  ForceArray force;
  AtomIDArray atomid;
  ParticleParameter2Array parameters;

  void clear()
  {
    poscharge.clear();
    atomtype.clear();
    force.clear();
    atomid.clear();
    parameters.clear();
  }
  
  void reserve(size_type n)
  {
    poscharge.reserve(n);
    atomtype.reserve(n);
    force.reserve(n);
    atomid.reserve(n);
    parameters.reserve(n);
  }

  void resize(size_type n)
  {
    poscharge.resize(n);
    atomtype.resize(n);
    force.resize(n);
    atomid.resize(n);
    parameters.resize(n);
  }

  void swap(CombinedParticleArray& from)
  {
    poscharge.swap(from.poscharge);
    atomtype.swap(from.atomtype);
    force.swap(from.force);
    atomid.swap(from.atomid);
    parameters.swap(from.parameters);
  }

  size_type begin(){
    return 0;
  }
  size_type end(){
    return poscharge.size();
  }

  CombinedParticleArray()
      : poscharge(), 
        atomtype(), 
        force(), atomid(), parameters()
  {
  }

  explicit CombinedParticleArray(int num)
    : poscharge(num), 
    atomtype(num), 
        force(num), atomid(num), parameters(num)
  {
  }

  CombinedParticleArray(const ParticleArray& pa)
  : poscharge(), 
    atomtype(), 
    force(), atomid(), parameters()
  {
    size_t num = pa.size();
    poscharge.resize(num);
    atomtype.resize(num);
    force.resize(num);
    atomid.resize(num);
    parameters.resize(num);
    for(size_t i=0;i<num;i++){
      poscharge[i].position = pa[i].position;
      /* poscharge[i].position.x = pa[i].position.x; */
      /* poscharge[i].position.y = pa[i].position.y; */
      /* poscharge[i].position.z = pa[i].position.z; */
      poscharge[i].charge = pa[i].charge;
      atomtype[i] = pa[i].atomtype;
      force[i] = pa[i].force;
      /* force[i].x = pa[i].force.x; */
      /* force[i].y = pa[i].force.y; */
      /* force[i].z = pa[i].force.z; */
      atomid[i] = pa[i].atomid;
      parameters[i].velocity = pa[i].velocity;
      /* parameters[i].velocity.x = pa[i].velocity.x; */
      /* parameters[i].velocity.y = pa[i].velocity.y; */
      /* parameters[i].velocity.z = pa[i].velocity.z; */
      parameters[i].mass = pa[i].mass;
      parameters[i].inv_mass = pa[i].inv_mass;
    }
  }
  void push_back(const Particle p)
  {
    ParticlePosCharge pc;
    pc.position = p.position;
    /* pc.position.x = p.position.x; */
    /* pc.position.y = p.position.y; */
    /* pc.position.z = p.position.z; */
    pc.charge = p.charge;
    atomtype.push_back(p.atomtype);
    poscharge.push_back(pc);
    force.push_back(p.force);
    atomid.push_back(p.atomid);
    ParticleParameter2 p2;
    p2.velocity = p.velocity;
    /* p2.velocity.x = p.velocity.x; */
    /* p2.velocity.y = p.velocity.y; */
    /* p2.velocity.z = p.velocity.z; */
    p2.mass = p.mass;
    p2.inv_mass = p.inv_mass;
    parameters.push_back(p2);
  }
  size_t size() const
  {
    return poscharge.size();
  }

  template<typename ID>
  inline Particle getparticle(const ID i) const
  {
    Particle p;
    p.position = poscharge[i].position;
    /* p.position.x = poscharge[i].position.x; */
    /* p.position.y = poscharge[i].position.y; */
    /* p.position.z = poscharge[i].position.z; */
    p.charge = poscharge[i].charge;
    p.force = force[i];
    /* p.force.x = force[i].x; */
    /* p.force.y = force[i].y; */
    /* p.force.z = force[i].z; */
    p.velocity = parameters[i].velocity;
    /* p.velocity.x = parameters[i].velocity.x; */
    /* p.velocity.y = parameters[i].velocity.y; */
    /* p.velocity.z = parameters[i].velocity.z; */
    p.mass = parameters[i].mass;
    p.inv_mass = parameters[i].inv_mass;
    p.atomtype = atomtype[i];
    p.atomid = atomid[i];
    return p;
  }
  template<typename ID>
  void setparticle(const Particle p, const ID i)
  {
    poscharge[i].position = p.position;
    /* poscharge[i].position.x = p.position.x; */
    /* poscharge[i].position.y = p.position.y; */
    /* poscharge[i].position.z = p.position.z; */
    poscharge[i].charge = p.charge;
    force[i] = p.force;
    /* force[i].x = p.force.x; */
    /* force[i].y = p.force.y; */
    /* force[i].z = p.force.z; */
    parameters[i].velocity = p.velocity;
    /* parameters[i].velocity.x = p.velocity.x; */
    /* parameters[i].velocity.y = p.velocity.y; */
    /* parameters[i].velocity.z = p.velocity.z; */
    parameters[i].mass = p.mass;
    parameters[i].inv_mass =  p.inv_mass;
    atomtype[i] = p.atomtype;
    atomid[i] =  p.atomid;
  }
};

struct PositionVelocity {
  Position position;
  Velocity velocity;
};

class PosVelArray : public std::vector<PositionVelocity>
{
 public:
  PosVelArray(int n)
  {
    resize(n);
  }
  PosVelArray(const ParticleArray& pa)
  {
    resize(pa.size());
    for(size_t i=0;i<size();i++){
      (*this)[i].position = pa[i].position;
      /* (*this)[i].position.x = pa[i].position.x; */
      /* (*this)[i].position.y = pa[i].position.y; */
      /* (*this)[i].position.z = pa[i].position.z; */
      (*this)[i].velocity = pa[i].velocity;
      /* (*this)[i].velocity.x = pa[i].velocity.x; */
      /* (*this)[i].velocity.y = pa[i].velocity.y; */
      /* (*this)[i].velocity.z = pa[i].velocity.z; */
    }
  }
  PosVelArray(const CombinedParticleArray& pa)
  {
    resize(pa.size());
    for(size_t i=0;i<size();i++){
      (*this)[i].position = pa.poscharge[i].position;
      /* (*this)[i].position.x = pa.poscharge[i].position.x; */
      /* (*this)[i].position.y = pa.poscharge[i].position.y; */
      /* (*this)[i].position.z = pa.poscharge[i].position.z; */
      (*this)[i].velocity = pa.parameters[i].velocity;
      /* (*this)[i].velocity.x = pa.parameters[i].velocity.x; */
      /* (*this)[i].velocity.y = pa.parameters[i].velocity.y; */
      /* (*this)[i].velocity.z = pa.parameters[i].velocity.z; */
    }
  }
};


//! particle data for force calculation
struct ParticleParameter1{
  Position position;
  Force force;
  double charge;
  Atomtype atomtype;
};

//! particle data for force calculation, Array of Struct form
typedef std::vector<ParticleParameter1> ParticleParameters1;

//! particle data, Struct of Array form
struct StructOfParticleDataArray{
  std::vector<Position> positions;
  std::vector<Velocity> velositys; 
  std::vector<Force> forces;
  std::vector<double> charges;
  std::vector<double> masses;
  std::vector<Atomtype> atomtypes;
};

template<typename ID>
inline const Particle getparticle(const CombinedParticleArray& pa, const ID i)
{
  return pa.getparticle(i);
}
template<typename ID>
inline Particle getparticle(CombinedParticleArray& pa, const ID i)
{
  return pa.getparticle(i);
}
template<typename ID>
inline const Particle& getparticle(const ParticleArray& pa, const ID i)
{
  return pa[i];
}
template<typename ID>
inline void setparticle(CombinedParticleArray& pa, const Particle p, const ID i)
{
  pa.setparticle(p,i);
}
template<typename ID>
inline void setparticle(ParticleArray& pa, const Particle p, const ID i)
{
  pa[i] = p;
}

template<class PA, typename ID>
inline const Position& getpos(const PA& pa, const ID i)
{
  return pa[i].position; 
}
template<typename ID>
inline const Position& getpos(const CombinedParticleArray& pa, const ID i)
{
  return pa.poscharge[i].position; 
}
template<typename ID>
inline const Position& getpos(const GhostParticleArray& pa, const ID i)
{
  return pa.poscharge[i].position; 
}
template <class GPA, typename ID>
inline Position& getpos(GPA &gpa, const ID i)
{
  return gpa.poscharge[i].position;
}
template <typename ID>
inline Position& getpos(ParticleArray &gpa, const ID i)
{
  return gpa[i].position;
}
template <typename ID>
inline Position& getpos(PosVelArray &gpa, const ID i)
{
  return gpa[i].position;
}
template <typename ID>
inline Position& getpos(PosChargeArray &gpa, const ID i)
{
  return gpa[i].position;
}

template <class GPA, typename ID>
inline double& getcharge(GPA &gpa, const ID i)
{
  return gpa.poscharge[i].charge;
}
template <class GPA, typename ID>
inline const double& getcharge(const GPA &gpa, const ID i)
{
  return gpa.poscharge[i].charge;
}
template <typename ID>
inline double& getcharge(ParticleArray &gpa, const ID i)
{
  return gpa[i].charge;
}
template <typename ID>
inline const double& getcharge(const ParticleArray &gpa, const ID i)
{
  return gpa[i].charge;
}
template <typename ID>
inline double& getcharge(PosChargeArray &gpa, const ID i)
{
  return gpa[i].charge;
}

template <class GPA, typename ID>
inline Atomtype& getatomtype(GPA &gpa, const ID i)
{
  return gpa.atomtype[i];
}
template <class GPA, typename ID>
inline const Atomtype& getatomtype(const GPA &gpa, const ID i)
{
  return gpa.atomtype[i];
}
template <typename ID>
inline Atomtype& getatomtype(ParticleArray &gpa, const ID i)
{
  return gpa[i].atomtype;
}
template <typename ID>
inline const Atomtype& getatomtype(const ParticleArray &gpa, const ID i)
{
  return gpa[i].atomtype;
}

template <class GPA, typename ID>
inline AtomID& getatomid(GPA &gpa, const ID i)
{
  return gpa.atomid[i];
}
template <class GPA, typename ID>
inline const AtomID& getatomid(const GPA &gpa, const ID i)
{
  return gpa.atomid[i];
}
template <typename ID>
inline AtomID& getatomid(ParticleArray &gpa, const ID i)
{
  return gpa[i].atomid;
}
template <typename ID>
inline const AtomID& getatomid(const ParticleArray &gpa, const ID i)
{
  return gpa[i].atomid;
}

template <class GPA, typename ID>
inline Velocity& getvelocity(GPA &gpa, const ID i)
{
  return gpa.parameters[i].velocity;
}
template <class GPA, typename ID>
inline const Velocity& getvelocity(const GPA &gpa, const ID i)
{
  return gpa.parameters[i].velocity;
}
template <typename ID>
inline Velocity& getvelocity(ParticleArray &gpa, const ID i)
{
  return gpa[i].velocity;
}
template <typename ID>
inline const Velocity& getvelocity(const ParticleArray &gpa, const ID i)
{
  return gpa[i].velocity;
}
template <typename ID>
inline Velocity& getvelocity(PosVelArray &gpa, const ID i)
{
  return gpa[i].velocity;
}
template <typename ID>
inline const Velocity& getvelocity(const PosVelArray &gpa, const ID i)
{
  return gpa[i].velocity;
}

template <class GPA, typename ID>
inline Force& getforce(GPA &gpa, const ID i)
{
  return gpa.force[i];
}
template <class GPA, typename ID>
inline const Force& getforce(const GPA &gpa, const ID i)
{
  return gpa.force[i];
}
template <typename ID>
inline Force& getforce(ParticleArray &gpa, const ID i)
{
  return gpa[i].force;
}
template <typename ID>
inline const Force& getforce(const ParticleArray &gpa, const ID i)
{
  return gpa[i].force;
}

template <class PA, typename ID>
inline double& getmass(PA &pa, const ID i)
{
  return pa.parameters[i].mass;
}
template <typename ID>
inline double& getmass(ParticleArray &pa, const ID i)
{
  return pa[i].mass;
}
template <class PA, typename ID>
inline const double& getmass(const PA &pa, const ID i)
{
  return pa.parameters[i].mass;
}
template <typename ID>
inline const double& getmass(const ParticleArray &pa, const ID i)
{
  return pa[i].mass;
}
template <class PA, typename ID>
inline double& getinvmass(PA &pa, const ID i)
{
  return pa.parameters[i].inv_mass;
}
template <typename ID>
inline double& getinvmass(ParticleArray &pa, const ID i)
{
  return pa[i].inv_mass;
}
template <class PA, typename ID>
inline const double& getinvmass(const PA &pa, const ID i)
{
  return pa.parameters[i].inv_mass;
}
template <typename ID>
inline const double& getinvmass(const ParticleArray &pa, const ID i)
{
  return pa[i].inv_mass;
}

template<class PA>
inline void append(PA& src, const PA& addpa)
{
}

template<>
inline void append(CombinedParticleArray& src, const CombinedParticleArray& addpa)
{
  src.poscharge.insert(src.poscharge.end(),addpa.poscharge.begin(),addpa.poscharge.end());
  src.atomtype.insert(src.atomtype.end(),addpa.atomtype.begin(),addpa.atomtype.end());
  src.force.insert(src.force.end(),addpa.force.begin(),addpa.force.end());
  src.atomid.insert(src.atomid.end(),addpa.atomid.begin(),addpa.atomid.end());
  src.parameters.insert(src.parameters.end(),addpa.parameters.begin(),addpa.parameters.end());
}

template<>
inline void append(ParticleArray& src, const ParticleArray& addpa)
{
  src.insert(src.end(),addpa.begin(),addpa.end());
}


struct WHPair{
  int h1;
  int h2;
  WHPair(){}
  WHPair(int const _h1, int const _h2)
      :h1(_h1), h2(_h2){}
};

struct H1List{
  int nh1;
  int h1[MAX_HBOND];
  int bondtype[MAX_HBOND];
  H1List(){}
  H1List(int const _nh1, int const _h1, int const _bondtype=-1){
    nh1 = 1;
    h1[0] = _h1;
    h1[1] = 0;
    h1[2] = 0;
    h1[3] = 0;
    bondtype[0] = _bondtype;
    bondtype[1] = -1;
    bondtype[2] = -1;
    bondtype[3] = -1;
  }
};

class WaterList : public std::map<int, WHPair>
{
 public:
  WaterList(){}
#if 1
  template<class PA>
  explicit WaterList(const PA& pa) {   // by AH. is this ok?
    for (size_t i = 0; i < pa.size(); ++i){
      if (getatomtype(pa,i) == ATOMTYPE_WO) {
        add_water(i, i+1, i+2);
      }
    }
  }
#endif
  virtual ~WaterList(){}

  std::map<int, int> reverse_list;

  void delete_water(const int o){
    iterator iter = find(o);
#if DEBUG
    if(iter == end()){
      printf("Error (delete_water): %d is no WaterO\n",o);
      exit(1);
    }
#endif
    erase(o);
    reverse_list.erase(iter->second.h1);
    reverse_list.erase(iter->second.h2);
  }

  void add_water(const int o, const int h1, const int h2){
#if DEBUG
    iterator o_iter = find(o);
    std::map<int, int>::iterator h1_iter = reverse_list.find(h1);
    std::map<int, int>::iterator h2_iter = reverse_list.find(h2);
    if(o_iter != end()){
      printf("Error (add_water): WaterO %d already exists\n",o);
      exit(1);
    }
    if(h1_iter != reverse_list.end()){
      printf("Error (add_water): WaterH1 %d already exists\n",h1);
      exit(1);
    }
    if(h2_iter != reverse_list.end()){
      printf("Error (add_water): WaterH2 %d already exists\n",h2);
      exit(1);
    }
#endif
    insert(std::make_pair(o, WHPair(h1, h2)));
    reverse_list.insert(std::make_pair(h1, o));
    reverse_list.insert(std::make_pair(h2, o));
  }
  
  void add_water(const WaterList::const_iterator& w){
    int o = w->first;
    int h1 = w->second.h1;
    int h2 = w->second.h2;
    insert(std::make_pair(o, WHPair(h1, h2)));
    reverse_list.insert(std::make_pair(h1,o));
    reverse_list.insert(std::make_pair(h2,o));
  }
  void append(const WaterList& wl){
    for(WaterList::const_iterator it=wl.begin();it!=wl.end();++it){
      add_water(it);
    }
  }

  void move_o(const int oold, const int onew){
    iterator oold_iter = find(oold);
#if DEBUG
    if(oold_iter == end()){
      printf("Error (WaterList::move_o): WaterO %d is not WaterO\n",oold);
      exit(1);
    }
#endif		
    int h1 = oold_iter->second.h1;
    int h2 = oold_iter->second.h2;

    // rewrite water_list
    (*this)[onew] = oold_iter->second;
    erase(oold);

    // rewrite reverse_list
    std::map<int,int>::iterator h1_iter = reverse_list.find(h1);
#if DEBUG
    if(h1_iter == reverse_list.end()){
      printf("Error (WaterList::move_o): old O's H1 is not in reverse_list: O:%d, H1:%d \n",oold, h1);
      exit(1);
    }			
#endif
    h1_iter->second = onew;

    std::map<int,int>::iterator h2_iter = reverse_list.find(h2);
#if DEBUG
    if(h2_iter == reverse_list.end()){
      printf("Error (WaterList::move_o): old O's H2 is not in reverse_list: O:%d, H2:%d \n",oold, h2);
      exit(1);
    }			
#endif
    h2_iter->second = onew;
  }

  void move_h(const int hold, const int hnew){
    std::map<int,int>::iterator hold_iter = reverse_list.find(hold);
#if DEBUG
    if(hold_iter == reverse_list.end()){
      printf("Error (Water_list::move_h): %d is not WaterH\n",hold);
      exit(1);
    }
#endif
    int o = hold_iter->second;

    // update reverse_list
    reverse_list[hnew] = o;
    reverse_list.erase(hold);
			
    // update water_list;
    iterator o_iter = find(o);
#if DEBUG
    if(o_iter == end()){
      printf("Error (WaterList::move_h): %d is not WaterO\n",o);
      exit(1);
    }
#endif
    if(o_iter->second.h1 == hold)
      o_iter->second.h1 = hnew;
    else if(o_iter->second.h2 == hold)
      o_iter->second.h2 = hnew;
    else
      printf("Error (Water_list::move_h): old H indicates wrong  O (not in water_list): H(%d),O(%d)\n",hold,o);
  }
};

class ShakeList : public std::map<int, H1List>
{
 public:
  ShakeList(){}
#if 1
  template<class PA>
  explicit ShakeList(const PA& pa ) {
    int i_ha = 1;  // if pa[0] = H atom, pa[1] = Heavy atom
    for(size_t i = 0; i < pa.size(); i++){
      if(getatomtype(pa,i) == ATOMTYPE_WO) {
        i+=2; // WH
      }else if(getatomtype(pa,i) < ATOMTYPE_WH){
        add_shake(i_ha, i);
      }else{
        i_ha = i;
      }
    }
  }
#endif
  virtual ~ShakeList(){}

  std::map<int, int> reverse_list;

  void delete_shake(const int ha){  // ha: heavy atom
    iterator iter = find(ha);
#if DEBUG
    if(iter == end()){
      printf("Error (delete_shake): %d is no ShakeHA\n",ha);
      exit(1);
    }
#endif
    // h1 is multiple
    for(int n1=0; n1<iter->second.nh1; n1++){
      reverse_list.erase(iter->second.h1[n1]);
    }
    erase(ha);
  }
  /*
  void add_shake(const int ha, const int h1){  // ha : heavy atom  h1 : H atom
    iterator ha_iter = find(ha);
    if(ha_iter != end()){  // already exists ha
      int nh1 = ha_iter->second.nh1;
      ha_iter->second.h1[nh1] = h1;
      ha_iter->second.nh1 = ++nh1;
#if DEBUG
      if(nh1 > MAX_HBOND){
        printf("Error (add_shake): nh1=%d more than %d\n",nh1, MAX_HBOND);
        exit(1);
      }
#endif
    }else{                 // not exists ha
      insert(std::make_pair(ha, H1List(1, h1)));
    }
    reverse_list.insert(std::make_pair(h1, ha));
  }
  */
  void add_shake(const int ha, const int h1, const int bondtype=-1){  // ha : heavy atom  h1 : H atom  bondtype : bond type
    iterator ha_iter = find(ha);
    if(ha_iter != end()){  // already exists ha
      int nh1 = ha_iter->second.nh1;
      ha_iter->second.h1[nh1] = h1;
      ha_iter->second.bondtype[nh1] = bondtype;
      ha_iter->second.nh1 = ++nh1;
#if DEBUG
      if(nh1 > MAX_HBOND){
        printf("Error (add_shake): nh1=%d more than %d\n",nh1, MAX_HBOND);
        exit(1);
      }
#endif
    }else{                 // not exists ha
      insert(std::make_pair(ha, H1List(1, h1, bondtype)));
    }
    reverse_list.insert(std::make_pair(h1, ha));
  }
  void add_shake(const int ha, const H1List& h1list){
    insert(std::make_pair(ha,h1list));
    for(int h=0;h<h1list.nh1;h++){
      reverse_list.insert(std::make_pair(h1list.h1[h], ha));
    }
  }
  void append(const ShakeList& sl){
    for(ShakeList::const_iterator it=sl.begin();it!=sl.end();++it){
      add_shake(it->first,it->second);
    }
  }

  void change_bondtype(const int ha, const int h1, const int bondtype){
    iterator ha_iter = find(ha);
    if(ha_iter != end()){ // found heavy atom
      int nh1 = ha_iter->second.nh1;
      int h;
      for(h=0;h<nh1;h++){
	if(ha_iter->second.h1[h] == h1)break;
      }
      if(h<nh1){ // found H
	ha_iter->second.bondtype[h] = bondtype;
      }else{
	printf("Error (change_bondtype): not found H atom %d for heavy atom %d\n",h1,ha);
      }
    }else{ // not found heavy atom
      printf("Error (change_bondtype): not found heavy atom %d\n",ha);
    }
  }

  void move_ha(const int haold, const int hanew){
    iterator haold_iter = find(haold);
    if(haold_iter == end()){
#if DEBUG
        printf("Check Error or Normal? (ShakeList::move_ha): old HA' is not in shakelist: HAold:%d HAnew:%d\n",haold,hanew);
#endif
      return;
    }
    int nh1 = haold_iter->second.nh1;

    // rewrite reverse_list
    for(int n1=0; n1<nh1; n1++){
      int h1 = haold_iter->second.h1[n1]; 
      std::map<int,int>::iterator h1_iter = reverse_list.find(h1);
#if DEBUG
      if(h1_iter == reverse_list.end()){
        printf("Error (ShakeList::move_ha): old HA's H1 is not in reverse_list: HAold:%d, H1:%d \n",haold, h1);
        exit(1);
      }
#endif
      if(h1_iter == reverse_list.end()){
        printf("Error (ShakeList::move_ha): old HA's H1 is not in reverse_list: HAold:%d, H1:%d \n",haold, h1);
      }
      
      h1_iter->second = hanew;
    }

    // rewrite shake_list
    (*this)[hanew] = haold_iter->second;
    erase(haold);
  }

  void move_h1(const int h1old, const int h1new){
    std::map<int,int>::iterator h1old_iter = reverse_list.find(h1old);
    if(h1old_iter == reverse_list.end()){
#if DEBUG
        printf("Check Error or Normal? (ShakeList::move_h1): old H1' is not in reverse_list: H1:%d\n",h1old);
#endif
      return;
    }
    int ha = h1old_iter->second;

    // update reverse_list
    reverse_list[h1new] = ha;
    reverse_list.erase(h1old);

    // update shake_list;
    iterator ha_iter = find(ha);
#if DEBUG
    if(ha_iter == end()){
      printf("Error (ShakeList::move_h1): %d is not ShakeHA\n",ha);
      exit(1);
    }
#endif
    bool move = false;
    int nh1 = ha_iter->second.nh1;
    for(int n1=0; n1<nh1; n1++){
      if(ha_iter->second.h1[n1] == h1old){
        ha_iter->second.h1[n1] = h1new;
        move = true;
      }
    }
    if(!move){
      printf("Error (Shake_list::move_h1): old H1 indicates wrong HA (not in shake_list): H1(%d),HA(%d)\n",h1old,ha);
    }
  }
};

template<typename PA>
inline
bool add_particle(PA& particle,
                  TypeRange& tr,
                  WaterList& waterlist,
                  ShakeList& shakelist,
                  const Particle& newparticle,
                  PotentialModel pm)
{
  bool stat=true;
  //  std::cout << "type range " << tr.begin << "--" << tr.end << std::endl;
  if(pm==OnlyLJPotential){
    if(tr.coulomb.end>tr.coulomb.begin){
      //      particle[tr.coulomb.end] = particle[tr.coulomb.begin];
      setparticle(particle,getparticle(particle,tr.coulomb.begin),tr.coulomb.end);
      if(getatomtype(particle,tr.coulomb.end) == ATOMTYPE_WH){
        waterlist.move_h(tr.coulomb.begin, tr.coulomb.end);
      }
    }
    if(tr.ljcoulomb.end>tr.ljcoulomb.begin){
      //      particle[tr.ljcoulomb.end] = particle[tr.ljcoulomb.begin];
      setparticle(particle,getparticle(particle,tr.ljcoulomb.begin),tr.ljcoulomb.end);
      if(getatomtype(particle,tr.ljcoulomb.begin) == ATOMTYPE_WO){
        waterlist.move_o(tr.ljcoulomb.begin, tr.ljcoulomb.end);
      }
      //#ifdef USE_SHAKE
#if 1
      else if(getatomtype(particle,tr.ljcoulomb.begin) < ATOMTYPE_WH){ // H atom, not WH atom
        shakelist.move_h1(tr.ljcoulomb.begin, tr.ljcoulomb.end);
      }else{
        shakelist.move_ha(tr.ljcoulomb.begin, tr.ljcoulomb.end);
      }
#endif
    }
    //    particle[tr.lj.end] = newparticle;
    setparticle(particle,newparticle,tr.lj.end);
    tr.end++;
    tr.lj.end++;
    tr.ljcoulomb.shift(1);
    tr.coulomb.shift(1); 
  }else if(pm==LJCoulombPotential){
    if(tr.coulomb.end>tr.coulomb.begin){
      //      particle[tr.coulomb.end] = particle[tr.coulomb.begin];
      setparticle(particle,getparticle(particle,tr.coulomb.begin),tr.coulomb.end);
      if(getatomtype(particle,tr.coulomb.end) == ATOMTYPE_WH){
        waterlist.move_h(tr.coulomb.begin, tr.coulomb.end);
      }
    }
    //    particle[tr.ljcoulomb.end] = newparticle;
    setparticle(particle,newparticle,tr.ljcoulomb.end);
    tr.end++;
    tr.ljcoulomb.end++;
    tr.coulomb.shift(1);
  }else if(pm==OnlyCoulombPotential){
    //    particle[tr.coulomb.end] = newparticle;
    setparticle(particle,newparticle,tr.coulomb.end);
    tr.end++;
    tr.coulomb.end++;
  }else{
    std::cout << " unknown type particle " << pm << std::endl;
    //    particle[tr.end] = newparticle;
    setparticle(particle,newparticle,tr.end);
    tr.end++;
    stat = false;
  }
  return stat;
}
template<typename PA>
inline
bool add_water(PA& particle,
               TypeRange& tr,
               WaterList& waterlist,
               const Particle& newparticle_o,
               const Particle& newparticle_h1,
               const Particle& newparticle_h2 )
{
  bool stat=true;
  if(tr.coulomb.end>tr.coulomb.begin){
    //    particle[tr.coulomb.end] = particle[tr.coulomb.begin];
    setparticle(particle,getparticle(particle,tr.coulomb.begin),tr.coulomb.end);
    // water
    if(getatomtype(particle,tr.coulomb.end) == ATOMTYPE_WH){
      waterlist.move_h(tr.coulomb.begin, tr.coulomb.end);
    }
  }
  //particle[tr.ljcoulomb.end] = newparticle_o;
  //particle[tr.coulomb.end+1] = newparticle_h1;
  //particle[tr.coulomb.end+2] = newparticle_h2;
  setparticle(particle,newparticle_o,tr.ljcoulomb.end);
  setparticle(particle,newparticle_h1,tr.coulomb.end+1);
  setparticle(particle,newparticle_h2,tr.coulomb.end+2);
	
  waterlist.add_water(tr.ljcoulomb.end, tr.coulomb.end+1, tr.coulomb.end+2);
  tr.end+=3;
  tr.ljcoulomb.end++;
  tr.coulomb.begin++;
  tr.coulomb.end+=3;

  return stat;
}

template<typename PA>
inline
bool add_shake(PA& particle,
               TypeRange& tr,
               WaterList& waterlist,
               ShakeList& shakelist,
               int nh1,
               const Particle& newparticle_ha,
               const Particle& newparticle_h1,
	       const int bondtype=-1)
{
  bool stat=true;
  // atoms of shake bond : LJCoulombPotential
  if(nh1==0){ // heavy atom at only first time
    if(tr.coulomb.end>tr.coulomb.begin){ // exist Coulomb only particle
      /* this case, 1st Coulomb only particle must be moved to last+1
       */
      //      particle[tr.coulomb.end] = particle[tr.coulomb.begin];
      setparticle(particle,getparticle(particle,tr.coulomb.begin),tr.coulomb.end);
      if(getatomtype(particle,tr.coulomb.end) == ATOMTYPE_WH){
        waterlist.move_h(tr.coulomb.begin, tr.coulomb.end);
      }
    }

    //    particle[tr.ljcoulomb.end] = newparticle_ha;
    setparticle(particle,newparticle_ha,tr.ljcoulomb.end);
    tr.end++;
    tr.ljcoulomb.end++;
    tr.coulomb.shift(1);
  }
  int ha_index = tr.ljcoulomb.end;

  // H atom
  int h1_index;
  if(tr.coulomb.end>tr.coulomb.begin){  // move 1st Coulomb only particle
    //    particle[tr.coulomb.end] = particle[tr.coulomb.begin];
    setparticle(particle,getparticle(particle,tr.coulomb.begin),tr.coulomb.end);
    if(getatomtype(particle,tr.coulomb.end) == ATOMTYPE_WH){
      waterlist.move_h(tr.coulomb.begin, tr.coulomb.end);
    }
  }
  //  particle[tr.ljcoulomb.end] = newparticle_h1;
  setparticle(particle,newparticle_h1,tr.ljcoulomb.end);
  h1_index = tr.ljcoulomb.end;
  tr.end++;
  tr.ljcoulomb.end++;
  tr.coulomb.shift(1);

  shakelist.add_shake(h1_index-1-nh1, h1_index, bondtype);

  return stat;
}

template<class PA>
inline
void make_atomid_map(const PA& particlearray,
                     const std::vector<TypeRange>& typerangearray,
                     std::map<AtomID,int>& atomid_to_index)
{
  for(size_t si=0;si<typerangearray.size();si++){
    for(int pi=typerangearray[si].begin;pi<typerangearray[si].end;pi++){
      atomid_to_index.insert(std::pair<AtomID,int>(getatomid(particlearray,pi),pi));
    }
  }
}

template<class SPA, class TPA>
inline 
void convert_particle(const SPA& spa,
                      const std::vector<TypeRange>& typerangearray,
                      TPA& tpa)
{
  for(size_t si=0;si<typerangearray.size();si++){
    for(int i=typerangearray[si].begin;i<typerangearray[si].end;i++){
      getpos(tpa,i) = getpos(spa,i);
      getcharge(tpa,i) = getcharge(spa,i);
      getatomtype(tpa,i) = getatomtype(spa,i);
      getatomid(tpa,i) = getatomid(spa,i);
      getforce(tpa,i) = getforce(spa,i);
      getvelocity(tpa,i) = getvelocity(spa,i);
      getmass(tpa,i) = getmass(spa,i);
      getinvmass(tpa,i) = getinvmass(spa,i);
    }
  }
}


#if 0

//! append  newparticle to last of range for type specified by rt,  original rage is tr and modified tr
template<typename PA>
inline
bool add_particle(PA& particle,
                  TypeRange& tr,
                  const Particle& newparticle,
                  PotentialModel pm)
{
  bool stat=true;
  //  std::cout << "type range " << tr.begin << "--" << tr.end << std::endl;
  if(pm==OnlyLJPotential){
    if(tr.coulomb.end>tr.coulomb.begin){
      particle[tr.coulomb.end] = particle[tr.coulomb.begin];
    }
    if(tr.ljcoulomb.end>tr.ljcoulomb.begin){
      particle[tr.ljcoulomb.end] = particle[tr.ljcoulomb.begin];
    }
    particle[tr.lj.end] = newparticle;
    tr.end++;
    tr.lj.end++;
    tr.ljcoulomb.shift(1);
    tr.coulomb.shift(1);
  }else if(pm==LJCoulombPotential){
    if(tr.coulomb.end>tr.coulomb.begin){
      particle[tr.coulomb.end] = particle[tr.coulomb.begin];
    }
    particle[tr.ljcoulomb.end] = newparticle;
    tr.end++;
    tr.ljcoulomb.end++;
    tr.coulomb.shift(1);
  }else if(pm==OnlyCoulombPotential){
    particle[tr.coulomb.end] = newparticle;
    tr.end++;
    tr.coulomb.end++;
  }else{
    std::cout << " unsupported potentail model " << std::endl;
    particle[tr.end] = newparticle;
    tr.end++;
    stat = false;
  }
  return stat;
}

#endif

#endif

