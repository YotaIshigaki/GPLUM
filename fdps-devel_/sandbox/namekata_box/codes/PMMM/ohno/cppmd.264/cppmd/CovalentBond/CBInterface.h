#ifndef CBINTERFACE_H
#define CBINTERFACE_H

#include <cassert>
#include "ParticleInfo.h"
#include "CovalentBondInfo.h"
#include "CBInterfaceFwd.h"
#include "CBObjects.h"
#include "ShortRangeInteraction.h"
#include "CBUnordered.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdio>

namespace CBModule {
struct CBContext {
  typedef std::vector<CovalentBondInfo::Bond> BondArray;
  typedef std::vector<CovalentBondInfo::Angle> AngleArray;
  typedef std::vector<CovalentBondInfo::Torsion> TorsionArray;
  typedef std::vector<CovalentBondInfo::Improper> ImproperArray;

#ifdef ORG_CB
  CBContext(const BondArray& bond, double& _energy, double& _shortenergy, double& _shortvirial)
  : pBond(&bond), pAngle(), pTorsion(), pImproper(), energy(_energy), shortenergy(_shortenergy), shortvirial(_shortvirial) {}
  CBContext(const AngleArray& angle, double& _energy, double& _shortenergy, double& _shortvirial)
    : pBond(), pAngle(&angle), pTorsion(), pImproper(), energy(_energy), shortenergy(_shortenergy), shortvirial(_shortvirial) {}
  CBContext(const TorsionArray& torsion, double& _energy, double& _shortenergy, double& _shortvirial)
    : pBond(), pAngle(), pTorsion(&torsion), pImproper(), energy(_energy), shortenergy(_shortenergy), shortvirial(_shortvirial) {}
  CBContext(const ImproperArray& improper, double& _energy, double& _shortenergy, double& _shortvirial)
    : pBond(), pAngle(), pTorsion(), pImproper(&improper), energy(_energy), shortenergy(_shortenergy), shortvirial(_shortvirial) {}

  const BondArray* const pBond;
  const AngleArray* const pAngle;
  const TorsionArray* const pTorsion;
  const ImproperArray* const pImproper;
#endif

  double& energy;
  double& shortenergy;
  double& shortvirial;
};

class CBInterface {
public:

  /* call from CovalentBond */
  typedef CBContext::BondArray BondArray;
  typedef CBContext::AngleArray AngleArray;
  typedef CBContext::TorsionArray TorsionArray;
  typedef CBContext::ImproperArray ImproperArray;

#if TUNE_CBMAP
  struct PI2 { CBInterface_::ParticleIndex index_of_atom[2]; };
  struct PI3 { CBInterface_::ParticleIndex index_of_atom[3]; };
  struct PI4 { CBInterface_::ParticleIndex index_of_atom[4]; };
  std::vector<PI2> bondIndexArray;
  std::vector<PI3> angleIndexArray;
  std::vector<PI4> torsionIndexArray;
  std::vector<PI4> improperIndexArray;
#endif  // TUNE_CBMAP


  CBInterface() : cbObjects(this), particleMap(), forceMap(), pParameterList(),
    scnb(2.0), scee(1.2), scnbInv(1.0/scnb), sceeInv(1.0/scee), 
    cancel_short(true), cancel_short_virial(false), cb14energy(0.0)
    {
      num_threads = 1;
    }

  void initialize()
  {
    //    printf("CBInterface() \n");
#ifdef _OPENMP
#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
    }
#endif      
    //    printf("CBInterface num_threads = %d \n",num_threads);
    tmpforce = new ForceArray[num_threads][2];
  }

  // DEBUG
  std::vector<int> c0;
  std::vector<int> c1;
  std::vector<double> ce;
  std::vector<int> ct;
  std::vector<int> cc;
  std::vector<int> c0_off;
  std::vector<int> c1_off;

  void calcel_debug_clear(){
    c0.clear();
    c1.clear();
    ce.clear();
    ct.clear();
    cc.clear();
    c0_off.clear();
    c1_off.clear();
  }
  void calcel_debug_dump(){
    if(c0.size()>0){
      for(std::vector<int>::size_type i=0;i<c0.size();i++){
	std::cout << " cancel " << ct[i] << " " << cc[i] << " (" << c0[i] << "," << c1[i] << ") " << ce[i] << std::endl;
      }
    }
    if(c0_off.size()>0){
      for(std::vector<int>::size_type i=0;i<c0_off.size();i++){
	std::cout << " cancel (" << c0_off[i] << "," << c1_off[i] << ") cut" << std::endl;
      }
    }
  }

  void setParameterList(const CovalentBondParameterList& param_list) {
    pParameterList = &param_list;
    if(DebugLog::verbose>1){
      std::cout << "set CovalentBondParameterList at " << pParameterList << std::endl;
    }
  } 
  const CovalentBondParameterList* getParameterList() {
    if(DebugLog::verbose>1){
      std::cout << "get CovalentBondParameterList at CBInterface.h" << std::endl;
    }
    return pParameterList;
  } 
  void clear_particle_map() {
    particleMap.clear();
    particleIndexMap.clear();
    forceMap.clear();
  }
  void insert_particle_map(ParticleArray& particleArray) {
    for (ParticleArray::iterator it = particleArray.begin();
         it != particleArray.end();++it) {
      Particle& particle = *it;
      particleMap[particle.atomid] = &particle;
    }
    if (DebugLog::verbose > 1){
      std::cout << "ParticleMap size (no typerangearray) " << particleMap.size() << std::endl;
    }
  }
  void insert_particle_map(ParticleArray& particleArray, ForceArray& forceArray) {
    ForceArray::iterator fit = forceArray.begin();
    for (ParticleArray::iterator it = particleArray.begin();
         it != particleArray.end();++it) {
      Particle& particle = *it;
      Force& force = *fit;
      particleMap[particle.atomid] = &particle;
      forceMap[particle.atomid] = &force;
    }
    if (DebugLog::verbose > 1){
      std::cout << "ParticleMap size (no typerangearray) " << particleMap.size() << std::endl;
    }
  }
  /*
    If particleArray has gap, use this version.
    add by ohno
  */
  void insert_particle_map(ParticleArray& particleArray,
                           const std::vector<TypeRange>& typerangearray) {
    int np=0;
    for (std::vector<TypeRange>::const_iterator tit = typerangearray.begin();
         tit != typerangearray.end();++tit) {
      for (int pi=tit->begin;pi<tit->end;pi++){
         Particle& particle = particleArray[pi];
         //         particleMap[particle.atomid] = &particle;
         particleMap.insert(std::make_pair(particle.atomid,&particle));
         np++;
      }
    }
  }
  template<class PA>
  void insert_particle_map(PA& particleArray,
                           ForceArray& forceArray,
                           const std::vector<TypeRange>& typerangearray) {
    int np=0;
    for (std::vector<TypeRange>::const_iterator tit = typerangearray.begin();
         tit != typerangearray.end();++tit) {
      for (int pi=tit->begin;pi<tit->end;pi++){
         Particle& particle = particleArray[pi];
         Force& force = forceArray[pi];
         //         particleMap[particle.atomid] = &particle;
         particleMap.insert(std::make_pair(particle.atomid,&particle));
         forceMap.insert(std::make_pair(particle.atomid,&force));
         np++;
      }
    }
  }
  void insert_particle_map(CombinedParticleArray& particleArray,
                           ForceArray& forceArray,
                           const std::vector<TypeRange>& typerangearray) {
    poscharge[0] = &(particleArray.poscharge);
    atomtype[0] = &(particleArray.atomtype);
    forcea[0] = &(particleArray.force);
    atomid[0] = &(particleArray.atomid);
#if 1
    for(int t=0;t<num_threads;t++){
      tmpforce[t][0].resize((*forcea[0]).size());
    }
#endif
    int np=0;
    for (std::vector<TypeRange>::const_iterator tit = typerangearray.begin();
         tit != typerangearray.end();++tit) {
      for (int pi=tit->begin;pi<tit->end;pi++){
         Force& force = forceArray[pi];
         CBInterface_::ParticleIndex pidx = {0,pi};
         AtomID aid = (*atomid[0])[pi];
         particleIndexMap.insert(std::make_pair(aid,pidx));
         forceMap.insert(std::make_pair(aid,&force));
         np++;
      }
    }
  }
  void insert_particle_map(ParticleArray& particleArray,
                           ForceArray& forceArray,
                           const std::vector<TypeRange>& typerangearray,
                           CBUnordered<AtomID, CBInterface_::ForceLocation>::Map& cbforcemap) {
    int np=0;
    for (std::vector<TypeRange>::const_iterator tit = typerangearray.begin();
         tit != typerangearray.end();++tit) {
      for (int pi=tit->begin;pi<tit->end;pi++){
         Particle& particle = particleArray[pi];
         Force& force = forceArray[pi];
         //         particleMap[particle.atomid] = &particle;
         particleMap.insert(std::make_pair(particle.atomid,&particle));
         forceMap.insert(std::make_pair(particle.atomid,&force));
         cbforcemap.insert(std::make_pair(particle.atomid,&force));
         np++;
      }
    }
  }
  template<class PA>
  void insert_particle_map(PA& particleArray,
                           ForceArray& forceArray,
                           const std::vector<TypeRange>& typerangearray,
                           const std::set<int>& bondatomidlist,
                           CBUnordered<AtomID, CBInterface_::ForceLocation>::Map& cbforcemap) {
    int np=0;
    for (std::vector<TypeRange>::const_iterator tit = typerangearray.begin();
         tit != typerangearray.end();++tit) {
      for (int pi=tit->begin;pi<tit->end;pi++){
         Particle& particle = particleArray[pi];
         int atomid = particle.atomid;
         std::set<int>::const_iterator bp;
         bp = bondatomidlist.find(atomid);
         if(bp!=bondatomidlist.end()){
           Force& force = forceArray[pi];
           //         particleMap[particle.atomid] = &particle;
           particleMap.insert(std::make_pair(atomid,&particle));
           forceMap.insert(std::make_pair(atomid,&force));
           cbforcemap.insert(std::make_pair(atomid,&force));
           np++;
         }
      }
    }
  }
  void insert_particle_map(ParticleArray& particleArray,
                           ForceArray& forceArray,
                           const std::vector<TypeRange>& typerangearray,
                           const std::set<int>& bondatomidlist,
                           std::vector<int>& cbindexarray) {
    int np=0;
    for (std::vector<TypeRange>::const_iterator tit = typerangearray.begin();
         tit != typerangearray.end();++tit) {
      for (int pi=tit->begin;pi<tit->end;pi++){
         Particle& particle = particleArray[pi];
         int atomid = particle.atomid;
         std::set<int>::const_iterator bp;
         bp = bondatomidlist.find(atomid);
         if(bp!=bondatomidlist.end()){
           Force& force = forceArray[pi];
           //         particleMap[particle.atomid] = &particle;
           particleMap.insert(std::make_pair(atomid,&particle));
           forceMap.insert(std::make_pair(atomid,&force));
           cbindexarray.push_back(pi);
           np++;
         }
      }
    }
  }
  void insert_particle_map(GhostParticleArray& particleArray,
                           ForceArray& forceArray,
                           const std::vector<TypeRange>& typerangearray,
                           const std::set<int>& bondatomidlist,
                           std::vector<int>& cbindexarray) {
    poscharge[1] = &(particleArray.poscharge);
    atomtype[1] = &(particleArray.atomtype);
    forcea[1] = &(particleArray.force);
    atomid[1] = &(particleArray.atomid);
#if 1
    for(int t=0;t<num_threads;t++){
      tmpforce[t][1].resize((*forcea[1]).size());
    }
#endif
    int np=0;
    for (std::vector<TypeRange>::const_iterator tit = typerangearray.begin();
         tit != typerangearray.end();++tit) {
      for (int pi=tit->begin;pi<tit->end;pi++){
         AtomID aid = (*atomid[1])[pi];
         std::set<int>::const_iterator bp;
         bp = bondatomidlist.find(aid);
         if(bp!=bondatomidlist.end()){
           Force& force = forceArray[pi];
           CBInterface_::ParticleIndex pidx = {1,pi};
           particleIndexMap.insert(std::make_pair(aid,pidx));
           forceMap.insert(std::make_pair(aid,&force));
           cbindexarray.push_back(pi);
           np++;
         }
      }
    }
  }
#ifdef TUNE_CBMAP
  void make_index_array(const CovalentBondList& bondlist) {
    // build (Bond|Angle|Torsion|Improper)IndexArray
    bondIndexArray.resize(bondlist.bond.size());
#pragma omp parallel for
    for (int i = 0; i < bondlist.bond.size(); ++i) {
      for (int j = 0; j < 2; ++j) {
        CBInterface::ParticleIndexMap::iterator p =
            particleIndexMap.find(bondlist.bond[i].id_of_atom[j]);
        bondIndexArray[i].index_of_atom[j] = p->second;
      }
    }
    angleIndexArray.resize(bondlist.angle.size());
#pragma omp parallel for
    for (int i = 0; i < bondlist.angle.size(); ++i) {
      for (int j = 0; j < 3; ++j) {
        CBInterface::ParticleIndexMap::iterator p =
            particleIndexMap.find(bondlist.angle[i].id_of_atom[j]);
        angleIndexArray[i].index_of_atom[j] = p->second;
      }
    }
    torsionIndexArray.resize(bondlist.torsion.size());
#pragma omp parallel for
    for (int i = 0; i < bondlist.torsion.size(); ++i) {
      for (int j = 0; j < 4; ++j) {
        CBInterface::ParticleIndexMap::iterator p =
            particleIndexMap.find(bondlist.torsion[i].id_of_atom[j]);
        torsionIndexArray[i].index_of_atom[j] = p->second;
      }
    }
    improperIndexArray.resize(bondlist.improper.size());
#pragma omp parallel for
    for (int i = 0; i < bondlist.improper.size(); ++i) {
      for (int j = 0; j < 4; ++j) {
        CBInterface::ParticleIndexMap::iterator p =
            particleIndexMap.find(bondlist.improper[i].id_of_atom[j]);
        improperIndexArray[i].index_of_atom[j] = p->second;
      }
    }
  }
#endif  // TUNE_CBMAP
  inline size_t sizeofparticleMap()
  {
    return particleMap.size();
  }

  void set_cancel_short() {
    cancel_short = true;
    cbObjects.set_cancel_short();
  }

  void unset_cancel_short() {
    cancel_short = false;
    cbObjects.unset_cancel_short();
  }

  void set_cancel_short_virial() {
    cancel_short_virial = true;
    cbObjects.set_cancel_short_virial();
  }

  void unset_cancel_short_virial() {
    cancel_short_virial = false;
    cbObjects.unset_cancel_short_virial();
  }
  
  template<typename PR>
  void subtractCBParticleForce(double& energy, double& shortenergy, double& shortvirial,
			       PR pi,
                               PR pj,
                               CBInterface_::ForceLocation fi,
                               CBInterface_::ForceLocation fj,
                               const Position& d);
  template<typename PR>
    void calcBond(const BondArray& bond, double& energy, double& shortenergy, double& shortvirial);

  template<typename PR>
    void calcAngle(const AngleArray& angle, double& energy, double& shortenergy, double& shortvirial);

  template<typename PR>
    void calcTorsion(const TorsionArray& torsion, double& energy, double& shortenergy, double& shortvirial);

  template<typename PR>
    void calcImproper(const ImproperArray& improper, double& energy, double& shortenergy, double& shortvirial);

  template <typename T, typename TI>
    void calcCBTorsionImproper(const std::vector<T>& torsion, double& energy, double& shortenergy, double& shortvirial);

# ifdef TUNE_CBMAP
  template<typename PR>
    void calcTorsion(const TorsionArray& torsion, double& energy, double& shortenergy, double& shortvirial, const std::vector<PI4>& torsionIndexArray);

  template<typename PR>
    void calcImproper(const ImproperArray& improper, double& energy, double& shortenergy, double& shortvirial, const std::vector<PI4>& improperIndexArray);

  template <typename T, typename TI>
    void calcCBTorsionImproper(const std::vector<T>& torsion, double& energy, double& shortenergy, double& shortvirial, const std::vector<PI4>& indexArray);
# endif  // TUNE_CBMAP


  void set_cancel_short(const double s0) {
    cbObjects.tmpshort = s0;
  }  
  double get_cancel_short() {
    return cbObjects.tmpshort;
  }
  /* end of call from CovalentBond */

  /* call from CBObjects */
  typedef CovalentBondInfo::Bond Bond;
  typedef CovalentBondInfo::Angle Angle;
  typedef CovalentBondInfo::Torsion Torsion;
  typedef CovalentBondInfo::Improper Improper;
  typedef std::vector<Bond>::const_iterator CBBondIterator;
  typedef std::vector<Angle>::const_iterator CBAngleIterator;
  typedef std::vector<Torsion>::const_iterator CBTorsionIterator;
  typedef std::vector<Improper>::const_iterator CBImproperIterator;
  typedef CBInterface_::ParticleLocation ParticleLocation;
  typedef CBInterface_::ParticleIndex ParticleIndex;
  typedef CBInterface_::Context Context;
  typedef CBUnordered<AtomID, ParticleLocation>::Map ParticleMap;
  typedef CBUnordered<AtomID, ParticleIndex>::Map ParticleIndexMap;
  typedef CBInterface_::ForceLocation ForceLocation;
  typedef CBUnordered<AtomID, ForceLocation>::Map ForceMap;

  inline void noneParticleLocation(ParticleLocation& np){
    np = NULL;
  }
  inline void noneParticleLocation(ParticleIndex& np){
    np.array = 0;
    np.index = 0;
  }

  ParticleIndex getParticleLocation(ParticleIndexMap& pmap, const int i);

  ForceLocation getForceLocation(int i) {
#ifdef HAVE_MAP_AT
    return forceMap.at(i);
#else
    ForceMap::iterator ff = forceMap.find(i);
    if(ff!=forceMap.end()){
      return ff->second;
    }else{
      std::cout << " not found atomid " << i << std::endl;
      return NULL;
    }
#endif
  }
  ParticleLocation scanParticleLocation(int i) {
    std::cout << "scan from particleMap size " << particleMap.size() << std::endl;
    /*
    for(ParticleMap::iterator it = particleMap.begin();
        it != particleMap.end(); ++it ) {
      std::cout << "(" << it->first;
      //      std::cout << "," << it->second;
      std::cout << ")";
    }
    std::cout << std::endl;
    */
    int c=0;
    for( ParticleMap::iterator it = particleMap.begin();
         it != particleMap.end(); ++it ) {
      //      std::cout << c << "(" << it->first << ")" ;
      if(it->first==i){
        return it->second;
      }
      c++;
    }
    return ParticleLocation(NULL);
  }

  template <typename T> void calcInteraction(const T& p,
                                             const double r,
                                             double& potential,
                                             double& force);
  template <typename T> void calcInteraction(const T& p, int iparm,
                                             const double r,
                                             double& potential,
                                             double& force);
  template <typename T> void calcInteraction(const T& p, int iparm,
                                             const double cos_theta,
					     const double sign,
                                             double& potential,
                                             double& force);

  template <typename T> int getParameterNum(const T& p) { return 1; }
  template <typename T> bool isCalc14Interaction(const T& p);

  bool inCutoffSphere(double r2) { return true; }

  void calcInteraction(ParticleIndex pi, ParticleIndex pj,
                       const double scnbInv,
                       const double sceeInv,
                       const Position& d,
                       double& r2,
                       double& shortPot,
                       double& longPot,
                       SpaceVector<double>& force) {
    double pot_nb=0.0, dp_nb=0.0;
    double pot_ee=0.0, dp_ee=0.0;
    int ai = pi.array;
    int aj = pj.array;
    AtomID ii = pi.index;
    AtomID ij = pj.index;
    /*
      This method is for new Particle. Short forces for new Particles are calculated by Pairlist. 
      Without longrange, Pairlist calculate force with shift-function.
    */
#ifdef SIMPLE_CUTOFF
    Interaction(r2,
                (*poscharge[ai])[ii].charge,
                (*poscharge[aj])[ij].charge,
                (*atomtype[ai])[ii],
                (*atomtype[aj])[ij],
                pot_nb, dp_nb,
                pot_ee, dp_ee);
#else
# ifdef CBCANCEL_SHIFTEDCL
    Interaction_LJShiftCoulombShift(r2,
                                    (*poscharge[ai])[ii].charge,
                                    (*poscharge[aj])[ij].charge,
                                    (*atomtype[ai])[ii],
                                    (*atomtype[aj])[ij],
                                    pot_nb, dp_nb,
                                    pot_ee, dp_ee);
# else
    Interaction_LJShiftCoulomb(r2,
			       (*poscharge[ai])[ii].charge,
			       (*poscharge[aj])[ij].charge,
			       (*atomtype[ai])[ii],
			       (*atomtype[aj])[ij],
			       pot_nb, dp_nb,
			       pot_ee, dp_ee);
# endif
#endif
#ifdef CHECK_ENERGY
    ShortRange::short_lj -= pot_nb;
    ShortRange::short_coulomb -= pot_ee;
    ShortRange::short_lj -= scnbInv * pot_nb;
    ShortRange::short_coulomb -= sceeInv * pot_ee;
#endif
    shortPot = scnbInv * pot_nb + sceeInv * pot_ee;
    longPot = 0.0;
    force = d * (scnbInv * dp_nb + sceeInv * dp_ee);
  }

  void calcInteraction(ParticleIndex pi, ParticleIndex pj,
                       const Position& d,
                       double& r2,
                       double& shortPot,
                       double& longPot,
                       SpaceVector<double>& force) {
    shortPot = 0.0;
    longPot = 0.0;
    double dp = 0.0;
    int ai = pi.array;
    int aj = pj.array;
    AtomID ii = pi.index;
    AtomID ij = pj.index;
    /*
      This method is for new Particle. Short forces for new Particles are calculated by Pairlist. 
      Without longrange, Pairlist calculate force with shift-function.
    */
#ifdef SIMPLE_CUTOFF
    Interaction(r2,
                (*poscharge[ai])[ii].charge,
                (*poscharge[aj])[ij].charge,
                (*atomtype[ai])[ii],
                (*atomtype[aj])[ij],
                shortPot, dp);
#else
# ifdef CBCANCEL_SHIFTEDCL
    Interaction_LJShiftCoulombShift(r2,
                                    (*poscharge[ai])[ii].charge,
                                    (*poscharge[aj])[ij].charge,
                                    (*atomtype[ai])[ii],
                                    (*atomtype[aj])[ij],
                                    shortPot, dp);
# else
    Interaction_LJShiftCoulomb(r2,
                                    (*poscharge[ai])[ii].charge,
                                    (*poscharge[aj])[ij].charge,
                                    (*atomtype[ai])[ii],
                                    (*atomtype[aj])[ij],
                                    shortPot, dp);
# endif
#endif
    force = d * dp;
  }

  void calcInteraction(ParticleIndex pi, ParticleIndex pj,
                       const Position& d,
                       double& r2,
                       double& shortPot,
                       double& longPot,
                       SpaceVector<double>& force,
		       double& virial) {
    shortPot = 0.0;
    longPot = 0.0;
    double dp = 0.0;
    int ai = pi.array;
    int aj = pj.array;
    AtomID ii = pi.index;
    AtomID ij = pj.index;
    /*
      This method is for new Particle. Short forces for new Particles are calculated by Pairlist. 
      Without longrange, Pairlist calculate force with shift-function.
    */
#ifdef SIMPLE_CUTOFF
    Interaction(r2,
                (*poscharge[ai])[ii].charge,
                (*poscharge[aj])[ij].charge,
                (*atomtype[ai])[ii],
                (*atomtype[aj])[ij],
                shortPot, dp);
#else
# ifdef CBCANCEL_SHIFTEDCL
    Interaction_LJShiftCoulombShift(r2,
                                    (*poscharge[ai])[ii].charge,
                                    (*poscharge[aj])[ij].charge,
                                    (*atomtype[ai])[ii],
                                    (*atomtype[aj])[ij],
                                    shortPot, dp);
# else
    Interaction_LJShiftCoulomb(r2,
                                    (*poscharge[ai])[ii].charge,
                                    (*poscharge[aj])[ij].charge,
                                    (*atomtype[ai])[ii],
                                    (*atomtype[aj])[ij],
                                    shortPot, dp);
# endif
#endif
    force = d * dp;
    virial = r2 * dp;
  }

  void dump_particleMap()
  {
    if(particleMap.size()>0){
      std::cout << "particleMap " << particleMap.size() << " ";
      CBUnordered<AtomID, ParticleLocation>::Map::iterator p;
      for(p=particleMap.begin();p!=particleMap.end();p++){
        std::cout << "(" << (*p).first << "," << (*p).second << ")"; 
      }
      std::cout << std::endl;
    }

  }

  void clear_tmpforce();
  void merge_tmpforce();


#ifdef ORG_CB
  template <typename T> typename std::vector<T>::const_iterator begin(Context);
  template <typename T> typename std::vector<T>::const_iterator end(Context);

  template <typename PR> PR getParticleLocation(Context pContext,
                                                const AtomID ai);

  template <typename T> ParticleLocation getParticleLocation(Context pContext,
                                                  const T& p, int i) {
    assert(i >= 0 || static_cast<size_t>(i) < sizeof(p.id_of_atom)/sizeof(p.id_of_atom[0]));
# ifdef HAVE_MAP_AT
    return particleMap.at(p.id_of_atom[i]);
# else
    ParticleMap::iterator pf = particleMap.find(p.id_of_atom[i]);
    if(pf!=particleMap.end()){
      return pf->second;
    }else{
      std::cout << " not found atomid " << p.id_of_atom[i] << std::endl;
      return NULL;
    }
# endif
  }
  template <typename T> ParticleIndex getParticleIndex(Context pContext,
                                                  const T& p, int i) {
    assert(i >= 0 || static_cast<size_t>(i) < sizeof(p.id_of_atom)/sizeof(p.id_of_atom[0]));
# ifdef HAVE_MAP_AT
    return particleIndexMap.at(p.id_of_atom[i]);
# else
    ParticleIndexMap::iterator pf = particleIndexMap.find(p.id_of_atom[i]);
    if(pf!=particleIndexMap.end()){
      return pf->second;
    }else{
      std::cout << " not found atomid " << p.id_of_atom[i] << std::endl;
      ParticleIndex p0 = {0,0};
      return p0;
    }
# endif
  }

  template <typename T> ForceLocation getForceLocation(Context pContext,
                                                       const T& p, int i) {
    assert(i >= 0 || static_cast<size_t>(i) < sizeof(p.id_of_atom)/sizeof(p.id_of_atom[0]));
# ifdef HAVE_MAP_AT
    return forceMap.at(p.id_of_atom[i]);
# else
    ForceMap::iterator ff = forceMap.find(p.id_of_atom[i]);
    if(ff!=forceMap.end()){
      return ff->second;
    }else{
      std::cout << " not found atomid " << p.id_of_atom[i] << std::endl;
      return NULL;
    }
# endif
  }

  void calcDistance(Context pContext, ParticleLocation pi, ParticleLocation pj,
                    Position& d) {
    d = pi->position - pj->position;
  }
  void calcDistance(Context pContext, ParticleIndex pi, ParticleIndex pj,
                    Position& d) {
    d = (*poscharge[pi.array])[pi.index].position - (*poscharge[pj.array])[pj.index].position;
  }
  bool isFixed(const Bond& pBond) {
    return pBond.shake;
  }

  template <typename T> void addPotentialEnergy(Context pContext,
                                                const double erg) {
    pContext->energy += erg;
  }
  void addForce(Context pContext, ParticleLocation p, const Force& force) {
    p->force += force;
  }
  void addForce(Context pContext, ParticleIndex p, const Force& force) {
    (*forcea[p.array])[p.index] += force;
  }
  void addshortForce(Context pContext, ForceLocation p, const Force& force) {
    (*p) += force;
  }
  void subForce(Context pContext, ParticleLocation p, const Force& force) {
    p->force -= force;
  }
  void subForce(Context pContext, ParticleIndex p, const Force& force) {
    (*forcea[p.array])[p.index] -= force;
  }
  void subshortForce(Context pContext, ForceLocation p, const Force& force) {
    (*p) -= force;
  }

  template <typename T> bool isCalc14Interaction(Context pContext, const T& p);

  bool inCutoffSphere(Context pContext, double r2) { return true; }
  void addShortPotentialEnergy(Context pContext, const double erg) {
    pContext->shortenergy += erg;
  }
  void addLongPotentialEnergy(Context pContext, const double erg) {
    pContext->energy += erg;
  }
  void addshortVirial(Context pContext, const double virial){
    pContext->shortvirial += virial;
  }

  void calcInteraction(Context pContext,
                       ParticleLocation pi, ParticleLocation pj,
                       const double scnbInv,
                       const double sceeInv,
                       const Position& d,
                       double& r2,
                       double& shortPot,
                       double& longPot,
                       SpaceVector<double>& force) {
    double pot_nb=0.0, dp_nb=0.0;
    double pot_ee=0.0, dp_ee=0.0;
#if !defined(USE_PAIRLIST) || defined(SIMPLE_CUTOFF)
    Interaction<ShortRange::LJ>(r2, *pi, *pj, pot_nb, dp_nb);
    Interaction<ShortRange::Coulomb>(r2, *pi, *pj, pot_ee, dp_ee);
#else
    Interaction_LJShiftCoulombShift(r2,
                                    pi->charge,
                                    pj->charge,
                                    pi->atomtype,
                                    pj->atomtype,
                                    pot_nb, dp_nb,
                                    pot_ee, dp_ee);
#endif
#ifdef CHECK_ENERGY
    ShortRange::short_lj -= pot_nb;
    ShortRange::short_coulomb -= pot_ee;
    ShortRange::short_lj -= scnbInv * pot_nb;
    ShortRange::short_coulomb -= sceeInv * pot_ee;
#endif
    shortPot = scnbInv * pot_nb + sceeInv * pot_ee;
    longPot = 0.0;
    force = d * (scnbInv * dp_nb + sceeInv * dp_ee);
  }
  void calcInteraction(Context pContext,
                       ParticleIndex pi, ParticleIndex pj,
                       const double scnbInv,
                       const double sceeInv,
                       const Position& d,
                       double& r2,
                       double& shortPot,
                       double& longPot,
                       SpaceVector<double>& force) {
    double pot_nb=0.0, dp_nb=0.0;
    double pot_ee=0.0, dp_ee=0.0;
    int ai = pi.array;
    int aj = pj.array;
    AtomID ii = pi.index;
    AtomID ij = pj.index;
    /*
      This method is for new Particle. Short forces for new Particles are calculated by Pairlist. 
      Without longrange, Pairlist calculate force with shift-function.
    */
#ifdef SIMPLE_CUTOFF
    Interaction(r2,
                (*poscharge[ai])[ii].charge,
                (*poscharge[aj])[ij].charge,
                (*atomtype[ai])[ii],
                (*atomtype[aj])[ij],
                pot_nb, dp_nb,
                pot_ee, dp_ee);
#else
    Interaction_LJShiftCoulombShift(r2,
                                    (*poscharge[ai])[ii].charge,
                                    (*poscharge[aj])[ij].charge,
                                    (*atomtype[ai])[ii],
                                    (*atomtype[aj])[ij],
                                    pot_nb, dp_nb,
                                    pot_ee, dp_ee);
#endif
#ifdef CHECK_ENERGY
    ShortRange::short_lj -= pot_nb;
    ShortRange::short_coulomb -= pot_ee;
    ShortRange::short_lj -= scnbInv * pot_nb;
    ShortRange::short_coulomb -= sceeInv * pot_ee;
#endif
    shortPot = scnbInv * pot_nb + sceeInv * pot_ee;
    longPot = 0.0;
    force = d * (scnbInv * dp_nb + sceeInv * dp_ee);
  }

  //! calculate LJC cancel force
  void calcInteraction(Context pContext,
                       ParticleLocation pi, ParticleLocation pj,
                       const Position& d,
                       double& r2,
                       double& shortPot,
                       double& longPot,
                       SpaceVector<double>& force) {
    shortPot = 0.0;
    longPot = 0.0;
    double dp = 0.0;
#if !defined(USE_PAIRLIST) || defined(SIMPLE_CUTOFF)
    Interaction<ShortRange::LJCoulomb>(r2, *pi, *pj, shortPot, dp);
#else
    Interaction_LJShiftCoulombShift(r2,
				    pi->charge,
				    pj->charge,
				    pi->atomtype,
				    pj->atomtype,
				    shortPot, dp);
#endif
    force = d * dp;
  }
  void calcInteraction(Context pContext,
                       ParticleIndex pi, ParticleIndex pj,
                       const Position& d,
                       double& r2,
                       double& shortPot,
                       double& longPot,
                       SpaceVector<double>& force) {
    shortPot = 0.0;
    longPot = 0.0;
    double dp = 0.0;
    int ai = pi.array;
    int aj = pj.array;
    AtomID ii = pi.index;
    AtomID ij = pj.index;
    /*
      This method is for new Particle. Short forces for new Particles are calculated by Pairlist. 
      Without longrange, Pairlist calculate force with shift-function.
    */
#ifdef SIMPLE_CUTOFF
    Interaction(r2,
                (*poscharge[ai])[ii].charge,
                (*poscharge[aj])[ij].charge,
                (*atomtype[ai])[ii],
                (*atomtype[aj])[ij],
                shortPot, dp);
#else
    Interaction_LJShiftCoulombShift(r2,
                                    (*poscharge[ai])[ii].charge,
                                    (*poscharge[aj])[ij].charge,
                                    (*atomtype[ai])[ii],
                                    (*atomtype[aj])[ij],
                                    shortPot, dp);
#endif
    force = d * dp;
  }

  void calcInteraction(Context pContext,
                       ParticleLocation pi, ParticleLocation pj,
                       const Position& d,
                       double& r2,
                       double& shortPot,
                       double& longPot,
                       SpaceVector<double>& force,
		       double& virial) {
    shortPot = 0.0;
    longPot = 0.0;
    double dp = 0.0;
#if !defined(USE_PAIRLIST) || defined(SIMPLE_CUTOFF)
    Interaction<ShortRange::LJCoulomb>(r2, *pi, *pj, shortPot, dp);
#else
    Interaction_LJShiftCoulombShift(r2,
				    pi->charge,
				    pj->charge,
				    pi->atomtype,
				    pj->atomtype,
				    shortPot, dp);
#endif
    force = d * dp;
    virial = r2*dp;
  }
  void calcInteraction(Context pContext,
                       ParticleIndex pi, ParticleIndex pj,
                       const Position& d,
                       double& r2,
                       double& shortPot,
                       double& longPot,
                       SpaceVector<double>& force,
		       double& virial) {
    shortPot = 0.0;
    longPot = 0.0;
    double dp = 0.0;
    int ai = pi.array;
    int aj = pj.array;
    AtomID ii = pi.index;
    AtomID ij = pj.index;
    /*
      This method is for new Particle. Short forces for new Particles are calculated by Pairlist. 
      Without longrange, Pairlist calculate force with shift-function.
    */
#ifdef SIMPLE_CUTOFF
    Interaction(r2,
                (*poscharge[ai])[ii].charge,
                (*poscharge[aj])[ij].charge,
                (*atomtype[ai])[ii],
                (*atomtype[aj])[ij],
                shortPot, dp);
#else
    Interaction_LJShiftCoulombShift(r2,
                                    (*poscharge[ai])[ii].charge,
                                    (*poscharge[aj])[ij].charge,
                                    (*atomtype[ai])[ii],
                                    (*atomtype[aj])[ij],
                                    shortPot, dp);
#endif
    force = d * dp;
    virial = r2 * dp;
  }
#endif
  /* end of call from CBObjects */


  void clear_cb14energy(){
    cb14energy = 0.0;
  }

  double get_cb14energy(){
    return cb14energy;
  }


private:
  typedef CovalentBondInfo::BondParameter BondParameter;
  typedef CovalentBondInfo::AngleParameter AngleParameter;
  typedef CovalentBondInfo::TorsionParameter TorsionParameter;
  typedef CovalentBondInfo::ImproperParameter ImproperParameter;

  template <typename TP> void calcInteractionTorsion(const TP& param,
                                                     const double r,
                                                     double& potential,
                                                     double& force);

  template <typename TP> void calcInteractionTorsion(const TP& param,
                                                     const double theta,
						     const double sign,
                                                     double& potential,
                                                     double& force);

  template<typename T>
    void convert_torsion_atomid_to_index(const ParticleIndexMap& pmap,
					 const std::vector<T>& torsion, 
					 std::vector< std::vector<ParticleIndex> >& torsion_index);
  
  void convert_torsion_improper_atomid_to_index(const TorsionArray& torsion, const ImproperArray& improper);

  CBObjects cbObjects;
  ParticleMap particleMap;
  ParticleIndexMap particleIndexMap;
  PosChargeArray* poscharge[2];
  AtomtypeArray* atomtype[2];
  ForceArray* forcea[2];
  AtomIDArray* atomid[2];
  ForceMap forceMap;
  const CovalentBondParameterList* pParameterList;
  double scnb,scee,scnbInv,sceeInv;
  bool cancel_short;
  bool cancel_short_virial;

  double cb14energy;

  std::vector< std::vector<ParticleIndex> > torsion_index;
  std::vector< std::vector<ParticleIndex> > improper_index;

  int num_threads;
  ForceArray (*tmpforce)[2];
};

}
#endif
