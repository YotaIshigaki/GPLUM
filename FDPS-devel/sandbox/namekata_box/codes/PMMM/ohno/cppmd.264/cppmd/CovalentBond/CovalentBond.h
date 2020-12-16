#ifndef COVALENTBOND_H
#define COVALENTBOND_H

#include "CovalentBondInfo.h"
#include "CBInterface.h"

//! Covalent Bond
template<class CBPARTICLE>
class CovalentBond {
public:
  void setParameterList(const CovalentBondParameterList& param_list) {
    cbInterface.setParameterList(param_list);
    if(DebugLog::verbose>1){
      param_list.print();
    }
  }

  const CovalentBondParameterList* getParameterList() {
    return cbInterface.getParameterList();
  }

  void set_cancel_short() {
    cbInterface.set_cancel_short();
  }

  void unset_cancel_short() {
    cbInterface.unset_cancel_short();
  }

  void set_cancel_short_virial() {
    cbInterface.set_cancel_short_virial();
  }

  void unset_cancel_short_virial() {
    cbInterface.unset_cancel_short_virial();
  }

  void calcForce(ParticleArray& particlearray, ParticleArray& ghost, CovalentBondList& bondlist, double& energy, double& shortenergy){
    cbInterface.clear_particle_map();
    cbInterface.insert_particle_map(particlearray);
    cbInterface.insert_particle_map(ghost);
    cbInterface.calcBond<CBPARTICLE>(bondlist.bond,energy, shortenergy);
    cbInterface.calcAngle<CBPARTICLE>(bondlist.angle,energy, shortenergy);
    cbInterface.calcTorsion<CBPARTICLE>(bondlist.torsion,energy, shortenergy);
    cbInterface.calcImproper<CBPARTICLE>(bondlist.improper,energy, shortenergy);
  }
  void map_particle(ParticleArray& particlearray, 
                    ForceArray& shortforcearray,
                    const std::vector<TypeRange>& typerangearray,
                    ParticleArray& ghost, 
                    ForceArray& ghostshortforcearray,
                    const std::vector<TypeRange>& ghosttyperange)
  {
    cbInterface.clear_particle_map();
    cbInterface.insert_particle_map(particlearray,shortforcearray,typerangearray);
    cbInterface.insert_particle_map(ghost,ghostshortforcearray,ghosttyperange);
    if (DebugLog::verbose > 1){
      std::cout << "map for CB particle " << cbInterface.sizeofparticleMap() << std::endl;
    }
  }
  void map_particle(ParticleArray& particlearray, 
                    ForceArray& shortforcearray,
                    const std::vector<TypeRange>& typerangearray,
                    ParticleArray& ghost, 
                    ForceArray& ghostshortforcearray,
                    const std::vector<TypeRange>& ghosttyperange,
                    CBModule::CBUnordered<AtomID, CBModule::CBInterface_::ForceLocation>::Map& ghostcbforcemap)
  {
    cbInterface.clear_particle_map();
    cbInterface.insert_particle_map(particlearray,shortforcearray,typerangearray);
    cbInterface.insert_particle_map(ghost,ghostshortforcearray,ghosttyperange,
                                    ghostcbforcemap);
    if (DebugLog::verbose > 1){
      std::cout << "map for CB particle " << cbInterface.sizeofparticleMap() << std::endl;
    }
  }
  template<class PA, class GPA>
  void map_particle(PA& particlearray, 
                    ForceArray& shortforcearray,
                    const std::vector<TypeRange>& typerangearray,
                    GPA& ghost, 
                    ForceArray& ghostshortforcearray,
                    const std::vector<TypeRange>& ghosttyperange,
                    const CovalentBondList& bondlist,
                    CBModule::CBUnordered<AtomID, CBModule::CBInterface_::ForceLocation>::Map& ghostcbforcemap)
  {
    cbInterface.clear_particle_map();
    cbInterface.insert_particle_map(particlearray,shortforcearray,typerangearray);
    cbInterface.insert_particle_map(ghost,ghostshortforcearray,ghosttyperange,
                                    bondlist.atomidset,
                                    ghostcbforcemap);
    if (DebugLog::verbose > 1){
      std::cout << "map for CB particle " << cbInterface.sizeofparticleMap() << std::endl;
    }
  }
  template<class PA, class GPA>
  void map_particle(PA& particlearray, 
                    ForceArray& shortforcearray,
                    const std::vector<TypeRange>& typerangearray,
                    GPA& ghost, 
                    ForceArray& ghostshortforcearray,
                    const std::vector<TypeRange>& ghosttyperange,
                    const CovalentBondList& bondlist,
                    std::vector<int>& ghostcbindexarray)
  {
    cbInterface.clear_particle_map();
    cbInterface.insert_particle_map(particlearray,shortforcearray,typerangearray);
    cbInterface.insert_particle_map(ghost,ghostshortforcearray,ghosttyperange,
                                    bondlist.atomidset,
                                    ghostcbindexarray);
#ifdef TUNE_CBMAP
    cbInterface.make_index_array(bondlist);
#endif  // TUNE_CBMAP
    if (DebugLog::verbose > 1){
      std::cout << "map for CB particle " << cbInterface.sizeofparticleMap() << std::endl;
    }
  }
  void calcForce(CovalentBondList& bondlist, double& energy, double& shortenergy){
    std::cout << " CovalentBond "<< std::endl;
    // DEBUG LJC cancel    cbInterface.set_cancel_short(0.0);
    cbInterface.clear_tmpforce();
    cbInterface.calcBond<CBPARTICLE>(bondlist.bond,energy, shortenergy);
    cbInterface.calcAngle<CBPARTICLE>(bondlist.angle,energy, shortenergy);
    cbInterface.calcTorsion<CBPARTICLE>(bondlist.torsion,energy, shortenergy);
    cbInterface.calcImproper<CBPARTICLE>(bondlist.improper,energy, shortenergy);
    cbInterface.merge_tmpforce();
  }
  void calcForce(CovalentBondList& bondlist, double& energy, double& shortenergy, double &shortvirial){
    //    std::cout << " CovalentBond "<< std::endl;
    // DEBUG LJC cancel    cbInterface.set_cancel_short(0.0);
    //    printf("clear before calc*\n");
    cbInterface.clear_tmpforce();
    cbInterface.calcBond<CBPARTICLE>(bondlist.bond,energy, shortenergy, shortvirial);
    cbInterface.calcAngle<CBPARTICLE>(bondlist.angle,energy, shortenergy, shortvirial);
#ifdef TUNE_CBMAP
    cbInterface.calcTorsion<CBPARTICLE>(bondlist.torsion,energy, shortenergy, shortvirial, cbInterface.torsionIndexArray);
    cbInterface.calcImproper<CBPARTICLE>(bondlist.improper,energy, shortenergy, shortvirial, cbInterface.improperIndexArray);
#else  // TUNE_CBMAP
    cbInterface.calcTorsion<CBPARTICLE>(bondlist.torsion,energy, shortenergy, shortvirial);
    cbInterface.calcImproper<CBPARTICLE>(bondlist.improper,energy, shortenergy, shortvirial);
#endif  // TUNE_CBMAP
    cbInterface.merge_tmpforce();
  }
  void calcForce(ParticleArray& particlearray, 
                 const std::vector<TypeRange>& typerangearray,
                 ParticleArray& ghost, 
                 const std::vector<TypeRange>& ghosttyperange,
                 CovalentBondList& bondlist, double& energy, double& shortenergy){
    std::cout << " CovalentBond "<< std::endl;
    cbInterface.clear_particle_map();
    cbInterface.insert_particle_map(particlearray,typerangearray);
    cbInterface.insert_particle_map(ghost,ghosttyperange);
    cbInterface.calcBond<CBPARTICLE>(bondlist.bond,energy, shortenergy);
    cbInterface.calcAngle<CBPARTICLE>(bondlist.angle,energy, shortenergy);
    cbInterface.calcTorsion<CBPARTICLE>(bondlist.torsion,energy, shortenergy);
    cbInterface.calcImproper<CBPARTICLE>(bondlist.improper,energy, shortenergy);
  }
  void dump_particleMap(){
    cbInterface.dump_particleMap();
  }

  void calcel_debug_clear(){
    cbInterface.calcel_debug_clear();
  }
  void calcel_debug_dump(){
    cbInterface.calcel_debug_dump();
  }

  double get_cb14energy(){
    return cbInterface.get_cb14energy();
  }

  void initialize(){
    cbInterface.initialize();
  }
private:
  CBModule::CBInterface cbInterface;
};

#endif

