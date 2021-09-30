//***************************************************************************************
//  This program is the data class and templete structure of "atsevb_main.cpp"
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

//--- base class for particle
class atom_pos{
public:
    //--- common properties
    MOL_TYPE   type;      // molecular type (enum MOL_TYPE)
    PS::S32    mol_id;    // molecular id
    PS::S32    atom_id;   // atom id
    PS::F64vec pos;
    
    //--- member functions
    MOL_TYPE   getType() const {
        return this->type;
    }
    PS::S32    getMolID() const {
        return this->mol_id;
    }
    PS::S32    getAtomID() const {
        return this->atom_id;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
};


//--- Full Particle class
//------ This class should contain all properties of the particle.
class FP :
  public atom_pos{
  public:
    //--- for intramolecular connection information
    PS::S32    bondPair_num;
    PS::S32    bondPair[4];   // atom_id of bonded pair atom (max 4 atom)
    
    //--- for coulomb interaction
//    static PS::F64 RSearch_coulomb;
    PS::F64    charge;
    PS::F64vec dipole;
    PS::F64    alpha;   // polarizability
    
    //--- for VDW interaction
//    static PS::F64 RSearch_VDW;
    PS::F64    VDW_R;
    PS::F64    VDW_D;
    
    //--- result
    PS::F64vec force_LJ;
    PS::F64    pot_LJ;
    PS::F64vec force_coulomb;
    PS::F64    pot_coulomb;
    PS::F64vec force_PM;
    
#ifdef VIRIAL_TEST
    PS::F64vec virt_LJ;
    PS::F64vec virt_coulomb;
#endif
    
    //--- member functions
    PS::S32    getBondPair_num() const {
        return this->bondPair_num;
    }
//    PS::F64    getRSearch() const {
//        return this->RSearch_coulomb;
//    }
    PS::F64    getCharge() const {
        return this->charge;
    }
    PS::F64vec getDipole() const {
        return this->dipole;
    }
    PS::F64    getAlpha() const {
        return this->alpha;
    }
    PS::F64    getVDW_R() const {
        return this->VDW_R;
    }
    PS::F64    getVDW_D() const {
        return this->VDW_D;
    }
    
    void copyFromForce(const ForceCoulomb & force) {
        this->force_coulomb = force.field_coulomb * this->charge;
        this->pot_coulomb   = force.pot_coulomb   * this->charge;
        
#ifdef VIRIAL_TEST
        this->virt_coulomb += force.virt_coulomb;
#endif
    }
    void copyFromForce(const ForceLJ & force) {
        this->force_LJ      = force.force_LJ;
        this->pot_LJ        = force.pot_LJ;
        
#ifdef VIRIAL_TEST
        this->virt_LJ      += force.virt_LJ;
#endif
    }
    void copyFromForce(const ForceScatter & force) {
        this->force_coulomb = force.field_coulomb * this->charge;
        this->pot_coulomb   = force.pot_coulomb   * this->charge;
        this->force_LJ      = force.force_LJ;
        this->pot_LJ        = force.pot_LJ;
        
#ifdef VIRIAL_TEST
        this->virt_LJ      += force.virt_LJ;
        this->virt_coulomb += force.virt_coulomb;
#endif
    }
    //--- for Particle Mesh
    PS::F64 getChargeParticleMesh() const {
        return this->charge;
    }
    void copyFromForceParticleMesh(const PS::F64vec & force_pm) {
        this->force_PM += force_pm * this->charge;
    }
    
    //--- clear dinamic data
    void clear() {
        this->dipole        = 0.0;
        this->force_coulomb = 0.0;
        this->pot_coulomb   = 0.0;
        this->force_LJ      = 0.0;
        this->pot_LJ        = 0.0;
        this->force_PM      = 0.0;
        
#ifdef VIRIAL_TEST
        this->virt_LJ       = 0.0;
        this->virt_coulomb  = 0.0;
        //std::cout << "virial is cleared in FP" << std::endl << std::endl;
#endif
    }
};


//--- Essential Particle class
//------ This class has the meta-data for calculate force.
//------ This class is the subset of Full Particle class.

//--------- for LJ interaction
class EPI_LJ :
  public atom_pos{
  public:
    static PS::F64 RSearch_LJ;
    PS::F64 VDW_R;
    PS::F64 VDW_D;
    
    //--- member functions
    PS::F64 getRSearch() const {
        return this->RSearch_LJ;
    }
    PS::F64 getVDW_R() const {
        return this->VDW_R;
    }
    PS::F64 getVDW_D() const {
        return this->VDW_D;
    }
    
    //--- copy data from FP class
    void copyFromFP(const FP & fp) {
        //--- in atom_pos
        this->type    = fp.type;
        this->atom_id = fp.atom_id;
        this->mol_id  = fp.mol_id;
        this->pos     = fp.pos;
        //--- added in EPI_VDW
        this->VDW_R   = fp.VDW_R;
        this->VDW_D   = fp.VDW_D;
    }
};

PS::F64 EPI_LJ::RSearch_LJ;

class EPJ_LJ :
  public EPI_LJ{
  public:
    void setPos(const PS::F64vec pos_new){
        this->pos = pos_new;
    }
};

//--------- for coulomb interaction (permanent point charge)
class EPI_coulomb :
  public atom_pos{
  public:
    static PS::F64 RSearch_coulomb;
    PS::F64 charge;
    
    //--- member functions
    PS::F64 getRSearch() const {
        return this->RSearch_coulomb;
    }
    PS::F64 getCharge() const {
        return this->charge;
    }
    
    //--- copy data from FP class
    void copyFromFP(const FP & fp) {
        //--- in atom_pos
        this->type    = fp.type;
        this->atom_id = fp.atom_id;
        this->mol_id  = fp.mol_id;
        this->pos     = fp.pos;
        //--- added in EPI_coulomb
        this->charge  = fp.charge;
    }
};

PS::F64 EPI_coulomb::RSearch_coulomb;

class EPJ_coulomb :
  public EPI_coulomb{
  public:
    void setPos(const PS::F64vec pos_new){
        this->pos = pos_new;
    }
};

//--------- for CalcForceScatter
//------------ contains LJ and direct part of Coulomb static
class EPI_scatter :
  public atom_pos{
  public:
    PS::F64 charge;
    PS::F64 VDW_R;
    PS::F64 VDW_D;
    
    PS::F64 getCharge() const {
        return this->charge;
    }
    PS::F64 getVDW_R() const {
        return this->VDW_R;
    }
    PS::F64 getVDW_D() const {
        return this->VDW_D;
    }
    
    //--- copy data from FP class
    void copyFromFP(const FP & fp) {
        //--- in atom_pos
        this->type    = fp.type;
        this->atom_id = fp.atom_id;
        this->mol_id  = fp.mol_id;
        this->pos     = fp.pos;
        //--- added in EPI_scatter
        this->charge  = fp.charge;
        this->VDW_R   = fp.VDW_R;
        this->VDW_D   = fp.VDW_D;
    }
  };
  
class EPJ_scatter :
  public EPI_scatter{
  public:
    void setPos(const PS::F64vec pos_new){
        this->pos = pos_new;
    }
  };


//--------- for coulomb dipole
class EPI_dipole :
  public atom_pos{
  public:
    PS::F64vec dipole;
    PS::F64    alpha;   // polarizability
    
    PS::F64vec getDipole() const {
        return this->dipole;
    }
    PS::F64    getAlpha() const {
        return this->alpha;
    }
    void copyFromFP(const FP & fp) {
        //--- in atom_pos
        this->type    = fp.type;
        this->atom_id = fp.atom_id;
        this->mol_id  = fp.mol_id;
        this->pos     = fp.pos;
        //--- added in EPI_dipole
        this->dipole  = fp.dipole;
        this->alpha   = fp.alpha;
    }
};

class EPJ_dipole :
  public EPI_dipole{
  public:
    void setPos(const PS::F64vec pos_new){
        this->pos = pos_new;
    }
};


//------ for intramolecular force
class EPI_intraforce{
public:
    //--- identifier
    MOL_TYPE   type;      // molecular type (enum MOL_TYPE)
    PS::S32    mol_id;    // molecular id
    PS::S32    atom_id;   // atom id
    
    //--- for intramolecular connection information
    //------ bond connection
    PS::S32    bondPair_num;
    PS::S32    bondPair[4];   // atom_id of bonded pair atom (max 4 atom)
    
    //--- member functions
    void copyFromFP(const FP & fp) {
        this->type    = fp.type;
        this->atom_id = fp.atom_id;
        this->mol_id  = fp.mol_id;
        
        this->bondPair_num = fp.bondPair_num;
        this->bondPair[0]  = fp.bondPair[0];
        this->bondPair[1]  = fp.bondPair[1];
        this->bondPair[2]  = fp.bondPair[2];
        this->bondPair[3]  = fp.bondPair[3];
    }
};
