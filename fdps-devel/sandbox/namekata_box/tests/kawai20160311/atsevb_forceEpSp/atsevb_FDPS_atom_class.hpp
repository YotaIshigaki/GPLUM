//***************************************************************************************
//  This program is the data class and templete structure of "atsevb_main.cpp"
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

//--- base classes for particle data ----------------------------------------------------
//------ identifier
class AtomID{
  private:
    PS::S32   atom_id;   // atom id
    PS::S32   mol_id;    // molecular id
    MOL_TYPE  mol_type;  // molecular type (enum class MOL_TYPE)
    
  public:
    //--- constructor
    AtomID(){
        this->atom_id  = -1;
        this->mol_id   = -1;
        this->mol_type = MOL_TYPE::ENUM_N_TOT;
    }
  
    //--- member functions
    void setAtomID(const PS::S32 & id){
        this->atom_id = id;
    }
    void setMolID(const PS::S32 & id){
        this->mol_id = id;
    }
    void setMolType(const MOL_TYPE & type){
        this->mol_type = type;
    }
    inline PS::S32 getAtomID() const {
        return this->atom_id;
    }
    inline PS::S32 getMolID() const {
        return this->mol_id;
    }
    inline MOL_TYPE getMolType() const {
        return this->mol_type;
    }
    
    //--- copy function
    template<class Tptcl>
    void copyAtomID(const Tptcl & fp){
        this->atom_id   = fp.getAtomID();
        this->mol_id    = fp.getMolID();
        this->mol_type  = fp.getMolType();
    }
};

//------ sub identifier
class AtomType{
  private:
    ATOM_TYPE atom_type;  // atom type (enum class ATOM_TYPE)

  public:
    //--- constructor
    AtomType(){
        this->atom_type = ATOM_TYPE::ENUM_N_TOT;  // illigal value.
    }
  
    void setAtomType(const ATOM_TYPE & type){
        this->atom_type = type;
    }
    inline ATOM_TYPE getAtomType() const {
        return this->atom_type;
    }
    
    template<class Tptcl>
    void copyAtomType(const Tptcl & fp){
        this->atom_type = fp.getAtomType();
    }
};

//--- connection data
class AtomConnect{
  public:
    //--- definition of data size
    MD_EXT::Connection<PS::S32, 4> bond;
    MD_EXT::Connection<PS::S32, 2> EVB;
    
    void getAllConnect(std::vector<PS::S32> & v_tmp) const {
        v_tmp.clear();
        v_tmp.reserve(  this->bond.getN()
                      + this->EVB.getN() );
        //--- add bond connection
        for(PS::S32 i=0; i<this->bond.getN(); i++){
            v_tmp.push_back(this->bond.get(i));
        }
        //-- add EVB connection
        for(PS::S32 i=0; i<this->EVB.getN(); i++){
            v_tmp.push_back(this->EVB.get(i));
        }
    }
    
    //--- copy from FP to EP
    template<class Tptcl>
    void copyAtomConnect(const Tptcl & fp){
        this->bond = fp.bond;
        this->EVB  = fp.EVB;
    }
    void clearAtomConnect(){
        this->bond.clear();
        this->EVB.clear();
    }
};

class IntraMask :
  //public MD_EXT::IntraMaskList<PS::S32>,
  public MD_EXT::IntraMaskList_fix<PS::S32, 16> {
  public:
};

//------ position
class AtomPos{
  private:
    PS::F64vec pos;
    
  public:
    //--- constructor
    AtomPos(){
        this->pos = 0.0;
    }
  
    //--- member functions
    void setPos(const PS::F64vec & pos_new){
        this->pos = pos_new;
    }
    inline PS::F64vec getPos() const {
        return this->pos;
    }
    
    //--- copy function
    template<class Tptcl>
    void copyAtomPos(const Tptcl & fp){
        this->pos = fp.getPos();
    }
};

//------ mass & velocity
class AtomVel{
  private:
    PS::F64    mass;  // mass of atom
    PS::F64    mass_inv;
    PS::F64vec vel;   // velocity of atom
      
  public:
    //--- constructor
    AtomVel(){
        this->mass = 0.0;  // this value will cause floating invalid error.
        this->vel  = 0.0;
    }
  
    void setMass(const PS::F64 & mass_new){
        this->mass     = mass_new;
        this->mass_inv = 1.0/mass_new;
    }
    void setVel(const PS::F64vec & vel_new){
        this->vel = vel_new;
    }
    inline void accel(const PS::F64vec & f,
                      const PS::F64    & dt){
        this->vel += f*dt*mass_inv;
    }
    inline PS::F64    getMass() const {
        return this->mass;
    }
    inline PS::F64vec getVel() const {
        return this->vel;
    }
};

//------ permanent charge
class AtomCharge{
  private:
    PS::F64 charge;
  
  public:
    //--- constructor
    AtomCharge(){
        this->charge = 0.0;
    }
  
    void setCharge(const PS::F64 & q){
        this->charge = q;
    }
    inline PS::F64 getCharge() const {
        return this->charge;
    }
    inline PS::F64 getChargeParticleMesh() const {
        return this->charge;
    }
    
    template<class Tptcl>
    void copyAtomCharge(const Tptcl & fp){
        this->charge = fp.getCharge();
    }
};

//------ VDW parameters
class AtomVDW{
  private:
    PS::F64 VDW_R;
    PS::F64 VDW_D;
  
  public:
    //--- constructor
    AtomVDW(){
        this->VDW_R = 0.0;
        this->VDW_D = 0.0;
    }
  
    void setVDW_R(const PS::F64 & vdw_r){
        this->VDW_R = vdw_r;
    }
    void setVDW_D(const PS::F64 & vdw_d){
        this->VDW_D = vdw_d;
    }
    inline PS::F64 getVDW_R() const {
        return this->VDW_R;
    }
    inline PS::F64 getVDW_D() const {
        return this->VDW_D;
    }
    
    template<class Tptcl>
    void copyAtomVDW(const Tptcl & fp){
        this->VDW_R = fp.getVDW_R();
        this->VDW_D = fp.getVDW_D();
    }
};

//------ induced dipole
class AtomDipole{
  private:
    PS::F64vec dipole;
    PS::F64    alpha;   // polarizability
  
  public:
    //--- constructor
    AtomDipole(){
        this->dipole = 0.0;
        this->alpha  = 0.0;
    }
  
    void setDipole(const PS::F64vec & d){
        this->dipole = d;
    }
    void setAlpha(const PS::F64 & a){
        this->alpha = a;
    }
    inline PS::F64vec getDipole() const {
        return this->dipole;
    }
    inline PS::F64    getAlpha() const {
        return this->alpha;
    }
    
    template<class Tptcl>
    void copyAtomDipole(const Tptcl & fp){
        this->dipole = fp.getDipole();
        this->alpha  = fp.getAlpha();
    }
};

//--- Force class -----------------------------------------------------------------------
//------ This class has the result of force. (used as temporary data)
//------ This class is the subset of Full Particle class.

//------ LJ interaction with cutoff length (simple cutoff)
class ForceLJ {
  private:
    PS::F64vec force_LJ;
    PS::F64vec  virt_LJ;
    PS::F64      pot_LJ;
    
  public:
    void clearForceLJ(){
        this->force_LJ = 0.0;
        this->virt_LJ  = 0.0;
        this->pot_LJ   = 0.0;
    }
    //--- constructor
    ForceLJ(){
        this->clearForceLJ();
    }
    
    PS::F64vec getForceLJ() const {
        return this->force_LJ;
    }
    PS::F64vec getVirtLJ() const {
        return this->virt_LJ;
    }
    PS::F64    getPotLJ() const {
        return this->pot_LJ;
    }
    inline void addForceLJ(const PS::F64vec & force){
        this->force_LJ += force;
    }
    inline void addVirtLJ(const PS::F64vec & virt){
        this->virt_LJ += virt;
    }
    inline void addPotLJ(const PS::F64 & pot){
        this->pot_LJ += pot;
    }
    
    void clear(){
        this->clearForceLJ();
    }
};

//------ coulomb interaction
class ForceCoulomb {
  private:
    PS::F64vec field_coulomb;
    PS::F64      pot_coulomb;
    
  public:
    void clearForceCoulomb(){
        field_coulomb = 0.0;
          pot_coulomb = 0.0;
    }
    //--- constructor
    ForceCoulomb(){
        clearForceCoulomb();
    }
    void clear(){
        clearForceCoulomb();
    }
  
    PS::F64vec getFieldCoulomb() const {
        return this->field_coulomb;
    }
    PS::F64    getPotCoulomb() const {
        return this->pot_coulomb;
    }

    void addFieldCoulomb(const PS::F64vec & f){
        this->field_coulomb += f;
    }
    void addPotCoulomb(const PS::F64 & p){
        this->pot_coulomb += p;
    }
};

//--- Intramoleecular interaction
class ForceIntra {
  private:
    PS::F64vec force;
    PS::F64vec virt;
    PS::F64    pot;
  public:
    void clearForceIntra(){
        this->force = 0.0;
        this->virt  = 0.0;
        this->pot   = 0.0;
    }
    //--- constructor
    ForceIntra(){
        this->clearForceIntra();
    }
  
    inline PS::F64vec getForceIntra() const {
        return this->force;
    }
    inline PS::F64vec getVirtIntra() const {
        return this->virt;
    }
    inline PS::F64    getPotIntra() const {
        return this->pot;
    }
    inline void addForceIntra(const PS::F64vec& f){
        this->force += f;
    }
    inline void addVirtIntra(const PS::F64vec& v){
        this->force += v;
    }
    inline void addPotIntra(const PS::F64& p){
        this->pot += p;
    }
    void clear(){
        this->clearForceIntra();
    }
};

//--- Full Particle class ---------------------------------------------------------------
//------ This class should contain all properties of the particle.
//------ This class is based on base classes.
class FP :
  public AtomID,
  public AtomType,
  public AtomConnect,
  public IntraMask,
  public AtomPos,
  public AtomVel,
  public AtomCharge,
  public AtomVDW,
  public AtomDipole,
  public ForceLJ,
  public ForceCoulomb,
  public ForceIntra{
  public:
    
    //--- recieve force calculation results
    //------ intramolecular force
    void copyFromForce(const ForceIntra & force) {
        this->addForceIntra(force.getForceIntra());
        this->addVirtIntra( force.getVirtIntra() );
        this->addPotIntra(  force.getPotIntra()  );
    }
    //------ coulomb force (permanent part, short range)
    void copyFromForce(const ForceCoulomb & force) {
        this->addFieldCoulomb( force.getFieldCoulomb() );
    }
    //------ LJ force
    void copyFromForce(const ForceLJ & force) {
        this->addForceLJ(force.getForceLJ());
        this->addVirtLJ( force.getVirtLJ() );
        this->addPotLJ(  force.getPotLJ()  );
    }
    //------ coulomb force from Particle Mesh (permanent part)
    void copyFromForceParticleMesh(const PS::F64vec & force_pm) {
        this->addFieldCoulomb( Normalize::realPMForce(force_pm) );
    }
    
    //--- output interaction result
    inline PS::F64vec getForce() const {
        return   this->getForceIntra()
               + this->getForceLJ()
               + this->getFieldCoulomb()*this->getCharge();
    }
    inline PS::F64vec getVirtCoulomb() const {
        PS::F64vec virt_tmp  = 0.0;
        PS::F64vec field_tmp = this->getFieldCoulomb();
        PS::F64vec pos_tmp   = Normalize::realPos( this->getPos() );
        virt_tmp[0] = field_tmp[0]*pos_tmp[0];
        virt_tmp[1] = field_tmp[1]*pos_tmp[1];
        virt_tmp[2] = field_tmp[2]*pos_tmp[2];
        virt_tmp = virt_tmp*this->getCharge();
        return virt_tmp;
    }
    
    //--- clear dinamic data
    void clear() {
        this->clearForceIntra();
        this->clearForceLJ();
        this->clearForceCoulomb();
    }
    
    //--- file I/O (overrided in "atsevb_FDPS_fileIO.hpp")
    void writeAscii(FILE * fp) const;
    void readAscii(FILE * fp);
};


//--- Essential Particle class ----------------------------------------------------------
//------ These classes have the meta-data for force calculation.
//------ They are the subset of Full Particle class.
//------ They are based on base classes.

//--------- for LJ interaction
class EPI_LJ :
  public AtomID,
  public IntraMask,
  public AtomPos,
  public AtomVDW{
  public:
    static PS::F64 rSearch;
    
    //--- member functions
    PS::F64 getRSearch() const {
        return this->rSearch;
    }
    
    //--- copy data from FP class
    void copyFromFP(const FP & fp) {
        copyAtomID(fp);
        copyIntraList(fp);
        copyAtomPos(fp);
        copyAtomVDW(fp);
    }
};
PS::F64 EPI_LJ::rSearch;

class EPJ_LJ :
  public AtomID,
  public AtomPos,
  public AtomVDW{
  public:
    static PS::F64 rSearch;
    
    //--- member functions
    PS::F64 getRSearch() const {
        return this->rSearch;
    }
    
    //--- copy data from FP class
    void copyFromFP(const FP & fp) {
        copyAtomID(fp);
        copyAtomPos(fp);
        copyAtomVDW(fp);
    }
};
PS::F64 EPJ_LJ::rSearch;

//--------- for coulomb interaction (permanent point charge)
class EPI_coulomb :
  public AtomID,
  public IntraMask,
  public AtomPos,
  public AtomCharge{
  public:
    static PS::F64 rSearch;
    
    //--- member functions
    PS::F64 getRSearch() const {
        return this->rSearch;
    }
    
    //--- copy data from FP class
    void copyFromFP(const FP & fp) {
        copyAtomID(fp);
        copyIntraList(fp);
        copyAtomPos(fp);
        copyAtomCharge(fp);
    }
};
PS::F64 EPI_coulomb::rSearch;

class EPJ_coulomb :
  public AtomID,
  public AtomPos,
  public AtomCharge{
  public:
    static PS::F64 rSearch;
    
    //--- member functions
    PS::F64 getRSearch() const {
        return this->rSearch;
    }
    
    //--- copy data from FP class
    void copyFromFP(const FP & fp) {
        copyAtomID(fp);
        copyAtomPos(fp);
        copyAtomCharge(fp);
    }
    //--- check same atom or not
    bool isSameAtom(const PS::S32 i_atom) const {
        return ( this->getAtomID() == i_atom );
    }
};
PS::F64 EPJ_coulomb::rSearch;

//class SPJ_coulomb {
//  private:
//    PS::F64    charge;
//    PS::F64vec pos;
//  public:
//    PS::F64vec getPos() const {
//        return this->pos;
//    }
//    void setPos(const PS::F64vec pos_new) {
//        this->pos = pos_new;
//    }
//    PS::F64    getCharge() const {
//        return this->charge;
//    }
//    void copyFromMoment(const PS::MomentMonopoleCutoff & mom) {
//        this->charge = mom.getCharge();
//        this->pos    = mom.getPos();
//    }
//    PS::MomentMonopoleCutoff convertToMoment() const {
//        return PS::MomentMonopoleCutoff(this->charge, this->pos);
//    }
//    void clear(){
//        this->charge = 0.0;
//        this->pos    = 0.0;
//    }
//    inline bool isSameAtom(const PS::S32 i_atom) const {
//        return false;
//    }
//};


//--------- for coulomb dipole
class EPI_dipole :
  public AtomID,
  public IntraMask,
  public AtomPos,
  public AtomDipole{
  private:
      
  public:
    static PS::F64 rSearch;
    
    PS::F64 getRSearch() const {
        return this->rSearch;
    }
  
    void copyFromFP(const FP & fp) {
        copyAtomID(fp);
        copyIntraList(fp);
        copyAtomPos(fp);
        copyAtomDipole(fp);
    }
};
PS::F64 EPI_dipole::rSearch;

class EPJ_dipole :
  public AtomID,
  public AtomPos,
  public AtomDipole{
  public:
    static PS::F64 rSearch;
    
    PS::F64 getRSearch() const {
        return this->rSearch;
    }
  
    void copyFromFP(const FP & fp) {
        copyAtomID(fp);
        copyAtomPos(fp);
        copyAtomDipole(fp);
    }
};
PS::F64 EPJ_dipole::rSearch;


//------ for intramolecular force
class EPI_intra :
  public AtomID,
  public AtomType,
  public AtomPos,
  public AtomConnect{
  private:
    //--- member valiables
    
  public:
    static PS::F64 rSearch;
    
    PS::F64 getRSearch() const {
        return this->rSearch;
    }
    
    //--- member functions
    void copyFromFP(const FP & fp) {
        copyAtomID(fp);
        copyAtomType(fp);
        copyAtomPos(fp);
        copyAtomConnect(fp);
    }
};
PS::F64 EPI_intra::rSearch;

class EPJ_intra :
  public EPI_intra {
  public:      
};

