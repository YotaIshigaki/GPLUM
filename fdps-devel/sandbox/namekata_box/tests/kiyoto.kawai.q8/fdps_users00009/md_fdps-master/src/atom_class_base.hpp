//***************************************************************************************
//  This code is base class for FullParticle and EssentialParticle.
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include <vector>
#include <string>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "md_enum.hpp"
#include "md_defs.hpp"


//--- base classes for particle data ----------------------------------------------------
//------ identifier
class AtomID{
protected:
    MD_DEFS::ID_type atom_id = -1;  // atom id
    MD_DEFS::ID_type mol_id  = -1;  // molecular id

public:
    //--- member functions
    void setAtomID(const MD_DEFS::ID_type &id){ this->atom_id = id; }
    void setMolID( const MD_DEFS::ID_type &id){ this->mol_id  = id; }
    inline MD_DEFS::ID_type getAtomID() const { return this->atom_id; }
    inline MD_DEFS::ID_type getMolID()  const { return this->mol_id;  }

    inline MD_DEFS::ID_type getId()     const { return this->atom_id; }

    //--- copy function
    template<class Tptcl>
    void copyAtomID(const Tptcl &fp){
        this->atom_id  = fp.getAtomID();
        this->mol_id   = fp.getMolID();
    }
};

//------ sub identifier
class AtomType{
protected:
    AtomName atom_type = static_cast<AtomName>(-1);  // illigal value
    MolName  mol_type  = static_cast<MolName>(-1);

public:
    void setAtomType(const AtomName &type){ this->atom_type = type; }
    void setAtomType(const std::string &str){
        this->setAtomType(ENUM::which_AtomName(str));
    }
    void setMolType(const MolName &type){ this->mol_type = type; }
    void setMolType(const std::string &str){
        this->setMolType(ENUM::which_MolName(str));
    }

    inline AtomName getAtomType() const { return this->atom_type; }
    inline MolName  getMolType()  const { return this->mol_type; }

    template<class Tptcl>
    void copyAtomType(const Tptcl &fp){
        this->atom_type = fp.getAtomType();
        this->mol_type  = fp.getMolType();
    }
};


//------ position
template <class Tf>
class AtomPos{
protected:
    PS::Vector3<Tf> pos = 0.0;

public:

    //--- member functions
    void setPos(const PS::Vector3<Tf> &pos_new){ this->pos = pos_new; }
    inline PS::Vector3<Tf> getPos() const { return this->pos; }

    //--- copy function
    template<class Tptcl>
    void copyAtomPos(const Tptcl &fp){
        this->pos = fp.getPos();
    }
};

//------ mass & velocity
template <class Tf>
class AtomVel{
protected:
    Tf              mass = 0.0;  // mass of atom
    PS::Vector3<Tf> vel  = 0.0;  // velocity of atom
    PS::Vector3<Tf> trj  = 0.0;  // trajectory of atom

public:
    void setMass(const Tf              &mass_new){ this->mass = mass_new; }
    void setVel( const PS::Vector3<Tf> &vel_new){  this->vel  = vel_new; }

    inline Tf              getMass() const { return this->mass; }
    inline PS::Vector3<Tf> getVel()  const { return this->vel; }

    inline void            clearTrj()                                { this->trj = 0.0;   }
    inline void            addTrj(const PS::Vector3<Tf> &move)       { this->trj += move; }
    inline PS::Vector3<Tf> getTrj()                            const { return this->trj;  }

    inline PS::F64 getKinetic() const {
        PS::F64 v2 = (this->vel*this->vel);
        return 0.5*v2/this->mass;
    }
};

//------ permanent charge
template <class Tf>
class AtomCharge{
protected:
    Tf charge = 0.0;

public:

    void setCharge(const Tf &q){ this->charge = q; }
    inline Tf getCharge()             const { return this->charge; }
    inline Tf getChargeParticleMesh() const { return this->charge; }

    template<class Tptcl>
    void copyAtomCharge(const Tptcl &fp){
        this->charge = fp.getCharge();
    }
};

//------ VDW parameters
template <class Tf>
class AtomVDW{
protected:
    Tf vdw_r = 0.0;
    Tf vdw_d = 0.0;

public:

    void setVDW_R(const Tf &vdw_r){ this->vdw_r = vdw_r; }
    void setVDW_D(const Tf &vdw_d){ this->vdw_d = vdw_d; }
    inline Tf getVDW_R() const { return this->vdw_r; }
    inline Tf getVDW_D() const { return this->vdw_d; }

    template<class Tptcl>
    void copyAtomVDW(const Tptcl &fp){
        this->vdw_r = fp.getVDW_R();
        this->vdw_d = fp.getVDW_D();
    }
};


//--- Force class -----------------------------------------------------------------------
//------ This class has the result of force. (used as temporary data)
//------ This class is the subset of Full Particle class.

//------ LJ interaction with cutoff length (simple cutoff)
template <class Tf>
class ForceLJ {
protected:
    PS::Vector3<Tf> force_LJ  = 0.0;
    PS::Vector3<Tf> virial_LJ = 0.0;
    Tf              pot_LJ    = 0.0;

public:
    void clearForceLJ(){
        this->force_LJ  = 0.0;
        this->virial_LJ = 0.0;
        this->pot_LJ    = 0.0;
    }
    void clear(){ this->clearForceLJ(); }

    inline PS::Vector3<Tf> getForceLJ()  const { return this->force_LJ;  }
    inline PS::Vector3<Tf> getVirialLJ() const { return this->virial_LJ; }
    inline Tf              getPotLJ()    const { return this->pot_LJ;    }
    inline void addForceLJ( const PS::Vector3<Tf> &force)  { this->force_LJ  += force;  }
    inline void addVirialLJ(const PS::Vector3<Tf> &virial) { this->virial_LJ += virial; }
    inline void addPotLJ(   const Tf              &pot)    { this->pot_LJ    += pot;    }

    template <class T>
    void copyForceLJ(const T &f){
        this->force_LJ  = f.getForceLJ();
        this->virial_LJ = f.getVirialLJ();
        this->pot_LJ    = f.getPotLJ();
    }
    template <class Tf_rhs>
    void copyFromForce(const ForceLJ<Tf_rhs> &f){
        this->copyForceLJ(f);
    }
};

//------ coulomb interaction
template <class Tf>
class ForceCoulomb {
protected:
    PS::Vector3<Tf> field_coulomb = 0.0;
    Tf              pot_coulomb   = 0.0;

public:
    void clearForceCoulomb(){
        field_coulomb = 0.0;
          pot_coulomb = 0.0;
    }
    void clear(){ this->clearForceCoulomb(); }

    inline PS::Vector3<Tf> getFieldCoulomb() const { return this->field_coulomb; }
    inline Tf              getPotCoulomb()   const { return this->pot_coulomb;   }
    inline void addFieldCoulomb(const PS::Vector3<Tf> &f){ this->field_coulomb += f; }
    inline void addPotCoulomb(  const Tf              &p){ this->pot_coulomb   += p; }

    template <class T>
    void copyForceCoulomb(const T &f){
        this->field_coulomb  = f.getFieldCoulomb();
        this->pot_coulomb    = f.getPotCoulomb();
    }
    template <class Tf_rhs>
    void copyFromForce(const ForceCoulomb<Tf_rhs> &f){
        this->copyForceCoulomb(f);
    }
};

//--- Intramoleecular interaction
template <class Tf>
class ForceIntra {
protected:
    PS::Vector3<Tf> force_intra  = 0.0;
    PS::Vector3<Tf> virial_intra = 0.0;
    Tf              pot_bond     = 0.0;
    Tf              pot_angle    = 0.0;
    Tf              pot_torsion  = 0.0;

public:
    void clearForceIntra(){
        this->force_intra  = 0.0;
        this->virial_intra = 0.0;
        this->pot_bond     = 0.0;
        this->pot_angle    = 0.0;
        this->pot_torsion  = 0.0;
    }
    void clear(){ this->clearForceIntra(); }

    inline PS::Vector3<Tf> getForceIntra()  const { return this->force_intra;  }
    inline PS::Vector3<Tf> getVirialIntra() const { return this->virial_intra; }
    inline Tf              getPotBond()     const { return this->pot_bond;     }
    inline Tf              getPotAngle()    const { return this->pot_angle;    }
    inline Tf              getPotTorsion()  const { return this->pot_torsion;  }
    inline void addForceIntra( const PS::Vector3<Tf> &f){ this->force_intra  += f; }
    inline void addVirialIntra(const PS::Vector3<Tf> &v){ this->virial_intra += v; }
    inline void addPotBond(    const Tf &p){ this->pot_bond    += p; }
    inline void addPotAngle(   const Tf &p){ this->pot_angle   += p; }
    inline void addPotTorsion( const Tf &p){ this->pot_torsion += p; }

    template <class T>
    void copyForceIntra(const T &f){
        this->force_intra  = f.getForceIntra();
        this->virial_intra = f.getVirialIntra();
        this->pot_bond     = f.getPotBond();
        this->pot_angle    = f.getPotAngle();
        this->pot_torsion  = f.getPotTorsion();
    }
    template <class Tf_rhs>
    void copyFromForce(const ForceIntra<Tf_rhs> &f){
        this->copyForceIntra(f);
    }
};
