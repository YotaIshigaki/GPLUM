//***************************************************************************************
//  This code is the FullParticle class and EssentialParticle class.
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include <sstream>
#include <tuple>
#include <unordered_map>
#include <stdexcept>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "unit.hpp"
#include "md_defs.hpp"
#include "md_enum.hpp"
#include "atom_class_base.hpp"
#include "file_IO_pos.hpp"
#include "file_IO_resume.hpp"
#include "file_Out_VMD.hpp"
#include "md_coef_table.hpp"


//------ connection
class AtomIntraMask :
  public AtomID {
  public:
    //--- static data for intra mask
    static std::unordered_map< MD_DEFS::ID_type,
                               MD_DEFS::MaskList > intra_mask_table;

    template <class Tptcl>
    void copyAtomIntraMask(const Tptcl &fp){
        this->copyAtomID(fp);
    }

    MD_DEFS::MaskList&       mask_list()       { return this->intra_mask_table[this->getId()]; }
    const MD_DEFS::MaskList& mask_list() const { return this->intra_mask_table[this->getId()]; }

    inline MD_DEFS::IntraMask find_mask(const MD_DEFS::ID_type id_j) const {
        return MD_DEFS::find_mask(this->mask_list(), id_j);
    }
    inline bool isFind_mask(const MD_DEFS::ID_type id_j) const {
        return MD_DEFS::isFind_mask(this->mask_list(), id_j);
    }
};
std::unordered_map< MD_DEFS::ID_type,
                    MD_DEFS::MaskList > AtomIntraMask::intra_mask_table;

class AtomConnect :
  public AtomIntraMask {
  public:
    MD_EXT::basic_connect<MD_DEFS::ID_type,
                          MD_DEFS::max_bond> bond;

    template <class Tptcl>
    void copyAtomConnect(const Tptcl &fp){
        this->copyAtomIntraMask(fp);
        this->bond = fp.bond;
    }

    void shift_pair_ID(const MD_DEFS::ID_type shift){
        auto bond_tmp = this->bond;
        this->bond.clear();
        for(const auto& b : bond_tmp){
            this->bond.add(shift + b);
        }
    }

    //--- static data for intra pair table
    static std::unordered_map< MD_DEFS::ID_type,
                               std::tuple<MD_DEFS::AngleList,
                                          MD_DEFS::TorsionList,
                                          MD_DEFS::TorsionList> > intra_pair_table;

    MD_DEFS::AngleList&   angle_list()   { return std::get<0>(this->intra_pair_table[this->getId()]); }
    MD_DEFS::TorsionList& dihedral_list(){ return std::get<1>(this->intra_pair_table[this->getId()]); }
    MD_DEFS::TorsionList& improper_list(){ return std::get<2>(this->intra_pair_table[this->getId()]); }

    const MD_DEFS::AngleList&   angle_list()    const { return std::get<0>(this->intra_pair_table[this->getId()]); }
    const MD_DEFS::TorsionList& dihedral_list() const { return std::get<1>(this->intra_pair_table[this->getId()]); }
    const MD_DEFS::TorsionList& improper_list() const { return std::get<2>(this->intra_pair_table[this->getId()]); }

    void clear_intra_list(){
        this->mask_list().clear();
        this->angle_list().clear();
        this->dihedral_list().clear();
        this->improper_list().clear();
    }
};
std::unordered_map< MD_DEFS::ID_type,
                    std::tuple<MD_DEFS::AngleList,
                               MD_DEFS::TorsionList,
                               MD_DEFS::TorsionList> > AtomConnect::intra_pair_table;


//--- derived force class
template <class Tf>
class ForceInter :
  public ForceLJ<Tf>,
  public ForceCoulomb<Tf> {
    public:

      //--- copy results of calculated force
      template <class Tforce>
      void copyFromForce(const Tforce &force){
          this->copyForceLJ(force);
          this->copyForceCoulomb(force);
      }

      //--- clear dinamic data
      void clearForceInter(){
          this->clearForceLJ();
          this->clearForceCoulomb();
      }
      void clear() {
          this->clearForceInter();
      }
};

template <class Tf>
class Force_FP :
  public ForceInter<Tf>,
  public ForceIntra<Tf> {
  public:

    //--- copy results of calculated force
    template <class Tforce>
    void copyFromForce(const Tforce &force){
        this->copyForceInter(force);
        this->copyForceIntra(force);
    }
    template <class Tf_rhs>
    void copyFromForce(const ForceInter<Tf_rhs> &force){
        this->copyForceInter(force);
    }
    template <class Tf_rhs>
    void copyFromForce(const ForceIntra<Tf_rhs> &force){
        this->copyForceIntra(force);
    }

    //--- clear dinamic data
    void clear() {
        this->clearForceInter();
        this->clearForceIntra();
    }
};


//--- Full Particle class ---------------------------------------------------------------
//------ This class should contain all properties of the particle.
//------ This class is based on base classes.
class Atom_FP :
  public AtomType,
  public AtomConnect,
  public AtomPos<PS::F32>,
  public AtomVel<PS::F32>,
  public AtomCharge<PS::F32>,
  public AtomVDW<PS::F32>,
  public Force_FP<PS::F32> {
  public:

    //--- output interaction result
    inline PS::F32vec getForce() const {
        return   this->getForceIntra()
               + this->getForceLJ()
               + this->getFieldCoulomb()*this->getCharge();
    }
    inline PS::F32vec getVirial() const {
        return   this->getVirialIntra()
               + this->getVirialLJ()
               + this->getPotCoulomb()*this->getCharge();
    }

    //--- copy model property from molecular model template (using FP class)
    void copyFromModelTemplate(const PS::S64 &atom_id_shift,
                               const PS::S64 &mol_id,
                               const Atom_FP &tmp){
        if(atom_id_shift < 0 || mol_id < 0){
            std::ostringstream oss;
            oss << "atom_id_shift & mol_id are must be > 0." << "\n"
                << "  atom_id_shift = " << atom_id_shift << "\n"
                << "  mol_id        = " << mol_id        << "\n";
            throw std::invalid_argument(oss.str());
        }
        *this = tmp;

        //--- id
        this->setAtomID(atom_id_shift + tmp.getAtomID());
        this->setMolID(mol_id);

        //--- connection pair
        this->shift_pair_ID(atom_id_shift);
    }

    std::string str() const {
        std::ostringstream oss;

        oss << "    AtomID   : " << this->getAtomID()   << "\n";
        oss << "    MolID    : " << this->getMolID()    << "\n";
        oss << "    AtomType : " << this->getAtomType() << "\n";
        oss << "    MolType  : " << this->getMolType()  << "\n";
        oss << "    Pos      : ("  << this->getPos().x
                           << ", " << this->getPos().y
                           << ", " << this->getPos().z << ")\n";
        oss << "    Mass     : "  << std::setw(12) << this->getMass()
            << "  | real value: " << std::setw(12) << this->getMass()*Unit::mass_C         << " [kg/mol]" << "\n";
        oss << "    charge   : "  << std::setw(12) << this->getCharge()
            << "  | real value: " << std::setw(12) << this->getCharge()/Unit::coef_coulomb << " [electron charge]" << "\n";
        oss << "    VDW_D    : "  << std::setw(12) << this->getVDW_D()
            << "  | real value: " << std::setw(12) << this->getVDW_D()*this->getVDW_D()    << " [kcal/mol]" << "\n";
        oss << "    VDW_R    : "  << std::setw(12) << this->getVDW_R()
            << "  | real value: " << std::setw(12) << this->getVDW_R()*2.0                 << " [angstrom]" << "\n";

        oss << "    bond: n=" << this->bond.size() << "\n";
        oss << "       id= ";
        bool isFirst = true;
        for(const auto& b : bond){
            if(!isFirst) oss << ", ";
            isFirst = false;

            oss << " " << b;
        }
        oss << "\n";

        return oss.str();
    }

    //--- interface for file I/O
    //------ for pos
    void write_pos_ascii(FILE *fp) const { FILE_IO::Pos::write_atom(fp, *this); }
    void read_pos_ascii( FILE *fp)       { FILE_IO::Pos::read_atom( fp, *this); }

    //------ for resume
    void write_resume_ascii(FILE *fp) const { FILE_IO::Resume::write_atom(fp, *this); }
    void read_resume_ascii( FILE *fp)       { FILE_IO::Resume::read_atom( fp, *this); }

    //------ for VMD
    std::string getResidue() const {
        const std::string res = MODEL::coef_table.residue.at( std::make_tuple( this->getMolType(),
                                                                               this->getAtomType() ) );
        return res.substr(0,3);  // length <= 3.
    }
    void write_VMD_ascii(FILE *fp) const { FILE_IO::VMD::write_atom(fp, *this); }
};

std::ostream& operator << (std::ostream &s, const Atom_FP &atom){
    s << atom.str();
    return s;
}


//--- Essential Particle class ----------------------------------------------------------
//------ These classes have the meta-data for force calculation.
//------ They are the subset of Full Particle class.
//------ They are based on base classes.

class EP_inter :
  public AtomIntraMask,
  public AtomType,
  public AtomPos<PS::F32>,
  public AtomCharge<PS::F32>,
  public AtomVDW<PS::F32> {
  private:
  //public:
    static PS::F32 Rcut;
    static PS::F32 Rcut_LJ;
    static PS::F32 Rcut_coulomb;

  public:
    static void setRcut_LJ(const PS::F32 &r){
        EP_inter::Rcut_LJ = r;
        EP_inter::Rcut    = std::max(r, EP_inter::Rcut_coulomb);
    }
    static void setRcut_coulomb(const PS::F32 &r){
        EP_inter::Rcut_coulomb = r;
        EP_inter::Rcut         = std::max(r, EP_inter::Rcut_LJ);
    }

    static PS::F32 getRSearch()      { return EP_inter::Rcut; }
    static PS::F32 getRcut_LJ()      { return EP_inter::Rcut_LJ; }
    static PS::F32 getRcut_coulomb() { return EP_inter::Rcut_coulomb; }

    template <class T>
    void copyFromFP(const T &fp){
        this->copyAtomIntraMask(fp);
        this->copyAtomType(fp);
        this->copyAtomPos(fp);
        this->copyAtomCharge(fp);
        this->copyAtomVDW(fp);
    }
};
PS::F32 EP_inter::Rcut         = 0.0;
PS::F32 EP_inter::Rcut_LJ      = 0.0;
PS::F32 EP_inter::Rcut_coulomb = 0.0;


class EP_intra :
  public AtomType,
  public AtomConnect,
  public AtomPos<PS::F32> {
  private:
    static PS::F32 Rcut;

  public:
    static void    setRcut(const PS::F32 r) { EP_intra::Rcut = r;    }
    static PS::F32 getRSearch()             { return EP_intra::Rcut; }

    template <class T>
    void copyFromFP(const T &fp){
        this->copyAtomType(fp);
        this->copyAtomConnect(fp);
        this->copyAtomPos(fp);
    }
};
PS::F32 EP_intra::Rcut = 0.0;
