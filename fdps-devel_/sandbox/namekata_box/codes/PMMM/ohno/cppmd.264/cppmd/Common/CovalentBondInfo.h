#ifndef COVALENTBONDINFO_H
#define COVALENTBONDINFO_H

#include <algorithm>
#include <vector>
#include <set>
#include "Common.h"
#include "ParticleInfo.h" 

namespace CovalentBondInfo{
extern int bond_assign_offset;
extern int angle_assign_offset;
extern int torsion_assign_offset;
extern int improper_assign_offset;

  struct Bond{
    union{
      int id_of_atom[2];
      int index_of_atom[2];
    };
    int typeofbond;
    int shake;                 // 0:not shake bond, 1: shake bond
    inline int get_assign_atomid() const
    {
      return id_of_atom[bond_assign_offset];
    }
  };
  struct BondParameter {
    double force_constant;      // RK: force constant, kcal/mol
    double equilibrium_length;  // REQ: the equilibrium bond length, angstroms
  };

  struct Angle{
    union{
      int id_of_atom[3];
      int index_of_atom[3];
    };
    int typeofangle;
    inline int get_assign_atomid() const
    {
      return id_of_atom[angle_assign_offset];
    }
  };
  struct AngleParameter {
    double force_constant;      // TK: force constant, kcal/mol A**2
    double equilibrium_angle;   // TEQ: the equilibrium angle, radians
  };
  
  struct Torsion{
    union{
      int id_of_atom[4];
      int index_of_atom[4];
    };
    int typeoftorsion;
    bool calc14interaction;
    inline int get_assign_atomid() const
    {
      return id_of_atom[torsion_assign_offset];
    }
  };
  struct TorsionParameter {
    double force_constant;      // PK: force constant, kcal/mol
    double periodicity;         // PN: periodicity of the dihedral
    double phase;               // PHASE: phase of the dihedral, radians
  };

  struct Improper{
    union{
      int id_of_atom[4];
      int index_of_atom[4];
    };
    int typeofimproper;
    inline int get_assign_atomid() const
    {
      return id_of_atom[improper_assign_offset];
    }
  };
  struct ImproperParameter {
    double force_constant;      // PK: force constant, kcal/mol
    double periodicity;         // PN: periodicity of the dihedral
    double phase;               // PHASE: phase of the dihedral, radians
  };


//! 
  struct BondList {
    std::vector<Bond> BondArray;
    std::vector<Angle> AngleArray;
    std::vector<Torsion> TorsionArray;
    std::vector<Improper> ImproperArray;

    BondList();
    BondList(int num_bond, int num_angle, 
             int num_torsion, int num_improper);
    BondList(const BondList& bl);
    
    void clear();
    void print();
    void nozero_print();
    //! calc size of packed int array
    /*!
      @return size of packed array in BYTE
      @note for estimate required size of MPI buffer
     */
    int size_of_packed_int() const;
    //! pack array of int
    /*!
      @param[out] bondpack int array which must be enough large
      @note for write to MPI buffer
      @note use size_of_packed_int() to get required size of bondpack 
     */
    int pack_int_array(int* bondpack) const;
    //! unpack array of int
    /*!
      @note for read from MPI buffer
     */
    void unpack_int_array(int bondpack[]);
    void push_back(BondList addlist);
    //! @return size in BYTE
    size_t size();
    template<typename CB>
    void pick_up(std::vector<CB>& destarray,
                 const std::vector<CB>& srcarray, const AtomID atomid)
    {
#if 0
      // g++ 4.4.5 does not take the following line. ???
      for(std::vector<CB>::size_type bi=0;bi<srcarray.size();bi++){
#else
      for(size_t bi=0;bi<srcarray.size();bi++){
#endif
        if(atomid==srcarray[bi].get_assign_atomid()){
          destarray.push_back(srcarray[bi]);
        }
      }
    }
    void pick_up_bondlist(CovalentBondInfo::BondList& bondlist,
                          const AtomID atomid);
    template<typename CB>
    class MachAtomID
    {
     public:
      MachAtomID(AtomID atomid){aid=atomid;}
      bool operator()(CB cb) const {return cb.get_assign_atomid()==aid;}
     private:
      AtomID aid;
    };
    template<typename CB>
    class MachAtomIDIndex
    {
     public:
      MachAtomIDIndex(AtomID atomid, int index){aid=atomid;idx=index;}
      bool operator()(CB cb) const {return cb.id_of_atom[idx]==aid;}
     private:
      AtomID aid;
      int idx;
    };
    void remove_atomid(const AtomID atomid);
  };

  struct BondParameterList {
    std::vector<BondParameter> BondParameterArray;
    std::vector<AngleParameter> AngleParameterArray;
    std::vector<TorsionParameter> TorsionParameterArray;
    std::vector<ImproperParameter> ImproperParameterArray;
  };

inline int bond_assign_atom(int assign_type)
{
  if(assign_type<0)assign_type=0;
  if(assign_type>1)assign_type=1;
  return assign_type;
}
inline int angle_assign_atom(int assign_type)
{
  if(assign_type<0)assign_type=0;
  if(assign_type>1)assign_type=1;
  return assign_type;
}
inline int torsion_assign_atom(int assign_type)
{
  if(assign_type<0)assign_type=0;
  if(assign_type>1)assign_type=1;
  return assign_type;
}
inline int improper_assign_atom(int assign_type)
{
  if(assign_type<0)assign_type=0;
  if(assign_type>2)assign_type=2;
  return assign_type;
}
}

/// Information of Covalent Bonds
class CovalentBondList{
public:
  std::vector<CovalentBondInfo::Bond> bond;
  std::vector<CovalentBondInfo::Angle> angle;
  std::vector<CovalentBondInfo::Torsion> torsion;
  std::vector<CovalentBondInfo::Improper> improper;
  std::set<int> atomidset;

  CovalentBondList();

  explicit CovalentBondList(int num);

  ~CovalentBondList();

  void clear();

  void merge_bondlistarray(const std::vector<CovalentBondInfo::BondList>& bondlistarray);
  
  int make_atomidset(){
    atomidset.clear();
    for(size_t b=0;b<bond.size();b++){
      atomidset.insert(bond[b].id_of_atom[0]);
      atomidset.insert(bond[b].id_of_atom[1]);
    }
    for(size_t b=0;b<angle.size();b++){
      atomidset.insert(angle[b].id_of_atom[0]);
      atomidset.insert(angle[b].id_of_atom[1]);
      atomidset.insert(angle[b].id_of_atom[2]);
    }
    for(size_t b=0;b<torsion.size();b++){
      atomidset.insert(torsion[b].id_of_atom[0]);
      atomidset.insert(torsion[b].id_of_atom[1]);
      atomidset.insert(torsion[b].id_of_atom[2]);
      atomidset.insert(torsion[b].id_of_atom[3]);
    }
    for(size_t b=0;b<improper.size();b++){
      atomidset.insert(improper[b].id_of_atom[0]);
      atomidset.insert(improper[b].id_of_atom[1]);
      atomidset.insert(improper[b].id_of_atom[2]);
      atomidset.insert(improper[b].id_of_atom[3]);
    }

    return atomidset.size();
  }

  void pick_up_bondlist(CovalentBondInfo::BondList& bondlist,
                        const ParticleArray& particlearray,
                        const TypeRange& typerange,
                        int assign_type);

  void print();
  
  int size_of_packed_int();
  //! pack array of int
  /*!
    @param[out] bondpack int array which must be enough large
    @note for write to MPI buffer
    @note use size_of_packed_int() to get required size of bondpack 
  */
  int pack_int_array(int* bondpack);
  //! unpack array of int
  /*!
    @note for read from MPI buffer
  */
  void unpack_int_array(int bondpack[]);

  CovalentBondList& append(const CovalentBondList& cbl);
  CovalentBondList& swap(CovalentBondList& cbl);

};

/// Information of Covalent Bonds
class CovalentBondParameterList{
public:
  std::vector<CovalentBondInfo::BondParameter> bond;
  std::vector<CovalentBondInfo::AngleParameter> angle;
  std::vector<CovalentBondInfo::TorsionParameter> torsion;
  std::vector<CovalentBondInfo::ImproperParameter> improper;

  CovalentBondParameterList();

  CovalentBondParameterList(int num);

  ~CovalentBondParameterList();

  void clear();

  void print() const;

  int size_of_packed_double(int nums[4]);
  //! pack array of int
  /*!
    @param[out] bondpack int array which must be enough large
    @note for write to MPI buffer
    @note use size_of_packed_int() to get required size of bondpack 
  */
  int pack_double_array(int nums[4], double* bondpack);
  //! unpack array of int
  /*!
    @note for read from MPI buffer
  */
  void unpack_double_array(int nums[4], double bondpack[]);
};

#endif

