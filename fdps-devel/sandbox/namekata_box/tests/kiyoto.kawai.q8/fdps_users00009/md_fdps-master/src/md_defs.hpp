/**************************************************************************************************/
/**
* @file  md_defs.hpp
* @brief definitions for data types, file name, ...etc.
*/
/**************************************************************************************************/
#pragma once

#include <vector>
#include <string>

#include <particle_simulator.hpp>

/**
* @brief definitions for data types, file name, ...etc.
*/
namespace MD_DEFS {

    const std::string condition_sequence_file{"condition_sequence.inp"};
    const std::string condition_molecule_file{"condition_molecule.inp"};

    const std::string model_dir{"./model/"};

    const std::string pos_data_dir{"./posdata"};
    const std::string resume_data_dir{"./resume"};
    const std::string VMD_data_dir{"./pdb"};

    constexpr size_t max_bond = 4;

    using ID_type  = PS::S64;

    struct IntraMask{
        ID_type id;
        PS::F32 scale_LJ;
        PS::F32 scale_coulomb;

    private:
        void _assign_impl(const ID_type id,
                          const PS::F32 scale_LJ,
                          const PS::F32 scale_coulomb){
            this->id            = id;
            this->scale_LJ      = scale_LJ;
            this->scale_coulomb = scale_coulomb;
        }
        void _assign_impl(const IntraMask &rhs){
            this->_assign_impl(rhs.id,
                               rhs.scale_LJ,
                               rhs.scale_coulomb);
        }

    public:
        IntraMask(){
            this->_assign_impl(-1, 1.0, 1.0);
        }
        IntraMask(const ID_type id,
                  const PS::F32 scale_LJ,
                  const PS::F32 scale_coulomb){
            this->_assign_impl(id,
                               scale_LJ,
                               scale_coulomb);
        }

        IntraMask& operator = (const IntraMask &rhs){
            this->_assign_impl(rhs);
            return *this;
        }

        IntraMask& setId(const ID_type id){
            this->id = id;
            return *this;
        }
    };
    bool operator == (const IntraMask &lhs, const IntraMask &rhs){
        return (lhs.id            == rhs.id &&
                lhs.scale_LJ      == rhs.scale_LJ &&
                lhs.scale_coulomb == rhs.scale_coulomb );
    }

    using AngleSet    = typename std::tuple<ID_type, ID_type, ID_type>;
    using TorsionSet  = typename std::tuple<ID_type, ID_type, ID_type, ID_type>;

    using MaskList    = typename std::vector<IntraMask>;
    using AngleList   = typename std::vector<AngleSet>;
    using TorsionList = typename std::vector<TorsionSet>;


    //--- mask interface
    //! @brief search interface for mask list.
    //! @return IntraMask.id = -1 means "not found". AtomID must be >= 0.
    IntraMask find_mask(const MaskList &list, const ID_type id){
        auto itr = std::find_if( list.begin(), list.end(),
                                 [id](const IntraMask &mask){ return mask.id == id; } );
        if(itr == list.end()){
            return IntraMask{-1, 1.0, 1.0};    // id = -1: means "not found", AtomID must be >= 0.
        } else {
            return *itr;
        }
    }
    //! @brief simple search interface for mask list.
    //! @return "true" means the id found in mask list.
    bool isFind_mask(const MaskList &list, const ID_type id){
        auto itr = std::find_if( list.begin(), list.end(),
                                 [id](const IntraMask &mask){ return mask.id == id; } );
        if(itr == list.end()){
            return false;
        } else {
            return true;
        }
    }
}


//--- output operators
inline std::ostream& operator << (std::ostream& s, const MD_DEFS::IntraMask &rhs){
    s << "(ID= "       << rhs.id
      << ", LJ= "      << rhs.scale_LJ
      << ", coulomb= " << rhs.scale_coulomb << ")";
    return s;
}
inline std::ostream& operator << (std::ostream& s, const MD_DEFS::AngleSet &rhs){
    s << " (" << std::get<0>(rhs)
      << ", " << std::get<1>(rhs)
      << ", " << std::get<2>(rhs) << ")";
    return s;
}
inline std::ostream& operator << (std::ostream& s, const MD_DEFS::TorsionSet &rhs){
    s << " (" << std::get<0>(rhs)
      << ", " << std::get<1>(rhs)
      << ", " << std::get<2>(rhs)
      << ", " << std::get<3>(rhs) << ")";
    return s;
}
