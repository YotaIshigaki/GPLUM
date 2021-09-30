//----------------------------------------------------------------------------------------
//  This file is enum class of coefficients for interaction.
//----------------------------------------------------------------------------------------
#pragma once

#include <string>
#include <stdexcept>


//--- enum indicator for intramolecular interaction
enum class IntraFuncForm : int {
    none,

    harmonic,
    anharmonic,

    cos,
    OPLS_3,
};

enum class TorsionShape : int {
    dihedral,
    improper,
};

namespace ENUM {

    //====================================
    //  enum interface for IntraFuncForm
    //====================================
    static const std::map<std::string, IntraFuncForm> table_str_IntraFuncForm{
        {"none"         , IntraFuncForm::none        },
        {"harmonic"     , IntraFuncForm::harmonic    },
        {"anharmonic"   , IntraFuncForm::anharmonic  },
        {"cos"          , IntraFuncForm::cos         },
        {"OPLS_3"       , IntraFuncForm::OPLS_3      },
    };

    static const std::map<IntraFuncForm, std::string> table_IntraFuncForm_str{
        {IntraFuncForm::none        , "none"         },
        {IntraFuncForm::harmonic    , "harmonic"     },
        {IntraFuncForm::anharmonic  , "anharmonic"   },
        {IntraFuncForm::cos         , "cos"          },
        {IntraFuncForm::OPLS_3      , "OPLS_3"       },
    };

    IntraFuncForm which_IntraFuncForm(const std::string &str){
        if(table_str_IntraFuncForm.find(str) != table_str_IntraFuncForm.end()){
            return table_str_IntraFuncForm.at(str);
        } else {
            std::cerr << "  IntraFuncForm: input = " << str << std::endl;
            throw std::out_of_range("undefined enum value in IntraFuncForm.");
        }
    }

    std::string what(const IntraFuncForm &e){
        if(table_IntraFuncForm_str.find(e) != table_IntraFuncForm_str.end()){
            return table_IntraFuncForm_str.at(e);
        } else {
            using type_base = typename std::underlying_type<IntraFuncForm>::type;
            std::cerr << "  IntraFuncForm: input = " << static_cast<type_base>(e) << std::endl;
            throw std::out_of_range("undefined enum value in IntraFuncForm.");
        }
    }

    //====================================
    //  enum interface for TorsionShape
    //====================================
    static const std::map<std::string, TorsionShape> table_str_TorsionShape{
        {"dihedral", TorsionShape::dihedral},
        {"improper", TorsionShape::improper},
    };

    static const std::map<TorsionShape, std::string> table_TorsionShape_str{
        {TorsionShape::dihedral, "dihedral"},
        {TorsionShape::improper, "improper"},
    };

    TorsionShape which_TorsionShape(const std::string &str){
        if(table_str_TorsionShape.find(str) != table_str_TorsionShape.end()){
            return table_str_TorsionShape.at(str);
        } else {
            std::cerr << "  TorsionShape: input = " << str << std::endl;
            throw std::out_of_range("undefined enum value in TorsionShape.");
        }
    }

    std::string what(const TorsionShape &e){
        if(table_TorsionShape_str.find(e) != table_TorsionShape_str.end()){
            return table_TorsionShape_str.at(e);
        } else {
            using type_base = typename std::underlying_type<TorsionShape>::type;
            std::cerr << "  TorsionShape: input = " << static_cast<type_base>(e) << std::endl;
            throw std::out_of_range("undefined enum value in TorsionShape.");
        }
    }

}

//--- specialize for std::to_string()
namespace std {
    inline string to_string(const IntraFuncForm &e){ return ENUM::what(e); }
    inline string to_string(const TorsionShape  &e){ return ENUM::what(e); }
}

//--- output function as "std::cout << (enum class::value)"
inline std::ostream& operator << (std::ostream& s, const IntraFuncForm &e){
    s << ENUM::what(e);
    return s;
}

inline std::ostream& operator << (std::ostream& s, const TorsionShape &e){
    s << ENUM::what(e);
    return s;
}

//--- specialized hash function
//        (support for gcc 4.x ~ 5.x. naturally supported by gcc 6.0 or later.)
//         ref: http://qiita.com/taskie/items/479d649ea1b20bacbe03
namespace std {
    template <>
    struct hash<IntraFuncForm> {
        size_t operator() (IntraFuncForm x) const noexcept {
            using type = typename underlying_type<IntraFuncForm>::type;
            return hash<type>{}(static_cast<type>(x));
        }
    };

    template <>
    struct hash<TorsionShape> {
        size_t operator() (TorsionShape x) const noexcept {
            using type = typename underlying_type<TorsionShape>::type;
            return hash<type>{}(static_cast<type>(x));
        }
    };

}
