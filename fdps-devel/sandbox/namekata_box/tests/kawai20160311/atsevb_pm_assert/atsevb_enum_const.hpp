//***************************************************************************************
//  This program is the definition of enumulator constant.
//    This code is used by "atsevb_main.cpp"
//***************************************************************************************
#pragma once

enum class MOL_TYPE : int {
    polymer,
    wall,
    solvent,
    cluster,
    
    test,
    
    ENUM_N_TOT    // this tag must be at last.
};

enum class ATOM_TYPE : int {
    H_water,
    O_water,
    
    H_oxo,
    O_oxo,
    
    O_so3,
    S_so3,
    C_sc,
    F_sc,
    
    test,
    
    ENUM_N_TOT    // this tag must be at last.
};

enum class IO_MODE : int {
    pos,
    resume,
    VMD,
};

namespace ENUM_TAG {
    std::array<std::string, static_cast<int>(ATOM_TYPE::ENUM_N_TOT)> name_tag;
    std::array<std::string, static_cast<int>(ATOM_TYPE::ENUM_N_TOT)> residue_tag;
    
    void Init(){
        //--- set name tag
        name_tag.at( static_cast<int>(ATOM_TYPE::H_water) ) = "HW ";
        name_tag.at( static_cast<int>(ATOM_TYPE::O_water) ) = "OW ";
        name_tag.at( static_cast<int>(ATOM_TYPE::H_oxo) )   = "HH ";
        name_tag.at( static_cast<int>(ATOM_TYPE::O_oxo) )   = "OH ";
        
        name_tag.at( static_cast<int>(ATOM_TYPE::test) )   = "dum";
        
        
        //--- set region tag
        residue_tag.at( static_cast<int>(ATOM_TYPE::H_water) ) = "WAT";
        residue_tag.at( static_cast<int>(ATOM_TYPE::O_water) ) = "WAT";
        residue_tag.at( static_cast<int>(ATOM_TYPE::H_oxo) )   = "OXO";
        residue_tag.at( static_cast<int>(ATOM_TYPE::O_oxo) )   = "OXO";
        
        residue_tag.at( static_cast<int>(ATOM_TYPE::test) )   = "dum";
    }
}

