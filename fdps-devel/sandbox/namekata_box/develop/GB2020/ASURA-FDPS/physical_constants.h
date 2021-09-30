#pragma once

namespace physical_constants {

    // Physical constants
    extern const double Ggrav;  // Gravitational constants
    extern const double kBoltz; // Boltzmann constant 
    
    // Mass units
    extern const double dalton; // Unified atomic mass unit
    extern const double Mproton; // proton mass
    extern const double Mhydrogen; // 1H + 2H
    extern const double Mhelium; // 3He + 4He
    extern const double Mcarbon; // 12C + 13C
    extern const double Mnitrogen; // 14N + 15N
    extern const double Moxygen; // 16O + 17O + 18O
    extern const double Mneon; // 20Ne + 21Ne + 22Ne
    extern const double Mmagnesium; // 24Mg + 25Mg + 26Mg
    extern const double Msilicon; // 28Si + 29Si + 30Si
    extern const double Msulfur; // 32S + 33S + 34S + 36S
    extern const double Mcalcium; // 40Ca + 42Ca + 43Ca + 44Ca + 46Ca + 48Ca
    extern const double Miron; // 54Fe + 56Fe + 57Fe + 58Fe
    extern const double Mnickel; // 58Ni + 60Ni + 61Ni + 62Ni + 64Ni
    extern const double Meuropium; // 151Eu + 153Eu
    extern const double Msolar; // Solar mass
    
    // Length units
    extern const double km; // kilo-meters
    extern const double AU; // Astronomical unit
    extern const double pc, kpc, Mpc, Gpc; // parsec
    
    // Time units
    extern const double yr, kyr, Myr, Gyr; // year

}
namespace phys_const = physical_constants;
