//***************************************************************************************
//  This program is the definition of enumulator constant.
//    This code is used by "atsevb_main.cpp"
//***************************************************************************************
#pragma once

#include<cmath>

namespace Unit {
    
    //--- fundamental parameters
    static constexpr PS::F64 pi    = 3.141592653589793;  // circular constant
    static constexpr PS::F64 angm  = 1.0e-10;            // [m] [angstrom] to [meter]
    
    //--- ref: The NIST reference (http://physics.nist.gov/cuu/Constants/index.html)
    static constexpr PS::F64 ec    = 1.6021766208e-19;   // [C] charge strength of electron
    static constexpr PS::F64 e0    = 8.854187817620e-12; // [C^2/J·m] permittivity of vacuum
    static constexpr PS::F64 avgdr = 6.022140857e+23;    // [/mol] Avogadro constant
    static constexpr PS::F64 blz   = 1.38064852e-23;     // [J/K] Boltzmann constant
    
    //--- ref: The elements, Theodore gray, 2010
    static constexpr PS::F64 mass_C        = 12.0107e-3;  // [kg/mol] mass of Carbon atom
    //--- ref: JSME text series "Thermodynamics"
    static constexpr PS::F64 mech_heat_equ = 4.1868e+3;   // [J/kcal] mechanical equivalent of heat
    
    //--- standard values                                               // J = [kg·m^2 / s^2]
    static constexpr PS::F64 norm_energy = mech_heat_equ/avgdr;         // [mol·J/kcal]
    static constexpr PS::F64 norm_length = angm;                        // [m]
    static constexpr PS::F64 norm_mass   = mass_C/avgdr;                // [kg]
    static constexpr PS::F64 norm_vel    = sqrt(norm_energy/norm_mass); // [m/s]
    static constexpr PS::F64 norm_time   = norm_length/norm_vel;        // [s]
    static constexpr PS::F64 norm_temp   = norm_energy/blz;             // [K·mol/kcal]
    static constexpr PS::F64 norm_press  = norm_energy/(norm_length
                                                        *norm_length
                                                        *norm_length);  // [mol·m^3·J/kcal]
    static constexpr PS::F64 norm_dens   = norm_mass/(norm_length
                                                      *norm_length
                                                      *norm_length);    // [kg/m^3]
    
    static constexpr PS::F64 coef_coulomb = sqrt(ec*ec/(4.0*pi*e0*norm_length*norm_energy));
    
}

