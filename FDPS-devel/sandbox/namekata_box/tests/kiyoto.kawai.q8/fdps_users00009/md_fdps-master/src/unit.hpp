/**************************************************************************************************/
/**
* @file  unit.hpp
* @brief definition of physical constants and normalized parameters.
*/
/**************************************************************************************************/
#pragma once

#include <cmath>
#include <iostream>

/**
* @brief definitions of physical constants and normalized parameters.
*/
namespace Unit {

    //--- fundamental parameters
    static constexpr PS::F64 pi    = 3.141592653589793;  // circular constant
    static constexpr PS::F64 angm  = 1.0e-10;            // [m] [angstrom] to [meter]

    //--- ref: The NIST reference (http://physics.nist.gov/cuu/Constants/index.html)
    static constexpr PS::F64 ec    = 1.6021766208e-19;   // [C] charge of electron
    static constexpr PS::F64 e0    = 8.854187817620e-12; // [C^2/J·m] permittivity of vacuum
    static constexpr PS::F64 avgdr = 6.022140857e+23;    // [/mol] Avogadro constant
    static constexpr PS::F64 blz   = 1.38064852e-23;     // [J/K] Boltzmann constant

    //--- ref: The elements, Theodore gray, 2010
    static constexpr PS::F64 mass_C        = 12.0107e-3;  // [kg/mol] mass of Carbon atom
    //--- ref: JSME text series "Thermodynamics"
    static constexpr PS::F64 mech_heat_equ = 4.1868e+3;   // [J/kcal] mechanical equivalent of heat

    //--- standard values                                                    // J = [kg·m^2 / s^2]
    static constexpr PS::F64 norm_energy = mech_heat_equ/avgdr;              // [J·mol/kcal]
    static constexpr PS::F64 norm_length = angm;                             // [m]
    static constexpr PS::F64 norm_mass   = mass_C/avgdr;                     // [kg]
    static constexpr PS::F64 norm_vel    = std::sqrt(norm_energy/norm_mass); // [m/s]
    static constexpr PS::F64 norm_time   = norm_length/norm_vel;             // [s]
    static constexpr PS::F64 norm_temp   = norm_energy/blz;                  // [K·mol/kcal]
    static constexpr PS::F64 norm_press  = norm_energy/( std::pow(norm_length, 3) );  // [J·mol/kcal·m^3]
    static constexpr PS::F64 norm_dens   = norm_mass/(   std::pow(norm_length, 3) );  // [kg/m^3]

    static constexpr PS::F64 femto_second = 1.e-15;                                                // [s]
    static constexpr PS::F64 coef_coulomb = std::sqrt(ec*ec/(4.0*pi*e0*norm_length*norm_energy));  // [sqrt(kcal/mol)]


    //! @brief displaying normalized units.
    void print_unit(){
		std::cout << " Normarized units:\n"
		          << "   energy       = " << std::setw(18) << norm_energy  << " [J·mol/kcal]\n"
				  << "   length       = " << std::setw(18) << norm_length  << " [m]\n"
				  << "   mass         = " << std::setw(18) << norm_mass    << " [kg]\n"
				  << "   velocity     = " << std::setw(18) << norm_vel     << " [m/s]\n"
				  << "   time         = " << std::setw(18) << norm_time    << " [s]\n"
				  << "   temperature  = " << std::setw(18) << norm_temp    << " [K·mol/kcal]\n"
				  << "   pressure     = " << std::setw(18) << norm_press   << " [J·mol/kcal·m^3]\n"
				  << "   density      = " << std::setw(18) << norm_dens    << " [kg/m^3]\n"
				  << "   coulomb_coef = " << std::setw(18) << coef_coulomb << " [sqrt(kcal/mol)]\n" << std::endl;
	}

    //! @brief convert normalized time into real time
    PS::F64 to_real_time(const PS::F64 norm_time) { return norm_time/Unit::femto_second*Unit::norm_time; }

    //! @brief convert real time into normalized time
    PS::F64 to_norm_time(const PS::F64 real_time) { return real_time*Unit::femto_second/Unit::norm_time; }
}
