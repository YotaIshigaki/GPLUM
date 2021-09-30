/* C++ headers */
#include "common.h"
/* User-defined headers */
#include "physical_constants.h"

namespace physical_constants {

    // Physical constants
    // Note that the values below are based on 2018CODATA.
    // (see https://physics.nist.gov/cuu/Constants/index.html for fundamental
    //  physical constants recommended by CODATA)
    const double Ggrav = 6.67430e-8;  // Gravitational constants [cm^{3}/g/s^{2}]
    const double kBoltz = 1.380649e-16; // Boltzmann constant [cm^{2} g/s^{2}/K]
    
    // Mass units
    // Note that the mass of atoms are based on 2018CODATA and CIAAW.
    // (i) dalton is the unified atomic mass unit.
    //     The value is taken from 2018CODATA.
    // (ii) The mass of each atom is taken from CIAAW (https://www.ciaaw.org).
    //      For CIAAW, read the section of Standard atomic weights in
    //      https://en.wikipedia.org/wiki/Relative_atomic_mass
    const double dalton = 1.66053906660e-24; // Unified atomic mass unit [g]
    const double Mhydrogen = 1.00797 * dalton; // 1H + 2H
    const double Mhelium = 4.002602 * dalton; // 3He + 4He
    const double Mcarbon = 12.0107 * dalton; // 12C + 13C
    const double Mnitrogen = 14.00685 * dalton; // 14N + 15N
    const double Moxygen = 15.99940 * dalton; // 16O + 17O + 18O
    const double Mneon = 20.1797 * dalton; // 20Ne + 21Ne + 22Ne
    const double Mmagnesium = 24.3055 * dalton; // 24Mg + 25Mg + 26Mg
    const double Msilicon = 28.085 * dalton; // 28Si + 29Si + 30Si
    const double Msulfur = 32.0675 * dalton; // 32S + 33S + 34S + 36S
    const double Mcalcium = 40.078 * dalton; // 40Ca + 42Ca + 43Ca + 44Ca + 46Ca + 48Ca
    const double Miron = 55.845 * dalton; // 54Fe + 56Fe + 57Fe + 58Fe
    const double Mnickel = 58.6934 * dalton; // 58Ni + 60Ni + 61Ni + 62Ni + 64Ni
    const double Meuropium = 151.964 * dalton; // 151Eu + 153Eu
    const double Msolar = 1.989e33; // Solar mass [g]
    
    // Length units
    const double km = 1.0e5; // kilo-meters [cm]
    const double AU = 1.49597870e13; // Astronomical unit [cm]
    const double pc = 3.0857e18; // parsec [cm]
    const double kpc = 1.0e3 * pc; // kilo-parsecs [cm]
    const double Mpc = 1.0e6 * pc; // mega-parsecs [cm]
    const double Gpc = 1.0e9 * pc; // giga-parsecs [cm]
    
    // Time units
    const double yr = 3.15576e7; // year [s]
    const double kyr = 1.0e3 * yr; // kilo-years [s]
    const double Myr = 1.0e6 * yr; // mega-years [s]
    const double Gyr = 1.0e9 * yr; // giga-years [s]

}

