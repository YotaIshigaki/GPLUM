#pragma once 

#define REAL    double

#define NEWTON_SI_CGS   1.e+05
#define LENGTH_SI_CGS   1.e+02
#define MASS_SI_CGS     1.e+03
#define TIME_SI_CGS     1.e+00
#define ENERGY_SI_CGS   1.e+07
#define CHARGE_SI_CGS   3.e+09

// fundermental constants
#define GRAVITY_CONSTANT 6.6725985e-11   // [Nm^2/kg^2] = [m^3/kg/s^2]
#define LIGHT_C          2.997924580e+08 // [m/s]
#define PLANK_HBAR       1.0545726663e-34 // [Js]

// fundermental constants [cgs]
#define GRAVITY_CONSTANT_CGS 6.6725985e-8     // [dyne m^2/kg^2] = [cm^3/g/s^2]
#define LIGHT_C_CGS          2.997924580e+10  // [cm/s]
#define PLANK_HBAR_CGS       1.0545726663e-27 // [ergs]

// constants about nucleus
#define ELEMENTAL_CHARGE 1.6021773349e-19 // [C]
#define ELECTRON_MASS    9.109389754e-31 // [kg]
#define ELECTRON_CHARGE  (-ELEMENTAL_CHARGE) // [C]
#define PROTON_MASS      1.672623110e-27 // [kg]
#define PROTON_CHARGE    (+ELEMENTAL_CHARGE) // [C]
#define NEUTRON_MASS     1.674928610e-27 // [kg]

// constants about nucleus [cgs]
#define ELEMENTAL_CHARGE_CGS 4.8065320047e-10    // [esu]
#define ELECTRON_MASS_CGS    9.109389754e-28     // [g]
#define ELECTRON_CHARGE_CGS  (-ELEMENTAL_CHARGE) // [esu]
#define PROTON_MASS_CGS      1.672623110e-24     // [g]
#define PROTON_CHARGE_CGS    (+ELEMENTAL_CHARGE) // [esu]
#define NEUTRON_MASS_CGS     1.674928610e-24     // [g]

// constants about chemical physics
#define AVOGADROS_NUMBER 6.022136736e+23 // [1/mol]
#define BOLTZMANN_CONSTANT 1.38065812e-23 // [J/K]
#define UNIVERSAL_GAS_CONSTANT (BOLTZMANN_CONSTANT/AVOGADROS_NUMBER) // [J/K/mol]

// constants about chemical physics [cgs]
#define AVOGADROS_NUMBER_CGS   6.022136736e+23 // [1/mol]
#define BOLTZMANN_CONSTANT_CGS 1.38065812e-16  // [erg/K]
#define UNIVERSAL_GAS_CONSTANT_CGS (BOLTZMANN_CONSTANT_CGS/AVOGADROS_NUMBER_CGS) // [erg/K/mol]

// constants about atom
#define ALPHA            (ELECTRON_CHARGE*ELECTRON_CHARGE\
                         *LIGHT_C/PLANK_HBAR/1e7)   // [1]
#define ALPHA_C          (ALPHA*LIGHT_C)            // [m/s]
#define ALPHA_C2         (ALPHA_C*ALPHA_C)          // [m^2/s^2]
#define HBAR_ALPHA_C     (PLANK_HBAR*ALPHA*LIGHT_C) // [Jm]
#define HBAR_C           (PLANK_HBAR*LIGHT_C)       // [Jm]

#define BOHR_RADIUS      (PLANK_HBAR/(ELECTRON_MASS*ALPHA_C))  // [m]
#define RYDBERG_CONSTANT (1.097373153413e7) // [1/m]
#define CLASSICAL_ELECTRON_RADIUS (HBAR_ALPHA_C/(ELECTRON_MASS*LIGHT_C*LIGHT_C)) // [m]

// constants about atom [cgs]
#define ALPHA_CGS            (ELECTRON_CHARGE_CGS*ELECTRON_CHARGE_CGS\
                              *LIGHT_C_CGS/PLANK_HBAR_CGS/1e7)   // [1]
#define ALPHA_C_CGS          (ALPHA_CGS*LIGHT_C_CGS)            // [cm/s]
#define ALPHA_C2_CGS         (ALPHA_C_CGS*ALPHA_C_CGS)          // [cm^2/s^2]
#define HBAR_ALPHA_C_CGS     (PLANK_HBAR_CGS*ALPHA_CGS*LIGHT_C_CGS) // [erg cm]
#define HBAR_C_CGS           (PLANK_HBAR_CGS*LIGHT_C_CGS)       // [erg cm]

#define BOHR_RADIUS_CGS      (PLANK_HBAR_CGS/(ELECTRON_MASS_CGS*ALPHA_C_CGS))  // [cm]
#define RYDBERG_CONSTANT_CGS (1.097373153413e5) // [1/cm]
#define CLASSICAL_ELECTRON_RADIUS_CGS (HBAR_ALPHA_C_CGS/(ELECTRON_MASS_CGS*LIGHT_C_CGS*LIGHT_C_CGS)) // [cm]


// base quantities of atomic unit
#define ENERGY_AU        (ELECTRON_MASS*ALPHA_C2)          // [J]
#define MOMENT_AU        (ELECTRON_MASS*ALPHA_C)           // [kgm/s]
#define LENGTH_AU        (PLANK_HBAR/(ELECTRON_MASS*ALPHA_C))  // [m]
#define TIME_AU          (PLANK_HBAR/(ELECTRON_MASS*ALPHA_C2)) // [s]

// base quantities of atomic unit [cgs]
#define ENERGY_AU_CGS        (ELECTRON_MASS_CGS*ALPHA_C2_CGS)          // [erg]
#define MOMENT_AU_CGS        (ELECTRON_MASS_CGS*ALPHA_C_CGS)           // [gcm/s]
#define LENGTH_AU_CGS        (PLANK_HBAR_CGS/(ELECTRON_MASS_CGS*ALPHA_C_CGS))  // [cm]
#define TIME_AU_CGS          (PLANK_HBAR_CGS/(ELECTRON_MASS_CGS*ALPHA_C2_CGS)) // [s]


// energy convertor
#define J_TO_EV          (1.0/1.602177335e-19) // [eV/J]
#define EV_TO_J          (1.602177335e-19)     // [J/eV]

// energy convertor [cgs]
#define ERG_TO_EV          (1.0/1.602177335e-12) // [eV/erg]
#define EV_TO_ERG          (1.602177335e-12)     // [erg/eV]


// astrophysical constants
#define MSUN_CGS    (1.98892e+33)       // [g]
#define MEARTH_CGS  (5.97258e+27)       // [g]
#define AU_CGS      (1.49597870e+13)    // [cm]
#define PC_CGS      (3.08568025e+18)    // [cm]
#define KPC_CGS     (PC_CGS*1.e+3)      // [cm]
#define MPC_CGS     (PC_CGS*1.e+6)      // [cm]
#define YEAR_CGS    (3.1556926e+7)      // [s]
#define MEGAYEAR_CGS   (YEAR_CGS*1.e+6) // [s]
#define GIGAYEAR_CGS   (YEAR_CGS*1.e+9) // [s]
#define HUBBLE_CGS     (100.0*MPC_CGS*1.e-5) // [Mpc/km/s] -> [1/s]
#define CMB_TEMPERATURE (2.735)         // [K] /* The CMD Temperature observed by COBE 2.735+-0.006 K */

// some convenient units
#define VELOCITY_KMS_CGS    (1.e+5) // [km/s] -> [cm/s]

// some convenient constants for mathematical operation 
#ifndef M_PI //pi 
#define M_PI 3.14159265358979323846
#endif

#ifndef M_E //e
#define M_E 2.7182818284590452354
#endif

#ifndef M_LOG2E //log_2 e
#define M_LOG2E 1.4426950408889634074
#endif

#ifndef M_LOG10E //log_10 e
#define M_LOG10E 0.43429448190325182765
#endif

#ifndef M_LN2 //log_e 2
#define M_LN2 0.69314718055994530942
#endif

#ifndef M_LN10 //log_e 10
#define M_LN10 2.30258509299404568402
#endif

#ifndef M_PI_2 //pi/2
#define M_PI_2 1.57079632679489661923
#endif

#ifndef M_PI_4 //pi/4
#define M_PI_4 0.78539816339744830962
#endif

#ifndef M_1_PI //1/pi
#define M_1_PI 0.31830988618379067154
#endif

#ifndef M_2_PI //2/pi
#define M_2_PI 0.63661977236758134308
#endif

#ifndef M_SQRTPI //sqrt(pi)
#define M_SQRTPI 1.77245385090551602729
#endif

#ifndef M_2_SQRTPI //2/sqrt(pi)
#define M_2_SQRTPI 1.12837916709551257390
#endif

#ifndef M_SQRT2 //sqrt(2)
#define M_SQRT2 1.41421356237309504880
#endif

#ifndef M_SQRT3 //sqrt(3)
#define M_SQRT3 1.73205080756887729352
#endif

#ifndef M_SQRT1_2 //1/sqrt(2)
#define M_SQRT1_2 0.70710678118654752440
#endif

#ifndef M_LNPI //log_e(pi)
#define M_LNPI 1.14472988584940017414
#endif

#ifndef M_EULER // eular constant
#define M_EULER 0.57721566490153286061
#endif

