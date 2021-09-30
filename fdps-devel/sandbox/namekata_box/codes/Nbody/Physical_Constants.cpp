/* Standard headers */
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
/* User-defined headers */
#include "Physical_Constants.h"

namespace physical_constants {

// Physical constants [cgs]
//   cLight    :=  Speed of light
//   Ggrav     :=  Gravitational constant
//   kBoltz    :=  Boltzman constant [erg/K]
//   NAvogadro :=  Avogadro constant [1/mol]
//   Rgasconst :=  gas constant [erg/K/mol]
//   hPlanck   :=  Planck constant
//   hbPlanck  :=  h/(2\pi)
//   aBohr     :=  Bohr radius
//   sigmaBohr :=  Bohr cross section (PI*aBohr**2)
//   sigmaSB   :=  Stefan-Boltmann constant
//   aRDC      :=  radiation (density) constant; =4*sigmaSB/cLight.
//   alphaFS   :=  Fine-structure constant
//   gam_ideal :=  Ratio of specific heats of ideal gas
const double cLight=3.0e10;
const double Ggrav=6.673e-8;
const double kBoltz=1.3806503e-16;
const double NAvogadro=6.02214129e23;
const double R_gas_const=kBoltz*NAvogadro;
const double hPlanck=6.626068e-27,hbPlanck=1.054571475e-27;
const double aBohr=5.291772108e-9;
const double sigmaBohr=M_PI*aBohr*aBohr;
const double sigmaSB=5.670373e-5;
const double aRDC=4.0e0*sigmaSB/cLight;
const double alphaFS=1.0e0/137.04e0;
const double gam_ideal=5.0e0/3.0e0;

// Length units 
const double micron=1.0e-4,Angstrom=1.0e-8;
const double km=1.0e5,Rsolar=6.95e10,AU=1.49597870e13,pc=3.0857e18;
const double kpc=1.0e3*pc,Mpc=1.0e6*pc,Gpc=1.0e9*pc;
// Time units
const double minute=6.0e1,hour=3.6e3,day=8.64e4,month=2.628e6;
const double yr=3.15576e7,kyr=3.15576e10,Myr=3.15576e13,Gyr=3.15576e16;
// Mass units
const double Mlunar=7.34581e25,Mearth=5.974e27,
             Msaturn=5.688e29,Mjupiter=1.8986e30,
             Msolar=1.989e33;

// Mass of elements
//  (1) For each element, name, atomic number, mass number are added as a comment.
//  (2) Elements with high atomic number have many isotopes with different half-lives.
//      Here, the mass of most stable isotope is set to the parameter M_*. 
//      Also, the mass numbers for them are often uncertein.
//  (3) dalton is the same as the unified atomic mass unit (CODATA2010)
const double dalton   = 1.660538921e-24;
const double M_elctrn = 9.10938188e-28;        // electron
const double M_prtn   = 1.67262158e-24;        // proton
const double M_ntrn   = 1.67492716e-24;        // neutron
const double M_H      = 1.008e0*dalton;        // Hydrogen (1,1)
const double M_Hm     = M_H + M_elctrn;        // H-
const double M_H2I    = 2.01588e0*dalton;      // Molecular hydrogen
const double M_H2p    = M_H2I - M_elctrn;      // H2+
const double M_D      = 2.0141017778e0*dalton; // Deuterium (1,2)
const double M_He     = 4.002602e0*dalton;     // Helium (2,4)
const double M_HeII   = M_He - M_elctrn;       // He+
const double M_HeIII  = M_He - 2.0e0*M_elctrn; // He++
const double M_Li     = 6.94e0*dalton;         // Lithium (3,7)
const double M_Be     = 9.012182e0*dalton;     // Beryllium (4,9)
const double M_B      = 10.81e0*dalton;        // Bron (5,11)
const double M_C      = 12.011e0*dalton;       // Carbon (6,12)
const double M_C13    = 13.0033548378e0*dalton;// Carbon-13 (6,13)
const double M_N      = 14.0067e0*dalton;      // Nitrogen (7,14)
const double M_N15    = 15.0001088982e0*dalton;// Nitrogen-15 (7,15)
const double M_O      = 15.999e0*dalton;       // Oxygen (8,16)
const double M_F      = 18.9984032e0*dalton;   // Fluorine (9,19)
const double M_Ne     = 20.1797e0*dalton;      // Neon (10,20)
const double M_Na     = 22.98976928e0*dalton;  // Sodium (11,23)
const double M_Mg     = 24.3050e0*dalton;      // Magnesium (12,24)
const double M_Al     = 26.9815386e0*dalton;   // Aluminium (Aluminum) (13,27)
const double M_Si     = 28.085e0*dalton;       // Silicon (14,28)
const double M_P      = 30.973762e0*dalton;    // Phosphorus (15,31)
const double M_S      = 32.06e0*dalton;        // Sulfur (16,32)
const double M_Cl     = 35.45e0*dalton;        // Chlorine (17,35)
const double M_Ar     = 39.948e0*dalton;       // Argon (18,40)
const double M_K      = 39.0983e0*dalton;      // Potassium (19,39)
const double M_Ca     = 40.078e0*dalton;       // Calcium (20,40) 
const double M_Sc     = 44.955912e0*dalton;    // Scandium (21,45)
const double M_Ti     = 47.867e0*dalton;       // Titanium (22,48)
const double M_V      = 50.9415e0*dalton;      // Vanadium (23,51)
const double M_Cr     = 51.9961e0*dalton;      // Chromium (24,52)
const double M_Mn     = 54.938045e0*dalton;    // Manganese (25,55)
const double M_Fe     = 55.845e0*dalton;       // Iron (26,56)
const double M_Co     = 58.933195e0*dalton;    // Cobalt (27,58)
const double M_Ni     = 58.6934e0*dalton;      // Nickel (28,58)
const double M_Cu     = 63.546e0*dalton;       // Copper (29,64)
const double M_Zn     = 65.38e0*dalton;        // Zinc (30,65)
const double M_Ga     = 69.723e0*dalton;       // Gallium (31,70)
const double M_Ge     = 72.63e0*dalton;        // Germanium (32,73)
const double M_As     = 74.92160e0*dalton;     // Arsenic (33,75)
const double M_Se     = 78.96e0*dalton;        // Selenium (34,79)
const double M_Br     = 79.904e0*dalton;       // Bromine (35,80)
const double M_Kr     = 83.798e0*dalton;       // Krypton (36,84)
const double M_Rb     = 85.4678e0*dalton;      // Rubidium (37,85)
const double M_Sr     = 87.62e0*dalton;        // Strontium (38,88)
const double M_Y      = 88.90585e0*dalton;     // Yttrium (39,89)
const double M_Zr     = 91.224e0*dalton;       // Zirconium (40,91)
const double M_Nb     = 92.90638e0*dalton;     // Niobium (41,93)
const double M_Mo     = 95.96e0*dalton;        // Molybdenum (42,96)
const double M_Tc     = 99.0e0*dalton;         // Technetium (43,98)
const double M_Ru     = 101.07e0*dalton;       // Ruthenium (44,101)
const double M_Rh     = 102.90550e0*dalton;    // Rhodium (45,103)
const double M_Pd     = 106.42e0*dalton;       // Palladium (46,106)
const double M_Ag     = 107.8682e0*dalton;     // Silver (47,108)
const double M_Cd     = 112.411e0*dalton;      // Cadmium (48,112)
const double M_In     = 114.818e0*dalton;      // Indium (49,115)
const double M_Sn     = 118.710e0*dalton;      // Tin (50,119)
const double M_Sb     = 121.760e0*dalton;      // Antimony (51,122)
const double M_Te     = 127.60e0*dalton;       // Tellurium (52,128)
const double M_I      = 126.90447e0*dalton;    // Iodine (53,127)
const double M_Xe     = 131.293e0*dalton;      // Xenon (54,131)
const double M_Cs     = 132.9054519e0*dalton;  // Caesium (Cesium) (55,133)
const double M_Ba     = 137.327e0*dalton;      // Barium (56,137)
const double M_La     = 138.90547e0*dalton;    // Lanthanum (57,139)
const double M_Ce     = 140.116e0*dalton;      // Cerium (58,140)
const double M_Pr     = 140.90765e0*dalton;    // Praseodymium (59,141)
const double M_Nd     = 144.242e0*dalton;      // Neodymium (60,144)
const double M_Pm     = 145.0e0*dalton;        // Promethium (61,145)
const double M_Sm     = 150.36e0*dalton;       // Samarium (62,150)
const double M_Eu     = 151.964e0*dalton;      // Europium (63,152)
const double M_Gd     = 157.25e0*dalton;       // Gadolinium (64,157)
const double M_Tb     = 158.92535e0*dalton;    // Terbium (65,159)
const double M_Dy     = 162.500e0*dalton;      // Dysprosium (66,163)
const double M_Ho     = 164.93032e0*dalton;    // Holmium (67,165)
const double M_Er     = 167.259e0*dalton;      // Erbium (68,167)
const double M_Tm     = 168.93421e0*dalton;    // Thulium (69,169)
const double M_Yb     = 173.054e0*dalton;      // Ytterbium (70,173)
const double M_Lu     = 174.9668e0*dalton;     // Lutetium (71,175)
const double M_Hf     = 178.49e0*dalton;       // Hafnium (72,177)
const double M_Ta     = 180.94788e0*dalton;    // Tantalum (73,181)
const double M_W      = 183.84e0*dalton;       // Tungsten (74,184)
const double M_Re     = 186.207e0*dalton;      // Rhenium (75,186)
const double M_Os     = 190.23e0*dalton;       // Osmium (76,190)
const double M_Ir     = 192.217e0*dalton;      // Iridium (77,192)
const double M_Pt     = 195.084e0*dalton;      // Platinum (78,195)
const double M_Au     = 196.966569e0*dalton;   // Gold (79,197)
const double M_Hg     = 200.59e0*dalton;       // Mercury (80,201)
const double M_Tl     = 204.38e0*dalton;       // Thallium (81,204)
const double M_Pb     = 207.2e0*dalton;        // Lead (82,207)
const double M_Bi     = 208.98040e0*dalton;    // Bismuth (83,209)
const double M_Po     = 209.0e0*dalton;        // Polonium (84,209)
const double M_At     = 210.0e0*dalton;        // Astatine (85,210)
const double M_Rn     = 222.0e0*dalton;        // Radon (86,222)
const double M_Fr     = 223.0e0*dalton;        // Francium (87,223)
const double M_Ra     = 226.0e0*dalton;        // Radium (88,226)
const double M_Ac     = 227.0e0*dalton;        // Actinium (89,227)
const double M_Th     = 232.03806e0*dalton;    // Thorium (90,232)
const double M_Pa     = 231.03588e0*dalton;    // Protactinium (91,231)
const double M_U      = 238.02891e0*dalton;    // Uranium (92,238)
const double M_Np     = 237.0e0*dalton;        // Neptunium (93,237)
const double M_Pu     = 244.0e0*dalton;        // Plutonium (94,244)
const double M_Am     = 243.0e0*dalton;        // Americium (95,243)
const double M_Cm     = 247.0e0*dalton;        // Curium (96,247)
const double M_Bk     = 247.0e0*dalton;        // Berkelium (97,247)
const double M_Cf     = 251.0e0*dalton;        // Californium (98,251)
const double M_Es     = 252.0e0*dalton;        // Einsteinium (99,252)
const double M_Fm     = 257.0e0*dalton;        // Fermium (100,257)
const double M_Md     = 258.0e0*dalton;        // Mendelevium (101,258) 
const double M_No     = 259.0e0*dalton;        // Nobelium (102,259)
const double M_Lr     = 266.0e0*dalton;        // Lawrencium (103,266)
const double M_Rf     = 267.0e0*dalton;        // Rutherfordium (104,267)
const double M_Db     = 268.0e0*dalton;        // Dubnium (105,268)
const double M_Sg     = 271.0e0*dalton;        // Seaborgium (106,271)
const double M_Bh     = 270.0e0*dalton;        // Bohrium (107,270)
const double M_Hs     = 269.0e0*dalton;        // Hassium (108,269) 
const double M_Mt     = 278.0e0*dalton;        // Meitnerium (109,278)
const double M_Ds     = 281.0e0*dalton;        // Darmstadtium (110,281)
const double M_Rg     = 281.0e0*dalton;        // Roentgenium (111,281)
const double M_Cn     = 285.0e0*dalton;        // Copernicium (112,285)
const double M_Uut    = 286.0e0*dalton;        // Ununtrium (113,286)
const double M_Fl     = 289.0e0*dalton;        // Flerovium (114,289)
const double M_Uup    = 289.0e0*dalton;        // Ununpentium (115,289)
const double M_Lv     = 293.0e0*dalton;        // Livermorium (116,293)
const double M_Uus    = 294.0e0*dalton;        // Ununseptium (117,294)
const double M_Uuo    = 294.0e0*dalton;        // Ununoctium (118,294)

// Energy units
const double eV=1.60217646e-12;
const double keV=1.0e3*eV,MeV=1.0e6*eV,GeV=1.0e9*eV,
             TeV=1.0e12*eV,PeV=1.0e15*eV;
// Charge units
const double e_cgs=4.8032e-10,e_mks=1.6022e-19;
// Luminosity units
const double Lbol_solar=3.85e33;

// Others
//   re_Cl         := Classical electron radius
//   sigmaThomson  := Thomson cross section
//   lambdaCompton := Compton wavelength
const double re_Cl=(e_cgs*e_cgs)/(M_elctrn*cLight*cLight);
const double sigmaThomson=8.0e0*M_PI*re_Cl*re_Cl/3.0e0;
const double lambdaCompton=hPlanck/(M_elctrn*cLight);

// Mass abundance
// [Ref.] recommended value in Asplund et al.(2009)[ARA+A,47,481].
const double XHydro_solar=0.7381e0;
const double YHelium_solar=0.2485e0;
const double Zmetal_solar=0.0134e0;

// Unit transform
const double Hz_to_eV=hPlanck/eV;
const double eV_to_Hz=eV/hPlanck;
const double Kelvin_to_eV=kBoltz/eV;
const double eV_to_Kelvin=eV/kBoltz;

/*-------------------------------------------------------------------*/
/////////////////////////   F U N C T I O N   /////////////////////////
///////////////////////// < E V _ T O _ C M > /////////////////////////
/*-------------------------------------------------------------------*/
inline double eV_to_cm(const double E){ 
   return cLight*hPlanck/(E*eV);
}

/*-------------------------------------------------------------------*/
////////////////////////    F U N C T I O N    ////////////////////////
//////////////////////// < E R G _ T O _ C M > ////////////////////////
/*-------------------------------------------------------------------*/
inline double erg_to_cm(const double E){
   return cLight*hPlanck/E;
}

/*-------------------------------------------------------------------*/
/////////////////////////   F U N C T I O N   /////////////////////////
///////////////////////// < H Z _ T O _ C M > /////////////////////////
/*-------------------------------------------------------------------*/
inline double Hz_to_cm(const double nu){
   return cLight/nu;
}

}
