#pragma once

namespace physical_constants {

// Physical constants
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
extern const double cLight;
extern const double Ggrav;
extern const double kBoltz;
extern const double NAvogadro;
extern const double R_gas_const;
extern const double hPlanck,hbPlanck;
extern const double aBohr;
extern const double sigmaBohr;
extern const double sigmaSB;
extern const double aRDC;
extern const double alphaFS;
extern const double gam_ideal;

// Length units 
extern const double micron,Angstrom;
extern const double km,Rsolar,AU,pc;
extern const double kpc,Mpc,Gpc;
// Time units
extern const double minute,hour,day,month;
extern const double yr,kyr,Myr,Gyr;
// Mass units
extern const double Mlunar,Mearth,Msaturn,Mjupiter,Msolar;

// Mass of elements
//  (1) For each element, name, atomic number, mass number are added as a comment.
//  (2) Elements with high atomic number have many isotopes with different half-lives.
//      Here, the mass of most stable isotope is set to the parameter M_*. 
//      Also, the mass numbers for them are often uncertein.
//  (3) dalton is the same as the unified atomic mass unit (CODATA2010)
extern const double dalton;
extern const double M_elctrn; // electron
extern const double M_prtn;   // proton
extern const double M_ntrn;   // neutron
extern const double M_H;      // Hydrogen (1,1)
extern const double M_Hm;     // H-
extern const double M_H2I;    // Molecular hydrogen
extern const double M_H2p;    // H2+
extern const double M_D;      // Deuterium (1,2)
extern const double M_He;     // Helium (2,4)
extern const double M_HeII;   // He+
extern const double M_HeIII;  // He++
extern const double M_Li;     // Lithium (3,7)
extern const double M_Be;     // Beryllium (4,9)
extern const double M_B;      // Bron (5,11)
extern const double M_C;      // Carbon (6,12)
extern const double M_C13;    // Carbon-13 (6,13)
extern const double M_N;      // Nitrogen (7,14)
extern const double M_N15;    // Nitrogen-15 (7,15)
extern const double M_O;      // Oxygen (8,16)
extern const double M_F;      // Fluorine (9,19)
extern const double M_Ne;     // Neon (10,20)
extern const double M_Na;     // Sodium (11,23)
extern const double M_Mg;     // Magnesium (12,24)
extern const double M_Al;     // Aluminium (Aluminum) (13,27)
extern const double M_Si;     // Silicon (14,28)
extern const double M_P;      // Phosphorus (15,31)
extern const double M_S;      // Sulfur (16,32)
extern const double M_Cl;     // Chlorine (17,35)
extern const double M_Ar;     // Argon (18,40)
extern const double M_K;      // Potassium (19,39)
extern const double M_Ca;     // Calcium (20,40) 
extern const double M_Sc;     // Scandium (21,45)
extern const double M_Ti;     // Titanium (22,48)
extern const double M_V;      // Vanadium (23,51)
extern const double M_Cr;     // Chromium (24,52)
extern const double M_Mn;     // Manganese (25,55)
extern const double M_Fe;     // Iron (26,56)
extern const double M_Co;     // Cobalt (27,58)
extern const double M_Ni;     // Nickel (28,58)
extern const double M_Cu;     // Copper (29,64)
extern const double M_Zn;     // Zinc (30,65)
extern const double M_Ga;     // Gallium (31,70)
extern const double M_Ge;     // Germanium (32,73)
extern const double M_As;     // Arsenic (33,75)
extern const double M_Se;     // Selenium (34,79)
extern const double M_Br;     // Bromine (35,80)
extern const double M_Kr;     // Krypton (36,84)
extern const double M_Rb;     // Rubidium (37,85)
extern const double M_Sr;     // Strontium (38,88)
extern const double M_Y;      // Yttrium (39,89)
extern const double M_Zr;     // Zirconium (40,91)
extern const double M_Nb;     // Niobium (41,93)
extern const double M_Mo;     // Molybdenum (42,96)
extern const double M_Tc;     // Technetium (43,98)
extern const double M_Ru;     // Ruthenium (44,101)
extern const double M_Rh;     // Rhodium (45,103)
extern const double M_Pd;     // Palladium (46,106)
extern const double M_Ag;     // Silver (47,108)
extern const double M_Cd;     // Cadmium (48,112)
extern const double M_In;     // Indium (49,115)
extern const double M_Sn;     // Tin (50,119)
extern const double M_Sb;     // Antimony (51,122)
extern const double M_Te;     // Tellurium (52,128)
extern const double M_I;      // Iodine (53,127)
extern const double M_Xe;     // Xenon (54,131)
extern const double M_Cs;     // Caesium (Cesium) (55,133)
extern const double M_Ba;     // Barium (56,137)
extern const double M_La;     // Lanthanum (57,139)
extern const double M_Ce;     // Cerium (58,140)
extern const double M_Pr;     // Praseodymium (59,141)
extern const double M_Nd;     // Neodymium (60,144)
extern const double M_Pm;     // Promethium (61,145)
extern const double M_Sm;     // Samarium (62,150)
extern const double M_Eu;     // Europium (63,152)
extern const double M_Gd;     // Gadolinium (64,157)
extern const double M_Tb;     // Terbium (65,159)
extern const double M_Dy;     // Dysprosium (66,163)
extern const double M_Ho;     // Holmium (67,165)
extern const double M_Er;     // Erbium (68,167)
extern const double M_Tm;     // Thulium (69,169)
extern const double M_Yb;     // Ytterbium (70,173)
extern const double M_Lu;     // Lutetium (71,175)
extern const double M_Hf;     // Hafnium (72,177)
extern const double M_Ta;     // Tantalum (73,181)
extern const double M_W;      // Tungsten (74,184)
extern const double M_Re;     // Rhenium (75,186)
extern const double M_Os;     // Osmium (76,190)
extern const double M_Ir;     // Iridium (77,192)
extern const double M_Pt;     // Platinum (78,195)
extern const double M_Au;     // Gold (79,197)
extern const double M_Hg;     // Mercury (80,201)
extern const double M_Tl;     // Thallium (81,204)
extern const double M_Pb;     // Lead (82,207)
extern const double M_Bi;     // Bismuth (83,209)
extern const double M_Po;     // Polonium (84,209)
extern const double M_At;     // Astatine (85,210)
extern const double M_Rn;     // Radon (86,222)
extern const double M_Fr;     // Francium (87,223)
extern const double M_Ra;     // Radium (88,226)
extern const double M_Ac;     // Actinium (89,227)
extern const double M_Th;     // Thorium (90,232)
extern const double M_Pa;     // Protactinium (91,231)
extern const double M_U;      // Uranium (92,238)
extern const double M_Np;     // Neptunium (93,237)
extern const double M_Pu;     // Plutonium (94,244)
extern const double M_Am;     // Americium (95,243)
extern const double M_Cm;     // Curium (96,247)
extern const double M_Bk;     // Berkelium (97,247)
extern const double M_Cf;     // Californium (98,251)
extern const double M_Es;     // Einsteinium (99,252)
extern const double M_Fm;     // Fermium (100,257)
extern const double M_Md;     // Mendelevium (101,258) 
extern const double M_No;     // Nobelium (102,259)
extern const double M_Lr;     // Lawrencium (103,266)
extern const double M_Rf;     // Rutherfordium (104,267)
extern const double M_Db;     // Dubnium (105,268)
extern const double M_Sg;     // Seaborgium (106,271)
extern const double M_Bh;     // Bohrium (107,270)
extern const double M_Hs;     // Hassium (108,269) 
extern const double M_Mt;     // Meitnerium (109,278)
extern const double M_Ds;     // Darmstadtium (110,281)
extern const double M_Rg;     // Roentgenium (111,281)
extern const double M_Cn;     // Copernicium (112,285)
extern const double M_Uut;    // Ununtrium (113,286)
extern const double M_Fl;     // Flerovium (114,289)
extern const double M_Uup;    // Ununpentium (115,289)
extern const double M_Lv;     // Livermorium (116,293)
extern const double M_Uus;    // Ununseptium (117,294)
extern const double M_Uuo;    // Ununoctium (118,294)

// Energy units
extern const double eV,keV,MeV,GeV,TeV,PeV;
// Charge units
extern const double e_cgs,e_mks;
// Luminosity units
extern const double Lbol_solar;

// Others
//   re_Cl         := Classical electron radius
//   sigmaThomson  := Thomson cross section
//   lambdaCompton := Compton wavelength
extern const double re_Cl;
extern const double sigmaThomson;
extern const double lambdaCompton;

// Mass abundance
// [Ref.] recommended value in Asplund et al.(2009)[ARA+A,47,481].
extern const double XHydro_solar;
extern const double YHelium_solar;
extern const double Zmetal_solar;

// Unit transform
extern const double Hz_to_eV;
extern const double eV_to_Hz;
extern const double Kelvin_to_eV;
extern const double eV_to_Kelvin;
extern double eV_to_cm(const double E);
extern double erg_to_cm(const double E);
extern double Hz_to_cm(const double nu);

}
