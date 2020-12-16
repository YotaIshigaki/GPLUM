# -*- coding: utf-8 -*-
#=================================
#   Module import
#=================================
try:
    import math
except ImportError:
    # ignore error when os not found.
    print("Module os,sys,math,struct,types not found.")
    quit()

#=======================================
#   Definitiion of physical constants
#=======================================
# Physical constants [cgs]
#   cLight    :=  Speed of light
#   Ggrav     :=  Gravitational constant
#   kBoltz    :=  Boltzman constant
#   hPlanck   :=  Planck constant
#   hbPlanck  :=  h/(2\pi)
#   aBohr     :=  Bohr radius
#   sigmaBohr :=  Bohr cross section (PI*aBohr**2)
#   sigmaSB   :=  Stefan-Boltmann constant
#   alphaFS   :=  Fine-structure constant
#   gam_ideal :=  Ratio of specific heats of ideal gas
cLight    = 3.0e10
Ggrav     = 6.673e-8
kBoltz    = 1.3806503e-16
hPlanck   = 6.626068e-27
hbPlanck  = 1.054571475e-27
aBohr     = 5.291772108e-9
sigmaBohr = math.pi*aBohr*aBohr
sigmaSB   = 5.670373e-5
alphaFS   = 1.0e0/137.04e0
gam_ideal = 5.0e0/3.0e0
# Length units
micron   = 1.0e-4
Angstrom = 1.0e-8
km       = 1.0e5
Rsolar   = 6.95e10
AU       = 1.49597870e13
pc       = 3.0857e18
kpc      = 3.0857e21
Mpc      = 3.0857e24
Gpc      = 3.0857e27     
# Time units
minute   = 6.0e1
hour     = 3.6e3
day      = 8.64e4
month    = 2.628e6
yr       = 3.15576e7
kyr      = 3.15576e10
Myr      = 3.15576e13
Gyr      = 3.15576e16
# Mass units
#   dalton := the same as the unified atomic mass unit (CODATA2010)
Mlunar    = 7.34581e25
Mearth    = 5.974e27
Msaturn   = 5.688e29
Mjupiter  = 1.8986e30
Msolar    = 1.989e33
dalton    = 1.660538921e-24
Melectron = 9.10938188e-28
Mproton   = 1.67262158e-24
Mneutron  = 1.67492716e-24
MHydro    = 1.673532497e-24
MHelium   = 4.002602e0*dalton
MHydroMol = 2.01588e0*dalton
# Energy units
eV = 1.60217646e-12
# Charge units
e_cgs = 4.8032e-10
e_mks = 1.6022e-19
# Luminosity units
Lbol_solar = 3.85e33

# Others
#   re_Cl         := Classical electron radius
#   sigmaThomson  := Thomson cross section
#   lambdaCompton := Compton wavelength
re_Cl = (e_cgs*e_cgs) \
      / (Melectron*cLight*cLight)
sigmaThomson  = 8.0e0*math.pi*re_Cl*re_Cl/3.0e0
lambdaCompton = hPlanck/(Melectron*cLight)

# Mass abundance
# [Ref.] recommended value in Asplund et al.(2009)[ARA+A,47,481].
XHydro_solar  = 0.7381e0
YHelium_solar = 0.2485e0
Zmetal_solar  = 0.0134e0
# Unit transform
Hz_to_eV     = hPlanck / eV
eV_to_Hz     = eV      / hPlanck
Kelvin_to_eV = kBoltz  / eV
eV_to_Kelvin = eV      / kBoltz

#=======================================
#   Definition of functions
#=======================================
def eV_to_cm(E):
    return cLight*hPlanck/(E*eV)

def cm_to_eV(lmd):
    return (hPlanck*cLight/lmd)/eV

def erg_to_cm(E): 
    return cLight*hPlanck/E

def Hz_to_cm(nu):
    return cLight/nu

def cm_to_Hz(lmd):
    return cLight/lmd
