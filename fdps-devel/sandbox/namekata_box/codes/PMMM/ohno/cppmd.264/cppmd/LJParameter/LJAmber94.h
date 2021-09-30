#ifndef LJAMBER94_H
#define LJAMBER94_H
#include "LJAmber.h"

/* LJ parameters from  AMBER 94 */
/*
  http://ambermd.org/
  http://ambermd.org/dbase.html
  http://ambermd.org/amber10.ffparms.tar.bz2
  parm/parm94.dat
 */

//! generic "extra point"
static LJAmberParameter LJAmberParameterEP = {"EP",0.0,0.0}; 

//! Lennard Jones Parameters in parm/parm94.dat of AMBER10
static LJAmberParameter ljamber94parameters[] = {
    {"H", 0.6000,0.0157},         // !Ferguson base pair geom.
    {"HO",0.0000,0.0000},         // OPLS Jorgensen, JACS,110,(1988),1657
    {"HS",0.6000,0.0157},         // W. Cornell CH3SH --> CH3OH FEP
    {"HC",1.4870,0.0157},         // OPLS
    {"H1",1.3870,0.0157},         // Veenstra et al JCC,8,(1992),963
    {"H2",1.2870,0.0157},         // Veenstra et al JCC,8,(1992),963
    {"H3",1.1870,0.0157},         // Veenstra et al JCC,8,(1992),963
    {"HP",1.1000,0.0157},         // Veenstra et al JCC,8,(1992),963
    {"HA",1.4590,0.0150},         // Spellmeyer
    {"H4",1.4090,0.0150},         // Spellmeyer, one electrowithdr. neighbor
    {"H5",1.3590,0.0150},         // Spellmeyer, two electrowithdr. neighbor
    {"HW",0.0000,0.0000},         // TIP3P water model
    {"O", 1.6612,0.2100},         // OPLS
    {"O2",1.6612,0.2100},         // OPLS
#ifdef DEBUG_COULOMB
    {"OW",0.0000,0.0000},         // TIP3P water model
#else
    {"OW",1.7683,0.1520},         // TIP3P water model
#endif
    {"OH",1.7210,0.2104},         // OPLS
    {"OS",1.6837,0.1700},         // OPLS ether
    {"CT",1.9080,0.1094},         // Spellmeyer
    {"CA",1.9080,0.0860},         // Spellmeyer
    {"CM",1.9080,0.0860},         // Spellmeyer
    {"C", 1.9080,0.0860},         // OPLS
    {"N", 1.8240,0.1700},         // OPLS
    {"N3",1.8240,0.1700},         // OPLS
    {"S", 2.0000,0.2500},         // W. Cornell CH3SH and CH3SCH3 FEP's
    {"SH",2.0000,0.2500},         // W. Cornell CH3SH and CH3SCH3 FEP's
    {"P", 2.1000,0.2000},         // JCC,7,(1986),230;
    {"IM",2.47,0.1},              // Cl- Smith & Dang, JCP 1994,100:5,3757
    {"Li",1.1370,0.0183},         // Li+ Aqvist JPC 1990,94,8021. (adapted)
    {"IP",1.8680,0.00277},        // Na+ Aqvist JPC 1990,94,8021. (adapted)
    {"K", 2.6580,0.000328},       // K+ Aqvist JPC 1990,94,8021. (adapted)
    {"Rb",2.9560,0.00017},        // Rb+ Aqvist JPC 1990,94,8021. (adapted)
    {"Cs",3.3950,0.0000806},      // Cs+ Aqvist JPC 1990,94,8021. (adapted)
    {"I", 2.35,0.40},             // JCC,7,(1986),230;
    {"F", 1.75,0.061},            // Gough et al. JCC 13,(1992),963.
    {"IB",5.0,0.1},               // solvated ion for vacuum approximation
    LJAmberParameterEP            // generic "extra point"
};

//! parm94 of AMBER10
class LJAmber94 : public LJAmber {
 public:
 LJAmber94() : LJAmber(ljamber94parameters,sizeof(ljamber94parameters)/sizeof(LJAmberParameter)) {}
};

#endif
