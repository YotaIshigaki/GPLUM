#ifndef AMBERFILE_AMBERPRMTOP_H
#define AMBERFILE_AMBERPRMTOP_H

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif  // HAVE_CONFIG_H

#include <cstdio>
#include <string>
#include <vector>
#include "cppmd-int.h"

namespace amberfile {

template <typename T>
struct AmberPrmtopT {
  typedef T                 int_type;

  AmberPrmtopT()
      : version(0000.000), date("00/00/00"), time("00:00:00"),
        natom(0), ntypes(0), nbonh(0), mbona(0), ntheth(0), mtheta(0),
        nphih(0), mphia(0), nhparm(0), nparm(0), nnb(0), nres(0), nbona(0),
        ntheta(0), nphia(0), numbnd(0), numang(0), nptra(0), natyp(0), nphb(0),
        ifpert(0), nbper(0), ngper(0), ndper(0), mbper(0), mgper(0), mdper(0),
        ifbox(0), nmxrs(0), ifcap(0), numextra(0), ncopy(0) {
  }

  float                     version;
  std::string               date;
  std::string               time;

  std::string               title;

  T                         natom, ntypes, nbonh, mbona, ntheth, mtheta, nphih,
                            mphia, nhparm, nparm, nnb, nres, nbona, ntheta,
                            nphia, numbnd, numang, nptra, natyp, nphb, ifpert,
                            nbper, ngper, ndper, mbper, mgper, mdper, ifbox,
                            nmxrs, ifcap, numextra, ncopy;

  std::vector<std::string>  graph;
  std::vector<double>       chrg;
  std::vector<double>       amass;
  std::vector<T>            iac;
  std::vector<T>            numex;  //!< Number of Excluded Atom
  std::vector<T>            ico;
  std::vector<std::string>  labres;
  std::vector<T>            ipres;

  std::vector<double>       rk, req, tk, teq, pk, pn, phase, solty, cn1, cn2;

  std::vector<T>            bh, b, th, t, ph, p, natex;  //!< natex: Excluded Atom List

  std::vector<double>       asol, bsol, hbcut;

  std::vector<std::string>  isymbl;
  std::vector<std::string>  itree;
  std::vector<T>            join;
  std::vector<T>            irotat;

#if 0
  T                         iptres, nspm, nspsol;
  std::vector<T>            nsp;
  std::vector<double>       beta;
  std::vector<double>       box;
#endif
};

// Define CPPMD_ENABLE_LARGEMODEL in config.h when you use 64-bit integer
// in Amber file I/O.
//
#ifdef CPPMD_ENABLE_LARGEMODEL
typedef AmberPrmtopT<int64_t> AmberPrmtop;
#else  // CPPMD_ENABLE_LARGEMODEL
typedef AmberPrmtopT<int32_t> AmberPrmtop;
#endif  // CPPMD_ENABLE_LARGEMODEL

int ReadFromFile(const std::string& file_name, AmberPrmtop* prm);
int ReadFromFile(FILE* fp, AmberPrmtop* prm);

}  // namespcae amberfile

#endif  // AMBERFILE_AMBERPRMTOP_H
