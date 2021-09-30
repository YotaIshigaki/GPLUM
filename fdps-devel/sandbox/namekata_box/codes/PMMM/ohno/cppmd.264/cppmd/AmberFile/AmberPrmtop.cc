#include "AmberPrmtop.h"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <map>
#include <sstream>
#include <string>

#include "ReadAmber-inl.h"

namespace {

void ReadPointers(const amberfile::readamb::VersionInfo& version,
                  const amberfile::readamb::SectionIndex& index,
                  amberfile::AmberPrmtop* prm, FILE* fp) {
  using namespace std;
  vector<amberfile::AmberPrmtop::int_type> intvec;
  ReadData(version, index, "POINTERS", 30, 'I', &intvec, fp);
  vector<amberfile::AmberPrmtop::int_type>::const_iterator it = intvec.begin();
  prm->natom = *it++;
  prm->ntypes = *it++;
  prm->nbonh = *it++;
  prm->mbona = *it++;
  prm->ntheth = *it++;
  prm->mtheta = *it++;
  prm->nphih = *it++;
  prm->mphia = *it++;
  prm->nhparm = *it++;
  prm->nparm = *it++;
  prm->nnb = *it++;
  prm->nres = *it++;
  prm->nbona = *it++;
  prm->ntheta = *it++;
  prm->nphia = *it++;
  prm->numbnd = *it++;
  prm->numang = *it++;
  prm->nptra = *it++;
  prm->natyp = *it++;
  prm->nphb = *it++;
  prm->ifpert = *it++;
  prm->nbper = *it++;
  prm->ngper = *it++;
  prm->ndper = *it++;
  prm->mbper = *it++;
  prm->mgper = *it++;
  prm->mdper = *it++;
  prm->ifbox = *it++;
  prm->nmxrs = *it++;
  prm->ifcap = *it++;
  if (intvec.size() >= 31) prm->numextra = *it++;
  if (intvec.size() >= 32) prm->ncopy = *it++;
}

}  // namespace

namespace amberfile {

int ReadFromFile(const std::string& file_name, AmberPrmtop* prm) {
  FILE* fp = fopen(file_name.c_str(), "r");
  if (fp == NULL) return 1;
  int ret = ReadFromFile(fp, prm);
  fclose(fp);
  return ret;
}

int ReadFromFile(FILE* fp, AmberPrmtop* prm) {
  using namespace std;

  readamb::VersionInfo version;
  readamb::SectionIndex index;

  readamb::BuildIndex(&index, &version, fp);

  prm->version = version.version;
  prm->date = version.date;
  prm->time = version.time;

  readamb::ReadTitle(version, index, &prm->title, fp);
  ReadPointers(version, index, prm, fp);

# define READ(f,n,v,t) do { \
  readamb::ReadData(version, index, (f), (n), (t), &prm->v, fp); } while (0)

  READ("ATOM_NAME", prm->natom, graph, 'A');
  READ("CHARGE", prm->natom, chrg, 'E');
  READ("MASS", prm->natom, amass, 'E');
  READ("ATOM_TYPE_INDEX", prm->natom, iac, 'I');
  READ("NUMBER_EXCLUDED_ATOMS", prm->natom, numex, 'I');
  READ("NONBONDED_PARM_INDEX", prm->ntypes * prm->ntypes, ico, 'I');
  READ("RESIDUE_LABEL", prm->nres, labres, 'A');
  READ("RESIDUE_POINTER", prm->nres, ipres, 'I');
  READ("BOND_FORCE_CONSTANT", prm->numbnd, rk, 'E');
  READ("BOND_EQUIL_VALUE", prm->numbnd, req, 'E');
  READ("ANGLE_FORCE_CONSTANT", prm->numang, tk, 'E');
  READ("ANGLE_EQUIL_VALUE", prm->numang, teq, 'E');
  READ("DIHEDRAL_FORCE_CONSTANT", prm->nptra, pk, 'E');
  READ("DIHEDRAL_PERIODICITY", prm->nptra, pn, 'E');
  READ("DIHEDRAL_PHASE", prm->nptra, phase, 'E');
  READ("SOLTY", prm->natyp, solty, 'E');
  READ("LENNARD_JONES_ACOEF", prm->ntypes * (prm->ntypes + 1) / 2, cn1, 'E');
  READ("LENNARD_JONES_BCOEF", prm->ntypes * (prm->ntypes + 1) / 2, cn2, 'E');
  READ("BONDS_INC_HYDROGEN", prm->nbonh * 3, bh, 'I');
  READ("BONDS_WITHOUT_HYDROGEN", prm->nbona * 3, b, 'I');
  READ("ANGLES_INC_HYDROGEN", prm->ntheth * 4, th, 'I');
  READ("ANGLES_WITHOUT_HYDROGEN", prm->ntheta * 4, t, 'I');
  READ("DIHEDRALS_INC_HYDROGEN", prm->nphih * 5, ph, 'I');
  READ("DIHEDRALS_WITHOUT_HYDROGEN", prm->nphia * 5, p, 'I');
  READ("EXCLUDED_ATOMS_LIST", prm->nnb, natex, 'I');
  READ("HBOND_ACOEF", prm->nphb, asol, 'E');
  READ("HBOND_BCOEF", prm->nphb, bsol, 'E');
  READ("HBCUT", prm->nphb, hbcut, 'E');
  READ("AMBER_ATOM_TYPE", prm->natom, isymbl, 'A');
  READ("TREE_CHAIN_CLASSIFICATION", prm->natom, itree, 'A');
  READ("JOIN_ARRAY", prm->natom, join, 'I');
  READ("IROTAT", prm->natom, irotat, 'I');
  if (prm->ifbox > 0) {}
  if (prm->ifcap > 0) {}
  if (prm->ifpert > 0) {}

#undef READ

  return 0;
}

}  // namespace amberfile
