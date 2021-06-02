#include <algorithm>
#include <vector>
#include "Common.h"
#include "CovalentBondInfo.h"


namespace CovalentBondInfo{
int bond_assign_offset = 0;
int angle_assign_offset = 1;
int torsion_assign_offset = 1;
int improper_assign_offset= 2;

//! 
BondList::BondList()
    : BondArray(0),
      AngleArray(0),
      TorsionArray(0),
      ImproperArray(0)
{}
BondList::BondList(int num_bond, int num_angle, 
                   int num_torsion, int num_improper)
    : BondArray(0),
      AngleArray(0),
      TorsionArray(0),
      ImproperArray(0)
{
  BondArray.resize(num_bond);
  for(int i=0;i<num_bond;i++) {
    BondArray[i].id_of_atom[0] = 0;
    BondArray[i].id_of_atom[1] = 0;
    BondArray[i].typeofbond = 0;
    BondArray[i].shake = 0;                 // int shake: +1
  }
  AngleArray.resize(num_angle);
  for(int i=0;i<num_angle;i++) {
    AngleArray[i].id_of_atom[0] = 0;
    AngleArray[i].id_of_atom[1] = 0;
    AngleArray[i].id_of_atom[2] = 0;
    AngleArray[i].typeofangle = 0;
  }
  TorsionArray.resize(num_torsion);
  for(int i=0;i<num_torsion;i++) {
    TorsionArray[i].id_of_atom[0] = 0;
    TorsionArray[i].id_of_atom[1] = 0;
    TorsionArray[i].id_of_atom[2] = 0;
    TorsionArray[i].id_of_atom[3] = 0;
    TorsionArray[i].typeoftorsion = 0;
    TorsionArray[i].calc14interaction = 0;
  }
  ImproperArray.resize(num_improper);
  for(int i=0;i<num_improper;i++) {
    ImproperArray[i].id_of_atom[0] = 0;
    ImproperArray[i].id_of_atom[1] = 0;
    ImproperArray[i].id_of_atom[2] = 0;
    ImproperArray[i].id_of_atom[3] = 0;
    ImproperArray[i].typeofimproper = 0;
  }
      
}
BondList::BondList(const BondList& bl)
    : BondArray(bl.BondArray),
      AngleArray(bl.AngleArray),
      TorsionArray(bl.TorsionArray),
      ImproperArray(bl.ImproperArray)
{
}
    
void BondList::clear()
{
  BondArray.clear();
  AngleArray.clear();
  TorsionArray.clear();
  ImproperArray.clear();
}
void BondList::print()
{
  std::cout << "Bond " << BondArray.size() << " ";
  for(std::vector<Bond>::size_type b=0;b<BondArray.size();b++){
    std::cout << "(" << BondArray[b].id_of_atom[0];
    std::cout << ":" << BondArray[b].id_of_atom[1];
    std::cout << ":" << BondArray[b].typeofbond;
    std::cout << ":" << BondArray[b].shake << ")";
  }
  std::cout << std::endl;
}
void BondList::nozero_print()
{
  if(BondArray.size()>0){
    print();
  }
}
//! calc size of packed int array
/*!
  @return size of packed array in BYTE
  @note for estimate required size of MPI buffer
*/
int BondList::size_of_packed_int () const
{
  return (5 
          + BondArray.size()*4              // int shake: +1 
          + AngleArray.size()*4 
          + TorsionArray.size()*6 
          + ImproperArray.size()*5)*sizeof(int);
}
//! pack array of int
/*!
  @param[out] bondpack int array which must be enough large
  @note for write to MPI buffer
  @note use size_of_packed_int() to get required size of bondpack 
*/
int BondList::pack_int_array(int* bondpack) const
{
  int ntotal, nb, na, nt, ni;
  nb = BondArray.size();
  na = AngleArray.size();
  nt = TorsionArray.size();
  ni = ImproperArray.size();
  ntotal = 5 + nb*4 + na*4 + nt*6 + ni*5;   // int shake: +1
  //      bondpack = new int[ntotal];
  bondpack[0] = ntotal;
  bondpack[1] = nb;
  bondpack[2] = na;
  bondpack[3] = nt;
  bondpack[4] = ni;
  int bp=5;
  for(int b=0;b<nb;b++){
    bondpack[bp] = BondArray[b].id_of_atom[0];
    bp++;
    bondpack[bp] = BondArray[b].id_of_atom[1];
    bp++;
    bondpack[bp] = BondArray[b].typeofbond;
    bp++;
    bondpack[bp] = BondArray[b].shake;      // int shake: +1
    bp++;
  }
  for(int a=0;a<na;a++){
    bondpack[bp] = AngleArray[a].id_of_atom[0];
    bp++;
    bondpack[bp] = AngleArray[a].id_of_atom[1];
    bp++;
    bondpack[bp] = AngleArray[a].id_of_atom[2];
    bp++;
    bondpack[bp] = AngleArray[a].typeofangle;
    bp++;
  }
  for(int t=0;t<nt;t++){
    bondpack[bp] = TorsionArray[t].id_of_atom[0];
    bp++;
    bondpack[bp] = TorsionArray[t].id_of_atom[1];
    bp++;
    bondpack[bp] = TorsionArray[t].id_of_atom[2];
    bp++;
    bondpack[bp] = TorsionArray[t].id_of_atom[3];
    bp++;
    bondpack[bp] = TorsionArray[t].typeoftorsion;
    bp++;
    bondpack[bp] = TorsionArray[t].calc14interaction ? 1:0;
    bp++;
  }
  for(int i=0;i<ni;i++){
    bondpack[bp] = ImproperArray[i].id_of_atom[0];
    bp++;
    bondpack[bp] = ImproperArray[i].id_of_atom[1];
    bp++;
    bondpack[bp] = ImproperArray[i].id_of_atom[2];
    bp++;
    bondpack[bp] = ImproperArray[i].id_of_atom[3];
    bp++;
    bondpack[bp] = ImproperArray[i].typeofimproper;
    bp++;
  }
  return ntotal;
}
//! unpack array of int
/*!
  @note for read from MPI buffer
*/
void BondList::unpack_int_array(int bondpack[])
{
  //int ntotal = bondpack[0];
  int nb = bondpack[1];
  int na = bondpack[2];
  int nt = bondpack[3];
  int ni = bondpack[4];
  int offset = BondArray.size();
  BondArray.resize(offset+nb);
  int bp = 5;
  for(int b=offset;b<offset+nb;b++){
    BondArray[b].id_of_atom[0] = bondpack[bp];
    bp++;
    BondArray[b].id_of_atom[1] = bondpack[bp];
    bp++;
    BondArray[b].typeofbond = bondpack[bp];
    bp++;
    BondArray[b].shake = bondpack[bp];      // int shake: +1
    bp++;
  }
  offset = AngleArray.size();
  AngleArray.resize(offset+na);
  for(int a=offset;a<offset+na;a++){
    AngleArray[a].id_of_atom[0] = bondpack[bp];
    bp++;
    AngleArray[a].id_of_atom[1] = bondpack[bp];
    bp++;
    AngleArray[a].id_of_atom[2] = bondpack[bp];
    bp++;
    AngleArray[a].typeofangle = bondpack[bp];
    bp++;
  }
  offset = TorsionArray.size();
  TorsionArray.resize(offset+nt);
  for(int t=offset;t<offset+nt;t++){
    TorsionArray[t].id_of_atom[0] = bondpack[bp];
    bp++;
    TorsionArray[t].id_of_atom[1] = bondpack[bp];
    bp++;
    TorsionArray[t].id_of_atom[2] = bondpack[bp];
    bp++;
    TorsionArray[t].id_of_atom[3] = bondpack[bp];
    bp++;
    TorsionArray[t].typeoftorsion = bondpack[bp];
    bp++;
    TorsionArray[t].calc14interaction = (bondpack[bp]==1);
    bp++;
  }
  offset = ImproperArray.size();
  ImproperArray.resize(offset+ni);
  for(int i=offset;i<offset+ni;i++){
    ImproperArray[i].id_of_atom[0] = bondpack[bp];
    bp++;
    ImproperArray[i].id_of_atom[1] = bondpack[bp];
    bp++;
    ImproperArray[i].id_of_atom[2] = bondpack[bp];
    bp++;
    ImproperArray[i].id_of_atom[3] = bondpack[bp];
    bp++;
    ImproperArray[i].typeofimproper = bondpack[bp];
    bp++;
  }
}
void BondList::push_back(BondList addlist)
{
  copy(addlist.BondArray.begin(),addlist.BondArray.end(),back_inserter(BondArray));
  copy(addlist.AngleArray.begin(),addlist.AngleArray.end(),back_inserter(AngleArray));
  copy(addlist.TorsionArray.begin(),addlist.TorsionArray.end(),back_inserter(TorsionArray));
  copy(addlist.ImproperArray.begin(),addlist.ImproperArray.end(),back_inserter(ImproperArray));
}
//! @return size in BYTE
size_t BondList::size()
{
  return sizeof(Bond)*BondArray.size()
      + sizeof(Angle)*AngleArray.size()
      + sizeof(Torsion)*TorsionArray.size()
      + sizeof(Improper)*ImproperArray.size();
}
void BondList::pick_up_bondlist(CovalentBondInfo::BondList& bondlist,
                                const AtomID atomid)
{
  pick_up<Bond>(bondlist.BondArray,BondArray,atomid);
  pick_up<Angle>(bondlist.AngleArray,AngleArray,atomid);
  pick_up<Torsion>(bondlist.TorsionArray,TorsionArray,atomid);
  pick_up<Improper>(bondlist.ImproperArray,ImproperArray,atomid);
}
void BondList::remove_atomid(const AtomID atomid)
{
  {
#if NOBUG_STL
    std::vector<Bond>::iterator end_it = std::remove_if(BondArray.begin(),BondArray.end(),MachAtomID<Bond>(atomid));
#else
    std::vector<Bond>::iterator end_it = std::remove_if(BondArray.begin(),BondArray.end(),MachAtomIDIndex<Bond>(atomid,bond_assign_offset));
#endif
    BondArray.erase(end_it,BondArray.end());
  }
  {
#if NOBUG_STL
    std::vector<Angle>::iterator end_it = std::remove_if(AngleArray.begin(),AngleArray.end(),MachAtomID<Angle>(atomid));
#else
    std::vector<Angle>::iterator end_it = std::remove_if(AngleArray.begin(),AngleArray.end(),MachAtomIDIndex<Angle>(atomid,angle_assign_offset));
#endif
    AngleArray.erase(end_it,AngleArray.end());
  }
  {
#if NOBUG_STL
    std::vector<Torsion>::iterator end_it = std::remove_if(TorsionArray.begin(),TorsionArray.end(),MachAtomID<Torsion>(atomid));
#else
    std::vector<Torsion>::iterator end_it = std::remove_if(TorsionArray.begin(),TorsionArray.end(),MachAtomIDIndex<Torsion>(atomid,torsion_assign_offset));
#endif
    TorsionArray.erase(end_it,TorsionArray.end());
  }
  {
#if NOBUG_STL
    std::vector<Improper>::iterator end_it = std::remove_if(ImproperArray.begin(),ImproperArray.end(),MachAtomID<Improper>(atomid));
#else
    std::vector<Improper>::iterator end_it = std::remove_if(ImproperArray.begin(),ImproperArray.end(),MachAtomIDIndex<Improper>(atomid,improper_assign_offset));
#endif
    ImproperArray.erase(end_it,ImproperArray.end());
  }
}

}


/// Information of Covalent Bonds
CovalentBondList::CovalentBondList() : bond(), angle(), torsion(), improper() 
{
}

CovalentBondList::CovalentBondList(int num) : bond(num), angle(num), torsion(num), improper(num) 
{
}

CovalentBondList::~CovalentBondList() {
  clear();
}

void CovalentBondList::clear() {
  if(!bond.empty()){
    bond.clear();
  }
  if(!angle.empty()){
    angle.clear();
  }
  if(!torsion.empty()){
    torsion.clear();
  }
  if(!improper.empty()){
    improper.clear();
  }
}

void CovalentBondList::merge_bondlistarray(const std::vector<CovalentBondInfo::BondList>& bondlistarray)
{
  for(std::vector<CovalentBondInfo::BondList>::size_type s=0;
      s<bondlistarray.size();s++){
    for(std::vector<CovalentBondInfo::Bond>::size_type b=0;
        b<bondlistarray[s].BondArray.size();b++){
      bond.push_back(bondlistarray[s].BondArray[b]);
    }
    for(std::vector<CovalentBondInfo::Angle>::size_type b=0;
        b<bondlistarray[s].AngleArray.size();b++){
      angle.push_back(bondlistarray[s].AngleArray[b]);
    }
    for(std::vector<CovalentBondInfo::Torsion>::size_type b=0;
        b<bondlistarray[s].TorsionArray.size();b++){
      torsion.push_back(bondlistarray[s].TorsionArray[b]);
    }
    for(std::vector<CovalentBondInfo::Improper>::size_type b=0;
        b<bondlistarray[s].ImproperArray.size();b++){
      improper.push_back(bondlistarray[s].ImproperArray[b]);
    }
  }
  if(DebugLog::verbose>1){
    std::cout << "total CovalentBondInfo Bond " << bond.size();
    std::cout << " Angle " << angle.size();
    std::cout << " Torsion " << torsion.size();
    std::cout << " Improper " << improper.size() << std::endl;
  }
}

void CovalentBondList::pick_up_bondlist(CovalentBondInfo::BondList& bondlist,
                                        const ParticleArray& particlearray,
                                        const TypeRange& typerange,
                                        int assign_type)
{
  int assign_id = CovalentBondInfo::bond_assign_atom(assign_type);
  for(std::vector<CovalentBondInfo::Bond>::size_type bi=0;
      bi<bond.size();bi++){
    int atomid = bond[bi].id_of_atom[assign_id];
    for(int pi=typerange.begin;pi<typerange.end;pi++){
      if(particlearray[pi].atomid==atomid){
        bondlist.BondArray.push_back(bond[bi]);
      }
    }
  }
  assign_id = CovalentBondInfo::angle_assign_atom(assign_type);
  for(std::vector<CovalentBondInfo::Angle>::size_type bi=0;
      bi<angle.size();bi++){
    int atomid = angle[bi].id_of_atom[assign_id];
    for(int pi=typerange.begin;pi<typerange.end;pi++){
      if(particlearray[pi].atomid==atomid){
        bondlist.AngleArray.push_back(angle[bi]);
      }
    }
  }
  assign_id = CovalentBondInfo::torsion_assign_atom(assign_type);
  for(std::vector<CovalentBondInfo::Torsion>::size_type bi=0;
      bi<torsion.size();bi++){
    int atomid = torsion[bi].id_of_atom[assign_id];
    for(int pi=typerange.begin;pi<typerange.end;pi++){
      if(particlearray[pi].atomid==atomid){
        bondlist.TorsionArray.push_back(torsion[bi]);
      }
    }
  }
  assign_id = CovalentBondInfo::improper_assign_atom(assign_type);
  for(std::vector<CovalentBondInfo::Improper>::size_type bi=0;
      bi<improper.size();bi++){
    int atomid = improper[bi].id_of_atom[assign_id];
    for(int pi=typerange.begin;pi<typerange.end;pi++){
      if(particlearray[pi].atomid==atomid){
        bondlist.ImproperArray.push_back(improper[bi]);
      }
    }
  }
}

void CovalentBondList::print() {
  std::cout << "Bond " << std::endl;
  for(size_t s=0;s<bond.size();s++){
    std::cout << "AtomID(" << bond[s].id_of_atom[0] << "," << bond[s].id_of_atom[1] << ") type " << bond[s].typeofbond << " shake " << bond[s].shake << std::endl;
  }
}

int CovalentBondList::size_of_packed_int()
{
  int ntotal, nb, na, nt, ni, nai;
  nb = bond.size();
  na = angle.size();
  nt = torsion.size();
  ni = improper.size();
  nai = atomidset.size();
  ntotal = 6 + nb*3 + na*4 + nt*6 + ni*5 + nai;
  return ntotal;
}

int CovalentBondList::pack_int_array(int* bondpack)
{
  int ntotal, nb, na, nt, ni, nai;
  nb = bond.size();
  na = angle.size();
  nt = torsion.size();
  ni = improper.size();
  nai = atomidset.size();
  ntotal = 6 + nb*3 + na*4 + nt*6 + ni*5 + nai;

  //  printf("CovalentBondList::pack_int_array\n");

  bondpack[0] = ntotal;
  bondpack[1] = nb;
  bondpack[2] = na;
  bondpack[3] = nt;
  bondpack[4] = ni;
  bondpack[5] = nai;
  int bp=6;
  for(int b=0;b<nb;b++){
    bondpack[bp] = bond[b].id_of_atom[0];
    bp++;
    bondpack[bp] = bond[b].id_of_atom[1];
    bp++;
    bondpack[bp] = bond[b].typeofbond;
    bp++;
  }
  for(int a=0;a<na;a++){
    bondpack[bp] = angle[a].id_of_atom[0];
    bp++;
    bondpack[bp] = angle[a].id_of_atom[1];
    bp++;
    bondpack[bp] = angle[a].id_of_atom[2];
    bp++;
    bondpack[bp] = angle[a].typeofangle;
    bp++;
  }
  for(int t=0;t<nt;t++){
    bondpack[bp] = torsion[t].id_of_atom[0];
    bp++;
    bondpack[bp] = torsion[t].id_of_atom[1];
    bp++;
    bondpack[bp] = torsion[t].id_of_atom[2];
    bp++;
    bondpack[bp] = torsion[t].id_of_atom[3];
    bp++;
    bondpack[bp] = torsion[t].typeoftorsion;
    bp++;
    bondpack[bp] = torsion[t].calc14interaction ? 1:0;
    bp++;
  }
  for(int i=0;i<ni;i++){
    bondpack[bp] = improper[i].id_of_atom[0];
    bp++;
    bondpack[bp] = improper[i].id_of_atom[1];
    bp++;
    bondpack[bp] = improper[i].id_of_atom[2];
    bp++;
    bondpack[bp] = improper[i].id_of_atom[3];
    bp++;
    bondpack[bp] = improper[i].typeofimproper;
    bp++;
  }
  for(std::set<int>::iterator it=atomidset.begin();it!=atomidset.end();it++){
    bondpack[bp] = (*it);
    bp++;
  }
  return ntotal;
 
}

void CovalentBondList::unpack_int_array(int bondpack[])
{
  int nb = bondpack[1];
  int na = bondpack[2];
  int nt = bondpack[3];
  int ni = bondpack[4];
  int nai = bondpack[5];
  int offset = bond.size();
  bond.resize(offset+nb);
  int bp = 6;
  for(int b=offset;b<offset+nb;b++){
    bond[b].id_of_atom[0] = bondpack[bp];
    bp++;
    bond[b].id_of_atom[1] = bondpack[bp];
    bp++;
    bond[b].typeofbond = bondpack[bp];
    bp++;
  }
  offset = angle.size();
  angle.resize(offset+na);
  for(int a=offset;a<offset+na;a++){
    angle[a].id_of_atom[0] = bondpack[bp];
    bp++;
    angle[a].id_of_atom[1] = bondpack[bp];
    bp++;
    angle[a].id_of_atom[2] = bondpack[bp];
    bp++;
    angle[a].typeofangle = bondpack[bp];
    bp++;
  }
  offset = torsion.size();
  torsion.resize(offset+nt);
  for(int t=offset;t<offset+nt;t++){
    torsion[t].id_of_atom[0] = bondpack[bp];
    bp++;
    torsion[t].id_of_atom[1] = bondpack[bp];
    bp++;
    torsion[t].id_of_atom[2] = bondpack[bp];
    bp++;
    torsion[t].id_of_atom[3] = bondpack[bp];
    bp++;
    torsion[t].typeoftorsion = bondpack[bp];
    bp++;
    torsion[t].calc14interaction = (bondpack[bp]==1);
    bp++;
  }
  offset = improper.size();
  improper.resize(offset+ni);
  for(int i=offset;i<offset+ni;i++){
    improper[i].id_of_atom[0] = bondpack[bp];
    bp++;
    improper[i].id_of_atom[1] = bondpack[bp];
    bp++;
    improper[i].id_of_atom[2] = bondpack[bp];
    bp++;
    improper[i].id_of_atom[3] = bondpack[bp];
    bp++;
    improper[i].typeofimproper = bondpack[bp];
    bp++;
  }
  for(int i=0;i<nai;i++){
    atomidset.insert(bondpack[bp]);
    bp++;
  }
}

CovalentBondList& CovalentBondList::append(const CovalentBondList& cbl) {
  this->bond.insert(this->bond.end(), cbl.bond.begin(), cbl.bond.end());
  this->angle.insert(this->angle.end(), cbl.angle.begin(), cbl.angle.end());
  this->torsion.insert(this->torsion.end(),
                       cbl.torsion.begin(), cbl.torsion.end());
  this->improper.insert(this->improper.end(),
                        cbl.improper.begin(), cbl.improper.end());
  return *this;
}

CovalentBondList& CovalentBondList::swap(CovalentBondList& cbl) {
  this->bond.swap(cbl.bond);
  this->angle.swap(cbl.angle);
  this->torsion.swap(cbl.torsion);
  this->improper.swap(cbl.improper);
  return *this;
}

/// Information of Covalent Bonds
CovalentBondParameterList::CovalentBondParameterList()
    : bond(),
      angle(),
      torsion(),
      improper() {
}

CovalentBondParameterList::CovalentBondParameterList(int num)
    : bond(num),
      angle(num),
      torsion(num),
      improper(num) {
}

CovalentBondParameterList::~CovalentBondParameterList() {
  clear();
}

void CovalentBondParameterList::clear() {
  if(!bond.empty()){
    bond.clear();
  }
  if(!angle.empty()){
    angle.clear();
  }
  if(!torsion.empty()){
    torsion.clear();
  }
  if(!improper.empty()){
    improper.clear();
  }
}

void CovalentBondParameterList::print() const {
  std::cout << "Bond Parameter";
  for(std::vector<CovalentBondInfo::BondParameter>::size_type b=0;
      b<bond.size();b++){
    std::cout << " " << b << "(" << bond[b].force_constant << "," << bond[b].equilibrium_length << ")";
  }
  std::cout << std::endl;
}

int CovalentBondParameterList::size_of_packed_double(int nums[4])
{
  int ntotal, nb, na, nt, ni;
  nb = bond.size();
  na = angle.size();
  nt = torsion.size();
  ni = improper.size();
  ntotal = nb*2 + na*2 + nt*3 + ni*3;
  nums[0] = nb;
  nums[1] = na;
  nums[2] = nt;
  nums[3] = ni;
  return ntotal;
}

int CovalentBondParameterList::pack_double_array(int nums[4], double* bondpack)
{
  int ntotal, nb, na, nt, ni;
  nb = bond.size();
  na = angle.size();
  nt = torsion.size();
  ni = improper.size();
  ntotal = nb*2 + na*2 + nt*3 + ni*3;
  nums[0] = nb;
  nums[1] = na;
  nums[2] = nt;
  nums[3] = ni;
  int bp=0;
  for(int b=0;b<nb;b++){
    bondpack[bp] = bond[b].force_constant;
    bp++;
    bondpack[bp] = bond[b].equilibrium_length;
    bp++;
  }
  for(int a=0;a<na;a++){
    bondpack[bp] = angle[a].force_constant;
    bp++;
    bondpack[bp] = angle[a].equilibrium_angle;
    bp++;
  }
  for(int t=0;t<nt;t++){
    bondpack[bp] = torsion[t].force_constant;
    bp++;
    bondpack[bp] = torsion[t].periodicity;
    bp++;
    bondpack[bp] = torsion[t].phase;
    bp++;
  }
  for(int i=0;i<ni;i++){
    bondpack[bp] = improper[i].force_constant;
    bp++;
    bondpack[bp] = improper[i].periodicity;
    bp++;
    bondpack[bp] = improper[i].phase;
    bp++;
  }
  return ntotal;
}

void CovalentBondParameterList::unpack_double_array(int nums[4], double bondpack[])
{
  int nb = nums[0];
  int na = nums[1];
  int nt = nums[2];
  int ni = nums[3];
  int offset = bond.size();
  bond.resize(offset+nb);
  int bp = 0;
  for(int b=offset;b<offset+nb;b++){
    bond[b].force_constant = bondpack[bp];
    bp++;
    bond[b].equilibrium_length = bondpack[bp];
    bp++;
  }
  offset = angle.size();
  angle.resize(offset+na);
  for(int a=offset;a<offset+na;a++){
    angle[a].force_constant = bondpack[bp];
    bp++;
    angle[a].equilibrium_angle = bondpack[bp];
    bp++;
  }
  offset = torsion.size();
  torsion.resize(offset+nt);
  for(int t=offset;t<offset+nt;t++){
    torsion[t].force_constant = bondpack[bp];
    bp++;
    torsion[t].periodicity = bondpack[bp];
    bp++;
    torsion[t].phase = bondpack[bp];
    bp++;
  }
  offset = improper.size();
  improper.resize(offset+ni);
  for(int i=offset;i<offset+ni;i++){
    improper[i].force_constant = bondpack[bp];
    bp++;
    improper[i].periodicity = bondpack[bp];
    bp++;
    improper[i].phase = bondpack[bp];
    bp++;
  }
}
