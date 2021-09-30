#ifndef H_USER_DEFINED_CLASS
#define H_USER_DEFINED_CLASS

#include <particle_mesh.hpp>
#include <param_fdps.h>
#include <param.h>

const PS::F64 RCUT = PS::PM::CUTOFF_RADIUS / (PS::F64)SIZE_OF_MESH;

// =================
// Force class
// =================
class Force{
public:
  PS::F64vec acc;
  PS::F64    pot;
  void clear(){
    acc = 0.0;
    pot = 0.0;
  }
};

// =========================================
// Full particle and essential particle class
// =========================================
class FP{
public:
  PS::S64 id;
  PS::S64 type;

  PS::F64vec pos, vel, acc;
  PS::F64 charge,pot;

  PS::F64 getRsearch() const {
    return RCUT;
  }
  PS::F64vec getPos() const {
    return pos;
  }
  void setPos(const PS::F64vec _pos){
    pos = _pos;
  }
  void copyFromForce(const Force& force){
    acc = force.acc;
    pot = force.pot;
  }
  PS::F64 getChargeParticleMesh() const {
    return charge;
  }
  void copyFromForceParticleMesh(const Force& f){
    acc = f.acc;
    pot = f.pot;
  }
};

class EP{
public:
#ifdef ENABLE_SINGLE
  PS::F32vec pos;
  PS::S32 charge;
#else
  PS::F64vec pos;
  PS::F64 charge;
#endif

  PS::F64 getRSearch()   const { return RCUT; }
  PS::F64vec getPos()    const { return pos; }
  PS::F64    getCharge() const { return charge; }
  void copyFromFP(const FP & fp){
    pos  = fp.pos;
    charge = fp.charge;
  }
  void setPos(const PS::F64vec3 _pos){
    pos = _pos;
  }
};

class MomentQuadrupoleCutoff{
public:
  PS::S32    n_ptcl;
  PS::F64vec pos;
  PS::F64    charge;
  PS::F64vec dipole;
  PS::F64mat quadrupole;
  PS::F64ort vertex_out_;
  MomentQuadrupoleCutoff(){
    n_ptcl = 0;
    pos = 0.0;
    charge = 0.0;
    dipole = 0.0;
    quadrupole = 0.0;
    vertex_out_.init();
  }
  MomentQuadrupoleCutoff(const PS::S32 n,
			 const PS::F64 c,
			 const PS::F64vec & p,
			 const PS::F64vec & d,
			 const PS::F64mat & q,
			 const PS::F64ort & v_out){
    n_ptcl = n;
    charge = c;
    pos = p;
    dipole = d;
    quadrupole = q;
    vertex_out_ = v_out;
  }
  void init(){
    n_ptcl = 0;
    charge = 0.0;
    pos = 0.0;
    dipole = 0.0;
    quadrupole = 0.0;
    vertex_out_.init();
  }
  PS::F64vec getPos() const {
    return pos;
  }
  PS::F64 getCharge() const {
    return charge;
  }
  PS::F64ort getVertexOut() const { return vertex_out_; }

  template<class Tepj>
  void accumulateAtLeaf(Tepj & epj){
    n_ptcl++;
    charge += epj.getCharge();
    pos += epj.getPos();
    vertex_out_.merge(epj.getPos(), epj.getRSearch());
  }
  template<class Tepj>
  void accumulateAtLeaf2(Tepj & epj){
    PS::F64 ctmp = epj.getCharge();
    PS::F64vec ptmp = epj.getPos() - pos;
    dipole += ctmp * ptmp;
    PS::F64 cx = ctmp * ptmp.x;
    PS::F64 cy = ctmp * ptmp.y;
    PS::F64 cz = ctmp * ptmp.z;
    quadrupole.xx += cx * ptmp.x;
    quadrupole.yy += cy * ptmp.y;
    quadrupole.zz += cz * ptmp.z;
    quadrupole.xy += cx * ptmp.y;
    quadrupole.xz += cx * ptmp.z;
    quadrupole.yz += cy * ptmp.z;
}
  void set(){
    pos = pos / (PS::F64)n_ptcl;
  }
  void accumulate(const MomentQuadrupoleCutoff & mom){
    n_ptcl += mom.n_ptcl;
    charge += mom.charge;
    pos += mom.n_ptcl * mom.pos;
    vertex_out_.merge(mom.vertex_out_);
  }
  void accumulate2(const MomentQuadrupoleCutoff & mom){
    PS::F64 ctmp = mom.charge;
    PS::F64vec ptmp = mom.pos - pos;
    dipole += ctmp * ptmp + mom.dipole;
    PS::F64 cx = ctmp * ptmp.x;
    PS::F64 cy = ctmp * ptmp.y;
    PS::F64 cz = ctmp * ptmp.z;
    quadrupole.xx += cx * ptmp.x + mom.quadrupole.xx;
    quadrupole.yy += cy * ptmp.y + mom.quadrupole.yy;
    quadrupole.zz += cz * ptmp.z + mom.quadrupole.zz;
    quadrupole.xy += cx * ptmp.y + mom.quadrupole.xy;
    quadrupole.xz += cx * ptmp.z + mom.quadrupole.xz;
    quadrupole.yz += cy * ptmp.z + mom.quadrupole.yz;
}
  // for DEBUG 
  void dump(std::ostream & fout = std::cout) const {
    fout<<"charge="<<charge<<std::endl;
    fout<<"pos="<<pos<<std::endl;
    fout<<"vertex_out_="<<vertex_out_<<std::endl;
  }
};

class SPJQuadrupoleCutoff{
public:
  PS::S32 n_ptcl;
#ifdef ENABLE_SINGLE
  PS::F32    charge;
  PS::F32vec pos;
  PS::F32vec dipole;
  PS::F32mat quadrupole;
#else
  PS::F64 charge;
  PS::F64vec pos;
  PS::F64vec dipole;
  PS::F64mat quadrupole;
#endif
  template<class Tmom>
  void copyFromMoment(const Tmom & mom){
    n_ptcl = mom.n_ptcl;
    charge = mom.charge;
    pos = mom.pos;
    dipole = mom.dipole;
    quadrupole = mom.quadrupole;
  }
  void clear(){
    n_ptcl = 0;
    charge = 0.0;
    pos = 0.0;
    dipole = 0.0;
    quadrupole = 0.0;
  }
  PS::F64 getCharge() const {
    return charge;
  }
  PS::F64vec getPos() const {
    return pos;
  }
  void setPos(const PS::F64vec & pos_new) {
    pos = pos_new;
  }
  MomentQuadrupoleCutoff convertToMoment() const {
    return MomentQuadrupoleCutoff(n_ptcl,charge, pos, dipole, quadrupole, PS::F64ort(0.0,0.0));
  }
};

#endif
