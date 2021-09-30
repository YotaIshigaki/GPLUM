#pragma once
class FileHeader{
public:
    PS::S64 n_body;
    PS::F64 time;
    PS::S32 readAscii(FILE * fp) {
        fscanf(fp, "%lf\n", &time);
        fscanf(fp, "%lld\n", &n_body);
        return n_body;
    }
    void writeAscii(FILE* fp) const {
        fprintf(fp, "%e\n", time);
        fprintf(fp, "%lld\n", n_body);
    }
};

class Force{
public:
  PS::F32vec acc;
  PS::F32    pot;
  void clear(){ acc=0.f; pot=0.f;}
};

class FPGrav{
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64    pot;    

    PS::F64 eps;

    PS::F64vec getPos() const {
        return pos;
    }

    PS::F64 getCharge() const {
        return mass;
    }

    void copyFromForce(const Force & force) {
        acc = force.acc;
        pot = force.pot;
    }

    void writeAscii(FILE* fp) const {
        fprintf(fp, "%lld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
                this->id, this->mass,
                this->pos.x, this->pos.y, this->pos.z,
                this->vel.x, this->vel.y, this->vel.z);
    }

    void readAscii(FILE* fp) {
        fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
               &this->id, &this->mass,
               &this->pos.x, &this->pos.y, &this->pos.z,
               &this->vel.x, &this->vel.y, &this->vel.z);
        }

};

class EPI{
public:
  PS::F64vec pos;
  PS::F32    eps2;
  PS::F64vec getPos() const { return pos; }
  void copyFromFP(const FPGrav &fp){
    pos  = fp.pos;
    eps2 = fp.eps*fp.eps;
  }
};
class EPJ {
public:
  PS::F64vec pos;
  PS::F32    mass;
  PS::F32    eps2;
  PS::F64vec getPos() const { return pos; }
  PS::F64 getCharge() const { return mass; }
  void copyFromFP(const FPGrav &fp){
    pos  = fp.pos;
    mass = fp.mass;
    eps2 = fp.eps*fp.eps;
  }
};

void CalcGravityEP(const EPI * ep_i,
		   const PS::S32 n_ip,
		   const EPJ * ep_j,
		   const PS::S32 n_jp,
		   Force * force) {
  PS::F32vec xiloc[n_ip];
  PS::F32vec xjloc[n_jp];
  for(int i=0;i<n_ip;i++) xiloc[i] = ep_i[i].pos - ep_i[0].pos;
  for(int j=0;j<n_jp;j++) xjloc[j] = ep_j[j].pos - ep_i[0].pos;

  for(PS::S32 i = 0; i < n_ip; i++){
    PS::F32vec xi = xiloc[i];
    PS::F32 eps2 = ep_i[i].eps2;
    PS::F32vec ai = 0.0;
    PS::F32 poti = 0.0;
    for(PS::S32 j = 0; j < n_jp; j++){
      PS::F32vec rij    = xi - xjloc[j];
      PS::F32    r3_inv = rij * rij + eps2 + ep_j[j].eps2;
      PS::F32    r_inv  = 1.0/sqrt(r3_inv);
      r3_inv  = r_inv * r_inv;
      r_inv  *= ep_j[j].getCharge();
      r3_inv *= r_inv;
      ai     -= r3_inv * rij;
      poti   -= r_inv;
    }
    force[i].acc += ai;
    force[i].pot += poti;
  }
}

void CalcGravitySP(const EPI * ep_i,
		   const PS::S32 n_ip,
		   const PS::SPJMonopole * ep_j,
		   const PS::S32 n_jp,
		   Force * force) {
  PS::F32vec xiloc[n_ip];
  PS::F32vec xjloc[n_jp];
  for(int i=0;i<n_ip;i++) xiloc[i] = ep_i[i].pos - ep_i[0].pos;
  for(int j=0;j<n_jp;j++) xjloc[j] = ep_j[j].pos - ep_i[0].pos;

  for(PS::S32 i = 0; i < n_ip; i++){
    PS::F32vec xi = xiloc[i];
    PS::F32 eps2 = ep_i[i].eps2;
    PS::F32vec ai = 0.0;
    PS::F32 poti = 0.0;
    for(PS::S32 j = 0; j < n_jp; j++){
      PS::F32vec rij    = xi - xjloc[j];
      PS::F32    r3_inv = rij * rij + eps2;
      PS::F32    r_inv  = 1.0/sqrt(r3_inv);
      r3_inv  = r_inv * r_inv;
      r_inv  *= ep_j[j].getCharge();
      r3_inv *= r_inv;
      ai     -= r3_inv * rij;
      poti   -= r_inv;
    }
    force[i].acc += ai;
    force[i].pot += poti;
  }
}

#include"kernel.h"
