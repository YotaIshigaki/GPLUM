#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<iostream>
#include<fstream>
#include<sys/time.h>
#include<particle_simulator.hpp>

/*Functions of Force, velocity, position and so on*/

//#define cdvfile_input
#define cdvfile_output

class FileHeader{
public:
  PS::F32 boxdh;
  PS::F32 readAscii(FILE *fp){
#ifdef cdvfile_input
    fscanf(fp,"'box_sx=%*f box_sy=%*f box_sz=%*f box_ex=%lf box_ey=%*f box_ez=%*f\n'box_wt=0.007 r1=0.5\n",
    	   &boxdh);
#endif
    return -1;
  }
  void writeAscii(FILE *fp)const{
#ifdef cdvfile_output
    fprintf(fp,"'box_sx=%lf box_sy=%lf box_sz=%lf box_ex=%lf box_ey=%lf box_ez=%lf\n'box_wt=0.007 r1=0.5\n",
	    -boxdh,-boxdh,-boxdh,boxdh,boxdh,boxdh);
#endif
  }
};

class ForceLJ{
public:
  PS::F64vec acc;
  PS::F32 pot;
  //PS::F32 virial;
  //int neighbor[15];
  void clear(){
    acc = 0.0;
    pot = 0.0;
    // for(int i=0; i<15; i++){
    //   neighbor[i] = -1;
    // }
  }
};

class FPLJ{
public:
  PS::S32 id;
  PS::F32vec pos;
  PS::F32vec vel;
  PS::F64vec acc;
  PS::F32 pot;
  //PS::F32 virial;
  PS::F32 search_radius;
  //  int neighbor[15];
  PS::F32 getRsearch() const{
    return this->search_radius;
  }
  void setPos(const PS::F32vec & p) { pos = p; }
  PS::F32vec getPos() const { return pos; } //必須メンバ関数
  void copyFromForce(const ForceLJ & force){ //必須メンバ関数
    acc = force.acc;
    pot = force.pot;
    //virial = force.virial;
    // for(int i=0; i<15; i++){
    //   neighbor[i] = force.neighbor[i];
    // }
  }
  //*  cdv file  out put  *//
#ifdef cdvfile_output
  void writeAscii(FILE* fp) const{
    if(pot<-1.5){
      fprintf(fp, "%ld %ld %17.16e %17.16e %17.16e\n",
	      this->id, this->pos.x, this->pos.y, this->pos.z); 
    }
  }
#else
  void writeAscii(FILE* fp) const{
    fprintf(fp, "%lld %lld %lf %lf %lf %lf %lf %lf %lf\n",
            this->id,
	    this->pos.x, this->pos.y, this->pos.z,
	    this->vel.x, this->vel.y, this->vel.z,
	    this->pot); 
  }
#endif
#ifdef cdvfile_input
  void readAscii(FILE* fp){
    fscanf(fp, "%lld %lld %lf %lf %lf\n",
	   &this->id, &this->pos.x, &this->pos.y, &this->pos.z);
  }
#else
  void readAscii(FILE* fp){
    // fscanf(fp, "%d %d %f %f %f %f %f %f\n",
    // 	   &this->id,
    // 	   &this->pos.x, &this->pos.y, &this->pos.z,
    // 	   &this->vel.x, &this->vel.y, &this->vel.z
    // 	   );
    fscanf(fp, "%d %f %f %f %f %f %f\n",
	   &id,
	   &pos.x, &pos.y, &pos.z,
	   &vel.x, &vel.y, &vel.z);
    // printf("%d %d %f %f %f %f %f %f\n",
    // 	   pos.x, pos.y, pos.z,
    // 	   vel.x, vel.y, vel.z);
  }
#endif
};

class EPILJ{
public:
  PS::S32 id;
  PS::F32vec pos;
  PS::F32vec getPos() const { return pos;}
  void copyFromFP(const FPLJ & fp){
    pos = fp.pos;
    id = fp.id;
  }
};

//PS::F64 EPILJ::eps = 1.0/32.0;
class EPJLJ{
public:

  PS::S32 id;
  PS::F32vec pos;
  PS::F32 search_radius;
  PS::F32 getRSearch() const{
    return this->search_radius;
  }
  void copyFromFP(const FPLJ & fp){
    pos = fp.pos;
    id = fp.id;
    search_radius = fp.search_radius;
  }
  PS::F32vec getPos() const { return pos; }
  void setPos(const PS::F32vec & pos_new){ pos = pos_new;}
  //PS::F64 getCharge() const { return mass; }
};

#define TUNING
#ifdef TUNING
static unsigned long flop = 0;
#endif
struct CalcForceEpEp{
  void operator () (const EPILJ * ep_i,
                    const PS::S32 n_ip,
                    const EPJLJ * ep_j,
                    const PS::S32 n_jp,
                    ForceLJ * force){
    for(PS::S32 i = 0; i < n_ip; i++){
      PS::F32vec pos_i = ep_i[i].pos;
      PS::F64vec acc_i = 0.0d;
      PS::F32    pot_i = 0.0f;
      PS::S32    id_i = ep_i[i].id;
      for(PS::S32 j = 0; j < n_jp; j++){
        PS::S32 id_j = ep_j[j].id;
	//if( id_i == id_j ) continue;
	PS::F32vec rij = pos_i - ep_j[j].pos;
	PS::F32 r2 = rij * rij;
#ifdef TUNING
	flop += 8;
#endif
        if(r2 <= cut_off2 && r2!=0.0){
	  const PS::F32 r2_inv = 1.f/r2;
	  const PS::F32 r6_inv =  r2_inv*r2_inv*r2_inv;
	  const PS::F32 r12_inv = r6_inv*r6_inv;
	  pot_i += 4.f*(r12_inv - r6_inv);
	  const PS::F32 dphi = (48.f*r12_inv - 24.f*r6_inv)*r2_inv;
	  acc_i += dphi* rij;
#ifdef TUNING
	  flop += 17;
#endif
	}
      }
      force[i].acc += acc_i;
      force[i].pot += pot_i;
#ifdef TUNING
      flop += 2;
#endif
      //force[i].virial += virial_i;
    }
  }
  CalcForceEpEp(const PS::F32 _cut_off2) : cut_off2(_cut_off2) {}
private:
  PS::F32 cut_off2;
};


template<class Tpsys>
void set_initial_pos(Tpsys & system,
		     const int unit_num,
		     const double unit_length,
		     const PS::F32vec sideh
		     ){
  const int np = system.getNumberOfParticleLocal();
  int i = 0;
  for(int ix = 0; ix<unit_num; ix++){
    for(int iy = 0; iy<unit_num; iy++){
      for(int iz = 0; iz<unit_num; iz++){
	system[i].id = i;
	system[i].pos.x = -sideh.x + unit_length*ix;
	system[i].pos.y = -sideh.y + unit_length*iy;
	system[i].pos.z = -sideh.z + unit_length*iz;
	i++;
	system[i].id = i;
	system[i].pos.x = -sideh.x + unit_length*(ix+0.5);
	system[i].pos.y = -sideh.y + unit_length*(iy+0.5);
	system[i].pos.z = -sideh.z + unit_length*iz;
	i++;
	system[i].id = i;
	system[i].pos.x = -sideh.x + unit_length*(ix+0.5);
	system[i].pos.y = -sideh.y + unit_length*iy;
	system[i].pos.z = -sideh.z + unit_length*(iz+0.5);
	i++;
	system[i].id = i;
	system[i].pos.x = -sideh.x + unit_length*ix;
	system[i].pos.y = -sideh.y + unit_length*(iy+0.5);
	system[i].pos.z = -sideh.z + unit_length*(iz+0.5);
	i++;
      }
    }
  }
}

template<class Tpsys>
void velocity_scaling(Tpsys &system,const PS::F64 T){
  const int np = system.getNumberOfParticleLocal();
  PS::F64vec sum = 0.0;
  for(int i = 0; i < np; i++) sum += system[i].vel * system[i].vel;
  const PS::F64 s = sqrt((3.*np*T)/(sum.x+sum.y+sum.z));
  for(int i=0;i<np;i++) system[i].vel *= s;

  sum = 0.0;
  for(int i = 0; i < np; i++) sum += system[i].vel;
  sum /= (PS::F64)np;
  for(int i = 0; i < np; i++) system[i].vel -= sum;
}

template<class Tpsys>
void set_initial_vel(Tpsys &system)
{
  const int np = system.getNumberOfParticleLocal();
  srand(1);
  /* 速度の設定 */
  for(int i = 0; i < np; i++)
    for(int j = 0; j < 3; j++)
      system[i].vel[j] = (float)rand() / (float)RAND_MAX;
  velocity_scaling(system,1.0);
}
template<class Tpsys>
void PBC(Tpsys &system,
	 const PS::F32vec sideh,
	 const PS::F32vec side){
  const int np = system.getNumberOfParticleLocal();
  for(int i=0; i<np ;i++){
    if(system[i].pos.x >= sideh.x) system[i].pos.x -= side.x;
    if(system[i].pos.y >= sideh.y) system[i].pos.y -= side.y;
    if(system[i].pos.z >= sideh.z) system[i].pos.z -= side.z;
    if(system[i].pos.x < -sideh.x) system[i].pos.x += side.x;
    if(system[i].pos.y < -sideh.y) system[i].pos.y += side.y;
    if(system[i].pos.z < -sideh.z) system[i].pos.z += side.z;
  }
}

template<class Tpsys>
void calc_energy(const Tpsys &system,PS::F64 &pot_tot,PS::F64 &kin_tot){
  const int np = system.getNumberOfParticleLocal();
  PS::F64 pot = 0.0,kin = 0.0;
  PS::F64vec sum = 0.0;
  pot_tot = 0.0,kin_tot = 0.0;

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for reduction(+:pot,kin)
#endif
  for(int i = 0; i < np; i++){
    pot += system[i].pot;
    kin += system[i].vel.x * system[i].vel.x + system[i].vel.y * system[i].vel.y + system[i].vel.z * system[i].vel.z;
  }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&pot, &pot_tot, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&kin, &kin_tot, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
#else
  pot_tot = 0.5*pot;
  kin_tot = 0.5*kin;
#endif
}

