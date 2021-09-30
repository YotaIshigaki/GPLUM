#ifndef __TREEPM_HPP__
#define __TREEPM_HPP__

#include <sys/times.h>
#include <sys/time.h>
#include <param_fdps.h>
#include <mpi.h>

#include "run_param.hpp"
#include "cosmology.hpp"

#define TINY (1.0e-30)

template<class T>
T reverse_endian(T value){
    char * first = reinterpret_cast<char*>(&value);
    char * last = first + sizeof(T);
    std::reverse(first, last);
    return value;
}

class Result_treepm {
public:
  PS::F32vec acc;
  PS::F32    pot;

  void clear() {
    acc = 0.0;
    pot = 0.0;
  }
};

class FPtreepm {
public:
  PS::S64    id;
  PS::F64    mass;
  static PS::F64    eps;
  PS::F64vec pos;
  PS::F64vec vel;
  PS::F64vec acc;
  PS::F64vec acc_pm;
    /*
  PS::S32    id;
  PS::F32    mass;
  static PS::F32    eps;
  PS::F32vec pos;
  PS::F32vec vel;
  PS::F32vec acc;
  PS::F32vec acc_pm;
    */
  PS::F64vec getPos() const {
    return pos;
  }

  PS::F64 getCharge() const {
    return mass;
  }

  void copyFromFP(const FPtreepm & fp) {
    this->mass = fp.mass;
    this->pos  = fp.pos;
  }

  void copyFromForce(const Result_treepm & force) {
    this->acc = force.acc;
  }

  constexpr PS::F64 getRSearch() const {
      PS::F64 rcut = 3.0/SIZE_OF_MESH;
      return rcut;
  }

  void setPos(const PS::F64vec pos_new) {
    this->pos = pos_new;
  }

  PS::F64 getChargeParticleMesh() const {
    return this->mass;
  }

  void copyFromForceParticleMesh(const PS::F64vec & acc_pm) {
    this->acc_pm = acc_pm;
  }

  void clear() {
    this->acc = 0.0;
  }

  void writeParticleBinary(FILE *fp) {
    int count;
    count = 0;
    count += fwrite(&mass,   sizeof(PS::F32),1,fp);
    count += fwrite(&eps,    sizeof(PS::F32),1,fp);
    count += fwrite(&pos[0], sizeof(PS::F64),1,fp);
    count += fwrite(&pos[1], sizeof(PS::F64),1,fp);
    count += fwrite(&pos[2], sizeof(PS::F64),1,fp);
    count += fwrite(&vel[0], sizeof(PS::F64),1,fp);
    count += fwrite(&vel[1], sizeof(PS::F64),1,fp);
    count += fwrite(&vel[2], sizeof(PS::F64),1,fp);
  }

  int readParticleBinary(FILE *fp) {
    int count;
    count += fread(&mass,   sizeof(PS::F32),1,fp);
    count += fread(&eps,    sizeof(PS::F32),1,fp);
    count += fread(&pos[0], sizeof(PS::F64),1,fp);
    count += fread(&pos[1], sizeof(PS::F64),1,fp);
    count += fread(&pos[2], sizeof(PS::F64),1,fp);
    count += fread(&vel[0], sizeof(PS::F64),1,fp);
    count += fread(&vel[1], sizeof(PS::F64),1,fp);
    count += fread(&vel[2], sizeof(PS::F64),1,fp);
    return count;
  }



    void readBinary(FILE *fp){
        PS::F32 x, y, z, vx, vy, vz, m;
        PS::S32 i;
        fread(&x, 4, 1, fp);
        fread(&vx, 4, 1, fp);
        fread(&y, 4, 1, fp);
        fread(&vy, 4, 1, fp);
        fread(&z, 4, 1, fp);
        fread(&vz, 4, 1, fp);
        fread(&m,   4, 1, fp);
        fread(&i,   4, 1, fp);

/*
        pos.x = reverse_endian(x) / (64.0/h);
        pos.y = reverse_endian(y) / (64.0/h);
        pos.z = reverse_endian(z) / (64.0/h);
        vel.x = reverse_endian(vx) / 64.0;
        vel.y = reverse_endian(vy) / 64.0;
        vel.z = reverse_endian(vz) / 64.0;
        mass = reverse_endian(m) / (6.097e13/h);
*/
/*
        pos.x = reverse_endian(x) / 64.0;
        pos.y = reverse_endian(y) / 64.0;
        pos.z = reverse_endian(z) / 64.0;
        vel.x = reverse_endian(vx) / 1600.0;
        vel.y = reverse_endian(vy) / 1600.0;
        vel.z = reverse_endian(vz) / 1600.0;
        mass = reverse_endian(m) / 1.914e16;
*/

        pos.x = reverse_endian(x) / 64.0;
        pos.y = reverse_endian(y) / 64.0;
        pos.z = reverse_endian(z) / 64.0;
        vel.x = reverse_endian(vx) / 3200.0;
        vel.y = reverse_endian(vy) / 3200.0;
        vel.z = reverse_endian(vz) / 3200.0;
        mass = reverse_endian(m) / 1.524e17;

        id = reverse_endian(i);
    }

    void writeBinary(FILE *fp){
        PS::F32vec x = pos;
        PS::F32vec v = vel;
        PS::S32 i = id;
        fwrite(&x.x,   sizeof(PS::F32), 1, fp);
        fwrite(&v.x,   sizeof(PS::F32), 1, fp);
        fwrite(&x.y,   sizeof(PS::F32), 1, fp);
        fwrite(&v.y,   sizeof(PS::F32), 1, fp);
        fwrite(&x.z,   sizeof(PS::F32), 1, fp);
        fwrite(&v.z,   sizeof(PS::F32), 1, fp);
        fwrite(&mass,   sizeof(PS::F32), 1, fp);
        fwrite(&i,   sizeof(PS::S32), 1, fp);
    }

    void writeAscii(FILE *fp){
        fprintf(fp, "%lld\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
                this->id, this->mass, this->pos.x, this->pos.y, this->pos.z, this->vel.x, this->vel.y, this->vel.z, 
                this->acc.x, this->acc.y, this->acc.z, this->acc_pm.x, this->acc_pm.y, this->acc_pm.z);
    }
    
    PS::F64 calcDtime(run_param &this_run) {
        PS::F64 dtime_v, dtime_a, dtime;
        PS::F64 vnorm, anorm;
#if 0
        vnorm = sqrt(SQR(this->vel))+TINY;
        return 0.5 * this->eps/vnorm;
#else
        vnorm = sqrt(SQR(this->vel))+TINY;
        anorm = sqrt(SQR(this->acc+this->acc_pm))+TINY;
        dtime_v = this->eps/vnorm;
        dtime_a = sqrt(this->eps/anorm)*CUBE(this_run.anow);
        dtime = fmin(0.5*dtime_v, dtime_a);
        return dtime;
#endif
  }
};

PS::F64 FPtreepm::eps;

class EPItreepm {
public:
  PS::S64    id;
  PS::F32    eps;
  PS::F64vec pos;

  PS::F64vec getPos() const {
    return this->pos;
  }

  void copyFromFP(const FPtreepm & fp) {
    this->id = fp.id;
    this->eps = fp.eps;
    this->pos = fp.pos;
  }
  
};

class EPJtreepm {
public:
  PS::S64    id;
  PS::F64vec pos;
  PS::F64    mass;
  //  PS::F64    rcut;

  PS::F64vec getPos() const {
    return this->pos;
  }

  PS::F64 getCharge() const {
    return this->mass;
  }

  void copyFromFP(const FPtreepm & fp) {
    this->id = fp.id;
    this->mass = fp.mass;
    this->pos = fp.pos;
  }

  PS::F64 getRSearch() const {
    PS::F64 rcut = 3.0/SIZE_OF_MESH;
    return rcut;
  }

  void setPos(const PS::F64vec pos_new) {
    this->pos = pos_new;
  }
};

inline PS::F64 gfactor_S2(const PS::F64 rad, const PS::F64 eps_pm) 
{
  PS::F64 R;
  PS::F64 g;
  PS::F64 S;

  R = 2.0*rad/eps_pm;
  R = (R > 2.0) ? 2.0 : R;
  S = R-1.0;
  S = (S > 0.0) ? S : 0.0;

  g = 1.0 + CUBE(R)*(-1.6+SQR(R)*(1.6+R*(-0.5+R*(0.15*R-12.0/35.0))))
    -CUBE(S)*CUBE(S)*(3.0/35.0+R*(18.0/35.0+0.2*R));

  return g;
}

template <class TPJ>
class calc_pp_force {
public:
  void operator () (EPItreepm *iptcl,
                    const PS::S32 ni,
                    TPJ *jptcl,
                    const PS::S32 nj,
                    Result_treepm *ppforce) 
    {
        for(PS::S32 i=0;i < ni;i++) {

            PS::F64 eps2 = SQR(iptcl[i].eps);
            for(PS::S32 j=0; j < nj;j++) {
                PS::F64vec dr = iptcl[i].pos - jptcl[j].pos;
                PS::F64 rsq = dr*dr;
                PS::F64 rad = sqrt(rsq+eps2);
                PS::F64 gfact = gfactor_S2(rad, 3.0/SIZE_OF_MESH);
                PS::F64 rinv  = 1.0/rad;
                PS::F64 mrinv3 = jptcl[j].mass*CUBE(rinv);
                ppforce[i].acc -= dr*gfact*mrinv3;
            }
        }
    }
};

void output_data(PS::ParticleSystem<FPtreepm> &ptcl, 
                 run_param &this_run, 
                 char *filename)
{
  FILE *output_fp = fopen(filename,"w");

  if(output_fp == NULL) {
    fprintf(stderr, "File %s cannot be written.",filename);
    exit(EXIT_FAILURE);
  }

  this_run.mpi_rank = PS::Comm::getRank();
  this_run.mpi_nproc = PS::Comm::getNumberOfProc();
  this_run.npart_local = ptcl.getNumberOfParticleLocal();
  this_run.npart_total = ptcl.getNumberOfParticleGlobal();
  this_run.write_header(output_fp);

  for(PS::S64 i=0;i<this_run.npart_local;i++) {
    ptcl[i].writeParticleBinary(output_fp);
  }

  fclose(output_fp);
}

void input_data(PS::ParticleSystem<FPtreepm> &ptcl, run_param &this_run, 
		char *filename)
{
  FILE *input_fp = fopen(filename,"r");

  if(input_fp == NULL) {
    fprintf(stderr, "File %s not found.\n", filename);
    exit(EXIT_FAILURE);
  }

  this_run.read_header(input_fp);

  ptcl.setNumberOfParticleLocal(this_run.npart_local);
  for(PS::S64 i=0;i<ptcl.getNumberOfParticleLocal();i++) {
    ptcl[i].readParticleBinary(input_fp);
  }

  fclose(input_fp);
}


void output_data_in_run(PS::ParticleSystem<FPtreepm> &ptcl, run_param &this_run)
{
    static char prefix[1024];
    sprintf(prefix,"snap_%05d",
            this_run.output_indx);
    if(this_run.znow < this_run.output_timing[this_run.output_indx]+0.001) {
        ptcl.writeParticleAscii(prefix, "%s_%05d_%05d.dat", this_run);
        this_run.output_indx++;
    }
/*
    static char filename[256], directory_name[256];
    sprintf(directory_name,"%s_%d",
            this_run.model_name, this_run.output_indx);
    sprintf(filename, "%s/%s_%d-%d", 
            directory_name, this_run.model_name, 
            this_run.output_indx, this_run.mpi_rank);
    
    make_directory(directory_name);
    if(this_run.znow < this_run.output_timing[this_run.output_indx]+0.001) {
        ptcl.writeParticleAscii(directory_name, "%s/snap_%05d_%05d.dat", this_run);
        this_run.output_indx++;
    }
*/
/*
    static char filename[256], directory_name[256];
    sprintf(directory_name,"%s_%d",
            this_run.model_name, this_run.output_indx);
    sprintf(filename, "%s/%s_%d-%d", 
            directory_name, this_run.model_name, 
            this_run.output_indx, this_run.mpi_rank);
    
    make_directory(directory_name);

    if(this_run.znow < this_run.output_timing[this_run.output_indx]+0.001) {
        output_data(ptcl, this_run, filename);
        this_run.output_indx++;
    }
*/
}

void drift_ptcl(PS::ParticleSystem<FPtreepm> &ptcl, PS::DomainInfo &dinfo, 
                const PS::F64 dtime)
{
    PS::S32 npart_local = ptcl.getNumberOfParticleLocal();
    for(PS::S64 i=0;i<npart_local;i++) 
        ptcl[i].pos += ptcl[i].vel*dtime;
    ptcl.adjustPositionIntoRootDomain(dinfo);
}

void reverse_ptcl_acc(PS::ParticleSystem<FPtreepm> &ptcl)
{
  PS::S64 npart_local = ptcl.getNumberOfParticleLocal();

  for(PS::S64 i=0;i<npart_local;i++) {
    ptcl[i].acc *= -1.0;
    ptcl[i].acc_pm *= -1.0;
  }

}

void kick_ptcl(PS::ParticleSystem<FPtreepm> &ptcl, const PS::F64 dtime, 
               run_param &this_run)
{
    PS::S64 npart_local = ptcl.getNumberOfParticleLocal();

    PS::F64 om = this_run.cosm.omegam;
    PS::F64 ov = this_run.cosm.omegav;

    PS::F64 anow = this_run.cosm.timetoa(this_run.tnow);
  
    //PS::F64 at = sqrt(1.e0+om*(1.e0/anow-1.e0)+ov*(SQR(anow)-1.e0))/anow; // is this right ?
    PS::F64 at = sqrt(1.e0+om*(1.e0/anow-1.e0)+ov*(1.0/(anow*anow)-1.e0))/anow;
    PS::F64 bt = 1.0/CUBE(anow);

    PS::F64 atdt1 = 1.0+at*dtime;
    PS::F64 vfact = (2.0-atdt1)/atdt1;
    PS::F64 afact = bt*dtime/atdt1;

    for(PS::S64 i=0;i<npart_local;i++) {
        ptcl[i].vel = vfact*ptcl[i].vel + afact*(ptcl[i].acc+ptcl[i].acc_pm);
        //    ptcl[i].vel = vfact*ptcl[i].vel - afact*(ptcl[i].acc+ptcl[i].acc_pm);
    }
}


PS::F64 calc_dtime(PS::ParticleSystem<FPtreepm> &ptcl, run_param &this_run){
    PS::F64 dtime;

    dtime = DBL_MAX;
    for(PS::S64 i=0;i<ptcl.getNumberOfParticleLocal();i++) {
        dtime = fmin(dtime, ptcl[i].calcDtime(this_run));
    }

    std::cerr<<"dtime="<<dtime<<std::endl;

    if(this_run.noutput >= 1) {
        COSM::REAL zred_next;
        COSM::REAL zred_next_output;
        COSM::REAL time_next_output;

        if(this_run.znow < this_run.output_timing[this_run.output_indx] && 
           this_run.output_indx+1 < this_run.noutput) {
            zred_next_output = this_run.output_timing[this_run.output_indx+1];
        }else{
            zred_next_output = this_run.output_timing[this_run.output_indx];
        }
        zred_next = this_run.cosm.timetoz(this_run.tnow+dtime);
        
        if(zred_next < zred_next_output){
            time_next_output = this_run.cosm.ztotime(zred_next_output/1.0001);
            dtime = time_next_output - this_run.tnow;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &dtime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    
    return dtime;
}


#endif /* __TREEPM_HPP__ */
