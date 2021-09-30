#include <particle_simulator.hpp>
using namespace PS;

#include <gen_pos.hpp>
#include <string>

#include "timer.h"
Profile timer;

class FORCE{
public:
  F64vec f;
  F64    p;
  void clear(){f = 0.0; p = 0.0;}
};

class FP{
public:
  F64vec r,v,f;
  F64    p;
  F64    rs;
  void clear(){r = v = f = 0.0;}
  void copyFromForce(const FORCE& fc){f = fc.f;p = fc.p;}
  F64vec getPos()const{return r;}
  void setPos(F64vec _r){r = _r;}
};

class EP{
public:
  F64vec r;
  F64    rs;
  void   clear(){r = 0.0;}
  F64    getRSearch() const {return rs;}
  F64vec getPos() const {return r;}
  void   setPos(const F64vec &_r){r = _r;};
  void   copyFromFP(const FP &fp){r = fp.r; rs = fp.rs;}
};

struct CalcForcePot{
public:
  const F64 rc,rc2;
  CalcForcePot(const F64 _rc):rc(_rc),rc2(_rc*_rc){}
  void operator()(const EP *epi,const S32 ni,
                  const EP *epj,const S32 nj,
                  FORCE *force){
    for(S32 i=0; i<ni; i++){
      F64vec ri = epi[i].r; F64vec fi = 0.0;
      F64 pi = 0.0;
      for(S32 j=0; j<nj; j++){
        F64vec rij = ri - epj[j].r;
	const F64 r2 = rij * rij;
        if(r2==0.0 || r2>rc2) continue;
        const F64 r2i = 1.0/r2;
        const F64 r6i = r2i * r2i * r2i;
        fi += r6i*(48.0*r6i-24.0)*r2i * rij;
	pi += 4.0*r6i*(r6i-1.0);
      }
      force[i].f = fi;
      force[i].p = pi;
    }
  }
};

struct CalcForce{
public:
  const F64 rc,rc2;
  CalcForce(const F64 _rc):rc(_rc),rc2(_rc*_rc){}
  void operator()(const EP *epi,const S32 ni,
                  const EP *epj,const S32 nj,
                  FORCE *force){
    for(S32 i=0; i<ni; i++){
      F64vec ri = epi[i].r; F64vec fi = 0.0;
      for(S32 j=0; j<nj; j++){
        F64vec rij = ri - epj[j].r;
	const F64 r2 = rij * rij;
        if(r2==0.0 || r2>rc2) continue;
        const F64 r2i = 1.0/r2;
        const F64 r6i = r2i * r2i * r2i;
        fi += r6i*(48.0*r6i-24.0)*r2i * rij;
      }
      force[i].f = fi;
    }
  }
};

int main(int argc,char **argv){
  Initialize(argc,argv);

  // default vaules
  S32 nx=4,ny=4,nz=4;
  S32 nstep = 100;
  F64 rho = 1.05;
  F64 rc = 3.5;
  S32 nskip = 1;
  F64 skin = 0.0;
  // read options
  for(int i=1;i<argc;i++){
    std::string tag = argv[i];
    if(tag == "-s" || tag == "--nstep"){
      nstep = atoi(argv[++i]);
      if(Comm::getRank()==0) fprintf(stderr,"nstep = %d\n",nstep);
      continue;
    }
    if(tag == "-x"){
      nx = atoi(argv[++i]);
      if(Comm::getRank()==0) fprintf(stderr,"nx = %d\n",nx);
      continue;
    }
    if(tag == "-y"){
      ny = atoi(argv[++i]);
      if(Comm::getRank()==0) fprintf(stderr,"ny = %d\n",ny);
      continue;
    }
    if(tag == "-z"){
      nz = atoi(argv[++i]);
      if(Comm::getRank()==0) fprintf(stderr,"nz = %d\n",nz);
      continue;
    }
    if(tag == "-r" || tag == "--density"){
      rho = atof(argv[++i]);
      if(Comm::getRank()==0) fprintf(stderr,"rho = %lf\n",rho);
      continue;
    }
    if(tag == "-c" || tag == "--cutoff"){
      rc = atof(argv[++i]);
      if(Comm::getRank()==0) fprintf(stderr,"rc = %lf\n",rc);
      continue;
    }

    if(tag == "--nskip"){
      nskip = atoi(argv[++i]);
      if(Comm::getRank()==0) fprintf(stderr,"nskip = %d\n",nskip);
      continue;
    }
    if(tag == "--skin"){
      skin = atof(argv[++i]);
      if(Comm::getRank()==0) fprintf(stderr,"skin = %lf\n",skin);
      continue;
    }
    if(Comm::getRank()==0) fprintf(stderr,"error: undefined option %s\n",tag.c_str());
    PS::Abort();
  }
  const PS::F64 rs = rc + skin;

  ParticleSystem<FP> ps;
  ps.initialize();
  DomainInfo di;
  di.initialize();
  generateFCC(ps,di,nx,ny,nz,rho);
  for(int i=0;i<ps.getNumberOfParticleLocal();i++) ps[i].rs = rs;
  //for(int i=0;i<ps.getNumberOfParticleLocal();i++) printf("%d %d %lf %lf %lf\n",i,Comm::getRank(),ps[i].r.x,ps[i].r.y,ps[i].r.z);

  di.decomposeDomainAll(ps);
  ps.exchangeParticle(di);
  TreeForForceShort<FORCE,EP,EP>::Scatter t;
  t.initialize(3*4*nx*ny*nz,0.0,64,256);
  t.calcForceAllAndWriteBack(CalcForce(rc),ps,di);

  const F64 dt = 0.0001;
  const F64 dth = 0.5*dt;
  S32 nl = ps.getNumberOfParticleLocal();
  INTERACTION_LIST_MODE isReuse = MAKE_LIST;

  timer.beg(Profile::TOTAL);
  for(int s=0;s<nstep;s++){
    timer.beg(Profile::INTEG);
    for(int i=0;i<nl;i++) ps[i].v+=ps[i].f*dth;
    for(int i=0;i<nl;i++){
      ps[i].r += ps[i].v*dt;
    }
    timer.end(Profile::INTEG);
    if(s%nskip == 0){
      ps.adjustPositionIntoRootDomain(di);
      timer.beg(Profile::DECOMP);
      di.decomposeDomainAll(ps);
      timer.end(Profile::DECOMP);
      timer.beg(Profile::EXCHANGE);
      ps.exchangeParticle(di);
      timer.end(Profile::EXCHANGE);
      nl = ps.getNumberOfParticleLocal();
      if(nskip == 0) isReuse = MAKE_LIST;
      else           isReuse = MAKE_LIST_FOR_REUSE;
    }else{
      isReuse = REUSE_LIST;
    }
    timer.beg(Profile::FORCE);
    t.calcForceAllAndWriteBack(
			       #ifdef CALC_ENG
			       CalcForcePot(rc),
			       #else
			       CalcForce(rc),
			       #endif
			       ps,di,true,
			       isReuse
			       );
    timer.end(Profile::FORCE);
    timer.beg(Profile::INTEG);
    for(int i=0;i<nl;i++) ps[i].v+=ps[i].f*dth;
    timer.end(Profile::INTEG);
    timer.beg(Profile::OTHER);
#ifdef CALC_ENG
    F64 pot = 0.0, kin = 0.0;
    for(int i=0;i<nl;i++){
      pot += ps[i].p;
      kin += 0.5*(ps[i].v*ps[i].v);
    }
    pot = Comm::getSum(pot);
    kin = Comm::getSum(kin);
#endif
    timer.end(Profile::OTHER);
    //if(Comm::getRank()==0) printf("%e %e %e\n",pot,kin,pot+kin);
  }
  timer.end(Profile::TOTAL);
  Finalize();
}
