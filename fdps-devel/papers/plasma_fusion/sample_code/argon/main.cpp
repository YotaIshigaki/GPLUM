#include <particle_simulator.hpp>
using namespace PS;
using namespace std;

#include <cstdio>
#include <cassert>
#include <random>

// global constant defenition
const F64 RCUT = 4.5; // cut-off length
const F64 NUM_DENSITY = 1.05; // number density
const F64 INIT_TEMP = 1.0; // initial temperature

// particle data class defenition
struct Force{ // using struct for class without private member
  F64vec f; // force
  F64    p; // potential
  void clear(){f = 0.0; p = 0.0;}
};

struct FP{
  F64vec r,v,f; // position, velocity, force
  F64    p;     // potential
  //void clear(){r = v = f = p = 0.0;}
  void copyFromForce(const Force& fc){f = fc.f;p = fc.p;}
  F64vec getPos() const {return r;};
  void setPos(const F64vec& _r){ r = _r;};
  void kick(const F64 dt){ v += f*dt; }
  void drift(const F64 dt){r += v*dt;}
};

struct EP{
  F64vec r; // position
  F64    getRSearch() const {return RCUT;}
  F64vec getPos() const {return r;}
  void   setPos(const F64vec &_r){r = _r;};
  void   copyFromFP(const FP &fp){r = fp.r;}
};

// interaction kernel function defenition
void Kernel(const EP *epi,const S32 ni,
	    const EP *epj,const S32 nj,
	    Force *force){
  const F64 rc2 = RCUT*RCUT;
  for(S32 i=0; i<ni; i++){
    F64vec ri = epi[i].r, fi = force[i].f;
    F64 pi = force[i].p;
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
    force[i].p = 0.5*pi;
  }
}

// user-defined functions for simulation
template <typename Tps>
F64 CalcKineticEnergy(const Tps& ps){
  F64 kin = 0.0;
  for(int i=0;i<ps.getNumberOfParticleLocal();i++) kin += ps[i].v * ps[i].v;
  return 0.5*Comm::getSum(kin);
}

template <typename Tps>
void ScaleVelocity(Tps& ps,F64 temp){
  const S64 nloc = ps.getNumberOfParticleLocal();
  const S64 ntot = Comm::getSum(nloc);
  const F64 kin = CalcKineticEnergy(ps);
  const F64 sfactor = std::sqrt(3.*temp*ntot/kin);
  for(int i=0;i<nloc;i++) ps[i].v *= sfactor;
}

template <typename Tps>
F64vec CalcTotalMomentum(const Tps& ps){
  F64vec mom = 0.0;
  for(int i=0;i<ps.getNumberOfParticleLocal();i++) mom += ps[i].v;
  return Comm::getSum(mom);
}

template <typename Tps>
void RemoveTotalMomentum(Tps& ps){
  const S64 nloc = ps.getNumberOfParticleLocal();
  const S64 ntot = Comm::getSum(nloc);
  const F64vec mom = CalcTotalMomentum(ps)/ntot;
  for(int i=0;i<nloc;i++) ps[i].v -= mom;
}

template <typename Tps>
F64 AccumulatePotential(const Tps& ps){
  const S64 nloc = ps.getNumberOfParticleLocal();
  const S64 ntot = Comm::getSum(nloc);
  F64 pot = 0.0;
  for(int i=0;i<nloc;i++) pot += ps[i].p;
  return Comm::getSum(pot);
}

// main function
int main(int argc,char **argv){
  assert(argc == 2);
  Initialize(argc,argv);
  ParticleSystem<FP> ps;
  ps.initialize();
  const S32 n = atoi(argv[1]);
  const F64 l = powf(n*n*n/NUM_DENSITY,1./3.); // length of cubic system
  const F64 lh = 0.5*l;
  const F64 u = l / (F64)n; // unit lattice size
  if(Comm::getRank()==0){ // generate particles on master process
    mt19937 rnd;
    ps.setNumberOfParticleLocal(n*n*n);
    S32 count = 0;
    for(S32 x=0;x<n;x++)
      for(S32 y=0;y<n;y++)
        for(S32 z=0;z<n;z++){
          ps[count].r = u*F64vec(x,y,z)-lh;
          ps[count].v = F64vec(rnd(),rnd(),rnd());
	  count++;
        }
  }else ps.setNumberOfParticleLocal(0); // set local number of particle 0 in the other processes
  DomainInfo di;
  di.initialize(0.3);
  di.decomposeDomainAll(ps);
  ps.exchangeParticle(di);
  TreeForForceShort<Force,EP,EP>::Scatter t;
  t.initialize(3*n*n*n,0.0,64,256);
  t.calcForceAllAndWriteBack(Kernel,ps,di);
  RemoveTotalMomentum(ps);
  ScaleVelocity(ps,INIT_TEMP);

  const F64 dt = 0.0005; const F64 dth = 0.5*dt;
  S64 nl = ps.getNumberOfParticleLocal();
  for(int s=0;s<10000;s++){
    for(int i=0;i<nl;i++){
      ps[i].kick(dth);
      ps[i].drift(dt);
    }
    ps.adjustPositionIntoRootDomain(di);
    di.decomposeDomainAll(ps);
    ps.exchangeParticle(di);
    nl = ps.getNumberOfParticleLocal();
    t.calcForceAllAndWriteBack(Kernel,ps,di);
    for(int i=0;i<nl;i++) ps[i].kick(dth);
    const F64 pot = AccumulatePotential(ps);
    const F64 kin = CalcKineticEnergy(ps);
    if(Comm::getRank()==0)
      cout << scientific << pot << " " << kin << " " << pot+kin << endl;
  }
  Finalize();
}
