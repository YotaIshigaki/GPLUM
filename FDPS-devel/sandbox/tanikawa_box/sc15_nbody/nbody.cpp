#include <particle_simulator.hpp>

using namespace PS;

class Nbody{
public:
    F64    mass, eps;
    F64vec pos, vel, acc;
    F64    pot;    
    F64vec getPos() const {return pos;}
    F64 getCharge() const {return mass;}
    void copyFromFP(const Nbody &in) { 
        mass = in.mass;
        pos  = in.pos;
        eps  = in.eps;
    }
    void copyFromForce(const Nbody &out) {
        acc = out.acc;
        pot = out.pot;
    }    
    void clear() {
        acc = 0.0;
        pot = 0.0;
    }
    void readAscii(FILE *fp) {
        fscanf(fp,
               "%lf%lf%lf%lf%lf%lf%lf%lf",
               &mass, &eps,
               &pos.x, &pos.y, &pos.z,
               &vel.x, &vel.y, &vel.z);
    }
	void writeAscii(FILE* fp) const {
		fprintf(fp, "%+e %+e %+e %+e %+e %+e %+e %+e\n", 
                mass, eps,
                pos.x, pos.y, pos.z,
                vel.x, vel.y, vel.z);
	}
    void predict(F64 dt) {
        vel += (0.5 * dt) * acc;
        pos += dt * vel;
    }
    void correct(F64 dt) {
        vel += (0.5 * dt) * acc;
    }
    F64 calcEnergy() const {
        return 0.5 * mass * (vel * vel + pot + mass / eps);
    }
};

template <class TPJ>
struct CalcGrav{
    void operator () (const Nbody * ip,
                      const S32 ni,
                      const TPJ * jp,
                      const S32 nj,
                      Nbody * force) {
        for(S32 i=0; i<ni; i++){
            F64vec xi  = ip[i].pos;
            F64    ep2 = ip[i].eps
                * ip[i].eps;
            F64vec ai = 0.0;
            F64    pi = 0.0;
            for(S32 j=0; j<nj;j++){
                F64vec xj = jp[j].pos;
                F64vec dr = xi - xj;
                F64 mj  = jp[j].mass;
                F64 dr2 = dr * dr + ep2;
                F64 dri = 1.0 / sqrt(dr2);                
                pi -= dri * mj;
                ai -= (dri * dri * dri
                       * mj) * dr;
            }
            force[i].acc += ai;
            force[i].pot += pi;
        }
    }
};

template<class Tpsys>
void predict(Tpsys &p,
             const F64 dt) {
    S32 n = p.getNumberOfParticleLocal();
    for(S32 i = 0; i < n; i++)
        p[i].predict(dt);
}

template<class Tpsys>
void correct(Tpsys &p,
             const F64 dt) {
    S32 n = p.getNumberOfParticleLocal();
    for(S32 i = 0; i < n; i++)
        p[i].correct(dt);
}

template<class Tpsys>
F64 calcEnergy(const Tpsys &system) {

    F64 etot = 0.0;
    F64 etot_loc = 0.0;

    const S32 n_loc = system.getNumberOfParticleLocal();
    for(S32 i = 0; i < n_loc; i++) {
        etot_loc += system[i].calcEnergy();
    }
    etot = Comm::getSum(etot_loc);

    return etot;
}

template <class TDI, class TPS, class TTFF>
void calcGravAllAndWriteBack(TDI &dinfo,
                             TPS &ptcl,
                             TTFF &tree) {
    dinfo.decomposeDomainAll(ptcl);
    ptcl.exchangeParticle(dinfo);    
    tree.calcForceAllAndWriteBack
        (CalcGrav<Nbody>(),
         CalcGrav<SPJMonopole>(),
         ptcl, dinfo);    
}

int main(int argc, char *argv[]) {
    F32 time  = 0.0;
    const F32 tend  = 10.0;
    const F32 dtime = 1.0 / 128.0;
    const F32 dtout = 1.0 / 8.0;
    PS::Initialize(argc, argv);
    if(argc != 3) {
        fprintf(stderr, "%s <ifile> <ofile>\n", argv[0]);
        Abort(1);
    }
    PS::DomainInfo dinfo;
    dinfo.initialize();
    PS::ParticleSystem<Nbody> ptcl;
    ptcl.initialize();
    PS::TreeForForceLong <Nbody, Nbody,
        Nbody>::Monopole grav;
    grav.initialize(0);
    ptcl.readParticleAscii(argv[1]);
    calcGravAllAndWriteBack(dinfo,
                            ptcl,
                            grav);
    F64 etot0 = calcEnergy(ptcl);
    if(Comm::getRank() == 0) {
        fprintf(stderr, "time: %10.7f energy: %+e energy error: %+e\n",
                time, etot0, (etot0 - etot0) / etot0);
    }
    while(time < tend) {
        predict(ptcl, dtime);        
        calcGravAllAndWriteBack(dinfo,
                                ptcl,
                                grav);
        correct(ptcl, dtime);        
        time += dtime;
        F64 etot1 = calcEnergy(ptcl);
        if(fmod(time, dtout) == 0.0 &&
           Comm::getRank() == 0) {
                fprintf(stderr, "time: %10.7f energy: %+e energy error: %+e\n",
                        time, etot1, (etot1 - etot0) / etot0);
        }
    }
    ptcl.writeParticleAscii(argv[2]);
    PS::Finalize();
    return 0;
}
