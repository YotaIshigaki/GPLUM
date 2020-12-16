#include <cmath>
#include "PairListInteraction.h"
#include "ShortRangeInteraction.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef CHECK_ENERGY
double shortpair_energy;
double ljpair_energy;
#endif

// pair list not holds value of pos, charge, LJ
// LJ force
void pairlistloopljf1(const int *jindex,
		      const int *lj,
		      const int npair,
		      const double posi[3],
		      const double chargei,
		      const double (*ljmp)[4],
		      const PosChargeArray& particlej,
		      const double cutoff2,
		      double force[3])
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  //  double e=energy;
#ifdef SIMPLE_CUTOFF
#else
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
#endif

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
#endif

#ifdef __FCC_VERSION
#pragma loop novrec 
#endif
  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

#ifdef K_SIMD
    dx = posi[0]-pacs[jid][0];
    dy = posi[1]-pacs[jid][1];
    dz = posi[2]-pacs[jid][2];
#else
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
#endif

    double r2 = dx*dx+dy*dy+dz*dz;

#ifdef SIMPLE_CUTOFF
    double inside;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
      inside=1.0;
    }
#else
    if(r2>cutoff2)r2=cutoff2;
#endif

    double _r = 1.0/sqrt(r2);

#ifdef SIMPLE_CUTOFF
    double _r2 = _r*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
    //    e += (ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6 )*inside;
    double dp = (ljmp[atj][2]*_r14 - ljmp[atj][3]*_r8 )*inside;
#else
    double r = r2*_r;
    double r3 = r2*r;
    double _r2 = _r*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
    //    e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
    //         - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
    //           );
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
                 );
#endif
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  //  energy = e;
}
// pair list not holds value of pos, charge, LJ
// LJ force and energy
void pairlistloopljfe1(const int *jindex,
                        const int *lj,
                        const int npair,
                        const double posi[3],
                        const double chargei,
                        const double (*ljmp)[4],
                        const PosChargeArray& particlej,
                        const double cutoff2,
                        double force[3],
                        double& energy)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
#ifdef SIMPLE_CUTOFF
#else
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
#endif

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
#endif

#ifdef __FCC_VERSION
#pragma loop novrec 
#endif
  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

#ifdef K_SIMD
    dx = posi[0]-pacs[jid][0];
    dy = posi[1]-pacs[jid][1];
    dz = posi[2]-pacs[jid][2];
#else
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
#endif

    double r2 = dx*dx+dy*dy+dz*dz;

#ifdef SIMPLE_CUTOFF
    double inside;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
      inside=1.0;
    }
#else
    if(r2>cutoff2)r2=cutoff2;
#endif

    double _r = 1.0/sqrt(r2);

#ifdef SIMPLE_CUTOFF
    double _r2 = _r*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
# ifdef CHECK_ENERGY
    double ljsingle = (ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6 )*inside;
    ljpair_energy += ljsingle;
    e += ljsingle;
# else
    e += (ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6 )*inside;
# endif
    double dp = (ljmp[atj][2]*_r14 - ljmp[atj][3]*_r8 )*inside;
#else
    double r = r2*_r;
    double r3 = r2*r;
    double _r2 = _r*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
    ljpair_energy += ljsingle;
    e += ljsingle;
# else
    e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
         - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
           );
# endif
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
                 );
#endif
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
}
// pair list not holds value of pos, charge, LJ
// LJ force, energy and virial
void pairlistloopljfe1(const int *jindex,
                        const int *lj,
                        const int npair,
                        const double posi[3],
                        const double chargei,
                        const double (*ljmp)[4],
                        const PosChargeArray& particlej,
                        const double cutoff2,
                        double force[3],
		       double& energy,
		       double& virial)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  double v=virial;
#ifdef SIMPLE_CUTOFF
#else
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
#endif

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
#endif

#ifdef __FCC_VERSION
#pragma loop novrec 
#endif
  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

#ifdef K_SIMD
    dx = posi[0]-pacs[jid][0];
    dy = posi[1]-pacs[jid][1];
    dz = posi[2]-pacs[jid][2];
#else
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
#endif

    double r2 = dx*dx+dy*dy+dz*dz;

#ifdef SIMPLE_CUTOFF
    double inside;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
      inside=1.0;
    }
#else
    if(r2>cutoff2)r2=cutoff2;
#endif

    double _r = 1.0/sqrt(r2);

#ifdef SIMPLE_CUTOFF
    double _r2 = _r*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
# ifdef CHECK_ENERGY
    double ljsingle = (ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6 )*inside;
    ljpair_energy += ljsingle;
    e += ljsingle;
# else
    e += (ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6 )*inside;
# endif
    double dp = (ljmp[atj][2]*_r14 - ljmp[atj][3]*_r8 )*inside;
    //    v += (ljmp[atj][2]*_r12 - ljmp[atj][3]*_r6 )*inside;
    v += dp*r2;
#else
    double r = r2*_r;
    double r3 = r2*r;
    double _r2 = _r*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
# ifdef CHECK_ENERGY
    double ljsingle = ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
			- ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
			);
    ljpair_energy += ljsingle;
    e += ljsingle;
# else
    e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
         - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
           );
# endif
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
                 );
    /*
    v += (ljmp[atj][2]*(_r12-r3*(r14co1-r14co2*r))
	  - ljmp[atj][3]*(_r6 -r3*(r8co1 -r8co2 *r))
	  );
    */
    v += dp*r2;
#endif
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
  virial = v;
}

// pair list not holds value of pos, charge, LJ
// Coulomb and LJ force
void pairlistloopljcf1(const int *jindex,
		       const int *lj,
		       const int npair,
		       const double posi[3],
		       const double chargei,
		       const double (*ljmp)[4],
		       const PosChargeArray& particlej,
		       const double cutoff2,
		       double force[3])
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
#ifdef SIMPLE_CUTOFF
#else
  //  double &r1co0 = ShortRange::r1co0, &r1co1=ShortRange::r1co1, &r1co2=ShortRange::r1co2;
  //  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  //  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r3co1=ShortRange::r3co1, &r3co2=ShortRange::r3co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
#endif

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
#endif

#ifdef __FCC_VERSION
#pragma loop novrec
#endif
  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

#ifdef K_SIMD
    dx = posi[0]-pacs[jid][0];
    dy = posi[1]-pacs[jid][1];
    dz = posi[2]-pacs[jid][2];
#else
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
#endif

    double r2 = dx*dx+dy*dy+dz*dz;
#ifdef K_SIMD
    double cij = chargei*pacs[jid][3];
#else
    double cij = chargei*particlej[jid].charge;
#endif

#ifdef SIMPLE_CUTOFF
    double inside;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
      inside=1.0;
    }
#else
    if(r2>cutoff2)r2=cutoff2;
#endif

    double _r = 1.0/sqrt(r2);

#ifdef SIMPLE_CUTOFF
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    //    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
    double dp = (ljmp[atj][2]*_r14 - ljmp[atj][3]*_r8 + cij*_r3)*inside;
#else
    double r = r2*_r;
    //    double r3 = r2*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    //    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3 -r*(r3co1 -r3co2 *r))
                 );
#endif
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
}
// pair list not holds value of pos, charge, LJ
// Coulomb and LJ force and energy
void pairlistloopljcfe1(const int *jindex,
                        const int *lj,
                        const int npair,
                        const double posi[3],
                        const double chargei,
                        const double (*ljmp)[4],
                        const PosChargeArray& particlej,
                        const double cutoff2,
                        double force[3],
                        double& energy)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
#ifdef SIMPLE_CUTOFF
#else
  double &r1co0 = ShortRange::r1co0, &r1co1=ShortRange::r1co1, &r1co2=ShortRange::r1co2;
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r3co1=ShortRange::r3co1, &r3co2=ShortRange::r3co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
#endif

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
#endif

#ifdef __FCC_VERSION
#pragma loop novrec
#endif
  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

#ifdef K_SIMD
    dx = posi[0]-pacs[jid][0];
    dy = posi[1]-pacs[jid][1];
    dz = posi[2]-pacs[jid][2];
#else
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
#endif

    double r2 = dx*dx+dy*dy+dz*dz;
#ifdef K_SIMD
    double cij = chargei*pacs[jid][3];
#else
    double cij = chargei*particlej[jid].charge;
#endif

#ifdef SIMPLE_CUTOFF
    double inside;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
      inside=1.0;
    }
#else
    if(r2>cutoff2)r2=cutoff2;
#endif

    double _r = 1.0/sqrt(r2);

#ifdef SIMPLE_CUTOFF
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6;
    double esingle = cij     *_r;
    ljpair_energy += ljsingle*inside;
    shortpair_energy += esingle*inside;
    e += (ljsingle + esingle)*inside;
# else
    e += (ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6 + cij*_r)*inside;
# endif
    double dp = (ljmp[atj][2]*_r14 - ljmp[atj][3]*_r8 + cij*_r3)*inside;
#else
    double r = r2*_r;
    double r3 = r2*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
    double esingle = cij     *(_r  +(r3*(r1co1 -r1co2 *r)-r1co0)) ;
    ljpair_energy += ljsingle;
    shortpair_energy += esingle;
    e += ( ljsingle + esingle );
# else
    e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
         - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
         + cij     *(_r  +(r3*(r1co1 -r1co2 *r)-r1co0 ))
           );
# endif
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3 -r*(r3co1 -r3co2 *r))
                 );
#endif
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
}
// pair list not holds value of pos, charge, LJ
// Coulomb and LJ force, energy and virial
void pairlistloopljcfe1(const int *jindex,
                        const int *lj,
                        const int npair,
                        const double posi[3],
                        const double chargei,
                        const double (*ljmp)[4],
                        const PosChargeArray& particlej,
                        const double cutoff2,
                        double force[3],
                        double& energy,
			double& virial)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  double v=virial;
#ifdef SIMPLE_CUTOFF
#else
  double &r1co0 = ShortRange::r1co0, &r1co1=ShortRange::r1co1, &r1co2=ShortRange::r1co2;
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r3co1=ShortRange::r3co1, &r3co2=ShortRange::r3co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
#endif

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
#endif

#ifdef __FCC_VERSION
#pragma loop novrec
#endif
  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

#ifdef K_SIMD
    dx = posi[0]-pacs[jid][0];
    dy = posi[1]-pacs[jid][1];
    dz = posi[2]-pacs[jid][2];
#else
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
#endif

    double r2 = dx*dx+dy*dy+dz*dz;
#ifdef K_SIMD
    double cij = chargei*pacs[jid][3];
#else
    double cij = chargei*particlej[jid].charge;
#endif

#ifdef SIMPLE_CUTOFF
    double inside;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
      inside=1.0;
    }
#else
    if(r2>cutoff2)r2=cutoff2;
#endif

    double _r = 1.0/sqrt(r2);

#ifdef SIMPLE_CUTOFF
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6;
    double esingle = cij     *_r;
    ljpair_energy += ljsingle*inside;
    shortpair_energy += esingle*inside;
    e += (ljsingle + esingle)*inside;
# else
    e += (ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6 + cij*_r)*inside;
# endif
    double dp = (ljmp[atj][2]*_r14 - ljmp[atj][3]*_r8 + cij*_r3)*inside;
    //    v += (ljmp[atj][2]*_r12 - ljmp[atj][3]*_r6 + cij*_r)*inside;
    v += r2*dp;
#else
    double r = r2*_r;
    double r3 = r2*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
    double esingle = cij     *(_r  +(r3*(r1co1 -r1co2 *r)-r1co0)) ;
    ljpair_energy += ljsingle;
    shortpair_energy += esingle;
    e += ( ljsingle + esingle );
# else
    e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
         - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
         + cij     *(_r  +(r3*(r1co1 -r1co2 *r)-r1co0 ))
           );
# endif
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3 -r*(r3co1 -r3co2 *r))
                 );
    /*
    v += (ljmp[atj][2]*(_r12-r3*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r6 -r3*(r8co1 -r8co2 *r))
               + cij     *(_r -r3*(r3co1 -r3co2 *r))
                 );
    */
    v += r2*dp;
#endif
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
  virial = v;
}

void pairlistloopljcfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const PosChargeArray& poschargei,
                       const AtomtypeArray& atomia,
                       const PosChargeArray& poschargej,
                       const int npl,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations)
{
  int i;
  double e=energy;
  double v=virial;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {poschargei[iid].position.x,
                      poschargei[iid].position.y,
                      poschargei[iid].position.z};
    double chargei = poschargei[iid].charge;
    int atomi;
    atomi = atomia[iid];
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doVirialcalculation){
      pairlistloopljcfe1(jindex[i],lj[i],npair[i],
			 posi, chargei, 
			 ljmp,
			 poschargej, 
			 cutoff2,
			 force[i],
			 e,
			 v);
    }else{
#ifdef CPPMD_ENABLE_FORTRAN_KERNEL
      pairlistloopljcfe1_i_posc_f_(jindex[i],lj[i],&npair[i],
				   posi, &chargei, 
				   ljmp,
				   &poschargej[0].position.x, 
				   &cutoff2,
				   force[i],
				   &e);
#else
      if(operations.doEnergycalculation){
	pairlistloopljcfe1(jindex[i],lj[i],npair[i],
			   posi, chargei, 
			   ljmp,
			   poschargej, 
			   cutoff2,
			   force[i],
			   e);
      }else{
	pairlistloopljcf1(jindex[i],lj[i],npair[i],
			  posi, chargei, 
			  ljmp,
			  poschargej, 
			  cutoff2,
			  force[i]);
      }
#endif
    }
#ifdef OVERLAP
# ifdef OVERLAP_PAIRLIST_THREAD
    {
      int tn = 0;
#ifdef _OPENMP
      tn = omp_get_thread_num();
#endif
      if(tn==0){
	MPICHECKER::mpi_checker();
      }
    }
# endif
#endif
  }
  energy = e;
}

// with offset(selfnpair) pairlist
void pairlistloopljcfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const PosChargeArray& poschargei,
                       const AtomtypeArray& atomia,
                       const PosChargeArray& poschargej,
                       const int npl,
                       const int *selfnpair,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations)
{
  int i;
  double e=energy;
  double v=virial;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {poschargei[iid].position.x,
                      poschargei[iid].position.y,
                      poschargei[iid].position.z};
    double chargei = poschargei[iid].charge;
    int atomi;
    atomi = atomia[iid];
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doVirialcalculation){
      pairlistloopljcfe1(&(jindex[i][selfnpair[i]]),
			 &(lj[i][selfnpair[i]]),npair[i],
			 posi, chargei, 
			 ljmp,
			 poschargej, 
			 cutoff2,
			 force[i],
			 e,
			 v);
    }else{
#ifdef CPPMD_ENABLE_FORTRAN_KERNEL
      pairlistloopljcfe1_i_posc_f_(&(jindex[i][selfnpair[i]]),
				   &(lj[i][selfnpair[i]]),&npair[i],
				   posi, &chargei, 
				   ljmp,
				   &poschargej[0].position.x, 
				   &cutoff2,
				   force[i],
				   &e);
#else
      if(operations.doEnergycalculation){
	pairlistloopljcfe1(&(jindex[i][selfnpair[i]]),
			   &(lj[i][selfnpair[i]]),npair[i],
			   posi, chargei, 
			   ljmp,
			   poschargej, 
			   cutoff2,
			   force[i],
			   e);
      }else{
	pairlistloopljcf1(&(jindex[i][selfnpair[i]]),
			  &(lj[i][selfnpair[i]]),npair[i],
			  posi, chargei, 
			  ljmp,
			  poschargej, 
			  cutoff2,
			  force[i]);
      }
#endif
    }
  }
  energy = e;
  virial = v;
}

// no offset pairlist
// LJ
template<typename PA, typename GPA>
void pairlistloopljfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const PA& particlei,
                       const GPA& particlej,
                       const int npl,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations)
{
  int i;
  double e=energy;
  double v=virial;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {particlei.poscharge[iid].position.x,
                      particlei.poscharge[iid].position.y,
                      particlei.poscharge[iid].position.z};
    double chargei = particlei.poscharge[iid].charge;
    int atomi = particlei.atomtype[iid];
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doVirialcalculation){
      pairlistloopljfe1(jindex[i],lj[i],npair[i],
			posi, chargei, 
			ljmp,
			particlej.poscharge, 
			cutoff2,
			force[i],
			e,
			v);
    }else{
#ifdef CPPMD_ENABLE_FORTRAN_KERNEL
      pairlistloopljfe1_i_posc_f_(jindex[i],lj[i],&npair[i],
				  posi, &chargei, 
				  ljmp,
				  &particlej.poscharge[0].position.x, 
				  &cutoff2,
				  force[i],
				  &e);
#else
      if(operations.doEnergycalculation){
	pairlistloopljfe1(jindex[i],lj[i],npair[i],
			  posi, chargei, 
			  ljmp,
			  particlej.poscharge, 
			  cutoff2,
			  force[i],
			  e);
      }else{
	pairlistloopljf1(jindex[i],lj[i],npair[i],
			 posi, chargei, 
			 ljmp,
			 particlej.poscharge, 
			 cutoff2,
			 force[i]);
      }
#endif
    }
#ifdef OVERLAP
# ifdef OVERLAP_PAIRLIST_THREAD
    {
      int tn = 0;
#ifdef _OPENMP
      tn = omp_get_thread_num();
#endif
      if(tn==0){
	MPICHECKER::mpi_checker();
      }
    }
# endif
#endif
  }
  energy = e;
  virial = v;
}
// no offset pairlist
// LJ Coulomb
template<typename PA, typename GPA>
void pairlistloopljcfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const PA& particlei,
                       const GPA& particlej,
                       const int npl,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations)
{
  int i;
  double e=energy;
  double v=virial;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {particlei.poscharge[iid].position.x,
                      particlei.poscharge[iid].position.y,
                      particlei.poscharge[iid].position.z};
    double chargei = particlei.poscharge[iid].charge;
    int atomi = particlei.atomtype[iid];
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doVirialcalculation){
      pairlistloopljcfe1(jindex[i],lj[i],npair[i],
			 posi, chargei, 
			 ljmp,
			 particlej.poscharge, 
			 cutoff2,
			 force[i],
			 e,
			 v);
    }else{
#ifdef CPPMD_ENABLE_FORTRAN_KERNEL
      pairlistloopljcfe1_i_posc_f_(jindex[i],lj[i],&npair[i],
				   posi, &chargei, 
				   ljmp,
				   &particlej.poscharge[0].position.x, 
				   &cutoff2,
				   force[i],
				   &e);
#else
      if(operations.doEnergycalculation){
	pairlistloopljcfe1(jindex[i],lj[i],npair[i],
			   posi, chargei, 
			   ljmp,
			   particlej.poscharge, 
			   cutoff2,
			   force[i],
			   e);
      }else{
	pairlistloopljcf1(jindex[i],lj[i],npair[i],
			  posi, chargei, 
			  ljmp,
			  particlej.poscharge, 
			  cutoff2,
			  force[i]);
      }
#endif
    }
#ifdef OVERLAP
# ifdef OVERLAP_PAIRLIST_THREAD
    {
      int tn = 0;
#ifdef _OPENMP
      tn = omp_get_thread_num();
#endif
      if(tn==0){
	MPICHECKER::mpi_checker();
      }
    }
# endif
#endif
  }
  energy = e;
  virial = v;
}

// with offset(selfnpair) pairlist
// jindex and lj are shared self and ghost
// LJ
template<typename PA, typename GPA>
void pairlistloopljfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const PA& particlei,
                       const GPA& particlej,
                       const int npl,
                       const int *selfnpair,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations)
{
  int i;
  double e=energy;
  double v=virial;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {particlei.poscharge[iid].position.x,
                      particlei.poscharge[iid].position.y,
                      particlei.poscharge[iid].position.z};
    double chargei = particlei.poscharge[iid].charge;
    int atomi = particlei.atomtype[iid];
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doVirialcalculation){
      pairlistloopljfe1(&(jindex[i][selfnpair[i]]),
			&(lj[i][selfnpair[i]]),npair[i],
			posi, chargei, 
			ljmp,
			particlej.poscharge, 
			cutoff2,
			force[i],
			e,
			v);
    }else{
#ifdef CPPMD_ENABLE_FORTRAN_KERNEL
      pairlistloopljfe1_i_posc_f_(&(jindex[i][selfnpair[i]]),
				  &(lj[i][selfnpair[i]]),&npair[i],
				  posi, &chargei, 
				  ljmp,
				  &particlej.poscharge[0].position.x, 
				  &cutoff2,
				  force[i],
				  &e);
#else
      if(operations.doEnergycalculation){
	pairlistloopljfe1(&(jindex[i][selfnpair[i]]),
			  &(lj[i][selfnpair[i]]),npair[i],
			  posi, chargei, 
			  ljmp,
			  particlej.poscharge, 
			  cutoff2,
			  force[i],
			  e);
      }else{
	pairlistloopljf1(&(jindex[i][selfnpair[i]]),
			 &(lj[i][selfnpair[i]]),npair[i],
			 posi, chargei, 
			 ljmp,
			 particlej.poscharge, 
			 cutoff2,
			 force[i]);
      }
#endif
    }
  }
  energy = e;
  virial = v;
}
// with offset(selfnpair) pairlist
// jindex and lj are shared self and ghost
// LJ Coulomb
template<typename PA, typename GPA>
void pairlistloopljcfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const PA& particlei,
                       const GPA& particlej,
                       const int npl,
                       const int *selfnpair,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations)
{
  int i;
  double e=energy;
  double v=virial;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {particlei.poscharge[iid].position.x,
                      particlei.poscharge[iid].position.y,
                      particlei.poscharge[iid].position.z};
    double chargei = particlei.poscharge[iid].charge;
    int atomi = particlei.atomtype[iid];
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doVirialcalculation){
      pairlistloopljcfe1(&(jindex[i][selfnpair[i]]),
			 &(lj[i][selfnpair[i]]),npair[i],
			 posi, chargei, 
			 ljmp,
			 particlej.poscharge, 
			 cutoff2,
			 force[i],
			 e,
			 v);
    }else{
#ifdef CPPMD_ENABLE_FORTRAN_KERNEL
      pairlistloopljcfe1_i_posc_f_(&(jindex[i][selfnpair[i]]),
				   &(lj[i][selfnpair[i]]),&npair[i],
				   posi, &chargei, 
				   ljmp,
				   &particlej.poscharge[0].position.x, 
				   &cutoff2,
				   force[i],
				   &e);
#else
      if(operations.doEnergycalculation){
	pairlistloopljcfe1(&(jindex[i][selfnpair[i]]),
			   &(lj[i][selfnpair[i]]),npair[i],
			   posi, chargei, 
			   ljmp,
			   particlej.poscharge, 
			   cutoff2,
			   force[i],
			   e);
      }else{
	pairlistloopljcf1(&(jindex[i][selfnpair[i]]),
			  &(lj[i][selfnpair[i]]),npair[i],
			  posi, chargei, 
			  ljmp,
			  particlej.poscharge, 
			  cutoff2,
			  force[i]);
      }
#endif
    }
  }
  energy = e;
  virial = v;
}

template
void pairlistloopljfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const CombinedParticleArray& particlei,
                       const CombinedParticleArray& particlej,
                       const int npl,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations);
template
void pairlistloopljfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const CombinedParticleArray& particlei,
                       const GhostParticleArray& particlej,
                       const int npl,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations);
template
void pairlistloopljfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const CombinedParticleArray& particlei,
                       const GhostParticleArray& particlej,
                       const int npl,
                       const int *selfnpair,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations);

template
void pairlistloopljcfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const CombinedParticleArray& particlei,
                       const CombinedParticleArray& particlej,
                       const int npl,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations);
template
void pairlistloopljcfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const CombinedParticleArray& particlei,
                       const GhostParticleArray& particlej,
                       const int npl,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations);
template
void pairlistloopljcfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const CombinedParticleArray& particlei,
                       const GhostParticleArray& particlej,
                       const int npl,
                       const int *selfnpair,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations);


void pairlistloopljcfe1_wolf_ljshift(const double alpha,
                                      const int *jindex,
                                      const int *lj,
                                      const int npair,
                                      const double posi[3],
                                      const double chargei,
                                      const double (*ljmp)[4],
                                      const PosChargeArray& particlej,
                                      const double cutoff2,
                                      double force[3],
				     double& energy,
				     double& virial)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  double v=virial;
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
  double ee, edp;
  const double alpharc2 = alpha*alpha*cutoff2;
  const double alpharc = alpha*sqrt(cutoff2);
  const double ec2 = M_2_SQRTPI*exp(-alpharc2);
  const double ec1 = 2.0*erfc(alpharc)/sqrt(cutoff2);
  const double ec = -ec1 - ec2;
  const double el = erfc(alpharc)/cutoff2+ec2/sqrt(cutoff2);

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
#endif

#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
  const double (*ftc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_force_table[0]));
  const double (*ptc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_pot_table[0]));
# endif
#endif

#ifdef __FCC_VERSION
#pragma loop novrec 
#endif
  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

#ifdef K_SIMD
    dx = posi[0]-pacs[jid][0];
    dy = posi[1]-pacs[jid][1];
    dz = posi[2]-pacs[jid][2];
#else
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
#endif

    double r2 = dx*dx+dy*dy+dz*dz;
#ifdef K_SIMD
    double cij = chargei*pacs[jid][3];
#else
    double cij = chargei*particlej[jid].charge;
#endif

    double inside = 1.0;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
    }

    double _r = 1.0/sqrt(r2);

    double r = r2*_r;
    double ar = alpha*r;
    double r3 = r2*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
    int iar = int((ar)*4.0);
    double dar = ((ar)*4.0-iar)*2.0-1.0;
    ee = ptc[iar][0]+dar*(ptc[iar][1]+dar*(ptc[iar][2]+dar*(ptc[iar][3]+dar*ptc[iar][4])));
    edp = ftc[iar][0]+dar*(ftc[iar][1]+dar*(ftc[iar][2]+dar*(ftc[iar][3]+dar*ftc[iar][4])));
# else
    ewaldrealforceandpottbl_(&ee,&edp,&ar);
# endif
#else
    ee = erfc(ar);
    edp = M_2_SQRTPI*ar*exp(-ar*ar)+ee;
#endif
# ifdef CHECK_ENERGY
    double ljsingle =ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
    double esingle = cij     *((_r*ee) + el*r + ec);
    ljpair_energy += ljsingle*inside;
    shortpair_energy += esingle*inside;
    e += ( ljsingle + esingle )*inside;
# else
    e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
         - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
	   + cij     *((_r*ee) + el*r + ec)
           )*inside;
# endif
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3*edp - el)
                 )*inside;
    v += dp*r2;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
  virial = v;
}
void pairlistloopljcfe1_wolf_ljshift(const double alpha,
                                      const int *jindex,
                                      const int *lj,
                                      const int npair,
                                      const double posi[3],
                                      const double chargei,
                                      const double (*ljmp)[4],
                                      const PosChargeArray& particlej,
                                      const double cutoff2,
                                      double force[3],
                                      double& energy)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
  double ee, edp;
  const double alpharc2 = alpha*alpha*cutoff2;
  const double alpharc = alpha*sqrt(cutoff2);
  const double ec2 = M_2_SQRTPI*exp(-alpharc2);
  const double ec1 = 2.0*erfc(alpharc)/sqrt(cutoff2);
  const double ec = -ec1 - ec2;
  const double el = erfc(alpharc)/cutoff2+ec2/sqrt(cutoff2);

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
#endif

#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
  const double (*ftc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_force_table[0]));
  const double (*ptc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_pot_table[0]));
# endif
#endif

#ifdef __FCC_VERSION
#pragma loop novrec 
#endif
  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

#ifdef K_SIMD
    dx = posi[0]-pacs[jid][0];
    dy = posi[1]-pacs[jid][1];
    dz = posi[2]-pacs[jid][2];
#else
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
#endif

    double r2 = dx*dx+dy*dy+dz*dz;
#ifdef K_SIMD
    double cij = chargei*pacs[jid][3];
#else
    double cij = chargei*particlej[jid].charge;
#endif

    double inside = 1.0;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
    }

    double _r = 1.0/sqrt(r2);

    double r = r2*_r;
    double ar = alpha*r;
    double r3 = r2*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
    int iar = int((ar)*4.0);
    double dar = ((ar)*4.0-iar)*2.0-1.0;
    ee = ptc[iar][0]+dar*(ptc[iar][1]+dar*(ptc[iar][2]+dar*(ptc[iar][3]+dar*ptc[iar][4])));
    edp = ftc[iar][0]+dar*(ftc[iar][1]+dar*(ftc[iar][2]+dar*(ftc[iar][3]+dar*ftc[iar][4])));
# else
    ewaldrealforceandpottbl_(&ee,&edp,&ar);
# endif
#else
    ee = erfc(ar);
    edp = M_2_SQRTPI*ar*exp(-ar*ar)+ee;
#endif
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
    double esingle = cij     *((_r*ee) + el*r + ec);
    ljpair_energy += ljsingle*inside;
    shortpair_energy += esingle*inside;
    e += (ljsingle + esingle)*inside;
# else
    e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
         - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
	   + cij     *((_r*ee) + el*r + ec)
           )*inside;
# endif
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3*edp - el)
                 )*inside;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
}
void pairlistloopljcf1_wolf_ljshift(const double alpha,
                                      const int *jindex,
                                      const int *lj,
                                      const int npair,
                                      const double posi[3],
                                      const double chargei,
                                      const double (*ljmp)[4],
                                      const PosChargeArray& particlej,
                                      const double cutoff2,
                                      double force[3])
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  //  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  //  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
  double ee, edp;
  const double alpharc2 = alpha*alpha*cutoff2;
  const double alpharc = alpha*sqrt(cutoff2);
  const double ec2 = M_2_SQRTPI*exp(-alpharc2);
  const double ec1 = 2.0*erfc(alpharc)/sqrt(cutoff2);
  //  const double ec = -ec1 - ec2;
  const double el = erfc(alpharc)/cutoff2+ec2/sqrt(cutoff2);

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
#endif

#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
  const double (*ftc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_force_table[0]));
  const double (*ptc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_pot_table[0]));
# endif
#endif

#ifdef __FCC_VERSION
#pragma loop novrec 
#endif
  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

#ifdef K_SIMD
    dx = posi[0]-pacs[jid][0];
    dy = posi[1]-pacs[jid][1];
    dz = posi[2]-pacs[jid][2];
#else
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
#endif

    double r2 = dx*dx+dy*dy+dz*dz;
#ifdef K_SIMD
    double cij = chargei*pacs[jid][3];
#else
    double cij = chargei*particlej[jid].charge;
#endif

    double inside = 1.0;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
    }

    double _r = 1.0/sqrt(r2);

    double r = r2*_r;
    double ar = alpha*r;
    //    double r3 = r2*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    //    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
    int iar = int((ar)*4.0);
    double dar = ((ar)*4.0-iar)*2.0-1.0;
    ee = ptc[iar][0]+dar*(ptc[iar][1]+dar*(ptc[iar][2]+dar*(ptc[iar][3]+dar*ptc[iar][4])));
    edp = ftc[iar][0]+dar*(ftc[iar][1]+dar*(ftc[iar][2]+dar*(ftc[iar][3]+dar*ftc[iar][4])));
# else
    ewaldrealforceandpottbl_(&ee,&edp,&ar);
# endif
#else
    ee = erfc(ar);
    edp = M_2_SQRTPI*ar*exp(-ar*ar)+ee;
#endif
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3*edp - el)
                 )*inside;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
}

// no offset pairlist
// LJ Wolf
template<typename PA, typename GPA>
void pairlistloopljcfe_wolf_ljshift(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const PA& particlei,
                       const GPA& particlej,
                       const int npl,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations)
{
  int i;
  double e=energy;
  double v=virial;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {particlei.poscharge[iid].position.x,
                      particlei.poscharge[iid].position.y,
                      particlei.poscharge[iid].position.z};
    double chargei = particlei.poscharge[iid].charge;
    int atomi = particlei.atomtype[iid];
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doVirialcalculation){
      pairlistloopljcfe1_wolf_ljshift(ShortRange::alpha,
				      jindex[i],lj[i],npair[i],
				      posi, chargei, 
				      ljmp,
				      particlej.poscharge, 
				      cutoff2,
				      force[i],
				      e,
				      v);
    }else{
      if(operations.doEnergycalculation){
	pairlistloopljcfe1_wolf_ljshift(ShortRange::alpha,
					jindex[i],lj[i],npair[i],
					posi, chargei, 
					ljmp,
					particlej.poscharge, 
					cutoff2,
					force[i],
					e);
      }else{
	pairlistloopljcf1_wolf_ljshift(ShortRange::alpha,
				       jindex[i],lj[i],npair[i],
				       posi, chargei, 
				       ljmp,
				       particlej.poscharge, 
				       cutoff2,
				       force[i]);
      }
    }
#ifdef OVERLAP
# ifdef OVERLAP_PAIRLIST_THREAD
    {
      int tn = 0;
#ifdef _OPENMP
      tn = omp_get_thread_num();
#endif
      if(tn==0){
	MPICHECKER::mpi_checker();
      }
    }
# endif
#endif
  }
  energy = e;
  virial = v;
}
// with offset(selfnpair) pairlist
// jindex and lj are shared self and ghost
// LJ Wolf
template<typename PA, typename GPA>
void pairlistloopljcfe_wolf_ljshift(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const PA& particlei,
                       const GPA& particlej,
                       const int npl,
				    const int *selfnpair,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations)
{
  int i;
  double e=energy;
  double v=virial;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {particlei.poscharge[iid].position.x,
                      particlei.poscharge[iid].position.y,
                      particlei.poscharge[iid].position.z};
    double chargei = particlei.poscharge[iid].charge;
    int atomi = particlei.atomtype[iid];
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doVirialcalculation){
      pairlistloopljcfe1_wolf_ljshift(ShortRange::alpha,
				      &(jindex[i][selfnpair[i]]),
				      &(lj[i][selfnpair[i]]),
				      npair[i],
				      posi, chargei, 
				      ljmp,
				      particlej.poscharge, 
				      cutoff2,
				      force[i],
				      e,
				      v);
    }else{
      if(operations.doEnergycalculation){
	pairlistloopljcfe1_wolf_ljshift(ShortRange::alpha,
					&(jindex[i][selfnpair[i]]),
					&(lj[i][selfnpair[i]]),
					npair[i],
					posi, chargei, 
					ljmp,
					particlej.poscharge, 
					cutoff2,
					force[i],
					e);
      }else{
	pairlistloopljcf1_wolf_ljshift(ShortRange::alpha,
				       &(jindex[i][selfnpair[i]]),
				       &(lj[i][selfnpair[i]]),
				       npair[i],
				       posi, chargei, 
				       ljmp,
				       particlej.poscharge, 
				       cutoff2,
				       force[i]);
      }
    }
  }
  energy = e;
  virial = v;
}

template
void pairlistloopljcfe_wolf_ljshift(const int (*jindex)[MAX_PAIR],
				    const int (*lj)[MAX_PAIR],
				    const int *npair,
				    const int *iindex,
				    const CombinedParticleArray& particlei,
				    const CombinedParticleArray& particlej,
				    const int npl,
				    const double cutoff2,
				    double (*force)[3],
				    double& energy,
				    double& virial,
				    const OperationSelector operations);
template
void pairlistloopljcfe_wolf_ljshift(const int (*jindex)[MAX_PAIR],
				    const int (*lj)[MAX_PAIR],
				    const int *npair,
				    const int *iindex,
				    const CombinedParticleArray& particlei,
				    const GhostParticleArray& particlej,
				    const int npl,
				    const double cutoff2,
				    double (*force)[3],
				    double& energy,
				    double& virial,
				    const OperationSelector operations);
template
void pairlistloopljcfe_wolf_ljshift(const int (*jindex)[MAX_PAIR],
				    const int (*lj)[MAX_PAIR],
				    const int *npair,
				    const int *iindex,
				    const CombinedParticleArray& particlei,
				    const GhostParticleArray& particlej,
				    const int npl,
				    const int *selfnpair,
				    const double cutoff2,
				    double (*force)[3],
				    double& energy,
				    double& virial,
				    const OperationSelector operations);

void pairlistloopljcfe1(const int *jindex,
                        const int *lj,
                        const int npair,
                        const double posi[3],
                        const double chargei,
                        const double (*ljmp)[4],
                        const ParticleArray& particlej,
                        const double cutoff2,
                        double force[3],
                        double& energy)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
#ifdef SIMPLE_CUTOFF
#else
  double &r1co0 = ShortRange::r1co0, &r1co1=ShortRange::r1co1, &r1co2=ShortRange::r1co2;
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r3co1=ShortRange::r3co1, &r3co2=ShortRange::r3co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
#endif

  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
    double r2 = dx*dx+dy*dy+dz*dz;
    double cij = chargei*particlej[jid].charge;
#ifdef SIMPLE_CUTOFF
    double inside = 1.0;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
      inside=1.0;
    }
#else
    if(r2>cutoff2)r2=cutoff2;
#endif

    double _r = 1.0/sqrt(r2);

#ifdef SIMPLE_CUTOFF
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6;
    double esingle = cij     *_r;
    ljpair_energy += ljsingle*inside;
    shortpair_energy += esingle*inside;
    e += (ljsingle + esingle)*inside;
# else
    e += (ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6 + cij*_r)*inside;
# endif
    double dp = (ljmp[atj][2]*_r14 - ljmp[atj][3]*_r8 + cij*_r3)*inside;
#else
    double r = r2*_r;
    double r3 = r2*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
    double esingle = cij     *(_r  +(r3*(r1co1 -r1co2 *r)-r1co0)) ;
    ljpair_energy += ljsingle;
    shortpair_energy += esingle;
    e += ( ljsingle + esingle );
# else
    e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
         - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
         + cij     *(_r  +(r3*(r1co1 -r1co2 *r)-r1co0 ))
           );
# endif
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3 -r*(r3co1 -r3co2 *r))
                 );
#endif
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
}

void pairlistloopljcfe1(const int *jindex,
                        const int *lj,
                        const int npair,
                        const double posi[3],
                        const double chargei,
                        const double (*ljmp)[4],
                        const ParticleArray& particlej,
                        const double cutoff2,
                        double force[3],
                        double& energy,
			double& virial)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  double v=virial;
#ifdef SIMPLE_CUTOFF
#else
  double &r1co0 = ShortRange::r1co0, &r1co1=ShortRange::r1co1, &r1co2=ShortRange::r1co2;
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r3co1=ShortRange::r3co1, &r3co2=ShortRange::r3co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
#endif

  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
    double r2 = dx*dx+dy*dy+dz*dz;
    double cij = chargei*particlej[jid].charge;
#ifdef SIMPLE_CUTOFF
    double inside = 1.0;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
      inside=1.0;
    }
#else
    if(r2>cutoff2)r2=cutoff2;
#endif

    double _r = 1.0/sqrt(r2);

#ifdef SIMPLE_CUTOFF
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6;
    double esingle = cij     *_r;
    ljpair_energy += ljsingle*inside;
    shortpair_energy += esingle*inside;
    e += (ljsingle + esingle)*inside;
# else
    e += (ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6 + cij*_r)*inside;
# endif
    double dp = (ljmp[atj][2]*_r14 - ljmp[atj][3]*_r8 + cij*_r3)*inside;
    //    v += (ljmp[atj][2]*_r12 - ljmp[atj][3]*_r6 + cij*_r)*inside;
    v += dp*r2;
#else
    double r = r2*_r;
    double r3 = r2*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
    double esingle = cij     *(_r  +(r3*(r1co1 -r1co2 *r)-r1co0)) ;
    ljpair_energy += ljsingle;
    shortpair_energy += esingle;
    e += ( ljsingle + esingle );
# else
    e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
         - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
         + cij     *(_r  +(r3*(r1co1 -r1co2 *r)-r1co0 ))
           );
# endif
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3 -r*(r3co1 -r3co2 *r))
                 );
    /*
    v += (ljmp[atj][2]*(_r12-r3*(r14co1-r14co2*r))
	  - ljmp[atj][3]*(_r6 -r3*(r8co1 -r8co2 *r))
	  + cij     *(_r -r3*(r3co1 -r3co2 *r))
	  );
    */
    v += dp*r2;
#endif
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
  virial = v;
}

template<>
void pairlistloopljcfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const ParticleArray& particlei,
                       const ParticleArray& particlej,
                       const int npl,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations)
{
  int i;
  double e=energy;
  double v=virial;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {particlei[iid].position.x,
                      particlei[iid].position.y,
                      particlei[iid].position.z};
    double chargei = particlei[iid].charge;
    int atomi;
    atomi = particlei[iid].atomtype;
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doVirialcalculation){
      pairlistloopljcfe1(jindex[i],lj[i],npair[i],
			 posi, chargei, 
			 ljmp,
			 particlej, cutoff2,
			 force[i],
			 e,
			 v);
    }else{
      pairlistloopljcfe1(jindex[i],lj[i],npair[i],
			 posi, chargei, 
			 ljmp,
			 particlej, cutoff2,
			 force[i],
			 e);
    }
#ifdef OVERLAP
# ifdef OVERLAP_PAIRLIST_THREAD
    {
      int tn = 0;
#ifdef _OPENMP
      tn = omp_get_thread_num();
#endif
      if(tn==0){
	MPICHECKER::mpi_checker();
      }
    }
# endif
#endif
  }
  energy = e;
  virial = v;
}

template<>
void pairlistloopljcfe(const int (*jindex)[MAX_PAIR],
                       const int (*lj)[MAX_PAIR],
                       const int *npair,
                       const int *iindex,
                       const ParticleArray& particlei,
                       const ParticleArray& particlej,
                       const int npl,
		       const int *selfnpair,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy,
		       double& virial,
		       const OperationSelector operations)
{
  int i;
  double e=energy;
  double v=virial;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {particlei[iid].position.x,
                      particlei[iid].position.y,
                      particlei[iid].position.z};
    double chargei = particlei[iid].charge;
    int atomi;
    atomi = particlei[iid].atomtype;
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doVirialcalculation){
      pairlistloopljcfe1(&(jindex[i][selfnpair[i]]),
			 &(lj[i][selfnpair[i]]),npair[i],
			 posi, chargei, 
			 ljmp,
			 particlej, cutoff2,
			 force[i],
			 e,
			 v);
    }else{
      pairlistloopljcfe1(&(jindex[i][selfnpair[i]]),
			 &(lj[i][selfnpair[i]]),npair[i],
			 posi, chargei, 
			 ljmp,
			 particlej, cutoff2,
			 force[i],
			 e);
    }
  }
  energy = e;
  virial = v;
}

void pairlistloopljcf1_ewald_ljshift(const double alpha,
				     const int *jindex,
				     const int *lj,
				     const int npair,
				     const double posi[3],
				     const double chargei,
				     const double (*ljmp)[4],
				     const PosChargeArray& particlej,
				     const double cutoff2,
				     double force[3] )
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  //  double e=energy;
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
  double ee, edp;

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
#endif

#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
  const double (*ftc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_force_table[0]));
  const double (*ptc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_pot_table[0]));
# endif
#endif

#ifdef __FCC_VERSION
#pragma loop novrec 
#endif
  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

#ifdef K_SIMD
    dx = posi[0]-pacs[jid][0];
    dy = posi[1]-pacs[jid][1];
    dz = posi[2]-pacs[jid][2];
#else
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
#endif

    double r2 = dx*dx+dy*dy+dz*dz;
#ifdef K_SIMD
    double cij = chargei*pacs[jid][3];
#else
    double cij = chargei*particlej[jid].charge;
#endif

    double inside = 1.0;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
    }

    double _r = 1.0/sqrt(r2);

    double r = r2*_r;
    double ar = alpha*r;
    double r3 = r2*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
    int iar = int((ar)*4.0);
    double dar = ((ar)*4.0-iar)*2.0-1.0;
    //    ee = ptc[iar][0]+dar*(ptc[iar][1]+dar*(ptc[iar][2]+dar*(ptc[iar][3]+dar*ptc[iar][4])));
    edp = ftc[iar][0]+dar*(ftc[iar][1]+dar*(ftc[iar][2]+dar*(ftc[iar][3]+dar*ftc[iar][4])));
# else
    ewaldrealforceandpottbl_(&ee,&edp,&ar);
# endif
#else
    ee = erfc(ar);
    edp = M_2_SQRTPI*ar*exp(-ar*ar)+ee;
#endif
    //    e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
    //         - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
    //         + cij     *(_r*ee)
    //           )*inside;
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3*edp)
                 )*inside;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  //  energy = e;
}
void pairlistloopljcfe1_ewald_ljshift(const double alpha,
                                      const int *jindex,
                                      const int *lj,
                                      const int npair,
                                      const double posi[3],
                                      const double chargei,
                                      const double (*ljmp)[4],
                                      const PosChargeArray& particlej,
                                      const double cutoff2,
                                      double force[3],
                                      double& energy)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
  double ee, edp;

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
#endif

#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
  const double (*ftc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_force_table[0]));
  const double (*ptc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_pot_table[0]));
# endif
#endif

#ifdef __FCC_VERSION
#pragma loop novrec 
#endif
  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

#ifdef K_SIMD
    dx = posi[0]-pacs[jid][0];
    dy = posi[1]-pacs[jid][1];
    dz = posi[2]-pacs[jid][2];
#else
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
#endif

    double r2 = dx*dx+dy*dy+dz*dz;
#ifdef K_SIMD
    double cij = chargei*pacs[jid][3];
#else
    double cij = chargei*particlej[jid].charge;
#endif

    double inside = 1.0;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
    }

    double _r = 1.0/sqrt(r2);

    double r = r2*_r;
    double ar = alpha*r;
    double r3 = r2*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
    int iar = int((ar)*4.0);
    double dar = ((ar)*4.0-iar)*2.0-1.0;
    ee = ptc[iar][0]+dar*(ptc[iar][1]+dar*(ptc[iar][2]+dar*(ptc[iar][3]+dar*ptc[iar][4])));
    edp = ftc[iar][0]+dar*(ftc[iar][1]+dar*(ftc[iar][2]+dar*(ftc[iar][3]+dar*ftc[iar][4])));
# else
    ewaldrealforceandpottbl_(&ee,&edp,&ar);
# endif
#else
    ee = erfc(ar);
    edp = M_2_SQRTPI*ar*exp(-ar*ar)+ee;
#endif
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
    double esingle = cij     *(_r*ee);
    ljpair_energy += ljsingle*inside;
    shortpair_energy += esingle*inside;
    e += ( ljsingle + esingle )*inside;
# else
    e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
         - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
         + cij     *(_r*ee)
           )*inside;
# endif
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3*edp)
                 )*inside;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
}

template<typename PA, typename GPA>
void pairlistloopljcfe_ewald_ljshift(const int (*jindex)[MAX_PAIR],
                                     const int (*lj)[MAX_PAIR],
                                     const int *npair,
                                     const int *iindex,
                                     const PA& particlei,
                                     const GPA& particlej,
                                     const int npl,
                                     const double cutoff2,
                                     double (*force)[3],
                                     double& energy,
				     double& virial,
				     const OperationSelector operations)
{
  int i;
  double e=energy;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {particlei.poscharge[iid].position.x,
                      particlei.poscharge[iid].position.y,
                      particlei.poscharge[iid].position.z};
    double chargei = particlei.poscharge[iid].charge;
    int atomi;
    atomi = particlei.atomtype[iid];
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doEnergycalculation){
      pairlistloopljcfe1_ewald_ljshift(ShortRange::alpha,
				       jindex[i],
				       lj[i],npair[i],
				       posi, chargei, 
				       ljmp,
				       particlej.poscharge, 
				       cutoff2,
				       force[i],
				       e);
    }else{
      pairlistloopljcf1_ewald_ljshift(ShortRange::alpha,
				      jindex[i],
				      lj[i],npair[i],
				      posi, chargei, 
				      ljmp,
				      particlej.poscharge, 
				      cutoff2,
				      force[i] );
    }
#ifdef OVERLAP
# ifdef OVERLAP_PAIRLIST_THREAD
    {
      int tn = 0;
#ifdef _OPENMP
      tn = omp_get_thread_num();
#endif
      if(tn==0){
	MPICHECKER::mpi_checker();
      }
    }
# endif
#endif
  }
  energy = e;
}

template<typename PA, typename GPA>
void pairlistloopljcfe_ewald_ljshift(const int (*jindex)[MAX_PAIR],
                                     const int (*lj)[MAX_PAIR],
                                     const int *npair,
                                     const int *iindex,
                                     const PA& particlei,
                                     const GPA& particlej,
                                     const int npl,
                                     const int *selfnpair,
                                     const double cutoff2,
                                     double (*force)[3],
                                     double& energy,
				     double& virial,
				     const OperationSelector operations)
{
  int i;
  double e=energy;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {particlei.poscharge[iid].position.x,
                      particlei.poscharge[iid].position.y,
                      particlei.poscharge[iid].position.z};
    double chargei =  particlei.poscharge[iid].charge;
    int atomi;
    atomi = particlei.atomtype[iid];
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doEnergycalculation){
      pairlistloopljcfe1_ewald_ljshift(ShortRange::alpha,
				       &(jindex[i][selfnpair[i]]),
				       &(lj[i][selfnpair[i]]),npair[i],
				       posi, chargei, 
				       ljmp,
				       particlej.poscharge,
				       cutoff2,
				       force[i],
				       e);
    }else{
      pairlistloopljcf1_ewald_ljshift(ShortRange::alpha,
				      &(jindex[i][selfnpair[i]]),
				      &(lj[i][selfnpair[i]]),npair[i],
				      posi, chargei, 
				      ljmp,
				      particlej.poscharge,
				      cutoff2,
				      force[i]);
    }
  }
  energy = e;
}

void pairlistloopljcf1_ewald(const double alpha,
			     const int *jindex,
			     const int *lj,
			     const int npair,
			     const double posi[3],
			     const double chargei,
			     const double (*ljmp)[4],
			     const PosChargeArray& particlej,
			     const double cutoff2,
			     double force[3])
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  //  double e=energy;
  double ee, edp;

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
#endif

#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
  const double (*ftc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_force_table[0]));
  //  const double (*ptc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_pot_table[0]));
# endif
#endif

  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

#ifdef K_SIMD
    dx = posi[0]-pacs[jid][0];
    dy = posi[1]-pacs[jid][1];
    dz = posi[2]-pacs[jid][2];
#else
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
#endif

    double r2 = dx*dx+dy*dy+dz*dz;
#ifdef K_SIMD
    double cij = chargei*pacs[jid][3];
#else
    double cij = chargei*particlej[jid].charge;
#endif

    double inside = 1.0;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
    }

    double _r = 1.0/sqrt(r2);

    double r = r2*_r;
    double ar = alpha*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
#ifdef K_SIMD
    int iar = int((ar)*4.0);
    double dar = ((ar)*4.0-iar)*2.0-1.0;
    //    ee = ptc[iar][0]+dar*(ptc[iar][1]+dar*(ptc[iar][2]+dar*(ptc[iar][3]+dar*ptc[iar][4])));
    edp = ftc[iar][0]+dar*(ftc[iar][1]+dar*(ftc[iar][2]+dar*(ftc[iar][3]+dar*ftc[iar][4])));
#else
    ewaldrealforceandpottbl_(&ee,&edp,&ar);
#endif
#else
    ee = erfc(ar);
    edp = M_2_SQRTPI*ar*exp(-ar*ar)+ee;
#endif
    //    e += ( ljmp[atj][0]*(_r12)
    //         - ljmp[atj][1]*(_r6)
    //         + cij     *(_r*ee)
    //           )*inside;
    double dp = (ljmp[atj][2]*(_r14)
               - ljmp[atj][3]*(_r8)
               + cij     *(_r3*edp)
                 )*inside;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  //  energy = e;
}

#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
void eeedp(double *pot, double *dp, const double *x)
{
  
  int ix = int((*x)*4.0);
  const double *c=ewaldRealInterpolation.ewald_real_force_table+ix*5;
  const double *pc=ewaldRealInterpolation.ewald_real_pot_table+ix*5;

  // Note: the following line may be unnecessary for normal situations.
  //  if ((x<0.0)||(x>=16.0)) return ret;

  double dx = ((*x)*4.0-ix)*2.0-1.0;
  *dp = c[0]+dx*(c[1]+dx*(c[2]+dx*(c[3]+dx*c[4])));
  *pot = pc[0]+dx*(pc[1]+dx*(pc[2]+dx*(pc[3]+dx*pc[4])));
};
#endif

void pairlistloopljcfe1_ewald(const double alpha,
                              const int *jindex,
                              const int *lj,
                              const int npair,
                              const double posi[3],
                              const double chargei,
                              const double (*ljmp)[4],
                              const PosChargeArray& particlej,
                              const double cutoff2,
                              double force[3],
                              double& energy)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  double ee, edp;

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
#endif

#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
  const double (*ftc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_force_table[0]));
  const double (*ptc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_pot_table[0]));
# else
#  define DIRECT_TABLE
#  ifdef DIRECT_TABLE
  const double (*ftc) = (&(ewaldRealInterpolation.ewald_real_force_table[0]));
  const double (*ptc) = (&(ewaldRealInterpolation.ewald_real_pot_table[0]));
#  endif
# endif
#endif

  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

#ifdef K_SIMD
    dx = posi[0]-pacs[jid][0];
    dy = posi[1]-pacs[jid][1];
    dz = posi[2]-pacs[jid][2];
#else
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
#endif

    double r2 = dx*dx+dy*dy+dz*dz;
#ifdef K_SIMD
    double cij = chargei*pacs[jid][3];
#else
    double cij = chargei*particlej[jid].charge;
#endif

    double inside = 1.0;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
    }

    double _r = 1.0/sqrt(r2);

    double r = r2*_r;
    double ar = alpha*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
    int iar = int((ar)*4.0);
    double dar = ((ar)*4.0-iar)*2.0-1.0;
    ee = ptc[iar][0]+dar*(ptc[iar][1]+dar*(ptc[iar][2]+dar*(ptc[iar][3]+dar*ptc[iar][4])));
    edp = ftc[iar][0]+dar*(ftc[iar][1]+dar*(ftc[iar][2]+dar*(ftc[iar][3]+dar*ftc[iar][4])));
# else
#  ifdef DIRECT_TABLE
#   if 0
    int iar = int((ar)*4.0);
    double dar = ((ar)*4.0-iar)*2.0-1.0;
    iar *= 5;
    ee = ptc[iar+0]+dar*(ptc[iar+1]+dar*(ptc[iar+2]+dar*(ptc[iar+3]+dar*ptc[iar+4])));
    edp = ftc[iar+0]+dar*(ftc[iar+1]+dar*(ftc[iar+2]+dar*(ftc[iar+3]+dar*ftc[iar+4])));
#   else
#    if 0
    int ix = int((ar)*4.0);
    const double *c=ewaldRealInterpolation.ewald_real_force_table+ix*5;
    const double *pc=ewaldRealInterpolation.ewald_real_pot_table+ix*5;
    double dar = ((ar)*4.0-ix)*2.0-1.0;
    edp = c[0]+dar*(c[1]+dar*(c[2]+dar*(c[3]+dar*c[4])));
    ee = pc[0]+dar*(pc[1]+dar*(pc[2]+dar*(pc[3]+dar*pc[4])));
#    else
    eeedp(&ee,&edp,&ar);
#    endif
#   endif
#  else
    ewaldrealforceandpottbl_(&ee,&edp,&ar);
#  endif
# endif
#else
    ee = erfc(ar);
    edp = M_2_SQRTPI*ar*exp(-ar*ar)+ee;
#endif
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6;
    double esingle = cij     *(_r*ee);
    ljpair_energy += ljsingle*inside;
    shortpair_energy += esingle*inside;
    e += ( ljsingle + esingle )*inside;
# else
    e += ( ljmp[atj][0]*(_r12)
         - ljmp[atj][1]*(_r6)
         + cij     *(_r*ee)
           )*inside;
# endif
    double dp = (ljmp[atj][2]*(_r14)
               - ljmp[atj][3]*(_r8)
               + cij     *(_r3*edp)
                 )*inside;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
}
void pairlistloopljcfe1_ewald(const double alpha,
                              const int *jindex,
                              const int *lj,
                              const int npair,
                              const double posi[3],
                              const double chargei,
                              const double (*ljmp)[4],
                              const PosChargeArray& particlej,
                              const double cutoff2,
                              double force[3],
                              double& energy,
			      double& virial)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  double ee, edp;
  double v=virial;

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
#endif

#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
  const double (*ftc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_force_table[0]));
  const double (*ptc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_pot_table[0]));
# endif
#endif

  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

#ifdef K_SIMD
    dx = posi[0]-pacs[jid][0];
    dy = posi[1]-pacs[jid][1];
    dz = posi[2]-pacs[jid][2];
#else
    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;
#endif

    double r2 = dx*dx+dy*dy+dz*dz;
#ifdef K_SIMD
    double cij = chargei*pacs[jid][3];
#else
    double cij = chargei*particlej[jid].charge;
#endif

    double inside = 1.0;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
    }

    double _r = 1.0/sqrt(r2);

    double r = r2*_r;
    double ar = alpha*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
#ifdef K_SIMD
    int iar = int((ar)*4.0);
    double dar = ((ar)*4.0-iar)*2.0-1.0;
    ee = ptc[iar][0]+dar*(ptc[iar][1]+dar*(ptc[iar][2]+dar*(ptc[iar][3]+dar*ptc[iar][4])));
    edp = ftc[iar][0]+dar*(ftc[iar][1]+dar*(ftc[iar][2]+dar*(ftc[iar][3]+dar*ftc[iar][4])));
#else
    ewaldrealforceandpottbl_(&ee,&edp,&ar);
#endif
#else
    ee = erfc(ar);
    edp = M_2_SQRTPI*ar*exp(-ar*ar)+ee;
#endif
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6;
    double esingle = cij     *(_r*ee);
    ljpair_energy += ljsingle*inside;
    shortpair_energy += esingle*inside;
    e += ( ljsingle + esingle )*inside;
# else
    e += ( ljmp[atj][0]*(_r12)
         - ljmp[atj][1]*(_r6)
         + cij     *(_r*ee)
           )*inside;
# endif
    double dp = (ljmp[atj][2]*(_r14)
               - ljmp[atj][3]*(_r8)
               + cij     *(_r3*edp)
                 )*inside;
    v += dp*r2;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
  virial = v;
}

template<typename PA, typename GPA>
void pairlistloopljcfe_ewald(const int (*jindex)[MAX_PAIR],
                             const int (*lj)[MAX_PAIR],
                             const int *npair,
                             const int *iindex,
                             const PA& particlei,
                             const GPA& particlej,
                             const int npl,
                             const double cutoff2,
                             double (*force)[3],
                             double& energy,
			     double& virial,
			     const OperationSelector operations)
{
  int i;
  double e=energy;
  double v=virial;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {particlei.poscharge[iid].position.x,
                      particlei.poscharge[iid].position.y,
                      particlei.poscharge[iid].position.z};
    double chargei = particlei.poscharge[iid].charge;
    int atomi;
    atomi = particlei.atomtype[iid];
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doVirialcalculation){
      pairlistloopljcfe1_ewald(ShortRange::alpha,
			       jindex[i],
			       lj[i],npair[i],
			       posi, chargei, 
			       ljmp,
			       particlej.poscharge, 
			       cutoff2,
			       force[i],
			       e,
			       v);
    }else{
#ifdef CPPMD_ENABLE_FORTRAN_KERNEL
      pairlistloopljewaldfe1_i_posc_f_(&ShortRange::alpha,
				       jindex[i],
				       lj[i],&npair[i],
				       posi, &chargei, 
				       ljmp,
				       &particlej.poscharge[0].position.x, 
				       &cutoff2,
				       force[i],
				       &e);
#else
      if(operations.doEnergycalculation){
	pairlistloopljcfe1_ewald(ShortRange::alpha,
				 jindex[i],
				 lj[i],npair[i],
				 posi, chargei, 
				 ljmp,
				 particlej.poscharge, 
				 cutoff2,
				 force[i],
				 e);
      }else{
	pairlistloopljcf1_ewald(ShortRange::alpha,
				jindex[i],
				lj[i],npair[i],
				posi, chargei, 
				ljmp,
				particlej.poscharge, 
				cutoff2,
				force[i] );
      }
#endif
    }
#ifdef OVERLAP
# ifdef OVERLAP_PAIRLIST_THREAD
    {
      int tn = 0;
#ifdef _OPENMP
      tn = omp_get_thread_num();
#endif
      if(tn==0){
	MPICHECKER::mpi_checker();
      }
    }
# endif
#endif
  }
  energy = e;
  virial = v;
}

template<typename PA, typename GPA>
void pairlistloopljcfe_ewald(const int (*jindex)[MAX_PAIR],
                             const int (*lj)[MAX_PAIR],
                             const int *npair,
                             const int *iindex,
                             const PA& particlei,
                             const GPA& particlej,
                             const int npl,
                             const int *selfnpair,
                             const double cutoff2,
                             double (*force)[3],
                             double& energy,
			     double& virial,
			     const OperationSelector operations)
{
  int i;
  double e=energy;
  double v=virial;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {particlei.poscharge[iid].position.x,
                      particlei.poscharge[iid].position.y,
                      particlei.poscharge[iid].position.z};
    double chargei = particlei.poscharge[iid].charge;
    int atomi;
    atomi = particlei.atomtype[iid];
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doVirialcalculation){
      pairlistloopljcfe1_ewald(ShortRange::alpha,
			       &(jindex[i][selfnpair[i]]),
			       &(lj[i][selfnpair[i]]),npair[i],
			       posi, chargei, 
			       ljmp,
			       particlej.poscharge, 
			       cutoff2,
			       force[i],
			       e,
			       v);
    }else{
#ifdef CPPMD_ENABLE_FORTRAN_KERNEL
      pairlistloopljewaldfe1_i_posc_f_(&ShortRange::alpha,
				       &(jindex[i][selfnpair[i]]),
				       &(lj[i][selfnpair[i]]),&npair[i],
				       posi, &chargei, 
				       ljmp,
				       &particlej.poscharge[0].position.x, 
				       &cutoff2,
				       force[i],
				       &e);
#else
      if(operations.doEnergycalculation){
	pairlistloopljcfe1_ewald(ShortRange::alpha,
				 &(jindex[i][selfnpair[i]]),
				 &(lj[i][selfnpair[i]]),npair[i],
				 posi, chargei, 
				 ljmp,
				 particlej.poscharge, 
				 cutoff2,
				 force[i],
				 e);
      }else{
	pairlistloopljcf1_ewald(ShortRange::alpha,
				&(jindex[i][selfnpair[i]]),
				&(lj[i][selfnpair[i]]),npair[i],
				posi, chargei, 
				ljmp,
				particlej.poscharge, 
				cutoff2,
				force[i] );
      }
#endif
    }
  }
  energy = e;
  virial = v;
}

template
void pairlistloopljcfe_ewald_ljshift(const int (*jindex)[MAX_PAIR],
                                     const int (*lj)[MAX_PAIR],
                                     const int *npair,
                                     const int *iindex,
                                     const CombinedParticleArray& particlei,
                                     const CombinedParticleArray& particlej,
                                     const int npl,
                                     const double cutoff2,
                                     double (*force)[3],
                                     double& energy,
				     double& virial,
				     const OperationSelector operations);
template
void pairlistloopljcfe_ewald_ljshift(const int (*jindex)[MAX_PAIR],
                                     const int (*lj)[MAX_PAIR],
                                     const int *npair,
                                     const int *iindex,
                                     const CombinedParticleArray& particlei,
                                     const GhostParticleArray& particlej,
                                     const int npl,
                                     const double cutoff2,
                                     double (*force)[3],
                                     double& energy,
				     double& virial,
				     const OperationSelector operations);
template
void pairlistloopljcfe_ewald_ljshift(const int (*jindex)[MAX_PAIR],
                                     const int (*lj)[MAX_PAIR],
                                     const int *npair,
                                     const int *iindex,
                                     const CombinedParticleArray& particlei,
                                     const GhostParticleArray& particlej,
                                     const int npl,
                                     const int *selfnpair,
                                     const double cutoff2,
                                     double (*force)[3],
                                     double& energy,
				     double& virial,
				     const OperationSelector operations);
template
void pairlistloopljcfe_ewald(const int (*jindex)[MAX_PAIR],
                             const int (*lj)[MAX_PAIR],
                             const int *npair,
                             const int *iindex,
                             const CombinedParticleArray& particlei,
                             const CombinedParticleArray& particlej,
                             const int npl,
                             const double cutoff2,
                             double (*force)[3],
                             double& energy,
			     double& virial,
		       const OperationSelector operations);
template
void pairlistloopljcfe_ewald(const int (*jindex)[MAX_PAIR],
                             const int (*lj)[MAX_PAIR],
                             const int *npair,
                             const int *iindex,
                             const CombinedParticleArray& particlei,
                             const GhostParticleArray& particlej,
                             const int npl,
                             const double cutoff2,
                             double (*force)[3],
                             double& energy,
			     double& virial,
		       const OperationSelector operations);
template
void pairlistloopljcfe_ewald(const int (*jindex)[MAX_PAIR],
                             const int (*lj)[MAX_PAIR],
                             const int *npair,
                             const int *iindex,
                             const CombinedParticleArray& particlei,
                             const GhostParticleArray& particlej,
                             const int npl,
                             const int *selfnpair,
                             const double cutoff2,
                             double (*force)[3],
                             double& energy,
			     double& virial,
		       const OperationSelector operations);

void pairlistloopljcf1_ewald(const double alpha,
			     const int *jindex,
			     const int *lj,
			     const int npair,
			     const double posi[3],
			     const double chargei,
			     const double (*ljmp)[4],
			     const ParticleArray& particlej,
			     const double cutoff2,
			     double force[3])
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  //  double e=energy;
  double ee, edp;

#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
  const double (*ftc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_force_table[0]));
  //  const double (*ptc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_pot_table[0]));
# endif
#endif

  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;

    double r2 = dx*dx+dy*dy+dz*dz;
    double cij = chargei*particlej[jid].charge;

    double inside = 1.0;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
    }

    double _r = 1.0/sqrt(r2);

    double r = r2*_r;
    double ar = alpha*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
    int iar = int((ar)*4.0);
    double dar = ((ar)*4.0-iar)*2.0-1.0;
    //    ee = ptc[iar][0]+dar*(ptc[iar][1]+dar*(ptc[iar][2]+dar*(ptc[iar][3]+dar*ptc[iar][4])));
    edp = ftc[iar][0]+dar*(ftc[iar][1]+dar*(ftc[iar][2]+dar*(ftc[iar][3]+dar*ftc[iar][4])));
# else
    ewaldrealforceandpottbl_(&ee,&edp,&ar);
# endif
#else
    ee = erfc(ar);
    edp = M_2_SQRTPI*ar*exp(-ar*ar)+ee;
#endif
    //    e += ( ljmp[atj][0]*(_r12)
    //         - ljmp[atj][1]*(_r6)
    //         + cij     *(_r*ee)
    //           )*inside;
    double dp = (ljmp[atj][2]*(_r14)
               - ljmp[atj][3]*(_r8)
               + cij     *(_r3*edp)
                 )*inside;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  //  energy = e;
}
void pairlistloopljcfe1_ewald(const double alpha,
                              const int *jindex,
                              const int *lj,
                              const int npair,
                              const double posi[3],
                              const double chargei,
                              const double (*ljmp)[4],
                              const ParticleArray& particlej,
                              const double cutoff2,
                              double force[3],
                              double& energy)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  double ee, edp;

#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
  const double (*ftc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_force_table[0]));
  const double (*ptc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_pot_table[0]));
# endif
#endif

  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;

    double r2 = dx*dx+dy*dy+dz*dz;
    double cij = chargei*particlej[jid].charge;

    double inside = 1.0;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
    }

    double _r = 1.0/sqrt(r2);

    double r = r2*_r;
    double ar = alpha*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
    int iar = int((ar)*4.0);
    double dar = ((ar)*4.0-iar)*2.0-1.0;
    ee = ptc[iar][0]+dar*(ptc[iar][1]+dar*(ptc[iar][2]+dar*(ptc[iar][3]+dar*ptc[iar][4])));
    edp = ftc[iar][0]+dar*(ftc[iar][1]+dar*(ftc[iar][2]+dar*(ftc[iar][3]+dar*ftc[iar][4])));
# else
    ewaldrealforceandpottbl_(&ee,&edp,&ar);
 # endif
#else
    ee = erfc(ar);
    edp = M_2_SQRTPI*ar*exp(-ar*ar)+ee;
#endif
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6;
    double esingle = cij     *(_r*ee);
    ljpair_energy += ljsingle*inside;
    shortpair_energy += esingle*inside;
    e += ( ljsingle + esingle )*inside;
# else
    e += ( ljmp[atj][0]*(_r12)
         - ljmp[atj][1]*(_r6)
         + cij     *(_r*ee)
           )*inside;
# endif
    double dp = (ljmp[atj][2]*(_r14)
               - ljmp[atj][3]*(_r8)
               + cij     *(_r3*edp)
                 )*inside;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
}
void pairlistloopljcfe1_ewald(const double alpha,
                              const int *jindex,
                              const int *lj,
                              const int npair,
                              const double posi[3],
                              const double chargei,
                              const double (*ljmp)[4],
                              const ParticleArray& particlej,
                              const double cutoff2,
                              double force[3],
                              double& energy,
			      double& virial)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  double ee, edp;
  double v=virial;

#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
  const double (*ftc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_force_table[0]));
  const double (*ptc)[5] = (double (*)[5])(&(ewaldRealInterpolation.ewald_real_pot_table[0]));
# endif
#endif

  for(j=0;j<npair;j++){
    int jid = jindex[j];
    int atj = lj[j];

    double dx, dy, dz;

    dx = posi[0]-particlej[jid].position.x;
    dy = posi[1]-particlej[jid].position.y;
    dz = posi[2]-particlej[jid].position.z;

    double r2 = dx*dx+dy*dy+dz*dz;
    double cij = chargei*particlej[jid].charge;

    double inside = 1.0;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
    }

    double _r = 1.0/sqrt(r2);

    double r = r2*_r;
    double ar = alpha*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
#ifdef CPPMD_ENABLE_TABLE_FOR_EWALD_REAL
# ifdef K_SIMD
    int iar = int((ar)*4.0);
    double dar = ((ar)*4.0-iar)*2.0-1.0;
    ee = ptc[iar][0]+dar*(ptc[iar][1]+dar*(ptc[iar][2]+dar*(ptc[iar][3]+dar*ptc[iar][4])));
    edp = ftc[iar][0]+dar*(ftc[iar][1]+dar*(ftc[iar][2]+dar*(ftc[iar][3]+dar*ftc[iar][4])));
# else
    ewaldrealforceandpottbl_(&ee,&edp,&ar);
# endif
#else
    ee = erfc(ar);
    edp = M_2_SQRTPI*ar*exp(-ar*ar)+ee;
#endif
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6;
    double esingle = cij     *(_r*ee);
    ljpair_energy += ljsingle*inside;
    shortpair_energy += esingle*inside;
    e += ( ljsingle + esingle )*inside;
# else
    e += ( ljmp[atj][0]*(_r12)
         - ljmp[atj][1]*(_r6)
         + cij     *(_r*ee)
           )*inside;
# endif
    double dp = (ljmp[atj][2]*(_r14)
               - ljmp[atj][3]*(_r8)
               + cij     *(_r3*edp)
                 )*inside;
    v += dp*r2;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
  virial = v;
}

template<>
void pairlistloopljcfe_ewald(const int (*jindex)[MAX_PAIR],
                             const int (*lj)[MAX_PAIR],
                             const int *npair,
                             const int *iindex,
                             const ParticleArray& particlei,
                             const ParticleArray& particlej,
                             const int npl,
                             const double cutoff2,
                             double (*force)[3],
                             double& energy,
			     double& virial,
			     const OperationSelector operations)
{
  int i;
  double e=energy;
  double v=virial;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {particlei[iid].position.x,
                      particlei[iid].position.y,
                      particlei[iid].position.z};
    double chargei = particlei[iid].charge;
    int atomi;
    atomi = particlei[iid].atomtype;
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doVirialcalculation){
      pairlistloopljcfe1_ewald(ShortRange::alpha,
			       jindex[i],
			       lj[i],npair[i],
			       posi, chargei, 
			       ljmp,
			       particlej, 
			       cutoff2,
			       force[i],
			       e,
			       v);
    }else{
      if(operations.doEnergycalculation){
	pairlistloopljcfe1_ewald(ShortRange::alpha,
				 jindex[i],
				 lj[i],npair[i],
				 posi, chargei, 
				 ljmp,
				 particlej, 
				 cutoff2,
				 force[i],
				 e);
      }else{
	pairlistloopljcf1_ewald(ShortRange::alpha,
				jindex[i],
				lj[i],npair[i],
				posi, chargei, 
				ljmp,
				particlej, 
				cutoff2,
				force[i] );
      }
    }
#ifdef OVERLAP
# ifdef OVERLAP_PAIRLIST_THREAD
    {
      int tn = 0;
#ifdef _OPENMP
      tn = omp_get_thread_num();
#endif
      if(tn==0){
	MPICHECKER::mpi_checker();
      }
    }
# endif
#endif
  }
  energy = e;
  virial = v;
}

template<>
void pairlistloopljcfe_ewald(const int (*jindex)[MAX_PAIR],
                             const int (*lj)[MAX_PAIR],
                             const int *npair,
                             const int *iindex,
                             const ParticleArray& particlei,
                             const ParticleArray& particlej,
                             const int npl,
                             const int *selfnpair,
                             const double cutoff2,
                             double (*force)[3],
                             double& energy,
			     double& virial,
			     const OperationSelector operations)
{
  int i;
  double e=energy;
  double v=virial;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,v) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
    int iid = iindex[i];
    double posi[3] = {particlei[iid].position.x,
                      particlei[iid].position.y,
                      particlei[iid].position.z};
    double chargei = particlei[iid].charge;
    int atomi;
    atomi = particlei[iid].atomtype;
    double (*ljmp)[4] = (double (*)[4])(&(ShortRange::ljmixparameters[atomi][0].potFirst));
    if(operations.doVirialcalculation){
      pairlistloopljcfe1_ewald(ShortRange::alpha,
			       &(jindex[i][selfnpair[i]]),
			       &(lj[i][selfnpair[i]]),npair[i],
			       posi, chargei, 
			       ljmp,
			       particlej, 
			       cutoff2,
			       force[i],
			       e,
			       v);
    }else{
      if(operations.doEnergycalculation){
	pairlistloopljcfe1_ewald(ShortRange::alpha,
				 &(jindex[i][selfnpair[i]]),
				 &(lj[i][selfnpair[i]]),npair[i],
				 posi, chargei, 
				 ljmp,
				 particlej, 
				 cutoff2,
				 force[i],
				 e);
      }else{
	pairlistloopljcf1_ewald(ShortRange::alpha,
				&(jindex[i][selfnpair[i]]),
				&(lj[i][selfnpair[i]]),npair[i],
				posi, chargei, 
				ljmp,
				particlej, 
				cutoff2,
				force[i] );
      }
    }
  }
  energy = e;
  virial = v;
}




// pair list holds value of pos, charge
void pairlistloopljcfe1(const double (*pos)[3],
                        const double *charge,
                        const double (*lj)[4],
                        const int npair,
                        const double posi[3],
                        const double chargei,
                        const double cutoff2,
                        double force[3],
                        double& energy)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
#ifdef SIMPLE_CUTOFF
#else
  double &r1co0 = ShortRange::r1co0, &r1co1=ShortRange::r1co1, &r1co2=ShortRange::r1co2;
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r3co1=ShortRange::r3co1, &r3co2=ShortRange::r3co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
#endif

  for(j=0;j<npair;j++){
    double dx, dy, dz;
    dx = posi[0]-pos[j][0];
    dy = posi[1]-pos[j][1];
    dz = posi[2]-pos[j][2];
    double r2 = dx*dx+dy*dy+dz*dz;
    double cij = chargei*charge[j];

#ifdef SIMPLE_CUTOFF
    double inside = 1.0;
    if(r2>=cutoff2){
      inside=0.0;
    }else{
      inside=1.0;
    }
#else
    if(r2>cutoff2)r2=cutoff2;
#endif

    double _r = 1.0/sqrt(r2);

#ifdef SIMPLE_CUTOFF
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
# ifdef CHECK_ENERGY
    double ljsingle = lj[j][0]*_r12 - lj[j][1]*_r6;
    double esingle = cij     *_r;
    ljpair_energy += ljsingle*inside;
    shortpair_energy += esingle*inside;
    e += (ljsingle + esingle)*inside;
# else
    e += (lj[j][0]*_r12 - lj[j][1]*_r6 + cij*_r)*inside;
# endif
    double dp = (lj[j][2]*_r14 - lj[j][3]*_r8 + cij*_r3)*inside;
#else
    double r = r2*_r;
    double r3 = r2*r;
    double _r2 = _r*_r;
    double _r3 = _r2*_r;
    double _r6 = _r2*_r2*_r2;
    double _r8 = _r6*_r2;
    double _r12 = _r6*_r6;
    double _r14 = _r8*_r6;
# ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
    double esingle = cij     *(_r  +(r3*(r1co1 -r1co2 *r)-r1co0)) ;
    ljpair_energy += ljsingle;
    shortpair_energy += esingle;
    e += ( ljsingle + esingle );
# else
    e += ( lj[j][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
         - lj[j][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
         + cij     *(_r  +(r3*(r1co1 -r1co2 *r)-r1co0 ))
           );
# endif
    double dp = (lj[j][2]*(_r14-r*(r14co1-r14co2*r))
               - lj[j][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3 -r*(r3co1 -r3co2 *r))
                 );
#endif
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
}

void pairlistloopljcfe(const double (*pos)[MAX_PAIR][3],
                       const double (*charge)[MAX_PAIR],
                       const double (*lj)[MAX_PAIR][4],
                       const int *npair,
                       const double (*posi)[3],
                       const double *chargei,
                       const int npl,
                       const double cutoff2,
                       double (*force)[3],
                       double& energy)
{
  int i;
  double e=energy;
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,ljpair_energy,shortpair_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,ljpair_energy,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,ljpair_energy,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e) schedule(dynamic)
#   endif
#  endif
# endif
#endif
  for(i=0;i<npl;i++){
#ifdef CPPMD_ENABLE_FORTRAN_KERNEL
    pairlistloopljcfe1_f_(pos[i],charge[i],lj[i],&npair[i],
                       posi[i], &chargei[i], &cutoff2,
                       force[i],
                       &e);
#else
    pairlistloopljcfe1(pos[i],charge[i],lj[i],npair[i],
                       posi[i], chargei[i], cutoff2,
                       force[i],
                       e);
#endif
#ifdef OVERLAP
# ifdef OVERLAP_PAIRLIST_THREAD
    {
      int tn = 0;
#ifdef _OPENMP
      tn = omp_get_thread_num();
#endif
      if(tn==0){
	MPICHECKER::mpi_checker();
      }
    }
# endif
#endif
  }
  energy = e;
}


void coulomb_jloop(const double (*pacs)[4], const int nj, const double posi[3], const double chargei, double& forcex, double& forcey, double& forcez, double& energy)
{
  double e = 0.0;
  double fx=0.0, fy=0.0, fz=0.0;
  int j;
#pragma loop norecurrence
#pragma loop noalias
#pragma loop reduction
#pragma loop novrec
#pragma loop simd
  for(j=0; j<nj;j++){
    double dx =  posi[0]-pacs[j][0];
    double dy =  posi[1]-pacs[j][1];
    double dz =  posi[2]-pacs[j][2];
    double cij = chargei*pacs[j][3];
    double r2 = dx*dx+dy*dy+dz*dz;
    double _r = 1.0/sqrt(r2);
    double cr = cij*_r;
    double dp = cr*_r*_r;
# ifdef CHECK_ENERGY
    shortpair_energy += cr;
# endif
    e += cr;
    fx += dp*dx;
    fy += dp*dy;
    fz += dp*dz;
  }
  energy += e;
  forcex += fx;
  forcey += fy;
  forcez += fz;
}

template<typename GPA>
void cellpairinteraction1(const ParticlePosCharge posci,
			  const GPA& particlej,
			  const std::vector<TypeRange>& typerangej,
			  const std::vector<int>& shorttarget_index,
			  Force& force,
			  const int icell, const int i,
			  const bool self)
{
  const double (*pacs)[4] = (double (*)[4])(&(particlej.poscharge[0].position.x));
  const double *posi = (double *)(&(posci.position.x));
  const double &chargei = posci.charge;
  double fx = 0.0, fy =0.0, fz=0.0;
  //  double e = 0.0;
  for(int cj=0;cj<shorttarget_index.size();cj++){
    int jcell = shorttarget_index[cj];
    int jmin = typerangej[jcell].begin;
    int jmax = typerangej[jcell].end;
    if(self&&(jcell==icell)){
#pragma loop unroll 4
      for(int j=jmin; j<i;j++){
	double dx =  posi[0]-pacs[j][0];
	double dy =  posi[1]-pacs[j][1];
	double dz =  posi[2]-pacs[j][2];
	double cij = chargei*pacs[j][3];
	double r2 = dx*dx+dy*dy+dz*dz;
	double _r = 1.0/sqrt(r2);
	double cr = cij*_r;
	double dp = cr*_r*_r;
	//	e += cr;
	fx += dp*dx;
	fy += dp*dy;
	fz += dp*dz;
      }
#pragma loop unroll 4
      for(int j=i+1; j<jmax;j++){
	double dx =  posi[0]-pacs[j][0];
	double dy =  posi[1]-pacs[j][1];
	double dz =  posi[2]-pacs[j][2];
	double cij = chargei*pacs[j][3];
	double r2 = dx*dx+dy*dy+dz*dz;
	double _r = 1.0/sqrt(r2);
	double cr = cij*_r;
	double dp = cr*_r*_r;
	//	e += cr;
	fx += dp*dx;
        fy += dp*dy;
	fz += dp*dz;
      }
    }else{
#pragma loop unroll 4
      for(int j=jmin; j<jmax;j++){
	double dx =  posi[0]-pacs[j][0];
	double dy =  posi[1]-pacs[j][1];
	double dz =  posi[2]-pacs[j][2];
	double cij = chargei*pacs[j][3];
	double r2 = dx*dx+dy*dy+dz*dz;
	double _r = 1.0/sqrt(r2);
	double cr = cij*_r;
	double dp = cr*_r*_r;
	//	e += cr;
	fx += dp*dx;
	fy += dp*dy;
	fz += dp*dz;
      }
    }
  }
  force.x += fx;
  force.y += fy;
  force.z += fz;
  //  energy += e;
}

template<typename GPA>
void cellpairinteraction1(const ParticlePosCharge posci,
			  const GPA& particlej,
			  const std::vector<TypeRange>& typerangej,
			  const std::vector<int>& shorttarget_index,
			  Force& force,
			  double& energy,
			  const int icell, const int i,
			  const bool self)
{
  const double (*pacs)[4] = (double (*)[4])(&(particlej.poscharge[0].position.x));
  const double *posi = (double *)(&(posci.position.x));
  const double &chargei = posci.charge;
  double fx = 0.0, fy =0.0, fz=0.0;
  double e = 0.0;
  for(int cj=0;cj<shorttarget_index.size();cj++){
    int jcell = shorttarget_index[cj];
    int jmin = typerangej[jcell].begin;
    int jmax = typerangej[jcell].end;
    if(self&&(jcell==icell)){
#pragma loop unroll 4
      for(int j=jmin; j<i;j++){
	double dx =  posi[0]-pacs[j][0];
	double dy =  posi[1]-pacs[j][1];
	double dz =  posi[2]-pacs[j][2];
	double cij = chargei*pacs[j][3];
	double r2 = dx*dx+dy*dy+dz*dz;
	double _r = 1.0/sqrt(r2);
	double cr = cij*_r;
	double dp = cr*_r*_r;
# ifdef CHECK_ENERGY
	shortpair_energy += cr;
# endif
	e += cr;
	fx += dp*dx;
	fy += dp*dy;
	fz += dp*dz;
      }
#pragma loop unroll 4
      for(int j=i+1; j<jmax;j++){
	double dx =  posi[0]-pacs[j][0];
	double dy =  posi[1]-pacs[j][1];
	double dz =  posi[2]-pacs[j][2];
	double cij = chargei*pacs[j][3];
	double r2 = dx*dx+dy*dy+dz*dz;
	double _r = 1.0/sqrt(r2);
	double cr = cij*_r;
	double dp = cr*_r*_r;
# ifdef CHECK_ENERGY
	shortpair_energy += cr;
# endif
	e += cr;
	fx += dp*dx;
        fy += dp*dy;
	fz += dp*dz;
      }
    }else{
#pragma loop unroll 4
      for(int j=jmin; j<jmax;j++){
	double dx =  posi[0]-pacs[j][0];
	double dy =  posi[1]-pacs[j][1];
	double dz =  posi[2]-pacs[j][2];
	double cij = chargei*pacs[j][3];
	double r2 = dx*dx+dy*dy+dz*dz;
	double _r = 1.0/sqrt(r2);
	double cr = cij*_r;
	double dp = cr*_r*_r;
# ifdef CHECK_ENERGY
	shortpair_energy += cr;
# endif
	e += cr;
	fx += dp*dx;
	fy += dp*dy;
	fz += dp*dz;
      }
    }
  }
  force.x += fx;
  force.y += fy;
  force.z += fz;
  energy += e;
}
		
#define CELLPAIR_OMP_STATIC
template<typename GPA>
void cellpairsinteraction(const CombinedParticleArray& particlei,
			  const std::vector<TypeRange>& typerangei,
			  const GPA& particlej,
			  const std::vector<TypeRange>& typerangej,
			  const std::vector<std::vector<int> >& shorttarget_index,
			  ForceArray& force,
			  double& energy,
			  const bool self,
			  const OperationSelector operations)
{
  //DEBUG
  /*
  {
    printf("target cell index for cell 0");
    for(int c=0;c<shorttarget_index[0].size();c++){
      printf(" %d",shorttarget_index[0][c]);
    }
    printf("\n");
  }
  */
  int cimax = typerangei.size();
  for(int ci=0;ci<cimax;ci++){
    int i;
    int imin, imax;
    double e = 0.0;
    imin = typerangei[ci].begin;
    imax = typerangei[ci].end;
    if(operations.doEnergycalculation){
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef CELLPAIR_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,shortpair_energy) schedule(runtime)
#  else
#   ifdef CELLPAIR_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef CELLPAIR_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e) schedule(runtime)
#  else
#   ifdef CELLPAIR_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e) schedule(dynamic)
#   endif
#  endif
# endif
#endif
      for(i=imin;i<imax;i++){
	cellpairinteraction1(particlei.poscharge[i],
			     particlej,
			     typerangej,
			     shorttarget_index[ci],
			     force[i],e,
			     ci,i,self);
      }
      energy += e;
    }else{
#ifdef _OPENMP
# ifdef CELLPAIR_OMP_RUNTIME
#pragma omp parallel for private(i) schedule(runtime)
# else
#  ifdef CELLPAIR_OMP_STATIC
#pragma omp parallel for private(i) schedule(static)
#  else
#pragma omp parallel for private(i) schedule(dynamic)
#  endif
# endif
#endif
      for(i=imin;i<imax;i++){
	cellpairinteraction1(particlei.poscharge[i],
			     particlej,
			     typerangej,
			     shorttarget_index[ci],
			     force[i],
			     ci,i,self);
      }
    }
#ifdef OVERLAP
    MPICHECKER::mpi_checker();
#endif
  }
}

template
void cellpairsinteraction(const CombinedParticleArray& particlei,
			  const std::vector<TypeRange>& typerangei,
			  const CombinedParticleArray& particlej,
			  const std::vector<TypeRange>& typerangej,
			  const std::vector<std::vector<int> >& shorttarget_index,
			  ForceArray& force,
			  double& energy,
			  const bool self,
			  const OperationSelector operations);

template
void cellpairsinteraction(const CombinedParticleArray& particlei,
			  const std::vector<TypeRange>& typerangei,
			  const GhostParticleArray& particlej,
			  const std::vector<TypeRange>& typerangej,
			  const std::vector<std::vector<int> >& shorttarget_index,
			  ForceArray& force,
			  double& energy,
			  const bool self,
			  const OperationSelector operations);





template<typename GPA>
void cellpairinteraction1(const ParticlePosCharge posci,
			  const std::vector<waterexclude> &waterexcludelist,
			  const GPA& particlej,
			  const std::vector<TypeRange>& typerangej,
			  const std::vector<int>& shorttarget_index,
			  Force& force,
			  const int icell, const int i,
			  const bool self)
{
  const double (*pacs)[4] = (double (*)[4])(&(particlej.poscharge[0].position.x));
  const double *posi = (double *)(&(posci.position.x));
  const double &chargei = posci.charge;
  double fx = 0.0, fy =0.0, fz=0.0;
  //  double e = 0.0;
  for(int cj=0;cj<shorttarget_index.size();cj++){
    int jcell = shorttarget_index[cj];
    int jmin = typerangej[jcell].begin;
    int jmax = typerangej[jcell].end;
    if(self&&(jcell==icell)){
#pragma loop unroll 4
      for(int j=jmin; j<i;j++){
	if((j==waterexcludelist[i].exid[0])||(j==waterexcludelist[i].exid[1]))continue;
	double dx =  posi[0]-pacs[j][0];
	double dy =  posi[1]-pacs[j][1];
	double dz =  posi[2]-pacs[j][2];
	double cij = chargei*pacs[j][3];
	double r2 = dx*dx+dy*dy+dz*dz;
	double _r = 1.0/sqrt(r2);
	double cr = cij*_r;
	double dp = cr*_r*_r;
	//	e += cr;
	fx += dp*dx;
	fy += dp*dy;
	fz += dp*dz;
      }
#pragma loop unroll 4
      for(int j=i+1; j<jmax;j++){
	if((j==waterexcludelist[i].exid[0])||(j==waterexcludelist[i].exid[1]))continue;
	double dx =  posi[0]-pacs[j][0];
	double dy =  posi[1]-pacs[j][1];
	double dz =  posi[2]-pacs[j][2];
	double cij = chargei*pacs[j][3];
	double r2 = dx*dx+dy*dy+dz*dz;
	double _r = 1.0/sqrt(r2);
	double cr = cij*_r;
	double dp = cr*_r*_r;
	//	e += cr;
	fx += dp*dx;
        fy += dp*dy;
	fz += dp*dz;
      }
    }else{
#pragma loop unroll 4
      for(int j=jmin; j<jmax;j++){
	double dx =  posi[0]-pacs[j][0];
	double dy =  posi[1]-pacs[j][1];
	double dz =  posi[2]-pacs[j][2];
	double cij = chargei*pacs[j][3];
	double r2 = dx*dx+dy*dy+dz*dz;
	double _r = 1.0/sqrt(r2);
	double cr = cij*_r;
	double dp = cr*_r*_r;
	//	e += cr;
	fx += dp*dx;
	fy += dp*dy;
	fz += dp*dz;
      }
    }
  }
  force.x += fx;
  force.y += fy;
  force.z += fz;
  //  energy += e;
}

template<typename GPA>
void cellpairinteraction1(const ParticlePosCharge posci,
			  const std::vector<waterexclude> &waterexcludelist,
			  const GPA& particlej,
			  const std::vector<TypeRange>& typerangej,
			  const std::vector<int>& shorttarget_index,
			  Force& force,
			  double& energy,
			  const int icell, const int i,
			  const bool self)
{
  const double (*pacs)[4] = (double (*)[4])(&(particlej.poscharge[0].position.x));
  const double *posi = (double *)(&(posci.position.x));
  const double &chargei = posci.charge;
  double fx = 0.0, fy =0.0, fz=0.0;
  double e = 0.0;
  for(int cj=0;cj<shorttarget_index.size();cj++){
    int jcell = shorttarget_index[cj];
    int jmin = typerangej[jcell].begin;
    int jmax = typerangej[jcell].end;
    if(self&&(jcell==icell)){
#pragma loop unroll 4
      for(int j=jmin; j<i;j++){
	double dx =  posi[0]-pacs[j][0];
	double dy =  posi[1]-pacs[j][1];
	double dz =  posi[2]-pacs[j][2];
	double cij = chargei*pacs[j][3];
	double r2 = dx*dx+dy*dy+dz*dz;
	double _r = 1.0/sqrt(r2);
	double cr = cij*_r;
	double dp = cr*_r*_r;
# ifdef CHECK_ENERGY
	shortpair_energy += cr;
# endif
	e += cr;
	if((j==waterexcludelist[i].exid[0])||(j==waterexcludelist[i].exid[1])){
	  continue;
	}
	fx += dp*dx;
	fy += dp*dy;
	fz += dp*dz;
      }
#pragma loop unroll 4
      for(int j=i+1; j<jmax;j++){
	double dx =  posi[0]-pacs[j][0];
	double dy =  posi[1]-pacs[j][1];
	double dz =  posi[2]-pacs[j][2];
	double cij = chargei*pacs[j][3];
	double r2 = dx*dx+dy*dy+dz*dz;
	double _r = 1.0/sqrt(r2);
	double cr = cij*_r;
	double dp = cr*_r*_r;
# ifdef CHECK_ENERGY
	shortpair_energy += cr;
# endif
	e += cr;
	if((j==waterexcludelist[i].exid[0])||(j==waterexcludelist[i].exid[1])){
	  continue;
	}
	fx += dp*dx;
        fy += dp*dy;
	fz += dp*dz;
      }
    }else{
#pragma loop unroll 4
      for(int j=jmin; j<jmax;j++){
	double dx =  posi[0]-pacs[j][0];
	double dy =  posi[1]-pacs[j][1];
	double dz =  posi[2]-pacs[j][2];
	double cij = chargei*pacs[j][3];
	double r2 = dx*dx+dy*dy+dz*dz;
	double _r = 1.0/sqrt(r2);
	double cr = cij*_r;
	double dp = cr*_r*_r;
# ifdef CHECK_ENERGY
	shortpair_energy += cr;
# endif
	e += cr;
	fx += dp*dx;
	fy += dp*dy;
	fz += dp*dz;
      }
    }
  }
  force.x += fx;
  force.y += fy;
  force.z += fz;
  energy += e;
}

template<typename GPA>
void cellpairsinteraction(const CombinedParticleArray& particlei,
			  const std::vector<TypeRange>& typerangei,
			  const std::vector<waterexclude> &waterexcludelist,
			  const GPA& particlej,
			  const std::vector<TypeRange>& typerangej,
			  const std::vector<std::vector<int> >& shorttarget_index,
			  ForceArray& force,
			  double& energy,
			  const bool self,
			  const OperationSelector operations)
{
  //DEBUG
  /*
  {
    printf("target cell index for cell 0");
    for(int c=0;c<shorttarget_index[0].size();c++){
      printf(" %d",shorttarget_index[0][c]);
    }
    printf("\n");
  }
  */
  int cimax = typerangei.size();
  for(int ci=0;ci<cimax;ci++){
    int i;
    int imin, imax;
    double e = 0.0;
    imin = typerangei[ci].begin;
    imax = typerangei[ci].end;
    if(operations.doEnergycalculation){
#ifdef _OPENMP
# ifdef CHECK_ENERGY
#  ifdef CELLPAIR_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e,shortpair_energy) schedule(runtime)
#  else
#   ifdef CELLPAIR_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,shortpair_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,shortpair_energy) schedule(dynamic)
#   endif
#  endif
# else
#  ifdef CELLPAIR_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: e) schedule(runtime)
#  else
#   ifdef CELLPAIR_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e) schedule(dynamic)
#   endif
#  endif
# endif
#endif
      for(i=imin;i<imax;i++){
	cellpairinteraction1(particlei.poscharge[i],
			     waterexcludelist,
			     particlej,
			     typerangej,
			     shorttarget_index[ci],
			     force[i],e,
			     ci,i,self);
      }
      energy += e;
    }else{
#ifdef _OPENMP
# ifdef CELLPAIR_OMP_RUNTIME
#pragma omp parallel for private(i) schedule(runtime)
# else
#  ifdef CELLPAIR_OMP_STATIC
#pragma omp parallel for private(i) schedule(static)
#  else
#pragma omp parallel for private(i) schedule(dynamic)
#  endif
# endif
#endif
      for(i=imin;i<imax;i++){
	cellpairinteraction1(particlei.poscharge[i],
			     waterexcludelist,
			     particlej,
			     typerangej,
			     shorttarget_index[ci],
			     force[i],
			     ci,i,self);
      }
    }
#ifdef OVERLAP
    MPICHECKER::mpi_checker();
#endif
  }
}

template
void cellpairsinteraction(const CombinedParticleArray& particlei,
			  const std::vector<TypeRange>& typerangei,
			  const std::vector<waterexclude> &waterexcludelist,
			  const CombinedParticleArray& particlej,
			  const std::vector<TypeRange>& typerangej,
			  const std::vector<std::vector<int> >& shorttarget_index,
			  ForceArray& force,
			  double& energy,
			  const bool self,
			  const OperationSelector operations);

template
void cellpairsinteraction(const CombinedParticleArray& particlei,
			  const std::vector<TypeRange>& typerangei,
			  const std::vector<waterexclude> &waterexcludelist,
			  const GhostParticleArray& particlej,
			  const std::vector<TypeRange>& typerangej,
			  const std::vector<std::vector<int> >& shorttarget_index,
			  ForceArray& force,
			  double& energy,
			  const bool self,
			  const OperationSelector operations);

