#include <cmath>
#include "PairListInteraction.h"
#include "ShortRangeInteraction.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/*
  Zero Dipole Sum
  Chemical Physics Letters 568-569 (2013) 26-32
  E = \sum_i \sum_j q_i q_j [U(r_ij)-u(r_c)] - [u(r_c)/2 + \alpha/sqrt(\pi)] \sum_i q_i^2  (3)

  u(r) = erfc(\alpha r)/r + [erfc(\alpha r)/(2r_c) + \alpha/sqrt(\pi) exp(-(\alpha r_c)^2)]r^2/r_c^2 (4)


  A = [erfc(\alpha r)/(2r_c) + \alpha/sqrt(\pi) exp(-(\alpha r_c)^2)]
  u(r) = erfc(\alpha r)/r + A r^2/r_c^2
       = [erfc(\alpha r) + A r^3/r_c^2]/r
  u(r_c) = [erfc(\alpha r_c) + A r_c]/r_c
         = erfc(\alpha r_c)/r_c + erfc(\alpha r)/(2r_c) + \alpha/sqrt(\pi) exp(-(\alpha r_c)^2)
         = 3/2 erfc(\alpha r_c)/r_c + \alpha/sqrt(\pi) exp(-(\alpha r_c)^2)

  -du(r)/dr = erfc(\alpha r)/r^2 + 2 \alpha/sqrt(\pi) exp(-(\alpha r)^2)/r - A2r/r_c^2
            = [erfc(\alpha r) + 2 \alpha r/sqrt(\pi) exp(-(\alpha r)^2) - 2Ar^3]/r^2


  When \alpha = 0,  erfc(\alpha r)=1, exp(-(\alpha r)^2)=1, A=1/2r_c, 
  u(r) = (1 + r^3/2r_c^3)/r
  -du(r)/dr = (1 - r^3/r_c^3)/r^2
 */

static bool calculated_self_energy_zerodipole = false;
static double self_energy_zerodipole = 0.0;

#ifdef CHECK_ENERGY
double zerodipole_energy;
#endif

template<typename PA>
double calc_self_energy_zerodipole(const int *iindex,
				   const PA& particlei,
				   const int npl,
				   const double cutoff2)
{
  if(calculated_self_energy_zerodipole==false){
    int i;
    double energy = self_energy_zerodipole;
    const double alpharc2 = ShortRange::alpha*ShortRange::alpha*cutoff2;
    const double alpharc = ShortRange::alpha*sqrt(cutoff2);
    const double ec2 = 0.5*ShortRange::alpha*M_2_SQRTPI*exp(-alpharc2);
    const double ec1 = 0.5*erfc(alpharc)/sqrt(cutoff2);
    const double el = (ec1+ec2)/cutoff2;
    const double ec = -3.0*ec1 - ec2;
    const double ese = (ec - ShortRange::alpha*M_2_SQRTPI*0.5)*0.5;

#ifdef _OPENMP
# ifdef PL_OMP_RUNTIME
#pragma omp parallel for private(i) reduction(+: energy) schedule(runtime)
# else
#  ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: energy) schedule(static)
#  else
#pragma omp parallel for private(i) reduction(+: energy) schedule(dynamic)
#  endif
# endif
#endif
    for(i=0;i<npl;i++){
      int iid = iindex[i];
      double chargei = particlei.poscharge[iid].charge;
      energy += ese*chargei*chargei;
    }
    calculated_self_energy_zerodipole = true;
    self_energy_zerodipole = energy;
#ifdef DUMP_ENERGY_CONTENT
    printf("selfenergy at ZD %f with %d atom factor %f alpha %f\n",energy,npl,ese,ShortRange::alpha);
#endif
  }
  return self_energy_zerodipole;
}

template
double calc_self_energy_zerodipole(const int *iindex,
				   const CombinedParticleArray& particlei,
				   const int npl,
				   const double cutoff2);


void pairlistloopljcf1_zerodipole0(const double alpha,
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
  //  double v=virial;
  const double ec1 = 0.5/sqrt(cutoff2);
  const double el = ec1/cutoff2;
  const double ec = -3.0*ec1;

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
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

    // e += ( ljmp[atj][0]*(_r12)
    //      - ljmp[atj][1]*(_r6)
    // 	   + cij     *(_r + el*r2 + ec)
    // 	   )*inside;
    double dp = (ljmp[atj][2]*(_r14)
               - ljmp[atj][3]*(_r8)
               + cij     *(_r3 - 2.0*el)
                 )*inside;
    //    v += dp*r2;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  //  energy = e;
  //  virial = v;
}

void pairlistloopljcfe1_zerodipole0(const double alpha,
				    const int *jindex,
				    const int *lj,
				    const int npair,
				    const double posi[3],
				    const double chargei,
				    const double (*ljmp)[4],
				    const PosChargeArray& particlej,
				    const double cutoff2,
				    double force[3],
				    double &energy)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  //  double v=virial;
  const double ec1 = 0.5/sqrt(cutoff2);
  const double el = ec1/cutoff2;
  const double ec = -3.0*ec1;

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
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

#ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6;
    double esingle = cij     *(_r + el*r2 + ec);
    ljpair_energy += ljsingle*inside;
    zerodipole_energy += esingle*inside;
    e += ( ljsingle + esingle )*inside;
#else
    e += ( ljmp[atj][0]*(_r12)
         - ljmp[atj][1]*(_r6)
	   + cij     *(_r + el*r2 + ec)
	   )*inside;
#endif
    double dp = (ljmp[atj][2]*(_r14)
               - ljmp[atj][3]*(_r8)
               + cij     *(_r3 - 2.0*el)
                 )*inside;
    //    v += dp*r2;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
  //  virial = v;
}

void pairlistloopljcfe1_zerodipole0(const double alpha,
				    const int *jindex,
				    const int *lj,
				    const int npair,
				    const double posi[3],
				    const double chargei,
				    const double (*ljmp)[4],
				    const PosChargeArray& particlej,
				    const double cutoff2,
				    double force[3],
				    double &energy,
				    double &virial)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  double v=virial;
  const double ec1 = 0.5/sqrt(cutoff2);
  const double el = ec1/cutoff2;
  const double ec = -3.0*ec1;

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
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

#ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6;
    double esingle = cij     *(_r + el*r2 + ec);
    ljpair_energy += ljsingle*inside;
    zerodipole_energy += esingle*inside;
    e += ( ljsingle + esingle )*inside;
#else
    e += ( ljmp[atj][0]*(_r12)
         - ljmp[atj][1]*(_r6)
	   + cij     *(_r + el*r2 + ec)
	   )*inside;
#endif
    double dp = (ljmp[atj][2]*(_r14)
               - ljmp[atj][3]*(_r8)
               + cij     *(_r3 - 2.0*el)
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


void pairlistloopljcf1_zerodipole0_ljshift(const double alpha,
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
  //  double v=virial;
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
  const double ec1 = 0.5/sqrt(cutoff2);
  const double el = ec1/cutoff2;
  const double ec = -3.0*ec1;

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
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

    // e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
    //      - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
    // 	   + cij     *(_r + el*r2 + ec)
    // 	   )*inside;
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3 - 2.0*el)
                 )*inside;
    //    v += dp*r2;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  //  energy = e;
  //  virial = v;
}

void pairlistloopljcfe1_zerodipole0_ljshift(const double alpha,
					   const int *jindex,
					   const int *lj,
					   const int npair,
					   const double posi[3],
					   const double chargei,
					   const double (*ljmp)[4],
					   const PosChargeArray& particlej,
					   const double cutoff2,
					   double force[3],
					   double &energy)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  //  double v=virial;
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
  const double ec1 = 0.5/sqrt(cutoff2);
  const double el = ec1/cutoff2;
  const double ec = -3.0*ec1;

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
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

#ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
    double esingle = cij     *(_r + el*r2 + ec);
    ljpair_energy += ljsingle*inside;
    zerodipole_energy += esingle*inside;
    e += ( ljsingle + esingle )*inside;
#else
    e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
         - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
	   + cij     *(_r + el*r2 + ec)
	   )*inside;
#endif
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3 - 2.0*el)
                 )*inside;
    //    v += dp*r2;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
  //  virial = v;
}

void pairlistloopljcfe1_zerodipole0_ljshift(const double alpha,
					   const int *jindex,
					   const int *lj,
					   const int npair,
					   const double posi[3],
					   const double chargei,
					   const double (*ljmp)[4],
					   const PosChargeArray& particlej,
					   const double cutoff2,
					   double force[3],
					   double &energy,
					   double &virial)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  double v=virial;
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
  const double ec1 = 0.5/sqrt(cutoff2);
  const double el = ec1/cutoff2;
  const double ec = -3.0*ec1;

#ifdef K_SIMD
  double (*pacs)[4] = (double (*)[4])(&(particlej[0].position.x));
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

#ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
    double esingle = cij     *(_r + el*r2 + ec);
    ljpair_energy += ljsingle*inside;
    zerodipole_energy += esingle*inside;
    e += ( ljsingle + esingle )*inside;
#else
    e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
         - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
	   + cij     *(_r + el*r2 + ec)
	   )*inside;
#endif
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3 - 2.0*el)
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

void pairlistloopljcf1_zerodipole(const double alpha,
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
  const double alpharc2 = alpha*alpha*cutoff2;
  const double alpharc = alpha*sqrt(cutoff2);
  const double ec2 = 0.5*M_2_SQRTPI*alpha*exp(-alpharc2);
  const double ec1 = 0.5*erfc(alpharc)/sqrt(cutoff2);
  const double el = (ec1+ec2)/cutoff2;
  const double ec = -3.0*ec1 - ec2;

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
    //         + cij     *(_r*ee + el*r2 + ec)
    //           )*inside;
    double dp = (ljmp[atj][2]*(_r14)
               - ljmp[atj][3]*(_r8)
               + cij     *(_r3*edp-2.0*el)
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

void pairlistloopljcfe1_zerodipole(const double alpha,
				   const int *jindex,
				   const int *lj,
				   const int npair,
				   const double posi[3],
				   const double chargei,
				   const double (*ljmp)[4],
				   const PosChargeArray& particlej,
				   const double cutoff2,
				   double force[3],
				   double &energy)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  double ee, edp;
  const double alpharc2 = alpha*alpha*cutoff2;
  const double alpharc = alpha*sqrt(cutoff2);
  const double ec2 = 0.5*M_2_SQRTPI*alpha*exp(-alpharc2);
  const double ec1 = 0.5*erfc(alpharc)/sqrt(cutoff2);
  const double el = (ec1+ec2)/cutoff2;
  const double ec = -3.0*ec1 - ec2;

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
#ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6;
    double esingle = cij     *(_r*ee + el*r2 + ec);
    ljpair_energy += ljsingle*inside;
    zerodipole_energy += esingle*inside;
    e += ( ljsingle + esingle
	   )*inside;
#else
    e += ( ljmp[atj][0]*(_r12)
	   - ljmp[atj][1]*(_r6)
	   + cij     *(_r*ee + el*r2 + ec)
	   )*inside;
#endif
    double dp = (ljmp[atj][2]*(_r14)
               - ljmp[atj][3]*(_r8)
               + cij     *(_r3*edp-2.0*el)
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

void pairlistloopljcfe1_zerodipole(const double alpha,
				   const int *jindex,
				   const int *lj,
				   const int npair,
				   const double posi[3],
				   const double chargei,
				   const double (*ljmp)[4],
				   const PosChargeArray& particlej,
				   const double cutoff2,
				   double force[3],
				   double &energy,
				   double &virial)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  double v=virial;
  double ee, edp;
  const double alpharc2 = alpha*alpha*cutoff2;
  const double alpharc = alpha*sqrt(cutoff2);
  const double ec2 = 0.5*M_2_SQRTPI*alpha*exp(-alpharc2);
  const double ec1 = 0.5*erfc(alpharc)/sqrt(cutoff2);
  const double el = (ec1+ec2)/cutoff2;
  const double ec = -3.0*ec1 - ec2;

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
#ifdef CHECK_ENERGY
    double esingle = cij     *(_r*ee + el*r2 + ec);
    double ljsingle = ljmp[atj][0]*_r12 - ljmp[atj][1]*_r6;
    ljpair_energy += ljsingle*inside;
    zerodipole_energy += esingle*inside;
    e += ( ljsingle + esingle )*inside;
#else
    e += ( ljmp[atj][0]*(_r12)
	   - ljmp[atj][1]*(_r6)
	   + cij     *(_r*ee + el*r2 + ec)
	   )*inside;
#endif
    double dp = (ljmp[atj][2]*(_r14)
               - ljmp[atj][3]*(_r8)
               + cij     *(_r3*edp-2.0*el)
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
void pairlistloopljcfe_zerodipole(const int (*jindex)[MAX_PAIR],
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
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(dynamic)
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
    if(ShortRange::alpha==0.0){
      if(operations.doVirialcalculation){
	pairlistloopljcfe1_zerodipole0(ShortRange::alpha,
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
	if(operations.doEnergycalculation){
	  pairlistloopljcfe1_zerodipole0(ShortRange::alpha,
					 jindex[i],
					 lj[i],npair[i],
					 posi, chargei, 
					 ljmp,
					 particlej.poscharge, 
					 cutoff2,
					 force[i],
					 e);
	}else{
	  pairlistloopljcf1_zerodipole0(ShortRange::alpha,
					jindex[i],
					lj[i],npair[i],
					posi, chargei, 
					ljmp,
					particlej.poscharge, 
					cutoff2,
					force[i] );
	}
      }
    }else{
      if(operations.doVirialcalculation){
	pairlistloopljcfe1_zerodipole(ShortRange::alpha,
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
	if(operations.doEnergycalculation){
	  pairlistloopljcfe1_zerodipole(ShortRange::alpha,
					jindex[i],
					lj[i],npair[i],
					posi, chargei, 
					ljmp,
					particlej.poscharge, 
					cutoff2,
					force[i],
					e);
	}else{
	  pairlistloopljcf1_zerodipole(ShortRange::alpha,
				       jindex[i],
				       lj[i],npair[i],
				       posi, chargei, 
				       ljmp,
				       particlej.poscharge, 
				       cutoff2,
				       force[i] );
	}
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

template<typename PA, typename GPA>
void pairlistloopljcfe_zerodipole(const int (*jindex)[MAX_PAIR],
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
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(dynamic)
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
    if(ShortRange::alpha==0.0){
      if(operations.doVirialcalculation){
	pairlistloopljcfe1_zerodipole0(ShortRange::alpha,
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
	if(operations.doEnergycalculation){
	  pairlistloopljcfe1_zerodipole0(ShortRange::alpha,
					 &(jindex[i][selfnpair[i]]),
					 &(lj[i][selfnpair[i]]),npair[i],
					 posi, chargei, 
					 ljmp,
					 particlej.poscharge, 
					 cutoff2,
					 force[i],
					 e);
	}else{
	  pairlistloopljcf1_zerodipole0(ShortRange::alpha,
					&(jindex[i][selfnpair[i]]),
					&(lj[i][selfnpair[i]]),npair[i],
					posi, chargei, 
					ljmp,
					particlej.poscharge, 
					cutoff2,
					force[i] );
	}
      }
    }else{
      if(operations.doVirialcalculation){
	pairlistloopljcfe1_zerodipole(ShortRange::alpha,
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
	if(operations.doEnergycalculation){
	  pairlistloopljcfe1_zerodipole(ShortRange::alpha,
					&(jindex[i][selfnpair[i]]),
					&(lj[i][selfnpair[i]]),npair[i],
					posi, chargei, 
					ljmp,
					particlej.poscharge, 
					cutoff2,
					force[i],
					e);
	}else{
	  pairlistloopljcf1_zerodipole(ShortRange::alpha,
				       &(jindex[i][selfnpair[i]]),
				       &(lj[i][selfnpair[i]]),npair[i],
				       posi, chargei, 
				       ljmp,
				       particlej.poscharge, 
				       cutoff2,
				       force[i] );
	}
      }
    }
  }
  energy = e;
  virial = v;
}



void pairlistloopljcf1_zerodipole_ljshift(const double alpha,
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
  //  double v=virial;
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
  double ee, edp;
  const double alpharc2 = alpha*alpha*cutoff2;
  const double alpharc = alpha*sqrt(cutoff2);
  const double ec2 = 0.5*M_2_SQRTPI*alpha*exp(-alpharc2);
  const double ec1 = 0.5*erfc(alpharc)/sqrt(cutoff2);
  const double el = (ec1+ec2)/cutoff2;
  const double ec = -3.0*ec1 - ec2;

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
    //    ee = ptc[iar][0]+dar*(ptc[iar][1]+dar*(ptc[iar][2]+dar*(ptc[iar][3]+dar*ptc[iar][4])));
    edp = ftc[iar][0]+dar*(ftc[iar][1]+dar*(ftc[iar][2]+dar*(ftc[iar][3]+dar*ftc[iar][4])));
#else
    ewaldrealforceandpottbl_(&ee,&edp,&ar);
#endif
#else
    ee = erfc(ar);
    edp = M_2_SQRTPI*ar*exp(-ar*ar)+ee;
#endif
    // e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
    //      - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
    // 	   + cij     *(_r*ee + el*r2 + ec)
    // 	   )*inside;
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3*edp-2.0*el)
                 )*inside;
    //    v += dp*r2;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  //  energy = e;
  //  virial = v;
}

void pairlistloopljcfe1_zerodipole_ljshift(const double alpha,
					   const int *jindex,
					   const int *lj,
					   const int npair,
					   const double posi[3],
					   const double chargei,
					   const double (*ljmp)[4],
					   const PosChargeArray& particlej,
					   const double cutoff2,
					   double force[3],
					   double &energy)
{
  int j;
  double fx=force[0], fy=force[1],fz=force[2];
  double e=energy;
  //  double v=virial;
  double &r6co0=ShortRange::r6co0, &r6co1=ShortRange::r6co1, &r6co2=ShortRange::r6co2;
  double &r12co0=ShortRange::r12co0, &r12co1=ShortRange::r12co1, &r12co2=ShortRange::r12co2;
  double &r8co1=ShortRange::r8co1, &r8co2=ShortRange::r8co2;
  double &r14co1=ShortRange::r14co1, &r14co2=ShortRange::r14co2;
  double ee, edp;
  const double alpharc2 = alpha*alpha*cutoff2;
  const double alpharc = alpha*sqrt(cutoff2);
  const double ec2 = 0.5*M_2_SQRTPI*alpha*exp(-alpharc2);
  const double ec1 = 0.5*erfc(alpharc)/sqrt(cutoff2);
  const double el = (ec1+ec2)/cutoff2;
  const double ec = -3.0*ec1 - ec2;

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
    double r3 = r2*r;
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
#ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
    double esingle = cij     *(_r*ee + el*r2 + ec);
    ljpair_energy += ljsingle*inside;
    zerodipole_energy += esingle*inside;
    e += ( ljsingle + esingle )*inside;
#else
    e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
         - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
	   + cij     *(_r*ee + el*r2 + ec)
	   )*inside;
#endif
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3*edp-2.0*el)
                 )*inside;
    //    v += dp*r2;
    fx += dx*dp;
    fy += dy*dp;
    fz += dz*dp;
  }
  force[0] = fx;
  force[1] = fy;
  force[2] = fz;
  energy = e;
  //  virial = v;
}

void pairlistloopljcfe1_zerodipole_ljshift(const double alpha,
					   const int *jindex,
					   const int *lj,
					   const int npair,
					   const double posi[3],
					   const double chargei,
					   const double (*ljmp)[4],
					   const PosChargeArray& particlej,
					   const double cutoff2,
					   double force[3],
					   double &energy,
					   double &virial)
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
  const double ec2 = 0.5*M_2_SQRTPI*alpha*exp(-alpharc2);
  const double ec1 = 0.5*erfc(alpharc)/sqrt(cutoff2);
  const double el = (ec1+ec2)/cutoff2;
  const double ec = -3.0*ec1 - ec2;

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
    double r3 = r2*r;
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
#ifdef CHECK_ENERGY
    double ljsingle = ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
      - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ));
    double esingle = cij     *(_r*ee + el*r2 + ec);
    ljpair_energy += ljsingle*inside;
    zerodipole_energy += esingle*inside;
    e += ( ljsingle + esingle )*inside;
#else
    e += ( ljmp[atj][0]*(_r12+(r3*(r12co1-r12co2*r)-r12co0))
         - ljmp[atj][1]*(_r6 +(r3*(r6co1 -r6co2 *r)-r6co0 ))
	   + cij     *(_r*ee + el*r2 + ec)
	   )*inside;
#endif
    double dp = (ljmp[atj][2]*(_r14-r*(r14co1-r14co2*r))
               - ljmp[atj][3]*(_r8 -r*(r8co1 -r8co2 *r))
               + cij     *(_r3*edp-2.0*el)
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
void pairlistloopljcfe_zerodipole_ljshift(const int (*jindex)[MAX_PAIR],
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
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(dynamic)
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
    if(ShortRange::alpha==0.0){
      if(operations.doVirialcalculation){
	pairlistloopljcfe1_zerodipole0_ljshift(ShortRange::alpha,
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
	if(operations.doEnergycalculation){
	  pairlistloopljcfe1_zerodipole0_ljshift(ShortRange::alpha,
						jindex[i],
						lj[i],npair[i],
						posi, chargei, 
						ljmp,
						particlej.poscharge, 
						cutoff2,
						force[i],
						e);
	}else{
	  pairlistloopljcf1_zerodipole0_ljshift(ShortRange::alpha,
					       jindex[i],
					       lj[i],npair[i],
					       posi, chargei, 
					       ljmp,
					       particlej.poscharge, 
					       cutoff2,
					       force[i] );
	}
      }
    }else{
      if(operations.doVirialcalculation){
	pairlistloopljcfe1_zerodipole_ljshift(ShortRange::alpha,
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
	if(operations.doEnergycalculation){
	  pairlistloopljcfe1_zerodipole_ljshift(ShortRange::alpha,
						jindex[i],
						lj[i],npair[i],
						posi, chargei, 
						ljmp,
						particlej.poscharge, 
						cutoff2,
						force[i],
						e);
	}else{
	  pairlistloopljcf1_zerodipole_ljshift(ShortRange::alpha,
					       jindex[i],
					       lj[i],npair[i],
					       posi, chargei, 
					       ljmp,
					       particlej.poscharge, 
					       cutoff2,
					       force[i] );
	}
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

template<typename PA, typename GPA>
void pairlistloopljcfe_zerodipole_ljshift(const int (*jindex)[MAX_PAIR],
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
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(dynamic)
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
    if(ShortRange::alpha==0.0){
      if(operations.doVirialcalculation){
	pairlistloopljcfe1_zerodipole0_ljshift(ShortRange::alpha,
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
	if(operations.doEnergycalculation){
	  pairlistloopljcfe1_zerodipole0_ljshift(ShortRange::alpha,
						&(jindex[i][selfnpair[i]]),
						&(lj[i][selfnpair[i]]),npair[i],
						posi, chargei, 
						ljmp,
						particlej.poscharge, 
						cutoff2,
						force[i],
						e);
	}else{
	  pairlistloopljcf1_zerodipole0_ljshift(ShortRange::alpha,
					       &(jindex[i][selfnpair[i]]),
					       &(lj[i][selfnpair[i]]),npair[i],
					       posi, chargei, 
					       ljmp,
					       particlej.poscharge, 
					       cutoff2,
					       force[i] );
	}
      }
    }else{
      if(operations.doVirialcalculation){
	pairlistloopljcfe1_zerodipole_ljshift(ShortRange::alpha,
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
	if(operations.doEnergycalculation){
	  pairlistloopljcfe1_zerodipole_ljshift(ShortRange::alpha,
						&(jindex[i][selfnpair[i]]),
						&(lj[i][selfnpair[i]]),npair[i],
						posi, chargei, 
						ljmp,
						particlej.poscharge, 
						cutoff2,
						force[i],
						e);
	}else{
	  pairlistloopljcf1_zerodipole_ljshift(ShortRange::alpha,
					       &(jindex[i][selfnpair[i]]),
					       &(lj[i][selfnpair[i]]),npair[i],
					       posi, chargei, 
					       ljmp,
					       particlej.poscharge, 
					       cutoff2,
					       force[i] );
	}
      }
    }
  }
  energy = e;
  virial = v;
}


template
void pairlistloopljcfe_zerodipole_ljshift(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe_zerodipole_ljshift(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe_zerodipole_ljshift(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe_zerodipole(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe_zerodipole(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe_zerodipole(const int (*jindex)[MAX_PAIR],
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



template<typename PA, typename GPA>
void pairlistloopljcfe_zerodipole0_ljshift(const int (*jindex)[MAX_PAIR],
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
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(dynamic)
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
      pairlistloopljcfe1_zerodipole0_ljshift(ShortRange::alpha,
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
      if(operations.doEnergycalculation){
	pairlistloopljcfe1_zerodipole0_ljshift(ShortRange::alpha,
					      jindex[i],
					      lj[i],npair[i],
					      posi, chargei, 
					      ljmp,
					      particlej.poscharge, 
					      cutoff2,
					      force[i],
					      e);
      }else{
	pairlistloopljcf1_zerodipole0_ljshift(ShortRange::alpha,
					     jindex[i],
					     lj[i],npair[i],
					     posi, chargei, 
					     ljmp,
					     particlej.poscharge, 
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

template<typename PA, typename GPA>
void pairlistloopljcfe_zerodipole0_ljshift(const int (*jindex)[MAX_PAIR],
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
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(dynamic)
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
      pairlistloopljcfe1_zerodipole0_ljshift(ShortRange::alpha,
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
      if(operations.doEnergycalculation){
	pairlistloopljcfe1_zerodipole0_ljshift(ShortRange::alpha,
					      &(jindex[i][selfnpair[i]]),
					      &(lj[i][selfnpair[i]]),npair[i],
					      posi, chargei, 
					      ljmp,
					      particlej.poscharge, 
					      cutoff2,
					      force[i],
					      e);
      }else{
	pairlistloopljcf1_zerodipole0_ljshift(ShortRange::alpha,
					     &(jindex[i][selfnpair[i]]),
					     &(lj[i][selfnpair[i]]),npair[i],
					     posi, chargei, 
					     ljmp,
					     particlej.poscharge, 
					     cutoff2,
					     force[i] );
      }
    }
  }
  energy = e;
  virial = v;
}


template
void pairlistloopljcfe_zerodipole0_ljshift(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe_zerodipole0_ljshift(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe_zerodipole0_ljshift(const int (*jindex)[MAX_PAIR],
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



template<typename PA, typename GPA>
void pairlistloopljcfe_zerodipole0(const int (*jindex)[MAX_PAIR],
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
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(dynamic)
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
      pairlistloopljcfe1_zerodipole0(ShortRange::alpha,
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
      if(operations.doEnergycalculation){
	pairlistloopljcfe1_zerodipole0(ShortRange::alpha,
				       jindex[i],
				       lj[i],npair[i],
				       posi, chargei, 
				       ljmp,
				       particlej.poscharge, 
				       cutoff2,
				       force[i],
				       e);
      }else{
	pairlistloopljcf1_zerodipole0(ShortRange::alpha,
				      jindex[i],
				      lj[i],npair[i],
				      posi, chargei, 
				      ljmp,
				      particlej.poscharge, 
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

template<typename PA, typename GPA>
void pairlistloopljcfe_zerodipole0(const int (*jindex)[MAX_PAIR],
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
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(runtime)
#  else
#   ifdef PL_OMP_STATIC
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(static)
#   else
#pragma omp parallel for private(i) reduction(+: e,v,ljpair_energy,zerodipole_energy) schedule(dynamic)
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
      pairlistloopljcfe1_zerodipole0(ShortRange::alpha,
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
      if(operations.doEnergycalculation){
	pairlistloopljcfe1_zerodipole0(ShortRange::alpha,
				       &(jindex[i][selfnpair[i]]),
				       &(lj[i][selfnpair[i]]),npair[i],
				       posi, chargei, 
				       ljmp,
				       particlej.poscharge, 
				       cutoff2,
				       force[i],
				       e);
      }else{
	pairlistloopljcf1_zerodipole0(ShortRange::alpha,
				      &(jindex[i][selfnpair[i]]),
				      &(lj[i][selfnpair[i]]),npair[i],
				      posi, chargei, 
				      ljmp,
				      particlej.poscharge, 
				      cutoff2,
				      force[i] );
      }
    }
  }
  energy = e;
  virial = v;
}


template
void pairlistloopljcfe_zerodipole0(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe_zerodipole0(const int (*jindex)[MAX_PAIR],
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
void pairlistloopljcfe_zerodipole0(const int (*jindex)[MAX_PAIR],
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
