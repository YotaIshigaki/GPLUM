#include "PMEInterface.h"
#include "PMELongRangeInteraction.h"
#include "EwaldInterface.h"
#include "EwaldLongRangeInteraction.h"
#include "ParticleInfo.h"
#include <cmath>
#include <cstdio>
#include <typeinfo>
using namespace std;

SpaceVector<double> L(62.0, 62.0, 62.0);
SpaceVector<int> grid_size(32, 32, 32);
const double bond = 0.957;
const double angle = 104.5;
const double theta = M_PI / 180.0 * angle / 2.0;
int nexternal = 2;
int nexternal_ewald = 1;

#define EWALD_REAL_FORCE(x) (M_2_SQRTPI*x*exp(-x*x) + erfc(x))
#define EWALD_REAL_POT(x) erfc(x)

void calcpv(std::vector<Position>& pv, int n) 
{
  pv.clear();
  for (int nx=-n;nx<=n;++nx) {
    for (int ny=-n;ny<=n;++ny) {
      for (int nz=-n;nz<=n;++nz) {
        if (nx == 0 && ny == 0 && nz == 0) continue;
        pv.push_back(Position(L.x*nx, L.y*ny, L.z*nz));
      }
    }
  }
}

template <typename TCA, typename TPS>
void calc(const string& msg,
          const LongRangeParameter& param,
          ParticleArray& particlearray, vector<ParticleRange>& selfrange,
          ParticleArray& ghost, vector<ParticleRange>& ghostrange)
{
  GridData gridcharge(param.grid_num);
  GridData gridpotential(param.grid_num);

  for (int i=0;i<particlearray.size();++i) {
    particlearray[i].force = 0.0;
  }
  TCA ca(0, param);
  TPS ps(0, param);
  double energy = 0.0;
  ca.assign(particlearray, selfrange, ghost, ghostrange, gridcharge);
  ps.solvePoisson(gridcharge, gridpotential, energy);
  cout.precision(12);
  cout << "energy " << energy << endl;
  ca.backinterpolate(particlearray, selfrange, ghost, ghostrange,
                     gridpotential);
  ca.addenergy(particlearray, selfrange, ghost, ghostrange, energy);
  cout << "energy(+self) " << energy << endl;
  cout.precision(6);
  for (int i=0;i<particlearray.size();++i) {
    cout << "wave " << i << particlearray[i].force << endl;
  }

  std::vector<Position> pv;
  calcpv(pv, nexternal_ewald);

  double ewaldrealpot = 0.0;
  double outer_ewaldrealpot = 0.0;
  for (int i=0;i<particlearray.size();++i) {
    Force f;
    Force f0;
    for (int j=0;j<particlearray.size();++j) {
      double qq = particlearray[i].charge*particlearray[j].charge;
      Position d = particlearray[i].position-particlearray[j].position;
      for (std::vector<Position>::const_iterator it=pv.begin();it != pv.end();
           ++it) {
        double r = (d+(*it)).norm();
        double dp = qq/(r*r*r);
        f += d*dp*EWALD_REAL_FORCE(param.alpha*r);
        outer_ewaldrealpot += qq/r*EWALD_REAL_POT(param.alpha*r);
      }
      if (i != j) {
        double r = d.norm();
        double dp = qq/(r*r*r);
        f += d*dp*EWALD_REAL_FORCE(param.alpha*r);
        f0 += d*dp;
        ewaldrealpot += qq/r*EWALD_REAL_POT(param.alpha*r);
        //cout << r << " " << EWALD_REAL_FORCE(param.alpha*r) << f << f0 << endl;
      }
    }
    //cout << i << f0 << endl;
    //cout << i << f << endl;
    cout << "real+wave " << i << f+particlearray[i].force << endl;
  }
  ewaldrealpot *= 0.5;
  cout << "energy(outerreal) " << outer_ewaldrealpot << endl;
  ewaldrealpot += outer_ewaldrealpot;
  cout.precision(12);
  cout << "energy(real+wave) " << ewaldrealpot+energy << endl;
  cout.precision(6);

  calcpv(pv, nexternal);
  //cout << pv.size() << endl;
  double innerpot = 0.0;
  double outerpot = 0.0;
  for (int i=0;i<particlearray.size();++i) {
    Force f;
    for (int j=0;j<particlearray.size();++j) {
      Position d = particlearray[i].position-particlearray[j].position;
      double qq = particlearray[i].charge*particlearray[j].charge;
      for (std::vector<Position>::const_iterator it=pv.begin();it != pv.end();
           ++it) {
        double r = (d+(*it)).norm();
        double dp = qq/(r*r*r);
        f += d*dp;
        outerpot += qq/r;
      }
    }
    //cout << i << f << endl;
    for (int j=0;j<particlearray.size();++j) {
      Position d = particlearray[i].position-particlearray[j].position;
      double qq = particlearray[i].charge*particlearray[j].charge;
      if (i != j) {
        double r = d.norm();
        double dp = qq/(r*r*r);
        f += d*dp;
        innerpot += qq/r;
      }
    }
    cout << "direct " << i << f << endl;
  }
  innerpot *= 0.5;
  cout.precision(12);
  cout << "energy(outer) " << outerpot << endl;
  cout << "energy(direct) " << innerpot+outerpot << endl;
  cout << msg << " denergy(real+wave-direct) "
       << (ewaldrealpot+energy)-(innerpot+outerpot) << endl;
  cout.precision(6);
}

int main(int argc, char **argv) {
  double alpha = 0.25;
  int surfaceDipole = 0;
  int direction = 0;
  int order = 12;
  int onlySmoothPME = 0;
  if (argc > 1) {
    alpha = atof(argv[1]);
  }
  if (argc > 2) {
    surfaceDipole = atoi(argv[2]);
  }
  if (argc > 3) {
    direction = atoi(argv[3]);
  }
  if (argc > 4) {
    L = SpaceVector<double>(atof(argv[4]));
  }
  if (argc > 5) {
    nexternal = atoi(argv[5]);
  }
  if (argc > 6) {
    grid_size = atoi(argv[6]);
  }
  if (argc > 7) {
    order = atoi(argv[7]);
  }
  if (argc > 8) {
    nexternal_ewald = atoi(argv[8]);
  }
  if (argc > 9) {
    onlySmoothPME = atoi(argv[9]);
  }
  ParticleArray particlearray(6);
  vector<ParticleRange> selfrange;
  ParticleArray ghost;
  vector<ParticleRange> ghostrange;
  GridData gridcharge(grid_size);
  GridData gridpotential(grid_size);

  const double bond_sin_theta = bond * sin(theta);
  const double bond_cos_theta = bond * cos(theta);

  // water 1 -- Oxygen
  particlearray[0].position = Position(L.x/2, L.y/2, L.z/2 + 1.0);
  particlearray[0].charge = -0.82;
  // water 1 -- Hydrogen
  particlearray[1].position = particlearray[0].position;
  particlearray[1].position.x += bond_sin_theta;
  particlearray[1].position.z += bond_cos_theta;
  particlearray[1].charge = 0.41;
  // water 1 -- Hydrogen
  particlearray[2].position = particlearray[0].position;
  particlearray[2].position.x -= bond_sin_theta;
  particlearray[2].position.z += bond_cos_theta;
  particlearray[2].charge = 0.41;

  // water 2 -- Oxygen
  particlearray[3].position = Position(L.x/2, L.y/2, L.z/2 - 1.0);
  particlearray[3].charge = -0.82;
  // water 2 -- Hydrogen
  particlearray[4].position = particlearray[3].position;
  if (direction == 0) {
    particlearray[4].position.x += bond_sin_theta;
    particlearray[4].position.z -= bond_cos_theta;
  }
  else {
    particlearray[4].position.x -= bond_sin_theta;
    particlearray[4].position.z += bond_cos_theta;
  }
  particlearray[4].charge = 0.41;
  // water 2 -- Hydrogen
  particlearray[5].position = particlearray[3].position;
  if (direction == 0) {
    particlearray[5].position.x -= bond_sin_theta;
    particlearray[5].position.z -= bond_cos_theta;
  }
  else {
    particlearray[5].position.x += bond_sin_theta;
    particlearray[5].position.z += bond_cos_theta;
  }
  particlearray[5].charge = 0.41;

  {
    ParticleRange pr;   pr.begin = 0; pr.end = 6;
    selfrange.push_back(pr);
  }
  {
    ParticleRange pr;   pr.begin = 0; pr.end = 0;
    ghostrange.push_back(pr);
  }

  //--------------------------------------------------------------------

  LongRangeParameter param;
  param.boxSize = L;
  param.grid_num = grid_size;
  param.order = order;
  param.alpha = alpha;
  param.surfaceDipole = surfaceDipole;
  cout << "* SmoothPME *" << endl;
  param.pmeType = SmoothPME;
  calc<PMEChargeAssign, PMEPoissonSolver>("SmoothPME", param,
                                          particlearray, selfrange,
                                          ghost, ghostrange);
  cout << "* SmoothPME end *" << endl;
  if (onlySmoothPME) return 0;
  cout << "* PME *" << endl;
  param.pmeType = PME;
  calc<PMEChargeAssign, PMEPoissonSolver>("PME", param,
                                          particlearray, selfrange,
                                          ghost, ghostrange);
  cout << "* PME end *" << endl;
  cout << "* Ewald *" << endl;
  calc<EwaldChargeAssign, EwaldPoissonSolver>("Ewald", param,
                                              particlearray, selfrange,
                                              ghost, ghostrange);
  cout << "* Ewald end *" << endl;

  return 0;
}
