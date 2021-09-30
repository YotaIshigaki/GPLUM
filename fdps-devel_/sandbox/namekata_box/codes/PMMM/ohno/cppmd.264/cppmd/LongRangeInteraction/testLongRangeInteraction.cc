#include "ParticleInfo.h"
#include "RealSpaceLongRangeInteraction.h"
#include <cmath>
#include <cstdio>
using namespace std;

const SpaceVector<double> L(62.0, 62.0, 62.0);
const SpaceVector<int> grid_size(32, 32, 32);
const double bond = 0.957;
const double angle = 104.5;
const double theta = M_PI / 180.0 * angle / 2.0;

int main(int argc, char **argv) {
  ParticleArray particlearray(6);
  vector<ParticleRange> selfrange;
  ParticleArray ghost;
  vector<ParticleRange> ghostrange;
  GridData gridcharge(grid_size);
  GridData gridpotential(grid_size);
  double energy = -1.0;

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
  particlearray[4].position.x += bond_sin_theta;
  particlearray[4].position.z -= bond_cos_theta;
  particlearray[4].charge = 0.41;
  // water 2 -- Hydrogen
  particlearray[5].position = particlearray[3].position;
  particlearray[5].position.x -= bond_sin_theta;
  particlearray[5].position.z -= bond_cos_theta;
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

  RealSpaceChargeAssign ca(0);
  RealSpacePoissonSolver ps(0);

  ca.assign(particlearray, selfrange, ghost, ghostrange, gridcharge);
  ps.solvePoisson(gridcharge, gridpotential, energy);
  ca.backinterpolate(particlearray, selfrange, ghost, ghostrange, gridpotential);

  printf("#realspace pot=%.15e\n", energy);

  return 0;
}
