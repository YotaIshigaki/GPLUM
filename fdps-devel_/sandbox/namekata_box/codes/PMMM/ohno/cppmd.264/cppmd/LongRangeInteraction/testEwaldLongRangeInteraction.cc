#include "EwaldInterface.h"
#include "EwaldLongRangeInteraction.h"

using namespace EwaldModule;

int main(int argc, char *argv[])
{
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif
  double tt=2.0;
  double w=10.0;
  LongRangeParameter param;
  param.cutoff = w/2;
  param.alpha = 1.0;
  param.surfaceDipole = 0;
  SpaceVector<int> gridSize(30, 30, 30);
  param.grid_num = gridSize;
  SpaceVector<double> box(w);
  param.boxSize = box;
  EwaldChargeAssign ca(0, param);
  EwaldPoissonSolver ps(0, param);
  GridData gridcharge(gridSize);
  GridData gridpotential(gridSize);
  EwaldLongRangeInteraction lr(0, param);

  ParticleArray particlearray(2);
  particlearray[0].position = Position(0.5*w+0.5*tt,0.0,0.0);
  particlearray[0].charge = 3.0;
  particlearray[1].position = Position(0.5*w-0.5*tt,0.0,0.0);
  particlearray[1].charge = -3.0;
  ParticleArray ghost;

  std::vector<ParticleRange> selfrange(1);
  std::vector<ParticleRange> ghostrange(1);
  selfrange[0].begin = 0;
  selfrange[0].end = particlearray.size();
  ghostrange[0].begin = 0;
  ghostrange[0].end = ghost.size();

  double energy = 0.0;
  ca.assign(particlearray, selfrange, ghost, ghostrange, gridcharge);
  ps.solvePoisson(gridcharge, gridpotential, energy);
  ca.backinterpolate(particlearray, selfrange,
                     ghost, ghostrange, gridpotential);
  std::cout << energy << std::endl;
  std::cout << particlearray[0].force << std::endl;
  std::cout << particlearray[1].force << std::endl;

  particlearray[0].force = 0.0;
  particlearray[1].force = 0.0;

  std::vector<TypeRange> tr(1);
  std::vector<TypeRange> ghosttr;
  tr[0].begin = 0;
  tr[0].end = particlearray.size();
  tr[0].lj.begin = tr[0].begin;
  tr[0].lj.end = tr[0].begin;
  tr[0].ljcoulomb.begin = tr[0].begin;
  tr[0].ljcoulomb.end = tr[0].begin;
  tr[0].coulomb.begin = tr[0].begin;
  tr[0].coulomb.end = tr[0].end;
  std::vector<int> self_longset_index(1,0);
  std::vector<int> ghost_longset_index;
  double potEnergy=0.0;
  lr.calcForce(particlearray, tr, self_longset_index,
               ghost, ghosttr, ghost_longset_index,
               potEnergy);
  std::cout << potEnergy << std::endl;
  std::cout << particlearray[0].force << std::endl;
  std::cout << particlearray[1].force << std::endl;
  return 0;
}
