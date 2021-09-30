//#include "UseMPI.h"
#include "PMEInterface.h"
#include "PMELongRangeInteraction.h"

int main(int argc, char *argv[])
{
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif
  typedef PMEModule::PMEInterface PMEInterface;

  double tt=2.0;
  double w=32.0;
  LongRangeParameter param;
  param.cutoff = 9.0;
  param.order = 12;
  //param.pmeType = PME;
  param.pmeType = SmoothPME;
  param.alpha = 1.0;
  param.surfaceDipole = 0;
  int gs = w*2;
  SpaceVector<int> gridSize(gs, gs, gs);
  param.grid_num = gridSize;
  SpaceVector<double> box(w);
  param.boxSize = box;
  PMEChargeAssign ca(0, param);
  PMEPoissonSolver ps(0, param);
  PMELongRangeInteraction lr(0, param);
  GridData gridcharge(gridSize);
  GridData gridpotential(gridSize);

#if 0
  ParticleArray particlearray(2);
  particlearray[0].position = Position(0.5*w+0.5*tt,0.0,0.0);
  particlearray[0].charge = 3.0;
  particlearray[1].position = Position(0.5*w-0.5*tt,0.0,0.0);
  particlearray[1].charge = -3.0;
#else  /// 1 Water
  ParticleArray particlearray(3);
  particlearray[0].position = Position(5,5,5);
  particlearray[1].position = Position(5.75695,5.58588,5);
  particlearray[2].position = Position(4.24305,5.58588,5);
  particlearray[0].charge = -0.834;
  particlearray[1].charge = 0.417;
  particlearray[2].charge = 0.417;
#endif
    
  ParticleArray ghost;

  std::vector<ParticleRange> selfrange(1);
  std::vector<ParticleRange> ghostrange(1);
  selfrange[0].begin = 0;
  selfrange[0].end = particlearray.size();
  ghostrange[0].begin = 0;
  ghostrange[0].end = ghost.size();
  ca.assign(particlearray, selfrange, ghost, ghostrange, gridcharge);
  double energy = 0.0;
  ps.solvePoisson(gridcharge, gridpotential, energy);
  ca.backinterpolate(particlearray, selfrange,
                     ghost, ghostrange, gridpotential);
  std::cout << energy << std::endl;
  Force fsum;
  for(int i=0;i<particlearray.size();i++){
    std::cout << particlearray[i].force << std::endl;
    fsum += particlearray[i].force;
  }
  std::cout << "Sum " << fsum << std::endl;


  std::cout << "use LongRangeInteraction" << std::endl;
  for(int i=0;i<particlearray.size();i++){
    particlearray[i].force = Force(0,0,0);
  }

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
  double potEnergy;
  lr.calcForce(particlearray, tr, self_longset_index,
               ghost, ghosttr, ghost_longset_index,
               potEnergy);
  std::cout << potEnergy << std::endl;
  fsum = Force(0.0,0.0,0.0);
  for(int i=0;i<particlearray.size();i++){
    std::cout << particlearray[i].force << std::endl;
    fsum += particlearray[i].force;
  }
  std::cout << "Sum " << fsum << std::endl;
  return 0;
}
