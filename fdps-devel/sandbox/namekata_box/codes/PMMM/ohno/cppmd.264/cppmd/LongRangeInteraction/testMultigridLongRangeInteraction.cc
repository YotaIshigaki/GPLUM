#include "UseMPI.h"
#include "MultigridLongRangeInteraction.h"

void dump(const SpaceVector<int>& gridSize, const SpaceVector<double>& h,
          GridData& gridcharge, GridData& gridpotential)
{
  int j0 = 0;
  double erecm = 0.0;
  int nyz = gridSize.y * gridSize.z;
  int nz  = gridSize.z;
  for (int ix=0,jx=j0;ix<gridSize.x;++ix,jx+=nyz) {
    for (int iy=0,jy=jx;iy<gridSize.y;++iy,jy+=nz) {
      for (int iz=0,jz=jy;iz<gridSize.z;++iz,++jz) {
        if (std::abs(gridcharge.gridvalue[jz]) > 0.2) {      
          erecm += gridcharge.gridvalue[jz]*gridpotential.gridvalue[jz];
          std::cout << "M " << ix << " " << iy << " " << iz << " " << gridcharge.gridvalue[jz] << " " << gridpotential.gridvalue[jz] << " " << erecm << std::endl;
        }
      }
    }
  }
  erecm *= h[0]*h[1]*h[2]*0.5;
  std::cout << "ErecM " << erecm << std::endl;
}

int main(int argc, char *argv[])
{
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif

  double tt=2.0;
  double w=50.0;
  LongRangeParameter param;
  param.cutoff = 9;
  param.alpha = 1.0;
  //param.alpha = sqrt(0.5);
  double mgAlpha = param.alpha * sqrt(2.0);
#if 0
  SpaceVector<int> gridSize(64, 64, 64);
#else
  SpaceVector<int> gridSize(128);
#endif
  param.grid_num = gridSize;
  SpaceVector<double> box(w);
  SpaceVector<double> h(box[0]/gridSize[0],
                        box[1]/gridSize[1],
                        box[2]/gridSize[2]);
  param.boxSize = box;
  param.multigridType = VCycle;
  param.multigridIteration = 5;
  param.multigridFFTlevel = -1;
  param.node_geometry.size = SpaceVector<int>(1,1,1);
  param.node_geometry.have_cell = SpaceVector<int>(1,1,1);
  std::cout << "# MultigridChargeAssign & MultigridPoissonSolver " << std::endl;
  MultigridChargeAssign ca(0, param);
  MultigridPoissonSolver ps(0, param);
  MultigridLongRangeInteraction lr(0, param);
  GridData gridcharge(gridSize);
  GridData gridpotential(gridSize);

#ifdef MG_TEST_CHARGE
  ParticleArray particlearray(1);
  particlearray[0].position = Position(0.5*w+0.5*tt,0.5*w,0.5*w);
  particlearray[0].charge = 1.0;
#else
  ParticleArray particlearray(2);
  particlearray[0].position = Position(0.5*w+0.5*tt,0.5*w,0.5*w);
  particlearray[0].charge = 3.0;
  particlearray[1].position = Position(0.5*w-0.5*tt,0.5*w,0.5*w);
  particlearray[1].charge = -3.0;
#endif
  ParticleArray ghost;

  std::vector<Position> shift;
  int nn = 1;
  shift.push_back(Position());
  for (int ix=-nn;ix<=nn;++ix) {
    for (int iy=-nn;iy<=nn;++iy) {
      for (int iz=-nn;iz<=nn;++iz) {
        if (ix == 0 && iy == 0 && iz == 0) continue;
        shift.push_back(Position(ix*box[0], iy*box[1], iz*box[2]));
      }
    }
  }

  std::vector<ParticleRange> selfrange(1);
  std::vector<ParticleRange> ghostrange(1);
  selfrange[0].begin = 0;
  selfrange[0].end = particlearray.size();
  ghostrange[0].begin = 0;
  ghostrange[0].end = ghost.size();
  ca.assign(particlearray, selfrange, ghost, ghostrange, gridcharge);
#ifdef MG_TEST_CHARGE
  double tot = 0;
  for (int i=0;i<gridcharge.size;++i) {
    tot += gridcharge.gridvalue[i];
  }
  std::cout << "AAA " << tot*h[0]*h[1]*h[2]/(4.0*M_PI) << std::endl;
  assert(std::abs(tot*h[0]*h[1]*h[2]/(4.0*M_PI)-particlearray[0].charge)<1e-6);
#else
  double energy = 0.0;
  ps.solvePoisson(gridcharge, gridpotential, energy);
  ca.backinterpolate(particlearray, selfrange,
                     ghost, ghostrange, gridpotential);
  std::cout << "Erec " << energy << std::endl;
  std::cout << particlearray[0].force << std::endl;
  std::cout << particlearray[1].force << std::endl;
  std::cout << std::endl;
#ifdef MG_DUMP
  dump(gridSize, h, gridcharge, gridpotential);
#endif

  std::cout << "# MultigridLongRangeInteraction " << std::endl;
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
  double potEnergy;
  lr.calcForce(particlearray, tr, self_longset_index,
               ghost, ghosttr, ghost_longset_index,
               potEnergy);
  double ewaldrealpot = particlearray[0].charge * particlearray[1].charge * 
                        erfc(param.alpha*tt)/tt;
  double coulomb = particlearray[0].charge * particlearray[1].charge/tt;
  std::cout << "Erec+Eself  " << potEnergy << std::endl;
  std::cout << "Ereal       " << ewaldrealpot << std::endl;
  std::cout << particlearray[0].force << std::endl;
  std::cout << particlearray[1].force << std::endl;
  std::cout << "Erec+Eself+Ereal  " << potEnergy+ewaldrealpot << std::endl;
  std::cout << "Ecoulomb          " << coulomb << std::endl;
#endif
  return 0;
}
