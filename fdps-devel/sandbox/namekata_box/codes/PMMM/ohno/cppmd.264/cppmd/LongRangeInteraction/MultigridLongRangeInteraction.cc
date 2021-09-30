#include "UseMPI.h"
#include "MultigridLongRangeInteraction.h"
#include "LGM.h"
#include "MGPoisson.h"
#include "Array3D.h"
#include "Geometry.h"

using namespace MultigridModule;

namespace {
void convert(const ParticleArray& particles,
             const std::vector<ParticleRange>& range,
             std::vector<Position>& cd,
             std::vector<double>& charge)
{
  for (std::vector<ParticleRange>::size_type it=0;it < range.size();
       ++it) {
    for (int i = range[it].begin; i < range[it].end; ++i) { 
      const Particle& pi = particles[i];
      cd.push_back(pi.position);
      charge.push_back(pi.charge);
    }
  }
}
void convertForce(ParticleArray& particles,
                  const std::vector<ParticleRange>& range,
                  std::vector<Force>::iterator& itfc)
{
  for (std::vector<ParticleRange>::size_type it=0;it < range.size();
       ++it) {
    for (int i = range[it].begin; i < range[it].end; ++i,++itfc) { 
      Particle& pi = particles[i];
      pi.force = *itfc;
    }
  }
}

}

MultigridChargeAssign::MultigridChargeAssign(int unitid,
                                             const LongRangeParameter& _param)
  : param(_param), nodeid(unitid),
    caMethod(new LGM(param.alpha, 0.0, 0.0, 0.0, param.cutoff, param.cutoff)),
    node_geometry(param.node_geometry),
    p_rho(), p_phi(), cd(), charge(), fc(),
    selfEnergyCalculated(false), selfEnergy(0.0)
{
  SpaceVector<int> pos = node_geometry.getAbsolutePosition(nodeid);
  SpaceVector<int> grid_num(param.grid_num);
  std::vector<int> localIndex;
  localIndex.push_back(pos[0]);
  localIndex.push_back(pos[1]);
  localIndex.push_back(pos[2]);
  caMethod->initialize(param.boxSize, node_geometry.size,
                       grid_num, localIndex);
  p_rho = R3D::createArray(caMethod->getLocalNumGrid()[0],
                           caMethod->getLocalNumGrid()[1],
                           caMethod->getLocalNumGrid()[2]);
  p_phi = R3D::createArray(caMethod->getLocalNumGrid()[0],
                           caMethod->getLocalNumGrid()[1],
                           caMethod->getLocalNumGrid()[2]);
}

MultigridChargeAssign::~MultigridChargeAssign()
{
  delete(caMethod);
  delete(p_rho);
}

void MultigridChargeAssign::assign(ParticleArray& particlearray,
                                   const std::vector<ParticleRange>& selfrange,
                                   ParticleArray& ghost,
                                   const std::vector<ParticleRange>& ghostrange,
                                   GridData& gridcharge)
{
  R3D::Array& rho = *p_rho;
  for (int ix=0;ix<caMethod->getLocalNumGrid()[0];++ix) {
    for (int iy=0;iy<caMethod->getLocalNumGrid()[1];++iy) {
      for (int iz=0;iz<caMethod->getLocalNumGrid()[2];++iz) {
        rho[ix][iy][iz] = 0.0;
      }
    }
  }
  cd.clear();
  charge.clear();
  convert(particlearray, selfrange, cd, charge);
  convert(ghost, ghostrange, cd, charge);
  caMethod->chargeAssignment(cd, charge, rho, 0);
  double localcharge = 0.0;
  for (int ix=0;ix<caMethod->getLocalNumGrid()[0];++ix) {
    for (int iy=0;iy<caMethod->getLocalNumGrid()[1];++iy) {
      for (int iz=0;iz<caMethod->getLocalNumGrid()[2];++iz) {
        localcharge += rho[ix][iy][iz];
      }
    }
  }
#ifdef USE_MPI
  double totalcharge = 0.0;
  MPI_Allreduce(&localcharge, &totalcharge, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
#else
  double totalcharge = localcharge;
#endif
  double fac = totalcharge / gridcharge.size;
#ifdef MG_TEST_CHARGE
  fac = 0.0;
#endif

  int nyz = caMethod->getNy() * caMethod->getNz();
  int nz = caMethod->getNz();
  int j0 = (caMethod->getLocalNumGridBegin()[0]*nyz) +
           (caMethod->getLocalNumGridBegin()[1]*nz) +
           (caMethod->getLocalNumGridBegin()[2]);
  for (int ix=0,jx=j0;ix<caMethod->getLocalNumGrid()[0];++ix,jx+=nyz) {
    for (int iy=0,jy=jx;iy<caMethod->getLocalNumGrid()[1];++iy,jy+=nz) {
      for (int iz=0,jz=jy;iz<caMethod->getLocalNumGrid()[2];++iz,++jz) {
        gridcharge.gridvalue[jz] = rho[ix][iy][iz] = rho[ix][iy][iz]-fac;
      }
    }
  }
}

void MultigridChargeAssign::backinterpolate(
                              ParticleArray& particlearray,
                              const std::vector<ParticleRange>& selfrange,
                              ParticleArray& ghost,
                              const std::vector<ParticleRange>& ghostrange,
                              GridData& gridpotential)
{
  R3D::Array& rho = *p_rho;
  R3D::Array& phi = *p_phi;
  cd.clear();
  charge.clear();
  convert(particlearray, selfrange, cd, charge);
  convert(ghost, ghostrange, cd, charge);
  int nyz = caMethod->getNy() * caMethod->getNz();
  int nz = caMethod->getNz();
  int j0 = (caMethod->getLocalNumGridBegin()[0]*nyz) +
           (caMethod->getLocalNumGridBegin()[1]*nz) +
           (caMethod->getLocalNumGridBegin()[2]);
  for (int ix=0,jx=j0;ix<caMethod->getLocalNumGrid()[0];++ix,jx+=nyz) {
    for (int iy=0,jy=jx;iy<caMethod->getLocalNumGrid()[1];++iy,jy+=nz) {
      for (int iz=0,jz=jy;iz<caMethod->getLocalNumGrid()[2];++iz,++jz) {
        phi[ix][iy][iz] = gridpotential.gridvalue[jz];
      }
    }
  }
  fc.resize(cd.size());
  caMethod->backInterpolation(cd, charge, fc, rho, phi, 0);
  std::vector<Force>::iterator itfc = fc.begin();
  convertForce(particlearray, selfrange, itfc);
  convertForce(ghost, ghostrange, itfc);
}

void MultigridChargeAssign::addenergy(
                              ParticleArray& particlearray,
                              const std::vector<ParticleRange>& selfrange,
                              ParticleArray& ghost,
                              const std::vector<ParticleRange>& ghostrange,
                              double& energy)
{
  if (!selfEnergyCalculated) {
    charge.resize(particlearray.size());
    for (ParticleArray::size_type i=0;i<particlearray.size();++i) {
      charge[i] = particlearray[i].charge;
    }
    selfEnergy = caMethod->calcSelfPotential(charge);
    selfEnergyCalculated = true;
  }
#ifdef MG_DUMP
  std::cout << "selfEnergy " << selfEnergy << std::endl;
#endif
  energy += selfEnergy;
}

MultigridPoissonSolver::MultigridPoissonSolver(int unitid,
                                               const LongRangeParameter& _param)
  : param(_param), nodeid(unitid), poisson(new MGPoisson()),
    node_geometry(param.node_geometry),
    pos(node_geometry.getAbsolutePosition(nodeid)),
    localNumGrid(param.grid_num[0]/node_geometry.size[0],
                 param.grid_num[1]/node_geometry.size[1],
                 param.grid_num[2]/node_geometry.size[2]),
    localSide(param.boxSize[0]/node_geometry.size[0],
              param.boxSize[1]/node_geometry.size[1],
              param.boxSize[2]/node_geometry.size[2]),
    p_rho(R3D::createArray(localNumGrid[0], localNumGrid[1], localNumGrid[2])),
    p_phi(R3D::createArray(localNumGrid[0], localNumGrid[1], localNumGrid[2]))
{
#ifdef USE_MPI
  int neighborRanks[3][3][3];
  for (int ix=0;ix<3;++ix) {
    for (int iy=0;iy<3;++iy) {
      for (int iz=0;iz<3;++iz) {
        SpaceVector<int> rpos(ix-1, iy-1, iz-1);
        neighborRanks[ix][iy][iz] = 
            node_geometry.getNodeIDofRelativePosition(rpos, nodeid);
      }
    }
  }
  poisson->Initialize(MPI_COMM_WORLD, &neighborRanks[0][0][0],
                      &localNumGrid[0], &localSide[0],
                      p_phi->data(), p_rho->data(),
                      param.multigridFFTlevel);
#else
  poisson->Initialize(&localNumGrid[0], &localSide[0],
                      p_phi->data(), p_rho->data(),
                      param.multigridFFTlevel);
#endif
}

MultigridPoissonSolver::~MultigridPoissonSolver()
{
  delete(poisson);
}

void MultigridPoissonSolver::solvePoisson(GridData& gridcharge,
                                          GridData& gridpotential,
                                          double& energy)
{
  R3D::Array& rho = *p_rho;
  R3D::Array& phi = *p_phi;
  int nyz = param.grid_num[1] * param.grid_num[2];
  int nz = param.grid_num[2];
  int j0 = (pos[0]*localNumGrid[0]*nyz) +
           (pos[1]*localNumGrid[1]*nz) +
           (pos[2]*localNumGrid[2]);
  SpaceVector<double> localSide(param.boxSize[0]/node_geometry.size[0],
                                param.boxSize[1]/node_geometry.size[1],
                                param.boxSize[2]/node_geometry.size[2]);
  const double FOUR_PI = M_PI*4.0;
  for (int ix=0,jx=j0;ix<localNumGrid[0];++ix,jx+=nyz) {
    for (int iy=0,jy=jx;iy<localNumGrid[1];++iy,jy+=nz) {
      for (int iz=0,jz=jy;iz<localNumGrid[2];++iz,++jz) {
        rho[ix][iy][iz] = gridcharge.gridvalue[jz] * FOUR_PI;
        phi[ix][iy][iz] = 0.0;
      }
    }
  }
  switch (param.multigridType) {
  case SOR:
    poisson->SORNIter(param.multigridIteration);
    break;
  case VCycle:
    poisson->VCycleNIter(param.multigridIteration);
    break;
  case FMG:
    poisson->FMG(param.multigridIteration);
    break;
  }
  energy = 0.0;
  double erecm = 0.0;
  for (int ix=0,jx=j0;ix<localNumGrid[0];++ix,jx+=nyz) {
    for (int iy=0,jy=jx;iy<localNumGrid[1];++iy,jy+=nz) {
      for (int iz=0,jz=jy;iz<localNumGrid[2];++iz,++jz) {
        gridpotential.gridvalue[jz] = phi[ix][iy][iz];
        energy += gridcharge.gridvalue[jz]*gridpotential.gridvalue[jz];
#ifdef MG_DUMP
if (std::abs(gridcharge.gridvalue[jz]) > 0.2) {      
  erecm += gridcharge.gridvalue[jz]*gridpotential.gridvalue[jz];
  std::cout << "M " << ix << " " << iy << " " << iz << " " << gridcharge.gridvalue[jz] << " " << gridpotential.gridvalue[jz] << " " << erecm << std::endl;
}
#endif
      }
    }
  }
  SpaceVector<double> h(param.boxSize[0]/param.grid_num[0],
                        param.boxSize[1]/param.grid_num[1],
                        param.boxSize[2]/param.grid_num[2]);
  energy *= h[0]*h[1]*h[2]*0.5;
  erecm *= h[0]*h[1]*h[2]*0.5;
#ifdef MG_DUMP
  std::cout << "Erec  " << energy << std::endl;
  std::cout << "ErecM " << erecm << std::endl;
#endif
}
