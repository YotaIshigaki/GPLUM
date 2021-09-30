#include "UseMPI.h"
#include "MultigridLongRangeInteraction.h"
#include "ParticleInfo.h"
#include <cmath>
#include <cstdio>
#include <typeinfo>
using namespace std;

SpaceVector<double> L(62.0, 62.0, 62.0);
//SpaceVector<int> grid_size(32, 32, 32);
SpaceVector<int> grid_size(128, 128, 128);
const double bond = 0.957;
const double angle = 104.5;
const double theta = M_PI / 180.0 * angle / 2.0;
int nexternal = 2;
int nexternal_ewald = 1;

#define EWALD_REAL_FORCE(x) (M_2_SQRTPI*x*exp(-x*x) + erfc(x))
#define EWALD_REAL_POT(x) erfc(x)


void end(std::string msg)
{
  cout << msg << endl;
  MPI_Finalize();
  exit(1);
}

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
#ifdef USE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  TCA ca(rank, param);
  TPS ps(rank, param);
#else
  TCA ca(0, param);
  TPS ps(0, param);
#endif
  double energy = 0.0;
  ca.assign(particlearray, selfrange, ghost, ghostrange, gridcharge);
  ps.solvePoisson(gridcharge, gridpotential, energy);
#ifdef USE_MPI
  double recvenergy;
  MPI_Allreduce(&energy, &recvenergy,
                1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  cout.precision(12);
#ifdef USE_MPI
  if (rank == 0) {
    cout << "energy " << recvenergy << endl;
  }
#else
  cout << "energy " << energy << endl;
#endif
  ca.backinterpolate(particlearray, selfrange, ghost, ghostrange,
                     gridpotential);
  ca.addenergy(particlearray, selfrange, ghost, ghostrange, energy);
#ifdef USE_MPI
  MPI_Allreduce(&energy, &recvenergy,
                1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  std::vector<Force> force;
  std::vector<Force> recvforce;
  for (int i=0;i<particlearray.size();++i) {
    force.push_back(particlearray[i].force + ghost[i].force);
    recvforce.push_back(Force());
  }
  MPI_Allreduce(&energy, &recvenergy,
                1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&force[0], &recvforce[0],
                3*force.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0) {
    cout << "energy(+self) " << recvenergy << endl;
    cout.precision(6);
    for (int i=0;i<particlearray.size();++i) {
      cout << "wave " << i << recvforce[i] << endl;
    }
  }
#else
  cout << "energy(+self) " << energy << endl;
  cout.precision(6);
  for (int i=0;i<particlearray.size();++i) {
    cout << "wave " << i << particlearray[i].force << endl;
  }
#endif
}

int main(int argc, char **argv) {
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif

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
#ifdef USE_MPI
  int rank;
  int nodenum;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nodenum);
  std::cout << "rank=" << rank << "/" << nodenum << std::endl;
#endif

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

#ifdef USE_MPI
  ghost.assign(particlearray.begin(), particlearray.end());
  {
    ParticleRange pr;
    pr.begin = particlearray.size() * rank / nodenum;
    pr.end = particlearray.size() * (rank+1) / nodenum;
    selfrange.push_back(pr);
  }
  {
    ParticleRange pr;
    pr.begin = 0;
    pr.end = particlearray.size() * rank / nodenum;
    ghostrange.push_back(pr);
    pr.begin = particlearray.size() * (rank+1) / nodenum;
    pr.end = particlearray.size();
    ghostrange.push_back(pr);
  }
#else
  {
    ParticleRange pr;   pr.begin = 0; pr.end = 6;
    selfrange.push_back(pr);
  }
  {
    ParticleRange pr;   pr.begin = 0; pr.end = 0;
    ghostrange.push_back(pr);
  }
#endif

  //--------------------------------------------------------------------

  LongRangeParameter param;
  param.boxSize = L;
  param.grid_num = grid_size;
  param.alpha = alpha;
  param.cutoff = 9.0;
  param.multigridType = VCycle;
  param.multigridIteration = 5;
  param.multigridFFTlevel = -1;
#ifdef USE_MPI
  int nn;
  for (nn=1;nn*nn*nn < nodenum;nn*=2) {
  }
  if (nn*nn*nn != nodenum) {
    if (rank == 0) cout << "error: nodenum must be 8**n" << endl;
    end("A");
  }
  param.node_geometry.size = SpaceVector<int>(nn,nn,nn);
  param.node_geometry.full_cell = SpaceVector<int>(nn,nn,nn);
  param.node_geometry.have_cell = SpaceVector<int>(1,1,1);
#else
  param.node_geometry.size = SpaceVector<int>(1,1,1);
  param.node_geometry.have_cell = SpaceVector<int>(1,1,1);
#endif

  cout << "* Multigrid *" << endl;
  calc<MultigridChargeAssign, MultigridPoissonSolver>("Multigrid", param,
                                              particlearray, selfrange,
                                              ghost, ghostrange);
  cout << "* Multigrid end *" << endl;

#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
