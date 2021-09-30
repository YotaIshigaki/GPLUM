#include <mpi.h>
#include "PME.h"
#include "ParticleInfo.h"
#include <cmath>
#include <cstdio>
#include <typeinfo>
#ifdef FFT_TIMER
#include "Timer.h"
#endif


#include <unistd.h>
extern char *optarg;
extern int optind, opterr, optopt;

using namespace std;

SpaceVector<double> L(62.0, 62.0, 62.0);
SpaceVector<int> grid_size(64, 64, 64);
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

template <typename PA, typename GPA>
void calc(const string& msg,
          const LongRangeParameter& param,
          PA& particlearray, vector<ParticleRange>& selfrange,
          GPA& ghost, vector<ParticleRange>& ghostrange,
	  const std::vector<ParticleRange>& self_selfenergy_range,
	  const std::vector<ParticleRange>& ghost_selfenergy_range,
	  int repeat)
{
  GridData gridcharge(param.grid_num);
  GridData gridpotential(param.grid_num);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //  cout << "calc rank " << rank<< endl;
  PMEModuleInterface *pmemoduleinterface = new PMEModuleInterface(rank, param, MPI_COMM_WORLD);
  PMELongRangeInteraction pmelongrange(rank, param, pmemoduleinterface);
  //  cout << "PMELongRangeInteraction construct rank " << rank<< endl;
  pmelongrange.initialize();
  //  cout << "PMELongRangeInteraction init rank " << rank<< endl;

  double energy =0.0;

  for (int r=0;r<repeat;r++){
    for (int i=0;i<particlearray.size();++i) {
      getforce(particlearray,i) = 0.0;
    }
    for (int i=0;i<ghost.size();++i) {
      getforce(ghost,i) = 0.0;
    }

    energy = 0.0;
    pmelongrange.chargeassign.assign(particlearray, selfrange, ghost, ghostrange, gridcharge, self_selfenergy_range, ghost_selfenergy_range);

    //  cout << " assign " << rank << endl;

    pmelongrange.poissonsolver.solvePoisson(gridcharge, gridpotential, energy);


    pmelongrange.chargeassign.backinterpolate(particlearray, selfrange, ghost, ghostrange,
					      gridpotential);
    pmelongrange.chargeassign.addenergy(particlearray, selfrange, ghost, ghostrange, energy);
  }


  double recvenergy;
  MPI_Allreduce(&energy, &recvenergy,
                1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  cout.precision(12);

  if (rank == 0) {
    cout << "energy " << recvenergy << endl;
  }

  MPI_Allreduce(&energy, &recvenergy,
		1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  std::vector<Force> force;
  std::vector<Force> recvforce;
  for (int i=0;i<particlearray.size();++i) {
    getforce(particlearray,i) += getforce(ghost,i);
    force.push_back(getforce(particlearray,i));
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
      cout << "force " << i << recvforce[i] << endl;
    }
  }
}

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  double alpha = 0.25;
  int surfaceDipole = 0;
  int direction = 0;
  int order = 12;
  int onlySmoothPME = 0;

  int repeat = 1;

  int opt;

  while ((opt = getopt(argc, argv, "r:a:s:d:L:n:g:o:e:S:")) != -1) {
    switch (opt) {
    case 'r':
      repeat = atoi(optarg);
      break;
    case 'a':
      alpha = atof(optarg);
      break;
    case 's':
      surfaceDipole = atoi(optarg);
      break;
    case 'd':
      direction = atoi(optarg);
      break;
    case 'L':
      L = SpaceVector<double>(atof(optarg));
      break;
    case 'n':
      nexternal = atoi(optarg);
      break;
    case 'g':
      grid_size = atoi(optarg);
      break;
    case 'o':
      order = atoi(optarg);
      break;
    case 'e':
      nexternal_ewald = atoi(optarg);
      break;
    case 'S':
      onlySmoothPME = atoi(optarg);
      break;
    default:
      printf("unknown\n");
    }
  }

  DebugLog::verbose = 1;
 
  CombinedParticleArray particlearray(6);
  vector<ParticleRange> selfrange;
  GhostParticleArray ghost;
  vector<ParticleRange> ghostrange;
  std::vector<ParticleRange> self_selfenergy_range;
  std::vector<ParticleRange> ghost_selfenergy_range;
  GridData gridcharge(grid_size);
  GridData gridpotential(grid_size);

  int rank;
  int nodenum;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nodenum);
  std::cout << "rank=" << rank << "/" << nodenum << std::endl;

  const double bond_sin_theta = bond * sin(theta);
  const double bond_cos_theta = bond * cos(theta);

  // water 1 -- Oxygen
  getpos(particlearray,0) = Position(L.x/2, L.y/2, L.z/2 + 1.0);
  getcharge(particlearray,0) = -0.82;
  // water 1 -- Hydrogen
  getpos(particlearray,1) = getpos(particlearray,0);
  getpos(particlearray,1).x += bond_sin_theta;
  getpos(particlearray,1).z += bond_cos_theta;
  getcharge(particlearray,1) = 0.41;
  // water 1 -- Hydrogen
  getpos(particlearray,2) =  getpos(particlearray,0);
  getpos(particlearray,2).x -= bond_sin_theta;
  getpos(particlearray,2).z += bond_cos_theta;
  getcharge(particlearray,2) = 0.41;

  // water 2 -- Oxygen
  getpos(particlearray,3) = Position(L.x/2, L.y/2, L.z/2 - 1.0);
  getcharge(particlearray,3) = -0.82;
  // water 2 -- Hydrogen
  getpos(particlearray,4) = getpos(particlearray,3);
  if (direction == 0) {
    getpos(particlearray,4).x += bond_sin_theta;
    getpos(particlearray,4).z -= bond_cos_theta;
  }
  else {
    getpos(particlearray,4).x -= bond_sin_theta;
    getpos(particlearray,4).z += bond_cos_theta;
  }
  getcharge(particlearray,4) = 0.41;
  // water 2 -- Hydrogen
  getpos(particlearray,5) = getpos(particlearray,3);
  if (direction == 0) {
    getpos(particlearray,5).x -= bond_sin_theta;
    getpos(particlearray,5).z -= bond_cos_theta;
  }
  else {
    getpos(particlearray,5).x += bond_sin_theta;
    getpos(particlearray,5).z += bond_cos_theta;
  }
  getcharge(particlearray,5) = 0.41;

  {
    ParticleRange pr;
    pr.begin = particlearray.size() * rank / nodenum;
    pr.end = particlearray.size() * (rank+1) / nodenum;
    selfrange.push_back(pr);
    self_selfenergy_range.push_back(pr);
  }
  {
    ghost.resize(particlearray.size());
    for(int i=0;i<particlearray.size();i++){
      getpos(ghost,i) = getpos(particlearray,i);
      getcharge(ghost,i) = getcharge(particlearray,i);
    }
    ParticleRange pr;
    pr.begin = 0;
    pr.end = particlearray.size() * rank / nodenum;
    ghostrange.push_back(pr);
    //    ghost_selfenergy_range.push_back(pr);
    pr.begin = particlearray.size() * (rank+1) / nodenum;
    pr.end = particlearray.size();
    ghostrange.push_back(pr);
    //    ghost_selfenergy_range.push_back(pr);
  }

  cout << "num cell " << selfrange.size() << std::endl;
  cout << "num gcell " << ghostrange.size() << std::endl;

  //--------------------------------------------------------------------

  LongRangeParameter param;
  param.boxSize = L;
  param.grid_num = grid_size;
  param.order = order;
  param.alpha = alpha;
  param.surfaceDipole = surfaceDipole;
  cout << "* SmoothPME *" << endl;
  param.pmeType = SmoothPME;
  calc("SmoothPME", param,
       particlearray, selfrange,
       ghost, ghostrange, 
       self_selfenergy_range, ghost_selfenergy_range,
       repeat);

  cout << "* SmoothPME end *" << endl;
  if (!onlySmoothPME) {
    cout << "* PME *" << endl;
    param.pmeType = PME;
    calc("PME", param,
	 particlearray, selfrange,
	 ghost, ghostrange, 
	 self_selfenergy_range, ghost_selfenergy_range,
	 repeat);
    cout << "* PME end *" << endl;
    cout << "* Ewald *" << endl;
  }

  MPI_Finalize();

#ifdef FFT_TIMER
  PerfCounter::print_time();
#endif
  return 0;
}
