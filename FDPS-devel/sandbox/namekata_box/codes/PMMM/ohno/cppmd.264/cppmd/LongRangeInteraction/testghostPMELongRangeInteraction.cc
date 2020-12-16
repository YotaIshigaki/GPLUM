#include "EwaldInterface.h"
#include "PMEInterface.h"
#include "PMELongRangeInteraction.h"
#include "EwaldLongRangeInteraction.h"
#include <cstring>
int main(int argc, char *argv[])
{
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif
  typedef PMEModule::PMEInterface PMEInterface;

  double tt=2.00;
  double w=64.0;
  LongRangeParameter param;
  param.cutoff = 9.0;
  double alpha = 0.0;
  double kCutoff = 0.0;
  EwaldModule::estimate_alpha_kCutoff(alpha,kCutoff,w,param.cutoff);
  std::cout << "default alpha, kCutoff " << alpha << " " << kCutoff << std::endl;
  std::cout << "estimate kCutoff for alpha=1.0 " 
            << EwaldModule::estimate_kCutoff(w,param.cutoff,1.0) << std::endl;
  std::cout << "estimate kCutoff for alpha=0.5 " 
            << EwaldModule::estimate_kCutoff(w,param.cutoff,0.5) << std::endl;
  param.order = 12;
  param.pmeType = PME;
  //  param.alpha = 1.0;
  param.alpha = alpha;
  param.kCutoff = kCutoff;
  param.surfaceDipole = 0;
  //  int gs = ceil(kCutoff*2.0);
  int gs = w;
  SpaceVector<int> gridSize(gs, gs, gs);
  param.grid_num = gridSize;
  SpaceVector<double> box(w);
  param.boxSize = box;
  /*
  PMEChargeAssign ca(0, param);
  PMEPoissonSolver ps(0, param);
  */
  LongRangeParameter sparam(param);
  sparam.pmeType = SmoothPME;
  PMEModule::PMEInterface* pmei = new PMEModule::PMEInterface(0,param);
  PMELongRangeInteraction lr(0, param, pmei);
  PMEModule::PMEInterface* pmeis = new PMEModule::PMEInterface(0,sparam);
  PMELongRangeInteraction slr(0, sparam, pmeis);
  EwaldModule::EwaldInterface *ei = new EwaldModule::EwaldInterface(0,param);
  EwaldLongRangeInteraction elr(0, param, ei);
  GridData gridcharge(gridSize);
  GridData gridpotential(gridSize);
  //  elr.chargeassign.setSide(w);
  double offset=0.0;

  enum Ptype{
    TwoAtom,
    Water
  };

  Ptype ptype=TwoAtom;
  
  if(argc>1){
    if((argv[1][0]=='W')||(argv[1][0]=='w')){
      ptype = Water;
    }
    if(argc>2){
      offset = strtod(argv[2],(char **)NULL);
    }
  }

  ParticleArray particlearray(1);
  ParticleArray ghost(1);

  if(ptype==Water){
    particlearray.resize(2);
    particlearray[0].position = Position(5.0+offset,5.0,5.0);
    particlearray[0].charge = -0.834;
    particlearray[1].position = Position(5.75695+offset,5.58588,5);
    particlearray[1].charge = 0.417;
    ghost[0].position = Position(4.24305+offset,5.58588,5);
    ghost[0].charge = 0.417;
  }else{
    particlearray[0].position = Position(0.5*w+0.5*tt+offset,0.0,0.0);
    particlearray[0].charge = 3.0;
    ghost[0].position = Position(0.5*w-0.5*tt+offset,0.0,0.0);
    ghost[0].charge = -3.0;
  }

  /*
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
  for(int i=0;i<particlearray.size();i++){
    std::cout << particlearray[i].force << std::endl;
    fsum += particlearray[i].force;
  }
  for(int i=0;i<ghost.size();i++){
    std::cout << ghost[i].force << std::endl;
    fsum += ghost[i].force;
  }
  std::cout << "Sum " << fsum << std::endl;
  */
  
  std::vector<TypeRange> tr(1);
  std::vector<TypeRange> ghosttr(1);
  tr[0].begin = 0;
  tr[0].end = particlearray.size();
  tr[0].lj.begin = tr[0].begin;
  tr[0].lj.end = tr[0].begin;
  tr[0].ljcoulomb.begin = tr[0].begin;
  tr[0].ljcoulomb.end = tr[0].begin;
  tr[0].coulomb.begin = tr[0].begin;
  tr[0].coulomb.end = tr[0].end;
  ghosttr[0].begin = 0;
  ghosttr[0].end = ghost.size();
  ghosttr[0].lj.begin = ghosttr[0].begin;
  ghosttr[0].lj.end = ghosttr[0].begin;
  ghosttr[0].ljcoulomb.begin = ghosttr[0].begin;
  ghosttr[0].ljcoulomb.end = ghosttr[0].begin;
  ghosttr[0].coulomb.begin = ghosttr[0].begin;
  ghosttr[0].coulomb.end = ghosttr[0].end;
  std::vector<int> self_longset_index(1,0);
  std::vector<int> ghost_longset_index(1,0);
  Force fsum;


  std::cout <<std::endl;
  std::cout << "use EwaldLongRangeInteraction" << std::endl;
  for(int i=0;i<particlearray.size();i++){
    particlearray[i].force = 0.0;
  }
  for(int i=0;i<ghost.size();i++){
    ghost[i].force = 0.0;
  }
  double EpotEnergy;
  elr.calcForce(particlearray, tr, self_longset_index,
               ghost, ghosttr, ghost_longset_index,
               EpotEnergy);
  std::cout << EpotEnergy << std::endl;
  fsum = Force(0.0,0.0,0.0);
  for(int i=0;i<particlearray.size();i++){
    std::cout << particlearray[i].force << std::endl;
    fsum += particlearray[i].force;
  }
  for(int i=0;i<ghost.size();i++){
    std::cout << ghost[i].force << std::endl;
    fsum += ghost[i].force;
  }
  std::cout << "Sum " << fsum << std::endl;

  std::vector<Force> ef(particlearray.size());
  std::vector<Force> egf(ghost.size());
  for(int i=0;i<ef.size();i++){
    ef[i] = particlearray[i].force;
  }
  for(int i=0;i<egf.size();i++){
    egf[i] = ghost[i].force;
  }

  std::cout <<std::endl;
  std::cout << "use PMELongRangeInteraction" << std::endl;
  for(int i=0;i<particlearray.size();i++){
    particlearray[i].force = 0.0;
  }
  for(int i=0;i<ghost.size();i++){
    ghost[i].force = 0.0;
  }
  double potEnergy;
  lr.calcForce(particlearray, tr, self_longset_index,
               ghost, ghosttr, ghost_longset_index,
               potEnergy);
  std::cout << potEnergy << " error " << potEnergy - EpotEnergy << " " << (potEnergy - EpotEnergy)/EpotEnergy << std::endl;
  fsum = Force(0.0,0.0,0.0);
  for(int i=0;i<particlearray.size();i++){
    std::cout << particlearray[i].force << " error " << particlearray[i].force - ef[i] << " " << (particlearray[i].force.x - ef[i].x)/ef[i].x << std::endl;
    fsum += particlearray[i].force;
  }
  for(int i=0;i<ghost.size();i++){
    std::cout << ghost[i].force  << " error " << ghost[i].force - egf[i] << " " << (ghost[i].force.x - egf[i].x)/egf[i].x <<  std::endl;
    fsum += ghost[i].force;
  }
  std::cout << "Sum " << fsum << std::endl;
  //  lr.pme_time_dump();

  std::cout <<std::endl;
  std::cout << "use SPMELongRangeInteraction" << std::endl;
  for(int i=0;i<particlearray.size();i++){
    particlearray[i].force = 0.0;
  }
  for(int i=0;i<ghost.size();i++){
    ghost[i].force = 0.0;
  }
  double SpotEnergy;
  slr.calcForce(particlearray, tr, self_longset_index,
                ghost, ghosttr, ghost_longset_index,
                SpotEnergy);
  std::cout << SpotEnergy <<  " error " << SpotEnergy - EpotEnergy << " " << (SpotEnergy - EpotEnergy)/EpotEnergy << std::endl;
  fsum = Force(0.0,0.0,0.0);
  for(int i=0;i<particlearray.size();i++){
    std::cout << particlearray[i].force << " error " << particlearray[i].force - ef[i] << " " << (particlearray[i].force.x - ef[i].x)/ef[i].x <<  std::endl;
    fsum += particlearray[i].force;
  }
  for(int i=0;i<ghost.size();i++){
    std::cout << ghost[i].force << " error " << ghost[i].force - egf[i] << " " << (ghost[i].force.x - egf[i].x)/egf[i].x << std::endl;
    fsum += ghost[i].force;
  }
  std::cout << "Sum " << fsum << std::endl;

  std::cout <<std::endl;
  std::cout << "Ewald " << std::endl;
  elr.lr_time_dump();
  std::cout <<std::endl;
  std::cout << "PME " << std::endl;
  lr.lr_time_dump();
  std::cout <<std::endl;
  std::cout << "SPME " << std::endl;
  slr.lr_time_dump();
  return 0;
}
