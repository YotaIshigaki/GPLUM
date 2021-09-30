#include <iostream>
#include <boost/timer.hpp>
#include "LJAmber94.h"
#ifndef RDTSCTIME
#define RDTSCTIME
#endif
#include "ShortRangeInteraction.h"

//using namespace std;

template<typename TP>
void zero_force(TP& p, int n)
{
  for(int j=0;j<n;j++){
    p[j].force(0.0,0.0,0.0);
  }
}

void setparticlearray(long n, ParticleParameters1& p)
{
  p.resize(n);
  for(long i=0;i<n;i++){
    p[i].position = (0.0,0.0,0.0);
    p[i].force = (0.0,0.0,0.0);
    p[i].charge = 0.0;
    p[i].atomtype = 0;
  }
}

template<typename TPSRC, typename TPDST>
void copyparticlearray(TPSRC& psrc, TPDST& pdst)
{
  long n = psrc.size();
  pdst.resize(n);
  for(long i=0;i<n;i++){
    pdst[i].position = psrc[i].position;
    pdst[i].force = psrc[i].force;
    pdst[i].charge = psrc[i].charge;
    pdst[i].atomtype = psrc[i].atomtype;
  }
}

void setijsets(long ni, long nj, PairRangeList& isets, std::vector<PairRangeList>& jsets)
{
  isets.resize(1);
  isets[0].begin=0;
  isets[0].end=ni;
  jsets.resize(1);
  jsets[0].resize(1);
  jsets[0][0].begin=0;
  jsets[0][0].end=nj;

}

int
main(int argc, char **argv)
{
  boost::timer tm;
  double time;

  LJAmber94 ljamber94;
  PairRangeList coulombisets;
  std::vector<PairRangeList> coulombjsets;
  PairRangeList ljisets;
  std::vector<PairRangeList> ljjsets;
  ParticleParameters1 particlei;
  ParticleParameters1 particlej;
  ParticleArray particlearrayi;
  ParticleArray particlearrayj;
  double energy;
  int n = 120;
  long ni=10, nj=120;
  int rl = 1000;
  int l;
  bool useEwald=false;
  double cutoff2 = 101.0;
  int testtype=0;
  int currenttest;

  int opt;
  while ((opt = getopt(argc, argv, "n:l:ec:t:"))!=-1){
    switch (opt) {
    case 'n':
      n = strtol(optarg,(char **)NULL,0);
      break;
    case 'l':
      rl = strtol(optarg,(char **)NULL,0);
      break;
    case 'e':
      useEwald = true;
      break;
    case 't':
      testtype = strtol(optarg,(char **)NULL,0);
      break;
    case 'c':
      double cutoff = strtod(optarg,(char **)NULL);
      cutoff2 = cutoff*cutoff;
      break;
    }
  }
  l = n*rl;
  ni = n;
  nj = n;
  
  LJ lj1pair;

  ShortRangeInteraction<LJEwaldReal, ParticleParameter1, ParticleParameter1> ljewald;
  ShortRangeInteraction<LJCoulomb, ParticleParameter1, ParticleParameter1> ljc;
  ShortRangeInteraction<Coulomb, ParticleParameter1, ParticleParameter1> cl;
  ShortRangeInteraction<EwaldReal, ParticleParameter1, ParticleParameter1> ewald;
  ShortRangeInteraction<LJ, ParticleParameter1, ParticleParameter1> lj;
  ShortRangeInteraction<LJ, Particle, Particle> lj0;



  ShortRangeInteractions<EwaldReal, LJEwaldReal, ParticleParameter1, ParticleParameter1> ewald_and_lje;
  ShortRangeInteractions<Coulomb, LJCoulomb, ParticleParameter1, ParticleParameter1> cl_and_ljc;

  EwaldAndLJEInteractions1 ewald_comb;
  EwaldAndLJEInteractions ewald_comb0;

  ljamber94.convertLJMixparameterArray(ShortRange::ljmixparameters);
  std::cout << ShortRange::ljmixparameters[0][0].potFirst << std::endl;
  std::cout << "Ewald alpha " << ShortRange::alpha <<  std::endl;
  std::cout << "Cutoff^2 " << cutoff2 << std::endl;


  setparticlearray(ni,particlei);
  for(long i=0;i<ni;i++){
    particlei[i].position.x = 3.0 + 7.0*((double)i/(double)ni);
  }
  setparticlearray(nj,particlej);

  copyparticlearray(particlei,particlearrayi);
  copyparticlearray(particlej,particlearrayj);

  /// set coulomb i,j-particlearray as all of i,j-particlearray
  setijsets(ni,nj,coulombisets,coulombjsets);
  /// change coulomb j-partilces as 2nd half of j-particlearray
  coulombjsets[0][0].begin=(nj>>1);
  /// set LJ i,j-particlearray as all of i,j-particlearray
  setijsets(ni,nj,ljisets,ljjsets);
  /// change LJ j-partilces as 1st half of j-particlearray
  ljjsets[0][0].end=(nj>>1);
  energy = 0.0;

  
  std::cout << "nj " << nj << " ni " << ni << " repeat " << rl << " total " << ni*nj*rl << " for each lj and coulomb" << std::endl;
  std::cout << "LJ i " << ljisets[0].end-ljisets[0].begin << " LJ j " << ljjsets[0][0].end-ljjsets[0][0].begin << " CL i " << coulombisets[0].end-coulombisets[0].begin <<  " CL j " << coulombjsets[0][0].end-coulombjsets[0][0].begin << std::endl;
  std::cout << "LJ " << (ljisets[0].end-ljisets[0].begin)*(ljjsets[0][0].end-ljjsets[0][0].begin) << std::endl;
  std::cout << "Coulomb " << (ljisets[0].end-ljisets[0].begin)*(coulombjsets[0][0].end-coulombjsets[0][0].begin) + (coulombisets[0].end-coulombisets[0].begin)*(ljjsets[0][0].end-ljjsets[0][0].begin + ljjsets[0][0].end-ljjsets[0][0].begin) << std::endl;

  if((testtype==0)||(testtype==1)){
    tm.restart();  
    for(int r=0;r<rl;r++){
      ewald_and_lje.loopijsets(ljisets,ljjsets,coulombisets,coulombjsets,particlei,particlej,energy);
    }
    time = tm.elapsed();
    std::cout << "ewald_and_lje  " << time << std::endl;
    std::cout << energy << std::endl;
  }

  if((testtype==0)||(testtype==2)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      ewald_comb.loopijsets(ljisets,ljjsets,coulombisets,coulombjsets,particlei,particlej,energy);
    }
    time = tm.elapsed();
    std::cout << "ewald_comb " << time << std::endl;
    std::cout << energy << std::endl;
  }

  if((testtype==0)||(testtype==3)){
    zero_force(particlei,ni);
  zero_force(particlej,nj);
  energy = 0.0;
  tm.restart();  
  for(int r=0;r<rl;r++){
    cl_and_ljc.loopijsets(ljisets,ljjsets,coulombisets,coulombjsets,particlei,particlej,energy);
  }
  time = tm.elapsed();
  std::cout << "cl_and_ljc " << time << std::endl;
  std::cout << energy << std::endl;
  }

  std::cout << "nj " << (nj>>1) << " ni " << (ni>>1) << " repeat " << rl << " total " << ni*nj*rl << std::endl;

  ljisets[0].end=(ni>>1);
  ljjsets[0][0].end=(nj>>1);

  if((testtype==0)||(testtype==4)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      cl.loopijsets(ljisets,ljjsets,particlei,particlej,energy);
    }
    time = tm.elapsed();
    std::cout << "cl " << time << std::endl;
    std::cout << energy << std::endl;
  }

  ljisets[0].end=ni;

  std::cout << "nj " << nj << " ni " << ni << " repeat " << rl << " total " << ni*nj*rl << std::endl;

  coulombjsets[0][0].begin=nj;
  ljjsets[0][0].end=nj;

  if((testtype==0)||(testtype==5)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      cl.loopijsets(ljisets,ljjsets,particlei,particlej,energy);
    }
    time = tm.elapsed();
    std::cout << "cl " << time << std::endl;
    std::cout << energy << std::endl;
  }

  if((testtype==0)||(testtype==6)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      ewald.loopijsets(ljisets,ljjsets,particlei,particlej,energy);
    }
    time = tm.elapsed();
    std::cout << "ewald " << time << std::endl;
    std::cout << energy << std::endl;
  }

  coulombjsets[0][0].begin=nj;
  ljjsets[0][0].end=nj;
  
  if((testtype==0)||(testtype==7)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      lj.loopijsets(ljisets,ljjsets,particlei,particlej,energy);
    }
    time = tm.elapsed();
    std::cout << "lj " << time << std::endl;
    std::cout << energy << std::endl;
  }

  if((testtype==0)||(testtype==8)){
    zero_force(particlearrayi,ni);
    zero_force(particlearrayj,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      lj0.loopijsets(ljisets,ljjsets,particlearrayi,particlearrayj,energy);
    }
    time = tm.elapsed();
    std::cout << "lj0 " << time << std::endl;
    std::cout << energy << std::endl;
  }

  if((testtype==0)||(testtype==9)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      ljc.loopijsets(ljisets,ljjsets,particlei,particlej,energy);
    }
    time = tm.elapsed();
    std::cout << "ljc " << time << std::endl;
    std::cout << energy << std::endl;
  }

  if((testtype==0)||(testtype==10)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      ljewald.loopijsets(ljisets,ljjsets,particlei,particlej,energy);
    }
    time = tm.elapsed();
    std::cout << "ljewald " << time << std::endl;
    std::cout << energy << std::endl;
  }
  
  coulombisets[0].end = 0;

  if((testtype==0)||(testtype==11)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      ewald_and_lje.loopijsets(ljisets,ljjsets,coulombisets,coulombjsets,particlei,particlej,energy);
    }
    time = tm.elapsed();
    std::cout << "ewald_and_lje (lje only) " << time << std::endl;
    std::cout << energy << std::endl;
  }

  if((testtype==0)||(testtype==12)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      ewald_and_lje.loopijsets(ljisets,ljjsets,coulombisets,coulombjsets,particlei,particlej,energy,cutoff2);
    }
    time = tm.elapsed();
    std::cout << "ewald_and_lje (lje only) cutoff " << time << std::endl;
    std::cout << energy << std::endl;
  }

  if((testtype==0)||(testtype==13)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      ewald_and_lje.ljc.loopijsets(ljisets,ljjsets,particlei,particlej,energy,cutoff2);
    }
    time = tm.elapsed();
    std::cout << "ewald_and_lje.ljc (lje only) cutoff " << time << std::endl;
    std::cout << energy << std::endl;
  }

  if((testtype==0)||(testtype==14)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      ewald_and_lje.lj.loopijsets(ljisets,ljjsets,particlei,particlej,energy,cutoff2);
    }
    time = tm.elapsed();
    std::cout << "ewald_and_lje.lj (lje only) cutoff " << time << std::endl;
    std::cout << energy << std::endl;
  }

  if((testtype==0)||(testtype==15)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      cl_and_ljc.lj.loopijsets(ljisets,ljjsets,particlei,particlej,energy,cutoff2);
    }
    time = tm.elapsed();
    std::cout << "cl_and_ljc.lj (ljc only) cutoff " << time << std::endl;
    std::cout << energy << std::endl;
  }

  /// all i particlearray are coulomb
  coulombisets[0].begin=0;
  coulombisets[0].end=ni;
  /// all j particlearray are coulomb
  coulombjsets[0][0].begin=0;
  coulombjsets[0][0].end=nj;
  ljjsets[0][0].begin=nj;
  ljjsets[0][0].end=nj;

  if((testtype==0)||(testtype==16)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      cl_and_ljc.coulomb.loopijsets(coulombisets,coulombjsets,particlei,particlej,energy,cutoff2);
    }
    time = tm.elapsed();
    std::cout << "cl_and_ljc.coulomb (cl only) cutoff " << time << std::endl;
    std::cout << energy << std::endl;
  }

  if((testtype==0)||(testtype==17)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      ewald_and_lje.coulomb.loopijsets(coulombisets,coulombjsets,particlei,particlej,energy,cutoff2);
    }
    time = tm.elapsed();
    std::cout << "ewald_and_lje.coulomb (cl only) cutoff " << time << std::endl;
    std::cout << energy << std::endl;
  }


  if((testtype==0)||(testtype==18)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      ewald_and_lje.lj.loopiself(particlei,energy,cutoff2);
    }
    time = tm.elapsed();
    std::cout << "ewald_and_lje.lj (lje only i self) cutoff " << time << std::endl;
    std::cout << energy << std::endl;
  }

  if((testtype==0)||(testtype==19)){
    zero_force(particlei,ni);
    zero_force(particlej,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      ewald_and_lje.coulomb.loopiself(particlei,energy,cutoff2);
    }
    time = tm.elapsed();
    std::cout << "ewald_and_lje.coulomb (cl only i self) cutoff " << time << std::endl;
    std::cout << energy << std::endl;
  }

  int nsub =10;
  int nset = ni/nsub;
  std::vector<TypeRange> tr(nset);
  for(int s=0;s<nset;s++){
    tr[s].begin = s*nsub;
    tr[s].end = (s+1)*nsub;
    if(tr[s].end>=ni)tr[s].end=ni;
    tr[s].lj.begin = tr[s].begin;
    tr[s].lj.end = tr[s].end;
    tr[s].ljcoulomb.begin = tr[s].end;
    tr[s].ljcoulomb.end = tr[s].end;
    tr[s].coulomb.begin = tr[s].end;
    tr[s].coulomb.end = tr[s].end;
  }
  std::vector< std::vector<int> > targetindex(nset);
  for(int s=0;s<nset;s++){
    for(int ts=0;ts<nset;ts++){
      if(ts!=s)targetindex[s].push_back(ts);
    }
  } 

  ShortRangeInteractionSet <ShortRange::LJ,LJType,Particle,Particle> ljset;
  ShortRangeInteractionSet <ShortRange::LJCoulomb,LJCoulombType,Particle,Particle> ljcset;
 
  if((testtype==0)||(testtype==20)){
    zero_force(particlearrayi,ni);
    zero_force(particlearrayj,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      ljset.loopself(particlearrayi,tr,targetindex,energy);
    }
    time = tm.elapsed();
    std::cout << "ljset (lj only i self)" << time << std::endl;
    std::cout << energy << std::endl;
  }

  int njset = nj/nsub;
  std::vector<TypeRange> gtr(nset);
  for(int s=0;s<njset;s++){
    gtr[s].begin = s*nsub;
    gtr[s].end = (s+1)*nsub;
    if(tr[s].end>=ni)gtr[s].end=ni;
    gtr[s].lj.begin = gtr[s].begin;
    gtr[s].lj.end = gtr[s].end;
    gtr[s].ljcoulomb.begin = gtr[s].end;
    gtr[s].ljcoulomb.end = gtr[s].end;
    gtr[s].coulomb.begin = gtr[s].end;
    gtr[s].coulomb.end = gtr[s].end;
  }
  std::vector< std::vector<int> > gtargetindex(nset);
  for(int s=0;s<nset;s++){
    for(int ts=0;ts<njset;ts++){
      gtargetindex[s].push_back(ts);
    }
  } 
  currenttest=21;
  if((testtype==0)||(testtype==currenttest)){

    std::cout << "ljset (lj only ij) short";
    for(int t=0;t<5;t++){
      zero_force(particlearrayi,ni);
      zero_force(particlearrayj,nj);
      energy = 0.0;
      tm.restart();  
      for(int r=0;r<rl;r++){
        ljset.loopisetjsetij(particlearrayi,tr,particlearrayj,gtr,gtargetindex,energy);
      }
      time = tm.elapsed();
      std::cout << " " << time;
    }
    std::cout << std::endl;
    std::cout << energy << std::endl;
  }

  currenttest++;
  if((testtype==0)||(testtype==currenttest)){

    std::cout << "ljset (loopijset) short";
    for(int t=0;t<5;t++){
      zero_force(particlearrayi,ni);
      zero_force(particlearrayj,nj);
      energy = 0.0;
      tm.restart();  
      for(int r=0;r<rl;r++){
        ljset.loopijset(particlearrayi,tr,particlearrayj,gtr,gtargetindex,energy);
      }
      time = tm.elapsed();
      std::cout << " " << time;
    }
    std::cout << std::endl;
    std::cout << energy << std::endl;
  }

  CoulombAndLJInteractionSet cljset;
  currenttest++;
  if((testtype==0)||(testtype==currenttest)){
    std::cout << "cljset () short";
    for(int t=0;t<5;t++){
      zero_force(particlearrayi,ni);
      zero_force(particlearrayj,nj);
      energy = 0.0;
      tm.restart();  
      for(int r=0;r<rl;r++){
        cljset.allloop(particlearrayi,tr,targetindex,particlearrayj,gtr,gtargetindex,energy);
      }
      time = tm.elapsed();
      std::cout << " " << time;
    }
    std::cout << std::endl;
    std::cout << energy << std::endl;
#ifdef RDTSCTIME
    std::cout << "lj " << cljset.ljtime << " ljc " << cljset.ljctime << " cl " << cljset.cltime << std::endl;
#endif
  }


  tr.resize(1);
  tr[0].begin=0;
  tr[0].end=ni;
  tr[0].lj.begin=tr[0].begin;
  tr[0].lj.end=tr[0].end;
  gtr.resize(1);
  gtr[0].begin=0;
  gtr[0].end=nj;
  gtr[0].lj.begin=gtr[0].begin;
  gtr[0].lj.end=gtr[0].end;
  gtargetindex.resize(1);
  gtargetindex[0].resize(1);
  gtargetindex[0][0] = 0;
  currenttest++;
  if((testtype==0)||(testtype==currenttest)){
    std::cout << "ljset (lj only ij) long";
    for(int t=0;t<5;t++){    
      zero_force(particlearrayi,ni);
      zero_force(particlearrayj,nj);
      energy = 0.0;
      tm.restart();  
      for(int r=0;r<rl;r++){
        ljset.loopisetjsetij(particlearrayi,tr,particlearrayj,gtr,gtargetindex,energy);
      }
      time = tm.elapsed();
      std::cout << " " << time;
    }
    std::cout << std::endl;
    std::cout << energy << std::endl;
  }
  currenttest++;
  if((testtype==0)||(testtype==currenttest)){
    std::cout << "ljset (loopijset) long";
    for(int t=0;t<5;t++){    
      zero_force(particlearrayi,ni);
      zero_force(particlearrayj,nj);
      energy = 0.0;
      tm.restart();  
      for(int r=0;r<rl;r++){
        ljset.loopijset(particlearrayi,tr,particlearrayj,gtr,gtargetindex,energy);
      }
      time = tm.elapsed();
      std::cout << " " << time;
    }
    std::cout << std::endl;
    std::cout << energy << std::endl;
  }

  tr.resize(nset);
  for(int s=0;s<nset;s++){
    tr[s].begin = s*nsub;
    tr[s].end = (s+1)*nsub;
    if(tr[s].end>=ni)tr[s].end=ni;
    tr[s].lj.begin = tr[s].begin;
    tr[s].lj.end = tr[s].begin;
    tr[s].ljcoulomb.begin = tr[s].begin;
    tr[s].ljcoulomb.end = tr[s].end;
    tr[s].coulomb.begin = tr[s].end;
    tr[s].coulomb.end = tr[s].end;
  }
  gtr.resize(njset);
  for(int s=0;s<njset;s++){
    gtr[s].begin = s*nsub;
    gtr[s].end = (s+1)*nsub;
    if(tr[s].end>=ni)gtr[s].end=ni;
    gtr[s].lj.begin = gtr[s].begin;
    gtr[s].lj.end = gtr[s].begin;
    gtr[s].ljcoulomb.begin = gtr[s].begin;
    gtr[s].ljcoulomb.end = gtr[s].end;
    gtr[s].coulomb.begin = gtr[s].end;
    gtr[s].coulomb.end = gtr[s].end;
  }
  gtargetindex.resize(nset);
  for(int s=0;s<nset;s++){
    gtargetindex[s].resize(njset);
    gtargetindex[s].clear();
    for(int ts=0;ts<njset;ts++){
      gtargetindex[s].push_back(ts);
    }
  }   
  currenttest++;
  if((testtype==0)||(testtype==currenttest)){
    zero_force(particlearrayi,ni);
    zero_force(particlearrayj,nj);
    energy = 0.0;
    tm.restart();  
    for(int r=0;r<rl;r++){
      ljcset.loopisetjsetij(particlearrayi,tr,particlearrayj,gtr,gtargetindex,energy);
    }
    time = tm.elapsed();
    std::cout << "ljcset (ljc only ij)" << time << std::endl;
    std::cout << energy << std::endl;
  }


  /*
  for(long i=0;i<ni;i++){
    std::cout <<  particlei[i].force.x << std::endl;
  }
  */
}
