#include <string>
#include <cstdio>
#include <iostream>
#include "SpaceVector.h"
#include "ParticleInfo.h"
#include "Common.h"

#include "NaClFCC.h"

#define NACLLATTICESPACEING 5.641
#define IP_AMBER94 28
#define IM_AMBER94 26

template<class PA>
NaClFCC<PA>::NaClFCC()
{
  this->latticeSpacing( NACLLATTICESPACEING, NACLLATTICESPACEING, NACLLATTICESPACEING );
  this->latticeNum(1);
  setLattice();
}

template<class PA>
NaClFCC<PA>::NaClFCC(int lnum)
{
  this->latticeSpacing( NACLLATTICESPACEING, NACLLATTICESPACEING, NACLLATTICESPACEING );
  this->latticeNum(lnum);
  this->side.x = this->latticeSpacing.x * this->latticeNum.x;
  this->side.y = this->latticeSpacing.y * this->latticeNum.y;
  this->side.z = this->latticeSpacing.z * this->latticeNum.z;
  setLattice();
}

template<class PA>
NaClFCC<PA>::NaClFCC(SpaceVector<int> lnum)
{
  this->latticeSpacing( NACLLATTICESPACEING, NACLLATTICESPACEING, NACLLATTICESPACEING );
  this->latticeNum = lnum;
  this->side.x = this->latticeSpacing.x * this->latticeNum.x;
  this->side.y = this->latticeSpacing.y * this->latticeNum.y;
  this->side.z = this->latticeSpacing.z * this->latticeNum.z;
  setLattice();
}

template<class PA>
int
NaClFCC<PA>::setLattice()
{
  int n = 0;
  Particle ip = Particle(); // Ion Plus (IP)
  Particle im = Particle(); // Ion Minus (IM)
  SpaceVector<double> pos;

  ip.velocity(0.0,0.0,0.0);
  ip.force(0.0,0.0,0.0);
  ip.mass=22.99/12.01;
  ip.inv_mass=1.0/ip.mass;
  ip.charge = 1.0;
  ip.atomtype = IP_AMBER94;

  im.velocity(0.0,0.0,0.0);
  im.force(0.0,0.0,0.0);
  im.mass=35.45/12.01;
  im.inv_mass=1.0/im.mass;
  im.charge = -1.0;
  im.atomtype = IM_AMBER94;

  this->side.x = this->latticeSpacing.x * this->latticeNum.x;
  this->side.y = this->latticeSpacing.y * this->latticeNum.y;
  this->side.z = this->latticeSpacing.z * this->latticeNum.z;
  this->particle.clear();
  for (int i = 0; i < this->latticeNum.x * 2; i++) {
    pos.x = i * (this->latticeSpacing.x/2.0);
    for (int j = 0; j < this->latticeNum.y * 2; j++) {
      pos.y = j * (this->latticeSpacing.y/2.0);
      for (int k = 0; k < this->latticeNum.z * 2; k++) {
          pos.z = k * (this->latticeSpacing.z/2.0);
          if ((i+j+k)%2 == 0) {
            ip.atomid = n;
            ip.position = pos;
            this->particle.push_back(ip);
            ++n;
          } else {
            im.atomid = n;
            im.position = pos;
            this->particle.push_back(im);
            ++n;
        }
      }
    }
  }
  this->potentialmodel.resize(n,LJCoulombPotential);
  this->natom = n;
  return n;
}


template<class PA>
void
NaClFCC<PA>::writePDB()
{
  int n = 0;
  char aname[2][5] = {"Na+ ","Cl- "};
  char rname[2][4] = {"Na+","Cl-"};

  size_t size;
  size = this->particle.size();
  for(n=0;size<n;n++){
    int atomtype = 0;
    if (getatomtype(this->particle,n) == IP_AMBER94) atomtype = 0;
    if (getatomtype(this->particle,n) == IM_AMBER94) atomtype = 1;
    printf("ATOM  %5d %4s %3s %5d    %8.3f%8.3f%8.3f\nTER\n",
           n+1, aname[atomtype], rname[atomtype], n+1,
           getpos(this->particle,n).x, getpos(this->particle,n).y, getpos(this->particle,n).z);
  }
  printf("END\n");
}

#ifdef OLDPARTICLE
template class
NaClFCC<ParticleArray>;
#else
template class
NaClFCC<CombinedParticleArray>;
#endif
