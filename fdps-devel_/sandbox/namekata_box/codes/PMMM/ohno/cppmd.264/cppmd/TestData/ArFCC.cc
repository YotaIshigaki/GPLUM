#include <string>
#include <cstdio>
#include <iostream>
#include "SpaceVector.h"
#include "ParticleInfo.h"
#include "Common.h"

#include "ArFCC.h"

#define ARLATTICESPACEING 5.77175

template<class PA>
ArFCC<PA>::ArFCC()
//  : NoInputSystem<PA>::latticeSpacing( ARLATTICESPACEING, ARLATTICESPACEING, ARLATTICESPACEING )
{
  this->latticeSpacing( ARLATTICESPACEING, ARLATTICESPACEING, ARLATTICESPACEING );
  this->latticeNum(1);
  setLattice();
}

template<class PA>
ArFCC<PA>::ArFCC(int lnum)
//  : latticeSpacing( ARLATTICESPACEING, ARLATTICESPACEING, ARLATTICESPACEING ),
//    latticeNum(lnum)
{
  this->latticeSpacing( ARLATTICESPACEING, ARLATTICESPACEING, ARLATTICESPACEING );
  this->latticeNum(lnum);
  this->side.x = this->latticeSpacing.x * this->latticeNum.x;
  this->side.y = this->latticeSpacing.y * this->latticeNum.y;
  this->side.z = this->latticeSpacing.z * this->latticeNum.z;
  setLattice();
}

template<class PA>
ArFCC<PA>::ArFCC(SpaceVector<int> lnum)
//  : latticeSpacing( ARLATTICESPACEING, ARLATTICESPACEING, ARLATTICESPACEING ),
//    latticeNum(lnum)
{
  this->latticeSpacing( ARLATTICESPACEING, ARLATTICESPACEING, ARLATTICESPACEING );
  this->latticeNum = lnum;
  this->side.x = this->latticeSpacing.x * this->latticeNum.x;
  this->side.y = this->latticeSpacing.y * this->latticeNum.y;
  this->side.z = this->latticeSpacing.z * this->latticeNum.z;
  setLattice();
}

template<class PA>
int
ArFCC<PA>::setLattice()
{
  int n = 0;
  Particle p = Particle();

  p.velocity(0.0,0.0,0.0);
  p.force(0.0,0.0,0.0);
  p.mass=3.32603;
  p.inv_mass=1.0/p.mass;
  p.atomtype = 0;

  this->side.x = this->latticeSpacing.x * this->latticeNum.x;
  this->side.y = this->latticeSpacing.y * this->latticeNum.y;
  this->side.z = this->latticeSpacing.z * this->latticeNum.z;
  this->particle.clear();
  for (int i = 0; i < this->latticeNum.x * 2; i++) {
    p.position.x = i * (this->latticeSpacing.x/2.0);
    for (int j = 0; j < this->latticeNum.y * 2; j++) {
      p.position.y = j * (this->latticeSpacing.y/2.0);
      for (int k = 0; k < this->latticeNum.z * 2; k++) {
        if ((i+j+k)%2 == 0) {
          p.position.z = k * (this->latticeSpacing.z/2.0);
          p.atomid = n;
          this->particle.push_back(p);
          ++n;
        }
      }
    }
  }
  this->potentialmodel.resize(n,OnlyLJPotential);
  this->natom = n;
  return n;
}

template<class PA>
void
ArFCC<PA>::writePDB()
{
  int n = 0;
  char aname[5] = "AR  ";
  char rname[4] = "AR ";
  
  size_t size;
  size = this->particle.size();
  for(n=0;size<n;n++){
    printf("ATOM  %5d %4s %3s %5d    %8.3f%8.3f%8.3f\nTER\n",
           n+1, aname, rname, n+1,
           getpos(this->particle,n).x, getpos(this->particle,n).y, getpos(this->particle,n).z);
  }
  printf("END\n");
}

#ifdef OLDPARTICLE
template class
ArFCC<ParticleArray>;
#else
template class
ArFCC<CombinedParticleArray>;
#endif
