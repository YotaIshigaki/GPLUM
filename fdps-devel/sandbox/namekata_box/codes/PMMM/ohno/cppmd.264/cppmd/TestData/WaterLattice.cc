#include <string>
#include <cstdio>
#include <iostream>
#include "SpaceVector.h"
#include "ParticleInfo.h"
#include "Common.h"
#include <stdlib.h>

#include "WaterLattice.h"

#define WATERLATTICESPACING 10
//#define OW_AMBER94 15
//#define HW_AMBER94 12

template<class PA>
WaterLattice<PA>:: WaterLattice()
{
  this->latticeSpacing( WATERLATTICESPACING, WATERLATTICESPACING, WATERLATTICESPACING );
  this->latticeNum(1);
  setLattice();
}

template<class PA>
WaterLattice<PA>::WaterLattice(int lnum)
{
  this->latticeSpacing( WATERLATTICESPACING, WATERLATTICESPACING, WATERLATTICESPACING );
  this->latticeNum(lnum);
  this->side.x = this->latticeSpacing.x * this->latticeNum.x;
  this->side.y = this->latticeSpacing.y * this->latticeNum.y;
  this->side.z = this->latticeSpacing.z * this->latticeNum.z;
  setLattice();
}

template<class PA>
WaterLattice<PA>::WaterLattice(SpaceVector<int> lnum)
{
  this->latticeSpacing( WATERLATTICESPACING, WATERLATTICESPACING, WATERLATTICESPACING );
  this->latticeNum= lnum;
  this->side.x = this->latticeSpacing.x * this->latticeNum.x;
  this->side.y = this->latticeSpacing.y * this->latticeNum.y;
  this->side.z = this->latticeSpacing.z * this->latticeNum.z;
  setLattice();
}

template<class PA>
int
WaterLattice<PA>::setLattice()
{
  int n = 0;

  Particle ow = Particle();
  Particle hw = Particle();
  SpaceVector<double> pos, vel;

  ow.velocity(0.0,0.0,0.0);
  ow.force(0.0,0.0,0.0);
  ow.mass=16.00/12.01;
  ow.inv_mass=1.0/ow.mass;
  ow.charge=-0.834;
  ow.atomtype=ATOMTYPE_WO;

  hw.velocity(0.0,0.0,0.0);
  hw.force(0.0,0.0,0.0);
  hw.mass=1.008/12.01;
  hw.inv_mass=1.0/hw.mass;
  hw.charge=0.417;
  hw.atomtype=ATOMTYPE_WH;

  this->side.x = this->latticeSpacing.x * this->latticeNum.x;
  this->side.y = this->latticeSpacing.y * this->latticeNum.y;
  this->side.z = this->latticeSpacing.z * this->latticeNum.z;
  this->particle.clear();
  this->potentialmodel.clear();

//  srand(2);
  srand(48);
  /*
    pos.x = 5.0;
    pos.y = 5.0;
    pos.z = 5.0;
    vel.x = 1.0;
    vel.y = 1.1;
    vel.z = 0.9;

    ow.atomid = 0;
    ow.position = pos;
    ow.velocity  = vel;
    this->particle.push_back(ow);
    potentialtype.push_back(LJCoulombParticle);

    SpaceVector<double> tran(0);

    hw.atomid = 1;
    tran.x = sin(104.520/2.0/180*M_PI) * 0.9572;
    tran.y = cos(104.520/2.0/180*M_PI) * 0.9572;
    tran.z = 0.0;
    hw.position = pos + tran;
    hw.velocity = 0.002 * vel;
    this->particle.push_back(hw);
    potentialtype.push_back(CoulombParticle);
	
    hw.atomid = 2;
    tran.x = -sin(104.520/2.0/180*M_PI) * 0.9572;
    tran.y = cos(104.520/2.0/180*M_PI) * 0.9572;
    tran.z = 0.0;
    hw.position = pos + tran;
    this->particle.push_back(hw);
    potentialtype.push_back(CoulombParticle);
  */

  ow.atomid=n;
  for(int i = 0; i < this->latticeNum.x; i++) {
    //pos.x = this->latticeSpacing.x * i + 30;
    pos.x = this->latticeSpacing.x * i + 5;
    vel.x = 1 - 2 * (rand()/(RAND_MAX+1.0));
    //vel.x =  (rand()/(RAND_MAX+1.0));
    for(int j = 0; j < this->latticeNum.y; j++) {
      //pos.y = this->latticeSpacing.y * j + 30;
      pos.y = this->latticeSpacing.y * j + 5;
      vel.y = 1 - 2 * (rand()/(RAND_MAX+1.0));
      //vel.y =  (rand()/(RAND_MAX+1.0));
      for(int k = 0; k < this->latticeNum.z; k++) {
        //for(int k = 0; k < 1; k++) {
        vel.z = 1 - 2 * (rand()/(RAND_MAX+1.0));
        //vel.z = (rand()/(RAND_MAX+1.0));
        SpaceVector<double> tran(0);
        //pos.z = this->latticeSpacing.z * k + 30;
        pos.z = this->latticeSpacing.z * k + 5;
        ow.atomid = n;
        ow.position = pos + tran;
        ow.velocity = vel * 0.0;
        //ow.velocity = vel * 0.1;
        //ow.velocity = vel * 0.02;
        //ow.velocity = vel * 0.2;
        //ow.velocity.x = 0.02;
        this->particle.push_back(ow);
        this->potentialmodel.push_back(LJCoulombPotential);
        ++n;
        hw.atomid = n;
        tran.x = sin(104.520/2.0/180*M_PI) * 0.9572;
        tran.y = cos(104.520/2.0/180*M_PI) * 0.9572;
        tran.z = 0.0;
        hw.position = pos + tran;
        //hw.velocity = vel * 0.02;
        this->particle.push_back(hw);
        this->potentialmodel.push_back(OnlyCoulombPotential);
        //potentialtype.push_back(LJCoulombParticle);
        ++n;
        hw.atomid = n;
        tran.x = -sin(104.520/2.0/180*M_PI) * 0.9572;
        tran.y = cos(104.520/2.0/180*M_PI) * 0.9572;
        tran.z = 0.0;
        hw.position = pos + tran;
//        hw.velocity = vel * 0.02;
        this->particle.push_back(hw);
        this->potentialmodel.push_back(OnlyCoulombPotential);
        //potentialtype.push_back(LJCoulombParticle);
        ++n;
      }
    }
  }

  this->natom = n;
  return n;
}

template<class PA>
void
WaterLattice<PA>::writePDB()
{
  int n = 0, res = 1;
  char aname[2][5] = {"OW  ","HW  "};
  char rname[2][4] = {"WAT","WAT"};

  size_t size;
  size = this->particle.size();
  for(n=0;size<n;n++){
    int atomtype = 0;
    if (getatomtype(this->particle,n) == ATOMTYPE_WO) atomtype = 0;
    if (getatomtype(this->particle,n) == ATOMTYPE_WH) atomtype = 1;
    //    printf("ATOM  %5d %4s %3s %5d    %8.3f%8.3f%8.3f\nTER\n",
    printf("ATOM  %5d %4s %3s %5d    %8.3f%8.3f%8.3f\n",
           n+1, aname[atomtype], rname[atomtype], res,
           getpos(this->particle,n).x, getpos(this->particle,n).y, getpos(this->particle,n).z);
    if(n%3==0) {
      printf("TER\n");
      res++;
    }
  }
  printf("END\n");
}

#ifdef OLDPARTICLE
template class
WaterLattice<ParticleArray>;
#else
template class
WaterLattice<CombinedParticleArray>;
#endif
