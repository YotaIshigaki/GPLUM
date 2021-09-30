#ifndef LOG_H
#define LOG_H

#include <iostream>
#include "Common.h"
#include "ParticleInfo.h"
#include "CovalentBondInfo.h"

void dump_typerange(std::vector<TypeRange>& typerangearray);

template<typename PA>
void dump_particle(PA& particlearray, 
          std::vector<TypeRange>& typerangearray, 
          std::vector<int>& setid,
          double energy);

template<typename PA>
void dump_atomid(const PA& particlearray, 
                 const std::vector<TypeRange>& typerangearray, 
                 const std::vector<int>& setid);

void dump_bondlistarray(std::vector<CovalentBondInfo::BondList> bondlistarray);

template<class PA>
void dump_shakelist(ShakeList& shakelist, PA& particlearray);

template<class PA>
void dump_shakelist_distance(ShakeList& shakelist, PA& particlearray);

#endif

