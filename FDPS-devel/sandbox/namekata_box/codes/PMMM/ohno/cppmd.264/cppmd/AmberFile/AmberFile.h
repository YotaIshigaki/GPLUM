#ifndef AMBERFILE_AMBERFILE_H
#define AMBERFILE_AMBERFILE_H

#include <string>

#include "ParticleInfo.h"
#include "CovalentBondInfo.h"

template<class PA>
int ReadAmberFile(const std::string& prmtop_file,
                  const std::string& restrt_file,
                  PA* particle_array,
                  CovalentBondList* covalent_bond_list,
                  CovalentBondParameterList* covalent_bond_parameter_list,
                  SpaceVector<double>* box,
                  SpaceVector<double>* skew,
                  int restart,
                  double *start_time);

#endif  // AMBERFILE_AMBERFILE_H
