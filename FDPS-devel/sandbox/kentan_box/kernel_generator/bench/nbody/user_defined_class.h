#ifndef H_USER_DEFINED_CLASS
#define H_USER_DEFINED_CLASS

#include <particle_simulator.hpp>

struct EPI{
  PS::F64vec pos;
  PS::F32    eps2;
};

struct EPJ{
  PS::F64vec pos;
  PS::F32    mass;
  PS::F32    eps2;
};

struct Force{
  PS::F32vec acc;
  PS::F32    pot;
};

#endif
