#pragma once
/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>

PS::F64 W(const PS::F64 r, const PS::F64 h);
PS::F64vec gradW(const PS::F64vec dr, const PS::F64 h);
PS::F64 dWdh(const PS::F64 r, const PS::F64 h);
