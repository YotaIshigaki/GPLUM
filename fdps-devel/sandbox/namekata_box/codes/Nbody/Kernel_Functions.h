#pragma once

#include <particle_simulator.hpp>

namespace kernel_functions {

extern PS::F64 W_SPH(const PS::F64 r, const PS::F64 h);
extern PS::F64 dWdr_SPH(const PS::F64 r, const PS::F64 h);
extern PS::F64 dWdh_SPH(const PS::F64 r, const PS::F64 h);
extern PS::F64 phi_SPH(const PS::F64 r, const PS::F64 h);
extern PS::F64 dphidr_SPH(const PS::F64 r, const PS::F64 h);
extern PS::F64 dphidh_SPH(const PS::F64 r, const PS::F64 h);
extern void check_SPH_kernels();

} // kernel_functions
