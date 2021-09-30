#pragma once
/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>

// SPH kernel used for the usual SPH calculation
PS::F64 W(const PS::F64 r, const PS::F64 h);
PS::F64vec gradW(const PS::F64vec dr, const PS::F64 h);
PS::F64 dWdh(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdrdh(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdh2(const PS::F64 r, const PS::F64 h);

// SPH kernel used for the constraint for Lagrangian
PS::F64 W_cstr(const PS::F64 r, const PS::F64 h);
PS::F64vec gradW_cstr(const PS::F64vec dr, const PS::F64 h);
PS::F64 dWdh_cstr(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdrdh_cstr(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdh2_cstr(const PS::F64 r, const PS::F64 h);

// W or W_cstr are choosen from the kernels listed below.
PS::F64 W_BSpline_M4(const PS::F64 r, const PS::F64 h);
PS::F64vec gradW_BSpline_M4(const PS::F64vec dr, const PS::F64 h);
PS::F64 dWdh_BSpline_M4(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdrdh_BSpline_M4(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdh2_BSpline_M4(const PS::F64 r, const PS::F64 h);

PS::F64 W_BSpline_M5(const PS::F64 r, const PS::F64 h);
PS::F64vec gradW_BSpline_M5(const PS::F64vec dr, const PS::F64 h);
PS::F64 dWdh_BSpline_M5(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdrdh_BSpline_M5(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdh2_BSpline_M5(const PS::F64 r, const PS::F64 h);

PS::F64 W_BSpline_M6(const PS::F64 r, const PS::F64 h);
PS::F64vec gradW_BSpline_M6(const PS::F64vec dr, const PS::F64 h);
PS::F64 dWdh_BSpline_M6(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdrdh_BSpline_M6(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdh2_BSpline_M6(const PS::F64 r, const PS::F64 h);

PS::F64 W_Wendland_C2(const PS::F64 r, const PS::F64 h);
PS::F64vec gradW_Wendland_C2(const PS::F64vec dr, const PS::F64 h);
PS::F64 dWdh_Wendland_C2(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdrdh_Wendland_C2(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdh2_Wendland_C2(const PS::F64 r, const PS::F64 h);

PS::F64 W_Wendland_C4(const PS::F64 r, const PS::F64 h);
PS::F64vec gradW_Wendland_C4(const PS::F64vec dr, const PS::F64 h);
PS::F64 dWdh_Wendland_C4(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdrdh_Wendland_C4(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdh2_Wendland_C4(const PS::F64 r, const PS::F64 h);

PS::F64 W_Wendland_C6(const PS::F64 r, const PS::F64 h);
PS::F64vec gradW_Wendland_C6(const PS::F64vec dr, const PS::F64 h);
PS::F64 dWdh_Wendland_C6(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdrdh_Wendland_C6(const PS::F64 r, const PS::F64 h);
PS::F64 d2Wdh2_Wendland_C6(const PS::F64 r, const PS::F64 h);

// check routine
void outputSPHKernelProfile(const PS::S32 root = 0);
