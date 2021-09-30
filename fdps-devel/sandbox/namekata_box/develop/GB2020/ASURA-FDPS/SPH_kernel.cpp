/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "mathematical_constants.h"
#include "SPH_kernel.h"


PS::F64 W(const PS::F64 r, const PS::F64 h) {
#if ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_BSPLINE_M4
    return W_BSpline_M4(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_BSPLINE_M5
    return W_BSpline_M5(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_BSPLINE_M6
    return W_BSpline_M6(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_WENDLAND_C2
    return W_Wendland_C2(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_WENDLAND_C4
    return W_Wendland_C4(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_WENDLAND_C6
    return W_Wendland_C6(r,h);
#else
#error The macro ASURA_FDPS_SPH_KERNEL is not defined or its value is incorrect.
#endif
}

PS::F64vec gradW(const PS::F64vec dr, const PS::F64 h) {
#if ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_BSPLINE_M4
    return gradW_BSpline_M4(dr,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_BSPLINE_M5
    return gradW_BSpline_M5(dr,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_BSPLINE_M6
    return gradW_BSpline_M6(dr,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_WENDLAND_C2
    return gradW_Wendland_C2(dr,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_WENDLAND_C4
    return gradW_Wendland_C4(dr,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_WENDLAND_C6
    return gradW_Wendland_C6(dr,h);
#else
#error The macro ASURA_FDPS_SPH_KERNEL is not defined or its value is incorrect.
#endif
}

PS::F64 dWdh(const PS::F64 r, const PS::F64 h) {
#if ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_BSPLINE_M4
    return dWdh_BSpline_M4(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_BSPLINE_M5
    return dWdh_BSpline_M5(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_BSPLINE_M6
    return dWdh_BSpline_M6(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_WENDLAND_C2
    return dWdh_Wendland_C2(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_WENDLAND_C4
    return dWdh_Wendland_C4(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_WENDLAND_C6
    return dWdh_Wendland_C6(r,h);
#else
#error The macro ASURA_FDPS_SPH_KERNEL is not defined or its value is incorrect.
#endif
}

PS::F64 d2Wdrdh(const PS::F64 r, const PS::F64 h) {
#if ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_BSPLINE_M4
    return d2Wdrdh_BSpline_M4(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_BSPLINE_M5
    return d2Wdrdh_BSpline_M5(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_BSPLINE_M6
    return d2Wdrdh_BSpline_M6(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_WENDLAND_C2
    return d2Wdrdh_Wendland_C2(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_WENDLAND_C4
    return d2Wdrdh_Wendland_C4(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_WENDLAND_C6
    return d2Wdrdh_Wendland_C6(r,h);
#else
#error The macro ASURA_FDPS_SPH_KERNEL is not defined or its value is incorrect.
#endif
}

PS::F64 d2Wdh2(const PS::F64 r, const PS::F64 h) {
#if ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_BSPLINE_M4
    return d2Wdh2_BSpline_M4(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_BSPLINE_M5
    return d2Wdh2_BSpline_M5(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_BSPLINE_M6
    return d2Wdh2_BSpline_M6(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_WENDLAND_C2
    return d2Wdh2_Wendland_C2(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_WENDLAND_C4
    return d2Wdh2_Wendland_C4(r,h);
#elif ASURA_FDPS_SPH_KERNEL == ASURA_FDPS_WENDLAND_C6
    return d2Wdh2_Wendland_C6(r,h);
#else
#error The macro ASURA_FDPS_SPH_KERNEL is not defined or its value is incorrect.
#endif
}

PS::F64 W_cstr(const PS::F64 r, const PS::F64 h) {
#ifndef ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT
    return W(r,h);
#else
#if ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_BSPLINE_M4
    return W_BSpline_M4(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_BSPLINE_M5
    return W_BSpline_M5(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_BSPLINE_M6
    return W_BSpline_M6(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_WENDLAND_C2
    return W_Wendland_C2(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_WENDLAND_C4
    return W_Wendland_C4(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_WENDLAND_C6
    return W_Wendland_C6(r,h);
#else
    return W(r,h);
#endif
#endif
}

PS::F64vec gradW_cstr(const PS::F64vec dr, const PS::F64 h) {
#ifndef ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT
    return gradW(dr,h);
#else
#if ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_BSPLINE_M4
    return gradW_BSpline_M4(dr,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_BSPLINE_M5
    return gradW_BSpline_M5(dr,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_BSPLINE_M6
    return gradW_BSpline_M6(dr,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_WENDLAND_C2
    return gradW_Wendland_C2(dr,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_WENDLAND_C4
    return gradW_Wendland_C4(dr,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_WENDLAND_C6
    return gradW_Wendland_C6(dr,h);
#else
    return gradW(dr,h);
#endif
#endif
}

PS::F64 dWdh_cstr(const PS::F64 r, const PS::F64 h) {
#ifndef ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT
    return dWdh(r,h);
#else
#if ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_BSPLINE_M4
    return dWdh_BSpline_M4(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_BSPLINE_M5
    return dWdh_BSpline_M5(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_BSPLINE_M6
    return dWdh_BSpline_M6(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_WENDLAND_C2
    return dWdh_Wendland_C2(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_WENDLAND_C4
    return dWdh_Wendland_C4(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_WENDLAND_C6
    return dWdh_Wendland_C6(r,h);
#else
    return dWdh(r,h);
#endif
#endif
}

PS::F64 d2Wdrdh_cstr(const PS::F64 r, const PS::F64 h) {
#ifndef ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT
    return d2Wdrdh(r,h);
#else
#if ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_BSPLINE_M4
    return d2Wdrdh_BSpline_M4(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_BSPLINE_M5
    return d2Wdrdh_BSpline_M5(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_BSPLINE_M6
    return d2Wdrdh_BSpline_M6(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_WENDLAND_C2
    return d2Wdrdh_Wendland_C2(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_WENDLAND_C4
    return d2Wdrdh_Wendland_C4(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_WENDLAND_C6
    return d2Wdrdh_Wendland_C6(r,h);
#else
    return d2Wdrdh(r,h);
#endif
#endif
}

PS::F64 d2Wdh2_cstr(const PS::F64 r, const PS::F64 h) {
#ifndef ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT
    return d2Wdh2(r,h);
#else
#if ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_BSPLINE_M4
    return d2Wdh2_BSpline_M4(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_BSPLINE_M5
    return d2Wdh2_BSpline_M5(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_BSPLINE_M6
    return d2Wdh2_BSpline_M6(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_WENDLAND_C2
    return d2Wdh2_Wendland_C2(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_WENDLAND_C4
    return d2Wdh2_Wendland_C4(r,h);
#elif ASURA_FDPS_SPH_KERNEL_FOR_CONSTRAINT == ASURA_FDPS_WENDLAND_C6
    return d2Wdh2_Wendland_C6(r,h);
#else
    return d2Wdh2(r,h);
#endif
#endif
}

PS::F64 W_BSpline_M4(const PS::F64 r, const PS::F64 h){
    // M4 Cubic spline kernel
    // (see Table.1 in Dehnen & Aly(2012)[MNRAS,425,1068]
    const PS::F64 u = r/h;
    const PS::F64 phu = std::max(0.0, 0.5 - u);
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 80.0 / (7.0 * math_const::pi * PWR2(h));
#else
    const PS::F64 coeff = 16.0 / (math_const::pi * PWR3(h));
#endif
    return coeff * (PWR3(p1u) - 4.0 * PWR3(phu));
}

PS::F64vec gradW_BSpline_M4(const PS::F64vec dr, const PS::F64 h){
    // Spatial gradient of M4 Cubic spline kernel.
    const PS::F64 r = std::sqrt(dr * dr);
    const PS::F64 u = r/h;
    const PS::F64 phu = std::max(0.0, 0.5-u);
    const PS::F64 p1u = std::max(0.0, 1.0-u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 240.0 / (7.0 * math_const::pi * PWR3(h));
#else
    const PS::F64 coeff = 48.0 / (math_const::pi * PWR4(h));
#endif
    // A term corresponding to Thomas & Couchman (1992) prescription
#if defined(ASURA_FDPS_USE_TC92_PRESCR)
    constexpr PS::F64 a = 1.0/3.0;
    const PS::F64 pau = std::max(0.0, a - u);
    const PS::F64 TC92 = 3.0 * PWR2(pau);
#else
    constexpr PS::F64 TC92 = 0.0;
#endif
    constexpr PS::F64 eps = 1.0e-6;
    return - coeff * (PWR2(p1u) - 4.0 * PWR2(phu) + TC92) * dr / (r + eps * h);
}

PS::F64 dWdh_BSpline_M4(const PS::F64 r, const PS::F64 h){
    // Partial derivative of M4 cubic spline kernel with h
    const PS::F64 u = r/h;
    const PS::F64 phu = std::max(0.0, 0.5 - u);
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 80.0 / (7.0 * math_const::pi * PWR3(h));
    return - coeff * (PWR2(p1u) * (2.0 - 5.0 * u) - 4.0 * PWR2(phu) * (1.0 - 5.0 * u));
#else
    const PS::F64 coeff = 16.0 / (math_const::pi * PWR4(h));
    return - coeff * (PWR2(p1u) * (3.0 - 6.0 * u) - 4.0 * PWR2(phu) * (1.5 - 6.0 * u));
#endif
}

PS::F64 d2Wdrdh_BSpline_M4(const PS::F64 r, const PS::F64 h){
    // Second partial derivative of M4 cubic spline kernel with h and r
    const PS::F64 u = r/h;
    const PS::F64 phu = std::max(0.0, 0.5 - u);
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 80.0 / (7.0 * math_const::pi * PWR4(h));
    return coeff * ( p1u * (9.0 - 15.0 * u)
                   - 4.0 * phu * (4.5 - 15.0 * u));
#else
    const PS::F64 coeff = 96.0 / (math_const::pi * PWR5(h));
    return coeff * ( p1u * (2.0 - 3.0 * u)
                   - 4.0 * phu * (1.0 - 3.0 * u));
#endif
}

PS::F64 d2Wdh2_BSpline_M4(const PS::F64 r, const PS::F64 h){
    // Second partial derivative of M4 cubic spline kernel with h
    const PS::F64 u = r/h;
    const PS::F64 phu = std::max(0.0, 0.5 - u);
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 480.0 / (7.0 * math_const::pi * PWR4(h));
    return coeff * ( p1u * (1.0 + u * (-5.0 + 5.0 * u))
                   - 4.0 * phu * (0.25 + u * (-2.5 + 5.0 * u)) );
#else
    const PS::F64 coeff = 96.0 / (math_const::pi * PWR5(h));
    return coeff * ( p1u * (2.0 + u * (-8.0 + 7.0 * u))
                   - 4.0 * phu * (0.5 + u * (-4.0 + 7.0 * u)) );
#endif
}

PS::F64 W_BSpline_M5(const PS::F64 r, const PS::F64 h) {
    // M5 quartic spline kernel
    const PS::F64 u = r/h;
    constexpr PS::F64 a = 1.0/5.0;
    constexpr PS::F64 b = 3.0/5.0;
    const PS::F64 pau = std::max(0.0, a - u);
    const PS::F64 pbu = std::max(0.0, b - u);
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 46875.0 / (2398.0 * math_const::pi * PWR2(h));
#else
    const PS::F64 coeff = 15625.0 / (512.0 * math_const::pi * PWR3(h));
#endif
    return coeff * (PWR4(p1u) - 5.0 * PWR4(pbu) + 10.0 * PWR4(pau));
}

PS::F64vec gradW_BSpline_M5(const PS::F64vec dr, const PS::F64 h) {
    // Spatial gradient of M5 quartic spline kernel
    const PS::F64 r = std::sqrt(dr * dr);
    const PS::F64 u = r/h;
    constexpr PS::F64 a = 1.0/5.0;
    constexpr PS::F64 b = 3.0/5.0;
    const PS::F64 pau = std::max(0.0, a - u);
    const PS::F64 pbu = std::max(0.0, b - u);
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 187500.0 / (2398.0 * math_const::pi * PWR3(h));
#else
    const PS::F64 coeff = 15625.0 / (128.0 * math_const::pi * PWR4(h));
#endif
    constexpr PS::F64 eps = 1.0e-6;
    return - coeff * (PWR3(p1u) - 5.0 * PWR3(pbu) + 10.0 * PWR3(pau)) * dr / (r + eps * h);
}

PS::F64 dWdh_BSpline_M5(const PS::F64 r, const PS::F64 h) {
    // Partial derivative of M5 quartic spline kernel with h
    const PS::F64 u = r/h;
    constexpr PS::F64 a = 1.0/5.0;
    constexpr PS::F64 b = 3.0/5.0;
    const PS::F64 pau = std::max(0.0, a - u);
    const PS::F64 pbu = std::max(0.0, b - u);
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 46875.0 / (2398.0 * math_const::pi * PWR3(h));
    return - coeff * ( PWR3(p1u) * (2.0 - 6.0 * u)
                     - 5.0 * PWR3(pbu) * (1.2 - 6.0 * u)
                     + 10.0 * PWR3(pau) * (0.4 - 6.0 * u));
#else
    const PS::F64 coeff = 15625.0 / (512.0 * math_const::pi * PWR4(h));
    return - coeff * ( PWR3(p1u) * (3.0 - 7.0 * u)
                     - 5.0 * PWR3(pbu) * (1.8 - 7.0 * u)
                     + 10.0 * PWR3(pau) * (0.6 - 7.0 * u) );
#endif
}

PS::F64 d2Wdrdh_BSpline_M5(const PS::F64 r, const PS::F64 h) {
    // Second partial derivative of M5 quartic spline kernel with h and r
    const PS::F64 u = r/h;
    constexpr PS::F64 a = 1.0/5.0;
    constexpr PS::F64 b = 3.0/5.0;
    const PS::F64 pau = std::max(0.0, a - u);
    const PS::F64 pbu = std::max(0.0, b - u);
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 281250.0 / (1199.0 * math_const::pi * PWR4(h));
    return coeff * ( PWR2(p1u) * (1.0 - 2.0 * u)
                   - 5.0 * PWR2(pbu) * (0.6 - 2.0 * u)
                   + 10.0 * PWR2(pau) * (0.2 - 2.0 * u) );
#else
    const PS::F64 coeff = 15625.0 / (128.0 * math_const::pi * PWR5(h));
    return coeff * ( PWR2(p1u) * (4.0 - 7.0 * u)
                   - 5.0 * PWR2(pbu) * (2.4 - 7.0 * u)
                   + 10.0 * PWR2(pau) * (0.8 - 7.0 * u) );
#endif
}

PS::F64 d2Wdh2_BSpline_M5(const PS::F64 r, const PS::F64 h) {
    // Second partial derivative of M5 quartic spline kernel with h
    const PS::F64 u = r/h;
    constexpr PS::F64 a = 1.0/5.0;
    constexpr PS::F64 b = 3.0/5.0;
    const PS::F64 pau = std::max(0.0, a - u);
    const PS::F64 pbu = std::max(0.0, b - u);
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 281250.0 / (2398.0 * math_const::pi * PWR4(h));
    return coeff * ( PWR2(p1u) * (1.0 + u * (-6.0 + 7.0 * u))
                   - 5.0 * PWR2(pbu) * (0.36 + u * (-3.6 + 7.0 * u))
                   + 10.0 * PWR2(pau) * (0.04 + u * (-1.2 + 7.0 * u)) );
#else
    const PS::F64 coeff = 15625.0 / (128.0 * math_const::pi * PWR5(h));
    return coeff * ( PWR2(p1u) * (3.0 + u * (-14.0 + 14.0 * u))
                   - 5.0 * PWR2(pbu) * (1.08 + u * (-8.4 + 14.0 * u))
                   + 10.0 * PWR2(pau) * (0.12 + u * (-2.8 + 14.0 * u)) );
#endif
}

PS::F64 W_BSpline_M6(const PS::F64 r, const PS::F64 h) {
    // M6 quintic spline kernel
    const PS::F64 u = r/h;
    constexpr PS::F64 a = 1.0/3.0;
    constexpr PS::F64 b = 2.0/3.0;
    const PS::F64 pau = std::max(0.0, a - u);
    const PS::F64 pbu = std::max(0.0, b - u);
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 15309.0 / (478.0 * math_const::pi * PWR2(h));
#else
    const PS::F64 coeff = 2187.0 / (40.0 * math_const::pi * PWR3(h));
#endif
    return coeff * (PWR5(p1u) - 6.0 * PWR5(pbu) + 15.0 * PWR5(pau));
}

PS::F64vec gradW_BSpline_M6(const PS::F64vec dr, const PS::F64 h) {
    // Spatial gradient of M6 quintic spline kernel
    const PS::F64 r = std::sqrt(dr * dr);
    const PS::F64 u = r/h;
    constexpr PS::F64 a = 1.0/3.0;
    constexpr PS::F64 b = 2.0/3.0;
    const PS::F64 pau = std::max(0.0, a - u);
    const PS::F64 pbu = std::max(0.0, b - u);
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 76545.0 / (478.0 * math_const::pi * PWR3(h));
#else
    const PS::F64 coeff = 2187.0 / (8.0 * math_const::pi * PWR4(h));
#endif
    constexpr PS::F64 eps = 1.0e-6;
    return - coeff * (PWR4(p1u) - 6.0 * PWR4(pbu) + 15.0 * PWR4(pau)) * dr / (r + eps * h);
}

PS::F64 dWdh_BSpline_M6(const PS::F64 r, const PS::F64 h) {
    // Partial derivative of M6 quintic spline kernel
    const PS::F64 u = r/h;
    constexpr PS::F64 a = 1.0/3.0;
    constexpr PS::F64 b = 2.0/3.0;
    const PS::F64 pau = std::max(0.0, a - u);
    const PS::F64 pbu = std::max(0.0, b - u);
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 15309.0 / (1434.0 * math_const::pi * PWR3(h));
    return - coeff * ( PWR4(p1u) * (6.0 - 21.0 * u)
                     - 6.0 * PWR4(pbu) * (4.0 - 21.0 * u)
                     + 15.0 * PWR4(pau) * (2.0 - 21.0 * u) );
#else
    const PS::F64 coeff = 2187.0 / (40.0 * math_const::pi * PWR4(h));
    return - coeff * ( PWR4(p1u) * (3.0 - 8.0 * u)
                     - 6.0 * PWR4(pbu) * (2.0 - 8.0 * u)
                     + 15.0 * PWR4(pau) * (1.0 - 8.0 * u) );
#endif
}

PS::F64 d2Wdrdh_BSpline_M6(const PS::F64 r, const PS::F64 h) {
    // Second partial derivative of M6 quintic spline kernel with h and r
    const PS::F64 u = r/h;
    constexpr PS::F64 a = 1.0/3.0;
    constexpr PS::F64 b = 2.0/3.0;
    const PS::F64 pau = std::max(0.0, a - u);
    const PS::F64 pbu = std::max(0.0, b - u);
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 76545.0 / (478.0 * math_const::pi * PWR4(h));
    return coeff * ( PWR3(p1u) * (3.0 - 7.0 * u)
                   - 6.0 * PWR3(pbu) * (2.0 - 7.0 * u)
                   + 15.0 * PWR3(pau) * (1.0 - 7.0 * u) );
#else
    const PS::F64 coeff = 2187.0 / (6.0 * math_const::pi * PWR5(h));
    return coeff * ( PWR3(p1u) * (3.0 - 6.0 * u)
                   - 6.0 * PWR3(pbu) * (2.0 - 6.0 * u)
                   + 15.0 * PWR3(pau) * (1.0 - 6.0 * u) );
#endif
}

PS::F64 d2Wdh2_BSpline_M6(const PS::F64 r, const PS::F64 h) {
    // Second partial derivative of M6 quintic spline kernel with h
    const PS::F64 u = r/h;
    constexpr PS::F64 a = 1.0/3.0;
    constexpr PS::F64 b = 2.0/3.0;
    const PS::F64 pau = std::max(0.0, a - u);
    const PS::F64 pbu = std::max(0.0, b - u);
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 15309.0 / (717.0 * math_const::pi * PWR4(h));
    return coeff * ( PWR3(p1u) * (9.0 + u * (-63.0 + 84.0 * u))
                   - 6.0 * PWR3(pbu) * (4.0 + u * (-42.0 + 84.0 * u))
                   + 15.0 * PWR3(pau) * (1.0 + u * (-21.0 + 84.0 * u)) );
#else
    const PS::F64 coeff = 2187.0 / (30.0 * math_const::pi * PWR5(h));
    return coeff * ( PWR3(p1u) * (9.0 + u * (-48.0 + 54.0 * u))
                   - 6.0 * PWR3(pbu) * (4.0 + u * (-32.0 + 54.0 * u))
                   + 15.0 * PWR3(pau) * (1.0 + u * (-16.0 + 54.0 * u)) );
#endif
}

PS::F64 W_Wendland_C2(const PS::F64 r, const PS::F64 h) {
    // Wendland C2 kernel
    // (see Table.1 in Dehnen & Aly (2012)[MNRAS,425,1068])
    const PS::F64 u = r/h;
    const PS::F64 p1u = std::max(0.0, 1.0 - u); // abbreviation of the positive part of 1-u
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 7.0 / (math_const::pi * PWR2(h));
#else
    const PS::F64 coeff = 21.0 / (2.0 * math_const::pi * PWR3(h));
#endif
    return coeff * PWR4(p1u) * (1.0 + 4.0*u);
}

PS::F64vec gradW_Wendland_C2(const PS::F64vec dr, const PS::F64 h) {
    // Spacial gradient of Wendland C2 kernel.
    const PS::F64 r = std::sqrt(dr * dr);
    const PS::F64 u = r/h;
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 140.0 / (math_const::pi * PWR4(h));
#else
    const PS::F64 coeff = 210.0 / (math_const::pi * PWR5(h));
#endif
    return - coeff * PWR3(p1u) * dr;
}

PS::F64 dWdh_Wendland_C2(const PS::F64 r, const PS::F64 h) {
    // Partial derivative of Wendland C2 kernel with h.
    const PS::F64 u = r/h;
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 14.0 / (math_const::pi * PWR3(h));
    return - coeff * PWR3(p1u) * (1.0 + u * (3.0 - 14.0*u));
#else
    const PS::F64 coeff = 21.0 / (2.0 * math_const::pi * PWR4(h));
    return - coeff * PWR3(p1u) * (3.0 + u * (9.0 - 32.0*u));
#endif
}

PS::F64 d2Wdrdh_Wendland_C2(const PS::F64 r, const PS::F64 h) {
    // Second partial derivative of Wendland C2 kernel with h and r.
    const PS::F64 u = r/h;
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 140.0 / (math_const::pi * PWR4(h));
    return coeff * PWR2(p1u) * u * (4.0 - 7.0 * u);
#else
    const PS::F64 coeff = 210.0 / (math_const::pi * PWR5(h));
    return coeff * PWR2(p1u) * u * (5.0 - 8.0 * u);
#endif
}

PS::F64 d2Wdh2_Wendland_C2(const PS::F64 r, const PS::F64 h) {
    // Second partial derivative of Wendland C2 kernel with h.
    const PS::F64 u = r/h;
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 14.0 / (math_const::pi * PWR4(h));
    return coeff * PWR2(p1u) * (3.0 + u * (6.0 + u * (-91.0 + 112.0 * u)));
#else
    const PS::F64 coeff = 126.0 / (math_const::pi * PWR5(h));
    return coeff * PWR2(p1u) * (1.0 + u * (2.0 + u * (-22.0 + 24.0 * u)));
#endif
}

PS::F64 W_Wendland_C4(const PS::F64 r, const PS::F64 h) {
    // Wendland C4 kernel
    // (see Table.1 in Dehnen & Aly (2012)[MNRAS,425,1068])
    const PS::F64 u = r/h;
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 9.0 / (math_const::pi * PWR2(h));
#else
    const PS::F64 coeff = 495.0 / (32.0 * math_const::pi * PWR3(h));
#endif
    constexpr PS::F64 a2 = 35.0/3.0;
    return coeff * PWR6(p1u) * (1.0 + u * (6.0 + a2 * u));
}

PS::F64vec gradW_Wendland_C4(const PS::F64vec dr, const PS::F64 h) {
    // Spatial gradiaent of Wendland C4 kernel
    const PS::F64 r = std::sqrt(dr * dr);
    const PS::F64 u = r/h;
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 168.0 / (math_const::pi * PWR4(h));
#else
    const PS::F64 coeff = 1155.0 / (4.0 * math_const::pi * PWR5(h));
#endif
    return - coeff * PWR5(p1u) * (1.0 + 5.0 * u) * dr;

}

PS::F64 dWdh_Wendland_C4(const PS::F64 r, const PS::F64 h) {
    // Partial derivative of Wendland C4 kernel with h
    const PS::F64 u = r/h;
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 6.0 / (math_const::pi * PWR3(h));
    return - coeff * PWR5(p1u) * (3.0 + u * (15.0 + u * (-11.0 - 175.0*u)));
#else
    const PS::F64 coeff = 165.0 / (32.0 * math_const::pi * PWR4(h));
    return - coeff * PWR5(p1u) * (9.0 + u * (45.0 + u * (-5.0 - 385.0*u)));
#endif
}

PS::F64 d2Wdrdh_Wendland_C4(const PS::F64 r, const PS::F64 h) {
    // Second partial derivative of Wendland C4 kernel with h and r
    const PS::F64 u = r/h;
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 336.0 / (math_const::pi * PWR4(h));
    return coeff * PWR4(p1u) * u * (2.0 + u * (8.0 - 25.0 * u));
#else
    const PS::F64 coeff = 5775.0 / (4.0 * math_const::pi * PWR5(h));
    return coeff * PWR4(p1u) * u * (1.0 + u * (4.0 - 11.0 * u));
#endif
}

PS::F64 d2Wdh2_Wendland_C4(const PS::F64 r, const PS::F64 h) {
    // Second partial derivative of Wendland C4 kernel with h
    const PS::F64 u = r/h;
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 6.0 / (math_const::pi * PWR4(h));
    return coeff * PWR4(p1u) * (9.0 + u * (36.0 + u * (-190.0 + u * (-940.0 + 1925.0 * u))));
#else
    const PS::F64 coeff = 495.0 / (8.0 * math_const::pi * PWR5(h));
    return coeff * PWR4(p1u) * (3.0 + u * (12.0 + u * (-40.0 + u * (-220.0 + 385.0 * u))));
#endif
}

PS::F64 W_Wendland_C6(const PS::F64 r, const PS::F64 h) {
    // Wendland C6 kernel
    // (see Table.1 in Dehnen & Aly (2012)[MNRAS,425,1068])
    const PS::F64 u = r/h;
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 78.0 / (7.0 * math_const::pi * PWR2(h));
#else
    const PS::F64 coeff = 1365.0 / (64.0 * math_const::pi * PWR3(h));
#endif
    return coeff * PWR8(p1u) * (1.0 + u * (8.0 + u * (25.0 + 32.0*u)));
}

PS::F64vec gradW_Wendland_C6(const PS::F64vec dr, const PS::F64 h) {
    // Spatial gradient of Wendland C6 kernel
    const PS::F64 r = std::sqrt(dr * dr);
    const PS::F64 u = r/h;
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 1716.0 / (7.0 * math_const::pi * PWR4(h));
#else
    const PS::F64 coeff = 15015.0 / (32.0 * math_const::pi * PWR5(h));
#endif
    return - coeff * PWR7(p1u) * (1.0 + u * (7.0 + 16.0*u)) * dr;
}

PS::F64 dWdh_Wendland_C6(const PS::F64 r, const PS::F64 h) {
    // Partial derivative of Wendland C6 kernel
    const PS::F64 u = r/h;
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 156.0 / (7.0 * math_const::pi * PWR3(h));
    return - coeff * PWR7(p1u) * (1.0 + u * (7.0 + u * (6.0 + u * (-70.0 - 208.0*u))));
#else
    const PS::F64 coeff = 1365.0 / (64.0 * math_const::pi * PWR4(h));
    return - coeff * PWR7(p1u) * (3.0 + u * (21.0 + u * (29.0 + u * (-133.0 - 448.0*u))));
#endif
}

PS::F64 d2Wdrdh_Wendland_C6(const PS::F64 r, const PS::F64 h) {
    // Second partial derivative of Wendland C6 kernel with h and r
    const PS::F64 u = r/h;
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 6864.0 / (7.0 * math_const::pi * PWR4(h));
    return coeff * PWR6(p1u) * u * (1.0 + u * (6.0 + u * (3.0 - 52.0 * u)));
#else
    const PS::F64 coeff = 15015.0 / (32.0 * math_const::pi * PWR5(h));
    return coeff * PWR6(p1u) * u * (5.0 + u * (30.0 + u * (21.0 - 224.0 * u)));
#endif
}

PS::F64 d2Wdh2_Wendland_C6(const PS::F64 r, const PS::F64 h) {
    // Second partial derivative of Wendland C6 kernel with h
    const PS::F64 u = r/h;
    const PS::F64 p1u = std::max(0.0, 1.0 - u);
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 coeff = 156.0 / (7.0 * math_const::pi * PWR4(h));
    return coeff * PWR6(p1u) * (3.0 + u * (18.0 + u * (-47.0 + u * (-492.0 + u * (-546.0 + 2912.0 * u)))));
#else
    const PS::F64 coeff = 4095.0 / (32.0 * math_const::pi * PWR5(h));
    return coeff * PWR6(p1u) * (2.0 + u * (12.0 + u * (-13.0 + u * (-218.0 + u * (-287.0 + 1120.0 * u)))));
#endif
}

void outputToFile(const std::string filename,
                  const std::vector<PS::F64vec> dr,
                  const PS::F64 h,
                  PS::F64 (*pfunc_W)(const PS::F64, const PS::F64),
                  PS::F64vec (*pfunc_gradW)(const PS::F64vec, const PS::F64),
                  PS::F64 (*pfunc_dWdh)(const PS::F64, const PS::F64),
                  PS::F64 (*pfunc_d2Wdrdh)(const PS::F64, const PS::F64),
                  PS::F64 (*pfunc_d2Wdh2)(const PS::F64, const PS::F64)) {
    std::ofstream ofs;
    ofs.open(filename.c_str(), std::ios::trunc);
    for (PS::S32 i = 0; i < dr.size(); i++) {
        const PS::F64 r = std::sqrt(dr[i] * dr[i]);
        ofs << r << "   "
            << pfunc_W(r, h) << "   "
            << dr[i] * pfunc_gradW(dr[i], h) << "    "
            << pfunc_dWdh(r, h) << "    "
            << pfunc_d2Wdrdh(r, h) << "    "
            << pfunc_d2Wdh2(r, h)
            << std::endl; 
    }
    ofs.close();
}

void outputSPHKernelProfile(const PS::S32 root) {
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    assert(0 <= root && root < n_proc);
    const PS::S32 my_rank = PS::Comm::getRank();
    if (my_rank != root) return;

    // Define sampling points
    constexpr PS::S32 N = 256;
    constexpr PS::F64 dx = 3.0 / N;
    std::vector<PS::F64vec> dr(N);
    for (PS::S32 i = 0; i < N; i++) {
        dr[i] = 0.0; // reset
        dr[i].x = dx * (0.5 + i);
    }

    // Output M4 Cubic spline kernel
    // In the following, h is set to the H/h values in Table.1
    // of Dehnen & Alley (2012)[MNRAS,425,1068] so that we can
    // compare the output data with their Figure 1 directly.
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    outputToFile("W_M4.txt", dr, 1.778002,
                 W_BSpline_M4, gradW_BSpline_M4, dWdh_BSpline_M4,
                 d2Wdrdh_BSpline_M4, d2Wdh2_BSpline_M4);
#else
    outputToFile("W_M4.txt", dr, 1.825742,
                 W_BSpline_M4, gradW_BSpline_M4, dWdh_BSpline_M4,
                 d2Wdrdh_BSpline_M4, d2Wdh2_BSpline_M4);
#endif

    // Output M5 Quatic spline kernel
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    outputToFile("W_M5.txt", dr, 1.977173,
                 W_BSpline_M5, gradW_BSpline_M5, dWdh_BSpline_M5,
                 d2Wdrdh_BSpline_M5, d2Wdh2_BSpline_M5);
#else
    outputToFile("W_M5.txt", dr, 2.018932,
                 W_BSpline_M5, gradW_BSpline_M5, dWdh_BSpline_M5,
                 d2Wdrdh_BSpline_M5, d2Wdh2_BSpline_M5);
#endif    

    // Output M6 Quintic spline kernel
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    outputToFile("W_M6.txt", dr, 2.158131,
                 W_BSpline_M6, gradW_BSpline_M6, dWdh_BSpline_M6,
                 d2Wdrdh_BSpline_M6, d2Wdh2_BSpline_M6);
#else
    outputToFile("W_M6.txt", dr, 2.195775,
                 W_BSpline_M6, gradW_BSpline_M6, dWdh_BSpline_M6,
                 d2Wdrdh_BSpline_M6, d2Wdh2_BSpline_M6);
#endif

    // Output Wendland C2 kernel
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    outputToFile("W_C2.txt", dr, 1.897367,
                 W_Wendland_C2, gradW_Wendland_C2, dWdh_Wendland_C2,
                 d2Wdrdh_Wendland_C2, d2Wdh2_Wendland_C2);
#else
    outputToFile("W_C2.txt", dr, 1.936492,
                 W_Wendland_C2, gradW_Wendland_C2, dWdh_Wendland_C2,
                 d2Wdrdh_Wendland_C2, d2Wdh2_Wendland_C2);
#endif

    // Output Wendland C4 kernel
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    outputToFile("W_C4.txt", dr, 2.171239,
                 W_Wendland_C4, gradW_Wendland_C4, dWdh_Wendland_C4,
                 d2Wdrdh_Wendland_C4, d2Wdh2_Wendland_C4);
#else
    outputToFile("W_C4.txt", dr, 2.207940,
                 W_Wendland_C4, gradW_Wendland_C4, dWdh_Wendland_C4,
                 d2Wdrdh_Wendland_C4, d2Wdh2_Wendland_C4);
#endif

    // Output Wendland C6 kernel
#if defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    outputToFile("W_C6.txt", dr, 2.415230,
                 W_Wendland_C6, gradW_Wendland_C6, dWdh_Wendland_C6,
                 d2Wdrdh_Wendland_C6, d2Wdh2_Wendland_C6);
#else
    outputToFile("W_C6.txt", dr, 2.449490,
                 W_Wendland_C6, gradW_Wendland_C6, dWdh_Wendland_C6,
                 d2Wdrdh_Wendland_C6, d2Wdh2_Wendland_C6);
#endif

    // Output gnuplot script
    const std::string filename = "plot_SPH_kernels.gp";
    std::ofstream ofs;
    ofs.open(filename.c_str(), std::ios::trunc);
    ofs << "set terminal postscript enhanced color eps" << std::endl;
    ofs << "set size square" << std::endl;
    ofs << "set grid" << std::endl;
    ofs << "set xrange [0.0:2.5]" << std::endl;
    ofs << "set xlabel \"r\"" << std::endl;
    ofs << "set ylabel \"W\"" << std::endl;
    ofs << "set key right top" << std::endl;
    ofs << "set output \"W.eps\"" << std::endl;
    ofs << "plot \"W_M4.txt\" u 1:2 w l lw 2 t \"M_{4} cubic spline\", \\" << std::endl;
    ofs << "     \"W_M5.txt\" u 1:2 w l lw 2 t \"M_{5} quartic spline\", \\" << std::endl;
    ofs << "     \"W_M6.txt\" u 1:2 w l lw 2 t \"M_{6} quintic spline\", \\" << std::endl;
    ofs << "     \"W_C2.txt\" u 1:2 w l lw 2 t \"Wendland C^{2}\", \\" << std::endl;
    ofs << "     \"W_C4.txt\" u 1:2 w l lw 2 t \"Wendland C^{4}\", \\" << std::endl;
    ofs << "     \"W_C6.txt\" u 1:2 w l lw 2 t \"Wendland C^{6}\"" << std::endl;
    ofs << "set ylabel \"r * dWdr\"" << std::endl;
    ofs << "set key right bottom" << std::endl;
    ofs << "set output \"rdWdr.eps\"" << std::endl;
    ofs << "plot \"W_M4.txt\" u 1:3 w l lw 2 t \"M_{4} cubic spline\", \\" << std::endl;
    ofs << "     \"W_M5.txt\" u 1:3 w l lw 2 t \"M_{5} quartic spline\", \\" << std::endl;
    ofs << "     \"W_M6.txt\" u 1:3 w l lw 2 t \"M_{6} quintic spline\", \\" << std::endl;
    ofs << "     \"W_C2.txt\" u 1:3 w l lw 2 t \"Wendland C^{2}\", \\" << std::endl;
    ofs << "     \"W_C4.txt\" u 1:3 w l lw 2 t \"Wendland C^{4}\", \\" << std::endl;
    ofs << "     \"W_C6.txt\" u 1:3 w l lw 2 t \"Wendland C^{6}\"" << std::endl;
    ofs << "set ylabel \"dWdh\"" << std::endl;
    ofs << "set key right bottom" << std::endl;
    ofs << "set output \"dWdh.eps\"" << std::endl;
    ofs << "plot \"W_M4.txt\" u 1:4 w l lw 2 t \"M_{4} cubic spline\", \\" << std::endl;
    ofs << "     \"W_M5.txt\" u 1:4 w l lw 2 t \"M_{5} quartic spline\", \\" << std::endl;
    ofs << "     \"W_M6.txt\" u 1:4 w l lw 2 t \"M_{6} quintic spline\", \\" << std::endl;
    ofs << "     \"W_C2.txt\" u 1:4 w l lw 2 t \"Wendland C^{2}\", \\" << std::endl;
    ofs << "     \"W_C4.txt\" u 1:4 w l lw 2 t \"Wendland C^{4}\", \\" << std::endl;
    ofs << "     \"W_C6.txt\" u 1:4 w l lw 2 t \"Wendland C^{6}\"" << std::endl;
    ofs << "set ylabel \"d2Wdrdh\"" << std::endl;
    ofs << "set key right top" << std::endl;
    ofs << "set output \"d2Wdrdh.eps\"" << std::endl;
    ofs << "plot \"W_M4.txt\" u 1:5 w l lw 2 t \"M_{4} cubic spline\", \\" << std::endl;
    ofs << "     \"W_M5.txt\" u 1:5 w l lw 2 t \"M_{5} quartic spline\", \\" << std::endl;
    ofs << "     \"W_M6.txt\" u 1:5 w l lw 2 t \"M_{6} quintic spline\", \\" << std::endl;
    ofs << "     \"W_C2.txt\" u 1:5 w l lw 2 t \"Wendland C^{2}\", \\" << std::endl;
    ofs << "     \"W_C4.txt\" u 1:5 w l lw 2 t \"Wendland C^{4}\", \\" << std::endl;
    ofs << "     \"W_C6.txt\" u 1:5 w l lw 2 t \"Wendland C^{6}\"" << std::endl;
    ofs << "set ylabel \"d2Wdh2\"" << std::endl;
    ofs << "set key right top" << std::endl;
    ofs << "set output \"d2Wdh2.eps\"" << std::endl;
    ofs << "plot \"W_M4.txt\" u 1:6 w l lw 2 t \"M_{4} cubic spline\", \\" << std::endl;
    ofs << "     \"W_M5.txt\" u 1:6 w l lw 2 t \"M_{5} quartic spline\", \\" << std::endl;
    ofs << "     \"W_M6.txt\" u 1:6 w l lw 2 t \"M_{6} quintic spline\", \\" << std::endl;
    ofs << "     \"W_C2.txt\" u 1:6 w l lw 2 t \"Wendland C^{2}\", \\" << std::endl;
    ofs << "     \"W_C4.txt\" u 1:6 w l lw 2 t \"Wendland C^{4}\", \\" << std::endl;
    ofs << "     \"W_C6.txt\" u 1:6 w l lw 2 t \"Wendland C^{6}\"" << std::endl;
    ofs.close();

}
