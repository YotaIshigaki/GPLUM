/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "mathematical_constants.h"
#include "SPH_kernel.h"

/* Kernel Function */
PS::F64 W(const PS::F64 r, const PS::F64 h){
    // M4 Cubic spline kernel
    // (see Eq. (4) in Springel (2005)[MNRAS,364,1105])
    const PS::F64 u = r/h;
    const PS::F64 cc=8.0/(math_const::pi*h*h*h);
    if (u <= 0.5) {
        const PS::F64 u2 = u*u;
        return cc*(1.0+u2*6.0*(-1.0+u));
    }
    else if ((0.5 < u) && (u <= 1.0)) {
        const PS::F64 s = 1.0-u;
        return cc*2.0*s*s*s;
    }
    else {
        return 0.0;
    }
}

PS::F64vec gradW(const PS::F64vec dr, const PS::F64 h){
    // This subroutine gives \nabla W(r,h), i.e.,
    // \dfrac{\partial W(r,h)}{\partial r}\dfrac{dr}{r}.
    const PS::F64 r = std::sqrt(dr * dr);
    const PS::F64 u=r/h;
    const PS::F64 cc = -48.0/(math_const::pi*h*h*h*h);
#if defined(USE_PRESCR_OF_THOMAS_COUCHMAN_1992)
    if (u <= 1.0/3.0) {
        return dr * cc*(1.0/3.0)/(r);
    }
    else if ((1.0/3.0 < u) && (u <= 0.5)) {
        return dr * cc*u*(2.0-3.0*u)/(r);
    }
    else if ((0.5 < u) && (u < 1.0)) {
        return dr * cc*(1.0-u)*(1.0-u)/(r);
    }
    else {
        return PS::F64vec(0.0, 0.0, 0.0);
    }
#else
    if ((0.0 < u) && (u <= 0.5)) {
        return dr * cc*u*(2.0-3.0*u)/(r);
    }
    else if ((0.5 < u) && (u < 1.0)) {
        return dr * cc*(1.0-u)*(1.0-u)/(r);
    }
    else {
        // r=0 case is included in this branch
        return PS::F64vec(0.0, 0.0, 0.0);
    }
#endif
}

PS::F64 dWdh(const PS::F64 r, const PS::F64 h){
   // This subroutine gives dW(r,h)/dh, i.e.,
   // \dfrac{\partial W(r,h)}{\partial h}.
   const PS::F64 u=r/h;
   const PS::F64 cc=-24.0/(math_const::pi*h*h*h*h);
   if (u <= 0.5) {
      const PS::F64 u2 = u*u;
      return cc*(1.0
                +u2*(-10.0
                     +12.0*u));
   }
   else if ((0.5 < u) && (u < 1.0)) {
      const PS::F64 s = 1.0-u;
      return cc*2.0*s*s*(1.0-2.0*u);
   }
   else {
      return 0.0;
   }

}
