#pragma once
/* C++ headers */
#include "common.h"
/* FDPS header */
#include <particle_simulator.hpp>
/* User-defined headers */

class HernquistModel {
private:
    PS::F64 M,rs,rt;

public:
    HernquistModel(const PS::F64 M, const PS::F64 rs, const PS::F64 rt);
    PS::F64 getDensity(const PS::F64 r) const;
    PS::F64 getDensityCutoffed(const PS::F64 r) const;
    PS::F64 getDensityXYZ(const PS::F64 x,
                          const PS::F64 y,
                          const PS::F64 z) const;
    PS::F64 getDensityCutoffedXYZ(const PS::F64 x,
                                  const PS::F64 y,
                                  const PS::F64 z) const;
    PS::F64 getEnclosedMass(const PS::F64 r) const;
    PS::F64 getRadiusFromMass(const PS::F64 Mr) const;
    PS::F64 getPotential(const PS::F64 r) const;
    PS::F64 getUnnormalizedDF(const PS::F64 E) const;
    PS::F64 getDF(const PS::F64 E) const;
    PS::F64 getVelocity(const PS::F64 phi) const;
};


class MiyamotoNagaiModel {
private:
    PS::F64 M,Rs,zs,Rt,zt;
public:

    MiyamotoNagaiModel(const PS::F64 M, const PS::F64 Rs, const PS::F64 zs,
                       const PS::F64 Rt, const PS::F64 zt);
    PS::F64 getDensity(const PS::F64 R, const PS::F64 z) const;
    PS::F64 getDensityCutoffed(const PS::F64 R, const PS::F64 z) const;
    PS::F64 getDensityXYZ(const PS::F64 x, const PS::F64 y, const PS::F64 z) const;
    PS::F64 getDensityCutoffedXYZ(const PS::F64 x, const PS::F64 y, const PS::F64 z) const;
    PS::F64 getPotential(const PS::F64 R, const PS::F64 z) const;
};

class ExponentialDiskModel {
private:
    PS::F64 M,Rs,zs,Rt,zt;
    PS::F64 rho0;
public:

    ExponentialDiskModel(const PS::F64 M, const PS::F64 Rs, const PS::F64 zs,
                         const PS::F64 Rt, const PS::F64 zt);
    PS::F64 getDensity(const PS::F64 R, const PS::F64 z) const;
    PS::F64 getDensityCutoffed(const PS::F64 R, const PS::F64 z) const;
    PS::F64 getDensityXYZ(const PS::F64 x, const PS::F64 y, const PS::F64 z) const;
    PS::F64 getDensityCutoffedXYZ(const PS::F64 x, const PS::F64 y, const PS::F64 z) const;
    PS::F64 getMassCutoffed() const;
};
