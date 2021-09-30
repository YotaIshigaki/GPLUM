/* C++ headers */
#include "common.h"
/* FDPS header */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "mathematical_constants.h"
#include "physical_constants.h"
#include "run_parameters.h"
#include "density_distribution_models.h"

//*** Hernquist model
HernquistModel::HernquistModel(const PS::F64 M,
                               const PS::F64 rs,
                               const PS::F64 rt) {
    this->M = M;
    this->rs = rs;
    this->rt = rt;
}

PS::F64 HernquistModel::getDensity(const PS::F64 r) const {
    return  (this->M * this->rs)
          / (2.0 * math_const::pi * r * CUBE(r + this->rs) );
}

PS::F64 HernquistModel::getDensityCutoffed(const PS::F64 r) const {
    if (r <= this->rt) return this->getDensity(r);
    else return 0.0;
}

PS::F64 HernquistModel::getDensityXYZ(const PS::F64 x,
                                      const PS::F64 y,
                                      const PS::F64 z) const {
    const PS::F64 r = std::sqrt(x*x + y*y + z*z);
    return this->getDensity(r);
}

PS::F64 HernquistModel::getDensityCutoffedXYZ(const PS::F64 x,
                                              const PS::F64 y,
                                              const PS::F64 z) const {
    const PS::F64 r = std::sqrt(x*x + y*y + z*z);
    if (r <= this->rt) return this->getDensity(r);
    else return 0.0;
}

PS::F64 HernquistModel::getEnclosedMass(const PS::F64 r) const {
    return this->M * SQ(r)/SQ(r + this->rs);
}

PS::F64 HernquistModel::getRadiusFromMass(const PS::F64 Menc) const {
    return  this->rs * (Menc + std::sqrt(Menc * this->M))
          / (this->M - Menc);
}

PS::F64 HernquistModel::getPotential(const PS::F64 r) const {
    return - (phys_const::Ggrav * this->M)/(r + this->rs);
}

PS::F64 HernquistModel::getUnnormalizedDF(const PS::F64 E) const {
    const PS::F64 q = std::sqrt(- this->rs * E/(phys_const::Ggrav * this->M));
    const PS::F64 p = 1.0 - SQ(q);
    return std::pow(p, -2.5) * (3*std::asin(q) + q * std::sqrt(p) * (1.0 - 2.0 * SQ(q)) * (8.0*PWR4(q) - 8.0*SQ(q) - 3.0));  
}

PS::F64 HernquistModel::getDF(const PS::F64 E) const {
    const PS::F64 vg = std::sqrt(phys_const::Ggrav * this->M / this->rs);
    const PS::F64 coeff = this->M / (8.0 * std::sqrt(2.0) * CUBE(math_const::pi * this->rs * vg));
    return coeff * this->getUnnormalizedDF(E);
}

PS::F64 HernquistModel::getVelocity(const PS::F64 phi) const {
    const PS::F64 vmax = std::sqrt(-2.0 * phi);
    const PS::F64 fmax = this->getUnnormalizedDF(phi) * SQ(vmax);
    // Note that we use DF * SQ(v) as the probabilify distribution of v
    // in order to prevent the resultant distibution from diverging at v=0.
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    PS::F64 v; 
    for (;;) {
        v = vmax * dist(run_param::prng::mt);
        const PS::F64 f = fmax * dist(run_param::prng::mt);
        const PS::F64 fv = this->getUnnormalizedDF(phi + 0.5*SQ(v)) * SQ(v);
        if (f < fv) break;
    } 
    return v;
}


//*** Miyamoto-Nagai model
MiyamotoNagaiModel::MiyamotoNagaiModel(const PS::F64 M,
                                       const PS::F64 Rs,
                                       const PS::F64 zs,
                                       const PS::F64 Rt,
                                       const PS::F64 zt) {
    this->M  = M;
    this->Rs = Rs;
    this->zs = zs;
    this->Rt = Rt;
    this->zt = zt;
}

PS::F64 MiyamotoNagaiModel::getDensity(const PS::F64 R,
                                       const PS::F64 z) const {
    const PS::F64 q = std::sqrt(SQ(z) + SQ(this->zs));
    return (SQ(this->zs) * this->M)/(4.0 * math_const::pi)
          * (this->Rs * SQ(R) + (this->Rs + 3.0*q) * SQ(this->Rs + q))
          / (std::pow(SQ(R)+SQ(this->Rs + q), 2.5) * CUBE(q));
}

PS::F64 MiyamotoNagaiModel::getDensityCutoffed(const PS::F64 R,
                                               const PS::F64 z) const {
    if (R <= this->Rt && std::fabs(z) <= this->zt) {
        return this->getDensity(R,z);
    } else {
        return 0.0;
    }
}

PS::F64 MiyamotoNagaiModel::getDensityXYZ(const PS::F64 x,
                                          const PS::F64 y,
                                          const PS::F64 z) const {
    const PS::F64 R = std::sqrt(x*x + y*y);
    return this->getDensity(R,z);
}

PS::F64 MiyamotoNagaiModel::getDensityCutoffedXYZ(const PS::F64 x,
                                                  const PS::F64 y,
                                                  const PS::F64 z) const {
    const PS::F64 R = std::sqrt(x*x + y*y);
    return this->getDensityCutoffed(R,z);
}

PS::F64 MiyamotoNagaiModel::getPotential(const PS::F64 R,
                                         const PS::F64 z) const {
    return - (phys_const::Ggrav * this->M)
           / (SQ(R) + SQ(this->Rs + std::sqrt(SQ(z) + SQ(this->zs))));
}


//*** Exponential disk model
ExponentialDiskModel::ExponentialDiskModel(const PS::F64 M,
                                           const PS::F64 Rs,
                                           const PS::F64 zs,
                                           const PS::F64 Rt,
                                           const PS::F64 zt) {
    this->M  = M;
    this->Rs = Rs;
    this->zs = zs;
    this->Rt = Rt;
    this->zt = zt;
    this->rho0 = M / (4.0 * math_const::pi * Rs * Rs * zs);
}

PS::F64 ExponentialDiskModel::getDensity(const PS::F64 R, const PS::F64 z) const {
    return this->rho0 * std::exp(-R/this->Rs) * std::exp(-std::abs(z)/this->zs);
}

PS::F64 ExponentialDiskModel::getDensityCutoffed(const PS::F64 R, const PS::F64 z) const {
    if (R <= this->Rt && std::fabs(z) <= this->zt) return this->getDensity(R,z);
    else return 0.0;
}

PS::F64 ExponentialDiskModel::getDensityXYZ(const PS::F64 x,
                                            const PS::F64 y,
                                            const PS::F64 z) const {
    const PS::F64 R = std::sqrt(x*x + y*y);
    return this->getDensity(R,z);
}


PS::F64 ExponentialDiskModel::getDensityCutoffedXYZ(const PS::F64 x,
                                                    const PS::F64 y,
                                                    const PS::F64 z) const {
    const PS::F64 R = std::sqrt(x*x + y*y);
    return this->getDensityCutoffed(R,z);
}

PS::F64 ExponentialDiskModel::getMassCutoffed() const {
    return this->M 
          * (1.0 - (1.0 + this->Rt/this->Rs) * std::exp(-this->Rt/this->Rs))
          * (1.0 - std::exp(-this->zt/this->zs));
            
}
