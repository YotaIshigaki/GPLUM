#pragma once
/* C++ headers */
#include "common.h"
/* FDPS header */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "user_defined.h"

void predictFeedbackRadius(PS::ParticleSystem<FP_star> & psys_star,
                           PS::DomainInfo & dinfo,
                           PS::TreeForForceShort<Force_hydro, EP_hydro, EP_hydro>::Symmetry & tree,
                           const PS::F64 t, const PS::F64 dt);

void calcKernelSize(PS::ParticleSystem<FP_gas> & psys_gas,
                    PS::ParticleSystem<FP_star> & psys_star,
                    PS::DomainInfo & dinfo, 
                    PS::TreeForForceShort<Force_knl_sz, EPI_knl_sz, EPJ_knl_sz>::Gather & tree,
                    const PS::F64 t, const PS::F64 dt,
                    const FunctionUseType use_type);

void calcDensityAndPressure(PS::ParticleSystem<FP_gas> & psys_gas,
                            PS::ParticleSystem<FP_star> & psys_star,
                            PS::DomainInfo & dinfo, 
                            PS::TreeForForceShort<Force_dens, EP_hydro, EP_hydro>::Gather & tree,
                            const PS::F64 t, const PS::F64 dt,
                            const FunctionUseType use_type);

void calcSoundSpeed(PS::ParticleSystem<FP_gas>& psys);

void calcArtificialViscosity(PS::ParticleSystem<FP_gas>& psys,
                             const PS::F64 dt);

void checkDensityFluctuation(PS::ParticleSystem<FP_gas>& psys);
