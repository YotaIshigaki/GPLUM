#pragma once
/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "user_defined.h"

#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_NORMAL_RUN
void readTipsyFile(std::string& input_file_name,
                   PS::ParticleSystem<FP_dm>& psys);

void GalaxyIC(PS::ParticleSystem<FP_dm>& psys_dm,
              PS::ParticleSystem<FP_star>& psys_star,
              PS::ParticleSystem<FP_gas>& psys_gas);

PS::F64ort getDomainInSphericalCoordinate(const PS::F64ort box_cart);

PS::F64ort getDomainInCylindricalCoordinate(const PS::F64ort box_cart);

void GordonBellIC(PS::ParticleSystem<FP_dm>& psys_dm,
                  PS::ParticleSystem<FP_star>& psys_star,
                  PS::ParticleSystem<FP_gas>& psys_gas,
                  PS::DomainInfo & dinfo);
#endif

#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_GLASS_DATA_GENERATION_MODE
void GlassDataGenerationModeIC(PS::ParticleSystem<FP_dm>& psys_dm,
                               PS::ParticleSystem<FP_star>& psys_star,
                               PS::ParticleSystem<FP_gas>& psys_gas);
#endif

#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_PARTICLE_COLLISION_TEST
void ParticleCollisionTestIC(PS::ParticleSystem<FP_dm>& psys_dm,
                             PS::ParticleSystem<FP_star>& psys_star,
                             PS::ParticleSystem<FP_gas>& psys_gas);
#endif

#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_SHOCK_TUBE_TEST
void ShockTubeTestIC(PS::ParticleSystem<FP_dm>& psys_dm,
                     PS::ParticleSystem<FP_star>& psys_star,
                     PS::ParticleSystem<FP_gas>& psys_gas);
#endif

#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_SURFACE_TENSION_TEST
void SurfaceTensionTestIC(PS::ParticleSystem<FP_dm>& psys_dm,
                          PS::ParticleSystem<FP_star>& psys_star,
                          PS::ParticleSystem<FP_gas>& psys_gas);
#endif

#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_POINT_EXPLOSION_TEST
void PointExplosionTestIC(PS::ParticleSystem<FP_dm>& psys_dm,
                          PS::ParticleSystem<FP_star>& psys_star,
                          PS::ParticleSystem<FP_gas>& psys_gas);
#endif
