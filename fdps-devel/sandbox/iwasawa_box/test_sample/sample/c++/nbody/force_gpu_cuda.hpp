#pragma once
#include <iostream>
#include <particle_simulator.hpp>
#include <helper_cuda.h>
#include "user-defined.hpp"

extern PS::F64 WTIME_KERNEL;
extern PS::F64 WTIME_SEND_EPJ;
extern PS::F64 WTIME_H2D;
extern PS::F64 WTIME_D2H;

PS::S32 DispatchKernelWithSP(const PS::S32 tag,
                             const PS::S32 n_walk,
                             const FPGrav ** epi,
                             const PS::S32 *  n_epi,
                             const FPGrav ** epj,
                             const PS::S32 *  n_epj,
                             const PS::SPJMonopole ** spj,
                             const PS::S32  * n_spj);

PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 * ni,
                       FPGrav      ** force);


PS::S32 DispatchKernelIndex(const PS::S32 tag,
                            const PS::S32 n_walk,
                            const FPGrav ** epi,
                            const PS::S32 *  n_epi,
                            const PS::S32 ** id_epj,
                            const PS::S32 *  n_epj,
                            const PS::S32 ** id_spj,
                            const PS::S32 *  n_spj,
                            const FPGrav * epj,
                            const PS::S32 n_epj_tot,
                            const PS::SPJMonopole * spj,
                            const PS::S32 n_spj_tot,
                            const bool send_flag);


