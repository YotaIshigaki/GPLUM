#pragma once
#include <iostream>
#include <particle_simulator.hpp>

#ifdef ENABLE_GPU_CUDA
#include <helper_cuda.h>
#elif ENABLE_PEZY
#include<pzcl/pzcl_ocl_wrapper.h>
#include<PZSDKHelper.h>
#endif

#include "user-defined.hpp"

void InitializeDevice();

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


