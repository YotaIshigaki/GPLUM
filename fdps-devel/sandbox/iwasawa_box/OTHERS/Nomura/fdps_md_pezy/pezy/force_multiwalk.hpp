#pragma once
#include <iostream>
#include <particle_simulator.hpp>

#ifdef ENABLE_PEZY
#include<pzcl/pzcl_ocl_wrapper.h>
#include<PZSDKHelper.h>
#endif

#include "./../main.hpp"

void InitializeDevice();

PS::S32 DispatchKernelWithSP(const PS::S32 tag,
                             const PS::S32 n_walk,
                             const EPILJ ** epi,
                             const PS::S32 *  n_epi,
                             const EPJLJ ** epj,
                             const PS::S32 *  n_epj
			     );

PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 * ni,
                       ForceLJ      ** force);


