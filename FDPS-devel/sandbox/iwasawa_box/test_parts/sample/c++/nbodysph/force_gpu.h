#pragma once
#include <iostream>
#include <particle_simulator.hpp>
#include <helper_cuda.h>
#include "header.h"

PS::S32 CalcDensDispatch(const PS::S32 tag,
                         const PS::S32 n_walk,
                         const EPI::Dens ** epi,
                         const PS::S32 *  n_epi,
                         const PS::S32 ** id_epj,
                         const PS::S32 *  n_epj,
                         const EPJ::Dens * epj,
                         const PS::S32 n_epj_tot,
                         const bool send_flag);

PS::S32 CalcDrvtDispatch(const PS::S32 tag,
                         const PS::S32 n_walk,
                         const EPI::Drvt ** epi,
                         const PS::S32 *  n_epi,
                         const PS::S32 ** id_epj,
                         const PS::S32 *  n_epj,
                         const EPJ::Drvt * epj,
                         const PS::S32 n_epj_tot,
                         const bool send_flag);

PS::S32 CalcHydroDispatch(const PS::S32 tag,
                          const PS::S32 n_walk,
                          const EPI::Hydro ** epi,
                          const PS::S32 *  n_epi,
                          const PS::S32 ** id_epj,
                          const PS::S32 *  n_epj,
                          const EPJ::Hydro * epj,
                          const PS::S32 n_epj_tot,
                          const bool send_flag);

PS::S32 CalcGravDispatch(const PS::S32 tag,
                         const PS::S32 n_walk,
                         const EPI::Grav ** epi,
                         const PS::S32 *  n_epi,
                         const PS::S32 ** id_epj,
                         const PS::S32 *  n_epj,
                         const PS::S32 ** id_spj,
                         const PS::S32 *  n_spj,
                         const EPJ::Grav * epj,
                         const PS::S32 n_epj_tot,
                         const PS::SPJMonopole * spj,
                         const PS::S32 n_spj_tot,
                         const bool send_flag);

template<class Tforce>
PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 * ni,
                       Tforce      ** force);



