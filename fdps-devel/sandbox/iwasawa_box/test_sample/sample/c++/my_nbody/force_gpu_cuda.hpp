#pragma once
#include <iostream>
#include <particle_simulator.hpp>
#include <helper_cuda.h>
#include "user-defined.hpp"

extern PS::F64 WTIME_KERNEL;
extern PS::F64 WTIME_SEND_IP;
extern PS::F64 WTIME_COPY_IP;
extern PS::F64 WTIME_SEND_ID;
extern PS::F64 WTIME_COPY_ID;
extern PS::F64 WTIME_SEND_JP;
extern PS::F64 WTIME_COPY_JP;
extern PS::F64 WTIME_RECV_FORCE;
extern PS::F64 WTIME_COPY_FORCE;
extern PS::F64 WTIME_COPY_IP_ID;
extern PS::F64 WTIME_COPY_IP_JP;

extern PS::F64 WTIME_H2D_IP;
extern PS::F64 WTIME_H2D_LIST;
extern PS::F64 WTIME_D2H_FORCE;
extern PS::F64 WTIME_H2D_ALL_PTCL;

PS::S32 DispatchKernelWithSP(const PS::S32 tag,
                             const PS::S32 n_walk,
                             const EPI ** epi,
                             const PS::S32 *  n_epi,
                             const EPJ ** epj,
                             const PS::S32 *  n_epj,
                             const PS::SPJMonopole ** spj,
                             const PS::S32  * n_spj);

PS::S32 DispatchKernelIndex(const PS::S32 tag,
                            const PS::S32 n_walk,
                            const EPI ** epi,
                            const PS::S32 *  n_epi,
                            const PS::S32 ** id_epj,
                            const PS::S32 *  n_epj,
                            const PS::S32 ** id_spj,
                            const PS::S32 *  n_spj,
                            const EPJ * epj,
                            const PS::S32 n_epj_tot,
                            const PS::SPJMonopole * spj,
                            const PS::S32 n_spj_tot,
                            const bool send_flag);

PS::S32 DispatchKernelIndex2(const PS::S32 tag,
                             const PS::S32 n_walk,
                             const EPIGrav ** epi,
                             const PS::S32 *  n_epi,
                             const PS::S32 ** id_epj,
                             const PS::S32 *  n_epj,
                             const PS::S32 ** id_spj,
                             const PS::S32 *  n_spj,
                             const EPJGrav * epj,
                             const PS::S32 n_epj_tot,
                             const PS::SPJMonopole * spj,
                             const PS::S32 n_spj_tot,
                             const bool send_flag);

PS::S32 DispatchKernelIndexStream(const PS::S32 tag,
                                  const PS::S32 n_walk,
                                  const EPI ** epi,
                                  const PS::S32 *  n_epi,
                                  const PS::S32 ** id_epj,
                                  const PS::S32 *  n_epj,
                                  const PS::S32 ** id_spj,
                                  const PS::S32 *  n_spj,
                                  const EPJ * epj,
                                  const PS::S32 n_epj_tot,
                                  const PS::SPJMonopole * spj,
                                  const PS::S32 n_spj_tot,
                                  const bool send_flag);

PS::S32 DispatchKernelIndexStream2(const PS::S32 tag,
                                   const PS::S32 n_walk,
                                   const EPI ** epi,
                                   const PS::S32 *  n_epi,
                                   const PS::S32 ** id_epj,
                                   const PS::S32 *  n_epj,
                                   const PS::S32 ** id_spj,
                                   const PS::S32 *  n_spj,
                                   const EPJ * epj,
                                   const PS::S32 n_epj_tot,
                                   const PS::SPJMonopole * spj,
                                   const PS::S32 n_spj_tot,
                                   const bool send_flag);

PS::S32 DispatchKernelStream(const PS::S32 tag,
                             const PS::S32 n_walk,
                             const EPI ** epi,
                             const PS::S32 *  n_epi,
                             const EPJ ** epj,
                             const PS::S32 *  n_epj,
                             const PS::SPJMonopole ** spj,
                             const PS::S32  * n_spj);

PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 * ni,
                       Force      ** force);

PS::S32 RetrieveKernel2(const PS::S32 tag,
                        const PS::S32 n_walk,
                        const PS::S32 * ni,
                        Force      ** force);

PS::S32 RetrieveKernelStream(const PS::S32 tag,
                             const PS::S32 n_walk,
                             const PS::S32 * ni,
                             Force      ** force);

PS::S32 RetrieveKernelStream2(const PS::S32 tag,
                              const PS::S32 n_walk,
                              const PS::S32 * ni,
                              Force      ** force);
