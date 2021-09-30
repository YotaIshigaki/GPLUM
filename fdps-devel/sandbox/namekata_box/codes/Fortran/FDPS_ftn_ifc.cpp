/* Standard headers */
#include <stdio.h>
#include <iostream>
#include <cmath>
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "FDPS_Manipulators.h"
using namespace FDPS_Manipulators;

extern "C" {

//----------------------
//  Basic
//----------------------
void fdps_initialize() {
   PS_Initialize();
}
void fdps_finalize() {
   PS_Finalize();
}
void fdps_abort(const int err_num) {
   PS_Abort(err_num);
}

//----------------------
//  Particle System 
//----------------------
void fdps_create_psys(int *psys_num,
                      char *psys_info) {
   create_psys(psys_num,psys_info);
}
void fdps_delete_psys(const int psys_num) {
   delete_psys(psys_num);
}
void fdps_init_psys(const int psys_num) {
   init_psys(psys_num);
}
void fdps_get_psys_info(const int psys_num,
                        char *psys_info,
                        size_t *charlen) {
   get_psys_info(psys_num,psys_info,charlen);
}
long long int fdps_get_psys_memsize(const int psys_num) {
   return get_psys_memsize(psys_num);
}
void fdps_get_psys_time_prof(const int psys_num,
                             PS::TimeProfile *prof) {
   get_psys_time_prof(psys_num,prof);
}
void fdps_clear_psys_time_prof(const int psys_num) {
   clear_psys_time_prof(psys_num);
}
void fdps_set_nptcl_smpl(const int psys_num,
                         const int nptcl) {
   set_nptcl_smpl(psys_num,nptcl);
}
void fdps_set_nptcl_loc(const int psys_num,
                        const int nptcl) {
   set_nptcl_loc(psys_num,nptcl);
}
int fdps_get_nptcl_loc(const int psys_num) {
   return get_nptcl_loc(psys_num);
}
int fdps_get_nptcl_glb(const int psys_num) {
   return get_nptcl_glb(psys_num);
}
void fdps_get_psys_cptr(const int psys_num,
                        void **cptr) {
   get_psys_cptr(psys_num,cptr);
}
void fdps_exchange_particle(const int psys_num,
                            const int dinfo_num) {
   exchange_particle(psys_num,dinfo_num);
}
//------------------------------------------------------------------
// [Comment] A place where fdps_add_particle_*() are generated.
// fdps-autogen:add_particle;
void fdps_add_particle00000(const int psys_num,
                            const full_particle *ptcl) {
   add_particle00000(psys_num,ptcl);
}

//------------------------------------------------------------------
void fdps_remove_particle(const int psys_num, 
                          const int nptcl,
                          int *ptcl_indx) {
   remove_particle(psys_num,nptcl,ptcl_indx);
}

//----------------------
//  Domain Info
//----------------------
void fdps_create_dinfo(int *dinfo_num) {
   create_dinfo(dinfo_num);
}
void fdps_delete_dinfo(const int dinfo_num) {
   delete_dinfo(dinfo_num);
}
void fdps_init_dinfo(const int dinfo_num,
                     const float coef_ema) {
   init_dinfo(dinfo_num,coef_ema);
}
void fdps_get_dinfo_time_prof(const int dinfo_num, 
                              PS::TimeProfile *prof) {
   get_dinfo_time_prof(dinfo_num,prof);
}
void fdps_clear_dinfo_time_prof(const int dinfo_num) {
   clear_dinfo_time_prof(dinfo_num);
}
void fdps_set_nums_domain(const int dinfo_num,
                          const int nx,
                          const int ny,
                          const int nz) {
   set_nums_domain(dinfo_num,nx,ny,nz);
}
void fdps_set_boundary_condition(const int dinfo_num,
                                 const enum PS_BOUNDARY_CONDITION bc) {
   set_boundary_condition(dinfo_num,bc);
}
void fdps_set_pos_root_domain(const int dinfo_num,
                              const PS::F32vec low,
                              const PS::F32vec high) {
   set_pos_root_domain(dinfo_num,low,high);
}
void fdps_collect_sample_particle(const int dinfo_num,
                                  const int psys_num, 
                                  const bool clear,
                                  const float weight) {
   collect_sample_particle(dinfo_num,psys_num,clear,weight);
}
void fdps_decompose_domain(const int dinfo_num) {
   decompose_domain(dinfo_num);
}
void fdps_decompose_domain_all(const int dinfo_num,
                               const int psys_num, 
                               const float weight) {
   decompose_domain_all(dinfo_num,psys_num,weight);
}

//----------------------
//  Tree 
//----------------------
void fdps_create_tree(int *tree_num,
                      char *tree_info) {
   create_tree(tree_num,tree_info);
}
void fdps_delete_tree(const int tree_num) {
   delete_tree(tree_num);
}
void fdps_init_tree(const int tree_num,
                    const int nptcl,
                    const float theta,
                    const int n_leaf_limit,
                    const int n_group_limit) {
   init_tree(tree_num,nptcl,theta,
             n_leaf_limit,n_group_limit);
}
void fdps_get_tree_info(const int tree_num,
                        char *tree_info,
                        size_t *charlen) {
   get_tree_info(tree_num,tree_info,charlen);
}
long long int fdps_get_tree_memsize(const int tree_num) {
   return get_tree_memsize(tree_num);
}
void fdps_get_tree_time_prof(const int tree_num,
                             PS::TimeProfile *prof) {
   get_tree_time_prof(tree_num,prof);
}
void fdps_clear_tree_time_prof(const int tree_num) {
   clear_tree_time_prof(tree_num);
}
long long int fdps_get_num_interact_ep_ep_loc(const int tree_num) {
   return get_num_interact_ep_ep_loc(tree_num);
}
long long int fdps_get_num_interact_ep_sp_loc(const int tree_num) {
   return get_num_interact_ep_sp_loc(tree_num);
}
long long int fdps_get_num_interact_ep_ep_glb(const int tree_num) {
   return get_num_interact_ep_ep_glb(tree_num);
}
long long int fdps_get_num_interact_ep_sp_glb(const int tree_num) {
   return get_num_interact_ep_sp_glb(tree_num);
}
void fdps_clear_num_interact(const int tree_num) {
   return clear_num_interact(tree_num);
}
long long int fdps_get_num_tree_walk_loc(const int tree_num) {
   return get_num_tree_walk_loc(tree_num);
}
long long int fdps_get_num_tree_walk_glb(const int tree_num) {
   return get_num_tree_walk_glb(tree_num);
}

//------------------------------------------------------------------
// [Comment] A place where fdps_calc_force_*() are generated.
// fdps-autogen:calc_force_all_and_write_back;
void fdps_calc_force_all_and_write_back_l00000(const int tree_num,
                                               void (*pfunc_ep_ep)(struct full_particle *,
                                                                   int ,
                                                                   struct full_particle *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct full_particle *,
                                                                   int ,
                                                                   PS::SPJMonopole *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00000(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00001(const int tree_num,
                                               void (*pfunc_ep_ep)(struct full_particle *,
                                                                   int ,
                                                                   struct full_particle *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct full_particle *,
                                                                   int ,
                                                                   PS::SPJQuadrupole *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00001(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00002(const int tree_num,
                                               void (*pfunc_ep_ep)(struct full_particle *,
                                                                   int ,
                                                                   struct full_particle *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct full_particle *,
                                                                   int ,
                                                                   PS::SPJMonopoleGeometricCenter *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00002(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00003(const int tree_num,
                                               void (*pfunc_ep_ep)(struct full_particle *,
                                                                   int ,
                                                                   struct full_particle *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct full_particle *,
                                                                   int ,
                                                                   PS::SPJDipoleGeometricCenter *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00003(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00004(const int tree_num,
                                               void (*pfunc_ep_ep)(struct full_particle *,
                                                                   int ,
                                                                   struct full_particle *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct full_particle *,
                                                                   int ,
                                                                   PS::SPJQuadrupoleGeometricCenter *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00004(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00005(const int tree_num,
                                               void (*pfunc_ep_ep)(struct full_particle *,
                                                                   int ,
                                                                   struct essential_particle_j *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct full_particle *,
                                                                   int ,
                                                                   PS::SPJMonopole *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00005(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00006(const int tree_num,
                                               void (*pfunc_ep_ep)(struct full_particle *,
                                                                   int ,
                                                                   struct essential_particle_j *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct full_particle *,
                                                                   int ,
                                                                   PS::SPJQuadrupole *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00006(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00007(const int tree_num,
                                               void (*pfunc_ep_ep)(struct full_particle *,
                                                                   int ,
                                                                   struct essential_particle_j *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct full_particle *,
                                                                   int ,
                                                                   PS::SPJMonopoleGeometricCenter *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00007(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00008(const int tree_num,
                                               void (*pfunc_ep_ep)(struct full_particle *,
                                                                   int ,
                                                                   struct essential_particle_j *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct full_particle *,
                                                                   int ,
                                                                   PS::SPJDipoleGeometricCenter *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00008(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00009(const int tree_num,
                                               void (*pfunc_ep_ep)(struct full_particle *,
                                                                   int ,
                                                                   struct essential_particle_j *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct full_particle *,
                                                                   int ,
                                                                   PS::SPJQuadrupoleGeometricCenter *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00009(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00010(const int tree_num,
                                               void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                   int ,
                                                                   struct full_particle *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                   int ,
                                                                   PS::SPJMonopole *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00010(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00011(const int tree_num,
                                               void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                   int ,
                                                                   struct full_particle *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                   int ,
                                                                   PS::SPJQuadrupole *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00011(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00012(const int tree_num,
                                               void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                   int ,
                                                                   struct full_particle *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                   int ,
                                                                   PS::SPJMonopoleGeometricCenter *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00012(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00013(const int tree_num,
                                               void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                   int ,
                                                                   struct full_particle *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                   int ,
                                                                   PS::SPJDipoleGeometricCenter *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00013(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00014(const int tree_num,
                                               void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                   int ,
                                                                   struct full_particle *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                   int ,
                                                                   PS::SPJQuadrupoleGeometricCenter *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00014(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00015(const int tree_num,
                                               void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                   int ,
                                                                   struct essential_particle_j *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                   int ,
                                                                   PS::SPJMonopole *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00015(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00016(const int tree_num,
                                               void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                   int ,
                                                                   struct essential_particle_j *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                   int ,
                                                                   PS::SPJQuadrupole *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00016(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00017(const int tree_num,
                                               void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                   int ,
                                                                   struct essential_particle_j *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                   int ,
                                                                   PS::SPJMonopoleGeometricCenter *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00017(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00018(const int tree_num,
                                               void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                   int ,
                                                                   struct essential_particle_j *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                   int ,
                                                                   PS::SPJDipoleGeometricCenter *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00018(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

void fdps_calc_force_all_and_write_back_l00019(const int tree_num,
                                               void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                   int ,
                                                                   struct essential_particle_j *,
                                                                   int , 
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                   int ,
                                                                   PS::SPJQuadrupoleGeometricCenter *,
                                                                   int ,
                                                                   struct full_particle *),
                                               const int psys_num,
                                               const int dinfo_num) {
   calc_force_all_and_write_back_l00019(tree_num, 
                                        pfunc_ep_ep,
                                        pfunc_ep_sp,
                                        psys_num,
                                        dinfo_num);
}

// fdps-autogen:calc_force_all;
void fdps_calc_force_all_l00000(const int tree_num,
                                void (*pfunc_ep_ep)(struct full_particle *,
                                                    int ,
                                                    struct full_particle *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct full_particle *,
                                                    int ,
                                                    PS::SPJMonopole *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00000(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00001(const int tree_num,
                                void (*pfunc_ep_ep)(struct full_particle *,
                                                    int ,
                                                    struct full_particle *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct full_particle *,
                                                    int ,
                                                    PS::SPJQuadrupole *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00001(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00002(const int tree_num,
                                void (*pfunc_ep_ep)(struct full_particle *,
                                                    int ,
                                                    struct full_particle *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct full_particle *,
                                                    int ,
                                                    PS::SPJMonopoleGeometricCenter *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00002(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00003(const int tree_num,
                                void (*pfunc_ep_ep)(struct full_particle *,
                                                    int ,
                                                    struct full_particle *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct full_particle *,
                                                    int ,
                                                    PS::SPJDipoleGeometricCenter *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00003(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00004(const int tree_num,
                                void (*pfunc_ep_ep)(struct full_particle *,
                                                    int ,
                                                    struct full_particle *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct full_particle *,
                                                    int ,
                                                    PS::SPJQuadrupoleGeometricCenter *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00004(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00005(const int tree_num,
                                void (*pfunc_ep_ep)(struct full_particle *,
                                                    int ,
                                                    struct essential_particle_j *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct full_particle *,
                                                    int ,
                                                    PS::SPJMonopole *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00005(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00006(const int tree_num,
                                void (*pfunc_ep_ep)(struct full_particle *,
                                                    int ,
                                                    struct essential_particle_j *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct full_particle *,
                                                    int ,
                                                    PS::SPJQuadrupole *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00006(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00007(const int tree_num,
                                void (*pfunc_ep_ep)(struct full_particle *,
                                                    int ,
                                                    struct essential_particle_j *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct full_particle *,
                                                    int ,
                                                    PS::SPJMonopoleGeometricCenter *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00007(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00008(const int tree_num,
                                void (*pfunc_ep_ep)(struct full_particle *,
                                                    int ,
                                                    struct essential_particle_j *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct full_particle *,
                                                    int ,
                                                    PS::SPJDipoleGeometricCenter *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00008(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00009(const int tree_num,
                                void (*pfunc_ep_ep)(struct full_particle *,
                                                    int ,
                                                    struct essential_particle_j *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct full_particle *,
                                                    int ,
                                                    PS::SPJQuadrupoleGeometricCenter *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00009(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00010(const int tree_num,
                                void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                    int ,
                                                    struct full_particle *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                    int ,
                                                    PS::SPJMonopole *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00010(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00011(const int tree_num,
                                void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                    int ,
                                                    struct full_particle *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                    int ,
                                                    PS::SPJQuadrupole *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00011(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00012(const int tree_num,
                                void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                    int ,
                                                    struct full_particle *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                    int ,
                                                    PS::SPJMonopoleGeometricCenter *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00012(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00013(const int tree_num,
                                void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                    int ,
                                                    struct full_particle *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                    int ,
                                                    PS::SPJDipoleGeometricCenter *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00013(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00014(const int tree_num,
                                void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                    int ,
                                                    struct full_particle *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                    int ,
                                                    PS::SPJQuadrupoleGeometricCenter *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00014(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00015(const int tree_num,
                                void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                    int ,
                                                    struct essential_particle_j *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                    int ,
                                                    PS::SPJMonopole *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00015(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00016(const int tree_num,
                                void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                    int ,
                                                    struct essential_particle_j *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                    int ,
                                                    PS::SPJQuadrupole *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00016(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00017(const int tree_num,
                                void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                    int ,
                                                    struct essential_particle_j *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                    int ,
                                                    PS::SPJMonopoleGeometricCenter *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00017(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00018(const int tree_num,
                                void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                    int ,
                                                    struct essential_particle_j *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                    int ,
                                                    PS::SPJDipoleGeometricCenter *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00018(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

void fdps_calc_force_all_l00019(const int tree_num,
                                void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                    int ,
                                                    struct essential_particle_j *,
                                                    int , 
                                                    struct full_particle *),
                                void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                    int ,
                                                    PS::SPJQuadrupoleGeometricCenter *,
                                                    int ,
                                                    struct full_particle *),
                                const int psys_num,
                                const int dinfo_num) {
   calc_force_all_l00019(tree_num, 
                         pfunc_ep_ep,
                         pfunc_ep_sp,
                         psys_num,
                         dinfo_num);
}

// fdps-autogen:calc_force_making_tree;
void fdps_calc_force_making_tree_l00000(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJMonopole *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00000(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00001(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJQuadrupole *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00001(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00002(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJMonopoleGeometricCenter *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00002(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00003(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJDipoleGeometricCenter *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00003(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00004(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJQuadrupoleGeometricCenter *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00004(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00005(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJMonopole *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00005(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00006(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJQuadrupole *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00006(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00007(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJMonopoleGeometricCenter *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00007(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00008(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJDipoleGeometricCenter *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00008(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00009(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJQuadrupoleGeometricCenter *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00009(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00010(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJMonopole *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00010(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00011(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJQuadrupole *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00011(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00012(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJMonopoleGeometricCenter *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00012(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00013(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJDipoleGeometricCenter *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00013(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00014(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJQuadrupoleGeometricCenter *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00014(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00015(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJMonopole *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00015(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00016(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJQuadrupole *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00016(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00017(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJMonopoleGeometricCenter *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00017(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00018(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJDipoleGeometricCenter *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00018(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00019(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct full_particle *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJQuadrupoleGeometricCenter *,
                                                            int ,
                                                            struct full_particle *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00019(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00020(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJMonopole *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00020(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00021(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJQuadrupole *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00021(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00022(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJMonopoleGeometricCenter *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00022(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00023(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJDipoleGeometricCenter *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00023(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00024(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJQuadrupoleGeometricCenter *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00024(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00025(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJMonopole *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00025(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00026(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJQuadrupole *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00026(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00027(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJMonopoleGeometricCenter *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00027(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00028(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJDipoleGeometricCenter *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00028(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00029(int const tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct full_particle *,
                                                            int ,
                                                            PS::SPJQuadrupoleGeometricCenter *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00029(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00030(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJMonopole *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00030(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00031(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJQuadrupole *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00031(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00032(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJMonopoleGeometricCenter *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00032(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00033(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJDipoleGeometricCenter *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00033(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00034(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct full_particle *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJQuadrupoleGeometricCenter *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00034(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00035(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJMonopole *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00035(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00036(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJQuadrupole *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00036(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00037(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJMonopoleGeometricCenter *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00037(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00038(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJDipoleGeometricCenter *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00038(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

void fdps_calc_force_making_tree_l00039(int const tree_num,
                                        void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                            int ,
                                                            struct essential_particle_j *,
                                                            int , 
                                                            struct force *),
                                        void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                            int ,
                                                            PS::SPJQuadrupoleGeometricCenter *,
                                                            int ,
                                                            struct force *),
                                        const int dinfo_num) {
   calc_force_making_tree_l00039(tree_num, 
                                 pfunc_ep_ep,
                                 pfunc_ep_sp,
                                 dinfo_num);
}

// fdps-autogen:calc_force_and_write_back;
void fdps_calc_force_and_write_back_l00000(const int tree_num,
                                           void (*pfunc_ep_ep)(struct full_particle *,
                                                               int ,
                                                               struct full_particle *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct full_particle *,
                                                               int ,
                                                               PS::SPJMonopole *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00000(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00001(const int tree_num,
                                           void (*pfunc_ep_ep)(struct full_particle *,
                                                               int ,
                                                               struct full_particle *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct full_particle *,
                                                               int ,
                                                               PS::SPJQuadrupole *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00001(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00002(const int tree_num,
                                           void (*pfunc_ep_ep)(struct full_particle *,
                                                               int ,
                                                               struct full_particle *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct full_particle *,
                                                               int ,
                                                               PS::SPJMonopoleGeometricCenter *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00002(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00003(const int tree_num,
                                           void (*pfunc_ep_ep)(struct full_particle *,
                                                               int ,
                                                               struct full_particle *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct full_particle *,
                                                               int ,
                                                               PS::SPJDipoleGeometricCenter *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00003(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00004(const int tree_num,
                                           void (*pfunc_ep_ep)(struct full_particle *,
                                                               int ,
                                                               struct full_particle *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct full_particle *,
                                                               int ,
                                                               PS::SPJQuadrupoleGeometricCenter *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00004(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00005(const int tree_num,
                                           void (*pfunc_ep_ep)(struct full_particle *,
                                                               int ,
                                                               struct essential_particle_j *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct full_particle *,
                                                               int ,
                                                               PS::SPJMonopole *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00005(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00006(const int tree_num,
                                           void (*pfunc_ep_ep)(struct full_particle *,
                                                               int ,
                                                               struct essential_particle_j *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct full_particle *,
                                                               int ,
                                                               PS::SPJQuadrupole *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00006(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00007(const int tree_num,
                                           void (*pfunc_ep_ep)(struct full_particle *,
                                                               int ,
                                                               struct essential_particle_j *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct full_particle *,
                                                               int ,
                                                               PS::SPJMonopoleGeometricCenter *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00007(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00008(const int tree_num,
                                           void (*pfunc_ep_ep)(struct full_particle *,
                                                               int ,
                                                               struct essential_particle_j *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct full_particle *,
                                                               int ,
                                                               PS::SPJDipoleGeometricCenter *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00008(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00009(const int tree_num,
                                           void (*pfunc_ep_ep)(struct full_particle *,
                                                               int ,
                                                               struct essential_particle_j *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct full_particle *,
                                                               int ,
                                                               PS::SPJQuadrupoleGeometricCenter *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00009(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00010(const int tree_num,
                                           void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                               int ,
                                                               struct full_particle *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                               int ,
                                                               PS::SPJMonopole *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00010(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00011(const int tree_num,
                                           void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                               int ,
                                                               struct full_particle *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                               int ,
                                                               PS::SPJQuadrupole *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00011(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00012(const int tree_num,
                                           void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                               int ,
                                                               struct full_particle *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                               int ,
                                                               PS::SPJMonopoleGeometricCenter *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00012(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00013(const int tree_num,
                                           void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                               int ,
                                                               struct full_particle *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                               int ,
                                                               PS::SPJDipoleGeometricCenter *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00013(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00014(const int tree_num,
                                           void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                               int ,
                                                               struct full_particle *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                               int ,
                                                               PS::SPJQuadrupoleGeometricCenter *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00014(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00015(const int tree_num,
                                           void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                               int ,
                                                               struct essential_particle_j *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                               int ,
                                                               PS::SPJMonopole *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00015(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00016(const int tree_num,
                                           void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                               int ,
                                                               struct essential_particle_j *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                               int ,
                                                               PS::SPJQuadrupole *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00016(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00017(const int tree_num,
                                           void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                               int ,
                                                               struct essential_particle_j *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                               int ,
                                                               PS::SPJMonopoleGeometricCenter *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00017(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00018(const int tree_num,
                                           void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                               int ,
                                                               struct essential_particle_j *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                               int ,
                                                               PS::SPJDipoleGeometricCenter *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00018(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

void fdps_calc_force_and_write_back_l00019(const int tree_num,
                                           void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                               int ,
                                                               struct essential_particle_j *,
                                                               int , 
                                                               struct full_particle *),
                                           void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                               int ,
                                                               PS::SPJQuadrupoleGeometricCenter *,
                                                               int ,
                                                               struct full_particle *),
                                           const int psys_num) {
   calc_force_and_write_back_l00019(tree_num, 
                                    pfunc_ep_ep,
                                    pfunc_ep_sp,
                                    psys_num);
}

//------------------------------------------------------------------

//------------------------------------------------------------------
// [Comment] A place where fdps_get_neighbor_list_*() are generated.
// fdps-autogen:get_neighbor_list;
//------------------------------------------------------------------

//----------------------
//  MPI comm. 
//----------------------
int fdps_get_rank() {
   return PS::Comm::getRank();
}
int fdps_get_num_procs() {
   return PS::Comm::getNumberOfProc();
}
int fdps_get_rank_multi_dim(const int id) {
   return PS::Comm::getRankMultiDim(id);
}
int fdps_get_num_procs_multi_dim(const int id) {
   return PS::Comm::getNumberOfProcMultiDim(id);
}
void fdps_get_logical_and(const bool in,
                          bool *out) {
   *out = PS::Comm::synchronizeConditionalBranchAND(in);
}
void fdps_get_logical_or(const bool in,
                         bool *out) {
   *out = PS::Comm::synchronizeConditionalBranchOR(in);
}

//------------------------------------------------------------------
// [Comment] A place where fdps_get_min_value*(), etc. are generated.
// fdps-autogen:get_min_value;
void fdps_get_min_value_i32(const int f_in, int *f_out) {
    //*f_out = PS::Comm::getMinValue(f_in);
    int tmp = f_in;
    *f_out = PS::Comm::getMinValue(tmp);
}

void fdps_get_min_value_i64(const PS::S64 f_in, PS::S64 *f_out) {
    //*f_out = PS::Comm::getMinValue(f_in);
    PS::S64 tmp = f_in;
    *f_out = PS::Comm::getMinValue(tmp);
}

void fdps_get_min_value_r32(const PS::F32 f_in, PS::F32 *f_out) {
    //*f_out = PS::Comm::getMinValue(f_in);
    PS::F32 tmp = f_in;
    *f_out = PS::Comm::getMinValue(tmp);
}

void fdps_get_min_value_r64(const PS::F64 f_in, PS::F64 *f_out) {
    //*f_out = PS::Comm::getMinValue(f_in);
    PS::F64 tmp = f_in;
    *f_out = PS::Comm::getMinValue(tmp);
}

void fdps_get_min_value_w_id_r32(const PS::F32 f_in,
                                 const int i_in,
                                 PS::F32 *f_out,
                                 int *i_out) {
    PS::Comm::getMinValue(f_in,i_in,*f_out,*i_out);
}

void fdps_get_min_value_w_id_r64(const PS::F64 f_in,
                                 const int i_in,
                                 PS::F64 *f_out,
                                 int *i_out) {
    PS::Comm::getMinValue(f_in,i_in,*f_out,*i_out);
}

// fdps-autogen:get_max_value;
void fdps_get_max_value_i32(const int f_in, int *f_out) {
    //*f_out = PS::Comm::getMaxValue(f_in);
    int tmp = f_in;
    *f_out = PS::Comm::getMaxValue(tmp);
}

void fdps_get_max_value_i64(const PS::S64 f_in, PS::S64 *f_out) {
    //*f_out = PS::Comm::getMaxValue(f_in);
    PS::S64 tmp = f_in;
    *f_out = PS::Comm::getMaxValue(tmp);
}

void fdps_get_max_value_r32(const PS::F32 f_in, PS::F32 *f_out) {
    //*f_out = PS::Comm::getMaxValue(f_in);
    PS::F32 tmp = f_in;
    *f_out = PS::Comm::getMaxValue(tmp);
}

void fdps_get_max_value_r64(const PS::F64 f_in, PS::F64 *f_out) {
    //*f_out = PS::Comm::getMaxValue(f_in);
    PS::F64 tmp = f_in;
    *f_out = PS::Comm::getMaxValue(tmp);
}

void fdps_get_max_value_w_id_r32(const PS::F32 f_in,
                                 const int i_in,
                                 PS::F32 *f_out,
                                 int *i_out) {
    PS::Comm::getMaxValue(f_in,i_in,*f_out,*i_out);
}

void fdps_get_max_value_w_id_r64(const PS::F64 f_in,
                                 const int i_in,
                                 PS::F64 *f_out,
                                 int *i_out) {
    PS::Comm::getMaxValue(f_in,i_in,*f_out,*i_out);
}

// fdps-autogen:get_sum;
void fdps_get_sum_i32(const int f_in, int *f_out) {
    //*f_out = PS::Comm::getSum(f_in);
    int tmp = f_in;
    *f_out = PS::Comm::getSum(tmp);
}

void fdps_get_sum_i64(const PS::S64 f_in, PS::S64 *f_out) {
    //*f_out = PS::Comm::getSum(f_in);
    PS::S64 tmp = f_in;
    *f_out = PS::Comm::getSum(tmp);
}

void fdps_get_sum_r32(const PS::F32 f_in, PS::F32 *f_out) {
    //*f_out = PS::Comm::getSum(f_in);
    PS::F32 tmp = f_in;
    *f_out = PS::Comm::getSum(tmp);
}

void fdps_get_sum_r64(const PS::F64 f_in, PS::F64 *f_out) {
    //*f_out = PS::Comm::getSum(f_in);
    PS::F64 tmp = f_in;
    *f_out = PS::Comm::getSum(tmp);
}

// fdps-autogen:broadcast;
void fdps_broadcast_scalar_i32(int *val, int n, int src) {
    PS::Comm::broadcast(val,n,src);
}

void fdps_broadcast_array_i32(int *val, int n, int src) {
    PS::Comm::broadcast(val,n,src);
}

void fdps_broadcast_scalar_i64(PS::S64 *val, int n, int src) {
    PS::Comm::broadcast(val,n,src);
}

void fdps_broadcast_array_i64(PS::S64 *val, int n, int src) {
    PS::Comm::broadcast(val,n,src);
}

void fdps_broadcast_scalar_r32(PS::F32 *val, int n, int src) {
    PS::Comm::broadcast(val,n,src);
}

void fdps_broadcast_array_r32(PS::F32 *val, int n, int src) {
    PS::Comm::broadcast(val,n,src);
}

void fdps_broadcast_scalar_r64(PS::F64 *val, int n, int src) {
    PS::Comm::broadcast(val,n,src);
}

void fdps_broadcast_array_r64(PS::F64 *val, int n, int src) {
    PS::Comm::broadcast(val,n,src);
}

//------------------------------------------------------------------

double fdps_get_wtime() {
   return PS::GetWtime();
}

//----------------------
//  Utility
//----------------------
void fdps_mt_init_genrand(const int s) {
   mt_init_genrand(s);
}
int fdps_mt_genrand_int31() {
   return mt_genrand_int31();
}
double fdps_mt_genrand_real1() {
   return mt_genrand_real1();
}
double fdps_mt_genrand_real2() {
   return mt_genrand_real2();
}
double fdps_mt_genrand_real3() {
   return mt_genrand_real3();
}
double fdps_mt_genrand_res53() {
   return mt_genrand_res53();
}

} // END of extern "C"
