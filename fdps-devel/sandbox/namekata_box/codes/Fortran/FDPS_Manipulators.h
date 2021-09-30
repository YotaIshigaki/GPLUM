#pragma once
/* Standard headers */
#include <string>
#include <vector>
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "user_defined.hpp"

namespace FDPS_Manipulators {
   //--------------------------------
   // Parameters
   //--------------------------------
   enum PS_BOUNDARY_CONDITION {
      BC_OPEN,
      BC_PERIODIC_X,
      BC_PERIODIC_Y,
      BC_PERIODIC_Z,
      BC_PERIODIC_XY,
      BC_PERIODIC_XZ,
      BC_PERIODIC_YZ,
      BC_PERIODIC_XYZ,
      BC_SHEARING_BOX,
      BC_USER_DEFINED
   };

   //--------------------------------
   // Initializer of FDPS_Manipulators
   //--------------------------------
   extern void Initialize(int argc, char *argv[]);
   //--------------------------------
   // FDPS Initializer/Finalizer
   //--------------------------------
   extern void PS_Initialize();
   extern void PS_Finalize();
   extern void PS_Abort(const int err_num);
   //--------------------------------
   // ParticleSystem manipulators
   //--------------------------------
   extern void create_psys(int *psys_num,
                           char *psys_info);
   extern void delete_psys(const int psys_num);
   extern void init_psys(const int psys_num);
   extern void get_psys_info(const int psys_num,
                             char *psys_info,
                             size_t *charlen);
   extern long long int get_psys_memsize(const int psys_num);
   extern void get_psys_time_prof(const int psys_num, 
                                  PS::TimeProfile *prof);
   extern void clear_psys_time_prof(const int psys_num);
   extern void set_nptcl_smpl(const int psys_num,
                              const int numPtcl);
   extern void set_nptcl_loc(const int psys_num,
                             const int numPtcl);
   extern int get_nptcl_loc(const int psys_num);
   extern int get_nptcl_glb(const int psys_num);
   extern void get_psys_cptr(const int psys_num,
                             void **cptr);
   extern void exchange_particle(const int psys_num,
                                 const int dinfo_num);
   // fdps-autogen:add_particle; 
   void add_particle00000(const int psys_num,
                          const full_particle *ptcl);
   
   extern void remove_particle(const int psys_num, 
                               const int numPtcl,
                               int *ptcl_indx);

   //----------------------------
   // DomainInfo manipulators
   //----------------------------
   extern void create_dinfo(int *dinfo_num);
   extern void delete_dinfo(const int dinfo_num);
   extern void init_dinfo(const int dinfo_num,
                          const float coef_ema);
   extern void get_dinfo_time_prof(const int dinfo_num, 
                                   PS::TimeProfile *prof);
   extern void clear_dinfo_time_prof(const int dinfo_num);
   extern void set_nums_domain(const int dinfo_num,
                               const int nx,
                               const int ny,
                               const int nz);
   extern void set_boundary_condition(const int dinfo_num, 
                                      const enum PS_BOUNDARY_CONDITION bc);
   extern void set_pos_root_domain(const int dinfo_num,
                                   const PS::F32vec low,
                                   const PS::F32vec high);
   extern void collect_sample_particle(const int dinfo_num,
                                       const int psys_num, 
                                       const bool clear,
                                       const float weight);
   extern void decompose_domain(const int dinfo_num);
   extern void decompose_domain_all(const int dinfo_num,
                                    const int psys_num,
                                    const float weight);

   //----------------------------
   // TreeForForce manipulators
   //----------------------------
   extern void create_tree(int *tree_num,
                           char *tree_info);
   extern void delete_tree(const int tree_num);
   extern void init_tree(const int tree_num,
                         const int numPtcl,
                         const float theta, 
                         const int n_leaf_limit,
                         const int n_group_limit);
   extern void get_tree_info(const int tree_num,
                             char *tree_info,
                             size_t *charlen);
   extern long long int get_tree_memsize(const int tree_num);
   extern void get_tree_time_prof(const int tree_num, 
                                  PS::TimeProfile *prof);
   extern void clear_tree_time_prof(const int tree_num);
   extern long long int get_num_interact_ep_ep_loc(const int tree_num);
   extern long long int get_num_interact_ep_sp_loc(const int tree_num);
   extern long long int get_num_interact_ep_ep_glb(const int tree_num);
   extern long long int get_num_interact_ep_sp_glb(const int tree_num);
   extern void clear_num_interact(const int tree_num);
   extern long long int get_num_tree_walk_loc(const int tree_num);
   extern long long int get_num_tree_walk_glb(const int tree_num);
   // fdps-autogen:calc_force_all_and_write_back;
   extern void calc_force_all_and_write_back_l00000(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00001(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00002(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00003(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00004(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00005(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00006(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00007(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00008(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00009(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00010(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00011(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00012(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00013(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00014(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00015(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00016(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00017(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00018(const int tree_num,
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
                                                    const int dinfo_num);
   
   extern void calc_force_all_and_write_back_l00019(const int tree_num,
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
                                                    const int dinfo_num);
   
   // fdps-autogen:calc_force_all;
   extern void calc_force_all_l00000(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00001(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00002(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00003(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00004(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00005(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00006(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00007(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00008(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00009(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00010(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00011(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00012(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00013(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00014(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00015(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00016(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00017(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00018(const int tree_num,
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
                                     const int dinfo_num);
   
   extern void calc_force_all_l00019(const int tree_num,
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
                                     const int dinfo_num);
   
   // fdps-autogen:calc_force_making_tree;
   extern void calc_force_making_tree_l00000(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00001(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00002(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00003(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00004(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00005(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00006(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00007(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00008(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00009(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00010(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00011(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00012(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00013(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00014(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00015(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00016(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00017(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00018(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00019(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00020(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00021(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00022(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00023(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00024(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00025(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00026(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00027(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00028(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00029(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00030(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00031(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00032(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00033(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00034(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00035(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00036(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00037(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00038(const int tree_num,
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
                                             const int dinfo_num);
   
   extern void calc_force_making_tree_l00039(const int tree_num,
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
                                             const int dinfo_num);
   
   // fdps-autogen:calc_force_and_write_back;
   extern void calc_force_and_write_back_l00000(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00001(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00002(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00003(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00004(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00005(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00006(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00007(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00008(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00009(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00010(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00011(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00012(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00013(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00014(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00015(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00016(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00017(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00018(const int tree_num,
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
                                                const int psys_num);
   
   extern void calc_force_and_write_back_l00019(const int tree_num,
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
                                                const int psys_num);
   
   // fdps-autogen:get_neighbor_list;

   //----------------------------
   // Utility functions
   //----------------------------
   extern void mt_init_genrand(int s);
   extern int mt_genrand_int31(void);
   extern double mt_genrand_real1(void);
   extern double mt_genrand_real2(void);
   extern double mt_genrand_real3(void);
   extern double mt_genrand_res53(void);

   //----------------------------
   // Other utility functions
   //----------------------------
   static std::vector<std::string> split(const std::string& s,
                                         const std::string& delim);
   static void check_split(void);
}
