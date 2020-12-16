#pragma once
/* Standard headers */
#include <string>
#include <vector>
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "user_defined.hpp"

namespace FDPS_Manipulators {
   extern void Initialize(int argc, char *argv[]);
   extern void PS_Initialize();
   extern void PS_Finalize();
   //--------------------------------
   // ParticleSystem manipulators
   //--------------------------------
   extern void create_psys(int *psys_num, char *psys_info);
   extern void delete_psys(int *psys_num);
   extern void init_psys(int *psys_num);
   extern void get_psys_info(int *psys_num, char *psys_info, size_t *charlen);
   extern void set_nptcl_loc(int *psys_num,int *numPtcl);
   extern int get_nptcl_loc(int *psys_num);
   extern void get_psys_cptr(int *psys_num, void **cptr);
   extern void exchange_particle(int *psys_num, int *dinfo_num);
   //----------------------------
   // DomainInfo manipulators
   //----------------------------
   extern void create_dinfo(int *dinfo_num);
   extern void delete_dinfo(int *dinfo_num);
   extern void init_dinfo(int *dinfo_num);
   extern void decompose_domain_all(int *dinfo_num, int *psys_num);
   //----------------------------
   // TreeForForce manipulators
   //----------------------------
   extern void create_tree(int *tree_num, char *tree_info);
   extern void delete_tree(int *tree_num);
   extern void init_tree(int *tree_num, int *numPtcl);
   extern void get_tree_info(int *tree_num, char *tree_info, size_t *charlen);
   extern void calc_force_all_and_write_back_s(int *tree_num,
                                               void (*pfunc_ep_ep)(struct full_particle *, 
                                                                   int ,
                                                                   struct full_particle *,
                                                                   int ,
                                                                   struct full_particle *),
                                               int *psys_num,
                                               int *dinfo_num);
   extern void calc_force_all_and_write_back_l(int *tree_num,
                                               void (*pfunc_ep_ep)(struct full_particle *,
                                                                   int,
                                                                   struct full_particle *,
                                                                   int ,
                                                                   struct full_particle *),
                                               void (*pfunc_ep_sp)(struct full_particle *,
                                                                   int ,
                                                                   PS::SPJMonopole *,
                                                                   int ,
                                                                   struct full_particle *),
                                               int *psys_num,
                                               int *dinfo_num);
   //----------------------------
   // Utility functions
   //----------------------------
   extern void mt_init_genrand(int *s);
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
