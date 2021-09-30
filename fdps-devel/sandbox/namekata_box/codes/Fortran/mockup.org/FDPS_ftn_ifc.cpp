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

void fdps_initialize() {
   PS_Initialize();
}

void fdps_finalize() {
   PS_Finalize();
}

//----------------------
//  Particle System 
//----------------------
void fdps_create_psys(int *psys_num, char *psys_info) {
   create_psys(psys_num,psys_info);
}
void fdps_delete_psys(int *psys_num) {
   delete_psys(psys_num);
}
void fdps_init_psys(int *psys_num) {
   init_psys(psys_num);
}
void fdps_get_psys_info(int *psys_num, char *psys_info, size_t *charlen) {
   get_psys_info(psys_num,psys_info,charlen);
}
void fdps_set_nptcl_loc(int *psys_num, int *nptcl) {
   set_nptcl_loc(psys_num,nptcl);
}
int fdps_get_nptcl_loc(int *psys_num) {
   return get_nptcl_loc(psys_num);
}
void fdps_get_psys_cptr(int *psys_num, void **cptr) {
   get_psys_cptr(psys_num,cptr);
}
void fdps_exchange_particle(int *psys_num, int *dinfo_num) {
   exchange_particle(psys_num,dinfo_num);
}

//----------------------
//  Domain Info
//----------------------
void fdps_create_dinfo(int *dinfo_num) {
   create_dinfo(dinfo_num);
}
void fdps_delete_dinfo(int *dinfo_num) {
   delete_dinfo(dinfo_num);
}
void fdps_init_dinfo(int *dinfo_num) {
   init_dinfo(dinfo_num);
}
void fdps_decompose_domain_all(int *dinfo_num, int *psys_num) {
   decompose_domain_all(dinfo_num,psys_num);
}

//----------------------
//  Tree 
//----------------------
void fdps_create_tree(int *tree_num, char *tree_info) {
   create_tree(tree_num,tree_info);
}
void fdps_delete_tree(int *tree_num) {
   delete_tree(tree_num);
}
void fdps_init_tree(int *tree_num, int *nptcl) {
   init_tree(tree_num,nptcl);
}
void fdps_get_tree_info(int *tree_num, char *tree_info, size_t *charlen) {
   get_tree_info(tree_num,tree_info,charlen);
}

// [Comment] A place where fdps_calc_force_all_and_write_back_*()
//           are generated.
// fdps-auto-gen-key:calcForceAllAndWriteBack():cpp_ifc
void fdps_calc_force_all_and_write_back_s(int *tree_num,
                                          void (*pfunc_ep_ep)(struct full_particle *,
                                                              int ,
                                                              struct full_particle *,
                                                              int , 
                                                              struct full_particle *),
                                           int *psys_num,
                                           int *dinfo_num) {
   calc_force_all_and_write_back_s(tree_num,
                                   pfunc_ep_ep,
                                   psys_num,
                                   dinfo_num);
}
void fdps_calc_force_all_and_write_back_l(int *tree_num,
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
                                          int *psys_num,
                                          int *dinfo_num) {
   calc_force_all_and_write_back_l(tree_num,
                                   pfunc_ep_ep,
                                   pfunc_ep_sp,
                                   psys_num,
                                   dinfo_num);
}

//----------------------
//  MPI comm. 
//----------------------
int fdps_comm_get_rank() {
   return PS::Comm::getRank();
}
int fdps_comm_get_num_procs() {
   return PS::Comm::getNumberOfProc();
}

//----------------------
//  Utility
//----------------------
void fdps_mt_init_genrand(int *s) {
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
