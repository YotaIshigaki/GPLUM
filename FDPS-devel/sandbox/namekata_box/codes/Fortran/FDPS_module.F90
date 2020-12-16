!==================================
!    MODULE: FDPS module
!==================================
module FDPS_module
   use, intrinsic :: iso_c_binding
   implicit none

   !* Private parameters
   !**** data types
   integer, parameter, private :: int16 = selected_int_kind(4)
   integer, parameter, private :: int32 = selected_int_kind(9)
   integer, parameter, private :: int64 = selected_int_kind(18)
   integer, parameter, private :: real32 = selected_real_kind(6)
   integer, parameter, private :: real64 = selected_real_kind(15)
   !**** space dimension
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
   integer, parameter, private :: space_dim=3
#else
   integer, parameter, private :: space_dim=2
#endif

   !* Enum types
   !**** PS::BOUNDARY_CONDITION
   enum, bind(c)
      enumerator :: fdps_bc_open
      enumerator :: fdps_bc_periodic_x
      enumerator :: fdps_bc_periodic_y
      enumerator :: fdps_bc_periodic_z
      enumerator :: fdps_bc_periodic_xy
      enumerator :: fdps_bc_periodic_xz
      enumerator :: fdps_bc_periodic_yz
      enumerator :: fdps_bc_periodic_xyz
      enumerator :: fdps_bc_shearing_box
      enumerator :: fdps_bc_user_defined
   end enum

   !**** FDPS controller
   type, public :: FDPS_controller
      ! この構造体は external 文をユーザに書かせないため必要.

      ! Members
      integer :: x 

   contains
      ! Methods
      !-(basic functions)
      procedure :: PS_initialize
      procedure :: PS_finalize
      procedure :: PS_abort
      !-(ParticleSystem functions)
      procedure :: create_psys
      procedure :: delete_psys
      procedure :: init_psys
      procedure :: get_psys_info
      procedure :: get_psys_memsize
      procedure :: get_psys_time_prof
      procedure :: clear_psys_time_prof
      procedure :: set_nptcl_smpl
      procedure :: set_nptcl_loc
      procedure :: get_nptcl_loc
      procedure :: get_nptcl_glb
      !-----------------------------------------------
      ! [Comment] A place where private procedures 
      !           and generic procedure of get_psys_fptr*() 
      !           are generated.
      ! fdps-autogen:get_psys_fptr:method;
      procedure, private :: get_psys_fptr00000
      generic :: get_psys_fptr => get_psys_fptr00000
      !-----------------------------------------------
      procedure :: exchange_particle
      !-----------------------------------------------
      ! [Comment] A place where private procedures 
      !           and generic procedure of add_particle*() 
      !           are generated.
      ! fdps-autogen:add_particle:method;
      procedure, private :: add_particle00000
      generic :: add_particle => add_particle00000
      !-----------------------------------------------
      procedure :: remove_particle
      !-(DomainInfo functions)
      procedure :: create_dinfo
      procedure :: delete_dinfo
      procedure :: init_dinfo
      procedure :: get_dinfo_time_prof
      procedure :: clear_dinfo_time_prof
      procedure :: set_nums_domain
      procedure :: set_boundary_condition
      procedure, private :: set_pos_root_domain_a32
      procedure, private :: set_pos_root_domain_a64
      procedure, private :: set_pos_root_domain_v32
      procedure, private :: set_pos_root_domain_v64
      generic :: set_pos_root_domain => set_pos_root_domain_a32, &
                                        set_pos_root_domain_a64, &
                                        set_pos_root_domain_v32, &
                                        set_pos_root_domain_v64
      procedure :: collect_sample_particle
      procedure :: decompose_domain
      procedure :: decompose_domain_all
      !-(Tree functions)
      procedure :: create_tree
      procedure :: delete_tree
      procedure :: init_tree
      procedure :: get_tree_info
      procedure :: get_tree_memsize
      procedure :: get_tree_time_prof
      procedure :: clear_tree_time_prof
      procedure :: get_num_interact_ep_ep_loc
      procedure :: get_num_interact_ep_sp_loc
      procedure :: get_num_interact_ep_ep_glb
      procedure :: get_num_interact_ep_sp_glb
      procedure :: clear_num_interact
      procedure :: get_num_tree_walk_loc
      procedure :: get_num_tree_walk_glb
      procedure, private :: calc_force_all_and_write_back_s
      procedure, private :: calc_force_all_and_write_back_l
      generic   :: calc_force_all_and_write_back => calc_force_all_and_write_back_s, &
                                                    calc_force_all_and_write_back_l
      procedure, private :: calc_force_all_s
      procedure, private :: calc_force_all_l
      generic   :: calc_force_all => calc_force_all_s, &
                                     calc_force_all_l
      procedure, private :: calc_force_making_tree_s
      procedure, private :: calc_force_making_tree_l
      generic   :: calc_force_making_tree => calc_force_making_tree_s, &
                                             calc_force_making_tree_l
      procedure, private :: calc_force_and_write_back_s
      procedure, private :: calc_force_and_write_back_l
      generic   :: calc_force_and_write_back => calc_force_and_write_back_s, &
                                                calc_force_and_write_back_l
      !-----------------------------------------------
      ! [Comment] A place where private procedures 
      !           and generic procedure of get_neighbor_list*() 
      !           are generated.
      ! fdps-autogen:get_neighbor_list:method;
      procedure, private :: get_neighbor_list00000
      procedure, private :: get_neighbor_list00001
      generic :: get_neighbor_list => get_neighbor_list00000, &
                                      get_neighbor_list00001
      !-----------------------------------------------
      !-(MPI comm. functions)
      procedure :: get_rank 
      procedure :: get_rank_multi_dim
      procedure :: get_num_procs
      procedure :: get_num_procs_multi_dim
      procedure :: get_logical_and
      procedure :: get_logical_or
      !-----------------------------------------------
      ! [Comment] A place where private procedures 
      !           and generic procedure of reduction
      !           routines such as get_min_value*
      !           are generated.
      ! fdps-autogen:get_min_value:method;
      procedure, private :: get_min_value_i32
      procedure, private :: get_min_value_i64
      procedure, private :: get_min_value_r32
      procedure, private :: get_min_value_r64
      procedure, private :: get_min_value_w_id_r32
      procedure, private :: get_min_value_w_id_r64
      generic :: get_min_value => get_min_value_i32, &
                                  get_min_value_i64, &
                                  get_min_value_r32, &
                                  get_min_value_r64, &
                                  get_min_value_w_id_r32, &
                                  get_min_value_w_id_r64
      ! fdps-autogen:get_max_value:method;
      procedure, private :: get_max_value_i32
      procedure, private :: get_max_value_i64
      procedure, private :: get_max_value_r32
      procedure, private :: get_max_value_r64
      procedure, private :: get_max_value_w_id_r32
      procedure, private :: get_max_value_w_id_r64
      generic :: get_max_value => get_max_value_i32, &
                                  get_max_value_i64, &
                                  get_max_value_r32, &
                                  get_max_value_r64, &
                                  get_max_value_w_id_r32, &
                                  get_max_value_w_id_r64
      ! fdps-autogen:get_sum:method;
      procedure, private :: get_sum_i32
      procedure, private :: get_sum_i64
      procedure, private :: get_sum_r32
      procedure, private :: get_sum_r64
      generic :: get_sum => get_sum_i32, &
                            get_sum_i64, &
                            get_sum_r32, &
                            get_sum_r64
      ! fdps-autogen:broadcast:method;
      procedure, private :: broadcast_scalar_i32
      procedure, private :: broadcast_array_i32
      procedure, private :: broadcast_scalar_i64
      procedure, private :: broadcast_array_i64
      procedure, private :: broadcast_scalar_r32
      procedure, private :: broadcast_array_r32
      procedure, private :: broadcast_scalar_r64
      procedure, private :: broadcast_array_r64
      generic :: broadcast => broadcast_scalar_i32, &
                              broadcast_array_i32, &
                              broadcast_scalar_i64, &
                              broadcast_array_i64, &
                              broadcast_scalar_r32, &
                              broadcast_array_r32, &
                              broadcast_scalar_r64, &
                              broadcast_array_r64
      !-----------------------------------------------
      procedure :: get_wtime

      !-(Utility functions)
      procedure :: MT_init_genrand
      procedure :: MT_genrand_int31
      procedure :: MT_genrand_real1
      procedure :: MT_genrand_real2
      procedure :: MT_genrand_res53
   end type FDPS_controller

   !* Private routines
   private :: PS_initialize
   private :: PS_finalize
   private :: PS_abort
   !-(ParticleSystem functions)
   private :: create_psys
   private :: delete_psys
   private :: init_psys
   private :: get_psys_info
   private :: get_psys_memsize
   private :: get_psys_time_prof
   private :: clear_psys_time_prof
   private :: set_nptcl_smpl
   private :: set_nptcl_loc
   private :: get_nptcl_loc
   private :: get_nptcl_glb
   !-------------------------------------------
   ! [Comment] A place where private procedures 
   !           get_psys_fptr*() are declared.
   ! fdps-autogen:get_psys_fptr:decl;
   private :: get_psys_fptr00000
   !-------------------------------------------
   private :: exchange_particle
   !-------------------------------------------
   ! [Comment] A place where private procedures 
   !           add_particle*() are declared.
   ! fdps-autogen:add_particle:decl;
   private :: add_particle00000
   !-------------------------------------------
   private :: remove_particle
   !-(DomainInfo functions)
   private :: create_dinfo
   private :: delete_dinfo
   private :: init_dinfo
   private :: get_dinfo_time_prof
   private :: clear_dinfo_time_prof
   private :: set_nums_domain
   private :: set_boundary_condition
   private :: set_pos_root_domain_a32
   private :: set_pos_root_domain_a64
   private :: set_pos_root_domain_v32
   private :: set_pos_root_domain_v64
   private :: collect_sample_particle
   private :: decompose_domain
   private :: decompose_domain_all
   !-(Tree functions)
   private :: create_tree
   private :: delete_tree
   private :: init_tree
   private :: get_tree_info
   private :: get_tree_memsize
   private :: get_tree_time_prof
   private :: clear_tree_time_prof
   private :: get_num_interact_ep_ep_loc
   private :: get_num_interact_ep_sp_loc
   private :: get_num_interact_ep_ep_glb
   private :: get_num_interact_ep_sp_glb
   private :: clear_num_interact
   private :: get_num_tree_walk_loc
   private :: get_num_tree_walk_glb
   private :: calc_force_all_and_write_back_s
   private :: calc_force_all_and_write_back_l
   private :: calc_force_all_s
   private :: calc_force_all_l
   private :: calc_force_making_tree_s
   private :: calc_force_making_tree_l
   private :: calc_force_and_write_back_s
   private :: calc_force_and_write_back_l
   !-------------------------------------------
   ! [Comment] A place where private procedures 
   !           get_neighbor_list*() are declared.
   ! fdps-autogen:get_neighbor_list:decl;
   private :: get_neighbor_list00000
   private :: get_neighbor_list00001
   !-------------------------------------------
   !-(MPI comm. functions)
   private :: get_rank 
   private :: get_rank_multi_dim
   private :: get_num_procs
   private :: get_num_procs_multi_dim
   private :: get_logical_and
   private :: get_logical_or
   !-------------------------------------------
   ! [Comment] A place where private procedures 
   !           for reduction operations are declared.
   ! fdps-autogen:get_min_value:decl;
   private :: get_min_value_i32
   private :: get_min_value_i64
   private :: get_min_value_r32
   private :: get_min_value_r64
   private :: get_min_value_w_id_r32
   private :: get_min_value_w_id_r64
   ! fdps-autogen:get_max_value:decl;
   private :: get_max_value_i32
   private :: get_max_value_i64
   private :: get_max_value_r32
   private :: get_max_value_r64
   private :: get_max_value_w_id_r32
   private :: get_max_value_w_id_r64
   ! fdps-autogen:get_sum:decl;
   private :: get_sum_i32
   private :: get_sum_i64
   private :: get_sum_r32
   private :: get_sum_r64
   ! fdps-autogen:broadcast:decl;
   private :: broadcast_scalar_i32
   private :: broadcast_array_i32
   private :: broadcast_scalar_i64
   private :: broadcast_array_i64
   private :: broadcast_scalar_r32
   private :: broadcast_array_r32
   private :: broadcast_scalar_r64
   private :: broadcast_array_r64
   !-------------------------------------------
   private :: get_wtime
   !-(Utility functions)
   private :: MT_init_genrand
   private :: MT_genrand_int31
   private :: MT_genrand_real1
   private :: MT_genrand_real2
   private :: MT_genrand_res53


   !* C++ function interfaces
   interface
      !--------------------------
      !  Initializer/Finalizer
      !--------------------------
      subroutine fdps_initialize() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
      end subroutine fdps_initialize

      subroutine fdps_finalize() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
      end subroutine fdps_finalize

      subroutine fdps_abort(err_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: err_num
      end subroutine fdps_abort

      !----------------------
      !  Particle System 
      !----------------------
      subroutine fdps_create_psys(psys_num,psys_info) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(inout) :: psys_num
         character(kind=c_char), dimension(*), intent(in) :: psys_info
      end subroutine fdps_create_psys

      subroutine fdps_delete_psys(psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num
      end subroutine fdps_delete_psys

      subroutine fdps_init_psys(psys_num)  bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num
      end subroutine fdps_init_psys

      subroutine fdps_get_psys_info(psys_num,psys_info,charlen)  bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num
         character(kind=c_char), dimension(*), intent(inout) :: psys_info
         integer(kind=c_size_t), intent(inout) :: charlen
      end subroutine fdps_get_psys_info

      function fdps_get_psys_memsize(psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_psys_memsize
         integer(kind=c_int), value, intent(in) :: psys_num
      end function fdps_get_psys_memsize

      subroutine fdps_get_psys_time_prof(psys_num,prof) bind(c)
         use, intrinsic :: iso_c_binding
         use fdps_time_profile
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num
         type(fdps_time_prof), intent(inout) :: prof
      end subroutine fdps_get_psys_time_prof

      subroutine fdps_clear_psys_time_prof(psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num
      end subroutine fdps_clear_psys_time_prof

      subroutine fdps_set_nptcl_smpl(psys_num,nptcl) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num,nptcl
      end subroutine fdps_set_nptcl_smpl

      subroutine fdps_set_nptcl_loc(psys_num,nptcl) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num,nptcl
      end subroutine fdps_set_nptcl_loc

      function fdps_get_nptcl_loc(psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_get_nptcl_loc
         integer(kind=c_int), value, intent(in) :: psys_num
      end function fdps_get_nptcl_loc

      function fdps_get_nptcl_glb(psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_get_nptcl_glb
         integer(kind=c_int), value, intent(in) :: psys_num
      end function fdps_get_nptcl_glb

      subroutine fdps_get_psys_cptr(psys_num,cptr) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num
         type(c_ptr), intent(inout) :: cptr
      end subroutine fdps_get_psys_cptr

      subroutine fdps_exchange_particle(psys_num,dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num,dinfo_num
      end subroutine fdps_exchange_particle

      !----------------------------------------------------------
      ! [Comment] A place where the interface statements for 
      !           fdps_add_particle*() are generated.
      ! fdps-autogen:add_particle:ifc;
   subroutine fdps_add_particle00000(psys_num,ptcl) bind(c)
      use, intrinsic :: iso_c_binding
      use user_defined_types
      implicit none
      integer(kind=c_int), value, intent(in) :: psys_num
      type(full_particle), intent(in) :: ptcl
   end subroutine fdps_add_particle00000
   
      !----------------------------------------------------------

      subroutine fdps_remove_particle(psys_num,nptcl,ptcl_indx) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: psys_num,nptcl
         integer(kind=c_int), dimension(*), intent(in) :: ptcl_indx
      end subroutine fdps_remove_particle

      !----------------------
      !  Domain Info
      !----------------------
      subroutine fdps_create_dinfo(dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(inout) :: dinfo_num 
      end subroutine fdps_create_dinfo

      subroutine fdps_delete_dinfo(dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num 
      end subroutine fdps_delete_dinfo

      subroutine fdps_init_dinfo(dinfo_num,coef_ema) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num 
         real(kind=c_float), value, intent(in) :: coef_ema
      end subroutine fdps_init_dinfo

      subroutine fdps_get_dinfo_time_prof(dinfo_num,prof) bind(c)
         use, intrinsic :: iso_c_binding
         use fdps_time_profile
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num
         type(fdps_time_prof), intent(inout) :: prof
      end subroutine fdps_get_dinfo_time_prof

      subroutine fdps_clear_dinfo_time_prof(dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num
      end subroutine fdps_clear_dinfo_time_prof

      subroutine fdps_set_nums_domain(dinfo_num,nx,ny,nz) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num,nx,ny,nz
      end subroutine fdps_set_nums_domain

      subroutine fdps_set_boundary_condition(dinfo_num,bc) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num
         integer(kind=c_int), value, intent(in) :: bc
      end subroutine fdps_set_boundary_condition

      subroutine fdps_set_pos_root_domain(dinfo_num,low,high) bind(c)
         use, intrinsic :: iso_c_binding
         use fdps_vector
         implicit none
         integer(kind=c_int), value, intent(in):: dinfo_num
         type(fdps_f32vec), value, intent(in) :: low,high
      end subroutine fdps_set_pos_root_domain

      subroutine fdps_collect_sample_particle(dinfo_num,psys_num, &
                                              clear,weight) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num,psys_num
         logical(kind=c_bool), value, intent(in) :: clear
         real(kind=c_float), value, intent(in) :: weight
      end subroutine fdps_collect_sample_particle

      subroutine fdps_decompose_domain(dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num
      end subroutine fdps_decompose_domain

      subroutine fdps_decompose_domain_all(dinfo_num,psys_num, &
                                           weight) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: dinfo_num,psys_num
         real(kind=c_float), value, intent(in) :: weight
      end subroutine fdps_decompose_domain_all

      !----------------------
      !  Tree 
      !----------------------
      subroutine fdps_create_tree(tree_num,tree_info) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(inout) :: tree_num
         character(kind=c_char), dimension(*), intent(in) :: tree_info
      end subroutine fdps_create_tree

      subroutine fdps_delete_tree(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num
      end subroutine fdps_delete_tree

      subroutine fdps_init_tree(tree_num,nptcl,theta, &
                                n_leaf_limit,n_group_limit) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,nptcl
         real(kind=c_float), value, intent(in) :: theta
         integer(kind=c_int), value, intent(in) :: n_leaf_limit,n_group_limit
      end subroutine fdps_init_tree

      subroutine fdps_get_tree_info(tree_num,tree_info,charlen) bind(c) 
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num
         character(kind=c_char), dimension(*), intent(inout) :: tree_info
         integer(kind=c_size_t), intent(inout) :: charlen
      end subroutine fdps_get_tree_info

      function fdps_get_tree_memsize(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_tree_memsize
         integer(kind=c_int), value, intent(in) :: tree_num
      end function fdps_get_tree_memsize

      subroutine fdps_get_tree_time_prof(tree_num,prof) bind(c)
         use, intrinsic :: iso_c_binding
         use fdps_time_profile
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num
         type(fdps_time_prof), intent(inout) :: prof
      end subroutine fdps_get_tree_time_prof

      subroutine fdps_clear_tree_time_prof(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num
      end subroutine fdps_clear_tree_time_prof

      function fdps_get_num_interact_ep_ep_loc(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_num_interact_ep_ep_loc
         integer(kind=c_int), value, intent(in) :: tree_num
      end function fdps_get_num_interact_ep_ep_loc

      function fdps_get_num_interact_ep_sp_loc(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_num_interact_ep_sp_loc
         integer(kind=c_int), value, intent(in) :: tree_num
      end function fdps_get_num_interact_ep_sp_loc

      function fdps_get_num_interact_ep_ep_glb(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_num_interact_ep_ep_glb
         integer(kind=c_int), value, intent(in) :: tree_num
      end function fdps_get_num_interact_ep_ep_glb

      function fdps_get_num_interact_ep_sp_glb(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_num_interact_ep_sp_glb
         integer(kind=c_int), value, intent(in) :: tree_num
      end function fdps_get_num_interact_ep_sp_glb

      subroutine fdps_clear_num_interact(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num
      end subroutine fdps_clear_num_interact

      function fdps_get_num_tree_walk_loc(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_num_tree_walk_loc
         integer(kind=c_int), value, intent(in) :: tree_num
      end function fdps_get_num_tree_walk_loc

      function fdps_get_num_tree_walk_glb(tree_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long) :: fdps_get_num_tree_walk_glb
         integer(kind=c_int), value, intent(in) :: tree_num
      end function fdps_get_num_tree_walk_glb

      !----------------------------------------------------------
      ! [Comment] A place where the interface statements for 
      !           fdps_calc_force_*() are generated.
      ! fdps-autogen:calc_force_all_and_write_back:ifc;
      subroutine fdps_calc_force_all_and_write_back_l00000(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00000
      
      subroutine fdps_calc_force_all_and_write_back_l00001(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00001
      
      subroutine fdps_calc_force_all_and_write_back_l00002(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00002
      
      subroutine fdps_calc_force_all_and_write_back_l00003(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00003
      
      subroutine fdps_calc_force_all_and_write_back_l00004(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00004
      
      subroutine fdps_calc_force_all_and_write_back_l00005(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00005
      
      subroutine fdps_calc_force_all_and_write_back_l00006(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00006
      
      subroutine fdps_calc_force_all_and_write_back_l00007(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00007
      
      subroutine fdps_calc_force_all_and_write_back_l00008(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00008
      
      subroutine fdps_calc_force_all_and_write_back_l00009(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00009
      
      subroutine fdps_calc_force_all_and_write_back_l00010(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00010
      
      subroutine fdps_calc_force_all_and_write_back_l00011(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00011
      
      subroutine fdps_calc_force_all_and_write_back_l00012(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00012
      
      subroutine fdps_calc_force_all_and_write_back_l00013(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00013
      
      subroutine fdps_calc_force_all_and_write_back_l00014(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00014
      
      subroutine fdps_calc_force_all_and_write_back_l00015(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00015
      
      subroutine fdps_calc_force_all_and_write_back_l00016(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00016
      
      subroutine fdps_calc_force_all_and_write_back_l00017(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00017
      
      subroutine fdps_calc_force_all_and_write_back_l00018(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00018
      
      subroutine fdps_calc_force_all_and_write_back_l00019(tree_num, &
                                                           pfunc_ep_ep, &
                                                           pfunc_ep_sp, &
                                                           psys_num,    &
                                                           dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l00019
      
      ! fdps-autogen:calc_force_all:ifc;
      subroutine fdps_calc_force_all_l00000(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00000
      
      subroutine fdps_calc_force_all_l00001(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00001
      
      subroutine fdps_calc_force_all_l00002(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00002
      
      subroutine fdps_calc_force_all_l00003(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00003
      
      subroutine fdps_calc_force_all_l00004(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00004
      
      subroutine fdps_calc_force_all_l00005(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00005
      
      subroutine fdps_calc_force_all_l00006(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00006
      
      subroutine fdps_calc_force_all_l00007(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00007
      
      subroutine fdps_calc_force_all_l00008(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00008
      
      subroutine fdps_calc_force_all_l00009(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00009
      
      subroutine fdps_calc_force_all_l00010(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00010
      
      subroutine fdps_calc_force_all_l00011(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00011
      
      subroutine fdps_calc_force_all_l00012(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00012
      
      subroutine fdps_calc_force_all_l00013(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00013
      
      subroutine fdps_calc_force_all_l00014(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00014
      
      subroutine fdps_calc_force_all_l00015(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00015
      
      subroutine fdps_calc_force_all_l00016(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00016
      
      subroutine fdps_calc_force_all_l00017(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00017
      
      subroutine fdps_calc_force_all_l00018(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00018
      
      subroutine fdps_calc_force_all_l00019(tree_num, &
                                            pfunc_ep_ep, &
                                            pfunc_ep_sp, &
                                            psys_num,    &
                                            dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_l00019
      
      ! fdps-autogen:calc_force_making_tree:ifc;
      subroutine fdps_calc_force_making_tree_l00000(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00000
      
      subroutine fdps_calc_force_making_tree_l00001(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00001
      
      subroutine fdps_calc_force_making_tree_l00002(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00002
      
      subroutine fdps_calc_force_making_tree_l00003(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00003
      
      subroutine fdps_calc_force_making_tree_l00004(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00004
      
      subroutine fdps_calc_force_making_tree_l00005(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00005
      
      subroutine fdps_calc_force_making_tree_l00006(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00006
      
      subroutine fdps_calc_force_making_tree_l00007(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00007
      
      subroutine fdps_calc_force_making_tree_l00008(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00008
      
      subroutine fdps_calc_force_making_tree_l00009(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00009
      
      subroutine fdps_calc_force_making_tree_l00010(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00010
      
      subroutine fdps_calc_force_making_tree_l00011(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00011
      
      subroutine fdps_calc_force_making_tree_l00012(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00012
      
      subroutine fdps_calc_force_making_tree_l00013(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00013
      
      subroutine fdps_calc_force_making_tree_l00014(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00014
      
      subroutine fdps_calc_force_making_tree_l00015(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00015
      
      subroutine fdps_calc_force_making_tree_l00016(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00016
      
      subroutine fdps_calc_force_making_tree_l00017(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00017
      
      subroutine fdps_calc_force_making_tree_l00018(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00018
      
      subroutine fdps_calc_force_making_tree_l00019(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00019
      
      subroutine fdps_calc_force_making_tree_l00020(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00020
      
      subroutine fdps_calc_force_making_tree_l00021(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00021
      
      subroutine fdps_calc_force_making_tree_l00022(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00022
      
      subroutine fdps_calc_force_making_tree_l00023(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00023
      
      subroutine fdps_calc_force_making_tree_l00024(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00024
      
      subroutine fdps_calc_force_making_tree_l00025(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00025
      
      subroutine fdps_calc_force_making_tree_l00026(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00026
      
      subroutine fdps_calc_force_making_tree_l00027(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00027
      
      subroutine fdps_calc_force_making_tree_l00028(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00028
      
      subroutine fdps_calc_force_making_tree_l00029(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00029
      
      subroutine fdps_calc_force_making_tree_l00030(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00030
      
      subroutine fdps_calc_force_making_tree_l00031(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00031
      
      subroutine fdps_calc_force_making_tree_l00032(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00032
      
      subroutine fdps_calc_force_making_tree_l00033(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00033
      
      subroutine fdps_calc_force_making_tree_l00034(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00034
      
      subroutine fdps_calc_force_making_tree_l00035(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00035
      
      subroutine fdps_calc_force_making_tree_l00036(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00036
      
      subroutine fdps_calc_force_making_tree_l00037(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00037
      
      subroutine fdps_calc_force_making_tree_l00038(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00038
      
      subroutine fdps_calc_force_making_tree_l00039(tree_num, &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_making_tree_l00039
      
      ! fdps-autogen:calc_force_and_write_back:ifc;
      subroutine fdps_calc_force_and_write_back_l00000(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00000
      
      subroutine fdps_calc_force_and_write_back_l00001(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00001
      
      subroutine fdps_calc_force_and_write_back_l00002(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00002
      
      subroutine fdps_calc_force_and_write_back_l00003(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00003
      
      subroutine fdps_calc_force_and_write_back_l00004(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00004
      
      subroutine fdps_calc_force_and_write_back_l00005(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00005
      
      subroutine fdps_calc_force_and_write_back_l00006(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00006
      
      subroutine fdps_calc_force_and_write_back_l00007(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00007
      
      subroutine fdps_calc_force_and_write_back_l00008(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00008
      
      subroutine fdps_calc_force_and_write_back_l00009(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00009
      
      subroutine fdps_calc_force_and_write_back_l00010(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00010
      
      subroutine fdps_calc_force_and_write_back_l00011(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00011
      
      subroutine fdps_calc_force_and_write_back_l00012(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00012
      
      subroutine fdps_calc_force_and_write_back_l00013(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00013
      
      subroutine fdps_calc_force_and_write_back_l00014(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00014
      
      subroutine fdps_calc_force_and_write_back_l00015(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00015
      
      subroutine fdps_calc_force_and_write_back_l00016(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00016
      
      subroutine fdps_calc_force_and_write_back_l00017(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00017
      
      subroutine fdps_calc_force_and_write_back_l00018(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00018
      
      subroutine fdps_calc_force_and_write_back_l00019(tree_num, &
                                                       pfunc_ep_ep, &
                                                       pfunc_ep_sp, &
                                                       psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: tree_num,psys_num
         type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_and_write_back_l00019
      
      !----------------------------------------------------------

      !----------------------------------------------------------
      ! [Comment] A place where the interface statements for
      !           fdps_get_neighbor_list*() are generated.
      ! fdps-autogen:get_neighbor_list:ifc;
      !----------------------------------------------------------

      !----------------------
      !  MPI comm. 
      !----------------------
      function fdps_get_rank() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_get_rank
      end function fdps_get_rank

      function fdps_get_rank_multi_dim(id) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_get_rank_multi_dim
         integer(kind=c_int), value, intent(in) :: id
      end function fdps_get_rank_multi_dim

      function fdps_get_num_procs() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_get_num_procs
      end function fdps_get_num_procs

      function fdps_get_num_procs_multi_dim(id) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_get_num_procs_multi_dim
         integer(kind=c_int), value, intent(in) :: id
      end function fdps_get_num_procs_multi_dim

      subroutine fdps_get_logical_and(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         logical(kind=c_bool), value, intent(in) :: f_in
         logical(kind=c_bool), intent(inout) :: f_out
      end subroutine fdps_get_logical_and

      subroutine fdps_get_logical_or(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         logical(kind=c_bool), value, intent(in) :: f_in
         logical(kind=c_bool), intent(inout) :: f_out
      end subroutine fdps_get_logical_or

      !-------------------------------------------
      ! [Comment] A place where the interface statements for 
      !           reduction functions are generated.
      ! fdps-autogen:get_min_value:ftn_ifc;
      subroutine fdps_get_min_value_i32(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: f_in
         integer(kind=c_int), intent(inout) :: f_out
      end subroutine fdps_get_min_value_i32
      
      subroutine fdps_get_min_value_i64(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long), value, intent(in) :: f_in
         integer(kind=c_long_long), intent(inout) :: f_out
      end subroutine fdps_get_min_value_i64
      
      subroutine fdps_get_min_value_r32(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_float), value, intent(in) :: f_in
         real(kind=c_float), intent(inout) :: f_out
      end subroutine fdps_get_min_value_r32
      
      subroutine fdps_get_min_value_r64(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double), value, intent(in) :: f_in
         real(kind=c_double), intent(inout) :: f_out
      end subroutine fdps_get_min_value_r64
      
      subroutine fdps_get_min_value_w_id_r32(f_in,i_in,f_out,i_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_float), value, intent(in) :: f_in
         real(kind=c_float), intent(inout) :: f_out
         integer(kind=c_int), value, intent(in) :: i_in
         integer(kind=c_int), intent(inout) :: i_out
      end subroutine fdps_get_min_value_w_id_r32
      
      subroutine fdps_get_min_value_w_id_r64(f_in,i_in,f_out,i_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double), value, intent(in) :: f_in
         real(kind=c_double), intent(inout) :: f_out
         integer(kind=c_int), value, intent(in) :: i_in
         integer(kind=c_int), intent(inout) :: i_out
      end subroutine fdps_get_min_value_w_id_r64
      
      ! fdps-autogen:get_max_value:ftn_ifc;
      subroutine fdps_get_max_value_i32(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: f_in
         integer(kind=c_int), intent(inout) :: f_out
      end subroutine fdps_get_max_value_i32
      
      subroutine fdps_get_max_value_i64(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long), value, intent(in) :: f_in
         integer(kind=c_long_long), intent(inout) :: f_out
      end subroutine fdps_get_max_value_i64
      
      subroutine fdps_get_max_value_r32(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_float), value, intent(in) :: f_in
         real(kind=c_float), intent(inout) :: f_out
      end subroutine fdps_get_max_value_r32
      
      subroutine fdps_get_max_value_r64(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double), value, intent(in) :: f_in
         real(kind=c_double), intent(inout) :: f_out
      end subroutine fdps_get_max_value_r64
      
      subroutine fdps_get_max_value_w_id_r32(f_in,i_in,f_out,i_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_float), value, intent(in) :: f_in
         real(kind=c_float), intent(inout) :: f_out
         integer(kind=c_int), value, intent(in) :: i_in
         integer(kind=c_int), intent(inout) :: i_out
      end subroutine fdps_get_max_value_w_id_r32
      
      subroutine fdps_get_max_value_w_id_r64(f_in,i_in,f_out,i_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double), value, intent(in) :: f_in
         real(kind=c_double), intent(inout) :: f_out
         integer(kind=c_int), value, intent(in) :: i_in
         integer(kind=c_int), intent(inout) :: i_out
      end subroutine fdps_get_max_value_w_id_r64
      
      ! fdps-autogen:get_sum:ftn_ifc;
      subroutine fdps_get_sum_i32(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: f_in
         integer(kind=c_int), intent(inout) :: f_out
      end subroutine fdps_get_sum_i32
      
      subroutine fdps_get_sum_i64(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long), value, intent(in) :: f_in
         integer(kind=c_long_long), intent(inout) :: f_out
      end subroutine fdps_get_sum_i64
      
      subroutine fdps_get_sum_r32(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_float), value, intent(in) :: f_in
         real(kind=c_float), intent(inout) :: f_out
      end subroutine fdps_get_sum_r32
      
      subroutine fdps_get_sum_r64(f_in,f_out) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double), value, intent(in) :: f_in
         real(kind=c_double), intent(inout) :: f_out
      end subroutine fdps_get_sum_r64
      
      ! fdps-autogen:broadcast:ftn_ifc;
      subroutine fdps_broadcast_scalar_i32(val,n,src) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(in) :: val
         integer(kind=c_int), value, intent(in) :: n,src
      end subroutine fdps_broadcast_scalar_i32
      
      subroutine fdps_broadcast_array_i32(vals,n,src) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: n
         integer(kind=c_int), dimension(n), intent(inout) :: vals
         integer(kind=c_int), value, intent(in) :: src
      end subroutine fdps_broadcast_array_i32
      
      subroutine fdps_broadcast_scalar_i64(val,n,src) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_long_long), intent(in) :: val
         integer(kind=c_int), value, intent(in) :: n,src
      end subroutine fdps_broadcast_scalar_i64
      
      subroutine fdps_broadcast_array_i64(vals,n,src) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: n
         integer(kind=c_long_long), dimension(n), intent(inout) :: vals
         integer(kind=c_int), value, intent(in) :: src
      end subroutine fdps_broadcast_array_i64
      
      subroutine fdps_broadcast_scalar_r32(val,n,src) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_float), intent(in) :: val
         integer(kind=c_int), value, intent(in) :: n,src
      end subroutine fdps_broadcast_scalar_r32
      
      subroutine fdps_broadcast_array_r32(vals,n,src) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: n
         real(kind=c_float), dimension(n), intent(inout) :: vals
         integer(kind=c_int), value, intent(in) :: src
      end subroutine fdps_broadcast_array_r32
      
      subroutine fdps_broadcast_scalar_r64(val,n,src) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double), intent(in) :: val
         integer(kind=c_int), value, intent(in) :: n,src
      end subroutine fdps_broadcast_scalar_r64
      
      subroutine fdps_broadcast_array_r64(vals,n,src) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: n
         real(kind=c_double), dimension(n), intent(inout) :: vals
         integer(kind=c_int), value, intent(in) :: src
      end subroutine fdps_broadcast_array_r64
      
      !-------------------------------------------

      function fdps_get_wtime() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double) :: fdps_get_wtime
      end function fdps_get_wtime

      !----------------------
      !  Utility
      !----------------------
      subroutine fdps_mt_init_genrand(s) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: s 
      end subroutine fdps_mt_init_genrand

      function fdps_mt_genrand_int31() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_mt_genrand_int31
      end function fdps_mt_genrand_int31

      function fdps_mt_genrand_real1() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double) :: fdps_mt_genrand_real1
      end function fdps_mt_genrand_real1

      function fdps_mt_genrand_real2() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double) :: fdps_mt_genrand_real2
      end function fdps_mt_genrand_real2

      function fdps_mt_genrand_real3() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double) :: fdps_mt_genrand_real3
      end function fdps_mt_genrand_real3

      function fdps_mt_genrand_res53() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         real(kind=c_double) :: fdps_mt_genrand_res53
      end function fdps_mt_genrand_res53

   end interface

   contains

   !----------------------------------------------------------
   subroutine PS_initialize(this)
      implicit none
      class(FDPS_controller) :: this

      call fdps_initialize()

   end subroutine PS_initialize 

   !----------------------------------------------------------
   subroutine PS_finalize(this)
      implicit none
      class(FDPS_controller) :: this

      call fdps_finalize()

   end subroutine PS_finalize 

   !----------------------------------------------------------
   subroutine PS_abort(this,err_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN), optional :: err_num

      if (present(err_num)) then
         call fdps_abort(err_num)
      else
         call fdps_abort(-1)
      end if

   end subroutine PS_abort

   !##########################################################
   subroutine create_psys(this,psys_num,psys_info_in)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(INOUT) :: psys_num
      character(len=*,kind=c_char), intent(IN) :: psys_info_in
      character(len=:,kind=c_char), allocatable :: psys_info

      psys_info = trim(psys_info_in) // c_null_char
      call fdps_create_psys(psys_num,psys_info)

   end subroutine create_psys

   !----------------------------------------------------------
   subroutine delete_psys(this,psys_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num

      call fdps_delete_psys(psys_num)

   end subroutine delete_psys

   !----------------------------------------------------------
   subroutine init_psys(this,psys_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num

      call fdps_init_psys(psys_num)

   end subroutine init_psys

   !----------------------------------------------------------
   subroutine get_psys_info(this,psys_num,psys_info)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num
      character(len=*,kind=c_char), intent(INOUT) :: psys_info
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize,kind=c_char) :: buf
      integer(kind=c_size_t) :: charlen

      call fdps_get_psys_info(psys_num,buf,charlen)
      psys_info = buf(1:charlen) ! copy
      ! [** Important **]
       !    You should use the intrinsic function trim() when
       !    you compare psys_info with an immediate value of
       !    a string of characters.
       !    cf.) 
       !    write(*,100)psys_info
       !    write(*,100)trim(psys_info)
       !    100 format(a,'|')

   end subroutine get_psys_info

   !----------------------------------------------------------
   function get_psys_memsize(this,psys_num)
      implicit none
      integer(kind=c_long_long) :: get_psys_memsize
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num

      get_psys_memsize =  fdps_get_psys_memsize(psys_num)
      
   end function get_psys_memsize

   !----------------------------------------------------------
   subroutine get_psys_time_prof(this,psys_num,prof)
      use fdps_time_profile
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num
      type(fdps_time_prof), intent(INOUT) :: prof

      call fdps_get_psys_time_prof(psys_num,prof)
      
   end subroutine get_psys_time_prof

   !----------------------------------------------------------
   subroutine clear_psys_time_prof(this,psys_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num

      call fdps_clear_psys_time_prof(psys_num)
      
   end subroutine clear_psys_time_prof

   !----------------------------------------------------------
   subroutine set_nptcl_smpl(this,psys_num,nptcl)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num,nptcl

      call fdps_set_nptcl_smpl(psys_num,nptcl)

   end subroutine set_nptcl_smpl

   !----------------------------------------------------------
   subroutine set_nptcl_loc(this,psys_num,nptcl)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num,nptcl

      call fdps_set_nptcl_loc(psys_num,nptcl)

   end subroutine set_nptcl_loc

   !----------------------------------------------------------
   function get_nptcl_loc(this,psys_num)
      implicit none
      integer :: get_nptcl_loc
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num

      get_nptcl_loc =  fdps_get_nptcl_loc(psys_num)
      
   end function get_nptcl_loc

   !----------------------------------------------------------
   function get_nptcl_glb(this,psys_num)
      implicit none
      integer :: get_nptcl_glb
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num

      get_nptcl_glb =  fdps_get_nptcl_glb(psys_num)
      
   end function get_nptcl_glb

   !----------------------------------------------------------
   ! [Comment] A place where the definitions or implementations
   !           get_psys_fptr*() are generated.
   ! fdps-autogen:get_psys_fptr:impl;
   subroutine get_psys_fptr00000(this,psys_num,fptr_to_FP)
      use user_defined_types
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num
      type(full_particle), dimension(:), pointer, intent(INOUT) :: fptr_to_FP
      !* Local variables
      integer(kind=c_int) :: nptcl_loc
      type(c_ptr) :: cptr_to_FP
   
      call fdps_get_psys_cptr(psys_num,cptr_to_FP)
      nptcl_loc = fdps_get_nptcl_loc(psys_num)
      call c_f_pointer(cptr_to_FP,fptr_to_FP,[nptcl_loc])
   
   end subroutine get_psys_fptr00000
   

   !----------------------------------------------------------
   subroutine exchange_particle(this,psys_num,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num,dinfo_num

      call fdps_exchange_particle(psys_num,dinfo_num)

   end subroutine exchange_particle

   !----------------------------------------------------------
   ! [Comment] A place where the definitions or implementations
   !           add_particle*() are generated.
   ! fdps-autogen:add_particle:impl;
   subroutine add_particle00000(this,psys_num,ptcl)
      use, intrinsic :: iso_c_binding
      use user_defined_types
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: psys_num
      type(full_particle), intent(in) :: ptcl
   
      call fdps_add_particle00000(psys_num,ptcl)
   
   end subroutine add_particle00000
   

   !----------------------------------------------------------
   subroutine remove_particle(this,psys_num,nptcl,ptcl_indx)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: psys_num
      integer(kind=c_int), intent(IN) :: nptcl
      integer(kind=c_int), dimension(nptcl), intent(IN) :: ptcl_indx

      call fdps_remove_particle(psys_num,nptcl,ptcl_indx)

   end subroutine remove_particle

   !##########################################################
   subroutine create_dinfo(this,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(INOUT) :: dinfo_num

      call fdps_create_dinfo(dinfo_num)

   end subroutine create_dinfo

   !----------------------------------------------------------
   subroutine delete_dinfo(this,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num

      call fdps_delete_dinfo(dinfo_num)

   end subroutine delete_dinfo

   !----------------------------------------------------------
   subroutine init_dinfo(this,dinfo_num,coef_ema)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num
      real(kind=c_float), intent(IN), optional :: coef_ema

      if (present(coef_ema)) then
         call fdps_init_dinfo(dinfo_num,coef_ema)
      else
         call fdps_init_dinfo(dinfo_num,1.0_c_float)
      end if

   end subroutine init_dinfo

   !----------------------------------------------------------
   subroutine get_dinfo_time_prof(this,dinfo_num,prof)
      use fdps_time_profile
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num
      type(fdps_time_prof), intent(INOUT) :: prof

      call fdps_get_dinfo_time_prof(dinfo_num,prof)
      
   end subroutine get_dinfo_time_prof

   !----------------------------------------------------------
   subroutine clear_dinfo_time_prof(this,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num

      call fdps_clear_dinfo_time_prof(dinfo_num)
      
   end subroutine clear_dinfo_time_prof

   !----------------------------------------------------------
   subroutine set_nums_domain(this,dinfo_num,nx,ny,nz)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num,nx,ny
      integer(kind=c_int), intent(IN), optional :: nz

      if (present(nz)) then
         call fdps_set_nums_domain(dinfo_num,nx,ny,nz)
      else
         call fdps_set_nums_domain(dinfo_num,nx,ny,1_c_int)
      end if

   end subroutine set_nums_domain

   !----------------------------------------------------------
   subroutine set_boundary_condition(this,dinfo_num,bc)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num,bc
      
      call fdps_set_boundary_condition(dinfo_num,bc)

   end subroutine set_boundary_condition

   !----------------------------------------------------------
   subroutine set_pos_root_domain_a32(this,dinfo_num,low,high)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num
      real(kind=c_float), dimension(space_dim), intent(in) :: low,high
      !* Local variables
      type(fdps_f32vec) :: low_,high_
      
      low_ = low; high_ = high
      call fdps_set_pos_root_domain(dinfo_num,low_,high_)

   end subroutine set_pos_root_domain_a32

   !----------------------------------------------------------
   subroutine set_pos_root_domain_a64(this,dinfo_num,low,high)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num
      real(kind=c_double), dimension(space_dim), intent(in) :: low,high
      !* Local variables
      type(fdps_f32vec) :: low_,high_
      
      low_ = low; high_ = high
      call fdps_set_pos_root_domain(dinfo_num,low_,high_)

   end subroutine set_pos_root_domain_a64

   !----------------------------------------------------------
   subroutine set_pos_root_domain_v32(this,dinfo_num,low,high)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num
      type(fdps_f32vec), intent(in) :: low,high 

      call fdps_set_pos_root_domain(dinfo_num,low,high)

   end subroutine set_pos_root_domain_v32

   !----------------------------------------------------------
   subroutine set_pos_root_domain_v64(this,dinfo_num,low,high)
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num
      type(fdps_f64vec), intent(in) :: low,high 
      !* Local variables
      type(fdps_f32vec) :: low_,high_
      
      low_ = low; high_ = high
      call fdps_set_pos_root_domain(dinfo_num,low_,high_)

   end subroutine set_pos_root_domain_v64

   !----------------------------------------------------------
   subroutine collect_sample_particle(this,dinfo_num,psys_num, &
                                      clear,weight)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num,psys_num
      logical(kind=c_bool), intent(IN), optional :: clear
      real(kind=c_float), intent(IN), optional :: weight
      !* Local variables
      logical(kind=c_bool) :: clear_
      real(kind=c_float) :: weight_

      if (present(clear)) then
         clear_ = clear
      else
         clear_ = .true. ! default value
      end if
      if (present(weight)) then
         weight_ = weight
      else
         weight_ = 1.0_c_float ! default value
      end if
      call fdps_collect_sample_particle(dinfo_num,psys_num, &
                                        clear_,weight_)

   end subroutine collect_sample_particle

   !----------------------------------------------------------
   subroutine decompose_domain(this,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num

      call fdps_decompose_domain(dinfo_num)

   end subroutine decompose_domain

   !----------------------------------------------------------
   subroutine decompose_domain_all(this,dinfo_num,psys_num, &
                                   weight)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: dinfo_num,psys_num
      real(kind=c_float), intent(IN), optional :: weight
      !* Local variables
      real(kind=c_float) :: wgh

      if (present(weight)) then
         call fdps_decompose_domain_all(dinfo_num,psys_num,weight)
      else
         wgh = real(fdps_get_nptcl_loc(psys_num),kind=c_float)
         call fdps_decompose_domain_all(dinfo_num,psys_num,wgh)
      end if

   end subroutine decompose_domain_all

   !##########################################################
   subroutine create_tree(this,tree_num,tree_info_in)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(INOUT) :: tree_num
      character(len=*,kind=c_char), intent(IN) :: tree_info_in
      character(len=:,kind=c_char), allocatable :: tree_info

      tree_info = trim(tree_info_in) // c_null_char
      call fdps_create_tree(tree_num,tree_info)

   end subroutine create_tree

   !----------------------------------------------------------
   subroutine delete_tree(this,tree_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      call fdps_delete_tree(tree_num)

   end subroutine delete_tree

   !----------------------------------------------------------
   subroutine init_tree(this,tree_num,nptcl, &
                        theta,n_leaf_limit,n_group_limit)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,nptcl
      real(kind=c_float), intent(IN), optional :: theta
      integer(kind=c_int), intent(IN), optional :: n_leaf_limit,n_group_limit
      !* Local variables
      real(kind=c_float) :: theta_
      integer(kind=c_int) :: n_leaf_limit_,n_group_limit_

      if (present(theta)) then
         theta_ = theta
      else
         theta_ = 0.7_c_float
      end if
      if (present(n_leaf_limit)) then
         n_leaf_limit_ = n_leaf_limit
      else
         n_leaf_limit_ = 8_c_int
      end if
      if (present(n_group_limit)) then
         n_group_limit_ = n_group_limit
      else
         n_group_limit_ = 64_c_int
      end if

      call fdps_init_tree(tree_num,      &
                          nptcl,         &
                          theta_,        &
                          n_leaf_limit_, &
                          n_group_limit_)

   end subroutine init_tree

   !----------------------------------------------------------
   subroutine get_tree_info(this,tree_num,tree_info)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num
      character(len=*,kind=c_char), intent(INOUT) :: tree_info
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize,kind=c_char) :: buf
      integer(kind=c_size_t) :: charlen

      call fdps_get_tree_info(tree_num,buf,charlen)
      tree_info = buf(1:charlen) ! copy
      ! [** Important **]
      !    You should use the intrinsic function trim() when
      !    you compare tree_info with an immediate value of
      !    a string of characters.
      !    cf.) 
      !    write(*,100)tree_info
      !    write(*,100)trim(tree_info)
      !    100 format(a,'|')

   end subroutine get_tree_info

   !----------------------------------------------------------
   function get_tree_memsize(this,tree_num)
      implicit none
      integer(kind=c_long_long) :: get_tree_memsize
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      get_tree_memsize =  fdps_get_tree_memsize(tree_num)
      
   end function get_tree_memsize

   !----------------------------------------------------------
   subroutine get_tree_time_prof(this,tree_num,prof)
      use fdps_time_profile
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num
      type(fdps_time_prof), intent(INOUT) :: prof

      call fdps_get_tree_time_prof(tree_num,prof)
      
   end subroutine get_tree_time_prof

   !----------------------------------------------------------
   subroutine clear_tree_time_prof(this,tree_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      call fdps_clear_tree_time_prof(tree_num)
      
   end subroutine clear_tree_time_prof

   !----------------------------------------------------------
   function get_num_interact_ep_ep_loc(this,tree_num)
      implicit none
      integer(kind=c_long_long) :: get_num_interact_ep_ep_loc
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      get_num_interact_ep_ep_loc = fdps_get_num_interact_ep_ep_loc(tree_num)
      
   end function get_num_interact_ep_ep_loc

   !----------------------------------------------------------
   function get_num_interact_ep_sp_loc(this,tree_num)
      implicit none
      integer(kind=c_long_long) :: get_num_interact_ep_sp_loc
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      get_num_interact_ep_sp_loc = fdps_get_num_interact_ep_sp_loc(tree_num)
      
   end function get_num_interact_ep_sp_loc

   !----------------------------------------------------------
   function get_num_interact_ep_ep_glb(this,tree_num)
      implicit none
      integer(kind=c_long_long) :: get_num_interact_ep_ep_glb
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      get_num_interact_ep_ep_glb = fdps_get_num_interact_ep_ep_glb(tree_num)
      
   end function get_num_interact_ep_ep_glb

   !----------------------------------------------------------
   function get_num_interact_ep_sp_glb(this,tree_num)
      implicit none
      integer(kind=c_long_long) :: get_num_interact_ep_sp_glb
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      get_num_interact_ep_sp_glb = fdps_get_num_interact_ep_sp_glb(tree_num)
      
   end function get_num_interact_ep_sp_glb

   !----------------------------------------------------------
   subroutine clear_num_interact(this,tree_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      call fdps_clear_num_interact(tree_num)

   end subroutine clear_num_interact

   !----------------------------------------------------------
   function get_num_tree_walk_loc(this,tree_num)
      implicit none
      integer(kind=c_long_long) :: get_num_tree_walk_loc
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      get_num_tree_walk_loc = fdps_get_num_tree_walk_loc(tree_num)
      
   end function get_num_tree_walk_loc

   !----------------------------------------------------------
   function get_num_tree_walk_glb(this,tree_num)
      implicit none
      integer(kind=c_long_long) :: get_num_tree_walk_glb
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num

      get_num_tree_walk_glb = fdps_get_num_tree_walk_glb(tree_num)
      
   end function get_num_tree_walk_glb

   !----------------------------------------------------------
   subroutine calc_force_all_and_write_back_s(this,        &
                                              tree_num,    &
                                              pfunc_ep_ep, &
                                              psys_num,    &
                                              dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,psys_num,dinfo_num
      type(c_funptr), intent(in) :: pfunc_ep_ep
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize,kind=c_char) :: psys_info,tree_info,info

      call get_psys_info(this,psys_num,psys_info)
      call get_tree_info(this,tree_num,tree_info)
      info = trim(psys_info) // ',' // trim(tree_info)

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_all_and_write_back:impl:short;
      !--------------
      case default
         write(*,'(a)')"something wrong occurs!!"
         call PS_finalize(this)
         stop 0
      end select

   end subroutine calc_force_all_and_write_back_s

   !----------------------------------------------------------
   subroutine calc_force_all_and_write_back_l(this,        &
                                              tree_num,    &
                                              pfunc_ep_ep, &
                                              pfunc_ep_sp, &
                                              psys_num,    &
                                              dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,psys_num,dinfo_num
      type(c_funptr), intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize,kind=c_char) :: psys_info,tree_info,info

      call get_psys_info(this,psys_num,psys_info)
      call get_tree_info(this,tree_num,tree_info) 
      info = trim(psys_info) // ',' // trim(tree_info)

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_all_and_write_back:impl:long;
      case("full_particle,Long,full_particle,full_particle,full_particle,Monopole")
         call fdps_calc_force_all_and_write_back_l00000(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,full_particle,Quadrupole")
         call fdps_calc_force_all_and_write_back_l00001(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter")
         call fdps_calc_force_all_and_write_back_l00002(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,full_particle,DipoleGeometricCenter")
         call fdps_calc_force_all_and_write_back_l00003(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter")
         call fdps_calc_force_all_and_write_back_l00004(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,essential_particle_j,Monopole")
         call fdps_calc_force_all_and_write_back_l00005(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,essential_particle_j,Quadrupole")
         call fdps_calc_force_all_and_write_back_l00006(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter")
         call fdps_calc_force_all_and_write_back_l00007(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter")
         call fdps_calc_force_all_and_write_back_l00008(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter")
         call fdps_calc_force_all_and_write_back_l00009(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,full_particle,Monopole")
         call fdps_calc_force_all_and_write_back_l00010(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,full_particle,Quadrupole")
         call fdps_calc_force_all_and_write_back_l00011(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter")
         call fdps_calc_force_all_and_write_back_l00012(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter")
         call fdps_calc_force_all_and_write_back_l00013(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter")
         call fdps_calc_force_all_and_write_back_l00014(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,essential_particle_j,Monopole")
         call fdps_calc_force_all_and_write_back_l00015(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole")
         call fdps_calc_force_all_and_write_back_l00016(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter")
         call fdps_calc_force_all_and_write_back_l00017(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter")
         call fdps_calc_force_all_and_write_back_l00018(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter")
         call fdps_calc_force_all_and_write_back_l00019(tree_num,    &
                                                        pfunc_ep_ep, &
                                                        pfunc_ep_sp, &
                                                        psys_num,    &
                                                        dinfo_num)
      
      !--------------
      case default
         write(*,'(a)')"something wrong occurs!"
         call PS_finalize(this)
         stop 0
      end select

   end subroutine calc_force_all_and_write_back_l

   !----------------------------------------------------------
   subroutine calc_force_all_s(this,        &
                               tree_num,    &
                               pfunc_ep_ep, &
                               psys_num,    &
                               dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,psys_num,dinfo_num
      type(c_funptr), intent(in) :: pfunc_ep_ep
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize,kind=c_char) :: psys_info,tree_info,info

      call get_psys_info(this,psys_num,psys_info)
      call get_tree_info(this,tree_num,tree_info)
      info = trim(psys_info) // ',' // trim(tree_info)

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_all:impl:short;
      !--------------
      case default
         write(*,'(a)')"something wrong occurs!!"
         call PS_finalize(this)
         stop 0
      end select

   end subroutine calc_force_all_s

   !----------------------------------------------------------
   subroutine calc_force_all_l(this,        &
                               tree_num,    &
                               pfunc_ep_ep, &
                               pfunc_ep_sp, &
                               psys_num,    &
                               dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,psys_num,dinfo_num
      type(c_funptr), intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize,kind=c_char) :: psys_info,tree_info,info

      call get_psys_info(this,psys_num,psys_info)
      call get_tree_info(this,tree_num,tree_info) 
      info = trim(psys_info) // ',' // trim(tree_info)

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_all:impl:long;
      case("full_particle,Long,full_particle,full_particle,full_particle,Monopole")
         call fdps_calc_force_all_l00000(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,full_particle,Quadrupole")
         call fdps_calc_force_all_l00001(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter")
         call fdps_calc_force_all_l00002(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,full_particle,DipoleGeometricCenter")
         call fdps_calc_force_all_l00003(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter")
         call fdps_calc_force_all_l00004(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,essential_particle_j,Monopole")
         call fdps_calc_force_all_l00005(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,essential_particle_j,Quadrupole")
         call fdps_calc_force_all_l00006(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter")
         call fdps_calc_force_all_l00007(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter")
         call fdps_calc_force_all_l00008(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter")
         call fdps_calc_force_all_l00009(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,full_particle,Monopole")
         call fdps_calc_force_all_l00010(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,full_particle,Quadrupole")
         call fdps_calc_force_all_l00011(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter")
         call fdps_calc_force_all_l00012(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter")
         call fdps_calc_force_all_l00013(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter")
         call fdps_calc_force_all_l00014(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,essential_particle_j,Monopole")
         call fdps_calc_force_all_l00015(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole")
         call fdps_calc_force_all_l00016(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter")
         call fdps_calc_force_all_l00017(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter")
         call fdps_calc_force_all_l00018(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter")
         call fdps_calc_force_all_l00019(tree_num,    &
                                         pfunc_ep_ep, &
                                         pfunc_ep_sp, &
                                         psys_num,    &
                                         dinfo_num)
      
      !--------------
      case default
         write(*,'(a)')"something wrong occurs!"
         call PS_finalize(this)
         stop 0
      end select

   end subroutine calc_force_all_l

   !----------------------------------------------------------
   subroutine calc_force_making_tree_s(this,        &
                                       tree_num,    &
                                       pfunc_ep_ep, &
                                       dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,dinfo_num
      type(c_funptr), intent(in) :: pfunc_ep_ep
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize,kind=c_char) :: info

      call get_tree_info(this,tree_num,info)

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_making_tree:impl:short;
      !--------------
      case default
         write(*,'(a)')"something wrong occurs!!"
         call PS_finalize(this)
         stop 0
      end select

   end subroutine calc_force_making_tree_s

   !----------------------------------------------------------
   subroutine calc_force_making_tree_l(this,        &
                                       tree_num,    &
                                       pfunc_ep_ep, &
                                       pfunc_ep_sp, &
                                       dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,dinfo_num
      type(c_funptr), intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize,kind=c_char) :: info

      call get_tree_info(this,tree_num,info) 

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_making_tree:impl:long;
      case("Long,full_particle,full_particle,full_particle,Monopole")
         call fdps_calc_force_making_tree_l00000(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,full_particle,full_particle,Quadrupole")
         call fdps_calc_force_making_tree_l00001(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter")
         call fdps_calc_force_making_tree_l00002(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,full_particle,full_particle,DipoleGeometricCenter")
         call fdps_calc_force_making_tree_l00003(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter")
         call fdps_calc_force_making_tree_l00004(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,full_particle,essential_particle_j,Monopole")
         call fdps_calc_force_making_tree_l00005(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,full_particle,essential_particle_j,Quadrupole")
         call fdps_calc_force_making_tree_l00006(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter")
         call fdps_calc_force_making_tree_l00007(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter")
         call fdps_calc_force_making_tree_l00008(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter")
         call fdps_calc_force_making_tree_l00009(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,essential_particle_i,full_particle,Monopole")
         call fdps_calc_force_making_tree_l00010(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,essential_particle_i,full_particle,Quadrupole")
         call fdps_calc_force_making_tree_l00011(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter")
         call fdps_calc_force_making_tree_l00012(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter")
         call fdps_calc_force_making_tree_l00013(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter")
         call fdps_calc_force_making_tree_l00014(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,essential_particle_i,essential_particle_j,Monopole")
         call fdps_calc_force_making_tree_l00015(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole")
         call fdps_calc_force_making_tree_l00016(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter")
         call fdps_calc_force_making_tree_l00017(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter")
         call fdps_calc_force_making_tree_l00018(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter")
         call fdps_calc_force_making_tree_l00019(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,full_particle,full_particle,Monopole")
         call fdps_calc_force_making_tree_l00020(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,full_particle,full_particle,Quadrupole")
         call fdps_calc_force_making_tree_l00021(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,full_particle,full_particle,MonopoleGeometricCenter")
         call fdps_calc_force_making_tree_l00022(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,full_particle,full_particle,DipoleGeometricCenter")
         call fdps_calc_force_making_tree_l00023(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,full_particle,full_particle,QuadrupoleGeometricCenter")
         call fdps_calc_force_making_tree_l00024(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,full_particle,essential_particle_j,Monopole")
         call fdps_calc_force_making_tree_l00025(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,full_particle,essential_particle_j,Quadrupole")
         call fdps_calc_force_making_tree_l00026(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,full_particle,essential_particle_j,MonopoleGeometricCenter")
         call fdps_calc_force_making_tree_l00027(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,full_particle,essential_particle_j,DipoleGeometricCenter")
         call fdps_calc_force_making_tree_l00028(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,full_particle,essential_particle_j,QuadrupoleGeometricCenter")
         call fdps_calc_force_making_tree_l00029(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,essential_particle_i,full_particle,Monopole")
         call fdps_calc_force_making_tree_l00030(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,essential_particle_i,full_particle,Quadrupole")
         call fdps_calc_force_making_tree_l00031(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,essential_particle_i,full_particle,MonopoleGeometricCenter")
         call fdps_calc_force_making_tree_l00032(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,essential_particle_i,full_particle,DipoleGeometricCenter")
         call fdps_calc_force_making_tree_l00033(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,essential_particle_i,full_particle,QuadrupoleGeometricCenter")
         call fdps_calc_force_making_tree_l00034(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,essential_particle_i,essential_particle_j,Monopole")
         call fdps_calc_force_making_tree_l00035(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,essential_particle_i,essential_particle_j,Quadrupole")
         call fdps_calc_force_making_tree_l00036(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,essential_particle_i,essential_particle_j,MonopoleGeometricCenter")
         call fdps_calc_force_making_tree_l00037(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,essential_particle_i,essential_particle_j,DipoleGeometricCenter")
         call fdps_calc_force_making_tree_l00038(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      case("Long,force,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter")
         call fdps_calc_force_making_tree_l00039(tree_num,    &
                                                 pfunc_ep_ep, &
                                                 pfunc_ep_sp, &
                                                 dinfo_num)
      
      !--------------
      case default
         write(*,'(a)')"something wrong occurs!"
         call PS_finalize(this)
         stop 0
      end select

   end subroutine calc_force_making_tree_l

   !----------------------------------------------------------
   subroutine calc_force_and_write_back_s(this,        &
                                          tree_num,    &
                                          pfunc_ep_ep, &
                                          psys_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,psys_num
      type(c_funptr), intent(in) :: pfunc_ep_ep
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize,kind=c_char) :: psys_info,tree_info,info

      call get_psys_info(this,psys_num,psys_info)
      call get_tree_info(this,tree_num,tree_info)
      info = trim(psys_info) // ',' // trim(tree_info)

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_and_write_back:impl:short;
      !--------------
      case default
         write(*,'(a)')"something wrong occurs!!"
         call PS_finalize(this)
         stop 0
      end select

   end subroutine calc_force_and_write_back_s

   !----------------------------------------------------------
   subroutine calc_force_and_write_back_l(this,        &
                                          tree_num,    &
                                          pfunc_ep_ep, &
                                          pfunc_ep_sp, &
                                          psys_num)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num,psys_num
      type(c_funptr), intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize,kind=c_char) :: psys_info,tree_info,info

      call get_psys_info(this,psys_num,psys_info)
      call get_tree_info(this,tree_num,tree_info) 
      info = trim(psys_info) // ',' // trim(tree_info)

      select case (trim(info))
      !--------------
      ! fdps-autogen:calc_force_and_write_back:impl:long;
      case("full_particle,Long,full_particle,full_particle,full_particle,Monopole")
         call fdps_calc_force_and_write_back_l00000(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,full_particle,full_particle,Quadrupole")
         call fdps_calc_force_and_write_back_l00001(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter")
         call fdps_calc_force_and_write_back_l00002(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,full_particle,full_particle,DipoleGeometricCenter")
         call fdps_calc_force_and_write_back_l00003(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter")
         call fdps_calc_force_and_write_back_l00004(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,full_particle,essential_particle_j,Monopole")
         call fdps_calc_force_and_write_back_l00005(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,full_particle,essential_particle_j,Quadrupole")
         call fdps_calc_force_and_write_back_l00006(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter")
         call fdps_calc_force_and_write_back_l00007(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter")
         call fdps_calc_force_and_write_back_l00008(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter")
         call fdps_calc_force_and_write_back_l00009(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,full_particle,Monopole")
         call fdps_calc_force_and_write_back_l00010(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,full_particle,Quadrupole")
         call fdps_calc_force_and_write_back_l00011(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter")
         call fdps_calc_force_and_write_back_l00012(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter")
         call fdps_calc_force_and_write_back_l00013(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter")
         call fdps_calc_force_and_write_back_l00014(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,essential_particle_j,Monopole")
         call fdps_calc_force_and_write_back_l00015(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole")
         call fdps_calc_force_and_write_back_l00016(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter")
         call fdps_calc_force_and_write_back_l00017(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter")
         call fdps_calc_force_and_write_back_l00018(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      case("full_particle,Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter")
         call fdps_calc_force_and_write_back_l00019(tree_num,    &
                                                    pfunc_ep_ep, &
                                                    pfunc_ep_sp, &
                                                    psys_num)
      
      !--------------
      case default
         write(*,'(a)')"something wrong occurs!"
         call PS_finalize(this)
         stop 0
      end select

   end subroutine calc_force_and_write_back_l

   !----------------------------------------------------------
   ! [Comment] A place where the implementations of
   !           get_neighbor_list*() are generated.
   ! fdps-autogen:get_neighbor_list:impl;
   subroutine get_neighbor_list00000(this,tree_num,pos,r_search,num_epj,fptr_to_EPJ)
      use user_defined_types
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num
      type(fdps_f64vec), intent(IN) :: pos
      real(kind=c_double), intent(IN) :: r_search
      integer(kind=c_int), intent(INOUT) :: num_epj
      type(full_particle), dimension(:), pointer, intent(INOUT) :: fptr_to_EPJ
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize,kind=c_char) :: info
      type(c_ptr) :: cptr_to_EPJ
   
      call get_tree_info(this,tree_num,info)
      select case (trim(info))
   
      case default
         write(*,100)"full_particle",trim(info)
         100 format("[Error]"/ &
                    "   The specified EssentialParticleJ does not have"/ &
                    "   getRSearch() method or the specified tree does"/ &
                    "   not support neighbor search. Please check the"/ &
                    "   definitions of the structure and the tree:"/ &
                    "   - EssentialParticleJ: ",a/ &
                    "   - TreeInfo: ",a)
         flush(6)
         call PS_abort(this)
         stop 1
      end select
   
      !* Convert C-pointer to Fortran-pointer
      call c_f_pointer(cptr_to_EPJ,fptr_to_EPJ,[num_epj])
   
   end subroutine get_neighbor_list00000
   
   subroutine get_neighbor_list00001(this,tree_num,pos,r_search,num_epj,fptr_to_EPJ)
      use user_defined_types
      use fdps_vector
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: tree_num
      type(fdps_f64vec), intent(IN) :: pos
      real(kind=c_double), intent(IN) :: r_search
      integer(kind=c_int), intent(INOUT) :: num_epj
      type(essential_particle_j), dimension(:), pointer, intent(INOUT) :: fptr_to_EPJ
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize,kind=c_char) :: info
      type(c_ptr) :: cptr_to_EPJ
   
      call get_tree_info(this,tree_num,info)
      select case (trim(info))
   
      case default
         write(*,100)"essential_particle_j",trim(info)
         100 format("[Error]"/ &
                    "   The specified EssentialParticleJ does not have"/ &
                    "   getRSearch() method or the specified tree does"/ &
                    "   not support neighbor search. Please check the"/ &
                    "   definitions of the structure and the tree:"/ &
                    "   - EssentialParticleJ: ",a/ &
                    "   - TreeInfo: ",a)
         flush(6)
         call PS_abort(this)
         stop 1
      end select
   
      !* Convert C-pointer to Fortran-pointer
      call c_f_pointer(cptr_to_EPJ,fptr_to_EPJ,[num_epj])
   
   end subroutine get_neighbor_list00001
   

   !##########################################################
   function get_rank(this)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int) :: get_rank

      get_rank = fdps_get_rank()
      
   end function get_rank

   !----------------------------------------------------------
   function get_rank_multi_dim(this,id)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int) :: get_rank_multi_dim
      integer(kind=c_int), intent(IN) :: id

      get_rank_multi_dim = fdps_get_rank_multi_dim(id)
      
   end function get_rank_multi_dim

   !----------------------------------------------------------
   function get_num_procs(this)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int) :: get_num_procs
      
      get_num_procs = fdps_get_num_procs()

   end function get_num_procs

   !----------------------------------------------------------
   function get_num_procs_multi_dim(this,id)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int) :: get_num_procs_multi_dim
      integer(kind=c_int), intent(IN) :: id
      
      get_num_procs_multi_dim = fdps_get_num_procs_multi_dim(id)

   end function get_num_procs_multi_dim

   !----------------------------------------------------------
   subroutine get_logical_and(this,f_in,f_out)
      implicit none
      class(FDPS_controller) :: this
      logical(kind=c_bool), intent(IN) :: f_in
      logical(kind=c_bool), intent(INOUT) :: f_out

      call fdps_get_logical_and(f_in,f_out)

   end subroutine get_logical_and

   !----------------------------------------------------------
   subroutine get_logical_or(this,f_in,f_out)
      implicit none
      class(FDPS_controller) :: this
      logical(kind=c_bool), intent(IN) :: f_in
      logical(kind=c_bool), intent(INOUT) :: f_out

      call fdps_get_logical_or(f_in,f_out)

   end subroutine get_logical_or

   !----------------------------------------------------------
   ! [Comment] A place where the implementations of
   !           reduction functions are generated.
   ! fdps-autogen:get_min_value:ftn_impl;
   subroutine get_min_value_i32(this,f_in,f_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: f_in
      integer(kind=c_int), intent(INOUT) :: f_out
   
       call fdps_get_min_value_i32(f_in,f_out)
   
   end subroutine get_min_value_i32
   
   subroutine get_min_value_i64(this,f_in,f_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_long_long), intent(IN) :: f_in
      integer(kind=c_long_long), intent(INOUT) :: f_out
   
       call fdps_get_min_value_i64(f_in,f_out)
   
   end subroutine get_min_value_i64
   
   subroutine get_min_value_r32(this,f_in,f_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_float), intent(IN) :: f_in
      real(kind=c_float), intent(INOUT) :: f_out
   
       call fdps_get_min_value_r32(f_in,f_out)
   
   end subroutine get_min_value_r32
   
   subroutine get_min_value_r64(this,f_in,f_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_double), intent(IN) :: f_in
      real(kind=c_double), intent(INOUT) :: f_out
   
       call fdps_get_min_value_r64(f_in,f_out)
   
   end subroutine get_min_value_r64
   
   subroutine get_min_value_w_id_r32(this,f_in,i_in,f_out,i_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_float), intent(IN) :: f_in
      real(kind=c_float), intent(INOUT) :: f_out
      integer(kind=c_int), intent(IN) :: i_in
      integer(kind=c_int), intent(INOUT) :: i_out
   
      call fdps_get_min_value_w_id_r32(f_in,i_in,f_out,i_out)
   
   end subroutine get_min_value_w_id_r32
   
   subroutine get_min_value_w_id_r64(this,f_in,i_in,f_out,i_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_double), intent(IN) :: f_in
      real(kind=c_double), intent(INOUT) :: f_out
      integer(kind=c_int), intent(IN) :: i_in
      integer(kind=c_int), intent(INOUT) :: i_out
   
      call fdps_get_min_value_w_id_r64(f_in,i_in,f_out,i_out)
   
   end subroutine get_min_value_w_id_r64
   
   ! fdps-autogen:get_max_value:ftn_impl;
   subroutine get_max_value_i32(this,f_in,f_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: f_in
      integer(kind=c_int), intent(INOUT) :: f_out
   
       call fdps_get_max_value_i32(f_in,f_out)
   
   end subroutine get_max_value_i32
   
   subroutine get_max_value_i64(this,f_in,f_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_long_long), intent(IN) :: f_in
      integer(kind=c_long_long), intent(INOUT) :: f_out
   
       call fdps_get_max_value_i64(f_in,f_out)
   
   end subroutine get_max_value_i64
   
   subroutine get_max_value_r32(this,f_in,f_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_float), intent(IN) :: f_in
      real(kind=c_float), intent(INOUT) :: f_out
   
       call fdps_get_max_value_r32(f_in,f_out)
   
   end subroutine get_max_value_r32
   
   subroutine get_max_value_r64(this,f_in,f_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_double), intent(IN) :: f_in
      real(kind=c_double), intent(INOUT) :: f_out
   
       call fdps_get_max_value_r64(f_in,f_out)
   
   end subroutine get_max_value_r64
   
   subroutine get_max_value_w_id_r32(this,f_in,i_in,f_out,i_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_float), intent(IN) :: f_in
      real(kind=c_float), intent(INOUT) :: f_out
      integer(kind=c_int), intent(IN) :: i_in
      integer(kind=c_int), intent(INOUT) :: i_out
   
      call fdps_get_max_value_w_id_r32(f_in,i_in,f_out,i_out)
   
   end subroutine get_max_value_w_id_r32
   
   subroutine get_max_value_w_id_r64(this,f_in,i_in,f_out,i_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_double), intent(IN) :: f_in
      real(kind=c_double), intent(INOUT) :: f_out
      integer(kind=c_int), intent(IN) :: i_in
      integer(kind=c_int), intent(INOUT) :: i_out
   
      call fdps_get_max_value_w_id_r64(f_in,i_in,f_out,i_out)
   
   end subroutine get_max_value_w_id_r64
   
   ! fdps-autogen:get_sum:ftn_impl;
   subroutine get_sum_i32(this,f_in,f_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: f_in
      integer(kind=c_int), intent(INOUT) :: f_out
   
       call fdps_get_sum_i32(f_in,f_out)
   
   end subroutine get_sum_i32
   
   subroutine get_sum_i64(this,f_in,f_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_long_long), intent(IN) :: f_in
      integer(kind=c_long_long), intent(INOUT) :: f_out
   
       call fdps_get_sum_i64(f_in,f_out)
   
   end subroutine get_sum_i64
   
   subroutine get_sum_r32(this,f_in,f_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_float), intent(IN) :: f_in
      real(kind=c_float), intent(INOUT) :: f_out
   
       call fdps_get_sum_r32(f_in,f_out)
   
   end subroutine get_sum_r32
   
   subroutine get_sum_r64(this,f_in,f_out)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_double), intent(IN) :: f_in
      real(kind=c_double), intent(INOUT) :: f_out
   
       call fdps_get_sum_r64(f_in,f_out)
   
   end subroutine get_sum_r64
   
   ! fdps-autogen:broadcast:ftn_impl;
   subroutine broadcast_scalar_i32(this,val,src)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(INOUT) :: val
      integer(kind=c_int), optional, intent(IN) :: src
   
       if (present(src)) then
          call fdps_broadcast_scalar_i32(val,1,src)
       else
          call fdps_broadcast_scalar_i32(val,1,0)
       end if
   
   end subroutine broadcast_scalar_i32
   
   subroutine broadcast_array_i32(this,vals,n,src)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: n
      integer(kind=c_int), dimension(n), intent(INOUT) :: vals
      integer(kind=c_int), optional, intent(IN) :: src
   
      if (present(src)) then
         call fdps_broadcast_array_i32(vals,n,src)
      else
         call fdps_broadcast_array_i32(vals,n,0)
      end if
   
   end subroutine broadcast_array_i32
   
   subroutine broadcast_scalar_i64(this,val,src)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_long_long), intent(INOUT) :: val
      integer(kind=c_int), optional, intent(IN) :: src
   
       if (present(src)) then
          call fdps_broadcast_scalar_i64(val,1,src)
       else
          call fdps_broadcast_scalar_i64(val,1,0)
       end if
   
   end subroutine broadcast_scalar_i64
   
   subroutine broadcast_array_i64(this,vals,n,src)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: n
      integer(kind=c_long_long), dimension(n), intent(INOUT) :: vals
      integer(kind=c_int), optional, intent(IN) :: src
   
      if (present(src)) then
         call fdps_broadcast_array_i64(vals,n,src)
      else
         call fdps_broadcast_array_i64(vals,n,0)
      end if
   
   end subroutine broadcast_array_i64
   
   subroutine broadcast_scalar_r32(this,val,src)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_float), intent(INOUT) :: val
      integer(kind=c_int), optional, intent(IN) :: src
   
       if (present(src)) then
          call fdps_broadcast_scalar_r32(val,1,src)
       else
          call fdps_broadcast_scalar_r32(val,1,0)
       end if
   
   end subroutine broadcast_scalar_r32
   
   subroutine broadcast_array_r32(this,vals,n,src)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: n
      real(kind=c_float), dimension(n), intent(INOUT) :: vals
      integer(kind=c_int), optional, intent(IN) :: src
   
      if (present(src)) then
         call fdps_broadcast_array_r32(vals,n,src)
      else
         call fdps_broadcast_array_r32(vals,n,0)
      end if
   
   end subroutine broadcast_array_r32
   
   subroutine broadcast_scalar_r64(this,val,src)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_double), intent(INOUT) :: val
      integer(kind=c_int), optional, intent(IN) :: src
   
       if (present(src)) then
          call fdps_broadcast_scalar_r64(val,1,src)
       else
          call fdps_broadcast_scalar_r64(val,1,0)
       end if
   
   end subroutine broadcast_scalar_r64
   
   subroutine broadcast_array_r64(this,vals,n,src)
      use, intrinsic :: iso_c_binding
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(in) :: n
      real(kind=c_double), dimension(n), intent(INOUT) :: vals
      integer(kind=c_int), optional, intent(IN) :: src
   
      if (present(src)) then
         call fdps_broadcast_array_r64(vals,n,src)
      else
         call fdps_broadcast_array_r64(vals,n,0)
      end if
   
   end subroutine broadcast_array_r64
   
   !----------------------------------------------------------

   !----------------------------------------------------------
   function get_wtime(this)
      implicit none
      class(FDPS_controller) :: this 
      real(kind=c_double) :: get_wtime

      get_wtime = fdps_get_wtime()

   end function get_wtime

   !##########################################################
   subroutine MT_init_genrand(this,s)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int), intent(IN) :: s

      call fdps_mt_init_genrand(s)
      
   end subroutine MT_init_genrand

   !----------------------------------------------------------
   function MT_genrand_int31(this)
      implicit none
      class(FDPS_controller) :: this
      integer(kind=c_int) :: MT_genrand_int31

      MT_genrand_int31 = fdps_mt_genrand_int31()

   end function MT_genrand_int31

   !----------------------------------------------------------
   function MT_genrand_real1(this)
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_double) :: MT_genrand_real1

      MT_genrand_real1 = fdps_mt_genrand_real1()

   end function MT_genrand_real1

   !----------------------------------------------------------
   function MT_genrand_real2(this)
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_double) :: MT_genrand_real2

      MT_genrand_real2 = fdps_mt_genrand_real2()

   end function MT_genrand_real2

   !----------------------------------------------------------
   function MT_genrand_real3(this)
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_double) :: MT_genrand_real3

      MT_genrand_real3 = fdps_mt_genrand_real3()

   end function MT_genrand_real3

   !----------------------------------------------------------
   function MT_genrand_res53(this)
      implicit none
      class(FDPS_controller) :: this
      real(kind=c_double) :: MT_genrand_res53

      MT_genrand_res53 = fdps_mt_genrand_res53()

   end function MT_genrand_res53

end module FDPS_module
