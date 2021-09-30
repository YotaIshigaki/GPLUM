!==================================
!    MODULE: FDPS module
!==================================
module FDPS_module
   use, intrinsic :: iso_c_binding
   implicit none

   !* Private parameters
   integer, parameter :: int16 = selected_int_kind(4)
   integer, parameter :: int32 = selected_int_kind(9)
   integer, parameter :: int64 = selected_int_kind(18)
   integer, parameter :: real32 = selected_real_kind(6)
   integer, parameter :: real64 = selected_real_kind(15)

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
      !-(ParticleSystem functions)
      procedure :: create_psys
      procedure :: delete_psys
      procedure :: init_psys
      procedure :: get_psys_info
      procedure :: set_nptcl_loc
      procedure :: get_nptcl_loc
      procedure :: get_psys_fptr
      procedure :: exchange_particle
      !-(DomainInfo functions)
      procedure :: create_dinfo
      procedure :: delete_dinfo
      procedure :: init_dinfo
      procedure :: decompose_domain_all
      !-(Tree functions)
      procedure :: create_tree
      procedure :: delete_tree
      procedure :: init_tree
      procedure :: get_tree_info
      procedure, private :: calc_force_all_and_write_back_s
      procedure, private :: calc_force_all_and_write_back_l
      generic   :: calc_force_all_and_write_back => calc_force_all_and_write_back_s, &
                                                    calc_force_all_and_write_back_l
      !-(MPI comm. functions)
      procedure :: comm_get_rank 
      procedure :: comm_get_num_procs
     !procedure :: comm_get_min_value
     !procedure :: comm_get_max_value
     !procedure :: comm_get_sum
     !procedure :: comm_broadcast
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
   !-(ParticleSystem functions)
   private :: create_psys
   private :: delete_psys
   private :: init_psys
   private :: get_psys_info
   private :: get_nptcl_loc
   private :: set_nptcl_loc
   private :: get_psys_fptr
   !-(DomainInfo functions)
   private :: create_dinfo
   private :: delete_dinfo
   private :: init_dinfo
   private :: decompose_domain_all
   !-(Tree functions)
   private :: create_tree
   private :: delete_tree
   private :: init_tree
   private :: get_tree_info
   private :: calc_force_all_and_write_back_s
   private :: calc_force_all_and_write_back_l
   !-(MPI comm. functions)
   private :: comm_get_rank 
   private :: comm_get_num_procs
  !private :: comm_get_min_value
  !private :: comm_get_max_value
  !private :: comm_get_sum
  !private :: comm_broadcast
   !-(Utility functions)
   private :: MT_init_genrand
   private :: MT_genrand_int31
   private :: MT_genrand_real1
   private :: MT_genrand_real2
   private :: MT_genrand_res53

   !* Polymorphism
!   private :: comm_get_min_value_i32
!   private :: comm_get_min_value_r32
!   private :: comm_get_min_value_r64
!   private :: comm_get_max_value_i32
!   private :: comm_get_max_value_r32
!   private :: comm_get_max_value_r64
!   private :: comm_get_sum_i32
!   private :: comm_get_sum_r32
!   private :: comm_get_sum_r64
!   private :: comm_broadcast_i32
!   private :: comm_broadcast_r32
!   private :: comm_broadcast_r64
!
!   interface comm_get_min_value
!      module procedure comm_get_min_value_i32
!      module procedure comm_get_min_value_r32
!      module procedure comm_get_min_value_r64
!   end interface comm_get_min_value
!
!   interface comm_get_max_value
!      module procedure comm_get_max_value_i32
!      module procedure comm_get_max_value_r32
!      module procedure comm_get_max_value_r64
!   end interface comm_get_max_value
!
!   interface comm_get_sum
!      module procedure comm_get_sum_i32
!      module procedure comm_get_sum_r32
!      module procedure comm_get_sum_r64
!   end interface comm_get_sum
!
!   interface comm_broadcast
!      module procedure comm_broadcast_i32
!      module procedure comm_broadcast_r32
!      module procedure comm_broadcast_r64
!   end interface comm_broadcast

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
         integer(kind=c_int), intent(in) :: psys_num
      end subroutine fdps_delete_psys

      subroutine fdps_init_psys(psys_num)  bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(in) :: psys_num
      end subroutine fdps_init_psys

      subroutine fdps_get_psys_info(psys_num,psys_info,charlen)  bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(in) :: psys_num
         character(kind=c_char), dimension(*), intent(inout) :: psys_info
         integer(kind=c_size_t), intent(inout) :: charlen
      end subroutine fdps_get_psys_info

      subroutine fdps_set_nptcl_loc(psys_num,nptcl) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(in) :: psys_num,nptcl
      end subroutine fdps_set_nptcl_loc

      function fdps_get_nptcl_loc(psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_get_nptcl_loc
         integer(kind=c_int), intent(in) :: psys_num
      end function fdps_get_nptcl_loc

      subroutine fdps_get_psys_cptr(psys_num,cptr) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(in) :: psys_num
         type(c_ptr), intent(inout) :: cptr
      end subroutine fdps_get_psys_cptr

      subroutine fdps_exchange_particle(psys_num,dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(in) :: psys_num,dinfo_num
      end subroutine fdps_exchange_particle

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
         integer(kind=c_int), intent(in) :: dinfo_num 
      end subroutine fdps_delete_dinfo

      subroutine fdps_init_dinfo(dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(in) :: dinfo_num 
      end subroutine fdps_init_dinfo

      subroutine fdps_decompose_domain_all(dinfo_num,psys_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(in) :: dinfo_num,psys_num
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
         integer(kind=c_int), intent(in) :: tree_num
      end subroutine fdps_delete_tree

      subroutine fdps_init_tree(tree_num,nptcl) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(in) :: tree_num,nptcl
      end subroutine fdps_init_tree

      subroutine fdps_get_tree_info(tree_num,tree_info,charlen) bind(c) 
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(in) :: tree_num
         character(kind=c_char), dimension(*), intent(inout) :: tree_info
         integer(kind=c_size_t), intent(inout) :: charlen
      end subroutine fdps_get_tree_info

      ! [Comment] A place where the interface statements for 
      !           fdps_calc_force_all_and_write_back_*() 
      !           are generated.
      ! fdps-auto-gen-key:calcForceAllAndWriteBack():ftn_ifc
      subroutine fdps_calc_force_all_and_write_back_s(tree_num,    &
                                                      pfunc_ep_ep, &
                                                      psys_num,    &
                                                      dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), intent(in), value :: pfunc_ep_ep
      end subroutine fdps_calc_force_all_and_write_back_s

      subroutine fdps_calc_force_all_and_write_back_l(tree_num,    &
                                                      pfunc_ep_ep, &
                                                      pfunc_ep_sp, &
                                                      psys_num,    &
                                                      dinfo_num) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(in) :: tree_num,psys_num,dinfo_num
         type(c_funptr), intent(in), value :: pfunc_ep_ep,pfunc_ep_sp
      end subroutine fdps_calc_force_all_and_write_back_l

      !----------------------
      !  MPI comm. 
      !----------------------
      function fdps_comm_get_rank() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_comm_get_rank
      end function fdps_comm_get_rank

      function fdps_comm_get_num_procs() bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int) :: fdps_comm_get_num_procs
      end function fdps_comm_get_num_procs

      !----------------------
      !  Utility
      !----------------------
      subroutine fdps_mt_init_genrand(s) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), intent(in) :: s 
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

   !##########################################################
   subroutine create_psys(this,psys_num,psys_info_in)
      implicit none
      class(FDPS_controller) :: this
      integer, intent(INOUT) :: psys_num
      character(len=*,kind=c_char), intent(IN) :: psys_info_in
      character(len=:,kind=c_char), allocatable :: psys_info

      psys_info = trim(psys_info_in) // c_null_char
      call fdps_create_psys(psys_num,psys_info)

   end subroutine create_psys

   !----------------------------------------------------------
   subroutine delete_psys(this,psys_num)
      implicit none
      class(FDPS_controller) :: this
      integer, intent(IN) :: psys_num

      call fdps_delete_psys(psys_num)

   end subroutine delete_psys

   !----------------------------------------------------------
   subroutine init_psys(this,psys_num)
      implicit none
      class(FDPS_controller) :: this
      integer, intent(IN) :: psys_num

      call fdps_init_psys(psys_num)

   end subroutine init_psys

   !----------------------------------------------------------
   subroutine get_psys_info(this,psys_num,psys_info)
      implicit none
      class(FDPS_controller) :: this
      integer, intent(IN) :: psys_num
      character(len=*), intent(INOUT) :: psys_info
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
   subroutine set_nptcl_loc(this,psys_num,nptcl)
      implicit none
      class(FDPS_controller) :: this
      integer, intent(IN) :: psys_num,nptcl

      call fdps_set_nptcl_loc(psys_num,nptcl)

   end subroutine set_nptcl_loc

   !----------------------------------------------------------
   function get_nptcl_loc(this,psys_num)
      implicit none
      integer :: get_nptcl_loc
      class(FDPS_controller) :: this
      integer, intent(IN) :: psys_num

      get_nptcl_loc =  fdps_get_nptcl_loc(psys_num)
      
   end function get_nptcl_loc

   !----------------------------------------------------------
   subroutine get_psys_fptr(this,psys_num,fptr_to_FP)
      use user_defined_types
      implicit none
      class(FDPS_controller) :: this
      integer, intent(IN) :: psys_num
      type(full_particle), dimension(:), pointer, intent(INOUT) :: fptr_to_FP
      !* Local variables
      integer :: nptcl_loc
      type(c_ptr) :: cptr_to_FP

      call fdps_get_psys_cptr(psys_num,cptr_to_FP)
      nptcl_loc = fdps_get_nptcl_loc(psys_num) 
      call c_f_pointer(cptr_to_FP,fptr_to_FP,[nptcl_loc])

   end subroutine get_psys_fptr

   !----------------------------------------------------------
   subroutine exchange_particle(this,psys_num,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer, intent(IN) :: psys_num,dinfo_num

      call fdps_exchange_particle(psys_num,dinfo_num)

   end subroutine exchange_particle

   !##########################################################
   subroutine create_dinfo(this,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer, intent(INOUT) :: dinfo_num

      call fdps_create_dinfo(dinfo_num)

   end subroutine create_dinfo

   !----------------------------------------------------------
   subroutine delete_dinfo(this,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer, intent(IN) :: dinfo_num

      call fdps_delete_dinfo(dinfo_num)

   end subroutine delete_dinfo

   !----------------------------------------------------------
   subroutine init_dinfo(this,dinfo_num)
      implicit none
      class(FDPS_controller) :: this
      integer, intent(IN) :: dinfo_num

      call fdps_init_dinfo(dinfo_num)

   end subroutine init_dinfo

   !----------------------------------------------------------
   subroutine decompose_domain_all(this,dinfo_num,psys_num)
      implicit none
      class(FDPS_controller) :: this
      integer, intent(IN) :: dinfo_num,psys_num

      call fdps_decompose_domain_all(dinfo_num,psys_num)

   end subroutine decompose_domain_all

   !##########################################################
   subroutine create_tree(this,tree_num,tree_info_in)
      implicit none
      class(FDPS_controller) :: this
      integer, intent(INOUT) :: tree_num
      character(len=*,kind=c_char), intent(IN) :: tree_info_in
      character(len=:,kind=c_char), allocatable :: tree_info

      tree_info = trim(tree_info_in) // c_null_char
      call fdps_create_tree(tree_num,tree_info)

   end subroutine create_tree

   !----------------------------------------------------------
   subroutine delete_tree(this,tree_num)
      implicit none
      class(FDPS_controller) :: this
      integer, intent(IN) :: tree_num

      call fdps_delete_tree(tree_num)

   end subroutine delete_tree

   !----------------------------------------------------------
   subroutine init_tree(this,tree_num,nptcl)
      implicit none
      class(FDPS_controller) :: this
      integer, intent(IN) :: tree_num,nptcl

      call fdps_init_tree(tree_num,nptcl)

   end subroutine init_tree

   !----------------------------------------------------------
   subroutine get_tree_info(this,tree_num,tree_info)
      implicit none
      class(FDPS_controller) :: this
      integer, intent(IN) :: tree_num
      character(len=*), intent(INOUT) :: tree_info
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
   subroutine calc_force_all_and_write_back_s(this,        &
                                              tree_num,    &
                                              pfunc_ep_ep, &
                                              psys_num,    &
                                              dinfo_num)
      use user_defined_types
      implicit none
      class(FDPS_controller) :: this
      integer, intent(IN) :: tree_num,psys_num,dinfo_num
      type(c_funptr), intent(in) :: pfunc_ep_ep
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize) :: tree_info

      call get_tree_info(this,tree_num,tree_info) 

      ! fdps-auto-gen-key:calcForceAllAndWriteBack():impl:short
      select case (trim(tree_info))
      case("Short,full_particle,full_particle,full_particle,Gather")
         call fdps_calc_force_all_and_write_back_s(tree_num,    &
                                                   pfunc_ep_ep, &
                                                   psys_num,    &
                                                   dinfo_num)
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
      use user_defined_types
      implicit none
      class(FDPS_controller) :: this
      integer, intent(IN) :: tree_num,psys_num,dinfo_num
      type(c_funptr), intent(in) :: pfunc_ep_ep,pfunc_ep_sp
      !* Local parameters
      integer, parameter :: bufsize=256
      !* Local variables
      character(len=bufsize) :: tree_info

      call get_tree_info(this,tree_num,tree_info) 

      ! fdps-auto-gen-key:calcForceAllAndWriteBack():impl:long
      select case (trim(tree_info))
      case ("Long,full_particle,full_particle,full_particle,Monopole")
         call fdps_calc_force_all_and_write_back_l(tree_num,    &
                                                   pfunc_ep_ep, &
                                                   pfunc_ep_sp, &
                                                   psys_num,    &
                                                   dinfo_num)
      case default
         write(*,'(a)')"something wrong occurs!"
         call PS_finalize(this)
         stop 0
      end select

   end subroutine calc_force_all_and_write_back_l


   !##########################################################
   function comm_get_rank(this)
      implicit none
      class(FDPS_controller) :: this
      integer :: comm_get_rank

      comm_get_rank = fdps_comm_get_rank()
      
   end function comm_get_rank

   !----------------------------------------------------------
   function comm_get_num_procs(this)
      implicit none
      class(FDPS_controller) :: this
      integer :: comm_get_num_procs
      
      comm_get_num_procs = fdps_comm_get_num_procs()

   end function comm_get_num_procs

!   !----------------------------------------------------------
!   function comm_get_min_value_i32(this,val)
!      implicit none
!      class(FDPS_controller) :: this
!      integer, intent(IN) :: val
!      integer :: comm_get_min_value_i32
!      !* External routines
!      integer, external :: fdps_comm_get_min_value_i32
!
!      comm_get_min_value_i32 = fdps_comm_get_min_value_i32()
!
!   end function comm_get_min_value_i32
!
!   !----------------------------------------------------------
!   function comm_get_min_value_r32(this,val)
!      implicit none
!      class(FDPS_controller) :: this
!      real, intent(IN) :: val
!      real :: comm_get_min_value_r32
!      !* External routines
!      real, external :: fdps_comm_get_min_value_r32
!
!      comm_get_min_value_r32 = fdps_comm_get_min_value_r32()
!
!   end function comm_get_min_value_r32
!
!   !----------------------------------------------------------
!   function comm_get_min_value_r64(this,val)
!      implicit none
!      class(FDPS_controller) :: this
!      double precision, intent(IN) :: val
!      double precision :: comm_get_min_value_r64
!      !* External routines
!      double precision, external :: fdps_comm_get_min_value_r64
!
!      comm_get_min_value_r64 = fdps_comm_get_min_value_r64()
!
!   end function comm_get_min_value_r64
!
!   !----------------------------------------------------------
!   function comm_get_max_value_i32(this,val)
!      implicit none
!      class(FDPS_controller) :: this
!      integer, intent(IN) :: val
!      integer :: comm_get_max_value_i32
!      !* External routines
!      integer, external :: fdps_comm_get_max_value_i32
!
!      comm_get_max_value_i32 = fdps_comm_get_max_value_i32()
!
!   end function comm_get_max_value_i32
!
!   !----------------------------------------------------------
!   function comm_get_max_value_r32(this,val)
!      implicit none
!      class(FDPS_controller) :: this
!      real, intent(IN) :: val
!      real :: comm_get_max_value_r32
!      !* External routines
!      real, external :: fdps_comm_get_max_value_r32
!
!      comm_get_max_value_r32 = fdps_comm_get_max_value_r32()
!
!   end function comm_get_max_value_r32
!
!   !----------------------------------------------------------
!   function comm_get_max_value_r64(this,val)
!      implicit none
!      class(FDPS_controller) :: this
!      double precision, intent(IN) :: val
!      double precision :: comm_get_max_value_r64
!      !* External routines
!      double precision, external :: fdps_comm_get_max_value_r64
!
!      comm_get_max_value_r64 = fdps_comm_get_max_value_r64()
!
!   end function comm_get_max_value_r64
!
!   !----------------------------------------------------------
!   function comm_get_sum_i32(this,val)
!      implicit none
!      class(FDPS_controller) :: this
!      integer, intent(IN) :: val
!      integer :: comm_get_sum_i32
!      !* External routines
!      integer, external :: fdps_comm_get_sum_i32
!
!      comm_get_sum_i32 = fdps_comm_get_sum_i32()
!
!   end function comm_get_sum_i32
!
!   !----------------------------------------------------------
!   function comm_get_sum_r32(this,val)
!      implicit none
!      class(FDPS_controller) :: this
!      real, intent(IN) :: val
!      real :: comm_get_sum_r32
!      !* External routines
!      real, external :: fdps_comm_get_sum_r32
!
!      comm_get_sum_r32 = fdps_comm_get_sum_r32()
!
!   end function comm_get_sum_r32
!
!   !----------------------------------------------------------
!   function comm_get_sum_r64(this,val)
!      implicit none
!      class(FDPS_controller) :: this
!      double precision, intent(IN) :: val
!      double precision :: comm_get_sum_r64
!      !* External routines
!      double precision, external :: fdps_comm_get_sum_r64
!
!      comm_get_sum_r64 = fdps_comm_get_sum_r64()
!
!   end function comm_get_sum_r64
!
!   !----------------------------------------------------------
!   function comm_broadcast_i32(this,val)
!      implicit none
!      class(FDPS_controller) :: this
!      integer, intent(IN) :: val
!      integer :: comm_broadcast_i32
!      !* External routines
!      integer, external :: fdps_comm_broadcast_i32
!
!      comm_broadcast_i32 = fdps_comm_broadcast_i32()
!
!   end function comm_broadcast_i32
!
!   !----------------------------------------------------------
!   function comm_broadcast_r32(this,val)
!      implicit none
!      class(FDPS_controller) :: this
!      real, intent(IN) :: val
!      real :: comm_broadcast_r32
!      !* External routines
!      real, external :: fdps_comm_broadcast_r32
!
!      comm_broadcast_r32 = fdps_comm_broadcast_r32()
!
!   end function comm_broadcast_r32
!
!   !----------------------------------------------------------
!   function comm_broadcast_r64(this,val)
!      implicit none
!      class(FDPS_controller) :: this
!      double precision, intent(IN) :: val
!      double precision :: comm_broadcast_r64
!      !* External routines
!      double precision, external :: fdps_comm_broadcast_r64
!
!      comm_broadcast_r64 = fdps_comm_broadcast_r64()
!
!   end function comm_broadcast_r64

   !##########################################################
   subroutine MT_init_genrand(this,s)
      implicit none
      class(FDPS_controller) :: this
      integer, intent(IN) :: s

      call fdps_mt_init_genrand(s)
      
   end subroutine MT_init_genrand

   !----------------------------------------------------------
   function MT_genrand_int31(this)
      implicit none
      class(FDPS_controller) :: this
      integer :: MT_genrand_int31

      MT_genrand_int31 = fdps_mt_genrand_int31()

   end function MT_genrand_int31

   !----------------------------------------------------------
   function MT_genrand_real1(this)
      implicit none
      class(FDPS_controller) :: this
      double precision :: MT_genrand_real1

      MT_genrand_real1 = fdps_mt_genrand_real1()

   end function MT_genrand_real1

   !----------------------------------------------------------
   function MT_genrand_real2(this)
      implicit none
      class(FDPS_controller) :: this
      double precision :: MT_genrand_real2

      MT_genrand_real2 = fdps_mt_genrand_real2()

   end function MT_genrand_real2

   !----------------------------------------------------------
   function MT_genrand_real3(this)
      implicit none
      class(FDPS_controller) :: this
      double precision :: MT_genrand_real3

      MT_genrand_real3 = fdps_mt_genrand_real3()

   end function MT_genrand_real3

   !----------------------------------------------------------
   function MT_genrand_res53(this)
      implicit none
      class(FDPS_controller) :: this
      double precision :: MT_genrand_res53

      MT_genrand_res53 = fdps_mt_genrand_res53()

   end function MT_genrand_res53

end module FDPS_module
