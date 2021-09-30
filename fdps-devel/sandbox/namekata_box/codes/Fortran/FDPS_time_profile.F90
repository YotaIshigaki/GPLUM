!================================
!   MODULE: FDPS time profile
!================================
module fdps_time_profile
   use, intrinsic :: iso_c_binding
   implicit none

   !**** PS::TimeProfile
   type, public, bind(c) :: fdps_time_prof
      real(kind=c_double) :: collect_sample_particle
      real(kind=c_double) :: decompose_domain
      real(kind=c_double) :: exchange_particle
      real(kind=c_double) :: make_local_tree
      real(kind=c_double) :: make_global_tree
      real(kind=c_double) :: calc_force
      real(kind=c_double) :: calc_moment_local_tree
      real(kind=c_double) :: calc_moment_global_tree
      real(kind=c_double) :: make_LET_1st
      real(kind=c_double) :: make_LET_2nd
      real(kind=c_double) :: exchange_LET_1st
      real(kind=c_double) :: exchange_LET_2nd
   end type fdps_time_prof

end module fdps_time_profile

