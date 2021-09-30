!===============================
!  MODULE: user_defined_types
!===============================
module user_defined_types
   use, intrinsic :: iso_c_binding
   implicit none

   !* FullParticle
   type, public, bind(c) :: full_particle
      real(c_double) :: mass
      real(c_double), dimension(3) :: pos
   end type full_particle

end module user_defined_types

!-----------------------------------------------------------------------
!////////////////////////  S U B R O U T I N E  ////////////////////////
!//////////////////////// < F T N _ C L E A R > ////////////////////////
!-----------------------------------------------------------------------
subroutine ftn_clear(ptcl) bind(c)
   use, intrinsic :: iso_c_binding
   use user_defined_types
   implicit none
   type(full_particle), intent(inout) :: ptcl

   write(*,*)'This is a test routine'

   ptcl%mass   = 1.0d0
   ptcl%pos(1) = 2.0d0
   ptcl%pos(2) = 3.0d0
   ptcl%pos(3) = 4.0d0 
   
end subroutine ftn_clear
