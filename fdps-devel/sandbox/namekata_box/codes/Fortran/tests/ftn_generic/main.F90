!========================
! MODULE: test module
!========================
module testmod

   !* Public routine
   public :: calc_x

   !* Private routines
   private :: calc_x_func
   private :: calc_x_sub

   !* Polymorphism
   interface calc_x
      module procedure calc_x_func
      module procedure calc_x_sub
   end interface calc_x

   contains

   function calc_x_func() 
      implicit none
      double precision :: calc_x_func
      calc_x_func = 1.0d0
   end function calc_x_func

   subroutine calc_x_sub(x)
      implicit none
      double precision, intent(INOUT) :: x
      x = 1.0d0
   end subroutine calc_x_sub

end module testmod

!-----------------------------------------------------------------------
!/////////////////////// M A I N   R O U T I N E ///////////////////////
!-----------------------------------------------------------------------
program main
   use testmod
   implicit none
   !* Local variables
   double precision :: x

   write(*,*)'x = ',calc_x()
   call calc_x(x)
   write(*,*)'x = ',x

   stop
end program main
