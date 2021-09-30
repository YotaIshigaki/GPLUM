!===============================
!  MODULE: test module
!===============================
module testmod
   use, intrinsic :: iso_c_binding
   implicit none

   ! define interface of c function.
   interface
      integer(kind=c_int) function call_it (func, arg) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_funptr), intent(in), value :: func
         integer(kind=c_int), intent(in), value :: arg
      end function call_it
   end interface

   contains

   ! define procedure passed to c function.
   ! it must be interoperable!
   integer(kind=c_int) function double_it (arg) bind(c)
     integer(kind=c_int), intent(in), value :: arg
     double_it = arg + arg
   end function double_it

   ! call c function.
   subroutine foobar ()
     type(c_funptr) :: cproc
     integer(kind=c_int) :: i

     ! get c procedure pointer.
     cproc = c_funloc (double_it)

     ! use it.
     do i = 1_c_int, 10_c_int
       print *, call_it (cproc, i)
     end do
   end subroutine foobar

end module testmod

!-----------------------------------------------------------------------
!/////////////////////// M A I N   R O U T I N E ///////////////////////
!-----------------------------------------------------------------------
program main
   use testmod
   implicit none

   call foobar()

   stop
end program main
