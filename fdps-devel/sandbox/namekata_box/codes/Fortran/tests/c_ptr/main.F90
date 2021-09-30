program main
   use, intrinsic :: iso_c_binding
   implicit none
   !* Local parameters
   integer, parameter :: size_of_array=7
   !* Local variables
   integer :: i
   type(c_ptr) :: cptr1, cptr2
   integer, target :: array(size_of_array), scalar
   integer, pointer :: pa(:), ps

   !* Obtain the addresses of variables in the C-lang style
   cptr1 = c_loc(array(1))
   ! The programmer needs to ensure that the
   ! array is contiguous if required by the C procedure.

   cptr2 = c_loc(scalar)

   !* Assign C pointers to Fortran pointer
   call c_f_pointer(cptr2, ps)
   call c_f_pointer(cptr1, pa, shape=[size_of_array])

   !* Set the substances
   scalar = 2
   do i=1,size_of_array
      array(i) = i
   end do

   !* Check the content of pointers
   write(*,*)'ps = ',ps
   write(*,*)'lbound(pa,1) = ',lbound(pa,1)
   write(*,*)'ubound(pa,1) = ',ubound(pa,1)
   write(*,*)'The content of pa(:) :'
   do i=1,size_of_array
      write(*,*)i,pa(i)
   end do

   stop
end program main
