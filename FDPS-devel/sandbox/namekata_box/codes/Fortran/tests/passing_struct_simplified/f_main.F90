!-----------------------------------------------------------------------
!/////////////////////// M A I N   R O U T I N E ///////////////////////
!-----------------------------------------------------------------------
subroutine f_main()
   use, intrinsic :: iso_c_binding
   implicit none
   !* Local structure
   type, bind(c) :: full_particle
      real(c_double) :: mass
      real(c_double), dimension(3) :: pos
   end type full_particle
   !* Local parameters
   integer, parameter :: nptcl=4
   !* Local variables
   integer :: i
   type(c_ptr) :: cptr_to_FP_array
   type(full_particle), dimension(:), pointer :: fptr_to_FP_array
   !* External routines
   external :: load_array_cpp

   !* Receive FP data from C++
   call load_array_cpp(cptr_to_FP_array)
   call c_f_pointer(cptr_to_FP_array,fptr_to_FP_array,[nptcl])
   do i=1,nptcl
      write(*,*)'-------------------------------'
      write(*,*)'i      = ',i
      write(*,*)'mass   = ',fptr_to_FP_array(i)%mass
      write(*,*)'pos(1) = ',fptr_to_FP_array(i)%pos(1)
      write(*,*)'pos(2) = ',fptr_to_FP_array(i)%pos(2)
      write(*,*)'pos(3) = ',fptr_to_FP_array(i)%pos(3)
   end do
   
   stop
end subroutine f_main

