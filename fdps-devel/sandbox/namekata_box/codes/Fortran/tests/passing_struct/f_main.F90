!==========================
!   MODULE: data types
!==========================
module data_types
   use, intrinsic :: iso_c_binding
   implicit none

   !**** Full particle
   ! [Case I]
   type, bind(c) :: full_particle
      real(c_double) :: mass
      real(c_double), dimension(3) :: pos,vel
   end type full_particle
   ! [Case II]
  !type :: full_particle
  !   double precision :: mass
  !   double precision, dimension(3) :: pos,vel
  !end type full_particle
   ! Note that both definitions work well at least in gcc/gfortran.

   !**** Public routines
   public :: load_scalar
   public :: load_array

   contains

   !-------------------------------------------------------------------
   subroutine load_scalar(cptr_to_ptcl)
      implicit none
      type(c_ptr) :: cptr_to_ptcl
      !* External routines
      external :: load_scalar_cpp

      call load_scalar_cpp(cptr_to_ptcl)
      
   end subroutine load_scalar

   !-------------------------------------------------------------------
   subroutine load_array(cptr_to_ptcl)
      implicit none
      type(c_ptr) :: cptr_to_ptcl
      !* External routines
      external :: load_array_cpp

      call load_array_cpp(cptr_to_ptcl)
      
   end subroutine load_array


end module data_types

!-----------------------------------------------------------------------
!/////////////////////// M A I N   R O U T I N E ///////////////////////
!-----------------------------------------------------------------------
subroutine f_main()
   use data_types
   implicit none
   !* Local parameters
   integer, parameter :: nptcl=4
   !* Local variables
   integer :: i,j
   !-(load_scalar)
   type(c_ptr) :: cptr_to_FP
   type(full_particle), pointer :: fptr_to_FP
   !-(load_array)
   type(c_ptr) :: cptr_to_FP_array
   type(full_particle), dimension(:), pointer :: fptr_to_FP_array

   !* Receive FP data from C++
   call load_scalar(cptr_to_FP)
   call c_f_pointer(cptr_to_FP,fptr_to_FP)

   !* Check
   write(*,*)''
   write(*,*)'mass   = ',fptr_to_FP%mass
   write(*,*)'pos(1) = ',fptr_to_FP%pos(1)
   write(*,*)'pos(2) = ',fptr_to_FP%pos(2)
   write(*,*)'pos(3) = ',fptr_to_FP%pos(3)
   write(*,*)'vel(1) = ',fptr_to_FP%vel(1)
   write(*,*)'vel(3) = ',fptr_to_FP%vel(2)
   write(*,*)'vel(3) = ',fptr_to_FP%vel(3)

   !* Receive FP data from C++
   call load_array(cptr_to_FP_array)
   call c_f_pointer(cptr_to_FP_array,fptr_to_FP_array,[nptcl])
   do i=1,nptcl
      write(*,*)'-------------------------------'
      write(*,*)'i      = ',i
      write(*,*)'mass   = ',fptr_to_FP_array(i)%mass
      write(*,*)'pos(1) = ',fptr_to_FP_array(i)%pos(1)
      write(*,*)'pos(2) = ',fptr_to_FP_array(i)%pos(2)
      write(*,*)'pos(3) = ',fptr_to_FP_array(i)%pos(3)
      write(*,*)'vel(1) = ',fptr_to_FP_array(i)%vel(1)
      write(*,*)'vel(3) = ',fptr_to_FP_array(i)%vel(2)
      write(*,*)'vel(3) = ',fptr_to_FP_array(i)%vel(3)
   end do
   
   stop
end subroutine f_main

