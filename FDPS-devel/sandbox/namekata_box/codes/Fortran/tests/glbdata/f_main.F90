!============================
!   MODULE: 
!============================
module link_to_c_vars
   use, intrinsic :: iso_c_binding
   implicit none

   ! A test for iso_c_binding taken from
   ! http://wwweic.eri.u-tokyo.ac.jp/computer/manual/altix/compile/Fortran/Intel_Fdoc100/main_for/mergedProjects/bldaps_for/common/bldaps_interopc.htm
   integer(c_int), bind(c) :: c_extern
   integer(c_long) :: c2
   bind(c, name='myVariable') :: c2
   common /com/ r,s
   real(c_float) :: r,s,t
   bind(c) :: /com/, /single/
   common /single/ t

   ! Additional test
   real(c_double), bind(c) :: d
   real(c_double), bind(c, name="D") :: cptl_d
   real(c_double), dimension(16), bind(c) :: darray
   real(c_double), dimension(4,4), bind(c) :: array2d

   type, bind(c) :: point
      real(c_double) :: x,y,z
   end type point
   type(point), bind(c) :: p
   type(point), dimension(4), bind(c) :: parray
    
end module link_to_c_vars

!-----------------------------------------------------------------------
!/////////////////////// M A I N   R O U T I N E ///////////////////////
!-----------------------------------------------------------------------
subroutine f_main()
   use link_to_c_vars
   implicit none
   !* Local variables
   integer :: i,j

   write(*,'(a)')'***********************************'

   write(*,*)'c_extern = ',c_extern
   write(*,*)'myVariable = ',c2
   write(*,*)'com.r = ',r
   write(*,*)'com.s = ',s
   write(*,*)'single.t = ',t
   write(*,*)'d = ',d
   write(*,*)'D (cptl_d) = ',cptl_d
   do i=1,16
      write(*,100)i,darray(i)
      100 format("(i,darray) = ",i3,1x,1es25.16e3)
   end do
   do j=1,4
      do i=1,4
         write(*,200)i,j,array2d(i,j)
         200 format("(i,j,array2d) = ",2i3,1x,1es25.16e3)
      end do
   end do

   write(*,*)'p%x = ',p%x
   write(*,*)'p%y = ',p%y
   write(*,*)'p%z = ',p%z

   do i=1,4
      write(*,*)'i = ',i
      write(*,*)'parray%x = ',parray(i)%x
      write(*,*)'parray%y = ',parray(i)%y
      write(*,*)'parray%z = ',parray(i)%z
   end do
   
   stop
end subroutine f_main
