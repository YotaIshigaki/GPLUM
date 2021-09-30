program main
   use, intrinsic :: iso_c_binding
   use fdps_matrix
   implicit none
   !* Local variables
   type(fdps_f32mat) :: m1,m2,m3

   m1%xx = 2.0d0
   m1%xy = 3.0d0
   m1%xz = 5.0d0
   m1%yy = 4.0d0
   m1%yz = 3.0d0
   m1%zz = 1.0d0

   m2%xx = 3.0d0
   m2%xy = 1.0d0
   m2%xz = 2.0d0
   m2%yy = 4.0d0
   m2%yz = 2.0d0
   m2%zz = 1.0d0

   m3 = m1 * m2

   write(*,*)'xx = ',m3%xx
   write(*,*)'xy = ',m3%xy
   write(*,*)'xz = ',m3%xz
   write(*,*)'yy = ',m3%yy
   write(*,*)'yz = ',m3%yz
   write(*,*)'zz = ',m3%zz

   stop 0
end program main
