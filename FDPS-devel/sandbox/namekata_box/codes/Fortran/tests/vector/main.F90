program main
   use, intrinsic :: iso_c_binding
   use fdps_vector
   implicit none
   !* Local variables
   type(fdps_f64vec) :: v1,v2,v3
   real(kind=c_double), dimension(3) :: a
   real(kind=c_double) :: norm

   a = 3.0d0

   v1 = a
   v2 = v1
   v3 = v1 + v2
   norm = v3 * v3

   write(*,*)'v1: ',v1%x,v1%y,v1%z
   write(*,*)'v2: ',v2%x,v2%y,v2%z
   write(*,*)'v3: ',v3%x,v3%y,v3%z
   write(*,*)'norm = ',norm
   write(*,*)'loc(v1) = ',loc(v1)
   write(*,*)'loc(v1%x) = ',loc(v1%x)
   write(*,*)'loc(v1%y) = ',loc(v1%y)
   write(*,*)'loc(v1%z) = ',loc(v1%z)

   stop 0
end program main
