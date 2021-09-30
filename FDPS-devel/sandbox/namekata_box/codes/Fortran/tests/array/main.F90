program main
   implicit none
   !* Local variables
   double precision, dimension(2) :: a(1),a2
   double precision, dimension(2) :: b(4)
   double precision, dimension(2) :: c(2,2)
   double precision, dimension(2,2) :: d(3)

   a  = 1.0d0
   a2 = 1.0d0
   b  = 4.0d0
   c  = 2.0d0
   d  = 3.0d0
   write(*,*)'--------'
   write(*,*)a
   write(*,*)'--------'
   write(*,*)a2
   write(*,*)'--------'
   write(*,*)b
   write(*,*)'--------'
   write(*,*)c
   write(*,*)'--------'
   write(*,*)d
   write(*,*)'--------'

   stop 
end program main
