program main
   implicit none
   !* Local variables
   integer :: x,y,z,x2,y2

   x = 1.0d0
   y = 2.0d0
   x2 = 3.0d0
   y2 = 3.0d0
   z = x & 
   & + y &
   & + x2 &
   & + y2 

   stop
end program main
