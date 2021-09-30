program main
   implicit none
   !* Local parameters
   integer, parameter :: real32 = kind(1.0)
   integer, parameter :: real64 = kind(1.0d0)
   double precision, parameter :: eps=1.0d-10
   !* Local variables
   double precision :: x,y,cen
   double precision :: xrel,yrel
   real :: xf,yf
   real :: xrelf,yrelf
   !* External function
   real :: cast

   x = 1.0d0
   y = x + eps
   xf = cast(x)
   yf = cast(y)
   cen = 0.5d0*(x+y)
   xrel = x - cen
   yrel = y - cen
   xrelf = cast(xrel)
   yrelf = cast(yrel)

   !* Output
   write(*,*)'-------'
   write(*,10)'x     =',x
   write(*,10)'y     =',y
   write(*,10)'diff. =',y-x
   write(*,*)'-------'
   write(*,10)'xf    =',xf
   write(*,10)'yf    =',yf
   write(*,10)'diff. =',yf-xf
   write(*,*)'-------'
   write(*,10)'xrel  =',xrel
   write(*,10)'yrel  =',yrel
   write(*,10)'diff. =',yrel-xrel
   write(*,*)'-------'
   write(*,10)'xrelf =',xrelf
   write(*,10)'yrelf =',yrelf
   write(*,10)'diff. =',yrelf-xrelf
   write(*,*)'-------'

   10 format(a,1x,1es25.15e3)

   stop 0
end program main

real function cast(x) result(ret)
   implicit none
   double precision, intent(in) :: x
   !* Local variables
   real :: xf
   xf = real(x)
   ret = xf
end function cast

