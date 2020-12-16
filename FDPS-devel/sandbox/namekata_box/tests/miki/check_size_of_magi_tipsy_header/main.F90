program main
   use iso_c_binding
   implicit none
   type, bind(c) :: magi_tipsy_header
      real(kind=c_double) :: time
      integer(kind=c_int) :: nbodies
      integer(kind=c_int) :: ndim
      integer(kind=c_int) :: nsph
      integer(kind=c_int) :: ndark
      integer(kind=c_int) :: nstar
   end type magi_tipsy_header
   !* Local variables
   type(magi_tipsy_header) :: header
  
   open(unit=9,file='magi_tipsy_header.dat',action='write', &
        form='unformatted',access='stream',status='replace')
      write(9)header
   close(unit=9)

   stop 0
end program main
