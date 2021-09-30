program main
   implicit none
   !* Local parameters
   integer, parameter :: n = 16
   !* Local variables
   integer, dimension(0:n-1) :: a
   integer, dimension(:), allocatable :: b

   allocate( b(0:n-1) )

   call sub(a)

   stop 0
end program main

subroutine sub(arr)
   implicit none
   integer, parameter :: n = 16
   integer, dimension(0:n-1), intent(in) :: arr
end subroutine sub
