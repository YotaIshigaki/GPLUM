program main
   use mpi
   implicit none
   !* Local variables
   integer :: ierr
   integer :: nprocs,myrank

   call MPI_Init(ierr)
   call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)


   write(*,*)'Hello from ',myrank

   call MPI_Finalize(ierr)


   stop 0
end program main
