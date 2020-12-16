program main
   use mpi
   implicit none
   integer :: ierr
   integer :: n_proc,my_rank
   integer :: parent_group,group,comm
   integer :: n_proc_X, my_rank_X
   integer :: tmp
   integer, dimension(0:1) :: ranks

   call MPI_Init(ierr)
   call MPI_Comm_size(MPI_COMM_WORLD,n_proc,ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
   if (n_proc /= 4) then
      if (my_rank == 0) write(*,*)'Run with 4 MPI processes!'
      call MPI_Finalize(ierr)
      stop 0
   end if

   ! This code assumes the configuration below:
   ! 
   !     ---------
   !     | 2 | 3 |
   !  Y   ---------
   !     | 0 | 1 |
   !     ---------
   !         X
   !

   call MPI_Comm_group(MPI_COMM_WORLD,parent_group,ierr)
   ! Create groups along the X direction
   select case (my_rank)
      case(0:1)
         ranks(0) = 0
         ranks(1) = 1
      case(2:3)
         ranks(0) = 2
         ranks(1) = 3
   end select
   call MPI_Group_incl(parent_group, 2, ranks, group, ierr)
   ! Create a new communicator
   call MPI_Comm_create(MPI_COMM_WORLD, group, comm, ierr)
   ! Check the new communicator
   ! (1) check n_proc, my_rank
   call MPI_Comm_size(comm,n_proc_X,ierr)
   call MPI_Comm_rank(comm,my_rank_X,ierr)
   write(*,*)'my_rank = ',my_rank,' n_proc_X = ',n_proc_X,' my_rank_X = ',my_rank_X
   ! (2) check MPI_Bcast
   tmp = my_rank
   call MPI_Bcast(tmp,1,MPI_Integer,0,comm,ierr)
   write(*,*)'[MPI_Bcast] my_rank = ',my_rank,' my_rank_X = ',my_rank_X,' tmp = ',tmp
   ! (3) check MPI_Allreduce
   tmp = my_rank
   call MPI_Allreduce(MPI_IN_PLACE,tmp,1,MPI_Integer,MPI_SUM,comm,ierr)
   write(*,*)'[MPI_Allreduce] my_rank = ',my_rank,' my_rank_X = ',my_rank_X,' tmp = ',tmp
   



   call MPI_Finalize(ierr)

   stop 0
end program main
