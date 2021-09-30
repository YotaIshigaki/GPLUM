!-----------------------------------------------------------------------
!//////////////////          F U N C T I O N          //////////////////
!////////////////// < G E T _ T I M E _ O F _ F F T > //////////////////
!-----------------------------------------------------------------------
function get_time_of_fft(nc,npuy,icommy,npuz,icommz) 
   implicit none
   integer, intent(in) :: nc,npuy,icommy,npuz,icommz
   double precision :: get_time_of_fft
   ! Local parameters
   integer, parameter :: n_trials = 32
   ! Local variables
   integer :: i,j,k,ierr
   integer :: cnt_start,cnt_end,cnt_per_sec,cnt_max,diff
   double precision, dimension(2) :: x
   complex(kind(0d0)), dimension(:,:,:), allocatable :: a,b
  
   ! Set input data
   allocate( a(nc, nc/npuy, nc/npuz), &
             b(nc, nc/npuy, nc/npuz) )
   do k=1,nc/npuz
      do j=1,nc/npuy
         do i=1,nc
            call random_number(x)
            a(i,j,k) = cmplx(x(1), x(2))
         end do
      end do
   end do  

   ! Perform DFT
   call system_clock(cnt_start) 
   do i=1,n_trials
       call pzfft3dv(a,b,nc,nc,nc,icommy,icommz,npuy,npuz,-1)
       call pzfft3dv(a,b,nc,nc,nc,icommy,icommz,npuy,npuz,1)
   end do
   call system_clock(cnt_end, cnt_per_sec, cnt_max)
   if ( cnt_end < cnt_start ) then
      diff = (cnt_max - cnt_start) + cnt_end + 1
   else
      diff = cnt_end - cnt_start
   end if
   get_time_of_fft = diff/dble(cnt_per_sec)/dble(n_trials)

   ! Free
   deallocate(a, b)

end function get_time_of_fft

!-----------------------------------------------------------------------
!/////////////////////// M A I N   R O U T I N E ///////////////////////
!-----------------------------------------------------------------------
program main
   use mpi
   implicit none
   ! Local parameters
   integer, parameter :: cap_np_list=64 
   integer, parameter :: cap_nc_list=64
   integer, parameter :: nc_start=128, nc_end=512
   ! Local variables
   integer :: i,j,k,n,m,ierr
   integer :: n_proc, n_proc_max, my_rank
   integer :: parent_group,groupy,groupz
   integer :: tmp,pf
   integer :: npuy,npuz,iy,iz
   integer :: icommy,icommz
   integer :: nc
   integer :: size_np_list, size_nc_list
   integer :: ptr
   double precision :: etime
   integer, dimension(2) :: npu_cand
   integer, dimension(cap_np_list) :: np_list
   integer, dimension(cap_nc_list) :: nc_list
   integer, dimension(:), allocatable :: ranksy,ranksz
   double precision, dimension(:,:), allocatable :: obs_data
   ! External
   double precision :: get_time_of_fft

   ! Initialize MPI
   call MPI_Init(ierr)
   call MPI_Comm_size(MPI_COMM_WORLD,n_proc_max,ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)

   ! Make a list of n_proc
   size_np_list = 0
   n_proc = 4
   do
      size_np_list = size_np_list + 1
      if (size_np_list > cap_np_list) then
         if (my_rank == 0) then
            write(*,*)'cap_np_list is too small!'
         end if
         call MPI_Abort(MPI_COMM_WORLD,1,ierr)
         stop 1
      end if
      np_list(size_np_list) = n_proc
      n_proc = 2*n_proc
      if (n_proc > n_proc_max) exit
   end do

   ! Make a list of nc
   size_nc_list = 0
   nc = nc_start
   do 
      size_nc_list = size_nc_list + 1
      if (size_nc_list > cap_nc_list) then
         if (my_rank == 0) then
            write(*,*)'cap_nc_list is too small!'
         end if
         call MPI_Abort(MPI_COMM_WORLD,1,ierr)
         stop 1
      end if
      nc_list(size_nc_list) = nc
      nc = 2*nc
      if (nc > nc_end) exit
   end do

   ! Allocate 
   allocate( obs_data(size_nc_list, size_np_list) )

   do n=1,size_np_list
      n_proc = np_list(n)
      ! Output
      if (my_rank == 0) then
         write(*,*)'Measuring the case of n_proc = ',n_proc
      end if
      ! Set npuy & npuz
      npu_cand(1) = 1
      npu_cand(2) = 1
      tmp = n_proc
      ptr = 1
      pf = 2
      pf_outer: do
         ! Repeat a division until tmp cannot be divided by pf
         pf_inner: do
            if (mod(tmp,pf) == 0) then
               ! Update tmp
               tmp = tmp/pf
               ! Update npu_cand()
               npu_cand(ptr) = npu_cand(ptr)*pf
               ! Update ptr (the next write address) 
               if (ptr == 1) then
                  ptr = 2
               else
                  ptr = 1
               end if
            else
               exit pf_inner
            end if
         end do pf_inner

         ! Update pf
         if (tmp > 1) then
            pf = pf + 1
            if (pf == tmp) then
               npu_cand(ptr) = npu_cand(ptr)*pf
               exit pf_outer
            end if
         else
            exit pf_outer
         end if
      end do pf_outer
      npuy = npu_cand(1)
      npuz = npu_cand(2)
      if (npuy * npuz /= n_proc) then
         write(*,*)'prime factorization is something wrong!'
         call MPI_Abort(MPI_COMM_WORLD,1,ierr)
         stop 1
      end if
      if (my_rank == 0) then
         write(*,*)'npuy = ',npuy
         write(*,*)'npuz = ',npuz
      end if
      ! Set iy & iz
      iz = my_rank / npuy
      iy = my_rank - npuy * iz
      !write(*,*)'my_rank = ',my_rank,' iy = ',iy,' iz = ',iz
      ! Get the parent group
      call MPI_Comm_group(MPI_COMM_WORLD,parent_group,ierr)
      ! Create groups for communication along y direction
      allocate( ranksy(0:npuy-1) )
      do j=0,npuy-1
         ranksy(j)=npuy*iz + j
      end do
      call MPI_Group_incl(parent_group, npuy, ranksy, groupy, ierr)
      if (my_rank < n_proc) then
         call MPI_Comm_create(MPI_COMM_WORLD, groupy, icommy, ierr)
      else
         call MPI_Comm_create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, icommy, ierr)
      end if
      ! Create a group for communication along y direction
      allocate( ranksz(0:npuz-1) )
      do k=0,npuz-1
         ranksz(k)=npuy*k + iy
      end do
      call MPI_Group_incl(parent_group, npuz, ranksz, groupz, ierr)
      if (my_rank < n_proc) then
          call MPI_Comm_create(MPI_COMM_WORLD, groupz, icommz, ierr)
      else
          call MPI_Comm_create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, icommz, ierr)
      end if
      ! Measure the speed of FFT
      if (my_rank < n_proc) then
         ! Perform measurement
         do m=1,size_nc_list
            nc = nc_list(m)
            if (my_rank == 0) write(*,*)'nc = ',nc
            if (npuy <= nc .and. npuz <= nc) then
               etime = get_time_of_FFT(nc,npuy,icommy,npuz,icommz)
               if (my_rank == 0) then
                  obs_data(m,n) = etime
               end if
            end if
         end do
      end if
      ! Free the communicators and the groups
      call MPI_Group_free(groupy,ierr)
      if (icommy /= MPI_COMM_NULL) call MPI_Comm_free(icommy,ierr)
      deallocate( ranksy )
      call MPI_Group_free(groupz,ierr)
      if (icommz /= MPI_COMM_NULL) call MPI_Comm_free(icommz,ierr)
      deallocate( ranksz )
   end do

   ! Output
   if (my_rank == 0) then
      do i=1,size_nc_list
         nc = nc_list(i)
         write(*,100)nc, nc*nc*nc
         100 format("# NC = ",i10,"    NC3 = ",i10/ &
                    "# (n_proc, etime)")
         do j=1,size_np_list
            n_proc = np_list(j)
            write(*,200)n_proc,obs_data(i,j)
            200 format(i10,4x,1es25.16e3)
         end do
      end do
   end if

   ! Free
   deallocate( obs_data )

   ! Finalize MPI
   call MPI_Finalize(ierr)

   stop 0
end program main
