!-----------------------------------------------------------------------
!/////////////////////// < M A I N  R O U T I N E > ////////////////////
!-----------------------------------------------------------------------
subroutine f_main()
   use fdps_module
   use user_defined_types
   use io_manager
   implicit none
   !* Local parameters
   !-(force parameters)
   real, parameter :: theta = 0.5
   integer, parameter :: n_leaf_limit = 8
   integer, parameter :: n_group_limit = 64
   !-(domain decomposition)
   real, parameter :: coef_ema=0.3
   !* Local variables
   integer :: i,j,k,num_loop,ierr
   integer :: psys_num,dinfo_num,tree_num
   integer :: nptcl_loc,nptcl_glb
   logical :: clear
   double precision :: ekin0,epot0,etot0
   double precision :: ekin1,epot1,etot1
   double precision :: r,acc
   type(fdps_controller) :: fdps_ctrl
   type(full_particle), dimension(:), pointer :: ptcl
   type(c_funptr) :: pfunc_ep_ep,pfunc_ep_sp
   !-(IO)
   character(len=64) :: fname
   !-(Time measurements)
   integer :: clk_sta,clk_end,clk_rate
   double precision :: wtime,wtime_per_step

   !* Initialize FDPS
   call fdps_ctrl%PS_Initialize()

   !* Output the summary of the model parameters
   if (fdps_ctrl%get_rank() == 0) then
      call print_model_params()
   end if

   !* Create domain info object
   call fdps_ctrl%create_dinfo(dinfo_num)
   call fdps_ctrl%init_dinfo(dinfo_num,coef_ema)

   !* Create particle system object
   call fdps_ctrl%create_psys(psys_num,'full_particle')
   call fdps_ctrl%init_psys(psys_num)

   !* Make an initial condition
   call setup_IC(fdps_ctrl,psys_num)
  !call fdps_ctrl%ps_finalize()
  !stop 0

   !* Domain decomposition and exchange particle
   call fdps_ctrl%decompose_domain_all(dinfo_num,psys_num)
   call fdps_ctrl%exchange_particle(psys_num,dinfo_num)

   !* Create tree object
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)
   call fdps_ctrl%create_tree(tree_num, &
                              "Long,full_particle,full_particle,full_particle,Monopole")
   call fdps_ctrl%init_tree(tree_num,3*nptcl_loc,theta, &
                            n_leaf_limit,n_group_limit)

   !* Compute force at the initial time
   pfunc_ep_ep = c_funloc(calc_gravity_pp)
   pfunc_ep_sp = c_funloc(calc_gravity_psp)
   call fdps_ctrl%calc_force_all_and_write_back(tree_num,    &
                                                pfunc_ep_ep, &
                                                pfunc_ep_sp, &
                                                psys_num,    &
                                                dinfo_num)
   !* Compute energies at the initial time
  !clear = .true.
  !call calc_energy(fdps_ctrl,psys_num,etot0,ekin0,epot0,clear)

   !* Time integration
   num_loop = 0
   call system_clock(clk_sta)
   do
#ifdef OUTPUT_DATA 
      !* Output
      if (fdps_ctrl%get_rank() == 0) then
         write(*,50)num_loop,time_sys/time_end
         50 format('(num_loop, time_sys/time_end) = ',i10,1x,1es25.16e3)
      end if
      if ( (time_sys >= time_snap) .or. &
           (((time_sys + dt) - time_snap) > (time_snap - time_sys)) ) then
         call output(fdps_ctrl,psys_num)
         time_snap = time_snap + dt_snap
      end if
#endif
      !* Termination
      if (time_sys >= time_end) then
         exit
      end if
     
!#ifdef CHECK_ENERGY_ERROR 
!      !* Compute energies and output the results
!      clear = .true.
!      call calc_energy(fdps_ctrl,psys_num,etot1,ekin1,epot1,clear)
!      if (fdps_ctrl%get_rank() == 0) then
!         if ( (time_sys >= time_diag) .or. &
!              (((time_sys + dt) - time_diag) > (time_diag - time_sys)) ) then
!            write(*,100)time_sys,(etot1-etot0)/etot0
!            100 format("time: ",1es20.10e3,", energy error: ",1es20.10e3)
!            time_diag = time_diag + dt_diag
!         end if
!      end if
!#endif

      !* Leapfrog: Kick-Drift
      call kick(fdps_ctrl,psys_num,0.5d0*dt)
      time_sys = time_sys + dt
      call drift(fdps_ctrl,psys_num,dt)

      !* Domain decomposition & exchange particle
      call fdps_ctrl%decompose_domain_all(dinfo_num,psys_num)
      call fdps_ctrl%exchange_particle(psys_num,dinfo_num)

      !* Force calculation
      pfunc_ep_ep = c_funloc(calc_gravity_pp)
      pfunc_ep_sp = c_funloc(calc_gravity_psp)
      call fdps_ctrl%calc_force_all_and_write_back(tree_num,    &
                                                   pfunc_ep_ep, &
                                                   pfunc_ep_sp, &
                                                   psys_num,    &
                                                   dinfo_num)
      !* Leapfrog: Kick
      call kick(fdps_ctrl,psys_num,0.5d0*dt)

      !* Update num_loop
      num_loop = num_loop + 1
     !if (num_loop == 32) exit
   end do
   call system_clock(clk_end,clk_rate)

   nptcl_glb = fdps_ctrl%get_nptcl_glb(psys_num)
   if (fdps_ctrl%get_rank() == 0) then
      wtime = (clk_end-clk_sta)/dble(clk_rate)
      wtime_per_step = wtime/dble(num_loop)
      write(*,*)'ntot     = ',nptcl_glb
      write(*,*)'num_loop = ',num_loop
      write(*,*)'wtime    = ',wtime,' [s]'
      write(*,*)'wtime    = ',wtime_per_step,' [s/step]'
   end if

   !* Finalize FDPS
   call fdps_ctrl%PS_Finalize()
   stop 0

end subroutine f_main

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!///////////////////////// < S E T U P _ I C > /////////////////////////
!-----------------------------------------------------------------------
subroutine setup_IC(fdps_ctrl,psys_num)
   use fdps_vector
   use fdps_module
   use user_defined_types
   use io_manager
   implicit none
   type(fdps_controller), intent(IN) :: fdps_ctrl
   integer, intent(IN) :: psys_num
   !* Local variables
   integer :: i,j,k,ierr
   integer :: io_stat
   integer :: nprocs,myrank
   integer :: nptcl_loc
   double precision :: rmin,rmax,r2min,r2max,r2,r
   double precision :: zmin,zmax,z
   double precision :: sigma_v
   type(full_particle), dimension(:), pointer :: ptcl
   character(len=5) :: file_num,proc_num
   character(len=64) :: fname
   !* External routines
   character(len=5), external :: fnum5

   !* Get # of MPI processes and rank number
   nprocs = fdps_ctrl%get_num_procs()
   myrank = fdps_ctrl%get_rank()

#ifdef RESTART_MODE
   !* Open file 
   write(file_num,"(i5.5)")snap_num_rst
   write(proc_num,"(i5.5)")myrank
   fname =  trim(root_dir) // "/" &
         // trim(file_prefix_1st) // file_num // "-" &
         // trim(file_prefix_2nd) // proc_num // ".dat"
   open(unit=9,file=trim(fname),action='read',form='unformatted', &
        access='stream',status='old',iostat=io_stat)
      if (io_stat /= 0) then
         write(*,'(a)')'cannot open the file ' // trim(fname)
         call fdps_ctrl%ps_abort(9)
         stop 1
      end if
      read(9)snap_num
      read(9)time_diag
      read(9)time_snap
      read(9)time_sys
      read(9)nptcl_loc
      !** Set # of local particles
      call fdps_ctrl%set_nptcl_loc(psys_num,nptcl_loc)
      !** Read particle data from the file
      call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
      do i=1,nptcl_loc
         read(9)ptcl(i)%id,ptcl(i)%mass, &
                ptcl(i)%pos%x,ptcl(i)%pos%y,ptcl(i)%pos%z, &
                ptcl(i)%vel%x,ptcl(i)%vel%y,ptcl(i)%vel%z
      end do
   close(unit=9)
#else
#ifdef TWO_BODY_COLLISION_TEST
   if (myrank == 0) then
      !* Set # of local particles
      nptcl_loc = 2
      call fdps_ctrl%set_nptcl_loc(psys_num,nptcl_loc)
      
      !* Place two particles
      !** get the pointer to full particle data
      call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
      do i=1,nptcl_loc
         ptcl(i)%id    = i
         ptcl(i)%mass  = M_particle
         if (i == 1) then
            !* Left particle
            ptcl(i)%pos%x = - 10.0d0*r_particle
            ptcl(i)%vel%x =   2.0d0*r_particle/Tdur
         else
            !* Right particle
            ptcl(i)%pos%x =   10.0d0*r_particle
            ptcl(i)%vel%x = - 2.0d0*r_particle/Tdur
         end if 
         ptcl(i)%pos%y = 0.0d0
         ptcl(i)%pos%z = 0.0d0
         ptcl(i)%vel%y = 0.0d0
         ptcl(i)%vel%z = 0.0d0
      end do
      write(*,*)'Initial sep.  = ',ptcl(2)%pos%x-ptcl(1)%pos%x,' [cm]'
      write(*,*)'Initial velc. = ',ptcl(1)%vel%x,' [cm/s]'
   else
      call fdps_ctrl%set_nptcl_loc(psys_num,0)
   end if
#elif defined(TWO_BODY_COLLAPSE_TEST)
   if (myrank == 0) then
      !* Set # of local particles
      nptcl_loc = 2
      call fdps_ctrl%set_nptcl_loc(psys_num,nptcl_loc)
      
      !* Place two particles
      !** get the pointer to full particle data
      call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
      do i=1,nptcl_loc
         ptcl(i)%id    = i
         ptcl(i)%mass  = M_particle
         if (i == 1) then
            !* Left particle
            ptcl(i)%pos%x = - 10.0d0*r_particle
         else
            !* Right particle
            ptcl(i)%pos%x =   10.0d0*r_particle
         end if 
         ptcl(i)%pos%y = 0.0d0
         ptcl(i)%pos%z = 0.0d0
         ptcl(i)%vel%x = 0.0d0
         ptcl(i)%vel%y = 0.0d0
         ptcl(i)%vel%z = 0.0d0
      end do
      write(*,*)'Initial sep.  = ',ptcl(2)%pos%x-ptcl(1)%pos%x,' [cm]'
   else
      call fdps_ctrl%set_nptcl_loc(psys_num,0)
   end if
#elif defined(KEPLER_ROTATION_TEST)
   if (myrank == 0) then
      !* Set # of local particles
      nptcl_loc = 2
      call fdps_ctrl%set_nptcl_loc(psys_num,nptcl_loc)
      
      !* Place two particles
      !** get the pointer to full particle data
      call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
      do i=1,nptcl_loc
         ptcl(i)%id    = i
         ptcl(i)%mass  = M_particle
         if (i == 1) then
            !* Left particle
            ptcl(i)%pos%x = - a_inner_ring
            ptcl(i)%vel%y = - dsqrt(Ggrav*M_chariklo/a_inner_ring)
         else
            !* Right particle
            ptcl(i)%pos%x =   a_inner_ring
            ptcl(i)%vel%y =   dsqrt(Ggrav*M_chariklo/a_inner_ring)
         end if 
         ptcl(i)%pos%y = 0.0d0
         ptcl(i)%pos%z = 0.0d0
         ptcl(i)%vel%x = 0.0d0
         ptcl(i)%vel%z = 0.0d0
      end do
   else
      call fdps_ctrl%set_nptcl_loc(psys_num,0)
   end if
#else
   !* Set # of local particles
   nptcl_loc = N_particle/nprocs
  !nptcl_loc = 16384/nprocs
   if (myrank == 0) write(*,*)'nptcl_loc = ',nptcl_loc
   call fdps_ctrl%set_nptcl_loc(psys_num,nptcl_loc)

   !* Create an uniform narrow ring of particles
   !** get the pointer to full particle data
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   !** initialize Mersenne twister
   call fdps_ctrl%MT_init_genrand(myrank)
   rmin = a_inner_ring - 0.5d0 * W_inner_ring
   rmax = a_inner_ring + 0.5d0 * W_inner_ring
   r2min = rmin * rmin
   r2max = rmax * rmax
  !zmin = - 0.01d0 * W_inner_ring
  !zmax =   0.01d0 * W_inner_ring
   zmin = - 2.0d0 * r_particle
   zmax =   2.0d0 * r_particle
   do i=1,nptcl_loc
      ptcl(i)%id   = i + nptcl_loc * myrank
      ptcl(i)%mass = M_particle
      do 
         ptcl(i)%pos%x = (2.0d0*fdps_ctrl%MT_genrand_res53()-1.0d0) * rmax
         ptcl(i)%pos%y = (2.0d0*fdps_ctrl%MT_genrand_res53()-1.0d0) * rmax
         ptcl(i)%pos%z = (2.0d0*fdps_ctrl%MT_genrand_res53()-1.0d0) * zmax
         r2 = ptcl(i)%pos%x * ptcl(i)%pos%x &
            + ptcl(i)%pos%y * ptcl(i)%pos%y
         z  = ptcl(i)%pos%z
         if ((r2min <= r2) .and. (r2 <= r2max) .and. &
             (zmin <= z) .and. (z <= zmax)) exit
      end do
      r = dsqrt(r2)
      ptcl(i)%vel%x = - dsqrt(Ggrav*M_chariklo/r) * (ptcl(i)%pos%y/r)
      ptcl(i)%vel%y =   dsqrt(Ggrav*M_chariklo/r) * (ptcl(i)%pos%x/r)
      ptcl(i)%vel%z = 0.0d0
   end do
   if (myrank == 0) then
      write(*,*)'zmin = ',zmin*1.0d-2,' [m]'
      write(*,*)'zmax = ',zmax*1.0d-2,' [m]'
      sigma_v = dsqrt((2.0d0*Ggrav*M_chariklo/M_particle) &
                     *(1.0d0/(a_inner_ring*a_inner_ring)  &
                      -1.0d0/(a_inner_ring*a_inner_ring &
                             +zmax*zmax)))
      write(*,*)'\sigma_{v} = ',sigma_v,' [cm/s]'
   end if
#endif
#endif

   !* Output
  !fname = 'IC' // fnum5(myrank) // '.dat'
  !open(unit=9,file=trim(fname),action='write',status='replace')
  !   do i=1,nptcl_loc
  !      write(9,'(3es25.16e3)')ptcl(i)%pos%x,ptcl(i)%pos%y,ptcl(i)%pos%z
  !   end do
  !close(unit=9)

   !* Release the pointer
   nullify( ptcl )

end subroutine setup_IC

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!/////////////////////////     < K I C K >     /////////////////////////
!-----------------------------------------------------------------------
subroutine kick(fdps_ctrl,psys_num,dt)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(IN) :: fdps_ctrl
   integer, intent(IN) :: psys_num
   double precision, intent(IN) :: dt
   !* Local variables
   integer :: i,nptcl_loc
   type(full_particle), dimension(:), pointer :: ptcl

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)

   !* Get the pointer to full particle data
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   do i=1,nptcl_loc
      ptcl(i)%vel%x = ptcl(i)%vel%x + ptcl(i)%acc%x * dt
      ptcl(i)%vel%y = ptcl(i)%vel%y + ptcl(i)%acc%y * dt
      ptcl(i)%vel%z = ptcl(i)%vel%z + ptcl(i)%acc%z * dt
   end do
  !write(*,*)'i = 1: ',ptcl(1)%acc%x 
  !write(*,*)'i = 2: ',ptcl(2)%acc%x 
  !write(*,*)'sep. = ',ptcl(2)%pos%x-ptcl(1)%pos%x
  !call fdps_ctrl%ps_finalize()
  !stop 0
   nullify(ptcl)


end subroutine kick

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!/////////////////////////    < D R I F T >    /////////////////////////
!-----------------------------------------------------------------------
subroutine drift(fdps_ctrl,psys_num,dt)
   use fdps_vector
   use fdps_module
   use user_defined_types
   implicit none
   type(fdps_controller), intent(IN) :: fdps_ctrl
   integer, intent(IN) :: psys_num
   double precision, intent(IN) :: dt
   !* Local variables
   integer :: i,nptcl_loc
   type(full_particle), dimension(:), pointer :: ptcl

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)

   !* Get the pointer to full particle data
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
   do i=1,nptcl_loc
      ptcl(i)%pos%x = ptcl(i)%pos%x + ptcl(i)%vel%x * dt
      ptcl(i)%pos%y = ptcl(i)%pos%y + ptcl(i)%vel%y * dt
      ptcl(i)%pos%z = ptcl(i)%pos%z + ptcl(i)%vel%z * dt
   end do
   nullify(ptcl)

end subroutine drift

! !-----------------------------------------------------------------------
! !//////////////////////    S U B R O U T I N E    //////////////////////
! !////////////////////// < C A L C _ E N E R G Y > //////////////////////
! !-----------------------------------------------------------------------
! subroutine calc_energy(fdps_ctrl,psys_num,etot,ekin,epot,clear)
!    use fdps_vector
!    use fdps_module
!    use user_defined_types
!    implicit none
!    type(fdps_controller), intent(IN) :: fdps_ctrl
!    integer, intent(IN) :: psys_num
!    double precision, intent(INOUT) :: etot,ekin,epot
!    logical, intent(IN) :: clear
!    !* Local variables
!    integer :: i,nptcl_loc
!    double precision :: etot_loc,ekin_loc,epot_loc
!    type(full_particle), dimension(:), pointer :: ptcl
! 
!    !* Clear energies
!    if (clear .eqv. .true.) then
!       etot = 0.0d0
!       ekin = 0.0d0
!       epot = 0.0d0
!    end if
! 
!    !* Get # of local particles
!    nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)
!    call fdps_ctrl%get_psys_fptr(psys_num,ptcl)
! 
!    !* Compute energies
!    ekin_loc = 0.0d0
!    epot_loc = 0.0d0 
!    do i=1,nptcl_loc
!       ekin_loc = ekin_loc + ptcl(i)%mass * ( ptcl(i)%vel%x * ptcl(i)%vel%x &
!                                            + ptcl(i)%vel%y * ptcl(i)%vel%y &
!                                            + ptcl(i)%vel%z * ptcl(i)%vel%z)
!       epot_loc = epot_loc + ptcl(i)%mass * (ptcl(i)%pot + ptcl(i)%mass/ptcl(i)%eps)
!    end do
!    ekin_loc = ekin_loc * 0.5d0
!    epot_loc = epot_loc * 0.5d0
!    etot_loc = ekin_loc + epot_loc
!    call fdps_ctrl%get_sum(ekin_loc,ekin)
!    call fdps_ctrl%get_sum(epot_loc,epot)
!    call fdps_ctrl%get_sum(etot_loc,etot)
! 
!    !* Release the pointer
!    nullify(ptcl)
! 
! end subroutine calc_energy

!-----------------------------------------------------------------------
!///////////////////////// S U B R O U T I N E /////////////////////////
!/////////////////////////   < O U T P U T >   /////////////////////////
!-----------------------------------------------------------------------
subroutine output(fdps_ctrl,psys_num)
   use fdps_vector
   use fdps_module
   use user_defined_types
   use io_manager
   implicit none
   type(fdps_controller), intent(IN) :: fdps_ctrl
   integer, intent(IN) :: psys_num
   !* Local variables
   integer :: i,nptcl_loc
   integer :: myrank
   character(len=5) :: file_num,proc_num
   character(len=64) :: cmd,sub_dir,fname
   type(full_particle), dimension(:), pointer :: ptcl

   !* Get the rank number
   myrank = fdps_ctrl%get_rank()

   !* Get # of local particles
   nptcl_loc = fdps_ctrl%get_nptcl_loc(psys_num)

   !* Get the pointer to full particle data 
   call fdps_ctrl%get_psys_fptr(psys_num,ptcl)

   !* Output
   write(file_num,"(i5.5)")snap_num
   write(proc_num,"(i5.5)")myrank
   fname =  trim(root_dir) // "/" &
         // trim(file_prefix_1st) // file_num // "-" &
         // trim(file_prefix_2nd) // proc_num // ".dat"
#if OUTPUT_DATA_IN_ASCII 
   open(unit=9,file=trim(fname),action='write',status='replace')
      do i=1,nptcl_loc
         write(9,100)ptcl(i)%id,ptcl(i)%mass, &
                     ptcl(i)%pos%x,ptcl(i)%pos%y,ptcl(i)%pos%z, &
                     ptcl(i)%vel%x,ptcl(i)%vel%y,ptcl(i)%vel%z
         100 format(i8,1x,7e25.16e3)
      end do
   close(unit=9)
#else
   open(unit=9,file=trim(fname),action='write',form='unformatted', &
        access='stream',status='replace')
      write(9)snap_num+1 ! for the next output
      write(9)time_diag
      write(9)time_snap+dt_snap ! for the next output
      write(9)time_sys
      write(9)nptcl_loc
      do i=1,nptcl_loc
         write(9)ptcl(i)%id,ptcl(i)%mass, &
                     ptcl(i)%pos%x,ptcl(i)%pos%y,ptcl(i)%pos%z, &
                     ptcl(i)%vel%x,ptcl(i)%vel%y,ptcl(i)%vel%z
      end do
   close(unit=9)
#endif
   nullify(ptcl)

   !* Update snap_num
   snap_num = snap_num + 1

end subroutine output

!-----------------------------------------------------------------------
!////////////////////////  S U B R O U T I N E  ////////////////////////
!////////////////////////     < F N U M 5 >     ////////////////////////
!-----------------------------------------------------------------------
character(len=5) function fnum5(n) result(ret)
   implicit none
   integer, intent(IN) :: n

   write(ret,"(i5.5)")n

end function fnum5
