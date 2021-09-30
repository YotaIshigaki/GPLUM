!===============================
!   MODULE: User defined types
!===============================
module user_defined_types
   use, intrinsic :: iso_c_binding
   use fdps_vector
   use fdps_super_particle
   implicit none

   !**** Full particle type
   type, public, bind(c) :: full_particle !$fdps FP,EPI,EPJ,Force
      !$fdps copyFromForce full_particle (pot,pot) (acc,acc)
      !$fdps copyFromFP full_particle (id,id) (mass,mass) (eps,eps) (pos,pos) 
      !$fdps clear id=keep, mass=keep, eps=keep, pos=keep, vel=keep
      integer(kind=c_long_long) :: id
      real(kind=c_double) :: mass !$fdps charge
      real(kind=c_double) :: eps
      type(fdps_f64vec) :: pos !$fdps position
      type(fdps_f64vec) :: vel !$fdps velocity
      real(kind=c_double) :: pot
      type(fdps_f64vec) :: acc
   end type full_particle

   contains

   !**** Interaction function (particle-particle)
   subroutine calc_gravity_pp(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      integer(c_int), intent(in), value :: n_ip,n_jp
      type(full_particle), dimension(n_ip), intent(in) :: ep_i
      type(full_particle), dimension(n_jp), intent(in) :: ep_j
      type(full_particle), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(c_int) :: i,j
      real(c_double) :: eps2,poti,r3_inv,r_inv
      type(fdps_f64vec) :: xi,ai,rij

      write(*,'(a)')'$'
      write(*,'(a)')'$ EPI-EPJ interaction [start]'
      write(*,'(a,i5)')'$ n_ip = ',n_ip
      write(*,'(a,i5)')'$ n_jp = ',n_jp
      write(*,'(a)')'$'

      !* Compute force
      do i=1,n_ip
         eps2 = ep_i(i)%eps * ep_i(i)%eps
         xi%x = ep_i(i)%pos%x
         xi%y = ep_i(i)%pos%y
         xi%z = ep_i(i)%pos%z
         ai%x = 0.0d0
         ai%y = 0.0d0
         ai%z = 0.0d0
         poti = 0.0d0
         !***** Debug (1)
        !if (ep_i(i)%id == 10) then
            write(*,'(a)')'########## EPJ int #####'
            write(*,'(a,i5)')'id = ',ep_i(i)%id
            write(*,'(a,1es25.16e3)')'x = ',xi%x
            write(*,'(a,1es25.16e3)')'y = ',xi%y
            write(*,'(a,1es25.16e3)')'z = ',xi%z
        !end if
         do j=1,n_jp
            rij%x  = xi%x - ep_j(j)%pos%x
            rij%y  = xi%y - ep_j(j)%pos%y
            rij%z  = xi%z - ep_j(j)%pos%z
            r3_inv = rij%x*rij%x &
                   + rij%y*rij%y &
                   + rij%z*rij%z &
                   + eps2
            r_inv  = 1.0d0/sqrt(r3_inv)
            r3_inv = r_inv * r_inv
            r_inv  = r_inv * ep_j(j)%mass
            r3_inv = r3_inv * r_inv
            ai%x   = ai%x - r3_inv * rij%x
            ai%y   = ai%y - r3_inv * rij%y
            ai%z   = ai%z - r3_inv * rij%z
            poti   = poti - r_inv
            !* Debug
           !if (ep_i(i)%id == 10) then
               write(*,'(a)')'---------- EPJ int -----'
               write(*,'(a,i5)')'j  = ',j
               write(*,'(a,i5)')'id = ',ep_j(j)%id
               write(*,'(a,1es25.16e3)')'mj     = ',ep_j(j)%mass
               write(*,'(a,1es25.16e3)')'xj     = ',ep_j(j)%pos%x
               write(*,'(a,1es25.16e3)')'yj     = ',ep_j(j)%pos%y
               write(*,'(a,1es25.16e3)')'zj     = ',ep_j(j)%pos%z
               write(*,'(a,1es25.16e3)')'rij%x  = ',rij%x
               write(*,'(a,1es25.16e3)')'rij%y  = ',rij%y
               write(*,'(a,1es25.16e3)')'rij%z  = ',rij%z
               write(*,'(a,1es25.16e3)')'r3_inv = ',r3_inv
           !end if
         end do
         f(i)%pot   = f(i)%pot   + poti
         f(i)%acc%x = f(i)%acc%x + ai%x
         f(i)%acc%y = f(i)%acc%y + ai%y
         f(i)%acc%z = f(i)%acc%z + ai%z
        !if (ep_i(i)%id == 10) then
            write(*,'(a,1es25.16e3)')'fp = ',f(i)%pot
            write(*,'(a,1es25.16e3)')'fx = ',f(i)%acc%x
            write(*,'(a,1es25.16e3)')'fy = ',f(i)%acc%y
            write(*,'(a,1es25.16e3)')'fz = ',f(i)%acc%z
        !end if
      end do

      write(*,'(a)')'$'
      write(*,'(a)')'$ EPI-SPJ interaction [end]'
      write(*,'(a)')'$'

   end subroutine calc_gravity_pp

   !**** Interaction function (particle-super particle)
   subroutine calc_gravity_psp(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      integer(c_int), intent(in), value :: n_ip,n_jp
      type(full_particle), dimension(n_ip), intent(in) :: ep_i
      type(fdps_spj_monopole), dimension(n_jp), intent(in) :: ep_j
      type(full_particle), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(c_int) :: i,j
      real(c_double) :: eps2,poti,r3_inv,r_inv
      type(fdps_f64vec) :: xi,ai,rij

      write(*,'(a)')'$'
      write(*,'(a)')'$ EPI-SPJ interaction [start]' 
      write(*,'(a,i5)')'$ n_ip = ',n_ip
      write(*,'(a,i5)')'$ n_jp = ',n_jp
      write(*,'(a)')'$'

      do i=1,n_ip
         eps2 = ep_i(i)%eps * ep_i(i)%eps
         xi%x = ep_i(i)%pos%x
         xi%y = ep_i(i)%pos%y
         xi%z = ep_i(i)%pos%z
         ai%x = 0.0d0
         ai%y = 0.0d0
         ai%z = 0.0d0
         poti = 0.0d0
         !***** Debug (1)
        !if (ep_i(i)%id == 10) then
            write(*,'(a)')'########## SPJ int #####'
            write(*,'(a,i5)')'id = ',ep_i(i)%id
            write(*,'(a,1es25.16e3)')'x = ',xi%x
            write(*,'(a,1es25.16e3)')'y = ',xi%y
            write(*,'(a,1es25.16e3)')'z = ',xi%z
        !end if
         do j=1,n_jp
            rij%x  = xi%x - ep_j(j)%pos%x
            rij%y  = xi%y - ep_j(j)%pos%y
            rij%z  = xi%z - ep_j(j)%pos%z
            r3_inv = rij%x*rij%x &
                   + rij%y*rij%y &
                   + rij%z*rij%z &
                   + eps2
            r_inv  = 1.0d0/sqrt(r3_inv)
            r3_inv = r_inv * r_inv
            r_inv  = r_inv * ep_j(j)%mass
            r3_inv = r3_inv * r_inv
            ai%x   = ai%x - r3_inv * rij%x
            ai%y   = ai%y - r3_inv * rij%y
            ai%z   = ai%z - r3_inv * rij%z
            poti   = poti - r_inv
            !**** Debug (2)
           !if (ep_i(i)%id == 10) then
               write(*,'(a)')'---------- SPJ int -----'
               write(*,'(a,i5)')'j  = ',j
               write(*,'(a,1es25.16e3)')'mj     = ',ep_j(j)%mass
               write(*,'(a,1es25.16e3)')'xj     = ',ep_j(j)%pos%x
               write(*,'(a,1es25.16e3)')'yj     = ',ep_j(j)%pos%y
               write(*,'(a,1es25.16e3)')'zj     = ',ep_j(j)%pos%z
               write(*,'(a,1es25.16e3)')'rij%x  = ',rij%x
               write(*,'(a,1es25.16e3)')'rij%y  = ',rij%y
               write(*,'(a,1es25.16e3)')'rij%z  = ',rij%z
               write(*,'(a,1es25.16e3)')'r3_inv = ',r3_inv
           !end if
!*******************************************
! xi%x and ep_j(j)%pos%x return NaN
            if ( rij%x /= rij%x ) then
               write(*,*) "Error in psp"
               write(*,'(a,1p2e12.3)') "x-pos", xi%x, ep_j(j)%pos%x
               write(*,'(a,1p2e12.3)') "y-pos", xi%y, ep_j(j)%pos%y
               write(*,'(a,1p2e12.3)') "z-pos", xi%z, ep_j(j)%pos%z
               stop
            end if
!*******************************************
         end do
         f(i)%pot   = f(i)%pot   + poti
         f(i)%acc%x = f(i)%acc%x + ai%x
         f(i)%acc%y = f(i)%acc%y + ai%y
         f(i)%acc%z = f(i)%acc%z + ai%z
        !if (ep_i(i)%id == 10) then
            write(*,'(a,1es25.16e3)')'fp = ',f(i)%pot
            write(*,'(a,1es25.16e3)')'fx = ',f(i)%acc%x
            write(*,'(a,1es25.16e3)')'fy = ',f(i)%acc%y
            write(*,'(a,1es25.16e3)')'fz = ',f(i)%acc%z
        !end if
      end do

      write(*,'(a)')'$'
      write(*,'(a)')'$ EPI-SPJ interaction [end]'
      write(*,'(a)')'$'

   end subroutine calc_gravity_psp

end module user_defined_types
