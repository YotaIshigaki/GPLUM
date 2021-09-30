!================================
!   MODULE : user_defined
!================================
module user_defined_types
   use, intrinsic :: iso_c_binding
   implicit none

   !* Full particle type
   type, bind(c) :: full_particle
      real(c_double) :: mass
      real(c_double) :: eps
      real(c_double), dimension(3) :: pos
      real(c_double), dimension(3) :: vel
      real(c_double) :: pot
      real(c_double), dimension(3) :: acc
   end type full_particle

   !* Essential particle I type
   type, bind(c) :: essential_particle_i
      real(c_double) :: mass
      real(c_double) :: eps
      real(c_double), dimension(3) :: pos
   end type essential_particle_i

   !* Essential particle J type
   type, bind(c) :: essential_particle_j
      real(c_double) :: mass
      real(c_double), dimension(3) :: pos
   end type essential_particle_j

   !* Force class
   type, bind(c) :: force
      real(c_double) :: pot
      real(c_double), dimension(3) :: acc
   end type force

end module user_defined_types


!===============================
!  MODULE: FDPS interface
!===============================
module fdps_ifc
   use, intrinsic :: iso_c_binding
   implicit none

   ! define interface of c function.
   interface
      subroutine calc_force_all(func,ep_i,n_ip,ep_j,n_jp,f) bind(c)
         use, intrinsic :: iso_c_binding
         use user_defined_types
         type(c_funptr), intent(in), value :: func
         integer(c_int), intent(in) :: n_ip,n_jp
         type(essential_particle_i), dimension(n_ip), intent(in) :: ep_i
         type(essential_particle_j), dimension(n_jp), intent(in) :: ep_j
         type(force), dimension(n_ip), intent(inout) :: f
      end subroutine calc_force_all
   end interface

   contains

   ! define procedure passed to c function.
   ! it must be interoperable!
   subroutine calc_gravity(ep_i,n_ip,ep_j,n_jp,f) bind(c)
      use user_defined_types
      integer(c_int), intent(in) :: n_ip,n_jp
      type(essential_particle_i), dimension(n_ip), intent(in) :: ep_i
      type(essential_particle_j), dimension(n_jp), intent(in) :: ep_j
      type(force), dimension(n_ip), intent(inout) :: f
      !* Local variables
      integer(c_int) :: i,j,k
      real(c_double) :: eps2,poti,r3_inv,r_inv
      real(c_double), dimension(3) :: xi,ai,rij

      do i=1,n_ip
         eps2 = ep_i(i)%eps * ep_i(i)%eps
         xi   = ep_i(i)%pos
         ai   = 0.0d0
         poti = 0.0d0
         do j=1,n_jp 
            do k=1,3
               rij(k) = xi(k) - ep_j(j)%pos(k)
            end do
            r3_inv = rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3) + eps2
            r_inv  = 1.0d0/sqrt(r3_inv)
            r3_inv = r_inv * r_inv
            r_inv  = r_inv * ep_j(j)%mass
            r3_inv = r3_inv * r_inv
            do k=1,3
               ai(k) = ai(k) - r3_inv * rij(k)
            end do
            poti = poti - r_inv
         end do
         do k=1,3
            f(i)%acc(k) = f(i)%acc(k) + ai(k)
         end do
         f(i)%pot = f(i)%pot + poti
      end do

   end subroutine calc_gravity

end module fdps_ifc

!-----------------------------------------------------------------------
!/////////////////////// M A I N   R O U T I N E ///////////////////////
!-----------------------------------------------------------------------
subroutine f_main()
   use fdps_ifc
   use user_defined_types
   implicit none
   !* Local parameters
   integer(c_int), parameter :: nptcl=4096
   real(c_double), parameter :: rmax=3.0d0,r2max=rmax*rmax
   !* Local variables
   integer(c_int) :: i,j,k
   integer(c_int) :: n_ip,n_jp
   double precision :: r,r2,a
   double precision, dimension(3) :: pos,acc
   type(essential_particle_i), dimension(:), allocatable :: ep_i
   type(essential_particle_j), dimension(:), allocatable :: ep_j
   type(force), dimension(:), allocatable :: f
   type(c_funptr) :: cproc

   !* Initialize ep_i 
   n_ip = nptcl
   allocate( ep_i(n_ip) ) 
   do i=1,n_ip
      ep_i(i)%mass = 1.0d0/n_ip
      ep_i(i)%eps  = 1.0d0/32.0d0
      do 
         call random_number(pos)
         pos = (2.0d0*pos - 1.0d0) * rmax
         r2 = pos(1)*pos(1) + pos(2)*pos(2) + pos(3)*pos(3)
         if (r2 < r2max) exit
      end do
      ep_i(i)%pos = pos
   end do
   !* Initialize ep_j
   n_jp = n_ip
   allocate( ep_j(n_jp) )
   do j=1,n_jp
      ep_j(j)%mass = ep_i(j)%mass
      ep_j(j)%pos  = ep_i(j)%pos
   end do
   !* Initialize f
   allocate( f(n_ip) )
   do i=1,n_ip
      f(i)%pot = 0.0d0
      f(i)%acc = 0.0d0
   end do

   !* Get c procedure pointer.
   cproc = c_funloc(calc_gravity)

   !* Compute force
   call calc_force_all(cproc,ep_i,n_ip,ep_j,n_jp,f)

   !* Output
   open(unit=9,file='result.txt',action='write',status='replace')
      do i=1,n_ip
         pos = ep_i(i)%pos
         acc = f(i)%acc
         r = sqrt(pos(1)*pos(1) + pos(2)*pos(2) + pos(3)*pos(3))
         a = sqrt(acc(1)*acc(1) + acc(2)*acc(2) + acc(3)*acc(3))
         write(9,'(2es25.16e3)')r,a
      end do
   close(unit=9)

   stop 0
end subroutine f_main
