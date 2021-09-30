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
      real(kind=c_double)  mass !$fdps charge
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
     !type(fdps_f64vec) :: xi,ai,rij
     !real(c_double), dimension(3) :: xi,xj,ai,rij
      real(c_double), dimension(3) :: xi,ai,rij
      real(c_double), dimension(n_jp) :: mj
      real(c_double), dimension(3,n_jp) :: xj

#ifdef INTERACTION_FUNCTION_IS_EMPTY
!     Do nothing
#else
      !* Compute force
     !do i=1,n_ip
     !   eps2 = ep_i(i)%eps * ep_i(i)%eps
     !   xi%x = ep_i(i)%pos%x
     !   xi%y = ep_i(i)%pos%y
     !   xi%z = ep_i(i)%pos%z
     !   ai%x = 0.0d0
     !   ai%y = 0.0d0
     !   ai%z = 0.0d0
     !   poti = 0.0d0
     !   do j=1,n_jp
     !      rij%x  = xi%x - ep_j(j)%pos%x
     !      rij%y  = xi%y - ep_j(j)%pos%y
     !      rij%z  = xi%z - ep_j(j)%pos%z
     !      r3_inv = rij%x*rij%x &
     !             + rij%y*rij%y &
     !             + rij%z*rij%z &
     !             + eps2
     !      r_inv  = 1.0d0/sqrt(r3_inv)
     !      r3_inv = r_inv * r_inv
     !      r_inv  = r_inv * ep_j(j)%mass
     !      r3_inv = r3_inv * r_inv
     !      ai%x   = ai%x - r3_inv * rij%x
     !      ai%y   = ai%y - r3_inv * rij%y
     !      ai%z   = ai%z - r3_inv * rij%z
     !      poti   = poti - r_inv
     !      ! [IMPORTANT NOTE]
     !      !   In the innermost loop, we use the components of vectors
     !      !   directly for vector operations because of the following
     !      !   reasion. Except for intel compilers with `-ipo` option,
     !      !   most of Fortran compilers use function calls to perform
     !      !   vector operations like rij = x - ep_j(j)%pos.
     !      !   This significantly slow downs the speed of the code.
     !      !   By using the components of vector directly, we can avoid 
     !      !   these function calls.
     !   end do
     !   f(i)%pot = f(i)%pot + poti
     !   f(i)%acc%x = f(i)%acc%x + ai%x
     !   f(i)%acc%y = f(i)%acc%y + ai%y
     !   f(i)%acc%z = f(i)%acc%z + ai%z
     !end do

      do j=1,n_jp
         mj(j)   = ep_j(j)%mass
         xj(1,j) = ep_j(j)%pos%x
         xj(2,j) = ep_j(j)%pos%y
         xj(3,j) = ep_j(j)%pos%z
      end do
      do i=1,n_ip
         eps2 = ep_i(i)%eps * ep_i(i)%eps
         xi(1) = ep_i(i)%pos%x
         xi(2) = ep_i(i)%pos%y
         xi(3) = ep_i(i)%pos%z
         ai = 0.0d0
         poti = 0.0d0
         do j=1,n_jp
            rij(1) = xi(1) - xj(1,j)
            rij(2) = xi(2) - xj(2,j)
            rij(3) = xi(3) - xj(3,j)
            r3_inv = rij(1)*rij(1) &
                   + rij(2)*rij(2) &
                   + rij(3)*rij(3) &
                   + eps2
            r_inv  = 1.0d0/sqrt(r3_inv)
            r3_inv = r_inv * r_inv
            r_inv  = r_inv * mj(j)
            r3_inv = r3_inv * r_inv
            ai(1)  = ai(1) - r3_inv * rij(1)
            ai(2)  = ai(2) - r3_inv * rij(2)
            ai(3)  = ai(3) - r3_inv * rij(3)
            poti   = poti - r_inv
         end do
         f(i)%pot = f(i)%pot + poti
         f(i)%acc%x = f(i)%acc%x + ai(1)
         f(i)%acc%y = f(i)%acc%y + ai(2)
         f(i)%acc%z = f(i)%acc%z + ai(3)
      end do
#endif

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
     !type(fdps_f64vec) :: xi,ai,rij
      real(c_double), dimension(3) :: xi,ai,rij
      real(c_double), dimension(n_jp) :: mj
      real(c_double), dimension(3,n_jp) :: xj

#ifdef INTERACTION_FUNCTION_IS_EMPTY
      ! Do nothing
#else
     !do i=1,n_ip
     !   eps2 = ep_i(i)%eps * ep_i(i)%eps
     !   xi%x = ep_i(i)%pos%x
     !   xi%y = ep_i(i)%pos%y
     !   xi%z = ep_i(i)%pos%z
     !   ai%x = 0.0d0
     !   ai%y = 0.0d0
     !   ai%z = 0.0d0
     !   poti = 0.0d0
     !   do j=1,n_jp
     !      rij%x  = xi%x - ep_j(j)%pos%x
     !      rij%y  = xi%y - ep_j(j)%pos%y
     !      rij%z  = xi%z - ep_j(j)%pos%z
     !      r3_inv = rij%x*rij%x &
     !             + rij%y*rij%y &
     !             + rij%z*rij%z &
     !             + eps2
     !      r_inv  = 1.0d0/sqrt(r3_inv)
     !      r3_inv = r_inv * r_inv
     !      r_inv  = r_inv * ep_j(j)%mass
     !      r3_inv = r3_inv * r_inv
     !      ai%x   = ai%x - r3_inv * rij%x
     !      ai%y   = ai%y - r3_inv * rij%y
     !      ai%z   = ai%z - r3_inv * rij%z
     !      poti   = poti - r_inv
     !   end do
     !   f(i)%pot = f(i)%pot + poti
     !   f(i)%acc%x = f(i)%acc%x + ai%x
     !   f(i)%acc%y = f(i)%acc%y + ai%y
     !   f(i)%acc%z = f(i)%acc%z + ai%z
     !end do

      do j=1,n_jp
         mj(j)   = ep_j(j)%mass
         xj(1,j) = ep_j(j)%pos%x
         xj(2,j) = ep_j(j)%pos%y
         xj(3,j) = ep_j(j)%pos%z
      end do
      do i=1,n_ip
         eps2 = ep_i(i)%eps * ep_i(i)%eps
         xi(1) = ep_i(i)%pos%x
         xi(2) = ep_i(i)%pos%y
         xi(3) = ep_i(i)%pos%z
         ai = 0.0d0
         poti = 0.0d0
         do j=1,n_jp
            rij(1) = xi(1) - xj(1,j)
            rij(2) = xi(2) - xj(2,j)
            rij(3) = xi(3) - xj(3,j)
            r3_inv = rij(1)*rij(1) &
                   + rij(2)*rij(2) &
                   + rij(3)*rij(3) &
                   + eps2
            r_inv  = 1.0d0/sqrt(r3_inv)
            r3_inv = r_inv * r_inv
            r_inv  = r_inv * mj(j)
            r3_inv = r3_inv * r_inv
            ai(1)  = ai(1) - r3_inv * rij(1)
            ai(2)  = ai(2) - r3_inv * rij(2)
            ai(3)  = ai(3) - r3_inv * rij(3)
            poti   = poti - r_inv
         end do
         f(i)%pot = f(i)%pot + poti
         f(i)%acc%x = f(i)%acc%x + ai(1)
         f(i)%acc%y = f(i)%acc%y + ai(2)
         f(i)%acc%z = f(i)%acc%z + ai(3)
      end do
#endif

   end subroutine calc_gravity_psp

end module user_defined_types
