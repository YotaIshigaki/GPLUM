!===============================
!   MODULE: Model parameters
!===============================
module model_parameters
   use, intrinsic :: iso_c_binding
   implicit none

   !* Mathematical constants
   real(kind=c_double), parameter :: pi=4.0d0*datan(1.0d0)

   !* Physical constants
   real(kind=c_double), parameter :: km=1.0d5
   real(kind=c_double), parameter :: yr=3.1556925168d7
   real(kind=c_double), parameter :: Ggrav=6.67408d-8 ! CODATA2014

   !* Model parameters in the cgs unit
   !-(chariklo)
   real(kind=c_double), parameter :: r_chariklo=125.0d0*km ! 125.0 [km] 
   real(kind=c_double), parameter :: rho_chariklo=1.0d0 ! 1.0 [g/cm^3]
   real(kind=c_double), parameter :: M_chariklo=(4.0d0*pi/3.0d0)*(r_chariklo**3.0d0)*rho_chariklo
   !-(inner ring)
   real(kind=c_double), parameter :: a_inner_ring=390.6d0*km ! 390.6[km]
   real(kind=c_double), parameter :: W_inner_ring=6.7d0*km ! 6.7[km]
   real(kind=c_double), parameter :: tau_inner_ring=0.38d0
   !-(outer ring)
   real(kind=c_double), parameter :: a_outer_ring=404.8d0*km ! 404.8[km]
   real(kind=c_double), parameter :: W_outer_ring=3.5d0*km ! 3.5[km]
   real(kind=c_double), parameter :: tau_outer_ring=0.06d0
   !-(particles)
   real(kind=c_double), parameter :: r_particle=5.0d2 ! 5[m]
   real(kind=c_double), parameter :: rho_particle=0.5d0*rho_chariklo
   real(kind=c_double), parameter :: M_particle=(4.0d0*pi/3.0d0)*(r_particle**3.0d0)*rho_particle
   integer(kind=c_long_long), parameter :: N_particle=int(tau_inner_ring*(2.0d0*a_inner_ring*W_inner_ring)/(r_particle*r_particle))
   !-(collision model)
   real(kind=c_double), parameter :: Tkep=2.0d0*pi*dsqrt(((a_inner_ring-0.5d0*W_inner_ring)**3.0d0)/(Ggrav*M_chariklo))
   real(kind=c_double), parameter :: Tdur=0.0025d0*Tkep
   real(kind=c_double), parameter :: e_CoR=0.1d0
   real(kind=c_double), parameter :: k_spg=(2.0d0*pi/Tdur)**2.0d0
   real(kind=c_double), parameter :: eta_dp=4.0d0*pi/(Tdur*dsqrt(1.0d0+(pi/(-dlog(e_CoR)))**2.0d0))

   !* Public routines
   public :: print_model_params

   contains

   !-------------------------------------------------------------------
   !/////////////           S U B R O U T I N E           /////////////
   !///////////// < P R I N T _ M O D E L _ P A R A M S > /////////////
   !-------------------------------------------------------------------
   subroutine print_model_params()
      implicit none
      !* Local variables
      double precision :: omega,sigma,sigma_ptcl
      double precision :: sigma_R,lambda_crit

      !* Quantities related to self-gravitational instability
      !  [Toomre 1964]
      !  For a stable stelar disk, 
      !  Q = \kappa * \sigma_{R}/(3.36*G*\Sigma) < 1
      !  should be satisfied.
      !  In other words, 
      !  \sigma_{R} > 3.36 * G * \Sigma / \kappa
      !  must be satisfied.
      !
      !  [Salo 1995]
      !  The critical wavelength for the self-gravitational
      !  instability can be written as
      !  \lambda_{crit} = 4\pi^{2} G \Sigma/\Omega^{2}
      !  for a Keplerian disk.
      !
      omega = dsqrt(Ggrav*M_chariklo/(a_inner_ring**3.0d0))
      sigma = (M_particle*N_particle) &
            / (2.0d0*a_inner_ring*W_inner_ring)
      sigma_ptcl = M_particle/(pi*r_particle*r_particle)
      sigma_R = 3.36d0*Ggrav*sigma/omega
      lambda_crit = 4.0d0*pi*pi*Ggrav*sigma/(omega*omega)

      write(*,100)r_chariklo/km, &
                  rho_chariklo, &
                  M_chariklo, &
                  a_inner_ring/km, &
                  W_inner_ring/km, &
                  tau_inner_ring, &
                  a_outer_ring/km, &
                  W_outer_ring/km, &
                  tau_outer_ring, &
                  r_particle/1.0d2, &
                  rho_particle, &
                  M_particle, &
                  N_particle, &
                  M_particle*N_particle/M_chariklo, &
                  Tkep/yr, &
                  Tdur/yr, &
                  e_CoR, &
                  k_spg, &
                  eta_dp, &
                  0.1d0/eta_dp/yr, &
                  sigma, &
                  sigma_ptcl, &
                  sigma_R, &
                  lambda_crit/W_inner_ring
      100 format("========================================"/          &
                 "(1) Chariklo"/                                      &
                 "    r_chariklo     = ",1es25.16e3," [km]"/          &
                 "    rho_chariklo   = ",1es25.16e3," [g/cm^3]"/      &
                 "    M_chariklo     = ",1es25.16e3," [g]"//          &
                 "(2) Inner ring"/                                    &
                 "    a_inner_ring   = ",1es25.16e3," [km]"/          &
                 "    W_inner_ring   = ",1es25.16e3," [km]"/          &
                 "    tau_inner_ring = ",1es25.16e3//                 &
                 "(3) Outer ring"/                                    &
                 "    a_outer_ring   = ",1es25.16e3," [km]"/          &
                 "    W_outer_ring   = ",1es25.16e3," [km]"/          &
                 "    tau_outer_ring = ",1es25.16e3//                 &
                 "(4) Particle system"/                               &
                 "    r_particle     = ",1es25.16e3," [m]"/           &
                 "    rho_particle   = ",1es25.16e3," [g/cm^3]"/      &
                 "    M_particle     = ",1es25.16e3," [g]"/           &
                 "    N_particle     = ",i12/                         &
                 "    Total mass     = ",1es25.16e3," [M_chariklo]"// &
                 "(5) Collision model"/                               &
                 "    Tkep           = ",1es25.16e3," [yr]"/          &
                 "    Tdur           = ",1es25.16e3," [yr]"/          &
                 "    e_CoR          = ",1es25.16e3/                  &
                 "    k_spg          = ",1es25.16e3/                  &
                 "    eta_dp         = ",1es25.16e3/                  &
                 "    0.1/eta_dp     = ",1es25.16e3," [yr]"//         &
                 "(6) Other useful quantities"/                       &
                 "    \Sigma_{ring}  = ",1es25.16e3," [g/cm^{2}]"/    &
                 "    \Sigma_{ptcl}  = ",1es25.16e3," [g/cm^{2}]"/    &
                 "    \sigma_{R}     = ",1es25.16e3," [cm/s]"/        &
                 "    \lambda_{crit} = ",1es25.16e3," [W_in_ring]"/   &
                 "========================================"/)

   end subroutine print_model_params

end module model_parameters

!===============================
!   MODULE: IO manager
!===============================
module io_manager
   use model_parameters
   implicit none
 
   !* Public parameters & variables
   !-(timing)
#ifdef TWO_BODY_COLLISION_TEST
   double precision, parameter :: time_end = 10.0d0*Tdur
   double precision, parameter :: dt = 1.0d-2*Tdur
   double precision, parameter :: dt_diag = 0.1d0*Tkep
   double precision, parameter :: dt_snap = 0.1d0*Tdur
#elif defined(TWO_BODY_COLLAPSE_TEST)
   double precision, parameter :: time_end = 10.0d0*Tkep
   double precision, parameter :: dt = 1.0d-2*Tdur
   double precision, parameter :: dt_diag = 0.1d0*Tkep
   double precision, parameter :: dt_snap = 0.1d0*Tkep
#elif defined(KEPLER_ROTATION_TEST)
   double precision, parameter :: time_end = 5.0d0*Tkep
   double precision, parameter :: dt = 1.0d-3*Tkep
   double precision, parameter :: dt_diag = 0.1d0*Tkep
   double precision, parameter :: dt_snap = 5.0d-2*Tkep
#else
   double precision, parameter :: time_end = 10.0d0*Tkep
   double precision, parameter :: dt = Tdur/16.0d0
   double precision, parameter :: dt_diag = 0.1d0*Tkep
   double precision, parameter :: dt_snap = Tkep
#endif
   double precision :: time_diag=0.0d0,time_snap=0.0d0,time_sys=0.0d0

   !-(I/O)
   character(len=16), parameter :: root_dir="result"
   character(len=16), parameter :: file_prefix_1st="snap"
   character(len=16), parameter :: file_prefix_2nd="proc"
   integer, parameter :: snap_num_rst=50 ! ID of the file used for the restart job
   integer :: snap_num=0


end module io_manager


!===============================
!   MODULE: User defined types
!===============================
module user_defined_types
   use, intrinsic :: iso_c_binding
   use fdps_vector
   use fdps_super_particle
   use model_parameters
   implicit none

   !**** Full particle type
   type, public, bind(c) :: full_particle !$fdps FP,EPI,EPJ,Force
      !$fdps copyFromForce full_particle (acc,acc)
      !$fdps copyFromFP full_particle (id,id) (mass,mass) (pos,pos) (vel,vel)
      !$fdps clear id=keep, mass=keep, pos=keep, vel=keep
      integer(kind=c_long_long) :: id
      real(kind=c_double)  mass !$fdps charge
      type(fdps_f64vec) :: pos !$fdps position
      type(fdps_f64vec) :: vel !$fdps velocity
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
      integer(c_long_long) :: idi
      real(c_double) :: r,r2,r_inv,r3_inv
      real(c_double) :: xij,vnij,f_ovlp,mj_eff
      real(c_double), dimension(3) :: ri,vi,ai,rij,vij
      integer(c_long_long), dimension(n_jp) :: idj
      real(c_double), dimension(n_jp) :: mj
      real(c_double), dimension(3,n_jp) :: rj,vj

#ifdef INTERACTION_FUNCTION_IS_EMPTY
!     Do nothing
#else
      !* Compute force
      do j=1,n_jp
         idj(j)  = ep_j(j)%id
         mj(j)   = ep_j(j)%mass
         rj(1,j) = ep_j(j)%pos%x
         rj(2,j) = ep_j(j)%pos%y
         rj(3,j) = ep_j(j)%pos%z
         vj(1,j) = ep_j(j)%vel%x
         vj(2,j) = ep_j(j)%vel%y
         vj(3,j) = ep_j(j)%vel%z
      end do
      do i=1,n_ip
         idi   = ep_i(i)%id
         ri(1) = ep_i(i)%pos%x
         ri(2) = ep_i(i)%pos%y
         ri(3) = ep_i(i)%pos%z
         vi(1) = ep_i(i)%vel%x
         vi(2) = ep_i(i)%vel%y
         vi(3) = ep_i(i)%vel%z
         ai = 0.0d0
         do j=1,n_jp
            if (idi == idj(j)) cycle
            rij(1) = ri(1) - rj(1,j)
            rij(2) = ri(2) - rj(2,j)
            rij(3) = ri(3) - rj(3,j)
            r2     = rij(1)*rij(1) &
                   + rij(2)*rij(2) &
                   + rij(3)*rij(3)
            r      = sqrt(r2)
            r_inv  = 1.0d0/r
            r3_inv = r_inv * r_inv * r_inv
            xij    = r - 2.0d0 * r_particle
#if !defined(KEPLER_ROTATION_TEST) && !defined(TWO_BODY_COLLAPSE_TEST)
            ! Spring-dashpot force
            ! (here, we assumed that mi=mj, ri,phys=rj,phys)
            if (xij < 0.0d0) then
               vij(1) = vi(1) - vj(1,j)
               vij(2) = vi(2) - vj(2,j)
               vij(3) = vi(3) - vj(3,j)
               vnij   = (vij(1)*rij(1)  &
                        +vij(2)*rij(2)  &
                        +vij(3)*rij(3)) &
                      * r_inv
               ai(1)  = ai(1) - 0.5d0 * (k_spg*xij + eta_dp*vnij) * r_inv * rij(1)
               ai(2)  = ai(2) - 0.5d0 * (k_spg*xij + eta_dp*vnij) * r_inv * rij(2)
               ai(3)  = ai(3) - 0.5d0 * (k_spg*xij + eta_dp*vnij) * r_inv * rij(3)
            end if
#endif
#if !defined(TWO_BODY_COLLISION_TEST) && !defined(KEPLER_ROTATION_TEST)
            ! Effective gravitational force
            f_ovlp = r/r_particle
            mj_eff = Ggrav * mj(j) * min(1.0d0, f_ovlp*f_ovlp*f_ovlp)
            ai(1)  = ai(1) - mj_eff * r3_inv * rij(1)
            ai(2)  = ai(2) - mj_eff * r3_inv * rij(2)
            ai(3)  = ai(3) - mj_eff * r3_inv * rij(3)
#endif
         end do
#if !defined(TWO_BODY_COLLISION_TEST) && !defined(TWO_BODY_COLLAPSE_TEST)
         ! Force from the Chariklo
         r3_inv = ri(1)*ri(1) + ri(2)*ri(2) + ri(3)*ri(3)
         r_inv  = 1.0d0/sqrt(r3_inv)
         r3_inv = Ggrav * M_chariklo * r_inv * r_inv * r_inv
         ai(1)  = ai(1) - r3_inv * ri(1)
         ai(2)  = ai(2) - r3_inv * ri(2)
         ai(3)  = ai(3) - r3_inv * ri(3)
#endif
         !* Substitute into FORCE type
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
      real(c_double) :: r3_inv,r_inv
      real(c_double), dimension(3) :: ri,ai,rij
      real(c_double), dimension(n_jp) :: mj
      real(c_double), dimension(3,n_jp) :: rj

#ifdef INTERACTION_FUNCTION_IS_EMPTY
      ! Do nothing
#else
#if !defined(TWO_BODY_COLLISION_TEST) && !defined(KEPLER_ROTATION_TEST)
      do j=1,n_jp
         mj(j)   = ep_j(j)%mass
         rj(1,j) = ep_j(j)%pos%x
         rj(2,j) = ep_j(j)%pos%y
         rj(3,j) = ep_j(j)%pos%z
      end do
      do i=1,n_ip
         ri(1) = ep_i(i)%pos%x
         ri(2) = ep_i(i)%pos%y
         ri(3) = ep_i(i)%pos%z
         ai = 0.0d0
         do j=1,n_jp
            rij(1) = ri(1) - rj(1,j)
            rij(2) = ri(2) - rj(2,j)
            rij(3) = ri(3) - rj(3,j)
            r3_inv = rij(1)*rij(1) &
                   + rij(2)*rij(2) &
                   + rij(3)*rij(3)
            r_inv  = 1.0d0/sqrt(r3_inv)
            r3_inv = r_inv * r_inv
            r_inv  = r_inv * mj(j)
            r3_inv = r3_inv * r_inv
            ai(1)  = ai(1) - r3_inv * rij(1)
            ai(2)  = ai(2) - r3_inv * rij(2)
            ai(3)  = ai(3) - r3_inv * rij(3)
         end do
         f(i)%acc%x = f(i)%acc%x + Ggrav*ai(1)
         f(i)%acc%y = f(i)%acc%y + Ggrav*ai(2)
         f(i)%acc%z = f(i)%acc%z + Ggrav*ai(3)
      end do
#endif
#endif

   end subroutine calc_gravity_psp

end module user_defined_types
