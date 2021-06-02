!===============
! FDPS enum
!===============
module fdps_enum_types

   !**** PS::BOUNDARY_CONDITION
   enum, bind(c)
      enumerator :: bc_open
      enumerator :: bc_periodic
   end enum 
  !enum, bind(c) :: fdps_boundary_condition
  !   enumerator :: bc_open
  !   enumerator :: bc_periodic
  !end enum fdps_boundary_condition

end module fdps_enum_types

!-----------------------------------------------------------------------
!/////////////////////// M A I N   R O U T I N E ///////////////////////
!-----------------------------------------------------------------------
subroutine f_main(bc)
   use, intrinsic :: iso_c_binding
   use fdps_enum_types
   implicit none
   !* Interface
   interface
      subroutine send_enum(bc) bind(c)
         use, intrinsic :: iso_c_binding
         implicit none
         integer(kind=c_int), value, intent(in) :: bc
      end subroutine send_enum
   end interface
   integer(kind=c_int), value, intent(in) :: bc

   !* Check received enum value
   if (bc == bc_open) then
      write(*,*)'bc_open'
   else if (bc == bc_periodic) then
      write(*,*)'bc_periodic'
   end if

   !* Send enum value
   call send_enum(bc_periodic)
   
   stop
end subroutine f_main

