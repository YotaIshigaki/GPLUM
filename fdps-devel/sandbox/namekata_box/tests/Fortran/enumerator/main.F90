module fdps_module
   use, intrinsic :: iso_c_binding
   implicit none 

   !***** PS::INTERACTION_LIST_MODE
  !enum, bind(c) :: fdps_interaction_list_mode
  !   enumerator :: make_list
  !   enumerator :: make_list_for_reuse
  !   enumerator :: reuse_list
  !end enum fdps_interaction_list_mode
   enum, bind(c) 
      enumerator :: make_list
      enumerator :: make_list_for_reuse
      enumerator :: reuse_list
   end enum 

   !* Public member functions 
   public :: test

   contains

   !--------------------------------------------------------------------
   subroutine test()
     !implicit none
     !type(fdps_interaction_list_mode), intent(in) :: val  
     !write(*,*)'val = ',val
   end subroutine test

end module fdps_module


program main
   use fdps_module
   implicit none

  !call test(make_list)

end program main
