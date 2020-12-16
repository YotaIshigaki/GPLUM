/* Standard headers */
#include <iostream>

enum fdps_boundary_condition {
   bc_open,
   bc_periodic
};
extern "C" void f_main_(enum fdps_boundary_condition bc);

int main(int argc, char *argv[]) {
   
   //* Call Fortran routine
   enum fdps_boundary_condition bc;
   bc = bc_periodic;
   f_main_(bc);

   return 0;
}

extern "C" {

void send_enum(enum fdps_boundary_condition bc){
   if (bc == bc_open) {
      std::cout << "bc_open in C++" << std::endl;
   } else if (bc == bc_periodic) {
      std::cout << "bc_periodic in C++" << std::endl;
   }
}

} // END of extern "C"
