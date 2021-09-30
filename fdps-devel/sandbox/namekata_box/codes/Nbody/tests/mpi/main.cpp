/* Standard headers */
#include <cstdio>
#include <cstdlib>
#include <iostream>
/* FDPS headers */
#include <particle_simulator.hpp>


int main(int argc, char* argv[]) {

   PS::Initialize(argc,argv);


   PS::Finalize();
   return 0;

}
