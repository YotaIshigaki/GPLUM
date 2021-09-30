#include <iostream>
#include "mpi.h"

int main(int argc, char *argv[]) {
   //* Local variables
   int nprocs,myrank;

   MPI::Init(argc,argv);
   nprocs = MPI::COMM_WORLD.Get_size();
   myrank = MPI::COMM_WORLD.Get_rank();

   std::cout << "Hello from " << myrank << std::endl;

   MPI::Finalize();

   return 0;
}
