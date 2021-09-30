#include <iostream>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>

int main(){
   //* Initialize
   MPI::Init();
   int myrank = MPI::COMM_WORLD.Get_rank();
   int nprocs = MPI::COMM_WORLD.Get_size();

   //* MPI-IO
   MPI::File file;
   file.Open(MPI::COMM_WORLD,"test.dat",
             MPI::MODE_CREATE | MPI::MODE_WRONLY,
             MPI::INFO_NULL);
   if (myrank == 0) {
      MPI::Offset offset = {0};
      double x = 1.0e0;
      file.Write_at(offset,&x,1,MPI::DOUBLE);
   }
   file.Close();

   MPI::Finalize();
   return 0;
}
