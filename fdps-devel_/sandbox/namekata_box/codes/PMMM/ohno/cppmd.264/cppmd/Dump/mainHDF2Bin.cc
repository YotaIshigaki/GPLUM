#include <cerrno>
#include <cstdlib>
#include <mpi.h>
#include "HDF2Bin.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  HDF2Bin hdf2bin(argc,argv, MPI_COMM_WORLD);
  
  MPI_Finalize();
  return EXIT_SUCCESS;
}
