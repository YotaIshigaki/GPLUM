#include <cerrno>
#include <cstdlib>
#include <mpi.h>
#include "HDF2RST.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  HDF2RST hdf2rst(argc,argv, MPI_COMM_WORLD);
  
  MPI_Finalize();
  return EXIT_SUCCESS;
}
