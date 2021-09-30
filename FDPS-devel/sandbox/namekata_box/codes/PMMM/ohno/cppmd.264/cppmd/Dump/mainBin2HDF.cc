#include <cerrno>
#include <cstdlib>
#include <mpi.h>
#include "Bin2HDF.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  Bin2HDF bin2hdf(argc,argv, MPI_COMM_WORLD);
  
  MPI_Finalize();
  return EXIT_SUCCESS;
}
