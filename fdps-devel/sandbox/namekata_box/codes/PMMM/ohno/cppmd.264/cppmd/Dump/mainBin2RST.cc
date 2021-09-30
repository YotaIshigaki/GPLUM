#include <cerrno>
#include <cstdlib>
#include <mpi.h>
#include "Bin2RST.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  Bin2RST bin2rst(argc,argv, MPI_COMM_WORLD);
  
  MPI_Finalize();
  return EXIT_SUCCESS;
}
