#include <cerrno>
#include <cstdlib>
#include <mpi.h>
#include "Cppmd.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  {
    Cppmd core(argc, argv, MPI_COMM_WORLD);
    core.Run();
  }
  MPI_Finalize();
  return EXIT_SUCCESS;
}
