#if PARALLELIZATION_TYPE == 0
#include "main-serial.cpp"
#elif PARALLELIZATION_TYPE == 1
#include "main-omp.cpp"
#elif PARALLELIZATION_TYPE == 2
#include "main-flat_mpi.cpp"
#elif PARALLELIZATION_TYPE == 3
#include "main-hybrid.cpp"
#else
#error macro ``PARALLELIZATION_TYPE'' has an invalid value.
#endif
