#ifndef MPIFFT3D_H
#define MPIFFT3D_H

#include "UseMPI.h"
#include "FFT3D.h"

namespace PMEModule {
struct Communicator {
#ifdef USE_MPI
  MPI_Comm communicator;
#endif
};
}

#endif
