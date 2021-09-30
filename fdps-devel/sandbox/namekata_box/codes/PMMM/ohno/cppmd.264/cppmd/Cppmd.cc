#include "Cppmd.h"

#include <cerrno>
#include <mpi.h>
#include <iostream>
#include <ostream>
#include "Config.h"

#include "simulator.h"

Cppmd::Cppmd(int argc, char *argv[], MPI_Comm comm):
    config_(new Config(argc, argv, comm))
{
  std::cout << *config_ << std::flush;
}

Cppmd::~Cppmd() {
  if (config_) delete config_;
}

void Cppmd::Run() {
  // TODO: integrate simulator class into Cppmd class
  cppmd::Simulator simulator(*config_);
  simulator.Run();
}
