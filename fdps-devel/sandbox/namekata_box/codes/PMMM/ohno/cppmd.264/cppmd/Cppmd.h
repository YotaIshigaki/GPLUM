#ifndef CPPMD_H_
#define CPPMD_H_ 1

#include <mpi.h>
#include <iosfwd>
#include <string>

class Config;

class Cppmd {
 public:
  Cppmd(int argc, char *argv[], MPI_Comm comm);
  ~Cppmd();

  void Run();

 private:
  Config *config_;
};

#endif  // CPPMD_H_
