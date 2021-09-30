#ifndef BIN2RST_H
#define BIN2RST_H
#include <mpi.h>
#include <string>

class Bin2RST {
public:
  char **argv;
  int argc;
  MPI_Comm comm;
  int comm_rank;
  int comm_size;
  std::string args;

  std::string filenamebase;
  int number_of_restorefiles;
  std::string crdname;
  std::string rstname;

  Bin2RST(int argc, char* argv[], MPI_Comm comm);

  template<typename T>
    void Read(T &t);

  int ParseArguments(int argc, char* argv[]);

};
#endif
