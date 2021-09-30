#include <mpi.h>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include "Bin2RST.h"
#include "Converter.h"

Bin2RST::Bin2RST(int argc, char* argv[], MPI_Comm comm):
  argv(argv), argc(argc), comm(comm), comm_rank(0), comm_size(1), args(""),
  filenamebase("binaryrestorefile"), number_of_restorefiles(1),
  crdname(""), rstname("")
{
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);
  ParseArguments(argc, argv);
  Converter<CombinedParticleArray> converter(filenamebase,number_of_restorefiles,comm);
  if(!crdname.empty()){
    converter.dump_amber_crd(crdname);
  }
  if(!rstname.empty()){
    converter.dump_amber_rst(rstname);
  }
}

template<typename T>
void Bin2RST::Read(T &t) {
  std::istringstream iss(::optarg);
  iss.exceptions(std::ios::badbit | std::ios::failbit);
  iss >> t;
}

int Bin2RST::ParseArguments(int argc, char* argv[])
{
  // record arguments
  args += argv[0];
  for (char **p = argv + 1; *p; ++p) { args += " "; args += *p; }

  enum {
    LONG_OPTION_START = 300,

    HELP = LONG_OPTION_START,

    FILENAME,
    NUMBER,
    AMBERFILE_MDCRD,
    AMBERFILE_MDRST
  };

  for (;;) {
    static struct option long_opts[] = {
      {"help",     no_argument,       0, 'h'},

      {"filename", required_argument, 0, 'f'},
      {"number",   required_argument, 0, 'n'},
      {"mdcrd",    required_argument, 0, AMBERFILE_MDCRD},
      {"mdrst",    required_argument, 0, AMBERFILE_MDRST},

      {0, 0, 0, 0}
    };
    int opt_index = 0;
    int opt = ::getopt_long(argc, argv,
			       ":hf:n:",
			       long_opts, &opt_index);
    if (opt == -1) break;
    try{
      switch (opt) {
      case 'f': Read(filenamebase); break;
      case 'n': Read(number_of_restorefiles); break;
      case AMBERFILE_MDCRD: Read(crdname); break;
      case AMBERFILE_MDRST: Read(rstname); break;
      case 'h'/*HELP*/: {
	{
	  std::cout <<
	    "options:\n"
	    "  -h,  --help\n"
	    "  -f,  --filename<s>\n"
	    "  -n,  --number<s>\n"
	    "       --mdcrd=<s>\n"
	    "       --mdrst=<s>\n"
		    << std::flush;
	}
	MPI_Finalize();
	exit(EXIT_SUCCESS);
      }
      case ':':  {  // missing option argument
	std::cout << "option requires an argument";
	if (0 < ::optopt && ::optopt < LONG_OPTION_START)
	  std::cout << " -- '" << static_cast<char>(::optopt) << "'";
	else if (::optopt >= LONG_OPTION_START)
	  std::cout << " -- '"
	       << long_opts[::optopt - LONG_OPTION_START].name << "'";
	std::cout << ".  try '-h' for help\n" << std::flush;
	MPI_Finalize();
	exit(EXIT_FAILURE);
      }
      default:  /* case '?': */ {  // unknown option
	std::cout << "cppmd: unknown option";
	if (0 < ::optopt && ::optopt < LONG_OPTION_START)
	  std::cout << " -- '" << static_cast<char>(::optopt) << "'";
	else
	  std::cout << " -- '" << &argv[::optind - 1][2] << "'";
	std::cout << ".  try '-h' for help\n" << std::flush;
	MPI_Finalize();
	exit(EXIT_FAILURE);
      }
      }
    }
    catch (...) {  // invalid option argument
      {
	std::cout << "invalid option argument";
        if (0 < ::optopt && ::optopt < LONG_OPTION_START)
	  std::cout << " -- '" << static_cast<char>(::optopt) << "'";
        else
	  std::cout << " -- '" << argv[::optind - 1] << "'";
	std::cout << ".  try '-h' for help\n" << std::flush;
      }
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }
}
