#ifndef CONFIG_H_
#define CONFIG_H_ 1

#include <mpi.h>
#include <iosfwd>
#include <string>
#include "Vec.h"

class Config {
 public:
  enum CellMethod {
    kFullCube, kHalfShell, kFullShell, kSmallBall
  };
  enum CoulombType {
    kEwald, kDirect, kZeroDipole
  };
  enum InputMode {
    kAMBER, kPDB, kRESTART, kWATER, kARGON, kNACL
  };
  enum TempCtrlMethod {
    kNONE, kVEL_SCALING, kNOSE_HOOVER, kLANGEVIN, kANDERSEN_HOOVER
  };

 public:
  Config(int argc, char *argv[], MPI_Comm comm);
  ~Config();

  friend std::ostream &operator<<(std::ostream &os, const Config &c);
  friend std::istream &operator>>(std::istream &is, Config &c);

 private:
  template<typename T> void Read(T &t);

  void NormalizeString(std::string &s);
  void ReadCellIndexType(int &i);
  void ReadCoulombType(int &i);
  void ReadInputMode(int &i);
  void ReadTempCtrlMethod(int &i);

  void ParseArguments(int argc, char *argv[]);
  void ValidateParams();

 public:
  char **argv;
  int argc;
  MPI_Comm comm;
  int comm_rank;
  int comm_size;
  std::string args;

  // md control
  struct Md {
    Md():
        delta_t(.5), cutoff(.0), bsize(.0), offset(.0), velocity(.0), tmax(0),
        input_mode(/*?*/0), citype(1), ci_update_interval(40),
        coulomb_type(/*?*/0), num_lattice(1), total_num_particle(16),
        withlong(0), withbond(0), verbose(0), particle_dump(0), print_interval(0),
        reduce_interval(0), dump_without_water(0), calc_space(0),
        io_node_number(16), move_comm_pattern(1), short_comm_pattern(1),
        celldiv(4), celldiv3d_in_node(0), copynum3d(1), nodediv3d(0),
	  unitnode3d(0), withcorrecttrans(0) {}
    double delta_t;
    double cutoff;
    double bsize;
    double offset;
    double velocity;
    long tmax;
    int input_mode;
    int citype;
    int ci_update_interval;
    int coulomb_type;
    int num_lattice;
    int total_num_particle;
    int withlong;
    int withbond;
    int verbose;
    int particle_dump;
    int print_interval;
    int reduce_interval;
    int dump_without_water;
    int calc_space;
    int io_node_number;
    int move_comm_pattern;
    int short_comm_pattern;
    int celldiv;
    V3i celldiv3d_in_node;
    V3i copynum3d;
    V3i nodediv3d;
    V3i unitnode3d;
    int withcorrecttrans;
  } md;

  // pme
  struct Pme {
    Pme():
    alpha(0.0), grid_length(0.0), grid_number(0), order(4) {}
    double alpha;
    double grid_length;
    int grid_number;
    int order;
  } pme;

  // shake
  struct Shake {
    Shake():
        tolerance(1.0e-10), type(2), max_iterate(3000) {}
    double tolerance;
    int type;
    int max_iterate;
  } shake;

  // temperature
  struct TempCtrl {
    TempCtrl():
        temperature(300.0), tau_t(1.0), tau_p(1.0), method(0), interval(50) {}
    double temperature;
    double tau_t;
    double tau_p;
    int method;
    int interval;
  } tempctrl;

  // amber file
  struct AmberFile {
    AmberFile():
        mdcrd_interval(0), mdrst_interval(0) {}
    std::string inpcrd;
    std::string prmtop;
    std::string restrt;
    std::string mdcrd;
    std::string mdrst;
    int mdcrd_interval;
    int mdrst_interval;
  } amberfile;
};

#endif  // CONFIG_H_
