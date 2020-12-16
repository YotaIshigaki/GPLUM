#include "Config.h"

#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <mpi.h>
#include <getopt.h>
#include <algorithm>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>

Config::Config(int argc, char *argv[], MPI_Comm comm):
    argv(argv), argc(argc), comm(comm), comm_rank(0), comm_size(1), args("")
{
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);
  ParseArguments(argc, argv);
  ValidateParams();
}

Config::~Config() {
}

std::ostream &operator<<(std::ostream &os, const Config &c) {
  if (c.comm_rank != 0) return os;
  os <<
      "# arguments:                   " << c.args << "\n"
      "# comm_size:                   " << c.comm_size << "\n"

      "# md.delta_t:                  " << c.md.delta_t << "\n"
      "# md.cutoff:                   " << c.md.cutoff << "\n"
      "# md.bsize:                    " << c.md.bsize << "\n"
      "# md.offset:                   " << c.md.offset << "\n"
      "# md.velocity:                 " << c.md.velocity << "\n"
      "# md.tmax:                     " << c.md.tmax << "\n"
      "# md.input_mode:               " << c.md.input_mode << "\n"
      "# md.citype:                   " << c.md.citype << "\n"
      "# md.ci_update_interval:       " << c.md.ci_update_interval << "\n"
      "# md.coulomb_type:             " << c.md.coulomb_type << "\n"
      "# md.num_lattice:              " << c.md.num_lattice << "\n"
      "# md.total_num_particle:       " << c.md.total_num_particle << "\n"
      "# md.withlong:                 " << c.md.withlong << "\n"
      "# md.withbond:                 " << c.md.withbond << "\n"
      "# md.verbose:                  " << c.md.verbose << "\n"
      "# md.particle_dump:            " << c.md.particle_dump << "\n"
      "# md.print_interval:           " << c.md.print_interval << "\n"
      "# md.reduce_interval:          " << c.md.reduce_interval << "\n"
      "# md.dump_without_water:       " << c.md.dump_without_water << "\n"
      "# md.calc_space:               " << c.md.calc_space << "\n"
      "# md.io_node_number:           " << c.md.io_node_number << "\n"
      "# md.move_comm_pattern:        " << c.md.move_comm_pattern << "\n"
      "# md.short_comm_pattern:       " << c.md.short_comm_pattern << "\n"
      "# md.celldiv:                  " << c.md.celldiv << "\n"
      "# md.celldiv3d_in_node:        " << c.md.celldiv3d_in_node << "\n"
      "# md.copynum3d:                " << c.md.copynum3d << "\n"
      "# md.nodediv3d:                " << c.md.nodediv3d << "\n"
      "# md.unitnode3d:               " << c.md.unitnode3d << "\n"
      "# md.withcorrecttrans:         " << c.md.withcorrecttrans << "\n"

      "# pme.alpha:                   " << c.pme.alpha << "\n"
      "# pme.grid_length:             " << c.pme.grid_length << "\n"
      "# pme.grid_number:             " << c.pme.grid_number << "\n"
      "# pme.order:                   " << c.pme.order << "\n"

      "# shake.tolerance:             " << c.shake.tolerance << "\n"
      "# shake.type:                  " << c.shake.type << "\n"
      "# shake.max_iterate:           " << c.shake.max_iterate << "\n"

      "# tempctrl.temperature:        " << c.tempctrl.temperature << "\n"
      "# tempctrl.tau_t:              " << c.tempctrl.tau_t << "\n"
      "# tempctrl.tau_p:              " << c.tempctrl.tau_p << "\n"
      "# tempctrl.method:             " << c.tempctrl.method << "\n"
      "# tempctrl.interval:           " << c.tempctrl.interval << "\n"

      "# amberfile.inpcrd:            " << c.amberfile.inpcrd << "\n"
      "# amberfile.prmtop:            " << c.amberfile.prmtop << "\n"
      "# amberfile.restrt:            " << c.amberfile.restrt << "\n"
      "# amberfile.mdcrd:             " << c.amberfile.mdcrd << "\n"
      "# amberfile.mdrst:             " << c.amberfile.mdrst << "\n"
      "# amberfile.mdcrd_interval:    " << c.amberfile.mdcrd_interval << "\n"
      "# amberfile.mdrst_interval:    " << c.amberfile.mdrst_interval << "\n"
  ;
  return os;
}

std::istream &operator>>(std::istream &is, Config &c) {
  return is;
}

template<typename T>
void Config::Read(T &t) {
  std::istringstream iss(::optarg);
  iss.exceptions(std::ios::badbit | std::ios::failbit);
  iss >> t;
}

void Config::NormalizeString(std::string &s) {
  // convert to lower case
  std::transform(s.begin(), s.end(), s.begin(),
                 static_cast<int(*)(int)>(std::tolower));
  // convert '_' to '-'
  std::replace(s.begin(), s.end(), '_', '-');
}

void Config::ReadCellIndexType(int &i) {
  std::string s(::optarg);
  NormalizeString(s);
  if (s == "0" || s == "fullcube")  { i = kFullCube; return; }
  if (s == "1" || s == "halfshell") { i = kHalfShell; return; }
  if (s == "2" || s == "fullshell") { i = kFullShell; return; }
  if (s == "3" || s == "smallball") { i = kSmallBall; return; }
  throw s;
}

void Config::ReadCoulombType(int &i) {
  std::string s(::optarg);
  NormalizeString(s);
  if (s == "0" || s == "ewald")  { i = kEwald; return; }
  if (s == "1" || s == "direct") { i = kDirect; return; }
  if (s == "2" || s == "zerodipole") { i = kZeroDipole; return; }
  throw s;
}

void Config::ReadInputMode(int &i) {
  std::string s(::optarg);
  NormalizeString(s);
  if (s == "0" || s == "amber")   { i = kAMBER; return; }
  //if (s == "1" || s == "pdb")     { i = kPDB; return; }
  //if (s == "2" || s == "restart") { i = kRESTART; return; }
  if (s == "3" || s == "water")   { i = kWATER; return; }
  if (s == "4" || s == "argon")   { i = kARGON; return; }
  if (s == "5" || s == "nacl")    { i = kNACL; return; }
  throw s;
}

void Config::ReadTempCtrlMethod(int &i) {
  std::string s(::optarg);
  NormalizeString(s);
  if (s == "0" || s == "none")            { i = kNONE; return; }
  if (s == "1" || s == "vscaling")        { i = kVEL_SCALING; return; }
  if (s == "2" || s == "nose-hoover")     { i = kNOSE_HOOVER; return; }
  if (s == "3" || s == "langevin")        { i = kLANGEVIN; return; }
  if (s == "4" || s == "andersen-hoover") { i = kANDERSEN_HOOVER; return; }
  throw s;
}

void Config::ParseArguments(int argc, char *argv[]) {
  using namespace std;

  // record arguments
  args += argv[0];
  for (char **p = argv + 1; *p; ++p) { args += " "; args += *p; }

  enum {
    LONG_OPTION_START = 300,

    HELP = LONG_OPTION_START,

    DELTA_T,
    CUTOFF,
    BSIZE,
    OFFSET,
    VELOCITY,
    TMAX,
    INPUT_MODE,
    CITYPE,
    CI_UPDATE_INTERVAL,
    COULOMB_TYPE,
    NUM_LATTICE,
    TOTAL_NUM_PARTICLE,
    WITHLONG,
    WITHBOND,
    VERBOSE,
    PARTICLE_DUMP,
    PRINT_INTERVAL,
    REDUCE_INTERVAL,
    DUMP_WITHOUT_WATER,
    CALC_SPACE,
    IO_NODE_NUMBER,
    REAL_NODE_NUMBER,
    MOVE_COMM_PATTERN,
    SHORT_COMM_PATTERN,
    CELLDIV,
    CELLDIV3D_IN_NODE,
    COPYNUM3D,
    NODEDIV3D,
    UNITNODE3D,
    WITHCORRECTTRANS,

    PME_ALPHA,
    PME_GRID_LENGTH,
    PME_GRID_NUMBER,
    PME_ORDER,

    SHAKE_TOLERANCE,
    SHAKE_TYPE,
    SHAKE_MAX_ITERATE,

    TEMPCTRL_TEMPERATURE,
    TEMPCTRL_TAU_T,
    TEMPCTRL_TAU_P,
    TEMPCTRL_METHOD,
    TEMPCTRL_INTERVAL,

    AMBERFILE_INPCRD,
    AMBERFILE_PRMTOP,
    AMBERFILE_RESTRT,
    AMBERFILE_MDCRD,
    AMBERFILE_MDRST,
    AMBERFILE_MDCRD_INTERVAL,
    AMBERFILE_MDRST_INTERVAL
  };

  for (;;) {
    static struct option long_opts[] = {
      {"help",                no_argument,       0, 'h'},

      {"delta-t",             required_argument, 0, DELTA_T},
      {"cutoff",              required_argument, 0, 'C'},
      {"bsize",               required_argument, 0, 's'},
      {"offset",              required_argument, 0, 'S'},
      {"velocity",            required_argument, 0, 'V'},
      {"tmax",                required_argument, 0, 'm'},
      {"input-mode",          required_argument, 0, INPUT_MODE},
      {"citype",              required_argument, 0, 'i'},
      {"ci-update-interval",  required_argument, 0, CI_UPDATE_INTERVAL},
      {"coulomb-type",        required_argument, 0, COULOMB_TYPE},
      {"num-lattice",         required_argument, 0, 'l'},
      {"total-num-particle",  required_argument, 0, 'n'},
      {"withlong",            no_argument,       0, 'L'},
      {"withbond",            no_argument,       0, 'B'},
      {"verbose",             required_argument, 0, 'v'},
      {"particle-dump",       required_argument, 0, 'd'},
      {"print-interval",      required_argument, 0, 'p'},
      {"reduce-interval",     required_argument, 0, 'r'},
      {"dump-without-water",  required_argument, 0, 'W'},
      {"calc-space",          required_argument, 0, CALC_SPACE},
      {"io-node-number",      required_argument, 0, IO_NODE_NUMBER},
      {"move-comm-pattern",   required_argument, 0, MOVE_COMM_PATTERN},
      {"short-comm-pattern",  required_argument, 0, SHORT_COMM_PATTERN},
      {"celldiv",             required_argument, 0, 'c'},
      {"celldiv3d-in-node",   required_argument, 0, CELLDIV3D_IN_NODE},
      {"copynum3d",           required_argument, 0, COPYNUM3D},
      {"nodediv3d",           required_argument, 0, NODEDIV3D},
      {"unitnode3d",          required_argument, 0, UNITNODE3D},
      {"withcorrecttrans",    no_argument,       0, WITHCORRECTTRANS},

      {"pme-alpha",               required_argument, 0, PME_ALPHA},
      {"pme-grid-length",         required_argument, 0, PME_GRID_LENGTH},
      {"pme-grid-number",         required_argument, 0, PME_GRID_NUMBER},
      {"pme-order",               required_argument, 0, PME_ORDER},

      {"shake-tolerance",     required_argument, 0, SHAKE_TOLERANCE},
      {"shake-type",          required_argument, 0, SHAKE_TYPE},
      {"shake-max-iterate",   required_argument, 0, SHAKE_MAX_ITERATE},

      {"temperature",         required_argument, 0, TEMPCTRL_TEMPERATURE},
      {"tau-t",               required_argument, 0, TEMPCTRL_TAU_T},
      {"tau-p",               required_argument, 0, TEMPCTRL_TAU_P},
      {"tempctrl-method",     required_argument, 0, TEMPCTRL_METHOD},
      {"tempctrl-interval",   required_argument, 0, TEMPCTRL_INTERVAL},

      {"inpcrd",              required_argument, 0, AMBERFILE_INPCRD},
      {"prmtop",              required_argument, 0, AMBERFILE_PRMTOP},
      {"restrt",              required_argument, 0, AMBERFILE_RESTRT},
      {"mdcrd",               required_argument, 0, AMBERFILE_MDCRD},
      {"mdrst",               required_argument, 0, AMBERFILE_MDRST},
      {"mdcrd-interval",      required_argument, 0, AMBERFILE_MDCRD_INTERVAL},
      {"mdrst-interval",      required_argument, 0, AMBERFILE_MDRST_INTERVAL},

      {0, 0, 0, 0}
    };
    int opt_index = 0;
    int opt = ::getopt_long(argc, argv,
                            ":hC:s:S:V:m:l:n:i:LBv:p:r:W:c:",
                            long_opts, &opt_index);
    if (opt == -1) break;
    try {
      switch (opt) {
        case DELTA_T:                    Read(md.delta_t); break;
        case 'C'/*CUTOFF*/:              Read(md.cutoff); break;
        case 's'/*BSIZE*/:               Read(md.bsize); break;
        case 'S'/*OFFSET*/:              Read(md.offset); break;
        case 'V'/*VELOCITY*/:            Read(md.velocity); break;
        case 'm'/*TMAX*/:                Read(md.tmax); break;
        case INPUT_MODE:                 ReadInputMode(md.input_mode); break;
        case 'i'/*CITYPE*/:              ReadCellIndexType(md.citype); break;
        case CI_UPDATE_INTERVAL:         Read(md.ci_update_interval); break;
        case COULOMB_TYPE:               ReadCoulombType(md.coulomb_type);
                                         break;
        case 'l'/*NUM_LATTICE*/:         Read(md.num_lattice); break;
        case 'n'/*TOTAL_NUM_PARTICLE*/:  Read(md.total_num_particle); break;
        case 'L'/*WITHLONG*/:            ++md.withlong; break;
        case 'B'/*WITHBOND*/:            ++md.withbond; break;
        case 'v'/*VERBOSE*/:             Read(md.verbose); break;
        case 'd'/*PARTICLE_DUMP*/:       Read(md.particle_dump); break;
        case 'p'/*PRINT_INTERVAL*/:      Read(md.print_interval); break;
        case 'r'/*REDUCE_INTERVAL*/:     Read(md.reduce_interval); break;
        case 'W'/*DUMP_WITHOUT_WATER*/:  Read(md.dump_without_water); break;
        case CALC_SPACE:                 Read(md.calc_space); break;
        case IO_NODE_NUMBER:             Read(md.io_node_number); break;
        case MOVE_COMM_PATTERN:          Read(md.move_comm_pattern); break;
        case SHORT_COMM_PATTERN:         Read(md.short_comm_pattern); break;
        case 'c'/*CELLDIV*/:             Read(md.celldiv); break;
        case CELLDIV3D_IN_NODE:          Read(md.celldiv3d_in_node); break;
        case COPYNUM3D:                  Read(md.copynum3d); break;
        case NODEDIV3D:                  Read(md.nodediv3d); break;
        case UNITNODE3D:                 Read(md.unitnode3d); break;
        case WITHCORRECTTRANS:           ++md.withcorrecttrans; break;

        case PME_ALPHA:                  Read(pme.alpha); break;
        case PME_GRID_LENGTH:            Read(pme.grid_length); break;
        case PME_GRID_NUMBER:            Read(pme.grid_number); break;
        case PME_ORDER:                  Read(pme.order); break;

        case SHAKE_TOLERANCE:            Read(shake.tolerance); break;
        case SHAKE_TYPE:                 Read(shake.type); break;
        case SHAKE_MAX_ITERATE:          Read(shake.max_iterate); break;

        case TEMPCTRL_TEMPERATURE:       Read(tempctrl.temperature); break;
        case TEMPCTRL_TAU_T:             Read(tempctrl.tau_t); break;
        case TEMPCTRL_TAU_P:             Read(tempctrl.tau_p); break;
        case TEMPCTRL_METHOD:            ReadTempCtrlMethod(tempctrl.method); break;
        case TEMPCTRL_INTERVAL:          Read(tempctrl.interval); break;

        case AMBERFILE_INPCRD:           Read(amberfile.inpcrd); break;
        case AMBERFILE_PRMTOP:           Read(amberfile.prmtop); break;
        case AMBERFILE_RESTRT:           Read(amberfile.restrt); break;
        case AMBERFILE_MDCRD:            Read(amberfile.mdcrd); break;
        case AMBERFILE_MDRST:            Read(amberfile.mdrst); break;
        case AMBERFILE_MDCRD_INTERVAL:   Read(amberfile.mdcrd_interval); break;
        case AMBERFILE_MDRST_INTERVAL:   Read(amberfile.mdrst_interval); break;

        case 'h'/*HELP*/: {
          if (comm_rank == 0) {
            cout <<
                "This is cppmd. A molecular dynamics core software.\n"
                "usage: cppmd [options]\n"
                "options:\n"
                "  -h  --help\n"
                "      --delta-t=<d>\n"
                "  -C, --cutoff=<d>\n"
                "  -s, --bsize=<d>\n"
                "  -S, --offset=<d>\n"
                "  -V, --velocity=<i>\n"
                "  -m, --tmax=<l>\n"
                "      --input-mode=<i>\n"
                "  -i, --citype=<i>\n"
                "      --ci-update-interval=<i>\n"
                "      --coulomb-type=<i>\n"
                "  -l, --num-lattice=<i>\n"
                "  -n, --total-num-particle=<i>\n"
                "  -L, --withlong=<i>\n"
                "  -B, --withbond=<i>\n"
                "  -v, --verbose=<i>\n"
                "  -d, --particle-dump=<i>\n"
                "  -p, --print-interval=<i>\n"
                "  -r, --reduce-interval=<i>\n"
                "  -W, --dump-without-water=<i>\n"
                "      --calc-space=<i>\n"
                "      --io-node-number=<i>\n"
                "      --move-comm-pattern=<i>\n"
                "      --short-comm-pattern=<i>\n"
                "  -c, --celldiv=<i>\n"
                "      --celldiv3d-in-node=<ixixi>\n"
                "      --copynum3d=<ixixi>\n"
                "      --nodediv3d=<ixixi>\n"
                "      --unitnode3d=<ixixi>\n"
                "      --withcorrecttrans\n"
                "\n"
                "      --pme-alpha=<d>\n"
                "      --pme-grid-length=<d>\n"
                "      --pme-grid-number=<i>\n"
                "      --pme-order=<i>\n"
                "\n"
                "      --shake-tolerance=<d>\n"
                "      --shake-type=<i>\n"
                "      --shake-max-iterate=<i>\n"
                "\n"
                "      --temperature=<d>\n"
                "      --tau-t=<d>\n"
                "      --tau-p=<d>\n"
                "      --tempctrl-method=<i>\n"
                "      --tempctrl-interval=<i>\n"
                "\n"
                "      --inpcrd=<s>\n"
                "      --prmtop=<s>\n"
                "      --restrt=<s>\n"
                "      --mdcrd=<s>\n"
                "      --mdrst=<s>\n"
                "      --mdcrd-interval=<i>\n"
                "      --mdrst-interval=<i>\n"
                << flush;
          }
          MPI_Finalize();
          exit(EXIT_SUCCESS);
        }
        case ':': {  // missing option argument
          if (comm_rank == 0) {
            cout << "cppmd: option requires an argument";
            if (0 < ::optopt && ::optopt < LONG_OPTION_START)
              cout << " -- '" << static_cast<char>(::optopt) << "'";
            else if (::optopt >= LONG_OPTION_START)
              cout << " -- '"
                  << long_opts[::optopt - LONG_OPTION_START].name << "'";
            cout << ".  try '-h' for help\n" << flush;
          }
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
        default: /* case '?': */ {  // unknown option
          if (comm_rank == 0) {
            cout << "cppmd: unknown option";
            if (0 < ::optopt && ::optopt < LONG_OPTION_START)
              cout << " -- '" << static_cast<char>(::optopt) << "'";
            else
              cout << " -- '" << &argv[::optind - 1][2] << "'";
            cout << ".  try '-h' for help\n" << flush;
          }
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
      }
    }
    catch (...) {  // invalid option argument
      if (comm_rank == 0) {
        cout << "cppmd: invalid option argument";
        if (0 < ::optopt && ::optopt < LONG_OPTION_START)
          cout << " -- '" << static_cast<char>(::optopt) << "'";
        else
          cout << " -- '" << argv[::optind - 1] << "'";
        cout << ".  try '-h' for help\n" << flush;
      }
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }
}

void Config::ValidateParams() {
}
