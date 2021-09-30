#ifndef AMBERFILE_AMBERRESTRT_H
#define AMBERFILE_AMBERRESTRT_H

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif  // HAVE_CONFIG_H

#include <cstdio>
#include <string>
#include <vector>
#include "cppmd-int.h"

namespace amberfile {

const double  kVelocityUnit = 20.455;     // [picosecond]

template <typename T>
struct AmberRestrtT {
  typedef T               int_type;

  AmberRestrtT()
      : version(0000.000), date("00/00/00"), time("00:00:00"),
        irest(0), natom(0), current_time(0.0) {
  }

  float                   version;
  std::string             date;
  std::string             time;

  T                       irest;

  std::string             title;

  T                       natom;
  double                  current_time;   // [picosecond]

  std::vector<double>     coordinates;
  std::vector<double>     velocities;     // [angstroms per 1/20.455 ps]
  std::vector<double>     box;
};

// Define CPPMD_ENABLE_LARGEMODEL in config.h when you use 64-bit integer
// in Amber file I/O.
//
#ifdef CPPMD_ENABLE_LARGEMODEL
typedef AmberRestrtT<int64_t> AmberRestrt;
#else  // CPPMD_ENABLE_LARGEMODEL
typedef AmberRestrtT<int32_t> AmberRestrt;
#endif  // CPPMD_ENABLE_LARGEMODEL

int ReadFromFile(const std::string& file_name, AmberRestrt* rst, int irest);
int ReadFromFile(FILE* fp, AmberRestrt* rst, int irest);

int WriteToFile(const AmberRestrt& rst, const std::string& file_name,
                int ntb = 1);
int WriteToFile(const AmberRestrt& rst, FILE* fp, int ntb = 1);
FILE* PrepareTitledFile(const AmberRestrt& rst, const std::string& file_name);
void CloseTitledFile(FILE* fp);
int AppendCoordinatesToFile(const AmberRestrt& rst, FILE* fp,
                            int ntb = 1, int iwrap = 0);
int AppendVelocitiesToFile(const AmberRestrt& rst, FILE* fp, int ntb = 1);

}  // namespace amberfile

#endif  // AMBERFILE_AMBERRESTRT_H
