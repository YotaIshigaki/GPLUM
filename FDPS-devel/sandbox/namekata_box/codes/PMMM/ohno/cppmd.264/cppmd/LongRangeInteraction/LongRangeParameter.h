#ifndef LONGRANGEPARAMETER_H
#define LONGRANGEPARAMETER_H

#include "Geometry.h"

typedef enum PMEType {
  PME,
  SmoothPME
} PMEType;

typedef enum MGType {
  SOR,
  VCycle,
  FMG,
} MGType;

//! How to Parallelize LongRange Calculation
/*!
  NoneLong : without LongRange Calculation
  Combine : Calculate Long on node that calculate Short
  Separate : Calculate Long only
  nD : parallelize n dimensionally
       0D : only one node
       1D : parallelize X-axis, gather atom info from YZ-Plate
       2D : parallelize XY-axis, gather atom info from Z-Tube
       3D : parallelize XYZ-axis, gather atom info form nearest neighbour
 */
typedef enum LongRangeMPIPlan {
  NoneLong,
  Combine0D,
  Combine1D,
  Combine2D,
  Combine3D,
  Separate0D,
  Separate1D,
  Separate2D,
  Separate3D,
} LongRangeMPIPlan;

struct LongRangeParameter {
  double cutoff;
  double alpha;
  double kCutoff;
  int surfaceDipole;
  SpaceVector<double> boxSize;
  SpaceVector<double> gridLengths;
  int order;
  PMEType pmeType;
  SpaceVector<int> grid_num;
  GeometryXYZ node_geometry;
  MGType multigridType;
  int multigridIteration;
  int multigridFFTlevel;
};

#endif
