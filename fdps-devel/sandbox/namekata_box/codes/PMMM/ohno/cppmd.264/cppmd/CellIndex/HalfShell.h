#ifndef HALFSHELL_H
#define HALFSHELL_H

#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <iostream>

#include <cassert>

#include "Common.h"
#include "ParticleInfo.h"
#include "CellIndex.h"

#include "Log.h"

size_t
halfshell_target(int celldiv, SpaceVector<int> ci, 
                 SpaceVector<double> cellsize, SpaceVector<double> margin,
                 double cutoff,
                 std::vector<int>& send_target_cell, 
                 std::vector<int>& recv_target_cell, 
                 ShiftCellArray& shift_cell);


#endif
