#ifndef CUBICCELL_H
#define CUBICCELL_H

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
cubic_target(int celldiv, SpaceVector<int> ci, int reach,
             std::vector<int>& send_target_cell, 
             std::vector<int>& recv_target_cell, 
             ShiftCellArray& shift_cell);

size_t
cubic_target(SpaceVector<int> celldiv, SpaceVector<int> ci, int reach,
             std::vector<int>& send_target_cell, 
             std::vector<int>& recv_target_cell, 
             ShiftCellArray& shift_cell);
#endif
