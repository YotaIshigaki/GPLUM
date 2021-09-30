#pragma once
#include <particle_simulator.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "param.h"
#include "mathfunc.h"
#include "kernel.h"
#include "class.h"
#include "force.h"
#include "prototype.h"

#ifdef ENABLE_GPU_CUDA
#include "force_gpu.h"
#endif
