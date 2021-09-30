#pragma once
#include <particle_simulator.hpp>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "param.h"
#include "mathfunc.h"
struct kernel_t{
	kernel_t(){}
	static PS::F64 supportRadius(){	return 3.5;	}
};
//#include "kernel.h"
#include "EoS.h"
#include "class.h"
#include "prototype.h"
