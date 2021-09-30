#include <stdio.h>
#include <stdlib.h>
#include <particle_simulator.hpp>

#include "treepm.hpp"
#include "run_param.hpp"
#include "mpi.h"

#include "prototype.h"

PS::F64 get_dtime(PS::ParticleSystem<FPtreepm> &ptcl, 
		  run_param &this_run) 
{
  return ptcl[0].pos[0];
}
