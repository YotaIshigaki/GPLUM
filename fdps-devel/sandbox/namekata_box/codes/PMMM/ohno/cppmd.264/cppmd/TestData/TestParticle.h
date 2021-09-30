#ifndef TESTPARTICLE_H
#define TESTPARTICLE_H
#include "CubicCell.h"

namespace TestParticle {
  // must be replaced data I/O or make initial data
  void makeparticle(ParticleArray& particlearray,
                    std::vector<int>& particle_setid,
                    std::vector<TypeRange>& typerangearray,
                    int total_num_particle, int total_num_set, 
                    int total_num_unit, int unit_id,
                    PostProcess& postprocess);
};

#endif
