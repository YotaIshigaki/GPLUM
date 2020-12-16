#ifndef CBINTERFACEFWD_H
#define CBINTERFACEFWD_H

namespace CBModule {
struct CBContext;

class CBInterface;
namespace CBInterface_ {
typedef struct {
  int array;
  AtomID index;
} ParticleIndex;
  typedef Particle* ParticleLocation;
  typedef Force* ForceLocation;
  typedef CBContext* Context;
};
}

#endif
