#ifndef CLASS_DEVICE
#define CLASS_DEVICE

#include <pzc_vector.h>

class EpiDev{
public:
  float3 pos;
  int    id;
  int    walk;
};

class EpjDev{
public:
  float3 pos;
  float  charge;
  int    id;
};

class ForceDev{
public:
  float3 f;
  float  u;
};

#endif
