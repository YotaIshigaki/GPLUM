#include <math.h>
typedef struct cvec{
    double x,y,z;
}Cvec_Float64, *PCvec_Float64;

Cvec_Float64 Cvec_Float64_zero();
    
typedef struct cvec_f32{
    float x, y, z;
}Cvec_Float32, *PCvec_Float32;

typedef struct cmat_f32{
    float xx, yy, zz, xy, xz, yz;
}Cmat_Float32, *PCmat_Float32;

typedef struct cmat_f64{
    double xx, yy, zz, xsy, xz, yz;
}Cmat_Float64, *PCmat_Float64;

