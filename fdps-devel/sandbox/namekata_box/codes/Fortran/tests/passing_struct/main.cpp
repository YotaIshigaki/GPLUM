/* Standard headers */
#include <iostream>
#include <vector>

extern "C" void f_main_();

class FullParticle {
   public:
      double mass;
      double pos[3];
      double vel[3];
};
FullParticle ps;
FullParticle pa[4];
std::vector<FullParticle> FPvec(4);

int main(int argc, char *argv[]) {
   
   //* Initialize ps
   ps.mass   = 2.4e0;
   ps.pos[0] = 1.0e0;
   ps.pos[1] = 2.0e0;
   ps.pos[2] = 3.0e0;
   ps.vel[0] = 4.0e0;
   ps.vel[1] = 5.0e0;
   ps.vel[2] = 6.0e0;

   //* Initialize pa
   for (int i=0; i<4; i++) {
      pa[i].mass   = (double)(i+1);
      pa[i].pos[0] = (double)(i+1);
      pa[i].pos[1] = (double)(i+1);
      pa[i].pos[2] = (double)(i+1);
      pa[i].vel[0] = (double)(i+1);
      pa[i].vel[1] = (double)(i+1);
      pa[i].vel[2] = (double)(i+1);
   }

   // Initialize FPvec
   for (int i=0; i<4; i++) {
      FPvec[i].mass   = (double)(i+1);
      FPvec[i].pos[0] = (double)(i+1);
      FPvec[i].pos[1] = (double)(i+1);
      FPvec[i].pos[2] = (double)(i+1);
      FPvec[i].vel[0] = (double)(i+1);
      FPvec[i].vel[1] = (double)(i+1);
      FPvec[i].vel[2] = (double)(i+1);
   }

   //* Call Fortran routine
   f_main_();

   return 0;
}

extern "C" {

void load_scalar_cpp_(void **p) {
   //*p = (void *)&ps; 
   *p = (void *)&FPvec[0]; 
}

void load_array_cpp_(void **p) {

   //*p = (void *)&pa;
   *p = (void *)&FPvec[0];
   
}


} // END of extern "C"
