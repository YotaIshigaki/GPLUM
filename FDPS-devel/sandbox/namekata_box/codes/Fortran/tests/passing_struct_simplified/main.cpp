/* Standard headers */
#include <iostream>
#include <vector>

extern "C" void f_main_();

class FullParticle {
   public:
      double mass;
      double pos[3];
};
FullParticle pa[4];

int main(int argc, char *argv[]) {
   
   //* Initialize pa
   for (int i=0; i<4; i++) {
      pa[i].mass   = (double)(i+1);
      pa[i].pos[0] = (double)(i+1);
      pa[i].pos[1] = (double)(i+1);
      pa[i].pos[2] = (double)(i+1);
   }

   //* Call Fortran routine
   f_main_();

   return 0;
}

extern "C" {
void load_array_cpp_(void **p) {
   *p = (void *)&pa;
}
} // END of extern "C"
