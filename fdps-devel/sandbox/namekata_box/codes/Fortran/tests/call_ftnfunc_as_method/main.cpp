/* Standard headers */
#include <iostream>
#include <vector>

extern "C" void f_main_();

//* Declaration
class FullParticle {
   public:
      double mass;
      double pos[3];
      double getCharge() const {
         return this->mass;
      }
      void clear();
};
//* Implementation
extern "C" void ftn_clear(FullParticle *ptcl);
void FullParticle::clear() {
   ftn_clear(this);
}

FullParticle ptcl;

int main(int argc, char *argv[]) {
   
   //* Call Fortran routine `ftn_clear`
   ptcl.clear();

   //* Check ptcl
   std::cout << ptcl.mass   << std::endl;
   std::cout << ptcl.pos[0] << std::endl;
   std::cout << ptcl.pos[1] << std::endl;
   std::cout << ptcl.pos[2] << std::endl;

   return 0;
}
