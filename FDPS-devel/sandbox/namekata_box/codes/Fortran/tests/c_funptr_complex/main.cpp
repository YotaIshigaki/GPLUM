//* Full Particle class
class full_particle {
   public:
      double mass;
      double eps;
      double pos[3];
      double vel[3];
      double pot;
      double acc[3];
};
//* Essential Particle I class
class essential_particle_i {
   public:
      double mass;
      double eps;
      double pos[3];
};
//* Essential Particle J class
class essential_particle_j {
   public:
      double mass;
      double pos[3];
};
//* Force class
class force {
   public:
      double pot;
      double acc[3];
};

//* Fortran main routine
extern "C" void f_main_();

//* Functions called from Fortran program
extern "C" {

void calc_force_all(void (*func)(struct essential_particle_i *,
                             int *n_ip, 
                             struct essential_particle_j *,
                             int *n_jp,
                             struct force *),
                struct essential_particle_i *ep_i, 
                int *n_ip,
                struct essential_particle_j *ep_j,
                int *n_jp,
                struct force *f) {
    func(ep_i,n_ip,ep_j,n_jp,f);
}

} // END of extern "C"

int main(int argc, char *argv[]) {

   //* Call Fortran main routine
   f_main_();

}
