#include <iostream>
#include <fstream>

typedef struct {
    double time;
    int nbodies;
    int ndim;
    int nsph;
    int ndark;
    int nstar;
} MAGI_Tipsy_Header;

typedef struct {
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    int idx;
} MAGI_Tipsy_Particle;

int main(int argc, char *argv[]) {

    MAGI_Tipsy_Header header;
    header.time = 5.0;
    header.nbodies = 4096;
    header.ndim = 3;
    header.nsph = 1024;
    header.ndark = 1024;
    header.nstar = 2048;

    MAGI_Tipsy_Particle ptcl;
    ptcl.mass = 1.3;
    ptcl.pos[0] = 1.0;
    ptcl.pos[1] = 2.0;
    ptcl.pos[2] = 3.0;
    ptcl.vel[0] = - 1.0;
    ptcl.vel[1] = - 2.0;
    ptcl.vel[2] = - 3.0;
    ptcl.eps    = 0.01;
    ptcl.idx    = 1;
 
    const std::string file_name = "data.dat"; 
    std::ofstream output_file;
    output_file.open(file_name.c_str(), std::ios::trunc | std::ios::binary); 
    output_file.write((char *) &header, sizeof(header));
    output_file.write((char *) &ptcl, sizeof(ptcl));
    output_file.close();

    return 0;
}
