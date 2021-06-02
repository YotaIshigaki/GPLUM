#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>

typedef class full_particle {
public:
    long long id;
    double mass;

    void readBinary(FILE *fp) {
        //fread(this, sizeof(class full_particle), 1, fp);
        //fread(this, sizeof(*this), 1, fp);
        fread(&this->id, sizeof(long long), 1, fp);
        fread(&this->mass, sizeof(double), 1, fp);
    }
} Full_particle;

int main(int argc, char *argv[]) {

    // Set particle data
    Full_particle ptcl_out;
    ptcl_out.id = 7;
    ptcl_out.mass = 7.7;
    std::cout << "sizeof(class full_particle) = " << sizeof(class full_particle) << std::endl;

    // Output particle data
    FILE *fp;
    char filename[64];
    strcpy(filename,"data.dat");
    if ((fp = fopen(filename,"wb")) == NULL) {
        std::cerr << "cannot open file " << filename << std::endl; 
        std::exit(-1);
    }
    fwrite(&ptcl_out, sizeof(Full_particle), 1, fp);
    fclose(fp);

    // Read particle data
    Full_particle ptcl_in;
    if ((fp = fopen(filename,"rb")) == NULL) {
        std::cerr << "cannot open file " << filename << std::endl; 
        std::exit(-1);
    }
    ptcl_in.readBinary(fp);
    long pos = std::ftell(fp);
    std::cout << "id   = " << ptcl_in.id << std::endl;
    std::cout << "mass = " << ptcl_in.mass << std::endl;
    std::cout << "pos  = " << pos  << std::endl;
    fclose(fp);
    
    return 0;
}
