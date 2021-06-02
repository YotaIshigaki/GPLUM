#include "particle_simulator.hpp"

class FP{
public:
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::S64 id;
    PS::F64 mass;
    void readAscii(FILE *fp) {
        fscanf(fp, "%lf%lf%lf",
               &pos[0], &pos[1], &pos[2]);
    }
    void dump(std::ostream & fout=std::cout) const {
        fout<<"pos= "<<pos<<std::endl;
    }
};

inline double frand()
{
    return (double) rand() / ((double) RAND_MAX + 1.);
}
