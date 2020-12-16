#include "particle_simulator.hpp"

class FP{
public:
    PS::F64vec pos;
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::S64 id;
    PS::F64 mass;
    void dump(std::ostream & fout=std::cout) const {
//        fout<<"id="<<id<<std::endl;
        fout<<"pos= "<<pos<<std::endl;
    }
};

inline double frand()
{
    return (double) rand() / ((double) RAND_MAX + 1.);
//    return genrand_res53();
}
