// Include the standard C++ headers
#include <cmath>
#include <math.h>
#include <cfloat>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <sys/stat.h>
#include <time.h>
// Include the header file of FDPS
#include <particle_simulator.hpp>

void initialize(const PS::S32 pfmm=5,
                const PS::S32vec n_cell=(8,8,8)) {
    std::cout << "pfmm = " << pfmm << std::endl;
    std::cout << "n_cell = " << n_cell << std::endl;
}

int main(int argc, char *argv[]) {

    initialize();
    return 0;
}
