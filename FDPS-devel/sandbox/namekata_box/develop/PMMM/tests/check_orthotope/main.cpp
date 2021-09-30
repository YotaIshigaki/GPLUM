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


int main(int argc, char *argv[]) {

    PS::F64ort target_box, tc_vertex;

    target_box = PS::F64ort(PS::F64vec(0.0, 0.0, 0.0), 
                            PS::F64vec(1.0, 1.0, 1.0));
    tc_vertex  = PS::F64ort(PS::F64vec(0.2, 0.2, 0.2),
                            PS::F64vec(0.8, 0.8, 0.8));
    std::cout << target_box.contained(tc_vertex) << std::endl; // (A)
    std::cout << tc_vertex.contained(target_box) << std::endl; // (B)

    std::cout << target_box.isContainedBy(tc_vertex) << std::endl; // (C)
    std::cout << tc_vertex.isContainedBy(target_box) << std::endl; // (D)

    return 0;
}
