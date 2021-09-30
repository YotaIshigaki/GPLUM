// Include the standard C++ headers
#include <cmath>
#include <math.h>
#include <cfloat>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <time.h>
#include <complex>
// Include the header file of FDPS
#include <particle_simulator.hpp>

int main(int argc, char *argv[]) {
    using real_t = double;
    using cplx_t = std::complex<real_t>;
    const int p = 5;
    const PS::F64 charge = 1.0;
    const PS::F64vec center(0.0);

    PS::MultipoleMoment<real_t, cplx_t> mm;
    mm.alloc(p);
    PS::F64vec pos;

    // Case 1: pos = (1,0,0)
    pos.x = 1.0;
    pos.y = 0.0;
    pos.z = 0.0;
    mm.clear();
    mm.assign_particle(center, pos, charge);
    std::cout << "[Case1] pos = " << pos << std::endl;
    std::cout << "dipole:" << std::endl;
    std::cout << "   mm[1] = " << mm.buf[1] << std::endl;
    std::cout << "   mm[2] = " << mm.buf[2] << std::endl;
    std::cout << "   mm[3] = " << mm.buf[3] << std::endl;
    std::cout << "quadrupole:" << std::endl;
    std::cout << "   mm[4] = " << mm.buf[4] << std::endl;
    std::cout << "   mm[5] = " << mm.buf[5] << std::endl;
    std::cout << "   mm[6] = " << mm.buf[6] << std::endl;
    std::cout << "   mm[7] = " << mm.buf[7] << std::endl;
    std::cout << "   mm[8] = " << mm.buf[8] << std::endl;
    std::cout << std::endl;

    // Case 2: pos = (0,1,0)
    pos.x = 0.0;
    pos.y = 1.0;
    pos.z = 0.0;
    mm.clear();
    mm.assign_particle(center, pos, charge);
    std::cout << "[Case2] pos = " << pos << std::endl;
    std::cout << "dipole:" << std::endl;
    std::cout << "   mm[1] = " << mm.buf[1] << std::endl;
    std::cout << "   mm[2] = " << mm.buf[2] << std::endl;
    std::cout << "   mm[3] = " << mm.buf[3] << std::endl;
    std::cout << "quadrupole:" << std::endl;
    std::cout << "   mm[4] = " << mm.buf[4] << std::endl;
    std::cout << "   mm[5] = " << mm.buf[5] << std::endl;
    std::cout << "   mm[6] = " << mm.buf[6] << std::endl;
    std::cout << "   mm[7] = " << mm.buf[7] << std::endl;
    std::cout << "   mm[8] = " << mm.buf[8] << std::endl;
    std::cout << std::endl;

    // Case 3: pos = (0,0,1)
    pos.x = 0.0;
    pos.y = 0.0;
    pos.z = 1.0;
    mm.clear();
    mm.assign_particle(center, pos, charge);
    std::cout << "[Case2] pos = " << pos << std::endl;
    std::cout << "dipole:" << std::endl;
    std::cout << "   mm[1] = " << mm.buf[1] << std::endl;
    std::cout << "   mm[2] = " << mm.buf[2] << std::endl;
    std::cout << "   mm[3] = " << mm.buf[3] << std::endl;
    std::cout << "quadrupole:" << std::endl;
    std::cout << "   mm[4] = " << mm.buf[4] << std::endl;
    std::cout << "   mm[5] = " << mm.buf[5] << std::endl;
    std::cout << "   mm[6] = " << mm.buf[6] << std::endl;
    std::cout << "   mm[7] = " << mm.buf[7] << std::endl;
    std::cout << "   mm[8] = " << mm.buf[8] << std::endl;


    return 0;
}
