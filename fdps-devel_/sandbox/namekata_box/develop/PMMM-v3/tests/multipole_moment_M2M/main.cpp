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
    const PS::F64vec pos(1.0, 0.0, 0.0); 

    PS::MultipoleMoment<real_t, cplx_t> mm_cell, mm_ptcl;
    mm_cell.alloc(p);
    mm_ptcl.alloc(p);

    mm_ptcl.clear();
    mm_ptcl.assign_particle(pos, pos, charge); // P2M
    mm_cell.clear();
    mm_cell.assign_from_MM(mm_ptcl, center, pos); // M2M
    std::cout << "Check mm_cell" << std::endl;
    std::cout << "dipole:" << std::endl;
    std::cout << "   mm[1] = " << mm_cell.buf[1] << std::endl;
    std::cout << "   mm[2] = " << mm_cell.buf[2] << std::endl;
    std::cout << "   mm[3] = " << mm_cell.buf[3] << std::endl;
    std::cout << "quadrupole:" << std::endl;
    std::cout << "   mm[4] = " << mm_cell.buf[4] << std::endl;
    std::cout << "   mm[5] = " << mm_cell.buf[5] << std::endl;
    std::cout << "   mm[6] = " << mm_cell.buf[6] << std::endl;
    std::cout << "   mm[7] = " << mm_cell.buf[7] << std::endl;
    std::cout << "   mm[8] = " << mm_cell.buf[8] << std::endl;
    std::cout << std::endl;

    return 0;
}
