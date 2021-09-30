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
#include <unordered_map>
#include <random>
#include <sys/stat.h>
#include <time.h>
// Include the header file of FDPS
#include <particle_simulator.hpp>

uint64_t separateBit(const uint64_t _s_in){
    uint64_t _s = _s_in;
    _s = (_s | _s<<32) & 0xffff00000000ffff; //11111111 11111111 00000000 00000000 00000000 00000000 11111111 11111111
    _s = (_s | _s<<16) & 0x00ff0000ff0000ff; //00000000 11111111 00000000 00000000 11111111 00000000 00000000 11111111
    _s = (_s | _s<<8) & 0xf00f00f00f00f00f;  //1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111
    _s = (_s | _s<<4) & 0x30c30c30c30c30c3;  //11 00 00 11 00 11 00 00 11
    return (_s | _s<<2) & 0x9249249249249249;  //1 0 0 1 0 0 1 0 0 1 0 0 1
}

uint64_t shrinkBit_failed(const uint64_t _s_in){
    uint64_t _s = _s_in;
    _s = (_s | _s>>2);
    std::cout << "step 1  = " << PS::GetBinString(_s) << std::endl;
    _s = _s & 0x9249249249249249;  //1 0 0 1 0 0 1 0 0 1 0 0 1
    std::cout << "step 2  = " << PS::GetBinString(_s) << std::endl;
    _s = (_s | _s>>4);
    std::cout << "step 3  = " << PS::GetBinString(_s) << std::endl;
    _s = _s & 0x30c30c30c30c30c3;  //11 00 00 11 00 11 00 00 11
    std::cout << "step 4  = " << PS::GetBinString(_s) << std::endl;
    _s = (_s | _s>>8);
    std::cout << "step 5  = " << PS::GetBinString(_s) << std::endl;
    _s = _s  & 0xf00f00f00f00f00f;  //1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111 0000 0000 1111
    std::cout << "step 6  = " << PS::GetBinString(_s) << std::endl;
    _s = (_s | _s>>16);
    std::cout << "step 7  = " << PS::GetBinString(_s) << std::endl;
    _s = _s & 0x00ff0000ff0000ff; //00000000 11111111 00000000 00000000 11111111 00000000 00000000 11111111
    std::cout << "step 8  = " << PS::GetBinString(_s) << std::endl;
    _s = (_s | _s>>32);
    std::cout << "step 9  = " << PS::GetBinString(_s) << std::endl;
    _s = _s & 0xffff00000000ffff; //11111111 11111111 00000000 00000000 00000000 00000000 11111111 11111111
    std::cout << "step 10 = " << PS::GetBinString(_s) << std::endl;
    return _s;
}

uint64_t shrinkBit(const uint64_t _s_in){
    uint64_t _s = _s_in;
    _s = (_s | _s>>2)  & 0x30c30c30c30c30c3; // Repeat of pattern 00 00 11
    _s = (_s | _s>>4)  & 0xf00f00f00f00f00f; // Repeat of pattern 0000 0000 1111
    _s = (_s | _s>>8)  & 0x00ff0000ff0000ff; // Repeat of pattern 00000000 00000000 11111111
    _s = (_s | _s>>16) & 0xffff00000000ffff; // Repeat of pattern 0^{16} 0^{16} 1^{16},
                                             // where x^{n} reprensets a state in which x is continuously arranged in n times.
    _s = (_s | _s>>32) & 0x00000000ffffffff; // Repaat of pattern 0^{32} 0^{32} 1^{32}
    return _s;
}

int main(int argc, char *argv[]) {

    uint64_t x = 0x1fffffULL; 

    std::cout << "x      = " << PS::GetBinString(x) << std::endl;
    uint64_t xsep = separateBit(x);
    std::cout << "xsep   = " << PS::GetBinString(xsep) << std::endl;
    uint64_t xshrnk = shrinkBit(xsep);
    std::cout << "xshrnk = " << PS::GetBinString(xshrnk) << std::endl;
    std::cout << "x      = " << PS::GetBinString(x) << std::endl;

    return 0;
}
