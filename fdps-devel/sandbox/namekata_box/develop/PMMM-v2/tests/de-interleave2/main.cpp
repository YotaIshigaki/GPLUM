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
    _s = (_s | _s<<32) &  0x00000000ffffffff;  //0000 0000 0000 0000 0000 0000 0000 0000 1111 1111 1111 1111 1111 1111 1111 1111 
    _s = (_s | _s<<16) &  0x0000ffff0000ffff;  //0000 0000 0000 0000 1111 1111 1111 1111 0000 0000 0000 0000 1111 1111 1111 1111 
    _s = (_s | _s<<8) &  0x00ff00ff00ff00ff;  //0000 0000 1111 1111 0000 0000 1111 1111     
    _s = (_s | _s<<4) &   0x0f0f0f0f0f0f0f0f;  //0000 1111 0000 1111 0000 1111              
    _s = (_s | _s<<2) &   0x3333333333333333;  //00 11 00 11 00 11                          
    return (_s | _s<<1) & 0x5555555555555555;  //0101 0101                                  
}     

uint64_t shrinkBit(const uint64_t _s_in){
    uint64_t _s = _s_in;
    _s = (_s | (_s >>  1)) & 0x3333333333333333ULL;
    _s = (_s | (_s >>  2)) & 0x0F0F0F0F0F0F0F0FULL;
    _s = (_s | (_s >>  4)) & 0x00FF00FF00FF00FFULL;
    _s = (_s | (_s >>  8)) & 0x0000FFFF0000FFFFULL;
    _s = (_s | (_s >> 16)) & 0x00000000FFFFFFFFULL;
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
