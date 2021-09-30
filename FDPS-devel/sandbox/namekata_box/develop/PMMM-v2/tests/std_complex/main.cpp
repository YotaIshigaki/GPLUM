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

int main(int argc, char *argv[]) {

    typedef std::complex<double> cplx_base_t;
    std::cout << "sizeof(std::complex<double>) = " << sizeof(std::complex<double>) << std::endl;
    std::cout << "sizeof(std::complex<std::complex<double> >) = " << sizeof(std::complex<std::complex<double> >) << std::endl;

    std::complex<cplx_base_t> z;
    std::cout << "sizeof(std::complex<std::complex<double> >.real()) = " << sizeof(z.real()) << std::endl;
    std::cout << "sizeof(std::complex<std::complex<double> >.imag()) = " << sizeof(z.imag()) << std::endl;
     

    return 0;
}
