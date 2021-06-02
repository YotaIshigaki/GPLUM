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
#include <sys/stat.h>
#include <time.h>
#include <typeinfo>
#include <complex>
#include <fftw3.h>
#include <complex.h>

int main(int argc, char *argv[]) {

#if 0
    std::cout << sizeof(std::complex<double>) << std::endl;
    std::cout << sizeof(fftw_complex) << std::endl;

    std::complex<double> z_std;
    fftw_complex z_fftw;

    z_std.real(1.0);
    z_std.imag(1.0);
    z_fftw[0] = z_std.real();
    z_fftw[1] = z_std.imag();

    double complex z_c99;
#else
    const double x = 1.0;
    const double y = 2.0;
    const std::complex<double> z_std_src(x, y);
    std::complex<double> z_std;
    z_std = z_std_src;
    std::cout << z_std << std::endl;
#endif

    return 0;
}
