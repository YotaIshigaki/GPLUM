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
#include <array>
#include <sys/stat.h>
#include <time.h>
#include <complex>

namespace mathematical_constants {

constexpr double pi = std::atan(1.0) * 4.0; // circular constant

}
namespace math_const = mathematical_constants;

template <int n>
class TestClass {
public:
    double x[n];
};

int main(int argc, char *argv[]) {

    constexpr int size = 8;
    TestClass<size> test;

    return 0;
}
