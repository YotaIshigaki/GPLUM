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


template <class T>
void output(T in);

template<>
void output(double in) {
    std::cout << in << std::endl;
}
template<>
void output(double *in) {
    std::cout << *in << std::endl;
}

int main(int argc, char *argv[]) {

    double x;
    x = 1.0;
    output(x);

    return 0;
}
