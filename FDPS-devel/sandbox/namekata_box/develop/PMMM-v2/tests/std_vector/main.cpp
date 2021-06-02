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

    std::vector<int> list_a;
    const int imax=8;
    for (int i=0; i<imax; i++) list_a.push_back(i);

    std::vector<int> list_b;
    list_b = list_a;

    for (int i=0; i<list_b.size(); i++)
        std::cout << "i = " << i << ", list_b[i] = " << list_b[i] << std::endl;

    return 0;
}
