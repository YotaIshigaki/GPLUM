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
#include <typeinfo>

void func(int sizes[3]) {
   for (int i=0; i<3; i++) std::cout << sizes[i] << std::endl;
}

int main(int argc, char *argv[]) {

    int x,y,z;
    x = 10;
    y = 11;
    z = 5;

    int sz[] = {x, y, z};
    func(sz);

    return 0;
}
