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

struct Value {
   int a;
   int b;
};


int main(int argc, char *argv[]) {

    std::unordered_map<int, Value> map;
    Value val;
    val.a = 8;
    map[1] = val;
    std::cout << "a = " << map[1].a << std::endl;
    map[1].b = 16;
    std::cout << "a = " << map[1].a << std::endl;
    std::cout << "b = " << map[1].b << std::endl;
    
    return 0;
}
