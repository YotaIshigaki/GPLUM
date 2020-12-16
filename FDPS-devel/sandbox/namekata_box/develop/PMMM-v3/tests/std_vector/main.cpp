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

    // Copy test
    std::vector<int> list_a;
    const int imax=8;
    for (int i=0; i<imax; i++) list_a.push_back(i);

    std::vector<int> list_b;
    list_b = list_a;

    for (int i=0; i<static_cast<int>(list_b.size()); i++)
        std::cout << "i = " << i << ", list_b[i] = " << list_b[i] << std::endl;

    // resize & push_back test
    std::vector<int> list_c;
    const int n_data_1st = 16;
    for (int i=0; i<n_data_1st; i++) list_c.push_back(i);
    std::cout << "list_c.size()     = " << list_c.size() << std::endl;
    std::cout << "list_c.capacity() = " << list_c.capacity() << std::endl;
    list_c.resize(0);
    std::cout << "list_c is resized!" << std::endl;
    std::cout << "list_c.size()     = " << list_c.size() << std::endl;
    std::cout << "list_c.capacity() = " << list_c.capacity() << std::endl;
    const int n_data_2nd = 8;
    for (int i=0; i<n_data_2nd; i++) list_c.push_back(i);
    std::cout << "list_c is reset." << std::endl;
    std::cout << "list_c.size()     = " << list_c.size() << std::endl;
    std::cout << "list_c.capacity() = " << list_c.capacity() << std::endl;

    return 0;
}
