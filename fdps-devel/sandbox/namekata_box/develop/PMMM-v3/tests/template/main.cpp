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

template <int p>
class MultipoleMoment {
public:
    enum {
        order = p,
        length = (p+1)*(p+1),
    };
    using value_type = double;
    value_type buf[length];
};

class FP_nbody {
public:
    double charge;
    double pos[3];
    double charge_sum;
    double pot;
    double grad[3];
};

class SPJMonopoleGeometricCenter{
public:
    double charge;
    double pos[3];
};

template <int p>
class SPJMultipoleGeometricCenter{
public:
    enum {
        order = p,
    };
    double charge;
    double pos[3];
    MultipoleMoment<p> mm;
    using value_type = typename MultipoleMoment<p>::value_type;
    constexpr int getMultipoleOrder() const {
        return order;
    }
};

// primary function template
template <class Tptcl>
void CalcForce(const FP_nbody * ep_i,
               const int n_ip,
               const Tptcl * ep_j,
               const int n_jp,
               FP_nbody * force) {
    for(int i = 0; i < n_ip; i++){
        force[i].charge_sum = 0.0;
        double sum {0.0};
        for(int j = 0; j < n_jp; j++){
            sum += ep_j[j].charge;
        }
        force[i].charge_sum += sum;
    }
}

// secondary function template (1)
template <>
void CalcForce<SPJMonopoleGeometricCenter>(
        const FP_nbody * ep_i,
        const int n_ip,
        const SPJMonopoleGeometricCenter * ep_j,
        const int n_jp,
        FP_nbody * force) {
    for(int i = 0; i < n_ip; i++){
        force[i].charge_sum = 0.0;
        double sum {0.0};
        for(int j = 0; j < n_jp; j++){
            sum += ep_j[j].charge;
        }
        force[i].charge_sum += sum;
    }
}

// function override
template <int p>
void CalcForce(const FP_nbody * ep_i,
               const int n_ip,
               const SPJMultipoleGeometricCenter<p> * ep_j,
               const int n_jp,
               FP_nbody * force) {
    for(int i = 0; i < n_ip; i++){
        force[i].charge_sum = 0.0;
        double sum {0.0};
        for(int j = 0; j < n_jp; j++){
            sum += ep_j[j].charge;
        }
        force[i].charge_sum += sum;
    }
}


template <int n, class real_t=double>
class TestClass {
public:
    using value_type = real_t;
    double x[n];
};


int main(int argc, char *argv[]) {

    std::cout << sizeof(TestClass<2>::value_type) << std::endl;
    std::cout << sizeof(SPJMultipoleGeometricCenter<2>::value_type) << std::endl;
    std::cout << SPJMultipoleGeometricCenter<3>::order << std::endl;
    SPJMultipoleGeometricCenter<3> spj;
    std::cout << spj.order << std::endl;
    constexpr int p_spj = spj.getMultipoleOrder();
    std::cout << sizeof(MultipoleMoment<p_spj>::value_type) << std::endl;
    return 0;
}
