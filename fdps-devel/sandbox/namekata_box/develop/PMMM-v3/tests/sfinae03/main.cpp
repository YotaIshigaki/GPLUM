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

// getMultipoleOrder
template<typename Tptcl>
struct HasgetMultipoleOrderMethod
{  
    template<typename U, int(*)() > struct SFINAE {};
    template<typename U> static char Test(SFINAE<U, &U::getMultipoleOrder> *);
    template<typename U> static int Test(...);
    static const bool value = sizeof(Test<Tptcl>(0)) == sizeof(char);
};
template<class Tptcl>
int GetMyMultipoleOrder(std::true_type)
{   
    return Tptcl::getMultipoleOrder();
}
template<class Tptcl>
int GetMyMultipoleOrder(std::false_type)
{   
    return -1;
}
template <class Tptcl>
int GetMyMultipoleOrder() { // wrapper function
    return GetMyMultipoleOrder<Tptcl>(std::integral_constant<bool, HasgetMultipoleOrderMethod<Tptcl>::value>());
}

// getMultipole
template<typename Tspj, typename Tbuf>
struct HasgetMultipoleMethod
{  
    template<typename U, void(U::*)(int, Tbuf *) > struct SFINAE0 {};
    template<typename U, void(U::*)(int, Tbuf *) const > struct SFINAE1 {};
    template<typename U> static char Test(SFINAE0<U, &U::getMultipole> *);
    template<typename U> static char Test(SFINAE1<U, &U::getMultipole> *);
    template<typename U> static int Test(...);
    static const bool value = sizeof(Test<Tspj>(0)) == sizeof(char);
};
template<class Tspj, class Tbuf>
void GetMyMultipole(Tspj & spj, int buflen, Tbuf buf[], std::true_type)
{   
    spj.getMultipole(buflen, buf);
}
template<class Tspj, class Tbuf>
void GetMyMultipole(Tspj & spj, int buflen, Tbuf buf[], std::false_type)
{   
    std::cout << "hoge" << std::endl;
    for (int i=0; i<buflen; i++) buf[i]=0.0;
}
template <class Tspj, class Tbuf>
void GetMyMultipole(Tspj & spj, int buflen, Tbuf buf[]) { // wrapper function
    GetMyMultipole(spj, buflen, buf, std::integral_constant<bool, HasgetMultipoleMethod<Tspj, Tbuf>::value>());
}

class FP{
public:
    long long int id;
    double r_search;
    long long int getId() const {
        return this->id;
    }
    double getRSearch() const{
        return this->r_search;
    }
};


template <int p>
class SPJMultipoleGeometricCenter{
public:
    enum {
        order = p,
        length = (p+1)*(p+1),
    };
    double charge;
    double pos[3];
    double mm[length];
    static int getMultipoleOrder() { return p; }
    template <class Tbuf>
    void getMultipole(int buflen, Tbuf *buf) {
        assert(buflen <= length);
        for (int i=0; i<buflen; i++) buf[i] = mm[i];
    }
};

class Func {
public:
    template <typename Tbuf>
    void func(int buflen, Tbuf *buf) {
    }
};


template <typename Tbuf>
void func(int buflen, Tbuf *buf){
}

int main(int argc, char *argv[]) {

    using SPJ_t = SPJMultipoleGeometricCenter<3>;
    // Test GetMyMultipoleOrder
    const int p_spj = GetMyMultipoleOrder<SPJ_t>();
    std::cout << "p_spj = " << p_spj << std::endl;

    // Test GetMyMultipole
    SPJ_t spj;
    for (int i=0; i<SPJ_t::length; i++) spj.mm[i] = 1.0;
    constexpr int buflen = 4;
    double *buf = new double[buflen];
    GetMyMultipole(spj, buflen, buf);
    for (int i=0; i<buflen; i++)
        std::cout << "i = " << i << " buf = " << buf[i] << std::endl;
    delete [] buf;


    return 0;
}
