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

// getId
template<typename Tptcl>
class HasgetIdMethod
{  
    template<class U,
        class = typename std::enable_if<
            (std::is_same<decltype(&U::getId), int(U::*)()>::value ||
             std::is_same<decltype(&U::getId), int(U::*)() const>::value ||
             std::is_same<decltype(&U::getId), long long int(U::*)()>::value ||   
             std::is_same<decltype(&U::getId), long long int(U::*)() const>::value
            ) &&
            std::is_member_pointer<decltype(&U::getId)>::value
                >::type
    >
    static std::true_type check(int);
    template <class U>
    static std::false_type check(...);
public:
    static constexpr bool value = decltype(check<Tptcl>(0))::value;
};
template<class Tptcl>
long long int GetMyId(Tptcl & ptcl, std::true_type)
{
   return ptcl.getId();
}
template<class Tptcl>
long long int GetMyId(Tptcl & ptcl, std::false_type) 
{
   return -1;
}
template <class Tptcl>
long long int GetMyId(Tptcl & ptcl) { // wrapper function
   return GetMyId(ptcl, std::integral_constant<bool, HasgetIdMethod<Tptcl>::value>());
}

// getRSearch
template<typename Tptcl>
class HasgetRSearchMethod
{
    template<class U,
        class = typename std::enable_if<
            (std::is_same<decltype(&U::getRSearch), float(U::*)()>::value ||
             std::is_same<decltype(&U::getRSearch), float(U::*)() const>::value ||
             std::is_same<decltype(&U::getRSearch), double(U::*)()>::value ||
             std::is_same<decltype(&U::getRSearch), double(U::*)() const>::value
            ) &&
            std::is_member_pointer<decltype(&U::getRSearch)>::value
                >::type
    >
    static std::true_type check(int);
    template <class U>
    static std::false_type check(...);
public:
    static constexpr bool value = decltype(check<Tptcl>(0))::value;
};
template<class Tptcl>
double GetMyRSearch(Tptcl ptcl, std::true_type)
{  
   return ptcl.getRSearch();
}
template<class Tptcl>
double GetMyRSearch(Tptcl ptcl, std::false_type)
{  
   return 0.0;
}
template <class Tptcl>
double GetMyRSearch(Tptcl ptcl) {
   return GetMyRSearch(ptcl, std::integral_constant<bool, HasgetRSearchMethod<Tptcl>::value>());
}

#if 0
// getDipole
template<typename Tspj>
struct HasgetDipoleMethod
{
   template<typename U, F32vec(U::*)() > struct SFINAE0 {};
   template<typename U, F32vec(U::*)() const > struct SFINAE1 {};
   template<typename U, F64vec(U::*)() > struct SFINAE2 {};
   template<typename U, F64vec(U::*)() const > struct SFINAE3 {};
   template<typename U> static char Test(SFINAE0<U, &U::getDipole> *);
   template<typename U> static char Test(SFINAE1<U, &U::getDipole> *);
   template<typename U> static char Test(SFINAE2<U, &U::getDipole> *);
   template<typename U> static char Test(SFINAE3<U, &U::getDipole> *);
   template<typename U> static int Test(...);
   static const bool value = sizeof(Test<Tspj>(0)) == sizeof(char);
};
template<class Tspj>
F64vec GetMyDipole(Tspj & spj, std::true_type)
{
   return spj.getDipole();
}
template<class Tspj>
F64vec GetMyDipole(Tspj & spj, std::false_type) 
{
   return F64vec(0.0);
}
template <class Tspj>
F64vec GetMyDipole(Tspj & spj) { // wrapper function
   return GetMyDipole(spj, std::integral_constant<bool, HasgetDipoleMethod<Tspj>::value>());
}

// getQuadrupole
template<typename Tspj>
struct HasgetQuadrupoleMethod
{
   template<typename U, F32mat(U::*)() > struct SFINAE0 {};
   template<typename U, F32mat(U::*)() const > struct SFINAE1 {};
   template<typename U, F64mat(U::*)() > struct SFINAE2 {};
   template<typename U, F64mat(U::*)() const > struct SFINAE3 {};
   template<typename U> static char Test(SFINAE0<U, &U::getQuadrupole> *);
   template<typename U> static char Test(SFINAE1<U, &U::getQuadrupole> *);
   template<typename U> static char Test(SFINAE2<U, &U::getQuadrupole> *);
   template<typename U> static char Test(SFINAE3<U, &U::getQuadrupole> *);
   template<typename U> static int Test(...);
   static const bool value = sizeof(Test<Tspj>(0)) == sizeof(char);
};
template<class Tspj>
F64mat GetMyQuadrupole(Tspj & spj, std::true_type)
{
   return spj.getQuadrupole();
}
template<class Tspj>
F64mat GetMyQuadrupole(Tspj & spj, std::false_type) 
{
   return F64mat(0.0);
}
template <class Tspj>
F64mat GetMyQuadrupole(Tspj & spj) { // wrapper function
   return GetMyQuadrupole(spj, std::integral_constant<bool, HasgetQuadrupoleMethod<Tspj>::value>());
}
#endif


// getMultipoleOrder
template<typename Tspj>
struct HasgetMultipoleOrderMethod
{ 
    template<class U, class = typename std::enable_if<!std::is_member_pointer<decltype(&U::getMultipoleOrder)>::value>::type>
    static std::true_type check(int);
    template <class U>
    static std::false_type check(...);
public:
    static constexpr bool value = decltype(check<Tspj>(0))::value;
};

template<class Tspj>
constexpr int GetMyMultipoleOrder(std::true_type)
{  
   return Tspj::getMultipoleOrder();
}
template<class Tspj>
constexpr int GetMyMultipoleOrder(std::false_type)
{  
   return -1;
}
template <class Tspj>
constexpr int GetMyMultipoleOrder() { // wrapper function
   return GetMyMultipoleOrder<Tspj>(std::integral_constant<bool, HasgetMultipoleOrderMethod<Tspj>::value>());
}

// getMultipole
#if 1
template<typename Tspj, typename Tbuf>
class HasgetMultipoleMethod
{  
    template<class U, class V,
        class = typename std::enable_if<
            std::is_same<decltype(&U::template getMultipole<V>), void(U::*)(int, V *)>::value &&
            std::is_member_pointer<decltype(&U::template getMultipole<V>)>::value
                >::type
    >
    static std::true_type check(int);
    template <class U, class V>
    static std::false_type check(...);
public:
    static constexpr bool value = decltype(check<Tspj,Tbuf>(0))::value;
};
template<class Tspj, class Tbuf>
void GetMyMultipole(Tspj & spj, int buflen, Tbuf buf[], std::true_type)
{   
    spj.getMultipole(buflen, buf);
}
template<class Tspj, class Tbuf>
void GetMyMultipole(Tspj & spj, int buflen, Tbuf buf[], std::false_type)
{   
    for (int i=0; i<buflen; i++) buf[i]=Tbuf(0);
}
template <class Tspj, class Tbuf>
void GetMyMultipole(Tspj & spj, int buflen, Tbuf buf[]) { // wrapper function
    GetMyMultipole(spj, buflen, buf, std::integral_constant<bool, HasgetMultipoleMethod<Tspj, Tbuf>::value>());
}
#else
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
#endif

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
    static constexpr int getMultipoleOrder() { return p; }
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

    // Test GetMyId, GetMyRSearch
    FP ptcl;
    ptcl.id = 1;
    ptcl.r_search = 1.0;
    std::cout << "GetMyId      = " << GetMyId(ptcl) << std::endl;
    std::cout << "GetMyRSearch = " << GetMyRSearch(ptcl) << std::endl;


    constexpr int p = 3;
    using SPJ_t = SPJMultipoleGeometricCenter<p>;
    // Test GetMyMultipoleOrder
    constexpr int p_spj = GetMyMultipoleOrder<SPJ_t>();
    std::cout << "p_spj = " << p_spj << std::endl;

    // Check 
    static_assert(std::is_same<decltype(&func<double>), void(*)(int, double *)>::value);
    std::cout << "typeid(void(*)(int, double*)).name()        = " << typeid(void(*)(int, double *)).name() << std::endl;
    std::cout << "typeid(&func<double>).name()                = " << typeid(&func<double>).name() << std::endl;
    std::cout << "typeid(&Func::func<double>).name()          = " << typeid(&Func::func<double>).name() << std::endl;
    std::cout << "typeid(&SPJ_t::getMultipole).name()         = " << typeid(&SPJ_t::getMultipole<double>).name() << std::endl;
    std::cout << "typeid(void(SPJ_t::*)(int, double*)).name() = " << typeid(void(SPJ_t::*)(int, double *)).name() << std::endl;
    static_assert(std::is_same<decltype(&SPJ_t::getMultipole<double>), void(SPJ_t::*)(int, double *)>::value);
    // static_assert(,) terminates comlilation with outputing
    // a string of characters given by the second argument
    // if a condition given by the 1st argument is false.

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
