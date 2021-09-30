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

static inline double factorial(const int i){
    assert(i >= 0);
    return i ? double(i) * factorial(i-1) : 1.0;
}
static inline double factinv(const int i){
    return 1.0 / factorial(i);
}

template <typename real_t=double>
class TableForRlm {
public:
    std::vector<real_t> tbl_inv;
    std::vector<real_t> tbl_factinv;

    TableForRlm() {}

    void init (int p) {
        for(int i=0; i<p+1; i++){
            tbl_inv.push_back(real_t(1.0) / real_t(i));
        }
        for(int i=0; i<2*p+1; i++){
            assert(factorial(i) > 0);
            tbl_factinv.push_back(1.0 / factorial(i));
        }
    }
};

template <typename Ttbl>
class TableManager {
public:
    std::unordered_map<int, int> map;
    std::vector<Ttbl> vec;

    TableManager() {}

    void set(int p) {
        const int adr = vec.size();
        Ttbl tmp; 
        tmp.init(p);
        vec.push_back(tmp);
        map[p] = adr;
    }

    Ttbl & get(int p){
        const int adr = map[p];
        return vec[adr];
    }

};


int main(int argc, char *argv[]) {

    TableManager<TableForRlm<double> > tbl_mgr;
  
#if 1 
    int p;

    p = 3; 
    tbl_mgr.set(p);
    TableForRlm<double> tbl;
    tbl = tbl_mgr.get(p);
    for (int i=0; i<2*p+1; i++) {
        std::cout << tbl.tbl_factinv[i] << std::endl;
    }
#endif

    return 0;
}
