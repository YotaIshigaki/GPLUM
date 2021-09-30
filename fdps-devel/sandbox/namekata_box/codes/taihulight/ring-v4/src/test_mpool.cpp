#include<iostream>
#include"memory_pool.hpp"
//#include"particle_simulator.hpp"

using namespace ParticleSimulator;

struct A{
    int x;
    double y;
};

struct B{
    int x;
    float y[5];
};

struct C{
    A x;
    B y;
};

int main(){
    MemoryPool::initialize(100000);
    MemoryPool::dump();
    TempArray<A> a;
    a.resizeNoInitialize(100);
    std::cerr<<"--- after reserve ---"<<std::endl;
    a.dump();
    std::cerr<<"--- before free ---"<<std::endl;
    a.free();
    std::cerr<<"--- after free ---"<<std::endl;
    a.dump();

    /*
    const int n = 10;
    TempArray<A> a[1000];
    TempArray<B> b[1000];
    TempArray<C> c[1000];
    for(size_t i=0; i<n; i++){
        a[i].reserve(10);
        b[i].reserve(10);
        c[i].reserve(10);
    }
    std::cerr<<"--- reserve ---"<<std::endl;
    MemoryPool::dump();
        
    for(size_t i=0; i<n; i++){
        b[i].free();
    }
    std::cerr<<"--- free b ---"<<std::endl;
    MemoryPool::dump();

    for(size_t i=n; i<2*n; i++){
        c[i].reserve(10);
    }
    std::cerr<<"--- reserve c ---"<<std::endl;
    MemoryPool::dump();
    
    for(size_t i=n; i<2*n; i++){
        a[i].reserve(10);
    }
    std::cerr<<"--- reserve a ---"<<std::endl;
    MemoryPool::dump();
    */

    
    return 0;
}
