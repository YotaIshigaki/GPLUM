#include <cassert>
#include <iostream>
#include <vector>
#include "particle_simulator.hpp"

/* FullParticle */
class FP {
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    std::vector<int> list;

    PS::F64vec getPos() const { return this->pos; }
    void setPos(const PS::F64vec & pos_) {
        this->pos = pos_; 
    }
    void copyFromFP(const FP & fp) {
        this->id   = fp.id;
        this->mass = fp.mass;
        this->pos  = fp.pos;
    }
    void copyFromForce(const FP & force) {
    }
};

template <class T>
void checkStatus(PS::ReallocatableArray<T> & ra) {

    std::cout << "size()       = " << ra.size() << std::endl;
    std::cout << "capacity()   = " << ra.capacity() << std::endl;
    std::cout << "getMemSize() = " << ra.getMemSize() << std::endl;
    std::cout << "getPointer() = 0x" 
              << std::hex << std::setw(16) << std::setfill('0')
              << reinterpret_cast<unsigned long long>(ra.getPointer())
              << std::dec << std::endl; 
    if (ra.size() > 0) {
        std::cout << "Output (i, ra[i].list.size(), ra[i].list, &ra[i].list)" << std::endl;
        for (int i=0; i<ra.size(); i++) {
           if (ra[i].list.size() > 0) {
               std::cout << i << " , " 
                         << ra[i].list.size() << " , ";
               for (int val : ra[i].list)
                  std::cout << val << " , ";
               std::cout << std::hex << std::setw(16) << std::setfill('0')
                         << reinterpret_cast<unsigned long long>(ra[i].list.data())
                         << std::dec << std::endl; 
           } else {
               std::cout << "... " << i << " skipped." << std::endl;
           }
        }
    }

}


int main(int argc, char **argv) {
   PS::Initialize(argc, argv);

   PS::ReallocatableArray<FP> fp_org;
   PS::ReallocatableArray<FP> fp_sorted;

   //=========================================
   // The list of public member functions
   //=========================================
   // size();
   // capacity();
   // reserve();
   // push_back();
   // resizeNoInitialize();
   // getMemSize();
   // getPointer();
   // pushBackNoCheck();
   // clearSize()
   // increaseSize()
   // decreaseSize()
   // reserveAtLeast()
   // reserveEmptyAreaAtLeast()
   // freeMem()
   // reallocMem()

   //* Check status
   std::cout << std::endl;
   std::cout << "*** check status of fp_org" << std::endl;
   checkStatus(fp_org);

   std::cout << std::endl;
   std::cout << "*** check status of fp_sorted" << std::endl;
   checkStatus(fp_sorted);

   //* Set fp_org
   std::cout << std::endl;
   std::cout << "*** Setup fp_org" << std::endl;
   PS::S32 numPtcl=4;
   fp_org.reserve(numPtcl);
   for (PS::S32 i=0; i<numPtcl; i++) {
       class FP fp;
       fp.id = i;
       for (PS::S32 j=0; j<i+1; j++) {
            fp.list.push_back(10*i + j);
       }
       fp_org.push_back(fp);
   }
   checkStatus(fp_org);

   //* Set fp_sorted
   std::cout << std::endl;
   std::cout << "*** Setup fp_sorted"  << std::endl;
   numPtcl = fp_org.capacity(); 
   std::cout << "numPtcl = " << numPtcl << std::endl;
   fp_sorted.reserve(numPtcl);
   fp_sorted.resizeNoInitialize(numPtcl);
   for (PS::S32 i=0; i<numPtcl; i++)
      fp_sorted[i] = fp_org[i];
   checkStatus(fp_sorted);


  //* Use the second reserve()
  //std::cout << "*** test 2" << std::endl;
  //num = 8;
  //fp_array.reserve(num);
  //checkStatus(fp_array);


  PS::Finalize();
  return 0;
}


