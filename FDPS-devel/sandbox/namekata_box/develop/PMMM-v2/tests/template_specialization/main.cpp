// Include the standard C++ headers
#include <cmath>
#include <math.h>
#include <cfloat>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <sys/stat.h>
#include <time.h>

struct TagKeyNormal{};
struct TagKeyMeshBased{};

template <typename T>
class MortonKey {
};

template <>
class MortonKey<TagKeyNormal> {
private:
    unsigned long long key_;
    MortonKey(){};
    ~MortonKey(){};
    MortonKey(const MortonKey &);
    MortonKey & operator = (const MortonKey &);
    static MortonKey & getInstance() {
        static MortonKey inst;
        return inst;
    }
public:
    static void initialize(unsigned long long key) {
        getInstance().key_ = key;
    }
    static unsigned long long getKey() {
        return getInstance().key_;
    }
};

template<>
class MortonKey<TagKeyMeshBased> {
private:
    unsigned long long key_;
    MortonKey(){};
    ~MortonKey(){};
    MortonKey(const MortonKey &);
    MortonKey & operator = (const MortonKey &);
    static MortonKey & getInstance() {
        static MortonKey inst;
        return inst;
    }
public:
    static void initialize(unsigned long long key) {
        getInstance().key_ = key;
    }
    static unsigned long long getKey() {
        return getInstance().key_;
    }
};


int main(int argc, char *argv[]) {
  
    unsigned long long val = 1; 
    MortonKey<TagKeyMeshBased>::initialize(val); 
    std::cout << MortonKey<TagKeyMeshBased>::getKey() << std::endl;

    return 0;
}
