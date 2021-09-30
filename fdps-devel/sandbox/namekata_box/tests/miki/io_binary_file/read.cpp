#include <iostream>
#include <fstream>
#include <array>
#include <cassert>
#include <cstdint>

typedef struct {
    double time;
    int nbodies;
    int ndim;
    int nsph;
    int ndark;
    int nstar;
} MAGI_Tipsy_Header;

typedef struct {
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    int idx;
} MAGI_Tipsy_Particle;

template <class T>
inline std::string GetBinString(const T val)
{
    constexpr int size = sizeof(T);
    if ( !val ) {
        std::string str = "0";
        constexpr int n_bit = 8 * size;
        for (int i = 1; i <= n_bit - 1; i++) str.append("0");
        return str;
    }
    constexpr int block_size = sizeof(uint16_t);
    constexpr int n_block = size / block_size;
    T T_tmp = val;
    uint16_t * head = reinterpret_cast<uint16_t *>(&T_tmp);
    std::string str;
    for (int n=0; n<n_block; n++) {
        uint16_t tmp = *(head + n_block - 1 - n);
        for (int i=0; i<16; i++) {
            if ( (tmp & 1) == 0 ) str.insert(str.begin(), '0');
            else str.insert(str.begin(), '1');
            tmp >>= 1;
        }
    }
    return str;
}

template <class T>
T byteswap(const T val) {
    constexpr int size = sizeof(T);
    constexpr int block_size = sizeof(uint16_t);
    if (size > block_size) {
        assert(size % block_size == 0);
        constexpr int n_block = size / block_size;
        std::array<uint16_t *, n_block> block;
        T T_tmp = val;
        uint16_t * head = reinterpret_cast<uint16_t *>(&T_tmp);
        for (int i = 0; i < n_block; i++) block[i] = head + i;
        for (int i = 0; i < n_block/2; i++) {
            uint16_t high = *block[i];
            uint16_t low  = *block[n_block - 1 - i];
            *block[n_block - 1 - i] = (high >> 8) | (high << 8);
            *block[i] = (low >> 8) | (low << 8);
        }
        return T_tmp;
    } else {
        return val;
    }
}

int main(int argc, char *argv[]) {

    MAGI_Tipsy_Header header;
    MAGI_Tipsy_Particle ptcl;
   
    const std::string file_name = "data.dat"; 
    std::ifstream input_file;
    input_file.open(file_name.c_str(), std::ios::in | std::ios::binary); 
    input_file.read((char *) &header, sizeof(header));
    input_file.read((char *) &ptcl, sizeof(ptcl));
    input_file.close();

    // byte swap
    header.time    = byteswap(header.time);
    header.nbodies = byteswap(header.nbodies);
    header.ndim    = byteswap(header.ndim);
    header.nsph    = byteswap(header.nsph);
    header.ndark   = byteswap(header.ndark);
    header.nstar   = byteswap(header.nstar);

    ptcl.mass   = byteswap(ptcl.mass);
    ptcl.pos[0] = byteswap(ptcl.pos[0]);
    ptcl.pos[1] = byteswap(ptcl.pos[1]);
    ptcl.pos[2] = byteswap(ptcl.pos[2]);
    ptcl.vel[0] = byteswap(ptcl.vel[0]);
    ptcl.vel[1] = byteswap(ptcl.vel[1]);
    ptcl.vel[2] = byteswap(ptcl.vel[2]);
    ptcl.eps    = byteswap(ptcl.eps);
    ptcl.idx    = byteswap(ptcl.idx);


    std::cout << "header.time    = " << header.time    << std::endl;
    std::cout << "header.nbodies = " << header.nbodies << std::endl;
    std::cout << "header.ndim    = " << header.ndim    << std::endl;
    std::cout << "header.nsph    = " << header.nsph    << std::endl;
    std::cout << "header.ndark   = " << header.ndark   << std::endl;
    std::cout << "header.nstar   = " << header.nstar   << std::endl;
    std::cout << "-----" << std::endl;
    std::cout << "ptcl.mass   = " << ptcl.mass   << std::endl;
    std::cout << "ptcl.pos[0] = " << ptcl.pos[0] << std::endl;
    std::cout << "ptcl.pos[1] = " << ptcl.pos[1] << std::endl;
    std::cout << "ptcl.pos[2] = " << ptcl.pos[2] << std::endl;
    std::cout << "ptcl.vel[0] = " << ptcl.vel[0] << std::endl;
    std::cout << "ptcl.vel[1] = " << ptcl.vel[1] << std::endl;
    std::cout << "ptcl.vel[2] = " << ptcl.vel[2] << std::endl;
    std::cout << "ptcl.eps    = " << ptcl.eps    << std::endl;
    std::cout << "ptcl.idx    = " << ptcl.idx    << std::endl;

    return 0;
}
