#pragma once
#include <cstdint>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <array>
#include <vector>

template <class T>
T byteswap(const T val) {
    constexpr int size = sizeof(T);
    constexpr int block_size = sizeof(uint8_t);
    if (size > block_size) {
        assert(size % block_size == 0);
        constexpr int n_block = size / block_size;
        std::array<uint8_t *, n_block> block;
        T T_tmp = val;
        uint8_t * head = reinterpret_cast<uint8_t *>(&T_tmp);
        for (int i = 0; i < n_block; i++) block[i] = head + i;
        for (int i = 0; i < n_block/2; i++) {
            uint8_t high = *block[i];
            uint8_t low  = *block[n_block - 1 - i];
            *block[n_block - 1 - i] = high;
            *block[i] = low;
        }
        return T_tmp;
    } else {
        return val;
    }
}

uint64_t get_file_size(const std::string & file_name) {
    std::ifstream ifs;
    ifs.open(file_name.c_str(), std::ios::binary);
    ifs.seekg(0, std::ios::end);
    auto eofpos = ifs.tellg();
    ifs.clear();
    ifs.seekg(0, std::ios::beg);
    auto begpos = ifs.tellg();
    auto size = eofpos - begpos;
    ifs.close();
    return size;
}

void readParticleData(PS::ParticleSystem<FP_nbody> & psys,
                      const std::string & file_name_base) {
    const PS::F64 time_offset = PS::GetWtime();
    constexpr PS::S32 ptcl_data_size = 66;

    // File name
    std::stringstream ss;
    ss << file_name_base << std::setfill('0') << std::setw(5) << PS::Comm::getRank();
    const std::string file_name = ss.str();

    // Calculate # of particles in the file
    PS::U64 file_size = get_file_size(file_name);
    PS::S64 n_loc = file_size/ptcl_data_size;
    if (PS::Comm::getRank() == 0) {
        std::cout << n_loc << " particles are stored in file " << file_name
                  << " (file size is " << file_size << " [B])" << std::endl;
    }
    psys.setNumberOfParticleLocal(n_loc);

    // Read particle data
    std::ifstream input_file;
    input_file.open(file_name.c_str(), std::ios::in | std::ios::binary);
    for (PS::S32 i = 0; i < n_loc; i++) {
        // < Data format >
        //  (short) Type, (double) x [kpc], (double) y [kpc], (double) z [kpc], (double) vx [km/s], (double) vy [km/s], (double) vz [km/s], (double) timestep [year], (double) h [kpc] 
        int16_t type;
        double pos[3];
        double vel[3];
        double dt;
        double h;
        input_file.read((char *)&type, sizeof(int16_t));
        input_file.read((char *)pos, 3 * sizeof(double));
        input_file.read((char *)vel, 3 * sizeof(double));
        input_file.read((char *)&dt, sizeof(double));
        input_file.read((char *)&h, sizeof(double));
        psys[i].id    = i;
#ifdef READ_DATA_WITH_BYTESWAP
        psys[i].type  = byteswap(type); 
        psys[i].pos.x = byteswap(pos[0]);
        psys[i].pos.y = byteswap(pos[1]);
        psys[i].pos.z = byteswap(pos[2]);
        psys[i].vel.x = byteswap(vel[0]);
        psys[i].vel.y = byteswap(vel[1]);
        psys[i].vel.z = byteswap(vel[2]);
        psys[i].dt    = byteswap(dt); 
        psys[i].h     = byteswap(h); 
#else
        psys[i].type  = type; 
        psys[i].pos.x = pos[0];
        psys[i].pos.y = pos[1];
        psys[i].pos.z = pos[2];
        psys[i].vel.x = vel[0];
        psys[i].vel.y = vel[1];
        psys[i].vel.z = vel[2];
        psys[i].dt    = dt; 
        psys[i].h     = h; 
#endif
#if 0
        std::cout << psys[i].type << " "
                  << psys[i].pos.x << " "
                  << psys[i].pos.y << " "
                  << psys[i].pos.z << " "
                  << psys[i].vel.x << " "
                  << psys[i].vel.y << " "
                  << psys[i].vel.z << " "
                  << psys[i].dt << " "
                  << psys[i].h << std::endl;
#endif
    }
    input_file.close();
    if (PS::Comm::getRank() == 0) {
        std::cout << PS::GetWtime() - time_offset << " [s] were required to read a file." << std::endl;
    }
    //PS::Finalize();
    //std::exit(0);
}
