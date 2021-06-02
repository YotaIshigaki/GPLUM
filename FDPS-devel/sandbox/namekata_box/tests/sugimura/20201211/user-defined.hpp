#pragma once
// Include FDPS header
#include <particle_simulator.hpp>

constexpr PS::S32 N_NGB_TYPICAL = 100;
constexpr PS::S32 N_NGB_MAX = 2 * N_NGB_TYPICAL;

class Force {
public:
    PS::S32 id[N_NGB_MAX]; // j-particle's ID
    PS::F64 val[N_NGB_MAX]; // off-diagonal values

    void clear() {
        for (PS::S32 i=0; i<N_NGB_MAX; i++) {
            id[i] = -1;
            val[i] = 0.0;
        }
    }
};

class FullParticle {
public:
    PS::S32 id; // i-particle's ID
    Force force; 
};

class CommBuffer {
public:
    PS::S32 i; // i-particle's ID
    PS::S32 j; // j-particle's ID
    PS::F64 val; // j-particle's force data
};
