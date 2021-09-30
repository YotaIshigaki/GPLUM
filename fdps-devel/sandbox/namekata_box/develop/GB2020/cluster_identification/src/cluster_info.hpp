#pragma once
#include <iostream>
#include <vector>

namespace ParticleSimulator {

    class ClusterInfo {
    public:
        U64 id;
        F64ort vertex;
        std::vector<S32> ranks;
    };

} // END of namespace of ParticleSimulator
