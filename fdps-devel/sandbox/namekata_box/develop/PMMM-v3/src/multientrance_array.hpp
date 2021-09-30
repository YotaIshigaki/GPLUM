#pragma once
#include <iostream>
#include <iomanip>
#include <vector>

namespace ParticleSimulator {

    template <class T>
    class MultientranceArray {
    public:
        std::vector<int> counts, displs;
        std::vector<T> data;

    };

} // END of namespace of ParticleSimulator
