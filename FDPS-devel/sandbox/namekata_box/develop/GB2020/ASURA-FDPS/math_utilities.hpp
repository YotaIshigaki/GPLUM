#pragma once
/* C++ headers */
#include "common.h"
/* FDPS header */
#include <particle_simulator.hpp>
/* user-defined headers */
#include "macro_defs.h"

namespace math_utilities {

    // non-MPI version
    template <class real_t>
    real_t qmc3d(std::function<real_t(const real_t, const real_t, const real_t)> func,
                 const real_t xmin, const real_t xmax,
                 const real_t ymin, const real_t ymax,
                 const real_t zmin, const real_t zmax,
                 const std::size_t N, 
                 std::mt19937_64 & mt) {
        const real_t dv = (xmax - xmin) * (ymax - ymin) * (zmax - zmin) / N;
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        PS::F64 sum {0.0};
        for (int i = 0; i < N; i++) {
            const real_t x = xmin + (xmax - xmin) * dist(mt);
            const real_t y = ymin + (ymax - ymin) * dist(mt);
            const real_t z = zmin + (zmax - zmin) * dist(mt);
            sum += func(x,y,z);
        }
        return sum * dv;
    }

    // MPI-parallel version
    template <class real_t>
    real_t qmc3d(std::function<real_t(const real_t, const real_t, const real_t)> func,
                 const real_t xmin, const real_t xmax,
                 const real_t ymin, const real_t ymax,
                 const real_t zmin, const real_t zmax,
                 const std::size_t N_glb, 
                 std::mt19937_64 & mt,
                 const PS::CommInfo & comm_info) {
        const int n_proc = comm_info.getNumberOfProc();
        const int my_rank = comm_info.getRank();
        const int rem = N_glb % n_proc; // remainder
        const int quo = N_glb / n_proc; // quotient
        const int N_loc = (my_rank < rem) ? (quo + 1) : quo;
        const real_t dv = (xmax - xmin) * (ymax - ymin) * (zmax - zmin) / N_glb;
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        PS::F64 sum {0.0};
        for (int i = 0; i < N_loc; i++) {
            const real_t x = xmin + (xmax - xmin) * dist(mt);
            const real_t y = ymin + (ymax - ymin) * dist(mt);
            const real_t z = zmin + (zmax - zmin) * dist(mt);
            sum += func(x,y,z);
        }
        return comm_info.getSum(sum) * dv;
    }

}
namespace math_util = math_utilities;
