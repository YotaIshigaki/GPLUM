/* Headers */
#include "mathematical_constants.h"
#include "physical_constants.h"

PS::F64 getPosCellCenter(const PS::F64 pos_left_bound,
                         const PS::F64 pos_right_bound,
                         const PS::S32 number_of_cells,
                         const PS::S32 i) {
    assert(0 <= i && i < number_of_cells);
    assert(pos_left_bound < pos_right_bound);
    const PS::F64 dx = (pos_right_bound - pos_left_bound) / number_of_cells;
    if ( number_of_cells % 2 == 0) {
        // # of cells is even.
        if (i < number_of_cells/2) {
            return pos_left_bound + dx * (i + 0.5);
        }
        else {
            return pos_right_bound - dx * (number_of_cells - i - 0.5);
        }
    } 
    else {
        // # of cells is odd.
        const PS::F64 center = 0.5 * (pos_left_bound + pos_right_bound);
        if (i < number_of_cells/2) {
            return center - dx * (number_of_cells/2 - i);
        }
        else {
            return center + dx * (i - number_of_cells/2);
        }
    }
}


void EvrardTestIC(PS::ParticleSystem<FP_sph>& psys,
                  PS::BOUNDARY_CONDITION& bc,
                  PS::F64ort& pos_root_domain,
                  PS::F64 & time_dump,
                  PS::F64 & dt_dump,
                  PS::F64 & time_end,
                  PS::S32 gen_mode) {
    // Place SPH particles
    PS::S64 N_sph = 0;
    const PS::F64 M_sph = 1.0;
    const PS::F64 radius_of_sphere = 1.0;
    if (gen_mode == 0) {
        // In this mode, we create an initial distribution of particles
        // by rescaling the positions of particles which are placed in a grid.
        // (1) Count # of particles in a sphere of radius 1
        const PS::S64 N_ptcl_per_side = 64;
        PS::F64 x,y,z;
        // x-loop
        for (PS::S32 i = 0; i < N_ptcl_per_side; i++) {
            x = getPosCellCenter(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, i);
            // y-loop
            for (PS::S32 j = 0; j < N_ptcl_per_side; j++) {
                y = getPosCellCenter(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, j);
                // z-loop
                for (PS::S32 k = 0; k < N_ptcl_per_side; k++) {
                    z = getPosCellCenter(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, k);
                    const PS::F64 r = std::sqrt(x*x + y*y + z*z);
                    if (r <= radius_of_sphere) N_sph++;
                }
            }
        }
        assert(N_sph > 0);
        psys.setNumberOfParticleLocal(N_sph);
        // (2) Actually place particles
        PS::S64 id = 0;
        // x-loop
        for (PS::S32 i = 0; i < N_ptcl_per_side; i++) {
            x = getPosCellCenter(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, i);
            // y-loop
            for (PS::S32 j = 0; j < N_ptcl_per_side; j++) {
                y = getPosCellCenter(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, j);
                // z-loop
                for (PS::S32 k = 0; k < N_ptcl_per_side; k++) {
                    z = getPosCellCenter(-radius_of_sphere, radius_of_sphere, N_ptcl_per_side, k);
                    const PS::F64 r = std::sqrt(x*x + y*y + z*z);
                    if (r <= radius_of_sphere) {
                        const PS::F64 r_new = radius_of_sphere * std::pow(r/radius_of_sphere,1.5);
                        psys[id].id = id;
                        psys[id].mass = M_sph / N_sph;
                        psys[id].pos.x = (r_new / r) * x;
                        psys[id].pos.y = (r_new / r) * y;
                        psys[id].pos.z = (r_new / r) * z;
                        psys[id].vel       = 0.0;
                        psys[id].acc_grav  = 0.0;
                        psys[id].pot_grav  = 0.0;
                        psys[id].acc_hydro = 0.0;
                        psys[id].eng = 0.05 * M_sph / radius_of_sphere;
                        psys[id].smth = radius_of_sphere * std::pow((PS::F64)N_neighbor/(PS::F64)N_sph, 1.0/3.0);
                        // Note that entropy is determined in setEntropy().
                        id++;
                    }
                }
            }
        }
    }
    else if (gen_mode == 1) {
        // In this mode, we set an initial distribution of particles by reading a file.
        std::vector<PS::F64vec> ptcl;
        // Read position data
        const std::string filename = "result/glass_data.txt";
        std::ifstream data_file; 
        data_file.open(filename.c_str(), std::ios::in | std::ios::binary);
        if (data_file) {
           while(true) {
              PS::F64vec pos;
              data_file.read((char *)&pos, sizeof(PS::F64vec));
              if (data_file.eof()) break;
              ptcl.push_back(pos);
           }
        }
        data_file.close();
        const PS::S64 n_ptcl_in = ptcl.size();
        if (PS::Comm::getRank() == 0)
            std::cout << n_ptcl_in << " particles are read." << std::endl;
        // Count # of particles in a sphere of radius 1
        N_sph = 0;
        for (PS::S32 i = 0; i < n_ptcl_in; i++) {
            const PS::F64 r2 = ptcl[i] * ptcl[i];
            if (r2 < 1.0) N_sph++;
        }
        if (PS::Comm::getRank() == 0)
            std::cout << N_sph << " particles of them will be used to make a Evrard sphere." << std::endl;
        // Place SPH particles
        psys.setNumberOfParticleLocal(N_sph);
        PS::S64 j = -1;
        for (PS::S64 i = 0; i < N_sph; i++) {
            psys[i].id = (i + 1);
            psys[i].mass = M_sph / N_sph;
            PS::F64 r2;
            for(;;) {
                j++;
                r2 = ptcl[j] * ptcl[j];
                if (r2 < 1.0) break;
            }
            const PS::F64 r = std::sqrt(r2);
            if (r > 0.0) {
                const PS::F64 r_new = radius_of_sphere * std::pow(r,1.5);
                psys[i].pos = (r_new / r) * ptcl[j];
            } else {
                psys[i].pos = ptcl[j];
            }
            psys[i].vel       = 0.0;
            psys[i].acc_grav  = 0.0;
            psys[i].pot_grav  = 0.0;
            psys[i].acc_hydro = 0.0;
            psys[i].eng = 0.05 * M_sph / radius_of_sphere;
            psys[i].smth = radius_of_sphere * std::pow((PS::F64)N_neighbor/(PS::F64)N_sph, 1.0/3.0);
            // Note that entropy is determined in setEntropy().
        }
    }
    else {
        if (PS::Comm::getRank() == 0) std::cout << "Given gen_mode is not supported." << std::endl;
        PS::Abort();
        std::exit(-1);
    }
    // Set boundary condition
    bc = PS::BOUNDARY_CONDITION_OPEN;
    // Set gravitational softening
    eps_grav = 1.0e-4 * radius_of_sphere;
    // Set I/O intervals
    const PS::F64 unit_dens = 3.0 * M_sph / (4.0 * math_const::pi * std::pow(radius_of_sphere, 3.0));
    const PS::F64 unit_time = std::sqrt(math_const::pi * math_const::pi / 8.0)
                            * std::pow(radius_of_sphere, 1.5)
                            / std::sqrt(M_sph);
    dt_dump = 0.1 * unit_time;
    time_dump = dt_dump;
    time_end = unit_time;
    if (PS::Comm::getRank() == 0) {
        std::cout << "unit_dens = " << unit_dens << std::endl;
        std::cout << "unit_time = " << unit_time << std::endl;
    }
    // Set maximum timestep
    dt_max = 1.0e-3;

    if (PS::Comm::getRank() == 0) {
        std::cout << "An initial condition for the Evrard test is made." << std::endl;
        std::cout << "(N_sph = " << N_sph << ")" << std::endl;
    }
}

void MakeGlassIC(PS::ParticleSystem<FP_sph>& psys,
                 PS::BOUNDARY_CONDITION& bc,
                 PS::F64ort& pos_root_domain,
                 PS::F64 & time_dump,
                 PS::F64 & dt_dump,
                 PS::F64 & time_end) {
    // Model parameters
    const PS::S64 N_sph = std::pow(2,18);
    // Initialize pseudorandom number generator
    PS::MTTS mt;
    mt.init_genrand(0);
    // Place SPH particles
    psys.setNumberOfParticleLocal(N_sph);
    for (PS::S32 i = 0; i < N_sph; i++) {
        psys[i].id = (i + 1);
        psys[i].mass = 8.0/N_sph;
        PS::F64 x,y,z;
        for(;;) {
            x = (2.0 * mt.genrand_res53() - 1.0);
            y = (2.0 * mt.genrand_res53() - 1.0);
            z = (2.0 * mt.genrand_res53() - 1.0);
            if ((-1.0 <= x) && (x < 1.0) &&
                (-1.0 <= y) && (y < 1.0) &&
                (-1.0 <= z) && (z < 1.0)) break;
        }
        psys[i].pos.x = x;
        psys[i].pos.y = y;
        psys[i].pos.z = z;
        psys[i].vel = 0.0;
        psys[i].eng = 1.0;  
        psys[i].smth = 2.0 * std::pow((PS::F64)N_neighbor/(PS::F64)N_sph, 1.0/3.0); // smoothing length
        // [Notes]
        //   (1) The value of the specific thermal energy is chosen
        //       so that the sound speed is nearly equal to 1.
        //   (2) the value of the entropy function is determined in setEntropy().
    }
    // Set boundary condition
    bc = PS::BOUNDARY_CONDITION_PERIODIC_XYZ;
    pos_root_domain.low_  = (PS::F64vec)(-1.0, -1.0, -1.0);
    pos_root_domain.high_ = (PS::F64vec)( 1.0,  1.0,  1.0);
    // Set gravitational softening
    eps_grav = 0.01;
    // Set I/O intervals
    const PS::F64 tcross = 2.0 * std::sqrt(3.0); 
    dt_dump   = 4.0 * tcross;
    time_dump = dt_dump;
    time_end  = 64.0 * tcross;
    if (PS::Comm::getRank() == 0) {
        std::cout << "The sound crossing time = " << tcross << std::endl;
    }
    // Set maximum timestep
    dt_max = 1.0e-2;

}

