#pragma once
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "vector3.h"

struct Particle{
    int    id;
    dvec3  pos;
    double mass;
    dvec3  acc_direct;
    double phi_direct;
    dvec3  acc_app;
    double phi_app;

    void clear_pot(){
        acc_direct = dvec3(0.0);
        phi_direct = 0.0;
        acc_app    = dvec3(0.0);
        phi_app    = 0.0;
    }
    void move_accp(){
        acc_app += acc_direct;
        phi_app += phi_direct;
        acc_direct = 0.0;
        phi_direct = 0.0;
    }
    dvec3 avecdiff_rel() const {
        return (acc_app - acc_direct) / acc_direct.abs();
    }
    double adiff_rel() const {
        return (acc_app - acc_direct).abs() / acc_direct.abs();
    }
    double adiff_abs() const {
        return (acc_app - acc_direct).abs();
    }
    double aabs() const {
        return acc_direct.abs();
    }
    double pdiff_rel() const {
        return fabs((phi_app - phi_direct) / phi_direct);
    }
    double pdiff_abs() const {
        return fabs(phi_app - phi_direct);
    }
    double pabs() const {
        return fabs(phi_direct);
    }

    static dvec3 rand_vec(){
        return dvec3(drand48(), drand48(), drand48());
    }

    static void gen_rand_dist(
            const int NP,
            Particle ptcl[],
            const long seed = 19810614)
    {
#if defined(READ_INITIAL_DATA)
        std::cout << "Reading initial particle data." << std::endl;
        {
            const std::string filename = "ic.dat";
            std::ifstream input_file;
            input_file.open(filename.c_str(), std::ios::binary);
            for (int i=0; i<NP; i++) {
                long long int id_tmp;
                double mass_tmp;
                double x_tmp, y_tmp, z_tmp;
                input_file.read((char *)&id_tmp, sizeof(long long int));
                input_file.read((char *)&mass_tmp, sizeof(double));
                input_file.read((char *)&x_tmp, sizeof(double));
                input_file.read((char *)&y_tmp, sizeof(double));
                input_file.read((char *)&z_tmp, sizeof(double));
                ptcl[i].id = id_tmp;
                ptcl[i].mass = mass_tmp;
                ptcl[i].pos.x = x_tmp;
                ptcl[i].pos.y = y_tmp;
                ptcl[i].pos.z = z_tmp;
            }
            input_file.close();
        }
#else
        srand48(seed);

        for(int i=0; i<NP; i++){
            ptcl[i].id  = i;
            ptcl[i].pos = rand_vec();
            ptcl[i].clear_pot();
        }
        const double pi = 4.0 * std::atan(1.0);
        double msum = 0.0;
        for(int i=0; i<NP; i++){
#if defined(CHECK_COSMO_SIM)
            msum += 
                (ptcl[i].mass  = 3.0/(8.0*pi*NP));
#else
            msum += 
                (ptcl[i].mass  = drand48() * (1.0/NP));
#endif
        }
#ifndef NON_CHARGE_NEUTRAL
        for(int i=0; i<NP; i++){
            ptcl[i].mass -= msum / NP;
        }
#endif
#endif
        // Output
        std::string filename = "IC.txt";
        std::ofstream output_file;
        output_file.open(filename.c_str(), std::ios::trunc);
        for (int i = 0; i < NP; i++) {
            output_file << ptcl[i].pos.x << " "
                        << ptcl[i].pos.y << " "
                        << ptcl[i].pos.z << " "
                        << ptcl[i].mass <<  std::endl;
        }
        output_file.close();
    }
};
