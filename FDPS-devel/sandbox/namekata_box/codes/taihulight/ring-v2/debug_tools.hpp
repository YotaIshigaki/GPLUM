#pragma once
/* Headers of std. libs */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <limits>
/* Header of FDPS */
#include <particle_simulator.hpp>
/* User-defined header */
#include "user-defined.hpp"

/* Macros */
#if defined(CHECK_NAN)
#define NAN_TEST(psys,msg) \
    do {                           \
        nan_test((psys),(msg));    \
    } while(0);
#else
#define NAN_TEST(psys,msg) \
    do {} while(0);
#endif
//-----------
#if defined(CHECK_ANGULAR_MOMENTUM)
#define ANGULAR_MOMENTUM_TEST0(psys) \
    do {                            \
        angmom_test0((psys));        \
    } while(0);
#else
#define ANGULAR_MOMENTUM_TEST0(psys) \
    do {} while(0);
#endif
//-----------
#if defined(CHECK_ANGULAR_MOMENTUM)
#define ANGULAR_MOMENTUM_TEST(psys) \
    do {                            \
        angmom_test((psys));        \
    } while(0);
#else
#define ANGULAR_MOMENTUM_TEST(psys) \
    do {} while(0);
#endif

/* Global variables */
extern PS::S32 N_SATELLITE;                         
extern PS::ReallocatableArray<Satellite> SATELLITE; 
static PS::F64vec L0; // Total angular momentum of the system at t=0


template<class Tpsys>
void nan_test(const Tpsys & system, std::string msg) {

    PS::S32 num_nan_loc = 0;
    for (PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++) {
         if ( (std::isfinite(system[i].pos.x) != true) ||
              (std::isfinite(system[i].pos.y) != true) ||
              (std::isfinite(system[i].pos.z) != true) ||
              (std::isfinite(system[i].mass)  != true) ||
              (std::isfinite(system[i].vel.x) != true) ||
              (std::isfinite(system[i].vel.y) != true) ||
              (std::isfinite(system[i].vel.z) != true)) {
             num_nan_loc++;
         }
    }
    PS::S32 num_nan;
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0)
        std::cout << msg << " " << num_nan << std::endl;
    
}

template<class Tpsys>
void nan_test(const Tpsys system[], const PS::S32 n, std::string msg) {

    PS::S32 num_nan_loc = 0;
    //for (PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++) {
    for (PS::S32 i=0; i<n; i++) {
         if ( (std::isfinite(system[i].pos.x) != true) ||
              (std::isfinite(system[i].pos.y) != true) ||
              (std::isfinite(system[i].pos.z) != true) ||
              (std::isfinite(system[i].mass)  != true) ||
              (std::isfinite(system[i].vel.x) != true) ||
              (std::isfinite(system[i].vel.y) != true) ||
              (std::isfinite(system[i].vel.z) != true)) {
             num_nan_loc++;
         }
    }
    PS::S32 num_nan;
    MPI::COMM_WORLD.Allreduce(&num_nan_loc,&num_nan,1,MPI::INT,MPI::SUM);
    if (PS::Comm::getRank() == 0)
        std::cout << msg << " " << num_nan << std::endl;
}

template<class Tpsys>
void angmom_test0(const Tpsys & system) {

    PS::F64vec L_loc = 0.0;
    PS::S32 n_loc = system.getNumberOfParticleLocal();
    for (PS::S32 i=0; i<n_loc; i++) {
        L_loc += system[i].mass * (system[i].pos ^ system[i].vel);
    }
    for (PS::S32 i=0; i<N_SATELLITE; i++) {
        L_loc += SATELLITE[i].mass * (SATELLITE[i].pos ^ SATELLITE[i].vel);
    }
    MPI::COMM_WORLD.Allreduce(&L_loc.x,&L0.x,1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&L_loc.y,&L0.y,1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&L_loc.z,&L0.z,1,MPI::DOUBLE,MPI::SUM);
    if (PS::Comm::getRank() == 0) {
       std::cout << "L0.x = " << L0.x << std::endl;
       std::cout << "L0.y = " << L0.y << std::endl;
       std::cout << "L0.z = " << L0.z << std::endl;
    }
}

template<class Tpsys>
void angmom_test0(const Tpsys & system, const PS::S32 n_loc) {

    PS::F64vec L_loc = 0.0;
    for (PS::S32 i=0; i<n_loc; i++) {
        L_loc += system[i].mass * (system[i].pos ^ system[i].vel);
    }
    for (PS::S32 i=0; i<N_SATELLITE; i++) {
        L_loc += SATELLITE[i].mass * (SATELLITE[i].pos ^ SATELLITE[i].vel);
    }
    MPI::COMM_WORLD.Allreduce(&L_loc.x,&L0.x,1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&L_loc.y,&L0.y,1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&L_loc.z,&L0.z,1,MPI::DOUBLE,MPI::SUM);
    if (PS::Comm::getRank() == 0) {
       std::cerr << "L0.x = " << L0.x << std::endl;
       std::cerr << "L0.y = " << L0.y << std::endl;
       std::cerr << "L0.z = " << L0.z << std::endl;
    }
}

template<class Tpsys>
void angmom_test(const Tpsys & system) {

    PS::F64vec L_loc = 0.0;
    PS::S32 n_loc = system.getNumberOfParticleLocal();
    for (PS::S32 i=0; i<n_loc; i++) {
        L_loc += system[i].mass * (system[i].pos ^ system[i].vel);
    }
    for (PS::S32 i=0; i<N_SATELLITE; i++) {
        L_loc += SATELLITE[i].mass * (SATELLITE[i].pos ^ SATELLITE[i].vel);
    }
    PS::F64vec L_tot;
    MPI::COMM_WORLD.Allreduce(&L_loc.x,&L_tot.x,1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&L_loc.y,&L_tot.y,1,MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&L_loc.z,&L_tot.z,1,MPI::DOUBLE,MPI::SUM);
    if (PS::Comm::getRank() == 0) {
        if(L0.x != 0.0){
            std::cout << "L.x = " << L_tot.x << " " 
                      << "rel.err. = " << std::abs((L_tot.x-L0.x)/L0.x) << std::endl;
        }
        if(L0.y != 0.0){
            std::cout << "L.y = " << L_tot.y << " " 
                      << "rel.err. = " << std::abs((L_tot.y-L0.y)/L0.y) << std::endl;
        }
        if(L0.z != 0.0){
            std::cout << "L.z = " << L_tot.z << " "
                      << "rel.err. = " << std::abs((L_tot.z-L0.z)/L0.z) << std::endl;
        }
    }

}

template<class Tpsys>
void checkLoadBalanceOfNumberOfParticle(const Tpsys & system) {

   PS::S32 n_loc = system.getNumberOfParticleLocal();
   PS::S32 n_loc_min,n_loc_max;
   MPI::COMM_WORLD.Allreduce(&n_loc,&n_loc_min,1,MPI::INT,MPI::MIN);
   MPI::COMM_WORLD.Allreduce(&n_loc,&n_loc_max,1,MPI::INT,MPI::MAX);
   if (PS::Comm::getRank() == 0) {
      std::cout << "n_loc_min = " << n_loc_min << std::endl;   
      std::cout << "n_loc_max = " << n_loc_max << std::endl;   
   }
   if (n_loc == n_loc_min)
       std::cout << "rank_min = " << PS::Comm::getRank() << std::endl;
   if (n_loc == n_loc_max)
       std::cout << "rank_max = " << PS::Comm::getRank() << std::endl;
}

template<class Tpsys>
void checkRingWidth(const Tpsys & system) {

    PS::F64 r2min_loc,r2max_loc;
    r2min_loc =  1.0e300;
    r2max_loc = -1.0e300;
    for (PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++) {
         PS::F64 r2 = system[i].pos * system[i].pos;
         if (r2 < r2min_loc) r2min_loc = r2;
         if (r2 > r2max_loc) r2max_loc = r2;
        
    }
    PS::F64 r2min,r2max;
    MPI::COMM_WORLD.Allreduce(&r2min_loc,&r2min,1,MPI::DOUBLE,MPI::MIN);
    MPI::COMM_WORLD.Allreduce(&r2max_loc,&r2max,1,MPI::DOUBLE,MPI::MAX);
    if (PS::Comm::getRank() == 0) {
        std::cout << "rmin = " << std::sqrt(r2min) << std::endl;
        std::cout << "rmax = " << std::sqrt(r2max) << std::endl;
    }
}

