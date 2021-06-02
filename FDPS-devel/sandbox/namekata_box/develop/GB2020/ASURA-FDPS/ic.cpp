/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "debug_utilities.h"
#include "mathematical_constants.h"
#include "physical_constants.h"
#include "run_parameters.h"
#include "user_defined.h"
#include "SPH_kernel.h"
#include "density_distribution_models.h"
#include "numerical_recipes.hpp"
#include "math_utilities.hpp"
#include "ic.h"

template<class T>
T reverseEndian(T value){
    char * first = reinterpret_cast<char*>(&value);
    char * last = first + sizeof(T);
    std::reverse(first, last);
    return value;
}

#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_NORMAL_RUN
/* Class definitions */
// The following two classes are used in function readTipsyFile,
// which is used to read particle data created by MAGI.
class MAGI_Tipsy_header {
public:
    double time;
    int nbodies;
    int ndim;
    int nsph;
    int ndark;
    int nstar;
};

class MAGI_Tipsy_particle {
public:
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    int idx;
};

void readTipsyFile(std::string& input_file_name,
                   PS::ParticleSystem<FP_dm>& psys) {
    static PS::S32 ONE = 1;
    static bool is_little_endian = *reinterpret_cast<char*>(&ONE) == ONE;
    // This function is used to read particle data created by
    // MAGI (https://bitbucket.org/ymiki/magi). The particle
    // data must be in the TIPSY format.

    // Read buffers
    MAGI_Tipsy_header header;
    std::vector<MAGI_Tipsy_particle> ptcl;

    // Read particle data
    std::ifstream input_file;
    input_file.open(input_file_name.c_str(), std::ios::in | std::ios::binary);
    if (input_file) {
        input_file.read((char *)&header, sizeof(MAGI_Tipsy_header));
        if (!is_little_endian) {
            std::cout << "nbodies = " << header.nbodies << std::endl;
            header.nbodies = reverseEndian(header.nbodies);
            std::cout << "nbodies = " << header.nbodies << std::endl;
        }
        ptcl.resize(header.nbodies);
        input_file.read((char *)&ptcl[0], sizeof(MAGI_Tipsy_particle)*header.nbodies);
    }
    else {
        std::cout << "cannot open the file " << input_file_name << std::endl;
        PS::Abort(-1);
        std::exit(1);
    }
    input_file.close();

    // Copy particle data
    psys.setNumberOfParticleLocal(header.nbodies);
    for (PS::S32 i=0; i<header.nbodies; i++) {
        if (is_little_endian) {
            psys[i].mass  = ptcl[i].mass;
            psys[i].pos.x = ptcl[i].pos[0];
            psys[i].pos.y = ptcl[i].pos[1];
            psys[i].pos.z = ptcl[i].pos[2];
            psys[i].vel.x = ptcl[i].vel[0];
            psys[i].vel.y = ptcl[i].vel[1];
            psys[i].vel.z = ptcl[i].vel[2];
        } else {
            psys[i].mass  = reverseEndian(ptcl[i].mass);
            psys[i].pos.x = reverseEndian(ptcl[i].pos[0]);
            psys[i].pos.y = reverseEndian(ptcl[i].pos[1]);
            psys[i].pos.z = reverseEndian(ptcl[i].pos[2]);
            psys[i].vel.x = reverseEndian(ptcl[i].vel[0]);
            psys[i].vel.y = reverseEndian(ptcl[i].vel[1]);
            psys[i].vel.z = reverseEndian(ptcl[i].vel[2]);
        }
    }

}

void GalaxyIC(PS::ParticleSystem<FP_dm>& psys_dm,
              PS::ParticleSystem<FP_star>& psys_star,
              PS::ParticleSystem<FP_gas>& psys_gas) {
    // Define the code units of MAGI
    //    [Important]
    //    (1) The values MUST BE consistent with "the computational units"
    //        written in the file doc/unit.txt, which is output by MAGI
    //        when we create a particle data with MAGI.
    //    (2) The MAGI's code units are DIFFERENT for unit systems
    //        a user choose. For detail, read Section "Unit systems in
    //        inc/constants.[c h]" in https://bitbucket.org/ymiki/magi.
    //    (3) In this sample code, "Galactic scale" unit is adopted.
    //        It is consistent with ./magi_data/cfg/Galaxy.cfg.
    const PS::F64 magi_unit_mass = 1.0e8 * phys_const::Msolar;
    const PS::F64 magi_unit_leng = phys_const::kpc;
    const PS::F64 magi_unit_time = 1.0e2 * phys_const::Myr;
    const PS::F64 magi_unit_velc = magi_unit_leng/magi_unit_time;
    // Initialize pseudorandom number generator
    PS::MTTS mt;
    mt.init_genrand(0);
    // Place Nbody particles
    std::string filename = "./magi_data/dat/Galaxy.tipsy";
    readTipsyFile(filename, psys_dm);
    const PS::S64 N_nbody = psys_dm.getNumberOfParticleLocal();
    for (PS::S32 i = 0; i < N_nbody; i++) {
        psys_dm[i].id = i;
        psys_dm[i].acc  = 0.0;
        psys_dm[i].pot  = 0.0;
    }
    // Place star particles [TODO: this should be merged to the above part]
    psys_star.setNumberOfParticleLocal(0); // tentative
    // Place SPH particles to form an exponential-disk
    const PS::S64 N_sph = (1<<20); // 2^{20}
    const PS::F64 M_gas = 1.0e10 * phys_const::Msolar;
    const PS::F64 Rs = 7.0 * phys_const::kpc; // scale radius
    const PS::F64 Rt = 12.5 * phys_const::kpc; // truncation radius
    const PS::F64 zd = 400.0 * phys_const::pc; // scale height
    const PS::F64 zt = 1.0 * phys_const::kpc; // truncation height
    const PS::F64 temp = 1.0e4; // gas temperature
    const PS::F64 mu = run_param::ism::mu; // mean molecular weight relative to the mass of hydrogen
    psys_gas.setNumberOfParticleLocal(N_sph);
    PS::S64 id = 0;
    for (PS::S32 i = 0; i < N_sph; i++) {
        // First make a uniform disk with a finite thickness
        PS::F64 x,y,z,R2;
        for(;;) {
            x = (2.0 * mt.genrand_res53() - 1.0) * Rt;
            y = (2.0 * mt.genrand_res53() - 1.0) * Rt;
            z = (2.0 * mt.genrand_res53() - 1.0) * zt;
            R2 = (x * x + y * y);
            if ((R2 < Rt * Rt) && (std::abs(z) < zt)) break;
        }
        const PS::F64 R = std::sqrt(R2);
        // Then re-scale particle position to generate an exponential disk
        PS::F64 R_low = 0.0, R_high = Rt, R_new;
        const PS::F64 eps = 1.0e-6;
        const PS::F64 val_trgt = (R/Rt)*(R/Rt);
        for(;;) {
            R_new = 0.5 * (R_low + R_high);
            const PS::F64 val = (1.0 - (R_new/Rs + 1.0)*std::exp(-R_new/Rs)) 
                              / (1.0 - (   Rt/Rs + 1.0)*std::exp(-   Rt/Rs));
            if (val < val_trgt) R_low  = R_new;
            if (val > val_trgt) R_high = R_new;
            const PS::F64 reldiff = 2.0 * std::abs(R_low - R_high)/(R_low + R_high);
            if (reldiff < eps) {
                R_new = 0.5 * (R_low + R_high);
                break;
            }
        }
        PS::F64 z_new;
        if (z >= 0.0) z_new = - zd * std::log(1.0 - (z/zt) * (1.0 - std::exp(-zt/zd)));
        else z_new = zd * std::log(1.0 + (z/zt) * (1.0 - std::exp(-zt/zd)));
        // Set  
        psys_gas[id].id = id + N_nbody;
        psys_gas[id].mass0 = M_gas / N_sph;
        psys_gas[id].mass = psys_gas[id].mass0;
        psys_gas[id].pos.x = (R_new / R) * x;
        psys_gas[id].pos.y = (R_new / R) * y;
        psys_gas[id].pos.z = z_new;
        psys_gas[id].vel       = 0.0;
        psys_gas[id].acc_grav  = 0.0;
        psys_gas[id].pot_grav  = 0.0;
        psys_gas[id].acc_hydro = 0.0;
        psys_gas[id].eng = (phys_const::kBoltz * temp)/((run_param::ism::gamma - 1.0) * mu * phys_const::Mproton);
        psys_gas[id].h = std::pow(Rt*Rt*zt, 1.0/3.0) * std::pow((PS::F64) run_param::sph::N_ngb / (PS::F64)N_sph, 1.0/3.0);
        psys_gas[id].alpha = run_param::sph::alpha_AV_ini;
        psys_gas[id].mabn[0]  = run_param::ism::Xhydrogen_solar;
        psys_gas[id].mabn[1]  = run_param::ism::Yhelium_solar;
        psys_gas[id].mabn[2]  = run_param::ism::Zcarbon_solar;
        psys_gas[id].mabn[3]  = run_param::ism::Znitrogen_solar;
        psys_gas[id].mabn[4]  = run_param::ism::Zoxygen_solar;
        psys_gas[id].mabn[5]  = run_param::ism::Zneon_solar;
        psys_gas[id].mabn[6]  = run_param::ism::Zmagnesium_solar;
        psys_gas[id].mabn[7]  = run_param::ism::Zsilicon_solar;
        psys_gas[id].mabn[8]  = run_param::ism::Zsulfur_solar;
        psys_gas[id].mabn[9]  = run_param::ism::Zcalcium_solar;
        psys_gas[id].mabn[10] = run_param::ism::Ziron_solar;
        psys_gas[id].mabn[11] = run_param::ism::Znickel_solar;
        psys_gas[id].mabn[12] = run_param::ism::Zeuropium_solar;
        // Note that entropy is determined in setEntropy().
        id++;
    }
    // Unit convertion (MAGI unit, CGS unit -> G=M=R=1 system)
    run_param::unit::mass = 0.0;
    run_param::unit::leng = 0.0;
    for (PS::S32 i = 0; i < N_nbody; i++) {
        psys_dm[i].mass *= magi_unit_mass;
        psys_dm[i].pos  *= magi_unit_leng;
        psys_dm[i].vel  *= magi_unit_velc;
        run_param::unit::mass += psys_dm[i].mass;
        const PS::F64 r = std::sqrt(psys_dm[i].pos * psys_dm[i].pos);
        if (r > run_param::unit::leng) run_param::unit::leng = r;
    }
    std::cout << "Total mass in N-body particles = "
              << run_param::unit::mass/phys_const::Msolar << " [Msolar]" << std::endl;
    for (PS::S32 i = 0; i< N_sph; i++) {
        run_param::unit::mass += psys_gas[i].mass;
        const PS::F64 r = std::sqrt(psys_gas[i].pos * psys_gas[i].pos);
        if (r > run_param::unit::leng) run_param::unit::leng = r;
    }
    run_param::unit::time = std::sqrt(CUBE(run_param::unit::leng)/(phys_const::Ggrav * run_param::unit::mass));
    run_param::unit::velc = run_param::unit::leng / run_param::unit::time;
    run_param::unit::eng  = run_param::unit::mass * SQ(run_param::unit::velc);
    run_param::unit::spen = run_param::unit::eng  / run_param::unit::mass;
    for (PS::S32 i = 0; i < N_nbody; i++) {
        psys_dm[i].mass /= run_param::unit::mass;
        psys_dm[i].pos  /= run_param::unit::leng;
        psys_dm[i].vel  /= run_param::unit::velc;
    }
    for (PS::S32 i = 0; i < N_sph; i++) {
        psys_gas[i].mass0 /= run_param::unit::mass;
        psys_gas[i].mass  /= run_param::unit::mass;
        psys_gas[i].pos   /= run_param::unit::leng;
        psys_gas[i].vel   /= run_param::unit::velc;
        psys_gas[i].h     /= run_param::unit::leng;
        psys_gas[i].eng   /= run_param::unit::spen;
    }
    // Set boundary condition
    run_param::basic::bc = PS::BOUNDARY_CONDITION_OPEN;
    // Set the parameters for gravity calculation
    run_param::grav::soft::theta = 0.5;
    run_param::grav::soft::n_group_limit = 64;
    run_param::grav::soft::n_leaf_limit = 8;
    run_param::grav::soft::eps_dm = 1.0e2 * phys_const::pc / run_param::unit::leng;
    run_param::grav::soft::eps_star = 3.0 * phys_const::pc / run_param::unit::leng;
    run_param::grav::soft::eps_gas = 3.0 * phys_const::pc / run_param::unit::leng;
    run_param::grav::soft::dt_fid = 1.0e5 * phys_const::yr / run_param::unit::time;
    run_param::grav::hard::eps = 0.01 * phys_const::pc / run_param::unit::leng;
    run_param::grav::hard::eta = 0.1;
    run_param::grav::hard::eta_s = 1.0e-3;
    run_param::grav::hard::dt_limit = run_param::grav::soft::dt_fid / 8;
    // Set the parameters for SPH calculation
    run_param::sph::n_group_limit = 64;
    run_param::sph::n_leaf_limit = 1;
    run_param::sph::n_jp_limit = 4096;
    // Set I/O intervals
    run_param::io::dt_dump = 0.001;
    run_param::io::time_dump = run_param::io::dt_dump;
    run_param::io::dt_dump_rst = 0.001;
    run_param::io::time_dump_rst = run_param::io::dt_dump;
    // Set the end time of the simulation
    run_param::basic::time_end = 10.0;

    if (PS::Comm::getRank() == 0) {
        std::cout << "An initial condition for isolated galaxy simulation is made." << std::endl;
        std::cout << "N_nbody = " << N_nbody << ", N_sph = " << N_sph << std::endl;
    }


#if 0
    // [for DEBUG]
    // Compute the surface gas density and output it.
    // Also check the distribution of particles w.r.t. the z coordinate.
    {
        // Set the resolution of bins, etc.
        const PS::S32 nbin = 64;
        const PS::F64 safety = 1.001;
        // Compute the maxium clyndrical radius of gas particles
        PS::F64 Rmax = 0.0;
        for (PS::S32 i = 0; i < N_sph; i++) {
            const PS::F64 R = std::sqrt(psys_gas[i].pos.x * psys_gas[i].pos.x
                                       +psys_gas[i].pos.y * psys_gas[i].pos.y);
            if (R > Rmax) Rmax = R;
        }
        // Compute the surface density
        const PS::F64 dR = (safety * Rmax)/nbin;
        PS::F64 Sigma[nbin] = {};
        for (PS::S32 i = 0; i < N_sph; i++) {
            const PS::F64 R = std::sqrt(psys_gas[i].pos.x * psys_gas[i].pos.x
                                       +psys_gas[i].pos.y * psys_gas[i].pos.y);
            const PS::S32 indx = R/dR;
            Sigma[indx] += psys_gas[i].mass;
        }
        for (PS::S32 n = 0; n < nbin; n++) 
            Sigma[n] /= 2.0 * math_const::pi * ((n+1)*(n+1)-n*n)*dR*dR;
        // Output the surface density
        std::string filename = "Sigma.txt";
        std::ofstream output_file;
        output_file.open(filename.c_str(),std::ios::trunc);
        output_file.setf(std::ios_base::scientific,
                         std::ios_base::floatfield);
        output_file << std::setprecision(15) << std::showpos;
        for (PS::S32 n = 0; n < nbin; n++) {
            output_file << (n + 0.5) * dR << " " << Sigma[n] << " " << std::endl;
        }
        output_file.close();
        // Compute the distribution function 
        const PS::F64 zmax = zt/run_param::unit::leng; 
        const PS::F64 dz = (safety * 2.0 * zmax)/nbin;
        PS::S32 dist_func[nbin] = {};
        for (PS::S32 i = 0; i < N_sph; i++) {
            const PS::S32 indx = (psys_gas[i].pos.z + safety * zmax)/dz;
            dist_func[indx]++;
        }
        // Output the distribution function
        filename = "dist_func.txt";
        output_file.open(filename.c_str(),std::ios::trunc);
        for (PS::S32 n = 0; n < nbin; n++) {
            const PS::F64 z = (n + 0.5) * dz - safety * zmax;
            output_file << z << " " << dist_func[n] << std::endl;
        }
        output_file.close();
    }
#endif
}

PS::F64ort getDomainInSphericalCoordinate(const PS::F64ort box_cart) {
    // Make a list of the coordinates of the vertex of box_cart
    PS::F64 x,y,z;
    std::vector<PS::F64vec> pos_list;
    x = box_cart.low_.x;  y = box_cart.low_.y;  z = box_cart.low_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
    x = box_cart.high_.x; y = box_cart.low_.y;  z = box_cart.low_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
    x = box_cart.low_.x;  y = box_cart.high_.y; z = box_cart.low_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
    x = box_cart.high_.x; y = box_cart.high_.y; z = box_cart.low_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
    x = box_cart.low_.x;  y = box_cart.low_.y;  z = box_cart.high_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
    x = box_cart.high_.x; y = box_cart.low_.y;  z = box_cart.high_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
    x = box_cart.low_.x;  y = box_cart.high_.y; z = box_cart.high_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
    x = box_cart.high_.x; y = box_cart.high_.y; z = box_cart.high_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
   
    // Calculate the ranges of the spherical coordinate corresponding to box_cart
    const PS::F64vec O(0.0); // origin
    const PS::F64 TWO_PI = 8.0 * std::atan(1.0); 
    const PS::F64 F64MAX = std::numeric_limits<PS::F64>::max();
    // Handle the case that box_cart contains the origin
    PS::F64 r_min = std::sqrt(box_cart.getDistanceMinSQ(O));
    PS::F64 r_max {0.0};
    // Handle the cases that box_cart crosses the z-axis
    PS::F64 mu_min { F64MAX};
    PS::F64 mu_max {-F64MAX}; 
    if ((box_cart.low_.x <= 0.0) && (0.0 <= box_cart.high_.x) &&
        (box_cart.low_.y <= 0.0) && (0.0 <= box_cart.high_.y)) {
        if (box_cart.low_.z <= 0.0)  mu_min = -1.0;
        if (0.0 <= box_cart.high_.z) mu_max =  1.0;
    }
    // Determine how to calculate the range of \phi
    bool use_0to2pi_range {false};
    PS::F64 phi_min { F64MAX};
    PS::F64 phi_max {-F64MAX};
    if (box_cart.low_.x <= 0.0) {
        // In this case, the \phi-range is given as a part of [0, 2\pi].
        use_0to2pi_range = true;
        if (0.0 <= box_cart.high_.x) {
            // In this case, box_cart contains the plane of x=0.
            phi_min = 0.0;
            phi_max = TWO_PI;
        }
    }
    // Shrink *_[min|max] using the coordinates of the vertex of box_cart
    for (PS::S32 i = 0; i < pos_list.size(); i++) {
        const PS::F64vec pos = pos_list[i];
        const PS::F64 r = std::sqrt(pos * pos);
        const PS::F64 mu = pos.z / r;
        PS::F64 phi = std::atan2(pos.y, pos.x);
        if (use_0to2pi_range && phi < 0.0) phi += TWO_PI;
        r_max = std::max(r_max, r);
        mu_min = std::min(mu_min, mu);
        mu_max = std::max(mu_max, mu);
        phi_min = std::min(phi_min, phi);
        phi_max = std::max(phi_max, phi);
    }

    PS::F64ort box_sph;
    box_sph.low_.x  = r_min;
    box_sph.high_.x = r_max;
    box_sph.low_.y  = mu_min;
    box_sph.high_.y = mu_max;
    box_sph.low_.z  = phi_min;
    box_sph.high_.z = phi_max;
    return box_sph;
}

PS::F64ort getDomainInCylindricalCoordinate(const PS::F64ort box_cart) {
    // Make a list of the coordinates of the vertex of box_cart
    PS::F64 x,y,z;
    std::vector<PS::F64vec> pos_list;
    x = box_cart.low_.x;  y = box_cart.low_.y;  z = box_cart.low_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
    x = box_cart.high_.x; y = box_cart.low_.y;  z = box_cart.low_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
    x = box_cart.low_.x;  y = box_cart.high_.y; z = box_cart.low_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
    x = box_cart.high_.x; y = box_cart.high_.y; z = box_cart.low_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
    x = box_cart.low_.x;  y = box_cart.low_.y;  z = box_cart.high_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
    x = box_cart.high_.x; y = box_cart.low_.y;  z = box_cart.high_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
    x = box_cart.low_.x;  y = box_cart.high_.y; z = box_cart.high_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
    x = box_cart.high_.x; y = box_cart.high_.y; z = box_cart.high_.z;
    pos_list.push_back(PS::F64vec(x,y,z));
   
    // Calculate the ranges of the cylindrical coordinate corresponding to box_cart
    const PS::F64vec O(0.0); // origin
    const PS::F64 TWO_PI = 8.0 * std::atan(1.0);
    const PS::F64 F64MAX = std::numeric_limits<PS::F64>::max();
    // Handle the case that box_cart contains the z-axis
    PS::Orthotope2<PS::F64> box_cart_xy;
    box_cart_xy.low_.x = box_cart.low_.x;
    box_cart_xy.low_.y = box_cart.low_.y;
    box_cart_xy.high_.x = box_cart.high_.x;
    box_cart_xy.high_.y = box_cart.high_.y;
    PS::F64 R_min = std::sqrt(box_cart_xy.getDistanceMinSQ(PS::Vector2<PS::F64>(0.0)));
    PS::F64 R_max {0.0};
    // Determine how to calculate the range of \phi
    bool use_0to2pi_range {false};
    PS::F64 phi_min { F64MAX};
    PS::F64 phi_max {-F64MAX};
    if (box_cart.low_.x <= 0.0) {
        // In this case, the \phi-range is given as a part of [0, 2\pi].
        use_0to2pi_range = true;
        if (0.0 <= box_cart.high_.x) {
            // In this case, box_cart contains the plane of x=0.
            phi_min = 0.0;
            phi_max = TWO_PI;
        }
    }
    PS::F64 z_min { F64MAX};
    PS::F64 z_max {-F64MAX};
    // Shrink _[min|max] using the coordinates of the vertex of box_cart
    for (PS::S32 i = 0; i < pos_list.size(); i++) {
        const PS::F64vec pos = pos_list[i];
        const PS::F64 R = std::sqrt(SQ(pos.x) + SQ(pos.y));
        PS::F64 phi = std::atan2(pos.y, pos.x);
        if (use_0to2pi_range && phi < 0.0) phi += TWO_PI;
        const PS::F64 z = pos.z;
        R_max   = std::max(R_max, R);
        phi_min = std::min(phi_min, phi);
        phi_max = std::max(phi_max, phi);
        z_min   = std::min(z_min, z);
        z_max   = std::max(z_max, z);
    }

    PS::F64ort box_cyl;
    box_cyl.low_.x  = R_min;
    box_cyl.high_.x = R_max;
    box_cyl.low_.y  = phi_min;
    box_cyl.high_.y = phi_max;
    box_cyl.low_.z  = z_min;
    box_cyl.high_.z = z_max;
    return box_cyl;
}

void defineLogarithmicGrid(std::vector<PS::F64> & r,
                           std::vector<PS::F64> & dr,
                           const PS::F64 rb_min,
                           const PS::F64 rb_max,
                           const PS::S32 N) {
    assert(rb_min > 0.0 && rb_max > 0.0 && N > 0);
    r.resize(N);
    dr.resize(N);
    const PS::F64 log_rb_min = std::log(rb_min);
    const PS::F64 log_rb_max = std::log(rb_max);
    const PS::F64 dlog_rb = (log_rb_max - log_rb_min) / N;
    for (PS::S32 i = 0; i < N; i++) {
        const PS::F64 log_rb_m = log_rb_min + dlog_rb * i;
        const PS::F64 log_rb_p = log_rb_min + dlog_rb * (i+1);
        const PS::F64 rb_m = std::exp(log_rb_m);
        const PS::F64 rb_p = std::exp(log_rb_p);
        dr[i] = rb_p - rb_m;
        r[i] = 0.5 * (rb_m + rb_p);
    }
}

void defineLogarithmicGrid(std::vector<PS::F64> & x,
                           std::vector<PS::F64> & dx,
                           const PS::F64 xb_lo, // the left boundary of domain
                           const PS::F64 xb_hi, // the right boundary of domain
                           const PS::F64 dx_min, // the minimum resolution
                           const PS::F64 x_org, // the location where the highest resolution is acheived
                           const PS::S32 N_inp) // A rough number of grid points
{
    if ((xb_lo - x_org)*(xb_hi - x_org) < 0.0) {
        // In this case, the specified domain contains x_org.
        // We use two logarithmic grids to cover [xb_lo, xb_hi].
        const PS::F64 D_max = std::max(std::fabs(xb_lo - x_org),
                                       std::fabs(xb_hi - x_org));
        const PS::S32 N_half = N_inp/2;
        std::vector<PS::F64> x_half, dx_half;
        // Calculate a logarithmic grid
        defineLogarithmicGrid(x_half, dx_half, dx_min, D_max, N_half);
        // Calculate the number of grids points required to cover [x_org, xb_hi]
        PS::S32 N_half_R {-1};
        for (PS::S32 i = 0; i < N_half; i++) {
            const PS::F64 xb = x_half[i] + 0.5 * dx_half[i];
            if (xb >= (xb_hi - x_org)) {
                N_half_R = i + 1;
                break;
            }
        }
        if (N_half_R == -1) N_half_R = N_half;
        // Calculate the number of grid points required to cover [xb_lo, x_org]
        PS::S32 N_half_L {-1};
        for (PS::S32 i = 0; i < N_half; i++) {
            const PS::F64 xb = x_half[i] + 0.5 * dx_half[i];
            if (xb >= (x_org - xb_lo)) {
                N_half_L = i + 1;
                break;
            }
        }
        if (N_half_L == -1) N_half_L = N_half;
        // Calculate x[], dx[]
        const PS::S32 N = N_half_L + 1 + N_half_R;
        x.resize(N);
        dx.resize(N);
        for (PS::S32 i = 0; i < N_half_L; i++) {
            x[N_half_L - 1 - i] = x_org - x_half[i];
            dx[N_half_L - 1 - i] = dx_half[i];
        }
        for (PS::S32 i = 0; i < N_half_R; i++) {
            x[i + N_half_L + 1] = x_org + x_half[i];
            dx[i + N_half_L + 1] = dx_half[i];
        }
        x[N_half_L] = x_org;
        dx[N_half_L] = (x[N_half_L+1] - 0.5 * dx[N_half_L+1])
                     - (x[N_half_L-1] + 0.5 * dx[N_half_L-1]);
        // Adjust the left- and right-most cells
        if (x[0] - 0.5 * dx[0] < xb_lo) {
            dx[0] = x[0] + 0.5 * dx[0]  - xb_lo;
            x[0] = xb_lo + 0.5 * dx[0];
        }
        if (x[N-1] + 0.5 * dx[N-1] > xb_hi) {
            dx[N-1] = xb_hi - (x[N-1] - 0.5 * dx[N-1]);
            x[N-1] = xb_hi - 0.5 * dx[N-1];
        }
    } else {
        // In this case, local domain does not contain x_org.
        // We use a single logarithmic grid to cover [xb_lo, xb_hi].
        const PS::F64 D_max = std::fabs(xb_hi - xb_lo);
        std::vector<PS::F64> x_tmp, dx_tmp;
        defineLogarithmicGrid(x_tmp, dx_tmp, dx_min, D_max, N_inp);
        // Calculate x[], dx[]
        const PS::S32 N = N_inp + 1;
        x.resize(N);
        dx.resize(N);
        if (xb_lo >= x_org) {
            // In this case, higher resolution is achieved at the side of xb_lo 
            for (PS::S32 i = 0; i < N_inp; i++) {
                x[i + 1]  = xb_lo + x_tmp[i];
                dx[i + 1] = dx_tmp[i];
            }
            // The leftmost cell
            const PS::F64 xb_2nd = x[1] - 0.5 * dx[1];
            dx[0] = xb_2nd - xb_lo;
            x[0] = 0.5 * (xb_lo + xb_2nd);
            // The rightmost cell
            if (x[N-1] + 0.5*dx[N-1] > xb_hi) {
                dx[N-1] = xb_hi - (x[N-1] - 0.5 * dx[N-1]);
                x[N-1] = xb_hi - 0.5 * dx[N-1];
            }
        } else {
            // In this case, higher resolution is achieved at the side of xb_hi
            for (PS::S32 i = 0; i < N_inp; i++) {
                x[N - 2 - i] = xb_hi - x_tmp[i];
                dx[N - 2 - i] = dx_tmp[i];
            }
            // the rightmost cell
            const PS::F64 xb_2nd = x[N-2] + 0.5 * dx[N-2];
            dx[N-1] = xb_hi - xb_2nd;
            x[N-1] = 0.5 * (xb_2nd + xb_hi);
            // the leftmost cell
            if (x[0] - 0.5 * dx[0] < xb_lo) {
                dx[0] = x[0] + 0.5 * dx[0] - xb_lo;
                x[0] = xb_lo + 0.5 * dx[0];
            }
        }
    }
}


void GordonBellIC(PS::ParticleSystem<FP_dm> & psys_dm,
                  PS::ParticleSystem<FP_star> & psys_star,
                  PS::ParticleSystem<FP_gas> & psys_gas,
                  PS::DomainInfo & dinfo) {
    // Local variables
    const PS::CommInfo comm_info = PS::Comm::getCommInfo();
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    const PS::S32 my_rank = PS::Comm::getRank();
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    PS::F64 time_start;
    // Set parameters for the Monte-Carlo integration
    constexpr PS::S32 N_smpl_glb = 100000000;
    // Define the parameters of the Hernquist model
    // (used to make the dark matter distribution)
    const PS::S64 N_dm = (1<<20);
    const PS::F64 M_dm = 1.0e12 * phys_const::Msolar;
    const PS::F64 rs_dm = 21.5 * phys_const::kpc; // scale radius
    const PS::F64 rt_dm = 1.0e3 * phys_const::kpc; // truncation radius
    HernquistModel model_dm(M_dm, rs_dm, rt_dm);
    const PS::F64 M_dm_cut = model_dm.getEnclosedMass(rt_dm);
    const PS::F64 M_dm_ptcl = M_dm_cut / N_dm;
    const PS::F64ort pos_domain_dm = PS::F64ort(PS::F64vec(-rt_dm, -rt_dm, -rt_dm),
                                                PS::F64vec( rt_dm,  rt_dm,  rt_dm));
    // Define the parameters of the Miyamoto-Nagai model
    // (used to make the star distribution)
    const PS::S64 N_star = (1<<20);
    const PS::F64 M_star = 4.0e10 * phys_const::Msolar;
    const PS::F64 Rs_star = 3.5 * phys_const::kpc; // scale radius
    const PS::F64 Rt_star = 15.0 * phys_const::kpc; // truncation radius
    const PS::F64 zs_star = 400.0 * phys_const::pc; // scale height
    const PS::F64 zt_star = 2.0 * phys_const::kpc; // truncation height
    MiyamotoNagaiModel model_star(M_star, Rs_star, zs_star, Rt_star, zt_star);
    const PS::F64 M_star_cut = math_util::qmc3d<PS::F64>(
        [model_star](const PS::F64 x, const PS::F64 y, const PS::F64 z)
            -> PS::F64 {return model_star.getDensityCutoffedXYZ(x,y,z);},
        -Rt_star, Rt_star,
        -Rt_star, Rt_star,
        -zt_star, zt_star,
        N_smpl_glb,
        run_param::prng::mt,
        comm_info);
    const PS::F64 M_star_ptcl = M_star_cut / N_star;
    const PS::F64ort pos_domain_star = PS::F64ort(PS::F64vec(-Rt_star, -Rt_star, -zt_star),
                                                  PS::F64vec( Rt_star,  Rt_star,  zt_star));
    // Define the parameters of the exponential disk model
    // (used to make the gas distribution)
    const PS::S64 N_gas = (1<<20); // 2^{20}
    const PS::F64 M_gas = 5.0e9 * phys_const::Msolar;
    const PS::F64 Rs_gas = 7.0 * phys_const::kpc; // scale radius
    const PS::F64 Rt_gas = 20.0 * phys_const::kpc; // truncation radius
    const PS::F64 zs_gas = 400.0 * phys_const::pc; // scale height
    const PS::F64 zt_gas = 1.0 * phys_const::kpc; // truncation height
    const PS::F64 temp = 1.0e4; // gas temperature
    const PS::F64 mu = run_param::ism::mu; // mean molecular weight relative to the proton mass
    ExponentialDiskModel model_gas(M_gas, Rs_gas, zs_gas, Rt_gas, zt_gas);
    const PS::F64 M_gas_cut = model_gas.getMassCutoffed();
    const PS::F64 M_gas_ptcl = M_gas_cut / N_gas;
    const PS::F64ort pos_domain_gas = PS::F64ort(PS::F64vec(-Rt_gas, -Rt_gas, -zt_gas),
                                                 PS::F64vec( Rt_gas,  Rt_gas,  zt_gas));
    // Get the number of subdomains in each direction and the coordinate of the local domain
    PS::S32 n_domain[3];
    PS::S32 id_domain[3];
    for (PS::S32 i = 0; i < 3; i++) {
        n_domain[i] = dinfo.getNDomain(i);
        id_domain[i] = dinfo.getRank1D(PS::Comm::getRank(), i);
        if (my_rank == 0) {
            std::cout << "n_domain[" << i << "] = " << n_domain[i] << std::endl; 
        }
    }
    // Create instances of CommInfo
    PS::CommInfo comm_info_x = PS::Comm::getCommInfo();
    std::vector<PS::S32> ranks_y, ranks_z;
    for(PS::S32 i = 0; i < comm_info_x.getNumberOfProc(); i++){
        if (id_domain[0] == dinfo.getRank1D(i,0)) { // having the same x
            ranks_y.push_back(i);
            if (id_domain[1] == dinfo.getRank1D(i,1)) { // having the same (x,y)
                ranks_z.push_back(i);
            }
        }
    }                                                     
    PS::CommInfo comm_info_y = comm_info_x.create(ranks_y.size(), &ranks_y[0]); 
    PS::CommInfo comm_info_z = comm_info_x.create(ranks_z.size(), &ranks_z[0]);
    if (my_rank == 0) {
        std::cout << "Instances of PS::CommInfo are created." << std::endl;
    }
    //PS::Finalize(); std::exit(0);
    // Determine the boundaries of the subdomains
    constexpr PS::F64 eps = 5.0e-3;
    constexpr PS::S32 itermax = 64;
    std::vector<PS::F64> xdiv, ydiv, zdiv;
    xdiv.resize(n_domain[0] + 1);
    ydiv.resize(n_domain[1] + 1);
    zdiv.resize(n_domain[2] + 1);
    std::vector<PS::F64> x_dm, dx_dm, y_dm, dy_dm, z_dm, dz_dm;
    std::vector<PS::F64> x_star, dx_star, y_star, dy_star, z_star, dz_star;
    std::vector<PS::F64> x_gas, dx_gas, y_gas, dy_gas, z_gas, dz_gas;
    std::vector<PS::F64> dN, N_cum;
    // (1) decompose domain along the X-direction
    if (my_rank == 0) std::cout << "Decomposing in X..." << std::endl;
    const PS::S32 n_proc_x = comm_info_x.getNumberOfProc();
    const PS::S32 my_rank_x = comm_info_x.getRank();
    xdiv[0]           = - rt_dm;
    xdiv[n_domain[0]] =   rt_dm;
    time_start = PS::GetWtime();
    if (n_domain[0] > 1) {
#if 1
        // Make grids
        // [FIXME] Current version of the implementation makes grids
        //         assuming that the system is not rotated and translated.
        constexpr PS::S32 nx = 1024;
        constexpr PS::S32 ny = 256;
        constexpr PS::S32 nz = ny;
        const PS::F64 rs_min = std::min({rs_dm, Rs_star, zs_star, Rs_gas, zs_gas});
        const PS::F64 rb_max = std::max({rt_dm, Rt_star, zt_star, Rt_gas, zt_gas});
        // (1) Make x grids
        //     (They must be the same because we're now trying to divide in x)
        defineLogarithmicGrid(x_dm, dx_dm, -rb_max, rb_max, 0.01 * rs_min, 0.0, nx);
        x_star.resize(x_dm.size());
        dx_star.resize(x_dm.size());
        x_gas.resize(x_dm.size());
        dx_gas.resize(x_dm.size());
        for (PS::S32 i = 0; i < x_dm.size(); i++) {
            x_star[i] = x_dm[i];
            dx_star[i] = dx_dm[i];
            x_gas[i] = x_dm[i];
            dx_gas[i] = dx_dm[i];
        }
        // (2) Make y & z grids
        //     (They CAN BE DIFFERENT for the type of particles)
        // (2-1) y & z grids for DM
        defineLogarithmicGrid(y_dm, dy_dm, -rb_max, rb_max, 0.01 * rs_min, 0.0, ny);
        z_dm.resize(y_dm.size());
        dz_dm.resize(y_dm.size());
        for (PS::S32 i = 0; i < y_dm.size(); i++) {
            z_dm[i] = y_dm[i];
            dz_dm[i] = dy_dm[i];
        }
        // (2-2) y & z grids for star
        const PS::F64 rs_min_star = std::min({Rs_star, zs_star});
        const PS::F64 rb_max_star = std::max({Rt_star, zt_star});
        defineLogarithmicGrid(y_star, dy_star, -rb_max_star, rb_max_star, 0.01 * rs_min_star, 0.0, ny);
        z_star.resize(y_star.size());
        dz_star.resize(y_star.size());
        for (PS::S32 i = 0; i < y_star.size(); i++) {
            z_star[i] = y_star[i];
            dz_star[i] = dy_star[i];
        }
        // (2-3) y & z grids for gas
        const PS::F64 rs_min_gas = std::min({Rs_gas, zs_gas});
        const PS::F64 rb_max_gas = std::max({Rt_gas, zt_gas});
        defineLogarithmicGrid(y_gas, dy_gas, -rb_max_gas, rb_max_gas, 0.01 * rs_min_gas, 0.0, ny);
        z_gas.resize(y_gas.size());
        dz_gas.resize(y_gas.size());
        for (PS::S32 i = 0; i < y_gas.size(); i++) {
            z_gas[i] = y_gas[i];
            dz_gas[i] = dy_gas[i];
        }
        // Calculate dN(x) in parallel
        dN.resize(x_dm.size());
        const PS::S32 rem = x_dm.size() % n_proc_x; // remainder
        const PS::S32 quo = x_dm.size() / n_proc_x; // quotient
        const PS::S32 sendcount = (my_rank_x < rem) ? (quo + 1) : quo;
        std::vector<PS::S32> recvcounts(n_proc_x), recvdispls(n_proc_x);
        comm_info_x.allGather(&sendcount, 1, &recvcounts[0]);
        recvdispls[0]=0;
        for (PS::S32 i = 1; i < n_proc_x; i++)
            recvdispls[i] = recvdispls[i-1] + recvcounts[i-1];
        const PS::S32 ista = recvdispls[my_rank_x];
        const PS::S32 iend = ista + sendcount - 1;
        for (PS::S32 i = ista; i <= iend; i++) {
            dN[i] = 0.0;
            const PS::F64 x = x_dm[i];
            const PS::F64 dx = dx_dm[i]; 
            // contribution from DM
            for (PS::S32 k = 0; k < z_dm.size(); k++) {
                const PS::F64 z = z_dm[k];
                const PS::F64 dz = dz_dm[k];
                for (PS::S32 j = 0; j < y_dm.size(); j++) {
                    const PS::F64 y = y_dm[j];
                    const PS::F64 dy = dy_dm[j];
                    PS::F64 r = std::sqrt(x*x + y*y + z*z);
                    PS::F64 rho;
                    if (r > 0.0) {
                        rho = model_dm.getDensityCutoffedXYZ(x,y,z);
                    } else {
                        // Handling the divergence at r = 0
                        const PS::F64 r_min_tmp = std::min(std::fabs(dx), std::min(std::fabs(dy), std::fabs(dz)));
                        const PS::F64 r_max_tmp = 0.5 * std::sqrt(dx*dx+dy*dy+dz*dz);
                        const PS::F64 rho_min = model_dm.getEnclosedMass(r_min_tmp)
                                              / (4.0 * math_const::pi * CUBE(r_min_tmp) / 3.0);
                        const PS::F64 rho_max = model_dm.getEnclosedMass(r_max_tmp)
                                              / (4.0 * math_const::pi * CUBE(r_max_tmp) / 3.0);
                        rho = 0.5 * (rho_min + rho_max);
                    }
                    dN[i] += rho * (dx * dy * dz) / M_dm_ptcl;
                    //dN[i] += rho * (dx * dy * dz); // for debug
                }
            }
            // contribution from star
            for (PS::S32 k = 0; k < z_star.size(); k++) {
                const PS::F64 z = z_star[k];
                const PS::F64 dz = dz_star[k];
                for (PS::S32 j = 0; j < y_star.size(); j++) {
                    const PS::F64 y = y_star[j];
                    const PS::F64 dy = dy_star[j];
                    const PS::F64 rho = model_star.getDensityCutoffedXYZ(x,y,z);
                    dN[i] += rho * (dx * dy * dz) / M_star_ptcl;
                    //dN[i] += rho * (dx * dy * dz); // for debug
                }
            }
            // contribution from gas
            for (PS::S32 k = 0; k < z_gas.size(); k++) {
                const PS::F64 z = z_gas[k];
                const PS::F64 dz = dz_gas[k];
                for (PS::S32 j = 0; j < y_gas.size(); j++) {
                    const PS::F64 y = y_gas[j];
                    const PS::F64 dy = dy_gas[j];
                    const PS::F64 rho = model_gas.getDensityCutoffedXYZ(x,y,z);
                    dN[i] += rho * (dx * dy * dz) / M_gas_ptcl;
                    //dN[i] += rho * (dx * dy * dz); // for debug
                }
            }
        }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                       &dN[0], &recvcounts[0], &recvdispls[0], PS::GetDataType<PS::F64>(),
                       comm_info_x.getCommunicator());
#endif
        // Calculate the cumulative mass from x = - xcut.
        N_cum.resize(x_dm.size());
        N_cum[0] = dN[0];
        for (PS::S32 i = 1; i < N_cum.size(); i++) {
            N_cum[i] = N_cum[i-1] + dN[i-1];
        }
        if (my_rank == 0) {
            std::cout << "N_cum      = " << N_cum.back() << std::endl;
            std::cout << "N_dm_cut   = " << M_dm_cut / M_dm_ptcl << std::endl;
            std::cout << "N_star_cut = " << M_star_cut / M_star_ptcl << std::endl;
            std::cout << "N_gas_cut  = " << M_gas_cut / M_gas_ptcl << std::endl;
            // for debug
            //std::cout << "M_cum      = " << N_cum.back() / phys_const::Msolar << std::endl;
            //std::cout << "M_dm_cut   = " << M_dm_cut / phys_const::Msolar << std::endl;
            //std::cout << "M_star_cut = " << M_star_cut / phys_const::Msolar << std::endl;
            //std::cout << "M_gas_cut  = " << M_gas_cut / phys_const::Msolar << std::endl;
        }
        // Find dividing points in x
        for (PS::S32 i = 0; i < n_domain[0] - 1; i++) {
            const PS::F64 N_trgt = (N_dm + N_star + N_gas) * (i+1) 
                                 / (PS::F64) n_domain[0];
            PS::S32 idx = -1;
            for (PS::S32 k = 0; k < x_dm.size(); k++) {
                // for debug
                //if (my_rank == 0) {
                //    std::cout << "i = " << i
                //              << " k = " << k
                //              << " N_cum[k] = " << N_cum[k]
                //              << " N_trgt = " << N_trgt
                //              << std::endl;
                //}
                if (N_cum[k] > N_trgt) {
                    idx = k;
                    break;
                }
            }
            assert(idx >= 0);
            PS::F64 x;
            if (idx > 0) {
                const PS::F64 slope = (N_cum[idx] - N_cum[idx-1])/dx_dm[idx];
                x = (x_dm[idx] - 0.5*dx_dm[idx]) + (N_trgt - N_cum[idx-1]) / slope;
            } else {
                const PS::F64 slope = N_cum[idx]/dx_dm[idx];
                x = (x_dm[idx] - 0.5*dx_dm[idx]) + N_trgt / slope;
            }
            xdiv[i+1] = x;
        }
#else
        const PS::S32 rem = N_smpl_glb % n_proc_x; // remainder
        const PS::S32 quo = N_smpl_glb / n_proc_x; // quotient
        const PS::S32 N_smpl_loc = (my_rank_x < rem) ? (quo + 1) : quo;
        for (PS::S32 i = 0; i < n_domain[0] - 1; i++) {
            // Calculate the x coordinate of the boundary i+1/2.
            const PS::F64 N_trgt = N_dm * (i+1) / (PS::F64) n_domain[0];
            PS::F64 x_min = xdiv[0];
            PS::F64 x_lo = xdiv[i];
            PS::F64 x_hi = xdiv[n_domain[0]];
            PS::F64 x_bisec = 0.5 * (x_lo + x_hi);
            PS::S32 iter = 0;
            for (;;) {
                // Estimate the number of particles in [x_min, x_bisec]
                const PS::F64 dv = SQ(2.0 * rt_dm) * (x_bisec - x_min) / N_smpl_glb;
                PS::F64 N_eval {0.0}; 
                for (PS::S32 k = 0; k < N_smpl_loc; k++) {
                    const PS::F64 x = x_min + (x_bisec - x_min) * dist(run_param::prng::mt);
                    const PS::F64 y = rt_dm * (2.0 * dist(run_param::prng::mt) - 1.0);
                    const PS::F64 z = rt_dm * (2.0 * dist(run_param::prng::mt) - 1.0);
                    N_eval += model_dm.getDensityCutoffedXYZ(x,y,z) / M_dm_ptcl;
                }
                N_eval = comm_info_x.getSum(N_eval) * dv;
                // Check if the termination condition is satisfied or not
                const PS::F64 relerr = std::fabs((N_eval - N_trgt)/N_trgt);
                if (relerr <= eps) break;
                // Calculate the next x_bisec and update h_L & h_U
                if ((N_eval < N_trgt) && (x_lo < x_bisec)) x_lo = x_bisec;
                if ((N_eval > N_trgt) && (x_bisec < x_hi)) x_hi = x_bisec;
                x_bisec = 0.5*(x_lo + x_hi);
                iter++;
                // Check an error
                if (iter >= itermax) {
                    if (my_rank == 0) {
                        std::cout << "Too many iterations in the decomposition along X!" << std::endl;
                    }
                    PS::Abort(-1); std::exit(-1);
                }
            }
            xdiv[i+1] = x_bisec;
        }
#endif
    }
    PS::Comm::barrier();
    if (my_rank == 0) std::cout << "took " << PS::GetWtime() - time_start << "[s]." << std::endl;
    // (2) decompose domain along the Y-direction
    if (my_rank == 0) std::cout << "Decomposing in Y..." << std::endl;
    const PS::S32 n_proc_y = comm_info_y.getNumberOfProc();
    const PS::S32 my_rank_y = comm_info_y.getRank();
    ydiv[0]           = - rt_dm;
    ydiv[n_domain[1]] =   rt_dm;
    time_start = PS::GetWtime();
    if (n_domain[1] > 1) {
#if 1
        // Make grids
        constexpr PS::S32 ny = 512;
        constexpr PS::S32 nx = 256;
        constexpr PS::S32 nz = nx;
        const PS::F64 rs_min = std::min({rs_dm, Rs_star, zs_star, Rs_gas, zs_gas});
        const PS::F64 rb_max = std::max({rt_dm, Rt_star, zt_star, Rt_gas, zt_gas});
        // (1) Make y grids
        //     (They must be the same because we're now trying to divide in y)
        defineLogarithmicGrid(y_dm, dy_dm, -rb_max, rb_max, 0.01 * rs_min, 0.0, ny);
        y_star.resize(y_dm.size());
        dy_star.resize(y_dm.size());
        y_gas.resize(y_dm.size());
        dy_gas.resize(y_dm.size());
        for (PS::S32 i = 0; i < y_dm.size(); i++) {
            y_star[i] = y_dm[i];
            dy_star[i] = dy_dm[i];
            y_gas[i] = y_dm[i];
            dy_gas[i] = dy_dm[i];
        }
        // (2) Make x & z grids
        //     (They CAN BE DIFFERENT for the type of particles)
        const PS::F64 xb_lo = xdiv[id_domain[0]];
        const PS::F64 xb_hi = xdiv[id_domain[0] + 1];
        // (2-1) x & z grids for DM
        defineLogarithmicGrid(x_dm, dx_dm, xb_lo, xb_hi, 0.01 * rs_min, 0.0, nx);
        defineLogarithmicGrid(z_dm, dz_dm, -rb_max, rb_max, 0.01 * rs_min, 0.0, nz);
        // (2-2) x & z grids for star
        const PS::F64 rs_min_star = std::min({Rs_star, zs_star});
        const PS::F64 rb_max_star = std::max({Rt_star, zt_star});
        const PS::F64 xb_lo_star = (xb_lo <= -rb_max_star && -rb_max_star <= xb_hi) ?  -rb_max_star : xb_lo;
        const PS::F64 xb_hi_star = (xb_lo <=  rb_max_star &&  rb_max_star <= xb_hi) ?   rb_max_star : xb_hi;
        defineLogarithmicGrid(x_star, dx_star, xb_lo_star, xb_hi_star, 0.01 * rs_min_star, 0.0, nx);
        defineLogarithmicGrid(z_star, dz_star, -rb_max_star, rb_max_star, 0.01 * rs_min_star, 0.0, nz);
        // (2-3) x & z grids for gas
        const PS::F64 rs_min_gas = std::min({Rs_gas, zs_gas});
        const PS::F64 rb_max_gas = std::max({Rt_gas, zt_gas});
        const PS::F64 xb_lo_gas = (xb_lo <= -rb_max_gas && -rb_max_gas <= xb_hi) ?  -rb_max_gas : xb_lo;
        const PS::F64 xb_hi_gas = (xb_lo <=  rb_max_gas &&  rb_max_gas <= xb_hi) ?   rb_max_gas : xb_hi;
        defineLogarithmicGrid(x_gas, dx_gas, xb_lo_gas, xb_hi_gas, 0.01 * rs_min_gas, 0.0, nx);
        defineLogarithmicGrid(z_gas, dz_gas, -rb_max_gas, rb_max_gas, 0.01 * rs_min_gas, 0.0, nz);
#if 0
        if (my_rank == 0) {
            std::ofstream ofs;
            ofs.open("x_dm.txt", std::ios::trunc);
            for (PS::S32 i = 0; i < x_dm.size(); i++) {
                ofs << x_dm[i] << "   " << dx_dm[i] << std::endl;
            }
            ofs.close();
            ofs.open("x_star.txt", std::ios::trunc);
            for (PS::S32 i = 0; i < x_star.size(); i++) {
                ofs << x_star[i] << "   " << dx_star[i] << std::endl;
            }
            ofs.close();
            ofs.open("x_gas.txt", std::ios::trunc);
            for (PS::S32 i = 0; i < x_gas.size(); i++) {
                ofs << x_gas[i] << "   " << dx_gas[i] << std::endl;
            }
            ofs.close();
        }
#endif
        // Calculate dN(y) in parallel
        dN.resize(y_dm.size());
        const PS::S32 rem = y_dm.size() % n_proc_y; // remainder
        const PS::S32 quo = y_dm.size() / n_proc_y; // quotient
        const PS::S32 sendcount = (my_rank_y < rem) ? (quo + 1) : quo;
        std::vector<PS::S32> recvcounts(n_proc_y), recvdispls(n_proc_y);
        comm_info_y.allGather(&sendcount, 1, &recvcounts[0]);
        recvdispls[0]=0;
        for (PS::S32 i = 1; i < n_proc_y; i++)
            recvdispls[i] = recvdispls[i-1] + recvcounts[i-1];
        const PS::S32 jsta = recvdispls[my_rank_y];
        const PS::S32 jend = jsta + sendcount - 1;
        for (PS::S32 j = jsta; j <= jend; j++) {
            dN[j] = 0.0;
            const PS::F64 y = y_dm[j];
            const PS::F64 dy = dy_dm[j];
            // contribution from DM
            for (PS::S32 k = 0; k < z_dm.size(); k++) {
                const PS::F64 z = z_dm[k];
                const PS::F64 dz = dz_dm[k];
                for (PS::S32 i = 0; i < x_dm.size(); i++) {
                    const PS::F64 x = x_dm[i];
                    const PS::F64 dx = dx_dm[i];
                    PS::F64 r = std::sqrt(x*x + y*y + z*z);
                    PS::F64 rho;
                    if (r > 0.0) {
                        rho = model_dm.getDensityCutoffedXYZ(x,y,z);
                    } else {
                        // Handling the divergence at r = 0
                        const PS::F64 r_min_tmp = std::min(std::fabs(dx), std::min(std::fabs(dy), std::fabs(dz)));
                        const PS::F64 r_max_tmp = 0.5 * std::sqrt(dx*dx+dy*dy+dz*dz);
                        const PS::F64 rho_min = model_dm.getEnclosedMass(r_min_tmp)
                                              / (4.0 * math_const::pi * CUBE(r_min_tmp) / 3.0);
                        const PS::F64 rho_max = model_dm.getEnclosedMass(r_max_tmp)
                                              / (4.0 * math_const::pi * CUBE(r_max_tmp) / 3.0);
                        rho = 0.5 * (rho_min + rho_max);
                    }
                    dN[j] += rho * (dx * dy * dz) / M_dm_ptcl;
                }
            }
            // contribution from star
            for (PS::S32 k = 0; k < z_star.size(); k++) {
                const PS::F64 z = z_star[k];
                const PS::F64 dz = dz_star[k];
                for (PS::S32 i = 0; i < x_star.size(); i++) {
                    const PS::F64 x = x_star[i];
                    const PS::F64 dx = dx_star[i];
                    const PS::F64 rho = model_star.getDensityCutoffedXYZ(x,y,z);
                    dN[j] += rho * (dx * dy * dz) / M_star_ptcl;
                }
            }
            // contribution from gas
            for (PS::S32 k = 0; k < z_gas.size(); k++) {
                const PS::F64 z = z_gas[k];
                const PS::F64 dz = dz_gas[k];
                for (PS::S32 i = 0; i < x_gas.size(); i++) {
                    const PS::F64 x = x_gas[i];
                    const PS::F64 dx = dx_gas[i];
                    const PS::F64 rho = model_gas.getDensityCutoffedXYZ(x,y,z);
                    dN[j] += rho * (dx * dy * dz) / M_gas_ptcl;
                }
            }
        }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                       &dN[0], &recvcounts[0], &recvdispls[0], PS::GetDataType<PS::F64>(),
                       comm_info_y.getCommunicator());
#endif
        // Calculate the cumulative mass from y = - ycut.
        N_cum.resize(y_dm.size());
        N_cum[0] = dN[0];
        for (PS::S32 i = 1; i < N_cum.size(); i++) {
            N_cum[i] = N_cum[i-1] + dN[i-1];
        }
        if (my_rank == 0) {
            const PS::F64 fac = n_domain[0];
            std::cout << "N_cum      = " << N_cum.back() << std::endl;
            std::cout << "N_dm_cut   = " << M_dm_cut / (M_dm_ptcl * fac) << std::endl;
            std::cout << "N_star_cut = " << M_star_cut / (M_star_ptcl * fac) << std::endl;
            std::cout << "N_gas_cut  = " << M_gas_cut / (M_gas_ptcl * fac) << std::endl;
        }
        // Find dividing points in y
        for (PS::S32 i = 0; i < n_domain[1] - 1; i++) {
            const PS::F64 N_trgt = (N_dm + N_star + N_gas) * (i+1) 
                                 / (PS::F64) (n_domain[0] * n_domain[1]);
            PS::S32 idx = -1;
            for (PS::S32 k = 0; k < y_dm.size(); k++) {
                if (N_cum[k] > N_trgt) {
                    idx = k;
                    break;
                }
            }
            assert(idx >= 0);
            PS::F64 y;
            if (idx > 0) {
                const PS::F64 slope = (N_cum[idx] - N_cum[idx-1])/dy_dm[idx];
                y = (y_dm[idx] - 0.5*dy_dm[idx]) + (N_trgt - N_cum[idx-1]) / slope;
            } else {
                const PS::F64 slope = N_cum[idx]/dx_dm[idx];
                y = (y_dm[idx] - 0.5*dy_dm[idx]) + N_trgt / slope;
            }
            ydiv[i+1] = y;
        }
#else
        const PS::S32 rem = N_smpl_glb % n_proc_y; // remainder
        const PS::S32 quo = N_smpl_glb / n_proc_y; // quotient
        const PS::S32 N_smpl_loc = (my_rank_y < rem) ? (quo + 1) : quo;
        for (PS::S32 i = 0; i < n_domain[1] - 1; i++) {
            // Calculate the x coordinate of the boundary i+1/2.
            const PS::F64 N_trgt = N_dm * (i+1) / (PS::F64) (n_domain[0] * n_domain[1]);
            const PS::F64 x_min = xdiv[id_domain[0]];
            const PS::F64 x_max = xdiv[id_domain[0] + 1];
            PS::F64 y_min = ydiv[0];
            PS::F64 y_lo = ydiv[i];
            PS::F64 y_hi = ydiv[n_domain[1]];
            PS::F64 y_bisec = 0.5 * (y_lo + y_hi);
            PS::S32 iter = 0;
            for (;;) {
                // Estimate the number of particles in [x_min, x_bisec]
                const PS::F64 dv = (2.0 * rt_dm) * (x_max - x_min) * (y_bisec - y_min) / N_smpl_glb;
                PS::F64 N_eval {0.0}; 
                for (PS::S32 k = 0; k < N_smpl_loc; k++) {
                    const PS::F64 x = x_min + (x_max - x_min) * dist(run_param::prng::mt);
                    const PS::F64 y = y_min + (y_bisec - y_min) * dist(run_param::prng::mt);
                    const PS::F64 z = rt_dm * (2.0 * dist(run_param::prng::mt) - 1.0);
                    N_eval += model_dm.getDensityCutoffedXYZ(x,y,z) / M_dm_ptcl;
                }
                N_eval = comm_info_y.getSum(N_eval) * dv;
                // Check if the termination condition is satisfied or not
                const PS::F64 relerr = std::fabs((N_eval - N_trgt)/N_trgt);
                if (relerr <= eps) break;
                //if (my_rank == 0) std::cout << "relerr = " << relerr << std::endl;
                // Calculate the next x_bisec and update h_L & h_U
                if ((N_eval < N_trgt) && (y_lo < y_bisec)) y_lo = y_bisec;
                if ((N_eval > N_trgt) && (y_bisec < y_hi)) y_hi = y_bisec;
                y_bisec = 0.5*(y_lo + y_hi);
                iter++;
                // Check an error
                if (iter >= itermax) {
                    if (my_rank == 0) {
                        std::cout << "Too many iterations in the decomposition along Y!" << std::endl;
                    }
                    PS::Abort(-1); std::exit(-1);
                }
            }
            ydiv[i+1] = y_bisec;
        }
#endif
    }
    PS::Comm::barrier();
    if (my_rank == 0) std::cout << "took " << PS::GetWtime() - time_start << "[s]." << std::endl;
    // (3) decompose domain along the Z-direction
    if (my_rank == 0) std::cout << "Decomposing in Z..." << std::endl;
    const PS::S32 n_proc_z = comm_info_z.getNumberOfProc();
    const PS::S32 my_rank_z = comm_info_z.getRank();
    zdiv[0]           = - rt_dm;
    zdiv[n_domain[2]] =   rt_dm;
    time_start = PS::GetWtime();
    if (n_domain[2] > 1) {
#if 1
        // Make grids
        // [FIXME] Current version of the implementation makes grids
        //         assuming that the system is not rotated and translated.
        constexpr PS::S32 nz = 512;
        constexpr PS::S32 nx = 256;
        constexpr PS::S32 ny = nx;
        const PS::F64 rs_min = std::min({rs_dm, Rs_star, zs_star, Rs_gas, zs_gas});
        const PS::F64 rb_max = std::max({rt_dm, Rt_star, zt_star, Rt_gas, zt_gas});
        // (1) Make z grids
        //     (They must be the same because we're now trying to divide in z)
        defineLogarithmicGrid(z_dm, dz_dm, -rb_max, rb_max, 0.01 * rs_min, 0.0, nz);
        z_star.resize(z_dm.size());
        dz_star.resize(z_dm.size());
        z_gas.resize(z_dm.size());
        dz_gas.resize(z_dm.size());
        for (PS::S32 i = 0; i < z_dm.size(); i++) {
            z_star[i] = z_dm[i];
            dz_star[i] = dz_dm[i];
            z_gas[i] = z_dm[i];
            dz_gas[i] = dz_dm[i];
        }
        // (2) Make x & y grids
        //     (They CAN BE DIFFERENT for the type of particles)
        const PS::F64 xb_lo = xdiv[id_domain[0]];
        const PS::F64 xb_hi = xdiv[id_domain[0] + 1];
        const PS::F64 yb_lo = ydiv[id_domain[1]];
        const PS::F64 yb_hi = ydiv[id_domain[1] + 1];
        // (2-1) x & y grids for DM
        defineLogarithmicGrid(x_dm, dx_dm, xb_lo, xb_hi, 0.01 * rs_min, 0.0, nx);
        defineLogarithmicGrid(y_dm, dy_dm, yb_lo, yb_hi, 0.01 * rs_min, 0.0, ny);
        // (2-2) x & y grids for star
        const PS::F64 rs_min_star = std::min({Rs_star, zs_star});
        const PS::F64 rb_max_star = std::max({Rt_star, zt_star});
        const PS::F64 xb_lo_star = (xb_lo <= -rb_max_star && -rb_max_star <= xb_hi) ?  -rb_max_star : xb_lo;
        const PS::F64 xb_hi_star = (xb_lo <=  rb_max_star &&  rb_max_star <= xb_hi) ?   rb_max_star : xb_hi;
        const PS::F64 yb_lo_star = (yb_lo <= -rb_max_star && -rb_max_star <= yb_hi) ?  -rb_max_star : yb_lo;
        const PS::F64 yb_hi_star = (yb_lo <=  rb_max_star &&  rb_max_star <= yb_hi) ?   rb_max_star : yb_hi;
        defineLogarithmicGrid(x_star, dx_star, xb_lo_star, xb_hi_star, 0.01 * rs_min_star, 0.0, nx);
        defineLogarithmicGrid(y_star, dy_star, yb_lo_star, yb_hi_star, 0.01 * rs_min_star, 0.0, ny);
        // (2-3) x & y grids for gas
        const PS::F64 rs_min_gas = std::min({Rs_gas, zs_gas});
        const PS::F64 rb_max_gas = std::max({Rt_gas, zt_gas});
        const PS::F64 xb_lo_gas = (xb_lo <= -rb_max_gas && -rb_max_gas <= xb_hi) ?  -rb_max_gas : xb_lo;
        const PS::F64 xb_hi_gas = (xb_lo <=  rb_max_gas &&  rb_max_gas <= xb_hi) ?   rb_max_gas : xb_hi;
        const PS::F64 yb_lo_gas = (yb_lo <= -rb_max_gas && -rb_max_gas <= yb_hi) ?  -rb_max_gas : yb_lo;
        const PS::F64 yb_hi_gas = (yb_lo <=  rb_max_gas &&  rb_max_gas <= yb_hi) ?   rb_max_gas : yb_hi;
        defineLogarithmicGrid(x_gas, dx_gas, xb_lo_gas, xb_hi_gas, 0.01 * rs_min_gas, 0.0, nx);
        defineLogarithmicGrid(y_gas, dy_gas, yb_lo_gas, yb_hi_gas, 0.01 * rs_min_gas, 0.0, ny);
        // Calculate dN(z) in parallel
        dN.resize(z_dm.size());
        const PS::S32 rem = z_dm.size() % n_proc_z; // remainder
        const PS::S32 quo = z_dm.size() / n_proc_z; // quotient
        const PS::S32 sendcount = (my_rank_z < rem) ? (quo + 1) : quo;
        std::vector<PS::S32> recvcounts(n_proc_z), recvdispls(n_proc_z);
        comm_info_z.allGather(&sendcount, 1, &recvcounts[0]);
        recvdispls[0]=0;
        for (PS::S32 i = 1; i < n_proc_z; i++)
            recvdispls[i] = recvdispls[i-1] + recvcounts[i-1];
        const PS::S32 ksta = recvdispls[my_rank_z];
        const PS::S32 kend = ksta + sendcount - 1;
        for (PS::S32 k = 0; k < z_dm.size(); k++) {
            dN[k] = 0.0;
            const PS::F64 z = z_dm[k];
            const PS::F64 dz = dz_dm[k];
            // contribution from DM
            for (PS::S32 j = 0; j < y_dm.size(); j++) {
                const PS::F64 y = y_dm[j];
                const PS::F64 dy = dy_dm[j];
                for (PS::S32 i = 0; i < x_dm.size(); i++) {
                    const PS::F64 x = x_dm[i];
                    const PS::F64 dx = dx_dm[i];
                    PS::F64 r = std::sqrt(x*x + y*y + z*z);
                    PS::F64 rho;
                    if (r > 0.0) {
                        rho = model_dm.getDensityCutoffedXYZ(x,y,z);
                    } else {
                        // Handling the divergence at r = 0
                        const PS::F64 r_min_tmp = std::min(std::fabs(dx), std::min(std::fabs(dy), std::fabs(dz)));
                        const PS::F64 r_max_tmp = 0.5 * std::sqrt(dx*dx+dy*dy+dz*dz);
                        const PS::F64 rho_min = model_dm.getEnclosedMass(r_min_tmp)
                                              / (4.0 * math_const::pi * CUBE(r_min_tmp) / 3.0);
                        const PS::F64 rho_max = model_dm.getEnclosedMass(r_max_tmp)
                                              / (4.0 * math_const::pi * CUBE(r_max_tmp) / 3.0);
                        rho = 0.5 * (rho_min + rho_max);
                    }
                    dN[k] += rho * (dx * dy * dz) / M_dm_ptcl;
                }
            }
            // contribution from star
            for (PS::S32 j = 0; j < y_star.size(); j++) {
                const PS::F64 y = y_star[j];
                const PS::F64 dy = dy_star[j];
                for (PS::S32 i = 0; i < x_star.size(); i++) {
                    const PS::F64 x = x_star[i];
                    const PS::F64 dx = dx_star[i];
                    const PS::F64 rho = model_star.getDensityCutoffedXYZ(x,y,z);
                    dN[k] += rho * (dx * dy * dz) / M_star_ptcl;
                }
            }
            // contribution from gas
            for (PS::S32 j = 0; j < y_gas.size(); j++) {
                const PS::F64 y = y_gas[j];
                const PS::F64 dy = dy_gas[j];
                for (PS::S32 i = 0; i < x_gas.size(); i++) {
                    const PS::F64 x = x_gas[i];
                    const PS::F64 dx = dx_gas[i];
                    const PS::F64 rho = model_gas.getDensityCutoffedXYZ(x,y,z);
                    dN[k] += rho * (dx * dy * dz) / M_gas_ptcl;
                }
            }
        }
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                       &dN[0], &recvcounts[0], &recvdispls[0], PS::GetDataType<PS::F64>(),
                       comm_info_z.getCommunicator());
#endif
        // Calculate the cumulative mass from z = - zcut.
        N_cum.resize(z_dm.size());
        N_cum[0] = dN[0];
        for (PS::S32 i = 1; i < N_cum.size(); i++) {
            N_cum[i] = N_cum[i-1] + dN[i-1];
        }
        if (my_rank == 0) {
            const PS::F64 fac = n_domain[0] * n_domain[1];
            std::cout << "N_cum      = " << N_cum.back() << std::endl;
            std::cout << "N_dm_cut   = " << M_dm_cut / (M_dm_ptcl * fac) << std::endl;
            std::cout << "N_star_cut = " << M_star_cut / (M_star_ptcl * fac) << std::endl;
            std::cout << "N_gas_cut  = " << M_gas_cut / (M_gas_ptcl * fac) << std::endl;
        }
        // Find dividing points in z
        for (PS::S32 i = 0; i < n_domain[2] - 1; i++) {
            const PS::F64 N_trgt = (N_dm + N_star + N_gas) * (i+1) 
                                 / (PS::F64) (n_domain[0] * n_domain[1] * n_domain[2]);
            PS::S32 idx = -1;
            for (PS::S32 k = 0; k < z_dm.size(); k++) {
                if (N_cum[k] > N_trgt) {
                    idx = k;
                    break;
                }
            }
            assert(idx >= 0);
            PS::F64 z;
            if (idx > 0) {
                const PS::F64 slope = (N_cum[idx] - N_cum[idx-1])/dz_dm[idx];
                z = (z_dm[idx] - 0.5*dz_dm[idx]) + (N_trgt - N_cum[idx-1]) / slope;
            } else {
                const PS::F64 slope = N_cum[idx]/dx_dm[idx];
                z = (z_dm[idx] - 0.5*dz_dm[idx]) + N_trgt / slope;
            }
            zdiv[i+1] = z;
        }
#else
        const PS::S32 rem = N_smpl_glb % n_proc_z; // remainder
        const PS::S32 quo = N_smpl_glb / n_proc_z; // quotient
        const PS::S32 N_smpl_loc = (my_rank_z < rem) ? (quo + 1) : quo;
        for (PS::S32 i = 0; i < n_domain[2] - 1; i++) {
            // Calculate the z coordinate of the boundary i+1/2.
            const PS::F64 N_trgt = N_dm * (i+1) / (PS::F64) (n_domain[0] * n_domain[1] * n_domain[2]);
            const PS::F64 x_min = xdiv[id_domain[0]];
            const PS::F64 x_max = xdiv[id_domain[0] + 1];
            const PS::F64 y_min = ydiv[id_domain[1]];
            const PS::F64 y_max = ydiv[id_domain[1] + 1];
            PS::F64 z_min = zdiv[0];
            PS::F64 z_lo = zdiv[i];
            PS::F64 z_hi = zdiv[n_domain[2]];
            PS::F64 z_bisec = 0.5 * (z_lo + z_hi);
            PS::S32 iter = 0;
            for (;;) {
                // Estimate the number of particles in [x_min, x_bisec]
                const PS::F64 dv = (x_max - x_min) * (y_max - y_min) * (z_bisec - z_min) / N_smpl_glb;
                PS::F64 N_eval {0.0}; 
                for (PS::S32 k = 0; k < N_smpl_loc; k++) {
                    const PS::F64 x = x_min + (x_max - x_min) * dist(run_param::prng::mt);
                    const PS::F64 y = y_min + (y_max - y_min) * dist(run_param::prng::mt);
                    const PS::F64 z = z_min + (z_bisec - z_min) * dist(run_param::prng::mt);
                    N_eval += model_dm.getDensityCutoffedXYZ(x,y,z) / M_dm_ptcl;
                }
                N_eval = comm_info_z.getSum(N_eval) * dv;
                // Check if the termination condition is satisfied or not
                const PS::F64 relerr = std::fabs((N_eval - N_trgt)/N_trgt);
                if (relerr <= eps) break;
                //if (my_rank == 0) std::cout << "relerr = " << relerr << std::endl;
                // Calculate the next x_bisec and update h_L & h_U
                if ((N_eval < N_trgt) && (z_lo < z_bisec)) z_lo = z_bisec;
                if ((N_eval > N_trgt) && (z_bisec < z_hi)) z_hi = z_bisec;
                z_bisec = 0.5*(z_lo + z_hi);
                iter++;
                // Check an error
                if (iter >= itermax) {
                    if (my_rank == 0) {
                        std::cout << "Too many iterations in the decomposition along X!" << std::endl;
                    }
                    PS::Abort(-1); std::exit(-1);
                }
            }
            zdiv[i+1] = z_bisec;
        }
#endif
    }
    PS::Comm::barrier();
    if (my_rank == 0) std::cout << "took " << PS::GetWtime() - time_start << "[s]." << std::endl;
    // Set local domain
    PS::F64ort pos_my_domain;
    pos_my_domain.low_.x  = xdiv[id_domain[0]];
    pos_my_domain.high_.x = xdiv[id_domain[0] + 1];
    pos_my_domain.low_.y  = ydiv[id_domain[1]];
    pos_my_domain.high_.y = ydiv[id_domain[1] + 1];
    pos_my_domain.low_.z  = zdiv[id_domain[2]];
    pos_my_domain.high_.z = zdiv[id_domain[2] + 1];
    const PS::F64 vol_my_domain = pos_my_domain.getVolume();
    const PS::F64vec len_my_domain = pos_my_domain.getFullLength();
    const PS::F64ort pos_my_domain_sph = getDomainInSphericalCoordinate(pos_my_domain);
    // Place DM particles
    if (my_rank == 0) {
        std::cout << "Generating the distribution of DM particles..." << std::endl;
    }
    PS::Comm::barrier(); time_start = PS::GetWtime();
    // (1) Calculate the number of local DM particles
    const PS::S32 N_dm_loc = (PS::S32) math_util::qmc3d<PS::F64>(
        [model_dm, M_dm_ptcl](const PS::F64 x, const PS::F64 y, const PS::F64 z)
            -> PS::F64 {return model_dm.getDensityCutoffedXYZ(x,y,z) / M_dm_ptcl;},
        pos_my_domain.low_.x, pos_my_domain.high_.x,
        pos_my_domain.low_.y, pos_my_domain.high_.y,
        pos_my_domain.low_.z, pos_my_domain.high_.z,
        N_smpl_glb,
        run_param::prng::mt);
    std::vector<PS::S32> list(n_proc);
    PS::Comm::allGather(&N_dm_loc, 1, &list[0]);
    PS::S64 id_offset {0};
    for (PS::S32 i = 0; i < my_rank; i++) id_offset += list[i];
    // (2) Set the number of local DM particles
    psys_dm.setNumberOfParticleLocal(N_dm_loc);
    const PS::S32 N_dm_glb = psys_dm.getNumberOfParticleGlobal();
    // (3) Generate a distribution of DM particles
    PS::F64ort pos_domain_gen;
    bool gen_dm_dist = pos_domain_dm.calcIntersection(pos_my_domain,
                                                      pos_domain_gen);
    if (N_dm_loc > 0 && gen_dm_dist == false) {
        std::cout << "No intersection although N_dm_loc > 0." << std::endl;
        std::cout << "N_dm_loc = " << N_dm_loc << " rank = " << my_rank << std::endl;
        PS::Abort(-1); std::exit(-1);
    }
    PS::F64 msum {0.0};
    PS::F64vec pos_cm(0.0), vel_cm(0.0);
    if (gen_dm_dist) {
        const PS::F64 r_min = pos_my_domain_sph.low_.x;
        const PS::F64 r_max = pos_my_domain_sph.high_.x;
        const PS::F64 M_min = std::min(M_dm_cut, model_dm.getEnclosedMass(r_min));
        const PS::F64 M_max = std::min(M_dm_cut, model_dm.getEnclosedMass(r_max));
        const PS::F64 mu_min = pos_my_domain_sph.low_.y;
        const PS::F64 mu_max = pos_my_domain_sph.high_.y;
        const PS::F64 phi_min = pos_my_domain_sph.low_.z;
        const PS::F64 phi_max = pos_my_domain_sph.high_.z;
        for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
            psys_dm[i].id = i + id_offset;
            // Mass
            psys_dm[i].mass = M_dm_ptcl;
            msum += M_dm_ptcl;
            // Position
            PS::F64 r, theta, phi;
            PS::F64vec pos;
            for (;;) {
                const PS::F64 M = M_min + (M_max - M_min) * dist(run_param::prng::mt);
                r = model_dm.getRadiusFromMass(M);
                const PS::F64 mu = mu_min + (mu_max - mu_min) * dist(run_param::prng::mt);
                theta = std::acos(mu);
                phi = phi_min + (phi_max - phi_min) * dist(run_param::prng::mt);
                pos.x = r * std::sin( theta ) * std::cos( phi );
                pos.y = r * std::sin( theta ) * std::sin( phi );
                pos.z = r * std::cos( theta );
                if (pos_my_domain.contains(pos)) break;
            }
            psys_dm[i].pos.x = pos.x;
            psys_dm[i].pos.y = pos.y;
            psys_dm[i].pos.z = pos.z;
            pos_cm += psys_dm[i].mass * psys_dm[i].pos;
            // Velocity
            const PS::F64 pot = model_dm.getPotential(r);
            const PS::F64 v = model_dm.getVelocity(pot);
            theta = std::acos(2.0 * dist(run_param::prng::mt) - 1.0);
            phi = 2.0 * math_const::pi * dist(run_param::prng::mt);
            psys_dm[i].vel.x = v * std::sin( theta ) * std::cos( phi );
            psys_dm[i].vel.y = v * std::sin( theta ) * std::sin( phi );
            psys_dm[i].vel.z = v * std::cos( theta );
            vel_cm += psys_dm[i].mass * psys_dm[i].vel;
            // Clear
            psys_dm[i].acc  = 0.0;
            psys_dm[i].pot  = 0.0;
        }
    }
    msum = PS::Comm::getSum(msum);
    pos_cm = PS::Comm::getSum(pos_cm) / msum;
    vel_cm = PS::Comm::getSum(vel_cm) / msum;
    for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
        psys_dm[i].pos -= pos_cm;
        psys_dm[i].vel -= vel_cm;
    }
    PS::Comm::barrier();
    if (my_rank == 0) std::cout << "took " << PS::GetWtime() - time_start << "[s]." << std::endl;
    // Place star particles 
    if (my_rank == 0) {
        std::cout << "Generating the distribution of star particles..." << std::endl;
    }
    PS::Comm::barrier(); time_start = PS::GetWtime();
    // (1) Calculate the number of local star particles
    const PS::S32 N_star_loc = (PS::S32) math_util::qmc3d<PS::F64>(
        [model_star, M_star_ptcl](const PS::F64 x, const PS::F64 y, const PS::F64 z)
            -> PS::F64 {return model_star.getDensityCutoffedXYZ(x,y,z) / M_star_ptcl;},
        pos_my_domain.low_.x, pos_my_domain.high_.x,
        pos_my_domain.low_.y, pos_my_domain.high_.y,
        pos_my_domain.low_.z, pos_my_domain.high_.z,
        N_smpl_glb,
        run_param::prng::mt);
    PS::Comm::allGather(&N_star_loc, 1, &list[0]);
    id_offset = 0;
    for (PS::S32 i = 0; i < my_rank; i++) id_offset += list[i];
    // (2) Set the number of local star particles
    psys_star.setNumberOfParticleLocal(N_star_loc);
    const PS::S64 N_star_glb = psys_star.getNumberOfParticleGlobal();
    // (3) Generate a distribution of star particles
    bool gen_star_dist = pos_domain_star.calcIntersection(pos_my_domain,
                                                          pos_domain_gen);
    if (N_star_loc > 0 && gen_star_dist == false) {
        std::cout << "No intersection although N_star_loc > 0." << std::endl;
        std::cout << "N_star_loc = " << N_star_loc << " rank = " << my_rank << std::endl;
        PS::Abort(-1); std::exit(-1);
    }
    msum = 0.0;
    pos_cm = 0.0;
    if (gen_star_dist) {
        for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
            psys_star[i].id = i + N_dm_glb + id_offset;
            // Mass
            psys_star[i].mass0 = M_star_ptcl;
            psys_star[i].mass  = M_star_ptcl;
            msum += M_star_ptcl;
            // Position
            const PS::F64vec len = pos_domain_gen.getFullLength();
            const PS::F64 fmax = model_star.getDensityCutoffedXYZ(0.0, 0.0, 0.0);
            PS::F64vec pos;
            for (;;) {
                pos.x = pos_domain_gen.low_.x + len.x * dist(run_param::prng::mt);
                pos.y = pos_domain_gen.low_.y + len.y * dist(run_param::prng::mt);
                pos.z = pos_domain_gen.low_.z + len.z * dist(run_param::prng::mt);
                const PS::F64 f = fmax * dist(run_param::prng::mt);
                const PS::F64 f_xyz = model_star.getDensityCutoffedXYZ(pos.x, pos.y, pos.z);
                if (f < f_xyz) break;
            }
            psys_star[i].pos.x = pos.x;
            psys_star[i].pos.y = pos.y;
            psys_star[i].pos.z = pos.z;
            pos_cm += psys_star[i].mass * psys_star[i].pos;
            // Velocity
            psys_star[i].vel = 0.0; // determined after the 1st gravity calculation
            // Acceleration, etc.
            psys_star[i].acc = 0.0;
            psys_star[i].pot = 0.0;
            // Quantities related to feed back
            psys_star[i].mabn[0]  = run_param::ism::Xhydrogen_solar;
            psys_star[i].mabn[1]  = run_param::ism::Yhelium_solar;
            psys_star[i].mabn[2]  = run_param::ism::Zcarbon_solar;
            psys_star[i].mabn[3]  = run_param::ism::Znitrogen_solar;
            psys_star[i].mabn[4]  = run_param::ism::Zoxygen_solar;
            psys_star[i].mabn[5]  = run_param::ism::Zneon_solar;
            psys_star[i].mabn[6]  = run_param::ism::Zmagnesium_solar;
            psys_star[i].mabn[7]  = run_param::ism::Zsilicon_solar;
            psys_star[i].mabn[8]  = run_param::ism::Zsulfur_solar;
            psys_star[i].mabn[9]  = run_param::ism::Zcalcium_solar;
            psys_star[i].mabn[10] = run_param::ism::Ziron_solar;
            psys_star[i].mabn[11] = run_param::ism::Znickel_solar;
            psys_star[i].mabn[12] = run_param::ism::Zeuropium_solar;
            psys_star[i].t_form = - 10 * phys_const::Gyr; 
            psys_star[i].FBcnt = 0;  
            psys_star[i].t_SNII = 0.0; 
            psys_star[i].t_SNIa = 0.0; 
            psys_star[i].t_AGB = 0.0;  
            psys_star[i].t_NSM = 0.0;  
            psys_star[i].FBrad = 0.0;  
        }
    }
    msum = PS::Comm::getSum(msum);
    pos_cm = PS::Comm::getSum(pos_cm) / msum;
    for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
        psys_star[i].pos -= pos_cm;
    }
    PS::Comm::barrier();
    if (my_rank == 0) std::cout << "took " << PS::GetWtime() - time_start << "[s]." << std::endl;
    // Place SPH particles 
    if (my_rank == 0) {
        std::cout << "Generating the distribution of gas particles..." << std::endl;
    }
    PS::Comm::barrier(); time_start = PS::GetWtime();
    // (1) Calculate the number of local gas particles
    const PS::S32 N_gas_loc = (PS::S32) math_util::qmc3d<PS::F64>(
        [model_gas, M_gas_ptcl](const PS::F64 x, const PS::F64 y, const PS::F64 z)
            -> PS::F64 {return model_gas.getDensityCutoffedXYZ(x,y,z) / M_gas_ptcl;},
        pos_my_domain.low_.x, pos_my_domain.high_.x,
        pos_my_domain.low_.y, pos_my_domain.high_.y,
        pos_my_domain.low_.z, pos_my_domain.high_.z,
        N_smpl_glb,
        run_param::prng::mt);
    PS::Comm::allGather(&N_gas_loc, 1, &list[0]);
    id_offset = 0;
    for (PS::S32 i = 0; i < my_rank; i++) id_offset += list[i];
    // (2) Set the number of local gas particles
    psys_gas.setNumberOfParticleLocal(N_gas_loc);
    const PS::S64 N_gas_glb = psys_gas.getNumberOfParticleGlobal();
    // (3) Generate a distribution of gas particles
    bool gen_gas_dist = pos_domain_gas.calcIntersection(pos_my_domain,
                                                        pos_domain_gen);
    if (N_gas_loc > 0 && gen_gas_dist == false) {
        std::cout << "No intersection although N_gas_loc > 0." << std::endl;
        std::cout << "N_gas_loc = " << N_gas_loc << " rank = " << my_rank << std::endl;
        PS::Abort(-1); std::exit(-1);
    }
    msum = 0.0;
    pos_cm = 0.0;
    if (gen_gas_dist) {
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            psys_gas[i].id = i + N_dm_glb + N_star_glb + id_offset;
            // Mass
            psys_gas[i].mass0 = M_gas_ptcl;
            psys_gas[i].mass  = M_gas_ptcl;
            msum += M_gas_ptcl; 
            // Position
            const PS::F64vec len = pos_domain_gen.getFullLength();
            const PS::F64 fmax = model_gas.getDensityCutoffedXYZ(0.0, 0.0, 0.0);
            PS::F64vec pos;
            for (;;) {
                pos.x = pos_domain_gen.low_.x + len.x * dist(run_param::prng::mt);
                pos.y = pos_domain_gen.low_.y + len.y * dist(run_param::prng::mt);
                pos.z = pos_domain_gen.low_.z + len.z * dist(run_param::prng::mt);
                const PS::F64 f = fmax * dist(run_param::prng::mt);
                const PS::F64 f_xyz = model_gas.getDensityCutoffedXYZ(pos.x, pos.y, pos.z);
                if (f < f_xyz) break;
            }
            psys_gas[i].pos.x = pos.x;
            psys_gas[i].pos.y = pos.y;
            psys_gas[i].pos.z = pos.z;
            pos_cm += psys_gas[i].mass * psys_gas[i].pos;
            // Velocity
            psys_gas[i].vel = 0.0; // determined later
            // Other quantities
            psys_gas[i].acc_grav  = 0.0;
            psys_gas[i].pot_grav  = 0.0;
            psys_gas[i].acc_hydro = 0.0;
            psys_gas[i].eng = (phys_const::kBoltz * temp)/((run_param::ism::gamma - 1.0) * mu * phys_const::Mproton);
            psys_gas[i].h = std::pow(Rt_gas*Rt_gas*zt_gas, 1.0/3.0) * std::pow((PS::F64) run_param::sph::N_ngb / (PS::F64)N_gas, 1.0/3.0);
            psys_gas[i].alpha = run_param::sph::alpha_AV_ini;
            psys_gas[i].mabn[0]  = run_param::ism::Xhydrogen_solar;
            psys_gas[i].mabn[1]  = run_param::ism::Yhelium_solar;
            psys_gas[i].mabn[2]  = run_param::ism::Zcarbon_solar;
            psys_gas[i].mabn[3]  = run_param::ism::Znitrogen_solar;
            psys_gas[i].mabn[4]  = run_param::ism::Zoxygen_solar;
            psys_gas[i].mabn[5]  = run_param::ism::Zneon_solar;
            psys_gas[i].mabn[6]  = run_param::ism::Zmagnesium_solar;
            psys_gas[i].mabn[7]  = run_param::ism::Zsilicon_solar;
            psys_gas[i].mabn[8]  = run_param::ism::Zsulfur_solar;
            psys_gas[i].mabn[9]  = run_param::ism::Zcalcium_solar;
            psys_gas[i].mabn[10] = run_param::ism::Ziron_solar;
            psys_gas[i].mabn[11] = run_param::ism::Znickel_solar;
            psys_gas[i].mabn[12] = run_param::ism::Zeuropium_solar;
        }
    }
    msum = PS::Comm::getSum(msum);
    pos_cm = PS::Comm::getSum(pos_cm) / msum;
    for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
        psys_gas[i].pos -= pos_cm;
    }
    PS::Comm::barrier();
    if (my_rank == 0) std::cout << "took " << PS::GetWtime() - time_start << "[s]." << std::endl;
    // Unit convertion 
    run_param::unit::mass = 0.0;
    run_param::unit::leng = 0.0;
    for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
        run_param::unit::mass += psys_dm[i].mass;
        const PS::F64 r = std::sqrt(SQ(psys_dm[i].pos));
        if (r > run_param::unit::leng) run_param::unit::leng = r;
    }
    for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
        run_param::unit::mass += psys_star[i].mass;
        const PS::F64 r = std::sqrt(SQ(psys_star[i].pos));
        if (r > run_param::unit::leng) run_param::unit::leng = r;
    }
    for (PS::S32 i = 0; i< psys_gas.getNumberOfParticleLocal(); i++) {
        run_param::unit::mass += psys_gas[i].mass;
        const PS::F64 r = std::sqrt(SQ(psys_gas[i].pos));
        if (r > run_param::unit::leng) run_param::unit::leng = r;
    }
    run_param::unit::mass = PS::Comm::getSum(run_param::unit::mass);
    run_param::unit::leng = PS::Comm::getMaxValue(run_param::unit::leng);
    if (my_rank == 0) {
        std::cout << "Total mass of the system = "
                  << run_param::unit::mass/phys_const::Msolar << " [Msolar]" << std::endl;
        std::cout << "The size of the system = "
                  << run_param::unit::leng/phys_const::kpc << " [kpc]" << std::endl;
    }
    run_param::unit::time = std::sqrt(CUBE(run_param::unit::leng)/(phys_const::Ggrav * run_param::unit::mass));
    run_param::unit::velc = run_param::unit::leng / run_param::unit::time;
    run_param::unit::eng  = run_param::unit::mass * SQ(run_param::unit::velc);
    run_param::unit::spen = run_param::unit::eng  / run_param::unit::mass;
    for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
        psys_dm[i].mass /= run_param::unit::mass;
        psys_dm[i].pos  /= run_param::unit::leng;
        psys_dm[i].vel  /= run_param::unit::velc;
    }
    for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
        psys_star[i].mass0  /= run_param::unit::mass;
        psys_star[i].mass   /= run_param::unit::mass;
        psys_star[i].pos    /= run_param::unit::leng;
        psys_star[i].vel    /= run_param::unit::velc;
        psys_star[i].t_form /= run_param::unit::time;
        psys_star[i].t_SNII /= run_param::unit::time;
        psys_star[i].t_SNIa /= run_param::unit::time;
        psys_star[i].t_AGB  /= run_param::unit::time;
        psys_star[i].t_NSM  /= run_param::unit::time;
    }
    for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
        psys_gas[i].mass0 /= run_param::unit::mass;
        psys_gas[i].mass  /= run_param::unit::mass;
        psys_gas[i].pos   /= run_param::unit::leng;
        psys_gas[i].vel   /= run_param::unit::velc;
        psys_gas[i].h     /= run_param::unit::leng;
        psys_gas[i].eng   /= run_param::unit::spen;
    }
    // Set boundary condition
    run_param::basic::bc = PS::BOUNDARY_CONDITION_OPEN;
    // Set the parameters for gravity calculation
    run_param::grav::soft::theta = 0.5;
    run_param::grav::soft::n_group_limit = 64;
    run_param::grav::soft::n_leaf_limit = 8;
    run_param::grav::soft::eps_dm = 1.0e2 * phys_const::pc / run_param::unit::leng;
    run_param::grav::soft::eps_star = 3.0 * phys_const::pc / run_param::unit::leng;
    run_param::grav::soft::eps_gas = 3.0 * phys_const::pc / run_param::unit::leng;
    run_param::grav::soft::dt_fid = 1.0e5 * phys_const::yr / run_param::unit::time;
    run_param::grav::hard::eps = 0.01 * phys_const::pc / run_param::unit::leng;
    run_param::grav::hard::eta = 0.1;
    run_param::grav::hard::eta_s = 1.0e-3;
    run_param::grav::hard::dt_limit = run_param::grav::soft::dt_fid / 8;
    // Set the parameters for SPH calculation
    run_param::sph::n_group_limit = 64;
    run_param::sph::n_leaf_limit = 1;
    run_param::sph::n_jp_limit = 4096;
    // Set I/O intervals
    run_param::io::dt_dump = 0.001;
    run_param::io::time_dump = run_param::io::dt_dump;
    run_param::io::dt_dump_rst = 0.001;
    run_param::io::time_dump_rst = run_param::io::dt_dump;
    // Set the end time of the simulation
    run_param::basic::time_end = 10.0;

    if (PS::Comm::getRank() == 0) {
        std::cout << "An initial condition for isolated galaxy simulation is made." << std::endl;
        std::cout << "N_dm_glb = " << N_dm_glb 
                  << ", N_star_glb = " << N_star_glb
                  << ", N_gas_glb = " << N_gas_glb << std::endl;
    }
    dbg_utils::fout << (N_dm_loc + N_star_loc + N_gas_loc) << std::endl;
    PS::Finalize(); std::exit(0);

#if 0
    // Compare the mass density profile due to the dark matter with the analytic one
    {
        // Set the resolution of cells
        const PS::S32 Nr = 8192;
        const PS::F64 safety = 1.001;
        // Compute the maximum spherical radius 
        PS::F64 rmax = 0.0;
        for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
            const PS::F64 r = std::sqrt(SQ(psys_dm[i].pos));
            rmax = std::max(rmax, r);
        }
        rmax = PS::Comm::getMaxValue(rmax);
        // Compute the mass density
        const PS::F64 dr = (safety * rmax)/Nr;
        PS::F64 rho_loc[Nr] = {};
        for (PS::S32 i = 0; i < psys_dm.getNumberOfParticleLocal(); i++) {
            const PS::F64 r = std::sqrt(SQ(psys_dm[i].pos));
            const PS::S32 ii = r/dr;
            rho_loc[ii] += psys_dm[i].mass;
        }
        PS::F64 rho_tot[Nr];
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Reduce(rho_loc, rho_tot, Nr, PS::GetDataType<PS::F64>(),
                   MPI_SUM, 0, comm_info.getCommunicator());
#else
        for (PS::S32 i = 0; i < Nr; i++) rho_tot[i]=rho_loc[i];
#endif
        // Output the mass density
        if (my_rank == 0) {
            for (PS::S32 i = 0; i < Nr; i++) {
                const PS::F64 r = (i + 0.5) * dr;
                rho_tot[i] /= (4.0 * math_const::pi * r * r) * dr;
            }
            const std::string file_name = "rho_dm.txt";
            std::ofstream ofs;
            ofs.open(file_name.c_str(), std::ios::trunc);
            ofs.setf(std::ios_base::scientific, std::ios_base::floatfield);
            ofs << std::setprecision(15) << std::showpos;
            for (PS::S32 i = 0; i < Nr; i++) {
                const PS::F64 r = (i + 0.5) * dr;
                const PS::F64 r_cgs = r * run_param::unit::leng;
                const PS::F64 rho_exact = model_dm.getDensity(r_cgs)
                                        / (run_param::unit::mass/CUBE(run_param::unit::leng));
                ofs << r << "   "
                    << rho_tot[i] << "   "
                    << rho_exact 
                    << std::endl;
            }
            ofs.close();
        }
    }
    // Compare the mass density profile due to the stars with the analytic one
    {
        // Set the resolution of cells
        const PS::S32 NR = 128, Nz = 128;
        const PS::F64 safety = 1.001;
        // Compute the maximul cylindrical radius and |z| of star particles
        PS::F64 Rmax = 0.0, zmax = 0.0;
        for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
            const PS::F64 R = std::sqrt(SQ(psys_star[i].pos.x)
                                       +SQ(psys_star[i].pos.y));
            const PS::F64 z = std::fabs(psys_star[i].pos.z);
            Rmax = std::max(Rmax, R);
            zmax = std::max(zmax, z);
        }
        Rmax = PS::Comm::getMaxValue(Rmax);
        zmax = PS::Comm::getMaxValue(zmax);
        // Compute the mass density
        const PS::F64 dR = (safety * Rmax)/NR;
        const PS::F64 dz = (safety * zmax)/Nz;
        PS::F64 rho_loc[NR * Nz] = {};
        for (PS::S32 i = 0; i < psys_star.getNumberOfParticleLocal(); i++) {
            const PS::F64 R = std::sqrt(SQ(psys_star[i].pos.x)
                                       +SQ(psys_star[i].pos.y));
            const PS::F64 z = std::fabs(psys_star[i].pos.z);
            const PS::S32 ii = R/dR;
            const PS::S32 jj = z/dz;
            rho_loc[ii + NR * jj] += psys_star[i].mass;
        }
        PS::F64 rho_tot[NR * Nz];
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Reduce(rho_loc, rho_tot, NR * Nz, PS::GetDataType<PS::F64>(),
                   MPI_SUM, 0, comm_info.getCommunicator());
#else
        for (PS::S32 i = 0; i < NR * Nz; i++) rho_tot[i]=rho_loc[i];
#endif
        // Output the mass density
        if (my_rank == 0) {
            for (PS::S32 j = 0; j < Nz; j++) {
                for (PS::S32 i = 0; i < NR; i++) {
                    const PS::F64 R = (i + 0.5) * dR;
                    rho_tot[i + NR * j] /= 2.0 * (2.0 * math_const::pi * R) * (dR * dz);
                    // The 1st factor of 2 of the RHS is necessary because we use
                    // particles at z<0 to calculate rho(R,|z|).
                }
            }
            const std::string file_name = "rho_star.txt";
            std::ofstream ofs;
            ofs.open(file_name.c_str(), std::ios::trunc);
            ofs.setf(std::ios_base::scientific, std::ios_base::floatfield);
            ofs << std::setprecision(15) << std::showpos;
            for (PS::S32 j = 0; j < Nz; j++) {
                const PS::F64 z = (j + 0.5) * dz;
                const PS::F64 z_cgs = z * run_param::unit::leng;
                for (PS::S32 i = 0; i < NR; i++) {
                    const PS::F64 R = (i + 0.5) * dR;
                    const PS::F64 R_cgs = R * run_param::unit::leng;
                    const PS::F64 rho_exact = model_star.getDensity(R_cgs,z_cgs)
                                            / (run_param::unit::mass/CUBE(run_param::unit::leng));
                    ofs << R << "   " << z << "   "
                        << rho_tot[i + NR * j] << "   "
                        << rho_exact 
                        << std::endl;
                }
            }
            ofs.close();
        }
    }
    // Compute the mass gas density and output it.
    {
        // Set the resolution of cells, etc.
        const PS::S32 NR=128, Nz = 128;
        const PS::F64 safety = 1.001;
        // Compute the maxium clyndrical radius of gas particles
        PS::F64 Rmax = 0.0, zmax = 0.0;
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            const PS::F64 R = std::sqrt(psys_gas[i].pos.x * psys_gas[i].pos.x
                                       +psys_gas[i].pos.y * psys_gas[i].pos.y);
            const PS::F64 z = std::fabs(psys_gas[i].pos.z);
            Rmax = std::max(Rmax, R);
            zmax = std::max(zmax, z);
        }
        Rmax = PS::Comm::getMaxValue(Rmax);
        zmax = PS::Comm::getMaxValue(zmax);
        // Compute the surface density
        const PS::F64 dR = (safety * Rmax)/NR;
        const PS::F64 dz = (safety * zmax)/Nz;
        PS::F64 rho_loc[NR * Nz] = {};
        for (PS::S32 i = 0; i < psys_gas.getNumberOfParticleLocal(); i++) {
            const PS::F64 R = std::sqrt(SQ(psys_gas[i].pos.x)
                                       +SQ(psys_gas[i].pos.y));
            const PS::F64 z = std::fabs(psys_gas[i].pos.z);
            const PS::S32 ii = R/dR;
            const PS::S32 jj = z/dz;
            rho_loc[ii + NR * jj] += psys_gas[i].mass;
        }
        PS::F64 rho_tot[NR * Nz];
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Reduce(rho_loc, rho_tot, NR * Nz, PS::GetDataType<PS::F64>(),
                   MPI_SUM, 0, comm_info.getCommunicator());
#else
        for (PS::S32 i = 0; i < NR * Nz; i++) rho_tot[i] = rho_loc[i];
#endif
        // Output the mass density
        if (my_rank == 0) {
            for (PS::S32 j = 0; j < Nz; j++) {
                for (PS::S32 i = 0; i < NR; i++) {
                    const PS::F64 R = (i + 0.5) * dR;
                    rho_tot[i + NR * j] /= 2.0 * (2.0 * math_const::pi * R) * (dR * dz);
                }
            }
            const std::string file_name = "rho_gas.txt";
            std::ofstream ofs;
            ofs.open(file_name.c_str(),std::ios::trunc);
            ofs.setf(std::ios_base::scientific, std::ios_base::floatfield);
            ofs << std::setprecision(15) << std::showpos;
            for (PS::S32 j = 0; j < Nz; j++) {
                const PS::F64 z = (j + 0.5) * dz;
                const PS::F64 z_cgs = z * run_param::unit::leng;
                for (PS::S32 i = 0; i < NR; i++) {
                    const PS::F64 R = (i + 0.5) * dR;
                    const PS::F64 R_cgs = R * run_param::unit::leng;
                    const PS::F64 rho_exact = model_gas.getDensity(R_cgs,z_cgs)
                                            / (run_param::unit::mass/CUBE(run_param::unit::leng));
                    ofs << R << "   " << z << "   "
                        << rho_tot[i + NR * j] << "   "
                        << rho_exact 
                        << std::endl;
                }
            }
            ofs.close();
        }
    }
#endif
}


#endif // ASURA_FDPS_NORMAL_RUN

#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_GLASS_DATA_GENERATION_MODE
void GlassDataGenerationModeIC(PS::ParticleSystem<FP_dm>& psys_dm,
                               PS::ParticleSystem<FP_star>& psys_star,
                               PS::ParticleSystem<FP_gas>& psys_gas) {
    // This test uses SPH particles only
    constexpr PS::S32 N_ptcl = (1<<18);
    psys_dm.setNumberOfParticleLocal(0);
    psys_star.setNumberOfParticleLocal(0); 
    psys_gas.setNumberOfParticleLocal(N_ptcl);
    // Set the size of the computational domain
    run_param::basic::pos_root_domain = PS::F64ort(PS::F64vec(0.0),
                                                   PS::F64vec(1.0));
    // Set the physical quantities of particles
    std::uniform_real_distribution<double> dist(0.0, 1.0);
#if !defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
    const PS::F64 h_ini = std::cbrt((PS::F64)run_param::sph::N_ngb/(PS::F64)N_ptcl);
#else
    const PS::F64 h_ini = std::sqrt((PS::F64)run_param::sph::N_ngb/(PS::F64)N_ptcl);
#endif
    for (PS::S32 i = 0; i < N_ptcl; i++) {
        psys_gas[i].id  = i;
        psys_gas[i].mass0 = 1.0 / N_ptcl;
        psys_gas[i].mass = 1.0 / N_ptcl;
        psys_gas[i].pos.x = dist(run_param::prng::mt);
        psys_gas[i].pos.y = dist(run_param::prng::mt);
#if !defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
        psys_gas[i].pos.z = dist(run_param::prng::mt);
#endif
        psys_gas[i].vel = 0.0;
        psys_gas[i].eng = 1.0;
        psys_gas[i].h = h_ini;
        psys_gas[i].alpha = run_param::sph::alpha_AV_ini;
        for (PS::S32 k = 0; k < CELibYield_Number; k++)
            psys_gas[i].mabn[k] = 0.0;
    }
    // Set the simulation unit 
    run_param::unit::mass = 1.0;
    run_param::unit::leng = 1.0;
    run_param::unit::time = 1.0;
    // Set boundary condition
    run_param::basic::bc = PS::BOUNDARY_CONDITION_PERIODIC_XYZ;
    // Set the parameters for gravity calculation
    // (Dummy values are set for the parameters that are not used actually)
    run_param::grav::soft::theta = 0.5;
    run_param::grav::soft::n_group_limit = 64;
    run_param::grav::soft::n_leaf_limit = 8;
    run_param::grav::soft::eps_dm = 1.0e-2 / run_param::unit::leng;
    run_param::grav::soft::eps_star = 1.0e-2 / run_param::unit::leng;
    run_param::grav::soft::eps_gas = 1.0e-2 / run_param::unit::leng;
    run_param::grav::soft::dt_fid = 0.01 / run_param::unit::time;
    run_param::grav::hard::eps = 1.0e-2 / run_param::unit::leng;
    run_param::grav::hard::eta = 0.1;
    run_param::grav::hard::eta_s = 1.0e-3;
    run_param::grav::hard::dt_limit = run_param::grav::soft::dt_fid / 8;
    // Set the parameters for SPH calculation
    run_param::sph::n_group_limit = 64;
    run_param::sph::n_leaf_limit = 8;
    run_param::sph::n_jp_limit = 4096;
    // Set I/O intervals
    const PS::F64 tcross = std::sqrt(3.0);
    run_param::io::dt_dump = tcross / run_param::unit::time;
    run_param::io::time_dump = run_param::io::dt_dump;
    run_param::io::dt_dump_rst = tcross / run_param::unit::time;
    run_param::io::time_dump_rst = run_param::io::dt_dump;
    // Set the end time of the simulation
    run_param::basic::time_end = (64.0 * tcross) / run_param::unit::time;

    if (PS::Comm::getRank() == 0) {
        std::cout << "An initial condition for glass data generation is made." << std::endl;
        std::cout << "N_sph = " << N_ptcl << std::endl;
        std::cout << "h_ini = " << h_ini << std::endl;
    }
}
#endif // ASURA_FDPS_GLASS_DATA_GENERATION_MODE

#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_PARTICLE_COLLISION_TEST
void ParticleCollisionTestIC(PS::ParticleSystem<FP_dm>& psys_dm,
                             PS::ParticleSystem<FP_star>& psys_star,
                             PS::ParticleSystem<FP_gas>& psys_gas) {
    // This test uses SPH particles only
    constexpr PS::S32 N_ptcl = 2;
    psys_dm.setNumberOfParticleLocal(0);
    psys_star.setNumberOfParticleLocal(0); 
    psys_gas.setNumberOfParticleLocal(N_ptcl);
    // Set the size of the computational domain
    run_param::basic::pos_root_domain = PS::F64ort(PS::F64vec(0.0),
                                                   PS::F64vec(1.0));
    // Set the physical quantities of particles
    constexpr PS::F64 h_ini = 0.3;
    constexpr PS::F64 v_ini = 0.1;
    // (1) Left particle
    {
        const PS::S32 idx = 0;
        psys_gas[idx].id  = idx;
        psys_gas[idx].mass0 = 1.0;
        psys_gas[idx].mass = 1.0;
        psys_gas[idx].pos = PS::F64vec(0.5-0.1, 0.0);
        psys_gas[idx].vel = PS::F64vec(v_ini, 0.0);
        psys_gas[idx].h = h_ini;
        psys_gas[idx].alpha = run_param::sph::alpha_AV_ini;
        psys_gas[idx].eng = 1.0;
        for (PS::S32 k = 0; k < CELibYield_Number; k++)
            psys_gas[idx].mabn[k] = 0.0;
    }
    // (2) Right particle
    {
        const PS::S32 idx = 1;
        psys_gas[idx].id  = idx;
        psys_gas[idx].mass0 = 1.0;
        psys_gas[idx].mass = 1.0;
        psys_gas[idx].pos = PS::F64vec(0.5+0.1, 0.0);
        psys_gas[idx].vel = PS::F64vec(-v_ini, 0.0);
        psys_gas[idx].h = h_ini;
        psys_gas[idx].alpha = run_param::sph::alpha_AV_ini;
        psys_gas[idx].eng = 1.0;
        for (PS::S32 k = 0; k < CELibYield_Number; k++)
            psys_gas[idx].mabn[k] = 0.0;
    }
    // Set the simulation unit 
    run_param::unit::mass = 1.0;
    run_param::unit::leng = 1.0;
    run_param::unit::time = 1.0;
    // Set boundary condition
    run_param::basic::bc = PS::BOUNDARY_CONDITION_PERIODIC_XYZ;
    // Set the parameters for gravity calculation
    // (Dummy values are set for the parameters that are not used actually)
    run_param::grav::soft::theta = 0.5;
    run_param::grav::soft::n_group_limit = 64;
    run_param::grav::soft::n_leaf_limit = 8;
    run_param::grav::soft::eps_dm = 1.0e-2 / run_param::unit::leng;
    run_param::grav::soft::eps_star = 1.0e-2 / run_param::unit::leng;
    run_param::grav::soft::eps_gas = 1.0e-2 / run_param::unit::leng;
    run_param::grav::soft::dt_fid = 0.01 / run_param::unit::time;
    run_param::grav::hard::eps = 1.0e-2 / run_param::unit::leng;
    run_param::grav::hard::eta = 0.1;
    run_param::grav::hard::eta_s = 1.0e-3;
    run_param::grav::hard::dt_limit = run_param::grav::soft::dt_fid / 8;
    // Set the parameters for SPH calculation
    run_param::sph::n_group_limit = 64;
    run_param::sph::n_leaf_limit = 8;
    run_param::sph::n_jp_limit = 4096;
    // Set I/O intervals
    run_param::io::dt_dump = 0.1 / run_param::unit::time;
    run_param::io::time_dump = run_param::io::dt_dump;
    run_param::io::dt_dump_rst = 0.1 / run_param::unit::time;
    run_param::io::time_dump_rst = run_param::io::dt_dump;
    // Set the end time of the simulation
    run_param::basic::time_end = 8.0 / run_param::unit::time;

    if (PS::Comm::getRank() == 0) {
        std::cout << "An initial condition for particle collision test is made." << std::endl;
        std::cout << "N_sph = " << N_ptcl << std::endl;
        std::cout << "h_ini = " << h_ini << std::endl;
    }
}
#endif // ASURA_FDPS_PARTICLE_COLLISION_TEST

#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_SHOCK_TUBE_TEST
void ShockTubeTestIC(PS::ParticleSystem<FP_dm>& psys_dm,
                     PS::ParticleSystem<FP_star>& psys_star,
                     PS::ParticleSystem<FP_gas>& psys_gas) {
    // Set the parameters of shock tube test
    // (Sod problem [Test1 in Toro's book; tstop=0.25])
    constexpr PS::F64 rho_L = 1.0;
    constexpr PS::F64 pres_L = 1.0;
    constexpr PS::F64 rho_R = 0.125;
    constexpr PS::F64 pres_R = 0.1;
    constexpr PS::F64 t_end = 0.25;
    // (common parameters)
    constexpr PS::F64 xlen_L = 1.0;
    constexpr PS::F64 xlen_R = 1.0;
    constexpr PS::F64 ylen   = (xlen_L + xlen_R)/16.0;
    constexpr PS::F64 zlen   = ylen;
    constexpr PS::F64 frac_L = (rho_L * xlen_L)/(rho_L * xlen_L + rho_R * xlen_R);
    constexpr PS::F64 frac_R = 1.0 - frac_L;
    // Set the simulation unit
    run_param::unit::mass = 1.0;
    run_param::unit::leng = 1.0;
    run_param::unit::time = 1.0;
    // Set the computational domain and the size of the dense region
    PS::F64ort pos_root_domain;
    pos_root_domain.low_.x  = - xlen_L;         
    pos_root_domain.high_.x = xlen_R;
    pos_root_domain.low_.y  = - 0.5 * ylen;
    pos_root_domain.high_.y =   0.5 * ylen;
    pos_root_domain.low_.z  = - 0.5 * zlen;
    pos_root_domain.high_.z =   0.5 * zlen;
    run_param::basic::pos_root_domain = pos_root_domain;
    // Read glass data
    PS::F64ort pos_glass_data;
    std::vector<PS::F64vec> ptcl;
    const std::string filename = "result/glass_data.txt";
    std::ifstream ifs;
    ifs.open(filename.c_str(), std::ios::in | std::ios::binary);
    if (ifs) {
        while(true) {
            PS::F64vec pos;
            ifs.read((char *)&pos, sizeof(PS::F64vec));
            if (ifs.eof()) break;
            ptcl.push_back(pos);
            pos_glass_data.merge(pos);
        }
    } else {
        std::cout << "Could not open file " << filename << std::endl;
        assert(false);
    }
    ifs.close();
    const PS::S64 N_ptcl = ptcl.size();
    std::cout << N_ptcl << " particles are read." << std::endl;
    // Allocate memory for particles
    psys_dm.setNumberOfParticleLocal(0);
    psys_star.setNumberOfParticleLocal(0);
    psys_gas.setNumberOfParticleLocal(N_ptcl);
    // Rescale the position and the mass of the particles
    const PS::F64 M_tot = (rho_L*xlen_L + rho_R*xlen_R)*ylen*zlen;
    const PS::F64 h_ini = 3.0 / std::cbrt(N_ptcl);
    const PS::F64vec lo = pos_glass_data.low_;
    const PS::F64vec hi = pos_glass_data.high_;
    const PS::F64vec cen = pos_glass_data.getCenter();
    const PS::F64vec len = pos_glass_data.getFullLength();
    const PS::F64vec lo_root_domain = pos_root_domain.low_;
    const PS::F64vec hi_root_domain = pos_root_domain.high_;
    const PS::F64vec cen_root_domain = pos_root_domain.getCenter();
    const PS::F64vec len_root_domain = pos_root_domain.getFullLength();
    for (PS::S32 i = 0; i < N_ptcl; i++) {
        // Rescale mass
        psys_gas[i].mass0 = M_tot/(PS::F64)N_ptcl;
        psys_gas[i].mass = psys_gas[i].mass0;
        // Rescale position
        // (1) First, we scale the computational box.
        PS::F64vec pos = ptcl[i];
        pos.x = (len_root_domain.x/len.x) * (pos.x - cen.x);
        pos.y = (len_root_domain.y/len.y) * (pos.y - cen.y);
        pos.z = (len_root_domain.z/len.z) * (pos.z - cen.z);
        // (2) Then, we stretch the particle distribution to give
        //     the required left/right states.
        if ((pos.x - lo_root_domain.x) < (frac_L * len_root_domain.x)) {
           const PS::F64 fac = xlen_L / (frac_L * len_root_domain.x);
           pos.x = fac * (pos.x - lo_root_domain.x) + lo_root_domain.x;
        } else {
           const PS::F64 fac = xlen_R / (frac_R * len_root_domain.x);
           pos.x = hi_root_domain.x - fac * (hi_root_domain.x - pos.x);
        }
        // (3) After that, we move some of the particle to be
        //     included in the specified region.
        while (pos.x < lo_root_domain.x) {
            pos.x += len_root_domain.x;
        }
        while (pos.x >= hi_root_domain.x) {
            pos.x -= len_root_domain.x;
        }
        while (pos.y < lo_root_domain.y) {
            pos.y += len_root_domain.y;
        }
        while (pos.y >= hi_root_domain.y) {
            pos.y -= len_root_domain.y;
        }
        while (pos.z < lo_root_domain.z) {
            pos.z += len_root_domain.z;
        }
        while (pos.z >= hi_root_domain.z) {
            pos.z -= len_root_domain.z;
        }
        psys_gas[i].pos = pos;
        // Set the other physical quantities
        psys_gas[i].id = i;
        psys_gas[i].vel = 0.0;
        psys_gas[i].h = h_ini;
        psys_gas[i].alpha = run_param::sph::alpha_AV_ini;
        if (psys_gas[i].pos.x <= 0.0) {
            psys_gas[i].eng = pres_L / (rho_L * (run_param::ism::gamma - 1.0));
        } else {
            psys_gas[i].eng = pres_R / (rho_R  * (run_param::ism::gamma - 1.0));
        }
        for (PS::S32 k = 0; k < CELibYield_Number; k++)
            psys_gas[i].mabn[0] = 0.0;
    }
    // Set boundary condition
    run_param::basic::bc = PS::BOUNDARY_CONDITION_PERIODIC_XYZ;
    // Set the parameters for gravity calculation
    // (Dummy values are set for the parameters that are not used actually)
    run_param::grav::soft::theta = 0.5;
    run_param::grav::soft::n_group_limit = 64;
    run_param::grav::soft::n_leaf_limit = 8;
    run_param::grav::soft::eps_dm = 1.0e-2 / run_param::unit::leng;
    run_param::grav::soft::eps_star = 1.0e-2 / run_param::unit::leng;
    run_param::grav::soft::eps_gas = 1.0e-2 / run_param::unit::leng;
    run_param::grav::soft::dt_fid = 0.01 / run_param::unit::time;
    run_param::grav::hard::eps = 1.0e-2 / run_param::unit::leng;
    run_param::grav::hard::eta = 0.1;
    run_param::grav::hard::eta_s = 1.0e-3;
    run_param::grav::hard::dt_limit = run_param::grav::soft::dt_fid / 8;
    // Set the parameters for SPH calculation
    run_param::sph::n_group_limit = 64;
    run_param::sph::n_leaf_limit = 8;
    run_param::sph::n_jp_limit = 4096;
    // Set I/O intervals
    run_param::io::dt_dump = 0.1 * t_end / run_param::unit::time;
    run_param::io::time_dump = run_param::io::dt_dump;
    run_param::io::dt_dump_rst = 0.1 * t_end / run_param::unit::time;
    run_param::io::time_dump_rst = run_param::io::dt_dump;
    // Set the end time of the simulation
    run_param::basic::time_end = 1.05 * t_end / run_param::unit::time;

    if (PS::Comm::getRank() == 0) {
        std::cout << "An initial condition for shock tube test is made." << std::endl;
        std::cout << "N_sph = " << N_ptcl << std::endl;
        std::cout << "h_ini = " << h_ini << std::endl;
    }

}
#endif // ASURA_FDPS_SHOCK_TUBE_TEST


#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_SURFACE_TENSION_TEST
void SurfaceTensionTestIC(PS::ParticleSystem<FP_dm>& psys_dm,
                          PS::ParticleSystem<FP_star>& psys_star,
                          PS::ParticleSystem<FP_gas>& psys_gas) {
    // Set the parameters of surface tension test
    constexpr PS::S32 dens_cont_param = 2; 
    constexpr PS::F64 dens_cont = SQ(dens_cont_param);
    constexpr PS::F64 rho_low = 1.0;
    constexpr PS::F64 rho_high = dens_cont * rho_low;
    constexpr PS::F64 pres = 2.5;
    constexpr PS::S32 N_grid = 64; // the number of grids per side in the low-dense region
    constexpr PS::F64 dx = 1.0 / N_grid; // the cell size in the low-dense region
    constexpr PS::F64 dx_s = dx / dens_cont_param; // the cell size in the high-dense region
    // Set the computational domain and the size of the dense region
    run_param::basic::pos_root_domain = PS::F64ort(PS::F64vec(0.0),
                                                   PS::F64vec(1.0));
    PS::F64ort pos_dense_rgn;
    pos_dense_rgn.low_.x = (N_grid/4) * dx;
    pos_dense_rgn.low_.y = (N_grid/4) * dx;
    pos_dense_rgn.high_.x = ((3*N_grid)/4) * dx;
    pos_dense_rgn.high_.y = ((3*N_grid)/4) * dx;
    const PS::F64 Area_high = (pos_dense_rgn.high_.x - pos_dense_rgn.low_.x)
                            * (pos_dense_rgn.high_.y - pos_dense_rgn.low_.y);
    const PS::F64 Area_low  = 1.0 - Area_high;
    // Choose how do we create the initial particle distribution 
    constexpr bool use_prng = false;
    PS::S32 N_ptcl {0};
    if (!use_prng) {
        // In this case, we place particles on grids.
        // First, we count the necessary number of particles
        for (PS::S32 i = 0; i < N_grid; i++) {
            const PS::F64 x = dx * (i + 0.5);
            for (PS::S32 j = 0; j < N_grid; j++) {
                const PS::F64 y = dx * (j + 0.5);
                PS::F64vec pos = PS::F64vec(x,y);
                if (pos_dense_rgn.contains(pos)) N_ptcl += SQ(dens_cont_param);
                else N_ptcl++;
            }
        }
        // Allocate memory for particles
        psys_dm.setNumberOfParticleLocal(0);
        psys_star.setNumberOfParticleLocal(0); 
        psys_gas.setNumberOfParticleLocal(N_ptcl);
        // Then, we place particles
        //constexpr PS::F64 eps = 0.0;
        constexpr PS::F64 eps = 1.0e-8;
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        PS::S32 idx {0};
        for (PS::S32 i = 0; i < N_grid; i++) {
            const PS::F64 x = dx * (i + 0.5);
            for (PS::S32 j = 0; j < N_grid; j++) {
                const PS::F64 y = dx * (j + 0.5);
                PS::F64vec pos = PS::F64vec(x,y);
                if (pos_dense_rgn.contains(pos)) {
                    for (PS::S32 ii = -dens_cont_param/2; ii < dens_cont_param/2; ii++) {
                        pos.x = x + dx_s * (ii + 0.5);
                        for (PS::S32 jj = -dens_cont_param/2; jj < dens_cont_param/2; jj++) {
                            pos.y = y + dx_s * (jj + 0.5);
                            PS::F64 p;
                            p = dist(run_param::prng::mt);
                            pos.x += eps * (2.0 * p - 1.0);
                            p = dist(run_param::prng::mt);
                            pos.y += eps * (2.0 * p - 1.0);
                            psys_gas[idx].pos = pos;
                            idx++;
                        }
                    }
                } else {
                    PS::F64 p;
                    p = dist(run_param::prng::mt);
                    pos.x += eps * (2.0 * p - 1.0);
                    p = dist(run_param::prng::mt);
                    pos.y += eps * (2.0 * p - 1.0);
                    psys_gas[idx].pos = pos;
                    idx++;
                }
            }
        }
    } else {
        // In this case, we make the initial particle distribution
        // using the pseudo-random number generator.
        // Determine the number of particles and allocate memory for particles
        N_ptcl = (1<<16);
        psys_dm.setNumberOfParticleLocal(0);
        psys_star.setNumberOfParticleLocal(0); 
        psys_gas.setNumberOfParticleLocal(N_ptcl);
        // Then, we place particles
        const PS::F64 prob_high = rho_high * Area_high
                                / (rho_high * Area_high + rho_low * Area_low);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (PS::S32 i = 0; i < N_ptcl; i++) {
            PS::F64vec pos;
            const PS::F64 p = dist(run_param::prng::mt);
            if (p <= prob_high) {
                // In this case, a particle belongs to the high dense region
                PS::F64 f;
                f = dist(run_param::prng::mt);
                pos.x = pos_dense_rgn.low_.x
                      + (pos_dense_rgn.high_.x - pos_dense_rgn.low_.x) * f;
                f = dist(run_param::prng::mt);
                pos.y = pos_dense_rgn.low_.y
                      + (pos_dense_rgn.high_.y - pos_dense_rgn.low_.y) * f;
            } else {
                for(;;) {
                    PS::F64 f;
                    f = dist(run_param::prng::mt);
                    pos.x = run_param::basic::pos_root_domain.low_.x
                          + (run_param::basic::pos_root_domain.high_.x
                            -run_param::basic::pos_root_domain.low_.x) * f;
                    f = dist(run_param::prng::mt);
                    pos.y = run_param::basic::pos_root_domain.low_.y
                          + (run_param::basic::pos_root_domain.high_.y
                            -run_param::basic::pos_root_domain.low_.y) * f;
                    if ((pos.x <  pos_dense_rgn.low_.x) ||
                        (pos.x >= pos_dense_rgn.high_.x) ||
                        (pos.y <  pos_dense_rgn.low_.y) ||
                        (pos.y >= pos_dense_rgn.high_.y)) break;
                }
            }
            psys_gas[i].pos = pos;
        }
        
    }
    // Set physical quantities
    run_param::unit::mass = 1.0;
    run_param::unit::leng = 1.0;
    run_param::unit::time = 1.0;
    const PS::F64 total_mass = (rho_high * Area_high + rho_low * Area_low);
    const PS::F64 h_ini = 3.0 * run_param::unit::leng / std::sqrt(N_ptcl);
    for (PS::S32 i = 0; i < N_ptcl; i++) {
        psys_gas[i].id = i;
        psys_gas[i].mass0 = total_mass / N_ptcl;
        psys_gas[i].mass = psys_gas[i].mass0;
        psys_gas[i].vel = 0.0;
        psys_gas[i].h = h_ini;
        psys_gas[i].alpha = run_param::sph::alpha_AV_ini;
        if (pos_dense_rgn.contains(psys_gas[i].pos)) {
            psys_gas[i].eng = pres / (rho_high * (run_param::ism::gamma - 1.0));
        } else {
            psys_gas[i].eng = pres / (rho_low  * (run_param::ism::gamma - 1.0));
        }
        for (PS::S32 k = 0; k < CELibYield_Number; k++)
            psys_gas[i].mabn[0] = 0.0;
    }
    // Set boundary condition
    run_param::basic::bc = PS::BOUNDARY_CONDITION_PERIODIC_XY;
    // Set the parameters for gravity calculation
    // (Dummy values are set for the parameters that are not used actually)
    run_param::grav::soft::theta = 0.5;
    run_param::grav::soft::n_group_limit = 64;
    run_param::grav::soft::n_leaf_limit = 8;
    run_param::grav::soft::eps_dm = 1.0e-2 / run_param::unit::leng;
    run_param::grav::soft::eps_star = 1.0e-2 / run_param::unit::leng;
    run_param::grav::soft::eps_gas = 1.0e-2 / run_param::unit::leng;
    run_param::grav::soft::dt_fid = 0.01 / run_param::unit::time;
    run_param::grav::hard::eps = 1.0e-2 / run_param::unit::leng;
    run_param::grav::hard::eta = 0.1;
    run_param::grav::hard::eta_s = 1.0e-3;
    run_param::grav::hard::dt_limit = run_param::grav::soft::dt_fid / 8;
    // Set the parameters for SPH calculation
    run_param::sph::n_group_limit = 64;
    run_param::sph::n_leaf_limit = 8;
    run_param::sph::n_jp_limit = 4096;
    // Set I/O intervals
    run_param::io::dt_dump = 0.1 / run_param::unit::time;
    run_param::io::time_dump = run_param::io::dt_dump;
    run_param::io::dt_dump_rst = 0.1 / run_param::unit::time;
    run_param::io::time_dump_rst = run_param::io::dt_dump;
    // Set the end time of the simulation
    run_param::basic::time_end = 8.0 / run_param::unit::time;

    if (PS::Comm::getRank() == 0) {
        std::cout << "An initial condition for surface tenstion test is made." << std::endl;
        std::cout << "N_sph = " << N_ptcl << std::endl;
        std::cout << "h_ini = " << h_ini << std::endl;
    }
}
#endif // ASURA_FDPS_SURFACE_TENSION_TEST


#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_POINT_EXPLOSION_TEST
void PointExplosionTestIC(PS::ParticleSystem<FP_dm>& psys_dm,
                          PS::ParticleSystem<FP_star>& psys_star,
                          PS::ParticleSystem<FP_gas>& psys_gas) {
    // Set the parameters of point explosion test
    constexpr PS::F64 M_tot = 1.0;
    constexpr PS::F64 E_inj = 1.0;
    constexpr PS::F64 t_end = 0.05;
    // Set the simulation unit
    run_param::unit::mass = 1.0;
    run_param::unit::leng = 1.0;
    run_param::unit::time = 1.0;
    // Set the computational domain and the size of the dense region
    PS::F64ort pos_root_domain;
    pos_root_domain.low_.x  = - 0.5;         
    pos_root_domain.high_.x =   0.5;
    pos_root_domain.low_.y  = - 0.5;
    pos_root_domain.high_.y =   0.5;
    pos_root_domain.low_.z  = - 0.5;
    pos_root_domain.high_.z =   0.5;
    run_param::basic::pos_root_domain = pos_root_domain;
    // Read glass data
    PS::F64ort pos_glass_data;
    std::vector<PS::F64vec> ptcl;
    const std::string filename = "result/glass_data.txt";
    std::ifstream ifs;
    ifs.open(filename.c_str(), std::ios::in | std::ios::binary);
    if (ifs) {
        while(true) {
            PS::F64vec pos;
            ifs.read((char *)&pos, sizeof(PS::F64vec));
            if (ifs.eof()) break;
            ptcl.push_back(pos);
            pos_glass_data.merge(pos);
        }
    } else {
        std::cout << "Could not open file " << filename << std::endl;
        assert(false);
    }
    ifs.close();
    const PS::S64 N_ptcl = ptcl.size();
    std::cout << N_ptcl << " particles are read." << std::endl;
    // Allocate memory for particles
    psys_dm.setNumberOfParticleLocal(0);
    psys_star.setNumberOfParticleLocal(0);
    psys_gas.setNumberOfParticleLocal(N_ptcl);
    // Rescale the position and the mass of the particles
    const PS::F64vec lo = pos_glass_data.low_;
    const PS::F64vec hi = pos_glass_data.high_;
    const PS::F64vec cen = pos_glass_data.getCenter();
    const PS::F64vec len = pos_glass_data.getFullLength();
    const PS::F64vec lo_root_domain = pos_root_domain.low_;
    const PS::F64vec hi_root_domain = pos_root_domain.high_;
    const PS::F64vec cen_root_domain = pos_root_domain.getCenter();
    const PS::F64vec len_root_domain = pos_root_domain.getFullLength();
    const PS::F64 fac = std::max(len_root_domain.x, std::max(len_root_domain.y, len_root_domain.z))
                      / std::max(len.x, std::max(len.y, len.z));
    for (PS::S32 i = 0; i < N_ptcl; i++) {
        // Rescale mass
        psys_gas[i].mass0 = M_tot/(PS::F64)N_ptcl;
        psys_gas[i].mass = psys_gas[i].mass0;
        // Rescale position
        // (1) First, we scale the computational box.
        PS::F64vec pos = ptcl[i];
        pos.x = fac * (pos.x - cen.x);
        pos.y = fac * (pos.y - cen.y);
        pos.z = fac * (pos.z - cen.z);
        // (2) After that, we move some of the particle to be
        //     included in the specified region.
        while (pos.x < lo_root_domain.x) {
            pos.x += len_root_domain.x;
        }
        while (pos.x >= hi_root_domain.x) {
            pos.x -= len_root_domain.x;
        }
        while (pos.y < lo_root_domain.y) {
            pos.y += len_root_domain.y;
        }
        while (pos.y >= hi_root_domain.y) {
            pos.y -= len_root_domain.y;
        }
        while (pos.z < lo_root_domain.z) {
            pos.z += len_root_domain.z;
        }
        while (pos.z >= hi_root_domain.z) {
            pos.z -= len_root_domain.z;
        }
        psys_gas[i].pos = pos;
    }
    // Determine the energy-input radius
    const PS::S32 N_ngb_inj = 32;
    PS::F64 r_inj = 3.0 / std::cbrt((PS::F64)N_ptcl);
    PS::F64 r_L = 0.0;
    PS::F64 r_U = 10.0 * r_inj;
    for (;;) {
        PS::S32 N_ngb = 0;
        for (PS::S32 i = 0; i < N_ptcl; i++) {
            const PS::F64vec pos = psys_gas[i].pos;
            const PS::F64 r = std::sqrt(pos * pos);
            if (r < r_inj) N_ngb++;
        }
        // Termination condition
        if (N_ngb == N_ngb_inj) break;
        // Update r_L, r_U
        if ((N_ngb < N_ngb_inj) && (r_L < r_inj)) r_L = r_inj;
        if ((N_ngb > N_ngb_inj) && (r_U > r_inj)) r_U = r_inj;
        // Update r_inj
        const PS::F64 s = std::cbrt((PS::F64)N_ngb_inj/(PS::F64)std::max(N_ngb,1));
        r_inj = 0.5 * r_inj * (1.0 + s); // Hernquist & Katz (1989)
        if ((r_inj <= r_L) || (r_inj == r_U)) {
           r_inj = 0.5 * (r_L + r_U); // switch to bisection search.
        } else if (r_U < r_inj) {
           r_inj = r_U;
        }
    }
    // Set the other physical quantities
    const PS::F64 h_ini = 3.0 / std::cbrt(N_ptcl);
    PS::F64 wtot {0.0};
    for (PS::S32 i = 0; i < N_ptcl; i++) {
        psys_gas[i].id = i;
        psys_gas[i].vel = 0.0;
        psys_gas[i].h = h_ini;
        psys_gas[i].alpha = run_param::sph::alpha_AV_ini;
        const PS::F64 r = std::sqrt(SQ(psys_gas[i].pos));
        const PS::F64 wj = W(r, r_inj);
        wtot += wj;
        psys_gas[i].eng = E_inj * wj;
        for (PS::S32 k = 0; k < CELibYield_Number; k++)
            psys_gas[i].mabn[0] = 0.0;
    }
    for (PS::S32 i = 0; i < N_ptcl; i++) {
        psys_gas[i].eng /= wtot; // To normalize
        psys_gas[i].eng += 1.0e-6; // Add tiny thermal energy
        psys_gas[i].eng /= psys_gas[i].mass; // To convert to specific energy
    }
    // Set boundary condition
    run_param::basic::bc = PS::BOUNDARY_CONDITION_PERIODIC_XYZ;
    // Set the parameters for gravity calculation
    // (Dummy values are set for the parameters that are not used actually)
    run_param::grav::soft::theta = 0.5;
    run_param::grav::soft::n_group_limit = 64;
    run_param::grav::soft::n_leaf_limit = 8;
    run_param::grav::soft::eps_dm = 1.0e-2 / run_param::unit::leng;
    run_param::grav::soft::eps_star = 1.0e-2 / run_param::unit::leng;
    run_param::grav::soft::eps_gas = 1.0e-2 / run_param::unit::leng;
    run_param::grav::soft::dt_fid = 0.01 / run_param::unit::time;
    run_param::grav::hard::eps = 1.0e-2 / run_param::unit::leng;
    run_param::grav::hard::eta = 0.1;
    run_param::grav::hard::eta_s = 1.0e-3;
    run_param::grav::hard::dt_limit = run_param::grav::soft::dt_fid / 8;
    // Set the parameters for SPH calculation
    run_param::sph::n_group_limit = 64;
    run_param::sph::n_leaf_limit = 8;
    run_param::sph::n_jp_limit = 4096;
    // Set I/O intervals
    run_param::io::dt_dump = 0.1 * t_end / run_param::unit::time;
    run_param::io::time_dump = run_param::io::dt_dump;
    run_param::io::dt_dump_rst = 0.1 * t_end / run_param::unit::time;
    run_param::io::time_dump_rst = run_param::io::dt_dump;
    // Set the end time of the simulation
    run_param::basic::time_end = 1.05 * t_end / run_param::unit::time;

    if (PS::Comm::getRank() == 0) {
        std::cout << "An initial condition for point explosion test is made." << std::endl;
        std::cout << "N_sph = " << N_ptcl << std::endl;
        std::cout << "M_tot = " << M_tot << std::endl;
        std::cout << "E_inj = " << E_inj << std::endl;
        std::cout << "N_ngb_inj = " << N_ngb_inj << std::endl;
        std::cout << "r_inj = " << r_inj << std::endl;
    }

}
#endif // ASURA_FDPS_POINT_EXPLOSION_TEST
