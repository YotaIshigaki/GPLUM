/* C++ headers */
#include "common.h"
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "macro_defs.h"
#include "mathematical_constants.h"
#include "physical_constants.h"
#include "run_parameters.h"
#include "user_defined.h"
#include "ic.h"

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
                   PS::ParticleSystem<FP_nbody>& psys) {
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
#ifdef READ_DATA_WITH_BYTESWAP
        std::cout << "nbodies = " << header.nbodies << std::endl;
        header.nbodies = byteswap(header.nbodies);
        std::cout << "nbodies = " << header.nbodies << std::endl;
#endif
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
#ifdef READ_DATA_WITH_BYTESWAP
        psys[i].mass  = byteswap(ptcl[i].mass);
        psys[i].pos.x = byteswap(ptcl[i].pos[0]);
        psys[i].pos.y = byteswap(ptcl[i].pos[1]);
        psys[i].pos.z = byteswap(ptcl[i].pos[2]);
        psys[i].vel.x = byteswap(ptcl[i].vel[0]);
        psys[i].vel.y = byteswap(ptcl[i].vel[1]);
        psys[i].vel.z = byteswap(ptcl[i].vel[2]);
#else
        psys[i].mass  = ptcl[i].mass;
        psys[i].pos.x = ptcl[i].pos[0];
        psys[i].pos.y = ptcl[i].pos[1];
        psys[i].pos.z = ptcl[i].pos[2];
        psys[i].vel.x = ptcl[i].vel[0];
        psys[i].vel.y = ptcl[i].vel[1];
        psys[i].vel.z = ptcl[i].vel[2];
#endif
    }

}

void GalaxyIC(PS::ParticleSystem<FP_nbody>& psys_nbody,
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
    readTipsyFile(filename, psys_nbody);
    const PS::S64 N_nbody = psys_nbody.getNumberOfParticleLocal();
    for (PS::S32 i = 0; i < N_nbody; i++) {
        psys_nbody[i].id = i;
        psys_nbody[i].acc  = 0.0;
        psys_nbody[i].pot  = 0.0;
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
        psys_gas[id].eng = (phys_const::kBoltz * temp)/((run_param::ism::gamma - 1.0) * mu * phys_const::Mhydrogen);
        psys_gas[id].smth = std::pow(Rt*Rt*zt, 1.0/3.0) * std::pow((PS::F64) run_param::sim::N_neighbor / (PS::F64)N_sph, 1.0/3.0);
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
        psys_nbody[i].mass *= magi_unit_mass;
        psys_nbody[i].pos  *= magi_unit_leng;
        psys_nbody[i].vel  *= magi_unit_velc;
        run_param::unit::mass += psys_nbody[i].mass;
        const PS::F64 r = std::sqrt(psys_nbody[i].pos * psys_nbody[i].pos);
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
        psys_nbody[i].mass /= run_param::unit::mass;
        psys_nbody[i].pos  /= run_param::unit::leng;
        psys_nbody[i].vel  /= run_param::unit::velc;
    }
    for (PS::S32 i = 0; i < N_sph; i++) {
        psys_gas[i].mass0 /= run_param::unit::mass;
        psys_gas[i].mass  /= run_param::unit::mass;
        psys_gas[i].pos   /= run_param::unit::leng;
        psys_gas[i].vel   /= run_param::unit::velc;
        psys_gas[i].smth  /= run_param::unit::leng;
        psys_gas[i].eng   /= run_param::unit::spen;
    }
    // Set boundary condition
    run_param::basic::bc = PS::BOUNDARY_CONDITION_OPEN;
    // Set the other parameters
    run_param::sim::eps_grav = 1.0e-3;
    run_param::sim::mass_avg = psys_gas[0].mass0;
    // Set I/O intervals
    run_param::io::dt_dump = 0.001;
    run_param::io::time_dump = run_param::io::dt_dump;
    run_param::io::dt_dump_rst = 0.001;
    run_param::io::time_dump_rst = run_param::io::dt_dump;
    // Set the end time of the simulation
    run_param::basic::time_end = 10.0;
    // Set maximum timestep
    run_param::sim::dt_max = 1.0e-4;

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

