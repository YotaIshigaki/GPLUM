#include <cassert>
#include <iostream>
#include <random>
#include <type_traits>
#include <vector>

#include "particle_simulator.hpp"

class forceNearEP {
  public:
    PS::S32 nb = 0;
    void clear() {}
};

/*FullParticle class for FDPS*/
class SphereFP {
  public:
    PS::S32 gid;
    PS::F64vec3 pos;
    PS::F64 RSearch;
    PS::S32 nb = 0;

    PS::F64vec getPos() const { return this->pos; }

    void copyFromForce(const forceNearEP &force) { nb = force.nb; }

    void setPos(const PS::F64vec3 &pos_) { this->pos = pos_; }

    void writeAscii(FILE *fp) const {
        const double radius = 1;
        fprintf(fp, "S\t%8d\t%.17g\t%.17g\t%.17g\t%.17g\n", this->gid, radius, this->pos[0], this->pos[1],
                this->pos[2]);
    }
};

class SphereNearForceEP {
  public:
    PS::S32 gid;
    PS::F64vec3 pos;
    PS::F64 RSearch;

    PS::F64vec getPos() const { return this->pos; }

    void setPos(const PS::F64vec3 &pos) { this->pos = pos; }

    PS::F64 getRSearch() const { return this->RSearch; }

    void copyFromFP(const SphereFP &fp) {
        gid = fp.gid;
        pos = fp.pos;
        RSearch = fp.RSearch;
    }
};

class SphereNearForceEPJ {
  public:
    PS::S32 gid;
    PS::F64vec3 pos;
    PS::F64 RSearch;

    PS::F64vec getPos() const { return this->pos; }

    void setPos(const PS::F64vec3 &pos) { this->pos = pos; }

    PS::F64 getRSearch() const { return this->RSearch; }

    void copyFromFP(const SphereFP &fp) {
        gid = fp.gid;
        pos = fp.pos;
        RSearch = fp.RSearch;
    }
};

class CalcNearForceEPIJ {
  public:
    // constructor
    CalcNearForceEPIJ() {}

    // copy constructor
    CalcNearForceEPIJ(const CalcNearForceEPIJ &obj) {}

    void operator()(const SphereNearForceEP *const ep_i, const PS::S32 Nip, const SphereNearForceEPJ *const ep_j,
                    const PS::S32 Njp, forceNearEP *const force) {

        for (PS::S32 i = 0; i < Nip; ++i) {
            force[i].clear();
            //auto &sphereI = ep_i[i];
            const SphereNearForceEP *ip = &ep_i[i];
            for (PS::S32 j = 0; j < Njp; ++j) {
                //auto &sphereJ = ep_j[j];
                const SphereNearForceEPJ *jp = &ep_j[j];
                force[i].nb++;
            }
        }
    }
};

class HogeWithgetRsearch
{
public:
   PS::F64 rsearchval;
   //PS::F64 getRSearch() const { return rsearchval; }
   PS::F64 getRSearch() { return rsearchval; }
};
class Hoge{};

int main(int argc, char **argv) {
    PS::Initialize(argc, argv);

#if 0
    class HogeWithgetRsearch p1;
    class Hoge p2;
    p1.rsearchval = 1.0;
    std::cout << PS::GetRSearch(p1) << std::endl;
    std::cout << PS::GetRSearch(p2) << std::endl;


    class SphereNearForceEP ptcl;
    ptcl.RSearch = 1.0;

    std::cout << PS::GetRSearch(ptcl) << std::endl;
    PS::Finalize();
    std::exit(0);
#endif

    PS::ParticleSystem<SphereFP> systemSP;
    systemSP.initialize();
    PS::DomainInfo dinfo;

    PS::TreeForForceShort<forceNearEP, SphereNearForceEP, SphereNearForceEPJ>::Symmetry treeSymmetry;
    PS::TreeForForceShort<forceNearEP, SphereNearForceEP, SphereNearForceEPJ>::Scatter treeScatter;
    PS::TreeForForceShort<forceNearEP, SphereNearForceEP, SphereNearForceEPJ>::Gather treeGather;
    CalcNearForceEPIJ calcNearForceFtr;

    //const int nLocal = 1000;
    const int nLocal = 200000;
    //const int nLocal = 2000000;
    //const int nLocal = 10000000;
    systemSP.setNumberOfParticleLocal(nLocal);
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        systemSP[i].pos[0] = sin(i) * 100;
        systemSP[i].pos[1] = cos(i) * 100;
        systemSP[i].pos[2] = sin(cos(i)) * 100;
        systemSP[i].RSearch = 10.0;
    }

    std::ofstream output_file;
    output_file.open("IC.txt",std::ios::trunc);
    for (int i = 0; i < nLocal; i++)
        output_file << systemSP[i].pos << std::endl;
    output_file.close();

    treeSymmetry.initialize(nLocal * 4);
    treeScatter.initialize(nLocal * 4);
    treeGather.initialize(nLocal * 4);

    PS::Comm::barrier();

    dinfo.initialize();
    //dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
    //dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_X);
    //dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_Y);
    //dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_Z);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XY);
    //dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_YZ);
    //dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XZ);
    //dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain(PS::F64vec3(0, 0, 0), PS::F64vec3(100, 100, 100));
    systemSP.adjustPositionIntoRootDomain(dinfo);

    PS::Comm::barrier();
    dinfo.decomposeDomainAll(systemSP);
    systemSP.exchangeParticle(dinfo);
    assert(systemSP.getNumberOfParticleLocal() > 0);

    PS::Comm::barrier();
    treeSymmetry.calcForceAll(calcNearForceFtr, systemSP, dinfo);
    printf("TreeSymmetry complete\n");

    PS::Comm::barrier();
    treeScatter.calcForceAll(calcNearForceFtr, systemSP, dinfo);
    printf("TreeScatter complete\n");

    PS::Comm::barrier();
    treeGather.calcForceAll(calcNearForceFtr, systemSP, dinfo);
    printf("TreeGather complete\n");

    PS::Comm::barrier();
    PS::Finalize();
    return 0;
}
