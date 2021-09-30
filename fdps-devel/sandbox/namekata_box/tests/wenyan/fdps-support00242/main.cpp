#include <cassert>
#include <iostream>
#include <type_traits>
#include <vector>

#include "particle_simulator.hpp"

class forceNearEP {
  public:
    void clear() {}
};

/*FullParticle class for FDPS*/
class SphereFP {
  public:
    PS::S32 gid;
    PS::F64vec3 pos;
    PS::F64 RSearch;

    PS::F64vec getPos() const { return this->pos; }

    void copyFromForce(const forceNearEP &force) {}

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
            auto &sphereI = ep_i[i];
            for (PS::S32 j = 0; j < Njp; ++j) {
                auto &sphereJ = ep_j[j];
            }
        }
    }
};

int main(int argc, char **argv) {

    //static_assert(std::is_trivially_destructible<PS::F64vec3>::value);
    //static_assert(std::is_trivially_copyable<PS::F64vec3>::value);
    //static_assert(std::is_trivially_copyable<PS::F64vec2>::value);

    PS::Initialize(argc, argv);

    PS::ParticleSystem<SphereFP> systemSP;
    systemSP.initialize();
    PS::DomainInfo dinfo;

    PS::TreeForForceShort<forceNearEP, SphereNearForceEP, SphereNearForceEPJ>::Symmetry treeSymmetry;
    PS::TreeForForceShort<forceNearEP, SphereNearForceEP, SphereNearForceEPJ>::Scatter treeScatter;
    CalcNearForceEPIJ calcNearForceFtr;

    if (PS::Comm::getRank() == 0) {
        systemSP.setNumberOfParticleLocal(5);
        // for the same particle, EPI and EPJ should have the same search radius
        systemSP[0].gid = 0;
        systemSP[0].RSearch = 8.22;
        systemSP[0].setPos(PS::F64vec3(35.2035, 0.83452, 31.2377));
        systemSP[1].gid = 1;
        systemSP[1].RSearch = 8.22;
        systemSP[1].setPos(PS::F64vec3(37.3486, 13.503, 24.2473));
        systemSP[2].gid = 2;
        systemSP[2].RSearch = 8.22;
        systemSP[2].setPos(PS::F64vec3(49.1056, 21.501, 31.0099));
        systemSP[3].gid = 3;
        systemSP[3].RSearch = 8.22;
        systemSP[3].setPos(PS::F64vec3(20.7771, 3.96282, 36.4081));
        systemSP[4].gid = 4;
        systemSP[4].RSearch = 8.22;
        systemSP[4].setPos(PS::F64vec3(40.3761, 3.27552, 3.38661));
        for (PS::S32 i = 0; i < systemSP.getNumberOfParticleLocal(); i++) {
            std::cout << systemSP[i].pos << std::endl;
        }
    } else {
        systemSP.setNumberOfParticleLocal(0);
    }
    PS::Comm::barrier();

    dinfo.initialize();
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain(PS::F64vec3(0, 0, 0), PS::F64vec3(100, 100, 100));
    systemSP.adjustPositionIntoRootDomain(dinfo);

    PS::Comm::barrier();
    dinfo.decomposeDomainAll(systemSP);
    systemSP.exchangeParticle(dinfo);
    assert(systemSP.getNumberOfParticleLocal() > 0);

#if 1
    if (PS::Comm::getRank() == 0) std::cout << "1st force calc." << std::endl;
    PS::Comm::barrier();
    treeSymmetry.initialize(20);
    treeSymmetry.calcForceAllAndWriteBack(calcNearForceFtr, systemSP, dinfo);
#endif

    if (PS::Comm::getRank() == 0) std::cout << "2nd force calc." << std::endl;
    PS::Comm::barrier();
    treeScatter.initialize(20);
    treeScatter.calcForceAllAndWriteBack(calcNearForceFtr, systemSP, dinfo);

    PS::Comm::barrier();
    PS::Finalize();
    return 0;
}
