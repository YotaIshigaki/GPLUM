#include <cassert>
#include <iostream>
#include <vector>

#include "particle_simulator.hpp"

class forceNearEP {
public:
  std::vector<PS::F64vec3> nb; // the number of neighbors found
  void clear() { nb.clear(); }
};

/*FullParticle class for FDPS*/
class SphereFP {
public:
  PS::S32 gid;
  PS::F64vec3 pos; // Tubule center or protein center
  PS::F64 RSearch;
  std::vector<PS::F64vec3> nb; // the number of neighbors found

  PS::F64vec getPos() const { return this->pos; }

  void copyFromForce(const forceNearEP &force) {
    for(int i=0;i<force.nb.size();i++){
      nb.push_back(force.nb[i]);
    }
  }

  void setPos(const PS::F64vec3 &pos_) { this->pos = pos_; }

  void writeAscii(FILE *fp) const {
    const double radius = 1;
    fprintf(fp, "S\t%8d\t%.17g\t%.17g\t%.17g\t%.17g\n", this->gid, radius,
            this->pos[0], this->pos[1], this->pos[2]);
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
    RSearch = fp.RSearch * 10; // test for EPI/EPJ rsearch differ greatly;
  }
};

class CalcNearForceEPIJ {
public:
  // constructor
  CalcNearForceEPIJ() {}

  // copy constructor
  CalcNearForceEPIJ(const CalcNearForceEPIJ &obj) {}

  void operator()(const SphereNearForceEP *const ep_i, const PS::S32 Nip,
                  const SphereNearForceEPJ *const ep_j, const PS::S32 Njp,
                  forceNearEP *const force) {
    for (PS::S32 i = 0; i < Nip; ++i) {
      force[i].clear();
      auto &sphereI = ep_i[i];
      for (PS::S32 j = 0; j < Njp; ++j) {
        auto &sphereJ = ep_j[j];
        force[i].nb.push_back(sphereJ.pos);
      }
    }
  }
};

int main(int argc, char **argv) {
  PS::Initialize(argc, argv);

  PS::ParticleSystem<SphereFP> systemSP;
  systemSP.initialize();
  PS::DomainInfo dinfo;

  PS::TreeForForceShort<forceNearEP, SphereNearForceEP,
                        SphereNearForceEPJ>::Symmetry treeSymmetry;
  PS::TreeForForceShort<forceNearEP, SphereNearForceEP,
                        SphereNearForceEPJ>::Scatter treeScatter;
  CalcNearForceEPIJ calcNearForceFtr;

  systemSP.setNumberOfParticleLocal(2); // 1 particle only
  systemSP[0].pos.x = 1;
  systemSP[0].pos.y = 10;
  systemSP[0].pos.z = 10;
  systemSP[0].gid = 0;
  systemSP[0].RSearch = 1.0; // EPI Rsearch = 1.0, EPJ Rsearch = 10.0
  systemSP[0].nb.clear();
  systemSP[1].pos.x = 15;
  systemSP[1].pos.y = 10;
  systemSP[1].pos.z = 10;
  systemSP[1].gid = 1;
  systemSP[1].RSearch = 1.0; // EPI Rsearch = 1.0, EPJ Rsearch = 10.0
  systemSP[1].nb.clear();

  dinfo.initialize(0.5);
  dinfo.setBoundaryCondition(
      PS::BOUNDARY_CONDITION_PERIODIC_X); // PERIODIC_X
  dinfo.setPosRootDomain(PS::F64vec3(0, 0, 0), PS::F64vec3(20, 20, 20));
  systemSP.adjustPositionIntoRootDomain(dinfo);
  dinfo.decomposeDomainAll(systemSP);
  systemSP.exchangeParticle(dinfo);

  treeSymmetry.initialize(10000);
  treeScatter.initialize(10000);

  treeSymmetry.calcForceAllAndWriteBack(calcNearForceFtr, systemSP, dinfo);
  int nbfound = systemSP[0].nb.size();
  std::cout << nbfound << "neighbors found with Symmetry search for Particle 0" << std::endl;
  for (int i = 0; i < nbfound; i++) {
    std::cout << systemSP[0].nb[i] << std::endl;
  }

  systemSP[0].nb.clear();
  systemSP[1].nb.clear();

  treeScatter.calcForceAllAndWriteBack(calcNearForceFtr, systemSP, dinfo);
  nbfound = systemSP[0].nb.size();
  std::cout << nbfound << "neighbors found with Scatter search for Particle 0" << std::endl;
  for (int i = 0; i < nbfound; i++) {
    std::cout << systemSP[0].nb[i] << std::endl;
  }


  PS::Finalize();
  return 0;
}
