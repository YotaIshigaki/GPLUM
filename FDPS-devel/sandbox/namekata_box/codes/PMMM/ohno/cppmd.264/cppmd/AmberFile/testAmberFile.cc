#include "AmberFile.h"
#include <cstdio>

using namespace std;

void printParticle(const vector<Particle>& particles) {
  printf("PPPPPPPPPPPPPP %d particles\n", particles.size());
  for (int i = 0; i < particles.size(); ++i) {
    const Particle& p = particles[i];
    printf(
        "%d\tP:(% f,% f,% f), V:(% f,% f,% f), F:(% f,% f,% f), C:% f, M:%f, T:%d, I:%d\n",
        i, p.position.x, p.position.y, p.position.z,
        p.velocity.x, p.velocity.y, p.velocity.z,
        p.force.x, p.force.y, p.force.z,
        p.charge, p.mass, p.atomtype, p.atomid);
  }
}

void printBond(const vector<CovalentBondInfo::Bond>& info,
               const vector<CovalentBondInfo::BondParameter>& param) {
  printf("BBBBBBBBBBBBBB %d bonds, %d params\n", info.size(), param.size());
  for (int i = 0; i < info.size(); ++i) {
    const int type = info[i].typeofbond;
    printf(
        "%d\tB:(% 5d,% 5d) [%d] %f, %f\n",
        i, info[i].id_of_atom[0], info[i].id_of_atom[1], info[i].typeofbond,
        param[type].force_constant, param[type].equilibrium_length);
  }
}

void printAngle(const vector<CovalentBondInfo::Angle>& info,
                const vector<CovalentBondInfo::AngleParameter>& param) {
  printf("AAAAAAAAAAAAAA %d angles, %d params\n", info.size(), param.size());
  for (int i = 0; i < info.size(); ++i) {
    const int type = info[i].typeofangle;
    printf(
        "%d\tA:(% 5d,% 5d,% 5d) [%d] %f, %f\n",
        i, info[i].id_of_atom[0], info[i].id_of_atom[1], info[i].id_of_atom[2],
        info[i].typeofangle,
        param[type].force_constant,
        param[type].equilibrium_angle);
  }
}

void printTorsion(const vector<CovalentBondInfo::Torsion>& info,
                  const vector<CovalentBondInfo::TorsionParameter>& param) {
  printf("TTTTTTTTTTTTTT %d torsions, %d params\n", info.size(), param.size());
  for (int i = 0; i < info.size(); ++i) {
    printf(
        "%d\tT:(% 5d,% 5d,% 5d,% 5d) [%d]\n",
        i, info[i].id_of_atom[0], info[i].id_of_atom[1], info[i].id_of_atom[2],
        info[i].id_of_atom[3], info[i].typeoftorsion);
    int type = info[i].typeoftorsion;
    do {
      printf( "\t%f, %f, %f, %f\n",
          param[type].force_constant, param[type].periodicity,
          param[type].phase);
    } while (param[type++].periodicity < 0);
  }
}

void printImproper(const vector<CovalentBondInfo::Improper>& info,
                  const vector<CovalentBondInfo::ImproperParameter>& param) {
  printf("IIIIIIIIIIIIII %d impropers, %d params\n", info.size(), param.size());
  for (int i = 0; i < info.size(); ++i) {
    printf(
        "%d\tI:(% 5d,% 5d,% 5d,% 5d) [%d]\n",
        i, info[i].id_of_atom[0], info[i].id_of_atom[1], info[i].id_of_atom[2],
        info[i].id_of_atom[3], info[i].typeofimproper);
    int type = info[i].typeofimproper;
    do {
      printf( "\t%f, %f, %f, %f\n",
          param[type].force_constant, param[type].periodicity,
          param[type].phase);
    } while (param[type++].periodicity < 0);
  }
}

void printBoxAndSkew(const SpaceVector<double>& box,
                     const SpaceVector<double>& skew) {
  printf("BSBSBSBSBSBSBS\n");
  printf("\tBOX:(% f,% f,% f)\n\tSKEW:(% f,% f, %f)\n",
         box.x, box.y, box.z, skew.x, skew.y, skew.z);
}

int main(int argc, char **argv) {
  if (argc < 3) {
    printf("usage: %s prmtop restrt\n", argv[0]);
    return 1;
  }

  string prmtop_file(argv[1]);
  string restrt_file(argv[2]);
  ParticleArray particles;
  CovalentBondList bond_list;
  CovalentBondParameterList param_list;
  SpaceVector<double> box;
  SpaceVector<double> skew;

  ReadAmberFile(prmtop_file, restrt_file, &particles, &bond_list, &param_list,
                &box, &skew);

  printParticle(particles);
  printBond(bond_list.bond, param_list.bond);
  printAngle(bond_list.angle, param_list.angle);
  printTorsion(bond_list.torsion, param_list.torsion);
  printImproper(bond_list.improper, param_list.improper);
  printBoxAndSkew(box, skew);

  return 0;
}
