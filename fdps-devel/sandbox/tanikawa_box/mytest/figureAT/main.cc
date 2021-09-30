#include <iostream>
#include <vector>
#include <cassert>

#include "myheader.hpp"

template <class T>
void loadInitialCondition(std::FILE *fp,
                          PS::ParticleSystem<T> & system)
{
    PS::S32 itmp;
    PS::F64 dtmp;

    PS::S32 nloc;

    std::fscanf(fp, "%lf", &dtmp);
    std::fscanf(fp, "%d", &nloc);
    
    system.setNumberOfParticleLocal(nloc);
    
    for(int i = 0; i < nloc; i++) {
        system[i].id     = i;
        std::fscanf(fp, "%d", &itmp);
        std::fscanf(fp, "%lf", &(system[i].mass));
        std::fscanf(fp, "%lf%lf%lf", &(system[i].pos[0]), &(system[i].pos[1]), &(system[i].pos[2]));
        std::fscanf(fp, "%lf%lf%lf", &(system[i].vel[0]), &(system[i].vel[1]), &(system[i].vel[2]));
    }

    return;
}

static PS::S32 nsampave = 10000;

int main(int argc, char **argv)
{
    PS::Initialize(argc, argv);

    PS::S32 rank = MPI::COMM_WORLD.Get_rank();
    PS::DomainInfo dinfo;
    PS::ParticleSystem<FP> system;

    srand(rank);

    dinfo.initialize();
    dinfo.setDomain(7, 6, 1);

    system.initialize();
    system.setAverageTargetNumberOfSampleParticlePerProcess(nsampave);
    system.setNumberOfParticleLocal(0);

    system.readParticleAscii(argv[1]);

    dinfo.decomposeDomainAll(system);
    system.exchangeParticle(dinfo);

    if(rank == 0) {
        char fname[1024];
        std::sprintf(fname, "out/domain_all.txt");
        std::FILE *fp = std::fopen(fname, "w");
        dinfo.getPosDomainTotal(fp);
        std::fclose(fp);
    }

    {
        char fname[1024];
        std::FILE *fp;
        std::sprintf(fname, "out/systeme_%04d.txt", rank);
        fp = std::fopen(fname, "w");
        for(int i = 0; i < system.getNumberOfParticleLocal(); i++) {
            fprintf(fp, "%+e %+e %+e\n", system[i].pos[0], system[i].pos[1], system[i].pos[2]);
        }    
        std::fclose(fp);
    }

    PS::Finalize();

    return 0;
}
