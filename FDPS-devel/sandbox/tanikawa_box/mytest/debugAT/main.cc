#include <iostream>
#include <vector>
#include <cassert>

#include "myheader.hpp"

template <class T>
void makeInitialCondition(PS::ParticleSystem<T> & system)
{
    int nloc = system.getNumberOfParticleLocal();

    for(int i = 0; i < nloc; i++) {
        system[i].id     = i;
        system[i].pos[0] = frand() - 0.5;
        system[i].pos[1] = 2. * frand() - 1.;
        system[i].pos[2] = 16. * frand();
    }
    
    return;
}

static const PS::S32 ntgt0 = 292;
static const PS::S32 ntgt1 = 401;

int main(int argc, char **argv)
{
    PS::Initialize(argc, argv);

    PS::S32 rank = MPI::COMM_WORLD.Get_rank();
    PS::S32 nloc0 =  512 * (rank + 1);
    PS::S32 nloc1 =  800 * (rank + 1);
    PS::DomainInfo dinfo;
    PS::ParticleSystem<FP> system0, system1;

    srand(rank);

    dinfo.initialize();
    dinfo.setDomain(2, 2, 2);

    if(rank == 0) {
        char fname[1024];
        std::sprintf(fname, "out/domain_all.txt");
        std::FILE *fp = std::fopen(fname, "w");
        dinfo.getRootDomain(fp);
        std::fclose(fp);
    }

    PS::S32 nloc0_all, nloc1_all;
    MPI::COMM_WORLD.Allreduce(&nloc0, &nloc0_all, 1, MPI::INT, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(&nloc1, &nloc1_all, 1, MPI::INT, MPI::SUM);

    system0.initialize();
    system0.setAverageTargetNumberOfSampleParticlePerProcess(37);
    system0.createParticle(nloc0_all);
    system0.setNumberOfParticleLocal(nloc0);
    system1.initialize();
    system1.setAverageTargetNumberOfSampleParticlePerProcess(50);
    system1.createParticle(nloc1_all);
    system1.setNumberOfParticleLocal(nloc1);

    makeInitialCondition(system0);
    makeInitialCondition(system1);

    {
        char fname[1024];
        std::sprintf(fname, "out/system0i_%04d.txt", rank);
        std::FILE *fp;
        fp = std::fopen(fname, "w");
        for(int i = 0; i < system0.getNumberOfParticleLocal(); i++) {
            fprintf(fp, "%+e %+e %+e\n", system0[i].pos[0], system0[i].pos[1], system0[i].pos[2]);
        }    
        std::fclose(fp);

        std::sprintf(fname, "out/system1i_%04d.txt", rank);
        fp = std::fopen(fname, "w");
        for(int i = 0; i < system1.getNumberOfParticleLocal(); i++) {
            fprintf(fp, "%+e %+e %+e\n", system1[i].pos[0], system1[i].pos[1], system1[i].pos[2]);
        }    
        std::fclose(fp);
    }
    
    dinfo.collectSampleParticle(system0); // **********************
    dinfo.collectSampleParticle(system1, 1.0, false); // **********************

    {
        char fname[1024];
        std::sprintf(fname, "out/psmp_%04d.txt", rank);
        std::FILE *fp = std::fopen(fname, "w");
        dinfo.getSampleParticleLocal(fp);
        std::fclose(fp);
    }
    
    {
        char fname[1024];
        std::sprintf(fname, "out/system0f_%04d.txt", rank);
        std::FILE *fp;
        fp = std::fopen(fname, "w");
        for(int i = 0; i < system0.getNumberOfParticleLocal(); i++) {
            fprintf(fp, "%+e %+e %+e\n", system0[i].pos[0], system0[i].pos[1], system0[i].pos[2]);
        }    
        std::fclose(fp);

        std::sprintf(fname, "out/system1f_%04d.txt", rank);
        fp = std::fopen(fname, "w");
        for(int i = 0; i < system1.getNumberOfParticleLocal(); i++) {
            fprintf(fp, "%+e %+e %+e\n", system1[i].pos[0], system1[i].pos[1], system1[i].pos[2]);
        }    
        std::fclose(fp);
    }

    dinfo.decomposeDomain();

    if(rank == 0) {
        char fname[1024];
        std::sprintf(fname, "out/psmp_all.txt");
        std::FILE *fp = std::fopen(fname, "w");
        dinfo.getSampleParticleTotal(fp);
        std::fclose(fp);
    }

    {
        char fname[1024];
        std::sprintf(fname, "out/posdomain_%04d.txt", rank);
        std::FILE *fp = std::fopen(fname, "w");
        dinfo.getPosDomainTotal(fp);
        std::fclose(fp);
        
    }

    system0.exchangeParticle(dinfo);
    system1.exchangeParticle(dinfo);

    {
        char fname[1024];
        std::FILE *fp;
        std::sprintf(fname, "out/system0e_%04d.txt", rank);
        fp = std::fopen(fname, "w");
        for(int i = 0; i < system0.getNumberOfParticleLocal(); i++) {
            fprintf(fp, "%+e %+e %+e\n", system0[i].pos[0], system0[i].pos[1], system0[i].pos[2]);
        }    
        std::fclose(fp);
        std::sprintf(fname, "out/system1e_%04d.txt", rank);
        fp = std::fopen(fname, "w");
        for(int i = 0; i < system1.getNumberOfParticleLocal(); i++) {
            fprintf(fp, "%+e %+e %+e\n", system1[i].pos[0], system1[i].pos[1], system1[i].pos[2]);
        }    
        std::fclose(fp);
    }

/*    
    dinfo.collectSampleParticle(system0); // **********************
    dinfo.collectSampleParticle(system1, 1.0, false); // **********************

    dinfo.decomposeDomain();

    {
        char fname[1024];
        std::sprintf(fname, "out/posdomain_%04d_2.txt", rank);
        std::FILE *fp = std::fopen(fname, "w");
        dinfo.getPosDomainTotal(fp);
        std::fclose(fp);
        
    }
*/

    PS::Finalize();

    return 0;
}
