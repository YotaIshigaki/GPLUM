#include <iostream>
#include <vector>
#include <cassert>

#ifdef LONG_CUTOFF
#include "myheader_cutoff.hpp"
#else
#include "myheader.hpp"
#endif

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
        system[i].nj = 0;
        system[i].njreal = 0;
    }

    return;
}

static PS::S32 nsampave = 100;
static PS::S32 nptclmax = 32768;
static PS::F32 half_len = 1e3;

int main(int argc, char **argv)
{
    PS::Initialize(argc, argv);

    PS::S32 rank = MPI::COMM_WORLD.Get_rank();
    PS::DomainInfo dinfo;
    PS::ParticleSystem<FP> system;
#ifdef LONG_CUTOFF
//    PS::TreeForForceLong<Force, EPI, EPJ, PS::MomentMonoPole, PS::SPJMonoPole>::WithCutoff treegrav;
    PS::TreeForForceLong<Force, EPI, EPJ, PS::MomentMonoPolePeriodic, PS::SPJMonoPolePeriodic>::WithCutoff treegrav;
#else
    PS::TreeForForceLong<Force, EPI, EPJ, Moment, SPJ>::Normal treegrav;
#endif

    srand(rank);

    dinfo.initialize();
//    dinfo.setDomain(1, 1, 1);
    dinfo.setDomain(2, 2, 2);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain(0., 1.);

    system.initialize();
    system.setAverageTargetNumberOfSampleParticlePerProcess(nsampave);
    system.createParticle(nptclmax);
    system.setNumberOfParticleLocal(0);

    if(rank == 0) {
        FILE *fp = fopen(argv[1], "r");
        loadInitialCondition(fp, system);
        fclose(fp);
    }

    dinfo.collectSampleParticle(system);
    dinfo.decomposeDomain();

    system.exchangeParticle(dinfo);

//    treegrav.initialize(nptclmax, 0.5);
    treegrav.initialize(nptclmax, 0.0);

#define DEBUG
#ifdef DEBUG
    treegrav.initializeLocalTree(half_len);
    treegrav.setParticleLocalTree(system);
    treegrav.mortonSortLocalTreeOnly();
    treegrav.linkCellLocalTreeOnly();    
    treegrav.calcMomentLocalTreeOnly();
    treegrav.exchangeLocalEssentialTree(dinfo);
    treegrav.setLocalEssentialTreeToGlobalTree();
    treegrav.mortonSortGlobalTreeOnly();
    treegrav.linkCellGlobalTreeOnly();
    treegrav.calcMomentGlobalTreeOnly();
    treegrav.makeIPGroup();
    treegrav.calcForceAndWriteBack(calcForceEpEp(), calcForceEpSp(), system);    
#else
    // It doesn't work. Because 'half_len_grav_glb' is too small!
    treegrav.calcForceAllAndWriteBack(calcForceEpEp(), calcForceEpSp(), system, dinfo);
#endif

    {
        char filename[1024];
        sprintf(filename, "dump%02d.log", rank);
        FILE *fp = fopen(filename, "w");
        PS::S32 nptcl = system.getNumberOfParticleLocal();
        for(PS::S32 i = 0; i < nptcl; i++) {
            fprintf(fp, "%6d", system[i].id);
            /*
            fprintf(fp, " %+e %+e %+e",
                    system[i].pos.x, system[i].pos.y, system[i].pos.z);
            */
            fprintf(fp, " %+e %+e %+e",
                    system[i].acc.x, system[i].acc.y, system[i].acc.z);
            fprintf(fp, " %+e", system[i].pot);
            fprintf(fp, " %6d", system[i].nj);
            fprintf(fp, " %6d\n", system[i].njreal);
        }
        fclose(fp);
    }

    PS::Finalize();

    return 0;
}
