#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <limits>

#include "particle_simulator.hpp"
#include "particle.hpp"
#include "force.hpp"

static FILE * FilePointer;

inline void debugFunction(const char * const literal) {
    if(PS::Comm::getRank() == 0) {
        printf("hogehogehoge %s\n", literal);
        fflush(stdout);
    }
}

template <class T>
void outputTimeProfile(PS::S64 nstep,
                       const char * const comment,
                       T & obj) {
    PS::TimeProfile time = obj.getTimeProfile();
    PS::F64         ttot = time.getTotalTime();

    PS::F64 collectMax = PS::Comm::getMaxValue(time.collect_sample_particle);
    PS::F64 decompoMax = PS::Comm::getMaxValue(time.decompose_domain);
    PS::F64 exchanpMax = PS::Comm::getMaxValue(time.exchange_particle);
    PS::F64 makelocMax = PS::Comm::getMaxValue(time.make_local_tree);
    PS::F64 makegloMax = PS::Comm::getMaxValue(time.make_global_tree);
    PS::F64 calcforMax = PS::Comm::getMaxValue(time.calc_force);
    PS::F64 calcmloMax = PS::Comm::getMaxValue(time.calc_moment_local_tree);
    PS::F64 calcmglMax = PS::Comm::getMaxValue(time.calc_moment_global_tree);
    PS::F64 makele1Max = PS::Comm::getMaxValue(time.make_LET_1st);
    PS::F64 makele2Max = PS::Comm::getMaxValue(time.make_LET_2nd);
    PS::F64 exchal1Max = PS::Comm::getMaxValue(time.exchange_LET_1st);
    PS::F64 exchal2Max = PS::Comm::getMaxValue(time.exchange_LET_2nd);
    PS::F64 tttotalMax = PS::Comm::getMaxValue(ttot);

    if(PS::Comm::getRank() == 0) {
        fprintf(FilePointer, "%-35s %e\n", comment, tttotalMax / nstep);
        fprintf(FilePointer, "Collect sample particles: %e\n", collectMax / nstep);
        fprintf(FilePointer, "Decompose domain:         %e\n", decompoMax / nstep);
        fprintf(FilePointer, "Exchange particle:        %e\n", exchanpMax / nstep);
        fprintf(FilePointer, "Make local tree:          %e\n", makelocMax / nstep);
        fprintf(FilePointer, "Make global tree:         %e\n", makegloMax / nstep);
        fprintf(FilePointer, "Calc force:               %e\n", calcforMax / nstep);
        fprintf(FilePointer, "Calc moment local tree:   %e\n", calcmloMax / nstep);
        fprintf(FilePointer, "Calc moment global tree:  %e\n", calcmglMax / nstep);
        fprintf(FilePointer, "Make LET 1st:             %e\n", makele1Max / nstep);
        fprintf(FilePointer, "Make LET 2nd:             %e\n", makele2Max / nstep);
        fprintf(FilePointer, "Exchange LET 1st:         %e\n", exchal1Max / nstep);
        fprintf(FilePointer, "Exchange LET 2nd:         %e\n", exchal2Max / nstep);
        fprintf(FilePointer, "\n");
    }
}

template <class Tptcl>
void generateParticle(PS::S64 nptcl,
                      Tptcl & ptcl) {

    PS::S64 irank = PS::Comm::getRank();
    PS::S64 nproc = PS::Comm::getNumberOfProc();
    PS::S64 nloc  = nptcl * (irank + 1) / nproc - nptcl * irank / nproc;
    ptcl.setNumberOfParticleLocal(nloc);
    PS::MT::init_genrand(irank);
    for(PS::S64 i = 0; i < nloc; i++) {
        ptcl[i].mass   = 1. / (PS::F64)nptcl;
        ptcl[i].size   = 3. / pow(nptcl, 1./3.);
        ptcl[i].pos[0] = PS::MT::genrand_res53();
        ptcl[i].pos[1] = PS::MT::genrand_res53();
        ptcl[i].pos[2] = PS::MT::genrand_res53();
        ptcl[i].vel[0] = PS::MT::genrand_res53();
        ptcl[i].vel[1] = PS::MT::genrand_res53();
        ptcl[i].vel[2] = PS::MT::genrand_res53();
    }

#if 0
    for(PS::S32 i = 0; i < nproc; i++) {
        if(i == irank) {
            fprintf(stdout, "Rank: %8d NoP: %10d\n", irank, nloc);
        }
        PS::Comm::barrier();
    }
#endif
}

template <class Tptcl>
void integrateOrbit(PS::F64 dt,
                    Tptcl & ptcl) {
    for(PS::S64 i = 0; i < ptcl.getNumberOfParticleLocal(); i++) {
        ptcl[i].integrateOrbit(dt);
    }
}

int main(int argc, char **argv)
{
    PS::S64 nptcl = atoi(argv[1]);

    PS::Initialize(argc, argv);

    PS::DomainInfo dinfo;
    dinfo.initialize();

    PS::ParticleSystem<FP> fptcl;
    fptcl.initialize();
    fptcl.createParticle(0);
    fptcl.setNumberOfParticleLocal(0);

    PS::TreeForForceShort<Force, FP, FP>::Gather   gather;
    PS::TreeForForceShort<Force, FP, FP>::Scatter  scatter;
    PS::TreeForForceShort<Force, FP, FP>::Symmetry symmetry;
    PS::TreeForForceLong<Force, FP, FP>::Monopole  longopen;
    PS::TreeForForceLong<Force, FP, FP>::MonopoleWithCutoff longperiodic;
    gather.initialize(0);
    scatter.initialize(0);
    symmetry.initialize(0);
    longopen.initialize(0);
    longperiodic.initialize(0);

    char filename[256];
    sprintf(filename, "bench/weakscaling_p%08lld_t%03lld_n%013lld.log",
            PS::Comm::getNumberOfProc() * omp_get_max_threads(), omp_get_max_threads(), nptcl);
    FilePointer = fopen(filename, "w");

    for(PS::S64 iop = 0; iop < 2; iop++) {

        if(iop != 0) {
            dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
            dinfo.setPosRootDomain(0., 1.);
        }
        
        generateParticle(nptcl, fptcl);
        fptcl.adjustPositionIntoRootDomain(dinfo);
        
        dinfo.decomposeDomainAll(fptcl);
        fptcl.exchangeParticle(dinfo);
        gather.calcForceAllAndWriteBack(calcForce(), fptcl, dinfo);
        scatter.calcForceAllAndWriteBack(calcForce(), fptcl, dinfo);
        symmetry.calcForceAllAndWriteBack(calcForce(), fptcl, dinfo);
        if(iop == 0) {
#error the number of arguments is 4!
            longopen.calcForceAllAndWriteBack(calcForce(), fptcl, dinfo);
        } else {
            longperiodic.calcForceAllAndWriteBack(calcForce(), fptcl, dinfo);
        }
        
        dinfo.clearTimeProfile();
        fptcl.clearTimeProfile();
        gather.clearTimeProfile();
        scatter.clearTimeProfile();
        symmetry.clearTimeProfile();
        if(iop == 0) {
            longopen.clearTimeProfile();
        } else {
            longperiodic.clearTimeProfile();
        }
        
        PS::S64 nrepeat = 10;
        PS::S64 irepeat = 0;
        PS::F64 dt      = 1. / 128.;
        while(irepeat < nrepeat) {
            integrateOrbit(dt, fptcl);
            fptcl.adjustPositionIntoRootDomain(dinfo);

            PS::Comm::barrier();
            dinfo.decomposeDomainAll(fptcl);
            PS::Comm::barrier();
            fptcl.exchangeParticle(dinfo);
            PS::Comm::barrier();
            gather.calcForceAllAndWriteBack(calcForce(), fptcl, dinfo);
            PS::Comm::barrier();
            scatter.calcForceAllAndWriteBack(calcForce(), fptcl, dinfo);
            PS::Comm::barrier();
            symmetry.calcForceAllAndWriteBack(calcForce(), fptcl, dinfo);
            PS::Comm::barrier();
            if(iop == 0) {
                longopen.calcForceAllAndWriteBack(calcForce(), fptcl, dinfo);
            } else {
                longperiodic.calcForceAllAndWriteBack(calcForce(), fptcl, dinfo);
            }
            PS::Comm::barrier();
            
            irepeat++;
        }

        if(PS::Comm::getRank() == 0) {
            fprintf(FilePointer, "#################################################\n");
            if(iop == 0) {
                fprintf(FilePointer, "Open boundary\n");
            } else {
                fprintf(FilePointer, "Periodic boundary\n");
            }
        }
        outputTimeProfile(nrepeat, "Domain Info Total", dinfo);
        outputTimeProfile(nrepeat, "Particle System Total", fptcl);
        outputTimeProfile(nrepeat, "Tree for Force (Gather) Total", gather);
        outputTimeProfile(nrepeat, "Tree for Force (Scatter) Total", scatter);
        outputTimeProfile(nrepeat, "Tree for Force (Symmetry) Total", symmetry);
        if(iop == 0) {
            outputTimeProfile(nrepeat, "Tree for Force (Long) Total", longopen);
        } else {
            outputTimeProfile(nrepeat, "Tree for Force (Long) Total", longperiodic);
        }

    }

    fclose(FilePointer);

    PS::Finalize();

    return 0;
}
