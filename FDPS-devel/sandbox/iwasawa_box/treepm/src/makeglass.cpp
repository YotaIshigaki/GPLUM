#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>
#include <float.h>
#include <cstdio>
#include <cstdlib>
#include <particle_simulator.hpp>
#include <particle_mesh.hpp>

#include "treepm.hpp"
#include "run_param.hpp"
#include "cosmology.hpp"
#include "mpi.h"

#include "prototype.h"

template <class Tp>
void setup_random_particle(Tp &ptcl, run_param &this_run, 
			   const PS::S32 npart_total)
{
  PS::S32 rank = PS::Comm::getRank();
  PS::S64 npart_local = (rank == 0) ? npart_total : 0;

  ptcl.setNumberOfParticleLocal(npart_local);

  this_run.npart_total = npart_total;
  this_run.mpi_nproc = PS::Comm::getNumberOfProc();
  this_run.mpi_rank = rank;

  PS::MT::init_genrand(0);

  for(PS::S32 i=0;i<npart_local;i++) {
    ptcl[i].mass = 3.0*this_run.cosm.omegam/(8.0*AC::PI*(PS::F32)npart_total);

    ptcl[i].pos[0] = PS::MT::genrand_res53();
    ptcl[i].pos[1] = PS::MT::genrand_res53();
    ptcl[i].pos[2] = PS::MT::genrand_res53();
    ptcl[i].vel[0] = 0.0;
    ptcl[i].vel[1] = 0.0;
    ptcl[i].vel[2] = 0.0;
    ptcl[i].eps = 0.1/pow(npart_total,1.0/3.0);
  }

}

int main(int argc, char **argv)
{
  std::cout << std::setprecision(15);
  std::cerr << std::setprecision(15);

  PS::Initialize(argc, argv);

  PS::PM::ParticleMesh pm;
  PS::ParticleSystem<FPtreepm> ptcl;
  PS::DomainInfo domain_info;

  run_param this_run;

  this_run.npart_total = 2097152;
  //  this_run.npart_total = 262144;
  //  this_run.npart_total = 32768;
  this_run.cosm.omegam = 1.0;
  this_run.cosm.omegav = 0.0;
  this_run.cosm.omegab = 0.0;
  this_run.cosm.omeganu = 0.0;

  this_run.step = 0;
  this_run.tnow = 0.05;
  this_run.update_expansion(this_run.tnow);

  ptcl.initialize();
  domain_info.initialize();

  setup_random_particle(ptcl, this_run, this_run.npart_total);

  this_run.input_params(argv[1]);

  domain_info.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  domain_info.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), 
			       PS::F64vec(1.0, 1.0, 1.0));

  domain_info.decomposeDomainAll(ptcl);

  ptcl.exchangeParticle(domain_info);
  this_run.npart_local = ptcl.getNumberOfParticleLocal();

  PS::TreeForForceLong<Result_treepm, EPItreepm, EPJtreepm>::MonopoleWithCutoff treepm_tree;

  treepm_tree.initialize(3*ptcl.getNumberOfParticleGlobal());

  treepm_tree.calcForceAllAndWriteBack
    (calc_pp_force<EPJtreepm>(),
     calc_pp_force<PS::SPJMonopoleCutoff>(),
     ptcl,
     domain_info);

  pm.calcForceAllAndWriteBack(ptcl, domain_info);
  PS::F64 dtime = calc_dtime(ptcl, this_run);
  PS::F64 dtime_prev, dtime_mid;

  this_run.output_diag(dtime);

  drift_ptcl(ptcl, domain_info, 0.5*dtime);

  domain_info.decomposeDomainAll(ptcl);
  ptcl.exchangeParticle(domain_info);
  this_run.npart_local = ptcl.getNumberOfParticleLocal();

  while(this_run.tnow < this_run.tend && 
	this_run.output_indx < this_run.noutput ) {

    treepm_tree.calcForceAllAndWriteBack
      (calc_pp_force<EPJtreepm>(),
       calc_pp_force<PS::SPJMonopoleCutoff>(),
       ptcl,
       domain_info);
    
    pm.calcForceAllAndWriteBack(ptcl, domain_info);

    this_run.tnow += 0.5*dtime;
    this_run.update_expansion(this_run.tnow);

    reverse_ptcl_acc(ptcl);
    kick_ptcl(ptcl, dtime, this_run);


    this_run.tnow += 0.5*dtime;
    this_run.update_expansion(this_run.tnow);

    dtime_prev = dtime;
    dtime = calc_dtime(ptcl, this_run);
    dtime_mid = 0.5*(dtime_prev + dtime);

    drift_ptcl(ptcl, domain_info, dtime_mid);

    domain_info.decomposeDomainAll(ptcl);
    ptcl.exchangeParticle(domain_info);
    this_run.npart_local = ptcl.getNumberOfParticleLocal();

    output_data_in_run(ptcl, this_run);

    this_run.step++;
    this_run.output_diag(dtime);
  }

  PS::Finalize();

  return(0);
}
