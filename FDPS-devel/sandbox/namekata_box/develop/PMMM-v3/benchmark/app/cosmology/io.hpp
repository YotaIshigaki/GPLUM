#pragma once
// Include header files of FDPS
#include <particle_simulator.hpp>
// Include user-defined headers
#include "make_directory.h"

template <class Tpsys>
void output_data(Tpsys & psys, run_param & this_run, char *filename) {
  FILE *fp = fopen(filename,"w");
  if (fp == NULL) {
      fprintf(stderr, "File %s cannot be written.",filename);
      exit(EXIT_FAILURE);
  }

  this_run.mpi_rank = PS::Comm::getRank();
  this_run.mpi_nproc = PS::Comm::getNumberOfProc();
  this_run.npart_local = psys.getNumberOfParticleLocal();
  this_run.npart_total = psys.getNumberOfParticleGlobal();
  this_run.write_header(fp);

  for (PS::S64 i=0; i < this_run.npart_local; i++) {
    psys[i].writeParticleBinary(fp);
  }
  fclose(fp);
}

template <class Tpsys>
void output_data_in_run(Tpsys & psys, run_param &this_run) {
    if (this_run.znow < this_run.output_timing[this_run.output_indx]+0.001) {
        if (PS::Comm::getRank() == 0) {
            std::cout << "outputing at z = "
                      << this_run.output_timing[this_run.output_indx]
                      << std::endl;
        }
        static char filename[256], directory_name[256];
        sprintf(directory_name,"%s_%d",
                this_run.model_name, this_run.output_indx);
        sprintf(filename, "%s/%s_%d-%d", 
                directory_name, this_run.model_name, 
                this_run.output_indx, this_run.mpi_rank);
        make_directory(directory_name);
        output_data(psys, this_run, filename);
        this_run.output_indx++;
    }
}
