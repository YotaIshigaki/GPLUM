/* Standard headers */
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#endif
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL 
#include <omp.h>
#endif
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "Particle_Class.h"
#include "SPH_Objects_Class.h"
#include "Simulation_Units.h"
#include "IO.h"

//================================
//* Class Defs.: IO_Controller
//================================
IO_Controller::IO_Controller(): output_dir("output") {
   time = 0.0e0;
   dt   = 0.0e0;
   output_time = 0.0e0;
   output_interval = 0.0e0;
   ndump = 0;
}

void IO_Controller::initialize() {
   std::string shell_cmd;
   if (PS::Comm::getRank() == 0) {
      shell_cmd = "mkdir -p " + output_dir;
      system(shell_cmd.c_str());
   }
}

void IO_Controller::set_config(PS::F64 stop_time,
                               PS::F64 output_interval) {
   this->stop_time = stop_time;
   this->output_interval = output_interval;
}

void IO_Controller::write_data(SPH_Objects& SPH_objs) {

   if (time >= output_time) {
      //* Write misc data
      write_misc_data();
      //* Write SPH data
      write_SPH_data(SPH_objs);
      //* Update ndump & output_time (setting for next output)
      ndump++;
      output_time += output_interval; 
      //* Output STDOUT
      if (PS::Comm::getRank() == 0){
         std::cout << "===== OUTPUT ========================" << std::endl;
         std::cout << " output is done at t = " << time  << std::endl;
         std::cout << " Next output will be done at t = " 
                   << output_time << std::endl;
         std::cout << " (stop_time = " << stop_time << ")" << std::endl;
         std::cout << "=====================================" << std::endl;
      }
   }

}

void IO_Controller::read_data(SPH_Objects& SPH_objs, PS::S32 file_num) {
    
   //* Read misc data
   read_misc_data(file_num);
   //* Read SPH data
   read_SPH_data(SPH_objs,file_num);
   //std::ofstream file;
   //file.open("check.txt",std::ios::trunc);
   //   for (PS::S32 i=0; i<SPH_objs.system.getNumberOfParticleLocal(); i++) {
   //      file << SPH_objs.system[i].x.x << " "
   //           << SPH_objs.system[i].x.y << " "
   //           << SPH_objs.system[i].x.z << std::endl;
   //   }
   //file.close();
   //std::exit(0);
   //* Initialize domain info
   SPH_objs.dinfo.initialize(); 
   //SPH_objs.dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
   //SPH_objs.dinfo.setPosRootDomain(PS::F64vec(-0.5e0*Lx,-0.5e0*Ly,-0.5e0*Lz),                                
   //                                PS::F64vec(+0.5e0*Lx,+0.5e0*Ly,+0.5e0*Lz)); 
   SPH_objs.dinfo.decomposeDomain();
   //* Exchange particles
   SPH_objs.system.exchangeParticle(SPH_objs.dinfo);
   //* Make trees and initialize them
   SPH_objs.init_tree();
   //* Compute density/pressure/accelerations
   SPH_objs.update_h();
   SPH_objs.calc_rho();
   SPH_objs.calc_pres();
   // Note that acceleration and the rate of change of specific internal
   // energy are read from the file.

}

void IO_Controller::write_misc_data() {
   //* Local variables
   std::string filename;
   std::ostringstream string_stream;
   std::ofstream output_file;

   if (PS::Comm::getRank() == 0) {
      string_stream << "misc" 
                    << std::setfill('0') 
                    << std::setw(5) 
                    << ndump << ".dat";
      filename = output_dir + "/" + string_stream.str();
      output_file.open(filename.c_str(),std::ios::trunc | std::ios::binary);
      //* Write the contents of IO_Controller
      output_file.write(reinterpret_cast<char*>(&this->time),
                        sizeof(this->time));
      output_file.write(reinterpret_cast<char*>(&this->dt),
                        sizeof(this->dt));
      output_file.write(reinterpret_cast<char*>(&this->stop_time),
                        sizeof(this->stop_time));
      output_file.write(reinterpret_cast<char*>(&this->output_time),
                        sizeof(this->output_time));
      output_file.write(reinterpret_cast<char*>(&this->output_interval),
                        sizeof(this->output_interval));
      output_file.write(reinterpret_cast<char*>(&this->ndump),
                        sizeof(this->ndump));
      //* Write the other parameters such as....
      output_file.close();
      //* Information
      std::cout << filename << " write complete." << std::endl;
   }

}

void IO_Controller::read_misc_data(PS::S32 file_num) {
   //* Local variables
   std::string filename;
   std::ostringstream string_stream;
   std::ifstream input_file;

   if (PS::Comm::getRank() == 0) {
      string_stream << "misc" 
                    << std::setfill('0') 
                    << std::setw(5) 
                    << file_num << ".dat";
      filename = output_dir + "/" + string_stream.str();
      input_file.open(filename.c_str(),std::ios::in | std::ios::binary);
      //* Read the contents of IO_Controller
      input_file.read(reinterpret_cast<char*>(&this->time),
                      sizeof(this->time));
      input_file.read(reinterpret_cast<char*>(&this->dt),
                      sizeof(this->dt));
      input_file.read(reinterpret_cast<char*>(&this->stop_time),
                      sizeof(this->stop_time));
      input_file.read(reinterpret_cast<char*>(&this->output_time),
                      sizeof(this->output_time));
      input_file.read(reinterpret_cast<char*>(&this->output_interval),
                      sizeof(this->output_interval));
      input_file.read(reinterpret_cast<char*>(&this->ndump),
                      sizeof(this->ndump));
      //* Read the other parameters ....
      input_file.close();
   }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
   //* MPI communication
   MPI::Datatype data_type;
   //** Broadcast IO_Controller 
   data_type = MPI::Datatype::Match_size(MPI::TYPECLASS_REAL,
                                         sizeof(this->time));
   MPI::COMM_WORLD.Bcast(&this->time,1,data_type,0);
   data_type = MPI::Datatype::Match_size(MPI::TYPECLASS_REAL,
                                         sizeof(this->dt));
   MPI::COMM_WORLD.Bcast(&this->dt,1,data_type,0);
   data_type = MPI::Datatype::Match_size(MPI::TYPECLASS_REAL,
                                         sizeof(this->stop_time));
   MPI::COMM_WORLD.Bcast(&this->stop_time,1,data_type,0);
   data_type = MPI::Datatype::Match_size(MPI::TYPECLASS_REAL,
                                         sizeof(this->output_time));
   MPI::COMM_WORLD.Bcast(&this->output_time,1,data_type,0);
   data_type = MPI::Datatype::Match_size(MPI::TYPECLASS_REAL,
                                         sizeof(this->output_interval));
   MPI::COMM_WORLD.Bcast(&this->output_interval,1,data_type,0);
   data_type = MPI::Datatype::Match_size(MPI::TYPECLASS_INTEGER,
                                         sizeof(this->output_interval));
   MPI::COMM_WORLD.Bcast(&this->output_interval,1,data_type,0);
   //** Broadcast the others...
#endif

   if (PS::Comm::getRank() == 0)
      std::cout << filename << " read complete" << std::endl;
}

void IO_Controller::write_SPH_data(SPH_Objects& SPH_objs) {
   //* Local variables
   std::string filename;
   std::ostringstream string_stream;
   PS::S32 numPtclGlbl  = SPH_objs.system.getNumberOfParticleGlobal();
   PS::S32 numPtclLocal = SPH_objs.system.getNumberOfParticleLocal();
   std::vector<SPH_IO> IO_buf;
   SPH_IO tmp;

   //* Set filename 
   string_stream << "hydro" 
                 << std::setfill('0') 
                 << std::setw(5) 
                 << ndump << ".dat";
   filename = output_dir + "/" + string_stream.str();

   //* Store into IO_buf
   for(PS::S32 i=0; i<numPtclLocal; i++){
      tmp.copyFromFP(SPH_objs.system[i]);
      IO_buf.push_back(tmp);
   }
#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL
   //* Thread-parallel or serial only
   std::ofstream output_file;
   output_file.open(filename.c_str(),std::ios::trunc | std::ios::binary);
   output_file.write(reinterpret_cast<const char*>(&numPtclLocal),
                     sizeof(numPtclLocal));
   output_file.write(reinterpret_cast<const char*>(&IO_buf[0]),
                     IO_buf.size() * sizeof(SPH_IO));
   output_file.close();
#else
   //* MPI-parallel
   int myrank = MPI::COMM_WORLD.Get_rank();
   int nprocs = MPI::COMM_WORLD.Get_size();
   MPI::Status status;
   MPI::File output_file = MPI::File::Open(MPI::COMM_WORLD, filename.c_str(),
                                           MPI::MODE_CREATE | MPI::MODE_WRONLY,
                                           MPI::INFO_NULL);
   //* Write # of ptcls at RANK 0
   MPI::Offset offset = {0};
   MPI::Datatype data_type = MPI::Datatype::Match_size(MPI::TYPECLASS_INTEGER,
                                                       sizeof(numPtclGlbl));
   if (myrank == 0) {
      output_file.Write_at(offset,&numPtclGlbl,1,data_type); 
   }
   offset += sizeof(numPtclGlbl);
   //* Write ptcls data at all RANK
   // (1) etype (element type)
   MPI::Datatype etype = MPI::BYTE.Create_contiguous(sizeof(SPH_IO));
   etype.Commit();
   // (2) filtype
   MPI::Datatype filetype = etype.Create_contiguous(numPtclLocal);
   filetype.Commit();
   // (3) disp
   std::vector<PS::S32> sendbuf(nprocs), ptcl_nums(nprocs);
   for (PS::S32 i=0; i<nprocs; i++) sendbuf[i]=numPtclLocal;
   MPI::COMM_WORLD.Alltoall(&sendbuf[0],1,MPI::INT,
                            &ptcl_nums[0],1,MPI::INT);
   MPI::Offset disp = {offset};
   if (myrank > 0) {
      for (int irank=0; irank<myrank; irank++) 
         disp += ptcl_nums[irank] * sizeof(SPH_IO); // in bytes
   }
   // (4) write
   output_file.Set_view(disp,etype,filetype,"native",MPI::INFO_NULL);
   output_file.Write_all(&IO_buf[0],numPtclLocal,etype,status);
   output_file.Close();
   // (5) free objects
   etype.Free();
   filetype.Free();
#endif
   //* Information
   if (PS::Comm::getRank() == 0)
      std::cout << filename << " write complete." << std::endl;

}

void IO_Controller::read_SPH_data(SPH_Objects& SPH_objs, PS::S32 file_num) {
   //* Local variables
   std::string filename;
   std::ostringstream string_stream;
   PS::S32 numPtclGlbl;
   PS::S32 numPtclLocal;
   PS::S32 numPtclRem;

   //* Set filename
   string_stream << "hydro" 
                 << std::setfill('0') 
                 << std::setw(5) 
                 << file_num << ".dat";
   filename = output_dir + "/" + string_stream.str();

#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL
   //* Thread-parallel or serial only
   std::ifstream input_file;
   input_file.open(filename.c_str(),std::ios::in | std::ios::binary);
   input_file.read(reinterpret_cast<char*>(&numPtclLocal), 
                   sizeof(numPtclLocal));
   std::vector<SPH_IO> IO_buf(numPtclLocal); // specify # of elements
   input_file.read(reinterpret_cast<char*>(&IO_buf[0]),
                   numPtclLocal * sizeof(SPH_IO)); 
   input_file.close();
#else
   //* MPI-parallel
   int myrank = MPI::COMM_WORLD.Get_rank();
   int nprocs = MPI::COMM_WORLD.Get_size();
   MPI::File input_file = MPI::File::Open(MPI::COMM_WORLD, filename.c_str(),
                                          MPI::MODE_RDONLY, MPI::INFO_NULL);
   //* Read # of ptcls at RANK 0
   MPI::Offset offset = {0};
   MPI::Datatype data_type = MPI::Datatype::Match_size(MPI::TYPECLASS_INTEGER,
                                                       sizeof(numPtclGlbl));
   if (myrank == 0) {
      input_file.Read_at(offset,&numPtclGlbl,1,data_type);
   }
   offset += sizeof(numPtclGlbl); // in bytes
   MPI::COMM_WORLD.Bcast(&numPtclGlbl,1,data_type,0);
   //* Determine numPtclLocal
   numPtclLocal = numPtclGlbl / nprocs;
   if ((numPtclRem = numPtclGlbl % nprocs) != 0) {
      if ((myrank+1) <= numPtclRem) numPtclLocal++;
   }
   //* Read ptcls data
   // (1) etype
   MPI::Datatype etype = MPI::BYTE.Create_contiguous(sizeof(SPH_IO));
   etype.Commit();
   // (2) filetype
   MPI::Datatype filetype = etype.Create_contiguous(numPtclLocal);
   filetype.Commit();
   // (3) disp
   std::vector<PS::S32> sendbuf(nprocs),ptcl_nums(nprocs);
   for (PS::S32 i=0; i<nprocs; i++) sendbuf[i]=numPtclLocal;
   MPI::COMM_WORLD.Alltoall(&sendbuf[0],1,MPI::INT,
                            &ptcl_nums[0],1,MPI::INT);
   MPI::Offset disp = {offset};
   if (myrank > 0) {
      for (int irank=0; irank<myrank; irank++) 
         disp += ptcl_nums[irank] * sizeof(SPH_IO); // in bytes
   }
   // (4) Read
   std::vector<SPH_IO> IO_buf(numPtclLocal);
   MPI::Status status;
   input_file.Set_view(disp,etype,filetype,"native",MPI::INFO_NULL);
   input_file.Read_all(&IO_buf[0],numPtclLocal,etype,status);
   input_file.Close();
   // (5) free objects
   etype.Free();
   filetype.Free();
#endif

   SPH_objs.system.setNumberOfParticleLocal(numPtclLocal);
   for(PS::S32 i=0; i<numPtclLocal; i++){
      SPH_IO tmp = IO_buf[i];
      tmp.copyToFP(SPH_objs.system[i]);
   }

   if (PS::Comm::getRank() == 0)
      std::cout << filename << " read complete" << std::endl;

}

IO_Controller IO_ctrl;

