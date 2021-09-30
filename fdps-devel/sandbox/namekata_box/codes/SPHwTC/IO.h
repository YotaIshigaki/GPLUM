#pragma once
/* Standard headers */
#include <string>
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "Particle_Class.h"
#include "SPH_Objects_Class.h"

//================================
//* Class Decl.: IO_Controller
//================================
class IO_Controller {
   public:
      //const std::string output_dir = {"output"}; 
      const std::string output_dir; 
      PS::F64 time,dt;
      PS::F64 stop_time;
      PS::F64 output_time;
      PS::F64 output_interval;
      PS::S32 ndump;
      // Constructors
      IO_Controller();
      // Members
      void initialize();
      void set_config(PS::F64 stop_time,PS::F64 output_interval);
      void write_data(SPH_Objects& SPH_objs);
      void read_data(SPH_Objects& SPH_objs, PS::S32 file_num);
   private:
      void write_misc_data();
      void read_misc_data(PS::S32 file_num);
      void write_SPH_data(SPH_Objects& SPH_objs);
      void read_SPH_data(SPH_Objects& SPH_objs, PS::S32 file_num);
};
extern IO_Controller IO_ctrl;
