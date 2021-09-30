/* Standard headers */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "preprocess_keywords.h"
#include "Physical_Constants.h"

namespace simulation_units {

//* Global variables
double unit_leng;
double unit_mass;
double unit_velc;
double unit_dens;
double unit_time;
double unit_engy;
double unit_pres;
double unit_temp;
double unit_rate;
double unit_accl;

/*-------------------------------------------------------------------*/
//////////////////        S U B R O U T I N E        //////////////////
////////////////// < S E T U P _ S I M _ U N I T S > //////////////////
/*-------------------------------------------------------------------*/
void setup_sim_units() {
   using namespace physical_constants;
   //* Local parameters
#if (Simulation_Dimension == 1)
   const double nu=1.0e0;
#elif (Simulation_Dimension == 2)
   const double nu=2.0e0;
#elif (Simulation_Dimension == 3)
   const double nu=3.0e0;
#endif

   //* [1] Define simulation units
   unit_leng = pc;
   unit_mass = Msolar;
   unit_velc = std::sqrt(Ggrav*unit_mass/unit_leng);
   unit_dens = unit_mass/std::pow(unit_leng,nu);
   unit_time = std::sqrt(std::pow(unit_leng,3.0e0)/unit_mass/Ggrav);
   unit_engy = unit_mass*unit_velc*unit_velc;
   unit_pres = unit_engy/std::pow(unit_leng,nu);
   unit_temp = (unit_engy/kBoltz)*(M_H/unit_mass);
   unit_rate = unit_pres/unit_time;
   unit_accl = unit_velc/unit_time;
   //* [Notes]
   //   In this units, the equation of state take the form of
   //     p = \rho*Tgas/\mu,
   //   where, \rho and Tgas is normalized quantities and
   //   \mu is mean molecular weight with respect to MHydro.

   //* [2] Output simulation units
   PS::S32 myrank = PS::Comm::getRank();
   if (myrank == 0) {
      std::cout << "**** The simulation units ***************"    << std::endl;
      std::cout << "unit_leng = " << unit_leng << " [cm]"         << std::endl;
      std::cout << "unit_mass = " << unit_mass << " [g]"          << std::endl;
      std::cout << "unit_velc = " << unit_velc << " [cm/s]"       << std::endl;
      std::cout << "unit_dens = " << unit_dens << " [g/cm^3]"     << std::endl;
      std::cout << "unit_time = " << unit_time << " [s]"          << std::endl;
      std::cout << "unit_engy = " << unit_engy << " [erg]"        << std::endl;
      std::cout << "unit_pres = " << unit_pres << " [erg/cm^3]"   << std::endl;
      std::cout << "unit_temp = " << unit_temp << " [K]"          << std::endl;
      std::cout << "unit_rate = " << unit_rate << " [erg/cm^3/s]" << std::endl;
      std::cout << "unit_accl = " << unit_accl << " [cm/s^2]"     << std::endl;
      std::cout << "*****************************************"    << std::endl;
      std::cout << std::endl;

      std::string filename;
      std::ofstream output_file;
      filename = "simulation_units.dat";
      output_file.open(filename.c_str(),std::ios::trunc | std::ios::binary);
         output_file << unit_leng << std::endl;
         output_file << unit_mass << std::endl;
         output_file << unit_velc << std::endl;
         output_file << unit_dens << std::endl;
         output_file << unit_time << std::endl;
         output_file << unit_engy << std::endl;
         output_file << unit_pres << std::endl;
         output_file << unit_temp << std::endl;
         output_file << unit_rate << std::endl;
         output_file << unit_accl << std::endl;
      output_file.close();
   }
}

/*-------------------------------------------------------------------*/
////////////////////       S U B R O U T I N E       //////////////////
//////////////////// < R E A D _ S I M _ U N I T S > //////////////////
/*-------------------------------------------------------------------*/
void read_sim_units() { 
   //* Local variables
   int myrank;
   std::string filename;
   std::ifstream input_file;

   myrank = PS::Comm::getRank();
   if (myrank == 0) {
      filename = "simulation_units.dat";
      input_file.open(filename.c_str(),std::ios::in | std::ios::binary);
         input_file >> unit_leng;
         input_file >> unit_mass;
         input_file >> unit_velc;
         input_file >> unit_dens;
         input_file >> unit_time;
         input_file >> unit_engy;
         input_file >> unit_pres;
         input_file >> unit_temp;
         input_file >> unit_rate;
         input_file >> unit_accl;
      input_file.close();
   }
   PS::Comm::broadcast(&unit_leng,1,0);
   PS::Comm::broadcast(&unit_mass,1,0);
   PS::Comm::broadcast(&unit_velc,1,0);
   PS::Comm::broadcast(&unit_dens,1,0);
   PS::Comm::broadcast(&unit_time,1,0);
   PS::Comm::broadcast(&unit_engy,1,0);
   PS::Comm::broadcast(&unit_pres,1,0);
   PS::Comm::broadcast(&unit_temp,1,0);
   PS::Comm::broadcast(&unit_rate,1,0);
   PS::Comm::broadcast(&unit_accl,1,0);

}

} // END of simulation_units
