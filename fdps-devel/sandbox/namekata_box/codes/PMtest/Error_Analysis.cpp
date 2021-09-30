/* Standard headers */
#include <cassert>
#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
/* FDPS headers */
#include <particle_simulator.hpp>
/* PMMM headers */
#include <vector3.h>
inline dvec3 minimum_image(const dvec3 &inp){
   return dvec3(
         inp.x - round(inp.x),
         inp.y - round(inp.y),
         inp.z - round(inp.z));
}
#include <ewald.h>
/* User-defined headers */
#include "preprocess_keywords.h"
#include "Particle_Class.h"
#include "NaCl.h"



void error_analysis(PS::ParticleSystem<Nbody_FP>& system){

   //* Get MPI parallel info.
   PS::S32 nprocs = PS::Comm::getNumberOfProc();
   PS::S32 myrank = PS::Comm::getRank();

   // /* ==================================================
   //    [Test] Computation of Madelung energy
   // ================================================== */
   // {
   //    //* Notification
   //    if (myrank == 0) std::cout << "Madelung energy test started!!" << std::endl;

   //    //* Define NaCl crystal  
   //    //const int NGrid1D = 2;
   //    //const int NGrid1D = 4;
   //    //const int NGrid1D = 6;
   //    const int NGrid1D = 8;
   //    const int NGrid3D = NGrid1D * NGrid1D * NGrid1D;
   //    NaCl<NGrid1D> NaCl_data;

   //    //* Define ptcl[]
   //    Particle ptcl[NGrid3D];
   //    for (int i=0; i<NGrid3D; i++) {
   //       ptcl[i].mass  = NaCl_data.atoms[i].q;
   //       ptcl[i].pos.x = NaCl_data.atoms[i].x;
   //       ptcl[i].pos.y = NaCl_data.atoms[i].y;
   //       ptcl[i].pos.z = NaCl_data.atoms[i].z;
   //       ptcl[i].phi_direct = 0.0;
   //       ptcl[i].acc_direct = 0.0;
   //    }

   //    //* Compute potential & force with the Ewald method
   //    double msum = 0.0;
   //    for(int i=0; i<NGrid3D; i++) msum += ptcl[i].mass;

   //    const double alpha = 2.4;
   //    eval_k_space<5>(NGrid3D, alpha, ptcl);
   //    eval_r_space<3>(NGrid3D, alpha, msum, ptcl);

   //    //* Output
   //    if (myrank == 0) {
   //       //* output total energy
   //       double Ptot{0.0};
   //       for (int i=0; i<NGrid3D; i++) {
   //          Ptot += ptcl[i].mass * ptcl[i].phi_direct;
   //       }
   //       std::cout << std::setprecision(15) 
   //                 << "Ptot = " << Ptot << std::endl;
   //       ////* output potential
   //       //std::string filename {"potential.dat"};
   //       //std::ofstream output_file;
   //       //output_file.open(filename.c_str(),std::ios::trunc);
   //       //for (int i=0; i<NGrid3D; i++) {
   //       //   output_file << ptcl[i].pos.x << " "
   //       //               << ptcl[i].pos.y << " "
   //       //               << ptcl[i].pos.z << " "
   //       //               << ptcl[i].phi_direct 
   //       //               << std::endl;
   //       //}
   //       //output_file.close();
   //    }

   //    PS::Finalize();
   //    std::exit(0);
   // }


   //* Construct a list of # of particles
   PS::S32 numPtclGlobal = system.getNumberOfParticleGlobal();
   PS::S32 numPtclLocal  = system.getNumberOfParticleLocal();
   std::vector<PS::S32> numPtclList(nprocs),sendbuf(nprocs);
   for (PS::S32 i=0; i<nprocs; i++) sendbuf[i]=numPtclLocal;  
   MPI::COMM_WORLD.Alltoall(&sendbuf[0],1,MPI::INT,
                            &numPtclList[0],1,MPI::INT);

   //* Setup ptcl
   PS::F64 m[numPtclGlobal] ,x[numPtclGlobal] ,y[numPtclGlobal] ,z[numPtclGlobal];
   PS::S32 i_start {0};
   for (PS::S32 irank=0; irank<nprocs; irank++) {
      //** Compute i_start
      if (irank > 0){
         i_start += numPtclList[irank-1];
      }
      //** Setup m,x,y,z if irank=myrank
      if (irank == myrank) {
         for (PS::S32 i=0; i<numPtclLocal; i++) {
            PS::S32 ii = i_start + (i+1) - 1;
            m[ii] = system[i].m;
            x[ii] = system[i].x.x;
            y[ii] = system[i].x.y;
            z[ii] = system[i].x.z;
         }
      }
      //** Broadcast
      MPI::COMM_WORLD.Bcast(&m[i_start],numPtclList[irank],MPI::DOUBLE,irank);
      MPI::COMM_WORLD.Bcast(&x[i_start],numPtclList[irank],MPI::DOUBLE,irank);
      MPI::COMM_WORLD.Bcast(&y[i_start],numPtclList[irank],MPI::DOUBLE,irank);
      MPI::COMM_WORLD.Bcast(&z[i_start],numPtclList[irank],MPI::DOUBLE,irank);
   }
   //** Substitute into ptcl[]
   Particle ptcl[numPtclGlobal];
   for (PS::S32 i=0; i<numPtclGlobal; i++) {
      ptcl[i].mass  = m[i];
      ptcl[i].pos.x = x[i];
      ptcl[i].pos.y = y[i];
      ptcl[i].pos.z = z[i];
      ptcl[i].phi_direct = 0.0;
      ptcl[i].acc_direct = 0.0;
   }
   //* Check ptcl[]
   if (myrank == 0) {
      std::string filename {"ptcl.txt"};
      std::ofstream output_file;
      output_file.open(filename.c_str(),std::ios::trunc);
      for (PS::S32 i=0; i<numPtclGlobal; i++) {
         output_file << std::setprecision(15)  
                     << std::showpos           
                     << ptcl[i].pos.x << " "   
                     << ptcl[i].pos.y << " "   
                     << ptcl[i].pos.z << " "   
                     << ptcl[i].mass  << " "   
                     << std::endl;             
      }
      output_file.close();
   }

   //* Compute the total mass
   double msum = 0.0;
   for(PS::S32 i=0; i<numPtclGlobal; i++) msum += ptcl[i].mass;
   if (myrank == 0) std::cout << "[check2] msum = " << msum << std::endl;

   //* Compute the Ewald summation
   FILE *fp = fopen("ewald.dat","r");
   if (fp) {
      for (PS::S32 i=0; i<numPtclGlobal; i++) {
         fscanf(fp,"%la %la %la %la\n",
                &ptcl[i].phi_direct,
                &ptcl[i].acc_direct.x,
                &ptcl[i].acc_direct.y,
                &ptcl[i].acc_direct.z);
      }
   } else {
      const double alpha = 2.4;
      eval_k_space<5>(numPtclGlobal, alpha, ptcl);
      eval_r_space<3>(numPtclGlobal, alpha, msum, ptcl);

      // We need to perform the sign inversion because ewald.h
      // assumes the potential shape of 1/r, instead of -1/r.
      for (PS::S32 i=0; i<numPtclGlobal; i++) {
         ptcl[i].phi_direct *= -1.0;
      }

      FILE *fp = fopen("ewald.dat","w");
      for (PS::S32 i=0; i<numPtclGlobal; i++) {
         fprintf(fp,"%a %a %a %a\n",
                 ptcl[i].phi_direct,
                 ptcl[i].acc_direct.x,
                 ptcl[i].acc_direct.y,
                 ptcl[i].acc_direct.z);
      }
      fclose(fp);
   }

   //* Compute RMS error
   PS::F64 RMS_p_error{0.0},RMS_f_error{0.0};
   PS::F64 Del_p_avrg{0.0};
   i_start = 0;
   if (myrank > 0)
      for (PS::S32 irank=0; irank<myrank; irank++)
         i_start += numPtclList[irank];
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      PS::S32 ii = i_start + (i+1) - 1;
      //* Potential error
      PS::F64 dp = system[i].pot - ptcl[ii].phi_direct;
      PS::F64 p2 = ptcl[ii].phi_direct * ptcl[ii].phi_direct;
      RMS_p_error += (dp*dp)/p2;
      Del_p_avrg += dp;

      //* Force error
      PS::F64vec df;
      df.x = system[i].agrv.x - ptcl[ii].acc_direct.x;
      df.y = system[i].agrv.y - ptcl[ii].acc_direct.y;
      df.z = system[i].agrv.z - ptcl[ii].acc_direct.z;
      PS::F64 f2 = ptcl[ii].acc_direct * ptcl[ii].acc_direct;
      RMS_f_error += (df*df)/f2;
      
      //* check
      PS::F64vec a = system[i].agrv;
      std::cout << "---------------------------- " << std::endl;
      std::cout << "ii          = " << ii << std::endl;
      std::cout << "m           = " << ptcl[ii].mass << std::endl;
      std::cout << "x           = " << ptcl[ii].pos.x << std::endl;
      std::cout << "y           = " << ptcl[ii].pos.y << std::endl;
      std::cout << "z           = " << ptcl[ii].pos.z << std::endl;
      std::cout << "p (P3M)     = " << system[i].pot << std::endl;
      std::cout << "p (Ewald)   = " << ptcl[ii].phi_direct << std::endl;
      std::cout << "delta p     = " << dp << std::endl;
      std::cout << "|f| (P3M)   = " << std::sqrt(a * a) << std::endl;
      std::cout << "|f| (Ewald) = " << std::sqrt(f2) << std::endl;
      std::cout << "|df|        = " << std::sqrt(df * df) << std::endl;
   }
   RMS_p_error = std::sqrt(PS::Comm::getSum(RMS_p_error)/numPtclGlobal);
   RMS_f_error = std::sqrt(PS::Comm::getSum(RMS_f_error)/numPtclGlobal);
   Del_p_avrg  = PS::Comm::getSum(Del_p_avrg)/numPtclGlobal;
   //* Compute Del_p_disp
   PS::F64 Del_p_disp{0.0};
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      PS::S32 ii = i_start + (i+1) - 1;
      PS::F64 dp = system[i].pot - ptcl[ii].phi_direct;
      PS::F64 disp = dp - Del_p_avrg;
      Del_p_disp += disp*disp;
   }
   Del_p_disp = std::sqrt(PS::Comm::getSum(Del_p_disp)/numPtclGlobal);

   //* Total potential
   PS::F64 Ptot_TreePM{0.0},Ptot_Ewald{0.0};
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      PS::S32 ii = i_start + (i+1) - 1;
      Ptot_TreePM += system[i].m * system[i].pot;
      Ptot_Ewald  += ptcl[ii].mass * ptcl[ii].phi_direct;
   }
   Ptot_TreePM = PS::Comm::getSum(Ptot_TreePM);
   Ptot_Ewald  = PS::Comm::getSum(Ptot_Ewald);

   //* Output 
   MPI::COMM_WORLD.Barrier();
   if (PS::Comm::getRank() == 0) {
      std::cout << "******** Summary report ********" << std::endl;
      std::cout << "[1] Potential error = " << RMS_p_error << std::endl;
      std::cout << "    - Del.P (avrg)  = " << Del_p_avrg  << std::endl; 
      std::cout << "    - Del.P (disp)  = " << Del_p_disp  << std::endl; 
      std::cout << "    - Ptot (TreePM) = " << Ptot_TreePM << std::endl;
      std::cout << "    - Ptot (Ewald)  = " << Ptot_Ewald  << std::endl;
      std::cout << "[2] Force error     = " << RMS_f_error << std::endl;
   }
   //* Output
   std::ostringstream string_stream;
   string_stream << "pot"
                 << std::setfill('0')
                 << std::setw(5)
                 << PS::Comm::getRank()
                 << ".dat";
   std::string filename = string_stream.str();
   std::ofstream output_file;
   output_file.open(filename.c_str(),std::ios::trunc);
   for (PS::S32 i=0; i<numPtclLocal; i++) {
      PS::S32 ii = i_start + (i+1) - 1;
      PS::F64 rij = std::sqrt(system[i].x * system[i].x);
      output_file << rij << " "
                  << system[i].pot << " "
                  << ptcl[ii].phi_direct << " "
                  << std::sqrt(std::sqrt(system[i].agrv * system[i].agrv)) << " "
                  << std::sqrt(std::sqrt(ptcl[ii].acc_direct * ptcl[ii].acc_direct)) << " "
                  << std::endl;
   }
   output_file.close();
   

   PS::Finalize();
   std::exit(0);

};

PS::F64 calc_NaCl_error(PS::ParticleSystem<Nbody_FP>& system,
                        PS::S32 nstep){

   //* Analytical solution for \sum q_{i}*\phi_{i}
   //PS::F64 E0 {1.747558e0}; // Wikipedia
   PS::F64 E0 {1.7475645946332}; // Ewald method
   PS::F64 Ebuf {0.0}, Esum {0.0};

   //* Get MPI parallel info.
   PS::S32 nprocs = PS::Comm::getNumberOfProc();
   PS::S32 myrank = PS::Comm::getRank();

   //* Particle num.
   PS::S32 numPtclLocal  = system.getNumberOfParticleLocal();
   PS::S32 numPtclGlobal = system.getNumberOfParticleGlobal(); 
   std::vector<PS::S32> numPtclList(nprocs),sendbuf(nprocs);
   for (PS::S32 i=0; i<nprocs; i++) sendbuf[i]=numPtclLocal;  
   MPI::COMM_WORLD.Alltoall(&sendbuf[0],1,MPI::INT,
                            &numPtclList[0],1,MPI::INT);

   //* Compute particle energy
   PS::S32 i_start {0};
   std::vector<PS::F64> E(numPtclGlobal);
   for (PS::S32 irank=0; irank<nprocs; irank++) {
      //** Compute i_start
      if (irank > 0){
         i_start += numPtclList[irank-1];
      }
      //** Setup m,x,y,z if irank=myrank
      if (irank == myrank) {
         for (PS::S32 i=0; i<numPtclLocal; i++) {
            PS::S32 ii = i_start + i;
            // rescaled potential energy of individual particles
            E[ii] = numPtclGlobal * system[i].m * system[i].pot;
            // total energy
            Ebuf += system[i].m * system[i].pot;
         }
      }
      //** Broadcast
      MPI::COMM_WORLD.Bcast(&E[i_start],numPtclList[irank],MPI::DOUBLE,irank);
   }
   MPI::COMM_WORLD.Allreduce(&Ebuf,&Esum,1,MPI::DOUBLE,MPI::SUM);
   std::sort(E.begin(),E.end());

   //* Output [1]
   if (myrank == 0) {
      std::cout << std::setprecision(15) 
                << "|Esum - E0| = "
                << std::abs(Esum-E0)
                << std::endl;
   }

   //* Output [2]
   //if (myrank == 0) {
   //   //* E.dat
   //   std::string filename {"E.dat"};
   //   std::ofstream output_file;
   //   output_file.open(filename.c_str(),std::ios::trunc);
   //   for (PS::S32 i=0; i<numPtclGlobal; i++) {
   //      output_file << std::setprecision(15)
   //                  << (i*1.0)/numPtclGlobal << " "
   //                  << E[i] << std::endl;
   //   }
   //   output_file.close();
   //   //* dE.dat
   //   filename = "dE.dat";
   //   output_file.open(filename.c_str(),std::ios::trunc);
   //   for (PS::S32 i=0; i<numPtclGlobal; i++) {
   //      output_file << std::setprecision(15)
   //                  << (i*1.0)/numPtclGlobal << " "
   //                  << E[i]-E0 << std::endl;
   //   }
   //   output_file.close();
   //}

   return std::abs(Esum-E0)/E0;
   
}
