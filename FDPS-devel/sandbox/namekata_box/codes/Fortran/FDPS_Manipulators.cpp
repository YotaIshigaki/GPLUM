#pragma once
/* Standard headers */
#include <string.h>
#include <string>
#include <cstring>
#include <vector>
#include <typeinfo>
/* FDPS headers */
#include <particle_simulator.hpp>
/* User-defined headers */
#include "FDPS_Manipulators.h"

namespace FDPS_Manipulators {
   //===================
   // Private classes 
   //===================
   class Psys_Data {
   public:
      int id;
      void * ptr;
      std::string info;
   };

   class Dinfo_Data {
   public:
      int id;
      void * ptr;
   };

   class Tree_Data {
   public:
      int id;
      void * ptr;
      std::string info;
      char * info_char;
      std::string type;
      std::string Force;
      std::string EPI;
      std::string EPJ;
      std::string mode;
   };

   class NBS_Target_Data {
   public:
      PS::F64vec pos;
      PS::F64 r_search;

      PS::F64vec getPos() const {
         return this->pos;
      }

      PS::F64 getRSearch() const {
         return this->r_search;
      }
   };

   //===================
   // Private variables 
   //===================
   // -(command-line arguments)
   static int argc_save;
   static char **argv_save;
   // -(ParticleSystem type)
   static int num_psys;
   static int num_psys_creation;
   static std::vector<Psys_Data> psys_vector;
   // -(DomainInfo type)
   static int num_dinfo;
   static int num_dinfo_creation;
   static std::vector<Dinfo_Data> dinfo_vector;
   // -(TreeForForce type)
   static int num_tree;
   static int num_tree_creation;
   static std::vector<Tree_Data> tree_vector;
   // -(MPI comm.)
   // -(Utility)

   //=======================
   // Functions
   //=======================
   //-(Initializer)
   void Initialize(int argc, char *argv[]) {
      argc_save = argc;
      argv_save = argv;
      num_psys  = 0;
      num_dinfo = 0;
      num_tree  = 0;
      num_psys_creation = 0;
      num_dinfo_creation = 0;
      num_tree_creation = 0;
   }
   //-(FDPS initilizer/finalizer)
   void PS_Initialize() {
      PS::Initialize(argc_save,argv_save);
      // Display # of MPI processes and threads
      PS::S32 nprocs = PS::Comm::getNumberOfProc();
      PS::S32 nthrds = PS::Comm::getNumberOfThread();
      PS::S32 myrank = PS::Comm::getRank();
      if (myrank == 0) {
      std::cout << "============================================================" << std::endl
                << " This is a mockup program of Fortran interface of FDPS!"      << std::endl
                << "   # of processes is " << nprocs                              << std::endl
                << "   # of thread is    " << nthrds                              << std::endl
                << "============================================================" << std::endl;
      }
   }
   void PS_Finalize() {
      PS::Finalize();
   }
   void PS_Abort(const int err_num) {
      PS::Abort(err_num);
   }
   //--------------------------------
   // ParticleSystem manipulators
   //--------------------------------
   void create_psys(int *psys_num,
                    char *psys_info) {
      std::string psys_info_ = psys_info;
      Psys_Data psys_data;
      //-----------------------------------------------
      // fdps-autogen:create_psys;
      if (psys_info_ == "full_particle") {
         PS::ParticleSystem<full_particle> *psys;
         psys = new PS::ParticleSystem<full_particle>;
         psys_data.id   = num_psys_creation;
         psys_data.ptr  = (void *) psys;
         psys_data.info = psys_info_;
      
      } else { 
         PS::S32 myrank = PS::Comm::getRank();
         if (myrank == 0) {
             std::cout << "FullParticle `" << psys_info_ << "` "
                       << "does not exist." << std::endl;
         }
         PS::Finalize();
         std::exit(EXIT_FAILURE);
      }
      
      //-----------------------------------------------
      psys_vector.push_back(psys_data);
      *psys_num = psys_data.id;
      num_psys_creation++;
      num_psys++;
      ////std::cout << "psys successfully created!" << std::endl;
   }
   void delete_psys(const int psys_num) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            it = psys_vector.erase(it);
            num_psys--;
         }
      }
   }
   void init_psys(const int psys_num) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:init_psys;
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               psys->initialize();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //-----------------------------------------------
         }
      }
   }
   void get_psys_info(const int psys_num,
                      char *psys_info,
                      size_t *charlen) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            std::strcpy(psys_info,it->info.c_str());
            *charlen  = std::strlen(psys_info);
         }
      }
   }
   long long int get_psys_memsize(const int psys_num) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:get_psys_memsize;
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               return (long long int) psys->getMemSizeUsed();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //-----------------------------------------------
         }
      }
   }
   void get_psys_time_prof(const int psys_num,
                           PS::TimeProfile *prof) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:get_psys_time_prof;
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               *prof = psys->getTimeProfile();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //-----------------------------------------------
         }
      }
   }
   void clear_psys_time_prof(const int psys_num) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:clear_psys_time_prof;
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               psys->clearTimeProfile();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //-----------------------------------------------
         }
      }
   }
   void set_nptcl_smpl(const int psys_num,
                       const int numPtcl){
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:set_nptcl_smpl;
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               psys->setAverageTargetNumberOfSampleParticlePerProcess(numPtcl);
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //-----------------------------------------------
         }
      }
   }
   void set_nptcl_loc(const int psys_num,
                      const int numPtcl){
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:set_nptcl_loc;
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               psys->setNumberOfParticleLocal(numPtcl);
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //-----------------------------------------------
         }
      }
   }
   int get_nptcl_loc(const int psys_num) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:get_nptcl_loc;
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               return psys->getNumberOfParticleLocal();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //-----------------------------------------------
         }
      }
   }
   int get_nptcl_glb(const int psys_num) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:get_nptcl_glb;
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               return psys->getNumberOfParticleGlobal();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //-----------------------------------------------
         }
      }
   }
   void get_psys_cptr(const int psys_num,
                      void **cptr) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:get_psys_cptr;
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               *cptr = (void *) &((*psys)[0]);// [!!IMPORTANT!!]//   The operator [0] is necessary to obtain the correct address of//   the particle system object.
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //-----------------------------------------------
         }
      }
   }
   void exchange_particle(const int psys_num,
                          const int dinfo_num) {
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Get psys and issue exchangeParticle()
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:exchange_particle;
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               psys->exchangeParticle(*dinfo);
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //-----------------------------------------------
         }
      }
   }
   //-----------------------------------------------------------------------------------
   // fdps-autogen:add_particle;
   void add_particle00000(const int psys_num,
                          const full_particle *ptcl) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            PS::ParticleSystem<full_particle> *psys;
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
            psys->addOneParticle(*ptcl);
         }
      }
   }
   
   //-----------------------------------------------------------------------------------
   void remove_particle(const int psys_num, 
                        const int numPtcl,
                        int *ptcl_indx) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:remove_particle;
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               psys->removeParticle(ptcl_indx,numPtcl);
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //-----------------------------------------------
         }
      }
   }
   //----------------------------
   // DomainInfo manipulators
   //----------------------------
   void create_dinfo(int *dinfo_num) {
      Dinfo_Data dinfo_data;
      PS::DomainInfo *dinfo;
      dinfo = new PS::DomainInfo;
      dinfo_data.id = num_dinfo_creation;
      dinfo_data.ptr = (void *) dinfo;
      dinfo_vector.push_back(dinfo_data);
      *dinfo_num = dinfo_data.id;
      num_dinfo_creation++;
      num_dinfo++;
      //std::cout << "domain info successfully created!" << std::endl;
   }
   void delete_dinfo(const int dinfo_num) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            it = dinfo_vector.erase(it);
            num_dinfo--;
         }
      }
   }
   void init_dinfo(const int dinfo_num,
                   const float coef_ema) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            PS::DomainInfo *dinfo;
            dinfo = (PS::DomainInfo *) it->ptr;
            dinfo->initialize(coef_ema);
         }
      }
   }
   void get_dinfo_time_prof(const int dinfo_num,
                            PS::TimeProfile *prof) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            PS::DomainInfo *dinfo;
            dinfo = (PS::DomainInfo *) it->ptr;
            *prof = dinfo->getTimeProfile();
         }
      }
   }
   void clear_dinfo_time_prof(const int dinfo_num) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            PS::DomainInfo *dinfo;
            dinfo = (PS::DomainInfo *) it->ptr;
            dinfo->clearTimeProfile();
         }
      }
   }
   void set_nums_domain(const int dinfo_num,
                        const int nx,
                        const int ny,
                        const int nz) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            PS::DomainInfo *dinfo;
            dinfo = (PS::DomainInfo *) it->ptr;
            dinfo->setDomain(nx,ny,nz);
         }
      }
   }
   void set_boundary_condition(const int dinfo_num,
                               const enum PS_BOUNDARY_CONDITION bc) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            PS::DomainInfo *dinfo;
            dinfo = (PS::DomainInfo *) it->ptr;
            if (bc == BC_OPEN) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
            } else if (bc == BC_PERIODIC_X) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_X);
            } else if (bc == BC_PERIODIC_Y) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_Y);
            } else if (bc == BC_PERIODIC_Z) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_Z);
            } else if (bc == BC_PERIODIC_XY) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XY);
            } else if (bc == BC_PERIODIC_XZ) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XZ);
            } else if (bc == BC_PERIODIC_YZ) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_YZ);
            } else if (bc == BC_PERIODIC_XYZ) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
            } else if (bc == BC_SHEARING_BOX) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_SHEARING_BOX);
            } else if (bc == BC_USER_DEFINED) {
               dinfo->setBoundaryCondition(PS::BOUNDARY_CONDITION_USER_DEFINED);
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) 
                  std::cout << "Unknown boundary condition is specified." << std::endl;
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
         }
      }
   }
   void set_pos_root_domain(const int dinfo_num,
                            const PS::F32vec low,
                            const PS::F32vec high) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            PS::DomainInfo *dinfo;
            dinfo = (PS::DomainInfo *) it->ptr;
            dinfo->setPosRootDomain(low,high);
         }
      }
   }
   void collect_sample_particle(const int dinfo_num,
                                const int psys_num,
                                const bool clear,
                                const float weight) {
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Get psys and issue collectSampleParticle()
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:collect_sample_particle;
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               dinfo->collectSampleParticle(*psys,clear,weight);
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //-----------------------------------------------
         }
      }
   }
   void decompose_domain(const int dinfo_num) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            PS::DomainInfo *dinfo;
            dinfo = (PS::DomainInfo *) it->ptr;
            dinfo->decomposeDomain();
         }
      }
   }
   void decompose_domain_all(const int dinfo_num,
                             const int psys_num,
                             const float weight) {
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Get psys and issue decomposeDomainAll()
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            //-----------------------------------------------
            // fdps-autogen:decompose_domain_all;
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               dinfo->decomposeDomainAll(*psys,weight);
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //-----------------------------------------------
         }
      }
   }
   //----------------------------
   // TreeForForce manipulators
   //----------------------------
   void create_tree(int *tree_num,
                    char *tree_info) {
      std::string tree_info_ = tree_info;
      Tree_Data tree_data;
      //------------------------------------------------------
      // fdps-autogen:create_tree;
      if (tree_info_ == "Long,full_particle,full_particle,full_particle,Monopole") {
         PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
         tree = new PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,full_particle,full_particle,Quadrupole") {
         PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
         tree = new PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter") {
         PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,full_particle,full_particle,DipoleGeometricCenter") {
         PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter") {
         PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,full_particle,essential_particle_j,Monopole") {
         PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
         tree = new PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,full_particle,essential_particle_j,Quadrupole") {
         PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
         tree = new PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter") {
         PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter") {
         PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
         PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,essential_particle_i,full_particle,Monopole") {
         PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
         tree = new PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,essential_particle_i,full_particle,Quadrupole") {
         PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
         tree = new PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter") {
         PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter") {
         PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
         PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,essential_particle_i,essential_particle_j,Monopole") {
         PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
         tree = new PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole") {
         PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
         tree = new PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
         PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
         PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
         PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,full_particle,full_particle,Monopole") {
         PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *tree;
         tree = new PS::TreeForForceLong<force,full_particle,full_particle>::Monopole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,full_particle,full_particle,Quadrupole") {
         PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *tree;
         tree = new PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,full_particle,full_particle,MonopoleGeometricCenter") {
         PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,full_particle,full_particle,DipoleGeometricCenter") {
         PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,full_particle,full_particle,QuadrupoleGeometricCenter") {
         PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,full_particle,essential_particle_j,Monopole") {
         PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *tree;
         tree = new PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,full_particle,essential_particle_j,Quadrupole") {
         PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *tree;
         tree = new PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,full_particle,essential_particle_j,MonopoleGeometricCenter") {
         PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,full_particle,essential_particle_j,DipoleGeometricCenter") {
         PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
         PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,essential_particle_i,full_particle,Monopole") {
         PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *tree;
         tree = new PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,essential_particle_i,full_particle,Quadrupole") {
         PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *tree;
         tree = new PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,essential_particle_i,full_particle,MonopoleGeometricCenter") {
         PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,essential_particle_i,full_particle,DipoleGeometricCenter") {
         PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
         PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,essential_particle_i,essential_particle_j,Monopole") {
         PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *tree;
         tree = new PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,essential_particle_i,essential_particle_j,Quadrupole") {
         PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *tree;
         tree = new PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
         PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
         PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else if (tree_info_ == "Long,force,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
         PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
         tree = new PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter;
         tree_data.id   = num_tree_creation;
         tree_data.ptr  = (void *) tree;
         tree_data.info = tree_info_;
      
      } else {
         PS::S32 myrank = PS::Comm::getRank();
         if (myrank == 0) {
            std::cout << "cannot create TreeForForce `" << tree_info_ << "` "
                      << std::endl;
         }
         PS::Finalize();
         std::exit(EXIT_FAILURE);
      }
      
      //------------------------------------------------------
      tree_vector.push_back(tree_data);
      *tree_num = tree_data.id;
      num_tree_creation++;
      num_tree++;
      //std::cout << "tree sucessfully created!" << std::endl;
   }
   void delete_tree(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            it = tree_vector.erase(it);
            num_tree--;
         }
      }
   }
   void init_tree(const int tree_num,
                  const int numPtcl,
                  const float theta, 
                  const int n_leaf_limit,
                  const int n_group_limit) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:init_tree;
            if (it->info == "Long,full_particle,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid TreeForForce number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //------------------------------------------------------
         }
      }
   }
   void get_tree_info(const int tree_num,
                      char *tree_info,
                      size_t *charlen) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            std::strcpy(tree_info,it->info.c_str());
            *charlen  = std::strlen(tree_info);
         }
      }
   }
   long long int get_tree_memsize(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_tree_memsize;
            if (it->info == "Long,full_particle,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getMemSizeUsed();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid TreeForForce number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //------------------------------------------------------
         }
      }
   }
   void get_tree_time_prof(const int tree_num,
                           PS::TimeProfile *prof) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_tree_time_prof;
            if (it->info == "Long,full_particle,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               *prof = tree->getTimeProfile();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid TreeForForce number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //------------------------------------------------------
         }
      }
   }
   void clear_tree_time_prof(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:clear_tree_time_prof;
            if (it->info == "Long,full_particle,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearTimeProfile();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid TreeForForce number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //------------------------------------------------------
         }
      }
   }
   long long int get_num_interact_ep_ep_loc(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_num_interact_ep_ep_loc;
            if (it->info == "Long,full_particle,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPLocal();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid TreeForForce number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //------------------------------------------------------
         }
      }
   }
   long long int get_num_interact_ep_sp_loc(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_num_interact_ep_sp_loc;
            if (it->info == "Long,full_particle,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPLocal();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid TreeForForce number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //------------------------------------------------------
         }
      }
   }
   long long int get_num_interact_ep_ep_glb(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_num_interact_ep_ep_glb;
            if (it->info == "Long,full_particle,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPEPGlobal();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid TreeForForce number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //------------------------------------------------------
         }
      }
   }
   long long int get_num_interact_ep_sp_glb(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_num_interact_ep_sp_glb;
            if (it->info == "Long,full_particle,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfInteractionEPSPGlobal();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid TreeForForce number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //------------------------------------------------------
         }
      }
   }
   void clear_num_interact(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:clear_num_interact;
            if (it->info == "Long,full_particle,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               tree->clearNumberOfInteraction();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid TreeForForce number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //------------------------------------------------------
         }
      }
   }
   long long int get_num_tree_walk_loc(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_num_tree_walk_loc;
            if (it->info == "Long,full_particle,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkLocal();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid TreeForForce number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //------------------------------------------------------
         }
      }
   }
   long long int get_num_tree_walk_glb(const int tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            //------------------------------------------------------
            // fdps-autogen:get_num_tree_walk_glb;
            if (it->info == "Long,full_particle,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,full_particle,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,full_particle,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,full_particle,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,full_particle,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Monopole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,Quadrupole") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,MonopoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,DipoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else if (it->info == "Long,force,essential_particle_i,essential_particle_j,QuadrupoleGeometricCenter") {
               PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
               tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
               return (long long int) tree->getNumberOfWalkGlobal();
            
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid TreeForForce number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
            
            //------------------------------------------------------
         }
      }
   }
   //-----------------------------------------------------------------------------------
   // fdps-autogen:calc_force_all_and_write_back;
   void calc_force_all_and_write_back_l00000(const int tree_num,
                                             void (*pfunc_ep_ep)(struct full_particle *,
                                                                 int ,
                                                                 struct full_particle *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct full_particle *,
                                                                 int ,
                                                                 PS::SPJMonopole *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00001(const int tree_num,
                                             void (*pfunc_ep_ep)(struct full_particle *,
                                                                 int ,
                                                                 struct full_particle *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct full_particle *,
                                                                 int ,
                                                                 PS::SPJQuadrupole *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00002(const int tree_num,
                                             void (*pfunc_ep_ep)(struct full_particle *,
                                                                 int ,
                                                                 struct full_particle *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct full_particle *,
                                                                 int ,
                                                                 PS::SPJMonopoleGeometricCenter *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00003(const int tree_num,
                                             void (*pfunc_ep_ep)(struct full_particle *,
                                                                 int ,
                                                                 struct full_particle *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct full_particle *,
                                                                 int ,
                                                                 PS::SPJDipoleGeometricCenter *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00004(const int tree_num,
                                             void (*pfunc_ep_ep)(struct full_particle *,
                                                                 int ,
                                                                 struct full_particle *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct full_particle *,
                                                                 int ,
                                                                 PS::SPJQuadrupoleGeometricCenter *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00005(const int tree_num,
                                             void (*pfunc_ep_ep)(struct full_particle *,
                                                                 int ,
                                                                 struct essential_particle_j *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct full_particle *,
                                                                 int ,
                                                                 PS::SPJMonopole *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00006(const int tree_num,
                                             void (*pfunc_ep_ep)(struct full_particle *,
                                                                 int ,
                                                                 struct essential_particle_j *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct full_particle *,
                                                                 int ,
                                                                 PS::SPJQuadrupole *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00007(const int tree_num,
                                             void (*pfunc_ep_ep)(struct full_particle *,
                                                                 int ,
                                                                 struct essential_particle_j *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct full_particle *,
                                                                 int ,
                                                                 PS::SPJMonopoleGeometricCenter *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00008(const int tree_num,
                                             void (*pfunc_ep_ep)(struct full_particle *,
                                                                 int ,
                                                                 struct essential_particle_j *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct full_particle *,
                                                                 int ,
                                                                 PS::SPJDipoleGeometricCenter *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00009(const int tree_num,
                                             void (*pfunc_ep_ep)(struct full_particle *,
                                                                 int ,
                                                                 struct essential_particle_j *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct full_particle *,
                                                                 int ,
                                                                 PS::SPJQuadrupoleGeometricCenter *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00010(const int tree_num,
                                             void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                 int ,
                                                                 struct full_particle *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                 int ,
                                                                 PS::SPJMonopole *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00011(const int tree_num,
                                             void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                 int ,
                                                                 struct full_particle *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                 int ,
                                                                 PS::SPJQuadrupole *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00012(const int tree_num,
                                             void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                 int ,
                                                                 struct full_particle *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                 int ,
                                                                 PS::SPJMonopoleGeometricCenter *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00013(const int tree_num,
                                             void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                 int ,
                                                                 struct full_particle *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                 int ,
                                                                 PS::SPJDipoleGeometricCenter *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00014(const int tree_num,
                                             void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                 int ,
                                                                 struct full_particle *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                 int ,
                                                                 PS::SPJQuadrupoleGeometricCenter *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00015(const int tree_num,
                                             void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                 int ,
                                                                 struct essential_particle_j *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                 int ,
                                                                 PS::SPJMonopole *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00016(const int tree_num,
                                             void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                 int ,
                                                                 struct essential_particle_j *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                 int ,
                                                                 PS::SPJQuadrupole *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00017(const int tree_num,
                                             void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                 int ,
                                                                 struct essential_particle_j *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                 int ,
                                                                 PS::SPJMonopoleGeometricCenter *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00018(const int tree_num,
                                             void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                 int ,
                                                                 struct essential_particle_j *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                 int ,
                                                                 PS::SPJDipoleGeometricCenter *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_and_write_back_l00019(const int tree_num,
                                             void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                                 int ,
                                                                 struct essential_particle_j *,
                                                                 int ,
                                                                 struct full_particle *),
                                             void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                                 int ,
                                                                 PS::SPJQuadrupoleGeometricCenter *,
                                                                 int ,
                                                                 struct full_particle *),
                                             const int psys_num,
                                             const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   //-----------------------------------------------------------------------------------
   // fdps-autogen:calc_force_all;
   void calc_force_all_l00000(const int tree_num,
                              void (*pfunc_ep_ep)(struct full_particle *,
                                                  int ,
                                                  struct full_particle *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct full_particle *,
                                                  int ,
                                                  PS::SPJMonopole *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00001(const int tree_num,
                              void (*pfunc_ep_ep)(struct full_particle *,
                                                  int ,
                                                  struct full_particle *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct full_particle *,
                                                  int ,
                                                  PS::SPJQuadrupole *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00002(const int tree_num,
                              void (*pfunc_ep_ep)(struct full_particle *,
                                                  int ,
                                                  struct full_particle *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct full_particle *,
                                                  int ,
                                                  PS::SPJMonopoleGeometricCenter *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00003(const int tree_num,
                              void (*pfunc_ep_ep)(struct full_particle *,
                                                  int ,
                                                  struct full_particle *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct full_particle *,
                                                  int ,
                                                  PS::SPJDipoleGeometricCenter *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00004(const int tree_num,
                              void (*pfunc_ep_ep)(struct full_particle *,
                                                  int ,
                                                  struct full_particle *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct full_particle *,
                                                  int ,
                                                  PS::SPJQuadrupoleGeometricCenter *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00005(const int tree_num,
                              void (*pfunc_ep_ep)(struct full_particle *,
                                                  int ,
                                                  struct essential_particle_j *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct full_particle *,
                                                  int ,
                                                  PS::SPJMonopole *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00006(const int tree_num,
                              void (*pfunc_ep_ep)(struct full_particle *,
                                                  int ,
                                                  struct essential_particle_j *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct full_particle *,
                                                  int ,
                                                  PS::SPJQuadrupole *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00007(const int tree_num,
                              void (*pfunc_ep_ep)(struct full_particle *,
                                                  int ,
                                                  struct essential_particle_j *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct full_particle *,
                                                  int ,
                                                  PS::SPJMonopoleGeometricCenter *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00008(const int tree_num,
                              void (*pfunc_ep_ep)(struct full_particle *,
                                                  int ,
                                                  struct essential_particle_j *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct full_particle *,
                                                  int ,
                                                  PS::SPJDipoleGeometricCenter *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00009(const int tree_num,
                              void (*pfunc_ep_ep)(struct full_particle *,
                                                  int ,
                                                  struct essential_particle_j *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct full_particle *,
                                                  int ,
                                                  PS::SPJQuadrupoleGeometricCenter *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00010(const int tree_num,
                              void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                  int ,
                                                  struct full_particle *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                  int ,
                                                  PS::SPJMonopole *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00011(const int tree_num,
                              void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                  int ,
                                                  struct full_particle *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                  int ,
                                                  PS::SPJQuadrupole *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00012(const int tree_num,
                              void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                  int ,
                                                  struct full_particle *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                  int ,
                                                  PS::SPJMonopoleGeometricCenter *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00013(const int tree_num,
                              void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                  int ,
                                                  struct full_particle *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                  int ,
                                                  PS::SPJDipoleGeometricCenter *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00014(const int tree_num,
                              void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                  int ,
                                                  struct full_particle *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                  int ,
                                                  PS::SPJQuadrupoleGeometricCenter *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00015(const int tree_num,
                              void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                  int ,
                                                  struct essential_particle_j *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                  int ,
                                                  PS::SPJMonopole *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00016(const int tree_num,
                              void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                  int ,
                                                  struct essential_particle_j *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                  int ,
                                                  PS::SPJQuadrupole *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00017(const int tree_num,
                              void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                  int ,
                                                  struct essential_particle_j *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                  int ,
                                                  PS::SPJMonopoleGeometricCenter *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00018(const int tree_num,
                              void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                  int ,
                                                  struct essential_particle_j *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                  int ,
                                                  PS::SPJDipoleGeometricCenter *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   void calc_force_all_l00019(const int tree_num,
                              void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                  int ,
                                                  struct essential_particle_j *,
                                                  int ,
                                                  struct full_particle *),
                              void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                  int ,
                                                  PS::SPJQuadrupoleGeometricCenter *,
                                                  int ,
                                                  struct full_particle *),
                              const int psys_num,
                              const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAll(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo);
   }
   
   //-----------------------------------------------------------------------------------
   // fdps-autogen:calc_force_making_tree;
   void calc_force_making_tree_l00000(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJMonopole *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00001(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJQuadrupole *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00002(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJMonopoleGeometricCenter *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00003(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJDipoleGeometricCenter *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00004(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJQuadrupoleGeometricCenter *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00005(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJMonopole *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00006(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJQuadrupole *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00007(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJMonopoleGeometricCenter *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00008(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJDipoleGeometricCenter *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00009(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJQuadrupoleGeometricCenter *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00010(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJMonopole *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00011(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJQuadrupole *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00012(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJMonopoleGeometricCenter *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00013(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJDipoleGeometricCenter *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00014(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJQuadrupoleGeometricCenter *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00015(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJMonopole *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00016(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJQuadrupole *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00017(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJMonopoleGeometricCenter *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00018(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJDipoleGeometricCenter *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00019(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct full_particle *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJQuadrupoleGeometricCenter *,
                                                          int ,
                                                          struct full_particle *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00020(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJMonopole *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Monopole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00021(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJQuadrupole *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,full_particle,full_particle>::Quadrupole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00022(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJMonopoleGeometricCenter *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00023(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJDipoleGeometricCenter *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00024(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJQuadrupoleGeometricCenter *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00025(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJMonopole *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Monopole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00026(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJQuadrupole *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00027(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJMonopoleGeometricCenter *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00028(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJDipoleGeometricCenter *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00029(const int tree_num,
                                      void (*pfunc_ep_ep)(struct full_particle *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct full_particle *,
                                                          int ,
                                                          PS::SPJQuadrupoleGeometricCenter *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00030(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJMonopole *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Monopole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00031(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJQuadrupole *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00032(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJMonopoleGeometricCenter *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00033(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJDipoleGeometricCenter *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00034(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct full_particle *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJQuadrupoleGeometricCenter *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00035(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJMonopole *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00036(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJQuadrupole *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00037(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJMonopoleGeometricCenter *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00038(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJDipoleGeometricCenter *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   void calc_force_making_tree_l00039(const int tree_num,
                                      void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                          int ,
                                                          struct essential_particle_j *,
                                                          int ,
                                                          struct force *),
                                      void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                          int ,
                                                          PS::SPJQuadrupoleGeometricCenter *,
                                                          int ,
                                                          struct force *),
                                      const int dinfo_num) {
      // Get tree
      PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<force,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceMakingTree(pfunc_ep_ep, pfunc_ep_sp, *dinfo);
   }
   
   //-----------------------------------------------------------------------------------
   // fdps-autogen:calc_force_and_write_back;
   void calc_force_and_write_back_l00000(const int tree_num,
                                         void (*pfunc_ep_ep)(struct full_particle *,
                                                             int ,
                                                             struct full_particle *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct full_particle *,
                                                             int ,
                                                             PS::SPJMonopole *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00001(const int tree_num,
                                         void (*pfunc_ep_ep)(struct full_particle *,
                                                             int ,
                                                             struct full_particle *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct full_particle *,
                                                             int ,
                                                             PS::SPJQuadrupole *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Quadrupole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00002(const int tree_num,
                                         void (*pfunc_ep_ep)(struct full_particle *,
                                                             int ,
                                                             struct full_particle *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct full_particle *,
                                                             int ,
                                                             PS::SPJMonopoleGeometricCenter *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00003(const int tree_num,
                                         void (*pfunc_ep_ep)(struct full_particle *,
                                                             int ,
                                                             struct full_particle *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct full_particle *,
                                                             int ,
                                                             PS::SPJDipoleGeometricCenter *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00004(const int tree_num,
                                         void (*pfunc_ep_ep)(struct full_particle *,
                                                             int ,
                                                             struct full_particle *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct full_particle *,
                                                             int ,
                                                             PS::SPJQuadrupoleGeometricCenter *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00005(const int tree_num,
                                         void (*pfunc_ep_ep)(struct full_particle *,
                                                             int ,
                                                             struct essential_particle_j *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct full_particle *,
                                                             int ,
                                                             PS::SPJMonopole *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Monopole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00006(const int tree_num,
                                         void (*pfunc_ep_ep)(struct full_particle *,
                                                             int ,
                                                             struct essential_particle_j *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct full_particle *,
                                                             int ,
                                                             PS::SPJQuadrupole *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::Quadrupole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00007(const int tree_num,
                                         void (*pfunc_ep_ep)(struct full_particle *,
                                                             int ,
                                                             struct essential_particle_j *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct full_particle *,
                                                             int ,
                                                             PS::SPJMonopoleGeometricCenter *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00008(const int tree_num,
                                         void (*pfunc_ep_ep)(struct full_particle *,
                                                             int ,
                                                             struct essential_particle_j *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct full_particle *,
                                                             int ,
                                                             PS::SPJDipoleGeometricCenter *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00009(const int tree_num,
                                         void (*pfunc_ep_ep)(struct full_particle *,
                                                             int ,
                                                             struct essential_particle_j *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct full_particle *,
                                                             int ,
                                                             PS::SPJQuadrupoleGeometricCenter *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,full_particle,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00010(const int tree_num,
                                         void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                             int ,
                                                             struct full_particle *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                             int ,
                                                             PS::SPJMonopole *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Monopole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00011(const int tree_num,
                                         void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                             int ,
                                                             struct full_particle *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                             int ,
                                                             PS::SPJQuadrupole *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::Quadrupole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00012(const int tree_num,
                                         void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                             int ,
                                                             struct full_particle *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                             int ,
                                                             PS::SPJMonopoleGeometricCenter *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00013(const int tree_num,
                                         void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                             int ,
                                                             struct full_particle *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                             int ,
                                                             PS::SPJDipoleGeometricCenter *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00014(const int tree_num,
                                         void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                             int ,
                                                             struct full_particle *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                             int ,
                                                             PS::SPJQuadrupoleGeometricCenter *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,full_particle>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00015(const int tree_num,
                                         void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                             int ,
                                                             struct essential_particle_j *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                             int ,
                                                             PS::SPJMonopole *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Monopole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00016(const int tree_num,
                                         void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                             int ,
                                                             struct essential_particle_j *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                             int ,
                                                             PS::SPJQuadrupole *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::Quadrupole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00017(const int tree_num,
                                         void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                             int ,
                                                             struct essential_particle_j *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                             int ,
                                                             PS::SPJMonopoleGeometricCenter *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::MonopoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00018(const int tree_num,
                                         void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                             int ,
                                                             struct essential_particle_j *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                             int ,
                                                             PS::SPJDipoleGeometricCenter *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::DipoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   void calc_force_and_write_back_l00019(const int tree_num,
                                         void (*pfunc_ep_ep)(struct essential_particle_i *,
                                                             int ,
                                                             struct essential_particle_j *,
                                                             int ,
                                                             struct full_particle *),
                                         void (*pfunc_ep_sp)(struct essential_particle_i *,
                                                             int ,
                                                             PS::SPJQuadrupoleGeometricCenter *,
                                                             int ,
                                                             struct full_particle *),
                                         const int psys_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == tree_num) {
            tree = (PS::TreeForForceLong<full_particle,essential_particle_i,essential_particle_j>::QuadrupoleGeometricCenter *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Issue calcForce*()
      tree->calcForceAndWriteBack(pfunc_ep_ep, pfunc_ep_sp, *psys);
   }
   
   //-----------------------------------------------------------------------------------
   // fdps-autogen:get_neighbor_list;

   //----------------------------
   // Utility functions
   //----------------------------
   void mt_init_genrand(const int s) {
      unsigned long seed;
      seed = (unsigned long) s;
      PS::MT::init_genrand(seed);
   }
   int mt_genrand_int31(void) {
      return (int) PS::MT::genrand_int31();
   }
   double mt_genrand_real1(void){
      return PS::MT::genrand_real1();
   } 
   double mt_genrand_real2(void) {
      return PS::MT::genrand_real2(); 
   } 
   double mt_genrand_real3(void) {
      return PS::MT::genrand_real3();
   }
   double mt_genrand_res53(void) {
      return PS::MT::genrand_res53();
   }

   //----------------------------
   // Other utility functions
   //----------------------------
   std::vector<std::string> split(const std::string& s,
                                  const std::string& delim) {
      // This function splits a string by a given delimiter.
      // [Refs.]
      //   1) http://qiita.com/_meki/items/4328c98964ea33b0db0d
      //   2) http://qiita.com/iseki-masaya/items/70b4ee6e0877d12dafa8 and
      std::vector<std::string> results;
      std::string::size_type pos = 0;
      while (pos != std::string::npos) {
         std::string::size_type pos_delim_head = s.find(delim,pos);
         if (pos_delim_head == std::string::npos) { // the end of the target string
            // Register characters between pos and npos as a string
            // return the result to the caller.
            results.push_back(s.substr(pos));
            return results;
         }
         std::string::size_type len = pos_delim_head - pos; 
         results.push_back(s.substr(pos,len));
         pos = pos_delim_head + delim.size();
      }
   }

    void check_split(void) {
       std::string input_string = "interact_type,Force,EPI,EPJ,mode";
       std::string delim = ",";
       std::vector<std::string> results = split(input_string,delim);
       for (std::vector<std::string>::iterator it = results.begin(); it != results.end(); ++it)
          std::cout << *it << delim;
       std::cout << std::endl;
    }


};
