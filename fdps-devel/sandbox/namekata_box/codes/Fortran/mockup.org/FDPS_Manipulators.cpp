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
   //--------------------------------
   // ParticleSystem manipulators
   //--------------------------------
   void create_psys(int *psys_num, char *psys_info) {
      std::string psys_info_ = psys_info;
      Psys_Data psys_data;
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
      psys_vector.push_back(psys_data);
      *psys_num = psys_data.id;
      num_psys_creation++;
      num_psys++;
      ////std::cout << "psys successfully created!" << std::endl;
   }
   void delete_psys(int *psys_num) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == *psys_num) {
            it = psys_vector.erase(it);
            num_psys--;
         }
      }
   }
   void init_psys(int *psys_num) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == *psys_num) {
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
         }
      }
   }
   void get_psys_info(int *psys_num, char *psys_info, size_t *charlen) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == *psys_num) {
            std::strcpy(psys_info,it->info.c_str());
            *charlen  = std::strlen(psys_info);
         }
      }
   }
   void set_nptcl_loc(int *psys_num, int *numPtcl){
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == *psys_num) {
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               psys->setNumberOfParticleLocal(*numPtcl);
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
         }
      }
   }
   int get_nptcl_loc(int *psys_num) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == *psys_num) {
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
         }
      }
   }
   void get_psys_cptr(int *psys_num, void **cptr) {
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == *psys_num) {
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               *cptr = (void *) &((*psys)[0]);
               // [!!IMPORTANT!!]
               //   The operator [0] is necessary to obtain the correct address of
               //   the particle system object.
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
         }
      }
   }
   void exchange_particle(int *psys_num, int *dinfo_num) {
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == *dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Get psys and issue exchangeParticle()
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == *psys_num) {
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
   void delete_dinfo(int *dinfo_num) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == *dinfo_num) {
            it = dinfo_vector.erase(it);
            num_dinfo--;
         }
      }
   }
   void init_dinfo(int *dinfo_num) {
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == *dinfo_num) {
            PS::DomainInfo *dinfo;
            dinfo = (PS::DomainInfo *) it->ptr;
            dinfo->initialize();
         }
      }
   }
   void decompose_domain_all(int *dinfo_num, int *psys_num) {
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == *dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Get psys and issue decomposeDomainAll()
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == *psys_num) {
            if (it->info == "full_particle") {
               PS::ParticleSystem<full_particle> *psys;
               psys = (PS::ParticleSystem<full_particle> *) it->ptr;
               dinfo->decomposeDomainAll(*psys);
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid ParticleSystem number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
         }
      }
   }
   //----------------------------
   // TreeForForce manipulators
   //----------------------------
   void create_tree(int *tree_num, char *tree_info) {
      std::string tree_info_ = tree_info;
      Tree_Data tree_data;
      if (tree_info_ == "Long,full_particle,full_particle,full_particle,Monopole") {
         PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
         tree = new PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole;
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
      tree_vector.push_back(tree_data);
      *tree_num = tree_data.id;
      num_tree_creation++;
      num_tree++;
      //std::cout << "tree sucessfully created!" << std::endl;
   }
   void delete_tree(int *tree_num) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == *tree_num) {
            it = tree_vector.erase(it);
            num_tree--;
         }
      }
   }
   void init_tree(int *tree_num, int *numPtcl) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == *tree_num) {
            if (it->info == "Long,full_particle,full_particle,full_particle,Monopole") {
               PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
               tree = (PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *) it->ptr;
               tree->initialize(*numPtcl);
            } else {
               PS::S32 myrank = PS::Comm::getRank();
               if (myrank == 0) {
                  std::cout << "An invalid TreeForForce number is received." << std::endl;
               }
               PS::Finalize();
               std::exit(EXIT_FAILURE);
            }
         }
      }
   }
   void get_tree_info(int *tree_num, char *tree_info, size_t *charlen) {
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == *tree_num) {
            std::strcpy(tree_info,it->info.c_str());
            *charlen  = std::strlen(tree_info);
         }
      }
   }

   void calc_force_all_and_write_back_s(int *tree_num,
                                        void (*pfunc_ep_ep)(struct full_particle *,
                                                            int ,
                                                            struct full_particle *,
                                                            int ,
                                                            struct full_particle *),
                                        int *psys_num,
                                        int *dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == *tree_num) {
            tree = (PS::TreeForForceLong<full_particle,
                                         full_particle,
                                         full_particle>::Monopole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == *psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == *dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForceAllAndWriteBack()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep,
                                     *psys,*dinfo);
   }

   void calc_force_all_and_write_back_l(int *tree_num,
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
                                        int *psys_num,
                                        int *dinfo_num) {
      // Get tree
      PS::TreeForForceLong<full_particle,full_particle,full_particle>::Monopole *tree;
      for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {
         if (it->id == *tree_num) {
            tree = (PS::TreeForForceLong<full_particle,
                                         full_particle,
                                         full_particle>::Monopole *) it->ptr;
         }
      }
      // Get psys
      PS::ParticleSystem<full_particle> *psys;
      for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {
         if (it->id == *psys_num) {
            psys = (PS::ParticleSystem<full_particle> *) it->ptr;
         }
      }
      // Get dinfo
      PS::DomainInfo *dinfo;
      for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {
         if (it->id == *dinfo_num) {
            dinfo = (PS::DomainInfo *) it->ptr;
         }
      }
      // Issue calcForceAllAndWriteBack()
      tree->calcForceAllAndWriteBack(pfunc_ep_ep,
                                     pfunc_ep_sp,
                                     *psys,*dinfo);
   }

   //----------------------------
   // Utility functions
   //----------------------------
   void mt_init_genrand(int *s) {
      unsigned long seed;
      seed = (unsigned long) *s;
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
