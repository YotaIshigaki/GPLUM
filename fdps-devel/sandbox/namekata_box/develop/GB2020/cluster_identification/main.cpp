// Include the standard C++ headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cstdint>
#include <vector>
// Include the header file of FDPS
#include <particle_simulator.hpp>
// Include user-defined headers
#include "user_defined.hpp"
#include "io.hpp"



int main(int argc, char *argv[]) {

    // Initialize FDPS
    PS::Initialize(argc, argv);

    // Make an instance of ParticleSystem and initialize it
    PS::ParticleSystem<FP_nbody> psys;
    psys.initialize();

    // Make an instance of DomainInfo and initialize it
    PS::DomainInfo dinfo;
    dinfo.initialize();

    // Read particle data
    const std::string file_name_base = "./ptcl_data_";
    readParticleData(psys, file_name_base);

    // Perform domain decomposition 
    dinfo.decomposeDomainAll(psys);
    
    // Perform particle exchange
    psys.exchangeParticle(dinfo);

    // Make tree structures
    const PS::S32 n_loc = psys.getNumberOfParticleLocal();
    PS::TreeForForceShort<Force_clus_id, EP_clus_id, EP_clus_id>::Symmetry tree;
    tree.initialize(3 * n_loc);


    // Perform force calculation
    tree.calcForceAllAndWriteBack(CalcForceEmpty,
                                  psys, dinfo);

    // Find cluster
    //const PS::F64 dt_crit = 3.0e4;
    //const PS::F64 dt_crit = 2.0e4;
    const PS::F64 dt_crit = 1.0e4;
    std::vector<PS::ClusterInfo> cluster_list = tree.getClusterList(dt_crit);

    // Sort cluster list
    std::sort(cluster_list.begin(), cluster_list.end(),  
              [](const PS::ClusterInfo & l, const PS::ClusterInfo & r)
                  -> bool {return l.id < r.id;} ); 


    // Output 
    if (PS::Comm::getRank() == 0) {
        const std::string file_name = "cluster_data";
        std::ofstream output_file;
        output_file.open(file_name.c_str(), std::ios::trunc);
        for (const auto & e : cluster_list) {
            output_file << e.vertex << "   "
                        << std::oct << e.id << "   ";
            for (PS::S32 i = 0; i < e.ranks.size(); i++)
                output_file << std::dec << e.ranks[i] << "   ";
            output_file << std::endl;
        }
        output_file.close();
    }

    // Test to create a communicator
    if (PS::Comm::getRank() == 0)
        std::cout << "start a test to create a communicator." << std::endl;
    MPI_Group parent_group;
    MPI_Comm_group(MPI_COMM_WORLD,&parent_group);
    for (const auto & e : cluster_list) {
        if (0615115255451604045676 == e.id) {
            MPI_Group group;
            MPI_Comm comm;
            bool is_contained {false};
            for (PS::S32 i = 0; i < e.ranks.size(); i++) 
                if (PS::Comm::getRank() == e.ranks[i]) is_contained = true;
            if (is_contained) {
                MPI_Group_incl(parent_group, e.ranks.size(), &e.ranks[0], &group);
                const PS::S32 tag = 0;
                const PS::S32 ret = MPI_Comm_create_group(MPI_COMM_WORLD, group, tag, &comm);
                if (ret) std::cout << "Something wrong occurs." << std::endl;
                else std::cout << "A communicator is created @ rank = " << PS::Comm::getRank() << std::endl;
            }
        }
    }

    // Output particle distribution whose dt < dt_crit
    {
        std::stringstream ss;
        ss << "pos_data_" << std::setw(5) << std::setfill('0') << PS::Comm::getRank() << ".txt";
        const std::string file_name = ss.str();
        std::ofstream output_file;
        output_file.open(file_name.c_str(), std::ios::trunc);
        for (PS::S32 i = 0; i < n_loc; i++) {
            if (psys[i].dt < dt_crit) {
                output_file << psys[i].pos << "   " << psys[i].type << std::endl;
            }
        }
        output_file.close();
    }

    // Finalize FDPS
    PS::Finalize();
    return 0;
}
