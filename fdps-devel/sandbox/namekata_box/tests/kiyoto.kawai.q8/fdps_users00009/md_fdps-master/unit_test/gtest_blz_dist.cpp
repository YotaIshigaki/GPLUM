//=======================================================================================
//  This is unit test of boltzmann distribution generator.
//=======================================================================================

#include <fstream>
#include <iostream>
#include <vector>

#include <gtest/gtest.h>
#include "boltzmann_dist.hpp"
#include "str_tool.hpp"

#include <random>



class BolzmannDist :
    public ::testing::Test{
    protected:
        MD_EXT::boltzmann_dist blz_dist;

        std::vector<PS::F64> cumulative_dist;
        std::vector<PS::S32> blz_count;
        std::vector<PS::F64> x_list;

        std::vector<PS::F64> cumulative_dist_ref;
        std::vector<PS::S32> blz_count_ref;
        std::vector<PS::F64> x_list_ref;

        PS::F64 range_min, dv;

        const int64_t N_data = 400000;
        const PS::F64 eps    = 1.e-10;

        virtual void SetUp(){
            std::mt19937 mt;
            std::uniform_real_distribution<> dist_real(0.0, 1.0); // must be [0.0, 1.0)
            mt.seed(54321);

            const std::string out_file_name = "test_bin/blz_dist_result.dat";
            const std::string ref_file_name = "unit_test/ref/blz_dist_result_ref.dat";

            //--- internal result
            //------ get internal table value
            blz_dist.get_cumulative_dist(cumulative_dist, range_min, dv);
            this->x_list.clear();
            for(size_t i=0; i<cumulative_dist.size(); ++i){
                this->x_list.push_back(range_min + PS::F64(i)*dv);
            }

            //------ count result
            blz_count.resize(cumulative_dist.size());
            for(int64_t i=0; i<N_data; ++i){
                PS::F64 blz_value = blz_dist( PS::F64(dist_real(mt)) );

                PS::S32 index    =  PS::S32( (blz_value - range_min)/dv );
                blz_count[index] += 1;
            }

            //------ record internal result
            std::ofstream file{out_file_name};
            file << "x "
                 << "blz_count "
                 << "blz_cumulative" << std::endl;
            for(size_t i=0; i<cumulative_dist.size(); ++i){
                file << std::scientific << std::setprecision(15) << x_list[i]          << " "
                     << std::dec        << std::setprecision(10) << blz_count[i]       << " "
                     << std::scientific << std::setprecision(15) << cumulative_dist[i] << std::endl;
            }
            file.close();

            //--- load reference result
            this->x_list_ref.clear();
            this->blz_count_ref.clear();
            this->cumulative_dist_ref.clear();

            std::ifstream file_ref{ref_file_name};
            if(file_ref.fail()) throw std::ios_base::failure("reference data: " + ref_file_name + " was not found.");

            std::string              line;
            std::vector<std::string> str_list;
            while( getline(file_ref, line) ){
                STR_TOOL::removeCR(line);
                str_list = STR_TOOL::split(line, " ");

                if(str_list.size() < 3) continue;
                if( not STR_TOOL::isNumeric(str_list[0]) ) continue;
                if( not STR_TOOL::isInteger(str_list[1]) ) continue;
                if( not STR_TOOL::isNumeric(str_list[2]) ) continue;

                this->x_list_ref.push_back(          std::stod(str_list[0]) );
                this->blz_count_ref.push_back(       std::stoi(str_list[1]) );
                this->cumulative_dist_ref.push_back( std::stod(str_list[2]) );
            }
        }
};

PS::F64 rel_diff(const PS::F64 lhs, const PS::F64 rhs){
    if(lhs == 0.0 && rhs == 0.0) return 0.0;

    PS::F64 diff = lhs - rhs;
    if(lhs == 0.0){
        return diff/rhs;
    } else {
        return diff/lhs;
    }
}

//--- unit test definition, CANNOT use "_" in test/test_case name.
TEST_F(BolzmannDist, result){
    ASSERT_EQ(x_list.size(),    x_list_ref.size());
    ASSERT_EQ(blz_count.size(), blz_count_ref.size());

    for(size_t i=0; i<x_list.size(); ++i){
        ASSERT_TRUE(rel_diff(x_list[i], x_list_ref[i]) < eps) << "result: " << x_list[i]
                                                              << " | ref: " << x_list_ref[i];
        ASSERT_EQ(blz_count[i], blz_count_ref[i]);
    }
}
TEST_F(BolzmannDist, internalTable){
    ASSERT_EQ(x_list.size(),          x_list_ref.size());
    ASSERT_EQ(cumulative_dist.size(), cumulative_dist_ref.size());

    for(size_t i=0; i<x_list.size(); ++i){
        ASSERT_TRUE(rel_diff(x_list[i], x_list_ref[i]) < eps)
            << "result: " << x_list[i]
            << " | ref: " << x_list_ref[i];
        ASSERT_TRUE(rel_diff(cumulative_dist[i], cumulative_dist_ref[i]) < eps)
            << "result: " << cumulative_dist[i]
            << " | ref: " << cumulative_dist_ref[i];
    }
}

#include "gtest_main.hpp"
