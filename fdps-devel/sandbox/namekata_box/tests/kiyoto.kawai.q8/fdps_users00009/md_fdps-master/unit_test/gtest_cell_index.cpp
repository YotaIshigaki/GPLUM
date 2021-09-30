//=======================================================================================
//  This is unit test of MD_EXT::CellIndex.
//     module location: ./generic_ext/cell_index.hpp
//=======================================================================================

#include <gtest/gtest.h>
#include <particle_simulator.hpp>

#include <random>

#include <cell_index.hpp>

//--- unit test definition, CANNOT use "_" in test/test_case name.
TEST(CellIndex, init){
    MD_EXT::CellIndex<int> cell_index;

    cell_index.initDomain( PS::F32vec{-1.0, -1.0, -1.0},
                           PS::F32vec{ 2.0,  2.0,  2.0} );
    EXPECT_EQ(cell_index.get_domain(), std::make_pair(PS::F32vec{-1.0, -1.0, -1.0},
                                                      PS::F32vec{ 2.0,  2.0,  2.0}));

    cell_index.initIndex(3, 6, 12);
    EXPECT_EQ(cell_index.get_grid_size().x, 3);
    EXPECT_EQ(cell_index.get_grid_size().y, 6);
    EXPECT_EQ(cell_index.get_grid_size().z, 12);
    EXPECT_FLOAT_EQ(cell_index.get_cell_size().x, 1.0);
    EXPECT_FLOAT_EQ(cell_index.get_cell_size().y, 0.5);
    EXPECT_FLOAT_EQ(cell_index.get_cell_size().z, 0.25);

    cell_index.init( PS::F32vec{0.0, 0.0, 0.0},
                     PS::F32vec{2.0, 4.0, 3.0},
                     0.24 );
    EXPECT_EQ(cell_index.get_grid_size().x, 10);
    EXPECT_EQ(cell_index.get_grid_size().y, 18);
    EXPECT_EQ(cell_index.get_grid_size().z, 14);
    EXPECT_FLOAT_EQ(cell_index.get_cell_size().x, 2.0/10.0);
    EXPECT_FLOAT_EQ(cell_index.get_cell_size().y, 4.0/18.0);
    EXPECT_FLOAT_EQ(cell_index.get_cell_size().z, 3.0/14.0);
}

TEST(CellIndex, neighborList){
    const int seed = 19937;
    const int n    = 8192;

    const double l_domain = 100.0;
    const double r_search = 9.0;

    std::mt19937 mt;
    std::uniform_real_distribution<> dist_real(0.0, l_domain);
    mt.seed(seed);

    //--- test target
    std::vector<PS::F32vec> pos;

    MD_EXT::CellIndex<size_t>        cell_index;
    std::vector<std::vector<size_t>> pair_list_ref;

    //--- make pos data
    pos.clear();
    for(int i=0; i<n; ++i){
        pos.push_back( PS::F32vec{ static_cast<PS::F32>(dist_real(mt)),
                                   static_cast<PS::F32>(dist_real(mt)),
                                   static_cast<PS::F32>(dist_real(mt)) } );
    }

    //--- make 1D pair list (direct method, reference data)
    pair_list_ref.clear();
    pair_list_ref.resize(n);
    const double r2_search = r_search*r_search;
    for(int i=0; i<n; ++i){
        pair_list_ref.at(i).clear();
        for(int j=0; j<n; ++j){
            if(j == i) continue;
            PS::F32vec r_pos = pos.at(j) - pos.at(i);
            PS::F32    r2    = r_pos*r_pos;
            if(r2 <= r2_search){
                pair_list_ref.at(i).push_back(j);
            }
        }
    }

    //--- init cell index
    cell_index.init( PS::F32vec{0.0, 0.0, 0.0},
                     PS::F32vec{l_domain, l_domain, l_domain},
                     r_search );

    //--- make index table in cell
    for(int i=0; i<n; ++i){
        cell_index.add(pos.at(i), i);
    }

    //--- compare result (internal buffer)
    for(int i=0; i<n; ++i){
        const auto& list_ref  = pair_list_ref.at(i);
        const auto  cell_data = cell_index.get_data_list(pos.at(i), r_search);
        for(const auto tgt : list_ref){
            size_t n_tgt = std::count(cell_data.begin(), cell_data.end(), tgt);
            EXPECT_EQ(n_tgt, 1) << " i= " << i << " tgt= " << tgt;
        }
    }

    //--- compare result (const interface for outside buffer)
    std::vector<MD_EXT::CellIndex<size_t>::index_type> cell_list;
    std::vector<MD_EXT::CellIndex<size_t>::value_type> cell_data;
    for(int i=0; i<n; ++i){
        const auto& list_ref = pair_list_ref.at(i);
        cell_list.clear();
        cell_data.clear();
        cell_index.get_data_list(pos.at(i), r_search, cell_list, cell_data);
        for(const auto tgt : list_ref){
            size_t n_tgt = std::count(cell_data.begin(), cell_data.end(), tgt);
            EXPECT_EQ(n_tgt, 1) << " i= " << i << " tgt= " << tgt;
        }
    }

}

#include "gtest_main.hpp"
