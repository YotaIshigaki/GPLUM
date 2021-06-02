//=======================================================================================
//  This is unit test of additional wrapper of PS::Comm::.
//     provides the interface for some STL container.
//     module location: ./generic_ext/comm_tool.hpp
//=======================================================================================

#include <gtest/gtest.h>
#include <particle_simulator.hpp>
#include "comm_tool_allGather.hpp"

#include <random>



//==========================================
// MPI allGather
//==========================================
struct DataBasic {
    std::vector<int>              vec_int;
    std::vector<std::vector<int>> vec_vec_int;
    std::vector<std::string>      vec_str;
    std::vector<std::pair<int, float>>  vec_pair_i_f;

    DataBasic() = default;
    ~DataBasic() = default;

    void clear(){
        this->vec_int.clear();
        this->vec_vec_int.clear();
        this->vec_str.clear();
        this->vec_pair_i_f.clear();
    }

    void generate(const int seed, const size_t N_data){
        std::mt19937 mt;
        std::uniform_int_distribution<>  dist_int(0,10);
        std::uniform_real_distribution<> dist_real(-99.9,99.9);
        mt.seed(seed);

        this->clear();

        std::vector<int> buff_vi;
        std::string      buff_str;

        //--- for basic pattern
        buff_vi.clear();
        buff_str = "";
        for(size_t i=0; i<N_data; ++i){
            int   data_int  = dist_int(mt);
            float data_real = dist_real(mt);

            this->vec_int.push_back(data_int);

            if(data_int == 10){
                this->vec_vec_int.push_back(buff_vi);
                this->vec_str.push_back(buff_str);
                buff_vi.clear();
                buff_str = "";
            } else {
                buff_vi.push_back(data_int);
                buff_str += std::to_string(data_int);
            }

            this->vec_pair_i_f.push_back( std::make_pair( data_int, data_real ) );
        }
    }
};

class AllGatherBasic :
    public ::testing::Test{
    protected:
        //--- for basic data pattern
        std::vector<DataBasic> data;
        int n_proc;
        int id_proc;

        virtual void SetUp(){
            n_proc  = PS::Comm::getNumberOfProc();
            id_proc = PS::Comm::getRank();

            size_t N_data = 10000;

            this->data.resize(n_proc);
            for(int i=0; i<n_proc; ++i){
                int seed = 19937*(1 + i);
                this->data[i].generate(seed, N_data);
            }
        }
};

//--- unit test definition, CANNOT use "_" in test/test_case name.
TEST_F(AllGatherBasic, VecInt){
    auto recv_vec_vi     = COMM_TOOL::allGather(data[id_proc].vec_int     );
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(data.at(i).vec_int     , recv_vec_vi.at(i)     ) << "source_proc = " << i;
    }
}

TEST_F(AllGatherBasic, VecVecInt){
    auto recv_vec_vvi    = COMM_TOOL::allGather(data[id_proc].vec_vec_int );
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(data.at(i).vec_vec_int , recv_vec_vvi.at(i)    ) << "source_proc = " << i;
    }
}

TEST_F(AllGatherBasic, VecStr){
    auto recv_vec_vs     = COMM_TOOL::allGather(data[id_proc].vec_str     );
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(data.at(i).vec_str     , recv_vec_vs.at(i)     ) << "source_proc = " << i;
    }
}

TEST_F(AllGatherBasic, VecPairIntFloat){
    auto recv_vec_vp_i_f = COMM_TOOL::allGather(data[id_proc].vec_pair_i_f);
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(data.at(i).vec_pair_i_f, recv_vec_vp_i_f.at(i) ) << "source_proc = " << i;
    }
}



struct DataRecursive {
    std::vector<std::vector<std::vector<int>>> vec_vec_vec_int;
    std::vector<std::vector<std::string>>      vec_vec_str;
    std::vector<std::pair<std::string, std::vector<int>>>  vec_pair_s_vi;

    DataRecursive() = default;
    ~DataRecursive() = default;

    void clear(){
        this->vec_vec_vec_int.clear();
        this->vec_vec_str.clear();
        this->vec_pair_s_vi.clear();
    }

    void generate(const int seed, const size_t N_data){
        std::mt19937 mt;
        std::uniform_int_distribution<>  dist_int(0,11);
        mt.seed(seed);

        this->clear();

        std::vector<std::vector<int>> buff_v_vi;
        std::vector<std::string>      buff_v_str;
        std::vector<int> buff_vi;
        std::string      buff_str;

        buff_v_vi.clear();
        buff_v_str.clear();
        buff_vi.clear();
        buff_str = "";
        for(size_t i=0; i<N_data; ++i){
            int data = dist_int(mt);

            if(       data == 11){
                this->vec_vec_vec_int.push_back(buff_v_vi);
                this->vec_vec_str.push_back(buff_v_str);

                buff_v_vi.clear();
                buff_v_str.clear();

            } else if(data == 10){
                buff_v_vi.push_back(buff_vi);
                buff_v_str.push_back(buff_str);

                this->vec_pair_s_vi.push_back( std::make_pair( buff_str,
                                                               buff_vi ) );

                buff_vi.clear();
                buff_str = "";

            } else {
                buff_vi.push_back(data);
                buff_str += std::to_string(data);
            }
        }
    }
};

class AllGatherRecursive :
    public ::testing::Test{
    protected:
        //--- for basic data pattern
        std::vector<DataRecursive> data;
        int n_proc;
        int id_proc;

        virtual void SetUp(){
            n_proc  = PS::Comm::getNumberOfProc();
            id_proc = PS::Comm::getRank();

            size_t N_data = 10000;

            this->data.resize(n_proc);
            for(int i=0; i<n_proc; ++i){
                int seed = 19937*(1 + i);
                this->data[i].generate(seed, N_data);
            }
        }
};

//--- unit test definition, CANNOT use "_" in test/test_case name.
TEST_F(AllGatherRecursive, VecVecVecInt){
    auto recv_vec_vvvi    = COMM_TOOL::allGather(data[id_proc].vec_vec_vec_int);
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(data.at(i).vec_vec_vec_int, recv_vec_vvvi.at(i)    ) << "source_proc = " << i;
    }
}

TEST_F(AllGatherRecursive, VecVecStr){
    auto recv_vec_vvs     = COMM_TOOL::allGather(data[id_proc].vec_vec_str    );
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(data.at(i).vec_vec_str    , recv_vec_vvs.at(i)     ) << "source_proc = " << i;
    }
}

TEST_F(AllGatherRecursive, VecPairStrVecInt){
    auto recv_vec_vp_s_vi = COMM_TOOL::allGather(data[id_proc].vec_pair_s_vi  );
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(data.at(i).vec_pair_s_vi  , recv_vec_vp_s_vi.at(i) ) << "source_proc = " << i;
    }
}

#include "gtest_main_mpi.hpp"
