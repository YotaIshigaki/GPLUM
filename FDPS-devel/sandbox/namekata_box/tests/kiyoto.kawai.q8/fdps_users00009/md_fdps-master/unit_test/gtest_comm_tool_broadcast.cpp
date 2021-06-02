//=======================================================================================
//  This is unit test of additional wrapper of PS::Comm::.
//     provides the interface for some STL container.
//     module location: ./generic_ext/comm_tool.hpp
//=======================================================================================

#include <gtest/gtest.h>
#include <particle_simulator.hpp>
#include "comm_tool_broadcast.hpp"

#include <random>



//==========================================
// MPI broadcast
//==========================================
struct DataBasic {
    std::vector<int>              vec_int;
    std::vector<std::vector<int>> vec_vec_int;
    std::vector<std::string>      vec_str;
    std::vector<std::pair<int, float>>  vec_pair_i_f;
    std::unordered_map<int, float>      map_i_f;
    std::unordered_multimap<int, float> m_map_i_f;

    DataBasic() = default;
    ~DataBasic() = default;

    void clear(){
        this->vec_int.clear();
        this->vec_vec_int.clear();
        this->vec_str.clear();
        this->vec_pair_i_f.clear();
        this->map_i_f.clear();
        this->m_map_i_f.clear();
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

            this->map_i_f[data_int]   = data_real;
            this->m_map_i_f.insert( std::make_pair(data_int, data_real) );
        }
    }
};

class BroadcastBasic :
    public ::testing::Test{
    protected:
        //--- for basic data pattern
        std::vector<DataBasic> data;
        std::vector<DataBasic> ref;

        int n_proc;

        virtual void SetUp(){
            n_proc  = PS::Comm::getNumberOfProc();
            size_t N_data = 10000;

            this->data.resize(n_proc);
            this->ref.resize(n_proc);
            for(int i=0; i<n_proc; ++i){
                int seed = 19937*(1 + i);
                this->data[i].generate(seed, N_data);
                this->ref[i].generate(seed, N_data);
            }
        }
};

//--- unit test definition, CANNOT use "_" in test/test_case name.
TEST_F(BroadcastBasic, VecInt){
    for(int i=0; i<n_proc; ++i){
        COMM_TOOL::broadcast(data[i].vec_int     , i);
    }
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(ref[i].vec_int     , data[i].vec_int     ) << "source_proc = " << i;
    }
}
TEST_F(BroadcastBasic, VecVecInt){
    for(int i=0; i<n_proc; ++i){
        COMM_TOOL::broadcast(data[i].vec_vec_int , i);
    }
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(ref[i].vec_vec_int , data[i].vec_vec_int ) << "source_proc = " << i;
    }
}
TEST_F(BroadcastBasic, VecStr){
    for(int i=0; i<n_proc; ++i){
        COMM_TOOL::broadcast(data[i].vec_str     , i);
    }
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(ref[i].vec_str     , data[i].vec_str     ) << "source_proc = " << i;
    }
}
TEST_F(BroadcastBasic, VecPairIntFloat){
    for(int i=0; i<n_proc; ++i){
        COMM_TOOL::broadcast(data[i].vec_pair_i_f, i);
    }
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(ref[i].vec_pair_i_f, data[i].vec_pair_i_f) << "source_proc = " << i;
    }
}
TEST_F(BroadcastBasic, MapIntFloat){
    for(int i=0; i<n_proc; ++i){
        COMM_TOOL::broadcast(data[i].map_i_f     , i);
    }
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(ref[i].map_i_f     , data[i].map_i_f     ) << "source_proc = " << i;
    }
}
TEST_F(BroadcastBasic, MultiMapIntFloat){
    for(int i=0; i<n_proc; ++i){
        COMM_TOOL::broadcast(data[i].m_map_i_f   , i);
    }
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(ref[i].m_map_i_f   , data[i].m_map_i_f   ) << "source_proc = " << i;
    }
}



struct DataRecursive {
    std::vector<std::vector<std::vector<int>>> vec_vec_vec_int;
    std::vector<std::vector<std::string>>      vec_vec_str;
    std::vector<std::pair<std::string, std::vector<int>>>  vec_pair_s_vi;
    std::unordered_map<std::string, std::vector<int>>      map_s_vi;
    std::unordered_multimap<std::string, std::vector<int>> m_map_s_vi;

    DataRecursive() = default;
    ~DataRecursive() = default;

    void clear(){
        this->vec_vec_vec_int.clear();
        this->vec_vec_str.clear();
        this->vec_pair_s_vi.clear();
        this->map_s_vi.clear();
        this->m_map_s_vi.clear();
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
                this->map_s_vi.insert( std::make_pair( buff_str,
                                                       buff_vi ) );
                this->m_map_s_vi.insert( std::make_pair( buff_str,
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

class BroadcastRecursive :
    public ::testing::Test{
    protected:
        //--- for basic data pattern
        std::vector<DataRecursive> data;
        std::vector<DataRecursive> ref;

        int n_proc;

        virtual void SetUp(){
            n_proc  = PS::Comm::getNumberOfProc();

            size_t N_data = 10000;

            this->data.resize(n_proc);
            this->ref.resize(n_proc);
            for(int i=0; i<n_proc; ++i){
                int seed = 19937*(1 + i);
                this->data[i].generate(seed, N_data);
                this->ref[i].generate(seed, N_data);
            }
        }
};

//--- unit test definition, CANNOT use "_" in test/test_case name.
TEST_F(BroadcastRecursive, VecVecVecInt){
    for(int i=0; i<n_proc; ++i){
        COMM_TOOL::broadcast(data[i].vec_vec_vec_int, i);
    }
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(ref[i].vec_vec_vec_int, data[i].vec_vec_vec_int) << "source_proc = " << i;
    }
}
TEST_F(BroadcastRecursive, VecVecStr){
    for(int i=0; i<n_proc; ++i){
        COMM_TOOL::broadcast(data[i].vec_vec_str    , i);
    }
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(ref[i].vec_vec_str    , data[i].vec_vec_str    ) << "source_proc = " << i;
    }
}
TEST_F(BroadcastRecursive, VecPairStrVecInt){
    for(int i=0; i<n_proc; ++i){
        COMM_TOOL::broadcast(data[i].vec_pair_s_vi  , i);
    }
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(ref[i].vec_pair_s_vi  , data[i].vec_pair_s_vi  ) << "source_proc = " << i;
    }
}
TEST_F(BroadcastRecursive, MapStrVecInt){
    for(int i=0; i<n_proc; ++i){
        COMM_TOOL::broadcast(data[i].map_s_vi       , i);
    }
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(ref[i].map_s_vi       , data[i].map_s_vi       ) << "source_proc = " << i;
    }
}
TEST_F(BroadcastRecursive, MultiMapStrVecInt){
    for(int i=0; i<n_proc; ++i){
        COMM_TOOL::broadcast(data[i].m_map_s_vi     , i);
    }
    for(int i=0; i<n_proc; ++i){
        ASSERT_EQ(ref[i].m_map_s_vi     , data[i].m_map_s_vi     ) << "source_proc = " << i;
    }
}

#include "gtest_main_mpi.hpp"
