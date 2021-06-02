#include<iostream>
#include<random>
#include<particle_simulator.hpp>
#include"gtest/gtest.h"

class CommInfoTest : public ::testing::Test{
protected:
    PS::S32 rank_0_;
    PS::S32 rank_1_;
    PS::S32 rank_2_;
    PS::S32 rank_3_;
    PS::S32 n_proc_0_;
    PS::S32 n_proc_1_;
    PS::S32 n_proc_2_;
    PS::S32 n_proc_3_;
    PS::CommInfo comm_info_0_;
    PS::CommInfo comm_info_1_;
    PS::CommInfo comm_info_2_;
    PS::CommInfo comm_info_3_;
    
    virtual void SetUp(){
        comm_info_0_ = PS::Comm::getCommInfo();
        rank_0_ = PS::Comm::getRank();
        n_proc_0_ = PS::Comm::getNumberOfProc();
        
        comm_info_1_ = PS::Comm::split(rank_0_%2, rank_0_);
        rank_1_ = comm_info_1_.getRank();
        n_proc_1_ = comm_info_1_.getNumberOfProc();

        comm_info_2_ = comm_info_1_.split(rank_1_%2, rank_1_);
        rank_2_ = comm_info_2_.getRank();
        n_proc_2_ = comm_info_2_.getNumberOfProc();

        comm_info_3_ = comm_info_2_.split(rank_2_%2, rank_2_);
        rank_3_ = comm_info_3_.getRank();
        n_proc_3_ = comm_info_3_.getNumberOfProc();
    }
    virtual void TearDown(){}
};

class CommInfoTest_create : public ::testing::Test{
protected:
    PS::S32 rank_0_;
    PS::S32 rank_1_;
    PS::S32 rank_2_;
    PS::S32 rank_3_;
    PS::S32 n_proc_0_;
    PS::S32 n_proc_1_;
    PS::S32 n_proc_2_;
    PS::S32 n_proc_3_;
    PS::CommInfo comm_info_0_;
    PS::CommInfo comm_info_1_;
    PS::CommInfo comm_info_2_;
    PS::CommInfo comm_info_3_;

    template<typename Tcomminfo>
    PS::CommInfo my_create(const Tcomminfo & comm_info_old){
        std::vector<PS::S32 >ranks;
        for(int i=0; i<comm_info_old.getNumberOfProc(); i++){
            if(i%2==0){
                ranks.push_back(i);
            }
        }
        return comm_info_old.create(ranks.size(), &ranks[0]);
    }
    virtual void SetUp(){
        comm_info_0_ = PS::Comm::getCommInfo();
        comm_info_1_ = my_create(PS::Comm::getCommInfo());
        if(comm_info_1_.isNotCommNull()){
            comm_info_2_ = my_create(comm_info_1_);
            if(comm_info_2_.isNotCommNull()){
                comm_info_3_ = my_create(comm_info_2_);
            }
        }

        n_proc_0_ = PS::Comm::getNumberOfProc();
        rank_0_ = PS::Comm::getRank();        
        n_proc_1_ = comm_info_1_.getNumberOfProc();
        rank_1_ = comm_info_1_.getRank();        
        n_proc_2_ = comm_info_2_.getNumberOfProc();
        rank_2_   = comm_info_2_.getRank();
        n_proc_3_ = comm_info_3_.getNumberOfProc();
        rank_3_   = comm_info_3_.getRank();
    }
    virtual void TearDown(){}
};


TEST_F(CommInfoTest, basic){
    EXPECT_EQ(rank_0_, PS::Comm::getRank());
    std::cerr<<"rank_0_= "<<rank_0_;
    if(comm_info_1_.isNotCommNull()){
        std::cerr<<" rank_1_= "<<rank_1_;
        if(comm_info_2_.isNotCommNull()){
            std::cerr<<" rank_2_= "<<rank_2_;
            if(comm_info_3_.isNotCommNull()){
                std::cerr<<" rank_3_= "<<rank_3_;
            }
        }
    }
    std::cerr<<std::endl;
    EXPECT_EQ(rank_1_, rank_0_/2);
    comm_info_0_.free();
    comm_info_1_.free();
    comm_info_2_.free();
    comm_info_3_.free();
}


TEST_F(CommInfoTest_create, basic){
    EXPECT_EQ(rank_0_, PS::Comm::getRank());
    if(comm_info_1_.isNotCommNull()){
        std::cerr<<" rank_1_= "<<rank_1_;
        EXPECT_EQ(rank_1_, rank_0_/2);
        if(comm_info_2_.isNotCommNull()){
            std::cerr<<" rank_2_= "<<rank_2_;
            EXPECT_EQ(rank_2_, rank_1_/2);
            if(comm_info_3_.isNotCommNull()){
                EXPECT_EQ(rank_3_, rank_2_/2);
                std::cerr<<" rank_3_= "<<rank_3_;
            }
        }
    }
    std::cerr<<std::endl;
    comm_info_0_.free();
    comm_info_1_.free();
    comm_info_2_.free();
    comm_info_3_.free();
}

TEST_F(CommInfoTest, logical){
    bool loc = false;
    if(rank_1_==0){ loc = true; }
    bool glb_or  = comm_info_1_.synchronizeConditionalBranchOR(loc);
    EXPECT_EQ(glb_or, true);
    bool glb_and  = comm_info_1_.synchronizeConditionalBranchAND(loc);
    EXPECT_EQ(glb_and, false);
}

TEST_F(CommInfoTest, bcast){
    PS::F64vec val(-1.0, -2.0, -3.0);
    if(rank_1_==0){ val = PS::F64vec(1.0, 2.0, 3.0); }
    comm_info_1_.broadcast(&val, 1);
    EXPECT_EQ(val.x, 1.0);
    EXPECT_EQ(val.y, 2.0);
    EXPECT_EQ(val.z, 3.0);
}

TEST_F(CommInfoTest, gather){
    const int dst = 0;
    const auto n_send = rank_1_ + 1;
    PS::ReallocatableArray<PS::F64vec> val_send(n_send, n_send, PS::MemoryAllocMode::Pool);
    for(auto i=0; i<rank_1_+1; i++) val_send[i] = PS::F64vec(rank_1_);
    PS::ReallocatableArray<PS::F64vec> val_recv(PS::MemoryAllocMode::Pool);
    comm_info_1_.gatherVAll(val_send, n_send, val_recv, dst);
    if(rank_1_ == dst){
        int n_cnt = 0;
        for(auto i=0; i<n_proc_1_; i++){
            for(auto j=0; j<i+1; j++){
                EXPECT_EQ(val_recv[n_cnt].x, (PS::F64)i);
                EXPECT_EQ(val_recv[n_cnt].y, (PS::F64)i);
                EXPECT_EQ(val_recv[n_cnt].z, (PS::F64)i);
                n_cnt++;
            }
        }
    }
}

TEST_F(CommInfoTest, scatter){
    const int src = 0;
    PS::ReallocatableArray<int> n_send(n_proc_1_, n_proc_1_, PS::MemoryAllocMode::Pool);
    int n_send_tot = 0;
    for(auto i=0; i<n_proc_1_; i++){
        n_send[i] = i+1;
        n_send_tot += n_send[i];
    }
    PS::ReallocatableArray<PS::F64vec> val_send(n_send_tot, n_send_tot, PS::MemoryAllocMode::Pool);
    int n_cnt = 0;
    for(auto i=0; i<n_proc_1_; i++){
        for(auto j=0; j<i+1; j++){
            val_send[n_cnt] = PS::F64vec(i);
            n_cnt++;
        }
    }
    PS::ReallocatableArray<PS::F64vec> val_recv(PS::MemoryAllocMode::Pool);
    comm_info_1_.scatterVAll(val_send, n_send.getPointer(), val_recv, src);
    const auto n_recv = val_recv.size();
    EXPECT_EQ(n_recv, rank_1_+1);
    EXPECT_EQ(n_recv, val_recv.size());
    for(auto i=0; i<n_recv; i++){
        EXPECT_EQ(val_recv[i].x, (PS::F64)rank_1_);
        EXPECT_EQ(val_recv[i].y, (PS::F64)rank_1_);
        EXPECT_EQ(val_recv[i].z, (PS::F64)rank_1_);
    }
}

TEST_F(CommInfoTest, allgather){
    const auto n_send = rank_1_ + 1;
    PS::ReallocatableArray<PS::F64vec> val_send(n_send, n_send, PS::MemoryAllocMode::Pool);
    for(auto i=0; i<n_send; i++) val_send[i] = PS::F64vec(rank_1_);
    PS::ReallocatableArray<PS::F64vec> val_recv(PS::MemoryAllocMode::Pool);
    comm_info_1_.allGatherVAll(val_send, n_send, val_recv);
    int n_cnt = 0;
    for(auto i=0; i<n_proc_1_; i++){
        for(auto j=0; j<i+1; j++){
            EXPECT_EQ(val_recv[n_cnt].x, (PS::F64)i);
            EXPECT_EQ(val_recv[n_cnt].y, (PS::F64)i);
            EXPECT_EQ(val_recv[n_cnt].z, (PS::F64)i);
            n_cnt++;
        }
    }
}

TEST_F(CommInfoTest, alltoall){
    PS::ReallocatableArray<PS::F64vec> val_send(PS::MemoryAllocMode::Pool);
    PS::ReallocatableArray<int> n_send(n_proc_1_, n_proc_1_, PS::MemoryAllocMode::Pool);
    PS::ReallocatableArray<int> n_disp_send(n_proc_1_+1, n_proc_1_+1, PS::MemoryAllocMode::Pool);
    n_disp_send[0] = 0;
    for(auto i=0; i<n_proc_1_; i++){
        n_send[i] = i+1;
        n_disp_send[i+1] = n_disp_send[i] + n_send[i];
        for(auto j=0; j<n_send[i]; j++){
            val_send.push_back(PS::F64vec(i));
        }
    }
    PS::ReallocatableArray<PS::F64vec> val_recv(PS::MemoryAllocMode::Pool);
    comm_info_1_.allToAllVAll(val_send, n_send.getPointer(), val_recv);
    for(auto i=0; i<val_recv.size(); i++){
        EXPECT_EQ(val_recv[i].x, (PS::F64)rank_1_);
        EXPECT_EQ(val_recv[i].y, (PS::F64)rank_1_);
        EXPECT_EQ(val_recv[i].z, (PS::F64)rank_1_);
    }
}

TEST_F(CommInfoTest, getmin){
    float loc_f = (float)rank_1_+10.0;
    auto glb_f = comm_info_1_.getMinValue(loc_f);
    EXPECT_EQ(glb_f, (float)0+10.0);

    PS::F64vec loc_v = PS::F64vec((PS::F64)(rank_1_+1)*10.0,
                                  -100.0*(PS::F64)(rank_1_+1),
                                  (PS::F64)(rank_1_+1)*1000.0);
    auto glb_v = comm_info_1_.getMinValue(loc_v);
    EXPECT_EQ(glb_v.x, 10.0);
    EXPECT_EQ(glb_v.y, -100.0*(PS::F64)(n_proc_1_));
    EXPECT_EQ(glb_v.z, 1000.0);
}

TEST_F(CommInfoTest, getsum){
    float loc_f = (float)rank_1_+10.0;
    float * loc_recv = new float[n_proc_1_];
    auto glb_f = comm_info_1_.getSum(loc_f);
    comm_info_1_.allGather(&loc_f, 1, loc_recv);
    float ans_f = 0.0;
    for(auto i=0; i<n_proc_1_; i++){
        ans_f += loc_recv[i];
    }
    if(rank_0_ == 3){
        std::cerr<<"glb_f= "<<glb_f<<" ans_f= "<<ans_f<<std::endl;
    }
    EXPECT_EQ(glb_f, ans_f);

    PS::F64vec loc_v = PS::F64vec((PS::F64)(rank_1_+1)*10.0,
                                  -100.0*(PS::F64)(rank_1_+1),
                                  (PS::F64)(rank_1_+1)*1000.0);
    auto glb_v = comm_info_1_.getSum(loc_v);
    PS::F64vec * loc_recv_v = new PS::F64vec[n_proc_1_];
    comm_info_1_.allGather(&loc_v, 1, loc_recv_v);
    PS::F64vec ans_v = 0.0;
    for(auto i=0; i<n_proc_1_; i++){
        ans_v += loc_recv_v[i];
    }
    if(rank_0_ == 3){
        std::cerr<<"glb_v= "<<glb_v<<" ans_v= "<<ans_v<<std::endl;
    }
    EXPECT_EQ(glb_v.x, ans_v.x);
    EXPECT_EQ(glb_v.y, ans_v.y);
    EXPECT_EQ(glb_v.z, ans_v.z);
}


int main(int argc, char *argv[]){
    std::cerr<<"start"<<std::endl;
    PS::Initialize(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    auto result = RUN_ALL_TESTS();

    PS::Finalize();
    return 0;
}
