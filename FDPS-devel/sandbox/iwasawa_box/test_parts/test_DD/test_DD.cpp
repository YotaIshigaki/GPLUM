#include<iostream>
#include<unistd.h>
#include<particle_simulator.hpp>
#include"gtest/gtest.h"

int N_TOT;
int N_SMP;

class FPGrav{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pot;
    PS::F64vec getPos() const { return pos; }
};

class DDTest : public ::testing::Test{
protected:
    PS::S32 rank_0_;
    PS::S32 rank_1_;
    PS::S32 n_proc_0_;
    PS::S32 n_proc_1_;
    PS::CommInfo comm_0_;
    PS::CommInfo comm_1_;
    PS::F64 coef_ema_;
    PS::ParticleSystem<FPGrav> system_grav_;
    PS::DomainInfo dinfo_;
    virtual void SetUp(){
        comm_0_ = PS::Comm::getCommInfo();
        auto my_rank = comm_0_.getRank();
        comm_1_ = PS::Comm::getCommInfo().split(my_rank%2, my_rank);
        rank_0_ = comm_0_.getRank();
        rank_1_ = comm_1_.getRank();
        n_proc_0_ = comm_0_.getNumberOfProc();
        n_proc_1_ = comm_1_.getNumberOfProc();
        coef_ema_ = 0.3;
        std::cerr<<"N_TOT= "<<N_TOT
                 <<" N_SMP= "<<N_SMP
                 <<std::endl;
    }
};

TEST_F(DDTest, basic){

}

/*
TEST(test_wheretogo_case, test_wheretogo){
    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    const PS::S32 n_smp = 100;
    const PS::S32 n_proc = 10;
    const PS::F64 dx = 1.0;
    const PS::F64 len_peri = dx*n_proc;
    const PS::S32 my_rank = 0;
    PS::F64 * coord_sep = new PS::F64[n_proc+1];
    coord_sep[0] = 0.0;
    for(auto i=0; i<n_proc; i++){
        coord_sep[i+1] = coord_sep[i] + dx;
    }
    PS::MTTS mt;
    mt.init_genrand(0);
    for(auto i=0; i<n_smp; i++){
        const auto coord_smp = ((double)(n_proc) * mt.genrand_res53() * 0.99999);
        const auto rank = dinfo.whereToGoOpen(coord_smp, coord_sep, my_rank);
        EXPECT_LE(coord_sep[rank], coord_smp);
        EXPECT_LT(coord_smp, coord_sep[rank+1]);
    }

    for(auto i=0; i<n_smp; i++){
        const auto coord_smp = ((double)(n_proc) * mt.genrand_res53() * 0.99999);
        const auto rank = dinfo.whereToGoPeriodic(coord_smp, coord_sep, my_rank, len_peri, n_proc);
        EXPECT_LE(coord_sep[rank], coord_smp);
        EXPECT_LT(coord_smp, coord_sep[rank+1]);
    }    
}
*/
/*
TEST(test_select_case, test_select){
    const PS::S32 n_smp = 1000;
    const PS::S32 n_proc = 100;
    const PS::F64 dx = 1.0;
    const PS::F64 len_peri = dx*n_proc;
    PS::F64 * coord_smp = new PS::F64[n_smp];
    PS::S32 * rank_send = new PS::S32[n_smp];
    PS::F64 * coord_sep = new PS::F64[n_proc+1];
    coord_sep[0] = 0.0;
    for(auto i=0; i<n_proc; i++){
        coord_sep[i+1] = coord_sep[i] + dx;
    }
    PS::MTTS mt;
    mt.init_genrand(0);
    for(auto i=0; i<n_smp; i++){
        coord_smp[i] = ((double)(n_proc) * mt.genrand_res53() * 0.99999);
    }
    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    for(auto my_rank=0; my_rank<n_proc; my_rank++){
        dinfo.selectExchangeSmpOMP(n_smp, coord_smp, rank_send, my_rank, coord_sep, n_proc, len_peri, false);
        for(auto i=0; i<n_smp; i++){
            EXPECT_LE(coord_sep[rank_send[i]], coord_smp[i]);
            EXPECT_LT(coord_smp[i], coord_sep[rank_send[i]+1]);
        }
        dinfo.selectExchangeSmpOMP(n_smp, coord_smp, rank_send, my_rank, coord_sep, n_proc, len_peri, true);
        for(auto i=0; i<n_smp; i++){
            EXPECT_LE(coord_sep[rank_send[i]], coord_smp[i]);
            EXPECT_LT(coord_smp[i], coord_sep[rank_send[i]+1]);
        }
    }
}
*/
/*
TEST(test_count_case, test_count){
    const PS::S32 n_bucket = 100;
    PS::S32 * n_obj_in_bucket = new PS::S32[n_bucket];
    const PS::S32 n_obj = 1000;
    PS::S32 * id_bucket = new PS::S32[n_obj];
    for(auto i=0; i<n_obj; i++){
        id_bucket[i] = i % n_bucket;
    }
    PS::CountNObjInBucketOMP(n_bucket, n_obj_in_bucket, n_obj, id_bucket);
    for(auto i=0; i<n_bucket; i++){
        EXPECT_LE(n_obj_in_bucket[i], 10);
    }
    delete [] n_obj_in_bucket;
    delete [] id_bucket;
}
*/

/*
TEST(test_dd_open_case, test_open_dd){
    const auto my_rank = PS::Comm::getRank();
    const auto n_proc  = PS::Comm::getNumberOfProc();
    PS::MTTS mt;
    mt.init_genrand(my_rank);
    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
    system_grav.setAverageTargetNumberOfSampleParticlePerProcess(N_SMP);
    PS::S32 n_loc = N_TOT / n_proc * (mt.genrand_int31()%3);
    //PS::S32 n_loc = N_TOT;
    //if(my_rank != 0) n_loc = 0;
    std::cerr<<"my_rank= "<<my_rank
             <<" n_loc= "<<n_loc
             <<std::endl;
    system_grav.setNumberOfParticleLocal(n_loc);

    for(auto i=0; i<n_loc; i++){
        system_grav[i].pos.x = mt.genrand_res53();
        system_grav[i].pos.y = mt.genrand_res53();
        system_grav[i].pos.z = mt.genrand_res53();
    }

    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.collectSampleParticle(system_grav);
    dinfo.decomposeDomainMultiStep2(false);
}
*/

/*
TEST(test_dd_prriodic_case, test_periodic_dd){
    const auto my_rank = PS::Comm::getRank();
    const auto n_proc  = PS::Comm::getNumberOfProc();
    PS::MTTS mt;
    mt.init_genrand(my_rank);
    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
    system_grav.setAverageTargetNumberOfSampleParticlePerProcess(N_SMP);
    PS::S32 n_loc = N_TOT / n_proc * (mt.genrand_int31()%3);
    system_grav.setNumberOfParticleLocal(n_loc);


    for(auto i=0; i<n_loc; i++){
        system_grav[i].pos.x = mt.genrand_res53();
        system_grav[i].pos.y = 2.0*mt.genrand_res53();
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        system_grav[i].pos.z = 3.0*mt.genrand_res53();
#endif
    }

    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XY);
#else
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
#endif

#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0), 
                           PS::F64vec(1.0, 2.0));
#else
    dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), 
                           PS::F64vec(1.0, 2.0, 3.0));
#endif
    dinfo.collectSampleParticle(system_grav);
    dinfo.decomposeDomainMultiStep2(false);

    if(my_rank==0){
        for(auto i=0; i<n_proc; i++){
            std::cerr<<"dinfo.getPosDomain(i)= "<<dinfo.getPosDomain(i)<<std::endl;
        }
    }

    system_grav.exchangeParticle2(dinfo);
    n_loc = system_grav.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_loc; i++){
        EXPECT_LE(dinfo.getPosDomain(my_rank).low_.x, system_grav[i].pos.x);
        EXPECT_LT(system_grav[i].pos.x, dinfo.getPosDomain(my_rank).high_.x);
        EXPECT_LE(dinfo.getPosDomain(my_rank).low_.y, system_grav[i].pos.y);
        EXPECT_LT(system_grav[i].pos.y, dinfo.getPosDomain(my_rank).high_.y);
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        EXPECT_LE(dinfo.getPosDomain(my_rank).low_.z, system_grav[i].pos.z);
        EXPECT_LT(system_grav[i].pos.z, dinfo.getPosDomain(my_rank).high_.z);
#endif
    }
    
}
*/

/*
TEST(test_exptcl_case, test_exptcl){
    const auto my_rank = PS::Comm::getRank();
    const auto n_proc  = PS::Comm::getNumberOfProc();
    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
    system_grav.setAverageTargetNumberOfSampleParticlePerProcess(N_SMP);
    PS::S32 n_loc = N_TOT / n_proc;
    system_grav.setNumberOfParticleLocal(n_loc);
    PS::MTTS mt;
    mt.init_genrand(my_rank);
    for(auto i=0; i<n_loc; i++){
        system_grav[i].pos.x = mt.genrand_res53();
        system_grav[i].pos.y = mt.genrand_res53();
        system_grav[i].pos.z = mt.genrand_res53();
    }

    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), 
                                 PS::F64vec(1.0, 1.0, 1.0));
    dinfo.decomposeDomainAll(system_grav);
    if(my_rank == 0){
        for(auto i=0; i<PS::Comm::getNumberOfProc(); i++){
            std::cerr<<"i= "<<i<<" pos= "<<dinfo.getPosDomain(i)<<std::endl;
        }
    }
    system_grav.exchangeParticle2(dinfo);
    n_loc = system_grav.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_loc; i++){
        EXPECT_LE(dinfo.getPosDomain(my_rank).low_.x, system_grav[i].pos.x);
        EXPECT_LT(system_grav[i].pos.x, dinfo.getPosDomain(my_rank).high_.x);
        EXPECT_LE(dinfo.getPosDomain(my_rank).low_.y, system_grav[i].pos.y);
        EXPECT_LT(system_grav[i].pos.y, dinfo.getPosDomain(my_rank).high_.y);
        EXPECT_LE(dinfo.getPosDomain(my_rank).low_.z, system_grav[i].pos.z);
        EXPECT_LT(system_grav[i].pos.z, dinfo.getPosDomain(my_rank).high_.z);
    }
}
*/


int main(int argc, char * argv[]){
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);
    
    ::testing::InitGoogleTest(&argc, argv);

    PS::S64 n_tot = 1024;
    PS::S32 n_smp = 30;
    PS::S32 c;
    while((c=getopt(argc,argv,"n:N:")) != -1){
        switch(c){
        case 'n':
            n_smp = atoi(optarg);
            if(PS::Comm::getRank() == 0) {
                std::cerr << "n_smp = " << n_smp << std::endl;
            }
            break;
        case 'N':
            n_tot = atoi(optarg);
            if(PS::Comm::getRank() == 0) {
                std::cerr << "n_tot = " << n_tot << std::endl;
            }
            break;
        default:
            if(PS::Comm::getRank() == 0) {
                std::cerr<<"No such option! Available options are here."<<std::endl;
            }
            PS::Abort();
        }
    }

    N_TOT = n_tot;
    N_SMP = n_smp;

    auto ret = RUN_ALL_TESTS();

    PS::Finalize();
    
    return ret;
    
}
