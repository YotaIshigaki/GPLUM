#include<iostream>
#include<fstream>
#include<unistd.h>
#include<particle_simulator.hpp>

#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
#include"particle_def_2d.hpp"
#else
#include"particle_def.hpp"
#endif

//#define OPEN_BOUNDARY
#define PERIODIC_BOUNDARY

//#define FORCE_TYPE_LONG
#define FORCE_TYPE_SHORT
//#define FORCE_TYPE_LONG_SEARCH

int main(int argc, char *argv[]){
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);

    char sinput[1024];
    int c;
    while((c=getopt(argc,argv,"i:h")) != -1){
        switch(c){
        case 'i':
            sprintf(sinput,optarg);
            break;
        case 'h':
            std::cerr<<"i: input file name (nemo ascii)"<<std::endl;
            return 0;
        }
    }

    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
    PS::S32 n_grav_glb, n_grav_loc;
    PS::F64 time_sys;
    ReadNemoAscii(system_grav, n_grav_glb, n_grav_loc, time_sys, sinput);

    PS::DomainInfo dinfo;
    dinfo.initialize();
#ifdef OPEN_BOUNDARY
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
    dinfo.collectSampleParticle(system_grav);
    dinfo.decomposeDomain();
    system_grav.exchangeParticle(dinfo);
    n_grav_loc = system_grav.getNumberOfParticleLocal();
#elif defined PERIODIC_BOUNDARY
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.collectSampleParticle(system_grav);
    dinfo.decomposeDomain();
    system_grav.exchangeParticle(dinfo);
    n_grav_loc = system_grav.getNumberOfParticleLocal();
#endif

#ifdef FORCE_TYPE_LONG
    PS::TreeForForceLong<ForceGrav, EPIGrav, EPJGrav>::Monopole tree_grav;
    PS::F32 theta_grav = 0.4;
    tree_grav.initialize(n_grav_glb, theta_grav);

#if 1
    tree_grav.calcForceAllAndWriteBackWithCheck
	(CalcForceEpEp<EPIGrav, EPJGrav, ForceGrav>(),
	 CalcForceSpEp<EPIGrav, PS::SPJMonopole, ForceGrav>(),
	 system_grav, dinfo, true);
    tree_grav.checkForce
	( CalcForceEpEp<EPIGrav, EPJGrav, ForceGrav>(),
	  CompareGrav<ForceGrav>(),
	  dinfo);
    std::cout<<"mem size used:"<<(double)tree_grav.getMemSizeUsed()*1e-9<<"[GByte]"<<std::endl;
#else
    tree_grav.calcForceAllAndWriteBack
	(CalcForceEpEp<EPIGrav, EPJGrav, ForceGrav>(),
	 CalcForceSpEp<EPIGrav, PS::SPJMonopole, ForceGrav>(),
	 system_grav, dinfo, true);
#endif

    //////////////////////////
    //// LONG CUTOFF MODE ////
    //#ifdef FORCE_TYPE_LONG_CUTOFF
    PS::TreeForForceLong<ForceGrav, EPIGravCutoff, EPJGravCutoff>::MonopoleWithCutoff tree_grav_cutoff;
    for(PS::S32 i=0; i<system_grav.getNumberOfParticleLocal(); i++){
	system_grav[i].r_cutoff = 0.5;
    }
    PS::F32 theta_grav_cutoff = 0.4;    
    tree_grav_cutoff.initialize(n_grav_glb, theta_grav_cutoff);
    
#if 1
    tree_grav_cutoff.calcForceAllAndWriteBackWithCheck
        (CalcForceEpEpCutoff<EPIGravCutoff, EPJGravCutoff, ForceGrav>(),
         CalcForceSpEpCutoff<EPIGravCutoff, PS::SPJMonopoleCutoff, ForceGrav>(),
         system_grav, dinfo, true);
    tree_grav_cutoff.checkForce
        ( CalcForceEpEpCutoff<EPIGravCutoff, EPJGravCutoff, ForceGrav>(),
          CompareGrav<ForceGrav>(),
          dinfo);
    std::cout<<"mem size used:"<<(double)tree_grav.getMemSizeUsed()*1e-9<<"[GByte]"<<std::endl;
#else
    tree_grav_cutoff.calcForceAllAndWriteBack
	(CalcForceEpEp<EPIGravCutoff, EPJGravCutoff, ForceGrav>(),
	 CalcForceSpEp<EPIGravCutoff, PS::SPJMonopoleCutoff, ForceGrav>(),
	 system_grav, dinfo, true);
#endif

#endif // FORCE_TYPE_LONG
    

///////////////////////
///////////////////////
//// LONG SHORT ///////
///////////////////////

#ifdef FORCE_TYPE_SHORT
    ///////////////
    //// SHORT ////
    PS::ParticleSystem<FPSPH> system_sph;
    system_sph.initialize();
    PS::S32 n_sph_glb, n_sph_loc;
    ReadNemoAscii(system_sph, n_sph_glb, n_sph_loc, time_sys, sinput);
    PS::F64 half_len_sph_glb = system_sph.getHalfLength();
    for(PS::S32 i=0; i<n_sph_loc; i++){
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
        system_sph[i].r_search = pow( (n_sph_glb/(half_len_sph_glb*half_len_sph_glb)), -0.5) * 3 * (1+(system_sph[i].id%10+1)*0.1);
#else //PARTICLE_SIMULATOR_TWO_DIMENSION
        //system_sph[i].r_search = pow( (n_sph_glb/(half_len_sph_glb*half_len_sph_glb*half_len_sph_glb)), -0.333333) * 3 * (1+(system_sph[i].id%10+1)*0.1);
        //system_sph[i].r_search = pow( (n_sph_glb/(1.0*1.0*1.0)), -0.333333) * 1 * (1+(system_sph[i].id%10+1)*0.1);
        system_sph[i].r_search = 0.1;
#endif // PARTICLE_SIMULATOR_TWO_DIMENSION
    }
    system_sph.exchangeParticle(dinfo);
    n_sph_loc = system_sph.getNumberOfParticleLocal();


    /////////////////////
    //// GATHER MODE ////
    PS::TreeForForceShort<ResultDens, EPIGather, EPJGather>::Gather tree_gather;
    tree_gather.initialize(n_sph_glb, 0.0, 1, 16);
#if 1
    tree_gather.calcForceAllAndWriteBackWithCheck(CalcDensityGather(), system_sph, dinfo);
    tree_gather.checkForce( CalcDensityGather(), CompareDensity(), dinfo);
#else
    tree_gather.calcForceAllAndWriteBack(CalcDensityGather(), system_sph, dinfo);
#endif

    
/*
    //////////////////////
    //// SCATTER MODE ////
    PS::TreeForForceShort<ResultDens, EPIScatter, EPJScatter>::Scatter tree_scatter;
    tree_scatter.initialize(n_sph_glb);
#if 1
    tree_scatter.calcForceAllAndWriteBackWithCheck(CalcDensityScatter(), system_sph, dinfo);
    tree_scatter.checkForce( CalcDensityScatter(), CompareDensity(), dinfo);
#else
    tree_scatter.calcForceAllAndWriteBack(CalcDensityScatter(), system_sph, dinfo);
#endif
    PS::TimeProfile tp = tree_scatter.getTimeProfile();
    std::cout<<"********************************************"<<std::endl;
    std::cout<<"********************************************"<<std::endl;
    std::cout<<"********************************************"<<std::endl;
    std::cout<<"tree_scatter.getTimeProfile().calc_force="<<tp.calc_force<<std::endl;
    std::cout<<"tree_scatter.getTimeProfile().make_LET_1st="<<tree_scatter.getTimeProfile().make_LET_1st<<std::endl;
    std::cout<<"tree_scatter.getTimeProfile().exchange_LET_1st="<<tree_scatter.getTimeProfile().exchange_LET_1st<<std::endl;
    std::cout<<"tree_scatter.getTimeProfile().make_LET_2nd="<<tree_scatter.getTimeProfile().make_LET_2nd<<std::endl;
    std::cout<<"tree_scatter.getTimeProfile().exchange_LET_2nd="<<tree_scatter.getTimeProfile().exchange_LET_2nd<<std::endl;
    std::cout<<"tree_scatter.getNumberOfInteractionEPEPLocal()="<<tree_scatter.getNumberOfInteractionEPEPLocal()<<std::endl;
    std::cout<<"tree_scatter.getNumberOfInteractionEPEPGlobal()="<<tree_scatter.getNumberOfInteractionEPEPGlobal()<<std::endl;
    std::cout<<"********************************************"<<std::endl;
    std::cout<<"********************************************"<<std::endl;
    std::cout<<"********************************************"<<std::endl;
*/
/*
    ///////////////////////
    //// SYMMETRY MODE ////
    PS::TreeForForceShort<ResultDens, EPISymmetry, EPJSymmetry>::Symmetry tree_symmetry;
    tree_symmetry.initialize(n_sph_glb);
#if 1
    tree_symmetry.calcForceAllAndWriteBackWithCheck(CalcDensitySymmetry(), system_sph, dinfo);
    tree_symmetry.checkForce( CalcDensitySymmetry(), CompareDensity(), dinfo);
#else
    tree_symmetry.calcForceAllAndWriteBack(CalcDensitySymmetry(), system_sph, dinfo);
#endif
*/
#endif //FORCE_TYPE_SHORT


////////////////////////
//// LONG SEARCH ///////
////////////////////////

#ifdef FORCE_TYPE_LONG_SEARCH
    ////////////////////////////////
    //// FORCE_TYPE_LONG_SEARCH ////
    PS::ParticleSystem<FPLongScatter> system_p3t;
    system_p3t.initialize();
    PS::S32 n_p3t_glb, n_p3t_loc;
    ReadNemoAscii(system_p3t, n_p3t_glb, n_p3t_loc, time_sys, sinput);
    PS::F64 half_len_p3t_glb = system_p3t.getHalfLength();
    for(PS::S32 i=0; i<n_p3t_loc; i++){
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
	system_p3t[i].r_search = pow( (n_p3t_glb/(half_len_p3t_glb*half_len_p3t_glb)), -0.5) * 3 * (1+(system_p3t[i].id%10+1)*0.1);
#else //PARTICLE_SIMULATOR_TWO_DIMENSION
	//system_p3t[i].r_search = pow( (n_p3t_glb/(half_len_p3t_glb*half_len_p3t_glb*half_len_p3t_glb)), -0.333333) * 3 * (1+(system_p3t[i].id%10+1)*0.1);
	system_p3t[i].r_search = pow( (n_p3t_glb/(1.0*1.0*1.0)), -0.333333) * 1 * (1+(system_p3t[i].id%10+1)*0.1);
	//system_p3t[i].r_search = 0.1;
	system_p3t[i].r_search_ep = system_p3t[i].r_search * 0.5;
#endif // PARTICLE_SIMULATOR_TWO_DIMENSION
    }
    system_p3t.exchangeParticle(dinfo);
    n_p3t_loc = system_p3t.getNumberOfParticleLocal();

#if 0
    ////////////////////////
    //// CUTOFF SCATTER ////
    PS::TreeForForceLong<ForceLongScatter, EPILongScatter, EPJLongScatter>::MonopoleWithCutoffScatterSearch tree_long_cutoff_scatter;
    tree_long_cutoff_scatter.initialize(n_p3t_glb, 0.4);
    tree_long_cutoff_scatter.setParticleLocalTree(system_p3t);
    tree_long_cutoff_scatter.setRootCell(dinfo);
    tree_long_cutoff_scatter.mortonSortLocalTreeOnly();
    tree_long_cutoff_scatter.linkCellLocalTreeOnly();
    tree_long_cutoff_scatter.calcMomentLocalTreeOnly();
    tree_long_cutoff_scatter.exchangeLocalEssentialTree(dinfo);


#else
    // scatter open
    /////////////////
    //// SCATTER ////
    std::cerr<<"scatter test"<<std::endl;
    PS::TreeForForceLong<ForceLongScatter, EPILongScatter, EPJLongScatter>::MonopoleWithScatterSearch tree_long_scatter;
    tree_long_scatter.initialize(n_p3t_glb, 0.4);
    tree_long_scatter.setParticleLocalTree(system_p3t);
    tree_long_scatter.setRootCell(dinfo);
    tree_long_scatter.mortonSortLocalTreeOnly();
    tree_long_scatter.linkCellLocalTreeOnly();
    tree_long_scatter.calcMomentLocalTreeOnly();
    tree_long_scatter.exchangeLocalEssentialTree(dinfo);
    tree_long_scatter.setLocalEssentialTreeToGlobalTree();
    tree_long_scatter.mortonSortGlobalTreeOnly();
    tree_long_scatter.linkCellGlobalTreeOnly();
    tree_long_scatter.calcMomentGlobalTreeOnly();
    tree_long_scatter.makeIPGroup();
    tree_long_scatter.calcForceAndWriteBack
        (CalcForceSearchEPEP<EPILongScatter, EPJLongScatter, ForceLongScatter>(), 
         CalcForceSearchEPSP<EPILongScatter, PS::SPJMonopoleScatter, ForceLongScatter>(), 
         system_p3t, true);
    tree_long_scatter.checkForce( CalcForceSearchEPEP<EPILongScatter, EPJLongScatter, 
                                  ForceLongScatter>(),
                                  CompareNNGB<ForceLongScatter>(), dinfo);
    //std::cout<<"tree_long_scatter.getNumberOfInteractionEPEPGlobal()="<<tree_long_scatter.getNumberOfInteractionEPEPGlobal()<<std::endl;
    //std::cout<<"tree_long_scatter.getNumberOfInteractionEPSPGlobal()="<<tree_long_scatter.getNumberOfInteractionEPSPGlobal()<<std::endl;
    FPLongScatter fp_ls;
    bool err_ls = false;
    for(PS::S32 i=0; i<system_p3t.getNumberOfParticleLocal(); i++){
        fp_ls.pos = system_p3t[i].getPos();
        fp_ls.r_search = system_p3t[i].getRSearch();
        
        EPJLongScatter * epj_ls = NULL;
        PS::S32 n_neighbor_ls = tree_long_scatter.getNeighborListOneParticle(fp_ls, epj_ls);
        if(system_p3t[i].n_ngb != n_neighbor_ls-1){
            std::cout<<"n_neighbor_ls="<<n_neighbor_ls<<" system_p3t[i].n_ngb="<<system_p3t[i].n_ngb<<std::endl;
            err_ls = true;
        }
        for(PS::S32 j=0; j<n_neighbor_ls; j++){
            if( epj_ls[j].getPos().getDistanceSQ(fp_ls.pos) > epj_ls[j].getRSearch()*epj_ls[j].getRSearch()){
                std::cout<<"epj_ls[j].getPos()="<<epj_ls[j].getPos()<<std::endl;
                std::cout<<"fp_ls.pos="<<fp_ls.pos<<std::endl;
                std::cout<<"epj_ls[j].getRSearch()="<<epj_ls[j].getRSearch()<<std::endl;
                err_ls = true;
            }
        }

    }
    if(err_ls) std::cout<<"FAIL: getNeighborListOneParticle"<<std::endl;
    else std::cout<<"PASS: getNeighborListOneParticle"<<std::endl;

#endif // scatter open

/*
#if 1
    tree_scatter.checkForce( CalcDensityScatter(), CompareDensity(), dinfo);
#else
    tree_scatter.calcForceAllAndWriteBack(CalcDensityScatter(), system_sph, dinfo);
#endif
*/

/*
    PS::TimeProfile tp = tree_scatter.getTimeProfile();
    std::cout<<"********************************************"<<std::endl;
    std::cout<<"********************************************"<<std::endl;
    std::cout<<"********************************************"<<std::endl;
    std::cout<<"tree_scatter.getTimeProfile().calc_force="<<tp.calc_force<<std::endl;
    std::cout<<"tree_scatter.getTimeProfile().make_LET_1st="<<tree_scatter.getTimeProfile().make_LET_1st<<std::endl;
    std::cout<<"tree_scatter.getTimeProfile().exchange_LET_1st="<<tree_scatter.getTimeProfile().exchange_LET_1st<<std::endl;
    std::cout<<"tree_scatter.getTimeProfile().make_LET_2nd="<<tree_scatter.getTimeProfile().make_LET_2nd<<std::endl;
    std::cout<<"tree_scatter.getTimeProfile().exchange_LET_2nd="<<tree_scatter.getTimeProfile().exchange_LET_2nd<<std::endl;
    std::cout<<"tree_scatter.getNumberOfInteractionEPEPLocal()="<<tree_scatter.getNumberOfInteractionEPEPLocal()<<std::endl;
    std::cout<<"tree_scatter.getNumberOfInteractionEPEPGlobal()="<<tree_scatter.getNumberOfInteractionEPEPGlobal()<<std::endl;
    std::cout<<"********************************************"<<std::endl;
    std::cout<<"********************************************"<<std::endl;
    std::cout<<"********************************************"<<std::endl;
*/

#endif //FORCE_TYPE_LONG_SEARCH





    PS::Finalize();
    return 0;
}
