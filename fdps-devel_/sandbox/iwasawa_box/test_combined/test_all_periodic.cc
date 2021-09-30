#include<iostream>
#include<fstream>
#include<unistd.h>
#include<particle_simulator.hpp>
#include<MT.hpp>
#include"particle_def.hpp"
//#include"particle_def_2d.hpp"

int main(int argc, char *argv[]){
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);

    PS::S32 my_rank = PS::Comm::getRank();
    
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

#if 0
    /////////////////////
    //// LONG CUTOFF ////
    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
    PS::S32 n_grav_glb, n_grav_loc;
    PS::F32 time_sys;
    ReadNemoAscii(system_grav, n_grav_glb, n_grav_loc, time_sys, sinput);


    PS::DomainInfo dinfo_grav;
    dinfo_grav.initialize();
    try{
	dinfo_grav.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XY);
	//dinfo_grav.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
    }
    catch(const char * str){
	std::cerr<<str<<std::endl;
	PS::Abort(-1);
    }
    dinfo_grav.setPosRootDomain(-1.0, 1.0); //  new 
    dinfo_grav.collectSampleParticle(system_grav);
    dinfo_grav.decomposeDomain();
    system_grav.exchangeParticle(dinfo_grav);
    n_grav_loc = system_grav.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_grav_loc; i++){
	//system_grav[i].r_cutoff = 1.0;
	system_grav[i].r_cutoff = 0.5;
	system_grav[i].eps = 1.0/128.0;
    }

    PS::TreeForForceLong<ForceGrav, EPIGravCutoff, EPJGravCutoff>::MonopoleWithCutoff tree_grav;
    PS::F32 theta_grav = 0.5;
    //PS::F32 theta_grav = 0.0;
    PS::F32 n_leaf_limit = 10;
    PS::F32 n_grp_limit = 30;
    tree_grav.initialize(n_grav_glb, theta_grav, n_leaf_limit, n_grp_limit);

    tree_grav.setParticleLocalTree(system_grav);

    tree_grav.setRootCell(dinfo_grav);
    tree_grav.mortonSortLocalTreeOnly();
    tree_grav.checkMortonSortLocalTreeOnly();

    tree_grav.linkCellLocalTreeOnly();
    tree_grav.checkMakeLocalTree();

    tree_grav.calcMomentLocalTreeOnly();
    tree_grav.checkCalcMomentLocalTree();
    
    tree_grav.exchangeLocalEssentialTree(dinfo_grav);
    tree_grav.checkExchangeLocalEssentialTree(dinfo_grav);


    tree_grav.setLocalEssentialTreeToGlobalTree();

    tree_grav.mortonSortGlobalTreeOnly();
    tree_grav.checkMortonSortGlobalTreeOnly();

    tree_grav.linkCellGlobalTreeOnly();
    tree_grav.checkMakeGlobalTree();
    
    tree_grav.calcMomentGlobalTreeOnly();
    tree_grav.checkCalcMomentGlobalTree();


    tree_grav.makeIPGroup();
    tree_grav.checkMakeIPGroup();

    PS::S32 n_ipg_grav = tree_grav.getNumberOfIPG();
    for(PS::S32 i=0; i<n_ipg_grav; i++){
        tree_grav.makeInteractionList(i);
        //tree_grav.checkMakeInteractionList(dinfo_grav);
        tree_grav.calcForceOnly
	    (CalcForceEpEpCutoff<EPIGravCutoff, EPJGravCutoff, ForceGrav>(),
	     CalcForceSpEpCutoff<EPIGravCutoff, PS::SPJMonopoleCutoff, ForceGrav>(),
	     i);
    }

    tree_grav.copyForceOriginalOrder();
    tree_grav.checkForce
	( CalcForceEpEpCutoff<EPIGravCutoff, EPJGravCutoff, ForceGrav>(),
	  CompareGrav<ForceGrav>(),
	  dinfo_grav);

#endif // LONG_FORCE
    
    //#ifdef SHORT_FORCE
#if 1
    ///////////////
    //// SHORT ////
    PS::ParticleSystem<FPSPH> system_sph;
    system_sph.initialize();
    PS::S32 n_sph_glb, n_sph_loc;
    PS::F32 time_sys;
    ReadNemoAscii(system_sph, n_sph_glb, n_sph_loc, time_sys, sinput);

    PS::DomainInfo dinfo;
    dinfo.initialize();
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain(-1.0, 1.0);
    dinfo.collectSampleParticle(system_sph);
    dinfo.decomposeDomain();
    bool pa[3];
    dinfo.getPeriodicAxis(pa);
    PS::F64 half_len_sph_glb = system_sph.getHalfLength();
    for(PS::S32 i=0; i<n_sph_loc; i++){
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
	system_sph[i].r_search = pow( (n_sph_glb/(half_len_sph_glb*half_len_sph_glb)), -0.5) * 3 * (1+(system_sph[i].id%10+1)*0.1);
#else
	system_sph[i].r_search = pow( (n_sph_glb/(half_len_sph_glb*half_len_sph_glb*half_len_sph_glb)), -0.333333) * 3 * (1+(system_sph[i].id%10+1)*0.1);
#endif
    }
    system_sph.exchangeParticle(dinfo);
    n_sph_loc = system_sph.getNumberOfParticleLocal();
    
    ///////////////
    /// SCATTER ///
#if 1
    PS::TreeForForceShort<ResultDens, EPIScatter, EPJScatter>::Scatter tree_scatter;
    tree_scatter.initialize(n_sph_glb);
#if 0
    tree_scatter.calcForceAllAndWriteBack(CalcDensityScatter(), system_sph, dinfo);
    tree_scatter.checkForce( CalcDensityScatter(), CompareDensity(), dinfo);
#else
    tree_scatter.setParticleLocalTree(system_sph);
    tree_scatter.setRootCell(dinfo);
    tree_scatter.mortonSortLocalTreeOnly();
    tree_scatter.checkMortonSortLocalTreeOnly();
    
    tree_scatter.linkCellLocalTreeOnly();
    tree_scatter.checkMakeLocalTree();

    tree_scatter.calcMomentLocalTreeOnly();
    tree_scatter.checkCalcMomentLocalTree();

    tree_scatter.exchangeLocalEssentialTree(dinfo);
    tree_scatter.checkExchangeLocalEssentialTree(dinfo);

    tree_scatter.setLocalEssentialTreeToGlobalTree();

    tree_scatter.mortonSortGlobalTreeOnly();
    tree_scatter.checkMortonSortGlobalTreeOnly();

    tree_scatter.linkCellGlobalTreeOnly();
    tree_scatter.checkMakeGlobalTree();

    tree_scatter.calcMomentGlobalTreeOnly();
    tree_scatter.checkCalcMomentGlobalTree();

    tree_scatter.makeIPGroup();
    tree_scatter.checkMakeIPGroup();

    PS::S32 n_ipg_sph_scatter = tree_scatter.getNumberOfIPG();
    for(PS::S32 i=0; i<n_ipg_sph_scatter; i++){
        tree_scatter.makeInteractionList(i);
        //tree_scatter.checkMakeInteractionList(dinfo, i);
        tree_scatter.calcForceOnly( CalcDensityScatter(), i);
    }
    tree_scatter.copyForceOriginalOrder();
    tree_scatter.checkForce( CalcDensityScatter(), CompareDensity(), dinfo);

#endif
#endif

    ///////////////
    /// GATHER ///
#if 0
    PS::TreeForForceShort<ResultDens, EPIGather, EPJGather>::Gather tree_gather;
    tree_gather.initialize(n_sph_glb);
#if 1
    tree_gather.calcForceAllAndWriteBackWithCheck(CalcDensityGather(), system_sph, dinfo);
    tree_gather.checkForce( CalcDensityGather(), CompareDensity(), dinfo);
    //tree_gather.calcForceAllAndWriteBack(CalcDensityGather(), system_sph, dinfo);
    //tree_gather.checkForce( CalcDensityGather(), CompareDensity(), dinfo);
#else
    tree_gather.setParticleLocalTree(system_sph);

    tree_gather.setRootCell(dinfo);
    tree_gather.mortonSortLocalTreeOnly();
    tree_gather.checkMortonSortLocalTreeOnly();


    tree_gather.linkCellLocalTreeOnly();
    tree_gather.checkMakeLocalTree();

    tree_gather.calcMomentLocalTreeOnly();
    tree_gather.checkCalcMomentLocalTree();

    tree_gather.exchangeLocalEssentialTree(dinfo);
    tree_gather.checkExchangeLocalEssentialTree(dinfo);

    tree_gather.setLocalEssentialTreeToGlobalTree();

    tree_gather.mortonSortGlobalTreeOnly();
    tree_gather.checkMortonSortGlobalTreeOnly();

    tree_gather.linkCellGlobalTreeOnly();
    tree_gather.checkMakeGlobalTree();

    tree_gather.calcMomentGlobalTreeOnly();
    tree_gather.checkCalcMomentGlobalTree();

    tree_gather.makeIPGroup();
    tree_gather.checkMakeIPGroup();

    PS::S32 n_ipg_sph_gather = tree_gather.getNumberOfIPG();
    for(PS::S32 i=0; i<n_ipg_sph_gather; i++){
        tree_gather.makeInteractionList(i);
        //tree_gather.checkMakeInteractionList(dinfo, i);
        tree_gather.calcForceOnly( CalcDensityGather(), i);
    }
    tree_gather.copyForceOriginalOrder();
    tree_gather.checkForce( CalcDensityGather(), CompareDensity(), dinfo);
#endif
#endif

    ////////////////
    /// SYMMETRY ///
#if 0
    PS::TreeForForceShort<ResultDens, EPISymmetry, EPJSymmetry>::Symmetry tree_symmetry;
    tree_symmetry.initialize(n_sph_glb);

#if 0
    tree_symmetry.calcForceAllAndWriteBack(CalcDensitySymmetry(), system_sph, dinfo);
    tree_symmetry.checkForce( CalcDensitySymmetry(), CompareDensity(), dinfo);
#else
    tree_symmetry.setParticleLocalTree(system_sph);
    tree_symmetry.setRootCell(dinfo);
    tree_symmetry.mortonSortLocalTreeOnly();
    tree_symmetry.checkMortonSortLocalTreeOnly();

    tree_symmetry.linkCellLocalTreeOnly();
    tree_symmetry.checkMakeLocalTree();

    tree_symmetry.calcMomentLocalTreeOnly();
    tree_symmetry.checkCalcMomentLocalTree();

    tree_symmetry.exchangeLocalEssentialTree(dinfo);
    tree_symmetry.checkExchangeLocalEssentialTree(dinfo, 1e-4);


    tree_symmetry.setLocalEssentialTreeToGlobalTree();

    tree_symmetry.mortonSortGlobalTreeOnly();
    tree_symmetry.checkMortonSortGlobalTreeOnly();

    tree_symmetry.linkCellGlobalTreeOnly();
    tree_symmetry.checkMakeGlobalTree();

    tree_symmetry.calcMomentGlobalTreeOnly();
    tree_symmetry.checkCalcMomentGlobalTree();

    tree_symmetry.makeIPGroup();
    tree_symmetry.checkMakeIPGroup();


    PS::S32 n_ipg_sph_symmetry = tree_symmetry.getNumberOfIPG();
    for(PS::S32 i=0; i<n_ipg_sph_symmetry; i++){
        tree_symmetry.makeInteractionList(i);
        //tree_symmetry.checkMakeInteractionList(i);
        tree_symmetry.calcForceOnly( CalcDensitySymmetry(), i);
    }

    tree_symmetry.copyForceOriginalOrder();
    tree_symmetry.checkForce( CalcDensitySymmetry(), CompareDensity(), dinfo);
#endif
#endif
#endif //SHORT_FORCE

    PS::Finalize();
    return 0;

}
