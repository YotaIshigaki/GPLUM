#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <map>
#include <dirent.h>
#include <ctime>
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include <mpi.h>
#endif

#define USE_RE_SEARCH_NEIGHBOR
#define PIKG_USE_FDPS_VECTOR

#include <particle_simulator.hpp>

#define GPLUMVERSION "2.4 (2021/09)"

#define PRC(x) std::cerr << #x << " = " << x << ", "
#define PRL(x) std::cerr << #x << " = " << x << "\n"

constexpr PS::F64 MY_PI = 3.14159265358979323846L;

#include "mathfunc.h"
#include "kepler.h"
#include "particle.h"

using EPI_t    = EPIGrav;
using EPJ_t    = EPJGrav;
#ifdef USE_QUAD
using SPJ_t    = MySPJQuadrupole;
using Moment_t = MyMomentQuadrupole;
#else
using SPJ_t    = MySPJMonopole;
using Moment_t = MyMomentMonopole;
#endif
using Force_t  = ForceGrav;
using FP_t  = FPGrav;
using FPH_t = FPHard;

#include "energy.h"
#include "neighbor.h"
#include "disk.h"
#include "gravity_hard.h"
#include "gravity_soft.h"
#include "collisionA.h"
#include "collisionB.h"
#include "hermite.h"
#include "hard.h"
#include "read.h"
#include "time.h"
#include "func.h"
#include "gravity_kernel_epep.hpp"
#include "gravity_kernel_epsp.hpp"

#ifdef USE_INDIVIDUAL_CUTOFF
using MY_SEARCH_MODE = PS::SEARCH_MODE_LONG_SYMMETRY;
#else
using MY_SEARCH_MODE = PS::SEARCH_MODE_LONG_SCATTER;
#endif
#ifdef USE_POLAR_COORDINATE
using Tree_t = PS::TreeForForce<MY_SEARCH_MODE, Force_t, EPI_t, EPJ_t, Moment_t, Moment_t, SPJ_t, PS::CALC_DISTANCE_TYPE_NEAREST_X>;
#else
using Tree_t = PS::TreeForForce<MY_SEARCH_MODE, Force_t, EPI_t, EPJ_t, Moment_t, Moment_t, SPJ_t, PS::CALC_DISTANCE_TYPE_NORMAL>;
#endif


int main(int argc, char *argv[])
{
    PS::Initialize(argc, argv); 
    showGplumVersion(GPLUMVERSION);
    //PS::Comm::barrier();
    time_t wtime_start_program = time(NULL);
    
    ////////////////////////
    /*   Set Parameters   */
    ////////////////////////
    // Set Default Parameters
    char param_file[256] = "parameter.dat";
    
    char init_file[256] = "INIT_3000.dat";
    char output_dir[256] = "OUTPUT";
    bool bHeader  = false;
    bool bRestart = false;

    bool makeInit = false;

    PS::F64 coef_ema = 0.3;
    PS::S32 nx = 0, ny = 0;
    
    PS::F64 theta         = 0.5;
    PS::S32 n_leaf_limit  = 8;
    PS::S32 n_group_limit = 256;
    PS::S32 n_smp_ave     = 100;
    
    PS::F64 t_end       = 1.;
    PS::F64 dt_snap     = pow2(-5);
    PS::F64 dt_snap_tmp = pow2(-5);

    PS::F64 r_max = 40.;
    PS::F64 r_min = 0.1;

    PS::S32 seed = 1;
    PS::F64 wtime_max = 0.;

    PS::S32 reset_step = 1024;

    //FP_t::setGamma(FP_t::gamma);
    
    // Read Parameter File
    opterr = 0;
    PS::S32 opt;
    char init_file_opt[256];
    char output_dir_opt[256];
    PS::S32 seed_opt;
    PS::F64 dt_opt;
    PS::F64 Rcut0_opt, Rcut1_opt;
    PS::S32 nx_opt, ny_opt;
    bool opt_r=false, opt_i=false, opt_s=false, opt_o=false,
        opt_D=false, opt_R=false, opt_S=false, opt_x=false, opt_y=false;
    while ((opt = getopt(argc, argv, "p:ri:s:e:o:D:R:S:x:y:")) != -1) {
        switch (opt) {
        case 'p':
            sprintf(param_file,"%s",optarg);
            break;
        case 'r':
            opt_r = true;
            break;
        case 'i':
            sprintf(init_file_opt,"%s",optarg);
            opt_i = true;
            break;
        case 's':
            seed_opt = std::atoi(optarg);
            opt_s = true;
            break;
        case 'e':
            wtime_max = std::atof(optarg)*60.*60.;
            break;
        case 'o':
            sprintf(output_dir_opt,"%s",optarg);
            opt_o = true;
            break;
        case 'D':
            dt_opt = pow2(-std::atoi(optarg));
            opt_D = true;
            break;
        case 'R':
            Rcut0_opt = std::atof(optarg);
            opt_R = true;
            break;
        case 'S':
            Rcut1_opt = std::atof(optarg);
            opt_S = true;
            break;
        case 'x':
            nx_opt = std::atoi(optarg);
            opt_x = true;
            break;
        case 'y':
            ny_opt = std::atoi(optarg);
            opt_y = true;
            break;
        default:
            std::cout << "Usage: "
                      <<  argv[0]
                      << " [-p argment] [-r] [-i argment] [-s argment] [-e argment] [-o argment] [-D argment] [-R argment] [-S argment] [-x argment] [-y argment] "
                      << std::endl;
            break;
        }
    }

    if ( readParameter(param_file, init_file, bHeader, output_dir, bRestart,  makeInit,
                       coef_ema, nx, ny,
                       theta, n_leaf_limit, n_group_limit, n_smp_ave,
                       t_end, dt_snap, dt_snap_tmp, r_max, r_min, seed, reset_step) ){
        PS::Comm::barrier();
        PS::Abort();
        PS::Comm::barrier();
        PS::Finalize();
        return 0;
        
    }
    
    if (opt_r) bRestart = true;
    if (opt_i) sprintf(init_file,"%s",init_file_opt);
    if (opt_s) seed = seed_opt;
    if (opt_o) sprintf(output_dir,"%s",output_dir_opt);
    if (opt_D) FP_t::dt_tree = dt_opt;
    if (opt_R) FP_t::R_cut0 = Rcut0_opt;
    if (opt_S) FP_t::R_cut1 = Rcut1_opt;
    if (opt_x) nx = nx_opt;
    if (opt_y) ny = ny_opt;

    if ( nx == 0 && ny == 0 ) {
        nx = (int)sqrt(PS::Comm::getNumberOfProc());
        while ( PS::Comm::getNumberOfProc()%nx != 0 ) nx++;
        ny = PS::Comm::getNumberOfProc()/nx;
    } else if ( nx == 0 ) {
        nx = PS::Comm::getNumberOfProc()/ny;
    } else if ( ny == 0 ) {
        ny = PS::Comm::getNumberOfProc()/nx;
    }

    if ( checkParameter(init_file, bHeader, output_dir, bRestart,  makeInit,
                        coef_ema, nx, ny,
                        theta, n_leaf_limit, n_group_limit, n_smp_ave,
                        t_end, dt_snap, dt_snap_tmp, r_max, r_min, seed, reset_step) ){
        PS::Comm::barrier();
        PS::Abort();
        PS::Comm::barrier();
        PS::Finalize();
        return 0;
    }

    FP_t::setGamma(FP_t::gamma);
    
    char dir_name[256];
    //sprintf(dir_name,"./%s",output_dir);
    sprintf(dir_name,"%s",output_dir);
    if ( makeOutputDirectory(dir_name) ){
        PS::Comm::barrier();
        PS::Abort();
        PS::Comm::barrier();
        PS::Finalize();
        return 0;
    }
    if ( bRestart ) sprintf(init_file, "%s/snap_tmp.dat", output_dir);
    
    srand48( seed + PS::Comm::getRank() );

    PS::F64 time_sys = 0.;
    PS::S32 n_tot  = 0;
    PS::F64 de_max = 0.;
    PS::F64 de_d_cum = 0.;
    PS::S32 istep     = 0;
    PS::S32 isnap     = 0;
    PS::S32 n_col_tot  = 0;
    PS::S32 n_frag_tot = 0;
    
    ////////////////////////////////
    /*   Create Particle System   */
    ////////////////////////////////
    // Soft System
    PS::ParticleSystem<FP_t> system_grav;
    system_grav.initialize();
    system_grav.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave);
    PS::S32 n_loc = 0;
    Energy e_init, e_now;
#ifdef OUTPUT_DETAIL
    PS::F64 ekin_before = 0., ekin_after = 0., edisp_d = 0.;
#endif
    PS::S64 id_next = 0;
    
    NeighborList NList;
    //ExPair::initialize();

    if ( !bRestart ) {
        if ( makeInit ){
            //Make Initial Condition
            SolidDisk::createInitialCondition(system_grav);
            istep = 0;
            isnap = 0;
            time_sys = 0.;
            sprintf(init_file, "NONE");
            bHeader = false;
        } else {
            // Read Initial File
            if ( bHeader ) {
                FileHeader header;
                PS::F64 dt_tree = FP_t::dt_tree;
                system_grav.readParticleAscii(init_file, header);
                if ( PS::Comm::getRank() == 0 ){
                    istep = (PS::S32)round(header.time/dt_tree);
                    isnap = (PS::S32)round(header.time/dt_snap);
                    e_init = header.e_init;
                    e_now = header.e_now;
                }
                PS::Comm::barrier();
                PS::Comm::broadcast(&istep, 1);
                PS::Comm::broadcast(&isnap, 1);
                PS::Comm::broadcast(&e_init, 1);
                time_sys = istep*dt_tree;
                id_next = header.id_next-1;
            } else {
                system_grav.readParticleAscii(init_file);
                istep = 0;
                isnap = 0;
                time_sys = 0.;
            }
        }
    } else {
        FileHeader header;
        PS::F64 dt_tree = FP_t::dt_tree;
        system_grav.readParticleBinary(init_file, header);
        if ( PS::Comm::getRank() == 0 ){
            istep = (PS::S32)round(header.time/dt_tree);
            isnap = (PS::S32)round(header.time/dt_snap);
            e_init = header.e_init;
            e_now = header.e_now;
        }
        PS::Comm::barrier();
        PS::Comm::broadcast(&istep, 1);
        PS::Comm::broadcast(&isnap, 1);
        PS::Comm::broadcast(&e_init, 1);
        time_sys = istep*dt_tree;
        id_next = header.id_next-1;
    }
    n_loc = system_grav.getNumberOfParticleLocal();
    n_tot = system_grav.getNumberOfParticleGlobal();
#pragma omp parallel for reduction (max: id_next)
    for ( PS::S32 i=0; i<n_loc; i++ ){
        //#pragma omp critical
        //{
        if ( system_grav[i].id > id_next ) id_next = system_grav[i].id;
        //}
        system_grav[i].time = time_sys;
        system_grav[i].neighbor = system_grav[i].n_cluster = 0;
        system_grav[i].isMerged = false;
        system_grav[i].isDead = false;
        if ( system_grav[i].r_planet <= 0. ) system_grav[i].setRPlanet();
        if ( system_grav[i].f == 0. ) system_grav[i].f = FP_t::increase_factor;
    }
    id_next = PS::Comm::getMaxValue(id_next);
    id_next ++;
    setCutoffRadii(system_grav);
#ifdef USE_POLAR_COORDINATE
    setPosPolar(system_grav);
#endif
        
    // Hard System
    HardSystem system_hard;
    ExParticleSystem<FP_t> system_ex;
    system_hard.clear();
    system_ex.initialize();
    PS::S32 n_ngb_tot        = 0;
    PS::S32 n_with_ngb       = 0;
    PS::S32 n_largestcluster = 0;
    PS::S32 n_cluster        = 0;
    PS::S32 n_isoparticle    = 0;
    

    ////////////////////
    /*   Set Domain   */
    ////////////////////
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.setNumberOfDomainMultiDimension(nx,ny,1);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
#ifdef USE_POLAR_COORDINATE
    dinfo.setPosRootDomainX(-MY_PI, MY_PI);
#endif
    dinfo.collectSampleParticle(system_grav, true);
    dinfo.decomposeDomain();
    
    system_grav.exchangeParticle(dinfo);
    n_loc = system_grav.getNumberOfParticleLocal();
    
    setIDLocalAndMyrank(system_grav, NList);

    /////////////////////
    /*   Create Tree   */
    /////////////////////
    Tree_t tree_grav;
    tree_grav.initialize(n_tot, theta, n_leaf_limit, n_group_limit);
#ifdef USE_P2P_FAST
    tree_grav.setExchangeLETMode(PS::EXCHANGE_LET_P2P_FAST);
#endif
    
#ifdef USE_RE_SEARCH_NEIGHBOR
    for (PS::S32 i=0; i<2; i++) {
#endif
        tree_grav.calcForceAllAndWriteBack(
#ifdef USE_INDIVIDUAL_CUTOFF
                                           CalcForceLongEPEP(FP_t::eps2),
#else
                                           CalcForceLongEPEP(FP_t::eps2, FP_t::r_out, FP_t::r_search),
#endif
                                           CalcForceLongEPSP(FP_t::eps2),
                                           system_grav,
                                           dinfo,
                                           true, PS::MAKE_LIST_FOR_REUSE, false);
        //NList.initializeList(system_grav);
        if ( bRestart ) {
            correctForceLong(system_grav, tree_grav, NList, n_ngb_tot, n_with_ngb);
        } else {
            correctForceLongInitial(system_grav, tree_grav, NList, n_ngb_tot, n_with_ngb);
        }
#ifdef USE_RE_SEARCH_NEIGHBOR
    }
#endif

#ifdef INDIRECT_TERM
    calcIndirectTerm(system_grav);
#endif

#ifdef GAS_DRAG
    GasDisk gas_disk;
    gas_disk.calcGasDrag(system_grav, time_sys);
#endif
    
    //////////////////
    /*  File Open   */
    //////////////////
    std::ofstream fout_eng;
    std::ofstream fout_col;
    std::ofstream fout_rem;
    char sout_eng[256];
    char sout_col[256];
    char sout_rem[256];

    if ( PS::Comm::getRank() == 0 ) {
        sprintf(sout_eng, "%s/energy.dat",    dir_name);
        //sprintf(sout_col, "%s/collision.dat", dir_name);
        //sprintf(sout_rem, "%s/remove.dat",    dir_name);
        sprintf(sout_col, "%s/collision%06d.dat", dir_name, isnap+1);
        sprintf(sout_rem, "%s/remove%06d.dat",    dir_name, isnap+1);
        if ( time_sys == 0. ) {
            fout_eng.open(sout_eng, std::ios::out);
            //fout_col.open(sout_col, std::ios::out);
            //fout_rem.open(sout_rem, std::ios::out);
        } else {
            fout_eng.open(sout_eng, std::ios::app);
            //fout_col.open(sout_col, std::ios::app);
            //fout_rem.open(sout_rem, std::ios::app);
        }
        fout_col.open(sout_col, std::ios::out);
        fout_rem.open(sout_rem, std::ios::out);
    }
    //PS::Comm::barrier();
    
    showParameter(init_file, dir_name, makeInit,
                  time_sys,
                  coef_ema, nx, ny,
                  theta, n_leaf_limit, n_group_limit, n_smp_ave,
                  t_end, dt_snap, dt_snap_tmp, r_max, r_min, seed, reset_step);

    ////////////////////////////////
    /*  Preparation Before Loop   */
    ////////////////////////////////
    e_now.calcEnergy(system_grav);
    if ( !bHeader ) e_init = e_now;
    PS::F64 de =  e_now.calcEnergyError(e_init);
    
    Wtime wtime;

    PS::Comm::barrier();
    wtime.init = wtime.now = PS::GetWtime();

    outputStep(system_grav, time_sys, e_init, e_now, de,
               n_col_tot, n_frag_tot, dir_name, isnap, id_next, fout_eng, 
               wtime, n_largestcluster, n_cluster, n_isoparticle, 
               (time_sys==0.) );
    istep ++;
    isnap ++;
    
    
    //////////////////////
    ///   Loop Start
    while(1){
        
        PS::S32 n_col = 0;
        PS::S32 n_frag = 0;
        PS::S32 n_remove = 0;
#ifdef GAS_DRAG
        PS::F64 edisp_gd = 0.;
#endif
        
        n_loc = system_grav.getNumberOfParticleLocal();

        PS::Comm::barrier();
        wtime.start_soft = PS::GetWtime();
#ifdef CALC_WTIME
        wtime.lap(wtime.start_soft);
#endif

        ////////////////////
        ///   Soft Part
        
        ///////////////////////////
        /*   1st Velocity kick   */
        ///////////////////////////
#ifdef GAS_DRAG
        correctEnergyForGas(system_grav, edisp_gd, false);
#endif
        velKick(system_grav);
#ifdef OUTPUT_DETAIL
        calcKineticEnergy(system_grav, ekin_before);
#endif
        
        ///   Soft Part
        ////////////////////
        

        PS::Comm::barrier();
        wtime.end_soft = wtime.start_hard = PS::GetWtime();
        wtime.soft += wtime.soft_step = wtime.end_soft - wtime.start_soft;
#ifdef CALC_WTIME
        wtime.calc_soft_force += wtime.lap(wtime.end_soft);
#endif
        
     
        /////////////////////
        ///   Hard Part

        ////////////////////////////
        /*   Create Hard System   */
        ////////////////////////////
        NList.createConnectedRankList();
        //NList.makeIdMap(system_grav);
        
        //PS::S32 ** ex_data_send = new PS::S32*[n_send];
        //PS::S32 ** ex_data_recv = new PS::S32*[n_send];
        //for ( PS::S32 ii=0; ii<n_send; ii++ ) {
        //    PS::S32 n_size = NList.getNumberOfPairConnected(ii) * ExPair::getSize();            
        //    ex_data_send[ii] = new PS::S32[n_size];
        //    ex_data_recv[ii] = new PS::S32[n_size];
        //}
        NList.resizeExDataBuffer();

        bool check = true;
        bool check_loc = true;
        PS::S32 TAG = 10;
        while ( check ) {
            NList.createNeighborCluster(system_grav);
            NList.inputExData(system_grav);
            check_loc = NList.exchangeExData(system_grav, TAG);
            check = PS::Comm::synchronizeConditionalBranchOR(check_loc);
            //PS::Comm::barrier();
            TAG ++ ;
        }
        
        // for ( PS::S32 ii=0; ii<n_send; ii++ ) {          
        //     delete [] ex_data_send[ii];
        //     delete [] ex_data_recv[ii];
        // }
        // delete [] ex_data_send;
        // delete [] ex_data_recv;

        NList.selectSendRecvParticle(system_grav);
        
#ifdef CALC_WTIME
        PS::Comm::barrier();
        wtime.create_cluster += wtime.create_cluster_step = wtime.lap(PS::GetWtime());
#endif

        system_ex.resize(NList.getNumberOfRankSend(), NList.getNumberOfRankRecv());

            
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        /////////////////////
        /*   Gather Hard   */
        /////////////////////
        // Send & Receive Particles Number
        system_ex.inputNumberOfExParticleSend(system_grav, NList);
        system_ex.sendRecvNumberOfExParticle(NList);

        // Send & Receive Particles
        system_ex.inputAdress();
        system_ex.inputExParticleSend(system_grav, NList);
        system_ex.sendRecvExParticle(NList);        
        system_ex.inputNeighborListOfExParticleRecv();
#endif
        
#ifdef CALC_WTIME
        PS::Comm::barrier();
        wtime.communication_step = wtime.lap(PS::GetWtime());
#endif
    
        ////////////////////////
        /*   Time Integrate   */
        ////////////////////////
        PS::S32 n_in = 0;
        PS::S32 n_out = 0;
        system_hard.clear();
        system_hard.reserve_hard(2 * (n_with_ngb + system_ex.getNumberOfParticleRecv()));
        n_in = system_hard.makeList(system_grav, system_ex);
        n_out = system_hard.timeIntegrate(system_grav, system_ex,
                                          NList.n_list, system_ex.n_list, istep);
        //system_hard.showParticleID();
        assert( n_in == n_out ); 
        assert( n_loc == system_hard.getNumberOfParticleLocal()-system_hard.getNumberOfFragmentLocal()
                +system_ex.getNumberOfParticleSend()-system_ex.getNumberOfParticleRecv() );
            
        //PS::Comm::barrier();
        n_col  = system_hard.getNumberOfCollisionGlobal();
        if ( n_col ) n_frag = system_hard.addFragment2ParticleSystem(system_grav, system_ex,
                                                                     id_next, fout_col);
        
        n_col_tot   += n_col;
        n_frag_tot  += n_frag;
        e_now.edisp += system_hard.getEnergyDissipationGlobal();
#ifdef OUTPUT_DETAIL
        edisp_d      = system_hard.getHardEnergyDissipationGlobal();
#endif
        if( time_sys+FP_t::dt_tree  == dt_snap*isnap ){
            n_largestcluster = system_hard.getNumberOfParticleInLargestClusterGlobal();
            n_cluster        = system_hard.getNumberOfClusterGlobal();
            n_isoparticle    = system_hard.getNumberOfIsolatedParticleGlobal();
        }

#ifdef CALC_WTIME
        PS::Comm::barrier();
        wtime.calc_hard_force += wtime.calc_hard_force_step = wtime.lap(PS::GetWtime());
#endif
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        //////////////////////
        /*   Scatter Hard   */
        //////////////////////
        system_ex.returnExParticle(NList);
        system_ex.outputExParticleSend(system_grav, NList);       
#endif
        ///   Hard Part
        ////////////////////              
   
        PS::Comm::barrier();
        wtime.end_hard = wtime.start_soft = PS::GetWtime();
        wtime.hard += wtime.hard_step = wtime.end_hard - wtime.start_hard;
#ifdef CALC_WTIME
        wtime.communication_step += wtime.lap(wtime.end_hard);
        wtime.communication += wtime.communication_step;
#endif


        ////////////////////
        ///   Soft Part

        ////////////////////////
        /*   Calculate Soft   */
        ////////////////////////
#ifdef USE_POLAR_COORDINATE
        setPosPolar(system_grav);
#endif
        system_grav.exchangeParticle(dinfo);
        setIDLocalAndMyrank(system_grav, NList);
        tree_grav.calcForceAllAndWriteBack(
#ifdef USE_INDIVIDUAL_CUTOFF
                                           CalcForceLongEPEP(FP_t::eps2),
#else
                                           CalcForceLongEPEP(FP_t::eps2, FP_t::r_out, FP_t::r_search),
#endif
                                           CalcForceLongEPSP(FP_t::eps2),
                                           system_grav,
                                           dinfo,
                                           true, PS::MAKE_LIST_FOR_REUSE, false);
#ifdef CALC_WTIME
        PS::Comm::barrier();
        wtime.calc_soft_force += wtime.calc_soft_force_step = wtime.lap(PS::GetWtime());
#endif
        //NList.initializeList(system_grav);
        correctForceLong(system_grav, tree_grav, NList, n_ngb_tot, n_with_ngb);
#ifdef INDIRECT_TERM
        calcIndirectTerm(system_grav);
#endif
#ifdef CALC_WTIME
        PS::Comm::barrier();
        wtime.neighbor_search += wtime.neighbor_search_step = wtime.lap(PS::GetWtime());
#endif
        
#ifdef GAS_DRAG
        gas_disk.calcGasDrag(system_grav, time_sys+FP_t::dt_tree);
#endif

        ///////////////////////////
        /*   2nd Velocity kick   */
        ///////////////////////////
#ifdef OUTPUT_DETAIL
        calcKineticEnergy(system_grav, ekin_after);
#endif
#ifndef CORRECT_NEIGHBOR
        velKick(system_grav);
#else
        velKick2nd(system_grav);
#endif
#ifdef GAS_DRAG
        correctEnergyForGas(system_grav, edisp_gd, true);
        e_now.edisp += edisp_gd;
#endif
        
        //////////////////////////
        /*   Calculate Energy   */
        //////////////////////////
#ifdef OUTPUT_DETAIL
        Energy e_tmp = e_now;
        PS::F64 dekin_d = ekin_after - ekin_before;
#endif
        e_now.calcEnergy(system_grav);
#ifdef OUTPUT_DETAIL
        PS::F64 dephi_d_d = e_now.ephi_d - e_tmp.ephi_d;
        PS::F64 dephi_s_d = e_now.ephi_sun - e_tmp.ephi_sun;
        PS::F64 de_d = ( dekin_d + dephi_d_d + dephi_s_d - edisp_d ) / e_init.etot;
        de_d_cum += de_d;
#endif

        ///////////////
        /*   Merge   */
        ///////////////
        if ( n_col ) {
            MergeParticle(system_grav, n_col, e_now.edisp);
        }

        // Remove Particle Out Of Boundary
        n_remove = removeParticlesOutOfBoundary(system_grav, e_now.edisp, r_max, r_min, fout_rem);

        ///////////////////////////
        /*   Re-Calculate Soft   */
        ///////////////////////////
        if ( n_col || n_remove || istep % reset_step == reset_step-1 ) {
            if( istep % reset_step == reset_step-1 ) {
#ifdef USE_POLAR_COORDINATE
                setPosPolar(system_grav);
#endif
                dinfo.decomposeDomainAll(system_grav);
                system_grav.exchangeParticle(dinfo);
 
                // Remove Particle Out Of Boundary
                //removeParticlesOutOfBoundary(system_grav, e_now.edisp, r_max, r_min, fout_rem);
            }
                
            // Reset Number Of Particles
            n_tot = system_grav.getNumberOfParticleGlobal();
            n_loc = system_grav.getNumberOfParticleLocal();
            
#ifdef CALC_WTIME
            PS::Comm::barrier();
            wtime.lap(PS::GetWtime());
#endif
            setIDLocalAndMyrank(system_grav, NList);
            setCutoffRadii(system_grav);
            tree_grav.calcForceAllAndWriteBack(
#ifdef USE_INDIVIDUAL_CUTOFF
                                               CalcForceLongEPEP(FP_t::eps2),
#else
                                               CalcForceLongEPEP(FP_t::eps2, FP_t::r_out, FP_t::r_search),
#endif
                                               CalcForceLongEPSP(FP_t::eps2),
                                               system_grav,
                                               dinfo,
                                               true, PS::MAKE_LIST_FOR_REUSE, false);
#ifdef CALC_WTIME
            PS::Comm::barrier();
            PS::F64 time_tmp = wtime.lap(PS::GetWtime());
            wtime.calc_soft_force_step += time_tmp;
            wtime.calc_soft_force += time_tmp;
#endif
            //NList.initializeList(system_grav);
            correctForceLongInitial(system_grav, tree_grav, NList, n_ngb_tot, n_with_ngb);
#ifdef INDIRECT_TERM
    calcIndirectTerm(system_grav);
#endif
#ifdef CALC_WTIME
            PS::Comm::barrier();
            time_tmp = wtime.lap(PS::GetWtime());
            wtime.neighbor_search_step += time_tmp;
            wtime.neighbor_search += time_tmp;
#endif
#ifdef GAS_DRAG
#pragma omp parallel for
            for(PS::S32 i=0; i<n_loc; i++) system_grav[i].acc += system_grav[i].acc_gd;
#endif
            e_now.calcEnergy(system_grav);
        }
   
        ///   Soft Part
        ////////////////////

        //for(PS::S32 i=0; i<n_loc; i++) system_grav[i].dump();

        PS::Comm::barrier();
        wtime.now = wtime.end_soft = PS::GetWtime();
        wtime.soft += wtime.end_soft - wtime.start_soft;
        wtime.soft_step +=  wtime.end_soft - wtime.start_soft;
#ifdef CALC_WTIME
        wtime.calc_soft_force += wtime.calc_soft_force_step = wtime.lap(wtime.now);
#endif
        
        ////////////////
        /*   Output   */
        ////////////////
        time_sys += FP_t::dt_tree;
        istep ++;

        PS::F64 de =  e_now.calcEnergyError(e_init);
        PS::F64 de_tmp = sqrt(de*de);
        if( de_tmp > de_max ) de_max = de_tmp;
        if ( PS::Comm::getRank() == 0 ) {
            std::cout << std::fixed << std::setprecision(8)
                      << "Time: " << time_sys
                      << "  NumberOfParticle: " << n_tot 
                      << "  NumberOfCol: " << n_col_tot
                      << "  NumberOfFrag: " << n_frag_tot << std::endl
                      << std::scientific << std::setprecision(15)
                      << "                  "
#ifdef OUTPUT_DETAIL
                      << "HardEnergyError: " << de_d_cum
                      << "  "
#endif
                      << "EnergyError: " << de
                      << "  MaxEnergyError: " << de_max << std::endl;
            
            //std::cerr << std::scientific<<std::setprecision(15);
            //PRC(etot1); PRL(ekin);
            //PRC(ephi); PRC(ephi_s); PRL(ephi_d);
        }
                
        if( time_sys  == dt_snap*isnap ){
            outputStep(system_grav, time_sys, e_init, e_now, de,
                       n_col_tot, n_frag_tot, dir_name, isnap, id_next, fout_eng,
                       wtime, n_largestcluster, n_cluster, n_isoparticle);
            makeSnapTmp(system_grav, time_sys, e_init, e_now, dir_name, id_next);
            isnap ++;

            if ( time_sys >= t_end || 
                 (wtime_max > 0. && wtime_max < difftime(time(NULL), wtime_start_program)) ) break;

            if ( PS::Comm::getRank() == 0 ) {
                fout_col.close();
                fout_rem.close();
                sprintf(sout_col, "%s/collision%06d.dat", dir_name, isnap);
                sprintf(sout_rem, "%s/remove%06d.dat",    dir_name, isnap);
                fout_col.open(sout_col, std::ios::out);
                fout_rem.open(sout_rem, std::ios::out);
            }
        }

        if( fmod(time_sys, dt_snap_tmp) == 0. ){
            makeSnapTmp(system_grav, time_sys, e_init, e_now, dir_name, id_next);
            
            if ( wtime_max > 0. && wtime_max < difftime(time(NULL), wtime_start_program) ) break;
        }
#ifdef CALC_WTIME
        PS::Comm::barrier();
        wtime.output += wtime.output_step = wtime.lap(PS::GetWtime());
#endif
    }
       
    ///   Loop End
    ////////////////////

    if ( PS::Comm::getRank() == 0 ) {
        fout_eng.close();
        fout_col.close();
        fout_rem.close();
    }
    
    PS::Comm::barrier();

    wtime.now = PS::GetWtime();
    wtime.showTime(dir_name, wtime_start_program, de_max, dinfo, system_grav, tree_grav);
    
    PS::Finalize();
    return 0;
}
