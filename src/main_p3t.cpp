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

#include <particle_simulator.hpp>

#define PRC(x) std::cerr << #x << " = " << x << ", "
#define PRL(x) std::cerr << #x << " = " << x << "\n"

#include "mathfunc.h"
#include "kepler.h"
#include "energy.h"
#include "particle.h"
#include "neighbor.h"
#include "disk.h"
#include "gravity.h"
#include "gravity_kernel.hpp"
#include "collisionA.h"
#include "collisionB.h"
#include "hermite.h"
#include "hard.h"
#include "read.h"
#include "time.h"
#include "func.h"

PS::F64 EPGrav::eps2   = 0.;
#ifndef USE_INDIVIDUAL_RADII
PS::F64 EPGrav::r_out;
PS::F64 EPGrav::r_out_inv;
PS::F64 EPGrav::r_search;
#endif
PS::F64 EPGrav::R_cut     = 1.;
PS::F64 EPGrav::R_search0 = 1.;
PS::F64 EPGrav::R_search1 = 1.;
PS::F64 EPGrav::gamma     = 0.5;
PS::F64 EPGrav::g2        = 0.25;   // gamma^2
PS::F64 EPGrav::g_1_inv   = -2.;    // 1/(gamma-1)
PS::F64 EPGrav::g2_1_inv  = -4./3.; // 1/(gamma^2-1)
PS::F64 EPGrav::g_1_inv7  = -128.; // 1/(gamma-1)^7
PS::F64 EPGrav::w_y;      // (g^6-9g^5+45*g^4-60g^3*log(g)-45g^2+9g-1)/3(g-1)^7

PS::F64 FPGrav::m_sun     = 1.;
PS::F64 FPGrav::dens      = 5.049667e6;
PS::F64 FPGrav::dt_tree   = pow2(-5);
PS::F64 FPGrav::dt_min    = pow2(-13);
PS::F64 FPGrav::eta       = 0.01;
PS::F64 FPGrav::eta_0     = 0.001;
PS::F64 FPGrav::alpha2    = 0.01;
PS::F64 FPGrav::rHill_min = 0;
PS::F64 FPGrav::rHill_max = 0;

PS::F64 HardSystem::f = 1.;
PS::F64 Collision0::f = 1.;
PS::F64 Collision0::m_min = 9.426627927538057e-12;


int main(int argc, char *argv[])
{
   PS::Initialize(argc, argv);
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
    PS::S32 nx = (int)sqrt(PS::Comm::getNumberOfProc());
    while ( PS::Comm::getNumberOfProc()%nx != 0 ) nx++;
    PS::S32 ny = PS::Comm::getNumberOfProc()/nx;
    
    PS::F64 theta         = 0.5;
    PS::S32 n_leaf_limit  = 8;
    PS::S32 n_group_limit = 256;
    PS::S32 n_smp_ave     = 100;
    
    PS::F64 t_end   = pow2(-2);
    PS::F64 dt_snap = pow2(-5);

    PS::F64 r_max = 40.;
    PS::F64 r_min = 0.1;

    PS::S32 seed = 1;
    PS::F64 wtime_max = 0.;

    EPGrav::setGamma(EPGrav::gamma);
    
    // Read Parameter File
    opterr = 0;
    PS::S32 opt;
    char init_file_opt[256];
    char output_dir_opt[256];
    PS::S32 seed_opt;
    PS::F64 dt_opt;
    PS::F64 Rcut_opt;
    bool opt_r=false, opt_i=false, opt_s=false, opt_o=false, opt_D=false, opt_R=false;
    while ((opt = getopt(argc, argv, "p:ri:s:e:o:D:R:")) != -1) {
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
            Rcut_opt = std::atof(optarg);
            opt_R = true;
            break;
        default:
            std::cout << "Usage: "
                      <<  argv[0]
                      << " [-p argment] [-r] [-i argment] [-s argment] [-e argment] [-o argment] [-D argment] [-R argment] arg1 ..."
                      << std::endl;
            break;
        }
    }
    if ( readParameter(param_file, init_file, bHeader, output_dir, bRestart,  makeInit,
                       coef_ema, nx, ny,
                       theta, n_leaf_limit, n_group_limit, n_smp_ave,
                       t_end, dt_snap, r_max, r_min, seed) ){
        std::cerr << "Parameter file has NOT successfully read" << std::endl;
        PS::Abort();
        return 0;
    }
    if (opt_r) bRestart = true;
    if (opt_i) sprintf(init_file,"%s",init_file_opt);
    if (opt_s) seed = seed_opt;
    if (opt_o) sprintf(output_dir,"%s",output_dir_opt);
    if (opt_D) FPGrav::dt_tree = dt_opt;
    if (opt_R) EPGrav::R_cut = EPGrav::R_search0 = Rcut_opt;
    
    char dir_name[256];
    sprintf(dir_name,"./%s",output_dir);
    if ( makeOutputDirectory(dir_name) ){
        PS::Abort();
        return 0;
    }
    if ( bRestart ) {
        if ( getLastSnap(dir_name, init_file) ) {
            PS::Abort();
            return 0;
        }
        bHeader = true;
        makeInit = false;
    }

    srand48( seed + PS::Comm::getRank() );

    PS::F64 time_sys = 0.;
    PS::S32 n_tot  = 0;
    PS::F64 de_max = 0.;
    PS::S32 istep = 0;
    PS::S32 isnap = 0;
    PS::S32 n_col_tot  = 0;
    PS::S32 n_frag_tot = 0;
    
    ////////////////////////////////
    /*   Create Particle System   */
    ////////////////////////////////
    // Soft System
    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
    system_grav.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave);
    PS::S32 n_loc = 0;
    Energy e_init, e_now;
#ifdef OUTPUT_DETAIL
    PS::F64 ekin_before = 0., ekin_after = 0., edisp_d = 0.;
#endif
    PS::S32 id_next = 0;
    
    NeighborList NList;
    NList.initialize();
    ExPair::initialize();

    if ( makeInit && !bRestart ){
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
            PS::F64 dt_tree = FPGrav::dt_tree;
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
        } else {
            system_grav.readParticleAscii(init_file);
            istep = 0;
            isnap = 0;
            time_sys = 0.;
        }
    }
    n_loc = system_grav.getNumberOfParticleLocal();
    n_tot = system_grav.getNumberOfParticleGlobal();
#pragma omp parallel for
    for ( PS::S32 i=0; i<n_loc; i++ ){
#pragma omp critical
        {
            if ( system_grav[i].id > id_next ) id_next = system_grav[i].id;
        }
        system_grav[i].time = time_sys;
        system_grav[i].neighbor = system_grav[i].n_cluster = 0;
    }
    id_next = PS::Comm::getMaxValue(id_next);
    id_next ++;
    setCutoffRadii(system_grav);
        
    // Hard System
    HardSystem system_hard;
    ExParticleSystem<FPGrav> system_ex;
    system_hard.clear();
    system_ex.initialize();
    PS::S32 nei_dist         = 0;
    PS::S32 nei_tot_loc      = 0;
    PS::S32 n_largestcluster = 0;
    PS::S32 n_cluster        = 0;
    PS::S32 n_isoparticle    = 0;

    ////////////////////
    /*   Set Domain   */
    ////////////////////
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.setNumberOfDomainMultiDimension(nx,ny,1);
    dinfo.collectSampleParticle(system_grav, true);
    dinfo.decomposeDomain();
    system_grav.exchangeParticle(dinfo);
    n_loc = system_grav.getNumberOfParticleLocal();
    
    inputIDLocalAndMyrank(system_grav);
    
    /////////////////////
    /*   Create Tree   */
    /////////////////////
#ifdef USE_INDIVIDUAL_RADII
#ifdef USE_QUAD
    PS::TreeForForceLong<ForceGrav, EPGrav, EPGrav>::QuadrupoleWithSymmetrySearch tree_grav;
#else
    PS::TreeForForceLong<ForceGrav, EPGrav, EPGrav>::MonopoleWithSymmetrySearch tree_grav;
#endif
#else
#ifdef USE_QUAD
    PS::TreeForForceLong<ForceGrav, EPGrav, EPGrav>::QuadrupoleWithScatterSearch tree_grav;
#else
    PS::TreeForForceLong<ForceGrav, EPGrav, EPGrav>::MonopoleWithScatterSearch tree_grav;
#endif
#endif
    tree_grav.initialize(n_tot, theta, n_leaf_limit, n_group_limit);

    tree_grav.calcForceAllAndWriteBack(CalcForceLongEP<EPGrav,EPGrav,ForceGrav>,
#ifdef USE_INDIVIDUAL_RADII
#ifdef USE_QUAD
                                       CalcForceLongSP<EPGrav,PS::SPJQuadrupoleInAndOut,ForceGrav>,
#else
                                       CalcForceLongSP<EPGrav,PS::SPJMonopoleInAndOut,ForceGrav>,
#endif //USE_QUAD
#else
#ifdef USE_QUAD
                                       CalcForceLongSP<EPGrav,PS::SPJQuadrupoleScatter,ForceGrav>,
#else
                                       CalcForceLongSP<EPGrav,PS::SPJMonopoleScatter,ForceGrav>,
#endif //USE_QUAD
#endif //USE_INDIVIDUAL_RADII
                                       system_grav,
                                       dinfo);
    NList.initializeList(system_grav);
    correctForceLongInitial(system_grav, tree_grav, NList, nei_dist, nei_tot_loc);

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
    sprintf(sout_eng, "%s/energy.dat",    dir_name);
    sprintf(sout_col, "%s/collision.dat", dir_name);
    sprintf(sout_rem, "%s/remove.dat",    dir_name);
    if ( time_sys == 0. ) {
        fout_eng.open(sout_eng, std::ios::out);
        fout_col.open(sout_col, std::ios::out);
        fout_rem.open(sout_rem, std::ios::out);
    } else {
        fout_eng.open(sout_eng, std::ios::app);
        fout_col.open(sout_col, std::ios::app);
        fout_rem.open(sout_rem, std::ios::app);
    }
    PS::Comm::barrier();
    showParameter(init_file, dir_name, makeInit,
                  time_sys,
                  coef_ema, nx, ny,
                  theta, n_leaf_limit, n_group_limit, n_smp_ave,
                  t_end, dt_snap, r_max, r_min, seed);

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
               n_col_tot, n_frag_tot, dir_name, isnap, fout_eng, 
               wtime, n_largestcluster, n_cluster, n_isoparticle, 
               (time_sys==0.) );
    istep ++;
    isnap ++;
    
    //////////////////////
    ///   Loop Start
    while(time_sys < t_end){
        
        PS::S32 n_col = 0;
        PS::S32 n_frag = 0;
        std::map<PS::S32, PS::S32> id2id_loc;
        id2id_loc.clear();
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
        NList.makeIdMap(system_grav);

        PS::S32 n_send = NList.getNumberOfRankConnected();
        PS::S32 ** ex_data_send = new PS::S32*[n_send];
        PS::S32 ** ex_data_recv = new PS::S32*[n_send];
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {
            PS::S32 n_size = NList.getNumberOfPairConnected(ii) * ExPair::getSize();            
            ex_data_send[ii] = new PS::S32[n_size];
            ex_data_recv[ii] = new PS::S32[n_size];
        }

        bool check = true;
        bool check_loc = true;
        PS::S32 TAG = 10;
        while ( check ) {
            NList.createNeighborCluster(system_grav);
            NList.inputExData(system_grav);
            check_loc = NList.exchangeExData(system_grav, TAG, ex_data_send, ex_data_recv);
            check = PS::Comm::synchronizeConditionalBranchOR(check_loc);
            //PS::Comm::barrier();
            TAG ++ ;
        }
        
        for ( PS::S32 ii=0; ii<n_send; ii++ ) {          
            delete [] ex_data_send[ii];
            delete [] ex_data_recv[ii];
        }
        delete [] ex_data_send;
        delete [] ex_data_recv;

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
        if( time_sys+FPGrav::dt_tree  == dt_snap*isnap ){
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
        system_grav.exchangeParticle(dinfo);
        inputIDLocalAndMyrank(system_grav);
        tree_grav.calcForceAllAndWriteBack(CalcForceLongEP<EPGrav,EPGrav,ForceGrav>,
#ifdef USE_INDIVIDUAL_RADII
#ifdef USE_QUAD
                                           CalcForceLongSP<EPGrav,PS::SPJQuadrupoleInAndOut,ForceGrav>,
#else
                                           CalcForceLongSP<EPGrav,PS::SPJMonopoleInAndOut,ForceGrav>,
#endif //USE_QUAD
#else
#ifdef USE_QUAD
                                           CalcForceLongSP<EPGrav,PS::SPJQuadrupoleScatter,ForceGrav>,
#else
                                           CalcForceLongSP<EPGrav,PS::SPJMonopoleScatter,ForceGrav>,
#endif //USE_QUAD
#endif //USE_INDIVIDUAL_RADII
                                           system_grav,
                                           dinfo);
#ifdef CALC_WTIME
        PS::Comm::barrier();
        wtime.calc_soft_force += wtime.calc_soft_force_step = wtime.lap(PS::GetWtime());
#endif
        NList.initializeList(system_grav);
        correctForceLong(system_grav, tree_grav, NList, nei_dist, nei_tot_loc);
#ifdef CALC_WTIME
        PS::Comm::barrier();
        wtime.neighbor_search += wtime.neighbor_search_step = wtime.lap(PS::GetWtime());
#endif
        
#ifdef GAS_DRAG
        gas_disk.calcGasDrag(system_grav, time_sys+FPGrav::dt_tree);
#endif

        ///////////////////////////
        /*   2nd Velocity kick   */
        ///////////////////////////
#ifdef OUTPUT_DETAIL
        calcKineticEnergy(system_grav, ekin_after);
#endif
        velKick(system_grav);
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
#endif

        ///////////////
        /*   Merge   */
        ///////////////
        if ( n_col ) {
            MergeParticle(system_grav, n_col, e_now.edisp);
        }
        
        ///////////////////////////
        /*   Re-Calculate Soft   */
        ///////////////////////////
        if ( n_col || istep % 1024 == 0 ) {
            if(istep % 1024 == 0) {
                dinfo.decomposeDomainAll(system_grav);
                system_grav.exchangeParticle(dinfo);
                
                // Remove Particle Out Of Boundary
                removeOutOfBoundaryParticle(system_grav, e_now.edisp, r_max, r_min, fout_rem);
            }
                
            // Reset Number Of Particles
            n_tot = system_grav.getNumberOfParticleGlobal();
            n_loc = system_grav.getNumberOfParticleLocal();
            
#ifdef CALC_WTIME
            PS::Comm::barrier();
            wtime.lap(PS::GetWtime());
#endif
            inputIDLocalAndMyrank(system_grav);
            setCutoffRadii(system_grav);
            tree_grav.calcForceAllAndWriteBack(CalcForceLongEP<EPGrav,EPGrav, ForceGrav>,
#ifdef USE_INDIVIDUAL_RADII
#ifdef USE_QUAD
                                               CalcForceLongSP<EPGrav,PS::SPJQuadrupoleInAndOut,ForceGrav>,
#else
                                               CalcForceLongSP<EPGrav,PS::SPJMonopoleInAndOut,ForceGrav>,
#endif //USE_QUAD
#else
#ifdef USE_QUAD
                                               CalcForceLongSP<EPGrav,PS::SPJQuadrupoleScatter,ForceGrav>,
#else
                                               CalcForceLongSP<EPGrav,PS::SPJMonopoleScatter,ForceGrav>,
#endif //USE_QUAD
#endif //USE_INDIVIDUAL_RADII
                                               system_grav,
                                               dinfo);
#ifdef CALC_WTIME
            PS::Comm::barrier();
            PS::F64 time_tmp = wtime.lap(PS::GetWtime());
            wtime.calc_soft_force_step += time_tmp;
            wtime.calc_soft_force += time_tmp;
#endif
            NList.initializeList(system_grav);
            correctForceLongInitial(system_grav, tree_grav, NList, nei_dist, nei_tot_loc);
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
        time_sys += FPGrav::dt_tree;
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
                      << "HardEnergyError: " << de_d
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
                       n_col_tot, n_frag_tot, dir_name, isnap, fout_eng,
                       wtime, n_largestcluster, n_cluster, n_isoparticle);
            isnap ++;

            if ( wtime_max > 0. && wtime_max < difftime(time(NULL), wtime_start_program) ) break; 
        }
#ifdef CALC_WTIME
        PS::Comm::barrier();
        wtime.output += wtime.output_step = wtime.lap(PS::GetWtime());
#endif
    }
       
    ///   Loop End
    ////////////////////
    
    fout_eng.close();
    fout_col.close();
    fout_rem.close();
    
    PS::Comm::barrier();

    wtime.now = PS::GetWtime();
    wtime.showTime(dir_name, wtime_start_program, de_max, dinfo, system_grav, tree_grav);
    
    PS::Finalize();
    return 0;
}
