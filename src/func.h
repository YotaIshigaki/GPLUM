#pragma once

void showGplumVersion(std::string version)
{
    if ( PS::Comm::getRank() == 0 ) {
        version.resize(16,' ');
        
        std::cout << "                                           \n "
                  << "     __________________________________   \n"
                  << "     /                                  \\  \n"
                  << "    |    ____ ____  _    _   _ __  __    | \n"
                  << "    |   / ___|  _ \\| |  | | | |  \\/  |   | \n"
                  << "    |  | |  _| |_) | |  | | | | |\\/| |   | \n"
                  << "    |  | |_| |  __/| |__| |_| | |  | |   | \n"
                  << "    |   \\____|_|   |_____\\___/|_|  |_|   | \n"
                  << "    |                                    | \n"
                  << "    |  Global Planetary Simulation Code  | \n"
                  << "    |    with Mass-dependent Cut-off     | \n"
                  << "    |       Version " << version << "     | \n"
                  << "     \\__________________________________/  \n"
                  << "                                           \n"
                  << "     Licence: MIT (see, https://github.com/YotaIshigaki/GPLUM/blob/master/LICENSE) \n"
                  << "                                           \n"
                  << "     Copyright (C) 2020                    \n"
                  << "       Yota Ishigaki, Junko Kominmi, Junichiro Makino, \n"
                  << "       Masaki Fujimoto and Masaki Iwasawa              \n"
                  << "                                                       \n";
    }
}

template <class Tpsys>
void calcMeanMass(Tpsys & pp,
                  PS::F64 & m_mean,
                  PS::F64 & m_max,
                  PS::F64 & nei_mean)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    const PS::S32 n_glb = pp.getNumberOfParticleGlobal();
    PS::F64 m_sum_loc = 0.;
    PS::F64 m_max_loc = 0.;
    PS::S32 nei_sum_loc = 0;
    
    for (PS::S32 i=0; i<n_loc; i++ ){
        m_sum_loc += pp[i].mass;
        if ( pp[i].mass > m_max_loc ) m_max_loc = pp[i].mass;
        nei_sum_loc += pp[i].neighbor;
    }
    m_mean = PS::Comm::getSum(m_sum_loc) / n_glb;
    m_max = PS::Comm::getMaxValue(m_max_loc);
    nei_mean = (PS::F64)PS::Comm::getSum(nei_sum_loc) / n_glb;
}

template <class Tpsys>
void makeSnap(Tpsys & pp,
              PS::F64 time_sys,
              Energy e_init,
              Energy e_now,
              const char * dir_name,
              const PS::S32 isnap,
              const PS::S64 id_next)
{
    FileHeader header(pp.getNumberOfParticleGlobal(), id_next, time_sys, e_init, e_now);
    char filename[256];
    sprintf(filename, "%s/snap%06d.dat", dir_name, isnap);
    pp.writeParticleAscii(filename, header);
}

template <class Tpsys>
void makeSnapTmp(Tpsys & pp,
                  PS::F64 time_sys,
                  Energy e_init,
                  Energy e_now,
                  const char * dir_name,
                  const PS::S64 id_next)
{
    FileHeader header(pp.getNumberOfParticleGlobal(), id_next, time_sys, e_init, e_now);
    char filename[256];
    sprintf(filename, "%s/snap_tmp.dat", dir_name);
    pp.writeParticleBinary(filename, header);
}

template <class Tpsys>
void outputStep(Tpsys & pp,
                PS::F64 time_sys,
                Energy e_init,
                Energy e_now,
                PS::F64 de,
                PS::S32 n_col_tot,
                PS::S32 n_frag_tot,
                const char * dir_name,
                const PS::S32 isnap,
                const PS::S64 id_next,
                std::ofstream & fout_eng,
                Wtime wtime,
                PS::S32 n_largestcluster,
                PS::S32 n_cluster,
                PS::S32 n_isoparticle,
                bool bSnap=true)
{
    const PS::S32 n_tot = pp.getNumberOfParticleGlobal();

    if ( bSnap ) makeSnap(pp, time_sys, e_init, e_now, dir_name, isnap, id_next);

#ifdef OUTPUT_DETAIL
    PS::F64 m_mean = 0.;
    PS::F64 m_max = 0.;
    PS::F64 nei_mean = 0.;
    calcMeanMass(pp, m_mean, m_max, nei_mean);
#endif

    if(PS::Comm::getRank() == 0 && bSnap){
        //PS::F64 de =  e_now.calcEnergyError(e_init);
        //PS::F64 de_tmp = sqrt(de*de);
        //if( de_tmp > de_max ) de_max = de_tmp;
        fout_eng  << std::fixed<<std::setprecision(8)
                  << time_sys << "\t" << n_tot << "\t"
                  << std::scientific<<std::setprecision(15)
                  << e_now.etot << "\t" << de << "\t"
                  << n_largestcluster << "\t" << n_cluster << "\t" << n_isoparticle
#ifdef OUTPUT_DETAIL
                  << "\t" << m_max << "\t" << m_mean << "\t" << nei_mean
#endif
#ifdef CALC_WTIME
                  << "\t" << wtime.soft_step << "\t" << wtime.hard_step << "\t"
                  << wtime.calc_soft_force_step << "\t" << wtime.neighbor_search_step << "\t"
                  << wtime.calc_hard_force_step << "\t"
                  << wtime.create_cluster_step << "\t" << wtime.communication_step << "\t"
                  << wtime.output_step 
#endif
                  << std::endl;
    }
}

template <class Tpsys>
void setIDLocalAndMyrank(Tpsys & pp,
                         NeighborList & NList)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    PS::S32 myrank = PS::Comm::getRank();
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        pp[i].id_local = i;
        pp[i].myrank = myrank;
        pp[i].inDomain = true;
        pp[i].isSent = false;
    }

    NList.makeIdMap(pp);
}

#ifdef USE_POLAR_COORDINATE
template <class Tpsys>
void setPosPolar(Tpsys & pp)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_loc; i++) pp[i].setPosPolar();
}
#endif

template <class Tpsys>
void MergeParticle(Tpsys & pp,
                   PS::S32 n_col,
                   PS::F64 & edisp)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    PS::S32 n_remove = 0;
    PS::S32 * remove = new PS::S32[n_col];
    PS::F64 edisp_loc = 0.;

#pragma omp parallel for reduction (-:edisp_loc)
    for ( PS::S32 i=0; i<n_loc; i++ ){
        if ( pp[i].isMerged ) {
            for ( PS::S32 j=0; j<n_loc; j++ ){              
                if ( pp[j].id == pp[i].id && i != j ){
                    PS::F64 mi = pp[i].mass;
                    PS::F64 mj = pp[j].mass;
                    PS::F64vec vrel = pp[j].vel - pp[i].vel;
                    pp[i].mass += mj;
                    pp[i].vel = ( mi*pp[i].vel + mj*pp[j].vel )/(mi+mj);
                    //pp[i].acc = ( mi*pp[i].acc + mj*pp[j].acc )/(mi+mj);
#ifdef GAS_DRAG
                    pp[i].acc_gd = ( mi*pp[i].acc_gd + mj*pp[j].acc_gd )/(mi+mj);
#endif
                    pp[i].phi   = ( mi*pp[i].phi   + mj*pp[j].phi   )/(mi+mj);
                    pp[i].phi_d = ( mi*pp[i].phi_d + mj*pp[j].phi_d )/(mi+mj);
                    
                    edisp_loc -= 0.5 * mi*mj/(mi+mj) * vrel*vrel;
#pragma omp critical
                    {
                        remove[n_remove] = j;
                        n_remove ++;
                    }
                    assert ( pp[i].pos == pp[j].pos );
                    assert ( pp[j].isDead );
                }
            }
            pp[i].isMerged = false;   
        }
    }
    PS::Comm::barrier();
    edisp += PS::Comm::getSum(edisp_loc);
    
    if ( n_remove ){
        pp.removeParticle(remove, n_remove);
    }
    delete [] remove;
}

template <class Tpsys>
PS::S32 removeParticlesOutOfBoundary(Tpsys & pp,
                                     PS::F64 & edisp,
                                     const PS::F64 r_max,
                                     const PS::F64 r_min,
                                     std::ofstream & fout_rem)
{
    const PS::F64 rmax2 = r_max*r_max;
    const PS::F64 rmin2 = r_min*r_min;
    PS::F64 edisp_loc = 0.;
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();

    static std::vector<PS::S32> n_remove_list;
    static std::vector<PS::S32> n_remove_adr;
    static std::vector<FP_t> remove_list_loc;
    static std::vector<FP_t> remove_list_glb;
    n_remove_list.resize(n_proc);
    n_remove_adr.resize(n_proc);

    static std::vector<PS::S32> remove_list;
    remove_list.clear();

#ifdef INDIRECT_TERM
    PS::F64 e_ind_before = 0.;
    PS::F64 e_ind_after  = 0.;
#endif
    
#pragma omp parallel for
    for ( PS::S32 i=0; i<n_loc; i++ ){
        PS::F64vec posi = pp[i].pos;
        PS::F64    pos2 = posi*posi;
        if ( pos2 > rmax2 || pos2 < rmin2 ){
#pragma omp critical
            {
                remove_list.push_back(i);
            }
        }
    }

    PS::S32 n_remove_loc = remove_list.size();
    PS::S32 n_remove_glb = PS::Comm::getSum(n_remove_loc);

    /*if ( n_remove_glb == 1 ){

        if ( n_remove_loc ) {
            PS::S32 i_remove = remove_list.at(0);
            
            PS::F64    massi = pp[i_remove].mass;
            PS::F64vec veli = pp[i_remove].vel;
            edisp_loc -= 0.5*massi* veli*veli;
            edisp_loc -= massi * pp[i_remove].phi_s;
            edisp_loc -= massi * pp[i_remove].phi_d;
            edisp_loc -= massi * pp[i_remove].phi;
        
            std::cerr << "Remove Particle " << pp[i_remove].id << std::endl
                      << "Position : " << std::setprecision(15) << pp[i_remove].pos << std::endl;
            fout_rem << std::fixed<<std::setprecision(8)
                     << pp[i_remove].time << "\t" << pp[i_remove].id << "\t"
                     << std::scientific << std::setprecision(15) << pp[i_remove].mass << "\t"
                     << pp[i_remove].pos.x << "\t" << pp[i_remove].pos.y << "\t" << pp[i_remove].pos.z << "\t"
                     << pp[i_remove].vel.x << "\t" << pp[i_remove].vel.y << "\t" << pp[i_remove].vel.z
                     << std::endl;
        }
        
        } else if ( n_remove_glb > 1 ){*/
    
    if ( n_remove_glb ){
        
        //PS::S32 * n_remove_list   = nullptr;
        //PS::S32 * n_remove_adr    = nullptr;
        //FP_t *  remove_list_loc = nullptr;
        //FP_t *  remove_list_glb = nullptr;
        
        if ( PS::Comm::getRank() == 0 ){
            //n_remove_list   = new PS::S32[n_proc];
            //n_remove_adr    = new PS::S32[n_proc];
            //remove_list_glb = new FP_t[n_remove_glb];
            remove_list_glb.resize(n_remove_glb);
        }
        //remove_list_loc = new FP_t[n_remove_loc];
        remove_list_loc.resize(n_remove_loc);
        
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Gather(&n_remove_loc, 1, PS::GetDataType(n_remove_loc),
                   &n_remove_list[0], 1, PS::GetDataType(n_remove_list[0]),  0, MPI_COMM_WORLD);
#else
        n_remove_list[0]  = n_remove_loc;
#endif
        //PS::Comm::gather(&n_remove_loc, 1, n_remove_list);
        
        if ( PS::Comm::getRank() == 0 ){
            PS::S32 tmp_remove = 0;
            for ( PS::S32 i=0; i<n_proc; i++ ){
                n_remove_adr[i]  = tmp_remove;
                tmp_remove += n_remove_list[i];
            }
            assert ( n_remove_glb == tmp_remove );
        }

        for ( PS::S32 i=0; i<n_remove_loc; i++ ) {
            remove_list_loc[i] = pp[remove_list.at(i)];
        }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        MPI_Gatherv(&remove_list_loc[0], n_remove_loc,                        PS::GetDataType(remove_list_loc[0]),
                    &remove_list_glb[0], &n_remove_list[0], &n_remove_adr[0], PS::GetDataType(remove_list_glb[0]), 0, MPI_COMM_WORLD);
#else
        for(PS::S32 i=0; i<n_remove_loc; i++) remove_list_glb[i] = remove_list_loc[i];
#endif
        //PS::Comm::gatherV(remove_list_loc, n_remove_loc, remove_list_glb, n_remove_list, n_remove_adr);
        
        if ( PS::Comm::getRank() == 0 ){
            for ( PS::S32 i=0; i<n_remove_glb; i++ ) {
            
                PS::F64    massi = remove_list_glb[i].mass;
                PS::F64vec veli  = remove_list_glb[i].vel;
                edisp_loc -= 0.5*massi* veli*veli;
                edisp_loc -= massi * remove_list_glb[i].phi_s;
                edisp_loc -= massi * remove_list_glb[i].phi_d;
                edisp_loc -= massi * remove_list_glb[i].phi;
            
                for ( PS::S32 j=0; j<i; j++ ) {
                    if ( remove_list_glb[i].id != remove_list_glb[j].id ) {
                        PS::F64    massj = remove_list_glb[j].mass;
                        PS::F64vec posi  = remove_list_glb[i].pos;
                        PS::F64vec posj  = remove_list_glb[j].pos;
                        PS::F64    eps2   = FP_t::eps2;
                        
                        PS::F64vec dr = posi - posj;
                        PS::F64    rinv = 1./sqrt(dr*dr + eps2);
                        
                        edisp_loc += - massi * massj * rinv;
                    }
                }
        
                std::cerr << "Remove Particle " << remove_list_glb[i].id << std::endl
                          << "Position : " << std::setprecision(15) << remove_list_glb[i].pos << std::endl;
                fout_rem << std::fixed<<std::setprecision(8)
                         << remove_list_glb[i].time << "\t" << remove_list_glb[i].id << "\t"
                         << std::scientific << std::setprecision(15) << remove_list_glb[i].mass << "\t"
                         << remove_list_glb[i].pos.x << "\t" << remove_list_glb[i].pos.y << "\t" << remove_list_glb[i].pos.z << "\t"
                         << remove_list_glb[i].vel.x << "\t" << remove_list_glb[i].vel.y << "\t" << remove_list_glb[i].vel.z
                         << std::endl;
            }
            
            //delete [] n_remove_list;
            //delete [] n_remove_adr;
            //delete [] remove_list_glb;
        
        }
        
        //delete [] remove_list_loc;

#ifdef INDIRECT_TERM
        e_ind_before = calcIndirectEnergy(pp);
#endif
    }
    
    if (n_remove_loc) pp.removeParticle(&remove_list[0], n_remove_loc);
            
    edisp += PS::Comm::getSum(edisp_loc);

#ifdef INDIRECT_TERM
    if (n_remove_glb) {
        e_ind_after = calcIndirectEnergy(pp);
        edisp += e_ind_after - e_ind_before;
    }
#endif

    return n_remove_glb;
}

template <class Tpsys>
void correctEnergyForGas(Tpsys & pp,
                         PS::F64 & edisp_gd,
                         bool second)
{// energy correction for gas drag
    PS::F64 edisp_gd_loc = 0.;
    PS::F64 coef = 0.25; if (second) coef *= -1.;
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    
#pragma omp parallel for reduction(+:edisp_gd_loc)
    for(PS::S32 i=0; i<n_loc; i++){
        edisp_gd_loc += pp[i].mass * pp[i].acc_gd
            * (pp[i].vel + coef * pp[i].acc_gd * FP_t::dt_tree);
    }
    edisp_gd += 0.5 * FP_t::dt_tree * PS::Comm::getSum(edisp_gd_loc);
}
