#pragma once

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
              const PS::S32 isnap)
{
    FileHeader header(pp.getNumberOfParticleGlobal(), time_sys, e_init, e_now);
    char filename[256];
    sprintf(filename, "%s/snap%06d.dat", dir_name, isnap);
    pp.writeParticleAscii(filename, header);
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
                std::ofstream & fout_eng,
                Wtime wtime,
                PS::S32 n_largestcluster,
                PS::S32 n_cluster,
                PS::S32 n_isoparticle,
                bool bSnap=true)
{
    const PS::S32 n_tot = pp.getNumberOfParticleGlobal();

    if ( bSnap ) makeSnap(pp, time_sys, e_init, e_now, dir_name, isnap);

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
                  <<std::endl;
    }
}

template <class Tpsys>
void inputIDLocalAndMyrank(Tpsys & pp)
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
}

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
void removeOutOfBoundaryParticle(Tpsys & pp,
                                 PS::F64 & edisp,
                                 const PS::F64 r_max,
                                 const PS::F64 r_min,
                                 std::ofstream & fout_rem)
{
    PS::F64 rmax2 = r_max*r_max;
    PS::F64 rmin2 = r_min*r_min;
    PS::F64 edisp_loc = 0.;
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();

    PS::S32 i_remove = -1;
    
#pragma omp parallel for
    for ( PS::S32 i=0; i<n_loc; i++ ){
        PS::F64vec posi = pp[i].pos;
        PS::F64    pos2 = posi*posi;
        if ( pos2 > rmax2 || pos2 < rmin2 ){
#pragma omp critical
            {
                if ( pos2 > rmax2 ) rmax2 = pos2;
                if ( pos2 < rmin2 ) rmin2 = pos2;
                i_remove = i;
            }
        }
    }

    if ( i_remove > -1 ){
        PS::F64    massi = pp[i_remove].mass;
        PS::F64vec veli = pp[i_remove].vel;
        edisp_loc -= 0.5*massi* veli*veli;
        edisp_loc -= massi * pp[i_remove].phi_s;
        edisp_loc -= massi * pp[i_remove].phi_d;
        edisp_loc -= massi * pp[i_remove].phi;

        std::cout << "Remove Particle " << i_remove << std::endl
                  << "Position : " << std::setprecision(15) << pp[i_remove].pos << std::endl;
        fout_rem << pp[i_remove].time << "\t" << pp[i_remove].id << "\t" << pp[i_remove].mass << "\t"
                 << pp[i_remove].pos.x << "\t" << pp[i_remove].pos.y << "\t" << pp[i_remove].pos.z << "\t"
                 << pp[i_remove].vel.x << "\t" << pp[i_remove].vel.y << "\t" << pp[i_remove].vel.z
                 << std::endl;
        pp.removeParticle(&i_remove, 1);
    }
    edisp += PS::Comm::getSum(edisp_loc);
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
            * (pp[i].vel + coef * pp[i].acc_gd * FPGrav::dt_tree);
    }
    edisp_gd += 0.5 * FPGrav::dt_tree * PS::Comm::getSum(edisp_gd_loc);
}
