/*
system_soft[i].n_epj_: without particle itself
pair_hadr_loc_: without particle itself
epj_ngb_arra_ has itself.
*/

/*
calcAcc0Acc1AndR2()
*/

/*
dead particle has zero-mass
 */
#ifdef USE_INTRINSIC_FOR_X86
#include<immintrin.h>
#endif

#include"kepler.hpp"

class SortAdr{
public:
    static PS::F64 * time;
    bool operator() (const PS::S32 & left, const PS::S32 & right) const {
        return time[left] < time[right];
    }
};

PS::F64 * SortAdr::time;

// in the case of merge, dt_limit must be changed. 
inline PS::F64 CalcDtLimit(const PS::F64 time_sys, 
                           const PS::F64 time_sync, 
                           const PS::F64 dt_limit_org){
    PS::F64 dt_limit_ret = dt_limit_org;
    PS::F64 s = time_sys / dt_limit_ret;
    //while(s != PS::F64(PS::S32(s))){
    while(s != PS::F64(PS::S64(s))){
        s *= 2.0;
        dt_limit_ret *= 0.5;
    }
    return dt_limit_ret;
}

#include"force.hpp"
class HardSystem{
private:
    void makeDictionary(std::vector<PTCLHard> & ptcl,
                        const std::vector< std::pair<PS::S64, PS::S64> > & pair_id,
                        std::unordered_map<PS::S64, PS::S64> & idx_to_adr,
                        std::vector< std::pair<PS::S64, PS::S64> > & pair_adr){
        idx_to_adr.clear();
        pair_adr.clear();
        const PS::S32 n_ptcl = ptcl.size();
        for(PS::S32 i=0; i<n_ptcl; i++){
            idx_to_adr.insert( std::pair<PS::S64, PS::S64>(ptcl[i].id, i) );
        }
        const PS::S32 n_interaction = pair_id.size();
        PS::S32 first_adr_old = -1;
        for(PS::S32 i=0; i<n_interaction; i++){
            std::unordered_map<PS::S64, PS::S64>::iterator itr = idx_to_adr.find(pair_id[i].second);
            const PS::S64 first_adr = idx_to_adr[pair_id[i].first];
            if(first_adr_old != first_adr){
                ptcl[first_adr].adr_pair = i;
                first_adr_old = first_adr;
            }
            if(itr == idx_to_adr.end()){
                pair_adr.push_back( std::pair<PS::S64, PS::S64>(first_adr, -1));
            }
            else{
                pair_adr.push_back( std::pair<PS::S64, PS::S64>(first_adr,
                                                                idx_to_adr[pair_id[i].second]));
            }
        }
    }

    PS::F64 mass_sun_;
    PS::F64vec pos_sun_;
    PS::F64vec vel_sun_;

    template<class Tptcl>
    void merge2body(Tptcl & ptcl0, 
		    Tptcl & ptcl1, 
		    PS::F64 & eng_disp, 
		    std::vector<MergeLog> & merge_history){
        const PS::F64 m0 = ptcl0.mass;
        const PS::F64 m1 = ptcl1.mass;
        const PS::S64 id0 = ptcl0.id;
        const PS::S64 id1 = ptcl1.id;

        Tptcl & ptcl_merge = (m0>m1) ? ptcl0 : ( (m1>m0) ? ptcl1 : ((id0<id1) ? ptcl0 : ptcl1) );
        Tptcl & ptcl_dead = (&ptcl_merge == &ptcl0) ? ptcl1 : ptcl0;

        const PS::F64vec r01 = ptcl0.pos - ptcl1.pos;
        const PS::F64vec v01 = ptcl0.vel - ptcl1.vel;
        const PS::F64vec r0 = ptcl0.pos - pos_sun_;
        const PS::F64vec r1 = ptcl1.pos - pos_sun_;
        const PS::F64 pot0 = -mass_sun_ * m0 / sqrt(r0*r0);
        const PS::F64 pot1 = -mass_sun_ * m1 / sqrt(r1*r1);
        const PS::F64 pot_old = pot0 + pot1 - m0*m1/sqrt(r01*r01);
        const PS::F64 m_red = (m0*m1) / (m0+m1);

        const PS::F64 merge_time = ptcl0.time + time_origin_;
        merge_history.push_back( MergeLog(merge_time, ptcl_merge, ptcl_dead) );
        ptcl_merge.merge(ptcl_dead);
        ptcl_dead.mass = 0.0;

        const PS::F64vec r_new = ptcl_merge.pos - pos_sun_;
        const PS::F64 pot_new = -mass_sun_ * ptcl_merge.mass / sqrt(r_new*r_new); 
        eng_disp += 0.5*m_red*v01*v01 - (pot_new - pot_old);
/*
        std::cerr<<"eng_disp="<<eng_disp<<" pot_new="<<pot_new
                 <<" pot_old="<<pot_old<<std::endl;
        std::cerr<<"ptcl_merge.n_ngb="<<ptcl_merge.n_ngb<<" ptcl_dead.n_ngb="<<ptcl_dead.n_ngb<<std::endl;
*/
    }

    template<class Tptcl>
    void merge2body(std::vector<Tptcl> & ptcl, 
                    const std::vector< std::pair<PS::S64, PS::S64> > pair_adr,
                    const PS::S64 adr0, 
		    const PS::S64 adr1, 
                    PS::F64 & eng_disp, 
		    std::vector<MergeLog> & merge_history){
        Tptcl & ptcl0 = ptcl[adr0];
        Tptcl & ptcl1 = ptcl[adr1];
        const PS::F64 m0 = ptcl0.mass;
        const PS::F64 m1 = ptcl1.mass;
        const PS::S64 id0 = ptcl0.id;
        const PS::S64 id1 = ptcl1.id;
        Tptcl & ptcl_merge = (m0>m1) ? ptcl0 : ( (m1>m0) ? ptcl1 : ((id0<id1) ? ptcl0 : ptcl1) );
        Tptcl & ptcl_dead = (&ptcl_merge == &ptcl0) ? ptcl1 : ptcl0;

        std::set<PS::S64> adr_jp;
        adr_jp.clear();

        const PS::S32 n0 = ptcl0.n_ngb;
        const PS::S32 adr_pair_ngb_0 = ptcl0.adr_pair;
        for(PS::S32 j=0; j<n0; j++){
            PS::S32 adr = pair_adr[adr_pair_ngb_0+j].second;
            const Tptcl & ptclj = ptcl[adr];
            if( ptcl1.id == ptclj.id) continue;
            adr_jp.insert(adr);
        }
        const PS::S32 n1 = ptcl1.n_ngb;
        const PS::S32 adr_pair_ngb_1 = ptcl1.adr_pair;
        for(PS::S32 j=0; j<n1; j++){
            PS::S32 adr = pair_adr[adr_pair_ngb_1+j].second;
            const Tptcl & ptclj = ptcl[adr];
            if( ptcl0.id == ptclj.id) continue;
            adr_jp.insert(adr);
        }
        PS::F64 pot_old = 0.0;
        const size_t n_ngb = adr_jp.size();
        std::set<PS::S64>::iterator adr_itr = adr_jp.begin();

#if 0
        for(size_t j=0; j<n_ngb; j++, adr_itr++){
            const Tptcl & ptclj = ptcl[*adr_itr];
            //std::cerr<<"---------------"<<std::endl;
            //std::cerr<<"---------------"<<std::endl;
            if(ptcl0.id != ptclj.id && ptcl1.id != ptclj.id){
                //std::cerr<<"ptcl0.id="<<ptcl0.id<<" ptcl1.id="<<ptcl1.id<<" ptclj.id="<<ptclj.id<<std::endl;
                //std::cerr<<"ptcl0.pos="<<ptcl0.pos<<" ptcl1.pos="<<ptcl1.pos<<" ptclj.pos="<<ptclj.pos<<std::endl;
                const PS::F64vec dr0 = ptcl0.pos - ptclj.pos;
                pot_old -= (ptcl0.mass * ptclj.mass) / sqrt(dr0*dr0);
                const PS::F64vec dr1 = ptcl1.pos - ptclj.pos;
                pot_old -= (ptcl1.mass * ptclj.mass) / sqrt(dr1*dr1);
                //std::cerr<<"dr0="<<dr0<<" dr1="<<dr1<<std::endl;
            }
            //std::cerr<<"---------------"<<std::endl;
            //std::cerr<<"---------------"<<std::endl;
        }
#endif

        const PS::F64vec r0 = ptcl0.pos - pos_sun_;
        pot_old -= mass_sun_ * m0 / sqrt(r0*r0);
        const PS::F64vec r1 = ptcl1.pos - pos_sun_;
        pot_old -= mass_sun_ * m1 / sqrt(r1*r1);
        const PS::F64vec r01 = ptcl0.pos - ptcl1.pos;
        pot_old -= m0*m1/sqrt(r01*r01);
        const PS::F64vec v01 = ptcl0.vel - ptcl1.vel;
        const PS::F64 m_red = (m0*m1) / (m0+m1);
        PS::F64 kin = 0.5*m_red*v01*v01;

        const PS::F64 merge_time = ptcl0.time + time_origin_;
        merge_history.push_back( MergeLog(merge_time, ptcl_merge, ptcl_dead) );
        ptcl_merge.merge(ptcl_dead);
        ptcl_dead.mass = 0.0;

        const PS::F64vec r_new = ptcl_merge.pos - pos_sun_;
        PS::F64 pot_new = -mass_sun_ * ptcl_merge.mass / sqrt(r_new*r_new); 
        adr_itr = adr_jp.begin();
#if 0
        //std::cerr<<"---------------"<<std::endl;
        //std::cerr<<"---------------"<<std::endl;
        for(size_t j=0; j<n_ngb; j++, adr_itr++){
            Tptcl & ptclj = ptcl[*adr_itr];
            if(ptcl0.id != ptclj.id && ptcl1.id != ptclj.id){
                const PS::F64vec dr = ptcl_merge.pos - ptclj.pos;
                pot_new -= (ptcl_merge.mass * ptclj.mass) / sqrt(dr*dr);
                //std::cerr<<"dr="<<dr<<std::endl;
            }
        }
#endif
        //std::cerr<<"---------------"<<std::endl;
        //std::cerr<<"---------------"<<std::endl;
        eng_disp += kin - (pot_new - pot_old);
/*
        std::cerr<<"eng_disp="<<eng_disp<<" pot_new="<<pot_new
                 <<" pot_old="<<pot_old<<" kin="<<kin<<std::endl;
        std::cerr<<"ptcl_merge.n_ngb="<<ptcl_merge.n_ngb<<" ptcl_dead.n_ngb="<<ptcl_dead.n_ngb<<std::endl;
*/
    }

    template<class Tptcl>
    void merge2body(Tptcl & ptcl0, Tptcl & ptcl1, std::vector<MergeLog> & merge_history){
        const PS::F64 m0 = ptcl0.mass;
        const PS::F64 m1 = ptcl1.mass;
        const PS::S64 id0 = ptcl0.id;
        const PS::S64 id1 = ptcl1.id;
        Tptcl & ptcl_merge = (m0>m1) ? ptcl0 : ( (m1>m0) ? ptcl1 : ((id0<id1) ? ptcl0 : ptcl1) );
        Tptcl & ptcl_dead = (&ptcl_merge == &ptcl0) ? ptcl1 : ptcl0;
        const PS::F64 merge_time = ptcl0.time + time_origin_;
        merge_history.push_back( MergeLog(merge_time, ptcl_merge, ptcl_dead) );
        ptcl_merge.merge(ptcl_dead);
        ptcl_dead.mass = 0.0;
    }

    template<class Tptcl>
    void merge2body(Tptcl & ptcl0, Tptcl & ptcl1){
        const PS::F64 m0 = ptcl0.mass;
        const PS::F64 m1 = ptcl1.mass;
        const PS::S64 id0 = ptcl0.id;
        const PS::S64 id1 = ptcl1.id;
        Tptcl & ptcl_merge = (m0>m1) ? ptcl0 : ( (m1>m0) ? ptcl1 : ((id0<id1) ? ptcl0 : ptcl1) );
        Tptcl & ptcl_dead = (&ptcl_merge == &ptcl0) ? ptcl1 : ptcl0;
        ptcl_merge.merge(ptcl_dead);
        ptcl_dead.mass = 0.0;
    }

    void calcAcc0AndAcc1FromSun(const PS::F64vec & pos,
                                const PS::F64vec & vel,
                                PS::F64vec & acc0,
                                PS::F64vec & acc1){
        const PS::F64vec rij = pos - pos_sun_;
        const PS::F64vec vij = vel - vel_sun_;
        const PS::F64 r_inv = 1.0 / sqrt(rij * rij);
        const PS::F64 r2_inv = r_inv * r_inv;
        const PS::F64 r3_inv = r2_inv * r_inv;
        const PS::F64 m_r3 = mass_sun_ * r3_inv;
        const PS::F64vec F0 = -m_r3*rij;
        const PS::F64vec F1 = -m_r3*vij - 3.0*rij*vij*r2_inv*F0;
        acc0 += F0;
        acc1 += F1;
    }

public:
    class Wtime{
    public:
        PS::F64 predict_;
        PS::F64 force_;
        PS::F64 force_kepler_;
        PS::F64 correct_;
        PS::F64 sort_;
        void clear(){
            predict_ = force_ = force_kepler_ = correct_ = sort_ = 0.0;
        }
        void dump(std::ostream & fout){
            fout<<"predict_= "<<predict_<<" force_= "<< force_<<" force_kepler_= "<<force_kepler_
                <<" correct_= "<<correct_<<" sort_= "<<sort_<<std::endl;
        }
    } wtime_profile_;

    void setSun(const PS::F64 m, const PS::F64vec & p, const PS::F64vec & v){
        mass_sun_ = m;
        pos_sun_ = p;
        vel_sun_ = v;
    }

    void clearEngDisp(){
        eng_disp_ = 0.0;
    }

    std::vector<PS::S32> adr_ptcl_loc_;
    std::vector< std::pair<PS::S64, PS::S64> > pair_id_loc_;
    //std::vector<EPJSoft*> epj_ngb_array_ptr_;
    std::vector<EPJSoft> epj_ngb_array_;
    std::vector<PS::S32> adr_head_epj_ngb_;
    std::vector<PTCLHard> ptcl_loc_;
    PS::S64 n_loop_;
    PS::F64 a0_offset_sq_;
    PS::F64 time_origin_;
    template<class Ttree, class Tsystem>
    void setUp(Ttree & tree_soft,
               Tsystem & system_soft,
               const PS::S32 n_tot_loc, 
               const PS::F64 r_out,
               const PS::F64 r_in,
               const PS::F64 time_origin=0.0){
        time_origin_ = time_origin;
        adr_ptcl_loc_.clear();
        for(PS::S32 i=0; i<n_tot_loc; i++){
            if(tree_soft.getForce(i).n_ngb > 0){
                adr_ptcl_loc_.push_back(i);
            }
        }
        const PS::S32 n_loc = adr_ptcl_loc_.size();
	/*
        if(epj_ngb_array_ptr_.capacity() < static_cast<size_t>(n_loc) ){
            epj_ngb_array_ptr_.reserve(n_loc + (n_loc+3)/3);
        }
        epj_ngb_array_ptr_.clear();
        epj_ngb_array_ptr_.resize(n_loc);
	*/
	if(epj_ngb_array_.capacity() < static_cast<size_t>(n_loc) ){
            epj_ngb_array_.reserve(n_loc + (n_loc+3)/3);
        }
        if(ptcl_loc_.capacity() < static_cast<size_t>(n_loc) ){
            ptcl_loc_.reserve(n_loc + (n_loc+3)/3);
        }
        ptcl_loc_.clear();
        ptcl_loc_.resize(n_loc);
	
        for(PS::S32 i=0; i<n_loc; i++){
            const PS::S32 adr = adr_ptcl_loc_[i];
            ptcl_loc_[i].id = system_soft[adr].id;
            ptcl_loc_[i].mass = system_soft[adr].mass;
            ptcl_loc_[i].pos = system_soft[adr].pos;
            ptcl_loc_[i].vel = system_soft[adr].vel;
            ptcl_loc_[i].n_ngb = tree_soft.getForce(adr).n_ngb;
        }
        pair_id_loc_.clear();
#if 1
	//std::cout<<"check a"<<std::endl;
	// 2016 02/05 new version of neighbour search
	adr_head_epj_ngb_.clear();
	epj_ngb_array_.clear();
	PS::S32 offset = 0;
	for(PS::S32 i=0; i<n_loc; i++){
	    adr_head_epj_ngb_.push_back(offset);
            const PS::S32 adr = adr_ptcl_loc_[i];
	    EPJSoft * nbl = NULL;
            PS::S32 n_ngb = tree_soft.getNeighborListOneParticle(system_soft[adr], nbl);
	    //std::cout<<"i="<<i<<" adr="<<adr<<" n_ngb="<<n_ngb<<" offset="<<offset<<std::endl;
	    PS::S32 j2 = 0;
	    for(PS::S32 j=0; j<n_ngb; j++){
		if(system_soft[adr].id == (nbl+j)->id) continue;
		epj_ngb_array_.push_back(nbl[j]);
		/*
		std::cout<<"epj_ngb_array_.size()="<<epj_ngb_array_.size()<<" offset+j2="<<offset+j2<<std::endl;
		std::cout<<"nbl[j].id="<<nbl[j].id<<std::endl;
		std::cout<<"epj_ngb_array_[offset+j2].id="<<epj_ngb_array_[offset+j2].id<<std::endl;
		std::cout<<"nbl[j].pos="<<nbl[j].pos<<std::endl;
		std::cout<<"epj_ngb_array_[offset+j2].pos="<<epj_ngb_array_[offset+j2].pos<<std::endl;
		*/
		pair_id_loc_.push_back(std::pair<PS::S64, PS::S64>( (PS::S64)system_soft[adr].id, (PS::S64)(nbl+j)->id));
		j2++;
	    }
	    offset += n_ngb - 1; // subtraction means reself interaction
	    //std::cout<<std::endl;
        } // end of for over i
	//std::cout<<"check b"<<std::endl;
#else	
        tree_soft.epj_neighbor_->clearSize();
        for(PS::S32 i=0; i<n_loc;){
            const PS::S32 adr = adr_ptcl_loc_[i];
            PS::S32 n_ngb = tree_soft.getNeighborListOneParticle(system_soft[adr], epj_ngb_array_ptr_[i] );
            if(n_ngb < 0){
                // do again
                std::cerr<<"do again"<<std::endl;
                std::cerr<<"tree_soft.epj_neighbor_->capacity()="<<tree_soft.epj_neighbor_->capacity()<<std::endl;
                pair_id_loc_.clear();
                tree_soft.epj_neighbor_->clearSize();
                i = 0;
            } 
            else{
                for(PS::S32 j=0; j<n_ngb; j++){
                    if(system_soft[adr].id == (epj_ngb_array_ptr_[i]+j)->id) continue;
                    pair_id_loc_.push_back(std::pair<PS::S64, PS::S64>( (PS::S64)system_soft[adr].id, 
                                                                        (PS::S64)(epj_ngb_array_ptr_[i]+j)->id));
                }
                i++; // improtant.
            }
        } // end of for over i
#endif
	
#ifdef DEBUG_PRINT_PLANET
	std::cerr<<"check 1"<<std::endl;
#endif
        /*
          std::cerr<<"pair_id_loc_.size()= "<<pair_id_loc_.size()<<std::endl;
          for(PS::S32 i=0; i<pair_id_loc_.size(); i++){
          std::cerr<<"pair_id_loc_[i].first="<<pair_id_loc_[i].first<<" pair_id_loc_[i].second="<<pair_id_loc_[i].second<<std::endl;
          }
        */

        const PS::F64 eps_sq = EPISoft::eps * EPISoft::eps;
        const PS::F64 r_search_sq = EPJSoft::r_search * EPJSoft::r_search;
	//std::cout<<"check c"<<std::endl;
        for(PS::S32 i=0; i<n_loc; i++){
            const PS::S32 adr = adr_ptcl_loc_[i];
            const PS::F64vec pos_i = system_soft[adr].getPos();
            PS::F64vec & acc_long_i = system_soft[adr].acc; // accumulate acc
            PS::F64vec & acc_short_i = ptcl_loc_[i].acc0; // accumulate acc
            PS::F64 & pot_tot_i = system_soft[adr].pot_tot; // accumulate pot
	    //EPJSoft * epj_tmp = epj_ngb_array_ptr_[i];
	    //const PS::S32 n_ngb = tree_soft.getForce(adr).n_ngb + 1;
	    const PS::S32 n_ngb = tree_soft.getForce(adr).n_ngb;
	    EPJSoft * epj_tmp = &(epj_ngb_array_[adr_head_epj_ngb_[i]]);
	    //std::cout<<"adr_head_epj_ngb_[i]="<<adr_head_epj_ngb_[i]<<" n_ngb="<<n_ngb<<std::endl;
            for(PS::S32 j=0; j<n_ngb; j++){
                if(system_soft[adr].id == (epj_tmp+j)->id) continue;
                calc_acc_split(pos_i,          acc_long_i, 
                               acc_short_i,    pot_tot_i, 
                               epj_tmp[j].pos, epj_tmp[j].mass,
                               eps_sq,         r_out, 
                               r_in,           r_search_sq);
            }
	    // for debug
	    system_soft[adr].acc_pla += acc_long_i + acc_short_i;
        }
	//std::cout<<"check d"<<std::endl;
#ifdef DEBUG_PRINT_PLANET
	std::cerr<<"check 2"<<std::endl;
#endif

	/// set a0 offset
	PS::F64 m_min = PS::LARGE_FLOAT;
	const PS::S32 n = ptcl_loc_.size();
	for(PS::S32 i=0; i<n; i++){
	    if(m_min > ptcl_loc_[i].mass && ptcl_loc_[i].mass > 0.0) m_min = ptcl_loc_[i].mass;
	}
#ifdef DEBUG_PRINT_PLANET
	std::cerr<<"m_min="<<m_min<<std::endl;
#endif

	PS::F64 a0_offset_sq_loc = 0.1 * m_min / (r_out*r_out);
	//PS::F64 a0_offset_sq_loc = 0.0;
	a0_offset_sq_ = PS::Comm::getMinValue(a0_offset_sq_loc);

#ifdef DEBUG_PRINT_PLANET
	std::cerr<<"a0_offset_sq_="<<a0_offset_sq_<<std::endl;
#endif

    }

    
    template<class Tsystem>
    void copyVelSoftToHard(const Tsystem & system){
	const PS::S32 n = adr_ptcl_loc_.size();
#pragma omp parallel for
	for(PS::S32 i=0; i<n; i++){
	    const PS::S32 adr = adr_ptcl_loc_[i];
	    ptcl_loc_[i].vel = system[adr].vel;
	}
    }
    
    std::unordered_map<PS::S64, PS::S64> idx_to_adr_loc_;
    std::vector< std::pair<PS::S64, PS::S64> > pair_adr_loc_;

    void makeDictionaryLocal(){
        makeDictionary(ptcl_loc_, pair_id_loc_, idx_to_adr_loc_, pair_adr_loc_);
    }


    std::vector<PS::S32> n_interaction_array_;
    std::vector<PS::S32> n_ptcl_array_;
    std::vector< std::pair<PS::S32, PS::S32> > n_int_n_ptcl_array_;
    std::vector<PS::S32> n_interaction_disp_;
    std::vector<PS::S32> n_ptcl_disp_;
    std::vector<PTCLHard> ptcl_multi_glb_;
    std::vector< std::pair<PS::S64, PS::S64> > pair_id_multi_glb_;
    void gatherData(){
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        n_int_n_ptcl_array_.reserve(n_proc);
        std::pair<PS::S32, PS::S32> n_int_n_ptcl(pair_id_multi_loc_.size(), ptcl_multi_loc_.size());
        PS::Comm::allGather(&n_int_n_ptcl, 1, &n_int_n_ptcl_array_[0]);
        n_interaction_array_.reserve(n_proc);
        n_ptcl_array_.reserve(n_proc);
        n_interaction_disp_.reserve(n_proc+1);
        n_ptcl_disp_.reserve(n_proc+1);
        n_interaction_disp_[0] = n_ptcl_disp_[0] = 0;
        for(PS::S32 i=0; i<n_proc; i++){
            n_interaction_array_[i] = n_int_n_ptcl_array_[i].first;
            n_ptcl_array_[i] = n_int_n_ptcl_array_[i].second;
            n_interaction_disp_[i+1] = n_interaction_disp_[i] + n_interaction_array_[i];
            n_ptcl_disp_[i+1] = n_ptcl_disp_[i] + n_ptcl_array_[i];
        }

        pair_id_multi_glb_.resize(n_interaction_disp_[n_proc]);
        ptcl_multi_glb_.resize(n_ptcl_disp_[n_proc]);

        PS::Comm::gatherV(&pair_id_multi_loc_[0], pair_id_multi_loc_.size(), &pair_id_multi_glb_[0],
                          &n_interaction_array_[0], &n_interaction_disp_[0]);

        PS::Comm::gatherV(&ptcl_multi_loc_[0], ptcl_multi_loc_.size(),  &ptcl_multi_glb_[0],
                          &n_ptcl_array_[0], &n_ptcl_disp_[0]);

    }

    std::unordered_map<PS::S64, PS::S64> idx_to_adr_multi_glb_;
    std::vector< std::pair<PS::S64, PS::S64> > pair_adr_multi_glb_;
    void makeDictionaryGlobal(){
        //makeDictionary(ptcl_glb_, pair_id_glb_, idx_to_adr_glb_, pair_adr_glb_);
        makeDictionary(ptcl_multi_glb_, pair_id_multi_glb_, idx_to_adr_multi_glb_, pair_adr_multi_glb_);
    }

    std::vector<PTCLPred> ptcl_pred_;
    std::vector<PTCLForce> ptcl_force_;
    std::vector<PS::F64> time_next_;
    std::vector<PS::S32> adr_sorted_;
    PS::F64 time_sys_;
    PS::F64 time_end_;
    PS::F64 dt_limit_;
    PS::F64 time_sync_;
    PS::S32 n_active_;

    std::vector<bool> merge_flag_glb_;
    bool merge_state_;
    void setUpGlb(const PS::F64 time_end, 
                  const PS::F64 dt_limit){
        const PS::S32 n = ptcl_multi_glb_.size();
        n_active_ = n;
        time_sys_ = 0.0;
        time_end_ = time_end;
        dt_limit_ = dt_limit;
        time_sync_ = time_sys_ + dt_limit_;
        adr_sorted_.clear();
        adr_sorted_.resize(n);
        ptcl_pred_.resize(n);
        ptcl_force_.resize(n);
        time_next_.resize(n);
        merge_flag_glb_.clear();
        merge_flag_glb_.resize(n);
        for(PS::S32 ip=0; ip<n; ip++){
            adr_sorted_[ip] = ip;
            ptcl_pred_[ip].pos = ptcl_multi_glb_[ip].pos;
            ptcl_pred_[ip].vel = ptcl_multi_glb_[ip].vel;
            ptcl_force_[ip].acc0 = ptcl_force_[ip].acc1 = 0.0;
            ptcl_multi_glb_[ip].time = time_next_[ip] = 0.0;
            ptcl_multi_glb_[ip].setRMerge();
            merge_flag_glb_[ip] = false;
        }
        merge_state_ = false;
        n_loop_ = 0;
        wtime_profile_.clear();
    }

    std::vector< std::pair<PS::S64, PS::S64> > pair_merge_adr_;
    void calcAcc0AndAcc1(const PS::F64 r_out,
                         const PS::F64 r_in){
        const PS::S32 n_thread_max = PS::Comm::getNumberOfThread();
        static std::vector< std::pair<PS::S64, PS::S64> > * pair_merge_adr_omp;
        static bool first = true;
        if(first){
            pair_merge_adr_omp = new std::vector< std::pair<PS::S64, PS::S64> >[n_thread_max];
            first = false;
        }
        for(PS::S32 i=0; i<n_thread_max; i++){
            pair_merge_adr_omp[i].clear();
        }
        pair_merge_adr_.clear();
        const PS::F64 eps_sq = EPISoft::eps * EPISoft::eps;
        const PS::S32 ni = n_active_;
#pragma omp parallel for
        for(PS::S32 ip=0; ip<ni; ip++){
            const PS::S32 ith = PS::Comm::getThreadNum();
            const PS::S32 adr_i = adr_sorted_[ip];
            ptcl_force_[adr_i].clear();
            const PS::S32 nj = ptcl_multi_glb_[adr_i].n_ngb;
            const PS::S32 adr_pair = ptcl_multi_glb_[adr_i].adr_pair;
            for(PS::S32 jp=0; jp<nj; jp++){
                PS::F64 r2 = 0.0;
                const PS::S32 adr_j = pair_adr_multi_glb_[adr_pair+jp].second;
                if(   ptcl_multi_glb_[adr_i].id != pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair].first
                      || ptcl_multi_glb_[adr_j].id != pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair+jp].second){
                    std::cerr<<"ptcl_multi_glb_[adr_i].id="<<ptcl_multi_glb_[adr_i].id<<" ptcl_multi_glb_[adr_j].id="<<ptcl_multi_glb_[adr_j].id<<std::endl;
                    std::cerr<<"pair_multi_id_glb_[ptcl_multi_glb_[adr_i].adr_pair].first="<<pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair].first
                             <<" pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair+jp].second="<<pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair+jp].second<<std::endl;
                    std::cerr<<std::endl;
                }
                assert(ptcl_multi_glb_[adr_i].id == pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair].first);
                assert(ptcl_multi_glb_[adr_j].id == pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair+jp].second);
#ifdef FORDEBUG
                CalcAcc0AndAcc1Cutoff(ptcl_pred_[adr_i].pos,   ptcl_pred_[adr_i].vel,
                                      ptcl_force_[adr_i].acc0, ptcl_force_[adr_i].acc1,
                                      ptcl_force_[adr_i].pot,  ptcl_pred_[adr_j].pos,
                                      ptcl_pred_[adr_j].vel,   ptcl_multi_glb_[adr_j].mass,
                                      eps_sq,                  r_out, r_in);		
#else // FORDEBUG
#ifdef MERGE
		/*
                CalcAcc0Acc1AndR2Cutoff(ptcl_pred_[adr_i].pos,   ptcl_pred_[adr_i].vel,
                                        ptcl_force_[adr_i].acc0, ptcl_force_[adr_i].acc1,
                                        r2, 
                                        ptcl_pred_[adr_j].pos,   ptcl_pred_[adr_j].vel,
                                        ptcl_multi_glb_[adr_j].mass,   eps_sq,
                                        r_out, r_in);
		*/
                CalcAcc0Acc1AndR2Cutoff(ptcl_pred_[adr_i].pos,   ptcl_pred_[adr_i].vel,
                                        ptcl_force_[adr_i].acc0_pla, ptcl_force_[adr_i].acc1_pla,
                                        r2, 
                                        ptcl_pred_[adr_j].pos,   ptcl_pred_[adr_j].vel,
                                        ptcl_multi_glb_[adr_j].mass,   eps_sq,
                                        r_out, r_in);
                const PS::F64 r_merge = (ptcl_multi_glb_[adr_i].r_merge + ptcl_multi_glb_[adr_j].r_merge);
                const PS::F64 r_merge_2 = r_merge * r_merge;
                if(r2 < r_merge_2){
                    //pair_merge_adr_.push_back( std::pair<PS::S64, PS::S64>(adr_i, adr_j) );
                    pair_merge_adr_omp[ith].push_back( std::pair<PS::S64, PS::S64>(adr_i, adr_j) );
                }
#else // MERGE
		/*
                CalcAcc0AndAcc1Cutoff(ptcl_pred_[adr_i].pos,   ptcl_pred_[adr_i].vel,
                                      ptcl_force_[adr_i].acc0, ptcl_force_[adr_i].acc1,
                                      ptcl_pred_[adr_j].pos,   ptcl_pred_[adr_j].vel,
                                      ptcl_multi_glb_[adr_j].mass,   eps_sq,
                                      r_out, r_in);
		*/
                CalcAcc0AndAcc1Cutoff(ptcl_pred_[adr_i].pos,   ptcl_pred_[adr_i].vel,
                                      ptcl_force_[adr_i].acc0_pla, ptcl_force_[adr_i].acc1_pla,
                                      ptcl_pred_[adr_j].pos,   ptcl_pred_[adr_j].vel,
                                      ptcl_multi_glb_[adr_j].mass,   eps_sq,
                                      r_out, r_in);
#endif // MERGE
#endif // FORDEBUG
            }
        }
#ifdef MERGE
        for(PS::S32 i=0; i<n_thread_max; i++){
            const size_t j_max = pair_merge_adr_omp[i].size();
            for(size_t j=0; j<j_max; j++){
                pair_merge_adr_.push_back(pair_merge_adr_omp[i][j]);
            }
        }
        const PS::S32 list_len = pair_merge_adr_.size();
        bool sort_again = false;
        for(PS::S32 i=0; i<list_len; i++){
            const PS::S32 adr_i = pair_merge_adr_[i].first;
            const PS::S32 adr_j = pair_merge_adr_[i].second;
            if( time_next_[adr_i] != time_next_[adr_j] && !ptcl_multi_glb_[adr_j].isDead()){
		// if the merging two particles are not integraed at the same time,
		// the particle non-integrated particle is forced to be integraed.
                std::cout<<"adr_i="<<adr_i
                         <<" adr_j="<<adr_j<<std::endl;
                std::cout<<"time_next_[adr_i]="<<time_next_[adr_i]
                         <<" time_next_[adr_j]="<<time_next_[adr_j]<<std::endl;
                std::cout<<"ptcl_multi_glb_[adr_i].n_ngb="<<ptcl_multi_glb_[adr_i].n_ngb<<std::endl;
                time_next_[adr_j] = time_next_[adr_i];
                ptcl_multi_glb_[adr_j].dt = time_next_[adr_j] - ptcl_multi_glb_[adr_j].time;
                ptcl_force_[adr_j].clear();
                const PS::S32 nk = ptcl_multi_glb_[adr_j].n_ngb;
                const PS::S32 adr_pair = ptcl_multi_glb_[adr_j].adr_pair;
                for(PS::S32 kp=0; kp<nk; kp++){
                    // Here, not care about second (j-th and k-th) merger.
                    const PS::S32 adr_k = pair_adr_multi_glb_[adr_pair+kp].second;
		    /*
                    CalcAcc0AndAcc1Cutoff(ptcl_pred_[adr_j].pos,   ptcl_pred_[adr_j].vel,
                                          ptcl_force_[adr_j].acc0, ptcl_force_[adr_j].acc1,
                                          ptcl_pred_[adr_k].pos,   ptcl_pred_[adr_k].vel,
                                          ptcl_multi_glb_[adr_k].mass,   eps_sq,
                                          r_out, r_in);
		    */
                    CalcAcc0AndAcc1Cutoff(ptcl_pred_[adr_j].pos,   ptcl_pred_[adr_j].vel,
                                          ptcl_force_[adr_j].acc0_pla, ptcl_force_[adr_j].acc1_pla,
                                          ptcl_pred_[adr_k].pos,   ptcl_pred_[adr_k].vel,
                                          ptcl_multi_glb_[adr_k].mass,   eps_sq,
                                          r_out, r_in);
                }
                n_active_++;
                sort_again = true;
            }
        }
        if(sort_again){
            const PS::S32 n = ptcl_multi_glb_.size();
            SortAdr::time = &time_next_[0];
            std::sort(&adr_sorted_[0], &adr_sorted_[0]+n, SortAdr());
        }
#endif
    }

    void calcAcc0AndAcc1Kepler(const PS::F64 mass_sun,
                               const PS::F64vec & pos_sun,
                               const PS::F64vec & vel_sun){
        const PS::S32 ni = n_active_;
#pragma omp parallel for
        for(PS::S32 ip=0; ip<ni; ip++){
            const PS::S32 adr = adr_sorted_[ip];
            PTCLPred & pred = ptcl_pred_[adr];
            PTCLForce & force = ptcl_force_[adr];
            const PS::F64vec rij = pred.pos - pos_sun;
            const PS::F64vec vij = pred.vel - vel_sun;
            const PS::F64 r_inv = 1.0 / sqrt(rij * rij);
            const PS::F64 r2_inv = r_inv * r_inv;
            const PS::F64 r3_inv = r2_inv * r_inv;
            const PS::F64 m_r3 = mass_sun * r3_inv;
            const PS::F64vec F0 = -m_r3*rij;
            const PS::F64vec F1 = -m_r3*vij - 3.0*rij*vij*r2_inv*F0;
            force.acc0 += F0;
            force.acc1 += F1;
        }
    }
    
    
    void setInitDt(const PS::F64 eta){
        const PS::S32 ni = n_active_;
        const PS::F64 dt_limit_new = CalcDtLimit(time_sys_, time_sync_, dt_limit_);
#pragma omp parallel for
        for(PS::S32 ip=0; ip<ni; ip++){
            const PS::S32 adr = adr_sorted_[ip];
            ptcl_multi_glb_[adr].setDt2nd(ptcl_force_[adr], eta, dt_limit_new, a0_offset_sq_);
        }
    }

    void copyForceToPtcl(){
        const PS::S32 ni = n_active_;
#pragma omp parallel for
        for(PS::S32 ip=0; ip<ni; ip++){
            const PS::S32 adr = adr_sorted_[ip];
            ptcl_multi_glb_[adr].acc0 = ptcl_force_[adr].acc0;
            ptcl_multi_glb_[adr].acc1 = ptcl_force_[adr].acc1;
            ptcl_multi_glb_[adr].acc0_pla = ptcl_force_[adr].acc0_pla;
            ptcl_multi_glb_[adr].acc1_pla = ptcl_force_[adr].acc1_pla;
#ifdef FORDEBUG
            ptcl_multi_glb_[adr].pot = ptcl_force_[adr].pot;
#endif
        }
    }

    void sortAndSelectIp(){
        const PS::S32 ni_old = n_active_;
#pragma omp parallel for
        for(PS::S32 ip=0; ip<ni_old; ip++){
            const PS::S32 adr = adr_sorted_[ip];
            time_next_[adr] += ptcl_multi_glb_[adr].dt;
        }
        SortAdr::time = &time_next_[0];
        if(merge_state_){
            PS::S32 n_tot = ptcl_multi_glb_.size();
            std::sort(&adr_sorted_[0], &adr_sorted_[0]+n_tot, SortAdr());
        }
        else{
            std::sort(&adr_sorted_[0], &adr_sorted_[0]+ni_old, SortAdr());
        }
        const PS::F64 t_n = time_next_[adr_sorted_[0]];
        const PS::S32 n = ptcl_multi_glb_.size();
        for(n_active_=1; n_active_<n; n_active_++){
            if(t_n < time_next_[adr_sorted_[n_active_]]) {
                break;
            }
        }
    }

#if 1
    void predictAll(){
        static const PS::F64 inv3 = 1.0 / 3.0;
        const PS::S32 n = ptcl_multi_glb_.size();
        const PS::F64 time_next = time_next_[adr_sorted_[0]];
#pragma omp parallel for
        for(PS::S32 i=0; i<n/2; i++){
            const PTCLHard & p_0 = ptcl_multi_glb_[i*2];
            const PTCLHard & p_1 = ptcl_multi_glb_[i*2+1];
            const PS::F64 dt_0 = time_next - p_0.time;
            const PS::F64 dt_1 = time_next - p_1.time;
            ptcl_pred_[i*2].pos = p_0.pos + dt_0*(p_0.vel  + 0.5*dt_0*(p_0.acc0 + inv3*dt_0*p_0.acc1));
            ptcl_pred_[i*2].vel = p_0.vel + dt_0*(p_0.acc0 + 0.5*dt_0*p_0.acc1);
            ptcl_pred_[i*2+1].pos = p_1.pos + dt_1*(p_1.vel  + 0.5*dt_1*(p_1.acc0 + inv3*dt_1*p_1.acc1));
            ptcl_pred_[i*2+1].vel = p_1.vel + dt_1*(p_1.acc0 + 0.5*dt_1*p_1.acc1);
#if 0
            const PS::F64 c1_0 = 0.5 * dt_0 * dt_0;
            const PS::F64 c2_0 = c1_0 * inv3 * dt_0;
            const PS::F64 c1_1 = 0.5 * dt_1 * dt_1;
            const PS::F64 c2_1 = c1_1 * inv3 * dt_1;
            ptcl_pred_[i*2].pos = p_0.pos + dt_0*p_0.vel  + c1_0*p_0.acc0 + c2_0*p_0.acc1;
            ptcl_pred_[i*2].vel = p_0.vel + dt_0*p_0.acc0 + c1_0*p_0.acc1;
            ptcl_pred_[i*2+1].pos = p_1.pos + dt_1*p_1.vel  + c1_1*p_1.acc0 + c2_1*p_1.acc1;
            ptcl_pred_[i*2+1].vel = p_1.vel + dt_1*p_1.acc0 + c1_1*p_1.acc1;
#endif
        }
        if(n%2 != 0){
            const PTCLHard & p_0 = ptcl_multi_glb_[n-1];
            const PS::F64 dt_0 = time_next - p_0.time;
            ptcl_pred_[n-1].pos = p_0.pos + dt_0*(p_0.vel  + 0.5*dt_0*(p_0.acc0 + inv3*dt_0*p_0.acc1));
            ptcl_pred_[n-1].vel = p_0.vel + dt_0*(p_0.acc0 + 0.5*dt_0*p_0.acc1);
        }
        /*
        //#pragma omp parallel for
        for(PS::S32 i=0; i<n; i++){
        const PTCLHard & p = ptcl_multi_glb_[i];
        const PS::F64 dt = time_next - p.time;
	    ptcl_pred_[i].pos = p.pos + dt*(p.vel  + 0.5*dt*(p.acc0 + inv3*dt*p.acc1));
        ptcl_pred_[i].vel = p.vel + dt*(p.acc0 + 0.5*dt*p.acc1);
        }
        */
    }
#else
#ifdef INTRINSIC_X86
    void predictAll(){
        asm("###predictAll");
        typedef double v2df __attribute__((vector_size(16)));
        typedef double v4df __attribute__((vector_size(32)));
        const PS::F64 tn = time_next_[adr_sorted_[0]];
        const PS::S32 n = ptcl_multi_glb_.size();
        for(PS::S32 i=0; i<n; i++){
            const PTCLHard & p_0 = ptcl_multi_glb_[i];
            const v4df pos{p_0.pos.x, p_0.pos.y, p_0.pos.z, 0.0};
            const v4df vel{p_0.vel.x, p_0.vel.y, p_0.vel.z, 0.0};
            const v4df acc0{p_0.acc0.x, p_0.acc0.y, p_0.acc0.z, 0.0};
            const v4df acc1{p_0.acc1.x, p_0.acc1.y, p_0.acc1.z, 0.0};
            const PS::F64 dt = tn - p_0.time;
            const v4df c0{dt, dt, dt, 0.0};
            const v4df c1 = c0 * (v4df){0.5, 0.5, 0.5, 0.0};
            const v4df c2 = c0 * (v4df){1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0};
            const v4df pos_pre = pos + c0*(vel  + c1*(acc0 + c2*acc1));
            const v4df vel_pre = vel + c0*(acc0 + c1*acc1);
            PS::F64vec & pos_pred = ptcl_pred_[i].pos;
            PS::F64vec & vel_pred = ptcl_pred_[i].vel;
            const v2df pos_pre_l = __builtin_ia32_vextractf128_pd256(pos_pre, 0);
            const v2df pos_pre_h = __builtin_ia32_vextractf128_pd256(pos_pre, 1);
            const v2df vel_pre_l = __builtin_ia32_vextractf128_pd256(vel_pre, 0);
            const v2df vel_pre_h = __builtin_ia32_vextractf128_pd256(vel_pre, 1);
            pos_pred.x = __builtin_ia32_vec_ext_v2df( pos_pre_l, 0 );
            pos_pred.y = __builtin_ia32_vec_ext_v2df( pos_pre_l, 1 );
            pos_pred.z = __builtin_ia32_vec_ext_v2df( pos_pre_h, 0 );
            vel_pred.x = __builtin_ia32_vec_ext_v2df( vel_pre_l, 0 );
            vel_pred.y = __builtin_ia32_vec_ext_v2df( vel_pre_l, 1 );
            vel_pred.z = __builtin_ia32_vec_ext_v2df( vel_pre_h, 0 );
        }
    }
#endif

#endif

    std::vector<MergeLog> merge_history_;
    PS::F64 eng_disp_;
    HardSystem(){
        merge_history_.clear();
        merge_history_.reserve(1000);
        eng_disp_ = 0.0;
    }
    
    void correctIp(const PS::F64 eta,
                   const PS::F64 dt_limit_org,
                   const PS::F64 r_out,
                   const PS::F64 r_in){
        const PS::S32 ni = n_active_;
        time_sys_ = time_next_[adr_sorted_[0]];
        if(time_sys_ == time_sync_) time_sync_ += dt_limit_org;
        const PS::F64 dt_limit_new = CalcDtLimit(time_sys_, time_sync_, dt_limit_org);
#pragma omp parallel for
        for(PS::S32 ip=0; ip<ni; ip++){
            const PS::S32 adr = adr_sorted_[ip];
            ptcl_multi_glb_[adr].correct(ptcl_force_[adr], eta, dt_limit_new, a0_offset_sq_);
        }
#ifdef MERGE
        merge_state_ = false;
        const PS::S32 n_int = pair_merge_adr_.size();
        for(PS::S32 i=0; i<n_int; i++){
            const PS::S32 adr_i = pair_merge_adr_[i].first;
            const PS::S32 adr_j = pair_merge_adr_[i].second;
            if(merge_flag_glb_[adr_i] == true || merge_flag_glb_[adr_j] == true) continue;
            //PTCLHard & ptcl_i = ptcl_multi_glb_[adr_i];
            //PTCLHard & ptcl_j = ptcl_multi_glb_[adr_j];
            //merge2body(ptcl_i, ptcl_j, merge_history_);
            //merge2body(ptcl_i, ptcl_j, eng_disp_, merge_history_);
            merge2body(ptcl_multi_glb_, pair_adr_multi_glb_, adr_i, adr_j, eng_disp_, merge_history_);

            PS::S64 adr_dead = (ptcl_multi_glb_[adr_i].mass == 0.0) ? adr_i : adr_j;
            merge_flag_glb_[adr_dead] = true;
            merge_state_ = true;
        }

        if(merge_state_){
            const PS::F64 eps_sq = EPISoft::eps * EPISoft::eps;
            const PS::S32 ni = n_active_;
            for(PS::S32 ip=0; ip<ni; ip++){
                const PS::S32 adr_i = adr_sorted_[ip];
                if(ptcl_multi_glb_[adr_i].isDead()){
                    ptcl_multi_glb_[adr_i].dt = PS::LARGE_FLOAT;
                }
                else{
                    ptcl_force_[adr_i].clear();
                    const PS::S32 nj = ptcl_multi_glb_[adr_i].n_ngb;
                    const PS::S32 adr_pair = ptcl_multi_glb_[adr_i].adr_pair;
                    for(PS::S32 jp=0; jp<nj; jp++){
                        const PS::S32 adr_j = pair_adr_multi_glb_[adr_pair+jp].second;
                        if(ptcl_multi_glb_[adr_j].isDead()) continue;
                        CalcAcc0AndAcc1Cutoff(ptcl_multi_glb_[adr_i].pos,    ptcl_multi_glb_[adr_i].vel,
                                              ptcl_force_[adr_i].acc0, ptcl_force_[adr_i].acc1,
                                              ptcl_multi_glb_[adr_j].pos,   ptcl_multi_glb_[adr_j].vel,
                                              ptcl_multi_glb_[adr_j].mass,   eps_sq,
                                              r_out, r_in);
                    }
                    calcAcc0AndAcc1FromSun(ptcl_multi_glb_[adr_i].pos, 
                                           ptcl_multi_glb_[adr_i].vel,
                                           ptcl_force_[adr_i].acc0,      
                                           ptcl_force_[adr_i].acc1);
                    ptcl_multi_glb_[adr_i].acc0 = ptcl_force_[adr_i].acc0;
                    ptcl_multi_glb_[adr_i].acc1 = ptcl_force_[adr_i].acc1;
                    ptcl_multi_glb_[adr_i].setDt2nd(ptcl_force_[adr_i], 0.01*eta, dt_limit_new, a0_offset_sq_);
                }
            }
        }
#endif
    }

    bool evolve(const PS::F64 r_out, const PS::F64 r_in,
                const PS::F64 eta,   const PS::F64 dt_limit_org){
        PS::F64 offset = PS::GetWtime();
        predictAll();
        wtime_profile_.predict_ += PS::GetWtime() - offset;

        offset = PS::GetWtime();
        calcAcc0AndAcc1(r_out, r_in);
        wtime_profile_.force_ += PS::GetWtime() - offset;

        offset = PS::GetWtime();
        correctIp(eta, dt_limit_org, r_out, r_in);
        wtime_profile_.correct_ += PS::GetWtime() - offset;
	
#ifdef FORDEBUG	
        if(n_active_ == ptcl_multi_glb_.size()){
            calc_eng();
            std::cout<<"time_sys_="<<time_sys_
		     <<" time_end_="<<time_end_
		     <<" loop="<< loop
		     <<" eng_tot_init_="<<eng_tot_init_
		     <<" eng_tot_="<<eng_tot_
		     <<" (eng_tot_init_ - eng_tot_)/eng_tot_init_="
                     <<(eng_tot_init_ - eng_tot_)/eng_tot_init_<<std::endl;
            std::cout<<"ptcl_multi_glb_[0].pos="<<ptcl_multi_glb_[0].pos<<std::endl;
            std::cout<<std::endl;
        }
#endif
        offset = PS::GetWtime();
        sortAndSelectIp();
        wtime_profile_.sort_ += PS::GetWtime() - offset;
	
        n_loop_++;
        if(time_sys_ == time_end_) return true;
        else return false;
    }

    bool evolveKepler(const PS::F64 r_out,        const PS::F64 r_in,
		      const PS::F64 eta,          const PS::F64 dt_limit_org,
		      const PS::F64 mass_sun=1.0, const PS::F64vec pos_sun=0.0,
		      const PS::F64vec vel_sun=0.0){
        PS::F64 offset = PS::GetWtime();
        predictAll();
        wtime_profile_.predict_ += PS::GetWtime() - offset;

        offset = PS::GetWtime();
        calcAcc0AndAcc1(r_out, r_in);
        wtime_profile_.force_ += PS::GetWtime() - offset;

        offset = PS::GetWtime();
        calcAcc0AndAcc1Kepler(mass_sun, pos_sun, vel_sun);
        wtime_profile_.force_kepler_ += PS::GetWtime() - offset;

	reduce_force();

        offset = PS::GetWtime();
        correctIp(eta, dt_limit_org, r_out, r_in);
        wtime_profile_.correct_ += PS::GetWtime() - offset;

        offset = PS::GetWtime();
        sortAndSelectIp();
        wtime_profile_.sort_ += PS::GetWtime() - offset;

        n_loop_++;
        if(time_sys_ == time_end_) return true;
        else return false;
    }
    
    void scatterData(){
        const PS::S32 n_loc = ptcl_multi_loc_.size();
        PS::Comm::scatterV(&ptcl_multi_glb_[0], &n_ptcl_array_[0], &n_ptcl_disp_[0],
                           &ptcl_multi_loc_[0], n_loc);

	
#if 0
        if(PS::Comm::getRank() == 1){
            for(PS::S32 i=0; i<n_loc; i++){
                ptcl_hard_loc[i].dump();
            }
        }
#endif
    }

    template<class Tsoft>
    void copyPtclHardToSoft(Tsoft & system_soft){
        const PS::S32 n = ptcl_multi_loc_.size();
#pragma omp parallel for
        for(PS::S32 i=0; i<n; i++){
            const PS::S32 adr = adr_ptcl_multi_loc_[i];
            assert(system_soft[adr].id == ptcl_multi_loc_[i].id);
            system_soft[adr].mass = ptcl_multi_loc_[i].mass;
            system_soft[adr].pos = ptcl_multi_loc_[i].pos;
            system_soft[adr].vel = ptcl_multi_loc_[i].vel;
        }
    }

    PS::F64 eng_tot_init_;
    PS::F64 eng_kin_init_;
    PS::F64 eng_pot_init_;

    PS::F64 eng_tot_;
    PS::F64 eng_kin_;
    PS::F64 eng_pot_;

    static void calc_eng_impl(const std::vector<PTCLHard> ptcl, PS::F64 & et, PS::F64 & ek, PS::F64 & ep){
        const PS::S32 n = ptcl.size();
        et = ek = ep = 0.0;
        for(PS::S32 i=0; i<n; i++){
            const PTCLHard & p = ptcl[i];
            ek += 0.5 * p.mass * p.vel * p.vel;
            ep += 0.5 * p.mass * p.pot;
        }
        et = ek + ep;
    }

    void calc_eng_init(){
        calc_eng_impl(ptcl_multi_glb_, eng_tot_init_, eng_kin_init_, eng_pot_init_);
    }

    void calc_eng(){
        calc_eng_impl(ptcl_multi_glb_, eng_tot_, eng_kin_, eng_pot_);
    }

    std::vector< std::pair<PS::S64, PS::S64> > pair_id_2body_loc_;
    std::vector< std::pair<PS::S64, PS::S64> > pair_id_multi_loc_;
    std::vector< std::pair<PS::S64, PS::S64> > pair_adr_2body_loc_; // point to ptcl_loc_
    std::vector<PTCLHard> ptcl_2body_loc_;
    std::vector<PTCLHard> ptcl_multi_loc_;
    std::vector<PS::S64> adr_ptcl_2body_loc_;
    std::vector<PS::S64> adr_ptcl_multi_loc_;
    void selectIsolatedSystem(){
        std::unordered_map<PS::S64, PS::S64> idx_to_adr_2body;
        makeDictionaryLocal();
        pair_id_multi_loc_.clear();
        pair_id_2body_loc_.clear();
        ptcl_2body_loc_.clear();
        ptcl_multi_loc_.clear();
        idx_to_adr_2body.clear();
        adr_ptcl_2body_loc_.clear();
        adr_ptcl_multi_loc_.clear();
        PS::S64 adr_i_old = -1;
        //const size_t n_pair = pair_adr_loc_.size();
        const size_t n_pair = pair_id_loc_.size();
        //std::cerr<<"n_pair="<<n_pair<<std::endl;
        PS::S64 n_cnt_2body = 0;
        for(size_t i=0; i<n_pair; i++){
            //const PS::S64 id_i = pair_id_loc_[i].first;
            //const PS::S64 id_j = pair_id_loc_[i].second;
            const PS::S64 adr_i = pair_adr_loc_[i].first;
            const PS::S64 adr_j = pair_adr_loc_[i].second;
            if( ptcl_loc_[adr_i].n_ngb == 1 && adr_j >= 0 && ptcl_loc_[adr_j].n_ngb == 1){
                pair_id_2body_loc_.push_back( pair_id_loc_[i] );
                ptcl_2body_loc_.push_back( ptcl_loc_[adr_i] );
                adr_ptcl_2body_loc_.push_back( adr_ptcl_loc_[adr_i] );
                idx_to_adr_2body.insert(std::pair<PS::S64, PS::S64>(ptcl_loc_[adr_i].id, n_cnt_2body));
                n_cnt_2body++;
            }
            else{
                pair_id_multi_loc_.push_back(pair_id_loc_[i]);
                if(adr_i != adr_i_old){
                    ptcl_multi_loc_.push_back( ptcl_loc_[adr_i] );
                    adr_ptcl_multi_loc_.push_back( adr_ptcl_loc_[adr_i] );
                    adr_i_old = adr_i;
                }
            }
        }
	if(pair_id_loc_.size() != pair_id_multi_loc_.size() + pair_id_2body_loc_.size()){
	    std::cerr<<"pair_id_loc_.size()="<<pair_id_loc_.size()
		     <<" pair_id_multi_loc_.size()="<<pair_id_multi_loc_.size()
		     <<" pair_id_2body_loc_.size()="<<pair_id_2body_loc_.size()<<std::endl;
	}
	if( ptcl_loc_.size() != ptcl_multi_loc_.size() + ptcl_2body_loc_.size() ){
	    std::cerr<<"ptcl_loc_.size()="<<ptcl_loc_.size()
		     <<" ptcl_multi_loc_.size()="<<ptcl_multi_loc_.size()
		     <<" ptcl_2body_loc_.size()="<<ptcl_2body_loc_.size()<<std::endl;
	}
        assert( pair_id_loc_.size() == pair_id_multi_loc_.size() + pair_id_2body_loc_.size() );
        assert( ptcl_loc_.size() == ptcl_multi_loc_.size() + ptcl_2body_loc_.size() );
        pair_adr_2body_loc_.clear();
        const size_t n_pair_2body = pair_id_2body_loc_.size();
        for(size_t i=0; i<n_pair_2body; i++){
            pair_adr_2body_loc_.push_back(std::pair<PS::S64, PS::S64>
                                          (idx_to_adr_2body[pair_id_2body_loc_[i].first],
                                           idx_to_adr_2body[pair_id_2body_loc_[i].second]));
        }
    }
    

    void predict2body(PTCLPred & pred_i,       PTCLPred & pred_j,
                      const PTCLHard & ptcl_i, const PTCLHard & ptcl_j, const PS::F64 dt){
        static const PS::F64 inv3 = 1.0 / 3.0;
        pred_i.pos = ptcl_i.pos + dt*(ptcl_i.vel  + 0.5*dt*(ptcl_i.acc0 + inv3*dt*ptcl_i.acc1));
        pred_i.vel = ptcl_i.vel + dt*(ptcl_i.acc0 + 0.5*dt* ptcl_i.acc1);
        pred_j.pos = ptcl_j.pos + dt*(ptcl_j.vel  + 0.5*dt*(ptcl_j.acc0 + inv3*dt*ptcl_j.acc1));
        pred_j.vel = ptcl_j.vel + dt*(ptcl_j.acc0 + 0.5*dt* ptcl_j.acc1);
    }

    void calcEng2body(const PTCLHard & pi,    const PTCLHard & pj,
                      const PS::F64 mass_sun, const PS::F64vec & pos_sun, const PS::F64vec & pos_vel,
                      const PS::F64 eps_sq,
                      PS::F64 & Ekin, PS::F64 & Epot, PS::F64 & Etot){
        Ekin = Epot = Etot = 0.0;
        Ekin = 0.5 * pi.mass * pi.vel * pi.vel
            + 0.5 * pj.mass * pj.vel * pj.vel;
        PS::F64vec rij = pi.pos - pj.pos;
        Epot -= pi.mass * pj.mass / sqrt(rij*rij+eps_sq);
        PS::F64vec dr = pi.pos - pos_sun;
        Epot -= pi.mass*mass_sun / sqrt(dr*dr);
        dr = pj.pos - pos_sun;
        Epot -= pj.mass*mass_sun / sqrt(dr*dr);
        Etot = Ekin + Epot;
    }
    
    std::vector<PTCLPred> ptcl_pred_2body_loc_;
    std::vector<PTCLForce> ptcl_force_2body_loc_;
    std::vector<bool> merge_flag_2body_loc_;
    //std::vector< std::pair<PS::S64, PS::S64> > pair_isolated_merge_adr_;
    template<class Tsoft>
    void evolveIsolatedSystem(Tsoft & system,
                              const PS::F64 r_out, const PS::F64 r_in,
                              const PS::F64 mass_sun, const PS::F64vec pos_sun,
                              const PS::F64vec vel_sun,
                              const PS::F64 eta_s, const PS::F64 eta,
                              const PS::F64 time_end, const PS::F64 dt_limit_org){
        const PS::S32 n_thread_max = PS::Comm::getNumberOfThread();
        static std::vector<MergeLog> * merge_log;
        static PS::F64 * eng_disp_omp;
        static bool first = true;
        if(first){
            merge_log = new std::vector<MergeLog>[n_thread_max];
            eng_disp_omp = new PS::F64[n_thread_max];
            first = false;
        }
        for(PS::S32 i=0; i<n_thread_max; i++){
            merge_log[i].clear();
            eng_disp_omp[i] = 0.0;
        }
        const size_t n_2body = ptcl_2body_loc_.size();
        //std::cerr<<"PS::Comm::getRank()="<<PS::Comm::getRank()<<" n="<<n<<std::endl;
        ptcl_pred_2body_loc_.clear();
        ptcl_pred_2body_loc_.resize(n_2body);
        ptcl_force_2body_loc_.clear();
        ptcl_force_2body_loc_.resize(n_2body);
        merge_flag_2body_loc_.clear();
        merge_flag_2body_loc_.resize(n_2body);
        //pair_isolated_merge_adr_.clear();
#pragma omp parallel for
        for(size_t ip=0; ip<n_2body; ip++){
            ptcl_pred_2body_loc_[ip].pos = ptcl_2body_loc_[ip].pos;
            ptcl_pred_2body_loc_[ip].vel = ptcl_2body_loc_[ip].vel;
            ptcl_2body_loc_[ip].time = 0.0;
            ptcl_2body_loc_[ip].setRMerge();
            merge_flag_2body_loc_[ip] = false;
        }
        const PS::S32 n_pair = pair_adr_2body_loc_.size();
        const PS::F64 eps_sq = EPISoft::eps * EPISoft::eps;
#pragma omp parallel for
        for(PS::S32 np=0; np<n_pair; np++){
            PS::S32 ith = PS::Comm::getThreadNum();
            PS::F64 time_sys_tmp = 0.0;
            PS::F64 time_sync_tmp = time_sys_tmp + dt_limit_org;
            PS::F64 time_end_tmp = time_end;

            const PS::S64 adr_i = pair_adr_2body_loc_[np].first;
            const PS::S64 adr_j = pair_adr_2body_loc_[np].second;
            if(adr_i > adr_j) continue;

            PTCLHard & ptcl_i = ptcl_2body_loc_[adr_i];
            PTCLHard & ptcl_j = ptcl_2body_loc_[adr_j];
            PTCLForce & force_i = ptcl_force_2body_loc_[adr_i];
            PTCLForce & force_j = ptcl_force_2body_loc_[adr_j];
            PTCLPred & pred_i = ptcl_pred_2body_loc_[adr_i];
            PTCLPred & pred_j = ptcl_pred_2body_loc_[adr_j];
/*
            PS::F64 Ekin0, Epot0, Etot0, Ekin1, Epot1, Etot1;
            calcEng2body(ptcl_i, ptcl_j,
                         mass_sun, pos_sun, vel_sun, eps_sq,
                         Ekin0, Epot0, Etot0);
*/
            PS::F64 r2 = 0.0;
            force_i.clear();
            force_j.clear();
	    /*
            CalcAcc0Acc1AndR2CutoffPair(ptcl_i.mass, pred_i.pos, pred_i.vel, force_i.acc0, force_i.acc1,
                                        ptcl_j.mass, pred_j.pos, pred_j.vel, force_j.acc0, force_j.acc1,
                                        r2,   eps_sq, r_out, r_in);
	    */
            CalcAcc0Acc1AndR2CutoffPair(ptcl_i.mass, pred_i.pos, pred_i.vel, force_i.acc0_pla, force_i.acc1_pla,
                                        ptcl_j.mass, pred_j.pos, pred_j.vel, force_j.acc0_pla, force_j.acc1_pla,
                                        r2,   eps_sq, r_out, r_in);
	    /*
            CalcAcc0Acc1AndR2(pred_i.pos,    pred_i.vel,
                              force_i.acc0,  force_i.acc1,
                              pos_sun, vel_sun, mass_sun, eps_sq);
	    
            CalcAcc0Acc1AndR2(pred_j.pos,     pred_j.vel,
                              force_j.acc0,   force_j.acc1,
                              pos_sun, vel_sun, mass_sun, eps_sq);
	    */

            CalcAcc0Acc1AndR2(pred_i.pos,    pred_i.vel,
                              force_i.acc0,  force_i.acc1,
                              pos_sun, vel_sun, mass_sun, 0.0);
	    
            CalcAcc0Acc1AndR2(pred_j.pos,     pred_j.vel,
                              force_j.acc0,   force_j.acc1,
                              pos_sun, vel_sun, mass_sun, 0.0);
	    force_i.reduce();
	    force_j.reduce();
            PS::F64 dt_limit_new = CalcDtLimit(time_sys_tmp, time_sync_tmp, dt_limit_org);
            ptcl_i.setDt2nd(force_i, eta_s, dt_limit_new, a0_offset_sq_);
            ptcl_j.setDt2nd(force_j, eta_s, dt_limit_new, a0_offset_sq_);
            PS::F64 dt = (ptcl_i.dt < ptcl_j.dt) ? ptcl_i.dt : ptcl_j.dt;
            ptcl_i.dt = ptcl_j.dt = dt;
            ptcl_i.acc0 = force_i.acc0;
            ptcl_i.acc1 = force_i.acc1;
            ptcl_j.acc0 = force_j.acc0;
            ptcl_j.acc1 = force_j.acc1;
            ptcl_i.acc0_pla = force_i.acc0_pla;
            ptcl_i.acc1_pla = force_i.acc1_pla;
            ptcl_j.acc0_pla = force_j.acc0_pla;
            ptcl_j.acc1_pla = force_j.acc1_pla;
            while(time_sys_tmp != time_end_tmp){
                predict2body(pred_i, pred_j,
                             ptcl_i, ptcl_j,
                             dt);
                force_i.clear();
                force_j.clear();
		/*
                CalcAcc0Acc1AndR2CutoffPair(ptcl_i.mass, pred_i.pos, pred_i.vel, force_i.acc0, force_i.acc1,
                                            ptcl_j.mass, pred_j.pos, pred_j.vel, force_j.acc0, force_j.acc1,
                                            r2,   eps_sq, r_out, r_in);
		*/
                CalcAcc0Acc1AndR2CutoffPair(ptcl_i.mass, pred_i.pos, pred_i.vel, force_i.acc0_pla, force_i.acc1_pla,
                                            ptcl_j.mass, pred_j.pos, pred_j.vel, force_j.acc0_pla, force_j.acc1_pla,
                                            r2,   eps_sq, r_out, r_in);
		/*
                CalcAcc0Acc1AndR2(pred_i.pos,	  pred_i.vel,
                                  force_i.acc0,   force_i.acc1,
                                  pos_sun, vel_sun, mass_sun, eps_sq);
                CalcAcc0Acc1AndR2(pred_j.pos,     pred_j.vel,
                                  force_j.acc0,   force_j.acc1,
                                  pos_sun, vel_sun, mass_sun, eps_sq);
		*/
		/*
                CalcAcc0Acc1AndR2(pred_i.pos,	  pred_i.vel,
                                  force_i.acc0,   force_i.acc1,
                                  pos_sun, vel_sun, mass_sun, 0.0);
                CalcAcc0Acc1AndR2(pred_j.pos,     pred_j.vel,
                                  force_j.acc0,   force_j.acc1,
                                  pos_sun, vel_sun, mass_sun, 0.0);
		*/
                CalcAcc0Acc1AndR2(pred_i.pos,	  pred_i.vel,
                                  force_i.acc0,   force_i.acc1,
                                  pos_sun, vel_sun, mass_sun, 0.0);
                CalcAcc0Acc1AndR2(pred_j.pos,     pred_j.vel,
                                  force_j.acc0,   force_j.acc1,
                                  pos_sun, vel_sun, mass_sun, 0.0);
		force_i.reduce();
		force_j.reduce();
#ifdef MERGE
                const PS::F64 r_merge = ptcl_i.r_merge + ptcl_j.r_merge;
                const PS::F64 r_merge_2 = r_merge * r_merge;
                bool merge = false;
                if(r2 < r_merge_2){
                    //pair_isolated_merge_adr_.push_back( std::pair<PS::S64, PS::S64>(adr_i, adr_j) );
                    merge = true;
                }
#endif //MERGE
                time_sys_tmp += dt;
                if(time_sys_tmp == time_sync_tmp) time_sync_ += dt_limit_org;
                dt_limit_new = CalcDtLimit(time_sys_tmp, time_sync_tmp, dt_limit_org);
                ptcl_i.correct(force_i, eta, dt_limit_new, a0_offset_sq_);
                ptcl_j.correct(force_j, eta, dt_limit_new, a0_offset_sq_);
                dt = (ptcl_2body_loc_[adr_i].dt < ptcl_2body_loc_[adr_j].dt) ? ptcl_2body_loc_[adr_i].dt : ptcl_2body_loc_[adr_j].dt;
                ptcl_i.dt = ptcl_j.dt = dt;
#ifdef MERGE
                if(merge){
                    //merge2body(ptcl_i, ptcl_j, merge_log[ith]);
                    merge2body(ptcl_i, ptcl_j, eng_disp_omp[ith], merge_log[ith]);
                    //merge2body(ptcl_2body_loc_, adr_i, adr_j, eng_disp_omp[ith], merge_log[ith]);
                    PTCLHard & ptcl_merge = (ptcl_i.mass == 0.0) ? ptcl_j : ptcl_i;
                    const PS::F64 dt_rem = time_end_tmp - time_sys_tmp;
                    DriveKepler(mass_sun, ptcl_merge.mass,
                                pos_sun, ptcl_merge.pos,
                                vel_sun, ptcl_merge.vel, dt_rem);
                    ptcl_merge.time = time_end_tmp;
                    time_sys_tmp = time_end_tmp;
                }
#endif //MERGE
            }
/*
            calcEng2body(ptcl_i,   ptcl_j,
                         mass_sun, pos_sun, vel_sun, eps_sq,
                         Ekin1, Epot1, Etot1);
*/
        }
        for(PS::S32 i=0; i<n_thread_max; i++){
            const size_t j_max = merge_log[i].size();
            for(size_t j=0; j<j_max; j++){
                merge_history_.push_back(merge_log[i][j]);
            }
            eng_disp_ += eng_disp_omp[i];
        }
        //std::cerr<<"PS::Comm::getRank()="<<PS::Comm::getRank()<<" check7"<<std::endl;
        for(size_t i=0; i<n_2body; i++){
            const PS::S64 adr = adr_ptcl_2body_loc_[i];
            assert(ptcl_2body_loc_[i].id == system[adr].id);
            system[adr].mass = ptcl_2body_loc_[i].mass;
            system[adr].pos = ptcl_2body_loc_[i].pos;
            system[adr].vel = ptcl_2body_loc_[i].vel;
        }
    }

    void reduce_force(){
        const PS::S32 ni = n_active_;
#pragma omp parallel for
        for(PS::S32 ip=0; ip<ni; ip++){
            const PS::S32 adr = adr_sorted_[ip];
            ptcl_force_[adr].reduce();
	}
    }

    void evolve_multi_sirial(const PS::F64 t_end, 
			     const PS::F64 dt_limit_hard,
			     const PS::F64 r_out,
			     const PS::F64 r_in,
			     const PS::F64 eta_s,
			     const PS::F64 eta){
	makeDictionaryGlobal();
	setUpGlb(t_end, dt_limit_hard);
	calcAcc0AndAcc1(r_out, r_in);
	calcAcc0AndAcc1Kepler(mass_sun_, pos_sun_, vel_sun_);
	reduce_force();
	setInitDt(eta_s);
	copyForceToPtcl();
	//calc_eng_init();
	sortAndSelectIp();
	bool end_hard_part = false;
	while( !end_hard_part ){
	    end_hard_part = evolveKepler(r_out, r_in, eta, dt_limit_hard);
	}
    }

    PS::S64 n_loop_cum_;
    size_t n_ptcl_2body_loc_;
    size_t n_pair_2body_loc_;
    size_t n_ptcl_2body_glb_;
    size_t n_pair_2body_glb_;
    size_t n_ptcl_multi_;
    size_t n_pair_multi_;
    void accumulate_counter(){
	n_loop_cum_ += n_loop_;
	n_ptcl_2body_loc_ += ptcl_2body_loc_.size();
	n_pair_2body_loc_ += pair_id_2body_loc_.size();
	n_ptcl_multi_ += ptcl_multi_glb_.size();
	n_pair_multi_ += pair_id_multi_glb_.size();
    }
    void dump_counter(std::ofstream & fout){
	n_ptcl_2body_glb_ = PS::Comm::getSum(n_ptcl_2body_loc_);
	n_pair_2body_glb_ = PS::Comm::getSum(n_pair_2body_loc_);
	fout<<"n_loop_cum="<<n_loop_cum_
	    <<" n_ptcl_2body_loc_="<<n_ptcl_2body_loc_
	    <<" n_pair_2body_loc_="<<n_pair_2body_loc_
	    <<" n_ptcl_2body_glb_="<<n_ptcl_2body_glb_
	    <<" n_pair_2body_glb_="<<n_pair_2body_glb_
	    <<" n_ptcl_multi_="<<n_ptcl_multi_
	    <<" n_pair_multi_="<<n_pair_multi_<<std::endl;
    }
    void clear_counter(){
	n_loop_cum_ = 0;
	n_ptcl_2body_loc_ = n_pair_2body_loc_ 
	    = n_ptcl_2body_glb_ = n_pair_2body_glb_
	    = n_ptcl_multi_ =  n_pair_multi_ = 0;
    }
};
