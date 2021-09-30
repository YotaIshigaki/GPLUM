//#define DEBUG_2017_10_30

#define PARTICLE_SIMULATOR_CHECK_SEARCH_MODE

#include<iostream>
#include<fstream>
#include<unistd.h>
#include<random>
#include<particle_simulator.hpp>
#include<map>

PS::F64 dEcrit = 1e-2;
//PS::EXCHANGE_LET_MODE EX_LET_MODE = PS::EXCHANGE_LET_A2A;
//PS::EXCHANGE_LET_MODE EX_LET_MODE = PS::EXCHANGE_LET_P2P_EXACT;
PS::EXCHANGE_LET_MODE EX_LET_MODE = PS::EXCHANGE_LET_P2P_FAST;

#define GENERATE_ENUM_ITERATOR(T) \
    inline T operator++(T& x) { return x = (T)(std::underlying_type<T>::type(x) + 1); } \
    inline T operator*(T c) { return c; } \
    inline T begin(T r) { return static_cast<T>(0); }

typedef PS::F64 REAL;
typedef PS::F64vec REALVEC;
typedef PS::F64ort REALORT;

template<class Tsys>
void SetRsearchRand(Tsys & sys){
    PS::S32 n_loc = sys.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_loc; i++){
        sys[i].r_search = sys[i].r_search_1;
    }
}

template<class Tsys>
void SetRsearchAve(Tsys & sys){
    PS::S32 n_loc = sys.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_loc; i++){
        sys[i].r_search = sys[i].r_search_0;
    }
}

enum class IcMode{
    uni_qube_grid,
    n_tot,
    uni_qube_rand,
    uni_qube_ill,
};

GENERATE_ENUM_ITERATOR(IcMode);
IcMode end(IcMode) { return IcMode::n_tot; }
//IcMode end(IcMode) { return IcMode::uni_qube_grid; }

typedef int BcMode;

typedef int FrMode;

typedef int IlMode;
const IlMode n_il_mode = 3;
IlMode il_mode_ar[n_il_mode];

typedef int IdMode;

std::map<int, std::string> MapIdMode = {
    {0, "multiwalk:no,  index: no"}, 
    {1, "multiwalk:yes, index: no"}, 
    {2, "multiwalk:yes, index: yes"}, 
};

std::map<int, std::string> MapSearchMode = {
    {0, "LONG_NO_CUTOFF"},
    {1, "LONG_CUTOFF"},
    {2, "LONG_SCATTER"},
    {3, "LONG_CUTOFF_SCATTER"},
    {4, "LONG_SYMMETRY"},
    {5, "SHORT_GATHER"},
    {6, "SHORT_SCATTER"},
    {7, "SHORT_SYMMETRY"},
};

std::map<int, std::string> MapInteractionListMode = {
    {0, "MAKE_LIST"},
    {1, "MAKE_LIST_FOR_REUSE"},
    {2, "REUSE_LIST"},
};

struct TestResult{
    bool pass;
    IcMode ic_mode;
    BcMode bc_mode;
    FrMode fr_mode;
    enum PS::INTERACTION_LIST_MODE il_mode;
    IdMode id_mode;
    PS::S64 n_err_ngb;
    PS::S64 n_err_mass;
    PS::S64 n_err_pot;
    TestResult(const TestResult & x) : pass(x.pass), ic_mode(x.ic_mode), bc_mode(x.bc_mode),
                                       fr_mode(x.fr_mode), il_mode(x.il_mode), id_mode(x.id_mode), n_err_ngb(x.n_err_ngb), n_err_mass(x.n_err_mass), n_err_pot(x.n_err_pot){}

    TestResult(const bool _pass, const IcMode _ic_mode, const BcMode _bc_mode,
               const FrMode _fr_mode, const PS::INTERACTION_LIST_MODE _il_mode,
               const IdMode _id_mode, const PS::S64 _n_err_ngb, const PS::S64 _n_err_mass, const PS::S64 _n_err_pot=0)
        : pass(_pass), ic_mode(_ic_mode), bc_mode(_bc_mode),
          fr_mode(_fr_mode), il_mode(_il_mode), id_mode(_id_mode), n_err_ngb(_n_err_ngb), n_err_mass(_n_err_mass), n_err_pot(_n_err_pot){}
    
    void dump(const std::ostream & fout=std::cout){
        if(pass) std::cerr<<"pass:";
        else{
            if( (fr_mode == static_cast<int>(PS::LONG_NO_CUTOFF)
                 || fr_mode == static_cast<int>(PS::LONG_CUTOFF) )
                && n_err_ngb > 0){
                std::cerr<<"not sure: ";
            }
            else{ std::cerr<<"fail:";}
        }
        std::cerr<<" ic_mode= "<<static_cast<int>(ic_mode)
                 <<" bc_mode= "<<static_cast<int>(bc_mode)
                 <<" fr_mode= "<<MapSearchMode[static_cast<int>(fr_mode)]
                 <<" il_mode= "<<MapInteractionListMode[static_cast<int>(il_mode)]
                 <<" id_mode= "<<MapIdMode[static_cast<int>(id_mode)]
                 <<" n_err_ngb= "<<n_err_ngb
                 <<" n_err_mass= "<<n_err_mass
                 <<" n_err_pot= "<<n_err_pot
                 <<std::endl;
    }
};

struct Force{
    PS::S32 n_ngb;
    PS::S32 id;
    REAL mass;
    REAL pot;
    void clear() {
        id = -1;
        n_ngb = 0.0;
        mass = 0.0;
        pot = 0.0;
    }
};

class FP{
public:
    PS::S64 id;
    REAL mass;
    REAL mass_tmp;
    REALVEC pos;
    REALVEC vel;
    PS::S32 n_ngb;
    REAL pot;
    REAL r_search;
    REAL r_search_0;
    REAL r_search_1;
    REALVEC getPos() const { return pos; }
    void setPos(const REALVEC & p) { pos = p; }
    void copyFromForce(const Force & force){
        n_ngb = force.n_ngb;
        mass_tmp = force.mass;
        pot = force.pot;
    }
    REAL getRSearch() const {
        return r_search;
    }
};

struct EpSimple{
    PS::S32 id = 0;
    REAL mass = 0.0;
    REALVEC pos = 0.0;
    REALVEC getPos() const {
        return pos;
    }
    REAL getCharge() const {
	return mass;
    }
    void setPos(const REALVEC & _pos){
	pos = _pos;
    }
    void copyFromFP(const FP & fp){
        id   = fp.id;
        pos  = fp.pos;
        mass = fp.mass;
    }
    template<class Tep>
    void copyFromEpInFdps(const Tep & ep){
        id       = ep.id;
        mass     = ep.mass;
        pos      = ep.pos;
    }
    //REAL getRSearch(){ return 0.0; }
};

struct EpAll{
    PS::S32 id = 0;
    REAL mass = 0.0;
    REALVEC pos = 0.0;
    REAL r_search = 0.0;
    REALVEC getPos() const {
        return pos;
    }
    REAL getCharge() const {
	return mass;
    }
    REAL getRSearch() const {
        return r_search;
    }
    void setPos(const REALVEC & _pos){
	pos = _pos;
    }
    void copyFromFP(const FP & fp){
        id   = fp.id;
        pos  = fp.pos;
        mass = fp.mass;
	r_search  = fp.r_search;
    }
    template<class Tep>
    void copyFromEpInFdps(const Tep & ep){
        id       = ep.id;
        mass     = ep.mass;
        pos      = ep.pos;
        r_search = ep.r_search;
    }
};
template<>
void EpAll::copyFromEpInFdps<EpSimple>(const EpSimple & ep){
    id       = ep.id;
    mass     = ep.mass;
    pos      = ep.pos;
}    


class MomentAll{
    public:
    REAL mass;
    REALVEC pos;
    REALORT vertex_out_;
    REALORT vertex_in_;
    MomentAll(){
        mass = 0.0;
        pos = 0.0;
        vertex_out_.init();
        vertex_in_.init();
    }
    MomentAll(const REAL m, const REALVEC & p){ 
        mass = m;
        pos = p;
    }
    REALORT getVertexOut() const { return vertex_out_; }
    REALORT getVertexIn() const { return vertex_in_; }
    void init(){
        mass = 0.0;
        pos = 0.0;
        vertex_out_.init();
        vertex_in_.init();
    }
    REALVEC getPos() const {
        return pos;
    }
    REAL getCharge() const {
        return mass;
    }
    template<class Tepj>
    void accumulateAtLeaf(const Tepj & epj){
        this->mass += epj.getCharge();
        this->pos += epj.getCharge() * epj.getPos();
        (this->vertex_out_).merge(epj.getPos(), epj.getRSearch());
        (this->vertex_in_).merge(epj.getPos());
    }
    template<class Tepj>
    void accumulateAtLeaf2(const Tepj & epj){}
    void set(){
        pos = pos / mass;
    }
    void accumulate(const MomentAll & mom){
        this->mass += mom.mass;
        this->pos += mom.mass * mom.pos;
        (this->vertex_out_).merge(mom.vertex_out_);
        (this->vertex_in_).merge(mom.vertex_in_);
    }
    void accumulate2(const MomentAll & mom){}
};

struct SpAll{
    int id;
    REAL mass;
    REALVEC pos;
    template<class Tmom>
    void copyFromMoment(const Tmom & mom){
        REAL mass = mom.mass;
        REALVEC pos = mom.pos;
        this->mass = mass;
        this->pos = pos;
        id = -1;
    }
    void clear(){
        mass = 0.0;
        pos = 0.0;
        id = -1;
    }
    REAL getCharge() const {
        return mass;
    }
    REALVEC getPos() const {
        return pos;
    }
    void setPos(const REALVEC & pos_new) {
        pos = pos_new;
    }
    MomentAll convertToMoment() const {
        return MomentAll(mass, pos);
    }

    template<class Tsp>
    void copyFromSpInFdps(const Tsp & sp){
        pos = sp.pos;
        mass = sp.mass;
    }
    
};

template<int id>
struct GetTreeType{ 
    typedef PS::TreeForForceLong<Force, EpAll, EpAll>::MonopoleWithSymmetrySearch 
    TreeType;
};

template<>
struct GetTreeType<0>{
    typedef PS::TreeForForceShort<Force, EpAll, EpAll>::Gather TreeType;
};


const int N_EPI_LIMIT = 100000;
const int N_EPJ_LIMIT = 1000000;
const int N_SPJ_LIMIT = 1000000;

static EpAll dev_epj[N_EPJ_LIMIT];
static SpAll dev_spj[N_SPJ_LIMIT];
static Force dev_force[N_EPI_LIMIT];

enum class SEARCH_R{
    EPI_R,
    EPJ_R,
    MAX_R,
    LARGE_R,
    ZERO_R,
    UNUSED,        
};

template<SEARCH_R Tsmode>
struct GetRSearchSq{
    template<class Tepi, class Tepj>
    REAL operator () (const Tepi & epi, const Tepj & epj){
        assert(0);
    }
};

// scatter
template<>
struct GetRSearchSq<SEARCH_R::EPJ_R>{
    template<class Tepi, class Tepj>
    REAL operator () (const Tepi & epi, const Tepj & epj){
        //if(epi.id==0){
        //    std::cerr<<"epj.id= "<<epj.id
        //             <<" epj.r_search= "<<epj.r_search<<std::endl;
        //}
        return epj.r_search*epj.r_search;
    }
};

// gather
template<>
struct GetRSearchSq<SEARCH_R::EPI_R>{
    template<class Tepi, class Tepj>
    REAL operator () (const Tepi & epi, const Tepj & epj){
        return epi.r_search*epi.r_search;
    }
};

// symmetry
template<>
struct GetRSearchSq<SEARCH_R::MAX_R>{
    template<class Tepi, class Tepj>
    REAL operator () (const Tepi & epi, const Tepj & epj){
        REAL r_search = std::max(epi.r_search, epj.r_search);
        return r_search*r_search;
    }
};

// Inf
template<>
struct GetRSearchSq<SEARCH_R::LARGE_R>{
    template<class Tepi, class Tepj>
    REAL operator () (const Tepi & epi, const Tepj & epj){
        return PS::LARGE_FLOAT;
    }
};

template<>
struct GetRSearchSq<SEARCH_R::ZERO_R>{
    template<class Tepi, class Tepj>
    REAL operator () (const Tepi & epi, const Tepj & epj){
        return 0.0;
        //return -0.00001; // in the case
    }
};

///////////
// SHORT //
// multiwalk:no index:no
template<SEARCH_R Tsearch_ngb_ep, SEARCH_R Tsearch_mass_ep, int Tindex_mode>
struct CalcForceShort{
    template<class Tepi, class Tepj>
    void operator () (const Tepi * ep_i,
                      const PS::S32 n_ip,
                      const Tepj * ep_j,
                      const PS::S32 n_jp,
                      Force * force){
        assert(0);
    }
};

template<SEARCH_R Tsearch_ngb_ep, SEARCH_R Tsearch_mass_ep>
struct CalcForceShort<Tsearch_ngb_ep, Tsearch_mass_ep, 0>{
    template<class Tepi, class Tepj>
    void operator () (const Tepi * ep_i,
                      const PS::S32 n_ip,
                      const Tepj * ep_j,
                      const PS::S32 n_jp,
                      Force * force){
	for(PS::S32 i=0; i<n_ip; i++){
	    const REALVEC xi = ep_i[i].getPos();
	    for(PS::S32 j=0; j<n_jp; j++){
		const REALVEC xj = ep_j[j].getPos();
                REALVEC rij = xi - xj;
                const REAL r_search_sq_ngb = GetRSearchSq<Tsearch_ngb_ep>()(ep_i[i], ep_j[j]);
		if(rij*rij <= r_search_sq_ngb){
                    //std::cerr<<"ep_i[i].id= "<<ep_i[i].id
                    //         <<" ep_j[j].id= "<<ep_j[j].id
                    //         <<"r= "<<sqrt(rij*rij)
                    //         <<" r_search= "<<ep_j[j].r_search
                    //         <<std::endl;
		    force[i].n_ngb++;
		}
                const REAL r_search_sq_mass = GetRSearchSq<Tsearch_mass_ep>()(ep_i[i], ep_j[j]);
		if(rij*rij <= r_search_sq_mass){
                    force[i].mass += ep_j[j].mass;
		}
	    }
	}
    }
};

// multiwalk:yes index:no
template<SEARCH_R Tsearch_ngb_ep, SEARCH_R Tsearch_mass_ep>
struct CalcForceShort<Tsearch_ngb_ep, Tsearch_mass_ep, 1>{
    template<class Tepi, class Tepj>
    PS::S32 operator () (const PS::S32   tag,
                         const PS::S32   n_walk,
                         const Tepi   ** epi,
                         const PS::S32 * n_epi,
                         const Tepj   ** epj,
                         const PS::S32 * n_epj){
        PS::S32 n_cnt = 0;
        for(PS::S32 iw=0; iw<n_walk; iw++){
            assert(n_cnt + n_epi[iw] < N_EPI_LIMIT);
            for(PS::S32 i=0; i<n_epi[iw]; i++, n_cnt++){
                const REALVEC xi = epi[iw][i].getPos();
                dev_force[n_cnt].clear();
                for(PS::S32 j=0; j<n_epj[iw]; j++){
                    const REALVEC xj = epj[iw][j].getPos();
                    REALVEC rij = xi - xj;
                    const REAL r_search_sq_ngb = GetRSearchSq<Tsearch_ngb_ep>()(epi[iw][i], epj[iw][j]);
                    if(rij*rij <= r_search_sq_ngb){
                        dev_force[n_cnt].n_ngb++;
                    }
                    const REAL r_search_sq_mass = GetRSearchSq<Tsearch_mass_ep>()(epi[iw][i], epj[iw][j]);
                    if(rij*rij <= r_search_sq_mass){
                        dev_force[n_cnt].mass += epj[iw][j].mass;
                    }
                }
            }
        }
        return 0;
    }
};

// multiwalk:yes index;yes
template<SEARCH_R Tsearch_ngb_ep, SEARCH_R Tsearch_mass_ep>
struct CalcForceShort<Tsearch_ngb_ep, Tsearch_mass_ep, 2>{
    template<class Tepi, class Tepj>
    PS::S32 operator () (const PS::S32    tag,
                         const PS::S32    n_walk,
                         const Tepi    ** epi,
                         const PS::S32  * n_epi,
                         const PS::S32 ** id_epj,
                         const PS::S32  * n_epj,
                         const Tepj     * epj,
                         const PS::S32    n_epj_tot,
                         const bool       send_flag){
        //if(n_epj_tot >= N_EPJ_LIMIT){
        //if(1){
        //    std::cerr<<" n_epj_tot= "<<n_epj_tot
        //             <<" N_EPJ_LIMIT= "<<N_EPJ_LIMIT
        //             //<<" n_epi[0]= "<<n_epi[0]
        //             <<" n_walk= "<<n_walk
        //             <<std::endl;
        //}
        assert(n_epj_tot < N_EPJ_LIMIT);
        if(send_flag==true){
            for(PS::S32 i=0; i<n_epj_tot; i++){
                //dev_epj[i] = epj[i];
                dev_epj[i].copyFromEpInFdps(epj[i]);
            }
        }
        else{
            PS::S32 n_cnt = 0;
            for(PS::S32 iw=0; iw<n_walk; iw++){
                for(PS::S32 i=0; i<n_epi[iw]; i++, n_cnt++){
                    assert(n_cnt + n_epi[iw] < N_EPI_LIMIT);
                    const REALVEC xi = epi[iw][i].getPos();
                    dev_force[n_cnt].clear();
                    for(PS::S32 j=0; j<n_epj[iw]; j++){
                        const PS::S32 id_j = id_epj[iw][j];
                        const REALVEC xj = dev_epj[id_j].getPos();
                        REALVEC rij = xi - xj;
                        const REAL r_search_sq_ngb = GetRSearchSq<Tsearch_ngb_ep>()(epi[iw][i], dev_epj[id_j]);
                        if(rij*rij <= r_search_sq_ngb){
                            dev_force[n_cnt].n_ngb++;
                        }
                        const REAL r_search_sq_mass = GetRSearchSq<Tsearch_mass_ep>()(epi[iw][i], dev_epj[id_j]);
                        if(rij*rij <= r_search_sq_mass){
                            dev_force[n_cnt].mass += epj[id_j].mass;
                        }
                    }
                }
            }
        }
        return 0;
    }
};

struct RetrieveShort{
    PS::S32 operator () (const PS::S32 tag,
                         const PS::S32 n_walk,
                         const PS::S32 * ni,
                         Force ** force){
        int n_cnt = 0;
        for(int iw=0; iw<n_walk; iw++){
            assert(n_cnt + ni[iw] < N_EPI_LIMIT);
            for(int i=0; i<ni[iw]; i++, n_cnt++){
                force[iw][i].n_ngb = dev_force[n_cnt].n_ngb;
                force[iw][i].mass = dev_force[n_cnt].mass;
            }
        }
        return 0;
    }
};



//////////
// LONG //

inline REAL CalcPot(const REALVEC & xi, const REALVEC & xj, const REAL mj){
    if(xi == xj) return 0.0;
    const REALVEC dx = xi - xj;
    const REAL dx_sq = dx*dx;
    return -mj / sqrt(dx_sq);
}


template<SEARCH_R Tsearch_ngb_ep, SEARCH_R Tsearch_mass_ep, SEARCH_R Tsearch_ngb_sp, SEARCH_R Tsearch_mass_sp, int Tindex_mode>
struct CalcForceLong{
    template<class Tepi, class Tepj>
    void operator () (const Tepi * ep_i,
                      const PS::S32 n_ip,
                      const Tepj * ep_j,
                      const PS::S32 n_jp,
                      Force * force){
        assert(0);
    }
};

// multiwalk:no index:no range:long
template<SEARCH_R Tsearch_ngb_ep, SEARCH_R Tsearch_mass_ep, SEARCH_R Tsearch_ngb_sp, SEARCH_R Tsearch_mass_sp>
struct CalcForceLong<Tsearch_ngb_ep, Tsearch_mass_ep, Tsearch_ngb_sp, Tsearch_mass_sp, 0>{
    template<class Tepi, class Tepj>
    void operator () (const Tepi * ep_i,
                      const PS::S32 n_ip,
                      const Tepj * ep_j,
                      const PS::S32 n_jp,
                      Force * force){
	for(PS::S32 i=0; i<n_ip; i++){
	    const REALVEC xi = ep_i[i].getPos();
            force[i].id = ep_i[i].id;
	    for(PS::S32 j=0; j<n_jp; j++){
		const REALVEC xj = ep_j[j].getPos();
		REALVEC rij = xi - xj;
                const REAL r_search_sq_ngb = GetRSearchSq<Tsearch_ngb_ep>()(ep_i[i], ep_j[j]);
		if(rij*rij <= r_search_sq_ngb){
                    //if(force[i].id==0){
                    //    std::cout<<"ep_j[j].id= "<<ep_j[j].id
                    //             <<" pos= "<<ep_j[j].pos
                    //             <<" mass= "<<ep_j[j].mass
                    //             <<std::endl;
                    //}
		    force[i].n_ngb++;
		}
                const REAL r_search_sq_mass = GetRSearchSq<Tsearch_mass_ep>()(ep_i[i], ep_j[j]);
                if(rij*rij <= r_search_sq_mass){
                    force[i].mass += ep_j[j].mass;
                    force[i].pot += CalcPot(ep_i[i].pos, ep_j[j].pos, ep_j[j].mass);
                    if(ep_j[j].mass == 0.0){
                        std::cerr<<"ep_j[j].pos= "<<ep_j[j].pos
                                 <<" r_search_sq_mass= "<<r_search_sq_mass
                                 <<std::endl;
                    }
                    assert(ep_j[j].mass != 0.0);
		}
	    }
            //if(force[i].id==0) std::cout<<std::endl;
	}
    }
};

// multiwalk:yes index:no
template<SEARCH_R Tsearch_ngb_ep, SEARCH_R Tsearch_mass_ep, SEARCH_R Tsearch_ngb_sp, SEARCH_R Tsearch_mass_sp>
struct CalcForceLong<Tsearch_ngb_ep, Tsearch_mass_ep, Tsearch_ngb_sp, Tsearch_mass_sp, 1>{
    template<class Tepi, class Tepj, class Tspj>
    PS::S32 operator () (const PS::S32   tag,
                         const PS::S32   n_walk,
                         const Tepi   ** epi,
                         const PS::S32 * n_epi,
                         const Tepj   ** epj,
                         const PS::S32 * n_epj,
                         const Tspj   ** spj,
                         const PS::S32 * n_spj){
        PS::S32 n_cnt = 0;
        for(PS::S32 iw=0; iw<n_walk; iw++){
            assert(n_cnt + n_epi[iw] < N_EPI_LIMIT);
            for(PS::S32 i=0; i<n_epi[iw]; i++, n_cnt++){
                const REALVEC xi = epi[iw][i].getPos();
                dev_force[n_cnt].clear();
                for(PS::S32 j=0; j<n_epj[iw]; j++){
                    const REALVEC xj = epj[iw][j].getPos();
                    REALVEC rij = xi - xj;
                    const REAL r_search_sq_ngb = GetRSearchSq<Tsearch_ngb_ep>()(epi[iw][i], epj[iw][j]);
                    if(rij*rij <= r_search_sq_ngb){
                        dev_force[n_cnt].n_ngb++;
                    }
                    const REAL r_search_sq_mass = GetRSearchSq<Tsearch_mass_ep>()(epi[iw][i], epj[iw][j]);
                    if(rij*rij <= r_search_sq_mass){
                        dev_force[n_cnt].mass += epj[iw][j].mass;
                        dev_force[n_cnt].pot += CalcPot(epi[iw][i].pos, epj[iw][j].pos, epj[iw][j].mass);
                    }
                }
                for(PS::S32 j=0; j<n_spj[iw]; j++){
                    const REALVEC xj = spj[iw][j].getPos();
                    const REAL r_search_sq_ngb = GetRSearchSq<Tsearch_ngb_sp>()(epi[iw][i], spj[iw][j]);
                    REALVEC rij = xi - xj;
                    if(rij*rij <= r_search_sq_ngb){
                        dev_force[n_cnt].n_ngb++;
                    }
                    const REAL r_search_sq_mass = GetRSearchSq<Tsearch_mass_sp>()(epi[iw][i], spj[iw][j]);
                    if(rij*rij <= r_search_sq_mass){
                        dev_force[n_cnt].mass += spj[iw][j].mass;
                        dev_force[n_cnt].pot += CalcPot(epi[iw][i].pos, spj[iw][j].pos, spj[iw][j].mass);
                    }
                }
            }
        }
        return 0;
    }
};

// multiwalk:yes index;yes
template<SEARCH_R Tsearch_ngb_ep, SEARCH_R Tsearch_mass_ep, SEARCH_R Tsearch_ngb_sp, SEARCH_R Tsearch_mass_sp>
struct CalcForceLong<Tsearch_ngb_ep, Tsearch_mass_ep, Tsearch_ngb_sp, Tsearch_mass_sp, 2>{
    template<class Tepi, class Tepj, class Tspj>
    PS::S32 operator () (const PS::S32    tag,
                         const PS::S32    n_walk,
                         const Tepi    ** epi,
                         const PS::S32  * n_epi,
                         const PS::S32 ** id_epj,
                         const PS::S32  * n_epj,
                         const PS::S32 ** id_spj,
                         const PS::S32  * n_spj,
                         const Tepj     * epj,
                         const PS::S32    n_epj_tot,
                         const Tspj     * spj,
                         const PS::S32    n_spj_tot,
                         const bool       send_flag){
        assert(n_epj_tot < N_EPJ_LIMIT);
        assert(n_spj_tot < N_SPJ_LIMIT);
        if(send_flag==true){
            for(PS::S32 i=0; i<n_epj_tot; i++){
                //dev_epj[i] = epj[i];
                dev_epj[i].copyFromEpInFdps(epj[i]);
            }
            for(PS::S32 i=0; i<n_spj_tot; i++){
                //dev_spj[i] = spj[i];
                dev_spj[i].copyFromSpInFdps(spj[i]);
            }
        }
        else{
            PS::S32 n_cnt = 0;
            for(PS::S32 iw=0; iw<n_walk; iw++){
                assert(n_cnt + n_epi[iw] < N_EPI_LIMIT);
                for(PS::S32 i=0; i<n_epi[iw]; i++, n_cnt++){
                    const REALVEC xi = epi[iw][i].getPos();
                    dev_force[n_cnt].clear();
                    for(PS::S32 j=0; j<n_epj[iw]; j++){
                        const PS::S32 id_j = id_epj[iw][j];
                        const REALVEC xj = dev_epj[id_j].getPos();
                        REALVEC rij = xi - xj;
                        const REAL r_search_sq_ngb = GetRSearchSq<Tsearch_ngb_ep>()(epi[iw][i], dev_epj[id_j]);
                        if(rij*rij <= r_search_sq_ngb){
                            dev_force[n_cnt].n_ngb++;
                        }
                        const REAL r_search_sq_mass = GetRSearchSq<Tsearch_mass_ep>()(epi[iw][i], dev_epj[id_j]);
                        if(rij*rij <= r_search_sq_mass){
                            dev_force[n_cnt].mass += dev_epj[id_j].mass;
                            dev_force[n_cnt].pot += CalcPot(epi[iw][i].pos, dev_epj[id_j].pos, dev_epj[id_j].mass);
                        }
                    }
                    for(PS::S32 j=0; j<n_spj[iw]; j++){
                        const PS::S32 id_j = id_spj[iw][j];
                        const REALVEC xj = dev_spj[id_j].getPos();
                        REALVEC rij = xi - xj;
                        const REAL r_search_sq_ngb = GetRSearchSq<Tsearch_ngb_sp>()(epi[iw][i], dev_spj[id_j]);
                        if(rij*rij <= r_search_sq_ngb){
                            dev_force[n_cnt].n_ngb++;
                        }
                        const REAL r_search_sq_mass = GetRSearchSq<Tsearch_mass_sp>()(epi[iw][i], dev_spj[id_j]);
                        if(rij*rij <= r_search_sq_mass){
                            dev_force[n_cnt].mass += dev_spj[id_j].mass;
                            dev_force[n_cnt].pot  += CalcPot(epi[iw][i].pos, dev_spj[id_j].pos, dev_spj[id_j].mass);
                        }
                    }
                }
            }
        }
        return 0;
    }    
};

struct RetrieveLong{
    PS::S32 operator () (const PS::S32 tag,
                         const PS::S32 n_walk,
                         const PS::S32 * ni,
                         Force ** force){
        int n_cnt = 0;
        for(int iw=0; iw<n_walk; iw++){
            assert(n_cnt + ni[iw] < N_EPI_LIMIT);
            for(int i=0; i<ni[iw]; i++, n_cnt++){
                force[iw][i].n_ngb = dev_force[n_cnt].n_ngb;
                force[iw][i].mass  = dev_force[n_cnt].mass;
                force[iw][i].pot   = dev_force[n_cnt].pot;
            }
        }
        return 0;
    }
};

void MakePlummerModel(const double mass_glb,
                      const long long int n_glb,
                      const long long int n_loc,
                      double *& mass,
                      PS::F64vec *& pos,
                      PS::F64vec *& vel,
                      const double eng = -0.25,
                      const int seed = 0){
    assert(eng < 0.0);
    static const double PI = atan(1.0) * 4.0;
    const double r_cutoff = 22.8 / (-3.0 * PI * mass_glb * mass_glb / (64.0 * -0.25)); // 22.8 is cutoff in Nbody units
    //const double r_cutoff = 22.8 * 0.25;
    mass = new double[n_loc];
    pos = new PS::F64vec[n_loc];
    vel = new PS::F64vec[n_loc];
    PS::MTTS mt;
    mt.init_genrand( PS::Comm::getRank() );
    for(int i=0; i<n_loc; i++){
        mass[i] = mass_glb / n_glb;
        double r_tmp = 9999.9;
        while(r_tmp > r_cutoff){ 
            double m_tmp = mt.genrand_res53();
            r_tmp = 1.0 / sqrt( pow(m_tmp, (-2.0/3.0)) - 1.0);
        }
        double phi = 2.0 * PI * mt.genrand_res53();
        double cth = 2.0 * (mt.genrand_real2() - 0.5);
        double sth = sqrt(1.0 - cth*cth);
        pos[i].x = r_tmp * sth * cos(phi);
        pos[i].y = r_tmp * sth * sin(phi);
        pos[i].z = r_tmp * cth;
        while(1){
            const double v_max = 0.1;
            const double v_try = mt.genrand_res53();
            const double v_crit = v_max * mt.genrand_res53();
            if(v_crit < v_try * v_try * pow( (1.0 - v_try * v_try), 3.5) ){
                const double ve = sqrt(2.0) * pow( (r_tmp*r_tmp + 1.0), -0.25 );
                phi = 2.0 * PI * mt.genrand_res53();
                cth = 2.0 * (mt.genrand_res53() - 0.5);
                sth = sqrt(1.0 - cth*cth);
                vel[i].x = ve * v_try * sth * cos(phi);
                vel[i].y = ve * v_try * sth * sin(phi);
                vel[i].z = ve * v_try * cth;
                break;
            }
        }
    }

    PS::F64vec cm_pos = 0.0;
    PS::F64vec cm_vel = 0.0;
    double  cm_mass = 0.0;
    for(int i=0; i<n_loc; i++){
        cm_pos += mass[i] * pos[i];
        cm_vel += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
    for(int i=0; i<n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }

    const double r_scale = -3.0 * PI * mass_glb * mass_glb / (64.0 * eng);
    const double coef = 1.0 / sqrt(r_scale);
    for(int i=0; i<n_loc; i++){
        pos[i] *= r_scale;
        vel[i] *= coef;
    }

    double r_max_sq = -1.0;
    for(int i=0; i<n_loc; i++){
        if(r_max_sq < pos[i] * pos[i]){
            r_max_sq = pos[i] * pos[i];
        }
    }
}


template<class Tpsys>
void SetParticlesPlummer(Tpsys & psys,
                         const PS::S64 n_glb,
                         PS::S32 & n_loc,  
                         REAL & t_sys){

    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);

    PS::F64 * mass;
    PS::F64vec * pos;
    PS::F64vec * vel;
    t_sys = 0.0;

    const PS::F64 m_tot = 1.0;
    const PS::F64 eng = -0.25;
    MakePlummerModel(m_tot, n_glb, n_loc, mass, pos, vel, eng);
    PS::S64 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    for(int i=0; i<n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}

void MakeUniformQubeModel(const double mass_glb,
                          const long long int n_glb,
                          const long long int n_loc,
                          double *& mass,
                          PS::F64vec *& pos,
                          PS::F64vec *& vel,
                          PS::F64 full_len = 1.0,
                          const PS::F64vec offset = PS::F64vec(0.0, 0.0, 0.0),
                          const int seed = 0){

    mass = new double[n_loc];
    pos = new PS::F64vec[n_loc];
    vel = new PS::F64vec[n_loc];
    PS::MTTS mt;
    mt.init_genrand( PS::Comm::getRank() );
    for(int i=0; i<n_loc; i++){
        mass[i] = mass_glb / n_glb;
        pos[i].x = full_len * mt.genrand_real2() - offset.x;
        pos[i].y = full_len * mt.genrand_real2() - offset.y;
        pos[i].z = full_len * mt.genrand_real2() - offset.z;
        vel[i].x = (mt.genrand_real2() - 0.5);
        vel[i].y = (mt.genrand_real2() - 0.5);
        vel[i].z = (mt.genrand_real2() - 0.5);
    }
}


template<class Tpsys>
void SetParticlesUniformQube(Tpsys & psys,
                             const PS::S64 n_glb,
                             PS::S32 & n_loc,  
                             REAL & t_sys,
                             const PS::F64 full_len = 1.0,
                             const PS::F64vec offset = PS::F64vec(0.0, 0.0, 0.0),
                             const IcMode ic=IcMode::uni_qube_rand){
    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    n_loc = n_glb / n_proc;
    if( n_glb % n_proc > my_rank) n_loc++;
    psys.setNumberOfParticleLocal(n_loc);
    PS::S64 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    
    PS::F64 * mass;
    PS::F64vec * pos;
    PS::F64vec * vel;
    t_sys = 0.0;
    const PS::F64 m_tot = 1.0;
    const PS::S32 seed = 0;
    if(ic==IcMode::uni_qube_rand){
        MakeUniformQubeModel(m_tot, n_glb, n_loc, mass, pos, vel, full_len, offset, seed);
    }
    else if(ic==IcMode::uni_qube_grid){
        PS::S32 n_glb_1d = cbrt(n_glb);
        assert(n_glb == (n_glb_1d*n_glb_1d*n_glb_1d) );
        PS::F64 dx = full_len / n_glb_1d;
        mass = new double[n_loc];
        pos = new PS::F64vec[n_loc];
        vel = new PS::F64vec[n_loc];
        PS::S64 i_e = i_h + n_loc;
        PS::S32 n_cnt = 0;
        for(PS::S32 i=i_h; i<i_e; i++, n_cnt++){
            PS::S32 ix = i / (n_glb_1d*n_glb_1d);
            PS::S32 iy = (i/n_glb_1d) % n_glb_1d;
            PS::S32 iz = i%n_glb_1d;
            pos[n_cnt].x = ix*dx - offset.x;
            pos[n_cnt].y = iy*dx - offset.y;
            pos[n_cnt].z = iz*dx - offset.z;
        }
    }
    else if(ic==IcMode::uni_qube_ill){
        PS::S32 n_glb_1d = cbrt(n_glb);
        assert(n_glb == (n_glb_1d*n_glb_1d*n_glb_1d) );
        PS::F64 dx = full_len / n_glb_1d;
        mass = new double[n_loc];
        pos = new PS::F64vec[n_loc];
        vel = new PS::F64vec[n_loc];
        PS::S64 i_e = i_h + n_loc;
        PS::S32 n_cnt = 0;
        for(PS::S32 i=i_h; i<i_e; i++, n_cnt++){
            PS::S32 ix = i / (n_glb_1d*n_glb_1d);
            //PS::S32 iy = (i/n_glb_1d) % n_glb_1d;
            PS::S32 iz = i%n_glb_1d;
            pos[n_cnt].x = ix*dx - offset.x;
            //pos[n_cnt].y = iy*dx - offset.y;
            pos[n_cnt].y = -offset.y;
            pos[n_cnt].z = iz*dx - offset.z;
        }
    }
    for(int i=0; i<n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}

void printHelp() {
    std::cerr<<"w: n_walk_limit (default: 10)"<<std::endl;
    std::cerr<<"t: theta (default: 0.5)"<<std::endl;
    std::cerr<<"l: n_leaf_limit (default: 8)"<<std::endl;
    std::cerr<<"n: n_group_limit (default: 64)"<<std::endl;
    std::cerr<<"N: n_tot (default: 512)"<<std::endl;
    std::cerr<<"h: help"<<std::endl;
}

template<class Tsys>
void ShiftPos(Tsys & sys){
    const PS::S32 n_loc = sys.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_loc; i++){
        sys[i].pos += 1.0;
    }
}

/*
template<class TNSM, class Tsys, class Tepj>
class GetNeighborClass{
    static void func(Tsys & system, Tepj & epj){}
};

template<class Tsys, class Tepj>
class GetNeighborClass<PS::TageNeighborSearchNo, Tsys, Tepj>{
    static void func(Tsys & system, Tepj & epj){}
};

template<class TNSM, class Tsys, class Tepj>
void GetNeighbor(Tsys & system,
                 Tepj & epj){
    GetNeighborClass<TNSM, Tsys, Tepj>::func(system, epj);
}
*/

template<class Tsys, class Ttree, class Tfunc_ep, class Tfunc_disp,
         class Tfunc_disp_id, class Tfunc_ret, class Tepj>
void CompareForceShort(Tsys & system, Ttree & tree, PS::DomainInfo & dinfo,
                       Tfunc_ep func_ep, Tfunc_disp func_disp,
                       Tfunc_disp_id func_disp_id, Tfunc_ret func_ret,
                       const int tag_max, const int n_walk_limit,
                       const int n_tot, const double theta,
                       const int n_leaf_limit, const int n_group_limit,
                       const IcMode ic_mode, const BcMode bc_mode, const FrMode fr_mode,
                       const std::vector<REALVEC> & pos_org,
                       std::vector<TestResult> & test_result,
                       Tepj){
    PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"Compare short force "<<std::endl;
    tree.initialize(n_tot, theta, n_leaf_limit, n_group_limit);
    tree.clearCounterAll();
    const PS::S32 n_loc = system.getNumberOfParticleLocal();
    Force * force_dir = new Force[n_loc];

    //tree.calcForceAll(func_ep, system, dinfo, true, PS::MAKE_LIST);
    //PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"start diretc force (Compare short force)"<<std::endl;
    tree.setParticleLocalTree(system, true);
    tree.setRootCell(dinfo);
    tree.clearCounterAll();
    tree.calcForceDirect(func_ep, force_dir, dinfo, true);
    PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"finish diretc force (Compare short force)"<<std::endl;

    PS::Comm::barrier();
    int n_interaction_list_mode = 3;

    ///////////////
    // multiwalk:no, index: no, id_mode: 0
    if(PS::Comm::getRank()==0){
        std::cerr<<"multiwalk:no, index:no, id_mode:0"<<std::endl;
    }
    for(int il=0; il<n_interaction_list_mode; il++){
    //for(int il=0; il<2; il++){
        PS::Comm::barrier();
        if(PS::Comm::getRank()==0){
            std::cerr<<"start: il= "<<il<<std::endl;
        }
        enum PS::INTERACTION_LIST_MODE il_mode = static_cast<enum PS::INTERACTION_LIST_MODE>(il_mode_ar[il]);
        if(il_mode != PS::REUSE_LIST) tree.freeMem();
        tree.clearCounterAll();
        tree.calcForceAll(func_ep, system, dinfo, true, il_mode);
        Force * force_tree = new Force[n_loc];
        for(int i=0; i<n_loc; i++) force_tree[i] = tree.getForce(i);
        PS::S64 n_err_ngb = 0;
        PS::S64 n_err_mass = 0;
        for(int i=0; i<n_loc; i++){
            if(force_tree[i].n_ngb != force_dir[i].n_ngb) n_err_ngb++;
            if(force_tree[i].mass != force_dir[i].mass) n_err_mass++;
            Tepj * next = NULL;
            PS::S32 n_ngb_tmp = tree.getNeighborListOneParticle(system[i], next);
            if(n_ngb_tmp != force_dir[i].n_ngb){
                //PS::S32 n_ngb_tmp_tmp = 0;
                //for(int j=0; j<n_loc; j++){
                //    PS::F64 dr = sqrt( (system[i].pos-system[j].pos)*(system[i].pos-system[j].pos) );
                //    if(dr < system[j].r_search){
                //        n_ngb_tmp_tmp++;
                //    }
                //}
                //std::cerr<<"system[i].id= "<<system[i].id<<std::endl;
                //std::cerr<<"system[i].pos= "<<system[i].pos
                //         <<" system[i].r_search= "<<system[i].r_search
                //         <<std::endl;
                //std::cerr<<"tree.epj_org_[i].id= "<<tree.epj_org_[i].id
                //         <<"tree.epj_org_[i].r_search= "<<tree.epj_org_[i].getRSearch()
                //         <<std::endl;
                //std::cerr<<"n_ngb_tmp= "<<n_ngb_tmp
                    //<<"n_ngb_tmp_tmp= "<<n_ngb_tmp_tmp
                //<<" force_dir[i].n_ngb= "<<force_dir[i].n_ngb
                //       <<std::endl;
                n_err_ngb++;
            }
        }
        PS::S64 n_err_ngb_glb = PS::Comm::getSum(n_err_ngb);
        PS::S64 n_err_mass_glb = PS::Comm::getSum(n_err_mass);
        
        bool pass_glb = true;
        if(n_err_ngb_glb != 0 || n_err_mass_glb != 0) pass_glb = false;
        test_result.push_back(TestResult(pass_glb, ic_mode, bc_mode, fr_mode, il_mode, 0, n_err_ngb_glb, n_err_mass_glb));
        // NOTE: if reusing is on and shift FPs, # of neighbors at the
        // reusing step is not correct, becuase at the reusing step
        // the information of vertex is not updated.
#if 0
        if(il_mode == static_cast<int>(PS::MAKE_LIST_FOR_REUSE)) ShiftPos(system);
#endif
        delete [] force_tree;
        if(PS::Comm::getRank()==0) std::cerr<<"finish: il= "<<il<<std::endl;
    }
    for(PS::S32 i=0; i<n_loc; i++) system[i].pos = pos_org[i];

    ///////////////
    // multiwalk:yes, index: no, id_mode: 1
    if(PS::Comm::getRank()==0){
        std::cerr<<"multiwalk:yes, index:no, id_mode:1"<<std::endl;
    }
    for(int il=0; il<n_interaction_list_mode; il++){
        PS::Comm::barrier();
        if(PS::Comm::getRank()==0){
            std::cerr<<"start: il= "<<il<<std::endl;
        }
        enum PS::INTERACTION_LIST_MODE il_mode = static_cast<enum PS::INTERACTION_LIST_MODE>(il_mode_ar[il]);
        if(il_mode != PS::REUSE_LIST) tree.freeMem();
        tree.clearCounterAll();
        tree.calcForceAllMultiWalk(func_disp, func_ret, tag_max, system, dinfo, n_walk_limit, true, il_mode);
        Force * force_tree = new Force[n_loc];
        for(int i=0; i<n_loc; i++) force_tree[i] = tree.getForce(i);
        PS::S64 n_err_ngb = 0;
        PS::S64 n_err_mass = 0;
        for(int i=0; i<n_loc; i++){
            if(force_tree[i].n_ngb != force_dir[i].n_ngb) n_err_ngb++;
            if(force_tree[i].mass != force_dir[i].mass) n_err_mass++;
            Tepj * next = NULL;
            PS::S32 n_ngb_tmp = tree.getNeighborListOneParticle(system[i], next);
            if(n_ngb_tmp != force_dir[i].n_ngb){
                n_err_ngb++;
            }
        }
        PS::S64 n_err_ngb_glb = PS::Comm::getSum(n_err_ngb);
        PS::S64 n_err_mass_glb = PS::Comm::getSum(n_err_mass);
        bool pass_glb = true;
        if(n_err_ngb_glb != 0 || n_err_mass_glb != 0) pass_glb = false;
        test_result.push_back(TestResult(pass_glb, ic_mode, bc_mode, fr_mode, il_mode, 1, n_err_ngb_glb, n_err_mass_glb));
#if 0
        if(il_mode == static_cast<int>(PS::MAKE_LIST_FOR_REUSE)) ShiftPos(system);
#endif
        delete [] force_tree;
        if(PS::Comm::getRank()==0) std::cerr<<"fhinish il= "<<il<<std::endl;
    }
    for(PS::S32 i=0; i<n_loc; i++) system[i].pos = pos_org[i];

    ///////////////
    // multiwalk:yes, index: no, id_mode: 2
    if(PS::Comm::getRank()==0){
        std::cerr<<"multiwalk:yes, index:no, id_mode:2"<<std::endl;
    }
    for(int il=0; il<n_interaction_list_mode; il++){
        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"start: il= "<<il<<std::endl;
        enum PS::INTERACTION_LIST_MODE il_mode = static_cast<enum PS::INTERACTION_LIST_MODE>(il_mode_ar[il]);
        if(il_mode != PS::REUSE_LIST) tree.freeMem();
        tree.clearCounterAll();
        tree.calcForceAllMultiWalkIndex(func_disp_id, func_ret, tag_max, system, dinfo, n_walk_limit, true, il_mode);
        Force * force_tree = new Force[n_loc];
        for(int i=0; i<n_loc; i++) force_tree[i] = tree.getForce(i);
        PS::S64 n_err_ngb = 0;
        PS::S64 n_err_mass = 0;
        for(int i=0; i<n_loc; i++){
            if(force_tree[i].n_ngb != force_dir[i].n_ngb) n_err_ngb++;
            if(force_tree[i].mass != force_dir[i].mass) n_err_mass++;
            Tepj * next = NULL;
            PS::S32 n_ngb_tmp = tree.getNeighborListOneParticle(system[i], next);
            if(n_ngb_tmp != force_dir[i].n_ngb){
                n_err_ngb++;
            }
        }
        PS::S64 n_err_ngb_glb = PS::Comm::getSum(n_err_ngb);
        PS::S64 n_err_mass_glb = PS::Comm::getSum(n_err_mass);
        bool pass_glb = true;
        if(n_err_ngb_glb != 0 || n_err_mass_glb != 0) pass_glb = false;
        test_result.push_back(TestResult(pass_glb, ic_mode, bc_mode, fr_mode, il_mode, 2, n_err_ngb_glb, n_err_mass_glb));
#if 0
        if(il_mode == static_cast<int>(PS::MAKE_LIST_FOR_REUSE)) ShiftPos(system);
#endif
        delete [] force_tree;
        if(PS::Comm::getRank()==0) std::cerr<<"fhinish il= "<<il<<std::endl;
    }
    for(PS::S32 i=0; i<n_loc; i++) system[i].pos = pos_org[i];

    delete [] force_dir;
}

template<class Tsys, class Ttree, class Tfunc_ep, class Tfunc_sp,
         class Tfunc_disp, class Tfunc_disp_id, class Tfunc_ret, class Tepj>
void CompareForceLong(Tsys & system, Ttree & tree, PS::DomainInfo & dinfo,
                      Tfunc_ep func_ep, Tfunc_sp func_sp, Tfunc_disp func_disp,
                      Tfunc_disp_id func_disp_id, Tfunc_ret func_ret,
                      const int tag_max, const int n_walk_limit,
                      const int n_tot, const double theta,
                      const int n_leaf_limit, const int n_group_limit,
                      const IcMode ic_mode, const BcMode bc_mode, const FrMode fr_mode,
                      const std::vector<REALVEC> & pos_org,
                      std::vector<TestResult> & test_result,
                      Tepj){

    std::cerr<<"0) CompareForceLong"<<std::endl;
        
    tree.initialize(n_tot, theta, n_leaf_limit, n_group_limit);
    const PS::S32 n_loc = system.getNumberOfParticleLocal();
    Force * force_dir = new Force[n_loc];
    //tree.calcForceAll(func_ep, func_sp, system, dinfo, true, PS::MAKE_LIST);
    tree.setParticleLocalTree(system, true);
    tree.setRootCell(dinfo);
    tree.clearCounterAll();
    tree.calcForceDirect(func_ep, force_dir, dinfo, true);
    int n_interaction_list_mode = 3; // original

    std::cerr<<"1) CompareForceLong"<<std::endl;
    // multiwalk:no, index: no, id_mode: 0
    for(int il=0; il<n_interaction_list_mode; il++){
        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"Long) start: il= "<<il<<std::endl;
        enum PS::INTERACTION_LIST_MODE il_mode = static_cast<enum PS::INTERACTION_LIST_MODE>(il_mode_ar[il]);
        if(il_mode != PS::REUSE_LIST) tree.freeMem();
        tree.clearCounterAll();
        tree.calcForceAll(func_ep, func_sp, system, dinfo, true, il_mode);
        Force * force_tree = new Force[n_loc];
        for(int i=0; i<n_loc; i++) force_tree[i] = tree.getForce(i);
        PS::S64 n_err_ngb  = 0;
        PS::S64 n_err_mass = 0;
        PS::S64 n_err_pot  = 0;
        for(int i=0; i<n_loc; i++){
            if(force_tree[i].n_ngb != force_dir[i].n_ngb){
                //std::cerr<<"force_tree[i].n_ngb= "<<force_tree[i].n_ngb
                //         <<" force_dir[i].n_ngb= "<<force_dir[i].n_ngb
                //         <<std::endl;
                n_err_ngb++;
            }
            if(force_tree[i].mass != force_dir[i].mass){
                n_err_mass++;
            }
            if( fabs( (force_tree[i].pot - force_dir[i].pot)/force_dir[i].pot) > dEcrit ){
                n_err_pot++;
            }
            Tepj * next = NULL;
            PS::S32 n_ngb_tmp = tree.getNeighborListOneParticle(system[i], next);
            if(fr_mode == static_cast<int>(PS::LONG_SCATTER) || (fr_mode == static_cast<int>(PS::LONG_SYMMETRY))){
                if(n_ngb_tmp != force_dir[i].n_ngb){
                    //std::cerr<<"n_ngb_tmp= "<<n_ngb_tmp
                    //         <<" force_dir[i].n_ngb= "<<force_dir[i].n_ngb
                    //         <<std::endl;
                    n_err_ngb++;
                }
            }
        }
        PS::S64 n_err_ngb_glb = PS::Comm::getSum(n_err_ngb);
        PS::S64 n_err_mass_glb = PS::Comm::getSum(n_err_mass);
        PS::S64 n_err_pot_glb = PS::Comm::getSum(n_err_pot);
        bool pass_glb = true;
        if(n_err_ngb_glb != 0 || n_err_mass_glb != 0 || n_err_pot_glb != 0) pass_glb = false;
        test_result.push_back(TestResult(pass_glb, ic_mode, bc_mode, fr_mode, il_mode, 0, n_err_ngb_glb, n_err_mass_glb, n_err_pot_glb));
        if(il_mode == static_cast<int>(PS::MAKE_LIST_FOR_REUSE)) ShiftPos(system);
        delete [] force_tree;
        if(PS::Comm::getRank()==0) std::cerr<<"Long) finish: il= "<<il<<std::endl;
    }
    for(PS::S32 i=0; i<n_loc; i++) system[i].pos = pos_org[i];


    
    std::cerr<<"2) CompareForceLong"<<std::endl;
    // multiwalk:yes, index: no, id_mode = 1
    for(int il=0; il<n_interaction_list_mode; il++){
        enum PS::INTERACTION_LIST_MODE il_mode = static_cast<enum PS::INTERACTION_LIST_MODE>(il_mode_ar[il]);
        if(il_mode != PS::REUSE_LIST) tree.freeMem();
        tree.clearCounterAll();
        tree.calcForceAllMultiWalk(func_disp, func_ret, tag_max, system, dinfo, n_walk_limit, true, il_mode);
        Force * force_tree = new Force[n_loc];
        for(int i=0; i<n_loc; i++) force_tree[i] = tree.getForce(i);
        PS::S64 n_err_ngb = 0;
        PS::S64 n_err_mass = 0;
        PS::S64 n_err_pot = 0;
        for(int i=0; i<n_loc; i++){
            if(force_tree[i].n_ngb != force_dir[i].n_ngb){
                n_err_ngb++;
            }
            if(force_tree[i].mass != force_dir[i].mass){
                n_err_mass++;
            }
            if( fabs( (force_tree[i].pot - force_dir[i].pot)/force_dir[i].pot) > dEcrit ){
                n_err_pot++;
            }
        }
        PS::S64 n_err_ngb_glb  = PS::Comm::getSum(n_err_ngb);
        PS::S64 n_err_mass_glb = PS::Comm::getSum(n_err_mass);
        PS::S64 n_err_pot_glb  = PS::Comm::getSum(n_err_pot);
        bool pass_glb = true;
        if(n_err_ngb_glb != 0 || n_err_mass_glb != 0 || n_err_pot_glb != 0) pass_glb = false;
        test_result.push_back(TestResult(pass_glb, ic_mode, bc_mode, fr_mode, il_mode, 1, n_err_ngb_glb, n_err_mass_glb, n_err_pot_glb));
        if(il_mode == static_cast<int>(PS::MAKE_LIST_FOR_REUSE)) ShiftPos(system);
        delete [] force_tree;
    }
    for(PS::S32 i=0; i<n_loc; i++) system[i].pos = pos_org[i];

    std::cerr<<"3) CompareForceLong"<<std::endl;
    // multiwalk:yes, index: yes, id_mode = 2
    for(int il=0; il<n_interaction_list_mode; il++){
        enum PS::INTERACTION_LIST_MODE il_mode = static_cast<enum PS::INTERACTION_LIST_MODE>(il_mode_ar[il]);
        if(il_mode != PS::REUSE_LIST) tree.freeMem();
        tree.clearCounterAll();
        tree.calcForceAllMultiWalkIndex(func_disp_id, func_ret, tag_max, system, dinfo, n_walk_limit, true, il_mode);
        Force * force_tree = new Force[n_loc];
        for(int i=0; i<n_loc; i++) force_tree[i] = tree.getForce(i);
        PS::S64 n_err_ngb  = 0;
        PS::S64 n_err_mass = 0;
        PS::S64 n_err_pot  = 0;
        for(int i=0; i<n_loc; i++){
            if(force_tree[i].n_ngb != force_dir[i].n_ngb){
                n_err_ngb++;
            }
            if(force_tree[i].mass != force_dir[i].mass){
                n_err_mass++;
            }
            if( fabs( (force_tree[i].pot - force_dir[i].pot)/force_dir[i].pot) > dEcrit ){
                n_err_pot++;
            }
        }
        PS::S64 n_err_ngb_glb = PS::Comm::getSum(n_err_ngb);
        PS::S64 n_err_mass_glb = PS::Comm::getSum(n_err_mass);
        PS::S64 n_err_pot_glb = PS::Comm::getSum(n_err_pot);
        bool pass_glb = true;
        if(n_err_ngb_glb != 0 || n_err_mass_glb != 0 || n_err_pot_glb != 0 ) pass_glb = false;
        test_result.push_back(TestResult(pass_glb, ic_mode, bc_mode, fr_mode, il_mode, 2, n_err_ngb_glb, n_err_mass_glb, n_err_pot_glb));
        if(il_mode == static_cast<int>(PS::MAKE_LIST_FOR_REUSE)) ShiftPos(system);
        delete [] force_tree;
    }
    for(PS::S32 i=0; i<n_loc; i++) system[i].pos = pos_org[i];
    std::cerr<<"4) CompareForceLong"<<std::endl;
    
    delete [] force_dir;

}

int main(int argc, char *argv[]) {
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
#if 0
    const int n_bc_mode = 1;
    int bc_mode_ar[n_bc_mode];
    bc_mode_ar[0] = PS::BOUNDARY_CONDITION_OPEN;
    //bc_mode_ar[0] = PS::BOUNDARY_CONDITION_PERIODIC_XYZ;
    //bc_mode_ar[0] = PS::BOUNDARY_CONDITION_PERIODIC_Y;
#else
    const int n_bc_mode = 3;
    int bc_mode_ar[n_bc_mode];
    bc_mode_ar[0] = PS::BOUNDARY_CONDITION_OPEN;
    bc_mode_ar[1] = PS::BOUNDARY_CONDITION_PERIODIC_XYZ;
    bc_mode_ar[2] = PS::BOUNDARY_CONDITION_PERIODIC_Y;
#endif

#if 0
    const FrMode n_fr_mode = 1;
    FrMode fr_mode_ar[n_fr_mode];
    //fr_mode_ar[0] = static_cast<int>(PS::SHORT_SYMMETRY);
    //fr_mode_ar[0] = static_cast<int>(PS::SHORT_SCATTER);
    //fr_mode_ar[0] = static_cast<int>(PS::LONG_NO_CUTOFF);
    //fr_mode_ar[0] = static_cast<int>(PS::LONG_SYMMETRY);
    fr_mode_ar[0] = static_cast<int>(PS::SHORT_GATHER);    
#elif 0
    const FrMode n_fr_mode = 4;
    FrMode fr_mode_ar[n_fr_mode];
    fr_mode_ar[0] = static_cast<int>(PS::LONG_NO_CUTOFF);
    fr_mode_ar[1] = static_cast<int>(PS::LONG_CUTOFF);
    fr_mode_ar[2] = static_cast<int>(PS::LONG_SYMMETRY);
    fr_mode_ar[3] = static_cast<int>(PS::LONG_SCATTER);
#elif 0
    const FrMode n_fr_mode = 3;
    FrMode fr_mode_ar[n_fr_mode];
    fr_mode_ar[0] = static_cast<int>(PS::SHORT_SCATTER);
    fr_mode_ar[1] = static_cast<int>(PS::SHORT_GATHER);
    fr_mode_ar[2] = static_cast<int>(PS::SHORT_SYMMETRY);
#else
    const FrMode n_fr_mode = 7;
    FrMode fr_mode_ar[n_fr_mode];
    fr_mode_ar[0] = static_cast<int>(PS::SHORT_SCATTER);
    fr_mode_ar[1] = static_cast<int>(PS::SHORT_GATHER);
    fr_mode_ar[2] = static_cast<int>(PS::SHORT_SYMMETRY);
    fr_mode_ar[3] = static_cast<int>(PS::LONG_NO_CUTOFF);
    fr_mode_ar[4] = static_cast<int>(PS::LONG_CUTOFF);
    fr_mode_ar[5] = static_cast<int>(PS::LONG_SYMMETRY);
    fr_mode_ar[6] = static_cast<int>(PS::LONG_SCATTER);
#endif
    il_mode_ar[0] = static_cast<int>(PS::MAKE_LIST);
    il_mode_ar[1] = static_cast<int>(PS::MAKE_LIST_FOR_REUSE);
    il_mode_ar[2] = static_cast<int>(PS::REUSE_LIST);
    
    PS::Initialize(argc, argv);
    
    PS::S32 tag_max = 1;
    PS::S32 n_walk_limit = 10;
    REAL theta = 0.5;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::S64 n_tot = 512;
    
    PS::S32 c;
    while((c=getopt(argc,argv,"w:t:l:n:N:hs:e:")) != -1){
        switch(c){
        case 'w':
            n_walk_limit = atoi(optarg);
            if(PS::Comm::getRank() == 0) std::cerr << "n_walk_limit= " << n_walk_limit << std::endl;
            break;
        case 't':
            theta = atof(optarg);
            if(PS::Comm::getRank() == 0) std::cerr << "theta= " << theta << std::endl;
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            if(PS::Comm::getRank() == 0) std::cerr << "n_leaf_limit= " << n_leaf_limit << std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            if(PS::Comm::getRank() == 0) std::cerr << "n_group_limit= " << n_group_limit << std::endl;
            break;
        case 'N':
            n_tot = atoi(optarg);
            if(PS::Comm::getRank() == 0) std::cerr << "n_tot= " << n_tot << std::endl;
            break;
        case 'e':
            dEcrit = atof(optarg);
            if(PS::Comm::getRank() == 0) std::cerr << "dEcrit= " << dEcrit << std::endl;
            break;
        case 'h':
            if(PS::Comm::getRank() == 0) {
                printHelp();
            }
            PS::Finalize();
            return 0;
        default:
            if(PS::Comm::getRank() == 0) {
                std::cerr<<"No such option! Available options are here."<<std::endl;
                printHelp();
            }
            PS::Abort();
        }
    }
    PS::F64 full_len = 1.0;
    PS::F64vec offset = PS::F64vec(0.0);

    PS::S32 n_loc    = 0;
    REAL time_sys = 0.0;
    std::vector<TestResult> test_result;
    for(auto ic_mode : IcMode()){
        std::cerr<<"--------------------"<<std::endl;
        std::cerr<<"INITIAL CONDITION: "<<static_cast<int>(ic_mode)<<std::endl;
        for(int bc=0; bc<n_bc_mode; bc++){
            std::cerr<<"--------------------"<<std::endl;
            std::cerr<<"BOUNDARY_CONDITION: ";
            //if(bc==0){
            //    std::cerr<<"OPEN"<<std::endl;
            //}
            //else if(bc==1){
            //    std::cerr<<"PERIODIC_XYZ"<<std::endl;
            //}
            
            PS::ParticleSystem<FP> system;
            system.initialize();
            system.setAverageTargetNumberOfSampleParticlePerProcess(200);
            SetParticlesUniformQube(system, n_tot, n_loc, time_sys, full_len, offset, ic_mode);
            PS::F64 dis_ave  = full_len / cbrt((PS::F64)(n_tot));
            //PS::F64 dis_ave  = full_len*0.5;
            PS::F64 mass_ave = 1.0 / n_tot;
            std::mt19937 mt;
            std::uniform_real_distribution<double> rand_real(0.1, 3.0);
            std::uniform_int_distribution<int> rand_int(1, 3);
            for(PS::S32 i=0; i<n_loc; i++){
                system[i].r_search_0 = dis_ave * 1.1;
                system[i].r_search_1 = dis_ave * rand_real(mt);
                //system[i].mass  = mass_ave * rand_int(mt);
                system[i].mass  = mass_ave;
            }
            PS::Comm::barrier();
            if(PS::Comm::getRank()==0){
                const int n_tmp = std::min(n_loc, 5);
                for(PS::S32 i=0; i<n_tmp; i++){
                    std::cerr<<"system[i].r_search_0= "<<system[i].r_search_0<<" r_search_1= "<<system[i].r_search_1<<" mass= "<<system[i].mass<<std::endl;
                }
            }
            PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"create domain info object"<<std::endl;

            PS::DomainInfo dinfo;
            PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"init domain info object"<<std::endl;
            dinfo.initialize();
            PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"collect domain info object"<<std::endl;
            dinfo.collectSampleParticle(system);
            enum PS::BOUNDARY_CONDITION bc_mode = static_cast<enum PS::BOUNDARY_CONDITION>(bc_mode_ar[bc]);

            dinfo.setBoundaryCondition(bc_mode);
            PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"set root domain"<<std::endl;
            dinfo.setPosRootDomain(offset, offset+PS::F64vec(full_len));

            PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"decompose domain"<<std::endl;

            dinfo.decomposeDomain();
            PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"n_loc(before) "<<n_loc<<std::endl;
            system.exchangeParticle(dinfo);
            n_loc = system.getNumberOfParticleLocal();
            PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"n_loc(after)= "<<n_loc<<std::endl;
            std::vector<REALVEC> pos_org;
            for(PS::S32 i=0; i<n_loc; i++) pos_org.push_back(system[i].pos);
            PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"start calc"<<std::endl;

            for(int fr=0; fr<n_fr_mode; fr++){
                enum PS::SEARCH_MODE fr_mode = static_cast<enum PS::SEARCH_MODE>(fr_mode_ar[fr]);

                if( static_cast<PS::SEARCH_MODE>(fr_mode) == PS::SHORT_SCATTER){
                    PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"start PS::SHORT_SCATTER"<<std::endl;
                    SetRsearchRand(system);
                    //PS::TreeForForceShort<Force, EpAll, EpAll>::Scatter tree;
                    PS::TreeForForceShort<Force, EpSimple, EpAll>::Scatter tree;
                    tree.setExchagneLETMode(EX_LET_MODE);
                    EpAll epj_tmp;
                    CompareForceShort(system, tree, dinfo,
                                      CalcForceShort<SEARCH_R::EPJ_R,  SEARCH_R::EPJ_R, 0>(),
                                      CalcForceShort<SEARCH_R::EPJ_R,  SEARCH_R::EPJ_R, 1>(), // scatter, index:no  multiwalk:yes
                                      CalcForceShort<SEARCH_R::EPJ_R,  SEARCH_R::EPJ_R, 2>(), // scatter, index:yes multiwalk:yes
                                      RetrieveShort(),
                                      tag_max, n_walk_limit,
                                      n_tot, theta, n_leaf_limit, n_group_limit,
                                      ic_mode, bc_mode, fr_mode, pos_org, test_result,
                                      epj_tmp);
                    PS::Comm::barrier();
                    if(PS::Comm::getRank()==0) std::cerr<<"finish PS::SHORT_SCATTER"<<std::endl;
                }
                
                if( static_cast<PS::SEARCH_MODE>(fr_mode) == PS::SHORT_GATHER){
                    PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"start PS::SHORT_GATHER"<<std::endl;
                    SetRsearchRand(system);
                    PS::TreeForForceShort<Force, EpAll, EpSimple>::Gather tree;
                    tree.setExchagneLETMode(EX_LET_MODE);
                    EpSimple epj_tmp;
                    CompareForceShort(system, tree, dinfo,
                                      CalcForceShort<SEARCH_R::EPI_R,  SEARCH_R::EPI_R, 0>(), 
                                      CalcForceShort<SEARCH_R::EPI_R,  SEARCH_R::EPI_R, 1>(), // index:no  multiwalk:yes
                                      CalcForceShort<SEARCH_R::EPI_R,  SEARCH_R::EPI_R, 2>(), // index:yes multiwalk:yes
                                      RetrieveShort(),
                                      tag_max, n_walk_limit,
                                      n_tot, theta, n_leaf_limit, n_group_limit,
                                      ic_mode, bc_mode, fr_mode, pos_org, test_result,
                                      epj_tmp);
                    PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"finish PS::SHORT_GATHER"<<std::endl;
                }
                

                if( static_cast<PS::SEARCH_MODE>(fr_mode) == PS::SHORT_SYMMETRY){
                    PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"start PS::SHORT_SYMMETRY"<<std::endl;
                    SetRsearchRand(system);
                    PS::TreeForForceShort<Force, EpAll, EpAll>::Symmetry tree;
                    tree.setExchagneLETMode(EX_LET_MODE);
                    EpAll epj_tmp;
                    CompareForceShort(system, tree, dinfo,
                                      CalcForceShort<SEARCH_R::MAX_R,  SEARCH_R::MAX_R, 0>(), 
                                      CalcForceShort<SEARCH_R::MAX_R,  SEARCH_R::MAX_R, 1>(), // index:no  multiwalk:yes
                                      CalcForceShort<SEARCH_R::MAX_R,  SEARCH_R::MAX_R, 2>(), // index:yes multiwalk:yes
                                      RetrieveShort(),
                                      tag_max, n_walk_limit,
                                      n_tot, theta, n_leaf_limit, n_group_limit,
                                      ic_mode, bc_mode, fr_mode, pos_org, test_result,
                                      epj_tmp);
                    PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"finish PS::SHORT_SYMMETRY"<<std::endl;
                }


                if( static_cast<PS::SEARCH_MODE>(fr_mode) == PS::LONG_NO_CUTOFF){
                    if(static_cast<PS::BOUNDARY_CONDITION>(bc_mode) == PS::BOUNDARY_CONDITION_OPEN) {
                        PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"start PS::LONG_NO_CUTOFF"<<std::endl;
                        SetRsearchAve(system);
                        //PS::TreeForForce<PS::SEARCH_MODE_LONG, Force, EpAll, EpAll, MomentAll, MomentAll, SpAll> tree;
                        PS::TreeForForceLong<Force, EpSimple, EpSimple>::Monopole tree;
                        tree.setExchagneLETMode(EX_LET_MODE);
                        EpSimple epj_tmp;
                        CompareForceLong(system, tree, dinfo,
                                         CalcForceLong<SEARCH_R::ZERO_R,  SEARCH_R::LARGE_R, SEARCH_R::UNUSED, SEARCH_R::UNUSED, 0>(),
                                         CalcForceLong<SEARCH_R::ZERO_R,  SEARCH_R::LARGE_R, SEARCH_R::UNUSED, SEARCH_R::UNUSED, 0>(), 
                                         CalcForceLong<SEARCH_R::ZERO_R,  SEARCH_R::LARGE_R, SEARCH_R::ZERO_R, SEARCH_R::LARGE_R, 1>(),
                                         CalcForceLong<SEARCH_R::ZERO_R,  SEARCH_R::LARGE_R, SEARCH_R::ZERO_R, SEARCH_R::LARGE_R, 2>(),
                                         RetrieveLong(),
                                         tag_max, n_walk_limit,
                                         n_tot, theta, n_leaf_limit, n_group_limit,
                                         ic_mode, bc_mode, fr_mode, pos_org, test_result,
                                         epj_tmp);
                        PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"finish PS::LONG_NO_CUTOFF"<<std::endl;
                    }
                }

                if( static_cast<PS::SEARCH_MODE>(fr_mode) == PS::LONG_CUTOFF){
                    PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"start PS::LONG_CUTOFF"<<std::endl;
                    SetRsearchAve(system);
                    PS::TreeForForceLong<Force, EpAll, EpAll>::MonopoleWithCutoff tree;
                    tree.setExchagneLETMode(EX_LET_MODE);
                    EpAll epj_tmp;
                    CompareForceLong(system, tree, dinfo,
                                     CalcForceLong<SEARCH_R::EPJ_R,  SEARCH_R::EPJ_R, SEARCH_R::UNUSED, SEARCH_R::UNUSED, 0>(),
                                     CalcForceLong<SEARCH_R::EPI_R,  SEARCH_R::EPI_R, SEARCH_R::UNUSED, SEARCH_R::UNUSED, 0>(), 
                                     CalcForceLong<SEARCH_R::EPJ_R,  SEARCH_R::EPJ_R, SEARCH_R::EPI_R, SEARCH_R::EPI_R, 1>(),
                                     CalcForceLong<SEARCH_R::EPJ_R,  SEARCH_R::EPJ_R, SEARCH_R::EPI_R, SEARCH_R::EPI_R, 2>(),
                                     RetrieveLong(),
                                     tag_max, n_walk_limit,
                                     n_tot, theta, n_leaf_limit, n_group_limit,
                                     ic_mode, bc_mode, fr_mode, pos_org, test_result,
                                     epj_tmp);
                    PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"finish PS::LONG_CUTOFF"<<std::endl;
                }


                if( static_cast<PS::SEARCH_MODE>(fr_mode) == PS::LONG_SCATTER
                    && static_cast<PS::BOUNDARY_CONDITION>(bc_mode) == PS::BOUNDARY_CONDITION_OPEN){
                    PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"start PS::LONG_SCATTER"<<std::endl;
                    SetRsearchRand(system);
                    //PS::TreeForForce<PS::SEARCH_MODE_LONG_SCATTER, Force, EpAll, EpAll, MomentAll, MomentAll, SpAll> tree;
                    PS::TreeForForceLong<Force, EpSimple, EpAll>::MonopoleWithScatterSearch tree;
                    tree.setExchagneLETMode(EX_LET_MODE);
                    EpAll epj_tmp;
                    CompareForceLong(system, tree, dinfo,
                                     CalcForceLong<SEARCH_R::EPJ_R,  SEARCH_R::LARGE_R, SEARCH_R::UNUSED, SEARCH_R::UNUSED, 0>(),
                                     CalcForceLong<SEARCH_R::ZERO_R, SEARCH_R::LARGE_R, SEARCH_R::UNUSED, SEARCH_R::UNUSED, 0>(), 
                                     CalcForceLong<SEARCH_R::EPJ_R,  SEARCH_R::LARGE_R, SEARCH_R::ZERO_R, SEARCH_R::LARGE_R, 1>(),
                                     CalcForceLong<SEARCH_R::EPJ_R,  SEARCH_R::LARGE_R, SEARCH_R::ZERO_R, SEARCH_R::LARGE_R, 2>(),
                                     RetrieveLong(),
                                     tag_max, n_walk_limit,
                                     n_tot, theta, n_leaf_limit, n_group_limit,
                                     ic_mode, bc_mode, fr_mode, pos_org, test_result,
                                     epj_tmp);
                    PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"finish PS::LONG_SCATTER"<<std::endl;
                }
                if( static_cast<PS::SEARCH_MODE>(fr_mode) == PS::LONG_SYMMETRY
                    && static_cast<PS::BOUNDARY_CONDITION>(bc_mode) == PS::BOUNDARY_CONDITION_OPEN){
                    PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"start PS::LONG_SYMMETRY"<<std::endl;
                    SetRsearchRand(system);
                    PS::TreeForForceLong<Force, EpAll, EpAll>::MonopoleWithSymmetrySearch tree;
                    tree.setExchagneLETMode(EX_LET_MODE);
                    EpAll epj_tmp;
                    CompareForceLong(system, tree, dinfo,
                                     CalcForceLong<SEARCH_R::MAX_R,  SEARCH_R::LARGE_R, SEARCH_R::UNUSED, SEARCH_R::UNUSED, 0>(),
                                     CalcForceLong<SEARCH_R::ZERO_R, SEARCH_R::LARGE_R, SEARCH_R::UNUSED, SEARCH_R::UNUSED, 0>(), 
                                     CalcForceLong<SEARCH_R::MAX_R,  SEARCH_R::LARGE_R, SEARCH_R::ZERO_R, SEARCH_R::LARGE_R, 1>(),
                                     CalcForceLong<SEARCH_R::MAX_R,  SEARCH_R::LARGE_R, SEARCH_R::ZERO_R, SEARCH_R::LARGE_R, 2>(),
                                     RetrieveLong(),
                                     tag_max, n_walk_limit,
                                     n_tot, theta, n_leaf_limit, n_group_limit,
                                     ic_mode, bc_mode, fr_mode, pos_org, test_result,
                                     epj_tmp);
                    PS::Comm::barrier(); if(PS::Comm::getRank()==0) std::cerr<<"finish PS::LONG_SYMMETRY"<<std::endl;
                }

            }
        }
    }
    if(PS::Comm::getRank() == 0){
        for(auto tr : test_result){
            //if(!tr.pass) tr.dump();
            if(1) tr.dump();
        }
    }
    
    PS::Finalize();
    return 0;
}
