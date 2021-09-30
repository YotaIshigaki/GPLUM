#include<particle_simulator.hpp>

////////////////////
//// LONG FORCE ////

////////////////////////////////////
// FORCE AND PARTICLE FOR GRAVITY //
class ForceGrav{
public:
    PS::F64vec acc;
    PS::F64 pot;
    PS::F64vec getAcc() const {return acc;}
    PS::F64 getPot() const {return pot;}
    void clear(){
        acc = 0.0;
        pot = 0.0;
    }
};

class FPGrav{
public:
    PS::F64vec pos;
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pot;
    static PS::F64 r_cutoff;
    static PS::F64 eps;
    PS::F64vec getPos() const {return pos;}    
    PS::S64 getID() const {return id;}
    PS::F64 getRSearch() const {return r_cutoff;}
    PS::F64 getCharge() const {return mass;}
    void copyFromForce(const ForceGrav & force){
        this->acc = force.acc;
        this->pot = force.pot;
    }
    void dump(std::ostream & fout=std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"vel="<<vel<<std::endl;
        fout<<"acc="<<acc<<std::endl;
        fout<<"pot="<<pot<<std::endl;
    }
};

PS::F64 FPGrav::r_cutoff;
PS::F64 FPGrav::eps;
//PS::F64 FPGrav::r_cutoff = 1.0/ 16.0;

////////////////////////
/// GRAVITY NO CUTOFF///
class EPIGrav{
public:
    PS::S64 id;
    PS::F64vec pos;
    static PS::F64 eps;
    PS::F64vec getPos() const {return pos;}
    PS::S64 getID() const {return id;}
    void copyFromFP(const FPGrav & fp){ 
        pos = fp.pos;
        id = fp.id;
    }
};

PS::F64 EPIGrav::eps;
//PS::F64 EPIGrav::eps = 1.0/32.0;

class EPJGrav{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    void copyFromFP(const FPGrav & fp){ 
        mass = fp.mass;
        pos = fp.pos;
        id = fp.id;
    }
    PS::F64vec getPos() const {return pos;}
    void setPos(const PS::F64vec & p){
        pos = p;
    }
    PS::F64 getCharge() const {return mass;}
    PS::S64 getID() const {return id;}
    void dump(std::ostream & fout = std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
    }
};

template<class Tepi, class Tepj, class Tforce>
struct CalcForceEpEp{
    void operator () (const Tepi * ep_i,
                      const PS::S32 n_ip,
                      const Tepj * ep_j,
                      const PS::S32 n_jp,
                      Tforce * force){
        PS::F64 eps2 = Tepi::eps * Tepi::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].getPos();
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            PS::S64 idi = ep_i[i].getID();
            for(PS::S32 j=0; j<n_jp; j++){
		PS::S64 idj = ep_j[j].getID();
                if( idi == idj ) continue;
                PS::F64vec rij = xi - ep_j[j].getPos();
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= ep_j[j].getCharge();
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};

template<class Tepi, class Tspj, class Tforce>
struct CalcForceSpEp{
    void operator () (const Tepi * ep_i,
                      const PS::S32 n_ip,
                      const Tspj * sp_j,
                      const PS::S32 n_jp,
                      Tforce * force){
        PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].getPos();
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                PS::F64vec rij = xi - sp_j[j].getPos();
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= sp_j[j].getCharge();
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};


//////////////////////////
/// GRAVITY WITH CUTOFF///
class EPIGravCutoff{
    PS::S64 id;
    PS::F64vec pos;
public:
    static PS::F64 eps;
    static PS::F64 r_cutoff;
    PS::F64vec getPos() const {return pos;}
    PS::S64 getID() const {return id;}
    void copyFromFP(const FPGrav & fp){
	id = fp.getID();
        pos = fp.getPos();
	r_cutoff = fp.getRSearch();
    }
};

//PS::F64 EPIGravCutoff::eps = 1.0/32.0;
PS::F64 EPIGravCutoff::eps;
PS::F64 EPIGravCutoff::r_cutoff;

class EPJGravCutoff{
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
public:
    static PS::F64 r_cutoff;
    void copyFromFP(const FPGrav & fp){ 
        mass = fp.getCharge();
        pos = fp.getPos();
        id = fp.getID();
	r_cutoff = fp.getRSearch();
    }
    PS::F64vec getPos() const { return pos;}
    void setPos(const PS::F64vec & p){ pos = p;}
    PS::F64 getCharge() const {return mass;}
    PS::F64 getRSearch() const {return r_cutoff;}
    PS::S64 getID() const {return id;}
    void dump(std::ostream & fout = std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
    }
};

PS::F64 EPJGravCutoff::r_cutoff;


template<class Tepi, class Tepj, class Tforce>
struct CalcForceEpEpCutoff{
    void operator () (const Tepi * ep_i,
                      const PS::S32 n_ip,
                      const Tepj * ep_j,
                      const PS::S32 n_jp,
                      Tforce * force){
        PS::F64 eps2 = Tepi::eps * Tepi::eps;
	PS::F64 r_cutoff_sq = Tepi::r_cutoff * Tepi::r_cutoff;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].getPos();
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            PS::S64 idi = ep_i[i].getID();
            for(PS::S32 j=0; j<n_jp; j++){
		PS::S64 idj = ep_j[j].getID();
                if( idi == idj ) continue;
                PS::F64vec rij = xi - ep_j[j].getPos();
		if(rij*rij > r_cutoff_sq) continue;
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= ep_j[j].getCharge();
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};

template<class Tepi, class Tspj, class Tforce>
struct CalcForceSpEpCutoff{
    void operator () (const Tepi * ep_i,
                      const PS::S32 n_ip,
                      const Tspj * sp_j,
                      const PS::S32 n_jp,
                      Tforce * force){
	PS::F64 eps2 = Tepi::eps * Tepi::eps;
	PS::F64 r_cutoff_sq = Tepi::r_cutoff * Tepi::r_cutoff;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].getPos();
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                PS::F64vec rij = xi - sp_j[j].getPos();
		if(rij*rij > r_cutoff_sq) continue;
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= sp_j[j].getCharge();
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};


template<class Tforce>
struct CompareGrav{
    void operator () (Tforce * grav_tree, 
                      Tforce * grav_direct, 
                      const PS::S32 n, 
                      std::ostream & fout){
	PS::S32 n_err = 0;
        bool err = false;
        for(PS::S32 i=0; i<n; i++){
            PS::F64 dpot = std::abs( (grav_tree[i].pot - grav_direct[i].pot) / grav_tree[i].pot );
            PS::F64vec dacc_vec = grav_tree[i].acc - grav_direct[i].acc;
            PS::F64 dacc = sqrt( (dacc_vec*dacc_vec) / (grav_tree[i].acc*grav_tree[i].acc) );
            //if( dpot > 1e-1 || dacc > 1e-1){
	    if( dpot > 1e-1){
		if(grav_tree[i].pot < grav_direct[i].pot){
		    fout<<"Compare Grav: FAIL"<<std::endl;
		    fout<<"grav_tree[i].pot="<<grav_tree[i].pot<<" grav_direct[i].pot="<<grav_direct[i].pot<<std::endl;
		    fout<<"grav_tree[i].acc="<<grav_tree[i].acc<<" grav_direct[i].acc="<<grav_direct[i].acc<<std::endl;
		}
                err = true;
		n_err++;
            }
        }
	if(err){
	    fout<<"n_err="<<n_err<<" n_err/n="<<(double)(n_err)/n<<std::endl;
	}
	else{
	    fout<<"Compare Grav: PASS"<<std::endl;
	}
    }
};


/////////////////////
//// SHORT FORCE ////


///////////
/// SPH ///

/////////////////////////////
/// FORCE & FULL PARTICLE ///

class ResultDens{
public:
    void clear(){
        dens = r_search_next = kernel_length = 0.0;
        n_neighbour = 0;
    }
    PS::F64 dens;
    PS::F64 r_search_next;
    PS::S32 n_neighbour;
    PS::F64 kernel_length;
};

class FPSPH{
public:
    void copyFromForce(const ResultDens & dens){
        this->dens = dens.dens;
        this->kernel_length = dens.kernel_length;
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    PS::F64 getRSearch() const {
        return this->r_search;
    }

    PS::F64 mass;
    PS::F64vec pos;
    PS::F64 r_search;
    PS::F64 kernel_length;
    PS::F64 dens;
    PS::F64 divv;
    PS::F64 pres;
    PS::F64 vel_sound;
    PS::F64vec acc;
    PS::F64 engdot;
    PS::S64 id;
    PS::F64vec vel;
    PS::F64 eng;
    PS::F64vec vel_half;
    PS::F64 eng_half;
    PS::S64 n_neighbour;
    void dump(std::ostream & fout=std::cout) const {
        fout<<"pos="<<pos<<std::endl;
        fout<<"dens="<<dens<<std::endl;
        fout<<"kernel_length="<<kernel_length<<std::endl;
    }
};

/////////////////////
//// SCATER MODE ////

class EPIScatter{
public:
    PS::S64 id;
    PS::F64vec pos;
    void copyFromFP(const FPSPH & fp){
        this->id = fp.id;
        this->pos = fp.getPos();
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
};

class EPJScatter{
public:
    enum{
        n_neighbour_crit = 7,
    };
    PS::F64 mass;
    PS::F64 r_search;
    PS::F64vec pos;
    PS::S64 id;
    PS::S64 getIndex() const { return this->id; }
    PS::F64 getCharge() const { return this->mass; }
    PS::F64vec getPos() const { return this->pos; }
    PS::F64 getRSearch() const { return this->r_search; }
    void copyFromFP(const FPSPH & fp){
        this->mass = fp.getCharge();
        this->r_search = fp.getRSearch();
        this->pos = fp.getPos();
        this->id = fp.id;
    }
    void setPos(const PS::F64vec & p){this->pos = p;}
    void dump(std::ostream & fout=std::cout) const {
        fout<<"id="<<id<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"r_search="<<r_search<<std::endl;
    }
};

inline PS::F64 CubicSpline(const PS::F64 r_sq,
			   const PS::F64 h_inv){
    const PS::F64 xi = sqrt(r_sq) * h_inv;
    const PS::F64 xi10 = (1.0-xi > 0.0) ? 1.0-xi : 0.0;
    const PS::F64 xi05 = (0.5-xi > 0.0) ? 0.5-xi : 0.0;
    return xi10*xi10*xi10 - 4.0*xi05*xi05*xi05;
}

inline PS::F64vec CubicSpline_ri(const PS::F64vec & rij,
				 const PS::F64 h_inv){
    PS::F64 r_sq = rij * rij;
    PS::F64 r = sqrt(r_sq);
    PS::F64 xi = r * h_inv;
    PS::F64 xi10 = (1.0-xi > 0.0) ? 1.0-xi : 0.0;
    PS::F64 xi05 = (0.5-xi > 0.0) ? 0.5-xi : 0.0;
    PS::F64 C = (3.0*xi10*xi10 - 12.0*xi05*xi05) * h_inv / r;
    return -C * rij;
}



struct CalcDensityScatter{
    void operator () (const EPIScatter * ep_i,
                      const PS::S32 n_ip,
                      const EPJScatter * ep_j,
                      const PS::S32 n_jp,
                      ResultDens * dens){
        //static PS::F64 Cnorm = 8.0/3.0; // 1D
        static PS::F64 Cnorm = 80.0/(7.0*M_PI); // 2D
        //static PS::F64 Cnorm = 16.0/M_PI; // 3D
        for(PS::S32 i=0; i<n_ip; i++){
            for(PS::S32 j=0; j<n_jp; j++){
                const PS::F64 r_crit_sq = ep_j[j].getRSearch() * ep_j[j].getRSearch();
                const PS::F64vec dr = ep_j[j].getPos() - ep_i[i].getPos();
                const PS::F64 dr_sq = dr * dr;
                if( r_crit_sq > dr_sq ){
                    const PS::F64 h_inv = 1.0 / ep_j[j].getRSearch();
		    //const PS::F64 h_inv_qb = h_inv * h_inv * h_inv;
		    const PS::F64 h_inv_sq = h_inv * h_inv;
                    dens[i].dens += ep_j[j].getCharge() * CubicSpline(dr_sq, h_inv) * h_inv_sq;
                }
            }
            dens[i].dens *= Cnorm;
        }
    }
};

/////////////////////
//// GATHER MODE ////

class EPIGather{
public:
    enum{
        n_neighbour_crit = 7,
    };
    PS::F64vec pos;
    PS::F64 r_search;
    PS::S64 id;
    PS::S64 getIndex() const { return this->id;}
    PS::F64vec getPos() const {
        return this->pos;
    }
    void setPos(const PS::F64vec & p){
        this->pos = p;
    }
    PS::F64 getRSearch() const {
        return this->r_search;
    }
    void copyFromFP(const FPSPH & fp){
        this->pos = fp.getPos();
        this->r_search = fp.getRSearch();
        id = fp.id;
    }
};

class EPJGather{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    void copyFromFP(const FPSPH & fp){
        this->mass = fp.getCharge();
        this->pos = fp.getPos();
        this->id = fp.id;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    void setPos(const PS::F64vec & p){
        this->pos = p;
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
};

struct CalcDensityGather{
    void operator () (const EPIGather * ep_i,
                      const PS::S32 n_ip,
                      const EPJGather * ep_j,
                      const PS::S32 n_jp,
                      ResultDens * dens){
        //static PS::F64 Cnorm = 8.0/3.0; // 1D
        static PS::F64 Cnorm = 80.0/(7.0*M_PI); // 2D
        //static PS::F64 Cnorm = 16.0/M_PI; // 3D
	//std::cerr<<"check"<<std::endl;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64 h_inv = 1.0 / ep_i[i].getRSearch();
            const PS::F64 r_crit_sq = ep_i[i].getRSearch() * ep_i[i].getRSearch();
	    //const PS::F64 h_inv_qb = h_inv * h_inv * h_inv; // for 3D
	    const PS::F64 h_inv_sq = h_inv * h_inv; // for 3D
            for(PS::S32 j=0; j<n_jp; j++){
                const PS::F64vec dr = ep_j[j].getPos() - ep_i[i].getPos();
                const PS::F64 dr_sq = dr * dr;
                if( r_crit_sq > dr_sq ){
                    dens[i].dens += ep_j[j].getCharge() * CubicSpline(dr_sq, h_inv);
                }
            }
	    //std::cerr<<"dens[i].dens0="<<dens[i].dens<<std::endl;
            dens[i].dens *= Cnorm * h_inv_sq;
	    //std::cerr<<"Cnorm="<<Cnorm<<std::endl;
	    //std::cerr<<"dens[i].dens1="<<dens[i].dens<<std::endl;
        }
    }
};

///////////////////////
//// SYMMETRY MODE ////

class EPISymmetry{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64 r_search;
    PS::F64vec pos;
    void copyFromFP(const FPSPH & fp){
        this->id = fp.id;
        this->mass = fp.getCharge();
        this->r_search = fp.getRSearch();
        this->pos = fp.getPos();
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    PS::F64 getRSearch() const {
        return this->r_search;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    void setPos(const PS::F32vec & p){
        this->pos = p;
    }
};
class EPJSymmetry{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64 r_search;
    PS::F64vec pos;
    void copyFromFP(const FPSPH & fp){
        this->id = fp.id;
        this->mass = fp.getCharge();
        this->r_search = fp.getRSearch();
        this->pos = fp.getPos();
    }
    PS::F64 getCharge() const {
        return this->mass;
    }
    PS::F64 getRSearch() const {
        return this->r_search;
    }
    PS::F64vec getPos() const {
        return this->pos;
    }
    void setPos(const PS::F32vec & p){
        this->pos = p;
    }
};
struct CalcDensitySymmetry{
    void operator () (const EPISymmetry * ep_i,
                      const PS::S32 n_ip,
                      const EPJSymmetry * ep_j,
                      const PS::S32 n_jp,
                      ResultDens * dens){
        //static PS::F64 Cnorm = 8.0/3.0; // 1D
        static PS::F64 Cnorm = 80.0/(7.0*M_PI); // 2D
        //static PS::F64 Cnorm = 16.0/M_PI; // 3D
        for(PS::S32 i=0; i<n_ip; i++){
            for(PS::S32 j=0; j<n_jp; j++){
                const PS::F64vec dr = ep_j[j].getPos() - ep_i[i].getPos();
                const PS::F64 dr_sq = dr * dr;
		const PS::F64 hij = (ep_i[i].getRSearch() + ep_j[j].getRSearch()) * 0.5;
		const PS::F64 h_inv = 1.0 / hij;
		//const PS::F64 h_inv_qb = h_inv * h_inv * h_inv;
		const PS::F64 h_inv_sq = h_inv * h_inv;
                if( hij*hij > dr_sq ){
                    //dens[i].dens += ep_j[j].getCharge() * CubicSpline(dr_sq, h_inv) * h_inv_qb;
		    dens[i].dens += ep_j[j].getCharge() * CubicSpline(dr_sq, h_inv) * h_inv_sq;
                }
            }
            dens[i].dens *= Cnorm;
        }
    }
};


///////////////////////////
//// COMPARE FUNCTIONS ////

struct CompareDensity{
    void operator () (ResultDens * dens0, ResultDens * dens1, 
                      const PS::S32 n, std::ostream & fout){
        bool err = false;
        for(PS::S32 i=0; i<n; i++){
            if( std::abs(dens0[i].dens - dens1[i].dens) > 1e-5){
                fout<<"CompareDensity: FAIL"<<std::endl;
                fout<<"desn0[i].dens="<<dens0[i].dens<<" dens1[i].dens="<<dens1[i].dens<<std::endl;
                err = true;
            }
        }
        if(!err) fout<<"CompareDensity: PASS"<<std::endl;
	else{
	    fout<<"CompareDensity: PASS"<<std::endl;
	}
	for(PS::S32 i=0; i<10; i++){
	    fout<<"dens0[i].dens="<<dens0[i].dens<<" dens1[i].dens="<<dens1[i].dens<<std::endl;
	}
    }
};



/////////////////////
/// I/O FUNCTIONS ///

template<class Tpsys>
void ReadNemoAscii(Tpsys & psys,
                   PS::S32 & n_glb,
                   PS::S32 & n_loc,  
                   PS::F32 & t_sys,
                   const char * ifile){
    std::ifstream finput;
    finput.open(ifile);
    assert(finput);
    PS::S32 dim;
    finput>>n_glb>>dim>>t_sys;
    std::cerr<<"ifile:"<<ifile<<std::endl;
    std::cerr<<"n_glb="<<n_glb<<std::endl;
    std::cerr<<"dim="<<dim<<std::endl;
    std::cerr<<"t_sys="<<t_sys<<std::endl;

    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    n_loc = n_glb/n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;

    psys.createParticle((n_glb/n_proc)*4+1000);
    psys.setNumberOfParticleLocal(n_loc);

    //PS::F32vec pos_shift(10.0, 20.0, 30.0); // for debug
    PS::F32vec pos_shift(0.0);

    PS::S32 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
    const PS::S32 i_t = i_h+n_loc;
    PS::F32 xf32;
    PS::F32vec vf32;

    for(PS::S32 i=i_h, n=0; i<i_t; i++, n++)psys[n].id = i;

    for(PS::S32 i=0; i<i_h; i++) finput>>xf32;
    for(PS::S32 i=i_h, n=0; i<i_t; i++, n++)finput>>psys[n].mass;
    for(PS::S32 i=i_t; i<n_glb; i++) finput>>xf32;

    for(PS::S32 i=0; i<i_h; i++){
	finput>>xf32;
	finput>>xf32;
	finput>>xf32;
    }
    for(PS::S32 i=i_h, n=0; i<i_t; i++, n++){
	finput>>psys[n].pos.x;
	finput>>psys[n].pos.y;
	finput>>xf32;
    }
    for(PS::S32 i=i_t; i<n_glb; i++){
	finput>>xf32;
	finput>>xf32;
	finput>>xf32;
    }
    for(PS::S32 i=0; i<i_h; i++){
	finput>>xf32;
	finput>>xf32;
	finput>>xf32;
    }
    for(PS::S32 i=i_h, n=0; i<i_t; i++, n++){
	finput>>psys[n].vel.x;
	finput>>psys[n].vel.y;
	finput>>xf32;
    }
    for(PS::S32 i=i_t; i<n_glb; i++){
	finput>>xf32;
	finput>>xf32;
	finput>>xf32;
    }
    finput.close();
}
