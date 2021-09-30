#include<iostream>
#include<fstream>
#include<unistd.h>
#include<particle_simulator.hpp>



class ForceGrav{
public:
    PS::F64vec acc;
    PS::F64 pot;
    void clear(){
        acc = 0.0;
        pot = 0.0;
    }
};

class FPGrav{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pot;
    PS::F64 r_search;
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & p) { pos = p; }
    void copyFromForce(const ForceGrav & force){
        acc = force.acc;
        pot = force.pot;
    }
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", 
                this->id, this->mass, this->pos.x, this->pos.y, this->pos.z, this->vel.x, this->vel.y, this->vel.z);
    }
    void readAscii(FILE* fp){
        fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
               &this->id, &this->mass, &this->pos.x, &this->pos.y, &this->pos.z, &this->vel.x, &this->vel.y, &this->vel.z);
    }
    void dump(std::ostream & fout){
	fout<<"id="<<id<<std::endl;
	fout<<"pos="<<pos<<std::endl;
	fout<<"r_search="<<r_search<<std::endl;
    }
};

class ForceRSearch{
public:
    PS::F64 r_cut;
    PS::F64 r_search;
    void clear(){
	r_search = 0.0;
	r_cut = 0.0;
    }
};

class EPIRSearch{
public:
    PS::S64 id;
    PS::F64 r_search;
    PS::F64vec pos;
    static const PS::S64 n_ngb_limit = 100;
    PS::F64vec getPos() const { return pos;}
    PS::F64 getRSearch() const { return r_search;}
    void copyFromFP(const FPGrav & fp){ 
        pos = fp.pos;
        id = fp.id;
	r_search = fp.r_search;
    }
};

class EPJRSearch{
public:
    PS::S64 id;
    PS::F64 r_search;
    PS::F64vec pos;
    void copyFromFP(const FPGrav & fp){ 
        pos = fp.pos;
        id = fp.id;
	r_search = fp.r_search;
    }
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & pos_new){ pos = pos_new;}
    PS::F64 getRSearch() const { return r_search;}
};

struct CalcRSearch{
    void operator () (const EPIRSearch * ep_i,
                      const PS::S32 n_ip,
                      const EPJRSearch * ep_j,
                      const PS::S32 n_jp,
                      ForceRSearch * force){
	const PS::S64 n_ngb_limit = EPIRSearch::n_ngb_limit;
	std::vector< std::pair<PS::F64, PS::S64> > dis2_id;
	for(PS::S32 i=0; i<n_ip; i++){
	    dis2_id.reserve(1000);
	    dis2_id.clear();
	    if(ep_i[i].r_search == 0.0) continue;
	    const PS::F64 r_search_sq = ep_i[i].r_search * ep_i[i].r_search;
	    force[i].clear();
	    for(PS::S32 j=0; j<n_jp; j++){
		PS::F64vec rij = ep_i[i].pos - ep_j[j].pos;
		PS::F64 dis2 = rij * rij;
		if( dis2 < r_search_sq ){
		    dis2_id.push_back( std::make_pair(dis2, ep_j[j].id) );
		}
	    }
	    if(dis2_id.size() >= n_ngb_limit){
		std::sort(dis2_id.begin(), dis2_id.end());
		force[i].r_cut = ( sqrt(dis2_id[n_ngb_limit-2].first) + sqrt(dis2_id[n_ngb_limit-1].first) ) * 0.5;
		force[i].r_search = 0.0;
		/*
		if(PS::Comm::getRank()==0){
		    std::cout<<"check a: dis2_id.size()="<<dis2_id.size()<<std::endl;
		    std::cout<<"force[i].r_cut="<<force[i].r_cut
			     <<" force[i].r_search="<<force[i].r_search<<std::endl;
		}
		*/
	    }
	    else{
		force[i].r_search = ep_i[i].r_search;
		/*
		if(PS::Comm::getRank()==0){
		    std::cout<<"check b: dis2_id.size()="<<dis2_id.size()<<std::endl;
		    std::cout<<"force[i].r_cut="<<force[i].r_cut
			     <<" force[i].r_search="<<force[i].r_search<<std::endl;
		}
		*/
	    }
	}
    }
};


template<class Tepi, class Tepj, class Tforce>
struct CalcForceDummy{
    void operator () (const Tepi * ep_i,
                      const PS::S32 n_ip,
                      const Tepj * ep_j,
                      const PS::S32 n_jp,
                      Tforce * force){}
};





class EPIGrav{
public:
    PS::S64 id;
    PS::F64 r_search;
    PS::F64vec pos;
    static PS::F64 eps;
    PS::F64vec getPos() const { return pos;}
    PS::F64 getRSearch() const { return r_search;}
    void copyFromFP(const FPGrav & fp){ 
        pos = fp.pos;
        id = fp.id;
	r_search = fp.r_search;
    }
};

PS::F64 EPIGrav::eps = 1e-3; // 1PC

class EPJGrav{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64 r_search;
    PS::F64vec pos;
    void copyFromFP(const FPGrav & fp){ 
        mass = fp.mass;
        pos = fp.pos;
        id = fp.id;
        r_search = fp.r_search;
    }
    PS::F64vec getPos() const { return pos; }
    PS::F64 getRSearch() const { return r_search;}
    void setPos(const PS::F64vec & pos_new){ pos = pos_new;}
    PS::F64 getCharge() const { return mass; }

};


#ifdef USEPHANTOMGRAPE
/*
struct CalcForceEpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const EPJGrav * ep_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){
        const PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        PhantomGrapeQuad pg;
        pg.set_eps2(eps2);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
        }
        for(PS::S32 i=0; i<n_jp; i++){
            const PS::F64 m_j = ep_j[i].getCharge();
            const PS::F64vec pos_j = ep_j[i].getPos();
            pg.set_epj_one(i, pos_j.x, pos_j.y, pos_j.z, m_j);
        }
        pg.run_epj(n_ip, n_jp);
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64 * p = &(force[i].pot);
            PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
            pg.get_accp_one(i, a[0], a[1], a[2], *p);
        }
    }
};
*/

struct CalcForceEpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const EPJGrav * ep_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){
        const PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        PhantomGrapeQuad pg;
        pg.set_eps2(eps2);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
        }
        PS::S32 loop_max = (n_jp-1) / PhantomGrapeQuad::NJMAX + 1;
        for(PS::S32 loop=0; loop<loop_max; loop++){
            const PS::S32 ih = PhantomGrapeQuad::NJMAX*loop;
            const PS::S32 n_jp_tmp = ( (n_jp - ih) < PhantomGrapeQuad::NJMAX) ? (n_jp - ih) : PhantomGrapeQuad::NJMAX;  
            const PS::S32 it =ih + n_jp_tmp;
            PS::S32 i_tmp = 0;
            for(PS::S32 i=ih; i<it; i++, i_tmp++){
                const PS::F64 m_j = ep_j[i].getCharge();
                const PS::F64vec pos_j = ep_j[i].getPos();
                pg.set_epj_one(i, pos_j.x, pos_j.y, pos_j.z, m_j);
            }
            pg.run_epj(n_ip, n_jp);
            for(PS::S32 i=0; i<n_ip; i++){
                PS::F64 * p = &(force[i].pot);
                PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
                pg.accum_accp_one(i, a[0], a[1], a[2], *p);
            }
        }
    }
};


#else
struct CalcForceEpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const EPJGrav * ep_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){

        PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            PS::S64 idi = ep_i[i].id;
            for(PS::S32 j=0; j<n_jp; j++){
                if( idi == ep_j[j].id ) continue;
                PS::F64vec rij = xi - ep_j[j].pos;
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= ep_j[j].mass;
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};
#endif // USEPHANTOMGRAPE



#ifdef USEPHANTOMGRAPE
#ifdef MONOPOLE
struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJMonopole * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){


        if(n_ip > PhantomGRAPE::PG_NIMAX){
            std::cout<<"n_ip(ep-sp)="<<n_ip<<std::endl;
            std::cout<<"PhantomGRAPE::PG_NIMAX(ep-sp)="<<PhantomGRAPE::PG_NIMAX<<std::endl;
        }
        if(n_jp > PhantomGRAPE::PG_NJMAX){
            std::cout<<"n_jp(ep-sp)="<<n_jp<<std::endl;
            std::cout<<"PhantomGRAPE::PG_NJMAX(ep-sp)="<<PhantomGRAPE::PG_NJMAX<<std::endl;
        }

        double xi[PhantomGRAPE::PG_NIMAX][3];
        double mxj[PhantomGRAPE::PG_NJMAX][4];
        double ai[PhantomGRAPE::PG_NIMAX][3];
        double pi[PhantomGRAPE::PG_NIMAX];
        const PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            xi[i][0] = ep_i[i].pos.x;
            xi[i][1] = ep_i[i].pos.y;
            xi[i][2] = ep_i[i].pos.z;
        }
        for(PS::S32 j=0; j<n_jp; j++){
            mxj[j][0] = sp_j[j].getPos().x;
            mxj[j][1] = sp_j[j].getPos().y;
            mxj[j][2] = sp_j[j].getPos().z;
            mxj[j][3] = sp_j[j].getCharge();
        }
        PhantomGRAPE pg;
        pg.set_eps2(eps2);
        pg.set_xj(n_jp, mxj);
        pg.set_xi(n_ip, xi);
        pg.run(n_ip, n_jp);
        pg.get_ai(n_ip, ai, pi);

        for(PS::S32 i=0; i<n_ip; i++){
            force[i].acc.x += ai[i][0];
            force[i].acc.y += ai[i][1];
            force[i].acc.z += ai[i][2];
            force[i].pot += pi[i];
        }

    }
};
#elif QUADRUPOLE
/*
struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJQuadrupole * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){

        const PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        PhantomGrapeQuad pg;
        pg.set_eps2(eps2);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
        }

        for(PS::S32 j=0; j<n_jp; j++){
            const PS::F64 m_j = sp_j[j].getCharge();
            const PS::F64vec pos_j = sp_j[j].getPos();
            const PS::F64mat q = sp_j[j].quad;
            pg.set_spj_one(j, pos_j.x, pos_j.y, pos_j.z, m_j,
                           q.xx, q.yy, q.zz, q.xy, q.yz, q.xz);
        }
        pg.run_spj(n_ip, n_jp);
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64 * p = &(force[i].pot);
            PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
            pg.accum_accp_one(i, a[0], a[1], a[2], *p);
        }
    }
};
*/

struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJQuadrupole * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){
        const PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        PhantomGrapeQuad pg;
        pg.set_eps2(eps2);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
        }
        PS::S32 loop_max = (n_jp-1) / PhantomGrapeQuad::NJMAX + 1;
        for(PS::S32 loop=0; loop<loop_max; loop++){
            const PS::S32 ih = PhantomGrapeQuad::NJMAX*loop;
            const PS::S32 n_jp_tmp = ( (n_jp - ih) < PhantomGrapeQuad::NJMAX) ? (n_jp - ih) : PhantomGrapeQuad::NJMAX;  
            const PS::S32 it = ih + n_jp_tmp;
            PS::S32 i_tmp = 0;
            for(PS::S32 i=ih; i<it; i++, i_tmp++){
                const PS::F64 m_j = sp_j[i].getCharge();
                const PS::F64vec pos_j = sp_j[i].getPos();
                const PS::F64mat q = sp_j[i].quad;
                pg.set_spj_one(i, pos_j.x, pos_j.y, pos_j.z, m_j,
                               q.xx, q.yy, q.zz, q.xy, q.yz, q.xz);
            }
            pg.run_spj(n_ip, n_jp_tmp);
            for(PS::S32 i=0; i<n_ip; i++){
                PS::F64 * p = &(force[i].pot);
                PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
                pg.accum_accp_one(i, a[0], a[1], a[2], *p);
            }
        }
    }
};

#endif // MONOPOLE QUADRUPOLE
#else // USEPHANTOMGRAPE

#ifdef MONOPOLE
struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJMonopole * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){
        PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                PS::F64vec rij = xi - sp_j[j].pos;
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= sp_j[j].mass;
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};
#elif QUADRUPOLE
struct CalcForceSpEp{
    void operator () (const EPIGrav * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJQuadrupole * sp_j,
                      const PS::S32 n_jp,
                      ForceGrav * force){



        PS::F64 eps2 = EPIGrav::eps * EPIGrav::eps;
        for(PS::S32 ip=0; ip<n_ip; ip++){
            PS::F64vec xi = ep_i[ip].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 jp=0; jp<n_jp; jp++){
                PS::F64 mj = sp_j[jp].mass;
                PS::F64vec xj= sp_j[jp].pos;
                PS::F64vec rij= xi - xj;
                PS::F64 r2 = rij * rij + eps2;
                PS::F64mat qj = sp_j[jp].quad;
                PS::F64 tr = qj.getTrace();
                PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                               (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                               (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
                PS::F64 qrr = qr * rij;
                PS::F64 r_inv = 1.0f/sqrt(r2);
                PS::F64 r2_inv = r_inv * r_inv;
                PS::F64 r3_inv = r2_inv * r_inv;
                PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
                PS::F64 qrr_r5 = r5_inv * qrr;
                PS::F64 qrr_r7 = r2_inv * qrr_r5;
                PS::F64 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
                PS::F64 B = -2.0*r5_inv;
                ai -= A*rij + B*qr;
                poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
            }
            force[ip].acc += ai;
            force[ip].pot += poti;
        }
    }
};
#endif // MONOPOLE
#endif // USEPHANTOMGRAPE
