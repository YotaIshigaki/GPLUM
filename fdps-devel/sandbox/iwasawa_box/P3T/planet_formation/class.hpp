#ifdef INTRINSIC_K
#define INTRINSIC
#include"phantomquad_for_p3t_k.hpp"
#endif
#ifdef INTRINSIC_X86
#define INTRINSIC
#include"phantomquad_for_p3t_x86.hpp"
#endif

class ForceSoft{
public:
    PS::F64vec acc;
    PS::F64 pot; // soft
    PS::S32 n_ngb;
    //PS::S32 n_ep;
    //PS::S32 n_sp;
    //PS::F64vec acc_pla; // for debug
    //PS::S32 id_ngb;
    //PS::F64vec acc_hard; // for debug
    void clear(){
        acc = 0.0;
        pot = 0.0;
        n_ngb = 0;
	//n_ep = 0;
	//n_sp = 0;
	//acc_pla = 0.0;
        //id_ngb = -1;
        //acc_hard = 0.0;
    }
};

class FPSoft{
public:
    PS::S64 id;
    PS::S64 adr;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc; // soft
    PS::F64 pot_tot; // soft + hard
    //PS::F64vec acc_hard; // for debug
    PS::S64 id_proc;
    PS::S32 n_ngb;
    PS::F64vec acc_pla; // for debug
    //PS::S32 n_ep;
    //PS::S32 n_sp;
    //PS::S32 n_ep_direct;
    //PS::S32 n_sp_direct;
    //PS::F64vec acc_direct;
    //PS::F64 pot_tot_direct; // soft + hard

    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & p) { pos = p; }
    void copyFromForce(const ForceSoft & force){
        acc = force.acc;
        pot_tot = force.pot;
	n_ngb = force.n_ngb;
	//n_ep = force.n_ep;
	//n_sp = force.n_sp;
        //acc_hard = force.acc_hard; // for debug
    }

    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %d \n", 
                this->id, this->mass, 
		this->pos.x, this->pos.y, this->pos.z,  // 3-5
		this->vel.x, this->vel.y, this->vel.z,  // 6-8
		this->acc.x, this->acc.y, this->acc.z,  // 9-11
		this->acc_pla.x, this->acc_pla.y, this->acc_pla.z, // 12-14
		this->pot_tot, this->n_ngb);
    }
    void readAscii(FILE* fp) {
        fscanf(fp, "%lld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d \n", 
                &this->id, &this->mass, 
		&this->pos.x, &this->pos.y, &this->pos.z,  // 3-5
		&this->vel.x, &this->vel.y, &this->vel.z,  // 6-8
		&this->acc.x, &this->acc.y, &this->acc.z,  // 9-11
		&this->acc_pla.x, &this->acc_pla.y, &this->acc_pla.z, // 12-14
		&this->pot_tot, &this->n_ngb);
    }

    /*
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", 
                this->id, this->mass, this->pos.x, this->pos.y, this->pos.z, this->vel.x, this->vel.y, this->vel.z);
    }
    */
    /*
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %d %d %d \n", 
                this->id, this->mass, 
		this->pos.x, this->pos.y, this->pos.z,  // 3-5
		this->vel.x, this->vel.y, this->vel.z,  // 6-8
		this->acc.x, this->acc.y, this->acc.z,  // 9-11
		this->acc_pla.x, this->acc_pla.y, this->acc_pla.z, // 12-14
		this->acc_direct.x, this->acc_direct.y, this->acc_direct.z, //15-17
		this->pot_tot, this->n_ngb, this->n_ep, this->n_sp);
    }
    void readAscii(FILE* fp) {
        fscanf(fp, "%lld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d \n", 
                &this->id, &this->mass, 
		&this->pos.x, &this->pos.y, &this->pos.z,  // 3-5
		&this->vel.x, &this->vel.y, &this->vel.z,  // 6-8
		&this->acc.x, &this->acc.y, &this->acc.z,  // 9-11
		&this->acc_pla.x, &this->acc_pla.y, &this->acc_pla.z, // 12-14
		&this->acc_direct.x, &this->acc_direct.y, &this->acc_direct.z, //15-17
		&this->pot_tot, &this->n_ngb, &this->n_ep, &this->n_sp);
    }
    */
    /*
    void readAscii(FILE* fp){
        fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
               &this->id, &this->mass, &this->pos.x, &this->pos.y, &this->pos.z, &this->vel.x, &this->vel.y, &this->vel.z);
    }
    */
    /*
    void readAscii(FILE* fp){
        fscanf(fp, "%lld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", 
               &this->id, &this->mass, &this->pos.x, &this->pos.y, &this->pos.z, &this->vel.x, &this->vel.y, &this->vel.z, &this->acc.x, &this->acc.y, &this->acc.z, &this->pot_tot, &this->n_ngb);
    }
    */
    /*
    void readAscii(FILE* fp){
        fscanf(fp, "%lld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", 
               &this->id, &this->mass, 
	       &this->pos.x, &this->pos.y, &this->pos.z, 
	       &this->vel.x, &this->vel.y, &this->vel.z, 
	       &this->acc.x, &this->acc.y, &this->acc.z, 
	       &this->acc_pla.x, &this->acc_pla.y, &this->acc_pla.z, 
	       &this->pot_tot, &this->n_ngb);
    }
    */



    void dump(std::ofstream & fout){
	fout<<"id= "<<id<<std::endl;
	fout<<"adr= "<<adr<<std::endl;
	fout<<"mass= "<<mass<<std::endl;
	fout<<"pos= "<<pos<<std::endl;
	fout<<"vel= "<<vel<<std::endl;
	fout<<"acc= "<<acc<<std::endl;
	fout<<"pot_tot= "<<pot_tot<<std::endl;
    }
};

class EPISoft{
public:
    PS::S64 id;
    PS::F64vec pos;
    static PS::F64 eps;
    PS::F64vec getPos() const { return pos;}
    void copyFromFP(const FPSoft & fp){ 
        pos = fp.pos;
        id = fp.id;
    }
    void dump(std::ostream & fout=std::cout) const {
	fout<<"id="<<id<<std::endl;
	fout<<"pos="<<pos<<std::endl;
	fout<<"eps="<<eps<<std::endl;
    }
};

PS::F64 EPISoft::eps = 1.0/1024.0;

class EPJSoft{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::S64 id_proc;
    static PS::F64 r_search;
    void copyFromFP(const FPSoft & fp){
        mass = fp.mass;
        pos = fp.pos;
        id = fp.id;
        vel = fp.vel;
	id_proc = fp.id_proc;
    }
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & pos_new){ pos = pos_new;}
    PS::F64 getCharge() const { return mass; }
    PS::F64 getRSearch() const { return r_search; }
    // FORDEBUG

    void dump(std::ostream & fout=std::cout) const {
	fout<<"id="<<id<<std::endl;
	fout<<"id_proc="<<id_proc<<std::endl;
	fout<<"mass="<<mass<<std::endl;
	fout<<"pos="<<pos<<std::endl;
	fout<<"vel="<<vel<<std::endl;
    }

#ifdef TMP2
    EPJSoft(){
	id = -1;
	mass = 0.0;
	pos = 0.0;
	vel = 0.0;
	id_proc = -1;
    }
#endif
#ifdef TMP3
    EPJSoft(){
	id = -100;
	mass = 100000.0;
	pos = 200000.0;
	vel = 300000.0;
	id_proc = -400;
    }
#endif

    void clear(){
	mass = 0.0;
	pos = vel = 0.0;
	id = id_proc = -1;
    }
    /*
    void messedup(){
	mass = 1000.0;
	pos = 10000.0;;
	vel = 0.0;
	id = id_proc = -10;
    }
    void messedup2(){
	mass = 1000.0;
	pos = 100.0;
	vel = 0.0;
	id = id_proc = -3;
    }
    void messedup3(){
	mass = 1000000.0;
	pos = 10000.0;
	vel = 0.0;
	id = id_proc = -4;
    }
    void messedup4(){
	mass = 1000000.0;
	pos = 0.0;
	vel = 10000.0;
	id = id_proc = -5;
    }
    */
    // FORDEBUG
};

PS::F64 EPJSoft::r_search;

#ifdef INTRINSIC

struct CalcForceEPEP{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 r_crit2 = EPJSoft::r_search * EPJSoft::r_search;
        static __thread PhantomGrapeQuad pg;
        pg.set_eps2(eps2);
        //pg.set_r_crit2(r_crit2);
	pg.set_r_crit2(r_crit2+eps2);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
#ifdef INTRINSIC_K
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
#else
            pg.set_xi_one_d(i, pos_i.x, pos_i.y, pos_i.z);
	    //pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
#endif
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
#ifdef INTRINSIC_K
		pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
#else
		pg.set_epj_one_d(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
		//pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
#endif
	    }
	    pg.run_epj_64bit_for_p3t(n_ip, n_jp_tmp);
            //pg.run_epj_for_p3t(n_ip, n_jp_tmp);
	    for(PS::S32 i=0; i<n_ip; i++){
		PS::F64 * p = &(force[i].pot);
		PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
		PS::F64 n_ngb = 0; 
		pg.accum_accp_one_d(i, a[0], a[1], a[2], *p, n_ngb);
		//pg.accum_accp_one(i, a[0], a[1], a[2], *p, n_ngb);
		force[i].n_ngb += (PS::S32)(n_ngb*1.0000001);
	    }
	}
        for(PS::S32 i=0; i<n_ip; i++){
	    force[i].n_ngb--;
	}
    }
};

struct CalcForceEPSP{
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      //const PS::SPJMonopoleScatter * sp_j,
		      const Tsp * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        static __thread PhantomGrapeQuad pg;
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
		pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
	    }
            pg.run_epj(n_ip, n_jp_tmp);
	    for(PS::S32 i=0; i<n_ip; i++){
		PS::F64 * p = &(force[i].pot);
		PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
		pg.accum_accp_one(i, a[0], a[1], a[2], *p);
	    }
	}
    }
};

#define CLAC_SP_64bit

struct CalcForceEPSPQuad{
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      //const PS::SPJMonopoleScatter * sp_j,
		      const Tsp * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        static __thread PhantomGrapeQuad pg;
        pg.set_eps2(eps2);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
#ifdef CLAC_SP_64bit
            pg.set_xi_one_d(i, pos_i.x, pos_i.y, pos_i.z);
#else
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
#endif
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
#ifdef CLAC_SP_64bit
                pg.set_spj_one_d(i, pos_j.x, pos_j.y, pos_j.z, m_j,
				 q.xx, q.yy, q.zz, q.xy, q.yz, q.xz);
#else
                pg.set_spj_one(i, pos_j.x, pos_j.y, pos_j.z, m_j,
                               q.xx, q.yy, q.zz, q.xy, q.yz, q.xz);
#endif
	    }
#ifdef CLAC_SP_64bit
            pg.run_spj_d(n_ip, n_jp_tmp);
#else
            pg.run_spj(n_ip, n_jp_tmp);
#endif
	    for(PS::S32 i=0; i<n_ip; i++){
		PS::F64 * p = &(force[i].pot);
		PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
#ifdef CLAC_SP_64bit
		pg.accum_accp_one_d(i, a[0], a[1], a[2], *p);
#else
		pg.accum_accp_one(i, a[0], a[1], a[2], *p);
#endif
	    }
	}
    }
};

#else //INTRINSIC

struct CalcForceEPEP{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 r_crit2 = EPJSoft::r_search * EPJSoft::r_search;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec xi = ep_i[i].pos;
            PS::S64 id_i = ep_i[i].id;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            PS::S32 & n_ngb_i = force[i].n_ngb;
	    //force[i].n_ep = n_jp;
            for(PS::S32 j=0; j<n_jp; j++){
                if(id_i == ep_j[j].id) continue;
                const PS::F64vec rij = xi - ep_j[j].pos;
                const PS::F64 r2 = rij * rij;
                const PS::F64 r2_eps = rij * rij + eps2;
                if(r2 < r_crit2){
                    if(id_i != ep_j[j].id){
                        n_ngb_i++;
                    }
                    continue;
                }
                const PS::F64 r_inv = 1.0/sqrt(r2_eps);
                const PS::F64 m_r = ep_j[j].mass * r_inv;
                const PS::F64 m_r3 = m_r * r_inv * r_inv;
                ai -= m_r3 * rij;
                poti -= m_r;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};


struct CalcForceEPSP{
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
		      const Tsp * sp_j,
                      //const PS::SPJMonopoleScatter * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
	    //force[i].n_sp = n_jp;
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


struct CalcForceEPSPQuad{
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
		      const Tsp * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
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




struct CalcForceEPSP2{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJMonopole * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
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

struct CalcForceEPEPSimple{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec xi = ep_i[i].pos;
            PS::S64 id_i = ep_i[i].id;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            PS::S32 & n_ngb_i = force[i].n_ngb;
            for(PS::S32 j=0; j<n_jp; j++){
                if(id_i == ep_j[j].id) continue;
                const PS::F64vec rij = xi - ep_j[j].pos;
                const PS::F64 r2_eps = rij * rij + eps2;
                const PS::F64 r_inv = 1.0/sqrt(r2_eps);
                const PS::F64 m_r = ep_j[j].mass * r_inv;
                const PS::F64 m_r3 = m_r * r_inv * r_inv;
                ai -= m_r3 * rij;
                poti -= m_r;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};


struct CalcForceEPSPSimple{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJMonopole * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                const PS::F64vec rij = xi - sp_j[j].pos;
                const PS::F64 r2_eps = rij * rij + eps2;
                const PS::F64 r_inv = 1.0/sqrt(r2_eps);
                const PS::F64 m_r = sp_j[j].mass * r_inv;
                const PS::F64 m_r3 = m_r * r_inv * r_inv;
                ai -= m_r3 * rij;
                poti -= m_r;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
	/*
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
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
	*/
    }
};

#endif //INTRINSIC

/*
struct CalcForceEPEP0{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 r_crit2 = EPJSoft::r_search * EPJSoft::r_search;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec xi = ep_i[i].pos;
            PS::S64 id_i = ep_i[i].id;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            PS::S32 & n_ngb_i = force[i].n_ngb;
            //PS::S32 & id_ngb_i = force[i].id_ngb;
            for(PS::S32 j=0; j<n_jp; j++){
                if(id_i == ep_j[j].id) continue;
                const PS::F64vec rij = xi - ep_j[j].pos;
                const PS::F64 r2 = rij * rij;
                const PS::F64 r2_eps = rij * rij + eps2;
                if(r2 < r_crit2){
                    //if(r2_eps < r_crit2){
                    if(id_i != ep_j[j].id){
                        //id_ngb_i = ep_j[j].id;
                        n_ngb_i++;
                    }
                    continue;
                }
                const PS::F64 r_inv = 1.0/sqrt(r2_eps);
                const PS::F64 m_r = ep_j[j].mass * r_inv;
                const PS::F64 m_r3 = m_r * r_inv * r_inv;
                ai -= m_r3 * rij;
                poti -= m_r;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};


struct CalcForceEPSP0{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJMonopoleScatter * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
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
*/

class PTCLPred{
public:
    PS::F64vec pos;
    PS::F64vec vel;
    void dump(std::ostream & fout=std::cout){
        fout<<"pos="<<pos<<std::endl;
        fout<<"vel="<<vel<<std::endl;
    }
};

class PTCLForce{
public:
    PS::F64vec acc0;
    PS::F64vec acc1;
    PS::F64vec acc0_pla;
    PS::F64vec acc1_pla;
    //PS::F64 pot;
    void clear(){
        acc0 = acc1 = 0.0;
        acc0_pla = acc1_pla = 0.0;
        //pot = 0.0;
    }
    void reduce(){
	acc0 += acc0_pla;
	acc1 += acc1_pla;
    }
    void dump(std::ostream & fout=std::cout){
        fout<<"acc0="<<acc0<<std::endl;
        fout<<"acc1="<<acc1<<std::endl;
    }
};

class PTCLHard{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64 time;
    PS::F64 dt;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc0;
    PS::F64vec acc1;
    PS::F64vec acc0_pla;
    PS::F64vec acc1_pla;
    PS::S32 n_ngb;
    PS::S32 adr_pair;
    PS::F64 pot; // debug
    PS::F64 r_merge;
    static PS::F64 r_factor;
    PTCLHard(){
        id = n_ngb = adr_pair = -1;
    }
    PTCLHard(const FPSoft & fp){
        id = fp.id;
        mass = fp.mass;
        pos = fp.pos;
        vel = fp.vel;
    }
    void setRMerge(){
        //static const PS::F64 rho_ave = 5.51 * ( (1.49597871e13*1.49597871e13*1.49597871e13) / 1.989e33); // [Msun/AU^3]
        //static const PS::F64 rho_ave = 2.0 * ( (1.49597871e13*1.49597871e13*1.49597871e13) / 1.989e33); // [Msun/AU^3]
	static const PS::F64 rho_ave = 3.0 * ( (1.49597871e13*1.49597871e13*1.49597871e13) / 1.989e33); // [Msun/AU^3]
        static const PS::F64 PI = 4.0*atan(1.0);
        static const PS::F64 C = 3.0/(4.0*PI*rho_ave);
        //r_merge = cbrt(C*mass) * 5.0;
        //r_merge = cbrt(C*mass);
	r_merge = cbrt(C*mass) * r_factor;
    }
    void copyFromFP(const FPSoft & fp){
        id = fp.id;
        mass = fp.mass;
        pos = fp.pos;
        vel = fp.vel;
    }
    void dump(std::ostream & fout=std::cout){
        fout<<"id="<<id<<std::endl;
        fout<<"n_ngb="<<n_ngb<<std::endl;
        fout<<"dt="<<dt<<std::endl;
        fout<<"time="<<time<<std::endl;
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"vel="<<vel<<std::endl;
        fout<<"acc0="<<acc0<<std::endl;
        fout<<"acc1="<<acc1<<std::endl;
        fout<<"n_ngb="<<n_ngb<<std::endl;
    }

    static PS::F64 calcDt2nd(const PS::F64vec a0,
			     const PS::F64vec a1,
                             const PS::F64 eta,
			     const PS::F64 a0_offset_sq=0.0){
        const PS::F64 s0 = a0 * a0 + a0_offset_sq;
        const PS::F64 s1 = a1 * a1;
	if(s0 == a0_offset_sq || s1 == 0.0){
            return PS::LARGE_FLOAT;
        }
        else{
	    return eta * sqrt(s0 / s1);
        }
    }

    void setDt2nd(const PTCLForce & force, 
                  const PS::F64 eta, 
                  const PS::F64 dt_limit,
		  const PS::F64 a0_offset_sq=0.0){
        //const PS::F64 dt_ref = calcDt2nd(force.acc0, force.acc1, eta, a0_offset_sq);
	//const PS::F64 dt_ref = calcDt2nd(force.acc0_pla, force.acc1_pla, eta, a0_offset_sq);
	const PS::F64 dt_ref = std::min( calcDt2nd(force.acc0_pla, force.acc1_pla, eta, a0_offset_sq), calcDt2nd(force.acc0, force.acc1, eta) );
        this->dt = dt_limit;
#if 1
        while(this->dt > dt_ref) this->dt *= 0.5;
#elif 1
	this->dt = 1.0/16384;
#else
	this->dt = dt_limit;
#endif
    }

    static PS::F64 calcDt4th(const PS::F64vec a0, const PS::F64vec a1,
                             const PS::F64vec a2, const PS::F64vec a3,
                             const PS::F64 eta,   const PS::F64 a0_offset_sq=0.0){
        const PS::F64 s0 = a0 * a0 + a0_offset_sq;
        const PS::F64 s1 = a1 * a1;
        const PS::F64 s2 = a2 * a2;
        const PS::F64 s3 = a3 * a3;
#if 1
        if(s0 == a0_offset_sq || s1 == 0.0){
            return PS::LARGE_FLOAT;
        }
        else{
            return eta * sqrt( (sqrt(s0*s2) + s1) / (sqrt(s1*s3) + s2) );
        }
#elif 1
	return 1.0/16384;
#else
	return PS::LARGE_FLOAT;
#endif
    }

    void merge(PTCLHard & ptcl_del){
        pos = mass*pos + ptcl_del.mass*ptcl_del.pos;
        vel = mass*vel + ptcl_del.mass*ptcl_del.vel;
        mass = mass + ptcl_del.mass;
        pos /= mass;
        vel /= mass;
        setRMerge();
    }
    
    void correct(const PTCLForce & force,
                 const PS::F64 eta,
                 const PS::F64 dt_limit,
		 const PS::F64 a0_offset_sq=0.0){
        static const PS::F64 inv3 = 1.0 / 3.0;
        const PS::F64 h = 0.5 * dt;
        const PS::F64 hinv = 2.0 / dt;
        const PS::F64vec A0p = (force.acc0 + this->acc0);
        const PS::F64vec A0m = (force.acc0 - this->acc0);
        const PS::F64vec A1p = (force.acc1 + this->acc1)*h;
        const PS::F64vec A1m = (force.acc1 - this->acc1)*h;
        const PS::F64vec vel_new = this->vel + h*( A0p - inv3*A1m );
        this->pos += h*( (this->vel + vel_new) + h*(-inv3*A0m));
        this->vel = vel_new;
        this->acc0 = force.acc0;
        this->acc1 = force.acc1;
#ifdef FORDEBUG
        this->pot = force.pot;
#endif
        this->time += dt;
        const PS::F64vec acc3 = (1.5*hinv*hinv*hinv) * (A1p - A0m);
        const PS::F64vec acc2 = (0.5*hinv*hinv) * A1m + h*acc3;
        //const PS::F64 dt_ref = calcDt4th(this->acc0, this->acc1, acc2, acc3, eta, a0_offset_sq);
        const PS::F64vec A0m_pla = (force.acc0_pla - this->acc0_pla);
        const PS::F64vec A1p_pla = (force.acc1_pla + this->acc1_pla)*h;
        const PS::F64vec A1m_pla = (force.acc1_pla - this->acc1_pla)*h;
        const PS::F64vec acc3_pla = (1.5*hinv*hinv*hinv) * (A1p_pla - A0m_pla);
        const PS::F64vec acc2_pla = (0.5*hinv*hinv) * A1m_pla + h*acc3_pla;
	//const PS::F64 dt_ref = calcDt4th(this->acc0_pla, this->acc1_pla, acc2_pla, acc3_pla, eta, a0_offset_sq);
	const PS::F64 dt_ref = std::min( calcDt4th(this->acc0_pla, this->acc1_pla, acc2_pla, acc3_pla, eta, a0_offset_sq), calcDt4th(this->acc0, this->acc1, acc2, acc3, eta) ); 
        this->acc0_pla = force.acc0_pla;
        this->acc1_pla = force.acc1_pla;
        const PS::F64 dt_old = this->dt;
        assert(dt_old != 0.0);
        this->dt = dt_limit;
        while(this->dt > dt_ref) this->dt *= 0.5;
        this->dt = dt_old*2 < this->dt ?  dt_old*2 : this->dt;
    }
    bool isDead(){
	return (mass > 0.0) ? false : true;
    }
};

PS::F64 PTCLHard::r_factor = 1.0;

class Energy{
public:
    PS::F64 kin;
    PS::F64 pot;
    PS::F64 pot_planet;
    PS::F64 tot;
    PS::F64 disp;
    Energy(){
        kin = pot = tot = disp = pot_planet = 0.0;
    }
    void clear(){
        kin = pot = tot = disp = pot_planet = 0.0;
    }
    void dump(std::ostream & fout=std::cout){
        fout<<"tot="<<tot<<" kin+pot="<<kin+pot<<" kin="<<kin<<" pot="<<pot
	    <<" disp="<<disp<<" pot_planet="<<pot_planet<<std::endl;
    }
};

class MergeLog{
public:
    PS::F64 time;
    PTCLHard ptcl_merged;
    PTCLHard ptcl_dead;
    MergeLog(){
	time = -1.0;
    }
    MergeLog(const PS::F64 t, const PTCLHard & p_m, const PTCLHard & p_d){
	time = t;
	ptcl_merged = p_m;
	ptcl_dead = p_d;
    }
    void dump(std::ostream & fout=std::cout){
	fout<<"time= "<<time<<std::endl;
	fout<<"ptcl_merged.id= "<<ptcl_merged.id<<std::endl;
	fout<<"ptcl_merged.mass= "<<ptcl_merged.mass<<std::endl;
	fout<<"ptcl_merged.pos= "<<ptcl_merged.pos<<std::endl;
	fout<<"ptcl_merged.vel= "<<ptcl_merged.vel<<std::endl;
	fout<<"ptcl_dead.id= "<<ptcl_dead.id<<std::endl;
	fout<<"ptcl_dead.mass= "<<ptcl_dead.mass<<std::endl;
	fout<<"ptcl_dead.pos= "<<ptcl_dead.pos<<std::endl;
	fout<<"ptcl_dead.vel= "<<ptcl_dead.vel<<std::endl;	
    }
};
