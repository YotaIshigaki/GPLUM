
/*
class FileHeader{
    FileHeader(const PS::S32 n_loc, const PS::S32 n_glb){
    }
};
*/

class Force_t{
public:
    PS::F64vec acc; // total
    PS::F64    pot;
    PS::F64vec acc_dash; // dashpot
    void clear(){
        acc = 0.0;
        pot = 0.0;
        acc_dash = 0.0;
    }
};

class FP_t{
public:
    PS::F64vec pos_car; // cartesian
    PS::F64vec pos_cyl; // cyl
    PS::F64    mass;
    PS::F64vec vel; // cartesian
    PS::F64vec vel_full; // for calculating disipation energy
    PS::S64    id;
    PS::F64vec acc;
    PS::F64 pot;
    PS::F64 r_coll;
    PS::F64vec acc_dash;
    static inline PS::F64 r_search;
    static inline PS::F64 eps;
    static inline PS::F64 kappa;
    static inline PS::F64 eta;
    PS::F64vec getPos() const {
        return pos_cyl;
    }
    PS::F64vec getPosCar() const {
        return pos_car;
    }
    void setPos(const PS::F64vec & pos_new){
        pos_cyl = pos_new;
    }

    PS::F64 getCharge() const {
        return mass;
    }
    PS::F64 getRSearch() const {
        return r_search;
    }
    void copyFromForce(const Force_t & force) {
        acc = force.acc;
        pot = force.pot;
        acc_dash = force.acc_dash;
    }
    void writeAscii(FILE* fp) const {
        fprintf(fp, "%lld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
                this->id, this->mass,
                this->pos_car.x, this->pos_car.y, this->pos_car.z,
                this->vel.x, this->vel.y, this->vel.z);
    }
    void readAscii(FILE* fp) {
        fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
	       &this->id, &this->mass,
               &this->pos_car.x, &this->pos_car.y, &this->pos_car.z,
               &this->vel.x, &this->vel.y, &this->vel.z);
    }
};

class EPI_t{
public:
    PS::F64vec pos_car; // cartesian
    PS::F64vec pos_cyl; // cyl
    PS::F64    mass;
    //PS::F64vec vel; // cartesian
    PS::F64vec vel_full; // for calculating disipation energy
    PS::S64    id;
    PS::F64 r_coll;
    PS::F64vec getPos() const {
        return pos_cyl;
    }
    PS::F64vec getPosCar() const {
        return pos_car;
    }
    void setPos(const PS::F64vec & pos_new){
        pos_cyl = pos_new;
    }
    PS::F64 getCharge() const {
        return mass;
    }
    PS::F64 getRSearch() const {
        return FP_t::r_search;
    }
    void copyFromFP(const FP_t & fp){ 
        pos_car  = fp.pos_car;
	pos_cyl  = fp.pos_cyl;
        mass = fp.mass;
        //vel  = fp.vel;
        vel_full  = fp.vel_full;
        id   = fp.id;
        r_coll = fp.r_coll;
    }
};

using EPJ_t = EPI_t;

class MyMomentMonopole{
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec pos_car;
    PS::F64ort boundary;
    MyMomentMonopole() : mass(0.0), pos(PS::F64vec(0.0)), pos_car(PS::F64vec(0.0)), boundary(PS::F64ort(PS::F64vec(-100.0), PS::F64vec(100.0))){}
    MyMomentMonopole(const PS::F64 m, const PS::F64vec & p, const PS::F64vec & p_car, const PS::F64ort & b) : mass(m), pos(p), pos_car(p_car), boundary(b){}
    void init(){
        mass = 0.0;
        pos = 0.0;
        pos_car = 0.0;
        boundary.init();
    }
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64 getCharge() const {
        return mass;
    }
    template<class Tepj>
    void accumulateAtLeaf(const Tepj & epj){
        mass += epj.getCharge();
        pos += epj.getCharge() * epj.getPos();
        pos_car += epj.getCharge() * epj.getPosCar();
        boundary.merge(epj.getPos());
    }
    template<class Tepj>
    void accumulateAtLeaf2(const Tepj & epj){}
    void set(){
        pos = pos / mass;
        pos_car = pos_car / mass;
    }
    void accumulate(const MyMomentMonopole & mom){
        mass += mom.mass;
        pos += mom.mass * mom.pos;
        pos_car += mom.mass * mom.pos_car;
        boundary.merge(mom.boundary);
    }
    void accumulate2(const MyMomentMonopole & mom){}
    // for DEBUG 
    void dump(std::ostream & fout = std::cout) const {
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"pos_car="<<pos_car<<std::endl;
    }
};

class MySPJMonopole{
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec pos_car;
    PS::F64ort boundary;
    template<class Tmom>
    void copyFromMoment(const Tmom & mom){
        this->mass     = mom.mass;
        this->pos      = mom.pos;
        this->pos_car  = mom.pos_car;
        this->boundary = mom.boundary;
    }
    void clear(){
        mass = 0.0;
        pos = 0.0;
        pos_car = 0.0;
        boundary.init();
    }
    PS::F64 getCharge() const {
        return mass;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64vec getPosCar() const {
        return pos_car;
    }
    MyMomentMonopole convertToMoment() const {
        return MyMomentMonopole(mass, pos, pos_car, boundary);
    }
};

class MyMomentQuadrupole{
public:
    PS::F64vec pos;
    PS::F64 mass;
    PS::F64mat quad;
    PS::F64vec pos_car;
    void init(){
        pos = 0.0;
        mass = 0.0;
        quad = 0.0;
        pos_car = 0.0;
    }
    MyMomentQuadrupole(){
        mass = 0.0;
        pos = 0.0;
        quad = 0.0;
        pos_car = 0.0;
    }
    MyMomentQuadrupole(const PS::F64 m, const PS::F64vec & p, const PS::F64mat & q, const PS::F64vec & p_car){
        mass = m;
        pos = p;
        quad = q;
        pos_car = p_car;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    template<class Tepj>
    void accumulateAtLeaf(const Tepj & epj){
        mass += epj.getCharge();
        pos  += epj.getCharge() * epj.getPos();
        pos_car += epj.getCharge() * epj.getPosCar();
    }
    template<class Tepj>
    void accumulateAtLeaf2(const Tepj & epj){
        PS::F64 ctmp = epj.getCharge();
        PS::F64vec ptmp = epj.getPosCar() - this->pos_car;
        PS::F64 cx = ctmp * ptmp.x;
        PS::F64 cy = ctmp * ptmp.y;
        PS::F64 cz = ctmp * ptmp.z;
        this->quad.xx += cx * ptmp.x;
        this->quad.yy += cy * ptmp.y;
        this->quad.zz += cz * ptmp.z;
        this->quad.xy += cx * ptmp.y;
        this->quad.xz += cx * ptmp.z;
        this->quad.yz += cy * ptmp.z;
    }
    void set(){
        pos = pos / mass;
        pos_car = pos_car / mass;
    }
    void accumulate(const MyMomentQuadrupole & mom){
        mass += mom.mass;
        pos += mom.mass * mom.pos;
        pos_car += mom.mass * mom.pos_car;
    }
    void accumulate2(const MyMomentQuadrupole & mom){
        PS::F64 mtmp = mom.mass;
        PS::F64vec ptmp = mom.pos_car - this->pos_car;
        PS::F64 cx = mtmp * ptmp.x;
        PS::F64 cy = mtmp * ptmp.y;
        PS::F64 cz = mtmp * ptmp.z;
        this->quad.xx += cx * ptmp.x + mom.quad.xx;
        this->quad.yy += cy * ptmp.y + mom.quad.yy;
        this->quad.zz += cz * ptmp.z + mom.quad.zz;
        this->quad.xy += cx * ptmp.y + mom.quad.xy;
        this->quad.xz += cx * ptmp.z + mom.quad.xz;
        this->quad.yz += cy * ptmp.z + mom.quad.yz;
    }
    void dump(std::ostream & fout = std::cout) const {
        fout<<"mass= "<<mass<<std::endl;
        fout<<"pos= "<<pos<<std::endl;
        fout<<"quad= "<<quad<<std::endl;
    }
};

class MySPJQuadrupole{
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec pos_car;
    PS::F64mat quad;
    PS::F64 getCharge() const {
        return mass;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64vec getPosCar() const {
        return pos_car;
    }
    void copyFromMoment(const MyMomentQuadrupole & mom){
        mass = mom.mass;
        pos = mom.pos;
        quad = mom.quad;
        pos_car = mom.pos_car;
    }
    MyMomentQuadrupole convertToMoment() const {
        return MyMomentQuadrupole(mass, pos, quad, pos_car);
    }
    void clear(){
        mass = 0.0;
        pos = 0.0;
        quad = 0.0;
        pos_car = 0.0;
    }
};


#if defined(USE_PIKG_KERNEL)

struct Epi0{
    PS::F32vec pos;
    PS::F32vec vel;
    PS::F32    mass;
    PS::F32    r_coll;
    PS::S64    id;
};
struct Epj0{
    PS::F32vec pos;
    PS::F32vec vel;
    PS::F32    mass;
    PS::F32    r_coll;
    PS::S64    id;
};
struct Force0{
    PS::F32vec acc;
    PS::F32vec acc_dash;
    PS::F32    pot;
};
struct Epi1{
    PS::F32vec pos;
};
struct Spj1{
    PS::F32vec pos;
    PS::F32    mass;
};
struct Force1{
    PS::F32vec acc;
    PS::F32    pot;
};

#include"kernel_ep.hpp"

template<typename Tpi, typename Tpj, typename Tfi>
void CalcForceEp(const Tpi * ep_i,
		 const PS::S32 n_ip,
		 const Tpj * ep_j,
		 const PS::S32 n_jp,
		 Tfi * force) {
    const PS::F32 eps2  = (PS::F32)(FP_t::eps*FP_t::eps);
    const PS::F32 kappa = (PS::F32)FP_t::kappa;
    const PS::F32 eta   = (PS::F32)FP_t::eta;  
    Epi0 epi[n_ip];
    Force0 f[n_ip];
    for(int i=0;i<n_ip;i++){
        epi[i].pos = ep_i[i].pos - ep_i[0].pos;
	epi[i].mass   = ep_i[i].mass;
	epi[i].r_coll = ep_i[i].r_coll;
	epi[i].vel    = ep_i[i].vel_full;
    }
    Epj0 epj[n_jp];
    for(int i=0;i<n_jp;i++){
	epj[i].pos = ep_j[i].pos - ep_i[0].pos;
	epj[i].mass   = ep_j[i].mass;
	epj[i].r_coll = ep_j[i].r_coll;
	epj[i].vel    = ep_j[i].vel_full;
    }
    CalcForceEpEpImpl(eps2, kappa, eta)(epi,n_ip,epj,n_jp,f);
    for(int i=0;i<n_ip;i++){
        force[i].acc      = f[i].acc;
        force[i].acc_dash = f[i].acc_dash;
	force[i].pot      = f[i].pot;
    }
}



#include"kernel_sp.hpp"

template<typename Tpi, typename Tpj, typename Tfi>
void CalcForceSpMono(const Tpi * ep_i,
		     const PS::S32 n_ip,
		     const Tpj * ep_j,
		     const PS::S32 n_jp,
		     Tfi * force) {
    const auto eps2  = FP_t::eps*FP_t::eps;
    Epi1 epi[n_ip];
    Force1 f[n_ip];
    for(int i=0;i<n_ip;i++){
        epi[i].pos = ep_i[i].pos - ep_i[0].pos;
    }
    Epj1 epj[n_jp];
    for(int i=0;i<n_jp;i++){
	epj[i].pos = ep_j[i].pos - ep_i[0].pos;
	epj[i].mass = ep_j[i].mass;
    }
    CalcForceEpSpImpl(eps2)(epi,n_ip,epj,n_jp,f);
    for(int i=0;i<n_ip;i++){
        force[i].acc = f[i].acc;
	force[i].pot = f[i].pot;
    }
}

#else

template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceEp{
    void operator ()(const Tpi * pi,
                     const PS::S32 ni,
                     const Tpj * pj,
                     const PS::S32 nj,
                     Tforce * force){
        const auto eps2  = FP_t::eps*FP_t::eps;
        const auto kappa = FP_t::kappa;
        const auto eta   = FP_t::eta;
	PS::F64vec xj[nj];
	for(auto j=0; j<nj; j++){
	  xj[j] = pj[j].getPosCar();;
	}
        for(auto i=0; i<ni; i++){
            const PS::F64vec xi = pi[i].getPosCar();
            PS::F64vec ai = 0.0;
            PS::F64vec ai_dash = 0.0;
            PS::F64 poti = 0.0;
            for(auto j=0; j<nj; j++){
                const auto r_coll    = (pi[i].r_coll + pj[j].r_coll);
                const auto r_coll_sq = r_coll*r_coll;
                //PS::F64vec rij       = xi - pj[j].getPosCar();
		PS::F64vec rij       = xi - xj[j];
                if(pi[i].id == pj[j].id) continue;
                PS::F64 r2_real = rij * rij + eps2;
                PS::F64 r2      = std::max(r2_real, r_coll_sq);
                PS::F64 r_inv   = 1.0/sqrt(r2);
                PS::F64 r2_inv  = r_inv * r_inv;
                PS::F64 pot = r_inv * pj[j].getCharge();
                if(r_coll_sq > r2_real){
                    ai     -= pj[j].getCharge() / (r_coll_sq*r_coll) * rij;
                    PS::F64 pot_offset = -1.5/r_coll;
                    poti   += 0.5*pj[j].getCharge()*(0.5*r2_real/(r_coll_sq*r_coll) + pot_offset);
                }
                else{
                    ai     -= pot * r2_inv * rij;
                    poti   -= 0.5 * pot;
                }
                //poti   -= 0.5 * pj[j].getCharge() * sqrt(r2_real) * r2_inv;
                if(r_coll_sq > r2_real){
                    PS::F64 m_r = pj[j].mass / (pi[i].mass+pj[j].mass);
                    PS::F64 r   = sqrt(r2_real);
                    PS::F64 dr  = r_coll-r ;
                    ai += kappa * m_r * dr/r * rij;
                    poti += 0.5*kappa*m_r*dr*dr * 0.5;
                    PS::F64vec vij = pi[i].vel_full - pj[j].vel_full;
                    PS::F64 rv = rij*vij;
                    //PS::F64vec a_eta = eta * m_r * rv * r2_inv * rij;
                    PS::F64vec a_eta = eta * m_r * rv / r2_real * rij;
                    ai_dash += a_eta;
                    ai += a_eta;
                }
            }
            force[i].acc += ai;
            force[i].acc_dash += ai_dash;
            force[i].pot += poti;
        }
    }
};

template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceSpMono{
    void operator ()(const Tpi * pi,
                     const PS::S32 ni,
                     const Tpj * pj,
                     const PS::S32 nj,
                     Tforce * force){
        const auto eps2 = FP_t::eps*FP_t::eps;
	PS::F64vec xj[nj];
	for(auto j=0; j<nj; j++){
	  xj[j] = pj[j].getPosCar();
	}
        for(auto i=0; i<ni; i++){
            const auto xi = pi[i].getPosCar();
            PS::F64vec ai = 0.0;
            PS::F64 poti  = 0.0;
            for(auto j=0; j<nj; j++){
	      //PS::F64vec rij    = xi - pj[j].getPosCar();
	        const auto rij    = xi - xj[j];
                auto r3_inv = rij * rij + eps2;
                auto r_inv  = 1.0/sqrt(r3_inv);
                r3_inv  = r_inv * r_inv;
                r_inv  *= pj[j].getCharge();
                r3_inv *= r_inv;
                ai     -= r3_inv * rij;
                poti   -= 0.5*r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};

template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceSpQuad{
    void operator ()(const Tpi * pi,
                     const PS::S32 ni,
                     const Tpj * pj,
                     const PS::S32 nj,
                     Tforce * force){
        const auto eps2 = FP_t::eps*FP_t::eps;
	PS::F64vec xj[nj];
	for(auto j=0; j<nj; j++){
	  xj[j] = pj[j].getPosCar();
	}
        for(auto i=0; i<ni; i++){
            PS::F64vec xi = pi[i].getPosCar();
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(auto j=0; j<nj; j++){
                PS::F64 mj = pj[j].getCharge();
                //PS::F64vec xj= pj[j].getPosCar();
                //PS::F64vec rij= xi - xj;
		PS::F64vec rij = xi - xj[j];
                PS::F64 r2 = rij * rij + eps2;
                PS::F64mat qj = pj[j].quad;
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
            force[i].acc += ai;
            force[i].pot += 0.5*poti;
        }
    }
};

#endif
