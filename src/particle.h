#pragma once

class ForceGrav{
public:
#ifdef CALC_EP_64bit
    PS::F64vec acc;
    PS::F64    phi;
#else
    PS::F32vec acc;
    PS::F32    phi;
#endif
    
    void clear(){
        acc = 0.;
        phi = 0.;
    }
};


class EPGrav{
  public:
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;

    PS::F64 mass;
    
    PS::S32 id;
    PS::S32 id_local;
    PS::S32 myrank;

#ifdef USE_INDIVIDUAL_RADII
    PS::F64 r_out;
    PS::F64 r_out_inv;
    PS::F64 r_search;
#else
    static PS::F64 r_out;
    static PS::F64 r_out_inv;
    static PS::F64 r_search;
#endif

    static PS::F64 eps2;
    static PS::F64 R_cut;
    static PS::F64 R_search0;
    static PS::F64 R_search1;
    static PS::F64 gamma;  
    static PS::F64 g2;       // gamma^2
    static PS::F64 g_1_inv;  // 1/(gamma-1)
    static PS::F64 g2_1_inv; // 1/(gamma^2-1)
    static PS::F64 g_1_inv7; // 1/(gamma-1)^7
    static PS::F64 w_y;      // (g^6-9g^5+45*g^4-60g^3*log(g)-45g^2+9g-1)/3(g-1)^7

    static void setGamma(PS::F64 g){
        gamma  = g;
        g2       = gamma*gamma;
        g_1_inv  = 1./(gamma - 1.);
        g2_1_inv = 1./(g2 - 1.);
        
        PS::F64 g_1_inv3 = g_1_inv * g_1_inv * g_1_inv;
        g_1_inv7 = g_1_inv3 * g_1_inv3 * g_1_inv;

        PS::F64 g4 = g2*g2;
        PS::F64 g6 = g4*g2;
        w_y = 7./3. * (g6 -9.*g4*gamma +45.*g4 -60.*g2*gamma*log(gamma) -45.*g2 +9.*gamma-1.)
            * g_1_inv7;
    }
    PS::F64 getGamma() const{ return gamma; }
    PS::F64 getGamma2() const{ return g2; }
    PS::F64 getGamma_1_inv() const{ return g_1_inv; }
    PS::F64 getGamma2_1_inv() const{ return g2_1_inv; }
    PS::F64 getEps2() const{ return eps2; }
    
    PS::F64 getROut() const { return r_out; }
    PS::F64 getROut_inv() const{ return r_out_inv; }
    PS::F64 getRSearch() const{ return r_search; }
    
    PS::F64vec getPos() const { return pos; }
    PS::F64 getCharge() const { return mass; }

    void copyFromFP(const EPGrav & fp){
        pos = fp.pos;
        vel = fp.vel;
        acc = fp.acc;
        
        mass = fp.mass;
        
        id = fp.id;
        id_local = fp.id_local;
        myrank = fp.myrank;
        
#ifdef USE_INDIVIDUAL_RADII
        r_out     = fp.r_out;
        r_out_inv = fp.r_out_inv;
        r_search  = fp.r_search;
#endif
    }
};


class FPGrav : public EPGrav {
 public:
    PS::F64vec acc_d;
    PS::F64vec acc_gd;
    PS::F64vec jerk;

    PS::F64 phi_s;
    PS::F64 phi_d;
    PS::F64 phi;
    static PS::F64 m_sun;
    static PS::F64 dens;
    
    PS::F64 time;
    PS::F64 dt;
    static PS::F64 dt_tree;
    static PS::F64 dt_min;
    static PS::F64 eta;
    static PS::F64 eta_0;
    static PS::F64 alpha2;
    
    PS::F64 r_planet;
    static PS::F64 rHill_min;
    static PS::F64 rHill_max;

    PS::S32 id_cluster;
    PS::S32 neighbor;
    PS::S32 n_cluster;
#ifdef CHECK_NEIGHBOR
    PS::S32 true_neighbor;
#endif
    
    bool inDomain;
    bool isSent;
    bool isDead;
    bool isMerged;

    PS::F64 getSemimajorAxis() const {
        return 1.0 / (2.0/sqrt(pos*pos) - vel*vel/m_sun);
    }
    PS::F64 getEccentricity(PS::F64 & ax) const {
        PS::F64 r = sqrt(pos*pos);
        PS::F64 rv = pos*vel;
        ax = 1.0 / (2.0/r - vel*vel/m_sun);
        PS::F64 ecccosu = 1. - r/ax;
        PS::F64 eccsinu2 = rv*rv/(m_sun*ax);
        return sqrt(ecccosu*ecccosu + eccsinu2);
    }
    PS::F64 getEccentricity() const {
        PS::F64 ax;
        return getEccentricity(ax);
    }
    PS::F64 getInclination(PS::F64vec & h) const {
        h.x = pos.y*vel.z - pos.z*vel.y;
        h.y = pos.z*vel.x - pos.x*vel.z;
        h.z = pos.x*vel.y - pos.y*vel.x;
        return atan2(sqrt(h.x*h.x + h.y*h.y), h.z);
    }
    PS::F64 getInclination() const {
        PS::F64vec h;
        return getInclination(h);
    }
    PS::F64 getRHill() const {
#ifndef ISOTROPIC
        PS::F64 ax = 1.0 / (2.0/sqrt(pos*pos) - vel*vel/m_sun);
        if ( ax > 0. ){
            return pow(mass/(3.*m_sun), 1./3.) * ax;
        } else {
            return pow(mass/(3.*m_sun), 1./3.) * sqrt(pos*pos);
        }
#else
        return sqrt(mass/m_sun) * sqrt(pos*pos);
#endif
    }
    //PS::F64 getRHill() const {
    //    return pow(mass/(3.*m_sun), 1./3.);
    //}
    PS::F64 getKeplerVelocity() const {
        PS::F64 r = sqrt(pos.x * pos.x + pos.y * pos.y);
        return sqrt(m_sun/r);
    }
#ifdef USE_INDIVIDUAL_RADII
    PS::F64 setROutRSearch(PS::F64 vdisp_k){
        //PS::F64 ax = 1.0 / (2.0/sqrt(pos*pos) - vel*vel/m_sun);
        PS::F64 rHill = std::max(getRHill(), rHill_min);
        PS::F64 v_kep = getKeplerVelocity();
        
        r_out     = R_cut    *rHill;
        r_out_inv = 1. / r_out;
#ifndef ISOTROPIC
        r_search  = R_search0*rHill + R_search1*vdisp_k*v_kep*dt_tree;
#else
        r_search  = R_search0*rHill + R_search1*vdisp_k*dt_tree;
#endif
        
        return rHill;
    }
#else
    static void setROutRSearch(PS::F64 rHill_glb,
#ifndef ISOTROPIC
                               PS::F64 v_kep_glb,
#endif
                               PS::F64 vdisp_k){
        r_out     = R_cut    *rHill_glb;
        r_out_inv = 1. / r_out;
#ifndef ISOTROPIC
        r_search  = R_search0*rHill_glb + R_search1*vdisp_k*v_kep_glb*dt_tree;
#else
        r_search  = R_search0*rHill_glb + R_search1*vdisp_k*dt_tree;
#endif
    }
#endif
    void setRPlanet(PS::F64 m) {
        r_planet = pow(0.75*m/(M_PI*dens), 1./3.);
    }
    void setRPlanet() {
        r_planet = pow(0.75*mass/(M_PI*dens), 1./3.);
    }
    
    void copyFromForce(const ForceGrav & force){
        acc = force.acc;
        phi = force.phi;
    }

    void writeAscii(FILE* fp) const {
        fprintf(fp, "%d\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%d\n",
        //fprintf(fp, "%d\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\n", 
                this->id, this->mass,
                this->pos.x, this->pos.y, this->pos.z,
                this->vel.x, this->vel.y, this->vel.z,
                this->neighbor);
        //this->r_out, this->r_search);
    }
    void readAscii(FILE* fp) {
        fscanf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",
        //fscanf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
               &this->id, &this->mass,
               &this->pos.x, &this->pos.y, &this->pos.z,
               &this->vel.x, &this->vel.y, &this->vel.z,
               &this->neighbor);
               //&this->r_out, &this->r_search);
    }

    void velKick(){
        vel += 0.5*dt_tree*acc;
    }
    
    void calcDeltatInitial(){
        PS::F64 dt_next = 0.5*dt_tree;
        
        PS::F64 a00 = alpha2 * mass / (r_out * r_out);
        PS::F64 dt_2 = (acc_d*acc_d + a00*a00)/(jerk*jerk);
        PS::F64 dt_1 = eta_0 * sqrt(dt_2);
        PS::F64 rem = fmod(time, dt_next);
    
        while( rem != 0.0 ){
            dt_next *= 0.5;
            rem = fmod(time, dt_next);
        }
        if ( dt > 0. ) 
            while( 2.*dt < dt_next ) dt_next *= 0.5;
        while( dt_1 < dt_next ) dt_next *= 0.5;
        if( dt_next < 2.*dt_min ) dt_next = dt_min;
        
        dt = dt_next;
    }
};



class FPHard : public FPGrav {
 public:
    PS::F64vec x0;
    PS::F64vec v0;
    PS::F64vec a0;
    PS::F64vec j0;
    PS::F64vec xp;
    PS::F64vec vp;
    PS::F64vec a2;
    PS::F64vec a3;

    std::vector<PS::S32> n_list;
    std::vector<PS::S32> n_hard_list;
    
    void clearList(){
        std::vector<PS::S32> tmp0, tmp1;
        tmp0.swap(n_list);
        tmp1.swap(n_hard_list);
    }
    void copyList(std::vector<PS::S32> list){
        n_list.resize(list.size());
        std::copy(list.begin(),list.end(),n_list.begin());
    }
    void copyList(PS::S32 * list){
        n_list.clear();
        n_list.reserve(neighbor);
        for ( PS::S32 i=0; i<neighbor; i++ ) n_list.push_back(list[i]);
    }
    void copyHardList(std::vector<PS::S32> list){
        n_hard_list.resize(list.size());
        std::copy(list.begin(),list.end(),n_hard_list.begin());
    }
    void copyHardList(PS::S32 * list){
        n_hard_list.clear();
        n_hard_list.reserve(neighbor);
        for ( PS::S32 i=0; i<neighbor; i++ ) n_hard_list.push_back(list[i]);
    }
    void makeHardList(std::map<PS::S32,PS::S32> & id_map){
        n_hard_list.clear();
        n_hard_list.reserve(neighbor);
        for ( PS::S32 i=0; i<neighbor; i++){
            n_hard_list.push_back(id_map.at(n_list.at(i)));
        }
    }
    FPHard(){
        jerk = 0.;
        x0 = v0 = a0 = j0 = 0.;
        xp = vp = 0.;
        a2 = a3 = 0.;
        
        clearList();
    }
    FPHard(const FPHard & fp) : FPGrav(fp){
        jerk = fp.jerk;
        x0 = fp.x0;
        v0 = fp.v0;
        a0 = fp.a0;
        j0 = fp.j0;
        xp = fp.xp;
        vp = fp.vp;
        a2 = fp.a2;
        a3 = fp.a3;

        copyList(fp.n_list);
        copyHardList(fp.n_hard_list);
    }
    FPHard(const FPGrav & fp) : FPGrav(fp){
        jerk = 0.;
        x0 = v0 = a0 = j0 = 0.;
        xp = vp = 0.;
        a2 = a3 = 0.;
        
        clearList();
    }
    FPHard &operator=(const FPHard & fp){
        FPGrav::operator=(fp);
        if ( this != &fp ){
            jerk = fp.jerk;
            x0 = fp.x0;
            v0 = fp.v0;
            a0 = fp.a0;
            j0 = fp.j0;
            xp = fp.xp;
            vp = fp.vp;
            a2 = fp.a2;
            a3 = fp.a3;
            
            copyList(fp.n_list);
            copyHardList(fp.n_hard_list);
        }
        return *this;
    }
    FPHard &operator=(const FPGrav & fp){
        FPGrav::operator=(fp);
        if ( this != &fp ){
            jerk = 0.;
            x0 = v0 = a0 = j0 = 0.;
            xp = vp = 0.;
            a2 = a3 = 0.;

            clearList();
        }
        return *this;
    }

    void predict(PS::F64 Dt){
        x0 = pos;
        v0 = vel;
        a0 = acc_d;
        j0 = jerk;
        xp = pos + Dt*vel + Dt*Dt*0.5*acc_d + Dt*Dt*Dt*jerk/6.0;
        vp = vel + Dt*acc_d + Dt*Dt*0.5*jerk;
    }
    void correct(PS::F64 Dt){
        PS::F64 Dt2 = Dt*Dt;
        PS::F64 Dt3 = Dt*Dt*Dt;
        PS::F64 Dt2inv = 1.0/Dt2;
        PS::F64 Dt3inv = 1.0/Dt3;
        a2 = -Dt2inv * (6.0*(a0-acc_d) + Dt*(4.0*j0+2.0*jerk));
        a3 = Dt3inv * (12.0*(a0-acc_d) + 6.0*Dt*(j0+jerk));
        pos = xp + a2*Dt2*Dt2/24.0 + 7.0*a3*Dt2*Dt3/720.0;
        vel = vp + a2*Dt3/6.0 + a3*Dt2*Dt2/24.0;
    }

#if 0
    void calcDeltatInitial(){
        PS::F64 dt_next = 0.5*dt_tree;

        PS::F64 a00 = alpha2 * mass / (r_out * r_out);
        PS::F64 dt_2 = (acc_d*acc_d + a00*a00)/(jerk*jerk);
        PS::F64 dt_1 = eta_0 * sqrt(dt_2);
        PS::F64 rem = fmod(time, dt_next);
    
        while( rem != 0.0 ){
            dt_next *= 0.5;
            rem = fmod(time, dt_next);
        }
        while( 2.*dt < dt_next ) dt_next *= 0.5;
        while( dt_1 < dt_next ) dt_next *= 0.5;
        if( dt_next < 2.*dt_min ) dt_next = dt_min;
        
        dt = dt_next;
    }
#endif

    void calcDeltat(){
        //PS::F64 dt_next = 2.*dt;
         PS::F64 dt_next = std::min(2.*dt, 0.5*dt_tree);

        PS::F64 a00 = alpha2 * mass / (r_out * r_out);
        PS::F64vec a1_2 = a2 + a3*dt;
        PS::F64vec a1_3 = a3;
        PS::F64 b1 = sqrt(acc_d*acc_d + a00*a00);
        PS::F64 b1_dot2 = jerk*jerk;
        PS::F64 b1_dot  = sqrt(b1_dot2);
        PS::F64 b1_22 = a1_2*a1_2;
        PS::F64 b1_2  = sqrt(b1_22);
        PS::F64 b1_3 = sqrt(a1_3*a1_3);
        PS::F64 dt_2 = (b1*b1_2 + b1_dot2)/(b1_dot*b1_3 + b1_22);
        PS::F64 dt_1 = eta * sqrt(dt_2);
        PS::F64 rem = fmod(time, dt_next);
        while(rem != 0.0){
            dt_next *= 0.5;
            rem = fmod(time, dt_next);
        }
        while(dt_1 < dt_next && 0.5*dt < dt_next) dt_next *= 0.5;
        if( dt_next < 2.*dt_min ) dt_next = dt_min;
        
        dt = dt_next;
    }
    
};


template <class Tpsys>
PS::F64 calcVelDisp(Tpsys & pp,
                    const PS::S32 n_tot,
                    const PS::S32 n_loc)
{
    //PS::F64 m_sun   = FPGrav::m_sun;
    PS::F64 ecc_rms_loc = 0.0;
    PS::F64 inc_rms_loc = 0.0;
    //PS::F64 v_kep_max_loc = 0.0;

#pragma omp parallel for reduction (+:ecc_rms_loc, inc_rms_loc)
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64 ax;
        PS::F64 ecc = pp[i].getEccentricity(ax);
        PS::F64vec h;
        PS::F64 inc = pp[i].getInclination(h);
        //PS::F64 v_kep = sqrt(m_sun/ax);

        ecc_rms_loc += ecc*ecc;
        inc_rms_loc += inc*inc;
        //#pragma omp critical
        //{
            //if ( v_kep > v_kep_max_loc ) v_kep_max_loc = v_kep;
        //}
    }
    
    PS::F64 ecc_rms2 = PS::Comm::getSum(ecc_rms_loc)/n_tot;
    PS::F64 inc_rms2 = PS::Comm::getSum(inc_rms_loc)/n_tot;
    //PS::F64 v_kep_max = PS::Comm::getMaxValue(v_kep_max_loc);
    
    //return (ecc_rms2 + inc_rms2)*v_kep_max;
    return sqrt(ecc_rms2 + inc_rms2);
}

#ifdef ISOTROPIC

template <class Tpsys>
PS::F64 calcRandomVel(Tpsys & pp,
                      const PS::S32 n_tot,
                      const PS::S32 n_loc)
{
    PS::F64    vel2_sum_loc = 0.;
    PS::F64vec vel_sum_loc  = 0.;

#pragma omp parallel for reduction (+:vel2_sum_loc, vel_sum_loc)
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64vec vel  = pp[i].vel;
        PS::F64    vel2 = vel*vel;

        vel2_sum_loc += vel2;
        vel_sum_loc  += vel;
    }

    PS::F64    vel2_ave = PS::Comm::getSum(vel2_sum_loc)/n_tot;
    PS::F64vec vel_ave  = PS::Comm::getSum(vel_sum_loc)/n_tot;
    
    return sqrt(vel2_ave - vel_ave*vel_ave);
}

#endif

template <class Tpsys>
void setCutoffRadii(Tpsys & pp)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    const PS::S32 n_tot = pp.getNumberOfParticleGlobal();
#ifndef ISOTROPIC
    PS::F64 v_disp_k = calcVelDisp(pp, n_tot, n_loc);
#else
    PS::F64 v_disp   = calcRandomVel(pp, n_tot, n_loc);
#endif
    
#ifdef USE_INDIVIDUAL_RADII
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
#ifndef ISOTROPIC
        pp[i].setROutRSearch(v_disp_k);
#else
        pp[i].setROutRSearch(v_disp);
#endif
        pp[i].setRPlanet();
    }
#else //USE_INDIVIDUAL_RADII
    PS::F64 rHill_loc = 0.;
    PS::F64 v_kep_loc = 0.;
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64 rHill = 0;
        if ( FPGrav::rHill_max <= 0 ) {
            rHill = std::max(pp[i].getRHill(), FPGrav::rHill_min);
        } else {
            rHill = std::min(FPGrav::rHill_max, std::max(pp[i].getRHill(), FPGrav::rHill_min));
        }
#ifndef ISOTROPIC
        PS::F64 v_kep = pp[i].getKeplerVelocity();
#endif
#pragma omp critical
        {
            if ( rHill > rHill_loc ) rHill_loc = rHill;
#ifndef ISOTROPIC
            if ( v_kep > v_kep_loc ) v_kep_loc = v_kep;
#endif
        }
        pp[i].setRPlanet();
    }
    PS::F64 rHill_glb = PS::Comm::getMaxValue(rHill_loc);
#ifndef ISOTROPIC
    PS::F64 v_kep_glb = PS::Comm::getMaxValue(v_kep_loc);
    FPGrav::setROutRSearch(rHill_glb, v_kep_glb, v_disp_k);
#else
    FPGrav::setROutRSearch(rHill_glb, v_disp);
#endif
#endif //USE_INDIVIDUAL_RADII
}
