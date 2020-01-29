#pragma once

#define SECONDORDER  5.e-1
#define THIRDORDER   1.6666666666666667e-1
#define FOURTHORDER  4.1666666666666667e-2
#define FIFTHORDER   8.3333333333333333e-3
#define SIXTHORDER   1.3888888888888889e-3
#define SEVENTHORDER 1.9841269841269841e-4
#define ALPHA        1.1666666666666667e-3
#define ONE_TWELFTH  8.3333333333333333e-2
#define ONE_SIXTIETH 1.6666666666666667e-2
#define ONE_420TH    2.3809523809523810e-2

    
inline PS::F64 calcDt2nd(PS::F64 eta,
                         PS::F64 alpha2,
                         PS::F64 acc0,
                         PS::F64vec acc,
                         PS::F64vec jerk){
    PS::F64 Acc2  = acc*acc + alpha2*acc0*acc0;
    PS::F64 Jerk2 = jerk*jerk;
    //PS::F64 dt2 = (Jerk2>0.) ? Acc2/Jerk2 : std::numeric_limits<double>::max();
    //return eta * sqrt(dt2);
    return (Jerk2>0.) ? eta*sqrt(Acc2/Jerk2) : std::numeric_limits<double>::max();
}

inline PS::F64 calcDt3rd(PS::F64 eta,
                         PS::F64 alpha2,
                         PS::F64 acc0,
                         PS::F64vec acc,
                         PS::F64vec snap){
    PS::F64 Acc2  = acc*acc + alpha2*acc0*acc0;
    PS::F64 Snap2 = snap*snap;
    return (Snap2>0.) ? eta*pow(Acc2/Snap2, 0.25) : std::numeric_limits<double>::max();
}

inline PS::F64 calcDt4th(PS::F64 eta,
                         PS::F64 alpha2,
                         PS::F64 acc0,
                         PS::F64vec acc,
                         PS::F64vec jerk,
                         PS::F64vec snap,
                         PS::F64vec crac){
    PS::F64 Acc   = sqrt(acc*acc + alpha2*acc0*acc0);
    PS::F64 Jerk2 = jerk*jerk;
    PS::F64 Jerk  = sqrt(Jerk2);
    PS::F64 Snap2 = snap*snap;
    PS::F64 Snap  = sqrt(Snap2);
    PS::F64 Crac  = sqrt(crac*crac);
    //PS::F64 dt2 = (Jerk>0.) ? (Acc*Snap + Jerk2)/(Jerk*Crac + Snap2) : std::numeric_limits<double>::max();
    //return eta * sqrt(dt2);
    return (Jerk>0.) ? eta*sqrt((Acc*Snap + Jerk2)/(Jerk*Crac + Snap2)) : std::numeric_limits<double>::max();
}

inline PS::F64 calcDt6th(PS::F64 eta,
                         PS::F64 alpha2,
                         PS::F64 acc0,
                         PS::F64vec acc,
                         PS::F64vec jerk,
                         PS::F64vec snap,
                         PS::F64vec crac,
                         PS::F64vec pop,
                         PS::F64vec a5){
    PS::F64 Acc   = sqrt(acc*acc + alpha2*acc0*acc0);
    PS::F64 Jerk2 = jerk*jerk;
    PS::F64 Snap  = sqrt(snap*snap);
    PS::F64 Crac  = sqrt(crac*crac);
    PS::F64 Pop2  = pop*pop;
    PS::F64 A5    = sqrt(a5*a5);
    return (Jerk2>0.) ? eta*pow((Acc*Snap + Jerk2)/(Crac*A5 + Pop2),THIRDORDER) : std::numeric_limits<double>::max();
}

class FPGrav : public EPGrav {
public:
#ifdef INTEGRATE_6TH_SUN
    PS::F64vec acc_;
    PS::F64vec snap_s;
#endif
    
    PS::F64vec acc_gd;
    PS::F64vec jerk_s;
    PS::F64vec jerk_d;
    
#ifdef INDIRECT_TERM
    static PS::F64vec acc_indirect;
    static PS::F64vec pos_g;
    static PS::F64vec vel_g;
    static PS::F64    mass_tot;
#endif

    PS::F64 phi_s;
    PS::F64 phi_d;
    PS::F64 phi;
    static PS::F64 m_sun;
    static PS::F64 dens;
    
    PS::F64 time;
    PS::F64 dt;
    PS::F64 acc0;
    static PS::F64 dt_tree;
    static PS::F64 dt_min;
    static PS::F64 eta;
    static PS::F64 eta_0;
    static PS::F64 eta_sun;
    static PS::F64 eta_sun0;
    static PS::F64 alpha2;
    
    PS::F64 r_planet;
    static PS::F64 r_cut_min;
    static PS::F64 r_cut_max;
    static PS::F64 p_cut;

    PS::S32 id_cluster;
    PS::S32 n_cluster;
#ifdef CHECK_NEIGHBOR
    PS::S32 true_neighbor;
#endif
    
    PS::S32 id_neighbor;
    
    bool inDomain;
    bool isSent;
    bool isDead;
    bool isMerged;

    static PS::F64 getSolarMass() { return m_sun; }
    
    PS::F64 getSemimajorAxis() const {
        return 1.0 / (2.0/sqrt(pos*pos) - vel*vel/m_sun);
    }
    PS::F64 getSemimajorAxis2() const {
        PS::F64 ax;
        if ( getEccentricity(ax) < 0.5 ) {
            return ax;
        } else {
            return sqrt(pos*pos);
        }
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
        PS::F64 ax = getSemimajorAxis2();
        return pow(mass/(3.*m_sun), 1./3.) * ax;
    }
    PS::F64 getKeplerVelocity() const {
        PS::F64 r = sqrt(pos.x * pos.x + pos.y * pos.y);
        return sqrt(m_sun/r);
    }

#ifdef INTEGRATE_6TH_SUN
    void setAcc_() { acc_ = acc_s + acc_d; }
#endif
    
#ifdef INDIRECT_TERM
    static PS::F64 getKineticEnergyOfSystem(){
        return 0.5 * mass_tot * vel_g * vel_g;
    }
#endif
    
#ifdef USE_INDIVIDUAL_CUTOFF
    
    PS::F64 setROutRSearch(PS::F64 vdisp_k){
        PS::F64 rHill = getRHill();
        PS::F64 v_kep = getKeplerVelocity();
        PS::F64 ax    = getSemimajorAxis2();
        
        if ( FPGrav::r_cut_max <= 0. ) {
            r_out = std::max(R_cut * pow(ax,-p_cut) * rHill, FPGrav::r_cut_min);
        } else {
            r_out = std::min(FPGrav::r_cut_max, std::max(R_cut * pow(ax,-p_cut) * rHill, FPGrav::r_cut_min) );
        }
        
#ifdef TEST_PTCL
        if ( r_out != 0. ) {
            r_out_inv = 1. / r_out;
        } else {
            r_out_inv = (PS::F64)std::numeric_limits<PS::F32>::max();
        }
#else
        r_out_inv = 1. / r_out;
#endif
            
#ifndef ISOTROPIC
        r_search  = R_search0*r_out + R_search1*vdisp_k*v_kep*dt_tree;
#else
        r_search  = R_search0*r_out + R_search1*vdisp_k*dt_tree;
#endif
#ifdef TEST_PTCL
        if ( r_out == 0. ) r_search = 0.;

        assert ( ( r_out > 0. && r_search > 0. ) ||
                 ( r_out == 0. && r_search == 0. && mass == 0. ) );
#else    
        assert ( r_out > 0. && r_search > 0. );
#endif
        
        return rHill;
    }
    
#else //USE_INDIVIDUAL_CUTOFF
    
    static void setROutRSearch(PS::F64 rHill_a_glb,
#ifndef ISOTROPIC
                               PS::F64 v_kep_glb,
#endif
                               PS::F64 vdisp_k){
        
        if ( FPGrav::r_cut_max <= 0 ) {
            r_out = std::max(R_cut * rHill_a_glb, FPGrav::r_cut_min);
        } else {
            r_out = std::min(FPGrav::r_cut_max, std::max(R_cut * rHill_a_glb, FPGrav::r_cut_min) );
        }

#ifdef TEST_PTCL
        if ( r_out != 0. ) {
            r_out_inv = 1. / r_out;
        } else {
            r_out_inv = (PS::F64)std::numeric_limits<PS::F32>::max();
        }
#else
        r_out_inv = 1. / r_out;
#endif
        
#ifndef ISOTROPIC
        r_search  = R_search0*r_out + R_search1*vdisp_k*v_kep_glb*dt_tree;
#else
        r_search  = R_search0*r_out + R_search1*vdisp_k*dt_tree;
#endif
#ifdef TEST_PTCL
        if ( r_out == 0. ) r_search = 0.;
        
        assert ( ( r_out > 0. && r_search > 0. ) ||
                 ( r_out == 0. && r_search == 0. ) );
#else    
        assert ( r_out > 0. && r_search > 0. );
#endif
    }
    
#endif //USE_INDIVIDUAL_CUTOFF
    
    void setRPlanet(PS::F64 m) {
        r_planet = pow(0.75*m/(M_PI*dens), 1./3.);
    }
    void setRPlanet() {
        r_planet = pow(0.75*mass/(M_PI*dens), 1./3.);
    }
    
    void copyFromForce(const ForceGrav & force){
        acc = force.acc;
        phi = force.phi;
        neighbor = force.neighbor;

        id_neighbor = force.id_neighbor;
    }

    void writeAscii(FILE* fp) const {
        if ( !fprintf(fp, "%d\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%d\n",
                      //!fprintf(fp, "%d\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\n", 
                      this->id, this->mass,
                      this->pos.x, this->pos.y, this->pos.z,
                      this->vel.x, this->vel.y, this->vel.z,
                      this->neighbor) ) {
            //this->r_out, this->r_search) ){
            errorMessage("The particle data has NOT been correctly written.");
            PS::Abort();
        }   
    }
    void readAscii(FILE* fp) {
        if ( !fscanf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",
                     //fscanf(fp, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                     &this->id, &this->mass,
                     &this->pos.x, &this->pos.y, &this->pos.z,
                     &this->vel.x, &this->vel.y, &this->vel.z,
                     &this->neighbor) ) {
            //&this->r_out, &this->r_search) ) {
            errorMessage("The particle data has NOT been correctly read.");
            PS::Abort();
        }
    }

    void velKick(){
#ifdef INDIRECT_TERM
        vel += 0.5*dt_tree*(acc + acc_indirect);
#else
        vel += 0.5*dt_tree*acc;
#endif
    }
    
    void calcDeltatInitial(){
        PS::F64 dt_next = 0.5*dt_tree;
        
#ifndef INTEGRATE_6TH_SUN
        PS::F64 dt_1 = std::min(calcDt2nd(eta_0,    alpha2, acc0, acc_d, jerk_d),
                                calcDt2nd(eta_sun0, alpha2,   0., acc_s, jerk_s));
        //PS::F64 dt_1 = calcDt2nd(eta_0, alpha2, acc0, acc_d, jerk_d);
#else
        PS::F64 dt_1 = std::min(calcDt2nd(eta_0,    alpha2, acc0, acc_d, jerk_d),
                                calcDt3rd(eta_sun0, alpha2,   0., acc_s, snap_s));
#endif
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
    PS::F64vec a0_s;
    PS::F64vec j0_s;
    PS::F64vec a0_d;
    PS::F64vec j0_d;
    PS::F64vec xp;
    PS::F64vec vp;
#ifdef INTEGRATE_6TH_SUN
    PS::F64vec s0_s;
    PS::F64vec ap;
#endif
    
#ifndef INTEGRATE_6TH_SUN
    PS::F64vec a2_s;
    PS::F64vec a3_s;
#else
    PS::F64vec a3_s;
    PS::F64vec a4_s;
    PS::F64vec a5_s;
#endif
    PS::F64vec a2_d;
    PS::F64vec a3_d;
    
    PS::F64 time_c;
    
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
    template <class Tpsys>
    void makeHardList(std::map<PS::S32,PS::S32> & id_map,
                      Tpsys & pp){
        n_hard_list.clear();
        for ( PS::S32 i=0; i<neighbor; i++){
            PS::S32 id_hard = id_map.at(n_list.at(i));
            if ( pp[id_hard].mass > 0. ) n_hard_list.push_back(id_hard);
        }
        neighbor = n_hard_list.size();
    }
    FPHard(){
        x0   = v0   = 0.;
        a0_s = j0_s = 0.;
        a0_d = j0_d = 0.;
        xp   = vp   = 0.;
#ifdef INTEGRATE_6TH_SUN
        s0_s   = 0.;
        ap     = 0.;
#endif
#ifndef INTEGRATE_6TH_SUN
        a2_s = a3_s = 0.;
#else
        a3_s = a4_s = a5_s = 0.;
#endif
        a2_d = a3_d = 0.;
        time_c = 0.;
        
        clearList();
    }
    FPHard(const FPHard & fp) : FPGrav(fp){
        x0   = fp.x0;
        v0   = fp.v0;
        a0_s = fp.a0_s;
        j0_s = fp.j0_s;
        a0_d = fp.a0_d;
        j0_d = fp.j0_d;
        xp   = fp.xp;
        vp   = fp.vp;
#ifdef INTEGRATE_6TH_SUN
        s0_s   = fp.s0_s;
        ap     = fp.ap;
#endif
#ifndef INTEGRATE_6TH_SUN
        a2_s = fp.a2_s;
        a3_s = fp.a3_s;
#else
        a3_s = fp.a3_s;
        a4_s = fp.a4_s;
        a5_s = fp.a5_s;
#endif
        a2_d = fp.a2_d;
        a3_d = fp.a3_d;
        
        time_c = fp.time_c;
        
        copyList(fp.n_list);
        copyHardList(fp.n_hard_list);
    }
    FPHard(const FPGrav & fp) : FPGrav(fp){
        x0   = v0   = 0.;
        a0_s = j0_s = 0.;
        a0_d = j0_d = 0.;
        xp   = vp   = 0.;
#ifdef INTEGRATE_6TH_SUN
        s0_s   = 0;
        ap     = 0;
#endif
#ifndef INTEGRATE_6TH_SUN
        a2_s = a3_s = 0.;
#else
        a3_s = a4_s = a5_s = 0.;
#endif
        a2_d = a3_d = 0.;
        
        time_c = fp.time;
        time = 0;
        
        clearList();
    }
    FPHard &operator=(const FPHard & fp){
        FPGrav::operator=(fp);
        if ( this != &fp ){
            x0   = fp.x0;
            v0   = fp.v0;
            a0_s = fp.a0_s;
            j0_s = fp.j0_s;
            a0_d = fp.a0_d;
            j0_d = fp.j0_d;
            xp   = fp.xp;
            vp   = fp.vp;
#ifdef INTEGRATE_6TH_SUN
            s0_s   = fp.s0_s;
            ap     = fp.ap;
#endif
#ifndef INTEGRATE_6TH_SUN
            a2_s = fp.a2_s;
            a3_s = fp.a3_s;
#else
            a3_s = fp.a3_s;
            a4_s = fp.a4_s;
            a5_s = fp.a5_s;
#endif
            a2_d = fp.a2_d;
            a3_d = fp.a3_d;
            
            time_c = fp.time_c;
            
            copyList(fp.n_list);
            copyHardList(fp.n_hard_list);
        }
        return *this;
    }
    FPHard &operator=(const FPGrav & fp){
        FPGrav::operator=(fp);
        if ( this != &fp ){
            x0   = v0   = 0.;
            a0_s = j0_s = 0.;
            a0_d = j0_d = 0.;
            xp   = vp   = 0.;
#ifdef INTEGRATE_6TH_SUN
            s0_s   = 0;
            ap     = 0;
#endif
#ifndef INTEGRATE_6TH_SUN
            a2_s = a3_s = 0.;
#else
            a3_s = a4_s = a5_s = 0.;
#endif
            a2_d = a3_d = 0.;
            
            time_c = fp.time;
            time = 0;

            clearList();
        }
        return *this;
    }

    void resetTime() {
        time = time_c + time;
        time_c = 0.;
    }
    PS::F64 getTime() const { return time_c + time; }

    void predict(PS::F64 Dt){
        x0   = pos;
        v0   = vel;
        a0_s = acc_s;
        j0_s = jerk_s;
#ifdef INTEGRATE_6TH_SUN
        s0_s = snap_s;
#endif
        j0_d = jerk_d;
        a0_d = acc_d;

#ifndef INTEGRATE_6TH_SUN
        PS::F64vec acc_sd  = acc_s  + acc_d;
        PS::F64vec jerk_sd = jerk_s + jerk_d;
        xp = pos + Dt*(vel    + Dt*(SECONDORDER*acc_sd + THIRDORDER*Dt*jerk_sd));
        vp = vel + Dt*(acc_sd + Dt* SECONDORDER*jerk_sd);
#else
        assert ( acc_ == acc_s + acc_d );
        PS::F64vec jerk_sd = jerk_s + jerk_d;
        PS::F64vec snap_sd = snap_s;
        xp = pos  + Dt*(vel     + Dt*(SECONDORDER*acc_    + Dt*(THIRDORDER*jerk_sd + Dt*FOURTHORDER*snap_sd)));
        vp = vel  + Dt*(acc_    + Dt*(SECONDORDER*jerk_sd + Dt* THIRDORDER*snap_sd));
        ap = acc_ + Dt*(jerk_sd + Dt* SECONDORDER*snap_sd);
#endif
    }
    void correct(PS::F64 Dt){
        PS::F64 Dt2 = Dt*Dt;
        PS::F64 Dt_inv  = 1./Dt;
        PS::F64 Dt_inv2 = Dt_inv *Dt_inv;
        PS::F64 Dt_inv3 = Dt_inv2*Dt_inv;
        
        PS::F64vec Am_s = acc_s - a0_s;
        PS::F64vec J0_s = Dt*j0_s;
        PS::F64vec J1_s = Dt*jerk_s;
#ifndef INTEGRATE_6TH_SUN 
        PS::F64vec A2_s =  3.*Am_s - (J1_s + 2.*J0_s);
        PS::F64vec A3_s = -2.*Am_s + (J1_s +    J0_s);

        a2_s = 2.*Dt_inv2*(A2_s + 3.*A3_s);
        a3_s = 6.*Dt_inv3* A3_s;
#else
        PS::F64vec S0_s = SECONDORDER*Dt2*s0_s;
        PS::F64vec S1_s = SECONDORDER*Dt2*snap_s;
        PS::F64vec A3_s =  10.*Am_s - (4.*J1_s + 6.*J0_s) + (   S1_s - 3.*S0_s);
        PS::F64vec A4_s = -15.*Am_s + (7.*J1_s + 8.*J0_s) - (2.*S1_s - 3.*S0_s);
        PS::F64vec A5_s =   6.*Am_s - 3.*(J1_s +    J0_s) + (   S1_s -    S0_s);

        PS::F64 Dt_inv4 = Dt_inv2*Dt_inv2;
        PS::F64 Dt_inv5 = Dt_inv3*Dt_inv2;

        a3_s =   3.*Dt_inv3*(A3_s + 4.*A4_s + 10.*A5_s);
        a4_s =  24.*Dt_inv4*(A4_s + 5.*A5_s);
        a5_s = 120.*Dt_inv5*A5_s;
#endif
        
        PS::F64vec Am_d = acc_d - a0_d;
        PS::F64vec J0_d = j0_d   * Dt;
        PS::F64vec J1_d = jerk_d * Dt;
        PS::F64vec A2_d =  3.*Am_d - (J1_d + 2.*J0_d);
        PS::F64vec A3_d = -2.*Am_d + (J1_d +    J0_d);

        a2_d = 2. * (A2_d + 3.*A3_d) * Dt_inv2;
        a3_d = 6. *  A3_d * Dt_inv3;

#ifndef INTEGRATE_6TH_SUN
        PS::F64vec A2_sd = A2_s+A2_d;
        PS::F64vec A3_sd = A3_s+A3_d;
        pos = xp + ONE_SIXTIETH*Dt2*(5.*A2_sd + 3.*ALPHA*A3_sd);
        vel = vp + ONE_TWELFTH *Dt *(4.*A2_sd + 3.      *A3_sd);
#else
        PS::F64vec A2_sd = A2_d;
        PS::F64vec A3_sd = A3_s+A3_d;
        PS::F64vec A4_sd = A4_s;
        PS::F64vec A5_sd = A5_s;
        pos  = xp + ONE_420TH   *Dt2*(35.*A2_sd + 21.*A3_sd + 14.*A4_sd + 10.*A5_sd);
        vel  = vp + ONE_SIXTIETH*Dt *(20.*A2_sd + 15.*A3_sd + 12.*A4_sd + 10.*A5_sd);
        assert( acc_ == acc_s + acc_d );
#endif
    }

    void calcDeltat(){
        PS::F64 dt_next = std::min(2.*dt, 0.5*dt_tree);

#ifndef INTEGRATE_6TH_SUN
        PS::F64 dt_1 = std::min(calcDt4th(eta,     alpha2, acc0, acc_d, jerk_d, a2_d, a3_d),
                                calcDt4th(eta_sun, alpha2,   0., acc_s, jerk_s, a2_s, a3_s));
        //PS::F64 dt_1 = calcDt4th(eta, alpha2, acc0, acc_d, jerk_d, a2_d+a3_d*dt, a3_d);
#else
        PS::F64 dt_1 = std::min(calcDt4th(eta,     alpha2, acc0, acc_d, jerk_d, a2_d,   a3_d),
                                calcDt6th(eta_sun, alpha2,   0., acc_s, jerk_s, snap_s, a3_s, a4_s, a5_s));
#endif
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
    PS::F64 ecc_rms_loc = 0.0;
    PS::F64 inc_rms_loc = 0.0;

#pragma omp parallel for reduction (+:ecc_rms_loc, inc_rms_loc)
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64 ax;
        PS::F64 ecc = pp[i].getEccentricity(ax);
        PS::F64vec h;
        PS::F64 inc = pp[i].getInclination(h);

        ecc_rms_loc += ecc*ecc;
        inc_rms_loc += inc*inc;
    }
    
    PS::F64 ecc_rms2 = PS::Comm::getSum(ecc_rms_loc)/n_tot;
    PS::F64 inc_rms2 = PS::Comm::getSum(inc_rms_loc)/n_tot;
    
    return sqrt(ecc_rms2 + inc_rms2);
}

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

#ifndef ISOTROPIC
#ifdef USE_INDIVIDUAL_CUTOFF

template <class Tpsys>
void setCutoffRadii(Tpsys & pp)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    const PS::S32 n_tot = pp.getNumberOfParticleGlobal();
    PS::F64 v_disp_k = calcVelDisp(pp, n_tot, n_loc);
    
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        pp[i].setROutRSearch(v_disp_k);
        pp[i].setRPlanet();
    }
}

#else //USE_INDIVIDUAL_CUTOFF

template <class Tpsys>
void setCutoffRadii(Tpsys & pp)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    const PS::S32 n_tot = pp.getNumberOfParticleGlobal();
    PS::F64 v_disp_k = calcVelDisp(pp, n_tot, n_loc);
    
    PS::F64 rHill_a_loc = 0.;
    PS::F64 v_kep_loc = 0.;
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64 ax      = pp[i].getSemimajorAxis2();
        PS::F64 rHill_a = pp[i].getRHill() * pow(ax,-FPGrav::p_cut);
        PS::F64 v_kep = pp[i].getKeplerVelocity();
        
#pragma omp critical
        {
            if ( rHill_a > rHill_a_loc ) rHill_a_loc = rHill_a;
            if ( v_kep > v_kep_loc ) v_kep_loc = v_kep;
        }
        pp[i].setRPlanet();
    }
    PS::F64 rHill_a_glb = PS::Comm::getMaxValue(rHill_a_loc);
    PS::F64 v_kep_glb = PS::Comm::getMaxValue(v_kep_loc);
    FPGrav::setROutRSearch(rHill_a_glb, v_kep_glb, v_disp_k);
}

#endif //USE_INDIVIDUAL_CUTOFF
#else //ISOTROPIC
#ifdef USE_INDIVIDUAL_CUTOFF

template <class Tpsys>
void setCutoffRadii(Tpsys & pp)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    const PS::S32 n_tot = pp.getNumberOfParticleGlobal();
    PS::F64 v_disp   = calcRandomVel(pp, n_tot, n_loc);
    
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        pp[i].setROutRSearch(v_disp);
        pp[i].setRPlanet();
    }
}

#else //USE_INDIVIDUAL_CUTOFF

void setCutoffRadii(Tpsys & pp)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    const PS::S32 n_tot = pp.getNumberOfParticleGlobal();
    PS::F64 v_disp   = calcRandomVel(pp, n_tot, n_loc);
    
    PS::F64 rHill_a_loc = 0.;
    PS::F64 v_kep_loc = 0.;
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64 ax      = pp[i].getSemimajorAxis2();
        PS::F64 rHill_a = pp[i].getRHill() * pow(ax,-FPGrav::p_cut);
#pragma omp critical
        {
            if ( rHill_a > rHill_a_loc ) rHill_a_loc = rHill_a;
        }
        pp[i].setRPlanet();
    }
    PS::F64 rHill_a_glb = PS::Comm::getMaxValue(rHill_a_loc);
    FPGrav::setROutRSearch(rHill_a_glb, v_disp);
}

#endif //USE_INDIVIDUAL_CUTOFF
#endif //ISOTROPIC
