#pragma once

#define SECONDORDER 5e-1
#define THIRDORDER  1.666666666666667e-1
#define FOURTHORDER 4.166666666666667e-2
#define FIFTHORDER  8.333333333333333e-3
#define ALPHA       1.166666666666667e-3

    
inline PS::F64 calcDt2nd(PS::F64 eta,
                         PS::F64 alpha2,
                         PS::F64 acc0,
                         PS::F64vec acc,
                         PS::F64vec jerk){
    PS::F64 Acc2  = acc*acc + alpha2*acc0*acc0;
    PS::F64 Jerk2 = jerk*jerk;
    PS::F64 dt2 = (Jerk2>0.) ? Acc2/Jerk2 : std::numeric_limits<double>::max();
    return eta * sqrt(dt2);
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
    PS::F64 dt2 = (Jerk>0.) ? (Acc*Snap + Jerk2)/(Jerk*Crac + Snap2) : std::numeric_limits<double>::max();
    return eta * sqrt(dt2);
}

class FPGrav : public EPGrav {
public:
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
        
        PS::F64 dt_1 = std::min(calcDt2nd(eta_0, alpha2, acc0, acc_d, jerk_d),
                                calcDt2nd(eta_0, alpha2,   0., acc_s, jerk_s));
        //PS::F64 dt_1 = calcDt2nd(eta_0, alpha2, acc0, acc_d, jerk_d);
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
    PS::F64vec a2_s;
    PS::F64vec a3_s;
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
        a0_s = j0_s = a0_d = j0_d = 0.;
        xp   = vp   = 0.;
        a2_s = a3_s = a2_d = a3_d = 0.;
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
        a2_s = fp.a2_s;
        a3_s = fp.a3_s;
        a2_d = fp.a2_d;
        a3_d = fp.a3_d;
        
        time_c = fp.time_c;
        
        copyList(fp.n_list);
        copyHardList(fp.n_hard_list);
    }
    FPHard(const FPGrav & fp) : FPGrav(fp){
        x0   = v0   = 0.;
        a0_s = j0_s = a0_d = j0_d = 0.;
        xp   = vp   = 0.;
        a2_s = a3_s = a2_d = a3_d = 0.;

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
            a2_s = fp.a2_s;
            a3_s = fp.a3_s;
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
            a0_s = j0_s = a0_d = j0_d = 0.;
            xp   = vp   = 0.;
            a2_s = a3_s = a2_d = a3_d = 0.;
            
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
        a0_d = acc_d;
        j0_s = jerk_s;
        j0_d = jerk_d;
        
        PS::F64vec acc_sd  = acc_s  + acc_d;
        PS::F64vec jerk_sd = jerk_s + jerk_d;
        xp = pos + Dt*(vel    + Dt*(SECONDORDER*acc_sd + THIRDORDER*Dt*jerk_sd));
        vp = vel + Dt*(acc_sd + Dt* SECONDORDER*jerk_sd);
    }
    void correct(PS::F64 Dt){
        PS::F64 Dt2 = Dt*Dt;
        PS::F64 Dt3 = Dt*Dt*Dt;
        PS::F64 Dt2inv = 1.0/Dt2;
        PS::F64 Dt3inv = 1.0/Dt3;
        a2_s = -Dt2inv * (6. *(a0_s-acc_s) +    Dt*(4.*j0_s+2.*jerk_s));
        a3_s =  Dt3inv * (12.*(a0_s-acc_s) + 6.*Dt*(   j0_s+   jerk_s));
        a2_d = -Dt2inv * (6. *(a0_d-acc_d) +    Dt*(4.*j0_d+2.*jerk_d));
        a3_d =  Dt3inv * (12.*(a0_d-acc_d) + 6.*Dt*(   j0_d+   jerk_d));
        
        PS::F64vec a2_sd = a2_s + a2_d;
        PS::F64vec a3_sd = a3_s + a3_d;
        pos = xp + Dt2*Dt2*(FOURTHORDER*a2_sd + ALPHA*FIFTHORDER *Dt*a3_sd);
        vel = vp +     Dt3*(THIRDORDER *a2_sd +       FOURTHORDER*Dt*a3_sd);
    }

    void calcDeltat(){
        PS::F64 dt_next = std::min(2.*dt, 0.5*dt_tree);

        PS::F64 dt_1 = std::min(calcDt4th(eta, alpha2, acc0, acc_d, jerk_d, a2_d+a3_d*dt, a3_d),
                                calcDt4th(eta, alpha2,   0., acc_s, jerk_s, a2_s+a3_s*dt, a3_s));
        //PS::F64 dt_1 = calcDt4th(eta, alpha2, acc0, acc_d, jerk_d, a2_d+a3_d*dt, a3_d);
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
