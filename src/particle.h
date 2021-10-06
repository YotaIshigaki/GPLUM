#pragma once

#define SAFTY_FACTOR 1.05


/////////////
/// Force ///
/////////////

class ForceGrav{
public:
    PS::F32vec acc;
    PS::F32    phi;
    PS::S32    neighbor;
    PS::S32    DUMMY_;

#ifdef FOR_PIKG01
    PS::S32    id_neighbor; 
    PS::S32    id_neighbor_dmmy; 
#else
    PS::S64    id_neighbor;        
#endif
    
    void clear(){
        acc         = 0.;
        phi         = 0.;
        neighbor    = 0;

        id_neighbor = -1;
    }
};


//////////////////////////
/// Essential Particle ///
//////////////////////////

class EPIGrav{
public:    
    PS::F64vec pos;     // position in cartesian
#ifdef USE_POLAR_COORDINATE
    PS::F64vec pos_pol; // position in polar
#endif

#ifdef FOR_PIKG01
    PS::S32 id;
    PS::S32 id_dmmy;
#else
    PS::S64 id;           // id number
#endif

#ifdef USE_INDIVIDUAL_CUTOFF
    PS::F64 r_out;        // cut-off radius
    PS::F64 r_search;     // search radius
#else
    static PS::F64 r_out;
    static PS::F64 r_search;
#endif

    PS::F64vec getPos() const {
#ifdef USE_POLAR_COORDINATE
        return pos_pol;
#else
        return pos;
#endif
    }
    PS::F64vec getPosCar() const { return pos; }
    PS::F64 getRSearch() const {
#ifdef USE_POLAR_COORDINATE
        return SAFTY_FACTOR * r_search / sqrt(pos*pos);
#else
        return SAFTY_FACTOR * r_search;
#endif
    }
    
    void copyFromFP(const EPIGrav & fp){
        pos      = fp.pos;
#ifdef USE_POLAR_COORDINATE
        pos_pol  = fp.pos_pol;
#endif

        id       = fp.id;

#ifdef USE_INDIVIDUAL_CUTOFF
        r_out    = fp.r_out;
        r_search = fp.r_search;
#endif
    }
};

#ifndef USE_INDIVIDUAL_CUTOFF
PS::F64 EPIGrav::r_out;
PS::F64 EPIGrav::r_search;
#endif


class EPJGrav : public EPIGrav {
public:
    PS::F64 mass;         // mass

    PS::F64vec vel;       // valocity
    PS::F64vec acc_d;     // acceleration for hard part

    PS::S32 id_local;     // local id number
    PS::S32 myrank;       // rank number

    PS::F64 getCharge() const { return mass; }

    void copyFromFP(const EPJGrav & fp){
        EPIGrav::copyFromFP(fp);
        mass     = fp.mass;
        
        vel      = fp.vel;
        acc_d    = fp.acc_d;

        id_local = fp.id_local;
        myrank   = fp.myrank;
    }
};


/////////////////////
/// Full Particle ///
/////////////////////

#define SECONDORDER  5.e-1
#define THIRDORDER   1.6666666666666667e-1
#define FOURTHORDER  4.1666666666666667e-2
#define FIFTHORDER   8.3333333333333333e-3
#define SIXTHORDER   1.3888888888888889e-3
#define SEVENTHORDER 1.9841269841269841e-4
#define ALPHA        1.1666666666666667e-3
#define ONE_TWELFTH  8.3333333333333333e-2
#define ONE_SIXTIETH 1.6666666666666667e-2
#define ONE_420TH    2.3809523809523810e-3

    
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


class FPGrav : public EPJGrav {
public:
    PS::F64vec acc;    // acceleration for soft part
    PS::F64vec acc_s;  // acceleration by sun
    PS::F64vec jerk_d; // jerk by planet
    PS::F64vec jerk_s; // jerk by sun
#ifdef INTEGRATE_6TH_SUN
    PS::F64vec acc_;
    PS::F64vec snap_s;
#endif
    PS::F64vec acc_gd; // acceleration by gas drag

    PS::F64 phi;       // potential for soft part
    PS::F64 phi_d;     // potential by planets
    PS::F64 phi_s;     // potential by sun
    static PS::F64 m_sun;
    static PS::F64 dens;

    static PS::F64 eps2;
    static PS::F64 eps2_sun;
    static PS::F64 R_cut0;
    static PS::F64 R_cut1;
    static PS::F64 R_search0;
    static PS::F64 R_search1;
#ifdef USE_RE_SEARCH_NEIGHBOR
    static PS::F64 R_search2;
    static PS::F64 R_search3;
#endif
    static PS::F64 gamma;
    static PS::F64 g_1_inv;  // 1/(g-1)
    static PS::F64 g_1_inv7; // 1/(g-1)^7
    static PS::F64 w_y;      // dW/dy if  y<g
    static PS::F64 f1;       // f(1;g)

    
#ifdef INDIRECT_TERM
    static PS::F64vec acc_indirect;
    static PS::F64vec pos_g;
    static PS::F64vec vel_g;
    static PS::F64    mass_tot;
#endif
    

#ifdef USE_INDIVIDUAL_CUTOFF
#ifndef CONSTANT_RANDOM_VELOCITY
    PS::F64 v_disp;
#else
    static PS::F64 v_disp;
#endif
#endif

#ifdef USE_INDIVIDUAL_CUTOFF
    PS::F64 r_out_inv;
#else
    static PS::F64 r_out_inv;
#endif
    
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
    PS::F64 f;
    static PS::F64 r_cut_min;
    static PS::F64 r_cut_max;
    static PS::F64 p_cut;
    static PS::F64 increase_factor;

    PS::S64 id_neighbor;
    
    PS::S64 id_cluster;
    PS::S32 n_cluster;

    PS::S32 neighbor;
    
    bool inDomain;
    bool isSent;
    bool isDead;
    bool isMerged;
#ifdef MERGE_BINARY
    bool isBinary;
    static PS::F64 R_merge;
#endif

    static void setGamma(PS::F64 g){
        gamma  = g;
        g_1_inv  = 1./(g - 1.);
        
        PS::F64 g2 = g*g;
        PS::F64 g_1_inv3 = g_1_inv * g_1_inv * g_1_inv;
        
        g_1_inv7 = g_1_inv3 * g_1_inv3 * g_1_inv;
        w_y = 7./3. * ((((((g- 9.)*g +45.)*g -60.*log(g))*g -45.)*g +9.)*g -1.) * g_1_inv7;
        f1 = (-10./3. + 14.*(g+1.) - 21.*((g+3.)*g+1.)
              + 35./3.*(((g+9.)*g+9.)*g+1.)
              - 70.*((g+3.)*g+1.)*g
              + 210.*(g+1.)*g2
              + (((g-7.)*g+21.)*g-35.)*g2*g2 ) * g_1_inv7;
    }
    PS::F64 getGamma() const{ return gamma; }
    PS::F64 getEps2() const{ return eps2;}
    PS::F64 getEps2_sun() const{ return eps2_sun;}

    PS::F64 getROut() const { return r_out; }
    PS::F64 getROut_inv() const { return r_out_inv; }

#ifdef USE_POLAR_COORDINATE
    void setPosPolar() {
        PS::F64 r = sqrt(pos*pos);
        pos_pol.x = atan2(pos.y, pos.x);
        pos_pol.y = log(r);
        pos_pol.z = asin(pos.z / r);
    }
#endif

    static PS::F64 getSolarMass() { return m_sun; }

#ifndef WITHOUT_SUN
    PS::F64 getSemimajorAxis() const {
#ifndef INDIRECT_TERM
        return 1.0 / (2.0/sqrt(pos*pos) - vel*vel/m_sun);
#else
        return 1.0 / (2.0/sqrt(pos*pos) - vel*vel/(m_sun+mass));
#endif
    }
    PS::F64 getSemimajorAxis2() const {
        PS::F64 ax;
        if ( getEccentricity(ax) < 0.6 ) {
            return ax;
        } else {
            return sqrt(pos*pos);
        }
    }
    PS::F64 getEccentricity(PS::F64 & ax) const {
        PS::F64 r = sqrt(pos*pos);
        PS::F64 rv = pos*vel;
#ifndef INDIRECT_TERM
        ax = 1.0 / (2.0/r - vel*vel/m_sun);
        PS::F64 ecccosu = 1. - r/ax;
        PS::F64 eccsinu2 = rv*rv/(m_sun*ax);
#else
        ax = 1.0 / (2.0/r - vel*vel/(m_sun+mass));
        PS::F64 ecccosu = 1. - r/ax;
        PS::F64 eccsinu2 = rv*rv/((m_sun+mass)*ax);
#endif
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
#endif // WITHOUT_SUN

#ifdef INTEGRATE_6TH_SUN
    void setAcc_() { acc_ = acc_s + acc_d; }
#endif
    
    
#ifdef USE_INDIVIDUAL_CUTOFF
    
    PS::F64 setROutRSearch(){
#ifndef WITHOUT_SUN
        PS::F64 rHill = getRHill();
        PS::F64 ax    = getSemimajorAxis2();

        PS::F64 r_out_i = std::max(R_cut0*pow(ax,-p_cut)*rHill, R_cut1*v_disp*dt_tree);
#else
        PS::F64 rHill = 0.;
        PS::F64 r_out_i = std::max(R_cut0*pow(mass * dt_tree * dt_tree, 1./3.), R_cut1*v_disp*dt_tree);
#endif
        
        r_out = std::max(r_out_i, r_cut_min);
        if ( r_cut_max > 0. ) r_out = std::min(r_cut_max, r_out);
        r_out_inv = 1. / r_out;    
        r_search  = R_search0*r_out + R_search1*v_disp*dt_tree;   
        assert ( r_out > 0. && r_search > 0. && r_search > r_out );
        
        return rHill;
    }
    
#else //USE_INDIVIDUAL_CUTOFF
    
    static void setROutRSearch(PS::F64 rHill_a_glb,
                               PS::F64 v_disp_glb){
      
#ifndef WITHOUT_SUN
        PS::F64 r_out_i = std::max(R_cut0*rHill_a_glb, R_cut1*v_disp_glb*dt_tree);
#else
        PS::F64 r_out_i = std::max(R_cut0*pow(rHill_a_glb * dt_tree * dt_tree, 1./3.), R_cut1*v_disp_glb*dt_tree);
#endif
        
        r_out = std::max(r_out_i, r_cut_min);
        if ( r_cut_max > 0. ) r_out = std::min(r_cut_max, r_out);
        r_out_inv = 1. / r_out;       
        r_search  = R_search0*r_out + R_search1*v_disp_glb*dt_tree;
        assert ( r_out > 0. && r_search > 0. );
    }
    
#endif //USE_INDIVIDUAL_CUTOFF
    
    void setRPlanet() {
        r_planet = pow(0.75*mass/(MY_PI*dens), 1./3.);
    }
    
    void copyFromForce(const ForceGrav & force){
        acc = force.acc;
        phi = force.phi;
        neighbor = force.neighbor;
        id_neighbor = force.id_neighbor;
    }

    void copy(const FPGrav & fp){
        EPJGrav::copyFromFP(fp);
        acc       = fp.acc;
        acc_s     = fp.acc_s;
        jerk_d    = fp.jerk_d;
        jerk_s    = fp.jerk_s;
#ifdef INTEGRATE_6TH_SUN
        acc_      = fp.acc_;
        snap_s    = fp.snap_s;
#endif
        acc_gd    = fp.acc_gd;
        
        phi       = fp.phi;
        phi_d     = fp.phi_d;
        phi_s     = fp.phi_s;

#ifdef USE_INDIVIDUAL_CUTOFF
#ifndef CONSTANT_RANDOM_VELOCITY
        v_disp    = fp.v_disp;
#endif
#endif

#ifdef USE_INDIVIDUAL_CUTOFF
        r_out_inv = fp.r_out_inv;
#endif

        time      = fp.time;
        dt        = fp.dt;
        acc0      = fp.acc0;
        r_planet  = fp.r_planet;
        f         = fp.f;

        id_neighbor = fp.id_neighbor;
        id_cluster  = fp.id_cluster;
        n_cluster   = fp.n_cluster;
        neighbor    = fp.neighbor;
    
        inDomain    = fp.inDomain;
        isSent      = fp.isSent;
        isDead      = fp.isDead;
        isMerged    = fp.isMerged;
#ifdef MERGE_BINARY
        isBinary    = fp.isBinary;
#endif
    }

    void dump(std::ostream & fout = std::cout) const {
        fout<<"id= "<<id<<std::endl;
        fout<<"mass= "<<mass<<std::endl;
        fout<<"pos= "<<pos<<std::endl;
#ifdef USE_POLAR_COORDINATE
        fout<<"pos_pol="<<pos_pol<<std::endl;
#endif
        fout<<"vel= "<<vel<<std::endl;
        fout<<"acc="<<acc<<std::endl;
        fout<<"phi="<<phi<<std::endl;
    }

    void readAscii(FILE* fp) {
        PS::S32 Flag;
        if ( !fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\n",
                     &this->id, &this->mass, &this->r_planet, &this->f, 
                     &this->pos.x, &this->pos.y, &this->pos.z,
                     &this->vel.x, &this->vel.y, &this->vel.z,
                     &this->neighbor,  &Flag) ) {
            //&this->r_out, &this->r_search) ) {
            errorMessage("The particle data have NOT been correctly read.");
            PS::Abort();
        }
#ifdef MERGE_BINARY
        isBinary = (bool)(Flag & (1<<0));
#endif
    }
    void writeAscii(FILE* fp) const {
        PS::S32 Flag = 0;
#ifdef MERGE_BINARY
        Flag |= ((PS::S32)isBinary)<<0;
#endif
        if ( !fprintf(fp, "%lld\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%20.15e\t%d\t%d\n",
                      this->id, this->mass, this->r_planet, this->f, 
                      this->pos.x, this->pos.y, this->pos.z,
                      this->vel.x, this->vel.y, this->vel.z,
                      this->neighbor, Flag) ) {
            //this->r_out, this->r_search) ){
            errorMessage("The particle data have NOT been correctly written.");
            PS::Abort();
        }   
    }
    void readBinary(FILE* fp) {
        FPGrav buf;
        if ( !fread(&buf, sizeof(buf), 1, fp) ) {
            errorMessage("The particle data have NOT been correctly read.");
            PS::Abort();
        }
        copy(buf);
    }
    void writeBinary(FILE* fp) const {
        FPGrav buf;
        buf.copy(*this);
        if ( !fwrite(&buf, sizeof(buf), 1, fp) ) {
            errorMessage("The particle data have NOT been correctly written.");
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
#ifdef AARSETH
        PS::F64 dt_1 = std::min(calcDt2nd(eta_0,    alpha2, acc0, acc_d, jerk_d),
                                calcDt2nd(eta_sun0, alpha2,   0., acc_s, jerk_s));
#else
        PS::F64 dt_1 = std::min(calcDt2nd(eta_0,    alpha2, acc0, acc_d, jerk_d),
                                calcDt3rd(eta_sun0, alpha2,   0., acc_s, snap_s));
#endif
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

PS::F64 FPGrav::m_sun     = 1.;
PS::F64 FPGrav::dens      = 5.049667e6;

PS::F64 FPGrav::eps2     = 0.;
PS::F64 FPGrav::eps2_sun = 0.;
PS::F64 FPGrav::R_cut0    = 2.;
PS::F64 FPGrav::R_cut1    = 8.;
PS::F64 FPGrav::R_search0 = 1.;
PS::F64 FPGrav::R_search1 = 4.;
#ifdef USE_RE_SEARCH_NEIGHBOR
PS::F64 FPGrav::R_search2 = 1.;
PS::F64 FPGrav::R_search3 = 4.;
#endif
PS::F64 FPGrav::gamma     = 0.5;
PS::F64 FPGrav::g_1_inv   = -2.;   // 1/(gamma-1)
PS::F64 FPGrav::g_1_inv7  = -128.; // 1/(gamma-1)^7
PS::F64 FPGrav::w_y;               // dW/dy when y<g
PS::F64 FPGrav::f1;                // f(1;g)

#ifdef INDIRECT_TERM
PS::F64vec FPGrav::acc_indirect = 0.;
PS::F64vec FPGrav::pos_g        = 0.;
PS::F64vec FPGrav::vel_g        = 0.;
PS::F64    FPGrav::mass_tot     = 0.;
#endif
#ifdef CONSTANT_RANDOM_VELOCITY
PS::F64 FPGrav::v_disp    = 0.;
#endif

#ifndef USE_INDIVIDUAL_CUTOFF
PS::F64 FPGrav::r_out_inv;
#endif

PS::F64 FPGrav::dt_tree   = pow2(-5);
PS::F64 FPGrav::dt_min    = pow2(-13);
PS::F64 FPGrav::eta       = 0.01;
PS::F64 FPGrav::eta_0     = 0.001;
PS::F64 FPGrav::eta_sun   = 0.01;
PS::F64 FPGrav::eta_sun0  = 0.001;
PS::F64 FPGrav::alpha2    = 1.;

PS::F64 FPGrav::r_cut_min = 0.;
PS::F64 FPGrav::r_cut_max = 0.;
PS::F64 FPGrav::p_cut     = 0.;
PS::F64 FPGrav::increase_factor = 1.;

#ifdef MERGE_BINARY
PS::F64 FPGrav::R_merge   = 0.2;
#endif



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
    
    std::vector<PS::S64> n_list;
    std::vector<PS::S32> n_hard_list;
    
    void clearList(){
        //std::vector<PS::S32> tmp0, tmp1;
        //tmp0.swap(n_list);
        //tmp1.swap(n_hard_list);
        n_list.clear();
        n_hard_list.clear();
    }
    void copyList(std::vector<PS::S64> list){
        n_list.resize(list.size());
        std::copy(list.begin(),list.end(),n_list.begin());
    }
    void copyList(PS::S64 * list){
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
    void makeHardList(std::map<PS::S64,PS::S32> & id_map){
        n_hard_list.clear();
        n_hard_list.reserve(neighbor);
        for ( PS::S32 i=0; i<neighbor; i++){
            n_hard_list.push_back(id_map.at(n_list.at(i)));
        }
    }
    template <class Tpsys>
    void makeHardList(std::map<PS::S64,PS::S32> & id_map,
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

        r_planet = f = 0.;
        
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

#ifndef WITHOUT_SUN
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
#else //WITHOUT_SUN
#ifndef INTEGRATE_6TH_SUN
        PS::F64vec A2_s = 0.;
        PS::F64vec A3_s = 0.;
        a2_s = a3_s = 0.;
#else
        PS::F64vec A3_s = 0.;
        PS::F64vec A4_s = 0.;
        PS::F64vec A5_s = 0.;
        a3_s = a4_s = a5_s = 0.;
#endif  
#endif //WITHOUT_SUN
        
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
#ifndef WITHOUT_SUN
        PS::F64 dt_1 = std::min(calcDt4th(eta,     alpha2, acc0, acc_d, jerk_d, a2_d, a3_d),
                                calcDt4th(eta_sun, alpha2,   0., acc_s, jerk_s, a2_s, a3_s));
        //PS::F64 dt_1 = calcDt4th(eta, alpha2, acc0, acc_d, jerk_d, a2_d+a3_d*dt, a3_d);
#else
        PS::F64 dt_1 = calcDt4th(eta,     alpha2, acc0, acc_d, jerk_d, a2_d, a3_d);
#endif
#else
#ifdef AARSETH
#ifndef WITHOUT_SUN
        PS::F64 dt_1 = std::min(calcDt4th(eta,     alpha2, acc0, acc_d, jerk_d, a2_d,   a3_d),
                                calcDt4th(eta_sun, alpha2,   0., acc_s, jerk_s, snap_s, a3_s));
#else
        PS::F64 dt_1 = calcDt4th(eta,     alpha2, acc0, acc_d, jerk_d, a2_d,   a3_d);
#endif
#else
#ifndef WITHOUT_SUN
        PS::F64 dt_1 = std::min(calcDt4th(eta,     alpha2, acc0, acc_d, jerk_d, a2_d,   a3_d),
                                calcDt6th(eta_sun, alpha2,   0., acc_s, jerk_s, snap_s, a3_s, a4_s, a5_s));
#else
        PS::F64 dt_1 = calcDt4th(eta,     alpha2, acc0, acc_d, jerk_d, a2_d,   a3_d);
#endif
#endif
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
void calcRandomVel(Tpsys & pp,
                   const PS::S32 n_tot,
                   const PS::S32 n_loc)
{
    PS::F64 r_max=0., r_min=-std::numeric_limits<PS::F64>::max();
#pragma omp parallel for reduction (max: r_max, r_min)
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64vec pos = pp[i].pos;
        PS::F64    r2  = pos.x*pos.x + pos.y*pos.y;
        r_max = std::max(r_max, r2);
        r_min = std::max(r_min, -r2);
    }
    r_max = sqrt(r_max)  * 1.01;
    r_min = sqrt(-r_min) * 0.99;
    r_max = PS::Comm::getMaxValue(r_max);
    r_min = PS::Comm::getMinValue(r_min);

    const PS::S32 N = 32;
    const PS::F64 dr    = (r_max - r_min) /N;
    const PS::F64 drinv = 1./dr;
#ifdef ISOTROPIC
    PS::F64vec v_ave_loc[N];
    PS::F64vec v_ave_glb[N];
    PS::F64    v_sq_loc[N];
    PS::F64    v_sq_glb[N];
#else
    PS::F64 v_disp_loc[N];
#endif
    PS::F64 v_disp_glb[N];
    PS::S32 n_ptcl_loc[N];
    PS::S32 n_ptcl_glb[N];

#ifdef ISOTROPIC
    for(PS::S32 i=0; i<N; i++) {
        v_ave_loc[i] = 0.;
        v_sq_loc[i] = 0.;
        n_ptcl_loc[i] = 0;
    }
    //#pragma omp parallel for reduction (+:v_ave_loc[:N], v_sq_loc[:N], n_ptcl_loc[:N])
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64vec pos = pp[i].pos;
        PS::F64    r2  = pos.x*pos.x + pos.y*pos.y;
        PS::S32 j = (PS::S32)((sqrt(r2) - r_min) * drinv);
        if ( j == N ) j = N - 1;
        assert ( 0<= j && j < N );
        PS::F64vec vec = pp[i].vel;
        v_ave_loc[j] += vec;
        v_sq_loc[j]  += vec*vec;
        n_ptcl_loc[j] ++;
    }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Allreduce(v_ave_loc,  v_ave_glb,  N, PS::GetDataType(*v_ave_loc),  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(v_sq_loc,   v_sq_glb,   N, PS::GetDataType(*v_sq_loc),   MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(n_ptcl_loc, n_ptcl_glb, N, PS::GetDataType(*n_ptcl_loc), MPI_SUM, MPI_COMM_WORLD);
#else
    for(PS::S32 i=0; i<N; i++) {
        v_ave_glb[i]  = v_ave_loc[i];
        v_sq_glb[i]   = v_sq_loc[i];
        n_ptcl_glb[i] = n_ptcl_loc[i];
    }
#endif

    PS::S32 n_tot0 = 0;
    for(PS::S32 i=0; i<N; i++) {
	v_disp_glb[i] = (n_ptcl_glb[i] > 1) ?
	    sqrt(v_sq_glb[i] / n_ptcl_glb[i]
		 - v_ave_glb[i]*v_ave_glb[i]/( n_ptcl_glb[i]* n_ptcl_glb[i]) ) : 0.;
        n_tot0 += n_ptcl_glb[i];
    }
    assert ( n_tot == n_tot0 ); 
    
#else //ISOTROPIC
    for(PS::S32 i=0; i<N; i++) {
        v_disp_loc[i] = 0.;
        n_ptcl_loc[i] = 0;
    }
    #pragma omp parallel for reduction (+:v_disp_loc[:N], n_ptcl_loc[:N])
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64vec pos = pp[i].pos;
        PS::F64    r2  = pos.x*pos.x + pos.y*pos.y;
        PS::F64    ri  = sqrt(r2);
        PS::S32 j = (PS::S32)((ri - r_min) * drinv);
        if ( j == N ) j = N - 1;
        assert ( 0<= j && j < N );
#if 1
        PS::F64vec v_kep;
        v_kep.x=-pos.y/ri ; v_kep.y=pos.x/ri; v_kep.z=0.;
        v_kep *= pp[i].getKeplerVelocity();
        PS::F64vec v_ran = pp[i].vel - v_kep;
#else
        PS::F64 ax;
        PS::F64 ecc = pp[i].getEccentricity(ax);
        PS::F64vec h;
        PS::F64 inc = pp[i].getInclination(h);
        PS::F64vec v_kep = pp[i].getKeplerVelocity();
        PS::F64    v_ran = (ecc*ecc + inc*inc) * v_kep*v_kep;
#endif
        v_disp_loc[j] += v_ran * v_ran;
        n_ptcl_loc[j] ++;
    }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Allreduce(v_disp_loc, v_disp_glb, N, PS::GetDataType(*v_disp_loc), MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(n_ptcl_loc, n_ptcl_glb, N, PS::GetDataType(*n_ptcl_loc), MPI_SUM, MPI_COMM_WORLD);
#else
    for(PS::S32 i=0; i<N; i++) {
        v_disp_glb[i] = v_disp_loc[i];
        n_ptcl_glb[i] = n_ptcl_loc[i];
    }
#endif

    PS::S32 n_tot0 = 0;
    for(PS::S32 i=0; i<N; i++) {
        v_disp_glb[i] = (n_ptcl_glb[i] > 0) ? sqrt(v_disp_glb[i] / n_ptcl_glb[i]) : 0.;
        n_tot0 += n_ptcl_glb[i];
    }
    assert ( n_tot == n_tot0 ); 
    
#endif //ISOTROPIC

#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64vec pos = pp[i].pos;
        PS::F64    r2  = pos.x*pos.x + pos.y*pos.y;
        PS::F64    v_dispi = 0.;
        PS::F64    ni      = 0.;
        for(PS::S32 j=0; j<N; j++){
            PS::F64 ddr  = r_min + (j+0.5)*dr - sqrt(r2);
            PS::F64 expr = exp(-ddr*ddr * drinv*drinv);
            v_dispi += v_disp_glb[j] * expr;
            ni      += (n_ptcl_glb[j] > 0) ? expr : 0.;
        }
        pp[i].v_disp = v_dispi / ni;
        assert ( pp[i].v_disp >= 0. ); 
    } 
}

template <class Tpsys>
PS::F64 calcRandomVelAll(Tpsys & pp,
                         const PS::S32 n_tot,
                         const PS::S32 n_loc)
{
#ifdef ISOTROPIC
    PS::F64vec v_ave_loc = 0.;
    PS::F64vec v_ave_glb = 0.;
    PS::F64    v_sq_loc = 0.;
    PS::F64    v_sq_glb = 0.;
#else
    PS::F64 v_disp_loc = 0.;
#endif
    PS::F64 v_disp_glb = 0.;
    PS::S32 n_ptcl_loc = 0;
    PS::S32 n_ptcl_glb = 0;

#ifdef ISOTROPIC
#pragma omp parallel for reduction (+:v_ave_loc, v_sq_loc, n_ptcl_loc)
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64vec vec = pp[i].vel;
        v_ave_loc += vec;
        v_sq_loc  += vec*vec;
        n_ptcl_loc ++;
    }

    v_ave_glb  = PS::Comm::getSum(v_ave_loc);
    v_sq_glb   = PS::Comm::getSum(v_sq_loc);
    n_ptcl_glb = PS::Comm::getSum(n_ptcl_loc);
    assert ( n_ptcl_glb == n_tot );
    
    v_disp_glb =  sqrt(v_sq_glb / n_ptcl_glb - v_ave_glb*v_ave_glb);
    
#else //ISOTROPIC  
#pragma omp parallel for reduction (+:v_disp_loc, n_ptcl_loc)
    for(PS::S32 i=0; i<n_loc; i++){
#if 1
        PS::F64vec pos = pp[i].pos;
        PS::F64    r2  = pos.x*pos.x + pos.y*pos.y;
        PS::F64    ri  = sqrt(r2);
        PS::F64vec v_kep;
        v_kep.x=-pos.y/ri ; v_kep.y=pos.x/ri; v_kep.z=0.;
        v_kep *= pp[i].getKeplerVelocity();
        PS::F64vec v_ran = pp[i].vel - v_kep;
#else
        PS::F64 ax;
        PS::F64 ecc = pp[i].getEccentricity(ax);
        PS::F64vec h;
        PS::F64 inc = pp[i].getInclination(h);
        PS::F64vec v_kep = pp[i].getKeplerVelocity();
        PS::F64    v_ran = (ecc*ecc + inc*inc) * v_kep*v_kep;
#endif
        v_disp_loc += v_ran * v_ran;
        n_ptcl_loc ++;
    }

    v_disp_glb = PS::Comm::getSum(v_disp_loc);
    n_ptcl_glb = PS::Comm::getSum(n_ptcl_loc);
    assert ( n_ptcl_glb == n_tot );
    
    v_disp_glb = sqrt(v_disp_glb / n_ptcl_glb);
    
#endif //ISOTROPIC
    
    return v_disp_glb;
}

#ifdef USE_INDIVIDUAL_CUTOFF

template <class Tpsys>
void setCutoffRadii(Tpsys & pp)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    const PS::S32 n_tot = pp.getNumberOfParticleGlobal();

#ifndef CONSTANT_RANDOM_VELOCITY
    calcRandomVel(pp, n_tot, n_loc);
#endif
    
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        pp[i].setROutRSearch();
        //pp[i].setRPlanet();
    }
}

#else //USE_INDIVIDUAL_CUTOFF

template <class Tpsys>
void setCutoffRadii(Tpsys & pp)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal();
    const PS::S32 n_tot = pp.getNumberOfParticleGlobal();

#ifndef CONSTANT_RANDOM_VELOCITY
    PS::F64 v_disp_glb = calcRandomVelAll(pp, n_tot, n_loc);
#else
    PS::F64 v_disp_glb = v_disp;
#endif
    
    PS::F64 rHill_a_loc = 0.;
#pragma omp parallel for reduction (max: rHill_a_loc)
    for(PS::S32 i=0; i<n_loc; i++){
#ifndef WITHOUT_SUN
        PS::F64 ax      = pp[i].getSemimajorAxis2();
        PS::F64 rHill_a = pp[i].getRHill() * pow(ax,-FPGrav::p_cut);
        
        rHill_a_loc = std::max(rHill_a, rHill_a_loc);
#else
        rHill_a_loc = std::max(pp[i].mass, rHill_a_loc);
#endif
        
        //pp[i].setRPlanet();
    }
    PS::F64 rHill_a_glb = PS::Comm::getMaxValue(rHill_a_loc);
    FPGrav::setROutRSearch(rHill_a_glb, v_disp_glb);
}

#endif //USE_INDIVIDUAL_CUTOFF



//////////////////////
/// Super Particle ///
//////////////////////

class MyMomentMonopole : public PS::MomentMonopole {
public:
#ifdef USE_POLAR_COORDINATE
    PS::F64vec pos_car;
    
    MyMomentMonopole() {
        pos_car = 0.;
    }
    MyMomentMonopole(const PS::F64 m, const PS::F64vec & p, const PS::F64vec & p_car) : PS::MomentMonopole(m, p) {
        pos_car = p_car;
    }
    void init(){
        PS::MomentMonopole::init();
        pos_car = 0.;
    }
    template<class Tepj>
    void accumulateAtLeaf(const Tepj & epj){
        PS::MomentMonopole::accumulateAtLeaf(epj);
        pos_car += epj.getCharge() * epj.getPosCar();
    }
    template<class Tepj>
    void accumulateAtLeaf2(const Tepj & epj){}
    void accumulate(const MyMomentMonopole & mom){
        PS::MomentMonopole::accumulate(mom);
        pos_car += mom.mass * mom.pos_car;
    }
    void accumulate2(const MyMomentMonopole & mom){}
    // for DEBUG 
    void dump(std::ostream & fout = std::cout) const {
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"pos_car="<<pos_car<<std::endl;
    }
#endif
    void set(){
        pos     = (mass != 0.) ? pos / mass : 0.;
#ifdef USE_POLAR_COORDINATE
        pos_car = (mass != 0.) ? pos_car / mass : 0.;
#endif
    }
};


class MySPJMonopole : public PS::SPJMonopole {
public:
#ifdef USE_POLAR_COORDINATE
    PS::F64vec pos_car;
    
    template<class Tmom>
    void copyFromMoment(const Tmom & mom){
        PS::SPJMonopole::copyFromMoment(mom);
        this->pos_car  = mom.pos_car;
    }
    void clear(){
        PS::SPJMonopole::clear();
        pos_car = 0.0;
    }
    PS::F64vec getPosCar() const { return pos_car; }
    MyMomentMonopole convertToMoment() const {
        return MyMomentMonopole(mass, pos, pos_car);
    }
#endif
};

class MyMomentQuadrupole : public PS::MomentQuadrupole {
public:
#ifdef USE_POLAR_COORDINATE
    PS::F64vec pos_car;
    
    MyMomentQuadrupole(){
        pos_car = 0.;
    }
    MyMomentQuadrupole(const PS::F64 m, const PS::F64vec & p, const PS::F64mat & q, const PS::F64vec & p_car) : PS::MomentQuadrupole(m, p, q) {
        pos_car = p_car;
    }
    void init(){
        PS::MomentQuadrupole::init();
        pos_car = 0.;
    }
    template<class Tepj>
    void accumulateAtLeaf(const Tepj & epj){
        PS::MomentQuadrupole::accumulateAtLeaf(epj);
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
    void accumulate(const MyMomentQuadrupole & mom){
        PS::MomentQuadrupole::accumulate(mom);
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
        fout<<"pos_car="<<pos_car<<std::endl;
        fout<<"quad= "<<quad<<std::endl;
    }
#endif
    void set(){
        pos     = (mass != 0.) ? pos / mass : 0.;
#ifdef USE_POLAR_COORDINATE
        pos_car = (mass != 0.) ? pos_car / mass : 0.;
#endif
    }
};

class MySPJQuadrupole : public PS::SPJQuadrupole {
public:
#ifdef USE_POLAR_COORDINATE
    PS::F64vec pos_car;

    PS::F64 getCharge() const { return mass; }
    PS::F64vec getPos() const { return pos; }
    PS::F64vec getPosCar() const { return pos_car; }
    void copyFromMoment(const MyMomentQuadrupole & mom){
        PS::SPJQuadrupole::copyFromMoment(mom);
        pos_car = mom.pos_car;
    }
    MyMomentQuadrupole convertToMoment() const {
        return MyMomentQuadrupole(mass, pos, quad, pos_car);
    }
    void clear(){
        PS::SPJQuadrupole::clear();
        pos_car = 0.0;
    }
#endif
};

