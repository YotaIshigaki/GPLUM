#pragma once

class ForceGrav{
public:
    PS::F64vec acc;
    PS::F64    phi;
    PS::S32 neighbor;

    PS::S32 id_neighbor;
    
    void clear(){
        acc = 0.;
        phi = 0.;
        neighbor = 0;

        id_neighbor = -1;
    }
};


class EPGrav{
public:
    PS::F64vec pos;   // position
    PS::F64vec vel;   // velocity
    PS::F64vec acc;   // acceleration for soft part
    PS::F64vec acc_s; // acceleration of sun for hard part
    PS::F64vec acc_d; // acceleration of planets for hard part

    PS::F64 mass;     // mass
    
    PS::S32 id;       // id number
    PS::S32 id_local; // local id number
    PS::S32 myrank;   // MPI rank number

    PS::S32 neighbor; // the number of neighbors

#ifdef USE_INDIVIDUAL_CUTOFF
    PS::F64 r_out;
    PS::F64 r_out_inv;
    PS::F64 r_search;
#else //USE_INDIVIDUAL_CUTOFF
    static PS::F64 r_out;
    static PS::F64 r_out_inv;
    static PS::F64 r_search;
#endif //USE_INDIVIDUAL_CUTOFF

    static PS::F64 eps2;
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
    
    PS::F64 getROut() const {
#ifdef TEST_PTCL
        static PS::F64 R_OUT_MIN = sqrt((PS::F64)std::numeric_limits<PS::F32>::min());
        return ( r_out == 0. ) ? R_OUT_MIN : r_out;    
#else
        return r_out;
#endif
    }
    PS::F64 getROut_inv() const { return r_out_inv; }
    PS::F64 getRSearch() const { return r_search; }
    
    PS::F64vec getPos() const { return pos; }
    PS::F64 getCharge() const {
#ifdef TEST_PTCL
        static PS::F64 MASS_MIN = (PS::F64)std::numeric_limits<PS::F32>::min();
        return ( mass == 0. ) ? MASS_MIN : mass;      
#else
        return mass;
#endif
    }
    PS::F64 getIdReal() const { return (PS::F64)id; }
    
    void copyFromFP(const EPGrav & fp){
        pos   = fp.pos;
        vel   = fp.vel;
        acc   = fp.acc;
        acc_s = fp.acc_s;
        acc_d = fp.acc_d;
    
        mass = fp.mass;
    
        id = fp.id;
        id_local = fp.id_local;
        myrank = fp.myrank;
    
#ifdef USE_INDIVIDUAL_CUTOFF
        r_out     = fp.r_out;
        r_out_inv = fp.r_out_inv;
        r_search  = fp.r_search;
#endif
    }
};
