#pragma once


inline PS::F64 cutoff_f(PS::F64 y){
    const PS::F64 g = EPGrav::gamma;
    const PS::F64 g2 = EPGrav::g2;
    
    PS::F64 y2 = y*y;
    PS::F64 y3 = y2*y;
    PS::F64 y4 = y2*y2;
    PS::F64 y5 = y3*y2;
    PS::F64 y6 = y3*y3;
    PS::F64 y7 = y4*y3;
    
    PS::F64 f = 0.0;
    f =( -10.0*y7/3.0 +14.0*(g+1.0)*y6 -21.0*(g2 +3.0*g +1.0)*y5
         +y4*35.0*(g2*g +9.0*g2 +9.0*g +1.0)/3.0 -70.0*y3*(g2*g +3.0*g2 + g)
         +210.0*(g +1.0)*g2*y2 -140.0*g2*g*y*log(y) +g2*g2*(g2*g -7.0*g2 +21.0*g -35.0) )
        * EPGrav::g_1_inv7;

    return f;   
}


inline PS::F64 cutoff_W(PS::F64 rij, PS::F64 r_out_inv){
    const PS::F64 g = EPGrav::gamma;
    PS::F64 y = rij * r_out_inv;
    
    PS::F64 W = 0.0;
    if ( 1.0 <= y ) {
        W = 1.0;
    } else if ( y <= g ){
        W = y * EPGrav::w_y;
    } else {
        PS::F64 f1 = cutoff_f(1.);
        PS::F64 f = cutoff_f(y);
        W = f +y*(1.-f1);
    }

    return W;
}

inline PS::F64 cutoff_K(PS::F64 rij, PS::F64 r_out_inv){
    const PS::F64 g = FPGrav::gamma;
    
    PS::F64 y = rij * r_out_inv;
    PS::F64 x = (g - y) * EPGrav::g_1_inv;
    
    PS::F64 K = 0.;
    if( x < 0. ) {
        K = 0.;
    } else if ( x >= 1.0 ) {
        K = 1.;
    } else {
        PS::F64 x2 = x*x;
        PS::F64 x4 = x2*x2;
        K = (((-20.0*x+70.0)*x-84.0)*x+35.0)*x4;
    }
    return K;
}

inline PS::F64 cutoff_dK(PS::F64 rij, PS::F64 rij_inv, PS::F64 rijvij, PS::F64 r_out_inv){
    static const PS::F64 g = FPGrav::gamma;
    PS::F64 y = rij * r_out_inv;
    PS::F64 x = (g - y) * EPGrav::g_1_inv;
    PS::F64 dx = - rijvij * rij_inv * r_out_inv * EPGrav::g_1_inv;

    PS::F64 dK = 0.;
    if( x < 0. || x >= 1.0 ) {
        dK = 0.;
    } else {
        PS::F64 x2 = x*x;
        PS::F64 x3 = x2*x;
        PS::F64 x4 = x2*x2;
        PS::F64 x5 = x4*x;
        PS::F64 x6 = x4*x2;
        dK = (-140.0*x6 + 420.0*x5 - 420.0*x4 + 140.0*x3) * dx;
    }
    return dK;
}

inline PS::F64 cutoff_W2(PS::F64 dr2, PS::F64 r_out_inv){
    return cutoff_W(sqrt(dr2), r_out_inv);
}

inline PS::F64 cutoff_W2(PS::F64 dr2, PS::F64 r_out_inv1, PS::F64 r_out_inv2){
    return cutoff_W(sqrt(dr2), std::min(r_out_inv1, r_out_inv2));
}
