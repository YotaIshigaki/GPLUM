#pragma once


inline PS::F64 cutoff_f(PS::F64 y){
    const PS::F64 g        = FPGrav::gamma;
    const PS::F64 g2       = g*g;
    
    return (((((((-10./3.*y + 14.*(g+1.))*y - 21.*((g+3.)*g+1.))*y
                + 35./3.*(((g+9.)*g+9.)*g+1.) )*y
               - 70.*((g+3.)*g+1.)*g )*y
              + 210.*(g+1.)*g2)*y - 140.*g2*g*log(y))*y
            + (((g-7.)*g+21.)*g-35.)*g2*g2 ) * FPGrav::g_1_inv7;
}


inline PS::F64 cutoff_W(PS::F64 rij, PS::F64 r_out_inv){
    const PS::F64 g = FPGrav::gamma;
    const PS::F64 y = rij * r_out_inv;
    
    if ( 1.0 <= y ) {
        return 1.0;
    } else if ( y <= g ){
        return y * FPGrav::w_y;
    } else {
        PS::F64 f = cutoff_f(y);
        return f +y*(1.-FPGrav::f1);
    }
}

inline PS::F64 cutoff_K(PS::F64 x){
    if( x < 0. ) {
        return 0.;
    } else if ( x >= 1.) {
        return 1.;
    } else {
        PS::F64 x2 = x*x;
        return (((-20.*x +70.)*x -84.)*x +35.)*x2*x2;
    }
}
inline PS::F64 cutoff_dKdr(PS::F64 x, PS::F64 r_out_inv){
    const PS::F64 x_1   = x - 1.;
    return ( x < 0. || x>= 1. ) ? 0.
        : ( 140.* x*x*x * x_1*x_1*x_1
            * r_out_inv*FPGrav::g_1_inv );
}
inline PS::F64 cutoff_d2Kdr2(PS::F64 x, PS::F64 r_out_inv){
    const PS::F64 x_1   = x - 1.;
    return ( x < 0. || x>= 1. ) ? 0.
        : ( 420.* x*x * x_1*x_1 * (1.-2.*x)
            * r_out_inv*r_out_inv*FPGrav::g_1_inv*FPGrav::g_1_inv );
}
inline PS::F64 cutoff_d3Kdr3(PS::F64 x, PS::F64 r_out_inv){
    const PS::F64 x_1   = x - 1.;
    return ( x < 0. || x>= 1. ) ? 0.
        : ( 840.* x * x_1 * (1.+5.*x*x_1)
            * r_out_inv*r_out_inv*r_out_inv*FPGrav::g_1_inv*FPGrav::g_1_inv*FPGrav::g_1_inv );
}

inline PS::F64 cutoff_K(PS::F64 rij, PS::F64 r_out_inv){
    const PS::F64 x = ( FPGrav::gamma - rij*r_out_inv ) * FPGrav::g_1_inv;
    return cutoff_K(x);
}
inline PS::F64 cutoff_dKdt(PS::F64 rij, PS::F64 r_out_inv, PS::F64 alpha){
    const PS::F64 x = ( FPGrav::gamma - rij*r_out_inv ) * FPGrav::g_1_inv;
    return alpha * rij * cutoff_dKdr(x, r_out_inv);
}
inline PS::F64 cutoff_d2Kdt2(PS::F64 rij, PS::F64 r_out_inv,
                             PS::F64 alpha, PS::F64 beta,
                             PS::F64 & dKdt){
    const PS::F64 x = ( FPGrav::gamma - rij*r_out_inv ) * FPGrav::g_1_inv;
    const PS::F64 dKdr   = cutoff_dKdr(x, r_out_inv);
    const PS::F64 alpha2 = alpha * alpha;
    dKdt = alpha * rij * dKdr;
    return rij * ( (beta + 4.*alpha2) * dKdr
                   + rij * alpha2 * cutoff_d2Kdr2(x, r_out_inv) );
}
inline PS::F64 cutoff_d3Kdt3(PS::F64 rij, PS::F64 r_out_inv,
                             PS::F64 alpha, PS::F64 beta, PS::F64 gamma,
                             PS::F64 & dKdt, PS::F64 & d2Kdt2){
    const PS::F64 x = ( FPGrav::gamma - rij*r_out_inv ) * FPGrav::g_1_inv;
    const PS::F64 dKdr   = cutoff_dKdr  (x, r_out_inv);
    const PS::F64 d2Kdr2 = cutoff_d2Kdr2(x, r_out_inv);
    const PS::F64 alpha2 = alpha * alpha;
    const PS::F64 beta_4alpha2  = beta + 4.*alpha2;
    dKdt   = alpha * rij * dKdr;
    d2Kdt2 = rij * ( beta_4alpha2 * dKdr + alpha2 * rij * d2Kdr2 );
    return rij * ( ( gamma + 4.*alpha*(3.*beta+7.*alpha2) ) * dKdr
                   + rij * ( 3.*alpha*beta_4alpha2 * d2Kdr2
                             + rij * alpha2*alpha * cutoff_d3Kdr3(x, r_out_inv) ) );
}

inline PS::F64 cutoff_W2(PS::F64 dr2, PS::F64 r_out_inv){
    return cutoff_W(sqrt(dr2), r_out_inv);
}

inline PS::F64 cutoff_W2(PS::F64 dr2, PS::F64 r_out_inv1, PS::F64 r_out_inv2){
    return cutoff_W(sqrt(dr2), std::min(r_out_inv1, r_out_inv2));
}
