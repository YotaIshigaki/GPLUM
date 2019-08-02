#pragma once


inline PS::F64 cutoff_f(PS::F64 y){
    const PS::F64 g        = EPGrav::gamma;
    const PS::F64 g2       = g*g;
    
    return (((((((-10./3.*y + 14.*(g+1.))*y - 21.*((g+3.)*g+1.))*y
                + 35./3.*(((g+9.)*g+9.)*g+1.) )*y
               - 70.*((g+3.)*g+1.)*g )*y
              + 210.*(g+1.)*g2)*y - 140.*g2*g*log(y))*y
            + (((g-7.)*g+21.)*g-35.)*g2*g2 ) * EPGrav::g_1_inv7;
}


inline PS::F64 cutoff_W(PS::F64 rij, PS::F64 r_out_inv){
    const PS::F64 g = EPGrav::gamma;
    const PS::F64 y = rij * r_out_inv;
    
    if ( 1.0 <= y ) {
        return 1.0;
    } else if ( y <= g ){
        return y * EPGrav::w_y;
    } else {
        PS::F64 f = cutoff_f(y);
        return f +y*(1.-EPGrav::f1);
    }
}

inline PS::F64 cutoff_K(PS::F64 rij, PS::F64 r_out_inv){
    const PS::F64 g = EPGrav::gamma;
    
    const PS::F64 y = rij * r_out_inv;
    const PS::F64 x = (g - y) * EPGrav::g_1_inv;
    
    if( x < 0. ) {
        return 0.;
    } else if ( x >= 1.0 ) {
        return 1.;
    } else {
        PS::F64 x2 = x*x;
        PS::F64 x4 = x2*x2;
        return (((-20.*x +70.)*x -84.)*x +35.)*x4;
    }
}

inline PS::F64 cutoff_dK(PS::F64 rij, PS::F64 rij_inv, PS::F64 rijvij, PS::F64 r_out_inv){
    const PS::F64 g = EPGrav::gamma;
    const PS::F64 y = rij * r_out_inv;
    const PS::F64 x = (g - y) * EPGrav::g_1_inv;
    const PS::F64 dx = - rijvij * rij_inv * r_out_inv * EPGrav::g_1_inv;

    if( x < 0. || x >= 1.0 ) {
        return 0.;
    } else {
        return (((-140.*x +420.)*x -420.)*x +140.)*x*x*x * dx;
    }
}

inline PS::F64 cutoff_W2(PS::F64 dr2, PS::F64 r_out_inv){
    return cutoff_W(sqrt(dr2), r_out_inv);
}

inline PS::F64 cutoff_W2(PS::F64 dr2, PS::F64 r_out_inv1, PS::F64 r_out_inv2){
    return cutoff_W(sqrt(dr2), std::min(r_out_inv1, r_out_inv2));
}
