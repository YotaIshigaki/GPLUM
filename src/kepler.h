#pragma once

PS::F64 KeplerEq(PS::F64 u,
                 PS::F64 ecc){ return u - ecc * sin(u); }
PS::F64 solveKeplerEq(PS::F64 l,
                      PS::F64 ecc)
{
    PS::F64 u;
    
    PS::F64 ecc2 = ecc*ecc;
    PS::F64 ecc3 = ecc2*ecc;
    PS::F64 ecc4 = ecc2*ecc2;
    PS::F64 ecc5 = ecc3*ecc2;
    PS::F64 ecc6 = ecc3*ecc3;
    u = l
        + (ecc     - ecc3/8. + ecc5/192.)*sin(l)
        + (ecc2/2. - ecc4/6. + ecc6/48. )*sin(2.*l)
        + (3.*ecc3/8. - 27.*ecc5/128.)*sin(3.*l)
        + (   ecc4/3. -  4.*ecc6/15. )*sin(4.*l)
        + 125.*ecc5/384.*sin(5.*l)
        +  27.*ecc6/ 80.*sin(6.*l);
#if 1
    if ( abs(KeplerEq(u,ecc)-l) > 1.e-15  ){
        PS::F64 u0;
        PS::S32 loop = 0;
        //u = l;
        do {
            u0 = u;
            PS::F64 sinu0 = sin(u0);
            PS::F64 cosu0 = cos(u0);
            u = u0 - ((u0 - ecc*sinu0 - l)/(1. - ecc*cosu0));
            loop++;
        } while( fabs(u - u0) > 1.e-15 && loop < 10 );
    }
#endif

    return u;
}
void posVel2OrbitalElement(PS::F64vec pos,
                           PS::F64vec vel,
                           PS::F64 mu,
                           PS::F64 & ax,
                           PS::F64 & ecc,
                           PS::F64 & n,
                           PS::F64 & u,
                           PS::F64vec & P,
                           PS::F64vec & Q)
{
    PS::F64 r2 = pos*pos;
    PS::F64 r = sqrt(r2);
    PS::F64 rinv = 1./r;
    PS::F64 v2 = vel*vel;
    PS::F64 rv = pos*vel;
    
    ax = 1.0 / (2.0*rinv - v2 / mu);
    PS::F64 ecccosu = 1. - r/ax;
    PS::F64 eccsinu = rv/sqrt(mu*ax);
    ecc = sqrt(ecccosu*ecccosu + eccsinu*eccsinu);
    n = sqrt(mu / (ax*ax*ax));
    
    u = atan2(eccsinu, ecccosu);
    
    PS::F64 cosu = ecccosu/ecc;
    PS::F64 sinu = eccsinu/ecc;
    PS::F64 aninv = sqrt(ax/mu);
    PS::F64 ecc_sq = sqrt(1.-ecc*ecc);
    P = rinv*cosu * pos - aninv*sinu * vel;
    Q = (rinv*sinu * pos + aninv*(cosu-ecc) * vel )/ecc_sq;
}
void orbitalElement2PosVel(PS::F64vec & pos,
                           PS::F64vec & vel,
                           PS::F64 mu,
                           PS::F64 ax,
                           PS::F64 ecc,
                           PS::F64 n,
                           PS::F64 u,
                           PS::F64vec P,
                           PS::F64vec Q)
{
    PS::F64 cosu = cos(u);
    PS::F64 sinu = sin(u);
    PS::F64 ecc_sq = sqrt(1.-ecc*ecc);
    pos = ax * ((cosu-ecc) * P + ecc_sq*sinu * Q);
    PS::F64 rinv = sqrt(1./(pos*pos));
    vel = ax*ax*n*rinv* (-sinu * P + ecc_sq*cosu * Q);
}
