#pragma once

class SolidDisk{
public:
    static PS::S32 n_init;
    static PS::F64 m_init;

    static PS::F64 p;
    //static PS::F64 f_in; 
    //static PS::F64 f_out;
    static PS::F64 f_dust;
    static PS::F64 eta_ice;
    static PS::F64 a_in;
    static PS::F64 a_out;
    static PS::F64 a_ice;

    static PS::F64 ecc_hill;
    static PS::F64 inc_hill;

    static PS::F64 calcDustMass(const PS::F64 a0,
                                const PS::F64 a1,
                                const bool inIce) {
        const PS::F64 L_CGS = 14959787070000;
        const PS::F64 M_CGS = 1.9884e33;
        
        if ( a1 < a0 ) return 0.;
        if ( inIce ) {
            const PS::F64 coef_in = 10. * f_dust /M_CGS*L_CGS*L_CGS;
            return 2.*M_PI*coef_in/(2.-p) * ( pow(a1, 2.-p) - pow(a0, 2.-p) );
        } else {
            const PS::F64 coef_out = 10. * f_dust * eta_ice /M_CGS*L_CGS*L_CGS;
            return 2.*M_PI*coef_out/(2.-p) * ( pow(a1, 2.-p) - pow(a0, 2.-p) );
        }
    }

    static PS::F64 getSemimajorAxis(const PS::F64 a0,
                                    const PS::F64 a1) {
        assert ( a0 < a1 );
        PS::F64 R = drand48();
        
        if ( p != 2 ){
            return pow( (pow(a1,2.-p) - pow(a0,2.-p)) * R + pow(a0,2.-p), 1./(2.-p) );
        } else {
            return exp( (log(a1) - log(a0)) * R + log(a0) );
        }
    }

    template <class Tpsys>
    static void createInitialCondition(Tpsys & pp){
        if ( PS::Comm::getRank() == 0 ){
            const PS::F64 m_sun = FPGrav::m_sun;
            PS::F64 m_in  = 0.;
            PS::F64 m_out = 0.;
            PS::S32 n_in  = 0;
            //PS::S32 n_out = 0;
        
            ////////////////////////////////////
            /*   Set Particle Mass & Number   */
            ////////////////////////////////////
            if ( a_out < a_ice ) {
                m_in = calcDustMass(a_in, a_out, true);
                m_out = 0.;
            } else if ( a_ice < a_in ) {
                m_in = 0.;
                m_out = calcDustMass(a_in, a_out, false);
            } else {
                m_in  = calcDustMass(a_in,  a_ice, true);
            m_out = calcDustMass(a_ice, a_out, false);
            }
            
            assert( n_init >= 0 );
            assert( m_init >= 0. );
            if ( m_init == 0. ) {
                assert( n_init > 0 );
                m_init = (m_in + m_out) / n_init;
            }
            if ( n_init == 0 ){
                assert( m_init > 0. );
                n_init = (m_in + m_out) / m_init;
            }
            
            n_in = (PS::S32)round(m_in/(m_in + m_out) * n_init);
            //n_out = n_init - n_in;
            
            ////////////////////////////////
            /*   Create Particle System   */
            ////////////////////////////////
            pp.setNumberOfParticleLocal(n_init);
            
            for ( PS::S32 i=0; i<n_init; i++ ){
                pp[i].id = i;
                pp[i].mass = m_init;
                
                // set orbital element
                PS::F64 ax;
                PS::F64 h = pow(pp[i].mass/(3.*m_sun), 1./3.);
                if ( a_out < a_ice || a_ice < a_in ) {
                    ax = getSemimajorAxis(a_in, a_out);
                } else {
                    if ( i < n_in ) {
                        ax = getSemimajorAxis(a_in, a_ice);
                    } else {
                        ax = getSemimajorAxis(a_ice, a_out);
                    }
                }
                PS::F64 ecc = getGaussian(ecc_hill*h);
                PS::F64 inc = getGaussian(inc_hill*h);
                
                PS::F64 l   = 2 * M_PI * drand48();
                PS::F64 u   = solveKeplerEq(l, ecc);
                PS::F64 omg = 2 * M_PI * drand48();
                PS::F64 OMG = 2 * M_PI * drand48();
                
                PS::F64 n = sqrt(m_sun / (ax*ax*ax));
                
                PS::F64vec P, Q;
                P.x =  cos(omg)*cos(OMG) - sin(omg)*sin(OMG)*cos(inc);
                P.y =  cos(omg)*sin(OMG) + sin(omg)*cos(OMG)*cos(inc);
                P.z =  sin(omg)*sin(inc);
                Q.x = -sin(omg)*cos(OMG) - cos(omg)*sin(OMG)*cos(inc);
                Q.y = -sin(omg)*sin(OMG) + cos(omg)*cos(OMG)*cos(inc);
                Q.z =  cos(omg)*sin(inc);
                
                orbitalElement2PosVel(pp[i].pos, pp[i].vel, m_sun,
                                      ax, ecc, n, u, P, Q);
            }
        } else {
            pp.setNumberOfParticleLocal(0);
        }
    }
};

PS::S32 SolidDisk::n_init = 0;
PS::F64 SolidDisk::m_init = 0.;
PS::F64 SolidDisk::p      = 1.5;
PS::F64 SolidDisk::f_dust  = 0.71;
PS::F64 SolidDisk::eta_ice = 30./7.1;
PS::F64 SolidDisk::a_in  = 0.98;
PS::F64 SolidDisk::a_out = 1.02;
PS::F64 SolidDisk::a_ice = 2.0;
PS::F64 SolidDisk::ecc_hill  = 2.0;
PS::F64 SolidDisk::inc_hill  = 1.0;


class GasDisk{
public:
    static PS::F64 alpha_gas;
    static PS::F64 beta_gas;
    static PS::F64 f_gas;
    static PS::F64 tau_gas;
    static PS::F64 C_d;
    static PS::F64 mu;

    PS::F64 coef_rho_gas;
    PS::F64 coef_cs_vk;
    PS::F64 coef_acc_gd;

    GasDisk(){
        const PS::F64 L_CGS = 14959787070000;
        const PS::F64 M_CGS = 1.9884e33;
        const PS::F64 T     = 365.25*24.*60.*60./(2.*M_PI);

        coef_rho_gas = 1.4e-9 * f_gas /M_CGS*L_CGS*L_CGS*L_CGS;

        const PS::F64 k_B = 1.380649e-16 /(M_CGS*L_CGS*L_CGS)*T*T;
        const PS::F64 N_A = 6.022140857e23;
        const PS::F64 m_H = 1./N_A /M_CGS;
        PS::F64 coef_cs = sqrt(k_B * 280 / (mu * m_H));
        PS::F64 coef_vk = sqrt(FPGrav::m_sun);
        coef_cs_vk = coef_cs / coef_vk;

        coef_acc_gd = 0.5*C_d*M_PI;

        if ( PS::Comm::getRank() == 0 ) {
            std::cout << "rho_gas at 1 AU = " << coef_rho_gas << std::endl
                      << "cs/vk at 1 AU   = " << coef_cs_vk << std::endl;
        }
    }

    template <class Tpsys>
    void calcGasDrag(Tpsys & pp,
                     PS::F64 time,
                     PS::F64 L=1.,
                     bool clear=true){

        const PS::S32 n_loc = pp.getNumberOfParticleLocal();

#pragma omp parallel for
        for(PS::S64 i=0; i<n_loc; i++){
            PS::F64 r2 = pp[i].pos.x*pp[i].pos.x + pp[i].pos.y*pp[i].pos.y;
            PS::F64 r_inv = 1./sqrt(r2);
            PS::F64 r = r2 * r_inv;
            PS::F64 rho_gas = coef_rho_gas * pow(r, -alpha_gas);
            if ( tau_gas != 0. ) rho_gas *= exp(-time / tau_gas);

            PS::F64 cs_vk = coef_cs_vk * sqrt(sqrt(r)) * pow(L, 1./8.);

            PS::F64vec ev(-pp[i].pos.y*r_inv, pp[i].pos.x*r_inv, 0.0);
            PS::F64vec vkep = sqrt(FPGrav::m_sun * r_inv) * ev;
            PS::F64 eta = 0.5 * (alpha_gas + beta_gas) * cs_vk * cs_vk;
            PS::F64vec vgas = (1.0 - eta)*vkep;
            PS::F64vec u = pp[i].vel - vgas;
            //PRL(eta);

            //PS::F64 rplanet = cbrt(0.75*pp[i].mass/(M_PI*FPGrav::dens));
            if (clear) pp[i].acc_gd = 0.;
            //pp[i].acc_gd += -coef_acc_gd * rplanet * rplanet * rho_gas * sqrt(u*u) * u / pp[i].mass;
            pp[i].acc_gd += -coef_acc_gd * pp[i].r_planet * pp[i].r_planet * rho_gas * sqrt(u*u) * u / pp[i].mass;
            pp[i].acc += pp[i].acc_gd;
        }
    }
};

PS::F64 GasDisk::alpha_gas = 11./4.;
PS::F64 GasDisk::beta_gas  = 0.5;
PS::F64 GasDisk::f_gas     = 0.71;
PS::F64 GasDisk::tau_gas   = 1.e6*2.*M_PI;
PS::F64 GasDisk::C_d       = 1.;
PS::F64 GasDisk::mu        = 2.34;

