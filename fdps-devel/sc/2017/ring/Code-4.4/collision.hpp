#pragma once

class Collision{
    static PS::F64 eps_tan;
    static PS::F64 eps_nor;
    static PS::F64 safty_factor_dr;
public:
    template<class Tpi, class Tpj>
    static void getNewPosAndVel(Tpi & pi, 
                                const Tpj & pj,
                                const PS::F64 _eps_tan,
                                const PS::F64 _eps_nor,
                                const PS::F64 _safty_factor_dr){
        const PS::F64vec rji = pj.pos  - pi.pos;
        const PS::F64vec vji = pj.vel  - pi.vel;
        const PS::F64 m_tot  = pi.mass + pj.mass;
        const PS::F64 r_sq   = rji * rji;
        const PS::F64 dis    = sqrt(rji*rji);
        PS::F64 dr     = 0.5*(2.0*Epj::r_coll - dis)*_safty_factor_dr;
        if(dr == 0.0) dr = Epj::r_coll*0.01;
        const PS::F64 C      = dr / dis;
        pi.pos -= C * rji;
        PS::F64vec nor = rji / dis;
        const PS::F64vec v_nor = (vji*nor)*nor;
        const PS::F64vec v_tan = vji - v_nor;
        pi.vel += pj.mass/m_tot*(1.0+_eps_nor)*v_nor;
        const PS::F64vec rji_tmp = (pj.pos+C*rji) - pi.pos;
    }
    template<class Tsys, class Ttree, class Tepj>
    static void correction(Tsys & system, 
                           const Ttree & tree,
                           const Tepj * epj_ngb,
                           const PS::F64 _eps_tan=eps_tan,
                           const PS::F64 _eps_nor=eps_nor,
                           const PS::F64 _safty_factor_dr=safty_factor_dr){
        const PS::S32 n = system.getNumberOfParticleLocal();
        for(PS::S32 i=0; i<n; i++){
            //std::cerr<<"system[i].n_coll="<<system[i].n_coll<<std::endl;
            if(system[i].n_coll > 0){
                //assert(system[i].n_coll < 2);
                //getNewPosAndVel(system[i], tree.getForce(i).ngb, _eps_tan, _eps_nor, _safty_factor_dr);
                const PS::S32 adr = tree.getForce(i).adr_ngb;
                getNewPosAndVel(system[i], epj_ngb[adr], _eps_tan, _eps_nor, _safty_factor_dr);
            }
        }
    }
};

PS::F64 Collision::eps_tan = 1.0;
PS::F64 Collision::eps_nor = 0.5;
PS::F64 Collision::safty_factor_dr = 1.0;
