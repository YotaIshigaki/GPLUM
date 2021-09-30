//***************************************************************************************
//  This is basic kick and drift routine for atom.
//***************************************************************************************
#pragma once

namespace ATOM_MOVE {

    template<class Tpsys>
    void kick(const PS::F64 &dt,
                    Tpsys   &psys){

        //--- kick atom
        const PS::S64 n = psys.getNumberOfParticleLocal();

        PS::F64vec v_barycentric = 0.0;
        PS::F64    mass_total    = 0.0;
        for(PS::S64 i=0; i<n; ++i){
            PS::F64vec v_new = psys[i].getVel()
                             + psys[i].getForce()*( dt/psys[i].getMass() );
            psys[i].setVel(v_new);

            v_barycentric += psys[i].getMass()*v_new;
            mass_total    += psys[i].getMass();
        }

        //--- cancel barycentric velocity
        mass_total    = PS::Comm::getSum(mass_total);
        v_barycentric = PS::Comm::getSum(v_barycentric);
        v_barycentric = v_barycentric*(1.0/mass_total);
        for(PS::S64 i=0; i<n; ++i){
            psys[i].setVel( psys[i].getVel() - v_barycentric );
        }
    }

    template<class Tpsys>
    PS::F64 drift(const PS::F64 &dt,
                        Tpsys   &psys){

        PS::F64 max_move = 0.0;
        const PS::S64 n  = psys.getNumberOfParticleLocal();
        for(PS::S64 i=0; i<n; ++i){
            PS::F64vec move    = psys[i].getVel()*dt;
            PS::F64vec pos_new = psys[i].getPos() + Normalize::normDrift(move);
            psys[i].addTrj( move );
            psys[i].setPos( Normalize::periodicAdjustNorm(pos_new) );

            //--- check largest move at step
            move     = Normalize::normDrift(move);
            max_move = std::max(max_move, move*move);
        }
        max_move = std::sqrt(max_move);

        if(max_move > 0.5){
            std::cerr << "  max_move = " << max_move << std::endl;
            throw std::logic_error("atoms speed runaway");
        }
        return max_move;
    }

}
