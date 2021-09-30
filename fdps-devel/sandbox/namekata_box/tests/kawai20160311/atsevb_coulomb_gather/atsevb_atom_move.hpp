//***************************************************************************************
//  This is kick and drift routine.
//    This code is used by "atsevb_main.cpp"
//***************************************************************************************
#pragma once

template<class Tpsys>
void Kick(Tpsys & system,
          const PS::F64 dt){
    PS::S32 n = system.getNumberOfParticleLocal();
    for(int i=0; i<n; i++){
    //    std::cout << "ID:" << system[i].getAtomID() << endl
    //              << " pos  ={" << system[i].getPos() << "}" << std::endl
    //              << " force={" << system[i].getForce() << "}" << std::endl << std::flush;
        system[i].accel(system[i].getForce(), dt);
    //    std::cout << "ID:" << system[i].getAtomID()
    //              << " vel={" << system[i].getVel() << "}" << std::endl;
    
    
        //--- test const force
      //  PS::F64vec f = PS::F64vec(1.0, 0.0, 0.0);
      //  if( (i % 2) == 0 ) f = -f;
      //  system[i].accel(f, dt);
    }
}

template<class Tpsys>
void Drift(Tpsys & system,
           const PS::F64 dt){
    PS::S32 n = system.getNumberOfParticleLocal();
    
    //return;  // skip move
    for(int i=0; i<n; i++){
        PS::F64vec pos_new = system[i].getPos() + Normalize::normDrift( system[i].getVel()*dt );
    //    std::cerr << "ID:" << system[i].getAtomID()
    //              << " pos_new={" << Normalize::realPos(pos_new) << "}" << std::endl;
        //std::cout << "[B] pos_new.x = " << pos_new.x 
        //          << "(" << PS::Comm::getRank() << ")" << std::endl;
        //std::cout << "[B] pos_new.x = " << pos_new.x 
        //          << "(" << PS::Comm::getRank() << ")" << std::endl;
        //std::cout << "[B] pos_new.x = " << pos_new.x 
        //          << "(" << PS::Comm::getRank() << ")" << std::endl;
        system[i].setPos(pos_new);
        //system[i].setPos( Normalize::periodicAdjust(pos_new) );

        //pos_new = system[i].getPos();
        //std::cout << "[A] pos_new.x = " << pos_new.x 
        //          << "(" << PS::Comm::getRank() << ")" << std::endl;
        //std::cout << "[A] pos_new.x = " << pos_new.x 
        //          << "(" << PS::Comm::getRank() << ")" << std::endl;
        //std::cout << "[A] pos_new.x = " << pos_new.x 
        //          << "(" << PS::Comm::getRank() << ")" << std::endl;
        
    //    std::cout << "  pos:" << system[i].getPos() << std::endl;
        
        //--- test move
    //    std::cout << "ID:" << system[i].getAtomID() << endl
    //              << " pos  ={" << system[i].getPos() << "}" << std::endl
    //              << " force={" << system[i].getForce() << "}" << std::endl << std::flush;
    //    PS::F64vec pos_new = Normalize::normDrift( PS::F64vec(0.5, 0.0, 0.0) );
    //    if( (system[i].getAtomID() % 2) == 0 ){
    //        pos_new = system[i].getPos() + pos_new;
    //    } else {
    //        pos_new = system[i].getPos() - pos_new;
    //    }
    //    system[i].setPos( Normalize::periodicAdjust(pos_new) );
    }
}

