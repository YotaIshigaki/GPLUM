#include<iostream>
#include<fstream>
#include<unistd.h>
#include<particle_simulator.hpp>
#include"kepler.hpp"

int main(int argc, char *argv[]){
    static const PS::F64 PI = 4.0 * atan(1.0);
    const PS::F64 ax0 = 1.0;
    const PS::F64 ecc0 = 0.7;
    const PS::F64 inc0 = 0.1 * 2.0 * PI;
    const PS::F64 OMG0 = 0.2 * 2.0 * PI;
    const PS::F64 omg0 = 0.3 * 2.0 * PI;
    const PS::F64 u0 = 0.4 * 2.0 * PI;
    const PS::F64 mass0 = 1.0;
    const PS::F64 mass1 = 0.0;
    const PS::F64 l0 = u0 - ecc0*sin(u0);
    
    PS::F64vec pos0_old, pos1_old, vel0_old, vel1_old;
    OrbParam2PosVel(pos0_old, pos1_old, vel0_old, vel1_old, mass0, mass1,
		    ax0, ecc0, inc0, OMG0, omg0, u0);

    PS::F64 ax1, ecc1, inc1, OMG1, omg1;
    PosVel2OrbParam(ax1, ecc1, inc1, OMG1, omg1,
		    pos0_old, pos1_old, vel0_old, vel1_old, mass0, mass1);

    const PS::F64 n = sqrt( (mass0+mass1) / (ax1*ax1*ax1) );
    const PS::F64 Tkep = 2.0 * PI / n;
    //const PS::F64 dt = Tkep * 0.5;
    const PS::F64 dt = Tkep;
    const PS::F64 l1 = n * dt + l0;
    PS::F64 u1 = solve_keplereq(l1, ecc1);
    const PS::F64 l2 = u1 - ecc1*sin(u1);
    std::cout<<"ax1= "<<ax1<<" ax0= "<<ax0<<std::endl;
    std::cout<<"ecc1= "<<ecc1<<" ecc0= "<<ecc0<<std::endl;
    std::cout<<"inc1= "<<inc1<<" inc0= "<<inc0<<std::endl;
    std::cout<<"OMG1= "<<OMG1<<" OMG0= "<<OMG0<<std::endl;
    std::cout<<"omg1= "<<omg1<<" omg0= "<<omg0<<std::endl;
    std::cout<<"l2= "<<l2<<" l1= "<<l1<<" l0= "<<l0<<std::endl;
    std::cout<<"u1= "<<u1<<" u0= "<<u0<<std::endl;

    PS::F64vec pos0_new, pos1_new, vel0_new, vel1_new;
    OrbParam2PosVel(pos0_new, pos1_new, vel0_new, vel1_new, mass0, mass1,
		    ax1, ecc1, inc1, OMG1, omg1, u1);
	
    std::cout<<"pos0_old="<<pos0_old<<" pos0_new="<<pos0_new<<std::endl;
    std::cout<<"pos1_old="<<pos1_old<<" pos1_new="<<pos1_new<<std::endl;
    std::cout<<"vel0_old="<<vel0_old<<" vel0_new="<<vel0_new<<std::endl;
    std::cout<<"vel1_old="<<vel1_old<<" vel1_new="<<vel1_new<<std::endl;
    
    return 0;
}
