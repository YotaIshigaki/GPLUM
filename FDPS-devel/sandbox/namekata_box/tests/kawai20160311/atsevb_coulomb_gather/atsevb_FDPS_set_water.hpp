//***************************************************************************************
//  This program is the intramolecular force calculation for "atsevb_main.cpp"
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include<cmath>
#include<stdexcept>

namespace SetWater {
    
    //--- coefficient table for intra force ---------------------------------------------
    //------ coefficient data class for intra force
    class CoefBond {
      public:
        PS::F64 a, d, r0, sign;
        
        template<class TcoefBond>
        inline void copyCoef(const TcoefBond & coef){
            this->a    = coef.a;
            this->d    = coef.d;
            this->r0   = coef.r0;
            this->sign = coef.sign;  // sign for force equation
        }
    };
    
    class CoefAngle {
      public:
        PS::F64 k, theta0;
        
        template<class TcoefAngle>
        inline void copyCoef(const TcoefAngle & coef){
            this->k      = coef.k;
            this->theta0 = coef.theta0;
        }
    };
    
    //------ grobal object of coefficient table
    MD_EXT::MultiTable<CoefBond,  ATOM_TYPE, ATOM_TYPE> coefTableBond;
    MD_EXT::MultiTable<CoefAngle, ATOM_TYPE, ATOM_TYPE, ATOM_TYPE> coefTableAngle;
    
    //------ initialize coefficient table
    void SetIntraParamWater(){
        //--- reserve table size
        PS::S32 key_size = static_cast<int>(ATOM_TYPE::ENUM_N_TOT);
        coefTableBond.reserve(key_size*key_size);
        coefTableAngle.reserve(key_size*key_size*key_size);
        
        //--- bond coefficient for water
        //------ H-O & O-H
        CoefBond tmp_bond;
        tmp_bond.a    = 2.287;
        tmp_bond.d    = 116.09;
        tmp_bond.r0   = 0.995;   // [angm]
        tmp_bond.sign = 1.0;
        coefTableBond.set(tmp_bond, ATOM_TYPE::H_water,
                                    ATOM_TYPE::O_water);  // H-O
        tmp_bond.sign = -1.0;
        coefTableBond.set(tmp_bond, ATOM_TYPE::O_water,
                                    ATOM_TYPE::H_water);  // O-H
                                   
        //--- angle coefficient for water
        //------ H-O-H & O<HH
        CoefAngle tmp_angle;
        tmp_angle.theta0 = 112.5*Unit::pi/180.0;        // [rad]
        tmp_angle.k      = 75.9/sin(tmp_angle.theta0);
        coefTableAngle.set(tmp_angle, ATOM_TYPE::H_water,
                                      ATOM_TYPE::O_water,
                                      ATOM_TYPE::H_water); // H-O-H
        coefTableAngle.set(tmp_angle, ATOM_TYPE::O_water,
                                      ATOM_TYPE::H_water,
                                      ATOM_TYPE::H_water); // O<HH
    }
    
    //------ setting data table for intermolecular force
    template<class T>
    void SetInterParamWater(T & coefTable_inter){
        InterForceParam coef_inter;
        //------ for H
        coef_inter.mass   = 1.0079e-3/Unit::mass_C;
        coef_inter.charge = 0.4175*Unit::coef_coulomb;
        coef_inter.VDW_R  = 0.0;
        coef_inter.VDW_D  = 0.0;
        coefTable_inter.set(coef_inter, ATOM_TYPE::H_water);
        //------ for O
        coef_inter.mass   = 16.000e-3/Unit::mass_C;
        coef_inter.charge = -0.8350*Unit::coef_coulomb;
        coef_inter.VDW_R  =  3.553145;
        coef_inter.VDW_D  =  0.1554253;
        coefTable_inter.set(coef_inter, ATOM_TYPE::O_water);
    }
    
    //--- check collistion or not between new position and all particles in "system".
    template<class Tpsys>
    bool isCollisionWater(const Tpsys & system,
                          const PS::S32 & n_increment,
                          const PS::F64vec & pos_new){
        
        constexpr PS::F64 range = 1.3*1.3;  // radius [angm^2]
        
        for(PS::S32 i=0; i<n_increment; i++){
            PS::F64vec r_ij = pos_new - Normalize::realPos( system[i].getPos() );
            if( range > r_ij*r_ij ) return true;
        }
        
        return false;
    }
    
    //--- installing water molecule into the system
    template<class Tpsys, class Ttable,typename Trand>
    void InstallWater(Tpsys & system,
                      PS::S32 & n_increment,
                      PS::S32 & mol_increment,
                      const Ttable & coefTable_inter,
                      const PS::S32 & n_add_mol,
                      Trand & mt){
        
        //--- proceed in PS::Comm::getRank() == 0
        if( PS::Comm::getRank() != 0 ) return;
        
        //--- message
        std::cerr << "  installing water:" << n_add_mol << " molecules." << std::endl;
        
        //--- data length check
        PS::S32 n_add_atom = n_add_mol*3;
        if( system.getNumberOfParticleLocal() < (n_increment + n_add_atom) ){
            throw std::length_error("ERROR: n_loc of system is too small.");
        }
        
        //--- install molecules
        for(PS::S32 add=0; add<n_add_mol; add++){
            PS::S32 count = 0;
            while(true){
                //--- finit loop condition
                count++;
                if(count > 100000){
                    throw std::domain_error("ERROR: cannot find installing position");
                }
                
                //--- searching install position
                PS::F64vec pos_O;
                pos_O[0] = mt.genrand_res53();
                pos_O[1] = mt.genrand_res53();
                pos_O[2] = mt.genrand_res53();
                pos_O = Normalize::realPos(pos_O);
                if( isCollisionWater(system, n_increment, pos_O) ) continue;
                
                //--- installing
                PS::S32 mol_id = mol_increment;
                //------ O atom of water
                PS::S32 O_id = n_increment;
                InterForceParam coef_inter = coefTable_inter.get( ATOM_TYPE::O_water );
                system[n_increment].setPos( Normalize::normPos(pos_O) );
                system[n_increment].setAtomID( O_id );
                system[n_increment].setMolID( mol_id );
                system[n_increment].setMolType( MOL_TYPE::solvent );
                system[n_increment].setAtomType( ATOM_TYPE::O_water );
                system[n_increment].setVel( PS::F64vec(0.0, 0.0, 0.0) );
                system[n_increment].setMass(   coef_inter.mass );
                system[n_increment].setCharge( coef_inter.charge );
                system[n_increment].setVDW_R(  coef_inter.VDW_R );
                system[n_increment].setVDW_D(  coef_inter.VDW_D );
                
            //    std::cerr << " installing: O  pos={" << system[n_increment].getPos() << "}" << std::endl
            //              << "           real_pos={" << Normalize::realPos(system[n_increment].getPos())
            //              << "}" << std::endl;
                
                n_increment++;
                //------ H atom of water (1)
                PS::S32 H_id_1 = n_increment;
                                coef_inter = coefTable_inter.get( ATOM_TYPE::H_water );
                CoefBond coef_bond = coefTableBond.get( ATOM_TYPE::H_water,
                                                        ATOM_TYPE::O_water );
                PS::F64vec pos_shift = PS::F64vec(coef_bond.r0, 0.0, 0.0);
                PS::F64vec pos_H_1 = Normalize::normPos(pos_O + pos_shift);
                           pos_H_1 = Normalize::periodicAdjust(pos_H_1);
                system[n_increment].setPos( pos_H_1 );
                system[n_increment].setAtomID( H_id_1 );
                system[n_increment].setMolID( mol_id );
                system[n_increment].setMolType( MOL_TYPE::solvent );
                system[n_increment].setAtomType( ATOM_TYPE::H_water );
                system[n_increment].setVel( PS::F64vec(0.0, 0.0, 0.0) );
                system[n_increment].setMass(   coef_inter.mass );
                system[n_increment].setCharge( coef_inter.charge );
                system[n_increment].setVDW_R(  coef_inter.VDW_R );
                system[n_increment].setVDW_D(  coef_inter.VDW_D );
                
            //    std::cerr << " installing: H1 pos={" << system[n_increment].getPos() << "}" << std::endl
            //              << "           real_pos={" << Normalize::realPos(system[n_increment].getPos())
            //              << "}" << std::endl;
                
                n_increment++;
                //------ H atom of water (2)
                PS::S32 H_id_2 = n_increment;
                CoefAngle coef_angle = coefTableAngle.get( ATOM_TYPE::H_water,
                                                           ATOM_TYPE::O_water,
                                                           ATOM_TYPE::H_water );
                PS::F64vec pos_H_2;
                           pos_H_2[0] = pos_shift.x*cos( coef_angle.theta0 );
                           pos_H_2[1] = pos_shift.x*sin( coef_angle.theta0 );
                           pos_H_2 = Normalize::normPos(pos_O + pos_H_2);
                           pos_H_2 = Normalize::periodicAdjust(pos_H_2);
                system[n_increment].setPos( pos_H_2 );
                system[n_increment].setAtomID( H_id_2 );
                system[n_increment].setMolID( mol_id );
                system[n_increment].setMolType( MOL_TYPE::solvent );
                system[n_increment].setAtomType( ATOM_TYPE::H_water );
                system[n_increment].setVel( PS::F64vec(0.0, 0.0, 0.0) );
                system[n_increment].setMass(   coef_inter.mass );
                system[n_increment].setCharge( coef_inter.charge );
                system[n_increment].setVDW_R(  coef_inter.VDW_R );
                system[n_increment].setVDW_D(  coef_inter.VDW_D );
                
            //    std::cerr << " installing: H2 pos={" << system[n_increment].getPos() << "}" << std::endl
            //              << "           real_pos={" << Normalize::realPos(system[n_increment].getPos())
            //              << "}" << std::endl;
                
                n_increment++;
                
                //--- connection setting
                //------ O atom of water
                system[O_id].bond.clear();
                system[O_id].bond.add( H_id_1 );
                system[O_id].bond.add( H_id_2 );
                //------ H atom of water (1)
                system[H_id_1].bond.clear();
                system[H_id_1].bond.add( O_id );
                //------ H atom of water (2)
                system[H_id_2].bond.clear();
                system[H_id_2].bond.add( O_id );
                
                
                //--- increment molID
                mol_increment++;
                break;
            }
        }
    }
    
    //--- calculate intramolecular potential --------------------------------------------
    //------ bond potential P(0)-P(1)
    inline void CalcBondForce(const IntraForceTarget<3> &target,
                              ForceIntra* force){
        const CoefBond coef = coefTableBond.get(target.getAtomType(0),
                                                target.getAtomType(1));
        
        //--- potential
        constexpr PS::F64 factor_pot = 7.0/(12.0*2.0);  // counteract double count
        PS::F64 ar  = coef.a*(target.getR(0) - coef.r0);
        PS::F64 ar2 = ar*ar;
        force[target.getIndex(0)].addPotIntra( coef.d*(ar2*(1.0-ar) + (factor_pot*ar2*ar2)) );
        
        //--- force
        constexpr PS::F64 factor_force = 7.0/3.0;
        PS::F64    ebp   = coef.d*coef.a*( (2.0*ar) - (3.0*ar2) + factor_force*ar*ar2 );
        PS::F64vec f_tmp = ebp*target.getR(0)*target.getRij(0);
        force[target.getIndex(0)].addForceIntra( coef.sign*f_tmp );
        
        //--- virial
        f_tmp[0] = f_tmp[0]*target.getRij(0)[0];
        f_tmp[1] = f_tmp[1]*target.getRij(0)[1];
        f_tmp[2] = f_tmp[2]*target.getRij(0)[2];
        force[target.getIndex(0)].addVirtIntra( -0.5*f_tmp );
    }
    //------ angle potential P(k=2)-P(i=0)-P(j=1)
    inline void CalcAngleForce_center(const IntraForceTarget<3> &target,
                                      ForceIntra* force){
        
        const CoefAngle  coef = coefTableAngle.get(target.getAtomType(0),
                                                   target.getAtomType(1),
                                                   target.getAtomType(2));
        const PS::F64vec R_ij = target.getRij(0);
        const PS::F64    rij  = sqrt(R_ij*R_ij);
        const PS::F64vec R_ik = target.getRij(1) + target.getRij(0);
        const PS::F64    rik  = sqrt(R_ik*R_ik);
        constexpr PS::F64 factor = 0.5/3.0;  // counteract triple count
        
        PS::F64 proi = R_ij*R_ik;
        PS::F64 cost = proi/( rij*rik );
        PS::F64 cost_diff = cost - cos(coef.theta0);
        
        //--- potential
        force[target.getIndex(0)].addPotIntra( factor*coef.k*cost_diff*cost_diff );
        
        //--- force
        PS::F64    r_inv = 1.0/( (rij*rij)*(rik*rik) );
        PS::F64vec rtri  = ( (rij*rij*R_ik + rik*rik*R_ij)*cost - rij*rik*(R_ij + R_ik) )*r_inv;
        force[target.getIndex(0)].addForceIntra( -coef.k*cost_diff*rtri );
        
        //--- virial
        // added in CalcAngleForce_linear().
    }
    //------ angle potential P(i=0)-P(j=1)-P(k=2)
    inline void CalcAngleForce_linear(const IntraForceTarget<3> &target,
                                      ForceIntra* force){
        
        const CoefAngle  coef = coefTableAngle.get(target.getAtomType(0),
                                                   target.getAtomType(1),
                                                   target.getAtomType(2));
        const PS::F64vec R_ij = target.getRij(0);
        const PS::F64    rij  = sqrt(R_ij*R_ij);
        const PS::F64vec R_jk = target.getRij(1);
        const PS::F64    rjk  = sqrt(R_jk*R_jk);
        constexpr PS::F64 factor = 0.5/3.0;  // counteract triple count
        
        PS::F64 proi = R_ij*R_jk;
        PS::F64 cost = proi/( rij*rjk );
        PS::F64 cost_diff = cost - cos(coef.theta0);
        
        //--- potential
        force[target.getIndex(0)].addPotIntra( factor*coef.k*cost_diff*cost_diff );
        
        //--- force
        PS::F64    r_inv = 1.0/( rij*rij*rjk );
        PS::F64vec rtri  = (rij*R_jk - rjk*R_ij*cost)*r_inv;
        force[target.getIndex(0)].addForceIntra( -coef.k*cost_diff*rtri );
        
        //--- virial
        rtri[0] = rtri[0]*R_ij[0];
        rtri[1] = rtri[1]*R_ij[1];
        rtri[2] = rtri[2]*R_ij[2];
        force[target.getIndex(0)].addVirtIntra( -coef.k*cost_diff*rtri );
    }
    
    //------ search intramolecular pair at root node
    template<class Tepi, class Tepj, class Ttable>
    void CalcForceIntra(const Tepi* EPI,
                        const PS::S32 root_i,
                        const Tepj* EPJ,
                        const PS::S32 n_EPJ,
                        const Ttable& IDtable,
                        ForceIntra* force){
        
        //--- data buffer
        IntraForceTarget<3> target;  // target atom list length = 3
        target.clear();
        target.add(EPI, root_i);  // input root node
        
        //--- bond potential  Pi-P(j)
        for(PS::S32 j=0; j<EPI[root_i].bond.getN(); j++){
            target.add(EPJ, IDtable.at( EPI[root_i].bond.get(j) ) );  // input j node
            CalcBondForce(target, force);
            target.pop();  // delete j node
        }
        //--- angle potential
        for(PS::S32 j=0; j<EPI[root_i].bond.getN(); j++){
            PS::S32 index_j = IDtable.at(EPI[root_i].bond.get(j) );
            target.add(EPJ, index_j);  // input j node
            //--- linear shape. Pi-P(j)-P(k)
            for(PS::S32 k=0; k<EPJ[index_j].bond.getN(); k++){
                if( target.checkRootNode( EPJ[index_j].bond.get(k) ) ) continue;  // skip Pi == P(k)
                PS::S32 index_k = IDtable.at( EPJ[index_j].bond.get(k) );
            
                target.add(EPJ, index_k);  // add k node
                CalcAngleForce_linear(target, force);
                target.pop();  //delete k node
            }
            //--- center shape. P(k)-Pi-P(j)
            for(PS::S32 k=j+1; k<EPI[root_i].bond.getN(); k++){
                
                target.add(EPJ, IDtable.at( EPI[root_i].bond.get(k) ) );  // input k node
                CalcAngleForce_center(target, force);
                target.pop();  // delete k node
            }
            target.pop();  // delete j node
        }
    }
}

