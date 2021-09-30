//***************************************************************************************
//  This routine is initializer for "atsevb_main.cpp".
//***************************************************************************************
#pragma once

class SystemSetting {
  private:
    //--- singleton rule
    SystemSetting() = default;
    ~SystemSetting() = default;
    
    
    //--- for time
    PS::S64 istep;
    PS::S64 nstep_st;
    PS::S64 nstep_ed;
    PS::F64 dt;
    
    //--- for Tree
    const PS::S32 n_leaf_limit = 8;
    const PS::F64 coef_ema = 0.3;
    PS::F64 theta;
    PS::S32 n_group_limit;
    PS::S32 cycle_dinfo;
    
  public:
    //--- unable object copy
    SystemSetting(const SystemSetting&) = delete;
    SystemSetting& operator=(const SystemSetting&) = delete;
    SystemSetting(SystemSetting&&) = delete;
    SystemSetting& operator=(SystemSetting&&) = delete;
        
    //--- call singleton object
    static SystemSetting* getInstance(void) {
        static SystemSetting singleton;
        return &singleton;
    }
  
    void setIstep(const PS::S64 i){       this->istep = i; }
    void setNstep_st(const PS::S64 i){    this->nstep_st = i; }
    void setNstep_ed(const PS::S64 i){    this->nstep_ed = i; }
    void set_dt(const PS::F64 t){         this->dt = t/Unit::norm_time; }
    void setTheta(const PS::F64 f){       this->theta = f; }
    void setNGroupLimit(const PS::S32 n){ this->n_group_limit = n; }
    void setCycleDinfo(const PS::S32 i){  this->cycle_dinfo = i; }
    
    PS::S64 getIstep() const {       return this->istep; }
    PS::S64 getNstep_st() const {    return this->nstep_st; }
    PS::S64 getNstep_ed() const {    return this->nstep_ed; }
    PS::F64 get_dt() const {         return this->dt; }
    PS::F64 getTheta() const {       return this->theta; }
    PS::S32 getNGroupLimit() const { return this->n_group_limit; }
    PS::S32 getCycleDinfo() const {  return this->cycle_dinfo; }
    
    //--- const value
    PS::S32 getNLeafLimit() const {  return this->n_leaf_limit; }
    PS::F64 getCoefEma() const {     return this->coef_ema; }
    
    //--- next step
    void StepNext(){
        this->istep++;
    }
    
    //--- dinfo timing
    bool isDinfoUpdate(){
        return ( ((this->istep - this->nstep_st) % this->cycle_dinfo) == 0 );
    }
    
    //--- judge continue lopp or not
    bool isLoopContinue(){
        return ( this->istep <= this->nstep_ed );
    }
};

template<class Tdinfo>
void Init_Condition(Tdinfo & dinfo){
    
    //--- get instance of SystemSetting
    SystemSetting* setting = SystemSetting::getInstance();
    
    //--- input settings
    //------ time
    setting->setNstep_st(0);    // test conditions
    setting->setIstep(setting->getNstep_st());
    setting->setNstep_ed(10000);
    //setting->set_dt(1.0e-15);   // [s]
    setting->set_dt(1.0e-15);   // [s]
    //------ tree parameters
    setting->setTheta(0.5);
    setting->setNGroupLimit(64);
    //------ update timing of domain info
    setting->setCycleDinfo(20);
    
    //------ file I/O timing
    FILE_IO::cycle_eng_log = 1;
    FILE_IO::cycle_pos     = 10;
    FILE_IO::cycle_resume  = 10000;
    FILE_IO::cycle_vmd     = 1000;
    
    //--- set dinfo
    dinfo.initialize(setting->getCoefEma());
    //------ set boundary condition
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    //------ set system box size
    Normalize::setBoxSize( PS::F64vec(100.0, 100.0, 100.0) ); // real space size
    dinfo.setPosRootDomain(PS::F64vec( 0.0, 0.0, 0.0),
                           PS::F64vec( 1.0, 1.0, 1.0) );   // normalized space (fixed value for PM)
    
    //--- initialize cutoff radius
    EPI_LJ::rSearch      = Normalize::normCutOff(12.0);  // 12 angm. normalized space
    EPI_intra::rSearch   = Normalize::normCutOff(12.0);  //  6 angm. normalized space
    EPI_coulomb::rSearch = 3.0/SIZE_OF_MESH;             //  3 mesh. normalized space (fixed)
    EPJ_LJ::rSearch      = EPI_LJ::rSearch;     
    EPJ_intra::rSearch   = EPI_intra::rSearch;  
    EPJ_coulomb::rSearch = EPI_coulomb::rSearch;
    
    if(PS::Comm::getRank() == 0){
        std::cerr << "  cutoff condition" << std::endl
                  << " LJ     : " << Normalize::realCutOff(EPI_LJ::rSearch) << std::endl
                  << " intra  : " << Normalize::realCutOff(EPI_intra::rSearch) << std::endl
                  << " coulomb: " << Normalize::realCutOff(EPI_coulomb::rSearch) << std::endl;
    }
}

template<class Tpsys>
void Init_Molecule(Tpsys & system){
    
    //--- number of molecule and atom
    PS::S32 n_wat    = 1;  // number of water molecule
    PS::S32 atom_wat = 3;  // atoms in water molecule
    
    //--- setting name tag
    ENUM_TAG::Init();
    
    //--- setting file I/O
    FILE_IO::Init();
    
    //--- initialize randum number generator
    PS::MTTS mt;
    mt.init_genrand(PS::Comm::getRank()*PS::Comm::getNumberOfThread()+PS::Comm::getThreadNum());
  
    //--- setting forcefield parameter
    SetWater::SetIntraParamWater();
    SetWater::SetInterParamWater( Grobal::coefTable_inter );
    
    
    //--- below routine is performed in process 0 only
    if( PS::Comm::getRank() != 0 ) return;
    
  //  //--- setting the total number of atom
  //  PS::S32 n_atom_tot = n_wat*atom_wat;
  //  system.setNumberOfParticleLocal( n_atom_tot );
  //  
  //  
  //  //--- install molecules
  //  PS::S32 n_increment   = 0;
  //  PS::S32 mol_increment = 0;
  //  //------ for polymer
  //  //------ for water molecule
  //  SetWater::InstallWater(system,
  //                         n_increment,
  //                         mol_increment,
  //                         Grobal::coefTable_inter,
  //                         n_wat,
  //                         mt);
  
    //--- test condition --- charge + VDW 2 body ------------------------------
    system.setNumberOfParticleLocal(2);
    
    system[0].setPos( Normalize::normPos( PS::F64vec(40.0, 10.0, 10.0) ) );
    system[0].setAtomID( 0 );
    system[0].setMolID( 0 );
    system[0].setMolType( MOL_TYPE::polymer );
    system[0].setAtomType( ATOM_TYPE::C_sc );
    system[0].setVel( PS::F64vec(0.0, 0.0, 0.0) );
    system[0].setMass( 0.1 );
    system[0].setCharge( +1.0 );
    system[0].setVDW_R(  3.5 );
    system[0].setVDW_D(  0.05 );
    
    
    system[1].setPos( Normalize::normPos( PS::F64vec(42.5, 10.0, 10.0) ) );
    system[1].setAtomID( 1 );
    system[1].setMolID( 1 );
    system[1].setMolType( MOL_TYPE::polymer );
    system[1].setAtomType( ATOM_TYPE::C_sc );
    system[1].setVel( PS::F64vec(0.0, 0.0, 0.0) );
    system[1].setMass( 0.1 );
    system[1].setCharge( -1.0 );
    system[1].setVDW_R(  3.5 );
    system[1].setVDW_D(  0.05 );
    //-------------------------------------------------------------------------
    
    //--- test condition --- VDW & charge nbody -------------------------------
  //  const PS::S32 n_test  = 20; // n
  //  const PS::F64 q_value = 1.0; // charge strength
  //  const PS::F64 T_real  = 600; // [K]
  //  const PS::F64 mass    = 1.0;
  //  PS::F64 r2_dist = Normalize::normXLen(10.0);
  //          r2_dist = r2_dist*r2_dist;
  //  system.setNumberOfParticleLocal(n_test);
  //  
  //  for(PS::S32 i=0; i<n_test; i++){
  //      
  //      PS::F64vec pos_tmp;
  //      PS::S32 pos_count = 0;
  //      while(true){
  //          pos_tmp[0] = mt.genrand_res53();  // [0,1) random number
  //          pos_tmp[1] = mt.genrand_res53();
  //          pos_tmp[2] = mt.genrand_res53();
  //                     
  //          bool flag = true;
  //          for(PS::S32 j=0; j<i; j++){
  //              PS::F64vec r_ij = system[i].getPos() - pos_tmp;
  //              PS::F64    r2   = r_ij*r_ij;
  //              if(r2_dist > r2) flag = false;
  //          }
  //          
  //          if(flag) break;       
  //          pos_count++;
  //          
  //          if(pos_count > 1000000){
  //              throw std::domain_error("ERROR: cannot find installing position");
  //          }
  //      }
  //      
  //      PS::F64vec vel_tmp;
  //      PS::F64    rdev    = sqrt( 2.0*T_real/(Unit::norm_temp*(mass/Unit::mass_C)) );
  //      
  //      PS::F64    random1 = sqrt( -log(mt.genrand_res53()) );
  //      PS::F64    random2 = 2.0*(Unit::pi)*(mt.genrand_res53());
  //                 vel_tmp[0] = rdev*random1*sin(random2);
  //                 
  //                 random1 = sqrt( -log(mt.genrand_res53()) );
  //                 random2 = 2.0*(Unit::pi)*(mt.genrand_res53());
  //                 vel_tmp[1] = rdev*random1*sin(random2);
  //                 
  //                 random1 = sqrt( -log(mt.genrand_res53()) );
  //                 random2 = 2.0*(Unit::pi)*(mt.genrand_res53());
  //                 vel_tmp[2] = rdev*random1*sin(random2);
  //                 
  //              //   vel_tmp = Normalize::normPos( vel_tmp );
  //                 
  //      PS::F64   charge_tmp;
  //      ATOM_TYPE type_tmp;
  //      if( (i % 2) == 0){
  //          charge_tmp = q_value*(Unit::coef_coulomb);
  //          type_tmp   = ATOM_TYPE::H_water;
  //      } else {
  //          charge_tmp = -q_value*(Unit::coef_coulomb);
  //          type_tmp   = ATOM_TYPE::O_water;
  //      }
  //      
  //      system[i].setPos( pos_tmp );
  //      system[i].setAtomID( i );
  //      system[i].setMolID( i );
  //      system[i].setMolType( MOL_TYPE::test );
  //      system[i].setAtomType( type_tmp );
  //      system[i].setVel( vel_tmp );
  //      system[i].setMass( mass/Unit::mass_C );
  //      system[i].setCharge( charge_tmp );
  //      system[i].setVDW_R(  3.5 );
  //      system[i].setVDW_D(  0.15 );
  //  }
    //-------------------------------------------------------------------------
    
    //--- initialzie intra list
    PS::S32 order = 3;
    MD_EXT::Init_IntraList<PS::S32>(system, order);
}

