/*************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
particle class for 3D!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
**************************************************/
#pragma once

class FileHeader{
public:
  PS::F64 time;
  PS::S32 Nbody;

  void readBinary(FILE* fp) {
    fread( &this->Nbody, sizeof(PS::S32), 1, fp );
  }
  void writeAscii(FILE* fp) const{
    fprintf(fp,"%e\n",time);
    fprintf(fp,"%d\n",Nbody);
    fprintf( fp, "itype,id,istate,dens,x,y,z,velx,vely,velz,sigx,y,z,xy,xz,yz,epsEx,y,z,xy,xz,yz,epsPx,y,z,xy,xz,yz,omg.xy,xz,yz,acc.x,y,z,one\n");
  }
  void writeBinary(FILE* fp) const{
    fwrite( &this->time , sizeof(PS::F64), 1, fp );
    fwrite( &this->Nbody, sizeof(PS::S32), 1, fp );
  }
};

class Hydro{
public:
  PS::F64    dens_rate;
  PS::F64vec acc;
  PS::F64mat eps_rate;
  PS::F64vec omg_rate;
  void clear(){
    dens_rate = 0.0;
    acc       = 0.0;
    eps_rate  = 0.0;
    omg_rate  = 0.0;
  }
};

struct FP{
  PS::S64    id;     //const particle number
  PS::U32    bdry;
  PS::U32    itype;  //const particle type
  PS::U32    istate; //2:elastic,3:plastic
  PS::U32    prank;
  // PS::S32 count;
  PS::F64    mass;   //const
  PS::F64    dens;
  PS::F64vec pos;
  PS::F64vec vel;
  PS::F64vec acc;    //force
  PS::F64mat sig;
  PS::F64mat eps;
  PS::F64mat epsp;   //plastic strain
  PS::F64mat epse;   //elastic strain

  PS::F64    dens_rate;
  PS::F64mat sig_rate;
  PS::F64mat eps_rate;
  PS::F64vec omg_rate;
  /*** ^antisymmetric tensor!!! = omg_rate.XY = -omg_rate.YX ***/
  //leap-frog
  PS::F64    dens_half;
  PS::F64vec vel_half;
  PS::F64mat sig_half;
  // phisical property
  PS::F64    al_phi;
  PS::F64    kc;
  PS::F64    coh;
  PS::F64    Ey;     // Young's modulus
  // get Prank
  // Eigen Solver
  PS::F64mat R_tensor;
  // displacement
  PS::F64vec disp;

  PS::F64vec getPos() const{
    return this->pos;//FP->FDPS
  }
  void copyFromForce(const Hydro & force){
    this->dens_rate = force.dens_rate;
    this->acc       = force.acc;
    this->eps_rate  = force.eps_rate;
    this->omg_rate  = force.omg_rate;
  }
  void setPos(const PS::F64vec& pos){
    this->pos       = pos;
  }
  PS::F64 getRSearch() const{
    return kappa * smth_len;// 2.0*1.2*d
  }

  void readBinary(FILE* fp) {
    fread(this, sizeof(FP), 1, fp);
      // fread( &this->bdry      , sizeof(PS::U32)   , 1, fp );
      // fread( &this->itype     , sizeof(PS::U32)   , 1, fp );
      // fread( &this->id        , sizeof(PS::U32)   , 1, fp );
      // fread( &this->istate    , sizeof(PS::U32)   , 1, fp );
      // fread( &this->mass      , sizeof(PS::F64)   , 1, fp );
      // fread( &this->dens      , sizeof(PS::F64)   , 1, fp );
      // fread( &this->pos       , sizeof(PS::F64vec), 1, fp );
      // fread( &this->vel       , sizeof(PS::F64vec), 1, fp );
      // fread( &this->acc       , sizeof(PS::F64vec), 1, fp );
      // fread( &this->sig       , sizeof(PS::F64mat), 1, fp );
      // fread( &this->eps       , sizeof(PS::F64mat), 1, fp );
      // fread( &this->epsp      , sizeof(PS::F64mat), 1, fp );
      // fread( &this->epse      , sizeof(PS::F64mat), 1, fp );
      // fread( &this->dens_rate , sizeof(PS::F64)   , 1, fp );
      // fread( &this->sig_rate  , sizeof(PS::F64mat), 1, fp );
      // fread( &this->eps_rate  , sizeof(PS::F64mat), 1, fp );
      // fread( &this->omg_rate  , sizeof(PS::F64vec), 1, fp );
      // fread( &this->dens_half , sizeof(PS::F64)   , 1, fp );
      // fread( &this->vel_half  , sizeof(PS::F64vec), 1, fp );
      // fread( &this->sig_half  , sizeof(PS::F64mat), 1, fp );
      // fread( &this->al_phi    , sizeof(PS::F64)   , 1, fp );
      // fread( &this->kc        , sizeof(PS::F64)   , 1, fp );
      // fread( &this->coh       , sizeof(PS::F64)   , 1, fp );
      // fread( &this->Ey        , sizeof(PS::F64)   , 1, fp );
      // fread( &this->prank     , sizeof(PS::U32)   , 1, fp );
      // fread( &this->R_tensor  , sizeof(PS::F64mat), 1, fp );
      // fread( &this->disp      , sizeof(PS::F64vec), 1, fp );
  }
  // void writeAscii(FILE* fp) const{
  //   fprintf(fp, "%d,%d,%lf,%lf,%lf,%lf\n",
  //   this->id, this->itype, this->pos.x, this->pos.y, this->pos.z, this->dens);
  // }
  void writeAscii(FILE* fp) const{
    if ( itype  < S_W )
    fprintf(fp,"%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
    this->itype, this->id, this->istate, this->dens, this->pos.x, this->pos.y, this->pos.z,this->vel.x, this->vel.y, this->vel.z, 
    this->sig.xx, this->sig.yy, this->sig.zz, this->sig.xy, this->sig.xz, this->sig.yz, 
    this->epse.xx, this->epse.yy, this->epse.zz, this->epse.xy, this->epse.xz, this->epse.yz,  
    this->epsp.xx, this->epsp.yy, this->epsp.zz, this->epsp.xy, this->epsp.xz, this->epsp.yz, 
    this->omg_rate.x, this->omg_rate.y, this->omg_rate.z, this->acc.x, this->acc.y, this->acc.z, this->disp.x);
  }
  void writeBinary(FILE* fp) const{
    if ( itype < S_W ){
      fwrite(this, sizeof(FP), 1, fp);
      // PS::U32 u32_lst[4] = {this->id, this->bdry, this->itype, this->istate}
      // PS::F64 f64_lst[3] =  {this->mass, this->dens, this->dens_rate}
      // fwrite( &this->id        , sizeof(PS::U32)   , 1, fp );
      // fwrite( &this->bdry      , sizeof(PS::U32)   , 1, fp );
      // fwrite( &this->itype     , sizeof(PS::U32)   , 1, fp );
      // fwrite( &this->istate    , sizeof(PS::U32)   , 1, fp );
      // fwrite( &this->mass      , sizeof(PS::F64)   , 1, fp );
      // fwrite( &this->dens      , sizeof(PS::F64)   , 1, fp );
      // fwrite( &this->pos       , sizeof(PS::F64vec), 1, fp );
      // fwrite( &this->vel       , sizeof(PS::F64vec), 1, fp );
      // fwrite( &this->acc       , sizeof(PS::F64vec), 1, fp );
      // fwrite( &this->sig       , sizeof(PS::F64mat), 1, fp );
      // fwrite( &this->eps       , sizeof(PS::F64mat), 1, fp );
      // fwrite( &this->epsp      , sizeof(PS::F64mat), 1, fp );
      // fwrite( &this->epse      , sizeof(PS::F64mat), 1, fp );
      // fwrite( &this->dens_rate , sizeof(PS::F64)   , 1, fp );
      // fwrite( &this->sig_rate  , sizeof(PS::F64mat), 1, fp );
      // fwrite( &this->eps_rate  , sizeof(PS::F64mat), 1, fp );
      // fwrite( &this->omg_rate  , sizeof(PS::F64vec), 1, fp );
      // fwrite( &this->dens_half , sizeof(PS::F64)   , 1, fp );
      // fwrite( &this->vel_half  , sizeof(PS::F64vec), 1, fp );
      // fwrite( &this->sig_half  , sizeof(PS::F64mat), 1, fp );
      // fwrite( &this->al_phi    , sizeof(PS::F64)   , 1, fp );
      // fwrite( &this->kc        , sizeof(PS::F64)   , 1, fp );
      // fwrite( &this->coh       , sizeof(PS::F64)   , 1, fp );
      // fwrite( &this->Ey        , sizeof(PS::F64)   , 1, fp );
      // fwrite( &this->prank     , sizeof(PS::U32)   , 1, fp );
      // fwrite( &this->R_tensor  , sizeof(PS::F64mat), 1, fp );
      // fwrite( &this->disp      , sizeof(PS::F64vec), 1, fp );
    }
  }


  void CalcStressRateDP_3d(const PS::F64 dt){
    /***** calc trial stress_rate (linear elastic) *****/
    PS::F64mat stress_rate;
    const PS::F64 tr_eps = eps_rate.getTrace();
    PS::F64 temp = Hook1*Ey * tr_eps;
    stress_rate = product_02_3d(temp, del_3d) + product_02_3d(Hook2*Ey, eps_rate);

    /***** calculate yield function: yieldf  *****/
    PS::F64mat sig_tr = sig + product_02_3d(dt, stress_rate);
    PS::F64 I1 = sig_tr.getTrace();
    temp = I1 * dim_1;
    PS::F64mat sij = sig_tr - product_02_3d(temp, del_3d);// stress deviation
    PS::F64 J2 = J2_3d(sij);
    PS::F64 yieldf = sqrt(J2) + al_phi * I1 - kc * coh;
    sig_rate = stress_rate;
    istate = ELASTIC;
    // if (yieldf < zero) {
    //   sig_rate = stress_rate;
    //   istate = ELASTIC;
    // }
    if (yieldf >= zero) {
      ////// plastic eps
      // ++count;
      istate = PLASTIC;
      I1 = sig.getTrace();
      temp = I1 * dim_1;
      sij = sig - product_02_3d(temp, del_3d);
      J2 = J2_3d(sij);
      if( J2 != 0.0 ){
        PS::F64 s_eps = sij.xx * eps_rate.xx + sij.yy * eps_rate.yy + sij.zz * eps_rate.zz 
                        + 2.0 * (sij.xy * eps_rate.xy + sij.xz * eps_rate.xz + sij.yz * eps_rate.yz);
        PS::F64 G_J2 = G*Ey / sqrt(J2);
        PS::F64 lamb_rate = (3.0 * al_phi * K*Ey * tr_eps + G_J2 * s_eps) / (27.0 * al_phi * K*Ey * sin_psi + G*Ey);
        temp = lamb_rate * 9.0 * K*Ey * sin_psi;
        PS::F64 temp2 = lamb_rate * G_J2;
        sig_rate = stress_rate - product_02_3d(temp, del_3d) - product_02_3d(temp2, sij);
      }
    }
    /***** Jaumann stress rate *****/
    sig_rate.xx += 2.0 * ( sig.xy * omg_rate.x + sig.xz * omg_rate.y);
    sig_rate.yy += 2.0 * (-sig.xy * omg_rate.x + sig.yz * omg_rate.z);
    sig_rate.zz += 2.0 * (-sig.xz * omg_rate.y - sig.yz * omg_rate.z);
    sig_rate.xy += omg_rate.x * (sig.yy - sig.xx) + sig.xz * omg_rate.z + sig.yz * omg_rate.y;
    sig_rate.xz += omg_rate.y * (sig.zz - sig.xx) + sig.yz * omg_rate.x - sig.xy * omg_rate.z;
    sig_rate.yz += omg_rate.z * (sig.zz - sig.yy) - sig.xy * omg_rate.y - sig.xz * omg_rate.x;
  }

  void ReturnMapping (const PS::F64 dt){
    PS::F64 I1 = sig.getTrace();
    PS::F64 temp = I1 * dim_1;
    PS::F64mat sij = sig - product_02_3d(temp, del_3d);
    PS::F64 J2 = J2_3d(sij);
    ////tension cracking treatment
    temp = al_phi * I1 - kc * coh;
    if( temp > zero ){
      sig -= product_02_3d(temp / al_phi * dim_1, del_3d);
      I1 = sig.getTrace();
      sij = sig - product_02_3d((I1 * dim_1), del_3d);
      J2 = J2_3d(sij);
    }

    temp = -al_phi * I1 + kc * coh;
    ////stress-scaling back procedure
    if( temp < sqrt(J2) ){
      PS::F64 rn = temp / sqrt(J2);
      //sig = -al_phi*I1+kc*coh0;
      sig = product_02_3d(rn, sij) + product_02_3d(I1 * dim_1, del_3d);
      I1 = sig.getTrace();
      sij = sig - product_02_3d((I1 * dim_1), del_3d);
      J2 = J2_3d(sij);
    }
    epse = product_02_3d(EComp1/Ey, sij) + product_02_3d(EComp2/Ey * I1, del_3d);
    epsp = eps - epse;
    //PS::F64 f = sqrt(J2)+al_phi*I1-kc*coh0;
    //std::cout<<f<<", "<<std::flush;
  }
};

struct EP{
  PS::U32    id;
  PS::U32    itype;
  PS::F64    mass;
  PS::F64    dens;
  PS::F64vec pos;
  PS::F64vec vel;
  PS::F64mat sig;
  PS::F64mat R_tensor;
  PS::F64vec disp;

  PS::F64vec getPos()const{
    return this->pos;
  }
  void copyFromFP(const FP& rp){
    this->pos      = rp.pos;
    this->vel      = rp.vel;
    this->mass     = rp.mass;
    this->dens     = rp.dens;
    this->id       = rp.id;
    this->itype    = rp.itype;
    this->sig      = rp.sig;
    this->R_tensor = rp.R_tensor;
    this->disp     = rp.disp;
  }
  void setPos(const PS::F64vec& pos){
    this->pos   = pos;
  }
  PS::F64 getRSearch() const{
    return kappa * smth_len;
  }
};

