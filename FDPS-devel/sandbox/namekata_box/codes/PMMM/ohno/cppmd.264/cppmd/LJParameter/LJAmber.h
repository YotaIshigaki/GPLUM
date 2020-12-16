#ifndef LJAMBER_H
#define LJAMBER_H

#include "UnitParameter.h"
#include "LJParameter.h"

/* LJ parameters of AMBER */
/* 
   2^(1/6) * (sigma_i + sigma_j)/2 = R_i + R_j
             (sigma_i + sigma_j)/2 = 2^(-1/6) * (R_i + R_j) 

   2^(1/6) * sigma_i / 2 = R_i
   sigma_i = R_i * 2.0 * 2^(-1/6)

   A = 4*epsilon*{(sigma_i+sigma_j)/2}^12
     =   epsilon*(R_i+R_j)^12
   B = 4*epsilon*{(sigma_i+sigma_j)/2}^6
     = 2*epsilon*(R_i+R_j)^6

   E = A*r^(-12) - B*r^(-6)

   dE/dr = -12*A*r^(-13) + 6*B*r^(-7)
         = -6*r^(-7)*{2*A*r^(-6)-B}
         = 0 at  r^6 = 2*A/B = 2*{(sigma_i+sigma_j)/2}^6 = (R_i + R_j)^6
 */
//! Lennard Jones parameters for Amber format
struct LJAmberParameter {
  std::string name;
  double r;
  double epsilon;
};

//! Lennard Jones for Amber format
class LJAmber {
 public:
  LJAmberParameter *ljamberparameters;
  int num_lj;
  
  LJAmber() {num_lj=1; ljamberparameters = new LJAmberParameter[num_lj];}
  LJAmber(LJAmberParameter* ljamber, int num) : ljamberparameters() {
    num_lj=num;
    ljamberparameters = ljamber;
  }

    bool get_lj(std::string name, double& r, double& epsilon){
    int i;
    for(i=0;i<num_lj;i++){
      if(name!=ljamberparameters[i].name)break;
    }
    if(i<num_lj){
      r = ljamberparameters[i].r;
      epsilon = ljamberparameters[i].epsilon;
      return true;
    }else{
      r = 0.0;
      epsilon = 0.0;
      return false;
    }
  }

  LJMixParameter getLJMixParameter(long i, long j){
    LJMixParameter ljmp;
    double r = ljamberparameters[i].r + ljamberparameters[j].r;
    double epsilon = UnitParameter::normalizeEnergy_kcal_mol(sqrt(ljamberparameters[i].epsilon*ljamberparameters[j].epsilon));
    double r2 = r*r;
    double r6 = r2*r2*r2;
    double r12 = r6*r6;
    ljmp.potFirst    =      epsilon*r12;
    ljmp.potSecond   =  2.0*epsilon*r6;
    ljmp.forceFirst  = 12.0*epsilon*r12;
    ljmp.forceSecond = 12.0*epsilon*r6;
    return ljmp;
  }
  
  long convertLJMixparameterArray(LJMixparameterArray& ljma){
    ljma.resize2d(num_lj,num_lj);
    for(long i=0;i<num_lj;i++){
      for(long j=0;j<num_lj;j++){
        ljma[i][j] = getLJMixParameter(i,j);
      }
    }
    return num_lj;
  }
  long convertLJMixName(LJMixName& ljmn){
    ljmn.resize2d(num_lj,num_lj);
    for(long i=0;i<num_lj;i++){
      for(long j=0;j<num_lj;j++){
        ljmn[i][j] = ljamberparameters[i].name + "-" + ljamberparameters[j].name;
      }
    }
    return num_lj;
  }
  LJParameter getLJParameter(int i=0)
  {
    if(i>=num_lj)i=num_lj-1;
    LJParameter ljp(ljamberparameters[i].name,
                    ljamberparameters[i].r*pow(2.0,5.0/6.0),
                    ljamberparameters[i].epsilon);
    return ljp;
  }
  LJParameter getmixParameter(int i=0, int j=0)
  {
    if(i>=num_lj)i=num_lj-1;
    if(j>=num_lj)j=num_lj-1;
    LJParameter ljmix = LJParameterStuff::mixParameter(getLJParameter(i),getLJParameter(j));
    return ljmix;
  }

};

#endif
