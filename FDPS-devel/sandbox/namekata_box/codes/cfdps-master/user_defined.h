//user_defined.c


typedef struct full_particle{ //!fdps FP,EPI,EPJ,Force
    //!fdps copyFromForce full_particle (pot,pot) (acc,acc)
    //!fdps copyFromFP full_particle (id,id) (mass,mass) (eps,eps) (pos,pos) 
    //!fdps clear id=keep, mass=keep, eps=keep, pos=keep, vel=keep
    long long id;
    double  mass;      //!fdps charge
    double  eps;
    Cvec_Float64 pos;  //!fdps position
    Cvec_Float64 vel;  //!fdps velocity
    double  pot;
    Cvec_Float64 acc;  
}Full_particle,*PFull_particle;

