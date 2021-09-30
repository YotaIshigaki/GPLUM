#include "config.h"

double PeriodicDistance(const double x1, const double x2, const int dim){
    double X = x1-x2;
#ifdef PERIODIC_RUN
    if(X>Pall.Lboxh[dim]){
        X -= Pall.Lbox[dim];
    } else if(X<-Pall.Lboxh[dim]){
        X += Pall.Lbox[dim];
    }
#endif
    return X;
}


void PeriodicWrappingi(double Pos[restrict]){

#ifdef PERIODIC_RUN
#if (DIMENSION >= 1)
    if(Pos[0] >= Pall.Lbox[0]) Pos[0] -= Pall.Lbox[0];
    if(Pos[0] < 0.e0)          Pos[0] += Pall.Lbox[0];
#endif

#if (DIMENSION >= 2)
    if(Pos[1] >= Pall.Lbox[1]) Pos[1] -= Pall.Lbox[1];
    if(Pos[1] < 0.e0)          Pos[1] += Pall.Lbox[1];
#endif

#if (DIMENSION >= 3)
    if(Pos[2] >= Pall.Lbox[2]) Pos[2] -= Pall.Lbox[2];
    if(Pos[2] < 0.e0)          Pos[2] += Pall.Lbox[2];
#endif
#endif
    return;
}

#if 0
void PeriodicWrappingi(double Pos[restrict]){

#ifdef PERIODIC_RUN
#if (DIMENSION >= 1)
    if(Pos[0] - Pall.BoxCenter[0] >= Pall.Lboxh[0]) Pos[0] -= Pall.Lbox[0];
    if(Pos[0] - Pall.BoxCenter[0] <  Pall.Lboxh[0]) Pos[0] += Pall.Lbox[0];
#endif

#if (DIMENSION >= 2)
    if(Pos[1] - Pall.BoxCenter[1] >= Pall.Lboxh[1]) Pos[1] -= Pall.Lbox[1];
    if(Pos[1] - Pall.BoxCenter[1] <  Pall.Lboxh[1]) Pos[1] += Pall.Lbox[1];
#endif

#if (DIMENSION >= 3)
    if(Pos[2] - Pall.BoxCenter[2] >= Pall.Lboxh[2]) Pos[2] -= Pall.Lbox[2];
    if(Pos[2] - Pall.BoxCenter[2] <  Pall.Lboxh[2]) Pos[2] += Pall.Lbox[2];
#endif
#endif
    return;
}
#endif

#if 0
void PeriodicWrappingi(double Pos[restrict], double Center[restrict]){

#ifdef PERIODIC_RUN
    while(Pos[0] > Center[0] + Pall.Lboxh[0]) Pos[0] -= Pall.Lbox[0];
    while(Pos[0] < Center[0] - Pall.Lboxh[0]) Pos[0] += Pall.Lbox[0];

    while(Pos[1] > Center[1] + Pall.Lboxh[1]) Pos[1] -= Pall.Lbox[1];
    while(Pos[1] < Center[1] - Pall.Lboxh[1]) Pos[1] += Pall.Lbox[1];

    while(Pos[2] > Center[2] + Pall.Lboxh[2]) Pos[2] -= Pall.Lbox[2];
    while(Pos[2] < Center[2] - Pall.Lboxh[2]) Pos[2] += Pall.Lbox[2];
    /*
    if(Pos[0] > Center[0] + Pall.Lboxh[0]) Pos[0] -= Pall.Lbox[0];
    if(Pos[0] < Center[0] - Pall.Lboxh[0]) Pos[0] += Pall.Lbox[0];

    if(Pos[1] > Center[1] + Pall.Lboxh[1]) Pos[1] -= Pall.Lbox[1];
    if(Pos[1] < Center[1] - Pall.Lboxh[1]) Pos[1] += Pall.Lbox[1];

    if(Pos[2] > Center[2] + Pall.Lboxh[2]) Pos[2] -= Pall.Lbox[2];
    if(Pos[2] < Center[2] - Pall.Lboxh[2]) Pos[2] += Pall.Lbox[2];
    */
#endif
    return;
}
#endif

void PeriodicWrapping(void){

#ifdef PERIODIC_RUN
    for(int i=0;i<Pall.Ntotal;i++){
        //PeriodicWrappingi(Pbody[i]->Pos,Pall.BoxCenter);
        PeriodicWrappingi(Pbody[i]->Pos);
        PeriodicWrappingi(Pbody[i]->PosP);
        if(Pbody[i]->Type == TypeHydro){
            PeriodicWrappingi(PbodyHydro(i)->PosP);
        }
    }
#endif
    return;
}

