#include "config.h"

static inline __attribute__((always_inline)) void *SendBuf(const bool target, void *sendbuf){
    return target
        ? MPI_IN_PLACE
        : sendbuf;
}


// BOUNDARY_CONDITION_TYPE=3  
static bool SphericalAbsorptionBoundary(const int index){

    double R = NORM(Pbody[index]->Pos);

    if(R > BOUNDARY_CONDITION_SPHERICAL_SHELL_EDGE/Pall.UnitLength){
        return true;
    } else {
        return false;
    }
}

// BOUNDARY_CONDITION_TYPE=10  
static void RadiusBaseDragBoundaryConditionForCosmologicalRun(void){

#ifndef COSMOLOGICAL_RUN //{
    return ;
#endif // COSMOLOGICAL_RUN //}

    double ScaleFactor = 1.0/(Pall.Redshift+1.0);
    double R_edge = BOUNDARY_CONDITION_COSMOLOGICAL_SHELL_EDGE*ScaleFactor/Pall.UnitLength;
    double R_edge2 = SQ(R_edge);

    double Hz = CalcHubbleParameterZ();

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type != TypeDM){
            double R2 = NORM2(Pbody[i]->Pos);
            if(R2 < R_edge2) continue;

            double v_hubble[] = {Hz*Pbody[i]->Pos[0],Hz*Pbody[i]->Pos[1],Hz*Pbody[i]->Pos[2]};
            double v_diff[] = {Pbody[i]->Vel[0]-v_hubble[0],
                               Pbody[i]->Vel[1]-v_hubble[1],
                               Pbody[i]->Vel[2]-v_hubble[2]};
            double Factor = BOUNDARY_CONDITION_COSMOLOGICAL_DRAG_FACTOR*(sqrt(R2)/R_edge);
            for(int k=0;k<3;k++){
                Pbody[i]->Vel[k] -= Factor*v_diff[k];
            }
        }
    }


    return ;
}

// BOUNDARY_CONDITION_TYPE=11  
static void VelocityBaseDragBoundaryConditionForCosmologicalRun(void){

#ifndef COSMOLOGICAL_RUN //{
    return ;
#endif // COSMOLOGICAL_RUN //}

    double ScaleFactor = 1.0/(Pall.Redshift+1.0);
    double VelocityEdge = BOUNDARY_CONDITION_COSMOLOGICAL_MAX_VELOCITY;

    double Hz = CalcHubbleParameterZ();
    double H0 = CalcHubbleParameterZ0();

    // fprintf(stderr,"Hz H0 = %g %g \n",Hz,H0);
    double VelConvert = (Pall.UnitLength/Pall.UnitTime/VELOCITY_KMS_CGS);
    double VelReConvert = VELOCITY_KMS_CGS*(Pall.UnitTime/Pall.UnitLength);
    // fprintf(stderr,"Vel Conv = %g, Vel Reconv %g \n",VelConvert,VelReConvert);

    int counter = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type != TypeDM){
            double v_hubble[] = {Hz*Pbody[i]->Pos[0]*Pall.UnitLength/MPC_CGS, // km/s
                                 Hz*Pbody[i]->Pos[1]*Pall.UnitLength/MPC_CGS,
                                 Hz*Pbody[i]->Pos[2]*Pall.UnitLength/MPC_CGS};
            double v_peculiar[] = {Pbody[i]->Vel[0]*VelConvert-v_hubble[0], // km/s
                                   Pbody[i]->Vel[1]*VelConvert-v_hubble[1],
                                   Pbody[i]->Vel[2]*VelConvert-v_hubble[2]};
            double v_norm = NORM(v_peculiar); // km/s
            if(v_norm > VelocityEdge){ // km/s
                double Factor = (1.0-BOUNDARY_CONDITION_COSMOLOGICAL_DRAG_FACTOR);

                fprintf(stderr,"Pre (%g) ",v_norm);
                for(int k=0;k<3;k++){
                    // fprintf(stderr,"[%d][%d] Vori %g ",counter, k, Pbody[i]->Vel[k]);
                    //Pbody[i]->Vel[k] -= Factor*v_diff[k]*VelReConvert;
                    Pbody[i]->Vel[k] = (v_peculiar[k]*Factor+v_hubble[k])*VelReConvert;

                    // fprintf(stderr," %g | %g \n",Pbody[i]->Vel[k],Factor*v_diff[k]*VelReConvert);
                    // fflush(NULL);
                }
                double v_up[] = {Pbody[i]->Vel[0]*VelConvert-v_hubble[0], // km/s
                                 Pbody[i]->Vel[1]*VelConvert-v_hubble[1],
                                 Pbody[i]->Vel[2]*VelConvert-v_hubble[2]};
                fprintf(stderr,"Post (%g) \n",NORM(v_up));

                counter ++;
            }
        }
    }
    

    MPI_Allreduce(MPI_IN_PLACE,&counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        dprintlmpi(counter);
    }
    fflush(NULL);
    // MPI_Barrier(MPI_COMM_WORLD);

    // MPI_Finalize();
    // exit(1);

    return ;
}


static int GetNearIntegerUpper(int IntNumber){
    if(IntNumber <= 0)
        return 0;
    IntNumber -= 1;
    IntNumber |= IntNumber >> 1;
    IntNumber |= IntNumber >> 2;
    IntNumber |= IntNumber >> 4;
    IntNumber |= IntNumber >> 8;
    IntNumber |= IntNumber >> 16;
    return IntNumber + 1;
}

static int GetNearIntegerLower(int IntNumber){
    if(IntNumber <= 0)
        return 0;
    IntNumber -= 1;
    IntNumber |= IntNumber >> 1;
    IntNumber |= IntNumber >> 2;
    IntNumber |= IntNumber >> 4;
    IntNumber |= IntNumber >> 8;
    IntNumber |= IntNumber >> 16;
    return (IntNumber+1)>>1;
}


static bool CheckDragStep(void){

#ifdef USE_DRAG_TERM //{

    int DragStep = (int)(fmin(Pall.dtmax*(Pall.UnitTime/YEAR_CGS),DRAG_TERM_INTERVAL)/(Pall.dtmin*(Pall.UnitTime/YEAR_CGS))+0.5);
    // fprintf(stderr,"++ %g %g %g %d\n",Pall.dtmax*(Pall.UnitTime/YEAR_CGS),DRAG_TERM_COSMOLOGICAL_RUN_INTERVAL,
            // Pall.dtmin*(Pall.UnitTime/YEAR_CGS),DragStep);

    DragStep = MAX(GetNearIntegerLower(DragStep),1);

    //fprintf(stderr,"-- %d %d\n",DragStep,Pall.TStep);
    
    if(Pall.TStep%DragStep == 0){ 
        return true;
    } else {
        return false;
    }

#else // USE_DRAG_TERM //}//{
    return true;
#endif // USE_DRAG_TERM //}
}


static void DragTermForCosmologicalRun(void){

#ifdef COSMOLOGICAL_RUN //{
#ifdef USE_DRAG_TERM //{

    if(CheckDragStep() == false) return;


    double ScaleFactor = 1.0/(Pall.Redshift+1.0);
    double VelocityEdge = BOUNDARY_CONDITION_COSMOLOGICAL_MAX_VELOCITY;

    double Hz = CalcHubbleParameterZ();
    double H0 = CalcHubbleParameterZ0();

    // fprintf(stderr,"Hz H0 = %g %g \n",Hz,H0);
    double VelConvert = (Pall.UnitLength/Pall.UnitTime/VELOCITY_KMS_CGS);
    double VelReConvert = VELOCITY_KMS_CGS*(Pall.UnitTime/Pall.UnitLength);
    // fprintf(stderr,"Vel Conv = %g, Vel Reconv %g \n",VelConvert,VelReConvert);

    int counter = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        if((Pbody[i]->Active)&&(Pbody[i]->Type != TypeDM)){
            double v_hubble[] = {Hz*Pbody[i]->Pos[0]*Pall.UnitLength/MPC_CGS, // km/s
                                 Hz*Pbody[i]->Pos[1]*Pall.UnitLength/MPC_CGS,
                                 Hz*Pbody[i]->Pos[2]*Pall.UnitLength/MPC_CGS};
            double v_peculiar[] = {Pbody[i]->Vel[0]*VelConvert-v_hubble[0], // km/s
                                   Pbody[i]->Vel[1]*VelConvert-v_hubble[1],
                                   Pbody[i]->Vel[2]*VelConvert-v_hubble[2]};
            double v_norm = NORM(v_peculiar); // km/s
            if(v_norm > VelocityEdge){ // km/s
                double Factor = (1.0-BOUNDARY_CONDITION_COSMOLOGICAL_DRAG_FACTOR);

                fprintf(stderr,"Pre (%g) ",v_norm);
                for(int k=0;k<3;k++){
                    // fprintf(stderr,"[%d][%d] Vori %g ",counter, k, Pbody[i]->Vel[k]);
                    //Pbody[i]->Vel[k] -= Factor*v_diff[k]*VelReConvert;
                    Pbody[i]->Vel[k] = (v_peculiar[k]*Factor+v_hubble[k])*VelReConvert;


                    // fprintf(stderr," %g | %g \n",Pbody[i]->Vel[k],Factor*v_diff[k]*VelReConvert);
                    // fflush(NULL);
                }
                double v_up[] = {Pbody[i]->Vel[0]*VelConvert-v_hubble[0], // km/s
                                 Pbody[i]->Vel[1]*VelConvert-v_hubble[1],
                                 Pbody[i]->Vel[2]*VelConvert-v_hubble[2]};
                fprintf(stderr,"Post (%g) \n",NORM(v_up));

                counter ++;
            }
        }
    }
    

    MPI_Allreduce(MPI_IN_PLACE,&counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Drag term affetcts %d particles\n",counter);
    }
    fflush(NULL);
    // MPI_Barrier(MPI_COMM_WORLD);

#endif // USE_DRAG_TERM //}
#endif // COSMOLOGICAL_RUN //}

    return ;
}

static void DragTermForNonCosmologicalRun(void){

#ifndef COSMOLOGICAL_RUN //{
#ifdef USE_DRAG_TERM //{

#ifdef USE_BARYON_COMM //{
    if(MPI_BARYON_COMM_WORLD == MPI_COMM_NULL){
        return ;
    }
#endif // USE_BARYON_COMM //{


    if(CheckDragStep() == false) return;

    const double VelocityEdge = BOUNDARY_CONDITION_COSMOLOGICAL_MAX_VELOCITY;

    const double VelConvert = (Pall.UnitLength/Pall.UnitTime/VELOCITY_KMS_CGS);
    const double VelReConvert = VELOCITY_KMS_CGS*(Pall.UnitTime/Pall.UnitLength);

    int counter = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        if((Pbody[i]->Active)&&(Pbody[i]->Type != TypeDM)){
            double v_peculiar[] = {Pbody[i]->Vel[0]*VelConvert, // km/s
                                   Pbody[i]->Vel[1]*VelConvert,
                                   Pbody[i]->Vel[2]*VelConvert};
            double v_norm = NORM(v_peculiar); // km/s
            if(v_norm > VelocityEdge){ // km/s
                double Factor = (1.0-BOUNDARY_CONDITION_COSMOLOGICAL_DRAG_FACTOR);

                fprintf(stderr,"Pre (%g) ",v_norm);
                for(int k=0;k<3;k++){
                    Pbody[i]->Vel[k] = (v_peculiar[k]*Factor)*VelReConvert;
                }
                double v_up[] = {Pbody[i]->Vel[0]*VelConvert, // km/s
                                 Pbody[i]->Vel[1]*VelConvert,
                                 Pbody[i]->Vel[2]*VelConvert};
                fprintf(stderr,"Post (%g) \n",NORM(v_up));

                counter ++;
            }
        }
    }
    
#ifdef USE_BARYON_COMM //{
    MPI_Reduce(SendBuf(MPIGetBaryonMyID()==MPI_ROOT_RANK,&counter),&counter,1,MPI_INT,MPI_SUM,MPI_ROOT_RANK,MPI_BARYON_COMM_WORLD);
    if(MPIGetBaryonMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Drag term affetcts %d particles\n",counter);
    }
#else // USE_BARYON_COMM //}//{
    MPI_Reduce(SendBuf(MPIGetMyID()==MPI_ROOT_RANK,&counter),&counter,1,MPI_INT,MPI_SUM,MPI_ROOT_RANK,MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Drag term affetcts %d particles\n",counter);
    }
#endif // USE_BARYON_COMM //{
    fflush(NULL);


#endif // USE_DRAG_TERM //}
#endif // COSMOLOGICAL_RUN //}

    return ;
}


void AddDragTerm(void){

#ifdef USE_DRAG_TERM //{

#   ifdef COSMOLOGICAL_RUN //{
    DragTermForCosmologicalRun();
#   endif // COSMOLOGICAL_RUN //}

#   ifndef COSMOLOGICAL_RUN //{
    DragTermForNonCosmologicalRun();
#   endif // COSMOLOGICAL_RUN //}

#endif // USE_DRAG_TERM //}
    return ;
}


void ImposeBoundaryCondition(const int mode){

#ifdef USE_BOUNDARY_CONDITION //{
    switch(mode){
        case 0: // Periodic boundary
            break;
        case 1: // Box boundary with reflecting walls.
            break;
        case 2: // Absorption boundary with a box shape.
            break;
        case 3: // Absorption boundary with a spherical shell.
            ParticleRemoverArbitraryCriteria(SphericalAbsorptionBoundary);
            break;
        case 10: // Add drag term for baryon particles out of a spherical region.
            RadiusBaseDragBoundaryConditionForCosmologicalRun();
            break;
        case 11: // Add drag term for baryon particles with unreasonably very high speed.
            VelocityBaseDragBoundaryConditionForCosmologicalRun();
            break;
    }
#endif // USE_BOUNDARY_CONDITION //}
    return ;
}

