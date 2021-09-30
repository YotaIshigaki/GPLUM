#include "config.h"
#include "NeighborSearch.h"

/*
 * This file includes routines which check the members of structrures.
 * These are usuful for debugging.
 */

/*
 * This is the core part of the subroutine CheckHydroStructures(const int mode).
 */
static void CheckHydroStructuresEngine(const int Index){

    for(int k=0;k<3;k++){
        if((fpclassify(PhydroBody(Index)->Pos[k]) == FP_INFINITE)||
           (fpclassify(PhydroBody(Index)->Pos[k]) == FP_NAN)){
            fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s[%d] is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                    PhydroBody(Index)->GlobalID,Index,"Pos",k,PhydroBody(Index)->Pos[k]);
            StructureReportPhydro(Index);
            fflush(NULL);
            assert(fpclassify(PhydroBody(Index)->Pos[k]) != FP_INFINITE);
            assert(fpclassify(PhydroBody(Index)->Pos[k]) != FP_NAN);
        }

        if((fpclassify(PhydroBody(Index)->Vel[k]) == FP_INFINITE)||
           (fpclassify(PhydroBody(Index)->Vel[k]) == FP_NAN)){
            fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s[%d] is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                    PhydroBody(Index)->GlobalID,Index,"Vel",k,PhydroBody(Index)->Vel[k]);
            StructureReportPhydro(Index);
            fflush(NULL);
            assert(fpclassify(PhydroBody(Index)->Vel[k]) != FP_INFINITE);
            assert(fpclassify(PhydroBody(Index)->Vel[k]) != FP_NAN);
        }

        if((fpclassify(Phydro[Index]->HydroAcc[k]) == FP_INFINITE)||
           (fpclassify(Phydro[Index]->HydroAcc[k]) == FP_NAN)){
            fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s[%d] is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                    PhydroBody(Index)->GlobalID,Index,"HydroAcc",k,Phydro[Index]->HydroAcc[k]);
            StructureReportPhydro(Index);
            fflush(NULL);
            assert(fpclassify(Phydro[Index]->HydroAcc[k]) != FP_INFINITE);
            assert(fpclassify(Phydro[Index]->HydroAcc[k]) != FP_NAN);
        }

        if((fpclassify(Phydro[Index]->Rot[k]) == FP_INFINITE)||
           (fpclassify(Phydro[Index]->Rot[k]) == FP_NAN)){
            fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s[%d] is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                    PhydroBody(Index)->GlobalID,Index,"Rot",k,Phydro[Index]->Rot[k]);
            StructureReportPhydro(Index);
            fflush(NULL);
            assert(fpclassify(Phydro[Index]->Rot[k]) != FP_INFINITE);
            assert(fpclassify(Phydro[Index]->Rot[k]) != FP_NAN);
        }
    }

    if((fpclassify(Phydro[Index]->Rho) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Rho) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Rho",Phydro[Index]->Rho);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Rho) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Rho) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->Kernel) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Kernel) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Kernel",Phydro[Index]->Kernel);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Kernel) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Kernel) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->EnergyDensity) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->EnergyDensity) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"EnergyDensity",Phydro[Index]->EnergyDensity);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->EnergyDensity) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->EnergyDensity) != FP_NAN);
    }

#ifdef USE_GRAD_H //{
    if((fpclassify(Phydro[Index]->Gradh) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Gradh) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Gradh",Phydro[Index]->Gradh);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Gradh) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Gradh) != FP_NAN);
    }
#ifdef USE_GRAD_N //{
    if((fpclassify(Phydro[Index]->NumberDensity) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->NumberDensity) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"NumberDensity",Phydro[Index]->NumberDensity);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->NumberDensity) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->NumberDensity) != FP_NAN);
    }
    if((fpclassify(Phydro[Index]->GradN) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->GradN) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"GradN",Phydro[Index]->GradN);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->GradN) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->GradN) != FP_NAN);
    }
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
    
    if((fpclassify(Phydro[Index]->Div) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Div) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Div",Phydro[Index]->Div);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Div) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Div) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->F) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->F) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"F",Phydro[Index]->F);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->F) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->F) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->U) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->U) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"U",Phydro[Index]->U);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->U) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->U) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->Du) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Du) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Du",Phydro[Index]->Du);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Du) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Du) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->DuCooling) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->DuCooling) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"DuCooling",Phydro[Index]->DuCooling);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->DuCooling) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->DuCooling) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->DQheat) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->DQheat) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"DQheat",Phydro[Index]->DQheat);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->DQheat) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->DQheat) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->Z) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Z) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Z",Phydro[Index]->Z);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Z) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Z) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->Vsig) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Vsig) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Vsig",Phydro[Index]->Vsig);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Vsig) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Vsig) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->Alpha) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->Alpha) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"Alpha",Phydro[Index]->Alpha);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->Alpha) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->Alpha) != FP_NAN);
    }

    if((fpclassify(Phydro[Index]->DAlpha) == FP_INFINITE)||
       (fpclassify(Phydro[Index]->DAlpha) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                PhydroBody(Index)->GlobalID,Index,"DAlpha",Phydro[Index]->DAlpha);
        StructureReportPhydro(Index);
        fflush(NULL);
        assert(fpclassify(Phydro[Index]->DAlpha) != FP_INFINITE);
        assert(fpclassify(Phydro[Index]->DAlpha) != FP_NAN);
    }

#ifdef USE_CELIB //{
#endif // USE_CELIB //}
    return ;
}

/*
 * This function checks members in the hydro structure. If there are inf or nan,
 * this function prints the result and calls ``assert''.
 */
void CheckHydroStructures(const int mode){

    if(mode == NONE){ // all hydro structures.
        for(int i=0;i<Pall.Nhydro;i++){
            CheckHydroStructuresEngine(i);
        }
    } else { // a special hydro structure.
        CheckHydroStructuresEngine(mode);
    }

    return ;
}


/*
 * This is the core part of the subroutine CheckBodyStructures(const int mode).
 */
static void CheckBodyStructuresEngine(const int Index){

    for(int k=0;k<3;k++){
        if((fpclassify(Pbody[Index]->Pos[k]) == FP_INFINITE)||
           (fpclassify(Pbody[Index]->Pos[k]) == FP_NAN)){
            fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s[%d] is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                    Pbody[Index]->GlobalID,Index,"Pos",k,Pbody[Index]->Pos[k]);
            StructureReportPbody(Index);
            fflush(NULL);

            assert(fpclassify(Pbody[Index]->Pos[k]) != FP_INFINITE);
            assert(fpclassify(Pbody[Index]->Pos[k]) != FP_NAN);
        }

        if((fpclassify(Pbody[Index]->Vel[k]) == FP_INFINITE)||
           (fpclassify(Pbody[Index]->Vel[k]) == FP_NAN)){
            fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s[%d] is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                    Pbody[Index]->GlobalID,Index,"Pos",k,Pbody[Index]->Vel[k]);
            StructureReportPbody(Index);
            fflush(NULL);
            assert(fpclassify(Pbody[Index]->Vel[k]) != FP_INFINITE);
            assert(fpclassify(Pbody[Index]->Vel[k]) != FP_NAN);
        }

        if((fpclassify(Pbody[Index]->Velh[k]) == FP_INFINITE)||
           (fpclassify(Pbody[Index]->Velh[k]) == FP_NAN)){
            fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s[%d] is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                    Pbody[Index]->GlobalID,Index,"Pos",k,Pbody[Index]->Velh[k]);
            StructureReportPbody(Index);
            fflush(NULL);
            assert(fpclassify(Pbody[Index]->Velh[k]) != FP_INFINITE);
            assert(fpclassify(Pbody[Index]->Velh[k]) != FP_NAN);
        }

        if((fpclassify(Pbody[Index]->Acc[k]) == FP_INFINITE)||
           (fpclassify(Pbody[Index]->Acc[k]) == FP_NAN)){
            fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s[%d] is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                    Pbody[Index]->GlobalID,Index,"Pos",k,Pbody[Index]->Acc[k]);
            StructureReportPbody(Index);
            fflush(NULL);
            assert(fpclassify(Pbody[Index]->Acc[k]) != FP_INFINITE);
            assert(fpclassify(Pbody[Index]->Acc[k]) != FP_NAN);
        }
    }


    if((fpclassify(Pbody[Index]->Mass) == FP_INFINITE)||
       (fpclassify(Pbody[Index]->Mass) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                Pbody[Index]->GlobalID,Index,"EnergyDensity",Pbody[Index]->Mass);
        StructureReportPbody(Index);
        fflush(NULL);
        assert(fpclassify(Pbody[Index]->Mass) != FP_INFINITE);
        assert(fpclassify(Pbody[Index]->Mass) != FP_NAN);
    }

    if((fpclassify(Pbody[Index]->Pot) == FP_INFINITE)||
       (fpclassify(Pbody[Index]->Pot) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                Pbody[Index]->GlobalID,Index,"EnergyDensity",Pbody[Index]->Pot);
        StructureReportPbody(Index);
        fflush(NULL);
        assert(fpclassify(Pbody[Index]->Pot) != FP_INFINITE);
        assert(fpclassify(Pbody[Index]->Pot) != FP_NAN);
    }


    if((fpclassify(Pbody[Index]->Eps) == FP_INFINITE)||
       (fpclassify(Pbody[Index]->Eps) == FP_NAN)){
        fprintf(stderr,"Check[%04d]: GID %ld, Leaf %d |  %s is FP_INFINITE|PF_NAN, %g\n",MPIGetMyID(),
                Pbody[Index]->GlobalID,Index,"EnergyDensity",Pbody[Index]->Eps);
        StructureReportPbody(Index);
        fflush(NULL);
        assert(fpclassify(Pbody[Index]->Eps) != FP_INFINITE);
        assert(fpclassify(Pbody[Index]->Eps) != FP_NAN);
    }


    return ;
}

/*
 * This function checks members in the body structure. If there are inf or nan,
 * this function prints the result and calls ``assert''.
 */
void CheckBodyStructures(const int mode){

    if(mode == NONE){ // all body structures.
        for(int i=0;i<Pall.Ntotal;i++){
            CheckBodyStructuresEngine(i);
        }
    } else { // a special body structure.
        CheckBodyStructuresEngine(mode);
    }

    return ;
}
