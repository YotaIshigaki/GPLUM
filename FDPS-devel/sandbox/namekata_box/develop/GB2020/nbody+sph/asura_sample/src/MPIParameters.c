#include "config.h"
#include "MPIParameters.h"
#include "CommunicationTable.h"

static int numprocs,myid,namelen,numgrapes,numprocspower;
static char processor_name[MPI_MAX_PROCESSOR_NAME];
static int Factors[3] = {0,0,0};

void InitializeMPIEnv(int *argc, char **argv){

    MPI_Init(argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);

    MPISetMyID(myid);
    MPISetNumProcs(numprocs);
#if DECOMPOSITION_TYPE == 0 //{
    MPISetNumProcsPower(numprocs);
#elif DECOMPOSITION_TYPE == 1 //}//{
    MPISetXYZFactors();
#endif // DECOMPOSITION_TYPE //}
    MPISetNumGrapes(MIN(MPIGetNumProcs(),4));
    // MPISetNameLen(namelen);
    // MPISetProcessorName(processor_name);

    return ;
}


void FinalizeMPIEnv(void){

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return ;
}


void MPISetXYZFactors(void){

#if DIMENSION == 1
    MPI_Dims_create(numprocs,1,Factors);
#elif DIMENSION == 2
    MPI_Dims_create(numprocs,2,Factors);
#elif DIMENSION == 3
    MPI_Dims_create(numprocs,3,Factors);
#else 
    fprintf(stderr,"Wrong dimension was used. DIMENSION == %d\n",DIMENSION);
    fprintf(stderr,"Need to change the DIMENSION flag\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
#endif 

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"The result of the node fractionization is ");
        for(int k=0;k<DIMENSION;k++){
            if(k==DIMENSION-1){
                fprintf(stderr,"%d.\n",Factors[k]);
            } else {
                fprintf(stderr,"%d*",Factors[k]);
            }
        }
        fflush(NULL);
    }

    return ;
}


void MPISetMyID(int current_id){
    myid = current_id;
    return;
}

void MPISetNumProcs(int current_numprocs){
    numprocs = current_numprocs;
    return;
}

void MPISetNumGrapes(int current_numgrapes){
    numgrapes = current_numgrapes;
    return;
}

void MPISetNameLen(int current_namelen){
    namelen = current_namelen;
    return;
}

void MPISetProcessorName(char *current_processor_name){
    strcpy(processor_name,current_processor_name);
    return;
}

void MPISetNumProcsPower(int current_numprocs){

    int ndummy = 1;

    numprocspower = 0;
    if(current_numprocs%2 != 0){
        numprocspower = 0;
        return;
    }

    while(ndummy != current_numprocs){
        ndummy *= 2;
        numprocspower ++;
    }

    return;
}

int MPIGetMyID(void){
    return myid;
}

int MPIGetNumProcs(void){
    return numprocs;
}

int MPIGetNumGrapes(void){
    return numgrapes;
}

int MPIGetNameLen(void){
    return namelen;
}

char *MPIGetProcessorName(void){
    return processor_name;
}

int MPIGetNumProcsPower(void){
    return numprocspower;
}

void MPIGetXYZFactors(int FactorsDist[]){
    FactorsDist[0] = Factors[0];
    FactorsDist[1] = Factors[1];
    FactorsDist[2] = Factors[2];
    return ;
}

int MPIGetFactor(const int direction){
    assert(direction <= DIMENSION);
    return Factors[direction];
}

int MPIGetFactorFromID(const int direction){
#if DIMENSION == 1 //{
    assert(direction <= DIMENSION);
    return myid;
#elif DIMENSION == 2//}//{
    if(direction == 0){
        return myid/Factors[1];
    }else if(direction == 1){
        return (myid%Factors[1]);
    } else {
        assert(direction <= DIMENSION);
        return NONE;
    }
#elif DIMENSION == 3//}//{
    if(direction == 0){
        return myid/(Factors[1]*Factors[2]);
    }else if(direction == 1){
        return (myid%(Factors[1]*Factors[2])/Factors[2]);
    }else if(direction == 2){
        return (myid%(Factors[1]*Factors[2])%Factors[2]);
    } else {
        assert(direction <= DIMENSION);
        return NONE;
    }
#endif
}


#ifdef USE_BARYON_COMM //{

static int HydroNumProcs,HydroMyID;
static int ActiveHydroNumProcs,ActiveHydroMyID;

int MPIGetHydroMyID(void){
    return HydroMyID;
}

int MPIGetHydroNumProcs(void){
    return HydroNumProcs;
}

int MPIGetActiveHydroMyID(void){
    return ActiveHydroMyID;
}

int MPIGetActiveHydroNumProcs(void){
    return ActiveHydroNumProcs;
}


void CreateHydroGroup(void){

    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();
 
    int HydroFlags[NProcs];
    for(int i=0;i<NProcs;i++) HydroFlags[i] = 0;

    if(Pall.Nhydro > 0)
        HydroFlags[MyID] = 1;

    MPI_Allreduce(MPI_IN_PLACE,HydroFlags,NProcs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    
    int Group[NProcs];

    int counter = 0;
    for(int i=0;i<NProcs;i++){
        if(HydroFlags[i] == 1){
            Group[counter] = i;
            counter ++;
        }
    }
    MPI_Group group_world;

    //get the group under MPI_COMM_WORLD
    MPI_Comm_group(MPI_COMM_WORLD,&group_world);

    // create the new group
    MPI_Group_incl(group_world,counter,Group,&hydro_group);

    // create the new communicator
    MPI_Comm_create(MPI_COMM_WORLD,hydro_group,&MPI_HYDRO_COMM_WORLD);

    if(MPI_HYDRO_COMM_WORLD != MPI_COMM_NULL){
        MPI_Comm_rank(MPI_HYDRO_COMM_WORLD,&HydroMyID);
        MPI_Comm_size(MPI_HYDRO_COMM_WORLD,&HydroNumProcs);
    } else {
        HydroMyID = NONE;
        HydroNumProcs = NONE;
    }

    InitializeHydroCommunicationOrder();

    return ;
}

void DeleteHydroGroup(void){

    if(MPI_HYDRO_COMM_WORLD != MPI_COMM_NULL){
        MPI_Comm_free(&MPI_HYDRO_COMM_WORLD);
        MPI_Group_free(&hydro_group);
    }
    return ;
}

//


static void CheckActiveDomain(int HydroFlags[], const int MyID, const int NProcs){

    if(Pall.NActivesHydro > 0)
        HydroFlags[MyID] = 1;

    if(Pall.Nhydro > 0)
        HydroFlags[MyID] = 1;

    // Check active hydro domain[local] and hydro domain[]


    // Check active hydro domain[local] and hydro domain[]


    return ;
}


void CreateHydroActiveGroup(void){

    if(MPI_HYDRO_COMM_WORLD == MPI_COMM_NULL) 
        return ;

    const int HydroMyID = MPIGetMyID();
    const int HydroNProcs = MPIGetNumProcs();
 
    int HydroFlags[HydroNProcs];
    for(int i=0;i<HydroNProcs;i++) HydroFlags[i] = 0;

    // List up mydomain and the contacted domains.
    CheckActiveDomain(HydroFlags,HydroMyID,HydroNProcs);
#warning need update.

    MPI_Allreduce(MPI_IN_PLACE,HydroFlags,HydroNProcs,MPI_INT,MPI_SUM,MPI_HYDRO_COMM_WORLD);
    
    int Group[HydroNProcs];

    int counter = 0;
    for(int i=0;i<HydroNProcs;i++){
        if(HydroFlags[i] > 1){
            Group[counter] = i;
            counter ++;
        }
    }
    MPI_Group group_world;

    //get the group under MPI_COMM_WORLD
    MPI_Comm_group(MPI_HYDRO_COMM_WORLD,&group_world);

    // create the new group
    MPI_Group_incl(group_world,counter,Group,&activehydro_group);

    // create the new communicator
    MPI_Comm_create(MPI_HYDRO_COMM_WORLD,activehydro_group,
            &MPI_ACTIVEHYDRO_COMM_WORLD);

    if(MPI_ACTIVEHYDRO_COMM_WORLD != MPI_COMM_NULL){
        MPI_Comm_rank(MPI_ACTIVEHYDRO_COMM_WORLD,&ActiveHydroMyID);
        MPI_Comm_size(MPI_ACTIVEHYDRO_COMM_WORLD,&ActiveHydroNumProcs);
    } else {
        ActiveHydroMyID = NONE;
        ActiveHydroNumProcs = NONE;
    }

    return ;
}

void DeleteActiveHydroGroup(void){

    if(MPI_ACTIVEHYDRO_COMM_WORLD != MPI_COMM_NULL){
        MPI_Comm_free(&MPI_ACTIVEHYDRO_COMM_WORLD);
        MPI_Group_free(&activehydro_group);
    }
    return ;
}

static int BaryonNumProcs,BaryonMyID;
static int ActiveBaryonNumProcs,ActiveBaryonMyID;

int MPIGetBaryonMyID(void){
    return BaryonMyID;
}

int MPIGetBaryonNumProcs(void){
    return BaryonNumProcs;
}

int MPIGetActiveBaryonMyID(void){
    return ActiveBaryonMyID;
}

int MPIGetActiveBaryonNumProcs(void){
    return ActiveBaryonNumProcs;
}

void CreateBaryonGroup(void){

    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();
 
    int HydroFlags[NProcs];
    for(int i=0;i<NProcs;i++) HydroFlags[i] = 0;

    if((Pall.Nhydro+Pall.Nstars+Pall.Nsink) > 0)
        HydroFlags[MyID] = 1;

    MPI_Allreduce(MPI_IN_PLACE,HydroFlags,NProcs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    
    int Group[NProcs];

    int counter = 0;
    for(int i=0;i<NProcs;i++){
        if(HydroFlags[i] == 1){
            Group[counter] = i;
            counter ++;
        }
    }
    MPI_Group group_world;

    //get the group under MPI_COMM_WORLD
    MPI_Comm_group(MPI_COMM_WORLD,&group_world);

    // create the new group
    MPI_Group_incl(group_world,counter,Group,&baryon_group);

    // create the new communicator
    MPI_Comm_create(MPI_COMM_WORLD,baryon_group,&MPI_BARYON_COMM_WORLD);

    if(MPI_BARYON_COMM_WORLD != MPI_COMM_NULL){
        MPI_Comm_rank(MPI_BARYON_COMM_WORLD,&BaryonMyID);
        MPI_Comm_size(MPI_BARYON_COMM_WORLD,&BaryonNumProcs);
    } else {
        BaryonMyID = NONE;
        BaryonNumProcs = NONE;
    }

    InitializeBaryonCommunicationOrder();

    return ;
}

void DeleteBaryonGroup(void){

    if(MPI_BARYON_COMM_WORLD != MPI_COMM_NULL){
        MPI_Comm_free(&MPI_BARYON_COMM_WORLD);
        MPI_Group_free(&baryon_group);
    }
    return ;
}

void CreateActiveBaryonGroup(void){

    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();
 
    int HydroFlags[NProcs];
    for(int i=0;i<NProcs;i++) HydroFlags[i] = 0;

    if((Pall.NActivesHydro+Pall.NActivesStars+Pall.NActivesSink) > 0)
        HydroFlags[MyID] = 1;

    MPI_Allreduce(MPI_IN_PLACE,HydroFlags,NProcs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    int Group[NProcs];

    int counter = 0;
    for(int i=0;i<NProcs;i++){
        if(HydroFlags[i] == 1){
            Group[counter] = i;
            counter ++;
        }
    }
    MPI_Group group_world;

    //get the group under MPI_COMM_WORLD
    MPI_Comm_group(MPI_COMM_WORLD,&group_world);

    // create the new group
    MPI_Group_incl(group_world,counter,Group,&activebaryon_group);

    // create the new communicator
    MPI_Comm_create(MPI_COMM_WORLD,activebaryon_group,&MPI_ACTIVEBARYON_COMM_WORLD);

    if(MPI_ACTIVEBARYON_COMM_WORLD != MPI_COMM_NULL){
        MPI_Comm_rank(MPI_ACTIVEBARYON_COMM_WORLD,&ActiveBaryonMyID);
        MPI_Comm_size(MPI_ACTIVEBARYON_COMM_WORLD,&ActiveBaryonNumProcs);
    } else {
        ActiveBaryonMyID = NONE;
        ActiveBaryonNumProcs = NONE;
    }

    InitializeActiveBaryonCommunicationOrder();

    return ;
}

void DeleteActiveBaryonGroup(void){

    if(MPI_ACTIVEBARYON_COMM_WORLD != MPI_COMM_NULL){
        MPI_Comm_free(&MPI_ACTIVEBARYON_COMM_WORLD);
        MPI_Group_free(&activebaryon_group);
    }
    return ;
}


int RankTranslator(MPI_Group Group1, int Rank1, MPI_Group Group2){
    int Rank2;
    MPI_Group_translate_ranks(Group1,1,&Rank1,Group2,&Rank2);
    return Rank2;
}
#endif // USE_BARYON_COMM //}

