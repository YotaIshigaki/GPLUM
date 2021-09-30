#include "config.h"
#include "StructureOperation.h"
#include "OrderingKey.h"

#define StructureStretch (ON)

static void SortStructPbody(void);
static void SortStructPhydro(void);
static void SortStructPstar(void);
static void SortStructPsink(void);
// static void ReconnectBodytoHydro(void);
// static void ReconnectBodytoStar(void);
// static void ReconnectBodytoSink(void);
// static void ReconnectHydrotoBody(void);
// static void ReconnectStartoBody(void);
// static void ReconnectSinktoBody(void);
static void StructureReportPbodyElement(StructPbodyptr Pb);
static void StructureReportPhydroElement(StructPhydroptr Ph);
static void StructureReportPstarElement(StructPstarptr Ps);
static void StructureReportPinkElement(StructPstarptr Ps);


void AddStructPbody(const int nadd){

    int Nsize = MAX(nadd,NAdditionUnit);

    StructPbodyptr Pb = PbodyElements;
    while(Pb->Next != NULL){
        Pb = Pb->Next;        
    }

    StructPbodyptr PbodyNew = (void *)malloc(Nsize*sizeof(StructPbody));
    memset(PbodyNew,0,Nsize*sizeof(StructPbody));

    Pb->Next = PbodyNew;
    for(int i=0;i<Nsize-1;i++)
        PbodyNew[i].Next = &(PbodyNew[i+1]);

    for(int i=0;i<Nsize;i++)
        PbodyNew[i].Use = OFF;

    return;
}


void AddStructPhydro(const int nadd){

    int Nsize = MAX(nadd,NAdditionUnit);

    StructPhydroptr Ph = PhydroElements;
    while(Ph->Next != NULL)
        Ph = Ph->Next;        
    
    StructPhydroptr PhydroNew = (void *)malloc(Nsize*sizeof(StructPhydro));
    memset(PhydroNew,0,Nsize*sizeof(StructPhydro));

    Ph->Next = PhydroNew;
    for(int i=0;i<Nsize-1;i++)
        PhydroNew[i].Next = &(PhydroNew[i+1]);
    PhydroNew[Nsize-1].Next = NULL;
   
    for(int i=0;i<Nsize;i++)
        PhydroNew[i].Use = OFF;

    return;
}

void AddStructPstar(const int nadd){

    int Nsize = MAX(nadd,NAdditionUnit);

    StructPstarptr Ps = PstarElements;
    while(Ps->Next != NULL)
        Ps = Ps->Next;        

    StructPstarptr PstarNew = (void *)malloc(Nsize*sizeof(StructPstar));
    memset(PstarNew,0,Nsize*sizeof(StructPstar));

    Ps->Next = PstarNew;
    for(int i=0;i<Nsize-1;i++)
        PstarNew[i].Next = &(PstarNew[i+1]);
    PstarNew[Nsize-1].Next = NULL;

    for(int i=0;i<Nsize;i++)
        PstarNew[i].Use = OFF;

    return;
}

void AddStructPsink(const int nadd){

    int Nsize = MAX(nadd,NAdditionUnit);

    StructPsinkptr Psk = PsinkElements;
    while(Psk->Next != NULL)
        Psk = Psk->Next;        

    StructPsinkptr PskNew = (void *)malloc(Nsize*sizeof(StructPsink));
    memset(PskNew,0,Nsize*sizeof(StructPsink));

    Psk->Next = PskNew;
    for(int i=0;i<Nsize-1;i++)
        PskNew[i].Next = &(PskNew[i+1]);
    PskNew[Nsize-1].Next = NULL;

    for(int i=0;i<Nsize;i++)
        PskNew[i].Use = OFF;

    return;
}

#define SORT_MORTON     (0)
#define SORT_TIMESTEP   (1)

#define SortStepSkip (100)
static int SortStep = SortStepSkip+1;
void SortAllStructures(void){

    if(SortStep > Pall.TStepTotal)
        return;

#if (StructureStretch)
#if PARTICLE_SORT_TYPE == 0  
    TimingResults.KeyGenerationThisStep = GetElapsedTime();
    MakeOrderingKey();
    TimingResults.KeyGenerationThisStep = GetElapsedTime()-TimingResults.KeyGenerationThisStep;
//#error !
#elif PARTICLE_SORT_TYPE == 1
//#error ?
#endif

    TimingResults.SortStructuresThisStep = GetElapsedTime();
    SortStructPbody();
    SortStructPhydro();
    SortStructPstar();
    SortStructPsink();
    ReConnectPointers();
    TimingResults.SortStructuresThisStep = GetElapsedTime()-TimingResults.SortStructuresThisStep;
#endif
    SortStep = SortStepSkip + Pall.TStepTotal;

    return;
}

void ShuffleStructures(void){

    int Size = PbodyElementsSize;

    for(int i=0;i<Size;i++){
        const int index_0 = gsl_rng_uniform_int(RandomGenerator,Size);
        const int index_1 = gsl_rng_uniform_int(RandomGenerator,Size);
        StructPbody Pb = PbodyElements[index_0];
        PbodyElements[index_0] = PbodyElements[index_1];
        PbodyElements[index_1] = Pb;
    }

    for(int i=0;i<PbodyElementsSize-1;i++)
        PbodyElements[i].Next = &(PbodyElements[i+1]);
    PbodyElements[PbodyElementsSize-1].Next = NULL;

    ReconnectBodytoHydro();
    ReconnectBodytoStar();
    ReconnectBodytoSink();

    return;
}


#if PARTICLE_SORT_TYPE == 0  
static int PbodyKeyCmp(const void *x, const void *y){
    const StructPbodyptr pointer1 = (StructPbodyptr)x;
    const StructPbodyptr pointer2 = (StructPbodyptr)y;

    if( pointer1->Use == OFF && pointer2->Use == OFF ){
        return 0;
    } else if ( pointer1->Use == OFF && pointer2->Use == ON ){
        return 1;
    } else if ( pointer1->Use == ON && pointer2->Use == OFF ){
        return -1;
    } else if ( pointer1->Type < pointer2->Type ){
        return -1;
    } else if ( pointer1->Type > pointer2->Type ){
        return 1;
    } else {
        if( pointer1->OrderingKey > pointer2->OrderingKey)
            return 1;
        else if( pointer1->OrderingKey < pointer2->OrderingKey)
            return -1;
        else
            return 0;
    }
}
#elif PARTICLE_SORT_TYPE == 1
static int PbodyKeyCmp(const void *x, const void *y){
    const StructPbodyptr pointer1 = (StructPbodyptr)x;
    const StructPbodyptr pointer2 = (StructPbodyptr)y;

    if( pointer1->Use == OFF && pointer2->Use == OFF ){
        return 0;
    } else if ( pointer1->Use == OFF && pointer2->Use == ON ){
        return 1;
    } else if ( pointer1->Use == ON && pointer2->Use == OFF ){
        return -1;
    } else if ( pointer1->Type < pointer2->Type ){
        return -1;
    } else if ( pointer1->Type > pointer2->Type ){
        return 1;
    } else {
        if( pointer1->dt > pointer2->dt)
            return 1;
        else if( pointer1->dt < pointer2->dt)
            return -1;
        else
            return 0;
    }
}
#endif

static void SortStructPbody(void){

    //sprintlmpi("Start sort body");
    qsort(PbodyElements, PbodyElementsSize,
        sizeof(StructPbody),(int(*)(const void*, const void*))PbodyKeyCmp);

    for(int i=0;i<PbodyElementsSize-1;i++)
        PbodyElements[i].Next = &(PbodyElements[i+1]);
    PbodyElements[PbodyElementsSize-1].Next = NULL;

    ReconnectBodytoHydro();
    ReconnectBodytoStar();
    ReconnectBodytoSink();
    return;
}

#if PARTICLE_SORT_TYPE == 0  
static int PhydroKeyCmp(const void *x, const void *y){
    const StructPhydroptr pointer1 = (StructPhydroptr)x;
    const StructPhydroptr pointer2 = (StructPhydroptr)y;
    if( pointer1->Use == OFF && pointer2->Use == OFF){
        return 0;
    }else if ( pointer1->Use == OFF && pointer2->Use == ON){
        return 1;
    }else if ( pointer1->Use == ON && pointer2->Use == OFF){
        return -1;
    }else if( pointer1->Body->OrderingKey > pointer2->Body->OrderingKey){
        return 1;
    }else if( pointer1->Body->OrderingKey < pointer2->Body->OrderingKey){
        return -1;
    }else{
        return 0;
    }
}
#elif PARTICLE_SORT_TYPE == 1
static int PhydroKeyCmp(const void *x, const void *y){
    const StructPhydroptr pointer1 = (StructPhydroptr)x;
    const StructPhydroptr pointer2 = (StructPhydroptr)y;
    if( pointer1->Use == OFF && pointer2->Use == OFF){
        return 0;
    }else if ( pointer1->Use == OFF && pointer2->Use == ON){
        return 1;
    }else if ( pointer1->Use == ON && pointer2->Use == OFF){
        return -1;
    }else if( pointer1->Body->dt > pointer2->Body->dt){
        return 1;
    }else if( pointer1->Body->dt < pointer2->Body->dt){
        return -1;
    }else{
        return 0;
    }
}
#endif

static void SortStructPhydro(void){

    if(Pall.Nhydro == 0)
        return;

    //sprintlmpi("Start sort hydro");
    qsort(PhydroElements, PhydroElementsSize,
        sizeof(StructPhydro),(int(*)(const void*, const void*))PhydroKeyCmp);

    for(int i=0;i<PhydroElementsSize-1;i++)
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    PhydroElements[PhydroElementsSize-1].Next = NULL;

    ReconnectHydrotoBody();

    return;
}

#if PARTICLE_SORT_TYPE == 0  
static int PstarKeyCmp(const void *x, const void *y){
    const StructPstarptr pointer1 = (StructPstarptr)x;
    const StructPstarptr pointer2 = (StructPstarptr)y;
    if( pointer1->Use == OFF && pointer2->Use == OFF){
        return 0;
    }else if ( pointer1->Use == OFF && pointer2->Use == ON){
        return 1;
    }else if ( pointer1->Use == ON && pointer2->Use == OFF){
        return -1;
    }else if( pointer1->Body->OrderingKey > pointer2->Body->OrderingKey){
        return 1;
    }else if( pointer1->Body->OrderingKey < pointer2->Body->OrderingKey){
        return -1;
    }else{
        return 0;
    }
}
#elif PARTICLE_SORT_TYPE == 1
static int PstarKeyCmp(const void *x, const void *y){
    const StructPstarptr pointer1 = (StructPstarptr)x;
    const StructPstarptr pointer2 = (StructPstarptr)y;
    if( pointer1->Use == OFF && pointer2->Use == OFF){
        return 0;
    }else if ( pointer1->Use == OFF && pointer2->Use == ON){
        return 1;
    }else if ( pointer1->Use == ON && pointer2->Use == OFF){
        return -1;
    }else if( pointer1->Body->dt > pointer2->Body->dt){
        return 1;
    }else if( pointer1->Body->dt < pointer2->Body->dt){
        return -1;
    }else{
        return 0;
    }
}
#endif

static void SortStructPstar(void){

    if(Pall.Nstars == 0)
        return;

    //sprintlmpi("Start sort star");
    qsort(PstarElements, PstarElementsSize,
        sizeof(StructPstar),(int(*)(const void*, const void*))PstarKeyCmp);

    for(int i=0;i<PstarElementsSize-1;i++)
        PstarElements[i].Next = &(PstarElements[i+1]);
    PstarElements[PstarElementsSize-1].Next = NULL;

    ReconnectStartoBody();

    return;
}

#if PARTICLE_SORT_TYPE == 0  
static int PsinkKeyCmp(const void *x, const void *y){
    const StructPsinkptr pointer1 = (StructPsinkptr)x;
    const StructPsinkptr pointer2 = (StructPsinkptr)y;
    if( pointer1->Use == OFF && pointer2->Use == OFF){
        return 0;
    }else if ( pointer1->Use == OFF && pointer2->Use == ON){
        return 1;
    }else if ( pointer1->Use == ON && pointer2->Use == OFF){
        return -1;
    }else if( pointer1->Body->OrderingKey > pointer2->Body->OrderingKey){
        return 1;
    }else if( pointer1->Body->OrderingKey < pointer2->Body->OrderingKey){
        return -1;
    }else{
        return 0;
    }
}
#elif PARTICLE_SORT_TYPE == 1
static int PsinkKeyCmp(const void *x, const void *y){
    const StructPsinkptr pointer1 = (StructPsinkptr)x;
    const StructPsinkptr pointer2 = (StructPsinkptr)y;
    if( pointer1->Use == OFF && pointer2->Use == OFF){
        return 0;
    }else if ( pointer1->Use == OFF && pointer2->Use == ON){
        return 1;
    }else if ( pointer1->Use == ON && pointer2->Use == OFF){
        return -1;
    }else if( pointer1->Body->dt > pointer2->Body->dt){
        return 1;
    }else if( pointer1->Body->dt < pointer2->Body->dt){
        return -1;
    }else{
        return 0;
    }
}
#endif

static void SortStructPsink(void){

    if(Pall.Nsink == 0)
        return;

    //sprintlmpi("Start sort star");
    qsort(PsinkElements, PsinkElementsSize,
        sizeof(StructPsink),(int(*)(const void*, const void*))PsinkKeyCmp);

    for(int i=0;i<PsinkElementsSize-1;i++)
        PsinkElements[i].Next = &(PsinkElements[i+1]);
    PsinkElements[PsinkElementsSize-1].Next = NULL;

    ReconnectSinktoBody();

    return;
}

int CountEmptyPbody(void){
    int counter = 0;
    for(StructPbodyptr Pb = PbodyElements; Pb; Pb = Pb->Next)
        counter += Pb->Use?0:1;
    return counter;
}

int CountEmptyPhydro(void){
    int counter = 0;
    for(StructPhydroptr Ph = PhydroElements; Ph; Ph = Ph->Next)
        counter += Ph->Use?0:1;
    return counter;
}

int CountEmptyPstar(void){
    int counter = 0;
    for(StructPstarptr Ps = PstarElements; Ps; Ps = Ps->Next)
        counter += Ps->Use?0:1;
    return counter;
}

int CountEmptyPsink(void){
    int counter = 0;
    for(StructPsinkptr Psk= PsinkElements; Psk; Psk = Psk->Next)
        counter += Psk->Use?0:1;
    return counter;
}

void StretchStructPbody(const int nadd){

#if (StructureStretch)
    int Nsize = MAX(nadd,NAdditionUnit);

    PbodyElements = realloc(PbodyElements,(PbodyElementsSize+Nsize)*sizeof(StructPbody));
    memset(PbodyElements+PbodyElementsSize,0,Nsize*sizeof(StructPbody));

    for(int i=0;i<(PbodyElementsSize+Nsize)-1;i++)
        PbodyElements[i].Next = &(PbodyElements[i+1]);
    PbodyElements[PbodyElementsSize+Nsize-1].Next = NULL;

    for(int i=PbodyElementsSize;i<PbodyElementsSize+Nsize;i++)
        PbodyElements[i].Use = OFF;

    PbodyElementsSize += Nsize;

    ReconnectBodytoHydro();
    ReconnectBodytoStar();
    ReconnectBodytoSink();
#else
    AddStructPbody(nadd);
    PbodyElementsSize += Nsize;
#endif

    return;
}

void StretchStructPhydro(const int nadd){

#if (StructureStretch)
    int Nsize = MAX(nadd,NAdditionUnit);

    PhydroElements = realloc(PhydroElements,(PhydroElementsSize+Nsize)*sizeof(StructPhydro));
    memset(PhydroElements+PhydroElementsSize,0,Nsize*sizeof(StructPhydro));

    for(int i=0;i<(PhydroElementsSize+Nsize)-1;i++)
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    PhydroElements[PhydroElementsSize+Nsize-1].Next = NULL;

    for(int i=PhydroElementsSize;i<PhydroElementsSize+Nsize;i++)
        PhydroElements[i].Use = OFF;

    PhydroElementsSize += Nsize;

    //dprintlmpi(PhydroElementsSize);
    //fflush(NULL);
    ReconnectHydrotoBody();
#else
    AddStructPhydro(nadd);
    PhydroElementsSize += Nsize;
#endif

    return;
}


void StretchStructPstar(const int nadd){

#if (StructureStretch)
    int Nsize = MAX(nadd,NAdditionUnit);

    PstarElements = realloc(PstarElements,(PstarElementsSize+Nsize)*sizeof(StructPstar));
    memset(PstarElements+PstarElementsSize,0,Nsize*sizeof(StructPstar));

    for(int i=0;i<(PstarElementsSize+Nsize)-1;i++)
        PstarElements[i].Next = &(PstarElements[i+1]);
    PstarElements[PstarElementsSize+Nsize-1].Next = NULL;

    for(int i=PstarElementsSize;i<PstarElementsSize+Nsize;i++)
        PstarElements[i].Use = OFF;

    PstarElementsSize += Nsize;

    ReconnectStartoBody();
#else
    AddStructPstar(nadd);
    PstarElementsSize += Nsize;
#endif

    return;
}

void StretchStructPsink(const int nadd){

#if (StructureStretch)
    int Nsize = MAX(nadd,NAdditionUnit);

    PsinkElements = realloc(PsinkElements,(PsinkElementsSize+Nsize)*sizeof(StructPsink));
    memset(PsinkElements+PsinkElementsSize,0,Nsize*sizeof(StructPsink));

    for(int i=0;i<(PsinkElementsSize+Nsize)-1;i++)
        PsinkElements[i].Next = &(PsinkElements[i+1]);
    PsinkElements[PsinkElementsSize+Nsize-1].Next = NULL;

    for(int i=PsinkElementsSize;i<PsinkElementsSize+Nsize;i++)
        PsinkElements[i].Use = OFF;

    PsinkElementsSize += Nsize;

    ReconnectSinktoBody();
#else
    AddStructPink(nadd);
    PsinkElementsSize += Nsize;
#endif

    return;
}

void GenerateStructPbody(const int Nallocate){

    int AllocationSize = MAX(Nallocate,NAdditionUnit);
    PbodySize = AllocationSize;
    PbodyElementsSize = AllocationSize;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(AllocationSize*sizeof(StructPbodyptr));

    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));

    for(int i=0;i<(AllocationSize-1);i++)
        PbodyElements[i].Next = &(PbodyElements[i+1]);
    PbodyElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++)
        Pbody[i] = PbodyElements+i;

    return;
}

void GenerateStructPhydro(const int Nallocate){

    int AllocationSize = MAX(Nallocate,NAdditionUnit);
    PhydroSize = AllocationSize;
    PhydroElementsSize = AllocationSize;

    PhydroElements = malloc(AllocationSize*sizeof(StructPhydro));
    Phydro = malloc(AllocationSize*sizeof(StructPhydroptr));

    memset(PhydroElements,0,AllocationSize*sizeof(StructPhydro));

    for(int i=0;i<(AllocationSize-1);i++)
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++)
        Phydro[i] = PhydroElements+i;

    return;
}

void GenerateStructPstar(const int Nallocate){

    int AllocationSize = MAX(Nallocate,NAdditionUnit);
    PstarSize = AllocationSize;
    PstarElementsSize = AllocationSize;
    
    PstarElements = malloc(AllocationSize*sizeof(StructPstar));
    Pstar = malloc(AllocationSize*sizeof(StructPstarptr));

    memset(PstarElements,0,AllocationSize*sizeof(StructPstar));

    for(int i=0;i<(AllocationSize-1);i++)
        PstarElements[i].Next = &(PstarElements[i+1]);
    PstarElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++)
        Pstar[i] = PstarElements+i;

    return;
}

void GenerateStructPsink(const int Nallocate){

    int AllocationSize = MAX(Nallocate,NAdditionUnit);
    PsinkSize = AllocationSize;
    PsinkElementsSize = AllocationSize;
    
    PsinkElements = malloc(AllocationSize*sizeof(StructPsink));
    Psink = malloc(AllocationSize*sizeof(StructPsinkptr));

    memset(PsinkElements,0,AllocationSize*sizeof(StructPsink));

    for(int i=0;i<(AllocationSize-1);i++)
        PsinkElements[i].Next = &(PsinkElements[i+1]);
    PsinkElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++)
        Psink[i] = PsinkElements+i;

    return;
}

void FreeStructPbody(void){

    PbodySize = PbodyElementsSize = 0;

    free(PbodyElements);
    free(Pbody);

    return;
}

void FreeStructPhydro(void){

    PhydroSize = PhydroElementsSize = 0;

    free(PhydroElements);
    free(Phydro);

    return;
}

void FreeStructPstar(void){

    PstarSize = PstarElementsSize = 0;
    
    free(PstarElements);
    free(Pstar);

    return;
}

void FreeStructPsink(void){

    PsinkSize = PsinkElementsSize = 0;
    
    free(PsinkElements);
    free(Psink);

    return;
}

StructPbodyptr ReturnEmptyBodyStructurePointer(void){

    if(PbodySize == 0){
        // sprintlmpi("Add New Structure!");
        GenerateStructPbody(1);
    }

    for(StructPbodyptr Pb = PbodyElements; Pb; Pb = Pb->Next){
        if(Pb->Use == OFF){
            //pprintlmpi(Pb);
            return Pb;
        }
    }

#if (StructureStretch)
    StretchStructPbody(1);
#else
    AddStructPbody(1);
#endif
    for(StructPbodyptr Pb = PbodyElements; Pb; Pb = Pb->Next){
        if(Pb->Use == OFF)
            return Pb;
    }
    return NULL;
}

StructPhydroptr ReturnEmptyHydroStructurePointer(void){

    // search unused hydro structure
    if(PhydroSize == 0){
        //sprintlmpi("Add New Structure!");
        GenerateStructPhydro(1);
    }

    for(StructPhydroptr Ph = PhydroElements; Ph; Ph = Ph->Next){
        if(Ph->Use == OFF){
            return Ph;
        } 
    }

#if (StructureStretch)
    StretchStructPhydro(1);
#else
    AddStructPhydro(1);
#endif
    for(StructPhydroptr Ph = PhydroElements; Ph; Ph = Ph->Next){
        if(Ph->Use == OFF)
            return Ph;
    }
    return NULL;
}

StructPstarptr ReturnEmptyStarStructurePointer(void){

    // search unused star structure
    if(PstarSize == 0){
        //sprintlmpi("Add New Structure!");
        GenerateStructPstar(1);
    }

    for(StructPstarptr Ps = PstarElements; Ps; Ps = Ps->Next){
        if(Ps->Use == OFF)
            return Ps;
    }

#if (StructureStretch)
    StretchStructPstar(1);
#else
    AddStructPstar(1);
#endif
    for(StructPstarptr Ps = PstarElements; Ps; Ps = Ps->Next){
        if(Ps->Use == OFF)
            return Ps;
    }
    return NULL;
}

StructPsinkptr ReturnEmptySinkStructurePointer(void){

    // search unused star structure
    if(PsinkSize == 0){
        //sprintlmpi("Add New Structure!");
        GenerateStructPsink(1);
    }

    for(StructPsinkptr Psk = PsinkElements; Psk; Psk = Psk->Next){
        if(Psk->Use == OFF)
            return Psk;
    }

#if (StructureStretch)
    StretchStructPsink(1);
#else
    AddStructPsink(1);
#endif
    for(StructPsinkptr Psk = PsinkElements; Psk; Psk = Psk->Next){
        if(Psk->Use == OFF)
            return Psk;
    }
    return NULL;
}

StructPbodyptr ReturnEmptyBodyStructurePointerWithCounter(int *IndexStart){

    if(PbodySize == 0){
        // sprintlmpi("Add New Structure!");
        GenerateStructPbody(1);
    }

    int counter = *IndexStart;
    for(StructPbodyptr Pb = PbodyElements + *IndexStart; Pb; Pb = Pb->Next){
        if(Pb->Use == OFF){
            *IndexStart = counter;
            return Pb;
        }
        counter ++;
    }
    *IndexStart = 0;

#if (StructureStretch)
    StretchStructPbody(1);
#else
    AddStructPbody(1);
#endif
    for(StructPbodyptr Pb = PbodyElements; Pb; Pb = Pb->Next){
        if(Pb->Use == OFF)
            return Pb;
    }
    return NULL;
}

StructPhydroptr ReturnEmptyHydroStructurePointerWithCounter(int *IndexStart){

    // search unused hydro structure
    if(PhydroSize == 0){
        //sprintlmpi("Add New Structure!");
        GenerateStructPhydro(1);
    }


    int counter = *IndexStart;
    for(StructPhydroptr Ph = PhydroElements + *IndexStart; Ph; Ph = Ph->Next){
        if(Ph->Use == OFF){
            *IndexStart = counter;
            return Ph;
        } 
        counter ++;
    }
    *IndexStart = 0;

#if (StructureStretch)
    StretchStructPhydro(1);
#else
    AddStructPhydro(1);
#endif
    for(StructPhydroptr Ph = PhydroElements; Ph; Ph = Ph->Next){
        if(Ph->Use == OFF)
            return Ph;
    }
    return NULL;
}

StructPstarptr ReturnEmptyStarStructurePointerWithCounter(int *IndexStart){

    // search unused star structure
    if(PstarSize == 0){
        //sprintlmpi("Add New Structure!");
        GenerateStructPstar(1);
    }

    int counter = *IndexStart;
    for(StructPstarptr Ps = PstarElements + *IndexStart; Ps; Ps = Ps->Next){
        if(Ps->Use == OFF){
            *IndexStart = counter;
            return Ps;
        } 
        counter ++;
    }
    *IndexStart = 0;

#if (StructureStretch)
    StretchStructPstar(1);
#else
    AddStructPstar(1);
#endif
    for(StructPstarptr Ps = PstarElements; Ps; Ps = Ps->Next){
        if(Ps->Use == OFF)
            return Ps;
    }
    return NULL;
}

StructPsinkptr ReturnEmptySinkStructurePointerWithCounter(int *IndexStart){

    // search unused star structure
    if(PsinkSize == 0){
        //sprintlmpi("Add New Structure!");
        GenerateStructPsink(1);
    }

    int counter = *IndexStart;
    for(StructPsinkptr Psk = PsinkElements + *IndexStart; Psk; Psk = Psk->Next){
        if(Psk->Use == OFF){
            *IndexStart = counter;
            return Psk;
        }
        counter ++;
    }
    *IndexStart = 0;

#if (StructureStretch)
    StretchStructPsink(1);
#else
    AddStructPsink(1);
#endif
    for(StructPsinkptr Psk = PsinkElements; Psk; Psk = Psk->Next){
        if(Psk->Use == OFF)
            return Psk;
    }
    return NULL;
}

void ReConnectPointers(void){

    // for all particles 
    if(Pall.Ntotal > PbodySize){
        PbodySize = (int)(ForAngelsShare*Pall.Ntotal);
        free(Pbody);
        Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
        //dprintlmpi(PbodySize);
    }
    int counter = 0;
    for(StructPbodyptr Pb = PbodyElements; Pb; Pb = Pb->Next){
        if(Pb->Use == ON){
            Pbody[counter] = Pb;
            counter ++;
        }
    }
    if(counter != Pall.Ntotal){
        //MPI_Finalize();
        fprintf(stderr,"[%03d] Element Number of Body is incorrect!\n",MPIGetMyID());
        fprintf(stderr,"MyID = %d, conuter = %d, Pall.Ntotal = %ld\n",MPIGetMyID(),counter,Pall.Ntotal);
        MPI_Abort(MPI_COMM_WORLD,ParticleDataArrayConsistencyError);
        exit(DecompositionElementNumberError);
    }

    // for hydro particles 
    if(Pall.Nhydro > PhydroSize){
        PhydroSize = (int)(ForAngelsShare*Pall.Nhydro);
        free(Phydro);
        Phydro = malloc(PhydroSize*sizeof(StructPhydroptr));
        //dprintlmpi(PhydroSize);
    }
    counter = 0;
    int ncounter = 0;
    for(StructPhydroptr Ph = PhydroElements; Ph; Ph = Ph->Next){
        if(Ph->Use == ON){
            Phydro[counter] = Ph;
            counter ++;
        } else {
            ncounter ++;
        }
    }
    //dprintlmpi(ncounter);
    if(counter != Pall.Nhydro){
        fprintf(stderr,"[%03d] Element Number of Hydro is incorrect!\n",MPIGetMyID());
        fprintf(stderr,"MyID = %d, conuter = %d, Pall.Nhydro = %ld\n",MPIGetMyID(),counter,Pall.Nhydro);
        MPI_Abort(MPI_COMM_WORLD,ParticleDataArrayConsistencyError);
        exit(DecompositionElementNumberError);
    }

    // for star particles 
    if(Pall.Nstars > PstarSize){
        PstarSize = (int)(ForAngelsShare*Pall.Nstars);
        free(Pstar);
        Pstar = malloc(PstarSize*sizeof(StructPstarptr));
        //dprintlmpi(PstarSize);
    }
    counter = 0;
    for(StructPstarptr Ps = PstarElements; Ps; Ps = Ps->Next){
        if(Ps->Use == ON){
            Pstar[counter] = Ps;
            counter ++;
        }
    }
    if(counter != Pall.Nstars){
        //MPI_Finalize();
        fprintf(stderr,"[%03d] Element Number of Star is incorrect!\n",MPIGetMyID());
        fprintf(stderr,"MyID = %d, conuter = %d, Pall.Nstars = %ld\n",MPIGetMyID(),counter,Pall.Nstars);
        MPI_Abort(MPI_COMM_WORLD,ParticleDataArrayConsistencyError);
        exit(DecompositionElementNumberError);
    }

    // for sink particles 
    if(Pall.Nsink > PsinkSize){
        PsinkSize = (int)(MAX(ForAngelsShare*Pall.Nsink,NAdditionUnit));
        free(Psink);
        Psink = malloc(PsinkSize*sizeof(StructPsinkptr));
    }
    counter = 0;
    for(StructPsinkptr Psk = PsinkElements; Psk; Psk = Psk->Next){
        if(Psk->Use == ON){
            Psink[counter] = Psk;
            counter ++;
        }
    }
    if(counter != Pall.Nsink){
        //MPI_Finalize();
        fprintf(stderr,"[%03d] Element Number of Sink is incorrect!\n",MPIGetMyID());
        fprintf(stderr,"MyID = %d, conuter = %d, Pall.Nsink = %ld\n",MPIGetMyID(),counter,Pall.Nsink);
        assert(counter == Pall.Nsink);
        MPI_Abort(MPI_COMM_WORLD,ParticleDataArrayConsistencyError);
        exit(DecompositionElementNumberError);
    }
    return;
}

void ReConnectPointersForNSM(void){

    // for all particles 
    if(Pall.Ntotal > PbodySize){
        PbodySize = (int)(ForAngelsShare*Pall.Ntotal);
        free(Pbody);
        Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
        //dprintlmpi(PbodySize);
    }
    int counter = 0;
    for(StructPbodyptr Pb = PbodyElements; Pb; Pb = Pb->Next){
        if(Pb->Use == ON){
            Pbody[counter] = Pb;
            counter ++;
        }
    }
    if(counter != Pall.Ntotal){
        //MPI_Finalize();
        fprintf(stderr,"[%03d] Element Number of Body is incorrect!\n",MPIGetMyID());
        fprintf(stderr,"MyID = %d, conuter = %d, Pall.Ntotal = %ld\n",MPIGetMyID(),counter,Pall.Ntotal);
        MPI_Abort(MPI_COMM_WORLD,ParticleDataArrayConsistencyError);
        exit(DecompositionElementNumberError);
    }

    // for hydro particles 
    if(Pall.Nhydro > PhydroSize){
        PhydroSize = (int)(ForAngelsShare*Pall.Nhydro);
        free(Phydro);
        Phydro = malloc(PhydroSize*sizeof(StructPhydroptr));
        //dprintlmpi(PhydroSize);
    }
    counter = 0;
    int ncounter = 0;
    for(StructPhydroptr Ph = PhydroElements; Ph; Ph = Ph->Next){
        if(Ph->Use == ON){
            Phydro[counter] = Ph;
            counter ++;
        } else {
            ncounter ++;
        }
    }
    //dprintlmpi(ncounter);
    if(counter != Pall.Nhydro){
        fprintf(stderr,"[%03d] Element Number of Hydro is incorrect!\n",MPIGetMyID());
        fprintf(stderr,"MyID = %d, conuter = %d, Pall.Nhydro = %ld\n",MPIGetMyID(),counter,Pall.Nhydro);
        MPI_Abort(MPI_COMM_WORLD,ParticleDataArrayConsistencyError);
        exit(DecompositionElementNumberError);
    }

    // for star particles 
    if(Pall.Nstars > PstarSize){
        PstarSize = (int)(ForAngelsShare*Pall.Nstars);
        free(Pstar);
        Pstar = malloc(PstarSize*sizeof(StructPstarptr));
        //dprintlmpi(PstarSize);
    }
    counter = 0;
    for(StructPstarptr Ps = PstarElements; Ps; Ps = Ps->Next){
        if(Ps->Use == ON){
            Pstar[counter] = Ps;
            counter ++;
        }
    }
    if(counter != Pall.Nstars){
        //MPI_Finalize();
        fprintf(stderr,"[%03d] Element Number of Star is incorrect!\n",MPIGetMyID());
        fprintf(stderr,"MyID = %d, conuter = %d, Pall.Nstars = %ld\n",MPIGetMyID(),counter,Pall.Nstars);
        MPI_Abort(MPI_COMM_WORLD,ParticleDataArrayConsistencyError);
        exit(DecompositionElementNumberError);
    }

    // for sink particles 
    if(Pall.Nsink > PsinkSize){
        PsinkSize = (int)(MAX(ForAngelsShare*Pall.Nsink,NAdditionUnit));
        free(Psink);
        Psink = malloc(PsinkSize*sizeof(StructPsinkptr));
    }
    counter = 0;
    for(StructPsinkptr Psk = PsinkElements; Psk; Psk = Psk->Next){
        if(Psk->Use == ON){
            Psink[counter] = Psk;
            counter ++;
        }
    }
    if(counter != Pall.Nsink){
        //MPI_Finalize();
        fprintf(stderr,"[%03d] Element Number of Sink is incorrect!\n",MPIGetMyID());
        fprintf(stderr,"MyID = %d, conuter = %d, Pall.Nsink = %ld\n",MPIGetMyID(),counter,Pall.Nsink);
        assert(counter == Pall.Nsink);
        MPI_Abort(MPI_COMM_WORLD,ParticleDataArrayConsistencyError);
        exit(DecompositionElementNumberError);
    }
    return;
}

void ReconnectBodytoHydro(void){

    for(int i=0;i<PbodyElementsSize;i++){
        if(PbodyElements[i].Use == ON){
            if(PbodyElements[i].Type == TypeHydro){
                StructPhydroptr Ph = PbodyElements[i].Baryon;
                Ph->Body = PbodyElements+i;
            }
        }
    }
    return;
}

void ReconnectBodytoStar(void){

    for(int i=0;i<PbodyElementsSize;i++){
        if(PbodyElements[i].Use == ON){
            if(PbodyElements[i].Type == TypeStar){
                StructPstarptr Ps = PbodyElements[i].Baryon;
                Ps->Body = PbodyElements+i;
            }
        }
    }
    return;
}

void ReconnectBodytoSink(void){

    for(int i=0;i<PbodyElementsSize;i++){
        if(PbodyElements[i].Use == ON){
            if(PbodyElements[i].Type == TypeSink){
                StructPsinkptr Psk = PbodyElements[i].Baryon;
                Psk->Body = PbodyElements+i;
            }
        }
    }
    return;
}

void ReconnectHydrotoBody(void){

    for(int i=0;i<PhydroElementsSize;i++){
        if(PhydroElements[i].Use == ON){
            StructPbodyptr Pb = PhydroElements[i].Body;
            Pb->Baryon = (void *)(PhydroElements+i);
        }
    }
    return;
}

void ReconnectStartoBody(void){

    for(int i=0;i<PstarElementsSize;i++){
        if(PstarElements[i].Use == ON){
            StructPbodyptr Pb = PstarElements[i].Body;
            Pb->Baryon = (void *)(PstarElements+i);
        }
    }
    return;
}

void ReconnectSinktoBody(void){

    for(int i=0;i<PsinkElementsSize;i++){
        if(PsinkElements[i].Use == ON){
            StructPbodyptr Pb = PsinkElements[i].Body;
            Pb->Baryon = (void *)(PsinkElements+i);
        }
    }
    return;
}

void StructureReportPbody(const int Index){

    fprintf(stderr,"=== Structure report of Pbody[%d] ===\n",Index);

    fprintf(stderr,"\n");

    fprintf(stderr,"Pointers\n");
    fprintf(stderr,"Pbody[%d]->Next = %p\n",Index,Pbody[Index]->Next);
    fprintf(stderr,"Pbody[%d]->Baryon = %p\n",Index,Pbody[Index]->Baryon);
    fprintf(stderr,"Pbody[%d]->GlobalID = %ld\n",Index,Pbody[Index]->GlobalID);

    fprintf(stderr,"\n");
    //fprintf(stderr,"Pbody[%d]->Use = %d\n",Index,Pbody[Index]->Use);
    if(Pbody[Index]->Use == true)   fprintf(stderr,"Pbody[%d]->Use = true\n",Index);
    else                            fprintf(stderr,"Pbody[%d]->Use = false\n",Index);
    if(Pbody[Index]->Active == true)    fprintf(stderr,"Pbody[%d]->Active = true\n",Index);
    else                                fprintf(stderr,"Pbody[%d]->Active = false\n",Index);

    fprintf(stderr,"\n");
    fprintf(stderr,"Pbody[%d]->Type = %d\n",Index,Pbody[Index]->Type);

    fprintf(stderr,"\n");
    fprintf(stderr,"Pbody[%d]->Pos = {%.19g, %.19g, %.19g}\n",Index,
            Pbody[Index]->Pos[0],Pbody[Index]->Pos[1],Pbody[Index]->Pos[2]);
    fprintf(stderr,"Pbody[%d]->PosP = {%.19g, %.19g, %.19g}\n",Index,
            Pbody[Index]->PosP[0],Pbody[Index]->PosP[1],Pbody[Index]->PosP[2]);
    fprintf(stderr,"Pbody[%d]->Vel = {%.19g, %.19g, %.19g}\n",Index,
            Pbody[Index]->Vel[0],Pbody[Index]->Vel[1],Pbody[Index]->Vel[2]);
    fprintf(stderr,"Pbody[%d]->Velh = {%.19g, %.19g, %.19g}\n",Index,
            Pbody[Index]->Velh[0],Pbody[Index]->Velh[1],Pbody[Index]->Velh[2]);
    fprintf(stderr,"Pbody[%d]->Acc = {%.19g, %.19g, %.19g}\n",Index,
            Pbody[Index]->Acc[0],Pbody[Index]->Acc[1],Pbody[Index]->Acc[2]);
    fprintf(stderr,"Pbody[%d]->Pot = %.19g\n",Index,Pbody[Index]->Pot);
    fprintf(stderr,"Pbody[%d]->Mass = %.19g\n",Index,Pbody[Index]->Mass);
    fprintf(stderr,"Pbody[%d]->Eps = %.19g\n",Index,Pbody[Index]->Eps);

    fprintf(stderr,"Pbody[%d]->InteractionList = %d\n",Index,Pbody[Index]->InteractionList);

    fprintf(stderr,"\n");
    fprintf(stderr,"Pbody[%d]->dt = %.19g\n",Index,Pbody[Index]->dt);
    fprintf(stderr,"Pbody[%d]->EraLocal = %.19g\n",Index,Pbody[Index]->EraLocal);

    return ;
}

void StructureReportPhydro(const int Index){

    fprintf(stderr,"=== Structure report of Phydro[%d] ===\n",Index);

    fprintf(stderr,"\n");

    fprintf(stderr,"Pointers\n");
    fprintf(stderr,"Phydro[%d]->Next = %p\n",Index,Phydro[Index]->Next);
    fprintf(stderr,"Phydro[%d]->Body = %p\n",Index,Phydro[Index]->Body);
    fprintf(stderr,"Phydro[%d]->GlobalID = %ld\n",Index,PhydroBody(Index)->GlobalID);
    fprintf(stderr,"Phydro[%d]->dt_hydro = %.19g\n",Index,Phydro[Index]->dt_hydro);

    fprintf(stderr,"\n");
    if(Phydro[Index]->Use == true)   fprintf(stderr,"Phydro[%d]->Use = true\n",Index);
    else                            fprintf(stderr,"Phydro[%d]->Use = false\n",Index);
    if(Phydro[Index]->CoolingFlag == true)  fprintf(stderr,"Phydro[%d]->CoolingFlag = true\n",Index);
    else                                    fprintf(stderr,"Phydro[%d]->CoolingFlag = false\n",Index);

    fprintf(stderr,"Phydro[%d]->Nlist = %u\n",Index,Phydro[Index]->Nlist);
    fprintf(stderr,"Phydro[%d]->NextLeaf = %lu\n",Index,Phydro[Index]->NextLeaf);

    fprintf(stderr,"Phydro[%d]->VelP = {%.19g, %.19g, %.19g}\n",Index,
            Phydro[Index]->VelP[0],Phydro[Index]->VelP[1],Phydro[Index]->VelP[2]);

    fprintf(stderr,"Phydro[%d]->Rho = %.19g\n",Index,Phydro[Index]->Rho);
    fprintf(stderr,"Phydro[%d]->RhoPred = %.19g\n",Index,Phydro[Index]->RhoPred);
    fprintf(stderr,"Phydro[%d]->Kernel = %.19g\n",Index,Phydro[Index]->Kernel);
    fprintf(stderr,"Phydro[%d]->KernelPred = %.19g\n",Index,Phydro[Index]->KernelPred);
    fprintf(stderr,"Phydro[%d]->Div = %.19g\n",Index,Phydro[Index]->Div);
    fprintf(stderr,"Phydro[%d]->Rot = {%.19g, %.19g, %.19g}\n",Index,
            Phydro[Index]->Rot[0],Phydro[Index]->Rot[1],Phydro[Index]->Rot[2]);
    fprintf(stderr,"Phydro[%d]->F = %.19g\n",Index,Phydro[Index]->F);
    fprintf(stderr,"Phydro[%d]->HydroAcc = {%.19g, %.19g, %.19g}\n",Index,
            Phydro[Index]->HydroAcc[0],Phydro[Index]->HydroAcc[1],Phydro[Index]->HydroAcc[2]);
    fprintf(stderr,"Phydro[%d]->U = %.19g\n",Index,Phydro[Index]->U);
    fprintf(stderr,"Phydro[%d]->UPred = %.19g\n",Index,Phydro[Index]->UPred);
    fprintf(stderr,"Phydro[%d]->Du = %.19g\n",Index,Phydro[Index]->Du);
    //fprintf(stderr,"Phydro[%d]->DuPrev = %.19g\n",Index,Phydro[Index]->DuPrev);
#ifdef USE_GRAD_H //{
    fprintf(stderr,"Phydro[%d]->Gradh = %.19g\n",Index,Phydro[Index]->Gradh);
#ifdef USE_GRAD_N //{
    fprintf(stderr,"Phydro[%d]->NumberDensity = %.19g\n",Index,Phydro[Index]->NumberDensity);
    fprintf(stderr,"Phydro[%d]->GradN = %.19g\n",Index,Phydro[Index]->GradN);
    fprintf(stderr,"Phydro[%d]->fij = %.19g\n",Index,Phydro[Index]->fij);
#endif // USE_GRAD_N //}
#endif // USE_GRAD_H //}
    fprintf(stderr,"Phydro[%d]->DuCooling = %.19g\n",Index,Phydro[Index]->DuCooling);
    fprintf(stderr,"Phydro[%d]->DQheat = %.19g\n",Index,Phydro[Index]->DQheat);
    fprintf(stderr,"Phydro[%d]->dMass = %.19g\n",Index,Phydro[Index]->dMass);
    fprintf(stderr,"Phydro[%d]->Z = %.19g\n",Index,Phydro[Index]->Z);
    fprintf(stderr,"Phydro[%d]->ZII = %.19g\n",Index,Phydro[Index]->ZII);
    fprintf(stderr,"Phydro[%d]->ZIa = %.19g\n",Index,Phydro[Index]->ZIa);
    fprintf(stderr,"Phydro[%d]->dZII = %.19g\n",Index,Phydro[Index]->dZII);
    fprintf(stderr,"Phydro[%d]->dZIa = %.19g\n",Index,Phydro[Index]->dZIa);
    fprintf(stderr,"Phydro[%d]->Vsig = %.19g\n",Index,Phydro[Index]->Vsig);
    fprintf(stderr,"\n");
#ifdef USE_VARIABLE_ALPHA
    fprintf(stderr,"Phydro[%d]->Alpha = %.19g\n",Index,Phydro[Index]->Alpha);
    fprintf(stderr,"Phydro[%d]->DAlpha = %.19g\n",Index,Phydro[Index]->DAlpha);
#endif // USE_VARIABLE_ALPHA
    fprintf(stderr,"\n");
#if (UseSFModelSpawn)
    fprintf(stderr,"Phydro[%d]->SpawnTimes = %hd\n",Index,Phydro[Index]->SpawnTimes);
    fprintf(stderr,"Phydro[%d]->SpawnMass = %.19g\n",Index,Phydro[Index]->SpawnMass);
#endif

    StructureReportPbodyElement(Phydro[Index]->Body);
    return ;
}

void StructureReportPstar(const int Index){

    fprintf(stderr,"=== Structure report of Pstar[%d] ===\n",Index);

    fprintf(stderr,"\n");

    fprintf(stderr,"Pointers\n");
    fprintf(stderr,"Pstar[%d]->Next = %p\n",Index,Pstar[Index]->Next);
    fprintf(stderr,"Pstar[%d]->Body = %p\n",Index,Pstar[Index]->Body);
    fprintf(stderr,"Pstar[%d]->GlobalID = %ld\n",Index,PstarBody(Index)->GlobalID);

    fprintf(stderr,"\n");
    if(Pstar[Index]->Use == true)   fprintf(stderr,"Pstar[%d]->Use = true\n",Index);
    else                            fprintf(stderr,"Pstar[%d]->Use = false\n",Index);
    if(Pstar[Index]->TypeII == true)   fprintf(stderr,"Pstar[%d]->TypeII = true\n",Index);
    else                            fprintf(stderr,"Pstar[%d]->TypeII = false\n",Index);
    if(Pstar[Index]->TypeIa == true)   fprintf(stderr,"Pstar[%d]->TypeIa = true\n",Index);
    else                            fprintf(stderr,"Pstar[%d]->TypeIa = false\n",Index);

    fprintf(stderr,"\n");
    fprintf(stderr,"Pstar[%d]->IMFTYPE = %hd\n",Index,Pstar[Index]->IMFTYPE);
    fprintf(stderr,"Pstar[%d]->NthChildren = %hd\n",Index,Pstar[Index]->NthChildren);
    fprintf(stderr,"Pstar[%d]->ParentGlobalID = %lu\n",Index,Pstar[Index]->ParentGlobalID);

    fprintf(stderr,"\n");
    fprintf(stderr,"Pstar[%d]->InitialMass = %.19g\n",Index,Pstar[Index]->InitialMass);
    fprintf(stderr,"Pstar[%d]->Mass = %.19g\n",Index,Pstar[Index]->Mass);
    fprintf(stderr,"Pstar[%d]->FormationTime = %.19g\n",Index,Pstar[Index]->FormationTime);
    fprintf(stderr,"Pstar[%d]->Z = %.19g\n",Index,Pstar[Index]->Z);
    fprintf(stderr,"Pstar[%d]->ZII = %.19g\n",Index,Pstar[Index]->ZII);
    fprintf(stderr,"Pstar[%d]->ZIa = %.19g\n",Index,Pstar[Index]->ZIa);
    fprintf(stderr,"Pstar[%d]->TempForm = %.19g\n",Index,Pstar[Index]->TempForm);
    fprintf(stderr,"Pstar[%d]->RhoForm = %.19g\n",Index,Pstar[Index]->RhoForm);

#ifdef USE_CELIB //{
    fprintf(stderr,"Pstar[%d]->SNIICount = %d\n",Index,Pstar[Index]->SNIICount);
    fprintf(stderr,"Pstar[%d]->EventTimeSNII = %.19g\n",Index,Pstar[Index]->EventTimeSNII);
    fprintf(stderr,"Pstar[%d]->SNIaCount = %d\n",Index,Pstar[Index]->SNIaCount);
    fprintf(stderr,"Pstar[%d]->EventTimeSNIa = %.19g\n",Index,Pstar[Index]->EventTimeSNIa);
#ifdef USE_CELIB_AGB //{
    fprintf(stderr,"Pstar[%d]->AGBCount = %d\n",Index,Pstar[Index]->AGBCount);
    fprintf(stderr,"Pstar[%d]->EventTimeAGB = %.19g\n",Index,Pstar[Index]->EventTimeAGB);
#ifdef USE_CELIB_NSM //{
    fprintf(stderr,"Pstar[%d]->NSMCount = %d\n",Index,Pstar[Index]->NSMCount);
    fprintf(stderr,"Pstar[%d]->EventTimeNSM = %.19g\n",Index,Pstar[Index]->EventTimeNSM);
#endif // USE_CELIB_NSM //}
#endif // USE_CELIB_AGB //}
    fprintf(stderr,"Pstar[%d]->Element[H]  = %.19g\n",Index,Pstar[Index]->Elements[CELibYield_H]);
    fprintf(stderr,"Pstar[%d]->Element[He] = %.19g\n",Index,Pstar[Index]->Elements[CELibYield_He]);
    fprintf(stderr,"Pstar[%d]->Element[C]  = %.19g\n",Index,Pstar[Index]->Elements[CELibYield_C]);
    fprintf(stderr,"Pstar[%d]->Element[N]  = %.19g\n",Index,Pstar[Index]->Elements[CELibYield_N]);
    fprintf(stderr,"Pstar[%d]->Element[O]  = %.19g\n",Index,Pstar[Index]->Elements[CELibYield_O]);
    fprintf(stderr,"Pstar[%d]->Element[Ne] = %.19g\n",Index,Pstar[Index]->Elements[CELibYield_Ne]);
    fprintf(stderr,"Pstar[%d]->Element[Mg] = %.19g\n",Index,Pstar[Index]->Elements[CELibYield_Mg]);
    fprintf(stderr,"Pstar[%d]->Element[Si] = %.19g\n",Index,Pstar[Index]->Elements[CELibYield_Si]);
    fprintf(stderr,"Pstar[%d]->Element[S]  = %.19g\n",Index,Pstar[Index]->Elements[CELibYield_S]);
    fprintf(stderr,"Pstar[%d]->Element[Ca] = %.19g\n",Index,Pstar[Index]->Elements[CELibYield_Ca]);
    fprintf(stderr,"Pstar[%d]->Element[Fe] = %.19g\n",Index,Pstar[Index]->Elements[CELibYield_Fe]);
    fprintf(stderr,"Pstar[%d]->Element[Ni] = %.19g\n",Index,Pstar[Index]->Elements[CELibYield_Ni]);
#endif // USE_CELIB //}

    fprintf(stderr,"\n");

    StructureReportPbodyElement(Pstar[Index]->Body);
    return ;
}

void StructureReportPsink(const int Index){

    fprintf(stderr,"=== Structure report of Psink[%d] ===\n",Index);

    fprintf(stderr,"\n");

    fprintf(stderr,"Pointers\n");
    fprintf(stderr,"Psink[%d]->Next = %p\n",Index,Psink[Index]->Next);
    fprintf(stderr,"Psink[%d]->Body = %p\n",Index,Psink[Index]->Body);
    fprintf(stderr,"Psink[%d]->GlobalID = %ld\n",Index,PsinkBody(Index)->GlobalID);

    fprintf(stderr,"\n");
    if(Psink[Index]->Use == true)   fprintf(stderr,"Psink[%d]->Use = true\n",Index);
    else                            fprintf(stderr,"Psink[%d]->Use = false\n",Index);

    fprintf(stderr,"\n");
    fprintf(stderr,"Psink[%d]->FormationTime = %.19g\n",Index,Psink[Index]->FormationTime);
    fprintf(stderr,"Psink[%d]->AM = %.19g, %.19g, %.19g\n",Index,
            Psink[Index]->AM[0],Psink[Index]->AM[1],Psink[Index]->AM[2]);
    fprintf(stderr,"Psink[%d]->Z = %.19g\n",Index,Psink[Index]->Z);
    fprintf(stderr,"Psink[%d]->ZII = %.19g\n",Index,Psink[Index]->ZII);
    fprintf(stderr,"Psink[%d]->ZIa = %.19g\n",Index,Psink[Index]->ZIa);

    fprintf(stderr,"\n");

    StructureReportPbodyElement(Psink[Index]->Body);
    return ;
}

static void StructureReportPbodyElement(StructPbodyptr Pb){

    fprintf(stderr,"=== Structure report of Pb at %p ===\n",Pb);

    fprintf(stderr,"\n");

    fprintf(stderr,"Pointers\n");
    fprintf(stderr,"Pb->Next = %p\n",Pb->Next);
    fprintf(stderr,"Pb->Baryon = %p\n",Pb->Baryon);
    fprintf(stderr,"Pb->GlobalID = %ld\n",Pb->GlobalID);

    fprintf(stderr,"\n");
    //fprintf(stderr,"Pbody[%d]->Use = %d\n",Index,Pbody[Index]->Use);
    if(Pb->Use == true)   fprintf(stderr,"Pb->Use = true\n");
    else                            fprintf(stderr,"Pb->Use = false\n");
    if(Pb->Active == true)    fprintf(stderr,"Pb->Active = true\n");
    else                                fprintf(stderr,"Pb->Active = false\n");

    fprintf(stderr,"\n");
    fprintf(stderr,"Pb->Type = %d\n",Pb->Type);

    fprintf(stderr,"\n");
    fprintf(stderr,"Pb->Pos = {%.19g, %.19g, %.19g}\n",Pb->Pos[0],Pb->Pos[1],Pb->Pos[2]);
    fprintf(stderr,"Pb->PosP = {%.19g, %.19g, %.19g}\n",Pb->PosP[0],Pb->PosP[1],Pb->PosP[2]);
    fprintf(stderr,"Pb->Vel = {%.19g, %.19g, %.19g}\n",Pb->Vel[0],Pb->Vel[1],Pb->Vel[2]);
    fprintf(stderr,"Pb->Velh = {%.19g, %.19g, %.19g}\n",Pb->Velh[0],Pb->Velh[1],Pb->Velh[2]);
    fprintf(stderr,"Pb->Acc = {%.19g, %.19g, %.19g}\n",Pb->Acc[0],Pb->Acc[1],Pb->Acc[2]);
    fprintf(stderr,"Pb->Pot = %.19g\n",Pb->Pot);
    fprintf(stderr,"Pb->Mass = %.19g\n",Pb->Mass);
    fprintf(stderr,"Pb->Eps = %.19g\n",Pb->Eps);

    fprintf(stderr,"Pb->InteractionList = %d\n",Pb->InteractionList);

    fprintf(stderr,"\n");
    fprintf(stderr,"Pb->dt = %.19g\n",Pb->dt);
    fprintf(stderr,"Pb->EraLocal = %.19g\n",Pb->EraLocal);

    return ;
}

static void StructureReportPhydroElement(StructPhydroptr Ph){

    fprintf(stderr,"=== Structure report of Ph at %p ===\n",Ph);

    fprintf(stderr,"\n");

    fprintf(stderr,"Pointers\n");
    fprintf(stderr,"Ph->Next = %p\n",Ph->Next);
    fprintf(stderr,"Ph->Body = %p\n",Ph->Body);

    fprintf(stderr,"\n");
    if(Ph->Use == true)   fprintf(stderr,"Ph->Use = true\n");
    else                 fprintf(stderr,"Ph->Use = false\n");
    if(Ph->CoolingFlag == true)  fprintf(stderr,"Ph->CoolingFlag = true\n");
    else                                    fprintf(stderr,"Ph->CoolingFlag = false\n");

    fprintf(stderr,"Ph->Nlist = %u\n",Ph->Nlist);
    fprintf(stderr,"Ph->NextLeaf = %lu\n",Ph->NextLeaf);

    fprintf(stderr,"Ph->VelP = {%.19g, %.19g, %.19g}\n",
            Ph->VelP[0],Ph->VelP[1],Ph->VelP[2]);

    fprintf(stderr,"Ph->Rho = %.19g\n",Ph->Rho);
    fprintf(stderr,"Ph->RhoPred = %.19g\n",Ph->RhoPred);
    fprintf(stderr,"Ph->Kernel = %.19g\n",Ph->Kernel);
    fprintf(stderr,"Ph->KernelPred = %.19g\n",Ph->KernelPred);
    fprintf(stderr,"Ph->Div = %.19g\n",Ph->Div);
    fprintf(stderr,"Ph->Rot = {%.19g, %.19g, %.19g}\n",Ph->Rot[0],Ph->Rot[1],Ph->Rot[2]);
    fprintf(stderr,"Ph->F = %.19g\n",Ph->F);
    fprintf(stderr,"Ph->HydroAcc = {%.19g, %.19g, %.19g}\n",Ph->HydroAcc[0],Ph->HydroAcc[1],Ph->HydroAcc[2]);
    fprintf(stderr,"Ph->U = %.19g\n",Ph->U);
    fprintf(stderr,"Ph->UPred = %.19g\n",Ph->UPred);
    fprintf(stderr,"Ph->Du = %.19g\n",Ph->Du);
    fprintf(stderr,"Ph->DuCooling = %.19g\n",Ph->DuCooling);
    fprintf(stderr,"Ph->DQheat = %.19g\n",Ph->DQheat);
    fprintf(stderr,"Ph->dMass = %.19g\n",Ph->dMass);
    fprintf(stderr,"Ph->Z = %.19g\n",Ph->Z);
    fprintf(stderr,"Ph->ZII = %.19g\n",Ph->ZII);
    fprintf(stderr,"Ph->ZIa = %.19g\n",Ph->ZIa);
    fprintf(stderr,"Ph->dZII = %.19g\n",Ph->dZII);
    fprintf(stderr,"Ph->dZIa = %.19g\n",Ph->dZIa);
    fprintf(stderr,"Ph->Vsig = %.19g\n",Ph->Vsig);
    fprintf(stderr,"\n");
#ifdef USE_VARIABLE_ALPHA
    fprintf(stderr,"Ph->Alpha = %.19g\n",Ph->Alpha);
    fprintf(stderr,"Ph->DAlpha = %.19g\n",Ph->DAlpha);
#endif // USE_VARIABLE_ALPHA
    fprintf(stderr,"\n");
#if (UseSFModelSpawn)
    fprintf(stderr,"Ph->SpawnTimes = %hd\n",Ph->SpawnTimes);
    fprintf(stderr,"Ph->SpawnMass = %.19g\n",Ph->SpawnMass);
#endif
    return ;
}

static void StructureReportPstarElement(StructPstarptr Ps){

    fprintf(stderr,"=== Structure report of Ps at %p ===\n",Ps);

    fprintf(stderr,"\n");

    fprintf(stderr,"Pointers\n");
    fprintf(stderr,"Ps->Next = %p\n",Ps->Next);
    fprintf(stderr,"Ps->Body = %p\n",Ps->Body);

    fprintf(stderr,"\n");
    if(Ps->Use == true)   fprintf(stderr,"Ps->Use = true\n");
    else                  fprintf(stderr,"Ps->Use = false\n");
    if(Ps->TypeII == true)  fprintf(stderr,"Ps->TypeII = true\n");
    else                    fprintf(stderr,"Ps->TypeII = false\n");
    if(Ps->TypeIa == true)  fprintf(stderr,"Ps->TypeIa = true\n");
    else                    fprintf(stderr,"Ps->TypeIa = false\n");

    fprintf(stderr,"\n");
    fprintf(stderr,"Ps->IMFTYPE = %hd\n",Ps->IMFTYPE);
    fprintf(stderr,"Ps->NthChildren = %hd\n",Ps->NthChildren);
    fprintf(stderr,"Ps->ParentGlobalID = %lu\n",Ps->ParentGlobalID);

    fprintf(stderr,"\n");
    fprintf(stderr,"Ps->InitialMass = %.19g\n",Ps->InitialMass);
    fprintf(stderr,"Ps->Mass = %.19g\n",Ps->Mass);
    fprintf(stderr,"Ps->FormationTime = %.19g\n",Ps->FormationTime);
    fprintf(stderr,"Ps->Z = %.19g\n",Ps->Z);
    fprintf(stderr,"Ps->ZII = %.19g\n",Ps->ZII);
    fprintf(stderr,"Ps->ZIa = %.19g\n",Ps->ZIa);
    fprintf(stderr,"Ps->TempForm = %.19g\n",Ps->TempForm);
    fprintf(stderr,"Ps->RhoForm = %.19g\n",Ps->RhoForm);

#ifdef USE_CELIB //{
    fprintf(stderr,"Ps->SNIICount = %d\n",Ps->SNIICount);
    fprintf(stderr,"Ps->EventTimeSNII = %.19g\n",Ps->EventTimeSNII);
    fprintf(stderr,"Ps->SNIaCount = %d\n",Ps->SNIaCount);
    fprintf(stderr,"Ps->EventTimeSNIa = %.19g\n",Ps->EventTimeSNIa);
#ifdef USE_CELIB_AGB //{
    fprintf(stderr,"Ps->AGBCount = %d\n",Ps->AGBCount);
    fprintf(stderr,"Ps->EventTimeAGB = %.19g\n",Ps->EventTimeAGB);
#ifdef USE_CELIB_NSM //{
    fprintf(stderr,"Ps->NSMCount = %d\n",Ps->NSMCount);
    fprintf(stderr,"Ps->EventTimeNSM = %.19g\n",Ps->EventTimeNSM);
#endif // USE_CELIB_NSM //}
#endif // USE_CELIB_AGB //}
    fprintf(stderr,"Ps->Element[H]  = %.19g\n",Ps->Elements[CELibYield_H]);
    fprintf(stderr,"Ps->Element[He] = %.19g\n",Ps->Elements[CELibYield_He]);
    fprintf(stderr,"Ps->Element[C]  = %.19g\n",Ps->Elements[CELibYield_C]);
    fprintf(stderr,"Ps->Element[N]  = %.19g\n",Ps->Elements[CELibYield_N]);
    fprintf(stderr,"Ps->Element[O]  = %.19g\n",Ps->Elements[CELibYield_O]);
    fprintf(stderr,"Ps->Element[Ne] = %.19g\n",Ps->Elements[CELibYield_Ne]);
    fprintf(stderr,"Ps->Element[Mg] = %.19g\n",Ps->Elements[CELibYield_Mg]);
    fprintf(stderr,"Ps->Element[Si] = %.19g\n",Ps->Elements[CELibYield_Si]);
    fprintf(stderr,"Ps->Element[S]  = %.19g\n",Ps->Elements[CELibYield_S]);
    fprintf(stderr,"Ps->Element[Ca] = %.19g\n",Ps->Elements[CELibYield_Ca]);
    fprintf(stderr,"Ps->Element[Fe] = %.19g\n",Ps->Elements[CELibYield_Fe]);
    fprintf(stderr,"Ps->Element[Ni] = %.19g\n",Ps->Elements[CELibYield_Ni]);
#endif // USE_CELIB //}


    fprintf(stderr,"\n");

    return ;
}

static void StructureReportPsinkElement(StructPsinkptr Psk){

    fprintf(stderr,"=== Structure report of Psink at %p ===\n",Psk);

    fprintf(stderr,"\n");

    fprintf(stderr,"Pointers\n");
    fprintf(stderr,"Psk->Next = %p\n",Psk->Next);
    fprintf(stderr,"Psk->Body = %p\n",Psk->Body);

    fprintf(stderr,"\n");
    if(Psk->Use == true)      fprintf(stderr,"Psk->Use = true\n");
    else                      fprintf(stderr,"Psk->Use = false\n");

    fprintf(stderr,"\n");
    fprintf(stderr,"Psk->FormationTime = %.19g\n",Psk->FormationTime);
    fprintf(stderr,"Psk->AM = %.19g, %.19g, %.19g\n",Psk->AM[0],Psk->AM[1],Psk->AM[2]);
    fprintf(stderr,"Psk->Z = %.19g\n",Psk->Z);
    fprintf(stderr,"Psk->ZII = %.19g\n",Psk->ZII);
    fprintf(stderr,"Psk->ZIa = %.19g\n",Psk->ZIa);

    fprintf(stderr,"\n");

    return ;
}


void StructureSizeReport(void){

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"=== Structure size report ===\n");
        fprintf(stderr,"===                       ===\n");
        fprintf(stderr,"===  Pbody  = %04d byte   ===\n",sizeof(StructPbody));
        fprintf(stderr,"===  Phydro = %04d byte   ===\n",sizeof(StructPhydro));
        fprintf(stderr,"===  Pstar  = %04d byte   ===\n",sizeof(StructPstar));
        fprintf(stderr,"===  Psink  = %04d byte   ===\n",sizeof(StructPsink));
        fprintf(stderr,"===                       ===\n");
        fprintf(stderr,"===  ===================  ===\n");
    }
    return; 
}
