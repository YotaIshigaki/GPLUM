#include "config.h"
#include "PlantGravityTree.h"
#include "NeighborSearch.h"

#define NFOFCrit    (100)

static double LinkingLength;
static double LinkingLength2;
static void MakeFOFGroups(const int index, const int nlist, const int Neighbors[static MaxNeighborSize]);
static void InitializeStructFOF(void);

static void AllocateStructFOFCatalog(const int NFOFGroups);
static void ReleaseStructFOFCatalog(void);
static void MakeFOFCatalog(void);
static void PickUpHaloData(const int index, const int FOFindex);
static void SortFOFCatalogByMass(void);

void AllocateStructFOF(void){

    FOF = malloc(sizeof(struct StructFOF)*Pall.Ntotal);
    if(FOF == NULL){
        fprintf(stderr,"Allcoation of FOF is a failture.\n");
        fprintf(stderr,"In line %d, the function name is %s, the file name is %s\n",
                __LINE__,__FUNCTION__,__FILE__);
        exit(EXIT_FAILURE);
    }
    FOFSize = Pall.Ntotal;
    return ;
}

void ReleaseStructFOF(void){

    if(FOFSize > 0)
        free(FOF);
    return ;
}

void SetFOFLinkingLength(const double LL){
    LinkingLength = LL;
    LinkingLength2 = SQ(LL);
    return;
}

static void InitializeStructFOF(void){

    for(int i=0;i<Pall.Ntotal;i++){
        FOF[i].Next = NONE;
        FOF[i].Head = i;
        FOF[i].Tail = i;
    }
    return ;
}

static int Neighbors[MaxNeighborSize];
void FOFHaloFinder(void){

    PlantGravityTree();

    InitializeStructFOF();

    int RootNodeID = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        int CurentNodeID = GravityNode[RootNodeID].Children;
        int nlist;
        do {
            nlist = 0;
            CurentNodeID = GetNeighborsIterativeApproachFOF(CurentNodeID,Pbody[i]->PosP,LinkingLength,&nlist,Neighbors);
            MakeFOFGroups(i,nlist,Neighbors);

        }while(CurentNodeID != RootNodeID);

    }

    MakeFOFCatalog();

    return ;
}

static void MakeFOFGroups(const int index, const int nlist, const int Neighbors[static MaxNeighborSize]){

    for(int i=0;i<nlist;i++){
        int leaf = Neighbors[i];
        //Exclude myself
        if(index == leaf)
            continue ;

        if(FOF[leaf].Head != FOF[index].Head){
            int new_head = FOF[index].Head;
            int new_tail = FOF[leaf].Tail;
            int prev_tail = FOF[index].Tail;

            FOF[prev_tail].Next = FOF[leaf].Head;

            int current_leaf = FOF[index].Head;
            do{
                FOF[current_leaf].Head = new_head;
                FOF[current_leaf].Tail = new_tail;
                current_leaf = FOF[current_leaf].Next;
            } while(current_leaf != NONE);

        }
    }

    return ;
}

static void AllocateStructFOFCatalog(const int NFOFGroups){

    if(FOFCatalogSize > 0)
        ReleaseStructFOFCatalog();

    FOFCatalog = malloc(sizeof(struct StructFOFCatalog)*NFOFGroups);
    if(FOFCatalog == NULL){
        fprintf(stderr,"Allcoation of FOFCatalog is a failture.\n");
        fprintf(stderr,"In line %d, the function name is %s, the file name is %s\n",
                __LINE__,__FUNCTION__,__FILE__);
        exit(EXIT_FAILURE);
    }
    FOFCatalogSize = NFOFGroups;
    return;
}

static void ReleaseStructFOFCatalog(void){

    if(FOFCatalogSize > 0)
        free(FOFCatalog);
    return ;
}

static void MakeFOFCatalog(void){

    NFOFGroups = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        if(i == FOF[i].Head){
            NFOFGroups ++;
        }
    }
    dprintlmpi(NFOFGroups);
    AllocateStructFOFCatalog(NFOFGroups);

    int counter = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        if(i == FOF[i].Head){
            PickUpHaloData(i,counter);
            counter ++;
        }
    }
    SortFOFCatalogByMass();

    return ;
}

static void PickUpHaloData(const int index, const int FOFindex){

    int current_leaf = FOF[index].Head;
    FOFCatalog[FOFindex].Head = current_leaf;

    FOFCatalog[FOFindex].Number = 0;
    double Mass = 0.0;
    double Position[3]={0.0,0.0,0.0};
    double Velocity[3]={0.0,0.0,0.0};
    do{
        Mass += Pbody[current_leaf]->Mass;
        Position[0] += Pbody[current_leaf]->Mass*Pbody[current_leaf]->PosP[0];
        Position[1] += Pbody[current_leaf]->Mass*Pbody[current_leaf]->PosP[1];
        Position[2] += Pbody[current_leaf]->Mass*Pbody[current_leaf]->PosP[2];

        Velocity[0] += Pbody[current_leaf]->Mass*Pbody[current_leaf]->Vel[0];
        Velocity[1] += Pbody[current_leaf]->Mass*Pbody[current_leaf]->Vel[1];
        Velocity[2] += Pbody[current_leaf]->Mass*Pbody[current_leaf]->Vel[2];
    
        FOFCatalog[FOFindex].Number ++;
        current_leaf = FOF[current_leaf].Next;
    }while(current_leaf != NONE);

    double iMass = 1.0/Mass; 
    for(int k=0;k<3;k++){
        Position[k] *= iMass;
        Velocity[k] *= iMass;
    }

    FOFCatalog[FOFindex].Mass = Mass;

    FOFCatalog[FOFindex].Pos[0] = Position[0];
    FOFCatalog[FOFindex].Pos[1] = Position[1];
    FOFCatalog[FOFindex].Pos[2] = Position[2];

    FOFCatalog[FOFindex].Vel[0] = Velocity[0];
    FOFCatalog[FOFindex].Vel[1] = Velocity[1];
    FOFCatalog[FOFindex].Vel[2] = Velocity[2];

    return ;
}

static int FOFCatalogMassCmp(const void *x, const void *y){
    const struct StructFOFCatalog *pointer1 = (struct StructFOFCatalog *)x;
    const struct StructFOFCatalog *pointer2 = (struct StructFOFCatalog *)y;
    if(pointer1->Mass > pointer2->Mass){
        return -1;
    } else if(pointer1->Mass < pointer2->Mass){
        return 1;
    } else {
        return 0;
    }
}

static void SortFOFCatalogByMass(void){

    qsort(FOFCatalog, NFOFGroups,
        sizeof(struct StructFOFCatalog),(int(*)(const void*, const void*))FOFCatalogMassCmp);
    return;
}

void ListUpFOFCatalog(void){

    fprintf(stderr,"Print out whole FOFCatalog members.\n");
    for(int i=0;i<NFOFGroups;i++){
        if(FOFCatalog[i].Number > NFOFCrit){
            fprintf(stderr,"ID = %d, Number = %d, Mass = %g, Pos = (%g ,%g,%g)\n",
                i,FOFCatalog[i].Number,FOFCatalog[i].Mass,
                FOFCatalog[i].Pos[0],FOFCatalog[i].Pos[1],FOFCatalog[i].Pos[2]);
        }
    }

    return ;
}

struct StructFOFIDSort{
    unsigned long int GlobalID;
    int Type;
    short NthChildren;
};

static int FOFIDSortCmp(const void *x, const void *y){

    const struct StructFOFIDSort *pointer1 = (struct StructFOFIDSort *)x;
    const struct StructFOFIDSort *pointer2 = (struct StructFOFIDSort *)y;
    if(pointer1->GlobalID > pointer2->GlobalID){
        return 1;
    } else if(pointer1->GlobalID < pointer2->GlobalID){
        return -1;
    } else {
        if(pointer1->Type > pointer2->Type){
            return 1;
        } else if(pointer1->Type < pointer2->Type){
            return -1;
        } else {
            if(pointer1->NthChildren > pointer2->NthChildren){
                return 1;
            } else if(pointer1->NthChildren < pointer2->NthChildren){
                return -1;
            } else {
                return 0;
            }
        }
    }
}

struct StructFOFHeader{
    /* Unit */
    double UnitLength;
    double UnitMass;
    double UnitTime;
    /* Unit */
    /* Time */
    double TCurrent;
    double Redshift;
    /* Time */
    /* Cosmological Parameters */
    double OmegaB;      // OmegaB = rho_b/rho_crit. 
    double OmegaCDM;    // OmegaCDM = rho_cdm/rho_crit.
    double OmegaM;      // OmegaM = OmegaB + OmegaCDM.
    double OmegaL;      // OmegaL = 1.0 - OmegaM. 
    double Hubble;      // The Hubble paramter at z=0. hubble x 100 kms/s/Mpc(in simulation unit)
    /* Cosmological Parameters */
    /* Number of particles */
    unsigned long int    Ntotal_t;
    unsigned long int    NDM_t;
    unsigned long int    Nhydro_t;
    unsigned long int    Nstars_t;
    /* Number of particles */
} FOFHeader;

static void CopyPallToFOFHeader(void){

    /* Unit */
    FOFHeader.UnitLength = Pall.UnitLength;
    FOFHeader.UnitMass = Pall.UnitMass;
    FOFHeader.UnitTime = Pall.UnitTime;
    /* Unit */

    /* Time */
    FOFHeader.TCurrent = Pall.TCurrent;
    FOFHeader.Redshift = Pall.Redshift;
    /* Time */

    /* Cosmological Parameters */
    FOFHeader.OmegaB = Pall.OmegaB;
    FOFHeader.OmegaCDM = Pall.OmegaCDM;
    FOFHeader.OmegaM = Pall.OmegaM;
    FOFHeader.OmegaL = Pall.OmegaL;
    FOFHeader.Hubble = Pall.Hubble;
    /* Cosmological Parameters */

    /* Number of particles */
    FOFHeader.Ntotal_t = Pall.Ntotal_t;
    FOFHeader.NDM_t = Pall.NDM_t;
    FOFHeader.Nhydro_t = Pall.Nhydro_t;
    FOFHeader.Nstars_t = Pall.Nstars_t;
    /* Number of particles */

    return ; 
}

char OutPutDir[] = "./data/FOF";
void WriteFOFCatalog(void){

    FILE *fp;
    char fname[MaxCharactersInLine];

    int NFOFWrite = 0;
    for(int i=0;i<NFOFGroups;i++)
        if(FOFCatalog[i].Number > NFOFCrit)
            NFOFWrite ++; 

    sprintf(fname,"%s/FOF.%04d",OutPutDir,Pall.OutPutFileNumber);
    fprintf(stderr,"%s\n",fname);

    MakeDir(OutPutDir);
    CopyPallToFOFHeader();

    FileOpen(fp,fname,"wb");

    fwrite(&FOFHeader,sizeof(struct StructFOFHeader),1,fp);

    fwrite(&NFOFWrite,sizeof(int),1,fp);
    fwrite(&LinkingLength,sizeof(double),1,fp);
    fwrite(FOFCatalog,sizeof(struct StructFOFCatalog),NFOFWrite,fp);

    struct StructFOFIDSort *FOFIDSort;
    FOFIDSort = malloc(sizeof(struct StructFOFIDSort)*FOFCatalog[0].Number);
    for(int i=0;i<NFOFWrite;i++){
        // sort based on global id
        int counter = 0;
        int current_leaf = FOFCatalog[i].Head;
        do{
            FOFIDSort[counter].GlobalID = Pbody[current_leaf]->GlobalID;
            FOFIDSort[counter].Type = Pbody[current_leaf]->Type;
            if(Pbody[current_leaf]->Type == TypeStar){
                FOFIDSort[counter].NthChildren = PbodyStar(current_leaf)->NthChildren;
            } else {
                FOFIDSort[counter].NthChildren = NONE;
            }

            counter ++;
            current_leaf = FOF[current_leaf].Next;
        }while(current_leaf != NONE);
        assert(counter == FOFCatalog[i].Number);

        qsort(FOFIDSort,FOFCatalog[i].Number,sizeof(struct StructFOFIDSort),
                (int(*)(const void*, const void*))FOFIDSortCmp);

        // write index
        fwrite(FOFIDSort,sizeof(struct StructFOFIDSort),FOFCatalog[i].Number,fp);

    }
    free(FOFIDSort);

    fclose(fp);

    return ;
}
