#include "config.h" 

#define KeyElementBit (8*sizeof(long long int)/3)    // is the bit length for each direction.
#define BitMax ((1<<KeyElementBit)-1)
#define BitMask 1

static long long int ReturnMortonKey(const long long int keyelement);

void MakeOrderingKey(void){

    double Width,InvWidth,PosMin[3];

    if(GravityRoot.NumberofLeaves > 0){
        Width = GravityRoot.Width;
        InvWidth = 1.0/Width;
        PosMin[0] = GravityRoot.PosMin[0]; 
        PosMin[1] = GravityRoot.PosMin[1];
        PosMin[2] = GravityRoot.PosMin[2]; 
    } else if(HydroRoot.NumberofLeaves > 0) {
        Width = HydroRoot.Width;
        InvWidth = 1.0/Width;
        PosMin[0] = HydroRoot.PosMin[0]; 
        PosMin[1] = HydroRoot.PosMin[1];
        PosMin[2] = HydroRoot.PosMin[2]; 
    } else {
        // self check!
        Width = InvWidth = 0.e0;
        PosMin[0] = Pbody[0]->PosP[0];
        PosMin[1] = Pbody[0]->PosP[1];
        PosMin[2] = Pbody[0]->PosP[2];

        double PosMax[3] = {Pbody[0]->PosP[0],Pbody[0]->PosP[1],Pbody[0]->PosP[2]};
        for(int i=1;i<Pall.Ntotal;i++){
            PosMin[0] = fmin(PosMin[0],Pbody[i]->PosP[0]);
            PosMin[1] = fmin(PosMin[1],Pbody[i]->PosP[1]);
            PosMin[2] = fmin(PosMin[2],Pbody[i]->PosP[2]);

            PosMax[0] = fmax(PosMax[0],Pbody[i]->PosP[0]);
            PosMax[1] = fmax(PosMax[1],Pbody[i]->PosP[1]);
            PosMax[2] = fmax(PosMax[2],Pbody[i]->PosP[2]);
        }
        Width = fmax(fmax(PosMax[0]-PosMin[0],PosMax[1]-PosMin[1]),PosMax[2]-PosMin[2]);
        InvWidth = 1.e0/Width;

        //fprintf(stderr,"Any root of tree is allocated.\n");
        //exit(ROOT_ALLOCATION_ERROR);
    }

    for(int i=0;i<Pall.Ntotal;i++){
        long long int keyx = (long long int)(BitMax*(Pbody[i]->PosP[0]-PosMin[0])*InvWidth);
        long long int keyy = (long long int)(BitMax*(Pbody[i]->PosP[1]-PosMin[1])*InvWidth);
        long long int keyz = (long long int)(BitMax*(Pbody[i]->PosP[2]-PosMin[2])*InvWidth);
        Pbody[i]->OrderingKey = (ReturnMortonKey(keyx)|(ReturnMortonKey(keyy)<<1)|(ReturnMortonKey(keyz)<<2));
        //dlprintlmpi(Pbody[i]->OrderingKey);
    }

    return ;
}

static long long int ReturnMortonKey(const long long int keyelement){
    const static int NumberofKeyElement = KeyElementBit;
    long long int key = 0;
    for(int i=0;i<NumberofKeyElement;i++)
        key |= (((keyelement>>i)&BitMask)<<(3*i));

    return key;
}

