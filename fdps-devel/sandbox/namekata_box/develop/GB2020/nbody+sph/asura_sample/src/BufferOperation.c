#include "config.h"
#include "BufferOperation.h"

static void ReleaseBufferExportSend(void);
static void ReleaseBufferExportRecv(void);
static void ReleaseBufferImportSend(void);
static void ReleaseBufferImportRecv(void);

void InitializeCommunicationBuffers(void){

    int NProcs = MPIGetNumProcs();
    NumberofBufferExportSendAllocated = malloc(sizeof(int)*(NProcs-1)+1);
    for(int i=0;i<NProcs-1;i++)
        NumberofBufferExportSendAllocated[i] = 0;
    NumberofBufferExportRecvAllocated = 0;

    NumberofBufferImportRecvAllocated = malloc(sizeof(int)*(NProcs-1)+1);
    for(int i=0;i<NProcs-1;i++)
        NumberofBufferImportRecvAllocated[i] = 0;
    NumberofBufferImportSendAllocated = 0;

    BufferExportSend = malloc(sizeof(void *)*(NProcs-1)+1);
    BufferImportRecv = malloc(sizeof(void *)*(NProcs-1)+1);

    return;
}

void ReportAllocatedMemorySizes(void){

    int NProcs = MPIGetNumProcs();

    size_t ExportSendByte = 0;
    for(int i=0;i<NProcs-1;i++)
        ExportSendByte += NumberofBufferExportSendAllocated[i];
    size_t ExportRecvByte = NumberofBufferExportRecvAllocated;

    size_t ImportSendByte = NumberofBufferImportSendAllocated;
    size_t ImportRecvByte = 0;
    for(int i=0;i<NProcs-1;i++)
        ImportRecvByte += NumberofBufferImportRecvAllocated[i];

    //mean/min/max
    long int LocalBuffersInKB[4] = {ExportSendByte/1024,ExportRecvByte/1024,ImportSendByte/1024,ImportRecvByte/1024};
    long int MeanBuffersInKB[4],MinBuffersInKB[4],MaxBuffersInKB[4];
    MPI_Allreduce(LocalBuffersInKB,MeanBuffersInKB,4,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(LocalBuffersInKB,MinBuffersInKB,4,MPI_LONG,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(LocalBuffersInKB,MaxBuffersInKB,4,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
    MeanBuffersInKB[0] /= NProcs; MeanBuffersInKB[1] /= NProcs;
    MeanBuffersInKB[2] /= NProcs; MeanBuffersInKB[3] /= NProcs;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"\n");
        fprintf(stderr,"=== Report: allocated communication buffer ===\n");
        fprintf(stderr,"    ExportSend: Mean %.3g MB, \tMin %.3g MB, \tMax %.3g MB\n",
                MeanBuffersInKB[0]/(1024.0),MinBuffersInKB[0]/(1024.0),MaxBuffersInKB[0]/(1024.0));
        fprintf(stderr,"    ExportRecv: Mean %.3g MB, \tMin %.3g MB, \tMax %.3g MB\n",
                MeanBuffersInKB[1]/(1024.0),MinBuffersInKB[1]/(1024.0),MaxBuffersInKB[1]/(1024.0));
        fprintf(stderr,"    ImportSend: Mean %.3g MB, \tMin %.3g MB, \tMax %.3g MB\n",
                MeanBuffersInKB[2]/(1024.0),MinBuffersInKB[2]/(1024.0),MaxBuffersInKB[2]/(1024.0));
        fprintf(stderr,"    ImportRecv: Mean %.3g MB, \tMin %.3g MB, \tMax %.3g MB\n",
                MeanBuffersInKB[3]/(1024.0),MinBuffersInKB[3]/(1024.0),MaxBuffersInKB[3]/(1024.0));
        //fprintf(stderr,"    ExportSend: %.3g MB\n",(ExportSendByte/(1024.0*1024.0)));
        //fprintf(stderr,"    ExportRecv: %.3g MB\n",(ExportRecvByte/(1024.0*1024.0)));
        //fprintf(stderr,"    ImportSend: %.3g MB\n",(ImportSendByte/(1024.0*1024.0)));
        //fprintf(stderr,"    ImportRecv: %.3g MB\n",(ImportRecvByte/(1024.0*1024.0)));
        fprintf(stderr,"=== Report: allocated communication buffer ===\n");
    }
    return;
}
    
void CheckSizeofBufferExportSend(const int NumberofExportSend[restrict], const size_t SizeofElement){
    int NProcs = MPIGetNumProcs();
    for(int i=0;i<NProcs-1;i++){
        size_t BufferSizeInByte = (MAX(NumberofExportSend[i],1))*SizeofElement;
        if(NumberofBufferExportSendAllocated[i] < BufferSizeInByte){
            BufferSizeInByte *= ForAngelsShare;
            if(NumberofBufferExportSendAllocated[i] == 0){
                BufferSizeInByte = MAX(FirstAllocationSize*SizeofElement,BufferSizeInByte);
                BufferExportSend[i] = malloc(BufferSizeInByte);
            } else {
                BufferExportSend[i] = realloc(BufferExportSend[i],BufferSizeInByte);
            }
            NumberofBufferExportSendAllocated[i] = BufferSizeInByte;
        }
    }
    return;
}

void CheckSizeofBufferExportSendIndex(const int NumberofExportSend, const size_t SizeofElement, const int Index){
    size_t BufferSizeInByte = (MAX(NumberofExportSend,1))*SizeofElement;
    if(NumberofBufferExportSendAllocated[Index] < BufferSizeInByte){
        BufferSizeInByte *= ForAngelsShare;
        if(NumberofBufferExportSendAllocated[Index] == 0){
            BufferSizeInByte = MAX(FirstAllocationSize*SizeofElement,BufferSizeInByte);
            BufferExportSend[Index] = malloc(BufferSizeInByte);
        } else {
            BufferExportSend[Index] = realloc(BufferExportSend[Index],BufferSizeInByte);
        }
        NumberofBufferExportSendAllocated[Index] = BufferSizeInByte;
    }
    return;
}

void CheckSizeofBufferExportRecv(const int NumberofExportRecv, const size_t SizeofElement){
    size_t BufferSizeInByte = (MAX(NumberofExportRecv,1))*SizeofElement;
    if(NumberofBufferExportRecvAllocated < BufferSizeInByte){
        BufferSizeInByte *= ForAngelsShare;
        if(NumberofBufferExportRecvAllocated == 0){
            BufferSizeInByte = MAX(FirstAllocationSize*SizeofElement,BufferSizeInByte);
            BufferExportRecv = malloc(BufferSizeInByte);
        } else {
            BufferExportRecv = realloc(BufferExportRecv,BufferSizeInByte);
        }
        NumberofBufferExportRecvAllocated = BufferSizeInByte;
    }
    return;
}

void CheckSizeofBufferImportSend(const int NumberofImportSend, const size_t SizeofElement){
    size_t BufferSizeInByte = (MAX(NumberofImportSend,1))*SizeofElement;
    if(NumberofBufferImportSendAllocated < BufferSizeInByte){
        BufferSizeInByte *= ForAngelsShare;
        if(NumberofBufferImportSendAllocated == 0){
            BufferSizeInByte = MAX(FirstAllocationSize*SizeofElement,BufferSizeInByte);
            BufferImportSend = malloc(BufferSizeInByte);
        } else {
            BufferImportSend = realloc(BufferImportSend,BufferSizeInByte);
        }
        NumberofBufferImportSendAllocated = BufferSizeInByte;
    }
    return;
}

void CheckSizeofBufferImportRecv(const int NumberofImportRecv[restrict], const size_t SizeofElement){
    int NProcs = MPIGetNumProcs();
    for(int i=0;i<NProcs-1;i++){
        size_t BufferSizeInByte = (MAX(NumberofImportRecv[i],1))*SizeofElement;
        if(NumberofBufferImportRecvAllocated[i] < BufferSizeInByte){
            BufferSizeInByte *= ForAngelsShare;
            if(NumberofBufferImportRecvAllocated[i] == 0){
                BufferSizeInByte = MAX(FirstAllocationSize*SizeofElement,BufferSizeInByte);
                BufferImportRecv[i] = malloc(BufferSizeInByte);
            } else {
                BufferImportRecv[i] = realloc(BufferImportRecv[i],BufferSizeInByte);
            }
            NumberofBufferImportRecvAllocated[i] = BufferSizeInByte;
        }
    }
    return;
}

void CheckSizeofBufferImportRecvIndex(const int NumberofImportRecv, const size_t SizeofElement, const int Index){
    size_t BufferSizeInByte = (MAX(NumberofImportRecv,1))*SizeofElement;
    if(NumberofBufferImportRecvAllocated[Index] < BufferSizeInByte){
        BufferSizeInByte *= ForAngelsShare;
        if(NumberofBufferImportRecvAllocated[Index] == 0){
            BufferSizeInByte = MAX(FirstAllocationSize*SizeofElement,BufferSizeInByte);
            BufferImportRecv[Index] = malloc(BufferSizeInByte);
        } else {
            BufferImportRecv[Index] = realloc(BufferImportRecv[Index],BufferSizeInByte);
        }
        NumberofBufferImportRecvAllocated[Index] = BufferSizeInByte;
    }
    return;
}

void ReleaseAllBuffers(void){
    ReleaseBufferExportSend();
    ReleaseBufferExportRecv();
    ReleaseBufferImportSend();
    ReleaseBufferImportRecv();
    return;
}

static void ReleaseBufferExportSend(void){
    int NProcs = MPIGetNumProcs();
    for(int i=0;i<NProcs-1;i++){
        if(NumberofBufferExportSendAllocated[i] > 0){
            free(BufferExportSend[i]);
            NumberofBufferExportSendAllocated[i] = 0;
        }
    }
    return;
}

static void ReleaseBufferExportRecv(void){
    if(NumberofBufferExportRecvAllocated > 0){
        free(BufferExportRecv);
        NumberofBufferExportRecvAllocated = 0;
    }
    return;
}

static void ReleaseBufferImportSend(void){
    if(NumberofBufferImportSendAllocated > 0){
        free(BufferImportSend);
        NumberofBufferImportSendAllocated = 0;
    }
    return;
}

static void ReleaseBufferImportRecv(void){
    int NProcs = MPIGetNumProcs();
    for(int i=0;i<NProcs-1;i++){
        if(NumberofBufferImportRecvAllocated[i] > 0){
            free(BufferImportRecv[i]);
            NumberofBufferImportRecvAllocated[i] = 0;
        }
    }
    return;
} 
