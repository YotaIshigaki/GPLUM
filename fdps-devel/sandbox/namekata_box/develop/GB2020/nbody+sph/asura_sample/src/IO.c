#include "config.h"
#include "Read.h"
#include "IO.h"
#include "StructuresForIO.h"
#include "NeighborSearch.h"
#include <getopt.h>
#include <ctype.h>

static void ParallelReadWriteAllData(void);
static void MergerAllData(void);
static void MakeLockFile(void);
static bool CheckLockFile(void);
static void RemoveLockFile(void);

/*
 * This function shows the usage of ASURA if the case that unknown parameters
 * are input from the command line.
 */
static void ShowUsage(void){

    fprintf(stderr,"Usage: ASURA Ver. %d.%d.%d\n",ASURA_MAJOR_VERSION,ASURA_MINOR_VERSION,ASURA_MICRO_VERSION);
    fprintf(stderr,"\n");

    fprintf(stderr,"New run: ");
    fprintf(stderr,"./asura.out\n");
    fprintf(stderr,"\n");

    fprintf(stderr,"Restart run: ");
    fprintf(stderr,"./asura.out -r [filename]\n");
    fprintf(stderr," Note that the number of processors and processor index should be removed from [filename].\n");
    fprintf(stderr," For instance, ./asura.out -r ./File.0001 \n");
    fprintf(stderr,"\n");

    fprintf(stderr,"./asura.out -R [filename]\n");
    fprintf(stderr," The same as -r, but this mode uses the latest outputfile and also cares the existence of the lockfile.\n");
    fprintf(stderr," For instance, ./asura.out -R ./File.0001 \n");
    fprintf(stderr,"\n");

    fprintf(stderr,"./asura.out -S [filename]\n");
    fprintf(stderr," The same as -R, but this mode seeks the latest outputfile from the given ID to 10000.\n");
    fprintf(stderr," For instance, ./asura.out -S ./File.0001 \n");
    fprintf(stderr," Note that if there are ./File.0100..../File.0500 and ASURA runs with `-S ./File.0001',\n");
    fprintf(stderr,"  ASURA adotps ./File.0500 as the restart file. \n");
    fprintf(stderr,"\n");

    fprintf(stderr,"File operations: ");
    fprintf(stderr,"./asura.out [option] [filename]\n");
    fprintf(stderr," The following commands are available.\n");
    fprintf(stderr," -d:\t\t Double the restart files.\n");
    fprintf(stderr," -h:\t\t Halve the restart files.\n");
    fprintf(stderr," -m:\t\t Merge the restart files into a single restart file.\n");
    fprintf(stderr,"\n");

    fprintf(stderr,"Other switches: \n");
    fprintf(stderr," -v:\t\t Show version.\n");

    return ;
}

/*
 * Routines for Header operations. 
 */
#ifdef USE_LEAN_IO_FORMAT
static bool LeanIOCondition(void){

    // fprintf(stderr,"=== %d %d %d\n",Pall.OutPutFileNumber,LEAN_IO_FORMAT_FULL_IO_INTERVAL,
            // Pall.OutPutFileNumber%LEAN_IO_FORMAT_FULL_IO_INTERVAL);
    if(Pall.OutPutFileNumber%LEAN_IO_FORMAT_FULL_IO_INTERVAL == 0){
        return false;
    } else {
        return true;
    }
}
#endif

static void WriteHeader(FILE *fp){

    int HeaderArray[HEADER_SIZE];
    for(int i=0;i<HEADER_SIZE;i++)
        HeaderArray[i] = 0;
    
    HeaderArray[HD_SizeofHeader] = HEADER_SIZE;
    HeaderArray[HD_SizeofInt] = sizeof(int);
#ifdef COMPACT_IO_FORMAT
    HeaderArray[HD_CompactFormatFlag] = 1;
#ifdef USE_LEAN_IO_FORMAT
    if(LeanIOCondition()){
        HeaderArray[HD_LeanFormatFlag] = 1;
    } else {
        HeaderArray[HD_LeanFormatFlag] = 0;
    }
#endif // USE_LEAN_IO_FORMAT
#else // COMPACT_IO_FORMAT
    HeaderArray[HD_CompactFormatFlag] = 0;
#endif // COMPACT_IO_FORMAT

#ifdef COMPACT_DOUBLE_IO_FORMAT //{
    HeaderArray[HD_CompactDoubleFormatFlag] = 1;
    HeaderArray[HD_CompactFormatFlag] = 0;
    HeaderArray[HD_CompactMixFormatFlag] = 0;
    HeaderArray[HD_LeanFormatFlag] = 0;
#endif // COMPACT_DOUBLE_IO_FORMAT //}

#ifdef COMPACT_MIX_IO_FORMAT //{
    HeaderArray[HD_CompactMixFormatFlag] = 1;
    HeaderArray[HD_CompactDoubleFormatFlag] = 0;
    HeaderArray[HD_CompactFormatFlag] = 0;
    HeaderArray[HD_LeanFormatFlag] = 0;
#endif // COMPACT_MIX_IO_FORMAT //}

    // write residual
    fwrite(HeaderArray,sizeof(int),HEADER_SIZE,fp);
    //fprintf(stderr,"Header size = %d\n",(int)(HEADER_SIZE*sizeof(int)));

    return;
}

static void ReadHeader(FILE *fp, int HeaderArray[]){

    fread(HeaderArray,sizeof(int),HEADER_SIZE,fp);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Header information : %d %d %d %d\n",
                HeaderArray[HD_SizeofHeader],HeaderArray[HD_SizeofInt],
                HeaderArray[HD_CompactFormatFlag],HeaderArray[HD_LeanFormatFlag]);
    }
    return;
}


#ifdef USE_LEAN_IO_FORMAT
static bool CheckLeanIOFormat(char *FileName){

    FILE *fp;
    fp = fopen(FileName,"r");
    int HeaderArray[HEADER_SIZE];
    ReadHeader(fp,HeaderArray);
    fclose(fp);

    if(HeaderArray[HD_LeanFormatFlag] == 1){
        return true;
    } else {
        return false;
    }
}
#endif


/*
 * Search the latest file ID. The search start from the given ID by the command
 * line parameter. If USE_LEAN_IO_FORMAT flag is defined, check
 * HD_LeanFormatFlag in the header part and choose adequate file ID.
 *
 * This function first picks up the initial file ID from the given commannd line
 * parameters, and then searches the restart file ID from the initial file ID.
 * The detailed procedure of the restard file ID search is as follows. (1)
 * increment file ID and check the file is exist or not.  (2) If the file is not
 * exist, stop incrementation. (3) Check lockfile. If it exists, decrement the
 * file ID. (4) If USE_LEAN_IO_FORMAT is defined, decrement the file ID till the
 * HD_LeanFormatFlag  == 0. (5) Set the restat file name to FilePath[].
 */

static void LatestFilesID(char *optarg, char FilePath[]){



#if 1
    if(FindWordFromString(optarg,"temp",".")){
        strcpy(FilePath,optarg);
        return ;
    }
#else
    if(!CheckFile("notemp")){
        if(FindWordFromString(optarg,"temp",".")){
            strcpy(FilePath,optarg);
            return ;
        }
    }
#endif

    char TmpOptarg[MaxCharactersInLine];
    strcpy(TmpOptarg,optarg);

    int PositionofSeparator;
#if 0
    for(int l=strlen(TmpOptarg)-1;l!=0;l--){
        fprintf(stderr,"%c\n",TmpOptarg[l]);
        if(isdigit(TmpOptarg[l]) != 0){
            PositionofSeparator = l;
        }
    }
#endif
    for(int l=strlen(TmpOptarg)-1;l!=0;l--){
        // fprintf(stderr,"--%c\n",TmpOptarg[l]);
        if(ispunct(TmpOptarg[l]) != 0){
            PositionofSeparator = l;
            // fprintf(stderr,"+++%c\n",TmpOptarg[l]);
            break;
        }
    }
    PositionofSeparator ++;
    // dprintlmpi(PositionofSeparator);
    // fprintf(stderr,"%c\n",TmpOptarg[PositionofSeparator]);
    // fprintf(stderr," aaaa %c\n",TmpOptarg[PositionofSeparator]);

    char ID[MaxCharactersInLine];
    //strcpy(ID,TmpOptarg+PositionofSeparator+1);
    strcpy(ID,TmpOptarg+PositionofSeparator);
    // fprintf(stderr,"ID = %s %d\n",ID,atoi(ID));
   
    char Path[MaxCharactersInLine];
    TmpOptarg[PositionofSeparator] = '\0';
    strcpy(Path,TmpOptarg);
    // fprintf(stderr,"Path = %s\n",Path);

    /// search file ID. Start from atoi(ID);
    
    //char FilePath[MaxCharactersInLine];
    int counter = atoi(ID);
    while(1){
        sprintf(FilePath,"%s%04d.%03d.%03d",Path,counter,MPIGetNumProcs(),0);
        //fprintf(stderr,"%s\n",FilePath);
        if(CheckFile(FilePath) != true){
            break;
        }
        if(counter > 100000000){
            fprintf(stderr,"I cannot find restart files. %s:%s:%d\n",__FILE__,__FUNCTION__,__LINE__);
            fflush(NULL);
            MPI_Abort(MPI_COMM_WORLD,RestartFileSearchError);
        }
        counter ++;
    }
    
    if(CheckLockFile()==true){
        counter --;
    }

    // The restart file ID.
    counter --;

#ifdef USE_LEAN_IO_FORMAT
    // back from the current file ID (= counter) in order to find a restartable
    // file.
    while(1){
        char fname[MaxCharactersInLine];
        sprintf(fname,"%s%04d.%03d.%03d",Path,counter,MPIGetNumProcs(),0);
        if(CheckLeanIOFormat(fname) == false){
            break; 
        }
        counter --;
        if(counter < 0){
            fprintf(stderr,"I cannot find restart files. %s:%s:%d\n",__FILE__,__FUNCTION__,__LINE__);
            fflush(NULL);
            MPI_Abort(MPI_COMM_WORLD,RestartFileSearchError);
        }
    }
#endif // USE_LEAN_IO_FORMAT

    sprintf(FilePath,"%s%04d",Path,counter);
    fprintf(stderr,"%s\n",FilePath);
    fflush(NULL);

    return ;
}

/*
 * Search the latest file ID. Unlike the previous function, LatestFilesID, this
 * function search from the given ID to 10000, and returns the latest ID.  Even the
 * case that there are missing files between them, this function seek the latest
 * file ID. This function also cares the lockfile.
 */
static void LatestFilesIDWithSkipMode(char *optarg, char FilePath[]){

#if 1
    if(FindWordFromString(optarg,"temp",".")){
        strcpy(FilePath,optarg);
        return ;
    }
#else
    if(!CheckFile("notemp")){
        if(FindWordFromString(optarg,"temp",".")){
            strcpy(FilePath,optarg);
            return ;
        }
    }
#endif

    char TmpOptarg[MaxCharactersInLine];
    strcpy(TmpOptarg,optarg);

    int PositionofSeparator;
    for(int l=strlen(TmpOptarg)-1;l!=0;l--){
        if(ispunct(TmpOptarg[l]) != 0){
            PositionofSeparator = l;
            break;
        }
    }
    PositionofSeparator ++;

    char ID[MaxCharactersInLine];
    strcpy(ID,TmpOptarg+PositionofSeparator);
   
    char Path[MaxCharactersInLine];
    TmpOptarg[PositionofSeparator] = '\0';
    strcpy(Path,TmpOptarg);

    
    int counter = atoi(ID);
    int lastone = NONE;

    do{
        sprintf(FilePath,"%s%04d.%03d.%03d",Path,counter,MPIGetNumProcs(),0);
        if(CheckFile(FilePath) == true){
            lastone = counter;
            dprintlmpi(lastone);
        }
        counter ++;
    }while(counter < 10000);

    if(lastone == NONE){
        fprintf(stderr,"I cannot find restart files. %s:%s:%d\n",__FILE__,__FUNCTION__,__LINE__);
        fflush(NULL);
        MPI_Abort(MPI_COMM_WORLD,RestartFileSearchError);
    }

    counter = lastone;

    if(CheckLockFile()==true){
        counter --;
    }

    // The restart file ID.
    counter --;

#ifdef USE_LEAN_IO_FORMAT
    // back from the current file ID (= counter) in order to find a restartable
    // file.
    while(1){
        char fname[MaxCharactersInLine];
        sprintf(fname,"%s%04d.%03d.%03d",Path,counter,MPIGetNumProcs(),0);
        if(CheckLeanIOFormat(fname) == false){
            break; 
        }
        counter --;
        if(counter < 0){
            fprintf(stderr,"I cannot find restart files. %s:%s:%d\n",__FILE__,__FUNCTION__,__LINE__);
            fflush(NULL);
            MPI_Abort(MPI_COMM_WORLD,RestartFileSearchError);
        }
    }
#endif // USE_LEAN_IO_FORMAT

    sprintf(FilePath,"%s%04d",Path,counter);
    fprintf(stderr,"%s\n",FilePath);
    fflush(NULL);

    return ;
}


/*
 * This function returns the restartable file ID if the case with
 * USE_LEAN_IO_FORMAT and -r option to restart a simulation.  
 */
static void FindRestartableFilesID(char *optarg, char FilePath[]){

    if(FindWordFromString(optarg,"temp",".")){
        strcpy(FilePath,optarg);
        return ;
    }

    char TmpOptarg[MaxCharactersInLine];
    strcpy(TmpOptarg,optarg);

    int PositionofSeparator;
#if 0 
    for(int l=strlen(TmpOptarg)-1;l!=0;l--){
        //fprintf(stderr,"%c\n",TmpOptarg[l]);
        if(isdigit(TmpOptarg[l]) != 0){
            PositionofSeparator = l;
        }
    }
#else
    for(int l=strlen(TmpOptarg)-1;l!=0;l--){
        if(ispunct(TmpOptarg[l]) != 0){
            PositionofSeparator = l;
            break;
        }
    }
    PositionofSeparator ++;
#endif
    // dprintlmpi(PositionofSeparator);
    // fprintf(stderr,"%c\n",TmpOptarg[PositionofSeparator]);

    char ID[MaxCharactersInLine];
    // strcpy(ID,TmpOptarg+PositionofSeparator+1);
     strcpy(ID,TmpOptarg+PositionofSeparator);
    // fprintf(stderr,"ID = %s %d\n",ID,atoi(ID));
   
    char Path[MaxCharactersInLine];
    TmpOptarg[PositionofSeparator] = '\0';
    strcpy(Path,TmpOptarg);
    // fprintf(stderr,"Path = %s\n",Path);

    /// search file ID. Start from atoi(ID);
    
    //char FilePath[MaxCharactersInLine];
    int counter = atoi(ID);
    /*
    while(1){
        sprintf(FilePath,"%s%04d.%03d.%03d",Path,counter,MPIGetNumProcs(),0);
        //fprintf(stderr,"%s\n",FilePath);
        if(CheckFile(FilePath) != true){
            break;
        }
        if(counter > 100000000){
            fprintf(stderr,"I cannot find restart files. %s:%s:%d\n",__FILE__,__FUNCTION__,__LINE__);
            fflush(NULL);
            MPI_Abort(MPI_COMM_WORLD,RestartFileSearchError);
        }
        counter ++;
    }
    
    if(CheckLockFile()==true){
        counter --;
    }
    */

    // The restart file ID.
    // counter --;

#ifdef USE_LEAN_IO_FORMAT
    // back from the current file ID (= counter) in order to find a restartable
    // file.
    while(1){
        char fname[MaxCharactersInLine];
        sprintf(fname,"%s%04d.%03d.%03d",Path,counter,MPIGetNumProcs(),0);
        if(CheckLeanIOFormat(fname) == false){
            break; 
        }
        counter --;
        if(counter < 0){
            fprintf(stderr,"I cannot find restart files. %s:%s:%d\n",__FILE__,__FUNCTION__,__LINE__);
            fflush(NULL);
            MPI_Abort(MPI_COMM_WORLD,RestartFileSearchError);
        }
    }
#endif // USE_LEAN_IO_FORMAT

    sprintf(FilePath,"%s%04d",Path,counter);
    fprintf(stderr,"%s\n",FilePath);
    fflush(NULL);

    return ;
}

/*
 * This function analyzes the parameters of the run. If the cases of a new and
 * a restart runs, this function set the status to Pall.RunStatus. While the
 * cases of date file operations, this function calls other functions in order
 * to do the required operations and then terminates all processes.
 */
static int RestartID = 0;
void GetRunStatus(const int argc, char **argv){

    extern char *optarg;

    int result;
    while((result=getopt(argc,argv,"vR:r:d:h:m:S:"))!=-1){
        switch(result){
            case 'd':
                Pall.RunStatus = DataFileOperation_Double;
                strcpy(Pall.RestartFileName,optarg);
                if(MPIGetMyID() == MPI_ROOT_RANK){
                    fprintf(stderr,"Data File Operation: Doubled\n");
                    fprintf(stderr," Data File number %d ->%d\n",
                            MPIGetNumProcs(),2*MPIGetNumProcs());
                    fprintf(stderr," Restart File Name : %s\n",Pall.RestartFileName);
                    fflush(NULL);
                }
                break;

            case 'h':
                Pall.RunStatus = DataFileOperation_Halve;
                strcpy(Pall.RestartFileName,optarg);
                if(MPIGetMyID() == MPI_ROOT_RANK){
                    fprintf(stderr,"Data File Operation: Halve\n");
                    fprintf(stderr," Data File number %d ->%d\n",
                            MPIGetNumProcs(),MPIGetNumProcs()/2);
                    fprintf(stderr," Restart File Name : %s\n",Pall.RestartFileName);
                    fflush(NULL);
                }
                break;

            case 'm':
                Pall.RunStatus = DataFileOperation_Halve;
                strcpy(Pall.RestartFileName,optarg);
                if(MPIGetMyID() == MPI_ROOT_RANK){
                    fprintf(stderr,"Data File Operation: Merge\n");
                    fprintf(stderr," Data File number %d ->%d\n",
                            MPIGetNumProcs(),MPIGetNumProcs()/2);
                    fprintf(stderr," Restart File Name : %s\n",Pall.RestartFileName);
                    fflush(NULL);
                }
                break;

            case 'r':
                Pall.RunStatus = RestartSimulation;
                //strcpy(Pall.RestartFileName,optarg);
                if(MPIGetMyID() == MPI_ROOT_RANK){
                    FindRestartableFilesID(optarg,Pall.RestartFileName);
                    //LatestFilesID(optarg,Pall.RestartFileName);
                }
                MPI_Bcast(Pall.RestartFileName,MaxCharactersInLine,MPI_CHAR,MPI_ROOT_RANK,MPI_COMM_WORLD);
                if(MPIGetMyID() == MPI_ROOT_RANK){
                    fprintf(stderr,"RunStatus: Restart\n");
                    fprintf(stderr," Restart File Name : %s\n",Pall.RestartFileName);
                    fflush(NULL);
                }
                break;

            case 'R':
                Pall.RunStatus = RestartSimulation;
                if(MPIGetMyID() == MPI_ROOT_RANK){
                    LatestFilesID(optarg,Pall.RestartFileName);
                }
                MPI_Bcast(Pall.RestartFileName,MaxCharactersInLine,MPI_CHAR,MPI_ROOT_RANK,MPI_COMM_WORLD);

                if(MPIGetMyID() == MPI_ROOT_RANK){
                    fprintf(stderr,"RunStatus: Restart\n");
                    fprintf(stderr," Restart File Name : %s\n",Pall.RestartFileName);
                    fflush(NULL);
                }
                break;

            case 'S':
                Pall.RunStatus = RestartSimulation;
                if(MPIGetMyID() == MPI_ROOT_RANK){
                    LatestFilesIDWithSkipMode(optarg,Pall.RestartFileName);
                }
                MPI_Bcast(Pall.RestartFileName,MaxCharactersInLine,MPI_CHAR,MPI_ROOT_RANK,MPI_COMM_WORLD);

                if(MPIGetMyID() == MPI_ROOT_RANK){
                    fprintf(stderr,"RunStatus: Restart\n");
                    fprintf(stderr," Restart File Name : %s\n",Pall.RestartFileName);
                    fflush(NULL);
                }
                break;

            case 'v':
                if(MPIGetMyID() == MPI_ROOT_RANK){
                    fprintf(stderr,"ASURA version %d.%d.%d\n",ASURA_MAJOR_VERSION,ASURA_MINOR_VERSION,ASURA_MICRO_VERSION);
                    ShowUsage();
                    fflush(NULL);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                exit(EXIT_SUCCESS);
                break;

            case '?':
                //fprintf(stderr,"Unknown Parameter!\n");
                if(MPIGetMyID() == MPI_ROOT_RANK)
                    ShowUsage();
                fflush(NULL);
                exit(UnknownParameter);
                break;

            default:
                Pall.RunStatus = NewSimulation;
                fprintf(stderr,"RunStatus: NewSimulation\n");
                fflush(NULL);
                break;
        }
    }
    int Flag;
    MPI_Allreduce(&(Pall.RunStatus),&Flag,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(Flag >0){
        // dprintlmpi(Flag);
        Pall.RunStatus = RestartSimulation;
        MPI_Bcast(Pall.RestartFileName,MaxCharactersInLine,MPI_CHAR,MPI_ROOT_RANK,MPI_COMM_WORLD);
        fprintf(stderr,"[%02d] Restart File Name : %s\n",MPIGetMyID(),Pall.RestartFileName);
    }

    CheckIOFileName();

    if(Pall.RunStatus == DataFileOperation_Double){
        ParallelReadWriteAllData();
        fflush(NULL);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(HeaderError);
    } else if(Pall.RunStatus == DataFileOperation_Halve){
        ParallelReadWriteAllData();
        fflush(NULL);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(HeaderError);
    } else if(Pall.RunStatus == DataFileOperation_Merge){
        MergerAllData();
        fflush(NULL);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(HeaderError);
    }

    return;
}

int GetRunStatusForSingleNodeAnalysis(const int argc, char **argv){

    extern char *optarg;

    int result;
    while((result=getopt(argc,argv,"r:"))!=-1){
        switch(result){
            case 'r':
                Pall.RunStatus = RestartSimulation;
                strcpy(Pall.RestartFileName,optarg);
                fprintf(stderr,"RunStatus: Restart\n");
                fprintf(stderr," Restart File Name : %s\n",Pall.RestartFileName);
                fflush(NULL);
                break;

            case '?':
                fprintf(stderr,"Unknown Parameter!\n");
                fflush(NULL);
                exit(UnknownParameter);
                break;

            default:
                Pall.RunStatus = NewSimulation;
                fprintf(stderr,"RunStatus: NewSimulation\n");
                fflush(NULL);
                break;
        }
    }
    int Flag;
    MPI_Allreduce(&(Pall.RunStatus),&Flag,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(Flag >0){
        dprintlmpi(Flag);
        Pall.RunStatus = RestartSimulation;
        MPI_Bcast(Pall.RestartFileName,MaxCharactersInLine,MPI_CHAR,MPI_ROOT_RANK,MPI_COMM_WORLD);
        fprintf(stderr,"[%02d] Restart File Name : %s\n",MPIGetMyID(),Pall.RestartFileName);
        //MPI_Barrier(MPI_COMM_WORLD);
        //MPI_Finalize();
        //exit(0);
    }

    int NProcs = 1;
    /*
    while((result=getopt(argc,argv,"r:"))!=-1){
        switch(result){
            case 'l':
                NProcs = atoi(optarg);
                strcpy(Pall.RestartFileName,optarg);
                fprintf(stderr," Total File Number : %d\n",NProcs);
                fflush(NULL);
                break;

            default:
                NProcs = 1;
                fprintf(stderr," Total File Number : %d\n",NProcs);
                fflush(NULL);
                break;
        }
    }
    MPI_Bcast(&NProcs,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    */

    CheckIOFileName();

    return NProcs;
}

void CheckIOFileName(void){

    if(Pall.BaseFileName == NULL){
        fprintf(stderr,"Base file name is undefined.\n");
        MPI_Finalize();
        exit(BaseFileNameUnDefined);
    }

    return;
}

/*
 * Routines for lockfile operations. 
 */

static bool FistFlagForLockFile = true;
static char LockFileName[MaxCharactersInLine];
static void MakeLockFile(void){

    if(FistFlagForLockFile == true){
        if(MPIGetMyID() == MPI_ROOT_RANK){
            sprintf(LockFileName,"./lockfile.%d",LOCK_FILE_ID);
            FistFlagForLockFile = false;
        }
    }
    if(MPIGetMyID() == MPI_ROOT_RANK){
        FILE *fp;
	    fp = fopen(LockFileName,"w");
	    if(fp == NULL){
		    fprintf(stderr,"Lock file cannot open.\n");
            fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD,LockFileMakeError);
	    }
        fclose(fp);
    }
    return ;
}

static void RemoveLockFile(void){

    int RemoveStatus;
    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(CheckFile(LockFileName) == true){
            RemoveStatus = remove(LockFileName);
            if(RemoveStatus == -1){ // error handling.
                fprintf(stderr,"Lock file cannot remove.");
                fflush(stderr);
                MPI_Abort(MPI_COMM_WORLD,LockFileRemoveError);
            }
        }
    }
    return ;
}

static bool CheckLockFile(void){

    if(FistFlagForLockFile == true){
        if(MPIGetMyID() == MPI_ROOT_RANK){
            sprintf(LockFileName,"./lockfile.%d",LOCK_FILE_ID);
            FistFlagForLockFile = false;
        }
    }

    bool LockFileCheckFlag;
    if(CheckFile(LockFileName) == true){
        fprintf(stderr,"Found lock file!\n");
        RemoveLockFile();
        LockFileCheckFlag = true;
        if(MPIGetMyID() == MPI_ROOT_RANK)
            fprintf(stderr,"Found lock file!\n");
    }  else {
        LockFileCheckFlag = false;
    }
    return LockFileCheckFlag;
}

/*
 * Routines for read/write operations. 
 */

void WriteAllData(void){

    FILE *fp;
    char fname[MaxCharactersInLine];
    int DataSize = 0;
    int NProcs = MPIGetNumProcs();

    sprintf(fname,"%s",Pall.BaseFileName);

    if(MPIGetMyID()==0){
        //FileOpen(fp,fname,"w");
        //fclose(fp);
        char ClearFile[MaxCharactersInLine];
        //sprintf(ClearFile,"rm -rf %s",Pall.BaseFileName);
        //system(ClearFile);
        sprintf(ClearFile,"%s",Pall.BaseFileName);
        int result = remove(ClearFile);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for(int i=0;i<NProcs;i++){
        if(MPIGetMyID()==i){
            if(i==0){ // file open mode "w"
                FileOpen(fp,fname,"w");
                fwrite(&Pall,sizeof(struct StructPall),1,fp);
                DataSize += sizeof(struct StructPall);
            } else { // file open mode "a"
                FileOpen(fp,fname,"a");
            }

            for(int k=0;k<Pall.Ntotal;k++){
                fwrite(Pbody[k],sizeof(StructPbody),1,fp);
                DataSize += sizeof(StructPbody);
                if(Pbody[k]->Type == TypeHydro){ // type hydro
                    fwrite(Pbody[k]->Baryon,sizeof(StructPhydro),1,fp);
                    DataSize += sizeof(StructPhydro);
                } else if (Pbody[k]->Type == TypeStar){ // type star
                    fwrite(Pbody[k]->Baryon,sizeof(StructPstar),1,fp);
                    DataSize += sizeof(StructPstar);
                } else if (Pbody[k]->Type != TypeDM){
                    fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                            MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
                    MPI_Finalize();
                    exit(WriteTypeUnDefined);
                }
            }
            fclose(fp);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    int GlobalDataSize;
    MPI_Allreduce(&DataSize,&GlobalDataSize,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(MPIGetMyID() == 0){
        fprintf(stderr,"Out put Data size = %d byte\n",GlobalDataSize);
        fprintf(stderr,"Out put Data size = %d MB\n",GlobalDataSize/1024/1024);
    }

    return;
}

void ReadAllData(void){

    FILE *fp;
    char fname[MaxCharactersInLine];

    int NProcs = MPIGetNumProcs();

    sprintf(fname,"%s",Pall.RestartFileName);
    for(int i=0;i<NProcs;i++){
        if(MPIGetMyID()==i){
            FileOpen(fp,fname,"rb");
            fread(&Pall,sizeof(struct StructPall),1,fp);
            fclose(fp);
        }
    }

    Pall.Ntotal = Pall.Ntotal_t/NProcs;
    if(MPIGetMyID() == 0)
        Pall.Ntotal += Pall.Ntotal_t%NProcs;

    int Offset = (Pall.Ntotal_t/NProcs)*MPIGetMyID();
    if(MPIGetMyID() != 0)
        Offset += Pall.Ntotal_t%NProcs;

    GenerateStructPbody(Pall.Ntotal);

    for(int i=0;i<NProcs;i++){
        if(MPIGetMyID()==i){
            FileOpen(fp,fname,"rb");
            fseek(fp,sizeof(struct StructPall),SEEK_CUR);
            int CountTypes[NTypes];
            memset(CountTypes,0,sizeof(int)*NTypes);

            //count offset byte
            int OffsetByte = sizeof(struct StructPall);
            for(int k=0;k<Offset;k++){
                //fseek(fp,sizeof(StructPbody),SEEK_CUR);
                StructPbody Pb;
                fread(&Pb,sizeof(StructPbody),1,fp);
                OffsetByte += sizeof(StructPbody);
                if(Pb.Type == TypeHydro){
                    fseek(fp,sizeof(StructPhydro),SEEK_CUR);
                    OffsetByte += sizeof(StructPhydro);
                } else if (Pb.Type == TypeStar){
                    fseek(fp,sizeof(StructPstar),SEEK_CUR);
                    OffsetByte += sizeof(StructPstar);
                } else if (Pb.Type != TypeDM){
                    fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                            MPIGetMyID(),Pb.Type,__FILE__,__LINE__);
                    MPI_Finalize();
                    exit(ReadTypeUnDefined);
                }
            }

            // read only body!
            for(int k=0;k<Pall.Ntotal;k++){
                StructPbodyptr Pb = PbodyElements[k].Next;
                fread(PbodyElements+k,sizeof(StructPbody),1,fp);
                PbodyElements[k].Next = Pb;
                CountTypes[0] ++;
                if(PbodyElements[k].Type == TypeHydro){ // type hydro
                    fseek(fp,sizeof(StructPhydro),SEEK_CUR);
                    CountTypes[PbodyElements[k].Type] ++;
                } else if (PbodyElements[k].Type == TypeStar){ // type hydro
                    fseek(fp,sizeof(StructPstar),SEEK_CUR);
                    CountTypes[PbodyElements[k].Type] ++;
                } else if (PbodyElements[k].Type != TypeDM){
                    fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                            MPIGetMyID(),PbodyElements[i].Type,__FILE__,__LINE__);
                    MPI_Finalize();
                    exit(ReadTypeUnDefined);
                }
            }
            Pall.Nhydro = CountTypes[TypeHydro];
            Pall.Nstars = CountTypes[TypeStar];

            dprintlmpi(CountTypes[0]);
            dprintlmpi(CountTypes[1]);
            dprintlmpi(CountTypes[2]);

            // read extra data!
            if(CountTypes[TypeHydro]>0)
                GenerateStructPhydro(CountTypes[TypeHydro]);
            if(CountTypes[TypeStar]>0)
                GenerateStructPstar(CountTypes[TypeStar]);

            memset(CountTypes,0,sizeof(int)*NTypes);
            //fseek(fp,sizeof(struct StructPall),SEEK_SET);
            fseek(fp,OffsetByte,SEEK_SET);

            for(int k=0;k<Pall.Ntotal;k++){
                fseek(fp,sizeof(StructPbody),SEEK_CUR);
                if(PbodyElements[k].Type == TypeHydro){ // type hydro
                    StructPhydroptr Ph = PhydroElements[CountTypes[TypeHydro]].Next;
                    fread(PhydroElements+CountTypes[TypeHydro],sizeof(StructPhydro),1,fp);
                    PhydroElements[CountTypes[TypeHydro]].Next = Ph;
                    PhydroElements[CountTypes[TypeHydro]].Body = PbodyElements+k;
                    PbodyElements[k].Baryon = (void*)(PhydroElements+CountTypes[TypeHydro]);
                    CountTypes[TypeHydro] ++;
                } else if (PbodyElements[k].Type == TypeStar){ // type hydro
                    StructPstarptr Ps = PstarElements[CountTypes[TypeStar]].Next;
                    fread(PstarElements+CountTypes[TypeStar],sizeof(StructPstar),1,fp);
                    PstarElements[CountTypes[TypeStar]].Next = Ps;
                    PhydroElements[CountTypes[TypeStar]].Body = PbodyElements+k;
                    PbodyElements[k].Baryon = (void*)(PhydroElements+CountTypes[TypeStar]);
                    CountTypes[TypeStar] ++;
                } else if (PbodyElements[k].Type != TypeDM){
                    fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                            MPIGetMyID(),PbodyElements[i].Type,__FILE__,__LINE__);
                    MPI_Finalize();
                    exit(ReadTypeUnDefined);
                }
            }

            fclose(fp);
        }
    }
    ReConnectPointers();

    return;
}

void DataDump(void){

    static int LastOutputStep = 0;
    char FileOperation[MaxCharactersInLine];

    if(Pall.TStepTotal - LastOutputStep < 10)
        return;

    static double LastDataDump = 0;
    double CurrentTime = GetElapsedTime();
    if(CurrentTime-LastDataDump < DATA_DUMP_INTERVAL*60.e0 )
        return;
    LastDataDump = CurrentTime;

    // move back up 
    char OldFile[MaxCharactersInLine];
    char NewFile[MaxCharactersInLine];
    sprintf(OldFile,"%s.temp.%03d.%03d",Pall.BaseFileName,MPIGetNumProcs(),MPIGetMyID());
    sprintf(NewFile,"%s.temp.back.%03d.%03d",Pall.BaseFileName,MPIGetNumProcs(),MPIGetMyID());
    int result = rename(OldFile,NewFile);

    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    sprintf(FileOperation,"%s.temp.%03d.%03d",
            Pall.BaseFileName,MPIGetNumProcs(),MPIGetMyID());
    for(int i=0;i<MPIGetNumProcs();i++){
        if(i == MPIGetMyID()){
            ParallelWriteAllData(FileOperation);
            fflush(NULL);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // remove back up 
    result = remove(NewFile);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    LastOutputStep = Pall.TStepTotal;
    LastDataDump = CurrentTime;

    return;
}

void DataFullDump(void){

    static double LastDataDumpEpoch = 0;
    if(LastDataDumpEpoch == 0.e0){
        LastDataDumpEpoch = GetElapsedTime();
        return ;
    }
    double CurrentTime = GetElapsedTime();
    int write_flag = 0;
    if(CurrentTime-LastDataDumpEpoch > DATA_DUMP_INTERVAL*60.e0)
        write_flag ++;
    int global_write_flag;
    MPI_Allreduce(&write_flag,&global_write_flag,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(2*global_write_flag < MPIGetNumProcs()){
        return;
    }
    MakeLockFile();

    char OldFile[MaxCharactersInLine];
    char NewFile[MaxCharactersInLine];
    sprintf(OldFile,"%s.temp.%03d.%03d",Pall.BaseFileName,MPIGetNumProcs(),MPIGetMyID());
    sprintf(NewFile,"%s.temp.back.%03d.%03d",Pall.BaseFileName,MPIGetNumProcs(),MPIGetMyID());
    int result = rename(OldFile,NewFile);

    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    char FileOperation[MaxCharactersInLine];
    sprintf(FileOperation,"%s.temp.%03d.%03d",Pall.BaseFileName,MPIGetNumProcs(),MPIGetMyID());
#ifdef WRITE_ALLDATA_AT_ONCE
    ParallelDumpAllData(FileOperation);
    MPI_Barrier(MPI_COMM_WORLD);
#else
    for(int i=0;i<MPIGetNumProcs();i++){
        if(i == MPIGetMyID()){
            ParallelDumpAllData(FileOperation);
            fflush(NULL);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
    result = remove(NewFile);
    fflush(NULL);

    MPI_Barrier(MPI_COMM_WORLD);
    RemoveLockFile();

    LastDataDumpEpoch = CurrentTime;

#ifdef STOP_AFTER_DATA_DUMP
    FinalizeMPIEnv();
    exit(EXIT_SUCCESS);
#endif // STOP_AFTER_DATA_DUMP

    return;
}

void FileOutPutConstantInterval(void){

    MakeLockFile();

    void (*BinaryDataWriteFunctionPointer)(char fname[]);

    BinaryDataWriteFunctionPointer = ParallelWriteAllData;
//void ParallelWriteAllDataASCIIFormat(char fname[]){
    // BinaryDataWriteFunctionPointer = ParallelWriteAllDataASCIIFormat;

#ifdef COSMOLOGICAL_RUN //{
    double Toffset = CalcZtoT(Pall.InitialRedshift);
    // fprintf(stderr,"%1.15g %1.15g %d %1.15g %1.15g %1.15g %s:%d\n",
            // Pall.TCurrent,Toffset,Pall.OutPutFileNumber,Pall.OutPutInterval,
            // Pall.TCurrent-Toffset,
            // Pall.TCurrent-Toffset-0.5*Pall.dtmin,
            // __FUNCTION__,__LINE__);
    if(fmax(0.e0,Pall.TCurrent-Toffset)+0.5*Pall.dtmin >= Pall.OutPutFileNumber*Pall.OutPutInterval){
        // fprintf(stderr,"%g %g %d %g %s:%d\n",
                // Pall.TCurrent,Toffset,Pall.OutPutFileNumber,Pall.OutPutInterval,
                // __FUNCTION__,__LINE__);
#else // COSMOLOGICAL_RUN //}//{
    if(Pall.TCurrent >= Pall.OutPutFileNumber*Pall.OutPutInterval){
#endif // COSMOLOGICAL_RUN //}
        Pall.OutPutFileNumber ++;
        char fname[MaxCharactersInLine];
        sprintf(fname,"%s.%04d.%03d.%03d",Pall.BaseFileName,Pall.OutPutFileNumber,
                MPIGetNumProcs(),MPIGetMyID());
#ifdef WRITE_ALLDATA_AT_ONCE
        BinaryDataWriteFunctionPointer(fname);
        fflush(NULL);
#else   // WRITE_ALLDATA_AT_ONCE
        for(int i=0;i<MPIGetNumProcs();i++){
            if(i == MPIGetMyID()){
                BinaryDataWriteFunctionPointer(fname);
                fflush(NULL);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
#endif // WRITE_ALLDATA_AT_ONCE
#ifdef WRITE_DISPH_COMPATIBLE_MODE
        WriteDISPHAnalysisToolsCompatibleMode();
#endif //WRITE_DISPH_COMPATIBLE_MODE

#ifdef USE_ON_THE_FLY_4D2U //{
        Write4D2U();
#endif //USE_ON_THE_FLY_4D2U //}

    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    RemoveLockFile();

    /////////////////////////////////////////
#if 0
    for(int k=0;k<MPIGetNumProcs();k++){
        if(k == MPIGetMyID()){
            FILE *fp_hydro;
            char fname_hydro[MaxCharactersInLine];

            sprintf(fname_hydro,"%s.%03d.%02d.%02d.0",Pall.ASCIIFileName,Pall.OutPutFileNumber,
                    MPIGetNumProcs(),MPIGetMyID());
            FileOpen(fp_hydro,fname_hydro,"w");

            for(int i=0;i<Pall.Ntotal;i++){
                if(Pbody[i]->Type == TypeHydro){
                if(PbodyHydro(i)->Tag == 0){
                    fprintf(fp_hydro,"%ld %g %g %g %g %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                            Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                            Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                            Pbody[i]->Mass,Pall.ConvertNumberDensityToCGS*PbodyHydro(i)->Rho,
                            PbodyHydroU(i)*Pall.ConvertUtoT,PbodyHydroKernel(i),
                            PbodyHydroZ(i),PbodyHydroZII(i),PbodyHydroZIa(i));
                } 
                }
            }
            fclose(fp_hydro);
            fflush(NULL);

            sprintf(fname_hydro,"%s.%03d.%02d.%02d.1",Pall.ASCIIFileName,Pall.OutPutFileNumber,
                    MPIGetNumProcs(),MPIGetMyID());
            FileOpen(fp_hydro,fname_hydro,"w");

            for(int i=0;i<Pall.Ntotal;i++){
                if(Pbody[i]->Type == TypeHydro){
                if(PbodyHydro(i)->Tag == 1){
                    fprintf(fp_hydro,"%ld %g %g %g %g %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                            Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                            Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                            Pbody[i]->Mass,Pall.ConvertNumberDensityToCGS*PbodyHydro(i)->Rho,
                            PbodyHydroU(i)*Pall.ConvertUtoT,PbodyHydroKernel(i),
                            PbodyHydroZ(i),PbodyHydroZII(i),PbodyHydroZIa(i));
                }
                } 
            }
            fclose(fp_hydro);
            fflush(NULL);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
    /////////////////////////////////////////


#ifdef TASK_GALACTIC_CENTER
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->dZII =
        Phydro[i]->dZIa = 0.e0;
        Phydro[i]->ZII =
        Phydro[i]->ZIa = 0.e0;
    }
#endif //TASK_GALACTIC_CENTER

    return;
}

void BinaryDump(void){

    Pall.OutPutFileNumber ++;
    char fname[MaxCharactersInLine];

    char fname_old[MaxCharactersInLine];
    sprintf(fname,"%s.dump.%03d.%03d",Pall.BaseFileName,MPIGetNumProcs(),MPIGetMyID());
    char fname_new[MaxCharactersInLine];
    sprintf(fname,"%s.dump.old.%03d.%03d",Pall.BaseFileName,MPIGetNumProcs(),MPIGetMyID());
    rename(fname_old,fname_new);
    MPI_Barrier(MPI_COMM_WORLD);

    sprintf(fname,"%s.dump.%03d.%03d",Pall.BaseFileName,MPIGetNumProcs(),MPIGetMyID());
    ParallelWriteAllData(fname);

    return;
}

static void ParallelWriteAllDataFullIO(FILE *fp){

    for(int k=0;k<Pall.Ntotal;k++){
        fwrite(Pbody[k],sizeof(StructPbody),1,fp);
        if(Pbody[k]->Type == TypeHydro){ // type hydro
            fwrite(Pbody[k]->Baryon,sizeof(StructPhydro),1,fp);
        } else if (Pbody[k]->Type == TypeStar){ // type star
            fwrite(Pbody[k]->Baryon,sizeof(StructPstar),1,fp);
        } else if (Pbody[k]->Type == TypeSink){ // type sink
            fwrite(Pbody[k]->Baryon,sizeof(StructPsink),1,fp);
        } else if (Pbody[k]->Type != TypeDM){
            fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                    MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
            MPI_Finalize();
            MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
        }
    }

    return ;
}

static void ParallelWriteAllDataCompactIO(FILE *fp){

    for(int k=0;k<Pall.Ntotal;k++){
        struct StructPbodyIOCompact PbodyIO = CopyPbodyToTemporalStructureCompact(k);
        fwrite(&PbodyIO,sizeof(struct StructPbodyIOCompact),1,fp);
        if(Pbody[k]->Type == TypeHydro){ // type hydro
            struct StructPhydroIOCompact PhydroIO = CopyPhydroToTemporalStructureCompactElement(Pbody[k]->Baryon);
            fwrite(&PhydroIO,sizeof(struct StructPhydroIOCompact),1,fp);
        } else if (Pbody[k]->Type == TypeStar){ // type star
            struct StructPstarIOCompact PstarIO = CopyPstarToTemporalStructureCompactElement(Pbody[k]->Baryon);
            fwrite(&PstarIO,sizeof(struct StructPstarIOCompact),1,fp);
        } else if (Pbody[k]->Type == TypeSink){ // type sink
            struct StructPsinkIOCompact PsinkIO = CopyPsinkToTemporalStructureCompactElement(Pbody[k]->Baryon);
            fwrite(&PsinkIO,sizeof(struct StructPsinkIOCompact),1,fp);
        } else if (Pbody[k]->Type != TypeDM){
            fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                    MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
            MPI_Finalize();
            MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
        }
    }

    return ;
}

static void ParallelWriteAllDataCompactDoubleIO(FILE *fp){

    for(int k=0;k<Pall.Ntotal;k++){
        struct StructPbodyIOCompactDouble PbodyIO = CopyPbodyToTemporalStructureCompactDouble(k);
        fwrite(&PbodyIO,sizeof(struct StructPbodyIOCompactDouble),1,fp);
        if(Pbody[k]->Type == TypeHydro){ // type hydro
            struct StructPhydroIOCompactDouble PhydroIO = CopyPhydroToTemporalStructureCompactDoubleElement(Pbody[k]->Baryon);
            fwrite(&PhydroIO,sizeof(struct StructPhydroIOCompactDouble),1,fp);
        } else if (Pbody[k]->Type == TypeStar){ // type star
            struct StructPstarIOCompactDouble PstarIO = CopyPstarToTemporalStructureCompactDoubleElement(Pbody[k]->Baryon);
            fwrite(&PstarIO,sizeof(struct StructPstarIOCompactDouble),1,fp);
        } else if (Pbody[k]->Type == TypeSink){ // type sink
            struct StructPsinkIOCompactDouble PsinkIO = CopyPsinkToTemporalStructureCompactDoubleElement(Pbody[k]->Baryon);
            fwrite(&PsinkIO,sizeof(struct StructPsinkIOCompactDouble),1,fp);
        } else if (Pbody[k]->Type != TypeDM){
            fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                    MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
            MPI_Finalize();
            MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
        }
    }

    return ;
}

static void ParallelWriteAllDataCompactMixIO(FILE *fp){

    for(int k=0;k<Pall.Ntotal;k++){
        struct StructPbodyIOCompactMix PbodyIO = CopyPbodyToTemporalStructureCompactMix(k);
        fwrite(&PbodyIO,sizeof(struct StructPbodyIOCompactMix),1,fp);
        if(Pbody[k]->Type == TypeHydro){ // type hydro
            struct StructPhydroIOCompactMix PhydroIO = CopyPhydroToTemporalStructureCompactMixElement(Pbody[k]->Baryon);
            fwrite(&PhydroIO,sizeof(struct StructPhydroIOCompactMix),1,fp);
        } else if (Pbody[k]->Type == TypeStar){ // type star
            struct StructPstarIOCompactMix PstarIO = CopyPstarToTemporalStructureCompactMixElement(Pbody[k]->Baryon);
            fwrite(&PstarIO,sizeof(struct StructPstarIOCompactMix),1,fp);
        } else if (Pbody[k]->Type == TypeSink){ // type sink
            struct StructPsinkIOCompactMix PsinkIO = CopyPsinkToTemporalStructureCompactMixElement(Pbody[k]->Baryon);
            fwrite(&PsinkIO,sizeof(struct StructPsinkIOCompactMix),1,fp);
        } else if (Pbody[k]->Type != TypeDM){
            fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                    MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
            MPI_Finalize();
            MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
        }
    }

    return ;
}

#ifdef USE_LEAN_IO_FORMAT
static void ParallelWriteAllDataLeanIO(FILE *fp){

    for(int k=0;k<Pall.Ntotal;k++){
        struct StructPbodyIOLean PbodyIO = CopyPbodyToTemporalStructureLean(k);
        fwrite(&PbodyIO,sizeof(struct StructPbodyIOLean),1,fp);
        if(Pbody[k]->Type == TypeHydro){ // type hydro
            struct StructPhydroIOLean PhydroIO = CopyPhydroToTemporalStructureLeanElement(Pbody[k]->Baryon);
            fwrite(&PhydroIO,sizeof(struct StructPhydroIOLean),1,fp);
        } else if (Pbody[k]->Type == TypeStar){ // type star
            struct StructPstarIOLean PstarIO = CopyPstarToTemporalStructureLeanElement(Pbody[k]->Baryon);
            fwrite(&PstarIO,sizeof(struct StructPstarIOLean),1,fp);
        } else if (Pbody[k]->Type == TypeSink){ // type sink
            struct StructPsinkIOLean PsinkIO = CopyPsinkToTemporalStructureLeanElement(Pbody[k]->Baryon);
            fwrite(&PsinkIO,sizeof(struct StructPsinkIOLean),1,fp);
        } else if (Pbody[k]->Type != TypeDM){
            fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                    MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
            MPI_Finalize();
            MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
        }
    }

    return ;
}
#endif

void ParallelWriteAllData(char fname[]){

    FILE *fp;

    FileOpen(fp,fname,"wb");
    // write header //
    WriteHeader(fp);
    // write header //
    fwrite(&Pall,sizeof(struct StructPall),1,fp);
    fwrite(&TimingResults,sizeof(struct StructTimingResults),1,fp);

#ifdef COMPACT_IO_FORMAT
#ifdef USE_LEAN_IO_FORMAT
    if(LeanIOCondition()){
        ParallelWriteAllDataLeanIO(fp);
    } else {
        ParallelWriteAllDataCompactIO(fp);
    }
#else // USE_LEAN_IO_FORMAT
    ParallelWriteAllDataCompactIO(fp);
#endif // USE_LEAN_IO_FORMAT
#elif defined(COMPACT_DOUBLE_IO_FORMAT) // COMPACT_IO_FORMAT
    ParallelWriteAllDataCompactDoubleIO(fp);
#elif defined(COMPACT_MIX_IO_FORMAT) // COMPACT_IO_FORMAT
    ParallelWriteAllDataCompactMixIO(fp);
#else // COMPACT_IO_FORMAT
    ParallelWriteAllDataFullIO(fp);
#endif // COMPACT_IO_FORMAT

    if(GSL_EFAILED == gsl_rng_fwrite(fp,RandomGenerator)){
        fprintf(stderr,"The random generator write error!\n");
        MPI_Abort(MPI_COMM_WORLD,GSL_RANDOM_GENERATOR_WRITE_ERROR);
    }
    fclose(fp);

    return;
}

void ParallelDumpAllData(char fname[]){

    FILE *fp;

    FileOpen(fp,fname,"wb");
    // write header //
    WriteHeader(fp);
    // write header //
    fwrite(&Pall,sizeof(struct StructPall),1,fp);
    fwrite(&TimingResults,sizeof(struct StructTimingResults),1,fp);

    for(int k=0;k<Pall.Ntotal;k++){
        fwrite(Pbody[k],sizeof(StructPbody),1,fp);
        if(Pbody[k]->Type == TypeHydro){ // type hydro
            fwrite(Pbody[k]->Baryon,sizeof(StructPhydro),1,fp);
        } else if (Pbody[k]->Type == TypeStar){ // type star
            fwrite(Pbody[k]->Baryon,sizeof(StructPstar),1,fp);
        } else if (Pbody[k]->Type == TypeSink){ // type sink
            fwrite(Pbody[k]->Baryon,sizeof(StructPsink),1,fp);
        } else if (Pbody[k]->Type != TypeDM){
            fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                    MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
            MPI_Finalize();
            MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
        }
    }
    if(GSL_EFAILED == gsl_rng_fwrite(fp,RandomGenerator)){
        fprintf(stderr,"The random generator write error!\n");
        MPI_Abort(MPI_COMM_WORLD,GSL_RANDOM_GENERATOR_WRITE_ERROR);
    }
    fclose(fp);

    return;
}

void ParallelWriteAllDataASCIIFormat(char fname[]){

    FILE *fp;
    FileOpen(fp,fname,"w");

    for(int i=0;i<Pall.Nhydro;i++){
#ifdef USE_GRAD_H //{
        double Gradh = Phydro[i]->Gradh;
#else
        double Gradh = 0.e0;
#endif // USE_GRAD_H //}

        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g %g %g\n",
                PhydroBody(i)->GlobalID,
                PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2],
                Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2],
                Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->Nlist,Gradh,
                Phydro[i]->U,Pall.Gm1*Phydro[i]->Rho*Phydro[i]->U,
                Phydro[i]->Du,Phydro[i]->EnergyDensity);
    }
    fclose(fp);
    fflush(NULL);

    return ;
}


/*
 * Read restart files. 
 */
void ParallelReadAllData(void){

    FILE *fp;
    char fname[MaxCharactersInLine];
    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int HeaderArray[HEADER_SIZE];

    sprintf(fname,"%s.%03d.%03d",Pall.RestartFileName,NProcs,MyID);
    FileOpen(fp,fname,"rb");
    ReadHeader(fp,HeaderArray);
    fread(&Pall,sizeof(struct StructPall),1,fp);
    fread(&TimingResults,sizeof(struct StructTimingResults),1,fp);

    GenerateStructPbody(Pall.Ntotal);
    GenerateStructPhydro(Pall.Nhydro);
    GenerateStructPstar(Pall.Nstars);
    GenerateStructPsink(Pall.Nsink);

    bool RestartFromDump = FindWordFromString(fname,"temp",".");

    int CountTypes[NTypes];
    for(int i=0;i<NTypes;i++)
        CountTypes[i] = 0;

    if( ((HeaderArray[HD_CompactFormatFlag] == OFF)&&(HeaderArray[HD_CompactDoubleFormatFlag] == OFF)&&(HeaderArray[HD_CompactMixFormatFlag] == OFF)) || (RestartFromDump==true)){
        for(int k=0;k<Pall.Ntotal;k++){
            StructPbodyptr Pb = PbodyElements[k].Next;
            fread(PbodyElements+k,sizeof(StructPbody),1,fp);
            PbodyElements[k].Next = Pb;

            if(PbodyElements[k].Type == TypeHydro){ // type hydro
                StructPhydroptr Ph = PhydroElements[CountTypes[TypeHydro]].Next;
                fread(PhydroElements+CountTypes[TypeHydro],sizeof(StructPhydro),1,fp);
                PhydroElements[CountTypes[TypeHydro]].Next = Ph;
                PhydroElements[CountTypes[TypeHydro]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PhydroElements+CountTypes[TypeHydro]);
                CountTypes[TypeHydro] ++;
            } else if (PbodyElements[k].Type == TypeStar){ // type star
                StructPstarptr Ps = PstarElements[CountTypes[TypeStar]].Next;
                fread(PstarElements+CountTypes[TypeStar],sizeof(StructPstar),1,fp);
                PstarElements[CountTypes[TypeStar]].Next = Ps;
                PstarElements[CountTypes[TypeStar]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PstarElements+CountTypes[TypeStar]);
                CountTypes[TypeStar] ++;
            } else if (PbodyElements[k].Type == TypeSink){ // type sink
                StructPsinkptr Psk = PsinkElements[CountTypes[TypeSink]].Next;
                fread(PsinkElements+CountTypes[TypeSink],sizeof(StructPsink),1,fp);
                PsinkElements[CountTypes[TypeSink]].Next = Psk;
                PsinkElements[CountTypes[TypeSink]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PsinkElements+CountTypes[TypeSink]);
                CountTypes[TypeSink] ++;
            } else if (PbodyElements[k].Type == TypeDM){
                CountTypes[TypeDM] ++;
            } else if (PbodyElements[k].Type != TypeDM){
                fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                        MPIGetMyID(),PbodyElements[k].Type,__FILE__,__LINE__);
                MPI_Finalize();
                assert(1==0);
                exit(ReadTypeUnDefined);
            }
        }
    } else if(HeaderArray[HD_CompactFormatFlag] == ON){
        for(int k=0;k<Pall.Ntotal;k++){
            StructPbodyptr Pb = PbodyElements[k].Next;
            struct StructPbodyIOCompact   PbodyIO;
            fread(&PbodyIO,sizeof(struct StructPbodyIOCompact),1,fp);

            PbodyElements[k] = CopyTemporalStructureCompactToPbodyCompact(PbodyIO);
            PbodyElements[k].Use = ON;
            PbodyElements[k].Next = Pb;
            /* predictor copy */
#if 1
            for(int l=0;l<3;l++)
                PbodyElements[k].PosP[l] = PbodyElements[k].Pos[l];
#endif
            /* predictor copy */

            if(PbodyElements[k].Type == TypeHydro){ // type hydro
                StructPhydroptr Ph = PhydroElements[CountTypes[TypeHydro]].Next;
                struct StructPhydroIOCompact  PhydroIO;
                fread(&PhydroIO,sizeof(struct StructPhydroIOCompact),1,fp);

                PhydroElements[CountTypes[TypeHydro]] = CopyTemporalStructureCompactToPhydroCompact(PhydroIO);
                PhydroElements[CountTypes[TypeHydro]].Use = ON;
                PhydroElements[CountTypes[TypeHydro]].Next = Ph;
                PhydroElements[CountTypes[TypeHydro]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PhydroElements+CountTypes[TypeHydro]);
                /* predictor copy */
                PhydroElements[CountTypes[TypeHydro]].Mass = PbodyElements[k].Mass;
#if 1
                for(int l=0;l<3;l++){
                    PhydroElements[CountTypes[TypeHydro]].PosP[l] = PbodyElements[k].Pos[l];
                    PhydroElements[CountTypes[TypeHydro]].VelP[l] = PbodyElements[k].Vel[l];
                }
                PhydroElements[CountTypes[TypeHydro]].RhoPred = PhydroElements[CountTypes[TypeHydro]].Rho;
                PhydroElements[CountTypes[TypeHydro]].KernelPred = PhydroElements[CountTypes[TypeHydro]].Kernel;
                PhydroElements[CountTypes[TypeHydro]].UPred = PhydroElements[CountTypes[TypeHydro]].U;
                PhydroElements[CountTypes[TypeHydro]].Active = PhydroElements[CountTypes[TypeHydro]].Active;
#endif
                /* predictor copy */
                CountTypes[TypeHydro] ++;
            } else if (PbodyElements[k].Type == TypeStar){ // type star
                StructPstarptr Ps = PstarElements[CountTypes[TypeStar]].Next;
                struct StructPstarIOCompact   PstarIO;
                fread(&PstarIO,sizeof(struct StructPstarIOCompact),1,fp);

                PstarElements[CountTypes[TypeStar]] = CopyTemporalStructureCompactToPstarCompact(PstarIO);
                PstarElements[CountTypes[TypeStar]].Use = ON;
                PstarElements[CountTypes[TypeStar]].Next = Ps;
                PstarElements[CountTypes[TypeStar]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PstarElements+CountTypes[TypeStar]);
                CountTypes[TypeStar] ++;
            } else if (PbodyElements[k].Type == TypeSink){ // type sink
                StructPsinkptr Psk = PsinkElements[CountTypes[TypeSink]].Next;
                struct StructPsinkIOCompact   PsinkIO;
                fread(&PsinkIO,sizeof(struct StructPsinkIOCompact),1,fp);

                PsinkElements[CountTypes[TypeSink]] = CopyTemporalStructureCompactToPsinkCompact(PsinkIO);
                PsinkElements[CountTypes[TypeSink]].Use = ON;
                PsinkElements[CountTypes[TypeSink]].Next = Psk;
                PsinkElements[CountTypes[TypeSink]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PsinkElements+CountTypes[TypeSink]);
                for(int l=0;l<3;l++){
                    PsinkElements[CountTypes[TypeSink]].PosP[l] = PbodyElements[k].Pos[l];
                    PsinkElements[CountTypes[TypeSink]].VelP[l] = PbodyElements[k].Vel[l];
                }

                CountTypes[TypeSink] ++;
                // dprintlmpi(CountTypes[TypeSink]);
                //fflush(NULL);
            } else if (PbodyElements[k].Type == TypeDM){
                CountTypes[TypeDM] ++;
            } else if (PbodyElements[k].Type != TypeDM){
                fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                        MPIGetMyID(),PbodyElements[k].Type,__FILE__,__LINE__);
                MPI_Finalize();
                exit(ReadTypeUnDefined);
            }
        }
    } else if(HeaderArray[HD_CompactDoubleFormatFlag] == ON){
        for(int k=0;k<Pall.Ntotal;k++){
            StructPbodyptr Pb = PbodyElements[k].Next;
            struct StructPbodyIOCompactDouble   PbodyIO;
            fread(&PbodyIO,sizeof(struct StructPbodyIOCompactDouble),1,fp);

            PbodyElements[k] = CopyTemporalStructureCompactDoubleToPbodyCompactDouble(PbodyIO);
            PbodyElements[k].Use = ON;
            PbodyElements[k].Next = Pb;
            /* predictor copy */
#if 1
            for(int l=0;l<3;l++)
                PbodyElements[k].PosP[l] = PbodyElements[k].Pos[l];
#endif
            /* predictor copy */

            if(PbodyElements[k].Type == TypeHydro){ // type hydro
                StructPhydroptr Ph = PhydroElements[CountTypes[TypeHydro]].Next;
                struct StructPhydroIOCompactDouble  PhydroIO;
                fread(&PhydroIO,sizeof(struct StructPhydroIOCompactDouble),1,fp);

                PhydroElements[CountTypes[TypeHydro]] = CopyTemporalStructureCompactDoubleToPhydroCompactDouble(PhydroIO);
                PhydroElements[CountTypes[TypeHydro]].Use = ON;
                PhydroElements[CountTypes[TypeHydro]].Next = Ph;
                PhydroElements[CountTypes[TypeHydro]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PhydroElements+CountTypes[TypeHydro]);
                /* predictor copy */
                PhydroElements[CountTypes[TypeHydro]].Mass = PbodyElements[k].Mass;
#if 1
                for(int l=0;l<3;l++){
                    PhydroElements[CountTypes[TypeHydro]].PosP[l] = PbodyElements[k].Pos[l];
                    PhydroElements[CountTypes[TypeHydro]].VelP[l] = PbodyElements[k].Vel[l];
                }
                PhydroElements[CountTypes[TypeHydro]].RhoPred = PhydroElements[CountTypes[TypeHydro]].Rho;
                PhydroElements[CountTypes[TypeHydro]].KernelPred = PhydroElements[CountTypes[TypeHydro]].Kernel;
                PhydroElements[CountTypes[TypeHydro]].UPred = PhydroElements[CountTypes[TypeHydro]].U;
                PhydroElements[CountTypes[TypeHydro]].Active = PhydroElements[CountTypes[TypeHydro]].Active;
#endif
                /* predictor copy */
                CountTypes[TypeHydro] ++;
            } else if (PbodyElements[k].Type == TypeStar){ // type star
                StructPstarptr Ps = PstarElements[CountTypes[TypeStar]].Next;
                struct StructPstarIOCompactDouble   PstarIO;
                fread(&PstarIO,sizeof(struct StructPstarIOCompactDouble),1,fp);

                PstarElements[CountTypes[TypeStar]] = CopyTemporalStructureCompactDoubleToPstarCompactDouble(PstarIO);
                PstarElements[CountTypes[TypeStar]].Use = ON;
                PstarElements[CountTypes[TypeStar]].Next = Ps;
                PstarElements[CountTypes[TypeStar]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PstarElements+CountTypes[TypeStar]);
                CountTypes[TypeStar] ++;
            } else if (PbodyElements[k].Type == TypeSink){ // type sink
                StructPsinkptr Psk = PsinkElements[CountTypes[TypeSink]].Next;
                struct StructPsinkIOCompactDouble   PsinkIO;
                fread(&PsinkIO,sizeof(struct StructPsinkIOCompactDouble),1,fp);

                PsinkElements[CountTypes[TypeSink]] = CopyTemporalStructureCompactDoubleToPsinkCompactDouble(PsinkIO);
                PsinkElements[CountTypes[TypeSink]].Use = ON;
                PsinkElements[CountTypes[TypeSink]].Next = Psk;
                PsinkElements[CountTypes[TypeSink]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PsinkElements+CountTypes[TypeSink]);
                for(int l=0;l<3;l++){
                    PsinkElements[CountTypes[TypeSink]].PosP[l] = PbodyElements[k].Pos[l];
                    PsinkElements[CountTypes[TypeSink]].VelP[l] = PbodyElements[k].Vel[l];
                }

                CountTypes[TypeSink] ++;
                // dprintlmpi(CountTypes[TypeSink]);
                //fflush(NULL);
            } else if (PbodyElements[k].Type == TypeDM){
                CountTypes[TypeDM] ++;
            } else if (PbodyElements[k].Type != TypeDM){
                fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                        MPIGetMyID(),PbodyElements[k].Type,__FILE__,__LINE__);
                MPI_Finalize();
                exit(ReadTypeUnDefined);
            }
        }
    } else if(HeaderArray[HD_CompactMixFormatFlag] == ON){
        for(int k=0;k<Pall.Ntotal;k++){
            StructPbodyptr Pb = PbodyElements[k].Next;
            struct StructPbodyIOCompactMix   PbodyIO;
            fread(&PbodyIO,sizeof(struct StructPbodyIOCompactMix),1,fp);

            PbodyElements[k] = CopyTemporalStructureCompactMixToPbodyCompactMix(PbodyIO);
            PbodyElements[k].Use = ON;
            PbodyElements[k].Next = Pb;
            /* predictor copy */
#if 1
            for(int l=0;l<3;l++)
                PbodyElements[k].PosP[l] = PbodyElements[k].Pos[l];
#endif
            /* predictor copy */

            if(PbodyElements[k].Type == TypeHydro){ // type hydro
                StructPhydroptr Ph = PhydroElements[CountTypes[TypeHydro]].Next;
                struct StructPhydroIOCompactMix  PhydroIO;
                fread(&PhydroIO,sizeof(struct StructPhydroIOCompactMix),1,fp);

                PhydroElements[CountTypes[TypeHydro]] = CopyTemporalStructureCompactMixToPhydroCompactMix(PhydroIO);
                PhydroElements[CountTypes[TypeHydro]].Use = ON;
                PhydroElements[CountTypes[TypeHydro]].Next = Ph;
                PhydroElements[CountTypes[TypeHydro]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PhydroElements+CountTypes[TypeHydro]);
                /* predictor copy */
                PhydroElements[CountTypes[TypeHydro]].Mass = PbodyElements[k].Mass;
#if 1
                for(int l=0;l<3;l++){
                    PhydroElements[CountTypes[TypeHydro]].PosP[l] = PbodyElements[k].Pos[l];
                    PhydroElements[CountTypes[TypeHydro]].VelP[l] = PbodyElements[k].Vel[l];
                }
                PhydroElements[CountTypes[TypeHydro]].RhoPred = PhydroElements[CountTypes[TypeHydro]].Rho;
                PhydroElements[CountTypes[TypeHydro]].KernelPred = PhydroElements[CountTypes[TypeHydro]].Kernel;
                PhydroElements[CountTypes[TypeHydro]].UPred = PhydroElements[CountTypes[TypeHydro]].U;
                PhydroElements[CountTypes[TypeHydro]].Active = PhydroElements[CountTypes[TypeHydro]].Active;
#endif
                /* predictor copy */
                CountTypes[TypeHydro] ++;
            } else if (PbodyElements[k].Type == TypeStar){ // type star
                StructPstarptr Ps = PstarElements[CountTypes[TypeStar]].Next;
                struct StructPstarIOCompactMix   PstarIO;
                fread(&PstarIO,sizeof(struct StructPstarIOCompactMix),1,fp);

                PstarElements[CountTypes[TypeStar]] = CopyTemporalStructureCompactMixToPstarCompactMix(PstarIO);
                PstarElements[CountTypes[TypeStar]].Use = ON;
                PstarElements[CountTypes[TypeStar]].Next = Ps;
                PstarElements[CountTypes[TypeStar]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PstarElements+CountTypes[TypeStar]);
                CountTypes[TypeStar] ++;
            } else if (PbodyElements[k].Type == TypeSink){ // type sink
                StructPsinkptr Psk = PsinkElements[CountTypes[TypeSink]].Next;
                struct StructPsinkIOCompactMix   PsinkIO;
                fread(&PsinkIO,sizeof(struct StructPsinkIOCompactMix),1,fp);

                PsinkElements[CountTypes[TypeSink]] = CopyTemporalStructureCompactMixToPsinkCompactMix(PsinkIO);
                PsinkElements[CountTypes[TypeSink]].Use = ON;
                PsinkElements[CountTypes[TypeSink]].Next = Psk;
                PsinkElements[CountTypes[TypeSink]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PsinkElements+CountTypes[TypeSink]);
                for(int l=0;l<3;l++){
                    PsinkElements[CountTypes[TypeSink]].PosP[l] = PbodyElements[k].Pos[l];
                    PsinkElements[CountTypes[TypeSink]].VelP[l] = PbodyElements[k].Vel[l];
                }

                CountTypes[TypeSink] ++;
                // dprintlmpi(CountTypes[TypeSink]);
                //fflush(NULL);
            } else if (PbodyElements[k].Type == TypeDM){
                CountTypes[TypeDM] ++;
            } else if (PbodyElements[k].Type != TypeDM){
                fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                        MPIGetMyID(),PbodyElements[k].Type,__FILE__,__LINE__);
                MPI_Finalize();
                exit(ReadTypeUnDefined);
            }
        }
    }else{
        fprintf(stderr,"Header flag error!\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(HeaderError);
    }

    AllocateRandomGenerator();
    if(GSL_EFAILED == gsl_rng_fread(fp,RandomGenerator)){
        fprintf(stderr,"The random generator read error!\n");
        MPI_Abort(MPI_COMM_WORLD,GSL_RANDOM_GENERATOR_READ_ERROR);
    }

    fclose(fp);

    fprintf(stderr," [%02d] Reload RestartFileData : (Ntotal,Ndm,Nhydro,Nstar,Nsink) = (%d,%d,%d,%d,%d)\n",MyID,
            CountTypes[TypeDM]+CountTypes[TypeHydro]+CountTypes[TypeStar]+CountTypes[TypeSink],
            CountTypes[TypeDM],CountTypes[TypeHydro],CountTypes[TypeStar],CountTypes[TypeSink]);

    // ReConnectPointers();
    return;
}

static void ParallelReadWriteAllData(void){

    FILE *fp;
    char fname[MaxCharactersInLine];
    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int HeaderArray[HEADER_SIZE];
    int RunStatusOperation = Pall.RunStatus;

    sprintf(fname,"%s.%03d.%03d",Pall.RestartFileName,NProcs,MyID);
    FileOpen(fp,fname,"rb");
    ReadHeader(fp,HeaderArray);
    fread(&Pall,sizeof(struct StructPall),1,fp);
    Pall.RunStatus = RunStatusOperation;
    fread(&TimingResults,sizeof(struct StructTimingResults),1,fp);

    GenerateStructPbody(Pall.Ntotal);
    if(Pall.Nhydro_t>0)
        GenerateStructPhydro(Pall.Nhydro);
    if(Pall.Nstars_t>0)
        GenerateStructPstar(Pall.Nstars);

    int CountTypes[NTypes];
    for(int i=0;i<NTypes;i++)
        CountTypes[i] = 0;

    if(HeaderArray[HD_CompactFormatFlag] == OFF){
        for(int k=0;k<Pall.Ntotal;k++){
            StructPbodyptr Pb = PbodyElements[k].Next;
            fread(PbodyElements+k,sizeof(StructPbody),1,fp);
            PbodyElements[k].Next = Pb;

            if(PbodyElements[k].Type == TypeHydro){ // type hydro
                StructPhydroptr Ph = PhydroElements[CountTypes[TypeHydro]].Next;
                fread(PhydroElements+CountTypes[TypeHydro],sizeof(StructPhydro),1,fp);
                PhydroElements[CountTypes[TypeHydro]].Next = Ph;
                PhydroElements[CountTypes[TypeHydro]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PhydroElements+CountTypes[TypeHydro]);
                CountTypes[TypeHydro] ++;
            } else if (PbodyElements[k].Type == TypeStar){ // type star
                StructPstarptr Ps = PstarElements[CountTypes[TypeStar]].Next;
                fread(PstarElements+CountTypes[TypeStar],sizeof(StructPstar),1,fp);
                PstarElements[CountTypes[TypeStar]].Next = Ps;
                PstarElements[CountTypes[TypeStar]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PstarElements+CountTypes[TypeStar]);
                CountTypes[TypeStar] ++;
            } else if (PbodyElements[k].Type == TypeDM){
                CountTypes[TypeDM] ++;
            } else if (PbodyElements[k].Type != TypeDM){
                fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                        MPIGetMyID(),PbodyElements[k].Type,__FILE__,__LINE__);
                MPI_Finalize();
                exit(ReadTypeUnDefined);
            }
        }
    } else if(HeaderArray[HD_CompactFormatFlag] == ON){
        for(int k=0;k<Pall.Ntotal;k++){
            StructPbodyptr Pb = PbodyElements[k].Next;
            struct StructPbodyIOCompact   PbodyIO;
            fread(&PbodyIO,sizeof(struct StructPbodyIOCompact),1,fp);

            PbodyElements[k] = CopyTemporalStructureCompactToPbodyCompact(PbodyIO);
            PbodyElements[k].Use = ON;
            PbodyElements[k].Next = Pb;
            /* predictor copy */
#if 1
            for(int l=0;l<3;l++)
                PbodyElements[k].PosP[l] = PbodyElements[k].Pos[l];
#endif
            /* predictor copy */

            if(PbodyElements[k].Type == TypeHydro){ // type hydro
                StructPhydroptr Ph = PhydroElements[CountTypes[TypeHydro]].Next;
                struct StructPhydroIOCompact  PhydroIO;
                fread(&PhydroIO,sizeof(struct StructPhydroIOCompact),1,fp);

                PhydroElements[CountTypes[TypeHydro]] = CopyTemporalStructureCompactToPhydroCompact(PhydroIO);
                PhydroElements[CountTypes[TypeHydro]].Use = ON;
                PhydroElements[CountTypes[TypeHydro]].Next = Ph;
                PhydroElements[CountTypes[TypeHydro]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PhydroElements+CountTypes[TypeHydro]);
                /* predictor copy */
                PhydroElements[CountTypes[TypeHydro]].Mass = PbodyElements[k].Mass;
#if 1
                for(int l=0;l<3;l++){
                    PhydroElements[CountTypes[TypeHydro]].PosP[l] = PbodyElements[k].Pos[l];
                    PhydroElements[CountTypes[TypeHydro]].VelP[l] = PbodyElements[k].Vel[l];
                }
                PhydroElements[CountTypes[TypeHydro]].RhoPred = PhydroElements[CountTypes[TypeHydro]].Rho;
                PhydroElements[CountTypes[TypeHydro]].KernelPred = PhydroElements[CountTypes[TypeHydro]].Kernel;
                PhydroElements[CountTypes[TypeHydro]].UPred = PhydroElements[CountTypes[TypeHydro]].U;
                PhydroElements[CountTypes[TypeHydro]].Active = PhydroElements[CountTypes[TypeHydro]].Active;
#endif
                /* predictor copy */
                CountTypes[TypeHydro] ++;
            } else if (PbodyElements[k].Type == TypeStar){ // type hydro
                StructPstarptr Ps = PstarElements[CountTypes[TypeStar]].Next;
                struct StructPstarIOCompact   PstarIO;
                fread(&PstarIO,sizeof(struct StructPstarIOCompact),1,fp);

                PstarElements[CountTypes[TypeStar]] = CopyTemporalStructureCompactToPstarCompact(PstarIO);
                PstarElements[CountTypes[TypeStar]].Use = ON;
                PstarElements[CountTypes[TypeStar]].Next = Ps;
                PstarElements[CountTypes[TypeStar]].Body = PbodyElements+k;
                PbodyElements[k].Baryon = (void*)(PstarElements+CountTypes[TypeStar]);
                CountTypes[TypeStar] ++;
            } else if (PbodyElements[k].Type == TypeDM){
                CountTypes[TypeDM] ++;
            } else if (PbodyElements[k].Type != TypeDM){
                fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                        MPIGetMyID(),PbodyElements[k].Type,__FILE__,__LINE__);
                MPI_Finalize();
                exit(ReadTypeUnDefined);
            }
        }
    }else{
        fprintf(stderr,"Header flag error!\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(HeaderError);
    }

    AllocateRandomGenerator();
    if(GSL_EFAILED == gsl_rng_fread(fp,RandomGenerator)){
        fprintf(stderr,"The random generator read error!\n");
        MPI_Abort(MPI_COMM_WORLD,GSL_RANDOM_GENERATOR_READ_ERROR);
    }

    fclose(fp);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr," [%02d] Reload RestartFileData : (Ntotal,Ndm,Nhydro,Nstar) = (%d,%d,%d,%d)\n",MyID,
                CountTypes[TypeDM]+CountTypes[TypeHydro]+CountTypes[TypeStar],
                CountTypes[TypeDM],CountTypes[TypeHydro],CountTypes[TypeStar]);
    }

    ReConnectPointers();

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"End data read!\n");

    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(Pall.RunStatus == DataFileOperation_Double){
            fprintf(stderr,"Double\n");
        } else if(Pall.RunStatus == DataFileOperation_Halve){
            fprintf(stderr,"Halve\n");
        }
    }


    /// WritePartON
    if(Pall.RunStatus == DataFileOperation_Double){
        for(int i=0;i<MPIGetNumProcs();i++){
            if(MPIGetMyID() == i){
                fprintf(stderr,"Start! %d\n",MPIGetMyID());

                struct  StructPall     TempPall_first = Pall;
                struct  StructPall     TempPall_second = Pall;

                int Nfirst = Pall.Ntotal/2;
                //int Nsecond = Pall.Ntotal-Nfirst;

                TempPall_first.Ntotal = TempPall_first.NDM = 
                TempPall_first.Nhydro = TempPall_first.Nstars = 0;
                for(int k=0;k<Nfirst;k++){ // count first.
                    if(Pbody[k]->Type == TypeHydro){
                        TempPall_first.Nhydro ++;
                    } else if(Pbody[k]->Type == TypeStar){
                        TempPall_first.Nstars ++;
                    } else if (Pbody[k]->Type == TypeDM){
                        TempPall_first.NDM ++;
                    } else {
                        exit(1);
                    }
                    TempPall_first.Ntotal ++;
                }
                TempPall_second.Ntotal = TempPall_second.NDM = 
                TempPall_second.Nhydro = TempPall_second.Nstars = 0;
                for(int k=Nfirst;k<Pall.Ntotal;k++){ // count second.
                    if(Pbody[k]->Type == TypeHydro){
                        TempPall_second.Nhydro ++;
                    } else if(Pbody[k]->Type == TypeStar){
                        TempPall_second.Nstars ++;
                    } else if (Pbody[k]->Type == TypeDM){
                        TempPall_second.NDM ++;
                    } else {
                        exit(1);
                    }
                    TempPall_second.Ntotal ++;
                }

                long int Ntotal,NDM,Nhydro,Nstars;
                Ntotal = NDM = Nhydro = Nstars = 0;
                for(int k=0;k<Pall.Ntotal;k++){ // count second.
                    if(Pbody[k]->Type == TypeHydro){
                        Nhydro ++;
                    } else if(Pbody[k]->Type == TypeStar){
                        Nstars ++;
                    } else if(Pbody[k]->Type == TypeDM){
                        NDM ++;
                    } else {
                        exit(1);
                    }
                    Ntotal ++;
                }

                fprintf(stderr,"[%d] -f- %ld %ld %ld %ld\n",MPIGetMyID(),
                        TempPall_first.Ntotal,TempPall_first.NDM,TempPall_first.Nhydro,TempPall_first.Nstars);
                fprintf(stderr,"[%d] -s- %ld %ld %ld %ld\n",MPIGetMyID(),
                        TempPall_second.Ntotal,TempPall_second.NDM,TempPall_second.Nhydro,TempPall_second.Nstars);
                fprintf(stderr,"[%d] -A- %ld %ld %ld %ld\n",MPIGetMyID(),Ntotal,NDM,Nhydro,Nstars);
                fprintf(stderr,"[%d] -o- %ld %ld %ld %ld\n",MPIGetMyID(),
                        Pall.Ntotal,Pall.NDM,Pall.Nhydro,Pall.Nstars);
                fprintf(stderr,"[%d] -S- %ld %ld %ld %ld\n",MPIGetMyID(),
                    TempPall_first.Ntotal+TempPall_second.Ntotal,TempPall_first.NDM+TempPall_second.NDM,
                    TempPall_first.Nhydro+TempPall_second.Nhydro,TempPall_first.Nstars+TempPall_second.Nstars);

                //assert(Pall.Ntotal == TempPall_first.Ntotal+TempPall_second.Ntotal);
                //assert(Pall.NDM    == TempPall_first.NDM+TempPall_second.NDM);
                //assert(Pall.Nhydro == TempPall_first.Nhydro+TempPall_second.Nhydro);
                //assert(Pall.Nstars == TempPall_first.Nstars+TempPall_second.Nstars);

                char fname[MaxCharactersInLine];

                // write first.
                sprintf(fname,"%s.%04d.%03d.%03d",Pall.BaseFileName,Pall.OutPutFileNumber,
                        2*MPIGetNumProcs(),2*MPIGetMyID());
                FileOpen(fp,fname,"wb");

                // write header //
                WriteHeader(fp);
                // write header //
                fwrite(&TempPall_first,sizeof(struct StructPall),1,fp);
                fwrite(&TimingResults,sizeof(struct StructTimingResults),1,fp);
#ifdef COMPACT_IO_FORMAT //{
                for(int k=0;k<Nfirst;k++){
                    struct StructPbodyIOCompact PbodyIO = CopyPbodyToTemporalStructureCompact(k);
                    fwrite(&PbodyIO,sizeof(struct StructPbodyIOCompact),1,fp);
                    if(Pbody[k]->Type == TypeHydro){ // type hydro
                        struct StructPhydroIOCompact PhydroIO = 
                            CopyPhydroToTemporalStructureCompactElement(Pbody[k]->Baryon);
                        fwrite(&PhydroIO,sizeof(struct StructPhydroIOCompact),1,fp);
                    } else if (Pbody[k]->Type == TypeStar){ // type star
                        struct StructPstarIOCompact PstarIO = 
                            CopyPstarToTemporalStructureCompactElement(Pbody[k]->Baryon);
                        fwrite(&PstarIO,sizeof(struct StructPstarIOCompact),1,fp);
                    } else if (Pbody[k]->Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s:line %d\n",
                                MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
                        MPI_Finalize();
                        MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
                    }
                }

#elif defined(COMPACT_DOUBLE_IO_FORMAT)  //}//{
                for(int k=0;k<Nfirst;k++){
                    struct StructPbodyIOCompactDouble PbodyIO = CopyPbodyToTemporalStructureCompactDouble(k);
                    fwrite(&PbodyIO,sizeof(struct StructPbodyIOCompactDouble),1,fp);
                    if(Pbody[k]->Type == TypeHydro){ // type hydro
                        struct StructPhydroIOCompactDouble PhydroIO = 
                            CopyPhydroToTemporalStructureCompactDoubleElement(Pbody[k]->Baryon);
                        fwrite(&PhydroIO,sizeof(struct StructPhydroIOCompactDouble),1,fp);
                    } else if (Pbody[k]->Type == TypeStar){ // type star
                        struct StructPstarIOCompactDouble PstarIO = 
                            CopyPstarToTemporalStructureCompactDoubleElement(Pbody[k]->Baryon);
                        fwrite(&PstarIO,sizeof(struct StructPstarIOCompactDouble),1,fp);
                    } else if (Pbody[k]->Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s:line %d\n",
                                MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
                        MPI_Finalize();
                        MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
                    }
                }
#elif defined(COMPACT_MIX_IO_FORMAT)  //}//{
                for(int k=0;k<Nfirst;k++){
                    struct StructPbodyIOCompactMix PbodyIO = CopyPbodyToTemporalStructureCompactMix(k);
                    fwrite(&PbodyIO,sizeof(struct StructPbodyIOCompactMix),1,fp);
                    if(Pbody[k]->Type == TypeHydro){ // type hydro
                        struct StructPhydroIOCompactMix PhydroIO = 
                            CopyPhydroToTemporalStructureCompactMixElement(Pbody[k]->Baryon);
                        fwrite(&PhydroIO,sizeof(struct StructPhydroIOCompactMix),1,fp);
                    } else if (Pbody[k]->Type == TypeStar){ // type star
                        struct StructPstarIOCompactMix PstarIO = 
                            CopyPstarToTemporalStructureCompactMixElement(Pbody[k]->Baryon);
                        fwrite(&PstarIO,sizeof(struct StructPstarIOCompactMix),1,fp);
                    } else if (Pbody[k]->Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s:line %d\n",
                                MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
                        MPI_Finalize();
                        MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
                    }
                }
#else //COMPACT_IO_FORMAT //}
                for(int k=0;k<Nfirst;k++){
                    fwrite(Pbody[k],sizeof(StructPbody),1,fp);
                    if(Pbody[k]->Type == TypeHydro){ // type hydro
                        fwrite(Pbody[k]->Baryon,sizeof(StructPhydro),1,fp);
                    } else if (Pbody[k]->Type == TypeStar){ // type star
                        fwrite(Pbody[k]->Baryon,sizeof(StructPstar),1,fp);
                    } else if (Pbody[k]->Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                                MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
                        MPI_Finalize();
                        MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
                    }
                }
#endif
                if(GSL_EFAILED == gsl_rng_fwrite(fp,RandomGenerator)){
                    fprintf(stderr,"The random generator write error!\n");
                    MPI_Abort(MPI_COMM_WORLD,GSL_RANDOM_GENERATOR_WRITE_ERROR);
                }
                fclose(fp);

                // write second.
                sprintf(fname,"%s.%04d.%03d.%03d",Pall.BaseFileName,Pall.OutPutFileNumber,
                        2*MPIGetNumProcs(),2*MPIGetMyID()+1);
                FileOpen(fp,fname,"wb");

                // write header //
                WriteHeader(fp);
                // write header //
                fwrite(&TempPall_second,sizeof(struct StructPall),1,fp);
                fwrite(&TimingResults,sizeof(struct StructTimingResults),1,fp);
#ifdef COMPACT_IO_FORMAT //{
                for(int k=Nfirst;k<Pall.Ntotal;k++){
                    struct StructPbodyIOCompact PbodyIO = CopyPbodyToTemporalStructureCompact(k);
                    fwrite(&PbodyIO,sizeof(struct StructPbodyIOCompact),1,fp);
                    if(Pbody[k]->Type == TypeHydro){ // type hydro
                        struct StructPhydroIOCompact PhydroIO = 
                            CopyPhydroToTemporalStructureCompactElement(Pbody[k]->Baryon);
                        fwrite(&PhydroIO,sizeof(struct StructPhydroIOCompact),1,fp);
                    } else if (Pbody[k]->Type == TypeStar){ // type star
                        struct StructPstarIOCompact PstarIO = 
                            CopyPstarToTemporalStructureCompactElement(Pbody[k]->Baryon);
                        fwrite(&PstarIO,sizeof(struct StructPstarIOCompact),1,fp);
                    } else if (Pbody[k]->Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s:line %d\n",
                                MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
                        MPI_Finalize();
                        MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
                    }
                }
#elif defined(COMPACT_DOUBLE_IO_FORMAT)  //}//{
                for(int k=Nfirst;k<Pall.Ntotal;k++){
                    struct StructPbodyIOCompactDouble PbodyIO = CopyPbodyToTemporalStructureCompactDouble(k);
                    fwrite(&PbodyIO,sizeof(struct StructPbodyIOCompactDouble),1,fp);
                    if(Pbody[k]->Type == TypeHydro){ // type hydro
                        struct StructPhydroIOCompactDouble PhydroIO = 
                            CopyPhydroToTemporalStructureCompactDoubleElement(Pbody[k]->Baryon);
                        fwrite(&PhydroIO,sizeof(struct StructPhydroIOCompactDouble),1,fp);
                    } else if (Pbody[k]->Type == TypeStar){ // type star
                        struct StructPstarIOCompactDouble PstarIO = 
                            CopyPstarToTemporalStructureCompactDoubleElement(Pbody[k]->Baryon);
                        fwrite(&PstarIO,sizeof(struct StructPstarIOCompactDouble),1,fp);
                    } else if (Pbody[k]->Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s:line %d\n",
                                MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
                        MPI_Finalize();
                        MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
                    }
                }
#elif defined(COMPACT_MIX_IO_FORMAT)  //}//{
                for(int k=Nfirst;k<Pall.Ntotal;k++){
                    struct StructPbodyIOCompactMix PbodyIO = CopyPbodyToTemporalStructureCompactMix(k);
                    fwrite(&PbodyIO,sizeof(struct StructPbodyIOCompactMix),1,fp);
                    if(Pbody[k]->Type == TypeHydro){ // type hydro
                        struct StructPhydroIOCompactMix PhydroIO = 
                            CopyPhydroToTemporalStructureCompactMixElement(Pbody[k]->Baryon);
                        fwrite(&PhydroIO,sizeof(struct StructPhydroIOCompactMix),1,fp);
                    } else if (Pbody[k]->Type == TypeStar){ // type star
                        struct StructPstarIOCompactMix PstarIO = 
                            CopyPstarToTemporalStructureCompactMixElement(Pbody[k]->Baryon);
                        fwrite(&PstarIO,sizeof(struct StructPstarIOCompactMix),1,fp);
                    } else if (Pbody[k]->Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s:line %d\n",
                                MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
                        MPI_Finalize();
                        MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
                    }
                }
#else
                for(int k=Nfirst;k<Pall.Ntotal;k++){
                    fwrite(Pbody[k],sizeof(StructPbody),1,fp);
                    if(Pbody[k]->Type == TypeHydro){ // type hydro
                        fwrite(Pbody[k]->Baryon,sizeof(StructPhydro),1,fp);
                    } else if (Pbody[k]->Type == TypeStar){ // type star
                        fwrite(Pbody[k]->Baryon,sizeof(StructPstar),1,fp);
                    } else if (Pbody[k]->Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s:line %d\n",
                                MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
                        MPI_Finalize();
                        MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
                    }
                }
#endif
                if(GSL_EFAILED == gsl_rng_fwrite(fp,RandomGenerator)){
                    fprintf(stderr,"The random generator write error!\n");
                    MPI_Abort(MPI_COMM_WORLD,GSL_RANDOM_GENERATOR_WRITE_ERROR);
                }
                fclose(fp);
                fflush(NULL);

                fprintf(stderr,"End! %d\n",MPIGetMyID());
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    } else if(Pall.RunStatus == DataFileOperation_Halve){
        MPI_Status  mpi_status;
        struct  StructPall     TempPall;
        int SendRank,RecvRank;
        if(MPIGetMyID()%2 == 0){
            SendRank = MPIGetMyID() + 1;
            RecvRank = MPIGetMyID() + 1;
        } else {
            SendRank = MPIGetMyID() - 1;
            RecvRank = MPIGetMyID() - 1;
        }
        MPI_Sendrecv(&Pall,sizeof(struct StructPall),MPI_BYTE,SendRank,1,
                     &TempPall,sizeof(struct StructPall),MPI_BYTE,RecvRank,1,
                          MPI_COMM_WORLD,&mpi_status);

        for(int i=0;i<MPIGetNumProcs();i++){
            if(MPIGetMyID() == i){
                fprintf(stderr,"Start! %d\n",MPIGetMyID());
                struct  StructPall     TempPall_write = Pall;
                TempPall_write.Ntotal += TempPall.Ntotal;
                TempPall_write.NDM += TempPall.NDM;
                TempPall_write.Nhydro += TempPall.Nhydro;
                TempPall_write.Nstars += TempPall.Nstars;
            
                fprintf(stderr,"[%d] -f- %ld %ld %ld %ld\n",MPIGetMyID(),
                        TempPall_write.Ntotal,TempPall_write.NDM,TempPall_write.Nhydro,TempPall_write.Nstars);

                char fname[MaxCharactersInLine];

                if(MPIGetMyID()%2 == 0){ // write header.
                    sprintf(fname,"%s.%04d.%03d.%03d",Pall.BaseFileName,Pall.OutPutFileNumber,
                            MPIGetNumProcs()/2,MPIGetMyID()/2);
                    FileOpen(fp,fname,"wb");

                    // write header //
                    WriteHeader(fp);
                    // write header //
                    fwrite(&TempPall_write,sizeof(struct StructPall),1,fp);
                    fwrite(&TimingResults,sizeof(struct StructTimingResults),1,fp);
                } else { // reopen.
                    sprintf(fname,"%s.%04d.%03d.%03d",Pall.BaseFileName,Pall.OutPutFileNumber,
                            MPIGetNumProcs()/2,MPIGetMyID()/2);
                    FileOpen(fp,fname,"ab");
                }
#ifdef COMPACT_IO_FORMAT
                for(int k=0;k<Pall.Ntotal;k++){
                    struct StructPbodyIOCompact PbodyIO = CopyPbodyToTemporalStructureCompact(k);
                    fwrite(&PbodyIO,sizeof(struct StructPbodyIOCompact),1,fp);
                    if(Pbody[k]->Type == TypeHydro){ // type hydro
                        struct StructPhydroIOCompact PhydroIO = 
                            CopyPhydroToTemporalStructureCompactElement(Pbody[k]->Baryon);
                        fwrite(&PhydroIO,sizeof(struct StructPhydroIOCompact),1,fp);
                    } else if (Pbody[k]->Type == TypeStar){ // type star
                        struct StructPstarIOCompact PstarIO = 
                            CopyPstarToTemporalStructureCompactElement(Pbody[k]->Baryon);
                        fwrite(&PstarIO,sizeof(struct StructPstarIOCompact),1,fp);
                    } else if (Pbody[k]->Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s:line %d\n",
                                MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
                        MPI_Finalize();
                        MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
                    }
                }
#elif defined(COMPACT_DOUBLE_IO_FORMAT)
                for(int k=0;k<Pall.Ntotal;k++){
                    struct StructPbodyIOCompactDouble PbodyIO = CopyPbodyToTemporalStructureCompactDouble(k);
                    fwrite(&PbodyIO,sizeof(struct StructPbodyIOCompactDouble),1,fp);
                    if(Pbody[k]->Type == TypeHydro){ // type hydro
                        struct StructPhydroIOCompactDouble PhydroIO = 
                            CopyPhydroToTemporalStructureCompactDoubleElement(Pbody[k]->Baryon);
                        fwrite(&PhydroIO,sizeof(struct StructPhydroIOCompactDouble),1,fp);
                    } else if (Pbody[k]->Type == TypeStar){ // type star
                        struct StructPstarIOCompactDouble PstarIO = 
                            CopyPstarToTemporalStructureCompactDoubleElement(Pbody[k]->Baryon);
                        fwrite(&PstarIO,sizeof(struct StructPstarIOCompactDouble),1,fp);
                    } else if (Pbody[k]->Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s:line %d\n",
                                MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
                        MPI_Finalize();
                        MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
                    }
                }
#elif defined(COMPACT_MIX_IO_FORMAT)
                for(int k=0;k<Pall.Ntotal;k++){
                    struct StructPbodyIOCompactMix PbodyIO = CopyPbodyToTemporalStructureCompactMix(k);
                    fwrite(&PbodyIO,sizeof(struct StructPbodyIOCompactMix),1,fp);
                    if(Pbody[k]->Type == TypeHydro){ // type hydro
                        struct StructPhydroIOCompactMix PhydroIO = 
                            CopyPhydroToTemporalStructureCompactMixElement(Pbody[k]->Baryon);
                        fwrite(&PhydroIO,sizeof(struct StructPhydroIOCompactMix),1,fp);
                    } else if (Pbody[k]->Type == TypeStar){ // type star
                        struct StructPstarIOCompactMix PstarIO = 
                            CopyPstarToTemporalStructureCompactMixElement(Pbody[k]->Baryon);
                        fwrite(&PstarIO,sizeof(struct StructPstarIOCompactMix),1,fp);
                    } else if (Pbody[k]->Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s:line %d\n",
                                MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
                        MPI_Finalize();
                        MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
                    }
                }
#else
                for(int k=0;k<Pall.Ntotal;k++){
                    fwrite(Pbody[k],sizeof(StructPbody),1,fp);
                    if(Pbody[k]->Type == TypeHydro){ // type hydro
                        fwrite(Pbody[k]->Baryon,sizeof(StructPhydro),1,fp);
                    } else if (Pbody[k]->Type == TypeStar){ // type star
                        fwrite(Pbody[k]->Baryon,sizeof(StructPstar),1,fp);
                    } else if (Pbody[k]->Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s:line %d\n",
                                MPIGetMyID(),Pbody[k]->Type,__FILE__,__LINE__);
                        MPI_Finalize();
                        MPI_Abort(MPI_COMM_WORLD,WriteTypeUnDefined);
                    }
                }
#endif
                if(MPIGetMyID()%2 == 1){ // write random generator.
                    if(GSL_EFAILED == gsl_rng_fwrite(fp,RandomGenerator)){
                        fprintf(stderr,"The random generator write error!\n");
                        MPI_Abort(MPI_COMM_WORLD,GSL_RANDOM_GENERATOR_WRITE_ERROR);
                    }
                }
                fclose(fp);
                fflush(NULL);
                fprintf(stderr,"End! %d\n",MPIGetMyID());
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"End data write!\n");

    return;
}

/*
 * This function merges all of the restart files into a single file.  Only the
 * first MPI process do the operation and others are sleeping during the
 * operation.  After the opration, all processes are terminated.
 */
static void MergerAllData(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int HeaderArray[HEADER_SIZE];
    int RunStatusOperation = Pall.RunStatus;

    // Run only the root MPI process.
    if(MyID==MPI_ROOT_RANK){
        // Open write file and read files.
        FILE *fp_w;
        char fname_w[MaxCharactersInLine];

        sprintf(fname_w,"%s.%03d.%03d",Pall.RestartFileName,0,0);
        FileOpen(fp_w,fname_w,"wb");

        for(int i=0;i<NProcs;i++){
            FILE *fp_r;
            char fname_r[MaxCharactersInLine];
            sprintf(fname_r,"%s.%03d.%03d",Pall.RestartFileName,NProcs,i);
            FileOpen(fp_r,fname_r,"rb");
            ReadHeader(fp_r,HeaderArray);
            fread(&Pall,sizeof(struct StructPall),1,fp_r);
            Pall.RunStatus = RunStatusOperation;
            fread(&TimingResults,sizeof(struct StructTimingResults),1,fp_r);

            if(i==0){ // Write header.
                struct  StructPall     TempPall = Pall;
                TempPall.Ntotal = Pall.Ntotal_t;
                TempPall.NDM = Pall.NDM_t;
                TempPall.Nhydro = Pall.Nhydro_t;
                TempPall.Nstars = Pall.Nstars_t;

                WriteHeader(fp_w);
                fwrite(&TempPall,sizeof(struct StructPall),1,fp_w);
                fwrite(&TimingResults,sizeof(struct StructTimingResults),1,fp_w);
            }

            // Allocate structures.
            GenerateStructPbody(Pall.Ntotal);
            if(Pall.Nhydro_t>0) GenerateStructPhydro(Pall.Nhydro);
            if(Pall.Nstars_t>0) GenerateStructPstar(Pall.Nstars);

            int CountTypes[NTypes];
            for(int k=0;k<NTypes;k++)
                CountTypes[i] = 0;

            if(HeaderArray[HD_CompactFormatFlag] == OFF){
                for(int k=0;k<Pall.Ntotal;k++){
                    fread(PbodyElements+k,sizeof(StructPbody),1,fp_r);
                    fwrite(PbodyElements+k,sizeof(StructPbody),1,fp_w); // write
                    if(PbodyElements[k].Type == TypeHydro){ // type hydro
                        fread(PhydroElements+CountTypes[TypeHydro],sizeof(StructPhydro),1,fp_r);
                        fwrite(PhydroElements+CountTypes[TypeHydro],sizeof(StructPhydro),1,fp_w); // write
                        CountTypes[TypeHydro] ++;
                    } else if (PbodyElements[k].Type == TypeStar){ // type star
                        fread(PstarElements+CountTypes[TypeStar],sizeof(StructPstar),1,fp_r);
                        fwrite(PstarElements+CountTypes[TypeStar],sizeof(StructPstar),1,fp_w); // write
                        CountTypes[TypeStar] ++;
                    } else if (PbodyElements[k].Type == TypeDM){
                        CountTypes[TypeDM] ++;
                    } else if (PbodyElements[k].Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s\n",
                                MPIGetMyID(),PbodyElements[k].Type,__FILE__);
                        MPI_Finalize();
                        exit(ReadTypeUnDefined);
                    }
                }
            } else if(HeaderArray[HD_CompactFormatFlag] == ON){
                for(int k=0;k<Pall.Ntotal;k++){
                    StructPbodyptr Pb = PbodyElements[k].Next;
                    struct StructPbodyIOCompact   PbodyIO;
                    fread(&PbodyIO,sizeof(struct StructPbodyIOCompact),1,fp_r);
                    fwrite(&PbodyIO,sizeof(struct StructPbodyIOCompact),1,fp_w); // write
                    if(PbodyIO.Type == TypeHydro){ // type hydro
                        struct StructPhydroIOCompact  PhydroIO;
                        fread(&PhydroIO,sizeof(struct StructPhydroIOCompact),1,fp_r);
                        fwrite(&PhydroIO,sizeof(struct StructPhydroIOCompact),1,fp_w); // write
                        CountTypes[TypeHydro] ++;
                    } else if (PbodyIO.Type == TypeStar){ // type hydro
                        struct StructPstarIOCompact   PstarIO;
                        fread(&PstarIO,sizeof(struct StructPstarIOCompact),1,fp_r);
                        fwrite(&PstarIO,sizeof(struct StructPstarIOCompact),1,fp_w); // write
                        CountTypes[TypeStar] ++;
                    } else if (PbodyElements[k].Type == TypeDM){
                        CountTypes[TypeDM] ++;
                    } else if (PbodyElements[k].Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s\n",
                                MPIGetMyID(),PbodyElements[k].Type,__FILE__);
                        MPI_Finalize();
                        exit(ReadTypeUnDefined);
                    }
                }
            } else if(HeaderArray[HD_CompactDoubleFormatFlag] == ON){
                for(int k=0;k<Pall.Ntotal;k++){
                    StructPbodyptr Pb = PbodyElements[k].Next;
                    struct StructPbodyIOCompactDouble   PbodyIO;
                    fread(&PbodyIO,sizeof(struct StructPbodyIOCompactDouble),1,fp_r);
                    fwrite(&PbodyIO,sizeof(struct StructPbodyIOCompactDouble),1,fp_w); // write
                    if(PbodyIO.Type == TypeHydro){ // type hydro
                        struct StructPhydroIOCompactDouble  PhydroIO;
                        fread(&PhydroIO,sizeof(struct StructPhydroIOCompactDouble),1,fp_r);
                        fwrite(&PhydroIO,sizeof(struct StructPhydroIOCompactDouble),1,fp_w); // write
                        CountTypes[TypeHydro] ++;
                    } else if (PbodyIO.Type == TypeStar){ // type hydro
                        struct StructPstarIOCompactDouble   PstarIO;
                        fread(&PstarIO,sizeof(struct StructPstarIOCompactDouble),1,fp_r);
                        fwrite(&PstarIO,sizeof(struct StructPstarIOCompactDouble),1,fp_w); // write
                        CountTypes[TypeStar] ++;
                    } else if (PbodyElements[k].Type == TypeDM){
                        CountTypes[TypeDM] ++;
                    } else if (PbodyElements[k].Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s\n",
                                MPIGetMyID(),PbodyElements[k].Type,__FILE__);
                        MPI_Finalize();
                        exit(ReadTypeUnDefined);
                    }
                }
            } else if(HeaderArray[HD_CompactMixFormatFlag] == ON){
                for(int k=0;k<Pall.Ntotal;k++){
                    StructPbodyptr Pb = PbodyElements[k].Next;
                    struct StructPbodyIOCompactMix   PbodyIO;
                    fread(&PbodyIO,sizeof(struct StructPbodyIOCompactMix),1,fp_r);
                    fwrite(&PbodyIO,sizeof(struct StructPbodyIOCompactMix),1,fp_w); // write
                    if(PbodyIO.Type == TypeHydro){ // type hydro
                        struct StructPhydroIOCompactMix  PhydroIO;
                        fread(&PhydroIO,sizeof(struct StructPhydroIOCompactMix),1,fp_r);
                        fwrite(&PhydroIO,sizeof(struct StructPhydroIOCompactMix),1,fp_w); // write
                        CountTypes[TypeHydro] ++;
                    } else if (PbodyIO.Type == TypeStar){ // type hydro
                        struct StructPstarIOCompactMix   PstarIO;
                        fread(&PstarIO,sizeof(struct StructPstarIOCompactMix),1,fp_r);
                        fwrite(&PstarIO,sizeof(struct StructPstarIOCompactMix),1,fp_w); // write
                        CountTypes[TypeStar] ++;
                    } else if (PbodyElements[k].Type == TypeDM){
                        CountTypes[TypeDM] ++;
                    } else if (PbodyElements[k].Type != TypeDM){
                        fprintf(stderr,"[%02d] Type (%d) Undefined in %s\n",
                                MPIGetMyID(),PbodyElements[k].Type,__FILE__);
                        MPI_Finalize();
                        exit(ReadTypeUnDefined);
                    }
                }
            }else{
                fprintf(stderr,"Header flag error!\n");
                MPI_Abort(MPI_COMM_WORLD,UnDefinedParticleType);
                exit(HeaderError);
            }

            // Release structures.
            FreeStructPbody();
            FreeStructPhydro();
            FreeStructPstar();

            if(i==NProcs-1){ // Write footer.
                if(GSL_EFAILED == gsl_rng_fwrite(fp_w,RandomGenerator)){
                    fprintf(stderr,"The random generator write error!\n");
                    MPI_Abort(MPI_COMM_WORLD,GSL_RANDOM_GENERATOR_WRITE_ERROR);
                }
            }
            fclose(fp_r);
        }
        fclose(fp_w);

    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    fprintf(stderr,"End %s!\n",__FUNCTION__);
    exit(EXIT_SUCCESS);

    return;
}

static int ReturnTotalFileNumber(void){

    int TotalFileNumber = 1;
    char fname[MaxCharactersInLine];
    do {
        Snprintf(fname,"%s.%03d.%03d",Pall.RestartFileName,TotalFileNumber,0);
        if(CheckFile(fname))
            return TotalFileNumber;
        //TotalFileNumber *= 2;
        TotalFileNumber += 1;
    } while(TotalFileNumber < 8192);

    return NONE;

}

static int ReadParallelDataOnSingleNodeFull(FILE *fp, const int Ntotal, int counter, int CountTypes[]){

    for(int k=0;k<Ntotal;k++){ 
        StructPbodyptr Pb = PbodyElements[counter].Next;
        fread(PbodyElements+counter,sizeof(StructPbody),1,fp);
        PbodyElements[counter].Next = Pb;

        if(PbodyElements[counter].Type == TypeHydro){ // type hydro
            StructPhydroptr Ph = PhydroElements[CountTypes[TypeHydro]].Next;
            fread(PhydroElements+CountTypes[TypeHydro],sizeof(StructPhydro),1,fp);
            PhydroElements[CountTypes[TypeHydro]].Next = Ph;
            PhydroElements[CountTypes[TypeHydro]].Body = PbodyElements+counter;
            PbodyElements[counter].Baryon = (void*)(PhydroElements+CountTypes[TypeHydro]);
            CountTypes[TypeHydro] ++;
        } else if (PbodyElements[counter].Type == TypeStar){ // type star
            StructPstarptr Ps = PstarElements[CountTypes[TypeStar]].Next;
            fread(PstarElements+CountTypes[TypeStar],sizeof(StructPstar),1,fp);
            PstarElements[CountTypes[TypeStar]].Next = Ps;
            PstarElements[CountTypes[TypeStar]].Body = PbodyElements+counter;
            PbodyElements[counter].Baryon = (void*)(PstarElements+CountTypes[TypeStar]);
            CountTypes[TypeStar] ++;
        } else if (PbodyElements[counter].Type == TypeSink){ // type sink
            StructPsinkptr Psk = PsinkElements[CountTypes[TypeSink]].Next;
            fread(PsinkElements+CountTypes[TypeSink],sizeof(StructPsink),1,fp);
            PsinkElements[CountTypes[TypeSink]].Next = Psk;
            PsinkElements[CountTypes[TypeSink]].Body = PbodyElements+counter;
            PbodyElements[counter].Baryon = (void*)(PsinkElements+CountTypes[TypeSink]);
            CountTypes[TypeSink] ++;
        } else if (PbodyElements[counter].Type == TypeDM){
            CountTypes[TypeDM] ++;
        } else if (PbodyElements[counter].Type != TypeDM){
            fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                    MPIGetMyID(),PbodyElements[k].Type,__FILE__,__LINE__);
            MPI_Finalize();
            exit(ReadTypeUnDefined);
        }
        counter ++;
    }

    return counter;
}

static int ReadParallelDataOnSingleNodeCompact(FILE *fp, const int Ntotal, int counter, int CountTypes[]){

    for(int k=0;k<Ntotal;k++){
        StructPbodyptr Pb = PbodyElements[counter].Next;
        struct StructPbodyIOCompact   PbodyIO;
        fread(&PbodyIO,sizeof(struct StructPbodyIOCompact),1,fp);

        PbodyElements[counter] = CopyTemporalStructureCompactToPbodyCompact(PbodyIO);
        PbodyElements[counter].Use = ON;
        PbodyElements[counter].Next = Pb;
        /* predictor copy */
#if 1
        for(int l=0;l<3;l++)
            PbodyElements[counter].PosP[l] = PbodyElements[counter].Pos[l];
#endif
        /* predictor copy */

        if(PbodyElements[counter].Type == TypeHydro){ // type hydro
            StructPhydroptr Ph = PhydroElements[CountTypes[TypeHydro]].Next;
            struct StructPhydroIOCompact  PhydroIO;
            fread(&PhydroIO,sizeof(struct StructPhydroIOCompact),1,fp);

            PhydroElements[CountTypes[TypeHydro]] = CopyTemporalStructureCompactToPhydroCompact(PhydroIO);
            PhydroElements[CountTypes[TypeHydro]].Use = ON;
            PhydroElements[CountTypes[TypeHydro]].Next = Ph;
            PhydroElements[CountTypes[TypeHydro]].Body = PbodyElements+counter;
            //PbodyElements[k].Baryon = (void*)(PhydroElements+CountTypes[TypeHydro]);
            PbodyElements[counter].Baryon = (void*)(PhydroElements+CountTypes[TypeHydro]);
            /* predictor copy */
            PhydroElements[CountTypes[TypeHydro]].Mass = PbodyElements[counter].Mass;
#if 1
            for(int l=0;l<3;l++){
                PhydroElements[CountTypes[TypeHydro]].PosP[l] = PbodyElements[counter].Pos[l];
                PhydroElements[CountTypes[TypeHydro]].VelP[l] = PbodyElements[counter].Vel[l];
            }
            PhydroElements[CountTypes[TypeHydro]].RhoPred = PhydroElements[CountTypes[TypeHydro]].Rho;
            PhydroElements[CountTypes[TypeHydro]].KernelPred = PhydroElements[CountTypes[TypeHydro]].Kernel;
            PhydroElements[CountTypes[TypeHydro]].UPred = PhydroElements[CountTypes[TypeHydro]].U;
            PhydroElements[CountTypes[TypeHydro]].Active = PhydroElements[CountTypes[TypeHydro]].Active;
#endif
            /* predictor copy */
            CountTypes[TypeHydro] ++;
        } else if (PbodyElements[counter].Type == TypeStar){ // type star
            StructPstarptr Ps = PstarElements[CountTypes[TypeStar]].Next;
            struct StructPstarIOCompact   PstarIO;
            fread(&PstarIO,sizeof(struct StructPstarIOCompact),1,fp);

            PstarElements[CountTypes[TypeStar]] = CopyTemporalStructureCompactToPstarCompact(PstarIO);
            PstarElements[CountTypes[TypeStar]].Use = ON;
            PstarElements[CountTypes[TypeStar]].Next = Ps;
            PstarElements[CountTypes[TypeStar]].Body = PbodyElements+counter;
            //PbodyElements[k].Baryon = (void*)(PstarElements+CountTypes[TypeStar]);
            PbodyElements[counter].Baryon = (void*)(PstarElements+CountTypes[TypeStar]);
            CountTypes[TypeStar] ++;
        } else if (PbodyElements[counter].Type == TypeSink){ // type sink
            StructPsinkptr Psk = PsinkElements[CountTypes[TypeSink]].Next;
            struct StructPsinkIOCompact   PsinkIO;
            fread(&PsinkIO,sizeof(struct StructPsinkIOCompact),1,fp);

            PsinkElements[CountTypes[TypeSink]] = CopyTemporalStructureCompactToPsinkCompact(PsinkIO);
            PsinkElements[CountTypes[TypeSink]].Use = ON;
            PsinkElements[CountTypes[TypeSink]].Next = Psk;
            PsinkElements[CountTypes[TypeSink]].Body = PbodyElements+counter;
            PbodyElements[counter].Baryon = (void*)(PsinkElements+CountTypes[TypeSink]);
            CountTypes[TypeSink] ++;
        } else if (PbodyElements[counter].Type == TypeDM){
            CountTypes[TypeDM] ++;
        } else if (PbodyElements[counter].Type != TypeDM){
            fprintf(stderr,"[%02d] %d %d %d\n",
                    MPIGetMyID(),CountTypes[TypeDM],CountTypes[TypeHydro],CountTypes[TypeStar]);
            fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                    MPIGetMyID(),PbodyElements[counter].Type,__FILE__,__LINE__);
            MPI_Finalize();
            exit(ReadTypeUnDefined);
        }
        counter ++;
    }

    return counter;
}

static int ReadParallelDataOnSingleNodeCompactDouble(FILE *fp, const int Ntotal, int counter, int CountTypes[]){

    for(int k=0;k<Ntotal;k++){
        StructPbodyptr Pb = PbodyElements[counter].Next;
        struct StructPbodyIOCompactDouble   PbodyIO;
        fread(&PbodyIO,sizeof(struct StructPbodyIOCompactDouble),1,fp);

        PbodyElements[counter] = CopyTemporalStructureCompactDoubleToPbodyCompactDouble(PbodyIO);
        PbodyElements[counter].Use = ON;
        PbodyElements[counter].Next = Pb;
        /* predictor copy */
#if 1
        for(int l=0;l<3;l++)
            PbodyElements[counter].PosP[l] = PbodyElements[counter].Pos[l];
#endif
        /* predictor copy */

        if(PbodyElements[counter].Type == TypeHydro){ // type hydro
            StructPhydroptr Ph = PhydroElements[CountTypes[TypeHydro]].Next;
            struct StructPhydroIOCompactDouble  PhydroIO;
            fread(&PhydroIO,sizeof(struct StructPhydroIOCompactDouble),1,fp);

            PhydroElements[CountTypes[TypeHydro]] = CopyTemporalStructureCompactDoubleToPhydroCompactDouble(PhydroIO);
            PhydroElements[CountTypes[TypeHydro]].Use = ON;
            PhydroElements[CountTypes[TypeHydro]].Next = Ph;
            PhydroElements[CountTypes[TypeHydro]].Body = PbodyElements+counter;
            //PbodyElements[k].Baryon = (void*)(PhydroElements+CountTypes[TypeHydro]);
            PbodyElements[counter].Baryon = (void*)(PhydroElements+CountTypes[TypeHydro]);
            /* predictor copy */
            PhydroElements[CountTypes[TypeHydro]].Mass = PbodyElements[counter].Mass;
#if 1
            for(int l=0;l<3;l++){
                PhydroElements[CountTypes[TypeHydro]].PosP[l] = PbodyElements[counter].Pos[l];
                PhydroElements[CountTypes[TypeHydro]].VelP[l] = PbodyElements[counter].Vel[l];
            }
            PhydroElements[CountTypes[TypeHydro]].RhoPred = PhydroElements[CountTypes[TypeHydro]].Rho;
            PhydroElements[CountTypes[TypeHydro]].KernelPred = PhydroElements[CountTypes[TypeHydro]].Kernel;
            PhydroElements[CountTypes[TypeHydro]].UPred = PhydroElements[CountTypes[TypeHydro]].U;
            PhydroElements[CountTypes[TypeHydro]].Active = PhydroElements[CountTypes[TypeHydro]].Active;
#endif
            /* predictor copy */
            CountTypes[TypeHydro] ++;
        } else if (PbodyElements[counter].Type == TypeStar){ // type star
            StructPstarptr Ps = PstarElements[CountTypes[TypeStar]].Next;
            struct StructPstarIOCompactDouble   PstarIO;
            fread(&PstarIO,sizeof(struct StructPstarIOCompactDouble),1,fp);

            PstarElements[CountTypes[TypeStar]] = CopyTemporalStructureCompactDoubleToPstarCompactDouble(PstarIO);
            PstarElements[CountTypes[TypeStar]].Use = ON;
            PstarElements[CountTypes[TypeStar]].Next = Ps;
            PstarElements[CountTypes[TypeStar]].Body = PbodyElements+counter;
            //PbodyElements[k].Baryon = (void*)(PstarElements+CountTypes[TypeStar]);
            PbodyElements[counter].Baryon = (void*)(PstarElements+CountTypes[TypeStar]);
            CountTypes[TypeStar] ++;
        } else if (PbodyElements[counter].Type == TypeSink){ // type sink
            StructPsinkptr Psk = PsinkElements[CountTypes[TypeSink]].Next;
            struct StructPsinkIOCompactDouble   PsinkIO;
            fread(&PsinkIO,sizeof(struct StructPsinkIOCompactDouble),1,fp);

            PsinkElements[CountTypes[TypeSink]] = CopyTemporalStructureCompactDoubleToPsinkCompactDouble(PsinkIO);
            PsinkElements[CountTypes[TypeSink]].Use = ON;
            PsinkElements[CountTypes[TypeSink]].Next = Psk;
            PsinkElements[CountTypes[TypeSink]].Body = PbodyElements+counter;
            PbodyElements[counter].Baryon = (void*)(PsinkElements+CountTypes[TypeSink]);
            CountTypes[TypeSink] ++;
        } else if (PbodyElements[counter].Type == TypeDM){
            CountTypes[TypeDM] ++;
        } else if (PbodyElements[counter].Type != TypeDM){
            fprintf(stderr,"[%02d] %d %d %d\n",
                    MPIGetMyID(),CountTypes[TypeDM],CountTypes[TypeHydro],CountTypes[TypeStar]);
            fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                    MPIGetMyID(),PbodyElements[counter].Type,__FILE__,__LINE__);
            MPI_Finalize();
            exit(ReadTypeUnDefined);
        }
        counter ++;
    }

    return counter;
}

static int ReadParallelDataOnSingleNodeCompactMix(FILE *fp, const int Ntotal, int counter, int CountTypes[]){

    for(int k=0;k<Ntotal;k++){
        StructPbodyptr Pb = PbodyElements[counter].Next;
        struct StructPbodyIOCompactMix   PbodyIO;
        fread(&PbodyIO,sizeof(struct StructPbodyIOCompactMix),1,fp);

        PbodyElements[counter] = CopyTemporalStructureCompactMixToPbodyCompactMix(PbodyIO);
        PbodyElements[counter].Use = ON;
        PbodyElements[counter].Next = Pb;
        /* predictor copy */
#if 1
        for(int l=0;l<3;l++)
            PbodyElements[counter].PosP[l] = PbodyElements[counter].Pos[l];
#endif
        /* predictor copy */

        if(PbodyElements[counter].Type == TypeHydro){ // type hydro
            StructPhydroptr Ph = PhydroElements[CountTypes[TypeHydro]].Next;
            struct StructPhydroIOCompactMix  PhydroIO;
            fread(&PhydroIO,sizeof(struct StructPhydroIOCompactMix),1,fp);

            PhydroElements[CountTypes[TypeHydro]] = CopyTemporalStructureCompactMixToPhydroCompactMix(PhydroIO);
            PhydroElements[CountTypes[TypeHydro]].Use = ON;
            PhydroElements[CountTypes[TypeHydro]].Next = Ph;
            PhydroElements[CountTypes[TypeHydro]].Body = PbodyElements+counter;
            //PbodyElements[k].Baryon = (void*)(PhydroElements+CountTypes[TypeHydro]);
            PbodyElements[counter].Baryon = (void*)(PhydroElements+CountTypes[TypeHydro]);
            /* predictor copy */
            PhydroElements[CountTypes[TypeHydro]].Mass = PbodyElements[counter].Mass;
#if 1
            for(int l=0;l<3;l++){
                PhydroElements[CountTypes[TypeHydro]].PosP[l] = PbodyElements[counter].Pos[l];
                PhydroElements[CountTypes[TypeHydro]].VelP[l] = PbodyElements[counter].Vel[l];
            }
            PhydroElements[CountTypes[TypeHydro]].RhoPred = PhydroElements[CountTypes[TypeHydro]].Rho;
            PhydroElements[CountTypes[TypeHydro]].KernelPred = PhydroElements[CountTypes[TypeHydro]].Kernel;
            PhydroElements[CountTypes[TypeHydro]].UPred = PhydroElements[CountTypes[TypeHydro]].U;
            PhydroElements[CountTypes[TypeHydro]].Active = PhydroElements[CountTypes[TypeHydro]].Active;
#endif
            /* predictor copy */
            CountTypes[TypeHydro] ++;
        } else if (PbodyElements[counter].Type == TypeStar){ // type star
            StructPstarptr Ps = PstarElements[CountTypes[TypeStar]].Next;
            struct StructPstarIOCompactMix   PstarIO;
            fread(&PstarIO,sizeof(struct StructPstarIOCompactMix),1,fp);

            PstarElements[CountTypes[TypeStar]] = CopyTemporalStructureCompactMixToPstarCompactMix(PstarIO);
            PstarElements[CountTypes[TypeStar]].Use = ON;
            PstarElements[CountTypes[TypeStar]].Next = Ps;
            PstarElements[CountTypes[TypeStar]].Body = PbodyElements+counter;
            //PbodyElements[k].Baryon = (void*)(PstarElements+CountTypes[TypeStar]);
            PbodyElements[counter].Baryon = (void*)(PstarElements+CountTypes[TypeStar]);
            CountTypes[TypeStar] ++;
        } else if (PbodyElements[counter].Type == TypeSink){ // type sink
            StructPsinkptr Psk = PsinkElements[CountTypes[TypeSink]].Next;
            struct StructPsinkIOCompactMix   PsinkIO;
            fread(&PsinkIO,sizeof(struct StructPsinkIOCompactMix),1,fp);

            PsinkElements[CountTypes[TypeSink]] = CopyTemporalStructureCompactMixToPsinkCompactMix(PsinkIO);
            PsinkElements[CountTypes[TypeSink]].Use = ON;
            PsinkElements[CountTypes[TypeSink]].Next = Psk;
            PsinkElements[CountTypes[TypeSink]].Body = PbodyElements+counter;
            PbodyElements[counter].Baryon = (void*)(PsinkElements+CountTypes[TypeSink]);
            CountTypes[TypeSink] ++;
        } else if (PbodyElements[counter].Type == TypeDM){
            CountTypes[TypeDM] ++;
        } else if (PbodyElements[counter].Type != TypeDM){
            fprintf(stderr,"[%02d] %d %d %d\n",
                    MPIGetMyID(),CountTypes[TypeDM],CountTypes[TypeHydro],CountTypes[TypeStar]);
            fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                    MPIGetMyID(),PbodyElements[counter].Type,__FILE__,__LINE__);
            MPI_Finalize();
            exit(ReadTypeUnDefined);
        }
        counter ++;
    }

    return counter;
}

static int ReadParallelDataOnSingleNodeLean(FILE *fp, const int Ntotal, int counter, int CountTypes[]){

    for(int k=0;k<Ntotal;k++){
        StructPbodyptr Pb = PbodyElements[counter].Next;
        struct StructPbodyIOLean   PbodyIO;
        fread(&PbodyIO,sizeof(struct StructPbodyIOLean),1,fp);

        PbodyElements[counter] = CopyTemporalStructureLeanToPbodyLean(PbodyIO);
        PbodyElements[counter].Use = ON;
        PbodyElements[counter].Next = Pb;
        /* predictor copy */
#if 1
        for(int l=0;l<3;l++)
            PbodyElements[counter].PosP[l] = PbodyElements[counter].Pos[l];
#endif
        /* predictor copy */

        if(PbodyElements[counter].Type == TypeHydro){ // type hydro
            StructPhydroptr Ph = PhydroElements[CountTypes[TypeHydro]].Next;
            struct StructPhydroIOLean  PhydroIO;
            fread(&PhydroIO,sizeof(struct StructPhydroIOLean),1,fp);

            PhydroElements[CountTypes[TypeHydro]] = CopyTemporalStructureLeanToPhydroLean(PhydroIO);
            PhydroElements[CountTypes[TypeHydro]].Use = ON;
            PhydroElements[CountTypes[TypeHydro]].Next = Ph;
            PhydroElements[CountTypes[TypeHydro]].Body = PbodyElements+counter;
            PbodyElements[counter].Baryon = (void*)(PhydroElements+CountTypes[TypeHydro]);
            /* predictor copy */
            PhydroElements[CountTypes[TypeHydro]].Mass = PbodyElements[counter].Mass;
#if 1
            for(int l=0;l<3;l++){
                PhydroElements[CountTypes[TypeHydro]].PosP[l] = PbodyElements[counter].Pos[l];
                PhydroElements[CountTypes[TypeHydro]].VelP[l] = PbodyElements[counter].Vel[l];
            }
            PhydroElements[CountTypes[TypeHydro]].RhoPred = PhydroElements[CountTypes[TypeHydro]].Rho;
            PhydroElements[CountTypes[TypeHydro]].KernelPred = PhydroElements[CountTypes[TypeHydro]].Kernel;
            PhydroElements[CountTypes[TypeHydro]].UPred = PhydroElements[CountTypes[TypeHydro]].U;
            PhydroElements[CountTypes[TypeHydro]].Active = PhydroElements[CountTypes[TypeHydro]].Active;
#endif
            /* predictor copy */
            CountTypes[TypeHydro] ++;
        } else if (PbodyElements[counter].Type == TypeStar){ // type star
            StructPstarptr Ps = PstarElements[CountTypes[TypeStar]].Next;
            struct StructPstarIOLean   PstarIO;
            fread(&PstarIO,sizeof(struct StructPstarIOLean),1,fp);

            PstarElements[CountTypes[TypeStar]] = CopyTemporalStructureLeanToPstarLean(PstarIO);
            PstarElements[CountTypes[TypeStar]].Use = ON;
            PstarElements[CountTypes[TypeStar]].Next = Ps;
            PstarElements[CountTypes[TypeStar]].Body = PbodyElements+counter;
            PbodyElements[counter].Baryon = (void*)(PstarElements+CountTypes[TypeStar]);
            CountTypes[TypeStar] ++;
        } else if (PbodyElements[counter].Type == TypeSink){ // type sink
            StructPsinkptr Psk = PsinkElements[CountTypes[TypeSink]].Next;
            struct StructPsinkIOLean   PsinkIO;
            fread(&PsinkIO,sizeof(struct StructPsinkIOLean),1,fp);

            PsinkElements[CountTypes[TypeSink]] = CopyTemporalStructureLeanToPsinkLean(PsinkIO);
            PsinkElements[CountTypes[TypeSink]].Use = ON;
            PsinkElements[CountTypes[TypeSink]].Next = Psk;
            PsinkElements[CountTypes[TypeSink]].Body = PbodyElements+counter;
            PbodyElements[counter].Baryon = (void*)(PsinkElements+CountTypes[TypeSink]);
            CountTypes[TypeSink] ++;
        } else if (PbodyElements[counter].Type == TypeDM){
            CountTypes[TypeDM] ++;
        } else if (PbodyElements[counter].Type != TypeDM){
            fprintf(stderr,"[%02d] %d %d %d\n",
                    MPIGetMyID(),CountTypes[TypeDM],CountTypes[TypeHydro],CountTypes[TypeStar]);
            fprintf(stderr,"[%02d] Type (%d) Undefined in %s:%d\n",
                    MPIGetMyID(),PbodyElements[counter].Type,__FILE__,__LINE__);
            MPI_Finalize();
            exit(ReadTypeUnDefined);
        }
        counter ++;
    }

    return counter;
}

/*
 * This function reads all date on a sinlge node. This is used for the deta
 * analyses.
 */

void ReadParallelDataOnSingleNode(void){

    FILE *fp;
    char fname[MaxCharactersInLine];
    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int HeaderArray[HEADER_SIZE];

    int TotalFileNumber = ReturnTotalFileNumber();

    if(NProcs > 1){
        fprintf(stderr,"Serial run only\n");
        exit(EXIT_FAILURE);
    }

    char FileName[MaxCharactersInLine];
    strcpy(FileName,Pall.RestartFileName);

    fprintf(stderr,"%s\n",Pall.RestartFileName);
    sprintf(fname,"%s.%03d.%03d",Pall.RestartFileName,TotalFileNumber,0);

    fp = fopen(fname,"rb");
    if(fp == NULL){
        Pall.RunStatus = NONE;
        return ;
    } else {
        fclose(fp);
    }

    FileOpen(fp,fname,"rb");
    ReadHeader(fp,HeaderArray);
    fread(&Pall,sizeof(struct StructPall),1,fp);
    fread(&TimingResults,sizeof(struct StructTimingResults),1,fp);
    fclose(fp);

    fprintf(stderr,"%s\n",Pall.RestartFileName);
    fprintf(stderr,"%s\n",Pall.BaseFileName);

    GenerateStructPbody(Pall.Ntotal_t);
    GenerateStructPhydro(Pall.Nhydro_t);
    GenerateStructPstar(Pall.Nstars_t);
    GenerateStructPsink(Pall.Nsink_t);

    fprintf(stderr," [%02d] Reload RestartFileData : (Ntotal,Ndm,Nhydro,Nstar,Nsink) = (%ld,%ld,%ld,%ld,%ld)\n",MyID,
            Pall.Ntotal_t,Pall.NDM_t,Pall.Nhydro_t,Pall.Nstars_t,Pall.Nsink_t);


    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(HeaderArray[HD_CompactFormatFlag] == ON){
            fprintf(stderr,"IO format : Compact mode ON\n");
            if(HeaderArray[HD_LeanFormatFlag] == ON){
                fprintf(stderr,"IO format : Lean mode ON\n");
            }
        }
    }

    AllocateRandomGenerator();

    int CountTypes[NTypes];
    for(int i=0;i<NTypes;i++)
        CountTypes[i] = 0;

    int counter = 0;
    for(int i=0;i<TotalFileNumber;i++){

        sprintf(fname,"%s.%03d.%03d",FileName,TotalFileNumber,i);
        fprintf(stderr,"%s\n",fname);
        FileOpen(fp,fname,"rb");

        struct  StructPall     TempPall;
        struct  StructTimingResults     TempTimingResults;
        int TempHeaderArray[HEADER_SIZE];

        ReadHeader(fp,TempHeaderArray);
        fread(&TempPall,sizeof(struct StructPall),1,fp);
        fread(&TempTimingResults,sizeof(struct StructTimingResults),1,fp);

        fprintf(stderr," [%02d] Reload RestartFileData : (Ntotal,Ndm,Nhydro,Nstar,Nsink) = (%ld,%ld,%ld,%ld,%ld)\n",MyID,
                TempPall.Ntotal,TempPall.NDM,TempPall.Nhydro,TempPall.Nstars,TempPall.Nsink);

        if((HeaderArray[HD_CompactFormatFlag] == false)
                &&(HeaderArray[HD_CompactDoubleFormatFlag] == false)
                &&(HeaderArray[HD_CompactMixFormatFlag] == false)){
            counter = ReadParallelDataOnSingleNodeFull(fp,TempPall.Ntotal,counter,CountTypes);
        } else if((HeaderArray[HD_CompactFormatFlag] == true)&&(HeaderArray[HD_LeanFormatFlag] == false)){
            counter = ReadParallelDataOnSingleNodeCompact(fp,TempPall.Ntotal,counter,CountTypes);
        } else if((HeaderArray[HD_CompactFormatFlag] == true)&&(HeaderArray[HD_LeanFormatFlag] == true)){
            counter = ReadParallelDataOnSingleNodeLean(fp,TempPall.Ntotal,counter,CountTypes);
        } else if(HeaderArray[HD_CompactDoubleFormatFlag] == true){
            counter = ReadParallelDataOnSingleNodeCompactDouble(fp,TempPall.Ntotal,counter,CountTypes);
        } else if(HeaderArray[HD_CompactMixFormatFlag] == true){
            counter = ReadParallelDataOnSingleNodeCompactMix(fp,TempPall.Ntotal,counter,CountTypes);
        }

        if(i == 0){
            if(GSL_EFAILED == gsl_rng_fread(fp,RandomGenerator)){
                fprintf(stderr,"The random generator read error!\n");
                MPI_Abort(MPI_COMM_WORLD,GSL_RANDOM_GENERATOR_READ_ERROR);
            }
        }
        fclose(fp);
    }
    Pall.Ntotal = Pall.Ntotal_t;
    Pall.NDM = Pall.NDM_t;
    Pall.Nhydro = Pall.Nhydro_t;
    Pall.Nstars = Pall.Nstars_t;
    Pall.Nsink = Pall.Nsink_t;


    fprintf(stderr," [%02d] Reload RestartFileData : (Ntotal,Ndm,Nhydro,Nstar,Nsink) = (%d,%d,%d,%d,%d)\n",MyID,
            CountTypes[TypeDM]+CountTypes[TypeHydro]+CountTypes[TypeStar]+CountTypes[TypeSink],
            CountTypes[TypeDM],CountTypes[TypeHydro],CountTypes[TypeStar],CountTypes[TypeSink]);

    ReConnectPointers();

    return;
}

/*
 * Routines for ASCII format files. 
 */

void OutPutASCIIDATA(void){

	FILE *fp_hydro,*fp_star,*fp_dm;

    char fname_hydro[MaxCharactersInLine],fname_star[MaxCharactersInLine],fname_dm[MaxCharactersInLine];
    sprintf(fname_hydro,"%s.Hydro.%02d.%02d",Pall.ASCIIFileName,MPIGetNumProcs(),MPIGetMyID());
    sprintf(fname_star,"%s.Star.%02d.%02d",Pall.ASCIIFileName,MPIGetNumProcs(),MPIGetMyID());
    sprintf(fname_dm,"%s.DM.%02d.%02d",Pall.ASCIIFileName,MPIGetNumProcs(),MPIGetMyID());

    //if(Pall.Nhydro > 0)
        FileOpen(fp_hydro,fname_hydro,"w");
    //if(Pall.Nstars > 0)
        FileOpen(fp_star,fname_star,"w");
    //if(Pall.NDM > 0)
        FileOpen(fp_dm,fname_dm,"w");

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeHydro){
            fprintf(fp_hydro,"%ld %g %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                    Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                    Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                    Pbody[i]->Mass,Pall.ConvertNumberDensityToCGS*PbodyHydro(i)->Rho,
                    PbodyHydroU(i)*Pall.ConvertUtoT,PbodyHydroKernel(i));
        } else if(Pbody[i]->Type == TypeStar){
            fprintf(fp_star,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                    Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                    Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                    Pbody[i]->Mass);
        } else {
            fprintf(fp_dm,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                    Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                    Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                    Pbody[i]->Mass);
        }
    }

    //if(Pall.Nhydro > 0)
        fclose(fp_hydro);
    //if(Pall.Nstars > 0)
        fclose(fp_star);
    //if(Pall.NDM > 0)
        fclose(fp_dm);

	return;
}




//#define WRITE_DISPH_COMPATIBLE_MODE

#ifdef WRITE_DISPH_COMPATIBLE_MODE
void WriteDISPHAnalysisToolsCompatibleMode(void){

    for(int i=0;i<MPIGetNumProcs();i++){
        if(MPIGetMyID() == i){
            FILE *fp;
            char fname[MaxCharactersInLine];
            sprintf(fname,"%s.%04d.%03d.data",Pall.BaseFileName,
                    Pall.OutPutFileNumber,MPIGetMyID());
            FileOpen(fp,fname,"w"); 

            for(int i=0;i<Pall.Nhydro;i++){
#ifdef USE_DISPH
                fprintf(fp,"%g %g %g %g %g %g %g %g %g %g %d %d\n",
                Phydro[i]->PosP[0], Phydro[i]->PosP[1],
                PhydroBody(i)->Vel[0], PhydroBody(i)->Vel[1],
                Phydro[i]->Mass,Phydro[i]->Kernel,Phydro[i]->Rho,Phydro[i]->U,
                Pall.Gm1*Phydro[i]->Rho*Phydro[i]->U,Phydro[i]->EnergyDensity,
                Phydro[i]->Nlist,
#ifdef USE_PARTICLE_TAG
                Phydro[i]->Tag);
#else // USE_PARTICLE_TAG
                0);
#endif // USE_PARTICLE_TAG
#else //USE_DISPH
                fprintf(fp,"%g %g %g %g %g %g %g %g %g %g %g %g %d %d\n",
                Phydro[i]->PosP[0], Phydro[i]->PosP[1],
                PhydroBody(i)->Vel[0], PhydroBody(i)->Vel[1],
                Phydro[i]->Mass,Phydro[i]->Kernel,Phydro[i]->Rho,Phydro[i]->U,
                Pall.Gm1*Phydro[i]->Rho*Phydro[i]->U,Phydro[i]->Rho*Phydro[i]->U,
                Pall.Gm1*Phydro[i]->Rho*Phydro[i]->U,Phydro[i]->Rho*Phydro[i]->U,
                Phydro[i]->Nlist,
#ifdef USE_PARTICLE_TAG
                Phydro[i]->Tag);
#else // USE_PARTICLE_TAG
                0);
#endif // USE_PARTICLE_TAG
#endif //USE_DISPH
            }
            fclose(fp);
            fflush(NULL);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if(MPIGetMyID() == MPI_ROOT_RANK){
        // Generate Header.
        FILE *fp;
        char fname[MaxCharactersInLine];
        sprintf(fname,"%s.head.%04d.data",Pall.BaseFileName,Pall.OutPutFileNumber);
#if 0
        FileOpen(fp,fname,"w"); 
#ifdef USE_DISPH
        fprintf(fp,"1\n");
#else // USE_DISPH
        fprintf(fp,"0\n");
#endif // USE_DISPH
        fprintf(fp,"%ld\n",Pall.Nhydro_t);
        fprintf(fp,"%g\n",Pall.TCurrent);
        fprintf(fp,"%g\n",Pall.Gamma);
        fflush(NULL);
        fclose(fp);
#endif

        for(int i=0;i<MPIGetNumProcs();i++){
            char run_command[MaxCharactersInLine];
            sprintf(run_command,"cat %s.head.%04d.data %s.%04d.%03d.data > %s.%04d",
                    Pall.BaseFileName,Pall.OutPutFileNumber,
                    Pall.BaseFileName,Pall.OutPutFileNumber,i,
                    Pall.BaseFileName,Pall.OutPutFileNumber);
            //fprintf(stderr,"************* %s\n",run_command);
            system(run_command);
            fflush(NULL);
            sprintf(run_command,"rm -rf %s.head.%04d.data",
                    Pall.BaseFileName,Pall.OutPutFileNumber);
            system(run_command);
            sprintf(run_command,"rm -rf %s.%04d.%03d.data",
                    Pall.BaseFileName,Pall.OutPutFileNumber,i);
            system(run_command);
            fflush(NULL);
            //exit(1);
        }
    }

    return ;
}

#endif // WRITE_DISPH_COMPATIBLE_MODE


/*
 * This function writes important data of hydro particles. There are two modes.
 * If mode == NONE, this function writes data of all hydro particles. On the
 * other hand, if mode == i(>=0), this function writes data of i-th particle's
 * data. Even in the parallel runs, an output file, which merges all data of all
 * nodes, is stored in the current directory.
 */
void WriteHydroDataASCIIFormat(const int mode){

    FILE *fp;
    char fname[MaxCharactersInLine];

    Snprintf(fname,"./HydroData.%04d.dat",MPIGetMyID());
    FileOpen(fp,fname,"w");

    if(mode == NONE){ // Write all.
        for(int i=0;i<Pall.Nhydro;i++){
#ifdef USE_GRAD_H //{
            fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g %g %g\n",
                    PhydroBody(i)->GlobalID,
                    PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                    PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2],
                    Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2],
                    Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->Nlist,Phydro[i]->Gradh,
                    Phydro[i]->U,Pall.ConvertUtoT*Phydro[i]->U,
                    Phydro[i]->Du,Phydro[i]->EnergyDensity);
#else
            fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g %g\n",
                    PhydroBody(i)->GlobalID,
                    PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                    PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2],
                    Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2],
                    Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->Nlist,
                    Phydro[i]->U,Pall.ConvertUtoT*Phydro[i]->U,
                    Phydro[i]->Du,Phydro[i]->EnergyDensity);
#endif // USE_GRAD_H //}
        }
    } else {
        for(int i=0;i<Pall.Nhydro;i++){
            if(PhydroBody(i)->GlobalID == mode){
#ifdef USE_GRAD_H //{
                fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g %g %g\n",
                        PhydroBody(i)->GlobalID,
                        PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                        PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2],
                        Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2],
                        Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->Nlist,Phydro[i]->Gradh,
                        Phydro[i]->U,Pall.ConvertUtoT*Phydro[i]->U,
                        Phydro[i]->Du,Phydro[i]->EnergyDensity);
#else
                fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g %g\n",
                        PhydroBody(i)->GlobalID,
                        PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                        PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2],
                        Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2],
                        Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->Nlist,
                        Phydro[i]->U,Pall.ConvertUtoT*Phydro[i]->U,
                        Phydro[i]->Du,Phydro[i]->EnergyDensity);
#endif // USE_GRAD_H //}
            }
        }
    }
    fclose(fp);
    fflush(NULL);

    // merger!
    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(mode == NONE){
            system("cat ./HydroData.????.dat | sort -n > ./HydroData.dat");
        } else {
            sprintf(fname,"cat ./HydroData.????.dat | sort -n > ./HydroData.ID%06d.dat",mode);
            system(fname);
        }
        fflush(NULL);
        system("rm -rf ./HydroData.????.dat");
        fflush(NULL);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    return ;
}

void ShowHydroDataASCIIFormat(const int mode){

    for(int k=0;k<MPIGetNumProcs();k++){
        if(k == MPIGetMyID()){
            if(mode == NONE){ // Write all.
                for(int i=0;i<Pall.Nhydro;i++){
#ifdef USE_GRAD_H //{
                    fprintf(stderr,"%ld %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g %g %g\n",
                            PhydroBody(i)->GlobalID,
                            PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                            PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2],
                            Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2],
                            Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->Nlist,Phydro[i]->Gradh,
                            Phydro[i]->U,Pall.ConvertUtoT*Phydro[i]->U,
                            Phydro[i]->Du,Phydro[i]->EnergyDensity);
#else 
                    fprintf(stderr,"%ld %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g %g\n",
                            PhydroBody(i)->GlobalID,
                            PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                            PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2],
                            Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2],
                            Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->Nlist,
                            Phydro[i]->U,Pall.ConvertUtoT*Phydro[i]->U,
                            Phydro[i]->Du,Phydro[i]->EnergyDensity);
#endif // USE_GRAD_H //}
                }
            } else {
                for(int i=0;i<Pall.Nhydro;i++){
                    if(PhydroBody(i)->GlobalID == mode){
                        fprintf(stderr,"GID P[] V[] A[] R h Nb Gh U T Q Du q\n");
#ifdef USE_GRAD_H //{
                        /*
                        fprintf(stderr,"%ld %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g %g %g\n",
                                PhydroBody(i)->GlobalID,
                                PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                                PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2],
                                Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2],
                                Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->Nlist,Phydro[i]->Gradh,
                                Phydro[i]->U,Pall.ConvertUtoT*Phydro[i]->U,
                                Phydro[i]->Du,Phydro[i]->EnergyDensity);
                        */
                        fprintf(stderr,"%ld | %g %g %g | %g %g %g | %g %g %g |%g %g %d %g %g %g %g %g %g\n",
                                PhydroBody(i)->GlobalID,
                                PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                                PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2],
                                Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2],
                                Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->Nlist,Phydro[i]->Gradh,
                                Phydro[i]->U,Pall.ConvertUtoT*Phydro[i]->U,Phydro[i]->DQheat,
                                Phydro[i]->Du,Phydro[i]->EnergyDensity);
#else
                        fprintf(stderr,"GID P[] V[] A[] Nb U Q Du q\n");
                        fprintf(stderr,"%ld %g %g %g %g %g %g %d %g %g %g %g %g\n",
                                PhydroBody(i)->GlobalID,
                                PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                                Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2],
                                Phydro[i]->Nlist,
                                Phydro[i]->U,Pall.ConvertUtoT*Phydro[i]->U,Phydro[i]->DQheat,
                                Phydro[i]->Du,Phydro[i]->EnergyDensity);
#endif // USE_GRAD_H //}
                    }
                }
            }
        }
        fflush(NULL);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    return ;
}

/*
 * This function writes all of neighbors' information in a file.
 */
void WriteNeighborInformation(const int index, const int mode){

    int Nlist;
    int Neighbors[MaxNeighborSize];
    if(mode == 0){
        Nlist = GetNeighborsLimited(PhydroBody(index)->Pos,2.e0*Phydro[index]->Kernel,Neighbors);
    } else if (mode == 1){
        Nlist = GetNeighborsPairsLimited(PhydroBody(index)->Pos,2.e0*Phydro[index]->Kernel,Neighbors);
    } else if (mode == 2){
        Nlist = GetNeighborsPairsLimitedImported(PhydroBody(index)->Pos,2.e0*Phydro[index]->Kernel,Neighbors);
    } else {
        assert(mode < 3);
    }

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"./NeighborLog.%06ld.%d",PhydroBody(index)->GlobalID,mode);
    FileOpen(fp,fname,"w");
    fprintf(fp,"#Neighbor information of particle with GID = %06ld\n",PhydroBody(index)->GlobalID);
#ifdef USE_GRAD_H //{
    fprintf(fp,"#GID#1, Pos[3]#2-4, Vel[3]#5-7, HAcc[3]#8-10, Rho#11, K#12, Nlist#13, Gradh#14, U#15, T#16, Du#17, DuCooling#18, q#19\n");
#else
    fprintf(fp,"#GID#1, Pos[3]#2-4, Vel[3]#5-7, HAcc[3]#8-10, Rho#11, K#12, Nlist#13, U#15, T#16, Du#17, DuCooling#18, q#19\n");
#endif // USE_GRAD_H //}
    for(int i=0;i<Nlist;i++){
        int leaf = Neighbors[i];
#ifdef USE_GRAD_H //{
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g %g %g %g %g\n",
                PhydroBody(leaf)->GlobalID,
                PhydroBody(leaf)->Pos[0],PhydroBody(leaf)->Pos[1],PhydroBody(leaf)->Pos[2],
                PhydroBody(leaf)->Vel[0],PhydroBody(leaf)->Vel[1],PhydroBody(leaf)->Vel[2],
                Phydro[leaf]->HydroAcc[0],Phydro[leaf]->HydroAcc[1],Phydro[leaf]->HydroAcc[2],
                Phydro[leaf]->Rho,Phydro[leaf]->Kernel,Phydro[leaf]->Nlist,Phydro[leaf]->Gradh,
                Phydro[leaf]->U,Pall.ConvertUtoT*Phydro[leaf]->U,Phydro[leaf]->DQheat,
                Phydro[leaf]->Du,Phydro[leaf]->DuCooling,Phydro[leaf]->EnergyDensity);
#else
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g %g %g %g\n",
                PhydroBody(leaf)->GlobalID,
                PhydroBody(leaf)->Pos[0],PhydroBody(leaf)->Pos[1],PhydroBody(leaf)->Pos[2],
                PhydroBody(leaf)->Vel[0],PhydroBody(leaf)->Vel[1],PhydroBody(leaf)->Vel[2],
                Phydro[leaf]->HydroAcc[0],Phydro[leaf]->HydroAcc[1],Phydro[leaf]->HydroAcc[2],
                Phydro[leaf]->Rho,Phydro[leaf]->Kernel,Phydro[leaf]->Nlist,
                Phydro[leaf]->U,Pall.ConvertUtoT*Phydro[leaf]->U,Phydro[leaf]->DQheat,
                Phydro[leaf]->Du,Phydro[leaf]->DuCooling,Phydro[leaf]->EnergyDensity);
#endif // USE_GRAD_H //}
    }
    fclose(fp);

    return ;
}

/*
 * This function shows all of neighbors' information in the stderr.
 */
void ShowNeighborInformation(const int index, const int mode){

    int Nlist;
    int Neighbors[MaxNeighborSize];
    if(mode == 0){
        Nlist = GetNeighborsLimited(PhydroBody(index)->Pos,2.e0*Phydro[index]->Kernel,Neighbors);
    } else if (mode == 1){
        Nlist = GetNeighborsPairsLimited(PhydroBody(index)->Pos,2.e0*Phydro[index]->Kernel,Neighbors);
    } else if (mode == 2){
        Nlist = GetNeighborsPairsLimitedImported(PhydroBody(index)->Pos,2.e0*Phydro[index]->Kernel,Neighbors);
    } else {
        assert(mode < 3);
    }

    fprintf(stderr,"#Neighbor information of particle with GID = %06ld\n",PhydroBody(index)->GlobalID);
#ifdef USE_GRAD_H //{
    fprintf(stderr,"#GID#1, Pos[3]#2-4, Vel[3]#5-7, HAcc[3]#8-10, Rho#11, K#12, Nlist#13, Gradh#14, U#15, T#16, Q#17 Du#18, DuCooling#19, q#20\n");
#else
    fprintf(stderr,"#GID#1, Pos[3]#2-4, Vel[3]#5-7, HAcc[3]#8-10, Rho#11, K#12, Nlist#13, U#14, T#15, Q#16 Du#17, DuCooling#18, q#19\n");
#endif // USE_GRAD_H //}
    for(int i=0;i<Nlist;i++){
        int leaf = Neighbors[i];
#ifdef USE_GRAD_H //{
        /*
        fprintf(stderr,"%ld %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g %g %g %g %g\n",
                PhydroBody(leaf)->GlobalID,
                PhydroBody(leaf)->Pos[0],PhydroBody(leaf)->Pos[1],PhydroBody(leaf)->Pos[2],
                PhydroBody(leaf)->Vel[0],PhydroBody(leaf)->Vel[1],PhydroBody(leaf)->Vel[2],
                Phydro[leaf]->HydroAcc[0],Phydro[leaf]->HydroAcc[1],Phydro[leaf]->HydroAcc[2],
                Phydro[leaf]->Rho,Phydro[leaf]->Kernel,Phydro[leaf]->Nlist,Phydro[leaf]->Gradh,
                Phydro[leaf]->U,Pall.ConvertUtoT*Phydro[leaf]->U,Phydro[leaf]->DQheat,
                Phydro[leaf]->Du,Phydro[leaf]->DuCooling,Phydro[leaf]->EnergyDensity);
                */
        /*
        fprintf(stderr,"%ld %g %g %d %g %g %g %g %g %g %g\n",
                PhydroBody(leaf)->GlobalID,
                Phydro[leaf]->Rho,Phydro[leaf]->Kernel,Phydro[leaf]->Nlist,Phydro[leaf]->Gradh,
                Phydro[leaf]->U,Pall.ConvertUtoT*Phydro[leaf]->U,Phydro[leaf]->DQheat,
                Phydro[leaf]->Du,Phydro[leaf]->DuCooling,Phydro[leaf]->EnergyDensity);
                */
        fprintf(stderr,"%ld %g %g %g %g %g %d %g %g %g %g %g %g %g\n",
                PhydroBody(leaf)->GlobalID,
                Phydro[leaf]->HydroAcc[0],Phydro[leaf]->HydroAcc[1],Phydro[leaf]->HydroAcc[2],
                Phydro[leaf]->Rho,Phydro[leaf]->Kernel,Phydro[leaf]->Nlist,Phydro[leaf]->Gradh,
                Phydro[leaf]->U,Pall.ConvertUtoT*Phydro[leaf]->U,Phydro[leaf]->DQheat,
                Phydro[leaf]->Du,Phydro[leaf]->DuCooling,Phydro[leaf]->EnergyDensity);
#else
        fprintf(stderr,"%ld %g %g %g %g %g %g %g %g %g %g %g %d %g %g %g %g %g %g\n",
                PhydroBody(leaf)->GlobalID,
                PhydroBody(leaf)->Pos[0],PhydroBody(leaf)->Pos[1],PhydroBody(leaf)->Pos[2],
                PhydroBody(leaf)->Vel[0],PhydroBody(leaf)->Vel[1],PhydroBody(leaf)->Vel[2],
                Phydro[leaf]->HydroAcc[0],Phydro[leaf]->HydroAcc[1],Phydro[leaf]->HydroAcc[2],
                Phydro[leaf]->Rho,Phydro[leaf]->Kernel,Phydro[leaf]->Nlist,
                Phydro[leaf]->U,Pall.ConvertUtoT*Phydro[leaf]->U,Phydro[leaf]->DQheat,
                Phydro[leaf]->Du,Phydro[leaf]->DuCooling,Phydro[leaf]->EnergyDensity);
#endif // USE_GRAD_H //}
    }

    return ;
}

/* 
 * This function stores the time-steps of particles and related information on files.
 */
void WriteTimeSteps(void){

    // Write All
    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"./dt.%05d",MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        double ParticleEpsSize = Pbody[i]->Eps;
        double Acc[3] = {Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2]};
        double normAcc = NORM(Acc);
        double DiffAcc[3] = {Acc[0]-Pbody[i]->AccOld[0],
                             Acc[1]-Pbody[i]->AccOld[1],
                             Acc[2]-Pbody[i]->AccOld[2]};
        double normDiffAcc = NORM(DiffAcc);
        double normVel = NORM(Pbody[i]->Vel);

        fprintf(fp,"%ld %g %g %g %g\n",Pbody[i]->GlobalID,Pbody[i]->dt,
                sqrt(ParticleEpsSize/NORM(Pbody[i]->Acc)),
                //normAcc/normDiffAcc*Pbody[i]->dt,
                NORM(Pbody[i]->Pos),NORM(Pbody[i]->Vel));
    }
    fclose(fp);
    fflush(NULL);

#if 0
    // merger!
    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat ./dt.????? | sort -n > ./dt.dat");
        fflush(NULL);
        system("rm -rf ./dt.?????");
        fflush(NULL);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // Write Hydro
    Snprintf(fname,"./dt.hydro.%05d",MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        double ParticleEpsSize = PhydroBody(i)->Eps;
        double Acc[3] = {PhydroBody(i)->Acc[0],
                         PhydroBody(i)->Acc[1],
                         PhydroBody(i)->Acc[2]};
        double normAcc = NORM(Acc);
        double DiffAcc[3] = {Acc[0]-PhydroBody(i)->AccOld[0],
                             Acc[1]-PhydroBody(i)->AccOld[1],
                             Acc[2]-PhydroBody(i)->AccOld[2]};
        double normDiffAcc = NORM(DiffAcc);
        double normVel = NORM(PhydroBody(i)->Vel);

        fprintf(fp,"%ld %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroBody(i)->dt,
                Phydro[i]->dt_hydro,
                sqrt(ParticleEpsSize/NORM(PhydroBody(i)->Acc)),
                //normAcc/normDiffAcc*PhydroBody(i)->dt,
                NORM(PhydroBody(i)->Pos),NORM(PhydroBody(i)->Vel));
    }
    fclose(fp);

    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

#if 0
    // merger!
    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat ./dt.hydro.????? | sort -n > ./dt.hydro.dat");
        fflush(NULL);
        system("rm -rf ./dt.hydro.?????");
        fflush(NULL);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

#endif

    return ;
}


void WriteCurrentActiveParticles(const long int FileID, char suffix[]){

    // Write All
    if(MPIGetMyID() == MPI_ROOT_RANK)
        MakeDir("./Log");
    MPI_Barrier(MPI_COMM_WORLD);

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"./Log/Actives.%s.%ld.%05d",suffix,FileID,MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active){
            fprintf(fp,"%lu %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d\n",Pbody[i]->GlobalID,Pbody[i]->dt,
                    PhydroBody(i)->PosP[0],PhydroBody(i)->PosP[1],PhydroBody(i)->PosP[2],
                    //PhydroBody(i)->VelP[0],PhydroBody(i)->VelP[1],PhydroBody(i)->VelP[2],
                    Phydro[i]->VelP[0],Phydro[i]->VelP[1],Phydro[i]->VelP[2],
                    Phydro[i]->HydroAcc[0],Phydro[i]->HydroAcc[1],Phydro[i]->HydroAcc[2],
                    Pall.ConvertNumberDensityToCGS*Phydro[i]->RhoPred,Phydro[i]->KernelPred,Pall.ConvertUtoT*Phydro[i]->UPred,
                    Phydro[i]->Gradh,
                    Phydro[i]->EnergyDensity, Phydro[i]->EnergyDensityPred,
#ifdef USE_SPSPH //{
                    Phydro[i]->PseudoDensityPred, 
                    Phydro[i]->ZwPred, 
#else // USE_SPSPH //} //{
                    0.0,0.0,
#endif // USE_SPSPH //}
                    Phydro[i]->Nlist);
        }

    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    // merger!
    if(MPIGetMyID() == MPI_ROOT_RANK){
        char cmd[MaxCharactersInLine];
        Snprintf(cmd,"cat ./Log/Actives.%s.%ld.????? | sort -n > ./Log/Actives.%s.%ld",
                suffix,FileID,suffix,FileID);
        system(cmd);

        fflush(NULL);
        Snprintf(cmd,"rm -rf ./Log/Actives.%s.%ld.?????",suffix,FileID);
        system(cmd);
        fflush(NULL);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    return ;
}



#ifdef USE_ON_THE_FLY_4D2U //{

static int Write4D2UExtraCounter = 1;
/*
 * This is the counter which points the timing to write the extra 4D2U IO.
 * This covers from 1 to ON_THE_FLY_4D2U_EXTRA_INFO_INTERVAL.
 * When Write4D2U() is called, the counter resets to 1.
 */
static inline __attribute__((always_inline)) void *SendBuf(const bool target, void *sendbuf){
    return target
        ? MPI_IN_PLACE
        : sendbuf;
}


static void ASR_Write4D2UFormatHeader(FILE *fp, const int Index){

    fprintf(fp,"%d\n",Index);
    fprintf(fp,"%g\n",Pall.TCurrent*Pall.UnitTime/YEAR_CGS);
    fprintf(fp,"Time unit : yr.\n");
    fprintf(fp,"Length unit : pc.\n");
    fprintf(fp,"Mass unit : M_sun.\n");
    fprintf(fp,"Velocity unit : km/s.\n");
    fprintf(fp,"Nhydro %ld\n",Pall.Nhydro_t);
    fprintf(fp,"Nstar %ld\n",Pall.Nstars_t);
    fprintf(fp,"Nsink %ld\n",Pall.Nsink_t);
    fprintf(fp,"NDM %ld\n",Pall.NDM_t);

    return ;
}

static void ASR_Write4D2UFormatHeaderExtra(FILE *fp, const int Index, const int Numbers[]){

    fprintf(fp,"%d\n",Index);
    fprintf(fp,"%g\n",Pall.TCurrent*Pall.UnitTime/YEAR_CGS);
    fprintf(fp,"Time unit : yr.\n");
    fprintf(fp,"Length unit : pc.\n");
    fprintf(fp,"Mass unit : M_sun.\n");
    fprintf(fp,"Velocity unit : km/s.\n");
    fprintf(fp,"Nhydro %d\n",Numbers[0]);
    fprintf(fp,"Nstar %d\n",Numbers[1]);
    fprintf(fp,"Nsink %d\n",Numbers[2]);
    fprintf(fp,"NDM %d\n",Numbers[3]);

    return ;
}

struct StructWrite4D2UBody{
    int Type;
    long int GlobalID;
    float Mass;
    float Pos[3];
    float Vel[3];

    int   Tag;   // Particle Tag
    float Data1; // 
    float Data2; // 
    float Data3; // 
    float Data4; // 
};


void InitWrite4D2U(void){

    char dname[MaxCharactersInLine] = "4D2U";
    MakeDir(dname);
    MPI_Barrier(MPI_COMM_WORLD);

    return ;
}


static void ASR_Pack4D2UFormatHydro(struct StructWrite4D2UBody Body[]){

    const double LengthUnitConverter = Pall.UnitLength/PC_CGS;
    const double VelocityUnitConverter = (Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;

    for(int i=0;i<Pall.Nhydro;i++){
        Body[i].Type = TypeHydro;
        Body[i].GlobalID = PhydroBody(i)->GlobalID;
#ifdef USE_PARTICLE_TAG
        Body[i].Tag = Phydro[i]->Tag;
#else
        Body[i].Tag = 0;
#endif

        Body[i].Pos[0] = LengthUnitConverter*PhydroBody(i)->Pos[0];
        Body[i].Pos[1] = LengthUnitConverter*PhydroBody(i)->Pos[1];
        Body[i].Pos[2] = LengthUnitConverter*PhydroBody(i)->Pos[2];
        Body[i].Vel[0] = VelocityUnitConverter*PhydroBody(i)->Vel[0];
        Body[i].Vel[1] = VelocityUnitConverter*PhydroBody(i)->Vel[1];
        Body[i].Vel[2] = VelocityUnitConverter*PhydroBody(i)->Vel[2];

        Body[i].Mass  = PhydroBody(i)->Mass*Pall.UnitMass/MSUN_CGS;
        Body[i].Data1 = LengthUnitConverter*Phydro[i]->Kernel;
        Body[i].Data2 = Pall.ConvertNumberDensityToCGS*Phydro[i]->Rho;
        Body[i].Data3 = Pall.ConvertUtoT*Phydro[i]->U;
        Body[i].Data4 = Phydro[i]->Z;
    }

    return ;
}

static void ASR_Pack4D2UFormatStar(struct StructWrite4D2UBody Body[]){

    const double LengthUnitConverter = Pall.UnitLength/PC_CGS;
    const double VelocityUnitConverter = (Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;

    for(int i=0;i<Pall.Nstars;i++){
        Body[i].Type = TypeStar;
        Body[i].GlobalID = PstarBody(i)->GlobalID;
        Body[i].Tag = 0;

        Body[i].Pos[0] = LengthUnitConverter*PstarBody(i)->Pos[0];
        Body[i].Pos[1] = LengthUnitConverter*PstarBody(i)->Pos[1];
        Body[i].Pos[2] = LengthUnitConverter*PstarBody(i)->Pos[2];
        Body[i].Vel[0] = VelocityUnitConverter*PstarBody(i)->Vel[0];
        Body[i].Vel[1] = VelocityUnitConverter*PstarBody(i)->Vel[1];
        Body[i].Vel[2] = VelocityUnitConverter*PstarBody(i)->Vel[2];

        Body[i].Mass  = PstarBody(i)->Mass*Pall.UnitMass/MSUN_CGS;
        Body[i].Data1 = LengthUnitConverter*PstarBody(i)->Eps;
        Body[i].Data2 = (Pall.TCurrent-Pstar[i]->FormationTime)*Pall.UnitTime/YEAR_CGS;//Age
        Body[i].Data3 = 0.e0;
        Body[i].Data4 = Pstar[i]->Z;// Metallicity
    }

    return ;
}

static void ASR_Pack4D2UFormatSink(struct StructWrite4D2UBody Body[]){

    const double LengthUnitConverter = Pall.UnitLength/PC_CGS;
    const double VelocityUnitConverter = (Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;

    for(int i=0;i<Pall.Nsink;i++){
        Body[i].Type = TypeSink;
        Body[i].GlobalID = PsinkBody(i)->GlobalID;
        Body[i].Tag = 0;

        Body[i].Pos[0] = LengthUnitConverter*PsinkBody(i)->Pos[0];
        Body[i].Pos[1] = LengthUnitConverter*PsinkBody(i)->Pos[1];
        Body[i].Pos[2] = LengthUnitConverter*PsinkBody(i)->Pos[2];
        Body[i].Vel[0] = VelocityUnitConverter*PsinkBody(i)->Vel[0];
        Body[i].Vel[1] = VelocityUnitConverter*PsinkBody(i)->Vel[1];
        Body[i].Vel[2] = VelocityUnitConverter*PsinkBody(i)->Vel[2];

        Body[i].Mass  = PsinkBody(i)->Mass*Pall.UnitMass/MSUN_CGS;
        Body[i].Data1 = LengthUnitConverter*PsinkBody(i)->Eps;
        Body[i].Data2 = Psink[i]->AccretionRadius*Pall.UnitLength/PC_CGS;
        Body[i].Data3 = 0.e0;
        Body[i].Data4 = Psink[i]->Z;
    }

    return ;
}

static void ASR_Pack4D2UFormatDM(struct StructWrite4D2UBody Body[]){

    const double LengthUnitConverter = Pall.UnitLength/PC_CGS;
    const double VelocityUnitConverter = (Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;

    int *Index;
    Index = malloc(sizeof(int)*Pall.NDM);

    int counter = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeDM){
            Index[counter] = i;
            counter ++;
        }
    }

    for(int k=0;k<Pall.NDM;k++){
        int i = Index[k];
        Body[k].Type = TypeDM;
        Body[k].GlobalID = Pbody[i]->GlobalID;
        Body[k].Tag = 0;

        Body[k].Pos[0] = LengthUnitConverter*Pbody[i]->Pos[0];
        Body[k].Pos[1] = LengthUnitConverter*Pbody[i]->Pos[1];
        Body[k].Pos[2] = LengthUnitConverter*Pbody[i]->Pos[2];
        Body[k].Vel[0] = VelocityUnitConverter*Pbody[i]->Vel[0];
        Body[k].Vel[1] = VelocityUnitConverter*Pbody[i]->Vel[1];
        Body[k].Vel[2] = VelocityUnitConverter*Pbody[i]->Vel[2];

        Body[k].Mass =  Pbody[i]->Mass*Pall.UnitMass/MSUN_CGS;
        Body[k].Data1 = LengthUnitConverter*Pbody[i]->Eps;
        Body[k].Data2 = 0.e0;
        Body[k].Data3 = 0.e0;
        Body[k].Data4 = 0.e0;
    }

    free(Index);

    return ;
}


static void ASR_Pack4D2UFormat(struct StructWrite4D2UBody Body[], const int Type){

    if(Type == 0){ // Hydro
        ASR_Pack4D2UFormatHydro(Body);
    }else if(Type == 1){ // Star
        ASR_Pack4D2UFormatStar(Body);
    }else if(Type == 2){ // Sink
        ASR_Pack4D2UFormatSink(Body);
    }else if(Type == 3){ // DM
        ASR_Pack4D2UFormatDM(Body);
    }

    return ;
}

static void ASR_Pack4D2UFormatExtraHydro(struct StructWrite4D2UBody Body[]){

    const double LengthUnitConverter = Pall.UnitLength/PC_CGS;
    const double VelocityUnitConverter = (Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;
    const double ExtraInterval = Pall.OutPutInterval/ON_THE_FLY_4D2U_EXTRA_INFO_INTERVAL;

    int counter = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->dt_hydro*ON_THE_FLY_4D2U_EXTRA_INFO_CRITERIA >= ExtraInterval) continue;

        Body[counter].Type = TypeHydro;
        Body[counter].GlobalID = PhydroBody(i)->GlobalID;
#ifdef USE_PARTICLE_TAG
        Body[counter].Tag = Phydro[i]->Tag;
#else
        Body[counter].Tag = 0;
#endif

        Body[counter].Pos[0] = LengthUnitConverter*PhydroBody(i)->Pos[0];
        Body[counter].Pos[1] = LengthUnitConverter*PhydroBody(i)->Pos[1];
        Body[counter].Pos[2] = LengthUnitConverter*PhydroBody(i)->Pos[2];
        Body[counter].Vel[0] = VelocityUnitConverter*PhydroBody(i)->Vel[0];
        Body[counter].Vel[1] = VelocityUnitConverter*PhydroBody(i)->Vel[1];
        Body[counter].Vel[2] = VelocityUnitConverter*PhydroBody(i)->Vel[2];

        Body[counter].Mass  = PhydroBody(i)->Mass*Pall.UnitMass/MSUN_CGS;
        Body[counter].Data1 = LengthUnitConverter*Phydro[i]->Kernel;
        Body[counter].Data2 = Pall.ConvertNumberDensityToCGS*Phydro[i]->Rho;
        Body[counter].Data3 = Pall.ConvertUtoT*Phydro[i]->U;
        Body[counter].Data4 = Phydro[i]->Z;
        counter ++;
    }

    return ;
}

static void ASR_Pack4D2UFormatExtraStar(struct StructWrite4D2UBody Body[]){

    const double LengthUnitConverter = Pall.UnitLength/PC_CGS;
    const double VelocityUnitConverter = (Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;
    const double ExtraInterval = Pall.OutPutInterval/ON_THE_FLY_4D2U_EXTRA_INFO_INTERVAL;

    int counter = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(PstarBody(i)->dt*ON_THE_FLY_4D2U_EXTRA_INFO_CRITERIA >= ExtraInterval) continue;

        Body[counter].Type = TypeStar;
        Body[counter].GlobalID = PstarBody(i)->GlobalID;
        Body[counter].Tag = 0;

        Body[counter].Pos[0] = LengthUnitConverter*PstarBody(i)->Pos[0];
        Body[counter].Pos[1] = LengthUnitConverter*PstarBody(i)->Pos[1];
        Body[counter].Pos[2] = LengthUnitConverter*PstarBody(i)->Pos[2];
        Body[counter].Vel[0] = VelocityUnitConverter*PstarBody(i)->Vel[0];
        Body[counter].Vel[1] = VelocityUnitConverter*PstarBody(i)->Vel[1];
        Body[counter].Vel[2] = VelocityUnitConverter*PstarBody(i)->Vel[2];

        Body[counter].Mass  = PstarBody(i)->Mass*Pall.UnitMass/MSUN_CGS;
        Body[counter].Data1 = LengthUnitConverter*PstarBody(i)->Eps;
        Body[counter].Data2 = (Pall.TCurrent-Pstar[i]->FormationTime)*Pall.UnitTime/YEAR_CGS;//Age
        Body[counter].Data3 = 0.e0;
        Body[counter].Data4 = Pstar[i]->Z;// Metallicity

        counter ++;
    }

    return ;
}

static void ASR_Pack4D2UFormatExtraSink(struct StructWrite4D2UBody Body[]){

    const double LengthUnitConverter = Pall.UnitLength/PC_CGS;
    const double VelocityUnitConverter = (Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;
    const double ExtraInterval = Pall.OutPutInterval/ON_THE_FLY_4D2U_EXTRA_INFO_INTERVAL;

    int counter = 0;
    for(int i=0;i<Pall.Nsink;i++){
        if(PsinkBody(i)->dt*ON_THE_FLY_4D2U_EXTRA_INFO_CRITERIA >= ExtraInterval) continue;

        Body[counter].Type = TypeSink;
        Body[counter].GlobalID = PsinkBody(i)->GlobalID;
        Body[counter].Tag = 0;

        Body[counter].Pos[0] = LengthUnitConverter*PsinkBody(i)->Pos[0];
        Body[counter].Pos[1] = LengthUnitConverter*PsinkBody(i)->Pos[1];
        Body[counter].Pos[2] = LengthUnitConverter*PsinkBody(i)->Pos[2];
        Body[counter].Vel[0] = VelocityUnitConverter*PsinkBody(i)->Vel[0];
        Body[counter].Vel[1] = VelocityUnitConverter*PsinkBody(i)->Vel[1];
        Body[counter].Vel[2] = VelocityUnitConverter*PsinkBody(i)->Vel[2];

        Body[counter].Mass  = PsinkBody(i)->Mass*Pall.UnitMass/MSUN_CGS;
        Body[counter].Data1 = LengthUnitConverter*PsinkBody(i)->Eps;
        Body[counter].Data2 = Psink[i]->AccretionRadius*Pall.UnitLength/PC_CGS;
        Body[counter].Data3 = 0.e0;
        Body[counter].Data4 = Psink[i]->Z;

        counter ++;
    }

    return ;
}

static void ASR_Pack4D2UFormatExtraDM(struct StructWrite4D2UBody Body[]){

    const double LengthUnitConverter = Pall.UnitLength/PC_CGS;
    const double VelocityUnitConverter = (Pall.UnitLength/Pall.UnitTime)/VELOCITY_KMS_CGS;
    const double ExtraInterval = Pall.OutPutInterval/ON_THE_FLY_4D2U_EXTRA_INFO_INTERVAL;

    int *Index;
    Index = malloc(sizeof(int)*Pall.NDM);

    int counter = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeDM){
            if(Pbody[i]->dt*ON_THE_FLY_4D2U_EXTRA_INFO_CRITERIA >= ExtraInterval) continue;
            Index[counter] = i;
            counter ++;
        }
    }

    for(int k=0;k<counter;k++){
        int i = Index[k];
        Body[k].Type = TypeDM;
        Body[k].GlobalID = Pbody[i]->GlobalID;
        Body[k].Tag = 0;

        Body[k].Pos[0] = LengthUnitConverter*Pbody[i]->Pos[0];
        Body[k].Pos[1] = LengthUnitConverter*Pbody[i]->Pos[1];
        Body[k].Pos[2] = LengthUnitConverter*Pbody[i]->Pos[2];
        Body[k].Vel[0] = VelocityUnitConverter*Pbody[i]->Vel[0];
        Body[k].Vel[1] = VelocityUnitConverter*Pbody[i]->Vel[1];
        Body[k].Vel[2] = VelocityUnitConverter*Pbody[i]->Vel[2];

        Body[k].Mass =  Pbody[i]->Mass*Pall.UnitMass/MSUN_CGS;
        Body[k].Data1 = LengthUnitConverter*Pbody[i]->Eps;
        Body[k].Data2 = 0.e0;
        Body[k].Data3 = 0.e0;
        Body[k].Data4 = 0.e0;
    }

    free(Index);

    return ;
}

static void ASR_Pack4D2UFormatExtra(struct StructWrite4D2UBody Body[], const int Type){

    if(Type == 0){ // Hydro
        ASR_Pack4D2UFormatExtraHydro(Body);
    }else if(Type == 1){ // Star
        ASR_Pack4D2UFormatExtraStar(Body);
    }else if(Type == 2){ // Sink
        ASR_Pack4D2UFormatExtraSink(Body);
    }else if(Type == 3){ // DM
        ASR_Pack4D2UFormatExtraDM(Body);
    }

    return ;
}

//#define WRITE4D2U_VELOCITY 

static void ASR_Write4D2UFormat(FILE *fp, struct StructWrite4D2UBody Body[], const int Size){

    for(int i=0;i<Size;i++){
#ifdef WRITE4D2U_VELOCITY //{
        fprintf(fp,"%d %ld %g %g %g %g %g %g %g %g %g %g %g %d\n",
                Body[i].Type,
                Body[i].GlobalID,
                Body[i].Pos[0],Body[i].Pos[1],Body[i].Pos[2],
                Body[i].Vel[0],Body[i].Vel[1],Body[i].Vel[2],
                Body[i].Mass,
                Body[i].Data1,Body[i].Data2,Body[i].Data3,Body[i].Data4,
                Body[i].Tag);
#else // WRITE4D2U_VELOCITY //}//{
        fprintf(fp,"%d %ld %g %g %g %g %g %g %g %g %d\n",
                Body[i].Type,
                Body[i].GlobalID,
                Body[i].Pos[0],Body[i].Pos[1],Body[i].Pos[2],
                Body[i].Mass,
                Body[i].Data1,Body[i].Data2,Body[i].Data3,Body[i].Data4,
                Body[i].Tag);
#endif // WRITE4D2U_VELOCITY //}
    }

    return ;
}

static void ASR_SetSizeArray(int SizeArray[], const int Type){

    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();

    int mycount;
    if(Type == 0){ // Hydro
        mycount = Pall.Nhydro;
    }else if(Type == 1){ // Star
        mycount = Pall.Nstars;
    }else if(Type == 2){ // Sink
        mycount = Pall.Nsink;
    }else if(Type == 3){ // DM
        mycount = Pall.NDM;
    }

#if 0
    MPI_Win win_data;
    MPI_Win_create(SizeArray,NProcs*sizeof(int),sizeof(int),MPI_INFO_NULL,MPI_COMM_WORLD,&win_data);

    MPI_Win_fence(0,win_data);
    MPI_Put(&mycount,1,MPI_INT,MPI_ROOT_RANK,MyID,1,MPI_INT,win_data);
    MPI_Win_fence(0,win_data);
    MPI_Win_free(&win_data);
    if(MyID != MPI_ROOT_RANK){
        SizeArray[MyID] = mycount;
    }
#else
    MPI_Gather(&mycount,1,MPI_INT,SizeArray,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    if(MyID != MPI_ROOT_RANK){
        SizeArray[MyID] = mycount;
    }
#endif

    return ;
}

static void ASR_SetSizeArrayExtra(int SizeArray[], const int Type, const int LocalNumbers[]){

    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();

    int mycount;
    if(Type == 0){ // Hydro
        mycount = LocalNumbers[0];
    }else if(Type == 1){ // Star
        mycount = LocalNumbers[1];
    }else if(Type == 2){ // Sink
        mycount = LocalNumbers[2];
    }else if(Type == 3){ // DM
        mycount = LocalNumbers[3];
    }

    MPI_Win win_data;
    MPI_Win_create(SizeArray,NProcs*sizeof(int),sizeof(int),MPI_INFO_NULL,MPI_COMM_WORLD,&win_data);

    MPI_Win_fence(0,win_data);
    MPI_Put(&mycount,1,MPI_INT,MPI_ROOT_RANK,MyID,1,MPI_INT,win_data);
    MPI_Win_fence(0,win_data);
    MPI_Win_free(&win_data);

    if(MyID != MPI_ROOT_RANK)
        SizeArray[MyID] = mycount;

    return ;
}

void Write4D2U(void){

    // timing checker

    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();

    // Prepare buffer
    int BufferSize = MAX(Pall.Nhydro,MAX(Pall.Nstars,MAX(Pall.Nsink,Pall.NDM)));

    struct StructWrite4D2UBody *Body,*RecvBody;
    Body = malloc(sizeof(struct StructWrite4D2UBody)*BufferSize);

    MPI_Allreduce(MPI_IN_PLACE,&BufferSize,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

#ifdef USE_ON_THE_FLY_4D2U_WRITE_PARALLEL //{
    {

    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"4D2U/4D2U.tmp.%04d",Pall.OutPutFileNumber);
    FileOpen(fp,fname,"w");
    if(MyID == MPI_ROOT_RANK){
        ASR_Write4D2UFormatHeader(fp,Pall.OutPutFileNumber);
    }

    // hydro
    for(int k=0;k<4;k++){
        int Size = 0;
        if(k==TypeHydro){
            Size = Pall.Nhydro;
        } else if(k==TypeStar){
            Size = Pall.Nstars;
        } else if(k==TypeSink){
            Size = Pall.Nsink;
        } else if(k==TypeDM){
            Size = Pall.NDM;
        }
        ASR_Pack4D2UFormat(Body,k);
        ASR_Write4D2UFormat(fp,Body,Size);
    }
    fclose(fp);
    fflush(fp);


    sprintf(fname,"4D2U/4D2U.tmp.%04d",Pall.OutPutFileNumber);

    }
#else // USE_ON_THE_FLY_4D2U_WRITE_PARALLEL //}//}
    {

    if(MyID == MPI_ROOT_RANK){
        RecvBody = malloc(sizeof(struct StructWrite4D2UBody)*BufferSize);
    }

    FILE *fp;
    if(MyID == MPI_ROOT_RANK){
        char fname[MaxCharactersInLine];
        sprintf(fname,"4D2U/4D2U.%04d",Pall.OutPutFileNumber);
        FileOpen(fp,fname,"w");
        // Write header
        ASR_Write4D2UFormatHeader(fp,Pall.OutPutFileNumber);
    }

    // hydro
    for(int k=0;k<4;k++){
        int SizeArray[NProcs];
        for(int i=0;i<NProcs;i++)
            SizeArray[i] = 0;
        ASR_SetSizeArray(SizeArray,k);
        // if(MPIGetMyID() == MPI_ROOT_RANK)
            // for(int i=0;i<NProcs;i++)
                // dprintlmpi(SizeArray[i]);

        ASR_Pack4D2UFormat(Body,k);

        for(int i=0;i<NProcs;i++){
            MPI_Status mpi_status_Recv;
            MPI_Request mpi_request_Recv;
            if(MyID == MPI_ROOT_RANK){
                MPI_Irecv(RecvBody,SizeArray[i]*sizeof(struct StructWrite4D2UBody),
                    MPI_BYTE,i,0,MPI_COMM_WORLD,&mpi_request_Recv);
            }
            MPI_Status mpi_status_Send;
            MPI_Request mpi_request_Send;
            if(MyID == i){
                // dprintlmpi(SizeArray[i]);fflush(NULL);
                MPI_Isend(Body,SizeArray[i]*sizeof(struct StructWrite4D2UBody),
                    MPI_BYTE,MPI_ROOT_RANK,0,MPI_COMM_WORLD,&mpi_request_Send);
                MPI_Wait(&mpi_request_Send,&mpi_status_Send);
            } 

            if(MyID == MPI_ROOT_RANK){
                MPI_Wait(&mpi_request_Recv,&mpi_status_Recv);
                // WriteData
                ASR_Write4D2UFormat(fp,RecvBody,SizeArray[i]);
                fflush(NULL);
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }
    }


    if(MyID == MPI_ROOT_RANK){
        fclose(fp);
    }

    free(Body);
    if(MyID == MPI_ROOT_RANK){
        free(RecvBody);
    }
    }
#endif // USE_ON_THE_FLY_4D2U_WRITE_PARALLEL //}

    // Reset counter.
    Write4D2UExtraCounter = 1;

    return ;
}

static int Evaluate4D2UExtraInfo(const int Mode, const int Type){

    const double ExtraInterval = Pall.OutPutInterval/ON_THE_FLY_4D2U_EXTRA_INFO_INTERVAL;
    //const double dt_criteria = ON_THE_FLY_4D2U_EXTRA_INFO_CRITERIA*ExtraInterval;

    if(Mode == 0){ // count
        int counter = 0;
        if(Type == TypeHydro){
            for(int i=0;i<Pall.Nhydro;i++){
                if(Phydro[i]->dt_hydro*ON_THE_FLY_4D2U_EXTRA_INFO_CRITERIA < ExtraInterval){
                    counter ++;
                }
            }
        } else {
            for(int i=0;i<Pall.Ntotal;i++){
                if(Pbody[i]->Type == Type){
                    if(Pbody[i]->dt*ON_THE_FLY_4D2U_EXTRA_INFO_CRITERIA < ExtraInterval){
                        counter ++;
                    }
                }
            }
        }
        return counter;
    } else {
    }


    return 0;
}

void ResetWrite4D2UExtraCounter(void){

    const double ExtraInterval = Pall.OutPutInterval/ON_THE_FLY_4D2U_EXTRA_INFO_INTERVAL;
    while(Pall.TCurrent > 
        (Pall.OutPutFileNumber-1)*Pall.OutPutInterval+Write4D2UExtraCounter*ExtraInterval){
        Write4D2UExtraCounter ++;
    }

    return ;
}


void Write4D2UExtraInfo(void){

#if 0
#ifdef COSMOLOGICAL_RUN //{
    double Toffset = CalcZtoT(Pall.InitialRedshift);
    if(Pall.TCurrent-Toffset+0.5*Pall.dtmin >= Pall.OutPutFileNumber*Pall.OutPutInterval)
#else // COSMOLOGICAL_RUN //}//{
    if(Pall.TCurrent >= Pall.OutPutFileNumber*Pall.OutPutInterval)
#endif // COSMOLOGICAL_RUN //}
#endif

    const int MyID = MPIGetMyID();
    const int NProcs = MPIGetNumProcs();

    const double ExtraInterval = Pall.OutPutInterval/ON_THE_FLY_4D2U_EXTRA_INFO_INTERVAL;
    if(Pall.TCurrent >= 
        (Pall.OutPutFileNumber-1)*Pall.OutPutInterval+Write4D2UExtraCounter*ExtraInterval){

        // Prepare buffer
        int LocalNumbers[] = {Evaluate4D2UExtraInfo(0,TypeHydro),
                              Evaluate4D2UExtraInfo(0,TypeStar),
                              Evaluate4D2UExtraInfo(0,TypeSink),
                              Evaluate4D2UExtraInfo(0,TypeDM)};
        int Numbers[4] = {LocalNumbers[0],LocalNumbers[1],
                          LocalNumbers[2],LocalNumbers[3]};

        int BufferSize = MAX(Numbers[0],MAX(Numbers[1],MAX(Numbers[2],Numbers[3])));

        struct StructWrite4D2UBody *Body,*RecvBody;
        Body = malloc(sizeof(struct StructWrite4D2UBody)*BufferSize);

        MPI_Allreduce(MPI_IN_PLACE,&BufferSize,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&Numbers,4,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        if(MyID == MPI_ROOT_RANK){
            RecvBody = malloc(sizeof(struct StructWrite4D2UBody)*BufferSize);
        }

        FILE *fp;
        if(MyID == MPI_ROOT_RANK){
            char fname[MaxCharactersInLine];
            sprintf(fname,"4D2U/4D2U.%04d.Ext.%03d",Pall.OutPutFileNumber,Write4D2UExtraCounter);
            FileOpen(fp,fname,"w");
            ASR_Write4D2UFormatHeaderExtra(fp,Pall.OutPutFileNumber,Numbers);
        }

        // hydro
        for(int k=0;k<4;k++){
            int SizeArray[NProcs];
            ASR_SetSizeArrayExtra(SizeArray,k,LocalNumbers);
            ASR_Pack4D2UFormatExtra(Body,k);

            for(int i=0;i<NProcs;i++){
                MPI_Status mpi_status_Recv;
                MPI_Request mpi_request_Recv;
                if(MyID == MPI_ROOT_RANK){
                    MPI_Irecv(RecvBody,SizeArray[i]*sizeof(struct StructWrite4D2UBody),
                        MPI_BYTE,i,0,MPI_COMM_WORLD,&mpi_request_Recv);
                }
                MPI_Status mpi_status_Send;
                MPI_Request mpi_request_Send;
                if(MyID == i){
                    MPI_Isend(Body,SizeArray[i]*sizeof(struct StructWrite4D2UBody),MPI_BYTE,
                            MPI_ROOT_RANK,0,MPI_COMM_WORLD,&mpi_request_Send);
                    MPI_Wait(&mpi_request_Send,&mpi_status_Send);
                } 

                if(MyID == MPI_ROOT_RANK){
                    MPI_Wait(&mpi_request_Recv,&mpi_status_Recv);
                    // WriteData
                    ASR_Write4D2UFormat(fp,RecvBody,SizeArray[i]);
                    fflush(NULL);
                }

                MPI_Barrier(MPI_COMM_WORLD);
            }
        }

        if(MyID == MPI_ROOT_RANK){
            fclose(fp);
        }

        free(Body);
        if(MyID == MPI_ROOT_RANK){
            free(RecvBody);
        }

        Write4D2UExtraCounter ++;
    }

    return ;
}
#endif // USE_ON_THE_FLY_4D2U //}
