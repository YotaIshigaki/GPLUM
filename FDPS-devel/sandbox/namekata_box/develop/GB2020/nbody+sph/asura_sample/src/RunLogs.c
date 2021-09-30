/*
 * This file includes "log output rountines".
 * "WriteMakeFlags" writes a file "makefiles", which includes all of parameters in makefiles
 * On the other hand, "" writes a file "", 
 *  which includes important members in the structure "Pall".
 * Both subroutines write "*.???" where "???" means times of runs.
 * In the first run, the label is "000".
 */

#include "config.h"

char makefiles[] = "makeflags";
char parameters[] = "Pall";

#include "RunLogsAG.h"

void WriteMakeflags(void){

    if(MPIGetMyID()==MPI_ROOT_RANK){
        char FileName[MaxCharactersInLine];

        for(int i=0;i<999;i++){
            sprintf(FileName,"%s.%03d",makefiles,i);
            if(CheckFile(FileName) == false ){ // open and write.
                WriteMakeflagsAG(FileName);
                return ;
            }
        }
    }
    return ;
}
