#include "config.h"
#include <sys/stat.h>

bool CheckFile(const char FileName[]){
    FILE *fp;
    fp = fopen(FileName,"r");
    if(fp == NULL){
        return false;
    }
    fclose(fp);
    return true;
}

bool CheckDir(const char DirName[]){
    FILE *fp;
    fp = fopen(DirName,"r");
    if(fp == NULL){
        return false;
    }
    fclose(fp);
    return true;
}

void MakeDir(const char DirName[]){

    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(CheckDir(DirName) == false){
            int checkflag = mkdir(DirName,0755);
            if(checkflag < 0){
                fprintf(stderr,"Directory [%s] creation error.\n",DirName);
                MPI_Finalize();
                exit(DirectoryCreationError);
            } else {
                fprintf(stderr,"Directory [%s] is created.\n",DirName);
            }
        } 
    }

    return ;
}

/*
 * This function returns "true" when the each string which is obtained from *src
 * splited by *delimiter has *word. 
 */
bool FindWordFromString(char *src, char *word, char *delimiter){

    char src_copy[MaxCharactersInLine];
    strcpy(src_copy,src);
    char *tok;
    tok = strtok(src_copy,delimiter);
    while(tok != NULL){
        //dbg("%s\n",tok);
        if(strcmp(tok,word) == 0)
            return true;
        tok = strtok(NULL,delimiter);
    }
    return false;
}

bool GetBaseFileName(char *src, char *word, char *delimiter, char *BaseFileName){

    char src_copy[MaxCharactersInLine];
    strcpy(src_copy,src);
    char *tok;
    char oldtok[MaxCharactersInLine];
    tok = strtok(src_copy,delimiter);
    while(tok != NULL){
        // fprintf(stderr,"tok = %s, %s\n",tok,oldtok);
        if(strcmp(tok,word) == 0){
            // fprintf("%s %s",tok,oldtok);
            strcpy(BaseFileName,oldtok);
            // fflush(NULL);
            return true;
        }
        strcpy(oldtok,tok);
        tok = strtok(NULL,delimiter);
    }
    return false;
}
