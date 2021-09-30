#include "config.h"

char BaseNameOfParameterLogs[] = "./ParameterLogs";

void OutputParameterLogs(void){

    char FileName[MaxCharactersInLine];

    int Counter = 0;
    do{
        Snprintf(FileName,"%s.%03d",BaseNameOfParameterLogs,Counter);
        Counter ++;
    } while (CheckFile(FileName) == true);

    FILE *fp;
    FileOpen(fp,FileName,"w");

    fclose(fp);

    return ;
}
