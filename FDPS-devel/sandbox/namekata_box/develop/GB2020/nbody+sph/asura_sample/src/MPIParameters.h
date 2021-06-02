#pragma once 

void InitializeMPIEnv(int *argc, char **argv);
void FinalizeMPIEnv(void);
void MPISetXYZFactors(void);
void MPISetMyID(int current_id);
void MPISetNumProcs(int current_numprocs);
void MPISetNumGrapes(int current_numgrapes);
void MPISetNameLen(int current_namelen);
void MPISetProcessorName(char *current_processor_name);
void MPISetNumProcsPower(int current_numprocs);
int MPIGetMyID(void);
int MPIGetNumProcs(void);
int MPIGetNumGrapes(void);
int MPIGetNameLen(void);
char *MPIGetProcessorName(void);
int MPIGetNumProcsPower(void);
void MPIGetXYZFactors(int FactorsDist[]);
int MPIGetFactor(const int direction);
int MPIGetFactorFromID(const int direction);


#ifdef USE_BARYON_COMM //{ 
int MPIGetHydroMyID(void);
int MPIGetHydroNumProcs(void);
int MPIGetActiveHydroMyID(void);
int MPIGetActiveHydroNumProcs(void);
void CreateHydroGroup(void);
void DeleteHydroGroup(void);
void CreateHydroActiveGroup(void);
void DeleteCreateHydroActiveGroup(void);
int MPIGetBaryonMyID(void);
int MPIGetBaryonNumProcs(void);
int MPIGetActiveBaryonMyID(void);
int MPIGetActiveBaryonNumProcs(void);
void CreateBaryonGroup(void);
void DeleteBaryonGroup(void);
void CreateActiveBaryonGroup(void);
void DeleteActiveBaryonGroup(void);
int RankTranslator(MPI_Group Group1, int Rank1, MPI_Group Group2);

#endif // USE_BARYON_COMM //;

#define bprintlmpi(f)       printf("*** " #f " = %d:Rank%d:%s:line%d:%s()\n", (int)f,MPIGetMyID(),__FILE__,__LINE__,__FUNCTION__)
#define dprintlmpi(f)       printf("*** " #f " = %d:Rank%d:%s:line%d:%s()\n", f,MPIGetMyID(),__FILE__,__LINE__,__FUNCTION__)
#define zdprintlmpi(f)      printf("*** " #f " = %zd:Rank%d:%s:line%d:%s()\n", f,MPIGetMyID(),__FILE__,__LINE__,__FUNCTION__)
#define dlprintlmpi(f)      printf("*** " #f " = %ld:Rank%d:%s:line%d:%s()\n", f,MPIGetMyID(),__FILE__,__LINE__,__FUNCTION__)
#define fprintlmpi(f)       printf("*** " #f " = %f:Rank%d:%s:line%d:%s()\n", f,MPIGetMyID(),__FILE__,__LINE__,__FUNCTION__)
#define eprintlmpi(f)       printf("*** " #f " = %e:Rank%d:%s:line%d:%s()\n", f,MPIGetMyID(),__FILE__,__LINE__,__FUNCTION__)
#define gprintlmpi(f)       printf("*** " #f " = %g:Rank%d:%s:line%d:%s()\n", f,MPIGetMyID(),__FILE__,__LINE__,__FUNCTION__)
#define sprintlmpi(f)       printf("*** " #f " %s:Rank%d:line%d:%s()\n", __FILE__,MPIGetMyID(),__LINE__,__FUNCTION__)
#define pprintlmpi(f)       printf("*** " #f " = %p:Rank%d:%s:line%d:%s()\n", f,MPIGetMyID(),__FILE__,__LINE__,__FUNCTION__)
