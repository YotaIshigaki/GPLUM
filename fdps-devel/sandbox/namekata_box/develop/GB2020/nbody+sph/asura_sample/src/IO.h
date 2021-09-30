#pragma once 

void GetRunStatus(const int argc, char **argv);
int GetRunStatusForSingleNodeAnalysis(const int argc, char **argv);
void CheckIOFileName(void);
void WriteAllData(void);
void ReadAllData(void);
void DataDump(void);
void DataFullDump(void);
void FileOutPutConstantInterval(void);
void BinaryDump(void);
void ParallelWriteAllData(char fname[]);
void ParallelDumpAllData(char fname[]);
void ParallelWriteAllDataASCIIFormat(char fname[]);
void ParallelReadAllData(void);
void ReadParallelDataOnSingleNode(void);
void OutPutASCIIDATA(void);
void WriteHydroDataASCIIFormat(const int mode);
void ShowHydroDataASCIIFormat(const int mode);
void WriteNeighborInformation(const int index, const int mode);
void ShowNeighborInformation(const int index, const int mode);
void WriteTimeSteps(void);
void WriteCurrentActiveParticles(const long int FileID, char suffix[]);

#ifdef USE_ON_THE_FLY_4D2U //{
void InitWrite4D2U(void);
void ResetWrite4D2UExtraCounter(void);
void Write4D2U(void);
void Write4D2UExtraInfo(void);
#endif // USE_ON_THE_FLY_4D2U //}
