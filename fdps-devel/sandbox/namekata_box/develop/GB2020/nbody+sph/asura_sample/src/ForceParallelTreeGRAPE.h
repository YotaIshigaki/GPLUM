/* ForceParallelTreeGRAPE.c */
void CCC(const int ID);
void DDD(const char *func, const int line);

void InitializeRootForLET(void);
void InitializeParallelTreeGRAPE(void);
void UpdateNumberofLeavesInGroupInLET(const int Ng);
void ForceParallelTreeGRAPE(void);

#ifdef TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR
void ForceParallelTreeGRAPEInsert(const int NImport, double Pos[restrict][3], double Mass[restrict], double Eps[restrict]);
void CalcGravityDirectSymmetrizedPotentialWithPhantom(void);
#endif
