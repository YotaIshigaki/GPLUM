/* NeighborSearch.c */
int GetNeighbors(double Pos[restrict], const double h, int list[restrict]);
int GetNeighborsLimited(double Pos[restrict], const double h, int list[restrict]);
int GetNeighborsLimitedImported(double Pos[restrict], const double h, int list[restrict]);
int GetActiveNeighborsLimited(double Pos[restrict], const double h, int list[restrict]);
int GetNeighborsIterativeApproach(const int StartNodeID, double Pos[restrict], const double h, int *nlist, int list[restrict]);
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
int GetNeighborsSmoothedNumberIterativeApproach(const int StartNodeID, double Pos[restrict], const double h, int *nlist, int list[restrict], double *SmoothedNumber);
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
int GetNeighborsIterativeApproachFOF(const int StartNodeID, double Pos[restrict], const double LinkingLength, int *nlist, int list[restrict]);
int GetNeighborNumbers(double Pos[restrict], const double h);
int GetNeighborsPairs(double Pos[restrict], const double h, int list[restrict]);
int GetNeighborsPairsLimited(double Pos[restrict], const double h, int list[restrict]);
int GetNeighborsPairsLimitedImported(double Pos[restrict], const double h, int list[restrict]);
int GetNeighborsPairsInteractiveApproach(const int StartNodeID, double Pos[restrict], const double h, int *nlist, int list[restrict]);
int GetNeighborsDirect(double Pos[restrict], const double h, int list[restrict]);
int GetNeighborsPairsDirect(double Pos[restrict], const double h, int list[restrict]);
int GetNeighborsPairsDirectLimitedImported(double Pos[restrict], const double h, int list[restrict]);
int ReturnNeighborNumber(const int CurrentNodeID, double Pos[restrict], const double h);
int GetNeighborsLimitedNBCacheIndex(double Pos[restrict], const double h, int list[restrict]);
int GetNeighborsPairsLimitedNBCacheIndex(double Pos[restrict], const double h, int list[restrict]);

int GetNeighborsFromStellarTree(double Pos[restrict], const double h, int list[restrict]);
int GetNeighborsLimitedFromStellarTree(double Pos[restrict], const double h, int list[restrict]);
int GetActiveNeighborsLimitedFromStellarTree(double Pos[restrict], const double h, int list[restrict]);
int GetNeighborsIterativeApproachFromStellarTree(const int StartNodeID, 
        double Pos[restrict], const double h, int *nlist, int list[restrict]);
